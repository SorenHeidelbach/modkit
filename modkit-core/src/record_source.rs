//! Sources of alignment records for region-based workers (pileup, threshold
//! sampling). A worker asks for the records overlapping a reference interval
//! and gets them in coordinate-sorted order, either from an indexed BAM/CRAM
//! or from an in-memory store of aligned records ([`AlignedRecordStore`]).
//! The store lets `modkit pileup --align` pile up reads straight after
//! aligning them, without writing, sorting, indexing and re-reading a BAM.

use std::path::PathBuf;
use std::sync::Arc;

use anyhow::bail;
use rayon::prelude::*;
use rust_htslib::bam::{self, FetchDefinition, Read};
use rust_htslib::htslib;

use crate::util::reader_is_cram;

pub type RecordIter<'a> = Box<
    dyn Iterator<Item = Result<bam::Record, rust_htslib::errors::Error>> + 'a,
>;

pub trait RegionRecordSource: Send {
    /// Name of the reference sequence with this target id.
    fn target_name(&self, tid: u32) -> String;
    /// Records overlapping `[start, end)` on `tid`, in the order a
    /// coordinate-sorted BAM would return them.
    fn fetch_records<'a>(
        &'a mut self,
        tid: u32,
        start: u32,
        end: u32,
    ) -> anyhow::Result<RecordIter<'a>>;
}

impl RegionRecordSource for bam::IndexedReader {
    fn target_name(&self, tid: u32) -> String {
        String::from_utf8_lossy(self.header().tid2name(tid)).to_string()
    }

    fn fetch_records<'a>(
        &'a mut self,
        tid: u32,
        start: u32,
        end: u32,
    ) -> anyhow::Result<RecordIter<'a>> {
        self.fetch(FetchDefinition::Region(
            tid as i32,
            start as i64,
            end as i64,
        ))?;
        Ok(Box::new(self.records()))
    }
}

/// Where workers get their records from.
#[derive(Clone)]
pub enum RecordSourceSpec {
    /// An indexed BAM or CRAM file.
    Hts(PathBuf),
    /// Aligned records held in memory.
    Store(Arc<AlignedRecordStore>),
}

impl RecordSourceSpec {
    /// Open a source for one worker.
    pub fn open(
        &self,
        reference_fp: Option<&PathBuf>,
        thread_pool: Option<&rust_htslib::tpool::ThreadPool>,
    ) -> anyhow::Result<Box<dyn RegionRecordSource>> {
        match self {
            Self::Hts(fp) => {
                let mut reader = bam::IndexedReader::from_path(fp)?;
                if let Some(pool) = thread_pool {
                    reader.set_thread_pool(pool)?;
                }
                if reader_is_cram(&reader) {
                    let Some(reference_fp) = reference_fp else {
                        bail!("CRAM input requires reference")
                    };
                    reader.set_reference(reference_fp)?;
                    reader.set_cram_options(
                        htslib::hts_fmt_option_CRAM_OPT_REQUIRED_FIELDS,
                        htslib::sam_fields_SAM_FLAG
                            | htslib::sam_fields_SAM_RNAME
                            | htslib::sam_fields_SAM_POS
                            | htslib::sam_fields_SAM_MAPQ
                            | htslib::sam_fields_SAM_CIGAR
                            | htslib::sam_fields_SAM_SEQ
                            | htslib::sam_fields_SAM_AUX,
                    )?;
                }
                Ok(Box::new(reader))
            }
            Self::Store(store) => Ok(Box::new(StoreSource(store.clone()))),
        }
    }
}

struct StoredRecord {
    pos: i64,
    /// End position as computed by htslib (`bam_endpos`), used for overlap
    /// queries exactly like an htslib region fetch.
    end: i64,
    reverse: bool,
    /// Order the record was produced in, used to break ties when sorting
    /// so the order matches `samtools sort` of the aligner output.
    seq: u64,
    record: bam::Record,
}

#[derive(Default)]
struct StoreContig {
    records: Vec<StoredRecord>,
    max_span: i64,
}

/// Aligned records kept in memory, sorted by coordinate per contig.
pub struct AlignedRecordStore {
    target_names: Vec<String>,
    header_text: Vec<u8>,
    contigs: Vec<StoreContig>,
}

impl AlignedRecordStore {
    pub fn header_text(&self) -> &[u8] {
        &self.header_text
    }

    pub fn n_records(&self, tid: u32) -> usize {
        self.contigs.get(tid as usize).map(|c| c.records.len()).unwrap_or(0)
    }

    pub fn total_records(&self) -> usize {
        self.contigs.iter().map(|c| c.records.len()).sum()
    }

    /// Approximate size of the record data held, in bytes.
    pub fn approx_bytes(&self) -> usize {
        self.contigs
            .iter()
            .flat_map(|c| c.records.iter())
            .map(|r| {
                r.record.inner().l_data as usize
                    + std::mem::size_of::<StoredRecord>()
            })
            .sum()
    }

    fn records_overlapping(
        &self,
        tid: u32,
        start: u32,
        end: u32,
    ) -> impl Iterator<Item = &bam::Record> + '_ {
        let (records, lo) = match self.contigs.get(tid as usize) {
            Some(contig) => {
                let lower = start as i64 - contig.max_span;
                let lo = contig.records.partition_point(|r| r.pos < lower);
                (&contig.records[..], lo)
            }
            None => (&[][..], 0),
        };
        let (start, end) = (start as i64, end as i64);
        records[lo..]
            .iter()
            .take_while(move |r| r.pos < end)
            .filter(move |r| r.end > start)
            .map(|r| &r.record)
    }
}

pub struct AlignedRecordStoreBuilder {
    target_names: Vec<String>,
    header_text: Vec<u8>,
    contigs: Vec<Vec<StoredRecord>>,
    next_seq: u64,
}

impl AlignedRecordStoreBuilder {
    pub fn new(target_names: Vec<String>, header_text: Vec<u8>) -> Self {
        let contigs = (0..target_names.len()).map(|_| Vec::new()).collect();
        Self { target_names, header_text, contigs, next_seq: 0 }
    }

    /// Add a record, records must be pushed in output order. Unmapped
    /// records are ignored.
    pub fn push(&mut self, record: bam::Record) -> anyhow::Result<()> {
        let seq = self.next_seq;
        self.next_seq += 1;
        if record.is_unmapped() || record.tid() < 0 {
            return Ok(());
        }
        let tid = record.tid() as usize;
        let Some(contig) = self.contigs.get_mut(tid) else {
            bail!("record has invalid target id {tid}")
        };
        let end = unsafe { htslib::bam_endpos(record.inner()) };
        contig.push(StoredRecord {
            pos: record.pos(),
            end,
            reverse: record.is_reverse(),
            seq,
            record,
        });
        Ok(())
    }

    pub fn finish(self) -> AlignedRecordStore {
        let contigs = self
            .contigs
            .into_par_iter()
            .map(|mut records| {
                // same order as `samtools sort`: position, then strand, ties
                // kept in input order
                records.sort_unstable_by_key(|r| (r.pos, r.reverse, r.seq));
                let max_span =
                    records.iter().map(|r| r.end - r.pos).max().unwrap_or(0);
                StoreContig { records, max_span }
            })
            .collect();
        AlignedRecordStore {
            target_names: self.target_names,
            header_text: self.header_text,
            contigs,
        }
    }
}

struct StoreSource(Arc<AlignedRecordStore>);

impl RegionRecordSource for StoreSource {
    fn target_name(&self, tid: u32) -> String {
        self.0.target_names.get(tid as usize).cloned().unwrap_or_default()
    }

    fn fetch_records<'a>(
        &'a mut self,
        tid: u32,
        start: u32,
        end: u32,
    ) -> anyhow::Result<RecordIter<'a>> {
        Ok(Box::new(
            self.0.records_overlapping(tid, start, end).map(|r| Ok(r.clone())),
        ))
    }
}
