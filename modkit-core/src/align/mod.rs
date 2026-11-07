//! Align reads in a modBAM to a reference with minimap2 while keeping the
//! base modification tags (MM/ML) valid.
//!
//! The usual approach (`samtools fastq -T MM,ML | minimap2 -y`) works, but it
//! has a few sharp edges: minimap2 hard-clips supplementary alignments while
//! copying the full-length MM/ML tags, stale alignment tags (NM, MD, ...) get
//! duplicated when copying all tags, legacy `Mm`/`Ml` tags are easy to drop,
//! and the `MN` tag is never set. This command runs minimap2 in-process
//! (via the `minimap2` crate) and builds the output records directly from the
//! input records so that:
//!
//! * every output record carries the full read sequence (supplementary and
//!   secondary alignments are soft-clipped), so MM/ML always apply to SEQ,
//! * all per-read tags (MM/ML, mv, RG, ...) are carried over, while tags that
//!   describe a previous alignment are dropped,
//! * `MN` is set to the sequence length whenever modification tags are
//!   present,
//! * reads that do not align are kept (flagged unmapped) unless requested
//!   otherwise.

use std::ffi::CStr;
use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};

use anyhow::{anyhow, bail, Context};
use bio::alphabets::dna::revcomp;
use clap::{Args, ValueEnum};
use crossbeam_channel::bounded;
use log::{debug, info, warn};
use minimap2::{Aligner, Built, Mapping, Preset, Strand};
use rayon::prelude::*;
use rust_htslib::bam::record::{Aux, Cigar, CigarString};
use rust_htslib::bam::{self, HeaderView, Read};

use modkit_logging::init_logging;
use safe_record::SafeRecord;

use crate::command_utils::{get_bam_writer, get_serial_reader};
use crate::util::{add_modkit_pg_records, get_ticker};

/// Tags that describe an alignment, these are removed from the input record
/// because they will be stale (or recomputed) after re-alignment.
const ALIGNMENT_TAGS: [&[u8; 2]; 21] = [
    b"NM", b"MD", b"AS", b"SA", b"XS", b"XA", b"MC", b"MQ", b"ms", b"nn",
    b"tp", b"cm", b"s1", b"s2", b"de", b"rl", b"zd", b"cs", b"cg", b"HP",
    b"PS",
];
const MN_TAG: &[u8; 2] = b"MN";
const MOD_TAGS: [&[u8; 2]; 4] = [b"MM", b"Mm", b"ML", b"Ml"];

const FLAG_REVERSE: u16 = 0x10;
const FLAG_UNMAPPED: u16 = 0x4;
const FLAG_SECONDARY: u16 = 0x100;
const FLAG_QC_FAIL: u16 = 0x200;
const FLAG_DUPLICATE: u16 = 0x400;
const FLAG_SUPPLEMENTARY: u16 = 0x800;

#[derive(Copy, Clone, Debug, ValueEnum)]
pub enum AlignPreset {
    /// Oxford Nanopore reads (minimap2 `map-ont`).
    #[value(name = "map-ont")]
    MapOnt,
    /// High-quality long reads, e.g. ONT Q20+ (minimap2 `lr:hq`).
    #[value(name = "lr:hq")]
    LrHq,
    /// PacBio HiFi/CCS reads (minimap2 `map-hifi`).
    #[value(name = "map-hifi")]
    MapHifi,
    /// PacBio CLR reads (minimap2 `map-pb`).
    #[value(name = "map-pb")]
    MapPb,
    /// Long-read spliced alignment (minimap2 `splice`).
    #[value(name = "splice")]
    Splice,
    /// High-quality long-read spliced alignment (minimap2 `splice:hq`).
    #[value(name = "splice:hq")]
    SpliceHq,
    /// Assembly to reference, ~0.1% divergence (minimap2 `asm5`).
    #[value(name = "asm5")]
    Asm5,
    /// Assembly to reference, ~1% divergence (minimap2 `asm10`).
    #[value(name = "asm10")]
    Asm10,
    /// Assembly to reference, ~5% divergence (minimap2 `asm20`).
    #[value(name = "asm20")]
    Asm20,
}

impl AlignPreset {
    fn to_minimap2(self) -> Preset {
        match self {
            Self::MapOnt => Preset::MapOnt,
            Self::LrHq => Preset::LrHq,
            Self::MapHifi => Preset::MapHifi,
            Self::MapPb => Preset::MapPb,
            Self::Splice => Preset::Splice,
            Self::SpliceHq => Preset::SpliceHq,
            Self::Asm5 => Preset::Asm5,
            Self::Asm10 => Preset::Asm10,
            Self::Asm20 => Preset::Asm20,
        }
    }
}

#[derive(Args)]
#[command(arg_required_else_help = true)]
pub struct AlignModBam {
    /// Input BAM/CRAM/SAM with base modification tags. May be unaligned or
    /// aligned (aligned reads are re-aligned from their original
    /// orientation; secondary and supplementary records are skipped). Can be
    /// a path to a file or one of `-` or `stdin` to specify a stream from
    /// standard input.
    in_bam: String,
    /// Output BAM, can be a path to a file or one of `-` or `stdout` to
    /// specify a stream to standard output. The output is not sorted.
    out_bam: String,
    /// Reference sequence to align to, a FASTA file (optionally gzipped) or
    /// a pre-built minimap2 index (.mmi).
    #[arg(short = 'r', long = "ref")]
    reference: PathBuf,
    /// Alignment preset, same as minimap2 `-x`.
    #[arg(short = 'x', long, value_enum, default_value_t = AlignPreset::MapOnt)]
    preset: AlignPreset,
    /// Number of threads to use.
    #[clap(help_heading = "Compute Options")]
    #[arg(short, long, default_value_t = 4)]
    threads: usize,
    /// Number of reads to align per batch.
    #[clap(help_heading = "Compute Options")]
    #[arg(long, default_value_t = 2_000, hide_short_help = true)]
    batch_size: usize,
    /// Output secondary alignments. Secondary records carry the full read
    /// sequence so that their MM/ML tags remain valid.
    #[clap(help_heading = "Output Options")]
    #[arg(long, default_value_t = false)]
    secondary: bool,
    /// Do not output reads that fail to align.
    #[clap(help_heading = "Output Options")]
    #[arg(long, default_value_t = false)]
    no_unmapped: bool,
    /// Add MD tags to aligned records.
    #[clap(help_heading = "Output Options")]
    #[arg(long, default_value_t = false)]
    md: bool,
    /// Output SAM format instead of BAM.
    #[clap(help_heading = "Output Options")]
    #[arg(long, default_value_t = false)]
    output_sam: bool,
    /// Output debug logs to file at this path.
    #[clap(help_heading = "Logging Options")]
    #[arg(long, alias = "log")]
    log_filepath: Option<PathBuf>,
}

#[derive(Default)]
struct Counts {
    reads_in: AtomicUsize,
    skipped_secondary_supplementary: AtomicUsize,
    skipped_no_seq: AtomicUsize,
    skipped_hard_clipped: AtomicUsize,
    failed: AtomicUsize,
    mapped: AtomicUsize,
    unmapped: AtomicUsize,
    records_out: AtomicUsize,
}

impl AlignModBam {
    pub fn run(&self) -> anyhow::Result<()> {
        let _handle = init_logging(self.log_filepath.as_ref());
        if self.threads == 0 {
            bail!("threads must be at least 1");
        }
        if self.batch_size == 0 {
            bail!("batch-size must be at least 1");
        }

        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.threads)
            .build()?;

        info!(
            "loading/building minimap2 index from {:?} with preset {:?}",
            self.reference, self.preset
        );
        let aligner =
            build_aligner(&self.reference, self.preset, self.threads)?;
        let targets = aligner_targets(&aligner);
        info!("index contains {} reference sequence(s)", targets.len());

        let mut reader = get_serial_reader(&self.in_bam)?;
        let io_threads = self.threads.min(4);
        reader.set_threads(io_threads)?;
        let mut header = build_output_header(reader.header(), &targets)?;
        add_modkit_pg_records(&mut header);
        let mut writer =
            get_bam_writer(&self.out_bam, &header, self.output_sam)?;
        writer.set_threads(io_threads)?;

        let counts = Counts::default();
        let ticker = get_ticker();
        ticker.set_message("reads aligned");

        let (batch_snd, batch_rcv) = bounded::<Vec<SafeRecord>>(2);
        let (out_snd, out_rcv) = bounded::<Vec<SafeRecord>>(2);
        let settings = MapSettings {
            secondary: self.secondary,
            unmapped: !self.no_unmapped,
            md: self.md,
        };

        std::thread::scope(|scope| -> anyhow::Result<()> {
            let counts = &counts;
            let batch_size = self.batch_size;
            let reader_handle = scope.spawn(move || -> anyhow::Result<()> {
                let mut batch = Vec::with_capacity(batch_size);
                let mut record = bam::Record::new();
                while let Some(res) = reader.read(&mut record) {
                    res.context("failed to read input record")?;
                    counts.reads_in.fetch_add(1, Ordering::Relaxed);
                    match make_template(&record) {
                        Ok(template) => batch.push(SafeRecord::from(template)),
                        Err(skip) => {
                            let counter = match skip {
                                Skip::SecondaryOrSupplementary => {
                                    &counts.skipped_secondary_supplementary
                                }
                                Skip::NoSequence => &counts.skipped_no_seq,
                                Skip::HardClipped => {
                                    &counts.skipped_hard_clipped
                                }
                            };
                            counter.fetch_add(1, Ordering::Relaxed);
                        }
                    }
                    if batch.len() >= batch_size {
                        let full = std::mem::replace(
                            &mut batch,
                            Vec::with_capacity(batch_size),
                        );
                        if batch_snd.send(full).is_err() {
                            // downstream hung up (error), stop reading
                            return Ok(());
                        }
                    }
                }
                if !batch.is_empty() {
                    let _ = batch_snd.send(batch);
                }
                Ok(())
            });

            let ticker = &ticker;
            let writer_handle = scope.spawn(move || -> anyhow::Result<()> {
                for batch in out_rcv {
                    for safe_record in batch {
                        let record: bam::Record = safe_record.into();
                        writer
                            .write(&record)
                            .context("failed to write output record")?;
                        counts.records_out.fetch_add(1, Ordering::Relaxed);
                    }
                }
                Ok(())
            });

            let aligner = &aligner;
            let settings = &settings;
            let mut align_result = Ok(());
            for batch in batch_rcv {
                let n_reads = batch.len() as u64;
                let aligned = pool.install(|| {
                    batch
                        .into_par_iter()
                        .map(|template| {
                            align_read(template.into(), aligner, settings)
                        })
                        .collect::<Vec<_>>()
                });
                let mut out = Vec::with_capacity(aligned.len());
                for res in aligned {
                    match res {
                        Ok(AlignOutcome::Mapped(records)) => {
                            counts.mapped.fetch_add(1, Ordering::Relaxed);
                            out.extend(
                                records.into_iter().map(SafeRecord::from),
                            );
                        }
                        Ok(AlignOutcome::Unmapped(record)) => {
                            counts.unmapped.fetch_add(1, Ordering::Relaxed);
                            if let Some(record) = record {
                                out.push(SafeRecord::from(record));
                            }
                        }
                        Err(e) => {
                            debug!("failed to align read, {e}");
                            counts.failed.fetch_add(1, Ordering::Relaxed);
                        }
                    }
                }
                ticker.inc(n_reads);
                if out_snd.send(out).is_err() {
                    align_result = Err(anyhow!("output writer stopped"));
                    break;
                }
            }
            drop(out_snd);
            let writer_result = writer_handle
                .join()
                .map_err(|_| anyhow!("writer thread panicked"))?;
            let reader_result = reader_handle
                .join()
                .map_err(|_| anyhow!("reader thread panicked"))?;
            // report the most specific error first
            writer_result?;
            reader_result?;
            align_result
        })?;
        ticker.finish_and_clear();

        let get = |c: &AtomicUsize| c.load(Ordering::Relaxed);
        info!(
            "processed {} input records, {} reads aligned, {} unaligned{}, \
             wrote {} records",
            get(&counts.reads_in),
            get(&counts.mapped),
            get(&counts.unmapped),
            if self.no_unmapped { " (not written)" } else { "" },
            get(&counts.records_out),
        );
        let skipped_secondary = get(&counts.skipped_secondary_supplementary);
        if skipped_secondary > 0 {
            info!(
                "skipped {skipped_secondary} secondary/supplementary input \
                 records, reads are re-aligned from their primary record"
            );
        }
        let skipped_no_seq = get(&counts.skipped_no_seq);
        if skipped_no_seq > 0 {
            warn!("skipped {skipped_no_seq} input records without a sequence");
        }
        let skipped_hard_clipped = get(&counts.skipped_hard_clipped);
        if skipped_hard_clipped > 0 {
            warn!(
                "skipped {skipped_hard_clipped} hard-clipped primary input \
                 records, these do not contain the full read sequence"
            );
        }
        let failed = get(&counts.failed);
        if failed > 0 {
            warn!(
                "{failed} reads failed to align due to errors, see debug log"
            );
        }
        Ok(())
    }
}

fn build_aligner(
    reference: &PathBuf,
    preset: AlignPreset,
    threads: usize,
) -> anyhow::Result<Aligner<Built>> {
    if !reference.exists() {
        bail!("reference {reference:?} does not exist");
    }
    let mut builder = Aligner::builder()
        .preset(preset.to_minimap2())
        .with_index_threads(threads)
        .with_cigar()
        .with_cigar_clipping();
    // load the whole reference into a single index part so that target ids
    // are global and line up with the @SQ lines of the output header.
    builder.idxopt.batch_size = u64::MAX / 2;
    let aligner = builder
        .with_index(reference, None)
        .map_err(|e| anyhow!("failed to load index for {reference:?}, {e}"))?;
    if aligner.idx_parts.len() > 1 {
        bail!(
            "index {reference:?} has {} parts, multi-part indices are not \
             supported. Re-build the index with a larger batch size (minimap2 \
             -I) or pass the FASTA directly",
            aligner.idx_parts.len()
        );
    }
    Ok(aligner)
}

fn aligner_targets(aligner: &Aligner<Built>) -> Vec<(String, u32)> {
    (0..aligner.n_seq() as usize)
        .filter_map(|i| aligner.get_seq(i))
        .map(|seq| {
            let name = unsafe { CStr::from_ptr(seq.name) }
                .to_string_lossy()
                .into_owned();
            (name, seq.len)
        })
        .collect()
}

/// Make the output header: a fresh @HD and @SQ lines from the index, any
/// read-group, program and comment lines are carried over from the input.
fn build_output_header(
    input: &HeaderView,
    targets: &[(String, u32)],
) -> anyhow::Result<bam::Header> {
    let mut text = String::from("@HD\tVN:1.6\tSO:unsorted\n");
    for (name, len) in targets {
        text.push_str(&format!("@SQ\tSN:{name}\tLN:{len}\n"));
    }
    let input_text = String::from_utf8_lossy(input.as_bytes());
    for line in input_text.lines() {
        if line.starts_with("@RG")
            || line.starts_with("@PG")
            || line.starts_with("@CO")
        {
            text.push_str(line);
            text.push('\n');
        }
    }
    let view = HeaderView::from_bytes(text.as_bytes());
    Ok(bam::Header::from_template(&view))
}

enum Skip {
    SecondaryOrSupplementary,
    NoSequence,
    HardClipped,
}

/// Make an unaligned record in the original read orientation carrying all of
/// the per-read tags from `record`.
fn make_template(record: &bam::Record) -> Result<bam::Record, Skip> {
    if record.is_secondary() || record.is_supplementary() {
        return Err(Skip::SecondaryOrSupplementary);
    }
    if record.seq_len() == 0 {
        return Err(Skip::NoSequence);
    }
    if record.cigar().iter().any(|op| matches!(op, Cigar::HardClip(_))) {
        return Err(Skip::HardClipped);
    }
    let mut seq = record.seq().as_bytes();
    let mut qual = record.qual().to_vec();
    if record.is_reverse() && !record.is_unmapped() {
        seq = revcomp(seq);
        qual.reverse();
    }

    let mut template = bam::Record::new();
    template.set(record.qname(), None, &seq, &qual);
    template.set_tid(-1);
    template.set_pos(-1);
    template.set_mtid(-1);
    template.set_mpos(-1);
    template.set_insert_size(0);
    template.set_mapq(0);
    template.set_flags(
        FLAG_UNMAPPED | (record.flags() & (FLAG_QC_FAIL | FLAG_DUPLICATE)),
    );

    let mut has_mod_tags = false;
    for (tag, value) in record.aux_iter().filter_map(|res| res.ok()) {
        if tag == MN_TAG || ALIGNMENT_TAGS.iter().any(|t| *t == tag) {
            continue;
        }
        has_mod_tags |= MOD_TAGS.iter().any(|t| *t == tag);
        // tags were read from a valid record, push cannot fail on duplicates
        // because we skip none of the input tags twice.
        let _ = template.push_aux(tag, value);
    }
    if has_mod_tags {
        let _ = template.push_aux(MN_TAG, Aux::I32(seq.len() as i32));
    }
    Ok(template)
}

struct MapSettings {
    secondary: bool,
    unmapped: bool,
    md: bool,
}

enum AlignOutcome {
    Mapped(Vec<bam::Record>),
    Unmapped(Option<bam::Record>),
}

fn is_secondary(mapping: &Mapping) -> bool {
    !mapping.is_primary && !mapping.is_supplementary
}

fn to_cigar_string(ops: &[(u32, u8)]) -> anyhow::Result<CigarString> {
    ops.iter()
        .map(|&(len, op)| {
            let cigar = match op {
                0 => Cigar::Match(len),
                1 => Cigar::Ins(len),
                2 => Cigar::Del(len),
                3 => Cigar::RefSkip(len),
                4 => Cigar::SoftClip(len),
                5 => Cigar::HardClip(len),
                6 => Cigar::Pad(len),
                7 => Cigar::Equal(len),
                8 => Cigar::Diff(len),
                _ => bail!("invalid cigar op {op}"),
            };
            Ok(cigar)
        })
        .collect::<anyhow::Result<Vec<Cigar>>>()
        .map(CigarString)
}

fn sa_entry(mapping: &Mapping) -> Option<String> {
    let alignment = mapping.alignment.as_ref()?;
    Some(format!(
        "{},{},{},{},{},{};",
        mapping.target_name.as_ref()?,
        mapping.target_start + 1,
        mapping.strand,
        alignment.cigar_str.as_ref()?,
        mapping.mapq,
        alignment.nm
    ))
}

fn align_read(
    template: bam::Record,
    aligner: &Aligner<Built>,
    settings: &MapSettings,
) -> anyhow::Result<AlignOutcome> {
    let seq = template.seq().as_bytes();
    let qual = template.qual().to_vec();
    let qname = template.qname().to_vec();

    let mappings = aligner
        .map(&seq, false, settings.md, None, None, Some(&qname))
        .map_err(|e| anyhow!("minimap2 failed, {e}"))?
        .into_iter()
        .filter(|m| m.alignment.is_some() && m.target_id >= 0)
        .filter(|m| settings.secondary || !is_secondary(m))
        .collect::<Vec<Mapping>>();

    if mappings.is_empty() {
        let unmapped = if settings.unmapped { Some(template) } else { None };
        return Ok(AlignOutcome::Unmapped(unmapped));
    }

    let rc = if mappings.iter().any(|m| m.strand == Strand::Reverse) {
        let mut rev_qual = qual.clone();
        rev_qual.reverse();
        Some((revcomp(&seq), rev_qual))
    } else {
        None
    };

    // SA tags list the other primary/supplementary alignments of a chimeric
    // read.
    let sa_entries = mappings
        .iter()
        .map(|m| if is_secondary(m) { None } else { sa_entry(m) })
        .collect::<Vec<Option<String>>>();

    let mut records = Vec::with_capacity(mappings.len());
    for (i, mapping) in mappings.iter().enumerate() {
        let alignment = mapping.alignment.as_ref().unwrap();
        let cigar = to_cigar_string(alignment.cigar.as_deref().unwrap_or(&[]))?;
        let (aln_seq, aln_qual) = match (mapping.strand, rc.as_ref()) {
            (Strand::Reverse, Some((rc_seq, rc_qual))) => {
                (rc_seq.as_slice(), rc_qual.as_slice())
            }
            _ => (seq.as_slice(), qual.as_slice()),
        };

        let mut record = template.clone();
        record.set(&qname, Some(&cigar), aln_seq, aln_qual);
        record.set_tid(mapping.target_id);
        record.set_pos(mapping.target_start as i64);
        record.set_mapq(mapping.mapq.min(255) as u8);
        let mut flags = template.flags() & !FLAG_UNMAPPED;
        if mapping.strand == Strand::Reverse {
            flags |= FLAG_REVERSE;
        }
        if mapping.is_supplementary {
            flags |= FLAG_SUPPLEMENTARY;
        } else if is_secondary(mapping) {
            flags |= FLAG_SECONDARY;
        }
        record.set_flags(flags);

        record.push_aux(b"NM", Aux::I32(alignment.nm))?;
        if let Some(score) = alignment.alignment_score {
            record.push_aux(b"AS", Aux::I32(score))?;
        }
        if let Some(md) = alignment.md.as_ref() {
            record.push_aux(b"MD", Aux::String(md))?;
        }
        let tp = if is_secondary(mapping) { b'S' } else { b'P' };
        record.push_aux(b"tp", Aux::Char(tp))?;
        if sa_entries[i].is_some() {
            let sa = sa_entries
                .iter()
                .enumerate()
                .filter(|(j, _)| *j != i)
                .filter_map(|(_, e)| e.as_deref())
                .collect::<String>();
            if !sa.is_empty() {
                record.push_aux(b"SA", Aux::String(&sa))?;
            }
        }
        records.push(record);
    }
    Ok(AlignOutcome::Mapped(records))
}
