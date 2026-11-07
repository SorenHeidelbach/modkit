use std::cmp::Ordering;
use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::{stdout, BufWriter, Cursor, Stdout, Write};
use std::path::{Path, PathBuf};

use anyhow::{anyhow, bail, Context, Result as AnyhowResult};
use bitvec::order::Lsb0;
use bitvec::view::BitView;
use charming::component::{
    Axis, DataZoom, DataZoomType, Feature, Legend, Restore, SaveAsImage, Title,
    Toolbox, ToolboxDataZoom,
};
use charming::element::{
    AxisPointer, AxisPointerType, AxisType, Color, Tooltip, Trigger,
};
use charming::series::Bar;
use charming::{Chart, HtmlRenderer};
use crossbeam_channel::Sender;
use derive_new::new;
use gzp::deflate::Bgzf;
use gzp::par::compress::{ParCompress, ParCompressBuilder};
use indicatif::{MultiProgress, ProgressBar};
use itertools::Itertools;
use log::{debug, info, warn};
use prettytable::format::FormatBuilder;
use prettytable::{row, Table};
use random_color::RandomColor;
use rust_htslib::bam::HeaderView;

use crate::mod_base_code::{
    BaseState, DnaBase, ModCodeRepr, ModifiedBasesOptions, ProbHistogram,
    DNA_BASE_COLORS, MOD_COLORS,
};
use crate::pileup::bedrmod::BedRModArgs;
use crate::pileup::prob_stats::{ProbStats, PROB_STATS_HEADER_FIELDS};
use crate::pileup::{ModBasePileup2, PileupFeatureCounts2};
use crate::sample_probs::PercentileThresholdQual;
use crate::util::{get_ticker_with_rate, TAB};

pub trait PileupWriter<T> {
    fn write(
        &mut self,
        item: T,
        motif_labels: &[String],
    ) -> anyhow::Result<u64>;
}

pub trait OutWriter<T> {
    fn write(&mut self, item: T) -> AnyhowResult<u64>;
}

pub struct BedMethylWriter2<T: Write> {
    buff: Cursor<Vec<u8>>,
    inner: RecordingWriter<T>,
    return_mem: Sender<ModBasePileup2>,
    bedrmod_spec: bool,
}

impl BedMethylWriter2<BufWriter<std::io::Stdout>> {
    pub(crate) fn new_stdout(
        with_header: bool,
        prob_stats: bool,
        bed_rmod_args: &BedRModArgs,
        header: &HeaderView,
        modified_bases_options: Option<&Vec<ModifiedBasesOptions>>,
        multi_progress: MultiProgress,
        return_mem: Sender<ModBasePileup2>,
    ) -> anyhow::Result<Self> {
        let write_pb = multi_progress.add(get_ticker_with_rate());
        write_pb.set_message(format!("B written stdout"));
        write_pb.set_position(0);
        let mut writer = RecordingWriter::new_stdout(write_pb);
        if with_header {
            writer.write(
                bedmethyl_header_with_prob_stats(prob_stats).as_bytes(),
            )?;
        } else if bed_rmod_args.enabled() {
            let modified_bases_options =
                modified_bases_options.ok_or_else(|| {
                    anyhow!("--modified-bases required for --bedrmod")
                })?;
            let bedrmod_header =
                bed_rmod_args.header(&header, modified_bases_options)?;
            writer.write(bedrmod_header.as_bytes())?;
        }

        let buff = Cursor::new(vec![0u8; 1 << 20]);
        Ok(Self {
            buff,
            inner: writer,
            return_mem,
            bedrmod_spec: bed_rmod_args.enabled(),
        })
    }
}

impl BedMethylWriter2<ParCompress<Bgzf>> {
    pub(crate) fn new_bgzf(
        out_path: &PathBuf,
        with_header: bool,
        prob_stats: bool,
        bed_rmod_args: &BedRModArgs,
        header: &HeaderView,
        modified_bases_options: Option<&Vec<ModifiedBasesOptions>>,
        multi_progress: MultiProgress,
        return_mem: Sender<ModBasePileup2>,
        bgzf_threads: usize,
    ) -> anyhow::Result<Self> {
        let fh = File::create(out_path)?;
        let write_pb = multi_progress.add(get_ticker_with_rate());
        let out_fn = out_path
            .file_name()
            .and_then(|x| x.to_str())
            .map(|x| x.to_string())
            .unwrap_or_else(|| "failed to parse filename".to_string());
        write_pb.set_message(format!("B written to output: {out_fn}"));
        write_pb.set_position(0);
        let mut writer = RecordingWriter::new_bgzf(fh, bgzf_threads, write_pb);
        if with_header {
            writer.write(
                bedmethyl_header_with_prob_stats(prob_stats).as_bytes(),
            )?;
        } else if bed_rmod_args.enabled() {
            let modified_bases_options =
                modified_bases_options.ok_or_else(|| {
                    anyhow!("--modified-bases required for --bedrmod")
                })?;
            let bedrmod_header =
                bed_rmod_args.header(&header, modified_bases_options)?;
            writer.write(bedrmod_header.as_bytes())?;
        }
        let buff = Cursor::new(vec![0u8; 1 << 20]);
        Ok(Self {
            buff,
            inner: writer,
            return_mem,
            bedrmod_spec: bed_rmod_args.enabled(),
        })
    }
}

impl BedMethylWriter2<BufWriter<File>> {
    pub(crate) fn new(
        out_path: &PathBuf,
        with_header: bool,
        prob_stats: bool,
        bed_rmod_args: &BedRModArgs,
        header: &HeaderView,
        modified_bases_options: Option<&Vec<ModifiedBasesOptions>>,
        multi_progress: MultiProgress,
        return_mem: Sender<ModBasePileup2>,
    ) -> anyhow::Result<Self> {
        let fh = File::create(out_path)?;
        let write_pb = multi_progress.add(get_ticker_with_rate());
        let out_fn = out_path
            .file_name()
            .and_then(|x| x.to_str())
            .map(|x| x.to_string())
            .unwrap_or_else(|| "failed to parse filename".to_string());
        write_pb.set_message(format!("B written to output: {out_fn}"));
        write_pb.set_position(0);
        let mut writer = RecordingWriter::new_file(fh, write_pb);
        if with_header {
            writer.write(
                bedmethyl_header_with_prob_stats(prob_stats).as_bytes(),
            )?;
        } else if bed_rmod_args.enabled() {
            let modified_bases_options =
                modified_bases_options.ok_or_else(|| {
                    anyhow!("--modified-bases required for --bedrmod")
                })?;
            let bedrmod_header =
                bed_rmod_args.header(&header, modified_bases_options)?;
            writer.write(bedrmod_header.as_bytes())?;
        }

        let buff = Cursor::new(vec![0u8; 1 << 20]);
        Ok(Self {
            buff,
            inner: writer,
            return_mem,
            bedrmod_spec: bed_rmod_args.enabled(),
        })
    }
}

impl<T: Write> PileupWriter<ModBasePileup2> for BedMethylWriter2<T> {
    fn write(
        &mut self,
        item: ModBasePileup2,
        _motif_labels: &[String],
    ) -> anyhow::Result<u64> {
        if item.preformatted {
            self.inner.write(&item.formatted)?;
            let n_rows = item.formatted_rows;
            let _ = self.return_mem.send(item);
            return Ok(n_rows);
        }
        self.buff.set_position(0);
        let n_rows = item.position_feature_counts.len() as u64;
        let chrom_name = &item.chrom_name;
        let prob_stats = (!item.prob_stats.is_empty()).then(|| {
            debug_assert_eq!(
                item.prob_stats.len(),
                item.position_feature_counts.len()
            );
            &item.prob_stats
        });
        for (i, pfc) in item.position_feature_counts.iter().enumerate() {
            format_feature_counts2(
                chrom_name,
                &mut self.buff,
                pfc,
                self.bedrmod_spec,
                prob_stats.map(|x| &x[i]),
            )
            .unwrap();
            let pos = self.buff.position() as usize;
            if pos >= 1 << 20 {
                self.inner.write(&self.buff.get_ref()[..pos]).unwrap();
                self.buff.set_position(0);
            }
        }
        let pos = self.buff.position() as usize;
        self.inner.write(&self.buff.get_ref()[..pos]).unwrap();
        let _ = self.return_mem.send(item);
        Ok(n_rows)
    }
}

pub fn bedmethyl_header() -> String {
    bedmethyl_header_op("valid_coverage")
}

/// bedMethyl header, with the `--prob-stats` columns appended when
/// requested.
pub(crate) fn bedmethyl_header_with_prob_stats(prob_stats: bool) -> String {
    let header = bedmethyl_header();
    if prob_stats {
        let extra = PROB_STATS_HEADER_FIELDS.join("\t");
        format!("{}\t{extra}\n", header.trim_end_matches('\n'))
    } else {
        header
    }
}

pub fn bedrmod_bedmethyl_header() -> String {
    bedmethyl_header_op("coverage")
}

fn bedmethyl_header_op(coverage_field: &'static str) -> String {
    let fields = [
        "chrom",
        "chromStart",
        "chromEnd",
        "name",
        "score",
        "strand",
        "thickStart",
        "thickEnd",
        "color",
        coverage_field,
        "percent_modified",
        "count_modified",
        "count_canonical",
        "count_other_mod",
        "count_delete",
        "count_fail",
        "count_diff",
        "count_nocall",
    ];
    let fields = fields.join("\t");
    format!("#{fields}\n")
}

pub struct MultipleMotifBedmethylWriter<T: Write> {
    writer: RecordingWriter<T>,
    return_mem: Sender<ModBasePileup2>,
    bedrmod_spec: bool,
}

impl MultipleMotifBedmethylWriter<BufWriter<std::io::Stdout>> {
    pub(crate) fn new_stdout(
        with_header: bool,
        prob_stats: bool,
        bed_rmod_args: &BedRModArgs,
        header: &HeaderView,
        modified_bases_options: Option<&Vec<ModifiedBasesOptions>>,
        return_mem: Sender<ModBasePileup2>,
        multi_progress: MultiProgress,
    ) -> anyhow::Result<Self> {
        let mut writer = BufWriter::new(stdout());
        if with_header {
            writer.write(
                bedmethyl_header_with_prob_stats(prob_stats).as_bytes(),
            )?;
        } else if bed_rmod_args.enabled() {
            let modified_bases_options =
                modified_bases_options.ok_or_else(|| {
                    anyhow!("--modified-bases required for --bedrmod")
                })?;
            let bedrmod_header =
                bed_rmod_args.header(&header, modified_bases_options)?;
            writer.write(bedrmod_header.as_bytes())?;
        }

        let write_pb = multi_progress.add(get_ticker_with_rate());
        write_pb.set_message(format!("B written to stdout"));
        write_pb.set_position(0);
        let writer = RecordingWriter { inner: writer, pb: write_pb };

        Ok(Self { writer, return_mem, bedrmod_spec: bed_rmod_args.enabled() })
    }
}

impl MultipleMotifBedmethylWriter<BufWriter<File>> {
    pub(crate) fn new_file(
        out_path: &PathBuf,
        with_header: bool,
        prob_stats: bool,
        bed_rmod_args: &BedRModArgs,
        header: &HeaderView,
        modified_bases_options: Option<&Vec<ModifiedBasesOptions>>,
        multi_progress: MultiProgress,
        return_mem: Sender<ModBasePileup2>,
    ) -> anyhow::Result<Self> {
        let fh = File::create(out_path)?;
        let write_pb = multi_progress.add(get_ticker_with_rate());
        let out_fn = out_path
            .file_name()
            .and_then(|x| x.to_str())
            .map(|x| x.to_string())
            .unwrap_or_else(|| "failed to parse filename".to_string());
        write_pb.set_message(format!("B written to output: {out_fn}"));
        write_pb.set_position(0);

        let mut writer = RecordingWriter::new_file(fh, write_pb);
        if with_header {
            writer.write(
                bedmethyl_header_with_prob_stats(prob_stats).as_bytes(),
            )?;
        } else if bed_rmod_args.enabled() {
            let modified_bases_options =
                modified_bases_options.ok_or_else(|| {
                    anyhow!("--modified-bases required for --bedrmod")
                })?;
            let bedrmod_header =
                bed_rmod_args.header(&header, modified_bases_options)?;
            writer.write(bedrmod_header.as_bytes())?;
        }

        Ok(Self { writer, return_mem, bedrmod_spec: bed_rmod_args.enabled() })
    }
}

impl MultipleMotifBedmethylWriter<ParCompress<Bgzf>> {
    pub(crate) fn new_bgzf(
        out_path: &PathBuf,
        with_header: bool,
        prob_stats: bool,
        bed_rmod_args: &BedRModArgs,
        header: &HeaderView,
        modified_bases_options: Option<&Vec<ModifiedBasesOptions>>,
        multi_progress: MultiProgress,
        return_mem: Sender<ModBasePileup2>,
        bgzf_threads: usize,
    ) -> anyhow::Result<Self> {
        let fh = File::create(out_path)?;
        let write_pb = multi_progress.add(get_ticker_with_rate());
        let out_fn = out_path
            .file_name()
            .and_then(|x| x.to_str())
            .map(|x| x.to_string())
            .unwrap_or_else(|| "failed to parse filename".to_string());
        write_pb.set_message(format!("B written to output: {out_fn}"));
        write_pb.set_position(0);

        let mut writer = RecordingWriter::new_bgzf(fh, bgzf_threads, write_pb);
        if with_header {
            writer.write(
                bedmethyl_header_with_prob_stats(prob_stats).as_bytes(),
            )?;
        } else if bed_rmod_args.enabled() {
            let modified_bases_options =
                modified_bases_options.ok_or_else(|| {
                    anyhow!("--modified-bases required for --bedrmod")
                })?;
            let bedrmod_header =
                bed_rmod_args.header(&header, modified_bases_options)?;
            writer.write(bedrmod_header.as_bytes())?;
        }

        Ok(Self { writer, return_mem, bedrmod_spec: bed_rmod_args.enabled() })
    }
}

impl<T: Write> PileupWriter<ModBasePileup2>
    for MultipleMotifBedmethylWriter<T>
{
    fn write(
        &mut self,
        item: ModBasePileup2,
        motif_labels: &[String],
    ) -> anyhow::Result<u64> {
        if item.preformatted {
            self.writer.write(&item.formatted)?;
            let n_rows = item.formatted_rows;
            let _ = self.return_mem.send(item);
            return Ok(n_rows);
        }
        let mut rows_written = 0u64;
        let chrom_name = &item.chrom_name;
        let mut buff = Cursor::new([0u8; 1024]);
        for (i, feature_count) in
            item.position_feature_counts.iter().enumerate()
        {
            let prob_stats = item.prob_stats.get(i);
            let idxs = feature_count.motif_idxs.view_bits::<Lsb0>().iter_ones();
            for idx in idxs {
                let motif_label = &motif_labels[idx];
                let pos = feature_count.position;
                let pp1 = pos + 1;
                let fraction_modified = feature_count.n_modified as f32
                    / feature_count.filtered_coverage as f32;
                write!(&mut buff, "{}{TAB}", chrom_name)?;
                write!(&mut buff, "{pos}{TAB}{pp1}{TAB}")?;
                write!(
                    &mut buff,
                    "{},{motif_label}{TAB}",
                    feature_count.mod_code
                )?;
                write!(&mut buff, "{}{TAB}", feature_count.filtered_coverage)?;
                write!(&mut buff, "{}{TAB}", feature_count.raw_strand)?;
                write!(&mut buff, "{pos}{TAB}{pp1}{TAB}")?;
                write!(&mut buff, "255,0,0{TAB}")?;
                let cov_val = if self.bedrmod_spec {
                    feature_count.filtered_coverage + feature_count.n_filtered
                } else {
                    feature_count.filtered_coverage
                };
                write!(&mut buff, "{}{TAB}", cov_val)?;
                write!(&mut buff, "{:.2}{TAB}", fraction_modified * 100f32)?;
                write!(&mut buff, "{}{TAB}", feature_count.n_modified)?;
                write!(&mut buff, "{}{TAB}", feature_count.n_canonical)?;
                write!(&mut buff, "{}{TAB}", feature_count.n_other_modified)?;
                write!(&mut buff, "{}{TAB}", feature_count.n_delete)?;
                write!(&mut buff, "{}{TAB}", feature_count.n_filtered)?;
                write!(&mut buff, "{}{TAB}", feature_count.n_diff)?;
                write!(&mut buff, "{}", feature_count.n_nocall)?;
                if let Some(prob_stats) = prob_stats {
                    prob_stats.write_columns(&mut buff)?;
                }
                write!(&mut buff, "\n")?;
                let pos = buff.position() as usize;

                self.writer
                    .write(&buff.get_ref()[..pos])
                    .with_context(|| "failed to write row")?;
                buff.set_position(0);
                rows_written = rows_written.saturating_add(1);
            }
        }
        let _ = self.return_mem.send(item);
        Ok(rows_written)
    }
}

fn format_feature_counts2(
    chrom_name: &String,
    buff: &mut Cursor<Vec<u8>>,
    feature_count: &PileupFeatureCounts2,
    bedrmod_spec: bool,
    prob_stats: Option<&ProbStats>,
) -> anyhow::Result<()> {
    let pos = feature_count.position;
    let pp1 = pos + 1;
    let fraction_modified = feature_count.n_modified as f32
        / feature_count.filtered_coverage as f32;
    write!(buff, "{}{TAB}", chrom_name)?;
    write!(buff, "{pos}{TAB}{pp1}{TAB}")?;
    write!(buff, "{}{TAB}", feature_count.mod_code)?;
    write!(buff, "{}{TAB}", feature_count.filtered_coverage)?;
    write!(buff, "{}{TAB}", feature_count.raw_strand)?;
    write!(buff, "{pos}{TAB}{pp1}{TAB}")?;
    write!(buff, "255,0,0{TAB}")?;
    let cov_val = if bedrmod_spec {
        feature_count.filtered_coverage + feature_count.n_filtered
    } else {
        feature_count.filtered_coverage
    };
    write!(buff, "{}{TAB}", cov_val)?;
    write!(buff, "{:.2}{TAB}", fraction_modified * 100f32)?;
    write!(buff, "{}{TAB}", feature_count.n_modified)?;
    write!(buff, "{}{TAB}", feature_count.n_canonical)?;
    write!(buff, "{}{TAB}", feature_count.n_other_modified)?;
    write!(buff, "{}{TAB}", feature_count.n_delete)?;
    write!(buff, "{}{TAB}", feature_count.n_filtered)?;
    write!(buff, "{}{TAB}", feature_count.n_diff)?;
    write!(buff, "{}", feature_count.n_nocall)?;
    if let Some(prob_stats) = prob_stats {
        prob_stats.write_columns(buff)?;
    }
    write!(buff, "\n")?;
    Ok(())
}

pub(crate) struct RecordingWriter<T: Write> {
    inner: T,
    pb: ProgressBar,
}

impl RecordingWriter<BufWriter<std::io::Stdout>> {
    pub(crate) fn new_stdout(pb: ProgressBar) -> Self {
        Self { inner: BufWriter::with_capacity(1 << 20, std::io::stdout()), pb }
    }
}

impl RecordingWriter<BufWriter<File>> {
    pub(crate) fn new_file(file: File, pb: ProgressBar) -> Self {
        Self { inner: BufWriter::with_capacity(1 << 20, file), pb }
    }
}

impl RecordingWriter<ParCompress<Bgzf>> {
    pub(crate) fn new_bgzf(
        file: File,
        compress_threads: usize,
        pb: ProgressBar,
    ) -> Self {
        let inner = ParCompressBuilder::<Bgzf>::new()
            .num_threads(compress_threads)
            .unwrap()
            .from_writer(file);
        Self { inner, pb }
    }
}
impl<T: Write> RecordingWriter<T> {
    pub(crate) fn write(&mut self, bulk: &[u8]) -> anyhow::Result<()> {
        let n = bulk.len();
        self.inner.write_all(bulk)?;
        self.pb.inc(n as u64);
        self.flush()?;
        Ok(())
    }

    pub(crate) fn flush(&mut self) -> anyhow::Result<()> {
        self.inner.flush()?;
        Ok(())
    }
}

impl<T: Write> Drop for RecordingWriter<T> {
    fn drop(&mut self) {
        self.pb.finish_and_clear();
        let _ = self.inner.flush();
    }
}

pub struct PhasedBedMethylWriter<T: Write> {
    hp1_writer: RecordingWriter<T>,
    hp2_writer: RecordingWriter<T>,
    combined_writer: RecordingWriter<T>,
    return_mem: Sender<ModBasePileup2>,
}

impl PhasedBedMethylWriter<BufWriter<File>> {
    fn make_writer(
        out_dir: &PathBuf,
        prefix: &str,
        hp_label: &str,
        force: bool,
        multi_progress: MultiProgress,
    ) -> anyhow::Result<RecordingWriter<BufWriter<File>>> {
        let filename = format!("{prefix}{hp_label}.bedmethyl");
        let path = out_dir.join(Path::new(&filename));
        let file =
            if force { File::create(path) } else { File::create_new(path) }?;
        let write_pb = multi_progress.add(get_ticker_with_rate());
        write_pb.set_message(format!("B written: {hp_label}"));
        write_pb.set_position(0);
        let writer = RecordingWriter::new_file(file, write_pb);
        Ok(writer)
    }

    pub fn new_file(
        out_dir: &PathBuf,
        prefix: Option<&String>,
        force: bool,
        multi_progress: MultiProgress,
        return_mem: Sender<ModBasePileup2>,
    ) -> anyhow::Result<Self> {
        let prefix = prefix.map(|x| format!("{x}_")).unwrap_or("".to_string());
        let hp1_writer = Self::make_writer(
            &out_dir,
            &prefix,
            "hp1",
            force,
            multi_progress.clone(),
        )?;
        let hp2_writer = Self::make_writer(
            &out_dir,
            &prefix,
            "hp2",
            force,
            multi_progress.clone(),
        )?;
        let combined_writer = Self::make_writer(
            &out_dir,
            &prefix,
            "combined",
            force,
            multi_progress.clone(),
        )?;
        Ok(Self { hp1_writer, hp2_writer, combined_writer, return_mem })
    }
}

impl PhasedBedMethylWriter<ParCompress<Bgzf>> {
    fn make_writer(
        out_dir: &PathBuf,
        prefix: &str,
        hp_label: &str,
        force: bool,
        compression_threads: usize,
        multi_progress: MultiProgress,
    ) -> anyhow::Result<RecordingWriter<ParCompress<Bgzf>>> {
        let filename = format!("{prefix}{hp_label}.bed.gz");
        let path = out_dir.join(Path::new(&filename));
        let file =
            if force { File::create(path) } else { File::create_new(path) }?;
        let write_pb = multi_progress.add(get_ticker_with_rate());
        write_pb.set_message(format!("B written: {hp_label}"));
        write_pb.set_position(0);
        let writer =
            RecordingWriter::new_bgzf(file, compression_threads, write_pb);
        Ok(writer)
    }

    pub fn new_bgzf(
        out_dir: &PathBuf,
        prefix: Option<&String>,
        force: bool,
        multi_progress: MultiProgress,
        return_mem: Sender<ModBasePileup2>,
        compression_threads: usize,
    ) -> anyhow::Result<Self> {
        let prefix = prefix.map(|x| format!("{x}_")).unwrap_or("".to_string());
        let hp1_writer = Self::make_writer(
            &out_dir,
            &prefix,
            "hp1",
            force,
            compression_threads,
            multi_progress.clone(),
        )?;
        let hp2_writer = Self::make_writer(
            &out_dir,
            &prefix,
            "hp2",
            force,
            compression_threads,
            multi_progress.clone(),
        )?;
        let combined_writer = Self::make_writer(
            &out_dir,
            &prefix,
            "combined",
            force,
            compression_threads,
            multi_progress.clone(),
        )?;
        Ok(Self { hp1_writer, hp2_writer, combined_writer, return_mem })
    }
}

impl<T> PileupWriter<ModBasePileup2> for PhasedBedMethylWriter<T>
where
    T: Write + Send,
{
    fn write(
        &mut self,
        item: ModBasePileup2,
        _motif_labels: &[String],
    ) -> anyhow::Result<u64> {
        let chrom_name = &item.chrom_name;
        let combined_counts = &item.position_feature_counts;
        let [hp1, hp2] = &item.phased_feature_counts;
        let total_rows = combined_counts.len() + hp1.len() + hp2.len();

        // TODO: make the "buff"s part of the object.
        std::thread::scope(|scope| {
            let hp1_handle = scope.spawn(|| {
                let mut buff = Cursor::new(vec![0u8; 1 << 20]);
                for pfc in hp1.iter().filter(|x| x.is_valid()) {
                    format_feature_counts2(
                        chrom_name, &mut buff, pfc, false, None,
                    )
                    .unwrap();
                    let pos = buff.position() as usize;
                    if pos >= 1 << 20 {
                        self.hp1_writer.write(&buff.get_ref()[..pos]).unwrap();
                        buff.set_position(0);
                    }
                }
                let pos = buff.position() as usize;
                self.hp1_writer.write(&buff.get_ref()[..pos]).unwrap();
            });
            let hp2_handle = scope.spawn(|| {
                let mut buff = Cursor::new(vec![0u8; 1 << 20]);
                for pfc in hp2.iter().filter(|x| x.is_valid()) {
                    format_feature_counts2(
                        chrom_name, &mut buff, pfc, false, None,
                    )
                    .unwrap();
                    let pos = buff.position() as usize;
                    if pos >= 1 << 20 {
                        self.hp2_writer.write(&buff.get_ref()[..pos]).unwrap();
                        buff.set_position(0);
                    }
                }
                let pos = buff.position() as usize;
                self.hp2_writer.write(&buff.get_ref()[..pos]).unwrap();
            });
            let combined_handle = scope.spawn(|| {
                let mut buff = Cursor::new(vec![0u8; 1 << 20]);
                for pfc in combined_counts.iter().filter(|x| x.is_valid()) {
                    format_feature_counts2(
                        &chrom_name,
                        &mut buff,
                        pfc,
                        false,
                        None,
                    )
                    .unwrap();
                    let pos = buff.position() as usize;
                    if pos >= 1 << 20 {
                        self.combined_writer
                            .write(&buff.get_ref()[..pos])
                            .unwrap();
                        buff.set_position(0);
                    }
                }
                let pos = buff.position() as usize;
                self.combined_writer.write(&buff.get_ref()[..pos]).unwrap();
            });
            let _ = hp1_handle.join().unwrap();
            let _ = hp2_handle.join().unwrap();
            let _ = combined_handle.join().unwrap();
        });
        let _ = self.return_mem.send(item);

        Ok(total_rows as u64)
    }
}
