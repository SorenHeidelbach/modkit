//! Fast bedMethyl row formatting.
//!
//! Rows are formatted by the pileup worker threads (in parallel) into byte
//! buffers, the writer thread then only has to copy bytes to the output.
//! The output is identical to the `write!` based formatting used before,
//! see the tests at the bottom of this module.

use bitvec::order::Lsb0;
use bitvec::view::BitView;

use crate::mod_base_code::ModCodeRepr;
use crate::pileup::prob_stats::ProbStats;
use crate::pileup::{ModBasePileup2, PileupFeatureCounts2};

/// How the rows of a [`ModBasePileup2`] are written.
#[derive(Debug, Clone)]
pub(crate) struct RowFormatter {
    /// bedRMod: column 10 is the total coverage (valid + failed).
    pub(crate) bedrmod_spec: bool,
    /// When set, one row is written per motif matched at the position, with
    /// the name column `<mod_code>,<motif_label>`.
    pub(crate) motif_labels: Option<Vec<String>>,
}

impl RowFormatter {
    /// Formats all rows of `item` into `item.formatted` and records the
    /// number of rows.
    pub(crate) fn format_into(&self, item: &mut ModBasePileup2) {
        let mut buf = std::mem::take(&mut item.formatted);
        buf.clear();
        let mut n_rows = 0u64;
        let chrom = item.chrom_name.as_bytes();
        for (i, row) in item.position_feature_counts.iter().enumerate() {
            let prob_stats = item.prob_stats.get(i);
            match &self.motif_labels {
                None => {
                    push_row(
                        &mut buf,
                        chrom,
                        row,
                        None,
                        self.bedrmod_spec,
                        prob_stats,
                    );
                    n_rows += 1;
                }
                Some(labels) => {
                    for idx in row.motif_idxs.view_bits::<Lsb0>().iter_ones() {
                        push_row(
                            &mut buf,
                            chrom,
                            row,
                            Some(labels[idx].as_bytes()),
                            self.bedrmod_spec,
                            prob_stats,
                        );
                        n_rows += 1;
                    }
                }
            }
        }
        item.formatted = buf;
        item.formatted_rows = n_rows;
        item.preformatted = true;
    }
}

#[inline]
pub(crate) fn push_uint(buf: &mut Vec<u8>, mut v: u64) {
    let mut tmp = [0u8; 20];
    let mut i = tmp.len();
    loop {
        i -= 1;
        tmp[i] = b'0' + (v % 10) as u8;
        v /= 10;
        if v == 0 {
            break;
        }
    }
    buf.extend_from_slice(&tmp[i..]);
}

#[inline]
fn push_char(buf: &mut Vec<u8>, c: char) {
    let mut tmp = [0u8; 4];
    buf.extend_from_slice(c.encode_utf8(&mut tmp).as_bytes());
}

#[inline]
fn push_mod_code(buf: &mut Vec<u8>, code: ModCodeRepr) {
    match code {
        ModCodeRepr::Code(c) => push_char(buf, c),
        ModCodeRepr::ChEbi(x) => push_uint(buf, x as u64),
    }
}

/// Same text as `format!("{:.N}", v)` for a finite, non-negative `f32`
/// where N = log10(scale). Any `f32` times 10^N (N <= 4) is exact in an
/// `f64` (24 + 14 bits of mantissa), so rounding the product with ties to
/// even matches the exact decimal rounding done by `{:.N}`.
#[inline]
pub(crate) fn push_fixed<const DECIMALS: u32>(buf: &mut Vec<u8>, v: f32) {
    debug_assert!(DECIMALS <= 4);
    if !(v.is_finite() && v >= 0f32 && v < 1e9) {
        use std::io::Write;
        let _ = write!(buf, "{:.*}", DECIMALS as usize, v);
        return;
    }
    let scale = 10u64.pow(DECIMALS);
    let x = (v as f64 * scale as f64).round_ties_even() as u64;
    push_uint(buf, x / scale);
    buf.push(b'.');
    let frac = x % scale;
    // left-pad the fractional digits with zeros
    let mut div = scale / 10;
    while div > 0 {
        buf.push(b'0' + ((frac / div) % 10) as u8);
        div /= 10;
    }
}

#[inline]
fn push_row(
    buf: &mut Vec<u8>,
    chrom: &[u8],
    row: &PileupFeatureCounts2,
    motif_label: Option<&[u8]>,
    bedrmod_spec: bool,
    prob_stats: Option<&ProbStats>,
) {
    let pos = row.position as u64;
    let pp1 = row.position.wrapping_add(1) as u64;
    // computed exactly as before, in f32
    let fraction_modified =
        row.n_modified as f32 / row.filtered_coverage as f32;
    buf.extend_from_slice(chrom);
    buf.push(b'\t');
    push_uint(buf, pos);
    buf.push(b'\t');
    push_uint(buf, pp1);
    buf.push(b'\t');
    push_mod_code(buf, row.mod_code);
    if let Some(label) = motif_label {
        buf.push(b',');
        buf.extend_from_slice(label);
    }
    buf.push(b'\t');
    push_uint(buf, row.filtered_coverage as u64);
    buf.push(b'\t');
    push_char(buf, row.raw_strand);
    buf.push(b'\t');
    push_uint(buf, pos);
    buf.push(b'\t');
    push_uint(buf, pp1);
    buf.extend_from_slice(b"\t255,0,0\t");
    let cov_val = if bedrmod_spec {
        // u16 arithmetic as before (release builds wrap on overflow)
        row.filtered_coverage.wrapping_add(row.n_filtered)
    } else {
        row.filtered_coverage
    };
    push_uint(buf, cov_val as u64);
    buf.push(b'\t');
    push_fixed::<2>(buf, fraction_modified * 100f32);
    for v in [
        row.n_modified,
        row.n_canonical,
        row.n_other_modified,
        row.n_delete,
        row.n_filtered,
        row.n_diff,
    ] {
        buf.push(b'\t');
        push_uint(buf, v as u64);
    }
    buf.push(b'\t');
    push_uint(buf, row.n_nocall as u64);
    if let Some(prob_stats) = prob_stats {
        prob_stats.push_columns(buf);
    }
    buf.push(b'\n');
}

#[cfg(test)]
mod format_tests {
    use std::io::Write;

    use super::*;
    use crate::mod_base_code::DnaBase;
    use crate::pileup::prob_stats::ProbStats;

    const TAB: char = '\t';

    /// The `write!` based implementation this module replaces.
    fn reference_row(
        chrom_name: &str,
        feature_count: &PileupFeatureCounts2,
        motif_label: Option<&str>,
        bedrmod_spec: bool,
        prob_stats: Option<&ProbStats>,
    ) -> Vec<u8> {
        let mut buff = Vec::new();
        let pos = feature_count.position;
        let pp1 = pos.wrapping_add(1);
        let fraction_modified = feature_count.n_modified as f32
            / feature_count.filtered_coverage as f32;
        write!(buff, "{}{TAB}", chrom_name).unwrap();
        write!(buff, "{pos}{TAB}{pp1}{TAB}").unwrap();
        match motif_label {
            Some(label) => {
                write!(buff, "{},{label}{TAB}", feature_count.mod_code)
            }
            None => write!(buff, "{}{TAB}", feature_count.mod_code),
        }
        .unwrap();
        write!(buff, "{}{TAB}", feature_count.filtered_coverage).unwrap();
        write!(buff, "{}{TAB}", feature_count.raw_strand).unwrap();
        write!(buff, "{pos}{TAB}{pp1}{TAB}").unwrap();
        write!(buff, "255,0,0{TAB}").unwrap();
        let cov_val = if bedrmod_spec {
            feature_count
                .filtered_coverage
                .wrapping_add(feature_count.n_filtered)
        } else {
            feature_count.filtered_coverage
        };
        write!(buff, "{}{TAB}", cov_val).unwrap();
        write!(buff, "{:.2}{TAB}", fraction_modified * 100f32).unwrap();
        write!(buff, "{}{TAB}", feature_count.n_modified).unwrap();
        write!(buff, "{}{TAB}", feature_count.n_canonical).unwrap();
        write!(buff, "{}{TAB}", feature_count.n_other_modified).unwrap();
        write!(buff, "{}{TAB}", feature_count.n_delete).unwrap();
        write!(buff, "{}{TAB}", feature_count.n_filtered).unwrap();
        write!(buff, "{}{TAB}", feature_count.n_diff).unwrap();
        write!(buff, "{}", feature_count.n_nocall).unwrap();
        if let Some(prob_stats) = prob_stats {
            prob_stats.write_columns(&mut buff).unwrap();
        }
        write!(buff, "\n").unwrap();
        buff
    }

    struct XorShift(u64);
    impl XorShift {
        fn next(&mut self) -> u64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            self.0
        }
        fn u16(&mut self) -> u16 {
            match self.next() % 4 {
                0 => (self.next() % 4) as u16,
                1 => (self.next() % 100) as u16,
                2 => u16::MAX - (self.next() % 3) as u16,
                _ => self.next() as u16,
            }
        }
    }

    #[test]
    fn test_percent_matches_std_exhaustive_small() {
        for cov in 1u16..=2000 {
            for n_mod in 0..=cov {
                let v = (n_mod as f32 / cov as f32) * 100f32;
                let mut buf = Vec::new();
                push_fixed::<2>(&mut buf, v);
                assert_eq!(
                    String::from_utf8(buf).unwrap(),
                    format!("{v:.2}"),
                    "{n_mod}/{cov}"
                );
            }
        }
    }

    #[test]
    fn test_push_uint() {
        for v in [0u64, 1, 9, 10, 99, 100, 65535, 4294967295, u64::MAX] {
            let mut buf = Vec::new();
            push_uint(&mut buf, v);
            assert_eq!(String::from_utf8(buf).unwrap(), v.to_string());
        }
    }

    #[test]
    fn test_rows_match_reference_formatting() {
        let mut rng = XorShift(0x9E3779B97F4A7C15);
        let codes = [
            ModCodeRepr::Code('m'),
            ModCodeRepr::Code('h'),
            ModCodeRepr::Code('a'),
            ModCodeRepr::Code('C'),
            ModCodeRepr::ChEbi(21839),
            ModCodeRepr::ChEbi(76792),
        ];
        let labels = ["CG,0", "GATC,1", "CCWGG,1"];
        for i in 0..300_000u64 {
            let filtered_coverage = rng.u16();
            let row = PileupFeatureCounts2::new(
                (rng.next() % (u32::MAX as u64 + 1)) as u32,
                ['+', '-', '.'][(rng.next() % 3) as usize],
                filtered_coverage,
                codes[(rng.next() % codes.len() as u64) as usize],
                rng.u16(),
                if filtered_coverage == 0 {
                    0
                } else {
                    (rng.next() % (filtered_coverage as u64 + 1)) as u16
                },
                rng.u16(),
                rng.u16(),
                rng.u16(),
                rng.u16(),
                rng.u16(),
                (rng.next() % 8) as u8,
            );
            let bedrmod = i % 3 == 0;
            let prob_stats = if i % 2 == 0 {
                let mut values = (0..(rng.next() % 5))
                    .map(|_| (rng.next() % 1000) as f32 / 1000f32)
                    .collect::<Vec<f32>>();
                Some(ProbStats::from_values_for_test(&mut values))
            } else {
                None
            };
            let label = if i % 5 == 0 {
                Some(labels[(rng.next() % 3) as usize])
            } else {
                None
            };
            let mut fast = Vec::new();
            push_row(
                &mut fast,
                b"contig_1",
                &row,
                label.map(|l| l.as_bytes()),
                bedrmod,
                prob_stats.as_ref(),
            );
            let expected = reference_row(
                "contig_1",
                &row,
                label,
                bedrmod,
                prob_stats.as_ref(),
            );
            assert_eq!(
                String::from_utf8(fast).unwrap(),
                String::from_utf8(expected).unwrap(),
                "{row:?}"
            );
        }
    }

    #[test]
    fn test_format_into_multiple_motifs() {
        let mut item = ModBasePileup2::new_empty();
        item.chrom_name = "chr1".to_string();
        let row = PileupFeatureCounts2::new(
            10,
            '+',
            4,
            ModCodeRepr::Code('m'),
            1,
            3,
            0,
            0,
            0,
            0,
            0,
            0b101,
        );
        item.position_feature_counts = vec![row];
        let formatter = RowFormatter {
            bedrmod_spec: false,
            motif_labels: Some(vec![
                "CG,0".to_string(),
                "X".to_string(),
                "GATC,1".to_string(),
            ]),
        };
        formatter.format_into(&mut item);
        assert_eq!(item.formatted_rows, 2);
        let text = String::from_utf8(item.formatted.clone()).unwrap();
        let lines = text.lines().collect::<Vec<_>>();
        assert_eq!(lines.len(), 2);
        assert!(lines[0].starts_with("chr1\t10\t11\tm,CG,0\t4\t+"));
        assert!(lines[1].starts_with("chr1\t10\t11\tm,GATC,1\t4\t+"));
        let _ = DnaBase::C;
    }
}
