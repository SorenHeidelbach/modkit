//! Per-position base modification probability statistics for `pileup
//! --prob-stats`.
//!
//! The counts in a bedMethyl record are made after thresholding each read's
//! base modification probabilities into a call. This module keeps the
//! probabilities themselves, for every read that contributes a canonical,
//! modified, or filtered (failed threshold) call at a position, and
//! summarizes them per output record.
//!
//! For a record with modification code `m`, the per-read value is P(m) as
//! reported in the ML tag (implicitly unmodified positions contribute 0).
//! When modification codes are combined (`--combine-mods`) the per-read value
//! is the summed probability of all modifications of the primary base. The
//! "confidence" of a read is the probability of the most likely class
//! (canonical or any modification), where P(canonical) = 1 - sum(P(mod)).

use std::io::Write;

use crate::mod_base_code::{DnaBase, ModCodeRepr};
use crate::pileup::PileupFeatureCounts2;

/// Column names appended to the bedMethyl header.
pub(crate) const PROB_STATS_HEADER_FIELDS: [&str; 7] = [
    "n_prob",
    "mean_prob",
    "sd_prob",
    "q25_prob",
    "median_prob",
    "q75_prob",
    "mean_confidence",
];

/// Summary of the probabilities for one bedMethyl record.
#[derive(Debug, Copy, Clone, PartialEq, Default)]
pub struct ProbStats {
    /// Number of reads that contributed a probability.
    pub n: u32,
    pub mean: f32,
    /// Population standard deviation.
    pub sd: f32,
    pub q25: f32,
    pub median: f32,
    pub q75: f32,
    /// Mean probability of the most likely class (canonical or modified).
    pub mean_confidence: f32,
}

impl ProbStats {
    /// Compute statistics from per-read probabilities, `values` is sorted
    /// in place. Quantiles use linear interpolation between order statistics
    /// (the default method in numpy and R type 7).
    fn from_values(values: &mut [f32], confidences: &[f32]) -> Self {
        let n = values.len();
        if n == 0 {
            return Self::default();
        }
        let nf = n as f64;
        let mean = values.iter().map(|&x| x as f64).sum::<f64>() / nf;
        let var = values
            .iter()
            .map(|&x| {
                let d = x as f64 - mean;
                d * d
            })
            .sum::<f64>()
            / nf;
        values.sort_unstable_by(|a, b| a.total_cmp(b));
        let mean_confidence = if confidences.is_empty() {
            0f64
        } else {
            confidences.iter().map(|&x| x as f64).sum::<f64>()
                / confidences.len() as f64
        };
        Self {
            n: n as u32,
            mean: mean as f32,
            sd: var.sqrt() as f32,
            q25: quantile_sorted(values, 0.25),
            median: quantile_sorted(values, 0.5),
            q75: quantile_sorted(values, 0.75),
            mean_confidence: mean_confidence as f32,
        }
    }

    /// Writes the extra columns, each preceded by a tab. Records without any
    /// probabilities get a count of 0 and '.' for the statistics.
    pub(crate) fn write_columns<W: Write>(
        &self,
        w: &mut W,
    ) -> std::io::Result<()> {
        if self.n == 0 {
            write!(w, "\t0\t.\t.\t.\t.\t.\t.")
        } else {
            write!(w, "\t{}", self.n)?;
            for v in [
                self.mean,
                self.sd,
                self.q25,
                self.median,
                self.q75,
                self.mean_confidence,
            ] {
                w.write_all(b"\t")?;
                write_fixed4(w, v)?;
            }
            Ok(())
        }
    }
}

impl ProbStats {
    /// Appends the extra columns, each preceded by a tab, same text as
    /// `write_columns`.
    #[inline]
    pub(crate) fn push_columns(&self, buf: &mut Vec<u8>) {
        use crate::pileup::format::{push_fixed, push_uint};
        if self.n == 0 {
            buf.extend_from_slice(b"\t0\t.\t.\t.\t.\t.\t.");
            return;
        }
        buf.push(b'\t');
        push_uint(buf, self.n as u64);
        for v in [
            self.mean,
            self.sd,
            self.q25,
            self.median,
            self.q75,
            self.mean_confidence,
        ] {
            buf.push(b'\t');
            push_fixed::<4>(buf, v);
        }
    }

    #[cfg(test)]
    pub(crate) fn from_values_for_test(values: &mut [f32]) -> Self {
        let confidences = values.iter().map(|x| x.max(1f32 - x)).collect::<Vec<_>>();
        Self::from_values(values, &confidences)
    }
}

/// Writes a finite, non-negative `f32` with 4 decimal places, producing the
/// same text as `format!("{:.4}", v)`. Any `f32` multiplied by 10^4 is exact
/// in an `f64` (24 + 14 bits of mantissa), so rounding the product with ties
/// to even matches the exact decimal rounding used by `{:.4}`. This avoids
/// the slow exact float formatting path, which is costly in the (serial)
/// writer.
#[inline]
fn write_fixed4<W: Write>(w: &mut W, v: f32) -> std::io::Result<()> {
    if !(v.is_finite() && v >= 0f32 && v < 1e9) {
        return write!(w, "{:.4}", v);
    }
    let x = (v as f64 * 10_000f64).round_ties_even() as u64;
    write!(w, "{}.{:04}", x / 10_000, x % 10_000)
}

fn quantile_sorted(sorted: &[f32], q: f64) -> f32 {
    debug_assert!(!sorted.is_empty());
    let h = (sorted.len() - 1) as f64 * q;
    let lo = h.floor() as usize;
    let hi = h.ceil() as usize;
    let lo_v = sorted[lo] as f64;
    let hi_v = sorted[hi] as f64;
    (lo_v + (h - lo as f64) * (hi_v - lo_v)) as f32
}

/// One read's call at one (position, strand). Its per-code probabilities
/// are `values[start..start + n]`.
#[derive(Debug, Clone, Copy)]
struct Obs {
    idx: u32,
    start: u32,
    base: u8,
    n: u8,
}

/// Accumulates per-read probabilities over an interval. Observations are
/// appended to flat arrays, and grouped by output record index with a
/// counting sort when statistics are computed. Records are indexed as
/// `relative_position * 2 + strand`, where the strand index is 0 for '+'
/// and for combined strands ('.'), and 1 for '-'.
#[derive(Debug, Default)]
pub(crate) struct ProbAccumulator {
    width: usize,
    obs: Vec<Obs>,
    /// (interned modification code, probability)
    values: Vec<(u8, f32)>,
    codes: Vec<ModCodeRepr>,
    // scratch space used when computing statistics
    group_starts: Vec<u32>,
    order: Vec<u32>,
    scratch_values: Vec<f32>,
    scratch_conf: Vec<f32>,
}

impl ProbAccumulator {
    pub(crate) fn reset(&mut self, width: usize) {
        self.width = width;
        self.obs.clear();
        self.values.clear();
    }

    #[inline]
    fn code_id(&mut self, code: ModCodeRepr) -> u8 {
        match self.codes.iter().position(|c| *c == code) {
            Some(i) => i as u8,
            None => {
                // more than 255 distinct codes is not realistic
                assert!(self.codes.len() < u8::MAX as usize);
                self.codes.push(code);
                (self.codes.len() - 1) as u8
            }
        }
    }

    /// Record one read's probabilities at a position. `mod_probs` are the
    /// probabilities of each modification code reported for `primary_base`
    /// by this read.
    #[inline]
    pub(crate) fn add(
        &mut self,
        relative_position: u32,
        negative_strand: bool,
        primary_base: DnaBase,
        mod_probs: &[(ModCodeRepr, f32)],
    ) {
        let idx = relative_position as usize * 2 + negative_strand as usize;
        if idx >= self.width * 2 {
            // can only happen if the caller's position mapping is off,
            // drop rather than panic
            debug_assert!(false, "position out of range {idx}");
            return;
        }
        let start = self.values.len() as u32;
        let n = mod_probs.len().min(u8::MAX as usize);
        for &(code, p) in &mod_probs[..n] {
            let code_id = self.code_id(code);
            self.values.push((code_id, p));
        }
        self.obs.push(Obs {
            idx: idx as u32,
            start,
            base: primary_base as u8,
            n: n as u8,
        });
    }

    /// Compute the statistics for each record, `out` is aligned with `rows`.
    pub(crate) fn stats_for_rows(
        &mut self,
        start_pos: u32,
        rows: &[PileupFeatureCounts2],
        combine_mods: bool,
        out: &mut Vec<ProbStats>,
    ) {
        out.clear();
        out.reserve(rows.len());

        // group observations by record index (counting sort)
        let n_idx = self.width * 2;
        self.group_starts.clear();
        self.group_starts.resize(n_idx + 1, 0);
        for o in &self.obs {
            self.group_starts[o.idx as usize + 1] += 1;
        }
        for i in 0..n_idx {
            self.group_starts[i + 1] += self.group_starts[i];
        }
        self.order.clear();
        self.order.resize(self.obs.len(), 0);
        // scatter using the group starts as cursors, afterwards
        // group_starts[i] is the end of group i, shift back by one
        for (i, o) in self.obs.iter().enumerate() {
            let c = &mut self.group_starts[o.idx as usize];
            self.order[*c as usize] = i as u32;
            *c += 1;
        }
        for i in (1..=n_idx).rev() {
            self.group_starts[i] = self.group_starts[i - 1];
        }
        self.group_starts[0] = 0;

        for row in rows {
            let idx = row
                .position
                .checked_sub(start_pos)
                .map(|rel| rel as usize * 2 + (row.raw_strand == '-') as usize)
                .filter(|&idx| idx < n_idx);
            let Some(idx) = idx else {
                out.push(ProbStats::default());
                continue;
            };
            let group = &self.order[self.group_starts[idx] as usize
                ..self.group_starts[idx + 1] as usize];

            self.scratch_values.clear();
            self.scratch_conf.clear();
            let combined_base = if combine_mods {
                match row.mod_code {
                    ModCodeRepr::Code(c) => {
                        DnaBase::parse(c).ok().map(|b| b as u8)
                    }
                    ModCodeRepr::ChEbi(_) => None,
                }
            } else {
                None
            };
            let row_code_id =
                self.codes.iter().position(|c| *c == row.mod_code);
            for &oi in group {
                let o = self.obs[oi as usize];
                let vals = &self.values
                    [o.start as usize..o.start as usize + o.n as usize];
                let mut mod_sum = 0f32;
                let mut max_mod = 0f32;
                let mut code_value = None;
                for &(code_id, p) in vals {
                    mod_sum += p;
                    max_mod = max_mod.max(p);
                    if Some(code_id as usize) == row_code_id {
                        code_value = Some(p);
                    }
                }
                let mod_sum = mod_sum.min(1f32);
                let value = if combine_mods {
                    (combined_base == Some(o.base)).then_some(mod_sum)
                } else {
                    code_value
                };
                if let Some(value) = value {
                    self.scratch_values.push(value);
                    self.scratch_conf.push((1f32 - mod_sum).max(max_mod));
                }
            }
            out.push(ProbStats::from_values(
                &mut self.scratch_values,
                &self.scratch_conf,
            ));
        }
    }
}

#[cfg(test)]
mod prob_stats_tests {
    use super::*;
    use crate::mod_base_code::{HYDROXY_METHYL_CYTOSINE, METHYL_CYTOSINE};

    fn row(
        position: u32,
        strand: char,
        code: ModCodeRepr,
    ) -> PileupFeatureCounts2 {
        PileupFeatureCounts2::new(
            position, strand, 1, code, 0, 0, 0, 0, 0, 0, 0, 0,
        )
    }

    #[test]
    fn test_quantiles_linear_interpolation() {
        // numpy.quantile([0.1, 0.2, 0.4, 0.8], [.25, .5, .75])
        // -> [0.175, 0.3, 0.5]
        let mut values = vec![0.8f32, 0.1, 0.4, 0.2];
        let stats = ProbStats::from_values(&mut values, &[1.0, 0.5]);
        assert_eq!(stats.n, 4);
        assert!((stats.mean - 0.375).abs() < 1e-6);
        assert!((stats.q25 - 0.175).abs() < 1e-6);
        assert!((stats.median - 0.3).abs() < 1e-6);
        assert!((stats.q75 - 0.5).abs() < 1e-6);
        // population sd
        let expected_sd = (((0.1f64 - 0.375).powi(2)
            + (0.2f64 - 0.375).powi(2)
            + (0.4f64 - 0.375).powi(2)
            + (0.8f64 - 0.375).powi(2))
            / 4f64)
            .sqrt();
        assert!((stats.sd as f64 - expected_sd).abs() < 1e-6);
        assert!((stats.mean_confidence - 0.75).abs() < 1e-6);
    }

    #[test]
    fn test_single_value() {
        let mut values = vec![0.9f32];
        let stats = ProbStats::from_values(&mut values, &[0.9]);
        assert_eq!(stats.n, 1);
        assert_eq!(stats.sd, 0f32);
        assert_eq!(stats.q25, 0.9);
        assert_eq!(stats.median, 0.9);
        assert_eq!(stats.q75, 0.9);
    }

    #[test]
    fn test_accumulator_per_code_strand_and_reset() {
        let mut acc = ProbAccumulator::default();
        acc.reset(10);
        // two reads on + strand at position 3, one on - strand
        acc.add(
            3,
            false,
            DnaBase::C,
            &[(HYDROXY_METHYL_CYTOSINE, 0.1), (METHYL_CYTOSINE, 0.7)],
        );
        acc.add(
            3,
            false,
            DnaBase::C,
            &[(HYDROXY_METHYL_CYTOSINE, 0.3), (METHYL_CYTOSINE, 0.1)],
        );
        acc.add(3, true, DnaBase::C, &[(METHYL_CYTOSINE, 0.5)]);

        let start = 100u32;
        let rows = vec![
            row(103, '+', HYDROXY_METHYL_CYTOSINE),
            row(103, '+', METHYL_CYTOSINE),
            row(103, '-', METHYL_CYTOSINE),
            row(103, '-', HYDROXY_METHYL_CYTOSINE),
            row(104, '+', METHYL_CYTOSINE),
        ];
        let mut out = Vec::new();
        acc.stats_for_rows(start, &rows, false, &mut out);
        assert_eq!(out.len(), rows.len());
        assert_eq!(out[0].n, 2);
        assert!((out[0].mean - 0.2).abs() < 1e-6);
        assert_eq!(out[1].n, 2);
        assert!((out[1].mean - 0.4).abs() < 1e-6);
        // confidences: read 1 max(0.2 can, 0.7) = 0.7, read 2 max(0.6,
        // 0.3) = 0.6
        assert!((out[1].mean_confidence - 0.65).abs() < 1e-6);
        assert_eq!(out[2].n, 1);
        assert!((out[2].mean - 0.5).abs() < 1e-6);
        assert_eq!(out[3].n, 0);
        assert_eq!(out[4].n, 0);

        // combined mods
        let rows = vec![row(103, '+', ModCodeRepr::Code('C'))];
        acc.stats_for_rows(start, &rows, true, &mut out);
        assert_eq!(out[0].n, 2);
        assert!((out[0].mean - 0.6).abs() < 1e-6);

        // nothing left over after reset
        acc.reset(10);
        let rows = vec![row(103, '+', METHYL_CYTOSINE)];
        acc.stats_for_rows(start, &rows, false, &mut out);
        assert_eq!(out[0].n, 0);
    }

    #[test]
    fn test_write_fixed4_matches_std_formatting() {
        let check = |v: f32| {
            let mut buf = Vec::new();
            write_fixed4(&mut buf, v).unwrap();
            assert_eq!(String::from_utf8(buf).unwrap(), format!("{v:.4}"));
        };
        // ties and near-ties: x/256 probabilities, x/2^k values
        for q in 0..=1024u32 {
            check(q as f32 / 1024f32);
            check((q as f32 + 0.5) / 256f32);
        }
        for v in [0f32, 1f32, 0.00005, 0.00015, 0.00025, 0.12345, 0.99995] {
            check(v);
        }
        // pseudo-random values in [0, 1]
        let mut x = 0x2545F491u32;
        for _ in 0..200_000 {
            x ^= x << 13;
            x ^= x >> 17;
            x ^= x << 5;
            check(f32::from_bits(x).abs().min(1f32).max(0f32));
            check((x as f64 / u32::MAX as f64) as f32);
        }
        check(f32::NAN);
    }

    #[test]
    fn test_write_columns() {
        let mut buf = Vec::new();
        ProbStats::default().write_columns(&mut buf).unwrap();
        assert_eq!(String::from_utf8(buf).unwrap(), "\t0\t.\t.\t.\t.\t.\t.");
        let mut buf = Vec::new();
        let mut values = vec![0.25f32, 0.75];
        ProbStats::from_values(&mut values, &[0.75, 0.75])
            .write_columns(&mut buf)
            .unwrap();
        assert_eq!(
            String::from_utf8(buf).unwrap(),
            "\t2\t0.5000\t0.2500\t0.3750\t0.5000\t0.6250\t0.7500"
        );
    }
}
