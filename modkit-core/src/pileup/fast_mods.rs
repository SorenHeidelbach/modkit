//! Allocation-free iteration over the explicit base modification calls of a
//! record, as a faster drop-in for htslib's `bam_mods_at_next_pos` (via
//! `rust_htslib::bam::Record::basemods_position_iter`).
//!
//! Only "simple" MM/ML tags are parsed here: an `MM` tag (not `Mm`) and a
//! `ML` array of u8, no `MN` tag, every MM record on the `+` strand with an
//! A, C, G, or T primary base, one ChEBI code or one or more single-letter
//! codes, an optional `?`/`.` mode, deltas that fit within the read, and
//! exactly as many ML values as calls, and at least one delta per MM record
//! (htslib reports a phantom call for records without any). Anything else is handed to htslib
//! unchanged, so behavior (including errors) is the same as before. Within
//! the simple subset the emitted positions, and per position the codes,
//! canonical bases, strands, and qualities (in order), are the same as
//! htslib's; this is checked by the differential tests at the bottom.

use rust_htslib::bam::{
    self, record::Aux, record::BaseModificationsPositionIter,
};

use crate::errs::{MkError, MkResult};
use crate::mod_base_code::ModCodeRepr;

/// One modification call reported at a position, mirrors the fields of
/// htslib's `hts_base_mod` that pileup uses.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct ModCall {
    pub(crate) mod_code: ModCodeRepr,
    /// The MM primary base (as in `hts_base_mod.canonical_base`)
    pub(crate) canonical_base: i32,
    /// 0 for '+', 1 for '-' (as in `hts_base_mod.strand`)
    pub(crate) strand: i32,
    /// ML value (as in `hts_base_mod.qual`)
    pub(crate) qual: i32,
}

#[derive(Debug, Clone, Copy)]
struct Entry {
    seq_pos: u32,
    qual: u8,
}

#[derive(Debug, Clone, Copy)]
struct MmRecord {
    base: u8,
    // index into `codes`
    codes_start: u16,
    n_codes: u16,
    // '.' or no mode (as `MmTagInfo::is_implicit`)
    implicit: bool,
    // this record's entries are entries[cursor..end], ascending SEQ
    // position, `n_codes` consecutive entries per position
    cursor: u32,
    end: u32,
}

/// Reusable buffers for parsing, owned by a worker.
#[derive(Debug, Default)]
pub(crate) struct FastModsBuffers {
    entries: Vec<Entry>,
    records: Vec<MmRecord>,
    codes: Vec<ModCodeRepr>,
    deltas: Vec<u32>,
    last_deltas: Vec<u32>,
    positions: Vec<u32>,
}

impl FastModsBuffers {
    /// The modification codes per primary base and whether the base has
    /// an implicit-mode MM record, the same as `get_base_codes_and_implicits`
    /// gives for `MmTagInfo::from_record` of a simple tag.
    pub(crate) fn base_codes_and_implicits(
        &self,
    ) -> ([Vec<ModCodeRepr>; 4], [bool; 4]) {
        let mut codes: [Vec<ModCodeRepr>; 4] = Default::default();
        let mut implicits = [false; 4];
        for r in self.records.iter() {
            let idx = match r.base {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => unreachable!("only ACGT records are parsed"),
            };
            codes[idx].extend_from_slice(
                &self.codes[r.codes_start as usize
                    ..(r.codes_start + r.n_codes) as usize],
            );
            if r.implicit {
                implicits[idx] = true;
            }
        }
        (codes, implicits)
    }
}

/// Iterator over explicit modification calls in SEQ coordinates, either
/// parsed here or delegated to htslib.
pub(crate) enum ModCallsIter<'a> {
    Fast { buffers: &'a mut FastModsBuffers },
    Hts(BaseModificationsPositionIter<'a>),
}

impl<'a> ModCallsIter<'a> {
    /// The parse buffers when the tags were parsed here (not by htslib).
    pub(crate) fn fast_buffers(&self) -> Option<&FastModsBuffers> {
        match self {
            Self::Fast { buffers } => Some(buffers),
            Self::Hts(_) => None,
        }
    }

    /// Equivalent to `record.basemods_position_iter()`, an error here is
    /// the same error htslib would give.
    pub(crate) fn new(
        record: &'a bam::Record,
        buffers: &'a mut FastModsBuffers,
    ) -> Result<Self, rust_htslib::errors::Error> {
        if parse_simple(record, buffers).is_some() {
            Ok(Self::Fast { buffers })
        } else {
            record.basemods_position_iter().map(Self::Hts)
        }
    }

    /// Next position with explicit calls, `calls` is cleared and filled.
    /// Returns None when done.
    #[inline]
    pub(crate) fn next_into(
        &mut self,
        calls: &mut Vec<ModCall>,
    ) -> Option<MkResult<usize>> {
        calls.clear();
        match self {
            Self::Fast { buffers, .. } => {
                // merge the (sorted) positions of each MM record, at equal
                // positions calls are reported in MM record order
                let entries = &buffers.entries;
                let mut pos = u32::MAX;
                for r in buffers.records.iter() {
                    if r.cursor < r.end {
                        pos = pos.min(entries[r.cursor as usize].seq_pos);
                    }
                }
                if pos == u32::MAX {
                    return None;
                }
                for r in buffers.records.iter_mut() {
                    if r.cursor < r.end
                        && entries[r.cursor as usize].seq_pos == pos
                    {
                        for k in 0..r.n_codes {
                            let e = entries[(r.cursor + k as u32) as usize];
                            calls.push(ModCall {
                                mod_code: buffers.codes
                                    [(r.codes_start + k) as usize],
                                canonical_base: r.base as i32,
                                strand: 0,
                                qual: e.qual as i32,
                            });
                        }
                        r.cursor += r.n_codes as u32;
                    }
                }
                Some(Ok(pos as usize))
            }
            Self::Hts(iter) => match iter.next()? {
                Ok((pos, codes)) => {
                    calls.extend(codes.iter().map(|c| ModCall {
                        mod_code: ModCodeRepr::from(c.modified_base),
                        canonical_base: c.canonical_base,
                        strand: c.strand,
                        qual: c.qual,
                    }));
                    Some(Ok(pos as usize))
                }
                Err(e) => Some(Err(MkError::HtsLibError(e))),
            },
        }
    }
}

#[inline]
fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        _ => b'N',
    }
}

#[inline]
fn seq_base(encoded: &[u8], i: usize) -> u8 {
    crate::pileup::base_mods_adapter::decode_base(encoded, i)
}

/// Parse a simple MM/ML pair into `buffers.entries`, sorted as htslib
/// reports them. Returns None (and the caller falls back to htslib) for
/// anything outside the simple subset.
fn parse_simple(
    record: &bam::Record,
    buffers: &mut FastModsBuffers,
) -> Option<()> {
    let seq_len = record.seq_len();
    if seq_len == 0 || seq_len >= u32::MAX as usize {
        return None;
    }
    if record.aux(b"MN").is_ok() {
        return None;
    }
    let mm = match record.aux(b"MM").ok()? {
        Aux::String(s) => s.as_bytes(),
        _ => return None,
    };
    let ml = match record.aux(b"ML").ok()? {
        Aux::ArrayU8(arr) => arr,
        _ => return None,
    };
    if mm.is_empty() {
        return None;
    }

    let FastModsBuffers {
        entries,
        records,
        codes,
        deltas,
        last_deltas,
        positions,
    } = buffers;
    entries.clear();
    records.clear();
    codes.clear();

    let reverse = record.is_reverse();
    let encoded = record.seq().encoded;

    let mut i = 0usize;
    let mut ml_idx = 0usize;
    let ml_len = ml.len();
    while i < mm.len() {
        // primary base and strand
        let base = mm[i];
        if !matches!(base, b'A' | b'C' | b'G' | b'T') {
            return None;
        }
        i += 1;
        if mm.get(i) != Some(&b'+') {
            return None;
        }
        i += 1;
        // codes
        let codes_start = codes.len();
        if mm.get(i)?.is_ascii_digit() {
            let mut v = 0u32;
            let mut n_digits = 0;
            while i < mm.len() && mm[i].is_ascii_digit() {
                v = v.checked_mul(10)?.checked_add((mm[i] - b'0') as u32)?;
                i += 1;
                n_digits += 1;
            }
            // htslib stores ChEBI codes as a negative int
            if n_digits > 9 || v == 0 || v > i32::MAX as u32 {
                return None;
            }
            codes.push(ModCodeRepr::ChEbi(v));
        } else {
            while i < mm.len() && mm[i].is_ascii_alphabetic() {
                codes.push(ModCodeRepr::Code(mm[i] as char));
                i += 1;
            }
        }
        let n_codes = codes.len() - codes_start;
        if n_codes == 0 || codes.len() > u16::MAX as usize {
            return None;
        }
        // mode
        let implicit = match mm.get(i)? {
            b'?' => {
                i += 1;
                false
            }
            b'.' => {
                i += 1;
                true
            }
            _ => true,
        };
        // deltas
        deltas.clear();
        loop {
            match mm.get(i)? {
                b';' => {
                    i += 1;
                    break;
                }
                b',' => {
                    i += 1;
                    let mut v = 0u32;
                    let mut n_digits = 0;
                    while i < mm.len() && mm[i].is_ascii_digit() {
                        v = v
                            .checked_mul(10)?
                            .checked_add((mm[i] - b'0') as u32)?;
                        i += 1;
                        n_digits += 1;
                    }
                    if n_digits == 0 {
                        return None;
                    }
                    deltas.push(v);
                }
                _ => return None,
            }
        }
        // htslib reports a (quality 0) call at the first matching base for
        // a record without any deltas, leave those to htslib
        if deltas.is_empty() {
            return None;
        }
        if records.len() >= u16::MAX as usize {
            return None;
        }

        // map deltas to SEQ positions; in read orientation the base is
        // `base`, for reverse reads that is the complement read from the end
        // of SEQ. MM records for the same base often share the deltas (e.g.
        // C+h?/C+m?), then the positions of the previous record are reused.
        let reuse = records.last().map(|r| r.base == base).unwrap_or(false)
            && last_deltas == deltas;
        if !reuse {
            positions.clear();
        }
        if !reuse && !deltas.is_empty() {
            let seq_base_target = if reverse { complement(base) } else { base };
            let mut d = 0usize;
            let mut skip = deltas[0];
            let mut take = |p: usize, positions: &mut Vec<u32>| -> bool {
                // returns true when all deltas are consumed
                if seq_base(encoded, p) == seq_base_target {
                    if skip == 0 {
                        positions.push(p as u32);
                        d += 1;
                        if d == deltas.len() {
                            return true;
                        }
                        skip = deltas[d];
                    } else {
                        skip -= 1;
                    }
                }
                false
            };
            let mut done = false;
            if reverse {
                for p in (0..seq_len).rev() {
                    if take(p, positions) {
                        done = true;
                        break;
                    }
                }
            } else {
                for p in 0..seq_len {
                    if take(p, positions) {
                        done = true;
                        break;
                    }
                }
            }
            if !done {
                // deltas run off the end of the read
                return None;
            }
            if reverse {
                positions.reverse();
            }
        }
        // ML values: for each position (in read order), one per code
        let n_ml = positions.len() * n_codes;
        if ml_idx + n_ml > ml_len {
            return None;
        }
        let entries_start = entries.len();
        let n_pos = positions.len();
        for (jj, &p) in positions.iter().enumerate() {
            // ML is in read order, for reverse reads that's descending SEQ
            let j = if reverse { n_pos - 1 - jj } else { jj };
            for k in 0..n_codes {
                let qual = ml.get(ml_idx + j * n_codes + k)?;
                entries.push(Entry {
                    seq_pos: p,
                    qual,
                });
            }
        }
        ml_idx += n_ml;
        std::mem::swap(last_deltas, deltas);
        records.push(MmRecord {
            base,
            codes_start: codes_start as u16,
            n_codes: n_codes as u16,
            implicit,
            cursor: entries_start as u32,
            end: entries.len() as u32,
        });
    }
    if ml_idx != ml_len {
        return None;
    }
    Some(())
}

#[cfg(test)]
mod fast_mods_tests {
    use rust_htslib::bam::{
        self,
        record::{Aux, Cigar, CigarString},
    };

    use super::*;

    struct Rng(u64);
    impl Rng {
        fn next(&mut self) -> u64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            self.0
        }
        fn below(&mut self, n: u64) -> u64 {
            self.next() % n
        }
    }

    fn make_record(seq: &[u8], reverse: bool, mm: &str, ml: &[u8]) -> bam::Record {
        let mut record = bam::Record::new();
        let qual = vec![30u8; seq.len()];
        record.set(
            b"read",
            Some(&CigarString(vec![Cigar::Match(seq.len() as u32)])),
            seq,
            &qual,
        );
        record.set_pos(100);
        record.set_tid(0);
        record.set_flags(if reverse { 16 } else { 0 });
        record.push_aux(b"MM", Aux::String(mm)).unwrap();
        record.push_aux(b"ML", Aux::ArrayU8(ml.into())).unwrap();
        record
    }

    fn hts_calls(record: &bam::Record) -> Option<Vec<(usize, Vec<ModCall>)>> {
        let iter = record.basemods_position_iter().ok()?;
        let mut out = Vec::new();
        for r in iter {
            let (pos, codes) = r.ok()?;
            out.push((
                pos as usize,
                codes
                    .iter()
                    .map(|c| ModCall {
                        mod_code: ModCodeRepr::from(c.modified_base),
                        canonical_base: c.canonical_base,
                        strand: c.strand,
                        qual: c.qual,
                    })
                    .collect(),
            ));
        }
        Some(out)
    }

    fn fast_calls(
        record: &bam::Record,
        buffers: &mut FastModsBuffers,
    ) -> (bool, Option<Vec<(usize, Vec<ModCall>)>>) {
        let is_fast = parse_simple(record, buffers).is_some();
        let Ok(mut iter) = ModCallsIter::new(record, buffers) else {
            return (is_fast, None);
        };
        let mut out = Vec::new();
        let mut calls = Vec::new();
        while let Some(r) = iter.next_into(&mut calls) {
            let Ok(pos) = r else {
                return (is_fast, None);
            };
            out.push((pos, calls.clone()));
        }
        (is_fast, Some(out))
    }

    /// Random reads and simple MM/ML tags (various bases, codes, modes, and
    /// numbers of records), compare against htslib.
    #[test]
    fn test_fast_mods_match_htslib() {
        let mut rng = Rng(0x2545F4914F6CDD1D);
        let mut buffers = FastModsBuffers::default();
        let code_choices: [&[&str]; 4] = [
            &["m", "h", "hm", "mh", "21839", "c", "C"],
            &["a", "17596", "A"],
            &["o", "g", "G"],
            &["e", "b", "T", "17802"],
        ];
        let bases = [b'C', b'A', b'G', b'T'];
        let mut n_fast = 0usize;
        let n_trials = 20_000usize;
        for _ in 0..n_trials {
            let len = 1 + rng.below(120) as usize;
            let seq = (0..len)
                .map(|_| b"ACGT"[rng.below(4) as usize])
                .collect::<Vec<u8>>();
            let reverse = rng.below(2) == 1;
            // read-orientation sequence, MM refers to this
            let read_seq = if reverse {
                seq.iter().rev().map(|&b| complement(b)).collect::<Vec<_>>()
            } else {
                seq.clone()
            };
            let n_records = 1 + rng.below(3) as usize;
            let mut mm = String::new();
            let mut ml = Vec::new();
            let mut used_bases = Vec::new();
            for _ in 0..n_records {
                let bi = rng.below(4) as usize;
                let base = bases[bi];
                // repeated primary bases are allowed but codes must differ
                let code = code_choices[bi][rng.below(
                    code_choices[bi].len() as u64,
                ) as usize];
                if used_bases.iter().any(|(b, c)| *b == base && *c == code) {
                    continue;
                }
                used_bases.push((base, code));
                let n_codes = if code.chars().all(|c| c.is_ascii_digit()) {
                    1
                } else {
                    code.len()
                };
                let mode = ["", "?", "."][rng.below(3) as usize];
                let occurrences =
                    read_seq.iter().filter(|&&b| b == base).count();
                let mut deltas = Vec::new();
                let mut remaining = occurrences;
                while remaining > 0 && rng.below(4) != 0 {
                    let d = rng.below(remaining as u64) as usize;
                    // d skipped, one called
                    if d + 1 > remaining {
                        break;
                    }
                    deltas.push(d);
                    remaining -= d + 1;
                }
                mm.push(base as char);
                mm.push('+');
                mm.push_str(code);
                mm.push_str(mode);
                for d in deltas.iter() {
                    mm.push_str(&format!(",{d}"));
                }
                mm.push(';');
                for _ in 0..(deltas.len() * n_codes) {
                    ml.push(rng.below(256) as u8);
                }
                // sometimes a second record for the same base with the same
                // deltas (like C+h?/C+m?)
                if rng.below(2) == 0 {
                    let code2 = code_choices[bi][rng.below(
                        code_choices[bi].len() as u64,
                    ) as usize];
                    if !used_bases.iter().any(|(b, c)| *b == base && *c == code2)
                    {
                        used_bases.push((base, code2));
                        let n_codes2 =
                            if code2.chars().all(|c| c.is_ascii_digit()) {
                                1
                            } else {
                                code2.len()
                            };
                        mm.push(base as char);
                        mm.push('+');
                        mm.push_str(code2);
                        mm.push_str(["", "?", "."][rng.below(3) as usize]);
                        for d in deltas.iter() {
                            mm.push_str(&format!(",{d}"));
                        }
                        mm.push(';');
                        for _ in 0..(deltas.len() * n_codes2) {
                            ml.push(rng.below(256) as u8);
                        }
                    }
                }
            }
            if mm.is_empty() {
                continue;
            }
            let record = make_record(&seq, reverse, &mm, &ml);
            let expected = hts_calls(&record);
            let (is_fast, observed) = fast_calls(&record, &mut buffers);
            n_fast += is_fast as usize;
            if is_fast {
                parse_simple(&record, &mut buffers).unwrap();
                let infos =
                    crate::mod_bam::MmTagInfo::from_record(&record).unwrap();
                assert_eq!(
                    buffers.base_codes_and_implicits(),
                    crate::pileup::pileup_processor::get_base_codes_and_implicits(
                        &infos
                    ),
                    "{mm}"
                );
            }
            assert_eq!(
                observed, expected,
                "mm={mm} ml={ml:?} reverse={reverse} seq={}",
                String::from_utf8_lossy(&seq)
            );
        }
        // most generated tags are "simple"
        assert!(n_fast > n_trials / 2, "{n_fast}");
    }

    /// Tags outside the simple subset fall back to htslib.
    #[test]
    fn test_fast_mods_fallback() {
        let mut buffers = FastModsBuffers::default();
        let seq = b"ACGTCCGACGATCG";
        for (mm, ml) in [
            // duplex
            ("C-m?,0;", vec![200u8]),
            // N primary base
            ("N+n?,0;", vec![200u8]),
            // too few ML values
            ("C+m?,0,1;", vec![200u8]),
            // too many ML values
            ("C+m?,0;", vec![200u8, 3]),
            // deltas past the end of the read
            ("C+m?,100;", vec![200u8]),
            // malformed
            ("C+m?0;", vec![200u8]),
            ("C+m?,0", vec![200u8]),
        ] {
            for reverse in [false, true] {
                let record = make_record(seq, reverse, mm, &ml);
                assert!(
                    parse_simple(&record, &mut buffers).is_none(),
                    "{mm}"
                );
                let expected = hts_calls(&record);
                let (_, observed) = fast_calls(&record, &mut buffers);
                assert_eq!(observed, expected, "{mm}");
            }
        }
    }
}
