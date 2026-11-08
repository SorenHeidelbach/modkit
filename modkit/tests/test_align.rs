use std::collections::HashMap;

use rust_htslib::bam::record::{Aux, Cigar};
use rust_htslib::bam::{self, Read};

use crate::common::run_modkit;

mod common;

const REF: &str = "../tests/resources/CGI_ladder_3.6kb_ref.fa";
// aligned with minimap2 (mappy) from the same reads, used as the expected
// alignments
const EXPECTED_ALIGNED: &str =
    "../tests/resources/bc_anchored_10_reads.sorted.bam";
const UNALIGNED: &str = "../tests/resources/bc_anchored_10_reads.unmapped.bam";

fn aux_string(record: &bam::Record, tag: &[u8]) -> Option<String> {
    match record.aux(tag) {
        Ok(Aux::String(s)) => Some(s.to_string()),
        _ => None,
    }
}

fn aux_u8_array(record: &bam::Record, tag: &[u8]) -> Option<Vec<u8>> {
    match record.aux(tag) {
        Ok(Aux::ArrayU8(arr)) => Some(arr.iter().collect()),
        _ => None,
    }
}

fn aux_int(record: &bam::Record, tag: &[u8]) -> Option<i64> {
    match record.aux(tag) {
        Ok(Aux::I8(x)) => Some(x as i64),
        Ok(Aux::U8(x)) => Some(x as i64),
        Ok(Aux::I16(x)) => Some(x as i64),
        Ok(Aux::U16(x)) => Some(x as i64),
        Ok(Aux::I32(x)) => Some(x as i64),
        Ok(Aux::U32(x)) => Some(x as i64),
        _ => None,
    }
}

fn qname(record: &bam::Record) -> String {
    String::from_utf8_lossy(record.qname()).to_string()
}

fn read_records(fp: &std::path::Path) -> Vec<bam::Record> {
    let mut reader = bam::Reader::from_path(fp).unwrap();
    reader.records().map(|r| r.unwrap()).collect()
}

/// Every aligned record must carry the full read sequence (no hard clips) and
/// MN must equal the sequence length, so MM/ML stay valid.
fn check_records_have_valid_mod_tags(records: &[bam::Record]) {
    for record in records {
        let name = qname(record);
        assert!(
            !record.cigar().iter().any(|op| matches!(op, Cigar::HardClip(_))),
            "{name} is hard-clipped"
        );
        assert!(aux_string(record, b"MM").is_some(), "{name} missing MM");
        assert!(aux_u8_array(record, b"ML").is_some(), "{name} missing ML");
        assert_eq!(
            aux_int(record, b"MN"),
            Some(record.seq_len() as i64),
            "{name} MN tag does not match sequence length"
        );
        // stale alignment tags from the input must not be duplicated
        for tag in [b"NM", b"MD"] {
            let n = record
                .aux_iter()
                .filter_map(|r| r.ok())
                .filter(|(t, _)| t == tag)
                .count();
            assert!(n <= 1, "{name} has {n} {:?} tags", tag);
        }
    }
}

fn check_primary_against_expected(records: &[bam::Record]) {
    let expected = read_records(std::path::Path::new(EXPECTED_ALIGNED))
        .into_iter()
        .filter(|r| !r.is_secondary() && !r.is_supplementary())
        .map(|r| (qname(&r), r))
        .collect::<HashMap<String, bam::Record>>();
    let observed = records
        .iter()
        .filter(|r| !r.is_secondary() && !r.is_supplementary())
        .map(|r| (qname(r), r))
        .collect::<HashMap<String, &bam::Record>>();
    assert_eq!(observed.len(), expected.len());
    for (name, exp) in expected.iter() {
        let obs = observed.get(name).expect("missing read");
        assert_eq!(obs.tid(), exp.tid(), "{name}");
        assert_eq!(obs.pos(), exp.pos(), "{name}");
        assert_eq!(obs.is_reverse(), exp.is_reverse(), "{name}");
        assert_eq!(obs.cigar().to_string(), exp.cigar().to_string(), "{name}");
        assert_eq!(obs.seq().as_bytes(), exp.seq().as_bytes(), "{name}");
        assert_eq!(obs.qual(), exp.qual(), "{name}");
        assert_eq!(aux_string(obs, b"MM"), aux_string(exp, b"MM"), "{name}");
        assert_eq!(
            aux_u8_array(obs, b"ML"),
            aux_u8_array(exp, b"ML"),
            "{name}"
        );
    }
}

#[test]
fn test_align_help() {
    run_modkit(&["modbam", "align", "--help"]).unwrap();
}

#[test]
fn test_align_unaligned_input() {
    let out_dir = tempfile::tempdir().unwrap();
    let out_bam = out_dir.path().join("aligned.bam");
    run_modkit(&[
        "modbam",
        "align",
        "-r",
        REF,
        UNALIGNED,
        out_bam.to_str().unwrap(),
    ])
    .unwrap();
    let records = read_records(&out_bam);
    assert!(records.iter().all(|r| !r.is_unmapped()));
    // one read is chimeric, it has a supplementary alignment
    let supplementary =
        records.iter().filter(|r| r.is_supplementary()).collect::<Vec<_>>();
    assert_eq!(supplementary.len(), 1);
    let supp = supplementary[0];
    assert!(aux_string(supp, b"SA").is_some());
    let primary = records
        .iter()
        .find(|r| qname(r) == qname(supp) && !r.is_supplementary())
        .unwrap();
    assert!(aux_string(primary, b"SA").is_some());
    assert_eq!(primary.seq_len(), supp.seq_len());

    check_records_have_valid_mod_tags(&records);
    check_primary_against_expected(&records);
}

#[test]
fn test_align_realign_aligned_input() {
    // re-aligning already aligned reads (including reverse strand reads)
    // should give the same result as aligning the unaligned reads.
    let out_dir = tempfile::tempdir().unwrap();
    let out_bam = out_dir.path().join("realigned.bam");
    run_modkit(&[
        "modbam",
        "align",
        "-r",
        REF,
        EXPECTED_ALIGNED,
        out_bam.to_str().unwrap(),
    ])
    .unwrap();
    let records = read_records(&out_bam);
    check_records_have_valid_mod_tags(&records);
    check_primary_against_expected(&records);
}

#[test]
fn test_align_no_unmapped() {
    // align to a reference that the reads do not match, all reads are
    // unmapped
    let out_dir = tempfile::tempdir().unwrap();
    let ref_fp = out_dir.path().join("ref.fa");
    std::fs::write(
        &ref_fp,
        ">poly\n".to_string() + &"ACGTTGCA".repeat(100) + "\n",
    )
    .unwrap();

    let out_bam = out_dir.path().join("with_unmapped.bam");
    run_modkit(&[
        "modbam",
        "align",
        "-r",
        ref_fp.to_str().unwrap(),
        UNALIGNED,
        out_bam.to_str().unwrap(),
    ])
    .unwrap();
    let records = read_records(&out_bam);
    assert_eq!(records.len(), 10);
    assert!(records.iter().all(|r| r.is_unmapped()));
    check_records_have_valid_mod_tags(&records);

    let out_bam = out_dir.path().join("no_unmapped.bam");
    run_modkit(&[
        "modbam",
        "align",
        "--no-unmapped",
        "-r",
        ref_fp.to_str().unwrap(),
        UNALIGNED,
        out_bam.to_str().unwrap(),
    ])
    .unwrap();
    assert!(read_records(&out_bam).is_empty());
}

/// `pileup --align` (align in-process, pile up from memory) must give the
/// same bedMethyl as `pileup` on the equivalent aligned, sorted BAM.
#[test]
fn test_pileup_align_matches_bam_pileup() {
    let out_dir = tempfile::tempdir().unwrap();
    for extra_args in [
        vec![],
        vec!["--no-filtering"],
        vec!["--cpg", "--modified-bases", "5mC", "5hmC", "--ref", REF],
        vec!["--motif", "CG", "0", "--combine-strands", "--ref", REF],
    ] {
        let from_bam = out_dir.path().join("from_bam.bed");
        let streamed = out_dir.path().join("streamed.bed");
        let mut args =
            vec!["pileup", EXPECTED_ALIGNED, from_bam.to_str().unwrap()];
        args.extend(extra_args.iter());
        run_modkit(&args).unwrap();
        let mut args = vec![
            "pileup",
            "--align",
            REF,
            UNALIGNED,
            streamed.to_str().unwrap(),
        ];
        args.extend(extra_args.iter());
        run_modkit(&args).unwrap();
        let expected = std::fs::read_to_string(&from_bam).unwrap();
        let observed = std::fs::read_to_string(&streamed).unwrap();
        assert!(!expected.is_empty(), "{extra_args:?}");
        assert_eq!(observed, expected, "{extra_args:?}");
    }
}
