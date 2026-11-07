use anyhow::Context;
use itertools::Itertools;
use rust_htslib::bam;
use std::cmp::Ordering;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

use common::{check_against_expected_text_file, run_modkit};
use mod_kit::mod_base_code::{ModCodeRepr, METHYL_CYTOSINE};

mod common;

#[test]
fn test_pileup_help() {
    let pileup_help_args = ["pileup", "--help"];
    let _out = run_modkit(&pileup_help_args).unwrap();
}

#[test]
fn test_pileup_no_filt() {
    let temp_file = std::env::temp_dir().join("test_pileup_nofilt.bed");
    let args = [
        "pileup",
        "-i",
        "25", // use small interval to make sure chunking works
        "--no-filtering",
        "../tests/resources/bc_anchored_10_reads.sorted.bam",
        temp_file.to_str().unwrap(),
    ];

    run_modkit(&args).unwrap();

    check_against_expected_text_file(
        temp_file.to_str().unwrap(),
        "../tests/resources/modbam.modpileup_nofilt.methyl.bed",
    );
}

#[test]
fn test_pileup_with_filt() {
    let temp_file = std::env::temp_dir().join("test_pileup_withfilt.bed");
    let args = [
        "pileup",
        "-i",
        "25", // use small interval to make sure chunking works
        "-f",
        "1.0",
        "-p",
        "0.25",
        "--seed",
        "42",
        "--include-unmapped",
        "../tests/resources/bc_anchored_10_reads.sorted.bam",
        temp_file.to_str().unwrap(),
    ];

    run_modkit(&args).unwrap();

    check_against_expected_text_file(
        temp_file.to_str().unwrap(),
        "../tests/resources/modbam.modpileup_filt025.methyl.bed",
    );
}

#[test]
fn test_pileup_combine() {
    let test_adjusted_bam = std::env::temp_dir().join("test_combined.bed");
    let pileup_args = [
        "pileup",
        "--combine-mods",
        "--no-filtering",
        "../tests/resources/bc_anchored_10_reads.sorted.bam",
        test_adjusted_bam.to_str().unwrap(),
    ];
    run_modkit(&pileup_args).unwrap();
    assert!(test_adjusted_bam.exists());

    check_against_expected_text_file(
        test_adjusted_bam.to_str().unwrap(),
        "../tests/resources/modbam.modpileup_combined.methyl.bed",
    );
}


#[test]
fn test_pileup_no_mod_calls() {
    let empty_bedfile =
        std::env::temp_dir().join("test_pileup_no_mod_calls_outbed.bed");
    let args = [
        "pileup",
        "--no-filtering",
        "../tests/resources/empty-tags.sorted.bam",
        empty_bedfile.to_str().unwrap(),
    ];

    run_modkit(&args).unwrap();

    let reader = BufReader::new(File::open(empty_bedfile).unwrap());
    let lines = reader.lines().collect::<Vec<Result<String, _>>>();
    assert_eq!(lines.len(), 0);
}

#[test]
fn test_pileup_old_tags() {
    // generated from HG002_small.ch20._other.sorted.bam with
    // `modkit update-tags --mode ambiguous --no-implicit-probs` (v0.6.4)
    let updated_file =
        "../tests/resources/HG002_small.ch20._other.updated_tags.sorted.bam";
    let out_file = std::env::temp_dir().join("test_pileup_old_tags.bed");
    run_modkit(&[
        "pileup",
        "--no-filtering",
        updated_file,
        out_file.to_str().unwrap(),
    ])
    .unwrap();
    assert!(out_file.exists());
    check_against_expected_text_file(
        out_file.to_str().unwrap(),
        "../tests/resources/pileup-old-tags-regressiontest.methyl.bed",
    );
}

#[test]
fn test_pileup_with_region() {
    let temp_file = std::env::temp_dir().join("test_pileup_with_region.bed");
    let args = [
        "pileup",
        "-i",
        "25", // use small interval to make sure chunking works
        "--no-filtering",
        "--region",
        "oligo_1512_adapters:0-50",
        "../tests/resources/bc_anchored_10_reads.sorted.bam",
        temp_file.to_str().unwrap(),
    ];

    run_modkit(&args).unwrap();

    check_against_expected_text_file(
        temp_file.to_str().unwrap(),
        "../tests/resources/modbam.modpileup_nofilt_oligo_1512_adapters_10_50.\
         bed",
    );
}

#[test]
fn test_pileup_duplex_reads() {
    let temp_file = std::env::temp_dir().join("test_pileup_duplex_reads.bed");
    run_modkit(&[
        "pileup",
        "../tests/resources/duplex_modbam.sorted.bam",
        temp_file.to_str().unwrap(),
        "--region",
        "chr17",
        "--no-filtering",
    ])
    .unwrap();

    check_against_expected_text_file(
        temp_file.to_str().unwrap(),
        "../tests/resources/duplex_modbam_pileup_nofilt.bed",
    );
}

#[test]
fn test_pileup_cpg_motif_filtering() {
    let temp_file = std::env::temp_dir().join("test_cpg_motif_filtering.bed");
    run_modkit(&[
        "pileup",
        "../tests/resources/bc_anchored_10_reads.sorted.bam",
        temp_file.to_str().unwrap(),
        "--no-filtering",
        "--cpg",
        "--modified-bases",
        "5mC",
        "5hmC",
        "--ref",
        "../tests/resources/CGI_ladder_3.6kb_ref.fa",
    ])
    .unwrap();
    check_against_expected_text_file(
        temp_file.to_str().unwrap(),
        "../tests/resources/bc_anchored_10_reads_nofilt_cg_motif.bed",
    );
}

#[test]
fn test_pileup_cpg_motif_filtering_compressed_ref() {
    let temp_file = std::env::temp_dir()
        .join("test_pileup_cpg_motif_filtering_compressed_ref.bed");
    run_modkit(&[
        "pileup",
        "../tests/resources/bc_anchored_10_reads.sorted.bam",
        temp_file.to_str().unwrap(),
        "--no-filtering",
        "--cpg",
        "--modified-bases",
        "5mC",
        "5hmC",
        "--ref",
        "../tests/resources/CGI_ladder_3.6kb_ref.fa.gz",
    ])
    .unwrap();
    check_against_expected_text_file(
        temp_file.to_str().unwrap(),
        "../tests/resources/bc_anchored_10_reads_nofilt_cg_motif.bed",
    );
}

#[test]
fn test_pileup_cpg_motif_filtering_strand_combine() {
    let temp_file = std::env::temp_dir()
        .join("test_cpg_motif_filtering_strand_combine.bed");
    let interval_sizes =
        ["10", "88", "89", "90", "91", "92", "93", "94", "10000"];
    for interval_size in interval_sizes {
        run_modkit(&[
            "pileup",
            "../tests/resources/bc_anchored_10_reads.sorted.bam",
            temp_file.to_str().unwrap(),
            "--no-filtering",
            "-i",
            interval_size,
            "--cpg",
            "--combine-strands",
            "--modified-bases",
            "5mC",
            "5hmC",
            "--ref",
            "../tests/resources/CGI_ladder_3.6kb_ref.fa",
        ])
        .unwrap();
        check_against_expected_text_file(
            temp_file.to_str().unwrap(),
            "../tests/resources/\
             bc_anchored_10_reads_nofilt_cg_motif_strand_combine.bed",
        );
    }
}







#[test]
fn test_pileup_with_filt_position_filter() {
    let temp_file =
        std::env::temp_dir().join("test_pileup_with_filt_position_filter.bed");
    run_modkit(&[
        "pileup",
        "-i",
        "25", // use small interval to make sure chunking works
        "--filter-threshold",
        "0.6894531",
        "--modified-bases",
        "5mC",
        "5hmC",
        "--ref",
        "../tests/resources/CGI_ladder_3.6kb_ref.fa",
        "--include-positions",
        "../tests/resources/CGI_ladder_3.6kb_ref_include_positions.bed",
        "../tests/resources/bc_anchored_10_reads.sorted.bam",
        temp_file.to_str().unwrap(),
    ])
    .unwrap();

    check_against_expected_text_file(
        temp_file.to_str().unwrap(),
        "../tests/resources/modbam.modpileup_filt_positions_fixed_thresh.\
         methyl.bed",
    );
    run_modkit(&[
        "pileup",
        "-i",
        "25", // use small interval to make sure chunking works
        "-p",
        "0.25",
        "--modified-bases",
        "5mC",
        "5hmC",
        "--ref",
        "../tests/resources/CGI_ladder_3.6kb_ref.fa",
        "--include-positions",
        "../tests/resources/CGI_ladder_3.6kb_ref_include_positions.bed",
        "../tests/resources/bc_anchored_10_reads.sorted.bam",
        temp_file.to_str().unwrap(),
    ])
    .unwrap();

    check_against_expected_text_file(
        temp_file.to_str().unwrap(),
        "../tests/resources/modbam.modpileup_filt_positions_025.methyl.bed",
    );
}



#[test]
fn test_pileup_motifs_cg0_cgcg2() {
    let temp_file =
        std::env::temp_dir().join("test_pileup_motifs_cg0_cgcg2.bed");
    run_modkit(&[
        "pileup",
        "../tests/resources/CG_5mC_20230207_1700_6A_PAG66026_3c0abf27_oligo_741_adapters_modcalls_0th_sort_10_reads.bam",
        temp_file.to_str().unwrap(),
        "--motif", "CG", "0",
        "--motif", "CGCG", "2",
        "--no-filtering",
        "--ref", "../tests/resources/CGI_ladder_3.6kb_ref.fa",
        "--region", "oligo_741_adapters:22-62",
    ])
        .unwrap();

    check_against_expected_text_file(
        temp_file.to_str().unwrap(),
        "../tests/resources/cgcg2_cg0_test1.bed",
    );

    run_modkit(&[
        "pileup",
        "../tests/resources/CG_5mC_20230207_1700_6A_PAG66026_3c0abf27_oligo_741_adapters_modcalls_0th_sort_10_reads-2.bam",
        temp_file.to_str().unwrap(),
        "--motif", "CG", "0",
        "--motif", "CGCG", "2",
        "--no-filtering",
        "--ref", "../tests/resources/CGI_ladder_3.6kb_ref.fa",
        "--region", "oligo_741_adapters:22-62",
    ])
        .unwrap();

    check_against_expected_text_file(
        temp_file.to_str().unwrap(),
        "../tests/resources/cgcg2_cg0_test2.bed",
    );
}


#[test]
fn test_pileup_chebi_code_same_output() {
    let sorted_fields = |lines: Vec<String>| -> Vec<Vec<String>> {
        lines
            .into_iter()
            .map(|l| l.split('\t').map(|x| x.to_string()).collect::<Vec<_>>())
            .sorted_by(|a, b| {
                (&a[0], a[1].parse::<u64>().unwrap(), &a[5], &a[3]).cmp(&(
                    &b[0],
                    b[1].parse::<u64>().unwrap(),
                    &b[5],
                    &b[3],
                ))
            })
            .collect()
    };
    let pileup =
        std::env::temp_dir().join("test_chebi_code_same_output_pileup.bed");
    let expected = sorted_fields(read_lines(&PathBuf::from(
        "../tests/resources/modbam.modpileup_nofilt.methyl.bed",
    )));
    // generated from bc_anchored_10_reads.sorted.bam with
    // `modkit adjust-mods --convert h <code>` (v0.6.4)
    for (to_code, adjusted_bam) in [
        (
            ModCodeRepr::ChEbi(76792),
            "../tests/resources/bc_anchored_10_reads.h_to_chebi76792.sorted.bam",
        ),
        (
            ModCodeRepr::Code('c'),
            "../tests/resources/bc_anchored_10_reads.h_to_c.sorted.bam",
        ),
    ] {
        run_modkit(&[
            "pileup",
            adjusted_bam,
            pileup.to_str().unwrap(),
            "-i",
            "25", // use small interval to make sure chunking works
            "--no-filtering",
        ])
        .context("failed to generate pileup")
        .unwrap();

        let observed = read_lines(&pileup)
            .into_iter()
            .map(|l| {
                let mut fields =
                    l.split('\t').map(|x| x.to_string()).collect::<Vec<_>>();
                if fields[3] != METHYL_CYTOSINE.to_string() {
                    assert_eq!(fields[3], to_code.to_string());
                    fields[3] = "h".to_string();
                }
                fields.join("\t")
            })
            .collect::<Vec<_>>();
        assert_eq!(expected, sorted_fields(observed));
    }
}

#[test]
fn test_pileup_with_header() {
    let temp_file = std::env::temp_dir().join("test_pileup_nofilt.bed");
    let args = [
        "pileup",
        "-i",
        "25", // use small interval to make sure chunking works
        "--no-filtering",
        "--with-header",
        "../tests/resources/bc_anchored_10_reads.sorted.bam",
        temp_file.to_str().unwrap(),
    ];

    run_modkit(&args).unwrap();

    check_against_expected_text_file(
        temp_file.to_str().unwrap(),
        "../tests/resources/pileup_with_header.bed",
    );
}

#[test]
fn test_pilep_phased() {
    let test_tmp_dir = std::env::temp_dir();
    let phased_output_dir = test_tmp_dir.join("test_pilep_phased_pileup");
    let normal_output_pileup =
        test_tmp_dir.join("test_pilep_regular_pileup.bedmethyl");
    let args = [
        "pileup",
        "../tests/resources/test.sorted.phased.bam",
        phased_output_dir.to_str().unwrap(),
        "--modified-bases",
        "5mC",
        "--phased",
        "--ref",
        "../tests/resources/genome_for_phased_test.fasta",
        "--no-filtering",
    ];
    run_modkit(&args).unwrap();
    let args = [
        "pileup",
        "../tests/resources/test.sorted.phased.bam",
        normal_output_pileup.to_str().unwrap(),
        "--modified-bases",
        "5mC",
        "--ref",
        "../tests/resources/genome_for_phased_test.fasta",
        "--no-filtering",
    ];
    run_modkit(&args).unwrap();
    let combined_phased = phased_output_dir.join("combined.bedmethyl");
    assert!(combined_phased.exists() && combined_phased.is_file());
    assert!(normal_output_pileup.exists() && normal_output_pileup.is_file());
    check_against_expected_text_file(
        combined_phased.to_str().unwrap(),
        normal_output_pileup.to_str().unwrap(),
    );
}

fn read_lines(fp: &PathBuf) -> Vec<String> {
    BufReader::new(File::open(fp).unwrap())
        .lines()
        .map(|l| l.unwrap())
        .collect()
}

/// Runs pileup with and without --prob-stats and checks that the standard
/// columns are unchanged, the extra columns are well formed, and that every
/// read counted at a position (valid coverage + failed calls) contributed a
/// probability.
fn check_prob_stats_against_counts(name: &str, args: &[&str]) -> Vec<String> {
    let without = std::env::temp_dir().join(format!("{name}_without.bed"));
    let with = std::env::temp_dir().join(format!("{name}_with.bed"));
    let mut a = vec!["pileup"];
    a.extend_from_slice(args);
    let mut a_without = a.clone();
    a_without.push(without.to_str().unwrap());
    run_modkit(&a_without).unwrap();
    let mut a_with = a.clone();
    a_with.extend_from_slice(&["--prob-stats", with.to_str().unwrap()]);
    run_modkit(&a_with).unwrap();

    let without = read_lines(&without);
    let with = read_lines(&with);
    assert!(!with.is_empty());
    assert_eq!(without.len(), with.len());
    for (a, b) in without.iter().zip(with.iter()) {
        let fa = a.split('\t').collect::<Vec<_>>();
        let fb = b.split('\t').collect::<Vec<_>>();
        assert_eq!(fa.len(), 18);
        assert_eq!(fb.len(), 25, "{b}");
        assert_eq!(fa[..], fb[..18], "standard columns changed");
        let valid = fa[9].parse::<u32>().unwrap();
        let fail = fa[15].parse::<u32>().unwrap();
        let n = fb[18].parse::<u32>().unwrap();
        assert_eq!(n, valid + fail, "{b}");
        let stats = fb[19..]
            .iter()
            .map(|x| x.parse::<f32>().unwrap())
            .collect::<Vec<f32>>();
        let (mean, sd, q25, median, q75, conf) =
            (stats[0], stats[1], stats[2], stats[3], stats[4], stats[5]);
        for x in [mean, q25, median, q75] {
            assert!((0f32..=1f32).contains(&x), "{b}");
        }
        assert!(sd >= 0f32, "{b}");
        assert!(q25 <= median && median <= q75, "{b}");
        // the most likely class has probability of at least 1/N classes
        assert!(conf >= 0.33 && conf <= 1f32, "{b}");
    }
    with
}

#[test]
fn test_pileup_prob_stats() {
    // generic worker, small interval to check chunking
    let rows = check_prob_stats_against_counts(
        "test_pileup_prob_stats_generic",
        &[
            "-i",
            "25",
            "--no-filtering",
            "../tests/resources/bc_anchored_10_reads.sorted.bam",
        ],
    );

    // independently compute n and mean probability per position, strand,
    // and modification code from the read-level probabilities
    // generated with `modkit extract full bc_anchored_10_reads.sorted.bam`
    // (v0.6.4)
    let extract_fp = PathBuf::from(
        "../tests/resources/bc_anchored_10_reads.extract_full.tsv",
    );
    let extract = read_lines(&extract_fp);
    let header = extract[0].split('\t').collect::<Vec<_>>();
    let col = |name: &str| header.iter().position(|x| *x == name).unwrap();
    let (chrom, pos, strand, code, prob) = (
        col("chrom"),
        col("ref_position"),
        col("ref_mod_strand"),
        col("mod_code"),
        col("mod_qual"),
    );
    let mut expected: HashMap<(String, i64, String, String), Vec<f64>> =
        HashMap::new();
    for line in extract.iter().skip(1) {
        let f = line.split('\t').collect::<Vec<_>>();
        let p = f[pos].parse::<i64>().unwrap();
        if p < 0 {
            continue;
        }
        expected
            .entry((
                f[chrom].to_string(),
                p,
                f[strand].to_string(),
                f[code].to_string(),
            ))
            .or_default()
            .push(f[prob].parse::<f64>().unwrap());
    }
    for row in rows {
        let f = row.split('\t').collect::<Vec<_>>();
        let key = (
            f[0].to_string(),
            f[1].parse::<i64>().unwrap(),
            f[5].to_string(),
            f[3].to_string(),
        );
        let probs = expected.get(&key).expect("missing position in extract");
        assert_eq!(f[18].parse::<usize>().unwrap(), probs.len(), "{row}");
        let mean = probs.iter().sum::<f64>() / probs.len() as f64;
        let observed = f[19].parse::<f64>().unwrap();
        assert!((observed - mean).abs() < 1e-4, "{row} expected {mean}");
    }

    // optimized CpG workers, combining strands
    check_prob_stats_against_counts(
        "test_pileup_prob_stats_cpg",
        &[
            "--no-filtering",
            "--ref",
            "../tests/resources/CGI_ladder_3.6kb_ref.fa",
            "--cpg",
            "--modified-bases",
            "5mC",
            "5hmC",
            "--combine-strands",
            "../tests/resources/bc_anchored_10_reads.sorted.bam",
        ],
    );
    // generic workers, combined modification codes
    check_prob_stats_against_counts(
        "test_pileup_prob_stats_combine_mods",
        &[
            "--combine-mods",
            "--filter-threshold",
            "0.7",
            "../tests/resources/bc_anchored_10_reads.sorted.bam",
        ],
    );

    // header has the extra fields
    let out = std::env::temp_dir().join("test_pileup_prob_stats_header.bed");
    run_modkit(&[
        "pileup",
        "--no-filtering",
        "--header",
        "--prob-stats",
        "../tests/resources/bc_anchored_10_reads.sorted.bam",
        out.to_str().unwrap(),
    ])
    .unwrap();
    let lines = read_lines(&out);
    let header = lines[0].split('\t').collect::<Vec<_>>();
    assert_eq!(header.len(), 25);
    assert_eq!(
        header[18..],
        [
            "n_prob",
            "mean_prob",
            "sd_prob",
            "q25_prob",
            "median_prob",
            "q75_prob",
            "mean_confidence"
        ]
    );
}
