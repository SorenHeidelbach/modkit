# Aligning a modBAM to a reference

`modkit modbam align` aligns the reads in a modBAM (unaligned, or previously aligned) to a reference with [minimap2](https://github.com/lh3/minimap2), which is built into modkit, so no external aligner or FASTQ conversion is needed.
The base modification tags (`MM`/`ML`) and all other per-read tags (e.g. `RG`, `mv`, `qs`) are carried over to every output record.

```bash
modkit modbam align --ref reference.fa calls.bam aligned.bam
samtools sort -o aligned.sorted.bam aligned.bam && samtools index aligned.sorted.bam
modkit pileup aligned.sorted.bam pileup.bed
```

The reference can be a FASTA (optionally gzipped) or a pre-built minimap2 index (`.mmi`).
Use `-x/--preset` to pick the minimap2 preset (default `map-ont`, use e.g. `map-hifi` for PacBio HiFi reads).
The output is not sorted.
To skip compressing the intermediate file, write uncompressed BAM with `-u` and pipe it into `samtools sort`:

```bash
modkit modbam align -u -t 8 --ref reference.mmi calls.bam - | samtools sort -@ 4 -o aligned.sorted.bam
```

For large references, build the index once with `minimap2 -x map-ont -d reference.mmi reference.fa` (use the same preset) and pass the `.mmi`, otherwise the index is rebuilt on every run.

## Compared to `samtools fastq -T MM,ML | minimap2 -y`

The command does the same thing as the common pipeline, but avoids its pitfalls:

* **Supplementary alignments are soft-clipped.** By default minimap2 hard-clips supplementary alignments but copies the full-length `MM`/`ML` tags, so the tags refer to bases that are no longer in `SEQ` and the modification calls on those records can't be used.
  Every record written by `modkit modbam align` (primary, supplementary and, with `--secondary`, secondary) holds the full read sequence.
* **`MN` is set** to the sequence length on every record with modification tags, so downstream tools can check the tags apply to the sequence.
* **Old alignment tags are dropped** (`NM`, `MD`, `AS`, `SA`, `HP`, `PS`, minimap2 tags, ...) instead of being duplicated, and new `NM`, `AS`, `tp` and `SA` tags are written (`MD` with `--md`).
* **Legacy `Mm`/`Ml` tags are kept** (they're easy to drop by accident with `samtools fastq -T MM,ML`).
* **Already aligned input works.** Reverse-strand reads are converted back to their original orientation before alignment, and secondary/supplementary input records are skipped so each read is aligned once.

Reads that don't align are written as unmapped records (with their tags) unless `--no-unmapped` is passed.
Paired-end reads aren't supported.

## Aligning and piling up in one step

If the goal is a bedMethyl table, `modkit pileup --align` aligns the reads and piles them up without writing, sorting, indexing and re-reading an aligned BAM:

```bash
modkit pileup --align reference.mmi --align-preset map-ont -t 8 calls.bam pileup.bed
```

All the usual `pileup` options work (`--cpg`, `--motif`, `--modified-bases`, `--combine-strands`, thresholds, ...).
The output is identical to running `modkit modbam align`, `samtools sort` and `samtools index`, and then `modkit pileup` on the sorted BAM, including the estimated filter threshold.
The aligned records are kept in memory (sorted the same way as `samtools sort`), roughly 1.5 bytes per aligned base plus the MM/ML tags, so this suits bacterial genomes, metagenomes, targeted data and single chromosomes better than whole human genomes at high coverage.
No aligned BAM is produced; run `modkit modbam align` if you also need one.

## Building for speed

If you build modkit yourself for one machine, compiling the embedded minimap2 for that CPU makes alignment about 4% faster with identical output (`-ffp-contract=off` keeps floating-point results unchanged):

```bash
CFLAGS="-O3 -march=native -ffp-contract=off" cargo build --release
```

The resulting binary may not run on other CPUs.
