## 01_Library_Strandedness_Assessment

### Purpose

Determine the strandedness (orientation) of paired-end RNA-Seq reads relative to the reference transcriptome. This is a required upstream QC step before running strand-aware tools (e.g., featureCounts, StringTie, HTSeq) with the correct `--stranded` / `-s` parameters.

### Tool

[`how_are_we_stranded_here`](https://github.com/signalbash/how_are_we_stranded_here) v1.0.1

> Signal, B., & Kahlke, T. (2022). *how_are_we_stranded_here: quick determination of RNA-Seq strandedness.* BMC Bioinformatics, 23(1), 49. <https://doi.org/10.1186/s12859-022-04572-7>

### Platform

Local (terminal), Python 3.13, installed via `pip`

### Inputs

| Input | Description |
|----|----|
| Reference GTF | Genome annotation file |
| Reference cDNA | Transcriptome FASTA (`.fa`) |
| Raw paired-end FASTQ | Forward (`R1`) and reverse (`R2`) reads, gzip-compressed |

### Output

Predicted read orientation / library type (e.g., unstranded, `fr-firststrand`/RF, or `fr-secondstrand`/FR), reported as fractions of reads explained by each strand configuration.

------------------------------------------------------------------------

### Installation

``` bash
pip install how_are_we_stranded_here
```

The package ships with a template script that needs to be copied into the installed package directory before use:

``` bash
cp check_strandedness.py <path_to_conda_or_venv>/lib/python3.13/site-packages/how_are_we_stranded_here/check_strandedness.py
```

> Replace `<path_to_conda_or_venv>` with the path to your active Python environment (find it with `python -c "import sys; print(sys.prefix)"`).

### Usage

``` bash
check_strandedness \
  --gtf <path_to_project>/StrandednessCheck/reference.gtf \
  --transcripts <path_to_project>/StrandednessCheck/Arabidopsis_thaliana.TAIR10.cdna.all.fa \
  --reads_1 <path_to_project>/StrandednessCheck/SRR7546174_forward.gz \
  --reads_2 <path_to_project>/StrandednessCheck/SRR7546174_reverse.gz
```

**Arguments**

- `--gtf` — reference annotation (GTF format)
- `--transcripts` — reference cDNA/transcriptome FASTA
- `--reads_1` / `--reads_2` — forward/reverse raw FASTQ files (gzip-compressed)

By default, the tool subsamples the first 200,000 reads, builds a `kallisto` index, pseudoaligns the subsampled reads, and infers strandedness from the resulting BAM.

> **Note:** Installation issues can vary depending on your OS, Python version, and environment setup. For system-specific installation problems, please refer to the tool's [GitHub repository](https://github.com/signalbash/how_are_we_stranded_here) or its associated [research paper](https://doi.org/10.1186/s12859-022-04572-7).

------------------------------------------------------------------------

### Sample Output

``` text
converting gtf to bed
using kallisto_index as kallisto index
creating fastq files with first 200000 reads
quantifying with kallisto
[quant] fragment length distribution will be estimated from the data
[index] k-mer length: 31
[index] number of targets: 48,359
[index] number of k-mers: 45,567,187
[index] number of equivalence classes: 93,973
Warning: 5654 transcripts were defined in GTF file, but not in the index
[quant] running in paired-end mode
[quant] will process pair 1: stranded_test_forward_2/forward_sample.fq
                             stranded_test_forward_2/reverse_sample.fq
[quant] finding pseudoalignments for the reads ... done
[quant] processed 200,000 reads, 192,292 reads pseudoaligned
[quant] estimated average fragment length: 180.505
[   em] quantifying the abundances ... done
[   em] the Expectation-Maximization algorithm ran for 853 rounds
[  bam] writing pseudoalignments to BAM format .. done
[  bam] sorting BAM files .. done
[  bam] indexing BAM file .. done
checking strandedness
Reading reference gene model stranded_test_forward_2/reference.bed ... Done
Loading SAM/BAM file ...  Total 200000 usable reads were sampled
This is PairEnd Data
Fraction of reads failed to determine: 0.0337
Fraction of reads explained by "1++,1--,2+-,2-+": 0.0037 (0.4% of explainable reads)
Fraction of reads explained by "1+-,1-+,2++,2--": 0.9626 (99.6% of explainable reads)
Over 90% of reads explained by "1+-,1-+,2++,2--"
Data is likely RF/fr-firststrand
```

### Interpretation

Over 90% of explainable reads matched the `"1+-,1-+,2++,2--"` configuration, indicating the library is **RF / `fr-firststrand`**. This strandedness setting should be used consistently in all downstream strand-aware steps (e.g., `kallisto quant --rf-stranded`, `featureCounts -s 2`, `STAR --outSAMstrandField`).

------------------------------------------------------------------------
