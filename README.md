### Description

This repository contains the scripts and documentation used to characterize long non-coding RNAs (lncRNAs) during compatible and incompatible pollination in *Arabidopsis thaliana*.

### Citation (bioRxiv)

**Nischay Patel**, **Nilesh D. Gawande**, **Subramanian Sankaranarayanan** *Characterization of long non-coding RNAs during compatible and incompatible pollination in Arabidopsis thaliana (2026)* [doi: 10.64898/2026.04.29.721561](https://doi.org/10.64898/2026.04.29.721561)

------------------------------------------------------------------------

### Repository Structure

The repository is organized into six numbered folders, one per pipeline stage. Each folder follows the same layout:

```         
0X_<Stage_Name>/
├── 0X_Inputs/
├── 0X_Outputs/
└── 0X_<Stage_Name>.md  (or .R)
```

- **`0X_Inputs/`** — the input file(s) required for that stage
- **`0X_Outputs/`** — the output file(s) produced by that stage
- **`0X_<Stage_Name>.md` / `.R`** — the script(s) or documentation for that stage, including tool versions, citations, input/output file names, installation notes, and (for the `.R` scripts) commented code

> **Important:** These scripts are written as sequential stages of a single pipeline, not standalone tools — the output of one stage is frequently the input of a later one. To reproduce the results reported in the paper above, follow the manuscript's Methods section from the beginning and read it carefully alongside these scripts, rather than running any single folder in isolation. The overall analysis was also carried out **cross-platform**; some steps were run on a **Galaxy server**, others locally in **RStudio** or the **terminal**, and a few coding-potential tools via their respective **public webservers**.

------------------------------------------------------------------------

### Pipeline Overview

#### 01 — Library Strandedness Assessment

Determines the strandedness (orientation) of paired-end RNA-Seq reads relative to the reference transcriptome, a required upstream QC step before running any strand-aware tool.

> `01_Library_Strandedness_Assessment.md`

------------------------------------------------------------------------

#### 02 — Coding Potential Prediction

Classifies filtered transcripts as coding or non-coding using three complementary tools (CPC2, CPAT, LncFinder), followed by a Pfam-A domain scan to remove any remaining transcripts with protein-coding signatures.

> `02_Coding_Potential_Prediction.md`

------------------------------------------------------------------------

#### 03 — Transcript/Gene ID Reconciliation

Reconciles transcript and gene identifiers between the TAIR10 (Araport11) reference annotation and the gffcompare-assembled transcriptome, producing a single consistent ID space used in downstream stages.

> `03_Transcript_Gene_Id_Reconciliation.R`

------------------------------------------------------------------------

#### 04 — FEELnc: Classification of lncRNAs by Genomic Location

Classifies identified lncRNAs by their genomic location relative to neighboring protein-coding genes (e.g. intergenic, antisense, intronic, sense-overlapping).

> `04_FEELnc.md`

------------------------------------------------------------------------

#### 05 — Transcript Expression Analysis

Builds per-timepoint TPM expression matrices from Kallisto quantification, defines replicate-aware "expressed" transcripts, identifies condition-specific transcripts, and extracts the corresponding lncRNA FASTA/GTF subsets.

> `05_Transcript_Expression_Analysis.R`

------------------------------------------------------------------------

#### 06 — Cis/Trans Target Prediction

Predicts candidate mRNA targets for differentially expressed lncRNAs via genomic proximity (cis, using BEDTools) and expression correlation (cis and trans, using Pearson/Spearman correlation in R).

> `06_Cis_Trans_Target_Prediction.md`

------------------------------------------------------------------------

Other technical details are provided inside each stage's own script/documentation file.

------------------------------------------------------------------------

### Repository Citation

Characterization of long non-coding RNAs during compatible and incompatible pollination in Arabidopsis thaliana [Internet]. [cited 2026 July 8]. Available from: <https://github.com/nischaypatel4/LncRNA-Analysis-Pollination-Arabidopsis>

### Contact

[nischaypatel4\@gmail.com](mailto:nischaypatel4@gmail.com){.email}
