# =============================================================================
# 05_Transcript_Expression_Analysis_and_Filtering
# =============================================================================
#
# PURPOSE:
#   From per-sample Kallisto TPM quantification files, this script:
#     1) Builds a merged TPM matrix across all replicates/conditions
#     2) Defines "expressed" transcripts using a replicate-aware presence rule
#     3) Identifies condition-specific (Compatible-only / Incompatible-only)
#        transcripts
#     4) Subsets the expressed lncRNA set (all / annotated / novel)
#     5) Writes lncRNA-filtered TSV, FASTA, and GTF outputs for downstream use
#
# TIMEPOINT GENERALIZATION:
#   This script is written for a single timepoint (default: t = 10 min).
#   To reuse for another timepoint (e.g. t = 60 min), only the variables in
#   the "USER-CONFIGURABLE PARAMETERS" section below need to change:
#     - `timepoint`      -> e.g. "60"
#     - `samples`        -> the corresponding sample/replicate file names for
#                           that timepoint (Compatible replicates first,
#                           Incompatible replicates second)
#   All output file names automatically inherit the `timepoint` prefix
#   (e.g. "10_expressed_TPM_matrix.tsv" -> "60_expressed_TPM_matrix.tsv"),
#   so the existing naming convention is preserved without manual renaming.
#
# INPUTS REQUIRED (per timepoint):
#   - Per-sample Kallisto `abundance.tsv`-style files, one per replicate,
#     each containing at minimum `target_id` and `tpm` columns
#   - Transcript ID lists: lncRNAs_all.txt,
#     lncRNAs_annot.txt, lncRNAs_novel.txt (from the genomic classification
#     step; see 03_lncRNA_Genomic_Classification.md)
#   - Reference transcriptome FASTA (all assembled transcript sequences)
#   - Merged/annotated GTF : "Final_GTF_annotation" containing transcript_id
#     attributes for all transcripts
#
# PACKAGES REQUIRED:
#   CRAN:         tidyverse, dplyr, data.table, stringr, Hmisc, reshape2,
#                 VennDiagram, gridExtra, grid, readxl, openxlsx
#   Bioconductor: DESeq2, tximport, GenomicFeatures, txdbmaker, Biostrings,
#                 SummarizedExperiment
#   (installed/loaded automatically in Section 0 below)
#
# NOTE ON PATHS:
#   All paths are built from a single `base_dir` variable defined below.
#   Update `base_dir` to point at your local project root before running;
#   no other path in the script should need manual editing.
# =============================================================================


# -----------------------------------------------------------------------------
# 0) Setup: Install & Load Packages
# -----------------------------------------------------------------------------

# ---- CRAN packages ----
cran_pkgs <- c(
  "tidyverse",
  "dplyr",
  "data.table",
  "stringr",
  "Hmisc",
  "reshape2",
  "VennDiagram",
  "gridExtra",
  "grid",
  "readxl",
  "openxlsx"
)

installed_cran <- rownames(installed.packages())
for (pkg in cran_pkgs) {
  if (!pkg %in% installed_cran) {
    install.packages(pkg, dependencies = TRUE)
  }
}

# ---- Bioconductor packages ----
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_pkgs <- c(
  "DESeq2",
  "tximport",
  "GenomicFeatures",
  "txdbmaker",
  "Biostrings",
  "SummarizedExperiment"
)

installed_bioc <- rownames(installed.packages())
for (pkg in bioc_pkgs) {
  if (!pkg %in% installed_bioc) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

# -----------------------------------------------------------------------------
# Load libraries
# -----------------------------------------------------------------------------

library(tidyverse)
library(dplyr)
library(data.table)
library(stringr)
library(Hmisc)
library(reshape2)
library(VennDiagram)
library(gridExtra)
library(grid)
library(readxl)
library(openxlsx)
library(DESeq2)
library(tximport)
library(GenomicFeatures)
library(txdbmaker)
library(Biostrings)
library(SummarizedExperiment)


# -----------------------------------------------------------------------------
# USER-CONFIGURABLE PARAMETERS
# -----------------------------------------------------------------------------
# Update these values to point at your local project and the timepoint being
# processed. Everything below this block references these variables only.

# Root directory of the local project (update this to your local path)
base_dir <- "<path_to_project>/Analysis"

# Timepoint label used as a prefix for all output file names
# (e.g. "10" for t=10 min; use "60" for t=60 min, etc.)
timepoint <- "10"

# Replicate/sample names for this timepoint.
# IMPORTANT: order matters — Compatible replicates must come first,
# followed by Incompatible replicates, since indices below (1:3 / 4:6)
# assume this ordering. Adjust the split indices if your replicate counts
# differ per condition.
samples <- c("10_compatible1", "10_compatible2", "10_compatible4",
             "10_incompatible1", "10_incompatible2", "10_incompatible3")


# -----------------------------------------------------------------------------
# 7) TPM Matrix Construction and Expressed Transcripts
# -----------------------------------------------------------------------------

# ------------------ Define per-sample Kallisto TPM file paths ------------------
paths <- file.path(base_dir, "Data", paste0("KallistoCounts_", timepoint),
                    "Original", paste0(samples, ".tsv"))

# ------------------ Read and merge TPMs across all replicates ------------------
tpm_list <- lapply(seq_along(paths), function(i){
  df <- read.delim(paths[i], stringsAsFactors = FALSE)[, c("target_id", "tpm")]
  colnames(df) <- c("target_id", samples[i])
  df
})

tpm_mat <- Reduce(function(a, b) full_join(a, b, by = "target_id"), tpm_list)
cat(">> TPM matrix dims before filtering:", dim(tpm_mat), "\n")

# Replace NA with 0 (important for replicate counting — a transcript missing
# from a sample's Kallisto output is treated as unexpressed, not missing data)
tpm_mat[, samples][is.na(tpm_mat[, samples])] <- 0

# ------------------ Build the numeric expression matrix ------------------
expr <- as.matrix(tpm_mat[, samples])
rownames(expr) <- tpm_mat$target_id

# ------------------ Replicate-aware expression calls ------------------
# A transcript is called "expressed" in a condition if it exceeds the TPM
# threshold (0.05) in at least 2 of that condition's 3 replicates.

# Compatible condition (samples 1–3)
comp_ge2 <- rowSums(expr[, 1:3] > 0.05) >= 2

# Incompatible condition (samples 4–6)
ic_ge2   <- rowSums(expr[, 4:6] > 0.05) >= 2

# ------------------ Expressed transcript definition ------------------
# A transcript is retained as "expressed" if it meets the >=2/3 replicate
# threshold in EITHER condition (Compatible OR Incompatible).
expressed_ids <- rownames(expr)[comp_ge2 | ic_ge2]

cat(">> Expressed transcripts at t=", timepoint,
    " (>=2/3 in either condition):", length(expressed_ids), "\n")

# Subset the full TPM matrix down to expressed transcripts only
tpm_expr <- tpm_mat[tpm_mat$target_id %in% expressed_ids, ]

cat(">> TPM matrix dims after expression filtering:",
    dim(tpm_expr), "\n\n")

out_dir <- file.path(base_dir, "Data", paste0("KallistoCounts_", timepoint),
                      "Expressed_01")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ------------------ Write expressed TPM matrix (all samples) ------------------
expressed_tpm_path <- file.path(
  out_dir,
  paste0(timepoint, "_expressed_TPM_matrix.tsv")
)

expressed_tpm_df <- tpm_mat[
  tpm_mat$target_id %in% expressed_ids,
  c("target_id", samples),
  drop = FALSE
]

expressed_tpm_df <- expressed_tpm_df[order(expressed_tpm_df$target_id), ]

readr::write_tsv(expressed_tpm_df, expressed_tpm_path)

cat(">> Written expressed TPM matrix (all 6 samples):",
    expressed_tpm_path, "\n")

# ------------------ Condition-specific transcripts ------------------
# A transcript is "Compatible-specific" if it is expressed (TPM > 0.1) in
# 2 or 3 of the Compatible replicates AND completely absent (TPM == 0) in
# all Incompatible replicates (and vice versa for Incompatible-specific).

# Compatible-specific
comp_spec_2 <- rownames(expr)[
  rowSums(expr[, 1:3] > 0.1) == 2 &
    rowSums(expr[, 4:6] > 0) == 0
]

comp_spec_3 <- rownames(expr)[
  rowSums(expr[, 1:3] > 0.1) == 3 &
    rowSums(expr[, 4:6] > 0) == 0
]

comp_spec_ids <- union(comp_spec_2, comp_spec_3)

# Incompatible-specific
ic_spec_2 <- rownames(expr)[
  rowSums(expr[, 4:6] > 0.1) == 2 &
    rowSums(expr[, 1:3] > 0) == 0
]

ic_spec_3 <- rownames(expr)[
  rowSums(expr[, 4:6] > 0.1) == 3 &
    rowSums(expr[, 1:3] > 0) == 0
]

ic_spec_ids <- union(ic_spec_2, ic_spec_3)

# ------------------ Summary counts ------------------
cat(">> Compatible-specific (2/3):", length(comp_spec_2), "\n")
cat(">> Compatible-specific (3/3):", length(comp_spec_3), "\n")
cat(">> Incompatible-specific (2/3):", length(ic_spec_2), "\n")
cat(">> Incompatible-specific (3/3):", length(ic_spec_3), "\n\n")

# -----------------------------------------------------------------------------
# Output: expressed transcript IDs and condition-specific TPM matrices
# -----------------------------------------------------------------------------

# ------------------ Save expressed transcript IDs ------------------
writeLines(expressed_ids,
           file.path(out_dir, paste0(timepoint, "_expressed_transcript_ids.txt")))

# ------------------ Helper function: subset a TPM matrix to a set of IDs ------------------
write_subset_matrix <- function(ids, mat, samples, out_path) {
  subset_df <- mat[mat$target_id %in% ids, c("target_id", samples), drop = FALSE]
  subset_df <- subset_df[order(subset_df$target_id), ]
  readr::write_tsv(subset_df, out_path)
  cat(">> Wrote", nrow(subset_df), "rows to:", out_path, "\n")
}

# ------------------ Write full expressed TPM matrix ------------------
write_subset_matrix(
  expressed_ids, tpm_mat, samples,
  file.path(out_dir, paste0(timepoint, "_expressed_transcripts.tsv"))
)

# ------------------ Write condition-specific TPM matrices ------------------
spec_dir <- file.path(out_dir, "Specific")
if (!dir.exists(spec_dir)) dir.create(spec_dir)

write_subset_matrix(
  comp_spec_ids, tpm_mat, samples,
  file.path(spec_dir, paste0(timepoint, "_comp_spec_transcripts.tsv"))
)

write_subset_matrix(
  ic_spec_ids, tpm_mat, samples,
  file.path(spec_dir, paste0(timepoint, "_incomp_spec_transcripts.tsv"))
)

# ------------------ Write filtered per-sample Kallisto outputs ------------------
# For each original per-sample file, retain only expressed transcripts and
# write a companion "_expressed.tsv" file alongside the expressed matrices.
for (i in seq_along(paths)) {
  df <- read.delim(paths[i], stringsAsFactors = FALSE)
  df_filt <- df[df$target_id %in% expressed_ids, ]
  out_path <- file.path(
    out_dir,
    sub(".tsv$", "_expressed.tsv", basename(paths[i]))
  )
  write.table(df_filt, out_path, sep = "\t",
              quote = FALSE, row.names = FALSE)
}

cat("Expression filtering complete.\n")


# -----------------------------------------------------------------------------
# 8) Expressed lncRNAs at this timepoint
# -----------------------------------------------------------------------------
# Intersects the expressed transcript ID set (from Section 7) with the
# lncRNA IDs lists (all / annotated / novel) to
# identify which lncRNAs are actually expressed at this timepoint.

# Directories
expr_dir <- out_dir
lnc_result_dir <- file.path(base_dir, "Results", "Intersection_raw")

# Input files
expressed_file <- file.path(expr_dir, paste0(timepoint, "_expressed_transcript_ids.txt"))

lnc_all_file   <- file.path(lnc_result_dir, "lncRNAs_all.txt")
lnc_annot_file <- file.path(lnc_result_dir, "lncRNAs_annot.txt")
lnc_novel_file <- file.path(lnc_result_dir, "lncRNAs_novel.txt")

# Output files (same directory)
out_all   <- file.path(lnc_result_dir, paste0(timepoint, "_lncRNAs_all_expressed.txt"))
out_annot <- file.path(lnc_result_dir, paste0(timepoint, "_lncRNAs_annot_expressed.txt"))
out_novel <- file.path(lnc_result_dir, paste0(timepoint, "_lncRNAs_novel_expressed.txt"))

expressed_ids <- readLines(expressed_file)

lnc_all   <- readLines(lnc_all_file)
lnc_annot <- readLines(lnc_annot_file)
lnc_novel <- readLines(lnc_novel_file)

lnc_all_expr   <- intersect(lnc_all, expressed_ids)
lnc_annot_expr <- intersect(lnc_annot, expressed_ids)
lnc_novel_expr <- intersect(lnc_novel, expressed_ids)

writeLines(sort(lnc_all_expr),   out_all)
writeLines(sort(lnc_annot_expr), out_annot)
writeLines(sort(lnc_novel_expr), out_novel)

cat(">> Expressed lncRNAs (all):",   length(lnc_all_expr),   "\n")
cat(">> Expressed lncRNAs (annot):", length(lnc_annot_expr), "\n")
cat(">> Expressed lncRNAs (novel):", length(lnc_novel_expr), "\n")

cat("All expressed lncRNA files written to:\n", lnc_result_dir, "\n")


# -----------------------------------------------------------------------------
# 9) LncRNAs by class — write the common and condition-specific transcript TSVs
# -----------------------------------------------------------------------------
# Subsets the condition-specific and common expressed TPM matrices (Section 7)
# down to lncRNA transcripts only, using the expressed-lncRNA ID list from
# Section 8.

# ---------------------- Input Directories ----------------------
spec_dir   <- file.path(out_dir, "Specific")
common_dir <- out_dir
lnc_dir    <- lnc_result_dir

# ---------------------- Output Directories ----------------------
spec_out_dir   <- file.path(out_dir, "Specific_LncRNAs")
common_out_dir <- out_dir

# Create directories if they do not already exist
if (!dir.exists(spec_out_dir)) dir.create(spec_out_dir, recursive = TRUE)
if (!dir.exists(common_out_dir)) dir.create(common_out_dir, recursive = TRUE)

# ---------------------- Input LncRNA List ----------------------
lnc_list_file <- file.path(lnc_dir, paste0(timepoint, "_lncRNAs_all_expressed.txt"))

lnc_ids <- read_lines(lnc_list_file)
lnc_ids <- unique(str_trim(lnc_ids))
lnc_ids <- lnc_ids[nzchar(lnc_ids)]  # drop empty lines
cat(">> Loaded", length(lnc_ids), "LncRNA IDs\n\n")

# ---------------------- Mapping of input files to output files ----------------------
io_spec <- setNames(
  c(paste0(timepoint, "_comp_spec_LncRNAs.tsv"),
    paste0(timepoint, "_incomp_spec_LncRNAs.tsv")),
  c(paste0(timepoint, "_comp_spec_transcripts.tsv"),
    paste0(timepoint, "_incomp_spec_transcripts.tsv"))
)

io_common <- setNames(
  paste0(timepoint, "_expressed_TPM_matrix_LncRNAs.tsv"),
  paste0(timepoint, "_expressed_TPM_matrix.tsv")
)

# ---------------------- Helper Function ----------------------
# Reads an input TSV, filters rows to lncRNA transcript IDs, and writes the
# result to the corresponding output path.
subset_to_lnc <- function(in_dir, out_dir, in_name, out_name) {
  in_path  <- file.path(in_dir, in_name)
  out_path <- file.path(out_dir, out_name)

  if (!file.exists(in_path)) {
    message(">> Skipping (not found): ", in_path)
    return(invisible(NULL))
  }

  df <- readr::read_tsv(in_path, col_types = cols())
  if (!"target_id" %in% names(df)) {
    stop("File ", in_path, " has no 'target_id' column.")
  }

  df_sub <- df %>%
    filter(target_id %in% lnc_ids) %>%
    arrange(target_id)

  readr::write_tsv(df_sub, out_path)
  message(">> Wrote ", nrow(df_sub), " LncRNA rows to: ", out_path)
}

# ---------------------- Process condition-specific files ----------------------
for (i in seq_along(io_spec)) {
  subset_to_lnc(spec_dir, spec_out_dir, names(io_spec)[i], io_spec[i])
}

# ---------------------- Process common (all-sample) file ----------------------
for (i in seq_along(io_common)) {
  subset_to_lnc(common_dir, common_out_dir, names(io_common)[i], io_common[i])
}

message("\nAll LncRNA-filtered TSVs generated successfully!")


# -----------------------------------------------------------------------------
# 10) LncRNA FASTA extraction
# -----------------------------------------------------------------------------
# Extracts sequences for the expressed lncRNA set (Section 8) from the full
# assembled transcriptome FASTA, using Biostrings.

# Reference transcriptome FASTA (fixed reference, not timepoint-specific;
# update path if your assembly file is named/located differently)
transcriptome_fasta <- file.path(base_dir, "Data", "Fasta", "transcriptome_200bp_final.fasta")

lnc_ids_file  <- file.path(lnc_result_dir, paste0(timepoint, "_lncRNAs_all_expressed.txt"))
lnc_fasta_out <- file.path(base_dir, "Data", "Fasta",
                            paste0(timepoint, "_lncRNAs_all_expressed.fasta"))

keep_ids <- trimws(readLines(lnc_ids_file))
all_seqs <- readDNAStringSet(transcriptome_fasta, format = "fasta")
head(names(all_seqs))

# Strip any FASTA header description after the first whitespace so IDs match
seq_names_simple <- sub("\\s.*$", "", names(all_seqs))
keep_idx <- seq_names_simple %in% keep_ids
filtered_seqs <- all_seqs[keep_idx]

writeXStringSet(filtered_seqs, lnc_fasta_out, format = "fasta")
message("Kept ", sum(keep_idx), " / ", length(all_seqs), " sequences.")


# -----------------------------------------------------------------------------
# 11) LncRNA GTF extraction
# -----------------------------------------------------------------------------
# Filters the merged/annotated GTF down to only the expressed lncRNA
# transcripts (Section 8), by matching on the `transcript_id` attribute.

# Reference merged GTF (e.g. Final_GTF_annotation.gtf)
gtf_file_1 <- file.path(base_dir, "Data", "GTF", "Final_GTF_annotation.gtf")
tx_list    <- file.path(lnc_result_dir, paste0(timepoint, "_lncRNAs_all_expressed.txt"))
output_gtf <- file.path(base_dir, "Results", paste0("t", timepoint),
                         paste0(timepoint, "_LncRNAs_all_expressed.gtf"))

# assuming one transcript_id per line (no header)
tx_ids <- read_tsv(tx_list, col_names = FALSE)$X1

gtf <- read_tsv(
  gtf_file_1,
  comment = "#",
  col_names = c("seqname", "source", "feature", "start", "end", "score",
                "strand", "frame", "attribute"),
  col_types = cols(.default = col_character())
)

gtf <- gtf %>%
  mutate(transcript_id = str_match(attribute, 'transcript_id "([^"]+)"')[, 2])

gtf_filtered <- gtf %>%
  filter(transcript_id %in% tx_ids)

write.table(gtf_filtered,
            output_gtf,
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

message("Filtered GTF written to: ", output_gtf)
