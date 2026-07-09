# =============================================================================
# Transcript / Gene ID Reconciliation
# =============================================================================
#
# PURPOSE:
#   Reconciles transcript and gene identifiers between the TAIR10 reference
#   annotation and the gffcompare-assembled transcriptome, so that every
#   transcript in the final GTF carries a single, consistent ID: its
#   original reference ID where the assembled transcript corresponds to an
#   already-annotated transcript, or a new assembly-derived ID (TCONS/XLOC)
#   where it does not. This reconciled GTF is the reference GTF used by all
#   downstream transcript expression filtering and lncRNA analysis steps.

#   This script defines the ID space used everywhere downstream.
#   The join keys and the class_code-based decision logic in Sections 3–4
#   must not be altered without re-validating every downstream step that
#   consumes `final_mapping.tsv` or `Final_GTF_annotation.gtf`.
#
# -----------------------------------------------------------------------------
# INPUTS (exact file names, read in the order used below):
#   1) tair10_gtf.gtf          — TAIR10 (Araport11) reference annotation GTF
#   2) gffcompare.gtf          — gffcompare-assembled transcriptome GTF
#                                (used twice: Section 2 and Section 4)
#   3) trackingfile.tsv        — gffcompare .tracking file (no header row)
#
# INTERMEDIATE / FINAL OUTPUTS (exact file names, written in this order):
#   1) protein_coding_geneids_TAIR10.txt        — unique protein-coding gene IDs
#   2) protein_coding_transcriptids_TAIR10.txt  — unique protein-coding transcript IDs
#   3) lncRNA_geneids_TAIR10.txt                — unique reference lncRNA gene IDs
#   4) lncRNA_transcriptids_TAIR10.txt          — unique reference lncRNA transcript IDs
#   5) gff_mapping.tsv                          — 5-column map: transcript_id,
#                                                  gene_id, ref_gene_id, cmp_ref,
#                                                  class_code (from gffcompare.gtf)
#   6) final_mapping.tsv                        — reconciled ID map: t_id_prev,
#                                                  g_id_prev, t_id_new, g_id_new,
#                                                  exon_count, length, class_code
#   7) Final_GTF_annotation.gtf                   — gffcompare.gtf with gene_id /
#                                                  transcript_id attributes
#                                                  replaced according to
#                                                  final_mapping.tsv
#                                                  (this is the reference GTF
#                                                  used everywhere downstream)
#
# -----------------------------------------------------------------------------
# LOGIC OVERVIEW:
#
#   Section 1 — Reference annotation split
#     Reads the TAIR10 GTF, keeps only "transcript" features, and extracts
#     gene_id / transcript_id / transcript_biotype from the attribute column.
#     Splits reference transcripts into two sets by biotype
#     (protein_coding vs. lncRNA) and writes their unique gene/transcript IDs.
#     This defines the set of already-known/annotated genes and transcripts.
#
#   Section 2 — gffcompare attribute extraction
#     Reads the gffcompare-assembled GTF, keeps only "transcript" features,
#     and extracts five attributes per transcript: transcript_id, gene_id
#     (the assembled/query IDs), ref_gene_id and cmp_ref (the reference gene
#     the assembled transcript was compared against, if any), and class_code
#     (gffcompare's relationship code between the assembled transcript and
#     the reference — e.g. "=" exact match, "c" contained within a reference
#     transcript, or other codes for novel/partial/antisense relationships).
#
#   Section 3 — Final ID reconciliation
#     Joins the Section 2 map (`file1`) with the gffcompare tracking file
#     (`file2`), which supplies the assembler's own TCONS_/XLOC_ IDs and a
#     pipe-delimited info string per transcript. The info string is parsed
#     assuming the format "q1:<sample>|<query_transcript>|<exon_count>|...|<length>"
#     (a leading "q1:" prefix, sample name, exon count, and transcript length
#     as the last field); only the query transcript ID, exon count, and
#     length are retained. The two tables are joined on
#     transcript_id (file1) == query_transcript (file2).
#
#     Decision rule (must be preserved exactly): for each transcript,
#       - if class_code is "=" (exact match) or "c" (contained in a reference
#         transcript) → the assembled transcript corresponds to an
#         already-annotated transcript, so the ORIGINAL reference-derived
#         transcript_id and ref_gene_id are kept as the final IDs.
#       - for any other class_code → the assembled transcript does not
#         correspond to a known reference transcript, so the assembler's own
#         TCONS_ (transcript) and XLOC_ (gene) IDs are used as the final IDs.
#     The output retains both the previous and new IDs (t_id_prev/g_id_prev
#     vs. t_id_new/g_id_new) alongside exon_count, length, and class_code.
#
#   Section 4 — Apply reconciled IDs to the GTF
#     Re-reads gffcompare.gtf, extracts its current gene_id/transcript_id,
#     and left-joins against final_mapping.tsv on
#     (transcript_id, gene_id) == (t_id_prev, g_id_prev). For each row,
#     transcript_id_final/gene_id_final take the mapped new ID if a match
#     was found, otherwise fall back to the transcript's existing ID
#     (`coalesce`). The attribute column's gene_id and transcript_id values
#     are then replaced in place with these final IDs (transcript_id is only
#     replaced when a non-missing final value exists), leaving every other
#     attribute untouched. The result is written as `Final_GTF_annotation.gtf`.
#
# -----------------------------------------------------------------------------
# PACKAGES REQUIRED:
#   CRAN:         tidyverse, dplyr, data.table, stringr, Hmisc, reshape2,
#                 VennDiagram, gridExtra, grid, readxl, openxlsx
#   Bioconductor: DESeq2, tximport, GenomicFeatures, txdbmaker, Biostrings,
#                 SummarizedExperiment
#   (installed/loaded automatically in Section 0 below; several of these,
#   e.g. DESeq2/tximport/Biostrings, are not directly used in this
#   particular script but are kept for environment consistency with the
#   rest of the pipeline)
#
# NOTE ON PATHS:
#   All paths are built from a single `base_dir` variable defined below.
#   Update `base_dir` to point at your local project root before running;
#   no other path in the script should need manual editing. File names
#   themselves are left unchanged from the original pipeline.
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
# Root directory of the local project (update this to your local path).
# All paths below are built from this single variable.

base_dir <- "<path_to_project>/Analysis"


# -----------------------------------------------------------------------------
# 1) Splitting the TAIR10 reference into protein-coding vs. lncRNA sets
# -----------------------------------------------------------------------------

# Input: TAIR10 reference annotation GTF
gtf_file <- file.path(base_dir, "Data", "ReferenceFiles", "tair10_gtf.gtf")

# Read GTF (ignore comment lines starting with "#")
gtf <- read.table(gtf_file, sep = "\t", header = FALSE, quote = "",
                   comment.char = "#", stringsAsFactors = FALSE)

# Assign column names according to GTF format
colnames(gtf) <- c("seqname", "source", "feature", "start", "end", "score",
                    "strand", "frame", "attribute")

colnames(gtf)

# Filter only transcript features
transcripts <- gtf %>%
  filter(feature == "transcript")

# Extract attributes: gene_id, transcript_id, transcript_biotype
transcripts <- transcripts %>%
  mutate(
    gene_id = str_extract(attribute, 'gene_id "[^"]+"') %>% str_replace('gene_id "', '') %>% str_replace('"', ''),
    transcript_id = str_extract(attribute, 'transcript_id "[^"]+"') %>% str_replace('transcript_id "', '') %>% str_replace('"', ''),
    transcript_biotype = str_extract(attribute, 'transcript_biotype "[^"]+"') %>% str_replace('transcript_biotype "', '') %>% str_replace('"', '')
  )

# --- Protein coding ---
protein_coding_genes <- transcripts %>%
  filter(transcript_biotype == "protein_coding") %>%
  pull(gene_id) %>%
  unique()

protein_coding_transcripts <- transcripts %>%
  filter(transcript_biotype == "protein_coding") %>%
  pull(transcript_id) %>%
  unique()

# --- lncRNA ---
lncRNA_genes <- transcripts %>%
  filter(transcript_biotype == "lncRNA") %>%
  pull(gene_id) %>%
  unique()

lncRNA_transcripts <- transcripts %>%
  filter(transcript_biotype == "lncRNA") %>%
  pull(transcript_id) %>%
  unique()

# --- Write outputs ---
known_dir <- file.path(base_dir, "Data", "Known_transcripts")

write.table(protein_coding_genes,
            file.path(known_dir, "protein_coding_geneids_TAIR10.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(protein_coding_transcripts,
            file.path(known_dir, "protein_coding_transcriptids_TAIR10.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(lncRNA_genes,
            file.path(known_dir, "lncRNA_geneids_TAIR10.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(lncRNA_transcripts,
            file.path(known_dir, "lncRNA_transcriptids_TAIR10.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# --- Report ---
cat("Extracted", length(protein_coding_genes), "unique protein-coding gene IDs\n")
cat("Extracted", length(protein_coding_transcripts), "unique protein-coding transcript IDs\n")
cat("Extracted", length(lncRNA_genes), "unique lncRNA gene IDs\n")
cat("Extracted", length(lncRNA_transcripts), "unique lncRNA transcript IDs\n")


# -----------------------------------------------------------------------------
# 2) Building the gffcompare attribute map:
#    transcript_id, gene_id, ref_gene_id, cmp_ref, class_code
# -----------------------------------------------------------------------------

# Input: gffcompare-assembled transcriptome GTF
gtf_file <- file.path(base_dir, "Data", "GffCompare", "gffcompare.gtf")
# Output: 5-column attribute map
out_file <- file.path(base_dir, "Data", "Mapping", "gff_mapping.tsv")

# Read GTF (tab-delimited, skip comments)
gtf <- read.delim(gtf_file, header = FALSE, comment.char = "#")

# Filter only transcript features
gtf_tx <- gtf %>% filter(V3 == "transcript")

# Extract attribute function (robust regex)
extract_attr <- function(attr, key) {
  match <- str_match(attr, paste0(key, " ([^;]+);"))
  return(match[, 2])
}

map_df <- data.frame(
  transcript_id = extract_attr(gtf_tx$V9, "transcript_id"),
  gene_id       = extract_attr(gtf_tx$V9, "gene_id"),
  ref_gene_id   = extract_attr(gtf_tx$V9, "ref_gene_id"),
  cmp_ref       = extract_attr(gtf_tx$V9, "cmp_ref"),
  class_code    = extract_attr(gtf_tx$V9, "class_code"),
  stringsAsFactors = FALSE
)

# Write to file WITH header
write.table(map_df, out_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

cat("Clean 5-column map written to:", out_file, "\n")

# -----------------------------------------------------------------------------
# 3) Reconciling final transcript/gene IDs (reference IDs vs. assembler IDs)
# -----------------------------------------------------------------------------

# ----------- Input files -----------
# file1: the 5-column map produced in Section 2
file1 <- read_tsv(file.path(base_dir, "Data", "Mapping", "gff_mapping.tsv"),
                   col_types = "ccccc") # transcript_id, gene_id, ref_gene_id, cmp_ref, class_code

# file2: gffcompare's own .tracking file (no header row)
# Expected columns: TCONS_id, XLOC_id, geneName|cmp_ref, class_code, info_string
# where info_string is assumed to have the form
# "q1:<sample>|<query_transcript>|<exon_count>|...|<length>"
file2 <- read_tsv(file.path(base_dir, "Data", "GffCompare", "trackingfile.tsv"),
                   col_names = FALSE, col_types = "ccccc")

# Rename columns in file2 for clarity
colnames(file2) <- c("t_id_new_file2", "g_id_new_file2", "geneName_cmpRef", "class_code2", "info")

file2_parsed <- file2 %>%
  mutate(
    info_clean = str_remove(info, "^q1:"),
    info_parts = str_split(info_clean, "\\|")
  ) %>%
  mutate(
    query_transcript = vapply(info_parts, `[`, character(1), 2),
    exon_count       = as.integer(vapply(info_parts, `[`, character(1), 3)),
    length           = as.integer(vapply(
      info_parts,
      function(x) x[length(x)],
      character(1)
    ))
  ) %>%
  dplyr::select(
    t_id_new_file2,
    g_id_new_file2,
    query_transcript,
    class_code2,
    exon_count,
    length
  )

# ----------- Join with file1 -----------
final <- file1 %>%
  dplyr::left_join(file2_parsed, by = c("transcript_id" = "query_transcript")) %>%
  dplyr::mutate(
    t_id_prev = transcript_id,
    g_id_prev = gene_id,
    # ---- Final ID decision rule ----
    # class_code "=" (exact match) or "c" (contained within a reference
    # transcript) => this assembled transcript IS an already-annotated
    # transcript, so keep the reference-derived IDs.
    # Any other class_code => this is a novel/partial/unmatched transcript,
    # so use the assembler's own TCONS_/XLOC_ IDs instead.
    t_id_new = dplyr::case_when(
      class_code %in% c("=", "c") ~ t_id_prev,       # keep transcript_id
      TRUE                        ~ t_id_new_file2   # use TCONS_xxx
    ),
    g_id_new = dplyr::case_when(
      class_code %in% c("=", "c") ~ ref_gene_id,     # use ref_gene_id
      TRUE                        ~ g_id_new_file2   # use XLOC_xxx
    )
  ) %>%
  dplyr::select(
    t_id_prev,
    g_id_prev,
    t_id_new,
    g_id_new,
    exon_count,
    length,
    class_code
  )

# ----------- Write output -----------
write_tsv(final, file.path(base_dir, "Data", "Mapping", "final_mapping.tsv"))

# -----------------------------------------------------------------------------
# 4) Applying reconciled IDs to the assembled GTF
# -----------------------------------------------------------------------------
# (Originally labeled "9)" in the working script — renumbered here since this
#  is the fourth and final step of this script; see NOTE ON SECTION NUMBERING
#  above.)

gtf_file     <- file.path(base_dir, "Data", "GTF", "gffcompare.gtf")
mapping_file <- file.path(base_dir, "Data", "Mapping", "final_mapping.tsv")
output_file  <- file.path(base_dir, "Data", "GTF", "Final_GTF_annotation.gtf")

# ------------------ Step 1: Read the GTF ------------------
gtf <- read_tsv(
  gtf_file,
  comment = "#",
  col_names = c("seqname", "source", "feature", "start", "end",
                "score", "strand", "frame", "attribute"),
  col_types = cols(.default = col_character())
)

# Extract gene_id and transcript_id from attribute column
gtf <- gtf %>%
  mutate(
    gene_id = str_match(attribute, 'gene_id "([^"]+)"')[, 2],
    transcript_id = str_match(attribute, 'transcript_id "([^"]+)"')[, 2]
  )

# ------------------ Step 2: Read Mapping ------------------
mapping <- read_tsv(
  mapping_file,
  col_types = cols(
    t_id_prev = col_character(),
    g_id_prev = col_character(),
    t_id_new  = col_character(),
    g_id_new  = col_character()
  )
) %>%
  mutate(across(c(t_id_prev, g_id_prev, t_id_new, g_id_new),
                ~ str_replace_all(., '^"|"$', '')))

any(str_detect(mapping$t_id_new, '^".*"$'))

# ------------------ Step 3: Join and Replace IDs ------------------
gtf_updated <- gtf %>%
  left_join(mapping, by = c("transcript_id" = "t_id_prev", "gene_id" = "g_id_prev")) %>%
  mutate(
    transcript_id_final = coalesce(t_id_new, transcript_id),
    gene_id_final       = coalesce(g_id_new, gene_id)
  )

# ------------------ Step 4: Update attribute field ------------------
# Replace only the values of gene_id and transcript_id, leave everything else intact
attr1 <- str_replace(
  gtf_updated$attribute,
  'gene_id "[^"]+"',
  paste0('gene_id "', gtf_updated$gene_id_final, '"')
)

# Only replace transcript_id if it was originally present
attribute_updated <- ifelse(
  !is.na(gtf_updated$transcript_id_final),
  str_replace(attr1,
              'transcript_id "[^"]+"',
              paste0('transcript_id "', gtf_updated$transcript_id_final, '"')),
  attr1
)

# Add back updated attribute
gtf_final <- gtf_updated %>%
  mutate(attribute = attribute_updated) %>%
  dplyr::select(seqname, source, feature, start, end, score, strand, frame, attribute)

# ------------------ Step 5: Write updated GTF ------------------
write.table(gtf_final, output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

message("Updated GTF written to: ", output_file)
