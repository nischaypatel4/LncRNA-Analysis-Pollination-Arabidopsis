## Cis / Trans lncRNA Target Prediction

### Purpose

Predict candidate mRNA targets for differentially expressed (DE) lncRNAs using two complementary strategies:

- **Cis targets** — protein-coding transcripts located near each DE lncRNA in the genome (identified via `bedtools closest`), restricted to those whose target mRNA is itself DE, then filtered by expression sparsity and tested for expression correlation.
- **Trans targets** — all DE protein-coding transcripts genome-wide (excluding a lncRNA's own cis/overlap partners), tested for expression correlation with each DE lncRNA.

This pipeline is run per timepoint. All commands and scripts below are shown for **t = 10 min**; the identical commands/scripts are used for **t = 60 min** by substituting the `10_` filename prefix with `60_` throughout (input/output file names, directory names, and the `samples` vector in the R scripts).

------------------------------------------------------------------------

## Part A — Genomic Proximity: Cis-Target Identification (Bash / BEDTools)

### Tool

[BEDTools](https://bedtools.readthedocs.io/) — `closest` subcommand Version: *(not specified — confirm your installed version with `bedtools --version`)*

> Quinlan, A.R. & Hall, I.M. (2010). *BEDTools: a flexible suite of utilities for comparing genomic features.* Bioinformatics, 26(6), 841–842. <https://doi.org/10.1093/bioinformatics/btq033>

### Platform

Local (terminal)

### Inputs

| File | Description |
|------------------------------------|------------------------------------|
| `10_DE_LncRNAs_BED12.bed12` | BED12 of DE lncRNA transcripts at t=10 (t=60 equivalent: `60_DE_LncRNAs_BED12.bed12`) |
| `pctranscripts.bed12` | BED12 of all protein-coding transcripts (reference; shared across timepoints) |

### Outputs (in order produced)

| File | Description |
|------------------------------------|------------------------------------|
| `10_DE_LncRNAs.sorted_BED12.bed12` | Coordinate-sorted DE lncRNA BED12 |
| `pctranscripts.sorted.bed12` | Coordinate-sorted protein-coding transcript BED12 |
| `10_DE_LncRNAs_upstream.bed` | Closest protein-coding transcript per lncRNA, restricted to one directional side (see note below) |
| `10_DE_LncRNAs_downstream.bed` | Closest protein-coding transcript per lncRNA, restricted to the other directional side (see note below) |
| `10_DE_LncRNAs.cis_targets.bed` | Upstream + downstream hits merged and sorted |
| `10_DE_LncRNA.cis_targets.filtered.bed` | Above, restricted to ±100 kb of the lncRNA |
| `10_DE_LncRNAs_overlap.bed` | Protein-coding transcripts directly overlapping (distance = 0) a lncRNA |
| `10_DE_LncRNAs.cis_targets_final.bed` | Distance-filtered + overlap hits combined and sorted |
| `10_cis.tsv` | Final 11-column cis-pair table (see column layout below) |
| `all_pc_ids.txt` | Sorted, unique list of all protein-coding transcript IDs |

**`10_cis.tsv` column layout:** `chromA, startA, endA, lnc_id, strandA, chromB, startB, endB, pc_id, strandB, dist`

> **Note on `-iu`/`-id` and the upstream/downstream file names:** per the BEDTools manual, `-iu` restricts `closest` to the **downstream** feature (it *ignores upstream* features), and `-id` restricts it to the **upstream** feature (it *ignores downstream* features). In the commands below, the file named `_upstream.bed` is produced using `-iu` and the file named `_downstream.bed` is produced using `-id` — i.e. the flag-to-filename mapping appears reversed relative to BEDTools' own terminology. This does not affect the final merged/filtered cis-target output (both directions are combined before filtering), but if you ever need the upstream and downstream hits *individually*, double-check which file actually corresponds to which direction before relying on the name alone. The code is left exactly as originally run.

### Script

``` bash
# ---------------------------------------------------------------------------
# Cis-target identification via genomic proximity (BEDTools closest)
# Shown for t=10; identical for t=60 with the "10_" prefix replaced by "60_"
# ---------------------------------------------------------------------------

# Quick sanity check of the two input BED12 files
head -n 5 10_DE_LncRNAs_BED12.bed12
head -n 5 pctranscripts.bed12
# ----------------------------------------------------------------------------

# BEDTools closest requires coordinate-sorted input (chrom, then start)
sort -k1,1 -k2,2n 10_DE_LncRNAs_BED12.bed12   > 10_DE_LncRNAs.sorted_BED12.bed12
sort -k1,1 -k2,2n pctranscripts.bed12         > pctranscripts.sorted.bed12
# ----------------------------------------------------------------------------

# Sanity check the sorted files
head -n 5 10_DE_LncRNAs.sorted_BED12.bed12
head -n 5 pctranscripts.sorted.bed12
# ----------------------------------------------------------------------------

# Closest protein-coding transcript per lncRNA, direction restricted via -iu
# (-iu: ignore upstream features -> keeps the downstream side; see note above)
# -io: ignore any protein-coding feature that overlaps the lncRNA
# -D a: report signed distance relative to feature A (the lncRNA)
# -k 10: report up to the 10 closest hits per lncRNA
bedtools closest \
-a 10_DE_LncRNAs.sorted_BED12.bed12 \
-b pctranscripts.sorted.bed12 \
-iu \
-io \
-D a \
-k 10 \
> 10_DE_LncRNAs_upstream.bed
# ----------------------------------------------------------------------------

# Closest protein-coding transcript per lncRNA, direction restricted via -id
# (-id: ignore downstream features -> keeps the upstream side; see note above)
bedtools closest \
-a 10_DE_LncRNAs.sorted_BED12.bed12 \
-b pctranscripts.sorted.bed12 \
-id \
-io \
-D a \
-k 10 \
> 10_DE_LncRNAs_downstream.bed
# ----------------------------------------------------------------------------

# Sanity check both directional outputs, and confirm expected column count
head -n 5 10_DE_LncRNAs_upstream.bed
head -n 5 10_DE_LncRNAs_downstream.bed
head -n 1 10_DE_LncRNAs_upstream.bed   | awk '{print NF}'
head -n 1 10_DE_LncRNAs_downstream.bed | awk '{print NF}'
# ----------------------------------------------------------------------------

# Merge both directional hit sets and re-sort by chromosome + signed distance
# (column 25 = 12 BED12 fields for A + 12 BED12 fields for B + 1 distance field)
cat 10_DE_LncRNAs_upstream.bed 10_DE_LncRNAs_downstream.bed \
| sort -k1,1 -k25,25n \
> 10_DE_LncRNAs.cis_targets.bed
# ----------------------------------------------------------------------------

# How many raw (unfiltered) cis candidate pairs were found
wc -l 10_DE_LncRNAs.cis_targets.bed
# ----------------------------------------------------------------------------

# Restrict candidate pairs to within +/- 100 kb of the lncRNA
awk '($NF >= -100000 && $NF <= 100000)' 10_DE_LncRNAs.cis_targets.bed > 10_DE_LncRNA.cis_targets.filtered.bed
# ----------------------------------------------------------------------------

# How many pairs remain after the distance filter
wc -l 10_DE_LncRNA.cis_targets.filtered.bed
# ----------------------------------------------------------------------------

# Separately identify directly overlapping pairs (distance == 0), searching
# in both directions this time (no -iu/-id restriction, still excluding
# overlaps at the closest-feature search level via nothing extra here since
# -io is intentionally omitted so that overlapping pairs ARE captured)
bedtools closest \
-a 10_DE_LncRNAs.sorted_BED12.bed12 \
-b pctranscripts.sorted.bed12 \
-D a \
-k 10 \
| awk '$NF == 0' \
> 10_DE_LncRNAs_overlap.bed
# ----------------------------------------------------------------------------

# Combine the distance-filtered (+/-100kb) pairs with the overlap pairs
cat 10_DE_LncRNA.cis_targets.filtered.bed 10_DE_LncRNAs_overlap.bed \
| sort -k1,1 -k25,25n \
> 10_DE_LncRNAs.cis_targets_final.bed
# ----------------------------------------------------------------------------

# Extract the final 11 columns needed downstream:
# $1,$2,$3   = chromA, startA, endA (lncRNA coordinates)
# $4         = lnc_id (BED12 "name" field of feature A)
# $6         = strandA
# $13,$14,$15 = chromB, startB, endB (protein-coding transcript coordinates)
# $16        = pc_id (BED12 "name" field of feature B)
# $18        = strandB
# $25        = signed distance
awk 'BEGIN{OFS="\t"} { print $1,$2,$3,$4,$6,$13,$14,$15,$16,$18,$25 }' \
 10_DE_LncRNAs.cis_targets_final.bed \
> 10_cis.tsv
# ----------------------------------------------------------------------------

# List of all unique protein-coding transcript IDs (used later to build the
# trans candidate pair set)
cut -f4 pctranscripts.sorted.bed12 | sort -u > all_pc_ids.txt
```

------------------------------------------------------------------------

## Part B — Cis Pair Filtering & Correlation (R)

### Purpose

Restricts the raw cis pairs (`10_cis.tsv`) to those whose protein-coding partner is itself DE, then tests lncRNA–mRNA expression correlation (Pearson + Spearman) for pairs passing a per-transcript TPM sparsity filter.

### Inputs

| File | Description |
|------------------------------------|------------------------------------|
| `10_cis.tsv` | Final cis-pair table produced in Part A |
| `10_DE_mRNAs_transid.txt` | DE protein-coding transcript IDs (produced by the differential expression step; not covered in this document) |
| `10_compatible1_expressed.tsv`, `10_compatible2_expressed.tsv`, `10_compatible4_expressed.tsv`, `10_incompatible1_expressed.tsv`, `10_incompatible2_expressed.tsv`, `10_incompatible3_expressed.tsv` | Per-replicate expressed-transcript Kallisto TPM files (from the transcript expression filtering step) |

### Outputs

| File | Description |
|------------------------------------|------------------------------------|
| `10_cis_DE.tsv` | Cis pairs retained after requiring the mRNA partner to be DE |
| `10_cis_results.tsv` | Correlation results (Pearson + Spearman, with BH-adjusted p-values) for cis pairs passing the sparsity filter |

### Script

``` r
# =============================================================================
# Cis Pair Filtering & Expression Correlation (t=10 min)
# =============================================================================
# Identical for t=60: replace all "10_" prefixes below (file names, sample
# names, and directory names) with "60_".
#
# PATH NOTE: `base_dir` and `timepoint` are the only values that should need
# editing to point this at your local project; original hardcoded personal
# paths have been replaced with paths built from these two variables. All
# other logic and file names are unchanged from the original script.
# =============================================================================

base_dir  <- "<path_to_project>/Analysis"
timepoint <- "10"

cat("\n--- STEP 1: LOAD AND FILTER CIS PAIRS ---\n")

# Read the cis interaction table produced in Part A (bedtools closest pipeline)
cis <- read.table(file.path(base_dir, "Results", paste0("t", timepoint),
                             paste0(timepoint, "_cis"), paste0(timepoint, "_cis.tsv")),
                   header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(cis) <- c("chromA","startA","endA","lnc_id","strandA",
                   "chromB","startB","endB","pc_id","strandB","dist")

initial_pairs <- nrow(cis)
cat("Total raw cis pairs initially loaded  :", initial_pairs, "\n")

# Read the DE mRNA transcript IDs (produced by the DE analysis step)
de_mrna <- read.table(file.path(base_dir, "Results", paste0("t", timepoint),
                                 paste0(timepoint, "_DE_mRNAs_transid.txt")),
                       header = FALSE, stringsAsFactors = FALSE)

# Keep only rows where the mRNA transcript (pc_id) is in the DE list
cis_filtered <- cis[cis$pc_id %in% de_mrna[[1]], ]
retained_de_pairs <- nrow(cis_filtered)

# Write filtered pairs
out_dir <- file.path(base_dir, "Results", paste0("t", timepoint),
                      paste0(timepoint, "_cis"), "New")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

write.table(cis_filtered,
            file = file.path(out_dir, paste0(timepoint, "_cis_DE.tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


# =========================================================================
# 2. LOAD TPM MATRICES
# =========================================================================
cat("\n--- STEP 2: LOAD EXPRESSION DATA ---\n")

# Per-replicate expressed-transcript TPM files for this timepoint
# (order: Compatible replicates first, Incompatible replicates second —
# see 05_Transcript_Expression_Analysis.R for how these are produced)
samples <- c(paste0(timepoint, "_compatible1_expressed"),
             paste0(timepoint, "_compatible2_expressed"),
             paste0(timepoint, "_compatible4_expressed"),
             paste0(timepoint, "_incompatible1_expressed"),
             paste0(timepoint, "_incompatible2_expressed"),
             paste0(timepoint, "_incompatible3_expressed"))

paths <- file.path(base_dir, "Data", paste0("KallistoCounts_", timepoint),
                    "Expressed_01", paste0(samples, ".tsv"))

tpm_list <- lapply(seq_along(paths), function(i){
  df <- read.delim(paths[i], stringsAsFactors=FALSE)[, c("target_id","tpm")]
  colnames(df) <- c("target_id", samples[i])
  df
})

# Merge into one TPM matrix
tpm_mat <- Reduce(function(a,b) full_join(a,b, by="target_id"), tpm_list)

# Setup expression matrix (rows = transcripts, columns = samples)
expr <- as.matrix(tpm_mat[ , samples])
rownames(expr) <- tpm_mat$target_id
cat("Expression matrix built with", nrow(expr), "unique transcripts.\n")

# =========================================================================
# 3. CIS CORRELATION FUNCTION
# =========================================================================
# For each candidate lncRNA-mRNA cis pair:
#   - both transcripts must exceed `tpm_thresh` TPM in at least
#     `min_samples` of the 6 samples (checked independently and jointly)
#   - if that sparsity condition is met, Pearson and Spearman correlation
#     are computed on log2(TPM + 1)-transformed values
# Runs as a per-pair loop (appropriate here since the cis pair set is small
# relative to the genome-wide trans pair set in Part C).
run_robust_corr <- function(expr_mat, cis_df, min_samples, tpm_thresh, out_file) {
  
  cat(sprintf("\n--- STEP 3: CORRELATION & SPARSITY FILTER ---\n"))
  cat(sprintf("Rule: BOTH transcripts must have TPM > %.2f in >= %d samples\n\n", tpm_thresh, min_samples))
  
  # Tracking variables
  total_input <- nrow(cis_df)
  failed_filter <- 0
  tested_pairs <- 0
  
  # Pre-allocate a list to store results for speed
  res_list <- list()
  
  for (i in 1:nrow(cis_df)) {
    lnc <- cis_df$lnc_id[i]
    pc  <- cis_df$pc_id[i]
    
    # Ensure both exist in expression matrix
    if (lnc %in% rownames(expr_mat) && pc %in% rownames(expr_mat)) {
      
      # Extract RAW TPM for filtering
      lnc_vec <- as.numeric(expr_mat[lnc, ])
      pc_vec  <- as.numeric(expr_mat[pc, ])
      
      # 1. Apply Sparsity Filter on RAW TPM
      lnc_pass <- sum(lnc_vec > tpm_thresh)
      pc_pass  <- sum(pc_vec > tpm_thresh)
      both_pass <- sum(lnc_vec > tpm_thresh & pc_vec > tpm_thresh)
      
      # ONLY proceed if BOTH genes meet the minimum expression threshold
      if (lnc_pass >= min_samples && pc_pass >= min_samples) {
        tested_pairs <- tested_pairs + 1
        
        # Transform data to log2(TPM + 1) for correlation matching
        lnc_vec_log <- log2(lnc_vec + 1)
        pc_vec_log  <- log2(pc_vec + 1)
        
        # 2. Run Rank-Based Spearman Correlation (on log2 data)
        ct_spearman <- suppressWarnings(cor.test(lnc_vec_log, pc_vec_log, method = "spearman", exact = FALSE))
        
        # 3. Run Pearson Correlation (on log2 data)
        ct_pearson <- suppressWarnings(cor.test(lnc_vec_log, pc_vec_log, method = "pearson"))
        
        res_list[[length(res_list) + 1]] <- data.frame(
          lnc_id = lnc,
          pc_id = pc,
          pearson_r = unname(ct_pearson$estimate),
          pearson_pval = ct_pearson$p.value,
          spearman_rho = unname(ct_spearman$estimate),
          spearman_pval = ct_spearman$p.value,
          lnc_expr_samples = lnc_pass,
          pc_expr_samples = pc_pass,
          concurrent_samples = both_pass, 
          stringsAsFactors = FALSE
        )
      } else {
        failed_filter <- failed_filter + 1
      }
    } else {
      # Failed because one or both transcripts were entirely missing from the kallisto counts
      failed_filter <- failed_filter + 1 
    }
  }
  
  # Combine results and calculate FDR
  if (length(res_list) > 0) {
    results <- bind_rows(res_list)
    
    # FDR correction for BOTH methods
    results$pearson_padj <- p.adjust(results$pearson_pval, method = "BH")
    results$spearman_padj <- p.adjust(results$spearman_pval, method = "BH")
    
    # Rearrange columns slightly for clean output
    results <- results %>% 
      relocate(pearson_padj, .after = pearson_pval) %>%
      relocate(spearman_padj, .after = spearman_pval)
    
    # Calculate how many were significant in BOTH tests
    sig_both <- results %>% 
      filter(pearson_padj < 0.05 & spearman_padj < 0.05) %>% 
      nrow()
    
    cat("Pairs significant in BOTH tests       :", sig_both, "\n")
    
    write.table(results, file = out_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
    cat(sprintf("\nSaved %d tested pairs to: %s\n", nrow(results), basename(out_file)))
  } else {
    cat("\nNo pairs survived this strict filtering threshold!\n")
  }
}

# =========================================================================
# 4. EXECUTE FINAL PIPELINE
# =========================================================================

# Run the single, chosen scenario: TPM > 0.05 in at least 4/6 samples
run_robust_corr(expr, cis_filtered, min_samples = 4, tpm_thresh = 0.05, 
                out_file = file.path(out_dir, paste0(timepoint, "_cis_results.tsv")))

cat("\nPipeline complete.\n")
```

------------------------------------------------------------------------

## Part C — Trans Pair Construction & Correlation (R)

### Purpose

Builds the trans candidate pair set — every DE lncRNA paired against every DE protein-coding transcript genome-wide, **excluding** each lncRNA's own cis and overlap partners from Part A/B — then applies the same sparsity filter and correlation testing as Part B, using a vectorized pre-filter for speed given the much larger number of candidate pairs.

### Inputs

| File | Description |
|------------------------------------|------------------------------------|
| `10_cis.tsv` | Final cis-pair table from Part A (used here only to identify each lncRNA's cis partners to exclude) |
| `10_DE_LncRNAs_overlap.bed` | Overlap-pair BED file from Part A (used here only to identify each lncRNA's overlap partners to exclude) |
| `10_DE_mRNAs_transid.txt` | DE protein-coding transcript IDs |
| `10_compatible1_expressed.tsv`, `10_compatible2_expressed.tsv`, `10_compatible4_expressed.tsv`, `10_incompatible1_expressed.tsv`, `10_incompatible2_expressed.tsv`, `10_incompatible3_expressed.tsv` | Per-replicate expressed-transcript Kallisto TPM files |

### Outputs

| File | Description |
|------------------------------------|------------------------------------|
| `10_trans_results` | Correlation results (Pearson + Spearman, raw p-values) for trans pairs passing the sparsity filter |

### Script

``` r
# =============================================================================
# Trans Pair Construction & Expression Correlation (t=10 min)
# =============================================================================
# Identical for t=60: replace all "10_" prefixes below (file names, sample
# names, and directory names) with "60_".
#
# PATH NOTE: `base_dir` and `timepoint` are the only values that should need
# editing to point this at your local project. The original script also had
# an inconsistent path in one place ("/users/..." lowercase vs "/Users/..."
# elsewhere) — this has been unified via `base_dir`. All other logic and file
# names are unchanged from the original script.
# =============================================================================

# Load required packages
library(dplyr)
library(readr)
library(data.table)

base_dir  <- "<path_to_project>/Analysis"
timepoint <- "10"

cat("\n######################################################################\n")
cat(" TRANS-CORRELATION PIPELINE (t=", timepoint, " min) \n")
cat("######################################################################\n")

# =========================================================================
# 1. READ CIS, OVERLAP, AND TARGET LISTS
# =========================================================================
cat("\n--- STEP 1: LOAD ANNOTATIONS & BUILD TRANS PAIRS ---\n")

# cis pairs from Part A/B (used here only to exclude each lncRNA's cis partners)
cis <- fread(file.path(base_dir, "Results", paste0("t", timepoint),
                       paste0(timepoint, "_cis"), paste0(timepoint, "_cis.tsv")),
             header = FALSE)
colnames(cis)[c(4,9)] <- c("lnc_id","pc_id")

# overlap pairs from Part A (used here only to exclude each lncRNA's overlap partners)
overlap <- fread(file.path(base_dir, "Results", paste0("t", timepoint),
                            paste0(timepoint, "_cis"),
                            paste0(timepoint, "_DE_LncRNAs_overlap.bed")),
                  header = FALSE)
colnames(overlap)[c(4,16)] <- c("lnc_id","pc_id")

# DE protein-coding transcript IDs
de_pc_ids <- fread(file.path(base_dir, "Results", paste0("t", timepoint),
                              paste0(timepoint, "_DE_mRNAs_transid.txt")),
                    header = FALSE)[[1]]

# Combine cis + overlap partners per lncRNA, to EXCLUDE from the trans set
cis_overlap <- bind_rows(cis[,c("lnc_id","pc_id")],
                         overlap[,c("lnc_id","pc_id")]) %>% distinct()

# list of cis/overlap partners for each lncRNA
cis_list <- split(cis_overlap$pc_id, cis_overlap$lnc_id)
all_lnc_ids <- unique(cis_overlap$lnc_id)

# Build combinations of ALL lncRNAs against ALL DE mRNAs (excluding cis/overlap)
make_trans_pairs <- function(pc_ids) {
  trans_pairs_list <- vector("list", length(all_lnc_ids))
  names(trans_pairs_list) <- all_lnc_ids
  
  for (lnc in all_lnc_ids) {
    cis_pcs <- cis_list[[lnc]]
    trans_pcs <- setdiff(pc_ids, cis_pcs)
    trans_pcs <- unique(trans_pcs)
    
    # Remove any pc_id that is the same as the current lnc OR any lnc
    trans_pcs <- setdiff(trans_pcs, all_lnc_ids)
    
    trans_pairs_list[[lnc]] <- data.frame(
      lnc_id = lnc,
      pc_id  = trans_pcs,
      stringsAsFactors = FALSE
    )
  }
  return(do.call(rbind, trans_pairs_list))
}

lnc_trans_pairs_de <- make_trans_pairs(de_pc_ids)
cat("Total raw trans pairs generated (excluding cis) :", nrow(lnc_trans_pairs_de), "\n")


# =========================================================================
# 2. LOAD TPM MATRICES
# =========================================================================
cat("\n--- STEP 2: LOAD EXPRESSION DATA ---\n")

# Per-replicate expressed-transcript TPM files for this timepoint
samples <- c(paste0(timepoint, "_compatible1_expressed"),
             paste0(timepoint, "_compatible2_expressed"),
             paste0(timepoint, "_compatible4_expressed"),
             paste0(timepoint, "_incompatible1_expressed"),
             paste0(timepoint, "_incompatible2_expressed"),
             paste0(timepoint, "_incompatible3_expressed"))

paths <- file.path(base_dir, "Data", paste0("KallistoCounts_", timepoint),
                    "Expressed_01", paste0(samples, ".tsv"))

tpm_list <- lapply(seq_along(paths), function(i){
  df <- read.delim(paths[i], stringsAsFactors=FALSE)[, c("target_id","tpm")]
  colnames(df) <- c("target_id", samples[i])
  df
})

# Merge into one TPM matrix
tpm_mat <- Reduce(function(a,b) full_join(a,b, by="target_id"), tpm_list)

# Global zero filter: drop transcripts with zero TPM in every sample
all_tpms <- as.matrix(tpm_mat[ , samples])
keep_any  <- rowSums(all_tpms > 0) > 0
tpm_mat <- tpm_mat[keep_any, ]

# Setup expression matrix (rows = transcripts, columns = samples)
expr <- as.matrix(tpm_mat[ , samples])
rownames(expr) <- tpm_mat$target_id
cat("Expression matrix built with", nrow(expr), "unique transcripts.\n")


# =========================================================================
# 3. TRANS CORRELATION FUNCTION
# =========================================================================
# Same sparsity + correlation logic as the cis pipeline (Part B), but the
# sparsity pre-filter is vectorized up front (rather than checked per-pair
# inside the loop) since the trans pair set is far larger genome-wide.
run_robust_trans_corr <- function(expr_mat, trans_df, min_samples, tpm_thresh, out_dir) {
  
  cat(sprintf("\n--- STEP 3: CORRELATION & SPARSITY FILTER ---\n"))
  cat(sprintf("Rule: BOTH transcripts must have TPM > %.2f in >= %d samples\n\n", tpm_thresh, min_samples))
  
  # 1. Ultra-fast Pre-filter (Vectorized evaluated strictly on raw TPM)
  pass_counts <- rowSums(expr_mat > tpm_thresh)
  valid_transcripts <- names(pass_counts[pass_counts >= min_samples])
  
  initial_pairs <- nrow(trans_df)
  
  # Filter the combinations table to ONLY those where both genes pass
  trans_df <- trans_df %>%
    mutate(
      lnc_valid = lnc_id %in% valid_transcripts,
      pc_valid = pc_id %in% valid_transcripts
    )
  
  tested_df <- trans_df %>% filter(lnc_valid & pc_valid)
  failed_filter <- initial_pairs - nrow(tested_df)
  
  if(nrow(tested_df) == 0) {
    cat("No pairs survived the strict filtering threshold!\n")
    return(NULL)
  }
  
  # 2. Pre-apply log2(x + 1) transformation to ensure Pearson values match original script
  cat("Applying log2(TPM + 1) transformation to expression matrix for correlation matching...\n")
  expr_mat_log <- log2(expr_mat + 1)
  
  n_pairs <- nrow(tested_df)
  pearson_r <- numeric(n_pairs)
  pearson_p <- numeric(n_pairs)
  spearman_rho <- numeric(n_pairs)
  spearman_p <- numeric(n_pairs)
  
  cat("Starting correlations (This may take a moment for large trans networks)...\n")
  for (i in seq_len(n_pairs)) {
    # Extract directly from the log2 transformed matrix for speed & consistency
    x <- as.numeric(expr_mat_log[tested_df$lnc_id[i], ])
    y <- as.numeric(expr_mat_log[tested_df$pc_id[i], ])
    
    ct_p <- suppressWarnings(cor.test(x, y, method = "pearson"))
    ct_s <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
    
    pearson_r[i] <- ct_p$estimate
    pearson_p[i] <- ct_p$p.value
    spearman_rho[i] <- ct_s$estimate
    spearman_p[i] <- ct_s$p.value
    
    if (i %% 10000 == 0) {
      cat(sprintf("  Processed %d / %d pairs (%.1f%%)...\n", i, n_pairs, (i/n_pairs)*100))
    }
  }
  
  # 3. Compile Results
  results <- data.frame(
    lnc_id = tested_df$lnc_id,
    pc_id  = tested_df$pc_id,
    pearson_r = pearson_r,
    pearson_pval = pearson_p,
    spearman_rho = spearman_rho,
    spearman_pval = spearman_p,
    lnc_expr_samples = pass_counts[tested_df$lnc_id],
    pc_expr_samples  = pass_counts[tested_df$pc_id],
    stringsAsFactors = FALSE
  )
  
  # Count significance hits
  sig_both <- sum(results$pearson_pval < 0.05 & results$spearman_pval < 0.05, na.rm = TRUE)
  cat("\nPairs significant (raw p < 0.05) in BOTH tests :", sig_both, "\n")
  
  # Save the file
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  out_file <- file.path(out_dir, paste0(timepoint, "_trans_results.tsv"))
  
  write.table(results, file = out_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  cat(sprintf("Saved %d tested pairs to: %s\n", nrow(results), basename(out_file)))
  
  return(results)
}

# =========================================================================
# 4. EXECUTE FINAL PIPELINE
# =========================================================================

out_dir_new <- file.path(base_dir, "Results", paste0("t", timepoint),
                          paste0(timepoint, "_trans"), "New")

res_trans <- run_robust_trans_corr(
  expr_mat = expr, 
  trans_df = lnc_trans_pairs_de, 
  min_samples = 4, 
  tpm_thresh = 0.05, 
  out_dir = out_dir_new
)

cat("\nTrans Pipeline complete.\n")
```

------------------------------------------------------------------------
