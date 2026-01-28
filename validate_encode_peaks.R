#!/usr/bin/env Rscript
# ENCODE Peak Validation Script
# Tests wavematch accuracy against ENCODE gold-standard peaks
# Automatically detects files in the encode_h3k27ac_k562 directory

library(wavematch)
library(GenomicRanges)
library(IRanges)

cat("=", paste(rep("=", 70), collapse=""), "=\n")
cat("ENCODE Peak Validation - wavematch vs ENCODE Gold Standard\n")
cat("=", paste(rep("=", 70), collapse=""), "=\n\n")

# ============================================================================
# STEP 1: Auto-detect files
# ============================================================================

encode_dir <- "/Users/away/Downloads/encode_h3k27ac_k562"
if (!dir.exists(encode_dir)) {
  stop("ENCODE directory not found: ", encode_dir)
}

cat("[STEP 1] Detecting files...\n")

# Find narrowPeak file
narrowpeak_file <- file.path(encode_dir, "ENCFF000BWW.narrowPeak")
if (!file.exists(narrowpeak_file)) {
  stop("ENCODE narrowPeak file not found: ", narrowpeak_file)
}
cat("  ✓ Found ENCODE peaks:", narrowpeak_file, "\n")

# Find BAM files (use the main ChIP sample - ENCFF000BWX is typically the ChIP)
bam_files <- list.files(encode_dir, pattern = "\\.bam$", full.names = TRUE)
if (length(bam_files) == 0) {
  stop("No BAM files found in ", encode_dir)
}

# Prefer ENCFF000BWX.bam as it's typically the ChIP sample matching ENCFF000BWW peaks
chip_bam <- bam_files[grepl("ENCFF000BWX", bam_files)]
if (length(chip_bam) == 0) {
  chip_bam <- bam_files[1]  # Fallback to first BAM
  cat("  ⚠ Using first available BAM:", chip_bam, "\n")
} else {
  chip_bam <- chip_bam[1]
  cat("  ✓ Found ChIP BAM:", chip_bam, "\n")
}

# Test chromosome (chr22 is small and fast)
test_chromosome <- "chr22"
genome <- "hg38"
binSize <- 25L

cat("  Test chromosome:", test_chromosome, "\n")
cat("  Genome:", genome, "\n")
cat("  Bin size:", binSize, "bp\n\n")

# ============================================================================
# STEP 2: Load ENCODE peaks
# ============================================================================

cat("[STEP 2] Loading ENCODE peaks...\n")

encode_peaks <- utils::read.table(
  narrowpeak_file,
  sep = "\t",
  header = FALSE,
  stringsAsFactors = FALSE,
  comment.char = "",
  quote = ""
)

# narrowPeak format: chrom, chromStart, chromEnd, name, score, strand, 
#                    signalValue, pValue, qValue, peak
if (ncol(encode_peaks) < 3) {
  stop("ENCODE peaks file must have at least 3 columns (chrom, start, end)")
}

colnames(encode_peaks)[1:3] <- c("chromosome", "start", "end")
encode_peaks$chromosome <- as.character(encode_peaks$chromosome)
encode_peaks$start <- as.integer(encode_peaks$start)
encode_peaks$end <- as.integer(encode_peaks$end)

cat("  ✓ Loaded ", nrow(encode_peaks), " total peaks\n")

# Filter to test chromosome
encode_peaks_chr <- encode_peaks[encode_peaks$chromosome == test_chromosome, ]
cat("  ✓ Filtered to ", nrow(encode_peaks_chr), " peaks on ", test_chromosome, "\n\n")

if (nrow(encode_peaks_chr) == 0) {
  stop("No ENCODE peaks found on ", test_chromosome, 
       "\nAvailable chromosomes: ", paste(unique(encode_peaks$chromosome)[1:10], collapse=", "))
}

# Convert to GRanges (BED is 0-based, GRanges is 1-based)
encode_gr <- GRanges(
  seqnames = encode_peaks_chr$chromosome,
  ranges = IRanges(start = encode_peaks_chr$start + 1L, end = encode_peaks_chr$end)
)

# ============================================================================
# STEP 3: Load track data and run wavematch
# ============================================================================

cat("[STEP 3] Loading track data from BAM...\n")
load_start <- Sys.time()

track_data <- tryCatch({
  wavematchFromBam(
    bamFile = chip_bam,
    genome = genome,
    chrom = test_chromosome,
    binSize = binSize
  )
}, error = function(e) {
  cat("  ✗ Error loading BAM: ", e$message, "\n")
  cat("  Trying alternative chromosome...\n")
  # Try chr1 as fallback
  track_data <- wavematchFromBam(
    bamFile = chip_bam,
    genome = genome,
    chrom = "chr1",
    binSize = binSize
  )
  test_chromosome <<- "chr1"
  return(track_data)
})

load_time <- as.numeric(difftime(Sys.time(), load_start, units = "secs"))
cat("  ✓ Loaded ", length(track_data$values), " bins in ", round(load_time, 2), " seconds\n")
cat("  ✓ Chromosome:", track_data$chromosome, "\n\n")

# Update test chromosome if it changed
if (track_data$chromosome != test_chromosome) {
  test_chromosome <- track_data$chromosome
  # Re-filter ENCODE peaks to match
  encode_peaks_chr <- encode_peaks[encode_peaks$chromosome == test_chromosome, ]
  if (nrow(encode_peaks_chr) == 0) {
    stop("No ENCODE peaks found on ", test_chromosome)
  }
  encode_gr <- GRanges(
    seqnames = encode_peaks_chr$chromosome,
    ranges = IRanges(start = encode_peaks_chr$start + 1L, end = encode_peaks_chr$end)
  )
  cat("  ⚠ Switched to chromosome:", test_chromosome, "\n")
  cat("  ✓ ENCODE peaks on ", test_chromosome, ": ", length(encode_gr), "\n\n")
}

cat("[STEP 4] Running wavematch (optimized algorithm)...\n")
match_start <- Sys.time()

detected_peaks <- wavematchFromTrack(
  chromosome = track_data$chromosome,
  intervals = track_data$intervals,
  values = track_data$values,
  template_names = "haar",
  cascade_levels = 3L,
  iters = 20000L,
  alpha = 0.05,
  use_split_halves = TRUE,
  recenter_at_point_source = TRUE,
  mergePeaks = TRUE,
  mergeGapBP = 150L,
  verbose = TRUE
)

match_time <- as.numeric(difftime(Sys.time(), match_start, units = "secs"))
cat("\n  ✓ Peak calling completed in ", round(match_time, 2), " seconds\n")
cat("  ✓ Detected ", nrow(detected_peaks), " peaks\n\n")

# ============================================================================
# STEP 5: Compute overlap metrics
# ============================================================================

cat("[STEP 5] Computing overlap metrics...\n")

# Convert detected peaks to GRanges
if (nrow(detected_peaks) > 0) {
  detected_gr <- asGRanges(detected_peaks, keep.mcols = TRUE)
} else {
  detected_gr <- GRanges()
}

# Find overlaps with minimum overlap fraction (10% required for match)
min_overlap_frac <- 0.1
overlaps <- findOverlaps(detected_gr, encode_gr, minoverlap = min_overlap_frac)

matched_detected <- unique(queryHits(overlaps))
matched_encode <- unique(subjectHits(overlaps))

TP <- length(matched_detected)  # True Positives: detected peaks that match ENCODE
FP <- length(detected_gr) - TP   # False Positives: detected peaks with no ENCODE match
FN <- length(encode_gr) - length(matched_encode)  # False Negatives: ENCODE peaks not detected

# Calculate metrics
precision <- if (TP + FP > 0) TP / (TP + FP) else 0.0
recall <- if (TP + FN > 0) TP / (TP + FN) else 0.0
F1 <- if (precision + recall > 0) 2 * precision * recall / (precision + recall) else 0.0

# Jaccard index: intersection / union
if (TP > 0) {
  detected_matched <- detected_gr[matched_detected]
  encode_matched <- encode_gr[matched_encode]
  
  intersection <- sum(width(intersect(detected_matched, encode_matched)))
  union <- sum(width(union(detected_matched, encode_matched)))
  jaccard <- if (union > 0) intersection / union else 0.0
} else {
  jaccard <- 0.0
}

cat("  ✓ Computed overlap metrics\n\n")

# ============================================================================
# STEP 6: Print results
# ============================================================================

cat(rep("=", 72), "\n", sep="")
cat("VALIDATION RESULTS\n")
cat(rep("=", 72), "\n\n", sep="")

cat("Configuration:\n")
cat("  BAM file:           ", basename(chip_bam), "\n")
cat("  ENCODE peaks:       ", basename(narrowpeak_file), "\n")
cat("  Test chromosome:    ", test_chromosome, "\n")
cat("  Genome:             ", genome, "\n")
cat("  Bin size:           ", binSize, " bp\n\n")

cat("Peak Counts:\n")
cat("  ENCODE peaks (gold standard): ", length(encode_gr), "\n")
cat("  Detected peaks:              ", length(detected_gr), "\n")
cat("  True Positives (TP):         ", TP, "\n")
cat("  False Positives (FP):        ", FP, "\n")
cat("  False Negatives (FN):        ", FN, "\n\n")

cat("Metrics:\n")
cat("  Precision:  ", sprintf("%.4f", precision), " (", round(precision * 100, 2), "%)\n")
cat("  Recall:     ", sprintf("%.4f", recall), " (", round(recall * 100, 2), "%)\n")
cat("  F1-score:   ", sprintf("%.4f", F1), "\n")
cat("  Jaccard:    ", sprintf("%.4f", jaccard), "\n\n")

cat("Performance:\n")
cat("  Load time:      ", round(load_time, 2), " seconds\n")
cat("  Peak calling:   ", round(match_time, 2), " seconds\n")
cat("  Total time:     ", round(load_time + match_time, 2), " seconds\n\n")

# Peak width statistics
if (nrow(detected_peaks) > 0) {
  detected_widths <- detected_peaks$end - detected_peaks$start
  encode_widths <- encode_peaks_chr$end - encode_peaks_chr$start
  
  cat("Peak Width Statistics:\n")
  cat("  Detected peaks: ", min(detected_widths), "-", max(detected_widths), 
      " bp (mean: ", round(mean(detected_widths), 1), ", median: ", 
      round(median(detected_widths), 1), ")\n")
  cat("  ENCODE peaks:   ", min(encode_widths), "-", max(encode_widths), 
      " bp (mean: ", round(mean(encode_widths), 1), ", median: ", 
      round(median(encode_widths), 1), ")\n\n")
}

# Success criteria
cat("Success Criteria:\n")
cat("  Precision ≥ 0.90: ", if (precision >= 0.90) "✓ PASS" else "✗ FAIL", "\n")
cat("  Recall ≥ 0.90:    ", if (recall >= 0.90) "✓ PASS" else "✗ FAIL", "\n")
cat("  F1 ≥ 0.90:        ", if (F1 >= 0.90) "✓ PASS" else "✗ FAIL", "\n")
cat("  Jaccard ≥ 0.70:   ", if (jaccard >= 0.70) "✓ PASS" else "✗ FAIL", "\n\n")

cat(rep("=", 72), "\n", sep="")
cat("Validation complete!\n")
cat(rep("=", 72), "\n", sep="")

# Save results to CSV
results_df <- data.frame(
  method = "wavematch",
  chromosome = test_chromosome,
  n_encode_peaks = length(encode_gr),
  n_detected_peaks = length(detected_gr),
  TP = TP,
  FP = FP,
  FN = FN,
  precision = precision,
  recall = recall,
  F1 = F1,
  jaccard = jaccard,
  load_time_sec = load_time,
  peak_calling_time_sec = match_time,
  total_time_sec = load_time + match_time,
  bam_file = basename(chip_bam),
  encode_peaks_file = basename(narrowpeak_file),
  stringsAsFactors = FALSE
)

results_file <- file.path(encode_dir, paste0("wavematch_validation_", test_chromosome, ".csv"))
write.csv(results_df, results_file, row.names = FALSE)
cat("\nResults saved to:", results_file, "\n")
