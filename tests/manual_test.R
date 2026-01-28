# Manual test script for wavematch package
# Run this after loading the package: devtools::load_all()

cat("=== wavematch Manual Test Script ===\n\n")

# Create test directory
test_dir <- file.path(getwd(), "test_output")
if (!dir.exists(test_dir)) {
  dir.create(test_dir, recursive = TRUE)
}

# Generate pseudo bedGraph data
generate_pseudo_bedgraph <- function(chromosome = "chr1", 
                                     n_bins = 1000, 
                                     bin_size = 25,
                                     signal_strength = 10,
                                     noise_level = 1) {
  starts <- seq(0, by = bin_size, length.out = n_bins)
  ends <- starts + bin_size
  
  # Create some peaks
  peak_positions <- c(100, 300, 500, 700, 900)
  peak_width <- 5  # bins
  
  values <- rnorm(n_bins, mean = noise_level, sd = noise_level)
  
  # Add peaks
  for (pos in peak_positions) {
    peak_start <- max(1, pos - peak_width)
    peak_end <- min(n_bins, pos + peak_width)
    values[peak_start:peak_end] <- values[peak_start:peak_end] + signal_strength
  }
  
  # Ensure non-negative
  values <- pmax(0, values)
  
  data.frame(
    chromosome = rep(chromosome, n_bins),
    start = starts,
    end = ends,
    value = values,
    stringsAsFactors = FALSE
  )
}

# Generate pseudo uncertainty bedGraph (Consenrich-style)
generate_pseudo_uncertainty <- function(bedgraph_df, 
                                        base_uncertainty = 0.1,
                                        peak_uncertainty = 0.5) {
  uncertainty <- rep(base_uncertainty, nrow(bedgraph_df))
  
  # Higher uncertainty at peaks (lower confidence)
  high_signal <- bedgraph_df$value > quantile(bedgraph_df$value, 0.9)
  uncertainty[high_signal] <- peak_uncertainty
  
  data.frame(
    chromosome = bedgraph_df$chromosome,
    start = bedgraph_df$start,
    end = bedgraph_df$end,
    uncertainty = uncertainty,
    stringsAsFactors = FALSE
  )
}

cat("1. Generating pseudo bedGraph files...\n")
set.seed(123)

# Create signal bedGraph
bg <- generate_pseudo_bedgraph("chr1", n_bins = 1000)
bg_file <- file.path(test_dir, "test_signal.bedGraph")
write.table(bg[, c("chromosome", "start", "end", "value")],
            file = bg_file,
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)
cat("   Created:", bg_file, "\n")

# Create uncertainty bedGraph
unc_bg <- generate_pseudo_uncertainty(bg)
unc_file <- file.path(test_dir, "test_uncertainty.bedGraph")
write.table(unc_bg[, c("chromosome", "start", "end", "uncertainty")],
            file = unc_file,
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)
cat("   Created:", unc_file, "\n")

# Create multiple sample files
cat("\n2. Creating multiple sample files...\n")
for (i in 1:3) {
  bg_i <- generate_pseudo_bedgraph("chr1", n_bins = 1000, signal_strength = 10 + i)
  bg_i_file <- file.path(test_dir, paste0("sample", i, ".bedGraph"))
  write.table(bg_i[, c("chromosome", "start", "end", "value")],
              file = bg_i_file,
              sep = "\t",
              quote = FALSE,
              row.names = FALSE,
              col.names = FALSE)
  cat("   Created:", bg_i_file, "\n")
}

cat("\n3. Testing wavematchFromTrack...\n")
tryCatch({
  peaks <- wavematchFromTrack(
    chromosome = bg$chromosome[1],
    intervals = bg$start,
    values = bg$value,
    iters = 2000L,
    alpha = 0.05,
    verbose = FALSE
  )
  cat("   ✓ Found", nrow(peaks), "peaks\n")
}, error = function(e) {
  cat("   ✗ Error:", e$message, "\n")
})

cat("\n4. Testing wavematchFromBedGraph...\n")
tryCatch({
  peaks <- wavematchFromBedGraph(
    bedGraphFile = bg_file,
    iters = 2000L,
    alpha = 0.05,
    verbose = FALSE
  )
  cat("   ✓ Found", nrow(peaks), "peaks\n")
}, error = function(e) {
  cat("   ✗ Error:", e$message, "\n")
})

cat("\n5. Testing wavematchFromBedGraph with weights...\n")
tryCatch({
  peaks <- wavematchFromBedGraph(
    bedGraphFile = bg_file,
    weightsBedGraph = unc_file,
    iters = 2000L,
    alpha = 0.05,
    verbose = FALSE
  )
  cat("   ✓ Found", nrow(peaks), "peaks (with uncertainty weighting)\n")
}, error = function(e) {
  cat("   ✗ Error:", e$message, "\n")
})

cat("\n6. Testing wavematchFromSamples...\n")
tryCatch({
  sample_data <- data.frame(
    sample = c("sample1", "sample2", "sample3"),
    bedGraphFile = file.path(test_dir, paste0("sample", 1:3, ".bedGraph")),
    condition = c("treatment", "treatment", "control"),
    stringsAsFactors = FALSE
  )
  
  peaks_list <- wavematchFromSamples(
    sampleData = sample_data,
    consenrichUncertaintyBedGraph = unc_file,
    inputType = "bedgraph",
    iters = 2000L,
    alpha = 0.05,
    verbose = TRUE,
    ncores = 1L
  )
  
  cat("   ✓ Processed", length(peaks_list), "samples\n")
  for (name in names(peaks_list)) {
    cat("     -", name, ":", nrow(peaks_list[[name]]), "peaks\n")
  }
}, error = function(e) {
  cat("   ✗ Error:", e$message, "\n")
})

cat("\n7. Testing asGRanges...\n")
tryCatch({
  peaks <- wavematchFromTrack(
    chromosome = bg$chromosome[1],
    intervals = bg$start,
    values = bg$value,
    iters = 2000L,
    alpha = 0.05,
    verbose = FALSE
  )
  
  peaks_gr <- asGRanges(peaks, keep.mcols = TRUE)
  cat("   ✓ Converted to GRanges:", length(peaks_gr), "peaks\n")
}, error = function(e) {
  cat("   ✗ Error:", e$message, "\n")
})

cat("\n=== Test Summary ===\n")
cat("Test files created in:", test_dir, "\n")
cat("You can manually inspect the files and run individual tests.\n")
