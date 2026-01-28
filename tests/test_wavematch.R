# Test script for wavematch package
# Creates pseudo bedGraph files and tests functionality

library(testthat)

# Create temporary directory for test files
test_dir <- tempdir()
setwd(test_dir)

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

# Test 1: Basic wavematchFromTrack
test_that("wavematchFromTrack works with basic input", {
  # Create pseudo data
  bg <- generate_pseudo_bedgraph("chr1", n_bins = 500)
  
  # Test peak calling
  peaks <- wavematchFromTrack(
    chromosome = bg$chromosome[1],
    intervals = bg$start,
    values = bg$value,
    iters = 1000L,  # Reduced for testing
    alpha = 0.05
  )
  
  expect_true(is.data.frame(peaks))
  expect_true(nrow(peaks) > 0)
  expect_true(all(c("chromosome", "start", "end", "score", "p_value", "q_value") %in% colnames(peaks)))
})

# Test 2: wavematchFromBedGraph
test_that("wavematchFromBedGraph works", {
  # Create pseudo bedGraph file
  bg <- generate_pseudo_bedgraph("chr1", n_bins = 500)
  bg_file <- file.path(test_dir, "test_signal.bedGraph")
  write.table(bg[, c("chromosome", "start", "end", "value")],
              file = bg_file,
              sep = "\t",
              quote = FALSE,
              row.names = FALSE,
              col.names = FALSE)
  
  # Test peak calling
  peaks <- wavematchFromBedGraph(
    bedGraphFile = bg_file,
    iters = 1000L,
    alpha = 0.05,
    verbose = FALSE
  )
  
  expect_true(is.data.frame(peaks))
  expect_true(nrow(peaks) > 0)
})

# Test 3: wavematchFromBedGraph with weights
test_that("wavematchFromBedGraph works with uncertainty weights", {
  # Create pseudo bedGraph files
  bg <- generate_pseudo_bedgraph("chr1", n_bins = 500)
  bg_file <- file.path(test_dir, "test_signal2.bedGraph")
  write.table(bg[, c("chromosome", "start", "end", "value")],
              file = bg_file,
              sep = "\t",
              quote = FALSE,
              row.names = FALSE,
              col.names = FALSE)
  
  # Create uncertainty bedGraph
  unc_bg <- generate_pseudo_uncertainty(bg)
  unc_file <- file.path(test_dir, "test_uncertainty.bedGraph")
  write.table(unc_bg[, c("chromosome", "start", "end", "uncertainty")],
              file = unc_file,
              sep = "\t",
              quote = FALSE,
              row.names = FALSE,
              col.names = FALSE)
  
  # Test peak calling with weights
  peaks <- wavematchFromBedGraph(
    bedGraphFile = bg_file,
    weightsBedGraph = unc_file,
    iters = 1000L,
    alpha = 0.05,
    verbose = FALSE
  )
  
  expect_true(is.data.frame(peaks))
  expect_true(nrow(peaks) > 0)
})

# Test 4: wavematchFromSamples with bedGraph input
test_that("wavematchFromSamples works with bedGraph input", {
  # Create sample data
  bg1 <- generate_pseudo_bedgraph("chr1", n_bins = 500)
  bg2 <- generate_pseudo_bedgraph("chr1", n_bins = 500)
  
  bg1_file <- file.path(test_dir, "sample1.bedGraph")
  bg2_file <- file.path(test_dir, "sample2.bedGraph")
  
  write.table(bg1[, c("chromosome", "start", "end", "value")],
              file = bg1_file,
              sep = "\t",
              quote = FALSE,
              row.names = FALSE,
              col.names = FALSE)
  
  write.table(bg2[, c("chromosome", "start", "end", "value")],
              file = bg2_file,
              sep = "\t",
              quote = FALSE,
              row.names = FALSE,
              col.names = FALSE)
  
  # Create uncertainty bedGraph
  unc_bg <- generate_pseudo_uncertainty(bg1)
  unc_file <- file.path(test_dir, "shared_uncertainty.bedGraph")
  write.table(unc_bg[, c("chromosome", "start", "end", "uncertainty")],
              file = unc_file,
              sep = "\t",
              quote = FALSE,
              row.names = FALSE,
              col.names = FALSE)
  
  # Create sample data frame
  sample_data <- data.frame(
    sample = c("sample1", "sample2"),
    bedGraphFile = c(bg1_file, bg2_file),
    condition = c("treatment", "control"),
    stringsAsFactors = FALSE
  )
  
  # Test multi-sample processing
  peaks_list <- wavematchFromSamples(
    sampleData = sample_data,
    consenrichUncertaintyBedGraph = unc_file,
    inputType = "bedgraph",
    iters = 1000L,
    alpha = 0.05,
    verbose = FALSE,
    ncores = 1L
  )
  
  expect_true(is.list(peaks_list))
  expect_true(length(peaks_list) == 2)
  expect_true(all(c("sample1", "sample2") %in% names(peaks_list)))
  expect_true(is.data.frame(peaks_list$sample1))
  expect_true(is.data.frame(peaks_list$sample2))
})

# Test 5: wavematchFromSamples with combineResults
test_that("wavematchFromSamples combines results correctly", {
  # Create sample data (reuse from previous test)
  bg1 <- generate_pseudo_bedgraph("chr1", n_bins = 500)
  bg2 <- generate_pseudo_bedgraph("chr1", n_bins = 500)
  
  bg1_file <- file.path(test_dir, "sample1_combined.bedGraph")
  bg2_file <- file.path(test_dir, "sample2_combined.bedGraph")
  
  write.table(bg1[, c("chromosome", "start", "end", "value")],
              file = bg1_file,
              sep = "\t",
              quote = FALSE,
              row.names = FALSE,
              col.names = FALSE)
  
  write.table(bg2[, c("chromosome", "start", "end", "value")],
              file = bg2_file,
              sep = "\t",
              quote = FALSE,
              row.names = FALSE,
              col.names = FALSE)
  
  sample_data <- data.frame(
    sample = c("sample1", "sample2"),
    bedGraphFile = c(bg1_file, bg2_file),
    stringsAsFactors = FALSE
  )
  
  # Test with combineResults = TRUE
  peaks_combined <- wavematchFromSamples(
    sampleData = sample_data,
    inputType = "bedgraph",
    combineResults = TRUE,
    iters = 1000L,
    alpha = 0.05,
    verbose = FALSE,
    ncores = 1L
  )
  
  expect_true(is.data.frame(peaks_combined))
  expect_true("sample" %in% colnames(peaks_combined))
  expect_true(all(c("sample1", "sample2") %in% peaks_combined$sample))
})

# Test 6: asGRanges conversion
test_that("asGRanges converts peaks correctly", {
  bg <- generate_pseudo_bedgraph("chr1", n_bins = 500)
  
  peaks <- wavematchFromTrack(
    chromosome = bg$chromosome[1],
    intervals = bg$start,
    values = bg$value,
    iters = 1000L,
    alpha = 0.05
  )
  
  peaks_gr <- asGRanges(peaks, keep.mcols = TRUE)
  
  expect_true(inherits(peaks_gr, "GRanges"))
  expect_true(length(peaks_gr) == nrow(peaks))
})

cat("\nAll tests completed!\n")
