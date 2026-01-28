#' Wavelet-based peak calling from multiple samples with shared uncertainty weights
#'
#' Process multiple samples individually, but use a shared uncertainty
#' bedGraph to improve peak calling for all samples.
#'
#' @param sampleData data.frame with columns:
#'   \describe{
#'     \item{sample}{character, unique sample identifiers}
#'     \item{bamFile}{character, paths to BAM files (if using BAM input)}
#'     \item{bedGraphFile}{character, paths to bedGraph files (if using bedGraph input)}
#'   }
#'   Optional columns:
#'   \describe{
#'     \item{batch}{character/factor, batch information}
#'     \item{condition}{character/factor, experimental conditions}
#'     \item{...}{any other metadata columns}
#'   }
#' @param consenrichUncertaintyBedGraph path to uncertainty bedGraph file.
#'   This will be used as weights for ALL samples to improve peak calling.
#' @param genome genome name or assembly object (required if using BAM files)
#' @param inputType character, either "bam" or "bedgraph" (default: auto-detect from sampleData)
#' @param binSize integer, bin size in bp (only used if inputType="bam", default: 25)
#' @param useUncertaintyWeights logical, use uncertainty values as weights (default: TRUE)
#' @param uncertaintyToWeights function to convert uncertainty to weights.
#'   Default: \code{function(u) 1 / (u + 0.01)}
#' @param validateFiles logical, check all files exist before processing (default: TRUE)
#' @param combineResults logical, combine all samples into single data.frame (default: FALSE)
#' @param returnGRanges logical, return GRanges objects (default: FALSE)
#' @param ncores integer, number of cores for parallel processing
#' @param verbose logical, print progress messages
#' @param ... arguments passed to wavematchFromTrack
#'
#' @return If combineResults=FALSE: named list of peak data.frames/GRanges (one per sample)
#'         If combineResults=TRUE: single data.frame/GRanges with sample column
#'
#' @export
#' @importFrom parallel mclapply detectCores
wavematchFromSamples <- function(
  sampleData,
  consenrichUncertaintyBedGraph = NULL,
  genome = NULL,
  inputType = c("auto", "bam", "bedgraph"),
  binSize = 25L,
  useUncertaintyWeights = TRUE,
  uncertaintyToWeights = function(u) 1 / (u + 0.01),
  validateFiles = TRUE,
  combineResults = FALSE,
  returnGRanges = FALSE,
  ncores = 1L,
  verbose = TRUE,
  ...
) {
  
  # Validate sampleData
  if (!is.data.frame(sampleData) || nrow(sampleData) == 0L) {
    stop("sampleData must be a non-empty data.frame")
  }
  
  required <- c("sample")
  missing <- setdiff(required, colnames(sampleData))
  if (length(missing)) {
    stop("sampleData missing required column: ", paste(missing, collapse = ", "))
  }
  
  # Auto-detect input type
  inputType <- match.arg(inputType)
  if (inputType == "auto") {
    if ("bamFile" %in% colnames(sampleData)) {
      inputType <- "bam"
    } else if ("bedGraphFile" %in% colnames(sampleData)) {
      inputType <- "bedgraph"
    } else {
      stop("sampleData must contain either 'bamFile' or 'bedGraphFile' column")
    }
  }
  
  # Validate required columns based on input type
  if (inputType == "bam") {
    if (!"bamFile" %in% colnames(sampleData)) {
      stop("sampleData must contain 'bamFile' column when inputType='bam'")
    }
    if (is.null(genome)) {
      stop("genome must be provided when using BAM files")
    }
  } else {
    if (!"bedGraphFile" %in% colnames(sampleData)) {
      stop("sampleData must contain 'bedGraphFile' column when inputType='bedgraph'")
    }
  }
  
  # Validate files exist
  if (validateFiles) {
    # Check sample files
    sample_files <- if (inputType == "bam") {
      sampleData$bamFile
    } else {
      sampleData$bedGraphFile
    }
    
    missing_files <- sample_files[!file.exists(sample_files)]
    if (length(missing_files)) {
      stop("Sample files not found:\n", paste(missing_files, collapse = "\n"))
    }
    
    # Check uncertainty file
    if (!is.null(consenrichUncertaintyBedGraph)) {
      if (!file.exists(consenrichUncertaintyBedGraph)) {
        stop("Uncertainty bedGraph not found: ", consenrichUncertaintyBedGraph)
      }
    }
    
    # Check for duplicate samples
    if (any(duplicated(sampleData$sample))) {
      stop("Duplicate sample names found in sampleData$sample")
    }
  }
  
  # Load shared uncertainty weights (once for all samples)
  shared_weights_bg <- NULL
  if (useUncertaintyWeights && !is.null(consenrichUncertaintyBedGraph)) {
    if (verbose) {
      message("Loading shared uncertainty weights from: ", 
              consenrichUncertaintyBedGraph)
    }
    
    unc_bg <- utils::read.table(
      consenrichUncertaintyBedGraph,
      sep = "\t",
      header = FALSE,
      stringsAsFactors = FALSE
    )
    
    if (ncol(unc_bg) < 4L) {
      stop("Uncertainty bedGraph must have at least 4 columns: chrom, start, end, value")
    }
    
    colnames(unc_bg)[1:4] <- c("chromosome", "start", "end", "uncertainty")
    unc_bg$chromosome <- as.character(unc_bg$chromosome)
    unc_bg$start <- as.numeric(unc_bg$start)
    unc_bg$end <- as.numeric(unc_bg$end)
    unc_bg$uncertainty <- as.numeric(unc_bg$uncertainty)
    
    # Convert uncertainty to weights
    unc_bg$weight <- uncertaintyToWeights(unc_bg$uncertainty)
    shared_weights_bg <- unc_bg
    
    if (verbose) {
      message("  Converted uncertainty to weights (range: ", 
              round(min(unc_bg$weight, na.rm = TRUE), 3), " - ",
              round(max(unc_bg$weight, na.rm = TRUE), 3), ")")
    }
  }
  
  # Process each sample
  process_sample <- function(i) {
    row <- sampleData[i, ]
    sample_name <- as.character(row$sample)
    
    if (verbose) {
      message("Processing sample ", i, "/", nrow(sampleData), ": ", sample_name)
    }
    
    # Get signal for this sample
    if (inputType == "bam") {
      # Process BAM file genome-wide
      peaks <- wavematchFromBamGenome(
        bamFile = as.character(row$bamFile),
        genome = genome,
        binSize = binSize,
        weightsBedGraph = consenrichUncertaintyBedGraph,
        returnGRanges = returnGRanges,
        verbose = verbose && nrow(sampleData) == 1L,  # Only verbose for single sample
        ...
      )
    } else {
      # Process bedGraph file
      peaks <- wavematchFromBedGraph(
        bedGraphFile = as.character(row$bedGraphFile),
        weightsBedGraph = consenrichUncertaintyBedGraph,
        returnGRanges = returnGRanges,
        verbose = verbose && nrow(sampleData) == 1L,
        ...
      )
    }
    
    # Add all metadata columns
    if (is.data.frame(peaks)) {
      meta_cols <- setdiff(colnames(sampleData), 
                          c("bamFile", "bedGraphFile"))
      for (col in meta_cols) {
        peaks[[col]] <- row[[col]]
      }
    } else if (inherits(peaks, "GRanges")) {
      meta_cols <- setdiff(colnames(sampleData),
                          c("bamFile", "bedGraphFile"))
      for (col in meta_cols) {
        GenomicRanges::mcols(peaks)[[col]] <- row[[col]]
      }
    }
    
    peaks
  }
  
  # Process samples (parallel or sequential)
  if (ncores > 1L && nrow(sampleData) > 1L) {
    max_cores <- parallel::detectCores()
    ncores <- min(ncores, max_cores, nrow(sampleData))
    
    if (verbose) {
      message("Using ", ncores, " cores for parallel processing")
    }
    
    results <- parallel::mclapply(
      seq_len(nrow(sampleData)),
      process_sample,
      mc.cores = ncores
    )
  } else {
    results <- lapply(seq_len(nrow(sampleData)), process_sample)
  }
  
  # Remove NULL results
  results <- results[!vapply(results, is.null, logical(1L))]
  
  if (length(results) == 0L) {
    if (returnGRanges) {
      return(GenomicRanges::GRanges())
    }
    return(data.frame(
      chromosome = character(0),
      start = numeric(0),
      end = numeric(0),
      stringsAsFactors = FALSE
    ))
  }
  
  # Name results
  names(results) <- sampleData$sample[seq_along(results)]
  
  # Combine if requested
  if (combineResults) {
    if (all(vapply(results, inherits, logical(1L), "GRanges"))) {
      result <- do.call(c, results)
    } else {
      result <- do.call(rbind, results)
      rownames(result) <- NULL
    }
    
    if (returnGRanges && !inherits(result, "GRanges")) {
      result <- asGRanges(result, keep.mcols = TRUE)
    }
    
    return(result)
  }
  
  # Convert to GRanges if requested
  if (returnGRanges) {
    results <- lapply(results, function(x) {
      if (!inherits(x, "GRanges")) {
        asGRanges(x, keep.mcols = TRUE)
      } else {
        x
      }
    })
  }
  
  results
}
