#' Wavelet-based peak calling from a BAM file (genome-wide)
#'
#' High-level function that processes all chromosomes in a BAM file and returns
#' peaks as a data frame or GenomicRanges object.
#'
#' @param bamFile path to BAM file
#' @param genome either a genome name such as "hg38" or an object with chromSizes
#' @param binSize bin size in base pairs (default: 25)
#' @param chromosomes optional character vector of chromosomes to process.
#'   If NULL, processes all chromosomes in the genome.
#' @param blacklist optional blacklist file path or genome name to filter peaks.
#'   If NULL, no blacklist filtering is applied.
#' @param weightsBedGraph optional path to uncertainty bedGraph file.
#'   If provided, uncertainty values will be converted to weights and used to
#'   improve peak calling in low-confidence regions.
#' @param returnGRanges logical, if TRUE return a GRanges object instead of data.frame
#' @param ncores integer, number of cores for parallel processing (default: 1)
#' @param verbose logical, print progress messages
#' @param ... additional arguments passed to \code{wavematchFromTrack}
#'
#' @return data.frame or GRanges object with detected peaks
#' @export
#' @importFrom parallel mclapply detectCores
#' @importFrom GenomicRanges GRanges
wavematchFromBamGenome <- function(bamFile,
                                   genome,
                                   binSize = 25L,
                                   chromosomes = NULL,
                                   blacklist = NULL,
                                   weightsBedGraph = NULL,
                                   returnGRanges = FALSE,
                                   ncores = 1L,
                                   verbose = TRUE,
                                   ...) {
  
  if (!file.exists(bamFile)) {
    stop("BAM file not found: ", bamFile)
  }
  
  # Load genome resources
  if (is.character(genome) || is.null(genome$chromSizes)) {
    assembly <- load_genome_resources(genome)
    genome_name <- assembly$name
  } else {
    assembly <- genome
    genome_name <- NULL
  }
  
  if (is.null(assembly$chromSizes) || is.null(names(assembly$chromSizes))) {
    stop("assembly$chromSizes must be a named numeric vector")
  }
  
  # Determine chromosomes to process
  if (is.null(chromosomes)) {
    chromosomes <- names(assembly$chromSizes)
  } else {
    chromosomes <- intersect(chromosomes, names(assembly$chromSizes))
  }
  
  if (length(chromosomes) == 0L) {
    stop("No valid chromosomes found")
  }
  
  # Load weights bedGraph if provided
  weights_bg <- NULL
  if (!is.null(weightsBedGraph) && file.exists(weightsBedGraph)) {
    if (verbose) {
      message("Reading weights from: ", weightsBedGraph)
    }
    weights_bg <- utils::read.table(
      weightsBedGraph,
      sep = "\t",
      header = FALSE,
      stringsAsFactors = FALSE
    )
    if (ncol(weights_bg) < 4L) {
      stop("Weights bedGraph must have at least 4 columns: chrom, start, end, value")
    }
    colnames(weights_bg)[1:4] <- c("chromosome", "start", "end", "uncertainty")
    weights_bg$chromosome <- as.character(weights_bg$chromosome)
    weights_bg$start <- as.numeric(weights_bg$start)
    weights_bg$uncertainty <- as.numeric(weights_bg$uncertainty)
    # Convert uncertainty to weights: 1 / (uncertainty + epsilon)
    weights_bg$weight <- 1 / (weights_bg$uncertainty + 0.01)
  }

  if (verbose) {
    message("Processing ", length(chromosomes), " chromosome(s) from ", bamFile)
  }
  
  # Process each chromosome
  process_chrom <- function(chrom) {
    if (verbose) {
      message("  Processing chromosome: ", chrom)
    }
    
    # Get binned track
    track_data <- wavematchFromBam(
      bamFile = bamFile,
      genome  = assembly,
      chrom   = chrom,
      binSize = binSize
    )
    
    weights_chr <- NULL
    if (!is.null(weights_bg)) {
      weights_chr_bg <- weights_bg[weights_bg$chromosome == chrom, , drop = FALSE]
      if (nrow(weights_chr_bg) > 0L) {
        weights_chr <- tryCatch({
          cpp_match_weights(
            intervals = track_data$intervals,
            weight_starts = weights_chr_bg$start,
            weight_ends = weights_chr_bg$end,
            weight_values = weights_chr_bg$weight,
            interpolate = TRUE
          )
        }, error = function(e) {
          # Fallback to R implementation
          warning("C++ weight matching failed, using R fallback. Error: ", e$message, immediate. = TRUE)
          matched <- match(track_data$intervals, weights_chr_bg$start)
          weights_chr <- rep(NA_real_, length(track_data$intervals))
          if (any(!is.na(matched))) {
            weights_chr[!is.na(matched)] <- weights_chr_bg$weight[matched[!is.na(matched)]]
            if (any(is.na(weights_chr))) {
              weights_chr[is.na(weights_chr)] <- median(weights_chr, na.rm = TRUE)
            }
          }
          weights_chr
        })
      }
    }
    
    # Call peak caller
    peaks <- wavematchFromTrack(
      chromosome = track_data$chromosome,
      intervals  = track_data$intervals,
      values     = track_data$values,
      weights    = weights_chr,
      ...
    )
    
    peaks
  }
  
  # Parallel or sequential processing
  ncores <- as.integer(ncores)
  if (ncores > 1L && length(chromosomes) > 1L) {
    max_cores <- parallel::detectCores()
    ncores <- min(ncores, max_cores, length(chromosomes))
    
    if (verbose && ncores > 1L) {
      message("Using ", ncores, " cores for parallel processing")
    }
    
    all_peaks <- parallel::mclapply(
      chromosomes,
      process_chrom,
      mc.cores = ncores
    )
  } else {
    all_peaks <- lapply(chromosomes, process_chrom)
  }
  
  # Combine results
  all_peaks <- all_peaks[!vapply(all_peaks, is.null, logical(1L))]
  all_peaks <- all_peaks[vapply(all_peaks, function(x) is.data.frame(x) && nrow(x) > 0L, logical(1L))]
  
  if (length(all_peaks) == 0L) {
    if (returnGRanges) {
      return(GenomicRanges::GRanges())
    }
    return(data.frame(
      chromosome = character(0),
      start      = numeric(0),
      end        = numeric(0),
      stringsAsFactors = FALSE
    ))
  }
  
  # Combine all chromosomes
  result <- do.call(rbind, all_peaks)
  rownames(result) <- NULL
  
  # Apply blacklist filtering if requested
  if (!is.null(blacklist)) {
    if (verbose) {
      message("Filtering peaks using blacklist: ", blacklist)
    }
    result <- filterBlacklist(result, blacklist_file = blacklist)
  }
  
  # Convert to GRanges if requested
  if (returnGRanges) {
    # Get seqinfo if genome name available
    seqinfo <- NULL
    if (!is.null(genome_name)) {
      tryCatch({
        sizes_file <- get_genome_resource_file(genome_name, "sizes")
        chrom_sizes <- get_chrom_sizes_dict(sizes_file)
        seqinfo <- GenomeInfoDb::Seqinfo(
          seqnames = names(chrom_sizes),
          seqlengths = chrom_sizes
        )
      }, error = function(e) {
        # Use default seqinfo
      })
    }
    
    result <- asGRanges(result, keep.mcols = TRUE, seqinfo = seqinfo)
  }
  
  if (verbose) {
    message("Found ", if (returnGRanges) length(result) else nrow(result), " peaks")
  }
  
  result
}
