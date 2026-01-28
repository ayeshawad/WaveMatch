
#' @useDynLib wavematch, .registration = TRUE
#' @importFrom R6 R6Class
#' @importFrom Rsubread featureCounts
#' @importFrom Rcpp evalCpp
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
NULL

#' Wavelet-based peak calling from a BAM file
#'
#' High level helper that:
#' 1) wraps a BAM file as an Alignment object
#' 2) extracts a binned coverage track for one chromosome
#' 3) returns a list with intervals and values
#'
#' @param bamFile path to BAM file
#' @param genome either a genome name such as "hg38" or
#'   an object with a named numeric vector chromSizes
#' @param chrom chromosome name, e.g. "chr1"
#' @param binSize bin size in base pairs
#' @return list with components intervals (numeric) and values (numeric)
#' @export
wavematchFromBam <- function(bamFile,
                                   genome,
                                   chrom,
                                   binSize = 25L) {
  # make sure the BAM exists
  if (!file.exists(bamFile)) {
    stop("BAM file not found: ", bamFile)
  }
  # normalize genome argument to an assembly with $chromSizes
  if (is.character(genome) || is.null(genome$chromSizes)) {
    assembly <- load_genome_resources(genome)
  } else {
    assembly <- genome
  }
  # sanity check for chromSizes
  if (is.null(assembly$chromSizes) || is.null(names(assembly$chromSizes))) {
    stop("assembly$chromSizes must be a named numeric vector")
  }
  if (!chrom %in% names(assembly$chromSizes)) {
    stop("Chromosome ", chrom, " not found in assembly$chromSizes")
  }
  
  aln <- Alignment$new(
    filePath       = bamFile,
    assembly       = assembly,
    pairedEnd      = NULL,
    fragmentLength = NULL,
    scaleFactor    = NULL
  )
  
  binSize <- as.integer(binSize)
  values <- aln$getRawTrack(chrom = chrom, binSize = binSize)

  # bin starts as intervals
  intervals <- seq.int(0L, by = binSize, length.out = length(values))

  list(
    chromosome = chrom,
    intervals  = intervals,
    values     = values
  )
}

#' Wavelet-based peak calling from an existing binned track
#'
#' Convenience wrapper around \code{matchTemplateEmpirical} for users
#' who already have a binned signal track.
#'
#' @param chromosome chromosome name, for example "chr1"
#' @param intervals numeric vector of bin starts in base pairs
#' @param values numeric vector of signal values, same length as \code{intervals}
#' @param base_template numeric vector giving the base template that will be
#'   resampled to the automatically chosen run length in bins.
#'   Defaults to a simple three point derivative shape \code{c(-1, 0, 1)}.
#' @param initLenBins integer lower bound in bins for the automatic minimum
#'   run length selection.
#' @param asinh_threshold optional numeric threshold on the asinh transformed
#'   signal at candidate maxima. If \code{NULL} no explicit threshold is applied.
#' @param minScore optional numeric. If not \code{NULL}, keep only peaks with
#'   template response greater than or equal to this value.
#' @param maxPeaks optional integer. If not \code{NULL}, keep at most this many
#'   peaks, preferring the largest scores.
#' @param iters integer, number of blocks to sample when estimating the
#'   empirical null from local block maxima.
#' @param alpha numeric, FDR threshold on \code{q_value} used to filter matches.
#' @param randSeed optional integer seed for reproducibility of the empirical null.
#' @param smooth_k odd integer kernel size for the running median used when
#'   locating candidate maxima in the response.
#' @param weights optional numeric vector of weights per bin
#' @param mergePeaks logical, if TRUE (default), merge overlapping or nearby peaks
#'   from different template scales using \code{mergeMatches()}.
#' @param mergeGapBP integer, maximum gap in base pairs between peaks to merge them.
#'   Defaults to 150bp, appropriate for broad marks like H3K27ac.
#' @param auto_tune_cascade logical, if TRUE, automatically tune cascade level
#'   based on signal characteristics. Defaults to FALSE.
#' @param verbose logical, if TRUE (default), print progress messages.
#'
#' @return
#' A data frame with one row per detected match and columns:
#' \describe{
#'   \item{chromosome}{chromosome name}
#'   \item{start}{start coordinate in base pairs}
#'   \item{end}{end coordinate in base pairs}
#'   \item{center_idx}{index of the center bin in the input track}
#'   \item{score}{template response at the center}
#'   \item{p_value}{empirical right tail p value}
#'   \item{q_value}{Benjamini Hochberg adjusted p value}
#'   \item{neglog10_p}{minus log10 p value}
#'   \item{neglog10_q}{minus log10 q value}
#'   \item{span_bins}{template span in bins}
#' }
#'
#' If no matches are found an empty data frame with the same columns is returned.
#'
#' @export
wavematchFromTrack <- function(chromosome,
                                     intervals,
                                     values,
                                     base_template   = c(-1, 0, 1),
                                     template_names  = "haar",  # Default: haar (best for H3K27ac and most histone marks)
                                     cascade_levels  = 3L,
                                     use_scaling_function = TRUE,
                                     template_grid_len    = NULL,
                                     initLenBins     = 3L,
                                     asinh_threshold = NULL,
                                     minMatchLengthBP = NULL,
                                     minScore        = NULL,
                                     maxPeaks        = 100000L,
                                     iters           = 20000L,
                                     alpha           = 0.05,
                                     randSeed        = NULL,
                                     smooth_k        = 3L,
                                     null_mask = NULL,
                                     exclude_mask = NULL,
                                     use_split_halves = TRUE,
                                     recenter_at_point_source = TRUE,
                                     add_point_source_bp = TRUE,
                                     weights = NULL,
                                     mergePeaks = TRUE,
                                     mergeGapBP = 150L,
                                     auto_tune_cascade = FALSE,
                                     verbose = TRUE
                                    ) {
  stopifnot(length(intervals) == length(values))

  # Auto-tune cascade level if requested
  if (isTRUE(auto_tune_cascade)) {
    if (verbose) {
      message("Auto-tuning cascade level from signal characteristics...")
    }
    
    # Calculate interval_bp
    interval_bp <- 25L
    if (length(intervals) > 1L) {
      diffs <- diff(intervals)
      diffs <- diffs[is.finite(diffs) & diffs > 0]
      if (length(diffs) > 0L) {
        interval_bp <- as.integer(round(median(diffs, na.rm = TRUE)))
        if (!is.finite(interval_bp) || interval_bp < 1L) {
          interval_bp <- 25L
        }
      }
    }
    
    # Handle NULL cascade_levels: auto-tune for all templates
    if (is.null(cascade_levels)) {
      if (is.null(template_names)) {
        # Single template case: use default template name
        cascade_levels <- autoTuneCascadeLevel(
          values = values,
          interval_bp = interval_bp,
          template_name = "db3"
        )
        if (verbose) {
          message("  Auto-tuned cascade level: ", cascade_levels)
        }
      } else {
        # Multiple templates: tune each one
        cascade_levels <- integer(length(template_names))
        for (i in seq_along(template_names)) {
          tuned_cascade <- autoTuneCascadeLevel(
            values = values,
            interval_bp = interval_bp,
            template_name = template_names[i]
          )
          cascade_levels[i] <- tuned_cascade
          if (verbose) {
            message("  Auto-tuned cascade level for ", template_names[i], ": ", tuned_cascade)
          }
        }
      }
    } else {
      # cascade_levels provided: only override defaults (3L)
      if (is.null(template_names) && length(cascade_levels) == 1L && cascade_levels == 3L) {
        cascade_levels <- autoTuneCascadeLevel(
          values = values,
          interval_bp = interval_bp,
          template_name = "db3"
        )
        if (verbose) {
          message("  Auto-tuned cascade level: ", cascade_levels)
        }
      } else if (!is.null(template_names) && length(cascade_levels) == length(template_names)) {
        # Multiple templates: tune each one that uses default
        for (i in seq_along(template_names)) {
          if (cascade_levels[i] == 3L) {  # Only override if using default
            tuned_cascade <- autoTuneCascadeLevel(
              values = values,
              interval_bp = interval_bp,
              template_name = template_names[i]
            )
            cascade_levels[i] <- tuned_cascade
            if (verbose) {
              message("  Auto-tuned cascade level for ", template_names[i], ": ", tuned_cascade)
            }
          }
        }
      }
    }
  }

  # apply optional per bin weights with a warning if it fails!
  vals_use <- as.numeric(values)

  if (!is.null(weights)) {
    w <- as.numeric(weights)
    if (length(w) != length(vals_use)) {
      warning("`weights` length (", length(w),
              ") does not match `values` length (", length(vals_use),
              "). Ignoring weights.")
    } else {
      vals_use <- vals_use * w
    }
  }

  # option A classic single base_template, no wavelet names
  if (is.null(template_names)) {
    tmpl <- resolve_template_spec(
      base_template         = base_template,
      template_name         = NULL,
      cascade_level         = 3L,
      use_scaling_function  = TRUE,
      grid_len              = NULL,  # Use natural cascade length (no resampling)
      normalize             = TRUE
    )

    out <- matchTemplateEmpirical(
      chromosome       = chromosome,
      intervals        = intervals,
      values           = vals_use,
      base_template    = tmpl,
      initLenBins      = as.integer(initLenBins),
      asinh_threshold  = asinh_threshold,
      minMatchLengthBP = minMatchLengthBP,
      iters            = as.integer(iters),
      alpha            = alpha,
      minScore         = minScore,
      maxPeaks         = maxPeaks,
      randSeed         = randSeed,
      smooth_k         = as.integer(smooth_k),
      null_mask        = null_mask,
      exclude_mask     = exclude_mask,
      use_split_halves = use_split_halves,
      template_name    = NULL,  # Single template, no name
      cascade_level    = 3L,    # Default cascade level
      weights          = weights,
      verbose          = verbose
    )

    if (isTRUE(recenter_at_point_source) && nrow(out) > 0L) {
      out <- recenterAtPointSource(
        peaks               = out,
        intervals           = intervals,
        values              = vals_use,
        add_point_source_bp = add_point_source_bp
      )
    }

    # Merge peaks using simple distance-based merging
    if (isTRUE(mergePeaks) && nrow(out) > 0L) {
      message("    Processing ", nrow(out), " peaks for merging...")
      mergeGapBP <- as.integer(mergeGapBP)
      if (!is.finite(mergeGapBP) || mergeGapBP < 1L) {
        mergeGapBP <- 150L
      }
      # Calculate interval_bp from intervals
      interval_bp <- 25L  # Default
      if (length(intervals) > 1L) {
        interval_bp <- as.integer(round(median(diff(intervals), na.rm = TRUE)))
        if (!is.finite(interval_bp) || interval_bp < 1L) {
          interval_bp <- 25L
        }
      }
      # Simple distance-based merging (no adaptive scale selection)
      out <- mergeMatches(
        out, 
        mergeGapBP = mergeGapBP, 
        adaptive = FALSE,  # No adaptive scale selection
        interval_bp = interval_bp
      )
    }

    return(out)
  }

  # option B wavelet library, possibly multiple templates
  
  template_names <- as.character(template_names)
  n_tmpl <- length(template_names)

  if (length(cascade_levels) == 1L) {
    cascade_levels <- rep.int(as.integer(cascade_levels), n_tmpl)
  }
  if (length(use_scaling_function) == 1L) {
    use_scaling_function <- rep.int(as.logical(use_scaling_function), n_tmpl)
  }

  if (length(cascade_levels) != n_tmpl) {
    stop("`cascade_levels` must have length 1 or match `template_names`.")
  }
  if (length(use_scaling_function) != n_tmpl) {
    stop("`use_scaling_function` must have length 1 or match `template_names`.")
  }

  all_list <- vector("list", n_tmpl)
  any_rows <- FALSE

  for (k in seq_len(n_tmpl)) {
    nm   <- template_names[[k]]
    lev  <- cascade_levels[[k]]
    useS <- use_scaling_function[[k]]

    tmpl_k <- resolve_template_spec(
      base_template         = NULL,
      template_name         = nm,
      cascade_level         = lev,
      use_scaling_function  = useS,
      grid_len              = NULL,  # Use natural cascade length (no resampling)
      normalize             = TRUE
    )

    res_k <- matchTemplateEmpirical(
      chromosome       = chromosome,
      intervals        = intervals,
      values           = vals_use,
      base_template    = tmpl_k,
      initLenBins      = as.integer(initLenBins),
      asinh_threshold  = asinh_threshold,
      minMatchLengthBP = minMatchLengthBP,
      iters            = as.integer(iters),
      alpha            = alpha,
      minScore         = minScore,
      maxPeaks         = maxPeaks,
      randSeed         = randSeed,
      smooth_k         = as.integer(smooth_k),
      null_mask        = null_mask,
      exclude_mask     = exclude_mask,
      use_split_halves = use_split_halves,
      template_name    = nm,  # Pass template name for FDR grouping
      cascade_level    = lev, # Pass cascade level for FDR grouping
      weights          = weights,
      verbose          = verbose
    )

    if (nrow(res_k) > 0L) {
      any_rows <- TRUE
      # template_name already added by matchTemplateEmpirical
      res_k$cascade_level  <- as.integer(lev)
    }

    all_list[[k]] <- res_k
  }

  if (!any_rows) {
    # construct an empty with template columns
    out0 <- data.frame(
      chromosome   = character(0),
      start        = numeric(0),
      end          = numeric(0),
      center_idx   = integer(0),
      score        = integer(0),
      score_raw    = numeric(0),
      p_value      = numeric(0),
      q_value      = numeric(0),
      neglog10_p   = numeric(0),
      neglog10_q   = numeric(0),
      span_bins    = integer(0),
      template_name = character(0),
      cascade_level = integer(0),
      stringsAsFactors = FALSE
    )
    return(out0)
  }

  out <- do.call(rbind, all_list)
  rownames(out) <- NULL

  if (isTRUE(recenter_at_point_source) && nrow(out) > 0L) {
    out <- recenterAtPointSource(
      peaks               = out,
      intervals           = intervals,
      values              = vals_use,
      add_point_source_bp = add_point_source_bp
    )
  }

  # Merge peaks using simple distance-based merging
  if (isTRUE(mergePeaks) && nrow(out) > 0L) {
    message("    Processing ", nrow(out), " peaks for merging...")
    mergeGapBP <- as.integer(mergeGapBP)
    if (!is.finite(mergeGapBP) || mergeGapBP < 1L) {
      mergeGapBP <- 150L
    }
    # Calculate interval_bp from intervals
    interval_bp <- 25L  # Default
    if (length(intervals) > 1L) {
      interval_bp <- as.integer(round(median(diff(intervals), na.rm = TRUE)))
      if (!is.finite(interval_bp) || interval_bp < 1L) {
        interval_bp <- 25L
      }
    }
    # Simple distance-based merging (no adaptive scale selection)
    out <- mergeMatches(
      out, 
      mergeGapBP = mergeGapBP, 
      adaptive = FALSE,  # No adaptive scale selection
      interval_bp = interval_bp
    )
  }

  out
}

#' Wavelet-based peak calling from an existing bedGraph file
#'
#' Convenience wrapper that runs \code{wavematchFromTrack} on each
#' chromosome present in a four column bedGraph file and optionally writes
#' per chromosome and all chromosome narrowPeak like files.
#'
#' The bedGraph is assumed to have columns:
#' \code{chrom}, \code{start}, \code{end}, \code{value},
#' with equal spacing of \code{start} positions within each chromosome.
#'
#' @param weightsBedGraph optional path to uncertainty bedGraph file.
#'   If provided, uncertainty values will be converted to weights and used to
#'   improve peak calling in low-confidence regions.
#' @return A data frame with one row per detected peak over all chromosomes,
#'   in the same format as \code{wavematchFromTrack}, with
#'   optional files written to disk.
#'
#' @export
wavematchFromBedGraph <- function(bedGraphFile,
                                  weightsBedGraph = NULL,
                                  outDir = ".",
                                  outPrefix = NULL,
                                  chromosomes = NULL,
                                  writePerChrom = FALSE,
                                  writeAllChrom = FALSE,
                                  perChromSuffix = "narrowPeak",
                                  allChromFile = NULL,
                                  verbose = TRUE,
                                  base_template   = c(-1, 0, 1),
                                  template_names  = NULL,
                                  cascade_levels  = 3L,
                                  use_scaling_function = TRUE,
                                  template_grid_len    = NULL,
                                  initLenBins     = 3L,
                                  asinh_threshold = NULL,
                                  minMatchLengthBP = NULL,
                                  minScore        = NULL,
                                  maxPeaks        = 100000L,
                                  iters           = 20000L,
                                  alpha           = 0.05,
                                  randSeed        = NULL,
                                  smooth_k        = 3L,
                                  null_mask       = NULL,
                                  exclude_mask    = NULL,
                                  use_split_halves = TRUE,
                                  recenter_at_point_source = TRUE,
                                  add_point_source_bp = TRUE,
                                  weights         = NULL) {

  if (!file.exists(bedGraphFile)) {
    stop("bedGraph file not found: ", bedGraphFile)
  }

  outDir <- normalizePath(outDir, mustWork = FALSE)
  if (!dir.exists(outDir)) {
    dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
  }

  if (is.null(outPrefix)) {
    base <- basename(bedGraphFile)
    outPrefix <- sub("\\.[^.]*$", "", base)
  }

  # read four column bedGraph
  bg <- utils::read.table(
    bedGraphFile,
    header = FALSE,
    sep    = "\t",
    stringsAsFactors = FALSE
  )

  if (ncol(bg) < 4L) {
    stop("bedGraph file must have at least four columns: chrom, start, end, value")
  }

  colnames(bg)[1:4] <- c("chrom", "start", "end", "value")

  bg$chrom <- as.character(bg$chrom)
  bg$start <- as.numeric(bg$start)
  bg$end   <- as.numeric(bg$end)
  bg$value <- as.numeric(bg$value)

  chroms_all <- sort(unique(bg$chrom))

  if (!is.null(chromosomes)) {
    chromosomes <- intersect(as.character(chromosomes), chroms_all)
    if (!length(chromosomes)) {
      stop("None of the requested chromosomes are present in the bedGraph file.")
    }
  } else {
    chromosomes <- chroms_all
  }

  all_list <- vector("list", length(chromosomes))
  names(all_list) <- chromosomes

  # Read weights bedGraph if provided
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
    colnames(weights_bg)[1:4] <- c("chrom", "start", "end", "uncertainty")
    weights_bg$chrom <- as.character(weights_bg$chrom)
    weights_bg$start <- as.numeric(weights_bg$start)
    weights_bg$uncertainty <- as.numeric(weights_bg$uncertainty)
    # Convert uncertainty to weights: 1 / (uncertainty + epsilon)
    weights_bg$weight <- 1 / (weights_bg$uncertainty + 0.01)
  }

  if (verbose) {
    message("wavematchFromBedGraph: processing ", length(chromosomes),
            " chromosome(s) from ", bedGraphFile)
  }

  for (i in seq_along(chromosomes)) {
    chr <- chromosomes[[i]]

    bg_chr <- bg[bg$chrom == chr, , drop = FALSE]
    if (!nrow(bg_chr)) {
      next
    }

    intervals_chr <- bg_chr$start
    values_chr    <- bg_chr$value

    # Get weights for this chromosome
    weights_chr <- NULL
    if (!is.null(weights_bg)) {
      weights_chr_bg <- weights_bg[weights_bg$chrom == chr, , drop = FALSE]
      if (nrow(weights_chr_bg) > 0L) {
        # Match by position (assuming same binning)
        matched <- match(intervals_chr, weights_chr_bg$start)
        if (any(!is.na(matched))) {
          weights_chr <- rep(NA_real_, length(intervals_chr))
          weights_chr[!is.na(matched)] <- weights_chr_bg$weight[matched[!is.na(matched)]]
          # Fill NA with median weight
          if (any(is.na(weights_chr))) {
            weights_chr[is.na(weights_chr)] <- median(weights_chr, na.rm = TRUE)
          }
        }
      }
    }

    if (verbose) {
      message("  chromosome ", chr, " with ", nrow(bg_chr), " bins")
    }

    res_chr <- wavematchFromTrack(
      chromosome           = chr,
      intervals            = intervals_chr,
      values               = values_chr,
      base_template        = base_template,
      template_names       = template_names,
      cascade_levels       = cascade_levels,
      use_scaling_function = use_scaling_function,
      template_grid_len    = template_grid_len,
      initLenBins          = initLenBins,
      asinh_threshold      = asinh_threshold,
      minMatchLengthBP     = minMatchLengthBP,
      minScore             = minScore,
      maxPeaks             = maxPeaks,
      iters                = iters,
      alpha                = alpha,
      randSeed             = randSeed,
      smooth_k             = smooth_k,
      null_mask            = null_mask,
      exclude_mask         = exclude_mask,
      use_split_halves     = use_split_halves,
      recenter_at_point_source = recenter_at_point_source,
      add_point_source_bp  = add_point_source_bp,
      weights              = weights_chr,
      mergePeaks           = mergePeaks,
      mergeGapBP           = mergeGapBP,
      verbose              = verbose
    )

    all_list[[i]] <- res_chr

    if (isTRUE(writePerChrom) && is.data.frame(res_chr) && nrow(res_chr) > 0L) {
      file_chr <- file.path(
        outDir,
        sprintf("%s.%s.%s", outPrefix, chr, perChromSuffix)
      )
      if (verbose) {
        message("    writing ", file_chr)
      }
      .sp_write_narrowPeak(res_chr, file_chr)
    }
  }

  # bind all chromosomes
  res_all <- do.call(rbind, all_list)
  if (is.null(res_all)) {
  if (is.null(template_names)) {
    res_all <- data.frame(
      chromosome   = character(0),
      start        = numeric(0),
      end          = numeric(0),
      center_idx   = integer(0),
      score        = integer(0),
      score_raw    = numeric(0),
      p_value      = numeric(0),
      q_value      = numeric(0),
      neglog10_p   = numeric(0),
      neglog10_q   = numeric(0),
      span_bins    = integer(0),
      stringsAsFactors = FALSE
    )
  } else {
    res_all <- data.frame(
      chromosome    = character(0),
      start         = numeric(0),
      end           = numeric(0),
      center_idx    = integer(0),
      score         = integer(0),
      score_raw     = numeric(0),
      p_value       = numeric(0),
      q_value       = numeric(0),
      neglog10_p    = numeric(0),
      neglog10_q    = numeric(0),
      span_bins     = integer(0),
      template_name = character(0),
      cascade_level = integer(0),
      stringsAsFactors = FALSE
    )
  }
} else {
  rownames(res_all) <- NULL
}
  if (isTRUE(writeAllChrom) && nrow(res_all) > 0L) {
    if (is.null(allChromFile)) {
      allChromFile <- file.path(
        outDir,
        sprintf("%s.allChroms.%s", outPrefix, perChromSuffix)
      )
    }
    if (verbose) {
      message("writing combined file ", allChromFile)
    }
    .sp_write_narrowPeak(res_all, allChromFile)
  }

  res_all
}
