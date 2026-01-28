#' Merge overlapping or nearby structured peaks with template-aware logic
#'
#' This merges peaks detected by the same template scale that are within
#' \code{mergeGapBP} base pairs of each other on the same chromosome.
#' Peaks from different template scales are kept separate unless they
#' significantly overlap, preserving multi-scale detection information.
#'
#' Input is expected to be the output of \code{wavematchFromTrack()}, i.e.
#' a data.frame with at least the columns:
#'   chromosome, start, end, score, neglog10_p, neglog10_q
#'
#' @param df data.frame of peak calls (for example the result of
#'   \code{wavematchFromTrack()}).
#' @param mergeGapBP integer, maximum allowed gap in base pairs between peaks
#'   to consider them part of the same merged region. Defaults to 150.
#'   Can also be a numeric vector for adaptive merging (see details).
#' @param adaptive logical, if TRUE use density-aware adaptive merging.
#'   When TRUE, mergeGapBP is treated as base gap, and actual gap varies
#'   based on local peak density. Defaults to TRUE.
#' @param densityWindowBP integer, window size in bp for calculating local
#'   peak density (used only if adaptive=TRUE). Defaults to 10000.
#' @param minGapBP integer, minimum gap for dense regions (used only if adaptive=TRUE).
#'   Defaults to 50.
#' @param maxGapBP integer, maximum gap for sparse regions (used only if adaptive=TRUE).
#'   Defaults to 300.
#' @param maxPeaksPerRegion integer, maximum peaks per merged region.
#'   If a cluster exceeds this, it will be split. Defaults to 50.
#' @param interval_bp integer, bin size in base pairs (default 25).
#'   Used for scale-dependent merge distance calculation.
#' @param require_same_template logical, if TRUE only merge peaks from same template.
#'   Cross-scale merging requires significant overlap. Defaults to TRUE.
#' @param use_summit_clustering logical, if TRUE use summit distance for clustering.
#'   Defaults to TRUE.
#' @param null_median numeric, median of null distribution (optional).
#' @param null_mad numeric, MAD of null distribution (optional).
#' @param null_mean numeric, mean of null distribution (optional).
#' @param null_sd numeric, SD of null distribution (optional).
#' @param scale_map integer vector, adaptive scale map (scale index per bin).
#'   If provided, peaks are filtered to keep only those from the "correct" scale.
#' @param scale_L_bins integer vector, template bin lengths for each adaptive scale.
#' @param intervals numeric vector, genomic intervals (bin start positions).
#' @param signal_values numeric vector, original signal values for valley detection.
#'   If provided, valley detection will be used instead of distance-based merging.
#'
#' @return
#' A data.frame with one row per merged region and columns:
#'   chromosome, start, end, name, score, strand, signal,
#'   neglog10_p, neglog10_q, pointSource, contributing_templates, n_templates.
#' @export
mergeMatches <- function(df, 
                        mergeGapBP = 150L,
                        adaptive = TRUE,
                        densityWindowBP = 10000L,
                        minGapBP = 50L,
                        maxGapBP = 300L,
                        maxPeaksPerRegion = 50L,
                        interval_bp = 25L,
                        require_same_template = TRUE,
                        use_summit_clustering = TRUE,
                        null_median = NULL,
                        null_mad = NULL,
                        null_mean = NULL,
                        null_sd = NULL,
                        scale_map = NULL,
                        scale_L_bins = NULL,
                        intervals = NULL,
                        signal_values = NULL) {
  if (missing(df) || !is.data.frame(df)) {
    stop("`df` must be a data.frame of peak calls")
  }

  required_cols <- c("chromosome", "start", "end", "score", "neglog10_p", "neglog10_q")
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols)) {
    stop(
      "Input `df` is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  if (nrow(df) == 0L) {
    return(
      data.frame(
        chromosome  = character(0),
        start       = numeric(0),
        end         = numeric(0),
        name        = character(0),
        score       = integer(0),
        strand      = character(0),
        signal      = numeric(0),
        neglog10_p  = numeric(0),
        neglog10_q  = numeric(0),
        pointSource = integer(0),
        contributing_templates = character(0),
        n_templates = integer(0),
        stringsAsFactors = FALSE
      )
    )
  }

  MAX_NEGLOGP <- 10.0
  MIN_NEGLOGP <- 1.0e-10

  df_local <- df
  
  # Use score as signal if signal column is missing
  if (!"signal" %in% colnames(df_local)) {
    df_local$signal <- df_local$score
  }

  # Add span_bins if not present (for scale-dependent merging)
  if (!"span_bins" %in% colnames(df_local)) {
    # Estimate from peak width
    df_local$span_bins <- as.integer(round((df_local$end - df_local$start) / interval_bp))
    df_local$span_bins <- pmax(3L, df_local$span_bins)  # Minimum 3 bins
  }

  # Compute null distribution stats from input data if not provided
  if (is.null(null_median) || is.null(null_mad)) {
    if ("score_raw" %in% colnames(df_local)) {
      raw_responses <- df_local$score_raw[is.finite(df_local$score_raw)]
      if (length(raw_responses) > 10L) {
        null_median <- median(raw_responses, na.rm = TRUE)
        null_mad <- mad(raw_responses, na.rm = TRUE)
        null_mean <- mean(raw_responses, na.rm = TRUE)
        null_sd <- sd(raw_responses, na.rm = TRUE)
      }
    } else if ("signal" %in% colnames(df_local)) {
      signal_vals <- df_local$signal[is.finite(df_local$signal)]
      if (length(signal_vals) > 10L) {
        null_median <- median(signal_vals, na.rm = TRUE)
        null_mad <- mad(signal_vals, na.rm = TRUE)
        null_mean <- mean(signal_vals, na.rm = TRUE)
        null_sd <- sd(signal_vals, na.rm = TRUE)
      }
    }
  }

  # STEP 2: Valley-based merging (replaces distance-based merging)
  # Use valley detection and signal-to-valley ratio to determine if peaks should merge
  if (nrow(df_local) == 0L) {
    return(df_local)
  }
  
  # For backward compatibility, keep mergeGapBP but it's only used as fallback
  # Primary decision is based on valley detection
  if (isTRUE(adaptive)) {
    mergeGapBP <- computeAdaptiveMergeGaps(
      df_local, 
      baseGapBP = as.integer(mergeGapBP[1L]),
      densityWindowBP = as.integer(densityWindowBP),
      minGapBP = as.integer(minGapBP),
      maxGapBP = as.integer(maxGapBP)
    )
  } else {
    mergeGapBP <- rep(as.integer(mergeGapBP[1L]), nrow(df_local))
  }

  # basic types
  df_local$chromosome  <- as.character(df_local$chromosome)
  df_local$start       <- as.numeric(df_local$start)
  df_local$end         <- as.numeric(df_local$end)
  df_local$score       <- as.numeric(df_local$score)
  df_local$signal      <- as.numeric(df_local$signal)
  df_local$neglog10_p  <- as.numeric(df_local$neglog10_p)
  df_local$neglog10_q  <- as.numeric(df_local$neglog10_q)
  df_local$span_bins   <- as.integer(df_local$span_bins)

  # sort by chromosome, then template_name (if available), then start
  if ("template_name" %in% colnames(df_local)) {
    ord <- order(df_local$chromosome, 
                 df_local$template_name,
                 df_local$start, 
                 df_local$end)
  } else {
  ord <- order(df_local$chromosome, df_local$start, df_local$end)
  }
  df_sorted <- df_local[ord, , drop = FALSE]
  mergeGapBP_sorted <- mergeGapBP[ord]
  rownames(df_sorted) <- NULL

  n <- nrow(df_sorted)

  # ROBUST CLUSTERING WITH TEMPLATE-AWARE LOGIC
  cluster_id <- integer(n)
  current_cluster <- 1L
  cluster_id[1L] <- current_cluster

  for (i in 2L:n) {
    peak1 <- df_sorted[i - 1L, , drop = FALSE]
    peak2 <- df_sorted[i, , drop = FALSE]
    
    # Check if peaks should be merged using valley detection
    should_merge <- shouldMergePeaks(
      peak1 = peak1,
      peak2 = peak2,
      mergeGapBP = mergeGapBP_sorted[i],
      interval_bp = interval_bp,
      require_same_template = require_same_template,
      use_summit_clustering = use_summit_clustering,
      signal_values = signal_values,
      intervals = intervals,
      null_median = null_median,
      null_mad = null_mad
    )
    
    if (!should_merge) {
      current_cluster <- current_cluster + 1L
    }
    cluster_id[i] <- current_cluster
  }

  df_sorted$cluster_id <- cluster_id

  # split by cluster and aggregate
  split_clusters <- split(df_sorted, df_sorted$cluster_id)
  n_clusters <- length(split_clusters)
  
  if (n_clusters > 1000L) {
    message("      Processing ", n_clusters, " clusters...")
  }

  out_list <- vector("list", length(split_clusters))
  idx_out <- 0L
  i_global <- 0L

  for (gdf in split_clusters) {
    if (nrow(gdf) == 0L) next

    # OPTIMIZATION: Split large clusters if they exceed maxPeaksPerRegion
    if (nrow(gdf) > maxPeaksPerRegion) {
      sub_clusters <- splitLargeCluster(gdf, maxPeaksPerRegion)
      for (sub_gdf in sub_clusters) {
        if (nrow(sub_gdf) == 0L) next
        
        i_global <- i_global + 1L
        merged_peak <- aggregateCluster(
          sub_gdf, 
          i_global, 
          mergeGapBP_sorted[1L], 
          MAX_NEGLOGP, 
          MIN_NEGLOGP,
          null_median = null_median,
          null_mad = null_mad,
          null_mean = null_mean,
          null_sd = null_sd
        )
        idx_out <- idx_out + 1L
        out_list[[idx_out]] <- merged_peak
      }
    } else {
    i_global <- i_global + 1L
      merged_peak <- aggregateCluster(
        gdf, 
        i_global, 
        mergeGapBP_sorted[1L], 
        MAX_NEGLOGP, 
        MIN_NEGLOGP,
        null_median = null_median,
        null_mad = null_mad,
        null_mean = null_mean,
        null_sd = null_sd
      )
      idx_out <- idx_out + 1L
      out_list[[idx_out]] <- merged_peak
    }
  }

  if (idx_out == 0L) {
    return(
      data.frame(
        chromosome  = character(0),
        start       = numeric(0),
        end         = numeric(0),
        name        = character(0),
        score       = integer(0),
        strand      = character(0),
        signal      = numeric(0),
        neglog10_p  = numeric(0),
        neglog10_q  = numeric(0),
        pointSource = integer(0),
        contributing_templates = character(0),
        n_templates = integer(0),
        stringsAsFactors = FALSE
      )
    )
  }

  out <- do.call(rbind, out_list[seq_len(idx_out)])

  out <- out[order(out$chromosome, out$start, out$end), , drop = FALSE]
  rownames(out) <- NULL
  out
}

# Helper function: Determine if two peaks should be merged
# Implements template-aware, scale-dependent, summit-based merging logic
shouldMergePeaks <- function(peak1, peak2, mergeGapBP, interval_bp = 25L,
                             require_same_template = TRUE,
                             use_summit_clustering = TRUE,
                             signal_values = NULL,
                             intervals = NULL,
                             null_median = NULL,
                             null_mad = NULL) {
  
  # 1. Same chromosome?
  if (peak1$chromosome[1L] != peak2$chromosome[1L]) {
    return(FALSE)
  }
  
  # 2. Template-aware merging
  template1 <- if ("template_name" %in% colnames(peak1)) peak1$template_name[1L] else NA
  template2 <- if ("template_name" %in% colnames(peak2)) peak2$template_name[1L] else NA
  
  same_template <- FALSE
  if (require_same_template) {
    if (!is.na(template1) && !is.na(template2)) {
      same_template <- template1 == template2
    } else if (is.na(template1) && is.na(template2)) {
      same_template <- TRUE  # Both NA = treat as same
    } else {
      same_template <- FALSE  # One NA, one not = different
    }
  } else {
    same_template <- TRUE  # If not requiring same template, allow merging
  }
  
  # 3. Cross-scale merging: require significant overlap
  if (!same_template && require_same_template) {
    overlap <- min(peak1$end[1L], peak2$end[1L]) - max(peak1$start[1L], peak2$start[1L])
    if (overlap <= 0) {
      return(FALSE)  # No overlap = separate peaks
    }
    
    # Require >80% overlap of smaller peak for cross-scale merging
    width1 <- peak1$end[1L] - peak1$start[1L]
    width2 <- peak2$end[1L] - peak2$start[1L]
    smaller_width <- min(width1, width2)
    
    if (overlap < 0.8 * smaller_width) {
      return(FALSE)  # Insufficient overlap = separate peaks
    }
    # If we get here, peaks overlap significantly - allow merge
    return(TRUE)
  }
  
  # 4. VALLEY-BASED MERGING (replaces distance-based)
  # Use valley detection and signal-to-valley ratio for mathematically principled merging
  if (!is.null(signal_values) && !is.null(intervals) && 
      length(signal_values) == length(intervals)) {
    
    valley_result <- checkValleyBetweenPeaks(
      peak1 = peak1,
      peak2 = peak2,
      signal_values = signal_values,
      intervals = intervals,
      null_median = null_median,
      null_mad = null_mad
    )
    
    # Valley detection is definitive: if valley exists, peaks are distinct
    if (!is.null(valley_result) && !valley_result$should_merge) {
      return(FALSE)  # Distinct peaks detected
    }
    
    # If valley check says merge, proceed (but still check distance as fallback)
    if (!is.null(valley_result) && valley_result$should_merge) {
      # Additional check: if peaks are very far apart, don't merge even if no valley
      gap <- peak2$start[1L] - peak1$end[1L]
      if (gap > mergeGapBP * 2.0) {
        return(FALSE)  # Too far apart
      }
      return(TRUE)  # No valley, close together = merge
    }
  }
  
  # 5. FALLBACK: Distance-based merging (if valley detection unavailable)
  # Only used if signal_values not provided
  gap <- peak2$start[1L] - peak1$end[1L]
  
  # If peaks are overlapping, merge
  if (gap <= 0) {
    return(TRUE)
  }
  
  # If peaks are too far apart, don't merge
  if (gap > mergeGapBP) {
    return(FALSE)
  }
  
  # For close peaks without valley info, be conservative: don't merge
  # (Better to have duplicate peaks than merge distinct ones)
  return(FALSE)
}

# Valley detection function: mathematically principled peak separation
# Uses fast C++ implementation with binary search for O(log n) performance
# Returns list with should_merge (logical) and valley_ratio (numeric)
checkValleyBetweenPeaks <- function(peak1, peak2, signal_values, intervals,
                                    null_median = NULL, null_mad = NULL) {
  
  # Extract peak boundaries
  peak1_start <- as.numeric(peak1$start[1L])
  peak1_end <- as.numeric(peak1$end[1L])
  peak2_start <- as.numeric(peak2$start[1L])
  peak2_end <- as.numeric(peak2$end[1L])
  
  # Get peak signals (use signal column if available, else score)
  peak1_signal <- if ("signal" %in% colnames(peak1)) as.numeric(peak1$signal[1L]) else as.numeric(peak1$score[1L])
  peak2_signal <- if ("signal" %in% colnames(peak2)) as.numeric(peak2$signal[1L]) else as.numeric(peak2$score[1L])
  
  # Convert null stats to numeric (handle NULL)
  null_median_val <- if (!is.null(null_median)) as.numeric(null_median) else NA_real_
  null_mad_val <- if (!is.null(null_mad)) as.numeric(null_mad) else NA_real_
  
  # Call fast C++ implementation with binary search
  result <- tryCatch({
    wavematch:::cpp_check_valley_between_peaks(
      peak1_start = peak1_start,
      peak1_end = peak1_end,
      peak2_start = peak2_start,
      peak2_end = peak2_end,
      peak1_signal = peak1_signal,
      peak2_signal = peak2_signal,
      signal_values = signal_values,
      intervals = intervals,
      null_median = null_median_val,
      null_mad = null_mad_val
    )
  }, error = function(e) {
    # Fallback to R implementation if C++ fails
    warning("C++ valley detection failed, using R fallback. Error: ", e$message, immediate. = TRUE)
    return(checkValleyBetweenPeaks_R(peak1, peak2, signal_values, intervals, null_median, null_mad))
  })
  
  return(result)
}

# R fallback implementation (kept for compatibility)
checkValleyBetweenPeaks_R <- function(peak1, peak2, signal_values, intervals,
                                     null_median = NULL, null_mad = NULL) {
  
  # Get peak boundaries
  peak1_start <- peak1$start[1L]
  peak1_end <- peak1$end[1L]
  peak2_start <- peak2$start[1L]
  peak2_end <- peak2$end[1L]
  
  # If peaks overlap significantly, they're likely the same
  overlap <- min(peak1_end, peak2_end) - max(peak1_start, peak2_start)
  if (overlap > 0) {
    width1 <- peak1_end - peak1_start
    width2 <- peak2_end - peak2_start
    smaller_width <- min(width1, width2)
    if (overlap > 0.5 * smaller_width) {
      return(list(should_merge = TRUE, valley_ratio = 0.0))
    }
  }
  
  # Get valley region
  valley_start_bp <- peak1_end
  valley_end_bp <- peak2_start
  if (valley_start_bp > valley_end_bp) {
    valley_start_bp <- peak2_start
    valley_end_bp <- peak1_end
  }
  
  # Use binary search for efficiency (O(log n) instead of O(n))
  valley_start_idx <- findInterval(valley_start_bp, intervals, rightmost.closed = FALSE)
  valley_end_idx <- findInterval(valley_end_bp, intervals, rightmost.closed = TRUE)
  
  if (valley_start_idx > valley_end_idx || valley_start_idx > length(intervals) || valley_end_idx < 1L) {
    gap <- peak2_start - peak1_end
    if (gap <= 0) {
      return(list(should_merge = TRUE, valley_ratio = 0.0))
    }
    return(list(should_merge = FALSE, valley_ratio = NA_real_))
  }
  
  # Extract only valley signals (efficient)
  valley_indices <- valley_start_idx:valley_end_idx
  valley_indices <- valley_indices[valley_indices >= 1L & valley_indices <= length(signal_values)]
  valley_signals <- signal_values[valley_indices]
  valley_signals <- valley_signals[is.finite(valley_signals)]
  
  if (length(valley_signals) == 0L) {
    return(list(should_merge = FALSE, valley_ratio = NA_real_))
  }
  
  valley_min <- min(valley_signals, na.rm = TRUE)
  valley_median <- median(valley_signals, na.rm = TRUE)
  
  peak1_signal <- if ("signal" %in% colnames(peak1)) peak1$signal[1L] else peak1$score[1L]
  peak2_signal <- if ("signal" %in% colnames(peak2)) peak2$signal[1L] else peak2$score[1L]
  
  peak_max_signal <- max(peak1_signal, peak2_signal, na.rm = TRUE)
  peak_avg_signal <- (peak1_signal + peak2_signal) / 2.0
  
  valley_depth_ratio <- if (peak_max_signal > 0 && valley_min < peak_max_signal) {
    (peak_max_signal - valley_min) / peak_max_signal
  } else {
    0.0
  }
  
  valley_above_null <- TRUE
  if (!is.null(null_median) && !is.null(null_mad) && null_mad > 1e-10) {
    null_threshold <- null_median - 2.0 * null_mad
    valley_above_null <- valley_min > null_threshold
  } else if (!is.null(null_median)) {
    valley_above_null <- valley_min > null_median
  }
  
  # Decision logic (same as C++ version)
  if (valley_depth_ratio > 0.5) {
    return(list(should_merge = FALSE, valley_ratio = valley_depth_ratio))
  }
  if (!valley_above_null && valley_depth_ratio > 0.3) {
    return(list(should_merge = FALSE, valley_ratio = valley_depth_ratio))
  }
  if (valley_depth_ratio < 0.3 && valley_above_null) {
    return(list(should_merge = TRUE, valley_ratio = valley_depth_ratio))
  }
  if (length(valley_signals) >= 3L && peak_avg_signal > 0) {
    valley_to_peak_ratio <- valley_median / peak_avg_signal
    if (valley_to_peak_ratio < 0.5) {
      return(list(should_merge = FALSE, valley_ratio = valley_depth_ratio))
    }
  }
  
  return(list(should_merge = FALSE, valley_ratio = valley_depth_ratio))
}

# Helper function: Get scale-dependent merge distance
getScaleDependentMergeDistance <- function(span_bins1, span_bins2, baseGapBP, interval_bp = 25L) {
  # Convert span_bins to bp
  width1_bp <- span_bins1 * interval_bp
  width2_bp <- span_bins2 * interval_bp
  
  # Check if scales are similar (within 2 bins)
  scale_diff <- abs(span_bins1 - span_bins2)
  
  if (scale_diff <= 2L) {
    # Same scale: use average width as merge distance (50% of average width)
    avg_width <- mean(c(width1_bp, width2_bp))
    return(max(baseGapBP, avg_width * 0.5))
  } else {
    # Different scales: require overlap (negative = overlap required)
    return(-1L)
  }
}

# Helper function to compute adaptive merge gaps based on local density
computeAdaptiveMergeGaps <- function(df, baseGapBP, densityWindowBP, minGapBP, maxGapBP) {
  n <- nrow(df)
  if (n == 0L) return(integer(0))
  
  ord <- order(df$chromosome, df$start)
  df_sorted <- df[ord, , drop = FALSE]
  
  gaps <- integer(n)
  
  for (i in seq_len(n)) {
    chrom <- df_sorted$chromosome[i]
    center <- (df_sorted$start[i] + df_sorted$end[i]) / 2.0
    
    window_start <- center - densityWindowBP / 2.0
    window_end <- center + densityWindowBP / 2.0
    
    in_window <- (df_sorted$chromosome == chrom) &
                 ((df_sorted$start + df_sorted$end) / 2.0 >= window_start) &
                 ((df_sorted$start + df_sorted$end) / 2.0 <= window_end)
    
    density <- sum(in_window)
    
    if (density <= 1L) {
      gaps[i] <- maxGapBP
    } else if (density >= 20L) {
      gaps[i] <- minGapBP
    } else {
      density_norm <- (density - 1.0) / (20.0 - 1.0)
      gaps[i] <- as.integer(round(maxGapBP - density_norm * (maxGapBP - minGapBP)))
    }
  }
  
  gaps <- pmax(minGapBP, pmin(maxGapBP, gaps))
  gaps[order(ord)]
}

# Helper function to split large clusters
splitLargeCluster <- function(gdf, maxPeaks) {
  n_peaks <- nrow(gdf)
  if (n_peaks <= maxPeaks) {
    return(list(gdf))
  }
  
  ord <- order(gdf$start)
  gdf_sorted <- gdf[ord, , drop = FALSE]
  
  gaps <- diff(gdf_sorted$start)
  median_gap <- median(gaps, na.rm = TRUE)
  large_gaps <- which(gaps > median_gap * 2.0)
  
  if (length(large_gaps) == 0L) {
    n_splits <- ceiling(n_peaks / maxPeaks)
    split_size <- ceiling(n_peaks / n_splits)
    splits <- list()
    for (i in seq_len(n_splits)) {
      start_idx <- (i - 1L) * split_size + 1L
      end_idx <- min(i * split_size, n_peaks)
      splits[[i]] <- gdf_sorted[start_idx:end_idx, , drop = FALSE]
    }
    return(splits)
  } else {
    split_points <- c(0L, large_gaps, n_peaks)
    splits <- list()
    for (i in seq_len(length(split_points) - 1L)) {
      start_idx <- split_points[i] + 1L
      end_idx <- split_points[i + 1L]
      if (start_idx <= end_idx) {
        splits[[i]] <- gdf_sorted[start_idx:end_idx, , drop = FALSE]
      }
    }
    return(splits)
  }
}

# Helper function to aggregate a cluster into a single merged peak
aggregateCluster <- function(gdf, i_global, mergeGapBP, MAX_NEGLOGP, MIN_NEGLOGP,
                             null_median=NULL, null_mad=NULL, null_mean=NULL, null_sd=NULL) {
    chrom <- gdf$chromosome[1L]
    sMin  <- min(gdf$start)
    eMax  <- max(gdf$end)

    scSum  <- sum(gdf$score,  na.rm = TRUE)
    sigSum <- sum(gdf$signal, na.rm = TRUE)
    n_peaks <- nrow(gdf)

    maxS    <- -Inf
    peakAbs <- -1
    
    # Check if input has point_source column (from matchTemplateEmpirical)
    # point_source is relative position within peak (relative to peak start)
    has_point_source <- "point_source" %in% colnames(gdf)
    if (has_point_source) {
      # Convert relative point_source to absolute genomic positions
      valid_ps <- is.finite(gdf$point_source) & is.finite(gdf$start)
      if (any(valid_ps)) {
        point_sources_abs <- gdf$start[valid_ps] + gdf$point_source[valid_ps]
        # Find point source with maximum signal
        valid_signal <- is.finite(gdf$signal[valid_ps])
        if (any(valid_signal)) {
          best_idx <- which.max(gdf$signal[valid_ps][valid_signal])
          peakAbs <- point_sources_abs[valid_signal][best_idx]
        }
      }
    }

    pMax    <- -Inf
    pTail   <- 0.0
    pHasInf <- FALSE

    qMax    <- -Inf
    qMin    <- Inf
    qTail   <- 0.0
    qHasInf <- FALSE
  
  # Track contributing templates
  contributing_templates <- character(0)
  if ("template_name" %in% colnames(gdf)) {
    unique_templates <- unique(gdf$template_name[!is.na(gdf$template_name)])
    contributing_templates <- paste(sort(unique_templates), collapse = ",")
    if (contributing_templates == "") contributing_templates <- "unknown"
  } else {
    contributing_templates <- "unknown"
  }
  n_templates <- length(unique(gdf$template_name[!is.na(gdf$template_name)]))
  if (n_templates == 0L) n_templates <- 1L

    for (j in seq_len(n_peaks)) {
      start_j <- gdf$start[j]
      end_j   <- gdf$end[j]
      score_j <- gdf$score[j]
      signal_j <- gdf$signal[j]
      pLog10  <- gdf$neglog10_p[j]
      qLog10  <- gdf$neglog10_q[j]

      # track best signal for point source
      if (is.finite(signal_j) && signal_j > maxS) {
        maxS <- signal_j
        peakAbs <- floor((start_j + end_j) / 2)
      }

      # p log10 combination
      if (is.infinite(pLog10) || (!is.na(pLog10) && pLog10 >= MAX_NEGLOGP)) {
        pHasInf <- TRUE
      } else if (is.finite(pLog10)) {
        if (pLog10 > pMax) {
          if (is.infinite(pMax)) {
            pTail <- 1.0
          } else {
            pTail <- pTail * (10^(pMax - pLog10)) + 1.0
          }
          pMax <- pLog10
        } else {
          pTail <- pTail + 10^(pLog10 - pMax)
        }
      }

      # q log10 combination
      if (is.infinite(qLog10) ||
          (!is.na(qLog10) && (qLog10 >= MAX_NEGLOGP || qLog10 <= MIN_NEGLOGP))) {
        qHasInf <- TRUE
      } else if (is.finite(qLog10)) {
        if (qLog10 < qMin) {
          qMin <- if (qLog10 < MIN_NEGLOGP) MIN_NEGLOGP else qLog10
        }
        if (qLog10 > qMax) {
          if (is.infinite(qMax)) {
            qTail <- 1.0
          } else {
            qTail <- qTail * (10^(qMax - qLog10)) + 1.0
          }
          qMax <- qLog10
        } else {
          qTail <- qTail + 10^(qLog10 - qMax)
        }
      }
    }

  # RECOMPUTE SCORE FROM MERGED STATISTICS
  if (qHasInf) {
    merged_q <- 1e-10
  } else {
    merged_q <- 10^(-qMax) * qTail
    merged_q <- pmax(1e-10, pmin(1.0, merged_q))
  }
  
  merged_neglog10_q <- -log10(merged_q)
  merged_significance <- pmin(10, merged_neglog10_q * 3.33)
  
  # Compute merged effect size from strongest signal (maxS)
  if (!is.null(null_median) && !is.null(null_mad) && null_mad > 1e-10) {
    merged_effect_size <- (maxS - null_median) / max(null_mad, 1e-10)
  } else if (!is.null(null_mean) && !is.null(null_sd) && null_sd > 1e-10) {
    merged_effect_size <- (maxS - null_mean) / max(null_sd, 1e-10)
  } else {
    merged_effect_size <- 2.0
  }
  
  merged_effect_norm <- pmax(0, pmin(10, merged_effect_size))
  merged_combined <- 0.7 * merged_significance + 0.3 * merged_effect_norm
  
  score_min_expected <- 2.0
  score_max_expected <- 8.0
  
  if (merged_combined < score_min_expected) score_min_expected <- merged_combined
  if (merged_combined > score_max_expected) score_max_expected <- merged_combined
  
  score_range_expected <- score_max_expected - score_min_expected
  if (score_range_expected > 1e-10) {
    avgScore <- 100.0 + 900.0 * (merged_combined - score_min_expected) / score_range_expected
  } else {
    avgScore <- 500.0
  }
  
  if (!is.finite(avgScore)) avgScore <- 500.0
  avgScore <- max(0.0, min(1000.0, avgScore))
    scoreInt <- as.integer(round(avgScore))

    sigAvg <- sigSum / max(n_peaks, 1L)
    if (!is.finite(sigAvg)) sigAvg <- 0

  # harmonic like aggregation for p and q in log10 scale
    if (pHasInf) {
      pHMLog10 <- MAX_NEGLOGP
    } else if (is.infinite(pMax) || !(pTail > 0.0) || is.na(pTail)) {
      pHMLog10 <- MIN_NEGLOGP
    } else {
      pHMLog10 <- -log10(n_peaks) + (pMax + log10(pTail))
      pHMLog10 <- max(MIN_NEGLOGP, min(pHMLog10, MAX_NEGLOGP))
    }

    if (qHasInf) {
      qHMLog10 <- MAX_NEGLOGP
  } else if (is.infinite(qMax) || !is.finite(qTail) || !(qTail > 0.0) || is.na(qTail)) {
      qHMLog10 <- MIN_NEGLOGP
    } else {
      qHMLog10 <- -log10(n_peaks) + (qMax + log10(qTail))
      qHMLog10 <- max(MIN_NEGLOGP, min(qHMLog10, MAX_NEGLOGP))
    }

    qMinLog10 <- qMin
    qMaxLog10 <- qMax
    if (is.finite(qMinLog10) && qMinLog10 < MIN_NEGLOGP) {
      qMinLog10 <- MIN_NEGLOGP
    }
    if (is.finite(qMaxLog10) && qMaxLog10 > MAX_NEGLOGP) {
      qMaxLog10 <- MAX_NEGLOGP
    } else if (!is.finite(qMaxLog10) || !is.finite(qMinLog10) || qMaxLog10 < MIN_NEGLOGP) {
      qMinLog10 <- 0.0
      qMaxLog10 <- 0.0
    }

    if (peakAbs >= 0) {
      pointSource <- as.integer(peakAbs - sMin)
    } else {
      pointSource <- as.integer((eMax - sMin) %/% 2)
    }

    name <- sprintf(
    "structuredPeak|i=%d|gap=%dbp|ct=%d|qRange=%.3f_%.3f|tmpls=%s",
      i_global,
      mergeGapBP,
      n_peaks,
      qMinLog10,
    qMaxLog10,
    contributing_templates
    )

  data.frame(
      chromosome  = chrom,
      start       = as.numeric(sMin),
      end         = as.numeric(eMax),
      name        = name,
      score       = as.integer(scoreInt),
      strand      = ".",
      signal      = as.numeric(sigAvg),
      neglog10_p  = as.numeric(pHMLog10),
      neglog10_q  = as.numeric(qHMLog10),
      pointSource = as.integer(pointSource),
    contributing_templates = contributing_templates,
    n_templates = as.integer(n_templates),
      stringsAsFactors = FALSE
    )
}
