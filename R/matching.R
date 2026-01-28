# C++ functions will be available after Rcpp::compileAttributes() is run
# They are defined in src/cpp_template_matching.cpp and exported via RcppExports

#' Relative maxima in a numeric sequence
#'
#' Find indices of local maxima such that each maximum is the highest
#' value within a window of +/- orderBins. Groups contiguous tied maxima
#' and selects the center of each group.
#'
#' @param x numeric vector
#' @param orderBins integer, half window in bins
#' @param eps numeric, tolerance for floating point comparisons (default: machine epsilon * 10)
#'
#' @return integer vector of indices (1-based) of local maxima
#' @export
relativeMaxima <- function(x, orderBins, eps = NULL) {
  x <- as.numeric(x)
  n <- length(x)
  if (n == 0L) return(integer(0L))
  
  orderBins <- max(as.integer(orderBins), 1L)
  
  # Default tolerance: machine epsilon * 10 for floating point comparisons
  if (is.null(eps)) {
    eps <- .Machine$double.eps * 10
  }
  eps <- as.numeric(eps)
  
  # Use optimized C++ implementation
  idx <- tryCatch({
    cpp_relative_maxima(x, orderBins, eps)
  }, error = function(e) {
    # Fallback to R implementation if C++ not available
    out <- integer(0L)
    for (i in seq_len(n)) {
      s <- max(1L, i - orderBins)
      e <- min(n, i + orderBins)
      v <- x[s:e]
      if (!is.finite(x[i])) next
      # Use tolerance-aware comparison
      if (x[i] >= (max(v, na.rm = TRUE) - eps) && x[i] > (min(v, na.rm = TRUE) + eps)) {
        out <- c(out, i)
      }
    }
    out
  })
  
  # Handle ties: group contiguous tied maxima and pick center
  # This ensures that when multiple consecutive positions have the same value,
  # we select the center of the tied group rather than multiple separate maxima
  if (length(idx) > 1L && eps > 0.0) {
    # Sort indices to ensure they're in order (C++ might return unsorted)
    idx <- sort(idx)
    
    # Fill in gaps: if C++ returned non-contiguous indices but intermediate
    # positions have the same value (within tolerance), include them
    idx_expanded <- idx
    for (i in 1L:(length(idx) - 1L)) {
      gap_start <- idx[i]
      gap_end <- idx[i + 1L]
      if (gap_end > gap_start + 1L) {
        # Check if values in gap are equal (within tolerance)
        val_start <- x[gap_start]
        val_end <- x[gap_end]
        if (abs(val_start - val_end) <= eps) {
          # Check intermediate positions
          gap_positions <- (gap_start + 1L):(gap_end - 1L)
          gap_values <- x[gap_positions]
          # If all gap values equal start/end value (within tolerance), include them
          if (all(abs(gap_values - val_start) <= eps, na.rm = TRUE)) {
            idx_expanded <- c(idx_expanded, gap_positions)
          }
        }
      }
    }
    idx <- sort(unique(idx_expanded))
    
    groups <- list()
    start <- idx[1L]
    prev <- idx[1L]
    
    # Group contiguous indices
    for (i in 2L:length(idx)) {
      if (idx[i] == prev + 1L) {
        # Still contiguous
        prev <- idx[i]
      } else {
        # Gap found: break off previous group
        groups[[length(groups) + 1L]] <- c(start, prev)
        start <- idx[i]
        prev <- idx[i]
      }
    }
    # Add final group
    groups[[length(groups) + 1L]] <- c(start, prev)
    
    # Pick center of each group
    centers <- integer(length(groups))
    for (i in seq_along(groups)) {
      s <- groups[[i]][1L]
      e <- groups[[i]][2L]
      if (s == e) {
        centers[i] <- s
      } else {
        # Center of tied group: use integer division to get center
        centers[i] <- (s + e) %/% 2L
      }
    }
    
    return(centers)
  }
  
  # If only one or zero maxima, return as-is
  return(idx)
}

#' Benjamini-Hochberg FDR adjustment
#'
#' Simple BH implementation returning adjusted p values (q values).
#' @param p numeric vector of raw p values
#' @return numeric vector of same length with FDR-adjusted p values
#' @keywords internal
bh_fdr <- function(p) {
  p <- as.numeric(p)
  n <- length(p)
  if (n == 0L) return(p)
  
  o <- order(p, na.last = TRUE)
  ro <- order(o)
  ranked <- seq_len(n)
  
  adj <- p[o] * n / ranked
  adj <- rev(cummin(rev(adj)))
  adj <- pmin(adj, 1)
  
  adj[ro]
}

#' Parse asinh threshold specification
#'
#' Helper to turn a user threshold specification into
#' asinh transformed values, and a numeric cutoff on that scale.
#'
#' @param values numeric vector of raw signal values
#' @param asinh_threshold NULL, numeric, or character
#'
#' @return list with components:
#'   asinh_vals: numeric vector asinh(values), or NULL if disabled
#'   thr: numeric threshold on asinh scale, or NULL if disabled
#' @keywords internal
parse_asinh_threshold <- function(values, asinh_threshold) {
  if (is.null(asinh_threshold)) {
    return(list(asinh_vals = NULL, thr = NULL))
  }
  
  vals <- as.numeric(values)
  asinh_vals <- asinh(vals)
  
  thr_num <- NULL
  
  if (is.character(asinh_threshold)) {
    thr_str <- asinh_threshold[[1L]]
    
    if (grepl("^q:", thr_str)) {
      # quantile of positive asinh values, e.g. "q:0.75"
      q <- suppressWarnings(as.numeric(sub("^q:", "", thr_str)))
      if (!is.finite(q) || q < 0 || q > 1) {
        stop("Quantile in asinh_threshold must be between 0 and 1, got: ", thr_str)
      }
      
      nz <- asinh_vals[asinh_vals > 0]
      if (!length(nz)) {
        # no positive signal, effectively disable by setting a very low threshold
        thr_num <- -1e6
      } else {
        thr_num <- as.numeric(stats::quantile(
          nz,
          probs = q,
          na.rm = TRUE,
          type  = 7
        ))
      }
      
    } else {
      # plain numeric in character form, interpret as value on original scale
      v <- suppressWarnings(as.numeric(thr_str))
      if (!is.finite(v)) {
        stop("Could not interpret asinh_threshold: ", thr_str)
      }
      if (v < 0) {
        # negative means "no threshold"
        thr_num <- -1e6
      } else {
        thr_num <- asinh(v)
      }
    }
    
  } else {
    # assume already on asinh scale
    v <- as.numeric(asinh_threshold)
    if (!is.finite(v)) {
      stop("asinh_threshold must be numeric or a 'q:0.75' style string.")
    }
    thr_num <- v
  }
  
  list(asinh_vals = asinh_vals, thr = thr_num)
}



# helper: sample block maxima for the null (now uses optimized C++ version)
null_block_maxima <- function(resp,
                              half_window_bins,
                              iters,
                              mask = NULL) {
  resp <- as.numeric(resp)
  n    <- length(resp)
  if (n == 0L) return(numeric(0L))
  
  half_window_bins <- max(as.integer(half_window_bins), 1L)
  iters <- max(1000L, as.integer(iters))
  
  # Use C++ implementation for speed
  if (!is.null(mask)) {
    if (length(mask) != n) {
      warning("`mask` length does not match `resp` length inside null_block_maxima. Ignoring mask.")
      mask <- NULL
    } else {
      mask <- as.logical(mask)
      mask[is.na(mask)] <- FALSE
    }
  }
  
  # Call C++ function 
  if (is.null(mask)) {
    mask <- rep(TRUE, n)
  }
  
  tryCatch({
    cpp_sample_block_maxima(
      response = resp,
      half_window_bins = half_window_bins,
      nsamp = iters,
      mask = mask,
      eps = .Machine$double.eps * 10
    )
  }, error = function(e) {
    # Fallback to original R implementation
    centers <- seq.int(half_window_bins + 1L, n - half_window_bins, by = 1L)
    if (length(mask) == n) {
      centers <- centers[mask[centers]]
    }
    if (!length(centers)) return(numeric(0L))
    
    block_vals <- numeric(iters)
    for (i in seq_len(iters)) {
      c_idx <- sample(centers, size = 1L)
      s <- c_idx - half_window_bins
      e <- c_idx + half_window_bins
      v <- resp[s:e]
      v <- v[is.finite(v)]
      block_vals[i] <- if (length(v)) max(v, na.rm = TRUE) else NA_real_
    }
    block_vals <- block_vals[is.finite(block_vals)]
    if (length(block_vals) >= 10L) {
      low  <- stats::quantile(block_vals, 0.001, na.rm = TRUE, type = 7)
      high <- stats::quantile(block_vals, 0.999, na.rm = TRUE, type = 7)
      block_vals <- block_vals[block_vals > low & block_vals < high]
    }
    block_vals
  })
}





# Template-based peak detection with empirical null distribution and FDR correction
#
# Detects peaks by convolving signal with normalized template, finding local maxima,
# and testing significance against a null distribution built from split-halves sampling.
matchTemplateEmpirical <- function(chromosome,
                                   intervals,
                                   values,
                                   base_template,
                                   initLenBins     = 3L,
                                   asinh_threshold = NULL,
                                   minMatchLengthBP = NULL,
                                   iters           = 20000L,
                                   alpha           = 0.05,
                                   minScore        = NULL,
                                   maxPeaks        = 100000L,
                                   randSeed        = NULL,
                                   smooth_k        = 3L,
                                   null_mask = NULL,
                                   exclude_mask = NULL,
                                   use_split_halves = TRUE,
                                   template_name = NULL,
                                   cascade_level = NULL,
                                   weights = NULL,
                                   verbose = TRUE
                                  ) {
  
  stopifnot(length(intervals) == length(values))
  
  # Helper function for conditional messaging
  vmsg <- function(...) if (isTRUE(verbose)) message(...)
  
  vmsg("[TIMEPOINT 1] Function start - ", Sys.time())
  
  vals <- as.numeric(values)
  ints <- as.numeric(intervals)
  n <- length(vals)
  
  if (n < 5L) {
    return(
      data.frame(
        chromosome = character(0),
        start      = numeric(0),
        end        = numeric(0),
        center_idx = integer(0),
        score      = numeric(0),
        p_value    = numeric(0),
        q_value    = numeric(0),
        neglog10_p = numeric(0),
        neglog10_q = numeric(0),
        span_bins  = integer(0),
        stringsAsFactors = FALSE
      )
    )
  }
  
  # Process exclude mask: regions to exclude from null sampling and candidate detection
  exclude_mask_global <- rep(FALSE, n)
  if (!is.null(exclude_mask)) {
    if (length(exclude_mask) != n) {
      warning("exclude_mask length does not match values length. Ignoring exclude_mask.")
    } else {
      exclude_mask_global <- as.logical(exclude_mask)
      exclude_mask_global[is.na(exclude_mask_global)] <- FALSE
    }
  }
  
  step <- unique(diff(ints))
  if (length(step) != 1L) {
    stop("`intervals` must be evenly spaced.")
  }
  interval_bp <- step[1]
  
  # Optional global seed
  if (!is.null(randSeed)) {
    set.seed(as.integer(randSeed))
  }
  
  # Apply weights if provided: multiply signal values by weights
  if (!is.null(weights)) {
    if (length(weights) != n) {
      warning("`weights` length does not match `values` length. Ignoring weights.")
    } else {
      vals <- vals * as.numeric(weights)
    }
  }
  
  # Normalize template to unit L2 norm (required for proper convolution)
  template <- as.numeric(base_template)
  template_norm <- sqrt(sum(template * template))
  template <- template / max(template_norm, .Machine$double.xmin)
  template_rev <- rev(template)
  
  vmsg("[TIMEPOINT 2] Template prepared - length: ", length(template), " - ", Sys.time())
  
  # Calculate minimum match length: default is template length in base pairs
  if (is.null(minMatchLengthBP) || minMatchLengthBP < 1) {
    minMatchLengthBP <- length(template) * interval_bp
  }
  # Round up to nearest interval_bp boundary
  if (minMatchLengthBP %% interval_bp != 0) {
    minMatchLengthBP <- minMatchLengthBP + interval_bp - (minMatchLengthBP %% interval_bp)
  }
  # Convert to half-window size in bins (used for relative maxima detection)
  rel_window_bins <- max(1L, as.integer(((minMatchLengthBP / interval_bp) / 2) + 1))
  
  vmsg("    Using minMatchLengthBP: ", minMatchLengthBP, " (rel_window_bins: ", rel_window_bins, ")")
  
  # Parse signal threshold: default is 50th quantile of all signal values
  parse_signal_threshold <- function(val) {
    if (is.null(val)) {
      # Default: 50th quantile of all values
      return(quantile(vals, probs = 0.50, na.rm = TRUE, type = 7))
    }
    if (is.character(val)) {
      if (grepl("^q:", val)) {
        q_val <- as.numeric(sub("^q:", "", val))
        if (!is.finite(q_val) || q_val < 0 || q_val > 1) {
          stop("Quantile must be between 0 and 1")
        }
        return(quantile(vals, probs = q_val, na.rm = TRUE, type = 7))
      }
      val <- suppressWarnings(as.numeric(val))
    }
    if (is.numeric(val) && is.finite(val)) {
      return(if (val < 0) -1e6 else val)
    }
    return(quantile(vals, probs = 0.50, na.rm = TRUE, type = 7))
  }
  
  signal_threshold <- parse_signal_threshold(asinh_threshold)
  if (signal_threshold < 0) {
    signal_threshold <- -1e6  # Effectively no threshold
  }
  
  vmsg("    Signal threshold: ", round(signal_threshold, 3))
  
  # Convolve signal with reversed template using optimized C++ implementation
  # Reversing template converts convolution to correlation
  # C++ function handles large arrays efficiently (chunked for >200k bins, FFT for medium)
  vmsg("[TIMEPOINT 3] Starting convolution (optimized C++) - ", Sys.time())
  conv_start <- Sys.time()
  
  resp <- tryCatch({
    # Use optimized C++ convolution (10-50x faster for large arrays)
    # Automatically uses chunked direct convolution for >200k bins
    # Falls back to R's FFT for medium arrays if C++ fails
    cpp_convolve_same(vals, template_rev)
  }, error = function(e) {
    # Fallback to R's FFT-based convolution if C++ not available
    warning("C++ convolution failed, using R FFT fallback. Error: ", e$message, immediate. = TRUE)
    resp_temp <- stats::convolve(vals, template_rev, type = "open")
    half_tmpl <- length(template_rev) %/% 2L
    as.numeric(resp_temp[half_tmpl:(half_tmpl + n - 1L)])
  })
  
  resp <- as.numeric(resp)
  
  conv_time <- as.numeric(difftime(Sys.time(), conv_start, units = "secs"))
  vmsg("    Convolution completed in ", round(conv_time, 2), " seconds")
  
  vmsg("[TIMEPOINT 4] Convolution complete - ", Sys.time())
  
  # Split-halves null sampling: divide chromosome into two halves
  # Test left half using null from right, test right half using null from left
  mid <- floor(n / 2)
  half_left_mask <- seq_len(n) <= mid
  half_right_mask <- !half_left_mask
  
  chrom_start <- as.integer(ints[1])
  chrom_end <- as.integer(ints[n] + interval_bp)
  
  # Collect results from both halves efficiently
  half_results <- list()
  
  # Process split-halves: iterate twice, once for each half
  for (null_half_idx in 1:2) {
    if (null_half_idx == 1) {
      # First pass: test right half using null from left half
      null_mask_use <- half_left_mask
      test_mask_use <- half_right_mask
      tag <- "R"
    } else {
      # Second pass: test left half using null from right half
      null_mask_use <- half_right_mask
      test_mask_use <- half_left_mask
      tag <- "L"
    }
    
    # Apply exclude mask to null sampling region
    null_mask_use <- null_mask_use & !exclude_mask_global
    
    vmsg("[TIMEPOINT 5] Processing ", tag, " half (null from ", if(null_half_idx==1) "left" else "right", ") - ", Sys.time())
    
    # Sample null distribution: random block maxima from null sampling region
    iters_use <- max(1000L, as.integer(iters))
    block_maxima <- null_block_maxima(
      resp = resp,
      half_window_bins = rel_window_bins,
      iters = iters_use,
      mask = null_mask_use
    )
    
    # Fallback: if insufficient samples in half, use pooled mask (all non-excluded regions)
    if (length(block_maxima) < 25L) {
      vmsg("    Warning: Only ", length(block_maxima), " null samples, using pooled mask")
      pooled_mask <- !exclude_mask_global
      block_maxima <- null_block_maxima(
        resp = resp,
        half_window_bins = rel_window_bins,
        iters = iters_use,
        mask = pooled_mask
      )
    }
    
    if (length(block_maxima) < 25L) {
      vmsg("    Warning: Insufficient null samples (", length(block_maxima), "), skipping ", tag, " half")
      next
    }
    
    # Outlier trimming already done in C++ function (0.001-0.999 quantiles)
    # No need to trim again - saves redundant computation
    
    vmsg("    Null distribution: ", length(block_maxima), " samples (median: ", 
         round(median(block_maxima), 4), ")")
    
    # Pre-compute ECDF for efficient p-value calculation
    ecdf_obj <- ecdf(block_maxima)
    
    # Find local maxima in response signal
    cand_idx <- relativeMaxima(resp, orderBins = rel_window_bins, eps = NULL)
    
    # Filter candidates: boundary checks, test region, exclude mask, signal threshold
    candidate_mask <- (
      (cand_idx >= rel_window_bins) &
      (cand_idx < n - rel_window_bins) &
      (test_mask_use[cand_idx]) &
      (!exclude_mask_global[cand_idx]) &
      (vals[cand_idx] >= signal_threshold)
    )
    
    cand_idx <- cand_idx[candidate_mask]
    
    if (length(cand_idx) == 0L) {
      vmsg("    No candidates found in ", tag, " half")
      next
    }
    
    # Cap candidates by maxPeaks: keep top N by signal value
    if (!is.null(maxPeaks) && length(cand_idx) > maxPeaks) {
      cand_signals <- vals[cand_idx]
      top_order <- order(cand_signals, decreasing = TRUE)[1:maxPeaks]
      cand_idx <- cand_idx[top_order]
      vmsg("    Capped candidates to ", maxPeaks, " (top by signal value)")
    }
    
    # Calculate p-values using pre-computed ECDF survival function
    # Survival function: P(X >= x) = 1 - P(X < x) = 1 - ECDF(x)
    resp_values <- resp[cand_idx]
    p_raw <- 1 - ecdf_obj(resp_values)
    p_raw <- pmax(p_raw, .Machine$double.xmin)
    p_raw <- pmin(p_raw, 1.0)
    
    vmsg("    Found ", length(cand_idx), " candidates, ", sum(p_raw < 0.05), " with p < 0.05")
    
    # Calculate initial peak boundaries around candidate centers
    starts_idx <- pmax(cand_idx - rel_window_bins, 1L)
    ends_idx <- pmin(cand_idx + rel_window_bins, n)
    
    # Find point sources: bin with maximum signal within each peak window
    # Use optimized C++ function for speed (much faster than R loop)
    span_bins_use <- rel_window_bins * 2L + 1L
    point_sources_idx <- tryCatch({
      cpp_find_point_sources(
        values = vals,
        center_indices = cand_idx,
        span_bins = rep(span_bins_use, length(cand_idx))
      )
    }, error = function(e) {
      # Fallback to R implementation if C++ fails
      warning("C++ point source finding failed, using R fallback. Error: ", e$message, immediate. = TRUE)
      point_sources_idx <- integer(length(cand_idx))
      for (i in seq_along(cand_idx)) {
        window <- starts_idx[i]:ends_idx[i]
        point_sources_idx[i] <- window[which.max(vals[window])]
      }
      point_sources_idx
    })
    
    # Convert indices to genomic coordinates
    starts <- ints[starts_idx]
    ends <- ints[ends_idx]
    point_sources_abs <- ints[point_sources_idx] + max(1L, interval_bp %/% 2L)
    
    # Recenter peaks at point sources: adjust boundaries symmetrically around point source
    starts <- point_sources_abs - (rel_window_bins * interval_bp)
    ends <- point_sources_abs + (rel_window_bins * interval_bp)
    
    # Clip to chromosome boundaries
    starts <- pmax(as.integer(starts), chrom_start)
    ends <- pmin(as.integer(ends), chrom_end)
    
    # Calculate relative point source position within peak (for output)
    point_sources_rel <- (ints[point_sources_idx] - starts) + max(1L, interval_bp %/% 2L)
    
    # Calculate scores: (1 + response)² normalized to [250, 1000] range
    sq_scores <- (1 + resp[cand_idx])^2
    min_r <- min(sq_scores)
    max_r <- max(sq_scores)
    range_r <- max(max_r - min_r, 1.0)
    scores <- as.integer(250 + 750 * (sq_scores - min_r) / range_r)
    
    # Build data.frame for this half efficiently (much faster than incremental list appending)
    half_df <- data.frame(
      chromosome = rep(chromosome, length(cand_idx)),
      start = starts,
      end = ends,
      center_idx = cand_idx,
      score = scores,
      p_value = p_raw,
      q_value = NA_real_,  # Will be calculated after grouping
      neglog10_p = -log10(pmax(p_raw, 1e-10)),
      neglog10_q = NA_real_,
      span_bins = rep(rel_window_bins * 2L + 1L, length(cand_idx)),
      signal = vals[cand_idx],
      point_source = point_sources_rel,
      template_name = if (!is.null(template_name)) rep(template_name, length(cand_idx)) else NA_character_,
      cascade_level = if (!is.null(cascade_level)) rep(cascade_level, length(cand_idx)) else NA_integer_,
      tag = rep(tag, length(cand_idx)),
      stringsAsFactors = FALSE
    )
    
    # Store this half's results
    half_results[[length(half_results) + 1L]] <- half_df
  }
  
  # Combine results from both halves and apply FDR correction
  if (length(half_results) == 0L) {
    vmsg("    No matches detected")
    return(
      data.frame(
        chromosome = character(0),
        start      = numeric(0),
        end        = numeric(0),
        center_idx = integer(0),
        score      = integer(0),
        p_value    = numeric(0),
        q_value    = numeric(0),
        neglog10_p = numeric(0),
        neglog10_q = numeric(0),
        span_bins  = integer(0),
        point_source = integer(0),
        stringsAsFactors = FALSE
      )
    )
  }
  
  # Combine all halves into single data.frame
  out <- do.call(rbind, half_results)
  
  # Apply FDR correction: group by chromosome, template_name, and cascade_level
  # This ensures FDR is computed separately for each template/cascade combination
  group_cols <- c("chromosome")
  if ("template_name" %in% colnames(out) && !all(is.na(out$template_name))) {
    group_cols <- c(group_cols, "template_name")
  }
  if ("cascade_level" %in% colnames(out) && !all(is.na(out$cascade_level))) {
    group_cols <- c(group_cols, "cascade_level")
  }
  
  out$q_value <- NA_real_
  for (grp_idx in split(seq_len(nrow(out)), interaction(out[group_cols], drop = TRUE))) {
    if (length(grp_idx) > 0L) {
      p_grp <- out$p_value[grp_idx]
      q_grp <- bh_fdr(p_grp)
      out$q_value[grp_idx] <- q_grp
    }
  }
  
  # Filter by alpha threshold
  keep_final <- out$q_value <= alpha
  if (!any(keep_final)) {
    vmsg("    No matches passed FDR filter (alpha=", alpha, ")")
    return(
      data.frame(
        chromosome = character(0),
        start      = numeric(0),
        end        = numeric(0),
        center_idx = integer(0),
        score      = integer(0),
        p_value    = numeric(0),
        q_value    = numeric(0),
        neglog10_p = numeric(0),
        neglog10_q = numeric(0),
        span_bins  = integer(0),
        point_source = integer(0),
        stringsAsFactors = FALSE
      )
    )
  }
  
  out <- out[keep_final, , drop = FALSE]
  
  # Add neglog10 values
  out$neglog10_p <- -log10(pmax(out$p_value, 1e-10))
  out$neglog10_q <- -log10(pmax(out$q_value, 1e-10))
  
  # Sort by chromosome, start, end
  out <- out[order(out$chromosome, out$start, out$end), , drop = FALSE]
  rownames(out) <- NULL
  
  vmsg("[TIMEPOINT 6] Function complete - ", nrow(out), " peaks detected - ", Sys.time())
  
  out
}


#' Recenter peaks at a local point source
#'
#' Given a set of peak calls with center indices and spans in bins,
#' recenter each peak around the bin with maximal signal within its span.
#' This is useful when peaks are detected via template matching but should
#' be centered at the actual signal maximum.
#'
#' @param peaks data.frame with columns: chromosome, start, end, center_idx, ...
#' @param intervals numeric vector of bin starts in base pairs
#' @param values numeric vector of signal values
#' @param add_point_source_bp logical, if TRUE add point_source column
#'
#' @return data.frame with recentered peaks
#' @export
recenterAtPointSource <- function(peaks,
                                  intervals,
                                  values,
                                  add_point_source_bp = TRUE) {
  if (!is.data.frame(peaks)) {
    stop("`peaks` must be a data.frame, e.g. wavematchFromTrack() output")
  }
  
  required_cols <- c("center_idx", "span_bins")
  missing_cols <- setdiff(required_cols, colnames(peaks))
  if (length(missing_cols)) {
    stop(
      "Input `peaks` is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  vals <- as.numeric(values)
  ints <- as.numeric(intervals)
  
  if (length(vals) != length(ints)) {
    stop("`intervals` and `values` must have the same length")
  }
  
  n_track <- length(vals)
  if (n_track == 0L) {
    return(peaks)
  }
  
  # determine bin width from intervals
  step <- unique(diff(ints))
  if (length(step) != 1L) {
    stop("`intervals` must be evenly spaced")
  }
  interval_bp <- step[1]
  if (!is.finite(interval_bp) || interval_bp <= 0) {
    stop("Invalid interval step size")
  }
  
  # ensure basic geometry columns exist in peaks
  if (!"start" %in% names(peaks)) peaks$start <- NA_real_
  if (!"end"   %in% names(peaks)) peaks$end   <- NA_real_
  
  df <- peaks
  n_peaks <- nrow(df)
  
  point_idx <- tryCatch({
    cpp_find_point_sources(
      values = vals,
      center_indices = df$center_idx,
      span_bins = df$span_bins
    )
  }, error = function(e) {
    warning("C++ point source finding failed, using R fallback. Error: ", e$message, immediate. = TRUE)
    # Fallback to R implementation
    point_idx <- integer(n_peaks)
    point_idx[] <- NA_integer_
    
    for (i in seq_len(n_peaks)) {
      L    <- as.integer(df$span_bins[i])
      c_ix <- as.integer(df$center_idx[i])
      
      if (!is.finite(L) || L < 1L || !is.finite(c_ix)) {
        next
      }
      
      half <- L %/% 2L
      s <- max(1L, c_ix - half)
      e <- min(n_track, c_ix + half)
      
      if (s <= e && e <= length(vals)) {
        window <- s:e
        max_idx <- window[which.max(vals[window])]
        point_idx[i] <- max_idx
      }
    }
    point_idx
  })
  
  # update center_idx and optionally add point_source column
  df$center_idx <- point_idx
  
  if (isTRUE(add_point_source_bp)) {
    df$point_source <- ints[point_idx] + max(1L, interval_bp %/% 2L)
  }
  
  # update start/end to be symmetric around point source
  for (i in seq_len(n_peaks)) {
    if (!is.finite(point_idx[i])) next
    
    center_bp <- ints[point_idx[i]] + max(1L, interval_bp %/% 2L)
    L <- as.integer(df$span_bins[i])
    half_span_bp <- (L * interval_bp) %/% 2L
    
    start_bp <- center_bp - half_span_bp
    end_bp <- center_bp + half_span_bp
    
    df$start[i] <- as.numeric(start_bp)
    df$end[i]   <- as.numeric(end_bp)
  }
  
  df
}
