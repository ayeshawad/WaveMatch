#' Auto-tuning functions for wavematch parameters
#' 
#' These functions estimate optimal parameter values from signal characteristics
#' without requiring expensive peak calling operations.

#' Auto-tune cascade level based on peak width distribution
#' 
#' Uses autoMinLengthCandidates to estimate dominant peak widths, then converts
#' to appropriate cascade level for the given wavelet template.
#' 
#' Memory-efficient: O(n log n) operations (running median), no peak calling required.
#'
#' @param values numeric vector of signal values
#' @param interval_bp integer, bin size in base pairs, default 25L
#' @param template_name character, wavelet template name ("haar", "db2", "db3")
#' @param target_width_bp numeric, target peak width in bp (optional override)
#' @return integer, optimal cascade level
#' @export
autoTuneCascadeLevel <- function(values,
                                 interval_bp = 25L,
                                 template_name = "db3",
                                 target_width_bp = NULL) {
  vals <- as.numeric(values)
  vals <- vals[is.finite(vals)]
  
  if (length(vals) < 100L) {
    # Default cascade levels by template
    defaults <- list("haar" = 3L, "db2" = 2L, "db3" = 2L)
    default_val <- defaults[[template_name]]
    if (is.null(default_val)) default_val <- 2L
    return(default_val)
  }
  
  # Estimate peak widths using existing autoMinLengthCandidates function
  # This analyzes signal to find dominant run lengths (peak widths)
  width_candidates_bins <- tryCatch({
    autoMinLengthCandidates(
      values = vals,
      initLen = 3L,
      quantiles = c(0.50, 0.75, 0.90),
      weight_by_intensity = TRUE
    )
  }, error = function(e) {
    # Fallback: use simple heuristic
    as.integer(round(median(vals > median(vals, na.rm = TRUE), na.rm = TRUE) * 10L))
  })
  
  if (length(width_candidates_bins) == 0L || all(is.na(width_candidates_bins))) {
    # Fallback to defaults
    defaults <- list("haar" = 3L, "db2" = 2L, "db3" = 2L)
    return(defaults[[template_name]] %||% 2L)
  }
  
  # Convert bins to base pairs
  width_candidates_bp <- width_candidates_bins * as.integer(interval_bp)
  
  # Use median width as estimate (most representative)
  median_width_bp <- median(width_candidates_bp, na.rm = TRUE)
  
  # If target width provided, use it (allows user override)
  if (!is.null(target_width_bp) && is.finite(target_width_bp) && target_width_bp > 0) {
    median_width_bp <- target_width_bp
  }
  
  # Estimate cascade level needed to match this width
  # Template length ≈ 2^(cascade_level) * base_filter_length
  
  base_lengths <- list(
    "haar" = 2L,
    "db2" = 4L,
    "db3" = 6L
  )
  
  base_len <- base_lengths[[template_name]]
  if (is.null(base_len) || !is.finite(base_len)) {
    base_len <- 6L  # Default to db3
  }
  
  # Estimate: template_length_bins ≈ median_width_bp / interval_bp
  target_length_bins <- max(3L, as.integer(round(median_width_bp / interval_bp)))
  
  # Find cascade level: 2^level * base_len ≈ target_length
  # level ≈ log2(target_length / base_len)
  if (target_length_bins <= base_len) {
    cascade_level <- 1L
  } else {
    cascade_level <- as.integer(ceiling(log2(target_length_bins / base_len)))
  }
  
  # Apply caps to prevent over-smoothing
  if (template_name == "db2") {
    cascade_level <- min(cascade_level, 3L)
  } else if (template_name == "db3") {
    cascade_level <- min(cascade_level, 2L)
  }
  
  # Ensure minimum
  cascade_level <- max(1L, cascade_level)
  
  cascade_level
}

