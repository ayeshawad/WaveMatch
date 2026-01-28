#' Resample a base template to a target length in bins
#'
#' Given a base 1D template, resample it to a target number of bins
#' and L2 normalize the result.
#'
#' @param base_template numeric vector representing a unit norm template
#' @param target_bins desired length in bins (will be coerced to integer,
#'   and clamped to at least 3)
#'
#' @return numeric vector of length \code{target_bins} with L2 norm 1,
#'   or all zeros if normalization fails.
#' @export
makeTemplateAtBins <- function(base_template, target_bins) {
  x <- as.numeric(base_template)

  if (length(x) == 0L) {
    return(numeric(0L))
  }

  tb <- as.integer(max(3L, target_bins))

  if (tb == length(x)) {
    t <- x
  } else {
    # simple linear resampling on [0, 1]
    old_idx <- seq(0, 1, length.out = length(x))
    new_idx <- seq(0, 1, length.out = tb)
    t <- stats::approx(old_idx, x, xout = new_idx, rule = 2)$y
  }

  nrm <- sqrt(sum(t * t))

  if (!is.finite(nrm) || nrm <= 0) {
    return(rep(0, tb))
  }

  t / nrm
}

#' Template response at a given center index with zero padding
#'
#' Compute the dot product between a reversed template and a window
#' of the signal centered at a given index, using zero padding at edges.
#'
#' @param values numeric vector of signal values
#' @param center_idx integer index in \code{values} (1 based)
#' @param tmpl_rev numeric vector, template already reversed
#'
#' @return numeric scalar, the dot product response
#' @export
responseAt <- function(values, center_idx, tmpl_rev) {
  v <- as.numeric(values)
  center_idx <- as.integer(center_idx)
  L <- length(tmpl_rev)

  if (L == 0L || length(v) == 0L) {
    return(0.0)
  }

  half <- L %/% 2L

  # outside bounds
  idx <- center_idx + (-half:half)

  inside <- idx >= 1L & idx <= length(v)
  w <- numeric(L)
  w[inside] <- v[idx[inside]]

  sum(w * tmpl_rev)
}
