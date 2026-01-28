#' Automatic minimum run length candidates
#'
#' Infer candidate minimum run lengths (in bins) from a signal.
#'
#' 1) asinh transform of the signal
#' 2) subtract a running median to get residuals
#' 3) keep positive residuals and threshold at a high quantile
#' 4) group contiguous runs that pass the threshold
#' 5) drop runs shorter than \code{initLen}
#' 6) optionally weight runs by area (width times mean residual)
#' 7) choose representative run lengths at given cumulative weight quantiles
#'
#' @param values numeric vector of signal values at uniform genomic bins
#' @param initLen integer scalar, lower bound on minimum run length in bins
#' @param quantiles numeric vector of probabilities in (0,1) used to pick
#'   representative run lengths from the weighted distribution of widths
#' @param weight_by_intensity logical, if true weight each run by
#'   \code{width times mean residual}; if false weight each run equally
#'
#' @return integer vector of candidate run lengths (bins), sorted and unique.
#'   Always returns at least \code{initLen}.
#' @export
autoMinLengthCandidates <- function(values,
                                    initLen = 3L,
                                    quantiles = c(0.50, 0.75, 0.90),
                                    weight_by_intensity = TRUE) {
  # coerce and guard
  if (length(values) == 0L) {
    return(as.integer(initLen))
  }

  x <- as.numeric(values)
  x[!is.finite(x)] <- 0.0

  initLen <- as.integer(initLen)
  if (initLen < 1L) {
    initLen <- 1L
  }

  n <- length(x)
  if (n < 1L) {
    return(as.integer(initLen))
  }

  # asinh transform
  tr <- asinh(x)

  # choose median filter kernel size
  ks_target <- max(2L * initLen + 1L,
                   2L * as.integer(n * 0.005) + 1L)
  # force odd
  ks_odd <- bitwOr(ks_target, 1L)
  # maximum odd kernel size that does not exceed n
  ks_max <- if (n %% 2L == 1L) n else (n - 1L)

  if (n >= 3L) {
    ksize <- min(ks_odd, ks_max)
  } else {
    ksize <- 3L
  }

  # running median
  if (n >= ksize && ksize >= 3L) {
    med <- stats::runmed(tr, k = ksize, endrule = "median")
  } else {
    med <- rep(stats::median(tr), n)
  }

  trValues <- tr - med

  # keep only positive residuals
  nz <- trValues[trValues > 0]
  if (length(nz) == 0L) {
    return(as.integer(initLen))
  }

  # high quantile threshold on positive residuals
  thr <- stats::quantile(
    nz,
    probs = 0.90,
    type = 8,
    names = FALSE,
    na.rm = TRUE
  )

  mask <- trValues >= thr
  if (!any(mask)) {
    return(as.integer(initLen))
  }

  # contiguous runs where mask is TRUE
  ext <- c(FALSE, mask, FALSE)
  d   <- diff(ext)
  idx <- which(d != 0L)

  if (length(idx) < 2L) {
    return(as.integer(initLen))
  }

  # pair up starts and ends (end is exclusive)
  starts <- idx[seq.int(1L, length(idx) - 1L, by = 2L)]
  ends_excl <- idx[seq.int(2L, length(idx), by = 2L)]

  widths <- ends_excl - starts
  if (length(widths) == 0L) {
    return(as.integer(initLen))
  }

  # drop very short runs
  keep <- widths >= initLen
  if (!any(keep)) {
    return(as.integer(initLen))
  }

  starts    <- starts[keep]
  ends_excl <- ends_excl[keep]
  widths    <- widths[keep]

  # weights per run
  if (weight_by_intensity) {
    means <- vapply(
      seq_along(starts),
      function(i) {
        s <- starts[i]
        e <- ends_excl[i] - 1L
        v <- trValues[s:e]
        mean(v, na.rm = TRUE)
      },
      numeric(1L)
    )
    wts <- widths * means
  } else {
    wts <- rep(1.0, length(widths))
  }

  # if all weights nonpositive, fall back to median width
  if (!any(is.finite(wts)) || max(wts, na.rm = TRUE) <= 0) {
    med_w <- stats::median(widths)
    cand  <- as.integer(round(med_w))
    cand  <- max(cand, initLen)
    return(as.integer(cand))
  }

  # sort widths and accumulate weights
  ord <- order(widths)
  w   <- wts[ord]
  w[!is.finite(w)] <- 0.0

  cw <- cumsum(w)
  tot <- cw[length(cw)]
  if (tot <= 0) {
    med_w <- stats::median(widths)
    cand  <- as.integer(round(med_w))
    cand  <- max(cand, initLen)
    return(as.integer(cand))
  }

  cw <- cw / tot

  out <- integer(0L)
  widths_sorted <- widths[ord]

  for (q in quantiles) {
    qf <- as.numeric(q)
    if (!is.finite(qf)) {
      next
    }
    qf <- max(0.0, min(1.0, qf))

    j <- which(cw >= qf)[1L]
    if (is.na(j)) {
      j <- length(widths_sorted)
    }

    out <- c(out, as.integer(widths_sorted[j]))
  }

  # enforce minimum length, unique, sorted ascending
  out <- unique(out[out >= initLen])
  out <- sort(out)
  # add a small baseline ladder of multiples of initLen
  ladder <- as.integer(initLen * c(1L, 2L, 3L))
  out <- sort(unique(c(out, ladder)))
  out <- out[out >= initLen]
  
  if (length(out) == 0L) {
    return(as.integer(initLen))
  }

  as.integer(out)
}
