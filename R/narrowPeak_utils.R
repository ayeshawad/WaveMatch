# internal helper to write a narrowPeak like file
#
# Expects a data.frame with at least:
#   chromosome, start, end, score, neglog10_p, neglog10_q
#
# We map to the ten narrowPeak columns as:
#   1 chrom        -> chromosome
#   2 chromStart   -> start
#   3 chromEnd     -> end
#   4 name         -> "peak_<index>"
#   5 score        -> integer in [0,1000] derived from neglog10_q
#   6 strand       -> "."
#   7 signalValue  -> score (template response)
#   8 pValue       -> neglog10_p
#   9 qValue       -> neglog10_q
#  10 peak         -> -1 (no summit)
#
.sp_write_narrowPeak <- function(peaks, file) {
  if (!is.data.frame(peaks) || nrow(peaks) == 0L) {
    return(invisible(NULL))
  }

  required <- c("chromosome", "start", "end", "score",
                "neglog10_p", "neglog10_q")
  miss <- setdiff(required, names(peaks))
  if (length(miss)) {
    stop("peaks table is missing required columns: ",
         paste(miss, collapse = ", "))
  }

  chrom  <- as.character(peaks$chromosome)
  start  <- as.integer(round(peaks$start))
  end    <- as.integer(round(peaks$end))
  name   <- sprintf("peak_%d", seq_len(nrow(peaks)))

  # Use normalized score if available, otherwise derive from neglog10_q
  if ("score" %in% names(peaks) && all(is.finite(peaks$score))) {
    score <- as.integer(pmax(0, pmin(1000, round(peaks$score))))
  } else {
    # Fallback: derive from neglog10_q
    score_raw <- peaks$neglog10_q
    score_raw[!is.finite(score_raw)] <- 0
    score <- as.integer(pmax(0, pmin(1000, round(10 * score_raw))))
  }

  strand      <- rep(".", length(chrom))
  # Use score_raw if available, otherwise use normalized score
  if ("score_raw" %in% names(peaks)) {
    signalValue <- as.numeric(peaks$score_raw)
  } else {
    signalValue <- as.numeric(peaks$score)
  }
  pValue      <- as.numeric(peaks$neglog10_p)
  qValue      <- as.numeric(peaks$neglog10_q)
  # Use point_source_bp if available for peak summit
  if ("point_source_bp" %in% names(peaks)) {
    peak <- as.integer(round(peaks$point_source_bp - peaks$start))
  } else {
    peak <- rep(-1L, length(chrom))
  }

  np <- data.frame(
    chrom, start, end, name, score, strand,
    signalValue, pValue, qValue, peak,
    stringsAsFactors = FALSE
  )

  utils::write.table(
    np,
    file      = file,
    sep       = "\t",
    quote     = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )

  invisible(NULL)
}
