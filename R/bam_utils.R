# BAM utilities for wavematch
#' Check that a BAM file exists and is indexed
#'
#' This helper normalizes the path, checks that the BAM exists,
#' and ensures a BAM index is present. If no index is found,
#' it calls [Rsamtools::indexBam()].
#'
#' @param bam_file Path to a BAM file.
#'
#' @return TRUE (invisibly) if the file exists and is indexed.
#' @export
check_bam_file <- function(bam_file) {
  bam_file <- normalizePath(bam_file, mustWork = TRUE)

  if (!file.exists(bam_file)) {
    stop("Could not find BAM file: ", bam_file)
  }

  bai1 <- paste0(bam_file, ".bai")
  bai2 <- sub("\\.bam$", ".bai", bam_file)

  if (!file.exists(bai1) && !file.exists(bai2)) {
    message(
      "No index for ", bam_file,
      "  creating with Rsamtools::indexBam()"
    )
    Rsamtools::indexBam(bam_file)
  }

  invisible(TRUE)
}

#' Estimate read length from a BAM file
#'
#' Uses [Rsamtools::scanBam()] to look at the first `num_reads`
#' alignments and returns the median query width.
#'
#' @param bam_file Path to a BAM file.
#' @param num_reads Maximum number of reads to inspect.
#'
#' @return Integer scalar, estimated read length.
#' @export
estimate_read_length <- function(bam_file, num_reads = 10000L) {
  bam_file <- normalizePath(bam_file, mustWork = TRUE)
  check_bam_file(bam_file)
  
  param <- Rsamtools::ScanBamParam(what = "qwidth")
  aln <- GenomicAlignments::readGAlignments(
    bam_file,
    param = param,
    use.names = FALSE
  )
  widths <- GenomicRanges::width(aln)
  if (length(widths) == 0L) {
    stop("No reads found in BAM file: ", bam_file)
  }
  widths <- head(widths, num_reads)
  tab <- table(widths)
  as.integer(names(sort(tab, decreasing = TRUE)[1L]))
}
  

#' Estimate read lengths for multiple BAM files
#'
#' Vector interface to [estimate_read_length()].
#'
#' @param bam_files Character vector of BAM file paths.
#' @param num_reads Maximum number of reads to inspect in each file.
#'
#' @return Named integer vector of read lengths.
#' @export
get_read_lengths <- function(bam_files, num_reads = 10000L) {
  if (length(bam_files) == 0L) {
    stop("bam_files must have length greater than zero")
  }

  bam_files <- normalizePath(bam_files, mustWork = TRUE)

  vals <- vapply(
    bam_files,
    estimate_read_length,
    integer(1L),
    num_reads = num_reads
  )
  stats::setNames(vals, bam_files)
}

#' Check whether BAM files are paired end
#'
#' Lightweight wrapper around the C++ helper [cpp_inferPairedEnd()],
#' which inspects up to `max_reads` alignments per BAM and returns
#' TRUE if any alignment is part of a proper pair.
#'
#' @param bam_files Character vector of BAM file paths.
#' @param max_reads Integer, maximum number of records to inspect
#'   in each BAM (default 10000).
#'
#' @return Logical named vector, TRUE for paired end BAMs, FALSE otherwise.
#' @export
bams_are_paired_end <- function(bam_files, max_reads = 10000L) {
  if (length(bam_files) == 0L) {
    stop("bam_files must have length greater than zero")
  }

  bam_files <- normalizePath(bam_files, mustWork = TRUE)
  max_reads <- as.integer(max_reads)

  vapply(
    bam_files,
    function(path) {
      check_bam_file(path)
      cpp_inferPairedEnd(path, maxReads = max_reads)
    },
    logical(1L),
    USE.NAMES = TRUE
  )
}
