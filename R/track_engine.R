#' wavematch track representation
#'
#' Simple list based representation of a fixed step genomic track.
#'
#' @param chromosome chromosome name
#' @param intervals integer vector of start positions for each bin
#' @param values numeric vector of values for each bin
#' @param step_bp bin size in base pairs
#' @param kind character label, for example "coverage"
#'
#' @return list with class "SPTrack"
#' @export
Track <- function(chromosome, intervals, values, step_bp, kind = "coverage") {
  if (length(intervals) != length(values)) {
    stop("intervals and values must have the same length")
  }
  structure(
    list(
      chromosome = as.character(chromosome),
      intervals  = as.integer(intervals),
      values     = as.numeric(values),
      step_bp    = as.integer(step_bp),
      kind       = as.character(kind)
    ),
    class = "SPTrack"
  )
}

#' Bin coverage for a single BAM on a chromosome
#'
#' @param bam_file path to BAM file
#' @param chromosome chromosome name
#' @param chrom_length length of the chromosome in base pairs
#' @param step_bp bin size in base pairs
#'
#' @return [Track] object
#' @keywords internal
bin_coverage_one_bam <- function(
  bam_file,
  chromosome,
  chrom_length,
  step_bp
) {
  check_bam_file(bam_file)

  region <- GenomicRanges::GRanges(
    seqnames = chromosome,
    ranges   = IRanges::IRanges(1L, chrom_length)
  )

  param <- Rsamtools::ScanBamParam(which = region)
  aln <- GenomicAlignments::readGAlignments(bam_file, param = param)
  cov <- GenomicRanges::coverage(aln)[[chromosome]]

  breaks <- seq.int(1L, chrom_length, by = step_bp)
  if (tail(breaks, 1L) < chrom_length) {
    breaks <- c(breaks, chrom_length + 1L)
  }

  means <- numeric(length(breaks) - 1L)
  for (i in seq_len(length(means))) {
    s <- breaks[i]
    e <- breaks[i + 1L] - 1L
    means[i] <- mean(cov[s:e])
  }

  Track(
    chromosome = chromosome,
    intervals  = breaks[-length(breaks)],
    values     = means,
    step_bp    = step_bp,
    kind       = "coverage"
  )
}

#' Build coverage tracks for multiple BAM files
#'
#' Given a set of BAMs and a genome definition, this constructs
#' fixed step coverage tracks for each chromosome and each BAM.
#' At this stage this is a simple engine that can be fed into
#' wavelet matching functions.
#'
#' @param bam_files character vector of BAM paths
#' @param genome genome string, see [resolve_genome_name]
#' @param step_bp bin size in base pairs
#' @param chrom_sizes_file optional explicit path to sizes file
#' @param chromosomes optional character vector of chromosomes to include
#' @param exclude_chroms optional vector of chromosomes to drop
#' @param max_chroms optional limit on number of chromosomes processed
#'
#' @return named list, one entry per chromosome, each is a list of [Track] objects
#' @export
build_tracks_from_bams <- function(
  bam_files,
  genome,
  step_bp = 25L,
  chrom_sizes_file = NULL,
  chromosomes = NULL,
  exclude_chroms = character(),
  max_chroms = Inf
) {
  if (is.null(chrom_sizes_file)) {
    sizes_file <- get_genome_resource_file(genome, "sizes")
  } else {
    sizes_file <- chrom_sizes_file
  }

  chrom_sizes <- get_chrom_sizes_dict(
    sizes_file,
    exclude_chroms = exclude_chroms
  )

  if (is.null(chromosomes)) {
    chroms <- names(chrom_sizes)
  } else {
    chroms <- intersect(chromosomes, names(chrom_sizes))
  }

  if (!is.finite(max_chroms)) {
    max_chroms <- length(chroms)
  }
  chroms <- head(chroms, max_chroms)

  tracks <- list()

  for (chrom in chroms) {
    message("Building coverage tracks for ", chrom)
    chrom_len <- chrom_sizes[[chrom]]

    bam_tracks <- lapply(
      bam_files,
      function(bam) bin_coverage_one_bam(
        bam_file     = bam,
        chromosome   = chrom,
        chrom_length = chrom_len,
        step_bp      = step_bp
      )
    )

    tracks[[chrom]] <- bam_tracks
  }

  tracks
}
