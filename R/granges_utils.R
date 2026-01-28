#' Convert wavematch output to GenomicRanges
#'
#' Convert a data frame of peaks (as returned by \code{wavematchFromTrack}
#' or \code{wavematchFromBam}) into a \code{GRanges} object for use with
#' other Bioconductor packages.
#'
#' @param peaks data.frame with columns \code{chromosome}, \code{start}, \code{end}
#'   and optionally other metadata columns
#' @param keep.mcols logical, if TRUE preserve all columns as metadata columns
#' @param seqinfo optional \code{Seqinfo} object for chromosome information
#'
#' @return A \code{GRanges} object
#' @export
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
asGRanges <- function(peaks, keep.mcols = TRUE, seqinfo = NULL) {
  if (!is.data.frame(peaks) || nrow(peaks) == 0L) {
    return(GenomicRanges::GRanges())
  }
  
  required <- c("chromosome", "start", "end")
  missing <- setdiff(required, colnames(peaks))
  if (length(missing)) {
    stop("peaks must contain columns: ", paste(missing, collapse = ", "))
  }
  
  gr <- GenomicRanges::GRanges(
    seqnames = as.character(peaks$chromosome),
    ranges   = IRanges::IRanges(
      start = as.integer(round(peaks$start)),
      end   = as.integer(round(peaks$end))
    ),
    seqinfo  = seqinfo
  )
  
  if (isTRUE(keep.mcols)) {
    # Add all other columns as metadata
    meta_cols <- setdiff(colnames(peaks), required)
    if (length(meta_cols) > 0L) {
      GenomicRanges::mcols(gr) <- peaks[, meta_cols, drop = FALSE]
    }
  }
  
  gr
}

#' Filter peaks using a blacklist BED file
#'
#' Remove peaks that overlap with regions in a blacklist BED file.
#' This is essential for quality control in ChIP-seq and related assays.
#'
#' @param peaks data.frame or GRanges object with peak coordinates
#' @param blacklist_file path to blacklist BED file, or genome name (e.g., "hg38")
#'   to use the built-in blacklist
#' @param min_overlap minimum overlap fraction (0-1) to consider a peak blacklisted
#'
#' @return Filtered peaks (same type as input)
#' @export
#' @importFrom GenomicRanges GRanges findOverlaps width
#' @importFrom IRanges IRanges
filterBlacklist <- function(peaks, blacklist_file, min_overlap = 0.5) {
  # Convert peaks to GRanges if needed
  if (is.data.frame(peaks)) {
    peaks_gr <- asGRanges(peaks, keep.mcols = TRUE)
    return_df <- TRUE
  } else if (inherits(peaks, "GRanges")) {
    peaks_gr <- peaks
    return_df <- FALSE
  } else {
    stop("peaks must be a data.frame or GRanges object")
  }
  
  # Load blacklist
  if (is.character(blacklist_file) && length(blacklist_file) == 1L) {
    # Check if it's a genome name
    if (file.exists(blacklist_file)) {
      bl_file <- blacklist_file
    } else {
      # Try to get built-in blacklist
      tryCatch({
        bl_file <- get_genome_resource_file(blacklist_file, "blacklist")
      }, error = function(e) {
        stop("Could not find blacklist file: ", blacklist_file)
      })
    }
  } else {
    stop("blacklist_file must be a single file path or genome name")
  }
  
  # Read blacklist BED file
  bl <- utils::read.table(
    bl_file,
    sep = "\t",
    header = FALSE,
    stringsAsFactors = FALSE
  )
  
  if (ncol(bl) < 3L) {
    stop("Blacklist BED file must have at least 3 columns")
  }
  
  bl_gr <- GenomicRanges::GRanges(
    seqnames = bl[[1L]],
    ranges   = IRanges::IRanges(
      start = as.integer(bl[[2L]]) + 1L,  # BED is 0-based, GRanges is 1-based
      end   = as.integer(bl[[3L]])
    )
  )
  
  # Find overlaps
  overlaps <- GenomicRanges::findOverlaps(
    peaks_gr,
    bl_gr,
    minoverlap = as.integer(min_overlap * GenomicRanges::width(peaks_gr))
  )
  
  # Remove overlapping peaks
  keep_idx <- setdiff(seq_along(peaks_gr), S4Vectors::queryHits(overlaps))
  
  filtered <- peaks_gr[keep_idx]
  
  # Convert back to data.frame if needed
  if (return_df) {
    df <- data.frame(
      chromosome = as.character(GenomicRanges::seqnames(filtered)),
      start      = GenomicRanges::start(filtered),
      end        = GenomicRanges::end(filtered),
      stringsAsFactors = FALSE
    )
    
    if (ncol(GenomicRanges::mcols(filtered)) > 0L) {
      df <- cbind(df, as.data.frame(GenomicRanges::mcols(filtered)))
    }
    
    return(df)
  }
  
  filtered
}
