#' Alignment Class
#'
#' @description
#' Simple wrapper around an alignment file (e.g., BAM)
#' for generating read coverage tracks, etc.
#' @param filePath path to alignment file (BAM)
#' @param assembly Assembly object corresponding to the alignment
#' @param pairedEnd logical, whether the alignment is paired-end
#' @param fragmentLength integer, fragment length (in bp)
#' @param scaleFactor numeric, scaling factor for normalization
#' @return An Alignment object
#' @export
Alignment <- R6::R6Class(
  "Alignment",
  public = list(
    filePath       = NULL,
    assembly       = NULL,
    pairedEnd      = NULL,
    fragmentLength = NULL,
    scaleFactor    = NULL,
    
    initialize = function(filePath,
                          assembly,
                          pairedEnd      = NULL,
                          fragmentLength = NULL,
                          scaleFactor    = NULL) {
      stopifnot(file.exists(filePath))
      self$filePath       <- normalizePath(filePath)
      self$assembly       <- assembly
      self$pairedEnd      <- pairedEnd
      self$fragmentLength <- fragmentLength
      self$scaleFactor    <- scaleFactor
      
      if (is.null(self$pairedEnd)) {
        self$pairedEnd <- self$inferPairedEnd()
      }
      if (is.null(self$fragmentLength)) {
        self$fragmentLength <- self$inferFragmentLength()
      }
      if (is.null(self$scaleFactor)) {
        self$scaleFactor <- self$inferScaleFactor(method = "coverage")
      }
      
      invisible(self)
    },
    
    #' @description
    #' Return a raw (unscaled) alignment count track for a chromosome
    #' @param chrom chromosome name, e.g. "chr1"
    #' @param binSize size of genomic intervals over which to measure aligned reads
    #' @return numeric vector of counts per bin
    getRawTrack = function(chrom, binSize = 25L) {
      chromSizes <- self$assembly$chromSizes
      if (is.null(chromSizes) || is.null(names(chromSizes))) {
        stop("assembly$chromSizes must contain a named numeric vector chromSizes")
      }
      if (!chrom %in% names(chromSizes)) {
        stop("Chromosome ", chrom, " not found in assembly$chromSizes")
      }
      
      chromLen <- as.integer(chromSizes[[chrom]])
      step_bp  <- as.integer(binSize)
      
      # call C++ binner
      vals <- cpp_binned_coverage_bam(
        bam_path     = self$filePath,
        chrom        = chrom,
        chrom_length = chromLen,
        step_bp      = step_bp
      )
      
      vals <- as.numeric(vals)
      
      # optional scaling (default 1.0)
      if (!is.null(self$scaleFactor)) {
        vals <- vals * self$scaleFactor
      }
      
      return(vals)
    },
    
    #' @description
    #' Infer whether the file is paired-end by scanning a subset of reads.
    inferPairedEnd = function(maxReads = 10000L) {
      cpp_inferPairedEnd(self$filePath, maxReads = as.integer(maxReads))
    },
    
    #' @description
    #' Infer fragment length (paired or single-end).
    inferFragmentLength = function() {
      return(0L)
    },
    
    #' @description
    #' Determine scaling factor for normalization between replicates.
    inferScaleFactor = function(method = c("coverage", "RPKM", "DESeq")) {
      method <- match.arg(method)
      return(1.0)
    }
  ),
  
  private = list(
    getChromLength = function(chrom) {
      cs <- self$assembly$chromSizes
      if (is.null(names(cs))) {
        stop("assembly$chromSizes must be a named vector.")
      }
      if (!chrom %in% names(cs)) {
        stop("Chromosome ", chrom, " not found in assembly$chromSizes.")
      }
      return(cs[[chrom]])
    }
  )
)
