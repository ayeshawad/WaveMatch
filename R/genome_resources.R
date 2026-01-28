#' Effective genome sizes and resource helpers
#'
#' Internal constants and helpers to resolve genome names and locate
#' sizes, blacklist, and sparse resource files shipped with the package.
#'
#' @keywords internal
#' @name genome_resources


.EFFECTIVE_GENOME_SIZES <- list(
  hg19 = c(
    "50"  = 2685511454,
    "75"  = 2736124898,
    "100" = 2776919708,
    "150" = 2827436883,
    "200" = 2855463800,
    "250" = 2855044784
  ),
  hg38 = c(
    "50"  = 2701495711,
    "75"  = 2747877702,
    "100" = 2805636231,
    "150" = 2862010428,
    "200" = 2887553103,
    "250" = 2898802627
  ),
  t2t = c(
    "50"  = 2725240337,
    "75"  = 2786136059,
    "100" = 2814334875,
    "150" = 2931551487,
    "200" = 2936403235,
    "250" = 2960856300
  ),
  mm10 = c(
    "50"  = 2308125299,
    "75"  = 2407883243,
    "100" = 2467481008,
    "150" = 2494787038,
    "200" = 2520868989,
    "250" = 2538590322
  ),
  mm39 = c(
    "50"  = 2309746861,
    "75"  = 2410055689,
    "100" = 2468088461,
    "150" = 2495461690,
    "200" = 2521902382,
    "250" = 2538633971
  ),
  dm3 = c(
    "50"  = 130428510,
    "75"  = 135004387,
    "100" = 139647132,
    "150" = 144307658,
    "200" = 148523810,
    "250" = 151901455
  ),
  dm6 = c(
    "50"  = 125464678,
    "75"  = 127324557,
    "100" = 129789773,
    "150" = 129940985,
    "200" = 132508963,
    "250" = 132900923
  ),
  ce11 = c(
    "50"  = 95159402,
    "75"  = 96945370,
    "100" = 98259898,
    "150" = 98721103,
    "200" = 98672558,
    "250" = 101271756
  )
)

#' Resolve a genome name into a canonical label
#'
#' @param genome character, user supplied genome string
#'
#' @return canonical genome name such as "hg38"
#' @export
resolve_genome_name <- function(genome) {
  g <- tolower(genome)
  if (g %in% c("hg19", "grch37")) return("hg19")
  if (g %in% c("hg38", "grch38")) return("hg38")
  if (g %in% c("t2t", "chm13", "t2t-chm13")) return("t2t")
  if (g %in% c("mm10", "grcm38")) return("mm10")
  if (g %in% c("mm39", "grcm39")) return("mm39")
  if (g %in% c("dm3")) return("dm3")
  if (g %in% c("dm6")) return("dm6")
  if (g %in% c("ce10", "ws220")) return("ce10")
  if (g %in% c("ce11", "wbcel235")) return("ce11")
  stop("Genome ", genome, " is not recognized")
}

#' Get effective genome size for a genome and read length
#'
#' @param genome character, see [resolve_genome_name]
#' @param read_length numeric, read length in base pairs
#'
#' @return integer effective genome size in base pairs
#' @export
get_effective_genome_size <- function(genome, read_length) {
  g <- resolve_genome_name(genome)
  sizes <- .EFFECTIVE_GENOME_SIZES[[g]]
  if (is.null(sizes)) {
    stop("Defaults not available for ", g)
  }
  rl <- as.numeric(names(sizes))
  nearest <- rl[which.min(abs(rl - read_length))]
  unname(sizes[as.character(nearest)])
}

#' Locate a genome resource file shipped with the package
#'
#' @param genome character, see [resolve_genome_name]
#' @param type character, one of "sizes", "blacklist", "sparse"
#'
#' @return absolute file path
#' @export
get_genome_resource_file <- function(genome, type) {
  g <- resolve_genome_name(genome)
  file_name <- switch(
    tolower(type),
    "sizes"     = paste0(g, ".sizes"),
    "blacklist" = paste0(g, "_blacklist.bed"),
    "sparse"    = paste0(g, "_sparse.bed"),
    stop("Unknown file type: ", type)
  )
  path <- system.file("extdata", file_name, package = "wavematch")
  if (!nzchar(path)) {
    stop("Resource file ", file_name, " not found in extdata")
  }
  path
}

#' Read chromosome sizes file into a named vector
#'
#' This mirrors the misc_util getChromSizesDict helper from the Python code.
#'
#' @param sizes_file path to a sizes file
#' @param exclude_regex regular expression for allowed chromosomes
#' @param exclude_chroms character vector of chromosomes to drop
#'
#' @return named integer vector of chromosome sizes
#' @export
get_chrom_sizes_dict <- function(
  sizes_file,
  exclude_regex = "^chr[A-Za-z0-9]+$",
  exclude_chroms = character()
) {
  df <- utils::read.table(
    sizes_file,
    sep = "\t",
    header = FALSE,
    col.names = c("chrom", "size"),
    stringsAsFactors = FALSE
  )
  keep <- grepl(exclude_regex, df$chrom) & !(df$chrom %in% exclude_chroms)
  out <- df$size[keep]
  names(out) <- df$chrom[keep]
  out
}

#' Load genome assembly resources (chromSizes)
#'
#' Converts genome string (e.g. "hg38") into an assembly list with $chromSizes
#' using the installed extdata *.sizes file.
#'
#' @param genome character, genome name such as "hg38"
#'
#' @return list with elements: name, chromSizes
#' @export
load_genome_resources <- function(genome) {
  g <- resolve_genome_name(genome)
  sizes_file <- get_genome_resource_file(g, "sizes")
  chromSizes <- get_chrom_sizes_dict(sizes_file)
  list(
    name = g,
    chromSizes = chromSizes
  )
}
