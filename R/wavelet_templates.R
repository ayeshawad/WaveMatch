# Internal wavelet filter catalog
#
# Each entry returns a list with components:
#   h: scaling filter coefficients
#   g: wavelet filter coefficients
#
# These are standard Daubechies filter coefficients.

.sp_get_wavelet_filters <- function(name) {
  name <- tolower(as.character(name))
  
  if (name == "haar" || name == "db1") {
    # Haar scaling and wavelet filters
    h <- c(1, 1) / sqrt(2)
    g <- c(1, -1) / sqrt(2)
    return(list(h = h, g = g))
  }
  
  if (name == "db2") {
    # Daubechies 2 (length 4)
    s3 <- sqrt(3)
    denom <- 4 * sqrt(2)
    h <- c(
      (1 + s3) / denom,
      (3 + s3) / denom,
      (3 - s3) / denom,
      (1 - s3) / denom
    )
    # standard quadrature mirror
    g <- c(
      (1 - s3) / denom,
      -(3 - s3) / denom,
      (3 + s3) / denom,
      -(1 + s3) / denom
    )
    return(list(h = h, g = g))
  }
  
  if (name == "db3") {
    # Daubechies 3 (length 6)
    # numeric approximations to keep it readable
    h <- c(
      0.3326705529500826,
      0.8068915093110928,
      0.4598775021184915,
      -0.1350110200102546,
      -0.08544127388202666,
      0.03522629188570953
    )
    g <- c(
      0.03522629188570953,
      0.08544127388202666,
      -0.1350110200102546,
      -0.4598775021184915,
      0.8068915093110928,
      -0.3326705529500826
    )
    return(list(h = h, g = g))
  }
  
  stop("Unknown wavelet '", name,
       "'. Supported names currently include 'haar', 'db1', 'db2', 'db3'.")
}

# One step of cascade refinement
.sp_refine_sequence <- function(seq, filt) {
  seq  <- as.numeric(seq)
  filt <- as.numeric(filt)
  if (!length(seq) || !length(filt)) {
    return(numeric(0L))
  }

  # upsample by 2 with zeros in between
  up <- numeric(length(seq) * 2L)
  up[seq(1L, by = 2L, length.out = length(seq))] <- seq

  # ensure time series is at least as long as the filter
  if (length(up) < length(filt)) {
    up <- c(up, numeric(length(filt) - length(up)))
  }

  out <- stats::filter(up, filt, sides = 1L, circular = FALSE)
  out[is.na(out)] <- 0
  as.numeric(out)
}

# Build a template from wavelet filters by iterative refinement
#
# name: wavelet family, e.g. "haar", "db2"
# cascade_level: integer >= 1
# use_scaling: if TRUE, build scaling function like template,
#              else build wavelet like template.
# grid_len: if not NULL, resample final template to this length
#
# Returns a unit norm numeric vector.
build_wavelet_template <- function(name,
                                   cascade_level = 3L,
                                   use_scaling = TRUE,
                                   grid_len = NULL) {
  wf <- .sp_get_wavelet_filters(name)
  h <- wf$h
  g <- wf$g
  
  L <- as.integer(cascade_level)
  if (!is.finite(L) || L < 1L) {
    L <- 1L
  }
  
  # start from a single impulse
  if (isTRUE(use_scaling)) {
    seq <- 1
    for (lvl in seq_len(L)) {
      seq <- .sp_refine_sequence(seq, h)
    }
  } else {
    seq <- 1
    for (lvl in seq_len(L)) {
      seq <- .sp_refine_sequence(seq, g)
    }
  }
  
  seq <- as.numeric(seq)
  seq[!is.finite(seq)] <- 0
  
  # optional resampling to a fixed grid length for stability between names
  if (!is.null(grid_len)) {
    grid_len <- as.integer(grid_len)
    if (grid_len > 2L && length(seq) > 1L) {
      old_idx <- seq(0, 1, length.out = length(seq))
      new_idx <- seq(0, 1, length.out = grid_len)
      seq <- stats::approx(old_idx, seq, xout = new_idx, rule = 2)$y
      seq[!is.finite(seq)] <- 0
    }
  }
  
  nrm <- sqrt(sum(seq * seq))
  if (!is.finite(nrm) || nrm <= 0) {
    return(rep(0, length(seq)))
  }
  seq / nrm
}

# Optional registry for caching templates by spec
.sp_template_cache <- new.env(parent = emptyenv())

.sp_template_key <- function(name, cascade_level, use_scaling, grid_len) {
  sprintf("%s|L=%d|S=%d|G=%s",
          tolower(as.character(name)),
          as.integer(cascade_level),
          if (isTRUE(use_scaling)) 1L else 0L,
          if (is.null(grid_len)) "none" else as.integer(grid_len))
}

# Resolve a template specification to a numeric base template
#
# If base_template is provided, this is returned (optionally normalized).
# If template_name is provided, we build or fetch the corresponding wavelet template.
#
# Automatically caps db2/db3 template lengths to prevent over-smoothing
# and false positives. Uses cascade level to control natural template length.

resolve_template_spec <- function(base_template = NULL,
                                  template_name = NULL,
                                  cascade_level = 3L,
                                  use_scaling_function = TRUE,
                                  grid_len = NULL,
                                  normalize = TRUE) {
  # case 1: user gave an explicit base_template
  if (!is.null(base_template) && is.null(template_name)) {
    tmpl <- as.numeric(base_template)
    tmpl[!is.finite(tmpl)] <- 0
    if (normalize) {
      nrm <- sqrt(sum(tmpl * tmpl))
      if (is.finite(nrm) && nrm > 0) {
        tmpl <- tmpl / nrm
      }
    }
    return(tmpl)
  }
  
  # case 2: wavelet by name
  if (!is.null(template_name)) {
    # Auto-cap db2/db3 template lengths to prevent over-smoothing
    # Uses natural cascade-determined lengths without resampling to very long templates
    template_name_lower <- tolower(as.character(template_name))
    
      # Use natural cascade-determined length (don't force resampling)
    # Only resample if explicitly requested (grid_len provided) or for adaptive scaling
    cascade_level <- as.integer(cascade_level)
    
    if (!is.null(grid_len)) {
      # User provided grid_len - cap it for db2/db3 to prevent over-smoothing
      grid_len <- as.integer(grid_len)
      if (template_name_lower == "db2" && grid_len > 50L) {
        warning("Capping db2 template grid_len from ", grid_len, " to 50 to prevent over-smoothing", immediate. = TRUE)
        grid_len <- 50L
      } else if (template_name_lower == "db3" && grid_len > 60L) {
        warning("Capping db3 template grid_len from ", grid_len, " to 60 to prevent over-smoothing", immediate. = TRUE)
        grid_len <- 60L
      }
    } else {
      # Use natural cascade length (no resampling)
      # Natural cascade lengths:
      # - haar level 1: ~2 points, level 2: ~4, level 3: ~8
      # - db2 level 1: ~4 points, level 2: ~8, level 3: ~16
      # - db3 level 1: ~6 points, level 2: ~12, level 3: ~24
      # Keep grid_len = NULL to use natural length
      # Only cap cascade level to prevent excessive length
      if (template_name_lower == "db2" && cascade_level > 3L) {
        warning("Capping db2 cascade_level from ", cascade_level, " to 3 to prevent over-smoothing", immediate. = TRUE)
        cascade_level <- 3L
      } else if (template_name_lower == "db3" && cascade_level > 2L) {
        warning("Capping db3 cascade_level from ", cascade_level, " to 2 to prevent over-smoothing", immediate. = TRUE)
        cascade_level <- 2L
      }
    }
    
    # Also limit cascade level for db2/db3 to prevent excessive length (even if grid_len provided)
    if (template_name_lower == "db2" && cascade_level > 3L) {
      warning("Capping db2 cascade_level from ", cascade_level, " to 3 (Consenrich-like approach)", immediate. = TRUE)
      cascade_level <- 3L
    } else if (template_name_lower == "db3" && cascade_level > 2L) {
      warning("Capping db3 cascade_level from ", cascade_level, " to 2 (Consenrich-like approach)", immediate. = TRUE)
      cascade_level <- 2L
    }
    
    # Also limit cascade level for db2/db3 to prevent excessive length
    cascade_level <- as.integer(cascade_level)
    if (template_name_lower == "db2" && cascade_level > 3L) {
      warning("Capping db2 cascade_level from ", cascade_level, " to 3 (Consenrich-like approach)", immediate. = TRUE)
      cascade_level <- 3L
    } else if (template_name_lower == "db3" && cascade_level > 2L) {
      warning("Capping db3 cascade_level from ", cascade_level, " to 2 (Consenrich-like approach)", immediate. = TRUE)
      cascade_level <- 2L
    }
    
    key <- .sp_template_key(template_name,
                            cascade_level,
                            use_scaling_function,
                            grid_len)
    if (exists(key, envir = .sp_template_cache, inherits = FALSE)) {
      return(get(key, envir = .sp_template_cache, inherits = FALSE))
    }
    
    tmpl <- build_wavelet_template(
      name          = template_name,
      cascade_level = cascade_level,
      use_scaling   = isTRUE(use_scaling_function),
      grid_len      = grid_len
    )
    tmpl <- as.numeric(tmpl)
    tmpl[!is.finite(tmpl)] <- 0
    
    if (normalize) {
      nrm <- sqrt(sum(tmpl * tmpl))
      if (is.finite(nrm) && nrm > 0) {
        tmpl <- tmpl / nrm
      }
    }
    
    assign(key, tmpl, envir = .sp_template_cache)
    return(tmpl)
  }
  
  stop("Either `base_template` or `template_name` must be provided.")
}
