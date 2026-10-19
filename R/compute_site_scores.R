#' Compute site-level scores following BAGFOOT2 logic
#'
#' BAGFOOT2-like definitions:
#' - window: +/- 200 bp around motif
#' - motif region: motifcenter +/- motif width, extended by motifbuffer
#' - flanking region: all positions outside motif region
#' - FA (independent): sum of observed flank cuts per site
#' - FPD: mean log2(obs/exp) in motif minus mean log2(obs/exp) in flanks
#'
#' @param sites data.frame with columns chr, st, ed, dir
#' @param S AggregationPlotData-like object
#' @param tdensity cut count list per chromosome
#' @param motifwidth integer, motif width
#' @param motifcenter integer, motif center position (0 = infer from motif width)
#' @param halfspan integer, half-window around motif center; BAGFOOT2-like default = 200
#' @param nuccodes nucleotide code list
#' @param ftable bias table
#' @param motifbuffer integer, motif offset around motif region; BAGFOOT2-like default = 2
#' @param pseudocount numeric, small value to avoid division by zero
#' @param trim numeric, trimming proportion for mean
#' @param outfile optional output file path
#'
#' @return data.frame with one row per site
#' @export
compute_site_scores <- function(
  sites,
  S,
  tdensity,
  motifwidth,
  motifcenter = 0,
  halfspan = 200,
  nuccodes = NA,
  ftable = NA,
  motifbuffer = 2,
  pseudocount = 1e-4,
  trim = 0.1,
  outfile = NULL
) {
  if (is.null(sites) || nrow(sites) == 0) {
    out <- data.frame()
    if (!is.null(outfile)) data.table::fwrite(out, outfile)
    return(out)
  }

  counts <- CalcObservedExpectedCuts(
    S = S,
    motif = sites,
    tdensity = tdensity,
    motifcenter = motifcenter,
    halfspan = halfspan,
    nuccodes = nuccodes,
    motifwidth = motifwidth,
    ftable = ftable
  )

  observed <- counts$observed
  expected <- counts$expected
  valid_sites <- counts$sites

  if (is.null(observed) || is.null(expected) || nrow(observed) == 0) {
    out <- data.frame()
    if (!is.null(outfile)) data.table::fwrite(out, outfile)
    return(out)
  }

  observed[is.na(observed)] <- 0
  expected[is.na(expected)] <- 0

  observed_pc <- observed + pseudocount
  expected_pc <- expected + pseudocount

  site_logratio <- log2(observed_pc / expected_pc)

  if (motifcenter == 0) {
    motifcenter_use <- floor(motifwidth / 2) + 1
  } else {
    motifcenter_use <- motifcenter
  }

  full_width <- ncol(site_logratio)
  halfpoint <- floor(full_width / 2) + 1

  xr <- seq_len(motifwidth) - motifcenter_use
  x1 <- min(xr)
  x2 <- max(xr)

  motif_start <- halfpoint + x1 - motifbuffer
  motif_end   <- halfpoint + x2 + motifbuffer

  motif_start <- max(1, motif_start)
  motif_end   <- min(full_width, motif_end)

  motif_idx <- motif_start:motif_end

  # TIGER-like flanking region = all positions outside motif region
  flank_idx <- setdiff(seq_len(full_width), motif_idx)

  if (length(motif_idx) == 0) stop("Empty motif region.")
  if (length(flank_idx) == 0) stop("Empty flanking region.")

  # FA independent: sum of observed cuts in flanks
  site_fa_sum <- rowSums(
    observed[, flank_idx, drop = FALSE],
    na.rm = TRUE
  )

  # Optional normalized flank sum
  site_fa_sum_norm <- rowSums(
    (observed_pc / expected_pc)[, flank_idx, drop = FALSE],
    na.rm = TRUE
  )

  # FPD components: mean log-ratio in motif and flanks
  site_motif_mean_logratio <- apply(
    site_logratio[, motif_idx, drop = FALSE],
    1,
    mean,
    trim = trim,
    na.rm = TRUE
  )

  site_flank_mean_logratio <- apply(
    site_logratio[, flank_idx, drop = FALSE],
    1,
    mean,
    trim = trim,
    na.rm = TRUE
  )

  site_fpd <- site_motif_mean_logratio - site_flank_mean_logratio

  out <- data.frame(
    chr = valid_sites$chr,
    st = valid_sites$st,
    ed = valid_sites$ed,
    dir = valid_sites$dir,
    site_fa_sum = site_fa_sum,
    site_fa_sum_norm = site_fa_sum_norm,
    site_motif_mean_logratio = site_motif_mean_logratio,
    site_flank_mean_logratio = site_flank_mean_logratio,
    site_fpd = site_fpd,
    observed_sum_total = rowSums(observed, na.rm = TRUE),
    expected_sum_total = rowSums(expected, na.rm = TRUE)
  )

  if (!is.null(outfile)) {
    ext <- tolower(tools::file_ext(outfile))
    if (ext == "tsv") {
      data.table::fwrite(out, outfile, sep = "\t")
    } else {
      data.table::fwrite(out, outfile)
    }
  }

  return(out)
}
