#' Generate EVT-based p-values for delta FA and delta FPD, export CSV and barplots
#'
#' @param dat Data frame (same format as gen_bagplot_bagfoot2 input).
#' @param dataname1 Name of condition 1 (e.g., control).
#' @param dataname2 Name of condition 2 (e.g., treatment).
#' @param pval_threshold P-value threshold for significance (default: 0.05).
#' @param include_homer_results Whether to include HOMER results (default: FALSE).
#' @param homer_results_file Path to HOMER knownResults.txt.
#' @param pval_threshold_homer P-value threshold for HOMER (default: 0.05).
#' @param top_n Max bars per panel (default: 10). Shows all if fewer significant.
#' @return Invisible list: $data, $plot_up, $plot_dn.
#' @export
gen_evt_barplot_bagfoot2 <- function(
  dat,
  dataname1             = "cond1",
  dataname2             = "cond2",
  pval_threshold        = 0.05,
  include_homer_results = FALSE,
  homer_results_file    = NULL,
  pval_threshold_homer  = 0.05,
  top_n                 = 10
) {

  # -- 0. Helpers -------------------------------------------------------------

  trimbp <- function(sss) {
    sapply(sss, function(s) {
      ttt <- gregexpr('_[0-9]+bp', s)[[1]][1]
      ifelse(ttt < 1, s, substring(s, 1, ttt - 1))
    })
  }

  gev_pvalue <- function(x) {
    if (!requireNamespace("extRemes", quietly = TRUE))
      stop("Please install 'extRemes': install.packages('extRemes')")
    fit    <- extRemes::fevd(abs(x), type = "GEV")
    params <- fit$results$par
    1 - extRemes::pevd(
      q     = abs(x),
      loc   = params["location"],
      scale = params["scale"],
      shape = params["shape"],
      type  = "GEV"
    )
  }

  # -- 1. Extract vectors -----------------------------------------------------

  delta_FA   <- dat$lrtotal
  delta_FPD  <- dat$fdepth1 - dat$fdepth2
  motif_name <- as.character(trimbp(dat$motif))

  # -- 2. EVT p-values --------------------------------------------------------

  pval_FA  <- gev_pvalue(delta_FA)
  pval_FPD <- gev_pvalue(delta_FPD)

  # -- 3. Build table and deduplicate by AC code (keep best pval_FA) ----------

  result <- data.frame(
    cluster_id = motif_name,
    delta_FA   = delta_FA,
    delta_FPD  = delta_FPD,
    pval_FA    = pval_FA,
    pval_FPD   = pval_FPD,
    stringsAsFactors = FALSE
  )

  result <- result %>%
    dplyr::group_by(cluster_id) %>%
    dplyr::slice_min(pval_FA, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()

  # -- 4. Merge with v2_ig.txt to get name_IG ---------------------------------

  file_path_src <- file.path(getwd(), "inst", "data", "v2_ig.txt")
  file_path_pkg <- system.file("data", "v2_ig.txt", package = "bagfoot2")

  file_path <- if (file.exists(file_path_src)) {
    file_path_src
  } else if (file_path_pkg != "" && file.exists(file_path_pkg)) {
    file_path_pkg
  } else {
    stop("v2_ig.txt not found in inst/data or installed package.")
  }

  v2_raw       <- read.table(file_path, header = FALSE, sep = "\t",
                             fill = TRUE, stringsAsFactors = FALSE)
  header_names <- as.character(v2_raw[1, ])
  v2           <- v2_raw[-1, ]
  n            <- min(length(header_names), ncol(v2))
  v2           <- v2[, 1:n]
  colnames(v2) <- header_names[1:n]
  v2[]         <- lapply(v2, as.character)

  result <- merge(result, v2, by.x = "cluster_id", by.y = "cluster", all.x = TRUE)

  result$motif_name <- ifelse(
    !is.na(result$name_IG) & result$name_IG != "",
    result$name_IG,
    result$cluster_id
  )

  result <- result %>%
    dplyr::group_by(motif_name) %>%
    dplyr::slice_min(pval_FA, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()

  # -- 5. HOMER integration ---------------------------------------------------

  result$Homer_significant <- FALSE

  if (include_homer_results) {
    if (is.null(homer_results_file)) {
      warning("homer_results_file is NULL: Homer_significant set to FALSE.")
    } else {
      homer           <- read.table(homer_results_file, header = FALSE,
                                    sep = "\t", fill = TRUE)
      colnames(homer) <- c("Motif.Name", "Consensus", "P.value", "Log.P.value",
                           "q.value", "nb_target", "perc_target",
                           "nb_background", "perc_background")
      homer           <- homer[-1, ]
      homer$P.value   <- as.numeric(homer$P.value)
      homer$ac_code   <- toupper(sapply(homer$Motif.Name,
                                        function(x) strsplit(x, ":")[[1]][1]))
      sig_homer                <- homer$ac_code[homer$P.value <= pval_threshold_homer]
      result$Homer_significant <- toupper(result$cluster_id) %in% sig_homer
    }
  }

  # -- 6. Export CSV ----------------------------------------------------------

  outputfilename <- sprintf("evt_barplot_%s_%s", dataname2, dataname1)

  csv_cols <- c("motif_name", "cluster_id",
                "delta_FA", "delta_FPD",
                "pval_FA",  "pval_FPD")
  if (include_homer_results) csv_cols <- c(csv_cols, "Homer_significant")

  write.csv(result[, csv_cols, drop = FALSE],
            file      = sprintf("%s_TIGER_evt.csv", outputfilename),
            row.names = FALSE)

  # -- 7. Select top_n for each direction ------------------------------------

  sig_up <- result[!is.na(result$pval_FA) &
                     result$pval_FA < pval_threshold &
                     result$delta_FA > 0, ]

  sig_dn <- result[!is.na(result$pval_FA) &
                     result$pval_FA < pval_threshold &
                     result$delta_FA < 0, ]

  sig_up <- sig_up[order(-sig_up$delta_FA), ]
  sig_up <- sig_up[seq_len(min(nrow(sig_up), top_n)), ]

  sig_dn <- sig_dn[order(sig_dn$delta_FA), ]
  sig_dn <- sig_dn[seq_len(min(nrow(sig_dn), top_n)), ]

  # Display order: largest |delta_FA| on top
  sig_up <- sig_up[order(sig_up$delta_FA), ]
  sig_dn <- sig_dn[order(abs(sig_dn$delta_FA)), ]

  # Actual number of bars = max across both panels (≤ top_n)
  # Both panels are padded to this value so bar heights are identical.
  n_actual <- max(nrow(sig_up), nrow(sig_dn), 1L)

  # -- 7b. GLOBAL scales shared by both panels --------------------------------

  sq_gap <- 0.06

  all_abs        <- c(abs(sig_up$delta_FA), abs(sig_dn$delta_FA))
  global_abs_max <- if (length(all_abs) > 0) max(all_abs, na.rm = TRUE) else 1

  y_sq_fpd_g    <- global_abs_max * (1 + sq_gap)
  y_sq_hom_g    <- global_abs_max * (1 + sq_gap * 2.5)
  y_lim_g       <- global_abs_max * (1 + sq_gap * 4.5)
  fill_limits_g <- c(0, pval_threshold)

  # -- 7c. Padding helper -----------------------------------------------------
  # Pads df to exactly n rows with invisible placeholder rows so that both
  # panels have the same number of factor levels → ggplot allocates the same
  # physical bar height regardless of how many real motifs were found.
  # IMPORTANT: yval must already exist in df before calling pad_to_n so that
  # the template row and the real rows have identical column sets (no rbind
  # mismatch).

  pad_to_n <- function(df, n) {
    n_real <- nrow(df)
    if (n_real >= n) return(df)
    n_pad    <- n - n_real
    template <- df[1, , drop = FALSE]
    template[1, ] <- NA
    pad_rows <- template[rep(1, n_pad), , drop = FALSE]
    # Unique invisible names (spaces of increasing length)
    pad_rows$motif_name <- strrep(" ", seq_len(n_pad))
    pad_rows$yval       <- 0
    rbind(pad_rows, df)   # padding at bottom = bottom of flipped plot
  }

  # -- 8. Build barplots ------------------------------------------------------

  fill_col      <- c("#99000D", "#EF3B2C", "#FC9272", "#FDBBA1", "grey")
  col_fpd       <- "#B8560A"   # dark orange  — delta FPD significant
  col_homer     <- "#B8960A"   # yellow ochre — Homer/MEME enriched
  label_size    <- 4
  label_border  <- 0.6
  bar_height_in <- 0.35
  margin_in     <- 1.8

  build_barplot <- function(df, title_str) {
    if (nrow(df) == 0) return(NULL)

    # yval must be set BEFORE pad_to_n so template row has the same columns
    df$yval <- abs(df$delta_FA)

    # Pad to n_actual — both panels now have identical number of factor levels
    df      <- pad_to_n(df, n_actual)
    df$yval[is.na(df$yval)] <- 0

    # Factor levels: first = bottom, last = top (coord_flip)
    df$motif_name <- factor(df$motif_name, levels = df$motif_name)

   # rows_fpd <- df[!is.na(df$pval_FPD) & df$pval_FPD < pval_threshold, , drop = FALSE]
    rows_fpd <- df[!is.na(df$pval_FPD) & df$pval_FPD < pval_threshold &
                     !is.na(df$delta_FPD) & !is.na(df$delta_FA) &
                     sign(df$delta_FPD) != sign(df$delta_FA), , drop = FALSE]

    rows_hom <- df[!is.na(df$Homer_significant) & df$Homer_significant == TRUE, , drop = FALSE]

    p <- ggplot2::ggplot(df, ggplot2::aes(x = motif_name, y = yval, fill = pval_FA)) +
      ggplot2::geom_col(width = 0.8) +
      ggplot2::coord_flip(clip = "off", ylim = c(0, y_lim_g)) +
      ggplot2::scale_fill_gradientn(
        colours  = fill_col,
        limits   = fill_limits_g,
        oob      = scales::squish,
        #name     = expression(italic(p)*"-value"),
        name     = expression(atop(Delta*"FA", italic(p)*"-value")),
        na.value = "transparent"   # padding rows: NA fill → invisible bar
      ) +
      ggplot2::scale_x_discrete(expand = ggplot2::expansion(add = c(0.5, 0.5)))

    # FPD label-box: hollow rectangle, colored border, "FPD" text
    if (nrow(rows_fpd) > 0) {
      p <- p + ggplot2::geom_label(
        data        = rows_fpd,
        mapping     = ggplot2::aes(x = motif_name, y = y_sq_fpd_g, label = "FPD"),
        color       = col_fpd,
        fill        = "white",
        label.size  = label_border,
        size        = label_size,
        fontface    = "bold",
        inherit.aes = FALSE
      )
    }

    # ME label-box: hollow rectangle, colored border, "ME" text
    if (nrow(rows_hom) > 0) {
      p <- p + ggplot2::geom_label(
        data        = rows_hom,
        mapping     = ggplot2::aes(x = motif_name, y = y_sq_hom_g, label = "ME"),
        color       = col_homer,
        fill        = "white",
        label.size  = label_border,
        size        = label_size,
        fontface    = "bold",
        inherit.aes = FALSE
      )
    }

    p +
      ggplot2::theme_classic() +
      ggplot2::labs(title = title_str, x = NULL,
                    y = expression("|" * Delta * "FA|")) +
      ggplot2::theme(
        plot.margin  = ggplot2::margin(5.5, 60, 5.5, 5.5),
        axis.text.y  = ggplot2::element_text(size = 14),
        axis.text.x  = ggplot2::element_text(size = 14),
        axis.title   = ggplot2::element_text(size = 15),
        plot.title   = ggplot2::element_text(size = 16, face = "bold"),
        legend.text  = ggplot2::element_text(size = 13),
        legend.title = ggplot2::element_text(size = 13),
        legend.key   = ggplot2::element_blank()
      )
  }

  p_up <- build_barplot(sig_up, dataname2)
  p_dn <- build_barplot(sig_dn, dataname1)

  # -- 9. Save combined PDF + TIFF --------------------------------------------

  # Height adapts to actual number of bars (max across both panels, ≤ top_n)
  plot_height <- n_actual * bar_height_in + margin_in

  if (!is.null(p_up) || !is.null(p_dn)) {
    if (!requireNamespace("patchwork", quietly = TRUE))
      stop("Please install 'patchwork': install.packages('patchwork')")

    left_panel  <- if (!is.null(p_up)) p_up else patchwork::plot_spacer()
    right_panel <- if (!is.null(p_dn)) p_dn else patchwork::plot_spacer()

    combined <- left_panel + right_panel +
      patchwork::plot_layout(ncol = 2, widths = c(1, 1), guides = "collect")

    ggplot2::ggsave(sprintf("%s_FA_TIGER.pdf",  outputfilename),
                    plot = combined, width = 18, height = plot_height)

    ggplot2::ggsave(sprintf("%s_FA_TIGER.tiff", outputfilename),
                    plot = combined, width = 18, height = plot_height,
                    dpi = 300, compression = "lzw")
  }

  message("Done. Files saved with prefix: ", outputfilename)

  return(invisible(list(
    data    = result[, csv_cols, drop = FALSE],
    plot_up = p_up,
    plot_dn = p_dn
  )))
}
