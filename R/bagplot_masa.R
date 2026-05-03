
### function including homer results

#' Generate a bagplot with optional HOMER results integration and point classification
#'
#' @param dat Data frame containing bagplot input data.
#' @param dataname1 Name of the first dataset (e.g., control).
#' @param dataname2 Name of the second dataset (e.g., treatment).
#' @param qvaluethreshold_bagplot Threshold for q-value significance (default: 0.05).
#' @param factor Factor for bagplot scaling (default: 3).
#' @param pdf Whether to save the plot as a PDF (default: TRUE).
#' @param include_homer_results Whether to include HOMER results (default: FALSE).
#' @param homer_results_file Path to the HOMER results file (e.g., knownResults.txt).
#' @param pval_threshold_homer P-value threshold for HOMER results (default: 0.05).
#' @param use_cluster_motif_db Whether to apply motif name extraction specific to the cluster_motif database (default: FALSE).
#' @return A ggplot object with the bagplot.
#' @export
gen_bagplot_tiger <- function(
  dat, dataname1 = dataname1,
  dataname2 = dataname2,
  qvaluethreshold_bagplot = 0.05,
  factor = 3,
  pdf = TRUE,
  html = FALSE,
  noBag = FALSE,
  noFence = FALSE,
  noMedianPoint = FALSE,
  pdfwidth = 16,
  pdfheight = 12,
  textsizefactor = 1.0,
  printmode = FALSE,
  include_homer_results = FALSE,
  homer_results_file = NULL,
  pval_threshold_homer = 0.05,      # Threshold for HOMER P-value significance
  use_cluster_motif_db = TRUE) {

  # Prepare data
  x = dat$lrtotal
  y = dat$fdepth1 - dat$fdepth2

  trimbp <- function(sss) {
    doit <- function(sss) {
      ttt= gregexpr('_[0-9]+bp', sss)[[1]][1];
      ifelse(ttt<1, sss, substring(sss, 1,ttt-1));
    }
    sapply(sss, doit);
  }

  mat <- data.frame(dhp=dat$lrtotal, dfd=dat$fdepth1-dat$fdepth2);
  name=trimbp(dat$motif);

  datan = data.frame(x=x, y=y, name=as.character(name),stringsAsFactors=F);
  ss <- aplpack::compute.bagplot(datan$x, datan$y, factor=factor, approx.limit = nrow(datan));

  outlier = ss$pxy.outlier;
  outer = ss$pxy.outer;
  bag = ss$pxy.bag;
  hdepths = ss$hdepths;
  dat = datan;

  findOutliersIndex<-function(data, outlier) {
    cols <- colnames(outlier);
    dat2 <- data[,cols];
    nr = nrow(outlier);
    rg = 1:nrow(dat2);
    matchingIdx = unique(sort(unlist(sapply(1:nr, function(ii) rg[(outlier[ii,2] == dat2[,2]) & (outlier[ii,1]==dat2[,1])]))));
    if (length(matchingIdx) != nr) {
      stop('cannot locate outliers in the data.');
    }
    matchingIdx;
  }

  nametoshow = rep("", nrow(datan));
  nametoshowInGray = rep("", nrow(datan));

  outlieridx <- integer(0)
  if (!is.null(outlier) && nrow(outlier) > 0) {
    outlieridx = findOutliersIndex(datan, outlier);
    nametoshow[outlieridx] = as.character(datan$name[outlieridx]);
  }

  outeridx = findOutliersIndex(datan, outer);
  inneridx = findOutliersIndex(datan, bag);
  nametoshowInGray[outeridx] = as.character(datan$name[outeridx]);

  category = vector("character", length=nrow(datan));
  category[inneridx] = "bag";
  category[outeridx] = "fence";
  category[outlieridx ] = "outlier";

  mx = mean(datan$x);
  my = mean(datan$y);

  Sx <- cov(mat);
  #D2 <- mahalanobis(mat, colMeans(mat), Sx);
  D2 <- mahalanobis(mat,apply(mat, 2, median), Sx);
  pvalue =   1-pchisq(D2, ncol(mat));
  qvalue = p.adjust(pvalue, method="BH");  # FDR adjusted p-value (BH)

  udist = sqrt((datan$x-mx)^2 + (datan$y-my)^2);

  bagplottable= data.frame(name=datan$name,
                           log2_access_ratio=datan$x,
                           deltafd = datan$y,
                           category=category,
                           hdepth= hdepths,
                           udist=udist,
                           pvalue=pvalue,
                           qvalue=qvalue,
                           D2 = D2);

  outputfilename=sprintf("bagplot_cutcount_diff_total_footprinting_depth_%s-%s",dataname2,dataname1);

  #############################################################################
  #### new classification increase/decrease/discordant ###################

  #  x/y  bag and fence
  window_x_min <- min(bagplottable$log2_access_ratio[bagplottable$category %in% c("bag", "fence")])
  window_x_max <- max(bagplottable$log2_access_ratio[bagplottable$category %in% c("bag", "fence")])
  window_y_min <- min(bagplottable$deltafd[bagplottable$category %in% c("bag", "fence")])
  window_y_max <- max(bagplottable$deltafd[bagplottable$category %in% c("bag", "fence")])

  bagplottable$direction <- "non_outlier"

  is_outlier <- bagplottable$category == "outlier"
  log2_access_ratio <- bagplottable$log2_access_ratio
  deltafd <- bagplottable$deltafd

  # increase
  increase_mask <- (
    (log2_access_ratio > window_x_min & log2_access_ratio < window_x_max & deltafd > window_y_max) | # above of the rectangle
      (log2_access_ratio > window_x_max & deltafd > window_y_min & deltafd < window_y_max) |  # right from the rectangle
      (log2_access_ratio > 0 & deltafd > 0)                                                  # quadrant ++
  )

  #  decrease
  decrease_mask <- (
    (log2_access_ratio > window_x_min & log2_access_ratio < window_x_max & deltafd < window_y_min) | # bellow the resctangle
      (log2_access_ratio < window_x_min & deltafd > window_y_min & deltafd < window_y_max) |  # left from the rectangle
      (log2_access_ratio < 0 & deltafd < 0)                                                  # quadrant --
  )

  #  discordant : out of rectangle and different signs
  outside_x <- (log2_access_ratio < window_x_min) | (log2_access_ratio > window_x_max)
  outside_y <- (deltafd < window_y_min) | (deltafd > window_y_max)
  discordant_mask <- (sign(log2_access_ratio) != sign(deltafd)) & outside_x & outside_y

  # Application
  bagplottable$direction[is_outlier & increase_mask] <- "increase"
  bagplottable$direction[is_outlier & decrease_mask] <- "decrease"
  bagplottable$direction[is_outlier & discordant_mask] <- "discordant"

  # keep original direction labels ("increase","decrease","discordant") for static plots
  # (HTML/interactive mapping will be handled separately)

  #############################################################################
  #write.csv(bagplottable, file=sprintf('%s_bagplot_output.csv', outputfilename));

  labx=sprintf("Normalized Cutcount Diff. (%s-%s)",dataname2,dataname1);
  laby=sprintf("-Footprinting Depth Diff. (%s-%s)",dataname2,dataname1);

  do_draw<-function() {

    ss <- aplpack::bagplot(datan$x, datan$y, factor=factor,show.looppoints = FALSE,show.bagpoints=FALSE,
                           show.baghull = !noBag,  show.loophull= !noFence,  show.whiskers=F,xlab=labx,ylab=laby, col.looppoint = '#000000',cex=2,pch=0, cex.lab=1.5,  cex.axis=1.5, lwd = 3);
    points(datan$x[outeridx], datan$y[outeridx], pch=16, col= "#332280");      # outer
    points(datan$x[inneridx], datan$y[inneridx], pch=18, col= "#000000");   # bag

    #PVALUETHRESHOLD = 0.05;
    points(datan$x[qvalue < qvaluethreshold_bagplot], datan$y[qvalue < qvaluethreshold_bagplot],cex=2, pch=15, col= "#FF0000");   #  p-value < 0.05;q

    text(datan$x, datan$y, labels=nametoshow,adj=1, cex=1.5);
    text(datan$x, datan$y, labels=nametoshowInGray,col="#000000A0", adj=1, cex=0.9);
    plotrange <- par("usr");
    xmin = plotrange[1];
    xmax = plotrange[2];
    ymin = plotrange[3];
    ymax = plotrange[4];
    lines(c(xmin, xmax),c(0,0));
    lines(c(0,0),c(ymin, ymax));

  }
  do_draw();

  #bagplot from aplpack, saved in file
  #outputfilename=sprintf("bagplot_cutcount_diff_total_footprinting_depth_%s_%s_TIGER",dataname2,dataname1,qvaluethreshold_bagplot);
  #tiff(filename= sprintf("%s.tiff",outputfilename),compression="zip", width=2000,height=1500, units="px", pointsize=20);
  #do_draw();
  #dev.off();

  #if (pdf) {
  # pdf(file= sprintf("%s.pdf",outputfilename), width=16,height=12, pointsize=10,useDingbats=FALSE);
  #  do_draw();
  # dev.off();
  #}


  # Generate ggplot2 visualization


  # Compute rectangle limits for bag and fence
  window_x_min <- min(bagplottable$log2_access_ratio[bagplottable$category %in% c("bag", "fence")])
  window_x_max <- max(bagplottable$log2_access_ratio[bagplottable$category %in% c("bag", "fence")])
  window_y_min <- min(bagplottable$deltafd[bagplottable$category %in% c("bag", "fence")])
  window_y_max <- max(bagplottable$deltafd[bagplottable$category %in% c("bag", "fence")])
  center<- ss$center

  #########################
  ### ADD GGPLOT + HOMER ##

  qvalue_bagplot <- if ("qvalue" %in% names(bagplottable)) {
    as.numeric(bagplottable$qvalue)
  } else {
    rep(NA_real_, nrow(bagplottable))
  }

  bagplottable$qvalue_bagplot<-qvalue_bagplot

  # Compute significance according to bagfoot qvalue threshold
  bagplottable$significant_bagplot <- !is.na(bagplottable$qvalue) & (bagplottable$qvalue < qvaluethreshold_bagplot)


  # Compute significance according to HOMER P.value threshold

  if (include_homer_results) {
    if (is.null(homer_results_file)) {
      warning("HOMER results file must be provided when include_homer_results is TRUE. Plot will be created without HOMER integration.")
      bagplottable$significant_homer <- FALSE
    } else {
      homer_results <- read.table(homer_results_file, header = FALSE, sep = "\t", fill=TRUE)
      colnames(homer_results) <- c(
        "Motif.Name", "Consensus", "P.value", "Log.P.value", "q.value",
        "nb_target_seq_with_motif", "perc_target_seq_with_motif",
        "nb_background_seq_with_motif", "perc_background_seq_with_motif"
      )
      homer_results <- homer_results[-1,]
      homer_results$P.value <- as.numeric(homer_results$P.value)
      homer_results$Motif.Name <- sapply(homer_results$Motif.Name, function(x) strsplit(x, ":")[[1]][1])
      homer_results$Motif.Name <- toupper(homer_results$Motif.Name)
      bagplottable$name <- toupper(bagplottable$name)
      significant_motifs_homer <- homer_results %>%
        filter(P.value <= pval_threshold_homer) %>%
        pull(Motif.Name)
      bagplottable$significant_homer <- bagplottable$name %in% significant_motifs_homer

    }
  } else {
    bagplottable$significant_homer <- FALSE

  }



  ########################
  ########################
  # PREPARATION  CSV
  ########################
  #merge with cluster_motif data
  file_path_src <- file.path(getwd(), "inst", "data", "v2_ig.txt")
  file_path_pkg <- system.file("data", "v2_ig.txt", package = "tiger")

  if (file.exists(file_path_src)) {
    file_path <- file_path_src
  } else if (file_path_pkg != "" && file.exists(file_path_pkg)) {
    file_path <- file_path_pkg
  } else {
    stop("ERROR: v2_ig.txt file not found in source tree (inst/data) or installed package. Searched:\n",
         file_path_src, "\n", file_path_pkg)
  }

  ########### Reading motif cluster info file #############

  # Now read and merge

  # without header pour to avoid crash
  v2_raw <- read.table(file_path, header = FALSE, sep = "\t", fill = TRUE, stringsAsFactors = FALSE)

  # first line as names
  header_names <- as.character(v2_raw[1, ])

  # keep only data
  v2 <- v2_raw[-1, ]

  # Cleaning: If header has one column in more
  if (length(header_names) > ncol(v2)) {
    # cut
    header_names <- header_names[1:ncol(v2)]
  } else if (ncol(v2) > length(header_names)) {
    # if more data, we reduce
    v2 <- v2[, 1:length(header_names)]
  }

  # Assignation of names
  colnames(v2) <- header_names

  ################################

  bagplottable$name <- as.character(bagplottable$name)
  v2$cluster       <- as.character(v2$cluster)
  v2[] <- lapply(v2, as.character)

  bagplottable <- merge(bagplottable, v2, by.x = "name", by.y = "cluster", all.x = TRUE)
  print(colnames(bagplottable))

  str(bagplottable)
  is_bad <- sapply(bagplottable, function(x) is.list(x) || is.data.frame(x))
  if (any(is_bad)) {
    warning("Non atomic columns in bagplottable: ",
            paste(names(bagplottable)[is_bad], collapse = ", "))
    bagplottable <- bagplottable[, !is_bad, drop = FALSE]
  }

  # rename columns
  # keep original direction values in bagplottable_export (used for CSV); interactive HTML will map colors separately
  bagplottable_export <- bagplottable



  # Basic column renames
  bagplottable_export <- bagplottable_export %>%
    dplyr::rename(
      motif_name = "name_IG",           # name_IG => motif_name
      delta_FA = "log2_access_ratio",  # log2_access_ratio => delta_FA
      delta_FPD = "deltafd",     # deltafd => delta_FPD
      halfspace_depth  = "hdepth",
      euclidean_distance = "udist",
      mahalanobis_distance_D2 = "D2",
      archetype_name = "name" ,          # name => archetype_name
      D2_qvalue = "qvalue_bagplot",
      D2_pvalue = "pvalue",
      D2_significant = "significant_bagplot"
    )
  if (include_homer_results) {
    bagplottable_export <- bagplottable_export %>% dplyr::rename(Homer_significant = "significant_homer")
  } else {
    if ("significant_homer" %in% names(bagplottable_export)) bagplottable_export$significant_homer <- NULL
    if ("Homer_significant" %in% names(bagplottable_export)) bagplottable_export$Homer_significant <- NULL
  }

  str(bagplottable_export)
  required_cols <- c("motif_name","delta_FA","delta_FPD","category",
                     "mahalanobis_distance_D2","halfspace_depth","euclidean_distance",
                     "D2_pvalue","D2_qvalue","direction","D2_significant",
                     "archetype_name","motif_id","source_id","tf_name","family_name",
                     "motif_type","PMID")
  if (include_homer_results) {
    required_cols <- c(required_cols, "Homer_significant")
  }
  missing <- setdiff(required_cols, names(bagplottable_export))
  if (length(missing)) {
    stop("Columns missing in bagplottable_export: ",
         paste(missing, collapse = ", "))
  }


  cols_to_select <- c(
    "motif_name",
    "delta_FA",
    "delta_FPD",
    "category",
    "mahalanobis_distance_D2",
    "halfspace_depth",
    "euclidean_distance",
    "D2_pvalue",
    "D2_qvalue",
    "direction",
    "D2_significant",
    "archetype_name",
    "motif_id",
    "source_id",
    "tf_name",
    "family_name",
    "motif_type",
    "PMID"
  )
  if (include_homer_results) {
    cols_to_select <- append(cols_to_select, "Homer_significant", after = 11)
  }
  bagplottable_export <- bagplottable_export %>% dplyr::select(all_of(cols_to_select))


  # Save results to CSV (with all qvalues)
  outputfilename = sprintf("bagplot_cutcount_diff_total_footprinting_depth_%s_%s", dataname2, dataname1)
  write.csv(bagplottable_export, file = sprintf('%s_bagplot_output_TIGER.csv', outputfilename), row.names=FALSE)

  window_x_min <- min(bagplottable$log2_access_ratio[bagplottable$category %in% c("bag", "fence")], na.rm=TRUE)
  window_x_max <- max(bagplottable$log2_access_ratio[bagplottable$category %in% c("bag", "fence")], na.rm=TRUE)
  window_y_min <- min(bagplottable$deltafd[bagplottable$category %in% c("bag", "fence")], na.rm=TRUE)
  window_y_max <- max(bagplottable$deltafd[bagplottable$category %in% c("bag", "fence")], na.rm=TRUE)


  label_color <- function(direction, category) {
    if (direction == "increase") return("tomato2")
    if (direction == "decrease") return("royalblue1")
    if (direction == "non_outlier" & category == "fence") return("grey")
    if (direction == "non_outlier" & category == "bag") return("black")
    return("black")
  }

  bagplottable$label_col <- mapply(label_color, bagplottable$direction, bagplottable$category)

  # Plot
  baseplot <- ggplot(bagplottable, aes(log2_access_ratio, deltafd)) +

    # Rectangle around bag+fence
    annotate(
      "rect",
      xmin = window_x_min, xmax = window_x_max,
      ymin = window_y_min, ymax = window_y_max,
      alpha = 0.15, fill = "grey46", color = "grey35",
      linetype = "dashed", linewidth = 1
    ) +

    # Bag polygon
    geom_polygon(
      data = {
        df <- bagplottable %>% filter(category == "bag")
        df[chull(df$log2_access_ratio, df$deltafd), ]
      },
      aes(x = log2_access_ratio, y = deltafd),
      fill = NA, color = "black", linewidth = 0.5
    ) +
    # Fence polygon
    geom_polygon(
      data = {
        df <- bagplottable %>% filter(category == "fence")
        df[chull(df$log2_access_ratio, df$deltafd), ]
      },
      aes(x = log2_access_ratio, y = deltafd),
      fill = NA, color = "black", linetype = "dashed", linewidth = 0.5
    ) +

    # Bag points
    geom_point(
      data = bagplottable %>% filter(category == "bag"),
      color = "gray50", fill = "gray50", pch = 20, size = 1
    ) +
    # Fence points
    geom_point(
      data = bagplottable %>% filter(category == "fence"),
      color = "gray50", fill = "gray50", pch = 20, size = 1
    ) +

    # Outlier - square red/blue/black orange cercle  if bagfoot and homer significant bagfoot
    geom_point(
      data = bagplottable %>% filter(category == "outlier" & significant_bagplot & significant_homer & direction == "increase"),
      color = "darkGoldenrod1", fill = "red3", pch = 22, size = 4, stroke = 2
    ) +
    geom_point(
      data = bagplottable %>% filter(category == "outlier" & significant_bagplot & significant_homer & direction == "decrease"),
      color = "darkGoldenrod1", fill = "darkblue", pch = 22, size = 4, stroke = 2
    ) +
    geom_point(
      data = bagplottable %>% filter(category == "outlier" & significant_bagplot & significant_homer & direction == "discordant"),
      color = "darkGoldenrod1", fill = "black", pch = 22, size = 2, stroke = 2
    ) +
    # Outlier - circle red/blue/black orange cercle  if bagfoot outlier and homer significant
    geom_point(
      data = bagplottable %>% filter(category == "outlier" & !significant_bagplot & significant_homer & direction == "increase"),
      color = "darkGoldenrod1", fill = "white", pch = 21, size = 3, stroke = 2
    ) +
    geom_point(
      data = bagplottable %>% filter(category == "outlier" & !significant_bagplot & significant_homer & direction == "decrease"),
      color = "darkGoldenrod1", fill = "white", pch = 21, size = 3, stroke = 2
    ) +
    geom_point(
      data = bagplottable %>% filter(category == "outlier" & !significant_bagplot & significant_homer & direction == "discordant"),
      color = "darkGoldenrod1", fill = "black", pch = 21, size = 2, stroke = 2
    ) +
    # Outlier - square red/blue/black
    geom_point(
      data = bagplottable %>% filter(category == "outlier" & significant_bagplot & !significant_homer & direction == "increase"),
      color = "red3", fill = "red3", pch = 22, size = 4, stroke = 1
    ) +
    geom_point(
      data = bagplottable %>% filter(category == "outlier" & significant_bagplot & !significant_homer & direction == "decrease"),
      color = "darkblue", fill = "darkblue", pch = 22, size = 4, stroke = 1
    ) +
    geom_point(
      data = bagplottable %>% filter(category == "outlier" & significant_bagplot & !significant_homer & direction == "discordant"),
      color = "black", fill = "black", pch = 22, size = 1, stroke = 1
    ) +

    # Outlier - white circle if non-significatif bagplot
    geom_point(
      data = bagplottable %>% filter(category == "outlier" & !significant_bagplot),
      aes(color = direction),
      fill = "white", pch = 21, size = 1, stroke = 1
    ) +

    # Labels for outliers
    geom_text(
      data = bagplottable %>% filter(category == "outlier" & direction!="discordant"),
      aes(label = name_IG, color = direction),
      hjust = 1.2, size = 5, show.legend = FALSE
    ) +
    geom_text(
      data = bagplottable %>% filter(category == "outlier" & direction=="discordant"),
      aes(label = name_IG, color = direction),
      hjust = 1.2, size = 3, show.legend = FALSE
    ) +

    # Median point
    geom_point(aes(x = center[1], y = center[2]), color = "brown", size = 4, shape = 18) +

    # axes at 0
    #geom_hline(yintercept = 0, color = "black", linewidth = 0.8) +
    #geom_vline(xintercept = 0, color = "black", linewidth = 0.8) +
    scale_color_manual(
      values = c("increase" = "tomato2", "decrease" = "royalblue1", "discordant" = "black"),
      guide = "none"
    ) +
    scale_fill_manual(
      values = c("increase" = "tomato2", "decrease" = "royalblue1", "discordant" = "black"),
      guide = "none"
    ) +
    labs(
      #title = bquote( .(dataname1) ~ "vs" ~ .(dataname2)),
      x = bquote(Delta ~ " FA : " * .(dataname1) * " minus " * .(dataname2)),
      y = bquote(-Delta ~ " FPD : " * .(dataname1) * " minus " * .(dataname2))
    )+
    theme_classic()+theme(
      axis.title.x = element_text(size = 40),
      axis.title.y = element_text(size = 40),
      plot.title = element_text(size = 54, face = "bold", hjust = 0.5),
      axis.text.x = element_text(size = 16),  # Increase size of x-axis tick labels
      axis.text.y = element_text(size = 16)   # Increase size of y-axis tick labels
    )

  # Marginal boxplots
  #bagplot <- ggMarginal(
  # baseplot,
  #type = "boxplot",
  #margins = "both",
  #size = 10,
  #xparams = list(fill = "purple", alpha = 0.3),
  #yparams = list(fill = "orange", alpha = 0.3)
  #)

  bagplot <- baseplot

  if (pdf) {

    ggsave(sprintf("%s_TIGER_zoomIn.pdf", outputfilename), plot = bagplot,width=36, height=24, pointsize=24, useDingbats=FALSE)
    ggsave(sprintf("%s_TIGER.pdf", outputfilename), plot = bagplot,width=18, height=14, pointsize=14, useDingbats=FALSE)

  }

  ggsave(sprintf("%s_TIGER.png", outputfilename),
         plot = bagplot,
         width = 18, height = 14
         , units = "in",
         dpi = 600)

  # saving SVG
  svg(sprintf("%s_TIGER.svg", outputfilename),
      width = 16, height = 12)
  print(bagplot)
  dev.off()

  # --- plot interactif & formatting for ggplotly ---
  if (include_homer_results) {
    bagplottable_export <- bagplottable_export %>%
      mutate(
        delta_FA          = as.numeric(delta_FA),
        delta_FPD         = as.numeric(delta_FPD),
        direction         = as.character(direction),
        motif_name        = as.character(motif_name),
        Homer_significant = as.logical(Homer_significant),
        D2_pvalue         = as.numeric(D2_pvalue),
        D2_qvalue         = as.numeric(D2_qvalue)
      )
  } else {
    bagplottable_export <- bagplottable_export %>%
      mutate(
        delta_FA          = as.numeric(delta_FA),
        delta_FPD         = as.numeric(delta_FPD),
        direction         = as.character(direction),
        motif_name        = as.character(motif_name),
        D2_pvalue         = as.numeric(D2_pvalue),
        D2_qvalue         = as.numeric(D2_qvalue)
      )
  }

  # tooltip prep : column text (include HOMER only if requested)
  if (include_homer_results) {
    bagplottable_export <- bagplottable_export %>%
      dplyr::mutate(
        tooltip = paste0(
          "Motif: ", motif_name,
          "<br>HOMER: ", as.character(Homer_significant),
          "<br>DeltaFPD: ", round(delta_FPD, 3),
          "<br>DeltaFA: ", round(delta_FA, 3),
          "<br>pvalue (chi2): ", signif(D2_pvalue, 3),
          "<br>qvalue (BH): ", signif(D2_qvalue, 3)
        )
      )
  } else {
    bagplottable_export <- bagplottable_export %>%
      dplyr::mutate(
        tooltip = paste0(
          "Motif: ", motif_name,
          "<br>DeltaFPD: ", round(delta_FPD, 3),
          "<br>DeltaFA: ", round(delta_FA, 3),
          "<br>pvalue (chi2): ", signif(D2_pvalue, 3),
          "<br>qvalue (BH): ", signif(D2_qvalue, 3)
        )
      )
  }

  # minimal version
  baseplot_interactive <- ggplot(bagplottable_export,
                                 aes(delta_FA, delta_FPD)) +
    geom_point(
      aes(color = direction, text = tooltip),
      size = 2
    ) +
    theme_classic() +
    scale_color_manual(
      values = c("increase" = "tomato2",
                 "decrease" = "royalblue1",
                 "discordant" = "black"),
      guide = "none"
    ) +
    labs(
      x = paste0("\u0394FA : ", dataname1, " minus ", dataname2),
      y = paste0("- \u0394FPD : ", dataname1, " minus ", dataname2)
    )+
    theme_classic() + theme(
      axis.title.x = element_text(size = 30),
      axis.title.y = element_text(size = 30),
      plot.title = element_text(size = 34, face = "bold", hjust = 0.5),
      axis.text.x = element_text(size = 16),
      axis.text.y = element_text(size = 16)
    ) +
    scale_color_manual(values=c("increase"="tomato2","decrease"="royalblue1","discordant"="black"), guide="none")
  #


  # Export HTML interactive (no visible labels)
  if (html) {
    library(plotly)
    library(htmlwidgets)
    interactive_plot <- ggplotly(baseplot_interactive, tooltip = "text")
    saveWidget(interactive_plot, sprintf("%s_TIGER_interactive.html", outputfilename))
  }


  return(bagplot)

}


