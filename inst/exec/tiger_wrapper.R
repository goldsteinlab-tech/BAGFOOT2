#!/usr/bin/env Rscript

# ==============================================================================
# TIGER WRAPPER - UNEQUAL REPLICATES & ROBUST HANDLING
# ==============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(tiger)
  library(hash)
  library(data.table)
  library(Cairo)
  library(aplpack)
  library(stringr)
  library(dplyr)
  library(ggplot2)
  library(plotly)
})

# ============================ 1. OPTIONS DEFINITION ===========================
option_list <- list(

  make_option(c("-h", "--help"), action="store_true", default=FALSE, help="Show this help message and exit"),
  
  # --- INPUT / OUTPUT ---
  make_option(c("-b", "--bam_dir"), type="character", help="Directory containing merged BAM files and the replicas if Homer is used"),
  make_option(c("-p", "--peak_dir"), type="character", help="Directory containing peak files"),
  make_option(c("-o", "--output_dir"), type="character", help="Output directory (will be created if not exists)"),
  make_option(c("-m", "--mm_files_dir"), type="character", help="Directory containing motif files"),

  # --- CONDITIONS & GENOME ---
  make_option(c("-c", "--control"), type="character", help="Control condition name (e.g., fed24h)"),
  make_option(c("-t", "--treatment"), type="character", help="Treatment condition name (e.g., fasted24h)"),
  make_option(c("-g", "--genome"), type="character", default="mm10", help="Genome version [default: mm10]"),
  make_option(c("--project"), type="character", default="tiger_proj", help="Project name directory containing bam_merged and macs2 folders"),

  # --- FILE SUFFIXES ---
  make_option(c("--peak_suffix"), type="character", default="_peaks.narrowPeak", help="Suffix for peak files"),
  make_option(c("--bam_suffix"), type="character", default="_merged.bam", help="Suffix for BAM files"),

  # --- TIGER PARAMETERS ---
  # Note: default is NULL to allow processing ALL motifs if not specified
  make_option(c("--max_motifs"), type="integer", default=NULL, help="Max number of motifs to process. If not set, ALL motifs are processed."),
  #make_option(c("--with_map"), type="character", default="FALSE", help="Use mappability (TRUE/FALSE) "),
  make_option(c("--cores"), type="integer", default=2, help="Number of CPU cores"),

  # --- HOMER / REPLICATES ---
  make_option(c("--homer"), type="character", default="FALSE", help="Run HOMER analysis (TRUE/FALSE) [default: %default]"),
  make_option(c("--homer_bin"), type="character", default=NULL, help="Path to HOMER binaries folder (required if --homer TRUE)"),

  # Separate options for replicates to handle unequal designs (e.g., 2 vs 4)
  make_option(c("--reps_control"), type="character", default="1", 
              help="List of replicate indices for Control (e.g. '1,2,3' or '4,5') [default: 1]"),
  
  make_option(c("--reps_treatment"), type="character", default="1", 
              help="List of replicate indices for Treatment (e.g. '1,2,3' or '4,5') [default: 1]"),

  make_option(c("--peak_l2fc_cutoff"), type="double", default=0.5, help="Log2 Fold Change threshold to select significant differential peaks [default: %default]"),
  make_option(c("--peak_padj_cutoff"), type="double", default=0.05, help="Adjusted p-value (FDR) threshold to select significant differential peaks [default: %default]")
)
opt_parser <- OptionParser(option_list=option_list, usage = "tiger [options]", add_help_option = FALSE)
opt <- parse_args(opt_parser, print_help_and_exit = FALSE)

# ============================ 2. SETUP & CLEANING =============================
# --- personnalized help ---
if (opt$help) {
  cat("TIGER: a computational tool to predict transcription factor activity from ATAC-seq\n")
  cat("\n")
  cat("Version:   1.0.1\n")
  cat("About:     developed in the Goldstein Lab.\n")
  cat("Code:      https://github.com/goldsteinlab-tech/TIGER\n")
  cat("\n")
  
  # --- hide secret options ---
  # filter before affichage
  visible_options <- option_list[ !sapply(option_list, function(x) "--max_motifs" %in% x@long_flag) ]
  
  # temp
  dummy_parser <- OptionParser(option_list=visible_options, add_help_option=FALSE)
  
  print_help(dummy_parser)
  quit(status = 0)
} 

# --- SUPPRIMEZ LES LIGNES QUI ETAIENT ICI (print_help, quit, et l'accolade en trop) ---

# Create Output Directory
if (!is.null(opt$output_dir)) {
    dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)
    bg_dir <- normalizePath(opt$output_dir)
    setwd(bg_dir)
} else {
    stop("Error: --output_dir is required.")
}

# Helper for boolean arguments
to_logical <- function(x) {
  if (is.null(x)) return(FALSE)
  x <- as.character(x)
  if (toupper(x) %in% c("TRUE", "T", "YES", "1")) return(TRUE)
  return(FALSE)
}

# Helper to parse replicate indices "4,5" -> c(4,5)
parse_reps <- function(x) {
  if (is.null(x) || x == "") return(1)
  # remove space and cut at the comma
  nums <- as.integer(unlist(strsplit(gsub(" ", "", as.character(x)), ",")))
  return(nums)
}

#use_map <- to_logical(opt$with_map)
run_homer_flag <- to_logical(opt$homer)
use_map <- ""

# System Options (SILENT MODE)
options("tiger.mc.cores" = opt$cores)
options("tiger.mc.cores.motif" = opt$cores)
options("tiger.debug" = FALSE) # Disable verbose debug

# Global Variables assignment
bam_dir <- normalizePath(opt$bam_dir)
peak_dir <- normalizePath(opt$peak_dir)
mm_files_dir <- opt$mm_files_dir
genome <- opt$genome
set1 <- opt$control
set2 <- opt$treatment

# Explicitly define 'mm' in Global Environment for TIGER internals
mm <- genome
assign("mm", genome, envir = .GlobalEnv)

# ============================ 3. BAM & PEAK PROCESSING ========================
cat("\n>>> Step 1: Preparing BAM and Peak files...\n")

all_set <- c(set1, set2)

for (treat in all_set) {
  # BAM Processing using --bam_suffix
  bamfile <- file.path(bam_dir, paste0(treat, opt$bam_suffix))

  if (!file.exists(bamfile)) {
      stop(paste("Error: Missing BAM file:", bamfile))
  }

  cat(paste("   Processing:", treat, "\n"))

  # TIGER CutCount
  assign(paste0("cc_", treat), countReadsBAM(bamfile, mm))

  # Use 'treat' as name, not the BAM filename
  cutcountfile <- makeCutCountBAMWithName(bamfile, name = treat, refgenome = mm)

  # Hexamer & Bias Correction
  map_path <- ""
  #if (use_map) {
   # map_path <- file.path(mm_files_dir, paste0(mm, "Mappability"))
  #}

  # Hexamer filename
  hex_outfile <- paste0("Hexamer_", treat, "_", mm, "_withoutMap.txt")

  # Silent execution
  invisible(capture.output({
    MakeBiasCorrectionTableBAM(
      bamfile = bamfile,
      outfile = hex_outfile,
      refgenome = mm,
      np = 6,
      atac = TRUE,
      mapdir = ""
    )
  }))

  # Peak Processing using --peak_suffix
  peak_file_in <- file.path(peak_dir, paste0(treat, opt$peak_suffix))
  if (!file.exists(peak_file_in)) stop(paste("Error: Peak file not found:", peak_file_in))

  macs <- read.table(peak_file_in, sep="\t", header=FALSE, stringsAsFactors=FALSE)
  # Basic MACS columns check
  if (ncol(macs) >= 3) {
      names(macs)[1:3] <- c('chr','st','ed')
  } else {
      stop(paste("Error: Invalid peak file format (not enough columns):", peak_file_in))
  }

  nr <- nrow(macs)

  hotspot <- data.frame(
    ID = 1:nr,
    chrom. = macs$chr,
    start = macs$st + 1,
    end = macs$ed,
    MaxD = rep(0, nr), AveD = rep(0, nr), Zscore = rep(0, nr), pvalue = rep(0, nr)
  )
  write.table(hotspot, file=file.path(bg_dir, paste0(treat, "_peaks.csv")), sep=",", row.names=FALSE)
}

# Clean temp files
temp_f <- list.files(bg_dir, pattern = glob2rx("temp*.dat"), full.names = TRUE)
if(length(temp_f) > 0) file.remove(temp_f)

# ============================ 4. MERGING COUPLES ==============================
cat("\n>>> Step 2: Merging couples...\n")

# Combine Hotspots
assign(paste0(set1, "_", set2, "_peaks"), combineTwoHotspots(
  file.path(bg_dir, paste0(set1, "_peaks.csv")),
  file.path(bg_dir, paste0(set2, "_peaks.csv")),
  set1, set2
))

pooled_csv <- file.path(bg_dir, paste0("pooled_", set1, "_", set2, "_hotspot.csv"))
pooled <- read.csv(pooled_csv)
if ("ID.1" %in% colnames(pooled)){
  pooled <- pooled %>% dplyr::select(-ID.1)
  write.csv(pooled, pooled_csv, row.names = FALSE)
}

# Merge Frequency Tables
assign(paste0("Hexamer_", set1, "_", set2), merge_two_frequency_tables(
  paste0("Hexamer_", set1, "_", mm, "_withoutMap.txt"),
  paste0("Hexamer_", set2, "_", mm, "_withoutMap.txt"),
  paste0("Hexamer_", set1, "_", set2, "_merged_withoutMap.csv")
))
# ============================ 5. HOMER (CONDITIONAL) ==========================
homer_results_file <- NULL

if (run_homer_flag) {
  homer_out_dir <- file.path(bg_dir, "homer_motif_res")
  known_results <- file.path(homer_out_dir, "knownResults.txt")

  if (dir.exists(homer_out_dir) && file.exists(known_results)) {
    cat("\n>>> Step 3: HOMER results already present. Skipping HOMER analysis.\n")
    homer_results_file <- known_results
  } else {
    cat("\n>>> Step 3: HOMER Analysis enabled...\n")

    if (is.null(opt$homer_bin)) stop("Error: --homer_bin required if --homer is TRUE")

    # 1. Parsing des indices de réplicats (ex: "4,5" -> c(4, 5))
    ctrl_indices  <- parse_reps(opt$reps_control)
    treat_indices <- parse_reps(opt$reps_treatment)

    cat(paste("   Control Indices:  ", paste(ctrl_indices, collapse=","), "\n"))
    cat(paste("   Treatment Indices:", paste(treat_indices, collapse=","), "\n"))

    # 2. Conversion du fichier de pics
    process_hotspot_file(
      pooled_csv,
      file.path(bg_dir, paste0("pooled_", set1, "_", set2, "_hotspot_modify.txt"))
    )

    # 3. Création des Tag Directories
    cat("   Creating Tag Directories in BAM folder (Pattern: TD_Condition_Index)...\n")
    
    make_tag_dirs_custom <- function(condition, indices, bam_dir_path, out_dir_path, homer_path) {
      dirs_created <- c()
      # On boucle sur les indices RÉELS (ex: 4, 5)
      for (i in indices) {
        bam_name <- paste0(condition, "_", i, ".bam")
        bam_path <- file.path(bam_dir_path, bam_name)
        
        # MODIFICATION 1 : Ajout du préfixe TD_
        tag_dir_name <- paste0("TD_", condition, "_", i)
        
        # Le chemin de sortie dépend de l'argument out_dir_path passé à la fonction
        tag_dir_path <- file.path(out_dir_path, tag_dir_name)
        
        if (!file.exists(bam_path)) {
          stop(paste("Error: Expected replicate BAM not found:", bam_path))
        }
        
        cmd <- paste(file.path(homer_path, "makeTagDirectory"), tag_dir_path, bam_path, "-format sam")
        
        if (!dir.exists(tag_dir_path)) {
          cat(paste("     Processing:", bam_name, "->", tag_dir_name, "\n"))
          system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
        } else {
          cat(paste("     [Exists]", tag_dir_name, "\n"))
        }
        dirs_created <- c(dirs_created, tag_dir_path)
      }
      return(dirs_created)
    }

    # MODIFICATION 2 : On passe 'bam_dir' comme dossier de sortie au lieu de 'bg_dir'
    # Cela créera les dossiers TD_... dans le répertoire des bams.
    td_dirs_ctrl  <- make_tag_dirs_custom(set1, ctrl_indices, bam_dir, bam_dir, opt$homer_bin)
    td_dirs_treat <- make_tag_dirs_custom(set2, treat_indices, bam_dir, bam_dir, opt$homer_bin)
    
    td_dirs <- c(td_dirs_ctrl, td_dirs_treat)

    # 4. Annotation
    annot_out <- file.path(bg_dir, paste0(set1, "_", set2, "_annot.txt"))
    run_annotate_peak(
      opt$homer_bin, bam_dir,
      file.path(bg_dir, paste0("pooled_", set1, "_", set2, "_hotspot_modify.txt")),
      genome, td_dirs, annot_out
    )

    # 5. Analyse Différentielle
    # On crée les labels basés sur la longueur des listes d'indices
    group1_labels <- rep(set1, length(ctrl_indices))
    group2_labels <- rep(set2, length(treat_indices))
    
    cat("   Differential Analysis configuration:\n")
    cat("     - Background (", set1, "): ", length(group1_labels), " samples\n", sep="")
    cat("     - Target     (", set2, "): ", length(group2_labels), " samples\n", sep="")

    diff_out <- file.path(bg_dir, paste0(set1, "_", set2, "_output.txt"))
    run_diff_expression(opt$homer_bin, bam_dir, annot_out, diff_out, group1_labels, group2_labels)

    # 6. Extraction des pics & Motifs
    create_diff_peak_files(
      input_path = bg_dir,
      input_filename = paste0(set1, "_", set2, "_output.txt"),
      sep = "\t", chr = "Chr", start = "Start", end = "End",
      padj_cutoff = opt$peak_padj_cutoff,
      l2fc_cutoff = opt$peak_l2fc_cutoff,
      output_prefix = "my_peaks"
    )

    target_bed <- file.path(bg_dir, "my_peaks_simple.bed")
    # Gestion fallback fichier bed
    if (!file.exists(target_bed)) {
      misplaced_bed <- file.path(bam_dir, "my_peaks_simple.bed")
      if (file.exists(misplaced_bed)) file.rename(misplaced_bed, target_bed)
      else if (file.exists("my_peaks_simple.bed")) file.rename("my_peaks_simple.bed", target_bed)
    }

    if (!file.exists(target_bed)) {
        cat("   [Warning] No differential peaks found (or file missing).\n")
    } else {
        cat("   Running Motif Finding (findMotifsGenome)...\n")
        motif_file_arg <- file.path(mm_files_dir, "all_cluster_motifs.motifs")
        if (!file.exists(motif_file_arg)) {
            motif_file_arg <- system.file("data", "all_cluster_motifs.motifs", package="tiger")
        }
        run_homer(target_bed, homer_out_dir, genome, file.path(opt$homer_bin, "findMotifsGenome.pl"), paste("-mknown", motif_file_arg))
        homer_results_file <- known_results
    }
  }
}
# ============================ 6. TIGER FOOTPRINTING (ROBUST BATCH MODE) ========
cat("\n>>> Step 4: Running TIGER Footprinting (Batch size: 25)...\n")

# 1. Chargement des CutCounts (Une seule fois)
cat("   Loading CutCounts...\n")
assign(set1, CutCount(file = file.path(bg_dir, paste0(set1, "_AC_cutcount.bgr.gz")), count = 1, name = set1))
assign(set2, CutCount(file = file.path(bg_dir, paste0(set2, "_AC_cutcount.bgr.gz")), count = 1, name = set2))

# 2. Chargement des motifs
motiflistfile <- file.path(mm_files_dir, paste0("motiflist_", mm, "_cluster.txt"))
motifdb_dir   <- file.path(mm_files_dir, paste0(mm, "_cluster_fimo"))
motifdb_mm    <- MotifDB(motiflistfile = motiflistfile, directory = motifdb_dir)

if (!file.exists(motiflistfile)) stop("Error: Motif list missing")

motif_lines <- readLines(motiflistfile)
motif_lines <- motif_lines[nzchar(motif_lines)]
motif_no    <- if (is.null(opt$max_motifs)) length(motif_lines) else opt$max_motifs

cat(paste("   Total Motifs to process:", motif_no, "\n"))

# 3. Configuration fichiers
c <- paste0(set1, "-", set2)
bias_file_abs <- file.path(bg_dir, paste0("Hexamer_", set1, "_", set2, "_merged_withoutMap.csv"))
site_file_abs <- file.path(bg_dir, paste0("pooled_", set1, "_", set2, "_hotspot.csv"))

# Vérification stricte
if (!file.exists(bias_file_abs)) {
     bias_file_abs <- file.path(bg_dir, paste0("Hexamer_", set1, "_", set2, "_merged_withoutMap.csv"))
     if (!file.exists(bias_file_abs)) stop(paste("CRITICAL: Bias file not found at", bias_file_abs))
}
if (!file.exists(site_file_abs)) stop(paste("CRITICAL: Peak file not found at", site_file_abs))

gfootoption_couple <- GFootOption(
    biasfile = bias_file_abs,
    hexamer_pattern = paste0("nuccode_", mm, "_6mer_{chr}.dat")
)

# 4. BOUCLE INTELLIGENTE (CHECK & SKIP)
# On parcourt tout de 1 à la fin, et on saute si le dossier existe déjà.

for (r in seq(1, motif_no, by = 25)) {
    
    # Nom du dossier de sortie pour ce batch
    batch_output_dir <- file.path(bg_dir, paste0("OUTPUT_", c, "_", r))
    
    # --- LA CORRECTION EST ICI ---
    # Si le dossier existe et contient des fichiers, on considère que c'est fait.
    if (dir.exists(batch_output_dir) && length(list.files(batch_output_dir)) > 0) {
        cat(paste("   [SKIP] Batch starting at", r, "already exists in", basename(batch_output_dir), "\n"))
        next # On passe directement au suivant !
    }
    
    cat(paste("\n   >>> Processing batch starting at motif:", r, "\n"))
    
    GFoot_obj <- GFoot(
      control   = get(set1),
      treatment = get(set2),
      sitefile  = site_file_abs,
      motifDB   = motifdb_mm,
      gfootoption = gfootoption_couple,
      outputdir = paste0("OUTPUT_", c, "_", r),
      cachedir  = paste0("CACHE_", c, "_", r)
    )

    tryCatch({
        run(
          GFoot_obj,
          graphout = TRUE,
          yrange   = c(-2.2, 1.0),
          mc.cores = opt$cores,
          # Le range est bien défini ici
          range    = r:ifelse(r + 24 < motif_no, r + 24, motif_no), 
          run_set  = c
        )
    }, error = function(e) {
        cat(paste("\n[CRITICAL ERROR] Batch starting at", r, "failed!\n"))
        cat("Error message:", e$message, "\n")
        # On n'arrête pas tout le script, on essaie de continuer au batch suivant
        # ou on stop si c'est trop grave. Ici on stop pour éviter de spammer les logs.
        stop("Stopping execution due to TIGER error.") 
    })
    
    gc() # Nettoyage mémoire
}

cat("\n>>> Footprinting Done.\n")
# ============================ 7. FINAL PLOT ===================================
cat("\n>>> Step 5: Generating Bagplot & Report...\n")
gc()
setwd(bg_dir)

merged_name <- paste0(
  set1, "_", set2, "_On_pooled_", set1, "_", set2,
  "_hotspot_footprint_depth_table_merged_80_rep2.csv"
)

# CORRECTION 3 : Fusion intelligente (ignore les headers répétés)
cat("   Merging CSV files...\n")
merge_cmd <- paste0(
  "find . -type f -name '",
  set1, "_", set2, "_On_pooled_", set1, "_", set2,
  "_hotspot_footprint_depth_table.csv' -exec tail -n +2 -q {} ';' > ",
  merged_name
)
system(merge_cmd)

if (file.exists(merged_name) && file.info(merged_name)$size > 0) {
  
  cat("   Reading merged dataset (using fread)...\n")
  # CORRECTION 4 : fread pour éviter le crash mémoire
  dat <- data.table::fread(
    file.path(bg_dir, merged_name),
    header = FALSE,
    sep = " ", 
    quote = "\"",
    nThread = opt$cores,
    data.table = FALSE 
  )
  
  # Fallback si séparateur tabulation
  if (ncol(dat) < 16) {
      dat <- data.table::fread(file.path(bg_dir, merged_name), header=FALSE, sep="\t", quote="\"", data.table=FALSE)
  }

  colnames(dat) <- c(
    "num", "filename", "center", "logo", "motiflength", "motif",
    "memeNo", "memeEntry", "lrmotif", "lrflanking",
    "lrtotal", "fdepth1", "fdepth2", "meancc1", "meancc2", "numsite"
  )
  
  num_rows <- nrow(dat)
  cat(paste("   Dataset contains", num_rows, "rows.\n"))

  # Sécurité PDF
  generate_pdf <- TRUE
  if (num_rows > 100000) {
      cat("   [WARNING] Dataset >100k rows. PDF disabled to prevent crash.\n")
      generate_pdf <- FALSE
  }

  use_homer_plot <- run_homer_flag && !is.null(homer_results_file) && file.exists(homer_results_file)
  
  tryCatch({
	suppressWarnings(
      gen_bagplot_tiger(
        dat,
        dataname1 = set1,
        dataname2 = set2,
        qvaluethreshold_bagplot = 0.05,
        factor = 1.5,
        pdf = generate_pdf,
        html = TRUE,
        include_homer_results = use_homer_plot,
        homer_results_file = homer_results_file,
        pval_threshold_homer = 0.05
      )
	 )
      cat(paste0("\n>>> Success! Files generated in: ", bg_dir, "\n"))
  }, error = function(e) {
      cat("\n[ERROR] Plot generation failed. Error message:", e$message, "\n")
  })
  
} else {
  cat("\n[Warning] Merged CSV file is empty or missing.\n")
}

if(!is.null(dev.list())) dev.off()
gc()

if (file.exists("Rplots.pdf")) {
    file.remove("Rplots.pdf")
}