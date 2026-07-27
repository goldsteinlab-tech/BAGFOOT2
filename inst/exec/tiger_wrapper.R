#!/usr/bin/env Rscript

.libPaths(R.home("library"))

# ==============================================================================
# BAGFOOT2 WRAPPER - UNEQUAL REPLICATES & ROBUST HANDLING
# ==============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(bagfoot2)
  library(hash)
  library(data.table)
  library(Cairo)
  library(aplpack)
  library(stringr)
  library(dplyr)
  library(ggplot2)
  library(plotly)
  library(patchwork)
  library(extRemes)
  library(robustbase)
  library(CerioliOutlierDetection)
})

# ============================ 1. OPTIONS DEFINITION ===========================
option_list <- list(

  make_option(c("-h", "--help"), action="store_true", default=FALSE, help="Show this help message and exit"),

  # --- INPUT / OUTPUT ---
  make_option(c("-b", "--bam_dir"), type="character", help="Directory containing merged BAM files and the replicas if Homer is used"),
  make_option(c("-p", "--peak_dir"), type="character", help="Directory containing peak files"),
  make_option(c("-o", "--output_dir"), type="character", help="Output directory (will be created if not exists)"),
  make_option(c("-m", "--mm_files_dir"), type="character", default=NULL,
            help="Directory containing motif files. Optional: if omitted, mm10/hg38 files are downloaded automatically from Zenodo."),

  # --- CONDITIONS & GENOME ---
  make_option(c("-c", "--control"), type="character", help="Control condition name (e.g., fed24h)"),
  make_option(c("-t", "--treatment"), type="character", help="Treatment condition name (e.g., fasted24h)"),
  make_option(c("-g", "--genome"), type="character", default="mm10", help="Genome version [default: mm10]"),
  make_option(c("--project"), type="character", default="bagfoot2_proj", help="Project name directory containing bam_merged and macs2 folders"),

  # --- FILE SUFFIXES ---
  make_option(c("--peak_suffix"), type="character", default="_peaks.narrowPeak", help="Suffix for peak files"),
  make_option(c("--bam_suffix"), type="character", default="_merged.bam", help="Suffix for BAM files"),

  # --- BAGFOOT2 PARAMETERS ---
  # Hidden: max number of motifs to process. If not set, a genome-specific default is used.
  make_option(c("--max_motifs"), type="integer", default=NULL, help="Max number of motifs to process. If not set, ALL motifs are processed."),

  # Hidden: also generate the bagplot figure itself (off by default)
  make_option(c("--bagplot"), type="character", default="FALSE",
              help="Also generate the bagplot figure (TRUE/FALSE) [default: FALSE]"),

  # --- MAHALANOBIS BARPLOT (robust MCD D2 ranking) ---
  make_option(c("--maha_barplot"), type="character", default="TRUE",
              help="Generate the Mahalanobis (robust MCD) barplots (TRUE/FALSE) [default: TRUE]"),
  make_option(c("--maha_top_n"), type="integer", default=50,
              help="Max bars per panel in the Mahalanobis barplots [default: %default]"),
  make_option(c("--maha_mode"), type="character", default="both",
              help="Mahalanobis barplot layout: 'padj' (bar = -log10(FDR), fill = D2), 'D2' (bar = D2, fill = FDR), or 'both' [default: %default]"),
  make_option(c("--maha_qvalue"), type="double", default=0.05,
              help="FDR threshold for the Mahalanobis barplots [default: %default]"),
  make_option(c("--mcd_alpha"), type="double", default=0.75,
              help="MCD 'clean subset' fraction for the robust Mahalanobis distance. Lower towards 0.5 (max breakdown) if many motifs are expected to be true positives [default: %default]"),

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
opt_parser <- OptionParser(option_list=option_list, usage = "bagfoot2 [options]", add_help_option = FALSE)
opt <- parse_args(opt_parser, print_help_and_exit = FALSE)

# ============================ 2. SETUP & CLEANING =============================
# --- personnalized help ---
if (opt$help) {
  cat("BAGFOOT2: a computational tool to predict transcription factor activity from ATAC-seq\n")
  cat("\n")
  cat("Version:   2.0.0\n")
  cat("About:     developed in the Goldstein Lab.\n")
  cat("Code:      https://github.com/goldsteinlab-tech/BAGFOOT2\n")
  cat("\n")

  # --- hide secret options (not shown in --help) ---
  hidden_flags <- c("--max_motifs", "--bagplot", "--mcd_alpha")
  visible_options <- option_list[ !sapply(option_list, function(x) any(x@long_flag %in% hidden_flags)) ]

  dummy_parser <- OptionParser(option_list=visible_options, add_help_option=FALSE)
  print_help(dummy_parser)
  quit(status = 0)
}

# Create Output Directory
if (!is.null(opt$output_dir)) {
  dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)
  bg_dir <- normalizePath(opt$output_dir)
  setwd(bg_dir)
  bg_dir <- paste0(bg_dir, "/")   # add separator
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
  nums <- as.integer(unlist(strsplit(gsub(" ", "", as.character(x)), ",")))
  return(nums)
}

run_homer_flag   <- to_logical(opt$homer)
run_bagplot_flag <- to_logical(opt$bagplot)
run_maha_flag    <- to_logical(opt$maha_barplot)
use_map <- ""

# --- Validate --maha_mode early, before hours of footprinting ---
if (!opt$maha_mode %in% c("padj", "D2", "both")) {
  stop(paste0("\n[BAGFOOT2 ERROR] --maha_mode must be one of 'padj', 'D2', 'both' (got: '",
              opt$maha_mode, "')\n"), call. = FALSE)
}

# System Options (SILENT MODE)
options("bagfoot2.mc.cores" = opt$cores)
options("bagfoot2.mc.cores.motif" = opt$cores)
options("bagfoot2.debug" = FALSE)

# Global Variables assignment
bam_dir      <- normalizePath(opt$bam_dir)
peak_dir     <- normalizePath(opt$peak_dir)
genome <- opt$genome

# Motif files directory: utilise --mm_files_dir si fourni, sinon
# resolution/telechargement automatique (env var > cache > Zenodo)
if (!is.null(opt$mm_files_dir)) {
  mm_files_dir <- normalizePath(opt$mm_files_dir, mustWork = FALSE)
} else {
  cat("\n>>> --mm_files_dir non fourni, resolution automatique des fichiers motifs...\n")
  mm_files_dir <- ensure_motif_db(genome)
}
set1 <- opt$control
set2 <- opt$treatment

# Explicitly define 'mm' in Global Environment for BAGFOOT2 internals
mm <- genome
assign("mm", genome, envir = .GlobalEnv)

# ============================ PRE-FLIGHT CHECKS ===============================
cat("\n>>> Pre-flight checks...\n")

# Clear, prefixed fatal error (no R stack trace noise)
fail <- function(msg) stop(paste0("\n[BAGFOOT2 ERROR] ", msg, "\n"), call. = FALSE)

# --- Genome supported? ---
supported_genomes <- c("mm9", "mm10", "hg18", "hg19", "hg38")
if (!genome %in% supported_genomes) {
  fail(paste0("Genome '", genome, "' is not supported. ",
              "Choose one of: ", paste(supported_genomes, collapse = ", "), "."))
}

# --- BSgenome package for this genome installed? ---
bsgenome_pkg <- switch(genome,
                       "mm9"  = "BSgenome.Mmusculus.UCSC.mm9",
                       "mm10" = "BSgenome.Mmusculus.UCSC.mm10",
                       "hg18" = "BSgenome.Hsapiens.UCSC.hg18",
                       "hg19" = "BSgenome.Hsapiens.UCSC.hg19",
                       "hg38" = "BSgenome.Hsapiens.UCSC.hg38")
if (!requireNamespace(bsgenome_pkg, quietly = TRUE)) {
  fail(paste0("Reference genome package '", bsgenome_pkg, "' is required for '", genome,
              "' but is not installed in this conda environment.\n",
              "         Install it with:\n",
              "         Rscript -e 'BiocManager::install(\"", bsgenome_pkg, "\")'"))
}

# --- Input directories exist? ---
if (!dir.exists(bam_dir))      fail(paste0("BAM directory not found: ", bam_dir))
if (!dir.exists(peak_dir))     fail(paste0("Peak directory not found: ", peak_dir))
if (!dir.exists(mm_files_dir)) fail(paste0("Motif files directory (--mm_files_dir) not found: ", mm_files_dir))

# --- Motif resources present for this genome? ---
req_motiflist <- file.path(mm_files_dir, paste0("motiflist_", mm, "_cluster.txt"))
req_fimo_dir  <- file.path(mm_files_dir, paste0(mm, "_cluster_fimo"))
if (!file.exists(req_motiflist))
  fail(paste0("Motif list not found: ", req_motiflist,
              "\n         Make sure --mm_files_dir contains the ", mm, " support files."))
if (!dir.exists(req_fimo_dir))
  fail(paste0("FIMO motif directory not found: ", req_fimo_dir,
              "\n         Make sure --mm_files_dir contains the ", mm, " support files."))

# --- samtools available (needed for BAM indexing)? ---
if (Sys.which("samtools") == "")
  fail("'samtools' not found in PATH. Activate the conda environment that provides samtools.")

cat("   All pre-flight checks passed.\n")


# ============================ 3. BAM & PEAK PROCESSING ========================
cat("\n>>> Step 1: Preparing BAM and Peak files...\n")

all_set <- c(set1, set2)

for (treat in all_set) {
  tryCatch({
    # BAM Processing using --bam_suffix
    bamfile <- file.path(bam_dir, paste0(treat, opt$bam_suffix))
    if (!file.exists(bamfile)) stop(paste("Error: Missing BAM file:", bamfile))

    cat(paste("   Processing:", treat, "\n"))

    # BAGFOOT2 CutCount
    assign(paste0("cc_", treat), countReadsBAM(bamfile, mm))

    # Use 'treat' as name, not the BAM filename
    cutcountfile <- makeCutCountBAMWithName(bamfile, name = treat, refgenome = mm)

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

  }, error = function(e) {
    fail(paste0("Failed while processing condition '", treat, "'.\n",
                "         R message: ", e$message, "\n",
                "         Check that ", file.path(bam_dir, paste0(treat, opt$bam_suffix)),
                " exists, is indexed (.bai), and uses 'chr'-prefixed chromosomes."))
  })
}

# Clean temp files
temp_f <- list.files(bg_dir, pattern = glob2rx("temp*.dat"), full.names = TRUE)
if(length(temp_f) > 0) file.remove(temp_f)

# ============================ 4. MERGING COUPLES ==============================
cat("\n>>> Step 2: Merging couples...\n")

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

    # --- HOMER pre-flight ---
    if (is.null(opt$homer_bin))      fail("--homer_bin is required when --homer is TRUE.")
    if (!dir.exists(opt$homer_bin))  fail(paste0("HOMER bin directory not found: ", opt$homer_bin))

    for (bin in c("makeTagDirectory", "annotatePeaks.pl",
                  "getDiffExpression.pl", "findMotifsGenome.pl")) {
      if (!file.exists(file.path(opt$homer_bin, bin)))
        fail(paste0("HOMER binary '", bin, "' not found in ", opt$homer_bin,
                    ". Check your --homer_bin path."))
    }

    homer_root <- dirname(normalizePath(opt$homer_bin))
    genome_dir <- file.path(homer_root, "data", "genomes", genome)
    cfg_file   <- file.path(homer_root, "config.txt")
    genome_installed <- dir.exists(genome_dir) ||
      (file.exists(cfg_file) &&
         any(grepl(paste0("^", genome, "\\b"), readLines(cfg_file, warn = FALSE))))
    if (!genome_installed) {
      fail(paste0("Genome '", genome, "' is not installed in HOMER.\n",
                  "         Install it with:\n",
                  "         perl ", file.path(homer_root, "configureHomer.pl"),
                  " -install ", genome))
    }

    ctrl_indices  <- parse_reps(opt$reps_control)
    treat_indices <- parse_reps(opt$reps_treatment)

    cat(paste("   Control Indices:  ", paste(ctrl_indices, collapse=","), "\n"))
    cat(paste("   Treatment Indices:", paste(treat_indices, collapse=","), "\n"))

    process_hotspot_file(
      pooled_csv,
      file.path(bg_dir, paste0("pooled_", set1, "_", set2, "_hotspot_modify.txt"))
    )

    cat("   Creating Tag Directories in BAM folder (Pattern: TD_Condition_Index)...\n")
    make_tag_dirs_custom <- function(condition, indices, bam_dir_path, out_dir_path, homer_path) {
      dirs_created <- c()
      for (i in indices) {
        bam_name <- paste0(condition, "_", i, ".bam")
        bam_path <- file.path(bam_dir_path, bam_name)
        tag_dir_name <- paste0("TD_", condition, "_", i)
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

    td_dirs_ctrl  <- make_tag_dirs_custom(set1, ctrl_indices, bam_dir, bam_dir, opt$homer_bin)
    td_dirs_treat <- make_tag_dirs_custom(set2, treat_indices, bam_dir, bam_dir, opt$homer_bin)
    td_dirs <- c(td_dirs_ctrl, td_dirs_treat)

    annot_out <- file.path(bg_dir, paste0(set1, "_", set2, "_annot.txt"))
    run_annotate_peak(
      opt$homer_bin, bam_dir,
      file.path(bg_dir, paste0("pooled_", set1, "_", set2, "_hotspot_modify.txt")),
      genome, td_dirs, annot_out
    )

    group1_labels <- rep(set1, length(ctrl_indices))
    group2_labels <- rep(set2, length(treat_indices))

    cat("   Differential Analysis configuration:\n")
    cat("     - Background (", set1, "): ", length(group1_labels), " samples\n", sep="")
    cat("     - Target     (", set2, "): ", length(group2_labels), " samples\n", sep="")

    diff_out <- file.path(bg_dir, paste0(set1, "_", set2, "_output.txt"))
    run_diff_expression(opt$homer_bin, bam_dir, annot_out, diff_out, group1_labels, group2_labels)

    # --- Extraction des pics & Motifs ---
    # If there are NO differential peaks: warn, continue, and keep HOMER OFF for the rest.
    no_diff_peaks <- FALSE
    tryCatch({
      create_diff_peak_files(
        input_path = bg_dir,
        input_filename = paste0(set1, "_", set2, "_output.txt"),
        sep = "\t", chr = "Chr", start = "Start", end = "End",
        padj_cutoff = opt$peak_padj_cutoff,
        l2fc_cutoff = opt$peak_l2fc_cutoff,
        output_prefix = "my_peaks"
      )
    }, error = function(e) {
      if (grepl("No differential peak", e$message, ignore.case = TRUE)) {
        cat(sprintf(
          "\n[WARNING] No differential peaks found (padj <= %s, |log2FC| >= %s).\n",
          opt$peak_padj_cutoff, opt$peak_l2fc_cutoff))
        cat("          Continuing WITHOUT HOMER results. HOMER will be turned OFF for the rest.\n")
        no_diff_peaks <<- TRUE
      } else {
        fail(paste0("create_diff_peak_files failed: ", e$message))
      }
    })

    if (no_diff_peaks) {
      homer_results_file <- NULL   # HOMER off downstream
    } else {
      target_bed <- file.path(bg_dir, "my_peaks_simple.bed")
      if (!file.exists(target_bed)) {
        misplaced_bed <- file.path(bam_dir, "my_peaks_simple.bed")
        if (file.exists(misplaced_bed)) file.rename(misplaced_bed, target_bed)
        else if (file.exists("my_peaks_simple.bed")) file.rename("my_peaks_simple.bed", target_bed)
      }

      if (!file.exists(target_bed)) {
        cat("   [WARNING] Differential peak file missing. Continuing WITHOUT HOMER results.\n")
        homer_results_file <- NULL
      } else {
        cat("   Running Motif Finding (findMotifsGenome)...\n")
        motif_file_arg <- file.path(mm_files_dir, "all_cluster_motifs.motifs")
        if (!file.exists(motif_file_arg)) {
          motif_file_arg <- system.file("data", "all_cluster_motifs.motifs", package = "bagfoot2")
        }
        run_homer(target_bed, homer_out_dir, genome,
                  file.path(opt$homer_bin, "findMotifsGenome.pl"),
                  paste("-mknown", motif_file_arg, "-p ",opt$cores))
        homer_results_file <- known_results
      }
    }
  }
}

# ============================ 6. BAGFOOT2 FOOTPRINTING (ROBUST BATCH MODE) =====
cat("\n>>> Step 4: Running BAGFOOT2 Footprinting (Batch size: 25)...\n")

cat("   Loading CutCounts...\n")
assign(set1, CutCount(file = file.path(bg_dir, paste0(set1, "_AC_cutcount.bgr.gz")), count = 1, name = set1))
assign(set2, CutCount(file = file.path(bg_dir, paste0(set2, "_AC_cutcount.bgr.gz")), count = 1, name = set2))

motiflistfile <- file.path(mm_files_dir, paste0("motiflist_", mm, "_cluster.txt"))
motifdb_dir   <- file.path(mm_files_dir, paste0(mm, "_cluster_fimo"))
motifdb_mm    <- MotifDB(motiflistfile = motiflistfile, directory = motifdb_dir)

if (!file.exists(motiflistfile)) stop("Error: Motif list missing")

motif_lines <- readLines(motiflistfile)
motif_lines <- motif_lines[nzchar(motif_lines)]

# --- Number of motifs to process ---
# Priority: explicit --max_motifs > genome-specific default > all motifs in the list.
if (!is.null(opt$max_motifs)) {
  motif_no <- opt$max_motifs
} else {
  motif_no <- switch(mm,
                     "mm10" = 577,
                     "hg38" = 636,
                     length(motif_lines))   # default for other genomes: all
}
# Safety: never exceed the number of available motifs
motif_no <- min(motif_no, length(motif_lines))

cat(paste("   Total Motifs to process:", motif_no, "\n"))

c <- paste0(set1, "-", set2)
bias_file_abs <- file.path(bg_dir, paste0("Hexamer_", set1, "_", set2, "_merged_withoutMap.csv"))
site_file_abs <- file.path(bg_dir, paste0("pooled_", set1, "_", set2, "_hotspot.csv"))

if (!file.exists(bias_file_abs)) stop(paste("CRITICAL: Bias file not found at", bias_file_abs))
if (!file.exists(site_file_abs)) stop(paste("CRITICAL: Peak file not found at", site_file_abs))

gfootoption_couple <- GFootOption(
    biasfile = bias_file_abs,
    hexamer_pattern = paste0("nuccode_", mm, "_6mer_{chr}.dat")
)

for (r in seq(1, motif_no, by = 25)) {

    batch_output_dir <- file.path(bg_dir, paste0("OUTPUT_", c, "_", r))

    if (dir.exists(batch_output_dir) && length(list.files(batch_output_dir)) > 0) {
        cat(paste("   [SKIP] Batch starting at", r, "already exists in", basename(batch_output_dir), "\n"))
        next
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
        invisible(capture.output(suppressMessages(suppressWarnings(
          run(
            GFoot_obj,
            graphout = TRUE,
            yrange   = c(-2.2, 1.0),
            mc.cores = opt$cores,
            range    = r:ifelse(r + 24 < motif_no, r + 24, motif_no),
            run_set  = c
          )
        ))))
    }, error = function(e) {
      fail(paste0("BAGFOOT2 footprinting failed on batch starting at motif ", r, ".\n",
                  "         R message: ", e$message, "\n",
                  "         Common causes: missing nuccode_", mm,
                  "_6mer_{chr}.dat files in --mm_files_dir, or a BAM/peak chromosome",
                  " naming mismatch (BAGFOOT2 needs 'chr'-prefixed names)."))
    })

    gc()
}

cat("\n>>> Footprinting Done.\n")

# ============================ 7. FINAL PLOTS ==================================
cat("\n>>> Step 5: Generating plots & Report...\n")
gc()
setwd(bg_dir)

merged_name <- paste0(
  set1, "_", set2, "_On_pooled_", set1, "_", set2,
  "_hotspot_footprint_depth_table_merged_80_rep2.csv"
)

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

  dat <- data.table::fread(
    file.path(bg_dir, merged_name),
    header = FALSE,
    sep = " ",
    quote = "\"",
    nThread = opt$cores,
    data.table = FALSE
  )

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

  # PDF safety
  generate_pdf <- TRUE
  if (num_rows > 100000) {
      cat("   [WARNING] Dataset >100k rows. PDF disabled to prevent crash.\n")
      generate_pdf <- FALSE
  }

  # HOMER plot only if HOMER ran AND produced results (auto FALSE if no diff peaks)
  use_homer_plot <- run_homer_flag && !is.null(homer_results_file) && file.exists(homer_results_file)

  tryCatch({
    # --- EVT barplot : MAIN output, always generated ---
    cat("   Generating EVT barplots...\n")
    suppressWarnings(
      gen_evt_barplot_bagfoot2(dat,
                               dataname1 = set1,
                               dataname2 = set2,
                               pval_threshold = 0.05,
                               include_homer_results = use_homer_plot,
                               homer_results_file = homer_results_file,
                               pval_threshold_homer = 0.05,
                               top_n = 50)
    )

    # --- Mahalanobis barplots (robust MCD D2) + bagplot CSV ---------------
    #
    # NOTE ON PANEL ORDER: lrtotal = log2(set2 / set1), so a positive delta FA
    # means higher accessibility in set2 (the treatment). Inside
    # gen_bagplot_bagfoot2() the "increase" panel is titled with dataname1,
    # therefore dataname1 must be set2 and dataname2 must be set1 here.
    # This is the same (inverted) order used in the .Rmd, and it is the
    # OPPOSITE of the EVT call above, where the increase panel is titled
    # with dataname2.
    #
    # gen_bagplot_bagfoot2() is what produces the Mahalanobis barplots, so it
    # must be called whenever --maha_barplot is TRUE, even if --bagplot FALSE.
    if (run_maha_flag || run_bagplot_flag) {
      cat("   Generating Mahalanobis barplots (robust MCD)...\n")
      suppressWarnings(
        gen_bagplot_bagfoot2(
          dat,
          dataname1 = set2,
          dataname2 = set1,
          qvaluethreshold_bagplot = opt$maha_qvalue,
          factor = 1.5,
          pdf = generate_pdf && run_bagplot_flag,
          html = FALSE,
          include_homer_results = use_homer_plot,
          homer_results_file = homer_results_file,
          pval_threshold_homer = 0.05,
          mcd_alpha    = opt$mcd_alpha,
          maha_barplot = run_maha_flag,
          maha_top_n   = opt$maha_top_n,
          maha_mode    = opt$maha_mode
        )
      )
    }

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