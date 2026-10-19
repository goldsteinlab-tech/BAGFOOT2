#!/usr/bin/env Rscript
# ==============================================================================
# BAGFOOT2 - command-line wrapper
# Differential transcription factor footprinting from ATAC-seq
# https://github.com/goldsteinlab-tech/BAGFOOT2
# ==============================================================================

.libPaths(R.home("library"))

suppressPackageStartupMessages({
  library(optparse)
  library(bagfoot2)
  library(dplyr)
  library(ggplot2)
})

# ============================ 1. OPTIONS ======================================
option_list <- list(

  make_option(c("-h", "--help"), action = "store_true", default = FALSE,
              help = "Show this help message and exit"),

  # --- INPUT / OUTPUT ---
  make_option(c("-b", "--bam_dir"), type = "character",
              help = "Directory containing the merged BAM files (and the replicate BAMs if HOMER is used)"),
  make_option(c("-p", "--peak_dir"), type = "character",
              help = "Directory containing the peak files"),
  make_option(c("-o", "--output_dir"), type = "character",
              help = "Output directory (created if it does not exist)"),
  make_option(c("-m", "--mm_files_dir"), type = "character", default = NULL,
              help = "Directory containing the motif files. Optional: if omitted, mm10/hg38 files are downloaded automatically from Zenodo."),

  # --- CONDITIONS & GENOME ---
  make_option(c("-c", "--control"), type = "character",
              help = "Control condition name (e.g. fed)"),
  make_option(c("-t", "--treatment"), type = "character",
              help = "Treatment condition name (e.g. fasted)"),
  make_option(c("-g", "--genome"), type = "character", default = "mm10",
              help = "Genome version [default: %default]"),

  # --- FILE SUFFIXES ---
  make_option(c("--peak_suffix"), type = "character", default = "_peaks.narrowPeak",
              help = "Suffix of the peak files [default: %default]"),
  make_option(c("--bam_suffix"), type = "character", default = "_merged.bam",
              help = "Suffix of the merged BAM files [default: %default]"),

  # --- MAHALANOBIS BARPLOTS ---
  make_option(c("--maha_barplot"), type = "character", default = "TRUE",
              help = "Generate the Mahalanobis (D2 / FDR) barplots (TRUE/FALSE) [default: %default]"),
  make_option(c("--maha_top_n"), type = "integer", default = 20,
              help = "Max bars per panel in the barplots [default: %default]"),
  make_option(c("--maha_mode"), type = "character", default = "both",
              help = "'padj' (bar = -log10(FDR), fill = D2), 'D2' (bar = D2, fill = FDR) or 'both' [default: %default]"),
  make_option(c("--maha_qvalue"), type = "double", default = 0.05,
              help = "FDR threshold [default: %default]"),

  make_option(c("--cores"), type = "integer", default = 2,
              help = "Number of CPU cores [default: %default]"),

  # --- HOMER ---
  make_option(c("--homer"), type = "character", default = "FALSE",
              help = "Run HOMER differential peaks + motif enrichment (TRUE/FALSE) [default: %default]"),
  make_option(c("--homer_bin"), type = "character", default = NULL,
              help = "Path to the HOMER bin folder (required if --homer TRUE)"),
  make_option(c("--reps_control"), type = "character", default = "1",
              help = "Replicate indices for the control, e.g. '1,2,3' [default: %default]"),
  make_option(c("--reps_treatment"), type = "character", default = "1",
              help = "Replicate indices for the treatment, e.g. '1,2' [default: %default]"),
  make_option(c("--peak_l2fc_cutoff"), type = "double", default = 0.5,
              help = "|log2FC| threshold for differential peaks [default: %default]"),
  make_option(c("--peak_padj_cutoff"), type = "double", default = 0.05,
              help = "FDR threshold for differential peaks [default: %default]")
)

opt_parser <- OptionParser(option_list = option_list,
                           usage = "bagfoot2 [options]",
                           add_help_option = FALSE)
opt <- parse_args(opt_parser, print_help_and_exit = FALSE)

if (opt$help) {
  cat("BAGFOOT2: differential transcription factor footprinting from ATAC-seq\n\n")
  cat("Version:   2.0.0\n")
  cat("About:     developed in the Goldstein Lab.\n")
  cat("Code:      https://github.com/goldsteinlab-tech/BAGFOOT2\n\n")
  print_help(opt_parser)
  quit(status = 0)
}

# ============================ 2. SETUP ========================================

# Clear, prefixed fatal error (no R stack trace)
fail <- function(msg) stop(paste0("\n[BAGFOOT2 ERROR] ", msg, "\n"), call. = FALSE)

to_logical <- function(x) {
  if (is.null(x)) return(FALSE)
  toupper(as.character(x)) %in% c("TRUE", "T", "YES", "1")
}

parse_reps <- function(x) {
  if (is.null(x) || x == "") return(1L)
  as.integer(unlist(strsplit(gsub(" ", "", as.character(x)), ",")))
}

# --- Required arguments ---
for (arg in c("bam_dir", "peak_dir", "output_dir", "control", "treatment")) {
  if (is.null(opt[[arg]])) fail(paste0("--", arg, " is required. Run 'bagfoot2 --help'."))
}
if (opt$control == opt$treatment) fail("--control and --treatment must be different.")
if (!opt$maha_mode %in% c("padj", "D2", "both"))
  fail(paste0("--maha_mode must be one of 'padj', 'D2', 'both' (got: '", opt$maha_mode, "')."))
if (is.na(opt$cores) || opt$cores < 1) fail("--cores must be a positive integer.")

run_homer_flag <- to_logical(opt$homer)
run_maha_flag  <- to_logical(opt$maha_barplot)

# --- Parallelism: --cores footprinting workers, motifs processed one at a time
#     (MCCORES_MOTIF = 1 avoids forking the graphics device / fontconfig).
assignInNamespace("MCCORES",       as.integer(opt$cores), ns = "bagfoot2")
assignInNamespace("MCCORES_MOTIF", 1L,                    ns = "bagfoot2")
assignInNamespace("BFOOT_DEBUG",   FALSE,                 ns = "bagfoot2")

# --- Output directory ---
dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)
bg_dir <- paste0(normalizePath(opt$output_dir), "/")
setwd(bg_dir)

# --- Global variables ---
bam_dir  <- normalizePath(opt$bam_dir,  mustWork = FALSE)
peak_dir <- normalizePath(opt$peak_dir, mustWork = FALSE)
genome   <- opt$genome
set1     <- opt$control
set2     <- opt$treatment

# BAGFOOT2 internals expect 'mm' in the global environment
mm <- genome
assign("mm", genome, envir = .GlobalEnv)

# ============================ 3. PRE-FLIGHT CHECKS ============================
cat("\n>>> Pre-flight checks...\n")

supported_genomes <- c("mm9", "mm10", "hg18", "hg19", "hg38")
if (!genome %in% supported_genomes) {
  fail(paste0("Genome '", genome, "' is not supported. Choose one of: ",
              paste(supported_genomes, collapse = ", "), "."))
}

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

if (!dir.exists(bam_dir))  fail(paste0("BAM directory not found: ", bam_dir))
if (!dir.exists(peak_dir)) fail(paste0("Peak directory not found: ", peak_dir))

if (Sys.which("samtools") == "")
  fail("'samtools' not found in PATH. Activate the conda environment that provides samtools.")

# --- Motif files: --mm_files_dir, otherwise cache / automatic Zenodo download ---
if (!is.null(opt$mm_files_dir)) {
  mm_files_dir <- normalizePath(opt$mm_files_dir, mustWork = FALSE)
} else {
  cat("   --mm_files_dir not given: looking for motif files in the cache (download if needed)...\n")
  mm_files_dir <- ensure_motif_db(genome)
}

motiflistfile <- file.path(mm_files_dir, paste0("motiflist_", mm, "_cluster.txt"))
motifdb_dir   <- file.path(mm_files_dir, paste0(mm, "_cluster_fimo"))

if (!dir.exists(mm_files_dir))
  fail(paste0("Motif files directory not found: ", mm_files_dir))
if (!file.exists(motiflistfile))
  fail(paste0("Motif list not found: ", motiflistfile,
              "\n         Make sure --mm_files_dir contains the ", mm, " support files."))
if (!dir.exists(motifdb_dir))
  fail(paste0("FIMO motif directory not found: ", motifdb_dir,
              "\n         Make sure --mm_files_dir contains the ", mm, " support files."))

cat("   All pre-flight checks passed.\n")

# ============================ 4. BAM & PEAK PROCESSING ========================
cat("\n>>> Step 1: Preparing BAM and peak files...\n")

for (treat in c(set1, set2)) {
  tryCatch({
    bamfile <- file.path(bam_dir, paste0(treat, opt$bam_suffix))
    if (!file.exists(bamfile)) stop(paste("Missing BAM file:", bamfile))

    cat(paste("   Processing:", treat, "\n"))

    countReadsBAM(bamfile, mm)
    makeCutCountBAMWithName(bamfile, name = treat, refgenome = mm)

    invisible(capture.output({
      MakeBiasCorrectionTableBAM(
        bamfile   = bamfile,
        outfile   = paste0("Hexamer_", treat, "_", mm, "_withoutMap.txt"),
        refgenome = mm,
        np        = 6,
        atac      = TRUE,
        mapdir    = ""
      )
    }))

    peak_file_in <- file.path(peak_dir, paste0(treat, opt$peak_suffix))
    if (!file.exists(peak_file_in)) stop(paste("Peak file not found:", peak_file_in))

    macs <- read.table(peak_file_in, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
    if (ncol(macs) < 3) stop(paste("Invalid peak file (less than 3 columns):", peak_file_in))
    names(macs)[1:3] <- c("chr", "st", "ed")

    nr <- nrow(macs)
    hotspot <- data.frame(
      ID = 1:nr, chrom. = macs$chr, start = macs$st + 1, end = macs$ed,
      MaxD = 0, AveD = 0, Zscore = 0, pvalue = 0
    )
    write.table(hotspot, file = file.path(bg_dir, paste0(treat, "_peaks.csv")),
                sep = ",", row.names = FALSE)

  }, error = function(e) {
    fail(paste0("Failed while processing condition '", treat, "'.\n",
                "         R message: ", e$message, "\n",
                "         Check that ", file.path(bam_dir, paste0(treat, opt$bam_suffix)),
                " exists, is indexed (.bai), and uses 'chr'-prefixed chromosomes."))
  })
}

temp_f <- list.files(bg_dir, pattern = glob2rx("temp*.dat"), full.names = TRUE)
if (length(temp_f) > 0) file.remove(temp_f)

# ============================ 5. MERGING THE PAIR =============================
cat("\n>>> Step 2: Merging control and treatment...\n")

combineTwoHotspots(
  file.path(bg_dir, paste0(set1, "_peaks.csv")),
  file.path(bg_dir, paste0(set2, "_peaks.csv")),
  set1, set2
)

pooled_csv <- file.path(bg_dir, paste0("pooled_", set1, "_", set2, "_hotspot.csv"))
pooled <- read.csv(pooled_csv)
if ("ID.1" %in% colnames(pooled)) {
  pooled <- pooled %>% dplyr::select(-ID.1)
  write.csv(pooled, pooled_csv, row.names = FALSE)
}

merge_two_frequency_tables(
  paste0("Hexamer_", set1, "_", mm, "_withoutMap.txt"),
  paste0("Hexamer_", set2, "_", mm, "_withoutMap.txt"),
  paste0("Hexamer_", set1, "_", set2, "_merged_withoutMap.csv")
)

# ============================ 6. HOMER (OPTIONAL) =============================
homer_results_file <- NULL

if (run_homer_flag) {
  homer_out_dir <- file.path(bg_dir, "homer_motif_res")
  known_results <- file.path(homer_out_dir, "knownResults.txt")

  if (file.exists(known_results)) {
    cat("\n>>> Step 3: HOMER results already present. Skipping HOMER.\n")
    homer_results_file <- known_results
  } else {
    cat("\n>>> Step 3: HOMER analysis...\n")

    # --- HOMER pre-flight ---
    if (is.null(opt$homer_bin))     fail("--homer_bin is required when --homer is TRUE.")
    if (!dir.exists(opt$homer_bin)) fail(paste0("HOMER bin directory not found: ", opt$homer_bin))

    for (bin in c("makeTagDirectory", "annotatePeaks.pl",
                  "getDiffExpression.pl", "findMotifsGenome.pl")) {
      if (!file.exists(file.path(opt$homer_bin, bin)))
        fail(paste0("HOMER binary '", bin, "' not found in ", opt$homer_bin, "."))
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
                  "         perl ", file.path(homer_root, "configureHomer.pl"), " -install ", genome))
    }

    ctrl_indices  <- parse_reps(opt$reps_control)
    treat_indices <- parse_reps(opt$reps_treatment)
    cat("   Control replicates:  ", paste(ctrl_indices,  collapse = ","), "\n")
    cat("   Treatment replicates:", paste(treat_indices, collapse = ","), "\n")

    peaks_homer <- file.path(bg_dir, paste0("pooled_", set1, "_", set2, "_hotspot_modify.txt"))
    process_hotspot_file(pooled_csv, peaks_homer)

    # One tag directory per replicate: <bam_dir>/TD_<condition>_<index>
    make_tag_dirs <- function(condition, indices) {
      dirs <- c()
      for (i in indices) {
        bam_path <- file.path(bam_dir, paste0(condition, "_", i, ".bam"))
        td_path  <- file.path(bam_dir, paste0("TD_", condition, "_", i))
        if (!file.exists(bam_path)) fail(paste("Replicate BAM not found:", bam_path))
        if (!dir.exists(td_path)) {
          cat("     Tag directory:", basename(td_path), "\n")
          system(paste(file.path(opt$homer_bin, "makeTagDirectory"), td_path, bam_path, "-format sam"),
                 ignore.stdout = TRUE, ignore.stderr = TRUE)
        } else {
          cat("     [Exists]", basename(td_path), "\n")
        }
        dirs <- c(dirs, td_path)
      }
      dirs
    }

    td_dirs <- c(make_tag_dirs(set1, ctrl_indices), make_tag_dirs(set2, treat_indices))

    annot_out <- file.path(bg_dir, paste0(set1, "_", set2, "_annot.txt"))
    diff_out  <- file.path(bg_dir, paste0(set1, "_", set2, "_output.txt"))

    run_annotate_peak(opt$homer_bin, bam_dir, peaks_homer, genome, td_dirs, annot_out)
    run_diff_expression(opt$homer_bin, bam_dir, annot_out, diff_out,
                        rep(set1, length(ctrl_indices)),
                        rep(set2, length(treat_indices)))

    # No differential peak -> warn and continue without HOMER
    no_diff_peaks <- FALSE
    tryCatch({
      create_diff_peak_files(
        input_path = bg_dir,
        input_filename = basename(diff_out),
        sep = "\t", chr = "Chr", start = "Start", end = "End",
        padj_cutoff = opt$peak_padj_cutoff,
        l2fc_cutoff = opt$peak_l2fc_cutoff,
        output_prefix = "my_peaks"
      )
    }, error = function(e) {
      if (grepl("No differential peak", e$message, ignore.case = TRUE)) {
        cat(sprintf("\n[WARNING] No differential peaks (padj <= %s, |log2FC| >= %s).\n",
                    opt$peak_padj_cutoff, opt$peak_l2fc_cutoff))
        cat("          Continuing WITHOUT HOMER results.\n")
        no_diff_peaks <<- TRUE
      } else {
        fail(paste0("create_diff_peak_files failed: ", e$message))
      }
    })

    if (!no_diff_peaks) {
      target_bed <- file.path(bg_dir, "my_peaks_simple.bed")
      if (!file.exists(target_bed)) {
        misplaced_bed <- file.path(bam_dir, "my_peaks_simple.bed")
        if (file.exists(misplaced_bed)) file.rename(misplaced_bed, target_bed)
      }

      if (!file.exists(target_bed)) {
        cat("   [WARNING] Differential peak file missing. Continuing WITHOUT HOMER results.\n")
      } else {
        cat("   Running motif enrichment (findMotifsGenome)...\n")
        motif_file_arg <- file.path(mm_files_dir, "all_cluster_motifs.motifs")
        if (!file.exists(motif_file_arg)) {
          motif_file_arg <- system.file("data", "all_cluster_motifs.motifs", package = "bagfoot2")
        }
        run_homer(target_bed, homer_out_dir, genome,
                  file.path(opt$homer_bin, "findMotifsGenome.pl"),
                  paste("-mknown", motif_file_arg, "-p", opt$cores))
        homer_results_file <- known_results
      }
    }
  }
}

# ============================ 7. FOOTPRINTING (BATCHES OF 25) =================
cat("\n>>> Step 4: Footprinting (batches of 25 motifs)...\n")

cat("   Loading cut counts...\n")
assign(set1, CutCount(file = file.path(bg_dir, paste0(set1, "_AC_cutcount.bgr.gz")),
                      count = 1, name = set1))
assign(set2, CutCount(file = file.path(bg_dir, paste0(set2, "_AC_cutcount.bgr.gz")),
                      count = 1, name = set2))

motifdb_mm <- MotifDB(motiflistfile = motiflistfile, directory = motifdb_dir)

motif_lines <- readLines(motiflistfile, warn = FALSE)
motif_lines <- motif_lines[nzchar(motif_lines)]

# Number of motifs: fixed per genome (mm10: 577, hg38: 636), all motifs otherwise
motif_no <- switch(mm,
                   "mm10" = 577,
                   "hg38" = 636,
                   length(motif_lines))
motif_no <- min(motif_no, length(motif_lines))
cat("   Motifs to process:", motif_no, "\n")

couple        <- paste0(set1, "-", set2)
bias_file_abs <- file.path(bg_dir, paste0("Hexamer_", set1, "_", set2, "_merged_withoutMap.csv"))
if (!file.exists(bias_file_abs)) fail(paste("Bias file not found:", bias_file_abs))
if (!file.exists(pooled_csv))    fail(paste("Pooled peak file not found:", pooled_csv))

gfootoption_couple <- GFootOption(
  biasfile        = bias_file_abs,
  hexamer_pattern = paste0("nuccode_", mm, "_6mer_{chr}.dat")
)

csv_pattern <- paste0(set1, "_", set2, "_On_pooled_", set1, "_", set2,
                      "_hotspot_footprint_depth_table.csv")
batch_starts <- seq(1, motif_no, by = 25)

for (r in batch_starts) {

  # A batch is complete only if its footprint table exists and is not empty
  batch_csv <- file.path(bg_dir, paste0("OUTPUT_", couple, "_", r), csv_pattern)
  if (file.exists(batch_csv) && file.info(batch_csv)$size > 0) {
    cat("   [SKIP] Batch", r, "already complete.\n")
    next
  }

  cat("   >>> Batch starting at motif", r, "(", format(Sys.time(), "%H:%M:%S"), ")\n")

  GFoot_obj <- GFoot(
    control     = get(set1),
    treatment   = get(set2),
    sitefile    = pooled_csv,
    motifDB     = motifdb_mm,
    gfootoption = gfootoption_couple,
    outputdir   = paste0("OUTPUT_", couple, "_", r),
    cachedir    = paste0("CACHE_", couple, "_", r)
  )

  tryCatch({
    run(GFoot_obj,
        graphout = TRUE,
        yrange   = c(-2.2, 1.0),
        mc.cores = opt$cores,
        range    = r:min(r + 24, motif_no),
        run_set  = couple)
  }, error = function(e) {
    fail(paste0("Footprinting failed on the batch starting at motif ", r, ".\n",
                "         R message: ", e$message, "\n",
                "         Common causes: missing nuccode_", mm, "_6mer_{chr}.dat files, ",
                "or a BAM/peak chromosome naming mismatch ('chr' prefix needed).\n",
                "         Delete OUTPUT_", couple, "_", r, " and CACHE_", couple, "_", r,
                " before relaunching."))
  })

  gc()
}

cat("\n>>> Footprinting done.\n")

# ============================ 8. STATISTICS & PLOTS ===========================
cat("\n>>> Step 5: Statistics and plots...\n")
setwd(bg_dir)

merged_name <- paste0(set1, "_", set2, "_On_pooled_", set1, "_", set2,
                      "_hotspot_footprint_depth_table_merged.csv")

# Check that every batch produced its table before merging
found_csv <- list.files(bg_dir, pattern = glob2rx(csv_pattern),
                        recursive = TRUE, full.names = TRUE)
cat(sprintf("   Batches expected: %d | tables found: %d\n",
            length(batch_starts), length(found_csv)))
if (length(found_csv) < length(batch_starts)) {
  cat("   [WARNING] Some batches produced no footprint table.\n")
  cat("             Delete the corresponding OUTPUT_* / CACHE_* folders and relaunch\n")
  cat("             before trusting the results.\n")
}

system(paste0("find . -type f -name '", csv_pattern,
              "' -exec tail -n +2 -q {} ';' > ", merged_name))

if (!file.exists(merged_name) || file.info(merged_name)$size == 0) {
  fail("Merged footprint table is empty or missing. No plot generated.")
}

dat <- data.table::fread(file.path(bg_dir, merged_name), header = FALSE, sep = " ",
                         quote = "\"", nThread = opt$cores, data.table = FALSE)
if (ncol(dat) < 16) {
  dat <- data.table::fread(file.path(bg_dir, merged_name), header = FALSE, sep = "\t",
                           quote = "\"", nThread = opt$cores, data.table = FALSE)
}
if (ncol(dat) != 16) {
  fail(paste0("Merged table has ", ncol(dat), " columns, expected 16. ",
              "The footprint tables are probably malformed."))
}

colnames(dat) <- c("num", "filename", "center", "logo", "motiflength", "motif",
                   "memeNo", "memeEntry", "lrmotif", "lrflanking", "lrtotal",
                   "fdepth1", "fdepth2", "meancc1", "meancc2", "numsite")
cat("   Motifs in the merged table:", nrow(dat), "\n")

use_homer_plot <- run_homer_flag && !is.null(homer_results_file) && file.exists(homer_results_file)

# Bagplot + D2 / FDR barplots + CSV table.
# lrtotal = log2(treatment / control): a positive delta FA means higher in the
# treatment. The "increase" panel is titled with dataname1, so dataname1 must be
# the treatment (set2) and dataname2 the control (set1).
tryCatch({
  suppressWarnings(
    gen_bagplot_bagfoot2(
      dat,
      dataname1               = set2,
      dataname2               = set1,
      qvaluethreshold_bagplot = opt$maha_qvalue,
      factor                  = 1.5,
      include_homer_results   = use_homer_plot,
      homer_results_file      = homer_results_file,
      pval_threshold_homer    = 0.05,
      maha_barplot            = run_maha_flag,
      maha_top_n              = opt$maha_top_n,
      maha_mode               = opt$maha_mode
    )
  )
  cat(paste0("\n>>> Success! Results in: ", bg_dir, "\n"))
}, error = function(e) {
  fail(paste0("Plot generation failed: ", e$message))
})

if (!is.null(dev.list())) invisible(dev.off())
if (file.exists("Rplots.pdf")) invisible(file.remove("Rplots.pdf"))
