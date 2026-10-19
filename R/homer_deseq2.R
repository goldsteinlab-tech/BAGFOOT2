


#' Process Hotspot File
#'
#' Processes a hotspot file by removing unnecessary columns, formatting, and writing a new file.
#'
#' @param input_file Character. Path to the input hotspot file.
#' @param output_file Character. Path to the output formatted file.
#' @param sep_in Character. Input file separator (default: ",").
#' @param sep_out Character. Output file separator (default: "\\t").
#' @param skip_lines Integer. Number of lines to skip at the start of the input file (default: 1).
#'
#' @return None. Writes a formatted file to \code{output_file}.
#' @export
#'
#' @examples
#' # process_hotspot_file("hotspot.csv", "hotspot_formatted.txt")
process_hotspot_file <- function(
  input_file,
  output_file,
  sep_in     = ",",
  sep_out    = "\t",
  skip_lines = 1,
  remove_cols_first = 1,
  remove_cols_after = 4:7,
  strip_quotes = TRUE
) {
  if (!file.exists(input_file)) {
    stop("Input file does not exist: ", input_file)
  }

  df <- tryCatch(
    read.table(
      input_file,
      header = FALSE,
      sep = sep_in,
      skip = skip_lines,
      quote = "",
      stringsAsFactors = FALSE,
      check.names = FALSE
    ),
    error = function(e) {
      stop("Failed to read input_file: ", conditionMessage(e))
    }
  )

  if (ncol(df) == 0) {
    stop("Input file has no columns after skipping lines, cannot process.")
  }

  ## 1)delete first column
  if (!is.null(remove_cols_first)) {
    if (any(remove_cols_first < 1 | remove_cols_first > ncol(df))) {
      stop("remove_cols_first contains invalid column indices.")
    }
    df <- df[, -remove_cols_first, drop = FALSE]
  }

  ## 2) remplace ""
  if (strip_quotes) {
    df[] <- lapply(df, function(x) {
      if (is.character(x)) {
        gsub('"', "", x, fixed = TRUE)
      } else {
        x
      }
    })
  }

  ## 3) delete columns 4–7 after first deletion
  if (!is.null(remove_cols_after)) {
    if (length(remove_cols_after) > 0) {
      if (any(remove_cols_after < 1 | remove_cols_after > ncol(df))) {
        stop(
          "remove_cols_after contains invalid column indices. ",
          "Data has ", ncol(df), " columns after first removal."
        )
      }
      df <- df[, -remove_cols_after, drop = FALSE]
    }
  }

  ## 4) writing
  dir_out <- dirname(output_file)
  if (!dir.exists(dir_out)) {
    dir.create(dir_out, recursive = TRUE)
  }

  write.table(
    df,
    file      = output_file,
    sep       = sep_out,
    quote     = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )

  message("File processed and saved to: ", output_file)
  invisible(output_file)
}


#######################################################################################################
#' Create HOMER Tag Directories
#'
#' Creates HOMER tag directories for each BAM file in a working directory.
#'
#' @param homer_bin Character. Path to HOMER binaries.
#' @param work_dir Character. Working directory containing BAM files.
#' @param condition_order Character vector or NULL. Optional order of conditions.
#' @param verbose Logical. If TRUE, prints progress messages (default: TRUE).
#'
#' @return Character vector of created tag directory names.
#' @export
#'
#' @examples
#' # create_tag_directories("/path/to/homer/bin", "/path/to/bam_dir")
create_tag_directories <- function(
  homer_bin,
  work_dir,
  condition_order = NULL,
  pattern_bam = "_[0-9]+\\.bam$",
  verbose = TRUE
) {
  if (!dir.exists(work_dir)) {
    stop("work_dir does not exist: ", work_dir)
  }
  if (!dir.exists(homer_bin)) {
    warning("homer_bin directory does not exist: ", homer_bin,
            " (makeTagDirectory may not be found in PATH).")
  }

  old_path <- Sys.getenv("PATH")
  on.exit(Sys.setenv(PATH = old_path), add = TRUE)
  Sys.setenv(PATH = paste(homer_bin, old_path, sep = .Platform$path.sep))

  if (verbose) {
    message("HOMER bin PATH: ", homer_bin)
    message("BAM directory: ", work_dir)
  }

  bam_files <- list.files(
    path       = work_dir,
    pattern    = pattern_bam,
    full.names = FALSE
  )

  if (length(bam_files) == 0) {
    if (verbose) message("No BAM replicate files found in directory (pattern: ", pattern_bam, ").")
    return(character(0))
  }

  conditions <- sub(pattern_bam, "", bam_files)
  uniq_cond  <- unique(conditions)

  if (!is.null(condition_order)) {
    missing <- setdiff(condition_order, uniq_cond)
    extra   <- setdiff(uniq_cond, condition_order)

    if (length(missing) > 0) {
      stop("Conditions missing in condition_order: ",
           paste(missing, collapse = ", "))
    }
    if (length(extra) > 0 && verbose) {
      warning("Conditions present in BAM files but not in condition_order: ",
              paste(extra, collapse = ", "))
    }

    cond_factor <- factor(conditions, levels = condition_order)
    rep_num <- suppressWarnings(as.integer(sub(".*_([0-9]+)\\.bam$", "\\1", bam_files)))
    if (any(is.na(rep_num))) {
      stop("Cannot extract replicate numbers from BAM filenames; ",
           "ensure they match the pattern 'condition_<number>.bam'.")
    }

    ord <- order(cond_factor, rep_num)
    bam_files  <- bam_files[ord]
    conditions <- conditions[ord]
  }

  if (verbose) {
    message("BAM -> Tag Directory order:")
    for (i in seq_along(bam_files)) {
      message(sprintf("  %d: %s (condition: %s)", i, bam_files[i], conditions[i]))
    }
  }

  for (cond in unique(conditions)) {
    replicates <- grep(paste0("^", cond, "_[0-9]+\\.bam$"),
                       bam_files, value = TRUE)
    if (length(replicates) < 2 && verbose) {
      message("Warning: No replicate BAM file for condition: ", cond)
    }
  }

  td_dirs <- character(length(bam_files))

  for (i in seq_along(bam_files)) {
    bam      <- bam_files[i]
    base     <- sub("\\.bam$", "", bam)
    td_dir   <- paste0("TD_", base)
    td_path  <- file.path(work_dir, td_dir)
    td_dirs[i] <- td_dir

    bam_path <- file.path(work_dir, bam)

    if (dir.exists(td_path) && length(list.files(td_path)) > 0L) {
      if (verbose) {
        message("Directory ", td_dir,
                " exists and is not empty, skipping creation.")
      }
      next
    }

    if (!file.exists(bam_path)) {
      warning("BAM file does not exist (skipping): ", bam_path)
      next
    }

    if (verbose) {
      message("Processing: ", td_dir)
    }

    cmd <- sprintf("makeTagDirectory %s %s",
                   shQuote(td_path), shQuote(bam_path))

    if (verbose) message("> ", cmd)

    status <- system(cmd)
    if (status != 0) {
      warning("makeTagDirectory failed for ", bam_path,
              " with status ", status)
    }
  }

  invisible(td_dirs)
}

# create_tag_directories <- function(homer_bin,
#                                    work_dir,
#                                    condition_order = NULL,  # nouveau
#                                    verbose = TRUE) {
#
#   Sys.setenv(PATH = paste(Sys.getenv("PATH"), homer_bin, sep = ":"))
#   if (verbose) cat("HOMER bin PATH :", homer_bin, "\n")
#   if (verbose) cat("bam directory :", work_dir, "\n")
#
#   bam_files <- list.files(path = work_dir,
#                           pattern = "_[0-9]+\\.bam$",
#                           full.names = FALSE)
#
#   if (length(bam_files) == 0) {
#     if (verbose) cat("No BAM replicate files found in directory.\n")
#     return(character(0))
#   }
#
#   # conditions déduites des noms
#   conditions <- sub("_[0-9]+\\.bam$", "", bam_files)
#   uniq_cond  <- unique(conditions)
#
#   # si un ordre est donné, réordonner
#   if (!is.null(condition_order)) {
#     missing <- setdiff(condition_order, uniq_cond)
#     extra   <- setdiff(uniq_cond, condition_order)
#     if (length(missing) > 0) {
#       stop("Conditions manquantes dans condition_order: ",
#            paste(missing, collapse = ", "))
#     }
#     if (length(extra) > 0) {
#       warning("Conditions présentes dans les BAM mais pas dans condition_order: ",
#               paste(extra, collapse = ", "))
#     }
#
#     # ordre: par condition_order puis par numéro de réplica
#     cond_factor <- factor(conditions, levels = condition_order)
#     rep_num <- as.integer(sub(".*_([0-9]+)\\.bam$", "\\1", bam_files))
#     ord <- order(cond_factor, rep_num)
#     bam_files <- bam_files[ord]
#     conditions <- conditions[ord]
#   }
#
#   if (verbose) {
#     cat("Ordre des BAM -> TD :\n")
#     for (i in seq_along(bam_files)) {
#       cat(sprintf("  %d: %s (condition: %s)\n",
#                   i, bam_files[i], conditions[i]))
#     }
#   }
#
#   # check réplicats
#   for (cond in unique(conditions)) {
#     replicates <- grep(paste0("^", cond, "_[0-9]+\\.bam$"),
#                        bam_files, value = TRUE)
#     if (length(replicates) < 2) {
#       if (verbose) cat("There is no replicate bam file for condition:", cond, "\n")
#     }
#   }
#
#   td_dirs <- character(length(bam_files))
#   for (i in seq_along(bam_files)) {
#     bam <- bam_files[i]
#     base <- sub("\\.bam$", "", bam)
#     td_dir <- paste0("TD_", base)
#     td_path <- file.path(work_dir, td_dir)
#     td_dirs[i] <- td_dir
#
#     bam_path <- file.path(work_dir, bam)
#
#     if (dir.exists(td_path)) {
#       td_files <- list.files(td_path)
#       if (length(td_files) > 0) {
#         if (verbose) cat("Directory", td_dir,
#                          "exists and is not empty, skipping creation.\n")
#         next
#       }
#     }
#
#     if (verbose) cat("Processing:", td_dir, "\n")
#     cmd <- sprintf("makeTagDirectory %s %s",
#                    shQuote(td_path), shQuote(bam_path))
#     if (verbose) cat(">", cmd, "\n")
#     system(cmd)
#   }
#
#   return(td_dirs)
# }


#' Run HOMER Annotate Peaks
#'
#' Runs HOMER's annotatePeaks.pl on  peaks using specified tag directories.
#'
#' @param homer_bin Character. Path to HOMER binaries.
#' @param work_dir Character. Working directory.
#' @param merged_peak Character. Path to merged peak file.
#' @param genome Character. Genome identifier (e.g., "hg38", "mm10").
#' @param td_dirs Character vector. Tag directories.
#' @param annotate_outfile Character. Output file for annotation results.
#' @param verbose Logical. If TRUE, prints progress messages (default: TRUE).
#'
#' @return None. Writes annotation results to \code{annotate_outfile}.
#' @export
#'
#' @examples
#' # run_annotate_peak("/path/to/homer/bin", "/work/dir", "merged_peaks.bed", "mm10", c("TD_sample1", "TD_sample2"), "annotate.txt")

run_annotate_peak <- function(homer_bin, work_dir, merged_peak, genome, td_dirs, annotate_outfile, verbose = TRUE) {
  #Add homer in the PATH
  Sys.setenv(PATH = paste(Sys.getenv("PATH"), homer_bin, sep = ":"))
  if (verbose) cat("HOMER bin has been added in PATH :", homer_bin, "\n")

  # Change directory
  setwd(work_dir)

  # Tag directories
  td_dirs_arg <- paste(td_dirs, collapse = " ")

  # build command annotatePeaks.pl
  cmd_annotate <- sprintf(
    "annotatePeaks.pl %s %s -noann -noadj -d %s > %s",
    merged_peak,
    genome,
    td_dirs_arg,
    annotate_outfile
  )

  # Print and execute the commande
  if (verbose) cat("Running annotatePeaks.pl...\n>", cmd_annotate, "\n")
  system(cmd_annotate)
}


#' Run Differential Expression Analysis with HOMER
#'
#' Runs HOMER's getDiffExpression.pl to perform differential expression analysis between two groups.
#'
#' @param homer_bin Character. Path to HOMER binaries.
#' @param work_dir Character. Working directory.
#' @param annotate_outfile Character. Path to annotation output file.
#' @param diffexp_outfile Character. Output file for differential expression results.
#' @param group1 Character vector. Sample names for group 1.
#' @param group2 Character vector. Sample names for group 2.
#' @param verbose Logical. If TRUE, prints progress messages (default: TRUE).
#'
#' @return None. Writes differential expression results to \code{diffexp_outfile}.
#' @export
#'
#' @examples
#' # run_diff_expression("/path/to/homer/bin", "/work/dir", "annotate.txt", "diffexp.txt", c("sample1"), c("sample2"))

run_diff_expression <- function(homer_bin, work_dir, annotate_outfile, diffexp_outfile, group1, group2, verbose = TRUE) {
  Sys.setenv(PATH = paste(Sys.getenv("PATH"), homer_bin, sep = ":"))
  setwd(work_dir)

  # Prepar group args
  group_args <- paste(
    annotate_outfile,
    paste(group1, collapse = " "),
    paste(group2, collapse = " "),
    "-peaks -DESeq2"
  )

  cmd_diffexp <- sprintf("getDiffExpression.pl %s > %s", group_args, diffexp_outfile)

  if (verbose) cat("Running getDiffExpression.pl...\n>", cmd_diffexp, "\n")
  system(cmd_diffexp)
}

#' Create Differential Peak Files
#'
#' Creates files for differential peaks based on log2 fold change and adjusted p-value cutoffs.
#'
#' @param input_path Character. Directory containing the input file.
#' @param input_filename Character. Name of the input file.
#' @param sep Character. Separator used in the input file (default: "\\t").
#' @param chr Character. Column name for chromosome (default: "Chr").
#' @param start Character. Column name for start position (default: "Start").
#' @param end Character. Column name for end position (default: "End").
#' @param padj_cutoff Numeric. Adjusted p-value cutoff (default: 0.05).
#' @param l2fc_cutoff Numeric. Log2 fold change cutoff (default: 1).
#' @param output_prefix Character. Prefix for output files (default: "diff_peaks").
#'
#' @return List with two data frames: \code{full} (detailed results) and \code{simple} (BED format).
#' @export
#'
#' @examples
#' # create_diff_peak_files("/input/dir", "diffexp.txt")
create_diff_peak_files <- function(
  input_path,
  input_filename,
  sep = "\t",
  chr = "Chr",
  start = "Start",
  end = "End",
  padj_cutoff = 0.05,
  l2fc_cutoff = 1,
  output_prefix = "diff_peaks"
) {
  # Build the path
  input_file <- file.path(input_path, input_filename)
  # Read file
  df <- read.table(input_file, sep = sep, stringsAsFactors = FALSE, fill=TRUE, header=TRUE, quote = "")
  colnames(df) <- gsub(" ", ".", colnames(df))

  # Find column indexes
  l2fc_col_idx <- grep("Log2.Fold.Change", colnames(df))
  padj_col_idx <- grep("adj..p.value", colnames(df))

  if (any(sapply(list(l2fc_col_idx, padj_col_idx), length) == 0)) {
    stop("Impossible to find columns of the differential analysis in input file.")
  }

  # Rename via index (robuste)
  colnames(df)[l2fc_col_idx] <- "log2FoldChange"
  colnames(df)[padj_col_idx] <- "padj"

  # Check coordinates columns existence
  required_cols <- c(chr, start, end, "log2FoldChange", "padj")
  missing <- setdiff(required_cols, colnames(df))
  if (length(missing) > 0) stop(paste("Missing columns:", paste(missing, collapse = ", ")))

  # Cut-off/filtres
  filt <- subset(df, !is.na(padj) & padj <= padj_cutoff & abs(log2FoldChange) >= l2fc_cutoff)
  if (nrow(filt) == 0) stop("No differential peak.")

  # Annotation UP/DOWN
  filt$regulation <- ifelse(filt$log2FoldChange > 0, "up", "down")
  # Coordinates
  filt$coord <- paste0(filt[[chr]], ":", filt[[start]], "-", filt[[end]])

  # Organization of output tables
  main_df <- filt[, c(chr, start, end, "coord", "regulation", "log2FoldChange", "padj")]
  simple_df <- main_df[, c(chr, start, end, "coord", "regulation")]

  # Write files
  write.table(main_df, paste0(output_prefix, "_full.txt"), sep="\t", row.names=FALSE, quote=FALSE)
  write.table(simple_df, paste0(output_prefix, "_simple.bed"), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

  cat("Created files:\n",
      paste0(output_prefix, "_full.txt\n"),
      paste0(output_prefix, "_simple.bed\n"))

  return(list(full = main_df, simple = simple_df))
}

# create_diff_peak_files <- function(
#   input_path,
#   input_filename,
#   sep = "\t",
#   chr = "Chr",
#   start = "Start",
#   end = "End",
#   padj_cutoff = 0.05,
#   l2fc_cutoff = 1,
#   output_prefix = "diff_peaks"
# ) {
#   # Build the path
#   input_file <- file.path(input_path, input_filename)
#
#
#
#   # read file
#   df <- read.table(input_file, sep = sep, stringsAsFactors = FALSE,fill=TRUE,header=TRUE,, quote = "")
#   colnames(df) <- gsub(" ", ".", colnames(df))
#
#   # Find column names
#   l2fc_col <- grep("Log2.Fold.Change", colnames(df), value = TRUE)
#   #pval_col <- grep("p.value", colnames(df), value = TRUE)
#   padj_col <- grep("adj..p.value", colnames(df), value = TRUE)
#
#   if (any(sapply(list(l2fc_col, padj_col), length) == 0)) {
#     stop("Impossible to find columns of the differential analysis in input file.")
#   }
#
#   # Rename
#   colnames(df)[colnames(df) == l2fc_col] <- "log2FoldChange"
#   #colnames(df)[colnames(df) == pval_col] <- "pvalue"
#   colnames(df)[colnames(df) == padj_col] <- "padj"
#
#   #
#   #   # Renomme les colonnes d'intérêt
#   # colnames(df)[l2fc_col_idx] <- "log2FoldChange"
#   # colnames(df)[pval_col_idx] <- "pvalue"
#   # colnames(df)[padj_col_idx] <- "padj"
#
#   # # Check  coordinates column
#   required_cols <- c(chr, start, end, "log2FoldChange", "padj")
#   missing <- setdiff(required_cols, colnames(df))
#   if(length(missing) > 0) stop(paste("missing columns:", paste(missing, collapse = ", ")))
#
#
#   # Cut-off
#   filt <- subset(df, !is.na(padj) & padj <= padj_cutoff & abs(log2FoldChange) >= l2fc_cutoff)
#
#   if(nrow(filt) == 0) stop("no differential peak.")
#
#
#   # Annotation up/down
#   filt$regulation <- ifelse(filt$log2FoldChange > 0, "up", "down")
#
#   # coordinates
#   filt$coord <- paste0(filt[[chr]], ":", filt[[start]], "-", filt[[end]])
#
#   # organisation of the file
#   main_df <- filt[, c(chr, start, end, "coord", "regulation", "log2FoldChange", "padj")]
#
#   # Simplified file
#   simple_df <- main_df[, c(chr, start, end, "coord", "regulation")]
#
#   # File saved
#   write.table(main_df, paste0(output_prefix, "_full.txt"), sep="\t", row.names = FALSE, quote = FALSE)
#   write.table(simple_df, paste0(output_prefix, "_simple.bed"), sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
#
#   cat("Created files :\n",
#       paste0(output_prefix, "_full.txt\n"),
#       paste0(output_prefix, "_simple.bed\n"))
#
#
#   return(list(full = main_df, simple = simple_df))
# }


