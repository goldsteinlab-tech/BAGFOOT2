#' Get the default cache directory for BAGFOOT2's motif database files
#'
#' Resolves where genome support files (motifDB) should be stored/looked up.
#' Priority order: the BAGFOOT2_MOTIF_DB_DIR environment variable (useful on a
#' shared cluster where files already exist at a fixed path), then the
#' "bagfoot2.motif_db_dir" R option, then a default per-user cache directory.
#'
#' @return Character string. Path to the resolved cache directory.
#' @export
bagfoot2_cache_dir <- function() {

  # 1. Environment variable override (e.g. set in .bashrc on a cluster)
  env_dir <- Sys.getenv("BAGFOOT2_MOTIF_DB_DIR", unset = "")
  if (nzchar(env_dir)) return(env_dir)

  # 2. R option override (set via options() in a script/session)
  opt_dir <- getOption("bagfoot2.motif_db_dir", default = NULL)
  if (!is.null(opt_dir)) return(opt_dir)

  # 3. Default: per-user cache directory (created if missing)
  dir <- tools::R_user_dir("bagfoot2", which = "data")
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  dir
}

#' Download Genome Support Files (mm10/hg38)
#'
#' Downloads the required genome-wide motif database (motifDB) files from
#' Zenodo. These files are genome-wide FIMO motif scan results and are
#' necessary to run BAGFOOT2 for a given genome build.
#'
#' @param genome Character string. The genome build to download ("mm10" or "hg38").
#' @param destination_dir Character string. Where to save/extract the files.
#'   Defaults to \code{bagfoot2_cache_dir()}.
#' @param force Logical. If TRUE, re-downloads even if files already exist locally.
#' @param verify_checksum Logical. If TRUE (default), verifies the downloaded
#'   archive's MD5 checksum against the known value before extracting.
#' @return Character string. Path to the directory containing the extracted files.
#' @export
download_bagfoot2_files <- function(genome = "mm10",
                                 destination_dir = bagfoot2_cache_dir(),
                                 force = FALSE,
                                 verify_checksum = TRUE) {

  # 1. Define URLs (public Zenodo deposit, DOI 10.5281/zenodo.20782782)
  zenodo_id <- "20782782"
  base_url  <- sprintf("https://zenodo.org/records/%s/files", zenodo_id)

  urls <- list(
    "mm10" = file.path(base_url, "mm10_files.tar.gz?download=1"),
    "hg38" = file.path(base_url, "hg38_files.tar.gz?download=1")
  )

  # Known MD5 checksums, used to detect corrupted/incomplete downloads
  md5sums <- list(
    "mm10" = "8e1c73690d23dc35cc8b3e2b61f0b107",
    "hg38" = "8adadc106e89dcdaa600dbda68e1960f"
  )

  if (!genome %in% names(urls)) {
    stop("Genome ", genome, " is not supported. Choose 'mm10' or 'hg38'.")
  }

  # 2. Make sure the destination directory exists
  if (!dir.exists(destination_dir)) dir.create(destination_dir, recursive = TRUE)

  target_dir <- file.path(destination_dir, paste0(genome, "_files"))

  # 3. Skip download if files are already present
  if (dir.exists(target_dir) && !force) {
    message("Files for ", genome, " already present in ", target_dir, ". Skipping.")
    return(target_dir)
  }

  # 4. Download (large files: bump the timeout)
  dest_file <- file.path(destination_dir, paste0(genome, "_files.tar.gz"))
  message("Downloading ", genome, " files from Zenodo...")
  options(timeout = max(7200, getOption("timeout")))

  tryCatch({
    download.file(urls[[genome]], destfile = dest_file, mode = "wb", method = "libcurl")
  }, error = function(e) {
    stop("Download failed (", conditionMessage(e), "). Check your connection or the Zenodo URL.")
  })

  # 5. Verify integrity before extracting
  if (verify_checksum) {
    message("Verifying checksum...")
    actual_md5 <- tools::md5sum(dest_file)
    if (actual_md5 != md5sums[[genome]]) {
      unlink(dest_file)
      stop("Checksum mismatch for ", genome, " -- corrupted file, please retry the download.")
    }
    message("Checksum OK.")
  }

  # 6. Extract and clean up the archive
  message("Extracting...")
  untar(dest_file, exdir = destination_dir)
  unlink(dest_file)

  message("Done: files ready in ", target_dir)
  target_dir
}

#' Ensure motifDB files exist locally, downloading them only if needed
#'
#' Looks for an already-resolved motifDB directory for the given genome
#' (via \code{bagfoot2_cache_dir()}, which checks shared/cluster locations
#' first) before falling back to downloading from Zenodo. Validity is
#' checked against the actual files BAGFOOT2 expects: a
#' "motiflist_<genome>_cluster.txt" file and a "<genome>_cluster_fimo"
#' subdirectory. This is the function called by both the command-line
#' wrapper and direct R usage.
#'
#' @param genome Character string. The genome build needed ("mm10" or "hg38").
#' @param cache_dir Character string. Root directory to look in / download to.
#'   Defaults to \code{bagfoot2_cache_dir()}.
#' @param force Logical. If TRUE, forces a re-download even if files exist.
#' @return Character string. Path to the motifDB directory for this genome
#'   (i.e. the directory containing "motiflist_<genome>_cluster.txt" and
#'   "<genome>_cluster_fimo/").
#' @export
ensure_motif_db <- function(genome, cache_dir = bagfoot2_cache_dir(), force = FALSE) {

  target_dir <- file.path(cache_dir, paste0(genome, "_files"))

  # 1. Check that the expected files are actually there (not just the folder)
  motiflist_ok <- file.exists(file.path(target_dir, paste0("motiflist_", genome, "_cluster.txt")))
  fimo_ok      <- dir.exists(file.path(target_dir, paste0(genome, "_cluster_fimo")))

  if (motiflist_ok && fimo_ok && !force) return(target_dir)

  # 2. Not found / incomplete -> warn and download (prompt only if interactive)
  message(sprintf("motifDB files for %s not found/incomplete locally, downloading from Zenodo...", genome))

  if (interactive()) {
    ans <- readline(sprintf("Download %s now? [y/n] ", genome))
    if (!tolower(ans) %in% c("y", "yes")) {
      stop("Cancelled. Run ensure_motif_db('", genome, "') again when ready.")
    }
  }
  # On a cluster (non-interactive Rscript), download proceeds without a prompt.

  download_bagfoot2_files(genome = genome, destination_dir = cache_dir, force = force)
}

#' Build a MotifDB object for a given genome, downloading files if needed
#'
#' Convenience wrapper for users running BAGFOOT2 directly in R (as opposed to
#' the command-line wrapper): resolves/downloads the motifDB files for the
#' requested genome and returns a ready-to-use \code{MotifDB} object, with
#' paths built the same way \code{bagfoot2_wrapper.R} builds them.
#'
#' @param genome Character string. The genome build needed ("mm10" or "hg38").
#' @param cache_dir Character string. Root directory to look in / download to.
#'   Defaults to \code{bagfoot2_cache_dir()}.
#' @param force Logical. If TRUE, forces a re-download even if files exist.
#' @return A \code{MotifDB} object pointing to the resolved directory.
#' @export
motifDB_for_genome <- function(genome, cache_dir = bagfoot2_cache_dir(), force = FALSE) {
  mm_files_dir <- ensure_motif_db(genome, cache_dir = cache_dir, force = force)
  MotifDB(
    motiflistfile = file.path(mm_files_dir, paste0("motiflist_", genome, "_cluster.txt")),
    directory     = file.path(mm_files_dir, paste0(genome, "_cluster_fimo"))
  )
}
