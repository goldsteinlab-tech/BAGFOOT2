#' Download Genome Support Files (mm10/hg19)
#'
#' Downloads the required genome annotation and motif database files from Zenodo.
#' These files are necessary to run TIGER.
#'
#' @param genome Character string. The genome build to download (e.g., "mm10" or "hg19").
#' @param destination_dir Character string. Where to save the files. Default is "./tiger_files".
#' @param force Logical. If TRUE, overwrites existing files.
#' @return The path to the downloaded directory.
#' @export
download_tiger_files <- function(genome = "mm10", destination_dir = "./tiger_files", force = FALSE) {

  # 1. Define URLs

  urls <- list(
    "mm10" = "https://doi.org/10.5281/zenodo.18182843/mm10_files.tar.gz?download=1",
    "hg19" = "https://doi.org/10.5281/zenodo.18182843/hg19_files.tar.gz?download=1"
  )

  if (!genome %in% names(urls)) {
    stop("Genome ", genome, " is not supported yet. Choose 'mm10' or 'hg19'.")
  }

  #
  if (!dir.exists(destination_dir)) {
    dir.create(destination_dir, recursive = TRUE)
  }

  target_dir <- file.path(destination_dir, paste0(genome, "_files"))

  # 3. Check if exist
  if (dir.exists(target_dir) && !force) {
    message("Genome files found in ", target_dir, ". Skipping download.")
    return(target_dir)
  }

  # 4. download
  dest_file <- file.path(destination_dir, paste0(genome, "_files.tar.gz"))
  message("Downloading ", genome, " support files (~X GB). Please wait...")

  # Increase timeout for big files
  options(timeout = max(3600, getOption("timeout")))

  tryCatch({
    download.file(urls[[genome]], destfile = dest_file, mode = "wb")
  }, error = function(e) {
    stop("Download failed. Check internet connection or Zenodo link.")
  })

  # 5. unzip
  message("Unzipping files...")
  untar(dest_file, exdir = destination_dir)
  unlink(dest_file) # Supprimer l'archive

  message("Success! Files are ready in ", target_dir)
  return(target_dir)
}
