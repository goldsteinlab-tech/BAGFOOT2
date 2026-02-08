## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----echo=FALSE---------------------------------------------------------------
knitr::kable(
  data.frame(
    Step = c("Preprocessing", "Hexamer", "Peaks"),
    Output_File = c("BedGraph", "nuccode + Hexamer table", "peaks reformatted as CSV"),
    Description = c("Genome-wide Accessibility profile", "hexamer genome index + hexamer bias / frequency table", "peaks in standardized CSV format")
  )
)

## ----echo=FALSE---------------------------------------------------------------

 knitr::kable(
  data.frame(
    Step = c(
      "Generate nuccode files ** (once per genome)",
      "Generate BedGraph (incl. Tn5 offset correction)",
      "Build hexamer table",
      "Peak reformatting"
    ),
    `Output (per condition)` = c(
      "Nuccode* files (hexamer index in the genome, shared by all conditions)",
      "BedGraph file (shifted)",
      "Hexamer bias file",
      "Hotspot table (MACS2 peaks reformatted to MASA format; one row per peak with placeholder footprint metrics)"
    ),
    `Output (per Couple)` = c(
      "",
      "",
      "Combined hexamer (merge hexamer tables for both conditions)",
      "Combined peaks / pooled hotspot  (merged hotspot table defining a common  peaks universe for both conditions)"
    )
  ),
  align = "lcc",
  caption = "Summary Table"
)



## ----eval=FALSE---------------------------------------------------------------
# 
# library(masa)
# library(hash)
# library(data.table)
# library(Cairo)
# library(aplpack)
# library(stringr)
# library(dplyr)
# library(ggplot2)
# library(plotly)
# 

## ----include=FALSE,eval=FALSE-------------------------------------------------
# numcores <- 2     # No. of CPU cores to be used
# options("masa.mc.cores"=numcores)
# options("masa.mc.cores.motif"=numcores)
# options("masa.debug"=T)
# 
# project<-"urf"
# 
# # the names of the folders must be as below
# bg_dir <-paste0("/home/leslieco/bg_files_fed_fasted/")
# 
# macs2_dir <- paste0("/home/leslieco/",project,"/macs2/")
# 
# bam_dir <- paste0("/home/leslieco/",project,"/bam_merged/")
# 
# # mm10_files directory should be downloaded and copied in the masa directory
# mm_files_dir <- "/home/leslieco/mm10_files/"
# 
# mm<-"mm10"
# genome <- "mm10"
# 
# 
# #⚠️ in couples: CONTROL is the FIRST, followed by - (not underscore) flowed by the second
# 
# if (project=="urf"){
#   all_set<- c("fed","fasted")
#   couples<- c("fed-fasted") # the control here is fed
# }
# 
# knitr::opts_knit$set(root.dir = bg_dir)

## ----eval=FALSE---------------------------------------------------------------
# setwd(bg_dir)
# 
# 
# for (treat in all_set) {
#   # Build BAM file path
#   bamfile <- file.path(bam_dir, paste0(treat, ".sorted.bam"))
# 
#   # Count the number of cuts in the BAM file
#   cc <- countReadsBAM(bamfile, mm)
#   print(paste("treat:", treat, "cc:", cc))
#   assign(paste0("cc_", treat), cc)
# 
#   # Generate a BedGraph file with cleavage counts
#    cutcountfile = makeCutCountBAMWithName(bamfile,
#                                          name =  basenameWithoutExt(bamfile),
#                                          refgenome=mm)
# 
#   # Tn5 insertion bias correction
#    hex_file <- file.path(bg_dir, paste0("Hexamer_", treat, "_", mm, ".txt"))
# message("Checking Hexamer file: ", hex_file)
# message("Exists? ", file.exists(hex_file))
#    if (!file.exists(hex_file)) {
#     message("Creating Hexamer table (and possibly nuccode) for ", treat, " ...")
# 
#   tab = MakeBiasCorrectionTableBAM(
#   bamfile=bamfile,
#   outfile=paste0("Hexamer_",treat,"_mm10.txt"),
#   refgenome="mm10",
#   np=6,
#   atac=T   )  # Set TRUE for ATAC data. The default value is FALSE
# 
#    } else {
#     message("Hexamer file already exists for ", treat, " -> skipping bias table generation.")
#    }
# 
# 
#   # Read MACS2 peaks and convert to CSV
#   macs <- read.csv(file.path(macs2_dir, paste0(treat, "_peaks.narrowPeak")),
#                    sep = "\t", header = FALSE, stringsAsFactors = FALSE)
#   names(macs)[1:3] <- c('chr', 'st', 'ed')
#   nr <- nrow(macs)
# 
#   # Create a hotspot data frame for each peak
#   hotspot <- data.frame(
#     ID = 1:nr, chrom. = macs$chr, start = macs$st + 1, end = macs$ed,
#     MaxD = rep(0, nr),
#     AveD = rep(0, nr),
#     Zscore = rep(0, nr),
#     pvalue = rep(0, nr)
#   )
# 
#   # Write the hotspot table to CSV
#   write.table(hotspot, file = file.path(bg_dir, paste0(treat, "_peaks.csv")), sep = ",", row.names = FALSE)
# }
# 
# # Remove temporary files
# temp <- list.files(bg_dir,
#                    pattern = glob2rx("temp*.dat"),
#                    recursive = TRUE,
#                    full.names = TRUE)
# file.remove(temp)
# 
# # Merge files for MASA run
# for (c in couples) {
#   print(c)
#   set1 <- strsplit(c, "-", fixed = TRUE)[[1]][1]
#   set2 <- strsplit(c, "-", fixed = TRUE)[[1]][2]
# 
#   # Combine hotspots for the two conditions
#   assign(
#     paste0(set1, "_", set2, "_peaks"),
#     combineTwoHotspots(
#       file.path(bg_dir, paste0(set1, "_peaks.csv")),
#       file.path(bg_dir, paste0(set2, "_peaks.csv")),
#       set1, set2
#     )
#   )
# 
#   # Read and clean pooled hotspot file
#   pooled_path <- file.path(bg_dir, paste0("pooled_", set1, "_", set2, "_hotspot.csv"))
#   pooled <- read.csv(pooled_path)
#   if ("ID.1" %in% colnames(pooled)) {
#     pooled <- pooled %>%
#       dplyr::select(-ID.1)
#     write.csv(pooled, pooled_path, row.names = FALSE)
#   }
# 
#   # Merge the hexamer frequency tables for the two conditions
#   assign(
#     paste0("Hexamer_", set1, "_", set2),
#     merge_two_frequency_tables(
#       file.path(bg_dir, paste0("Hexamer_", set1, "_", mm, ".txt")),
#       file.path(bg_dir, paste0("Hexamer_", set2, "_", mm, ".txt")),
#       file.path(bg_dir, paste0("Hexamer_", set1, "_", set2, "_merged.csv"))
#     )
#   )
# }
# 
# 

## ----echo=TRUE, message=FALSE, warning=FALSE ,results = 'hide',eval=FALSE-----
# 
# 
# # ============================
# # SETTINGS
# # ============================
# 
# # Path to the HOMER binary directory
# homer_bin <- "/home/leslieco/software/homer/bin/"
# 
# # Directory containing all BAM files (replicates)
# work_dir <- bam_dir
# 
# # Peak file merged from previous MASA steps
# merged_peak <- file.path(bg_dir, "pooled_fasted_fed_hotspot_modify.txt")
# 
# # Directory for HOMER output
# output_dir <- file.path(bg_dir, paste0("homer_res_", project))
# 
# # Output file for annotatePeak.pl
# annotate_outfile <- file.path(bg_dir, "fasted_annot.txt")
# 
# # Output file for getDiffExpression.pl
# diffexp_outfile <- file.path(bg_dir, "fasted_output.txt")
# 
# # Reference genome (e.g., mm10, hg38)
# genome <- "mm10"
# 
# # Replicate groups for getDiffExpression.pl (adapt as needed)
# group1 <- c("fed", "fed", "fed") # 3 replicas
# group2 <- c("fasted", "fasted", "fasted") # 3 replicas
# 
# # Reformat the merged peak file for HOMER
# input_file1 <- file.path(bg_dir, "pooled_fed_fasted_hotspot.csv")
# output_file1 <- file.path(bg_dir, "pooled_fed_fasted_hotspot_modify.txt")
# 
# process_hotspot_file(
#   input_file1,
#   output_file1
# )
# 
# library("DESeq2")
# 
# # 1. Create tag directories for each BAM file (replicate)
# # This will create tag directories in the merged BAM directory
# # BAM files must be named with the pattern "_[0-9].bam"
# # Can take more than 30 min
# 
# # PATH expected for the Tag Directories
# 
# bam_files_all <- list.files(work_dir, pattern = "\\.bam$", full.names = TRUE)
# bam_files <- bam_files_all[grepl("[0-9]\\.bam$", bam_files_all)]
# td_paths <-  file.path(work_dir, paste0("TD_", sub("\\.bam$", "", basename(bam_files))))
# 
# # Check of TDs
# if (all(dir.exists(td_paths))) {
#   message("All Tag Directories already exist, skipping creation.")
#   td_dirs <- td_paths
# } else {
#   td_dirs <- create_tag_directories(
#   homer_bin ,
#   work_dir ,
#   condition_order = c("fed", "fasted"),  # control first, treatment second !!!
#   verbose = TRUE
# )
# }
# 
# 
# # 2. Peak annotation using HOMER
# # This will generate an annotation file in bg_dir
# run_annotate_peak(
#   homer_bin,
#   work_dir,
#   output_file1,
#   genome,
#   td_dirs,
#   annotate_outfile
# )
# 
# # 3. Differential analysis with HOMER
# # This will generate a diff expression file in bg_dir
# run_diff_expression(
#   homer_bin,
#   work_dir,
#   annotate_outfile,
#   diffexp_outfile,
#   group1,
#   group2
# )
# 
# # 4. Create a file with only differential peaks (BED format and table with l2fc, padj, pval) in work_dir(bam_merged)
# 
# diffexp_outfile1<-"fasted_output.txt"
# 
# create_diff_peak_files(
#   input_path = bg_dir,
#   input_filename = diffexp_outfile1,
#   sep = "\t",
#   chr = "Chr",
#   start = "Start",
#   end = "End",
#   padj_cutoff = 0.05, # adjust as needed
#   l2fc_cutoff = 0.5,    # adjust as needed
#   output_prefix = "my_peaks"
# )
# 
# # copy the BED file (not usually necessary in R, but shown for completeness) to bg_dir
# file.copy(
#   from = file.path(work_dir, "my_peaks_simple.bed"),
#   to = file.path(bg_dir, "my_peaks_simple.bed"),
#   overwrite = TRUE
# )
# 
# # 5. Run HOMER motif enrichment
# input_bed <- file.path(bg_dir, "my_peaks_simple.bed")
# output_folder <- file.path(bg_dir, "homer_motif_res")
# genome_build <- "mm10"
# homer_exe <- file.path(homer_bin, "findMotifsGenome.pl")
# 
# MASA_DIR <- "/home/leslieco/masa/"  # <-- User edits this line for their install
# motif_file <- file.path(MASA_DIR, "/inst/data/all_cluster_motifs.motifs")
# 
# extra_args1 <- paste0("-mknown ", motif_file)
# 
# # homer enrichment , can take more than 30min
# run_homer(
#   input_bed,
#   output_folder,
#   genome_build,
#   homer_exe,
#   extra_args1
# )
# 

## ----echo=TRUE,message=FALSE, warning=FALSE,results='hide', fig.show='hide',eval=FALSE----
# 
# # Set working directory to bg_dir
# 
# knitr::opts_knit$set(root.dir = "/home/leslieco/bg_files_fed_fasted/")
# 
# # Load libraries
# 
# library(masa)
# library(dplyr,tidyr)
# library(stringr)
# library(ggplot2)
# library(ggExtra)
# library(plotly)
# 
# # Base directories
# 
# bg_dir <- "/home/leslieco/bg_files_fed_fasted/"
# mm_files_dir <- "/home/leslieco/mm10_files/"
# mm <- "mm10"
# 
# # Analysis parameters
# 
# numcores <- 2
# motif_no <- 577 # total number of motifs in motifDB. This number should not be changed
# 
# 
# # Treatments and couples
# 
# all_set <- c("fed", "fasted")
# couples <- c("fed-fasted")
# 
# # HOMER motif results
# 
# homer_results_file <- file.path(bg_dir, "homer_motif_res", "knownResults.txt")
# 
# # Load motif database (FIMO-based archetype motif scan results)
# 
# motifdb_mm <- MotifDB(
# motiflistfile = file.path(mm_files_dir, paste0("motiflist_", mm, "_cluster.txt")),
# directory = file.path(mm_files_dir, paste0(mm, "_cluster_fimo"))
# )
# 
# # Initialize CutCount object from BedGraph files
# # Each CutCount object stores ATAC-seq cut profiles for a given condition
# 
# for (treat in all_set) {
#   assign(
#   treat,
#   CutCount(
#   file = file.path(bg_dir, paste0(treat, "_AC_cutcount.bgr.gz")),
#   count = 1,
#   name = treat))
# }
# 
# # Run MASA to generate footprint_depth_table files per couple
# 
# for (c in couples) {
#   message("Processing couple: ", c)
# 
#   set1 <- strsplit(c, "-", fixed = TRUE)[[1]][1] # control condition
#   set2 <- strsplit(c, "-", fixed = TRUE)[[1]][2] # treatment condition
# 
# # Define bias-correction options using the merged hexamer table
#   gfootoption_couple <- GFootOption(
#   biasfile = file.path(bg_dir, paste0("Hexamer_", set1, "_", set2, "_merged.csv")),
#   hexamer_pattern = paste0("nuccode", mm, "6mer{chr}.dat")
#   )
# 
#   setwd(bg_dir)
# 
# # Check whether footprint_depth_table outputs already exist (In case a restart was forced)
# 
#   output_files<- list.files(bg_dir,
#                             pattern= glob2rx(paste0(set1,"_",set2,"_On_pooled_",set1,"_",set2,
#                                                     "_hotspot_footprint_depth_table",".csv")),
#                             recursive = TRUE)
#   if (length(output_files)==0) {
#     start_motif<-1
#   } else{
#     # Resume: infer starting motif index from the last output folder name
#     start_motif<-dirname(output_files) %>%  # directory name
#       str_split(.,"\\_", simplify = TRUE) %>%   # separate by underscore
#       as.data.frame (stringsAsFactors=F)%>%
#       dplyr::select(tail(names(.), 1)) %>%  # keep only the last column
#       mutate_if(is.character,as.numeric) %>%  #convert to numeric
#       max+25 #Identify the last processed index and increment by 25 to start the next batch
#   }
# 
# # Run MASA in batches of 25 motifs
#   if (start_motif<motif_no) {
# 
#     for (r in seq(start_motif, motif_no, by = 25)) {
#       GFoot_obj <- GFoot(
#       control = get(set1),
#       treatment = get(set2),
#       sitefile = paste0("pooled_",set1,"_",set2,"_hotspot.csv"),
#       motifDB = motifdb_mm,
#       gfootoption = gfootoption_couple,
#       outputdir = file.path(bg_dir, paste0("OUTPUT_", c, "_", r)),
#       cachedir = file.path(bg_dir, paste0("CACHE_", c, "_", r))
#       )
# 
# 
#     run(
#       GFoot_obj,
#       graphout = T,
#       yrange = c(-2.2, 1.0),
#       mc.cores = numcores,
#       range=r:ifelse(r+24<motif_no,r+24,motif_no),
#       run_set=c
#     )
#     gc() #garbage collector
#   }
# 
#   setwd(bg_dir)
# 
# # Merge footprint depth tables produced for all motif batches
# # (central + flanking accessibility information for all motifs)
# 
#   merge_output<-paste0("find  . -type f -name '",
#                         set1,"_",set2,"_On_pooled_",set1,"_",set2,
#                         "_hotspot_footprint_depth_table.csv' -exec tail -n +2 -q {} ';'",
#                         ">",set1,"_",set2,"_On_pooled_",
#                          set1,"_",set2,"_hotspot_footprint_depth_table_merged.csv")
# 
#     system(merge_output)
# 
#     #read the merged footprint depth table
# # This table includes both central (motif) and flanking accessibility metrics
# 
#     dat= read.table(paste0(bg_dir,set1,"_",set2,"_On_pooled_",set1,"_",set2,"_hotspot_footprint_depth_table_merged.csv"), quote="\"", comment.char="")
# 
#     #add column names
#     colnames(dat)<-c("num","filename","center","logo","motiflength","motif","memeNo","memeEntry","lrmotif","lrflanking",
#                      "lrtotal","fdepth1","fdepth2","meancc1","meancc2","numsite")
# 
#  # Generate MASA plots and interactive report
# # The bagplot summarizes TF activity changes based on flanking + motif-core accessibility
# 
#    test_gen<-gen_bagplot_masa(dat,
#                               dataname1=set1,
#                               dataname2=set2,
#                               qvaluethreshold_bagplot=0.05,
#                               factor=1.5,
#                               pdf = TRUE,
#                               html=TRUE,
#                               include_homer_results = TRUE,
#                               homer_results_file = homer_results_file,
#                               pval_threshold_homer = 0.05)
# 
#      dev.off()
#   }
# }

## ----eval=TRUE----------------------------------------------------------------
# Session info for reproducibility/debugging
sessionInfo()


