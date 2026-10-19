#' Merge Two Frequency Tables
#'
#' Merges two frequency tables with identical genomic position counts, sums observed cuts, and calculates ratios and correction factors.
#'
#' @param file1 Character. Path to the first frequency table file.
#' @param file2 Character. Path to the second frequency table file.
#' @param outfile Character. Path to the output file where the merged table will be written.
#'
#' @return A data frame containing the merged frequency table with columns: ObCuts, GenomicPositionCount, ObCutRatio, GPRatio, and CorrectionFactor.
#' The merged table is also written to the specified output file.
#' @export
#'
merge_two_frequency_tables<- function(file1, file2, outfile) {
  t1<-read.table(file1);
  t2<-read.table(file2);

  if (all(t1$GenomicPositionCount==t2$GenomicPositionCount)) {
    nr= nrow(t1);
    ObCuts <- t1$ObCuts + t2$ObCuts;
    GenomicPositionCount= t1$GenomicPositionCount
    GPRatio <- c(t1$GenomicPositionCount[-nr] / sum(as.double(t1$GenomicPositionCount[-nr])), NA);
    ObCutRatio <- c(ObCuts[-nr] / sum(as.double(ObCuts[-nr])),NA);
    CorrectionFactor=c(GPRatio[-nr]/ObCutRatio[-nr], 1);

    t<-data.frame(ObCuts=ObCuts,GenomicPositionCount=GenomicPositionCount , ObCutRatio=ObCutRatio,GPRatio=GPRatio,CorrectionFactor=CorrectionFactor);

    row.names(t) <- row.names(t1);
    write.table(t, file=outfile);
  } else {
    stop('Different ref. genomes were used for Table 1 and Table 2.');
  }
}

merged_two_frequency_table = merge_two_frequency_tables;
