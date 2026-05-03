

#' Load Chromosome Ranges for Reference Genome
#'
#' Returns a list of chromosome ranges for a specified reference genome.
#'
#' @param refgenome Character. Reference genome name (e.g., "mm10", "mm9", "hg19", "hg18", "mm8","hg38").
#'
#' @return Named list of chromosome ranges (start and end positions).
#' @export
#'
#' @examples
#' #loadChromosomeRange("mm10")
loadChromosomeRange<-function(refgenome) {

  #ES added info from library(GenomicFeatures)   getChromInfoFromUCSC("mm10")
  hg38range=list('chr1'=c(1,195471971), 'chr2'=c(1,182113224),'chr3'=c(1,160039680),'chr4'=c(1,156508116),
                 'chr5'=c(1,151834684),'chr6'=c(1,149736546),'chr7'=c(1,145441459),'chr8'=c(1,129401213),
                 'chr9'=c(1,124595110),'chr10'=c(1,130694993),'chr11'=c(1,122082543),'chr12'=c(1,120129022),
                 'chr13'=c(1,120421639),'chr14'=c(1,124902244),'chr15'=c(1,104043685),'chr16'=c(1,98207768),
                 'chr17'=c(1,94987271),'chr18'=c(1,90702639),'chr19'=c(1,61431566),'chrX'=c(1,171031299),
                 'chrY'=c(1,91744698));

  mm10range=list('chr1'=c(1,195471971),  'chr2'=c(1,182113224),  'chr3'=c(1,160039680),  'chr4'=c(1,156508116),
                 'chr5'=c(1,151834684),  'chr6'=c(1,149736546),  'chr7'=c(1,145441459),  'chr8'=c(1,129401213),
                 'chr9'=c(1,124595110),  'chr10'=c(1,130694993),  'chr11'=c(1,122082543),  'chr12'=c(1,120129022),
                 'chr13'=c(1,120421639),  'chr14'=c(1,124902244),  'chr15'=c(1,104043685),  'chr16'=c(1,98207768),
                 'chr17'=c(1,94987271),  'chr18'=c(1,90702639),  'chr19'=c(1,61431566),  'chrX'=c(1,171031299),
                 'chrY'=c(1,91744698));

  mm9range=list('chr1'=c(1,197195432),  'chr2'=c(1,181748087),  'chr3'=c(1,159599783),  'chr4'=c(1,155630120),
                'chr5'=c(1,152537259),  'chr6'=c(1,149517037),  'chr7'=c(1,152524553),  'chr8'=c(1,131738871),
                'chr9'=c(1,124076172),  'chr10'=c(1,129993255),  'chr11'=c(1,121843856),  'chr12'=c(1,121257530),
                'chr13'=c(1,120284312),  'chr14'=c(1,125194864),  'chr15'=c(1,103494974),  'chr16'=c(1,98319150),
                'chr17'=c(1,95272651),  'chr18'=c(1,90772031),  'chr19'=c(1,61342430),  'chrX'=c(1,166650296),
                'chrY'=c(1,15902555));

  hg19range=list('chr1'=c(1,249250621),'chr2'=c(1,243199373),
                 'chr3'=c(1,198022430),'chr4'=c(1,191154276),
                 'chr5'=c(1,180915260),'chr6'=c(1,171115067),
                 'chr7'=c(1,159138663),'chr8'=c(1,146364022),
                 'chr9'=c(1,141213431), 'chr10'=c(1,135534747),
                 'chr11'=c(1,135006516), 'chr12'=c(1,133851895),
                 'chr13'=c(1,115169878), 'chr14'=c(1,107349540),
                 'chr15'=c(1,102531392), 'chr16'=c(1,90354753),
                 'chr17'=c(1,81195210),   'chr18'=c(1,78077248),
                 'chr19'=c(1,59128983),   'chr20'=c(1,63025520),
                 'chr21'=c(1,48129895),   'chr22'=c(1,51304566),
                 'chrX'=c(1,155270560),   'chrY'=c(1,59373566));

  hg18range = list(
    'chr1'=c(1,247249719),  'chr2'=c(1,242951149),
    'chr3'=c(1,199501827),  'chr4'=c(1,191273063),
    'chr5'=c(1,180857866),  'chr6'=c(1,170899992),
    'chr7'=c(1,158821424),  'chr8'=c(1,146274826),
    'chr9'=c(1,140273252),  'chr10'=c(1,135374737),
    'chr11'=c(1,134452384),  'chr12'=c(1,132349534),
    'chr13'=c(1,114142980),  'chr14'=c(1,106368585),
    'chr15'=c(1,100338915),  'chr16'=c(1,88827254),
    'chr17'=c(1,78774742),  'chr18'=c(1,76117153),
    'chr19'=c(1,63811651),  'chr20'=c(1,62435964),
    'chr21'=c(1,46944323),  'chr22'=c(1,49691432),
    'chrX'=c(1,154913754),  'chrY'=c(1,57772954));

  mm8range= list('chr1'=c(1,197069962),  'chr2'=c(1,181976762),
                 'chr3'=c(1,159872112),  'chr4'=c(1,155029701),
                 'chr5'=c(1,152003063),  'chr6'=c(1,149525685),
                 'chr7'=c(1,145134094),  'chr8'=c(1,132085098),
                 'chr9'=c(1,124000669),  'chr10'=c(1,129959148),
                 'chr11'=c(1,121798632),  'chr12'=c(1,120463159),
                 'chr13'=c(1,120614378),  'chr14'=c(1,123978870),
                 'chr15'=c(1,103492577),  'chr16'=c(1,98252459),
                 'chr17'=c(1,95177420),  'chr18'=c(1,90736837),
                 'chr19'=c(1,61321190),  'chrX'=c(1,165556469),
                 'chrY'=c(1,16029404))

  genomeinfo=list("mm8"=mm8range, "mm9"=mm9range, "mm10"=mm10range,"hg18"=hg18range, "hg19"=hg19range, "hg38"=hg38range);

  if (!refgenome %in% names(genomeinfo)) {
    stop(sprintf('%s : not supported yet.',refgenome));
  }

  genomeinfo[[refgenome]];

}


# call the function writeBGRperChromosomeInt in the code bagfoot.so compiled from a C code bagfoot_calc.c

#' @param bgrfilename
#'
#' @param chr
#' @param count
#' @param binsize
#' @param threshold
#'
# @useDynLib bagfoot writeBGRperChromosomeInt
#' @useDynLib tiger writeBGRperChromosomeInt
writeBGRperChromosomeR<-function(bgrfilename, chr, count, binsize=20, threshold=2) {
  countsize= length(count);

  out <- .C("writeBGRperChromosomeInt", filename = as.character(bgrfilename),chrom=as.character(chr),
            count=as.integer(count), span=as.integer(countsize), binsize=as.integer(binsize),
            thr=as.numeric(threshold));
}


#' Read Nucleotide Codes at Cut Sites
#'
#' Reads nucleotide codes from a file and returns codes at cut sites, shifted by half the n-mer size.
#'
#' @param nuccodefile Character. Path to nucleotide code RDS file.
#' @param np Integer. n-mer size.
#'
#' @return Integer vector of nucleotide codes at cut sites.
#' @export
#'
#' @examples
#' # readNucleotideCodeForChromosomeForCuts("nuccode_chr1_4mer.dat", 4)
readNucleotideCodeForChromosomeForCuts<-function(nuccodefile,np) {

	nuccode= readRDS(nuccodefile);
	nuccodeAtCuts = shiftarray(nuccode$code, np/2);
	nuccodeAtCuts[1:np]=4^np+1;
	nuccode={};
#	print(nuccodefile);
	as.integer(nuccodeAtCuts);
}

#' Calculate Expected Cleavage Rates
#'
#' Calculates expected cleavage rates for genomic regions based on nucleotide codes and a frequency table.
#'
#' @param tregion Data frame. Regions with chromosome, start, and end positions.
#' @param dir Data frame or vector. Direction information for each region.
#' @param nuccodes List. Nucleotide codes for each chromosome.
#' @param ftable Data frame. Cleavage rates table.
#'
#' @return Matrix of expected cleavage rates for each region.
#' @export
#'
#' @examples
#' # CalcExpectedRates(tregion, dir, nuccodes, ftable)
CalcExpectedRates<- function(tregion, dir=NA,nuccodes=NA,ftable=NA) {

  # ADD tiger : no hexamer table
  if (is.null(ftable) || length(ftable) == 0 || all(is.na(ftable))) {
    nrow_result <- nrow(tregion)
    if (nrow_result == 0) return(matrix(1, nrow=0, ncol=1))
    ncol_result <- max(tregion$ed - tregion$st) + 1
    return(matrix(1, nrow=nrow_result, ncol=ncol_result))
  }
  ###
  # Get unique sorted list of chromosomes from tregion
  locchrom = as.character(unique(sort(tregion$chr)));
  numchrom = length(locchrom);
  Plist = {};
  nregion = nrow(tregion);
  range=seq(from=1, to=nregion);
  chrs=unique(as.character(tregion$chr));

  # Function to calculate expected rates for each chromosome
  calcPerchrom<-function(ch) {
    nuccode = nuccodes[[ch]];  # Get nucleotide codes for the chromosome
    region = subset(tregion, chr==ch);  # Subset regions for the chromosome
    dirchr = subset(dir, tregion$chr==ch);  # Subset directions for the chromosome
    numregion = nrow(region);
    R={};
    inRange= region$ed < length(nuccode);  # Check if region end is within nucleotide code length
    outRange = !inRange;
    regionInRange = region[inRange,];  # Subset regions within range

    # Function to calculate Pj for each site
    calcPj<-function(ii) {
      site= regionInRange[ii,];
      Aj_CleavageRate=ftable$Aj[nuccode[site$st:site$ed]];  # Get cleavage rates for the site
      Pj=Aj_CleavageRate / sum(Aj_CleavageRate);  # Normalize cleavage rates
      if (dirchr[ii]==1) {
        returnvalue=Pj;  # Return Pj if direction is 1
      } else {
        returnvalue = rev(Pj);  # Reverse Pj if direction is not 1
      }
      returnvalue;
    }

    # Apply calcPj function to each region in range
    if (MCCORES==1) {
      QQ<-lapply(1:nrow(regionInRange) ,calcPj);
    } else {
      QQ<-parallel::mclapply(1:nrow(regionInRange) ,calcPj, mc.cores=MCCORES);
    }
    Qj= do.call("rbind",QQ)  # Combine results into a matrix
    QC =array(0,dim=c(length(inRange), ncol(Qj)));  # Initialize QC array
    QC[inRange,] = Qj;  # Fill QC with Qj values for regions in range
    QC;
  }

  # Apply calcPerchrom function to each chromosome
  Pj = parallel::mclapply(locchrom, calcPerchrom, mc.cores=MCCORES);

  # Function to filter out NaN results
  filterNANresult <- function(SS) {
    if (is.numeric(SS)) {
      nullrow = is.nan(apply(SS,1,sum));  # Identify rows with NaN values
      SS[nullrow,] = array(0, dim=c(sum(nullrow), ncol(SS)));  # Replace NaN rows with zeros
      result = SS;
    } else {
      result = NA;
    }
    result;
  }

  Pjj = lapply(Pj, filterNANresult);  # Apply filterNANresult to each element in Pj
  ExpectedCuts=do.call(rbind, Pjj);  # Combine results into a single data frame
  ExpectedCuts;
}
