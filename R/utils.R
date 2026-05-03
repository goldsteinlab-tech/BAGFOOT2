# Utilitaires used in this package

#' Print Formatted Output
#'
#' Prints formatted output to the console, similar to printf in C.
#'
#' @param ... Arguments passed to sprintf.
#'
#' @return None.
#' @export
#'
#' @examples
#' #printf("Hello %s\n", "world")
printf <- function(...) { cat(sprintf(...)); };

#' Extract Base Name Without Extension
#'
#' Returns the base name of a file, removing its extension.
#'
#' @param filename Character. Path to the file.
#'
#' @return Character. File name without its extension.
#' @export
#'
#' @examples
#' #extractBaseName("path/to/file.txt")
extractBaseName<- function(filename) {
  bname = basename(filename);
  ext=regexpr("\\.",bname)[[1]];
  bnameNoExt = substring(bname, 1,ext-1);
  bnameNoExt;
}

#' Count Cut Events on Sites
#'
#' Counts the number of cut events within specified genomic sites, grouped by chromosome.
#'
#' @param cutcount List. Cut count data per chromosome.
#' @param sites Data frame or character. Genomic sites or path to sites file.
#'
#' @return Numeric. Total cut count across all sites.
#' @export
#'
#' @examples
#' #countCutcountOnSites(cutcount, sites)
countCutcountOnSites <- function(cutcount, sites) {
  print("countCutcountOnSites")
  print(paste("sites:",sites,"is.character(sites):",is.character(sites)))
  if (is.character(sites)) {
    hs <- readcsv(sites);
  } else {
    hs <- sites;
  }

  hs$chr <- factor(hs$chr)
  chroms = levels(hs$chr);
  sum_for_chr<- function(chr) {
    tdchr = cutcount[[chr]];
    hschr = hs[hs$chr==chr,];
    cc = vector("numeric", length=max(hs$ed+1));
    cc[tdchr$st] = tdchr$val;
    cc[tdchr$ed+1] = cc[tdchr$ed+1] - tdchr$val;
    cc = cumsum(cc);
    s=sum(sapply(seq_along(hschr$st),function(ii) sum(cc[hschr$st[ii]:hschr$ed[ii]])));
    rm(cc);
    s;
  }
  sum_for_chr_es<-sum(unlist(lapply(chroms, sum_for_chr)))
  print(paste("sum_for_chr_es:",sum_for_chr_es))
  sum_for_chr_es;
}

#' Read Cut Count Data from BedGraph
#'
#' Reads cut count data from a BedGraph file and splits it by chromosome.
#'
#' @param cutcountfile Character. Path to the BedGraph file.
#'
#' @return List of data frames, one per chromosome.
#' @export
readcutcount2 <- function(cutcountfile) {
  cutcount = readbgr(cutcountfile);
  rl = rle(as.vector(cutcount$chr));
  chrs = rl$values;
  ed=cumsum(rl$lengths);
  st= c(1,ed[-length(ed)]+1);
  cutcountPerChrom=lapply(1:length(chrs), function(ii) {cutcount[st[ii]:ed[ii],];});
  names(cutcountPerChrom) = chrs;
  cutcount={};
  cutcountPerChrom;
}

#' Get Output Directory and Filename
#'
#' Constructs output directory and filename for results based on input parameters.
#'
#' @param outputdir_base Character. Base output directory.
#' @param dataname Character. Data name.
#' @param sitefilename Character. Site file name.
#' @param halfwidth Integer. Half-width value.
#'
#' @return List with elements 'dir' and 'filename'.
#' @export
getOutputDirAndFilename<- function(outputdir_base, dataname,  sitefilename, halfwidth) {
  hotspotname = gsub(".csv","",basename(sitefilename));
  od = file.path(outputdir_base, sprintf("%s_on%s_%dbp", dataname,hotspotname, halfwidth));
  list(dir = od,
       filename = file.path(od, sprintf('../%s_on%s_%dbp_out.txt', dataname, hotspotname,halfwidth)));
}

#' Read BAM Index
#'
#' Reads BAM index statistics and returns a data frame for chromosomes in the reference genome.
#'
#' @param bamFile Character. Path to BAM file.
#' @param refgenome Character. Reference genome ("mm9", "mm10", "hg19","hg38).
#'
#' @return Data frame with columns chr, maxloc, num, etc.
#' @export
readBAMIndex<-function(bamFile, refgenome="mm9") {
  #bamFile<-bamfile
  #if (system("samtools", intern=F)!=1) {
  #		stop('Cannot run "samtools".  Install "samtools" first!');
  #	};

  if (refgenome=="mm9" | refgenome=="mm10") {
    chroms = paste("chr", c(1:19,"X","Y"), sep="");
  }

  if (refgenome =="hg19" | refgenome=="hg38") {
    chroms = paste("chr", c(1:22,"X","Y"), sep="");
  }

  baiFile = sprintf('%s.bai',bamFile);
  if (!isBigEnough(baiFile)) {
    dlog(sprintf("Generating the index file of %s\n", bamFile));
    system(sprintf('samtools index %s', bamFile));
  }
  #idxout= system(sprintf('samtools idxstats %s', bamFile), intern=T);
  #con <- textConnection(idxout);
  dat <- read.csv(sprintf('%s.csv',bamFile),sep='\t',header=F);
  names(dat) <- c("chr","maxloc","num","etc");
  dat[dat$chr %in% chroms,];
}


#' Count Reads in BAM File
#'
#' Counts the total number of reads in a BAM file for the specified reference genome.
#'
#' @param bamFile Character. Path to BAM file.
#' @param refgenome Character. Reference genome (default: "mm10").
#'
#' @return Integer. Total read count.
#' @export
countReadsBAM<-function(bamFile, refgenome = "mm10") {
  #bamFile = fastedbam;
  dat = readBAMIndex(bamFile, refgenome= refgenome);
  ncount = sum(dat$num);
  ncount;
}


#' Check if File is Big Enough
#'
#' Checks if a file exists and its size is greater than a threshold.
#'
#' @param filename Character. Path to the file.
#'
#' @return Logical. TRUE if file exists and is large enough, FALSE otherwise.
#' @export
isBigEnough<-function(filename) {
  allowablesize = 300;
  if (file.exists(filename)) {
    ifelse(file.info(filename)$size > allowablesize , T,F);
  } else {
    F;
  }

}



#' Make CutCount BedGraph File from BAM
#'
#' Creates a BedGraph file containing cut counts from a BAM file.
#'
#' @param bamfile Character. Path to the input BAM file.
#' @param refgenome Character. Reference genome (e.g., "mm9", "hg19").
#'
#' @return Character. BedGraph filename.
#' @export
makeCutCountBAM<- function(bamFile, refgenome = "mm10") {
  makeCutCountBAMWithName(bamFile, name=basenameWithoutExt(bamFile), refgenome=refgenome);
}


#' Get Base Name Without Extension
#'
#' Returns the base name of a file, removing its extension.
#'
#' @param filename Character. Path to the file.
#'
#' @return Character. File name without its extension.
#' @export
basenameWithoutExt<-function(filename)  {
  sub("^([^.]*).*", "\\1",basename(filename))
}



#' Generate n-mer Bias Tables from BAM
#'
#' Generates n-mer bias tables from a BAM file, optionally using mappability files.
#'
#' @param bamfile Character. Path to the input BAM file.
#' @param outfile Character. Output filename.
#' @param refgenome Character. Reference genome.
#' @param np Integer. Number of base pairs (default: 6).
#' @param mapdir Character. Directory to mappability files.
#' @param atac Logical. TRUE for ATAC-seq data.
#'
#' @return Character. Output filename.
#' @export
MakeBiasCorrectionTableBAM<- function(bamfile="", outfile="", refgenome="", np=6, mapdir='',atac=F) {
  #examples:
  #tab=MakeBiasCorrectionTableBAM(bamfile="/mnt/Data1/MA/sjbaek/project/ghostfootprint/data/nakedDNA/NakedDNA_SRR769954_hg19_sorted.bam", outfile="Hexamer_FT_NakedDNA_hg19_withMap.txt", refgenome="hg19", np=6, mappability=T);
  #tab2=MakeBiasCorrectionTableBAM(bamfile="/mnt/Data1/MA/sjbaek/project/ghostfootprint/data/nakedDNA/NakedDNA_SRR769954_hg19_sorted.bam", outfile="Hexamer_FT_NakedDNA_hg19_withoutMap.txt", refgenome="hg19", np=6, mappability=F)
  #tab3=MakeBiasCorrectionTableBAM(bamfile= "/home/sjbaek/Data/data/stamseq/DS8497/DS8497_sorted_merged.bam", outfile="Hexamer_FT_DS8497_mm9_withMap.txt", refgenome="mm9", np=6, mappability=T);
  setwd(bg_dir)

  if (!file.exists(outfile)) {
    if (!atac) {
      shifts = c(0,0);
    } else {
      shifts = c(4,-5);
      cat(sprintf("The reads are shifted by (%d,%d) base-pairs for ATAC-seq data.\n", shifts[1],shifts[2]));
    }

    freq=calcFreqeuncyTableBAM(bamfile =bamfile,  refgenome=refgenome,  np=np,mapdir=mapdir,shifts=shifts);
    #	browser();
    colnames(freq) <- c("ObCuts", "GenomicPositionCount");
    gpsum = sum(as.numeric(freq$GenomicPositionCount[-nrow(freq)]));
    totalcuts =sum(freq$ObCuts[-nrow(freq)]);   # Both sum exclude the last element (others)
    ObCutRatio = freq$ObCuts / totalcuts;
    GPRatio = freq$GenomicPositionCount / gpsum;
    #AjRatio = freq$ObCuts / freq$GenomicPositionCount;
    ObCutRatio[nrow(freq)] = NA;
    GPRatio[nrow(freq)] = NA;
    CorrectionFactor= GPRatio/ObCutRatio;
    CorrectionFactor[nrow(freq)] = 1;
    tab=cbind(freq, ObCutRatio, GPRatio, CorrectionFactor);
    write.table(tab, file=outfile);
  }
  outfile;
}

#' Make CutCount BedGraph File from BAM with Name
#'
#' Creates a BedGraph file containing cut counts from a BAM file, using a sample name.
#'
#' @param bamFile Character. Path to the input BAM file.
#' @param name Character. Sample name.
#' @param refgenome Character. Reference genome (default: "mm10").
#' @param atac Logical. TRUE for ATAC-seq data (default: FALSE).
#'
#' @return Character. Path to the compressed BedGraph file (.bgr.gz).
#' @export
makeCutCountBAMWithName<- function(bamFile, name="", refgenome = "mm10",atac=F) {
  #OutputDir =  dirname(bamFile);
  OutputDir = getwd();
  bgrfilename = file.path(OutputDir, sprintf('%s_AC_cutcount.bgr',name));
  bgrgzfilename = paste(bgrfilename,'.gz', sep='');

  bamindex = readBAMIndex(bamFile);  # Make sure that the index file is generated before creating cutcounts
  # browser();
  if (!isBigEnough(bgrgzfilename)) {
    BuildCutCountProfileBAM(bamfile = bamFile, bgrfilename=bgrfilename, ftablefile = '',  refgenome=refgenome, compress=T,atac=atac);
  };
  bgrgzfilename;
}


#' Combine Two Hotspot Files
#'
#' Combines two hotspot (peak) files in CSV format by taking the union of hotspots.
#'
#' @param csvfile1 Character. Path to first hotspot file.
#' @param csvfile2 Character. Path to second hotspot file.
#' @param name1 Character. Name for first file.
#' @param name2 Character. Name for second file.
#'
#' @return Character. Path to the pooled hotspot file.
#' @export
combineTwoHotspots <- function(csvfile1, csvfile2, name1, name2) {
  #Outputdir = dirname(csvfile1);
  Outputdir = getwd();
  pooledhotspotfilename = file.path(Outputdir, sprintf('pooled_%s_%s_hotspot.csv',name1,name2));
  hotspotfiles = list(csvfile1, csvfile2);
  hotspotlist = lapply(hotspotfiles, readcsv);
  phs = pool_hotspots3(hotspotlist);
  phs = data.frame(ID=1:length(phs$chr),phs);
  write.table(phs, file=pooledhotspotfilename, sep=",", row.names=F);
  pooledhotspotfilename;
}


#' Read CSV File with Standardized Columns
#'
#' Reads a CSV file and renames columns 'chr', 'st', and 'ed' as needed.
#'
#' @param filename Character. Path to the CSV file.
#'
#' @return Data frame with standardized column names.
#' @export
readcsv <- function(filename)
{
  R=read.table(filename, sep=",", header=TRUE);
  #  header=c("chr","st","ed","MaxD","AveD","Zscore");
  idx=sapply(c('chr','st'),function(x) {specialgrep(x, names(R)); });
  names(R)[idx[["chr"]]] = "chr";
  names(R)[idx[["st"]]] = "st";
  names(R)[idx[["st"]]+1] = "ed";
  # browser();
  #   nc = min(5,ncol(R)+1);
  #    names(R)[2:nc] = header[1:(nc-1)];
  #readcsv$st = readcsv$st+1;
  R;
}


#' Special Grep for Column Names
#'
#' Searches for column names starting with a given pattern, escaping special characters.
#'
#' @param x Character. Pattern to search for.
#' @param y Character vector. Vector of names to search in.
#' @param ... Additional arguments passed to grep.
#'
#' @return Integer vector of matching indices.
#' @export
specialgrep <- function(x,y,...){
  grep(
    paste("^",
          gsub("([].^+?|[#\\-])","\\\\\\1",x)
          ,sep=""),
    y,...)
}


#' Pool Hotspots from Multiple Files
#'
#' Merges two or more MACS2 peak files into a single hotspot file.
#'
#' @param hotspotdata List of data frames. Hotspot data from multiple files.
#'
#' @return Data frame. Merged hotspot data.
#' @export
pool_hotspots3<- function(hotspotdata) {

  chrs = as.character(sort(unique(unlist(lapply(hotspotdata, function(x) {unique(x$chr);})))));  numData = length(hotspotdata);
  ihotspots = vector("list", length=numData);

  R={};
  if (length(hotspotdata)<2)
  {
    R=hotspotdata[[1]];
    R;
    return;
  }
  ihotspots = vector("list", length=numData);
  maxnumhotspot = 0;
  for (j in 1:length(hotspotdata))
  {
    ihotspots[[j]]= 1:length(hotspotdata[[j]]$st);
    maxnumhotspot = maxnumhotspot + length(hotspotdata[[j]]$st);
  }

  R$chr = vector("character", length=maxnumhotspot);
  R$st = vector("numeric", length=maxnumhotspot);
  R$ed = vector("numeric", length=maxnumhotspot);

  range_st= 1;
  for (i in 1:length(chrs))
  {
    jj=1;
    hotspots=vector("list");
    range_hotspots =  vector("list", length= numData);
    for (j in 1:numData) {
      hotspots[[j]]=vector("list");
      range_hotspots[[j]] = ihotspots[[j]][hotspotdata[[j]]$chr== chrs[i]];
      hotspots[[j]]$st= hotspotdata[[j]]$st[range_hotspots[[j]]];
      hotspots[[j]]$ed= hotspotdata[[j]]$ed[range_hotspots[[j]]];
    }
    st=hotspots[[1]]$st;
    ed=hotspots[[1]]$ed;

    for (j in 2:numData) {
      st=c(st,hotspots[[j]]$st);
      ed=c(ed,hotspots[[j]]$ed);
    }
    if (length(st)>0)
    {
      ord= order(st);
      st1 = st[ord];  ed1 = ed[ord];
      leng = length(st);
      sst = st1[1];  cmp = ed1[1];

      chr2= vector("character", length = leng);
      st2= vector("numeric", length = leng);
      ed2= vector("numeric", length = leng);

      ii=2;
      while (ii<=leng)
      {
        if (cmp < st1[ii])
        {
          chr2[jj]= chrs[i];  st2[jj] = sst;     ed2[jj] = cmp;

          jj=jj+1;
          cmp = ed1[ii];    sst = st1[ii];
        } else {
          cmp = max(cmp,ed1[ii]);
        }
        if (ii==leng)
        {
          chr2[jj] = chrs[i];
          st2[jj] = sst;
          ed2[jj] = cmp;
        }
        ii=ii+1;
      }

      range_ed = range_st+jj-1;
      R$chr[range_st:range_ed] = chr2[1:jj];
      R$st[range_st:range_ed] = st2[1:jj];
      R$ed[range_st:range_ed] = ed2[1:jj];
      range_st = range_st+jj;
    }
  } #i
  R$chr=R$chr[1:range_ed];
  R$st=R$st[1:range_ed];
  R$ed=R$ed[1:range_ed];

  hhs =  as.data.frame(matrix(0, nrow=range_ed, ncol=ncol(hotspotdata[[1]])));
  names(hhs) <- names(hotspotdata[[1]]);
  if ("ID" %in% names(hhs)) {
    hhs$ID = 1:range_ed;
  }
  hhs$chr = R$chr;
  hhs$st = R$st;
  hhs$ed = R$ed;
  hhs;
}

#' Merge Two Frequency Tables
#'
#' Merges two frequency tables with identical genomic positions, sums observed cuts, and computes ratios and correction factors.
#'
#' @param file1 Character. Path to the first frequency table.
#' @param file2 Character. Path to the second frequency table.
#' @param outfile Character. Path to the output file.
#'
#' @return None. Writes the merged table to \code{outfile}.
#' @export
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

merged_two_frequency_table = merge_two_frequency_tables



### Input
#' Build CutCount Profile from BAM
#'
#' Builds a cut count profile from a BAM file and writes it to a BedGraph file.
#'
#' @param bamfile Character. Path to BAM file.
#' @param dataname Character. Data name.
#' @param bgrfilename Character. Output BedGraph filename.
#' @param ftablefile Character. Frequency table file.
#' @param refgenome Character. Reference genome.
#' @param chroms Character vector. Chromosome names.
#' @param compress Logical. Compress output file (default: TRUE).
#' @param atac Logical. TRUE for ATAC-seq data (default: FALSE).
#'
#' @return None.
#' @export
BuildCutCountProfileBAM <- function (bamfile='', dataname='', bgrfilename='', ftablefile='', refgenome='hg19', chroms=NA, compress=T, atac=F) {

  if (!atac) {
    shifts = c(0,0);
  } else {
    shifts = c(4,-5);
    cat(sprintf("The reads are being shifted by (%d,%d) base-pairs for ATAC-seq data.\n", shifts[1],shifts[2]));
  }

  if (dataname=='') {
    dataname= sub('.bam','',basename(bamfile),ignore.case=T);
  }

  if (bgrfilename=='') {
    if (any(shifts==0)) {
      bgrfilename= file.path(sprintf("%s_adjusted_by_%d_%d_cutcount%s.bgr",dataname, shifts[1],shifts[2], ifelse(ftablefile=='','',sprintf('_%s', basename(ftablefile)))));
    } else {
      bgrfilename= file.path(sprintf("%s_cutcount%s.bgr",dataname, ifelse(ftablefile=='','',sprintf('_%s', basename(ftablefile)))));
    }
  }

  np=NA;

  if (ftablefile != '') {
    biascorrection=T;
    ftable=read.table(ftablefile);
    np= nchar(rownames(ftable)[1]);
  } else {
    biascorrection=F;
  }


  heightPixel=c(128,80,0);
  binsize=20;

  trackcolor='black';


  bgrdir = dirname(bgrfilename);
  bgrbase = basename(bgrfilename);
  bgrbasenoext = sub('.bgr','', bgrbase, ignore.case = T);
  filetype = '';


  gboption = sprintf("visibility=full maxHeightPixels=%d:%d:%d autoScale=on",
                     heightPixel[1],heightPixel[2],heightPixel[3]);
  description = sprintf('%s %s',
                        dataname, getDateStamp());

  chromrange = loadChromosomeRange(refgenome);

  if (is.na(chroms)) {
    chroms = names(chromrange);
  }

  tempfiles=paste(file.path(bgrdir, bgrbasenoext),'_',chroms,'.bgr', sep='');
  bgrtempfiles = hash::hash(chroms,tempfiles);

  writeBGRheader(bgrfilename, dataname, description, gboption, col2ColorCode(trackcolor));

  genome <- readGenome(refgenome);

  MakeCutCountProfilesForChromosome<-function(chr)   {
    cuts<-readCutSitesPerChromFromBAM(bamfile, chr, shifts);

    cuts <- cuts[cuts>0]; # delete negative indices

    #	cat(sprintf('Reads were shifted by %d and %d bps.\n', shifts[1],shifts[2]));
    #		browser();
    cutdensity = rep(0, chromrange[[chr]][2]);
    rl = rle(cuts);  # cuts is already sorted.
    if (biascorrection==T) {
      nuccode=readNucleotideCodeForChromosome(chr, nmer=np, genome=genome, nuccodefileDir='.', mappability=T);
      nuccodeAtCuts = shiftarray(nuccode$code, np/2);
      nuccodeAtCuts[1:np]=4^np+1;
      nuccode={};
      bias = ftable$CorrectionFactor;
      rl$lengths  = rl$lengths[rl$values != 0];
      rl$values=rl$values[rl$values != 0];
      cutdensity[rl$values] = rl$lengths * bias[nuccodeAtCuts[rl$values]];
    } else {
      cutdensity[rl$values] = rl$lengths;
    }

    if (any(is.na(cutdensity))) {
      cutdensity[is.na(cutdensity)] = 0;
    }
    #
    writeBGRperChromosomeR(bgrtempfiles[[chr]], chr, cutdensity, binsize=20, threshold=0);
    #	cat(sprintf('creating cut count profiles for %s...\n', chr));
    cutdensity={};
  }

  #	l<-parallel::mclapply(sample(chroms), MakeCutCountProfilesForChromosome, mc.cores=6);
  l<-lapply(sample(chroms), MakeCutCountProfilesForChromosome);

  for (chr in chroms) {
    tempfile = bgrtempfiles[[chr]];
    if (file.exists(tempfile)) {
      cmd1 = sprintf('cat %s >>%s',tempfile, bgrfilename);
      cmd2 = sprintf('rm %s',tempfile);
      system(cmd1, intern=TRUE);
      system(cmd2, intern=TRUE);
    }
  }

  if (compress==TRUE) {
    system(sprintf("gzip -f %s", bgrfilename), wait=FALSE);
  }
  #	browser();

}

#' Get Date Stamp
#'
#' Returns the current date as a formatted string.
#'
#' @return Character. Date stamp in "mm-dd-yy" format.
#' @export
getDateStamp <- function() {
  format(Sys.time(), "%m-%d-%y");
}

#' Write BedGraph Header
#'
#' Writes a header line to a BedGraph file.
#'
#' @param filename Character. Path to the BedGraph file.
#' @param name Character. Track name.
#' @param desc Character. Track description.
#' @param option Character. Track options.
#' @param color Character. Track color code.
#'
#' @return None.
#' @export
writeBGRheader<-function(filename, name, desc, option, color) {

  filep = file(filename,"w+");
  str=sprintf( "track type=bedGraph name=\"%s\" description=\"%s\" %s color=%s\n",
               name, desc, option, color);
  cat(str, file=filep);
  close(filep);
}

##' Convert Color Name to RGB Code
#'
#' Converts a color name to an RGB string for use in track headers.
#'
#' @param col Character. Color name.
#'
#' @return Character. RGB color code as "R,G,B".
#' @export
col2ColorCode<-function(col) {
  cc=col2rgb(col)
  sprintf('%d,%d,%d', cc[1],cc[2],cc[3]);
}

#' Read Reference Genome
#'
#' Loads a reference genome as a BSgenome object.
#'
#' @param refgenome Character. Reference genome name.
#'
#' @return BSgenome object.
#' @export
readGenome <- function(refgenome) {
  if (refgenome == 'mm9') {
    genome <- BSgenome.Mmusculus.UCSC.mm9::BSgenome.Mmusculus.UCSC.mm9
  } else if (refgenome == 'mm10') {
    genome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
  } else if (refgenome == 'hg19') {
    genome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
  } else if (refgenome == 'hg38') {
    genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  } else {
    stop(sprintf("Unsupported genome: %s", refgenome))
  }
  genome
}

#' Read BedGraph File
#'
#' Reads a BedGraph file and returns a data frame with columns chr, st, ed, and value.
#'
#' @param filename Character. Path to the BedGraph file.
#'
#' @return Data frame with columns chr, st, ed, value.
#' @export
readbgr <- function(filename)
{
  if (filename=="")
  {
    R=NULL;
    R;
  } else {
    # Read the first 7 lines of the file
    R= readLines(filename, n=7);
    i=1;
    # Find the first line that starts with "chr"
    while (i < 7)
    {
      if (substr(R[[i]],1,3)=="chr") break;
      i=i+1;
    }
    # If no such line is found, return an empty list
    if (i >= 7)
    {
      R={};
      sprintf("%s is not a readable bed file.", A);
      return;
    } else
    {
      # Read the file into a data frame, skipping lines until the "chr" line
      R=read.csv(filename, sep="", header=FALSE,comment.char="", skip=i-1, colClasses=c("factor","numeric","numeric","numeric"), stringsAsFactors=FALSE);
      # Set column names
      names(R) <- c("chr","st","ed","value");
      # Increment the start position by 1
      R$st = R$st+1;
      R;
    }
  }
}

# @useDynLib bagfoot intersect_intervals
#' @useDynLib tiger intersect_intervals
intersect_intervals_fastInC <- function(X1,X2,Y1,Y2)
{
	Isorted1= order(X1);
	Isorted2=  order(Y1);
	.Call("intersect_intervals", as.numeric(X1),as.numeric(X2),as.numeric(Y1),as.numeric(Y2),as.integer(Isorted1), as.integer(Isorted2));
}


# @useDynLib bagfoot intersect_intervals_index
#' @useDynLib tiger intersect_intervals_index
intersect_intervals_fast_indexInC <- function(X1,X2,Y1,Y2) {
	Isorted1= order(X1);
	Isorted2=  order(Y1);
	.Call("intersect_intervals_index", as.numeric(X1),as.numeric(X2),as.numeric(Y1),as.numeric(Y2),as.integer(Isorted1), as.integer(Isorted2));
}


intersect_intervals_fast_index <- intersect_intervals_fast_indexInC;
intersect_intervals_fast <- intersect_intervals_fastInC;
intersect_intervals_really_fast = intersect_intervals_fastInC;


#' @export
compare_intervals <- function(x1,x2,y1,y2) ifelse(y1>x2, -1, ifelse(y2<x1,1,0));

#' @export
intersect_intervals_fast_parallel <- function(P) {

    R={};
    if (length(P$Ast) ==0  || length(P$Bst)==0) {
      return(R);
    }
    R=intersect_intervals_fast(P$Ast, P$Aed, P$Bst,P$Bed);
    R;
}

#' @export
intersect_intervals_really_fast_parallel <- function(P) {
    R={};
    if (length(P$Ast) ==0  || length(P$Bst)==0) {
      return(R);
    }
    R=intersect_intervals_really_fast(P$Ast, P$Aed, P$Bst,P$Bed);
    R;
}

#' @export
csv2bed <- function(csvfilename, bedfilename) {
	dat<- readcsv(filename);
	D=data.frame(chr=as.character(dat$chr), st=dat$st-1, ed=dat$ed);
	writebed(bedfilename, D, header=F);
}

#' Read BED File
#'
#' Reads a BED file and returns a data frame with columns chr, st, ed.
#'
#' @param filename Character. Path to BED file.
#'
#' @return Data frame with columns chr, st, ed.
#' @export
readbed <- function(filename)
{
  # Read the first 5 lines of the file
  readbed = readLines(filename, n=5)

  # Initialize the line counter
  i = 1

  # Loop through the first 5 lines to find the line starting with "chr"
  while (i < 5)
  {
    if (substr(readbed[[i]], 1, 3) == "chr") break
    i = i + 1
  }

  # If no line starting with "chr" is found, return an empty list and print a message
  if (i >= 5)
  {
    readbed = list()
    sprintf("%s is not a readable bed file.", filename)
    return
  }
  else
  {
    # Read the file from the line starting with "chr" and set column names
    readbed = read.csv(filename, sep="", skip=i-1, header=FALSE, colClasses=c("character", "numeric", "numeric"))
    names(readbed) <- c("chr", "st", "ed")

    # Increment the start positions by 1
    readbed$st = readbed$st + 1
  }

  # Return the data frame
  readbed
}


#' Convert BedGraph to RDensity Format
#'
#' Converts a BedGraph file to RDensity format for density calculation.
#'
#' @param inputfilename Character. Input BedGraph file name.
#' @param outputfilename Character. Output RDensity file name.
#' @param binsize Integer. Bin size for density calculation (default: 20).
#'
#' @return None.
#' @export
#' @examples
#' if (file.exists("input.bgr")) {
#'  bgr2rdensity("input.bgr", "output.rdensity", 20)
#'} else {
#'  cat("The file 'input.bgr' does not exist.\n")}
bgr2rdensity <- function(inputfilename, outputfilename=NULL, binsize=20) {

  # Print version information
  cat(sprintf("bgr2rdensity v.1.0.0    May 18, 2011\n"));

  # If output file name is not provided, generate it from input file name
  if (is.null(outputfilename)) {
    outputfilename= sub('.bgr','.rdensity', inputfilename);
  }

  # Read the .bgr file
  b1= readbgr(inputfilename);
  chrlist = levels(b1$chr);

  # Initialize the density list and binsize
  rd = list(density=vector("list", length(chrlist)),binsize=binsize);
  rd$read <- function(chr, addr) { rd$density[[chr]][1 + floor(addr/ binsize)]; };
  names(rd$density) <- chrlist;

  # Helper functions to get and inverse address
  getaddr <- function(x) { r= 1+floor((x-1)/binsize) };
  invaddr <- function(r) { r= 1+(r-1)*binsize };

  # Loop through each chromosome
  for (ch in chrlist) {

    # Select rows for the current chromosome
    sel = (b1$chr == ch);
    maxaddr= max(b1$ed[sel]);
    maxx = 1 + floor(maxaddr/ binsize);

    st = b1$st[sel];
    ed = b1$ed[sel];
    dd = b1$value[sel];
    rrst = getaddr(st);
    rred = getaddr(ed);
    mst = invaddr(rrst);
    med = invaddr(rred);
    bb=length(rrst);

    # Initialize the density array
    density.raw = vector("numeric", length=maxx);

    # Calculate the density
    for (i in 1:bb) {
      if (rred[i] > rrst[i]) {
        density.raw[rrst[i]] = density.raw[rrst[i]] + dd[i] * (mst[i] + binsize - st[i]);
        density.raw[rred[i]] = density.raw[rred[i]] + dd[i] * (ed[i] - med[i]+1);
        if (rred[i] > rrst[i]+1) {
          density.raw[(rrst[i]+1):(rred[i]-1)] = density.raw[(rrst[i]+1):(rred[i]-1)] + dd[i]*binsize;
        }
      } else {
        density.raw[rrst[i]] = density.raw[rrst[i]] + dd[i] * (ed[i] - st[i]+1);
      }
    }

    # Print the sum of the density
    cat(paste(ch,"-",bb,":",sum(density.raw)," ", sum(dd[1:bb] * (ed[1:bb]-st[1:bb]+1)),"\n"));

    # Store the density in the list
    rd$density[[ch]]= as.integer(round(density.raw/20));
  }

  # Save the density list to the output file
  save(rd, file=outputfilename, compress=TRUE);

}


#' @export
readcsv2 <- function(filename)
{
	library(tools)
	filebase = basename(filename);
	ext = tolower(file_ext(filebase));
	if (ext == "bed") {
		R=readbed(filename);
	} else {
		R=readcsv(filename);
	}
    R;
}

#' @export
readhotspot<-function(filename)
{
    readcsv(filename);
}

#' @export
readannot<-function(filename)
{
    read.csv(filename, sep=",", header=TRUE);
}


#' @export
extendHotspot <- function(hotspot, ext) {
	if (ext >=0) {
		hotspot$st = hotspot$st - ext;
		hotspot$ed = hotspot$ed + ext;
		st = sapply(hotspot$st, function(x) { ifelse(x>0,x,1);});
		hotspot$st = st;
	}
	hotspot;
}

#' @export
writebed <- function(filename, D, header=T)
{
    trackname=basename(filename);
    description= basename(filename);
    fid = file(filename, 'w');
    oheader= sprintf('browser hide all\ntrack name="%s" description="%s"\n', trackname, description);

	if (header) {
	    cat(paste(oheader),file=fid);
	}
    for (i in 1: length(D$st))
    {
        cat(sprintf('%s\t%d\t%d\n',as.character(D$chr[i]), D$st[i]-1, D$ed[i]),file=fid);
    }
    close(fid);
}


#' @export
writeBGRHeader <- function(fp, dataname, datadescription)
{


}



#' Write BedGraph Per Chromosome
#'
#' Writes BedGraph data for a single chromosome to a file.
#'
#' @param fp File handle to write to.
#' @param chr Character. Chromosome name.
#' @param count Numeric vector. Cut counts.
#' @param binsize Integer. Bin size (default: 20).
#' @param thr Numeric. Threshold (default: 2).
#'
#' @return None.
#' @export
writeBGRPerChrom <- function(fp, chr, count,  binsize = 20, thr=2)
{
  # Initialize variables
  span = length(count);
  sum = 0;
  c = 0;
  cc = 0;
  prevc = 0;
  meanc = 0;
  prevj = 0;
  sz = 0;
  MINTHR = 1.0;
  ii = 1;
  jj = 0;

  # Loop through the count array
  while (ii <= span) {
    c = count[ii];

    # If the count is greater than the threshold
    if (c > thr) {
      jj = ii + 1;

      # Find the end of the current segment with the same count
      while (count[jj] == c && jj <= span) {
        jj = jj + 1;
      }

      # Write the segment to the file
      if (c > 0 && jj <= span) {
        if (c < 30.0) {
          cat(sprintf("%s %d %d %.1f\n", chr, ii - 1, jj - 1, c), file = fp);
        } else {
          cat(sprintf("%s %d %d %.0f\n", chr, ii - 1, jj - 1, c), file = fp);
        }
      }
      ii = jj;
    } else {
      jj = ii;
      cc = 0;
      prevj = jj;
      prevc = 0;

      # Process segments with counts less than or equal to the threshold
      while (jj <= span) {
        if (count[jj] > thr)
          break;
        sum = 0;
        sz = 1;

        # Sum counts within the bin size
        while (jj <= span) {
          if (count[jj] > thr || sz > binsize)
            break;
          sum = sum + count[jj];
          jj = jj + 1;
          sz = sz + 1;
        }
        cc = cc + 1;
        meanc = sum / (sz - 1);

        # Write the previous segment to the file if the mean count changes
        if (cc > 1 && (round(meanc) != round(prevc))) {
          if (prevc > MINTHR && prevj <= span) {
            if (prevc < 30.0) {
              cat(sprintf("%s %d %d %.1f\n", chr, ii - 1, prevj - 1, prevc), file = fp);
            } else {
              cat(sprintf("%s %d %d %.0f\n", chr, ii - 1, prevj - 1, prevc), file = fp);
            }
          }
          ii = prevj;
        }
        prevj = jj;
        prevc = meanc;
      }

      # Write the last segment to the file
      if (meanc > MINTHR && jj <= span) {
        if (meanc < 30.0) {
          cat(sprintf("%s %d %d %.1f\n", chr, ii - 1, jj - 1, meanc), file = fp);
        } else {
          cat(sprintf("%s %d %d %.0f\n", chr, ii - 1, jj - 1, meanc), file = fp);
        }
      }
      ii = jj;
    }
  }
}



#@useDynLib bagfoot readcutcount
#' @useDynLib tiger readcutcount
readCutCount <- function(filename, chr) {
    maxaddr = 250000000;
    out <- .C("readcutcount", filename = as.character(datafilepath),chr= as.character(chr), span=as.integer(0), count= as.integer(vector("numeric",length=maxaddr)));
    cutcount <- out$count[1:out$span];
    cutcount;
}


#' @export
file_check <- function(filelist)  {
    res = file.exists(filelist);
    ox<-ifelse(res,"O","X");
    for (i in 1:length(ox)) {
        cat(sprintf('[%s] %s\n', ox[i], filelist[[i]]));
    }
    return (all(res));
}

#' @export
getMaxLocBAM<-function(bamFile) {   # returns the maximum locations of each chromosome
	dat = readBAMIndex(bamFile);
	ll = list();
	lapply(1:nrow(dat), function(ii) { ll[[as.character(dat$chr[ii])]] <<-  dat$maxloc[ii]; });
	ll;
}
