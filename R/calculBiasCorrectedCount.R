

#' Read Cut Sites Per Chromosome from BAM
#'
#' Extracts cut site positions from a BAM file for a specified chromosome, handling both paired-end and single-end sequencing. Optionally applies position shifts.
#'
#' @param bamfile Character. Path to the BAM file.
#' @param chrom Character. Chromosome name.
#' @param shifts Numeric vector of length 2. Position shifts for start and end (default: c(0, 0)).
#'
#' @return Numeric vector of cut site positions.
#' @export
#'
#' @examples
#' # readCutSitesPerChromFromBAM("sample.bam", "chr1")
readCutSitesPerChromFromBAM <- function(bamfile, chrom,shifts=c(0,0)) {

  #  cat('readCutSitesPerChromFromBAM\n');

  TAGLENGTHLIMIT=1000;

  param <- Rsamtools::ScanBamParam(what=c('rname', 'pos','qwidth','mpos','isize','strand'),flag=Rsamtools::scanBamFlag(isSecondMateRead=FALSE ),
                                   which = GenomicRanges::GRanges(chrom, IRanges::IRanges(1,536870912)));   # T

  bam <- data.frame(Rsamtools::scanBam(bamfile,param=param)[[1]]);
  bam=bam[abs(bam$isize) < TAGLENGTHLIMIT,];

  paired = any(bam$isize>0);
  if (paired) { # Paired-end sequencing
    st = bam$pos + ifelse(bam$strand=='+', 0, bam$qwidth + bam$isize);
    ed = bam$pos+  ifelse(bam$strand=='+', bam$isize-1, bam$qwidth-1);
  } else { # Single-end sequencing
    st = bam$pos;
    ed = bam$pos+bam$qwidth-1;

  }

  if (paired) {
    cuts = sort(c(st-1+shifts[1],ed+shifts[2]));
  } else {
    cuts= sort(ifelse(bam$strand == '+', st-1+shifts[1], ed+shifts[2]));

  }

  cuts;
}


# @useDynLib bagfoot intersect_intervals
#' @useDynLib bagfoot2 intersect_intervals
intersect_intervals_fast_index2C <- function(X1,X2,Y1,Y2)
{
  Isorted1= order(X1);
  Isorted2=  order(Y1);
  .Call("intersect_intervals", as.numeric(X1),as.numeric(X2),as.numeric(Y1),as.numeric(Y2),as.integer(Isorted1), as.integer(Isorted2));
}


#' @export
readFrequencyTable<-function(ftablefile) {
	cat(sprintf("reading %s ......\n", ftablefile));
	ftable = read.table(ftablefile);
	nr = nrow(ftable);
	ftable[nr,c("GPRatio","ObCutRatio")] = 0;
	ftable = cbind(ftable, Aj = ftable$ObCuts/ ftable$GenomicPositionCount);
	ftable$Aj[nr] = 0;
	# cat(sprintf("Num of cuts =  %d ......\n", sum(ftable$ObCuts)));
	ftable;
}



#' Fast Parallel Interval Intersection
#'
#' Finds intersections between two sets of intervals in parallel.
#'
#' @param P List. Contains Ast, Aed, Bst, Bed vectors.
#'
#' @return List of intersection indices.
#' @export
#'
#' @examples
#' # intersect_intervals_fast_parallel_index2(list(Ast = ..., Aed = ..., Bst = ..., Bed = ...))
intersect_intervals_fast_parallel_index2 <- function(P) {

  R={};
  if (length(P$Ast) ==0  || length(P$Bst)==0) {
    return(R);
  }
  R=intersect_intervals_fast_index2C(P$Ast, P$Aed, P$Bst,P$Bed);
  R;
}

#' Find Overlapping Indices in A (Parallel)
#'
#' Finds indices of intervals in A that overlap with intervals in B, processing by chromosome and optionally in parallel.
#'
#' @param A Data frame. Columns: 'chr', 'st', 'ed'.
#' @param B Data frame. Columns: 'chr', 'st', 'ed'.
#' @param mc.cores Integer. Number of cores for parallel processing.
#'
#' @return Integer vector of indices in A that overlap with B.
#' @export
#'
#' @examples
#' #A <- data.frame(chr = c("chr1", "chr1"), st = c(100, 200), ed = c(150, 250))
#' #B <- data.frame(chr = c("chr1", "chr1"), st = c(120, 220), ed = c(170, 270))
#' #OverlappingIndexAParallel(A, B, mc.cores = 2)
OverlappingIndexAParallel <- function(A, B, mc.cores=MCCORES) {
  # Function to process each chromosome
  OverlappingIndexAParallelPerChrom <- function(i) {
    result = NULL
    iA = 1:length(A$st)  # Indices for A
    iB = 1:length(B$st)  # Indices for B
    rangeA = iA[A$chr == chrs[i]]  # Indices of A for current chromosome
    rangeB = iB[B$chr == chrs[i]]  # Indices of B for current chromosome
    if ((length(rangeA) != 0) && (length(rangeB) != 0)) {
      P = list(Ast = A$st[rangeA], Aed = A$ed[rangeA], Bst = B$st[rangeB], Bed = B$ed[rangeB])  # Create list of start and end positions
      RR = intersect_intervals_fast_parallel_index2(P)  # Find overlapping intervals
      if (!is.null(RR)) {
        result = rangeA[RR]  # Get indices of overlapping intervals
      }
    }
    result
  }

  chrs = levels(factor(A$chr))  # Get unique chromosomes

  if (mc.cores == 1) {
    kk = lapply(seq_along(chrs), OverlappingIndexAParallelPerChrom)  # Process chromosomes sequentially
  } else {
    kk = parallel::mclapply(seq_along(chrs), OverlappingIndexAParallelPerChrom, mc.cores = mc.cores)  # Process chromosomes in parallel
  }

  unlist(kk)  # Combine results
}

#' @export
CompareHotspotsFaster<-function(A,B, mc.cores=MCCORES) {
	idx = OverlappingIndexAParallel(A,B, mc.cores=mc.cores);
	A[idx,];
}

#' @export
filesizenotzero <- function(filename) {
	if (file.exists(filename)) {
		r = file.info(filename)$size > 0;
	} else {
		r = FALSE;
	}
	r;
}



#' Read BedGraph File
#'
#' Reads a bedGraph file and returns a data frame with columns: chr, st, ed, value.
#'
#' @param filename Character. Path to the bedGraph file.
#'
#' @return Data frame with columns: chr, st, ed, value.
#' @export
#'
#' @examples
#' # readbgr2("sample.bedgraph")
readbgr2 <- function(filename)
{
    if (filename=="")
    {
        R=NULL;
        R;
    } else {
        # browser();
        R= readLines(filename, n=7);
        i=1;
        while (i < 7)
        {
            if (substr(R[[i]],1,3)=="chr") break;
            i=i+1;
        }
        if (i >= 7)
        {
            R={};
            sprintf("%s is not a readable bed file.", A);
            return;
        } else
        {
            #R=read.csv(filename, sep="", header=FALSE,comment.char="", skip=i-1, colClasses=c("factor","numeric","numeric","numeric"), stringsAsFactors=FALSE);
            if (length(grep(".gz", filename))==1) {
				freadfilename = paste('zcat', filename);
			} else {
				freadfilename = filename;
			}

            R=data.frame(data.table::fread(freadfilename,  header=FALSE, skip=i-1, colClasses=c("factor","numeric","numeric","numeric"), stringsAsFactors=FALSE));
            names(R) <- c("chr","st","ed","value");
            R$chr <- factor(R$chr);
            R$st = R$st+1;
            R;
        }
    }
}

#' @export
readcutcount3 <- function(cutcountfile) {

	cutcount = readbgr2(cutcountfile);
	cutcountPerChrom = split(cutcount, cutcount$chr);
	cutcount={};
	cutcountPerChrom;
}

#' @export
readBGRAsList <- readcutcount3


#' Convert Bias Correction Table
#'
#' Converts a bias correction table to a different n-mer size and writes the result to a new file.
#'
#' @param ftablefile Character. Path to the frequency table file.
#' @param np Integer. Desired n-mer size.
#'
#' @return Data frame with converted bias correction table.
#' @export
#'
#' @examples
#' # convertBiasCorrectionTable("Hexamer_FT_DS8497_mm9_withMap.txt", 4)
convertBiasCorrectionTable<- function(ftablefile, np) {

  #examples:
  #tab4=convertBiasCorrectionTable('Hexamer_FT_DS8497_mm9_withMap.txt', 4);
  #tab2=convertBiasCorrectionTable('Hexamer_FT_DS8497_mm9_withMap.txt', 2);

  ftable=read.table(ftablefile);
  #if (np==2 || np==4) {
  subseq= getNucleotideString(np);
  repeatN = function(x) { paste(rep('N',x),collapse='');};
  seqNseqN = paste(repeatN((6-np)/2),subseq,repeatN((6-np)/2), sep='');

  tab=data.frame(t(sapply(1:length(seqNseqN), function(ii)  {
    seq=seqNseqN[ii];
    sel=CalcSeqCodeForSeq(seq)
    subtable = ftable[sel,];
    w=apply(subtable,2,sum);
  })));
  row.names(tab) = subseq;
  tab$CorrectionFactor = tab$GPRatio / tab$ObCutRatio;
  tab = rbind(tab, ftable[nrow(ftable),]);

  if (np == 4) {
    outputfile = sub('Hexamer','Tetramer', ftablefile);
  } else if (np == 2) {
    outputfile = sub('Hexamer','Dimer', ftablefile);
  }
  if (outputfile == ftablefile) {
    outputfile = paste('NucBias_',np,'_',ftablefile, sep='');
  }
  write.table(tab, file=outputfile);
  tab;
}

#' @export
FindSeqCodeInBase2ContainingSeqInBase1<-function(code, base1=4, base2=6) {  # default
	n=base1; m=base2;

	if (code <1 | code > 4^n) {
		cat(sprintf('%d is out of range.\n', code));
	}

	diff = m-n;
	hd=diff/2;
	loop=seq(from=0, to=4^hd-1);
	as.vector(sapply(sapply(loop, function(x) x*4^(m-hd)+(code-1)*4^(hd) ), function(y) y+loop)+1);

}


#' @export
CalcSeqCodeForSeq<- function(seq) {

	np = nchar(toupper(seq));
	seqchar <- strsplit(seq, "")[[1]]
	nucleotide= c('A','C','G','T');

	code = c(0);
	mult = 1;
	for (ii in np:1) {
		ch = substring(seq,ii,ii);
		if (ch == 'A')  {
#			code = code * 4;  #no change
		} else if (ch =='C') {
			code = code  + 1 * mult;
		} else if (ch =='G') {
			code = code  + 2 * mult;
		} else if (ch =='T') {
			code = code  + 3 * mult;
		} else if (ch =='N') {
			code =  c(code, code+mult,code+2*mult,code+3*mult);
		} else {
			code = NA
		}
		mult = mult*4;
	}
	code = code+1;
	code;
}
