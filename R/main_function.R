# Fonction runGfoot
#' Run GFoot BAGFOOT2 Analysis
#'
#' Runs the main GFoot BAGFOOT2 analysis, generating footprinting results and plots.
#'
#' @param Gfoot GFoot object containing analysis parameters and data.
#' @param range Optional. Range of motifs to analyze.
#' @param graphout Logical. If TRUE, output graphs.
#' @param halfwidths Numeric vector. List of half-widths for aggregation (default: c(50)).
#' @param normalization Character. Normalization method (default: "default").
#' @param yrange Numeric vector. Y-axis range for plots (default: c(-2.2, 0.5)).
#' @param run_set Character. Optional run set identifier.
#'
#' @export

runGfoot<-function(Gfoot, range=NA, graphout=F, halfwidths=c(50),normalization="default", yrange=c(-2.2,0.5),run_set="") {
  print(paste("run_set:",run_set))
  drawMotifAggPlotOnMotifSetsForMultipleRangesAndWithComparisons(
    sitefiledir=Gfoot@motifDB@directory,
    cache= Gfoot@cachedir,
    outputdir=Gfoot@outputdir,
    datanamelist = list(Gfoot@control@name, Gfoot@treatment@name),
    cutcountdatalist = list(Gfoot@control@cutcount, Gfoot@treatment@cutcount),
    ncutcounts = c(Gfoot@control@count,Gfoot@treatment@count),
    normalization = normalization,
    sitefilelist = list(Gfoot@sitefile, Gfoot@sitefile),
    motiflistlist = list(Gfoot@motifDB@motiflistfile,Gfoot@motifDB@motiflistfile),
    ftablefile = Gfoot@gfootoption@biasfile,
    nuccodefilePattern = Gfoot@gfootoption@hexamer_pattern,
    halfwidthlist = halfwidths,
    comparisonIndices = list(c(1,2)),
    range=range,
    graphout=graphout,
    yrange=yrange,
    run_set= run_set
  );

  scatterDataOutGroups(
    outputdir=Gfoot@outputdir,
    datanamelist = list(Gfoot@control@name, Gfoot@treatment@name),
    motiflistlist =  list(Gfoot@motifDB@motiflistfile,Gfoot@motifDB@motiflistfile),
    sitefilelist = list(Gfoot@sitefile, Gfoot@sitefile),
    halfwidthlist =halfwidths,
    comparisonIndices = list(c(1,2)),
    range=range
  );
}




#' Draw Motif Aggregation Plots for Multiple Ranges and Comparisons
#'
#' Draws motif aggregation plots for multiple ranges and comparisons, saving results and plots.
#'
#' @param sitefiledir Character. Directory containing site files.
#' @param cache Character. Cache directory.
#' @param outputdir Character. Output directory.
#' @param datanamelist List of data names.
#' @param cutcountdatalist List of cut count data.
#' @param normalization Character. Normalization method.
#' @param ncutcounts Numeric vector. Number of cut counts per dataset.
#' @param sitefilelist List of site files.
#' @param motiflistlist List of motif list files.
#' @param ftablefile Character. Frequency table file.
#' @param nuccodefilePattern Character. Pattern for nucleotide code files.
#' @param halfwidthlist Numeric vector. List of half-widths.
#' @param comparisonIndices List of index pairs for comparison.
#' @param range Optional. Range of motifs to analyze.
#' @param yrange Numeric vector. Y-axis range for plots.
#' @param graphout Logical. If TRUE, output graphs.
#' @param run_set Character. Optional run set identifier.
#'
#' @return None. Writes output files and plots.
#' @export
#'
#' @examples
#' # drawMotifAggPlotOnMotifSetsForMultipleRangesAndWithComparisons(...)
drawMotifAggPlotOnMotifSetsForMultipleRangesAndWithComparisons <- function(
    sitefiledir='',
    cache='',
    outputdir='',
    datanamelist = list(),
    cutcountdatalist = list(),
    normalization="default",
    ncutcounts = NA,
    sitefilelist = list(),
    motiflistlist = list(),
    ftablefile = Fed_mm9_withMap,
    nuccodefilePattern = nuccode_hexamer_withMappability_35mer_mm9,
    halfwidthlist = c(50,100),
    comparisonIndices = list(c(1,2), c(3,2)),
    range= NA,
    yrange=c(-2.2,0.5),
    graphout=F,
    run_set="") {
  # Display a message indicating that the function has been called
  print("drawMotifAggPlotOnMotifSetsForMultipleRangesAndWithComparisons")
    # Record the start time of the execution
  ptm <- proc.time()
    # Change the working directory
  setwd(bg_dir)
    # Initialize the threshold to 0
  threshold = 0;
    # Save the base output directory
  outputdir_base = outputdir;
    # Display the current working directory and the frequency table file
  print(paste(getwd(),"ftablefile:",ftablefile))
    # Read the frequency table if it does not already exist

  # Read the frequency table only if ftablefile is not empty and file exists
  if (ftablefile != "" && file.exists(ftablefile)) {
    freqtable <- readFrequencyTable(ftablefile)
    np <- round(log(nrow(freqtable)-1)/log(4))
  } else {
    np <- 6
    freqtable <- NULL  # table by default with 1 everywhere
  }

    # Display the elapsed time since the start of execution
  proc.time() - ptm
  print(paste("1"))
    # Record the start time of the execution
  ptm <- proc.time()
    # Define the chromosomes to use
  chrs = paste("chr", c(1:22, "X","Y"),sep='');
    # Generate the nucleotide code filenames for each chromosome
  nuccodefiles= sapply(chrs,function(x) gsub("{chr}",x,nuccodefilePattern,fixed=T));
   # Check the existence of the nucleotide code files
  doexist=file.exists(nuccodefiles);
  chrs= chrs[doexist];
  nuccodefiles=nuccodefiles[doexist];

    # Load or generate the nucleotide codes for each chromosome
  if(!file.exists(paste0(bg_dir,run_set,"_nuccodes.Rdata"))){
    nuccodes= parallel::mclapply(nuccodefiles, function(nuccodefile) {
      nuccode=readNucleotideCodeForChromosomeForCuts(nuccodefile,np);
      nuccode; }, mc.cores=MCCORES,mc.preschedule=F);

    save(nuccodes,file = paste0(bg_dir,run_set,"_nuccodes.Rdata"))
  }else{
    load(paste0(bg_dir,run_set,"_nuccodes.Rdata"))
  }
    # Display the elapsed time since the start of execution
  print(proc.time() - ptm )
  print(paste("2"))

    # Record the start time of the execution
  ptm <- proc.time()

    # Normalize the break counts if necessary
  if (normalization=="default") {
    oldncutcounts= ncutcounts;
    ncutcounts = sapply(1:length(cutcountdatalist), function(kk) countCutcountOnSites(cutcountdatalist[[kk]],sitefilelist[[kk]]));
    print(paste("ncutcounts2:",ncutcounts))
  }

    # Display the elapsed time since the start of execution
  print(proc.time() - ptm )
  print(paste("3") )

    # Record the start time of the execution
  ptm <- proc.time()


    # Check that the break counts are not zero
  if (is.na(ncutcounts[1])) {
    stop("Error in drawMotifAggPlotOnMotifSetsForMultipleRangesAndWithComparisons: ncutcounts should not be zero.")
  }

    # For each specified half-width
  for (halfwidth in halfwidthlist) {
    ll <- lapply(seq_along(datanamelist), function(ii) {
      print(ii);
      dataname = datanamelist[[ii]];
      print(paste("dataname:",dataname))
      print(paste("ncutcounts[ii]:",ncutcounts[ii]))
            # Generate the output directory and filename
      out=getOutputDirAndFilename(outputdir_base, dataname, sitefilelist[[ii]], halfwidth);
      print(paste("out:",out))

            # Appelle la fonction pour dessiner le plot d'agr?gation de motifs
            run<-drawMotifAggPlotOnMotifSets(dataname = dataname,
                                       cutcount = cutcountdatalist[[ii]],
                                       ncutcount = ncutcounts[ii],
                                       hotspotfile = sitefilelist[[ii]],
                                       motiflist = read.table(motiflistlist[[ii]], stringsAsFactors=F),
                                       freqtablefile = ftablefile,
                                       nuccodes = nuccodes,
                                       sitefiledir = sitefiledir,
                                       cache = cache,
                                       outputdir = out$dir,
                                       halfwidth = halfwidth,
                                       logyrange = yrange,
                                       range = range,
                                       graphout = graphout);
      gc();
    });
  }

    # Display the elapsed time since the start of execution
  print(proc.time() - ptm)
  print(paste("5"))

  # For each pair of comparison indices
  for (index in comparisonIndices) {
    for (halfwidth in halfwidthlist) {

      i1 = index[1];
      i2 = index[2];

      dataname1 = datanamelist[[i1]];
      dataname2 = datanamelist[[i2]];
      motiflistfile1 = motiflistlist[[i1]];
      motiflistfile2 = motiflistlist[[i2]];
      out1 = getOutputDirAndFilename(outputdir_base, dataname1, sitefilelist[[i1]], halfwidth);
      out2 = getOutputDirAndFilename(outputdir_base, dataname2, sitefilelist[[i2]], halfwidth);

      title = sprintf('Footprinting Aggregation Plot (%s vs. %s)', dataname1, dataname2);

      datadir1 = sprintf("%s_%dbp", dataname1, halfwidth);
      datadir2 = sprintf("%s_%dbp", dataname2, halfwidth);
      MAoutputdir = file.path(outputdir, 'comparison');

      dlog(sprintf("\n\n\n\n\n\\nMAoutputdir=%s\n", MAoutputdir))

            # Call the function to draw the multiple aggregation plot
      MultipleAggregationPlot(
        dataname1 = dataname1,
        dataname2 = dataname2,
        motiflistfile1 = motiflistfile1,
        motiflistfile2 = motiflistfile2,
        outputfile1 = out1$filename,
        outputfile2 = out2$filename,
        title = title,
        threshold = threshold,
        datadir1 = datadir1,
        datadir2 = datadir2,
        outputdir = MAoutputdir,
        range = range,
        logyrange = yrange,
        graphout = graphout
      );
    }
  }
}


# Classe GFoot
GFoot <- setClass("GFoot",
                  representation(control = "CutCount",
                                 treatment = "CutCount",
                                 sitefile = "character",
                                 motifDB = "MotifDB",
                                 gfootoption = "GFootOption",
                                 outputdir = "character",
                                 cachedir = "character",
                                 outputfiles="character"),
                  prototype = list(
                    sitefile='',
                    motifDB = MotifDB(),
                    gfootoption= GFootOption(),
                    outputdir = getwd(),
                    cachedir = file.path(getwd(), "cache"),
                    outputfiles='')
)


#' Draw Motif Aggregation Plot on Motif Sets
#'
#' Draws motif aggregation plots for a single motif set, saving results and plots.
#'
#' @param dataname Character. Data name.
#' @param cutcount Cut count data.
#' @param ncutcount Numeric. Number of cut counts.
#' @param hotspotfile Character. Path to hotspot file.
#' @param motiflist Data frame. Motif list.
#' @param freqtablefile Character. Frequency table file.
#' @param nuccodes List. Nucleotide codes.
#' @param sitefiledir Character. Directory containing site files.
#' @param cache Character. Cache directory.
#' @param outputdir Character. Output directory.
#' @param halfwidth Integer. Half-width for aggregation.
#' @param strandwisePlot Logical. If TRUE, plot strandwise.
#' @param logyrange Numeric vector. Y-axis range for log ratio plot.
#' @param range Optional. Range of motifs to analyze.
#' @param graphout Logical. If TRUE, output graphs.
#'
#' @return Data frame with results for each motif.
#' @export
#'
#' @examples
#' # drawMotifAggPlotOnMotifSets(...)
drawMotifAggPlotOnMotifSets<- function(dataname = '',
                                       cutcount = '',
                                       ncutcount= 0,
                                       hotspotfile = '',
                                       motiflist = '',
                                       freqtablefile=NA,
                                       nuccodes = NA,
                                       sitefiledir='',
                                       cache='.',
                                       outputdir = '',
                                       halfwidth=50,
                                       strandwisePlot=F,
                                       logyrange=c(-2.2,0.5),
                                       range=NA,
                                       graphout=F
) {
  # Print the function name
  print("drawMotifAggPlotOnMotifSets")

  # Read the frequency table from the file ONLY if it exists and is not empty
  if (!is.na(freqtablefile) && freqtablefile != "" && file.exists(freqtablefile)) {
    ftable <- readFrequencyTable(freqtablefile)
  } else {
    ftable <- NA
  }

  # Create cache directory if it does not exist
  if (!file.exists(cache) && cache != '') {
    dir.create(cache, recursive=T);
  }

  # Set cache to current working directory if it is empty
  if (cache=='') {
    cache = getwd();
  }

  # Create output directory if it does not exist
  if (!file.exists(outputdir) && cache != '') {
    dir.create(outputdir, recursive=T);
  }

  # Set output directory to current working directory if it is empty
  if (outputdir=='') {
    outputdir = getwd();
  }

  # Record the start time
  begintime = Sys.time();

  # Set range to the length of motif list if it is not provided
  if (is.na(range[1])) {
    range=1:length(motiflist$motif);
  }

  # Initialize result list
  result = vector("list", length=length(motiflist$motif));

  # Stop execution if ncutcount is zero
  if (ncutcount == 0) {
    stop("Error in drawMotifAggPlotOnMotifSets: ncutcount should not be zero.")
  }

  # Calculate scale factor
  scalefactor = 100000000/ncutcount;
  cat(sprintf('num cutcount = %g,  scale factor = %g\n', ncutcount, scalefactor));


  # Function to run calculation for each motif
  runCalc<-function(ii, graphout=F) {
    outputfilename = make.names(sprintf('%s_%s_hw%dbp',motiflist$motif[ii], dataname,halfwidth));
    jpgfile1= file.path(outputdir,sprintf("aggregation_%s_observed_expected_cuts.jpg",outputfilename));

    cat(sprintf('No.%3d Motif:%s (%s %d bp)\n', ii, motiflist$motif[ii], motiflist$logo[ii], motiflist$motiflength[ii]));

    motif=list(
      name = as.character(motiflist$motif[ii]),
      sitefile=file.path(sitefiledir, motiflist$filename[ii]),
      motiflogo=as.character(motiflist$logo[ii]),
      motifcenter=motiflist$center[ii]);

    motifname= motiflist$motif[ii];
    HS=list(
      name='Hotspot',
      hotspotfile=hotspotfile,
      cutcountfile=cutcount,
      cutcountsignature=digest::digest(cutcount));

    ChIP=list(
      name='ChIP',
      hotspotfile=hotspotfile);

    AggregationPlotData= list(motifinfo=motif,
                              dhsinfo =HS,
                              chipinfo = ChIP,
                              ftablefile= freqtablefile,
                              ftable = ftable,
                              nuccodes =nuccodes);

    ss=DrawStamMotifAggPlotOnChIPBoundRegions(AggregationPlotData,
                                              output=outputfilename,
                                              halfwidth=halfwidth,
                                              title=sprintf('%s (%s)',dataname, motifname),
                                              yrange=logyrange,
                                              cache=cache,
                                              outputdir=outputdir,
                                              strandwisePlot=strandwisePlot,
                                              scalefactor = scalefactor,
                                              graphout=graphout);
    ss;
  }

  # Run calculations in parallel if MCCORES_MOTIF is greater than 1
  if (MCCORES_MOTIF == 1) {
    result<-lapply(range, function(x) { runCalc(x, graphout=graphout); });
  } else {
    result<-parallel::mclapply(range, function(x) { runCalc(x, graphout=graphout); }
                               ,mc.cores=MCCORES_MOTIF, mc.preschedule=F
    );
  }

  # Combine results and write to output file
  rout<-do.call("rbind", result);
  hotspotname = gsub(".csv","",basename(hotspotfile));
  outputfilename = file.path(outputdir, sprintf('../%s_on%s_%dbp_out.txt', dataname,hotspotname,  halfwidth));

  write.table(rout,file=outputfilename,sep='\t');
  rout;
}



#' @export
setMethod(f="run",		#GFoot::run
          signature="GFoot",
          definition=function(obj, range=NA, graphout=F, yrange=c(-2.2,0.5), mc.cores= MCCORES, run_set=NA) {
            #obj <- loadCutcount(obj);
            if (length(obj@control@cutcount)==0) {
              #	browser();
              obj@control<-readCutCountSites(obj@control, obj@sitefile);   #read cutcount
            }
            if (length(obj@treatment@cutcount)==0) {
              obj@treatment<-readCutCountSites(obj@treatment, obj@sitefile); #read cutcount
            }
            #obj@treatment<-readCutCountSites(obj@treatment, obj@sitefile); #read cutcount
            cat(sprintf("Num. of CPU cores to use: %d\n",mc.cores));
            MCCORES = mc.cores;
            obj@outputfiles=runGfoot(
              obj,              # GFoot obj
              range =  range,   # Range of CREB
              graphout = graphout,
              yrange= yrange,
              run_set= run_set
            );
            return(obj);
          })




#' @export
setGeneric(name="run",    #GFoot::run
           def = function(obj,range=NA,graphout=F,yrange=c(-2.2,0.5), mc.cores=MCCORES, run_set=NA) {
             standardGeneric("run");
           });

#' @export
setMethod(f="run",		#GFoot::run
          signature="GFoot",
          definition=function(obj, range=NA, graphout=F, yrange=c(-2.2,0.5), mc.cores= MCCORES, run_set=NA) {
            #obj <- loadCutcount(obj);
            if (length(obj@control@cutcount)==0) {
              #	browser();
              obj@control<-readCutCountSites(obj@control, obj@sitefile);   #read cutcount
            }
            if (length(obj@treatment@cutcount)==0) {
              obj@treatment<-readCutCountSites(obj@treatment, obj@sitefile); #read cutcount
            }
            #obj@treatment<-readCutCountSites(obj@treatment, obj@sitefile); #read cutcount
            cat(sprintf("Num. of CPU cores to use: %d\n",mc.cores));
            MCCORES = mc.cores;
            obj@outputfiles=runGfoot(
              obj,              # GFoot obj
              range =  range,   # Range of CREB
              graphout = graphout,
              yrange= yrange,
              run_set= run_set
            );
            return(obj);
          })


#' Read BedGraph File
#'
#' Reads a BedGraph file and returns a data frame with columns "chr", "st", "ed", and "value".
#'
#' @param filename Character. Path to the BedGraph file.
#'
#' @return Data frame with columns "chr", "st", "ed", "value".
#' @export
#'
#' @examples
#' #readbgr("example.bed")
readbgr <- function(filename)
{
  if (filename=="")
  {
    # If the filename is empty, return NULL
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
    # If no line starts with "chr", return an empty list and a message
    if (i >= 7)
    {
      R={};
      sprintf("%s is not a readable bed file.", A);
      return;
    } else
    {
      # Read the file starting from the line with "chr"
      R=read.csv(filename, sep="", header=FALSE,comment.char="", skip=i-1, colClasses=c("factor","numeric","numeric","numeric"), stringsAsFactors=FALSE);
      # Set column names
      names(R) <- c("chr","st","ed","value");
      # Increment the start position by 1
      R$st = R$st+1;
      R;
    }
  }
}


#' Read Hotspot File
#'
#' Reads a hotspot file in CSV format.
#'
#' @param filename Character. Path to the hotspot file.
#'
#' @return Data frame containing hotspot data.
#' @export
#'
#' @examples
#' #readhotspot("hotspot.csv")
readhotspot <- function(filename) {
  readcsv(filename)
}


#' Read Annotation File
#'
#' Reads an annotation file in CSV format.
#'
#' @param filename Character. Path to the annotation file.
#'
#' @return Data frame containing annotation data.
#' @export
#'
#' @examples
#' #readannot("annotation.csv")
readannot <- function(filename) {
  read.csv(filename, sep = ",", header = TRUE)
}

#' Read BED File as CSV
#'
#' Reads a BED file and returns a data frame with additional columns.
#'
#' @param filename Character. Path to the BED file.
#'
#' @return Data frame with columns "chr", "st", "ed", "MaxD", "AveD", "Zscore".
#' @export
#'
#' @examples
#' #readbed_as_csv("bedfile.bed")
readbed_as_csv <- function(filename) {
  hs = readbed(filename)
  lnth = length(hs$st)
  chs = data.frame(
    chr = hs$chr,
    st = hs$st + 1,
    ed = hs$ed,
    MaxD = vector("numeric", length = lnth),
    AveD = vector("numeric", length = lnth),
    Zscore = vector("numeric", length = lnth)
  )
  chs
}

#' Filter Hotspot Data
#'
#' Filters hotspot data by category and threshold.
#'
#' @param hotspot Data frame. Hotspot data.
#' @param category Character. Category to filter by.
#' @param threshold Numeric. Threshold value.
#'
#' @return Data frame containing filtered hotspot data.
#' @export
#'
#' @examples
#' #filterHotspot(hotspot, "MaxD", 5)
filterHotspot <- function(hotspot, category, threshold) {
  test = hotspot[[category]] >= threshold
  R = hotspot[test, ]
  R
}


#'@export
readCutCountOnSites <- function(cutcount_file = "", site_file = "") {
  cutcount = NULL
  sites = NULL

  reduced_bgr_filename = file.path(dirname(cutcount_file), sprintf("%s.bgr", paste(extractBaseName(cutcount_file), 'on', extractBaseName(site_file), sep = "_")))
  reduced_bgr_filename_zip = paste(reduced_bgr_filename, ".gz", sep = "")

  overlapping <- function(A, B) {
    if (as.character(A$chr[1]) != as.character(B$chr[1])) {
      return(NULL)
    } else {
      R = vector("logical", length = nrow(A))
      rangeA = 1:length(A$st)
      rangeB = 1:length(B$st)
      if ((length(rangeA) != 0) && (length(rangeB) != 0)) {
        P = list(Ast = A$st[rangeA], Aed = A$ed[rangeA], Bst = B$st[rangeB], Bed = B$ed[rangeB])
        RR = intersect_intervals_really_fast_parallel(P)
        I1 = unique(RR)
        R[rangeA[I1]] = TRUE
      }
      return(R)
    }
  }


#' @export
  filterCutcountByHotspot <- function(sites, chr) {
    cat(sprintf("reducing cutcount data for %s\n", chr))
    hotspot = sites[sites$chr == chr, c("chr", "st", "ed")]
    extendHotspot <- function(hotspot, extensionbp = 200) {
      hotspot$st = hotspot$st - extensionbp
      hotspot$ed = hotspot$ed + extensionbp
      hotspot = hotspot[hotspot$st > 0, ]
      return(hotspot)
    }
    if (nrow(hotspot) == 0) {
      return(NULL)
    } else {
      exh = extendHotspot(hotspot, extensionbp = 200)
      cutcountchr = cutcount[[chr]]
      if (is.null(cutcountchr)) {
        return(NULL)
      }
      ovlp = overlapping(cutcountchr, exh)
      subcutcount = cutcountchr[ovlp, ]
      return(subcutcount)
    }
  }

  if (file.exists(reduced_bgr_filename_zip)) {
    cat(sprintf("reading the cutcount data file: %s\n", reduced_bgr_filename_zip))
    cutcount = readcutcount2(reduced_bgr_filename_zip)
  } else {
    cat(sprintf("converting the cutcount data file: %s into %s\n", cutcount_file, reduced_bgr_filename))
    cutcount = readcutcount2(cutcount_file)
    sites = readhotspot(site_file)

    reduced_cutcount = parallel::mclapply(names(cutcount), function(x) {
      filterCutcountByHotspot(sites, x)
    }, mc.cores = MCCORES)
    names(reduced_cutcount) <- names(cutcount)

    for (ii in seq_along(reduced_cutcount)) {
      chr = names(reduced_cutcount)[[ii]]
      dat = reduced_cutcount[[ii]]
      if (!is.null(dat)) {
        cat(sprintf("saving cutcount data (%s) for %s\n", reduced_bgr_filename, chr))
        dat$st = dat$st - 1  # zero-based format
        if (ii == 1) {
          write.table(dat, file = reduced_bgr_filename, sep = " ", append = FALSE, quote = FALSE, row.names = FALSE, col.names = FALSE)
        } else {
          write.table(dat, file = reduced_bgr_filename, sep = " ", append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
        }
      }
    }
    system(sprintf('gzip %s', reduced_bgr_filename))  # gzip
    cutcount = reduced_cutcount
  }
  return(cutcount)
}



#' @export
quicksummaryListOrString<-function(obj) {
#	obj = cutcountfile

	if (is.character(obj)) {
		r=digest::digest(obj);
	} else {
		r=digest::digest(c(unlist(lapply(obj, nrow))));
	}
	r;
}

#' Draw Motif Aggregation Plot on ChIP Bound Regions
#'
#' Draws motif aggregation plots on ChIP-bound regions and saves results.
#'
#' @param S AggregationPlotData object.
#' @param tag Character. Optional tag for output.
#' @param motifdir Numeric vector. Motif directions.
#' @param output Character. Output file name.
#' @param halfwidth Integer. Half-width for aggregation.
#' @param movingaverage Logical. Apply smoothing.
#' @param plot1yrange Numeric. Y-axis range for observed/expected plot.
#' @param yrange Numeric. Y-axis range for log ratio plot.
#' @param title Character. Plot title.
#' @param cache Character. Cache directory.
#' @param outputdir Character. Output directory.
#' @param strandwisePlot Logical. Plot strandwise.
#' @param scalefactor Numeric. Normalization factor.
#' @param graphout Logical. Output graphs.
#' @export
DrawStamMotifAggPlotOnChIPBoundRegions <- function(S, tag="", motifdir=c(1,-1),output='',  halfwidth=20,movingaverage=F, plot1yrange=0,  yrange=c(-2.5,1.5),title='', cache='.', outputdir='.',strandwisePlot=F, scalefactor=1, graphout=F) {

  # Set various fields in the S object from its subfields
  S$cutcountname = S$dhsinfo$name;
  S$DHS = S$dhsinfo$hotspotfile;
  S$cutcount = S$dhsinfo$cutcountfile;
  S$cutcountsignature = S$dhsinfo$cutcountsignature
  S$CHIP = S$chipinfo$hotspotfile;
  S$chipname = S$chipinfo$name;
  S$motif=  S$motifinfo$sitefile;
  S$output = paste(S$output,tag, sep="");

  # Print the motif file path
  print(S$motif);

  # Create a save file path based on the cache directory and a hash of various inputs
  savefile = file.path(cache,sprintf('%s_%s_motifplotresult.dat',make.names(S$motifinfo$name), digest::digest(c(S$motif,S$CHIP,S$DHS,S$cutcountsignature,halfwidth))));

  # Function to read hotspot regions from a file
  readHotspotRegions <- function(filename) {
    dhs  =data.frame(data.table::fread(filename));
    names(dhs) <- c("ID","chr","st","ed", "MaxD","AveD","Zscore","pvalue");
    dhs = dhs[dhs$chr != "chrM",];   # eliminate chrM regions from the analysis
    dhs;
  }

  # Check if the save file exists or if debugging is enabled
  if (!file.exists( savefile) || BFOOT_DEBUG) {

    motif_chip_dhs=NULL;

    # Read the motif data and process it
    motif=data.frame(data.table::fread(S$motif))[,-1];
    istart = match('start',names(motif));
    if (istart == 2) {
      names(motif)[(istart-1):(istart+3)] = c("chr","st","ed", "dir","value");
      direction = motifdir;
      direction[direction==1] = '+';
      direction[direction==-1] = '-';
      motif = motif[is.element(motif$dir, direction),];
    } else {
      names(motif)=c("chr","st","ed", "value","dir");
      motif = motif[is.element(motif$dir, motifdir),];
    }

    motifwidth = motif$ed[1]-motif$st[1]+1;

    # Read the DHS regions
    dhs = readHotspotRegions(S$DHS);

    # Read the CHIP regions if different from DHS
    if (S$DHS == S$CHIP) {
      chip = dhs;
    } else {
      chip =data.frame(data.table::fread(S$CHIP));
      if (ncol(chip)>7) {
        names(chip)=  c("ID","chr","st","ed", "MaxD","AveD","Zscore","pvalue");
      } else if (ncol(chip)==3) {
        names(chip)=  c("chr","st","ed");
      }
    }

    result=NULL;
    if (nrow(motif)>0) {

      # Compare motif with CHIP regions
      motif_chip= CompareHotspotsFaster(motif, chip, mc.cores=MCCORES);
      if (S$DHS == S$CHIP) {
        motif_chip_dhs = motif_chip;
      } else {
        motif_chip_dhs= CompareHotspotsFaster(motif_chip, dhs, mc.cores=MCCORES);
      }
    }

    # Limit the number of motif sites if necessary
    if (nrow(motif_chip_dhs)> MAX_MOTIF_SITES) {
      set.seed(1);
      motif_chip_dhs = motif_chip_dhs[sample(nrow(motif_chip_dhs), MAX_MOTIF_SITES),];
    }

    # Generate the plot and save the result
    result=StamMotifAggPlot(motif_chip_dhs, S,output=output, halfwidth=halfwidth,movingaverage=movingaverage,yrange=yrange,plot1yrange=plot1yrange, title=title,cache=cache,outputdir=outputdir,strandwisePlot=strandwisePlot,scalefactor=scalefactor, graphout=graphout);
    save(result, file=savefile);
  } else {
    # Load the result from the save file if it exists
    cat(sprintf('loading %s\n', savefile));
    load(savefile);
  }
  result;
}   ###### DrawStamMotifAggPlotOnChIPBoundRegions



#' calcNormalizedExpectedObservedCounts
#'
#' This function calculates the normalized observed and expected cut counts around motif sites.
#' It is a core step for footprinting analysis.
#'
#' Main steps:
#'   - For each motif site, extract observed and expected cut counts in a window around the motif.
#'   - Normalize these counts using a scale factor.
#'   - Compute the log2 ratio of observed/expected counts (accessibility ratio) for the motif and flanking regions.
#'   - Summarize the accessibility signal inside the motif and in the flanking regions.
#'
#' Outputs:
#'   - A list containing:
#'       * mobs: mean observed cut counts profile
#'       * mexp: mean expected cut counts profile
#'       * logratio: log2(observed/expected) profile (accessibility ratio)
#'       * baselogratio: mean logratio in flanking regions (background accessibility)
#'       * motiflogratio: mean logratio inside the motif (footprint accessibility)
#'
#' These results are used to quantify motif accessibility and footprint depth for downstream analysis.
#'
#' @param sites Motif sites to analyze
#' @param S AggregationPlotData object with required data
#' @param cache Cache directory
#' @param basehalfwidth Base window size
#' @param halfwidth Window size around motif
#' @param scalefactor Normalization factor
#' @param motifbuffer Buffer size around motif
#' @param halfflankingregionwidth Flanking region size
#' @param heatmapout Output file for heatmap (optional)
#' @return List with normalized observed/expected counts and accessibility metrics
#' @export
calcNormalizedExpectedObservedCounts <- function (sites, S, cache, basehalfwidth, halfwidth, scalefactor, motifbuffer=2,halfflankingregionwidth=50, heatmapout='') {
  #print(paste("calcNormalizedExpectedObservedCounts sites:",nrow(sites)))
  getSubCount <- function(from, to) {
    subrange_st = motifcenter+((basehalfwidth)+from);
    subrange_ed = motifcenter+((basehalfwidth)+to);

    st = ifelse(counts$sites$dir==1,counts$region$st + subrange_st - 1,counts$region$ed - subrange_ed + 1)
    ed = ifelse(counts$sites$dir==1,counts$region$st + subrange_ed - 1,counts$region$ed - subrange_st + 1)

    isForward = counts$sites$dir==1;
    subobserved = counts$observed[,(subrange_st:subrange_ed)-1];

    subobserved[isForward,] = counts$observed[isForward,(subrange_st:subrange_ed)+1];
    #print(paste("subrange_st:",subrange_st,"subrange_ed:",subrange_ed))
    #browser()

    print(paste("subrange_st:", subrange_st, "subrange_ed:", subrange_ed))
    print((counts$expected))
    #print(paste("expected head:",head(expected) ))
    #print(paste("summary de expected:",    summary(as.vector(expected))))

    subexpected = counts$expected[,subrange_st:subrange_ed];


    subcounts=list(sites=counts$sites, dir=counts$sites$dir, region =  data.frame(chr = counts$region$chr, st = st,  ed = ed), observed=subobserved*scalefactor, expected=subexpected*scalefactor);

  }




# Function to calculate average aggregation of a heatmap
# average_aggregation
#
# @param heatmap
#
# @return
# @export
average_aggregation<-function(heatmap){
    sum1 = apply(heatmap, 1,sum)
    notnull = !is.na(sum1)
    num= sum(notnull)
    apply(heatmap[notnull,], 2, sum)/num
  }


# Function to calculate log ratio
calc_logratio<-function(cc) {
    rowNA1 = is.na(apply(cc$expected,1,sum))
    rowNA2 = is.na(apply(cc$observed,1,sum))
    mexp = average_aggregation(cc$expected[!rowNA1,])
    mobs = average_aggregation(cc$observed[!rowNA2,])
    mexp[mexp<0.0001] = 0.0001
    mobs[mobs<0.0001] = 0.0001
    logratio = log(mobs/mexp,2)
    list(mobs=mobs, mexp=mexp,logratio= logratio)
  }

  # Extracting necessary information from S
  cutcountfile=S$cutcount
  motifname= S$motifinfo$name
  motifcenter=S$motifinfo$motifcenter
  motiflogo= S$motifinfo$motiflogo
  motifwidth=nchar(motiflogo)
  cutcountname = S$cutcountname
  ftable= S$ftable
  ftablefile=S$ftablefile
  nuccodes= S$nuccodes

  nsite = nrow(sites)
  print(paste("nsite:",nsite))
  if (motifcenter==0) {
    motifcenter = floor(nchar(motiflogo)/2)+1
  }

  savefile = file.path(cache,sprintf('%s_%s_2.dat',make.names(S$motifinfo$name), digest::digest(c(sites,quicksummaryListOrString(cutcountfile) , S$dhsinfo$name, halfwidth, motifname,ftablefile))))

  if (!file.exists( savefile) || heatmapout != '' || BFOOT_DEBUG) {
    if (is.list(cutcountfile)) {
      cutcount = cutcountfile
    } else {
      cutcount = readcutcount3(cutcountfile)
    }

    counts=CalcObservedExpectedCuts(S, sites, cutcount,motifcenter=1, halfspan = basehalfwidth, nuccodes=nuccodes,motifwidth=motifwidth, ftable=ftable)
    cutcount = {}
    subcounts <- getSubCount(-halfwidth, halfwidth)
    FRcounts <- getSubCount(-motifbuffer-(halfflankingregionwidth-1)-motifcenter, motifwidth+motifbuffer+halfflankingregionwidth-motifcenter)
    flankingregionwidth= motifwidth+2*(motifbuffer+halfflankingregionwidth)
    flankingcolumn=c(1:halfflankingregionwidth, (flankingregionwidth-halfflankingregionwidth+1):flankingregionwidth)

    R2<- calc_logratio(FRcounts)
    baselogratio = mean(R2$logratio[flankingcolumn], trim=0.1)
    motiflogratio = mean(R2$logratio[-flankingcolumn], trim=0.1)

    Rfwd={}
    Rbwd={}
    R=calc_logratio(subcounts)
    result=list(Rboth=R, Rfwd=Rfwd, Rbwd=Rbwd, baselogratio=baselogratio, motiflogratio=motiflogratio)

    if (heatmapout!='') {
      heatmapfile1 = file.path(sprintf('%s.expected.csv',heatmapout))
      heatmapfile2 = file.path(sprintf('%s.observed.csv',heatmapout))
      cat(sprintf('saving %s\nsaving %s\n', heatmapfile1, heatmapfile2))
      write.table(subcounts$expected, file= heatmapfile1, sep=',', row.names=F, col.names=F)
      write.table(subcounts$observed, file= heatmapfile2,sep=',',  row.names=F, col.names=F)
    }
    save(result, file=savefile)
  } else {
    cat(sprintf('loading %s\n', savefile))
    load(savefile)
  }
  result
}

#' Draw Aggregated Observed/Expected Cut Count Plot
#'
#' Draws a plot of aggregated observed and expected cut counts around motif sites.
#'
#' @param mobs Numeric vector. Mean observed cut counts.
#' @param mexp Numeric vector. Mean expected cut counts.
#' @param motifcenter Integer. Motif center position.
#' @param halfwidth Integer. Half-width for aggregation.
#' @param title Character. Plot title.
#' @param motifwidth Integer. Motif width.
#' @param linea Integer. Start of motif region.
#' @param lineb Integer. End of motif region.
#' @param motiflogo Character. Motif logo.
#' @param plot1yrange Numeric. Y-axis range for plot.
#'
#' @return None.
#' @export
drawAggregatedObservedExpectedCutcountPlot <- function (mobs, mexp, motifcenter, halfwidth, title,  motifwidth, linea, lineb,motiflogo,  plot1yrange) {

  # Define color codes for the plot lines
  colorcode=c('red','darkolivegreen','darkblue');

  # Determine the minimum and maximum values for the y-axis
  miny=min(mobs,mexp);
  maxy=max(mobs,mexp);

  # Calculate the x-axis range based on the motif logo and center
  xr =  seq(from=1, to= nchar(motiflogo)) - motifcenter;

  # Set the y-axis range based on the provided plot1yrange or calculated min and max values
  if (length(plot1yrange)<2) {
    yrange = c(miny,maxy);
  } else {
    yrange= plot1yrange;
  }

  hgt = 0;

  # Set the plot margins
  op <- par(mar=c(5,6,4,2) + 0.1);

  # Initialize the plot with the specified parameters
  plot((-halfwidth:halfwidth), mobs, "n", main = title, col = colorcode[1], ylim=c(miny,maxy),
       xlab = '(BP)',  ylab='', lwd=2, cex.axis=1.5, cex.lab=1.5, cex.sub=1.5,
       xaxs='i',yaxs='i', las=1);

  # Add a y-axis label
  mtext('Observed vs. Expected\nCounts', side=2, line=3, cex=1.5);

  # Add the motif logo text to the plot
  text(xr, rep(maxy- (maxy-miny)*0.05, length(xr)), labels= unlist(strsplit(motiflogo,'')),cex=0.9);

  # Plot the observed counts line
  lines((-halfwidth:halfwidth), mobs, col= colorcode[1],lwd=2.5,lty=1);    #observed

  # Plot the expected counts line
  lines((-halfwidth:halfwidth), mexp, col= colorcode[2],lwd=2.5,lty=1);    #expected

  # Add vertical lines to indicate the motif region if motifwidth is greater than 0
  if (motifwidth >0) {
    lines(c(linea, linea), yrange, col='green', lwd=1.5,lty=2);
    lines(c(lineb, lineb), yrange, col='green', lwd=1.5,lty=2);
  }

  # Add a legend to the plot
  legend(-halfwidth+1, y = yrange[2] - (yrange[2]-yrange[1])*0.05, bty='n',c("Observed", "Expected"), lwd=3, col = colorcode,cex=1.3);

}

#' Draw Log Ratio Plot
#'
#' Draws a plot of log2(observed/expected) cut counts around motif sites.
#'
#' @param logratio Numeric vector. Log ratio profile.
#' @param yrange Numeric vector. Y-axis range for plot.
#' @param halfwidth Integer. Half-width for aggregation.
#' @param title Character. Plot title.
#' @param motiflogratio Numeric. Mean log ratio inside motif.
#' @param baselogratio Numeric. Mean log ratio in flanking regions.
#' @param linea Integer. Start of motif region.
#' @param lineb Integer. End of motif region.
#' @param motiflogo Character. Motif logo.
#' @param motifcenter Integer. Motif center position.
#'
#' @return None.
#' @export
drawLogratioPlot <- function (logratio, yrange, halfwidth, title, motiflogratio,baselogratio, linea, lineb, motiflogo, motifcenter) {

  # Determine the minimum and maximum values for the y-axis
  miny=min(logratio);
  maxy=max(logratio);

  # Set the y-axis range based on the provided yrange or calculated min and max values
  if (length(yrange)<2) {
    yrange = c(miny,maxy);
  }

  # Set the plot margins
  op <- par(mar=c(5,6,4,2) + 0.1);

  linewidth= 3;

  # Initialize the plot with the specified parameters
  plot((-halfwidth:halfwidth), logratio,"n", main = title,
       xlab = '(BP)', ylab='', ylim=yrange, lwd=2, cex.axis=1.5, cex.lab=1.5, cex.sub=1.5,
       xaxs='i',yaxs='i',las=1);

  # Add baseline and mean lines
  lines(c(-halfwidth, halfwidth), c(motiflogratio,motiflogratio), col= 'darkred', lwd=2, lty=2);
  lines(c(-halfwidth, halfwidth), c(baselogratio,baselogratio), col= 'black', lwd=2, lty=3);

  # Plot the log ratio line
  lines((-halfwidth:halfwidth), logratio,"l", col = 'darkblue',lwd=linewidth);

  # Add a y-axis label
  mtext('Log Ratio',side=2, line=3,cex=1.5);
  par(op);

  # Add vertical lines to indicate the motif region
  lines(c(linea, linea), yrange, col='green', lwd=1.5,lty=2);
  lines(c(lineb, lineb), yrange, col='green', lwd=1.5,lty=2);

  # Calculate the x-axis range based on the motif logo and center
  xr =  seq(from=1, to= nchar(motiflogo)) - motifcenter;

  # Add the motif logo text to the plot
  text(xr, rep(maxy, length(xr)), labels= unlist(strsplit(motiflogo,'')),cex=0.9);

}




#' Stam Motif Aggregation Plot
#'
#' Performs core calculation for motif footprinting analysis, computing observed/expected cut counts and footprint depth.
#'
#' @param sites Data frame. Motif sites to analyze.
#' @param S AggregationPlotData object.
#' @param output Character. Output file name base.
#' @param halfwidth Integer. Half-width for aggregation.
#' @param movingaverage Logical. Apply smoothing.
#' @param yrange Numeric. Y-axis range for log ratio plot.
#' @param plot1yrange Numeric. Y-axis range for observed/expected plot.
#' @param title Character. Plot title.
#' @param motifoffset Integer. Offset for motif region.
#' @param trim Numeric. Trim parameter for mean calculation.
#' @param strandwisePlot Logical. Plot strandwise.
#' @param basehalfwidth Integer. Base window size.
#' @param cache Character. Cache directory.
#' @param outputdir Character. Output directory.
#' @param scalefactor Numeric. Normalization factor.
#' @param ttest Logical. Perform statistical test.
#' @param graphout Logical. Output graphs.
#'
#' @return List with footprinting results for the motif.
#' @export
StamMotifAggPlot<- function(sites, S,output='', halfwidth=20,movingaverage=F, yrange=0,plot1yrange=0, title='', motifoffset=2, trim=0, strandwisePlot=F,basehalfwidth=200, cache='.',outputdir='.',scalefactor=1,ttest=F, graphout=F) {
  # Print the number of sites
  #print(paste("StamMotifAggPlot sites:",nrow(sites)))

  # Extract motif information
  motifcenter=S$motifinfo$motifcenter;
  motiflogo= S$motifinfo$motiflogo;
  motifwidth=nchar(motiflogo);

  # Create output directory if it doesn't exist
  if (!file.exists(outputdir) && outputdir != '') {
    dir.create(outputdir, recursive=T);
  }

  # Define the data file path
  datafile = file.path(outputdir,sprintf("aggregation_%s.dat",output));

  # Check if data file exists and is not empty or if debugging is enabled
  if (!filesizenotzero(datafile) || BFOOT_DEBUG) {

    # Calculate normalized expected and observed counts
    RR <- calcNormalizedExpectedObservedCounts(sites, S, cache, basehalfwidth, halfwidth, scalefactor)

    # Extract results
    mobs = RR$Rboth$mobs;
    mexp = RR$Rboth$mexp;
    logratio = RR$Rboth$logratio;
    baselogratio=RR$baselogratio;
    motiflogratio=RR$motiflogratio;
    footprintdepth=motiflogratio-baselogratio; # calculating FPD

    # Define motif region
    xr =  seq(from=1, to= nchar(motiflogo)) - motifcenter;
    x1 = min(xr); x2 = max(xr);

    linea = x1 - motifoffset;
    lineb = x2 + motifoffset;

    halfpoint= floor(length(logratio)/2)+1;
    motifregion=seq(from=halfpoint+x1-motifoffset, to=halfpoint+x2+motifoffset);
    motifinside = rep(F, length(logratio));
    motifinside[seq(from=halfpoint+x1-motifoffset, to=halfpoint+x2+motifoffset)] = T;

    # Perform t-test if required
    if (ttest) {
      res<-wilcox.test(logratio~motifinside,  alternative="greater");

      mu = 0;
      xbar = mean(logratio[motifinside],trim=0.1);
      nn = sum(motifinside);
      ssd= sd(logratio);
      ttt= (xbar-mu) / (ssd/sqrt(nn));

      alpha = 0.05;
      t.alpha = qt(1-alpha, df=nn-1);

      ttestoutputfile = file.path(outputdir,sprintf('%s_ttest.txt', output));
      sink(ttestoutputfile);
      cat(sprintf('t=%f\n', ttt));
      cat(sprintf('critical value at %f = %f\n', alpha, t.alpha));
      cat(sprintf('Univariate t-test pvalue = %f\n', pt(ttt, df=nn-1)));
      print(res);
      sink();
    }

    # Calculate trimmed mean and save results
    trimmean=mean(logratio[motifregion],trim=trim);
    save(mobs, mexp, logratio,baselogratio,motiflogratio, footprintdepth,linea,lineb, file=datafile);

  } else {
    # Load existing data file
    cat(sprintf("loading datafile:%s\n", datafile));
    load(datafile);
  }

  # Generate graphs if required
  if (graphout) {
    if (length(yrange)<2) {
      rangestring = "";
    } else {
      rangestring = paste(yrange, collapse="to");
    }
    jpgfile1= file.path(outputdir,sprintf("aggregation_%s_observed_expected_cuts%s.jpg",output,rangestring));
    jpgfile2= file.path(outputdir,sprintf("log_ratio_%s_observed_expected_cuts%s.jpg",output,rangestring));
    pdffilename1 =  file.path(outputdir,sprintf("aggregation_%s_observed_expected_cuts%s.pdf",output,rangestring));
    pdffilename2 =  file.path(outputdir,sprintf("log_ratio_%s_observed_expected_cuts%s.pdf",output,rangestring));

    # Plot observed vs expected cut counts
    jpeg(jpgfile1,width = 600, height = 600);
    drawAggregatedObservedExpectedCutcountPlot(mobs, mexp, motifcenter, halfwidth, title,  motifwidth, linea, lineb,motiflogo,  plot1yrange);
    dev.off();

    # Plot log ratio
    jpeg(jpgfile2,width = 600, height = 600);
    drawLogratioPlot(logratio, yrange, halfwidth, title, motiflogratio,baselogratio,  linea, lineb, motiflogo, motifcenter);
    dev.off();
  }

  # Return result as a list
  result = list(name= title, numSites=nrow(sites), meanlog=motiflogratio, baselog=baselogratio, footprintdepth=footprintdepth, graphdatafile=datafile);
  result;
}

#' Calculate Observed and Expected Cuts
#'
#' Calculates observed and expected cut counts for motif sites.
#'
#' @param S AggregationPlotData object.
#' @param motif Data frame. Motif sites.
#' @param tdensity List. Cut count data per chromosome.
#' @param motifcenter Integer. Motif center position.
#' @param halfspan Integer. Window size around motif.
#' @param nuccodes List. Nucleotide codes.
#' @param motifwidth Integer. Motif width.
#' @param ftable Data frame. Frequency table.
#'
#' @return List with observed and expected cut counts.
#' @export
CalcObservedExpectedCuts <- function(S, motif, tdensity,motifcenter=0, halfspan = 50, nuccodes=NA,motifwidth=0, ftable=NA) {
  # Calculate the midpoint of the motif based on its direction and center
  if (motifcenter==0) {
    midpoint= ifelse(motif$dir==1,motif$st + floor(motifwidth/2), motif$ed-floor(motifwidth/2));
  } else {
    midpoint= ifelse(motif$dir==1, motif$st + motifcenter -1, motif$ed-(motifcenter-1));
  }

  # Define the target region around the midpoint
  tregion = data.frame(chr= motif$chr, st = midpoint - halfspan + ifelse(motif$dir==-1,-1,0), ed = midpoint + halfspan+ ifelse(motif$dir==-1,-1,0));  #AP1

  # Initialize a logical vector to check if regions are within valid chromosome ranges
  chroms = names(tdensity);
  isValid = vector("logical", length=nrow(tregion));

  # Check if the end of the target region is within the maximum coordinate of the chromosome
  for (chr in chroms) {
    maxcoord = max(tdensity[[chr]]$ed);
    isValid[tregion$chr==chr] = tregion$ed[tregion$chr==chr] <= maxcoord;
  };

  # Filter out invalid motifs and regions
  motif = motif[isValid,];
  tregion = tregion[isValid,];

  # Calculate observed cuts in the target region
  ObservedCuts = CalcHeatmapFastMatrixOutputC(tregion, tdensity, dir=motif$dir);
  observedsum= apply(ObservedCuts, 1, sum);

  # Calculate expected rates of cuts
  # bfoot version
 # ExpectedRates = CalcExpectedRates(tregion, dir=motif$dir,nuccodes=nuccodes,ftable=ftable);

  # Initialize the matrix for expected cuts
  #ExpectedCuts = matrix(0, nrow=nrow(ExpectedRates), ncol = ncol(ExpectedRates));

  # Calculate expected cuts based on observed sums and expected rates
  #for (ii in 1:nrow(ExpectedCuts)) {
  #  ExpectedCuts[ii,] = observedsum[ii]* ExpectedRates[ii,];
  #}

#tiger version
  # Gestion du expected : si pas de bias_table, bypass en mettant à 1
  if (is.null(ftable) || length(ftable) == 0 || all(is.na(ftable))) {
    ExpectedCuts <- matrix(1, nrow = nrow(ObservedCuts), ncol = ncol(ObservedCuts))
  } else {
    # Calculate expected rates of cuts (avec bias)
    ExpectedRates <- CalcExpectedRates(tregion, dir=motif$dir, nuccodes=nuccodes, ftable=ftable)
    # Initialize the matrix for expected cuts
    ExpectedCuts <- matrix(0, nrow=nrow(ExpectedRates), ncol=ncol(ExpectedRates))
    # Calculate expected cuts based on observed sums and expected rates
    for (ii in 1:nrow(ExpectedCuts)) {
      ExpectedCuts[ii, ] <- observedsum[ii] * ExpectedRates[ii, ]
    }
  }
#####
  # Return the results as a list
  list(sites=motif, region=tregion, observed=ObservedCuts, expected=ExpectedCuts);
}  #CalcAggregationMotif


#' Multiple Aggregation Plot
#'
#' Draws multiple aggregation plots for motif comparisons and saves results.
#'
#' @param dataname1 Character. Name of first dataset.
#' @param dataname2 Character. Name of second dataset.
#' @param motiflistfile1 Character. Path to first motif list file.
#' @param motiflistfile2 Character. Path to second motif list file.
#' @param outputfile1 Character. Path to first output data file.
#' @param outputfile2 Character. Path to second output data file.
#' @param title Character. Plot title.
#' @param threshold Numeric. Threshold for filtering.
#' @param datadir1 Character. Directory for first dataset.
#' @param datadir2 Character. Directory for second dataset.
#' @param outputdir Character. Output directory.
#' @param scalefactor Numeric vector. Scale factors for data.
#' @param range Optional. Range of motifs to analyze.
#' @param logyrange Numeric vector. Y-axis range for log ratio plot.
#' @param graphout Logical. Output graphs.
#'
#' @return None.
#' @export

MultipleAggregationPlot <- function(dataname1="", dataname2="", motiflistfile1="", motiflistfile2="",
                                    outputfile1="",outputfile2="",title = "",threshold=0,datadir1='',datadir2='',
                                    outputdir=' ',scalefactor= c(1,1),range=NA, logyrange=0, graphout=F) {

  # Create output directory if it does not exist
  if (!file.exists(outputdir) && outputdir != '') {
    dir.create(outputdir, recursive=T);
  }

  # Read motif lists from files
  motiflist1=read.table(motiflistfile1,stringsAsFactors=F);
  motiflist2=read.table(motiflistfile2,stringsAsFactors=F);

  # Extract motifs from the lists
  mfs1 = as.character(motiflist1$motif);
  mfs2 = as.character(motiflist2$motif);

  # Find common motifs between the two lists
  commonmotifs = intersect(mfs1,mfs2);

  # Filter motif lists to keep only common motifs
  common1= motiflist1$motif %in% commonmotifs;
  common2= motiflist2$motif %in% commonmotifs;

  motiflist1 = subset(motiflist1,common1);
  motiflist2 = subset(motiflist2,common2);

  # Read output data files
  printf("reading %s..\n", outputfile1);
  printf("reading %s..\n", outputfile2);

  output1 = read.table(outputfile1,sep='\t');
  output2 = read.table(outputfile2,sep='\t');

  # Filter output data to keep only common motifs
  output1 = output1[common1,];
  output2 = output2[common2,];

  # Filter data based on threshold or range
  if (is.na(range[1])) {
    selection= (output1$numSites>=threshold) | (output2$numSites>=threshold);
    output1 = output1[selection,];
    output2 = output2[selection,];
    motiflist1 = motiflist1[selection,];
    motiflist2 = motiflist2[selection,];
    motiflist = motiflist1;
  } else {
    motiflist1 = motiflist1[range,];
    motiflist2 = motiflist2[range,];
    motiflist = motiflist1;
  }


  # Function to draw double aggregated observed expected cut count plot
#' @export
  drawDoubleAggregatedObservedExpectedCutcountPlot <- function (mobs1, mobs2,motifcenter,
                                                                title, motiflogo,  plot1yrange=0) {

    halfwidth= (length(mobs1)-1)/2;
    motifwidth = nchar(motiflogo);
    colorcode=c('darkolivegreen','purple');

    miny=min(mobs1,mobs2);
    miny=0;

    maxy=max(mobs1,mobs2);

    xr =  seq(from=1, to= nchar(motiflogo)) - motifcenter;
    if (length(plot1yrange)<2) {
      yrange = c(miny,maxy);
    } else {
      yrange= plot1yrange;
    };

    hgt = 0;

    op <- par(mar=c(5,6,4,2) + 0.1);

    xr =  seq(from=1, to= nchar(motiflogo)) - motifcenter;
    x1 = min(xr); x2 = max(xr);

    linea = x1 - motifoffset;
    lineb = x2 + motifoffset;

    plot((-halfwidth:halfwidth), mobs1,"n", main = title, col = colorcode[1], ylim=c(miny,maxy),
         xlab = '(BP)',  ylab='', lwd=2, cex.axis=1.5, cex.lab=1.5, cex.sub=1.5,
         xaxs='i',yaxs='i', las=1);
    mtext('Normalized Cut Counts', side=2, line=3, cex=1.5);

    if (length(mobs1)!=length(mobs2)) {
      print("length of mobs1 and mobs2 differ.");
    }
    lines((-halfwidth:halfwidth), mobs1, col= colorcode[1],lwd=2.5,lty=1);
    lines((-halfwidth:halfwidth), mobs2, col= colorcode[2],lwd=2.5,lty=1);

    if (motifwidth >0) {
      lines(c(linea, linea), yrange, col='green', lwd=1.5,lty=2);
      lines(c(lineb, lineb), yrange, col='green', lwd=1.5,lty=2);
    }

    legend(-halfwidth+1, y = yrange[2] - (yrange[2]-yrange[1])*0.05, bty='n',c(dataname1, dataname2), lwd=3, col = colorcode,cex=1.3);
  }

  # Function to draw double log ratio plot
  drawDoubleLogratioPlot <- function (logratio1, logratio2, meanlog1, meanlog2, motifcenter, title, motiflogo, yrange=0) {

    motifwidth = nchar(motiflogo);
    colorcode=c('darkolivegreen','purple');

    miny=min(min(logratio1),min(logratio2));
    maxy=max(max(logratio1),max(logratio2));

    xr =  seq(from=1, to= nchar(motiflogo)) - motifcenter;
    x1 = min(xr); x2 = max(xr);

    linea = x1 - motifoffset;
    lineb = x2 + motifoffset;

    if (length(yrange)<2) {
      yrange = c(miny,maxy);
    }

    op <- par(mar=c(5,6,4,2) + 0.1);

    linewidth= 3;
    plot((-halfwidth:halfwidth), logratio1,"n", main = title,
         xlab = '(BP)', ylab='', ylim=yrange, lwd=2, cex.axis=1.5, cex.lab=1.5, cex.sub=1.5,
         xaxs='i',yaxs='i',las=1);

    lines(c(-halfwidth, halfwidth), c(meanlog1,meanlog1), col= colorcode[1], lwd=2, lty=2);
    lines(c(-halfwidth, halfwidth), c(meanlog2,meanlog2), col= colorcode[2], lwd=2, lty=2);

    lines(c(-halfwidth, halfwidth), c(0,0), col= 'black', lwd=1.5, lty=1);

    lines((-halfwidth:halfwidth), logratio1,"l", col = colorcode[1],lwd=linewidth);
    lines((-halfwidth:halfwidth), logratio2,"l", col = colorcode[2],lwd=linewidth);

    mtext('Footprinting Depth',side=2, line=3,cex=1.5);
    par(op);

    lines(c(linea, linea), yrange, col='green', lwd=1.5,lty=2);
    lines(c(lineb, lineb), yrange, col='green', lwd=1.5,lty=2);

    xr =  seq(from=1, to= nchar(motiflogo)) - motifcenter;
    text(xr, rep(maxy*0.95, length(xr)), labels= unlist(strsplit(motiflogo,'')),cex=0.9);

    legend(-halfwidth+1, y = yrange[1] + (yrange[2]-yrange[1])*0.15, bty='n',c(dataname1, dataname2), lwd=3, col = colorcode,cex=1.3);
  }

  # Set range if not provided
  if (is.na(range[1])) {
    range = 1:nrow(motiflist);
  }

  # Loop through each motif in the range
  for (ii in 1:length(range)) {
    motif = motiflist[ii,];
    name = as.character(motif$motif);
    print(name);

    sites1 = output1[ii,];
    sites2 = output2[ii,];

    numsite1 = sites1$numSites;
    numsite2 = sites2$numSites;

    name = as.character(motif$motif);
    motiflogo = as.character(motif$logo);
    motifcenter =motif$center;

    load(as.character(output1$graphdatafile[ii]));
    mobs1 = mobs;
    mexp1 = mexp;
    logratio1 =  logratio - baselogratio;
    meanlog1 = motiflogratio - baselogratio;

    load(as.character(output2$graphdatafile[ii]));
    mobs2 = mobs;
    mexp2 = mexp;
    logratio2 =  logratio - baselogratio;
    meanlog2 = motiflogratio - baselogratio;

    motifoffset=2;
    plot1yrange=0;

    halfwidth= (length(logratio1)-1)/2;
    filenamebase = make.names(name);

    if (length(logyrange)<2) {
      rangestring = "";
    } else {
      rangestring = paste(logyrange, collapse="to");
    }

    jpgfile1=file.path(outputdir,sprintf("aggregation_%s_%s_vs_%s_%dbp%s.jpg",filenamebase,dataname1, dataname2,halfwidth,rangestring));
    jpgfile2=file.path(outputdir,sprintf("logratio_%s_%s_vs_%s_%dbp%s.jpg",filenamebase,dataname1, dataname2,halfwidth,rangestring));

    pdffilename1 =  file.path(outputdir,sprintf("aggregation_%s_%s_vs_%s_%dbp%s.pdf",filenamebase,dataname1, dataname2,halfwidth,rangestring));
    pdffilename2 =  file.path(outputdir,sprintf("logratio_%s_%s_vs_%s_%dbp%s.pdf",filenamebase,dataname1, dataname2,halfwidth,rangestring));

    if (graphout) {
      jpeg(jpgfile1,width = 600, height = 600);
      drawDoubleAggregatedObservedExpectedCutcountPlot(mobs1, mobs2,motifcenter,name, motiflogo);
      dev.off();

      jpeg(jpgfile2,width = 600, height = 600);
      drawDoubleLogratioPlot(logratio1, logratio2, meanlog1, meanlog2, motifcenter, title, motiflogo, yrange=logyrange);
      dev.off();
    }
  }
}
