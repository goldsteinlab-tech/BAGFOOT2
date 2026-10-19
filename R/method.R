# Definition of methods for Gfoot objects
#' Read Cut Count Sites (S4 Generic)
#'
#' Generic method to read cut count data for specified sites.
#'
#' @param obj S4 object of class CutCount.
#' @param site_filename Character. Path to the site file.
#'
#' @return S4 object of class CutCount with cutcount slot populated.
#' @export
#'
#' @examples
#' # readCutCountSites(obj, "sites.txt")
setGeneric(name="readCutCountSites",
           def = function(obj, site_filename) {
             standardGeneric("readCutCountSites");
           });
#' @describeIn readCutCountSites Method for CutCount objects.
setMethod(f="readCutCountSites",
          signature="CutCount",
          definition=function(obj,site_filename) {
            if (length(obj@cutcount)==0) {
              cat(sprintf("reading %s ...\n", obj@file));
              obj@cutcount <- readCutCountOnSites(cutcount_file = obj@file,site_file = site_filename);
            }
            return(obj);
          });

#' Load Cut Count Data for GFoot Object (S4 Generic)
#'
#' Loads cut count data for control and treatment in a GFoot object.
#'
#' @param obj S4 object of class GFoot.
#'
#' @return S4 object of class GFoot with cutcount slots populated.
#' @export
#'
#' @examples
#' # loadCutcount(gfoot_obj)
setGeneric(name="loadCutcount",
           def = function(obj) {
             standardGeneric("loadCutcount");
           });

#' @describeIn loadCutcount Method for GFoot objects.
setMethod(f="loadCutcount",
          signature="GFoot",
          definition=function(obj) {
            if (length(obj@control@cutcount)==0) {
              obj@control<-readCutCountSites(obj@control, obj@sitefile);
            }
            if (length(obj@treatment@cutcount)==0) {
              obj@treatment<-readCutCountSites(obj@treatment, obj@sitefile);
            }
            return(obj);
          });

#' Run GFoot Analysis (S4 Generic)
#'
#' Runs GFoot analysis for a GFoot object.
#'
#' @param obj S4 object of class GFoot.
#' @param range Optional. Range for analysis.
#' @param graphout Logical. If TRUE, output graphs.
#' @param yrange Numeric vector. Y-axis range for plots.
#' @param mc.cores Integer. Number of CPU cores to use.
#' @param run_set Optional. Additional run set parameter.
#'
#' @return S4 object of class GFoot with output files.
#' @export
#'
#' @examples
#' # run(gfoot_obj, range = 1:10, graphout = TRUE)
setGeneric(name="run",
           def = function(obj,range=NA,graphout=F,yrange=c(-2.2,0.5), mc.cores=MCCORES, run_set=NA) {
             standardGeneric("run");
           });

#' @describeIn run Method for GFoot objects.
setMethod(f="run",
          signature="GFoot",
          definition=function(obj, range=NA, graphout=F, yrange=c(-2.2,0.5), mc.cores= MCCORES, run_set=NA) {
            if (length(obj@control@cutcount)==0) {
              obj@control<-readCutCountSites(obj@control, obj@sitefile);
            }
            if (length(obj@treatment@cutcount)==0) {
              obj@treatment<-readCutCountSites(obj@treatment, obj@sitefile);
            }
            cat(sprintf("Num. of CPU cores to use: %d\n",mc.cores));
            MCCORES = mc.cores;
            obj@outputfiles=runGfoot(
              obj,
              range =  range,
              graphout = graphout,
              yrange= yrange,
              run_set= run_set
            );
            return(obj);
          });

#' Scatter Data Out Groups
#'
#' Generates scatter plot data for group comparisons and writes output files.
#'
#' @param outputdir Character. Directory for output files.
#' @param datanamelist List of data names.
#' @param motiflistlist List of motif lists.
#' @param sitefilelist List of site files.
#' @param halfwidthlist Numeric vector of half-width values.
#' @param threshold Numeric. Threshold for filtering.
#' @param comparisonIndices List of index pairs for comparison.
#' @param range Optional. Range for scatter plot.
#'
#' @return Character vector of output filenames.
#' @export
#'
#' @examples
#' # scatterDataOutGroups("outdir", list("A", "B"), list("motifs.txt"), list("sites.txt"), c(50,100), 0, list(c(1,2)), NA)
scatterDataOutGroups<-function(
    outputdir='',
    datanamelist = list(),
    motiflistlist = list(),
    sitefilelist = list(),
    halfwidthlist = c(50,100),threshold = 0,
    comparisonIndices = list(c(1,2), c(3,2)),
    range=NA) {

  out=c(); # Initialize an empty vector to store output filenames
  for (index in comparisonIndices) { # Loop through each pair of comparison indices
    for (halfwidth in halfwidthlist) { # Loop through each half-width value

      i1 = index[1]; i2 = index[2]; # Extract the indices for comparison

      dataname1 = datanamelist[[i1]]; # Get the first data name
      dataname2 = datanamelist[[i2]]; # Get the second data name

      sitefilename = extractBaseName(basename(sitefilelist[[i1]])); # Extract the base name of the site file

      # Determine the motif list file based on the length of motiflistlist
      if (length(motiflistlist)==1) {
        motiflistfile1 =  motiflistlist
      } else {
        motiflistfile1 = motiflistlist[[i1]]
      }

      # Get the output directory and filename for the first and second data names
      out1=getOutputDirAndFilename(outputdir, datanamelist[[i1]], sitefilelist[[i1]], halfwidth);
      out2=getOutputDirAndFilename(outputdir, datanamelist[[i2]], sitefilelist[[i2]], halfwidth);
      outputfile1=out1$filename;
      outputfile2=out2$filename;

      # Create the title for the plot
      title = sprintf('Footprinting Aggregation Plot (%s vs. %s)', dataname1, dataname2);

      # Create the data directories for the first and second data names
      datadir1 =  sprintf("%s_%dbp",dataname1, halfwidth);
      datadir2 =  sprintf("%s_%dbp",dataname2, halfwidth);
      MAoutputdir = file.path(outputdir, 'comparison');

      # Create the output CSV file name
      outputcsvfile = sprintf("%s_%s_On_%s_footprint_depth_table.csv", dataname1, dataname2, sitefilename);

      # Call the ScatterDataOut function and store the output filename
      outfilename=ScatterDataOut(dataname1=dataname1,
                                 dataname2=dataname2,
                                 motiflistfile=motiflistfile1,
                                 outputfile1=outputfile1,
                                 outputfile2=outputfile2,
                                 title = title,threshold=threshold,datadir1=datadir1,datadir2=datadir2, outputcsv=
                                   outputcsvfile,
                                 range=range,
                                 #added to be able to write to different location
                                 outputdir=outputdir);
      out=c(out, outfilename); # Append the output filename to the list
    }
  }
  out; # Return the list of output filenames
}

#' Scatter Data Out
#'
#' Processes motif and cut count data, calculates statistics, and writes results to a CSV file.
#'
#' @param dataname1 Character. Name of first dataset.
#' @param dataname2 Character. Name of second dataset.
#' @param motiflistfile Character. Path to motif list file.
#' @param outputfile1 Character. Path to first output file.
#' @param outputfile2 Character. Path to second output file.
#' @param title Character. Plot title.
#' @param threshold Numeric. Threshold for filtering.
#' @param datadir1 Character. Directory for first dataset.
#' @param datadir2 Character. Directory for second dataset.
#' @param outputcsv Character. Output CSV file name.
#' @param range Optional. Range for motif selection.
#' @param outputdir Character. Output directory.
#'
#' @return Character. Path to output CSV file.
#' @export
#'
#' @examples
#' # ScatterDataOut("A", "B", "motifs.txt", "out1.txt", "out2.txt", "Title", 0, "dir1", "dir2", "output.csv")
ScatterDataOut <- function(dataname1="", dataname2="", motiflistfile="",
                           outputfile1="",outputfile2="",title = "",threshold=0,datadir1='',datadir2='', outputcsv=' ',range=NA, outputdir='') {

  # Read the motif list from the file
  motiflist=read.table(motiflistfile);

  # Function to calculate observed means
  calcObsMeans <- function (mobs1, mobs2,motifcenter, motiflogo, meanlog1, meanlog2, numsite) {
    motifoffset=2;
    halfwidth= (length(mobs1)-1)/2;
    motifwidth = nchar(motiflogo);
    xr =  seq(from=1, to= nchar(motiflogo)) - motifcenter;
    hgt = 0;
    x1 = min(xr); x2 = max(xr);
    linea = x1 - motifoffset;
    lineb = x2 + motifoffset;
    motifregion = linea:lineb  + halfwidth + 1;
    flanking = setdiff(1:length(mobs1), motifregion);
    c1=log(sum(mobs2[motifregion])/sum(mobs1[motifregion]),2);
    c2=log(sum(mobs2[flanking])/sum(mobs1[flanking]),2);
    c3=log(sum(mobs2)/sum(mobs1),2);
    meancc1 = mean(mobs1);
    meancc2 = mean(mobs2);
    data.frame(lrmotif=c1, lrflanking=c2, lrtotal=c3, fdepth1= meanlog1, fdepth2= meanlog2,meancc1= meancc1, meancc2=meancc2, numsite);
  }

  # Function to process each motif
  calcprog<-function(ii) {
    motif = motiflist[ii,];
    name = as.character(motif$motif);
    print(name);

    sites1 = output1[ii,];
    sites2 = output2[ii,];

    numsite1 = sites1$numSites;
    numsite2 = sites2$numSites;

    motif = motiflist[ii,];
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

    plot1yrange=0;

    halfwidth= (length(logratio1)-1)/2;
    filenamebase = make.names(name);

    calcObsMeans(mobs1, mobs2,motifcenter,motiflogo, meanlog1, meanlog2, numsite1);
  }

  # Filter motif list if range is specified
  if (!is.na(range[1])) {
    motiflist <- motiflist[range,];
  }

  # Read output files
  output1 = read.table(outputfile1,sep='\t');  #200
  output2 = read.table(outputfile2,sep='\t');  #199

  # Select motifs based on threshold
  selection= (output1$numSites>=threshold) & (output2$numSites>=threshold);

  output1 = output1[selection,];
  output2 = output2[selection,];
  motiflist = motiflist[selection,];

  # Calculate and combine results
  cal<-do.call("rbind", lapply(1:nrow(motiflist), calcprog));
  tabl = data.frame(motiflist, cal);

  # Write results to output file
  write.table(tabl, file=paste0(outputdir,"/",outputcsv))
  outputcsv;
}



#' Count Cut Count on Sites
#'
#' Counts cut events within specified genomic sites.
#'
#' @param cutcount List. Cut count data per chromosome.
#' @param sites Data frame or character. Genomic sites or path to sites file.
#'
#' @return Numeric. Total cut count across all sites.
#' @export
#'
#' @examples
#' # countCutcountOnSites(cutcount, sites)
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
  #print(paste("chroms:",chroms))
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
  sum_for_chr_es<-sum(unlist(parallel::mclapply(chroms, sum_for_chr,mc.cores=MCCORES))) #modif LC
  # sum_for_chr_es<-sum(unlist(lapply(chroms, sum_for_chr))) #modif LC
  print(paste("sum_for_chr_es:",sum_for_chr_es))
  sum_for_chr_es;
}

#' Read Cut Count on Sites
#'
#' Reads and reduces cut count data based on overlapping genomic sites.
#'
#' @param cutcount_file Character. Path to cut count file.
#' @param site_file Character. Path to site file.
#'
#' @return List. Reduced cut count data per chromosome.
#' @export
#'
#' @examples
#' # readCutCountOnSites("cutcount.bgr", "sites.txt")
readCutCountOnSites<-function(cutcount_file="", site_file="") {

  cutcount = NULL;
  sites = NULL;

  reduced_bgr_filename = file.path(dirname(cutcount_file),sprintf("%s.bgr", paste(extractBaseName(cutcount_file),'on',extractBaseName(site_file), sep="_")));

  reduced_bgr_filename_zip = paste(reduced_bgr_filename,".gz", sep="");

  overlapping<-function(A,B)
  {
    if (as.character(A$chr[1]) != as.character(B$chr[1])) {
      null;
    } else {
      R= vector("logical", length=nrow(A));

      rangeA = 1:length(A$st);
      rangeB = 1:length(B$st);

      if ((length(rangeA) !=0) && (length(rangeB) !=0))
      {
        P=list(Ast=A$st[rangeA],Aed=A$ed[rangeA],Bst= B$st[rangeB],Bed= B$ed[rangeB]);
        RR=  intersect_intervals_really_fast_parallel(P);
        I1 = unique(RR);
        R[rangeA[I1]]=TRUE;
      }
      # print(c(i,cend, cendcomp));
      R;
    }
  }


  filterCutcountByHotspot<-function(sites, chr) {
    cat(sprintf("reducing cutcount data for %s\n", chr));
    hotspot = sites[sites$chr==chr,c("chr","st","ed")];
    extendHotspot<-function(hotspot, extensionbp = 200)  {
      hotspot$st = hotspot$st - extensionbp;
      hotspot$ed = hotspot$ed + extensionbp;
      hotspot = hotspot[hotspot$st > 0,];
      hotspot;
    }
    if (nrow(hotspot)==0) {
      subcutcount=NULL;
    } else {
      #	browser();
      exh = extendHotspot(hotspot, extensionbp = 200);
      cutcountchr = cutcount[[chr]];

      ovlp = overlapping(cutcountchr, exh);
      subcutcount = cutcountchr[ovlp,];
    }
    subcutcount;
  };

  if (file.exists(reduced_bgr_filename_zip)) {
    cat(sprintf("reading the cutcount data file: %s\n", reduced_bgr_filename_zip));
    cutcount = readcutcount2(reduced_bgr_filename_zip);
  } else {
    cat(sprintf("converting the cutcount data file: %s into %s\n", cutcount_file, reduced_bgr_filename));
    cutcount = readcutcount2(cutcount_file);
    sites = readhotspot(site_file);

    #	browser();
    reduced_cutcount= parallel::mclapply(names(cutcount), function(x) { filterCutcountByHotspot(sites, x); }, mc.cores=MCCORES);
    names(reduced_cutcount)<- names(cutcount);

    for (ii in seq_along(reduced_cutcount)) {
      chr = names(reduced_cutcount)[[ii]];
      dat = reduced_cutcount[[ii]];
      if (!is.null(dat)) {
        cat(sprintf("saving cutcount data (%s) for %s\n",reduced_bgr_filename, chr));
        dat$st = dat$st - 1;     # zero-based format
        if (ii==1) {
          write.table(dat, file=reduced_bgr_filename, sep=" ",   append=F, quote=F, row.names=F, col.names=F);
        } else {
          write.table(dat, file=reduced_bgr_filename, sep=" ",   append=T, quote=F, row.names=F, col.names=F);
        }
      }
    }
    #browser();
    system(sprintf('gzip %s', reduced_bgr_filename));  #gzip
    cutcount = reduced_cutcount;
  }
  cutcount;
}

