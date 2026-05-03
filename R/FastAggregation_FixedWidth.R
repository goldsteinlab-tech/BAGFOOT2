
#' @useDynLib tiger CalcAggregationPerChromFixedWidthInC
CalcAggregationPerChromFastFixedWidthInC<- function(P) {
  oldw <- getOption("warn");
  options(warn = -1);

  region =  P$region;     # Array of regions
  density= P$density;      # Array of Tag Density

  x1 = as.integer(density$st);
  x2 = as.integer(density$ed);
  y1 = as.integer(region$st);
  y2 = as.integer(region$ed);

  if (all(is.na(P$dir))) {
    dir = seq(from=1,to=1, length=length(y1));
  } else {
    dir = P$dir;
  }

  width=y2[1]-y1[1]+1;
  Result=.Call("CalcAggregationPerChromFixedWidthInC", x1,x2,y1,y2,as.integer(dir), density$value);
  options(warn = oldw);
  array(data=Result, dim=c(length(y1),width));
}

#' @export
CalcHeatmapFastFixedWidthC<- function(tregion, tdensity, dir=NULL) {
  #browser();
  chrom = levels(factor(tregion$chr ));
  locchrom = as.character(tregion$chr);

  if (is.null(dir)) {
    dir = rep(NA, nrow(tregion));
  }


  numchrom = length(chrom);
  Plist = {};
  nregion = nrow(tregion);
  range=seq(from=1, to=nregion);

  runCalcAggregationPerChromFastPerChrom<- function(ch) {
    P = {};
    #print(ch);

    P$density = tdensity[[ch]];
    maxcoord = max(P$density$ed);

    P$region = subset(tregion, chr==ch ); # & ed < maxcoord
    #	browser();
    P$dir = subset(dir, tregion$chr==ch);
    withinRange = P$region$ed < maxcoord;
    outRange = !withinRange;
    P$region = P$region[withinRange,];
    P$dir = P$dir[withinRange];
    R = CalcAggregationPerChromFastFixedWidthInC(P);
    RC =array(0,dim=c(length(withinRange), ncol(R)));
    RC[withinRange,] = R;
    RC[outRange,] = array(0, dim=c(sum(outRange),ncol(R)));
    #browser(text="line#1623");
    RC;
  }

  chrs=unique(as.character(tregion$chr));

  result = lapply(chrs, runCalcAggregationPerChromFastPerChrom);
  do.call("rbind", result);
}

#' @export
CalcHeatmapFastMatrixOutputC <- function(tregion, tdensity,dir=NULL) {
  out = CalcHeatmapFastFixedWidthC(tregion, tdensity, dir=dir);
  out;
}
