#' CutCount Class
#'
#' S4 class to store cut count data and related metadata.
#'
#' @slot file Character. Path to the cut count file.
#' @slot cutcount List. Cut count data per chromosome.
#' @slot count Numeric. Total cut count.
#' @slot name Character. Name or label for the dataset.
#'
#' @return A CutCount object.
#' @export
#'
#' @examples
#' # cc <- CutCount(file = "sample.bgr", count = 1000, name = "Sample1")
CutCount <- setClass("CutCount",
                     representation(file = "character",
                                    cutcount = "list",
                                    count = "numeric",
                                    name = "character"),
                     prototype = list(file="",
                                      cutcount = list(),
                                      count = 0,
                                      name=""),
                     validity = function(object) {
                       if (file.exists(object@file) && object@count > 0) TRUE else "File not exist or count==0";
                     });



#' MotifDB Class
#'
#' S4 class to store motif database information.
#'
#' @slot motiflistfile Character. Path to motif list file.
#' @slot directory Character. Directory containing motif files.
#'
#' @return A MotifDB object.
#' @export
#'
#' @examples
#' # mdb <- MotifDB(motiflistfile = "motifs.txt", directory = "motif_dir")
MotifDB <- setClass("MotifDB",
                    representation(
                      motiflistfile = "character",
                      directory = "character"
                    ),
                    prototype = list(
                      motiflistfile ='',
                      directory='')
);


#' GFootOption Class
#'
#' S4 class to store GFoot analysis options.
#'
#' @slot biasfile Character. Path to bias correction file.
#' @slot hexamer_pattern Character. Pattern for hexamer code files.
#'
#' @return A GFootOption object.
#' @export
#'
GFootOption <- setClass("GFootOption",
                        representation(
                          biasfile = "character",
                          hexamer_pattern = "character"),
                        prototype = list(
                          biasfile="",
                          hexamer_pattern='')
);


#' GFoot Class
#'
#' S4 class to store all data and options for a GFoot analysis.
#'
#' @slot control CutCount. Control cut count data.
#' @slot treatment CutCount. Treatment cut count data.
#' @slot sitefile Character. Path to site file.
#' @slot motifDB MotifDB. Motif database object.
#' @slot gfootoption GFootOption. GFoot analysis options.
#' @slot outputdir Character. Output directory.
#' @slot cachedir Character. Cache directory.
#' @slot outputfiles Character. Output files generated.
#'
#' @return A GFoot object.
#' @export
#'
#' @examples
#' #gfoot <- GFoot(control = cc1, treatment = cc2, sitefile = "sites.txt", motifDB = mdb, gfootoption = opt)

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



#' Run GFoot Analysis (S4 Generic)
#'
#' Generic method to run GFoot analysis.
#'
#' @param obj GFoot object.
#' @param range Optional. Range of motifs to analyze.
#' @param graphout Logical. If TRUE, output graphs.
#' @param yrange Numeric vector. Y-axis range for plots.
#' @param mc.cores Integer. Number of CPU cores to use.
#' @param run_set Character. Optional run set identifier.
#'
#' @return GFoot object with output files.
#' @export
setGeneric(name="run",    #GFoot::run
           def = function(obj,range=NA,graphout=F,yrange=c(-2.2,0.5), mc.cores=MCCORES, run_set=NA) {
             standardGeneric("run");
           });

#' @describeIn run Method for GFoot objects.
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
          });
