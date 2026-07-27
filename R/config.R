# configuration
#' All the classes that we need to use in the package
#'
#' @slot file character.
#' @slot cutcount list.
#' @slot count numeric.
#' @slot name character.
#'
#' @returns nothing
#' @export
#'
#' @examples

# Configuration

#MCCORES=2;	         # Number of CPU cores for parallel computation
MCCORES=max(2, parallel::detectCores() - 1) #modif LC

#MCCORES_MOTIF = 1;     # Number of motifs calculated at the same time
MCCORES_MOTIF = 4 #modif LC

FAST_AGGREGATION_MC_CORES = 4; # Number of CPU cores for less memory-dependent computation

MAX_MOTIF_SITES = 10000000;    # Maximum number of FIMO sites for each motif

BFOOT_DEBUG = T              # Display the debug messages
