#' To occurrences
#' @title Transform abundance vector or matrix into occurrences
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords OTU ASV table abundance hill occurrence
#' @description Transform an absolute or relative abundance vector (one sample) or matrix (multipla samples) into an occurrence vector or matrix.
#' @param abund A vector or a matrix/data.frame indicating the absolute or relative abundances of multiple samples. Columns must refer to samples and rows to OTUs/ASVs.
#' @return A vector (one sample) or matrix (multiple samples) of occurrence data.
#' @seealso \code{\link{to.incidence}}
#' @examples
#' data(bat.diet.otutable)
#' to.occurrences(bat.diet.otutable)
#' to.occurrences(bat.diet.otutable[,1])
#' @references
#' Alberdi, A., Gilbert, M.T.P. (2019). A guide to the application of Hill numbers to DNA-based diversity analyses. Molecular Ecology Resources, 19, 804-817.
#' @export

to.occurrences <- function(abund){
if(missing(abund)) stop("Abundance vector or table is missing")
abund[abund != 0] <- 1
return(abund)
}
