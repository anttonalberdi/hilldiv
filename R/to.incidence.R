#' To incidence
#' @title Transform to incidence
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords OTU ASV table abundance hill
#' @description Transform a count (OTU/ASV) table from abundance to incidence.
#' @param otutable A matrix/data.frame indicating the (relative) abundances of multiple samples. Columns must refer to samples and rows to OTUs/ASVs.
#' @param hierarchy A two-column matrix indicating the relation between samples (first column) and groups (second column).
#' @param relative Whether to transform the incidence vector or matrix to relative (0-1) values. Default: relative=FALSE.
#' @return A vector of incidence data of a single system if no hierarchy table is specified and a matrix of incidence data of multiple systems if a hierarchy table is specified.
#' @seealso \code{\link{hill_div}}, \code{\link{div_part}}
#' @examples
#' data(bat.diet.otutable)
#' data(bat.diet.hierarchy)
#' to.incidence(bat.diet.otutable)
#' to.incidence(bat.diet.otutable,bat.diet.hierarchy)
#' to.incidence(bat.diet.otutable,bat.diet.hierarchy,relative=TRUE)
#' to.incidence(otutable=bat.diet.otutable,hierarchy=bat.diet.hierarchy,relative=TRUE)
#' @references
#' Alberdi, A., Gilbert, M.T.P. (2019). A guide to the application of Hill numbers to DNA-based diversity analyses. Molecular Ecology Resources, 19, 804-817.
#' @export

to.incidence <- function(otutable,hierarchy,relative){

if(missing(otutable)) stop("Count table is missing")
if(is.null(dim(otutable)) == TRUE) stop("The count table needs to be a bi-dimensional matrix")
if(missing(relative)){relative=FALSE}

if(missing(hierarchy)){
  #One system
  inci.vector <- rowSums(otutable != 0)
  if(relative == FALSE){
    return(inci.vector)
    }
  if(relative == TRUE){
    inci.vector.rel <- tss(inci.vector)
    return(inci.vector.rel)
    }
}else{
  #Multiple systems
  if(ncol(hierarchy) != 2) stop("The hierarchy table needs to have two columns")
  if(length(unique(hierarchy[,1])) == length(unique(hierarchy[,2]))) stop("The number of groups needs to be smaller than the number of samples")
  systems <- unique(hierarchy[,2])
  inci.otutable <- c()
  for(s in systems){
    samples <- as.character(hierarchy[which(hierarchy[,2] == s),1])
    otutable.subset <- otutable[,samples]
    inci.vector <- rowSums(otutable.subset != 0)
    inci.otutable <- cbind(inci.otutable,inci.vector)
  }
  colnames(inci.otutable) <- systems
  if(relative == FALSE){
    return(inci.otutable)
    }
  if(relative == TRUE){
    inci.otutable.rel <- tss(inci.otutable)
    return(inci.otutable.rel)
  }

}
}
