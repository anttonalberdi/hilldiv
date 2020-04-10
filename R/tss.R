#' Total Sum Scaling normalisation
#' @title Total Sum Scaling normalisation
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords normalisation Hill
#' @description Normalise a vector or count matrix to the range of 0-1.
#' @param abund A vector or a matrix/data.frame indicating the relative abundances of one or multiple samples, respectively. If a matrix/data.frame is provided, columns must refer to samples and rows to OTUs.
#' @return Normalised vector or matrix.
#' @seealso \code{\link{hill_div}}, \code{\link{index_div}}
#' @examples
#' data(bat.diet.otutable)
#' tss(bat.diet.otutable)
#' bat.diet.sample <- bat.diet.otutable[,1]
#' tss(bat.diet.sample)
#' @export

tss <- function(abund){
  #If input data is a vector
  if(is.null(dim(abund)) == TRUE){
  abund.norm <- abund/sum(abund)}
  #If input data is an OTU table
  if(is.null(dim(abund)) == FALSE){
  abund.norm <- sweep(abund, 2, colSums(abund), FUN="/")}
  return(abund.norm)
}
