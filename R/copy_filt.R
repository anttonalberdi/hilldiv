#' OTU/ASV copy number filtering
#' @title OTU/ASV copy number filtering
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords OTU ASV threshold
#' @description As DNA sequencing data include PCR and sequencing errors, copy number thresholds are commonly applied to discard the OTUs with low number of sequence copies. This threshold can be absolute or (ideally) relative to the sequencing depth of each sample.
#' @param abund A vector or a matrix/data.frame indicating the relative abundances of one or multiple samples, respectively. If a matrix/data.frame is provided, columns must refer to samples and rows to OTUs.
#' @param threshold An integer or a decimal number indicating the cut-off threshold. If an integer is provided, an absolute threshold is used (same threshold for all samples). If a decimal number is provided a relative copy number threshold is applied (dependent on the sequencing depth of each sample).
#' @param filter Whether to remove the OTUs/ASVs with no read counts. Default=TRUE.
#' @seealso \code{\link{depth_cov}}, \code{\link{tss}}
#' @examples
#' data(bat.diet.otutable)
#' #Remove singletons from all samples
#' copy_filt(bat.diet.otutable,2)
#' #Remove OTUs represented by less than 0.01% of the total reads per sample.
#' copy_filt(bat.diet.otutable,0.0001)
#' @references
#' Alberdi A, Aizpurua O, Bohmann K, Gopalakrishnan S, Lynggaard C, Nielsen M, Gilbert MTP. 2019. Promises and pitfalls of using high-throughput sequencing for diet analysis. Molecular Ecology Resources, 19(2), 327-348.
#' @export

copy_filt <- function(abund,threshold,filter){
  if(missing(filter)){filter=TRUE}
  if(is.null(dim(abund)) == FALSE){
  #It is an OTU table (multiple samples)
  sapply(1:ncol(abund), function(colnum){temp = abund[,colnum]
    if(threshold == round(threshold)){
    #Absolute
    rownums = which(temp < threshold)
    }else{
    #Relative
    rownums = which(temp < sum(temp)*threshold)
    }
  abund[rownums, colnum] <<- 0})
  if(filter == TRUE){
  abund <- abund[!apply(abund, 1, function(x) all(x == 0)), ]
  }
  }else{
  #It is a vector (single sample)
    if(threshold == round(threshold)){
    #Absolute
    abund[abund < threshold] <- 0
    }else{
    #Relative
    abund[abund < sum(abund)*threshold] <- 0
    }
  if(filter == TRUE){
  abund <- abund[abund != 0]
  }
  }
  return(abund)
}
