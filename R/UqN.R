#' Jaccard-type overlap
#' @title Jaccard-type overlap
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords beta dissimilarity similarity
#' @description The Jaccard-type overlap quantifies the effective proportion of OTUs or lineages in a system that are shared across all subsystems. Hence, this metric quantifies overlap from the perspective of the overall system. Its corresponding dissimilarity (1 - UqN) quantifies the effective proportion of nonshared OTUs or lineages in the overall system. UqN is integrated in the functions beta_dis() and pair_dis().
#' @param beta A beta diversity value based on Hill numbers.
#' @param qvalue The q value used to compute the beta diversity. It needs to be a positive number, usually between 0 and 5, but most commonly 0, 1 or 2. It can be an integer or contain decimals.
#' @param N An integer indicating sample size, the number of sampling units to be used to compute the similarity measure.
#' @return A Jaccard-type overlap value
#' @seealso \code{\link{div_part}}, \code{\link{beta_dis}}
#' @examples
#' UqN(beta=1.24,qvalue=1,N=2)
#' UqN(1.24,1,2)
#' @references
#' Alberdi, A., Gilbert, M.T.P. (2019). A guide to the application of Hill numbers to DNA-based diversity analyses. Molecular Ecology Resources, 19, 804-817.\cr\cr
#' Chao, A., Chiu, C.H., & Hsieh, T. C. (2012). Proposing a resolution to debates on diversity partitioning. Ecology, 93, 2037-2051.\cr\cr
#' Jost, L. (2007). Partitioning diversity into independent alpha and beta components. Ecology, 88, 2427-2439.
#' @export

UqN <- function(beta,qvalue,N){
if(qvalue==1){qvalue=0.99999}
value = ((1/beta)^(1-qvalue) - (1/N)^(1-qvalue)) / (1 - (1/N)^(1-qvalue))
return(value)
}
