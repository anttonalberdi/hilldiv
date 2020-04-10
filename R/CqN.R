#' Sørensen-type overlap
#' @title Sørensen-type overlap
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords beta dissimilarity similarity
#' @description The Sørensen-type overlap quantifies the effective average proportion of a sub-systems OTUs (or lineages in the case of phylodiversities) that is shared across all subsystems. This is thus a metric that quantifies overlap from the subsystems perspective. Its corresponding dissimilarity measure (1 - CqN) quantifies the effective average proportion of nonshared OTUs or lineages in a system. CqN is integrated in the functions beta_dis() and pair_dis().
#' @param beta A beta diversity value based on Hill numbers.
#' @param qvalue The q value used to compute the beta diversity. It needs to be a positive number, usually between 0 and 5, but most commonly 0, 1 or 2. It can be an integer or contain decimals.
#' @param N An integer indicating sample size, the number of sampling units to be used to compute the similarity measure.
#' @return A Sørensen-type overlap value
#' @seealso \code{\link{div_part}}, \code{\link{beta_dis}}
#' @examples
#' CqN(beta=1.24,qvalue=1,N=3)
#' CqN(1.24,1,3)
#' @references
#' Alberdi, A., Gilbert, M.T.P. (2019). A guide to the application of Hill numbers to DNA-based diversity analyses. Molecular Ecology Resources, 19, 804-817.\cr\cr
#' Chao, A., Chiu, C.H., & Hsieh, T. C. (2012). Proposing a resolution to debates on diversity partitioning. Ecology, 93, 2037-2051.\cr\cr
#' Jost, L. (2007). Partitioning diversity into independent alpha and beta components. Ecology, 88, 2427-2439.
#' @export

CqN <- function(beta,qvalue,N){
if(qvalue==1){qvalue=0.99999}
value = ((1/beta)^(qvalue-1) - (1/N)^(qvalue-1)) / (1 - (1/N)^(qvalue-1))
return(value)
}
