#' Sørensen-type turnover-complement
#' @title Sørensen-type turnover-complement
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords beta dissimilarity similarity
#' @description The Sørensen-type turnover-complement is the complement of the Sørensen-type turnover, which quantifies the normalized OTU turnover rate with respect to the average subsystem (i.e., alpha), thus provides the proportion of a typical subsystem that changes across subsystems. VqN is integrated in the functions beta_dis() and pair_dis().
#' @param beta A beta diversity value based on Hill numbers.
#' @param N An integer indicating sample size, the number of sampling units to be used to compute the similarity measure.
#' @return A Sørensen-type turnover-complement value
#' @seealso \code{\link{div_part}}, \code{\link{beta_dis}}
#' @examples
#' VqN(beta=1.24,N=2)
#' VqN(1.24,2)
#' @references
#' Alberdi, A., Gilbert, M.T.P. (2019). A guide to the application of Hill numbers to DNA-based diversity analyses. Molecular Ecology Resources, 19, 804-817.\cr\cr
#' Chao, A., Chiu, C.H., & Hsieh, T. C. (2012). Proposing a resolution to debates on diversity partitioning. Ecology, 93, 2037-2051.\cr\cr
#' Jost, L. (2007). Partitioning diversity into independent alpha and beta components. Ecology, 88, 2427-2439.
#' @export

VqN <- function(beta,N){
value = (N - beta)/(N-1)
return(value)
}
