#' Alpha diversity computation (based on Hill numbers)
#' @title Alpha diversity
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords alpha partitioning hill
#' @description Compute alpha diversity of a system comprised of multiple samples from a count (OTU/ASV) table. If a tree object is provided, the computed alpha diversity accounts for the phylogenetic relations across OTUs/ASVs.
#' @param countable A count table (matrix/data.frame) indicating the absolute or relative OTU/ASV abundances of multiple samples. Columns must refer to samples and rows to OTUs/ASVs.
#' @param qvalue A positive number, usually between 0 and 5, but most commonly 0, 1 or 2. It can be an integer or contain decimals.
#' @param tree A phylogenetic tree of class 'phylo'. The tip labels must match the row names in the count table. Use the function match_data() if the count table and tree names do not match.
#' @param weight A vector indicating the relative weight of each sample. The order needs to be identical to the order of the samples in the OTU table. The values need to sum up to 1. If empty, all samples are weighed the same.
#' @usage alpha_div(countable,qvalue,tree,weight)
#' @return An alpha diversity value.
#' @seealso \code{\link{div_part}}, \code{\link{gamma_div}}, \code{\link{match_data}}
#' @importFrom ape is.ultrametric
#' @importFrom geiger tips
#' @examples
#' data(bat.diet.otutable)
#' data(bat.diet.tree)
#' alpha_div(countable=bat.diet.otutable,qvalue=1)
#' alpha_div(countable=bat.diet.otutable,qvalue=1,tree=bat.diet.tree)
#' weight.vector = rep(1/ncol(bat.diet.otutable),ncol(bat.diet.otutable))
#' alpha_div(bat.diet.otutable,1,bat.diet.tree,weight.vector)
#' @references
#' Alberdi, A., Gilbert, M.T.P. (2019). A guide to the application of Hill numbers to DNA-based diversity analyses. Molecular Ecology Resources, 19, 804-817.\cr\cr
#' Chao, A., Chiu, C.H., & Hsieh, T. C. (2012). Proposing a resolution to debates on diversity partitioning. Ecology, 93, 2037-2051.\cr\cr
#' Jost, L. (2007). Partitioning diversity into independent alpha and beta components. Ecology, 88, 2427-2439.
#' @export

alpha_div <- function(countable,qvalue,tree,weight){

#Quality-check and warnings
if(missing(countable)) stop("OTU table is missing")
if(is.null(dim(countable)) == TRUE) stop("The OTU table is not a matrix")
if(dim(countable)[1] < 2) stop("The OTU table only less than 2 OTUs")
if(dim(countable)[2] < 2) stop("The OTU table contains less than 2 samples")
if(sum(colSums(countable)) != ncol(countable)) {countable <- tss(countable)}
if(missing(qvalue)) stop("q value is missing")
if(qvalue < 0) stop("q value needs to be possitive (equal or higher than zero)")
if (qvalue==1) {qvalue=0.99999}
if(missing(weight)) { weight= rep(1/ncol(countable),ncol(countable))}

if(missing(tree)){
        #Neutral
        pi <- as.data.frame(countable[apply(countable, 1, function(z) !all(z==0)),])
        pi.w <- sweep(pi,2,weight,"*")
        pi.w.q <- pi.w^qvalue
        pi.w.q[!pi] <- 0
        N <- length(weight)
        div <- sum(rowSums(pi.w.q))^(1/(1-qvalue))/N
        return(div)
}else{
        #Phylogenetic
        if(class(tree) != "phylo") stop("Tree needs to be an object of class Phylo")
        if(is.ultrametric(tree) == FALSE) stop("Tree needs to be ultrametric")
        if(identical(sort(rownames(countable)),sort(tree$tip.label)) == FALSE) stop("OTU/ASV names in the count table and tree do not match. Use match_data() to solve this issue.")
        countable <- as.data.frame(countable)
        wj <- weight
        N <- ncol(countable)
        Li <- tree$edge.length
        ltips <- sapply(tree$edge[, 2], function(node) geiger::tips(tree, node))
        aij <- matrix(unlist(lapply(ltips, function(TipVector) colSums(countable[TipVector,]))), ncol = N, byrow = TRUE)
        aij.wj <- sweep(aij, 2, wj, "*")
        T <- sum(sweep(aij.wj, 1, Li, "*"))
        L <- matrix(rep(Li, N), ncol = N)
        wm <-  matrix(rep(wj, length(Li)), ncol = N,byrow=TRUE)
        i <-  which(aij > 0)
        phylodiv <- sum(L[i] * (aij[i]*wm[i]/T)^qvalue)^(1/(1 - qvalue))/(N*T)
        return(phylodiv)
}
}
