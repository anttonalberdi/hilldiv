#' Tree depth
#' @title Tree depth
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords tree phylogeny
#' @description Computes phylogenetic tree depth based from a phylogenetic tree and a vector of (relative) abundances.
#' @param tree A phylogenetic tree of class 'phylo'. The tip labels must match the row names in the OTU table. Use the function match_data() if the OTU names do not match.
#' @param abund A vector or a matrix/data.frame indicating the relative abundances of one or multiple samples, respectively. If a matrix/data.frame is provided, columns must refer to samples and rows to OTUs.
#' @return A tree depth value
#' @seealso \code{\link{div_part}}, \code{\link{gamma_div}}, \code{\link{match_data}}
#' @importFrom geiger tips
#' @examples
#' data(bat.diet.otutable)
#' data(bat.diet.tree)
#' tree_depth(tree=bat.diet.tree,abund=bat.diet.otutable)
#' tree_depth(bat.diet.tree,bat.diet.otutable)
#' @export

tree_depth <- function(tree,abund){
if(class(tree) != "phylo") stop("Tree needs to be an object of class Phylo")
if(is.null(dim(abund)) == TRUE){
  vector <- tss(abund)
  }
if(is.null(dim(abund)) == FALSE){
  vector <- tss(abund[,1])
  names(vector) <- rownames(abund)
  }
if(identical(sort(names(vector)),sort(tree$tip.label)) == FALSE) stop("OTU names in the vector and tree do not match")
vector <- tss(vector)
Li <- tree$edge.length
ltips <- sapply(tree$edge[, 2], function(node) tips(tree, node))
ai <- unlist(lapply(ltips, function(TipVector) sum(vector[TipVector])))
T <- sum(Li * ai)
return(T)
}
