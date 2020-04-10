#' Diversity index computation
#' @title Diversity index computation
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords OTU ASV diversity Hill
#' @description Computes common diversity indices related to Hill numbers. If the input is a vector, the function computes the indices of a single sample, while if the input is a matrix (OTU table), the function computes individual diversity indices for each sample (column). An ultrametic OTU tree is required for computing phylogenetic diversity indices (Faith's PD, Allen's H and Rao's Q). If the relative abundances of each sample (vector or each column of the matrix) do not sum to 1, TSS normalisation is applied.
#' @param abund A vector or a matrix/data.frame indicating the relative abundances of one or multiple samples, respectively. If a matrix/data.frame is provided, columns must refer to samples and rows to OTUs.
#' @param tree An ultrametic tree of class 'phylo'. The tip labels must match the names of the vector values (if one sample) or matrix rows (if multiple samples). Use the function match_data() if the OTU names do not match.
#' @param index Diversity index to be computed ("richness", "shannon", "simpson", "faith", "allen", "rao"). Default without tree argument: index="richness". Default with tree argument: index="faith".
#' @seealso \code{\link{hill_div}}, \code{\link{div_part}}
#' @importFrom ape cophenetic.phylo
#' @examples
#' data(bat.diet.otutable)
#' data(bat.diet.tree)
#' data(bat.diet.hierarchy)
#' #One sample
#' bat.diet.sample <- bat.diet.otutable[,1]
#' index_div(bat.diet.sample)
#' index_div(bat.diet.sample,index="shannon")
#' #Multiple samples
#' index_div(bat.diet.otutable)
#' index_div(bat.diet.otutable,tree=bat.diet.tree,index="faith")
#' #Incidence-based
#' bat.diet.otutable.incidence <- to.incidence(bat.diet.otutable,bat.diet.hierarchy)
#' index_div(bat.diet.otutable.incidence)
#' index_div(bat.diet.otutable.incidence,index="simpson")
#' index_div(to.incidence(bat.diet.otutable,bat.diet.hierarchy),tree=bat.diet.tree)
#' @references
#' Alberdi, A., Gilbert, M.T.P. (2019). A guide to the application of Hill numbers to DNA-based diversity analyses. Molecular Ecology Resources, 19, 804-817.\cr\cr
#' Jost, L. (2006). Entropy and diversity. Oikos, 113, 363-375.\cr\cr
#' Rao, C. R. (1982). Diversity and dissimilarity coefficients: A unified approach. Theoretical Population Biology, 21, 24-43.\cr\cr
#' Shannon, C. E. (1948). A mathematical theory of communication. The Bell System Technical Journal, 27, 379-423.\cr\cr
#' @export

index_div <- function(abund,tree,index){

#index_div(bat.diet.otutable)
#index_div(bat.diet.otutable,index="shannon")
#index_div(bat.diet.otutable,tree=bat.diet.tree,index="faith")

#Data input control
if(missing(abund)) stop("Abundance data is missing")
if(missing(index) & missing(tree)){index="richness"}
if(missing(index) & !missing(tree)){index="faith"}

#Input data type identification
if(is.null(dim(abund)) == TRUE){
inputtype="onesample"
}else{
inputtype="multiplesamples"
}

#INDICES

#Richness
if(index == "richness"){
  if(!missing(tree)) warning("Phylogenetic tree was not used to compute richness.")
  #Function
  richness <- function(vector){sum(vector != 0)}

  #One sample
  if(inputtype == "onesample"){
  div <- richness(abund)
  return(div)
  }

  #Multiple samples
  if(inputtype == "multiplesamples"){
  divs <- apply(abund, 2, function(x) richness(x))
  return(divs)
  }
}

#Shannon index
if(index == "shannon"){
  if(!missing(tree)) warning("Phylogenetic tree was not used to compute Shannon index.")
  #Function
  shannon <- function(vector){
    pi <- tss(vector[vector != 0])
    div <- -sum(pi*log(pi))
    }

  #One sample
  if(inputtype == "onesample"){
  div <- shannon(abund)
  return(div)
  }

  #Multiple samples
  if(inputtype == "multiplesamples"){
  divs <- apply(abund, 2, function(x) shannon(x))
  return(divs)
  }

}

#Simpson index
if(index == "simpson"){
  if(!missing(tree)) warning("Phylogenetic tree was not used to compute Simpson index.")
  #Function
  simpson <- function(vector){
    pi <- tss(vector[vector != 0])
    div <- 1 - (sum(pi^2))
    }

  #One sample
  if(inputtype == "onesample"){
  div <- simpson(abund)
  return(div)
  }

  #Multiple samples
  if(inputtype == "multiplesamples"){
  divs <- apply(abund, 2, function(x) simpson(x))
  return(divs)
  }
}

#Faith's PD
if(index == "faith"){
if(missing(tree)) stop("Faith's PD cannot be computed without a phylogenetic tree")
  #Function
  faith <- function(vector,tree){
  Li <- tree$edge.length
  ltips <- sapply(tree$edge[, 2], function(node) geiger::tips(tree, node))
  ai <- unlist(lapply(ltips, function(TipVector) sum(vector[TipVector])))
  Li <- Li[ai != 0]
  phylodiv <- sum(Li)
  }

  #One sample
  if(inputtype == "onesample"){
  div <- faith(abund,tree)
  return(div)
  }

  #Multiple samples
  if(inputtype == "multiplesamples"){
  divs <- apply(abund, 2, function(x) faith(x,tree))
  return(divs)
  }

}

#Allen's H
if(index == "allen"){
if(missing(tree)) stop("Allen's H cannot be computed without a phylogenetic tree")
  #Function
  allen <- function(vector,tree){
  vector <- tss(vector)
  Li <- tree$edge.length
  ltips <- sapply(tree$edge[, 2], function(node) geiger::tips(tree, node))
  ai <- unlist(lapply(ltips, function(TipVector) sum(vector[TipVector])))
  T <- sum(Li * ai)
  Li <- Li[ai != 0]
  ai <- ai[ai != 0]
  phylodiv <- -sum(Li*ai*log(ai))
  }

  #One sample
  if(inputtype == "onesample"){
  phylodiv <- allen(abund,tree)
  return(phylodiv)
  }

  #Multiple samples
  if(inputtype == "multiplesamples"){
  phylodivs <- apply(abund, 2, function(x) allen(x,tree))
  return(phylodivs)
  }

}

#Rao's Q
if(index == "rao"){
if(missing(tree)) stop("Rao's Q cannot be computed without a phylogenetic tree")
  #Function
  rao <- function(vector,tree){
  vector <- tss(vector)
  pij <- outer(vector, vector)
  phylodist <- cophenetic.phylo(tree)/2
  phylodist <- phylodist[rownames(pij),colnames(pij)]
  phylodiv <- sum(phylodist * pij)
  }

  #One sample
  if(inputtype == "onesample"){
  phylodiv <- rao(abund,tree)
  return(phylodiv)
  }

  #Multiple samples
  if(inputtype == "multiplesamples"){
  phylodivs <- apply(abund, 2, function(x) rao(x,tree))
  return(phylodivs)
  }

}

}
