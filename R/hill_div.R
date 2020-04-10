#' Hill numbers computation
#' @title Hill numbers computation
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords OTU ASV diversity
#' @description Compute neutral or phylogenetic Hill numbers from a single sample (vector) or count table (matrix). Hill numbers or numbers equivalents of diversity indices are diversity measures that compute diversity in effective number of OTUs, i.e. the number of equally abundant OTUs that would be needed to give the same value of diversity.
#' @param count A vector or a matrix/data.frame indicating the (relative) counts of one or multiple samples, respectively. If a matrix/data.frame is provided, columns must refer to samples and rows to OTUs.
#' @param qvalue A positive integer or decimal number (>=0), usually between 0 and 3.
#' @param tree An ultrametic tree of class 'phylo'. The tip labels must match the names of the vector values (if one sample) or matrix rows (if multiple samples). Use the function match_data() if the OTU names do not match.
#' @param dist A dist object indicating the pairwise distances between samples. NOT implemented yet
#' @seealso \code{\link{index_div}}, \code{\link{div_part}}
#' @importFrom ape is.ultrametric
#' @importFrom geiger tips
#' @examples
#' data(bat.diet.otutable)
#' data(bat.diet.tree)
#' data(bat.diet.hierarchy)
#' #One sample
#' bat.diet.sample <- bat.diet.otutable[,1]
#' hill_div(bat.diet.sample,0)
#' hill_div(bat.diet.sample,qvalue=1)
#' #One sample (phylogenetic)
#' names(bat.diet.sample) <- rownames(bat.diet.otutable)
#' hill_div(bat.diet.sample,1,bat.diet.tree)
#' #Multiple samples
#' hill_div(bat.diet.otutable,0)
#' #Incidence-based
#' bat.diet.otutable.incidence <- to.incidence(bat.diet.otutable,bat.diet.hierarchy)
#' hill_div(bat.diet.otutable.incidence,qvalue=1)
#' hill_div(to.incidence(bat.diet.otutable,bat.diet.hierarchy),1)
#' @references
#' Alberdi, A., Gilbert, M.T.P. (2019). A guide to the application of Hill numbers to DNA-based diversity analyses. Molecular Ecology Resources, 19, 804-817.\cr\cr
#' Jost, L. (2006). Entropy and diversity. Oikos, 113, 363-375.\cr\cr
#' Hill, M. O. (1973). Diversity and evenness: a unifying notation and its consequences. Ecology, 54, 427-432.
#' @export

hill_div <- function(count,qvalue,tree,dist){

#Quality-check and warnings
if(missing(count)) stop("Count data is missing")
if(missing(qvalue)) stop("q value is missing")
if(qvalue==1){qvalue=0.99999}
if(qvalue < 0) stop("q value needs to be possitive (equal or higher than zero)")

#Hill numbers type definition
if(missing(tree) && missing(dist)){hilltype="neutral"}
if(!missing(tree) && missing(dist)){hilltype="phylogenetic"}
if(missing(tree) && !missing(dist)){hilltype="functional"}
if(!missing(tree) && !missing(dist)) stop("Phylogenetic and trait information cannot be added at once. Use either phylogenetic (tree) or functional (traits) information.")

#Input data type definition
if(is.null(dim(count)) == TRUE){
inputtype="onesample"
}else{
inputtype="multiplesamples"
}

###
# Neutral Hill numbers
###
if(hilltype == "neutral"){

#Function
neutral.Hill <- function(vector,qvalue){
    pi <- vector[vector!=0]
    sum(pi^qvalue)^(1/(1-qvalue))
    }

#One sample
if(inputtype == "onesample"){
div <- neutral.Hill(tss(count),qvalue)
return(div)
}

#Multiple samples
if(inputtype == "multiplesamples"){
divs <- apply(tss(count), 2, function(x) neutral.Hill(x,qvalue))
return(divs)
}

}

###
# Phylogenetic Hill numbers
###
if(hilltype == "phylogenetic"){

if(class(tree) != "phylo") stop("Tree needs to be an object of class Phylo")
if(is.ultrametric(tree) == FALSE) stop("Tree needs to be ultrametric")

#Function
phylogenetic.Hill <- function(vector,qvalue,tree){
  Li <- tree$edge.length
  ltips <- sapply(tree$edge[, 2], function(node) tips(tree, node))
  ai <- unlist(lapply(ltips, function(TipVector) sum(vector[TipVector])))
  T <- sum(Li * ai)
  Li <- Li[ai != 0]
  ai <- ai[ai != 0]
  sum(Li/T * ai^qvalue)^(1/(1-qvalue))
  }

  #One sample
  if(inputtype == "onesample"){
  if(identical(sort(names(count)),sort(tree$tip.label)) == FALSE) stop("OTU names in the vector and tree do not match. Use match_data() to solve this issue.")
  div <- phylogenetic.Hill(tss(count),qvalue,tree)
  return(div)
  }

  #Multiple samples
  if(identical(sort(rownames(count)),sort(tree$tip.label)) == FALSE) stop("OTU names in the vector and tree do not match. Use match_data() to solve this issue.")
  if(inputtype == "multiplesamples"){
  divs <- apply(tss(count), 2, function(x) phylogenetic.Hill(x,qvalue,tree))
  return(divs)
  }

}

###
# Functional Hill numbers
###
if(hilltype == "functional") stop("Functional Hill numbers have not been implemented yet in the function.")

#End of function
}
