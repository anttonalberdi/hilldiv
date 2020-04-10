#' Diversity profile
#' @title Diversity profile
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords alpha gamma beta hill
#' @description Create diversity profile vectors (single sample or system) or tables (multiple samples or groups) from count tables.
#' @param count A vector or a matrix indicating the (relative) OTU/ASV counts of one or multiple samples. If a matrix is provided, columns must refer to samples and rows to OTUs.
#' @param qvalues A vector of sequential orders of diversity (default from 0 to 5). qvalue=seq(from = 0, to = 5, by = (0.1))
#' @param tree A tree of class 'phylo'. The tip labels must match the names of the vector values (if one sample) or matrix rows (if multiple samples).
#' @param hierarchy A two-column matrix indicating the relation between samples (first column) and groups (second column).
#' @param level Whether to compute alpha or gamma diversities of the system or the groups specified in the hierarchy table.
#' @return A vector or matrix containing diversity values at different orders of diversity (as specified in qvalues).
#' @seealso \code{\link{div_profile_plot}}, \code{\link{hill_div}}
#' @importFrom ape drop.tip is.ultrametric
#' @importFrom geiger tips
#' @examples
#' data(bat.diet.otutable)
#' data(bat.diet.tree)
#' data(bat.diet.hierarchy)
#' #One sample example
#' bat.diet.sample <- bat.diet.otutable[,1]
#' div_profile(count=bat.diet.sample,qvalues=seq(from = 0, to = 5, by = (0.1)))
#' #One sample example (phylogenetic Hill numbers)
#' names(bat.diet.sample) <- rownames(bat.diet.otutable)
#' div_profile(count=bat.diet.sample,qvalues=seq(from = 0, to = 5, by = (0.1)),tree=bat.diet.tree)
#' #Multiple samples
#' div_profile(bat.diet.otutable)
#' #Multiple groups (gamma diversity)
#' div_profile(bat.diet.otutable,hierarchy=bat.diet.hierarchy,level="gamma")
#' #Multiple groups (alpha diversity)
#' div_profile(bat.diet.otutable,hierarchy=bat.diet.hierarchy,level="alpha")
#' @references
#' Alberdi, A., Gilbert, M.T.P. (2019). A guide to the application of Hill numbers to DNA-based diversity analyses. Molecular Ecology Resources, 19, 804-817.\cr\cr
#' Chao, A., Chiu, C.H., & Jost, L. (2014). Unifying species diversity, phylogenetic diversity, functional diversity, and related similarity and differentiation measures through hill numbers. Annual Review of Ecology Evolution and Systematics, 45, 297-324.
#' @export

div_profile <- function(count,qvalues,tree,hierarchy,level){

#Quality-check and warnings
if(missing(count)) stop("The countance data is missing")
if(missing(qvalues)) {qvalues= seq(from = 0, to = 5, by = (0.1))}
if(missing(level)) {level= "NA"}

#If input data is a single sample (vector)
if(is.null(dim(count)) == TRUE){
    profile <- c()
    for (o in qvalues){
    if(missing(tree)){
        div.value <- hill_div(count,o)
        }else{
        div.value <- hill_div(count,o,tree)
    }
    profile <- c(profile,div.value)
    }
    names(profile) <- qvalues
}

#Declare the fast version of hill_div if tree is inputed
if(!missing(tree)){

hill_div_fast <- function(count,qvalue,tree,dist){
    if(qvalue==1){qvalue=0.99999}
    phylogenetic.Hill.fast <- function(vector,qvalue,tree){
    Li <- tree$edge.length
    ai <- unlist(lapply(ltips, function(TipVector) sum(vector[TipVector])))
    T <- sum(Li * ai)
    Li <- Li[ai != 0]
    ai <- ai[ai != 0]
    sum(Li/T * ai^qvalue)^(1/(1-qvalue))
    }
    divs <- apply(tss(count), 2, function(x) phylogenetic.Hill.fast(x,qvalue,tree))
    return(divs)
}

alpha_div_fast <- function(otutable,qvalue,tree,weight){
if(missing(otutable)) stop("OTU table is missing")
if(is.null(dim(otutable)) == TRUE) stop("The OTU table is not a matrix")
if(dim(otutable)[1] < 2) stop("The OTU table only less than 2 OTUs")
if(dim(otutable)[2] < 2) stop("The OTU table contains less than 2 samples")
if(sum(colSums(otutable)) != ncol(otutable)) {otutable <- tss(otutable)}
if(missing(qvalue)) stop("q value is missing")
if(qvalue < 0) stop("q value needs to be possitive (equal or higher than zero)")
if(identical(sort(rownames(otutable)),sort(tree$tip.label)) == FALSE) stop("OTU names in the OTU table and tree do not match")
if(is.ultrametric(tree) == FALSE) stop("Tree needs to be ultrametric")
if (qvalue==1) {qvalue=0.99999}
if(missing(weight)) { weight= rep(1/ncol(otutable),ncol(otutable))}
otutable <- as.data.frame(otutable)
wj <- weight
N <- ncol(otutable)
Li <- tree$edge.length
aij <- matrix(unlist(lapply(ltips, function(TipVector) colSums(otutable[TipVector,]))), ncol = N, byrow = TRUE)
aij.wj <- sweep(aij, 2, wj, "*")
T <- sum(sweep(aij.wj, 1, Li, "*"))
L <- matrix(rep(Li, N), ncol = N)
wm <-  matrix(rep(wj, length(Li)), ncol = N,byrow=TRUE)
i <-  which(aij > 0)
phylodiv <- sum(L[i] * (aij[i]*wm[i]/T)^qvalue)^(1/(1 - qvalue))/(N*T)
return(phylodiv)
}

gamma_div_fast <- function(otutable,qvalue,tree,weight){
if(missing(otutable)) stop("OTU table is missing")
if(is.null(dim(otutable)) == TRUE) stop("The OTU table is not a matrix")
if(dim(otutable)[1] < 2) stop("The OTU table only less than 2 OTUs")
if(dim(otutable)[2] < 2) stop("The OTU table contains less than 2 samples")
if(sum(colSums(otutable)) != ncol(otutable)) {otutable <- tss(otutable)}
if(missing(qvalue)) stop("q value is missing")
if(qvalue < 0) stop("q value needs to be possitive (equal or higher than zero)")
if (qvalue==1) {qvalue=0.99999} # change q to the limit of the unity (0.99999) if q=1
if(is.ultrametric(tree) == FALSE) stop("Tree needs to be ultrametric")
if(identical(sort(rownames(otutable)),sort(tree$tip.label)) == FALSE) stop("OTU names in the OTU table and tree do not match")
if(missing(weight)) { weight= rep(1/ncol(otutable),ncol(otutable))}
otutable <- as.data.frame(otutable)
wj <- weight
N <- ncol(otutable)
Li <- tree$edge.length
aij <- matrix(unlist(lapply(ltips, function(TipVector) colSums(otutable[TipVector,]))), ncol = N, byrow = TRUE)
aij.wj <- sweep(aij, 2, wj, "*")
ai <- rowSums(aij.wj)
T <- sum(sweep(aij.wj, 1, Li, "*"))
L <- matrix(rep(Li, N), ncol = N)
Li <- Li[ai != 0] #Remove zeros
ai <- ai[ai != 0] #Remove zeros
wm <-  matrix(rep(wj, length(Li)), ncol = N, byrow=TRUE)
phylodiv <- (sum(Li * (ai/T)^qvalue)^(1/(1 - qvalue)))/T
return(phylodiv)
}

}

#If input data is a count table (matrix)
if(is.null(dim(count)) == FALSE){

    if(dim(count)[1] < 2) stop("The OTU table only less than 2 OTUs")
    if(dim(count)[2] < 2) stop("The OTU table contains less than 2 samples")

    #Without hierarchy
    if(missing(hierarchy)){
    profile <- c()
        for (o in qvalues){
            if(missing(tree)){
            if(level == "NA"){div.values <- hill_div(count,o)}
            if(level == "gamma"){div.values <- gamma_div(count,o)}
            if(level == "alpha"){div.values <- alpha_div(count,o)}
            }else{
            ltips <- sapply(tree$edge[, 2], function(node) tips(tree, node))
            if(level == "NA"){div.values <- hill_div_fast(count,o,tree)}
            if(level == "gamma"){div.values <- gamma_div(count,o,tree)}
            if(level == "alpha"){div.values <- alpha_div(count,o,tree)}
            }
        profile <- rbind(profile,div.values)
        }
        if(level == "NA"){
          rownames(profile) <- qvalues
          }
        if(level == "gamma"){
          profile <- c(profile)
          names(profile) <- qvalues
          }
        if(level == "alpha"){
          profile <- c(profile)
          names(profile) <- qvalues
          }

    }

    #With hierarchy
    if(!missing(hierarchy)){
    if(ncol(hierarchy) != 2) stop("The hierarchy table must contain two columns.")
    colnames(hierarchy) <- c("Sample","Group")
    groups <- as.character(sort(unique(hierarchy$Group)))
    profile <- c()
        for (g in groups){
            samples <- as.character(hierarchy[which(hierarchy$Group == g),1])
            count.subset <- count[,samples]
            count.subset <- as.data.frame(count.subset[apply(count.subset, 1, function(z) !all(z==0)),])

            #Subset tree
            if(!missing(tree)){
            missing.otus <- setdiff(tree$tip.label,rownames(count.subset))
            tree.subset <- drop.tip(tree,missing.otus)
            }

            for (o in qvalues){
                  if(missing(tree)){
                        if(level == "NA"){div.value <- gamma_div(count.subset,o)}
                        if(level == "gamma"){div.value <- gamma_div(count.subset,o)}
                        if(level == "alpha"){div.value <- alpha_div(count.subset,o)}
                  }else{
                        ltips <- sapply(tree.subset$edge[, 2], function(node) tips(tree.subset, node))
                        if(level == "NA"){div.value <- gamma_div_fast(count.subset,o,tree.subset)}
                        if(level == "gamma"){div.value <- gamma_div_fast(count.subset,o,tree.subset)}
                        if(level == "alpha"){div.value <- alpha_div_fast(count.subset,o,tree.subset)}
                        }
            profile <- rbind(profile,as.numeric(div.value))
            }
      }

    profile <- matrix(profile,nrow=length(qvalues))
    colnames(profile) <- as.character(groups)
    rownames(profile) <- as.character(qvalues)
  }
}
return(profile)
}
