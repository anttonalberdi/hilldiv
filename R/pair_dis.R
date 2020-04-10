#' Pairwise dissimilarity
#' @title Pairwise dissimilarity
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords hill numbers diversity partitioning beta
#' @description Computation of pairwise dissimilarities based on Hill numbers diversity partitioning
#' @param countable A matrix indicating the relative abundances of multiple samples. Columns should be samples and rows OTUs.
#' @param qvalue A positive integer or decimal number (>=0), usually between 0 and 3.
#' @param tree A phylogenetic tree of class 'phylo'. The tip labels must match the row names in the OTU table. Use the function match_data() if the OTU names do not match.
#' @param hierarchy A matrix indicating the relation between samples (first column) and groups.
#' @param metric A vector containing any combination of "C", "U", "V" or "S". If not provided, all metrics will be computed. metric="U", metric=c("U","S").
#' @return A list of matrices containing pairwise beta diversities and dissimilarity metrics.
#' @seealso \code{\link{hill_div}}, \code{\link{div_part}}, \code{\link{beta_dis}}
#' @examples
#' data(bat.diet.otutable)
#' data(bat.diet.tree)
#' data(bat.diet.hierarchy)
#' pair_dis(bat.diet.otutable,qvalue=1)
#'\donttest{
#' pair_dis(bat.diet.otutable,qvalue=1,tree=bat.diet.tree,metric="V")
#'}
#' @references
#' Alberdi, A., Gilbert, M.T.P. (2019). A guide to the application of Hill numbers to DNA-based diversity analyses. Molecular Ecology Resources, 19, 804-817.\cr\cr
#' Chao, A., Chiu, C.H., & Hsieh, T. C. (2012). Proposing a resolution to debates on diversity partitioning. Ecology, 93, 2037-2051.\cr\cr
#' Jost, L. (2007). Partitioning diversity into independent alpha and beta components. Ecology, 88, 2427-2439.
#' @export

pair_dis <- function(countable,qvalue,tree,hierarchy,metric){

#Quality-check and warnings
if(missing(countable)) stop("OTU table is missing")
if(is.null(dim(countable)) == TRUE) stop("The OTU table is not a matrix")
if(dim(countable)[1] < 2) stop("The OTU table only less than 2 OTUs")
if(dim(countable)[2] < 2) stop("The OTU table contains less than 2 samples")
if(sum(colSums(countable)) != ncol(countable)) {countable <- tss(countable)}
if(missing(qvalue)) stop("q value is missing")
if(qvalue < 0) stop("q value needs to be possitive (equal or higher than zero)")
if (qvalue==1) {qvalue=0.99999}
if(missing(metric)) { metric= c("C","U","V","S")}

#Declare fast alpha and gamma phylodiversities (without ltips, as it is the same for all combinations)

alpha_div_fast <- function(otutable,qvalue,tree){
if(missing(otutable)) stop("OTU table is missing")
if(is.null(dim(otutable)) == TRUE) stop("The OTU table is not a matrix")
if(dim(otutable)[1] < 2) stop("The OTU table only less than 2 OTUs")
if(dim(otutable)[2] < 2) stop("The OTU table contains less than 2 samples")
if(missing(qvalue)) stop("q value is missing")
if(qvalue < 0) stop("q value needs to be possitive (equal or higher than zero)")
if(identical(sort(rownames(otutable)),sort(tree$tip.label)) == FALSE) stop("OTU names in the OTU table and tree do not ")
if(ape::is.ultrametric(tree) == FALSE) stop("Tree needs to be ultrametric")
if (qvalue==1) {qvalue=0.99999}
otutable <- as.data.frame(otutable)
wj <- rep(1/ncol(otutable),ncol(otutable))
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

gamma_div_fast <- function(otutable,qvalue,tree){
if(missing(otutable)) stop("OTU table is missing")
if(is.null(dim(otutable)) == TRUE) stop("The OTU table is not a matrix")
if(dim(otutable)[1] < 2) stop("The OTU table only less than 2 OTUs")
if(dim(otutable)[2] < 2) stop("The OTU table contains less than 2 samples")
if(sum(colSums(otutable)) != ncol(otutable)) {otutable <- tss(otutable)}
if(missing(qvalue)) stop("q value is missing")
if(qvalue < 0) stop("q value needs to be possitive (equal or higher than zero)")
if (qvalue==1) {qvalue=0.99999} # change q to the limit of the unity (0.99999) if q=1
if(ape::is.ultrametric(tree) == FALSE) stop("Tree needs to be ultrametric")
if(identical(sort(rownames(otutable)),sort(tree$tip.label)) == FALSE) stop("OTU names in the OTU table and tree do not match")
otutable <- as.data.frame(otutable)
wj <- rep(1/ncol(otutable),ncol(otutable))
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

#Generate ltips
if(!missing(tree)){
ltips <- sapply(tree$edge[, 2], function(node) geiger::tips(tree, node))
}

#Count number of levels
if(!missing(hierarchy)){
leveln <- ncol(hierarchy)
}else{
leveln <- 1
}
levels <- paste(rep("L",leveln),seq(1:leveln),sep="")
if(!missing(hierarchy)){
  hierarchy[] <- lapply(hierarchy, as.character)
  colnames(hierarchy) <- levels
}

#Generate aggregated OTU tables
count.tables <- list()
count.tables[[1]] <- countable
if(leveln > 1){
  countable.sub <- countable
  for(i in c(2:leveln)){
    countable.sub <- merge(t(countable.sub),unique(hierarchy[,c(i-1,i)]), by.x="row.names",by.y=as.character(levels[i-1]))
    countable.sub <- countable.sub[,-1]
    countable.sub <- aggregate(subset(countable.sub, select=rownames(countable)), by=list(countable.sub[,as.character(levels[i])]), FUN=mean)
    rownames(countable.sub) <- countable.sub[,1]
    countable.sub <- t(countable.sub[,-1])
    count.tables[[i]] <- countable.sub
  }
}

#Generate results
results <- list()
names <- c()
for (i in c(1:leveln)){
#Generate matrices
count.table.sub <- count.tables[[i]]
indices <- sort(colnames(count.table.sub))

beta.matrix <- matrix(rep(NA,length(indices)^2), nrow = length(indices), ncol = length(indices))
colnames(beta.matrix) <- indices
rownames(beta.matrix) <- indices
if('C' %in% metric){ CqN.matrix <- beta.matrix }
if('U' %in% metric){ UqN.matrix <- beta.matrix }
if('V' %in% metric){ VqN.matrix <- beta.matrix }
if('S' %in% metric){ SqN.matrix <- beta.matrix }

#Populate matrices
for (x in indices){
for (y in indices){
if(is.na(beta.matrix[x,y])){ #to avoid repeating mirror operations
combination <- count.table.sub[,c(y,x)]

if(identical(x,y) == TRUE){
    beta <- NA
}else{
    if(missing(tree)){
    alpha <- alpha_div(combination,qvalue)
    gamma <- gamma_div(combination,qvalue)
    }else{
    alpha <- alpha_div_fast(combination,qvalue,tree)
    gamma <- gamma_div_fast(combination,qvalue,tree)
    }
    beta <- gamma/alpha

    beta.matrix[y,x] <- beta

    if('C' %in% metric){
    CqN.matrix[y,x] <- beta_dis(beta=beta,qvalue=qvalue,N=2,metric="C",type="dissimilarity")$CqN
    }

    if('U' %in% metric){
    UqN.matrix[y,x] <- beta_dis(beta=beta,qvalue=qvalue,N=2,metric="U",type="dissimilarity")$UqN
    }

    if('V' %in% metric){
    VqN.matrix[y,x] <- beta_dis(beta=beta,qvalue=qvalue,N=2,metric="V",type="dissimilarity")$VqN
    }

    if('S' %in% metric){
    SqN.matrix[y,x] <- beta_dis(beta=beta,qvalue=qvalue,N=2,metric="S",type="dissimilarity")$SqN
    }

}
}
}
}

#Append matrices to results
results <- append(results, list(Beta=beta.matrix))
if('C' %in% metric){results <- append(results, list(CqN=CqN.matrix))}
if('U' %in% metric){results <- append(results, list(UqN=UqN.matrix))}
if('V' %in% metric){results <- append(results, list(VqN=VqN.matrix))}
if('S' %in% metric){results <- append(results, list(SqN=SqN.matrix))}

#Append matrix names
names <- c(names,paste(paste("L",i,sep=""),"beta",sep="_"))
if('C' %in% metric){names <- c(names,paste(paste("L",i,sep=""),"CqN",sep="_"))}
if('U' %in% metric){names <- c(names,paste(paste("L",i,sep=""),"UqN",sep="_"))}
if('V' %in% metric){names <- c(names,paste(paste("L",i,sep=""),"VqN",sep="_"))}
if('S' %in% metric){names <- c(names,paste(paste("L",i,sep=""),"SqN",sep="_"))}
}

#Modify names
names(results) <- names

return(results)

}
