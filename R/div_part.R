#' Multi-level diversity partitioning (based on Hill numbers)
#' @title Multi-level diversity partitioning
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords alpha gamma beta hill
#' @description Multi-level diversity partitioning following the multiplicative definition based on Hill numbers. Hierarchical levels are defined from L1 (minimum, sample) to Ln (maximum, whole system), and as many intermediate levels as wanted can be defined in between. The hierarchical structure of the system is defined with the hierarchy table. If no hierarchy table is inputed, the function yields a simple two-level partitioning between alpha (L1), beta and gamma (L2).
#' @param countable A count table (matrix/data.frame) indicating the absolute or relative OTU/ASV abundances of multiple samples. Columns must refer to samples and rows to OTUs/ASVs.
#' @param qvalue A positive number, usually between 0 and 5, but most commonly 0, 1 or 2. It can be an integer or contain decimals.
#' @param tree A phylogenetic tree of class 'phylo'. The tip labels must match the row names in the OTU table. Use the function match_data() if the OTU names do not match.
#' @param hierarchy A matrix indicating the relation between samples (first column) and parent group(s).
#' @return A list object containing details of hierarchical diversity partitioning.
#' @seealso \code{\link{div_part}}, \code{\link{gamma_div}}, \code{\link{match_data}}
#' @importFrom stats aggregate
#' @examples
#' data(bat.diet.otutable)
#' data(bat.diet.tree)
#' data(bat.diet.hierarchy)
#' #Two level examples (L1=sample (alpha diversity), L2=whole system (gamma diversity))
#' div_part(bat.diet.otutable,qvalue=1)
#' div_part(bat.diet.otutable,qvalue=0,tree=bat.diet.tree)
#' #Three-level example (L1=sample, L2=species, L3=whole system)
#' div_part(bat.diet.otutable,qvalue=0,hierarchy=bat.diet.hierarchy)
#' @references
#' Alberdi, A., Gilbert, M.T.P. (2019). A guide to the application of Hill numbers to DNA-based diversity analyses. Molecular Ecology Resources, 19, 804-817.\cr\cr
#' Chao, A., Chiu, C.H., & Hsieh, T. C. (2012). Proposing a resolution to debates on diversity partitioning. Ecology, 93, 2037-2051.\cr\cr
#' Jost, L. (2007). Partitioning diversity into independent alpha and beta components. Ecology, 88, 2427–2439.
#' @export

div_part <- function(countable,qvalue,tree,hierarchy) {

#Quality-check and warnings
if(missing(countable)) stop("OTU table is missing")
if(is.null(dim(countable)) == TRUE) stop("The OTU table is not a matrix")
if(dim(countable)[1] < 2) stop("The OTU table only less than 2 OTUs")
if(dim(countable)[2] < 2) stop("The OTU table contains less than 2 samples")
if(sum(colSums(countable)) != ncol(countable)) {countable <- tss(countable)}
if(missing(qvalue)) stop("q value is missing")
if(qvalue < 0) stop("q value needs to be possitive (equal or higher than zero)")
if(qvalue==1) {qvalue=0.99999}

#####
# 2-level hierarchy (no hierarchy table)
#####

#Compute L1 and L2 diversities
if(missing(hierarchy)){
if(missing(tree)){
  L1_div <- alpha_div(countable,qvalue)
  L2_div <- gamma_div(countable,qvalue)
  }else{
  L1_div <- alpha_div(countable,qvalue,tree)
  L2_div <- gamma_div(countable,qvalue,tree)
}
div.vector <- c(L1_div,L2_div)
names(div.vector) <- c("L1","L2")

#Beta
beta <- L2_div/L1_div

#Sample size
N <- c(N1=ncol(countable),N2=1)

#Return values
if (qvalue==0.99999) {qvalue=1}
results <- list("Hierarchical_levels" = 2, "Order_diversity" = qvalue, "Hill_numbers" = div.vector, "Sample_size" = N, "Beta" = beta)
return(results)
}

#####
# Multi-level hierarchy (with hierarchy table)
#####
if(!missing(hierarchy)){

#Check nestedness of hierarchy
if(is.nested(hierarchy) == FALSE) stop("The groups in the hierarchy table are not nested.")

#Count number of levels
leveln <- ncol(hierarchy)
levels <- paste(rep("L",leveln+1),seq(1:(leveln+1)),sep="")

#Convert hierarchy columns to character
hierarchy[] <- lapply(hierarchy, as.character)
colnames(hierarchy) <- levels[-length(levels)]

#Generate aggregated OTU tables
count.tables <- list()
count.tables[[1]] <- countable
countable.sub <- countable
for(i in c(2:leveln)){
  countable.sub <- merge(t(countable.sub),unique(hierarchy[,c(i-1,i)]), by.x="row.names",by.y=as.character(levels[i-1]))
  countable.sub <- countable.sub[,-1]
  countable.sub <- aggregate(subset(countable.sub, select=rownames(countable)), by=list(countable.sub[,as.character(levels[i])]), FUN=sum)
  rownames(countable.sub) <- countable.sub[,1]
  countable.sub <- t(countable.sub[,-1])
  count.tables[[i]] <- countable.sub
}

#Generate vector of diversities at different hierarchical levels
div.vector <- c()
for (i in c(1:(leveln+1))){
if(i == 1){
  #If lowest level
  if(missing(tree)){
    div.vector <- c(div.vector,alpha_div(count.tables[[1]],qvalue))
    }else{
    div.vector <- c(div.vector,alpha_div(count.tables[[1]],qvalue,tree))
  }
}else if(i == leveln+1){
  #If highest level
  if(missing(tree)){
  div.vector <- c(div.vector,gamma_div(count.tables[[1]],qvalue))
  }else{
  div.vector <- c(div.vector,gamma_div(count.tables[[1]],qvalue,tree))
  }
}else{
  #Intermediate level
  if(missing(tree)){
  div.vector <- c(div.vector,alpha_div(count.tables[[i]],qvalue))
  }else{
  div.vector <- c(div.vector,alpha_div(count.tables[[i]],qvalue,tree))
  }
}
}

#Name levels
names(div.vector) <- levels

#Get beta values
beta.vector <- c()
N.vector <- c()
for(b in c(1:(leveln))){
beta <- div.vector[b+1]/div.vector[b]
names(beta) <- paste("B",paste(b,b+1,sep="_"),sep="")
beta.vector <- c(beta.vector,beta)
N <- ncol(count.tables[[b]])
N.vector <- c(N.vector,N)
}

N.vector <- c(N.vector,1)
names(N.vector) <- paste(rep("N",leveln+1),seq(1:(leveln+1)),sep="")

#Return values
if (qvalue==0.99999) {qvalue=1}
results <- list("Hierarchical_levels" = (leveln+1), "Order_diversity" = qvalue, "Hill_numbers"=div.vector, "Sample_size"=N.vector, "Beta"=beta.vector)
return(results)

}
}
