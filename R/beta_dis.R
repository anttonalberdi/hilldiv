#' (Dis)similarity computation from beta diversities based on Hill numbers
#' @title Beta dissimilarity
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords alpha beta gamma partitioning hill similarity dissimilarity
#' @description Compute dissimilarity or similarity values based on beta diversities (neutral or phylogenetic) and sample size.
#' @param beta A numeric beta diversity value or an object outputted by function div_part() (which contains all the information to compute (dis)similarities).
#' @param qvalue A positive number, usually between 0 and 5, but most commonly 0, 1 or 2. It can be an integer or contain decimals.
#' @param N An integer indicating sample size, the number of sampling units to be used to compute the (dis)similarity measures. The argument is ovewritten if a 'div_part' object is used.
#' @param metric A vector containing "C", "U", "V" or "S". C: Sørensen-type overlap or complement. U: Jaccard-type overlap or complement. V: Sørensen-type turnover or complement. S: Jaccard-type turnover or complement. See hilldiv wiki for further information.
#' @param type A character object containing either "similarity" or "dissimilarity". If 'similarity' is used, similarity metrics (0: completely different composition - 1: identical composition) are returned. If 'dissimilarity' is used, dissimilarity metrics (0: identical composition - 1:completely different composition) are returned.
#' @seealso \code{\link{div_part}}, \code{\link{gamma_div}}, \code{\link{pair_dis}}
#' @examples
#' data(bat.diet.otutable)
#' data(bat.diet.tree)
#' #Manually indicating beta diversity, order of diversity and sample size
#' beta_dis(beta=4.5,qvalue=1,N=8)
#' beta_dis(beta=4.5,qvalue=1,N=8,metric="C",type="similarity")
#' #Using an object created with the function div_part()
#' divpartobject <- div_part(bat.diet.otutable,qvalue=0,tree=bat.diet.tree)
#' beta_dis(divpartobject)
#' beta_dis(divpartobject,metric="S",type="similarity")
#' @references
#' Alberdi, A., Gilbert, M.T.P. (2019). A guide to the application of Hill numbers to DNA-based diversity analyses. Molecular Ecology Resources, 19, 804-817.\cr\cr
#' Chao, A., Chiu, C.H., & Hsieh, T. C. (2012). Proposing a resolution to debates on diversity partitioning. Ecology, 93, 2037-2051.\cr\cr
#' Jost, L. (2007). Partitioning diversity into independent alpha and beta components. Ecology, 88, 2427-2439.
#' @export

beta_dis <- function(beta,qvalue,N,metric,type){

#Quality-check and warnings
if(missing(beta)) stop("Beta diversity value or div_part output object is missing")
if(class(beta) == "numeric"){
  betan <- length(beta)
  betas <- beta
  Ns <- N
  if(missing(qvalue)) stop("The order of diversity (q) is missing")
  if (qvalue==1) {qvalue=0.99999}
  if(missing(N)) stop("The number of samples or groups (N) is missing")
}

if(class(beta) == "list"){
  if(names(beta[2]) != "Order_diversity") stop("The input object is not valid")
  betan <- length(beta$Beta)
  qvalue <- beta$Order_diversity
  if (qvalue==1) {qvalue=0.99999}
  betas <- beta$Beta
  Ns <- beta$Sample_size[-length(beta$Sample_size)]
}
if(missing(metric)) {metric= c("C","U","V","S")}
if(missing(type)) {type="dissimilarity"}

###### MULTIPLE HIERARCHIES NEED TO BE ADDED, and similarity functions updated
results <- list()

#Sørensen-type overlap (CqN, 1-CqN)
if ('C' %in% metric){
  CqNs <- c()
  for(i in c(1:betan)){
  CqNs <- c(CqNs,CqN(betas[i],qvalue,Ns[i]))
  names(CqNs)[i] <- paste("L",paste(i,i+1,sep="_"),sep="")
  }
  if (type == "dissimilarity"){
    rCqNs <- 1 - CqNs
    results <- append(results, list(CqN=rCqNs))
  }else{
    results <- append(results, list(CqN=CqNs))
  }
}

#Jaccard-type overlap (UqN, 1-UqN)
if ('U' %in% metric){
  UqNs <- c()
  for(i in c(1:betan)){
  UqNs <- c(UqNs,UqN(betas[i],qvalue,Ns[i]))
  names(UqNs)[i] <- paste("L",paste(i,i+1,sep="_"),sep="")
  }
  if (type == "dissimilarity"){
    rUqNs <- 1 - UqNs
    results <- append(results, list(UqN=rUqNs))
  }else{
    results <- append(results, list(UqN=UqNs))
  }
}

#Sørensen-type turnover-complement (VqN, 1-VqN)
if ('V' %in% metric){
  VqNs <- c()
  for(i in c(1:betan)){
  VqNs <- c(VqNs,VqN(betas[i],Ns[i]))
  names(VqNs)[i] <- paste("L",paste(i,i+1,sep="_"),sep="")
  }
  if (type == "dissimilarity"){
    rVqNs <- 1 - VqNs
    results <- append(results, list(VqN=rVqNs))
  }else{
    results <- append(results, list(VqN=VqNs))
  }
}

#Jaccard-type turnover-complement (SqN, 1-SqN)
if ('S' %in% metric){
  SqNs <- c()
  for(i in c(1:betan)){
  SqNs <- c(SqNs,SqN(betas[i],Ns[i]))
  names(SqNs)[i] <- paste("L",paste(i,i+1,sep="_"),sep="")
  }
  if (type == "dissimilarity"){
    rSqNs <- 1 - SqNs
    results <- append(results, list(SqN=rSqNs))
  }else{
    results <- append(results, list(SqN=SqNs))
  }
}

return(results)

}
