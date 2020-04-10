#' Depth coverage assessment
#' @title Depth coverage assessment
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords coverage sequencing
#' @description Coverage of the estimated Hill numbers at different orders of diversity.
#' @param abund A vector or a matrix/data.frame indicating the relative abundances of one or multiple samples, respectively. If a matrix/data.frame is provided, columns must refer to samples and rows to OTUs.
#' @param qvalue A positive integer or decimal number (>=0), usually between 0 and 3.
#' @return A matrix with observed diversity,  estimated diversities and coverage
#' @seealso \code{\link{hill_div}}, \code{\link{depth_filt}}
#' @examples
#' data(bat.diet.otutable)
#' depth_cov(bat.diet.otutable,0)
#' depth_cov(bat.diet.otutable,qvalue=1)
#' @references
#' Chao, A. & Jost, L. (2015) Estimating diversity and entropy profiles via discovery rates of new species. Methods in Ecology and Evolution, 6, 873-882.\cr\cr
#' Jost, L. (2006). Entropy and diversity. Oikos, 113, 363-375.\cr\cr
#' Hill, M. O. (1973). Diversity and evenness: a unifying notation and its consequences. Ecology, 54, 427-432.
#' @export

depth_cov <- function(abund,qvalue){
if(missing(abund)) stop("Abundance data is missing")
if(missing(qvalue)) stop("q value is missing")
if(qvalue < 0) stop("q value needs to be possitive (equal or higher than zero)")
if(qvalue > 3) warning("Diversity estimations are not accurate when q>3")

#Obtain sequencing depth
if(is.null(dim(abund)) == TRUE){
Depth <- sum(abund)
}else{
Depth <- apply(abund, 2, FUN=sum)
}

#Compute Hill number
Observed <- hill_div(abund,qvalue)

#Estimate Hill number (diversity profile estimator derived by Chao and Jost 2015)
#Chao, A. & Jost, L. (2015) Estimating diversity and entropy profiles via discovery rates of new species. Methods in Ecology and Evolution, 6, 873-882.

ChaoJost = function(x,q){
  x = x[x>0]
  n = sum(x)
  f1 = sum(x==1)
  f2 = sum(x==2)
  p1 = ifelse(f2>0,2*f2/((n-1)*f1+2*f2),ifelse(f1>0,2/((n-1)*(f1-1)+2),1))
  Sub <- function(q){
    if(q==0){
      sum(x>0) + (n-1)/n*ifelse(f2>0, f1^2/2/f2, f1*(f1-1)/2)
    }
    else if(q==1){
      r <- 1:(n-1)
      A <- sum(x/n*(digamma(n)-digamma(x)))
      B <- ifelse(f1==0|p1==1,0,f1/n*(1-p1)^(1-n)*(-log(p1)-sum((1-p1)^r/r)))
      exp(A+B)
    }else if(abs(q-round(q))==0){
      A <- sum(exp(lchoose(x,q)-lchoose(n,q)))
      ifelse(A==0,NA,A^(1/(1-q)))
    }else {
      sort.data = sort(unique(x))
      tab = table(x)
      term = sapply(sort.data,function(z){
        k=0:(n-z)
        sum(choose(k-q,k)*exp(lchoose(n-k-1,z-1)-lchoose(n,z)))
      })
      r <- 0:(n-1)
      A = sum(tab*term)
      B = ifelse(f1==0|p1==1,0,f1/n*(1-p1)^(1-n)*(p1^(q-1)-sum(choose(q-1,r)*(p1-1)^r)))
      (A+B)^(1/(1-q))
    }
  }
  sapply(q, Sub)
}

if(is.null(dim(abund)) == TRUE){
Estimated <- ChaoJost(abund,qvalue)
}else{
Estimated <- apply(abund, 2, function(x) ChaoJost(x,qvalue))
}

#Compute coverage
Coverage <- Observed/Estimated*100
Coverage[Coverage > 100] <- 100

#Combine and print
result <- cbind(Depth,round(Observed,2),round(Estimated,2),round(Coverage,2))
colnames(result) <- c("Depth","Observed","Estimated","Coverage")
return(result)
}
