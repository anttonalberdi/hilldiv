#' Diversity test
#' @title Diversity test
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords hill numbers comparison
#' @description Diversity comparison test between groups of samples. The function automatically assesses whether the data meets the properties for parametric statistics and performs the appropriate test accordingly: Students' T, ANOVA, Wilcoxon or Kruskal-Wallis. If the posthoc argument is set as TRUE, multiple group comparisons are complemented with post hoc pairwise tests, either Tukey test (parametric) or Dunn test with Benjamini-Hochberg correction (non-parametric).
#' @param countable A matrix indicating the relative abundances of multiple samples. Columns should be samples and rows OTUs.
#' @param qvalue A positive integer or decimal number (>=0), usually between 0 and 3.
#' @param hierarchy A two-column matrix indicating the relation between samples (first column) and groups (second column).
#' @param tree A phylogenetic tree of class 'phylo'. The tip labels must match the row names in the OTU table. Use the function match_data() if the OTU names do not match.
#' @param posthoc Whether to run post hoc pairwise analyses or not. If TRUE, an ANOVA will be complemented with a Tukey test and a Kruskal-Wallis test will be complemented with a Dunn test.
#' @return A statistical test output.
#' @seealso \code{\link{hill_div}}, \code{\link{div_part}}
#' @import stats
#' @importFrom FSA dunnTest
#' @examples
#' data(bat.diet.otutable)
#' data(bat.diet.tree)
#' data(bat.diet.hierarchy)
#' div_test(bat.diet.otutable,qvalue=0,hierarchy=bat.diet.hierarchy)
#' div_test(bat.diet.otutable,qvalue=1,hierarchy=bat.diet.hierarchy,tree=bat.diet.tree)
#' div_test(bat.diet.otutable,2,bat.diet.hierarchy,bat.diet.tree)
#' div_test(bat.diet.otutable,qvalue=1,hierarchy=bat.diet.hierarchy,posthoc=TRUE)
#' @references
#' Alberdi, A., Gilbert, M.T.P. (2019). A guide to the application of Hill numbers to DNA-based diversity analyses. Molecular Ecology Resources, 19, 804-817.\cr\cr
#' Chao, A., Chiu, C.H., & Jost, L. (2014). Unifying species diversity, phylogenetic diversity, functional diversity, and related similarity and differentiation measures through hill numbers. Annual Review of Ecology Evolution and Systematics, 45, 297-324.
#' @export

div_test <- function(countable,qvalue,hierarchy,tree,posthoc){
if(missing(countable)) stop("Count table is missing")
if(is.null(dim(countable)) == TRUE) stop("The OTU table is not a matrix")
if(dim(countable)[1] < 2) stop("The count table contains less than 2 OTUs/ASVs")
if(dim(countable)[2] < 2) stop("The count table contains less than 2 samples")
if(missing(qvalue)) stop("q value is missing")
if(qvalue < 0) stop("q value needs to be possitive (equal or higher than zero)")
if(missing(hierarchy)) stop("Hierarchy table is necessary to contrast groups of samples")
if(length(unique(hierarchy[,2])) == 1) stop("The hierarchy table contains a single group. The function requires at least two groups to contrast diversities.")
if(min(table(hierarchy[,2])) == 1) stop("All contrasting groups need to have at least 2 samples")
if(missing(posthoc)){posthoc=FALSE}

#Compute diversity values
if(missing(tree)){
div.values <- hill_div(countable,qvalue)
}else{
if(class(tree) != "phylo") stop("Tree needs to be an object of class Phylo")
if(identical(sort(rownames(countable)),sort(tree$tip.label)) == FALSE) stop("OTU/ASV names in the count table and tree do not match. Use match_data() to address this issue.")
div.values <- hill_div(countable,qvalue,tree)
}
colnames(hierarchy) <- c("Sample","Group")
div.values.groups <- merge(hierarchy,t(t(div.values)),by.y="row.names",by.x="Sample", sort = FALSE)
colnames(div.values.groups) <- c("Sample","Group","Value")
div.values.groups$Group <- as.factor(div.values.groups$Group)

#Data distribution (normality and homogeneity) assessment
shapiro <- shapiro.test(div.values.groups$Value)
barlett <- bartlett.test(Value ~ Group, data= div.values.groups)
if((shapiro$p.value >= 0.05) & (barlett$p.value >= 0.05)){
norm.homo=TRUE
}else{
norm.homo=FALSE
}

#Run statistical tests and output
if(length(unique(div.values.groups$Group)) == 2){
    if(norm.homo == TRUE){
    method <- "Student's t-Test"

    test <- t.test(Value ~ Group, data = div.values.groups)
    results <- list(data=div.values.groups,normality.pvalue=shapiro$p.value,homogeneity.pvalue=barlett$p.value,groups=length(unique(div.values.groups$Group)),method = "Student's t-Test",test=c(t=unname(test$statistic),df=unname(test$parameter),p.value=unname(test$p.value)))
    }else{
    test <- wilcox.test(Value ~ Group, data = div.values.groups)
    results <- list(data=div.values.groups,normality.pvalue=shapiro$p.value,homogeneity.pvalue=barlett$p.value,groups=length(unique(div.values.groups$Group)),method = "Wilcoxon Rank Sum Test",test=c(W=unname(test$statistic),df=unname(test$parameter),p.value=unname(test$p.value)))
    }
}else{
    if(norm.homo == TRUE){
    anova <- aov(Value ~ Group, data = div.values.groups)
    test <- summary(anova)
      if(posthoc == FALSE){
      results <- list(data=div.values.groups,normality.pvalue=shapiro$p.value,homogeneity.pvalue=barlett$p.value,groups=length(unique(div.values.groups$Group)),method = "ANOVA",test=test)
      }else{
      test.ph <- TukeyHSD(anova)$Group
      test.ph <- test.ph[,c(1,4)]
      colnames(test.ph)[2] <- c("P.adj")
      results <- list(data=div.values.groups,normality.pvalue=shapiro$p.value,homogeneity.pvalue=barlett$p.value,groups=length(unique(div.values.groups$Group)),method = "ANOVA",test=test, posthoc.method = "Tukey post-hoc test", posthoc=test.ph)
      }
    }else{
    test <- kruskal.test(Value ~ Group, data = div.values.groups)
      if(posthoc == FALSE){
      results <- list(data=div.values.groups,normality.pvalue=shapiro$p.value,homogeneity.pvalue=barlett$p.value,groups=length(unique(div.values.groups$Group)),method = "Kruskal-Wallis Test",test=c(chi.squared=unname(test$statistic),df=unname(test$parameter),p.value=unname(test$p.value)))
      }else{
      test.ph <- dunnTest(Value ~ Group, data= div.values.groups, method="bh")$res
      test.ph <- test.ph[,c(1,2,4)]
      rownames(test.ph) <- gsub(" - ","-",test.ph[,1])
      test.ph <- test.ph[,-1]
      results <- list(data=div.values.groups,normality.pvalue=shapiro$p.value,homogeneity.pvalue=barlett$p.value,groups=length(unique(div.values.groups$Group)),method = "Kruskal-Wallis Test",test=c(chi.squared=unname(test$statistic),df=unname(test$parameter),p.value=unname(test$p.value)),posthoc.method = "Dunn test with Benjamini-Hochberg correction",posthoc=test.ph)
      }
    }
}
return(results)
}
