#' Merge samples
#' @title Merge samples
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords OTU ASV table combine
#' @description Combines samples into groups defined by the hierarchy table, with the possibility to convert abundances into incidence data.
#' @param countable A count table (matrix/data.frame) indicating the absolute or relative OTU/ASV abundances of multiple samples. Columns must refer to samples and rows to OTUs/ASVs.
#' @param hierarchy A two-column matrix indicating the relation between samples (first column) and groups (second column).
#' @param relative Whether to output relative values or not. Default=TRUE.
#' @param incidence Whether to transform abundance into incidence data when merging. Default=FALSE.
#' @usage merge_samples(countable,hierarchy,incidence)
#' @return A count table
#' @seealso \code{\link{to.incidence}}
#' @examples
#' data(bat.diet.otutable)
#' data(bat.diet.hierarchy)
#' merge_samples(countable=bat.diet.otutable,hierarchy=bat.diet.hierarchy)
#' merge_samples(bat.diet.otutable,bat.diet.hierarchy)
#' merge_samples(countable=bat.diet.otutable,hierarchy=bat.diet.hierarchy, incidence=TRUE)
#' @references
#' Alberdi, A., Gilbert, M.T.P. (2019). A guide to the application of Hill numbers to DNA-based diversity analyses. Molecular Ecology Resources, 19, 804-817.\cr\cr
#' @export

merge_samples <- function(countable,hierarchy,relative,incidence){
if(missing(countable)) stop("Count table is missing")
if(is.null(dim(countable)) == TRUE) stop("The OTU table is not a matrix")
if(missing(hierarchy)) stop("A hierarchy table is necessary to merge samples")
if(length(unique(hierarchy[,2])) == 1) stop("The hierarchy table contains a single group. The function requires at least two groups to merge samples.")
if(missing(relative)){relative=TRUE}
if(missing(incidence)){incidence=FALSE}

colnames(hierarchy) <- c("Sample","Group")

groups <- as.character(unique(hierarchy[,2]))

countable.groups <- c()
for(g in groups){
group.samples <- as.character(hierarchy[hierarchy$Group == g,1])
countable.group <- countable[,group.samples]
if(incidence == FALSE){countable.group.agg <- rowMeans(countable.group)}
if(incidence == TRUE){countable.group.agg <- to.incidence(countable.group)}
countable.groups <- cbind(countable.groups,countable.group.agg)
}

colnames(countable.groups) <- groups
if(relative==TRUE){countable.groups <- tss(countable.groups)}

return(countable.groups)
}
