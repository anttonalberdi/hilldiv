#' Check if hierachy is nested
#' @title Check if hierachy is nested
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords hierarchy levels partitioning nestedness
#' @description Multi-level diversity partitioning requires the groups at different hierarchical levels to be nested. i.e. two samples that belong to a common parent group cannot have different grandparent groups. The best example of nested hierarchy is taxonomy: e.g. two species that belong to the same genus cannot belong to different families. This function checks whether the groups specified in a hierarchy table have a nested structure.
#' @param hierarchy A matrix indicating the relation between samples (first column) and parent groups.
#' @return A logical value (TRUE/FALSE).
#' @examples
#' data(bat.diet.hierarchy)
#' is.nested(bat.diet.hierarchy)
#' @export

is.nested <- function(hierarchy){
if(is.null(dim(hierarchy)) == TRUE) stop("The hierarchy object is not a two-dimensional table.")
leveln <- ncol(hierarchy)
logic.vector <- c()
if(leveln == 2){return(TRUE)}
if(leveln > 2){
for (i in c((leveln-1):2)){
levelonly <- length(unique(hierarchy[,i]))
levelparent <- nrow(unique(hierarchy[,c(i,i+1)]))
logic <- levelonly == levelparent
logic.vector <- c(logic.vector,logic)
}
return(all(logic.vector))
}
}
