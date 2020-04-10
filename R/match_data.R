#' Match data
#' @title Match data
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords OTU ASV tree names
#' @description Filter count tables and OTU/ASV phylogenetic trees to match OTUs/ASVs present in both data files..
#' @param countable A count table (matrix/data.frame) indicating the absolute or relative OTU/ASV abundances of multiple samples. Columns must refer to samples and rows to OTUs/ASVs.
#' @param tree An ultrametic tree of class 'phylo'.
#' @param output Whether to output a filtered count table ('countable') or a filtered OTU tree ('tree'). Default is empty, which only yields a message.
#' @param silent Whether to stop printing text on screen. Default=FALSE.
#' @seealso \code{\link{hill_div}}, \code{\link{index_div}}
#' @examples
#' data(bat.diet.otutable)
#' data(bat.diet.tree)
#' match_data(bat.diet.otutable,bat.diet.tree,output="countable")
#' match_data(bat.diet.otutable,bat.diet.tree,output="tree")
#' @export

match_data <- function(countable,tree,output,silent){
if(missing(silent)){silent=FALSE}

#Obtain data
  #OTUs in OTU table
  countable.otus <- rownames(countable)
  #Samples in OTU table
  countable.samples <- colnames(countable)
  #OTUs in tree
  tree.otus <- tree$tip.label

if(missing(output)){
  if((length(setdiff(tree.otus,countable.otus)) == 0) & (length(setdiff(countable.otus,tree.otus)) == 0)){message("OTUs in the OTU table and OTU tree match perfectly.")}
  if((length(setdiff(tree.otus,countable.otus)) > 0) & (length(setdiff(countable.otus,tree.otus)) == 0)){message("The OTU tree contains OTUs absent in the OTU table. Filter the OTU tree")}
  if((length(setdiff(tree.otus,countable.otus)) == 0) & (length(setdiff(countable.otus,tree.otus)) > 0)){message("The OTU table contains OTUs absent in the OTU tree. Filter the OTU table")}
  if((length(setdiff(tree.otus,countable.otus)) > 0) & (length(setdiff(countable.otus,tree.otus)) > 0)){message("The OTU table contains OTUs absent in the OTU tree and the OTU tree contains OTUs absent in the OTU table. Filter both files")}
  output="NA"
}

#Output OTU table
if(output == "countable"){
  if(length(setdiff(tree.otus,countable.otus)) > 0){
    if(silent == FALSE){
    message("All OTUs/ASVs present in the count table are present in the tree, but the tree contains OTUs/ASVs absent in the count table. Remember to filter the tree.")
    }
    return(countable)
  }
  if(length(setdiff(countable.otus,tree.otus)) > 0){
    OTUs.to.drop <- setdiff(countable.otus,tree.otus)
    countable.filt <- countable[!(row.names(countable) %in% OTUs.to.drop), ]
    OTUs.to.drop.string <- paste(OTUs.to.drop,collapse=", ")
    if(silent == FALSE){
    message("The following OTUs/ASVs were removed from the count table for being absent in the tree: ",OTUs.to.drop.string)
    }
    return(countable.filt)
  }
  if((length(setdiff(tree.otus,countable.otus)) == 0) & (length(setdiff(countable.otus,tree.otus)) == 0)){
    if(silent == FALSE){
    message("OTUs/ASVs in the count table and tree match perfectly. No new count table was created.")
    }
    return(countable)
  }

}

#Output tree
if(output == "tree"){
  if(length(setdiff(tree.otus,countable.otus)) > 0){
    OTUs.to.drop <- setdiff(tree.otus,countable.otus)
    tree.filt <- drop.tip(tree,OTUs.to.drop)
    OTUs.to.drop.string <- paste(OTUs.to.drop,collapse=", ")
    if(silent == FALSE){
    message("The following OTUs/ASVs were removed from the tree for being absent in the count table: ",OTUs.to.drop.string)
    }
    return(tree.filt)
  }
  if(length(setdiff(countable.otus,tree.otus)) > 0){
    if(silent == FALSE){
    message("All OTUs/ASVs present in the tree are present in the count table, but the count table contains OTUs/ASVs absent in the tree. Remember to filter the count table.")
    }
    return(tree)
  }
  if((length(setdiff(tree.otus,countable.otus)) == 0) & (length(setdiff(countable.otus,tree.otus)) == 0)){
    if(silent == FALSE){
    message("OTUs/ASVs in the tree and count table match perfectly. No new tree was created.")
    }
    return(tree)
  }
}

}
