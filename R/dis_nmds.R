#' Dissimilarity NMDS plot
#' @title Dissimilarity NMDS plot
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords hill numbers diversity partitioning pairwise ordination plot
#' @description Visualisation of pairwise dissimilarities
#' @param distance Matrix of pairwise dissimilarities, usually one of the matrices listed in the output object of the pair_dis() function.
#' @param hierarchy The first column containing the sample names while the second containing the groups names. If provided, dots are coloured according to groups and group centroids can also be visualised.
#' @param plot Whether to plot a NMDS or a Shepard plot. Default: "NMDS".
#' @param colour The number of vector items (colours, e.g. '#34k235'), must be of length one if no hierarchy table is added, or must equal the number of groups if the hierarchy table is provided.
#' @param centroids Whether to link sample dots with group centroids or not. A hierarchy table is necessary to draw centroids. Default: FALSE
#' @param labels Whether to print sample or group labels or both. A hierarchy table is necessary to plot grpup names. Default: "none".
#' @param legend Whether to print the legend or not. Default: TRUE.
#' @param runs Number of iterations for the NMDS function. Default: 100.
#' @return An NMDS or Shepard plot.
#' @seealso \code{\link{pair_dis}}, \code{\link{beta_dis}}
#' @import ggplot2
#' @import ggrepel
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @importFrom stats as.dist
#' @importFrom vegan metaMDS stressplot
#' @examples
#' data(bat.diet.otutable)
#' data(bat.diet.tree)
#' data(bat.diet.hierarchy)
#' pairdisres <- pair_dis(bat.diet.otutable,qvalue=0,hierarchy=bat.diet.hierarchy)
#' dis_nmds(pairdisres$L1_CqN)
#' dis_nmds(pairdisres$L1_CqN,hierarchy=bat.diet.hierarchy, centroids=TRUE)
#' dis_nmds(pairdisres$L1_CqN,hierarchy=bat.diet.hierarchy, centroids=TRUE, labels="group")
#' @references
#' Alberdi, A., Gilbert, M.T.P. (2019). A guide to the application of Hill numbers to DNA-based diversity analyses. Molecular Ecology Resources, 19, 804-817.\cr\cr
#' Chao, A., Chiu, C.H., & Hsieh, T. C. (2012). Proposing a resolution to debates on diversity partitioning. Ecology, 93, 2037-2051.\cr\cr
#' Jost, L. (2007). Partitioning diversity into independent alpha and beta components. Ecology, 88, 2427-2439.
#' @export

dis_nmds <- function(distance,hierarchy,colour,plot,centroids,labels,legend,runs){

if(missing(plot)){plot = "NMDS"}
if(missing(centroids)){centroids = FALSE}
if(missing(legend)){legend = TRUE}
if(missing(labels)){labels = "none"}
if((centroids == TRUE) & missing(hierarchy))stop("A hierarchy table is necessary to draw centroids.")
if((labels == "group") & missing(hierarchy))stop("A hierarchy table is necessary to print group names.")
if((labels == "both") & missing(hierarchy))stop("A hierarchy table is necessary to print group names.")
if(missing(runs)){runs = 100}


######
# Colours
######

if(!missing(hierarchy)){
	if(setequal(colnames(distance),hierarchy[,1]) == FALSE) stop("Distance matrix and hierarchy table do not match.")
	colnames(hierarchy) <- c("Sample","Group")
	if(missing(colour) || (length(colour) != length(unique(hierarchy[,2])))){
	getPalette <- colorRampPalette(brewer.pal(9, "Paired"))
	colour <- getPalette(length(unique(hierarchy[,2])))
	colourmap <- cbind(Group=as.character(unique(hierarchy[,2])),Colour=colour)
	hierarchy.colour <- merge(hierarchy,colourmap,by="Group")
	hierarchy.colour <- hierarchy.colour[,c(2,1,3)]
	hierarchy.colour[,3] <- as.character(hierarchy.colour[,3])
	}
}else{
  if(missing(colour)){
	colour <- "#999999"
  }
}

if(missing(hierarchy) & (length(colour) != 1))stop("The colour vector can only contain one value when not specifying a hierarchy table")

######
# NMDS plot
######

values.NMDS <- metaMDS(as.dist(distance), k = 2, trymax = runs)

if(missing(hierarchy)){
NMDS=data.frame(x=values.NMDS$point[,1],y=values.NMDS$point[,2])
NMDS$Sample <- rownames(NMDS)
}
if(!missing(hierarchy)){
NMDS=data.frame(x=values.NMDS$point[,1],y=values.NMDS$point[,2],Group=as.factor(hierarchy[,2]))
NMDS$Sample <- rownames(NMDS)

#Find centroids
NMDS.centroids=aggregate(NMDS[,c(1:2)],by=list(NMDS[,3]),FUN=mean)
colnames(NMDS.centroids) <- c("Group","x_cen","y_cen")
NMDS=merge(NMDS,NMDS.centroids,by="Group")
}

if(plot == "NMDS"){
  if(missing(hierarchy)){
    nmds.plot <- ggplot(NMDS, aes(x,y)) +
    geom_point(size=1,colour=colour) +
    theme(panel.background = element_rect(fill = 'white', colour = 'grey'))
    }
  if(!missing(hierarchy)){
    nmds.plot <- ggplot(NMDS, aes(x,y,colour=Group)) +
      geom_point(size=1) +
      scale_colour_manual(values = colour) +
      theme(panel.background = element_rect(fill = 'white', colour = 'grey'))
  }
  if(centroids == TRUE){
    nmds.plot <- nmds.plot +
    geom_point(aes(x=x_cen,y=y_cen),size=1,alpha=0) +
    geom_segment(aes(x=x_cen, y=y_cen, xend=x, yend=y), alpha=0.2)
      if(labels == "group" || labels == "both"){
        nmds.plot <- nmds.plot + geom_text_repel(aes(x=x_cen,y=y_cen,label=Group),size=3,vjust=0,force=0)
      }
  }

  if(labels == "sample" || labels == "both"){
  nmds.plot <- nmds.plot + geom_text_repel(aes(x=x,y=y,label=Sample),size=2,vjust=0,segment.size=0,box.padding=0.1,point.padding=0.1)
  }

	if(legend == FALSE){
	nmds.plot <- nmds.plot + theme(legend.position = "none")
	}

  return(nmds.plot)

}

######
# Shepard plot
######

if(plot == "Shepard"){
stressplot(values.NMDS)
}

}
