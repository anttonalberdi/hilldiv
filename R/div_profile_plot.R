#' Diversity profile plot
#' @title Diversity profile plot
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords alpha gamma beta hill
#' @description Plot diversity profiles from objects generated with the function div_profile().
#' @param profile A div_profile() object or a vector/matrix containg diversity profile(s), with columns indicating samples/groups and rows indicating orders of diversity (q-values).
#' @param colour A vector of RGB colours e.g. c("#34k235","#99cc00"). The number of vector items, must equal the number of samples or groups that are intended to plot.
#' @param log Whether to transform Hill numbers to logarithmic scale (TRUE) or not (FALSE). This is useful when there are large differences between q values (e.g. sharp drop from q=0 to q=1), which might complicate visualization. Default: log=FALSE
#' @param legend Whether to display the legend (TRUE) or not (FALSE) in diversity profiles containing multiple samples/groups. Default TRUE in multi-sample charts.
#' @return A diversity profile plot.
#' @seealso \code{\link{div_profile}}, \code{\link{hill_div}}
#' @import ggplot2
#' @importFrom data.table melt
#' @examples
#' data(bat.diet.otutable)
#' data(bat.diet.hierarchy)
#' #One sample example
#' bat.diet.sample <- bat.diet.otutable[,1]
#' profile.onesample <- div_profile(count=bat.diet.sample,qvalues=seq(from = 0, to = 5, by = (0.1)))
#' div_profile_plot(profile.onesample)
#' #Multiple samples
#' profile.multiplesamples <- div_profile(bat.diet.otutable)
#' div_profile_plot(profile.multiplesamples)
#' #Multiple groups (gamma diversity)
#' profile.multiplegroups <- div_profile(bat.diet.otutable,hierarchy=bat.diet.hierarchy,level="gamma")
#' div_profile_plot(profile.multiplegroups)
#' @references
#' Alberdi, A., Gilbert, M.T.P. (2019). A guide to the application of Hill numbers to DNA-based diversity analyses. Molecular Ecology Resources, 19, 804-817.\cr\cr
#' Chao, A., Chiu, C.H., & Jost, L. (2014). Unifying species diversity, phylogenetic diversity, functional diversity, and related similarity and differentiation measures through hill numbers. Annual Review of Ecology Evolution and Systematics, 45, 297-324.
#' @export

div_profile_plot <- function(profile,colour,log,legend){

#Data check
if(missing(profile)) stop("The diversity profile vector/table is missing")
if(is.null(dim(profile)) == TRUE){
  inputtype="onesample"
}else{
  inputtype="multiplesamples"
  }
if(missing(log)){log=FALSE}
if(missing(legend)){legend=TRUE}

#Declare colours
if(missing(colour) || (length(colour) != ncol(profile))){
colour <- colorRampPalette(brewer.pal(12, "Paired"))(ncol(profile))
}

#Single sample/group
if(inputtype == "onesample"){
  if(length(colour != 1)){colour="#99cc00"}
  profile.melted <- as.data.frame(cbind(names(profile),profile))
  colnames(profile.melted) <- c("Order","Profile")
  profile.melted[,1] <- as.numeric(as.character(profile.melted[,1]))
  profile.melted[,2] <- as.numeric(as.character(profile.melted[,2]))
  if(log == TRUE){profile.melted[,2] <- log(profile.melted[,2])}

  plot <- ggplot(profile.melted , aes_(x=~Order, y=~Profile, col=colour)) +
         geom_line() +
         xlab("Order of diversity") +
         ylab(if(log == TRUE){"Effective number of OTUs (log-transformed)" }else{"Effective number of OTUs"}) +
         theme_minimal() +
         theme(legend.position = "none")
  print(plot)
}

#Multiple samples/groups
if(inputtype == "multiplesamples"){
  profile.melted <- as.data.frame(melt(profile))
  colnames(profile.melted) <- c("Order","Sample","Value")
  profile.melted[,1] <- as.numeric(as.character(profile.melted[,1]))
  profile.melted[,3] <- as.numeric(as.character(profile.melted[,3]))
  if(log == TRUE){profile.melted[,3] <- log(profile.melted[,3])}

  #Plot
  plot <- ggplot(profile.melted , aes_(x=~Order, y=~Value, group=~Sample, colour=~Sample)) +
  geom_line() +
  xlab("Order of diversity") +
  ylab(if(log == TRUE){"Effective number of OTUs (log-transformed)" }else{"Effective number of OTUs"}) +
  scale_colour_manual(values = colour) +
  theme_minimal()
  if(legend == TRUE){
  plot <- plot + theme(legend.position = "right")}
  if(legend == FALSE){
  plot <- plot + theme(legend.position = "none")}
  print(plot)
}

}
