#' Diversity test plotting
#' @title Diversity test plotting
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords hill numbers comparison chart
#' @description Plot of diversity comparison between groups of samples
#' @param divtest Object outputed by the div_test() function
#' @param chart Chart type, either 'box' for boxplot, 'jitter' for jitter plot or 'violin' for violin plot. chart="box"
#' @param colour The number of vector items (colours, e.g. '#34k235'), must equal the number of groups that are intended to plot.
#' @param posthoc If 'TRUE' pairwise p-values of the posthoc analyses will be ploted. It requires the div_test() object to contain posthoc results.
#' @param threshold Maximum p-value to show in pairwise posthoc results (usually 0.05, but could be any other number between 0 an 1). P-values above the threshold will not be showed.
#' @return Chart of (mean) diversities of contrasting groups with optional posthoc results.
#' @seealso \code{\link{div_test}}, \code{\link{hill_div}}, \code{\link{div_part}}
#' @import ggplot2
#' @import ggpubr
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @importFrom stats as.dist
#' @examples
#' data(bat.diet.otutable)
#' data(bat.diet.hierarchy)
#' divtestres <- div_test(bat.diet.otutable,qvalue=0,hierarchy=bat.diet.hierarchy)
#' div_test_plot(divtestres,chart="box")
#' div_test_plot(divtestres,chart="violin")
#' divtest.res.ph <- div_test(bat.diet.otutable,qvalue=0,hierarchy=bat.diet.hierarchy,posthoc=TRUE)
#' div_test_plot(divtest.res.ph,chart="jitter",posthoc=TRUE,threshold=0.5)
#' @export

div_test_plot <- function(divtest,chart,colour,posthoc,threshold){
if(missing(chart)){chart="box"}
if(missing(posthoc)){posthoc=FALSE}
if((names(divtest)[1] != "data") & (names(divtest)[2] != "normality.pvalue")) stop("The input object does not seem to be a div_test output.")

#Get data table
divtestdata <- divtest$data
divtestdata$Group <- as.factor(divtestdata$Group)
#Sort factor levels
divtestdata$Group <- factor(divtestdata$Group, levels = as.character(unique(divtestdata$Group)))

#Declare colours
if(missing(colour) || (length(colour) < divtest$groups)){
  if(divtest$groups == 2){
    colour=c("#96c3dc","#1b63a5")
  }else{
    getPalette <- colorRampPalette(brewer.pal(divtest$groups, "Paired"))
    colour <- getPalette(divtest$groups)
  }
}

if(posthoc == TRUE){
  if(is.na(names(divtest)[7])) stop("The input div_test object does not seem to contain pairwise posthoc data. Re-run div_test() using 'posthoc=TRUE' argument.")

  #Prepare pairwisetable from posthoc data
  if(names(divtest)[7] == "posthoc.method"){
  combinations <- matrix(gsub(" $","",gsub("^ ","",unlist(strsplit(as.character(rownames(divtest$posthoc)), "-", fixed = TRUE)))),ncol=2,byrow=TRUE)
  pvalue <- round(divtest$posthoc[,2],3)
  pairwisetable <- as.data.frame(cbind(combinations,pvalue))
  colnames(pairwisetable) <- c("group1","group2","p")
  }
  pairwisetable[,1] <- as.character(pairwisetable[,1])
  pairwisetable[,2] <- as.character(pairwisetable[,2])
  pairwisetable[,3] <- as.numeric(as.character(pairwisetable[,3]))

  #Filter pairwisetable
  if(!missing(threshold)){
  pairwisetable <- pairwisetable[which(pairwisetable$p < threshold),]
  }

  #Set y values
  sortedgroups <- unique(sort(c(pairwisetable$group1,pairwisetable$group2)))
  datamax <- round(max(divtest$data[which(divtest$data$Group %in% sortedgroups),3]))
  datamin <- round(min(divtest$data[which(divtest$data$Group %in% sortedgroups),3]))
  datarange <- datamax - datamin
  by <- datarange * 0.1
  min <- datamax
  max <- min + (by*nrow(pairwisetable))
  ypos <- seq(min,max,by)[-1]
  pairwisetable$ypos <- ypos
}

#Plot
if(chart == "box"){
plot <- ggboxplot(divtestdata, x='Group', y='Value', color='Group', fill='Group', x.text.angle = 45) +
      ylab("Effective number of OTUs") + xlab("Groups") +
      scale_colour_manual(values=scales::alpha(colour, 1)) +
      scale_fill_manual(values=scales::alpha(colour, 0.5))
  if(posthoc == TRUE){
    plot <- suppressWarnings(plot + stat_pvalue_manual(pairwisetable, label = "p", y.position = "ypos"))
  }
print(plot)
}

if(chart == "jitter"){
plot <- ggboxplot(divtestdata, x='Group', y='Value', color='Group', add = "jitter", width = 0, x.text.angle = 45) +
      ylab("Effective number of OTUs") + xlab("Groups") +
      scale_colour_manual(values=scales::alpha(colour, 0))
    if(posthoc == TRUE){
      plot <- suppressWarnings(plot + stat_pvalue_manual(pairwisetable, label = "p", y.position = "ypos"))
    }
  print(plot)
  }

if(chart == "violin"){
  plot <- ggviolin(divtestdata, x='Group', y='Value', color='Group', fill='Group', x.text.angle = 45) +
        ylab("Effective number of OTUs") + xlab("Groups") +
        scale_fill_manual(values=scales::alpha(colour, 0.1)) +
        scale_colour_manual(values=scales::alpha(colour, 1))
    if(posthoc == TRUE){
      plot <- suppressWarnings(plot + stat_pvalue_manual(pairwisetable, label = "p", y.position = "ypos"))
    }
  print(plot)
  }

}
