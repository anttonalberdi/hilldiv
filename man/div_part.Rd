% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/div_part.R
\name{div_part}
\alias{div_part}
\title{Multi-level diversity partitioning}
\usage{
div_part(countable, qvalue, tree, hierarchy)
}
\arguments{
\item{countable}{A count table (matrix/data.frame) indicating the absolute or relative OTU/ASV abundances of multiple samples. Columns must refer to samples and rows to OTUs/ASVs.}

\item{qvalue}{A positive number, usually between 0 and 5, but most commonly 0, 1 or 2. It can be an integer or contain decimals.}

\item{tree}{A phylogenetic tree of class 'phylo'. The tip labels must match the row names in the OTU table. Use the function match_data() if the OTU names do not match.}

\item{hierarchy}{A matrix indicating the relation between samples (first column) and parent group(s).}
}
\value{
A list object containing details of hierarchical diversity partitioning.
}
\description{
Multi-level diversity partitioning following the multiplicative definition based on Hill numbers. Hierarchical levels are defined from L1 (minimum, sample) to Ln (maximum, whole system), and as many intermediate levels as wanted can be defined in between. The hierarchical structure of the system is defined with the hierarchy table. If no hierarchy table is inputed, the function yields a simple two-level partitioning between alpha (L1), beta and gamma (L2).
}
\details{
Multi-level diversity partitioning (based on Hill numbers)
}
\examples{
data(bat.diet.otutable)
data(bat.diet.tree)
data(bat.diet.hierarchy)
#Two level examples (L1=sample (alpha diversity), L2=whole system (gamma diversity))
div_part(bat.diet.otutable,qvalue=1)
div_part(bat.diet.otutable,qvalue=0,tree=bat.diet.tree)
#Three-level example (L1=sample, L2=species, L3=whole system)
div_part(bat.diet.otutable,qvalue=0,hierarchy=bat.diet.hierarchy)
}
\references{
Alberdi, A., Gilbert, M.T.P. (2019). A guide to the application of Hill numbers to DNA-based diversity analyses. Molecular Ecology Resources, 19, 804-817.\cr\cr
Chao, A., Chiu, C.H., & Hsieh, T. C. (2012). Proposing a resolution to debates on diversity partitioning. Ecology, 93, 2037-2051.\cr\cr
Jost, L. (2007). Partitioning diversity into independent alpha and beta components. Ecology, 88, 2427–2439.
}
\seealso{
\code{\link{div_part}}, \code{\link{gamma_div}}, \code{\link{match_data}}
}
\author{
Antton Alberdi, \email{anttonalberdi@gmail.com}
}
\keyword{alpha}
\keyword{beta}
\keyword{gamma}
\keyword{hill}
