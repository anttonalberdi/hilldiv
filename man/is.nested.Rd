% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is.nested.R
\name{is.nested}
\alias{is.nested}
\title{Check if hierachy is nested}
\usage{
is.nested(hierarchy)
}
\arguments{
\item{hierarchy}{A matrix indicating the relation between samples (first column) and parent groups.}
}
\value{
A logical value (TRUE/FALSE).
}
\description{
Multi-level diversity partitioning requires the groups at different hierarchical levels to be nested. i.e. two samples that belong to a common parent group cannot have different grandparent groups. The best example of nested hierarchy is taxonomy: e.g. two species that belong to the same genus cannot belong to different families. This function checks whether the groups specified in a hierarchy table have a nested structure.
}
\details{
Check if hierachy is nested
}
\examples{
data(bat.diet.hierarchy)
is.nested(bat.diet.hierarchy)
}
\author{
Antton Alberdi, \email{anttonalberdi@gmail.com}
}
\keyword{hierarchy}
\keyword{levels}
\keyword{nestedness}
\keyword{partitioning}
