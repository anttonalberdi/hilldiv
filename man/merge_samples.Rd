% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/merge_samples.R
\name{merge_samples}
\alias{merge_samples}
\title{Merge samples}
\usage{
merge_samples(countable,hierarchy,incidence)
}
\arguments{
\item{countable}{A count table (matrix/data.frame) indicating the absolute or relative OTU/ASV abundances of multiple samples. Columns must refer to samples and rows to OTUs/ASVs.}

\item{hierarchy}{A two-column matrix indicating the relation between samples (first column) and groups (second column).}

\item{relative}{Whether to output relative values or not. Default=TRUE.}

\item{incidence}{Whether to transform abundance into incidence data when merging. Default=FALSE.}
}
\value{
A count table
}
\description{
Combines samples into groups defined by the hierarchy table, with the possibility to convert abundances into incidence data.
}
\details{
Merge samples
}
\examples{
data(bat.diet.otutable)
data(bat.diet.hierarchy)
merge_samples(countable=bat.diet.otutable,hierarchy=bat.diet.hierarchy)
merge_samples(bat.diet.otutable,bat.diet.hierarchy)
merge_samples(countable=bat.diet.otutable,hierarchy=bat.diet.hierarchy, incidence=TRUE)
}
\references{
Alberdi, A., Gilbert, M.T.P. (2019). A guide to the application of Hill numbers to DNA-based diversity analyses. Molecular Ecology Resources, 19, 804-817.\cr\cr
}
\seealso{
\code{\link{to.incidence}}
}
\author{
Antton Alberdi, \email{anttonalberdi@gmail.com}
}
\keyword{ASV}
\keyword{OTU}
\keyword{combine}
\keyword{table}
