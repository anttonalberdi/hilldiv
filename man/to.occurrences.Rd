% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/to.occurrences.R
\name{to.occurrences}
\alias{to.occurrences}
\title{Transform abundance vector or matrix into occurrences}
\usage{
to.occurrences(abund)
}
\arguments{
\item{abund}{A vector or a matrix/data.frame indicating the absolute or relative abundances of multiple samples. Columns must refer to samples and rows to OTUs/ASVs.}
}
\value{
A vector (one sample) or matrix (multiple samples) of occurrence data.
}
\description{
Transform an absolute or relative abundance vector (one sample) or matrix (multipla samples) into an occurrence vector or matrix.
}
\details{
To occurrences
}
\examples{
data(bat.diet.otutable)
to.occurrences(bat.diet.otutable)
to.occurrences(bat.diet.otutable[,1])
}
\references{
Alberdi, A., Gilbert, M.T.P. (2019). A guide to the application of Hill numbers to DNA-based diversity analyses. Molecular Ecology Resources, 19, 804-817.
}
\seealso{
\code{\link{to.incidence}}
}
\author{
Antton Alberdi, \email{anttonalberdi@gmail.com}
}
\keyword{ASV}
\keyword{OTU}
\keyword{abundance}
\keyword{hill}
\keyword{occurrence}
\keyword{table}
