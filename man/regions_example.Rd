\name{regions_example}
\alias{regions_example}
\docType{data}

\title{Example of regions}

\description{
  Example of regions used to illustrate usage of \code{"\link{IWTomicsData}"} object constructor. 
  It contains four different region datasets: "elem1", "elem2", "elem3" and "control".
}

\usage{data(regions_example)}

\value{
  \code{"\link[GenomicRanges]{GRangesList}"} object.
}

\format{
  A \code{"\link[GenomicRanges]{GRangesList}"} object with genomic regions of each region dataset: 
  "elem1", "elem2", "elem3" and "control". 
}

\examples{
data(regions_example)
}
\keyword{datasets}