\name{features_example}
\alias{features_example}
\docType{data}

\title{Example of features}

\description{
  Example of features used to illustrate usage of \code{"\link{IWTomicsData}"} object constructor. 
  It contains measurements of two different genomic features, "ftr1" and "ftr2", corresponding to 
  four different region datasets: "elem1", "elem2", "elem3" and "control".
}

\usage{data(features_example)}

\value{
  list of two elements, one for each of the two features "ftr1" and "ftr2". 
}

\format{
  Each element of the list is a list of matrices with measurements corresponding to four different 
  region datasets: "elem1", "elem2", "elem3" and "control". Columns correspond to measurements in a region.
}

\examples{
data(features_example)
}
\keyword{datasets}
