\name{plotTest}
\alias{plotTest}

\title{Detailed plots of Interval-Wise Testing result}

\description{
  The function \code{plotTest} creates graphical representations of the Interval-Wise 
  Testing results present in the \code{test} slot of a \code{"\link{IWTomicsData}"} object. 
  In particular, it draws a heatmap of the adjusted p-value functions at each possible 
  scale (with the different maximum interval lengths considered), a plot of the adjusted 
  p-value curve at the chosen scale threshold and a plot of the feature measurements in 
  the region datasets tested (aligned curves or pointwise quantile curves). 
}

\usage{
plotTest(regionsFeatures, alpha=0.05, scale_threshold=NULL, nlevel=100, type="boxplot", 
  N_regions=pmin(lengthRegions(regionsFeatures),10), 
  probs=c(0.25,0.5,0.75), average=TRUE, size=TRUE, 
  id_features_subset=idFeatures(regionsFeatures), 
  col=1+seq_len(nRegions(regionsFeatures)), 
  ask=TRUE, xlab="Windows", ylim=NULL, ...)
}

\arguments{
  \item{regionsFeatures}{\code{"IWTomicsData"} object with \code{test} slot.}
  \item{alpha}{level of the hypothesis test used to select and display significant 
    results. Default \code{alpha=0.05}.}
  \item{scale_threshold}{threshold on the test scale (maximum interval length for 
    the p-value adjustment) to be used in the adjusted p-value plot. Can be either 
    a scalar (the same length for all features) or a vector (a length for each feature) 
    or a list of vectors (a vector for each test). If \code{NULL}, the maximum possible 
    interval length is used.}
  \item{nlevel}{number of desired color levels for the adjusted p-value heatmap.}
  \item{type}{type of plot to be drawn in addition to the adjusted p-value heatmap 
    and plot. Possible types are:
    \itemize{
      \item \code{"boxplot"} for the pointwise quantiles curves (default), 
      \item \code{"curves"} for the (aligned) curves.
    }
  }
  \item{N_regions}{number of regions to be randomly chosen (for each region dataset) in 
    \code{"curves"} type plot. Default plots a maximum of 10 curves.}
  \item{probs}{probabilities corresponding to the quantiles to be drawn in \code{"boxplot"} 
    type plot. Default plots quartile curves.}
  \item{average}{if \code{TRUE}, plot the mean curves.}
  \item{size}{if \code{TRUE}, plot the sample size in each position.}
  \item{id_features_subset}{vector with the identifiers of the features to be plotted.}
  \item{col}{vector of plotting colors for the different region datasets.}
  \item{ask}{if \code{TRUE} (default) the user is prompted before a new plot is drawn.}
  \item{xlab}{a title for the x axis.}
  \item{ylim}{the y limits of the feature measurement plot.}
  \item{\dots}{additional plot parameters.}
}

\value{No value returned. The function produces a graphical output.}

\references{
  A Pini and S Vantini (2017). Interval-Wise Testing for functional data. 
  \emph{Journal of Nonparametric Statistics}.
}

\author{
  Marzia A Cremona, Alessia Pini, Francesca Chiaromonte, Simone Vantini
}

\seealso{
  \code{\link{IWTomicsData}} for \code{"IWTomicsData"} class, constructors, accessors and methods; 
  \code{\link{IWTomicsTest}} function to perform the Interval-wise Testing;
  \code{\link{plotSummary}} to draw a summary plot of the test results.
}


\examples{
data(regionsFeatures_center)

## One sample test for 'control' regions and feature 'ftr1'
regionsFeatures_test=IWTomicsTest(regionsFeatures_center,
                                  id_region1='control',id_features_subset='ftr1')

## Plotting the results of the one sample test
plotTest(regionsFeatures_test,col=5)

## Set the scale threshold to 25
plotTest(regionsFeatures_test,col=5,scale_threshold=25)

## Plot curves instead of boxplots, do not report sample size in each position
plotTest(regionsFeatures_test,col=5,type='curves',size=FALSE)


## Two sample test for 'elem1', 'elem2' and 'elem3' vs 'control' regions
regionsFeatures_test=IWTomicsTest(regionsFeatures_center,
                                  id_region1=c('elem1','elem2','elem3'),
                                  id_region2=c('control','control','control'))

## Plotting the results of the two sample test
plotTest(regionsFeatures_test)

## Set scale thresholds for the different features
plotTest(regionsFeatures_test,scale_threshold=c(25,30))

## Plot only results regarding feature 'ftr2', setting scale thresholds for each test
plotTest(regionsFeatures_test,id_features_subset='ftr2',
         scale_threshold=list(test1=10,test2=20,test3=30))

## Report also 0.05 and 0.95 quantiles in the boxplots
plotTest(regionsFeatures_test,id_features_subset='ftr2',
         scale_threshold=list(test1=10,test2=20,test3=30),
         probs=c(0.05,0.25,0.5,0.75,0.9))
}