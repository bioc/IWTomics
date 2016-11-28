\name{plotSummary}
\alias{plotSummary}

\title{Summary plots of Interval-Wise Testing result}

\description{
  The function creates a graphical summary of the Interval-Wise Testing results. In particular, 
  it draws a heatmap of the adjusted p-value curves at the chosen scale thresholds, grouped by 
  the region datasets tested (test) or by the features tested (feature). Different raws in the heatmap correspond 
  to different features (when grouping by locations) or to different tests (when grouping by feature). 
  Color intensity is proportional to -log(p-value), i.e. darker colors correspond to lower 
  p-values. Red means that the test statistics is higher in the first dataset tested than in the 
  second one (or is positive in one sample test), while blue means that the test statistics is 
  lower in the first datset tested than in the second one (or is negative in one sample test).
}

\usage{
plotSummary(regionsFeatures, alpha=0.05, only_significant=FALSE, 
  scale_threshold=NULL, nlevel=100, groupby='test', 
  test=1:nTests(regionsFeatures), gaps_tests=NULL,
  id_features_subset=idFeaturesTest(regionsFeatures), gaps_features=NULL,
  ask=TRUE, filenames=NA, align_lab=NA, cellwidth=NA, cellheight=NA,
  xlab='Windows', ylab=ifelse(groupby=='test','Features','Tests'), ...)
}

\arguments{
  \item{regionsFeatures}{\code{"IWTomicsData"} object with \code{test} slot.}
  \item{alpha}{level of the hypothesis test used to select and display significant 
    results. Default \code{alpha=0.05}.}
  \item{only_significant}{if \code{TRUE}, only the significant tests (p-value<=alpha in some 
    position) will be plotted.}
  \item{scale_threshold}{threshold on the test scale (maximum interval length for 
    the p-value adjustment) to be used in the adjusted p-value plot. Can be either 
    a scalar (the same length for all features) or a vector (a length for each feature) 
    or a list of vectors (a vector for each test). If \code{NULL}, the maximum possible 
    interval length is used.}
  \item{nlevel}{number of desired color levels for the adjusted p-value heatmap.}
  \item{groupby}{how tests should be grouped for the summary plot. Possible types are:
    \itemize{
      \item \code{"test"} to draw a summary plot for all the tests about the same region 
        dataset(s) comparison,
      \item \code{"feature"} to draw a summary plot for all the tests about the same feature.
    }
  }
  \item{test}{vector of indices of the tests to be plotted.}
  \item{gaps_tests}{vector of \code{test} indices that show where to put gaps in the heatmap 
    between tests to be plotted. Only used if \code{groupby} is \code{"feature"}.}
  \item{id_features_subset}{vector with the identifiers of the features to be plotted.}
  \item{gaps_features}{vector of id_features_subset indices that show where to put gaps 
    in the heatmap between features. Only used if \code{groupby} is \code{"test"}.}
  \item{ask}{if \code{TRUE} (default) the user is prompted before a new plot is drawn.}
  \item{filenames}{file paths where to save the plots (one for each group of tests as defined 
    by \code{groupby}). Filetypes are decided by the extension in the paths. Currently the 
    following formats are supported: png, pdf, tiff, bmp, jpeg.}
  \item{align_lab}{a title for the alignment point. Ignored if region alignment is \code{"scale"}.}
  \item{cellwidth}{individual cell width in points.}
  \item{cellheight}{individual cell height in points.}
  \item{xlab}{a title for the x axis.}
  \item{ylab}{a title for the y axis.}
  \item{\dots}{additional plot parameters.}
}

\value{No value returned. The function produces a graphical output.}

\note{
  This function uses a modified version of \pkg{pheatmap} package 
  (\url{https://cran.r-project.org/package=pheatmap}).
}

\author{
  Marzia A Cremona, Alessia Pini, Francesca Chiaromonte, Simone Vantini
}

\seealso{
  \code{\link{IWTomicsData}} for \code{"IWTomicsData"} class, constructors, accessors and methods; 
  \code{\link{IWTomicsTest}} function to perform the Interval-wise Testing;
  \code{\link{plotTest}} to draw detailed plots of the test results.
}

\examples{
## -------------------------------------------------------------------------------------------
## CURVE ALIGNMENT CENTER
## -------------------------------------------------------------------------------------------
data(regionsFeatures_center)

## One sample test for different regions and features
regionsFeatures_test=IWTomicsTest(regionsFeatures_center,
                                  id_region1=c('elem1','elem2','elem3','control'),mu=1)

## Summary plots grouped by feature
## x11(10,5)
plotSummary(regionsFeatures_test,groupby='feature')

## Set scale thresholds for the different features
## x11(10,5)
plotSummary(regionsFeatures_test,groupby='feature',
            scale_threshold=c(25,30))

## Add a title for the alignment point
## x11(10,5)
plotSummary(regionsFeatures_test,groupby='feature',
            scale.threshold=c(25,30),align_lab='Center')

## Plot only significant tests
## x11(10,5)
plotSummary(regionsFeatures_test,groupby='feature',
            scale.threshold=c(25,30),align_lab='Center',
            only_significant=TRUE)

## Summary plots grouped by test
## x11(10,5)
plotSummary(regionsFeatures_test,groupby='test')


## Two sample test for all possible region comparisons (mu=0), 
## and one sample test for 'elem3' (mu=1) for feature 'ftr1'
regionsFeatures_test=IWTomicsTest(regionsFeatures_center,
                       id_region1=c('elem1','elem2','elem3',
                                    'elem1','elem1','elem2','elem1'),
                       id_region2=c('control','control','control',
                                    'elem2','elem3','elem3',''),
                       id_features_subset='ftr1')

## Summary plots grouped by feature
## x11(10,5)
plotSummary(regionsFeatures_test,groupby='feature')

## Put gaps between different types of test and add a title for the alignment point
## x11(10,5)
plotSummary(regionsFeatures_test,groupby='feature',
            gaps_tests=c(3,6),align_lab='Center')

## Plot only significant tests
## x11(10,5)
plotSummary(regionsFeatures_test,groupby='feature',
            gaps_tests=c(3,6),align_lab='Center',only_significant=TRUE)

## Plot only the first three tests
## x11(10,5)
plotSummary(regionsFeatures_test,groupby='feature',
            test=1:3,align_lab='Center')

## Change scale threshold for the first test
## x11(10,5)
plotSummary(regionsFeatures_test,groupby='feature',
            test=1:3,align_lab='Center',scale_threshold=list(t1=8,t2=50,t3=50))
}