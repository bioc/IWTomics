\name{IWTomicsTest}
\alias{IWTomicsTest}

\title{Interval-Wise Testing}

\description{
  The function implements the Interval-Wise Testing for omics data (both one sample 
  and two sample tests), an extended version of the Interval-Wise Testing for functional 
  data presented in Pini and Vantini (2017). This inferential procedure tests for 
  differences in feature measurements between two region datasets (two sample test) 
  or between a region dataset and a reference curve (one sample test).
}

\usage{
IWTomicsTest(regionsFeatures, 
  id_region1=idRegions(regionsFeatures)[1], id_region2=NULL, 
  id_features_subset=idFeatures(regionsFeatures), mu=0, 
  statistics="mean", probs=0.5, max_scale=NULL, paired=FALSE, B=1000)
}

\arguments{
  \item{regionsFeatures}{\code{"IWTomicsData"} object.}
  \item{id_region1}{identifier(s) of the region dataset(s) to be tested. If a vector 
    is provided, a test will be performed for each element of the vector.}
  \item{id_region2}{identifier(s) of the region dataset(s) to be tested for two sample 
    test. If \code{NULL} or empty string, one sample test is performed.}
  \item{id_features_subset}{vector with the identifiers of the features to be tested.}
  \item{mu}{the reference curve (center of symmetry) under the null hypothesis 
    in one sample test, or the difference between the two populations in two 
    sample test. Can be either a constant (the same constant curve for all features), 
    a vector of constants (a constant curve for each feature), or a list of vectors 
    containing its measurements for each feature (on the same grid as the features). 
    Default \code{mu=0}.}
  \item{statistics}{test statistics to be used in the test. Possible test statistics 
    are \code{"mean"}, \code{"median"}, \code{"variance"} and \code{"quantile"}.}
  \item{probs}{probabilities corresponding to the quantiles in test statistics 
    \code{"quantile"}. If multiple quantiles are selected, the test statistics is the sum 
    of the statistics on the different quantiles.}
  \item{max_scale}{the maximum interval length to be used for the p-value adjustment, 
    i.e. the maximum number of consecutive windows to be employed (can range from 1 to 
    the maximum number of consecutive measurements present for the feature tested). 
    Can be either a scalar (the same length for all features) or a vector (a length 
    for each feature) or a list of vectors (a vector for each test). If \code{NULL}, 
    the maximum possible interval length is used.}
  \item{paired}{if \code{TRUE}, a paired (two sample) test is performed.}
  \item{B}{number of iterations of the MC algorithm to evaluate the p-values of the 
    permutation tests. If B is greater than the number of possible permutations, 
    exact p-values are computed.}
}

\value{
  The function \code{IWTomicsTest} returns an object of class \code{"\link{IWTomicsData}"}
  containing the region datasets and feature datasets tested, and with the test 
  input and result in the optional slot \code{test}. 
}

\note{
  In this implementation of Interval-Wise Testing, the smallest scale considered 
  corresponds to the measurement resolution, i.e. the univariate (unadjusted) tests 
  are done on each measurement window. To change the smallest scale considered the 
  method \code{\link[IWTomics]{smooth}} can be employed. 

  If the region alignment is \code{"scale"}, the function 
  \code{\link[IWTomics]{smooth}} must be used before applying the 
  Interval-Wise Testing in order to measure the features over the grid in all the regions. 
}

\references{
  A Pini and S Vantini (2017). Interval-Wise Testing for functional data. 
  \emph{Journal of Nonparametric Statistics}.
}

\author{
  Marzia A Cremona, Alessia Pini, Francesca Chiaromonte, Simone Vantini
}

\seealso{
  \code{\link{IWTomicsData}} for \code{"IWTomicsData"} class, constructors, accessors and methods; 
  \code{\link[IWTomics]{plot}} method to plot \code{"IWTomicsData"} objects;
  \code{\link[IWTomics]{smooth}} to smooth curves before testing;
  \code{\link{plotTest}} to plot the test results;
  \code{\link{plotSummary}} to draw a summary plot of the test results.
}

\examples{
## -------------------------------------------------------------------------------------------
## -------------------------------------------------------------------------------------------
## EXAMPLE ON REAL DATA
## -------------------------------------------------------------------------------------------
## -------------------------------------------------------------------------------------------
## ETn Recombination hotspots data
## Two region datasets:
## ETns fixed 64-kb flanking regions and 64-kb control regions in mouse
## One feature measured in 1-kb windows:
## Recombination hotspots content
## ?ETn_example for details on the dataset
data(ETn_example)
ETn_example

## Two sample test to compare recombination hotspots 
## in ETn regions vs control regions
ETn_test=IWTomicsTest(ETn_example,
                     id_region1='ETn_fixed',id_region2='Control')
## Adjusted p-value
adjusted_pval(ETn_test)

## Adjusted p-value lowering the scale of the test
adjusted_pval(ETn_test,scale_threshold=10)


## -------------------------------------------------------------------------------------------
## -------------------------------------------------------------------------------------------
## EXAMPLE ON SIMULATED DATA
## -------------------------------------------------------------------------------------------
## -------------------------------------------------------------------------------------------

## -------------------------------------------------------------------------------------------
## -------------------------------------------------------------------------------------------
## CURVE ALIGNMENT CENTER
## -------------------------------------------------------------------------------------------
## -------------------------------------------------------------------------------------------
data(regionsFeatures_center)

## One sample test for 'control' regions and feature 'ftr1'
regionsFeatures_test=IWTomicsTest(regionsFeatures_center,
                                  id_region1='control',id_features_subset='ftr1')
adjusted_pval(regionsFeatures_test)

## Plotting the results of the one sample test
plotTest(regionsFeatures_test,col=5)


## Two sample test for 'elem1', 'elem2' and 'elem3' vs 'control' regions and feature 'ftr1'
regionsFeatures_test=IWTomicsTest(regionsFeatures_center,
                                  id_region1=c('elem1','elem2','elem3'),
                                  id_region2=c('control','control','control'),
                                  id_features_subset='ftr1')
adjusted_pval(regionsFeatures_test)

## Plotting the results of the two sample test
plotTest(regionsFeatures_test)

## Summary plot of the two sample test
## x11(10,5)
plotSummary(regionsFeatures_test,groupby='feature',align_lab='Center')

#############################################################################################
## Not run: 
#############################################################################################
## Using 'quantile' test statistics with multiple quantiles
regionsFeatures_test=IWTomicsTest(regionsFeatures_center,
                                  id_region1=c('elem1','elem2','elem3'),
                                  id_region2=c('control','control','control'),
                                  id_features_subset='ftr1',
                                  statistics='quantile',probs=c(0.25,0.75))
adjusted_pval(regionsFeatures_test)

## Plotting the results of the two sample test
plotTest(regionsFeatures_test)

## Summary plot of the two sample test
## x11(10,5)
plotSummary(regionsFeatures_test,groupby='feature',align_lab='Center')
#############################################################################################
## End(Not run) 
#############################################################################################


## -------------------------------------------------------------------------------------------
## -------------------------------------------------------------------------------------------
## CURVE ALIGNMENT SCALE
## -------------------------------------------------------------------------------------------
## -------------------------------------------------------------------------------------------
data(regionsFeatures_scale)

## Smooth the curves to have measurements on the same grid
regionsFeatures_smooth=smooth(regionsFeatures_scale,type='locpoly',scale_grid=30)

## Two sample test for 'elem1', 'elem2' and 'elem3' vs 'control' regions and feature 'ftr2'
regionsFeatures_test=IWTomicsTest(regionsFeatures_smooth,
                                  id_region1=c('elem1','elem2','elem3'),
                                  id_region2=c('control','control','control'),
                                  id_features_subset='ftr2')
adjusted_pval(regionsFeatures_test)

## Plotting the results of the two sample test
plotTest(regionsFeatures_test)

## Summary plot of the two sample test
## x11(10,5)
plotSummary(regionsFeatures_test,groupby='feature')
}