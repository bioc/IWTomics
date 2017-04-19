\name{IWTomics-package}
\alias{IWTomics-package}
\alias{IWTomics}
\docType{package}

\title{Interval-Wise Testing for Omics Data}

\description{
  Implementation of the Interval-Wise Testing for "Omics" data, an extended version 
  of the Interval-Wise Testing for functional data presented in Pini and Vantini (2017). 
  This inferential procedure tests for differences in "Omics" data between two groups 
  of genomic regions (or between a group of genomic regions and a reference curve), 
  and does not require fixing location and scale at the outset. 
}

\author{
  Marzia A Cremona, Alessia Pini, Francesca Chiaromonte, Simone Vantini
  
  Maintainer: Marzia A Cremona <mac78@psu.edu>
}
\references{
  A Pini and S Vantini (2017). Interval-Wise Testing for functional data. 
  \emph{Journal of Nonparametric Statistics}.
}
\keyword{package}

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

## -------------------------------------------------------------------------------------------
## PLOT DATA
## -------------------------------------------------------------------------------------------
## Plot the pointwise boxplot and averages of the curves in ETn and control regions
## (note that the box of the pointwise boxplot is zero, but the average is not)
plot(ETn_example)

## -------------------------------------------------------------------------------------------
## PERFORM THE TEST
## -------------------------------------------------------------------------------------------
## Two sample test to compare recombination hotspots 
## in ETn regions vs control regions
ETn_test=IWTomicsTest(ETn_example,
                     id_region1='ETn_fixed',id_region2='Control')
## Adjusted p-value
adjusted_pval(ETn_test)

## Plotting the test results
plotTest(ETn_test)

## Adjusted p-value lowering the scale of the test
adjusted_pval(ETn_test,scale_threshold=10)

## Plotting the test results, lowering the scale of the test
plotTest(ETn_test,scale_threshold=10)

## Summary plot of the two sample test
## x11(12,2)
plotSummary(ETn_test,groupby='feature',scale_threshold=10,
            align_lab='Integration site')


## -------------------------------------------------------------------------------------------
## -------------------------------------------------------------------------------------------
## EXAMPLE ON SIMULATED DATA
## -------------------------------------------------------------------------------------------
## -------------------------------------------------------------------------------------------
examples_path <- system.file("extdata",package="IWTomics")
## Four region datasets: 
## three different types of elements and one control
datasets=read.table(file.path(examples_path,"datasets.txt"),
                    sep="\t",header=TRUE,stringsAsFactors=FALSE)
datasets
## Two different features measured in all four types regions
features_datasetsTable=read.table(file.path(examples_path,"features_datasetsTable.txt"),
                                  sep="\t",header=TRUE,stringsAsFactors=FALSE)
features_datasetsTable

## -------------------------------------------------------------------------------------------
## GET DATA AND PLOT
## -------------------------------------------------------------------------------------------
## Get genomic regions for the four region datasets, 
## and the two features from Table files for each region dataset
regionsFeatures=IWTomicsData(datasets$regionFile,features_datasetsTable[,3:6],'center',
                            datasets$id,datasets$name,
                            features_datasetsTable$id,features_datasetsTable$name,
                            path=file.path(examples_path,'files'))

## Plot the pointwise boxplot of the curves in the different region datasets
plot(regionsFeatures)

## -------------------------------------------------------------------------------------------
## PERFORM THE TEST
## -------------------------------------------------------------------------------------------
## Two sample test for the two features in the comparisons
## 'elem1' vs 'control', 'elem2' vs 'control' 
## and 'elem3' vs 'control'
regionsFeatures_test=IWTomicsTest(regionsFeatures,id_region1=c('elem1','elem2','elem3'),
                                                  id_region2=c('control','control','control'))
## Adjusted p-value for each comparison and each feature
adjusted_pval(regionsFeatures_test)

## Plotting the results of the two sample test
plotTest(regionsFeatures_test)

## Summary plot of the two sample tests
## x11(10,5)
plotSummary(regionsFeatures_test,groupby='feature',align_lab='Center')
}