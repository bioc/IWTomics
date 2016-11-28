\name{plot-IWTomicsData}
\docType{methods}
\alias{plot}
\alias{plot,IWTomicsData-method}

\title{Plot \code{"IWTomicsData"} object}

\description{
  Method to plot objects of class \code{"\link{IWTomicsData}"}. The function create graphical 
  representations of the feature measurments in each region datasets, such as aligned curves, 
  pointwise quantile curves, scatterplot and smoothed scatterplot.
}

\usage{
\S4method{plot}{IWTomicsData}(x, type='boxplot',method='pearson',
  N_regions=pmin(lengthRegions(x),
                 ifelse(type=='curves',10,ifelse(type=='pairs',1000,+Inf))),
  probs=c(0.25,0.5,0.75), average=TRUE, size=TRUE,
  id_regions_subset=idRegions(x), id_features_subset=idFeatures(x),
  log_scale=FALSE, log_shift=0, col=1+seq_along(id_regions_subset),
  plot=TRUE, ask=TRUE, xlab='Windows', ylim=NULL,...)
}

\arguments{
  \item{x}{\code{"IWTomicsData"} object.}
  \item{type}{type of plot to be drawn. Possible types are:
    \itemize{
      \item \code{"boxplot"} for the pointwise quantiles curves (default),
      \item \code{"curves"} for the (aligned) curves, 
      \item \code{"pairs"} for the scatterplot matrix of the different features 
        (same resolution needed), 
      \item \code{"pairsSmooth"} for the smoothed scatterplot matrix of the different features 
        (same resolution needed).
    }
  }
  \item{method}{correlation coefficient to be computed and plotted if type="pairs" 
    or type="pairsSmooth". Possible types are "pearson" (default), "kendall" and "spearman". 
    See \code{\link[stats]{cor}} for details about the different methods. 
  }
  \item{N_regions}{number of regions to be randomly chosen (for each region dataset) in 
    \code{"curves"} and \code{"pairs"} type plot. Default plots a maximum of 10 and 1000 curves 
    for \code{"curves"} and \code{"pairs"}, respectively.}
  \item{probs}{probabilities corresponding to the quantiles to be drawn in \code{"boxplot"} 
    type plot. Default plots quartile curves.}
  \item{average}{if \code{TRUE}, plot the mean curves.}
  \item{size}{if \code{TRUE}, plot the sample size in each position.}
  \item{id_regions_subset}{vector with the identifiers of the region datasets to be plotted.}
  \item{id_features_subset}{vector with the identifiers of the features to be plotted.}
  \item{log_scale}{plot the logarithm of the measurements.}
  \item{log_shift}{scalar with the shift to be applied before applying the logorithm, when 
    \code{log_scale=TRUE}.}
  \item{col}{vector of plotting colors for the different region datasets. If type is 
    \code{"pairsSmooth"}, only the first element of \code{col} is used.}
  \item{plot}{it \code{TRUE} (default) a plot is drawn, otherwise a list of plot data is 
    returned.}
  \item{ask}{if \code{TRUE} (default) the user is prompted before a new plot is drawn.}
  \item{xlab}{a title for the x axis.}
  \item{ylim}{the y limits of the plot.}
  \item{\dots}{additional plot parameters.}
}

\value{
  \code{plot} returns a list with components depending on the plot type. 
  
  If type is \code{"pairs"} or \code{"pairsSmooth"} the list has components:
  \item{features_plot}{matrix with the data plotted, one column for each feature considered.}
  \item{features_cor}{correlation matrix.}
  \item{type}{type of plot drawn.}
  
  If type is \code{"boxplot"} or \code{"curves"} the list has components:
  \item{x_plot}{a list of vectors with the abscissa for each feature.}
  \item{features_plot}{if type is \code{"boxplot"}, a list of matrix lists with 
    the data plotted for each feature and each region datasets. Each column 
    represents one quantile drawn. If type is \code{"curves"}, a list of matrices 
    with the data plotted for each feature. Each column represents one curve drawn.}
  \item{features_average}{if \code{average} is \code{TRUE}, a list of matrices with the 
    average curves for each feature. Each column represents the average curve for a 
    region dataset.}
  \item{features_position_size}{if \code{size} is \code{TRUE}, a list of matrices with 
    the sample size for each feature. Each column represents the sample size in each 
    position for a region dataset.}
  \item{type}{type of plot drawn.}
  \item{col}{vector of plotting colors for the different region datasets.}
  \item{col_plot}{if  type is \code{"boxplot"}, list of plotting colors for each 
  feature (quantile curves colors and shaded band colors). If type is \code{"curves"}, 
  vector of plotting colors correspondent to the different curves in \code{features_plot}.}
}

\author{
  Marzia A Cremona, Alessia Pini, Francesca Chiaromonte, Simone Vantini
}

\seealso{
  \code{\link{IWTomicsData}} for \code{"IWTomicsData"} class, constructors, accessors 
    and methods;
  \code{\link[IWTomics]{smooth}} method to smooth curves in \code{"IWTomicsData"} objects;
  \code{\link{IWTomicsTest}} for the Interval-Wise Testing.
}
\examples{
data(regionsFeatures_center)

## Plot the pointwise boxplot of the curves in the different region datasets
plot(regionsFeatures_center,type='boxplot')

## Plot all the curves in the different region datasets
plot(regionsFeatures_center,type='curves',N_regions=lengthRegions(regionsFeatures_center))

## Scatterplot of the measurements in the different region datasets
plot(regionsFeatures_center,type='pairs',N_regions=lengthRegions(regionsFeatures_center))

## Smooth scatterplot of the measurements in the 'control' region datasets
plot(regionsFeatures_center,type='pairsSmooth',id_regions_subset='control',
     N_regions=lengthRegions(regionsFeatures_center)['control'],col=5)
}