\name{smooth-IWTomicsData}
\docType{methods}
\alias{smooth}
\alias{smooth,IWTomicsData-method}

\title{Smooth curves of a \code{"IWTomicsData"} object}

\description{
  The function allows to smooth the curves in a \code{"\link{IWTomicsData}"} object, 
  to fill in gaps in the measurements, to change measurement resolutions and, when 
  region alignment is \code{"scale"}, to measure the features in the different regions 
  on the same grid.
}

\usage{
\S4method{smooth}{IWTomicsData}(x, type="locpoly", 
  id_regions_subset=idRegions(x), id_features_subset=idFeatures(x),
  resolution_new=resolution(x)[id_features_subset], 
  scale_grid=unlist(lapply(lengthFeatures(x)[id_features_subset],
    function(length) max(unlist(length[id_regions_subset])))),
  fill_gaps=TRUE, bandwidth=5, degree=3, dist_knots=10, parallel=FALSE)
}

\arguments{
  \item{x}{\code{"IWTomicsData"} object.}
  \item{type}{type of smoothing used. Possible types are:
    \itemize{
      \item \code{"locpoly"} for local polynomials (see help of 
        \code{\link[KernSmooth]{locpoly}} function),
      \item \code{"kernel"} for Nadaraya-Watson kernel smoothing with Gaussian 
        kernel (see help of \code{\link{ksmooth}} function),
      \item \code{"splines"} for regression b-splines (see 
        \code{\link[fda]{create.bspline.basis}} and \code{\link[fda]{smooth.basis}).}
      Note that \code{"splines"} is very computational expensive when regions have 
      different length and/or gaps.
    }
  }
  \item{id_regions_subset}{vector with the identifiers of the region datasets to be smoothed.}
  \item{id_features_subset}{vector with the identifiers of the features to be smoothed.}
  \item{resolution_new}{resolution used to evaluate smoothed curves. As default, measurement 
    resolution is mainteined after smoothing and smoothed curves are evaluated on the 
    same grid as the raw curves. Can be either a constant (the same resolution for all 
    features) or a vector of constants (a resolution for each feature). Not used if 
    region alignment is \code{"scale"}.}
  \item{scale_grid}{number of equally-spaced grid points over which the smoothed curves 
    are evaluated. Can be either a constant (the same grid for all features) or a vector 
    of constants (a grid for each feature). Only used if region alignment is \code{"scale"}.}
  \item{fill_gaps}{if \code{TRUE} (default), smoothing is used to fill gaps 
    (\code{NA} measurements). Can be set to \code{FALSE} only if region alignment 
    is not \code{"scale"} and measurement resolution is mainteined.}
  \item{bandwidth}{bandwidth smoothing parameter used in \code{"locpoly"} or 
    \code{"kernel"} type smoothing.}
  \item{degree}{degree of local polynomial or b-splines used in \code{"locpoly"} 
    or \code{"splines"} type smoothing.}
  \item{dist_knots}{(approximate) distance between knots defining the b-splines used 
    in \code{"splines"} type smoothing.}
  \item{parallel}{if \code{TRUE} smoothing is run in parallel on the different curves.}
}

\note{
  If the regions have different lengths and/or gaps are present in the feature 
  measurements (\code{NA} measurements), \code{"splines"} type smoothing is very 
  computational expensive. In this case we suggest to use \code{"locpoly"} 
  or \code{"kernel"} type smoothing to smooth a large number of curves.
}

\value{
  \code{smooth} returns a \code{"IWTomicsData"} object with smoothed curves and 
  updated resolutions and lengths (without the optional slot \code{test}).
}

\author{
  Marzia A Cremona, Alessia Pini, Francesca Chiaromonte, Simone Vantini
}

\seealso{
  \code{\link{IWTomicsData}} for \code{"IWTomicsData"} class, constructors, accessors and methods; 
  \code{\link[IWTomics]{plot}} method to plot \code{"IWTomicsData"} objects;
  \code{\link{IWTomicsTest}} for the Interval-Wise Testing.
}
\examples{
## -------------------------------------------------------------------------------------------
## -------------------------------------------------------------------------------------------
## CURVE ALIGNMENT CENTER
## -------------------------------------------------------------------------------------------
## -------------------------------------------------------------------------------------------
data(regionsFeatures_center)
dev.new()
plot(regionsFeatures_center,type='curves',N_regions=lengthRegions(regionsFeatures_center))

## Smooth all the curves with local polynomials
regionsFeatures_smooth=smooth(regionsFeatures_center,type='locpoly')
plot(regionsFeatures_smooth,type='curves',N_regions=lengthRegions(regionsFeatures_smooth))

## Smooth only feature 'ftr1', only region datasets 'elem3' and 'control'
regionsFeatures_smooth=smooth(regionsFeatures_center,type='locpoly',id_features_subset='ftr1',
                              id_regions_subset=c('elem3','control'))
plot(regionsFeatures_smooth,type='curves',N_regions=lengthRegions(regionsFeatures_smooth))

## Smooth only feature 'ftr1' and change its resolution
regionsFeatures_smooth=smooth(regionsFeatures_center,type='locpoly',id_features_subset='ftr1',
                              resolution_new=4000)
plot(regionsFeatures_smooth,type='curves',N_regions=lengthRegions(regionsFeatures_smooth))

## Smooth without filling gaps
regionsFeatures_smooth=smooth(regionsFeatures_center,type='locpoly',fill_gaps=FALSE)
plot(regionsFeatures_smooth,type='curves',N_regions=lengthRegions(regionsFeatures_smooth))


## -------------------------------------------------------------------------------------------
## -------------------------------------------------------------------------------------------
## CURVE ALIGNMENT SCALE
## -------------------------------------------------------------------------------------------
## -------------------------------------------------------------------------------------------
data(regionsFeatures_scale)
dev.new()
plot(regionsFeatures_scale,type='curves',N_regions=lengthRegions(regionsFeatures_scale),
                           average=FALSE,size=FALSE)

## Smooth all the curves with local polynomials
regionsFeatures_smooth=smooth(regionsFeatures_scale,type='locpoly')

## Plot of the smoothed curves
plot(regionsFeatures_smooth,type='curves',N_regions=lengthRegions(regionsFeatures_smooth),
     average=FALSE)
}