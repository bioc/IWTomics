\name{IWTomicsData-class}
\docType{class}
\alias{IWTomicsData-class}
\alias{IWTomicsData}

% Constructor
\alias{IWTomicsData,GRangesList,list-method}
\alias{IWTomicsData,character,data.frame-method}
\alias{IWTomicsData,character,character-method}
\alias{IWTomicsData,character,matrix-method}
% Accessors
\alias{nRegions}
\alias{nFeatures}
\alias{lengthRegions}
\alias{lengthFeatures}
\alias{resolution}
\alias{metadata}
\alias{regions}
\alias{features}
\alias{idRegions}
\alias{idFeatures}
\alias{nameRegions}
\alias{nameFeatures}
\alias{alignment}
\alias{testInput}
\alias{nTests}
\alias{idRegionsTest}
\alias{idFeaturesTest}
\alias{adjusted_pval}
% Subsetting
\alias{[,IWTomicsData,ANY,ANY,ANY-method}
% Combining
\alias{c,IWTomicsData-method}
\alias{merge,IWTomicsData,IWTomicsData-method}
\alias{rbind,IWTomicsData-method}
\alias{cbind,IWTomicsData-method}
% show method
\alias{show,IWTomicsData-method}

\title{Class \code{"IWTomicsData"}}

\description{
  The class \code{"IWTomicsData"} defines a container for storing a collection 
  of aligned genomic region datasets, and their associated feature measurements, 
  to be used as input for the Interval-Wise Testing of "Omics" data. 
}

\section{Details}{
  An object of class \code{"IWTomicsData"} is a list of genomic locations organized 
  in different region datasets and aligned (e.g. around their center). Multiple 
  genomic feature measurements are associated to each location. In particular, each 
  feature is measured in windows of fixed size inside each location. As a consequence, 
  a vector of measurement is associated to each pair of locations and features. This 
  information is stored in the following slots:
  \describe{
    \item{\code{metadata}:}{list with \code{region_datasets} and \code{feature_datasets} 
      components. The component \code{region_datasets} is a data frame with names, file 
      names and size of each region dataset. The component \code{feature_datasets} is a 
      data frame with names, file names and resolution of each feature.}
    \item{\code{regions}:}{\code{?"\link[GenomicRanges]{GRangesList}"} object containing 
      the genomic locations of each region dataset.}
    \item{\code{alignment}:}{string indicating the region alignment type. Can be 
      \code{"left"}, \code{"right"}, \code{"center"} or \code{"scale"}.}
    \item{\code{features}:}{list of matrix lists, with columns of aligned feature 
      measurments corresponding to each feature in each region dataset.}
    \item{\code{length_features}:}{list of vector lists, with the number of measurments 
      corresponding to each feature in each region dataset.}
    \item{\code{test}:}{(optional) list with \code{input} and \code{result}, containing 
      test input and results. In particular, \code{input} is a list with components:
      \itemize{
        \item{\code{id_region1}: }{identifier(s) of the region dataset(s) tested.}
        \item{\code{id_region2}: }{identifier(s) of the region dataset(s) tested for two 
          sample test.}
        \item{\code{id_features_subset}: }{vector with the identifiers of the features tested.}
        \item{\code{mu}: }{the center of symmetry under the null hypothesis in one sample test, 
          or the difference between the two populations in two sample test.}
        \item{\code{statistics}: }{test statistics used in the test.}
        \item{\code{probs}: }{probabilities corresponding to the quantiles in test statistics 
          \code{"quantile"}.}
        \item{\code{max_scale}: }{the maximum interval length used for the p-value adjustment.}
        \item{\code{paired}: }{if \code{TRUE}, the test was paired.}
        \item{\code{B}: }{number of permutation used in the test.}
      }
      Each element of the list \code{result} is a list of test results for the features 
      tested. Each test result is a a list with components:
      \itemize{
        \item{\code{test}: }{string vector indicating the type of test performed, \code{"1pop"} 
          or \code{"2pop"} for one sample and two sample tests, respectively.}
        \item{\code{mu}: }{the center of symmetry under the null hypothesis in one sample test, 
          or the difference between the two populations in two sample test, for the 
          particular test considered.}
        \item{\code{max_scale} }{the maximum interval length used for the p-value adjustment.}
        \item{\code{T0_plot}: }{value of the test statistics without squaring (used by 
          \code{\link{plotSummary}} to draw the summary plot).}
        \item{\code{adjusted_pval}: }{adjusted p-value curve, i.e. adjusted p-values for each 
          point of the curves. The adjustment is done considering \code{max_scale} as length.}
        \item{\code{adjusted_pval_matrix}: }{matrix of size the number of points in the curves 
          with the adjusted p-value curves for each possible scale up to \code{max_scale}. 
          Row \code{i} of the matrix contains the adjusted p-value curve with correction 
          done up to scale \code{p-i+1} (the matrix contains \code{NA} for scale greater 
          than \code{max_scale}).}
        \item{\code{unadjusted_pval}: }{p-value curve, i.e. raw p-values for each point of the curves.}
        \item{\code{pval_matrix}: }{matrix of size the number of points in the curves with the 
          raw p-values of the multivariate tests. The element \code{(i,j)} of the matrix 
          contains the p-value of the joint NPC test of the components \code{j,j+1,...,j+(p-i)} 
          (the matrix contains \code{NA} for scale greater than \code{max_scale}).}
        \item{\code{exact}: }{logical value indicating whether the exact p-values have been computed.}
        \item{\code{notNA}: }{vector of logical vectors indicated the points of the curves where 
          the test was performed (used by \code{\link{plotSummary}} to draw the summary plot).}
      }
    }
    Objects of class \code{"IWTomicsData"} can be initialized from BED or Table files, 
    or they can be directly created by supplying a \code{"GRangesList"} of genomic region 
    datasets and a \code{"list"} with the aligned feature measurements (see Contructors).
    The optional slot \code{test} is filled by the function \code{\link{IWTomicsTest}} 
    that performs the Interval-Wise Testing. 
  }
}

\section{Constructors}{
  \describe{
    \item{}{\code{IWTomicsData(x, y, alignment='center',}
            \code{id_regions=NULL,} \code{name_regions=NULL,} \code{id_features=NULL,}
            \code{name_features=NULL,} \code{path=NULL,} \code{start.are.0based=TRUE,}
            \code{header=FALSE, ...)}:
            creates a \code{"IWTomicsData"} object from BED or Table files.}
    \item{}{\code{IWTomicsData(x, y, alignment='center',}
            \code{id_regions=NULL,} \code{name_regions=NULL,} \code{id_features=NULL,}
            \code{name_features=NULL,} \code{length_features=NULL)}:
            creates a \code{"IWTomicsData"} object from genomic regions datasets and 
            feature measurements.
      \describe{
        \item{\code{x}}{vector with the names of the region files containing 
          the region datasets to be loaded. BED and Table formats currently supported.
          \emph{Alternative constructor}: \code{"\link[GenomicRanges]{GRangesList}"} object 
          with genomic locations of each region dataset.}
        \item{\code{y}}{vector with the names of the feature files, or 
          dataframe with columns of feature file names corresponding to the different 
          region datasets. Each feature must be measured in windows of a fixed size inside 
          all the regions. BED and Table formats currently supported: either a row (with 
          4 columns chr, start, end, measure) for each window, or a row for each region 
          (with columns chr, start, end, value1, ..., valueN). Note that all files must 
          be sorted. \emph{Alternative constructor}: list of matrix lists, with columns 
          of aligned feature measurments corresponding to each feature in each region dataset.}
        \item{\code{alignment}}{region alignment. Possible types are: 
          \itemize{
            \item "\code{left}" for alignment of the starting positions, 
            \item "\code{right}" for alignment of the ending positions,
            \item "\code{center}" for alignment of the central positions (default), 
            \item "\code{scale}" for scaling all regions to the same length.}}
        \item{\code{id_regions}}{vector with the identifiers of the region datasets. 
          If \code{NULL}, \code{file_regions} or \code{names(regions)} are used.}
        \item{\code{name_regions}}{vector with the names of the region datasets to be used 
          in the output. If \code{NULL}, the identifiers \code{id_regions} are used.}
        \item{\code{id_features}}{vector with the identifiers of the features. If \code{NULL}, 
          \code{file_features} or \code{names(features)} are used.}
        \item{\code{name_features}}{vector with the names of the features to be used in 
          the output plots. If \code{NULL}, the identifiers \code{id_features} are used.}
        \item{\code{path}}{the directory that contains the files. If \code{NULL}, the 
          current working directory is used.}
        \item{\code{start.are.0based}}{if \code{TRUE} (default) the start position in the 
          region files are considered to be 0-based, and converted to 1-based in the 
          \code{"IWTomicsData"} object in output.}
        \item{\code{header}}{\code{TRUE} or \code{FALSE} (default) indicating if the files contain 
          the names of the variables as their first lines.}
        \item{\code{length_features}}{list of vector lists, with the number of measurments 
          corresponding to each feature in each region dataset.}
        \item{\dots}{additional parameters in input to \code{\link[utils]{read.delim}}.}
      }
    }
  }
}

\section{Accessors}{
  In the following code, \code{x} is a \code{"IWTomicsData"} object.
  \describe{
    \item{}{
      \code{nRegions(x)}: get the number of region datasets.
    }
    \item{}{
      \code{nFeatures(x)}: get the number of features.
    }
    \item{}{
      \code{dim(x)}: get the dimension of the object (number of region datasets, number of
      features).
    }
    \item{}{
      \code{lengthRegions(x)}: get the number of locations in each region dataset.
    }
    \item{}{
      \code{lengthFeatures(x)}: get a list of vector list, with the number of measurements 
      corresponding to each feature in each region dataset.
    }
    \item{}{
      \code{resolution(x)}: get the measurement resolution for each feature.
    }
    \item{}{
      \code{metadata(x)}: get the metadata associated with the object, i.e. a list with 
      \code{region_datasets} and \code{feature_datasets} components.
    }
    \item{}{
      \code{regions(x)}: get the \code{"GRangesList"} object containing the genomic 
      locations of each region dataset.
    }
    \item{}{
      \code{features(x)}: get a list of matrix lists, with columns of aligned feature 
      measurements corresponding to each feature in each region dataset.
    }
    \item{}{
      \code{idRegions(x)}: get the identifiers of the region datasets.
    }
    \item{}{
      \code{idFeatures(x)}: get the identifiers of the features.
    }
    \item{}{
      \code{nameRegions(x)}: get the names of the region datasets.
    }
    \item{}{
      \code{nameFeatures(x)}: get the names of the features.
    }
    \item{}{
      \code{alignment(x)}: get the region alignment.
    }
    \item{}{
      \code{testInput(x)}: get the test input (if present).
    }
    \item{}{
      \code{nTests(x)}: get the number of tests present.
    }
    \item{}{
      \code{idRegionsTest(x)}, \code{idRegionsTest(x,test)}: get the identifiers of the 
      region datasets in the different tests. The (optional) argument \code{test} 
      indicates the indices of the tests to be considered.
    }
    \item{}{
      \code{idFeatuersTest(x)}: get the identifiers of the features tested.
    }
    \item{}{
      \code{adjusted_pval(x)}, \code{adjusted_pval(x,test,id_features_subset,scale_threshold)}: 
      get the adjusted p-values of the different tests. The (optional) argument \code{test} 
      indicates the indices of the tests to be considered. The (optional) argument 
      \code{id_features_subset} is a vector with the identifiers of the features to be 
      cosidered. The (optional) argument \code{scale_threshold} is the threshold on 
      the test scale (maximum interval length for the p-value adjustment) for the 
      adjusted p-value computation. Can be either a scalar (the same length for all features) 
      or a vector (a length for each feature) or a list of vectors (a vector for each test). 
      See \code{\link{IWTomicsTest}} for more details.
    }
  }
}

\section{Subsetting}{
  In the following code, \code{x} is a \code{"IWTomicsData"} object. The optional slot 
  \code{test}, if present in \code{x}, is deleted when using subsetting methods.
  \describe{
    \item{}{
      \code{x[i,j]}: extract region dataset \code{i} and feature \code{j} in a new 
      \code{"IWTomicsData"} object. Both \code{i} and \code{j} can be logical vectors, 
      numeric vectors, character vectors (with region dataset and feature identifiers, 
      respectively), or missing.
    }
  }
}

\section{Combining}{
  In the following code, \code{x} is a \code{"IWTomicsData"} object. The optional slot 
  \code{test}, if present in \code{x}, is deleted when using combining methods.
  \describe{
    \item{}{
      \code{c(x,...)} and \code{merge(x,...)}: create a new \code{"IWTomicsData"} 
      object combining \code{x} with the \code{"IWTomicsData"} objects in \code{...}. 
      Any object in \code{...} must have the same region alignment as \code{x}, 
      and region datasets and features present in multiple objects must coincide. 
    }
    \item{}{
      \code{rbind(x,...)}: create a new \code{"IWTomicsData"} object combining the 
      features in \code{x} with the features in the \code{"IWTomicsData"} objects \code{...}. 
      Region datasets in \code{x} and any object in \code{...} must coincide and have 
      the same region alignment. 
    }
    \item{}{
      \code{cbind(x,...)}: create a new \code{"IWTomicsData"} object combining the 
      region datasets in \code{x} with the region datasets in the \code{"IWTomicsData"} 
      objects \code{...}. Features in \code{x} and any object in \code{...} must coincide 
      and the \code{"IWTomicsData"} objects must have the same region alignment. 
    }
	}
}

\section{Other methods}{
  \describe{
    \item{}{
      \code{show(x)}:
      The \code{show} method prints the number of region datasets, their alignment type 
      and the number of features in the \code{"IWTomicsData"} object. It also displays 
      names and size of the region datasets, and names and resolution of the features. 
      If the slot \code{test} is present in \code{x}, the \code{show} method prints also 
      the comparisons present. 
    }
  }
}

\author{
  Marzia A Cremona, Alessia Pini, Francesca Chiaromonte, Simone Vantini
}

\seealso{
  \code{\link[IWTomics]{plot}} method to plot \code{"IWTomicsData"} objects;
  \code{\link[IWTomics]{smooth}} method to smooth curves in \code{"IWTomicsData"} objects;
  \code{\link{IWTomicsTest}} for the Interval-Wise Testing.
}

\examples{
examples_path <- system.file("extdata",package="IWTomics")
datasets=read.table(file.path(examples_path,"datasets.txt"),
                    sep="\t",header=TRUE,stringsAsFactors=FALSE)
features_datasetsBED=read.table(file.path(examples_path,"features_datasetsBED.txt"),
                                sep="\t",header=TRUE,stringsAsFactors=FALSE)
features_datasetsTable=read.table(file.path(examples_path,"features_datasetsTable.txt"),
                                  sep="\t",header=TRUE,stringsAsFactors=FALSE)
data(regions_example)
data(features_example)

## -------------------------------------------------------------------------------------------
## CONSTRUCTION
## -------------------------------------------------------------------------------------------
## Get genomic regions for four region datasets, 
## and two features for each region dataset

## From BED files (check for consistency, time consuming) 
regionsFeaturesBED=IWTomicsData(datasets$regionFile,features_datasetsBED[,3:6],
                                'center',datasets$id,datasets$name,
                                features_datasetsBED$id,features_datasetsBED$name,
                                path=file.path(examples_path,'files'))
regionsFeaturesBED

## From Table files (less checks for consistency, more efficient)
regionsFeaturesTable=IWTomicsData(datasets$regionFile,features_datasetsTable[,3:6],
                                  'center',datasets$id,datasets$name,
                                  features_datasetsTable$id,features_datasetsTable$name,
                                  path=file.path(examples_path,'files'))
regionsFeaturesTable

## From genomic regions datasets and feature measurements.
regionsFeatures=IWTomicsData(regions_example,features_example,alignment='center')
regionsFeatures

## -------------------------------------------------------------------------------------------
## SUBSETTING 
## -------------------------------------------------------------------------------------------
## Extract a subset of region datasets and/or of features

## Get the first region dataset and the second features
regionsFeaturesBED[1,2]

## Get the first region dataset and the second features, using identifiers
regionsFeaturesBED['elem1','ftr2']

## Get the first two region datasets for all the features
regionsFeaturesBED[1:2,]

## Get all region datasets for the first feature
regionsFeaturesBED[,1]

## -------------------------------------------------------------------------------------------
## COMBINING 
## -------------------------------------------------------------------------------------------
data1=regionsFeaturesBED[1:2,1]
data2=regionsFeaturesBED[1:2,2]
data3=regionsFeaturesBED[2:3,]
data4=regionsFeaturesBED[4,]

## Merge different objects
data1
data2
c(data1,data2)
merge(data1,data2)

## Combine different features
data1
data2
cbind(data1,data2)

## Combine different regions
data3
data4
rbind(data3,data4)

## Combine methods together
rbind(cbind(data1,data2),data3,data4)
}
\keyword{classes}