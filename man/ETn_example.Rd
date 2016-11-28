\name{ETn_example}
\alias{ETn_example}
\docType{data}

\title{ETn Recombination hotspots data}

\description{
  Recombination hotspots data in the regions surrounding fixed ETns (elements of the 
  Early Transposon family of active Endogenous Retroviruses in mouse) and in control regions. 
}

\usage{data(ETn_example)}

\format{
  \code{ETn_example} is an object of class \code{"IWTomicsData"}, with two region datasets 
  "ETn fixed" and "Control" with center alignment and one feature "Recombination hotspots content". 
  
  In particular, the region dataset "ETn fixed" contains 1296 genomic regions of 64 kb 
  surrounding fixed ETns elements (32-kb flanking sequences upstream and 32-kb flanking 
  sequences downstream of each element). The region dataset "Control" contains 1142 regions
  of 64 kb without elements, used as control in the test. The regions are aligned around 
  their center (i.e. around the ETn integration sites). 
  
  Recombination hotspots measurements are associated to each "ETn fixed" and 
  "Control" region. In particular, this feature is measured in 1-kb windows, so that each region 
  is associated to a recombination hotspots curve made of 64 values. The measurement used is 
  the feature content, i.e. the fraction of the 1-kb window that is covered by recombination 
  hotspots. 
}

\value{
  \code{"IWTomicsData"} object.
}

\source{
  Data have been collected and pre-processed by: 
  R Campos-Sanchez, MA Cremona, A Pini, F Chiaromonte and KD Makova (2016). 
  Integration and fixation preferences of human and mouse endogenous retroviruses 
  uncovered with Functional Data Analysis. \emph{PLoS Computational Biology}. 12(6): 1-41.
  
  Fixed ETn positions come from:
  Y Zhang, IA Maksakova, L Gagnier, LN van de Lagemaat, DL Mager (2008). Genome-wide assessments
  reveal extremely high levels of polymorphism of two active families of mouse endogenous retroviral
  elements. \emph{PLoS Genetics}. 4: e1000007.
  
  Recombination hotspots data come from: 
  H Brunschwig, L Levi, E Ben-David, RW Williams, B Yakir, S Shifman (2012). Fine-scale maps 
  of recombination rates and hotspots in the mouse genome. \emph{Genetics}. 191: 757-764.
}

\examples{
data(ETn_example)
ETn_example

plot(ETn_example,type='boxplot')
}
\keyword{datasets}
