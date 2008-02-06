
%%% $Id: utilsForExamples.Rd,v 1.5 2008/02/06 05:40:31 goswami Exp $

\name{utilsForExamples}

\alias{partitionRep}
\alias{clusterPlot}
\alias{KMeansFuncGenerator1}
% \alias{KMeansFuncGenerator2}

\title{The utility function(s) for examples}

\description{
  The utility function(s) that are used in the example sections of the
  exported functions in this package.
}

\usage{
partitionRep(clusterInd)
clusterPlot(clusterInd,                    
            data,                  
            main              = '',     
            sub               = '',     
            knownClusterMeans = NULL,   
            \dots) 
KMeansFuncGenerator1(seed, plotIt = TRUE)
% KMeansFuncGenerator2(seed, plotIt = TRUE)
}

\arguments{

  \item{clusterInd}{vector of cluster indicators.}

  \item{data}{a matrix with two columns representing the two-dimensional
    data clustered by clusterInd.}

  \item{main}{the title of the plot.}

  \item{sub}{the sub-title of the plot.}

  \item{knownClusterMeans}{a matrix with two columns (for the two
    dimensions), the rows containing the cluster means. These are
    plotted when provided.}

  \item{seed}{the seed for random number generation.}

  \item{plotIt}{\code{logical}, controls the plotting of the generated
    data.}

  \item{\dots}{optional arguments to be passed to \code{plot}; see its
    documentation.}

}

\value{

  \item{partitionRep}{this function returns a list with two components,
    namely, \code{clusterLabels} (the unique cluster identifiers) and
    \code{clusters} (the partitioning of the cluster identifiers), as a
    list.}

  \item{KMeansFuncGenerator1}{this function returns a list containing
    the objects to be used as arguments to the exported functions in the
    respective example sections of this package.}

}

\author{Gopi Goswami \email{goswami@stat.harvard.edu}}

\seealso{\code{\link{evolMonteCarloClustering}},
\code{\link{findMaxTemper}}, \code{\link{placeTempers}}}

\keyword{datagen}

