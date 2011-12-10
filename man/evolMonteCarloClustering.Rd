
%%% $Id: evolMonteCarloClustering.Rd,v 1.13 2008/02/06 05:40:31 goswami Exp $

\name{evolMonteCarloClustering}

\alias{evolMonteCarloClustering}

\title{evolutionary Monte Carlo clustering algorithm}

\description{
  Given a possibly multi-modal and multi-dimensional clustering target
  density function and a temperature ladder this function produces samples
  from the target using the evolutionary Monte Carlo clustering (EMCC) 
  algorithm.

  Below \code{sampDim} refers to the dimension of the sample space,
  \code{temperLadderLen} refers to the length of the temperature ladder,
  and \code{levelsSaveSampForLen} refers to the length of the
  \code{levelsSaveSampFor}.
}

\usage{
evolMonteCarloClustering(nIters,              
                         temperLadder,              
                         startingVals,              
                         logTarDensFunc,            
                         MHMergeProb       = 0.5,  
                         moveProbsList     = NULL,  
                         moveNTimesList    = NULL,  
                         levelsSaveSampFor = NULL,  
                         saveFitness       = FALSE, 
                         verboseLevel      = 0,     
                         \dots)                       
}

\arguments{
  \item{nIters}{\code{integer} \eqn{>} 0.}

  \item{temperLadder}{\code{double} vector with all \emph{positive}
    entries, in \emph{decreasing} order.}

  \item{startingVals}{\code{double} matrix of dimension
    \code{temperLadderLen} \eqn{\times}{x} \code{sampDim} or vector of
    length \code{sampDim}, in which case the same starting values are
    used for every temperature level.}

  \item{logTarDensFunc}{\code{function} of two arguments
    \code{(draw, \dots)}
    that returns the target density evaluated in the log scale.}
  
  \item{MHMergeProb}{\code{double} in (0, 1). \emph{See details below
      for the use of this argument.}}
  
  \item{moveProbsList}{named \code{list} of probabilities adding upto 1.}

  \item{moveNTimesList}{named \code{list} of integers \eqn{\ge}{>=} 0.}

  \item{levelsSaveSampFor}{\code{integer} vector with \emph{positive}
    entries.}

  \item{saveFitness}{\code{logical}.}

  \item{verboseLevel}{\code{integer}, a value \eqn{\ge}{>=} 2 produces a
    lot of output.}

  \item{\dots}{optional arguments to be passed to
    \code{logTarDensFunc}.}

}

\details{
  \describe{

    \item{The EMCC algorithm}{The evolutionary Monte Carlo clustering
      (EMCC; Goswami and Liu, 2007) algorithm is composed of the
      following moves:}
    
    \tabular{rl}{
      \acronym{MH}
      \tab Metropolis-Hastings or mutation
      \cr

      \acronym{SCSC_ONE_NEW}
      \tab sub-cluster swap crossover: one new
      \cr

      \acronym{SCSC_TWO_NEW}
      \tab sub-cluster swap crossover: two new
      \cr

      \acronym{SCRC}
      \tab sub-cluster reallocation crossover
      \cr
      
      \acronym{RE}
      \tab (random) exchange 
    }

    The current function could be used to run the EMCC algorithm by
    specifying what moves to employ using the following variables.

    \describe{
      
      \item{\code{moveProbsList} and \code{moveNTimesList}}{The allowed
        names for components of \code{moveProbsList} and
        \code{moveNTimesList} come from the abbreviated names of the
        moves above. For example, the following specifications are
        valid:}

      \preformatted{
moveProbsList = list(MH           = 0.5,
                     SCSC_TWO_NEW = 0.25,
                     SCRC         = 0.25)
      }

      \preformatted{
moveNTimesList = list(MH           = 1, 
                      SCSC_TWO_NEW = floor(temperLadderLen / 2),
                      SCRC         = floor(temperLadderLen / 2),
                      RE           = temperLadderLen)           
      }
    }

    \item{\code{MHMergeProb}}{In the \acronym{MH} or the mutation step,
      each of the \code{sampDim}-many objects are proposed to either
      merge with an existing cluster or split to form its own cluster
      with probability \code{MHMergeProb} and (1 - \code{MHMergeProb}),
      respectively (see Goswami and Liu, 2007).}

    \item{\code{levelsSaveSampFor}}{By default, samples are saved and
      returned for temperature level \code{temperLadderLen}. The
      \code{levelsSaveSampFor} could be used to save samples from other
      temperature levels as well (e.g., \code{levelsSaveSampFor =
      1:temperLadderLen} saves samples from all levels).}

    \item{\code{saveFitness}}{The term \emph{fitness} refers to the
      function \eqn{H(x)}, where the target density of interest is given
      by:

      \deqn{g(x) \propto \exp[ -H(x) / \tau_{min} ]}

      \eqn{H(x)} is also known as the \emph{energy} function. By default,
      the fitness values are not saved, but one can do so by setting
      \code{saveFitness = TRUE}.}
    
    }
  }

\value{
  This function returns a list with the following components:

  \item{draws}{\code{array} of dimension \code{nIters} \eqn{\times}{x}
    \code{sampDim} \eqn{\times}{x} \code{levelsSaveSampForLen}, if
    \code{saveFitness = FALSE}. If \code{saveFitness = TRUE}, then the
    returned array is of dimension \code{nIters} \eqn{\times}{x}
    \code{(sampDim + 1)} \eqn{\times}{x} \code{levelsSaveSampForLen};
    i.e., each of the \code{levelsSaveSampForLen} matrices contain the
    fitness values in their last column.}

  \item{acceptRatios}{\code{matrix} of the acceptance rates for various
    moves used.}

  \item{detailedAcceptRatios}{\code{list} of matrices with detailed
    summary of the acceptance rates for various moves used.}

  \item{nIters}{the \code{nIters} argument.}

  \item{temperLadder}{the \code{temperLadder} argument.}

  \item{startingVals}{the \code{startingVals} argument.}

  \item{moveProbsList}{the \code{moveProbsList} argument.}

  \item{moveNTimesList}{the \code{moveNTimesList} argument.}  

  \item{levelsSaveSampFor}{the \code{levelsSaveSampFor} argument.}

  \item{time}{the time taken by the run.}
  
}

\note{
  The effect of leaving the default value \code{NULL} for some of the
  arguments above are as follows:

  \tabular{rl}{
    \code{moveProbsList}
    \tab \code{list(MH = 0.5, RC = 0.25, 'SCSC_TWO_NEW' = 0.25)}.
    \cr

    \code{moveNTimesList}
    \tab \code{list(MH = 1, RC = mm, 'SCSC_TWO_NEW' = mm, RE = nn)},  
    \cr

    \tab
    where \code{mm <- floor(nn / 2)} and \code{nn <- temperLadderLen}.
    \cr

    \code{levelsSaveSampFor}
    \tab \code{temperLadderLen}.

  }
}

\author{Gopi Goswami \email{goswami@stat.harvard.edu}}

\references{
  \cite{Gopi Goswami, Jun S. Liu and Wing H. Wong (2007). Evolutionary
    Monte Carlo Methods for Clustering. \emph{Journal of Computational
    and Graphical Statistics, 16:4:855-876}.}
}

\examples{
\dontrun{
## The following example is a simple stochastic optimization problem,
## the set up is same as that of findMaxTemper and placeTempers. Here
## no "heating up" is necessary, and hence the maximum temprature is
## the coldest one, namely, 0.5.
##
## However, we run evolMonteCarloClustering on this example with a
## temperature ladder that is the output of placeTempers, which
## assumes that the maximum temperature is 5.
KMeansObj  <- KMeansFuncGenerator1(-97531)
samplerObj <-
    with(KMeansObj,
     {
         temperLadder    <- c(5.0000000, 1.5593974, 1.1028349, 0.9220684,
                              0.7900778, 0.6496648, 0.5135825, 0.5000000)
         nLevels         <- length(temperLadder)
         sampDim         <- nrow(yy)
         startingVals    <- sample(c(0, 1),
                                   size    = nLevels * sampDim,
                                   replace = TRUE)
         startingVals    <- matrix(startingVals, nrow = nLevels, ncol = sampDim)
         moveProbsList   <- list(MH             = 0.4,
                                 RC             = 0.3,
                                 'SCSC_TWO_NEW' = 0.3)
         mm              <- floor(nLevels / 2)
         moveNTimesList  <- list(MH             = 1,
                                 RC             = mm,
                                 'SCSC_TWO_NEW' = mm,
                                 RE             = nLevels)
         evolMonteCarloClustering(nIters            = 5000,
                                  temperLadder      = temperLadder,
                                  startingVals      = startingVals,
                                  logTarDensFunc    = logTarDensFunc,
                                  moveProbsList     = moveProbsList,
                                  moveNTimesList    = moveNTimesList,
                                  levelsSaveSampFor = seq_len(nLevels),
                                  saveFitness       = TRUE,
                                  verboseLevel      = 1)
     })
print(samplerObj)
print(names(samplerObj))
with(c(samplerObj, KMeansObj),
 {
     print(acceptRatios)
     print(detailedAcceptRatios)
     print(dim(draws))
     fitnessCol <- ncol(draws[ , , 1])     
     sub        <- paste('uniform prior on # of clusters: DU[',
                         priorMinClusters, ', ',
                         priorMaxClusters, ']', sep = '')
     for (ii in rev(seq_along(levelsSaveSampFor))) {
         main <- paste('EMCC (MAP) clustering (temper = ',
                       round(temperLadder[levelsSaveSampFor[ii]], 3), ')',
                       sep = '')
         MAPRow <- which.min(draws[ , fitnessCol, ii])
         clusterPlot(clusterInd        = draws[MAPRow, -fitnessCol, ii],
                     data              = yy,
                     main              = main,
                     sub               = sub,
                     knownClusterMeans = knownClusterMeans)
     }
 })
}
}
\keyword{methods}

