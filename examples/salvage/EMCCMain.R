
###  $Id: EMCCMain.R,v 1.7 2008/02/03 04:18:56 goswami Exp $
###  
###     File: EMCCMain.R
###  Package: EMC
###  
###  Copyright (C) 2006-present Gopi Goswami
###
###  This program is free software; you can redistribute it and/or modify
###  it under the terms of the GNU General Public License as published by
###  the Free Software Foundation; either version 2 of the License, or
###  (at your option) any later version.
###
###  This program is distributed in the hope that it will be useful,
###  but WITHOUT ANY WARRANTY; without even the implied warranty of
###  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
###  GNU General Public License for more details.
###
###  For a copy of the GNU General Public License please write to the
###  Free Software Foundation, Inc.
###  59 Temple Place, Suite 330.
###  Boston, MA  02111-1307 USA.
###
###  For bugs in the code please contact:
###  <goswami@stat.harvard.edu>
###
###  SYNOPSIS
###
###
###
###  DESCRIPTION
###
###
###



if ('package:EMCC' %in% search( ))
    detach('package:EMCC')
library(EMCC, lib.loc='~/Rlib')


KMeansFuncGenerator <-
    function (seed = 13579, plotIt = TRUE)
{
    temperLadder <- c(40, 20, 10, 3, 1)
    nLevels <- length(temperLadder)
    tmp <- unname(as.matrix(read.table('../data/KMeans-data.txt', header = FALSE)))
    ##    tmp <- tmp[(tmp[ ,1] == 4 | tmp[ ,1] == 3), ]
    origLabels <- tmp[ , 1]
    yy <- tmp[ , -1]
    clusterPlot <-
        function (id, data, plotFname = NULL, main = "")
        {
            if (!is.null(plotFname)) 
                pdf(plotFname)

            uu <- unique(id)
            plot(data, type = 'n',
                 main = main,
                 xlab = as.expression(substitute(x[ll], list(ll = 1))),
                 ylab = as.expression(substitute(x[ll], list(ll = 2))))
            for (ii in 1:length(uu)) {
                cc <- as.character(ii)
                points(data[(id == uu[ii]), ], pch = cc, col = ii) 
            }
            if (!is.null(plotFname))
                dev.off( )
        }

    if (plotIt) {
        par(mfcol = c(2, 2))
        clusterPlot(origLabels, yy, main = "Original data")
    }


    ## helper function for the logTarDens func
    partitionRep <-
        function (cl)
        {
            stopifnot(is.vector(cl))
            nn <- length(cl)
            stopifnot(all(cl >= 0))
            stopifnot(all(cl < nn))
            cl <- floor(cl)
            uu <- unique(cl)
            luu <- length(uu)
            id <- 1:nn
            ll <- sapply(1:luu, FUN =
                         function (ii)
                     {
                         list(id[cl == uu[ii]])
                     })
            return(list(clusterLabels = uu, clusters = ll))
        }

    logTarDensFunc <-
        function (xx)
        {
            xx <- as.integer(xx)
   #         print(xx)
            scr <- partitionRep(xx)
            sum <- 0
            nClusters <- length(scr$clusterLabels)
            sapply(1:nClusters, FUN =
                   function (ii)
               {
                   clusterVals <- as.matrix(yy[scr$clusters[[ii]], ])
                   clusterMean <- c(apply(clusterVals, 2, mean))
                   apply(clusterVals, 1, FUN =
                         function (xx)
                     {
                         tmp <- xx - clusterMean
                         sum <<- sum + c(tmp %*% tmp)
                         return(0)
                     })
                   return(0)
               })
            priorMinClusters <- 2; priorMaxClusters <- 2
            sum <- -sum +
                ifelse(((nClusters < priorMinClusters) ||
                        (nClusters > priorMaxClusters)),
                       -Inf, -log(priorMaxClusters - priorMinClusters + 1))
    #         print(sum)
            return(sum)
        }
    sampDim <- length(origLabels)
    ss <- matrix(rep(c(0:1, rep(0, sampDim - 2)), nLevels),
                 nrow = nLevels, byrow = TRUE)
    moveProbsList <- list(MH = 0.3, 'SCSC_ONE_NEW' = 0.0,
                          'SCSC_TWO_NEW' = 0.7, SCRC = 0.0,
                          RC = 0.0, BCE = 0.0, BIRE = 0.0, BSE = 0.0)
    moveNTimesList <- list('SCSC_ONE_NEW' = 0, 'SCSC_TWO_NEW' = 10,
                           SCRC = 0, RE = 0, BCE = 0, BIRE = 0,
                           BSE = 0, CE = 0)
    return(list(temperLadder = temperLadder,
                startingVals = ss,
                origLabels = origLabels,
                yy = yy,
                logTarDensFunc = logTarDensFunc,
                moveProbsList = moveProbsList,
                moveNTimesList = moveNTimesList))
}

KMeansObj <- KMeansFuncGenerator(plotIt = TRUE)
KMeansEMCC <-
  evolMonteCarloClustering(nIters = 100,
                           temperLadder = KMeansObj$temperLadder,
                           startingVals = KMeansObj$startingVals,
                           logTarDensFunc = KMeansObj$logTarDensFunc,
                           moveProbsList = KMeansObj$moveProbsList,
                           moveNTimesList = KMeansObj$moveNTimesList,
                           levelsSaveSampFor = c(5),
                           saveFitness = TRUE,
                           verboseLevel = 1)
sampDim <- ncol(KMeansEMCC$draws[ , , 1]) - 1
fitness <- KMeansEMCC$draws[ , sampDim + 1, 1]
oo <- order(fitness, decreasing = FALSE)
bestClustering <- KMeansEMCC$draws[oo[1], 1:sampDim, 1]
clusterPlot(bestClustering, KMeansData[ , -1], main = "Clustered data")




