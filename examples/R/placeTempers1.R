
### $Id: placeTempers1.R,v 1.3 2008/02/05 20:29:20 goswami Exp $

if (!('package:EMCC' %in% search( )))
    library(EMCC, lib.loc='~/Rlib')

pdf(file.path('plots', 'placeTempers1.pdf'))

## The following example is a simple stochastic optimization problem,
## the set up is same as that of findMaxTemper. Here no "heating up"
## is necessary, and hence the maximum temprature is the coldest one,
## namely, 0.5.
##
## However, we do the temperature placement to show how placeTempers
## works, assuming the maximum temperature is 5.
KMeansObj       <- KMeansFuncGenerator1(-97531)
placeTempersObj <-
    with(KMeansObj,
     {
         nLevels      <- 15
         sampDim      <- nrow(yy)
         startingVals <- sample(c(0, 1),
                                size    = nLevels * sampDim,
                                replace = TRUE)
         startingVals <- matrix(startingVals, nrow = nLevels, ncol = sampDim)
         placeTempers(nIters            = 5000,
                      acceptRatioLimits = c(0.5, 0.6),
                      ladderLenMax      = 50,
                      startingVals      = startingVals,
                      logTarDensFunc    = logTarDensFunc,
                      temperLimits      = c(0.5, 5),
                      ladderLen         = nLevels,
                      scheme            = 'geometric',
                      levelsSaveSampFor = seq_len(nLevels),
                      saveFitness       = TRUE,
                      verboseLevel      = 1)
     })
print(placeTempersObj)
print(names(placeTempersObj))
with(c(placeTempersObj, KMeansObj),
 {
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

dev.off( )
