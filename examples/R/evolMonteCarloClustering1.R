
### $Id: evolMonteCarloClustering1.R,v 1.4 2008/02/06 05:11:10 goswami Exp $

if (!('package:EMCC' %in% search( )))
    library(EMCC, lib.loc='~/Rlib')

pdf(file.path('plots', 'evolMonteCarloClustering1.pdf'))

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

dev.off( )

