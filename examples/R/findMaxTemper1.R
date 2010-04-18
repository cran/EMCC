
### $Id: findMaxTemper1.R,v 1.4 2008/02/06 04:22:03 goswami Exp $

if (!('package:EMCC' %in% search( )))
    library(EMCC, lib.loc='~/Rlib')

pdf(file.path('plots', 'findMaxTemper1.pdf'))

## The following example is a simple stochastic optimization problem,
## and thus it does not require any "heating up", and hence the
## maximum temperature turns out to be the coldest one, i.e, 0.5.
adjMatSum <-
    function (xx)
{
    xx     <- as.integer(xx)
    adjMat <- outer(xx, xx, function (id1, id2) { id1 == id2 })
    sum(adjMat)
}
modeSensitive1 <-
    function (xx)
{
    with(partitionRep(xx),
     {
         rr   <- 1 + seq_along(clusterLabels)
         freq <- sapply(clusters, length)
         oo   <- order(freq, decreasing = TRUE)
         sum(sapply(clusters[oo], sum) * log(rr))
     })
}
entropy <-
    function (xx)
{
    yy <- table(as.vector(xx, mode = "numeric"))
    zz <- yy / length(xx)
    -sum(zz * log(zz))
}
maxProp <- 
    function (xx)
{
    yy <- table(as.vector(xx, mode = "numeric"))
    oo <- order(yy, decreasing = TRUE)
    yy[oo][1] / length(xx)
}
statsFuncList  <- list(adjMatSum, modeSensitive1, entropy, maxProp)
KMeansObj      <- KMeansFuncGenerator1(-97531)
maxTemperObj   <-
    with(KMeansObj,
     {
         temperLadder <- c(20, 10, 5, 1, 0.5)
         nLevels      <- length(temperLadder)
         sampDim      <- nrow(yy)
         startingVals <- sample(c(0, 1),
                                size    = nLevels * sampDim,
                                replace = TRUE)
         startingVals <- matrix(startingVals, nrow = nLevels, ncol = sampDim)
         findMaxTemper(nIters            = 5000,
                       statsFuncList     = statsFuncList,
                       temperLadder      = temperLadder,
                       startingVals      = startingVals,
                       logTarDensFunc    = logTarDensFunc,
                       levelsSaveSampFor = seq_len(nLevels),
                       doFullAnal        = TRUE,
                       saveFitness       = TRUE,
                       verboseLevel      = 1)
     })
print(maxTemperObj)
print(names(maxTemperObj))
with(c(maxTemperObj, KMeansObj),
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
