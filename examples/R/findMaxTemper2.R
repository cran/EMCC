
### $Id: findMaxTemper2.R,v 1.2 2008/02/06 04:23:02 goswami Exp $

if (!('package:EMCC' %in% search( )))
    library(EMCC, lib.loc='~/Rlib')

pdf(file.path('plots', 'findMaxTemper2.pdf'))

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
statsFuncList1 <- list(adjMatSum, modeSensitive1, entropy, maxProp)
KMeansObj      <- KMeansFuncGenerator2(-13579)
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

         ## Multimodality here would correspond to, say, the i-th
         ## coordiante pairing up with different objects, e.g., j and
         ## k
         statsFuncList2 <- 
             lapply(seq_len(sampDim), function (ii)
                {
                    function (xx)
                    {
                        ## sums the object numbers of the cluster
                        ## memebers of xx[ii]
                        sum(seq_along(xx)[abs(xx - xx[ii]) <= 0])
                    }
                })
         statsFuncList  <- c(statsFuncList1, statsFuncList2)
         findMaxTemper(nIters            = 5000,
                       statsFuncList     = statsFuncList,
                       startingVals      = startingVals,
                       logTarDensFunc    = logTarDensFunc,
                       temperLadder      = temperLadder,
                       temperLimits      = c(0.5, 20),
                       ladderLen         = nLevels,
                       levelsSaveSampFor = seq_len(nLevels),
                       doFullAnal        = TRUE,
                       saveFitness       = TRUE,
                       verboseLevel      = 2)
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
