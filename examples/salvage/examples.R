
examplesEvolMonteCarloClustering <-
    function ( )
{
    if ('package:EMCC' %in% search( ))
        detach('package:EMCC')
    library(EMCC, lib.loc='~/Rlib')

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

    ## The KMeans data set example 
    KMeansFuncGenerator <-
        function (seed = 13579, plotIt = TRUE)
        {
            set.seed(seed)                
            temperLadder <- c(20, 10, 3, 1)
            nLevels <- length(temperLadder)
            data('KMeansData')
            origLabels <- KMeansData[ , 1]
            yy <- KMeansData[ , -1]

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
                    priorMinClusters <- 3; priorMaxClusters <- 3
                    sum <- -sum +
                        ifelse(((nClusters < priorMinClusters) ||
                                (nClusters > priorMaxClusters)),
                               -Inf, -log(priorMaxClusters - priorMinClusters + 1))
                    return(sum)
                }
            
            sampDim <- length(origLabels)
            ss <- matrix(rep(c(0:2, rep(0, sampDim - 3)), nLevels),
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
        evolMonteCarloClustering(nIters = 300,
                                 temperLadder = KMeansObj$temperLadder,
                                 startingVals = KMeansObj$startingVals,
                                 logTarDensFunc = KMeansObj$logTarDensFunc,
                                 moveProbsList = KMeansObj$moveProbsList,
                                 moveNTimesList = KMeansObj$moveNTimesList,
                                 levelsSaveSampFor = length(KMeansObj$temperLadder),
                                 saveFitness = TRUE,
                                 verboseLevel = 1)
    sampDim <- ncol(KMeansEMCC$draws[ , , 1]) - 1
    fitness <- KMeansEMCC$draws[ , sampDim + 1, 1]
    oo <- order(fitness, decreasing = FALSE)
    bestClustering <- KMeansEMCC$draws[oo[1], 1:sampDim, 1]
    clusterPlot(bestClustering, KMeansData[ , -1], main = "Clustered data")    
}


examplesFindMaxTemper <-
    function ( )
{
    if ('package:EMCC' %in% search( ))
        detach('package:EMCC')
    library(EMCC, lib.loc='~/Rlib')

    ## The KMeans data set example 
    KMeansFuncGenerator <-
        function (seed = 13579, plotIt = TRUE)
        {
            set.seed(seed)                
            temperLadder <- c(20,  10,  5,  1)
            nLevels <- length(temperLadder)
            data('KMeansData')
            origLabels <- KMeansData[ , 1]
            yy <- KMeansData[ , -1]

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
                    priorMinClusters <- 3; priorMaxClusters <- 3
                    sum <- -sum +
                        ifelse(((nClusters < priorMinClusters) ||
                                (nClusters > priorMaxClusters)),
                               -Inf, -log(priorMaxClusters - priorMinClusters + 1))
                    return(sum)
                }            
            entropy <-
                function (xx)
                {
                    yy <- table(as.vector(xx, mode = "numeric"))
                    zz <- yy / sum(yy)
                    return(-sum(zz * log(zz)))
                }
            statsFuncList <- list(entropy)
            sampDim <- length(origLabels)
            ss <- matrix(rep(c(0:2, rep(0, sampDim - 3)), nLevels),
                         nrow = nLevels, byrow = TRUE)
            return(list(temperLadder = temperLadder,
                        startingVals = ss,
                        origLabels = origLabels,
                        yy = yy,
                        logTarDensFunc = logTarDensFunc,
                        statsFuncList = statsFuncList))
        }

    KMeansObj <- KMeansFuncGenerator( )
    maxTemperObj <- findMaxTemper(nIters = 300,
                                  statsFuncList = KMeansObj$statsFuncList,
                                  temperLadder = KMeansObj$temperLadder,
                                  startingVals = KMeansObj$startingVals,
                                  logTarDensFunc = KMeansObj$logTarDensFunc,
                                  levelsSaveSampFor = 1:4,
                                  guideMe = FALSE,
                                  verboseLevel = 1)
    print(maxTemperObj)    
}


examplesPlaceTempers <-
    function ( )
{
    if ('package:EMCC' %in% search( ))
        detach('package:EMCC')
    library(EMCC, lib.loc='~/Rlib')

    ## The KMeans data set example 
    KMeansFuncGenerator <-
        function (seed = 13579, plotIt = TRUE)
        {
            set.seed(seed)                
            data('KMeansData')
            origLabels <- KMeansData[ , 1]
            yy <- KMeansData[ , -1]

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
                    priorMinClusters <- 3; priorMaxClusters <- 3
                    sum <- -sum +
                        ifelse(((nClusters < priorMinClusters) ||
                                (nClusters > priorMaxClusters)),
                               -Inf, -log(priorMaxClusters - priorMinClusters + 1))
                    return(sum)
                }            
            sampDim <- length(origLabels)
            ss <- c(0:2, rep(0, sampDim - 3))
            return(list(startingVal = ss,
                        origLabels = origLabels,
                        yy = yy,
                        logTarDensFunc = logTarDensFunc))
        }

    KMeansObj <- KMeansFuncGenerator( )
    placeTempersObj <- placeTempers(nIters = 300,
                                    temperLimits = c(0.5, 10),
                                    acceptRatioLimits = c(0.3, 0.5),
                                    ladderLenMax = 50,
                                    startingVals = KMeansObj$startingVal,
                                    logTarDensFunc = KMeansObj$logTarDensFunc,
                                    scheme = "log",
                                    ladderLen = 10,
                                    verboseLevel = 1)
    print(placeTempersObj)    
}






