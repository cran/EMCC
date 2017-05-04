
###  $Id: funcGenerators.R,v 1.2 2008/02/05 19:30:08 goswami Exp $
###  
###     File: funcGenerators.R
###  Package: EMCC
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


.genMixNormalData <-
    function (nObs,
              knownClusterMeans,
              SigmaArr,
              mixProps = NULL)
    ## Note:
    ##
    ## the mean vectors of the mixture components are in the rows of
    ## knownClusterMeans
    ##
    ## the dispersion matrices of the mixture components are collected
    ## in the array SigmaArr
{
    nMixComps <- nrow(knownClusterMeans)
    mm        <- ncol(knownClusterMeans)
    if (is.null(mixProps))
        mixProps <- rep(1, nMixComps) / nMixComps
    
    ret       <- matrix(0, nObs, mm + 1)
    ret[ , 1] <- sample(seq_len(nMixComps),
                        size    = nObs,
                        prob    = mixProps,
                        replace = TRUE)
    for (ii in seq_len(nObs)) {
        comp                        <- ret[ii, 1]
        ret[ii, seq.int(2, mm + 1)] <- mvrnorm(n     = 1,
                                               mu    = knownClusterMeans[comp, ],
                                               Sigma = SigmaArr[ , , comp])
    }
    
    ret
}

partitionRep <-
    function (clusterInd)
{
    stopifnot(is.vector(clusterInd))
    nn <- length(clusterInd)
    stopifnot(all(clusterInd >= 0))
    stopifnot(all(clusterInd < nn))

    uu <- unique(clusterInd)
    id <- seq_len(nn)
    ll <- sapply(seq_along(uu), function (ii)
             {
                 list(id[clusterInd == uu[ii]])
             })
    list(clusterLabels = uu,
         clusters      = ll)
}

clusterPlot <-
    function (clusterInd,
              data,
              main              = '',
              sub               = '',
              knownClusterMeans = NULL,
              ...)                   
{
    uu <- unique(clusterInd)
    plot(data,
         type     = 'n',
         main     = as.expression(main),
         cex.main = 0.8,
         sub      = as.expression(sub),
         cex.sub  = 0.8,
         xlab     = as.expression(substitute(x[ll], list(ll = 1))),
         ylab     = as.expression(substitute(x[ll], list(ll = 2))),
         ...)
    for (ii in seq_along(uu)) {
        points(data[clusterInd == uu[ii], , drop = FALSE],
               pch = as.character(ii),
               col = 1 + ii) 
    }
    if (!is.null(knownClusterMeans))
        points(knownClusterMeans, col = 1, pch = 16)
}

KMeansFuncGenerator1 <-
    function (seed, plotIt = TRUE)
{
    set.seed(seed)    
    mm                     <- 2
    nMixComps              <- 2
    knownClusterMeans      <- matrix(0, nMixComps, mm)
    knownClusterMeans[1, ] <- c(0, 0)
    knownClusterMeans[2, ] <- c(1, 1)
    nObs                   <- 40
    sigma                  <- sqrt(2) / 8
    SigmaArr               <- array(dim = c(mm, mm, nMixComps))
    for (ii in seq_len(nMixComps))
        SigmaArr[ , , ii] <- diag(rep(sigma, mm))
    
    data       <- .genMixNormalData(nObs              = nObs,
                                    knownClusterMeans = knownClusterMeans,
                                    SigmaArr          = SigmaArr)
    origLabels <- data[ , 1]
    yy         <- data[ , -1]

    if (plotIt) {
        par(mfcol = c(2, 2))
        clusterPlot(clusterInd        = origLabels,
                    data              = yy,
                    main              = 'Original data',                    
                    knownClusterMeans = knownClusterMeans)
        sub <- paste('# of clusters = ', nMixComps, ', assumed known', sep = '')
        kmObj <- kmeans(yy,
                        centers  = nMixComps,
                        iter.max = 1000,
                        nstart   = 1000)
        cat('\nThe K-Means clustering summary:\n')
        print(kmObj)
        clusterPlot(clusterInd        = kmObj$cluster,
                    data              = yy,
                    main              = 'K-Means clusterting',
                    sub               = sub,
                    knownClusterMeans = knownClusterMeans)
        mcObj <- Mclust(yy, G = nMixComps)
        cat('\nThe MCLUST clustering summary:\n')
        print(mcObj)
        clusterPlot(clusterInd        = mcObj$classification,
                    data              = yy,
                    main              = 'MCLUST clusterting',
                    sub               = sub, 
                    knownClusterMeans = knownClusterMeans)
    }

    priorMinClusters <- 1
    priorMaxClusters <- 2
    priorNClusters   <-
        function (pr)
        {
            nClusters <- length(pr$clusterLabels)
            ifelse(((nClusters < priorMinClusters) ||
                    (nClusters > priorMaxClusters)),
                   -Inf,
                   -log(priorMaxClusters - priorMinClusters + 1))
        }
    
    logTarDensFunc <-
        function (draw, ...)
        {
            pr <- partitionRep(as.integer(draw))
            ss <- 0
            for (cluster in pr$clusters) {
                YY        <- yy[cluster, , drop = FALSE]
                YYColMean <- apply(YY, 2, mean)
                YYSwept   <- sweep(YY, 2, YYColMean)
                ss        <- ss + sum(YYSwept^2)
            }
##             print(pr$clusters)
##             print(-ss + priorNClusters(pr))
            -ss + priorNClusters(pr)
        }
    
    list(knownClusterMeans = knownClusterMeans,
         yy                = yy,
         origLabels        = origLabels,
         priorMinClusters  = priorMinClusters,
         priorMaxClusters  = priorMaxClusters,
         logTarDensFunc    = logTarDensFunc)
}

KMeansFuncGenerator2 <-
    function (seed, plotIt = TRUE)
{
    set.seed(seed)    
    mm            <- 2
    nMixComps     <- 4
    knownClusterMeans      <- matrix(0, nMixComps, mm)
    knownClusterMeans[1, ] <- c(-3, 0.75)
    knownClusterMeans[2, ] <- c(-3, -0.75)
    knownClusterMeans[3, ] <- c(3, 1)
    knownClusterMeans[4, ] <- c(3, -1)
    sigmaSq       <- c(1, 1, 0.2, 0.2)
    rho           <- c(0.95, 0.95, 0, 0)
    ## mixProps      <- c(0.3, 0.3, 0.2, 0.2)
    mixProps      <- NULL
    nObs          <- 60

    ARDisp <-
        function (rho)
        {
            tmp <- rep(1, mm)
            diag((1 - rho) * tmp) + rho * tmp %*% t(tmp)
        }

    SigmaArr <- array(dim = c(mm, mm, nMixComps))
    for (ii in seq_len(nMixComps))
        SigmaArr[ , , ii] <- sigmaSq[ii] * ARDisp(rho[ii])
    
    data       <- .genMixNormalData(nObs              = nObs,
                                    knownClusterMeans = knownClusterMeans,
                                    SigmaArr          = SigmaArr,
                                    mixProps          = mixProps)
    origLabels <- data[ , 1]
    yy         <- data[ , -1]

    if (plotIt) {
        par(mfcol = c(2, 2))
        clusterPlot(clusterInd        = origLabels,
                    data              = yy,
                    main              = 'Original data',                    
                    knownClusterMeans = knownClusterMeans)
        sub <- paste('# of clusters = ', nMixComps, ', assumed known', sep = '')
        kmObj <- kmeans(yy,
                        centers  = nMixComps,
                        iter.max = 10000,
                        nstart   = 10000)
        cat('\nThe K-Means clustering summary:\n')
        print(kmObj)
        clusterPlot(clusterInd        = kmObj$cluster,
                    data              = yy,
                    main              = 'K-Means clusterting',
                    sub               = sub,
                    knownClusterMeans = knownClusterMeans)
        mcObj <- Mclust(yy, G = nMixComps)
        cat('\nThe MCLUST clustering summary:\n')
        print(mcObj)
        clusterPlot(clusterInd        = mcObj$classification,
                    data              = yy,
                    main              = 'MCLUST clusterting',
                    sub               = sub, 
                    knownClusterMeans = knownClusterMeans)
    }

    priorMinClusters <- 1
    priorMaxClusters <- 5
    priorNClusters   <-
        function (pr)
        {
            nClusters <- length(pr$clusterLabels)
            ifelse(((nClusters < priorMinClusters) ||
                    (nClusters > priorMaxClusters)),
                   -Inf,
                   -log(priorMaxClusters - priorMinClusters + 1))
        }
    
    ## See section 5.4 of Goswami, Liu and Wong (2007) for details on
    ## the following the density
    nu0                <- 6
    kappa0             <- 0.05
    Lambda0            <- diag(mm)
    logDetLambda0Const <- nu0 * log(det(Lambda0)) / 2
    mu0                <- numeric(mm)
    const1             <- -mm * log(pi) / 2
    logTarDensFunc     <-
        function (draw, ...)
        {
            pr <- partitionRep(as.integer(draw))
            ss <- 0
            for (cluster in pr$clusters) {
                nn     <- length(cluster)
                nuN    <- nu0 + nn
                kappaN <- kappa0 + nn
                tmp1   <- kappa0 / kappaN
                
                ss <- ss + nn * const1 + mm * log(tmp1) / 2
                rr <- seq_len(mm)
                ss <- ss + sum(lgamma((nuN + 1 - rr) / 2) -
                               lgamma((nu0 + 1 - rr) / 2))

                YY        <- yy[cluster, , drop = FALSE]
                YYColMean <- apply(YY, 2, mean)
                YYSwept   <- sweep(YY, 2, YYColMean)
                tmp2      <- YYColMean - mu0
                LambdaN   <- (Lambda0 + (t(YYSwept) %*% YYSwept) +
                              (nn * tmp1) * (tmp2 %*% t(tmp2)))
                
                ss <- ss + logDetLambda0Const - (nuN * log(det(LambdaN)) / 2)
            }

            ss + priorNClusters(pr)
        }
    
    list(knownClusterMeans = knownClusterMeans,
         yy                = yy,
         origLabels        = origLabels,
         priorMinClusters  = priorMinClusters,
         priorMaxClusters  = priorMaxClusters,
         logTarDensFunc    = logTarDensFunc)
}

