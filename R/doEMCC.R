
###  $Id: doEMCC.R,v 1.11 2008/02/06 05:11:09 goswami Exp $
###  
###     File: doEMCC.R
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


### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### The following is the main 'workhorse' function.

TOEMCCMain <-
    function (nIters,              
              temperLadder,
              startingVals,
              logTarDensFunc,
              MHBlocks                 = NULL,
              MHBlockNTimes            = NULL,
              MHMoveType               = NULL,
              MHMergeProb              = 0.5,
              moveProbsList            = NULL,
              moveNTimesList           = NULL,
              moveSelectionCodesList   = NULL,
              moveSelectionTempersList = NULL,
              levelsSaveSampFor        = NULL,
              saveFitness              = FALSE,
              saveAcceptRatiosList     = FALSE,             
              timeInSecs               = -1,
              verboseLevel             = 0,
              ...)   
{
    ptm <- proc.time( )
    ## BEGIN: Error checks
    nIters       <- EMC:::.check.numericWithLLim(nIters, 0)
    temperLadder <- EMC:::.check.temperLadder(temperLadder)
    nLevels      <- as.integer(length(temperLadder))
    startingVals <- EMC:::.check.startingVals(startingVals, nLevels)
    startingVals <- apply(startingVals, c(1, 2), as.integer)
    sampDim      <- as.integer(ncol(startingVals))
    
    logTarDensFunc <- EMC:::.check.logTarDensFunc(logTarDensFunc)

    MHBlocks      <- EMC:::.check.MHBlocks(MHBlocks, sampDim)
    MHNBlocks     <- as.integer(length(MHBlocks))
    MHBlockNTimes <- EMC:::.check.MHBlockNTimes(MHBlockNTimes, MHNBlocks)

    moveProbsList  <- .check.moveProbsList(moveProbsList, nLevels)
    moveNTimesList <- .check.moveNTimesList(moveNTimesList,
                                            moveProbsList = moveProbsList,
                                            nLevels       = nLevels)
    
    MHMoveType    <- .check.MHMoveType(MHMoveType)
    MHMergeProb   <- .check.MHMergeProb(MHMergeProb)

    moveSelectionCodesList   <- .check.mSCList(moveSelectionCodesList,
                                               moveNTimesList = moveNTimesList)
    moveSelectionTempersList <- .check.mSTList(moveSelectionTempersList,
                                               moveNTimesList = moveNTimesList,
                                               temperColdest  = temperLadder[nLevels])
    
    levelsSaveSampFor <- EMC:::.check.levelsSaveSamplesFor(levelsSaveSampFor, nLevels)
    saveFitness       <- EMC:::.check.logical(saveFitness)    
    timeInSecs        <- EMC:::.check.timeInSecs(timeInSecs)      
    verboseLevel      <- EMC:::.check.numericWithLLim(verboseLevel, NA)  
    procTimeFunc      <- as.function(proc.time)             
    procTimeFuncEnv   <- new.env( )
    doCallFunc        <- as.function(EMC:::doCall)             
    doCallFuncEnv     <- new.env( )
    dotsList          <- list(...)
    argsList          <- EMC:::collectVarnames(ls( ))
    ## E N D: Error checks

    if (argsList$verboseLevel >= 3) {
        cat('The processed arguments\n')
        print(argsList)
    }
    if (argsList$verboseLevel >= 1) cat('\nBEGIN: EMCC\n')        
    res <- .Call('EMCCMainC', argsList)
    if (argsList$verboseLevel >= 1) cat('E N D: EMCC\n')

    res$nIters               <- nIters
    res$levelsSaveSampFor    <- levelsSaveSampFor
    res$temperLadder         <- temperLadder
    res$startingVals         <- startingVals
    res$moveProbsList        <- moveProbsList
    res$moveNTimesList       <- moveNTimesList
    res$detailedAcceptRatios <- EMC:::.get.detailedAcceptRatios(res$acceptRatiosList)

    if (!saveAcceptRatiosList)
        res$acceptRatiosList <- NULL

    if (!any(is.na(ptm))) res$time <- proc.time( ) - ptm
    class(res) <- 'EMCC'
    res
}

### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

evolMonteCarloClustering <-
  function (nIters,              
            temperLadder,
            startingVals,
            logTarDensFunc,
            MHMergeProb       = 0.5,
            moveProbsList     = NULL,
            moveNTimesList    = NULL,
            levelsSaveSampFor = NULL,
            saveFitness       = FALSE,
            verboseLevel      = 0,
            ...)              
{    
    if (is.null(moveNTimesList)) {
        temperLadder   <- EMC:::.check.temperLadder(temperLadder)
        nLevels        <- length(temperLadder)
        moveNTimesList <- list(MH = 1,
                               RE = nLevels)
    }

    ## THINKME: Below we are passing MHBlocks = NULL to TOEMCCMain, so
    ## that each dimension (i.e., each object to be clustered) is a
    ## block of its own. Is there way to implement a generalized
    ## blocking scheme in this context?
    res <- TOEMCCMain(nIters            = nIters,              
                      temperLadder      = temperLadder,
                      startingVals      = startingVals,
                      logTarDensFunc    = logTarDensFunc,
                      MHBlocks          = NULL,
                      MHBlockNTimes     = NULL,
                      MHMergeProb       = MHMergeProb,
                      moveProbsList     = moveProbsList,
                      moveNTimesList    = moveNTimesList,
                      levelsSaveSampFor = levelsSaveSampFor,
                      saveFitness       = saveFitness,
                      verboseLevel      = verboseLevel,
                      ...)
    res
}

### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

print.EMCC <-
    function (x, ...)
{
    if (!is.null(x$temperLadder)) {
        cat('\nThe temperature ladder:\n')
        print(x$temperLadder, ...)
    }
    cat('\nThe overall acceptance rate summary:\n')
    print(x$acceptRatios, ...)
    cat('\n')    
}







