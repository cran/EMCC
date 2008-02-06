
###  $Id: ladder.R,v 1.3 2008/02/03 04:18:54 goswami Exp $
###  
###     File: ladder.R
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


### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### The following is the rewrite of the package.

findMaxTemper <-
    function (nIters,
              statsFuncList,
              startingVals,
              logTarDensFunc,
              temperLadder      = NULL,
              temperLimits      = NULL,
              ladderLen         = 10,
              scheme            = 'exponential',
              schemeParam       = 0.5,
              cutoffDStats      = 1.96,
              cutoffESS         = 50,
              guideMe           = TRUE,
              levelsSaveSampFor = NULL,
              saveFitness       = FALSE,
              doFullAnal        = TRUE,
              verboseLevel      = 0,
              ...)   
{
    getEMCOutFunc <-
        function (ladder, ...)
        {
            EMCCOut <- TOEMCCMain(nIters            = nIters,
                                  temperLadder      = ladder,
                                  startingVals      = startingVals,
                                  logTarDensFunc    = logTarDensFunc,
                                  moveProbsList     = list(MH = 1.0),
                                  moveNTimesList    = list(MH = 1),
                                  levelsSaveSampFor = seq_along(ladder),
                                  saveFitness       = TRUE,
                                  verboseLevel      = verboseLevel,
                                  ...)
            
            if (verboseLevel >= 1) {
                cat('\nThe detailed MH acceptance ratios:\n')
                print(EMCCOut$detailedAcceptRatios$MH)
            }
            
            EMCCOut
        }

    ret <- EMC::findMaxTemper(nIters            = nIters,
                              statsFuncList     = statsFuncList,
                              startingVals      = startingVals,
                              logTarDensFunc    = logTarDensFunc,
                              MHPropNewFunc     = NULL,
                              logMHPropDensFunc = NULL,
                              temperLadder      = temperLadder,
                              temperLimits      = temperLimits,
                              ladderLen         = ladderLen,
                              scheme            = scheme,
                              schemeParam       = schemeParam,
                              cutoffDStats      = cutoffDStats,
                              cutoffESS         = cutoffESS,
                              guideMe           = guideMe,
                              levelsSaveSampFor = levelsSaveSampFor,
                              saveFitness       = saveFitness,
                              doFullAnal        = doFullAnal,
                              verboseLevel      = verboseLevel,
                              getEMCOutFunc     = getEMCOutFunc,
                              ...)   

    class(ret) <- 'EMCCMaxTemper'
    ret
}

### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

print.EMCCMaxTemper <-
    function (x, ...)
{
    EMC:::print.EMCMaxTemper(x, ...)
}

### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

placeTempers <-
    function (nIters,
              acceptRatioLimits,
              ladderLenMax,              
              startingVals,
              logTarDensFunc,
              temperLadder      = NULL,
              temperLimits      = NULL,
              ladderLen         = 15,
              scheme            = 'exponential',
              schemeParam       = 1.5,
              guideMe           = TRUE,
              levelsSaveSampFor = NULL,
              saveFitness       = FALSE,
              verboseLevel      = 0,
              ...)                 
{
    getEMCOutFunc <-
        function (ladder, ...)
        {
            nLevels <- length(ladder)            
            TOEMCCMain(nIters            = nIters,
                       temperLadder      = ladder,
                       startingVals      = startingVals,
                       logTarDensFunc    = logTarDensFunc,
                       moveProbsList     = list(MH = 0.3, RE = 0.7),
                       moveNTimesList    = list(MH = 1, RE = nLevels),
                       levelsSaveSampFor = seq_len(nLevels),
                       saveFitness       = TRUE,
                       verboseLevel      = verboseLevel,
                       ...)                        
        }

    ret <- EMC::placeTempers(nIters            = nIters,
                             acceptRatioLimits = acceptRatioLimits,
                             ladderLenMax      = ladderLenMax,
                             startingVals      = startingVals,
                             logTarDensFunc    = logTarDensFunc,
                             MHPropNewFunc     = NULL,
                             logMHPropDensFunc = NULL,
                             temperLadder      = temperLadder,
                             temperLimits      = temperLimits,
                             ladderLen         = ladderLen,
                             scheme            = scheme,
                             schemeParam       = schemeParam,
                             guideMe           = guideMe,
                             levelsSaveSampFor = levelsSaveSampFor,
                             saveFitness       = saveFitness,
                             verboseLevel      = verboseLevel,
                             getEMCOutFunc     = getEMCOutFunc,
                             ...)

    class(ret) <- 'EMCCPlaceTempers'
    ret
}

### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

print.EMCCPlaceTempers <-
    function (x, ...)
{
    EMC:::print.EMCPlaceTempers(x, ...)
}





