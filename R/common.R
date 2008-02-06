
###  $Id: common.R,v 1.1 2008/02/05 19:30:08 goswami Exp $
###  
###  File:    common.R
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


### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### The following functions are either new to this package or are
### re-definitions (due to relevance issues) of their same-named
### counter parts to be found in the package EMC

### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Some utility functions 

fatal <-
    function (msg, var, formatMsg = TRUE,
              stopMsg = 'Fix the above problem first')
{
    if (formatMsg)
        cat(strsplit(msg, ' ')[[1]], '\n', fill = TRUE)
    else
        cat(msg, '\n')
    if (!missing(var))
        print(var)
    stop(stopMsg, call. = FALSE)
}

### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.check.MHMoveType <-
    function (MHMoveType)
{
    if (is.null(MHMoveType))
        return(as.character('SPLIT_MERGE'))

    MHMoveType <- toupper(MHMoveType)
    if (MHMoveType %in% c('SPLIT_MERGE'))
        return(as.character(MHMoveType))

    if (MHMoveType %in% c('GIBBS'))
        msg <- 'not implemented yet'
      else
        msg <- paste('Please provide a valid MHMoveType :: it should be in',
                     '("SPLIT_MERGE", "GIBBS")')
    
    fatal(msg)
}

.check.MHMergeProb <-
    function (MHMergeProb)
{
    if (!is.numeric(MHMergeProb) ||
        !((0 < MHMergeProb) && (MHMergeProb < 1.0))) {
        msg <- paste('Please provide a valid MHMergeProb :: it should be',
                     'a probability, given object:')
        fatal(msg, MHMergeProb)
    }

    as.double(MHMergeProb)
}

### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.get.implementedMovenames <-
    function ( )
    c('MH', 'RC', 'SCSC_ONE_NEW', 'SCSC_TWO_NEW', 'SCRC',
      'RE', 'BCE', 'BIRE', 'BSE', 'CE') 

.get.specialMovenames <-
    function ( )
    c('BCE', 'BIRE', 'BSE', 'CE') 

.get.implementedMovesList <-
    function (func, implementedMovenames = .get.implementedMovenames( ))
    sapply(implementedMovenames, function (cc) func(0), simplify = FALSE)

.check.moveList.do <-
    function (aList,
              listname             = substitute(aList),
              implementedMovenames = .get.implementedMovenames( ))
{
    msg <- paste('Please provide a list for', listname)

    if (!is.list(aList)) {
        msg <- paste(msg, ', given object:', sep = '')
        fatal(msg, aList)
    }

    given <- names(aList)
    if (!('MH' %in% given)) {
        msg <- paste(msg, 'having a component called MH, given list:')
        fatal(msg, aList)
    }

    if (length(notKnown <- setdiff(given, implementedMovenames)) > 0) {
        msg <- paste('Some of the components of the list', listname, 'are not',
                     'implemented move names, allowable names are:')
        fatal(msg, implementedMovenames)
    }

    if (any(unlist(aList) < 0)) {
        msg <- paste('All the components of', listname, 'are not positive,',
                     'given list:')
        fatal(msg, aList)
    }

    if (aList$MH <= 0) {
        msg <- paste('The "MH" component of', listname, 'should be positive,',
                     'given list:')
        fatal(msg, aList)
    }
}

.check.moveList.fill <-
    function (aList, func)
{
    ll <- .get.implementedMovesList(func)
    for (nn in names(aList))
        ll[[nn]] <- (func)(aList[[nn]])

    ll
}    

.check.moveProbsList <-
    function (moveProbsList, nLevels)
{
    if (is.null(moveProbsList)) {
        if (nLevels == 1) {
            moveProbsList <- list(MH = as.double(1))
        }
        else {
            moveProbsList <- list(MH             = as.double(0.5),
                                  RC             = as.double(0.25),
                                  'SCSC_TWO_NEW' = as.double(0.25))
        }
        return(.check.moveList.fill(moveProbsList, as.double))
    }

    .check.moveList.do(moveProbsList)
    if (abs(sum(unlist(moveProbsList)) - 1) > (getOption('ts.eps') / 2)) {
        msg <- paste('The probabilities should add up to 1 in moveProbsList,',
                     'given list:')
        fatal(msg, moveProbsList)
    }

    .check.moveList.fill(moveProbsList, as.double)
}

.check.moveNTimesList <-
    function (moveNTimesList, moveProbsList, nLevels)
{
    nmpl     <- names(moveProbsList[moveProbsList > 0])
    defaults <- list(MH             = as.integer(1),
                     RC             = as.integer(floor(nLevels / 2)),
                     'SCSC_TWO_NEW' = as.integer(floor(nLevels / 2)),
                     'SCSC_ONE_NEW' = as.integer(floor(nLevels / 2)),
                     SCRC           = as.integer(floor(nLevels / 2)),
                     RE             = as.integer(nLevels),
                     BCE            = as.integer(max(1, nLevels - 3)),
                     BIRE           = as.integer(max(1, nLevels - 3)),
                     BSE            = as.integer(max(1, nLevels - 3)),
                     CE             = as.integer(max(1, nLevels - 3)))
    
    if (is.null(moveNTimesList)) {
        if (nLevels == 1) {
            moveNTimesList <- list(MH = as.integer(1))
        }
        else {
            moveNTimesList <- defaults[nmpl] 
        }
        return(.check.moveList.fill(moveNTimesList, as.integer))
    }

    .check.moveList.do(moveNTimesList)
    ll <- .check.moveList.fill(moveNTimesList, as.integer)

    nmntl <- names(moveNTimesList[moveNTimesList > 0])    
    if (length(forgot <- setdiff(nmpl, nmntl)) > 0) {
        warning('Some of the components (',  toString(forgot), ') of ',
                'moveNTimesList were made 1 since the corresponding ',
                'components of moveProbsList were positive', call. = FALSE)
        for (nn in forgot)
            ll[[nn]] <- as.integer(1)
    }
    
    for (nn in intersect(names(moveNTimesList), .get.specialMovenames( ))) {
        if (moveNTimesList[[nn]] > nLevels - 2) {
            msg <- paste('Please provide a valid moveNTimes for move ', nn,
                         ' in moveNTimesList :: it should have non-negative',
                         ' number <= ', nLevels - 2, ' since nLevels = ',
                         nLevels, '. The given moveNTimes:', sep = '')
            fatal(msg, moveNTimesList[[nn]])
        }
    }
    
    ll
}

### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.get.selectionCodes <-
    function ( )
    c('random', 'best', 'worst')

.get.defaults.selectionCodes <-
    function ( )
{
    list(RC             = c('best', 'best'),
         'SCSC_TWO_NEW' = c('best', 'best'),
         'SCSC_ONE_NEW' = c('best', 'best'),
         SCRC           = c('random'))
}

.check.mSCList.do <-
    function (moveSelectionCodesList, moveNTimesList, movename, nVals,
              selectionCodes = .get.selectionCodes( ))
{
    vv   <- moveSelectionCodesList[[movename]]
    nn   <- moveNTimesList[[movename]]
    msg1 <- paste('The', movename, 'component of moveSelectionCodesList ')
    
    if (is.null(vv)) {
        if (nn > 0) {
            msg <- paste(msg1, 'is not provided but the corresponding component ',
                         'of moveNTimesList is positive: ', nn)
            fatal(msg)
        }

        return(invisible(0))
    }
    
    if (length(vv) != nVals) {
        msg <- paste(msg1, ' should be of length ', nVals, ', given vector:',
                     sep = '')
        fatal(msg, vv)
    }

    if (length(invalid <- setdiff(tolower(vv), selectionCodes)) > 0) {
        msg <- paste(msg1, ' has some invalid selection codes: ',
                     toString(invalid), ', allowed codes are:', sep = '')
        fatal(msg, selectionCodes)
    }
}

.check.mSCList <-
    function (moveSelectionCodesList, moveNTimesList)
{
    defaults <- .get.defaults.selectionCodes( )

    if (is.null(moveSelectionCodesList))
        return(sapply(defaults, as.character, simplify = FALSE))

    if (!is.list(moveSelectionCodesList)) {
        msg <- paste('Please provide a list for moveSelectionCodesList,',
                     'given object:')
        fatal(msg, moveSelectionCodesList)
    }    

    for (vv in names(defaults)) {
        .check.mSCList.do(moveSelectionCodesList = moveSelectionCodesList,
                          moveNTimesList         = moveNTimesList,
                          movename               = vv,
                          nVals                  = length(defaults[[vv]]))
    }
    
    needComps  <- names(defaults)
    nmscl      <- names(moveSelectionCodesList)
    
    if (length(extraComps <- setdiff(nmscl, needComps)) > 0) 
        warning('The following extra components of moveSelectionCodesList ',
                'will be ignored: ', toString(extraComps), call. = FALSE)

    missingComps           <- setdiff(needComps, nmscl)
    moveSelectionCodesList <- c(moveSelectionCodesList, defaults[missingComps])
    sapply(moveSelectionCodesList[needComps], as.character, simplify = FALSE)
}

### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.get.defaults.selectionTempers <-
    function (temperColdest)
{
    list(RC             = rep(temperColdest, 2),
         'SCSC_TWO_NEW' = rep(temperColdest, 2),
         'SCSC_ONE_NEW' = rep(temperColdest, 2),
         SCRC           = temperColdest,
         BCE            = temperColdest)
}

.check.mSTList.do <-
    function (moveSelectionTempersList, moveNTimesList, movename, nVals)
{
    vv   <- moveSelectionTempersList[[movename]]
    nn   <- moveNTimesList[[movename]]
    msg1 <- paste('The', movename, 'component of moveSelectionTempersList')

    if (is.null(vv)) {
        if (nn > 0) {
            msg <- paste(msg1, 'is not provided but the corresponding component ',
                         'of moveNTimesList is positive: ', nn)
            fatal(msg)
        }
        
        return(invisible(0))
    }
    
    if (length(vv) != nVals) {
        msg <- paste(msg1, ' should be of length ', nVals, ', given vector:',
                     sep = '')
        fatal(msg, vv)
    }
    
    if (any(bad <- (vv < 0))) {
        msg <- paste(msg1, 'has some non-negative selectionTempers:')
        fatal(msg, vv[bad])
    }
}

.check.mSTList <-
    function (moveSelectionTempersList, moveNTimesList, temperColdest)
{
    defaults <- .get.defaults.selectionTempers(temperColdest)
    
    if (is.null(moveSelectionTempersList))
        return(sapply(defaults, as.double, simplify = FALSE))

    if (!is.list(moveSelectionTempersList)) {
        msg <- paste('Please provide a list for moveSelectionTempersList,',
                     'given object:')
        fatal(msg, moveSelectionTempersList)
    }

    for (vv in names(defaults)) {
        .check.mSTList.do(moveSelectionTempersList = moveSelectionTempersList,
                          moveNTimesList           = moveNTimesList,
                          movename                 = vv,
                          nVals                    = length(defaults[[vv]]))
    }    

    needComps <- names(defaults)
    nmstl     <- names(moveSelectionTempersList)
    
    if (length(extraComps <- setdiff(nmstl, needComps)) > 0)
        warning('The following extra components of moveSelectionTempersList ',
                'will be ignored: ', toString(extraComps), call. = FALSE)

    missingComps             <- setdiff(needComps, nmstl)
    moveSelectionTempersList <- c(moveSelectionTempersList, defaults[missingComps])
    sapply(moveSelectionTempersList[needComps], as.double, simplify = FALSE)
}





