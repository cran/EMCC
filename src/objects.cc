
/*
 *  $Id: objects.cc,v 1.5 2010/04/18 17:04:09 goswami Exp $
 *  
 *  File: objects.C
 *  Copyright (C) 2006-present Gopi Goswami 
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  For a copy of the GNU General Public License please write to the
 *  Free Software Foundation, Inc.
 *  59 Temple Place, Suite 330.
 *  Boston, MA  02111-1307 USA.
 *
 *  For bugs in the code please contact:
 *  <goswami@stat.harvard.edu>
 *
 *  SYNOPSIS
 *
 *
 *
 *  DESCRIPTION
 *
 *
 *
 */ 

#include <iostream>
#include "objects.h"

#ifndef DEBUG_OBJECTS
#define DEBUG_OBJECTS 0
#endif 

#if (!DEBUG_OBJECTS)
#define DEBUG(_x) ((void) 0)
#endif

/*****************************************************************************/
// The public functions.
/*****************************************************************************/

Sampler::Sampler (SEXP opts)
{
        SEXP SEXPTmp, SEXPCodes;
        int ii, jj, count, *intsTmp, *startingVals;
        double *doublesTmp;
        char *name;

        nIters_ = INTEGER(getListElement(opts, "nIters"))[0];
        timeInSecs_ = REAL(getListElement(opts, "timeInSecs"))[0];
        if (timeInSecs_ <= 0) nItersActual_ = nIters_;
        else                 nItersActual_ = -1;
        
        verboseLevel_ = INTEGER(getListElement(opts, "verboseLevel"))[0];
        if (nIters_ <= 1000) printEstTimeAt_ = 10;
        else                 printEstTimeAt_ = 100;
        printEstTimeNTimes_ = 10;
        printInitialDotsWhen_ = printEstTimeAt_ / 10;
        printDotAt_ = 0; nDotsPerLine_ = 20;
        eachDotWorth_ = (int) ceil((nIters_ - printEstTimeAt_ + 1.0) / \
                                  (printEstTimeNTimes_ * nDotsPerLine_));
        
        SEXPTmp = getListElement(opts, "temperLadder");
        nLevels_ = length(SEXPTmp);
        doublesTmp = REAL(SEXPTmp);
        temperLadder_ = (double *) R_alloc(nLevels_, sizeof(double));
        invTemperLadder_ = (double *) R_alloc(nLevels_, sizeof(double));
        for (ii = 0; ii < nLevels_; ++ii) {
                temperLadder_[ii] = doublesTmp[ii];
                invTemperLadder_[ii] = 1.0 / doublesTmp[ii];
        }
        logDensities_ = (double *) R_alloc(nLevels_, sizeof(double));
        logDensitiesStore_ = (double *) R_alloc(nLevels_, sizeof(double));
        scratch_SLC_ = (SampleLevelContext *) R_alloc(1, sizeof(struct SampleLevelContext));
        scratch_SLC_->logWeights_ = (double *) R_alloc(nLevels_, sizeof(double));
        scratch_SLC_->adjWeights_ = (double *) R_alloc(nLevels_, sizeof(double));
        scratch_SLC_->partialSum_ = (double *) R_alloc(nLevels_, sizeof(double));
        
        SEXPTmp = getListElement(opts, "startingVals");
        startingVals = INTEGER(SEXPTmp);        
        sampDim_ = ncols(SEXPTmp);
        SEXPTmp = getListElement(opts, "MHBlocks");
        MHNBlocks_ = length(SEXPTmp);
        
        SEXPTmp = getListElement(opts, "MHMoveType");
        if (strcmp(CHAR(STRING_ELT(SEXPTmp, 0)), "SPLIT_MERGE") == 0)
                MHMoveCode_ = SPLIT_MERGE;
        else if (strcmp(CHAR(STRING_ELT(SEXPTmp, 0)), "GIBBS") == 0)
                MHMoveCode_ = GIBBS;
        
        MHMergeProb_ = REAL(getListElement(opts, "MHMergeProb"))[0];

        MH_OLC_ = (ObjLabelContext *) R_alloc(1, sizeof(struct ObjLabelContext));
        boolMat_ = (bool **) R_alloc(sampDim_, sizeof(bool *));
        for (ii = 0; ii < sampDim_; ++ii)
                boolMat_[ii] = (bool *) R_alloc(sampDim_, sizeof(bool));
        for (ii = 0; ii < N_SUB_CLUSTERS; ++ii)
                SCXX_intersection_[ii] = (int *) R_alloc(sampDim_, sizeof(int));
        nSCXX_nonEmptyInRow_ = (int *) R_alloc(sampDim_, sizeof(int));
        SCXX_nonEmptyInRow_ = (int **) R_alloc(sampDim_, sizeof(int *));
        for (ii = 0; ii < sampDim_; ++ii)
                SCXX_nonEmptyInRow_[ii] = (int *) R_alloc(sampDim_, sizeof(int));
        SCXX_rows2OrMore_ = (int *) R_alloc(sampDim_, sizeof(int));
        labelDiffs_ = (int *) R_alloc(sampDim_, sizeof(int));
        labelDiffsPositive_ = (int *) R_alloc(sampDim_, sizeof(int));
        labelDiffsNegative_ = (int *) R_alloc(sampDim_, sizeof(int));
        SCRC_isAlreadySampled_ = (bool *) R_alloc(sampDim_, sizeof(bool));
        SCRC_union_ = (int *) R_alloc(sampDim_, sizeof(int));
        SCRC_subSampleIndices_ = (int *) R_alloc(sampDim_, sizeof(int));

        for (ii = 0; ii < N_IMPLEMENTED_MOVES; ++ii)
                moveObjs_[ii] = (MoveObject *) R_alloc(1, sizeof(struct MoveObject));
        strcpy(moveObjs_[MH]->name_, "MH");
        strcpy(moveObjs_[RC]->name_, "RC");
        strcpy(moveObjs_[SCSC_ONE_NEW]->name_, "SCSC_ONE_NEW");
        strcpy(moveObjs_[SCSC_TWO_NEW]->name_, "SCSC_TWO_NEW");
        strcpy(moveObjs_[SCRC]->name_, "SCRC");
        strcpy(moveObjs_[RE]->name_, "RE");
        strcpy(moveObjs_[BCE]->name_, "BCE");
        strcpy(moveObjs_[BIRE]->name_, "BIRE");
        strcpy(moveObjs_[BSE]->name_, "BSE");
        strcpy(moveObjs_[CE]->name_, "CE");
        
        SEXPTmp = getListElement(opts, "moveProbsList");
        for (ii = MH; ii <= CE; ++ii) {
                name = moveObjs_[ii]->name_;
                moveProbs_[ii] = REAL(getListElement(SEXPTmp, name))[0];
                if (moveProbs_[ii] > 0.0) moveProbsPositive_[ii] = true;
        }
        cumsumProbs_[MH] = moveProbs_[MH];
        for (ii = RC; ii <= CE; ++ii)
                cumsumProbs_[ii] = cumsumProbs_[ii - 1] + moveProbs_[ii];
        
        SEXPTmp = getListElement(opts, "moveNTimesList");
        for (ii = MH; ii <= CE; ++ii) {
                name = moveObjs_[ii]->name_;
                moveNTimes_[ii] = INTEGER(getListElement(SEXPTmp, name))[0];
        }
        
        SEXPTmp = getListElement(opts, "moveSelectionCodesList");
        SEXPCodes = getListElement(SEXPTmp, "RC");
        SEXPTmp = getListElement(opts, "moveSelectionTempersList");
        doublesTmp = REAL(getListElement(SEXPTmp, "RC"));
        for (ii = 0; ii < 2; ++ii) {
                moveSelectionCodes_[RC][ii] = \
                getSelectionCode(CHAR(STRING_ELT(SEXPCodes, ii)));
                moveSelectionTempers_[RC][ii] = doublesTmp[ii];
        }
        
        SEXPTmp = getListElement(opts, "moveSelectionCodesList");
        SEXPCodes = getListElement(SEXPTmp, "SCSC_ONE_NEW");
        SEXPTmp = getListElement(opts, "moveSelectionTempersList");
        doublesTmp = REAL(getListElement(SEXPTmp, "SCSC_ONE_NEW"));
        for (ii = 0; ii < 2; ++ii) {
                moveSelectionCodes_[SCSC_ONE_NEW][ii] = \
                getSelectionCode(CHAR(STRING_ELT(SEXPCodes, ii)));
                moveSelectionTempers_[SCSC_ONE_NEW][ii] = doublesTmp[ii];
        }
        
        SEXPTmp = getListElement(opts, "moveSelectionCodesList");
        SEXPCodes = getListElement(SEXPTmp, "SCSC_TWO_NEW");
        SEXPTmp = getListElement(opts, "moveSelectionTempersList");
        doublesTmp = REAL(getListElement(SEXPTmp, "SCSC_TWO_NEW"));
        for (ii = 0; ii < 2; ++ii) {
                moveSelectionCodes_[SCSC_TWO_NEW][ii] = \
                getSelectionCode(CHAR(STRING_ELT(SEXPCodes, ii)));
                moveSelectionTempers_[SCSC_TWO_NEW][ii] = doublesTmp[ii];
        }
        
        SEXPTmp = getListElement(opts, "moveSelectionCodesList");
        SEXPCodes = getListElement(SEXPTmp, "SCRC");
        SEXPTmp = getListElement(opts, "moveSelectionTempersList");
        doublesTmp = REAL(getListElement(SEXPTmp, "SCRC"));
        ii = 0;
        if (strcmp(CHAR(STRING_ELT(SEXPCodes, ii)), "same") == 0)
                moveSelectionCodes_[SCRC][ii] = SAME_SIZE;
        else if (strcmp(CHAR(STRING_ELT(SEXPCodes, ii)), "random") == 0)
                moveSelectionCodes_[SCRC][ii] = RANDOM_SIZE;
        moveSelectionTempers_[SCRC][ii] = doublesTmp[ii];
        
        SEXPTmp = getListElement(opts, "moveSelectionTempersList");
        doublesTmp = REAL(getListElement(SEXPTmp, "BCE"));
        ii = 0;        
        moveSelectionTempers_[BCE][ii] = doublesTmp[ii];
        SCSC_ONE_NEW_notPossible_ = 0;
        SCSC_TWO_NEW_notPossible_ = 0;
        SCRC_notPossible_ = 0;
        
        movePropCtrs_[MH] = PCMatNew(nLevels_, sampDim_);
        for (ii = RC; ii <= CE; ++ii)
                movePropCtrs_[ii] = PCMatNew(nLevels_, nLevels_);
        
        SEXPTmp = getListElement(opts, "levelsSaveSampFor");
        nLevelsSaveSampFor_ = length(SEXPTmp);
        intsTmp = INTEGER(SEXPTmp);
        levelsSaveSampFor_ = (int *) R_alloc(nLevelsSaveSampFor_, sizeof(int));
        for (ii = 0; ii < nLevelsSaveSampFor_; ++ii)
                levelsSaveSampFor_[ii] = intsTmp[ii] - 1;
        saveFitness_ = (Rboolean) LOGICAL(getListElement(opts, "saveFitness"))[0];        

        /* The user provided functions */
        logTarDensFunc_     = getListElement(opts, "logTarDensFunc");
        logTarDensArgsList_ = (ArgsList1 *) R_alloc(1, sizeof(struct ArgsList1));
        
        SEXPTmp = getListElement(opts, "doCallFunc");
        PROTECT(doCallFuncCall_ = lang4(SEXPTmp, R_NilValue,
                                        R_NilValue, R_NilValue));
        CountNProtected::incrementNProtected( );
        doCallFuncEnv_ = getListElement(opts, "doCallFuncEnv");
        
        SEXPTmp = getListElement(opts, "procTimeFunc");
        if (SEXPTmp == R_NilValue) {
                procTimeFuncCall_ = NULL; procTimeFuncEnv_ = NULL;
        } else {
                PROTECT(procTimeFuncCall_ = lang1(SEXPTmp));
                CountNProtected::incrementNProtected( );
                procTimeFuncEnv_ = getListElement(opts, "procTimeFuncEnv");
        }
        timeDetails_ = (TimeDetails *) R_alloc(1, sizeof(struct TimeDetails));
        
        clusterCurrDraws_      = (Cluster **) R_alloc(nLevels_, sizeof(Cluster *));
        clusterCurrDrawsStore_ = (Cluster **) R_alloc(nLevels_, sizeof(Cluster *));
        clusterPropDraws_      = (Cluster **) R_alloc(nLevels_, sizeof(Cluster *));
        for (ii = 0; ii < nLevels_; ++ii) {
                clusterCurrDraws_[ii] = new Cluster(sampDim_);
                clusterPropDraws_[ii] = new Cluster(sampDim_);
        }
        PROTECT(SEXPCurrDraws_ = allocVector(VECSXP, nLevels_));
        SEXPCurrDrawsStore_    = (SEXP *) R_alloc(nLevels_, sizeof(SEXP));
        PROTECT(SEXPPropDraws_ = allocVector(VECSXP, nLevels_));
        CountNProtected::incrementNProtected(2);
        /*
         * Making the SEXPCurrDraws and SEXPPropDraws components point to the 
         * SEXP data contained in clusterCurrDraws_ and clusterPropDraws_, 
         * respectively.         
         */
        for (ii = 0; ii < nLevels_; ++ii) {
                SEXPTmp = clusterCurrDraws_[ii]->getObjLabelsSEXPData( );
                SET_VECTOR_ELT(SEXPCurrDraws_, ii, SEXPTmp);
                SEXPTmp = clusterPropDraws_[ii]->getObjLabelsSEXPData( );
                SET_VECTOR_ELT(SEXPPropDraws_, ii, SEXPTmp);
                PHONY(PRINT_STUB_SEXP_IARRAY(VECTOR_ELT(SEXPCurrDraws_, ii), ", "););
                PHONY(PRINT_STUB_SEXP_IARRAY(VECTOR_ELT(SEXPPropDraws_, ii), ", "););
        }
        /* filling the starting values in */
        for (ii = 0; ii < nLevels_; ++ii) {
                for (jj = 0; jj < sampDim_; ++jj) {
                        count = jj * nLevels_ + ii;
                        clusterCurrDraws_[ii]->setObjLabel(jj, startingVals[count]);
                        clusterPropDraws_[ii]->setObjLabel(jj, startingVals[count]);
                }                
        }
        
        initArgsList1(logTarDensArgsList_);
        dotsList_ = getListElement(opts, "dotsList");

        /* fixed iter case */
        if (timeInSecs_ <= 0) { drawsArr_ = NULL; }
        
        /* fixed time case */
        drawsArr_ = (double ***) R_alloc(nLevelsSaveSampFor_, sizeof(double **));
        for (ii = 0; ii < nLevelsSaveSampFor_; ++ii) {
                drawsArr_[ii] = (double **) R_alloc(sampDim_, sizeof(double *));
                for (jj = 0; jj < sampDim_; ++jj) 
                        drawsArr_[ii][jj] = (double *) R_alloc(nItersActual_, sizeof(double));
        }
        init( );
}

Sampler::~Sampler (void)
{
        PHONY(Rprintf("Sampler destructor was called\n"););
        for (int ii = 0; ii < nLevels_; ++ii) {
                delete clusterCurrDraws_[ii];
                delete clusterPropDraws_[ii];
        }
}

SEXP
Sampler::run (void)
{       
        int nComps = 0, comp = 0;
        SEXP draws, acceptRatios, acceptRatiosList;
        SEXP samplerObj = R_NilValue, names;
        
        draws                = makeDraws( ); ++nComps;
        acceptRatios         = makeAcceptRatios( ); ++nComps;
        acceptRatiosList     = makeAcceptRatiosList( ); ++nComps;
        PROTECT(samplerObj   = allocVector(VECSXP, nComps));
        CountNProtected::incrementNProtected( );
        PROTECT(names        = allocVector(STRSXP, nComps)); 
        
        /* fill in the samplerObj */
        SET_VECTOR_ELT(samplerObj, comp, draws);
        SET_STRING_ELT(names, comp, mkChar("draws")); ++comp;
        SET_VECTOR_ELT(samplerObj, comp, acceptRatios);
        SET_STRING_ELT(names, comp, mkChar("acceptRatios")); ++comp;
        SET_VECTOR_ELT(samplerObj, comp, acceptRatiosList);
        SET_STRING_ELT(names, comp, mkChar("acceptRatiosList")); ++comp;
        setAttrib(samplerObj, R_NamesSymbol, names);
        PHONY(print( );
              PRINT_STUB_INT(CountNProtected::getNProtected( )););
        UNPROTECT(1 + CountNProtected::getNProtected( ));
        CountNProtected::resetNProtected( );
        return samplerObj;
}

/*****************************************************************************/
// The private functions: the MH move and its helpers.
/*****************************************************************************/

int 
Sampler::MH_SPLIT_MERGE_propNew (int obj, Cluster &curr, Cluster &prop)
{
        int freq = prop.getClusterFreqFromObj(obj);
        int nClusters = prop.getNClusters( );
        // freq == 1: merge is the only possibility
        if (freq == 1) {
                prop.merge(obj, prop, MH_OLC_);
        } else {
                // freq > 1 && nCluster == 1: split is the only possibility
                if (nClusters == 1) {
                        prop.split(obj, prop, MH_OLC_);
                } else {
                        // merge with probability MHMergeProb_
                        if (runif(0, 1) <= MHMergeProb_)
                                prop.merge(obj, prop, MH_OLC_);
                        // and split otherwise
                        else
                                prop.split(obj, prop, MH_OLC_);
                }
        }
        return 0;
}

double
Sampler::MH_SPLIT_MERGE_propDens (int obj, Cluster &curr, Cluster &prop)
{
        int propFreq = prop.getClusterFreqFromObj(obj);
        int currNClusters = curr.getNClusters( );        
        if (propFreq == 1) {
                /* 
                 * propFreq == 1 and currNClusters == 1:  split was the only
                 * possibility, thus return log(1.0)
                 */
                if (currNClusters == 1) {
                        return 0.0;
                } 
                /* 
                 * propFreq == 1 and currNClusters > 1:  both merge and split
                 * were the options and split was chosen
                 */
                else {
                        return log(1.0 - MHMergeProb_);
                }
        } else {
                int currFreq = curr.getClusterFreqFromObj(obj);
                /* 
                 * propFreq > 1 and currNClusters == 1:  merge was the only
                 * possibility
                 */
                if (currFreq == 1) {
                        return -log(static_cast<double>(currNClusters - 1));
                }
                /* 
                 * propFreq > 1 and currNClusters > 1:  both merge and split
                 * were the options and merge was chosen
                 */
                else {
                        return log(MHMergeProb_ / (currNClusters - 1));
                }
        }
        RAISE("should not reach here");
        return 0;
}

int
Sampler::MH_SPLIT_MERGE_move (void)
{
        int ii, ll = thisLevel_;
        double logPropDens, logNumerators[2], logDenominators[2];
        double logAlpha, alpha;
        ProposalCounter **ppc = movePropCtrs_[MH][ll];
        Cluster *currPtr = clusterCurrDraws_[ll], *propPtr = clusterPropDraws_[ll];
        SEXP SEXPCurr = VECTOR_ELT(SEXPCurrDraws_, ll);
        SEXP SEXPProp = VECTOR_ELT(SEXPPropDraws_, ll);        

        currPtr->tabulate( );
        propPtr->copy(*currPtr);
        for (int bb = 0; bb < MHNBlocks_; ++bb) {
                MH_SPLIT_MERGE_propNew(bb, *currPtr, *propPtr);
                logPropDens = (this->*logTarDensFP_)(SEXPProp);
                logDenominators[0] = logDensities_[ll] / temperLadder_[ll];
                logNumerators[0] = logPropDens / temperLadder_[ll];
                logDenominators[1] = MH_SPLIT_MERGE_propDens(bb, *currPtr, *propPtr);
                logNumerators[1] = MH_SPLIT_MERGE_propDens(bb, *propPtr, *currPtr);
                for (logAlpha = 0.0, ii = 0; ii < 2; ++ii)
                        logAlpha += logNumerators[ii] - logDenominators[ii];
                alpha = exp(logAlpha); alpha = MIN(1.0, alpha);
                if (R_FINITE(alpha) == FALSE) {
                        Rprintf("MH: level: %d | block: %d | iter: %d | " \
                                "alpha: %5.4g\n", ll, bb, thisIter_, alpha);
                        Rprintf("The current draw: logDensity: %g\n", logDenominators[0]);
                        utils_SEXP_iarray_print(SEXPCurr, ", ");
                        Rprintf("The proposal draw: logDensity: %g\n", logNumerators[0]);
                        utils_SEXP_iarray_print(SEXPProp, ", ");
                        error("alpha is non finite, inspect the above numbers\n");                
                }
                ppc[bb]->proposed_ += 1;        
                if (verboseLevel_ >= 100) 
                        Rprintf("MH: level: %d | block: %d | iter: %d | " \
                                "alpha: %5.4g\n", ll, bb, thisIter_, alpha);
        
                /* MH acceptance rejection step */
                if (runif(0, 1) <= alpha) {
                        if (verboseLevel_ >= 10) 
                                Rprintf("MH: level: %d | block: %d | iter: %d | " \
                                        "alpha: %5.4g [*** accepted]\n",
                                        ll, bb, thisIter_, alpha);
                        currPtr->registerChangeOfLabel(MH_OLC_->obj,
                                                       MH_OLC_->oldLabel,
                                                       MH_OLC_->newLabel);
                        logDensities_[ll] = logPropDens;
                        ppc[bb]->accepted_ += 1;
                } else {
                        propPtr->registerChangeOfLabel(MH_OLC_->obj,
                                                       MH_OLC_->newLabel,
                                                       MH_OLC_->oldLabel);       
                }
        }
        return 0;
}

int
Sampler::MH_GIBBS_move (void)
{
        // WRITE ME
        return 0;
}


/*****************************************************************************/
// The private functions: the RC move and its helpers.
/*****************************************************************************/
/*
 * Note: Both in RC_RANDOM_RANDOM_probAndLevels and
 * RC_RANDOM_RANDOM_probGivenLevels, we are setting *prob = 1.0, since
 * it's going to get cancelled in the expression
 * '(childrenToParentsProb / parentsToChildrenProb)' below.
 */
int 
Sampler::RC_RANDOM_RANDOM_probAndLevels (double *prob, int *sl)
{
        sl[0] = (int) floor(runif(0, nLevels_));
        int pos = (int) floor(runif(0, nLevels_ - 1));
        for (int ii = 0; ii < nLevels_; ++ii) {
                if (ii == pos) {
                        *prob = 1.0; sl[1] = ii + 1; return 0;
                }
        }
        RAISE("should not reach here, check your C++ code");
        return 0;
}

int 
Sampler::RC_RANDOM_RANDOM_probGivenLevels (double *prob, int *sl, double *logPropDens)
{
        *prob = 1.0; return 0;
}

int 
Sampler::RC_BEST_BEST_probAndLevels (double *prob, int *sl)
{
        int ll;
        double *ld = logDensities_, selTemper = moveSelectionTempers_[RC][0];
        double *lw = scratch_SLC_->logWeights_, *aw = scratch_SLC_->adjWeights_;
        double *ps = scratch_SLC_->partialSum_;
        /* A: All, N: Numerator, S: Sum, B: But, O: One */
        double mlwA, parentNA[2], parentSA, mlwBO[2], parentNBO[2], parentSBO[2], uu;

        /* choose the first parent */
        for (mlwA = R_NegInf, ll = 0; ll < nLevels_; ++ll) {
                lw[ll] = ld[ll] / selTemper; mlwA = MAX(mlwA, lw[ll]);
        }
        aw[0] = exp(lw[0] - mlwA); ps[0] = aw[0];
        for (ll = 1; ll < nLevels_; ++ll) {
                aw[ll] = exp(lw[ll] - mlwA); ps[ll] = ps[ll - 1] + aw[ll];
        }
        parentSA = ps[nLevels_ - 1]; uu = runif(0, parentSA);
        for (ll = 0; ll < nLevels_; ++ll)
                if (uu <= ps[ll]) { sl[0] = ll; break; }
        parentNA[0] = aw[sl[0]];

        /* choose the second parent */
        for (mlwBO[1] = R_NegInf, ll = 0; ll < nLevels_; ++ll)
                if (ll != sl[0]) mlwBO[1] = MAX(mlwBO[1], lw[ll]);
        if (sl[0] == 0) ps[0] = 0.0;
        else            ps[0] = aw[0];
        for (ll = 1; ll < nLevels_; ++ll) 
                if (ll == sl[0]) ps[ll] = ps[ll - 1];
                else             ps[ll] = ps[ll - 1] + exp(lw[ll] - mlwBO[1]);
        parentSBO[1] = ps[nLevels_ - 1]; uu = runif(0, parentSBO[1]);
        for (ll = 0; ll < nLevels_; ++ll)
                if (uu <= ps[ll]) { sl[1] = ll; break; }
        parentNA[1] = aw[sl[1]];
        parentNBO[1] = exp(lw[sl[1]] - mlwBO[1]);
        
        /* do the rest of the probability computation */
        for (mlwBO[0] = R_NegInf, ll = 0; ll < nLevels_; ++ll)
                if (ll != sl[1]) mlwBO[0] = MAX(mlwBO[0], lw[ll]);
        parentNBO[0] = exp(lw[sl[0]] - mlwBO[0]);
        for (parentSBO[0] = 0.0, ll = 0; ll < nLevels_; ++ll) 
                if (ll != sl[1]) parentSBO[0] += exp(lw[ll] - mlwBO[0]);
        *prob = (parentNA[0] / parentSA) * (parentNBO[1] / parentSBO[1]) + \
                (parentNA[1] / parentSA) * (parentNBO[0] / parentSBO[0]);
        return 0;
}

int 
Sampler::RC_BEST_BEST_probGivenLevels (double *prob, int *sl, double *logPropDens)
{
        int ll;
        double *ld = logDensities_, selTemper = moveSelectionTempers_[RC][0];
        double *lw = scratch_SLC_->logWeights_;
        /* A: All, N: Numerator, S: Sum, B: But, O: One */
        double mlwA, mlwBO[2] = { 0.0, 0.0 };
        double childNA[2] = { 0.0, 0.0 }, childSA;
        double childNBO[2] = { 0.0, 0.0 }, childSBO[2] = { 0.0, 0.0 };

        mlwA = R_NegInf; mlwBO[0] = R_NegInf; mlwBO[1] = R_NegInf;
        for (ll = 0; ll < nLevels_; ++ll) {
                if (ll == sl[0]) {
                        lw[ll] = logPropDens[0] / selTemper;
                        mlwBO[0] = MAX(mlwBO[0], lw[ll]);
                } else if (ll == sl[1]) {
                        lw[ll] = logPropDens[1] / selTemper;
                        mlwBO[1] = MAX(mlwBO[1], lw[ll]);
                } else {
                        lw[ll] = ld[ll] / selTemper;
                        mlwBO[0] = MAX(mlwBO[0], lw[ll]);
                        mlwBO[1] = MAX(mlwBO[1], lw[ll]);
                }
                mlwA = MAX(mlwA, lw[ll]);                
        }
        childSA = 0.0; childSBO[0] = 0.0; childSBO[1] = 0.0;
        for (ll = 0; ll < nLevels_; ++ll) {
                if (ll == sl[0]) {
                        childNA[0] = exp(lw[ll] - mlwA); childSA += childNA[0];
                        childNBO[0] = exp(lw[ll] - mlwBO[0]); childSBO[0] += childNBO[0];
                } else if (ll == sl[1]) {
                        childNA[1] = exp(lw[ll] - mlwA); childSA += childNA[1];
                        childNBO[1] = exp(lw[ll] - mlwBO[1]); childSBO[1] += childNBO[1];
                } else {
                        childSA += exp(lw[ll] - mlwA);
                        childNBO[0] += exp(lw[ll] - mlwBO[0]);
                        childNBO[1] += exp(lw[ll] - mlwBO[1]);
                }
        }
        *prob = (childNA[0] / childSA) * (childNBO[1] / childSBO[1]) + \
                (childNA[1] / childSA) * (childNBO[0] / childSBO[0]);
        return 0;
}

int
Sampler::RC_move (void)
{
        int ii, jj, sl[N_SELECTIONS], crosspoint;
        double *ld = logDensities_, *tl = temperLadder_;
        SEXP parents[N_SELECTIONS], children[N_SELECTIONS];
        int *parentVals[N_SELECTIONS], *childVals[N_SELECTIONS];
        double parentsToChildrenProb, childrenToParentsProb, logPropDens[N_SELECTIONS];
        double sum, alpha;

        (this->*RC_probAndLevelsFP_)(&parentsToChildrenProb, sl);
        for (ii = 0; ii < N_SELECTIONS; ++ii) {
                parents[ii] = VECTOR_ELT(SEXPCurrDraws_, sl[ii]);
                parentVals[ii] = INTEGER(parents[ii]);
                children[ii] = VECTOR_ELT(SEXPPropDraws_, sl[ii]);
                childVals[ii] = INTEGER(children[ii]);
        }
        /*
         * Here we choose a crosspoint in \{ 2:d \}, making it a
         * proper crosssover as opposed to exchange, as is the case
         * when coordinate 1 and d are chosen.
         */
        crosspoint = (int) floor(runif(1, sampDim_));
        for (jj = 0; jj < crosspoint; ++jj) {
                childVals[0][jj] = parentVals[0][jj];
                childVals[1][jj] = parentVals[1][jj];
        }
        for (jj = crosspoint; jj < sampDim_; ++jj) {
                childVals[0][jj] = parentVals[1][jj];
                childVals[1][jj] = parentVals[0][jj];
        }
        for (sum = 0.0, ii = 0; ii < 2; ++ii) {
                logPropDens[ii] = (this->*logTarDensFP_)(children[ii]);
                sum += (logPropDens[ii] - ld[sl[ii]]) / tl[sl[ii]];
        }
        (this->*RC_probGivenLevelsFP_)(&childrenToParentsProb, sl, logPropDens);
        alpha = exp(sum) * (childrenToParentsProb / parentsToChildrenProb);
        alpha = MIN(1.0, alpha);
        ProposalCounter *pc = movePropCtrs_[RC][sl[0]][sl[1]];
        pc->proposed_ += 1;
        if (verboseLevel_ >= 100)
                Rprintf("RC: levels: %d, %d | iter: %d | " \
                        "alpha: %5.4g\n", sl[0], sl[1], thisIter_, alpha);
        /* the acceptance rejection step */
        if (runif(0, 1) <= alpha) {
                pc->accepted_ += 1;
                if (verboseLevel_ >= 10)
                        Rprintf("RC: levels: %d, %d | iter: %d | " \
                                "alpha: %5.4g [*** accepted]\n",
                                sl[0], sl[1], thisIter_, alpha);
                for (ii = 0; ii < N_SELECTIONS; ++ii) {
                        SEXP SEXPTmp = VECTOR_ELT(SEXPCurrDraws_, sl[ii]);
                        SET_VECTOR_ELT(SEXPCurrDraws_, sl[ii],
                                       VECTOR_ELT(SEXPPropDraws_, sl[ii]));
                        SET_VECTOR_ELT(SEXPPropDraws_, sl[ii], SEXPTmp);
                        SWAP(Cluster *, clusterCurrDraws_[sl[ii]], clusterPropDraws_[sl[ii]]);
                        ld[sl[ii]] = logPropDens[ii];
                        clusterCurrDraws_[sl[ii]]->forceTabulate( );
                }
        }
        return 0;        
}

/*****************************************************************************/
// The private functions: the SCSC move and its helpers.
/*****************************************************************************/

int
Sampler::SCSC_ONE_NEW_RANDOM_RANDOM_probAndLevels (double *probs, int *sl)
{
        sl[0] = (int) floor(runif(0, nLevels_));
        if (sl[0] == 0)                   sl[1] = 1;
        else if (sl[0] == (nLevels_ - 1)) sl[1] = nLevels_ - 2;
        else {
                if (runif(0, 1) <= 0.5) sl[1] = sl[0] - 1;
                else                    sl[1] = sl[0] + 1;
        }
        // choose the (survivor, nonSurvivor) parents
        if (runif(0, 1) <= 0.5) {
                sl[SURV_PARENT_POS] = sl[0];
                sl[NON_SURV_PARENT_POS] = sl[1];
        } else {
                sl[SURV_PARENT_POS] = sl[1];
                sl[NON_SURV_PARENT_POS] = sl[0];
        }
        // CHECK ME: WILL THE FOLLOWING PROBABILITIES CANCEL OUT?
        probs[0] = 1.0; probs[1] = 1.0;
        return 0;
}       

int
Sampler::SCSC_ONE_NEW_RANDOM_RANDOM_probGivenLevels (double *probs, int *sl, double *logPropDens)
{
        // CHECK ME: WILL THE FOLLOWING PROBABILITIES CANCEL OUT?
        probs[0] = 1.0; probs[1] = 1.0;
        return 0;
}

/*
 * Note: The probAndLevels and probGivenLevels are a little different
 * for SCSC_ONE_NEW than for RC and SCSC_TWO_NEW. Here 'probs' is of
 * length 2 and 'sl' is of length 4.
 */
int
Sampler::SCSC_ONE_NEW_BEST_BEST_probAndLevels (double *probs, int *sl)
{
        int ll, ii, slHyp[2 * N_SELECTIONS];
        double *ld = logDensities_, selTemper = moveSelectionTempers_[SCSC_ONE_NEW][0];
        double *lw = scratch_SLC_->logWeights_, *aw = scratch_SLC_->adjWeights_;
        double *ps = scratch_SLC_->partialSum_;
        /* A: All, S: Sum */
        double mlwA, parentSA, uu;
        
        // following initialization avoids g++ warnings
        for (ii = 0; ii < 2 * N_SELECTIONS; ++ii)
                slHyp[ii] = 0;
        /* choose the first parent */
        for (mlwA = R_NegInf, ll = 0; ll < nLevels_; ++ll) {
                lw[ll] = ld[ll] / selTemper; mlwA = MAX(mlwA, lw[ll]);
        }
        aw[0] = exp(lw[0] - mlwA); ps[0] = aw[0];
        for (ll = 1; ll < nLevels_; ++ll) {
                aw[ll] = exp(lw[ll] - mlwA); ps[ll] = ps[ll - 1] + aw[ll];
        }
        parentSA = ps[nLevels_ - 1]; uu = runif(0, parentSA);
        for (ll = 0; ll < nLevels_; ++ll)
                if (uu <= ps[ll]) { sl[0] = ll; break; }

        /* choose the second parent */
        if (sl[0] == 0)                   sl[1] = 1;
        else if (sl[0] == (nLevels_ - 1)) sl[1] = nLevels_ - 2;
        else {
                if (runif(0, 1) < 0.5) sl[1] = sl[0] - 1;
                else                   sl[1] = sl[0] + 1;
        }

        // choose the (survivor, nonSurvivor) parents
        double survivorProb = aw[sl[0]] / (aw[sl[0]] + aw[sl[1]]);
        if (runif(0, 1) <= survivorProb) {
                sl[SURV_PARENT_POS] = sl[0];
                sl[NON_SURV_PARENT_POS] = sl[1];
        } else {
                sl[SURV_PARENT_POS] = sl[1];
                sl[NON_SURV_PARENT_POS] = sl[0];
                survivorProb = 1 - survivorProb;
        }        
        
        // the probability computation
        probs[0] = ((aw[sl[0]] / parentSA) * getNeighbourProb(sl[1], sl[0]) + \
                   (aw[sl[1]] / parentSA) * getNeighbourProb(sl[0], sl[1])) * \
                survivorProb;

        // the hypothetical survivor parent computation
        if ((sl[NON_SURV_PARENT_POS] == 0) ||
            (sl[NON_SURV_PARENT_POS] == (nLevels_ - 1)))
                return 0;
        // only if it's needed
        if (sl[SURV_PARENT_POS] == (sl[NON_SURV_PARENT_POS] - 1)) {
                slHyp[0] = sl[NON_SURV_PARENT_POS] + 1;
                slHyp[1] = sl[NON_SURV_PARENT_POS];
                slHyp[SURV_PARENT_POS] = slHyp[0];
                slHyp[NON_SURV_PARENT_POS] = slHyp[1];
        } else if (sl[SURV_PARENT_POS] == (sl[NON_SURV_PARENT_POS] + 1)) {
                slHyp[0] = sl[NON_SURV_PARENT_POS] - 1;
                slHyp[1] = sl[NON_SURV_PARENT_POS];
                slHyp[SURV_PARENT_POS] = slHyp[0];
                slHyp[NON_SURV_PARENT_POS] = slHyp[1];
        } else {
                RAISE("should not reach here");
        }
        survivorProb = aw[slHyp[SURV_PARENT_POS]] / (aw[slHyp[0]] + aw[slHyp[1]]);
        probs[1] = ((aw[slHyp[0]] / parentSA) * getNeighbourProb(slHyp[1], slHyp[0]) + \
                    (aw[slHyp[1]] / parentSA) * getNeighbourProb(slHyp[0], slHyp[1])) * \
                survivorProb;
        return 0;        
}

int
Sampler::SCSC_ONE_NEW_BEST_BEST_probGivenLevels (double *probs, int *sl, double *logPropDens)
{
        int ll, ii, slHyp[2 * N_SELECTIONS];
        double *ld = logDensities_, selTemper = moveSelectionTempers_[SCSC_ONE_NEW][0];
        double *lw = scratch_SLC_->logWeights_, *aw = scratch_SLC_->adjWeights_;
        /* A: All, S: Sum */
        double mlwA, childSA, survivorProb;

        // following initialization avoids g++ warnings
        for (ii = 0; ii < 2 * N_SELECTIONS; ++ii)
                slHyp[ii] = 0;
        for (mlwA = R_NegInf, ll = 0; ll < nLevels_; ++ll) {
                if (ll == sl[NON_SURV_PARENT_POS]) 
                        lw[ll] = logPropDens[0] / selTemper;
                else
                        lw[ll] = ld[ll] / selTemper;
                mlwA = MAX(mlwA, lw[ll]);
        }
        for (childSA = 0.0, ll = 0; ll < nLevels_; ++ll) {
                aw[ll] = exp(lw[ll] - mlwA); childSA += aw[ll];
        }
        
        survivorProb = aw[sl[SURV_PARENT_POS]] / (aw[sl[0]] + aw[sl[1]]);
        probs[0] = ((aw[sl[0]] / childSA) * getNeighbourProb(sl[1], sl[0]) + \
                    (aw[sl[1]] / childSA) * getNeighbourProb(sl[0], sl[1])) * \
                survivorProb;

        // the hypothetical survivor parent computation
        if ((sl[NON_SURV_PARENT_POS] == 0) ||
            (sl[NON_SURV_PARENT_POS] == (nLevels_ - 1)))
                return 0;
        // only if it's needed
        if (sl[SURV_PARENT_POS] == (sl[NON_SURV_PARENT_POS] - 1)) {
                slHyp[0] = sl[NON_SURV_PARENT_POS] + 1;
                slHyp[1] = sl[NON_SURV_PARENT_POS];
                slHyp[SURV_PARENT_POS] = slHyp[0];
                slHyp[NON_SURV_PARENT_POS] = slHyp[1];
        } else if (sl[SURV_PARENT_POS] == (sl[NON_SURV_PARENT_POS] + 1)) {
                slHyp[0] = sl[NON_SURV_PARENT_POS] - 1;
                slHyp[1] = sl[NON_SURV_PARENT_POS];
                slHyp[SURV_PARENT_POS] = slHyp[0];
                slHyp[NON_SURV_PARENT_POS] = slHyp[1];
        } else 
                RAISE("should not reach here");
        survivorProb = aw[slHyp[SURV_PARENT_POS]] / (aw[slHyp[0]] + aw[slHyp[1]]);
        probs[1] = ((aw[slHyp[0]] / childSA) * getNeighbourProb(slHyp[1], slHyp[0]) + \
                    (aw[slHyp[1]] / childSA) * getNeighbourProb(slHyp[0], slHyp[1])) * \
                survivorProb;
        return 0;
}

bool 
Sampler::SCSC_ONE_NEW_isParent (Cluster *plausibleParent, Cluster *nonSurvParent,
                                Cluster *modChild)
{
        int ii, nObjs = plausibleParent->getNObjs( );
        int *plausibleParentLabels = plausibleParent->getObjLabelsData( );
        int *nonSurvParentLables = nonSurvParent->getObjLabelsData( );
        int *modChildLabels = modChild->getObjLabelsData( );
        int count = 0, nPositive = 0, nNegative = 0, parentLabel = 0;
        
        if (nObjs != sampDim_)
                RAISE("we are in trouble");
        for (ii = 0; ii < nObjs; ++ii) {
                labelDiffs_[ii] = nonSurvParentLables[ii] - modChildLabels[ii];
                if (labelDiffs_[ii] == 0) continue;

                if (labelDiffs_[ii] > 0) {
                        labelDiffsPositive_[nPositive] = labelDiffs_[ii]; ++nPositive;
                } else if (labelDiffs_[ii] < 0) {
                        labelDiffsNegative_[nNegative] = labelDiffs_[ii]; ++nNegative;
                }
                
                if (count == 0) {
                        parentLabel = plausibleParentLabels[ii]; ++count;
                } else {
                        /*
                         * Plausible parent should have the same label
                         * responsible for non-zero label differences.
                         */
                        if (parentLabel != plausibleParentLabels[ii]) return false;
                }
        }

        if ((nPositive == 0) || (nNegative == 0)) return false;
        
        qsort(labelDiffsPositive_, nPositive, sizeof(int), intcmp);
        qsort(labelDiffsNegative_, nNegative, sizeof(int), intcmp);
        
        int commonVal = labelDiffsPositive_[0];
        for (ii = 1; ii < nPositive; ++ii)
                if (labelDiffsPositive_[ii] != commonVal) return false;

        for (ii = 0; ii < nNegative; ++ii)
                if (labelDiffsNegative_[ii] != -commonVal) return false;        
        return true;
}

int
Sampler::SCSC_ONE_NEW_propNew (Cluster *survParent, Cluster *nonSurvParent,
                               Cluster *modChild)
{
        int ii, jj, survParentNClusters, nonSurvParentNClusters;
        int survClusterSampled, nonSurvClustersSampled[N_SUB_CLUSTERS];
        int nn, nn1, *iarr1, nn2, *iarr2;
        double incr, sum, uu;
        
        survParentNClusters = survParent->getNClusters( );
        nonSurvParentNClusters = nonSurvParent->getNClusters( );
        // SCSC_ONE_NEW cannot be performed
        if (nonSurvParentNClusters < 2) return 1;

        for (nSCXX_rows2OrMore_ = 0, ii = 0; ii < survParentNClusters; ++ii) {
                nn1 = survParent->getClusterFreq(ii);
                iarr1 = survParent->getClusterMembersData(ii);
                int count = 0;
                for (jj = 0; jj < nonSurvParentNClusters; ++jj) {
                        nn2 = nonSurvParent->getClusterFreq(jj);
                        iarr2 = nonSurvParent->getClusterMembersData(jj);
                        bool tmp = utils_are_disjoint_iarray_iarray(nn1, iarr1, nn2, iarr2);
                        if (tmp == false) {
                                SCXX_nonEmptyInRow_[ii][count] = jj; ++count;
                        }
                }
                nSCXX_nonEmptyInRow_[ii] = count;
                if (nSCXX_nonEmptyInRow_[ii] >= 2) {
                        SCXX_rows2OrMore_[nSCXX_rows2OrMore_] = ii;
                        ++nSCXX_rows2OrMore_;
                }
        }
        
        PHONY(PRINT_STUB_INT(nSCXX_rows2OrMore_);
              utils_iarray_print(SCXX_rows2OrMore_, nSCXX_rows2OrMore_, ", "););

        // SCSC_ONE_NEW cannot be performed
        if (nSCXX_rows2OrMore_ == 0) return 1;
        
        // choose the survParent sampled cluster
        ii = (int) floor(nSCXX_rows2OrMore_ * runif(0, 1));
        survClusterSampled = SCXX_rows2OrMore_[ii];

        // choose the nonSurvParent sampled clusters
        nn = nSCXX_nonEmptyInRow_[survClusterSampled];
        jj = (int) floor(nn * runif(0, 1));
        nonSurvClustersSampled[0] = SCXX_nonEmptyInRow_[survClusterSampled][jj];
        
        incr = 1.0 / (nn - 1); sum = 0.0; uu = runif(0, 1);
        for (jj = 0; jj < nn; ++jj) {
                int kk = SCXX_nonEmptyInRow_[survClusterSampled][jj];
                if (kk == nonSurvClustersSampled[0]) continue;
                sum += incr;
                if (uu <= sum) {
                        nonSurvClustersSampled[1] = kk; break;
                }
        }

        // change the labels of the two SC intersections
        nn1 = survParent->getClusterFreq(survClusterSampled);
        iarr1 = survParent->getClusterMembersData(survClusterSampled);
        for (ii = 0; ii < N_SUB_CLUSTERS; ++ii) {
                nn2 = nonSurvParent->getClusterFreq(nonSurvClustersSampled[ii]);
                iarr2 = nonSurvParent->getClusterMembersData(nonSurvClustersSampled[ii]);
                utils_intersect_iarray_iarray(nn1, iarr1, nn2, iarr2,
                                              nSCXX_intersection_ + ii,
                                              SCXX_intersection_[ii]);

                PHONY(Rprintf("SCSC_ONE_NEW_propNew: the sampled clusters: %d, %d\n",
                              survClusterSampled, nonSurvClustersSampled[ii]);
                      Rprintf("SCXX_intersection_ [length: %d]:\n", nSCXX_intersection_[ii]);
                      survParent->prettyPrint( );
                      nonSurvParent->prettyPrint( );
                      utils_iarray_print(SCXX_intersection_[ii], nSCXX_intersection_[ii], ", "););

                int label = nonSurvParent->getClusterLabel(nonSurvClustersSampled[1 - ii]);
                modChild->setObjLabelsSliceFromOne(nSCXX_intersection_[ii],
                                                   SCXX_intersection_[ii], label);
        }
        return 0;
}

int 
Sampler::SCSC_ONE_NEW_move (void)
{
        int ii, sl[2 * N_SELECTIONS];
        double *ld = logDensities_, *tl = temperLadder_;
        Cluster *survParent, *nonSurvParent, *modChild;
        double parentsToChildrenProb, childrenToParentsProb, logPropDens;
        double probsPTC[2], probsCTP[2], sum, alpha;
        bool isParent;
        
        (this->*SCSC_ONE_NEW_probAndLevelsFP_)(probsPTC, sl);
        parentsToChildrenProb = probsPTC[0];        
        survParent = clusterCurrDraws_[sl[SURV_PARENT_POS]];
        survParent->tabulate( );
        nonSurvParent = clusterCurrDraws_[sl[NON_SURV_PARENT_POS]];
        nonSurvParent->tabulate( );
        modChild = clusterPropDraws_[sl[NON_SURV_PARENT_POS]];
        modChild->copy(*nonSurvParent);
        int err = SCSC_ONE_NEW_propNew(survParent, nonSurvParent, modChild);
        if (err) {
                PHONY(Rprintf("---> Could not do SCSC_ONE_NEW_propNew\n"););
                ++SCSC_ONE_NEW_notPossible_; return 1;
        }
        SEXP SEXPModChild = VECTOR_ELT(SEXPPropDraws_, sl[NON_SURV_PARENT_POS]);
        logPropDens = (this->*logTarDensFP_)(SEXPModChild);
        sum = (logPropDens - ld[sl[NON_SURV_PARENT_POS]]) / tl[sl[NON_SURV_PARENT_POS]];
        (this->*SCSC_ONE_NEW_probGivenLevelsFP_)(probsCTP, sl, &logPropDens);
        childrenToParentsProb = probsCTP[0];

        // checking the 'isParent' stuff
        if ((sl[NON_SURV_PARENT_POS] > 0) &&
            (sl[NON_SURV_PARENT_POS] < (nLevels_ - 1))) {
                Cluster *plausibleParent = NULL;
                
                if (sl[SURV_PARENT_POS] == (sl[NON_SURV_PARENT_POS] - 1))
                        plausibleParent = clusterCurrDraws_[sl[NON_SURV_PARENT_POS] + 1];
                else if (sl[SURV_PARENT_POS] == (sl[NON_SURV_PARENT_POS] + 1))
                        plausibleParent = clusterCurrDraws_[sl[NON_SURV_PARENT_POS] - 1];
                else
                        RAISE("should not reach here");
                plausibleParent->tabulate( );
                isParent = SCSC_ONE_NEW_isParent(plausibleParent, nonSurvParent, modChild);
                if (isParent == true) {
                        PHONY(Rprintf("---> SCSC_ONE_NEW isParent is true\n"););
                        parentsToChildrenProb += probsPTC[1];
                        childrenToParentsProb += probsCTP[1];
                }
        }
        
        alpha = exp(sum) * (childrenToParentsProb / parentsToChildrenProb);
        alpha = MIN(1.0, alpha);
        ProposalCounter *pc = movePropCtrs_[SCSC_ONE_NEW][sl[0]][sl[1]];
        pc->proposed_ += 1;
        if (verboseLevel_ >= 100)
                Rprintf("SCSC_ONE_NEW: levels: %d, %d | iter: %d | " \
                        "alpha: %5.4g\n", sl[0], sl[1], thisIter_, alpha);
        /* the acceptance rejection step */
        if (runif(0, 1) <= alpha) {
                pc->accepted_ += 1;
                if (verboseLevel_ >= 10)
                        Rprintf("SCSC_ONE_NEW: levels: %d, %d | iter: %d | " \
                                "alpha: %5.4g [*** accepted]\n",
                                sl[0], sl[1], thisIter_, alpha);
                
                PHONY(Rprintf("SCSC_ONE_NEW changes -fitness values as:\n");
                      ii = NON_SURV_PARENT_POS;
                      Rprintf("%d: %g ---> %g\n", sl[ii], ld[sl[ii]], logPropDens);
                      Rprintf("SCSC_ONE_NEW SCXX_intersection:\n");
                      for (ii = 0; ii < N_SUB_CLUSTERS; ++ii) {
                              PRINT_STUB_INT(nSCXX_intersection_[ii]);
                              utils_iarray_print(SCXX_intersection_[ii], nSCXX_intersection_[ii], ", ");
                      });
                
                // only swap things for the non-surviving parent.
                ii = NON_SURV_PARENT_POS;
                SEXP SEXPTmp = VECTOR_ELT(SEXPCurrDraws_, sl[ii]);
                SET_VECTOR_ELT(SEXPCurrDraws_, sl[ii],
                               VECTOR_ELT(SEXPPropDraws_, sl[ii]));
                SET_VECTOR_ELT(SEXPPropDraws_, sl[ii], SEXPTmp);
                SWAP(Cluster *, clusterCurrDraws_[sl[ii]], clusterPropDraws_[sl[ii]]);
                ld[sl[ii]] = logPropDens;
                clusterCurrDraws_[sl[ii]]->tabulate( );
        }        
        return 0;
}

int
Sampler::SCSC_TWO_NEW_RANDOM_RANDOM_probAndLevels (double *prob, int *sl)
{
        sl[0] = (int) floor(runif(0, nLevels_));
        int pos = (int) floor(runif(0, nLevels_ - 1));
        for (int ii = 0; ii < nLevels_; ++ii) {
                if (ii == pos) {
                        *prob = 1.0; sl[1] = ii + 1; return 0;
                }
        }
        RAISE("should not reach here, check your C++ code");
        return 0;
}

int
Sampler::SCSC_TWO_NEW_RANDOM_RANDOM_probGivenLevels (double *prob, int *sl, double *logPropDens)
{
        *prob = 1.0; return 0;
}

int
Sampler::SCSC_TWO_NEW_BEST_BEST_probAndLevels (double *prob, int *sl)
{
        int ll;
        double *ld = logDensities_, selTemper = moveSelectionTempers_[RC][0];
        double *lw = scratch_SLC_->logWeights_, *aw = scratch_SLC_->adjWeights_;
        double *ps = scratch_SLC_->partialSum_;
        /* A: All, N: Numerator, S: Sum, B: But, O: One */
        double mlwA, parentNA[2], parentSA, mlwBO[2], parentNBO[2], parentSBO[2], uu;

        /* choose the first parent */
        for (mlwA = R_NegInf, ll = 0; ll < nLevels_; ++ll) {
                lw[ll] = ld[ll] / selTemper; mlwA = MAX(mlwA, lw[ll]);
        }
        aw[0] = exp(lw[0] - mlwA); ps[0] = aw[0];
        for (ll = 1; ll < nLevels_; ++ll) {
                aw[ll] = exp(lw[ll] - mlwA); ps[ll] = ps[ll - 1] + aw[ll];
        }
        parentSA = ps[nLevels_ - 1]; uu = runif(0, parentSA);
        for (ll = 0; ll < nLevels_; ++ll)
                if (uu <= ps[ll]) { sl[0] = ll; break; }
        parentNA[0] = aw[sl[0]];

        /* choose the second parent */
        for (mlwBO[1] = R_NegInf, ll = 0; ll < nLevels_; ++ll)
                if (ll != sl[0]) mlwBO[1] = MAX(mlwBO[1], lw[ll]);
        if (sl[0] == 0) ps[0] = 0.0;
        else            ps[0] = aw[0];
        for (ll = 1; ll < nLevels_; ++ll) 
                if (ll == sl[0]) ps[ll] = ps[ll - 1];
                else             ps[ll] = ps[ll - 1] + exp(lw[ll] - mlwBO[1]);
        parentSBO[1] = ps[nLevels_ - 1]; uu = runif(0, parentSBO[1]);
        for (ll = 0; ll < nLevels_; ++ll)
                if (uu <= ps[ll]) { sl[1] = ll; break; }
        parentNA[1] = aw[sl[1]];
        parentNBO[1] = exp(lw[sl[1]] - mlwBO[1]);
        
        /* do the rest of the probability computation */
        for (mlwBO[0] = R_NegInf, ll = 0; ll < nLevels_; ++ll)
                if (ll != sl[1]) mlwBO[0] = MAX(mlwBO[0], lw[ll]);
        parentNBO[0] = exp(lw[sl[0]] - mlwBO[0]);
        for (parentSBO[0] = 0.0, ll = 0; ll < nLevels_; ++ll) 
                if (ll != sl[1]) parentSBO[0] += exp(lw[ll] - mlwBO[0]);
        *prob = (parentNA[0] / parentSA) * (parentNBO[1] / parentSBO[1]) + \
                (parentNA[1] / parentSA) * (parentNBO[0] / parentSBO[0]);
        return 0;
}

int
Sampler::SCSC_TWO_NEW_BEST_BEST_probGivenLevels (double *prob, int *sl, double *logPropDens)
{
        int ll;
        double *ld = logDensities_, selTemper = moveSelectionTempers_[RC][0];
        double *lw = scratch_SLC_->logWeights_;
        /* A: All, N: Numerator, S: Sum, B: But, O: One */
        double mlwA, mlwBO[2] = { 0.0, 0.0 };
        double childNA[2] = { 0.0, 0.0 }, childSA = 0.0;
        double childNBO[2] = { 0.0, 0.0 }, childSBO[2] = { 0.0, 0.0 };

        mlwA = R_NegInf; mlwBO[0] = R_NegInf; mlwBO[1] = R_NegInf;
        for (ll = 0; ll < nLevels_; ++ll) {
                if (ll == sl[0]) {
                        lw[ll] = logPropDens[0] / selTemper;
                        mlwBO[0] = MAX(mlwBO[0], lw[ll]);
                } else if (ll == sl[1]) {
                        lw[ll] = logPropDens[1] / selTemper;
                        mlwBO[1] = MAX(mlwBO[1], lw[ll]);
                } else {
                        lw[ll] = ld[ll] / selTemper;
                        mlwBO[0] = MAX(mlwBO[0], lw[ll]);
                        mlwBO[1] = MAX(mlwBO[1], lw[ll]);
                }
                mlwA = MAX(mlwA, lw[ll]);                
        }
        childSA = 0.0; childSBO[0] = 0.0; childSBO[1] = 0.0;
        for (ll = 0; ll < nLevels_; ++ll) {
                if (ll == sl[0]) {
                        childNA[0] = exp(lw[ll] - mlwA); childSA += childNA[0];
                        childNBO[0] = exp(lw[ll] - mlwBO[0]); childSBO[0] += childNBO[0];
                } else if (ll == sl[1]) {
                        childNA[1] = exp(lw[ll] - mlwA); childSA += childNA[1];
                        childNBO[1] = exp(lw[ll] - mlwBO[1]); childSBO[1] += childNBO[1];
                } else {
                        childSA += exp(lw[ll] - mlwA);
                        childNBO[0] += exp(lw[ll] - mlwBO[0]);
                        childNBO[1] += exp(lw[ll] - mlwBO[1]);
                }
        }
        *prob = (childNA[0] / childSA) * (childNBO[1] / childSBO[1]) + \
                (childNA[1] / childSA) * (childNBO[0] / childSBO[0]);
        return 0;
}

int
Sampler::SCSC_TWO_NEW_propNew (Cluster **parents, Cluster **children)
{
        int ii, jj, nIntersections, parentsNClusters[N_SELECTIONS];
        int clustersSampled[N_SELECTIONS][N_SUB_CLUSTERS];
        int nn1, *iarr1, nn2, *iarr2;
        
        for (ii = 0; ii < N_SELECTIONS; ++ii) {
                parentsNClusters[ii] = parents[ii]->getNClusters( );
                // SCSC_TWO_NEW cannot be performed
                if (parentsNClusters[ii] < 2) return 1;
        }

        nIntersections = 0;
        for (ii = 0; ii < parentsNClusters[0]; ++ii) {
                nn1 = parents[0]->getClusterFreq(ii);
                iarr1 = parents[0]->getClusterMembersData(ii);
                for (jj = 0; jj < parentsNClusters[1]; ++jj) {
                        nn2 = parents[1]->getClusterFreq(jj);
                        iarr2 = parents[1]->getClusterMembersData(jj);
                        bool tmp = utils_are_disjoint_iarray_iarray(nn1, iarr1, nn2, iarr2);
                        if (tmp == false) {
                                boolMat_[ii][jj] = true;
                                ++nIntersections;
                        }
                        else boolMat_[ii][jj] = false;
                }
        }
        // sample the first intersection at random
        utils_sample_bool_mat_cell(parentsNClusters[0], parentsNClusters[1],
                                   boolMat_, nIntersections,
                                   clustersSampled[0], clustersSampled[1]);
        // "remove" the row
        ii = clustersSampled[0][0];
        for (jj = 0; jj < parentsNClusters[1]; ++jj) boolMat_[ii][jj] = false;
        // "remove" the column
        jj = clustersSampled[1][0];
        for (ii = 0; ii < parentsNClusters[0]; ++ii) boolMat_[ii][jj] = false;
        // count remaining cells
        nIntersections = 0;
        for (ii = 0; ii < parentsNClusters[0]; ++ii)
                for (jj = 0; jj < parentsNClusters[1]; ++jj)
                        if (boolMat_[ii][jj] == true) ++nIntersections;
        // SCSC_TWO_NEW cannot be performed
        if (nIntersections == 0) return 1;
        
        // sample the second intersection at random from the rest
        utils_sample_bool_mat_cell(parentsNClusters[0], parentsNClusters[1],
                                   boolMat_, nIntersections,
                                   clustersSampled[0] + 1, clustersSampled[1] + 1);

        for (ii = 0; ii < N_SUB_CLUSTERS; ++ii) {
                int kk = clustersSampled[0][ii], ll = clustersSampled[1][ii];
                nn1 = parents[0]->getClusterFreq(kk);
                iarr1 = parents[0]->getClusterMembersData(kk);
                nn2 = parents[1]->getClusterFreq(ll);
                iarr2 = parents[1]->getClusterMembersData(ll);
                utils_intersect_iarray_iarray(nn1, iarr1, nn2, iarr2,
                                              nSCXX_intersection_ + ii,
                                              SCXX_intersection_[ii]);

                PHONY(Rprintf("SCSC_TWO_NEW_propNew: the sampled clusters: %d, %d\n",
                              kk, ll);
//                       parents[0]->printCluster(kk);
//                       parents[1]->printCluster(ll);
                      Rprintf("SCXX_intersection_ [length: %d]:\n", nSCXX_intersection_[ii]);
                      utils_iarray_print(SCXX_intersection_[ii], nSCXX_intersection_[ii], ", "););

                for (jj = 0; jj < N_SELECTIONS; ++jj) {
                        int label = parents[jj]->getClusterLabel(clustersSampled[jj][1 - ii]);
                        children[jj]->setObjLabelsSliceFromOne(nSCXX_intersection_[ii],
                                                               SCXX_intersection_[ii],
                                                               label);
                }
        }        
        return 0;
}

int 
Sampler::SCSC_TWO_NEW_move (void)
{
        int ii, sl[N_SELECTIONS];
        double *ld = logDensities_, *tl = temperLadder_;
        Cluster *parents[N_SELECTIONS], *children[N_SELECTIONS];
        SEXP SEXPProps[N_SELECTIONS];
        double parentsToChildrenProb, childrenToParentsProb, logPropDens[2];
        double sum, alpha;

        (this->*SCSC_TWO_NEW_probAndLevelsFP_)(&parentsToChildrenProb, sl);
        for (ii = 0; ii < N_SELECTIONS; ++ii) {
                parents[ii] = clusterCurrDraws_[sl[ii]];
                children[ii] = clusterPropDraws_[sl[ii]];
                parents[ii]->tabulate( );
                children[ii]->copy(*(parents[ii]));
                SEXPProps[ii] = VECTOR_ELT(SEXPPropDraws_, sl[ii]);
        }
        int err = SCSC_TWO_NEW_propNew(parents, children);
        if (err) {
                PHONY(Rprintf("---> Could not do SCSC_TWO_NEW_propNew\n"););
                ++SCSC_TWO_NEW_notPossible_; return 1;
        }
        for (sum = 0.0, ii = 0; ii < 2; ++ii) {
                logPropDens[ii] = (this->*logTarDensFP_)(SEXPProps[ii]);
                sum += (logPropDens[ii] - ld[sl[ii]]) / tl[sl[ii]];
        }
        (this->*SCSC_TWO_NEW_probGivenLevelsFP_)(&childrenToParentsProb, sl, logPropDens);
        alpha = exp(sum) * (childrenToParentsProb / parentsToChildrenProb);
        alpha = MIN(1.0, alpha);
        ProposalCounter *pc = movePropCtrs_[SCSC_TWO_NEW][sl[0]][sl[1]];
        pc->proposed_ += 1;
        if (verboseLevel_ >= 100)
                Rprintf("SCSC_TWO_NEW: levels: %d, %d | iter: %d | " \
                        "alpha: %5.4g\n", sl[0], sl[1], thisIter_, alpha);
        /* the acceptance rejection step */
        if (runif(0, 1) <= alpha) {
                pc->accepted_ += 1;
                if (verboseLevel_ >= 10)
                        Rprintf("SCSC_TWO_NEW: levels: %d, %d | iter: %d | " \
                                "alpha: %5.4g [*** accepted]\n",
                                sl[0], sl[1], thisIter_, alpha);
                PHONY(Rprintf("SCSC_TWO_NEW changes -fitness values as:\n");
                      for (ii = 0; ii < N_SELECTIONS; ++ii)
                      Rprintf("%d: %g ---> %g\n", sl[ii], ld[sl[ii]], logPropDens[ii]);
                      Rprintf("SCSC_TWO_NEW SCXX_intersection:\n");
                      for (ii = 0; ii < N_SUB_CLUSTERS; ++ii) {
                              PRINT_STUB_INT(nSCXX_intersection_[ii]);
                              utils_iarray_print(SCXX_intersection_[ii], nSCXX_intersection_[ii], ", ");
                      });
                
                for (ii = 0; ii < N_SELECTIONS; ++ii) {
                        SEXP SEXPTmp = VECTOR_ELT(SEXPCurrDraws_, sl[ii]);
                        SET_VECTOR_ELT(SEXPCurrDraws_, sl[ii],
                                       VECTOR_ELT(SEXPPropDraws_, sl[ii]));
                        SET_VECTOR_ELT(SEXPPropDraws_, sl[ii], SEXPTmp);
                        SWAP(Cluster *, clusterCurrDraws_[sl[ii]], clusterPropDraws_[sl[ii]]);
                        ld[sl[ii]] = logPropDens[ii];
                        clusterCurrDraws_[sl[ii]]->tabulate( );
                }
        }        
        return 0;
}


/*****************************************************************************/
// The private functions: the SCRC move and its helpers.
/*****************************************************************************/

int
Sampler::SCRC_RANDOM_SIZE_sizeComputation (int *sizeIn, int *sizeOut,
                                           double *logDensRatioSize)
{
        int sumSize = sizeIn[0] + sizeIn[1];
        sizeOut[0] = (int) floor(runif(1, sumSize));
        sizeOut[1] = sumSize - sizeOut[0];
        if ((sizeOut[0] <= 0) || (sizeOut[1] <= 0))
                RAISE("should not reach here");                
        *logDensRatioSize = 0.0;
        for (int ii = 0; ii < N_SUB_CLUSTERS; ++ii)
                *logDensRatioSize += \
                        (lgamma(sizeIn[ii] + 1) - lgamma(sizeOut[ii] + 1));
        return 0;
}

int
Sampler::SCRC_SAME_SIZE_sizeComputation (int *sizeIn, int *sizeOut,
                                         double *logDensRatioSize)
{
        for (int ii = 0; ii < N_SUB_CLUSTERS; ++ii) sizeOut[ii] = sizeIn[ii];
        *logDensRatioSize = 0.0;
        return 0;
}

bool 
Sampler::SCRC_isParent (Cluster *plausibleParent, Cluster *nonSurvParent,
                        Cluster *modChild)
{
        int ii, nObjs = plausibleParent->getNObjs( );
        int *plausibleParentLabels = plausibleParent->getObjLabelsData( );
        int *nonSurvParentLables = nonSurvParent->getObjLabelsData( );
        int *modChildLabels = modChild->getObjLabelsData( );
        int count = 0, nPositive = 0, nNegative = 0, parentLabel = 0;
        
        if (nObjs != sampDim_)
                RAISE("we are in trouble");
        for (ii = 0; ii < nObjs; ++ii) {
                labelDiffs_[ii] = nonSurvParentLables[ii] - modChildLabels[ii];
                if (labelDiffs_[ii] == 0) continue;

                if (labelDiffs_[ii] > 0) {
                        labelDiffsPositive_[nPositive] = labelDiffs_[ii]; ++nPositive;
                } else if (labelDiffs_[ii] < 0) {
                        labelDiffsNegative_[nNegative] = labelDiffs_[ii]; ++nNegative;
                }
                
                if (count == 0) {
                        parentLabel = plausibleParentLabels[ii]; ++count;
                } else {
                        /*
                         * Plausible parent should have the same label
                         * responsible for non-zero label differences.
                         */
                        if (parentLabel != plausibleParentLabels[ii]) return false;
                }
        }

        /*
         * Note: We cannot have deleted one sub-cluster intersection
         * i.e. cannot choose h_1 = 0. The first two cases below takes
         * care of that. The last case arises when we end up proposing
         * the mod_child to be same as the child.
         */
        if ((nPositive == 0) && (nNegative > 0))       return false;
        else if ((nPositive > 0) && (nNegative == 0))  return false;
        else if ((nPositive == 0) || (nNegative == 0)) return true;
        
        qsort(labelDiffsPositive_, nPositive, sizeof(int), intcmp);
        qsort(labelDiffsNegative_, nNegative, sizeof(int), intcmp);
        
        int commonVal = labelDiffsPositive_[0];
        for (ii = 1; ii < nPositive; ++ii)
                if (labelDiffsPositive_[ii] != commonVal) return false;

        for (ii = 0; ii < nNegative; ++ii)
                if (labelDiffsNegative_[ii] != -commonVal) return false;        
        return true;
}

int
Sampler::SCRC_propNew (Cluster *survParent, Cluster *nonSurvParent,
                       Cluster *modChild, double *logDensRatioSize)
{
        int ii, jj, survParentNClusters, nonSurvParentNClusters;
        int survClusterSampled, nonSurvClustersSampled[N_SUB_CLUSTERS];
        int nn, nn1, *iarr1, nn2, *iarr2;
        int sizeIn[N_SUB_CLUSTERS], sizeOut[N_SUB_CLUSTERS];
        double incr, sum, uu;
        
        survParentNClusters = survParent->getNClusters( );
        nonSurvParentNClusters = nonSurvParent->getNClusters( );
        // SCRC cannot be performed
        if (nonSurvParentNClusters < 2) return 1;
        
        for (nSCXX_rows2OrMore_ = 0, ii = 0; ii < survParentNClusters; ++ii) {
                nn1 = survParent->getClusterFreq(ii);
                iarr1 = survParent->getClusterMembersData(ii);
                int count = 0;
                for (jj = 0; jj < nonSurvParentNClusters; ++jj) {
                        nn2 = nonSurvParent->getClusterFreq(jj);
                        iarr2 = nonSurvParent->getClusterMembersData(jj);
                        bool tmp = utils_are_disjoint_iarray_iarray(nn1, iarr1, nn2, iarr2);
                        if (tmp == false) {
                                SCXX_nonEmptyInRow_[ii][count] = jj; ++count;
                        }
                }
                nSCXX_nonEmptyInRow_[ii] = count;
                if (nSCXX_nonEmptyInRow_[ii] >= 2) {
                        SCXX_rows2OrMore_[nSCXX_rows2OrMore_] = ii;
                        ++nSCXX_rows2OrMore_;
                }                
        }
        
        PHONY(PRINT_STUB_INT(nSCXX_rows2OrMore_);
              utils_iarray_print(SCXX_rows2OrMore_, nSCXX_rows2OrMore_, ", "););

        // SCRC cannot be performed
        if (nSCXX_rows2OrMore_ == 0) return 1;
        
        // choose the survParent sampled cluster
        ii = (int) floor(nSCXX_rows2OrMore_ * runif(0, 1));
        survClusterSampled = SCXX_rows2OrMore_[ii];

        // choose the nonSurvParent sampled clusters
        nn = nSCXX_nonEmptyInRow_[survClusterSampled];
        jj = (int) floor(nn * runif(0, 1));
        nonSurvClustersSampled[0] = SCXX_nonEmptyInRow_[survClusterSampled][jj];
        
        incr = 1.0 / (nn - 1); sum = 0.0; uu = runif(0, 1);
        for (jj = 0; jj < nn; ++jj) {
                int kk = SCXX_nonEmptyInRow_[survClusterSampled][jj];
                if (kk == nonSurvClustersSampled[0]) continue;
                sum += incr;
                if (uu <= sum) {
                        nonSurvClustersSampled[1] = kk; break;
                }
        }
        
        // form the union of the two sub-cluster intersections
        nSCRC_union_ = 0;
        nn1 = survParent->getClusterFreq(survClusterSampled);
        iarr1 = survParent->getClusterMembersData(survClusterSampled);
        for (ii = 0; ii < N_SUB_CLUSTERS; ++ii) {
                nn2 = nonSurvParent->getClusterFreq(nonSurvClustersSampled[ii]);
                iarr2 = nonSurvParent->getClusterMembersData(nonSurvClustersSampled[ii]);
                utils_intersect_iarray_iarray(nn1, iarr1, nn2, iarr2,
                                              nSCXX_intersection_ + ii,
                                              SCXX_intersection_[ii]);
                sizeIn[ii] = nSCXX_intersection_[ii];

                for (jj = 0; jj < nSCXX_intersection_[ii]; ++jj)
                        SCRC_union_[nSCRC_union_ + jj] = SCXX_intersection_[ii][jj];
                nSCRC_union_ += nSCXX_intersection_[ii];
                
                PHONY(Rprintf("SCRC_propNew: the sampled clusters: %d, %d\n",
                              survClusterSampled, nonSurvClustersSampled[ii]);
                      Rprintf("SCXX_intersection_ [length: %d]:\n", nSCXX_intersection_[ii]);
                      survParent->prettyPrint( );
                      nonSurvParent->prettyPrint( );
                      utils_iarray_print(SCXX_intersection_[ii], nSCXX_intersection_[ii], ", "););
        }

        // take a sub-sample from the union
        (this->*SCRC_sizeComputationFP_)(sizeIn, sizeOut, logDensRatioSize);
        int minOutSize = MIN(sizeOut[0], sizeOut[1]);
        utils_sample_indices_SRSWOR(nSCRC_union_, minOutSize,
                                    SCRC_isAlreadySampled_,
                                    SCRC_subSampleIndices_);

        // label all the members of the union
        int label = nonSurvParent->getClusterLabel(nonSurvClustersSampled[0]);
        modChild->setObjLabelsSliceFromOne(nSCRC_union_, SCRC_union_, label);
        // label the smaller subset differently
        label = nonSurvParent->getClusterLabel(nonSurvClustersSampled[1]);
        for (jj = 0; jj < minOutSize; ++jj) {
                int kk = SCRC_union_[SCRC_subSampleIndices_[jj]];
                modChild->setObjLabel(kk, label);
        }
        return 0;
}

int 
Sampler::SCRC_move (void)
{
        int ii, sl[2 * N_SELECTIONS];
        double *ld = logDensities_, *tl = temperLadder_;
        Cluster *survParent, *nonSurvParent, *modChild;
        double parentsToChildrenProb, childrenToParentsProb, logPropDens;
        double probsPTC[2], probsCTP[2], logDensRatioSize, sum, alpha;
        bool isParent;
        
        SCSC_ONE_NEW_BEST_BEST_probAndLevels(probsPTC, sl);
        parentsToChildrenProb = probsPTC[0];        
        survParent = clusterCurrDraws_[sl[SURV_PARENT_POS]];
        survParent->tabulate( );
        nonSurvParent = clusterCurrDraws_[sl[NON_SURV_PARENT_POS]];
        nonSurvParent->tabulate( );
        modChild = clusterPropDraws_[sl[NON_SURV_PARENT_POS]];
        modChild->copy(*nonSurvParent);
        int err = SCRC_propNew(survParent, nonSurvParent, modChild,
                               &logDensRatioSize);
        if (err) {                
                PHONY(Rprintf("---> Could not do SCRC_propNew\n"););
                ++SCRC_notPossible_; return 1;
        }
        SEXP SEXPModChild = VECTOR_ELT(SEXPPropDraws_, sl[NON_SURV_PARENT_POS]);
        logPropDens = (this->*logTarDensFP_)(SEXPModChild);
        sum = (logPropDens - ld[sl[NON_SURV_PARENT_POS]]) / tl[sl[NON_SURV_PARENT_POS]];
        SCSC_ONE_NEW_BEST_BEST_probGivenLevels(probsCTP, sl, &logPropDens);
        childrenToParentsProb = probsCTP[0];

        // checking the 'isParent' stuff
        if ((sl[NON_SURV_PARENT_POS] > 0) &&
            (sl[NON_SURV_PARENT_POS] < (nLevels_ - 1))) {
                Cluster *plausibleParent = NULL;
                
                if (sl[SURV_PARENT_POS] == (sl[NON_SURV_PARENT_POS] - 1))
                        plausibleParent = clusterCurrDraws_[sl[NON_SURV_PARENT_POS] + 1];
                else if (sl[SURV_PARENT_POS] == (sl[NON_SURV_PARENT_POS] + 1))
                        plausibleParent = clusterCurrDraws_[sl[NON_SURV_PARENT_POS] - 1];
                else
                        RAISE("should not reach here");
                plausibleParent->tabulate( );
                isParent = SCRC_isParent(plausibleParent, nonSurvParent, modChild);
                if (isParent == true) {
                        PHONY(Rprintf("---> SCRC isParent is true\n"););
                        parentsToChildrenProb += probsPTC[1];
                        childrenToParentsProb += probsCTP[1];
                }
        }

        sum += logDensRatioSize;
        alpha = exp(sum) * (childrenToParentsProb / parentsToChildrenProb);
        alpha = MIN(1.0, alpha);
        ProposalCounter *pc = movePropCtrs_[SCRC][sl[0]][sl[1]];
        pc->proposed_ += 1;
        if (verboseLevel_ >= 100)
                Rprintf("SCRC: levels: %d, %d | iter: %d | " \
                        "alpha: %5.4g\n", sl[0], sl[1], thisIter_, alpha);
        /* the acceptance rejection step */
        if (runif(0, 1) <= alpha) {
                pc->accepted_ += 1;
                if (verboseLevel_ >= 10)
                        Rprintf("SCRC: levels: %d, %d | iter: %d | " \
                                "alpha: %5.4g [*** accepted]\n",
                                sl[0], sl[1], thisIter_, alpha);
                PHONY(Rprintf("SCRC changes -fitness values as:\n");
                      ii = NON_SURV_PARENT_POS;
                      Rprintf("%d: %g ---> %g\n", sl[ii], ld[sl[ii]], logPropDens);
                      Rprintf("SCRC SCXX_intersection:\n");
                      for (ii = 0; ii < N_SUB_CLUSTERS; ++ii) {
                              PRINT_STUB_INT(nSCXX_intersection_[ii]);
                              utils_iarray_print(SCXX_intersection_[ii], nSCXX_intersection_[ii], ", ");
                      });
                
                // only swap things for the non-surviving parent.
                ii = NON_SURV_PARENT_POS;
                SEXP SEXPTmp = VECTOR_ELT(SEXPCurrDraws_, sl[ii]);
                SET_VECTOR_ELT(SEXPCurrDraws_, sl[ii],
                               VECTOR_ELT(SEXPPropDraws_, sl[ii]));
                SET_VECTOR_ELT(SEXPPropDraws_, sl[ii], SEXPTmp);
                SWAP(Cluster *, clusterCurrDraws_[sl[ii]], clusterPropDraws_[sl[ii]]);
                ld[sl[ii]] = logPropDens;
                clusterCurrDraws_[sl[ii]]->tabulate( );
        }        
        return 0;
}


/*****************************************************************************/
// The private functions: the exchange family of moves.
/*****************************************************************************/

int 
Sampler::RE_move (void)
{
        int sl[N_SELECTIONS];

        /* choose the levels to exchange */
        sl[0] = (int) floor(nLevels_ * runif(0, 1));
        if (sl[0] == 0)                   sl[1] = sl[0] + 1;
        else if (sl[0] == (nLevels_ - 1)) sl[1] = sl[0] - 1;
        else {
                if (runif(0, 1) <= 0.5) sl[1] = sl[0] + 1;
                else                    sl[1] = sl[0] - 1;
        }
        exchange(sl, movePropCtrs_[RE][sl[0]][sl[1]], "RE");
        return 0;
}

int 
Sampler::BCE_move (void)
{
        int sl[N_SELECTIONS];
        double *ld = logDensities_, *itl = invTemperLadder_;
        double *lw = scratch_SLC_->logWeights_;
        double selTemper = moveSelectionTempers_[BCE][0];
        double mlw, denom, CC, logR1, logR2, diff, alpha;
        
        sl[0] = thisStep_ + 2; mlw = R_NegInf; 
        for (int ll = 0; ll < sl[0]; ++ll) {
                lw[ll] = ld[ll] / selTemper; mlw = MAX(mlw, lw[ll]);
        }
        scratch_SLC_->maxLogWeights_ = mlw; 
        scratch_SLC_->endLevel_ = sl[0];
        sampleWithDetails( ); 
        sl[1] = scratch_SLC_->samp_;
        denom = exp(ld[sl[0]] / selTemper - mlw);
        CC = scratch_SLC_->sum_ - scratch_SLC_->numerator_;
        logR2 = log(scratch_SLC_->sum_ / (denom + CC));
        diff = (ld[sl[0]] - ld[sl[1]]);
        logR2 +=  diff / selTemper;
        logR1 = diff * (itl[sl[0]] - itl[sl[1]]);
        alpha = exp(logR1 + logR2); alpha = MIN(1.0, alpha);
        ProposalCounter *pc = movePropCtrs_[BCE][sl[0]][sl[1]];
        pc->proposed_ += 1;
        exchangeGivenProb(sl, pc, "BCE", alpha);
        return 0;
}

int 
Sampler::BIRE_move (void)
{
        int sl[N_SELECTIONS];
        double *ld = logDensities_, *tl = temperLadder_;
        double *lw = scratch_SLC_->logWeights_, mlw, denom, CC, alpha;

        sl[0] = thisStep_ + 2; mlw = R_NegInf;
        for (int ll = 0; ll < sl[0]; ++ll) {
                lw[ll] = (ld[ll] / tl[sl[0]]) - (ld[ll] / tl[ll]);
                mlw = MAX(mlw, lw[ll]);
        }
        scratch_SLC_->maxLogWeights_ = mlw;
        scratch_SLC_->endLevel_ = sl[0];
        sampleWithDetails( );
        sl[1] = scratch_SLC_->samp_;
        denom = exp((ld[sl[0]] / tl[sl[0]]) - (ld[sl[0]] / tl[sl[1]]) - mlw);
        CC = scratch_SLC_->sum_ - scratch_SLC_->numerator_;
        alpha = scratch_SLC_->sum_ / (denom + CC); alpha = MIN(1.0, alpha);
        ProposalCounter *pc = movePropCtrs_[BIRE][sl[0]][sl[1]];
        pc->proposed_ += 1;
        exchangeGivenProb(sl, pc, "BIRE", alpha);
        return 0;
}

int 
Sampler::BSE_move (void)
{
        int sl[N_SELECTIONS];
        double *ld = logDensities_, *tl = temperLadder_;
        double *lw = scratch_SLC_->logWeights_, mlw, denom, CC, alpha;

        sl[0] = thisStep_ + 2; mlw = R_NegInf;
        for (int ll = 0; ll < sl[0]; ++ll) {
                lw[ll] = (ld[ll] / tl[sl[0]]) + (ld[sl[0]] / tl[ll]);
                mlw = MAX(mlw, lw[ll]);
        }
        scratch_SLC_->maxLogWeights_ = mlw;
        scratch_SLC_->endLevel_ = sl[0];
        sampleWithDetails( );
        sl[1] = scratch_SLC_->samp_;
        denom = exp((ld[sl[0]] / tl[sl[0]]) + (ld[sl[1]] / tl[sl[1]]) - mlw);
        CC = scratch_SLC_->sum_ - scratch_SLC_->numerator_;
        alpha = scratch_SLC_->sum_ / (denom + CC); alpha = MIN(1.0, alpha);
        ProposalCounter *pc = movePropCtrs_[BSE][sl[0]][sl[1]];
        pc->proposed_ += 1;
        exchangeGivenProb(sl, pc, "BSE", alpha);
        return 0;
}

int
Sampler::cyclicShiftDraws (int ladderLength, int selLength)
{
        int jj, shiftBy = selLength + 1, rem;
        Cluster **ccds = clusterCurrDrawsStore_, **ccd = clusterCurrDraws_;
        SEXP *cds = SEXPCurrDrawsStore_, cd = SEXPCurrDraws_;
        double *lds = logDensitiesStore_, *ld = logDensities_;

        for (jj = 0; jj < ladderLength; ++jj) {
                ccds[jj] = ccd[jj];
                cds[jj] = VECTOR_ELT(cd, jj);
                lds[jj] = ld[jj];
        }
        for (jj = 0; jj < ladderLength; ++jj) {
                rem = (jj + shiftBy) % ladderLength;
                ccd[jj] = ccds[rem];
                SET_VECTOR_ELT(cd, jj, cds[rem]);
                ld[jj] = lds[rem];
        }
        return 0;
}

int 
Sampler::CE_move (void)
{
        int ladderLength, selLength, shiftBy, rem;
        double *ld = logDensities_, *tl = temperLadder_;
        double *lw = scratch_SLC_->logWeights_, mlw, sum, logDenom, CC, alpha;

        ladderLength = thisStep_ + 3; mlw = R_NegInf;
        for (int ll = 0; ll < (ladderLength - 1); ++ll) {
                shiftBy = ll + 1; lw[ll] = 0.0;
                for (int jj = 0; jj < ladderLength; ++jj) {
                        rem = (jj + shiftBy) % ladderLength;
                        lw[ll] += ld[rem] / tl[jj];
                }
                mlw = MAX(mlw, lw[ll]);
        }
        scratch_SLC_->maxLogWeights_ = mlw;
        scratch_SLC_->endLevel_ = ladderLength - 1;
        sampleWithDetails( );
        selLength = scratch_SLC_->samp_;
        sum = scratch_SLC_->sum_;
        CC = sum - scratch_SLC_->numerator_;
        logDenom = 0.0;
        for (int jj = 0; jj < ladderLength; ++jj)
                logDenom += ld[jj] / tl[jj];
        alpha = sum / (exp(logDenom - mlw) + CC); alpha = MIN(alpha, 1.0);
        ProposalCounter *pc = movePropCtrs_[CE][ladderLength - 1][selLength];
        pc->proposed_ += 1;
        if (runif(0, 1) <= alpha) {
                pc->accepted_ += 1;
                cyclicShiftDraws(ladderLength, selLength);
        }
        return 0;
}


/*****************************************************************************/
// The private functions: the utility functions.
/*****************************************************************************/
/*
 * The following has been borrowed from src/main/deriv.c from the R
 * distribution.
 */

SEXP
Sampler::lang5 (SEXP s, SEXP t, SEXP u, SEXP v, SEXP w)
{
        PROTECT(s);
        s = LCONS(s, list4(t, u, v, w));
        UNPROTECT(1);
        return s;
}

/*
 * The following has been taken from the "Manual" section of R
 * website, it gets the list element named str, or returns NULL, if
 * not found. It has been modified a little bit.
 */
SEXP
Sampler::getListElement (SEXP list, char *str)
{
        SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
        int ii, found = 0, nn = length(list);
        
        for (ii = 0; ii < nn; ii++) {
                if(strcmp(CHAR(STRING_ELT(names, ii)), str) == 0) {
                        elmt = VECTOR_ELT(list, ii); ++found; break;
                }
        }
        if (found == 0) {
                char *errMsg1, errMsg2[MAX_LINE_LENGTH];

                errMsg1 = (char *) R_alloc(MAX_LINE_LENGTH, sizeof(errMsg1));
                for (ii = 0; ii < nn; ii++) {
                        errMsg1 = strcat(errMsg1, CHAR(STRING_ELT(names, ii)));
                        errMsg1 = strcat(errMsg1,
                                         ((ii == (nn - 1)) ? "" : ", "));
                }
                        
                sprintf(errMsg2,
                        "No element called \"%s\" found in the SEXP list, " \
                        "check your C code; given elements are:\n%s\n",
                        str, errMsg1);
                error(errMsg2);
        }
        return elmt;
}

int
Sampler::gatherTimeDetails (void)
{
        SEXP SEXPTmp;
        double *doublesTmp;
        
        PROTECT(SEXPTmp = eval(procTimeFuncCall_, procTimeFuncEnv_));
        doublesTmp = REAL(SEXPTmp);
        timeDetails_->usr = doublesTmp[0];
        timeDetails_->sys = doublesTmp[1];
        UNPROTECT(1);
        return 0;
}

Sampler::ProposalCounter ***
Sampler::PCMatNew (int dim1, int dim2)
{
        ProposalCounter ***pppc;
        int ii, jj;
        
        pppc = (ProposalCounter ***) R_alloc(dim1, sizeof(ProposalCounter **));
        for (ii = 0; ii < dim1; ++ii) {
                pppc[ii] = (ProposalCounter **) R_alloc(dim2, sizeof(ProposalCounter *));
                for (jj = 0; jj < dim2; ++jj) {
                        pppc[ii][jj] = (ProposalCounter *) R_alloc(1, sizeof(struct ProposalCounter));
                        (pppc[ii][jj])->accepted_ = 0.0;
                        (pppc[ii][jj])->proposed_ = 0.0;
                }
        }
        return pppc;
}

enum Sampler::SelectionCode
Sampler::getSelectionCode (char const *codeName)
{
        SelectionCode code;
        
        if (strcmp(codeName, "random") == 0) code = RANDOM;
        else if (strcmp(codeName, "best") == 0) code = BEST;
        else code = WORST;
        return code;
}

double 
Sampler::logTarDensFuncUserRfunc (SEXP argDraw)
{
        ArgsList1 *al = logTarDensArgsList_;
        SEXP SEXPTmp;
        double res;
            
        SET_VECTOR_ELT(al->argsList_, al->posDraw_, argDraw);
        
        SETCADR(doCallFuncCall_, logTarDensFunc_);
        SETCADDR(doCallFuncCall_, al->argsList_);
        SETCADDDR(doCallFuncCall_, dotsList_);        
        PROTECT(SEXPTmp = eval(doCallFuncCall_, doCallFuncEnv_));
        res = REAL(SEXPTmp)[0];
        UNPROTECT(1);
        return res;
}

inline double
Sampler::getNeighbourProb (int jjGiven, int ii)
{
        if ((ii == 0) || (ii == (nLevels_ - 1))) return 1.0;
        return 0.5;
}

int
Sampler::exchangeGivenProb (int *sl, ProposalCounter *pc, char *moveName, 
                            double prob)
{        
        if (runif(0, 1) > prob) return 0;
        pc->accepted_ += 1; 
        if (verboseLevel_ >= 10) 
                Rprintf("%s: levels: %d, %d | iter: %d | " \
                        "alpha: %5.4g [*** accepted]\n",
                        moveName, sl[0], sl[1], thisIter_, prob);
        SEXP SEXPTmp = VECTOR_ELT(SEXPCurrDraws_, sl[0]);
        SET_VECTOR_ELT(SEXPCurrDraws_, sl[0], 
                       VECTOR_ELT(SEXPCurrDraws_, sl[1]));
        SET_VECTOR_ELT(SEXPCurrDraws_, sl[1], SEXPTmp);
        SWAP(Cluster *, clusterCurrDraws_[sl[0]], clusterCurrDraws_[sl[1]]);
        SWAP(double, logDensities_[sl[0]], logDensities_[sl[1]]);
        return 0;
}

int 
Sampler::exchange (int *sl, ProposalCounter *pc, char *moveName)
{
        double *ld = logDensities_, *itl = invTemperLadder_, alpha;
        
        pc->proposed_ += 1;
        alpha = exp((ld[sl[0]] - ld[sl[1]]) * (itl[sl[1]] - itl[sl[0]]));
        alpha = MIN(1.0, alpha);
        if (verboseLevel_ >= 100) 
                Rprintf("%s: levels: %d, %d | iter: %d | " \
                        "alpha: %5.4g\n", moveName, sl[0], sl[1], thisIter_, alpha);
        exchangeGivenProb(sl, pc, moveName, alpha);
        return 0;
}

int
Sampler::initArgsList1 (Sampler::ArgsList1 *al)
{
        int nComps = 0, comp = 0, nProtected = 0;
        SEXP names;

        al->posDraw_ = nComps++;

        PROTECT(al->argsList_ = allocVector(VECSXP, nComps)); 
        CountNProtected::incrementNProtected( );                
        PROTECT(names         = allocVector(STRSXP, nComps)); ++nProtected;

        for (comp = 0; comp < nComps; ++comp)
                SET_VECTOR_ELT(al->argsList_, comp, R_NilValue);
        comp = 0;
        SET_STRING_ELT(names, comp, mkChar("draw")); ++comp;

        setAttrib(al->argsList_, R_NamesSymbol, names);
        UNPROTECT(nProtected);
        return 0;
}

int
Sampler::init (void)
{
        if (logTarDensFunc_ != R_NilValue)
                logTarDensFP_ = &Sampler::logTarDensFuncUserRfunc;
        
        switch (MHMoveCode_) {
                case SPLIT_MERGE:
                        moveObjs_[MH]->FP_ = &Sampler::MH_SPLIT_MERGE_move;
                        break;
                        
                case  GIBBS:
                        moveObjs_[MH]->FP_ = &Sampler::MH_GIBBS_move;
                        break;
                        
                default:
                        RAISE("invalid MHMoveCode code");
        }
        if (nLevels_ == 1) {
                oneIterFP_ = &Sampler::oneIterWithOneLevel;
        } else if (nLevels_ == 2) {
                oneIterFP_ = &Sampler::oneIterWithTwoLevels;
                /*
                 * In this case the BEST-BEST case yields to proposal
                 * probability cacellations and hence we just use the
                 * RANDOM-RANDOM case a proxy for the BEST-BEST case.
                 */
                RC_probAndLevelsFP_             = &Sampler::RC_RANDOM_RANDOM_probAndLevels;
                RC_probGivenLevelsFP_           = &Sampler::RC_RANDOM_RANDOM_probGivenLevels;
                SCSC_ONE_NEW_probAndLevelsFP_   = &Sampler::SCSC_ONE_NEW_RANDOM_RANDOM_probAndLevels;
                SCSC_ONE_NEW_probGivenLevelsFP_ = &Sampler::SCSC_ONE_NEW_RANDOM_RANDOM_probGivenLevels;
                SCSC_TWO_NEW_probAndLevelsFP_   = &Sampler::SCSC_TWO_NEW_RANDOM_RANDOM_probAndLevels;
                SCSC_TWO_NEW_probGivenLevelsFP_ = &Sampler::SCSC_TWO_NEW_RANDOM_RANDOM_probGivenLevels;
        } else {
                oneIterFP_ = &Sampler::oneIterWithManyLevels;
                if (moveSelectionCodes_[RC][0] == RANDOM) {
                        RC_probAndLevelsFP_ = &Sampler::RC_RANDOM_RANDOM_probAndLevels;
                        RC_probGivenLevelsFP_ = &Sampler::RC_RANDOM_RANDOM_probGivenLevels;
                }
                else {
                        RC_probAndLevelsFP_ = &Sampler::RC_BEST_BEST_probAndLevels;
                        RC_probGivenLevelsFP_ = &Sampler::RC_BEST_BEST_probGivenLevels;
                }
                if (moveSelectionCodes_[SCSC_TWO_NEW][0] == RANDOM) {
                        SCSC_TWO_NEW_probAndLevelsFP_ = &Sampler::SCSC_TWO_NEW_RANDOM_RANDOM_probAndLevels;
                        SCSC_TWO_NEW_probGivenLevelsFP_ = &Sampler::SCSC_TWO_NEW_RANDOM_RANDOM_probGivenLevels;
                }
                else {
                        SCSC_TWO_NEW_probAndLevelsFP_ = &Sampler::SCSC_TWO_NEW_BEST_BEST_probAndLevels;
                        SCSC_TWO_NEW_probGivenLevelsFP_ = &Sampler::SCSC_TWO_NEW_BEST_BEST_probGivenLevels;
                }
                if (moveSelectionCodes_[SCSC_ONE_NEW][0] == RANDOM) {
                        SCSC_ONE_NEW_probAndLevelsFP_ = &Sampler::SCSC_ONE_NEW_RANDOM_RANDOM_probAndLevels;
                        SCSC_ONE_NEW_probGivenLevelsFP_ = &Sampler::SCSC_ONE_NEW_RANDOM_RANDOM_probGivenLevels;
                }
                else {
                        SCSC_ONE_NEW_probAndLevelsFP_ = &Sampler::SCSC_ONE_NEW_BEST_BEST_probAndLevels;
                        SCSC_ONE_NEW_probGivenLevelsFP_ = &Sampler::SCSC_ONE_NEW_BEST_BEST_probGivenLevels;
                }
        }
        if (moveSelectionCodes_[SCRC][0] == RANDOM)
                SCRC_sizeComputationFP_ = &Sampler::SCRC_RANDOM_SIZE_sizeComputation;
        else
                SCRC_sizeComputationFP_ = &Sampler::SCRC_SAME_SIZE_sizeComputation;
        moveObjs_[RC]->FP_           = &Sampler::RC_move;
        moveObjs_[SCSC_ONE_NEW]->FP_ = &Sampler::SCSC_ONE_NEW_move;
        moveObjs_[SCSC_TWO_NEW]->FP_ = &Sampler::SCSC_TWO_NEW_move;
        moveObjs_[SCRC]->FP_         = &Sampler::SCRC_move;
        moveObjs_[RE]->FP_           = &Sampler::RE_move;
        moveObjs_[BCE]->FP_          = &Sampler::BCE_move;
        moveObjs_[BIRE]->FP_         = &Sampler::BIRE_move;
        moveObjs_[BSE]->FP_          = &Sampler::BSE_move;
        moveObjs_[CE]->FP_           = &Sampler::CE_move;
        
        for (int ii = 0; ii < nLevels_; ++ii) {
                SEXP curr = VECTOR_ELT(SEXPCurrDraws_, ii);
                logDensities_[ii] = (this->*logTarDensFP_)(curr);
                if (R_FINITE(logDensities_[ii]) == FALSE) {
                        char errMsg[MAX_LINE_LENGTH];
                        sprintf(errMsg,
                                "logTarDens evaluation for level [%d] " \
                                "gives [%f], please use different starting " \
                                "values for this level", ii, logDensities_[ii]);
                        RAISE(errMsg);
                }
        }
        if (timeInSecs_ <= 0) {
                if (saveFitness_ == TRUE)
                        registerIterFP_ = &Sampler::registerIterFixedIterWithFitness;
                else
                        registerIterFP_ = &Sampler::registerIterFixedIter;
        }
        else {
                if (saveFitness_ == TRUE)
                        registerIterFP_ = &Sampler::registerIterFixedTimeWithFitness;
                else
                        registerIterFP_ = &Sampler::registerIterFixedTime;
        }
        return 0;       
}

int
Sampler::print (void)
{
        // WRITE ME
        PRINT_STUB_INT(nIters_);
        PRINT_STUB_DOUBLE(timeInSecs_);
        PRINT_STUB_INT(verboseLevel_);
        PRINT_STUB_INT(nLevels_);
        PRINT_STUB_DARRAY(temperLadder_, nLevels_, ", ");
        PRINT_STUB_DARRAY(invTemperLadder_, nLevels_, ", ");
        PRINT_STUB_DARRAY(logDensities_, nLevels_, ", ");
        PRINT_STUB_INT(sampDim_);
        PRINT_STUB_INT(MHNBlocks_);
        PRINT_STUB_DOUBLE(MHMergeProb_);
        PRINT_STUB_DARRAY(moveProbs_, N_IMPLEMENTED_MOVES, ", ");
        PRINT_STUB_DARRAY(cumsumProbs_, N_IMPLEMENTED_MOVES, ", ");
        PRINT_STUB_IARRAY(moveNTimes_, N_IMPLEMENTED_MOVES, ", ");
        PRINT_STUB_IARRAY(levelsSaveSampFor_, nLevelsSaveSampFor_, ", ");
        for (int ii = 0; ii < nLevels_; ++ii) {
                Rprintf("This level: %d with log density: %g\n", ii, logDensities_[ii]);
                clusterCurrDraws_[ii]->tabulate( );
                clusterCurrDraws_[ii]->prettyPrint( );
        }
        PRINT_STUB_INT(SCSC_ONE_NEW_notPossible_);
        PRINT_STUB_INT(SCSC_TWO_NEW_notPossible_);
        PRINT_STUB_INT(SCRC_notPossible_);
        return 0;
}

int
Sampler::sampleWithDetails (void)
{
        double *lw = scratch_SLC_->logWeights_, *aw = scratch_SLC_->adjWeights_;
        double *ps = scratch_SLC_->partialSum_, mlw = scratch_SLC_->maxLogWeights_;
        double uu, sum;
        int ll, nn = scratch_SLC_->endLevel_;
        
        aw[0] = exp(lw[0] - mlw); ps[0] = aw[0];
        for (ll = 1; ll < nn; ++ll) {
                aw[ll] = exp(lw[ll] - mlw); ps[ll] = ps[ll - 1] + aw[ll];
        }
        sum = ps[nn - 1]; uu = runif(0, sum);
        for (ll = 0; ll < nn; ++ll) {                
                if (uu <= ps[ll]) {
                        scratch_SLC_->samp_ = ll; 
                        scratch_SLC_->prob_ = aw[ll] / sum;
                        scratch_SLC_->numerator_ = aw[ll]; 
                        scratch_SLC_->sum_ = sum;
                        break;
                }
        }
        return 0;
}

int
Sampler::sampleWithDetailsExcludingOne (int sl)
{
        double *lw = scratch_SLC_->logWeights_, *aw = scratch_SLC_->adjWeights_;
        double *ps = scratch_SLC_->partialSum_, mlw = scratch_SLC_->maxLogWeights_;
        double uu, sum;
        int ll, nn = scratch_SLC_->endLevel_;
        
        if (sl == 0) ps[0] = 0.0;
        else         aw[0] = exp(lw[0] - mlw); ps[0] = aw[0];
        for (ll = 1; ll < nn; ++ll) {
                if (sl == ll) { ps[ll] = ps[ll - 1]; continue; }
                aw[ll] = exp(lw[ll] - mlw); ps[ll] = ps[ll - 1] + aw[ll];
        }
        sum = ps[nn - 1]; uu = runif(0, sum);
        for (ll = 0; ll < nn; ++ll) {
                if ((sl != ll) && (uu <= ps[ll])) {
                        scratch_SLC_->samp_ = ll; 
                        scratch_SLC_->prob_ = aw[ll] / sum;
                        scratch_SLC_->numerator_ = aw[ll]; 
                        scratch_SLC_->sum_ = sum;
                        break;
                }
        }
        return 0;
}

int 
Sampler::oneIterWithOneLevel (SEXP notRequired)
{
        for (int ii = 0; ii < moveNTimes_[MH]; ++ii) {
                thisLevel_ = 0;
                (this->*(moveObjs_[MH]->FP_))( );
        }
        return 0;
}       

int 
Sampler::setMoveNTimesAtIter (void)
{
        double uu = runif(0, 1);
        bool found = false;
        
        for (int jj = MH; jj <= CE; ++jj) {
                if (moveProbsPositive_[jj] == true) {
                        if ((found == false) && (uu <= cumsumProbs_[jj])) {
                                moveNTimesAtIter_[jj] = moveNTimes_[jj];
                                found = true;
                        }
                        else { moveNTimesAtIter_[jj] = 0; }
                }
                else moveNTimesAtIter_[jj] = moveNTimes_[jj];   
        }
        return 0;
}

int 
Sampler::oneIterWithTwoLevels (SEXP notRequired)
{
        int ii, sl[N_SELECTIONS];
        ProposalCounter *pc;
        char *moveName;
        
        setMoveNTimesAtIter( );
        for (ii = 0; ii < moveNTimesAtIter_[MH]; ++ii) {
                for (int ll = 0; ll < nLevels_; ++ll) {
                        thisLevel_ = ll;
                        (this->*(moveObjs_[MH]->FP_))( );
                }
        }
        for (ii = 0; ii < moveNTimesAtIter_[RC]; ++ii)
                (this->*(moveObjs_[RC]->FP_))( );
        
        // WRITE ME: WHAT HAPPENS FOR SCSC AND SCRC FAMILY OF MOVES?
        
        for (int jj = RE; jj <= CE; ++jj) {
                for (ii = 0; ii < moveNTimesAtIter_[jj]; ++ii) {
                        thisStep_ = ii; sl[0] = 0; sl[1] = 1;
                        pc = movePropCtrs_[jj][sl[0]][sl[1]];
                        moveName = moveObjs_[jj]->name_;
                        exchange(sl, pc, moveName);
                }
        }
        return 0;
}

int 
Sampler::oneIterWithManyLevels (SEXP notRequired)
{
        int  ii;
        
        setMoveNTimesAtIter( );
        PHONY(PRINT_STUB_STRING("we enter");
              PRINT_STUB_IARRAY(moveNTimesAtIter_, N_IMPLEMENTED_MOVES, ","););
        for (ii = 0; ii < moveNTimesAtIter_[MH]; ++ii) {
                for (int ll = 0; ll < nLevels_; ++ll) {
                        thisLevel_ = ll; (this->*(moveObjs_[MH]->FP_))( );
                }
        }        
        for (int jj = RC; jj <= CE; ++jj) {
                for (ii = 0; ii < moveNTimesAtIter_[jj]; ++ii) {
                        thisStep_ = ii; (this->*(moveObjs_[jj]->FP_))( );
                }
        }
        PHONY(PRINT_STUB_STRING("we leave"));
        return 0;
}

int 
Sampler::registerIterFixedIter (SEXP SEXPDraws)
{
        double *dest = REAL(SEXPDraws);
        
        for (int kk = 0; kk < nLevelsSaveSampFor_; ++kk) {
                int ll = levelsSaveSampFor_[kk];
                int *src = INTEGER(VECTOR_ELT(SEXPCurrDraws_, ll));
                int mm = (kk * nIters_ * sampDim_); 
                for (int jj = 0; jj < sampDim_; ++jj) {
                        int ii = mm + thisIter_ + jj * nIters_;
                        dest[ii] = src[jj];
                }
                PHONY(utils_darray_print(src, sampDim_, ", "););
        }
        return 0;
}

int 
Sampler::registerIterFixedIterWithFitness (SEXP SEXPDraws)
{
        double *dest = REAL(SEXPDraws);
        
        for (int kk = 0; kk < nLevelsSaveSampFor_; ++kk) {
                int ll = levelsSaveSampFor_[kk];
                int *src = INTEGER(VECTOR_ELT(SEXPCurrDraws_, ll));
                int mm = (kk * nIters_ * (sampDim_ + 1)), ii, jj;                 
                for (jj = 0; jj < sampDim_; ++jj) {
                        ii = mm + thisIter_ + jj * nIters_;
                        dest[ii] = src[jj];
                }
                ii = mm + thisIter_ + jj * nIters_;
                dest[ii] = -logDensities_[ll];
                PHONY(utils_darray_print(src, sampDim_, ", "););
        }
        PHONY(PRINT_STUB_INT(thisIter_);
              utils_darray_print(logDensities_, nLevels_, ", "););
        return 0;
}

int 
Sampler::registerIterFixedTime (SEXP SEXPDraws)
{
        RAISE("not implmented yet");
        return 0;
}

int 
Sampler::registerIterFixedTimeWithFitness (SEXP SEXPDraws)
{
        RAISE("not implmented yet");
        return 0;
}

int
Sampler::makeDrawsFixedIter (SEXP SEXPDraws)
{
        int ii, iter = 0;
        bool initialEstDone = false;
        double timeStartUsr, timeEndUsr, timeToFinish;
        double timeStartSys, timeEndSys;    
        
        gatherTimeDetails( );
        timeStartUsr = timeDetails_->usr;
        timeStartSys = timeDetails_->sys;
        GetRNGstate( ); 
        for (ii = 0; ii < nIters_; ++ii) {
                thisIter_ = ii; (this->*oneIterFP_)(R_NilValue);
                PHONY(PRINT_STUB_INT(thisIter_));
                (this->*registerIterFP_)(SEXPDraws);
                // WRITE ME: CHANGE THE FOLLOWING IMPLEMENTATION
                if (verboseLevel_ >= 1) {
                        iter = ii + 1;
                        if (initialEstDone == false) {
                                if ((iter % printInitialDotsWhen_) == 0) Rprintf(".");
                        }
                        else {
                                if (iter == printDotAt_) {
                                        Rprintf("."); printDotAt_ += eachDotWorth_;
                                }
                        }
                        if (iter == printEstTimeAt_) {
                                gatherTimeDetails( );
                                timeEndUsr = timeDetails_->usr;
                                timeToFinish = ((timeEndUsr - timeStartUsr) / iter * \
                                                (nIters_ - iter));
                                Rprintf("\n[Time to finish (est): %8.0f secs, " \
                                        "this iter: %10d]", ceil(timeToFinish), iter);
                                if (initialEstDone == false) initialEstDone = true;
                                printEstTimeAt_ += eachDotWorth_ * nDotsPerLine_;
                                printDotAt_ = iter + 1;
                        }                        
                }
        }
        PutRNGstate( );
        if (verboseLevel_ >= 1) {
                gatherTimeDetails( );
                timeEndUsr = timeDetails_->usr;
                timeEndSys = timeDetails_->sys;
                Rprintf("\n[Total time: %10.0f secs (usr), %5.0f secs (sys), " \
                        "this iter: %10d]\n", (timeEndUsr - timeStartUsr),
                        (timeEndSys - timeStartSys), iter);
        }
        PHONY(utils_SEXP_darray_print(SEXPDraws, ", "););
        return 0;     
}       

int
Sampler::makeDrawsFixedTime (SEXP SEXPDraws)
{
        RAISE("not implemented yet");
        return 0;
}

SEXP
Sampler::makeDraws (void)
{
        /*
         * draws is an array of dimension:
         * c(nIters_, sampDim_ (+ 1), nLevelsSaveSampFor_)
         */
        SEXP draws, drawsDim;
        int ii, nn, *intsTmp, dims[3];
        
        dims[0] = nIters_; dims[2] = nLevelsSaveSampFor_;
        if (saveFitness_ == TRUE) dims[1] = sampDim_ + 1;
        else                      dims[1] = sampDim_;
        nn = dims[0] * dims[1] * dims[2];
        PROTECT(draws = allocVector(REALSXP, nn));
        CountNProtected::incrementNProtected( );        
        PROTECT(drawsDim = allocVector(INTSXP, 3)); 
        intsTmp = INTEGER(drawsDim);
        for (ii = 0; ii < 3; ++ii) intsTmp[ii] = dims[ii];
        setAttrib(draws, R_DimSymbol, drawsDim);
        UNPROTECT(1);
        if (timeInSecs_ <= 0) makeDrawsFixedIter(draws);
        else                  makeDrawsFixedTime(draws);
        return draws;
}

SEXP
Sampler::makeAcceptRatios (void)
{
        int ii, bb, ll1, ll2, nVarious, row, nProtected = 0;
        double accepted, proposed, *doublesTmp;
        ProposalCounter ***pppc;
        SEXP acceptRatios, aRRowNames, aRColNames, aRDimNames;
        
        for (nVarious = 0, ii = MH; ii <= CE; ++ii) 
                if (moveNTimes_[ii] > 0) ++nVarious;        
        PROTECT(acceptRatios = allocMatrix(REALSXP, nVarious, 4));
        CountNProtected::incrementNProtected( );        
        doublesTmp = REAL(acceptRatios); 
        /* The acceptRatios for MH */
        accepted = 0; proposed = 0; pppc = movePropCtrs_[MH];
        for (ll1 = 0; ll1 < nLevels_; ++ll1) {
                for (bb = 0; bb < MHNBlocks_; ++bb) {
                        accepted += (pppc[ll1][bb])->accepted_;
                        proposed += (pppc[ll1][bb])->proposed_;
                }
        }
        if (fabs(proposed - 0) <= 0) doublesTmp[MH] = R_NaReal;
        else                         doublesTmp[MH] = accepted / proposed;
        doublesTmp[MH + nVarious] = accepted;
        doublesTmp[MH + 2 * nVarious] = proposed;        
        doublesTmp[MH + 3 * nVarious] = moveNTimes_[MH];        
        /* The acceptRatios for RC through CE */
        row = 1;        
        for (ii = RC; ii <= CE; ++ii) {
                /* no row to be filled in if the move wasn't used */
                if (moveNTimes_[ii] == 0) continue;                
                accepted = 0; proposed = 0; pppc = movePropCtrs_[ii];
                for (ll1 = 0; ll1 < nLevels_; ++ll1) {
                        for (ll2 = 0; ll2 < nLevels_; ++ll2) {
                                accepted += (pppc[ll1][ll2])->accepted_;
                                proposed += (pppc[ll1][ll2])->proposed_;
                        }
                }
                if (fabs(proposed - 0) <= 0) doublesTmp[row] = R_NaReal;
                else                         doublesTmp[row] = accepted / proposed;
                doublesTmp[row + nVarious] = accepted;
                doublesTmp[row + 2 * nVarious] = proposed;
                doublesTmp[row + 3 * nVarious] = moveNTimes_[ii];        
                ++row;
        }
        PROTECT(aRRowNames = allocVector(STRSXP, nVarious)); ++nProtected;
        row = 0;
        for (ii = MH; ii <= CE; ++ii) {
                if (moveNTimes_[ii] > 0) {                
                        SET_STRING_ELT(aRRowNames, row, mkChar(moveObjs_[ii]->name_));
                        ++row;
                }
        }
        PROTECT(aRColNames = allocVector(STRSXP, 4)); ++nProtected;
        SET_STRING_ELT(aRColNames, 0, mkChar("ratio"));
        SET_STRING_ELT(aRColNames, 1, mkChar("accepted"));
        SET_STRING_ELT(aRColNames, 2, mkChar("proposed"));
        SET_STRING_ELT(aRColNames, 3, mkChar("moveNTimes"));
        PROTECT(aRDimNames = allocVector(VECSXP, 2)); ++nProtected;
        SET_VECTOR_ELT(aRDimNames, 0, aRRowNames);
        SET_VECTOR_ELT(aRDimNames, 1, aRColNames);
        setAttrib(acceptRatios, R_DimNamesSymbol, aRDimNames);
        UNPROTECT(nProtected);
        return acceptRatios;
}

SEXP
Sampler::makeAcceptRatiosListMH (void)
{
        int bb, ll, jj, nProtected = 0;
        double *doublesTmp, aR;
        ProposalCounter **ppc;
        SEXP aRMat, aRRowNames, aRColNames, aRDimNames;

        PROTECT(aRMat = allocMatrix(REALSXP, nLevels_, MHNBlocks_)); ++nProtected;
        doublesTmp = REAL(aRMat); 
        for (ll = 0; ll < nLevels_; ++ll) {
                ppc = movePropCtrs_[MH][ll];
                for (bb = 0; bb < MHNBlocks_; ++bb) {
                        aR = R_NaReal;
                        if (ppc[bb]->proposed_ > 0)
                                aR = ppc[bb]->accepted_ / ppc[bb]->proposed_;
                        
                        jj             = (nLevels_ * bb + ll);
                        doublesTmp[jj] = aR;
                }
        }

        PROTECT(aRRowNames = allocVector(STRSXP, nLevels_)); ++nProtected;
        for (ll = 0; ll < nLevels_; ++ll) {
                sprintf(charsTmp_, "level%d", ll + 1);
                SET_STRING_ELT(aRRowNames, ll, mkChar(charsTmp_));
        }
        
        PROTECT(aRColNames = allocVector(STRSXP, MHNBlocks_)); ++nProtected;
        for (bb = 0; bb < MHNBlocks_; ++bb) {
                sprintf(charsTmp_, "block%d", bb + 1);
                SET_STRING_ELT(aRColNames, bb, mkChar(charsTmp_));
        }

        PROTECT(aRDimNames = allocVector(VECSXP, 2)); ++nProtected;
        SET_VECTOR_ELT(aRDimNames, 0, aRRowNames);
        SET_VECTOR_ELT(aRDimNames, 1, aRColNames);
        setAttrib(aRMat, R_DimNamesSymbol, aRDimNames);

        UNPROTECT(nProtected);
        return aRMat;        
}

SEXP
Sampler::makeAcceptRatiosListMove (int moveComp)
{
        int ll0, ll1, jj, nProtected = 0;
        double *doublesTmp, proposed, aR;
        ProposalCounter ***pppc = movePropCtrs_[moveComp];
        ProposalCounter *pc01, *pc10;
        SEXP aRMat, aRRowNames, aRColNames, aRDimNames;

        PROTECT(aRMat = allocMatrix(REALSXP, nLevels_, nLevels_)); ++nProtected;
        doublesTmp = REAL(aRMat); 
        for (ll0 = 0; ll0 < nLevels_; ++ll0) {
                for (ll1 = 0; ll1 < nLevels_; ++ll1) {
                        aR       = R_NaReal;
                        pc01     = pppc[ll0][ll1];
                        pc10     = pppc[ll1][ll0]; 
                        proposed = pc01->proposed_ + pc10->proposed_;
                        if (proposed > 0)
                                aR = (pc01->accepted_ + pc10->accepted_) / proposed;

                        jj             = (nLevels_ * ll1 + ll0);
                        doublesTmp[jj] = aR;
                }
        }

        PROTECT(aRRowNames = allocVector(STRSXP, nLevels_)); ++nProtected;
        PROTECT(aRColNames = allocVector(STRSXP, nLevels_)); ++nProtected;
        for (ll0 = 0; ll0 < nLevels_; ++ll0) {
                sprintf(charsTmp_, "level%d", ll0 + 1);
                SET_STRING_ELT(aRRowNames, ll0, mkChar(charsTmp_));
                SET_STRING_ELT(aRColNames, ll0, mkChar(charsTmp_));
        }
        
        PROTECT(aRDimNames = allocVector(VECSXP, 2)); ++nProtected;
        SET_VECTOR_ELT(aRDimNames, 0, aRRowNames);
        SET_VECTOR_ELT(aRDimNames, 1, aRColNames);
        setAttrib(aRMat, R_DimNamesSymbol, aRDimNames);

        UNPROTECT(nProtected);
        return aRMat;        
}

SEXP
Sampler::makeAcceptRatiosListCE (void)
{
        int ii, ll, sl, jj, nLadderLengths, nProtected = 0;
        double *doublesTmp, aR;
        ProposalCounter **ppc;
        SEXP aRMat, aRRowNames, aRColNames, aRDimNames;

        nLadderLengths = moveNTimes_[CE];
        PROTECT(aRMat = allocMatrix(REALSXP, nLadderLengths, nLevels_)); ++nProtected;
        doublesTmp = REAL(aRMat); 
        for (ii = 0; ii < nLadderLengths; ++ii) {
                /* initialize the accept ratio matrix with NAs */
                for (ll = 0; ll < nLevels_; ++ll) {
                        jj             = (nLadderLengths * ll + ii);
                        doublesTmp[jj] = R_NaReal;
                }
                
                ll = ii + 3; ppc = movePropCtrs_[CE][ll - 1];
                for (sl = 0; sl < (ll - 1); ++sl) {
                        aR = R_NaReal;
                        if (ppc[sl]->proposed_ > 0) 
                                aR = ppc[sl]->accepted_ / ppc[sl]->proposed_;

                        jj             = (nLadderLengths * sl + ii);
                        doublesTmp[jj] = aR;
                }
        }
        
        PROTECT(aRRowNames = allocVector(STRSXP, nLadderLengths)); ++nProtected;
        for (ii = 0; ii < nLadderLengths; ++ii) {
                sprintf(charsTmp_, "ladderLength%d", ii + 3);
                SET_STRING_ELT(aRRowNames, ii, mkChar(charsTmp_));
        }
        
        PROTECT(aRColNames = allocVector(STRSXP, nLevels_)); ++nProtected;
        for (ll = 0; ll < nLevels_; ++ll) {
                sprintf(charsTmp_, "level%d", ll + 1);
                SET_STRING_ELT(aRColNames, ll, mkChar(charsTmp_));
        }

        PROTECT(aRDimNames = allocVector(VECSXP, 2)); ++nProtected;
        SET_VECTOR_ELT(aRDimNames, 0, aRRowNames);
        SET_VECTOR_ELT(aRDimNames, 1, aRColNames);
        setAttrib(aRMat, R_DimNamesSymbol, aRDimNames);

        UNPROTECT(nProtected);
        return aRMat;
}

SEXP
Sampler::makeAcceptRatiosList (void)
{
        int ii, nComps, comp = 0, nProtected = 0;
        SEXP aRMat, aRList, aRListNames;

        /* MH will always be used */
        nComps = 1;
        for (ii = RC; ii <= CE; ++ii)
                if (moveNTimes_[ii] > 0) ++nComps;

        PROTECT(aRList      = allocVector(VECSXP, nComps));
        CountNProtected::incrementNProtected( );
        PROTECT(aRListNames = allocVector(STRSXP, nComps)); ++nProtected;

        /* the MH component: should always be included */
        PROTECT(aRMat = makeAcceptRatiosListMH( )); ++nProtected;
        SET_VECTOR_ELT(aRList, comp, aRMat); 
        SET_STRING_ELT(aRListNames, comp, mkChar("MH")); ++comp;

        for (ii = RC; ii <= BSE; ++ii) {
                /* no component to be filled in if the move wasn't used */
                if (moveNTimes_[ii] == 0) continue;                
                
                PROTECT(aRMat = makeAcceptRatiosListMove(ii)); ++nProtected;
                SET_VECTOR_ELT(aRList, comp, aRMat); 
                SET_STRING_ELT(aRListNames, comp, mkChar(moveObjs_[ii]->name_));
                ++comp;
        }

        /* the CE component: if the move was used */
        if (moveNTimes_[CE] > 0) {        
                PROTECT(aRMat = makeAcceptRatiosListCE( )); ++nProtected;
                SET_VECTOR_ELT(aRList, comp, aRMat); 
                SET_STRING_ELT(aRListNames, comp, mkChar("CE")); ++comp;
        }
        
        setAttrib(aRList, R_NamesSymbol, aRListNames);
        UNPROTECT(nProtected);
        return aRList;        
}
        
SEXP
EMCCMainC (SEXP argsList)
{
        Sampler sampler(argsList);
        SEXP samplerObj = R_NilValue;

        PROTECT(samplerObj = sampler.run( ));
        UNPROTECT(1);
        return samplerObj;
}

