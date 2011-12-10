
/*
 *  $Id: objects.h,v 1.1 2008/02/05 21:44:53 goswami Exp $
 *  
 *  File: objects.H
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


#ifndef OBJECTS_H
#define OBJECTS_H

#include <iostream>
#include "utils.h"
#include "cluster.h"

#ifdef __cplusplus
extern "C" {
#endif

#include <R.h>
#include <Rinternals.h>        
#include <Rmath.h>

        extern SEXP
        EMCCMainC (SEXP argsList);

#ifdef __cplusplus
}
#endif        

#define N_IMPLEMENTED_MOVES 10
#define N_SELECTIONS 2
#define N_SUB_CLUSTERS 2
#define MAX_SAMPLER_MOVE_NAME 100
#define SURV_PARENT_POS 2
#define NON_SURV_PARENT_POS 3

class Sampler {
public:
        Sampler (SEXP argsList);
        ~Sampler(void);
        SEXP run (void);
        
private:
        // The typedefs and enums.
        typedef int    (Sampler::*Util) (SEXP);
        typedef int    (Sampler::*Move) (void);
        typedef double (Sampler::*FuncPtr1) (SEXP);
        typedef int    (Sampler::*ProbAndLevels) (double *, int *);
        typedef int    (Sampler::*ProbGivenLevels) (double *, int *, double *);
        typedef int    (Sampler::*SCRC_sizeComputation) (int *, int *, double *);
        enum MoveTable {
                MH = 0,
                RC = 1,
                SCSC_ONE_NEW = 2,
                SCSC_TWO_NEW = 3,
                SCRC = 4,
                RE = 5,
                BCE = 6,
                BIRE = 7,
                BSE = 8,
                CE = 9
        };
        enum SelectionCode {
                RANDOM = 1,
                BEST,
                WORST,
                SAME_SIZE,
                RANDOM_SIZE
        };
        enum MHMoveCode {
                SPLIT_MERGE = -321,
                GIBBS
        };
        struct TimeDetails {
                double usr, sys;
        };
        struct MoveObject {
                char name_[MAX_SAMPLER_MOVE_NAME];
                Move FP_;
        };
        struct SampleLevelContext {
                double *logWeights_, *adjWeights_, *partialSum_;
                int samp_, endLevel_;
                double prob_, numerator_, sum_, maxLogWeights_;
        };
        struct ProposalCounter {
                double accepted_;
                double proposed_;
        };
        struct ArgsList1 {
                int posDraw_;
                SEXP argsList_;
        };        
        // Some scratch storage
        char charsTmp_[MAX_LINE_LENGTH];
        // The variables
        int nIters_, thisIter_;
        
        double timeInSecs_;
        int nItersActual_;
        
        int verboseLevel_, printEstTimeAt_, printEstTimeNTimes_;
        int printInitialDotsWhen_, printDotAt_, eachDotWorth_, nDotsPerLine_;
        
        int nLevels_, thisLevel_;
        double *temperLadder_, *invTemperLadder_;
        double *logDensities_, *logDensitiesStore_; 
        SampleLevelContext *scratch_SLC_;
        
        int sampDim_, MHNBlocks_, thisBlock_, thisStep_;
        MHMoveCode MHMoveCode_;
        double MHMergeProb_;
        ObjLabelContext *MH_OLC_;
        bool **boolMat_;
        int nSCXX_intersection_[N_SUB_CLUSTERS], *SCXX_intersection_[N_SUB_CLUSTERS];
        int *nSCXX_nonEmptyInRow_, **SCXX_nonEmptyInRow_;
        int nSCXX_rows2OrMore_, *SCXX_rows2OrMore_;
        int *labelDiffs_, *labelDiffsPositive_, *labelDiffsNegative_;
        bool *SCRC_isAlreadySampled_;
        int nSCRC_union_, *SCRC_union_, *SCRC_subSampleIndices_;

        MoveObject *moveObjs_[N_IMPLEMENTED_MOVES];
        double moveProbs_[N_IMPLEMENTED_MOVES], cumsumProbs_[N_IMPLEMENTED_MOVES];
        bool moveProbsPositive_[N_IMPLEMENTED_MOVES];
        int moveNTimes_[N_IMPLEMENTED_MOVES], moveNTimesAtIter_[N_IMPLEMENTED_MOVES];
        SelectionCode moveSelectionCodes_[N_IMPLEMENTED_MOVES][N_SELECTIONS];
        double moveSelectionTempers_[N_IMPLEMENTED_MOVES][N_SELECTIONS];
        int SCSC_ONE_NEW_notPossible_, SCSC_TWO_NEW_notPossible_;
        int SCRC_notPossible_;

        /*
         * MH: nLevels \times sampDim
         * all other moves: nLevels \times nLevels
         */
        ProposalCounter ***movePropCtrs_[N_IMPLEMENTED_MOVES];
        
        int nLevelsSaveSampFor_, *levelsSaveSampFor_;
        Rboolean saveFitness_;

        // The following will be used for logTarDensFunc
        SEXP logTarDensFunc_;
        ArgsList1 *logTarDensArgsList_;
        FuncPtr1 logTarDensFP_;
        
        SEXP doCallFuncCall_, doCallFuncEnv_;
        SEXP procTimeFuncCall_, procTimeFuncEnv_;
        TimeDetails *timeDetails_;        
        
        SEXP dotsList_;
        
        Util oneIterFP_, registerIterFP_;
        ProbAndLevels RC_probAndLevelsFP_, SCSC_TWO_NEW_probAndLevelsFP_;
        ProbAndLevels SCSC_ONE_NEW_probAndLevelsFP_;
        ProbGivenLevels RC_probGivenLevelsFP_, SCSC_TWO_NEW_probGivenLevelsFP_;
        ProbGivenLevels SCSC_ONE_NEW_probGivenLevelsFP_;
        SCRC_sizeComputation SCRC_sizeComputationFP_;
        
        /*
         * nLevels \times sampDim
         */        
        Cluster **clusterCurrDraws_, **clusterCurrDrawsStore_, **clusterPropDraws_;
        SEXP SEXPCurrDraws_, *SEXPCurrDrawsStore_, SEXPPropDraws_;
        double ***drawsArr_;
        
        // The MH move and its helpers.
        int MH_SPLIT_MERGE_propNew (int obj, Cluster &curr, Cluster &prop);
        double MH_SPLIT_MERGE_propDens (int obj, Cluster &curr, Cluster &prop);
        int MH_SPLIT_MERGE_move (void);
        int MH_GIBBS_move (void);
        // The RC move and its helpers.
        int RC_RANDOM_RANDOM_probAndLevels (double *prob, int *sl);
        int RC_RANDOM_RANDOM_probGivenLevels (double *prob, int *sl,
                                              double *logPropDens);
        int RC_BEST_BEST_probAndLevels (double *prob, int *sl);
        int RC_BEST_BEST_probGivenLevels (double *prob, int *sl,
                                          double *logPropDens);
        int RC_move (void);
        // The SCSC_ONE_NEW moves and its helpers.
        int SCSC_ONE_NEW_RANDOM_RANDOM_probAndLevels (double *probs, int *sl);
        int SCSC_ONE_NEW_RANDOM_RANDOM_probGivenLevels (double *probs, int *sl,
                                                        double *logPropDens);
        int SCSC_ONE_NEW_BEST_BEST_probAndLevels (double *probs, int *sl);
        int SCSC_ONE_NEW_BEST_BEST_probGivenLevels (double *probs, int *sl,
                                                    double *logPropDens);
        bool SCSC_ONE_NEW_isParent (Cluster *plausibleParent, Cluster *nonSurvParent,
                                    Cluster *modChild);
        int SCSC_ONE_NEW_propNew (Cluster *survParent, Cluster *nonSurvParent,
                                  Cluster *modChild);
        int SCSC_ONE_NEW_move (void);
        // The SCSC_TWO_NEW moves and its helpers.
        int SCSC_TWO_NEW_RANDOM_RANDOM_probAndLevels (double *prob, int *sl);
        int SCSC_TWO_NEW_RANDOM_RANDOM_probGivenLevels (double *prob, int *sl,
                                                        double *logPropDens);
        int SCSC_TWO_NEW_BEST_BEST_probAndLevels (double *prob, int *sl);
        int SCSC_TWO_NEW_BEST_BEST_probGivenLevels (double *prob, int *sl,
                                                    double *logPropDens);
        int SCSC_TWO_NEW_propNew (Cluster **parents, Cluster **children);
        int SCSC_TWO_NEW_move (void);
        // The SCRC moves and their helpers.
        bool SCRC_isParent (Cluster *plausibleParent, Cluster *nonSurvParent,
                            Cluster *modChild);
        int SCRC_RANDOM_SIZE_sizeComputation (int *sizeIn, int *sizeOut,
                                              double *logDensRatioSize);
        int SCRC_SAME_SIZE_sizeComputation (int *sizeIn, int *sizeOut,
                                            double *logDensRatioSize);
        int SCRC_propNew (Cluster *survParent, Cluster *nonSurvParent,
                          Cluster *modChild, double *logDensRatioSize);        
        int SCRC_move (void);
        // The exchange family of functions.
        int RE_move (void);
        int BCE_move (void);
        int BIRE_move (void);
        int BSE_move (void);
        int CE_move (void);
        // The utility functions.
        SEXP lang5 (SEXP s, SEXP t, SEXP u, SEXP v, SEXP w);
        SEXP getListElement (SEXP list, const char *str);
        int gatherTimeDetails (void);
        ProposalCounter *** PCMatNew (int dim1, int dim2);
        enum SelectionCode getSelectionCode (char const *codeName);
        double logTarDensFuncUserRfunc (SEXP argXX);
        inline double getNeighbourProb (int jjGiven, int ii);
        int exchangeGivenProb (int *sl, 
                               ProposalCounter *pc, 
                               const char *moveName, 
                               double prob);
        int exchange (int *sl, ProposalCounter *pc, const char *moveName);
        int cyclicShiftDraws (int ladderLength, int selLength);
        int initArgsList1 (ArgsList1 *al);
        int init (void);
        int print (void); 
        int sampleWithDetails (void);
        int sampleWithDetailsExcludingOne (int sl);
        int oneIterWithOneLevel (SEXP notRequired);
        int setMoveNTimesAtIter (void);
        int oneIterWithTwoLevels (SEXP notRequired);
        int oneIterWithManyLevels (SEXP notRequired);
        int registerIterFixedIter (SEXP SEXPDraws);
        int registerIterFixedIterWithFitness (SEXP SEXPDraws);
        int registerIterFixedTime (SEXP SEXPDraws);
        int registerIterFixedTimeWithFitness (SEXP SEXPDraws);
        int makeDrawsFixedIter (SEXP SEXPDraws);
        int makeDrawsFixedTime (SEXP SEXPDraws);
        SEXP makeDraws (void);
        SEXP makeAcceptRatios (void);
        SEXP makeAcceptRatiosList (void);
        SEXP makeAcceptRatiosListMH (void);
        SEXP makeAcceptRatiosListMove (int moveComp);
        SEXP makeAcceptRatiosListCE (void);
};              


#endif /* OBJECTS_H */
