
/*
 *  $Id: cluster.h,v 1.1 2008/02/05 21:44:52 goswami Exp $
 *  
 *  File: cluster.H
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


#ifndef CLUSTER_H
#define CLUSTER_H

#include <iostream>

#ifdef __cplusplus
extern "C" {
#endif

#include <R.h>  
#include <Rinternals.h>        
#include "utils.h"

        static int 
        intcmp(void const *v1, void const *v2)
        {
                return (*(int *)v1 - *(int *)v2);
        }
        
#ifdef __cplusplus
}
#endif

typedef struct ObjLabelContext {
        int obj, oldLabel, newLabel;
} ObjLabelContext;

class Cluster {
public:
	Cluster (int nObjs);
	~Cluster (void);
        // This function provides the R functions a handle on the private data
        inline SEXP getObjLabelsSEXPData (void) const { return objLabelsSEXP_; }
        inline int *getObjLabelsData (void) const { return objLabels_; }
	inline int getNObjs (void) const { return nObjs_; }
	inline int getNClusters (void) const { return nClusters_; }
        int * getClusterMembersData (int clst);
	int getObjLabel (int obj) const;
        int getClusterLabel (int clst) const;
        int getClusterFreq (int clst) const;
        int getClusterFreqFromObj (int obj);
        int getClusterPosFromObj (int obj);
        int setObjLabel (int obj, int label);
        int setObjLabelsAll (int const *labels);
        int setObjLabelsSliceFromOne (int nObjs, int const *objs, int label);
        int setObjLabelAndTabulate (int obj, int label);        
        int generateNewLabel (void);
	int tabulate (void);
        int forceTabulate (void);
        int printCluster (int clst);
	int prettyPrint (void);
        int copy (Cluster const &src);
	bool operator== (Cluster const &cl);
        int registerChangeOfLabel (int obj, int oldLabel, int newLabel);
        int merge (int obj, Cluster &changed, ObjLabelContext *olc);
        int split (int obj, Cluster &changed, ObjLabelContext *olc);
        
private:
        typedef int * IntPtr;
        /*
         * Convention: clusterLabels are restricted to be integers
         * between 0 and n_objects just to make it easier to check
         * that we don't, by mistake, have more clusters than objects.
         * 
         * Note: objLabels_ is just a pointer to the data contained in 
         * objLabelsSEXP_.
         */
        SEXP objLabelsSEXP_;    
        int nObjs_, *objLabels_;
	bool isTabulated_;
	int nClusters_, *clusterLabels_, *clusterFreqs_;
        IntPtr *clusterMembers_;
        /*
         * Scratch spaces.
         */
        int *storeLabels_;
        char errMsg_[MAX_LINE_LENGTH];
};


#endif /* CLUSTER_H */




