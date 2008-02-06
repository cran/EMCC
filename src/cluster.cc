
/*
 *  $Id: cluster.cc,v 1.3 2008/02/05 21:44:52 goswami Exp $
 *  
 *  File: cluster.C
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
#include "utils.h"
#include "cluster.h"

extern "C"  {
#include <R.h>  
#include <Rinternals.h>        
#include <Rmath.h>        
}

#ifndef DEBUG_CLUSTER
#define DEBUG_CLUSTER 0
#endif 

#if (!DEBUG_CLUSTER)
#define DEBUG(_x) ((void) 0)
#endif

Cluster::Cluster (int nObjs)
{        
        if (nObjs <= 0) RAISE("cannot create cluster of non-positive size");
        nObjs_ = nObjs;
        // FIXME: How to 'delete' objLabelsSEXP_?
        PROTECT(objLabelsSEXP_ = allocVector(INTSXP, nObjs_));
        CountNProtected::incrementNProtected( );
        objLabels_ = INTEGER(objLabelsSEXP_);
        clusterLabels_ = new int[nObjs_];
        clusterFreqs_ = new int[nObjs_];
        clusterMembers_ = new IntPtr[nObjs_];
        for (int ii = 0; ii < nObjs_; ++ii) {
                // zeroing out the clusterFreqs_'s
                clusterFreqs_[ii] = 0;
                clusterMembers_[ii] = new int[nObjs_];
        }
        // tabulating
        for (int ii = 0; ii < nObjs_; ++ii)
                objLabels_[ii] = 0;
        nClusters_ = 1;
        clusterLabels_[0] = 0;
        clusterFreqs_[0] = nObjs_;
        for (int ii = 0; ii < clusterFreqs_[0]; ++ii)
                clusterMembers_[0][ii] = ii;
        isTabulated_ = true;        
        // The scratch spaces.
        storeLabels_ = new int[nObjs_];
}


Cluster::~Cluster (void)
{
        DEBUG(Rprintf("Cluster destructor was called\n"););
        delete [] clusterLabels_;
        delete [] clusterFreqs_;
        for (int ii = 0; ii < nObjs_; ++ii) delete [] clusterMembers_[ii];
        delete [] clusterMembers_;
        // The scratch spaces.
        delete [] storeLabels_;
}

int *
Cluster::getClusterMembersData (int clst)
{
        if ((clst < 0) || (clst >= nClusters_)) 
                RAISE("cluster members data requested for invalid cluster number");
        return clusterMembers_[clst];
}


int
Cluster::getObjLabel (int obj) const {
        if ((obj < 0) || (obj >= nObjs_)) 
                RAISE("object label requested for invalid object number");
        return objLabels_[obj];
}


int
Cluster::getClusterLabel (int clst) const {
        if ((clst < 0) || (clst >= nClusters_)) 
                RAISE("cluster label requested for invalid cluster number");
        return clusterLabels_[clst];
}


int
Cluster::getClusterFreq (int clst) const {
        if ((clst < 0) || (clst >= nClusters_)) 
                RAISE("cluster frequency requested for invalid cluster number");
        return clusterFreqs_[clst];
}


int
Cluster::getClusterFreqFromObj (int obj) {
        tabulate( ); int label = getObjLabel(obj);
        for (int ii = 0; ii < nClusters_; ++ii)
                if (clusterLabels_[ii] == label) 
                        return clusterFreqs_[ii];
        sprintf(errMsg_, "cluster frequency for obj: %d could not be found", obj);
        RAISE(errMsg_);
        return 0;
}


int
Cluster::getClusterPosFromObj (int obj) {
        tabulate( ); int label = getObjLabel(obj);
        for (int ii = 0; ii < nClusters_; ++ii)
                if (clusterLabels_[ii] == label) return ii;
        RAISE("should not reach here, check your C++ code");
        return 0;
}

int
Cluster::setObjLabel (int obj, int label) {
        if ((obj < 0) || (obj >= nObjs_)) {
                sprintf(errMsg_, "object label cannot be set for invalid object " \
                        "number: %d", obj);
                RAISE(errMsg_);
        }
        if ((label < 0) || (label >= nObjs_)) {
                sprintf(errMsg_, "invalid object label: %d cannot be set for " \
                        "object: %d", label, obj);
                RAISE(errMsg_);
        }
        objLabels_[obj] = label;
        isTabulated_ = false;
        return 0;
}


int
Cluster::setObjLabelsAll (int const *labels) {
        for (int ii = 0; ii < nObjs_; ++ii) {
                int label = labels[ii];
                if ((label < 0) || (label >= nObjs_)) {
                        sprintf(errMsg_, "invalid object label: %d cannot be set for " \
                                "object: %d", label, ii);
                        RAISE(errMsg_);
                }
                objLabels_[ii] = label;
        }	
        isTabulated_ = false;
        return 0;
}


int
Cluster::setObjLabelsSliceFromOne (int nObjs, int const *objs, int label)
{
        if ((nObjs <= 0) || (nObjs > nObjs_)) {
                sprintf(errMsg_, "invalid nObjs: %d, it should be in [1, %d]",
                        nObjs, nObjs_);
                RAISE(errMsg_);
        }        
        if ((label < 0) || (label >= nObjs_)) {
                sprintf(errMsg_, "invalid object label: %d ", label);
                RAISE(errMsg_);
        }
        for (int ii = 0; ii < nObjs; ++ii) {
                int obj = objs[ii];
                if ((obj < 0) || (obj >= nObjs_)) {
                        sprintf(errMsg_, "object label cannot be set for invalid " \
                                "object number: %d", obj);
                        RAISE(errMsg_);
                }
                objLabels_[obj] = label;
        }
        isTabulated_ = false;
        return 0;
}


int
Cluster::setObjLabelAndTabulate (int obj, int label) {
        setObjLabel(obj, label);
        return tabulate( );
}


/*
 * Note: This returns a new label if one is found, othewise returns -1.
 */
int
Cluster::generateNewLabel (void) {
        tabulate( );
        // no new label could be assigned, since all objects form their 
        // own cluster, and hence return 
        if (nClusters_ == nObjs_) {
                RAISE("cannot generate a new label"); return -1;
        }
        for (int ii = 0; ii < nClusters_; ++ii) 
                        storeLabels_[ii] = clusterLabels_[ii];
        qsort(storeLabels_, nClusters_, sizeof(storeLabels_[0]), intcmp);
        if (storeLabels_[0] > 0) return 0;
        if (storeLabels_[nClusters_ - 1] + 1 < nObjs_) 
                return storeLabels_[nClusters_ - 1] + 1;
        for (int ii = 0; ii < nClusters_ - 1; ++ii)
                if (storeLabels_[ii] + 1 < storeLabels_[ii + 1])
                        return storeLabels_[ii] + 1;
        RAISE("should not reach here, check your C++ code");
        return 0;
}

int
Cluster::tabulate (void)
{
        if (isTabulated_ == true) return 0;
        int label = objLabels_[0];
        nClusters_ = 1; clusterLabels_[0] = label; clusterFreqs_[0] = 1;
        clusterMembers_[0][0] = 0;
        bool openNewCluster;
        for (int ii = 1; ii < nObjs_; ++ii) {
                // zeroing out the clusterFreqs_'s
                clusterFreqs_[ii] = 0;
                label = objLabels_[ii]; openNewCluster = true;
                for (int jj = 0; jj < nClusters_; ++jj) {
                        if (label == clusterLabels_[jj]) {
                                clusterMembers_[jj][clusterFreqs_[jj]] = ii;
                                ++(clusterFreqs_[jj]);
                                openNewCluster = false; break;
                        }
                }
                if (openNewCluster == true) {
                        clusterLabels_[nClusters_] = label;
                        clusterFreqs_[nClusters_] = 1;
                        clusterMembers_[nClusters_][0] = ii;
                        ++nClusters_;
                }
        }
        isTabulated_ = true;        
        return 0;
}


int
Cluster::forceTabulate (void) 
{
        isTabulated_ = false; tabulate( ); return 0;
}


int 
Cluster::printCluster (int clst)
{
        if ((clst < 0) || (clst >= nClusters_)) 
                RAISE("cluster printing requested for invalid cluster number");
        Rprintf("%d: %d: ", clusterLabels_[clst], clusterFreqs_[clst]);
        int jj;
        for (jj = 0; jj < clusterFreqs_[clst] - 1; ++jj)	
                Rprintf("%d, ", clusterMembers_[clst][jj]);	
        Rprintf("%d\n", clusterMembers_[clst][jj]);        
        return 0;
}


int
Cluster::prettyPrint (void)
{
	Rprintf("This cluster: [nObjs: %d, tabulated: %s]\n" \
		"The cluster labels [nClusters: %d]:\n", 
		nObjs_, (isTabulated_ == true) ? "YES" : "NO",
                nClusters_);
	int ii;
	for (ii = 0; ii < nClusters_ - 1; ++ii) 
		Rprintf("%d, ", clusterLabels_[ii]);			
	Rprintf("%d\n", clusterLabels_[ii]);			
	if (isTabulated_ == false) return 0;	
	Rprintf("The table for this cluster [label: frequency: members]\n");
	int jj;
	for (ii = 0; ii < nClusters_; ++ii) {
		Rprintf("%d: %d: ", clusterLabels_[ii], clusterFreqs_[ii]);
		for (jj = 0; jj < clusterFreqs_[ii] - 1; ++jj)	
			Rprintf("%d, ", clusterMembers_[ii][jj]);	
		Rprintf("%d\n", clusterMembers_[ii][jj]);
	}
        Rprintf("The cluster labels:\n");
        for (ii = 0; ii < nObjs_ - 1; ++ii)
		Rprintf("%d, ", objLabels_[ii]);			
	Rprintf("%d\n", objLabels_[ii]);
	return 0;		
}


int 
Cluster::copy (Cluster const &src)
{
        if (nObjs_ != src.nObjs_) 
                RAISE("cannot copy two cluster objects of different sizes");
        int ii;
        for (ii = 0; ii < nObjs_; ++ii)
                objLabels_[ii] = src.objLabels_[ii];
        isTabulated_ = src.isTabulated_;
        if (isTabulated_ == false) return 0;
        // if src is tabulated then copy the rest of the members
        nClusters_ = src.nClusters_;
        for (ii = 0; ii < nClusters_; ++ii) {
                clusterLabels_[ii] = src.clusterLabels_[ii];
                clusterFreqs_[ii] = src.clusterFreqs_[ii];
                for (int jj = 0; jj < src.clusterFreqs_[ii]; ++jj)
                        clusterMembers_[ii][jj] = src.clusterMembers_[ii][jj];
        }
        return 0;
}


bool
Cluster::operator== (Cluster const &cl)
{
        if (nObjs_ != cl.nObjs_) return false;
        if (nClusters_ != cl.nClusters_) return false;
        for (int ii = 0; ii < nClusters_; ++ii) {
                if (clusterFreqs_[ii] != cl.clusterFreqs_[ii]) return false;
                for (int jj = 0; jj < clusterFreqs_[ii]; ++jj) {
                        if (clusterMembers_[ii][jj] != cl.clusterMembers_[ii][jj]) 
                                return false;
                }
        }
        return true;        
}


int
Cluster::registerChangeOfLabel (int obj, int oldLabel, int newLabel)
{
        if (isTabulated_ == false)
                RAISE("cannot register change of label, do tabulation first");
        if (getObjLabel(obj) != oldLabel)
                RAISE("oldLabel does not match existing oldLabel");
        
        int ii, oldCluster = -1, newCluster = -1;
        for (ii = 0; ii < nClusters_; ++ii) {
                if (clusterLabels_[ii] == oldLabel) {
                        oldCluster = ii; break;
                }
        }
        if (oldCluster == -1)
                RAISE("oldCluster was not found, corrupted previous tabulation");
        for (ii = 0; ii < nClusters_; ++ii) {
                if (clusterLabels_[ii] == newLabel) {
                        newCluster = ii; break;
                }
        }

        objLabels_[obj] = newLabel;
        if (newCluster == -1) {
                // add obj to a new cluster
                clusterLabels_[nClusters_] = newLabel;
                clusterFreqs_[nClusters_] = 1;
                clusterMembers_[nClusters_][0] = obj;
                ++nClusters_;
        } else {
                // merge obj with an existing cluster
                ++(clusterFreqs_[newCluster]);
                int ff = clusterFreqs_[newCluster];
                clusterMembers_[newCluster][ff - 1] = obj;
        }

        // remove obj from the oldCluster
        utils_unique_iarray_remove_item(clusterMembers_[oldCluster],
                                        clusterFreqs_[oldCluster], obj);
        --(clusterFreqs_[oldCluster]);
        // delete the oldCluster if it's empty
        if (clusterFreqs_[oldCluster] == 0) {
                int jj = nClusters_ - 1;

                // copy the contents of the present last cluster to
                // this empty cluster
                clusterLabels_[oldCluster] = clusterLabels_[jj];
                clusterFreqs_[oldCluster] = clusterFreqs_[jj];
                for (ii = 0; ii < clusterFreqs_[jj]; ++ii) 
                        clusterMembers_[oldCluster][ii] = clusterMembers_[jj][ii];
                --nClusters_;
        }
        return 0;
}

int 
Cluster::merge (int obj, Cluster &changed, ObjLabelContext *olc)
{
        tabulate( ); 
        if (nClusters_ == 1) RAISE("merging is useless here");        
        int label = getObjLabel(obj);
        double uu = runif(0, 1), sum = 0.0, incr = 1.0 / (nClusters_ - 1);
        for (int ii = 0; ii < nClusters_; ++ii) {
                int clusterLabel = getClusterLabel(ii);
                if (clusterLabel == label) continue;
                sum += incr;
                if (uu <= sum) {
                        olc->obj = obj;
                        olc->oldLabel = label;
                        olc->newLabel = clusterLabel;
                        DEBUG(Rprintf("merging obj: %d: %d ==> %d\n", obj,
                                      label, clusterLabel););
                        changed.registerChangeOfLabel(olc->obj, olc->oldLabel, olc->newLabel);
                        return 0;
                }
        }
        RAISE("should not reach here");
        return 0;
}


int 
Cluster::split (int obj, Cluster &changed, ObjLabelContext *olc)
{
        int newLabel = generateNewLabel( );
        olc->obj = obj;
        olc->oldLabel = getObjLabel(obj);
        olc->newLabel = newLabel;
        DEBUG(Rprintf("splitting obj: %d: %d ==> %d\n", obj,
                      getObjLabel(obj), newLabel););      
        changed.registerChangeOfLabel(olc->obj, olc->oldLabel, olc->newLabel);
        return 0;
}
