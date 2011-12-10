
/*
 *  $Id: utils.cc,v 1.3 2008/02/05 21:44:54 goswami Exp $
 *  
 *  File: utils.C
 *
 *  Created by Gopi Goswami on Wed May 22 2006
 *  Copyright (C) 2006 Gopika R. Goswami
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
 *  goswami@stat.harvard.edu
 *
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

int
utils_iarray_print (int *arr, int nn, const char *sep)
{
        int ii;

        if ((arr == NULL) || (nn < 0))
                Rprintf("MALFORMED array\n");
        else if (nn == 0)
                Rprintf("EMPTY array\n");
        else {
                for (ii = 0; ii < (nn - 1); ++ii) 
                        Rprintf("%d%s", arr[ii], sep);
                Rprintf("%d\n", arr[ii]);
        }
        return 0;
}


int
utils_darray_print (double *arr, int nn, const char *sep)
{
        int ii;

        if ((arr == NULL) || (nn < 0))
                Rprintf("MALFORMED array\n");
        else if (nn == 0)
                Rprintf("EMPTY array\n");
        else {
                for (ii = 0; ii < (nn - 1); ++ii) 
                        Rprintf("%g%s", arr[ii], sep);
                Rprintf("%g\n", arr[ii]);
        }
        return 0;
}


int
utils_sarray_print (char **arr, int nn, const char *sep)
{
        int ii;

        if ((arr == NULL) || (nn < 0))
                Rprintf("MALFORMED array\n");
        else if (nn == 0)
                Rprintf("EMPTY array\n");
        else {
                for (ii = 0; ii < (nn - 1); ++ii) 
                        Rprintf("%s%s", arr[ii], sep);
                Rprintf("%s\n", arr[ii]);
        }
        return 0;
}


int
utils_SEXP_iarray_print (SEXP arr, const char *sep)
{
        int ii, nn = length(arr);

        if ((arr == NULL) || (nn < 0))
                Rprintf("MALFORMED array\n");
        else if (nn == 0)
                Rprintf("EMPTY array\n");
        else {
                for (ii = 0; ii < (nn - 1); ++ii) 
                        Rprintf("%d%s", INTEGER(arr)[ii], sep);
                Rprintf("%d\n", INTEGER(arr)[ii]);
        }
        return 0;
}


int
utils_SEXP_darray_print (SEXP arr, const char *sep)
{
        int ii, nn = length(arr);
        
        if ((arr == NULL) || (nn < 0))
                Rprintf("MALFORMED array\n");
        else if (nn == 0)
                Rprintf("EMPTY array\n");
        else {
                for (ii = 0; ii < (nn - 1); ++ii) 
                        Rprintf("%g%s", REAL(arr)[ii], sep);
                Rprintf("%g\n", REAL(arr)[ii]);
        }
        return 0;
}


bool
utils_is_int_in_iarray (int elem, int nn, int *iarr)
{
        int ii;
        
        if (nn == 0) return false;
        for (ii = 0; ii < nn; ++ii)
                if (elem == iarr[ii]) return true;
        return false;
}


int
utils_unique_iarray_remove_item (int *arr, int nn, int item)
{
        int ii, pos = -1;

        for (ii = 0; ii < nn; ++ii) {
                if (arr[ii] == item) {
                        pos = ii; break;
                }
        }
        if (pos == -1) {
                char errMsg[MAX_LINE_LENGTH];

                sprintf(errMsg,
                        "item %d is not present in the given array of length %d",
                        item, nn);
                RAISE(errMsg);
        }
        if (pos == nn - 1) {
                arr[pos] = -1; return 0;
        }
        for (ii = pos; ii < nn - 1; ++ii)
                arr[ii] = arr[ii + 1];
        arr[nn - 1] = -1;
        return 0;
}


bool
utils_are_disjoint_iarray_iarray (int nn1, int *iarr1, int nn2, int *iarr2)
{
        int ii, jj;

        for (ii = 0; ii < nn1; ++ii) {
                for (jj = 0; jj < nn2; ++jj) {
                        if (iarr1[ii] == iarr2[jj]) return false;
                }
        }
        return true;
}


int
utils_intersect_iarray_iarray (int nn1, int *iarr1, int nn2, int *iarr2,
                               int *nIntersection, int *intersection)
{
        int ii, jj, elem, count = 0;

        for (ii = 0; ii < nn1; ++ii) {                
                for (elem = iarr1[ii], jj = 0; jj < nn2; ++jj) {
                        if (iarr2[jj] == elem) {
                                intersection[count] = elem; ++count;
                                break;
                        }
                }
        }
        *nIntersection = count;
        return 0;                        
}


int
utils_sample_bool_mat_cell (int nrows, int ncols, bool **mat, int ntruths,
                            int *row, int *col)
{
        int ii, jj, ntruths_found = 0;
        double incr = 1.0 / ntruths, sum = 0.0, uu = runif(0, 1);
        char err_msg[MAX_LINE_LENGTH];
        
        for (ii = 0; ii < nrows; ++ii) {
                for (jj = 0; jj < ncols; ++jj) {
                        if (mat[ii][jj] == false) continue;
                        ++ntruths_found; sum += incr;
                        if (uu <= sum) {
                                *row = ii; *col = jj; return 0;
                        }
                }
        }
        sprintf(err_msg,
                "promised ntruths: %d, but ntruths_found: %d, check your C++ code",
                ntruths, ntruths_found);
        RAISE(err_msg);
        return 0;
}

int
utils_sample_indices_SRSWOR (int nPop, int nSubPop, bool *isAlreadySampled,
                             int *indices)
{
	int ii, jj, count, index;
	
	for (ii = 0; ii < nPop; ++ii) isAlreadySampled[ii] = false;
        // the first sample
	index = (int) floor(runif(0, nPop));
	indices[0] = index; isAlreadySampled[index] = true;

        // the rest of the samples
	for (jj = 1; jj < nSubPop; ++jj) {
		index = (int) floor(runif(0, (nPop - jj)));
		count = 0;
		for (ii = 0; ii < nPop; ++ii) {
			if (isAlreadySampled[ii] == true) continue;
                        ++count;
			if (count == (index + 1)) {
                                indices[jj] = ii; isAlreadySampled[ii] = true;
                        }
                }
	}
	return 0;
}

CountNProtected CountNProtected::theOne = CountNProtected( );

