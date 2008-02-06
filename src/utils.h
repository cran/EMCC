
/*
 *  $Id: utils.h,v 1.1 2008/02/05 21:44:54 goswami Exp $
 *  
 *  File: utils.H
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

#ifndef UTILS_H
#define UTILS_H

#include <iostream>

#ifdef __cplusplus
extern "C" {
#endif
        
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

#define MAX_WORD_LENGTH 1024
#define MAX_LINE_LENGTH 4096

#define DEBUG(_xx) \
do { \
      Rprintf("===== BEGIN: DEBUG block =====\n" \
               "file:  %s,  line:  %d\n", \
               __FILE__, __LINE__); \
       _xx \
       Rprintf("===== E N D: DEBUG block =====\n"); \
} while (0)

#define PHONY(_xx) ((void) 0)

#define SWAP(_type, _obj1, _obj2) \
do { \
_type _objTmp; \
\
_objTmp = _obj2; \
_obj2 = _obj1; \
_obj1 = _objTmp; \
} while (0)

#define RAISE(_ss) \
do { \
        Rprintf("error: %s [%s:%d]", _ss, __FILE__, __LINE__);   \
        abort( ); \
} while (0)

// #define RAISE(_ss) error((_ss))
#define SQR(_xx) ((_xx) * (_xx))
#define MIN(_aa, _bb) ((_aa) < (_bb) ? (_aa) : (_bb));
#define MAX(_aa, _bb) ((_aa) > (_bb) ? (_aa) : (_bb));
#define PRINT_STUB_INT(_var) Rprintf(#_var " = %d\n", _var);
#define PRINT_STUB_DOUBLE(_var) Rprintf(#_var " = %g\n", _var);
#define PRINT_STUB_STRING(_var) Rprintf(#_var " = %s\n", _var);
#define PRINT_STUB_POINTER(_var) Rprintf(#_var " = %p\n", _var);

#define PRINT_STUB_IARRAY(arr_, nn_, ss_) \
do { \
      Rprintf(#arr_ " [length = %d]:\n", nn_); \
      utils_iarray_print(arr_, nn_, ss_); \
} while (0)

#define PRINT_STUB_DARRAY(arr_, nn_, ss_) \
do { \
      Rprintf(#arr_ " [length = %d]:\n", nn_); \
      utils_darray_print(arr_, nn_, ss_); \
} while (0)

#define PRINT_STUB_SEXP_IARRAY(arr_, ss_) \
do { \
      Rprintf(#arr_ " [length = %d]:\n", length(arr_)); \
      utils_SEXP_iarray_print(arr_, ss_); \
} while (0)

#define PRINT_STUB_SEXP_DARRAY(arr_, ss_) \
do { \
      Rprintf(#arr_ " [length = %d]:\n", length(arr_)); \
      utils_SEXP_darray_print(arr_, ss_); \
} while (0)

        extern int
        utils_iarray_print (int *arr, int nn, char *sep);

        extern int
        utils_darray_print (double *arr, int nn, char *sep);

        extern int
        utils_sarray_print (char **arr, int nn, char *sep);

        extern int
        utils_SEXP_iarray_print (SEXP arr, char *sep);

        extern int
        utils_SEXP_darray_print (SEXP arr, char *sep);

        extern bool
        utils_is_int_in_iarray (int elem, int nn, int *iarr);

        extern int
        utils_unique_iarray_remove_item (int *arr, int nn, int item);
        
        extern bool
        utils_are_disjoint_iarray_iarray (int nn1, int *iarr1, int nn2, int *iarr2);

        extern int
        utils_intersect_iarray_iarray (int nn1, int *iarr1, int nn2, int *iarr2,
                                       int *nIntersection, int *intersection);

        extern int
        utils_sample_bool_mat_cell (int nrows, int ncols, bool **mat, int ntruths,
                                    int *row, int *col);
        
        extern int
        utils_sample_indices_SRSWOR (int nPop, int nSubPop, bool *isAlreadySampled,
                                     int *indices);
        
#ifdef __cplusplus
}
#endif

class CountNProtected {
public:
        CountNProtected (void) : nProtected_(0) { };
        static int incrementNProtected (int by = 1) {
                theOne.nProtected_ += by; return 0;
        }
        static int getNProtected (void) { return theOne.nProtected_; }
        static int resetNProtected (void) { theOne.nProtected_ = 0; return 0; }
                        
private:
        static CountNProtected theOne;
        int nProtected_;
};

#endif /* UTILS_H */

