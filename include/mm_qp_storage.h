/************************************************************************ *
* Goma - Multiphysics finite element software                             *
* Sandia National Laboratories                                            *
*                                                                         *
* Copyright (c) 2014 Sandia Corporation.                                  *
*                                                                         *
* Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,  *
* the U.S. Government retains certain rights in this software.            *
*                                                                         *
* This software is distributed under the GNU General Public License.      *
\************************************************************************/
 

#ifndef GOMA_MM_QP_STORAGE_H
#define GOMA_MM_QP_STORAGE_H
#include "rf_allo.h"
#include "rf_bc_const.h"
/*
 *$Id: mm_qp_storage.h,v 5.1 2007-09-18 18:53:46 prschun Exp $
 */

/*
 * QP_Storage:
 *
 *    This structure is a linked list, containing a persistent storage
 *  location that a boundary condition or set of boundary conditions
 *  can employ throughout the element assembly process. The storage
 *  is wiped out at the end of a residual or Jacobian fill operation.
 *
 *  The storage is associated with a guass point on an element
 *  side.  The storage is identified primary via an integer number,
 *  StorageType, usually associated with the boundary condition integer
 *  id. The storage consists of a pointer to void, which the boundary
 *  condition can use to store the pointer to its own malloced data.  
 *  Note that it must supply a routine in the function
 *  qp_storage_list_destroy() to destroy that data at the end of a
 *  fill operation.
 */

struct QP_Storage {
/*
 * The following are unused fields, that will be added in later
 * when we try to link qp's together from different elements.
 *
    double Real_x;
    double Real_y;
    double Real_z;
    int    ElemNum;
    int    LocalSideNum;
*/
    int    LocalQPNum;
    int    StorageType;
    void  *Storage;
    struct QP_Storage *next;
};
typedef struct QP_Storage QP_STORAGE_STRUCT;

/*
 * Prototypes for functions in mm_qp_storage.c
 */
EXTERN void elem_qp_storage_free(ELEM_SIDE_BC_STRUCT *);
extern void qp_storage_list_destroy(QP_STORAGE_STRUCT **);
extern void **side_qp_storage_findalloc(const int, const int,
					ELEM_SIDE_BC_STRUCT *);
extern void **qp_storage_findalloc(const int, const int, 
				   QP_STORAGE_STRUCT **);
extern void global_qp_storage_destroy(void);

#endif
