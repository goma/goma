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
 

/*
 *$Id: mm_qp_storage.c,v 5.1 2007-09-18 18:53:46 prschun Exp $
 */

/* Standard include files */
 
#include <stdio.h>
 
/* GOMA include files */
 
#include "el_elm.h"
#include "rf_allo.h"
#include "rf_bc_const.h"
#include "rf_bc.h"
#include "mm_as_structs.h"
#include "mm_as.h"
#include "mm_qp_storage.h"
#include "mm_interface.h"
#include "rd_mesh.h"
#include "el_elm_info.h"
#include "exo_struct.h"

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void
elem_qp_storage_free(ELEM_SIDE_BC_STRUCT *elem_side_bc)

    /*************************************************************************
     *
     * elem_qp_storage_free():
     *
     *  Frees all of the malloced memory associated with storage of 
     *  bc calculations at quadrature points on the sides of elements. 
     *
     *************************************************************************/
{
  int i, ielem_type, num_qp;
  QP_STORAGE_STRUCT **qps;
  if (elem_side_bc == NULL) return;
  do {
    qps = elem_side_bc->Side_QP_Storage;
    if (qps) {
      ielem_type = Elem_Type(EXO_ptr, elem_side_bc->ielem);
      num_qp = elem_info(NQUAD_SURF, ielem_type);
      for (i = 0; i < num_qp; i++) {
	qp_storage_list_destroy(qps + i);
      }
      safer_free((void **) &elem_side_bc->Side_QP_Storage);
    }
  } while ( (elem_side_bc = elem_side_bc->next_side_bc) != NULL );
  return;
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void
qp_storage_list_destroy(QP_STORAGE_STRUCT **qps_hdl)

    /*************************************************************************
     *
     * qp_storage_list_destroy()
     *
     *    Frees storage associated with a linked list of quadrature point
     * storage structures. It does this by first going to the end of the list
     * and working backwards through recursive calls.
     *
     * This routine needs to know the name of a function that destroys the
     * malloced memory for each type of malloced memory structure. If 
     * it doesn't know about a storage type, then the default is to try to
     * free the address of the storage.
     *************************************************************************/
{
  QP_STORAGE_STRUCT *qps = *qps_hdl;
  if (!qps) return;
  if (qps->next) {
    qp_storage_list_destroy(&(qps->next));
  }
  switch (qps->StorageType) {
  case VL_EQUIL_PRXN_BC:
  case IS_EQUIL_PRXN_BC:
  case SDC_SURFRXN_BC:
      interface_source_destroy((INTERFACE_SOURCE_STRUCT **) &(qps->Storage));
      break;
  default:
      safer_free(&(qps->Storage));
      break;
  }
  safer_free((void **) qps_hdl);
  return;
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void **side_qp_storage_findalloc(const int storageType, const int iquad, 
				 ELEM_SIDE_BC_STRUCT *elem_side_bc)

    /**************************************************************************
     * 
     * side_qp_storage_findalloc()
     *
     *   Given an elem_side_bc structure, this routine will return the 
     *   address of the pointer to the storage (i.e., the handle)
     *   associated with a boundary condition's calculations at a gauss 
     *   point on the surface. The boundary condition is associated via a
     *   storage type integer, storageType. The quadrature point is associated
     *   via its number.
     *   If it can't find the storage, then this function mallocs space for
     *   the storage structure on the surface and returns
     *   the handle to the storage location.
     *
     **************************************************************************/

{
  QP_STORAGE_STRUCT **qpsv;
  int ielem_type, nqp;
  qpsv = elem_side_bc->Side_QP_Storage;
  /*
   * Check to see if we have malloced the vector of pointers for this
   * element side that will base the linked lists on
   */
  if (!qpsv) {
    ielem_type = Elem_Type(EXO_ptr, elem_side_bc->ielem);
    nqp = elem_info(NQUAD_SURF, ielem_type);
    elem_side_bc->Side_QP_Storage = qpsv = 
	(QP_STORAGE_STRUCT **) alloc_ptr_1(nqp);
  }
  return (qp_storage_findalloc(storageType, iquad, qpsv + iquad));
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void **qp_storage_findalloc(const int storageType, const int iquad, 
			    QP_STORAGE_STRUCT **qpsh)

    /**************************************************************************
     * 
     * qp_storage_findalloc()
     *
     *   Given a handle to the start of QP_STORAGE linked list,
     *   this routine will return the 
     *   address of the pointer to the storage (i.e., the handle)
     *   associated with a boundary condition's calculations at a gauss 
     *   point on the surface. The boundary condition is associated via a
     *   storage type integer, storageType. The quadrature point is associated
     *   via its number.
     *   If it can't find the storage, then this function mallocs space for
     *   the storage structure at the end of the linked list and returns
     *   the address of the pointer to the storage (i.e., the handle).
     *
     **************************************************************************/
{
  QP_STORAGE_STRUCT *qps;
  if (! *qpsh) { 
    qps = *qpsh = alloc_struct_1(QP_STORAGE_STRUCT, 1);
  } else {
    for (qps = *qpsh; qps != NULL; qps = qps->next) {
      if (qps->StorageType == storageType) {
	return (&(qps->Storage));
      }
    }
    qps = alloc_struct_1(QP_STORAGE_STRUCT, 1);
  }
  qps->StorageType = storageType;
  qps->LocalQPNum = iquad;
  return  (&(qps->Storage));
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void
global_qp_storage_destroy(void)

   /**************************************************************************
    * 
    * global_qp_storage_destroy()
    *
    *   This routine will deallocate all memory for qp storage in all 
    *   elements.
    *
    **************************************************************************/
{
  int ielem, e_start, e_end;  
  e_start = EXO_ptr->eb_ptr[0];
  e_end   = EXO_ptr->eb_ptr[EXO_ptr->num_elem_blocks];
  for (ielem = e_start; ielem < e_end; ielem++) {
    if (First_Elem_Side_BC_Array[pg->imtrx][ielem]) {
      elem_qp_storage_free(First_Elem_Side_BC_Array[pg->imtrx][ielem]);
    }
  }
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
