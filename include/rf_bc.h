/************************************************************************ *
* Goma - Multiphysics finite element software                             *
* Sandia National Laboratories                                            *
*                                                                         *
* Copyright (c) 2022 Goma Developers, National Technology & Engineering   *
*               Solutions of Sandia, LLC (NTESS)                          *
*                                                                         *
* Under the terms of Contract DE-NA0003525, the U.S. Government retains   *
* certain rights in this software.                                        *
*                                                                         *
* This software is distributed under the GNU General Public License.      *
* See LICENSE file.                                                       *
\************************************************************************/

/*
 *$Id: rf_bc.h,v 5.2 2007-09-18 18:53:46 prschun Exp $
 */

#ifndef GOMA_RF_BC_H
#define GOMA_RF_BC_H

/*
 *
 *  Author:          Scott Hutchinson (1421)
 *  Date:            20 January 1993
 *  Revised:         22 January 1993
 *
 */

/* BC_Types[]
 *
 *     This is a vector of pointers to type struct Boundary_Condition.
 *  The structure holds the boundary condition information entered from the
 *  input file.  It is dynamically allocated in rf_input.c.
 *     A description of struct Boundary_Condition is contained in rf_bc_const.h.
 *
 *
 */

extern struct Boundary_Condition *BC_Types;
extern struct Rotation_Specs *ROT_Types;
extern int **ROT_list;
extern int **SS_list;
extern int Use_2D_Rotation_Vectors;

/*
 *  First_Elem_Side_BC_Array:
 *
 *     This is an array of pointers to type (struct elem_side_bc_struct).  The
 *  length of the array is equal to the number of elements defined on the
 *  local processor.  A NULL pointer in this array indicates that the current
 *  element does not have an integral boundary condition applied on any of its sides.
 *
 *
 */

extern struct elem_side_bc_struct ***First_Elem_Side_BC_Array;

/*
 *  First_Elem_Edge_BC_Array:
 *
 *     This is an array of pointers to type (struct elem_edge_bc_struct).  The
 *  length of the array is equal to the number of elements defined on the
 *  local processor.  A NULL pointer in this array indicates that the current
 *  element does not have an integral boundary condition applied on any of its edges.
 *
 *
 */

extern struct elem_edge_bc_struct ***First_Elem_Edge_BC_Array;

#endif
