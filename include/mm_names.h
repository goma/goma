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
 *$Id: mm_names.h,v 5.23 2010-07-21 16:39:27 hkmoffa Exp $
 */
/*
 * Here's a RECIPE for adding new boundary conditions so you don't have any
 * excuses not to add new ones.  The changes should be made in at least
 * four files (rf_bc_const.h, mm_names.h, mm_input.c, and bc_[method].c)
 * for some boundary conditions you may want to make additions elsewhere also.
 * One example of extra additions is in el_exoII_io.c where the ss_dup_list
 * is created - you may want to adapt the logic for rotation of new conditions.
 * (note that these lines are repeated at each place where you need to 
 * make changes):
 *  Boundary Condition  Step 1: Add Macro Substitution for your new boundary
 *                              condition in rf_bc_const.h - this is how you
 *      rf_bc_const.h           will refer to the new boundary condition 
 *                              throughout the code.  Make sure the integer
 *                              you choose is unique from all the other BC
 *                              types.
 *  Boundary Condition  Step 2: Add a description of your boundary condition
 *                              to the BC_Desc structure array in mm_names.h.
 *      mm_names.h              This structure includes information about the 
 *                              type of boundary condition, which equation it
 *                              applies to, what variables it is sensitive to,
 *                              whether to rotate mesh or momentum equations, 
 *                              etc.  It is very important that you fill out
 *                              this structure carefully, otherwise the code
 *                              won't know what to do.
 *  Boundary Condition  Step 3: Add your BC case to the correct format listing 
 *                              for reading the necessary arguments from the 
 *				goma input file in mm_input_bc.c
 *
 *  Boundary Condition  Step 3a: Some BC's may require changes in mm_bc.c
 *                              similar to what is done for PLANE, GEOM, CA, 
 *                              and all dirichlet conditions 
 *
 *  Boundary Condition  Step 4: Add a function call (and a function) in the 
 *                              correct routine for evaluating your boundary
 *      bc_colloc.c             condition.  This will probably be located in bc_colloc.c
 *      bc_integ.c              for collocated conditions or bc_integ.c for 
 *                              strong or weak integrated conditions.
 *  Boundary Condition  Step 5: Use and enjoy your new boundary condition
 *
 * Step 2 should be done below:
 *
 ***********************************************************************************
 * 
 * BC_descriptions:
 *
 * Structure to control input of boundary conditions:
 * contains three entries: 
 *   name1:    string for name of this bc in input deck
 *   name2:    alternate string for name of this bc in input deck
 *   method:   Descriptor of the method by which this condition is applied
 *   BC_Name:  integer which corresponds to this bc name (listed in rf_bc_const.h)
 *   equation: equation to which this condition is applied
 *   vector:   flag to indicate if this is a vector condition 
 *   rotate:   flag to indicate if the corresponding equations should be rotated
 *   sens:     flags to indicate which variables this bc is sensitive to
 *   i_apply:  flag which indicates if this condition only applies on one side of the  
 *             boundary, or if it should be evaluated from both sides
 *   DV_Index_Default: Default methodology for application of the variable
 *             in the case of discontinuous variables at an interface.
 *
 * NOTE: it is important to put in all of the variables for each BC 
 *       (or set them to default values) otherwise, the initialization of this
 *        array will be screwed up
 *
 * Note: For Dirichlet conditions, the sensitivity vector is not used, so you can
 *       just fill it with zeros if you want!
 *
 * This structure is initialized here so it doesn't need to be initialized anywhere
 *  else
 *
 * {
 *   char *name1, *name2;
 *   int method;
 *   int  BC_Name; 
 *   int equation;
 *   int vector;
 *   int rotate;
 *   int sens[MAX_VARIABLE_TYPES];
 *   int i_apply;
 *   int DV_Index_Default;
 * } 
 */
#ifndef _MM_NAMES_H
#define _MM_NAMES_H

#include "rf_bc_const.h"	/* just in case it has not yet been included */

/*   Here is a ruler for setting up sensitivity settings for the boundary 
condition.  Just copy the line of zeros below the first line of zeros and 
change the zeros to 1 under the appropriate var value. */
/*                                                                                                                                                                                                        1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1*/
/*                    1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6,6,6,6,6,7,7,7,7,7,7,7,7,7,7,8,8,8,8,8,8,8,8,8,8,9,9,9,9,9,9,9,9,9,9,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6,6,6,6,6,7,7,7,7,7,7,7,7,7,7,8,8,8,8,8,8,8,8,8,8,9,9,9,9,9,9,9,9,9,9*/
/*0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9*/

/*0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,*/
/**/
struct BC_descriptions  BC_Desc[] = 
{

  { "DX", "DX_BC", DIRICHLET, DX_BC, R_MESH1, SCALAR, NO_ROT, {0, 0, 0, 0, 0, 1, 0, 0, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB },
  { "DY", "DY_BC", DIRICHLET, DY_BC, R_MESH2, SCALAR, NO_ROT, {0, 0, 0, 0, 0, 0, 1, 0, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB },
  { "DZ", "DZ_BC", DIRICHLET, DZ_BC, R_MESH3, SCALAR, NO_ROT, {0, 0, 0, 0, 0, 0, 0, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB },
  { "DX_RS", "DX_RS_BC", DIRICHLET, DX_RS_BC, R_SOLID1, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB },
  { "DY_RS", "DY_RS_BC", DIRICHLET, DY_RS_BC, R_SOLID2, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "DZ_RS", "DZ_RS_BC", DIRICHLET, DZ_RS_BC, R_SOLID3, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB },
  { "DX_NOTHING", "DX_NOTHING_BC", DIRICHLET, DX_NOTHING_BC, R_MESH1, SCALAR, NO_ROT, {0, 0, 0, 0, 0, 1, 0, 0, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB },
  { "DY_NOTHING", "DY_NOTHING_BC", DIRICHLET, DY_NOTHING_BC, R_MESH2, SCALAR, NO_ROT, {0, 0, 0, 0, 0, 0, 1, 0, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB },
  { "DZ_NOTHING", "DZ_NOTHING_BC", DIRICHLET, DZ_NOTHING_BC, R_MESH3, SCALAR, NO_ROT, {0, 0, 0, 0, 0, 0, 0, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB },

  {  "DISTNG",   "DISTNG_BC",   COLLOCATE_SURF, DISTNG_BC, R_MESH_NORMAL, SCALAR, R_MESH1, {0, 0, 0, 1, 0, 0, 0, 0, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "DXDISTNG", "DXDISTNG_BC", COLLOCATE_SURF, DXDISTNG_BC, R_MESH1, SCALAR, NO_ROT,  {0, 0, 0, 1, 0, 0, 0, 0, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "DYDISTNG", "DYDISTNG_BC", COLLOCATE_SURF, DYDISTNG_BC, R_MESH2, SCALAR, NO_ROT,  {0, 0, 0, 1, 0, 0, 0, 0, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  {  "DZDISTNG", "DZDISTNG_BC", COLLOCATE_SURF, DZDISTNG_BC, R_MESH3, SCALAR, NO_ROT,  {0, 0, 0, 1, 0, 0, 0, 0, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,

  { "SPLINE",  "SPLINE_BC",  COLLOCATE_SURF, SPLINE_BC, R_MESH_NORMAL,SCALAR, R_MESH1, {0, 0, 0, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "SPLINEX", "SPLINEX_BC", COLLOCATE_SURF, SPLINEX_BC, R_MESH1, SCALAR, NO_ROT,  {0, 0, 0, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "SPLINEY", "SPLINEY_BC", COLLOCATE_SURF, SPLINEY_BC, R_MESH2, SCALAR, NO_ROT,  {0, 0, 0, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "SPLINEZ", "SPLINEZ_BC", COLLOCATE_SURF, SPLINEZ_BC, R_MESH3, SCALAR, NO_ROT,  {0, 0, 0, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "SPLINE_RS",  "SPLINE_RS_BC",  COLLOCATE_SURF, SPLINE_RS_BC, R_SOLID_NORMAL,SCALAR, R_SOLID1, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "SPLINEX_RS",  "SPLINEX_RS_BC",  COLLOCATE_SURF, SPLINEX_RS_BC, R_SOLID1,SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "SPLINEY_RS",  "SPLINEY_RS_BC",  COLLOCATE_SURF, SPLINEY_RS_BC, R_SOLID2,SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "SPLINEZ_RS",  "SPLINEZ_RS_BC",  COLLOCATE_SURF, SPLINEZ_RS_BC, R_SOLID3,SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "GEOM",  "GEOM_BC",  COLLOCATE_SURF, GEOM_BC, R_MESH_NORMAL,SCALAR, R_MESH1, {0, 0, 0, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "GEOMX", "GEOMX_BC", COLLOCATE_SURF, GEOMX_BC, R_MESH1, SCALAR, NO_ROT,  {0, 0, 0, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "GEOMY", "GEOMY_BC", COLLOCATE_SURF, GEOMY_BC, R_MESH2, SCALAR, NO_ROT,  {0, 0, 0, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "GEOMZ", "GEOMZ_BC", COLLOCATE_SURF, GEOMZ_BC, R_MESH3, SCALAR, NO_ROT,  {0, 0, 0, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,

  { "PLANE",  "PLANE_BC",    COLLOCATE_SURF, PLANE_BC, R_MESH_NORMAL, SCALAR, R_MESH1, {0, 0, 0, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "PLANEX", "PLANEX_BC",   COLLOCATE_SURF, PLANEX_BC,  R_MESH1, SCALAR, NO_ROT,  {0, 0, 0, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "PLANEY", "PLANEY_BC",   COLLOCATE_SURF, PLANEY_BC,  R_MESH2, SCALAR, NO_ROT,  {0, 0, 0, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "PLANEZ", "PLANEZ_BC",   COLLOCATE_SURF, PLANEZ_BC,  R_MESH3, SCALAR, NO_ROT,  {0, 0, 0, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "MOVING_PLANE",  "MOVING_PLANE_BC",    COLLOCATE_SURF, MOVING_PLANE_BC, R_MESH_NORMAL, SCALAR, R_MESH1, {0, 0, 0, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "MOVING_PLANE_ETCH",  "MOVING_PLANE_ETCH_BC",    COLLOCATE_SURF, MOVING_PLANE_ETCH_BC, R_MESH_NORMAL, SCALAR, R_MESH1, {0, 0, 0, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "MESH_CONSTRAINT",  "MESH_CONSTRAINT_BC",  COLLOCATE_SURF, MESH_CONSTRAINT_BC,  R_MESH_NORMAL, SCALAR, R_MESH1, {0, 0, 0, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB },
  { "FILLET",  "FILLET_BC",    COLLOCATE_SURF, FILLET_BC, R_MESH_NORMAL, SCALAR, R_MESH1, {0, 0, 0, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "DOUBLE_RAD",  "DOUBLE_RAD_BC",    COLLOCATE_SURF, DOUBLE_RAD_BC, R_MESH_NORMAL, SCALAR, R_MESH1, {0, 0, 0, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "ROLL_FLUID",  "ROLL_FLUID_BC",    COLLOCATE_SURF, ROLL_FLUID_BC, R_MESH_NORMAL, SCALAR, R_MESH1, {1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "P", "P_BC", DIRICHLET, P_BC, R_PRESSURE, SCALAR, NO_ROT, {0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,

  { "PSPG", "PSPG_BC", WEAK_INT_SURF, PSPG_BC, R_PRESSURE, SCALAR, NO_ROT, {1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0} , SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "FILL_INLET", "FILL_INLET_BC", SPECIAL, FILL_INLET_BC, R_FILL, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1} , SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
#ifdef COUPLED_FILL
  { "FILL_CA", "FILL_CA_BC", WEAK_INT_SURF, FILL_CA_BC, R_MOMENTUM1, VECTOR, NO_ROT, {0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1} , SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
#else /* COUPLED_FILL */
  { "FILL_CA", "FILL_CA_BC", WEAK_INT_SURF, FILL_CA_BC, R_MOMENTUM1, VECTOR, NO_ROT, {0, 0, 0, 0, 0, 1, 1, 1, 0, 0} , SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
#endif /* COUPLED_FILL */
  { "SHARP_CA_2D", "SHARP_CA_2D_BC", WEAK_SHARP_INT, SHARP_CA_2D_BC, R_MOMENTUM1, VECTOR, NO_ROT, {0, 0, 0, 0, 0, 1, 1, 1, 0, 0} , SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "SHARP_WETLIN_VELOCITY", "SHARP_WETLIN_VELOCITY_BC", WEAK_SHARP_INT, SHARP_WETLIN_VELOCITY_BC, R_MOMENTUM1, VECTOR, NO_ROT, {1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1} , SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "SHARP_BLAKE_VELOCITY", "SHARP_BLAKE_VELOCITY_BC", WEAK_SHARP_INT, SHARP_BLAKE_VELOCITY_BC, R_MOMENTUM1, VECTOR, NO_ROT, {1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1} , SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "SHARP_HOFFMAN_VELOCITY", "SHARP_HOFFMAN_VELOCITY_BC", WEAK_SHARP_INT, SHARP_HOFFMAN_VELOCITY_BC, R_MOMENTUM1, VECTOR, NO_ROT, {1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1} , SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "SHARP_COX_VELOCITY", "SHARP_COX_VELOCITY_BC", WEAK_SHARP_INT, SHARP_COX_VELOCITY_BC, R_MOMENTUM1, VECTOR, NO_ROT, {1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1} , SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "SHARP_SHIK_VELOCITY", "SHARP_SHIK_VELOCITY_BC", WEAK_SHARP_INT, SHARP_SHIK_VELOCITY_BC, R_MOMENTUM1, VECTOR, NO_ROT, {1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1} , SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "WETTING_TENSION", "WETTING_TENSION_BC", WEAK_INT_SURF, WETTING_TENSION_BC, R_MOMENTUM1, VECTOR, NO_ROT, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0} , SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "STRONG_FILL_CA", "STRONG_FILL_CA_BC", STRONG_INT_SURF, STRONG_FILL_CA_BC, R_FILL, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1} , SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "LS_INLET", "LS_INLET_BC", COLLOCATE_SURF, LS_INLET_BC, R_FILL, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1} , SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "QCONV", "QCONV_BC", WEAK_INT_SURF, QCONV_BC, R_ENERGY, SCALAR, NO_ROT, {0, 0, 0, 1, 0, 0, 0, 0, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "QUSER", "QUSER_BC", WEAK_INT_SURF, QUSER_BC, R_ENERGY, SCALAR, NO_ROT, {1, 1, 1, 1, 1, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "Q_LASER_WELD", "Q_LASER_WELD_BC", WEAK_INT_SURF, Q_LASER_WELD_BC, R_ENERGY, SCALAR, NO_ROT, {1, 1, 1, 1, 1, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "Q_RAIL", "Q_RAIL_BC", WEAK_INT_SURF, Q_RAIL_BC, R_ENERGY, SCALAR, NO_ROT, {1, 1, 1, 1, 1, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "Q_VAPOR", "Q_VAPOR_BC", WEAK_INT_SURF, Q_VAPOR_BC, R_ENERGY, SCALAR, NO_ROT, {1, 1, 1, 1, 1, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "QRAD",  "QRAD_BC",  WEAK_INT_SURF, QRAD_BC,  R_ENERGY, SCALAR, NO_ROT, {0, 0, 0, 1, 0, 0, 0, 0, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "QRAD_REPULSE_ROLL",  "QRAD_REPULSE_ROLL_BC",  WEAK_INT_SURF, QRAD_REPULSE_ROLL_BC,  R_ENERGY, SCALAR, NO_ROT, {0, 0, 0, 1, 1, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "QNOBC",  "QNOBC_BC",  WEAK_INT_SURF, QNOBC_BC,  R_ENERGY, SCALAR, NO_ROT, {0,0,0,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "POTENTIAL_NOBC",  "POTENTIAL_NOBC_BC",  WEAK_INT_SURF, POTENTIAL_NOBC_BC,  R_POTENTIAL, SCALAR, NO_ROT, {0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB},
  { "QSIDE", "QSIDE_BC", WEAK_INT_SURF, QSIDE_BC, R_ENERGY, SCALAR, NO_ROT, {0, 0, 0, 1, 0, 0, 0, 0, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "QSIDE_LS", "QSIDE_LS_BC", WEAK_INT_SURF, QSIDE_LS_BC, R_ENERGY, SCALAR, NO_ROT, {0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, CROSS_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "QSIDE_DIRECTION", "QSIDE_DIR_BC", WEAK_INT_SURF, QSIDE_DIR_BC, R_ENERGY, SCALAR, NO_ROT, {0, 0, 0, 1, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "Q_VELO_SLIP", "Q_VELO_SLIP_BC", WEAK_INT_SURF, Q_VELO_SLIP_BC, R_ENERGY, SCALAR, NO_ROT, {1, 1, 1, 1, 0, 0, 0, 0, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "T_MELT", "T_MELT_BC", STRONG_INT_SURF, T_MELT_BC, R_MESH_NORMAL, SCALAR, R_MESH1, {0, 0, 0, 1, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  
  { "N1", "N1_BC", DIRICHLET, N1_BC, R_NORMAL1, SCALAR, NO_ROT, {0, 0, 0, 0, 0, 1, 0, 0, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB },
  { "N2", "N2_BC", DIRICHLET, N2_BC, R_NORMAL2, SCALAR, NO_ROT, {0, 0, 0, 0, 0, 0, 1, 0, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB },
  { "N3", "N3_BC", DIRICHLET, N3_BC, R_NORMAL3, SCALAR, NO_ROT, {0, 0, 0, 0, 0, 0, 0, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB },

  { "S11", "S11_BC", DIRICHLET, S11_BC, R_STRESS11, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "S12", "S12_BC", DIRICHLET, S12_BC, R_STRESS12, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "S13", "S13_BC", DIRICHLET, S13_BC, R_STRESS13, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "S22", "S22_BC", DIRICHLET, S22_BC, R_STRESS22, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "S23", "S23_BC", DIRICHLET, S23_BC, R_STRESS23, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "S33", "S33_BC", DIRICHLET, S33_BC, R_STRESS33, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,

  { "S11_1", "S11_1_BC", DIRICHLET, S11_1_BC, R_STRESS11_1, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "S12_1", "S12_1_BC", DIRICHLET, S12_1_BC, R_STRESS12_1, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "S13_1", "S13_1_BC", DIRICHLET, S13_1_BC, R_STRESS13_1, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "S22_1", "S22_1_BC", DIRICHLET, S22_1_BC, R_STRESS22_1, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "S23_1", "S23_1_BC", DIRICHLET, S23_1_BC, R_STRESS23_1, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,


  { "S11_2", "S11_2_BC", DIRICHLET, S11_2_BC, R_STRESS11_2, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "S12_2", "S12_2_BC", DIRICHLET, S12_2_BC, R_STRESS12_2, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "S13_2", "S13_2_BC", DIRICHLET, S13_2_BC, R_STRESS13_2, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "S22_2", "S22_2_BC", DIRICHLET, S22_2_BC, R_STRESS22_2, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "S23_2", "S23_2_BC", DIRICHLET, S23_2_BC, R_STRESS23_2, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,

  { "S11_3", "S11_3_BC", DIRICHLET, S11_3_BC, R_STRESS11_3, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "S12_3", "S12_3_BC", DIRICHLET, S12_3_BC, R_STRESS12_3, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "S13_3", "S13_3_BC", DIRICHLET, S13_3_BC, R_STRESS13_3, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "S22_3", "S22_3_BC", DIRICHLET, S22_3_BC, R_STRESS22_3, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "S23_3", "S23_3_BC", DIRICHLET, S23_3_BC, R_STRESS23_3, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,

  { "S11_4", "S11_4_BC", DIRICHLET, S11_4_BC, R_STRESS11_4, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "S12_4", "S12_4_BC", DIRICHLET, S12_4_BC, R_STRESS12_4, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "S13_4", "S13_4_BC", DIRICHLET, S13_4_BC, R_STRESS13_4, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "S22_4", "S22_4_BC", DIRICHLET, S22_4_BC, R_STRESS22_4, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "S23_4", "S23_4_BC", DIRICHLET, S23_4_BC, R_STRESS23_4, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,

  { "S11_5", "S11_5_BC", DIRICHLET, S11_5_BC, R_STRESS11_5, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "S12_5", "S12_5_BC", DIRICHLET, S12_5_BC, R_STRESS12_5, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "S13_5", "S13_5_BC", DIRICHLET, S13_5_BC, R_STRESS13_5, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "S22_5", "S22_5_BC", DIRICHLET, S22_5_BC, R_STRESS22_5, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "S23_5", "S23_5_BC", DIRICHLET, S23_5_BC, R_STRESS23_5, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,

  { "S11_6", "S11_6_BC", DIRICHLET, S11_6_BC, R_STRESS11_6, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "S12_6", "S12_6_BC", DIRICHLET, S12_6_BC, R_STRESS12_6, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "S13_6", "S13_6_BC", DIRICHLET, S13_6_BC, R_STRESS13_6, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "S22_6", "S22_6_BC", DIRICHLET, S22_6_BC, R_STRESS22_6, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "S23_6", "S23_6_BC", DIRICHLET, S23_6_BC, R_STRESS23_6, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,

  { "S11_7", "S11_7_BC", DIRICHLET, S11_7_BC, R_STRESS11_7, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "S12_7", "S12_7_BC", DIRICHLET, S12_7_BC, R_STRESS12_7, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "S13_7", "S13_7_BC", DIRICHLET, S13_7_BC, R_STRESS13_7, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "S22_7", "S22_7_BC", DIRICHLET, S22_7_BC, R_STRESS22_7, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "S23_7", "S23_7_BC", DIRICHLET, S23_7_BC, R_STRESS23_7, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,

  { "STRESS_DEVELOPED", "STRESS_DEVELOPED_BC", STRONG_INT_SURF, STRESS_DEVELOPED_BC, R_STRESS11, STRESS, NO_ROT, {0,0,0,0,0,0,0,0,0},
  SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,

  { "G11", "G11_BC", DIRICHLET, G11_BC, R_GRADIENT11, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "G12", "G12_BC", DIRICHLET, G12_BC, R_GRADIENT12, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0},
  SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "G21", "G21_BC", DIRICHLET, G21_BC, R_GRADIENT21, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0},
  SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "G22", "G22_BC", DIRICHLET, G22_BC, R_GRADIENT22, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0},
  SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "G13", "G13_BC", DIRICHLET, G13_BC, R_GRADIENT13, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0},
  SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  {  "G23", "G23_BC", DIRICHLET, G23_BC, R_GRADIENT23, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0},
  SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "G31", "G31_BC", DIRICHLET, G31_BC, R_GRADIENT31, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0},
  SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "G32", "G32_BC", DIRICHLET, G32_BC, R_GRADIENT32, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0},
  SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "G33", "G33_BC", DIRICHLET, G33_BC, R_GRADIENT33, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1},
  SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,

  { "VOLT", "VOLT_BC", DIRICHLET, VOLT_BC, R_POTENTIAL, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0},
  SINGLE_PHASE },
  { "VOLT_USER", "VOLT_USER_BC", STRONG_INT_SURF, VOLT_USER_BC, R_POTENTIAL, SCALAR, NO_ROT, {0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0},
  SINGLE_PHASE },
  { "CURRENT", "CURRENT_BC", WEAK_INT_SURF, CURRENT_BC, R_POTENTIAL, SCALAR, NO_ROT, {0, 0, 0, 1, 0, 0, 0, 0, 0, 0},SINGLE_PHASE },
  { "CURRENT_SIC", "CURRENT_SIC_BC", STRONG_INT_SURF, CURRENT_SIC_BC, R_POTENTIAL, SCALAR, NO_ROT, {0, 0, 0, 1, 0, 0, 0, 0, 0, 0},SINGLE_PHASE },
  { "CURRENT_BV", "CURRENT_BV_BC", WEAK_INT_SURF, CURRENT_BV_BC, R_POTENTIAL, SCALAR, NO_ROT, {0, 0, 0, 1, 0, 0, 0, 0, 0, 0},SINGLE_PHASE },
  { "CURRENT_ORR", "CURRENT_ORR_BC", WEAK_INT_SURF, CURRENT_ORR_BC, R_POTENTIAL, SCALAR, NO_ROT, {0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0}, SINGLE_PHASE },
  { "CURRENT_HOR", "CURRENT_HOR_BC", WEAK_INT_SURF, CURRENT_HOR_BC, R_POTENTIAL, SCALAR, NO_ROT, {0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0}, SINGLE_PHASE },
  { "CURRENT_BV2", "CURRENT_BV2_BC", WEAK_INT_SURF, CURRENT_BV2_BC, R_POTENTIAL, SCALAR, NO_ROT, {0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0}, SINGLE_PHASE },  /* RSL 1/15/01 */
  { "CURRENT_NI", "CURRENT_NI_BC", WEAK_INT_SURF, CURRENT_NI_BC, R_POTENTIAL, SCALAR, NO_ROT, {0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0}, SINGLE_PHASE },  /* RSL 3/9/01 */
  { "CURRENT_USER", "CURRENT_USER_BC", WEAK_INT_SURF, CURRENT_USER_BC, R_POTENTIAL, SCALAR, NO_ROT, {0, 0, 0, 1, 0, 0, 0, 0, 0, 0},SINGLE_PHASE },
  { "CURRENT_USER_SIC", "CURRENT_USER_SIC_BC", STRONG_INT_SURF, CURRENT_USER_SIC_BC, R_POTENTIAL, SCALAR, NO_ROT, {0, 0, 0, 1, 0, 0, 0, 0, 0, 0},SINGLE_PHASE },
  { "QS", "QS_BC", DIRICHLET, QS_BC, R_SURF_CHARGE, SCALAR, NO_ROT, {0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB},
  { "SH_X", "SH_X_BC", DIRICHLET, SH_X_BC, R_SHELL_X, SCALAR, NO_ROT, {0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB},
  { "SH_Y", "SH_Y_BC", DIRICHLET, SH_Y_BC, R_SHELL_Y, SCALAR, NO_ROT, {0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB},
  { "SH_K", "SH_K_BC", DIRICHLET, SH_K_BC, R_SHELL_CURVATURE, SCALAR, NO_ROT, {0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB},
  { "SH_TENS", "SH_TENS_BC", DIRICHLET, SH_TENS_BC, R_SHELL_TENSION, SCALAR, NO_ROT, {0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB},
  { "SH_FLUID_STRESS", "SH_FLUID_STRESS_BC", COLLOCATE_SURF, SH_FLUID_STRESS_BC, R_SHELL_CURVATURE, VECTOR, NO_ROT,{0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1}, CROSS_PHASE, DVI_SINGLE_PHASE_DB},
  { "SH_FLUID_SHEAR", "SHEAR_TO_SHELL", WEAK_INT_SURF, SHEAR_TO_SHELL_BC, R_SHELL_TENSION, SCALAR, NO_ROT,{1,1,1,0,0,1,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, CROSS_PHASE, DVI_SINGLE_PHASE_DB},
 { "SH_SLOPE_X",  "SH_SLOPE_X_BC", SPECIAL, SH_SLOPE_X_BC, R_MESH1, SCALAR, NO_ROT, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
 { "SH_SLOPE_Y",  "SH_SLOPE_Y_BC", SPECIAL, SH_SLOPE_Y_BC, R_MESH2, SCALAR, NO_ROT, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "SH_USER", "SH_USER_BC", DIRICHLET, SH_USER_BC, R_SHELL_USER, SCALAR, NO_ROT, {0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB},
  { "SH_LUBP", "SH_LUBP_BC", DIRICHLET, SH_LUBP_BC, R_SHELL_LUBP, SCALAR, NO_ROT, {0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB},
  { "LUBP_SH_FP_MATCH", "LUBP_SH_FP_MATCH_BC", STRONG_INT_SURF , LUBP_SH_FP_MATCH_BC, R_LUBP, SCALAR, NO_ROT, { 0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1}, CROSS_PHASE, DVI_SIDTIE},
  { "LUBP_SH_FP_FLUX", "LUBP_SH_FP_FLUX_BC", COLLOCATE_SURF, LUBP_SH_FP_FLUX_BC, R_SHELL_FILMP, SCALAR, NO_ROT, { 0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1} , CROSS_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "SH_LUBP_SOLID", "SH_LUBP_SOLID_BC", WEAK_INT_SURF, SH_LUBP_SOLID_BC, R_MESH1, VECTOR, NO_ROT, {0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB},
 { "SH_LUBP_SOLID_RS", "SH_LUBP_SOLID_RS_BC", WEAK_INT_SURF, SH_LUBP_SOLID_RS_BC, R_SOLID1, VECTOR, NO_ROT, {0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB},

 { "SH_S11_WEAK", "SH_S11_WEAK_BC", WEAK_INT_SURF, SH_S11_WEAK_BC, R_MESH2, SCALAR, NO_ROT,
{0,0,0,0,0,1,1,1,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,
 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,
 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,
 0,0,0,0,0,0,0,0,1,1, 1,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB},

 { "SH_S22_WEAK", "SH_S22_WEAK_BC", WEAK_INT_SURF, SH_S22_WEAK_BC, R_MESH1, SCALAR, NO_ROT,
{0,0,0,0,0,1,1,1,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,
 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,
 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,
 0,0,0,0,0,0,0,0,1,1, 1,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB},

  /*
   *  147 - SH_SURF_DIV_V
   *  147 R_SHELL_SURF_DIV_V
   */
  { "GAMMA1", "SH_GAMMA1_BC", DIRICHLET, SH_GAMMA1_BC, R_SHELL_SURF_DIV_V, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB},

  { "GAMMA1_DERIV_SYMM",  "SH_GAMMA1_DERIV_SYMM_BC", COLLOCATE_SURF, SH_GAMMA1_DERIV_SYMM_BC, R_SHELL_SURF_DIV_V, SCALAR, NO_ROT, {0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,

  { "GAMMA2_DERIV_SYMM",  "SH_GAMMA2_DERIV_SYMM_BC", COLLOCATE_SURF, SH_GAMMA2_DERIV_SYMM_BC, R_SHELL_SURF_CURV, SCALAR, NO_ROT, {0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,

  { "DVZDR_ZERO",  "DVZDR_ZERO_BC", COLLOCATE_SURF, DVZDR_ZERO_BC, R_MOMENTUM1, SCALAR, NO_ROT, {1,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,

  { "LUB_PRESS", "LUB_PRESS_BC", DIRICHLET, LUB_PRESS_BC, R_LUBP, SCALAR, NO_ROT, {0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB},

  { "LUB_PRESS_2", "LUB_PRESS_2_BC", DIRICHLET, LUB_PRESS_2_BC, R_LUBP_2, SCALAR, NO_ROT, {0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB},
  { "GRAD_LUB_PRESS", "GRAD_LUB_PRESS_BC", STRONG_INT_SURF, GRAD_LUB_PRESS_BC, R_LUBP, SCALAR, NO_ROT,
{0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1},  
SINGLE_PHASE, DVI_SINGLE_PHASE_DB},

  { "SHELL_TEMP", "SHELL_TEMP_BC", DIRICHLET, SHELL_TEMP_BC, R_SHELL_ENERGY, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB},
  { "SHELL_GRAD_TEMP", "SHELL_GRAD_TEMP_BC", STRONG_INT_SURF, SHELL_GRAD_TEMP_BC, R_SHELL_ENERGY, SCALAR, NO_ROT, 
{0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,1}, 
SINGLE_PHASE, DVI_SINGLE_PHASE_DB},
  { "SHELL_GRAD_TEMP_NOBC", "SHELL_GRAD_TEMP_NOBC_BC", WEAK_INT_SURF, SHELL_GRAD_TEMP_NOBC_BC, R_SHELL_ENERGY, SCALAR, NO_ROT, 
{0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,1}, 
SINGLE_PHASE, DVI_SINGLE_PHASE_DB},

  { "SHELL_FILMP", "SHELL_FILMP_BC", DIRICHLET, SHELL_FILMP_BC, R_SHELL_FILMP, SCALAR, NO_ROT, 
{0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0
}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB},


  { "SHELL_FILMH", "SHELL_FILMH_BC", DIRICHLET, SHELL_FILMH_BC, R_SHELL_FILMH, SCALAR, NO_ROT, 
{0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0
}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB},
 
  { "SHELL_PARTC", "SHELL_PARTC_BC", DIRICHLET, SHELL_PARTC_BC, R_SHELL_PARTC, SCALAR, NO_ROT, 
{0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1 
},SINGLE_PHASE, DVI_SINGLE_PHASE_DB},

 { "SHELL_GRAD_FP", "SHELL_GRAD_FP_BC", STRONG_INT_SURF, SHELL_GRAD_FP_BC, R_SHELL_FILMP, SCALAR, NO_ROT, 
{0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1 
}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB},

 { "SHELL_FLOW_DEVELOPED", "SHELL_FLOW_DEVELOPED_BC", STRONG_INT_SURF, SHELL_FLOW_DEVELOPED_BC, R_SHELL_FILMP, SCALAR, NO_ROT, 
{0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1 
}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB},


 { "SHELL_GRAD_FP_NOBC", "SHELL_GRAD_FP_NOBC_BC", WEAK_INT_SURF, SHELL_GRAD_FP_NOBC_BC, R_SHELL_FILMP, SCALAR, NO_ROT, 
{0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1 
}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB},

 { "SHELL_GRAD_FH", "SHELL_GRAD_FH_BC", STRONG_INT_SURF, SHELL_GRAD_FH_BC, R_SHELL_FILMH, SCALAR, NO_ROT, 
{0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0 
}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB},

 { "SHELL_GRAD_FH_NOBC", "SHELL_GRAD_FH_NOBC_BC", WEAK_INT_SURF, SHELL_GRAD_FH_NOBC_BC, R_SHELL_FILMH, SCALAR, NO_ROT, 
{0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0 
}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB},

 { "SHELL_GRAD_PC", "SHELL_GRAD_PC_BC", STRONG_INT_SURF, SHELL_GRAD_PC_BC, R_SHELL_PARTC, SCALAR, NO_ROT, 
{0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1 
}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB},

 { "SHELL_GRAD_PC_NOBC", "SHELL_GRAD_PC_NOBC_BC", WEAK_INT_SURF, SHELL_GRAD_PC_NOBC_BC, R_SHELL_PARTC, SCALAR, NO_ROT, 
{0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1 
}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB},

  { "SHELL_OPEN_PRESS", "SHELL_OPEN_PRESS_BC", DIRICHLET, SHELL_OPEN_PRESS_BC, R_SHELL_SAT_OPEN, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB},

  { "SHELL_OPEN_PRESS_2", "SHELL_OPEN_PRESS_2_BC", DIRICHLET, SHELL_OPEN_PRESS_2_BC, R_SHELL_SAT_OPEN_2, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB},

  {  "F", "F_BC", DIRICHLET, F_BC, R_FILL, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1},SINGLE_PHASE, DVI_SINGLE_PHASE_DB },
  {  "F_DIODE", "F_DIODE_BC", DIRICHLET, F_DIODE_BC, R_FILL, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1},SINGLE_PHASE, DVI_SINGLE_PHASE_DB },
  {  "H", "H_BC", DIRICHLET, H_BC, R_CURVATURE, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
  SINGLE_PHASE, DVI_SINGLE_PHASE_DB },
  {  "SH", "SH_BC", DIRICHLET, SH_BC, R_SHEAR_RATE, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB },
  
  {  "F1", "F1_BC", DIRICHLET, F1_BC, R_PHASE1, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},SINGLE_PHASE, DVI_SINGLE_PHASE_DB },
  {  "F2", "F2_BC", DIRICHLET, F2_BC, R_PHASE2, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},SINGLE_PHASE, DVI_SINGLE_PHASE_DB },
  {  "F3", "F3_BC", DIRICHLET, F3_BC, R_PHASE3, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},SINGLE_PHASE, DVI_SINGLE_PHASE_DB },
  {  "F4", "F4_BC", DIRICHLET, F4_BC, R_PHASE4, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},SINGLE_PHASE, DVI_SINGLE_PHASE_DB },
  {  "F5", "F5_BC", DIRICHLET, F5_BC, R_PHASE5, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},SINGLE_PHASE, DVI_SINGLE_PHASE_DB },

  { "APR", "APR_BC", DIRICHLET, APR_BC, R_ACOUS_PREAL, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1} , SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "API", "API_BC", DIRICHLET, API_BC, R_ACOUS_PIMAG, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1} , SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "APR_PLANE_TRANS", "APR_PLANE_TRANS_BC", WEAK_INT_SURF, APR_PLANE_TRANS_BC, R_ACOUS_PREAL, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1} , SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "API_PLANE_TRANS", "API_PLANE_TRANS_BC", WEAK_INT_SURF, API_PLANE_TRANS_BC, R_ACOUS_PIMAG, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1} , SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "INTP", "INTP_BC", DIRICHLET, INTP_BC, R_LIGHT_INTP, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1} , SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "INTM", "INTM_BC", DIRICHLET, INTM_BC, R_LIGHT_INTM, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1} , SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "INTD", "INTD_BC", DIRICHLET, INTD_BC, R_LIGHT_INTD, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1} , SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "RESTIME", "RESTIME_BC", DIRICHLET, RESTIME_BC, R_RESTIME, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1} , SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,  
  { "LIGHTP_TRANS", "LIGHTP_TRANS_BC", WEAK_INT_SURF, LIGHTP_TRANS_BC, R_LIGHT_INTP, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1} , SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "LIGHTM_TRANS", "LIGHTM_TRANS_BC", WEAK_INT_SURF, LIGHTM_TRANS_BC, R_LIGHT_INTM, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1} , SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "LIGHTD_TRANS", "LIGHTD_TRANS_BC", WEAK_INT_SURF, LIGHTD_TRANS_BC, R_LIGHT_INTD, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1} , SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "LIGHTP_JUMP",   "LIGHTP_JUMP_BC", STRONG_INT_SURF, LIGHTP_JUMP_BC,  R_LIGHT_INTP,       SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1} , CROSS_PHASE_DISCONTINUOUS,  DVI_SIDTIE } ,
  { "LIGHTM_JUMP", "LIGHTM_JUMP_BC", STRONG_INT_SURF, LIGHTM_JUMP_BC, R_LIGHT_INTM, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1} , CROSS_PHASE_DISCONTINUOUS, DVI_SIDTIE} ,
  { "LIGHTP_JUMP_2",   "LIGHTP_JUMP_2_BC", STRONG_INT_SURF, LIGHTP_JUMP_2_BC,  R_LIGHT_INTP,       SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1} , CROSS_PHASE_DISCONTINUOUS,  DVI_SIDTIE } ,
  { "LIGHTM_JUMP_2", "LIGHTM_JUMP_2_BC", WEAK_INT_SURF, LIGHTM_JUMP_2_BC, R_LIGHT_INTM, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1} , CROSS_PHASE_DISCONTINUOUS, DVI_SIDTIE} ,
  { "APR_NOBC", "APR_NOBC_BC", WEAK_INT_SURF, APR_NOBC_BC, R_ACOUS_PREAL, SCALAR, NO_ROT, {0,0,0,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1} , SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "API_NOBC", "API_NOBC_BC", WEAK_INT_SURF, API_NOBC_BC, R_ACOUS_PIMAG, SCALAR, NO_ROT, {0,0,0,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1} , SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "APR_VELOCITY", "APR_VELOCITY_BC", WEAK_INT_SURF, APR_VELOCITY_BC, R_ACOUS_PREAL, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1} , SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "API_VELOCITY", "API_VELOCITY_BC", WEAK_INT_SURF, API_VELOCITY_BC, R_ACOUS_PIMAG, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1} , SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "EM_E1R", "EM_E1R_BC", DIRICHLET, EM_E1R_BC, R_EM_E1_REAL, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1} , SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "EM_E2R", "EM_E2R_BC", DIRICHLET, EM_E2R_BC, R_EM_E2_REAL, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1} , SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "EM_E3R", "EM_E3R_BC", DIRICHLET, EM_E3R_BC, R_EM_E3_REAL, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1} , SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "EM_E1I", "EM_E1I_BC", DIRICHLET, EM_E1I_BC, R_EM_E1_IMAG, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1} , SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "EM_E2I", "EM_E2I_BC", DIRICHLET, EM_E2I_BC, R_EM_E2_IMAG, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1} , SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "EM_E3I", "EM_E3I_BC", DIRICHLET, EM_E3I_BC, R_EM_E3_IMAG, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1} , SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "EM_H1R", "EM_H1R_BC", DIRICHLET, EM_H1R_BC, R_EM_H1_REAL, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1} , SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "EM_H2R", "EM_H2R_BC", DIRICHLET, EM_H2R_BC, R_EM_H2_REAL, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1} , SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "EM_H3R", "EM_H3R_BC", DIRICHLET, EM_H3R_BC, R_EM_H3_REAL, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1} , SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "EM_H1I", "EM_H1I_BC", DIRICHLET, EM_H1I_BC, R_EM_H1_IMAG, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1} , SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "EM_H2I", "EM_H2I_BC", DIRICHLET, EM_H2I_BC, R_EM_H2_IMAG, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1} , SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "EM_H3I", "EM_H3I_BC", DIRICHLET, EM_H3I_BC, R_EM_H3_IMAG, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1} , SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,

  { "NN", "NN_BC", DIRICHLET, NN_BC, R_BOND_EVOLUTION, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1} , SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "T", "T_BC", DIRICHLET, T_BC, R_ENERGY, SCALAR, NO_ROT, {0, 0, 0, 1, 0, 0, 0, 0, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "U", "U_BC", DIRICHLET, U_BC, R_MOMENTUM1, SCALAR, NO_ROT, {1, 0, 0, 0, 0, 0, 0, 0, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "V", "V_BC", DIRICHLET, V_BC, R_MOMENTUM2, SCALAR, NO_ROT, {0, 1, 0, 0, 0, 0, 0, 0, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "W", "W_BC", DIRICHLET, W_BC, R_MOMENTUM3, SCALAR, NO_ROT, {0, 0, 1, 0, 0, 0, 0, 0, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,

  { "EXT_V", "EXT_V_BC", DIRICHLET, EXT_V_BC, R_EXT_VELOCITY, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0, 0, 1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "LM1", "LM1_BC", DIRICHLET, LM1_BC, R_LAGR_MULT1, SCALAR, NO_ROT, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "LM2", "LM2_BC", DIRICHLET, LM2_BC, R_LAGR_MULT2, SCALAR, NO_ROT, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "LM3", "LM3_BC", DIRICHLET, LM3_BC, R_LAGR_MULT3, SCALAR, NO_ROT, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,

  { "UVARY", "UVARY_BC", COLLOCATE_SURF, UVARY_BC, R_MOMENTUM1, SCALAR, NO_ROT, {1, 0, 0, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "VVARY", "VVARY_BC", COLLOCATE_SURF, VVARY_BC, R_MOMENTUM2, SCALAR, NO_ROT, {0, 1, 0, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "WVARY", "WVARY_BC", COLLOCATE_SURF, WVARY_BC, R_MOMENTUM3, SCALAR, NO_ROT, {0, 0, 1, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "U_PARABOLA", "U_PARABOLA_BC", COLLOCATE_SURF, U_PARABOLA_BC, R_MOMENTUM1, SCALAR, NO_ROT, {1, 0, 0, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "V_PARABOLA", "V_PARABOLA_BC", COLLOCATE_SURF, V_PARABOLA_BC, R_MOMENTUM2, SCALAR, NO_ROT, {0, 1, 0, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "W_PARABOLA", "W_PARABOLA_BC", COLLOCATE_SURF, W_PARABOLA_BC, R_MOMENTUM3, SCALAR, NO_ROT, {0, 0, 1, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,

  { "T_USER", "T_USER_BC", COLLOCATE_SURF, T_USER_BC, R_ENERGY, SCALAR, NO_ROT, {0, 0, 0, 1, 0, 0, 0, 0, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,

  { "UUSER", "UUSER_BC", STRONG_INT_SURF, UUSER_BC, R_MOMENTUM1, SCALAR, NO_ROT, {1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "VUSER", "VUSER_BC", STRONG_INT_SURF, VUSER_BC, R_MOMENTUM2, SCALAR, NO_ROT, {1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "WUSER", "WUSER_BC", STRONG_INT_SURF, WUSER_BC, R_MOMENTUM3, SCALAR, NO_ROT, {1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "UUSER_COLLOC", "UUSER_COLLOC_BC", COLLOCATE_SURF, UUSER_COLLOC_BC, R_MOMENTUM1, SCALAR, NO_ROT, {1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "VUSER_COLLOC", "VUSER_COLLOC_BC", COLLOCATE_SURF, VUSER_COLLOC_BC, R_MOMENTUM2, SCALAR, NO_ROT, {1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "WUSER_COLLOC", "WUSER_COLLOC_BC", COLLOCATE_SURF, WUSER_COLLOC_BC, R_MOMENTUM3, SCALAR, NO_ROT, {1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,

  { "DX_USER", "DX_USER_BC", COLLOCATE_SURF, DX_USER_BC, R_MESH1, SCALAR, NO_ROT, {1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "DY_USER", "DY_USER_BC", COLLOCATE_SURF, DY_USER_BC, R_MESH2, SCALAR, NO_ROT, {1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "DZ_USER", "DZ_USER_BC", COLLOCATE_SURF, DZ_USER_BC, R_MESH3, SCALAR, NO_ROT, {1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "DX_USER_NODE", "DX_USER_NODE_BC", DIRICHLET, DX_USER_NODE_BC, R_MESH1, SCALAR, NO_ROT, {0, 0, 0, 0, 0, 1, 0, 0, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "DY_USER_NODE", "DY_USER_NODE_BC", DIRICHLET, DY_USER_NODE_BC, R_MESH2, SCALAR, NO_ROT, {0, 0, 0, 0, 0, 0, 1, 0, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "DZ_USER_NODE", "DZ_USER_NODE_BC", DIRICHLET, DZ_USER_NODE_BC, R_MESH3, SCALAR, NO_ROT, {0, 0, 0, 0, 0, 0, 0, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "TENSION_SHEET",    "TENSION_SHEET_BC", STRONG_INT_SURF, TENSION_SHEET_BC,   R_MESH_NORMAL, SCALAR, R_MESH1, {1, 1, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "KIN_LEAK",     "KIN_LEAK_BC", STRONG_INT_SURF, KIN_LEAK_BC,    R_MESH_NORMAL, SCALAR, R_MESH1, {1, 1, 1, 1, 1, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "KIN_LEAK_HEAT",     "KIN_LEAK_HEAT_BC", STRONG_INT_SURF, KIN_LEAK_HEAT_BC,    R_MESH_NORMAL, SCALAR, R_MESH1, {1, 1, 1, 1, 1, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
 { "KIN_CHEM",     "KIN_CHEM_BC",     STRONG_INT_SURF, KIN_CHEM_BC,    R_MESH_NORMAL, SCALAR, R_MESH1,     {1, 1, 1, 1, 1, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB  } ,
  { "KIN_ELECTRODEPOSITION",     "KIN_ELECTRODEPOSITION_BC",     STRONG_INT_SURF, KIN_ELECTRODEPOSITION_BC,    R_MESH_NORMAL, SCALAR, R_MESH1,     {0,0,0,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0}, SINGLE_PHASE } ,  /*  RSL 5/27/02  */
  { "KINEMATIC",    "KINEMATIC_BC", STRONG_INT_SURF, KINEMATIC_BC,   R_MESH_NORMAL, SCALAR, R_MESH1, {1, 1, 1, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "KINEMATIC_COLLOC",    "KINEMATIC_COLLOC_BC", COLLOCATE_SURF, KINEMATIC_COLLOC_BC,   R_MESH_NORMAL, SCALAR, R_MESH1, {1, 1, 1, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "KINEMATIC_PETROV",    "KINEMATIC_PETROV_BC", STRONG_INT_SURF, KINEMATIC_PETROV_BC,   R_MESH_NORMAL, SCALAR, R_MESH1, {1, 1, 1, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "KINEMATIC_SPECIES",    "KINEMATIC_SPECIES_BC", WEAK_INT_SURF, KINEMATIC_SPECIES_BC,   R_MASS, SCALAR, NO_ROT, {1, 1, 1, 1, 1, 1, 1, 1, 0, 0}, CROSS_PHASE_DISCONTINUOUS, DVI_DVVSIG } ,
  { "KINEMATIC_DISC",    "KINEMATIC_DISC_BC", STRONG_INT_SURF, KINEMATIC_DISC_BC,   R_MESH_NORMAL, SCALAR, R_MESH1, {1, 1, 1, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "KINEMATIC_EDGE",    "KINEMATIC_EDGE_BC", STRONG_INT_EDGE, KINEMATIC_EDGE_BC,   R_MESH_NORMAL, SCALAR, R_MESH1, {1, 1, 1, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "VNORM_LEAK",   "VNORM_LEAK_BC", STRONG_INT_SURF, VNORM_LEAK_BC,   R_MOM_NORMAL, SCALAR, R_MOMENTUM1, {1, 1, 1, 1, 1, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "VNORM_ELECTRODEPOSITION",   "VNORM_ELECTRODEPOSITION_BC",   STRONG_INT_SURF, VNORM_ELECTRODEPOSITION_BC,   R_MOM_NORMAL, SCALAR, R_MOMENTUM1, {1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0}, SINGLE_PHASE } ,  /*  RSL 5/30/02  */
  { "VELO_NORMAL",  "VELO_NORMAL_BC",  STRONG_INT_SURF, VELO_NORMAL_BC,  R_MOM_NORMAL, SCALAR, R_MOMENTUM1, {1, 1, 1, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
 { "VELO_NORMAL_LS",  "VELO_NORMAL_LS_BC",  STRONG_INT_SURF, VELO_NORMAL_LS_BC,  R_MOM_NORMAL, SCALAR, R_MOMENTUM1, {1, 1, 1, 0, 0, 1, 1, 1, 0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
 { "VELO_NORMAL_LS_COLLOC",  "VELO_NORMAL_LS_COLLOC_BC",  COLLOCATE_SURF, VELO_NORMAL_LS_COLLOC_BC,  R_MOM_NORMAL, SCALAR, R_MOMENTUM1, {1, 1, 1, 0, 0, 1, 1, 1, 0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
 { "VELO_NORMAL_LS_PETROV",  "VELO_NORMAL_LS_PETROV_BC",  STRONG_INT_SURF, VELO_NORMAL_LS_PETROV_BC,  R_MOM_NORMAL, SCALAR, R_MOMENTUM1, {1, 1, 1, 0, 0, 1, 1, 1, 0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "LS_ATTACH",  "LS_SUCKS",  STRONG_INT_SURF, LS_ATTACH_BC,  R_MOM_NORMAL, SCALAR, R_MOMENTUM1, {1, 1, 1, 0, 0, 1, 1, 1, 0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "VELO_NORM_COLLOC",  "VELO_NORM_COLLOC_BC",  COLLOCATE_SURF, VELO_NORM_COLLOC_BC,  R_MOM_NORMAL, SCALAR, R_MOMENTUM1, {1, 1, 1, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "VELO_NORMAL_EDGE_INT",  "VELO_NORMAL_EDGE_INT_BC",  STRONG_INT_EDGE, VELO_NORMAL_EDGE_INT_BC,   R_MOM_NORMAL, SCALAR, R_MOMENTUM1, {1, 1, 1, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "VELO_NORMAL_DISC",  "VELO_NORMAL_DISC_BC",  STRONG_INT_SURF, VELO_NORMAL_DISC_BC,  R_MOM_NORMAL, SCALAR, R_MOMENTUM1, {1, 1, 1, 0, 0, 1, 1, 1, 0, 0}, CROSS_PHASE_DISCONTINUOUS, DVI_SIDTIE } ,
  { "VELO_NORMAL_SOLID", "VELO_NORMAL_SOLID_BC", STRONG_INT_SURF, VELO_NORMAL_SOLID_BC, R_MOM_NORMAL,  SCALAR, R_MOMENTUM2, {1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 1,1,1,1 }, CROSS_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "VELO_NORMAL_EDGE",  "VELO_NORMAL_EDGE_BC",  COLLOCATE_EDGE, VELO_NORMAL_EDGE_BC,   R_MOM_NORMAL, SCALAR, R_MOMENTUM1, {1, 1, 1, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "VELO_TANGENT", "VELO_TANGENT_BC", STRONG_INT_SURF, VELO_TANGENT_BC, R_MOM_TANG1,  SCALAR, R_MOMENTUM1, {1, 1, 1, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "VELO_TANGENT_LS",  "VELO_TANGENT_LS_BC",  STRONG_INT_SURF, VELO_TANGENT_LS_BC,  R_MOM_TANG1, SCALAR, R_MOMENTUM1, {1, 1, 1, 0, 0, 1, 1, 1, 0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "VELO_TANGENT_3D", "VELO_TANGENT_3D_BC", STRONG_INT_SURF, VELO_TANGENT_3D_BC, R_MOM_TANG1,  SCALAR, R_MOMENTUM1, {1, 1, 1, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "VELO_TANGENT_SOLID", "VELO_TANGENT_SOLID_BC", STRONG_INT_SURF, VELO_TANGENT_SOLID_BC, R_MOM_TANG1,  SCALAR, R_MOMENTUM1, {1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 1,1,1,1 }, CROSS_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "VELO_TANGENT_EDGE_INT", "VELO_TANGENT_EDGE_INT_BC", STRONG_INT_EDGE, VELO_TANGENT_EDGE_INT_BC, R_MOM_TANG1,  SCALAR, R_MOMENTUM1, {1, 1, 1, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "VELO_TANGENT_EDGE", "VELO_TANGENT_EDGE_BC", COLLOCATE_EDGE, VELO_TANGENT_EDGE_BC, R_MOM_TANG1,  SCALAR, R_MOMENTUM1, {1, 1, 1, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "VELO_TANGENT_USER", "VELO_TANGENT_USER_BC", STRONG_INT_SURF, VELO_TANGENT_USER_BC, R_MOM_TANG1,  SCALAR, R_MOMENTUM1, {1, 1, 1, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "VELO_SLIP",    "VELO_SLIP_BC",    WEAK_INT_SURF, VELO_SLIP_BC,    R_MOMENTUM1,  VECTOR, NO_ROT, {1, 1, 1, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
{ "VELO_SLIP_POWER_CARD",    "VELO_SLIP_POWER_CARD_BC",    WEAK_INT_SURF, VELO_SLIP_POWER_CARD_BC,    R_MOM_TANG1,  VECTOR, R_MOMENTUM1, {1, 1, 1, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,

{ "VELO_SLIP_POWER",    "VELO_SLIP_POWER_BC",    WEAK_INT_SURF, VELO_SLIP_POWER_BC,    R_MOMENTUM1,  VECTOR, NO_ROT, {1, 1, 1, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "AIR_FILM",    "AIR_FILM_BC",    WEAK_INT_SURF, AIR_FILM_BC,    R_MOMENTUM1,  VECTOR, NO_ROT, {1, 1, 1, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "VELO_SLIP_FLUID",    "VELO_SLIP_FLUID_BC",    WEAK_INT_SURF, VELO_SLIP_FLUID_BC,    R_MOMENTUM1,  VECTOR, NO_ROT, {1, 1, 1, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
#ifdef COUPLED_FILL
  { "VELO_SLIP_FILL",    "VELO_SLIP_FILL_BC",    WEAK_INT_SURF, VELO_SLIP_FILL_BC,    R_MOMENTUM1,  VECTOR, NO_ROT, {1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "VELO_SLIP_ROT_FILL",    "VELO_SLIP_ROT_FILL_BC",    WEAK_INT_SURF, VELO_SLIP_ROT_FILL_BC,    R_MOMENTUM1,  VECTOR, NO_ROT, {1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "VELO_SLIP_LS",    "VELO_SLIP_LS_BC",    WEAK_INT_SURF, VELO_SLIP_LEVEL_BC,    R_MOMENTUM1,  VECTOR, NO_ROT, {1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "VELO_SLIP_LS_SIC",    "VELO_SLIP_LS_SIC_BC",    STRONG_INT_SURF, VELO_SLIP_LEVEL_SIC_BC,    R_MOMENTUM1,  VECTOR, R_MOMENTUM1, {1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
/*  { "VELO_SLIP_LS_SIC",    "VELO_SLIP_LS_SIC_BC",    STRONG_INT_SURF, VELO_SLIP_LEVEL_SIC_BC,    R_MOM_TANG1,  SCALAR, R_MOMENTUM1, {1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,  */
  { "VELO_SLIP_LS_ROT",    "VELO_SLIP_LS_ROT_BC",    WEAK_INT_SURF, VELO_SLIP_LS_ROT_BC,    R_MOMENTUM1,  VECTOR, NO_ROT, {1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
#else /* COUPLED_FILL */
  { "VELO_SLIP_FILL",    "VELO_SLIP_FILL_BC",    WEAK_INT_SURF, VELO_SLIP_FILL_BC,    R_MOMENTUM1,  VECTOR, NO_ROT, {1, 1, 1, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "VELO_SLIP_LS",    "VELO_SLIP_LS_BC",    WEAK_INT_SURF, VELO_SLIP_LEVEL_BC,    R_MOMENTUM1,  VECTOR, NO_ROT, {1, 1, 1, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "VELO_SLIP_LS_SIC",    "VELO_SLIP_LS_SIC_BC",    STRONG_INT_SURF, VELO_SLIP_LEVEL_SIC_BC,    R_MOMENTUM1,  VECTOR, NO_ROT, {1, 1, 1, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "VELO_SLIP_LS_ROT",    "VELO_SLIP_LS_ROT_BC",    WEAK_INT_SURF, VELO_SLIP_LS_ROT_BC,    R_MOMENTUM1,  VECTOR, NO_ROT, {1, 1, 1, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "VELO_SLIP_ROT_FILL",    "VELO_SLIP_ROT_FILL_BC",    WEAK_INT_SURF, VELO_SLIP_ROT_FILL_BC,    R_MOMENTUM1,  VECTOR, NO_ROT, {1, 1, 1, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
#endif /* COUPLED_FILL */
  { "VELO_SLIP_ROT",    "VELO_SLIP_ROT_BC",   WEAK_INT_SURF, VELO_SLIP_ROT_BC,    R_MOMENTUM1,  VECTOR, NO_ROT, {1, 1, 1, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "AIR_FILM_ROT",    "AIR_FILM_ROT_BC",   WEAK_INT_SURF, AIR_FILM_ROT_BC,    R_MOMENTUM1,  VECTOR, NO_ROT, {1, 1, 1, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "VELO_SLIP_ROT_FLUID",    "VELO_SLIP_ROT_FLUID_BC",   WEAK_INT_SURF, VELO_SLIP_ROT_FLUID_BC,    R_MOMENTUM1,  VECTOR, NO_ROT, {1, 1, 1, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "VELO_SLIP_ELECTROKINETIC", "VELO_SLIP_EK_BC", STRONG_INT_SURF, VELO_SLIP_EK_BC, R_MOM_TANG1,  SCALAR, R_MOMENTUM1, {1, 1, 1, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "VELO_SLIP_ELECTROKINETIC3D", "VELO_EK_3D_BC", STRONG_INT_SURF, VELO_EK_3D_BC, R_MOM_TANG1,  SCALAR, R_MOMENTUM1, {1, 1, 1, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "VELO_SLIP_SOLID", "VELO_SLIP_SOLID_BC", WEAK_INT_SURF, VELO_SLIP_SOLID_BC, R_MOMENTUM1,  VECTOR, NO_ROT, {1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 1,1,1,1 }, CROSS_PHASE, DVI_SINGLE_PHASE_DB } ,
{ "WETTING_SPEED_LINEAR",    "WETTING_SPEED_LINEAR_BC",    WEAK_INT_SURF, WETTING_SPEED_LIN_BC,    R_MOMENTUM1,  VECTOR, NO_ROT, {1,1,1,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
{ "WETTING_SPEED_BLAKE",    "WETTING_SPEED_BLAKE_BC",    WEAK_INT_SURF, WETTING_SPEED_BLAKE_BC,    R_MOMENTUM1,  VECTOR, NO_ROT, {1,1,1,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
{ "WETTING_SPEED_HOFFMAN",    "WETTING_SPEED_HOFFMAN_BC",    WEAK_INT_SURF, WETTING_SPEED_HOFFMAN_BC,    R_MOMENTUM1,  VECTOR, NO_ROT, {1,1,1,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
{ "WETTING_SPEED_COX",    "WETTING_SPEED_COX_BC",    WEAK_INT_SURF, WETTING_SPEED_COX_BC,    R_MOMENTUM1,  VECTOR, NO_ROT, {1,1,1,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
{ "WETTING_SPEED_SHIK",    "WETTING_SPEED_SHIK_BC",    WEAK_INT_SURF, WETTING_SPEED_SHIK_BC,    R_MOMENTUM1,  VECTOR, NO_ROT, {1,1,1,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
{ "BLAKE_DIRICHLET",    "BLAKE_DIRICHLET_BC",    WEAK_INT_SURF, BLAKE_DIRICHLET_BC,    R_MOMENTUM1,  VECTOR, NO_ROT, {1,1,1,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
{ "HOFFMAN_DIRICHLET",    "HOFFMAN_DIRICHLET_BC",    WEAK_INT_SURF, HOFFMAN_DIRICHLET_BC,    R_MOMENTUM1,  VECTOR, NO_ROT, {1,1,1,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
{ "COX_DIRICHLET",    "COX_DIRICHLET_BC",    WEAK_INT_SURF, COX_DIRICHLET_BC,    R_MOMENTUM1,  VECTOR, NO_ROT, {1,1,1,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
{ "SHIK_DIRICHLET",    "SHIK_DIRICHLET_BC",    WEAK_INT_SURF, SHIK_DIRICHLET_BC,    R_MOMENTUM1,  VECTOR, NO_ROT, {1,1,1,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
{ "BLAKE_DIRICH_ROLL",    "BLAKE_DIRICH_ROLL_BC",    WEAK_INT_SURF, BLAKE_DIRICH_ROLL_BC,    R_MOMENTUM1,  VECTOR, NO_ROT, {1,1,1,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
{ "HOFFMAN_DIRICH_ROLL",    "HOFFMAN_DIRICH_ROLL_BC",    WEAK_INT_SURF, HOFFMAN_DIRICH_ROLL_BC,    R_MOMENTUM1,  VECTOR, NO_ROT, {1,1,1,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
{ "COX_DIRICH_ROLL",    "COX_DIRICH_ROLL_BC",    WEAK_INT_SURF, COX_DIRICH_ROLL_BC,    R_MOMENTUM1,  VECTOR, NO_ROT, {1,1,1,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
{ "SHIK_DIRICH_ROLL",    "SHIK_DIRICH_ROLL_BC",    WEAK_INT_SURF, SHIK_DIRICH_ROLL_BC,    R_MOMENTUM1,  VECTOR, NO_ROT, {1,1,1,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
{ "LINEAR_WETTING_SIC",    "LINEAR_WETTING_SIC_BC",    WEAK_INT_SURF, LINEAR_WETTING_SIC_BC,    R_MOMENTUM1,  VECTOR, NO_ROT, {1,1,1,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
{ "HYSTERESIS_WETTING",    "HYSTERESIS_WETTING_BC",    WEAK_INT_SURF, HYSTERESIS_WETTING_BC,    R_MOMENTUM1,  VECTOR, NO_ROT, {1,1,1,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
{ "VELO_STREAMING", "VELO_STREAMING_BC", STRONG_INT_SURF, VELO_STREAMING_BC, R_MOM_TANG1,  SCALAR, R_MOMENTUM1, {1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB} ,

  { "Y", "Y_BC", DIRICHLET, Y_BC, R_MASS, SCALAR, NO_ROT, {0, 0, 0, 0, 1, 0, 0, 0, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,

  { "YFLUX",       "YFLUX_BC",       WEAK_INT_SURF, YFLUX_BC,       R_MASS, SCALAR, NO_ROT, {1, 1, 1, 1, 1, 1, 1, 1, 0, 1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
 { "YFLUX_BV",  "YFLUX_BV_BC",  WEAK_INT_SURF, YFLUX_BV_BC, R_MASS, SCALAR, NO_ROT, {1, 1, 1, 1, 1, 1, 1, 1, 0, 1},SINGLE_PHASE, DVI_SINGLE_PHASE_DB} ,
 { "YFLUX_HOR",  "YFLUX_HOR_BC",  WEAK_INT_SURF, YFLUX_HOR_BC, R_MASS, SCALAR, NO_ROT, {0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0},SINGLE_PHASE } ,
 { "YFLUX_ORR",  "YFLUX_ORR_BC",  WEAK_INT_SURF, YFLUX_ORR_BC, R_MASS, SCALAR, NO_ROT, {0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0},SINGLE_PHASE } ,
 { "YFLUX_H2O_ANODE",  "YFLUX_H2O_ANODE_BC",  WEAK_INT_SURF, YFLUX_H2O_ANODE_BC, R_MASS, SCALAR, NO_ROT, {0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0},SINGLE_PHASE } ,
 { "YFLUX_H2O_CATHODE",  "YFLUX_H2O_CATHODE_BC",  WEAK_INT_SURF, YFLUX_H2O_CATHODE_BC, R_MASS, SCALAR, NO_ROT, {0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0},SINGLE_PHASE } ,
 { "YFLUX_SULFIDATION",  "YFLUX_SULFIDATION_BC",  WEAK_INT_SURF, YFLUX_SULFIDATION_BC, R_MASS, SCALAR, NO_ROT, {1, 1, 1, 1, 1, 1, 1, 1, 0, 1},SINGLE_PHASE, DVI_SINGLE_PHASE_DB} ,
  { "YFLUX_BV2",  "YFLUX_BV2_BC",  WEAK_INT_SURF, YFLUX_BV2_BC, R_MASS, SCALAR, NO_ROT, {1,1,1,1,1,1,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0},SINGLE_PHASE } ,  /* RSL 3/9/01 */
  { "YFLUX_NI",  "YFLUX_NI_BC",  WEAK_INT_SURF, YFLUX_NI_BC, R_MASS, SCALAR, NO_ROT, {1,1,1,1,1,1,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0},SINGLE_PHASE } ,  /* RSL 3/9/01 */
  { "YTOTALFLUX_CONST", "YTOTALFLUX_CONST_BC", WEAK_INT_SURF, YTOTALFLUX_CONST_BC, R_MASS, SCALAR, NO_ROT, {1, 1, 1, 1, 1, 1, 1, 1, 0, 1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "YFLUX_EQUIL",  "YFLUX_EQUIL_BC", WEAK_INT_SURF, YFLUX_EQUIL_BC, R_MASS, SCALAR, NO_ROT, {1, 1, 1, 1, 1, 1, 1, 1, 0, 1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "YFLUX_CONST", "YFLUX_CONST_BC", WEAK_INT_SURF, YFLUX_CONST_BC, R_MASS, SCALAR, NO_ROT, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "YFLUX_USER",  "YFLUX_USER_BC",  WEAK_INT_SURF, YFLUX_USER_BC, R_MASS, SCALAR, NO_ROT, {1, 1, 1, 1, 1, 1, 1, 1, 0, 1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "YFLUX_ALLOY",  "YFLUX_ALLOY_BC",  WEAK_INT_SURF, YFLUX_ALLOY_BC, R_MASS, SCALAR, NO_ROT, {1, 1, 1, 1, 1, 1, 1, 1, 0, 1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "YFLUX_ETCH", "YFLUX_ETCH_BC", STRONG_INT_SURF, YFLUX_ETCH_BC, R_MASS, SCALAR, NO_ROT, {1, 1, 1, 1, 1, 1, 1, 1, 0, 1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "YUSER", "YUSER_BC", STRONG_INT_SURF, YUSER_BC, R_MASS, SCALAR, NO_ROT, {1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "FCD_CONC_GRAD","FCD_CONC_GRAD_BC", WEAK_INT_SURF, FICK_CHRGD_SURF_GRAD_BC,  R_MASS, SCALAR, NO_ROT, {0, 0, 0, 1, 1, 1, 1, 1, 0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,  0,0,1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "SURFACE_CHARGE", "SURFACE_CHARGE_BC", STRONG_INT_SURF, SURFACE_CHARGE_BC, R_MASS, SCALAR, NO_ROT, {1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,

 { "POROUS_FLUX", "POROUS_FLUX_BC", WEAK_INT_SURF, POROUS_FLUX_BC, R_POR_LIQ_PRES, SCALAR, NO_ROT, {1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 1,1,1,1 },SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "POROUS_CONV", "POROUS_CONV_BC", WEAK_INT_SURF, POROUS_CONV_BC, R_POR_LIQ_PRES, SCALAR, NO_ROT, {1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 1,1,1,1 },SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "POROUS_PRESSURE",  "POROUS_PRESSURE_BC", COLLOCATE_SURF, POROUS_PRESSURE_BC, R_POR_LIQ_PRES, SCALAR, NO_ROT, {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 1,1,1,1 }, CROSS_PHASE, DVI_SINGLE_PHASE_DB } ,
{ "POROUS_PRESSURE_LUB",  "POROUS_PRESSURE_LUB_BC", COLLOCATE_SURF, POROUS_PRESSURE_LUB_BC, R_POR_LIQ_PRES, SCALAR, NO_ROT, {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 1,1,1,1,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 1,1 }, CROSS_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "DARCY_CONTINUOUS",  "DARCY_CONTINUOUS_BC", STRONG_INT_SURF, DARCY_CONTINUOUS_BC, R_MOM_NORMAL, SCALAR, R_MOMENTUM1, {1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 1,1,1,1 }, CROSS_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "DARCY_LUB",  "DARCY_LUB_BC", WEAK_INT_SURF, DARCY_LUB_BC, R_LUBP, SCALAR, NO_ROT, {1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 1,1,1,1 }, CROSS_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "POROUS_LIQ_PRESSURE", "POR_LIQ_PRESSURE_BC", DIRICHLET, POROUS_LIQ_PRESSURE_BC, R_POR_LIQ_PRES, SCALAR, NO_ROT, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "POROUS_LIQ_PRESSURE_FILL", "POROUS_LIQ_PRESSURE_FILL_BC", COLLOCATE_SURF, POROUS_LIQ_PRESSURE_FILL_BC, R_POR_LIQ_PRES, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0, 0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,1,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 1,0,0,0 }, CROSS_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "POROUS_GAS_PRESSURE", "POR_GAS_PRESSURE_BC", DIRICHLET, POROUS_GAS_PRESSURE_BC, R_POR_GAS_PRES, SCALAR, NO_ROT, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "POROUS_TEMPERATURE", "POR_TEMP_BC", DIRICHLET, POROUS_TEMP_BC, R_POR_ENERGY, SCALAR, NO_ROT, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "POROUS_LIQ_FLUX_CONST", "POROUS_LIQ_FLUX_CONST_BC", WEAK_INT_SURF, POROUS_LIQ_FLUX_CONST_BC, R_POR_LIQ_PRES, SCALAR, NO_ROT, {1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 1,1,1,1 }, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "POROUS_GAS_FLUX_CONST", "POROUS_GAS_FLUX_CONST_BC", WEAK_INT_SURF, POROUS_GAS_FLUX_CONST_BC, R_POR_GAS_PRES, SCALAR, NO_ROT, {1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 1,1,1,1 }, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "POR_LIQ_FLUX_FILL", "POR_LIQ_FLUX_FILL_BC", WEAK_INT_SURF, POR_LIQ_FLUX_FILL_BC, R_POR_LIQ_PRES, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0, 0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,1,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 1,0,0,0 }, CROSS_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "P_LIQ_USER", "P_LIQ_USER_BC", COLLOCATE_SURF, P_LIQ_USER_BC, R_POR_LIQ_PRES, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0, 0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,1,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 1,0,0,0 } , SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "POROUS_SINK", "POROUS_SINK_BC", DIRICHLET, POROUS_SINK_BC, R_POR_SINK_MASS, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0, 0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0 } , SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
 { "SH_P_OPEN_USER", "SH_P_OPEN_USER_BC", COLLOCATE_SURF, SH_P_OPEN_USER_BC, R_SHELL_SAT_OPEN, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "H_FREE", "H_FREE_BC", WEAK_INT_SURF, H_FREE_BC, R_CURVATURE, SCALAR, NO_ROT, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0 }, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "LS_CA_H", "LS_CA_H_BC", WEAK_INT_SURF, LS_CA_H_BC, R_CURVATURE, SCALAR, NO_ROT, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0 }, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "YFLUX_SUS", "YFLUX_SUS_BC", WEAK_INT_SURF, YFLUX_SUS_BC, R_MASS,SCALAR, NO_ROT, {1,1,1,1,1,1,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,

  { "CAPILLARY", "CAPILLARY_BC", WEAK_INT_SURF, CAPILLARY_BC, R_MOMENTUM1, VECTOR, NO_ROT, {0, 0, 0, 1, 1, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "CAPILLARY_SHEAR_VISC", "CAPILLARY_SHEAR_VISC_BC", WEAK_INT_SURF, CAPILLARY_SHEAR_VISC_BC, R_MOMENTUM1, VECTOR, NO_ROT, {0, 0, 0, 1, 1, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "ELEC_TRACTION", "ELEC_TRACTION_BC", WEAK_INT_SURF, ELEC_TRACTION_BC, R_MOMENTUM1, VECTOR, NO_ROT, {0,0,0,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1}, CROSS_PHASE, DVI_CROSS_PHASE_CONJUGATE }, 
  { "ELEC_TRACTION_SOLID", "ELEC_TRACTION_SOLID_BC", WEAK_INT_SURF, ELEC_TRACTION_SOLID_BC, R_MESH1, VECTOR, NO_ROT, {0,0,0,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1}, CROSS_PHASE, DVI_CROSS_PHASE_CONJUGATE }, 
  { "CAP_REPULSE", "CAP_REPULSE_BC", WEAK_INT_SURF, CAP_REPULSE_BC, R_MOMENTUM1, VECTOR, NO_ROT, {0, 0, 0, 1, 1, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
 { "CAP_RECOIL_PRESS", "CAP_RECOIL_PRESS_BC", WEAK_INT_SURF, CAP_RECOIL_PRESS_BC, R_MOMENTUM1, VECTOR, NO_ROT, {0, 0, 0, 1, 1, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "CAP_REPULSE_ROLL", "CAP_REPULSE_ROLL_BC", WEAK_INT_SURF, CAP_REPULSE_ROLL_BC, R_MOMENTUM1, VECTOR, NO_ROT, {0, 0, 0, 1, 1, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "CAP_REPULSE_USER", "CAP_REPULSE_USER_BC", WEAK_INT_SURF, CAP_REPULSE_USER_BC, R_MOMENTUM1, VECTOR, NO_ROT, {0, 0, 0, 1, 1, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "CAP_REPULSE_TABLE", "CAP_REPULSE_TABLE_BC", WEAK_INT_SURF, CAP_REPULSE_TABLE_BC, R_MOMENTUM1, VECTOR, NO_ROT, {0, 0, 0, 1, 1, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "CAPILLARY_TABLE", "CAPILLARY_TABLE_BC", WEAK_INT_SURF, CAPILLARY_TABLE_BC, R_MOMENTUM1, VECTOR, NO_ROT, {0, 0, 0, 1, 1, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "PRESSURE_USER", "PRESSURE_USER_BC", WEAK_INT_SURF, PRESSURE_USER_BC, R_MOMENTUM1, VECTOR, NO_ROT, {1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "FLOW_PRESSURE",    "FLOW_PRESSURE_BC",    WEAK_INT_SURF, FLOW_PRESSURE_BC, R_MOMENTUM1, VECTOR,   
      NO_ROT, {0, 0, 0, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  {  "FLOW_HYDROSTATIC", "FLOW_HYDROSTATIC_BC", WEAK_INT_SURF, FLOW_HYDROSTATIC_BC, R_MOMENTUM1, VECTOR, 
      NO_ROT, {0, 0, 0, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  {  "FLOW_PRESSURE_VAR", "FLOW_PRESSURE_VAR_BC", WEAK_INT_SURF, FLOW_PRESSURE_VAR_BC, R_MOMENTUM1, VECTOR, 
      NO_ROT, {0, 0, 0, 1, 1, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  {  "FLOW_PRESS_USER",  "FLOW_PRESS_USER_BC",  WEAK_INT_SURF, FLOW_PRESS_USER_BC,  R_MOMENTUM1, VECTOR, 
      NO_ROT, {1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  {  "FLOW_STRESSNOBC",  "FLOW_STRESSNOBC_BC",  WEAK_INT_SURF, FLOW_STRESSNOBC_BC,  R_MOMENTUM1, VECTOR, 
      NO_ROT, {1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  {  "FLOW_GRADV",  "FLOW_GRADV_BC",  WEAK_INT_SURF, FLOW_GRADV_BC,  R_MOMENTUM1, VECTOR, 
      NO_ROT, {1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  {  "FLOW_GRADV_SIC",  "FLOW_GRADV_SIC_BC",  STRONG_INT_SURF, FLOW_GRADV_SIC_BC,  R_MOMENTUM1, VECTOR, 
      NO_ROT, {1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  {  "HYDROSTATIC_SYMM", "HYDROSTATIC_SYMM_BC", WEAK_INT_SURF, HYDROSTATIC_SYMM_BC, R_MOMENTUM1, VECTOR , 
      NO_ROT, {0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,

  { "SHEET_ENDSLOPE", "IDLER_LOC", SPECIAL, SHEET_ENDSLOPE_BC, R_MESH_NORMAL, SCALAR, R_MESH1, {0, 0, 0, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB },
  { "CA", "CA_BC", SPECIAL, CA_BC, R_MESH_NORMAL, SCALAR, R_MESH1, {0, 0, 0, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "CA_MOMENTUM", "CA_MOMENTUM_BC", SPECIAL, CA_MOMENTUM_BC, R_MOM_TANG1, SCALAR, R_MOMENTUM1, {1, 1, 1, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "CA_OR_FIX", "CA_OR_FIX_BC", SPECIAL, CA_OR_FIX_BC, R_MESH_NORMAL, SCALAR, R_MESH1, {0, 0, 0, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "MOVING_CA", "MOVING_CA_BC", SPECIAL, MOVING_CA_BC, R_MOM_NORMAL, SCALAR, R_MOMENTUM1, {1, 1, 1, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "VELO_THETA_TPL", 		/* name1 */
    "VELO_THETA_TPL_BC", 	/* name2 */
    SPECIAL, 			/* method from rf_bc_const.h */
    VELO_THETA_TPL_BC, 		/* BC_Name from rf_bc_const.h */
    R_MOM_TANG1, 		/* equation from rf_fem_const.h */
    SCALAR, 			/* vector (either VECTOR or SCALAR) */
    R_MOMENTUM1, 		/* rotate (meaningful if vector) */
    {1, 1, 1, 0, 0, 1, 1, 1, 0, 0}, /* sensitivity booleans (velocity, mesh) */
    SINGLE_PHASE, 		/* evaluate only from one phase */
    DVI_SINGLE_PHASE_DB 	/* special flags for discontinuous vbl (no) */
  },
  { "VELO_THETA_HOFFMAN", 		/* name1 */
    "VELO_THETA_HOFFMAN_BC", 	/* name2 */
    SPECIAL, 			/* method from rf_bc_const.h */
    VELO_THETA_HOFFMAN_BC, 		/* BC_Name from rf_bc_const.h */
    R_MOM_TANG1, 		/* equation from rf_fem_const.h */
    SCALAR, 			/* vector (either VECTOR or SCALAR) */
    R_MOMENTUM1, 		/* rotate (meaningful if vector) */
    {1, 1, 1, 0, 0, 1, 1, 1, 0, 0}, /* sensitivity booleans (velocity, mesh) */
    SINGLE_PHASE, 		/* evaluate only from one phase */
    DVI_SINGLE_PHASE_DB 	/* special flags for discontinuous vbl (no) */
  },
  { "VELO_THETA_COX", 		/* name1 */
    "VELO_THETA_COX_BC", 	/* name2 */
    SPECIAL, 			/* method from rf_bc_const.h */
    VELO_THETA_COX_BC, 		/* BC_Name from rf_bc_const.h */
    R_MOM_TANG1, 		/* equation from rf_fem_const.h */
    SCALAR, 			/* vector (either VECTOR or SCALAR) */
    R_MOMENTUM1, 		/* rotate (meaningful if vector) */
    {1, 1, 1, 0, 0, 1, 1, 1, 0, 0}, /* sensitivity booleans (velocity, mesh) */
    SINGLE_PHASE, 		/* evaluate only from one phase */
    DVI_SINGLE_PHASE_DB 	/* special flags for discontinuous vbl (no) */
  },
  { "VELO_THETA_SHIK", 		/* name1 */
    "VELO_THETA_SHIK_BC", 	/* name2 */
    SPECIAL, 			/* method from rf_bc_const.h */
    VELO_THETA_SHIK_BC, 		/* BC_Name from rf_bc_const.h */
    R_MOM_TANG1, 		/* equation from rf_fem_const.h */
    SCALAR, 			/* vector (either VECTOR or SCALAR) */
    R_MOMENTUM1, 		/* rotate (meaningful if vector) */
    {1, 1, 1, 0, 0, 1, 1, 1, 0, 0}, /* sensitivity booleans (velocity, mesh) */
    SINGLE_PHASE, 		/* evaluate only from one phase */
    DVI_SINGLE_PHASE_DB 	/* special flags for discontinuous vbl (no) */
  },

  { "SURFTANG",        "SURFTANG_BC",        SPECIAL, SURFTANG_BC,        R_MOMENTUM1, VECTOR, NO_ROT, {0, 0, 0, 1, 1, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "SURFTANG_SCALAR", "SURFTANG_SCALAR_BC", SPECIAL, SURFTANG_SCALAR_BC, R_MOMENTUM1, VECTOR, NO_ROT, {0, 0, 0, 1, 1, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,

  { "CAP_ENDFORCE",  "CAP_ENDFORCE_BC", SPECIAL, CAP_ENDFORCE_BC, R_MOMENTUM1, VECTOR, NO_ROT, {0, 0, 0, 1, 1, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "CAP_ENDFORCE_SCALAR", "CAP_ENDFORCE_SCALAR_BC", SPECIAL, CAP_ENDFORCE_SCALAR_BC, R_MOMENTUM1, VECTOR, NO_ROT, {0, 0, 0, 1, 1, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,

  { "CA_EDGE", "CA_EDGE_BC", COLLOCATE_EDGE, CA_EDGE_BC, R_MESH_NORMAL, SCALAR, R_MESH1, {0, 0, 0, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "CA_EDGE_OR_FIX", "CA_EDGE_OR_FIX_BC", COLLOCATE_EDGE, CA_EDGE_OR_FIX_BC, R_MESH_NORMAL, SCALAR, R_MESH1, {0, 0, 0, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "CA_EDGE_INT", "CA_EDGE_INT_BC", STRONG_INT_EDGE, CA_EDGE_INT_BC, R_MESH_NORMAL, SCALAR, R_MESH1, {0, 0, 0, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "VAR_CA_EDGE", "VAR_CA_EDGE_BC", STRONG_INT_EDGE, VAR_CA_EDGE_BC, R_MESH_NORMAL, SCALAR, R_MESH1, {0, 0, 0, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "VAR_CA_USER", "VAR_CA_USER_BC", STRONG_INT_EDGE, VAR_CA_USER_BC, R_MESH_NORMAL, SCALAR, R_MESH1, {0, 0, 0, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,  
  { "CA_EDGE_CURVE", "CA_EDGE_CURVE_BC", COLLOCATE_EDGE, CA_EDGE_CURVE_BC, R_MESH_NORMAL, SCALAR, R_MESH1, {0, 0, 0, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "CA_EDGE_CURVE_INT", "CA_EDGE_CURVE_INT_BC", STRONG_INT_EDGE, CA_EDGE_CURVE_INT_BC, R_MESH_NORMAL, SCALAR, R_MESH1, {0, 0, 0, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,  
  { "SURFTANG_EDGE", "SURFTANG_EDGE_BC", WEAK_INT_EDGE, SURFTANG_EDGE_BC, R_MOMENTUM1, VECTOR, NO_ROT, {0, 0, 0, 1, 1, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "SURFTANG_SCALAR_EDGE", "SURFTANG_SCALAR_EDGE_BC", WEAK_INT_EDGE, SURFTANG_SCALAR_EDGE_BC, R_MOMENTUM1, VECTOR, NO_ROT, {0, 0, 0, 1, 1, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,

  { "FIX", "FIX_BC", DIRICHLET, FIX_BC, R_ANYTHING, SCALAR, NO_ROT, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "PERIODIC", "PERIODIC_BC", DIRICHLET, PERIODIC_BC, R_ANYTHING, SCALAR, NO_ROT, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,

  { "GD_CONST",  "GD_CONST_BC",  COLLOCATE_SURF, GD_CONST_BC,  R_ANYTHING, SCALAR, NO_ROT, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "GD_LINEAR", "GD_LINEAR_BC", COLLOCATE_SURF, GD_LINEAR_BC, R_ANYTHING, SCALAR, NO_ROT, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "GD_INVERSE", "GD_INVERSE_BC", COLLOCATE_SURF, GD_INVERSE_BC, R_ANYTHING, SCALAR, NO_ROT, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "GD_PARAB",  "GD_PARAB_BC",  COLLOCATE_SURF, GD_PARAB_BC,  R_ANYTHING, SCALAR, NO_ROT, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "GD_PARAB_OFFSET",  "GD_PARAB_OFFSET_BC",  COLLOCATE_SURF, GD_PARAB_OFFSET_BC,  R_ANYTHING, SCALAR, NO_ROT, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "GD_CIRC",   "GD_CIRC_BC",   COLLOCATE_SURF, GD_CIRC_BC,   R_ANYTHING, SCALAR, NO_ROT, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "GD_POLYN",  "GD_POLYN_BC",  COLLOCATE_SURF, GD_POLYN_BC,  R_ANYTHING, SCALAR, NO_ROT, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "GD_TABLE",  "GD_TABLE_BC",  COLLOCATE_SURF, GD_TABLE_BC,  R_ANYTHING, SCALAR, NO_ROT, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "GD_TIME",  "GD_TIME_BC",  COLLOCATE_SURF, GD_TIME_BC,  R_ANYTHING, SCALAR, NO_ROT, {1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,

  { "FORCE",      "FORCE_BC",      WEAK_INT_SURF, FORCE_BC,      R_MESH1, VECTOR, NO_ROT, {0, 0, 0, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "FORCE_SIC",      "FORCE_SIC_BC",      STRONG_INT_SURF, FORCE_SIC_BC,      R_MESH1, VECTOR, NO_ROT, {0, 0, 0, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "FORCE_RS",      "FORCE_RS_BC",      WEAK_INT_SURF, FORCE_RS_BC,      R_SOLID1, VECTOR, NO_ROT, {0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,1,1,1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "FORCE_USER",      "FORCE_USER_BC",      WEAK_INT_SURF, FORCE_USER_BC,      R_MESH1, VECTOR, NO_ROT, {0, 0, 0, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "FORCE_USER_SIC",      "FORCE_USER_SIC_BC",      STRONG_INT_SURF, FORCE_USER_SIC_BC,      R_MESH1, VECTOR, NO_ROT, {0, 0, 0, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "FORCE_USER_RS",      "FORCE_USER_RS_BC",      WEAK_INT_SURF, FORCE_USER_RS_BC,      R_SOLID1, VECTOR, NO_ROT, {0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,1,1,1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "REP_FORCE",      "REP_FORCE_BC",      WEAK_INT_SURF, REP_FORCE_BC,      R_MESH1, VECTOR, NO_ROT, {0, 0, 0, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "REP_FORCE_RS",      "REP_FORCE_RS_BC",      WEAK_INT_SURF, REP_FORCE_RS_BC,      R_SOLID1, VECTOR, NO_ROT,  {0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,1,1,1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "ATTR_FORCE",      "ATTR_FORCE_BC",      WEAK_INT_SURF, ATTR_FORCE_BC,      R_MESH1, VECTOR, NO_ROT, {0, 0, 0, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "ATTR_FORCE_RS",      "ATTR_FORCE_RS_BC",      WEAK_INT_SURF, ATTR_FORCE_RS_BC,      R_SOLID1, VECTOR, NO_ROT,  {0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,1,1,1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "REP_FORCE_ROLL",      "REP_FORCE_ROLL_BC",      WEAK_INT_SURF, REP_FORCE_ROLL_BC,      R_MESH1, VECTOR, NO_ROT, {0, 0, 0, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "REP_FORCE_ROLL_RS",      "REP_FORCE_ROLL_RS_BC",      WEAK_INT_SURF, REP_FORCE_ROLL_RS_BC,      R_SOLID1, VECTOR, NO_ROT,  {0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,1,1,1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "REP_FORCE_SHU", "REP_FORCE_SHU_BC", WEAK_INT_SURF, REP_FORCE_SHU_BC, R_MESH1, VECTOR, NO_ROT, {0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB},
  { "REP_FORCE_SHU_SIC", "REP_FORCE_SHU_SIC_BC", STRONG_INT_SURF, REP_FORCE_SHU_SIC_BC, R_MESH1, VECTOR, NO_ROT, {0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB},
  { "NORM_FORCE", "NORM_FORCE_BC", WEAK_INT_SURF, NORM_FORCE_BC, R_MESH1, VECTOR, NO_ROT, {0, 0, 0, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "NORM_FORCE_RS", "NORM_FORCE_RS_BC", WEAK_INT_SURF, NORM_FORCE_RS_BC, R_SOLID1, VECTOR, NO_ROT, {0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,1,1,1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "FRICTION",      "FRICTION_BC",      STRONG_INT_SURF, FRICTION_BC,      R_MESH_TANG1, SCALAR, R_MESH1, {0,0,0,0,0,1,1,1,0,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "FRICTION_ACOUSTIC",      "FRICTION_ACOUSTIC_BC",      STRONG_INT_SURF, FRICTION_ACOUSTIC_BC,      R_MESH_TANG1, SCALAR, R_MESH1, {0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "FRICTION_RS",      "FRICTION_RS_BC",      STRONG_INT_SURF, FRICTION_RS_BC,      R_SOLID_TANG1, SCALAR, R_SOLID1,  {0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,1,1,1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "FRICTION_ACOUSTIC_RS",      "FRICTION_ACOUSTIC_RS_BC",      STRONG_INT_SURF, FRICTION_ACOUSTIC_RS_BC,      R_SOLID_TANG1, SCALAR, R_SOLID1, {0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "SLOPE", "SLOPE_BC", STRONG_INT_SURF, SLOPE_BC, R_MESH_NORMAL, SCALAR, R_MESH1, {0, 0, 0, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "SLOPEX", "SLOPEX_BC", STRONG_INT_SURF, SLOPEX_BC, R_MESH1, SCALAR, NO_ROT, {0, 0, 0, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "SLOPEY", "SLOPEY_BC", STRONG_INT_SURF, SLOPEY_BC, R_MESH2, SCALAR, NO_ROT, {0, 0, 0, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "SLOPEZ", "SLOPEZ_BC", STRONG_INT_SURF, SLOPEZ_BC, R_MESH3, SCALAR, NO_ROT, {0, 0, 0, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,

  { "FLUID_SOLID", "FLUID_SOLID_BC", COLLOCATE_SURF, FLUID_SOLID_BC, R_MOMENTUM1, VECTOR, NO_ROT, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, CROSS_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "SOLID_FLUID", "SOLID_FLUID_BC", COLLOCATE_SURF, SOLID_FLUID_BC, R_MESH1,     VECTOR, NO_ROT, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, CROSS_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "FLUID_SOLID_RS", "FLUID_SOLID_RS_BC", WEAK_SHIFT, FLUID_SOLID_RS_BC, R_MOMENTUM1, VECTOR, NO_ROT, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, CROSS_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "SOLID_FLUID_RS", "SOLID_FLUID_RS_BC", COLLOCATE_SURF, SOLID_FLUID_RS_BC, R_SOLID1,     VECTOR, NO_ROT, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, CROSS_PHASE, DVI_CROSS_PHASE_CONJUGATE } ,
  { "SOLID_FLUID_CONTACT", "SOLID_FLUID_CONTACT_BC", CONTACT_SURF, SOLID_FLUID_CONTACT_BC, R_MESH1,     VECTOR, NO_ROT, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, SINGLE_PHASE, DVI_CROSS_PHASE_CONJUGATE } ,
  { "SOLID_LAGRANGE_MULT", "SOLID_LAGRANGE_MULT", CONTACT_SURF, SOLID_LAGRANGE_MULT_BC, R_LAGR_MULT1,     VECTOR, NO_ROT, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, SINGLE_PHASE, DVI_CROSS_PHASE_CONJUGATE } ,
  { "LAGRANGE_NO_SLIP", "LAGRANGE_NO_SLIP_BC", CONTACT_SURF, LAGRANGE_NO_SLIP_BC, R_LAGR_MULT1,     VECTOR, NO_ROT, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, SINGLE_PHASE, DVI_CROSS_PHASE_CONJUGATE } ,
  { "BAAIJENS_SOLID_FLUID", "BAAIJENS_SOLID_FLUID_BC", CONTACT_SURF, BAAIJENS_SOLID_FLUID_BC, R_MESH1,     VECTOR, NO_ROT, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, SINGLE_PHASE, DVI_CROSS_PHASE_CONJUGATE } ,
  { "BAAIJENS_FLUID_SOLID", "BAAIJENS_FLUID_SOLID_BC", EMBEDDED_SURF, BAAIJENS_FLUID_SOLID_BC, R_MOMENTUM1,     VECTOR, NO_ROT, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, SINGLE_PHASE, DVI_CROSS_PHASE_CONJUGATE } ,
  { "FLUID_SOLID_CONTACT", "FLUID_SOLID_CONTACT_BC", EMBEDDED_SURF, FLUID_SOLID_CONTACT_BC, R_MOMENTUM1,     VECTOR, NO_ROT, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, SINGLE_PHASE, DVI_CROSS_PHASE_CONJUGATE } ,
  { "KIN_DISPLACEMENT", "KIN_DISPLACEMENT_BC", STRONG_INT_SURF, KIN_DISPLACEMENT_BC, R_MESH_NORMAL, SCALAR, R_MESH1,  {0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,1,1,1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
{ "KIN_DISPLACEMENT_RS", "KIN_DISPLACEMENT_RS_BC", STRONG_INT_SURF, KIN_DISPLACEMENT_RS_BC, R_MESH_NORMAL, SCALAR, R_MESH1,   {0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "KIN_DISPLACEMENT_PETROV", "KIN_DISPLACEMENT_PETROV_BC", STRONG_INT_SURF, KIN_DISPLACEMENT_PETROV_BC, R_MESH_NORMAL, SCALAR, R_MESH1,  {0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,1,1,1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "KIN_DISPLACEMENT_COLLOC", "KIN_DISPLACEMENT_COLLOC_BC", COLLOCATE_SURF, KIN_DISPLACEMENT_COLLOC_BC, R_MESH_NORMAL, SCALAR, R_MESH1,  {0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,1,1,1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "SHELL_SURFACE_CHARGE", "SHELL_SURFACE_CHARGE_BC", WEAK_INT_SURF, SHELL_SURFACE_CHARGE_BC, R_POTENTIAL, SCALAR, NO_ROT, {0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB},
  { "SHELL_SURFACE_CHARGE_SIC", "SHELL_SURFACE_CHARGE_SIC_BC", STRONG_INT_SURF, SHELL_SURFACE_CHARGE_SIC_BC, R_POTENTIAL, SCALAR, NO_ROT, {0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB},
  { "SURFACE_ELECTRIC_FIELD", "SURFACE_ELECTRIC_FIELD_BC", WEAK_SHELL_GRAD, SURFACE_ELECTRIC_FIELD_BC, R_POTENTIAL, SCALAR, NO_ROT, {0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB},
  { "SURFACE_ACOUSTIC_VELOCITY", "SURFACE_ACOUSTIC_VELOCITY_BC", WEAK_SHELL_GRAD, SURFACE_ACOUSTIC_VELOCITY_BC, R_ACOUS_PREAL, SCALAR, NO_ROT, {0,0,0,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB},
 { "SHELL_DIFF_KINEMATIC", "SHELL_DIFF_KINEMATIC_BC", STRONG_INT_SURF, SHELL_DIFF_KINEMATIC_BC, R_MESH_NORMAL, SCALAR, R_MESH1, {0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB},
  { "SURFACE_USER_SHELL", "SURFACE_USER_SHELL_BC", WEAK_SHELL_GRAD, SURFACE_USER_SHELL_BC, R_SHELL_USER, SCALAR, NO_ROT, {0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB},
  { "SURFACE_LUBRICATION", "SURFACE_LUBRICATION_BC", WEAK_SHELL_GRAD, SURFACE_LUBRICATION_BC, R_SHELL_LUBP, SCALAR, NO_ROT, {0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB},
  { "NO_SLIP", "NO_SLIP_BC", STRONG_INT_SURF, NO_SLIP_BC, R_MOMENTUM1, VECTOR, NO_ROT,   {1, 1, 1, 1, 1, 1, 1, 1, 0, 1}, CROSS_PHASE, DVI_CROSS_PHASE_CONJUGATE } ,
  { "NO_SLIP_RS", "NO_SLIP_RS_BC", STRONG_INT_SURF, NO_SLIP_RS_BC, R_MOMENTUM1, VECTOR, NO_ROT,   {1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0,0,0,0,0,1,1,1}, CROSS_PHASE, DVI_CROSS_PHASE_CONJUGATE } ,
  { "P_EQUIL", "P_EQUIL_BC", STRONG_INT_SURF, P_EQUIL_BC, R_MASS, SCALAR, NO_ROT,    {0, 0, 0, 0, 1, 0, 0, 0, 0, 1}, CROSS_PHASE,DVI_CROSS_PHASE_CONJUGATE } ,
  { "VP_EQUIL",   "VP_EQUIL_BC", STRONG_INT_SURF, VP_EQUIL_BC,  R_ENERGY,     SCALAR, NO_ROT, {0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 1, 1, 1, 1 }, CROSS_PHASE, DVI_CROSS_PHASE_CONJUGATE } ,
  { "VN_POROUS", "VN_POROUS_BC", STRONG_INT_SURF, VN_POROUS_BC, R_MOM_NORMAL, SCALAR, R_MOMENTUM1, {1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 1, 1, 1, 1 }, CROSS_PHASE, DVI_CROSS_PHASE_CONJUGATE } ,
  { "POROUS_GAS", "POROUS_GAS_BC", STRONG_INT_SURF, POROUS_GAS_BC, R_POR_LIQ_PRES, SCALAR, NO_ROT, {1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 1, 1, 1, 1 }, CROSS_PHASE, DVI_CROSS_PHASE_CONJUGATE } ,
  { "CONT_TANG_VEL", "CONT_TANG_VEL_BC", STRONG_INT_SURF, CONT_TANG_VEL_BC, R_MOMENTUM2, SCALAR,  NO_ROT, {1, 1, 1, 0, 0, 1, 1, 1, 0, 0}, CROSS_PHASE_DISCONTINUOUS , DVI_SIDTIE} ,
  { "CONT_NORM_VEL", "CONT_NORM_VEL_BC", STRONG_INT_SURF, CONT_NORM_VEL_BC, R_MOMENTUM1, SCALAR,  NO_ROT, {1, 1, 1, 0, 0, 1, 1, 1, 0, 0}, CROSS_PHASE_DISCONTINUOUS, DVI_SIDTIE } ,
  { "DISCONTINUOUS_VELO", "DISCONTINUOUS_VELO_BC", STRONG_INT_SURF, DISCONTINUOUS_VELO_BC, R_MOMENTUM1, SCALAR,  NO_ROT, {1, 1, 1, 1, 1, 1, 1, 1, 0, 0}, CROSS_PHASE_DISCONTINUOUS, DVI_SIDTIE } ,
  { "VL_EQUIL",   "VL_EQUIL_BC", STRONG_INT_SURF, VL_EQUIL_BC,  R_MASS,       SCALAR, NO_ROT, {0, 0, 0, 1, 1, 0, 0, 0, 0, 1}, CROSS_PHASE_DISCONTINUOUS, DVI_SIDTIE } ,
  { "Y_DISCONTINUOUS",   "Y_DISCONTINUOUS_BC", DIRICHLET, Y_DISCONTINUOUS_BC,  R_MASS,       SCALAR, NO_ROT, {0, 0, 0, 0, 1, 0, 0, 0, 0, 0}, CROSS_PHASE_DISCONTINUOUS, DVI_SIDTIE } ,
  { "VL_POLY",   "VL_POLY_BC", STRONG_INT_SURF, VL_POLY_BC,  R_MASS, SCALAR, NO_ROT, {0, 0, 0, 1, 1, 0, 0, 0, 0, 1}, CROSS_PHASE_DISCONTINUOUS, DVI_SIDTIE } ,
  { "LATENT_HEAT", "LATENT_HEAT_BC", WEAK_INT_SURF, LATENT_HEAT_BC, R_ENERGY, SCALAR, NO_ROT, {0, 0, 0, 1, 1, 0, 0, 0, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "T_CONTACT_RESIS",   "T_CONTACT_RESIS_BC", WEAK_INT_SURF, T_CONTACT_RESIS_BC,  R_ENERGY,       SCALAR, NO_ROT, {0, 0, 0, 1, 1, 0, 0, 0, 0, 0}, CROSS_PHASE_DISCONTINUOUS,  DVI_DVVSIG } ,
  { "T_CONTACT_RESIS_2",   "T_CONTACT_RESIS_2_BC", WEAK_INT_SURF, T_CONTACT_RESIS_2_BC,  R_ENERGY,       SCALAR, NO_ROT, {0, 0, 0, 1, 1, 0, 0, 0, 0, 0}, CROSS_PHASE_DISCONTINUOUS,  DVI_DVVSIG } ,


  { "PU", "PU_BC", DIRICHLET, PU_BC, R_PMOMENTUM1, SCALAR, NO_ROT, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "PV", "PV_BC", DIRICHLET, PV_BC, R_PMOMENTUM2, SCALAR, NO_ROT, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "PW", "PW_BC", DIRICHLET, PW_BC, R_PMOMENTUM3, SCALAR, NO_ROT, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "LATENT_HEAT_INTERNAL", "LATENT_HEAT_INTERNAL_BC", WEAK_INT_SURF, LATENT_HEAT_INTERNAL_BC, R_ENERGY, SCALAR, NO_ROT, {1, 1, 1, 1, 1, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "TABLE", "TABLE_BC", COLLOCATE_SURF, TABLE_BC, R_ANYTHING, SCALAR, NO_ROT, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "TABLE_WICV", "TABLE_WICV_BC", WEAK_INT_SURF, TABLE_WICV_BC, R_ANYTHING, VECTOR , NO_ROT, {0, 0, 0, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "TABLE_WICS", "TABLE_WICS_BC", WEAK_INT_SURF, TABLE_WICS_BC, R_ANYTHING, SCALAR , NO_ROT, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  /* Lagrange multiplier boundary condtions */
  { "FLOWRATE",    "FLOWRATE_BC",    WEAK_INT_SURF, LGR_FLOWRATE_BC, R_MOMENTUM1, VECTOR, NO_ROT, {0, 0, 0, 0, 0, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  /* level set "boundary" conditions */
  { "LS_CAPILLARY", "LS_CAPILLARY_BC", WEAK_INT_SURF, LS_CAPILLARY_BC, R_MOMENTUM1, VECTOR, NO_ROT,   {1, 1, 1, 1, 1, 1, 1, 1, 0, 1}, CROSS_PHASE, DVI_CROSS_PHASE_CONJUGATE } ,
  { "LS_CAP_CURVE", "LS_CAP_CURVE_BC", WEAK_INT_SURF, LS_CAP_CURVE_BC, R_MOMENTUM1, VECTOR, NO_ROT,   {1, 1, 1, 1, 1, 1, 1, 1, 0, 1}, CROSS_PHASE, DVI_CROSS_PHASE_CONJUGATE } ,
  { "LS_CAP_DIV_N", "LS_CAP_DIV_N_BC", WEAK_INT_SURF, LS_CAP_DIV_N_BC, R_MOMENTUM1, VECTOR, NO_ROT,   {1, 1, 1, 1, 1, 1, 1, 1, 0, 1}, CROSS_PHASE, DVI_CROSS_PHASE_CONJUGATE } ,
  { "LS_CAP_DIV_S_N", "LS_CAP_DIV_S_N_BC", WEAK_INT_SURF, LS_CAP_DIV_S_N_BC, R_MOMENTUM1, VECTOR, NO_ROT,   {1, 1, 1, 1, 1, 1, 1, 1, 0, 1}, CROSS_PHASE, DVI_CROSS_PHASE_CONJUGATE } ,
  { "LS_NO_SLIP", "LS_NO_SLIP_BC", EMBEDDED_SURF, LS_NO_SLIP_BC, R_LAGR_MULT1, VECTOR, NO_ROT,   {1, 1, 1, 1, 1, 1, 1, 1, 0, 1}, SINGLE_PHASE, DVI_CROSS_PHASE_CONJUGATE } ,
  { "LS_CAPILLARY_GHOST", "LS_CAPILLARY_GHOST_BC", EMBEDDED_SURF, LS_CAPILLARY_GHOST_BC, R_PRESSURE, SCALAR, NO_ROT,   {1, 1, 1, 1, 1, 1, 1, 1, 0, 1}, SINGLE_PHASE, DVI_CROSS_PHASE_CONJUGATE } ,
  { "LS_CAP_HYSING", "LS_CAP_HYSING_BC", WEAK_INT_SURF, LS_CAP_HYSING_BC, R_MOMENTUM1, VECTOR, NO_ROT,   {1, 1, 1, 1, 1, 1, 1, 1, 0, 1}, CROSS_PHASE, DVI_CROSS_PHASE_CONJUGATE } ,
  { "LS_CAP_DENNER_DIFF", "LS_CAP_DENNER_DIFF_BC", WEAK_INT_SURF, LS_CAP_DENNER_DIFF_BC, R_MOMENTUM1, VECTOR, NO_ROT,   {1, 1, 1, 1, 1, 1, 1, 1, 0, 1}, CROSS_PHASE, DVI_CROSS_PHASE_CONJUGATE } ,
  { "LS_T", "LS_T_BC", EMBEDDED_SURF, LS_T_BC, R_ENERGY, SCALAR, NO_ROT,   {1, 1, 1, 1, 1, 1, 1, 1, 0, 1}, CROSS_PHASE, DVI_CROSS_PHASE_CONJUGATE } ,
  { "LS_U", "LS_U_BC", EMBEDDED_SURF, LS_U_BC, R_MOMENTUM1, SCALAR, NO_ROT,   {1, 1, 1, 1, 1, 1, 1, 1, 0, 1}, CROSS_PHASE, DVI_CROSS_PHASE_CONJUGATE } ,
  { "LS_V", "LS_V_BC", EMBEDDED_SURF, LS_V_BC, R_MOMENTUM2, SCALAR, NO_ROT,   {1, 1, 1, 1, 1, 1, 1, 1, 0, 1}, CROSS_PHASE, DVI_CROSS_PHASE_CONJUGATE } ,
  { "LS_W", "LS_W_BC", EMBEDDED_SURF, LS_W_BC, R_MOMENTUM3, SCALAR, NO_ROT,   {1, 1, 1, 1, 1, 1, 1, 1, 0, 1}, CROSS_PHASE, DVI_CROSS_PHASE_CONJUGATE } ,
  { "LS_CONT_T", "LS_CONT_T_BC", EMBEDDED_SURF, LS_CONT_T_BC, R_ENERGY, SCALAR, NO_ROT,   {1, 1, 1, 1, 1, 1, 1, 1, 0, 1}, CROSS_PHASE, DVI_CROSS_PHASE_CONJUGATE } ,
  { "LS_CONT_FLUX", "LS_CONT_FLUX_BC", EMBEDDED_SURF, LS_CONT_FLUX_BC, R_ENERGY, SCALAR, NO_ROT,   {1, 1, 1, 1, 1, 1, 1, 1, 0, 1}, CROSS_PHASE, DVI_CROSS_PHASE_CONJUGATE } ,
  { "LS_CONT_VEL", "LS_CONT_VEL_BC", EMBEDDED_SURF, LS_CONT_VEL_BC, R_MOMENTUM1, VECTOR, NO_ROT,   {1, 1, 1, 1, 1, 1, 1, 1, 0, 1}, CROSS_PHASE, DVI_CROSS_PHASE_CONJUGATE } ,
  { "LS_CONT_TRACTION", "LS_CONT_TRACTION_BC", EMBEDDED_SURF, LS_CONT_TRACTION_BC, R_MOMENTUM1, VECTOR, NO_ROT,   {1, 1, 1, 1, 1, 1, 1, 1, 0, 1}, CROSS_PHASE, DVI_CROSS_PHASE_CONJUGATE } ,
  { "LS_FLOW_PRESSURE", "LS_FLOW_PRESSURE_BC", EMBEDDED_SURF, LS_FLOW_PRESSURE_BC, R_MOMENTUM1, VECTOR, NO_ROT,   {1, 1, 1, 1, 1, 1, 1, 1, 0, 1}, CROSS_PHASE, DVI_CROSS_PHASE_CONJUGATE } ,
  { "LS_ACOUSTIC_SOURCE", "LS_ACOUSTIC_SOURCE_BC", EMBEDDED_SURF, LS_ACOUSTIC_SOURCE_BC, R_MOMENTUM1, VECTOR, NO_ROT,   {1,1,1,1,1,1,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1} ,CROSS_PHASE, DVI_CROSS_PHASE_CONJUGATE } ,
  { "LS_Q", "LS_Q_BC", EMBEDDED_SURF, LS_Q_BC, R_ENERGY, SCALAR, NO_ROT,   {1, 1, 1, 1, 1, 1, 1, 1, 0, 1}, CROSS_PHASE, DVI_CROSS_PHASE_CONJUGATE } ,
  { "LS_QRAD", "LS_QRAD_BC", EMBEDDED_SURF, LS_QRAD_BC, R_ENERGY, SCALAR, NO_ROT,   {1, 1, 1, 1, 1, 1, 1, 1, 0, 1}, CROSS_PHASE, DVI_CROSS_PHASE_CONJUGATE } ,
 { "LS_RECOIL_PRESSURE", "LS_RECOIL_PRESSURE_BC", EMBEDDED_SURF, LS_RECOIL_PRESSURE_BC, R_MOMENTUM1, VECTOR, NO_ROT,   {1, 1, 1, 1, 1, 1, 1, 1, 0, 1}, CROSS_PHASE, DVI_CROSS_PHASE_CONJUGATE } ,
  { "LS_Q", "LS_Q_BC", EMBEDDED_SURF, LS_Q_BC, R_ENERGY, SCALAR, NO_ROT,   {1, 1, 1, 1, 1, 1, 1, 1, 0, 1}, CROSS_PHASE, DVI_CROSS_PHASE_CONJUGATE } ,
  { "LS_QLASER", "LS_QLASER_BC", EMBEDDED_SURF, LS_QLASER_BC, R_ENERGY, SCALAR, NO_ROT,   {1, 1, 1, 1, 1, 1, 1, 1, 0, 1}, CROSS_PHASE, DVI_CROSS_PHASE_CONJUGATE } ,
  { "LS_VAPOR", "LS_QVAPOR_BC", EMBEDDED_SURF, LS_QVAPOR_BC, R_ENERGY, SCALAR, NO_ROT,   {1, 1, 1, 1, 1, 1, 1, 1, 0, 1}, CROSS_PHASE, DVI_CROSS_PHASE_CONJUGATE } ,
  { "LS_SOLID_FLUID", "LS_SOLID_FLUID_BC", STRONG_INT_SURF, LS_SOLID_FLUID_BC, R_SOLID1, VECTOR, NO_ROT,   {1, 1, 1, 1, 1, 1, 1, 1, 0, 1}, CROSS_PHASE, DVI_CROSS_PHASE_CONJUGATE } ,
  { "LS_EXTV_FLUID_SIC", "LS_EXTV_FLUID_SIC_BC", EMBEDDED_SURF, LS_EXTV_FLUID_SIC_BC, R_EXT_VELOCITY, SCALAR, NO_ROT,   {1, 1, 1, 1, 1, 1, 1, 1, 0, 1}, CROSS_PHASE, DVI_CROSS_PHASE_CONJUGATE } ,
  { "LS_EIK_KINEMATIC", "LS_EIK_KINEMATIC_BC", EMBEDDED_SURF, LS_EIK_KINEMATIC_BC, R_FILL, SCALAR, NO_ROT,   {1, 1, 1, 1, 1, 1, 1, 1, 0, 1}, CROSS_PHASE, DVI_CROSS_PHASE_CONJUGATE } ,
  { "LS_EXTV_KINEMATIC", "LS_EXTV_KINEMATIC_BC", EMBEDDED_SURF, LS_EXTV_KINEMATIC_BC, R_FILL, SCALAR, NO_ROT,   {1, 1, 1, 1, 1, 1, 1, 1, 0, 1}, CROSS_PHASE, DVI_CROSS_PHASE_CONJUGATE } ,
  { "LS_EIK_KIN_LEAK", "LS_EIK_KIN_LEAK_BC", EMBEDDED_SURF, LS_EIK_KIN_LEAK_BC, R_FILL, SCALAR, NO_ROT,   {1, 1, 1, 1, 1, 1, 1, 1, 0, 1}, CROSS_PHASE, DVI_CROSS_PHASE_CONJUGATE } ,
  { "LS_EXTV_KIN_LEAK", "LS_EXTV_KIN_LEAK_BC", EMBEDDED_SURF, LS_EXTV_KIN_LEAK_BC, R_FILL, SCALAR, NO_ROT,   {1, 1, 1, 1, 1, 1, 1, 1, 0, 1}, CROSS_PHASE, DVI_CROSS_PHASE_CONJUGATE } ,
  { "LS_EXTV_LATENT", "LS_EXTV_LATENT_BC", EMBEDDED_SURF, LS_EXTV_LATENT_BC, R_FILL, SCALAR, NO_ROT,   {1, 1, 1, 1, 1, 1, 1, 1, 0, 1}, CROSS_PHASE, DVI_CROSS_PHASE_CONJUGATE } ,
  { "LS_YFLUX", "LS_YFLUX_BC", EMBEDDED_SURF, LS_YFLUX_BC, R_MASS, SCALAR, NO_ROT, {1, 1, 1, 1, 1, 1, 1, 1, 0, 1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "LS_LATENT_HEAT", "LS_LATENT_HEAT_BC", EMBEDDED_SURF, LS_LATENT_HEAT_BC, R_ENERGY, SCALAR, NO_ROT, {1, 1, 1, 1, 1, 1, 1, 1, 0, 1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "LS_ADC_MESH",    "LS_ADC_MESH_BC", LS_SPECIAL, LS_ADC_OLD_BC,  R_LEVEL_SET,  SCALAR, NO_ROT, {0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "LS_ADC",    "LS_ADC_BC", LS_SPECIAL, LS_ADC_BC,  R_LEVEL_SET,  SCALAR, NO_ROT, {0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "PF_CAPILLARY", "PF_CAPILLARY_BC", WEAK_INT_SURF, PF_CAPILLARY_BC, R_MOMENTUM1, VECTOR, NO_ROT,   {1, 1, 1, 1, 1, 1, 1, 1, 0, 1}, CROSS_PHASE, DVI_CROSS_PHASE_CONJUGATE } ,

  { "YFLUX_DISC_RXN", "YFLUX_DISC_RXN_BC", WEAK_INT_SURF, YFLUX_DISC_RXN_BC, R_MASS, SCALAR, NO_ROT, {1,1,1,1,1,1,1,1,0,1}, CROSS_PHASE_DISCONTINUOUS, DVI_MULTI_PHASE_VD } ,

  /* HKM -> Chemkin Boundary conditions , need ones on momentum, heat, and mass transfer */

  { "VL_EQUIL_PSEUDORXN", "VL_EQUIL_PRXN_BC", WEAK_INT_SURF, VL_EQUIL_PRXN_BC, R_MASS, SCALAR, NO_ROT, {1,1,1,1,1,1,1,1,0,1}, CROSS_PHASE_DISCONTINUOUS, DVI_MULTI_PHASE_VD } ,
  { "IS_EQUIL_PSEUDORXN", "IS_EQUIL_PRXN_BC", WEAK_INT_SURF, IS_EQUIL_PRXN_BC, R_MASS, SCALAR, NO_ROT, {1,1,1,1,1,1,1,1,0,1}, CROSS_PHASE_DISCONTINUOUS, DVI_MULTI_PHASE_VD } ,
  { "SPECIESFLUX_CHEM", "SF_CHEM_BC", WEAK_INT_SURF, SF_CHEM_BC, R_MASS, SCALAR, NO_ROT, {1,1,1,1,1,1,1,1,0,1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "KINEMATIC_SPECIESCHEM", "KINEMATIC_SC_BC", STRONG_INT_SURF, KINEMATIC_SC_BC, R_MESH_NORMAL, SCALAR, R_MESH1, {0, 0, 0, 1, 1, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB } ,
  { "MASSFLUX_STEFAN_FLOW", "MF_STEFANFLOW_BC", STRONG_INT_SURF, MF_STEFANFLOW_BC, R_MOM_NORMAL, SCALAR, R_MOMENTUM1, {1,1,1,1,1,1,1,1,0,1}, CROSS_PHASE_DISCONTINUOUS, DVI_SIDTIE } ,
  { "VNORM_STEFAN_FLOW", "VN_STEFANFLOW_BC", STRONG_INT_SURF, VN_STEFANFLOW_BC, R_MOM_NORMAL, SCALAR, R_MOMENTUM1, {1,1,1,0,1,1,1,1,0,0}, CROSS_PHASE_DISCONTINUOUS, DVI_SIDTIE } ,
  { "SURFDOMAINCHEMKIN_SURFACERXN", "SDC_SURFRXN_BC", WEAK_INT_SURF, SDC_SURFRXN_BC, R_MASS, SCALAR, NO_ROT, {1,1,1,1,1,1,1,1,1,1}, CROSS_PHASE_DISCONTINUOUS, DVI_MULTI_PHASE_VD } ,
  { "SURFDOMAINCHEMKIN_KIN_CHEM", "SDC_KIN_CHEM_BC", STRONG_INT_SURF, SDC_KIN_CHEM_BC, R_MESH_NORMAL, SCALAR, R_MESH1, {0, 0, 0, 1, 1, 1, 1, 1, 0, 0}, CROSS_PHASE_DISCONTINUOUS, DVI_MULTI_PHASE_SINGLE } ,
  { "SURFDOMAINCHEMKIN_KIN_STEFAN_FLOW", "SDC_KIN_SF_BC", STRONG_INT_SURF, SDC_KIN_SF_BC, R_MESH_NORMAL, SCALAR, R_MESH1, {0, 0, 0, 1, 1, 1, 1, 1, 0, 0}, SINGLE_PHASE, DVI_MULTI_PHASE_SINGLE } ,
  { "SURFDOMAINCHEMKIN_KIN_STEFAN_VOL_FLOW", "SDC_KIN_SFV_BC", STRONG_INT_SURF, SDC_KIN_SFV_BC, R_MESH_NORMAL, SCALAR, R_MESH1, {0, 0, 0, 1, 1, 1, 1, 1, 0, 0}, CROSS_PHASE_DISCONTINUOUS, DVI_MULTI_PHASE_SINGLE } ,
  { "SURFDOMAINCHEMKIN_HEATRXN", "SDC_HEATRXN_BC", WEAK_INT_SURF, SDC_HEATRXN_BC, R_ENERGY, SCALAR, NO_ROT, {1,1,1,1,1,1,1,1,1,1}, CROSS_PHASE_DISCONTINUOUS, DVI_VSIG } ,
  { "SURFDOMAINCHEMKIN_STEFAN_FLOW", "SDC_STEFANFLOW_BC", STRONG_INT_SURF, SDC_STEFANFLOW_BC, R_MOM_NORMAL, SCALAR, NO_ROT, {1,1,1,1,1,1,1,1,1,1}, CROSS_PHASE_DISCONTINUOUS, DVI_SIDTIE } ,
  { "SHELL_TFMP_PRES_BC", "SHELL_TFMP_PRES_BC", DIRICHLET, SHELL_TFMP_PRES_BC, R_TFMP_BOUND, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB},
  { "SHELL_TFMP_SAT_BC", "SHELL_TFMP_SAT_BC", DIRICHLET, SHELL_TFMP_SAT_BC, R_TFMP_MASS, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB},
  { "SHELL_TFMP_FREE_LIQ_BC", "SHELL_TFMP_FREE_LIQ_BC", WEAK_INT_SURF, SHELL_TFMP_FREE_LIQ_BC, R_TFMP_MASS, SCALAR, NO_ROT, {0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB},
  { "SHELL_TFMP_NUM_DIFF_BC", "SHELL_TFMP_NUM_DIFF_BC", WEAK_INT_SURF, SHELL_TFMP_NUM_DIFF_BC, R_TFMP_MASS, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB},
  { "SHELL_TFMP_GRAD_S_BC", "SHELL_TFMP_GRAD_S_BC", STRONG_INT_SURF, SHELL_TFMP_GRAD_S_BC, R_TFMP_MASS, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB},
  { "SHELL_TFMP_FREE_GAS_BC", "SHELL_TFMP_FREE_GAS_BC", WEAK_INT_SURF, SHELL_TFMP_FREE_GAS_BC, R_TFMP_BOUND, SCALAR, NO_ROT, {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB},
  { "SH_SDET_BC", "SH_SDET_BC", WEAK_INT_SURF, SH_SDET_BC, R_MESH1, SCALAR, NO_ROT, {0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB},
  { "SH_MESH2_WEAK_BC", "SH_MESH2_WEAK_BC", WEAK_INT_SURF, SH_MESH2_WEAK_BC, R_MESH2, SCALAR, NO_ROT, {0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB},
  { "SHELL_LUBRICATION_OUTFLOW_BC", "SHELL_LUBRICATION_OUTFLOW_BC", WEAK_INT_SURF, SHELL_LUBRICATION_OUTFLOW_BC, R_TFMP_BOUND, SCALAR, NO_ROT, {0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,}, SINGLE_PHASE, DVI_SINGLE_PHASE_DB},
};

/* order of sensitivities in list (as of 5/15/95): 
   vel1, vel2, vel3, temp, conc, mesh1, mesh2, mesh3, surf, press, stress tensor, velocity gradient */

int Num_BC_Names = sizeof(BC_Desc) / sizeof(struct BC_descriptions);


/*
 * EQ_Name []
 *  structure for IO processing of equations
 *  EQ_Name[] provides an association between an equation number
 *  and two strings (one long, one short) associated with that
 *  equation
 *    -> the numbers refer to equation numbers listed in
 *       rf_fem_const.h
 *
 *  Requirements:
 * ---------------
 *   The long and short strings must be unique.
 *   The numbers must be unique
 *   The first V_LAST entries must be contiguous and associated
 *     with equations listed in rf_fem_const.h
 *   After V_LAST, there must be a concentration section. This section
 *   must have more entries than there are species equations in the
 *   problem.
 */

struct Equation_Names EQ_Name[] = {
    { "R_MOMENTUM1", "VX", R_MOMENTUM1 } ,      /* 0 */
    { "R_MOMENTUM2", "VY", R_MOMENTUM2 } ,
    { "R_MOMENTUM3", "VZ", R_MOMENTUM3 } ,
    { "R_ENERGY",    "T",  R_ENERGY } ,
    { "R_MASS",      "Y",  R_MASS } ,
    { "R_MESH1",     "DMX", R_MESH1 } ,
    { "R_MESH2",     "DMY", R_MESH2 } ,
    { "R_MESH3",     "DMZ", R_MESH3 } ,
    { "R_MASS_SURF", "S",  R_MASS_SURF } ,
    { "R_PRESSURE",  "P",  R_PRESSURE } ,       /* 9 */

    { "R_STRESS11", "S11", R_STRESS11 } ,        
    { "R_STRESS12", "S12", R_STRESS12 } ,
    { "R_STRESS22", "S22", R_STRESS22 } ,
    { "R_STRESS13", "S13", R_STRESS13 } ,
    { "R_STRESS23", "S23", R_STRESS23 } ,
    { "R_STRESS33", "S33", R_STRESS33 } ,

    { "R_SOLID1",   "DMX_RS", R_SOLID1 } ,
    { "R_SOLID2",   "DMY_RS", R_SOLID2 } ,
    { "R_SOLID3",   "DMZ_RS", R_SOLID3 } ,

    { "R_GRADIENT11", "G11", R_GRADIENT11 } ,   /* 19 */
    { "R_GRADIENT12", "G12", R_GRADIENT12 } ,      
    { "R_GRADIENT21", "G21", R_GRADIENT21 } ,
    { "R_GRADIENT22", "G22", R_GRADIENT22 } ,
    { "R_GRADIENT13", "G13", R_GRADIENT13 } ,
    { "R_GRADIENT23", "G23", R_GRADIENT23 } ,
    { "R_GRADIENT31", "G31", R_GRADIENT31 } ,
    { "R_GRADIENT32", "G32", R_GRADIENT32 } ,
    { "R_GRADIENT33", "G33", R_GRADIENT33 } ,

    { "R_POTENTIAL",  "V",  R_POTENTIAL } ,
    { "R_FILL",       "F",  R_FILL } ,          /* 29 */
    { "R_SHEAR_RATE", "SH", R_SHEAR_RATE } ,      

    { "R_PMOMENTUM1", "PVX", R_PMOMENTUM1 } ,
    { "R_PMOMENTUM2", "PVY", R_PMOMENTUM2 } ,
    { "R_PMOMENTUM3", "PVZ", R_PMOMENTUM3 } ,   /* 33 */

    { "R_STRESS11_1", "S11_1", R_STRESS11_1 } ,
    { "R_STRESS12_1", "S12_1", R_STRESS12_1 } ,
    { "R_STRESS22_1", "S22_1", R_STRESS22_1 } ,
    { "R_STRESS13_1", "S13_1", R_STRESS13_1 } ,
    { "R_STRESS23_1", "S23_1", R_STRESS23_1 } ,
    { "R_STRESS33_1", "S33_1", R_STRESS33_1 } , /* 39 */

    { "R_STRESS11_2", "S11_2", R_STRESS11_2 } ,    
    { "R_STRESS12_2", "S12_2", R_STRESS12_2 } ,
    { "R_STRESS22_2", "S22_2", R_STRESS22_2 } , 
    { "R_STRESS13_2", "S13_2", R_STRESS13_2 } ,
    { "R_STRESS23_2", "S23_2", R_STRESS23_2 } ,
    { "R_STRESS33_2", "S33_2", R_STRESS33_2 } , 

    { "R_STRESS11_3", "S11_3", R_STRESS11_3 } ,
    { "R_STRESS12_3", "S12_3", R_STRESS12_3 } ,
    { "R_STRESS22_3", "S22_3", R_STRESS22_3 } , 
    { "R_STRESS13_3", "S13_3", R_STRESS13_3 } , /* 49 */
    { "R_STRESS23_3", "S23_3", R_STRESS23_3 } ,   
    { "R_STRESS33_3", "S33_3", R_STRESS33_3 } , 

    { "R_STRESS11_4", "S11_4", R_STRESS11_4 } ,
    { "R_STRESS12_4", "S12_4", R_STRESS12_4 } ,
    { "R_STRESS22_4", "S22_4", R_STRESS22_4 } , 
    { "R_STRESS13_4", "S13_4", R_STRESS13_4 } ,
    { "R_STRESS23_4", "S23_4", R_STRESS23_4 } ,
    { "R_STRESS33_4", "S33_4", R_STRESS33_4 } , 

    { "R_STRESS11_5", "S11_5", R_STRESS11_5 } ,
    { "R_STRESS12_5", "S12_5", R_STRESS12_5 } , /* 59 */
    { "R_STRESS22_5", "S22_5", R_STRESS22_5 } ,   
    { "R_STRESS13_5", "S13_5", R_STRESS13_5 } ,
    { "R_STRESS23_5", "S23_5", R_STRESS23_5 } ,
    { "R_STRESS33_5", "S33_5", R_STRESS33_5 } , 

    { "R_STRESS11_6", "S11_6", R_STRESS11_6 } ,
    { "R_STRESS12_6", "S12_6", R_STRESS12_6 } ,
    { "R_STRESS22_6", "S22_6", R_STRESS22_6 } , 
    { "R_STRESS13_6", "S13_6", R_STRESS13_6 } ,
    { "R_STRESS23_6", "S23_6", R_STRESS23_6 } ,
    { "R_STRESS33_6", "S33_6", R_STRESS33_6 } , /* 69 */

    { "R_STRESS11_7", "S11_7", R_STRESS11_7 } ,  
    { "R_STRESS12_7", "S12_7", R_STRESS12_7 } ,  
    { "R_STRESS22_7", "S22_7", R_STRESS22_7 } , 
    { "R_STRESS13_7", "S13_7", R_STRESS13_7 } ,
    { "R_STRESS23_7", "S23_7", R_STRESS23_7 } ,
    { "R_STRESS33_7", "S33_7", R_STRESS33_7 } , /* 75 */

    { "R_Species_0",   "Sp_0",   SPECIES_UNK_0  } ,  
    { "R_Species_1",   "Sp_1",   SPECIES_UNK_1  } ,  
    { "R_Species_2",   "Sp_2",   SPECIES_UNK_2  } ,  
    { "R_Species_3",   "Sp_3",   SPECIES_UNK_3  } ,  /* 79 */
    { "R_Species_4",   "Sp_4",   SPECIES_UNK_4  } ,  
    { "R_Species_5",   "Sp_5",   SPECIES_UNK_5  } , 
    { "R_Species_6",   "Sp_6",   SPECIES_UNK_6  } ,  
    { "R_Species_7",   "Sp_7",   SPECIES_UNK_7  } ,  
    { "R_Species_8",   "Sp_8",   SPECIES_UNK_8  } ,  
    { "R_Species_9",   "Sp_9",   SPECIES_UNK_9  } , 
    { "R_Species_10",  "Sp_10",  SPECIES_UNK_10 } , 
    { "R_Species_11",  "Sp_11",  SPECIES_UNK_11 } , 
    { "R_Species_12",  "Sp_12",  SPECIES_UNK_12 } ,  
    { "R_Species_13",  "Sp_13",  SPECIES_UNK_13 } , /* 89 */
    { "R_Species_14",  "Sp_14",  SPECIES_UNK_14 } ,  
    { "R_Species_15",  "Sp_15",  SPECIES_UNK_15 } ,  
    { "R_Species_16",  "Sp_16",  SPECIES_UNK_16 } , 
    { "R_Species_17",  "Sp_17",  SPECIES_UNK_17 } ,  
    { "R_Species_18",  "Sp_18",  SPECIES_UNK_18 } ,  
    { "R_Species_19",  "Sp_19",  SPECIES_UNK_19 } ,  
    { "R_Species_20",  "Sp_20",  SPECIES_UNK_20 } ,  
    { "R_Species_21",  "Sp_21",  SPECIES_UNK_21 } ,  
    { "R_Species_22",  "Sp_22",  SPECIES_UNK_22 } ,  
    { "R_Species_23",  "Sp_23",  SPECIES_UNK_23 } , /* 99 */
    { "R_Species_24",  "Sp_24",  SPECIES_UNK_24 } ,  
    { "R_Species_25",  "Sp_25",  SPECIES_UNK_25 } ,  
    { "R_Species_26",  "Sp_26",  SPECIES_UNK_26 } ,  
    { "R_Species_27",  "Sp_27",  SPECIES_UNK_27 } ,  
    { "R_Species_28",  "Sp_28",  SPECIES_UNK_28 } ,  
    { "R_Species_29",  "Sp_29",  SPECIES_UNK_29 } ,  /* 105 */

    { "R_VolFracPh_0",  "VFP_0", VOLF_PHASE_0   } ,  
    { "R_VolFracPh_1",  "VFP_1", VOLF_PHASE_1   } , 
    { "R_VolFracPh_2",  "VFP_2", VOLF_PHASE_2   } , 
    { "R_VolFracPh_3",  "VFP_3", VOLF_PHASE_3   } ,  /* 109 */
    { "R_VolFracPh_4",  "VFP_4", VOLF_PHASE_4   } ,  
    
    { "R_POR_LIQ_PRES",   "P_LIQ", POR_LIQ_PRES   } ,  
    { "R_POR_GAS_PRES",   "P_GAS", POR_GAS_PRES   } ,  
    { "R_POR_POROSITY",   "P_POR", POR_POROSITY   } ,  
    { "R_POR_ENERGY",     "P_TEMP", POR_TEMP      } , 
    { "R_POR_SATURATION", "P_SAT", POR_SATURATION } ,

    { "R_POR_SINK_MASS", "P_SINK_MASS", POR_SINK_MASS},  /*116 */

    { "R_VORT_DIR1",   "VDX",     R_VORT_DIR1 } ,    
    { "R_VORT_DIR2",   "VDX",     R_VORT_DIR1 } ,
    { "R_VORT_DIR3",   "VDX",     R_VORT_DIR1 } ,
    { "R_VORT_LAMBDA", "VLAMBDA", R_VORT_LAMBDA } ,  /* 120 */

    { "R_CURVATURE", "H",    R_CURVATURE},   

    { "R_LAGR_MULT1",   "LM1", R_LAGR_MULT1 }, 
    { "R_LAGR_MULT2",   "LM2", R_LAGR_MULT2 },
    { "R_LAGR_MULT3",   "LM3", R_LAGR_MULT3 },

    { "R_BOND_EVOLUTION",   "NN", R_BOND_EVOLUTION },

    { "R_SURF_CHARGE", "QS", R_SURF_CHARGE} ,

    { "R_EXT_VELOCITY", "V_EXT", R_EXT_VELOCITY},
    
    { "R_EFIELD1", "E1", R_EFIELD1},
    { "R_EFIELD2", "E2", R_EFIELD2},
    { "R_EFIELD3", "E3", R_EFIELD3},                /* 130 */

    { "R_ENORM", "ENORM", R_ENORM }, 

    { "R_NORMAL1", "N1", R_NORMAL1 },  
    { "R_NORMAL2", "N2", R_NORMAL2 },  
    { "R_NORMAL3", "N3", R_NORMAL3 },  

    { "R_SHELL_CURVATURE", "K", R_SHELL_CURVATURE}, 
    { "R_SHELL_CURVATURE2", "K2", R_SHELL_CURVATURE2}, 
    { "R_SHELL_TENSION", "TENS", R_SHELL_TENSION},
    { "R_SHELL_X", "SH_X", R_SHELL_X},                
    { "R_SHELL_Y", "SH_Y", R_SHELL_Y},
    { "R_SHELL_USER", "SH_U", R_SHELL_USER},        
    { "R_PHASE1", "F1", R_PHASE1},                  /* 141 */  
    { "R_PHASE2", "F2", R_PHASE2},  
    { "R_PHASE3", "F3", R_PHASE3},  
    { "R_PHASE4", "F4", R_PHASE4},  
    { "R_PHASE5", "F5", R_PHASE5}, 
    { "R_SHELL_ANGLE1", "SH_ANG1", R_SHELL_ANGLE1},
    { "R_SHELL_ANGLE2", "SH_ANG2", R_SHELL_ANGLE2},  

    { "R_SHELL_SURF_DIV_V", "SH_DIV_V", R_SHELL_SURF_DIV_V},
    { "R_SHELL_SURF_CURV", "SH_CURV", R_SHELL_SURF_CURV},
    { "R_N_DOT_CURL_V", "SH_NCURLV", R_N_DOT_CURL_V},
    { "R_GRAD_S_V_DOT_N1", "SH_GRADSVNX", R_GRAD_S_V_DOT_N1}, /* 151 */
    { "R_GRAD_S_V_DOT_N2", "SH_GRADSVNY", R_GRAD_S_V_DOT_N2},
    { "R_GRAD_S_V_DOT_N3", "SH_GRADSVNZ", R_GRAD_S_V_DOT_N3},
    { "R_ACOUS_PREAL", "ACOUS_PREAL", R_ACOUS_PREAL},
    { "R_ACOUS_PIMAG", "ACOUS_PIMAG", R_ACOUS_PIMAG},
    { "R_SHELL_DIFF_FLUX", "SH_J", R_SHELL_DIFF_FLUX},   /* 156 */
    { "R_SHELL_DIFF_CURVATURE", "SH_KD", R_SHELL_DIFF_CURVATURE},
    { "R_SHELL_NORMAL1", "SH_NX", R_SHELL_NORMAL1},
    { "R_SHELL_NORMAL2", "SH_NY", R_SHELL_NORMAL2},      
    { "R_SHELL_NORMAL3", "SH_NZ", R_SHELL_NORMAL3},      /* 160 */
    { "R_ACOUS_REYN_STRESS", "ACOUS_REYN_STRESS", R_ACOUS_REYN_STRESS},
    { "R_SHELL_BDYVELO", "SHELL_BDYVELO", R_SHELL_BDYVELO},
    { "R_SHELL_LUBP", "SHELL_LUBP", R_SHELL_LUBP},
    { "R_LUBP", "LUBP", R_LUBP},
    { "R_SHELL_FILMP", "SHELL_FILMP", R_SHELL_FILMP},
    { "R_SHELL_FILMH", "SHELL_FILMH", R_SHELL_FILMH},
    { "R_SHELL_PARTC", "SHELL_PARTC", R_SHELL_PARTC},   
    { "R_SHELL_SAT_CLOSED", "SHELL_SAT_CLOSED", R_SHELL_SAT_CLOSED},   /*  168 - SAR  */
    { "R_SHELL_SAT_OPEN", "SHELL_PRESS_OPEN", R_SHELL_SAT_OPEN},   /*  169 - SAR  */
    { "R_SHELL_ENERGY", "SH_T", R_SHELL_ENERGY},   /*  170 - PRS  */
    { "R_SHELL_DELTAH", "SH_DH", R_SHELL_DELTAH},   /*  171 - PRS  */
    { "R_SHELL_LUB_CURV", "SH_L_CURV", R_SHELL_LUB_CURV},  /*  172 - SAR  */
    { "R_SHELL_SAT_GASN", "SH_SAT_GASN", R_SHELL_SAT_GASN},
    { "R_SHELL_SHEAR_TOP", "SH_SHEAR_TOP", R_SHELL_SHEAR_TOP},   /*  174 */
    { "R_SHELL_SHEAR_BOT", "SH_SHEAR_BOT", R_SHELL_SHEAR_BOT},   /*  175  */
    { "R_SHELL_CROSS_SHEAR", "SH_CROSS_SHEAR", R_SHELL_CROSS_SHEAR},   /*  176  */
    { "R_MAX_STRAIN", "MAX_STRAIN", R_MAX_STRAIN},   /*  177 - SAR */
    { "R_CUR_STRAIN", "CUR_STRAIN", R_CUR_STRAIN},   /*  178 - SAR */
    { "R_LUBP_2", "LUBP_2", R_LUBP_2},               /*  179 - PRS */
    { "R_SHELL_SAT_OPEN_2", "SHELL_PRESS_OPEN_2", R_SHELL_SAT_OPEN_2},   /*  180 - PRS */
    { "R_SHELL_LUB_CURV_2", "SH_L_CURV_2", R_SHELL_LUB_CURV_2},  /*  181 - PRS  */
    { "R_LIGHT_INTP", "LIGHT_INTP", R_LIGHT_INTP}, /* 182 */
    { "R_LIGHT_INTM", "LIGHT_INTM", R_LIGHT_INTM}, /* 183 */
    { "R_LIGHT_INTD", "LIGHT_INTD", R_LIGHT_INTD},  /*   184  */
    { "R_TFMP_MASS", "TFMP_MASS", R_TFMP_MASS},     /*   185  */
    { "R_TFMP_BOUND", "TFMP_BOUND", R_TFMP_BOUND},     /*   186  */
    { "R_RESTIME", "RESTIME", R_RESTIME},     /*   187  */

    { "R_EM_E1_REAL", "EM_E1_REAL", R_EM_E1_REAL},
    { "R_EM_E1_IMAG", "EM_E1_IMAG", R_EM_E1_IMAG},
    { "R_EM_E2_REAL", "EM_E2_REAL", R_EM_E2_REAL},
    { "R_EM_E2_IMAG", "EM_E2_IMAG", R_EM_E2_IMAG}, /*   191  */
    { "R_EM_E3_REAL", "EM_E3_REAL", R_EM_E3_REAL},
    { "R_EM_E3_IMAG", "EM_E3_IMAG", R_EM_E3_IMAG},
    { "R_EM_H1_REAL", "EM_H1_REAL", R_EM_H1_REAL},
    { "R_EM_H1_IMAG", "EM_H1_IMAG", R_EM_H1_IMAG},
    { "R_EM_H2_REAL", "EM_H2_REAL", R_EM_H2_REAL},/*   196  */
    { "R_EM_H2_IMAG", "EM_H2_IMAG", R_EM_H2_IMAG},
    { "R_EM_H3_REAL", "EM_H3_REAL", R_EM_H3_REAL},
    { "R_EM_H3_IMAG", "EM_H3_IMAG", R_EM_H3_IMAG},/*   199  */
    /*
     *  Note -> these entries must remain until we get rid
     *          of putting the species unknowns after V_LAST
     *          There must be at least as many of these as there
     *          are species in the problem
     */

    { "R_Y0", "Y0", V_LAST + 0} ,                    /* 200 */
    { "R_Y1", "Y1", V_LAST + 1} ,
    { "R_Y2", "Y2", V_LAST + 2} ,
    { "R_Y3", "Y3", V_LAST + 3} ,
    { "R_Y4", "Y4", V_LAST + 4} ,
    { "R_Y5", "Y5", V_LAST + 5} ,
    { "R_Y6", "Y6", V_LAST + 6} ,
    { "R_Y7", "Y7", V_LAST + 7} ,
    { "R_Y8", "Y8", V_LAST + 8} ,
    { "R_Y9", "Y9", V_LAST + 9} ,                    
    { "R_Y10", "Y10", V_LAST + 10} ,                 /* 210 */
    { "R_Y11", "Y11", V_LAST + 11} ,                   
    { "R_Y12", "Y12", V_LAST + 12} ,                  
    { "R_Y13", "Y13", V_LAST + 13} ,                  
    { "R_Y14", "Y14", V_LAST + 14} ,                 
    { "R_Y15", "Y15", V_LAST + 15} ,                  
    { "R_Y16", "Y16", V_LAST + 16} ,                  
    { "R_Y17", "Y17", V_LAST + 17} ,                 
    { "R_Y18", "Y18", V_LAST + 18} ,                   
    { "R_Y19", "Y19", V_LAST + 19} ,                 
    { "R_Y20", "Y20", V_LAST + 20} ,                 /* 220 */  
    { "R_Y21", "Y21", V_LAST + 21} ,
    { "R_Y22", "Y22", V_LAST + 22} ,
    { "R_Y23", "Y23", V_LAST + 23} ,
    { "R_Y24", "Y24", V_LAST + 24} ,
    { "R_Y25", "Y25", V_LAST + 25} ,
    { "R_Y26", "Y26", V_LAST + 26} ,
    { "R_Y27", "Y27", V_LAST + 27} ,   
    { "R_Y28", "Y28", V_LAST + 28} ,
    { "R_Y29", "Y29", V_LAST + 29} ,                 /* 229 */

    /*
     * Add extra equation names for vector fields that can be rotated
     */
    
    { "R_MOM_NORMAL",  "DN",  R_MOM_NORMAL } ,       /* 230 */
    { "R_MOM_TANG1",   "DT1", R_MOM_TANG1 } ,
    { "R_MOM_TANG2",   "DT2", R_MOM_TANG2 } ,
    { "R_MESH_NORMAL", "VN",  R_MESH_NORMAL } ,
    { "R_MESH_TANG1",  "VT1", R_MESH_TANG1 } ,
    { "R_MESH_TANG2",  "VT2", R_MESH_TANG2 } ,       /* 235 */
    { "R_SOLID_NORMAL", "SN",  R_SOLID_NORMAL } ,
    { "R_SOLID_TANG1",  "ST1", R_SOLID_TANG1 } ,
    { "R_SOLID_TANG2",  "ST2", R_SOLID_TANG2 }         /* 238 */
};
int Num_EQ_Names = sizeof(EQ_Name) / sizeof(struct Equation_Names);  

/* structure to control input of variable names 
 * using values in rf_fem_const.h */

struct Equation_Names Var_Name[] =  {
    { "VELOCITY1",          "VX",  VELOCITY1 } ,	/* 0 */   
    { "VELOCITY2",          "VY",  VELOCITY2 } ,
    { "VELOCITY3",          "VZ",  VELOCITY3 } ,
    { "TEMPERATURE",        "T",  TEMPERATURE } ,
    { "MASS_FRACTION",      "Y",  MASS_FRACTION } ,
    { "MESH_DISPLACEMENT1", "DMX", MESH_DISPLACEMENT1 } ,
    { "MESH_DISPLACEMENT2", "DMY", MESH_DISPLACEMENT2 } ,
    { "MESH_DISPLACEMENT3", "DMZ", MESH_DISPLACEMENT3 } ,
    { "SURFACE",            "S", SURFACE } ,
    { "PRESSURE",           "P", PRESSURE } ,		/* 09 */

    { "POLYMER_STRESS11", "S11", POLYMER_STRESS11 } ,
    { "POLYMER_STRESS12", "S12", POLYMER_STRESS12 } ,
    { "POLYMER_STRESS22", "S22", POLYMER_STRESS22 } ,
    { "POLYMER_STRESS13", "S13", POLYMER_STRESS13 } ,
    { "POLYMER_STRESS23", "S23", POLYMER_STRESS23 } ,
    { "POLYMER_STRESS33", "S33", POLYMER_STRESS33 } ,

    { "SOLID_DISPLACEMENT1", "DMX_RS", SOLID_DISPLACEMENT1 } ,
    { "SOLID_DISPLACEMENT2", "DMY_RS", SOLID_DISPLACEMENT2 } ,
    { "SOLID_DISPLACEMENT3", "DMZ_RS", SOLID_DISPLACEMENT3 } ,
		    
    { "VELOCITY_GRADIENT11", "G11", VELOCITY_GRADIENT11 } , /* 19 */
    { "VELOCITY_GRADIENT12", "G12", VELOCITY_GRADIENT12 } ,
    { "VELOCITY_GRADIENT21", "G21", VELOCITY_GRADIENT21 } ,
    { "VELOCITY_GRADIENT22", "G22", VELOCITY_GRADIENT22 } ,
    { "VELOCITY_GRADIENT13", "G13", VELOCITY_GRADIENT13 } ,
    { "VELOCITY_GRADIENT23", "G23", VELOCITY_GRADIENT23 } ,
    { "VELOCITY_GRADIENT31", "G31", VELOCITY_GRADIENT31 } ,
    { "VELOCITY_GRADIENT32", "G32", VELOCITY_GRADIENT32 } ,
    { "VELOCITY_GRADIENT33", "G33", VELOCITY_GRADIENT33 } ,

    { "VOLTAGE",      "V", VOLTAGE } ,
    { "FILL",         "F", FILL } ,			/* 29 */
    { "SHEAR_RATE",   "SH", SHEAR_RATE } ,
    { "PVELOCITY1",          "PVX",  PVELOCITY1 } ,
    { "PVELOCITY2",          "PVY",  PVELOCITY2 } ,
    { "PVELOCITY3",          "PVZ",  PVELOCITY3 } ,

    { "POLYMER_STRESS11_1", "S11_1", POLYMER_STRESS11_1 } ,
    { "POLYMER_STRESS12_1", "S12_1", POLYMER_STRESS12_1 } ,
    { "POLYMER_STRESS22_1", "S22_1", POLYMER_STRESS22_1 } ,
    { "POLYMER_STRESS13_1", "S13_1", POLYMER_STRESS13_1 } ,
    { "POLYMER_STRESS23_1", "S23_1", POLYMER_STRESS23_1 } ,
    { "POLYMER_STRESS33_1", "S33_1", POLYMER_STRESS33_1 } , /* 39 */

    { "POLYMER_STRESS11_2", "S11_2", POLYMER_STRESS11_2 } ,
    { "POLYMER_STRESS12_2", "S12_2", POLYMER_STRESS12_2 } ,
    { "POLYMER_STRESS22_2", "S22_2", POLYMER_STRESS22_2 } ,
    { "POLYMER_STRESS13_2", "S13_2", POLYMER_STRESS13_2 } ,
    { "POLYMER_STRESS23_2", "S23_2", POLYMER_STRESS23_2 } ,
    { "POLYMER_STRESS33_2", "S33_2", POLYMER_STRESS33_2 } ,

    { "POLYMER_STRESS11_3", "S11_3", POLYMER_STRESS11_3 } ,
    { "POLYMER_STRESS12_3", "S12_3", POLYMER_STRESS12_3 } ,
    { "POLYMER_STRESS22_3", "S22_3", POLYMER_STRESS22_3 } ,
    { "POLYMER_STRESS13_3", "S13_3", POLYMER_STRESS13_3 } , /* 49 */
    { "POLYMER_STRESS23_3", "S23_3", POLYMER_STRESS23_3 } ,
    { "POLYMER_STRESS33_3", "S33_3", POLYMER_STRESS33_3 } ,

    { "POLYMER_STRESS11_4", "S11_4", POLYMER_STRESS11_4 } ,
    { "POLYMER_STRESS12_4", "S12_4", POLYMER_STRESS12_4 } ,
    { "POLYMER_STRESS22_4", "S22_4", POLYMER_STRESS22_4 } ,
    { "POLYMER_STRESS13_4", "S13_4", POLYMER_STRESS13_4 } ,
    { "POLYMER_STRESS23_4", "S23_4", POLYMER_STRESS23_4 } ,
    { "POLYMER_STRESS33_4", "S33_4", POLYMER_STRESS33_4 } ,

    { "POLYMER_STRESS11_5", "S11_5", POLYMER_STRESS11_5 } ,
    { "POLYMER_STRESS12_5", "S12_5", POLYMER_STRESS12_5 } , /* 59 */
    { "POLYMER_STRESS22_5", "S22_5", POLYMER_STRESS22_5 } ,
    { "POLYMER_STRESS13_5", "S13_5", POLYMER_STRESS13_5 } ,
    { "POLYMER_STRESS23_5", "S23_5", POLYMER_STRESS23_5 } ,
    { "POLYMER_STRESS33_5", "S33_5", POLYMER_STRESS33_5 } ,

    { "POLYMER_STRESS11_6", "S11_6", POLYMER_STRESS11_6 } ,
    { "POLYMER_STRESS12_6", "S12_6", POLYMER_STRESS12_6 } ,
    { "POLYMER_STRESS22_6", "S22_6", POLYMER_STRESS22_6 } ,
    { "POLYMER_STRESS13_6", "S13_6", POLYMER_STRESS13_6 } ,
    { "POLYMER_STRESS23_6", "S23_6", POLYMER_STRESS23_6 } ,
    { "POLYMER_STRESS33_6", "S33_6", POLYMER_STRESS33_6 } , /* 69 */

    { "POLYMER_STRESS11_7", "S11_7", POLYMER_STRESS11_7 } ,
    { "POLYMER_STRESS12_7", "S12_7", POLYMER_STRESS12_7 } ,
    { "POLYMER_STRESS22_7", "S22_7", POLYMER_STRESS22_7 } ,
    { "POLYMER_STRESS13_7", "S13_7", POLYMER_STRESS13_7 } ,
    { "POLYMER_STRESS23_7", "S23_7", POLYMER_STRESS23_7 } ,
    { "POLYMER_STRESS33_7", "S33_7", POLYMER_STRESS33_7 } ,  /* 75 */  

    { "Species_Conc_0",   "Sp_0",   SPECIES_UNK_0  } ,  
    { "Species_Conc_1",   "Sp_1",   SPECIES_UNK_1  } ,  
    { "Species_Conc_2",   "Sp_2",   SPECIES_UNK_2  } ,  
    { "Species_Conc_3",   "Sp_3",   SPECIES_UNK_3  } ,	/* 79 */
    { "Species_Conc_4",   "Sp_4",   SPECIES_UNK_4  } ,  
    { "Species_Conc_5",   "Sp_5",   SPECIES_UNK_5  } , 
    { "Species_Conc_6",   "Sp_6",   SPECIES_UNK_6  } ,  
    { "Species_Conc_7",   "Sp_7",   SPECIES_UNK_7  } ,  
    { "Species_Conc_8",   "Sp_8",   SPECIES_UNK_8  } ,  
    { "Species_Conc_9",   "Sp_9",   SPECIES_UNK_9  } , 
    { "Species_Conc_10",  "Sp_10",  SPECIES_UNK_10 } , 
    { "Species_Conc_11",  "Sp_11",  SPECIES_UNK_11 } , 
    { "Species_Conc_12",  "Sp_12",  SPECIES_UNK_12 } ,  
    { "Species_Conc_13",  "Sp_13",  SPECIES_UNK_13 } ,	/* 89 */
    { "Species_Conc_14",  "Sp_14",  SPECIES_UNK_14 } ,  
    { "Species_Conc_15",  "Sp_15",  SPECIES_UNK_15 } ,  
    { "Species_Conc_16",  "Sp_16",  SPECIES_UNK_16 } , 
    { "Species_Conc_17",  "Sp_17",  SPECIES_UNK_17 } ,  
    { "Species_Conc_18",  "Sp_18",  SPECIES_UNK_18 } ,  
    { "Species_Conc_19",  "Sp_19",  SPECIES_UNK_19 } ,  
    { "Species_Conc_20",  "Sp_20",  SPECIES_UNK_20 } ,  
    { "Species_Conc_21",  "Sp_21",  SPECIES_UNK_21 } ,  
    { "Species_Conc_22",  "Sp_22",  SPECIES_UNK_22 } ,  
    { "Species_Conc_23",  "Sp_23",  SPECIES_UNK_23 } ,	/* 99 */
    { "Species_Conc_24",  "Sp_24",  SPECIES_UNK_24 } ,  
    { "Species_Conc_25",  "Sp_25",  SPECIES_UNK_25 } ,  
    { "Species_Conc_26",  "Sp_26",  SPECIES_UNK_26 } ,  
    { "Species_Conc_27",  "Sp_27",  SPECIES_UNK_27 } ,  
    { "Species_Conc_28",  "Sp_28",  SPECIES_UNK_28 } ,  
    { "Species_Conc_29",  "Sp_29",  SPECIES_UNK_29 } ,	/* 105 */
    
    { "VolFracPh_0",  "VFP_0", VOLF_PHASE_0   } ,  
    { "VolFracPh_1",  "VFP_1", VOLF_PHASE_1   } , 
    { "VolFracPh_2",  "VFP_2", VOLF_PHASE_2   } , 
    { "VolFracPh_3",  "VFP_3", VOLF_PHASE_3   } ,  	/* 109 */
    { "VolFracPh_4",  "VFP_4", VOLF_PHASE_4   } ,  
    
    { "POR_LIQ_PRES",   "P_LIQ", POR_LIQ_PRES   } ,  
    { "POR_GAS_PRES",   "P_GAS", POR_GAS_PRES   } ,           
    { "POR_POROSITY",    "P_POR", POR_POROSITY   } ,  
    { "POR_TEMP",       "P_TEMP", POR_TEMP  } ,  
    { "POR_SATURATION",  "P_SAT", POR_SATURATION } ,

    { "POR_SINK_MASS", "SINK_MASS", POR_SINK_MASS} , 
  
    { "VORT_DIR1",   "VDX",  VORT_DIR1 } ,         
    { "VORT_DIR2",   "VDY",  VORT_DIR2 } ,
    { "VORT_DIR3",   "VDZ",  VORT_DIR3 } , 
    { "VORT_LAMBDA", "VLAMBDA",  VORT_LAMBDA } ,	/* 120 */

    { "CURVATURE", "H", CURVATURE},       

    { "LAGR_MULT1",   "LM1",  LAGR_MULT1 } ,         
    { "LAGR_MULT2",   "LM2",  LAGR_MULT2 } ,
    { "LAGR_MULT3",   "LM3",  LAGR_MULT3 } , 

    { "BOND_EVOLUTION",   "NN",  BOND_EVOLUTION } ,  

    { "SURF_CHARGE", "QS", SURF_CHARGE } ,

    { "EXT_VELOCITY", "V_EXT", EXT_VELOCITY},

    { "EFIELD1", "E1", EFIELD1},
    { "EFIELD2", "E2", EFIELD2},
    { "EFIELD3", "E3", EFIELD3},			/* 130 */

    { "ENORM", "ENORM", ENORM }, 

    { "NORMAL1", "N1", NORMAL1 },  
    { "NORMAL2", "N2", NORMAL2 },  
    { "NORMAL3", "N3", NORMAL3 },  
 
    { "SHELL_CURVATURE", "K", SHELL_CURVATURE}, 
    { "SHELL_CURVATURE2", "K2", SHELL_CURVATURE2}, 
    { "SHELL_TENSION", "TENS", SHELL_TENSION},
    { "SHELL_X", "SH_X", SHELL_X},                 
    { "SHELL_Y", "SH_Y", SHELL_Y},
    { "SHELL_USER", "SH_U", SHELL_USER},

    { "PHASE1", "F1", PHASE1},				/* 141 */
    { "PHASE2", "F2", PHASE2},
    { "PHASE3", "F3", PHASE3},
    { "PHASE4", "F4", PHASE4},
    { "PHASE5", "F5", PHASE5},
    
    { "SHELL_ANGLE1", "SH_ANG1", SHELL_ANGLE1},
    { "SHELL_ANGLE2", "SH_ANG2", SHELL_ANGLE2},

    { "SHELL_SURF_DIV_V", "SH_DIV_V", SHELL_SURF_DIV_V},
    { "SHELL_SURF_CURV", "SH_CURV", SHELL_SURF_CURV},
    { "N_DOT_CURL_V", "SH_NCURLV", N_DOT_CURL_V},
    { "GRAD_S_V_DOT_NX", "SH_GRADSVNX", GRAD_S_V_DOT_N1},	/* 151 */
    { "GRAD_S_V_DOT_NY", "SH_GRADSVNY", GRAD_S_V_DOT_N2},
    { "GRAD_S_V_DOT_NZ", "SH_GRADSVNZ", GRAD_S_V_DOT_N3},
    { "ACOUS_PREAL", "APR", ACOUS_PREAL},
    { "ACOUS_PIMAG", "API", ACOUS_PIMAG},
    { "SHELL_DIFF_FLUX", "SH_J", R_SHELL_DIFF_FLUX},
    { "SHELL_DIFF_CURVATURE", "SH_KD", R_SHELL_DIFF_CURVATURE},
    { "SHELL_NORMAL1", "SH_NX", R_SHELL_NORMAL1},
    { "SHELL_NORMAL2", "SH_NY", R_SHELL_NORMAL2},	
    { "SHELL_NORMAL3", "SH_NZ", R_SHELL_NORMAL3},	/* 160 */
    { "ACOUS_REYN_STRESS", "ARS", ACOUS_REYN_STRESS},
    { "SHELL_BDYVELO", "SH_BV", SHELL_BDYVELO},
    { "SHELL_LUBP", "SH_P", SHELL_LUBP},
    { "LUBP", "LUB_P", LUBP},
    { "SHELL_FILMP", "SH_FP", SHELL_FILMP},
    { "SHELL_FILMH", "SH_FH", SHELL_FILMH},
    { "SHELL_PARTC", "SH_PC", SHELL_PARTC},
    { "SHELL_SAT_CLOSED", "SH_SAT_CLOSED", SHELL_SAT_CLOSED},  /*  168 - SAR  */
    { "SHELL_PRESS_OPEN", "SH_P_OPEN", SHELL_PRESS_OPEN},        /*  169 - SAR  */
    { "SHELL_TEMPERATURE", "SH_TEMP", SHELL_TEMPERATURE},       /*  170 - PRS  */
    { "SHELL_DELTAH", "SH_DH", SHELL_DELTAH},                  /*  171 - PRS  */
    { "SHELL_LUB_CURV", "SH_L_CURV", SHELL_LUB_CURV},          /* 172 - SAR */
    { "SHELL_SAT_GASN", "SH_SAT_GASN", SHELL_SAT_GASN},
    { "SHELL_SHEAR_TOP", "SH_SHEAR_TOP", SHELL_SHEAR_TOP},
    { "SHELL_SHEAR_BOT", "SH_SHEAR_BOT", SHELL_SHEAR_BOT},
    { "SHELL_CROSS_SHEAR", "SH_CROSS_SHEAR", SHELL_CROSS_SHEAR},
    { "MAX_STRAIN", "MAX_STRAIN", MAX_STRAIN},  /* 177 - SAR */
    { "CUR_STRAIN", "CUR_STRAIN", CUR_STRAIN},  /* 178 - SAR */
    { "LUBP_2", "LUB_P_2", LUBP_2},             /* 179 - PRS */
    { "SHELL_PRESS_OPEN_2", "SH_P_OPEN_2", SHELL_PRESS_OPEN_2},        /*  180 - PRS  */
    { "SHELL_LUB_CURV_2", "SH_L_CURV_2", SHELL_LUB_CURV_2},          /* 181 - PRS */
    { "LIGHT_INTP", "INTP", LIGHT_INTP},
    { "LIGHT_INTM", "INTM", LIGHT_INTM},
    { "LIGHT_INTD", "INTD", LIGHT_INTD},                     /* 184 */
    { "TFMP_SAT", "SAT", TFMP_SAT},                     /* 185 */
    { "TFMP_PRES", "PRES", TFMP_PRES},                     /* 186 */
    { "RESTIME", "RST", RESTIME},                     /* 187 */

    { "EM_E1_REAL", "E1R", EM_E1_REAL},
    { "EM_E1_IMAG", "E1I", EM_E1_IMAG},
    { "EM_E2_REAL", "E2R", EM_E2_REAL},
    { "EM_E2_IMAG", "E2I", EM_E2_IMAG},
    { "EM_E3_REAL", "E3R", EM_E3_REAL},
    { "EM_E3_IMAG", "E3I", EM_E3_IMAG},
    { "EM_H1_REAL", "H1R", EM_H1_REAL},
    { "EM_H1_IMAG", "H1I", EM_H1_IMAG},
    { "EM_H2_REAL", "H2R", EM_H2_REAL},
    { "EM_H2_IMAG", "H2I", EM_H2_IMAG},
    { "EM_H3_REAL", "H3R", EM_H3_REAL},
    { "EM_H3_IMAG", "H3I", EM_H3_IMAG},

    { "MESH_POSITION1", "X",  MESH_POSITION1 } ,
    { "MESH_POSITION2", "Y",  MESH_POSITION2 } ,	/* 189 */
    { "MESH_POSITION3", "Z",  MESH_POSITION3 } ,

    { "VEL_NORM",       "VN", VEL_NORM } ,
    { "SPEED",       "VMAG", SPEED } ,

    { "D_VEL1_DT", "UDOT", D_VEL1_DT } ,
    { "D_VEL2_DT", "VDOT", D_VEL2_DT } ,
    { "D_VEL3_DT", "WDOT", D_VEL3_DT } ,
    { "D_T_DT",    "TDOT", D_T_DT } ,
    { "D_C_DT",    "CDOT", D_Y_DT } ,

    { "D_X1_DT",   "XDOT", D_X1_DT } ,
    { "D_X2_DT",   "YDOT", D_X2_DT } ,
    { "D_X3_DT",   "ZDOT", D_X3_DT } ,			/* 200 */
    { "D_S_DT",    "SDOT", D_S_DT } ,

    { "D_P_DT",    "PDOT", D_P_DT } ,

    { "SOLID_POSITION1", "X_RS",  SOLID_POSITION1 } ,  
    { "SOLID_POSITION2", "Y_RS",  SOLID_POSITION2 } ,
    { "SOLID_POSITION3", "Z_RS",  SOLID_POSITION3 } 	/* 205 */
};

int Num_Var_Names = sizeof(Var_Name) / sizeof(struct Equation_Names);  

struct Equation_Names Exo_Var_Names[] = 
{
  { "Velocity, x component",           "VX",     VELOCITY1 } ,
  { "Velocity, y component",           "VY",     VELOCITY2 } ,
  { "Velocity, z component",           "VZ",     VELOCITY3 } ,
  { "Temperature",                     "T",      TEMPERATURE } ,
  { "Species concentration",           "Y",      MASS_FRACTION } ,
  { "Mesh displacement, x component",  "DMX",    MESH_DISPLACEMENT1 } ,
  { "Mesh displacement, y component",  "DMY",    MESH_DISPLACEMENT2 } ,
  { "Mesh displacement, z component",  "DMZ",    MESH_DISPLACEMENT3 } ,
  { "SURFACE",                         "S",      SURFACE } ,
  { "Pressure",                        "P",      PRESSURE } ,

  { "Polymer stress, xx component",    "S11",    POLYMER_STRESS11 } ,
  { "Polymer stress, xy component",    "S12",    POLYMER_STRESS12 } ,
  { "Polymer stress, xz component",    "S13",    POLYMER_STRESS13 } ,
  { "Polymer stress, yy component",    "S22",    POLYMER_STRESS22 } ,
  { "Polymer stress, yz component",    "S23",    POLYMER_STRESS23 } ,
  { "Polymer stress, zz component",    "S33",    POLYMER_STRESS33 } ,

  { "SOLID displacement, x component", "DMX_RS", SOLID_DISPLACEMENT1 } ,
  { "SOLID displacement, y component", "DMY_RS", SOLID_DISPLACEMENT2 } ,
  { "SOLID displacement, z component", "DMZ_RS", SOLID_DISPLACEMENT3 } ,
		    
  { "Velocity Gradient, xx component", "G11",    VELOCITY_GRADIENT11 } ,
  { "Velocity Gradient, xy component", "G12",    VELOCITY_GRADIENT12 } ,
  { "Velocity Gradient, xz component", "G13",    VELOCITY_GRADIENT13 } ,
  { "Velocity Gradient, yx component", "G21",    VELOCITY_GRADIENT21 } ,
  { "Velocity Gradient, yy component", "G22",    VELOCITY_GRADIENT22 } ,
  { "Velocity Gradient, yz component", "G23",    VELOCITY_GRADIENT23 } ,
  { "Velocity Gradient, zx component", "G31",    VELOCITY_GRADIENT31 } ,
  { "Velocity Gradient, zy component", "G32",    VELOCITY_GRADIENT32 } ,
  { "Velocity Gradient, zz component", "G33",    VELOCITY_GRADIENT33 } ,

  { "Voltage",                         "VOLT",   VOLTAGE } ,
  { "Fill",                            "F",      FILL } ,
  { "Shear Rate",                      "SH",     SHEAR_RATE },
  { "Particle velocity, x-comp",        "PVX",     PVELOCITY1 },
  { "Particle velocity, y-comp",        "PVY",     PVELOCITY2 },
  { "Particle velocity, z-comp",        "PVZ",    PVELOCITY3 },

  { "Polymer stress component 11 mode 1", "S11_1", POLYMER_STRESS11_1 } ,
  { "Polymer stress component 12 mode 1", "S12_1", POLYMER_STRESS12_1 } ,
  { "Polymer stress component 22 mode 1", "S22_1", POLYMER_STRESS22_1 } ,
  { "Polymer stress component 13 mode 1", "S13_1", POLYMER_STRESS13_1 } ,
  { "Polymer stress component 23 mode 1", "S23_1", POLYMER_STRESS23_1 } ,
  { "Polymer stress component 33 mode 1", "S33_1", POLYMER_STRESS33_1 } ,

  { "Polymer stress component 11 mode 2", "S11_2", POLYMER_STRESS11_2 } ,
  { "Polymer stress component 12 mode 2", "S12_2", POLYMER_STRESS12_2 } ,
  { "Polymer stress component 22 mode 2", "S22_2", POLYMER_STRESS22_2 } ,
  { "Polymer stress component 13 mode 2", "S13_2", POLYMER_STRESS13_2 } ,
  { "Polymer stress component 23 mode 2", "S23_2", POLYMER_STRESS23_2 } ,
  { "Polymer stress component 33 mode 2", "S33_2", POLYMER_STRESS33_2 } ,

  { "Polymer stress component 11 mode 3", "S11_3", POLYMER_STRESS11_3 } ,
  { "Polymer stress component 12 mode 3", "S12_3", POLYMER_STRESS12_3 } ,
  { "Polymer stress component 22 mode 3", "S22_3", POLYMER_STRESS22_3 } ,
  { "Polymer stress component 13 mode 3", "S13_3", POLYMER_STRESS13_3 } ,
  { "Polymer stress component 23 mode 3", "S23_3", POLYMER_STRESS23_3 } ,
  { "Polymer stress component 33 mode 3", "S33_3", POLYMER_STRESS33_3 } ,

  { "Polymer stress component 11 mode 4", "S11_4", POLYMER_STRESS11_4 } ,
  { "Polymer stress component 12 mode 4", "S12_4", POLYMER_STRESS12_4 } ,
  { "Polymer stress component 22 mode 4", "S22_4", POLYMER_STRESS22_4 } ,
  { "Polymer stress component 13 mode 4", "S13_4", POLYMER_STRESS13_4 } ,
  { "Polymer stress component 23 mode 4", "S23_4", POLYMER_STRESS23_4 } ,
  { "Polymer stress component 33 mode 4", "S33_4", POLYMER_STRESS33_4 } ,

  { "Polymer stress component 11 mode 5", "S11_5", POLYMER_STRESS11_5 } ,
  { "Polymer stress component 12 mode 5", "S12_5", POLYMER_STRESS12_5 } ,
  { "Polymer stress component 22 mode 5", "S22_5", POLYMER_STRESS22_5 } ,
  { "Polymer stress component 13 mode 5", "S13_5", POLYMER_STRESS13_5 } ,
  { "Polymer stress component 23 mode 5", "S23_5", POLYMER_STRESS23_5 } ,
  { "Polymer stress component 33 mode 5", "S33_5", POLYMER_STRESS33_5 } ,

  { "Polymer stress component 11 mode 6", "S11_6", POLYMER_STRESS11_6 } ,
  { "Polymer stress component 12 mode 6", "S12_6", POLYMER_STRESS12_6 } ,
  { "Polymer stress component 22 mode 6", "S22_6", POLYMER_STRESS22_6 } ,
  { "Polymer stress component 13 mode 6", "S13_6", POLYMER_STRESS13_6 } ,
  { "Polymer stress component 23 mode 6", "S23_6", POLYMER_STRESS23_6 } ,
  { "Polymer stress component 33 mode 6", "S33_6", POLYMER_STRESS33_6 } ,

  { "Polymer stress component 11 mode 7", "S11_7", POLYMER_STRESS11_7 } ,
  { "Polymer stress component 12 mode 7", "S12_7", POLYMER_STRESS12_7 } ,
  { "Polymer stress component 22 mode 7", "S22_7", POLYMER_STRESS22_7 } ,
  { "Polymer stress component 13 mode 7", "S13_7", POLYMER_STRESS13_7 } ,
  { "Polymer stress component 23 mode 7", "S23_7", POLYMER_STRESS23_7 } ,
  { "Polymer stress component 33 mode 7", "S33_7", POLYMER_STRESS33_7 } ,   /*76*/ 

  { "species concentration comp 0",   "Sp_0",   SPECIES_UNK_0  } ,  
  { "species concentration comp 1",   "Sp_1",   SPECIES_UNK_1  } ,  
  { "species concentration comp 2",   "Sp_2",   SPECIES_UNK_2  } ,  
  { "species concentration comp 3",   "Sp_3",   SPECIES_UNK_3  } ,  
  { "species concentration comp 4",   "Sp_4",   SPECIES_UNK_4  } ,  
  { "species concentration comp 5",   "Sp_5",   SPECIES_UNK_5  } , 
  { "species concentration comp 6",   "Sp_6",   SPECIES_UNK_6  } ,  
  { "species concentration comp 7",   "Sp_7",   SPECIES_UNK_7  } ,  
  { "species concentration comp 8",   "Sp_8",   SPECIES_UNK_8  } ,  
  { "species concentration comp 9",   "Sp_9",   SPECIES_UNK_9  } , 
  { "species concentration comp 10",  "Sp_10",  SPECIES_UNK_10 } , 
  { "species concentration comp 11",  "Sp_11",  SPECIES_UNK_11 } , 
  { "species concentration comp 12",  "Sp_12",  SPECIES_UNK_12 } ,  
  { "species concentration comp 13",  "Sp_13",  SPECIES_UNK_13 } , 
  { "species concentration comp 14",  "Sp_14",  SPECIES_UNK_14 } ,  
  { "species concentration comp 15",  "Sp_15",  SPECIES_UNK_15 } ,  
  { "species concentration comp 16",  "Sp_16",  SPECIES_UNK_16 } , 
  { "species concentration comp 17",  "Sp_17",  SPECIES_UNK_17 } ,  
  { "species concentration comp 18",  "Sp_18",  SPECIES_UNK_18 } ,  
  { "species concentration comp 19",  "Sp_19",  SPECIES_UNK_19 } ,  
  { "species concentration comp 20",  "Sp_20",  SPECIES_UNK_20 } ,  
  { "species concentration comp 21",  "Sp_21",  SPECIES_UNK_21 } ,  
  { "species concentration comp 22",  "Sp_22",  SPECIES_UNK_22 } ,  
  { "species concentration comp 23",  "Sp_23",  SPECIES_UNK_23 } , 
  { "species concentration comp 24",  "Sp_24",  SPECIES_UNK_24 } ,  
  { "species concentration comp 25",  "Sp_25",  SPECIES_UNK_25 } ,  
  { "species concentration comp 26",  "Sp_26",  SPECIES_UNK_26 } ,  
  { "species concentration comp 27",  "Sp_27",  SPECIES_UNK_27 } ,  
  { "species concentration comp 28",  "Sp_28",  SPECIES_UNK_28 } ,  
  { "species concentration comp 29",  "Sp_29",  SPECIES_UNK_29 } ,  /*106*/
  
  { "VolFracPh_0",  "VFP_0", VOLF_PHASE_0   } ,  
  { "VolFracPh_1",  "VFP_1", VOLF_PHASE_1   } , 
  { "VolFracPh_2",  "VFP_2", VOLF_PHASE_2   } , 
  { "VolFracPh_3",  "VFP_3", VOLF_PHASE_3   } ,  
  { "VolFracPh_4",  "VFP_4", VOLF_PHASE_4   } ,  
    
  { "Liquid phase pressue, porous",   "P_LIQ", POR_LIQ_PRES   } ,  
  { "Gas phase pressure, porous",   "P_GAS", POR_GAS_PRES   } ,           
  { "Porosity of medium, porous",    "P_POR", POR_POROSITY   } ,  
  { "Temperature, porous", "P_TEMP", POR_TEMP } ,  
  { "Saturation of medium, porous",  "P_SAT", POR_SATURATION } ,  /*116*/

  { "Sink mass accumulation, porous", "SINK_MASS", POR_SINK_MASS}, 

  { "Vorticity direction, x component", "VDX", VORT_DIR1 } ,
  { "Vorticity direction, y component", "VDY", VORT_DIR2 } ,
  { "Vorticity direction, z component", "VDZ", VORT_DIR3 } ,
  { "Vorticity eigenvalue",             "VLAMBDA", VORT_LAMBDA } ,

  { "Curvature of level set field", "H", CURVATURE },

  { "Volume Lagrange Multiplier, x component", "LMX", LAGR_MULT1 } ,
  { "Volume Lagrange Multiplier, y component", "LMY", LAGR_MULT2 } ,
  { "Volume Lagrange Multiplier, z component", "LMZ", LAGR_MULT3 } ,

  { "Bond Structure", "NN", BOND_EVOLUTION },

  { "Surface charge density for shell elements", "QS", SURF_CHARGE },

  { "Extension velocity, normal component", "EXT_V", EXT_VELOCITY },

  { "Electric field, x component", "E1", EFIELD1 },  /*128*/
  { "Electric field, y component", "E2", EFIELD2 },
  { "Electric field, z component", "E3", EFIELD3 },

  { "Norm of potential field", "ENORM", ENORM },

  { "Normal to LS function, x component", "N1", NORMAL1 },
  { "Normal to LS function, y component", "N2", NORMAL2 },
  { "Normal to LS function, z component", "N3", NORMAL3 },

  { "Shell element curvature", "K", SHELL_CURVATURE},
  { "Shell element second curvature", "K2", SHELL_CURVATURE2},
  { "Shell element tension", "TENS", SHELL_TENSION},
  { "Shell element X coordinate", "SH_X", SHELL_X},
  { "Shell element Y coordinate", "SH_Y", SHELL_Y},
  { "Shell element user", "SH_U", SHELL_USER},
  { "First  phase function", "F1", PHASE1},
  { "Second phase function", "F2", PHASE2},
  { "Third  phase function", "F3", PHASE3},
  { "Fourth phase function", "F4", PHASE4},
  { "Fifth  phase function", "F5", PHASE5},
  
  { "Shell element orientation angle, theta", "SH_ANG1", SHELL_ANGLE1},
  { "Shell element orientation angle, phi", "SH_ANG2", SHELL_ANGLE2},

  { "(I-NN).del.v", "SH_DIV_V", SHELL_SURF_DIV_V},
  { "Mean curvature of surface", "SH_CURV", SHELL_SURF_CURV},
  { "normal component of vorticity", "SH_NCURLV", N_DOT_CURL_V},
  { "Normal x component of velocity gradient_s tensor", "SH_GRADSVNX", GRAD_S_V_DOT_N1},
  { "Normal y component of velocity gradient_s tensor", "SH_GRADSVNY", GRAD_S_V_DOT_N2},
  { "Normal z component of velocity gradient_S tensor", "SH_GRADSVNZ", GRAD_S_V_DOT_N3},
  { "Acoustic harmonic pressure - real part", "APR", ACOUS_PREAL},
  { "Acoustic harmonic pressure - imag part", "API", ACOUS_PIMAG},
  { "Shell surface diffusion flux", "SH_J", SHELL_DIFF_FLUX},
  { "Shell surface curvature", "SH_KD", SHELL_DIFF_CURVATURE},
  { "Shell normal vector, X component", "SH_NX", SHELL_NORMAL1},
  { "Shell normal vector, Y component", "SH_NY", SHELL_NORMAL2},
  { "Shell normal vector, Z component", "SH_NZ", SHELL_NORMAL3},
  { "Acoustic Reynolds Stress", "ARS", ACOUS_REYN_STRESS},
  { "Acoustic Boundary Velocity", "SH_BV", SHELL_BDYVELO},
  { "Shell Lubrication Pressure", "SH_P", SHELL_LUBP},
  { "Lubrication Pressure", "LUBP", LUBP},
  { "Shell Film Pressure", "SH_FP", SHELL_FILMP},
  { "Shell Film Thickness","SH_FH", SHELL_FILMH},
  { "Shell Particles Concentration", "SH_PC", SHELL_PARTC}, 
  { "Structured Porous Shell Saturation - Closed Cells", "SH_SAT_CLOSED", SHELL_SAT_CLOSED },
  { "Structured Porous Shell Pressure - Open Cells", "SH_P_OPEN", SHELL_PRESS_OPEN },
  { "Shell Temperature", "SH_TEMP", SHELL_TEMPERATURE },
  { "Shell delta_gap", "SH_DH", SHELL_DELTAH },
  { "Shell Lubrication Curvature", "SH_L_CURV", SHELL_LUB_CURV },
  { "Structured Porous Shell Saturation - Gas Compression", "SH_SAT_GASN", SHELL_SAT_GASN },
  { "Top wall shear rate", "SH_SHEAR_TOP", SHELL_SHEAR_TOP},
  { "Bottom wall shear rate", "SH_SHEAR_BOT", SHELL_SHEAR_BOT},
  { "Cross stream shear stress", "SH_CROSS_SHEAR", SHELL_CROSS_SHEAR},
  { "Maximum Von Mises strain", "MAX_STRAIN", MAX_STRAIN},
  { "Current Von Mises strain", "CUR_STRAIN", CUR_STRAIN},
  { "Lubrication Pressure 2", "LUBP_2", LUBP_2},
  { "Structured Porous Shell Pressure - Open Cells 2", "SH_P_OPEN_2", SHELL_PRESS_OPEN_2 },
  { "Shell Lubrication Curvature_2", "SH_L_CURV_2", SHELL_LUB_CURV_2 },
  { "Plus Propagating Intensity", "INTP", LIGHT_INTP },
  { "Minus Propagating Intensity", "INTM", LIGHT_INTM },
  { "Dispersive Scattering Intensity", "INTD", LIGHT_INTD },
  { "Thin Film Multiphase Lubrication Pressure", "TFMP_PRES", TFMP_PRES },
  { "Thin Film Multiphase Saturation", "TFMP_SAT", TFMP_SAT },
  { "Residence Time Function", "RST", RESTIME },
};

int Num_Exo_Var_Names = sizeof(Exo_Var_Names) / sizeof(struct Equation_Names);  
struct Equation_Names Post_Var_Name[] = 
{
  {"STREAM", "stream", -1 },
  {"STREAM_NORMAL_STRESS", "STREAM_NORMAL_STRESS", -1},
  {"STREAM_SHEAR_STRESS", "STREAM_SHEAR_STRESS", -1},
  {"CROSS_STREAM_SHEAR", "CROSS_STREAM_SHEAR", -1},
  {"MEAN_SHEAR", "MEAN_SHEAR", -1},
  {"PRESSURE_CONT", "PRESSURE_CONT", -1},
  {"FILL_CONT", "FILL_CONT", -1},
  {"FIRST_INVAR_STRAIN", "FIRST_INVAR_STRAIN", -1},
  {"SEC_INVAR_STRAIN", "SEC_INVAR_STRAIN", -1},
  {"THIRD_INVAR_STRAIN", "THIRD_INVAR_STRAIN", -1},
  {"DIV_VELOCITY", "DIV_VELOCITY", -1},
  {"VISCOSITY", "VISCOSITY", -1},
  {"DENSITY", "DENSITY", -1},
  {"NS_RESIDUALS", "NS_RESIDUALS", -1},
  {"MM_RESIDUALS", "MM_RESIDUALS", -1},
  {"DIFFUSION_VECTORS", "DIFFUSION_VECTORS", -1},
  {"FLUXLINES", "FLUXLINES", -1},
  {"CONDUCTION_VECTORS", "CONDUCTION_VECTORS", -1},
  {"ENERGY_FLUXLINES", "ENERGY_FLUXLINES", -1},
  {"TIME_DERIVATIVES", "TIME_DERIVATIVES", -1},
  {"STRESS_TENSOR", "STRESS_TENSOR", -1},
  {"REAL_STRESS_TENSOR", "REAL_STRESS_TENSOR", -1},
  {"STRAIN_TENSOR", "STRAIN_TENSOR", -1},
  {"LAGRANGE_CONVECTION", "LAGRANGE_CONVECTION", -1},
  {"SURFACE_VECTORS", "SURFACE_VECTORS", -1},
  {"SHELL_NORMALS", "SHELL_NORMALS", -1},
  {"ERROR_ZZ_VEL", "ERROR_ZZ_VEL", -1},
  {"ERROR_ZZ_Q", "ERROR_ZZ_Q", -1},
  {"ERROR_ZZ_P", "ERROR_ZZ_P", -1},
  {"USER_POST", "USER_POST", -1},
  {"THETA", "THETA", -1},
  {"SPEED", "SPEED", -1},
  {"AC_PRES", "AC_PRES", -1},
  {"TIMESTEPSIZE", "TIMESTEPSIZE", -1},
  {"WALLCLOCKTIME", "WALLCLOCKTIME", -1},
  {"CPUTIME", "CPUTIME", -1},
  {"LIGHT_COMP", "LIGHT_COMP", -1},
  {"UNTRACKED_SPEC", "UNTRACKED_SPEC", -1},
  {"EXTERNAL_FIELD", "EXTERNAL_FIELD", -1},
  {"SH_DIV_S_V_CONT", "SH_DIV_S_V_CONT", -1},
  {"VON_MISES_STRAIN", "VON_MISES_STRAIN", -1},
  {"VON_MISES_STRESS", "VON_MISES_STRESS", -1},
  {"LUB_HEIGHT", "LUB_HEIGHT", -1},
  {"VELO_SPEED", "VELO_SPEED", -1},
  {"GIES_CRIT", "GIES_CRIT", -1},
};

int Num_Post_Var_Names = sizeof(Post_Var_Name) / sizeof(struct Equation_Names);  

struct Equation_Names Var_Units[] = 
{
  { "Velocity, x component",           "[1]",     VELOCITY1 } ,
  { "Velocity, y component",           "[1]",     VELOCITY2 } ,
  { "Velocity, z component",           "[1]",     VELOCITY3 } ,
  { "Temperature",                     "[1]",     TEMPERATURE } ,
  { "Species concentration",           "[1]",     MASS_FRACTION } ,
  { "Mesh displacement, x component",  "[1]",     MESH_DISPLACEMENT1 } ,
  { "Mesh displacement, y component",  "[1]",     MESH_DISPLACEMENT2 } ,
  { "Mesh displacement, z component",  "[1]",     MESH_DISPLACEMENT3 } ,
  { "SURFACE",                         "[1]",     SURFACE } ,
  { "Pressure",                        "[1]",     PRESSURE } ,

  { "Polymer stress, xx component",    "[1]",     POLYMER_STRESS11 } ,
  { "Polymer stress, xy component",    "[1]",     POLYMER_STRESS12 } ,
  { "Polymer stress, xz component",    "[1]",     POLYMER_STRESS13 } ,
  { "Polymer stress, yy component",    "[1]",     POLYMER_STRESS22 } ,
  { "Polymer stress, yz component",    "[1]",     POLYMER_STRESS23 } ,
  { "Polymer stress, zz component",    "[1]",     POLYMER_STRESS33 } ,

  { "SOLID displacement, x component", "[1]",    SOLID_DISPLACEMENT1 } ,
  { "SOLID displacement, y component", "[1]",    SOLID_DISPLACEMENT2 } ,
  { "SOLID displacement, z component", "[1]",    SOLID_DISPLACEMENT3 } ,
		    
  { "Velocity Gradient, xx component", "[1]",     VELOCITY_GRADIENT11 } ,
  { "Velocity Gradient, xy component", "[1]",     VELOCITY_GRADIENT12 } ,
  { "Velocity Gradient, xz component", "[1]",     VELOCITY_GRADIENT13 } ,
  { "Velocity Gradient, yx component", "[1]",     VELOCITY_GRADIENT21 } ,
  { "Velocity Gradient, yy component", "[1]",     VELOCITY_GRADIENT22 } ,
  { "Velocity Gradient, yz component", "[1]",     VELOCITY_GRADIENT23 } ,
  { "Velocity Gradient, zx component", "[1]",     VELOCITY_GRADIENT31 } ,
  { "Velocity Gradient, zy component", "[1]",     VELOCITY_GRADIENT32 } ,
  { "Velocity Gradient, zz component", "[1]",     VELOCITY_GRADIENT33 } ,

  { "Voltage",                         "[1]",     VOLTAGE } ,
  { "Fill",                            "[1]",     FILL } ,
  { "Shear Rate",                      "[1]",     SHEAR_RATE },
  { "Particle velocity, x-comp",        "[1]",      PVELOCITY1 },
  { "Particle velocity, y-comp",        "[1]",      PVELOCITY2 },
  { "Particle velocity, z-comp",        "[1]",     PVELOCITY3 },

  { "Polymer stress component 11 mode 1", "[1]",   POLYMER_STRESS11_1 } ,
  { "Polymer stress component 12 mode 1", "[1]",   POLYMER_STRESS12_1 } ,
  { "Polymer stress component 22 mode 1", "[1]",   POLYMER_STRESS22_1 } ,
  { "Polymer stress component 13 mode 1", "[1]",   POLYMER_STRESS13_1 } ,
  { "Polymer stress component 23 mode 1", "[1]",   POLYMER_STRESS23_1 } ,
  { "Polymer stress component 33 mode 1", "[1]",   POLYMER_STRESS33_1 } ,

  { "Polymer stress component 11 mode 2", "[1]",   POLYMER_STRESS11_2 } ,
  { "Polymer stress component 12 mode 2", "[1]",   POLYMER_STRESS12_2 } ,
  { "Polymer stress component 22 mode 2", "[1]",   POLYMER_STRESS22_2 } ,
  { "Polymer stress component 13 mode 2", "[1]",   POLYMER_STRESS13_2 } ,
  { "Polymer stress component 23 mode 2", "[1]",   POLYMER_STRESS23_2 } ,
  { "Polymer stress component 33 mode 2", "[1]",   POLYMER_STRESS33_2 } ,

  { "Polymer stress component 11 mode 3", "[1]",   POLYMER_STRESS11_3 } ,
  { "Polymer stress component 12 mode 3", "[1]",   POLYMER_STRESS12_3 } ,
  { "Polymer stress component 22 mode 3", "[1]",   POLYMER_STRESS22_3 } ,
  { "Polymer stress component 13 mode 3", "[1]",   POLYMER_STRESS13_3 } ,
  { "Polymer stress component 23 mode 3", "[1]",   POLYMER_STRESS23_3 } ,
  { "Polymer stress component 33 mode 3", "[1]",   POLYMER_STRESS33_3 } ,

  { "Polymer stress component 11 mode 4", "[1]",   POLYMER_STRESS11_4 } ,
  { "Polymer stress component 12 mode 4", "[1]",   POLYMER_STRESS12_4 } ,
  { "Polymer stress component 22 mode 4", "[1]",   POLYMER_STRESS22_4 } ,
  { "Polymer stress component 13 mode 4", "[1]",   POLYMER_STRESS13_4 } ,
  { "Polymer stress component 23 mode 4", "[1]",   POLYMER_STRESS23_4 } ,
  { "Polymer stress component 33 mode 4", "[1]",   POLYMER_STRESS33_4 } ,

  { "Polymer stress component 11 mode 5", "[1]",   POLYMER_STRESS11_5 } ,
  { "Polymer stress component 12 mode 5", "[1]",   POLYMER_STRESS12_5 } ,
  { "Polymer stress component 22 mode 5", "[1]",   POLYMER_STRESS22_5 } ,
  { "Polymer stress component 13 mode 5", "[1]",   POLYMER_STRESS13_5 } ,
  { "Polymer stress component 23 mode 5", "[1]",   POLYMER_STRESS23_5 } ,
  { "Polymer stress component 33 mode 5", "[1]",   POLYMER_STRESS33_5 } ,

  { "Polymer stress component 11 mode 6", "[1]",   POLYMER_STRESS11_6 } ,
  { "Polymer stress component 12 mode 6", "[1]",   POLYMER_STRESS12_6 } ,
  { "Polymer stress component 22 mode 6", "[1]",   POLYMER_STRESS22_6 } ,
  { "Polymer stress component 13 mode 6", "[1]",   POLYMER_STRESS13_6 } ,
  { "Polymer stress component 23 mode 6", "[1]",   POLYMER_STRESS23_6 } ,
  { "Polymer stress component 33 mode 6", "[1]",   POLYMER_STRESS33_6 } ,

  { "Polymer stress component 11 mode 7", "[1]",   POLYMER_STRESS11_7 } ,
  { "Polymer stress component 12 mode 7", "[1]",   POLYMER_STRESS12_7 } ,
  { "Polymer stress component 22 mode 7", "[1]",   POLYMER_STRESS22_7 } ,
  { "Polymer stress component 13 mode 7", "[1]",   POLYMER_STRESS13_7 } ,
  { "Polymer stress component 23 mode 7", "[1]",   POLYMER_STRESS23_7 } ,
  { "Polymer stress component 33 mode 7", "[1]",   POLYMER_STRESS33_7 } ,   /*76*/ 

  { "species concentration comp 0",   "[1]",     SPECIES_UNK_0  } ,  
  { "species concentration comp 1",   "[1]",     SPECIES_UNK_1  } ,  
  { "species concentration comp 2",   "[1]",     SPECIES_UNK_2  } ,  
  { "species concentration comp 3",   "[1]",     SPECIES_UNK_3  } ,  
  { "species concentration comp 4",   "[1]",     SPECIES_UNK_4  } ,  
  { "species concentration comp 5",   "[1]",     SPECIES_UNK_5  } , 
  { "species concentration comp 6",   "[1]",     SPECIES_UNK_6  } ,  
  { "species concentration comp 7",   "[1]",     SPECIES_UNK_7  } ,  
  { "species concentration comp 8",   "[1]",     SPECIES_UNK_8  } ,  
  { "species concentration comp 9",   "[1]",     SPECIES_UNK_9  } , 
  { "species concentration comp 10",  "[1]",     SPECIES_UNK_10 } , 
  { "species concentration comp 11",  "[1]",     SPECIES_UNK_11 } , 
  { "species concentration comp 12",  "[1]",     SPECIES_UNK_12 } ,  
  { "species concentration comp 13",  "[1]",     SPECIES_UNK_13 } , 
  { "species concentration comp 14",  "[1]",     SPECIES_UNK_14 } ,  
  { "species concentration comp 15",  "[1]",     SPECIES_UNK_15 } ,  
  { "species concentration comp 16",  "[1]",     SPECIES_UNK_16 } , 
  { "species concentration comp 17",  "[1]",     SPECIES_UNK_17 } ,  
  { "species concentration comp 18",  "[1]",     SPECIES_UNK_18 } ,  
  { "species concentration comp 19",  "[1]",     SPECIES_UNK_19 } ,  
  { "species concentration comp 20",  "[1]",     SPECIES_UNK_20 } ,  
  { "species concentration comp 21",  "[1]",     SPECIES_UNK_21 } ,  
  { "species concentration comp 22",  "[1]",     SPECIES_UNK_22 } ,  
  { "species concentration comp 23",  "[1]",     SPECIES_UNK_23 } , 
  { "species concentration comp 24",  "[1]",     SPECIES_UNK_24 } ,  
  { "species concentration comp 25",  "[1]",     SPECIES_UNK_25 } ,  
  { "species concentration comp 26",  "[1]",     SPECIES_UNK_26 } ,  
  { "species concentration comp 27",  "[1]",     SPECIES_UNK_27 } ,  
  { "species concentration comp 28",  "[1]",     SPECIES_UNK_28 } ,  
  { "species concentration comp 29",  "[1]",     SPECIES_UNK_29 } ,  
  
  { "VolFracPh_0",  "[1]",   VOLF_PHASE_0   } ,  
  { "VolFracPh_1",  "[1]",   VOLF_PHASE_1   } , 
  { "VolFracPh_2",  "[1]",   VOLF_PHASE_2   } , 
  { "VolFracPh_3",  "[1]",   VOLF_PHASE_3   } ,  
  { "VolFracPh_4",  "[1]",   VOLF_PHASE_4   } ,  
    
  { "Liquid phase pressue, porous",   "[1]",   POR_LIQ_PRES   } ,  
  { "Gas phase pressure, porous",   "[1]",   POR_GAS_PRES   } ,           
  { "Porosity of medium, porous",    "[1]",   POR_POROSITY   } ,  
  { "Temperature, porous",    "[1]",   POR_TEMP   } ,  
  { "Saturation of medium, porous",  "[1]",   POR_SATURATION } ,  /*116*/

  { "Sink mass accumulation, porous", "[1]", POR_SINK_MASS }, 

  { "Vorticity direction, x component", "[1]",  VORT_DIR1 } ,
  { "Vorticity direction, y component", "[1]",  VORT_DIR2 } ,
  { "Vorticity direction, z component", "[1]",  VORT_DIR3 } ,
  { "Vorticity eigenvalue",             "[1]",  VORT_LAMBDA } ,

  { "Curvature of level set field", "H", CURVATURE },

  { "Volume Lagrange Multiplier, x component", "[1]",  LAGR_MULT1 } ,
  { "Volume Lagrange Multiplier, y component", "[1]",  LAGR_MULT2 } ,
  { "Volume Lagrange Multiplier, z component", "[1]",  LAGR_MULT3 } ,

  { "Bond Structure", "[1]", BOND_EVOLUTION },

  { "Surface charge density",           "[1]", SURF_CHARGE},
 
  { "Extension velocity, normal component", "[1]", EXT_VELOCITY},

  { "Electric field, x component", "[1]", EFIELD1 },           /*129*/
  { "Electric field, y component", "[1]", EFIELD2 },
  { "Electric field, z component", "[1]", EFIELD3 },

  { "Norm of potential field", "[1]", ENORM },

  { "Normal to LS function, x component", "[1]", NORMAL1 },
  { "Normal to LS function, y component", "[1]", NORMAL2 },
  { "Normal to LS function, z component", "[1]", NORMAL3 },

  { "Shell element curvature", "[1]", SHELL_CURVATURE},
  { "Shell element second curvature", "[1]", SHELL_CURVATURE2},
  { "Shell element tension", "[1]", SHELL_TENSION},
  { "Shell element x coordinate", "[1]", SHELL_X},
  { "Shell element t coordinate", "[1]", SHELL_Y},
  { "Shell element user", "[1]", SHELL_USER},

  { "First  phase function", "[1]", PHASE1},
  { "Second phase function", "[1]", PHASE2},
  { "Third  phase function", "[1]", PHASE3},
  { "Fourth phase function", "[1]", PHASE4},
  { "Fifth  phase function", "[1]", PHASE5},
  
  { "Shell element orientation angle, theta", "[1]", SHELL_ANGLE1},
  { "Shell element orientation angle, phi", "[1]", SHELL_ANGLE2},

  { "(I-NN).del.v", "[1]", SHELL_SURF_DIV_V},
  { "Mean curvature of surface", "[1]", SHELL_SURF_CURV},
  { "normal component of vorticity", "[1]", N_DOT_CURL_V},
  { "Normal x component of velocity gradient_s tensor", "[1]", GRAD_S_V_DOT_N1},
  { "Normal y component of velocity gradient_s tensor", "[1]", GRAD_S_V_DOT_N2},
  { "Normal z component of velocity gradient_s tensor", "[1]", GRAD_S_V_DOT_N3},
  { "Acoustic harmonic pressure - real part", "[1]", ACOUS_PREAL},
  { "Acoustic harmonic pressure - imag part", "[1]", ACOUS_PIMAG},
  { "Shell surface diffusion flux", "[1]", SHELL_DIFF_FLUX},
  { "Shell surface curvature", "[1]", SHELL_DIFF_CURVATURE},
  { "Shell normal vector, X component", "[1]", SHELL_NORMAL1},
  { "Shell normal vector, Y component", "[1]", SHELL_NORMAL2},
  { "Shell normal vector, Z component", "[1]", SHELL_NORMAL3},
  { "Acoustic Reynolds Stress", "[1]", ACOUS_REYN_STRESS},
  { "Acoustic Boundary Velocity", "[1]", SHELL_BDYVELO},
  { "Shell Lubrication Pressure", "[1]", SHELL_LUBP}, 
  { "Lubrication Pressure", "[1]", LUBP},
  { "Shell Film Pressure", "[1]", SHELL_FILMP},
  { "Shell Film Thickness", "[1]", SHELL_FILMH},
  { "Shell Particles Concentration", "[1]", SHELL_PARTC},
  { "Structured Porous Shell Saturation - Closed Cells", "[1]", SHELL_SAT_CLOSED},
  { "Structured Porous Shell Pressure - Open Cells", "[1]", SHELL_PRESS_OPEN},
  { "Shell temperature", "[1]", SHELL_TEMPERATURE},
  { "Shell delta H", "[1]", SHELL_DELTAH},
  { "Shell Lubrication Curvature", "[1]", SHELL_LUB_CURV},
  { "Structured Porous Shell Saturation - Gas Compression", "[1]", SHELL_SAT_GASN},
  { "Top wall shear rate", "[1]", SHELL_SHEAR_TOP},
  { "Bottom wall shear rate", "[1]", SHELL_SHEAR_BOT},
  { "Cross stream shear stress", "[1]", SHELL_CROSS_SHEAR},
  { "Maximum Von Mises strain", "[1]", MAX_STRAIN},
  { "Von Mises strain", "[1]", CUR_STRAIN},
  { "Lubrication Pressure 2", "[1]", LUBP_2},
  { "Structured Porous Shell Pressure - Open Cells 2", "[1]", SHELL_PRESS_OPEN_2},
  { "Shell Lubrication Curvature 2", "[1]", SHELL_LUB_CURV_2},
  { "Plus Intensity", "[1]", LIGHT_INTP},
  { "Minus Intensity", "[1]", LIGHT_INTM},
  { "Dispersive Intensity", "[1]", LIGHT_INTD},
  { "Thin Film Multiphase Lubrication Pressure", "[1]", TFMP_PRES},
  { "Thin Film Multiphase Saturation", "[1]", TFMP_SAT},
  { "Residence Time Function", "[1]", RESTIME},  
};

int Num_Var_Units = sizeof(Var_Units) / sizeof(struct Equation_Names);  

/*
 * Array of pointers is more convenient access to dynamically-allocated
 * BC Descriptions than otherwise requiring a search through BC_Types...
 */

struct BC_descriptions  **new_BC_Desc;

/*
 * This variable counts up any dynamically allocated BC_Description structures
 * so that we can send them to other processors.
 */

int num_new_BC_Desc;

/*
 * This variable is similar in function num_new_BC_Desc in that it counts
 * up the dynamically allocated table structures associated with BCS
 * so they can be sent to other processors
 */

int num_BC_Tables;

struct Data_Table *BC_Tables[MAX_BC_TABLES];
				
int num_MP_Tables;

struct Data_Table *MP_Tables[MAX_MP_TABLES];

int num_ext_Tables;

struct Data_Table *ext_Tables[MAX_EXT_TABLES];
				
int num_AC_Tables;

struct Data_Table *AC_Tables[MAX_AC_TABLES];
				
#endif
