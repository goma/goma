/************************************************************************ *
* Goma - Multiphysics finite element software                             *
* Sandia National Laboratories                                            *
*                                                                         *
* Copyright (c) 2015 Sandia Corporation.                                  *
*                                                                         *
* Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,  *
* the U.S. Government retains certain rights in this software.            *
*                                                                         *
* This software is distributed under the GNU General Public License.      *
\************************************************************************/
 

#ifndef _MM_FILL_TERMS_SEGREGATED_H
#define _MM_FILL_TERMS_SEGREGATED_H
#endif

#ifdef EXTERN
#undef EXTERN
#endif

#ifdef _MM_FILL_TERMS_SEGREGATED_C
#define EXTERN /* do nothing */
#endif

#ifndef _MM_FILL_TERMS_SEGREGATED_C
#define EXTERN extern
#endif

EXTERN void
load_splitb_esp(int ielem, Exo_DB *exo);

EXTERN void
load_splitb_fv(int ielem);

EXTERN void
load_splitb_fv_grads(int ielem);

EXTERN int assemble_aux_u
(double,   // Current time
 double);  // Current time step


EXTERN int assemble_press_poisson
(double,   // Current time
 double);  // Current time step


EXTERN int assemble_press_proj
(double,   // Current time
 double);  // Current time step

EXTERN int assemble_press_update
(double,   // Current time
 double);  // Current time step
