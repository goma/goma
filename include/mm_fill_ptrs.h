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
 
#ifndef _MM_FILL_PTRS_H
#define _MM_FILL_PTRS_H

#ifdef EXTERN
#undef EXTERN
#endif

#ifdef _MM_FILL_PTRS_C
#define EXTERN /* do nothing */
#endif

#ifndef _MM_FILL_PTRS_C
#define EXTERN extern
#endif


/*
 * Externals for the solution vector and its derivative at the
 * current and old times. These are filled in every time an element
 * dofptrs are set up.
 */
extern double *x_static;
extern double *x_old_static;
extern double *xdot_static;
extern double *xdot_old_static;
extern double *x_dbl_dot_static;
extern double *x_dbl_dot_old_static;
extern double *x_pred_static;

EXTERN int load_ei              /* mm_fill_ptrs.c                            */
PROTO((const int ,		/* element index                             */
       const Exo_DB *,	        /* exo data base                             */
       struct Element_Indices *,
       int)); /* Pointer to the structure to be filled
	                         * up. If 0, this is the globa ei and we are
				 * in master mode. If 1, this is a related ei
				 * and we are in slave mode                  */

EXTERN int load_elem_dofptr	/* mm_fill_ptrs.c                            */
PROTO((const int,		/* element index                             */
       const Exo_DB *,		/* Exodus database pointer                   */
       dbl * const,		/* x - unknown vector                        */
       dbl * const,		/* x_old - unknown vector, previous timestep */
       dbl * const,		/* xdot - difference approx to dx/dt         */
       dbl * const,             /* xdot_old - time derivative at t = t_n     */
       dbl * const,		/* resid_vector - residual RHS               */
       const int ));	       	/* if early_return is true we just get 
				 * some of the element information for 
				 * routines that need it like 
				 * global_h_elem_siz, but not all the 
				 * computationally expensive pointers        */

int 
load_elem_dofptr_all(const int ielem,
                     const Exo_DB * exo);

EXTERN int load_elem_aijaptr	/* mm_fill_ptrs.c                            */
PROTO((int [],			/* ija - column indeces                      */
       dbl []));		/* a - matrix nonzero values                 */

#endif /* _MM_FILL_PTRS_H */
