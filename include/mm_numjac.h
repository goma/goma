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
 * mm_numjac.h -- prototype declarations for mm_numjac.c
 */

#ifndef GOMA_MM_NUMJAC_H
#define GOMA_MM_NUMJAC_H

#include "dpi.h"
#include "el_elm.h"
#include "exo_struct.h"
#include "mm_more_utils.h"

struct Aztec_Linear_Solver_System;
struct elem_side_bc_struct;
#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_MM_NUMJAC_C
#define EXTERN
#else
#define EXTERN extern
#endif

#define ELEM_LIST_SIZE (MAX_SUR_ELEM_3D*MAX_SUR_ELEM_3D)

/* Good for Debug_Flag == -2: 1.0e-3, 1.0e-4, 1.0e-2, 1.0e-8, 1.0e-4, 1.0e-12 */
/* Good for Debug_Flag == -3: 1.0e-3, 1.0e-3, 1.0e-2, 1.0e-4, 1.0e-3, 1.0e-8 */
/*
 * The relative differences used for obtaining the most accurate jacobians
 * is normally set to the square root of machine precision, according to
 * most text books on the subject. This is roughly 1.0E-7.
 */
#define DELTA_UNKNOWN (1.0e-7)	           /* for infinitesimal finite difference */
//#define FD_DELTA_UNKNOWN (1.0e-4)	   /* for finite difference */
#define FD_DELTA_UNKNOWN (1.0e-6)	   /* for finite difference */

#define RESIDUAL_TOLERANCE (1.0e-3)        /* absolute error required for output */
#define SCALED_RESIDUAL_TOLERANCE (1.0e-2) /* scaled error required for output */
#define MIXED_RESIDUAL_TOLERANCE (1.0e-4)  /* If the scaled error is too large, but the
					    * aboslute error is smaller than this,
					    * then don't report.  This should never be
					    * larger than RESIDUAL_TOLERANCE. */
#define MIXED_SCALED_RESIDUAL_TOLERANCE (1.0e-3) /* If the aboslute error is too large,
						  * but the scaled error is smaller than
						  * this, then don't report.  This should
						  * never be larger than
						  * SCALED_RESIDUAL_TOLERANCE. */
#define SCALED_RESIDUAL_TOLERANCE_CUTOFF (1.0e-8) /* don't report scaled error if values
						    * are this small */
EXTERN void
numerical_jacobian_compute_stress(struct Aztec_Linear_Solver_System *ams,
		   double x[],	/* Solution vector for the current processor */
		   double resid_vector[],   /* Residual vector for the current
					     * processor */
		   double delta_t, /* time step size */
		   double theta, /* parameter to vary time integration from
				    explicit (theta = 1) to
				    implicit (theta = 0) */
		   double x_old[], /* Value of the old solution vector */
		   double x_older[], /* Value of the real old soln vect */

		   double xdot[], /* Value of xdot predicted for new solution */
		   double xdot_old[], /* Value of xdot at previous time */

		   double x_update[],
		   int num_total_nodes,

		   struct elem_side_bc_struct *first_elem_side_BC_array[],
				/* This is an array of pointers to the first
				   surface integral defined for each element.
				   It has a length equal to the total number
				   of elements defined on the current proc */
		   int Debug_Flag, /* flag for calculating numerical jacobian
				      -1 == calc num jac w/o rescaling
				      -2 == calc num jac w/  rescaling */
		   double time_value, /* current value of time */
		   Exo_DB *exo,	    /* ptr to whole fe mesh */
		   Dpi *dpi,        /* any distributed processing info */
		   double *h_elem_avg,
  double *U_norm);


EXTERN void numerical_jacobian	/* mm_numjac.c                               */
(struct Aztec_Linear_Solver_System *, /* ams                           */
       double [],		/* x - soln vector for current processor     */
       double [],		/* resid_vector -for current processor       */
       double ,			/* delta_t - time step size                  */
       double ,			/* theta - parameter varies time integration *
				 * from explicit (theta = 1) to              *
				 * implicit (theta = 0)                      */
       double [],		/* x_old - Value of the old solution vector  */
       double [],		/* x_older - Value of the real old soln vect */
       double [],		/* xdot - predicted for new solution         */
       double [],		/* xdot_old - Value of xdot at previous time */
       double [],		/* x_update                                  */
       int ,			/* num_total_nodes                           */
       struct elem_side_bc_struct *[], /* first_elem_side_BC_array           *
					* This is an array of pointers to    *
					* the first surface integral defined *
					* for each element.  It has a length *
					* equal to the total number of       *
					* elements on the current proc       */
       int ,			/* Debug_Flag - flag for calculating         *
				 * numerical jacobian:                       *
				 *    -1 == calc num jac w/o rescaling       *
				 *    -2 == calc num jac w/  rescaling       *
                                 *    -3 == calc num jac w/ diagonal scaling */
       double ,			/* time_value - current time                 */
       Exo_DB *,		/* exo - ptr to whole fe mesh                */
       Dpi *,			/* dpi - any distributed processing info     */
       double *,
       double *);

#ifndef COUPLED_FILL
EXTERN void numerical_jacobian_fill /* mm_numjac.c                           */
(int [],			/* ijaf - column pointers into fill matrix   */
       double [],		/* afill - non-zero entries in fill matrix   */
       double [],		/* xf - fill Solution vector                 */
       double [],		/* rf - Residual vector for fill eqns        */
       double ,			/* delta_t - time step size                  */
       double ,			/* theta - parameter varies time integration *
				 * from explicit (theta = 1) to              *
				 * implicit (theta = 0)                      */
       double [],		/* x - current big soln vector (everything)  */
       double [],		/* x_old  - old solution vector              */
       double [],		/* xdot - predicted for new solution         */
       int ,			/* Debug_Flag - flag for calculating         *
				 * numerical jacobian  -4 == calc num jac w/ *
				 * rescaling                                 */
       int [],			/* node_to_fill - this is a map from the     */
       Exo_DB *,		/* exo - ptr to whole fe mesh                */
       Dpi *);			/* dpi - ptr to parallel info                */
#endif /* not COUPLED_FILL */
extern double calc_numerical_delta(double);
extern void   AF_assemble_Residual_Only(void);
extern void   AF_restore_Jacobian_Flag(void);

#endif /* GOMA_MM_NUMJAC_H */
