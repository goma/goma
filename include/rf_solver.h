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
 *$Id: rf_solver.h,v 5.1 2007-09-18 18:53:47 prschun Exp $
 */

#ifndef _RF_SOLVER_H
#define _RF_SOLVER_H

#include "rf_io_const.h"	/* for definition of MAX_CHAR_IN_INPUT */
#include "rf_solver_const.h"	/* for kinds of direct solvers */

/*
 * These strings are read from the input deck and processed later as needed
 * on their way to be sacrificed to Aztec. The meanings of these strings
 * follow the Aztec manual description on pages 2-6, SAND 95-1559.
 */

#ifndef String_line_is_defined
#define String_line_is_defined
typedef char String_line [MAX_CHAR_IN_INPUT];
#endif

String_line Matrix_Solver;

String_line Matrix_Format;

String_line Matrix_Scaling;

String_line Matrix_Preconditioner;

String_line Matrix_Subdomain_Solver; /* new Aztec 2.x option */

String_line Matrix_Residual_Norm_Type;

String_line Matrix_Output_Type;

String_line Matrix_Factorization_Reuse;

String_line Matrix_Graph_Fillin; /* new Aztec 2.x option */

String_line Matrix_Maximum_Iterations;

String_line Matrix_Polynomial_Order;

String_line Matrix_Factor_Overlap;

String_line Matrix_Overlap_Type;

String_line Matrix_Krylov_Subspace;

String_line Matrix_Orthogonalization;

String_line Matrix_Auxiliary_Vector;

String_line Matrix_Convergence_Tolerance;

String_line Matrix_Drop_Tolerance;

String_line Matrix_Factorization_Save; /* Aztec 2 */

String_line Matrix_ILUT_Fill_Factor; /* Aztec 2 */

String_line Matrix_RILU_Relax_Factor; /* Aztec 2 */

String_line Matrix_BILU_Threshold; /* Trilinos 1 */

String_line Matrix_Relative_Threshold; /* Trilinos 2 */

String_line Matrix_Absolute_Threshold; /* Trilinos 2 */

String_line Amesos_Package;

String_line AztecOO_Solver;

/*
 * A new Aztec 2.0 option. There are more and difft options and our
 * previous options probably ought to be revised to reflect the newer
 * Aztec 2.0 capability specifications.
 */

String_line Matrix_Reorder;

int Linear_Solver;		/* Aztec, Sparse, MA28, UMFPACK */

int UMFPACK_IDIM;		/* UMFPACK STORAGE CONSTANT */
int UMFPACK_XDIM;		/* UMFPACK STORAGE CONSTANT */
int LOCA_UMF_ID;                /* UMFPACK SYSTEM ID */

int Max_Newton_Steps; 	/* Maximum number of Newton steps to take.     */
int Guess_Flag;		/* Indicates the type of initial guess         */

double damp_factor;
double damp_factor1;		/* Relaxation factor for Newton iteration */
				/* damp_factor1 = 1.0 is full Newton */
				/* damp_factor1 = 0.0 is not updating our */
				/*                   solution estimate */
double damp_factor2,            /* Additional damping factors for custom */
       damp_factor3,            /* schemes for automatic control with    */
       custom_tol1,             /* NORM(0,0) tolerances                  */
       custom_tol2,
       custom_tol3;
double var_damp[MAX_VARIABLE_TYPES];  /* variable specific damp factors */

int Newt_Jacobian_Reformation_stride;/*Stride for reformation of jacobian for
                                  modified newton scheme               */
int Time_Jacobian_Reformation_stride; 
int modified_newton;            /*boolean flag for modified Newton */
int save_old_A;                 /*boolean flag for saving old A matrix
				 for resolve reasons with AZTEC.   There
				 are at least four reasons, that you
				 can see in sl_util.c */
double convergence_rate_tolerance; /* tolerance for jacobian reformation
                                       based on convergence rate */
double modified_newt_norm_tol; /* tolerance for jacobian reformation 
                                       based on residual norm */

double Epsilon[3];	/* Used for determining stopping criteria.     */

int NZeros;             /* Number of nonzeros in this procs matrix     */


int nnz_own;			/* number of nonzeroes in the matrix that
				 * belong to this processor
				 */

int nnz_total;			/* total number of nonzeroes in the matrix
				 * that this processor sees, including the
				 * entries for external unknowns that will
				 * not be updated by this processor  */

int GNZeros;            /* Number of nonzeros in global matrix         */

int fill_zeros;         /* number of nonzeros in fill matrix for this
			  processor */

int Gfill_zeros;        /* number of nonzeros in fill matrix for the
			  global problem */

int PSPG;               /* 1 means pressure stabilized Petrov-Galerkin is used */
int PSPP;               /* 1 means pressure stabilized polynomial projection is used */
double PS_scaling;      /* This term is a constant scaling for the PSPG or PSPP term */

int Filter_Species, filter_species_material_number;
double c_min, c_max;

int Include_Visc_Sens, Visc_Sens_Copy;
				 /* 1 means to include the sensitivities of the 
				viscosity functions in the jacobian matrix.  
				0 indicates that these sensitivities are not to 
				be included.  The latter case is useful for highly
				 shear thinning viscosity models.  By disabling 
				the viscosity sensititives convergence can be 
				achieved albeit at a less than quadratic rate 
		       		*/

#endif
