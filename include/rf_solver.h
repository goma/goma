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
 *$Id: rf_solver.h,v 5.1 2007-09-18 18:53:47 prschun Exp $
 */

#ifndef GOMA_RF_SOLVER_H
#define GOMA_RF_SOLVER_H

#include "rf_fem_const.h"    /* MAX_NUM_MATRICES */
#include "rf_io_const.h"     /* for definition of MAX_CHAR_IN_INPUT */
#include "rf_solver_const.h" /* for kinds of direct solvers */

/*
 * These strings are read from the input deck and processed later as needed
 * on their way to be sacrificed to Aztec. The meanings of these strings
 * follow the Aztec manual description on pages 2-6, SAND 95-1559.
 */

#ifndef String_line_is_defined
#define String_line_is_defined
typedef char String_line[MAX_CHAR_IN_INPUT];
#endif

extern String_line Matrix_Solver;

extern String_line Matrix_Format;

extern String_line Matrix_Scaling;

extern String_line Matrix_Preconditioner;

extern String_line Matrix_Subdomain_Solver; /* new Aztec 2.x option */

extern String_line Matrix_Residual_Norm_Type;

extern String_line Matrix_Output_Type;

extern String_line Matrix_Factorization_Reuse;

extern String_line Matrix_Graph_Fillin; /* new Aztec 2.x option */

extern String_line Matrix_Maximum_Iterations;

extern String_line Matrix_Polynomial_Order;

extern String_line Matrix_Factor_Overlap;

extern String_line Matrix_Overlap_Type;

extern String_line Matrix_Krylov_Subspace;

extern String_line Matrix_Orthogonalization;

extern String_line Matrix_Auxiliary_Vector;

extern String_line Matrix_Convergence_Tolerance;

extern String_line Matrix_Drop_Tolerance;

extern String_line Matrix_Factorization_Save; /* Aztec 2 */

extern String_line Matrix_ILUT_Fill_Factor; /* Aztec 2 */

extern String_line Matrix_RILU_Relax_Factor; /* Aztec 2 */

extern String_line Matrix_BILU_Threshold; /* Trilinos 1 */

extern String_line Matrix_Relative_Threshold; /* Trilinos 2 */

extern String_line Matrix_Absolute_Threshold; /* Trilinos 2 */

extern String_line Amesos_Package;
extern String_line Amesos2_Package;

extern String_line AztecOO_Solver;

extern String_line Stratimikos_File[MAX_NUM_MATRICES];
extern String_line Amesos2_File[MAX_NUM_MATRICES];

/*
extern  * A new Aztec 2.0 option. There are more and difft options and our
extern  * previous options probably ought to be revised to reflect the newer
extern  * Aztec 2.0 capability specifications.
extern  */

extern String_line Matrix_Reorder;

extern int Linear_Solver; /* Aztec, Sparse, MA28, UMFPACK */

extern int UMFPACK_IDIM; /* UMFPACK STORAGE CONSTANT */
extern int UMFPACK_XDIM; /* UMFPACK STORAGE CONSTANT */
extern int LOCA_UMF_ID;  /* UMFPACK SYSTEM ID */

extern int Max_Newton_Steps;  /* Maximum number of Newton steps to take.     */
extern int Guess_Flag;        /* Indicates the type of initial guess         */
extern int Conformation_Flag; /* Indicates mapping from stress to log-conformation tensor */
extern int Print3DBCDup;

extern double damp_factor;
extern double damp_factor1; /* Relaxation factor for Newton iteration */
                            /* damp_factor1 = 1.0 is full Newton */
                            /* damp_factor1 = 0.0 is not updating our */
                            /*                   solution estimate */
extern double damp_factor2, /* Additional damping factors for custom */
    damp_factor3,           /* schemes for automatic control with    */
    custom_tol1,            /* NORM(0,0) tolerances                  */
    custom_tol2, custom_tol3;
extern double var_damp[MAX_VARIABLE_TYPES]; /* variable specific damp factors */

extern int Newt_Jacobian_Reformation_stride; /*Stride for reformation of jacobian for
                                   modified newton scheme               */
extern int Time_Jacobian_Reformation_stride;
extern int Newton_Line_Search_Type;
extern int modified_newton;               /*boolean flag for modified Newton */
extern int save_old_A;                    /*boolean flag for saving old A matrix
                                    for resolve reasons with AZTEC.   There
                                    are at least four reasons, that you
                                    can see in sl_util.c */
extern double convergence_rate_tolerance; /* tolerance for jacobian reformation
                                       based on convergence rate */
extern double modified_newt_norm_tol;     /* tolerance for jacobian reformation
                                           based on residual norm */

extern double Epsilon[MAX_NUM_MATRICES][3]; /* Used for determining stopping criteria.     */
extern int Solver_Output_Format;            /* Bitmap for Solver output columns    */
extern int Output_Variable_Stats;           /* Toggle for Variable Statistics    */
extern int Output_Variable_Regression;      /* Toggle for Variable Regression    */

extern int NZeros; /* Number of nonzeros in this procs matrix     */

extern int GNZeros; /* Number of nonzeros in global matrix         */

extern int fill_zeros; /* number of nonzeros in fill matrix for this
                  processor */

extern int Gfill_zeros; /* number of nonzeros in fill matrix for the
                   global problem */

extern int PSPG;          /* 1 means pressure stabilized Petrov-Galerkin is used */
extern int PSPP;          /* 1 means pressure stabilized polynomial projection is used */
extern double PS_scaling; /* This term is a constant scaling for the PSPG or PSPP term */
extern int Cont_GLS;      /* 1 means continuity stabilization is used */

extern int Filter_Species, filter_species_material_number;
extern double c_min, c_max;

extern int Include_Visc_Sens, Visc_Sens_Copy, Visc_Sens_Factor;
/* 1 means to include the sensitivities of the
viscosity functions in the jacobian matrix.
0 indicates that these sensitivities are not to
be included.  The latter case is useful for highly
shear thinning viscosity models.  By disabling
the viscosity sensititives convergence can be
achieved albeit at a less than quadratic rate
*/

#endif
