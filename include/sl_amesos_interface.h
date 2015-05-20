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
 
#ifndef _SL_AMESOS_INTERFACE_H
#define _SL_AMESOS_INTERFACE_H
#endif

#ifdef EXTERN
#undef EXTERN
#endif

#ifdef _SL_AMESOS_INTERFACE_CC
#define EXTERN extern "C"
#else
#define EXTERN extern
#endif


#ifdef ENABLE_AMESOS

/* Use this version when Amesos is linked in through Trilinos */
#ifdef TRILINOS
EXTERN void
amesos_solve_msr ( char *,
		   struct Aztec_Linear_Solver_System *,
		   double *,
		   double *,
		   int ,
		   int );
EXTERN int
amesos_solve_epetra( char *choice,
                     struct Aztec_Linear_Solver_System *ams,
		     double *x_,
		     double *resid_vector,
		     int imtrx );

EXTERN void
trilinos_solve_ls (double *,
		   int *,
		   int *,
		   double *,
		   double *,
		   double *,
		   int, 
		   int, 
		   int,
		   int);
#endif

EXTERN void *
construct_Epetra_CrsMatrix ( struct Aztec_Linear_Solver_System * );

EXTERN void *
construct_Epetra_CrsMatrix_ls ( int, int, int );



#endif
