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
 

#ifndef GOMA_SL_UMF_H
#define GOMA_SL_UMF_H

#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_SL_UMF_C
#define EXTERN /* do nothing */
#endif

#ifndef GOMA_SL_UMF_C
#define EXTERN extern
#endif

extern int SL_UMF
( int,
	 int *,
	 int *,
	 int *,
	 int *,
	 int *,
	 int *,
	 int *,
	 double *,
	 double *,
	 double *  );

struct UMF_Linear_Solver_System
{
	int n, nnz;
        int *ap, *ai, *atp, *ati;
	double *ax, *atx;
        void *symbolic, *numeric;
};

#endif /* GOMA_SL_UMF_H */


