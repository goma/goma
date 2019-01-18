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
 * $Id: mm_dil_viscosity.h,v 5.3 2008-03-13 01:12:31 hkmoffa Exp $
 */

#ifndef _MM_DIL_VISCOSITY_H
#define _MM_DIL_VISCOSITY_H

#ifdef EXTERN
#undef EXTERN
#endif
#ifdef _MM_DIL_VISCOSITY_C
#define EXTERN
#else
#define EXTERN extern
#endif


EXTERN double dil_viscosity		        /* mm_dil_viscosity.c                */
(struct Generalized_Newtonian *gn_local,  /* gn_local                          */
       dbl gamma[DIM][DIM],	        	/* gamma - strain rate tensor    */
       const dbl mu,
       const VISCOSITY_DEPENDENCE_STRUCT *d_mu, /* d_mu - viscosity dependence       */
       DILVISCOSITY_DEPENDENCE_STRUCT *d_dilMu);/* d_dilMu - dil_viscosity dependence */

#endif /* _MM_DIL_VISCOSITY_H */
