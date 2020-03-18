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

#ifndef GOMA_MM_DIL_VISCOSITY_H
#define GOMA_MM_DIL_VISCOSITY_H

#include "mm_as_structs.h"
#include "mm_mp_structs.h"
#include "mm_viscosity.h"
#include "std.h"
#ifdef EXTERN
#undef EXTERN
#endif
#ifdef GOMA_MM_DIL_VISCOSITY_C
#define EXTERN
#else
#define EXTERN extern
#endif


EXTERN double dil_viscosity                        /* mm_dil_viscosity.c                */
    (GEN_NEWT_STRUCT *gn_local, const dbl muValue, const VISCOSITY_DEPENDENCE_STRUCT *d_mu, DILVISCOSITY_DEPENDENCE_STRUCT *d_dilMu);/* d_dilMu - dil_viscosity dependence */

int
ls_modulate_dilviscosity ( double *kappa1,
                           double  kappa2,
                           double width,
                           double pm_minus,
                           double pm_plus,
                           DILVISCOSITY_DEPENDENCE_STRUCT *d_dilMu);

#endif /* GOMA_MM_DIL_VISCOSITY_H */
