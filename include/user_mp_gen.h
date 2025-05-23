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
 * $Id: user_mp_gen.h,v 5.2 2010-04-05 16:46:05 prschun Exp $
 */

#ifndef GOMA_USER_MP_GEN_H
#define GOMA_USER_MP_GEN_H

#include "el_elm.h"
#include "std.h"
#include "user_mp.h"
#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_USER_MP_GEN_C
#define EXTERN
#
#endif

#ifndef GOMA_USER_MP_GEN_C
#define EXTERN extern
#endif

EXTERN int usr_heat_source_gen(dbl *,              /* h - volumetric heat source                */
                               dbl[MDE],           /* dhdT - temperature dependence.            */
                               dbl[DIM][MDE],      /* dhdX - spatial dependence.                */
                               dbl[DIM][MDE],      /* dhdV - velocity dependence.               */
                               dbl[MAX_CONC][MDE], /* dhdC - concentration dependence.          */
                               dbl[MDE],           /* dhdVolt - voltabe dependence.          */
                               dbl *,              /* param - user-defined parameter list       */
                               dbl);               /* time */

EXTERN int usr_viscosity_gen(dbl *,              /* mu - viscosity                            */
                             dbl[DIM][DIM],      /* gamma_dot - strain rate tensor            */
                             dbl *,              /* d_mu_dgd - deriv of viscosity wrt strain  *
                                                  * rate inv                                  */
                             dbl[DIM][MDE],      /* d_mu_dv                                   */
                             dbl[DIM][MDE],      /* d_mu_dmesh                                */
                             dbl[MDE],           /* d_mu_dT                                   */
                             dbl[MDE],           /* d_mu_dp                                   */
                             dbl[MAX_CONC][MDE], /* d_mu_dC                                   */
                             dbl *param);        /* user-defined parameter list               */

#endif /* GOMA_USER_MP_GEN_H */
