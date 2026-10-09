/************************************************************************ *
* Goma - Multiphysics finite element software                             *
* Sandia National Laboratories                                            *
*                                                                         *
* Copyright (c) Goma Developers, National Technology & Engineering        *
*               Solutions of Sandia, LLC (NTESS)                          *
*                                                                         *
* Under the terms of Contract DE-NA0003525, the U.S. Government retains   *
* certain rights in this software.                                        *
*                                                                         *
* This software is distributed under the GNU General Public License.      *
* See LICENSE file.                                                       *
\************************************************************************/
#ifndef GOMA_AD_VISCOSITY_H
#define GOMA_AD_VISCOSITY_H
#ifdef GOMA_ENABLE_SACADO

#include "mm_as_structs.h"
#include "std.h"

GOMA_EXTERN dbl ad_viscosity_wrap(struct Generalized_Newtonian *gn_local);
GOMA_EXTERN dbl ad_sa_viscosity(struct Generalized_Newtonian *gn_local,
                                VISCOSITY_DEPENDENCE_STRUCT *d_mu);

#ifdef __cplusplus
#include "ad/structs.h"
ADType sst_viscosity(const ADType &Omega, const ADType &F2);
ADType ad_sa_viscosity(struct Generalized_Newtonian *gn_local);
ADType ad_numerical_viscosity(ADType s[DIM][DIM], /* total stress */
                              ADType gamma_cont[DIM][DIM],
                              int sdim); /* continuous shear rate */

ADType ad_arrhenius_simple_viscosity(struct Generalized_Newtonian *gn_local,
                                     ADType gamma_dot[DIM][DIM]);

ADType ad_carreau_arrhenius_viscosity(struct Generalized_Newtonian *gn_local,
                                      ADType gamma_dot[DIM][DIM]);
ADType ad_arrhenius_viscosity(struct Generalized_Newtonian *gn_local, ADType gamma_dot[DIM][DIM]);

ADType ad_ls_modulate_property(
    const ADType &p1, const ADType &p2, double width, double pm_minus, double pm_plus);
int ad_ls_modulate_viscosity(
    ADType &mu1, double mu2, double width, double pm_minus, double pm_plus, const int model);

ADType ad_bingham_viscosity(struct Generalized_Newtonian *gn_local, ADType gamma_dot[DIM][DIM]);
ADType ad_viscosity(struct Generalized_Newtonian *gn_local, ADType gamma_dot[DIM][DIM]);
#endif

#endif // GOMA_ENABLE_SACADO
#endif