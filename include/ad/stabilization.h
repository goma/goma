/************************************************************************ *
* Goma - Multiphysics finite element software                             *
* Sandia National Laboratories                                            *
*                                                                         *
* Copyright (c) 2026 Goma Developers, National Technology & Engineering   *
*               Solutions of Sandia, LLC (NTESS)                          *
*                                                                         *
* Under the terms of Contract DE-NA0003525, the U.S. Government retains   *
* certain rights in this software.                                        *
*                                                                         *
* This software is distributed under the GNU General Public License.      *
* See LICENSE file.                                                       *
\************************************************************************/

#ifndef GOMA_AD_STABILIZATION_H
#define GOMA_AD_STABILIZATION_H
#include "std.h"
#include "mm_fill_stabilization.h"
#ifdef __cplusplus
#include "ad/structs.h"
void ad_supg_tau_shakib(ADType &supg_tau, int dim, dbl dt, ADType diffusivity, int interp_eqn);
void ad_get_metric_tensor(ADType B[DIM][DIM], int dim, int element_type, ADType G[DIM][DIM]);
void ad_only_tau_momentum_shakib(ADType &tau, int dim, dbl dt, int pspg_scale);
#endif
#endif