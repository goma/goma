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
#include "ad/stabilization.h"
#include "ad/viscosity.h"
#include "mm_fill_stabilization.h"
/* GOMA include files */
#include "density.h"
#include "el_elm.h"
#include "el_elm_info.h"
#include "el_geom.h"
#include "mm_as.h"
#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "mm_eh.h"
#include "mm_fill_energy.h"
#include "mm_fill_stress.h"
#include "mm_mp.h"
#include "mm_mp_structs.h"
#include "rf_fem.h"
#include "rf_fem_const.h"
#include "std.h"
#ifdef GOMA_ENABLE_SACADO

void ad_supg_tau_shakib(ADType &supg_tau, int dim, dbl dt, ADType diffusivity, int interp_eqn) {
  ADType G[DIM][DIM];

  ad_get_metric_tensor(ad_fv->B, dim, ei[pg->imtrx]->ielem_type, G);

  ADType v_d_gv = 0;
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      v_d_gv += fabs(ad_fv->v[i] * G[i][j] * ad_fv->v[j]);
    }
  }

  ADType diff_g_g = 0;
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      diff_g_g += G[i][j] * G[i][j];
    }
  }
  diff_g_g *= 9 * diffusivity * diffusivity;

  if (dt > 0) {
    supg_tau = 1.0 / (sqrt(4 / (dt * dt) + v_d_gv + diff_g_g));
  } else {
    supg_tau = 1.0 / (sqrt(v_d_gv + diff_g_g) + 1e-14);
  }
}

void ad_get_metric_tensor(ADType B[DIM][DIM], int dim, int element_type, ADType G[DIM][DIM]) {
  dbl adjustment[DIM][DIM] = {{0}};
  const dbl invroot3 = 0.577350269189626;
  const dbl tetscale = 0.629960524947437; // 0.5 * cubroot(2)

  switch (element_type) {
  case LINEAR_TRI:
    adjustment[0][0] = (invroot3) * 2;
    adjustment[0][1] = (invroot3) * -1;
    adjustment[1][0] = (invroot3) * -1;
    adjustment[1][1] = (invroot3) * 2;
    break;
  case LINEAR_TET:
    adjustment[0][0] = tetscale * 2;
    adjustment[0][1] = tetscale * 1;
    adjustment[0][2] = tetscale * 1;
    adjustment[1][0] = tetscale * 1;
    adjustment[1][1] = tetscale * 2;
    adjustment[1][2] = tetscale * 1;
    adjustment[2][0] = tetscale * 1;
    adjustment[2][1] = tetscale * 1;
    adjustment[2][2] = tetscale * 2;
    break;
  default:
    adjustment[0][0] = 1.0;
    adjustment[1][1] = 1.0;
    adjustment[2][2] = 1.0;
    break;
  }

  // G = B * adjustment * B^T where B = J^-1

  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      G[i][j] = 0;
      for (int k = 0; k < dim; k++) {
        for (int m = 0; m < dim; m++) {
          G[i][j] += B[i][k] * adjustment[k][m] * B[j][m];
        }
      }
    }
  }
}
void ad_only_tau_momentum_shakib(ADType &tau, int dim, dbl dt, int pspg_scale) {
  ADType G[DIM][DIM];
  dbl inv_rho = 1.0;
  DENSITY_DEPENDENCE_STRUCT d_rho_struct;
  DENSITY_DEPENDENCE_STRUCT *d_rho = &d_rho_struct;
  ADType gamma[DIM][DIM];
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      gamma[i][j] = ad_fv->grad_v[i][j] + ad_fv->grad_v[j][i];
    }
  }

  if (pspg_scale) {
    dbl rho = density(d_rho, dt);
    if (rho > 0.0) {
      inv_rho = 1.0 / rho;
    }
  }

  ad_get_metric_tensor(ad_fv->B, dim, ei[pg->imtrx]->ielem_type, G);

  ADType v_d_gv = 0;
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      v_d_gv += fabs((ad_fv->v[i] - ad_fv->x_dot[i]) * G[i][j] * (ad_fv->v[j] - ad_fv->x_dot[j]));
    }
  }

  ADType mu = ad_viscosity(gn, gamma);
  for (int mode = 0; mode < vn->modes; mode++) {
    ADType mup = ad_viscosity(ve[mode]->gn, gamma);
    mu += mup;
  }

  ADType coeff = (12.0 * mu * mu);

  ADType diff_g_g = 0;
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      diff_g_g += coeff * G[i][j] * G[i][j];
    }
  }

  if (pd->TimeIntegration != STEADY) {
    tau = inv_rho / (sqrt(4 / (dt * dt) + v_d_gv + diff_g_g));
  } else {
    tau = inv_rho / (sqrt(v_d_gv + diff_g_g) + 1e-14);
  }
}
#endif // GOMA_ENABLE_SACADO