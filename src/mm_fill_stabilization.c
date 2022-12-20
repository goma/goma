#include <float.h>
#include <math.h>
#include <string.h>

#include "mm_as.h"
#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "mm_eh.h"
#include "mm_fill_stabilization.h"
#include "mm_fill_stress.h"
#include "mm_fill_terms.h"
#include "mm_mp.h"
#include "mm_mp_const.h"
#include "mm_mp_structs.h"
#include "mm_qtensor_model.h"
#include "mm_viscosity.h"
#include "rf_fem.h"
#include "rf_fem_const.h"
#include "rf_solver.h"
#include "std.h"
#include "user_mp.h"

static dbl yzbeta1(dbl scale,
                   int dim,
                   dbl Y,
                   dbl Z,
                   dbl d_Z[MDE],
                   dbl beta,
                   dbl u,
                   dbl d_u[MDE],
                   dbl grad_u[DIM],
                   dbl d_grad_u[MDE][DIM],
                   dbl h_elem,
                   int interp_eqn,
                   dbl deriv[MDE]);

static dbl
yzbeta2(dbl scale, dbl Y, dbl Z, dbl d_Z[MDE], dbl deriv[MDE], dbl h_elem, int interp_eqn);

static const dbl DIFFUSION_EPSILON = 1e-8;

void get_metric_tensor(dbl B[DIM][DIM], int dim, int element_type, dbl G[DIM][DIM]) {
  dbl adjustment[DIM][DIM] = {{0}};
  const dbl invroot3 = 0.577350269189626;
  const dbl tetscale = 0.629960524947437; // 0.5 * cubroot(2)

  switch (element_type) {
  case LINEAR_TRI:
    adjustment[0][0] = (invroot3)*2;
    adjustment[0][1] = (invroot3) * -1;
    adjustment[1][0] = (invroot3) * -1;
    adjustment[1][1] = (invroot3)*2;
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

void get_metric_tensor_deriv(dbl B[DIM][DIM],
                             dbl dB[DIM][DIM][DIM][MDE],
                             int dim,
                             int interp_base,
                             int element_type,
                             dbl dG[DIM][DIM][DIM][MDE]) {
  dbl adjustment[DIM][DIM] = {{0}};
  const dbl invroot3 = 0.577350269189626;
  const dbl tetscale = 0.629960524947437; // 0.5 * cubroot(2)

  switch (element_type) {
  case LINEAR_TRI:
    adjustment[0][0] = (invroot3)*2;
    adjustment[0][1] = (invroot3) * -1;
    adjustment[1][0] = (invroot3) * -1;
    adjustment[1][1] = (invroot3)*2;
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

  for (int a = 0; a < dim; a++) {
    for (int b = 0; b < ei[pg->imtrx]->dof[interp_base + a]; b++) {
      for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
          dG[i][j][a][b] = 0;
          for (int k = 0; k < dim; k++) {
            for (int m = 0; m < dim; m++) {
              dG[i][j][a][b] += dB[i][k][a][b] * adjustment[k][m] * B[j][m];
              dG[i][j][a][b] += B[i][k] * adjustment[k][m] * dB[j][m][a][b];
            }
          }
        }
      }
    }
  }
}

void supg_tau_momentum_shakib(SUPG_momentum_terms *supg_terms, int dim, dbl dt) {
  dbl G[DIM][DIM];
  dbl gamma[DIM][DIM];
  dbl mu;
  VISCOSITY_DEPENDENCE_STRUCT d_mu_struct;
  VISCOSITY_DEPENDENCE_STRUCT *d_mu = &d_mu_struct;

  int interp_eqn = VELOCITY1;
  get_metric_tensor(bf[interp_eqn]->B, dim, ei[pg->imtrx]->ielem_type, G);

  dbl v_d_gv = 0;
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      v_d_gv += fabs(fv->v[i] * G[i][j] * fv->v[j]);
    }
  }

  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      gamma[i][j] = fv->grad_v[i][j] + fv->grad_v[j][i];
    }
  }

  mu = viscosity(gn, gamma, d_mu);

  dbl coeff = (12.0 * mu * mu);
  dbl diff_g_g = 0;
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      diff_g_g += coeff * G[i][j] * G[i][j];
    }
  }

  dbl d_v_d_gv[DIM][MDE];
  for (int a = 0; a < dim; a++) {
    for (int k = 0; k < ei[pg->imtrx]->dof[VELOCITY1]; k++) {
      d_v_d_gv[a][k] = 0.0;
      for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
          d_v_d_gv[a][k] += delta(a, i) * bf[VELOCITY1 + a]->phi[k] * G[i][j] * fv->v[j] +
                            delta(a, j) * fv->v[i] * G[i][j] * bf[VELOCITY1 + a]->phi[k];
        }
      }
    }
  }

  if (pd->TimeIntegration != STEADY) {
    supg_terms->supg_tau = 1.0 / (sqrt(4 / (dt * dt) + v_d_gv + diff_g_g));
  } else {
    supg_terms->supg_tau = 1.0 / (sqrt(v_d_gv + diff_g_g) + 1e-14);
  }

  dbl d_diff_g_g_dmu = 2.0 * diff_g_g / mu;
  dbl supg_tau_cubed = supg_terms->supg_tau * supg_terms->supg_tau * supg_terms->supg_tau;

  for (int a = 0; a < dim; a++) {
    for (int k = 0; k < ei[pg->imtrx]->dof[VELOCITY1]; k++) {
      supg_terms->d_supg_tau_dv[a][k] =
          -0.5 * (d_v_d_gv[a][k] + (d_diff_g_g_dmu * d_mu->v[a][k])) * supg_tau_cubed;
    }
  }

  if (pd->e[pg->imtrx][MESH_DISPLACEMENT1]) {
    dbl dG[DIM][DIM][DIM][MDE];
    get_metric_tensor_deriv(bf[MESH_DISPLACEMENT1]->B, bf[MESH_DISPLACEMENT1]->dB, dim,
                            MESH_DISPLACEMENT1, ei[pg->imtrx]->ielem_type, dG);
    for (int a = 0; a < dim; a++) {
      for (int k = 0; k < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1 + a]; k++) {
        dbl v_d_gv_dx = 0;
        for (int i = 0; i < dim; i++) {
          for (int j = 0; j < dim; j++) {
            v_d_gv_dx += fv->v[i] * dG[i][j][a][k] * fv->v[j];
          }
        }

        dbl diff_g_g_dx = 0;
        for (int i = 0; i < dim; i++) {
          for (int j = 0; j < dim; j++) {
            diff_g_g_dx += 2 * coeff * dG[i][j][a][k] * G[i][j];
          }
        }
        supg_terms->d_supg_tau_dX[a][k] =
            -0.5 * (v_d_gv_dx + diff_g_g_dx + d_mu->X[a][k] * d_diff_g_g_dmu) * supg_tau_cubed;
      }
    }
  }
  if (pd->e[pg->imtrx][TEMPERATURE]) {
    for (int k = 0; k < ei[pg->imtrx]->dof[TEMPERATURE]; k++) {
      supg_terms->d_supg_tau_dT[k] = -0.5 * (d_mu->T[k] * d_diff_g_g_dmu) * supg_tau_cubed;
    }
  }
  if (pd->e[pg->imtrx][PRESSURE]) {
    for (int k = 0; k < ei[pg->imtrx]->dof[PRESSURE]; k++) {
      supg_terms->d_supg_tau_dP[k] = -0.5 * (d_mu->P[k] * d_diff_g_g_dmu) * supg_tau_cubed;
    }
  }
  if (pd->e[pg->imtrx][FILL]) {
    for (int k = 0; k < ei[pg->imtrx]->dof[FILL]; k++) {
      supg_terms->d_supg_tau_dF[k] = -0.5 * (d_mu->F[k] * d_diff_g_g_dmu) * supg_tau_cubed;
    }
  }
  if (pd->e[pg->imtrx][BOND_EVOLUTION]) {
    for (int k = 0; k < ei[pg->imtrx]->dof[BOND_EVOLUTION]; k++) {
      supg_terms->d_supg_tau_dnn[k] = -0.5 * (d_mu->nn[k] * d_diff_g_g_dmu) * supg_tau_cubed;
    }
  }
  if (pd->e[pg->imtrx][MASS_FRACTION]) {
    for (int w = 0; w < pd->Num_Species_Eqn; w++) {
      for (int k = 0; k < ei[pg->imtrx]->dof[MASS_FRACTION]; k++) {
        supg_terms->d_supg_tau_dC[w][k] = -0.5 * (d_mu->C[w][k] * d_diff_g_g_dmu) * supg_tau_cubed;
      }
    }
  }
}

void supg_tau_shakib(SUPG_terms *supg_terms, int dim, dbl dt, dbl diffusivity, int interp_eqn) {
  dbl G[DIM][DIM];

  get_metric_tensor(bf[interp_eqn]->B, dim, ei[pg->imtrx]->ielem_type, G);

  dbl v_d_gv = 0;
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      v_d_gv += fabs(fv->v[i] * G[i][j] * fv->v[j]);
    }
  }

  dbl diff_g_g = 0;
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      diff_g_g += G[i][j] * G[i][j];
    }
  }
  diff_g_g *= 9 * diffusivity * diffusivity;

  dbl d_v_d_gv[DIM][MDE];
  for (int a = 0; a < dim; a++) {
    for (int k = 0; k < ei[pg->imtrx]->dof[VELOCITY1]; k++) {
      d_v_d_gv[a][k] = 0.0;
      for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
          d_v_d_gv[a][k] += delta(a, i) * bf[VELOCITY1 + a]->phi[k] * G[i][j] * fv->v[j] +
                            delta(a, j) * fv->v[i] * G[i][j] * bf[VELOCITY1 + a]->phi[k];
        }
      }
    }
  }

  if (dt > 0) {
    supg_terms->supg_tau = 1.0 / (sqrt(4 / (dt * dt) + v_d_gv + diff_g_g));
  } else {
    supg_terms->supg_tau = 1.0 / (sqrt(v_d_gv + diff_g_g) + 1e-14);
  }

  for (int a = 0; a < dim; a++) {
    for (int k = 0; k < ei[pg->imtrx]->dof[VELOCITY1]; k++) {
      supg_terms->d_supg_tau_dv[a][k] = -0.5 * d_v_d_gv[a][k] * supg_terms->supg_tau *
                                        supg_terms->supg_tau * supg_terms->supg_tau;
    }
  }

  if (pd->e[pg->imtrx][MESH_DISPLACEMENT1]) {
    dbl dG[DIM][DIM][DIM][MDE];
    get_metric_tensor_deriv(bf[MESH_DISPLACEMENT1]->B, bf[MESH_DISPLACEMENT1]->dB, dim,
                            MESH_DISPLACEMENT1, ei[pg->imtrx]->ielem_type, dG);
    for (int a = 0; a < dim; a++) {
      for (int k = 0; k < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1 + a]; k++) {
        dbl v_d_gv_dx = 0;
        for (int i = 0; i < dim; i++) {
          for (int j = 0; j < dim; j++) {
            v_d_gv_dx += fv->v[i] * dG[i][j][a][k] * fv->v[j];
          }
        }

        dbl diff_g_g_dx = 0;
        for (int i = 0; i < dim; i++) {
          for (int j = 0; j < dim; j++) {
            diff_g_g_dx += 2 * dG[i][j][a][k] * G[i][j];
          }
        }
        diff_g_g_dx *= 9 * diffusivity * diffusivity;
        supg_terms->d_supg_tau_dX[a][k] = -0.5 * (v_d_gv_dx + diff_g_g_dx) * supg_terms->supg_tau *
                                          supg_terms->supg_tau * supg_terms->supg_tau;
      }
    }
  }
}

void supg_tau_gauss_point(SUPG_terms *supg_terms,
                          int dim,
                          dbl diffusivity,
                          const PG_DATA *pg_data) {
  dbl vnorm = 0;

  for (int i = 0; i < VIM; i++) {
    vnorm += fv->v[i] * fv->v[i];
  }
  vnorm = sqrt(vnorm);

  dbl hk = 0;
  for (int a = 0; a < ei[pg->imtrx]->ielem_dim; a++) {
    hk += pg_data->hsquared[a];
  }
  /* This is the size of the element */
  hk = sqrt(hk / ((dbl)ei[pg->imtrx]->ielem_dim));

  dbl D = diffusivity;

  dbl hk_dX[DIM][MDE];
  for (int a = 0; a < dim; a++) {
    for (int j = 0; j < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1 + a]; j++) {
      dbl tmp = 0;
      for (int b = 0; b < dim; b++) {
        tmp +=
            (2 * pg_data->hhv[b][a] * pg_data->dhv_dxnode[b][j]) / (2 * sqrt(pg_data->hsquared[b]));
      }
      hk_dX[a][j] = tmp / dim;
    }
  }

  dbl Pek = 0.5 * vnorm * hk / (D + DBL_EPSILON);

  dbl eta = Pek;
  dbl eta_dX[DIM][MDE];
  dbl eta_dV[DIM][MDE];
  if (Pek > 1) {
    eta = 1;
    for (int i = 0; i < DIM; i++) {
      for (int j = 0; j < MDE; j++) {
        eta_dX[i][j] = 0;
        eta_dV[i][j] = 0;
      }
    }
  } else {
    for (int i = 0; i < DIM; i++) {
      for (int j = 0; j < MDE; j++) {
        if (pd->e[pg->imtrx][VELOCITY1 + i]) {
          eta_dV[i][j] =
              0.5 * 0.5 * hk * fv->v[i] * bf[VELOCITY1 + i]->phi[j] / (vnorm * D + 1e-16);
        }

        if (pd->e[pg->imtrx][MESH_DISPLACEMENT1 + i]) {
          eta_dX[i][j] = 0.5 * vnorm * hk_dX[i][j] / D;
        }
      }
    }
  }

  if (vnorm > 0) {
    supg_terms->supg_tau = 0.5 * hk * eta / vnorm;

    for (int a = 0; a < VIM; a++) {
      int var = VELOCITY1 + a;
      for (int j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        supg_terms->d_supg_tau_dv[a][j] =
            0.5 * hk * eta * fv->v[a] * bf[var]->phi[j] / (-vnorm * vnorm * vnorm) +
            0.5 * hk * eta_dV[a][j] / vnorm;
      }

      var = MESH_DISPLACEMENT1 + a;
      for (int j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        supg_terms->d_supg_tau_dX[a][j] =
            0.5 * hk_dX[a][j] * eta / vnorm + 0.5 * hk * eta_dX[a][j] / vnorm;
      }
    }

  } else {
    supg_terms->supg_tau = 0;
    for (int i = 0; i < DIM; i++) {
      for (int j = 0; j < MDE; j++) {
        supg_terms->d_supg_tau_dv[i][j] = 0.0;
      }
      for (int j = 0; j < MDE; j++) {
        supg_terms->d_supg_tau_dX[i][j] = 0.0;
      }
    }
  }
}

void supg_tau(SUPG_terms *supg_terms,
              int dim,
              dbl diffusivity,
              const PG_DATA *pg_data,
              dbl dt,
              int shakib,
              int interp_eqn) {
  if (shakib) {
    supg_tau_shakib(supg_terms, dim, dt, diffusivity, interp_eqn);
  } else {
    supg_tau_gauss_point(supg_terms, dim, diffusivity, pg_data);
  }
}

static dbl yzbeta1(dbl scale,
                   int dim,
                   dbl Y,
                   dbl Z,
                   dbl d_Z[MDE],
                   dbl beta,
                   dbl u,
                   dbl d_u[MDE],
                   dbl grad_u[DIM],
                   dbl d_grad_u[MDE][DIM],
                   dbl h_elem,
                   int interp_eqn,
                   dbl deriv[MDE]) {

  //  static const dbl EPSILON = 1e-10;
  dbl Y_inv = 1.0 / Y;
  //  dbl resid_scale = Y_inv * Z + EPSILON;
  dbl inner = 0;
  for (int i = 0; i < dim; i++) {
    inner += Y_inv * Y_inv * grad_u[i] * grad_u[i];
  }
  for (int i = 0; i < MDE; i++) {
    deriv[i] = 0;
  }
  return 1 * fabs(Z) * (1.0 / (sqrt(inner + 1e-12))) * h_elem * 0.5;
}

static dbl
yzbeta2(dbl scale, dbl Y, dbl Z, dbl d_Z[MDE], dbl deriv[MDE], dbl h_elem, int interp_eqn) {
  //  static const dbl EPSILON = 1e-10;
  for (int k = 0; k < ei[pg->imtrx]->dof[interp_eqn]; k++) {
    deriv[k] = 0;
  }
  dbl yzbeta = 1.0 * fabs(Z) * h_elem * h_elem * 0.025;
  return yzbeta;
}

dbl yzbeta(dbl scale,
           int dim,
           dbl Y,
           dbl Z,
           dbl d_Z[MDE],
           dbl beta,
           dbl u,
           dbl d_u[MDE],
           dbl grad_u[DIM],
           dbl d_grad_u[MDE][DIM],
           dbl h_elem,
           int interp_eqn,
           dbl deriv[MDE]) {
  dbl Y_inv = 1.0 / Y;

  dbl gradunit[DIM];
  dbl grad_u_norm = 0;

  for (int i = 0; i < dim; i++) {
    grad_u_norm += grad_u[i] * grad_u[i];
  }
  grad_u_norm = sqrt(grad_u_norm) + DBL_EPSILON;
  dbl inv_grad_u_norm = 1 / grad_u_norm;

  for (int j = 0; j < ei[pd->mi[interp_eqn]]->dof[interp_eqn]; j++) {
    dbl inner = 0;
    for (int i = 0; i < dim; i++) {
      inner += Y_inv * grad_u[i] * Y_inv * grad_u[i];
    }
    inner += DIFFUSION_EPSILON;

    dbl d_inner = 0;
    for (int i = 0; i < dim; i++) {
      d_inner += 2 * Y_inv * d_grad_u[j][i] * Y_inv * grad_u[i];
    }

    //    dbl scalar_part = Y_inv * u + DIFFUSION_EPSILON;

    //    dbl p = 1-beta;
    //    dbl q = (beta/2.0 - 1);

    dbl d_grad_u_norm = 0;
    for (int i = 0; i < dim; i++) {
      d_grad_u_norm += bf[interp_eqn]->grad_phi[j][i] * grad_u[i];
    }
    d_grad_u_norm *= inv_grad_u_norm;

    dbl d_inv_grad_u_norm = -inv_grad_u_norm * inv_grad_u_norm * d_grad_u_norm;

    for (int i = 0; i < dim; i++) {
      gradunit[i] = grad_u[i] * inv_grad_u_norm;
    }

    dbl h_dc = 0;
    for (int i = 0; i < ei[pd->mi[interp_eqn]]->dof[interp_eqn]; i++) {
      for (int j = 0; j < dim; j++) {
        h_dc += fabs(gradunit[j] * d_grad_u[i][j]);
      }
    }

    // h_dc = 2 / h_dc;

    //    dbl d_h_dc = 0;

    // deriv[j] = (Z+DIFFUSION_EPSILON) / (fabs(Z+DIFFUSION_EPSILON)) * d_Z[j] * pow(inner,q);
    // deriv[j] += fabs(Z+DIFFUSION_EPSILON) * d_inner * q * pow(inner, q-1);

    deriv[j] = d_inv_grad_u_norm;
  }
  return inv_grad_u_norm;
}

dbl yzbeta_model(int model,
                 dbl scale,
                 dbl beta,
                 int dim,
                 dbl Y,
                 dbl Z,
                 dbl d_Z[MDE],
                 dbl u,
                 dbl d_u[MDE],
                 dbl grad_u[DIM],
                 dbl d_grad_u[MDE][DIM],
                 dbl h_elem,
                 int interp_eqn,
                 dbl deriv[MDE]) {
  dbl dc = 0;
  switch (model) {
  case YZBETA_ONE:
    dc = yzbeta1(scale, dim, Y, Z, d_Z, 1.0, u, d_u, grad_u, d_grad_u, h_elem, interp_eqn, deriv);
    break;
  case YZBETA_TWO:
    dc = yzbeta2(scale, Y, Z, d_Z, deriv, h_elem, interp_eqn);
    break;
  case YZBETA_MIXED: {
    dbl deriv1[MDE] = {0};
    dbl deriv2[MDE] = {0};
    dbl dc1, dc2;
    dc1 = yzbeta1(scale, dim, Y, Z, d_Z, 1.0, u, d_u, grad_u, d_grad_u, h_elem, interp_eqn, deriv);
    dc2 = yzbeta2(scale, Y, Z, d_Z, deriv, h_elem, interp_eqn);
    for (int j = 0; j < ei[pd->mi[interp_eqn]]->dof[interp_eqn]; j++) {
      deriv[j] = 0.5 * (deriv1[j] + deriv2[j]);
    }
    dc = 0.5 * (dc1 + dc2);
  } break;
  case YZBETA_CUSTOM:
    dc = yzbeta(scale, dim, Y, Z, d_Z, beta, u, d_u, grad_u, d_grad_u, h_elem, interp_eqn, deriv);
    break;
  default:
    GOMA_EH(GOMA_ERROR, "Unknown YZBETA Model");
    break;
  }

  return dc;
}
void get_supg_tau(SUPG_terms *supg_terms, int dim, dbl diffusivity, PG_DATA *pg_data) {
  dbl vnorm = 0;

  for (int i = 0; i < VIM; i++) {
    vnorm += fv->v[i] * fv->v[i];
  }
  vnorm = sqrt(vnorm);

  dbl hk = 0;
  for (int i = 0; i < dim; i++) {
    hk += sqrt(pg_data->hsquared[i]);
  }

  hk /= (dbl)dim;

  dbl D = diffusivity;

  dbl hk_dX[DIM][MDE];
  for (int a = 0; a < dim; a++) {
    for (int j = 0; j < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1 + a]; j++) {
      dbl tmp = 0;
      for (int b = 0; b < dim; b++) {
        tmp +=
            (2 * pg_data->hhv[b][a] * pg_data->dhv_dxnode[b][j]) / (2 * sqrt(pg_data->hsquared[b]));
      }
      hk_dX[a][j] = tmp / dim;
    }
  }

  dbl Pek = 0.5 * vnorm * hk / D;

  dbl eta = Pek;
  dbl eta_dX[DIM][MDE];
  dbl eta_dV[DIM][MDE];
  if (Pek > 1) {
    eta = 1;
    for (int i = 0; i < DIM; i++) {
      for (int j = 0; j < MDE; j++) {
        eta_dX[i][j] = 0;
        eta_dV[i][j] = 0;
      }
    }
  } else {
    for (int i = 0; i < DIM; i++) {
      for (int j = 0; j < MDE; j++) {
        if (pd->e[pg->imtrx][VELOCITY1 + i]) {
          eta_dV[i][j] = 0.5 * 0.5 * hk * fv->v[i] * bf[VELOCITY1 + i]->phi[j] / (vnorm * D);
        }

        if (pd->e[pg->imtrx][MESH_DISPLACEMENT1 + i]) {
          eta_dX[i][j] = 0.5 * vnorm * hk_dX[i][j] / D;
        }
      }
    }
  }

  if (vnorm > 0) {
    supg_terms->supg_tau = 0.5 * hk * eta / vnorm;

    for (int a = 0; a < VIM; a++) {
      int var = VELOCITY1 + a;
      for (int j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        supg_terms->d_supg_tau_dv[a][j] =
            0.5 * hk * eta * fv->v[a] * bf[var]->phi[j] / (-vnorm * vnorm * vnorm) +
            0.5 * hk * eta_dV[a][j] / vnorm;
      }

      var = MESH_DISPLACEMENT1 + a;
      for (int j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        supg_terms->d_supg_tau_dX[a][j] =
            0.5 * hk_dX[a][j] * eta / vnorm + 0.5 * hk * eta_dX[a][j] / vnorm;
      }
    }

  } else {
    supg_terms->supg_tau = 0;
    for (int i = 0; i < DIM; i++) {
      for (int j = 0; j < MDE; j++) {
        supg_terms->d_supg_tau_dv[i][j] = 0.0;
      }
      for (int j = 0; j < MDE; j++) {
        supg_terms->d_supg_tau_dX[i][j] = 0.0;
      }
    }
  }
}

int calc_pspg(dbl pspg[DIM],
              PSPG_DEPENDENCE_STRUCT *d_pspg,
              dbl time_value, /* current time */
              dbl tt,         /* parameter to vary time integration from
                                                 explicit (tt = 1) to implicit (tt = 0)    */
              dbl dt,         /* current time step size                    */
              const PG_DATA *pg_data) {
  const dbl h_elem_avg = pg_data->h_elem_avg;
  const dbl *hsquared = pg_data->hsquared; /* element size information for PSPG         */
  const dbl U_norm = pg_data->U_norm;      /* global velocity norm for PSPG calcs       */
  const dbl mu_avg = pg_data->mu_avg;      /* element viscosity for PSPG calculations   */
  const dbl rho_avg = pg_data->rho_avg;    /* element density for PSPG calculations   */
  const dbl *v_avg = pg_data->v_avg;       /* element velocity for PSPG calculations   */

  int dim;
  int p, a, b, c;
  int var;
  int w, j;

  int pspg_global;
  int pspg_local;

  dbl *v = fv->v; /* Velocity field. */
  dbl *grad_P = fv->grad_P;

  /*
   * Variables for vicosity and derivative
   */
  dbl mu;
  VISCOSITY_DEPENDENCE_STRUCT d_mu_struct; /* viscosity dependence */
  VISCOSITY_DEPENDENCE_STRUCT *d_mu = &d_mu_struct;

  /*
   * density and sensitivity terms
   */
  dbl rho;
  dbl rho_t;
  DENSITY_DEPENDENCE_STRUCT d_rho_struct; /* density dependence */
  DENSITY_DEPENDENCE_STRUCT *d_rho = &d_rho_struct;

  dbl f[DIM];                                  /* Body force. */
  MOMENTUM_SOURCE_DEPENDENCE_STRUCT df_struct; /* Body force dependence */
  MOMENTUM_SOURCE_DEPENDENCE_STRUCT *df = &df_struct;

  /*
   * Interpolation functions...
   */

  dbl phi_j;

  /*
   * Variables for Pressure Stabilization Petrov-Galerkin...
   */
  int meqn, meqn1;
  int r;
  int v_s[MAX_MODES][DIM][DIM], v_g[DIM][DIM];

  dbl mass, advection, diffusion, source, porous;

  dbl momentum[DIM]; /* momentum residual for PSPG */
  dbl x_dot[DIM];
  dbl v_dot[DIM];
  dbl *grad_v[DIM];
  dbl div_s[DIM];
  dbl div_G[DIM];

  /* variables for Brinkman porous flow */
  dbl por = 0, por2 = 0, per = 0, vis = 0, dvis_dT[MDE], sc = 0, speed = 0;

  dbl h_elem = 0;
  dbl Re;
  dbl tau_pspg = 0.;
  dbl tau_pspg1 = 0.;
  dbl d_tau_pspg_dX[DIM][MDE];
  dbl d_tau_pspg_dv[DIM][MDE];

  dbl hh_siz, vv_speed;

  dbl gamma[DIM][DIM]; /* shrearrate tensor based on velocity */

  int w0 = 0;

  int mode;

  /* For particle momentum model.
   */
  int species;              /* species number for particle phase,  */
  dbl ompvf;                /* 1 - partical volume fraction */
  int particle_momentum_on; /* boolean. */
  /* particle stress for suspension balance model*/
  dbl tau_p[DIM][DIM];
  dbl d_tau_p_dv[DIM][DIM][DIM][MDE];
  dbl d_tau_p_dvd[DIM][DIM][DIM][MDE];
  dbl d_tau_p_dy[DIM][DIM][MAX_CONC][MDE];
  dbl d_tau_p_dmesh[DIM][DIM][DIM][MDE];
  dbl d_tau_p_dp[DIM][DIM][MDE];

  static int is_initialized = FALSE;

  dim = pd->Num_Dim;

  /* initialize */
  for (a = 0; a < DIM; a++)
    pspg[a] = 0.;

  if (d_pspg == NULL) {
    /* I guess we won't be needing these! */
    d_mu = NULL;
    d_rho = NULL;
    df = NULL;
  } else if (!is_initialized) {

    memset(d_pspg->v, 0, sizeof(dbl) * DIM * DIM * MDE);
    memset(d_pspg->X, 0, sizeof(dbl) * DIM * DIM * MDE);
    memset(d_pspg->T, 0, sizeof(dbl) * DIM * MDE);
    memset(d_pspg->P, 0, sizeof(dbl) * DIM * MDE);
    memset(d_pspg->C, 0, sizeof(dbl) * DIM * MAX_CONC * MDE);
    memset(d_pspg->S, 0, sizeof(dbl) * DIM * MAX_MODES * DIM * DIM * MDE);
    memset(d_pspg->g, 0, sizeof(dbl) * DIM * DIM * DIM * MDE);
  }

  /* This is the flag for the standard global PSPG */
  if (PSPG == 1) {
    pspg_global = TRUE;
    pspg_local = FALSE;
  }
  /* This is the flag for the standard local PSPG */
  else if (PSPG == 2) {
    pspg_global = FALSE;
    pspg_local = TRUE;
  } else if (PSPG == 3) { // Shakib
    pspg_global = FALSE;
    pspg_local = FALSE;
  } else {
    return 0;
  }

  if (pd->gv[POLYMER_STRESS11]) {
    stress_eqn_pointer(v_s);

    v_g[0][0] = VELOCITY_GRADIENT11;
    v_g[0][1] = VELOCITY_GRADIENT12;
    v_g[1][0] = VELOCITY_GRADIENT21;
    v_g[1][1] = VELOCITY_GRADIENT22;
    v_g[0][2] = VELOCITY_GRADIENT13;
    v_g[1][2] = VELOCITY_GRADIENT23;
    v_g[2][0] = VELOCITY_GRADIENT31;
    v_g[2][1] = VELOCITY_GRADIENT32;
    v_g[2][2] = VELOCITY_GRADIENT33;
  }

  /* initialize dependencies */
  memset(d_tau_pspg_dv, 0, sizeof(dbl) * DIM * MDE);
  memset(d_tau_pspg_dX, 0, sizeof(dbl) * DIM * MDE);

  if (cr->MassFluxModel == DM_SUSPENSION_BALANCE && PSPG) {
    w0 = gn->sus_species_no;
    /* This is the divergence of the particle stress  */
    /* divergence_particle_stress(div_tau_p, d_div_tau_p_dgd, d_div_tau_p_dy,
         d_div_tau_p_dv, d_div_tau_p_dX, w0); */
  }

  if (pd->e[pg->imtrx][R_PMOMENTUM1]) {
    particle_momentum_on = 1;
    species = (int)mp->u_density[0];
    ompvf = 1.0 - fv->c[species];
  } else {
    particle_momentum_on = 0;
    species = -1;
    ompvf = 1.0;
  }

  // Global average for pspg_global's element size
  h_elem = h_elem_avg;

  /*** Density ***/
  rho = density(d_rho, time_value);

  if (pspg_global) {

    /* Now calculate the element Reynolds number based on a global
     * norm of the velocity and determine tau_pspg discretely from Re
     * The global version has no Jacobian dependencies
     */
    Re = rho * U_norm * h_elem / (2.0 * mu_avg);

    if (Re <= 3.0) {
      tau_pspg = PS_scaling * h_elem * h_elem / (12.0 * mu_avg);
    } else if (Re > 3.0) {
      tau_pspg = PS_scaling * h_elem / (2.0 * rho * U_norm);
    }
  } else if (pspg_local) {
    hh_siz = 0.;
    for (p = 0; p < dim; p++) {
      hh_siz += hsquared[p];
    }
    // Average value of h**2 in the element
    hh_siz = hh_siz / ((dbl)dim);

    // Average value of v**2 in the element
    vv_speed = 0.0;
    for (a = 0; a < WIM; a++) {
      vv_speed += v_avg[a] * v_avg[a];
    }

    // Use vv_speed and hh_siz for tau_pspg, note it has a continuous dependence on Re
    tau_pspg1 = rho_avg * rho_avg * vv_speed / hh_siz + (9.0 * mu_avg * mu_avg) / (hh_siz * hh_siz);
    if (pd->TimeIntegration != STEADY) {
      tau_pspg1 += 4.0 / (dt * dt);
    }
    tau_pspg = PS_scaling / sqrt(tau_pspg1);

    // tau_pspg derivatives wrt v from vv_speed
    if (d_pspg != NULL && pd->v[pg->imtrx][VELOCITY1]) {
      for (b = 0; b < dim; b++) {
        var = VELOCITY1 + b;
        if (pd->v[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            d_tau_pspg_dv[b][j] = -tau_pspg / tau_pspg1;
            d_tau_pspg_dv[b][j] *= rho_avg * rho_avg / hh_siz * v_avg[b] * pg_data->dv_dnode[b][j];
          }
        }
      }
    }

    // tau_pspg derivatives wrt mesh from hh_siz
    if (d_pspg != NULL && pd->v[pg->imtrx][MESH_DISPLACEMENT1] && pspg_local) {
      for (b = 0; b < dim; b++) {
        var = MESH_DISPLACEMENT1 + b;
        if (pd->v[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            d_tau_pspg_dX[b][j] = tau_pspg / tau_pspg1;
            d_tau_pspg_dX[b][j] *=
                (rho_avg * rho_avg * vv_speed + 18.0 * mu_avg * mu_avg / hh_siz) /
                (hh_siz * hh_siz);
            d_tau_pspg_dX[b][j] *= pg_data->hhv[b][b] * pg_data->dhv_dxnode[b][j] / ((dbl)dim);
          }
        }
      }
    }
  } else if (PSPG == 3) { // shakib

    for (a = 0; a < VIM; a++) {
      for (b = 0; b < VIM; b++) {
        gamma[a][b] = fv_old->grad_v[a][b] + fv_old->grad_v[b][a];
      }
    }

    // check if rho exists, otherwise set to 1 for pspg_tau
    dbl pspg_rho = rho;
    if (!DOUBLE_NONZERO(rho)) {
      pspg_rho = 1.0;
    }

    mu = viscosity(gn, gamma, NULL);

    dbl mup[MAX_MODES];
    if (pd->gv[POLYMER_STRESS11]) {
      for (int mode = 0; mode < vn->modes; mode++) {
        mup[mode] = viscosity(ve[mode]->gn, gamma, NULL);
      }
    }

    dbl G[DIM][DIM];
    get_metric_tensor(bf[pd->ShapeVar]->B, pd->Num_Dim, ei[pg->imtrx]->ielem_type, G);

    dbl tau_time = 0;
    // time term
    if (pd->TimeIntegration != STEADY) {
      tau_time += 4 * pspg_rho * pspg_rho / (dt * dt);
    }

    // advection
    dbl tau_adv = 0;
    for (int i = 0; i < dim; i++) {
      for (int j = 0; j < dim; j++) {
        // tau_adv += rho * rho * fv->v[i] * G[i][j] * fv->v[j];
        tau_adv += pspg_rho * pspg_rho * fv_old->v[i] * G[i][j] * fv_old->v[j];
      }
    }

    // diffusion
    dbl tau_diff = 0;
    dbl mu_total = mu;
    for (int mode = 0; mode < vn->modes; mode++) {
      mu_total += mup[mode];
    }
    dbl coeff = 12 * (mu_total * mu_total);
    for (int i = 0; i < dim; i++) {
      for (int j = 0; j < dim; j++) {
        tau_diff += coeff * G[i][j] * G[i][j];
      }
    }

    tau_pspg = PS_scaling * pspg_rho / sqrt(tau_time + tau_adv + tau_diff);

    // d/dx 1/sqrt(f(x)) => - f'(x) / (2 * f(x)^(3/2))
    if (d_pspg != NULL && pd->v[pg->imtrx][VELOCITY1]) {
      for (b = 0; b < dim; b++) {
        var = VELOCITY1 + b;
        if (pd->v[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            dbl tau_adv_dv = 0;
            for (int a = 0; a < dim; a++) {
              tau_adv_dv += pspg_rho * pspg_rho * bf[var]->phi[j] * G[b][a] * fv->v[a];
              tau_adv_dv += pspg_rho * pspg_rho * bf[var]->phi[j] * G[a][b] * fv->v[a];
            }

            // d_tau_pspg_dv[b][j] = -PS_scaling * pspg_rho * 0.5 * tau_pspg * tau_pspg * tau_pspg *
            // tau_adv_dv;
            d_tau_pspg_dv[b][j] =
                -PS_scaling * pspg_rho * 0.5 * (tau_adv_dv)*tau_pspg * tau_pspg * tau_pspg;
          }
        }
      }
    }
    if (d_pspg != NULL && pd->v[pg->imtrx][MESH_DISPLACEMENT1]) {
      dbl dG[DIM][DIM][DIM][MDE];
      get_metric_tensor_deriv(bf[MESH_DISPLACEMENT1]->B, bf[MESH_DISPLACEMENT1]->dB, dim,
                              MESH_DISPLACEMENT1, ei[pg->imtrx]->ielem_type, dG);
      for (b = 0; b < dim; b++) {
        var = MESH_DISPLACEMENT1 + b;
        if (pd->v[pg->imtrx][var]) {
          for (int k = 0; k < ei[pg->imtrx]->dof[var]; k++) {
            dbl tau_adv_dx = 0;
            for (int i = 0; i < dim; i++) {
              for (int j = 0; j < dim; j++) {
                tau_adv_dx += pspg_rho * pspg_rho * fv->v[i] * dG[i][j][b][k] * fv->v[j];
              }
            }
            dbl tau_diff_dx = 0;
            for (int i = 0; i < dim; i++) {
              for (int j = 0; j < dim; j++) {
                tau_diff_dx += coeff * 2 * dG[i][j][b][k] * G[i][j];
              }
            }
            d_tau_pspg_dX[b][k] = -PS_scaling * pspg_rho * 0.5 * (tau_adv_dx + tau_diff_dx) *
                                  tau_pspg * tau_pspg * tau_pspg;
          }
        }
      }
    }
  }

  for (a = 0; a < VIM; a++)
    grad_v[a] = fv->grad_v[a];

  /* load up shearrate tensor based on velocity */
  for (a = 0; a < VIM; a++) {
    for (b = 0; b < VIM; b++) {
      gamma[a][b] = grad_v[a][b] + grad_v[b][a];
    }
  }

  /*
   * get viscosity for velocity second derivative/diffusion
   * term in PSPG stuff
   */
  mu = viscosity(gn, gamma, d_mu);

  /* get variables we will need for momentum residual */

  for (a = 0; a < WIM; a++) {
    if (pd->TimeIntegration != STEADY && pd->v[pg->imtrx][MESH_DISPLACEMENT1 + a]) {
      x_dot[a] = fv_dot->x[a];
    } else {
      x_dot[a] = 0.;
    }

    if (pd->TimeIntegration != STEADY) {
      v_dot[a] = fv_dot->v[a];
    } else {
      v_dot[a] = 0.;
    }
  }

  for (p = 0; p < WIM; p++)
    div_s[p] = 0.;

  if (pd->gv[POLYMER_STRESS11]) {
    if (vn->evssModel == LOG_CONF || vn->evssModel == LOG_CONF_GRADV ||
        vn->evssModel == LOG_CONF_TRANSIENT || vn->evssModel == LOG_CONF_TRANSIENT_GRADV) {
      for (mode = 0; mode < vn->modes; mode++) {
        dbl lambda = 0.0;
        if (ve[mode]->time_constModel == CONSTANT) {
          lambda = ve[mode]->time_const;
        }
        dbl mup = viscosity(ve[mode]->gn, gamma, NULL);
        int dofs = ei[upd->matrix_index[v_s[mode][0][0]]]->dof[v_s[mode][0][0]];
        dbl grad_S[DIM][DIM][DIM] = {{{0.0}}};
        dbl s[MDE][DIM][DIM];
        dbl exp_s[MDE][DIM][DIM] = {{{0.0}}};
        dbl eig_values[DIM];
        dbl R[DIM][DIM];
        for (int k = 0; k < dofs; k++) {
          if (pg->imtrx == upd->matrix_index[POLYMER_STRESS11] &&
              (vn->evssModel == LOG_CONF_TRANSIENT_GRADV || vn->evssModel == LOG_CONF_TRANSIENT)) {
            for (int i = 0; i < VIM; i++) {
              for (int j = 0; j < VIM; j++) {
                if (j >= i) {
                  s[k][i][j] = *esp_old->S[mode][i][j][k];
                } else {
                  s[k][i][j] = *esp_old->S[mode][j][i][k];
                }
              }
            }
          } else {
            for (int i = 0; i < VIM; i++) {
              for (int j = 0; j < VIM; j++) {
                if (j >= i) {
                  s[k][i][j] = *esp->S[mode][i][j][k];
                } else {
                  s[k][i][j] = *esp->S[mode][j][i][k];
                }
              }
            }
          }
          compute_exp_s(s[k], exp_s[k], eig_values, R);
        }
        for (int p = 0; p < VIM; p++) {
          for (int q = 0; q < VIM; q++) {
            for (int r = 0; r < VIM; r++) {
              grad_S[r][p][q] = 0.;
              for (int i = 0; i < dofs; i++) {
                if (p <= q) {
                  grad_S[r][p][q] += exp_s[i][p][q] * bf[POLYMER_STRESS11]->grad_phi[i][r];
                } else {
                  grad_S[r][p][q] += exp_s[i][q][p] * bf[POLYMER_STRESS11]->grad_phi[i][r];
                }
              }
            }
          }
        }
        dbl div_exp_s[DIM];
        for (int r = 0; r < dim; r++) {
          div_exp_s[r] = 0.0;

          for (int q = 0; q < dim; q++) {
            div_exp_s[r] += grad_S[q][q][r];
          }
        }
        for (p = 0; p < WIM; p++) {
          div_s[p] += (mup / lambda) * div_exp_s[p];
        }
      }
    } else {
      for (p = 0; p < WIM; p++) {
        for (mode = 0; mode < vn->modes; mode++) {
          div_s[p] += fv->div_S[mode][p];
        }
      }
    }
  }

  if (!is_initialized) {
    memset(tau_p, 0, sizeof(dbl) * DIM * DIM);
    memset(d_tau_p_dv, 0, sizeof(dbl) * DIM * DIM * DIM * MDE);
    memset(d_tau_p_dvd, 0, sizeof(dbl) * DIM * DIM * DIM * MDE);
    memset(d_tau_p_dy, 0, sizeof(dbl) * DIM * DIM * MAX_CONC * MDE);
    memset(d_tau_p_dmesh, 0, sizeof(dbl) * DIM * DIM * DIM * MDE);
    memset(d_tau_p_dp, 0, sizeof(dbl) * DIM * DIM * MDE);
  }

  if (cr->MassFluxModel == DM_SUSPENSION_BALANCE || cr->MassFluxModel == HYDRODYNAMIC_QTENSOR_OLD ||
      cr->MassFluxModel == HYDRODYNAMIC_QTENSOR)
    particle_stress(tau_p, d_tau_p_dv, d_tau_p_dvd, d_tau_p_dy, d_tau_p_dmesh, d_tau_p_dp, w0);

  if (pd->gv[VELOCITY_GRADIENT11]) {
    for (p = 0; p < WIM; p++) {
      div_G[p] = fv->div_G[p];
    }
  } else {
    for (p = 0; p < WIM; p++) {
      div_G[p] = 0.;
    }
  }

  if (pd->e[upd->matrix_index[R_MOMENTUM1]][R_MOMENTUM1] & T_POROUS_BRINK) {
    if (mp->PorousMediaType != POROUS_BRINKMAN)
      GOMA_WH(-1, "Set Porous term multiplier in continuous medium");
    /* Short-hand notation for the four parameters in the Brinkman Equation. */
    por = mp->porosity;
    por2 = por * por;
    per = mp->permeability;
    if (mp->FlowingLiquidViscosityModel == CONSTANT) {
      /* Do nothing */
    } else if (mp->FlowingLiquidViscosityModel == MOLTEN_GLASS) {
      molten_glass_viscosity(&(mp->FlowingLiquid_viscosity), dvis_dT,
                             mp->u_FlowingLiquid_viscosity);
    } else if (mp->FlowingLiquidViscosityModel == USER) {
      usr_FlowingLiquidViscosity(mp->u_FlowingLiquid_viscosity);
      var = TEMPERATURE;
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        dvis_dT[j] = mp->d_FlowingLiquid_viscosity[var] * bf[var]->phi[j];
      }
    } else {
      GOMA_EH(GOMA_ERROR, "Don't recognize your FlowingLiquidViscosity model");
    }
    vis = mp->FlowingLiquid_viscosity;
    sc = mp->Inertia_coefficient;
  } else {
    por = 1.;
    por2 = 1.;
    per = 1.;
    vis = mp->viscosity;
    sc = 0.;
  }

  /* for porous media stuff */
  speed = 0.0;
  for (a = 0; a < WIM; a++) {
    speed += v[a] * v[a];
  }
  speed = sqrt(speed);

  /* get momentum source term */
  momentum_source_term(f, df, time_value);

  if (pd->gv[R_PMOMENTUM1]) {
    rho_t = ompvf * rho;
    meqn1 = R_PMOMENTUM1;
  } else {
    rho_t = rho;
    meqn1 = R_MOMENTUM1;
  }

  for (a = 0; a < WIM; a++) {
    meqn = meqn1 + a;

    mass = 0.;
    if ((pd->e[upd->matrix_index[meqn]][meqn] & T_MASS) && (pd->TimeIntegration != STEADY)) {
      mass = rho_t * v_dot[a] / por;
      mass *= pd->etm[upd->matrix_index[meqn]][meqn][(LOG2_MASS)];
    }

    advection = 0.;
    if (pd->e[upd->matrix_index[meqn]][meqn] & T_ADVECTION) {
      for (p = 0; p < WIM; p++) {
        advection += rho_t * (v[p] - x_dot[p]) * grad_v[p][a] / por2;
      }
      advection *= pd->etm[upd->matrix_index[meqn]][meqn][(LOG2_ADVECTION)];
    }

    diffusion = 0.;
    if (pd->e[upd->matrix_index[meqn]][meqn] & T_DIFFUSION) {
      diffusion = grad_P[a] - div_s[a];
      /*diffusion  -= div_tau_p[a]  */
      diffusion -= mu * div_G[a];
      diffusion *= pd->etm[upd->matrix_index[meqn]][meqn][(LOG2_DIFFUSION)];
    }

    source = 0.;
    if (pd->e[upd->matrix_index[meqn]][meqn] & T_SOURCE) {
      source = -f[a];
      source *= pd->etm[upd->matrix_index[meqn]][meqn][(LOG2_SOURCE)];
    }

    porous = 0.;
    if (pd->e[upd->matrix_index[meqn]][meqn] & T_POROUS_BRINK) {
      porous = v[a] * (rho_t * sc * speed / sqrt(per) + vis / per);
      porous *= pd->etm[upd->matrix_index[meqn]][meqn][(LOG2_POROUS_BRINK)];
    }

    momentum[a] = mass + advection + diffusion + source + porous;
    pspg[a] = tau_pspg * momentum[a];

    if (d_pspg != NULL && pd->v[pg->imtrx][VELOCITY1]) {
      for (b = 0; b < WIM; b++) {
        var = VELOCITY1 + b;
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          phi_j = bf[var]->phi[j];

          mass = 0.;
          if ((pd->e[upd->matrix_index[meqn]][meqn] & T_MASS) && (pd->TimeIntegration != STEADY)) {
            mass = rho_t / por * (1. + 2. * tt) * phi_j / dt * (dbl)delta(a, b);
            mass *= pd->etm[upd->matrix_index[meqn]][meqn][(LOG2_MASS)];
          }

          advection = 0.;
          if (pd->e[upd->matrix_index[meqn]][meqn] & T_ADVECTION) {
            advection = phi_j * grad_v[b][a];
            for (p = 0; p < WIM; p++) {
              advection += (v[p] - x_dot[p]) * bf[var]->grad_phi_e[j][b][p][a];
            }
            advection *= rho_t / por2;
            advection *= pd->etm[upd->matrix_index[meqn]][meqn][(LOG2_ADVECTION)];
          }

          diffusion = 0.;
          if (pd->e[upd->matrix_index[meqn]][meqn] & T_DIFFUSION) {
            diffusion -= d_mu->v[b][j] * div_G[a];
            /*diffusion -= d_div_tau_p_dv[a][b][j];*/
            diffusion *= pd->etm[upd->matrix_index[meqn]][meqn][(LOG2_DIFFUSION)];
          }

          source = 0.;
          if (pd->e[upd->matrix_index[meqn]][meqn] & T_SOURCE) {
            source -= df->v[a][b][j] * pd->etm[upd->matrix_index[meqn]][meqn][(LOG2_SOURCE)];
          }

          porous = 0.;
          if (pd->e[upd->matrix_index[meqn]][meqn] & T_POROUS_BRINK) {
            porous = (rho_t * sc * speed / sqrt(per) + vis / per) * phi_j;
            for (p = 0; p < WIM; p++) {
              porous += rho_t * sc / sqrt(per) * 2. * v[p] * v[a];
            }
            porous *= pd->etm[upd->matrix_index[meqn]][meqn][(LOG2_POROUS_BRINK)];
          }
          d_pspg->v[a][b][j] = tau_pspg * (mass + advection + diffusion + source + porous) +
                               d_tau_pspg_dv[b][j] * momentum[a];
        }
      }
    }

    if (d_pspg != NULL && pd->v[pg->imtrx][MESH_DISPLACEMENT1]) {
      for (b = 0; b < WIM; b++) {
        var = MESH_DISPLACEMENT1 + b;
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          phi_j = bf[var]->phi[j];

          advection = 0.;
          if ((pd->e[upd->matrix_index[meqn]][meqn] & T_ADVECTION)) {
            if (pd->TimeIntegration != STEADY) {
              advection = -(1. + 2. * tt) * phi_j / dt * grad_v[b][a];
            }
            for (p = 0; p < WIM; p++) {
              advection += (v[p] - x_dot[p]) * fv->d_grad_v_dmesh[p][a][b][j];
            }
            advection *= rho_t / por2;
            advection *= pd->etm[upd->matrix_index[meqn]][meqn][(LOG2_ADVECTION)];
          }

          diffusion = 0.;
          if (pd->e[upd->matrix_index[meqn]][meqn] & T_DIFFUSION) {
            diffusion = fv->d_grad_P_dmesh[a][b][j] - d_mu->X[b][j] * div_G[a] -
                        mu * fv->d_div_G_dmesh[a][b][j];
            /* diffusion -= d_div_tau_p_dX[a][b][j];*/
            for (mode = 0; mode < vn->modes; mode++) {
              diffusion -= fv->d_div_S_dmesh[mode][a][b][j];
            }
            diffusion *= pd->etm[upd->matrix_index[meqn]][meqn][(LOG2_DIFFUSION)];
          }

          source = 0.;
          if (pd->e[upd->matrix_index[meqn]][meqn] & T_SOURCE) {
            source -= df->X[a][b][j] * pd->etm[upd->matrix_index[meqn]][meqn][(LOG2_SOURCE)];
          }

          d_pspg->X[a][b][j] =
              tau_pspg * (advection + diffusion + source) + d_tau_pspg_dX[b][j] * momentum[a];
        }
      }
    }

    var = TEMPERATURE;
    if (d_pspg != NULL && pd->v[pg->imtrx][var]) {
      dbl d_rho_t_dT;

      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        phi_j = bf[var]->phi[j];

        d_rho_t_dT = ompvf * d_rho->T[j];

        mass = 0.;
        if ((pd->e[upd->matrix_index[meqn]][meqn] & T_MASS) && (pd->TimeIntegration != STEADY)) {
          mass = d_rho_t_dT / por * v_dot[a];
          mass *= pd->etm[upd->matrix_index[meqn]][meqn][(LOG2_MASS)];
        }

        advection = 0.;
        if (pd->e[upd->matrix_index[meqn]][meqn] & T_ADVECTION) {
          advection = 0.;
          for (p = 0; p < WIM; p++) {
            advection += (v[p] - x_dot[p]) * grad_v[p][a];
          }
          advection *= d_rho_t_dT / por2;
          advection *= pd->etm[upd->matrix_index[meqn]][meqn][(LOG2_ADVECTION)];
        }

        diffusion = 0.;
        if (pd->e[upd->matrix_index[meqn]][meqn] & T_DIFFUSION) {
          diffusion -= d_mu->T[j] * div_G[a];
          diffusion *= pd->etm[upd->matrix_index[meqn]][meqn][(LOG2_DIFFUSION)];
        }

        source = 0.;
        if (pd->e[upd->matrix_index[meqn]][meqn] & T_SOURCE) {
          source -= df->T[a][j] * pd->etm[upd->matrix_index[meqn]][meqn][(LOG2_SOURCE)];
        }

        porous = 0.;
        if (pd->e[upd->matrix_index[meqn]][meqn] & T_POROUS_BRINK) {
          porous = v[a] * (d_rho_t_dT * sc * speed / sqrt(per) + dvis_dT[j] / per);
          porous *= pd->etm[upd->matrix_index[meqn]][meqn][(LOG2_POROUS_BRINK)];
        }
        d_pspg->T[a][j] = tau_pspg * (mass + advection + diffusion + source + porous);
      }
    }

    var = PRESSURE;
    if (d_pspg != NULL && pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        diffusion = 0.;
        if (pd->e[upd->matrix_index[meqn]][meqn] & T_DIFFUSION) {
          diffusion = bf[var]->grad_phi[j][a] - d_mu->P[j] * div_G[a];
          diffusion *= pd->etm[upd->matrix_index[meqn]][meqn][(LOG2_DIFFUSION)];
        }

        d_pspg->P[a][j] = tau_pspg * (diffusion);
      }
    }

    var = MASS_FRACTION;
    if (d_pspg != NULL && pd->v[pg->imtrx][var]) {
      dbl d_rho_t_dC;

      for (w = 0; w < pd->Num_Species_Eqn; w++) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          phi_j = bf[var]->phi[j];

          d_rho_t_dC = ompvf * d_rho->C[w][j];
          if (particle_momentum_on && w == species)
            d_rho_t_dC -= phi_j * rho;

          mass = 0.;
          if ((pd->e[upd->matrix_index[meqn]][meqn] & T_MASS) && (pd->TimeIntegration != STEADY)) {
            mass = d_rho_t_dC / por * v_dot[a];
            mass *= pd->etm[upd->matrix_index[meqn]][meqn][(LOG2_MASS)];
          }

          advection = 0.;
          if (pd->e[upd->matrix_index[meqn]][meqn] & T_ADVECTION) {
            advection = 0.;
            for (p = 0; p < WIM; p++) {
              advection += (v[p] - x_dot[p]) * grad_v[p][a];
            }
            advection *= d_rho_t_dC / por2;
            advection *= pd->etm[upd->matrix_index[meqn]][meqn][(LOG2_ADVECTION)];
          }

          diffusion = 0.;
          if (pd->e[upd->matrix_index[meqn]][meqn] & T_DIFFUSION) {
            diffusion -= d_mu->C[w][j] * div_G[a];
            /*diffusion  -= d_div_tau_p_dy[a][w][j]; */
            diffusion *= pd->etm[upd->matrix_index[meqn]][meqn][(LOG2_DIFFUSION)];
          }

          source = 0.;
          if (pd->e[upd->matrix_index[meqn]][meqn] & T_SOURCE) {
            source -= df->C[a][w][j] * pd->etm[upd->matrix_index[meqn]][meqn][(LOG2_SOURCE)];
          }

          porous = 0.;
          if (pd->e[upd->matrix_index[meqn]][meqn] & T_POROUS_BRINK) {
            porous = v[a] * (d_rho_t_dC * sc * speed / sqrt(per));
            porous *= pd->etm[upd->matrix_index[meqn]][meqn][(LOG2_POROUS_BRINK)];
          }
          d_pspg->C[a][w][j] = tau_pspg * (mass + advection + diffusion + source + porous);
        }
      }
    }

    var = POLYMER_STRESS11;
    if (d_pspg != NULL && pd->v[pg->imtrx][var]) {
      for (mode = 0; mode < vn->modes; mode++) {
        for (b = 0; b < VIM; b++) {
          for (c = 0; c < VIM; c++) {
            var = v_s[mode][b][c];
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              phi_j = bf[var]->phi[j];

              diffusion = 0.;
              if (pd->e[upd->matrix_index[meqn]][meqn] & T_DIFFUSION) {
                diffusion = -(dbl)delta(a, c) * bf[var]->grad_phi[j][b];

                if (pd->CoordinateSystem != CARTESIAN) {
                  for (r = 0; r < VIM; r++) {
                    diffusion -= (dbl)delta(a, c) * phi_j * fv->grad_e[b][r][c];
                  }
                  for (r = 0; r < WIM; r++) {
                    diffusion -= (dbl)delta(a, r) * phi_j * fv->grad_e[c][b][r];
                  }
                }

                diffusion *= pd->etm[upd->matrix_index[meqn]][meqn][(LOG2_DIFFUSION)];
              }

              d_pspg->S[a][mode][b][c][j] = tau_pspg * diffusion;
            }
          }
        }
      }
    }

    var = VELOCITY_GRADIENT11;
    if (d_pspg != NULL && pd->v[pg->imtrx][var]) {
      for (b = 0; b < VIM; b++) {
        for (c = 0; c < VIM; c++) {
          var = v_g[b][c];
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];

            diffusion = 0.;
            if (pd->e[upd->matrix_index[meqn]][meqn] & T_DIFFUSION) {
              diffusion = -(dbl)delta(a, c) * bf[var]->grad_phi[j][b];

              if (pd->CoordinateSystem != CARTESIAN) {
                for (r = 0; r < VIM; r++) {
                  diffusion -= (dbl)delta(a, c) * phi_j * fv->grad_e[b][r][c];
                }
                for (r = 0; r < WIM; r++) {
                  diffusion -= (dbl)delta(a, r) * phi_j * fv->grad_e[c][b][r];
                }
              }
              diffusion *= mu;
              /* diffusion -= d_div_tau_p_dgd0[a][b][c][j]; */
              diffusion *= pd->etm[upd->matrix_index[meqn]][meqn][(LOG2_DIFFUSION)];
            }

            d_pspg->g[a][b][c][j] = tau_pspg * diffusion;
          }
        }
      }
    }
  }
  is_initialized = TRUE;
  return 0;
}

int calc_cont_gls(dbl *cont_gls,
                  CONT_GLS_DEPENDENCE_STRUCT *d_cont_gls,
                  dbl time_value,         /*Current time value, needed for density */
                  const PG_DATA *pg_data) /*Petrov-Galerkin data, needed for tau_cont   */
{
  dbl div_v = fv->div_v;
  dbl tau_cont, d_tau_dmesh[DIM][MDE], d_tau_dv[DIM][MDE];
  dbl advection_etm, advection, Re, div_phi_j_e_b, div_v_dmesh;
  int eqn;
  int var, dim, b, j, p;
  int advection_on = 0;
  static int is_initialized = FALSE;

  // Density terms
  dbl rho;
  DENSITY_DEPENDENCE_STRUCT d_rho_struct;
  DENSITY_DEPENDENCE_STRUCT *d_rho = &d_rho_struct;

  // Petrov-Galerkin values
  dbl h_elem = 0, U, hh_siz, vv;
  const dbl h_elem_avg = pg_data->h_elem_avg;
  const dbl *hsquared = pg_data->hsquared;
  const dbl U_norm = pg_data->U_norm;
  const dbl *v_avg = pg_data->v_avg;
  const dbl mu_avg = pg_data->mu_avg;

  // Initialize and Define
  dim = pd->Num_Dim;

  tau_cont = 0;
  *cont_gls = 0.0;
  eqn = R_PRESSURE;
  advection_on = pd->e[pg->imtrx][eqn] & T_ADVECTION;
  advection_etm = pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];

  if (d_cont_gls == NULL) {
    d_rho = NULL;
  } else if (!is_initialized) {
    memset(d_cont_gls->v, 0, sizeof(dbl) * DIM * MDE);
    memset(d_cont_gls->X, 0, sizeof(dbl) * DIM * MDE);
  }

  /* Calculate stabilization parameter tau_cont
   * From Wall:
   * tau_cont = h_elem*u_norm*max{Re,1}/2;
   * The Reynolds number, Re, is defined as in calc_pspg
   */
  if (Cont_GLS == 1) {
    h_elem = h_elem_avg;
    U = U_norm;
  } else {
    hh_siz = 0.0;
    vv = 0.0;
    for (p = 0; p < dim; p++) {
      hh_siz += hsquared[p] / ((dbl)dim);
      vv += v_avg[p] * v_avg[p];
    }
    h_elem = sqrt(hh_siz);
    U = sqrt(vv);
  }

  tau_cont = U * h_elem / 2.0;

  // We need density for the Reynolds number
  rho = density(d_rho, time_value);
  Re = rho * U * h_elem / (2.0 * mu_avg);
  if (Re > 1.0) {
    tau_cont *= Re;
  }

  // d_tau terms for Jacobian
  if (d_cont_gls != NULL && pd->v[pg->imtrx][VELOCITY1]) {
    for (b = 0; b < dim; b++) {
      var = VELOCITY1 + b;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          if (Cont_GLS == 1 || U == 0) {
            d_tau_dv[b][j] = 0.0;
          } else if (Re > 1.0) {
            d_tau_dv[b][j] =
                rho / (2.0 * mu_avg) * h_elem * h_elem * v_avg[b] * pg_data->dv_dnode[b][j];
          } else {
            d_tau_dv[b][j] = h_elem / (2.0 * U) * v_avg[b] * pg_data->dv_dnode[b][j];
          }
        }
      }
    }
  }

  if (d_cont_gls != NULL && pd->v[pg->imtrx][MESH_DISPLACEMENT1]) {
    for (b = 0; b < dim; b++) {
      var = MESH_DISPLACEMENT1 + b;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          if (Cont_GLS == 1 || h_elem == 0) {
            d_tau_dmesh[b][j] = 0.0;
          } else if (Re > 1.0) {
            d_tau_dmesh[b][j] = rho / (2.0 * mu_avg) * U * U * pg_data->hhv[b][b] *
                                pg_data->dhv_dxnode[b][j] / ((dbl)dim);
          } else {
            d_tau_dmesh[b][j] =
                U / (2.0 * h_elem) * pg_data->hhv[b][b] * pg_data->dhv_dxnode[b][j] / ((dbl)dim);
          }
        }
      }
    }
  }

  /*
   * Calculate residual
   * This term refers to the standard del dot v .
   */
  advection = 0.0;
  if (advection_on) {
    if (pd->v[pg->imtrx][VELOCITY1]) {
      advection = div_v;
      advection *= advection_etm;
    }
  }

  /*Multiply tau_cont by residual resulting in cont_gls
   *This eventually then gets combined with the stabilization functional
   *in assemble_momentum to form the GLS stabilization for continuity
   */
  *cont_gls = tau_cont * advection;

  // Determine Jacobian terms

  // J_v
  for (b = 0; b < WIM; b++) {
    var = VELOCITY1 + b;
    if (pd->v[pg->imtrx][var] && d_cont_gls != NULL) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        advection = 0.;
        if (advection_on) {
          div_phi_j_e_b = 0.;
          for (p = 0; p < VIM; p++) {
            div_phi_j_e_b += bf[var]->grad_phi_e[j][b][p][p];
          }
          advection = div_phi_j_e_b * advection_etm;
        }
        d_cont_gls->v[b][j] = tau_cont * advection + d_tau_dv[b][j] * div_v * advection_etm;
      }
    }
  }

  // J_d
  for (b = 0; b < dim; b++) {
    var = MESH_DISPLACEMENT1 + b;
    if (pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        advection = 0.0;
        if (advection_on) {
          if (pd->v[pg->imtrx][VELOCITY1] && d_cont_gls != NULL) {
            div_v_dmesh = fv->d_div_v_dmesh[b][j];
            advection += div_v_dmesh;
          }
        }
        advection *= advection_etm;
        d_cont_gls->X[b][j] = tau_cont * advection + d_tau_dmesh[b][j] * div_v * advection_etm;
      }
    }
  }

  is_initialized = TRUE;
  return 0;
}
