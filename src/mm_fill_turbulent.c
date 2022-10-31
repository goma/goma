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

#include <stdio.h>

/* GOMA include files */
#define GOMA_MM_FILL_TURBULENT_C
#include "el_elm.h"
#include "mm_as.h"
#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "mm_fill_terms.h"
#include "mm_fill_turbulent.h"
#include "mm_fill_util.h"
#include "mm_mp.h"
#include "mm_viscosity.h"
#include "rf_fem.h"
#include "rf_fem_const.h"
#include "std.h"

/*  _______________________________________________________________________  */

static int calc_sa_S(dbl *S,              /* strain rate invariant */
                     dbl omega[DIM][DIM], /* strain rate tensor */
                     dbl d_S_dv[DIM][MDE],
                     dbl d_S_dmesh[DIM][MDE]) {
  int mdofs = 0;
  int p, q, a, b;
  int vdofs, i, j, v;

  dbl grad_phi_e_omega[MDE][DIM][DIM][DIM]; /* transpose of grad(phi_i ea) tensor
                                             + grad(phi_i ea) tensor */
  dbl d_omega_dmesh[DIM][DIM][DIM][MDE];    /* d/dmesh(grad_v)T */

  int status = 1;

  /* Zero out sensitivities */

  if (d_S_dv != NULL)
    memset(d_S_dv, 0, sizeof(double) * DIM * MDE);
  if (d_S_dmesh != NULL)
    memset(d_S_dmesh, 0, sizeof(double) * DIM * MDE);

  *S = 0.;
  /* get gamma_dot invariant for viscosity calculations */
  for (a = 0; a < VIM; a++) {
    for (b = 0; b < VIM; b++) {
      *S += omega[a][b] * omega[a][b];
    }
  }

  *S = sqrt(0.5 * (*S));

  /* get stuff for Jacobian entries */
  v = VELOCITY1;
  vdofs = ei[pg->imtrx]->dof[v];

  if (d_S_dmesh != NULL || d_S_dv != NULL) {
    if (pd->v[pg->imtrx][R_MESH1]) {
      mdofs = ei[pg->imtrx]->dof[R_MESH1];
    }

    for (p = 0; p < VIM; p++) {
      for (q = 0; q < VIM; q++) {
        for (a = 0; a < VIM; a++) {
          for (i = 0; i < vdofs; i++) {
            grad_phi_e_omega[i][a][p][q] =
                bf[v]->grad_phi_e[i][a][p][q] - bf[v]->grad_phi_e[i][a][q][p];
          }
        }
      }
    }
  }

  /*
   * d( gamma_dot )/dmesh
   */

  if (pd->v[pg->imtrx][R_MESH1] && d_S_dmesh != NULL) {

    for (p = 0; p < VIM; p++) {
      for (q = 0; q < VIM; q++) {
        for (b = 0; b < VIM; b++) {
          for (j = 0; j < mdofs; j++) {

            d_omega_dmesh[p][q][b][j] =
                fv->d_grad_v_dmesh[p][q][b][j] - fv->d_grad_v_dmesh[q][p][b][j];
          }
        }
      }
    }

    /*
     * d( gammadot )/dmesh
     */

    if (*S != 0.) {
      for (b = 0; b < VIM; b++) {
        for (j = 0; j < mdofs; j++) {
          d_S_dmesh[b][j] = 0.;
          for (p = 0; p < VIM; p++) {
            for (q = 0; q < VIM; q++) {
              d_S_dmesh[b][j] += 0.5 * d_omega_dmesh[p][q][b][j] * omega[p][q] / *S;
            }
          }
        }
      }
    }
  }

  /*
   * d( gammadot )/dv
   */

  if (*S != 0. && d_S_dv != NULL) {
    for (a = 0; a < VIM; a++) {
      for (i = 0; i < vdofs; i++) {
        d_S_dv[a][i] = 0.;
        for (p = 0; p < VIM; p++) {
          for (q = 0; q < VIM; q++) {
            d_S_dv[a][i] += 0.5 * grad_phi_e_omega[i][a][p][q] * omega[p][q] / *S;
          }
        }
      }
    }
  }
  return (status);
}

/* assemble_spalart_allmaras -- assemble terms (Residual & Jacobian) for conservation
 *                              of eddy viscosity for Spalart Allmaras turbulent flow model
 *
 * in:
 *      time value
 *      theta (time stepping parameter, 0 for BE, 0.5 for CN)
 *      time step size
 *      Streamline Upwind Petrov Galerkin (PG) data structure
 *
 * out:
 *      lec -- gets loaded up with local contributions to resid, Jacobian
 *
 * Created:     August 2022 kristianto.tjiptowidjojo@averydennison.com
 *
 */
int assemble_spalart_allmaras(dbl time_value, /* current time */
                              dbl tt,         /* parameter to vary time integration from
                                                 explicit (tt = 1) to implicit (tt = 0)    */
                              dbl dt,         /* current time step size                    */
                              const PG_DATA *pg_data) {

  //! WIM is the length of the velocity vector
  int i, j, a, b;
  int eqn, var, peqn, pvar;
  int *pdv = pd->v[pg->imtrx];

  double mass, adv, src, diff;
  double src_1, src_2, diff_1, diff_2;

  int status = 0;

  eqn = EDDY_MU;
  double d_area = fv->wt * bf[eqn]->detJ * fv->h3;

  /* Get Eddy viscosity at Gauss point */
  double mu_e = fv->eddy_mu;

  /* Get fluid viscosity */
  // double mu;
  // VISCOSITY_DEPENDENCE_STRUCT d_mu_struct;
  // VISCOSITY_DEPENDENCE_STRUCT *d_mu = &d_mu_struct;
  // double gamma[DIM][DIM];
  // for (a = 0; a < VIM; a++) {
  //   for (b = 0; b < VIM; b++) {
  //     gamma[a][b] = fv->grad_v[a][b] + fv->grad_v[b][a];
  //   }
  // }
  // mu = viscosity(gn, gamma, d_mu);
  double mu_newt = mp->viscosity;

  /* Rate of rotation tensor  */
  double omega[DIM][DIM];
  for (a = 0; a < VIM; a++) {
    for (b = 0; b < VIM; b++) {
      omega[a][b] = (fv->grad_v[a][b] - fv->grad_v[b][a]);
    }
  }

  /* Vorticity */
  double S = 0.0;
  double dS_dvelo[DIM][MDE];
  double dS_dmesh[DIM][MDE];
  calc_sa_S(&S, omega, dS_dvelo, dS_dmesh);

  /* Get distance from nearest wall - use external field for now */
  int i_d_wall = mp->dist_wall_ext_field_index;
  double d = fv->external_field[i_d_wall];
  if (d < 1.0e-6)
    d = 1.0e-6;

  /* Model coefficients (constants) */
  double cb1 = 0.1355;
  double cb2 = 0.622;
  double cv1 = 7.1;
  double sigma = (2.0 / 3.0);
  double cw2 = 0.3;
  double cw3 = 2.0;
  double kappa = 0.41;
  double cw1 = (cb1 / kappa / kappa) + (1.0 + cb2) / sigma;

  /* More model coefficients (depends on mu_e) */
  double chi = mu_e / mu_newt;
  double fv1 = pow(chi, 3) / (pow(chi, 3) + pow(cv1, 3));
  double fv2 = 1.0 - (chi) / (1.0 + chi * fv1);
  double S_e = S + (mu_e * fv2) / (kappa * kappa * d * d);
  double r_max = 10.0;
  double r = 0.0;
  if (fabs(S_e) > 1.0e-6) {
    r = mu_e / (kappa * kappa * d * d * S_e);
  } else {
    r = r_max;
  }
  if (r > r_max) {
    r = r_max;
  }
  double g = r + cw2 * (pow(r, 6) - r);
  double fw_inside = (1.0 + pow(cw3, 6)) / (pow(g, 6) + pow(cw3, 6));
  double fw = g * pow(fw_inside, (1.0 / 6.0));

  /* Model coefficients sensitivity w.r.t. mu_e*/
  double dchi_dmu_e = 1.0 / mu_newt;
  double dfv1_dchi = 3 * pow(cv1, 3) * pow(chi, 2) / (pow((pow(chi, 3) + pow(cv1, 3)), 2));
  double dfv1_dmu_e = dfv1_dchi * dchi_dmu_e;
  double dfv2_dmu_e = (chi * chi * dfv1_dmu_e - dchi_dmu_e) / (pow(1 + fv1 * chi, 2));
  double dS_e_dmu_e = (fv2 + mu_e * dfv2_dmu_e) / (kappa * kappa * d * d);
  double dr_dmu_e = 0.0;
  if (r < r_max) {
    dr_dmu_e = 1.0 / (kappa * kappa * d * d * S_e) - (mu_e * kappa * kappa * d * d * dS_e_dmu_e) /
                                                         (kappa * kappa * d * d * S_e) /
                                                         (kappa * kappa * d * d * S_e);
  }
  double dg_dr = 1.0 + cw2 * (6.0 * pow(r, 5) - 1.0);
  double dg_dmu_e = dg_dr * dr_dmu_e;
  double dfw_inside_dg = -6. * (pow(g, 5) * (pow(cw3, 6) + 1)) / pow((pow(cw3, 6) + pow(g, 6)), 2);
  double dfw_dg = pow(fw_inside, (1.0 / 6.0)) +
                  g * dfw_inside_dg * (1.0 / 6.0) * pow(fw_inside, (1.0 / 6.0) - 1.0);
  double dfw_dmu_e = dfw_dg * dg_dmu_e;

  /* Model coefficients sensitivity w.r.t. velocity*/
  double dr_dvelo[DIM][MDE];
  double dg_dvelo[DIM][MDE];
  double dfw_dvelo[DIM][MDE];
  for (b = 0; b < VIM; b++) {
    for (j = 0; j < ei[pg->imtrx]->dof[VELOCITY1]; j++) {
      dr_dvelo[b][j] = -(mu_e * kappa * kappa * d * d * dS_dvelo[b][j]) /
                       (kappa * kappa * d * d * S_e) / (kappa * kappa * d * d * S_e);
      if (r == r_max)
        dr_dvelo[b][j] = 0.0;
      dg_dvelo[b][j] = dg_dr * dr_dvelo[b][j];
      dfw_dvelo[b][j] = dfw_dg * dg_dvelo[b][j];
    }
  }

  /*
   * Residuals_________________________________________________________________
   */

  if (af->Assemble_Residual) {
    /*
     * Assemble residual for eddy viscosity
     */
    eqn = EDDY_MU;
    peqn = upd->ep[pg->imtrx][eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {

      /* Assemble mass term */
      mass = 0.0;
      if (pd->TimeIntegration != STEADY) {
        if (pd->e[pg->imtrx][eqn] & T_MASS) {
          mass += fv_dot->eddy_mu * bf[eqn]->phi[i] * d_area;
          mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
        }
      }

      /* Assemble advection term */
      adv = 0;
      for (int p = 0; p < VIM; p++) {
        adv += fv->v[p] * fv->grad_eddy_mu[p];
      }
      adv *= bf[eqn]->phi[i] * d_area;
      adv *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];

      /* Assemble source terms */
      src_1 = cb1 * S_e * mu_e;
      src_2 = cw1 * fw * (mu_e * mu_e) / (d * d);
      src = -src_1 + src_2;
      src *= bf[eqn]->phi[i] * d_area;
      src *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

      /* Assemble diffusion terms */
      diff_1 = 0.0;
      diff_2 = 0.0;
      for (int p = 0; p < VIM; p++) {
        diff_1 += bf[eqn]->grad_phi[i][p] * (mu_newt + mu_e) * fv->grad_eddy_mu[p];
        diff_2 += bf[eqn]->phi[i] * cb2 * fv->grad_eddy_mu[p] * fv->grad_eddy_mu[p];
      }
      diff = (1.0 / sigma) * (diff_1 - diff_2);
      diff *= d_area;
      diff *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];

      lec->R[LEC_R_INDEX(peqn, i)] += mass + adv + src + diff;
    } /* end of for (i=0,ei[pg->imtrx]->dofs...) */
  }   /* end of if assemble residual */

  /*
   * Jacobian terms...
   */

  if (af->Assemble_Jacobian) {
    eqn = EDDY_MU;
    peqn = upd->ep[pg->imtrx][eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {

      /* Sensitivity w.r.t. eddy viscosity */
      var = EDDY_MU;
      if (pdv[var]) {
        pvar = upd->vp[pg->imtrx][var];

        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {

          /* Assemble mass term */
          mass = 0.0;
          if (pd->TimeIntegration != STEADY) {
            if (pd->e[pg->imtrx][eqn] & T_MASS) {
              mass += (1.0 + 2.0 * tt) / dt * bf[eqn]->phi[i] * d_area;
              mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
            }
          }

          /* Assemble advection term */
          adv = 0.0;
          for (int p = 0; p < VIM; p++) {
            adv += fv->v[p] * bf[eqn]->grad_phi[j][p];
          }
          adv *= bf[eqn]->phi[i] * d_area;
          adv *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];

          /* Assemble source term */
          src_1 = cb1 * (dS_e_dmu_e * bf[eqn]->phi[j] * mu_e + S_e * bf[eqn]->phi[j]);
          src_2 = (cw1 / d / d) * bf[var]->phi[j] * (dfw_dmu_e * mu_e * mu_e + fw * 2.0 * mu_e);
          src = -src_1 + src_2;
          src *= bf[eqn]->phi[i] * d_area;
          src *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

          /* Assemble diffusion terms */
          diff_1 = 0.0;
          diff_2 = 0.0;
          for (int p = 0; p < VIM; p++) {
            diff_1 += bf[eqn]->grad_phi[i][p] * (bf[eqn]->phi[j] * fv->grad_eddy_mu[p] +
                                                 (mu_newt + mu_e) * bf[eqn]->grad_phi[j][p]);
            diff_2 += bf[eqn]->phi[i] * cb2 * 2.0 * fv->grad_eddy_mu[p] * bf[var]->grad_phi[j][p];
          }
          diff = (1.0 / sigma) * (diff_1 - diff_2);
          diff *= d_area;
          diff *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];

          lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += mass + adv + src + diff;
        } /* End of loop over j */
      }   /* End of if the variable is active */

      /* Sensitivity w.r.t. velocity */
      for (b = 0; b < VIM; b++) {
        var = VELOCITY1 + b;
        if (pdv[var]) {
          pvar = upd->vp[pg->imtrx][var];

          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {

            /* Assemble advection term */
            adv = bf[var]->phi[j] * fv->grad_eddy_mu[b];
            adv *= bf[eqn]->phi[i] * d_area;
            adv *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];

            /* Assemble source term */
            src_1 = cb1 * dS_dvelo[b][j] * mu_e;
            src_2 = cw1 * dfw_dvelo[b][j] * (mu_e / d) * (mu_e / d);
            src = -src_1 + src_2;
            src *= bf[eqn]->phi[i] * d_area;
            src *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

            lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += adv + src;

          } /* End of loop over j */
        }   /* End of if the variale is active */
      }     /* End of loop over velocity components */

    } /* End of loop over i */
  }   /* End of if assemble Jacobian */
  return (status);
}
