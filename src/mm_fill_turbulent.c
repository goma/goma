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

#include "mm_fill_stabilization.h"
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
 *  Kessels, P. C. J. "Finite element discretization of the Spalart-Allmaras
 *  turbulence model." (2016).
 *
 *  Spalart, Philippe, and Steven Allmaras. "A one-equation turbulence model for
 *  aerodynamic flows." 30th aerospace sciences meeting and exhibit. 1992.
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
 * Modified:    June 2023 Weston Ortiz
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

  eqn = EDDY_NU;
  double d_area = fv->wt * bf[eqn]->detJ * fv->h3;

  /* Get Eddy viscosity at Gauss point */
  double mu_e = fv->eddy_nu;

  int negative_sa = false;
  // Previous workaround, see comment for negative_Se below
  //
  //   int transient_run = FALSE;
  //   if (pd->TimeIntegration != STEADY) {
  //     transient_run = true;
  //   }
  //   // Use old values for equation switching for transient runs
  //   // Seems to work reasonably well.
  //   if (transient_run && (fv_old->eddy_nu < 0)) {
  //     negative_sa = true;
  //   } else if (!transient_run && (mu_e < 0)) {
  //     // Kris thinks it might work with switching equations in steady state
  //     negative_sa = true;
  //   }
  if (mu_e < 0) {
    negative_sa = true;
  }

  /* Get fluid viscosity */
  double mu_newt = mp->viscosity;

  /* Rate of rotation tensor  */
  double omega[DIM][DIM];
  double omega_old[DIM][DIM];
  for (a = 0; a < VIM; a++) {
    for (b = 0; b < VIM; b++) {
      omega[a][b] = (fv->grad_v[a][b] - fv->grad_v[b][a]);
      omega_old[a][b] = (fv_old->grad_v[a][b] - fv_old->grad_v[b][a]);
    }
  }

  /* Vorticity */
  double S = 0.0;
  double S_old = 0;
  double dS_dvelo[DIM][MDE];
  double dS_dmesh[DIM][MDE];
  calc_sa_S(&S, omega, dS_dvelo, dS_dmesh);
  calc_sa_S(&S_old, omega_old, NULL, NULL);

  double d = fv->wall_distance;
  /* Get distance from nearest wall */
  if (d < 1.0e-6)
    d = 1.0e-6;

  /* Model coefficients (constants) */
  double cb1 = 0.1355;
  double cb2 = 0.622;
  double cv1 = 7.1;
  double cv2 = 0.7;
  double cv3 = 0.9;
  double sigma = (2.0 / 3.0);
  double cw2 = 0.3;
  double cw3 = 2.0;
  double cn1 = 16;
  double kappa = 0.41;
  double cw1 = (cb1 / kappa / kappa) + (1.0 + cb2) / sigma;

  /* More model coefficients (depends on mu_e) */
  double chi = mu_e / mu_newt;
  double fv1 = pow(chi, 3) / (pow(chi, 3) + pow(cv1, 3));
  double fv2 = 1.0 - (chi) / (1.0 + chi * fv1);
  double fn = 1.0;
  if (negative_sa) {
    fn = (cn1 + pow(chi, 3.0)) / (cn1 - pow(chi, 3));
  }
  double Sbar = (mu_e * fv2) / (kappa * kappa * d * d);
  int negative_Se = false;
  // I tried to use the old values for equation switching for transient runs but
  // I end up getting floating point errors. because of Se going very far below
  // zero I'm trying to use current values instead and hope that Newton's method
  // will converge with the switching equations
  // previously:
  // . double Sbar_old = (fv_old->eddy_nu * fv2) / (kappa * kappa * d * d);
  //   if (transient_run && (Sbar_old < -cv2 * S_old)) {
  //     negative_Se = true;
  //   } else if (!transient_run && (Sbar < -cv2 * S)) {
  // .   negative_Se = true;
  //   }
  if (Sbar < -cv2 * S) {
    negative_Se = true;
  }
  double S_e = S + Sbar;
  if (negative_Se) {
    S_e = S + S * (cv2 * cv2 * S + cv3 * Sbar) / ((cv3 - 2 * cv2) * S - Sbar);
  }
  double r_max = 10.0;
  double r = 0.0;
  if (fabs(S_e) > 1.0e-6) {
    r = mu_e / (kappa * kappa * d * d * S_e);
  } else {
    r = r_max;
  }
  if (r >= r_max) {
    r = r_max;
  }
  // Arbitrary limit to avoid floating point errors should only hit this when
  // S_e is very small and either mu_e or S_e are negative.  Which means we are
  // already trying to alleviate the issue.
  if (r < -100) {
    r = -100;
  }
  double g = r + cw2 * (pow(r, 6) - r);
  double fw_inside = (1.0 + pow(cw3, 6)) / (pow(g, 6) + pow(cw3, 6));
  double fw = g * pow(fw_inside, (1.0 / 6.0));

  /* Model coefficients sensitivity w.r.t. mu_e*/
  double dchi_dmu_e = 1.0 / mu_newt;
  double dfv1_dchi = 3 * pow(cv1, 3) * pow(chi, 2) / (pow((pow(chi, 3) + pow(cv1, 3)), 2));
  double dfv1_dmu_e = dfv1_dchi * dchi_dmu_e;
  double dfv2_dmu_e = (chi * chi * dfv1_dmu_e - dchi_dmu_e) / (pow(1 + fv1 * chi, 2));
  double dfn_dmu_e = 0.0;
  if (negative_sa) {
    dfn_dmu_e = 6.0 * (cn1 * chi * chi * dchi_dmu_e) / (pow(cn1 - pow(chi, 3.0), 2.0));
  }
  double dSbar_dmu_e = (fv2 + mu_e * dfv2_dmu_e) / (kappa * kappa * d * d);
  double dS_e_dmu_e = dSbar_dmu_e;
  if (negative_Se) {
    dS_e_dmu_e =
        (pow(cv2 - cv3, 2.0) * S * S * dSbar_dmu_e) / (pow(2 * cv2 * S - cv3 * S + Sbar, 2.0));
  }
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
  double dr_dS = 0;
  double dg_dS = 0;
  double dfw_dS = 0;
  if (fabs(S_e) > 1.0e-6) {
    dr_dS = -(mu_e * kappa * kappa * d * d) / (kappa * kappa * d * d * S_e) /
            (kappa * kappa * d * d * S_e);
  }
  if (r == r_max)
    dr_dS = 0.0;
  if (r < -100)
    dr_dS = 0.0;
  dg_dS = dg_dr * dr_dS;
  dfw_dS = dfw_dg * dg_dS;
  dbl dS_e_dS = 1.0;
  if (negative_Se) {
    S_e = S + S * (cv2 * cv2 * S + cv3 * Sbar) / ((cv3 - 2 * cv2) * S - Sbar);
    dS_e_dS += (cv2 * cv2 * S) / ((cv3 - 2 * cv2) * S - Sbar) +
               (cv2 * cv2 * S + cv3 * Sbar) / ((cv3 - 2 * cv2) * S - Sbar) -
               (cv3 - 2 * cv2) * S * (cv2 * cv2 * S + cv3 * Sbar) /
                   (pow(((cv3 - 2 * cv2) * S - Sbar), 2.0));
  }

  dbl supg = 1.;
  SUPG_terms supg_terms;
  if (mp->Mwt_funcModel == GALERKIN) {
    supg = 0.;
  } else if (mp->Mwt_funcModel == SUPG || mp->Mwt_funcModel == SUPG_GP ||
             mp->Mwt_funcModel == SUPG_SHAKIB) {
    supg = mp->Mwt_func;
    supg_tau_shakib(&supg_terms, pd->Num_Dim, dt, mu_newt, EDDY_NU);
  }

  /*
   * Residuals_________________________________________________________________
   */

  if (af->Assemble_Residual) {
    /*
     * Assemble residual for eddy viscosity
     */
    eqn = EDDY_NU;
    peqn = upd->ep[pg->imtrx][eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      dbl wt_func = bf[eqn]->phi[i];

      if (supg > 0) {
        if (supg != 0.0) {
          for (int p = 0; p < VIM; p++) {
            wt_func += supg * supg_terms.supg_tau * fv->v[p] * bf[eqn]->grad_phi[i][p];
          }
        }
      }

      /* Assemble mass term */
      mass = 0.0;
      if (pd->TimeIntegration != STEADY) {
        if (pd->e[pg->imtrx][eqn] & T_MASS) {
          mass += fv_dot->eddy_nu * wt_func * d_area;
          mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
        }
      }

      /* Assemble advection term */
      adv = 0;
      for (int p = 0; p < VIM; p++) {
        adv += fv->v[p] * fv->grad_eddy_nu[p];
      }
      adv *= wt_func * d_area;
      adv *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];

      double neg_c = 1.0;
      if (negative_sa) {
        neg_c = -1.0;
      }

      /* Assemble source terms */
      src_1 = cb1 * S_e * mu_e;
      src_2 = neg_c * cw1 * fw * (mu_e * mu_e) / (d * d);
      src = -src_1 + src_2;
      src *= wt_func * d_area;
      src *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

      /* Assemble diffusion terms */
      diff_1 = 0.0;
      diff_2 = 0.0;
      for (int p = 0; p < VIM; p++) {
        diff_1 += bf[eqn]->grad_phi[i][p] * (mu_newt + mu_e * fn) * fv->grad_eddy_nu[p];
        diff_2 += wt_func * cb2 * fv->grad_eddy_nu[p] * fv->grad_eddy_nu[p];
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
    eqn = EDDY_NU;
    peqn = upd->ep[pg->imtrx][eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {

      dbl wt_func = bf[eqn]->phi[i];

      dbl r_mass = 0;
      dbl r_adv = 0;
      dbl r_src = 0;
      dbl r_diff = 0;
      if (supg > 0) {
        if (supg != 0.0) {
          for (int p = 0; p < VIM; p++) {
            wt_func += supg * supg_terms.supg_tau * fv->v[p] * bf[eqn]->grad_phi[i][p];
          }
        }

        /* Assemble mass term */
        r_mass = 0.0;
        if (pd->TimeIntegration != STEADY) {
          if (pd->e[pg->imtrx][eqn] & T_MASS) {
            r_mass += fv_dot->eddy_nu;
            mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
          }
        }

        /* Assemble advection term */
        r_adv = 0;
        for (int p = 0; p < VIM; p++) {
          r_adv += fv->v[p] * fv->grad_eddy_nu[p];
        }
        r_adv *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];

        double neg_c = 1.0;
        if (negative_sa) {
          neg_c = -1.0;
        }

        /* Assemble source terms */
        src_1 = cb1 * S_e * mu_e;
        src_2 = neg_c * cw1 * fw * (mu_e * mu_e) / (d * d);
        r_src = -src_1 + src_2;
        r_src *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

        /* Assemble diffusion terms */
        diff_1 = 0.0;
        diff_2 = 0.0;
        for (int p = 0; p < VIM; p++) {
          diff_2 += wt_func * cb2 * fv->grad_eddy_nu[p] * fv->grad_eddy_nu[p];
        }
        r_diff = (1.0 / sigma) * (diff_1 - diff_2);
        r_diff *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
      }

      /* Sensitivity w.r.t. eddy viscosity */
      var = EDDY_NU;
      if (pdv[var]) {
        pvar = upd->vp[pg->imtrx][var];

        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {

          /* Assemble mass term */
          mass = 0.0;
          if (pd->TimeIntegration != STEADY) {
            if (pd->e[pg->imtrx][eqn] & T_MASS) {
              mass += (1.0 + 2.0 * tt) / dt * bf[eqn]->phi[j] * wt_func * d_area;
              mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
            }
          }

          /* Assemble advection term */
          adv = 0.0;
          for (int p = 0; p < VIM; p++) {
            adv += fv->v[p] * bf[eqn]->grad_phi[j][p];
          }
          adv *= wt_func * d_area;
          adv *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];

          double neg_c = 1.0;
          if (negative_sa) {
            neg_c = -1.0;
          }
          /* Assemble source term */
          src_1 = cb1 * (dS_e_dmu_e * bf[eqn]->phi[j] * mu_e + S_e * bf[eqn]->phi[j]);
          src_2 = neg_c *
                  ((cw1 / d / d) * bf[var]->phi[j] * (dfw_dmu_e * mu_e * mu_e + fw * 2.0 * mu_e));
          src = -src_1 + src_2;
          src *= wt_func * d_area;
          src *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

          /* Assemble diffusion terms */
          diff_1 = 0.0;
          diff_2 = 0.0;
          for (int p = 0; p < VIM; p++) {
            diff_1 += bf[eqn]->grad_phi[i][p] *
                      ((fn + mu_e * dfn_dmu_e) * bf[eqn]->phi[j] * fv->grad_eddy_nu[p] +
                       (mu_newt + mu_e * fn) * bf[eqn]->grad_phi[j][p]);
            diff_2 += wt_func * cb2 * 2.0 * fv->grad_eddy_nu[p] * bf[var]->grad_phi[j][p];
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
            if (supg > 0) {
              dbl d_wt_func =
                  supg * bf[var]->phi[j] * supg_terms.supg_tau * bf[eqn]->grad_phi[i][b];
              for (int p = 0; p < VIM; p++) {
                d_wt_func +=
                    supg * fv->v[p] * supg_terms.d_supg_tau_dv[p][j] * bf[eqn]->grad_phi[i][p];
              }
              lec->J[LEC_J_INDEX(peqn, pvar, i, j)] +=
                  (r_mass + r_adv + r_src + r_diff) * d_area * d_wt_func;
            }

            /* Assemble advection term */
            adv = bf[var]->phi[j] * fv->grad_eddy_nu[b];
            adv *= wt_func * d_area;
            adv *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];

            double neg_c = 1.0;
            if (negative_sa) {
              neg_c = -1.0;
            }
            /* Assemble source term */
            src_1 = cb1 * dS_e_dS * dS_dvelo[b][j] * mu_e;
            src_2 = neg_c * cw1 * dfw_dS * dS_dvelo[b][j] * (mu_e / d) * (mu_e / d);
            src = -src_1 + src_2;
            src *= wt_func * d_area;
            src *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

            lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += adv + src;

          } /* End of loop over j */
        }   /* End of if the variale is active */
      }     /* End of loop over velocity components */

      /* Sensistivity w.r.t. mesh */
      for (b = 0; b < pd->Num_Dim; b++) {
        var = MESH_DISPLACEMENT1 + b;
        if (pdv[var]) {
          pvar = upd->vp[pg->imtrx][var];
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            double d_area_dmesh = +fv->wt * bf[eqn]->d_det_J_dm[b][j] * fv->h3 +
                                  fv->wt * bf[eqn]->detJ * fv->dh3dmesh[b][j];
            if (supg > 0) {
              dbl d_wt_func = 0;
              for (int p = 0; p < VIM; p++) {
                d_wt_func +=
                    supg * fv->v[p] * supg_terms.d_supg_tau_dX[p][j] * bf[eqn]->grad_phi[i][p];
                d_wt_func +=
                    supg * fv->v[p] * supg_terms.supg_tau * bf[eqn]->d_grad_phi_dmesh[i][p][b][j];
              }
              lec->J[LEC_J_INDEX(peqn, pvar, i, j)] +=
                  (r_mass + r_adv + r_src + r_diff) * d_area * d_wt_func +
                  (r_mass + r_adv + r_src + r_diff) * d_area_dmesh * wt_func;
            }

            /* Assemble mass term */
            mass = 0.0;
            if (pd->TimeIntegration != STEADY) {
              if (pd->e[pg->imtrx][eqn] & T_MASS) {
                mass += fv_dot->eddy_nu * wt_func * d_area_dmesh;
                mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
              }
            }

            /* Assemble advection term */
            adv = 0;
            for (int p = 0; p < VIM; p++) {
              adv += fv->v[p] * fv->grad_eddy_nu[p];
            }
            dbl d_adv = 0;
            for (int p = 0; p < VIM; p++) {
              d_adv += fv->v[p] * fv->d_grad_eddy_nu_dmesh[p][b][j];
            }
            adv *= wt_func * d_area_dmesh;
            adv += d_adv * wt_func * d_area;
            adv *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];

            double neg_c = 1.0;
            if (negative_sa) {
              neg_c = -1.0;
            }

            /* Assemble source terms */
            dbl d_src_1 = cb1 * dS_e_dS * dS_dmesh[b][j] * mu_e;
            dbl d_src_2 = neg_c * cw1 * dfw_dS * dS_dmesh[b][j] * (mu_e * mu_e) / (d * d);
            src_1 = cb1 * S_e * mu_e;
            src_2 = neg_c * cw1 * fw * (mu_e * mu_e) / (d * d);
            src = -src_1 + src_2;
            dbl d_src = -d_src_1 + d_src_2;
            src *= wt_func * d_area_dmesh;
            src += d_src * wt_func * d_area;
            src *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

            /* Assemble diffusion terms */
            dbl d_diff_1 = 0;
            dbl d_diff_2 = 0;
            diff_1 = 0.0;
            diff_2 = 0.0;
            for (int p = 0; p < VIM; p++) {
              diff_1 += bf[eqn]->grad_phi[i][p] * (mu_newt + mu_e * fn) * fv->grad_eddy_nu[p];
              diff_2 += wt_func * cb2 * fv->grad_eddy_nu[p] * fv->grad_eddy_nu[p];

              d_diff_1 += bf[eqn]->d_grad_phi_dmesh[i][p][b][j] * (mu_newt + mu_e * fn) *
                          fv->grad_eddy_nu[p];
              d_diff_1 += bf[eqn]->grad_phi[i][p] * (mu_newt + mu_e * fn) *
                          fv->d_grad_eddy_nu_dmesh[p][b][j];
              d_diff_2 +=
                  wt_func * cb2 * 2.0 * fv->d_grad_eddy_nu_dmesh[p][b][j] * fv->grad_eddy_nu[p];
            }
            diff = (1.0 / sigma) * (diff_1 - diff_2);
            dbl d_diff = (1.0 / sigma) * (d_diff_1 - d_diff_2);
            diff *= d_area_dmesh;
            diff += d_diff * d_area;
            diff *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];

            lec->R[LEC_R_INDEX(peqn, i)] += mass + adv + src + diff;
          } /* End of loop over j */
        }
      }

    } /* End of loop over i */
  }   /* End of if assemble Jacobian */
  return (status);
}
