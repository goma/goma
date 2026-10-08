
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

#include <cstddef>
#ifdef GOMA_ENABLE_SACADO

#include "ad/stabilization.h"
#include "ad/structs.h"
#include "ad/turbulence.h"
#include "ad/viscosity.h"
#include "density.h"
#include "el_elm.h"
#include "mm_as.h"
#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "mm_eh.h"
#include "mm_fill_energy.h"
#include "mm_mp.h"
#include "mm_mp_structs.h"
#include "rf_fem.h"
#include "rf_fem_const.h"
#include "std.h"

/*  _______________________________________________________________________  */

static int calc_vort_mag(ADType &vort_mag, ADType omega[DIM][DIM]) {

  vort_mag = 0.;
  /* get gamma_dot invariant for viscosity calculations */
  for (int a = 0; a < VIM; a++) {
    for (int b = 0; b < VIM; b++) {
      vort_mag += omega[a][b] * omega[a][b];
    }
  }

  vort_mag = sqrt(0.5 * (vort_mag) + 1e-14);
  return 0;
}

ADType ad_dcdd(int dim, int eqn, const ADType grad_U[DIM], ADType gs_inner_dot[DIM]) {
  ADType tmp = 0.0;
  ADType s[DIM] = {0.0};
  ADType r[DIM] = {0.0};
  for (int w = 0; w < dim; w++) {
    tmp += (ad_fv->v[w] - ad_fv->x_dot[w]) * (ad_fv->v[w] - ad_fv->x_dot[w]);
  }
  tmp = 1.0 / (sqrt(tmp + 1e-32));
  for (int w = 0; w < dim; w++) {
    s[w] = (ad_fv->v[w] - ad_fv->x_dot[w]) * tmp;
  }
  ADType mags = 0;
  for (int w = 0; w < dim; w++) {
    mags += (grad_U[w] * grad_U[w]);
  }
  mags = 1.0 / (sqrt(mags + 1e-32));
  for (int w = 0; w < dim; w++) {
    r[w] = grad_U[w] * mags;
  }

  ADType he = 0.0;
  for (int q = 0; q < ei[pg->imtrx]->dof[eqn]; q++) {
    ADType tmp = 0;
    for (int w = 0; w < dim; w++) {
      tmp += ad_fv->basis[eqn].grad_phi[q][w] * ad_fv->basis[eqn].grad_phi[q][w];
    }
    he += 1.0 / sqrt(std::max(tmp, 1e-20));
  }

  tmp = 0;
  for (int q = 0; q < ei[pg->imtrx]->dof[eqn]; q++) {
    for (int w = 0; w < dim; w++) {
      tmp += fabs(r[w] * ad_fv->basis[eqn].grad_phi[q][w]);
    }
  }
  ADType hrgn = 1.0 / (tmp + 1e-14);

  ADType magv = 0.0;
  for (int q = 0; q < VIM; q++) {
    magv += ad_fv->v[q] * ad_fv->v[q];
  }
  magv = sqrt(magv + 1e-32);

  ADType tau_dcdd = 0.5 * he * (1.0 / (mags + 1e-16)) * hrgn * hrgn;
  // ADType tau_dcdd = he * (1.0 / mags) * hrgn * hrgn / lambda;
  // printf("%g %g ", supg_tau.val(), tau_dcdd.val());
  // tau_dcdd = 1 / sqrt(1.0 / (supg_tau * supg_tau + 1e-32) +
  //                     1.0 / (tau_dcdd * tau_dcdd + 1e-32));
  // printf("%g \n ", tau_dcdd.val());
  ADType ss[DIM][DIM] = {{0.0}};
  ADType rr[DIM][DIM] = {{0.0}};
  ADType rdots = 0.0;
  for (int w = 0; w < dim; w++) {
    for (int z = 0; z < dim; z++) {
      ss[w][z] = s[w] * s[z];
      rr[w][z] = r[w] * r[z];
    }
    rdots += r[w] * s[w];
  }

  ADType inner_tensor[DIM][DIM] = {{0.0}};
  for (int w = 0; w < dim; w++) {
    for (int z = 0; z < dim; z++) {
      inner_tensor[w][z] = rr[w][z] - rdots * rdots * ss[w][z];
    }
  }

  for (int w = 0; w < dim; w++) {
    ADType tmp = 0.;
    for (int z = 0; z < dim; z++) {
      tmp += grad_U[w] * inner_tensor[w][z];
    }
    gs_inner_dot[w] = tmp;
    // gs_inner_dot[w] = grad_s[w][ii][jj];
  }
  return tau_dcdd;
}

template <typename scalar>
static int ad_calc_sa_S(scalar &S,                /* strain rate invariant */
                        scalar omega[DIM][DIM]) { /* strain rate tensor */
  S = 0.;
  /* get gamma_dot invariant for viscosity calculations */
  for (int a = 0; a < VIM; a++) {
    for (int b = 0; b < VIM; b++) {
      S += omega[a][b] * omega[a][b];
    }
  }

  S = sqrt(0.5 * (S) + 1e-16);

  return 0;
}

extern "C" void ad_sa_wall_func(double func[DIM],
                                double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE]) {
  // kind of hacky, near wall velocity is velocity at central node
  ADType unw[DIM];
  ADType eddy_nw;
  unw[0] = ADType(ad_fv->total_ad_variables, ad_fv->offset[VELOCITY1] + 8, *esp->v[0][8]);
  unw[1] = ADType(ad_fv->total_ad_variables, ad_fv->offset[VELOCITY2] + 8, *esp->v[1][8]);
  eddy_nw = ADType(ad_fv->total_ad_variables, ad_fv->offset[EDDY_NU] + 8, *esp->eddy_nu[8]);

  ADType mu = 0;
  dbl scale = 1.0;
  DENSITY_DEPENDENCE_STRUCT d_rho;
  if (gn->ConstitutiveEquation == TURBULENT_SA_DYNAMIC) {
    scale = density(&d_rho, tran->time_value);
  }
  int negative_mu_e = FALSE;
  if (fv_old->eddy_nu < 0) {
    negative_mu_e = TRUE;
  }

  double mu_newt = mp->viscosity;
  ADType fv1 = 1.0;
  if (negative_mu_e) {
    mu = 0;
  } else {

    ADType mu_e = eddy_nw;
    ADType cv1 = 7.1;
    ADType chi = mu_e / mu_newt;
    fv1 = pow(chi, 3) / (pow(chi, 3) + pow(cv1, 3));

    mu = scale * (mu_e * fv1);
    if (mu > 1e3 * mu_newt) {
      mu = 1e3 * mu_newt;
    }
  }

  ADType ut = 0;
  for (int i = 0; i < 2; i++) {
    ut += unw[i] * fv->stangent[0][i];
  }

  ADType normgv = 0;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      normgv += ad_fv->grad_v[i][j] * ad_fv->grad_v[i][j];
    }
  }
  normgv = std::sqrt(normgv);

  ADType mu_t = scale * ut * ut / (normgv + 1e-12);

  ADType mu_tt = std::max(mu_t - mu, 0);

  ADType omega[DIM][DIM];
  for (int a = 0; a < VIM; a++) {
    for (int b = 0; b < VIM; b++) {
      omega[a][b] = fv_old->grad_v[a][b] + fv_old->grad_v[b][a];
    }
  }

  ADType S = 0.;
  /* get gamma_dot invariant for viscosity calculations */
  for (int a = 0; a < VIM; a++) {
    for (int b = 0; b < VIM; b++) {
      S += omega[a][b] * omega[a][b];
    }
  }

  S = sqrt(0.5 * (S) + 1e-16);

  ADType eddy_nu = std::min(std::min(5, 1e-1 * scale * S * mp->viscosity), fv_old->eddy_nu * 1.5);

  func[0] = fv->eddy_nu - eddy_nu.val();

  // printf("eddynu = %g %g %g\n", eddy_nu.val(), fv->eddy_nu, func[0]);

  // for (int p = 0; p < WIM; p++) {
  //   for (int i = 0; i < ei[pg->imtrx]->dof[VELOCITY1 + p]; i++) {
  //     d_func[0][VELOCITY1 + p][i] = eddy_nu.dx(ad_fv->v_offset[p] + i);
  //   }
  // }

  // for (int i = 0; i < ei[pg->imtrx]->dof[EDDY_NU]; i++) {
  //   d_func[0][EDDY_NU][i] = eddy_nu.dx(ad_fv->eddy_nu_offset + i);
  // }
}

extern "C" void ad_omega_wall_func(double func[DIM],
                                   double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE]) {

  ADType d = std::max(fv->wall_distance, 0.0004);
  dbl nu = 1.5e-5;
  dbl beta1 = 0.075;
  ADType r = ad_fv->turb_omega - 6 * nu / (beta1 * d * d);
  func[0] = r.val();
  for (int j = 0; j < ei[pg->imtrx]->dof[TURB_OMEGA]; j++) {
    d_func[0][TURB_OMEGA][j] = r.dx(ad_fv->offset[TURB_OMEGA] + j);
  }
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

extern "C" int ad_assemble_spalart_allmaras(dbl time_value, /* current time */
                                            dbl tt, /* parameter to vary time integration from
                                                       explicit (tt = 1) to implicit (tt = 0)    */
                                            dbl dt, /* current time step size                    */
                                            const PG_DATA *pg_data) {

  //! WIM is the length of the velocity vector
  int i, j, a, b;
  int eqn, var, peqn, pvar;
  int *pdv = pd->v[pg->imtrx];

  int status = 0;

  eqn = EDDY_NU;
  double d_area = fv->wt * bf[eqn]->detJ * fv->h3;

  /* Get Eddy viscosity at Gauss point */
  ADType mu_e = ad_fv->eddy_nu;

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
  ADType omega[DIM][DIM];
  double omega_old[DIM][DIM];
  for (a = 0; a < VIM; a++) {
    for (b = 0; b < VIM; b++) {
      omega[a][b] = (ad_fv->grad_v[a][b] - ad_fv->grad_v[b][a]);
      omega_old[a][b] = (fv_old->grad_v[a][b] - fv_old->grad_v[b][a]);
    }
  }

  /* Vorticity */
  ADType S = 0.0;
  double S_old = 0;
  ad_calc_sa_S(S, omega);
  ad_calc_sa_S(S_old, omega_old);

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
  ADType chi = mu_e / mu_newt;
  ADType fv1 = pow(chi, 3) / (pow(chi, 3) + pow(cv1, 3));
  ADType fv2 = 1.0 - (chi) / (1.0 + chi * fv1);
  ADType fn = 1.0;
  if (negative_sa) {
    fn = (cn1 + pow(chi, 3.0)) / (cn1 - pow(chi, 3));
  }
  ADType Sbar = (mu_e * fv2) / (kappa * kappa * d * d);
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
  ADType S_e = S + Sbar;
  if (negative_Se) {
    S_e = S + S * (cv2 * cv2 * S + cv3 * Sbar) / ((cv3 - 2 * cv2) * S - Sbar);
  }
  double r_max = 10.0;
  ADType r = 0.0;
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
  ADType g = r + cw2 * (pow(r, 6) - r);
  ADType fw_inside = (1.0 + pow(cw3, 6)) / (pow(g, 6) + pow(cw3, 6));
  ADType fw = g * pow(fw_inside, (1.0 / 6.0));

  dbl supg = 1.;
  ADType supg_tau = 0;
  if (mp->Mwt_funcModel == GALERKIN) {
    supg = 0.;
  } else if (mp->Mwt_funcModel == SUPG || mp->Mwt_funcModel == SUPG_GP ||
             mp->Mwt_funcModel == SUPG_SHAKIB) {
    supg = mp->Mwt_func;
    ad_supg_tau_shakib(supg_tau, pd->Num_Dim, dt, mu_newt, EDDY_NU);
  }

  /*
   * Residuals_________________________________________________________________
   */

  std::vector<ADType> resid(ei[pg->imtrx]->dof[eqn]);
  for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
    resid[i] = 0;
  }

  if (af->Assemble_Residual) {
    /*
     * Assemble residual for eddy viscosity
     */
    eqn = EDDY_NU;
    peqn = upd->ep[pg->imtrx][eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      ADType wt_func = bf[eqn]->phi[i];

      if (supg > 0) {
        if (supg != 0.0) {
          for (int p = 0; p < VIM; p++) {
            wt_func += supg * supg_tau * ad_fv->v[p] * bf[eqn]->grad_phi[i][p];
          }
        }
      }

      /* Assemble mass term */
      ADType mass = 0.0;
      if (pd->TimeIntegration != STEADY) {
        if (pd->e[pg->imtrx][eqn] & T_MASS) {
          mass += ad_fv->eddy_nu_dot * wt_func * d_area;
          mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
        }
      }

      /* Assemble advection term */
      ADType adv = 0;
      for (int p = 0; p < VIM; p++) {
        adv += ad_fv->v[p] * ad_fv->grad_eddy_nu[p];
      }
      adv *= wt_func * d_area;
      adv *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];

      double neg_c = 1.0;
      if (negative_sa) {
        neg_c = -1.0;
      }

      /* Assemble source terms */
      ADType src_1 = cb1 * S_e * mu_e;
      ADType src_2 = neg_c * cw1 * fw * (mu_e * mu_e) / (d * d);
      ADType src = -src_1 + src_2;
      src *= wt_func * d_area;
      src *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

      /* Assemble diffusion terms */
      ADType diff_1 = 0.0;
      ADType diff_2 = 0.0;
      for (int p = 0; p < VIM; p++) {
        diff_1 += bf[eqn]->grad_phi[i][p] * (mu_newt + mu_e * fn) * ad_fv->grad_eddy_nu[p];
        diff_2 += wt_func * cb2 * ad_fv->grad_eddy_nu[p] * ad_fv->grad_eddy_nu[p];
      }
      ADType diff = (1.0 / sigma) * (diff_1 - diff_2);
      diff *= d_area;
      diff *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];

      resid[i] += mass + adv + src + diff;
      lec->R[LEC_R_INDEX(peqn, i)] += mass.val() + adv.val() + src.val() + diff.val();
    } /* end of for (i=0,ei[pg->imtrx]->dofs...) */
  } /* end of if assemble residual */

  /*
   * Jacobian terms...
   */

  if (af->Assemble_Jacobian) {
    eqn = EDDY_NU;
    peqn = upd->ep[pg->imtrx][eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {

      /* Sensitivity w.r.t. eddy viscosity */
      var = EDDY_NU;
      if (pdv[var]) {
        pvar = upd->vp[pg->imtrx][var];

        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += resid[i].dx(ad_fv->offset[EDDY_NU] + j);
        } /* End of loop over j */
      } /* End of if the variable is active */

      /* Sensitivity w.r.t. velocity */
      for (b = 0; b < VIM; b++) {
        var = VELOCITY1 + b;
        if (pdv[var]) {
          pvar = upd->vp[pg->imtrx][var];

          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += resid[i].dx(ad_fv->offset[var] + j);

          } /* End of loop over j */
        } /* End of if the variale is active */
      } /* End of loop over velocity components */

    } /* End of loop over i */
  } /* End of if assemble Jacobian */
  return (status);
}

ADType turb_omega_wall_bc(void) {
  double rho = density(NULL, tran->time_value);
  double beta1 = 0.075;
  double nu = mp->viscosity / rho;
  double Dy = 0.0;

  if (upd->turbulent_info->use_internal_wall_distance) {
    fv->wall_distance = 0.;
    if (pd->gv[pd->ShapeVar]) {
      int dofs = ei[upd->matrix_index[pd->ShapeVar]]->dof[pd->ShapeVar];
      for (int i = 0; i < dofs; i++) {
        dbl d =
            upd->turbulent_info
                ->wall_distances[ei[upd->matrix_index[pd->ShapeVar]]->gnn_list[pd->ShapeVar][i]];
        if (d > Dy) {
          Dy = d;
        }
      }
    }
  } else {
    GOMA_EH(GOMA_ERROR, "Unimplemented wall distance for turb_omega_wall_bc\n");
  }
  double omega_wall = 10.0 * (6 * nu) / (beta1 * Dy * Dy);
  return omega_wall;
}

static void
calc_blending_functions(dbl rho, DENSITY_DEPENDENCE_STRUCT *d_rho, ADType &F1, ADType &F2) {
  dbl d = fv->wall_distance;
  ADType k = ad_fv->turb_k;
  ADType omega = ad_fv->turb_omega;
  dbl nu = mp->viscosity / rho;
  dbl beta_star = 0.09;
  dbl sigma_omega2 = 0.856;
  ADType ddkdomega = 0;
  for (int i = 0; i < pd->Num_Dim; i++) {
    ddkdomega += ad_fv->grad_turb_k[i] * ad_fv->grad_turb_omega[i];
  }
  ADType CD_komega = std::max(2 * rho * sigma_omega2 * ddkdomega / (omega + 1e-20), 1e-20);
  ADType C500 = 500 * nu / (d * d * omega + 1e-20);

  ADType arg11 = 0;
  if (k > 0) {
    ADType arg11 = sqrt(k) / (beta_star * omega * d + 1e-20);
  }
  ADType arg12 = 4 * rho * sigma_omega2 * k / (CD_komega * d * d + 1e-20);
  ADType arg1 = std::min(std::max(arg11, C500), arg12);
  // ADType arg1 = std::min(std::max(sqrt(std::max(k+1e-20, 0)) / (beta_star * omega * d + 1e-20),
  // C500),
  //                        4 * rho * sigma_omega2 * k / (CD_komega * d * d + 1e-20));
  ADType arg2 = std::max(2 * arg11, C500);

  F1 = tanh(arg1 * arg1 * arg1 * arg1);
  F2 = tanh(arg2 * arg2);
#if 0
  if (d_F1 != NULL) {
    dbl d_dot_k_omega_dk[MDE];
    dbl d_dot_k_omega_domega[MDE];
    for (int i = 0; i < pd->Num_Dim; i++) {
      for (int j = 0; j < ei[pg->imtrx]->dof[TURB_K]; j++) {
        d_dot_k_omega_dk[j] = bf[TURB_K]->grad_phi[j][i] * fv->grad_turb_omega[i];
      }
      for (int j = 0; j < ei[pg->imtrx]->dof[TURB_OMEGA]; j++) {
        d_dot_k_omega_domega[j] = bf[TURB_OMEGA]->grad_phi[j][i] * fv->grad_turb_k[i];
      }
    }

    dbl d_CD_komega_dk[MDE] = {0.};
    dbl d_CD_komega_domega[MDE] = {0.};
    if (CD_komega > 1e-10) {
      for (int j = 0; j < ei[pg->imtrx]->dof[TURB_K]; j++) {
        d_CD_komega_dk[j] = 2 * rho * sigma_omega2 * d_dot_k_omega_dk[j] / omega;
      }
      for (int j = 0; j < ei[pg->imtrx]->dof[TURB_K]; j++) {
        d_CD_komega_domega[j] =
            2 * rho * sigma_omega2 * d_dot_k_omega_domega[j] / omega -
            2 * rho * sigma_omega2 * dot_k_omega * bf[TURB_OMEGA]->phi[j] / (omega * omega);
      }
    }

    
  }
#endif
}

/* assemble_turb_k -- assemble terms (Residual & Jacobian) for conservation
 *
 * SST-2003m turbulence model
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
 * Created:    July 2023 Weston Ortiz
 *
 */
extern "C" int ad_assemble_turb_k(dbl time_value, /* current time */
                                  dbl tt,         /* parameter to vary time integration from
                                                     explicit (tt = 1) to implicit (tt = 0)    */
                                  dbl dt,         /* current time step size                    */
                                  const PG_DATA *pg_data) {

  //! WIM is the length of the velocity vector
  int i;
  int eqn, peqn;
  int *pdv = pd->v[pg->imtrx];

  int status = 0;

  eqn = TURB_K;

  dbl d_area = fv->wt * bf[eqn]->detJ * fv->h3;

  dbl mu = mp->viscosity;
  DENSITY_DEPENDENCE_STRUCT d_rho_struct;
  DENSITY_DEPENDENCE_STRUCT *d_rho = &d_rho_struct;
  dbl rho = density(d_rho, time_value);
  ADType F1 = 0;
  ADType F2 = 0;
  calc_blending_functions(rho, d_rho, F1, F2);

  ADType SI;
  ADType gamma_dot[DIM][DIM];
  for (int i = 0; i < DIM; i++) {
    for (int j = 0; j < DIM; j++) {
      gamma_dot[i][j] = (ad_fv->grad_v[i][j] + ad_fv->grad_v[j][i]);
    }
  }
  ad_calc_shearrate(SI, gamma_dot);
  dbl a1 = 0.31;

  /* Rate of rotation tensor  */
  ADType omega[DIM][DIM];
  for (int a = 0; a < VIM; a++) {
    for (int b = 0; b < VIM; b++) {
      omega[a][b] = (ad_fv->grad_v[a][b] - ad_fv->grad_v[b][a]);
    }
  }

  /* Vorticity */
  ADType Omega = 0.0;
  calc_vort_mag(Omega, omega);

  dbl x_dot[DIM] = {0.};
  if (pd->gv[R_MESH1]) {
    for (int i = 0; i < DIM; i++) {
      x_dot[i] = fv_dot->x[i];
    }
  }

  dbl supg = 1.;
  ADType supg_tau = 0;
  if (mp->SAwt_funcModel == GALERKIN) {
    supg = 0.;
  } else if (mp->SAwt_funcModel == SUPG || mp->SAwt_funcModel == SUPG_GP ||
             mp->SAwt_funcModel == SUPG_SHAKIB) {
    supg = mp->SAwt_func;
    ad_supg_tau_shakib(supg_tau, pd->Num_Dim, dt, mu, TURB_K);
  }

  dbl beta_star = 0.09;
  dbl sigma_k1 = 0.85;
  dbl sigma_k2 = 1.0;

  // blended values
  ADType sigma_k = F1 * sigma_k1 + (1 - F1) * sigma_k2;

  ADType mu_t = rho * a1 * fv_old->turb_k / (std::max(a1 * ad_fv->turb_omega, Omega * F2) + 1e-16);
  ADType P = mu_t * SI * SI;
  ADType Plim = std::min(P, 20 * beta_star * rho * ad_fv->turb_omega * fv_old->turb_k);

  std::vector<ADType> resid(ei[pg->imtrx]->dof[eqn]);
  for (int i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
    resid[i] = 0;
  }
  /*
   * Residuals_________________________________________________________________
   */
  if (af->Assemble_Residual) {
    /*
     * Assemble residual for eddy viscosity
     */
    eqn = TURB_K;
    peqn = upd->ep[pg->imtrx][eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      ADType wt_func = bf[eqn]->phi[i];

      if (supg > 0) {
        if (supg != 0.0) {
          for (int p = 0; p < VIM; p++) {
            wt_func += supg * supg_tau * ad_fv->v[p] * bf[eqn]->grad_phi[i][p];
          }
        }
      }

      /* Assemble mass term */
      ADType mass = 0.0;
      if (pd->TimeIntegration != STEADY) {
        if (pd->e[pg->imtrx][eqn] & T_MASS) {
          mass += rho * ad_fv->turb_k_dot * wt_func * d_area;
          mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
        }
      }

      /* Assemble advection term */
      ADType adv = 0;
      for (int p = 0; p < VIM; p++) {
        adv += rho * (ad_fv->v[p] - x_dot[p]) * ad_fv->grad_turb_k[p];
      }
      adv *= wt_func * d_area;
      adv *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];

      /* Assemble source terms */
      ADType src = Plim - beta_star * rho * ad_fv->turb_omega * fv_old->turb_k;
      src *= -wt_func * d_area;
      src *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

      /* Assemble diffusion terms */
      ADType diff = 0.0;
      for (int p = 0; p < VIM; p++) {
        diff += bf[eqn]->grad_phi[i][p] * (mu + mu_t * sigma_k) * ad_fv->grad_turb_k[p];
      }
      diff *= d_area;
      diff *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];

      resid[i] += mass + adv + src + diff;
      lec->R[LEC_R_INDEX(peqn, i)] += mass.val() + adv.val() + src.val() + diff.val();
    } /* end of for (i=0,ei[pg->imtrx]->dofs...) */
  } /* end of if assemble residual */

  /*
   * Jacobian terms...
   */

  if (af->Assemble_Jacobian) {
    eqn = TURB_K;
    peqn = upd->ep[pg->imtrx][eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      for (int var = V_FIRST; var < V_LAST; var++) {

        /* Sensitivity w.r.t. velocity */
        if (pdv[var]) {
          int pvar = upd->vp[pg->imtrx][var];

          for (int j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            // J = &(lec->J[LEC_J_INDEX(peqn, pvar, ii, 0)]);
            lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += resid[i].dx(ad_fv->offset[var] + j);

          } /* End of loop over j */
        } /* End of if the variale is active */
      }
    }
  } /* End of if assemble Jacobian */
  return (status);
}

ADType ad_only_turb_k_omega_sst_viscosity(void) {
  ADType mu = 0;
  double mu_newt = mp->viscosity;
  dbl rho;
  rho = density(NULL, tran->time_value);
  ADType F1 = 0;
  ADType F2 = 0;
  calc_blending_functions(rho, NULL, F1, F2);

  ADType omega[DIM][DIM];
  for (int a = 0; a < VIM; a++) {
    for (int b = 0; b < VIM; b++) {
      omega[a][b] = (ad_fv->grad_v[a][b] - ad_fv->grad_v[b][a]);
    }
  }
  /* Vorticity */
  ADType Omega = 0.0;
  calc_vort_mag(Omega, omega);
  dbl a1 = 0.31;
  ADType mu_t = rho * a1 * ad_fv->turb_k / (std::max(a1 * ad_fv->turb_omega, Omega * F2) + 1e-16);

  mu = mu_newt + mu_t;
  return mu;
}

extern "C" dbl ad_turb_k_omega_sst_viscosity(VISCOSITY_DEPENDENCE_STRUCT *d_mu) {
  ADType mu = 0;
  double mu_newt = mp->viscosity;
  dbl rho;
  DENSITY_DEPENDENCE_STRUCT d_rho_struct;
  DENSITY_DEPENDENCE_STRUCT *d_rho = &d_rho_struct;
  rho = density(d_rho, tran->time_value);
  ADType F1 = 0;
  ADType F2 = 0;
  calc_blending_functions(rho, d_rho, F1, F2);

  ADType omega[DIM][DIM];
  for (int a = 0; a < VIM; a++) {
    for (int b = 0; b < VIM; b++) {
      omega[a][b] = (ad_fv->grad_v[a][b] - ad_fv->grad_v[b][a]);
    }
  }
  /* Vorticity */
  ADType Omega = 0.0;
  calc_vort_mag(Omega, omega);
  dbl a1 = 0.31;
  ADType mu_t = rho * a1 * ad_fv->turb_k / (std::max(a1 * ad_fv->turb_omega, Omega * F2) + 1e-16);

  mu = mu_newt + mu_t;
  if (d_mu != NULL) {
    for (int j = 0; j < ei[pg->imtrx]->dof[TURB_OMEGA]; j++) {
      d_mu->turb_omega[j] = mu.dx(ad_fv->offset[TURB_OMEGA] + j);
    }
    for (int j = 0; j < ei[pg->imtrx]->dof[TURB_K]; j++) {
      d_mu->turb_k[j] = mu.dx(ad_fv->offset[TURB_K] + j);
    }
    for (int b = 0; b < VIM; b++) {
      int var = VELOCITY1 + b;
      for (int j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_mu->v[b][j] = mu.dx(ad_fv->offset[var] + j);
      }
    }
  }
  return mu.val();
}

/* assemble_turb_omega -- assemble terms (Residual & Jacobian) for conservation
 *
 * k-omega SST
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
 * Created:    July 2023 Weston Ortiz
 *
 */
extern "C" int ad_assemble_turb_omega(dbl time_value, /* current time */
                                      dbl tt,         /* parameter to vary time integration from
                                                         explicit (tt = 1) to implicit (tt = 0)    */
                                      dbl dt, /* current time step size                    */
                                      const PG_DATA *pg_data) {

  //! WIM is the length of the velocity vector
  int i;
  int eqn, peqn;
  int *pdv = pd->v[pg->imtrx];

  int status = 0;

  eqn = TURB_OMEGA;

  dbl d_area = fv->wt * bf[eqn]->detJ * fv->h3;

  dbl mu = mp->viscosity;
  DENSITY_DEPENDENCE_STRUCT d_rho_struct;
  DENSITY_DEPENDENCE_STRUCT *d_rho = &d_rho_struct;
  dbl rho = density(d_rho, time_value);
  ADType F1 = 0;
  ADType F2 = 0;
  calc_blending_functions(rho, d_rho, F1, F2);

  ADType SI;
  ADType gamma_dot[DIM][DIM];
  for (int i = 0; i < DIM; i++) {
    for (int j = 0; j < DIM; j++) {
      gamma_dot[i][j] = (ad_fv->grad_v[i][j] + ad_fv->grad_v[j][i]);
    }
  }

  ADType omega[DIM][DIM];
  for (int a = 0; a < VIM; a++) {
    for (int b = 0; b < VIM; b++) {
      omega[a][b] = (ad_fv->grad_v[a][b] - ad_fv->grad_v[b][a]);
    }
  }
  /* Vorticity */
  ADType Omega = 0.0;
  calc_vort_mag(Omega, omega);

  ad_calc_shearrate(SI, gamma_dot);
  dbl a1 = 0.31;

  dbl x_dot[DIM] = {0.};
  if (pd->gv[R_MESH1]) {
    for (int i = 0; i < DIM; i++) {
      x_dot[i] = fv_dot->x[i];
    }
  }

  dbl supg = 1.;
  ADType supg_tau;
  if (mp->SAwt_funcModel == GALERKIN) {
    supg = 0.;
  } else if (mp->SAwt_funcModel == SUPG || mp->SAwt_funcModel == SUPG_GP ||
             mp->SAwt_funcModel == SUPG_SHAKIB) {
    supg = mp->SAwt_func;
    ad_supg_tau_shakib(supg_tau, pd->Num_Dim, dt, mu, TURB_OMEGA);
  }

  dbl beta_star = 0.09;
  dbl sigma_k1 = 0.85;
  dbl sigma_k2 = 1.0;
  dbl sigma_omega1 = 0.5;
  dbl sigma_omega2 = 0.856;
  dbl beta1 = 0.075;
  dbl beta2 = 0.0828;
  dbl gamma1 = beta1 / beta_star;
  dbl gamma2 = beta2 / beta_star;

  // blended values
  ADType gamma = F1 * gamma1 + (1 - F1) * gamma2;
  ADType sigma_k = F1 * sigma_k1 + (1 - F1) * sigma_k2;
  ADType sigma_omega = F1 * sigma_omega1 + (1 - F1) * sigma_omega2;
  ADType beta = F1 * beta1 + (1 - F1) * beta2;

  ADType mu_t = rho * a1 * ad_fv->turb_k / (std::max(a1 * fv_old->turb_omega, Omega * F2) + 1e-16);
  ADType P = mu_t * SI * SI;
  ADType Plim = std::min(P, 10 * beta_star * rho * fv_old->turb_omega * ad_fv->turb_k);

  std::vector<ADType> resid(ei[pg->imtrx]->dof[eqn]);
  for (int i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
    resid[i] = 0;
  }
  /*
   * Residuals_________________________________________________________________
   */
  if (af->Assemble_Residual) {
    /*
     * Assemble residual for eddy viscosity
     */
    eqn = TURB_OMEGA;
    peqn = upd->ep[pg->imtrx][eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      ADType wt_func = bf[eqn]->phi[i];

      if (supg > 0) {
        if (supg != 0.0) {
          for (int p = 0; p < VIM; p++) {
            wt_func += supg * supg_tau * ad_fv->v[p] * bf[eqn]->grad_phi[i][p];
          }
        }
      }

      /* Assemble mass term */
      ADType mass = 0.0;
      if (pd->TimeIntegration != STEADY) {
        if (pd->e[pg->imtrx][eqn] & T_MASS) {
          mass += rho * ad_fv->turb_omega_dot * wt_func * d_area;
          mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
        }
      }

      /* Assemble advection term */
      ADType adv = 0;
      for (int p = 0; p < VIM; p++) {
        adv += rho * (ad_fv->v[p] - x_dot[p]) * ad_fv->grad_turb_omega[p];
      }
      adv *= wt_func * d_area;
      adv *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];

      /* Assemble source terms */
      ADType src1 = (gamma * rho / (mu_t + 1e-16)) * Plim -
                    beta * rho * fv_old->turb_omega * fv_old->turb_omega;
      ADType src2 = 0;
      for (int p = 0; p < pd->Num_Dim; p++) {
        src2 += ad_fv->grad_turb_k[p] * fv_old->grad_turb_omega[p];
      }
      src2 *= 2 * (1 - F1) * rho * sigma_omega2 / (fv_old->turb_omega + 1e-16);
      ADType src = src1 + src2;
      src *= -wt_func * d_area;
      src *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

      /* Assemble diffusion terms */
      ADType diff = 0.0;
      for (int p = 0; p < VIM; p++) {
        diff += bf[eqn]->grad_phi[i][p] * (mu + mu_t * sigma_omega) * ad_fv->grad_turb_omega[p];
      }
      diff *= d_area;
      diff *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];

      resid[i] += mass + adv + src + diff;
      lec->R[LEC_R_INDEX(peqn, i)] += mass.val() + adv.val() + src.val() + diff.val();
    } /* end of for (i=0,ei[pg->imtrx]->dofs...) */
  } /* end of if assemble residual */

  /*
   * Jacobian terms...
   */

  if (af->Assemble_Jacobian) {
    eqn = TURB_OMEGA;
    peqn = upd->ep[pg->imtrx][eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      for (int var = V_FIRST; var < V_LAST; var++) {

        /* Sensitivity w.r.t. velocity */
        if (pdv[var]) {
          int pvar = upd->vp[pg->imtrx][var];

          for (int j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            // J = &(lec->J[LEC_J_INDEX(peqn, pvar, ii, 0)]);
            lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += resid[i].dx(ad_fv->offset[var] + j);

          } /* End of loop over j */
        } /* End of if the variale is active */
      }
    }
  } /* End of if assemble Jacobian */
  return (status);
}

ADType clipping_func(ADType x, ADType x_s, ADType x_e) {
  ADType psi = 0;
  ADType theta = M_PI_2 * (2 * x - (x_s + x_e)) / (x_s - x_e);

  if (x < x_s) {
    psi = 0;
  } else if (x > x_e) {
    psi = 1;
  } else {
    psi = 0.5 * (1 + std::sin(theta));
  }

  return psi;
}

ADType ad_only_turb_k_omega_viscosity(void) {
  ADType mu = 0;
  double mu_newt = mp->viscosity;
  dbl k_inf = upd->turbulent_info->k_inf;
  ADType psi = clipping_func(ad_fv->turb_k, 0, 10 * k_inf);
  ADType psi_neg = clipping_func(-ad_fv->turb_k, 0, 10 * k_inf);
  dbl rho = density(NULL, tran->time_value);

  ADType mu_t = rho * psi * ad_fv->turb_k / (std::exp(ad_fv->turb_omega));
  mu = mu_newt + mu_t;
  return mu;
}

extern "C" int ad_assemble_turb_k_modified(dbl time_value, /* current time */
                                           dbl tt, /* parameter to vary time integration from
                                                      explicit (tt = 1) to implicit (tt = 0)    */
                                           dbl dt, /* current time step size                    */
                                           const PG_DATA *pg_data) {

  //! WIM is the length of the velocity vector
  int i;
  int eqn, peqn;
  int *pdv = pd->v[pg->imtrx];

  int status = 0;

  eqn = TURB_K;

  dbl d_area = fv->wt * bf[eqn]->detJ * fv->h3;

  dbl omega_inf = std::exp(upd->turbulent_info->omega_inf);
  dbl k_inf = upd->turbulent_info->k_inf;

  // logarithmic formulation of k-omega model
  ADType omega = std::exp(ad_fv->turb_omega);
  ADType psi = clipping_func(ad_fv->turb_k, 0, 10 * k_inf);
  ADType psi_neg = clipping_func(-fv_old->turb_k, 0, 10 * k_inf);
  ADType k = psi * ad_fv->turb_k;
  // ADType k = psi * fv_old->turb_k;

  dbl sigma_k = 0.5;
  dbl beta_star = 9.0 / 100.0;

  ADType SI;
  ADType gamma_dot[DIM][DIM];
  for (int i = 0; i < DIM; i++) {
    for (int j = 0; j < DIM; j++) {
      gamma_dot[i][j] = (ad_fv->grad_v[i][j] + ad_fv->grad_v[j][i]);
    }
  }

  ADType Omega_tens[DIM][DIM];
  for (int a = 0; a < VIM; a++) {
    for (int b = 0; b < VIM; b++) {
      Omega_tens[a][b] = (ad_fv->grad_v[a][b] - ad_fv->grad_v[b][a]);
    }
  }
  /* Vorticity */
  ADType Omega = 0.0;
  calc_vort_mag(Omega, Omega_tens);

  ad_calc_shearrate(SI, gamma_dot);

  dbl rho = density(NULL, time_value);

  ADType mu_t = rho * k / (omega + 1e-16);
  // ADType mu_t_kdiff = rho * k / (omega + 1e-16) - rho * psi_neg * ad_fv->turb_k / omega_inf;
  ADType mu_t_kdiff = rho * k / (omega + 1e-16) - rho * psi_neg * fv_old->turb_k / omega_inf;
  dbl mu = mp->viscosity;

  ADType P = mu_t * Omega * Omega;

  P = std::min(P, 10 * beta_star * rho * omega * std::max(k, 0));

  dbl supg = 1.;
  ADType supg_tau;
  if (mp->SAwt_funcModel == GALERKIN) {
    supg = 0.;
  } else if (mp->SAwt_funcModel == SUPG || mp->SAwt_funcModel == SUPG_GP ||
             mp->SAwt_funcModel == SUPG_SHAKIB) {
    supg = mp->SAwt_func;
    ad_supg_tau_shakib(supg_tau, pd->Num_Dim, dt, mu_t, TURB_OMEGA);
  }

  std::vector<ADType> resid(ei[pg->imtrx]->dof[TURB_K]);
  for (int i = 0; i < ei[pg->imtrx]->dof[TURB_K]; i++) {
    resid[i] = 0;
  }
  /*
   * Residuals_________________________________________________________________
   */
  if (af->Assemble_Residual) {
    eqn = TURB_K;
    peqn = upd->ep[pg->imtrx][eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      ADType wt_func = bf[eqn]->phi[i];

      if (supg > 0) {
        if (supg != 0.0) {
          for (int p = 0; p < VIM; p++) {
            wt_func += supg * supg_tau * ad_fv->v[p] * bf[eqn]->grad_phi[i][p];
          }
        }
      }

      /* Assemble mass term */
      ADType mass = 0.0;
      if (pd->TimeIntegration != STEADY) {
        if (pd->e[pg->imtrx][eqn] & T_MASS) {
          mass += rho * ad_fv->turb_k_dot * wt_func * d_area;
          mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
        }
      }

      /* Assemble advection term */
      ADType adv = 0;
      for (int p = 0; p < VIM; p++) {
        adv += rho * (ad_fv->v[p] - ad_fv->x_dot[p]) * ad_fv->grad_turb_k[p];
      }
      adv *= wt_func * d_area;
      adv *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];

      /* Assemble source terms */
      ADType src1 = P - beta_star * rho * omega * std::max(k, 0);
      ADType src2 = 0;

      // ADType src2 = -beta_star * rho * omega_inf * psi_neg * ad_fv->turb_k;
      // ADType src2 = -beta_star * rho * omega_inf * psi_neg * ad_fv->turb_k;
      ADType src = (src1 + src2);
      // src = src1;
      src *= -wt_func * d_area;
      src *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

      /* Assemble diffusion terms */
      ADType diff = 0.0;
      for (int p = 0; p < VIM; p++) {
        diff += bf[eqn]->grad_phi[i][p] * (mu + mu_t_kdiff * sigma_k) * ad_fv->grad_turb_k[p];
      }
      diff *= d_area;
      diff *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];

      resid[i] -= mass + adv + src + diff;
      lec->R[LEC_R_INDEX(peqn, i)] -= mass.val() + adv.val() + src.val() + diff.val();
    } /* end of for (i=0,ei[pg->imtrx]->dofs...) */
  } /* end of if assemble residual */

  /*
   * Jacobian terms...
   */

  if (af->Assemble_Jacobian) {
    eqn = TURB_K;
    peqn = upd->ep[pg->imtrx][eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      for (int var = V_FIRST; var < V_LAST; var++) {

        /* Sensitivity w.r.t. velocity */
        if (pdv[var]) {
          int pvar = upd->vp[pg->imtrx][var];

          for (int j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            // J = &(lec->J[LEC_J_INDEX(peqn, pvar, ii, 0)]);
            lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += resid[i].dx(ad_fv->offset[var] + j);

          } /* End of loop over j */
        } /* End of if the variale is active */
      }
    }
  } /* End of if assemble Jacobian */
  return (status);
}
extern "C" int ad_assemble_turb_omega_modified(dbl time_value, /* current time */
                                               dbl tt, /* parameter to vary time integration from
                                                          explicit (tt = 1) to implicit (tt = 0) */
                                               dbl dt, /* current time step size */
                                               const PG_DATA *pg_data) {

  //! WIM is the length of the velocity vector
  int i;
  int eqn, peqn;
  int *pdv = pd->v[pg->imtrx];

  int status = 0;

  eqn = TURB_OMEGA;

  dbl d_area = fv->wt * bf[eqn]->detJ * fv->h3;

  dbl omega_inf = std::exp(upd->turbulent_info->omega_inf);
  dbl k_inf = upd->turbulent_info->k_inf;

  // logarithmic formulation of k-omega model
  ADType omega = std::exp(ad_fv->turb_omega);
  ADType psi = clipping_func(ad_fv->turb_k, 0, 10 * k_inf);
  ADType psi_neg = clipping_func(-ad_fv->turb_k, 0, 10 * k_inf);
  ADType k = psi * ad_fv->turb_k;

  dbl sigma_omega = 0.5;
  dbl beta_star = 9.0 / 100.0;
  dbl beta = 3.0 / 40.;
  dbl gamma = 5.0 / 9.0;

  ADType SI;
  ADType gamma_dot[DIM][DIM];
  for (int i = 0; i < DIM; i++) {
    for (int j = 0; j < DIM; j++) {
      gamma_dot[i][j] = (ad_fv->grad_v[i][j] + ad_fv->grad_v[j][i]);
    }
  }

  ADType Omega_tens[DIM][DIM];
  for (int a = 0; a < VIM; a++) {
    for (int b = 0; b < VIM; b++) {
      Omega_tens[a][b] = (ad_fv->grad_v[a][b] - ad_fv->grad_v[b][a]);
    }
  }
  /* Vorticity */
  ADType Omega = 0.0;
  calc_vort_mag(Omega, Omega_tens);

  ad_calc_shearrate(SI, gamma_dot);

  dbl rho = density(NULL, time_value);

  ADType mu_t = rho * k / (omega + 1e-16);
  ADType mu_t_kdiff = rho * k / (omega + 1e-16) - rho * psi_neg * ad_fv->turb_k / omega_inf;
  dbl mu = mp->viscosity;

  // ADType P = mu_t * SI * SI;
  ADType P = mu_t * Omega * Omega;
  P = std::min(P, 10 * beta_star * rho * omega * k);

  dbl supg = 1.;
  ADType supg_tau;
  if (mp->SAwt_funcModel == GALERKIN) {
    supg = 0.;
  } else if (mp->SAwt_funcModel == SUPG || mp->SAwt_funcModel == SUPG_GP ||
             mp->SAwt_funcModel == SUPG_SHAKIB) {
    supg = mp->SAwt_func;
    ad_supg_tau_shakib(supg_tau, pd->Num_Dim, dt, mu, TURB_OMEGA);
  }

  std::vector<ADType> resid(ei[pg->imtrx]->dof[TURB_OMEGA]);
  for (int i = 0; i < ei[pg->imtrx]->dof[TURB_OMEGA]; i++) {
    resid[i] = 0;
  }
  /*
   * Residuals_________________________________________________________________
   */
  if (af->Assemble_Residual) {
    /*
     * Assemble residual for eddy viscosity
     */
    eqn = TURB_OMEGA;
    peqn = upd->ep[pg->imtrx][eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      ADType wt_func = bf[eqn]->phi[i];

      if (supg > 0) {
        if (supg != 0.0) {
          for (int p = 0; p < VIM; p++) {
            wt_func += supg * supg_tau * ad_fv->v[p] * bf[eqn]->grad_phi[i][p];
          }
        }
      }

      /* Assemble mass term */
      ADType mass = 0.0;
      if (pd->TimeIntegration != STEADY) {
        if (pd->e[pg->imtrx][eqn] & T_MASS) {
          mass += rho * ad_fv->turb_omega_dot * wt_func * d_area;
          mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
        }
      }

      /* Assemble advection term */
      ADType adv = 0;
      for (int p = 0; p < VIM; p++) {
        adv += rho * (ad_fv->v[p] - ad_fv->x_dot[p]) * ad_fv->grad_turb_omega[p];
      }
      adv *= wt_func * d_area;
      adv *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];

      /* Assemble source terms */
      ADType src1 = (gamma / (k + 1e-16)) * P - beta * rho * omega;

      ADType src2 = 0;
      for (int p = 0; p < pd->Num_Dim; p++) {
        src2 += (mu + sigma_omega * mu_t) * ad_fv->grad_turb_omega[p] * ad_fv->grad_turb_omega[p];
      }
      ADType src = src1 + src2;
      src *= -wt_func * d_area;
      src *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

      /* Assemble diffusion terms */
      ADType diff = 0.0;
      for (int p = 0; p < VIM; p++) {
        diff += bf[eqn]->grad_phi[i][p] * (mu + mu_t * sigma_omega) * ad_fv->grad_turb_omega[p];
      }
      diff *= d_area;
      diff *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];

      resid[i] -= mass + adv + src + diff;
      lec->R[LEC_R_INDEX(peqn, i)] -= mass.val() + adv.val() + src.val() + diff.val();
    } /* end of for (i=0,ei[pg->imtrx]->dofs...) */

  } /* end of if assemble residual */

  /*
   * Jacobian terms...
   */

  if (af->Assemble_Jacobian) {
    eqn = TURB_OMEGA;
    peqn = upd->ep[pg->imtrx][eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      for (int var = V_FIRST; var < V_LAST; var++) {

        /* Sensitivity w.r.t. velocity */
        if (pdv[var]) {
          int pvar = upd->vp[pg->imtrx][var];

          for (int j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            // J = &(lec->J[LEC_J_INDEX(peqn, pvar, ii, 0)]);
            lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += resid[i].dx(ad_fv->offset[var] + j);

          } /* End of loop over j */
        } /* End of if the variale is active */
      }
    }
  } /* End of if assemble Jacobian */
  return (status);
}
extern "C" int
ad_assemble_turb_k_omega_modified(dbl time_value, /* current time */
                                  dbl tt,         /* parameter to vary time integration from
                                                     explicit (tt = 1) to implicit (tt = 0)    */
                                  dbl dt,         /* current time step size                    */
                                  const PG_DATA *pg_data) {

  //! WIM is the length of the velocity vector
  int i;
  int eqn, peqn;
  int *pdv = pd->v[pg->imtrx];

  int status = 0;

  eqn = TURB_OMEGA;

  dbl d_area = fv->wt * bf[eqn]->detJ * fv->h3;

  dbl omega_inf = upd->turbulent_info->omega_inf;
  dbl k_inf = upd->turbulent_info->k_inf;

  // logarithmic formulation of k-omega model
  ADType omega = std::exp(ad_fv->turb_omega);
  ADType psi = clipping_func(ad_fv->turb_k, 0, 10 * k_inf);
  ADType psi_neg = clipping_func(-ad_fv->turb_k, 0, 10 * k_inf);
  ADType k = psi * ad_fv->turb_k;

  dbl sigma_k = 0.5;
  dbl sigma_omega = 0.5;
  dbl beta_star = 9.0 / 100.0;
  dbl beta = 3.0 / 40.;
  dbl gamma = 5.0 / 9.0;

  ADType SI;
  ADType gamma_dot[DIM][DIM];
  for (int i = 0; i < DIM; i++) {
    for (int j = 0; j < DIM; j++) {
      gamma_dot[i][j] = (ad_fv->grad_v[i][j] + ad_fv->grad_v[j][i]);
    }
  }

  ADType Omega_tens[DIM][DIM];
  for (int a = 0; a < VIM; a++) {
    for (int b = 0; b < VIM; b++) {
      Omega_tens[a][b] = (ad_fv->grad_v[a][b] - ad_fv->grad_v[b][a]);
    }
  }
  /* Vorticity */
  ADType Omega = 0.0;
  calc_vort_mag(Omega, Omega_tens);

  ad_calc_shearrate(SI, gamma_dot);

  dbl rho = density(NULL, time_value);

  ADType mu_t = rho * k / (omega);
  // ADType mu_t_kdiff = rho * k / (omega);// - rho * psi_neg * ad_fv->turb_k / omega_inf;
  ADType mu_t_kdiff = rho * k / (omega)-rho * psi_neg * ad_fv->turb_k / omega_inf;
  dbl mu = mp->viscosity;

  ADType P = mu_t * Omega * Omega;

  dbl supg = 1.;
  // if (mp->SAwt_funcModel == GALERKIN) {
  // supg = 0.;
  // } else if (mp->SAwt_funcModel == SUPG || mp->SAwt_funcModel == SUPG_GP ||
  //  mp->SAwt_funcModel == SUPG_SHAKIB) {
  // supg = mp->SAwt_func;
  // }

  std::vector<std::vector<ADType>> resid(2);
  resid[0].resize(ei[pg->imtrx]->dof[TURB_OMEGA]);
  resid[1].resize(ei[pg->imtrx]->dof[TURB_K]);
  for (int i = 0; i < ei[pg->imtrx]->dof[TURB_OMEGA]; i++) {
    resid[0][i] = 0;
  }
  for (int i = 0; i < ei[pg->imtrx]->dof[TURB_K]; i++) {
    resid[1][i] = 0;
  }
  /*
   * Residuals_________________________________________________________________
   */
  if (af->Assemble_Residual) {
    ADType gs_inner_dot[DIM];
    ADType supg_tau_w, supg_tau_k;
    ad_supg_tau_shakib(supg_tau_w, pd->Num_Dim, dt, mu + sigma_omega * mu_t, TURB_OMEGA);
    ad_supg_tau_shakib(supg_tau_k, pd->Num_Dim, dt, mu + sigma_k * mu_t, TURB_K);
    ADType vshock_k = ad_dcdd(pd->Num_Dim, TURB_K, ad_fv->grad_turb_k, gs_inner_dot);
    vshock_k = std::min(supg_tau_k, vshock_k);
    ADType vshock_w = ad_dcdd(pd->Num_Dim, TURB_OMEGA, ad_fv->grad_turb_omega, gs_inner_dot);
    vshock_k = std::min(supg_tau_w, vshock_w);
    /*
     * Assemble residual for eddy viscosity
     */
    eqn = TURB_OMEGA;
    peqn = upd->ep[pg->imtrx][eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      ADType wt_func = bf[eqn]->phi[i];

      if (supg > 0) {
        if (supg != 0.0) {
          for (int p = 0; p < VIM; p++) {
            wt_func += supg * supg_tau_w * ad_fv->v[p] * bf[eqn]->grad_phi[i][p];
          }
        }
      }

      /* Assemble mass term */
      ADType mass = 0.0;
      if (pd->TimeIntegration != STEADY) {
        if (pd->e[pg->imtrx][eqn] & T_MASS) {
          mass += rho * ad_fv->turb_omega_dot * wt_func * d_area;
          mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
        }
      }

      /* Assemble advection term */
      ADType adv = 0;
      for (int p = 0; p < VIM; p++) {
        adv += rho * (ad_fv->v[p] - ad_fv->x_dot[p]) * ad_fv->grad_turb_omega[p];
      }
      adv *= wt_func * d_area;
      adv *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];

      /* Assemble source terms */
      ADType src1 = (gamma / (k + 1e-16)) * std::min(P, 10 * beta_star * rho * omega * k) -
                    beta * rho * omega;

      ADType src2 = 0;
      for (int p = 0; p < pd->Num_Dim; p++) {
        src2 += (mu + sigma_omega * mu_t) * ad_fv->grad_turb_omega[p] * ad_fv->grad_turb_omega[p];
      }
      ADType src = src1 + src2;
      src *= -wt_func * d_area;
      src *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

      /* Assemble diffusion terms */
      ADType diff = 0.0;
      for (int p = 0; p < VIM; p++) {
        diff += bf[eqn]->grad_phi[i][p] * (mu + mu_t * sigma_omega + vshock_w) *
                ad_fv->grad_turb_omega[p];
      }
      diff *= d_area;
      diff *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];

      resid[0][i] -= mass + adv + src + diff;
      lec->R[LEC_R_INDEX(peqn, i)] -= mass.val() + adv.val() + src.val() + diff.val();
    } /* end of for (i=0,ei[pg->imtrx]->dofs...) */

    eqn = TURB_K;
    peqn = upd->ep[pg->imtrx][eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      ADType wt_func = bf[eqn]->phi[i];

      if (supg > 0) {
        if (supg != 0.0) {
          for (int p = 0; p < VIM; p++) {
            wt_func += supg * supg_tau_k * ad_fv->v[p] * bf[eqn]->grad_phi[i][p];
          }
        }
      }

      /* Assemble mass term */
      ADType mass = 0.0;
      if (pd->TimeIntegration != STEADY) {
        if (pd->e[pg->imtrx][eqn] & T_MASS) {
          mass += rho * ad_fv->turb_k_dot * wt_func * d_area;
          mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
        }
      }

      /* Assemble advection term */
      ADType adv = 0;
      for (int p = 0; p < VIM; p++) {
        adv += rho * (ad_fv->v[p] - ad_fv->x_dot[p]) * ad_fv->grad_turb_k[p];
      }
      adv *= wt_func * d_area;
      adv *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];

      /* Assemble source terms */
      ADType src1 = std::min(P, 10 * beta_star * rho * omega * k) - beta_star * rho * omega * k;
      // ADType src2 = 0;

      ADType src2 = -beta_star * rho * omega_inf * psi_neg * ad_fv->turb_k;
      ADType src = src1 + src2;
      src *= -wt_func * d_area;
      src *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

      /* Assemble diffusion terms */
      ADType diff = 0.0;
      for (int p = 0; p < VIM; p++) {
        diff += bf[eqn]->grad_phi[i][p] * (mu + mu_t_kdiff * sigma_k + vshock_k) *
                ad_fv->grad_turb_k[p];
      }
      diff *= d_area;
      diff *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];

      resid[1][i] -= mass + adv + src + diff;
      lec->R[LEC_R_INDEX(peqn, i)] -= mass.val() + adv.val() + src.val() + diff.val();
    } /* end of for (i=0,ei[pg->imtrx]->dofs...) */
  } /* end of if assemble residual */

  /*
   * Jacobian terms...
   */

  if (af->Assemble_Jacobian) {
    eqn = TURB_OMEGA;
    peqn = upd->ep[pg->imtrx][eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      for (int var = V_FIRST; var < V_LAST; var++) {

        /* Sensitivity w.r.t. velocity */
        if (pdv[var]) {
          int pvar = upd->vp[pg->imtrx][var];

          for (int j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            // J = &(lec->J[LEC_J_INDEX(peqn, pvar, ii, 0)]);
            lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += resid[0][i].dx(ad_fv->offset[var] + j);

          } /* End of loop over j */
        } /* End of if the variale is active */
      }
    }
    eqn = TURB_K;
    peqn = upd->ep[pg->imtrx][eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      for (int var = V_FIRST; var < V_LAST; var++) {

        /* Sensitivity w.r.t. velocity */
        if (pdv[var]) {
          int pvar = upd->vp[pg->imtrx][var];

          for (int j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            // J = &(lec->J[LEC_J_INDEX(peqn, pvar, ii, 0)]);
            lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += resid[1][i].dx(ad_fv->offset[var] + j);

          } /* End of loop over j */
        } /* End of if the variale is active */
      }
    }
  } /* End of if assemble Jacobian */
  return (status);
}

void compute_sst_blending(ADType &F1, ADType &F2) {
  dbl mu = mp->viscosity;
  dbl rho = density(NULL, tran->time_value);

  ADType grad_k_dot_grad_omega = 0.;
  for (int i = 0; i < DIM; i++) {
    grad_k_dot_grad_omega += ad_fv->grad_turb_k[i] * ad_fv->grad_turb_omega[i];
  }
  ADType omega = std::max(1e-20, ad_fv->turb_omega);
  ADType k = std::max(1e-20, ad_fv->turb_k);
  dbl sigma_omega2 = 0.856;
  dbl beta_star = 0.09;
  dbl d = fv->wall_distance;

  ADType CD_kw = std::max(2 * rho * sigma_omega2 * (1 / omega) * grad_k_dot_grad_omega, 1e-10);

  ADType arg1 =
      std::min(std::max(std::sqrt(k) / (beta_star * omega * d), 500 * (mu / rho) / (d * d * omega)),
               4 * rho * sigma_omega2 * k / (CD_kw * d * d));

  ADType arg2 =
      std::max(2 * std::sqrt(k) / (beta_star * omega * d), 500 * (mu / rho) / (d * d * omega));

  F1 = std::tanh(arg1 * arg1 * arg1 * arg1);

  F2 = std::tanh(arg2 * arg2);
}

ADType ad_yzbeta(int dim, int eqn, ADType Y, ADType Z, const ADType grad_U[DIM]) {

  ADType invY = 1.0 / Y;
  ADType yzbeta = 0.0;

  ADType js = 0.0;
  for (int i = 0; i < dim; i++) {
    js += grad_U[i] * grad_U[i];
  }
  js = sqrt(std::max(js, 1e-20));

  ADType hdc = 0.0;
  for (int i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
    for (int a = 0; a < dim; a++) {
      hdc += std::abs(grad_U[a] / js) * ad_fv->basis[eqn].grad_phi[i][a];
    }
  }
  hdc = 2.0 / std::max(hdc, 1e-20);

  ADType inner = 0.0;
  for (int i = 0; i < dim; i++) {
    inner += invY * invY * grad_U[i] * grad_U[i];
  }
  inner = 1.0 / sqrt(std::max(inner, 1e-20));

  ADType yzbeta1 = std::abs(invY * Z) * inner;
  ADType yzbeta2 = std::abs(invY * Z) * std::pow(hdc / 2.0, 2.0);
  return 0.5 * (yzbeta1 + yzbeta2);
}

#if 1
extern "C" int ad_assemble_k_omega_sst_modified(dbl time_value, /* current time */
                                                dbl tt, /* parameter to vary time integration from
                                                           explicit (tt = 1) to implicit (tt = 0) */
                                                dbl dt, /* current time step size */
                                                const PG_DATA *pg_data) {

  //! WIM is the length of the velocity vector
  int i;
  int eqn, peqn;
  int *pdv = pd->v[pg->imtrx];

  int status = 0;

  eqn = TURB_K;

  dbl supg = 1.0;

  ADType G[DIM][DIM];

  ad_get_metric_tensor(ad_fv->B, pd->Num_Dim, ei[pg->imtrx]->ielem_type, G);

  ADType v_d_gv = 0;
  for (int i = 0; i < pd->Num_Dim; i++) {
    for (int j = 0; j < pd->Num_Dim; j++) {
      v_d_gv += fabs((ad_fv->v[i] - ad_fv->x_dot[i]) * G[i][j] * (ad_fv->v[j] - ad_fv->x_dot[j]));
    }
  }

  dbl rho = density(NULL, time_value);

  ADType W[DIM][DIM];
  for (int i = 0; i < DIM; i++) {
    for (int j = 0; j < DIM; j++) {
      W[i][j] = 0.5 * (ad_fv->grad_v[i][j] - ad_fv->grad_v[j][i]);
    }
  }

  ADType Omega = 0.0;
  for (int i = 0; i < DIM; i++) {
    for (int j = 0; j < DIM; j++) {
      Omega += W[i][j] * W[i][j];
    }
  }
  Omega = sqrt(std::max(Omega, 1e-20));

  dbl mu = mp->viscosity;

  ADType grad_k_dot_grad_omega = 0.;
  for (int i = 0; i < pd->Num_Dim; i++) {
    grad_k_dot_grad_omega += ad_fv->grad_turb_k[i] * ad_fv->grad_turb_omega[i];
  }

  ADType omega = std::max(1e-20, ad_fv->turb_omega);
  ADType k = std::max(1e-20, ad_fv->turb_k);
  dbl sigma_omega1 = 0.5;
  dbl sigma_omega2 = 0.856;
  dbl sigma_k1 = 0.85;
  dbl sigma_k2 = 1.0;
  dbl beta1 = 0.075;
  dbl beta2 = 0.0828;
  dbl beta_star = 0.09;
  // dbl kappa = 0.41;

  dbl gamma1 = 5.0 / 9.0; // beta1 / beta_star - sigma_omega1 * kappa * kappa / sqrt(beta_star);
  dbl gamma2 = 0.44;      // beta2 / beta_star - sigma_omega2 * kappa * kappa / sqrt(beta_star);

  ADType F1, F2;
  compute_sst_blending(F1, F2);

  ADType mu_turb = sst_viscosity(Omega, F2);

  ADType sigma_k = sigma_k1 * F1 + sigma_k2 * (1 - F1);
  ADType sigma_omega = sigma_omega1 * F1 + sigma_omega2 * (1 - F1);
  ADType beta = beta1 * F1 + beta2 * (1 - F1);

  ADType gamma = gamma1 * F1 + gamma2 * (1 - F1);

  ADType P = mu_turb * Omega * Omega;
  ADType Plim = std::min(P, 10 * beta_star * rho * omega * k);

  ADType d_area = fv->wt * ad_fv->detJ * fv->h3;
  std::vector<std::vector<ADType>> resid(2);
  resid[0].resize(ei[pg->imtrx]->dof[TURB_K]);
  resid[1].resize(ei[pg->imtrx]->dof[TURB_OMEGA]);
  for (int i = 0; i < ei[pg->imtrx]->dof[TURB_K]; i++) {
    resid[0][i] = 0;
  }
  for (int i = 0; i < ei[pg->imtrx]->dof[TURB_OMEGA]; i++) {
    resid[1][i] = 0;
  }
  /*
   * Residuals_________________________________________________________________
   */
  if (af->Assemble_Residual) {
    /*
     * Assemble residual for eddy viscosity
     */
    eqn = TURB_K;
    peqn = upd->ep[pg->imtrx][eqn];

    ADType Z_k = rho * ad_fv->turb_k_dot;
    for (int i = 0; i < pd->Num_Dim; i++) {
      Z_k += rho * ad_fv->v[i] * ad_fv->grad_turb_k[i];
    }
    ADType coeff = 12 * (mu + mu_turb * sigma_k) * (mu + mu_turb * sigma_k);

    ADType diff_g_g = 0;
    for (int i = 0; i < pd->Num_Dim; i++) {
      for (int j = 0; j < pd->Num_Dim; j++) {
        diff_g_g += coeff * G[i][j] * G[i][j];
      }
    }

    ADType sugn1 = 0;
    for (int i = 0; i < pd->Num_Dim; i++) {
      for (int j = 0; j < ei[pg->imtrx]->dof[eqn]; j++) {
        sugn1 += std::abs(ad_fv->v[i] * bf[eqn]->grad_phi[j][i]);
      }
    }
    sugn1 = 1.0 / std::max(sugn1, 1e-20);

    ADType sugn2 = dt / 2.0;

    ADType r[DIM] = {0};
    ADType norm_grad_k = 0;
    for (int i = 0; i < pd->Num_Dim; i++) {
      norm_grad_k = ad_fv->grad_turb_k[i] * ad_fv->grad_turb_k[i];
    }
    norm_grad_k = sqrt(norm_grad_k);
    for (int i = 0; i < pd->Num_Dim; i++) {
      r[i] = std::abs(ad_fv->grad_turb_k[i]) / norm_grad_k;
    }

    ADType h_rgn = 0;
    for (int i = 0; i < pd->Num_Dim; i++) {
      for (int j = 0; j < ei[pg->imtrx]->dof[eqn]; j++) {
        h_rgn += std::abs(r[i] * bf[eqn]->grad_phi[j][i]);
      }
    }
    h_rgn = 2.0 / std::max(h_rgn, 1e-20);

    ADType sugn3 = h_rgn * h_rgn * rho / (4 * (mu + mu_turb * sigma_k));
    sugn3 = 0;

    ADType tau_supg_k = 1.0 / (sqrt(4 / (dt * dt) + v_d_gv + diff_g_g));
    // ADType tau_supg_k = 1.0 / sqrt(1 / (sugn1 * sugn1) + 1 / (sugn2 * sugn2) + 1 / (sugn3 *
    // sugn3));

    // ADType vshock_k =
    // ad_yzbeta(pd->Num_Dim, TURB_K, upd->turbulent_info->k_inf, Z_k, ad_fv->grad_turb_k);

    ADType gs_inner_dot[DIM];
    ADType grad_turb_k[DIM] = {fv_old->grad_turb_k[0], fv_old->grad_turb_k[1],
                               fv_old->grad_turb_k[2]};
    ADType vshock_k = 0 * ad_dcdd(pd->Num_Dim, TURB_K, grad_turb_k, gs_inner_dot);

    vshock_k = std::min(tau_supg_k, vshock_k);
    // if (fv->wall_distance < 1) {
    //   // larger diffusion near wall
    //   vshock_k = std::max(vshock_k, 1e-2*std::exp(-fv->wall_distance * 20));
    // }

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      ADType wt_func = bf[eqn]->phi[i];

      if (supg > 0) {
        if (supg != 0.0) {
          for (int p = 0; p < VIM; p++) {
            wt_func += supg * tau_supg_k * ad_fv->v[p] * bf[eqn]->grad_phi[i][p];
          }
        }
      }

      /* Assemble mass term */
      ADType mass = 0.0;
      if (pd->TimeIntegration != STEADY) {
        if (pd->e[pg->imtrx][eqn] & T_MASS) {
          mass += rho * ad_fv->turb_k_dot;
          mass *= wt_func * d_area;
          mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
        }
      }

      /* Assemble advection term */
      ADType adv = 0;
      for (int p = 0; p < VIM; p++) {
        adv += rho * (ad_fv->v[p] - ad_fv->x_dot[p]) * ad_fv->grad_turb_k[p];
      }
      adv *= wt_func * d_area;
      adv *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];

      /* Assemble source terms */
      // Production - Destruction
      ADType src = Plim - beta_star * rho * omega * fv_old->turb_k;
      src *= -wt_func * d_area;
      src *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

      /* Assemble diffusion terms */
      ADType diff = 0.0;
      for (int p = 0; p < VIM; p++) {
        diff +=
            bf[eqn]->grad_phi[i][p] * (mu + mu_turb * sigma_k + vshock_k) * ad_fv->grad_turb_k[p];
        // + bf[eqn]->grad_phi[i][p] * vshock_k * gs_inner_dot[p];
      }
      diff *= d_area;
      diff *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];

      resid[0][i] -= mass + adv + src + diff;
      lec->R[LEC_R_INDEX(peqn, i)] -= mass.val() + adv.val() + src.val() + diff.val();
    } /* end of for (i=0,ei[pg->imtrx]->dofs...) */

    eqn = TURB_OMEGA;
    if (0 && fv->wall_distance < 0.005) {
      dbl d = std::max(0.0004, fv->wall_distance);
      dbl beta1 = 0.075;
      dbl nu = 1.57e-5;
      dbl ow = 10 * 6 * nu / (beta1 * d * d);
      peqn = upd->ep[pg->imtrx][eqn];
      for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
        ADType wt_func = bf[eqn]->phi[i];

        /* Assemble mass term */
        ADType mass = 0.0;
        if (pd->TimeIntegration != STEADY) {
          if (pd->e[pg->imtrx][eqn] & T_MASS) {
            mass += 1e1 * ad_fv->turb_omega * wt_func * d_area;
            mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
          }
        }

        /* Assemble advection term */
        ADType src = 0;
        for (int p = 0; p < VIM; p++) {
          src -= 1e1 * ow;
        }
        src *= wt_func * d_area;
        src *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

        /* Assemble diffusion terms */
        ADType diff = 0.0;
        // for (int p = 0; p < VIM; p++) {
        //   diff +=
        //       bf[eqn]->grad_phi[i][p] * (mu + mu_turb * sigma_omega) * ad_fv->grad_turb_omega[p];
        //   // + bf[eqn]->grad_phi[i][p] * vshock_w * gs_inner_dot[p];
        // }
        diff *= d_area;
        diff *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];

        resid[1][i] -= mass + src + diff;
        lec->R[LEC_R_INDEX(peqn, i)] -= mass.val() + src.val() + diff.val();
      } /* end of for (i=0,ei[pg->imtrx]->dofs...) */
    }
    {
      peqn = upd->ep[pg->imtrx][eqn];
      ADType Z_w = rho * ad_fv->turb_omega_dot;
      for (int i = 0; i < pd->Num_Dim; i++) {
        Z_w += rho * ad_fv->v[i] * ad_fv->grad_turb_omega[i];
      }

      coeff = 12 * (mu + mu_turb * sigma_omega) * (mu + mu_turb * sigma_omega);

      diff_g_g = 0;
      for (int i = 0; i < pd->Num_Dim; i++) {
        for (int j = 0; j < pd->Num_Dim; j++) {
          diff_g_g += coeff * G[i][j] * G[i][j];
        }
      }

      ADType norm_grad_omega = 0;
      for (int i = 0; i < pd->Num_Dim; i++) {
        norm_grad_omega = ad_fv->grad_turb_omega[i] * ad_fv->grad_turb_omega[i];
      }
      norm_grad_omega = sqrt(norm_grad_omega);
      for (int i = 0; i < pd->Num_Dim; i++) {
        r[i] = std::abs(ad_fv->grad_turb_omega[i]) / norm_grad_omega;
      }

      h_rgn = 0;
      for (int i = 0; i < pd->Num_Dim; i++) {
        for (int j = 0; j < ei[pg->imtrx]->dof[eqn]; j++) {
          h_rgn += std::abs(r[i] * bf[eqn]->grad_phi[j][i]);
        }
      }
      h_rgn = 2.0 / std::max(h_rgn, 1e-20);

      sugn3 = h_rgn * h_rgn / (4 * (mu + mu_turb * sigma_omega));
      // ADType tau_supg_w =
      // 1.0 / sqrt(1 / (sugn1 * sugn1) + 1 / (sugn2 * sugn2) + 1 / (sugn3 * sugn3));
      ADType tau_supg_w = 1.0 / (sqrt(4 / (dt * dt) + v_d_gv + diff_g_g));

      // ADType vshock_w = ad_yzbeta(pd->Num_Dim, TURB_OMEGA, upd->turbulent_info->omega_inf, Z_w,
      // ad_fv->grad_turb_omega);

      ADType grad_turb_w[DIM] = {fv_old->grad_turb_omega[0], fv_old->grad_turb_omega[1],
                                 fv_old->grad_turb_omega[2]};
      ADType vshock_w = 1 * ad_dcdd(pd->Num_Dim, TURB_OMEGA, grad_turb_w, gs_inner_dot);

      vshock_w = std::min(tau_supg_w, vshock_w);

      // if (fv->wall_distance < 1) {
      //   // larger diffusion near wall
      //   vshock_w = std::max(vshock_w, 1e-2*std::exp(-fv->wall_distance * 20));
      // }

      for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
        ADType wt_func = bf[eqn]->phi[i];

        if (supg > 0) {
          if (supg != 0.0) {
            for (int p = 0; p < VIM; p++) {
              wt_func += supg * tau_supg_w * ad_fv->v[p] * bf[eqn]->grad_phi[i][p];
            }
          }
        }

        /* Assemble mass term */
        ADType mass = 0.0;
        if (pd->TimeIntegration != STEADY) {
          if (pd->e[pg->imtrx][eqn] & T_MASS) {
            mass += rho * ad_fv->turb_omega_dot * wt_func * d_area;
            mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
          }
        }

        /* Assemble advection term */
        ADType adv = 0;
        for (int p = 0; p < VIM; p++) {
          adv += rho * (ad_fv->v[p] - ad_fv->x_dot[p]) * ad_fv->grad_turb_omega[p];
        }
        adv *= wt_func * d_area;
        adv *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];

        /* Assemble source terms */
        ADType src = (gamma * rho / mu_turb) * Plim - beta * rho * omega * omega +
                     2 * (1 - F1) * (rho * sigma_omega2 / omega) * grad_k_dot_grad_omega;
        // penalty for near wall
        // if (d < 0.4) {
        //   if (ad_fv->turb_omega > 2.422) {
        //   ADType red = 100*(ad_fv->turb_omega - 2.422);
        //   if (std::abs(red) > 10*std::abs(src)) {
        //     red = 10*std::abs(src) * red/std::abs(red);
        //   }
        //   src = -red;
        //   }
        // }
        src *= -wt_func * d_area;
        src *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

        /* Assemble diffusion terms */
        ADType diff = 0.0;
        for (int p = 0; p < VIM; p++) {
          diff += bf[eqn]->grad_phi[i][p] * (mu + mu_turb * sigma_omega + vshock_w) *
                  ad_fv->grad_turb_omega[p];
          // + bf[eqn]->grad_phi[i][p] * vshock_w * gs_inner_dot[p];
        }
        diff *= d_area;
        diff *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];

        resid[1][i] -= mass + adv + src + diff;
        lec->R[LEC_R_INDEX(peqn, i)] -= mass.val() + adv.val() + src.val() + diff.val();
      } /* end of for (i=0,ei[pg->imtrx]->dofs...) */
    }
  } /* end of if assemble residual */

  /*
   * Jacobian terms...
   */

  if (af->Assemble_Jacobian) {
    eqn = TURB_OMEGA;
    peqn = upd->ep[pg->imtrx][eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      for (int var = V_FIRST; var < V_LAST; var++) {

        /* Sensitivity w.r.t. velocity */
        if (pdv[var]) {
          int pvar = upd->vp[pg->imtrx][var];

          for (int j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            // J = &(lec->J[LEC_J_INDEX(peqn, pvar, ii, 0)]);
            lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += resid[1][i].dx(ad_fv->offset[var] + j);

          } /* End of loop over j */
        } /* End of if the variale is active */
      }
    }
    eqn = TURB_K;
    peqn = upd->ep[pg->imtrx][eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      for (int var = V_FIRST; var < V_LAST; var++) {

        /* Sensitivity w.r.t. velocity */
        if (pdv[var]) {
          int pvar = upd->vp[pg->imtrx][var];

          for (int j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            // J = &(lec->J[LEC_J_INDEX(peqn, pvar, ii, 0)]);
            lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += resid[0][i].dx(ad_fv->offset[var] + j);

          } /* End of loop over j */
        } /* End of if the variale is active */
      }
    }
  } /* End of if assemble Jacobian */
  return (status);
}
int ad_assemble_invariant(double tt, /* parameter to vary time integration from
                                      * explicit (tt = 1) to implicit (tt = 0)    */
                          double dt) /*  time step size                          */
{
  int dim;
  int p;

  int eqn;
  int peqn;
  int i;
  int status;

  dbl h3 = fv->h3; /* Volume element (scale factors). */
  /*
   * Galerkin weighting functions for i-th and a-th momentum residuals
   * and some of their derivatives...
   */

  ADType wt_func;

  /*
   * Interpolation functions for variables and some of their derivatives.
   */

  dbl wt;

  status = 0;

  /*
   * Unpack variables from structures for local convenience...
   */

  dim = pd->Num_Dim;

  /*
   * Bail out fast if there's nothing to do...
   */

  if (!pd->e[pg->imtrx][eqn = R_SHEAR_RATE]) {
    return (status);
  }

  peqn = upd->ep[pg->imtrx][eqn];

  wt = fv->wt; /* Numerical integration weight */

  ADType det_J = ad_fv->detJ; /* Really, ought to be mesh eqn. */

  ADType omega[DIM][DIM];
  for (int a = 0; a < VIM; a++) {
    for (int b = 0; b < VIM; b++) {
      omega[a][b] = ad_fv->grad_v[a][b] + ad_fv->grad_v[b][a];
    }
  }

  ADType S = 0.;
  /* get gamma_dot invariant for viscosity calculations */
  for (int a = 0; a < VIM; a++) {
    for (int b = 0; b < VIM; b++) {
      S += omega[a][b] * omega[a][b];
    }
  }

  if (S > 1e-20) {
    S = sqrt(0.5 * S);
  }

  /*
   * Residuals_________________________________________________________________
   */

  std::vector<ADType> resid(ei[pg->imtrx]->dof[eqn]);
  for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
    resid[i] = 0;
  }
  if (af->Assemble_Residual) {
    /*
     * Assemble the second_invariant equation
     */

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {

      wt_func = bf[eqn]->phi[i];

      ADType advection = 0.;

      if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
        advection = -S;
        advection *= wt_func * det_J * wt * h3;
        advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
      }

      /*
       * Diffusion Term..
       */

      /* OK this really isn't a diffusion term.  Its really a
       * filtering term.  But it looks like a diffusion operator.No?
       */

      ADType diffusion = 0.;

      if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {
        for (p = 0; p < dim; p++) {
          diffusion += ad_fv->basis[eqn].grad_phi[i][p] * ad_fv->grad_SH[p];
        }

        diffusion *= det_J * wt * h3;
        diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
      }

      /*
       * Source term...
       */

      ADType source = 0;

      if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
        source += ad_fv->SH;
        source *= wt_func * det_J * h3 * wt;
        source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
      }

      resid[i] += advection + source + diffusion;
      lec->R[LEC_R_INDEX(peqn, i)] += advection.val() + source.val() + diffusion.val();
    }
  }

  /*
   * Jacobian terms_________________________________________________________________
   */

  if (af->Assemble_Jacobian) {

    eqn = SHEAR_RATE;
    peqn = upd->ep[pg->imtrx][eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      for (int var = V_FIRST; var < V_LAST; var++) {

        /* Sensitivity w.r.t. velocity */
        if (pd->v[pg->imtrx][var]) {
          int pvar = upd->vp[pg->imtrx][var];

          for (int j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            // J = &(lec->J[LEC_J_INDEX(peqn, pvar, ii, 0)]);
            lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += resid[i].dx(ad_fv->offset[var] + j);

          } /* End of loop over j */
        } /* End of if the variale is active */
      }
    }
  } /* end of if(af,, */

  return (status);

} /* END of assemble_invariant */

extern "C" dbl visc_diss_heat_source_film_use_ad(HEAT_SOURCE_DEPENDENCE_STRUCT *d_h, dbl scale) {
  ADType h = 0;
  ADType gamma_dot[DIM][DIM];
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      gamma_dot[i][j] = ad_fv->grad_v[i][j] + ad_fv->grad_v[j][i];
    }
  }
  gamma_dot[2][2] = 2.0 * (-ad_fv->grad_v[0][0] - ad_fv->grad_v[1][1]);

  ADType gammadot;
  ad_calc_shearrate(gammadot, gamma_dot);

  ADType mu = ad_viscosity(gn, gamma_dot);

  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      h += mu * gamma_dot[j][i] * ad_fv->grad_v[i][j];
    }
  }
  h += mu * gamma_dot[2][2] * (-ad_fv->grad_v[0][0] - ad_fv->grad_v[1][1]);
  h = mu * gammadot * gammadot;
  h *= scale; // ad_fv->film_height;

  dbl alpha = mp->u_heat_source[0];
  dbl T_alpha = mp->u_heat_source[1];
  h += -alpha * (ad_fv->T - T_alpha) / ad_fv->film_height;

  for (int j = 0; j < ei[pg->imtrx]->dof[TEMPERATURE]; j++) {
    d_h->T[j] = h.dx(ad_fv->offset[TEMPERATURE] + j);
  }

  for (int j = 0; j < ei[pg->imtrx]->dof[FILM_HEIGHT]; j++) {
    d_h->film_height[j] = h.dx(ad_fv->offset[FILM_HEIGHT] + j);
  }

  for (int b = 0; b < pd->Num_Dim; b++) {
    for (int j = 0; j < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1 + b]; j++) {
      d_h->X[b][j] = h.dx(ad_fv->offset[MESH_DISPLACEMENT1 + b] + j);
    }
    for (int j = 0; j < ei[pg->imtrx]->dof[VELOCITY1 + b]; j++) {
      d_h->v[b][j] = h.dx(ad_fv->offset[VELOCITY1 + b] + j);
    }
  }

  return h.val();
}
#endif
#endif
