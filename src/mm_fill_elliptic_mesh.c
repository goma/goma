/************************************************************************ *
* Goma - Multiphysics finite element software                             *
* Sandia National Laboratories                                            *
*                                                                         *
* Copyright (c) 2023 Goma Developers, National Technology & Engineering   *
*               Solutions of Sandia, LLC (NTESS)                          *
*                                                                         *
* Under the terms of Contract DE-NA0003525, the U.S. Government retains   *
* certain rights in this software.                                        *
*                                                                         *
* This software is distributed under the GNU General Public License.      *
* See LICENSE file.                                                       *
\************************************************************************/

#include "mm_fill_elliptic_mesh.h"

#include "el_elm.h"
#include "mm_as.h"
#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "mm_eh.h"
#include "mm_fill_ptrs.h"
#include "mm_mp.h"
#include "rf_bc_const.h"
#include "rf_fem_const.h"
#include "std.h"
#include <math.h>

double elliptic_simple_abs_model(dbl xi, dbl xi_0, dbl a, dbl b, dbl c, dbl d, dbl e) {
  dbl x = fabs(xi - xi_0);
  return e + a / (b + c * x) + d * x;
}

double elliptic_dual_abs_model(dbl xi, dbl xi_0, dbl a, dbl b, dbl c, dbl d) {
  dbl x = fabs(xi) - fabs(xi_0);
  return d + a / (b + c * x);
}

int assemble_elliptic_mesh(void) {
  const int dim = pd->Num_Dim;

  int status = 0;
  /*
   * Unpack variables from structures for local convenience...
   */

  int eqn = R_MESH1; /* Well, yes, there really are 3, */
                     /* but for now just use the first */
                     /* to tell whether there is anything */
                     /* to do...*/
  dbl wt = fv->wt;

  if (!pd->e[pg->imtrx][eqn] || (ei[pg->imtrx]->ielem_dim < pd->Num_Dim)) {
    return -1;
  }

  // J[0][0] = dx/dxi
  // J[1][0] = dx/deta
  // J[2][0] = dx/dzeta
  // J[0][1] = dy/dxi
  // J[1][1] = dy/deta
  // J[2][1] = dy/dzeta
  // J[0][2] = dz/dxi
  // J[1][2] = dz/deta
  // J[2][2] = dz/dzeta

  dbl S[DIM] = {0};
  dbl d_S_dmesh[DIM][DIM][MDE] = {{{0}}};
  dbl T[DIM] = {0};
  dbl d_T_dmesh[DIM][DIM][MDE] = {{{0}}};
  if (dim == 2 || dim == 3) {
    dbl Snum = (SQUARE(bf[eqn]->J[0][0]) + SQUARE(bf[eqn]->J[0][1]));
    dbl Sden = (SQUARE(bf[eqn]->J[1][0]) + SQUARE(bf[eqn]->J[1][1]));
    S[0] = sqrt(Snum / Sden);
    S[1] = sqrt(Sden / Snum);
    T[0] = log(Snum);
    T[1] = log(Sden);
    for (int b = 0; b < dim; b++) {
      for (int m = 0; m < ei[pg->imtrx]->dof[R_MESH1 + b]; m++) {
        dbl dSnum = 2.0 * (bf[eqn]->dJ[0][0][b][m] * bf[eqn]->J[0][0] +
                           bf[eqn]->dJ[0][1][b][m] * bf[eqn]->J[0][1]);
        dbl dSden = 2.0 * (bf[eqn]->dJ[1][0][b][m] * bf[eqn]->J[1][0] +
                           bf[eqn]->dJ[1][1][b][m] * bf[eqn]->J[1][1]);
        d_S_dmesh[0][b][m] = S[0] * dSnum / (2 * Snum) - S[0] * dSden / (2 * Sden);
        d_S_dmesh[1][b][m] = S[1] * dSden / (2 * Sden) - S[1] * dSnum / (2 * Snum);
        d_T_dmesh[0][b][m] = dSnum / Snum;
        d_T_dmesh[1][b][m] = dSden / Sden;
      }
    }
  } else {
    GOMA_EH(GOMA_ERROR, "Unknown mesh dimension for elliptic mesh");
  }
  // } else if (dim == 3) {
  //   dbl Snum = (SQUARE(bf[eqn]->J[0][0]) + SQUARE(bf[eqn]->J[0][1]));
  //   dbl Sden = (SQUARE(bf[eqn]->J[1][0]) + SQUARE(bf[eqn]->J[1][1]));
  //   dbl Sxy = sqrt(Snum / Sden);
  //   dbl d_Sxy_dmesh[DIM][MDE] = {{0.}};
  //   T[0] = log((SQUARE(bf[eqn]->J[0][0]) + SQUARE(bf[eqn]->J[0][1]) + SQUARE(bf[eqn]->J[0][2])));
  //   for (int b = 0; b < dim; b++) {
  //     for (int m = 0; m < ei[pg->imtrx]->dof[R_MESH1 + b]; m++) {
  //       dbl dSnum = 2.0 * (bf[eqn]->dJ[0][0][b][m] * bf[eqn]->J[0][0] +
  //                          bf[eqn]->dJ[0][1][b][m] * bf[eqn]->J[0][1]);
  //       dbl dSden = 2.0 * (bf[eqn]->dJ[1][0][b][m] * bf[eqn]->J[1][0] +
  //                          bf[eqn]->dJ[1][1][b][m] * bf[eqn]->J[1][1]);
  //       d_Sxy_dmesh[b][m] = Sxy * dSnum / (2 * Snum) - Sxy * dSden / (2 * Sden);
  //       d_T_dmesh[0][b][m] = 2.0 *
  //                            (bf[eqn]->dJ[0][0][b][m] * bf[eqn]->J[0][0] +
  //                             bf[eqn]->dJ[0][1][b][m] * bf[eqn]->J[0][1] +
  //                             bf[eqn]->dJ[0][2][b][m] * bf[eqn]->J[0][2]) /
  //                            T[0];
  //     }
  //   }
  //
  //   Snum = (SQUARE(bf[eqn]->J[1][2]) + SQUARE(bf[eqn]->J[1][1]));
  //   Sden = (SQUARE(bf[eqn]->J[2][1]) + SQUARE(bf[eqn]->J[2][2]));
  //   dbl Syz = sqrt(Snum / Sden);
  //   T[1] = log((SQUARE(bf[eqn]->J[1][0]) + SQUARE(bf[eqn]->J[1][1])) + SQUARE(bf[eqn]->J[1][2]));
  //   dbl d_Syz_dmesh[DIM][MDE] = {{0.}};
  //   for (int b = 0; b < dim; b++) {
  //     for (int m = 0; m < ei[pg->imtrx]->dof[R_MESH1 + b]; m++) {
  //       dbl dSnum = 2.0 * (bf[eqn]->dJ[1][2][b][m] * bf[eqn]->J[1][2] +
  //                          bf[eqn]->dJ[1][1][b][m] * bf[eqn]->J[1][1]);
  //       dbl dSden = 2.0 * (bf[eqn]->dJ[2][1][b][m] * bf[eqn]->J[2][1] +
  //                          bf[eqn]->dJ[2][2][b][m] * bf[eqn]->J[2][2]);
  //       d_Syz_dmesh[b][m] = Syz * dSnum / (2 * Snum) - Syz * dSden / (2 * Sden);
  //       d_T_dmesh[1][b][m] = 2.0 *
  //                            (bf[eqn]->dJ[1][0][b][m] * bf[eqn]->J[1][0] +
  //                             bf[eqn]->dJ[1][1][b][m] * bf[eqn]->J[1][1] +
  //                             bf[eqn]->dJ[1][2][b][m] * bf[eqn]->J[1][2]) /
  //                            T[1];
  //     }
  //   }
  //
  //   Snum = (SQUARE(bf[eqn]->J[2][2]) + SQUARE(bf[eqn]->J[2][0]));
  //   Sden = (SQUARE(bf[eqn]->J[0][2]) + SQUARE(bf[eqn]->J[0][0]));
  //   dbl Szx = sqrt(Snum / Sden);
  //   T[2] = log((SQUARE(bf[eqn]->J[2][0]) + SQUARE(bf[eqn]->J[2][1]) + SQUARE(bf[eqn]->J[2][2])));
  //   dbl d_Szx_dmesh[DIM][MDE] = {{0.}};
  //   for (int b = 0; b < dim; b++) {
  //     for (int m = 0; m < ei[pg->imtrx]->dof[R_MESH1 + b]; m++) {
  //       dbl dSnum = 2.0 * (bf[eqn]->dJ[2][2][b][m] * bf[eqn]->J[2][2] +
  //                          bf[eqn]->dJ[2][0][b][m] * bf[eqn]->J[2][0]);
  //       dbl dSden = 2.0 * (bf[eqn]->dJ[0][2][b][m] * bf[eqn]->J[0][2] +
  //                          bf[eqn]->dJ[0][0][b][m] * bf[eqn]->J[0][0]);
  //       d_Szx_dmesh[b][m] = Szx * dSnum / (2 * Snum) - Szx * dSden / (2 * Sden);
  //       d_T_dmesh[2][b][m] = 2.0 *
  //                            (bf[eqn]->dJ[2][0][b][m] * bf[eqn]->J[2][0] +
  //                             bf[eqn]->dJ[2][1][b][m] * bf[eqn]->J[2][1] +
  //                             bf[eqn]->dJ[2][2][b][m] * bf[eqn]->J[2][2]) /
  //                            T[2];
  //     }
  //   }
  //
  //   // assemble components
  //   S[0] = Sxy + Szx;
  //   S[1] = 1 / Sxy + Syz;
  //   S[2] = 1 / Szx + 1 / Syz;
  //
  //   for (int b = 0; b < dim; b++) {
  //     for (int m = 0; m < ei[pg->imtrx]->dof[R_MESH1 + b]; m++) {
  //       d_S_dmesh[0][b][m] = d_Sxy_dmesh[b][m] + d_Szx_dmesh[b][m];
  //       d_S_dmesh[1][b][m] = d_Syz_dmesh[b][m];
  //       if (fabs(d_Sxy_dmesh[b][m]) > 0) {
  //         d_S_dmesh[1][b][m] += 1 / (d_Sxy_dmesh[b][m] + 1e-32);
  //       }
  //       d_S_dmesh[2][b][m] = 0;
  //       if (fabs(d_Szx_dmesh[b][m]) > 0) {
  //         d_S_dmesh[2][b][m] += 1 / (d_Szx_dmesh[b][m] + 1e-32);
  //       }
  //       if (fabs(d_Syz_dmesh[b][m]) > 0) {
  //         d_S_dmesh[2][b][m] += 1 / (d_Syz_dmesh[b][m]+1e-32);
  //       }
  //     }
  //   }
  //

  dbl eps_s = 1.2;

  dbl fxi = 1.0;
  if (elc_glob[ei[pg->imtrx]->mn]->fxi_model == CONSTANT) {
    fxi = elc_glob[ei[pg->imtrx]->mn]->fxi;
  } else if (elc_glob[ei[pg->imtrx]->mn]->fxi_model == ELLIPTIC_SIMPLE_ABS) {
    fxi = elliptic_simple_abs_model(
        fv->x0[0], elc_glob[ei[pg->imtrx]->mn]->u_fxi[0], elc_glob[ei[pg->imtrx]->mn]->u_fxi[1],
        elc_glob[ei[pg->imtrx]->mn]->u_fxi[2], elc_glob[ei[pg->imtrx]->mn]->u_fxi[3],
        elc_glob[ei[pg->imtrx]->mn]->u_fxi[4], elc_glob[ei[pg->imtrx]->mn]->u_fxi[5]);
  } else if (elc_glob[ei[pg->imtrx]->mn]->fxi_model == ELLIPTIC_DUAL_ABS) {
    fxi = elliptic_dual_abs_model(
        fv->x0[0], elc_glob[ei[pg->imtrx]->mn]->u_fxi[0], elc_glob[ei[pg->imtrx]->mn]->u_fxi[1],
        elc_glob[ei[pg->imtrx]->mn]->u_fxi[2], elc_glob[ei[pg->imtrx]->mn]->u_fxi[3],
        elc_glob[ei[pg->imtrx]->mn]->u_fxi[4]);
  } else {
    GOMA_EH(GOMA_ERROR, "Unknown Elliptic fxi model");
  }
  dbl geta = 1.0;
  if (elc_glob[ei[pg->imtrx]->mn]->geta_model == CONSTANT) {
    geta = elc_glob[ei[pg->imtrx]->mn]->geta;
  } else if (elc_glob[ei[pg->imtrx]->mn]->geta_model == ELLIPTIC_SIMPLE_ABS) {
    geta = elliptic_simple_abs_model(
        fv->x0[1], elc_glob[ei[pg->imtrx]->mn]->u_geta[0], elc_glob[ei[pg->imtrx]->mn]->u_geta[1],
        elc_glob[ei[pg->imtrx]->mn]->u_geta[2], elc_glob[ei[pg->imtrx]->mn]->u_geta[3],
        elc_glob[ei[pg->imtrx]->mn]->u_geta[4], elc_glob[ei[pg->imtrx]->mn]->u_geta[5]);
  } else if (elc_glob[ei[pg->imtrx]->mn]->geta_model == ELLIPTIC_DUAL_ABS) {
    geta = elliptic_dual_abs_model(
        fv->x0[1], elc_glob[ei[pg->imtrx]->mn]->u_geta[0], elc_glob[ei[pg->imtrx]->mn]->u_geta[1],
        elc_glob[ei[pg->imtrx]->mn]->u_geta[2], elc_glob[ei[pg->imtrx]->mn]->u_geta[3],
        elc_glob[ei[pg->imtrx]->mn]->u_geta[4]);
  } else {
    GOMA_EH(GOMA_ERROR, "Unknown Elliptic geta model");
  }
  dbl hzeta = 1.0;
  if (elc_glob[ei[pg->imtrx]->mn]->hzeta_model == CONSTANT) {
    hzeta = elc_glob[ei[pg->imtrx]->mn]->hzeta;
  } else if (elc_glob[ei[pg->imtrx]->mn]->hzeta_model == ELLIPTIC_SIMPLE_ABS) {
    hzeta = elliptic_simple_abs_model(
        fv->x0[2], elc_glob[ei[pg->imtrx]->mn]->u_hzeta[0], elc_glob[ei[pg->imtrx]->mn]->u_hzeta[1],
        elc_glob[ei[pg->imtrx]->mn]->u_hzeta[2], elc_glob[ei[pg->imtrx]->mn]->u_hzeta[3],
        elc_glob[ei[pg->imtrx]->mn]->u_hzeta[4], elc_glob[ei[pg->imtrx]->mn]->u_hzeta[5]);
  } else if (elc_glob[ei[pg->imtrx]->mn]->hzeta_model == ELLIPTIC_DUAL_ABS) {
    hzeta = elliptic_dual_abs_model(
        fv->x0[2], elc_glob[ei[pg->imtrx]->mn]->u_hzeta[0], elc_glob[ei[pg->imtrx]->mn]->u_hzeta[1],
        elc_glob[ei[pg->imtrx]->mn]->u_hzeta[2], elc_glob[ei[pg->imtrx]->mn]->u_hzeta[3],
        elc_glob[ei[pg->imtrx]->mn]->u_hzeta[4]);
  } else {
    GOMA_EH(GOMA_ERROR, "Unknown Elliptic hzeta model");
  }

  dbl sc[DIM] = {fxi, geta, hzeta};

  /*
   * Residuals_________________________________________________________________
   */

  if (af->Assemble_Residual) {
    /*
     * Assemble each component "a" of the momentum equation...
     */
    for (int aa = 0; aa < dim; aa++) {
      int a = (aa + 1) % dim;
      a = aa;
      int eqn = R_MESH1 + a;
      int peqn = upd->ep[pg->imtrx][eqn];

      int diffusion_on = pd->e[pg->imtrx][eqn] & T_DIFFUSION;
      int source_on = pd->e[pg->imtrx][eqn] & T_SOURCE;

      int diffusion_etm = pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
      int source_etm = pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

      for (int i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
        dbl diffusion = 0.;
        if (diffusion_on) {
          for (int p = 0; p < dim; p++) {
            diffusion += (S[a] + eps_s) * bf[eqn]->d_phi[i][p] * bf[eqn]->B[p][a];
          }
          diffusion *= bf[eqn]->detJ * wt;
          diffusion *= diffusion_etm;
        }

        dbl source = 0.;
        if (source_on) {
          source = -sc[a] * T[a] * bf[eqn]->dphidxi[i][a] * wt;

          source *= source_etm;
        }

        /*
         * porous term removed for mesh equation
         *  - the additional effects due  to porosity are entered
         *    into the consitutive equation for stress
         */

        lec->R[LEC_R_INDEX(peqn, i)] += diffusion + source;
      }
    }
  }

  /*
   * Jacobian terms...
   */

  if (af->Assemble_Jacobian) {
    for (int aa = 0; aa < dim; aa++) {
      int a = (aa + 1) % dim;
      a = aa;
      int eqn = R_MESH1 + a;
      int peqn = upd->ep[pg->imtrx][eqn];

      int diffusion_on = pd->e[pg->imtrx][eqn] & T_DIFFUSION;
      int source_on = pd->e[pg->imtrx][eqn] & T_SOURCE;

      int diffusion_etm = pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
      int source_etm = pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

      for (int i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
        /*
         * Set up some preliminaries that are needed for the (a,i)
         * equation for bunches of (b,j) column variables...
         */

        /*
         * J_d_d
         */
        for (int bb = 0; bb < dim; bb++) {
          int b = (bb + 1) % dim;
          b = bb;
          int var = MESH_DISPLACEMENT1 + b;
          if (pd->v[pg->imtrx][var]) {
            int pvar = upd->vp[pg->imtrx][var];
            for (int j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              dbl diffusion = 0.;
              if (diffusion_on) {
                dbl diff_a = 0.;
                dbl diff_b = 0.;

                for (int p = 0; p < dim; p++) {
                  diff_a +=
                      (d_S_dmesh[a][b][j]) * (bf[eqn]->d_phi[i][p] * bf[eqn]->B[p][a]) +
                      (S[a] + eps_s) * ((bf[eqn]->d_d_phi_dmesh[i][p][b][j] * bf[eqn]->B[p][a]) +
                                        (bf[eqn]->d_phi[i][p] * bf[eqn]->dB[p][a][b][j]));
                  diff_b += (S[a] + eps_s) * (bf[eqn]->d_phi[i][p] * bf[eqn]->B[p][a]);
                }

                diff_a *= bf[eqn]->detJ * wt;
                diff_b *= bf[eqn]->d_det_J_dm[b][j] * wt;

                diffusion = diff_a + diff_b;
                diffusion *= diffusion_etm;
              }

              dbl source = 0.;

              if (source_on) {
                source = -sc[a] * d_T_dmesh[a][b][j] * bf[eqn]->dphidxi[i][a] * wt;
                source *= source_etm;
              }

              lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += diffusion + source;
            }
          }
        }

      } /* end of loop over equations i  */
    } /* end of loop over equation directions a */
  } /* end of if jacobian */

  return (status);
}

void assemble_essential_elliptic_mesh(dbl func[DIM],
                                      dbl d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                                      int bc_name,
                                      dbl M) {
  int eqn = R_MESH1;

  dbl fxi = 1.0;
  dbl geta = 1.0;
  dbl hzeta = 1.0;

  dbl T[DIM] = {0};
  dbl d_T_dmesh[DIM][DIM][MDE] = {{{0}}};
  if (pd->Num_Dim == 2) {
    dbl Snum = (SQUARE(bf[eqn]->J[0][0]) + SQUARE(bf[eqn]->J[0][1]));
    dbl Sden = (SQUARE(bf[eqn]->J[1][0]) + SQUARE(bf[eqn]->J[1][1]));
    T[0] = log(Snum);
    T[1] = log(Sden);
    for (int b = 0; b < pd->Num_Dim; b++) {
      for (int m = 0; m < ei[pg->imtrx]->dof[R_MESH1 + b]; m++) {
        dbl dSnum = 2.0 * (bf[eqn]->dJ[0][0][b][m] * bf[eqn]->J[0][0] +
                           bf[eqn]->dJ[0][1][b][m] * bf[eqn]->J[0][1]);
        dbl dSden = 2.0 * (bf[eqn]->dJ[1][0][b][m] * bf[eqn]->J[1][0] +
                           bf[eqn]->dJ[1][1][b][m] * bf[eqn]->J[1][1]);
        d_T_dmesh[0][b][m] = dSnum / Snum;
        d_T_dmesh[1][b][m] = dSden / Sden;
      }
    }
  } else if (pd->Num_Dim == 3) {
    for (int a = 0; a < pd->Num_Dim; a++) {
      dbl inner = 0;
      for (int b = 0; b < pd->Num_Dim; b++) {
        inner += SQUARE(bf[eqn]->J[a][b]);
      }

      T[a] = log(inner);
      dbl inv = 1 / T[a];
      for (int b = 0; b < pd->Num_Dim; b++) {
        for (int r = 0; r < pd->Num_Dim; r++) {
          for (int m = 0; m < ei[pg->imtrx]->dof[R_MESH1 + b]; m++) {
            d_T_dmesh[0][b][m] += inv * 2.0 * bf[eqn]->dJ[a][b][r][m] * bf[eqn]->J[a][b];
          }
        }
      }
    }
  }

  switch (bc_name) {
  case ELLIPTIC_XI_REGULARIZATION_BC: {
    func[0] = -M * fxi * T[0] / fv->sdet;

    for (int b = 0; b < pd->Num_Dim; b++) {
      int var = MESH_DISPLACEMENT1 + b;
      for (int j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_func[0][var][j] = -M * fxi * d_T_dmesh[0][b][j] / fv->sdet +
                            M * fxi * T[0] * fv->dsurfdet_dx[b][j] / (fv->sdet * fv->sdet);
      }
    }
  } break;
  case ELLIPTIC_ETA_REGULARIZATION_BC: {

    func[0] = -M * geta * T[1] / fv->sdet;
    for (int b = 0; b < pd->Num_Dim; b++) {
      int var = MESH_DISPLACEMENT1 + b;
      for (int j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_func[0][var][j] = -M * geta * d_T_dmesh[1][b][j] / fv->sdet +
                            M * geta * T[1] * fv->dsurfdet_dx[b][j] / (fv->sdet * fv->sdet);
      }
    }
  } break;
  case ELLIPTIC_ZETA_REGULARIZATION_BC: {

    func[0] = -M * hzeta * T[2] / fv->sdet;
    for (int b = 0; b < pd->Num_Dim; b++) {
      int var = MESH_DISPLACEMENT1 + b;
      for (int j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_func[0][var][j] = -M * hzeta * d_T_dmesh[2][b][j] / fv->sdet +
                            M * hzeta * T[2] * fv->dsurfdet_dx[b][j] / (fv->sdet * fv->sdet);
      }
    }
  } break;
  default:
    GOMA_EH(GOMA_ERROR, "Unknown elliptic mesh regularization boundary condition");
    break;
  }
}
