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

#include "std.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "az_aztec.h"
#include "el_elm.h"
#include "el_geom.h"
#include "mm_as.h"
#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "mm_eh.h"
#include "mm_fill_ptrs.h"
#include "mm_fill_rs.h"
#include "mm_fill_shell.h"
#include "mm_fill_solid.h"
#include "mm_fill_species.h"
#include "mm_fill_terms.h"
#include "mm_fill_util.h"
#include "mm_mp.h"
#include "mm_mp_const.h"
#include "mm_mp_structs.h"
#include "mm_post_def.h"
#include "mm_shell_util.h"
#include "mm_std_models.h"
#include "mm_std_models_shell.h"
#include "mm_viscosity.h"
#include "rf_allo.h"
#include "rf_bc.h"
#include "rf_bc_const.h"
#include "rf_fem.h"
#include "rf_fem_const.h"
#include "rf_fill_const.h"
#include "rf_io.h"
#include "rf_io_const.h"
#include "rf_masks.h"
#include "rf_mp.h"
#include "rf_node_const.h"
#include "rf_solver.h"
#include "rf_solver_const.h"
#include "rf_vars_const.h"
#include "sl_util.h"
#include "std.h"
#include "user_mp.h"

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
  // J[0][1] = dy/dxi
  // J[1][1] = dy/deta

  dbl S[DIM] = {0};
  dbl d_S_dmesh[DIM][DIM][MDE] = {{{0}}};
  dbl T[DIM] = {0};
  dbl d_T_dmesh[DIM][DIM][MDE] = {{{0}}};
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

  dbl eps_s = 1.0;

  dbl fxi = 1.0;
  dbl geta = 1.0;

  dbl sc[DIM] = {fxi, geta};

  /*
   * Residuals_________________________________________________________________
   */

  if (af->Assemble_Residual) {
    /*
     * Assemble each component "a" of the momentum equation...
     */
    for (int aa = 0; aa < dim; aa++) {
      int a = (aa+1) % dim;
      int eqn = R_MESH1 + a;
      int peqn = upd->ep[pg->imtrx][eqn];

      int diffusion_on = pd->e[pg->imtrx][eqn] & T_DIFFUSION;
      int source_on = pd->e[pg->imtrx][eqn] & T_SOURCE;

      int diffusion_etm = pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
      int source_etm = pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

      for (int i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
        dbl diffusion = 0.;
        if (diffusion_on) {
          /*
           * use pseudo cartesian arbitrary mesh motion
           */
          for (int p = 0; p < dim; p++) {
            diffusion += (S[a] + eps_s) * bf[eqn]->d_phi[i][p] * bf[eqn]->B[p][a];
          }
          diffusion *= bf[eqn]->detJ * wt;

          diffusion *= diffusion_etm;
        }

        dbl source = 0.;
        if (source_on) {
          source = -sc[a] * T[a] * bf[eqn]->dphidxi[i][a] * wt;
          /* Source term only applies in Lagrangian mesh motion */

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
      int a = (aa+1) % dim;
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
          int b =  (bb+1) % dim;
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
    }   /* end of loop over equation directions a */
  }     /* end of if jacobian */

  return (status);
}

void assemble_essential_elliptic_mesh(dbl func[DIM],
                                      dbl d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                                      int bc_name,
                                      dbl M) {
  int eqn = R_MESH1;

  dbl fxi = 1.0;
  dbl geta = 1.0;

  dbl T[DIM] = {0};
  dbl d_T_dmesh[DIM][DIM][MDE] = {{{0}}};
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
                            M * fxi * T[1] * fv->dsurfdet_dx[b][j] / (fv->sdet * fv->sdet);
      }
    }
  } break;
  default:
    GOMA_EH(GOMA_ERROR, "Unknown elliptic mesh regularization boundary condition");
    break;
  }
}