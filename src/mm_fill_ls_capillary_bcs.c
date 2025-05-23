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
#include "mm_fill_ls_capillary_bcs.h"
#include "ac_particles.h"
#include "az_aztec.h"
#include "bc_colloc.h"
#include "bc_contact.h"
#include "density.h"
#include "el_elm.h"
#include "el_elm_info.h"
#include "el_geom.h"
#include "exo_struct.h"
#include "load_field_variables.h"
#include "mm_as.h"
#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "mm_dil_viscosity.h"
#include "mm_eh.h"
#include "mm_fill_aux.h"
#include "mm_fill_common.h"
#include "mm_fill_fill.h"
#include "mm_fill_ls.h"
#include "mm_fill_momentum.h"
#include "mm_fill_population.h"
#include "mm_fill_ptrs.h"
#include "mm_fill_rs.h"
#include "mm_fill_shell.h"
#include "mm_fill_solid.h"
#include "mm_fill_species.h"
#include "mm_fill_stabilization.h"
#include "mm_fill_stress.h"
#include "mm_fill_terms.h"
#include "mm_fill_util.h"
#include "mm_flux.h"
#include "mm_mp.h"
#include "mm_mp_const.h"
#include "mm_mp_structs.h"
#include "mm_ns_bc.h"
#include "mm_post_def.h"
#include "mm_qtensor_model.h"
#include "mm_shell_util.h"
#include "mm_species.h"
#include "mm_unknown_map.h"
#include "mm_viscosity.h"
#include "rf_allo.h"
#include "rf_bc.h"
#include "rf_bc_const.h"
#include "rf_fem.h"
#include "rf_fem_const.h"
#include "rf_solver.h"
#include "rf_vars_const.h"
#include "sl_util.h"
#include "sl_util_structs.h"
#include "std.h"
#include "stdbool.h"
#include "user_mp.h"
#include "user_mp_gen.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int continuous_surface_tension_old(double st,
                                          double csf[DIM][DIM],
                                          struct Level_Set_Interface *lsi_old) {
  int a, b;
  int status = 0;

  int var;

  var = ls->var;

  if (var != FILL) {
    GOMA_EH(GOMA_ERROR, "Unknown fill variable");
  }

  /* Fetch the level set interface functions. */
  load_lsi_old(ls->Length_Scale, lsi_old);

  /* If we're not near the zero level set, then csf is a zero tensor. */
  if (!lsi_old->near) {
    for (a = 0; a < VIM; a++) {
      for (b = 0; b < VIM; b++) {
        csf[a][b] = 0;
      }
    }

    return (status);
  }

  /* Calculate the CSF tensor. */
  for (a = 0; a < VIM; a++) {
    for (b = 0; b < VIM; b++) {
      csf[a][b] =
          st * lsi_old->delta * (((double)delta(a, b)) - lsi_old->normal[a] * lsi_old->normal[b]);
    }
  }

  return status;
}
int assemble_csf_tensor(void) {
  int i, j, a, b, p, q, ii, ledof;
  int eqn, peqn, var, pvar;

  struct Basis_Functions *bfm;
  dbl(*grad_phi_i_e_a)[DIM] = NULL;

  double wt, det_J, h3, phi_j, d_area;

  double csf[DIM][DIM];
  double d_csf_dF[DIM][DIM][MDE];
  double d_csf_dX[DIM][DIM][DIM][MDE];

  double source, source_etm;

  eqn = R_MOMENTUM1;
  if (!pd->e[pg->imtrx][eqn]) {
    return (0);
  }

  wt = fv->wt;
  h3 = fv->h3;
  if (ls->on_sharp_surf) /* sharp interface */
  {
    det_J = fv->sdet;
  } else /* diffuse interface */
  {
    det_J = bf[eqn]->detJ;
  }

  d_area = wt * det_J * h3;

  memset(csf, 0, sizeof(double) * DIM * DIM);
  memset(d_csf_dF, 0, sizeof(double) * DIM * DIM * MDE);
  memset(d_csf_dX, 0, sizeof(double) * DIM * DIM * DIM * MDE);

  continuous_surface_tension(mp->surface_tension, csf, d_csf_dF, d_csf_dX);

  /* finite difference calculation of path dependencies for
          subelement integration
          */
  if (ls->CalcSurfDependencies) {
    for (a = 0; a < WIM; a++) {
      eqn = R_MOMENTUM1 + a;
      peqn = upd->ep[pg->imtrx][eqn];
      bfm = bf[eqn];

      for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {

        ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][i];

        if (ei[pg->imtrx]->active_interp_ledof[ledof]) {

          ii = ei[pg->imtrx]->lvdof_to_row_lvdof[eqn][i];

          grad_phi_i_e_a = bfm->grad_phi_e[i][a];

          source = 0.;

          for (p = 0; p < VIM; p++) {
            for (q = 0; q < VIM; q++) {
              source += grad_phi_i_e_a[p][q] * csf[q][p];
            }
          }

          source *= -det_J * wt * h3;
          source *= pd->etm[pg->imtrx][eqn][LOG2_SOURCE];

          /* J_m_F
           */
          var = LS;
          pvar = upd->vp[pg->imtrx][var];

          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];
            lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source * phi_j;
          }
        }
      }
    }
    return (0);
  }

  /*
   * Residuals ____________________________________________________________________________
   */

  if (af->Assemble_Residual) {

    for (a = 0; a < WIM; a++) {
      eqn = R_MOMENTUM1 + a;
      peqn = upd->ep[pg->imtrx][eqn];
      bfm = bf[eqn];

      for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {

        ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][i];

        if (ei[pg->imtrx]->active_interp_ledof[ledof]) {

          ii = ei[pg->imtrx]->lvdof_to_row_lvdof[eqn][i];

          grad_phi_i_e_a = bfm->grad_phi_e[i][a];

          source = 0.;

#ifdef DO_NO_UNROLL

          for (p = 0; p < VIM; p++) {
            for (q = 0; q < VIM; q++) {
              source += grad_phi_i_e_a[p][q] * csf[q][p];
            }
          }
#else
          source += grad_phi_i_e_a[0][0] * csf[0][0];
          source += grad_phi_i_e_a[1][1] * csf[1][1];
          source += grad_phi_i_e_a[0][1] * csf[1][0];
          source += grad_phi_i_e_a[1][0] * csf[0][1];
          if (VIM == 3) {
            source += grad_phi_i_e_a[2][2] * csf[2][2];
            source += grad_phi_i_e_a[2][1] * csf[1][2];
            source += grad_phi_i_e_a[2][0] * csf[0][2];
            source += grad_phi_i_e_a[1][2] * csf[2][1];
            source += grad_phi_i_e_a[0][2] * csf[2][0];
          }
#endif

          source *= -det_J * wt * h3;
          source *= pd->etm[pg->imtrx][eqn][LOG2_SOURCE];

          lec->R[LEC_R_INDEX(peqn, ii)] += source;
        }
      }
    }
  }

  /*
   * Yacobian terms...
   */

  if (af->Assemble_Jacobian) {
    for (a = 0; a < WIM; a++) {
      eqn = R_MOMENTUM1 + a;
      peqn = upd->ep[pg->imtrx][eqn];
      bfm = bf[eqn];
      source_etm = pd->etm[pg->imtrx][eqn][LOG2_SOURCE];

      for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {

        ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][i];

        if (ei[pg->imtrx]->active_interp_ledof[ledof]) {

          ii = ei[pg->imtrx]->lvdof_to_row_lvdof[eqn][i];

          grad_phi_i_e_a = bfm->grad_phi_e[i][a];

          /* J_m_F
           */
          var = LS;
          if (pd->v[pg->imtrx][var]) {
            pvar = upd->vp[pg->imtrx][var];

            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              phi_j = bf[var]->phi[j];

              source = 0.;

#ifdef DO_NO_UNROLL

              for (p = 0; p < VIM; p++) {
                for (q = 0; q < VIM; q++) {
                  source += grad_phi_i_e_a[p][q] * d_csf_dF[q][p][j];
                }
              }

#else
              source += grad_phi_i_e_a[0][0] * d_csf_dF[0][0][j];
              source += grad_phi_i_e_a[0][1] * d_csf_dF[1][0][j];
              source += grad_phi_i_e_a[1][0] * d_csf_dF[0][1][j];
              source += grad_phi_i_e_a[1][1] * d_csf_dF[1][1][j];

              if (VIM == 3) {
                source += grad_phi_i_e_a[2][0] * d_csf_dF[0][2][j];
                source += grad_phi_i_e_a[2][1] * d_csf_dF[1][2][j];
                source += grad_phi_i_e_a[2][2] * d_csf_dF[2][2][j];
                source += grad_phi_i_e_a[0][2] * d_csf_dF[2][0][j];
                source += grad_phi_i_e_a[1][2] * d_csf_dF[2][1][j];
              }
#endif

              source *= -d_area;
              source *= source_etm;

              lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
            }
          }
          var = MESH_DISPLACEMENT1;
          if (pd->v[pg->imtrx][var]) {
            for (b = 0; b < VIM; b++) {
              var = MESH_DISPLACEMENT1 + b;
              pvar = upd->vp[pg->imtrx][var];

              for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                phi_j = bf[var]->phi[j];

                source = 0.;
#ifdef DO_NO_UNROLL
                for (p = 0; p < VIM; p++) {
                  for (q = 0; q < VIM; q++) {
                    source += grad_phi_i_e_a[p][q] * d_csf_dX[q][p][b][j];
                  }
                }
#else

                source += grad_phi_i_e_a[0][0] * d_csf_dX[0][0][b][j];
                source += grad_phi_i_e_a[0][1] * d_csf_dX[1][0][b][j];
                source += grad_phi_i_e_a[1][0] * d_csf_dX[0][1][b][j];
                source += grad_phi_i_e_a[1][1] * d_csf_dX[1][1][b][j];

                if (VIM == 3) {
                  source += grad_phi_i_e_a[2][0] * d_csf_dX[0][2][b][j];
                  source += grad_phi_i_e_a[2][1] * d_csf_dX[1][2][b][j];
                  source += grad_phi_i_e_a[2][2] * d_csf_dX[2][2][b][j];
                  source += grad_phi_i_e_a[0][2] * d_csf_dX[2][0][b][j];
                  source += grad_phi_i_e_a[1][2] * d_csf_dX[2][1][b][j];
                }
#endif

                source *= -d_area;
                source *= source_etm;

                lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
              }
            }
          }
        }
      }
    }
  }

  return (1);
}
int assemble_div_n_source(void) {
  int i, j, a, b, p, ii, ledof;
  int eqn, peqn, var, pvar;
  int dim;

  struct Basis_Functions *bfm;

  double wt, det_J, h3, phi_i, phi_j;
  double source, source1;

  eqn = R_MOMENTUM1;
  if (!pd->e[pg->imtrx][eqn]) {
    return (-1);
  }

  wt = fv->wt;
  h3 = fv->h3;

  if (ls->on_sharp_surf) /* sharp interface */
  {
    det_J = fv->sdet;
  } else /* diffuse interface */
  {
    det_J = bf[eqn]->detJ;
  }

  dim = pd->Num_Dim;

  /* finite difference calculation of path dependencies for
     subelement integration
  */
  if (ls->CalcSurfDependencies) {
    for (a = 0; a < WIM; a++) {
      eqn = R_MOMENTUM1 + a;
      peqn = upd->ep[pg->imtrx][eqn];
      bfm = bf[eqn];

      for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {

        ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][i];

        if (ei[pg->imtrx]->active_interp_ledof[ledof]) {

          ii = ei[pg->imtrx]->lvdof_to_row_lvdof[eqn][i];

          phi_i = bfm->phi[i];

          source = -mp->surface_tension * fv->div_n * fv->n[a] * lsi->delta;

          source *= phi_i * wt * det_J * h3;
          source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

          /* J_m_F
           */
          var = LS;
          pvar = upd->vp[pg->imtrx][var];

          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];
            lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source * phi_j;
          }
        }
      }
    }
    return (0);
  }

  /*
   * Residuals
   * ________________________________________________________________________________
   */

  if (af->Assemble_Residual) {
    for (a = 0; a < WIM; a++) {
      eqn = R_MOMENTUM1 + a;
      peqn = upd->ep[pg->imtrx][eqn];
      bfm = bf[eqn];

      for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {

        ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][i];

        if (ei[pg->imtrx]->active_interp_ledof[ledof]) {

          ii = ei[pg->imtrx]->lvdof_to_row_lvdof[eqn][i];

          phi_i = bfm->phi[i];

          source = 0.0;

          source = -mp->surface_tension * fv->div_n * fv->n[a] * lsi->delta;

          source *= phi_i * wt * det_J * h3;
          source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

          lec->R[LEC_R_INDEX(peqn, ii)] += source;
        }
      }
    }
  }

  /*
   * Jacobian terms...
   */

  if (af->Assemble_Jacobian) {
    for (a = 0; a < WIM; a++) {
      eqn = R_MOMENTUM1 + a;
      peqn = upd->ep[pg->imtrx][eqn];
      bfm = bf[eqn];

      for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {

        ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][i];

        if (ei[pg->imtrx]->active_interp_ledof[ledof]) {

          ii = ei[pg->imtrx]->lvdof_to_row_lvdof[eqn][i];

          phi_i = bfm->phi[i];

          /*
           * J_m_n
           */

          for (b = 0; b < dim; b++) {
            var = NORMAL1 + b;
            if (pd->v[pg->imtrx][var]) {

              pvar = upd->vp[pg->imtrx][var];

              for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {

                source = fv->div_n * bf[var]->phi[j] * delta(a, b);

                for (p = 0; p < VIM; p++) {
                  source += fv->n[a] * bf[var]->grad_phi_e[j][b][p][p];
                }
                source *= -mp->surface_tension * lsi->delta;
                source *= phi_i * wt * det_J * h3;
                source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

                lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
              }
            }
          }

          /*
           * J_m_F
           */

          var = LS;
          if (pd->v[pg->imtrx][var]) {
            pvar = upd->vp[pg->imtrx][var];

            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              source = 0.0;

              source = -mp->surface_tension * fv->div_n * fv->n[a];
              source *= lsi->d_delta_dF[j];
              source *= phi_i * wt * det_J * h3;
              source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

              lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
            }
          }

          /*
           * J_m_d
           */

          var = MESH_DISPLACEMENT1;

          if (pd->v[pg->imtrx][var]) {
            for (b = 0; b < dim; b++) {
              var = MESH_DISPLACEMENT1 + b;
              pvar = upd->vp[pg->imtrx][var];

              for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                source = -mp->surface_tension * fv->div_n * fv->n[a] * lsi->delta;

                source *= bf[var]->d_det_J_dm[b][j] * h3 + fv->dh3dmesh[b][j] * det_J;
                source *= wt * phi_i;
                source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

                source1 = -mp->surface_tension * fv->d_div_n_dmesh[b][j] * fv->n[a] * lsi->delta;
                source1 *= wt * phi_i * det_J * h3;
                source1 *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

                lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source + source1;
              }
            }
          }
        }
      }
    }
  }
  return (0);
}
int assemble_div_s_n_source(void) {
  int i, j, a, b, p, q, ii, ledof;
  int eqn, peqn, var, pvar;
  int dim;

  struct Basis_Functions *bfm;

  double wt, det_J, h3, phi_i, phi_j;
  double source, source1, source2, source3;

  eqn = R_MOMENTUM1;
  if (!pd->e[pg->imtrx][eqn]) {
    return (-1);
  }

  wt = fv->wt;
  h3 = fv->h3;

  if (ls->on_sharp_surf) /* sharp interface */
  {
    det_J = fv->sdet;
  } else /* diffuse interface */
  {
    det_J = bf[eqn]->detJ;
  }

  dim = pd->Num_Dim;

  /* finite difference calculation of path dependencies for
     subelement integration
   */
  if (ls->CalcSurfDependencies) {
    for (a = 0; a < WIM; a++) {
      eqn = R_MOMENTUM1 + a;
      peqn = upd->ep[pg->imtrx][eqn];
      bfm = bf[eqn];

      for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {

        ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][i];

        if (ei[pg->imtrx]->active_interp_ledof[ledof]) {

          ii = ei[pg->imtrx]->lvdof_to_row_lvdof[eqn][i];

          phi_i = bfm->phi[i];

          source = -mp->surface_tension * fv->div_s_n * fv->n[a] * lsi->delta;

          source *= phi_i * wt * det_J * h3;
          source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

          /* J_m_F
           */
          var = LS;
          pvar = upd->vp[pg->imtrx][var];

          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];
            lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source * phi_j;
          }
        }
      }
    }
    return (0);
  }

  /*
   * Wesiduals
   * ________________________________________________________________________________
   */

  if (af->Assemble_Residual) {
    for (a = 0; a < WIM; a++) {
      eqn = R_MOMENTUM1 + a;
      peqn = upd->ep[pg->imtrx][eqn];
      bfm = bf[eqn];

      for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {

        ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][i];

        if (ei[pg->imtrx]->active_interp_ledof[ledof]) {

          ii = ei[pg->imtrx]->lvdof_to_row_lvdof[eqn][i];

          phi_i = bfm->phi[i];

          source = 0.0;

          source = -mp->surface_tension * fv->div_s_n * fv->n[a] * lsi->delta;

          source *= phi_i * wt * det_J * h3;
          source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

          lec->R[LEC_R_INDEX(peqn, ii)] += source;
        }
      }
    }
  }

  /*
   * Yacobian terms...
   */

  if (af->Assemble_Jacobian) {
    for (a = 0; a < WIM; a++) {
      eqn = R_MOMENTUM1 + a;
      peqn = upd->ep[pg->imtrx][eqn];
      bfm = bf[eqn];

      for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {

        ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][i];

        if (ei[pg->imtrx]->active_interp_ledof[ledof]) {

          ii = ei[pg->imtrx]->lvdof_to_row_lvdof[eqn][i];

          phi_i = bfm->phi[i];

          /*
           * J_m_n
           */

          for (b = 0; b < dim; b++) {
            var = NORMAL1 + b;
            if (pd->v[pg->imtrx][var]) {

              pvar = upd->vp[pg->imtrx][var];

              for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                source1 = 0.0;
                for (p = 0; p < VIM; p++) {
                  source1 -= fv->n[p] * fv->grad_n[p][b] + fv->n[p] * fv->grad_n[b][p];
                }
                source1 *= fv->n[a] * bf[var]->phi[j];

                source2 = 0.0;
                for (p = 0; p < VIM; p++) {
                  for (q = 0; q < VIM; q++) {
                    source2 +=
                        (delta(p, q) - fv->n[p] * fv->n[q]) * bf[var]->grad_phi_e[j][b][p][q];
                  }
                }
                source2 *= fv->n[a];

                source3 = delta(a, b) * fv->div_s_n * bf[var]->phi[j];

                source = source1 + source2 + source3;
                source *= -mp->surface_tension * lsi->delta;
                source *= phi_i * wt * det_J * h3;
                source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

                lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
              }
            }
          }
          /*
           * J_m_F
           */
          var = FILL;
          if (pd->v[pg->imtrx][var]) {
            pvar = upd->vp[pg->imtrx][var];

            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              source = 0.0;

              source = -mp->surface_tension * fv->div_s_n * fv->n[a];
              source *= lsi->d_delta_dF[j];
              source *= phi_i * wt * det_J * h3;
              source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

              lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
            }
          }
        }
      }
    }
  }
  return (0);
}
int assemble_cap_hysing(double dt, double scale) {
  int i, j, a, b, p, q, k, ii, ledof;
  int eqn, peqn, var, pvar;

  struct Level_Set_Interface lsi_old_struct;
  struct Level_Set_Interface *lsi_old = &lsi_old_struct;

  struct Basis_Functions *bfm;
  dbl(*grad_phi_i_e_a)[DIM] = NULL;

  double wt, det_J, h3;

  double csf[DIM][DIM];

  double grad_s_v[DIM][DIM];
  //  double grad_s_phi_i_e_a[DIM][DIM];
  double d_grad_s_v_dv[DIM][DIM];

  double source;
  double diffusion;

  eqn = R_MOMENTUM1;
  if (!pd->e[pg->imtrx][eqn]) {
    return (0);
  }

  wt = fv->wt;
  h3 = fv->h3;
  if (ls->on_sharp_surf) /* sharp interface */
  {
    det_J = fv->sdet;
  } else /* diffuse interface */
  {
    det_J = bf[eqn]->detJ;
  }

  memset(csf, 0, sizeof(double) * DIM * DIM);

  continuous_surface_tension_old(mp->surface_tension, csf, lsi_old);

  for (p = 0; p < VIM; p++) {
    for (q = 0; q < VIM; q++) {
      grad_s_v[p][q] = fv->grad_v[p][q];
      for (k = 0; k < VIM; k++) {
        grad_s_v[p][q] -= (lsi_old->normal[p] * lsi_old->normal[k]) * fv->grad_v[k][q];
      }
    }
  }

  /* finite difference calculation of path dependencies for
     subelement integration
  */
  if (ls->CalcSurfDependencies) {
    GOMA_EH(GOMA_ERROR, "Calc surf dependencies not implemented");
  }

  /*
   * Residuals ____________________________________________________________________________
   */

  if (af->Assemble_Residual) {

    for (a = 0; a < WIM; a++) {
      eqn = R_MOMENTUM1 + a;
      peqn = upd->ep[pg->imtrx][eqn];
      bfm = bf[eqn];

      for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {

        ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][i];

        if (ei[pg->imtrx]->active_interp_ledof[ledof]) {

          ii = ei[pg->imtrx]->lvdof_to_row_lvdof[eqn][i];

          grad_phi_i_e_a = bfm->grad_phi_e[i][a];

          source = 0.;

          for (p = 0; p < VIM; p++) {
            for (q = 0; q < VIM; q++) {
              source += grad_phi_i_e_a[p][q] * csf[q][p];
            }
          }

          source *= -det_J * wt * h3;
          source *= pd->etm[pg->imtrx][eqn][LOG2_SOURCE];

          diffusion = 0;
          /*
           * 			grad_v[a][b] = d v_b
           *				       -----
           *				       d x_a
           */
          for (p = 0; p < VIM; p++) {
            for (q = 0; q < VIM; q++) {
              diffusion += bf[eqn]->grad_phi_e[i][a][p][q] * grad_s_v[p][q];
            }
          }

          diffusion *= -det_J * wt * h3 * (dt * mp->surface_tension * lsi_old->delta * scale);

          lec->R[LEC_R_INDEX(peqn, ii)] += source + diffusion;
        }
      }
    }
  }

  /*
   * Yacobian terms...
   */

  if (af->Assemble_Jacobian) {
    for (a = 0; a < WIM; a++) {
      eqn = R_MOMENTUM1 + a;
      peqn = upd->ep[pg->imtrx][eqn];
      bfm = bf[eqn];

      for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {

        ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][i];

        if (ei[pg->imtrx]->active_interp_ledof[ledof]) {

          ii = ei[pg->imtrx]->lvdof_to_row_lvdof[eqn][i];

          grad_phi_i_e_a = bfm->grad_phi_e[i][a];

          var = VELOCITY1;
          if (pd->v[pg->imtrx][var]) {
            for (b = 0; b < VIM; b++) {
              var = VELOCITY1 + b;
              pvar = upd->vp[pg->imtrx][var];

              for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                for (p = 0; p < VIM; p++) {
                  for (q = 0; q < VIM; q++) {
                    d_grad_s_v_dv[p][q] = bf[var]->grad_phi_e[j][b][p][q];
                    for (k = 0; k < VIM; k++) {
                      d_grad_s_v_dv[p][q] -=
                          lsi_old->normal[p] * lsi_old->normal[k] * bf[var]->grad_phi_e[j][b][k][q];
                    }
                  }
                }

                diffusion = 0;

                for (p = 0; p < VIM; p++) {
                  for (q = 0; q < VIM; q++) {
                    diffusion += bf[eqn]->grad_phi_e[i][a][p][q] * d_grad_s_v_dv[p][q];
                  }
                }

                diffusion *= -det_J * wt * h3 * (dt * mp->surface_tension * lsi_old->delta * scale);
                lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += diffusion;
              }
            }
          }

          var = MESH_DISPLACEMENT1;
          if (pd->v[pg->imtrx][var]) {
            GOMA_EH(GOMA_ERROR, "Jacobian terms for hysing capillary wrt mesh not implemented");
          }
        }
      }
    }
  }

  return 0;
}
int assemble_cap_denner_diffusion(double dt, double scale) {
  int i, j, a, b, p, q, k, ii, ledof;
  int eqn, peqn, var, pvar;

  double wt, det_J, h3;

  double csf[DIM][DIM];

  double grad_s_v[DIM][DIM];
  //  double grad_s_phi_i_e_a[DIM][DIM];
  double d_grad_s_v_dv[DIM][DIM];
  double d_grad_s_v_dF[DIM][DIM];
  //  double d_grad_s_v_dn[DIM][DIM];

  double diffusion;

  eqn = R_MOMENTUM1;
  if (!pd->e[pg->imtrx][eqn]) {
    return (0);
  }

  load_lsi(ls->Length_Scale);

  wt = fv->wt;
  h3 = fv->h3;
  if (ls->on_sharp_surf) /* sharp interface */
  {
    det_J = fv->sdet;
  } else /* diffuse interface */
  {
    det_J = bf[eqn]->detJ;
  }

  memset(csf, 0, sizeof(double) * DIM * DIM);

  for (p = 0; p < VIM; p++) {
    for (q = 0; q < VIM; q++) {
      grad_s_v[p][q] = fv->grad_v[p][q];
      for (k = 0; k < VIM; k++) {
        grad_s_v[p][q] -= (lsi->normal[p] * lsi->normal[k]) * fv->grad_v[k][q];
      }
    }
  }

  /* finite difference calculation of path dependencies for
     subelement integration
  */
  if (ls->CalcSurfDependencies) {
    GOMA_EH(GOMA_ERROR, "Calc surf dependencies not implemented");
  }

  /*
   * Residuals ____________________________________________________________________________
   */

  if (af->Assemble_Residual) {

    for (a = 0; a < WIM; a++) {
      eqn = R_MOMENTUM1 + a;
      peqn = upd->ep[pg->imtrx][eqn];

      for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {

        ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][i];

        if (ei[pg->imtrx]->active_interp_ledof[ledof]) {

          ii = ei[pg->imtrx]->lvdof_to_row_lvdof[eqn][i];

          diffusion = 0;

          for (p = 0; p < VIM; p++) {
            for (q = 0; q < VIM; q++) {
              diffusion += bf[eqn]->grad_phi_e[i][a][p][q] * grad_s_v[p][q];
            }
          }

          diffusion *= -det_J * wt * h3 * (dt * mp->surface_tension * lsi->delta * scale);

          lec->R[LEC_R_INDEX(peqn, ii)] += diffusion;
        }
      }
    }
  }

  /*
   * Yacobian terms...
   */

  if (af->Assemble_Jacobian) {
    for (a = 0; a < WIM; a++) {
      eqn = R_MOMENTUM1 + a;
      peqn = upd->ep[pg->imtrx][eqn];

      for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {

        ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][i];

        if (ei[pg->imtrx]->active_interp_ledof[ledof]) {

          ii = ei[pg->imtrx]->lvdof_to_row_lvdof[eqn][i];

          var = LS;
          if (pd->v[pg->imtrx][var]) {
            pvar = upd->vp[pg->imtrx][var];

            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {

              /* grouping lsi->delta into this now */
              for (p = 0; p < VIM; p++) {
                for (q = 0; q < VIM; q++) {
                  d_grad_s_v_dF[p][q] =
                      lsi->d_delta_dF[j] * fv->grad_v[p][q] + lsi->delta * fv->grad_v[p][q];
                  for (k = 0; k < VIM; k++) {
                    d_grad_s_v_dF[p][q] -=
                        lsi->d_delta_dF[j] * (lsi->normal[p] * lsi->normal[k]) * fv->grad_v[k][q];
                    d_grad_s_v_dF[p][q] -=
                        lsi->delta * (lsi->d_normal_dF[p][j] * lsi->normal[k]) * fv->grad_v[k][q];
                    d_grad_s_v_dF[p][q] -=
                        lsi->delta * (lsi->normal[p] * lsi->d_normal_dF[k][j]) * fv->grad_v[k][q];
                  }
                }
              }

              diffusion = 0;

              for (p = 0; p < VIM; p++) {
                for (q = 0; q < VIM; q++) {
                  diffusion += bf[eqn]->grad_phi_e[i][a][p][q] * d_grad_s_v_dF[p][q];
                }
              }

              diffusion *= -det_J * wt * h3 * (dt * mp->surface_tension * scale);
              lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += diffusion;
            }
          }

          var = VELOCITY1;
          if (pd->v[pg->imtrx][var]) {
            for (b = 0; b < VIM; b++) {
              var = VELOCITY1 + b;
              pvar = upd->vp[pg->imtrx][var];

              for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {

                for (p = 0; p < VIM; p++) {
                  for (q = 0; q < VIM; q++) {
                    d_grad_s_v_dv[p][q] = bf[var]->grad_phi_e[j][b][p][q];
                    for (k = 0; k < VIM; k++) {
                      d_grad_s_v_dv[p][q] -=
                          lsi->normal[p] * lsi->normal[k] * bf[var]->grad_phi_e[j][b][k][q];
                    }
                  }
                }

                diffusion = 0;

                for (p = 0; p < VIM; p++) {
                  for (q = 0; q < VIM; q++) {
                    diffusion += bf[eqn]->grad_phi_e[i][a][p][q] * d_grad_s_v_dv[p][q];
                  }
                }

                diffusion *= -det_J * wt * h3 * (dt * mp->surface_tension * lsi->delta * scale);
                lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += diffusion;
              }
            }
          }

          var = MESH_DISPLACEMENT1;
          if (pd->v[pg->imtrx][var]) {
            GOMA_EH(GOMA_ERROR,
                    "Jacobian terms for denner capillary diffusion wrt mesh not implemented");
          }
        }
      }
    }
  }

  return 0;
}
int assemble_cap_denner_diffusion_n(double dt, double scale) {
  int i, j, a, b, p, q, k, ii, ledof;
  int eqn, peqn, var, pvar;

  double wt, det_J, h3;

  double csf[DIM][DIM];

  double grad_s_v[DIM][DIM];
  //  double grad_s_phi_i_e_a[DIM][DIM];
  double d_grad_s_v_dv[DIM][DIM];
  double d_grad_s_v_dn[DIM][DIM];
  //  double d_grad_s_v_dn[DIM][DIM];

  double diffusion;

  eqn = R_MOMENTUM1;
  if (!pd->e[pg->imtrx][eqn]) {
    return (0);
  }

  load_lsi(ls->Length_Scale);

  wt = fv->wt;
  h3 = fv->h3;
  if (ls->on_sharp_surf) /* sharp interface */
  {
    det_J = fv->sdet;
  } else /* diffuse interface */
  {
    det_J = bf[eqn]->detJ;
  }

  memset(csf, 0, sizeof(double) * DIM * DIM);

  for (p = 0; p < VIM; p++) {
    for (q = 0; q < VIM; q++) {
      grad_s_v[p][q] = fv->grad_v[p][q];
      for (k = 0; k < VIM; k++) {
        grad_s_v[p][q] -= (fv->n[p] * fv->n[k]) * fv->grad_v[k][q];
      }
    }
  }

  /* finite difference calculation of path dependencies for
     subelement integration
  */
  if (ls->CalcSurfDependencies) {
    GOMA_EH(GOMA_ERROR, "Calc surf dependencies not implemented");
  }

  /*
   * Residuals ____________________________________________________________________________
   */

  if (af->Assemble_Residual) {

    for (a = 0; a < WIM; a++) {
      eqn = R_MOMENTUM1 + a;
      peqn = upd->ep[pg->imtrx][eqn];

      for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {

        ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][i];

        if (ei[pg->imtrx]->active_interp_ledof[ledof]) {

          ii = ei[pg->imtrx]->lvdof_to_row_lvdof[eqn][i];

          diffusion = 0;

          for (p = 0; p < VIM; p++) {
            for (q = 0; q < VIM; q++) {
              diffusion += bf[eqn]->grad_phi_e[i][a][p][q] * grad_s_v[p][q];
            }
          }

          diffusion *= -det_J * wt * h3 * (dt * mp->surface_tension * lsi->delta * scale);

          lec->R[LEC_R_INDEX(peqn, ii)] += diffusion;
        }
      }
    }
  }

  /*
   * Yacobian terms...
   */

  if (af->Assemble_Jacobian) {
    for (a = 0; a < WIM; a++) {
      eqn = R_MOMENTUM1 + a;
      peqn = upd->ep[pg->imtrx][eqn];

      for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {

        ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][i];

        if (ei[pg->imtrx]->active_interp_ledof[ledof]) {

          ii = ei[pg->imtrx]->lvdof_to_row_lvdof[eqn][i];

          var = LS;
          if (pd->v[pg->imtrx][var]) {
            pvar = upd->vp[pg->imtrx][var];

            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              diffusion = 0;

              for (p = 0; p < VIM; p++) {
                for (q = 0; q < VIM; q++) {
                  diffusion += bf[eqn]->grad_phi_e[i][a][p][q] * grad_s_v[p][q];
                }
              }

              diffusion *= lsi->d_delta_dF[j];

              diffusion *= -det_J * wt * h3 * (dt * mp->surface_tension * scale);
              lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += diffusion;
            }
          }

          var = VELOCITY1;
          if (pd->v[pg->imtrx][var]) {
            for (b = 0; b < WIM; b++) {
              var = VELOCITY1 + b;
              pvar = upd->vp[pg->imtrx][var];

              for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {

                for (p = 0; p < WIM; p++) {
                  for (q = 0; q < WIM; q++) {
                    d_grad_s_v_dv[p][q] = bf[var]->grad_phi_e[j][b][p][q];
                    for (k = 0; k < WIM; k++) {
                      d_grad_s_v_dv[p][q] -= fv->n[p] * fv->n[k] * bf[var]->grad_phi_e[j][b][k][q];
                    }
                  }
                }

                diffusion = 0;

                for (p = 0; p < WIM; p++) {
                  for (q = 0; q < WIM; q++) {
                    diffusion += bf[eqn]->grad_phi_e[i][a][p][q] * d_grad_s_v_dv[p][q];
                  }
                }

                diffusion *= -det_J * wt * h3 * (dt * mp->surface_tension * lsi->delta * scale);
                lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += diffusion;
              }
            }
          }

          var = R_NORMAL1;
          if (pd->v[pg->imtrx][var]) {
            for (b = 0; b < VIM; b++) {
              var = R_NORMAL1 + b;
              pvar = upd->vp[pg->imtrx][var];

              for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                for (p = 0; p < VIM; p++) {
                  for (q = 0; q < VIM; q++) {
                    d_grad_s_v_dn[p][q] = 0;
                    for (k = 0; k < VIM; k++) {
                      d_grad_s_v_dn[p][q] -= bf[var]->phi[j] * fv->n[k] * fv->grad_v[k][q];
                      d_grad_s_v_dn[p][q] -= fv->n[p] * bf[var]->phi[j] * fv->grad_v[k][q];
                    }
                  }
                }

                diffusion = 0;

                for (p = 0; p < VIM; p++) {
                  for (q = 0; q < VIM; q++) {
                    diffusion += bf[eqn]->grad_phi_e[i][a][p][q] * d_grad_s_v_dn[p][q];
                  }
                }

                diffusion *= -det_J * wt * h3 * (dt * mp->surface_tension * lsi->delta * scale);
                lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += diffusion;
              }
            }
          }

          var = MESH_DISPLACEMENT1;
          if (pd->v[pg->imtrx][var]) {
            GOMA_EH(GOMA_ERROR,
                    "Jacobian terms for denner capillary diffusion wrt mesh not implemented");
          }
        }
      }
    }
  }

  return 0;
}