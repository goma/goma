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

/*
 *   add user components to Residual and Jacobian below
 *
 */

if (af->Assemble_Residual) {
  for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
    phi_i = bf[eqn]->phi[i];

    /* Mass term */
    if (pd0->e[0][eqn] & T_MASS) {
      mass = phi_i * (fv->sh_u - 2 * fv->qs);
      mass *= pd0->etm[pg->imtrx][eqn][(LOG2_MASS)];
      diffusion = 0.0;
      diffusion *= pd0->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
    }
    res[i] += (mass + diffusion) * wt * h3 * det_J;
    lec->R[LEC_R_INDEX(peqn, i)] += res[i];
  }
}
/* Include Jacobian contributions from shell variables */
if (af->Assemble_Jacobian) {
  for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
    phi_i = bf[eqn]->phi[i];

    /* J_qs_qs:  Shell sensitivity */
    var = SHELL_USER;
    pvar = upd->vp[pg->imtrx][var];
    for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
      phi_j = bf[var]->phi[j];
      phi_t = phi_j;

      /* Mass term only */
      mass = 0.0;
      if (pd0->e[0][eqn] & T_MASS) {
        mass += phi_i * phi_t;
        mass *= pd0->etm[pg->imtrx][eqn][(LOG2_MASS)];
      }
      diffusion = 0.0;
      if (pd0->e[0][eqn] & T_MASS) {
        diffusion *= pd0->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
      }
      jac[i][pvar][j] += (mass + diffusion) * wt * h3 * det_J;
      lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += jac[i][pvar][j];
    }
    /* J_qs_qs:  Shell sensitivity */
    var = SURF_CHARGE;
    pvar = upd->vp[pg->imtrx][var];
    for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
      phi_j = bf[var]->phi[j];
      phi_t = phi_j;

      /* Mass term only */
      mass = 0.0;
      if (pd0->e[0][eqn] & T_MASS) {
        mass += phi_i * (-2 * phi_t);
        mass *= pd0->etm[pg->imtrx][eqn][(LOG2_MASS)];
      }
      diffusion = 0.0;
      if (pd0->e[0][eqn] & T_MASS) {
        diffusion *= pd0->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
      }
      jac[i][pvar][j] += (mass + diffusion) * wt * h3 * det_J;
      lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += jac[i][pvar][j];
    }
  }
}
