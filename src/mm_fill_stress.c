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

/* mm_fill_stress -- auxiliary routines helpful in matrix & residual assembly
 * for polymer stress equations
 *
 *   Copyright (C) 1993, 1994  Sandia National Laboratories
 */

/* Standard include files */

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>

/* GOMA include files */
#include "ac_stability.h"
#include "ac_stability_util.h"
#include "az_aztec.h"
#include "bc_colloc.h"
#include "el_elm.h"
#include "el_elm_info.h"
#include "el_geom.h"
#include "exo_struct.h"
#include "load_field_variables.h"
#include "mm_as.h"
#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "mm_bc.h"
#include "mm_eh.h"
#include "mm_fill_aux.h"
#include "mm_fill_fill.h"
#include "mm_fill_ls.h"
#include "mm_fill_solid.h"
#include "mm_fill_species.h"
#include "mm_fill_stabilization.h"
#include "mm_fill_terms.h"
#include "mm_fill_util.h"
#include "mm_mp.h"
#include "mm_mp_const.h"
#include "mm_mp_structs.h"
#include "mm_post_def.h"
#include "mm_unknown_map.h"
#include "mm_viscosity.h"
#include "polymer_time_const.h"
#include "rf_allo.h"
#include "rf_bc.h"
#include "rf_bc_const.h"
#include "rf_fem.h"
#include "rf_fem_const.h"
#include "rf_node_const.h"
#include "rf_solver.h"
#include "rf_solver_const.h"
#include "rf_vars_const.h"
#include "sl_util_structs.h"
#include "std.h"

#define GOMA_MM_FILL_STRESS_C
#include "mm_fill_stress.h"

extern struct Boundary_Condition *inlet_BC[MAX_VARIABLE_TYPES + MAX_CONC];

// direct call to a fortran LAPACK eigenvalue routine
extern FSUB_TYPE dsyev_(char *JOBZ,
                        char *UPLO,
                        int *N,
                        double *A,
                        int *LDA,
                        double *W,
                        double *WORK,
                        int *LWORK,
                        int *INFO,
                        int len_jobz,
                        int len_uplo);

/*  _______________________________________________________________________  */

/* assemble_stress -- assemble terms (Residual &| Jacobian) for polymer stress eqns
 *
 * in:
 * 	ei -- pointer to Element Indeces	structure
 *	pd -- pointer to Problem Description	structure
 *	af -- pointer to Action Flag		structure
 *	bf -- pointer to Basis Function		structure
 *	fv -- pointer to Field Variable		structure
 *	cr -- pointer to Constitutive Relation	structure
 *	md -- pointer to Mesh Derivative	structure
 *	me -- pointer to Material Entity	structure
 * 	ija -- vector of pointers into the a matrix
 * 	a  -- global Jacobian matrix
 * 	R  -- global residual vector
 *
 * out:
 *	a   -- gets loaded up with proper contribution
 *	lec -- gets loaded up with local contributions to resid, Jacobian
 * 	r  -- residual RHS vector
 *
 * Created:	Wed Dec  8 14:03:06 MST 1993 pasacki@sandia.gov
 *
 * Revised:	Sun Feb 27 06:53:12 MST 1994 pasacki@sandia.gov
 *
 *
 * Note: currently we do a "double load" into the addresses in the global
 *       "a" matrix and resid vector held in "esp", as well as into the
 *       local accumulators in "lec".
 *
 */

int assemble_stress(dbl tt, /* parameter to vary time integration from
                               explicit (tt = 1) to implicit (tt = 0) */
                    dbl dt, /* current time step size */
                    dbl h[DIM],
                    dbl hh[DIM][DIM],
                    dbl dh_dxnode[DIM][MDE],
                    dbl vcent[DIM], /* average element velocity, which is the
                                     * centroid velocity for Q2 and the average of
                                     * the vertices for Q1. It comes from the
                                     * routine "element_velocity."               */
                    dbl dvc_dnode[DIM][MDE]) {
  int dim, p, q, r, a, b, w = -1;

  int mode; /*counter for viscoelastic modes */

  int eqn, var;
  int peqn, pvar;

  int i, j, status, imtrx;

  dbl v[DIM];      /* Velocity field. */
  dbl x_dot[DIM];  /* current position field derivative wrt time. */
  dbl h3;          /* Volume element (scale factors). */
  dbl dh3dmesh_pj; /* Sensitivity to (p,j) mesh dof. */

  dbl grad_v[DIM][DIM];
  dbl gamma[DIM][DIM]; /* Shear-rate tensor based on velocity */
  dbl det_J;           /* determinant of element Jacobian */

  dbl d_det_J_dmesh_pj; /* for specific (p,j) mesh dof */

  dbl mass; /* For terms and their derivatives */
  dbl mass_a, mass_b, mass_c;
  dbl advection;
  dbl advection_a, advection_b, advection_c, advection_d;
  dbl diffusion;
  dbl source;
  dbl source1;
  dbl source_a, source_b, source_c;

  /*
   *
   * Note how carefully we avoid refering to d(phi[i])/dx[j] and refer instead
   * to the j-th component of grad_phi[j][i] so that this vector can be loaded
   * up with components that may be different in non Cartesian coordinate
   * systems.
   *
   * We will, however, insist on *orthogonal* coordinate systems, even if we
   * might permit them to be curvilinear.
   *
   * Assume all components of velocity are interpolated with the same kind
   * of basis function.
   */

  /*
   * Petrov-Galerkin weighting functions for i-th and ab-th stress residuals
   * and some of their derivatives...
   */

  dbl wt_func;

  /*
   * Interpolation functions for variables and some of their derivatives.
   */

  dbl phi_j;
  dbl wt;

  /* Variables for stress */

  int R_s[MAX_MODES][DIM][DIM];
  int v_s[MAX_MODES][DIM][DIM];
  int v_g[DIM][DIM];

  dbl s[DIM][DIM];     /* stress tensor */
  dbl s_dot[DIM][DIM]; /* stress tensor from last time step */
  dbl grad_s[DIM][DIM][DIM];
  dbl d_grad_s_dmesh[DIM][DIM][DIM][DIM]
                    [MDE]; /* derivative of grad of stress tensor for mode ve_mode */

  dbl g[DIM][DIM];      /* velocity gradient tensor */
  dbl gt[DIM][DIM];     /* transpose of velocity gradient tensor */
  dbl g_dot[DIM][DIM];  /* velocity gradient tensor time derivative */
  dbl gt_dot[DIM][DIM]; /* transpose of velocity gradient tensor time derivative */

  /* dot product tensors */

  dbl g_dot_g[DIM][DIM];
  dbl g_dot_gt[DIM][DIM];
  dbl gt_dot_g[DIM][DIM];
  dbl gt_dot_gt[DIM][DIM];
  dbl s_dot_s[DIM][DIM];
  dbl s_dot_g[DIM][DIM];
  dbl s_dot_gt[DIM][DIM];
  dbl g_dot_s[DIM][DIM];
  dbl gt_dot_s[DIM][DIM];

  /* polymer viscosity and derivatives */
  dbl mup;
  VISCOSITY_DEPENDENCE_STRUCT d_mup_struct;
  VISCOSITY_DEPENDENCE_STRUCT *d_mup = &d_mup_struct;
  POLYMER_TIME_CONST_DEPENDENCE_STRUCT d_lam_struct;
  POLYMER_TIME_CONST_DEPENDENCE_STRUCT *d_lam = &d_lam_struct;

  dbl d_mup_dv_pj;
  dbl d_mup_dmesh_pj;

  /* constitutive equation parameters */
  dbl alpha;      /* This is the Geisekus mobility parameter */
  dbl lambda = 0; /* polymer relaxation constant */

  /* advective terms are precalculated */
  dbl v_dot_del_g[DIM][DIM];
  dbl x_dot_del_g[DIM][DIM];
  dbl v_dot_del_s[DIM][DIM];
  dbl x_dot_del_s[DIM][DIM];

  dbl d_xdotdelg_dm;
  dbl d_xdotdelgt_dm;
  dbl d_xdotdels_dm;

  dbl d_vdotdelg_dm;
  dbl d_vdotdelgt_dm;
  dbl d_vdotdels_dm;

  /* SUPG variables */
  dbl h_elem = 0, h_elem_inv = 0, h_elem_deriv = 0;
  dbl supg = 0;

  status = 0;

  eqn = R_STRESS11;

  /*
   * Bail out fast if there's nothing to do...
   */

  if (!pd->e[pg->imtrx][eqn]) {
    return (status);
  }

  /*
   * Unpack variables from structures for local convenience...
   */

  dim = pd->Num_Dim;

  wt = fv->wt;

  det_J = bf[eqn]->detJ; /* Really, ought to be mesh eqn. */

  h3 = fv->h3; /* Differential volume element (scales). */

  v_g[0][0] = VELOCITY_GRADIENT11;
  v_g[0][1] = VELOCITY_GRADIENT12;
  v_g[1][0] = VELOCITY_GRADIENT21;
  v_g[1][1] = VELOCITY_GRADIENT22;
  v_g[0][2] = VELOCITY_GRADIENT13;
  v_g[1][2] = VELOCITY_GRADIENT23;
  v_g[2][0] = VELOCITY_GRADIENT31;
  v_g[2][1] = VELOCITY_GRADIENT32;
  v_g[2][2] = VELOCITY_GRADIENT33;

  /*
   * Field variables...
   */

  for (a = 0; a < WIM; a++) {
    v[a] = fv->v[a];

    /* note, these are zero for steady calculations */
    x_dot[a] = 0.0;
    for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
      if (pd->TimeIntegration != STEADY && pd->v[imtrx][MESH_DISPLACEMENT1 + a]) {
        x_dot[a] = fv_dot->x[a];
      }
    }
  }

  /*
   * In Cartesian coordinates, this velocity gradient tensor will
   * have components that are...
   *
   * 			grad_v[a][b] = d v_b
   *				       -----
   *				       d x_a
   */

  for (a = 0; a < VIM; a++) {
    for (b = 0; b < VIM; b++) {
      grad_v[a][b] = fv->grad_v[a][b];
    }
  }

  /* load up shearrate tensor based on velocity */
  for (a = 0; a < VIM; a++) {
    for (b = 0; b < VIM; b++) {
      gamma[a][b] = grad_v[a][b] + grad_v[b][a];
    }
  }

  for (a = 0; a < VIM; a++) {
    for (b = 0; b < VIM; b++) {
      g[a][b] = fv->G[a][b];
      gt[b][a] = g[a][b];
      if (pd->TimeIntegration != STEADY) {
        g_dot[a][b] = fv_dot->G[a][b];
        gt_dot[b][a] = g_dot[a][b];
      } else {
        g_dot[a][b] = 0.0;
        gt_dot[a][b] = 0.0;
      }
    }
  }

  if (vn->wt_funcModel == GALERKIN) {
    supg = 0.;
  } else if (vn->wt_funcModel == SUPG) {
    supg = vn->wt_func;
  }

  if (supg != 0.) {
    h_elem = 0.;
    for (p = 0; p < dim; p++) {
      h_elem += vcent[p] * vcent[p] * h[p];
    }
    h_elem = sqrt(h_elem) / 2.;
    if (h_elem == 0.) {
      h_elem_inv = 1.;
    } else {
      h_elem_inv = 1. / h_elem;
    }
  }
  /* end Petrov-Galerkin addition */

  /* get tensor dot products for future use */

  (void)tensor_dot(g, g, g_dot_g, VIM);
  (void)tensor_dot(g, gt, g_dot_gt, VIM);
  (void)tensor_dot(gt, g, gt_dot_g, VIM);
  (void)tensor_dot(gt, gt, gt_dot_gt, VIM);

  (void)stress_eqn_pointer(v_s);
  (void)stress_eqn_pointer(R_s);

  /* Begin loop over modes */
  for (mode = 0; mode < vn->modes; mode++) {

    load_modal_pointers(mode, tt, dt, s, s_dot, grad_s, d_grad_s_dmesh);

    /* precalculate advective terms of form (v dot del tensor)*/

    for (a = 0; a < VIM; a++) {
      for (b = 0; b < VIM; b++) {
        v_dot_del_g[a][b] = 0.;
        x_dot_del_g[a][b] = 0.;
        v_dot_del_s[a][b] = 0.;
        x_dot_del_s[a][b] = 0.;
        for (q = 0; q < WIM; q++) {
          v_dot_del_g[a][b] += v[q] * fv->grad_G[q][a][b];
          x_dot_del_g[a][b] += x_dot[q] * fv->grad_G[q][a][b];
          v_dot_del_s[a][b] += v[q] * grad_s[q][a][b];
          x_dot_del_s[a][b] += x_dot[q] * grad_s[q][a][b];
        }
      }
    }

    /*
     * Stress tensor...(Note "anti-BSL" sign convention on deviatoric stress)
     */

    /* get polymer viscosity */
    mup = viscosity(ve[mode]->gn, gamma, d_mup);

    /* get Geisekus mobility parameter */
    alpha = ve[mode]->alpha;

    /* get time constant */
    lambda = polymer_time_const(ve[mode]->time_const_st, gamma, d_lam);

    /* get tensor dot products for future use */

    (void)tensor_dot(s, s, s_dot_s, VIM);
    (void)tensor_dot(s, g, s_dot_g, VIM);
    (void)tensor_dot(s, gt, s_dot_gt, VIM);
    (void)tensor_dot(gt, s, gt_dot_s, VIM);
    (void)tensor_dot(g, s, g_dot_s, VIM);

    /*
     * Residuals_________________________________________________________________
     */

    if (af->Assemble_Residual) {
      /*
       * Assemble each component "ab" of the polymer stress equation...
       */
      for (a = 0; a < VIM; a++) {
        for (b = 0; b < VIM; b++) {

          if (a <= b) /* since the stress tensor is symmetric, only assemble the upper half */
          {
            eqn = R_s[mode][a][b];

            peqn = upd->ep[pg->imtrx][eqn];

            /*
             * In the element, there will be contributions to this many equations
             * based on the number of degrees of freedom...
             */

            for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
              wt_func = bf[eqn]->phi[i];
              /* add Petrov-Galerkin terms as necessary */
              if (supg != 0.) {
                for (w = 0; w < dim; w++) {
                  wt_func += supg * h_elem * v[w] * bf[eqn]->grad_phi[i][w];
                }
              }

              mass = 0.;

              if (pd->TimeIntegration != STEADY) {
                if (pd->e[pg->imtrx][eqn] & T_MASS) {
                  mass = s_dot[a][b] + mup * (g_dot[a][b] + gt_dot[a][b]);
                  mass *= wt_func * lambda * det_J * wt;
                  mass *= h3;
                  mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
                }
              }

              advection = 0.;
              if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
                if (DOUBLE_NONZERO(lambda)) {

                  advection += (v_dot_del_s[a][b] - x_dot_del_s[a][b] +
                                mup * (v_dot_del_g[a][b] - x_dot_del_g[a][b] + v_dot_del_g[b][a] -
                                       x_dot_del_g[b][a]));
                  advection -= (mup * (2. * gt_dot_g[a][b] + gt_dot_gt[a][b] + g_dot_g[a][b]) +
                                gt_dot_s[a][b] + s_dot_g[a][b]);

                  advection *= wt_func * lambda * det_J * wt * h3;
                  advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
                }
              }

              diffusion = 0.;
              if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {
                /* add SU term in here when appropriate */
                diffusion += 0.;
                diffusion *= -wt_func * det_J * wt * h3;
                diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
              }

              /*
               * Source term...
               */

              source = 0.;
              source1 = 0.;

              if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
                source += s[a][b];
                if (DOUBLE_NONZERO(alpha)) {
                  source1 +=
                      (s_dot_g[a][b] + s_dot_gt[a][b] + g_dot_s[a][b] + gt_dot_s[a][b] +
                       s_dot_s[a][b] / mup +
                       mup * (g_dot_g[a][b] + gt_dot_gt[a][b] + gt_dot_g[a][b] + g_dot_gt[a][b]));

                  source1 *= alpha * lambda;
                  source += source1;
                }

                source *= wt_func * det_J * h3 * wt;

                source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
              }

              /*
               * Add contributions to this residual (globally into Resid, and
               * locally into an accumulator)
               */

              lec->R[LEC_R_INDEX(peqn, i)] += mass + advection + diffusion + source;
            }
          }
        }
      }
    }

    /*
     * Jacobian terms...
     */

    if (af->Assemble_Jacobian) {
      for (a = 0; a < VIM; a++) {
        for (b = 0; b < VIM; b++) {
          if (a <= b) /* since the stress tensor is symmetric, only assemble the upper half */
          {

            eqn = R_s[mode][a][b];
            peqn = upd->ep[pg->imtrx][eqn];

            for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {

              wt_func = bf[eqn]->phi[i];
              /* add Petrov-Galerkin terms as necessary */
              if (supg != 0.) {
                for (w = 0; w < dim; w++) {
                  wt_func += supg * h_elem * v[w] * bf[eqn]->grad_phi[i][w];
                }
              }

              /*
               * Set up some preliminaries that are needed for the (a,i)
               * equation for bunches of (b,j) column variables...
               */

              /*
               * J_S_T
               */

              var = TEMPERATURE;
              if (pd->v[pg->imtrx][var]) {
                pvar = upd->vp[pg->imtrx][var];
                for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                  phi_j = bf[var]->phi[j];

                  mass = 0.;

                  if (pd->TimeIntegration != STEADY) {
                    if (pd->e[pg->imtrx][eqn] & T_MASS) {
                      mass = d_mup->T[j] * (g_dot[a][b] + gt_dot[a][b]);
                      mass *= wt_func;
                      mass *= h3 * lambda * det_J * wt;
                      mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
                    }
                  }

                  advection = 0.;
                  if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
                    if (DOUBLE_NONZERO(lambda)) {
                      advection += (v_dot_del_g[a][b] - x_dot_del_g[a][b] + v_dot_del_g[b][a] -
                                    x_dot_del_g[b][a]) -
                                   2. * gt_dot_g[a][b] - gt_dot_gt[a][b] - g_dot_g[a][b];
                      advection *= wt_func * d_mup->T[j] * lambda * det_J * wt * h3;
                      advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
                    }
                  }

                  diffusion = 0.;
                  if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {
                    diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
                  }

                  source = 0.;
                  source1 = 0.;
                  if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
                    if (DOUBLE_NONZERO(alpha)) {
                      source1 += (-s_dot_s[a][b] / (mup * mup) + (g_dot_g[a][b] + gt_dot_gt[a][b] +
                                                                  gt_dot_g[a][b] + g_dot_gt[a][b]));
                      source1 *= lambda * alpha * d_mup->T[j];
                      source += source1;
                    }
                    source *= wt_func * det_J * wt * h3;
                    source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
                  }

                  lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += mass + advection + diffusion + source;
                }
              }

              /*
               * J_S_v
               */
              for (p = 0; p < WIM; p++) {
                var = VELOCITY1 + p;
                if (pd->v[pg->imtrx][var]) {
                  pvar = upd->vp[pg->imtrx][var];
                  for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                    phi_j = bf[var]->phi[j];
                    d_mup_dv_pj = d_mup->v[p][j];

                    mass = 0.;

                    if (pd->TimeIntegration != STEADY) {
                      if (pd->e[pg->imtrx][eqn] & T_MASS) {
                        mass = d_mup_dv_pj * (g_dot[a][b] + gt_dot[a][b]) * wt_func;

                        if (supg != 0.) {
                          mass_a = supg * h_elem * phi_j * bf[eqn]->grad_phi[i][p];

                          for (w = 0; w < dim; w++) {
                            mass_a += supg * vcent[p] * dvc_dnode[p][j] * h[p] * h_elem_inv / 4. *
                                      v[w] * bf[eqn]->grad_phi[i][w];
                          }

                          mass_a *= s_dot[a][b] + mup * (g_dot[a][b] + gt_dot[a][b]);
                          mass += mass_a;
                        }

                        mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)] * lambda * det_J * wt * h3;
                      }
                    }

                    advection = 0.;

                    if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
                      if (DOUBLE_NONZERO(lambda)) {
                        advection_a = phi_j * (grad_s[p][a][b] +
                                               mup * (fv->grad_G[p][a][b] + fv->grad_G[p][b][a]));

                        advection_a *= wt_func;
                        advection_b = (v_dot_del_g[a][b] - x_dot_del_g[a][b] + v_dot_del_g[b][a] -
                                       x_dot_del_g[b][a]);

                        advection_b -= 2. * gt_dot_g[a][b] + gt_dot_gt[a][b] + g_dot_g[a][b];

                        advection_b *= wt_func * d_mup_dv_pj;

                        advection_c = 0.;
                        /* Petrov-Galerkin term */
                        if (supg != 0.) {
                          advection_c = supg * h_elem * phi_j * bf[eqn]->grad_phi[i][p];
                          for (w = 0; w < dim; w++) {
                            advection_c += supg * vcent[p] * h[p] * dvc_dnode[p][j] * h_elem_inv /
                                           4. * v[w] * bf[eqn]->grad_phi[i][w];
                          }

                          advection_c *=
                              (v_dot_del_s[a][b] - x_dot_del_s[a][b] +
                               mup * (v_dot_del_g[a][b] - x_dot_del_g[a][b] + v_dot_del_g[b][a] -
                                      x_dot_del_g[b][a]) -
                               mup * (2. * gt_dot_g[a][b] + gt_dot_gt[a][b] + g_dot_g[a][b]) -
                               gt_dot_s[a][b] - s_dot_g[a][b]);
                        }

                        advection = advection_a + advection_b + advection_c;
                        advection *= lambda * det_J * wt * h3;
                        advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
                      }
                    }

                    diffusion = 0.;

                    if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {
                      /* add SU term in here when appropriate */

                      diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
                    }

                    source = 0.;

                    if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
                      source_a = 0.;
                      if (DOUBLE_NONZERO(alpha)) {
                        source_a = -s_dot_s[a][b] / (mup * mup) + g_dot_g[a][b] + gt_dot_gt[a][b] +
                                   gt_dot_g[a][b] + g_dot_gt[a][b];
                        source_a *= wt_func * alpha * lambda * d_mup_dv_pj;
                      }

                      source_b = 0.;
                      if (supg != 0.) {
                        source_b = supg * h_elem * phi_j * bf[eqn]->grad_phi[i][p];

                        for (w = 0; w < dim; w++) {
                          source_b += supg * vcent[p] * dvc_dnode[p][j] * h[p] * h_elem_inv / 4. *
                                      v[w] * bf[eqn]->grad_phi[i][w];
                        }

                        source_b *= s[a][b] + (s_dot_g[a][b] + s_dot_gt[a][b] + g_dot_s[a][b] +
                                               gt_dot_s[a][b] + s_dot_s[a][b] / mup +
                                               mup * (g_dot_g[a][b] + gt_dot_gt[a][b] +
                                                      gt_dot_g[a][b] + g_dot_gt[a][b])) *
                                                  alpha * lambda;
                      }

                      source = source_a + source_b;
                      source *= det_J * wt * h3;
                      source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
                    }

                    lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += mass + advection + diffusion + source;
                  }
                }
              }

              /*
               * J_S_c
               */
              var = MASS_FRACTION;
              if (pd->v[pg->imtrx][var]) {
                pvar = MAX_PROB_VAR + w;
                for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                  phi_j = bf[var]->phi[j];

                  for (w = 0; w < pd->Num_Species_Eqn; w++) {

                    mass = 0.;

                    if (pd->TimeIntegration != STEADY) {
                      if (pd->e[pg->imtrx][eqn] & T_MASS) {
                        mass = d_mup->C[w][j] * (g_dot[a][b] + gt_dot[a][b]);
                        mass *= wt_func * lambda * det_J * wt;
                        mass *= h3;
                        mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
                      }
                    }

                    advection = 0.;
                    if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
                      if (DOUBLE_NONZERO(lambda)) {
                        advection += (v_dot_del_g[a][b] - x_dot_del_g[a][b] + v_dot_del_g[b][a] -
                                      x_dot_del_g[b][a]) -
                                     2. * gt_dot_g[a][b] - gt_dot_gt[a][b] - g_dot_g[a][b];
                        advection *= d_mup->C[w][j] * wt_func * lambda * det_J * wt * h3;
                        advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
                      }
                    }

                    diffusion = 0.;
                    if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {
                      diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
                    }

                    source = 0.;

                    if (pd->e[pg->imtrx][eqn] & T_SOURCE) {

                      if (DOUBLE_NONZERO(alpha)) {
                        source += -s_dot_s[a][b] / (mup * mup) + g_dot_g[a][b] + gt_dot_gt[a][b] +
                                  gt_dot_g[a][b] + g_dot_gt[a][b];
                        source *= alpha * lambda * d_mup->C[w][j];
                      }

                      source *= wt_func * det_J * wt * h3;
                      source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
                    }

                    if (w > 1) {
                      GOMA_EH(GOMA_ERROR, "Need more arrays for each species.");
                    }

                    lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += mass + advection + diffusion + source;
                  }
                }
              }

              /*
               * J_S_P
               */
              var = PRESSURE;
              if (pd->v[pg->imtrx][var]) {
                pvar = upd->vp[pg->imtrx][var];
                for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                  phi_j = bf[var]->phi[j];

                  mass = 0.;
                  if (pd->TimeIntegration != STEADY) {
                    if (pd->e[pg->imtrx][eqn] & T_MASS) {
                      mass = d_mup->P[j] * (g_dot[a][b] + gt_dot[a][b]);
                      mass *= wt_func * lambda * det_J * wt;
                      mass *= h3;
                      mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
                    }
                  }

                  advection = 0.;
                  if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
                    if (DOUBLE_NONZERO(lambda)) {

                      advection += (v_dot_del_g[a][b] - x_dot_del_g[a][b] + v_dot_del_g[b][a] -
                                    x_dot_del_g[b][a]) -
                                   2. * gt_dot_g[a][b] - gt_dot_gt[a][b] - g_dot_g[a][b];
                      advection *= wt_func * lambda * det_J * wt * h3 * d_mup->P[j];
                      advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
                    }
                  }

                  diffusion = 0.;

                  if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {
                    diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
                  }

                  source = 0.;

                  if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
                    if (DOUBLE_NONZERO(alpha)) {
                      source += (-s_dot_s[a][b] / (mup * mup) + (g_dot_g[a][b] + gt_dot_gt[a][b] +
                                                                 gt_dot_g[a][b] + g_dot_gt[a][b]));
                      source *= d_mup->P[j] * alpha * lambda;
                    }
                    source *= wt_func * det_J * wt * h3;
                    source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
                  }

                  lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += mass + advection + diffusion + source;
                }
              }

              /*
               * J_S_d
               */
              for (p = 0; p < dim; p++) {
                var = MESH_DISPLACEMENT1 + p;
                if (pd->v[pg->imtrx][var]) {
                  pvar = upd->vp[pg->imtrx][var];
                  for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                    phi_j = bf[var]->phi[j];
                    d_det_J_dmesh_pj = bf[eqn]->d_det_J_dm[p][j];
                    dh3dmesh_pj = fv->dh3dq[p] * bf[var]->phi[j];
                    d_mup_dmesh_pj = d_mup->X[p][j];
                    if (supg != 0.) {
                      h_elem_deriv = 0.;
                      for (q = 0; q < dim; q++) {
                        h_elem_deriv +=
                            hh[q][p] * vcent[q] * vcent[q] * dh_dxnode[q][j] * h_elem_inv / 4.;
                      }
                    }

                    mass = 0.;
                    mass_a = 0.;
                    mass_b = 0.;
                    if (pd->TimeIntegration != STEADY) {
                      if (pd->e[pg->imtrx][eqn] & T_MASS) {
                        mass_a = s_dot[a][b] + mup * (g_dot[a][b] + gt_dot[a][b]);
                        mass_a *= wt_func * (d_det_J_dmesh_pj * h3 + det_J * dh3dmesh_pj);

                        mass_b = d_mup_dmesh_pj * (g_dot[a][b] + gt_dot[a][b]);
                        mass_b *= wt_func * h3 * det_J;

                        mass_c = 0.;
                        if (supg != 0.) {
                          for (w = 0; w < dim; w++) {
                            mass_c +=
                                supg * (h_elem * v[w] * bf[eqn]->d_grad_phi_dmesh[i][w][p][j] +
                                        h_elem_deriv * v[w] * bf[eqn]->grad_phi[i][w]);
                          }
                          mass_c *= (s_dot[a][b] + mup * (g_dot[a][b] + gt_dot[a][b])) * h3 * det_J;
                        }

                        mass = mass_a + mass_b + mass_c;
                        mass *= lambda * wt;
                        mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
                      }
                    }

                    advection = 0.;

                    if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
                      if (DOUBLE_NONZERO(lambda)) {
                        /*
                         * Four parts:
                         *    advection_a =
                         *    	Int ( ea.(v-xdot).d(Vv)/dmesh h3 |Jv| )
                         *
                         *    advection_b =
                         *  (i)	Int ( ea.(v-xdot).Vv h3 d(|Jv|)/dmesh )
                         *  (ii)  Int ( ea.(v-xdot).d(Vv)/dmesh h3 |Jv| )
                         *  (iii) Int ( ea.(v-xdot).Vv dh3/dmesh |Jv|   )
                         *
                         * For unsteady problems, we have an
                         * additional term
                         *
                         *    advection_c =
                         *    	Int ( ea.d(v-xdot)/dmesh.Vv h3 |Jv| )
                         */

                        advection_a = 0.;

                        advection_a += (v_dot_del_s[a][b] - x_dot_del_s[a][b] +
                                        mup * (v_dot_del_g[a][b] - x_dot_del_g[a][b] +
                                               v_dot_del_g[b][a] - x_dot_del_g[b][a]));
                        advection_a -=
                            (mup * (2. * gt_dot_g[a][b] + gt_dot_gt[a][b] + g_dot_g[a][b]) +
                             gt_dot_s[a][b] + s_dot_g[a][b]);

                        advection_a *= wt_func * (d_det_J_dmesh_pj * h3 + det_J * dh3dmesh_pj);

                        advection_b = 0.;
                        d_vdotdelg_dm = 0.;
                        d_vdotdelgt_dm = 0.;
                        d_vdotdels_dm = 0.;
                        for (q = 0; q < WIM; q++) {
                          d_vdotdels_dm += (v[q] - x_dot[q]) * d_grad_s_dmesh[q][a][b][p][j];
                          d_vdotdelg_dm += (v[q] - x_dot[q]) * fv->d_grad_G_dmesh[q][a][b][p][j];
                          d_vdotdelgt_dm += (v[q] - x_dot[q]) * fv->d_grad_G_dmesh[q][b][a][p][j];
                        }

                        advection_b += (d_vdotdels_dm + mup * (d_vdotdelg_dm + d_vdotdelgt_dm));
                        advection_b += d_mup_dmesh_pj *
                                       ((v_dot_del_g[a][b] - x_dot_del_g[a][b] + v_dot_del_g[b][a] -
                                         x_dot_del_g[b][a]) -
                                        2. * gt_dot_g[a][b] - gt_dot_gt[a][b] - g_dot_g[a][b]);
                        advection_b *= wt_func * det_J * h3;

                        advection_c = 0.;
                        if (pd->TimeIntegration != STEADY) {
                          if (pd->e[pg->imtrx][eqn] & T_MASS) {
                            d_xdotdels_dm = (1. + 2. * tt) * phi_j / dt * grad_s[p][a][b];
                            d_xdotdelg_dm = (1. + 2. * tt) * phi_j / dt * fv->grad_G[p][a][b];
                            d_xdotdelgt_dm = (1. + 2. * tt) * phi_j / dt * fv->grad_G[p][b][a];

                            advection_c -= (d_xdotdels_dm + mup * (d_xdotdelg_dm + d_xdotdelgt_dm));

                            advection_c *= wt_func * h3 * det_J;
                          }
                        }

                        advection_d = 0.;
                        if (supg != 0.) {
                          for (w = 0; w < dim; w++) {
                            advection_d +=
                                supg * (h_elem * v[w] * bf[eqn]->d_grad_phi_dmesh[i][w][p][j] +
                                        h_elem_deriv * v[w] * bf[eqn]->grad_phi[i][w]);
                          }

                          advection_d *=
                              (v_dot_del_s[a][b] - x_dot_del_s[a][b] +
                               mup * (v_dot_del_g[a][b] - x_dot_del_g[a][b] + v_dot_del_g[b][a] -
                                      x_dot_del_g[b][a]) -
                               mup * (2. * gt_dot_g[a][b] + gt_dot_gt[a][b] + g_dot_g[a][b]) -
                               gt_dot_s[a][b] - s_dot_g[a][b]) *
                              det_J * h3;
                        }

                        advection = advection_a + advection_b + advection_c + advection_d;

                        advection *= wt * lambda * pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
                      }
                    }

                    /*
                     * Diffusion...
                     */

                    diffusion = 0.;

                    if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {
                      diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
                    }

                    /*
                     * Source term...
                     */

                    source = 0.;

                    if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
                      source_a = s[a][b];
                      source_b = 0.;

                      if (DOUBLE_NONZERO(alpha)) {
                        source_a += (s_dot_g[a][b] + s_dot_gt[a][b] + g_dot_s[a][b] +
                                     gt_dot_s[a][b] + s_dot_s[a][b] / mup +
                                     mup * (g_dot_g[a][b] + gt_dot_gt[a][b] + gt_dot_g[a][b] +
                                            g_dot_gt[a][b])) *
                                    alpha * lambda;
                        source_b +=
                            -s_dot_s[a][b] / (mup * mup) +
                            (g_dot_g[a][b] + gt_dot_gt[a][b] + gt_dot_g[a][b] + g_dot_gt[a][b]);
                        source_b *= d_mup_dmesh_pj * alpha * lambda;
                      }

                      source_a *= wt_func * (d_det_J_dmesh_pj * h3 + det_J * dh3dmesh_pj);

                      source_b *= wt_func * det_J * h3;

                      source_c = 0.;
                      if (supg != 0.) {

                        for (w = 0; w < dim; w++) {
                          source_c +=
                              supg * (h_elem * v[w] * bf[eqn]->d_grad_phi_dmesh[i][w][p][j] +
                                      h_elem_deriv * v[w] * bf[eqn]->grad_phi[i][w]);
                        }
                        source_c *= (s[a][b] + alpha * lambda *
                                                   (s_dot_g[a][b] + s_dot_gt[a][b] + g_dot_s[a][b] +
                                                    gt_dot_s[a][b] + s_dot_s[a][b] / mup +
                                                    mup * (g_dot_g[a][b] + gt_dot_gt[a][b] +
                                                           gt_dot_g[a][b] + g_dot_gt[a][b]))) *
                                    det_J * h3;
                      }

                      source = source_a + source_b + source_c;

                      source *= wt * pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
                    }

                    lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += mass + advection + diffusion + source;
                  }
                }
              }

              /*
               * J_S_G
               */
              for (p = 0; p < VIM; p++) {
                for (q = 0; q < VIM; q++) {
                  var = v_g[p][q];

                  if (pd->v[pg->imtrx][var]) {
                    pvar = upd->vp[pg->imtrx][var];
                    for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                      phi_j = bf[var]->phi[j];
                      mass = 0.;
                      if (pd->TimeIntegration != STEADY) {
                        if (pd->e[pg->imtrx][eqn] & T_MASS) {
                          mass = mup * (1. + 2. * tt) * phi_j / dt *
                                 ((double)delta(a, p) * (double)delta(b, q) +
                                  (double)delta(b, p) * (double)delta(a, q));
                          mass *= h3 * det_J;

                          mass *= wt_func * lambda * wt * pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
                        }
                      }

                      advection = 0.;
                      advection_a = 0.;

                      if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
                        if (DOUBLE_NONZERO(lambda)) {

                          for (r = 0; r < WIM; r++) {
                            advection_a += mup * (v[r] - x_dot[r]) * bf[var]->grad_phi[j][r];
                          }

                          advection += advection_a * ((double)delta(a, p) * (double)delta(b, q) +
                                                      (double)delta(b, p) * (double)delta(a, q));

                          advection -=
                              mup * phi_j *
                              (2. * (g[p][b] * (double)delta(a, q) +
                                     gt[a][p] * (double)delta(b, q)) +
                               (gt[p][b] * (double)delta(a, q) + gt[a][q] * (double)delta(b, p)) +
                               (g[q][b] * (double)delta(a, p) + g[a][p] * (double)delta(b, q)));

                          advection -= phi_j * (s[p][b] * (double)delta(a, q) +
                                                s[a][p] * (double)delta(b, q));

                          advection *= h3 * det_J;

                          advection *=
                              wt_func * wt * lambda * pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
                        }
                      }

                      /*
                       * Diffusion...
                       */

                      diffusion = 0.;

                      if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {
                        diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
                      }

                      /*
                       * Source term...
                       */

                      source = 0.;

                      if (pd->e[pg->imtrx][eqn] & T_SOURCE) {

                        if (DOUBLE_NONZERO(alpha)) {
                          source +=
                              phi_j *
                                  (s[a][p] * (double)delta(b, q) + s[a][q] * (double)delta(b, p) +
                                   s[q][b] * (double)delta(a, p) + s[p][b] * (double)delta(a, q)) +
                              mup * phi_j *
                                  (g[q][b] * (double)delta(a, p) + g[a][p] * (double)delta(b, q) +
                                   gt[p][b] * (double)delta(a, q) + gt[a][q] * (double)delta(b, p) +
                                   g[p][b] * (double)delta(a, q) + gt[a][p] * (double)delta(b, q) +
                                   g[a][q] * (double)delta(b, p) + gt[q][b] * (double)delta(a, p));
                          source *= alpha * lambda;
                        }

                        source *=
                            det_J * h3 * wt_func * wt * pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
                      }

                      lec->J[LEC_J_INDEX(peqn, pvar, i, j)] +=
                          mass + advection + diffusion + source;
                    }
                  }
                }
              }

              /*
               * J_S_S
               */
              for (p = 0; p < VIM; p++) {
                for (q = 0; q < VIM; q++) {
                  var = v_s[mode][p][q];

                  if (pd->v[pg->imtrx][var]) {
                    pvar = upd->vp[pg->imtrx][var];
                    for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                      phi_j = bf[var]->phi[j];
                      mass = 0.;
                      if (pd->TimeIntegration != STEADY) {
                        if (pd->e[pg->imtrx][eqn] & T_MASS) {
                          mass = (1. + 2. * tt) * phi_j / dt * (double)delta(a, p) *
                                 (double)delta(b, q);
                          mass *= h3 * det_J;
                          mass *= wt_func * lambda * wt * pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
                        }
                      }

                      advection = 0.;

                      if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
                        if (DOUBLE_NONZERO(lambda)) {
                          if ((a == p) && (b == q)) {
                            for (r = 0; r < WIM; r++) {
                              advection += (v[r] - x_dot[r]) * bf[var]->grad_phi[j][r];
                            }
                          }

                          advection -= phi_j * (gt[a][p] * (double)delta(b, q) +
                                                g[q][b] * (double)delta(a, p));

                          advection *= h3 * det_J;

                          advection *=
                              wt_func * wt * lambda * pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
                        }
                      }

                      /*
                       * Diffusion...
                       */

                      diffusion = 0.;

                      if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {
                        diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
                      }

                      /*
                       * Source term...
                       */

                      source = 0.;

                      if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
                        source_a = phi_j * (double)delta(a, p) * (double)delta(b, q);

                        source_b = 0.;
                        if (DOUBLE_NONZERO(alpha)) {
                          source_b =
                              phi_j *
                                  (g[q][b] * (double)delta(a, p) + gt[q][b] * (double)delta(a, p) +
                                   g[a][p] * (double)delta(b, q) + gt[a][p] * (double)delta(b, q)) +
                              phi_j *
                                  (s[q][b] * (double)delta(a, p) + s[a][p] * (double)delta(b, q)) /
                                  mup;

                          source_b *= alpha * lambda;
                        }

                        source = source_a + source_b;

                        source *=
                            det_J * h3 * wt_func * wt * pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
                      }

                      lec->J[LEC_J_INDEX(peqn, pvar, i, j)] +=
                          mass + advection + diffusion + source;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  return (status);
}

/* this stress routine does the EVSS formulation according to Fortin, 1995
 *who uses the regular stress equation and converts stress in the momentum
 *equation by adding the divergence of (g + gT).
 */

int assemble_stress_fortin(dbl tt, /* parameter to vary time integration from
                                    * explicit (tt = 1) to implicit (tt = 0) */
                           dbl dt, /* current time step size */
                           PG_DATA *pg_data) {
  int dim, p, q, r, a, b, w, k;

  int eqn, var;
  int peqn, pvar;
  int evss_gradv = 0;

  int i, j, status, mode;
  dbl v[DIM];      /* Velocity field. */
  dbl x_dot[DIM];  /* current position field derivative wrt time. */
  dbl h3;          /* Volume element (scale factors). */
  dbl dh3dmesh_pj; /* Sensitivity to (p,j) mesh dof. */

  dbl grad_v[DIM][DIM];
  dbl gamma[DIM][DIM]; /* Shear-rate tensor based on velocity */
  dbl det_J;           /* determinant of element Jacobian */

  dbl d_det_J_dmesh_pj; /* for specific (p,j) mesh dof */

  dbl mass; /* For terms and their derivatives */
  dbl mass_a, mass_b;
  dbl advection;
  dbl advection_a, advection_b, advection_c, advection_d;
  dbl diffusion;
  dbl source;
  dbl source1;
  dbl source_a = 0, source_b = 0, source_c = 0;
  int err;
  dbl alpha = 0;      /* This is the Geisekus mobility parameter */
  dbl lambda1 = 0;    /* polymer relaxation constant */
  dbl lambda2 = 0;    /* 2nd polymer relaxation constant -- for modified Jeffreys model */
  dbl elasticMod = 0; /* elastic modulus -- needed for the modified Jeffreys model */
  dbl lambda =
      0; /* lambda1 + lambda2 -- this is just lambda1 unless using the modified Jeffreys model */
  double xi;
  double d_xi_dF[MDE];
  dbl ucwt, lcwt; /* Upper convected derviative weight, Lower convected derivative weight */
  dbl eps = 0;    /* This is the PTT elongation parameter */
  double d_eps_dF[MDE];
  /*
   *
   * Note how carefully we avoid refering to d(phi[i])/dx[j] and refer instead
   * to the j-th component of grad_phi[j][i] so that this vector can be loaded
   * up with components that may be different in non Cartesian coordinate
   * systems.
   *
   * We will, however, insist on *orthogonal* coordinate systems, even if we
   * might permit them to be curvilinear.
   *
   * Assume all components of velocity are interpolated with the same kind
   * of basis function.
   */

  /*
   * Petrov-Galerkin weighting functions for i-th and ab-th stress residuals
   * and some of their derivatives...
   */

  dbl wt_func;

  /*
   * Interpolation functions for variables and some of their derivatives.
   */

  dbl phi_j;
  dbl wt;

  /* Variables for stress */

  int R_s[MAX_MODES][DIM][DIM];
  int v_s[MAX_MODES][DIM][DIM];
  int v_g[DIM][DIM];

  dbl s[DIM][DIM];     /* stress tensor */
  dbl s_dot[DIM][DIM]; /* stress tensor from last time step */
  dbl grad_s[DIM][DIM][DIM];
  dbl d_grad_s_dmesh[DIM][DIM][DIM][DIM]
                    [MDE]; /* derivative of grad of stress tensor for mode ve_mode */

  dbl g[DIM][DIM];      /* velocity gradient tensor */
  dbl gt[DIM][DIM];     /* transpose of velocity gradient tensor */
  dbl g_dot[DIM][DIM];  /* velocity gradient tensor time derivative */
  dbl gt_dot[DIM][DIM]; /* transpose of velocity gradient tensor time derivative */

  /* dot product tensors */

  dbl s_dot_s[DIM][DIM];
  dbl s_dot_g[DIM][DIM];
  dbl g_dot_s[DIM][DIM];
  dbl s_dot_gt[DIM][DIM];
  dbl gt_dot_s[DIM][DIM];

  dbl g_dot_g[DIM][DIM];
  dbl gt_dot_g[DIM][DIM];
  dbl gt_dot_gt[DIM][DIM];

  /* polymer viscosity and derivatives */
  dbl mup;
  VISCOSITY_DEPENDENCE_STRUCT d_mup_struct;
  VISCOSITY_DEPENDENCE_STRUCT *d_mup = &d_mup_struct;

  POLYMER_TIME_CONST_DEPENDENCE_STRUCT d_lam_struct;
  POLYMER_TIME_CONST_DEPENDENCE_STRUCT *d_lam = &d_lam_struct;

  SARAMITO_DEPENDENCE_STRUCT d_saramito_struct;
  SARAMITO_DEPENDENCE_STRUCT *d_saramito = &d_saramito_struct;

  /* 2nd polymer viscosity -- for modified Jefreys model*/
  dbl mupJeff;

  const bool jeffreysEnabled = (vn->ConstitutiveEquation == MODIFIED_JEFFREYS);

  // todo: will want to parse necessary parameters... for now hard code
  const bool saramitoEnabled =
      (vn->ConstitutiveEquation == SARAMITO_OLDROYDB || vn->ConstitutiveEquation == SARAMITO_PTT ||
       vn->ConstitutiveEquation == SARAMITO_GIESEKUS);

  dbl saramitoCoeff = 1.;

  dbl d_mup_dv_pj;
  dbl d_mup_dmesh_pj;

  /*  shift function */
  dbl at = 0.0;
  dbl d_at_dT[MDE];
  dbl wlf_denom;

  /* constitutive equation parameters */
  dbl Z = 1.0; /* This is the factor appearing in front of the stress tensor in PTT */
  dbl dZ_dtrace = 0.0;

  /* advective terms are precalculated */
  dbl v_dot_del_s[DIM][DIM];
  dbl x_dot_del_s[DIM][DIM];

  dbl v_dot_del_g[DIM][DIM];
  dbl v_dot_del_gt[DIM][DIM];

  dbl x_dot_del_g[DIM][DIM];
  dbl x_dot_del_gt[DIM][DIM];

  dbl d_xdotdels_dm;

  dbl d_vdotdels_dm;

  dbl trace = 0.0; /* trace of the stress tensor */

  /* SUPG variables */
  dbl supg = 0;

  if (vn->evssModel == EVSS_GRADV) {
    evss_gradv = 1;
  }

  status = 0;

  eqn = R_STRESS11;

  /*
   * Bail out fast if there's nothing to do...
   */

  if (!pd->e[pg->imtrx][eqn]) {
    return (status);
  }

  /*
   * Unpack variables from structures for local convenience...
   */

  dim = pd->Num_Dim;

  wt = fv->wt;

  det_J = bf[eqn]->detJ; /* Really, ought to be mesh eqn. */

  h3 = fv->h3; /* Differential volume element (scales). */

  /* load eqn and variable number in tensor form */
  (void)stress_eqn_pointer(v_s);
  (void)stress_eqn_pointer(R_s);

  v_g[0][0] = VELOCITY_GRADIENT11;
  v_g[0][1] = VELOCITY_GRADIENT12;
  v_g[1][0] = VELOCITY_GRADIENT21;
  v_g[1][1] = VELOCITY_GRADIENT22;
  v_g[0][2] = VELOCITY_GRADIENT13;
  v_g[1][2] = VELOCITY_GRADIENT23;
  v_g[2][0] = VELOCITY_GRADIENT31;
  v_g[2][1] = VELOCITY_GRADIENT32;
  v_g[2][2] = VELOCITY_GRADIENT33;

  /*
   * Field variables...
   */
  for (a = 0; a < WIM; a++) {
    v[a] = fv->v[a];

    /* note, these are zero for steady calculations */
    x_dot[a] = 0.0;
    if (pd->TimeIntegration != STEADY && pd->gv[MESH_DISPLACEMENT1 + a]) {
      x_dot[a] = fv_dot->x[a];
    }
  }

  /*
   * In Cartesian coordinates, this velocity gradient tensor will
   * have components that are...
   *
   * 			grad_v[a][b] = d v_b
   *				       -----
   *				       d x_a
   */

  for (a = 0; a < VIM; a++) {
    for (b = 0; b < VIM; b++) {
      grad_v[a][b] = fv->grad_v[a][b];
    }
  }

  /* load up shearrate tensor based on velocity */
  for (a = 0; a < VIM; a++) {
    for (b = 0; b < VIM; b++) {
      gamma[a][b] = grad_v[a][b] + grad_v[b][a];
    }
  }

  for (a = 0; a < VIM; a++) {
    for (b = 0; b < VIM; b++) {
      if (evss_gradv) {
        g[a][b] = fv->grad_v[a][b];
        gt[a][b] = fv->grad_v[b][a];
      } else {
        g[a][b] = fv->G[a][b];
        gt[b][a] = g[a][b];
      }
    }
  }

  if (vn->wt_funcModel == GALERKIN) {
    supg = 0.;
  } else if (vn->wt_funcModel == SUPG) {
    supg = vn->wt_func;
  }

  SUPG_terms supg_terms;
  if (supg != 0.) {
    supg_tau(&supg_terms, dim, 0.0, pg_data, dt, TRUE, eqn);
  }
  /* end Petrov-Galerkin addition */

  /*  shift factor  */
  if (pd->gv[TEMPERATURE]) {
    if (vn->shiftModel == CONSTANT) {
      at = vn->shift[0];
      for (j = 0; j < ei[pg->imtrx]->dof[TEMPERATURE]; j++) {
        d_at_dT[j] = 0.;
      }
    } else if (vn->shiftModel == MODIFIED_WLF) {
      wlf_denom = vn->shift[1] + fv->T - mp->reference[TEMPERATURE];
      if (wlf_denom != 0.) {
        at = exp(vn->shift[0] * (mp->reference[TEMPERATURE] - fv->T) / wlf_denom);
        for (j = 0; j < ei[pg->imtrx]->dof[TEMPERATURE]; j++) {
          d_at_dT[j] =
              -at * vn->shift[0] * vn->shift[1] / (wlf_denom * wlf_denom) * bf[TEMPERATURE]->phi[j];
        }
      } else {
        at = 1.;
      }
      for (j = 0; j < ei[pg->imtrx]->dof[TEMPERATURE]; j++) {
        d_at_dT[j] = 0.;
      }
    }
  } else {
    at = 1.;
  }

  // if a modified Jeffreys model is being run, load time derivative of velocity gradient
  if (jeffreysEnabled) {
    (void)tensor_dot(g, g, g_dot_g, VIM);
    (void)tensor_dot(gt, gt, gt_dot_gt, VIM);
    (void)tensor_dot(gt, g, gt_dot_g, VIM);

    for (a = 0; a < VIM; a++) {
      for (b = 0; b < VIM; b++) {
        v_dot_del_g[a][b] = 0.;
        v_dot_del_gt[a][b] = 0.;
        x_dot_del_g[a][b] = 0.;
        x_dot_del_gt[a][b] = 0.;

        for (q = 0; q < dim; q++) {
          v_dot_del_g[a][b] += v[q] * fv->grad_G[q][a][b];
          v_dot_del_gt[a][b] += v[q] * fv->grad_G[q][b][a];
          x_dot_del_g[a][b] += x_dot[q] * fv->grad_G[q][a][b];
          x_dot_del_gt[a][b] += x_dot[q] * fv->grad_G[q][b][a];
        }

        if (pd->TimeIntegration != STEADY) {
          g_dot[a][b] = fv_dot->G[a][b];
        } else {
          g_dot[a][b] = 0.;
        }
        gt_dot[b][a] = g_dot[a][b];
      }
    }
  }

  /* Begin loop over modes */
  for (mode = 0; mode < vn->modes; mode++) {

    load_modal_pointers(mode, tt, dt, s, s_dot, grad_s, d_grad_s_dmesh);

    /* precalculate advective terms of form (v dot del tensor)*/

    trace = 0.0;

    for (a = 0; a < VIM; a++) {
      trace += s[a][a];
      for (b = 0; b < VIM; b++) {
        v_dot_del_s[a][b] = 0.;
        x_dot_del_s[a][b] = 0.;
        for (q = 0; q < WIM; q++) {
          v_dot_del_s[a][b] += v[q] * grad_s[q][a][b];
          x_dot_del_s[a][b] += x_dot[q] * grad_s[q][a][b];
        }
      }
    }

    /*
     * Stress tensor...(Note "anti-BSL" sign convention on deviatoric stress)
     */

    /* get polymer viscosity */
    mup = viscosity(ve[mode]->gn, gamma, d_mup);

    if (saramitoEnabled == TRUE) {
      compute_saramito_model_terms(&saramitoCoeff, d_saramito, s, ve[mode]->gn, FALSE);
    } else {
      saramitoCoeff = 1.;
      d_saramito->tau_y = 0;

      for (int i = 0; i < VIM; ++i) {
        for (int j = 0; j < VIM; ++j) {
          d_saramito->s[i][j] = 0;
        }
      }
    }

    double d_alpha_dF[MDE];
    /* get Geisekus mobility parameter */
    if (ve[mode]->alphaModel == CONSTANT) {
      alpha = ve[mode]->alpha;
    } else if (ls != NULL && ve[mode]->alphaModel == VE_LEVEL_SET) {
      double pos_alpha = ve[mode]->pos_ls.alpha;
      double neg_alpha = ve[mode]->alpha;
      double width = ls->Length_Scale;
      err = level_set_property(neg_alpha, pos_alpha, width, &alpha, d_alpha_dF);
      GOMA_EH(err, "level_set_property() failed for mobility parameter.");
    } else {
      GOMA_EH(GOMA_ERROR, "Unknown mobility parameter model");
    }

    /* get time constant */
    lambda = polymer_time_const(ve[mode]->time_const_st, gamma, d_lam);

    xi = 0;
    if (ve[mode]->xiModel == CONSTANT) {
      xi = ve[mode]->xi;
    } else if (ls != NULL && ve[mode]->xiModel == VE_LEVEL_SET) {
      double pos_xi = ve[mode]->pos_ls.xi;
      double neg_xi = ve[mode]->xi;
      double width = ls->Length_Scale;
      err = level_set_property(neg_xi, pos_xi, width, &xi, d_xi_dF);
      GOMA_EH(err, "level_set_property() failed for ptt xi parameter.");
    } else {
      GOMA_EH(GOMA_ERROR, "Unknown PTT Xi parameter model");
    }

    ucwt = 1.0 - xi / 2.0;
    lcwt = xi / 2.0;

    if (ve[mode]->epsModel == CONSTANT) {
      eps = ve[mode]->eps;
    } else if (ls != NULL && ve[mode]->epsModel == VE_LEVEL_SET) {
      double pos_eps = ve[mode]->pos_ls.eps;
      double neg_eps = ve[mode]->eps;
      double width = ls->Length_Scale;
      err = level_set_property(neg_eps, pos_eps, width, &eps, d_eps_dF);
      GOMA_EH(err, "level_set_property() failed for ptt epsilon parameter.");
    } else {
      GOMA_EH(GOMA_ERROR, "Unknown PTT Epsilon parameter model");
    }

    if (jeffreysEnabled) {
      mupJeff = ve[mode]->muJeffreys;
      // if the modified Jeffreys model is used, the parsed value of lambda is the
      // elastic modulus rather than the time consant
      elasticMod = lambda;
      lambda1 = mup / elasticMod; // mup/G
      lambda2 = mupJeff / elasticMod;
      lambda = lambda1 + lambda2;
    }

    Z = 1.0;
    dZ_dtrace = 0;
    if (vn->ConstitutiveEquation == PTT) {
      if (vn->ptt_type == PTT_LINEAR) {
        Z = 1 + eps * lambda * trace / mup;
        dZ_dtrace = eps * lambda / mup;
      } else if (vn->ptt_type == PTT_EXPONENTIAL) {
        Z = exp(eps * lambda * trace / mup);
        dZ_dtrace = Z * eps * lambda / mup;
      } else {
        GOMA_EH(GOMA_ERROR, "Unrecognized PTT Form %d", vn->ptt_type);
      }
    }

    /* get tensor dot products for future use */

    if (DOUBLE_NONZERO(alpha))
      (void)tensor_dot(s, s, s_dot_s, VIM);

    if (ucwt != 0.) {
      (void)tensor_dot(s, g, s_dot_g, VIM);
      (void)tensor_dot(gt, s, gt_dot_s, VIM);
    }

    if (lcwt != 0.) {
      (void)tensor_dot(s, gt, s_dot_gt, VIM);
      (void)tensor_dot(g, s, g_dot_s, VIM);
    }
    /*
     * Residuals_________________________________________________________________
     */

    if (af->Assemble_Residual) {
      /*
       * Assemble each component "ab" of the polymer stress equation...
       */
      for (a = 0; a < VIM; a++) {
        for (b = 0; b < VIM; b++) {

          if (a <= b) /* since the stress tensor is symmetric, only assemble the upper half */
          {
            eqn = R_s[mode][a][b];

            /*
             * In the element, there will be contributions to this many equations
             * based on the number of degrees of freedom...
             */

            for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
              wt_func = bf[eqn]->phi[i];
              /* add Petrov-Galerkin terms as necessary */
              if (supg != 0.) {
                for (w = 0; w < dim; w++) {
                  wt_func += supg * supg_terms.supg_tau * v[w] * bf[eqn]->grad_phi[i][w];
                }
              }

              mass = 0.;

              if (pd->TimeIntegration != STEADY) {
                if (pd->e[pg->imtrx][eqn] & T_MASS) {
                  mass = s_dot[a][b];
                  mass *= wt_func * at * lambda * det_J * wt;
                  mass *= h3;
                  mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
                }
              }

              advection = 0.;
              if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
                if (DOUBLE_NONZERO(lambda)) {

                  advection += v_dot_del_s[a][b] - x_dot_del_s[a][b];
                  if (ucwt != 0.)
                    advection -= ucwt * (gt_dot_s[a][b] + s_dot_g[a][b]);
                  if (lcwt != 0.)
                    advection += lcwt * (s_dot_gt[a][b] + g_dot_s[a][b]);

                  advection *= wt_func * at * lambda * det_J * wt * h3;
                  advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
                }
              }

              diffusion = 0.;
              if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {
                diffusion *= det_J * wt * h3;
                diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
              }

              /*
               * Source term...
               */

              source = 0.;
              if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
                // consider whether saramitoCoeff should multiply here
                source += saramitoCoeff * Z * s[a][b] - at * mup * (g[a][b] + gt[a][b]);

                if (DOUBLE_NONZERO(alpha)) {
                  source1 = (s_dot_s[a][b] / mup);

                  source1 *= alpha * lambda * saramitoCoeff;
                  source += source1;
                }

                if (jeffreysEnabled) {
                  source -= mup * lambda2 *
                            (g_dot[a][b] + gt_dot[a][b] + v_dot_del_g[a][b] + v_dot_del_gt[a][b] -
                             (g_dot_g[a][b] + 2 * gt_dot_g[a][b] + gt_dot_gt[a][b]));
                }

                source *= wt_func * det_J * h3 * wt;

                source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
              }

              /*
               * Add contributions to this residual (globally into Resid, and
               * locally into an accumulator)
               */

              lec->R[LEC_R_INDEX(upd->ep[pg->imtrx][eqn], i)] +=
                  mass + advection + diffusion + source;
            }
          }
        }
      }
    }

    /*
     * Jacobian terms...
     */

    if (af->Assemble_Jacobian) {
      dbl R_source, R_advection; /* Places to put the raw residual portions
                                    instead of constantly recalcing them */
      for (a = 0; a < VIM; a++) {
        for (b = 0; b < VIM; b++) {
          if (a <= b) /* since the stress tensor is symmetric, only assemble the upper half */
          {
            eqn = R_s[mode][a][b];
            peqn = upd->ep[pg->imtrx][eqn];

            R_advection = v_dot_del_s[a][b] - x_dot_del_s[a][b];
            if (ucwt != 0.)
              R_advection -= ucwt * (gt_dot_s[a][b] + s_dot_g[a][b]);
            if (lcwt != 0.)
              R_advection += lcwt * (s_dot_gt[a][b] + g_dot_s[a][b]);

            R_source = Z * s[a][b];

            if (DOUBLE_NONZERO(alpha))
              R_source += alpha * lambda * (s_dot_s[a][b] / mup);
            R_source *= saramitoCoeff;
            R_source += -at * mup * (g[a][b] + gt[a][b]);

            if (jeffreysEnabled) {
              R_source -= at * lambda2 * mup *
                          (g_dot[a][b] + gt_dot[a][b] + v_dot_del_g[a][b] - x_dot_del_g[a][b] +
                           v_dot_del_gt[a][b] - x_dot_del_gt[b][a] -
                           (g_dot_g[a][b] + 2 * gt_dot_g[a][b] + gt_dot_gt[a][b]));
            }

            for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {

              wt_func = bf[eqn]->phi[i];
              /* add Petrov-Galerkin terms as necessary */
              if (supg != 0.) {
                for (w = 0; w < dim; w++) {
                  wt_func += supg * supg_terms.supg_tau * v[w] * bf[eqn]->grad_phi[i][w];
                }
              }

              /*
               * Set up some preliminaries that are needed for the (a,i)
               * equation for bunches of (b,j) column variables...
               */

              /*
               * J_S_T
               */

              var = TEMPERATURE;
              if (pd->v[pg->imtrx][var]) {
                pvar = upd->vp[pg->imtrx][var];
                for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                  phi_j = bf[var]->phi[j];

                  mass = 0.;
                  dbl d_lambda_dT = d_lam->T[j];

                  if (pd->TimeIntegration != STEADY) {
                    if (pd->e[pg->imtrx][eqn] & T_MASS) {
                      mass = s_dot[a][b] * d_at_dT[j] * lambda;

                      if (jeffreysEnabled)
                        mass += s_dot[a][b] * at * d_lambda_dT;

                      mass *= wt_func * det_J * wt;
                      mass *= h3;
                      mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
                    }
                  }

                  advection = 0.;
                  if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
                    if (DOUBLE_NONZERO(lambda)) {

                      advection += v_dot_del_s[a][b] - x_dot_del_s[a][b];
                      if (ucwt != 0.)
                        advection -= ucwt * (gt_dot_s[a][b] + s_dot_g[a][b]);
                      if (lcwt != 0.)
                        advection += lcwt * (s_dot_gt[a][b] + g_dot_s[a][b]);

                      advection *= d_at_dT[j] * lambda + (jeffreysEnabled ? d_lambda_dT * at : 0);
                      advection *=
                          pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)] * wt_func * det_J * wt * h3;
                    }
                  }

                  source = 0.;
                  source1 = 0.;
                  if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
                    source = -(g[a][b] + gt[a][b]) * (at * d_mup->T[j] + mup * d_at_dT[j]);

                    if (jeffreysEnabled) {
                      source -= d_mup->T[j] * lambda2 *
                                (g_dot[a][b] + gt_dot[a][b] + v_dot_del_g[a][b] -
                                 x_dot_del_g[a][b] + v_dot_del_gt[a][b] - v_dot_del_gt[b][a] -
                                 (g_dot_g[a][b] + 2 * gt_dot_g[a][b] + gt_dot_gt[a][b]));
                    }
                    if (DOUBLE_NONZERO(alpha)) {
                      source1 -= s_dot_s[a][b] / (mup * mup) * d_mup->T[j];
                      source1 *= lambda * alpha * saramitoCoeff;
                      source += source1;
                    }
                    source *= wt_func * det_J * wt * h3;
                    source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
                  }

                  lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += mass + advection + source;
                }
              }

              /*
               * J_S_v
               */
              for (p = 0; p < WIM; p++) {
                var = VELOCITY1 + p;
                if (pd->v[pg->imtrx][var]) {
                  pvar = upd->vp[pg->imtrx][var];
                  for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                    phi_j = bf[var]->phi[j];
                    d_mup_dv_pj = d_mup->v[p][j];
                    dbl d_lambda_dv_pj = d_lam->v[p][j];
                    if (jeffreysEnabled) {
                      d_lambda_dv_pj = d_mup_dv_pj / elasticMod;
                    }

                    mass = 0.;

                    if (pd->TimeIntegration != STEADY) {
                      if (pd->e[pg->imtrx][eqn] & T_MASS) {
                        if (supg != 0.) {
                          mass = supg * supg_terms.supg_tau * phi_j * bf[eqn]->grad_phi[i][p];

                          for (w = 0; w < dim; w++) {
                            mass += supg * supg_terms.d_supg_tau_dv[p][j] * v[w] *
                                    bf[eqn]->grad_phi[i][w];
                          }

                          mass *= s_dot[a][b];
                        }

                        mass *=
                            pd->etm[pg->imtrx][eqn][(LOG2_MASS)] * at * lambda * det_J * wt * h3;
                        mass += s_dot[a][b] * wt_func * pd->etm[pg->imtrx][eqn][(LOG2_MASS)] * at *
                                d_lambda_dv_pj * det_J * wt * h3;
                      }
                    }

                    advection = 0.;

                    if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
                      if (DOUBLE_NONZERO(lambda)) {
                        advection_a = phi_j * (grad_s[p][a][b]);

                        advection_a *= wt_func;

                        advection_b = 0.;
                        /* Petrov-Galerkin term */
                        if (supg != 0.) {

                          advection_b =
                              supg * supg_terms.supg_tau * phi_j * bf[eqn]->grad_phi[i][p];
                          for (w = 0; w < dim; w++) {
                            advection_b += supg * supg_terms.d_supg_tau_dv[p][j] * v[w] *
                                           bf[eqn]->grad_phi[i][w];
                          }

                          advection_b *= R_advection;
                        }

                        advection_c = 0.;
                        if (evss_gradv) {
                          if (pd->CoordinateSystem != CYLINDRICAL) {
                            if (ucwt != 0) {
                              for (k = 0; k < VIM; k++) {
                                advection_c -=
                                    ucwt * (bf[VELOCITY1 + a]->grad_phi_e[j][p][k][a] * s[k][b] +
                                            bf[VELOCITY1 + b]->grad_phi_e[j][p][k][b] * s[a][k]);
                              }
                            }
                            if (lcwt != 0.) {
                              for (k = 0; k < VIM; k++) {
                                advection_c +=
                                    lcwt * (bf[VELOCITY1 + b]->grad_phi_e[j][p][b][k] * s[a][k] +
                                            bf[VELOCITY1 + a]->grad_phi_e[j][p][a][k] * s[k][b]);
                              }
                            }
                          } else {
                            if (ucwt != 0) {
                              for (k = 0; k < VIM; k++) {
                                advection_c -=
                                    ucwt * (bf[VELOCITY1]->grad_phi_e[j][p][k][a] * s[k][b] +
                                            bf[VELOCITY1]->grad_phi_e[j][p][k][b] * s[a][k]);
                              }
                            }
                            if (lcwt != 0.) {
                              for (k = 0; k < VIM; k++) {
                                advection_c +=
                                    lcwt * (bf[VELOCITY1]->grad_phi_e[j][p][b][k] * s[a][k] +
                                            bf[VELOCITY1]->grad_phi_e[j][p][a][k] * s[k][b]);
                              }
                            }
                          }
                          advection_c *= wt_func;
                        }

                        advection = advection_a + advection_b + advection_c;
                        advection *= at * lambda * det_J * wt * h3;
                        advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];

                        advection += R_advection * d_lambda_dv_pj * wt_func * det_J * wt * h3 *
                                     pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
                      }
                    }

                    diffusion = 0.;
                    if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {
                      diffusion *= det_J * wt * h3;
                      diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
                    }

                    source = 0.;

                    if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
                      source_c = -at * d_mup_dv_pj * (g[a][b] + gt[a][b]);

                      if (jeffreysEnabled) {
                        source_c -= d_mup_dv_pj * at * lambda2 *
                                    (g_dot[a][b] + gt_dot[a][b] + v_dot_del_g[a][b] -
                                     x_dot_del_g[a][b] + v_dot_del_gt[a][b] - x_dot_del_gt[b][a] -
                                     (g_dot_g[a][b] + 2 * gt_dot_g[a][b] + gt_dot_gt[a][b]));

                        source_c -= lambda2 * mup *
                                    (phi_j * fv->grad_G[p][a][b] + phi_j * fv->grad_G[p][b][a]);
                      }
                      if (evss_gradv) {
                        if (pd->CoordinateSystem != CYLINDRICAL) {
                          source_c -= at * mup *
                                      (bf[VELOCITY1 + a]->grad_phi_e[j][p][a][b] +
                                       bf[VELOCITY1 + b]->grad_phi_e[j][p][b][a]);
                        } else {
                          source_c -= at * mup *
                                      (bf[VELOCITY1]->grad_phi_e[j][p][a][b] +
                                       bf[VELOCITY1]->grad_phi_e[j][p][b][a]);
                        }
                      }
                      source_c *= wt_func;

                      source_a = 0.;
                      if (DOUBLE_NONZERO(alpha)) {
                        source_a = -s_dot_s[a][b] / (mup * mup);
                        source_a *= wt_func * saramitoCoeff * alpha *
                                    (d_lambda_dv_pj * mup + lambda * d_mup_dv_pj);
                      }

                      source_b = 0.;
                      if (supg != 0.) {
                        source_b = supg * supg_terms.supg_tau * phi_j * bf[eqn]->grad_phi[i][p];

                        for (w = 0; w < dim; w++) {
                          source_b += supg * supg_terms.d_supg_tau_dv[p][j] * v[w] *
                                      bf[eqn]->grad_phi[i][w];
                        }

                        source_b *= R_source;
                      }

                      source = source_a + source_b + source_c;
                      source *= det_J * wt * h3;
                      source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
                    }

                    lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += mass + advection + diffusion + source;
                  }
                }
              }

              /*
               * J_S_c
               */
              var = MASS_FRACTION;
              if (pd->v[pg->imtrx][var]) {
                for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                  phi_j = bf[var]->phi[j];

                  for (w = 0; w < pd->Num_Species_Eqn; w++) {

                    source = 0.;

                    if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
                      source_a = -at * d_mup->C[w][j] * (g[a][b] + gt[a][b]);

                      if (jeffreysEnabled) {
                        source_a -= at * d_mup->C[w][j] * lambda2 *
                                    (g_dot[a][b] + gt_dot[a][b] + v_dot_del_g[a][b] -
                                     x_dot_del_g[a][b] + v_dot_del_gt[a][b] - v_dot_del_gt[b][a] -
                                     (g_dot_g[a][b] + 2 * gt_dot_g[a][b] + gt_dot_gt[a][b]));
                      }

                      source_b = 0.;
                      if (DOUBLE_NONZERO(alpha)) {
                        source_b -= s_dot_s[a][b] / (mup * mup);
                        source_b *= alpha * lambda * saramitoCoeff * d_mup->C[w][j];
                      }
                      source = source_a + source_b;
                      source *= wt_func * det_J * wt * h3;
                      source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
                    }

                    if (w > 1) {
                      GOMA_EH(GOMA_ERROR, "Need more arrays for each species.");
                    }

                    lec->J[LEC_J_INDEX(peqn, MAX_PROB_VAR + w, i, j)] += source;
                  }
                }
              }

              /*
               * J_S_P
               */
              var = PRESSURE;
              if (pd->v[pg->imtrx][var]) {
                pvar = upd->vp[pg->imtrx][var];
                for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                  phi_j = bf[var]->phi[j];

                  source = 0.;
                  if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
                    source_a += -at * d_mup->P[j] * (g[a][b] + gt[a][b]);

                    if (jeffreysEnabled) {
                      source_a -= at * d_mup->P[j] * lambda2 *
                                  (g_dot[a][b] + gt_dot[a][b] + v_dot_del_g[a][b] -
                                   x_dot_del_g[a][b] + v_dot_del_gt[a][b] - v_dot_del_gt[b][a] -
                                   (g_dot_g[a][b] + 2 * gt_dot_g[a][b] + gt_dot_gt[a][b]));
                    }

                    source_b = 0.;
                    if (DOUBLE_NONZERO(alpha)) {
                      source_b -= (s_dot_s[a][b] / (mup * mup));
                      source_b *= d_mup->P[j] * alpha * lambda * saramitoCoeff;
                    }
                    source = source_a + source_b;
                    source *= wt_func * det_J * wt * h3;
                    source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
                  }

                  lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += source;
                }
              }

              /*
               * J_S_d
               */
              for (p = 0; p < dim; p++) {
                var = MESH_DISPLACEMENT1 + p;
                if (pd->v[pg->imtrx][var]) {
                  pvar = upd->vp[pg->imtrx][var];
                  for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                    phi_j = bf[var]->phi[j];
                    d_det_J_dmesh_pj = bf[eqn]->d_det_J_dm[p][j];
                    dh3dmesh_pj = fv->dh3dq[p] * bf[var]->phi[j];
                    d_mup_dmesh_pj = d_mup->X[p][j];

                    mass = 0.;
                    mass_a = 0.;
                    mass_b = 0.;
                    if (pd->TimeIntegration != STEADY) {
                      if (pd->e[pg->imtrx][eqn] & T_MASS) {
                        mass_a = s_dot[a][b];
                        mass_a *= wt_func * (d_det_J_dmesh_pj * h3 + det_J * dh3dmesh_pj);

                        if (supg != 0.) {
                          for (w = 0; w < dim; w++) {
                            mass_b += supg * (supg_terms.supg_tau * v[w] *
                                                  bf[eqn]->d_grad_phi_dmesh[i][w][p][j] +
                                              supg_terms.d_supg_tau_dX[p][j] * v[w] *
                                                  bf[eqn]->grad_phi[i][w]);
                          }
                          mass_b *= s_dot[a][b] * h3 * det_J;
                        }

                        mass = mass_a + mass_b;
                        mass *= at * lambda * wt * pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
                      }
                    }

                    advection = 0.;

                    if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
                      if (DOUBLE_NONZERO(lambda)) {
                        /*
                         * Four parts:
                         *    advection_a =
                         *    	Int ( ea.(v-xdot).d(Vv)/dmesh h3 |Jv| )
                         *
                         *    advection_b =
                         *  (i)	Int ( ea.(v-xdot).Vv h3 d(|Jv|)/dmesh )
                         *  (ii)  Int ( ea.(v-xdot).d(Vv)/dmesh h3 |Jv| )
                         *  (iii) Int ( ea.(v-xdot).Vv dh3/dmesh |Jv|   )
                         *
                         * For unsteady problems, we have an
                         * additional term
                         *
                         *    advection_c =
                         *    	Int ( ea.d(v-xdot)/dmesh.Vv h3 |Jv| )
                         */

                        advection_a = R_advection;

                        advection_a *= wt_func * (d_det_J_dmesh_pj * h3 + det_J * dh3dmesh_pj);

                        d_vdotdels_dm = 0.;
                        for (q = 0; q < WIM; q++) {
                          d_vdotdels_dm += (v[q] - x_dot[q]) * d_grad_s_dmesh[q][a][b][p][j];
                        }

                        advection_b = d_vdotdels_dm;
                        advection_b *= wt_func * det_J * h3;

                        advection_c = 0.;
                        if (pd->TimeIntegration != STEADY) {
                          if (pd->e[pg->imtrx][eqn] & T_MASS) {
                            d_xdotdels_dm = (1. + 2. * tt) * phi_j / dt * grad_s[p][a][b];

                            advection_c -= d_xdotdels_dm;

                            advection_c *= wt_func * h3 * det_J;
                          }
                        }

                        advection_d = 0.;
                        if (supg != 0.) {
                          for (w = 0; w < dim; w++) {
                            advection_d += supg * (supg_terms.supg_tau * v[w] *
                                                       bf[eqn]->d_grad_phi_dmesh[i][w][p][j] +
                                                   supg_terms.d_supg_tau_dX[p][j] * v[w] *
                                                       bf[eqn]->grad_phi[i][w]);
                          }

                          advection_d *= (R_advection)*det_J * h3;
                        }

                        advection = advection_a + advection_b + advection_c + advection_d;

                        advection *= wt * at * lambda * pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
                      }
                    }

                    diffusion = 0.;
                    if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {
                      diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
                    }

                    /*
                     * Source term...
                     */

                    source = 0.;

                    if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
                      source_a = R_source;
                      source_b = -at * (g[a][b] + gt[a][b]);

                      if (DOUBLE_NONZERO(alpha)) {
                        source_b += -s_dot_s[a][b] / (mup * mup) * alpha * lambda * saramitoCoeff;
                      }

                      source_a *= wt_func * (d_det_J_dmesh_pj * h3 + det_J * dh3dmesh_pj);

                      source_b *= wt_func * det_J * h3 * d_mup_dmesh_pj;

                      source_c = 0.;
                      if (supg != 0.) {
                        for (w = 0; w < dim; w++) {
                          source_c +=
                              supg *
                              (supg_terms.supg_tau * v[w] * bf[eqn]->d_grad_phi_dmesh[i][w][p][j] +
                               supg_terms.d_supg_tau_dX[p][j] * v[w] * bf[eqn]->grad_phi[i][w]);
                        }
                        source_c *= R_source * det_J * h3;
                      }

                      dbl source_jeff = 0;
                      if (jeffreysEnabled) {

                        dbl source_ja = 0;
                        dbl source_jb = 0;
                        dbl source_jc = 0;
                        dbl source_jd = 0;

                        source_ja -= mup * lambda2 * at *
                                     (g_dot[a][b] + gt_dot[a][b] + v_dot_del_g[a][b] -
                                      x_dot_del_g[a][b] + v_dot_del_gt[a][b] - v_dot_del_gt[b][a] -
                                      (g_dot_g[a][b] + 2 * gt_dot_g[a][b] + gt_dot_gt[a][b]));

                        source_ja *= wt_func * (d_det_J_dmesh_pj * h3 + det_J * dh3dmesh_pj);

                        dbl d_vdotdelg_dm = 0.;
                        dbl d_vdotdelgt_dm = 0.;

                        for (q = 0; q < dim; q++) {
                          d_vdotdels_dm += (v[q] - x_dot[q]) * d_grad_s_dmesh[q][a][b][p][j];
                          d_vdotdelg_dm += (v[q] - x_dot[q]) * fv->d_grad_G_dmesh[q][a][b][p][j];
                          d_vdotdelgt_dm += (v[q] - x_dot[q]) * fv->d_grad_G_dmesh[q][b][a][p][j];
                        }

                        source_jb += (lambda2 * at * mup * (d_vdotdelg_dm + d_vdotdelgt_dm));
                        source_jb -= d_mup_dmesh_pj * lambda2 *
                                     (g_dot[a][b] + gt_dot[a][b] + v_dot_del_g[a][b] -
                                      x_dot_del_g[a][b] + v_dot_del_gt[a][b] - v_dot_del_gt[b][a] -
                                      (g_dot_g[a][b] + 2 * gt_dot_g[a][b] + gt_dot_gt[a][b]));

                        source_jb *= wt_func * det_J * h3;

                        if (pd->TimeIntegration != STEADY) {
                          const dbl d_xdotdelg_dm =
                              (1. + 2. * tt) * phi_j / dt * fv->grad_G[p][a][b];
                          const dbl d_xdotdelgt_dm =
                              (1. + 2. * tt) * phi_j / dt * fv->grad_G[p][b][a];

                          source_jc -= lambda2 * mup * (d_xdotdelg_dm + d_xdotdelgt_dm);

                          source_jc *= wt_func * h3 * det_J;
                        }

                        if (supg != 0.) {
                          for (w = 0; w < dim; w++) {
                            source_jd += supg * (supg_terms.supg_tau * v[w] *
                                                     bf[eqn]->d_grad_phi_dmesh[i][w][p][j] +
                                                 supg_terms.d_supg_tau_dX[p][j] * v[w] *
                                                     bf[eqn]->grad_phi[i][w]);
                          }
                          source_jd *= R_source * det_J * h3;

                          // advection_d *=
                          //     -lambda2 * at * det_J * h3 *
                          //     (g_dot[a][b] + gt_dot[a][b] + v_dot_del_g[a][b] - x_dot_del_g[a][b]
                          //     +
                          //      v_dot_del_gt[a][b] - v_dot_del_gt[b][a] -
                          //      (g_dot_g[a][b] + 2 * gt_dot_g[a][b] + gt_dot_gt[a][b]));
                        }
                        source_jeff = source_ja + source_jb + source_jc + source_jb + source_jd;
                      }

                      source = source_a + source_b + source_c + source_jeff;

                      source *= wt * pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
                    }

                    lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += mass + advection + diffusion + source;
                  }
                }
              }

              /*
               * J_S_G
               */
              if (evss_gradv == 0) {
                for (p = 0; p < VIM; p++) {
                  for (q = 0; q < VIM; q++) {
                    var = v_g[p][q];

                    if (pd->v[pg->imtrx][var]) {
                      pvar = upd->vp[pg->imtrx][var];
                      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                        phi_j = bf[var]->phi[j];
                        dbl dG[DIM][DIM] = {{0.0}};
                        dbl dGt[DIM][DIM] = {{0.0}};
                        dG[p][q] = phi_j;
                        dGt[q][p] = phi_j;
                        advection = 0.;
                        advection_a = 0.;
                        if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
                          if (DOUBLE_NONZERO(lambda)) {

                            advection -= ucwt * (s[p][b] * (double)delta(a, q) +
                                                 s[a][p] * (double)delta(b, q));
                            advection += lcwt * (s[a][q] * (double)delta(p, b) +
                                                 s[q][b] * (double)delta(a, p));

                            advection *= phi_j * h3 * det_J;

                            advection *= wt_func * wt * at * lambda *
                                         pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
                          }
                        }

                        /*
                         * Diffusion...
                         */

                        diffusion = 0.;

                        if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {
                          diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
                        }

                        /*
                         * Source term...
                         */

                        source = 0.;

                        if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
                          source = -at * mup * phi_j *
                                   ((double)delta(a, p) * (double)delta(b, q) +
                                    (double)delta(b, p) * (double)delta(a, q));
                          if (jeffreysEnabled) {
                            dbl source_jeffrey = 0;
                            // g_dot
                            source_jeffrey += (1. + 2. * tt) * phi_j / dt * (double)delta(a, p) *
                                              (double)delta(b, q);
                            // gt_dot
                            source_jeffrey += (1. + 2. * tt) * phi_j / dt * (double)delta(a, q) *
                                              (double)delta(b, p);
                            //// v_dot_del_g
                            if ((a == p) && (b == q)) {
                              for (r = 0; r < dim; r++) {
                                source_jeffrey += (v[r] - x_dot[r]) * bf[var]->grad_phi[j][r];
                              }
                            }
                            //// v_dot_del_gt
                            if ((a == q) && (b == p)) {
                              for (r = 0; r < dim; r++) {
                                source_jeffrey += (v[r] - x_dot[r]) * bf[var]->grad_phi[j][r];
                              }
                            }
                            //// g_dot_g
                            dbl dG_dot_g[DIM][DIM] = {{0.}};
                            dbl g_dot_dG[DIM][DIM] = {{0.}};
                            tensor_dot(dG, g, dG_dot_g, VIM);
                            tensor_dot(g, dG, g_dot_dG, VIM);
                            source_jeffrey -= g_dot_dG[a][b] + dG_dot_g[a][b];
                            //// g_dot_gt
                            // source_jeffrey -=  (g[a][p] * (double)delta(b, q) + gt[q][b] *
                            // (double)delta(a, p)) *bf[var]->phi[j];
                            //// 2*gt_dot_g
                            dbl dGt_dot_g[DIM][DIM] = {{0.}};
                            dbl gt_dot_dG[DIM][DIM] = {{0.}};
                            tensor_dot(dGt, g, dGt_dot_g, VIM);
                            tensor_dot(gt, dG, gt_dot_dG, VIM);
                            source_jeffrey -= 2. * (dGt_dot_g[a][b] + gt_dot_dG[a][b]);
                            // source_jeffrey -= 2.*(gt[a][p] * (double)delta(b, q) + g[q][b] *
                            // (double)delta(a, p)) * bf[var]->phi[j];
                            //// gt_dot_gt
                            source_jeffrey -= g_dot_dG[b][a] + dG_dot_g[b][a];
                            // source_jeffrey -=  (gt[a][p] * (double)delta(b, q) + g[q][b] *
                            // (double)delta(a, p)) * bf[var]->phi[j];

                            source_jeffrey *= -mup * lambda2;
                            source += source_jeffrey;
                          }
                          source *=
                              det_J * h3 * wt_func * wt * pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
                        }

                        lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += advection + diffusion + source;
                      }
                    }
                  }
                }
              }

              /*
               * J_S_F
               */
              var = FILL;
              if (pd->v[pg->imtrx][var]) {
                pvar = upd->vp[pg->imtrx][var];
                for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                  phi_j = bf[var]->phi[j];

                  mass = 0.;

                  if (pd->TimeIntegration != STEADY) {
                    if (pd->e[pg->imtrx][eqn] & T_MASS) {

                      mass = s_dot[a][b];
                      mass *= d_lam->F[j];
                      mass *= wt_func * at * det_J * wt;
                      mass *= h3;
                      mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
                    }
                  }

                  advection = 0.;

                  if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
                    if (d_lam->F[j] != 0.) {

                      advection += v_dot_del_s[a][b] - x_dot_del_s[a][b];
                      if (ucwt != 0.)
                        advection -= ucwt * (gt_dot_s[a][b] + s_dot_g[a][b]);
                      if (lcwt != 0.)
                        advection += lcwt * (s_dot_gt[a][b] + g_dot_s[a][b]);

                      advection *= d_lam->F[j];
                      advection *= wt_func * at * det_J * wt * h3;
                      advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
                    }
                  }

                  diffusion = 0.;

                  if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {
                    /* add SU term in here when appropriate */

                    diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
                  }

                  source = 0.;

                  if (pd->e[pg->imtrx][eqn] & T_SOURCE) {

                    double invmup = 1 / mup;
                    // PTT
                    if (eps != 0) {
                      // product rule + exponential
                      source += Z *
                                ((lambda * trace * d_eps_dF[j] * invmup) +
                                 (d_lam->F[j] * trace * eps * invmup) -
                                 (lambda * trace * eps * d_mup->F[j] * invmup * invmup)) *
                                s[a][b];
                    }

                    source += -at * d_mup->F[j] * (g[a][b] + gt[a][b]);

                    // Giesekus
                    if (alpha != 0.) {
                      source += s_dot_s[a][b] *
                                (-alpha * lambda * d_mup->F[j] * invmup * invmup +
                                 d_alpha_dF[j] * lambda * invmup + alpha * d_lam->F[j] * invmup);
                    }

                    source *= wt_func * det_J * h3 * wt;

                    source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
                  }

                  lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += mass + advection + diffusion + source;
                }
              }

              /*
               * J_S_S
               */
              for (p = 0; p < VIM; p++) {
                for (q = 0; q < VIM; q++) {
                  var = v_s[mode][p][q];

                  if (pd->v[pg->imtrx][var]) {
                    pvar = upd->vp[pg->imtrx][var];
                    for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                      phi_j = bf[var]->phi[j];
                      mass = 0.;
                      if (pd->TimeIntegration != STEADY) {
                        if (pd->e[pg->imtrx][eqn] & T_MASS) {
                          mass = (1. + 2. * tt) * phi_j / dt * (double)delta(a, p) *
                                 (double)delta(b, q);
                          mass *= h3 * det_J;
                          mass *= wt_func * at * lambda * wt * pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
                        }
                      }

                      advection = 0.;

                      if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
                        if (DOUBLE_NONZERO(lambda)) {
                          if ((a == p) && (b == q)) {
                            for (r = 0; r < WIM; r++) {
                              advection += (v[r] - x_dot[r]) * bf[var]->grad_phi[j][r];
                            }
                          }
                          advection -=
                              phi_j * ucwt *
                              (gt[a][p] * (double)delta(b, q) + g[q][b] * (double)delta(a, p));
                          advection +=
                              phi_j * lcwt *
                              (gt[q][b] * (double)delta(p, a) + g[a][p] * (double)delta(q, b));

                          advection *= h3 * det_J;

                          advection *= wt_func * wt * at * lambda *
                                       pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
                        }
                      }

                      /*
                       * Diffusion...
                       */

                      diffusion = 0.;

                      if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {
                        diffusion *= det_J * wt * h3;
                        diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
                      }

                      /*
                       * Source term...
                       */

                      source = 0.;

                      if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
                        source_a = Z * (double)delta(a, p) * (double)delta(b, q);
                        if (p == q)
                          source_a += s[a][b] * dZ_dtrace;
                        source_a *= saramitoCoeff;
                        // sensitivities for saramito model:
                        if (p <= q) {
                          source_a += d_saramito->s[p][q] * s[a][b] * Z;
                        }

                        source_b = 0.;
                        if (DOUBLE_NONZERO(alpha)) {
                          source_b =
                              alpha * lambda * saramitoCoeff *
                              (s[q][b] * (double)delta(a, p) + s[a][p] * (double)delta(b, q)) / mup;
                          if (p <= q) {
                            source_b +=
                                d_saramito->s[p][q] * alpha * lambda * (s_dot_s[a][b] / mup);
                          }
                        }

                        source = source_a + source_b;
                        source *= phi_j * det_J * h3 * wt_func * wt *
                                  pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
                      }

                      lec->J[LEC_J_INDEX(peqn, pvar, i, j)] +=
                          mass + advection + diffusion + source;
                    }
                  }
                }
              }
            }
          }
        }
      }
    } /* End Assemble Jacobian */
  }   /* End loop over modes */

  return (status);
}

/*
 * This routine assembles the stress with a log-conformation tensor formulation.
 */
int assemble_stress_log_conf(dbl tt, dbl dt, PG_DATA *pg_data) {
  int dim, q, a, b, w;
  int eqn, siz;

  int i, j, status, mode;
  int logc_gradv = 0;
  dbl v[DIM];
  dbl x_dot[DIM];
  dbl h3;

  dbl grad_v[DIM][DIM];
  dbl gamma[DIM][DIM];
  dbl det_J;

  dbl mass;
  dbl advection;
  dbl source;
  dbl source_term1[DIM][DIM];

  dbl wt_func;
  dbl wt;
  dbl tmp1[DIM][DIM], tmp2[DIM][DIM], tmp3[DIM][DIM];
  dbl advection_term1[DIM][DIM];

  // Variables for stress, velocity gradient
  int R_s[MAX_MODES][DIM][DIM];
  int v_s[MAX_MODES][DIM][DIM];
  dbl s[DIM][DIM], exp_s[DIM][DIM];
  dbl s_dot[DIM][DIM];
  dbl grad_s[DIM][DIM][DIM];
  dbl d_grad_s_dmesh[DIM][DIM][DIM][DIM][MDE];
  dbl gt[DIM][DIM];

  // Temperature shift
  dbl at = 0.0;
  dbl wlf_denom;

  // Consitutive prameters
  dbl alpha;
  dbl lambda = 0;
  dbl d_lambda;
  POLYMER_TIME_CONST_DEPENDENCE_STRUCT d_lam_struct;
  POLYMER_TIME_CONST_DEPENDENCE_STRUCT *d_lam = &d_lam_struct;
  dbl eps;
  dbl Z = 1.0;

  // Decomposition of velocity vector
  dbl M1[DIM][DIM];
  dbl eig_values[DIM];
  dbl R1[DIM][DIM];
  dbl R1_T[DIM][DIM];
  dbl Rt_dot_gradv[DIM][DIM];
  dbl D[DIM][DIM];
  dbl D_dot_D[DIM][DIM];

  // Advective terms
  dbl v_dot_del_s[DIM][DIM];
  dbl x_dot_del_s[DIM][DIM];

  // Trace of stress
  dbl trace = 0.0;

  // SUPG terms
  dbl supg = 0;

  status = 0;
  if (vn->evssModel == LOG_CONF_GRADV) {
    logc_gradv = 1;
  }

  eqn = R_STRESS11;
  // Check if we are actually needed
  if (!pd->e[pg->imtrx][eqn]) {
    return (status);
  }

  dim = pd->Num_Dim;
  wt = fv->wt;
  det_J = bf[eqn]->detJ;
  h3 = fv->h3;

  // Load pointers
  (void)stress_eqn_pointer(v_s);
  (void)stress_eqn_pointer(R_s);

  memset(s, 0, sizeof(double) * DIM * DIM);
  memset(exp_s, 0, sizeof(double) * DIM * DIM);

  // Load up field variables
  for (a = 0; a < WIM; a++) {
    // Velocity
    v[a] = fv->v[a];
    //
    if (pd->TimeIntegration != STEADY && pd->gv[MESH_DISPLACEMENT1 + a]) {
      x_dot[a] = fv_dot->x[a];
    } else {
      x_dot[a] = 0.0;
    }
  }

  // Velocity gradient
  for (a = 0; a < VIM; a++) {
    for (b = 0; b < VIM; b++) {
      grad_v[a][b] = fv->grad_v[a][b];
    }
  }

  // Shear rate
  for (a = 0; a < VIM; a++) {
    for (b = 0; b < VIM; b++) {
      gamma[a][b] = grad_v[a][b] + grad_v[b][a];
    }
  }

  // Velocity gradient projection
  for (a = 0; a < VIM; a++) {
    for (b = 0; b < VIM; b++) {
      gt[a][b] = fv->G[b][a];
    }
  }

  if (vn->wt_funcModel == GALERKIN) {
    supg = 0.0;
  } else if (vn->wt_funcModel == SUPG) {
    supg = vn->wt_func;
  }

  SUPG_terms supg_terms;
  if (supg != 0.0) {
    supg_tau(&supg_terms, dim, 0, pg_data, dt, true, eqn);
  }

  dbl dcdd_factor = 0.0;
  if (vn->shockcaptureModel == SC_DCDD) {
    dcdd_factor = vn->shockcapture;
  } else if (vn->shockcaptureModel != SC_NONE) {
    GOMA_EH(GOMA_ERROR, "Unknown shock capture model, only DCDD supported for LOG_CONF");
  }

  // Shift factor
  if (pd->gv[TEMPERATURE]) {
    if (vn->shiftModel == CONSTANT) {
      at = vn->shift[0];
    } else if (vn->shiftModel == MODIFIED_WLF) {
      wlf_denom = vn->shift[1] + fv->T - mp->reference[TEMPERATURE];
      if (wlf_denom != 0.0) {
        at = exp(vn->shift[0] * (mp->reference[TEMPERATURE] - fv->T) / wlf_denom);
      } else {
        at = 1.0;
      }
    }
  } else {
    at = 1.0;
  }

  // Loop over modes
  for (mode = 0; mode < vn->modes; mode++) {
    // Load up constants and some pointers
    load_modal_pointers(mode, tt, dt, s, s_dot, grad_s, d_grad_s_dmesh);

    // Giesekus mobility parameter
    alpha = ve[mode]->alpha;

    // Polymer time constant
    lambda = polymer_time_const(ve[mode]->time_const_st, gamma, d_lam);

    dbl xi = 0;
    if (ve[mode]->xiModel == CONSTANT) {
      xi = ve[mode]->xi;
    } else if (ls != NULL && ve[mode]->xiModel == VE_LEVEL_SET) {
      double pos_xi = ve[mode]->pos_ls.xi;
      double neg_xi = ve[mode]->xi;
      double width = ls->Length_Scale;
      int err = level_set_property(neg_xi, pos_xi, width, &xi, NULL);
      GOMA_EH(err, "level_set_property() failed for ptt xi parameter.");
    } else {
      GOMA_EH(GOMA_ERROR, "Unknown PTT Xi parameter model");
    }

    if (lambda <= 0.) {
      GOMA_WH(-1, "Trouble: Zero relaxation time with LOG_CONF");
      return -1;
    }

#ifdef ANALEIG_PLEASE
    analytical_exp_s(s, exp_s, eig_values, R1, NULL);
#else
    compute_exp_s(s, exp_s, eig_values, R1);
#endif

    /* Check to make sure eigenvalues are positive (negative eigenvalues will not
       work for log-conformation formulation). These eigenvalues are for the
       conformation tensor, not the log-conformation tensor. */
    if (eig_values[0] < 0. || eig_values[1] < 0. || (VIM > 2 && eig_values[2] < 0.)) {
      GOMA_WH(-1, "Error: Negative eigenvalue for conformation tensor");
      return -1;
    }

    memset(D, 0, sizeof(double) * DIM * DIM);
    D[0][0] = eig_values[0];
    D[1][1] = eig_values[1];
    if (VIM > 2) {
      D[2][2] = eig_values[2];
    }
    (void)tensor_dot(D, D, D_dot_D, VIM);

    // Decompose velocity gradient

    memset(M1, 0, sizeof(double) * DIM * DIM);
    memset(R1_T, 0, sizeof(double) * DIM * DIM);

    for (i = 0; i < VIM; i++) {
      for (j = 0; j < VIM; j++) {
        R1_T[i][j] = R1[j][i];
      }
    }

    for (i = 0; i < VIM; i++) {
      for (j = 0; j < VIM; j++) {
        Rt_dot_gradv[i][j] = 0.;
        for (w = 0; w < VIM; w++) {
          if (DOUBLE_NONZERO(xi)) {
            if (logc_gradv) {
              Rt_dot_gradv[i][j] +=
                  R1_T[i][w] * (grad_v[j][w] - 0.5 * xi * (grad_v[j][w] + grad_v[w][j]));
            } else {
              Rt_dot_gradv[i][j] += R1_T[i][w] * (gt[w][j] - 0.5 * xi * (gt[j][w] + gt[w][j]));
            }
          } else {
            if (logc_gradv) {
              Rt_dot_gradv[i][j] += R1_T[i][w] * grad_v[j][w];
            } else {
              Rt_dot_gradv[i][j] += R1_T[i][w] * gt[w][j];
            }
          }
        }
      }
    }

    for (i = 0; i < VIM; i++) {
      for (j = 0; j < VIM; j++) {
        M1[i][j] = 0.;
        for (w = 0; w < VIM; w++) {
          M1[i][j] += Rt_dot_gradv[i][w] * R1[w][j];
        }
      }
    }

    // Predetermine advective terms
    trace = exp_s[0][0] + exp_s[1][1];
    if (VIM > 2) {
      trace += exp_s[2][2];
    }

    for (a = 0; a < VIM; a++) {
      for (b = 0; b < VIM; b++) {
        v_dot_del_s[a][b] = 0.0;
        x_dot_del_s[a][b] = 0.0;
        for (q = 0; q < dim; q++) {
          v_dot_del_s[a][b] += v[q] * grad_s[q][a][b];
          x_dot_del_s[a][b] += x_dot[q] * grad_s[q][a][b];
        }
      }
    }

    // PTT exponent
    eps = ve[mode]->eps;

    // PTT
    Z = 1;
    if (vn->ConstitutiveEquation == PTT) {
      if (vn->ptt_type == PTT_LINEAR) {
        Z = 1 + eps * (trace - (double)VIM);
      } else if (vn->ptt_type == PTT_EXPONENTIAL) {
        Z = exp(eps * (trace - (double)VIM));
      } else {
        GOMA_EH(GOMA_ERROR, "Unrecognized PTT Form %d", vn->ptt_type);
      }
    }

    siz = sizeof(double) * DIM * DIM;
    memset(tmp1, 0, siz);
    memset(tmp2, 0, siz);
    memset(tmp3, 0, siz);
    memset(advection_term1, 0, siz);
    memset(source_term1, 0, siz);

    for (a = 0; a < VIM; a++) {
      for (b = 0; b < VIM; b++) {
        if (a != b) {
          d_lambda = eig_values[b] - eig_values[a];
          if (DOUBLE_NONZERO(d_lambda)) {
            double eiglog_a = log(DBL_SMALL), eiglog_b = log(DBL_SMALL);
            if (DOUBLE_NONZERO(eig_values[b])) {
              eiglog_b = fmax(eiglog_b, log(eig_values[b]));
            }
            if (DOUBLE_NONZERO(eig_values[a])) {
              eiglog_a = fmax(eiglog_a, log(eig_values[a]));
            }
            tmp1[a][b] += (eiglog_b - eiglog_a) / d_lambda;
            tmp1[a][b] *= (eig_values[a] * M1[b][a] + eig_values[b] * M1[a][b]);
          } else {
            tmp1[a][b] += eig_values[b] * (M1[a][b] + M1[b][a]);
          }
        }
        if (a == b) {
          source_term1[a][b] += Z * (1.0 - D[a][a]) / lambda;
          if (DOUBLE_NONZERO(alpha)) {
            source_term1[a][b] += alpha * (2.0 * D[a][a] - 1.0 - D_dot_D[a][a]) / lambda;
          }
          source_term1[a][b] /= fmax(DBL_SMALL, eig_values[a]);
          source_term1[a][b] += 2.0 * M1[a][a];
        }
      }
    }

    (void)tensor_dot(R1, tmp1, tmp2, VIM);
    (void)tensor_dot(tmp2, R1_T, advection_term1, VIM);
    (void)tensor_dot(R1, source_term1, tmp3, VIM);
    (void)tensor_dot(tmp3, R1_T, source_term1, VIM);

    if (af->Assemble_Residual) {
      for (a = 0; a < VIM; a++) {
        for (b = 0; b < VIM; b++) {
          if (a <= b) {
            eqn = R_s[mode][a][b];

            for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
              wt_func = bf[eqn]->phi[i];

              // SUPG weighting, this is SUPG with s, not e^s
              if (DOUBLE_NONZERO(supg)) {
                for (w = 0; w < dim; w++) {
                  wt_func += supg * supg_terms.supg_tau * v[w] * bf[eqn]->grad_phi[i][w];
                }
              }

              mass = 0.0;
              if (pd->TimeIntegration != STEADY) {
                if (pd->e[pg->imtrx][eqn] & T_MASS) {
                  mass = s_dot[a][b];
                  mass *= wt_func * at * det_J * wt * h3;
                  mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
                }
              }

              advection = 0.;
              if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
                advection += v_dot_del_s[a][b] - x_dot_del_s[a][b];
                advection -= advection_term1[a][b];
                advection *= wt_func * at * det_J * wt * h3;
                advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
              }

              dbl diffusion = 0.;
              if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {
                if (DOUBLE_NONZERO(dcdd_factor)) {
                  dbl tmp = 0.0;
                  dbl s[DIM] = {0.0};
                  dbl r[DIM] = {0.0};
                  for (int w = 0; w < dim; w++) {
                    tmp += (v[w] - x_dot[w]) * (v[w] - x_dot[w]);
                  }
                  tmp = 1.0 / (sqrt(tmp) + 1e-16);
                  for (int w = 0; w < dim; w++) {
                    s[w] = (v[w] - x_dot[w]) * tmp;
                  }
                  dbl mags = 0;
                  for (int w = 0; w < dim; w++) {
                    mags += (grad_s[w][a][b] * grad_s[w][a][b]);
                  }
                  mags = 1.0 / (sqrt(mags) + 1e-14);
                  for (int w = 0; w < dim; w++) {
                    r[w] = grad_s[w][a][b] * mags;
                  }

                  dbl he = 0.0;
                  for (int q = 0; q < ei[pg->imtrx]->dof[eqn]; q++) {
                    dbl tmp = 0;
                    for (int w = 0; w < dim; w++) {
                      tmp += bf[eqn]->grad_phi[q][w] * bf[eqn]->grad_phi[q][w];
                    }
                    he += 1.0 / sqrt(tmp);
                  }

                  dbl G[DIM][DIM];
                  get_metric_tensor(bf[eqn]->B, dim, ei[pg->imtrx]->ielem_type, G);

                  dbl hrgn = 0.0;
                  // for (int w = 0; w < dim; w++) {
                  //   for (int z = 0; z < dim; z++) {
                  //     //tmp += fabs(r[w] * G[w][z] * r[z]);
                  //   }
                  // }
                  tmp = 0;
                  for (int q = 0; q < ei[pg->imtrx]->dof[eqn]; q++) {
                    for (int w = 0; w < dim; w++) {
                      tmp += fabs(r[w] * bf[eqn]->grad_phi[q][w]);
                    }
                  }
                  hrgn = 1.0 / (tmp + 1e-14);

                  dbl magv = 0.0;
                  for (int q = 0; q < VIM; q++) {
                    magv += v[q] * v[q];
                  }
                  magv = sqrt(magv);

                  dbl tau_dcdd = 0.5 * he * magv * magv * pow((1.0 / (mags + 1e-16)) * hrgn, 1.0);
                  // dbl tau_dcdd = (1.0 / mags) * hrgn * hrgn;
                  tau_dcdd = 1 / sqrt(1.0 / (supg_terms.supg_tau * supg_terms.supg_tau + 1e-32) +
                                      (supg_terms.supg_tau * supg_terms.supg_tau *
                                           supg_terms.supg_tau * supg_terms.supg_tau +
                                       1e-32) +
                                      1.0 / (tau_dcdd * tau_dcdd + 1e-32));
                  dbl ss[DIM][DIM] = {{0.0}};
                  dbl rr[DIM][DIM] = {{0.0}};
                  dbl rdots = 0.0;
                  for (int w = 0; w < dim; w++) {
                    for (int z = 0; z < dim; z++) {
                      ss[w][z] = s[w] * s[z];
                      rr[w][z] = r[w] * r[z];
                    }
                    rdots += r[w] * s[w];
                  }

                  dbl inner_tensor[DIM][DIM] = {{0.0}};
                  for (int w = 0; w < dim; w++) {
                    for (int z = 0; z < dim; z++) {
                      inner_tensor[w][z] = rr[w][z] - rdots * rdots * ss[w][z];
                    }
                  }

                  dbl gs_inner_dot[DIM] = {0.0};
                  for (int w = 0; w < dim; w++) {
                    dbl tmp = 0.;
                    for (int z = 0; z < dim; z++) {
                      tmp += grad_s[w][a][b] * inner_tensor[w][z];
                    }
                    gs_inner_dot[w] = tmp;
                    // gs_inner_dot[w] = grad_s[w][a][b];
                  }

                  for (int w = 0; w < dim; w++) {
                    diffusion += tau_dcdd * gs_inner_dot[w] * bf[eqn]->grad_phi[i][w];
                  }
                  diffusion *= dcdd_factor * det_J * wt * h3;
                }
                diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
              }

              source = 0.0;
              if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
                source -= source_term1[a][b];
                source *= wt_func * det_J * h3 * wt;
                source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
              }
              lec->R[LEC_R_INDEX(upd->ep[pg->imtrx][eqn], i)] +=
                  mass + advection + diffusion + source;
            } // i loop
          }   // if a<=b
        }     // b loop
      }       // a loop
    }         // if Residual
  }
  return (status);
}

int assemble_stress_log_conf_transient(dbl tt, dbl dt, PG_DATA *pg_data) {
  int dim, p, q, a, b, w;
  int eqn, siz;

  int i, j, status, mode;
  int logc_gradv = 0;
  dbl v[DIM];
  dbl x_dot[DIM];
  dbl h3;

  dbl grad_v[DIM][DIM];
  dbl gamma[DIM][DIM];
  dbl det_J;

  dbl mass;
  dbl advection;
  dbl source;
  dbl source_term1[DIM][DIM];

  dbl wt_func;
  dbl wt;
  dbl tmp1[DIM][DIM], tmp2[DIM][DIM], tmp3[DIM][DIM];
  dbl advection_term1[DIM][DIM];

  // Variables for stress, velocity gradient
  int R_s[MAX_MODES][DIM][DIM];
  int v_s[MAX_MODES][DIM][DIM];
  dbl s[DIM][DIM], exp_s[DIM][DIM];
  dbl s_dot[DIM][DIM];
  dbl grad_s[DIM][DIM][DIM];
  dbl d_grad_s_dmesh[DIM][DIM][DIM][DIM][MDE];
  dbl gt[DIM][DIM];

  // Polymer viscosity
  dbl mup;
  VISCOSITY_DEPENDENCE_STRUCT d_mup_struct;
  VISCOSITY_DEPENDENCE_STRUCT *d_mup = &d_mup_struct;

  POLYMER_TIME_CONST_DEPENDENCE_STRUCT d_lam_struct;
  POLYMER_TIME_CONST_DEPENDENCE_STRUCT *d_lam = &d_lam_struct;
  // Temperature shift
  dbl at = 0.0;
  dbl wlf_denom;

  // Consitutive prameters
  dbl alpha;
  dbl lambda = 0;
  dbl d_lambda;
  dbl eps;
  dbl Z = 1.0;

  // Decomposition of velocity vector
  dbl M1[DIM][DIM];
  dbl eig_values[DIM];
  dbl R1[DIM][DIM];
  dbl R1_T[DIM][DIM];
  dbl Rt_dot_gradv[DIM][DIM];
  dbl D[DIM][DIM];
  dbl D_dot_D[DIM][DIM];

  // Advective terms
  dbl v_dot_del_s[DIM][DIM];
  dbl x_dot_del_s[DIM][DIM];

  // Trace of stress
  dbl trace = 0.0;

  // SUPG terms
  dbl supg = 0;

  status = 0;
  if (vn->evssModel == LOG_CONF_TRANSIENT_GRADV) {
    logc_gradv = 1;
  }

  eqn = R_STRESS11;
  // Check if we are actually needed
  if (!pd->e[pg->imtrx][eqn]) {
    return (status);
  }

  dim = pd->Num_Dim;
  wt = fv->wt;
  det_J = bf[eqn]->detJ;
  h3 = fv->h3;

  // Load pointers
  (void)stress_eqn_pointer(v_s);
  (void)stress_eqn_pointer(R_s);

  memset(s, 0, sizeof(double) * DIM * DIM);
  memset(exp_s, 0, sizeof(double) * DIM * DIM);

  // Load up field variables
  for (a = 0; a < dim; a++) {
    // Velocity
    v[a] = fv->v[a];
    //
    if (pd->TimeIntegration != STEADY && pd->gv[MESH_DISPLACEMENT1 + a]) {
      x_dot[a] = fv_dot->x[a];
    } else {
      x_dot[a] = 0.0;
    }
  }

  // Velocity gradient
  for (a = 0; a < VIM; a++) {
    for (b = 0; b < VIM; b++) {
      grad_v[a][b] = fv->grad_v[a][b];
    }
  }

  // Shear rate
  for (a = 0; a < VIM; a++) {
    for (b = 0; b < VIM; b++) {
      gamma[a][b] = grad_v[a][b] + grad_v[b][a];
    }
  }

  // Velocity gradient projection
  for (a = 0; a < VIM; a++) {
    for (b = 0; b < VIM; b++) {
      gt[a][b] = fv->G[b][a];
    }
  }

  if (vn->wt_funcModel == GALERKIN) {
    supg = 0.0;
  } else if (vn->wt_funcModel == SUPG) {
    supg = vn->wt_func;
  }

  const bool saramitoEnabled =
      (vn->ConstitutiveEquation == SARAMITO_OLDROYDB || vn->ConstitutiveEquation == SARAMITO_PTT ||
       vn->ConstitutiveEquation == SARAMITO_GIESEKUS);

  dbl saramitoCoeff = 1.;

  SUPG_terms supg_terms;
  if (supg != 0.0) {
    supg_tau_shakib(&supg_terms, dim, dt, 1e-6, POLYMER_STRESS11);
  }

  // Shift factor
  if (pd->gv[TEMPERATURE]) {
    if (vn->shiftModel == CONSTANT) {
      at = vn->shift[0];
    } else if (vn->shiftModel == MODIFIED_WLF) {
      wlf_denom = vn->shift[1] + fv->T - mp->reference[TEMPERATURE];
      if (wlf_denom != 0.0) {
        at = exp(vn->shift[0] * (mp->reference[TEMPERATURE] - fv->T) / wlf_denom);
      } else {
        at = 1.0;
      }
    }
  } else {
    at = 1.0;
  }

  // Loop over modes
  for (mode = 0; mode < vn->modes; mode++) {
    // Load up constants and some pointers
    load_modal_pointers(mode, tt, dt, s, s_dot, grad_s, d_grad_s_dmesh);

    // Polymer viscosity
    mup = viscosity(ve[mode]->gn, gamma, d_mup);

    // Giesekus mobility parameter
    alpha = ve[mode]->alpha;

    // Polymer time constant
    lambda = polymer_time_const(ve[mode]->time_const_st, gamma, d_lam);

#ifdef ANALEIG_PLEASE
    analytical_exp_s(fv_old->S[mode], exp_s, eig_values, R1, NULL);
#else
    compute_exp_s(fv_old->S[mode], exp_s, eig_values, R1);
#endif

    dbl tau[DIM][DIM] = {{0.0}};
    if (saramitoEnabled == TRUE) {
      for (int i = 0; i < VIM; i++) {
        for (int j = 0; j < VIM; j++) {
          tau[i][j] = mup / lambda * (exp_s[i][j] - delta(i, j));
        }
      }
      compute_saramito_model_terms(&saramitoCoeff, NULL, tau, ve[mode]->gn, FALSE);
    } else {
      saramitoCoeff = 1.;
    }

    /* Check to make sure eigenvalues are positive (negative eigenvalues will not
       work for log-conformation formulation). These eigenvalues are for the
       conformation tensor, not the log-conformation tensor. */
    if (eig_values[0] < 0. || eig_values[1] < 0. || (VIM > 2 && eig_values[2] < 0.)) {
      GOMA_WH(-1, "Error: Negative eigenvalue for conformation tensor");
      return -1;
    }

    memset(D, 0, sizeof(double) * DIM * DIM);
    D[0][0] = eig_values[0];
    D[1][1] = eig_values[1];
    if (VIM > 2) {
      D[2][2] = eig_values[2];
    }
    (void)tensor_dot(D, D, D_dot_D, VIM);

    // Decompose velocity gradient

    memset(M1, 0, sizeof(double) * DIM * DIM);
    memset(R1_T, 0, sizeof(double) * DIM * DIM);

    for (i = 0; i < VIM; i++) {
      for (j = 0; j < VIM; j++) {
        R1_T[i][j] = R1[j][i];
      }
    }

    for (i = 0; i < VIM; i++) {
      for (j = 0; j < VIM; j++) {
        Rt_dot_gradv[i][j] = 0.;
        for (w = 0; w < VIM; w++) {
          if (logc_gradv) {
            Rt_dot_gradv[i][j] += R1_T[i][w] * grad_v[j][w];
          } else {
            Rt_dot_gradv[i][j] += R1_T[i][w] * gt[w][j];
          }
        }
      }
    }

    for (i = 0; i < VIM; i++) {
      for (j = 0; j < VIM; j++) {
        M1[i][j] = 0.;
        for (w = 0; w < VIM; w++) {
          M1[i][j] += Rt_dot_gradv[i][w] * R1[w][j];
        }
      }
    }

    // Predetermine advective terms
    trace = exp_s[0][0] + exp_s[1][1];
    if (VIM > 2) {
      trace += exp_s[2][2];
    }

    for (a = 0; a < VIM; a++) {
      for (b = 0; b < VIM; b++) {
        v_dot_del_s[a][b] = 0.0;
        x_dot_del_s[a][b] = 0.0;
        for (q = 0; q < dim; q++) {
          v_dot_del_s[a][b] += v[q] * grad_s[q][a][b];
          x_dot_del_s[a][b] += x_dot[q] * grad_s[q][a][b];
        }
      }
    }

    // PTT exponent
    eps = ve[mode]->eps;

    // Exponential term for PTT
    Z = exp(eps * (trace - (double)VIM));

    siz = sizeof(double) * DIM * DIM;
    memset(tmp1, 0, siz);
    memset(tmp2, 0, siz);
    memset(tmp3, 0, siz);
    memset(advection_term1, 0, siz);
    memset(source_term1, 0, siz);

    for (a = 0; a < VIM; a++) {
      for (b = 0; b < VIM; b++) {
        if (a != b) {
          d_lambda = eig_values[b] - eig_values[a];
          if (DOUBLE_NONZERO(d_lambda)) {
            double eiglog_a = log(DBL_SMALL), eiglog_b = log(DBL_SMALL);
            if (DOUBLE_NONZERO(eig_values[b])) {
              eiglog_b = fmax(eiglog_b, log(eig_values[b]));
            }
            if (DOUBLE_NONZERO(eig_values[a])) {
              eiglog_a = fmax(eiglog_a, log(eig_values[a]));
            }
            tmp1[a][b] += (eiglog_b - eiglog_a) / d_lambda;
            tmp1[a][b] *= (eig_values[a] * M1[b][a] + eig_values[b] * M1[a][b]);
          } else {
            tmp1[a][b] += eig_values[b] * (M1[a][b] + M1[b][a]);
          }
        }
        if (a == b) {
          source_term1[a][b] += saramitoCoeff * Z * (1.0 - D[a][a]) / lambda;
          if (DOUBLE_NONZERO(alpha)) {
            source_term1[a][b] += alpha * (2.0 * D[a][a] - 1.0 - D_dot_D[a][a]) / lambda;
          }
          source_term1[a][b] /= fmax(DBL_SMALL, eig_values[a]);
          source_term1[a][b] += 2.0 * M1[a][a];
        }
      }
    }

    (void)tensor_dot(R1, tmp1, tmp2, VIM);
    (void)tensor_dot(tmp2, R1_T, advection_term1, VIM);
    (void)tensor_dot(R1, source_term1, tmp3, VIM);
    (void)tensor_dot(tmp3, R1_T, source_term1, VIM);

    if (af->Assemble_Residual) {
      for (a = 0; a < VIM; a++) {
        for (b = 0; b < VIM; b++) {
          if (a <= b) {
            eqn = R_s[mode][a][b];

            for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
              wt_func = bf[eqn]->phi[i];

              // SUPG weighting, this is SUPG with s, not e^s
              if (DOUBLE_NONZERO(supg)) {
                for (w = 0; w < dim; w++) {
                  wt_func += supg * supg_terms.supg_tau * v[w] * bf[eqn]->grad_phi[i][w];
                }
              }

              mass = 0.0;
              if (pd->TimeIntegration != STEADY) {
                if (pd->e[pg->imtrx][eqn] & T_MASS) {
                  mass = fv_dot->S[mode][a][b];
                  mass *= wt_func * at * det_J * wt * h3;
                  mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
                }
              }

              advection = 0.;
              if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
                advection += v_dot_del_s[a][b] - x_dot_del_s[a][b];
                advection -= advection_term1[a][b];
                advection *= wt_func * at * det_J * wt * h3;
                advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
              }

              source = 0.0;
              if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
                source -= source_term1[a][b];
                source *= wt_func * det_J * h3 * wt;
                source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
              }
              lec->R[LEC_R_INDEX(upd->ep[pg->imtrx][eqn], i)] += mass + advection + source;
            } // i loop
          }   // if a<=b
        }     // b loop
      }       // a loop
    }         // if Residual
    if (af->Assemble_Jacobian) {
      for (a = 0; a < VIM; a++) {
        for (b = 0; b < VIM; b++) {
          if (a <= b) {
            eqn = R_s[mode][a][b];
            int peqn = upd->ep[pg->imtrx][eqn];

            for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
              wt_func = bf[eqn]->phi[i];

              // SUPG weighting, this is SUPG with s, not e^s
              if (DOUBLE_NONZERO(supg)) {
                for (w = 0; w < dim; w++) {
                  wt_func += supg * supg_terms.supg_tau * v[w] * bf[eqn]->grad_phi[i][w];
                }
              }

              /*
               * J_S_S
               */
              for (p = 0; p < VIM; p++) {
                for (q = 0; q < VIM; q++) {
                  int var = v_s[mode][p][q];

                  if (pd->v[pg->imtrx][var]) {
                    int pvar = upd->vp[pg->imtrx][var];
                    for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                      double phi_j = bf[var]->phi[j];
                      mass = 0.;
                      if (pd->TimeIntegration != STEADY) {
                        if (pd->e[pg->imtrx][eqn] & T_MASS) {
                          mass = (1. + 2. * tt) * phi_j / dt * (double)delta(a, p) *
                                 (double)delta(b, q);
                          mass *= h3 * det_J;
                          mass *= wt_func * wt * pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
                        }
                      }

                      advection = 0.;

                      if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
                        if ((a == p) && (b == q)) {
                          for (int r = 0; r < dim; r++) {
                            advection += (v[r] - x_dot[r]) * bf[var]->grad_phi[j][r];
                          }
                        }

                        advection *= h3 * det_J;

                        advection *= wt_func * wt * pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
                      }

                      source = 0;
                      lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += mass + advection + source;
                    }
                  }
                }
              }
            } // i loop
          }   // if a<=b
        }     // b loop
      }       // a loop
    }         // if Residual
  }
  return (status);
}

/* this stress routine uses the viscoelastic equations to do a solid-fluid
 * interaction problem in an Eulerian context with the level set denote the
 * solid-fluid interface.
 */

int assemble_stress_level_set(dbl tt,           /* parameter to vary time integration from
                                                 * explicit (tt = 1) to implicit (tt = 0) */
                              dbl dt,           /* current time step size */
                              dbl h[DIM],       /* coordinate scale factors */
                              dbl hh[DIM][DIM], /* coordinate scale factors */
                              dbl dh_dxnode[DIM][MDE],
                              dbl vcent[DIM], /* Average element velocity, which is
                                               * the centroid velocity for Q2 and the
                                               * average of the vertices for Q1. It
                                               * comes from the routine
                                               * "element_velocity." */
                              dbl dvc_dnode[DIM][MDE]) {
  int dim, p, q, r, a, b, w;

  int eqn, var;
  int peqn, pvar;

  int i, j, mode;
  dbl v[DIM];      /* Velocity field. */
  dbl x_dot[DIM];  /* current position field derivative wrt time. */
  dbl h3;          /* Volume element (scale factors). */
  dbl dh3dmesh_pj; /* Sensitivity to (p,j) mesh dof. */

  dbl grad_v[DIM][DIM];
  dbl gamma[DIM][DIM]; /* Shear-rate tensor based on velocity */
  dbl det_J;           /* determinant of element Jacobian */

  dbl d_det_J_dmesh_pj; /* for specific (p,j) mesh dof */

  dbl mass; /* For terms and their derivatives */
  dbl mass_a, mass_b;
  dbl advection;
  dbl advection_a, advection_b, advection_c, advection_d;
  dbl source;
  dbl source_a = 0, source_b = 0, source_c = 0;

  /*
   *
   * Note how carefully we avoid refering to d(phi[i])/dx[j] and refer instead
   * to the j-th component of grad_phi[j][i] so that this vector can be loaded
   * up with components that may be different in non Cartesian coordinate
   * systems.
   *
   * We will, however, insist on *orthogonal* coordinate systems, even if we
   * might permit them to be curvilinear.
   *
   * Assume all components of velocity are interpolated with the same kind
   * of basis function.
   */

  /*
   * Petrov-Galerkin weighting functions for i-th and ab-th stress residuals
   * and some of their derivatives...
   */

  dbl wt_func;

  /*
   * Interpolation functions for variables and some of their derivatives.
   */

  dbl phi_j;
  dbl wt;

  /* Variables for stress */

  int R_s[MAX_MODES][DIM][DIM];
  int v_s[MAX_MODES][DIM][DIM];
  int v_g[DIM][DIM];

  dbl s[DIM][DIM];     /* stress tensor */
  dbl s_dot[DIM][DIM]; /* stress tensor from last time step */
  dbl grad_s[DIM][DIM][DIM];
  dbl d_grad_s_dmesh[DIM][DIM][DIM][DIM]
                    [MDE]; /* derivative of grad of stress tensor for mode ve_mode */

  dbl g[DIM][DIM];  /* velocity gradient tensor */
  dbl gt[DIM][DIM]; /* transpose of velocity gradient tensor */

  /* dot product tensors */

  dbl s_dot_g[DIM][DIM];
  dbl gt_dot_s[DIM][DIM];

  /* polymer viscosity and derivatives */
  dbl mup;
  VISCOSITY_DEPENDENCE_STRUCT d_mup_struct;
  VISCOSITY_DEPENDENCE_STRUCT *d_mup = &d_mup_struct;

  dbl d_mup_dv_pj;
  dbl d_mup_dmesh_pj;

  /* advective terms are precalculated */
  dbl v_dot_del_s[DIM][DIM];
  dbl x_dot_del_s[DIM][DIM];

  dbl d_xdotdels_dm;

  dbl d_vdotdels_dm;

  /* SUPG variables */
  dbl h_elem = 0, h_elem_inv = 0, h_elem_deriv = 0;
  dbl supg = 0;

  int status = 0;

  double H_ls, delta_ls, normal_ls[MAX_PDIM];
  int near_ls;
  dbl Gmod;

  eqn = R_STRESS11;

  /*
   * Bail out fast if there's nothing to do...
   */

  if (!pd->e[pg->imtrx][eqn]) {
    return (status);
  }

  /*
   * Unpack variables from structures for local convenience...
   */

  dim = pd->Num_Dim;

  wt = fv->wt;

  det_J = bf[eqn]->detJ; /* Really, ought to be mesh eqn. */

  h3 = fv->h3; /* Differential volume element (scales). */

  /* load eqn and variable number in tensor form */
  (void)stress_eqn_pointer(v_s);
  (void)stress_eqn_pointer(R_s);

  v_g[0][0] = VELOCITY_GRADIENT11;
  v_g[0][1] = VELOCITY_GRADIENT12;
  v_g[1][0] = VELOCITY_GRADIENT21;
  v_g[1][1] = VELOCITY_GRADIENT22;
  v_g[0][2] = VELOCITY_GRADIENT13;
  v_g[1][2] = VELOCITY_GRADIENT23;
  v_g[2][0] = VELOCITY_GRADIENT31;
  v_g[2][1] = VELOCITY_GRADIENT32;
  v_g[2][2] = VELOCITY_GRADIENT33;

  /*
   * Field variables...
   */

  for (a = 0; a < WIM; a++) {
    v[a] = fv->v[a];

    /* note, these are zero for steady calculations */
    if (pd->TimeIntegration != STEADY && pd->v[pg->imtrx][MESH_DISPLACEMENT1 + a]) {
      x_dot[a] = fv_dot->x[a];
    } else {
      x_dot[a] = 0.;
    }
  }

  /*
   * In Cartesian coordinates, this velocity gradient tensor will
   * have components that are...
   *
   * 			grad_v[a][b] = d v_b
   *				       -----
   *				       d x_a
   */

  for (a = 0; a < VIM; a++) {
    for (b = 0; b < VIM; b++) {
      grad_v[a][b] = fv->grad_v[a][b];
    }
  }

  /* load up shearrate tensor based on velocity */
  for (a = 0; a < VIM; a++) {
    for (b = 0; b < VIM; b++) {
      gamma[a][b] = grad_v[a][b] + grad_v[b][a];
    }
  }

  for (a = 0; a < VIM; a++) {
    for (b = 0; b < VIM; b++) {
      g[a][b] = fv->G[a][b];
      gt[b][a] = g[a][b];
    }
  }

  if (vn->wt_funcModel == GALERKIN) {
    supg = 0.;
  } else if (vn->wt_funcModel == SUPG) {
    supg = vn->wt_func;
  }

  if (supg != 0.) {
    h_elem = 0.;
    for (p = 0; p < dim; p++) {
      h_elem += vcent[p] * vcent[p] * h[p];
    }
    h_elem = sqrt(h_elem) / 2.;
    if (h_elem == 0.) {
      h_elem_inv = 1.;
    } else {
      h_elem_inv = 1. / h_elem;
    }
  }
  /* end Petrov-Galerkin addition */

  /* Get level set length scales, Heaviside function etc */
  level_set_interface(fv->F, fv->grad_F, ls->Length_Scale, 0, &near_ls, &H_ls, NULL, NULL,
                      &delta_ls, NULL, NULL, normal_ls, NULL, NULL);

  /* Begin loop over modes */
  for (mode = 0; mode < vn->modes; mode++) {

    load_modal_pointers(mode, tt, dt, s, s_dot, grad_s, d_grad_s_dmesh);

    /* precalculate advective terms of form (v dot del tensor)*/

    for (a = 0; a < VIM; a++) {
      for (b = 0; b < VIM; b++) {
        v_dot_del_s[a][b] = 0.;
        x_dot_del_s[a][b] = 0.;
        for (q = 0; q < dim; q++) {
          v_dot_del_s[a][b] += v[q] * grad_s[q][a][b];
          x_dot_del_s[a][b] += x_dot[q] * grad_s[q][a][b];
        }
      }
    }

    /*
     * Stress tensor...(Note "anti-BSL" sign convention on deviatoric stress)
     */

    /* get polymer viscosity */
    mup = viscosity(ve[mode]->gn, gamma, d_mup);

    /* hardwire shear modulus for now */
    /* Gmod = 1000000.; */
    Gmod = mup / ve[mode]->time_const_st->lambda0;

    /* get tensor dot products for future use */
    (void)tensor_dot(s, g, s_dot_g, VIM);
    (void)tensor_dot(gt, s, gt_dot_s, VIM);

    /* get tensor dot products for future use */

    /*
     * Residuals_________________________________________________________________
     */

    if (af->Assemble_Residual) {
      /*
       * Assemble each component "ab" of the polymer stress equation...
       */
      for (a = 0; a < VIM; a++) {
        for (b = 0; b < VIM; b++) {

          if (a <= b) /* since the stress tensor is symmetric, only assemble the upper half */
          {
            eqn = R_s[mode][a][b];

            /*
             * In the element, there will be contributions to this many equations
             * based on the number of degrees of freedom...
             */

            for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
              wt_func = bf[eqn]->phi[i];
              /* add Petrov-Galerkin terms as necessary */
              if (supg != 0.) {
                for (w = 0; w < dim; w++) {
                  wt_func += supg * h_elem * v[w] * bf[eqn]->grad_phi[i][w];
                }
              }

              /* The mass and advection terms constitute the substantial
               * derivative and the upper convective derivative of the
               * stress tensor. These terms are divided by the modulus,
               * G, and form the elastic solid response.
               */
              mass = 0.;
              if (pd->TimeIntegration != STEADY) {
                if (pd->e[pg->imtrx][eqn] & T_MASS) {
                  mass = s_dot[a][b] * H_ls / Gmod;
                  mass *= wt_func * det_J * wt * h3;
                  mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
                }
              }

              advection = 0.;
              if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
                advection += v_dot_del_s[a][b] - x_dot_del_s[a][b];
                advection -= (gt_dot_s[a][b] + s_dot_g[a][b]);

                advection *= wt_func * H_ls / Gmod * det_J * wt * h3;
                advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
              }

              /*
               * Source term...
               */

              source = 0.;
              if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
                /* This source term involves the Newtonian fluid
                 * viscosity expressed in the mixed formulation.
                 */
                source = s[a][b] * (1. - H_ls) / mup;

                /* This source term involves the rate of deformation
                 * tensor, which is used for both the solid and fluid
                 * formulations.
                 */

                source -= (g[a][b] + gt[a][b]);
                source *= wt_func * det_J * h3 * wt;
                source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
              }

              /*
               * Add contributions to this residual (globally into Resid, and
               * locally into an accumulator)
               */

              lec->R[LEC_R_INDEX(upd->ep[pg->imtrx][eqn], i)] += mass + advection + source;
            }
          }
        }
      }
    }

    /*
     * Jacobian terms...
     */

    if (af->Assemble_Jacobian) {
      dbl R_source, R_advection; /* Places to put the raw residual portions
                                            instead of constantly recalcing them */
      for (a = 0; a < VIM; a++) {
        for (b = 0; b < VIM; b++) {
          if (a <= b) /* since the stress tensor is symmetric, only assemble the upper half */
          {
            eqn = R_s[mode][a][b];
            peqn = upd->ep[pg->imtrx][eqn];

            R_advection = v_dot_del_s[a][b] - x_dot_del_s[a][b];
            R_advection -= (gt_dot_s[a][b] + s_dot_g[a][b]);
            R_advection *= H_ls / Gmod;

            R_source = s[a][b] * (1. - H_ls) / mup - (g[a][b] + gt[a][b]);

            for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {

              wt_func = bf[eqn]->phi[i];
              /* add Petrov-Galerkin terms as necessary */
              if (supg != 0.) {
                for (w = 0; w < dim; w++) {
                  wt_func += supg * h_elem * v[w] * bf[eqn]->grad_phi[i][w];
                }
              }

              /*
               * Set up some preliminaries that are needed for the (a,i)
               * equation for bunches of (b,j) column variables...
               */

              /*
               * J_S_T
               */

              var = TEMPERATURE;
              if (pd->v[pg->imtrx][var]) {
                pvar = upd->vp[pg->imtrx][var];
                for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                  phi_j = bf[var]->phi[j];

                  source = 0.;
                  if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
                    source -= s[a][b] * (1. - H_ls) / (mup * mup) * d_mup->T[j];
                    source *= wt_func * det_J * wt * h3;
                    source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
                  }

                  lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += source;
                }
              }

              /*
               * J_S_v
               */
              for (p = 0; p < WIM; p++) {
                var = VELOCITY1 + p;
                if (pd->v[pg->imtrx][var]) {
                  pvar = upd->vp[pg->imtrx][var];
                  for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                    phi_j = bf[var]->phi[j];
                    d_mup_dv_pj = d_mup->v[p][j];

                    mass = 0.;

                    if (pd->TimeIntegration != STEADY) {
                      if (pd->e[pg->imtrx][eqn] & T_MASS) {
                        if (supg != 0.) {
                          mass = supg * h_elem * phi_j * bf[eqn]->grad_phi[i][p];

                          for (w = 0; w < dim; w++) {
                            mass += supg * vcent[p] * dvc_dnode[p][j] * h[p] * h_elem_inv / 4. *
                                    v[w] * bf[eqn]->grad_phi[i][w];
                          }

                          mass *= s_dot[a][b] * H_ls / Gmod;
                        }

                        mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)] * det_J * wt * h3;
                      }
                    }

                    advection = 0.;

                    if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
                      advection_a = phi_j * (grad_s[p][a][b]) * H_ls / Gmod;

                      advection_a *= wt_func;

                      advection_b = 0.;
                      /* Petrov-Galerkin term */
                      if (supg != 0.) {
                        advection_b = supg * h_elem * phi_j * bf[eqn]->grad_phi[i][p];
                        for (w = 0; w < dim; w++) {
                          advection_b += supg * vcent[p] * h[p] * dvc_dnode[p][j] * h_elem_inv /
                                         4. * v[w] * bf[eqn]->grad_phi[i][w];
                        }

                        advection_b *= R_advection;
                      }

                      advection = advection_a + advection_b;
                      advection *= det_J * wt * h3;
                      advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
                    }

                    source = 0.;

                    if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
                      source_a = -d_mup_dv_pj * (s[a][b] * (1. - H_ls) / (mup * mup));
                      source_a *= wt_func;

                      source_b = 0.;
                      if (supg != 0.) {
                        source_b = supg * h_elem * phi_j * bf[eqn]->grad_phi[i][p];

                        for (w = 0; w < dim; w++) {
                          source_b += supg * vcent[p] * dvc_dnode[p][j] * h[p] * h_elem_inv / 4. *
                                      v[w] * bf[eqn]->grad_phi[i][w];
                        }

                        source_b *= R_source;
                      }

                      source = source_a + source_b;
                      source *= det_J * wt * h3;
                      source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
                    }

                    lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += mass + advection + source;
                  }
                }
              }

              /*
               * J_S_c
               */
              var = MASS_FRACTION;
              if (pd->v[pg->imtrx][var]) {
                for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                  phi_j = bf[var]->phi[j];

                  for (w = 0; w < pd->Num_Species_Eqn; w++) {

                    source = 0.;

                    if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
                      source = -d_mup->C[w][j] * (s[a][b] * (1. - H_ls) / (mup * mup));
                      source *= wt_func * det_J * wt * h3;
                      source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
                    }

                    if (w > 1) {
                      GOMA_EH(GOMA_ERROR, "Need more arrays for each species.");
                    }

                    lec->J[LEC_J_INDEX(peqn, MAX_PROB_VAR + w, i, j)] += source;
                  }
                }
              }

              /*
               * J_S_P
               */
              var = PRESSURE;
              if (pd->v[pg->imtrx][var]) {
                pvar = upd->vp[pg->imtrx][var];
                for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                  phi_j = bf[var]->phi[j];

                  source = 0.;
                  if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
                    source = -d_mup->P[j] * (s[a][b] * (1. - H_ls) / (mup * mup));
                    source *= wt_func * det_J * wt * h3;
                    source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
                  }

                  lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += source;
                }
              }

              /*
               * J_S_d
               */
              for (p = 0; p < dim; p++) {
                var = MESH_DISPLACEMENT1 + p;
                if (pd->v[pg->imtrx][var]) {
                  pvar = upd->vp[pg->imtrx][var];
                  for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                    phi_j = bf[var]->phi[j];
                    d_det_J_dmesh_pj = bf[eqn]->d_det_J_dm[p][j];
                    dh3dmesh_pj = fv->dh3dq[p] * bf[var]->phi[j];
                    d_mup_dmesh_pj = d_mup->X[p][j];
                    if (supg != 0.) {
                      h_elem_deriv = 0.;
                      for (q = 0; q < dim; q++) {
                        h_elem_deriv +=
                            hh[q][p] * vcent[q] * vcent[q] * dh_dxnode[q][j] * h_elem_inv / 4.;
                      }
                    }

                    mass = 0.;
                    mass_a = 0.;
                    mass_b = 0.;
                    if (pd->TimeIntegration != STEADY) {
                      if (pd->e[pg->imtrx][eqn] & T_MASS) {
                        mass_a = s_dot[a][b] * H_ls / Gmod;
                        mass_a *= wt_func * (d_det_J_dmesh_pj * h3 + det_J * dh3dmesh_pj);

                        if (supg != 0.) {
                          for (w = 0; w < dim; w++) {
                            mass_b +=
                                supg * (h_elem * v[w] * bf[eqn]->d_grad_phi_dmesh[i][w][p][j] +
                                        h_elem_deriv * v[w] * bf[eqn]->grad_phi[i][w]);
                          }
                          mass_b *= s_dot[a][b] * H_ls / Gmod * h3 * det_J;
                        }

                        mass = mass_a + mass_b;
                        mass *= wt * pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
                      }
                    }

                    advection = 0.;

                    if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
                      /*
                       * Four parts:
                       *    advection_a =
                       *    	Int ( ea.(v-xdot).d(Vv)/dmesh h3 |Jv| )
                       *
                       *    advection_b =
                       *  (i)	Int ( ea.(v-xdot).Vv h3 d(|Jv|)/dmesh )
                       *  (ii)  Int ( ea.(v-xdot).d(Vv)/dmesh h3 |Jv| )
                       *  (iii) Int ( ea.(v-xdot).Vv dh3/dmesh |Jv|   )
                       *
                       * For unsteady problems, we have an
                       * additional term
                       *
                       *    advection_c =
                       *    	Int ( ea.d(v-xdot)/dmesh.Vv h3 |Jv| )
                       */

                      advection_a = R_advection;

                      advection_a *= wt_func * (d_det_J_dmesh_pj * h3 + det_J * dh3dmesh_pj);

                      d_vdotdels_dm = 0.;
                      for (q = 0; q < dim; q++) {
                        d_vdotdels_dm += (v[q] - x_dot[q]) * d_grad_s_dmesh[q][a][b][p][j];
                      }

                      advection_b = d_vdotdels_dm;
                      advection_b *= H_ls / Gmod * wt_func * det_J * h3;

                      advection_c = 0.;
                      if (pd->TimeIntegration != STEADY) {
                        if (pd->e[pg->imtrx][eqn] & T_MASS) {
                          d_xdotdels_dm = (1. + 2. * tt) * phi_j / dt * grad_s[p][a][b];

                          advection_c -= d_xdotdels_dm;

                          advection_c *= H_ls / Gmod * wt_func * h3 * det_J;
                        }
                      }

                      advection_d = 0.;
                      if (supg != 0.) {
                        for (w = 0; w < dim; w++) {
                          advection_d +=
                              supg * (h_elem * v[w] * bf[eqn]->d_grad_phi_dmesh[i][w][p][j] +
                                      h_elem_deriv * v[w] * bf[eqn]->grad_phi[i][w]);
                        }

                        advection_d *= (R_advection)*det_J * h3;
                      }

                      advection = advection_a + advection_b + advection_c + advection_d;

                      advection *= wt * pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
                    }

                    /*
                     * Source term...
                     */

                    source = 0.;

                    if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
                      source_a = R_source;
                      source_a *= wt_func * (d_det_J_dmesh_pj * h3 + det_J * dh3dmesh_pj);

                      source_b = -s[a][b] * (1. - H_ls) / (mup * mup);
                      source_b *= wt_func * det_J * h3 * d_mup_dmesh_pj;

                      source_c = 0.;
                      if (supg != 0.) {
                        for (w = 0; w < dim; w++) {
                          source_c +=
                              supg * (h_elem * v[w] * bf[eqn]->d_grad_phi_dmesh[i][w][p][j] +
                                      h_elem_deriv * v[w] * bf[eqn]->grad_phi[i][w]);
                        }
                        source_c *= R_source * det_J * h3;
                      }

                      source = source_a + source_b + source_c;

                      source *= wt * pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
                    }

                    lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += mass + advection + source;
                  }
                }
              }

              /*
               * J_S_G
               */
              for (p = 0; p < VIM; p++) {
                for (q = 0; q < VIM; q++) {
                  var = v_g[p][q];

                  if (pd->v[pg->imtrx][var]) {
                    pvar = upd->vp[pg->imtrx][var];
                    for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                      phi_j = bf[var]->phi[j];
                      advection = 0.;
                      if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
                        advection -=
                            (s[p][b] * (double)delta(a, q) + s[a][p] * (double)delta(b, q));
                        advection *= phi_j * h3 * det_J;

                        advection *=
                            wt_func * wt * H_ls / Gmod * pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
                      }

                      /*
                       * Source term...
                       */

                      source = 0.;

                      if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
                        source = -phi_j * ((double)delta(a, p) * (double)delta(b, q) +
                                           (double)delta(b, p) * (double)delta(a, q));
                        source *=
                            det_J * h3 * wt_func * wt * pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
                      }

                      lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += advection + source;
                    }
                  }
                }
              }

              /*
               * J_S_S
               */
              for (p = 0; p < VIM; p++) {
                for (q = 0; q < VIM; q++) {
                  var = v_s[mode][p][q];

                  if (pd->v[pg->imtrx][var]) {
                    pvar = upd->vp[pg->imtrx][var];
                    for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                      phi_j = bf[var]->phi[j];
                      mass = 0.;
                      if (pd->TimeIntegration != STEADY) {
                        if (pd->e[pg->imtrx][eqn] & T_MASS) {
                          mass = (1. + 2. * tt) * phi_j / dt * (double)delta(a, p) *
                                 (double)delta(b, q);
                          mass *= h3 * det_J;
                          mass *= wt_func * H_ls / Gmod * wt * pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
                        }
                      }

                      advection = 0.;

                      if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
                        if ((a == p) && (b == q)) {
                          for (r = 0; r < WIM; r++) {
                            advection += (v[r] - x_dot[r]) * bf[var]->grad_phi[j][r];
                          }
                        }

                        advection -= phi_j * (gt[a][p] * (double)delta(b, q) +
                                              g[q][b] * (double)delta(a, p));

                        advection *= h3 * det_J;

                        advection *=
                            wt_func * wt * H_ls / Gmod * pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
                      }

                      /*
                       * Source term...
                       */

                      source = 0.;

                      if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
                        source =
                            phi_j * (double)delta(a, p) * (double)delta(b, q) * (1. - H_ls) / mup;

                        source *=
                            det_J * h3 * wt_func * wt * pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
                      }

                      lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += mass + advection + source;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  return (status);
}

int assemble_rate_of_strain(dbl tt, /* parameter to vary time integration from
                                     * explicit (tt = 1) to implicit (tt = 0) */
                            dbl dt) /* current time step size */
{
  int dim;
  int p, q, a, b;

  int eqn, var;
  int peqn, pvar;
  int i, j;
  int status;

  dbl h3;          /* Volume element (scale factors). */
  dbl dh3dmesh_pj; /* Sensitivity to (p,j) mesh dof. */

  dbl det_J; /* determinant of element Jacobian */

  dbl d_det_J_dmesh_pj; /* for specific (p,j) mesh dof */

  dbl advection;
  dbl advection_a, advection_b;
  dbl source;

  /*
   *
   * Note how carefully we avoid refering to d(phi[i])/dx[j] and refer instead
   * to the j-th component of grad_phi[j][i] so that this vector can be loaded
   * up with components that may be different in non Cartesian coordinate
   * systems.
   *
   * We will, however, insist on *orthogonal* coordinate systems, even if we
   * might permit them to be curvilinear.
   *
   * Assume all components of velocity are interpolated with the same kind
   * of basis function.
   */

  /*
   * Galerkin weighting functions for i-th and a-th momentum residuals
   * and some of their derivatives...
   */

  /*
   * Petrov-Galerkin weighting functions for i-th and ab-th stress residuals
   * and some of their derivatives...
   */

  dbl wt_func;

  /*
   * Interpolation functions for variables and some of their derivatives.
   */

  dbl phi_j;

  dbl wt;

  /* Variables for stress */

  int R_g[DIM][DIM];
  int v_g[DIM][DIM];

  status = 0;

  /*
   * Unpack variables from structures for local convenience...
   */

  dim = pd->Num_Dim;

  eqn = R_GRADIENT11;

  /*
   * Bail out fast if there's nothing to do...
   */

  if (!pd->e[pg->imtrx][eqn]) {
    return (status);
  }

  wt = fv->wt;

  det_J = bf[eqn]->detJ; /* Really, ought to be mesh eqn. */

  h3 = fv->h3; /* Differential volume element (scales). */

  /* load eqn and variable number in tensor form */

  v_g[0][0] = VELOCITY_GRADIENT11;
  v_g[0][1] = VELOCITY_GRADIENT12;
  v_g[1][0] = VELOCITY_GRADIENT21;
  v_g[1][1] = VELOCITY_GRADIENT22;
  v_g[0][2] = VELOCITY_GRADIENT13;
  v_g[1][2] = VELOCITY_GRADIENT23;
  v_g[2][0] = VELOCITY_GRADIENT31;
  v_g[2][1] = VELOCITY_GRADIENT32;
  v_g[2][2] = VELOCITY_GRADIENT33;

  R_g[0][0] = R_GRADIENT11;
  R_g[0][1] = R_GRADIENT12;
  R_g[1][0] = R_GRADIENT21;
  R_g[1][1] = R_GRADIENT22;
  R_g[0][2] = R_GRADIENT13;
  R_g[1][2] = R_GRADIENT23;
  R_g[2][0] = R_GRADIENT31;
  R_g[2][1] = R_GRADIENT32;
  R_g[2][2] = R_GRADIENT33;

  for (a = 0; a < VIM; a++) {
    for (b = 0; b < VIM; b++) {
      if (b < a && pd->e[pg->imtrx][v_g[a][b]]) {
        GOMA_WH(-1, "Rate of Strain expects symmetric gradient elements only");
      }
    }
  }
  /*
   * Field variables...
   */

  dbl gamma[DIM][DIM];
  for (a = 0; a < VIM; a++) {
    for (b = 0; b < VIM; b++) {
      gamma[a][b] = fv->grad_v[a][b] + fv->grad_v[b][a];
    }
  }

  dbl d_gamma_dv[DIM][DIM][DIM][MDE] = {{{{0.0}}}};
  for (a = 0; a < VIM; a++) {
    for (b = 0; b < VIM; b++) {
      for (p = 0; p < VIM; p++) {
        for (int k = 0; k < ei[pg->imtrx]->dof[VELOCITY1 + p]; k++) {
          d_gamma_dv[a][b][p][k] =
              bf[VELOCITY1]->grad_phi_e[k][p][a][b] + bf[VELOCITY1]->grad_phi_e[k][p][b][a];
        }
      }
    }
  }

  dbl d_gamma_dX[DIM][DIM][DIM][MDE] = {{{{0.0}}}};
  for (a = 0; a < VIM; a++) {
    for (b = 0; b < VIM; b++) {
      for (p = 0; p < VIM; p++) {
        for (int k = 0; k < ei[pg->imtrx]->dof[VELOCITY1 + p]; k++) {
          d_gamma_dX[a][b][p][k] = fv->d_grad_v_dmesh[a][b][p][k] + fv->d_grad_v_dmesh[b][a][p][k];
        }
      }
    }
  }

  dbl gammadot;
  dbl d_gd_dv[DIM][MDE];    /* derivative of strain rate invariant
                               wrt velocity */
  dbl d_gd_dmesh[DIM][MDE]; /* derivative of strain rate invariant
                               wrt mesh */
  calc_shearrate(&gammadot, gamma, d_gd_dv, d_gd_dmesh);

  /*
   * Residuals_________________________________________________________________
   */

  if (af->Assemble_Residual) {
    /*
     * Assemble each component "ab" of the velocity gradient equation...
     */
    for (a = 0; a < VIM; a++) {
      for (b = 0; b < VIM; b++) {
        if (b >= a) {
          eqn = R_g[a][b];
          /*
           * In the element, there will be contributions to this many equations
           * based on the number of degrees of freedom...
           */

          for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {

            wt_func = bf[eqn]->phi[i]; /* add Petrov-Galerkin terms as necessary */

            advection = 0.;

            if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
              advection -= gamma[a][b];
              advection *= wt_func * det_J * wt * h3;
              advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
            }

            /*
             * Source term...
             */

            source = 0;

            if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
              dbl Du_eps = sqrt(gammadot + gn->epsilon * gn->epsilon);
              source += Du_eps * fv->G[a][b];
              source *= wt_func * det_J * h3 * wt;
              source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
            }

            lec->R[LEC_R_INDEX(upd->ep[pg->imtrx][eqn], i)] += advection + source;
          }
        }
      }
    }
  }

  /*
   * Jacobian terms...
   */

  if (af->Assemble_Jacobian) {
    for (a = 0; a < VIM; a++) {
      for (b = 0; b < VIM; b++) {
        if (b >= a) {
          eqn = R_g[a][b];
          peqn = upd->ep[pg->imtrx][eqn];

          for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
            wt_func = bf[eqn]->phi[i]; /* add Petrov-Galerkin terms as necessary */

            /*
             * J_G_v
             */
            for (p = 0; p < WIM; p++) {
              var = VELOCITY1 + p;
              if (pd->v[pg->imtrx][var]) {
                pvar = upd->vp[pg->imtrx][var];
                for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                  phi_j = bf[var]->phi[j];

                  advection = 0.;

                  if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
                    advection -= d_gamma_dv[a][b][p][j];
                    advection *= wt_func * det_J * wt * h3;
                    advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
                  }

                  source = 0.;
                  if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
                    dbl Du_eps = sqrt(gammadot + gn->epsilon * gn->epsilon);
                    dbl deriv = d_gd_dv[p][j] * 0.5 * (1.0 / Du_eps);
                    source += deriv * fv->G[a][b];
                    source *= wt_func * det_J * h3 * wt;
                    source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
                  }

                  lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += advection + source;
                }
              }
            }

            /*
             * J_G_d
             */
            for (p = 0; p < dim; p++) {
              var = MESH_DISPLACEMENT1 + p;
              if (pd->v[pg->imtrx][var]) {
                pvar = upd->vp[pg->imtrx][var];
                for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                  phi_j = bf[var]->phi[j];

                  d_det_J_dmesh_pj = bf[eqn]->d_det_J_dm[p][j];

                  dh3dmesh_pj = fv->dh3dq[p] * bf[var]->phi[j];

                  advection = 0.;

                  if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
                    /*
                     * three parts:
                     *    advection_a =
                     *    	Int ( d(Vv)/dmesh h3 |Jv| )
                     *
                     *    advection_b =
                     *  (i)	Int ( Vv h3 d(|Jv|)/dmesh )
                     *  (ii)      Int ( Vv dh3/dmesh |Jv|   )
                     */

                    advection_a = -gamma[a][b];

                    advection_a *= (d_det_J_dmesh_pj * h3 + det_J * dh3dmesh_pj);

                    advection_b = -d_gamma_dX[a][b][p][j];

                    advection_b *= det_J * h3;

                    advection = advection_a + advection_b;

                    advection *= wt_func * wt * pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
                  }

                  source = 0.;

                  if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
                    dbl Du_eps = sqrt(gammadot + gn->epsilon * gn->epsilon);
                    dbl deriv = d_gd_dv[p][j] * 0.5 * (1.0 / Du_eps);
                    dbl source_a = 0.0;
                    source_a += deriv * fv->G[a][b];
                    source_a *= det_J * h3 * wt;
                    dbl source_b = 0.0;
                    source_b += fv->G[a][b];
                    source_b *= d_det_J_dmesh_pj * h3 + det_J * dh3dmesh_pj;
                    source = source_a + source_b;
                    source *= wt_func * wt * pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
                  }

                  lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += advection + source;
                }
              }
            }

            /*
             * J_G_G
             */

            for (p = 0; p < VIM; p++) {
              for (q = 0; q < VIM; q++) {
                var = v_g[p][q];

                if (pd->v[pg->imtrx][var]) {
                  pvar = upd->vp[pg->imtrx][var];
                  for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                    phi_j = bf[var]->phi[j];

                    source = 0.;

                    if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
                      if ((a == p) && (b == q)) {
                        dbl Du_eps = sqrt(gammadot + gn->epsilon * gn->epsilon);
                        source = Du_eps * phi_j * det_J * h3 * wt_func * wt *
                                 pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
                      }
                    }

                    lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += source;
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  return (status);
}

int assemble_gradient(dbl tt, /* parameter to vary time integration from
                               * explicit (tt = 1) to implicit (tt = 0) */
                      dbl dt) /* current time step size */
{
  int dim;
  int p, q, a, b;

  int eqn, var;
  int peqn, pvar;
  int i, j;
  int status;

  dbl h3;          /* Volume element (scale factors). */
  dbl dh3dmesh_pj; /* Sensitivity to (p,j) mesh dof. */

  dbl grad_v[DIM][DIM];
  dbl g[DIM][DIM]; /* velocity gradient tensor */

  dbl det_J; /* determinant of element Jacobian */

  dbl d_det_J_dmesh_pj; /* for specific (p,j) mesh dof */

  dbl advection;
  dbl advection_a, advection_b;
  dbl source;

  /*
   *
   * Note how carefully we avoid refering to d(phi[i])/dx[j] and refer instead
   * to the j-th component of grad_phi[j][i] so that this vector can be loaded
   * up with components that may be different in non Cartesian coordinate
   * systems.
   *
   * We will, however, insist on *orthogonal* coordinate systems, even if we
   * might permit them to be curvilinear.
   *
   * Assume all components of velocity are interpolated with the same kind
   * of basis function.
   */

  /*
   * Galerkin weighting functions for i-th and a-th momentum residuals
   * and some of their derivatives...
   */

  /*
   * Petrov-Galerkin weighting functions for i-th and ab-th stress residuals
   * and some of their derivatives...
   */

  dbl wt_func;

  /*
   * Interpolation functions for variables and some of their derivatives.
   */

  dbl phi_j;

  dbl wt;

  /* Variables for stress */

  int R_g[DIM][DIM];
  int v_g[DIM][DIM];

  status = 0;

  /*
   * Unpack variables from structures for local convenience...
   */

  dim = pd->Num_Dim;

  eqn = R_GRADIENT11;

  /*
   * Bail out fast if there's nothing to do...
   */

  if (!pd->e[pg->imtrx][eqn]) {
    return (status);
  }

  wt = fv->wt;

  det_J = bf[eqn]->detJ; /* Really, ought to be mesh eqn. */

  h3 = fv->h3; /* Differential volume element (scales). */

  /* load eqn and variable number in tensor form */

  v_g[0][0] = VELOCITY_GRADIENT11;
  v_g[0][1] = VELOCITY_GRADIENT12;
  v_g[1][0] = VELOCITY_GRADIENT21;
  v_g[1][1] = VELOCITY_GRADIENT22;
  v_g[0][2] = VELOCITY_GRADIENT13;
  v_g[1][2] = VELOCITY_GRADIENT23;
  v_g[2][0] = VELOCITY_GRADIENT31;
  v_g[2][1] = VELOCITY_GRADIENT32;
  v_g[2][2] = VELOCITY_GRADIENT33;

  R_g[0][0] = R_GRADIENT11;
  R_g[0][1] = R_GRADIENT12;
  R_g[1][0] = R_GRADIENT21;
  R_g[1][1] = R_GRADIENT22;
  R_g[0][2] = R_GRADIENT13;
  R_g[1][2] = R_GRADIENT23;
  R_g[2][0] = R_GRADIENT31;
  R_g[2][1] = R_GRADIENT32;
  R_g[2][2] = R_GRADIENT33;

  /*
   * Field variables...
   */

  /*
   * In Cartesian coordinates, this velocity gradient tensor will
   * have components that are...
   *
   * 			grad_v[a][b] = d v_b
   *				       -----
   *				       d x_a
   */

  for (a = 0; a < VIM; a++) {
    for (b = 0; b < VIM; b++) {
      grad_v[a][b] = fv->grad_v[a][b];
    }
  }

  for (a = 0; a < VIM; a++) {
    for (b = 0; b < VIM; b++) {
      g[a][b] = fv->G[a][b];
    }
  }

  /*
   * Residuals_________________________________________________________________
   */

  if (af->Assemble_Residual) {
    /*
     * Assemble each component "ab" of the velocity gradient equation...
     */
    for (a = 0; a < VIM; a++) {
      for (b = 0; b < VIM; b++) {
        eqn = R_g[a][b];
        /*
         * In the element, there will be contributions to this many equations
         * based on the number of degrees of freedom...
         */

        for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {

          wt_func = bf[eqn]->phi[i]; /* add Petrov-Galerkin terms as necessary */

          advection = 0.;

          if (upd->devss_traceless_gradient) {
            if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
              advection -= grad_v[a][b] - fv->div_v * delta(a, b) / ((dbl)VIM);
              advection *= wt_func * det_J * wt * h3;
              advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
            }
          } else {
            if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
              advection -= grad_v[a][b];
              advection *= wt_func * det_J * wt * h3;
              advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
            }
          }

          /*
           * Source term...
           */

          source = 0;

          if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
            source += g[a][b];
            source *= wt_func * det_J * h3 * wt;
            source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
          }

          lec->R[LEC_R_INDEX(upd->ep[pg->imtrx][eqn], i)] += advection + source;
        }
      }
    }
  }

  /*
   * Jacobian terms...
   */

  if (af->Assemble_Jacobian) {
    for (a = 0; a < VIM; a++) {
      for (b = 0; b < VIM; b++) {
        eqn = R_g[a][b];
        peqn = upd->ep[pg->imtrx][eqn];

        for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
          wt_func = bf[eqn]->phi[i]; /* add Petrov-Galerkin terms as necessary */

          /*
           * J_G_v
           */
          for (p = 0; p < WIM; p++) {
            var = VELOCITY1 + p;
            if (pd->v[pg->imtrx][var]) {
              pvar = upd->vp[pg->imtrx][var];
              for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                phi_j = bf[var]->phi[j];

                advection = 0.;

                if (upd->devss_traceless_gradient) {
                  dbl div_phi_j_e_p = 0.;
                  for (int b = 0; b < VIM; b++) {
                    div_phi_j_e_p += bf[var]->grad_phi_e[j][p][b][b];
                  }
                  if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
                    advection -=
                        bf[var]->grad_phi_e[j][p][a][b] - div_phi_j_e_p * delta(a, b) / ((dbl)VIM);
                    advection *= wt_func * det_J * wt * h3;
                    advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
                  }
                } else {
                  if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
                    advection -= bf[var]->grad_phi_e[j][p][a][b];
                    advection *= wt_func * det_J * wt * h3;
                    advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
                  }
                }

                source = 0.;

                lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += advection + source;
              }
            }
          }

          /*
           * J_G_d
           */
          for (p = 0; p < dim; p++) {
            var = MESH_DISPLACEMENT1 + p;
            if (pd->v[pg->imtrx][var]) {
              pvar = upd->vp[pg->imtrx][var];
              for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                phi_j = bf[var]->phi[j];

                d_det_J_dmesh_pj = bf[eqn]->d_det_J_dm[p][j];

                dh3dmesh_pj = fv->dh3dq[p] * bf[var]->phi[j];

                advection = 0.;

                if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
                  /*
                   * three parts:
                   *    advection_a =
                   *    	Int ( d(Vv)/dmesh h3 |Jv| )
                   *
                   *    advection_b =
                   *  (i)	Int ( Vv h3 d(|Jv|)/dmesh )
                   *  (ii)      Int ( Vv dh3/dmesh |Jv|   )
                   */

                  advection_a = -grad_v[a][b];

                  advection_a *= (d_det_J_dmesh_pj * h3 + det_J * dh3dmesh_pj);

                  advection_b = -fv->d_grad_v_dmesh[a][b][p][j];

                  advection_b *= det_J * h3;

                  advection = advection_a + advection_b;

                  advection *= wt_func * wt * pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
                }

                source = 0.;

                if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
                  source += g[a][b];

                  source *= d_det_J_dmesh_pj * h3 + det_J * dh3dmesh_pj;

                  source *= wt_func * wt * pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
                }

                lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += advection + source;
              }
            }
          }

          /*
           * J_G_G
           */

          for (p = 0; p < VIM; p++) {
            for (q = 0; q < VIM; q++) {
              var = v_g[p][q];

              if (pd->v[pg->imtrx][var]) {
                pvar = upd->vp[pg->imtrx][var];
                for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                  phi_j = bf[var]->phi[j];

                  source = 0.;

                  if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
                    if ((a == p) && (b == q)) {
                      source = phi_j * det_J * h3 * wt_func * wt *
                               pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
                    }
                  }

                  lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += source;
                }
              }
            }
          }
        }
      }
    }
  }
  return (status);
}

int tensor_dot(dbl t1[DIM][DIM], dbl t2[DIM][DIM], dbl t1_dot_t2[DIM][DIM], const int dim) {
  int i, j, k;
  int status;
  dbl v1[DIM];
  dbl v2[DIM];

  for (k = 0; k < dim; k++) {
    for (i = 0; i < dim; i++) {
      v1[i] = t1[k][i];
    }
    for (j = 0; j < dim; j++) {
      for (i = 0; i < dim; i++) {
        v2[i] = t2[i][j];
      }
      t1_dot_t2[k][j] = vec_dot(dim, v1, v2);
    }
  }

  status = 1;
  return (status);
}

dbl vec_dot(const int n1, dbl *v1, dbl *v2) {
  int i;
  dbl rc = 0.0;

  for (i = 0; i < n1; i++) {
    rc += *v1 * *v2;
    v1++;
    v2++;
  }
  return (rc);
}

void load_modal_pointers(int ve_mode, /* mode number */
                         dbl tt,
                         dbl dt,
                         dbl s[DIM][DIM],     /* stress tensor for mode ve_mode */
                         dbl s_dot[DIM][DIM], /* stress tensor time derivative for mode ve_mode */
                         dbl grad_s[DIM][DIM][DIM], /* grad of stress tensor for mode ve_mode */
                         dbl d_grad_s_dm[DIM][DIM][DIM][DIM][MDE]) /* derivative of grad of stress
                                                                      tensor for mode ve_mode */

{
  int a, b, p, q; /* indeces for dimensions */
  int j;          /* indeces for dofs */
  int var;
  /* load up things we need in the assembly routine for each mode in turn*/

  /* put stress in a nice working array */

  for (a = 0; a < VIM; a++) {
    for (b = 0; b < VIM; b++) {
      s[a][b] = fv->S[ve_mode][a][b];
      if (pd->TimeIntegration != STEADY) {
        s_dot[a][b] = fv_dot->S[ve_mode][a][b];
      } else {
        s_dot[a][b] = 0.;
      }
    }
  }

  for (p = 0; p < VIM; p++) {
    for (a = 0; a < VIM; a++) {
      for (b = 0; b < VIM; b++) {
        grad_s[p][a][b] = fv->grad_S[ve_mode][p][a][b];
      }
    }
  }

  var = MESH_DISPLACEMENT1;
  if (pd->v[pg->imtrx][var]) {
    for (p = 0; p < VIM; p++) {
      for (a = 0; a < VIM; a++) {
        for (b = 0; b < VIM; b++) {
          for (q = 0; q < ei[pg->imtrx]->ielem_dim; q++) {
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              d_grad_s_dm[p][a][b][q][j] = fv->d_grad_S_dmesh[ve_mode][p][a][b][q][j];
            }
          }
        }
      }
    }
  }
}

/******************************************************************************/
/* END of routine modal_esp_alloc */
/******************************************************************************/

int assemble_surface_stress(Exo_DB *exo, /* ptr to basic exodus ii mesh information */
                            double x[],
                            struct GomaLinearSolverData *ams,
                            dbl x_update[],      /* last update for x vector */
                            double delta_t,      /* current time step size */
                            double t_,           /* parameter to vary time integration from
                                                  * explicit (tt = 1) to implicit (tt = 0) */
                            int ielem_type,      /* element type  */
                            int ielem_type_fill, /* element type for fill function */
                            int id_side,         /* id number of current side according to
                                                  * EXODUS convention  */
                            int neighbor,        /* element neighboring this side */
                            int ielem,           /* current element */
                            int num_local_nodes) /* number of nodes per element */
{
  /*    TAB certifies that this function conforms to the exo/patran side numbering convention
   * 11/9/98. */

  /* LOCAL VARIABLES */
  int ip, ip1, i, j, dim; /* counters */
  int a, b, p, q;         /* more counters */
  int nodes_per_side;
  int local_elem_node_id[MAX_NODES_PER_SIDE];

  int eqn, peqn;
  int var, pvar;
  int err; /* status variable for functions */
  int ip_total, ip_total_fill;
  int found_it;

  /* Variables for stress */
  int R_s[MAX_MODES][DIM][DIM];
  int v_s[MAX_MODES][DIM][DIM];
  int S_map[MAX_MODES][DIM][DIM]; /* map var index to stress mode component */

  int ibc, ins, side_index;
  int table_ibc[MAX_MODES][DIM][DIM]; /* maps table boundary condition index for each stress mode */
  int num_zeros;

  int doit = 1;

  int mode;

  dbl s[DIM][DIM], s_dot[DIM][DIM], grad_s[DIM][DIM][DIM];
  dbl stress_neighbor[MAX_SURF_GP][MAX_MODES][DIM][DIM];
  dbl stress_update_v[MAX_SURF_GP][MAX_MODES][DIM][DIM];
  dbl s_n[MAX_MODES][DIM][DIM];
  dbl d_grad_s_dmesh[DIM][DIM][DIM][DIM][MDE];
  dbl **x_neighbor;

  dbl *J_S_S_v[MAX_MODES][DIM][DIM][MDE][DIM][DIM][MDE];
  dbl phi_neighbor[MAX_SURF_GP][MDE];

  dbl advection;

  double phi_j, phi_i;
  dbl rhs;
  double ss, tt, uu; /* Gaussian quadrature point locations  */
  double xi[DIM];    /* Local element coordinates of Gauss point. */
  dbl vdotn, vdotn_avg = 0;
  dbl vdotn_norm = 0;
  double wt; /* Quadrature weights units - ergs/(sec*cm*K) = g*cm/(sec^3*K)     */

  dbl *phi_v = NULL;

  dbl alpha = 0.5;

  /***************************************************************************/

  /* load eqn and variable number in tensor form */
  err = stress_eqn_pointer(S_map);

  err = stress_eqn_pointer(v_s);
  err = stress_eqn_pointer(R_s);

  /* initialization of the neighbor stress pointers array */

  memset(J_S_S_v, 0, MAX_MODES * DIM * DIM * MDE * DIM * DIM * MDE);

  /********************************************************************************/
  /*     START OF SURFACE LOOPS THAT REQUIRE INTEGRATION (WEAK SENSE)             */
  /*                AND REQUIRE ROTATION IN TO N-T FORM                           */
  /********************************************************************************/
  /* Find out the number of surface quadrature points
     -this is assumed independent of the surface */
  ip_total = elem_info(NQUAD_SURF, ielem_type);

  dim = pd->Num_Dim;

  /* allocate space for x_neighbor */

  /*  x_neighbor = (double **) array_alloc(2, ip_total, DIM, sizeof(double)); */
  /*  manually allocate space to correct (seeming) misalignment for HPUX */

  x_neighbor = (double **)smalloc(ip_total * (int)sizeof(double *));
  for (i = 0; i < DIM; i++) {
    x_neighbor[i] = (double *)smalloc(DIM * sizeof(double));
  }

  /* If no neighbor element found, check for table boundary condition on inlet */

  for (mode = 0; mode < vn->modes; mode++) {
    for (a = 0; a < VIM; a++) {
      for (b = 0; b < VIM; b++) {
        table_ibc[mode][a][b] = -1;
      }
    }
  }

  for (ibc = 0; ibc < Num_BC; ibc++) {
    if (BC_Types[ibc].BC_Name == TABLE_BC) {
      /*Loop over all side sets to find a match */
      for (ins = 0; ins < exo->num_side_sets; ins++) {
        if (Proc_SS_Ids[ins] == BC_Types[ibc].BC_ID) {

          /* Does it contain the element? */
          for (side_index = 0; side_index < exo->ss_num_sides[ins]; side_index++) {

            if (ielem == exo->ss_elem_list[exo->ss_elem_index[ins] + side_index]) {
              /* which variable is the table for? */
              switch (BC_Types[ibc].table->f_index) {
              case POLYMER_STRESS11: {
                table_ibc[0][0][0] = ibc;
                break;
              }
              case POLYMER_STRESS12: {
                table_ibc[0][0][1] = ibc;
                break;
              }
              case POLYMER_STRESS22: {
                table_ibc[0][1][1] = ibc;
                break;
              }
              case POLYMER_STRESS11_1: {
                table_ibc[1][0][0] = ibc;
                break;
              }
              case POLYMER_STRESS12_1: {
                table_ibc[1][0][1] = ibc;
                break;
              }
              case POLYMER_STRESS22_1: {
                table_ibc[1][1][1] = ibc;
                break;
              }
              case POLYMER_STRESS11_2: {
                table_ibc[2][0][0] = ibc;
                break;
              }
              case POLYMER_STRESS12_2: {
                table_ibc[2][0][1] = ibc;
                break;
              }
              case POLYMER_STRESS22_2: {
                table_ibc[2][1][1] = ibc;
                break;
              }
              case POLYMER_STRESS11_3: {
                table_ibc[3][0][0] = ibc;
                break;
              }
              case POLYMER_STRESS12_3: {
                table_ibc[3][0][1] = ibc;
                break;
              }
              case POLYMER_STRESS22_3: {
                table_ibc[3][1][1] = ibc;
                break;
              }
              case POLYMER_STRESS11_4: {
                table_ibc[4][0][0] = ibc;
                break;
              }
              case POLYMER_STRESS12_4: {
                table_ibc[4][0][1] = ibc;
                break;
              }
              case POLYMER_STRESS22_4: {
                table_ibc[4][1][1] = ibc;
                break;
              }
              case POLYMER_STRESS11_5: {
                table_ibc[5][0][0] = ibc;
                break;
              }
              case POLYMER_STRESS12_5: {
                table_ibc[5][0][1] = ibc;
                break;
              }
              case POLYMER_STRESS22_5: {
                table_ibc[5][1][1] = ibc;
                break;
              }
              case POLYMER_STRESS11_6: {
                table_ibc[6][0][0] = ibc;
                break;
              }
              case POLYMER_STRESS12_6: {
                table_ibc[6][0][1] = ibc;
                break;
              }
              case POLYMER_STRESS22_6: {
                table_ibc[6][1][1] = ibc;
                break;
              }
              case POLYMER_STRESS11_7: {
                table_ibc[7][0][0] = ibc;
                break;
              }
              case POLYMER_STRESS12_7: {
                table_ibc[7][0][1] = ibc;
                break;
              }
              case POLYMER_STRESS22_7: {
                table_ibc[7][1][1] = ibc;
                break;
              }
              default: {
                break;
              }
              }
            } /* end check for right element */
          }   /* end loop over sides of ins */
        }     /*close if loop for BC_ID check */
      }       /*close loop over all side sets */
    }         /* end if loop over tables */
  }           /* end loop over all bc's */

  num_zeros = 0;

  for (mode = 0; mode < vn->modes; mode++) {
    for (a = 0; a < VIM; a++) {
      for (b = 0; b < VIM; b++) {
        if ((a <= b) && (table_ibc[mode][a][b] == -1)) {
          num_zeros++;
        }
      }
    }
  }

  /* if num_zeros == 0, then table bc's exist for all stress modes
     if num_zeros == vn->modes*(VIM!), then no table bc's exist
     if num_zeros in between - there is a problem...  */

  /* Use one point surface quadrature integration
     to get the sign of v*n */

  /*   ip_total_fill = elem_info(NQUAD_SURF, ielem_type_fill); */

  ip_total_fill = 1;

  /* Surface integration over element */

  for (ip = 0; ip < ip_total_fill; ip++) {
    /* find the quadrature point locations for current ip */

    find_surf_st(ip, P1_QUAD, id_side, dim, xi, &ss, &tt, &uu);

    /* find the quadrature weight for current ip */
    wt = Gq_surf_weight(ip, ielem_type_fill);

    /* ****************************************/
    err = load_basis_functions(xi, bfd);
    GOMA_EH(err, "problem from load_basis_functions");

    err = beer_belly();
    GOMA_EH(err, "beer_belly");

    /* precalculate variables at  current integration pt.*/

    err = load_fv();
    GOMA_EH(err, "load_fv");

    /* calculate the determinant of the surface jacobian and the normal to
     * the surface all at one time
     */

    err = get_side_info(ielem_type, id_side, &nodes_per_side, local_elem_node_id);
    GOMA_EH(err, "get_side_info");

    surface_determinant_and_normal(ielem, ei[pg->imtrx]->iconnect_ptr, num_local_nodes,
                                   ei[pg->imtrx]->ielem_dim - 1, id_side, nodes_per_side,
                                   local_elem_node_id);

    do_LSA_mods(LSA_SURFACE);

    vdotn_avg = 0.;
    vdotn_norm = 0.;
    for (a = 0; a < WIM; a++) {
      vdotn_avg += fv->v[a] * fv->snormal[a];
      vdotn_norm += fv->v[a] * fv->v[a];
    }

    if (vdotn_avg < 0. && vdotn_avg * vdotn_avg / vdotn_norm > 1.e-12) {
      if (neighbor != -1) {
        err = neighbor_stress(exo, x, x_update, ielem, neighbor, stress_neighbor, stress_update_v,
                              phi_neighbor, num_local_nodes, nodes_per_side, local_elem_node_id,
                              ielem_type, ielem_type_fill, x_neighbor, S_map);
        GOMA_EH(err, "neighbor_stress");
      } else if ((neighbor == -1) && (num_zeros == 0)) {
        /* inlet table boundary consitions exist for the stress components */

        err =
            neighbor_stress_table(exo, x, x_update, ielem, stress_neighbor, stress_update_v,
                                  phi_neighbor, num_local_nodes, nodes_per_side, local_elem_node_id,
                                  ielem_type, ielem_type_fill, x_neighbor, S_map, table_ibc);
        GOMA_EH(err, "neighbor_stress_table");
      } else {
        /* if there is no neighbor and no tables, set this value to zero
         * I am assuming we will only get here for inflow
         * boundaries... later we will do something better.*/
        for (a = 0; a < MAX_SURF_GP; a++) {
          for (b = 0; b < MDE; b++) {
            phi_neighbor[a][b] = 0.0;
          }
        }

        phi_v = phi_neighbor[0];

        for (mode = 0; mode < vn->modes; mode++) {
          for (a = 0; a < VIM; a++) {
            for (b = 0; b < VIM; b++) {
              s_n[mode][a][b] = 0.;
            }
          }
        }
      }
    }
  }

  /* Surface integration over element */
  if (vdotn_avg < 0. && vdotn_avg * vdotn_avg / vdotn_norm > 1.e-12 &&
      ((neighbor != -1) || (num_zeros == 0))) {
    for (ip = 0; ip < ip_total; ip++) {
      /* find the quadrature point locations for current ip */

      find_surf_st(ip, ielem_type, id_side, dim, xi, &ss, &tt, &uu);

      /* find the quadrature weight for current ip */
      wt = Gq_surf_weight(ip, ielem_type);

      /* ****************************************/
      err = load_basis_functions(xi, bfd);
      GOMA_EH(err, "problem from load_basis_functions");

      err = beer_belly();
      GOMA_EH(err, "beer_belly");

      /* precalculate variables at  current integration pt.*/
      err = load_fv();
      GOMA_EH(err, "load_fv");

      /* calculate the determinant of the surface jacobian and the normal to
       * the surface all at one time */

      err = get_side_info(ielem_type, id_side, &nodes_per_side, local_elem_node_id);
      GOMA_EH(err, "get_side_info");

      surface_determinant_and_normal(ielem, ei[pg->imtrx]->iconnect_ptr, num_local_nodes,
                                     ei[pg->imtrx]->ielem_dim - 1, id_side, nodes_per_side,
                                     local_elem_node_id);

      do_LSA_mods(LSA_SURFACE);

      if ((neighbor != -1) || (num_zeros == 0)) {
        found_it = 0;
        for (ip1 = 0; ip1 < ip_total && (!found_it); ip1++) {
          if ((fabs(fv->x0[0] - x_neighbor[ip1][0]) < 1.e-7) &&
              (fabs(fv->x0[1] - x_neighbor[ip1][1]) < 1.e-7) &&
              (fabs(fv->x0[2] - x_neighbor[ip1][2]) < 1.e-7)) {
            found_it = 1;
            phi_v = phi_neighbor[ip1];

            for (mode = 0; mode < vn->modes; mode++) {
              for (a = 0; a < VIM; a++) {
                for (b = 0; b < VIM; b++) {
                  /* since the stress tensor is symmetric, only assemble the upper half */
                  if (a <= b) {
                    s_n[mode][a][b] = stress_neighbor[ip1][mode][a][b];
                  }
                }
              }
            }
          }
        }
      }

      vdotn = 0.;
      for (a = 0; a < WIM; a++) {
        vdotn += fv->v[a] * fv->snormal[a];
      }

      /* add alpha for upwinding and possible stability */
      vdotn *= alpha;

      vdotn_norm = sqrt(vdotn * vdotn);

      if ((vdotn < 0.) && (vdotn_norm > 1.e-7)) {
        /*
         * Put local contributions into global right-hand side
         * if it is not a right-hand side variable-it won't get added in (contribution
         * is zero)
         */
        for (mode = 0; mode < vn->modes; mode++) {
          load_modal_pointers(mode, t_, delta_t, s, s_dot, grad_s, d_grad_s_dmesh);

          if (vn->dg_J_model == FULL_DG && Linear_Solver != FRONT) {
            load_neighbor_pointers(exo, ams, neighbor, ielem_type, mode, R_s, v_s, J_S_S_v[mode]);
          }

          if (af->Assemble_Residual) {
            for (a = 0; a < VIM; a++) {
              for (b = 0; b < VIM; b++) {
                /* since the stress tensor is symmetric, only assemble the upper half */
                if (a <= b) {
                  eqn = R_s[mode][a][b];
                  if (pd->e[pg->imtrx][eqn]) {
                    peqn = upd->ep[pg->imtrx][eqn];

                    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
                      phi_i = bf[eqn]->phi[i];

                      rhs = phi_i * wt * fv->sdet * ve[mode]->time_const_st->lambda0 * vdotn *
                            (s[a][b] - s_n[mode][a][b]);

                      lec->R[LEC_R_INDEX(peqn, i)] -= rhs;
                    }
                  }
                }
              }
            }
          }

          if (af->Assemble_Jacobian && doit) {
            for (a = 0; a < VIM; a++) {
              for (b = 0; b < VIM; b++) {
                /* since the stress tensor is symmetric, only assemble the upper half */
                if (a <= b) {
                  eqn = R_s[mode][a][b];
                  if (pd->e[pg->imtrx][eqn]) {
                    peqn = upd->ep[pg->imtrx][eqn];

                    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
                      phi_i = bf[eqn]->phi[i];

                      /*
                       * J_S_S
                       */
                      var = eqn;
                      pvar = upd->vp[pg->imtrx][var];
                      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                        phi_j = bf[var]->phi[j];

                        advection = wt * fv->sdet * vdotn * ve[mode]->time_const_st->lambda0 *
                                    phi_j * phi_i;

                        /* Work better without this?????, see PRS concern above
                           Or is this a correction ???*/
                        lec->J[LEC_J_INDEX(peqn, pvar, i, j)] -= advection;

                        advection = 0;
                        if (vn->dg_J_model == FULL_DG && found_it) {
                          advection = wt * fv->sdet * vdotn * ve[mode]->time_const_st->lambda0 *
                                      phi_v[j] * phi_i;

                          if (Linear_Solver != FRONT) {
                            *J_S_S_v[mode][a][b][i][a][b][j] += advection;
                          } else {

                            /* Notice how this is diagonal only */
                            /* i.e., T_12_i is only depending on T_12_j on other face element*/
                            lec->J_stress_neighbor[LEC_J_STRESS_INDEX(id_side - 1, i, peqn, j)] +=
                                advection;
                          }
                        }
                      }

                      /*
                       * J_S_v  sensitivity of stress equation w.r.t. velocity
                       */

                      for (p = 0; p < dim; p++) {
                        var = VELOCITY1 + p;
                        if (pd->v[pg->imtrx][var]) {
                          pvar = upd->vp[pg->imtrx][var];
                          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                            phi_j = bf[var]->phi[j];

                            advection = phi_i * wt * fv->sdet * ve[mode]->time_const_st->lambda0 *
                                        phi_j * fv->snormal[p] * alpha *
                                        (s[a][b] - s_n[mode][a][b]);

                            lec->J[LEC_J_INDEX(peqn, pvar, i, j)] -= advection;
                          }
                        }
                      }

                      /*
                       * J_S_d  sensitivity of stress equation w.r.t. mesh displacement
                       */

                      for (p = 0; p < dim; p++) {
                        var = MESH_DISPLACEMENT1 + p;
                        if (pd->v[pg->imtrx][var]) {
                          pvar = upd->vp[pg->imtrx][var];
                          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                            phi_j = bf[var]->phi[j];

                            advection = 0.;
                            for (q = 0; q < dim; q++) {
                              advection += fv->sdet * fv->v[q] * fv->dsnormal_dx[q][p][j];
                            }

                            advection += fv->dsurfdet_dx[p][j] * vdotn;

                            advection *= phi_i * wt * ve[mode]->time_const_st->lambda0 *
                                         (s[a][b] - s_n[mode][a][b]);

                            lec->J[LEC_J_INDEX(peqn, pvar, i, j)] -= advection;
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        } /* end mode loop */
      }
    }
  }

  for (i = 0; i < ip_total; i++) {
    safer_free((void **)&x_neighbor[i]);
  }

  safer_free((void **)&x_neighbor);

  return 0;
}
/*****************************************************************************/
/* END of routine assemble_surface_stress */
/*****************************************************************************/

int neighbor_stress(Exo_DB *exo, /* ptr to basic exodus ii mesh information */
                    dbl x[],
                    dbl x_update[],
                    int current_elem,
                    int neighbor_elem,
                    dbl stress_neighbor[][MAX_MODES][DIM][DIM],
                    dbl snv[][MAX_MODES][DIM][DIM],
                    dbl phi_v[][MDE],
                    int num_local_nodes,
                    int nodes_per_side,
                    int local_elem_node_id[],
                    int ielem_type,
                    int ielem_type_fill,
                    dbl **x_n,
                    int v_s[MAX_MODES][DIM][DIM])
/*
 *   This function take the current element and side and
 *   knowing who the neighboring element is, finds the value
 *   of the fill function at the current gauss point location
 *   in the neighboring element.
 *

 *  TAB certifies that this function conforms to the exo/patran side numbering convention 11/9/98.
 */

{
  int gnn, ledof, i, j, mode;
  int iconnect_ptr;
  int id_side, id;
  int id_local_elem_coord[MAX_NODES_PER_SIDE];
  int ie;
  int ielem_shape;
  int inode[MAX_NODES_PER_SIDE];
  int ip, ip_total;
  int a, b, p; /* counters */
  const int dim = pd->Num_Dim;
  int nvdof;
  int status = 0;
  int v;

  dbl phi[MDE], phi_map[MDE], arg_j, s, t, u;
  dbl xi[DIM];

  ielem_shape = type2shape(ielem_type);
  ip_total = elem_info(NQUAD_SURF, ielem_type);

  for (ip = 0; ip < ip_total; ip++) {
    for (p = 0; p < DIM; p++) {
      x_n[ip][p] = 0.;
    }

    for (j = 0; j < MDE; j++) {
      phi_v[ip][j] = 0.;
    }

    for (mode = 0; mode < vn->modes; mode++) {
      for (a = 0; a < VIM; a++) {
        for (b = 0; b < VIM; b++) {
          {
            stress_neighbor[ip][mode][a][b] = 0.;
            snv[ip][mode][a][b] = 0.;
          }
        }
      }
    }
  }

  /* first get global node numbers for current element
     and side */
  iconnect_ptr = Proc_Connect_Ptr[current_elem]; /* find pointer to beginning */

  for (i = 0; i < nodes_per_side; i++) {
    id = local_elem_node_id[i];
    /* load up global node numbers on this side */
    inode[i] = Proc_Elem_Connect[iconnect_ptr + id];
  }

  /* find localside number for neighbor element from
     global node numbers of current element */

  id_side = find_id_side(neighbor_elem, nodes_per_side, inode, id_local_elem_coord, exo);

  for (ip = 0; ip < ip_total; ip++) {

    /* find the quadrature point locations for current ip */
    find_surf_st(ip, ielem_type, id_side, pd->Num_Dim, xi, &s, &t, &u);

    /*
     *  we are cheating here and hoping that the "ei[pg->imtrx]->dof[v]" for
     *  the current element will work on the neighbor element that
     *  we are trying to get information for
     */
    v = POLYMER_STRESS11;

    /* first load phi for the fill function */
    for (i = 0; i < ei[pg->imtrx]->dof[v]; i++) {
      phi[i] = newshape(xi, ielem_type, PSI, ei[pg->imtrx]->dof_list[v][i], ielem_shape,
                        pd->i[pg->imtrx][v], i);

      phi_v[ip][i] = phi[i];
    }

    v = pd->ShapeVar;
    /*
     *  we are cheating here and hoping that the "ei[pg->imtrx]->dof[v]" for
     *  the current element will work on the neighbor element that
     *  we are trying to get information for
     */

    iconnect_ptr = Proc_Connect_Ptr[neighbor_elem]; /* find pointer to beginning */
    for (i = 0; i < ei[pg->imtrx]->dof[v]; i++) {
      phi_map[i] = newshape(xi, ielem_type, PSI, ei[pg->imtrx]->dof_list[v][i], ielem_shape,
                            pd->i[pg->imtrx][v], i);
    }

    iconnect_ptr = Proc_Connect_Ptr[neighbor_elem]; /* find pointer to beginning */
    for (i = 0; i < ei[pg->imtrx]->dof[v]; i++) {
      gnn = Proc_Elem_Connect[iconnect_ptr + i];

      for (p = 0; p < dim; p++) {
        x_n[ip][p] += Coor[p][gnn] * phi_map[i];
      }
    }

    for (mode = 0; mode < vn->modes; mode++) {
      for (a = 0; a < VIM; a++) {
        for (b = 0; b < VIM; b++) {
          if (a <= b) {
            v = v_s[mode][a][b];
            if (pd->v[pg->imtrx][v]) {
              for (i = 0; i < num_local_nodes; i++) {
                gnn = Proc_Elem_Connect[iconnect_ptr + i];
                nvdof = Dolphin[pg->imtrx][gnn][v];
                for (j = 0; j < nvdof; j++) {
                  ledof = ei[pg->imtrx]->lvdof_to_ledof[v][j];
                  ie = Index_Solution(gnn, v, 0, j, ei[pg->imtrx]->matID_ledof[ledof], pg->imtrx);
                  GOMA_EH(ie, "Could not find vbl in sparse matrix.");
                  if (vn->dg_J_model == EXPLICIT_DG || vn->dg_J_model == SEGREGATED) {
                    arg_j = x[ie] - vn->dg_J_model_wt[0] * x_update[ie];
                  } else {
                    arg_j = x[ie];
                  }
                  stress_neighbor[ip][mode][a][b] += arg_j * phi[j];
                }
              }
            }
          }
        }
      }
    }
  }

  status = 1;
  return (status);
}
/*****************************************************************************/
/* END of routine neighbor_stress */
/*****************************************************************************/

int neighbor_stress_table(Exo_DB *exo, /* ptr to basic exodus ii mesh information */
                          dbl x[],
                          dbl x_update[],
                          int current_elem,
                          dbl stress_neighbor[][MAX_MODES][DIM][DIM],
                          dbl snv[][MAX_MODES][DIM][DIM],
                          dbl phi_v[][MDE],
                          int num_local_nodes,
                          int nodes_per_side,
                          int local_elem_node_id[],
                          int ielem_type,
                          int ielem_type_fill,
                          dbl **x_n,
                          int v_s[MAX_MODES][DIM][DIM],
                          int table_ibc[MAX_MODES][DIM][DIM])
/*
 *   This function take the current element and side and
 *   knowing where table boundary conditions are stored,
 *   finds the value of the fill function at the current
 *   gauss point location in the current element.
 *
 */

{
  int gnn, i, j, mode;
  int iconnect_ptr;
  int id_side, id;
  int id_local_elem_coord[MAX_NODES_PER_SIDE];
  int ielem_shape;
  int inode[MAX_NODES_PER_SIDE];
  int ip, ip_total;
  int a, b, p; /* counters */
  const int dim = pd->Num_Dim;
  int status = 0;
  int v;

  dbl phi_map[MDE], s, t, u;
  dbl xi[DIM];
  dbl x_ip[DIM], interp_val, slope;

  ielem_shape = type2shape(ielem_type);
  ip_total = elem_info(NQUAD_SURF, ielem_type);

  for (p = 0; p < DIM; p++) {
    x_ip[p] = 0.;
  }

  for (ip = 0; ip < ip_total; ip++) {
    for (p = 0; p < DIM; p++) {
      x_n[ip][p] = 0.;
    }

    for (j = 0; j < MDE; j++) {
      phi_v[ip][j] = 0.;
    }

    for (mode = 0; mode < vn->modes; mode++) {
      for (a = 0; a < VIM; a++) {
        for (b = 0; b < VIM; b++) {
          {
            stress_neighbor[ip][mode][a][b] = 0.;
            snv[ip][mode][a][b] = 0.;
          }
        }
      }
    }
  }

  /* first get global node numbers for current element and side */
  iconnect_ptr = Proc_Connect_Ptr[current_elem]; /* find pointer to beginning */

  for (i = 0; i < nodes_per_side; i++) {
    id = local_elem_node_id[i];
    /* load up global node numbers on this side */
    inode[i] = Proc_Elem_Connect[iconnect_ptr + id];
  }

  /* find localside number for neighbor element from
     global node numbers of current element */

  id_side = find_id_side(current_elem, nodes_per_side, inode, id_local_elem_coord, exo);

  for (ip = 0; ip < ip_total; ip++) {

    /* find the quadrature point locations for current ip */
    find_surf_st(ip, ielem_type, id_side, pd->Num_Dim, xi, &s, &t, &u);

    v = pd->ShapeVar;

    iconnect_ptr = Proc_Connect_Ptr[current_elem]; /* find pointer to beginning */
    for (i = 0; i < ei[pg->imtrx]->dof[v]; i++) {
      phi_map[i] = newshape(xi, ielem_type, PSI, ei[pg->imtrx]->dof_list[v][i], ielem_shape,
                            pd->i[pg->imtrx][v], i);
    }

    for (i = 0; i < ei[pg->imtrx]->dof[v]; i++) {
      gnn = Proc_Elem_Connect[iconnect_ptr + i];

      for (p = 0; p < dim; p++) {
        x_n[ip][p] += Coor[p][gnn] * phi_map[i];
      }

    } /* sets global positions for gaussian points */

    for (p = 0; p < dim; p++) {
      x_ip[p] = x_n[ip][p];
    }

    for (mode = 0; mode < vn->modes; mode++) {
      for (a = 0; a < VIM; a++) {
        for (b = 0; b < VIM; b++) {
          if (a <= b) {

            if (table_ibc[mode][a][b] != 0) {
              /*       printf("index %d, x_ip %f,%f
               * \n",BC_Types[table_ibc[mode][a][b]].table->t_index[0],x_ip[0],x_ip[1]); */

              /*set first coordinate of x_ip to interpolation abscissa (assumes only one) */
              if (BC_Types[table_ibc[mode][a][b]].table->t_index[0] == 1) {
                x_ip[0] = x_ip[1];
              }
              interp_val =
                  interpolate_table(BC_Types[table_ibc[mode][a][b]].table, x_ip, &slope, NULL);
              stress_neighbor[ip][mode][a][b] = interp_val;

              /* printf("mode %d a %d b %d: x %f interp_val %f \n",mode,a,b,x_ip[0],interp_val); */

            } else {
              stress_neighbor[ip][mode][a][b] = 0.;
            }
          }
        }
      }
    } /* close mode loop*/
  }   /* close ip loop */

  status = 1;
  return (status);
}

/*****************************************************************************/
/* END of routine neighbor_stress_table */
/*****************************************************************************/

void load_neighbor_pointers(Exo_DB *exo,
                            struct GomaLinearSolverData *ams,
                            int ielem, /* neighbor element */
                            int etype, /* element type */
                            int mode,  /* stress mode */
                            int R_s[MAX_MODES][DIM][DIM],
                            /* Equation number for mode ve_mode */
                            int v_s[MAX_MODES][DIM][DIM],
                            /* Variable number for mode ve_mode */
                            dbl *J_S_S[DIM][DIM][MDE][DIM][DIM][MDE])

/**********************************************************************
 *
 * load_neighbor_pointers:
 *
 *  This routine calculates the pointer array, J_S_S, defined
 *  below:
 *
 * J_S_S: Pointer array. This is a array of pointers to
 *                   d_R_S_ab,i/d_S_v_pq,j
 *        where R_S_ab,i is the residual equation to the
 *        ith dof of the ab stress component of the current element
 *        and S_v_pq,j is the jth dof of the pq stress
 *        component in the current upstream neighbor of the
 *        current element.
 *********************************************************************/
{
  int iconnect_ptr, ldof, v, ln, nunks, nnodes;
  int *enl;
  int dof_list[MAX_VARIABLE_TYPES + MAX_CONC][MDE];
  int gun_list[MAX_VARIABLE_TYPES + MAX_CONC][MDE];
  int dof[MAX_VARIABLE_TYPES + MAX_CONC];
  int gnn, i, j, ie, je, ja, meqn1, meqn2, mvar1, mvar2, eqn, var;
  int I, row_dofs, blk_row, J, K, blk_col;
  int *rpntr, *bpntr, *bindx, *indx, *ija;
  double *a;
  NODAL_VARS_STRUCT *nv;

  (void)memset((void *)dof, 0, sizeof(int) * (MAX_VARIABLE_TYPES + MAX_CONC));
  iconnect_ptr = exo->elem_node_pntr[ielem];
  enl = exo->elem_node_list + iconnect_ptr;
  nnodes = exo->elem_node_pntr[ielem + 1] - iconnect_ptr;

  /*
   *  Formulate a list of the number of degrees of freedom, dof[varType],
   *  in element, ielem. Restrict the list to vartypes of type, stress
   *  mode.  Also form the following maps:
   *     dof_list[v][ldof] -> local dof to local node map
   *     gun_list[v][ldof] -> local dof to proc unknown index.
   */
  for (v = v_s[mode][0][0]; v <= v_s[mode][2][2]; v++) {
    if (Num_Var_In_Type[pg->imtrx][v]) {
      ldof = 0;
      for (ln = 0; ln < nnodes; ln++) {
        /*
         * For this local node "ln", what is the global node number, "gnn"?
         */
        gnn = *(enl + ln);

        /*
         * For this variable at this local node, how many dofs are needed?
         * (according to this element) Note: this can be zero...
         *
         */
        nv = Nodes[gnn]->Nodal_Vars_Info[pg->imtrx];
        nunks = get_nv_ndofs_modMF(nv, v);
        dof[v] += nunks;
        GOMA_EH(nunks, "problem with nun for this var.");

        for (i = 0; i < nunks; i++) {
          dof_list[v][ldof] = ln;
          gun_list[v][ldof] = Index_Solution(gnn, v, 0, i, -1, pg->imtrx);
          ldof++;
        }
      }
    }
  }

  if (strcmp(Matrix_Format, "msr") == 0) {
    ija = ams->bindx;
    a = ams->val;
    for (meqn1 = 0; meqn1 < VIM; meqn1++) {
      for (meqn2 = 0; meqn2 < VIM; meqn2++) {
        if (meqn1 <= meqn2) {
          eqn = R_s[mode][meqn1][meqn2];
          if (pd->e[pg->imtrx][eqn]) {
            for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
              ie = ei[pg->imtrx]->gun_list[eqn][i];
              GOMA_EH(ie, "Bad eqn index.");
              mvar1 = meqn1;
              mvar2 = meqn2;
              if (mvar1 <= mvar2) {
                var = v_s[mode][mvar1][mvar2];
                if (pd->v[pg->imtrx][var]) {
                  for (j = 0; j < dof[var]; j++) {
                    je = gun_list[var][j];
                    GOMA_EH(je, "Bad var index.");
                    ja = (ie == je) ? ie : in_list(je, ija[ie], ija[ie + 1], ija);
                    GOMA_EH(ja, "Could not find vbl in sparse matrix.");
                    J_S_S[meqn1][meqn2][i][mvar1][mvar2][j] = a + ja;
                  }
                }
              }
            }
          }
        }
      }
    }
  } else if (strcmp(Matrix_Format, "vbr") == 0) {
    rpntr = ams->rpntr;
    bpntr = ams->bpntr;
    bindx = ams->bindx;
    indx = ams->indx;
    a = ams->val;
    for (meqn1 = 0; meqn1 < VIM; meqn1++) {
      for (meqn2 = 0; meqn2 < VIM; meqn2++) {
        if (meqn1 <= meqn2) {
          eqn = R_s[mode][meqn1][meqn2];
          if (pd->e[pg->imtrx][eqn]) {
            for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
              ie = ei[pg->imtrx]->gun_list[eqn][i];
              I = ei[pg->imtrx]->gnn_list[eqn][i];
              row_dofs = rpntr[I + 1] - rpntr[I];
              blk_row = Local_Offset[pg->imtrx][I][eqn];
              GOMA_EH(ie, "Bad eqn index.");
              mvar1 = meqn1;
              mvar2 = meqn2;
              if (mvar1 <= mvar2) {
                var = v_s[mode][mvar1][mvar2];
                if (pd->v[pg->imtrx][var]) {
                  for (j = 0; j < dof[var]; j++) {
                    J = *(enl + dof_list[var][j]);
                    K = in_list(J, bpntr[I], bpntr[I + 1], bindx);
                    blk_col = Local_Offset[pg->imtrx][J][var];
                    J_S_S[meqn1][meqn2][i][mvar1][mvar2][j] =
                        a + indx[K] + row_dofs * (blk_col + j) + blk_row + i;
                  }
                }
              }
            }
          }
        }
      }
    }
  } else if (strcmp(Matrix_Format, "epetra") == 0) {
    GOMA_EH(GOMA_ERROR, "load_neighbor_pointers unsupported by epetra");
  }
}
/***************************************************************************/
/***************************************************************************/
/***************************************************************************/
int segregate_stress_update(double x_update[]) {
  int a, b, ieqn, pvar, p, q;
  int eqn, var;
  int ldof_eqn, ldof_var;
  int ln;
  int R_s[DIM][DIM];
  double lump = 0;
  int status = 0;

  if (vn->dg_J_model != SEGREGATED)
    return (0);

  R_s[0][0] = POLYMER_STRESS11;
  R_s[0][1] = POLYMER_STRESS12;
  R_s[0][2] = POLYMER_STRESS13;
  R_s[1][0] = POLYMER_STRESS12;
  R_s[1][1] = POLYMER_STRESS22;
  R_s[1][2] = POLYMER_STRESS23;
  R_s[2][0] = POLYMER_STRESS13;
  R_s[2][1] = POLYMER_STRESS23;
  R_s[2][2] = POLYMER_STRESS33;

  for (ln = 0; ln < ei[pg->imtrx]->num_local_nodes; ln++) {
    for (a = 0; a < DIM; a++) {
      for (b = 0; b < DIM; b++) {
        if (a <= b) {
          eqn = R_s[a][b];

          if (pd->e[pg->imtrx][eqn] && (ldof_eqn = ei[pg->imtrx]->ln_to_first_dof[eqn][ln] != -1)) {
            ieqn = upd->ep[pg->imtrx][eqn];

            while (ldof_eqn <= ei[pg->imtrx]->ln_to_dof[eqn][ln]) {
              for (p = 0; p < DIM; p++) {
                for (q = 0; q < DIM; q++) {
                  if (p <= q && delta(a, p) && delta(q, b)) {
                    var = R_s[p][q];

                    lump = 0.0;

                    if (pd->v[pg->imtrx][var] &&
                        (ldof_var = ei[pg->imtrx]->ln_to_first_dof[var][ln] != -1)) {
                      pvar = upd->vp[pg->imtrx][var];

                      while (ldof_var <= ei[pg->imtrx]->ln_to_dof[var][ln]) {
                        lump += lec->J[LEC_J_INDEX(ieqn, pvar, ldof_eqn, ldof_var)];
                        ldof_var++;
                      }
                    }
                  }
                }
              }
              x_update[ei[pg->imtrx]->gun_list[eqn][ldof_eqn]] =
                  lec->R[LEC_R_INDEX(ieqn, ldof_eqn)] / lump;
              ldof_eqn++;
            }
          }
        }
      }
    }
  }
  return (status);
}

/* This routine calculates the adaptive viscosity from Sun et al., 1999.
 * The adaptive viscosity term multiplies the continuous and discontinuous
 * shear-rate, so it should cancel out and not affect the
 * solution, other than increasing the stability of the
 * algorithm in areas of high shear and stress.
 */
dbl numerical_viscosity(dbl s[DIM][DIM],                        /* total stress */
                        dbl gamma_cont[DIM][DIM],               /* continuous shear rate */
                        dbl d_mun_dS[MAX_MODES][DIM][DIM][MDE], /* derivative of mun wrt S */
                        dbl d_mun_dG[DIM][DIM][MDE])            /* derivative of mun wrt G */
{
  int a, b, j;
  int var, mode;
  int v_s[MAX_MODES][DIM][DIM];
  int v_g[DIM][DIM];
  dbl s_dbl_dot_s;
  dbl g_dbl_dot_g;
  dbl eps2;
  dbl eps; /* should migrate this to input deck */

  dbl mun;

  /* load pointers into equation/variable number */
  (void)stress_eqn_pointer(v_s);

  v_g[0][0] = VELOCITY_GRADIENT11;
  v_g[0][1] = VELOCITY_GRADIENT12;
  v_g[1][0] = VELOCITY_GRADIENT21;
  v_g[1][1] = VELOCITY_GRADIENT22;
  v_g[0][2] = VELOCITY_GRADIENT13;
  v_g[1][2] = VELOCITY_GRADIENT23;
  v_g[2][0] = VELOCITY_GRADIENT31;
  v_g[2][1] = VELOCITY_GRADIENT32;
  v_g[2][2] = VELOCITY_GRADIENT33;

  eps = vn->eps;

  eps2 = eps / 2.;

  s_dbl_dot_s = 0.;
  for (a = 0; a < VIM; a++) {
    for (b = 0; b < VIM; b++) {
      s_dbl_dot_s += s[a][b] * s[b][a];
    }
  }

  g_dbl_dot_g = 0.;
  for (a = 0; a < VIM; a++) {
    for (b = 0; b < VIM; b++) {
      g_dbl_dot_g += gamma_cont[a][b] * gamma_cont[b][a];
    }
  }

  mun = (sqrt(1. + eps2 * s_dbl_dot_s)) / sqrt(1. + eps2 * g_dbl_dot_g);

  for (mode = 0; mode < vn->modes; mode++) {
    for (a = 0; a < VIM; a++) {
      for (b = 0; b < VIM; b++) {
        var = v_s[mode][a][b];

        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_mun_dS[mode][a][b][j] =
              eps2 / (sqrt(1. + eps2 * s_dbl_dot_s) * sqrt(1. + eps2 * g_dbl_dot_g)) * s[a][b] *
              bf[var]->phi[j];
        }
      }
    }
  }

  for (a = 0; a < VIM; a++) {
    for (b = 0; b < VIM; b++) {
      var = v_g[a][b];

      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        /* d_mun_dG[a][b][j] = -eps*(sqrt(1.+ eps2*s_dbl_dot_s))/pow((1.+eps2*g_dbl_dot_g),1.5) */
        /* 	*gamma_cont[a][b]*bf[var]->phi[j]; */
        d_mun_dG[a][b][j] = -eps * mun / (1. + eps2 * g_dbl_dot_g)

                            * gamma_cont[a][b] * bf[var]->phi[j];
      }
    }
  }

  return (mun);
}

/* This routine sets a handy pointer to give the
 * correct equation number for a given mode
 * and a and b component: v[mode][a][b].
 */
int stress_eqn_pointer(int v_s[MAX_MODES][DIM][DIM]) {
  int status = 1;
  /* mode 0, polymer stress */
  v_s[0][0][0] = POLYMER_STRESS11;
  v_s[0][0][1] = POLYMER_STRESS12;
  v_s[0][0][2] = POLYMER_STRESS13;
  v_s[0][1][0] = POLYMER_STRESS12;
  v_s[0][1][1] = POLYMER_STRESS22;
  v_s[0][1][2] = POLYMER_STRESS23;
  v_s[0][2][0] = POLYMER_STRESS13;
  v_s[0][2][1] = POLYMER_STRESS23;
  v_s[0][2][2] = POLYMER_STRESS33;

  /* mode 1, polymer stress */
  v_s[1][0][0] = POLYMER_STRESS11_1;
  v_s[1][0][1] = POLYMER_STRESS12_1;
  v_s[1][0][2] = POLYMER_STRESS13_1;
  v_s[1][1][0] = POLYMER_STRESS12_1;
  v_s[1][1][1] = POLYMER_STRESS22_1;
  v_s[1][1][2] = POLYMER_STRESS23_1;
  v_s[1][2][0] = POLYMER_STRESS13_1;
  v_s[1][2][1] = POLYMER_STRESS23_1;
  v_s[1][2][2] = POLYMER_STRESS33_1;

  /* mode 2, polymer stress */
  v_s[2][0][0] = POLYMER_STRESS11_2;
  v_s[2][0][1] = POLYMER_STRESS12_2;
  v_s[2][0][2] = POLYMER_STRESS13_2;
  v_s[2][1][0] = POLYMER_STRESS12_2;
  v_s[2][1][1] = POLYMER_STRESS22_2;
  v_s[2][1][2] = POLYMER_STRESS23_2;
  v_s[2][2][0] = POLYMER_STRESS13_2;
  v_s[2][2][1] = POLYMER_STRESS23_2;
  v_s[2][2][2] = POLYMER_STRESS33_2;

  /* mode 3, polymer stress */
  v_s[3][0][0] = POLYMER_STRESS11_3;
  v_s[3][0][1] = POLYMER_STRESS12_3;
  v_s[3][0][2] = POLYMER_STRESS13_3;
  v_s[3][1][0] = POLYMER_STRESS12_3;
  v_s[3][1][1] = POLYMER_STRESS22_3;
  v_s[3][1][2] = POLYMER_STRESS23_3;
  v_s[3][2][0] = POLYMER_STRESS13_3;
  v_s[3][2][1] = POLYMER_STRESS23_3;
  v_s[3][2][2] = POLYMER_STRESS33_3;

  /* mode 4, polymer stress */
  v_s[4][0][0] = POLYMER_STRESS11_4;
  v_s[4][0][1] = POLYMER_STRESS12_4;
  v_s[4][0][2] = POLYMER_STRESS13_4;
  v_s[4][1][0] = POLYMER_STRESS12_4;
  v_s[4][1][1] = POLYMER_STRESS22_4;
  v_s[4][1][2] = POLYMER_STRESS23_4;
  v_s[4][2][0] = POLYMER_STRESS13_4;
  v_s[4][2][1] = POLYMER_STRESS23_4;
  v_s[4][2][2] = POLYMER_STRESS33_4;

  /* mode 5, polymer stress */
  v_s[5][0][0] = POLYMER_STRESS11_5;
  v_s[5][0][1] = POLYMER_STRESS12_5;
  v_s[5][0][2] = POLYMER_STRESS13_5;
  v_s[5][1][0] = POLYMER_STRESS12_5;
  v_s[5][1][1] = POLYMER_STRESS22_5;
  v_s[5][1][2] = POLYMER_STRESS23_5;
  v_s[5][2][0] = POLYMER_STRESS13_5;
  v_s[5][2][1] = POLYMER_STRESS23_5;
  v_s[5][2][2] = POLYMER_STRESS33_5;

  /* mode 6, polymer stress */
  v_s[6][0][0] = POLYMER_STRESS11_6;
  v_s[6][0][1] = POLYMER_STRESS12_6;
  v_s[6][0][2] = POLYMER_STRESS13_6;
  v_s[6][1][0] = POLYMER_STRESS12_6;
  v_s[6][1][1] = POLYMER_STRESS22_6;
  v_s[6][1][2] = POLYMER_STRESS23_6;
  v_s[6][2][0] = POLYMER_STRESS13_6;
  v_s[6][2][1] = POLYMER_STRESS23_6;
  v_s[6][2][2] = POLYMER_STRESS33_6;

  /* mode 7, polymer stress */
  v_s[7][0][0] = POLYMER_STRESS11_7;
  v_s[7][0][1] = POLYMER_STRESS12_7;
  v_s[7][0][2] = POLYMER_STRESS13_7;
  v_s[7][1][0] = POLYMER_STRESS12_7;
  v_s[7][1][1] = POLYMER_STRESS22_7;
  v_s[7][1][2] = POLYMER_STRESS23_7;
  v_s[7][2][0] = POLYMER_STRESS13_7;
  v_s[7][2][1] = POLYMER_STRESS23_7;
  v_s[7][2][2] = POLYMER_STRESS33_7;

  return (status);
}

void compute_exp_s(double s[DIM][DIM],
                   double exp_s[DIM][DIM],
                   double eig_values[DIM],
                   double R[DIM][DIM]) {

  int N = VIM;
  int LDA = N;
  int i, j, k;

  int INFO;
  int LWORK = 20;
  double WORK[LWORK];
  double A[VIM * VIM];
  double EIGEN_MAX = sqrt(sqrt(DBL_MAX));
  double eig_S[DIM];
  memset(eig_values, 0.0, sizeof(double) * DIM);
  memset(eig_S, 0.0, sizeof(double) * DIM);
  memset(WORK, 0, sizeof(double) * LWORK);

  // convert to column major
  for (i = 0; i < VIM; i++) {
    for (j = 0; j < VIM; j++) {
      A[i * VIM + j] = s[j][i];
    }
  }

  // eig solver
  dsyev_("V", "U", &N, A, &LDA, eig_S, WORK, &LWORK, &INFO, 1, 1);
  if (INFO > 0)
    fprintf(stderr, "eigsolver not converged %d\n", INFO);
  if (INFO < 0)
    fprintf(stderr, "eigsolver illegal entry %d\n", INFO);

  // transpose (revert to row major)
  for (i = 0; i < VIM; i++) {
    for (j = 0; j < VIM; j++) {
      R[i][j] = A[j * VIM + i];
    }
  }

  // exponentiate diagonal
  for (i = 0; i < VIM; i++) {
    eig_values[i] = MIN(exp(eig_S[i]), EIGEN_MAX);
  }

  memset(exp_s, 0, sizeof(double) * DIM * DIM);
  for (i = 0; i < VIM; i++) {
    for (j = 0; j < VIM; j++) {
      for (k = 0; k < VIM; k++) {
        exp_s[i][j] += R[i][k] * eig_values[k] * R[j][k];
      }
    }
  }

} // End compute_exp_s

void analytical_exp_s(double s[DIM][DIM],
                      double exp_s[DIM][DIM],
                      double eig_values[DIM],
                      double R[DIM][DIM],
                      double d_exp_s_ds[DIM][DIM][DIM][DIM]) {

  double I_S, II_S, disc, off_diag, q, p2, r, p, phi;
  int i, j, k, l, m;
  double tmp, tmp1, tmp2;
  double d_eig_dS[DIM][DIM][DIM], d_R_dS[DIM][DIM][DIM][DIM];

  double B[DIM][DIM], Q1[DIM][DIM], Q2[DIM][DIM];
  double eig_S[DIM] = {0.0, 0.0, 0.0};
  double EIGEN_MAX = sqrt(sqrt(DBL_MAX));
  memset(eig_S, 0, sizeof(double) * DIM);
  memset(d_eig_dS, 0, sizeof(double) * DIM * DIM * DIM);
  memset(d_R_dS, 0, sizeof(double) * DIM * DIM * DIM * DIM);
  if (d_exp_s_ds != NULL)
    memset(d_exp_s_ds, 0, sizeof(double) * DIM * DIM * DIM * DIM);

  /* Use Eigenvalue algorithm from Wikipedia - https://en.wikipedia.org/wiki/
     Eigenvalue_algorithm#Normal%2C_Hermitian%2C_and_real-symmetric_matrices */

  /*  Make sure stress is symmetric  */
  s[1][0] = s[0][1];
  s[2][0] = s[0][2];
  s[2][1] = s[1][2];

  if ((VIM == 2 || pd->CoordinateSystem == CYLINDRICAL)) {
    eig_S[2] = s[2][2];
    d_eig_dS[2][2][2] = 1.0;
    I_S = s[0][0] + s[1][1];
    II_S = s[0][0] * s[1][1] - SQUARE(s[0][1]);
    disc = sqrt(SQUARE(I_S) - 4 * II_S);
    eig_S[1] = 0.5 * (I_S + disc);
    eig_S[0] = 0.5 * (I_S - disc);

    if (d_exp_s_ds != NULL) {
      d_eig_dS[1][0][0] = 0.5 + 0.5 * (s[0][0] - s[1][1]) / disc;
      d_eig_dS[1][1][1] = 0.5 - 0.5 * (s[0][0] - s[1][1]) / disc;
      d_eig_dS[1][0][1] = d_eig_dS[1][1][0] = 4. * s[0][1] / disc;
      d_eig_dS[0][0][0] = 0.5 - 0.5 * (s[0][0] - s[1][1]) / disc;
      d_eig_dS[0][1][1] = 0.5 + 0.5 * (s[0][0] - s[1][1]) / disc;
      d_eig_dS[0][0][1] = d_eig_dS[0][1][0] = -4. * s[0][1] / disc;
    }
    for (i = 0; i < VIM; i++) {
      for (j = 0; j < VIM; j++) {
        R[i][j] = delta(i, j);
      }
    }
    if (DOUBLE_NONZERO(II_S)) {
      for (i = 0; i < 2; i++) {
        tmp1 = sqrt(SQUARE(s[0][1]) + SQUARE(s[0][0] - eig_S[i]));
        tmp2 = sqrt(SQUARE(s[0][1]) + SQUARE(s[1][1] - eig_S[i]));
        if (DOUBLE_NONZERO(tmp1)) {
          R[i][0] = s[0][1] / tmp1;
          R[i][1] = (eig_S[i] - s[0][0]) / tmp1;
          if (d_exp_s_ds != NULL) {
            d_R_dS[i][0][0][0] =
                -s[0][1] * (s[0][0] - eig_S[i]) * (1. - d_eig_dS[i][0][0]) / CUBE(tmp1);
            d_R_dS[i][0][0][1] = d_R_dS[i][0][1][0] =
                (SQUARE(s[0][0] - eig_S[i]) + s[0][1] * (s[0][0] - eig_S[i]) * d_eig_dS[i][0][1]) /
                CUBE(tmp1);
            d_R_dS[i][0][1][1] = s[0][1] * (s[0][0] - eig_S[i]) * d_eig_dS[i][1][1] / CUBE(tmp1);
            d_R_dS[i][1][0][0] = -SQUARE(s[0][1]) * (1. - d_eig_dS[i][0][0]) / CUBE(tmp1);
            d_R_dS[i][1][0][1] = d_R_dS[i][1][1][0] =
                (s[0][1] * (s[0][0] - eig_S[i]) + SQUARE(s[0][1]) * d_eig_dS[i][0][1]) / CUBE(tmp1);
            d_R_dS[i][1][1][1] = SQUARE(s[0][1]) * d_eig_dS[i][1][1] / CUBE(tmp1);
          }
        } else if (DOUBLE_NONZERO(tmp2)) {
          R[i][0] = (eig_S[i] - s[1][1]) / tmp2;
          R[i][1] = s[0][1] / tmp2;
          if (d_exp_s_ds != NULL) {
            d_R_dS[i][0][1][1] = SQUARE(s[0][1]) * (1. - d_eig_dS[i][1][1]) / CUBE(tmp2);
            d_R_dS[i][0][0][1] = d_R_dS[i][0][1][0] =
                -(s[0][1] * (s[1][1] - eig_S[i]) + SQUARE(s[0][1]) * d_eig_dS[i][0][1]) /
                CUBE(tmp2);
            d_R_dS[i][0][0][0] = -SQUARE(s[0][1]) * d_eig_dS[i][0][0] / CUBE(tmp2);
            d_R_dS[i][1][1][1] =
                s[0][1] * (s[1][1] - eig_S[i]) * (1. - d_eig_dS[i][1][1]) / CUBE(tmp2);
            d_R_dS[i][1][0][1] = d_R_dS[i][1][1][0] =
                -(SQUARE(s[1][1] - eig_S[i]) + s[0][1] * (s[1][1] - eig_S[i]) * d_eig_dS[i][0][1]) /
                CUBE(tmp2);
            d_R_dS[i][1][0][0] = -s[0][1] * (s[1][1] - eig_S[i]) * d_eig_dS[i][0][0] / CUBE(tmp2);
          }
        } else {
          R[i][i] = 1.0;
        }
      }
    } else {
      for (i = 0; i < VIM; i++) {
        R[i][i] = 1.0;
      }
    }
  } else {
    off_diag = SQUARE(s[0][1]) + SQUARE(s[0][2]) + SQUARE(s[1][2]);
    I_S = s[0][0] + s[1][1] + s[2][2];
    II_S = s[0][0] * s[1][1] + s[0][0] * s[2][2] + s[1][1] * s[2][2] -
           (SQUARE(s[0][1]) + SQUARE(s[0][2]) + SQUARE(s[1][2]));
    if (DOUBLE_NONZERO(off_diag)) {
      q = I_S / 3.;
      p2 = 2 * off_diag;
      for (i = 0; i < VIM; i++) {
        p2 += SQUARE(s[i][i] - q);
      }
      p = sqrt(p2 / 6.);
      for (i = 0; i < VIM; i++) {
        for (j = 0; j < VIM; j++) {
          B[i][j] = (s[i][j] - q * delta(i, j)) / p;
        }
      }
      r = B[0][0] * B[1][1] * B[2][2] + 2. * (B[0][1] * B[1][2] * B[0][2]) -
          B[0][0] * SQUARE(B[1][2]) - B[1][1] * SQUARE(B[0][2]) - B[2][2] * SQUARE(B[0][1]);
      r *= 0.5;
      if (r <= -1.) {
        phi = M_PIE / 3.;
      } else if (r >= 1.) {
        phi = 0.;
      } else {
        phi = acos(r) / 3.;
      }
      for (i = 0; i < VIM; i++) {
        eig_S[i] = q + 2 * p * cos(phi + i * 2 * M_PIE / 3.);
      }
    } else {
      for (i = 0; i < DIM; i++) {
        eig_S[i] = s[i][i];
      }
    }

    /*   c
    c       solve 3x3 linear system
    c
            deter=xj(1,1)*(xj(2,2)*xj(3,3)-xj(3,2)*xj(2,3))
         1          -xj(2,1)*(xj(1,2)*xj(3,3)-xj(3,2)*xj(1,3))
         2          +xj(3,1)*(xj(2,3)*xj(1,2)-xj(2,2)*xj(1,3))
            detinv=1./deter
            u(1)=detinv*((-r(1))*(xj(3,3)*xj(2,2)-xj(3,2)*xj(2,3))
         1          +(-r(2))*(xj(3,2)*xj(1,3)-xj(3,3)*xj(1,2))
         2          +(-r(3))*(xj(2,3)*xj(1,2)-xj(2,2)*xj(1,3)))
            u(2)=detinv*((-r(1))*(xj(3,1)*xj(2,3)-xj(3,3)*xj(2,1))
         1          +(-r(2))*(xj(3,3)*xj(1,1)-xj(3,1)*xj(1,3))
         2          +(-r(3))*(xj(2,1)*xj(1,3)-xj(2,3)*xj(1,1)))
            u(3)=detinv*((-r(1))*(xj(3,2)*xj(2,1)-xj(3,1)*xj(2,2))
         1          +(-r(2))*(xj(3,1)*xj(1,2)-xj(3,2)*xj(1,1))
         2          +(-r(3))*(xj(2,2)*xj(1,1)-xj(2,1)*xj(1,2)))
    */
    if (DOUBLE_NONZERO(I_S) || DOUBLE_NONZERO(II_S)) {
      for (i = 0; i < VIM; i++) {
        for (j = 0; j < VIM; j++) {
          for (k = 0; k < VIM; k++) {
            Q1[j][k] = s[j][k] - eig_S[(i + 1) % 3] * delta(j, k);
            Q2[j][k] = s[j][k] - eig_S[(i + 2) % 3] * delta(j, k);
          }
        }
        for (j = 0; j < VIM; j++) {
          for (k = 0; k < VIM; k++) {
            R[j][i] = Q1[j][k] * Q2[k][j];
          }
        }
      }
      /*  Normalize eigenvectors   */
      for (j = 0; j < VIM; j++) {
        tmp = sqrt(SQUARE(R[0][j]) + SQUARE(R[1][j]) + SQUARE(R[2][j]));
        if (DOUBLE_NONZERO(tmp)) {
          for (i = 0; i < VIM; i++) {
            R[i][j] /= tmp;
          }
        }
      }
    } else {
      for (j = 0; j < VIM; j++) {
        for (k = 0; k < VIM; k++) {
          R[j][k] = delta(j, k);
        }
      }
    }
  }

  // exponentiate diagonal
  for (i = 0; i < DIM; i++) {
    eig_values[i] = MIN(exp(eig_S[i]), EIGEN_MAX);
  }

  memset(exp_s, 0, sizeof(double) * DIM * DIM);
  for (i = 0; i < VIM; i++) {
    for (j = 0; j < VIM; j++) {
      for (k = 0; k < VIM; k++) {
        exp_s[i][j] += R[i][k] * eig_values[k] * R[j][k];
      }
    }
  }
  if (d_exp_s_ds != NULL) {
    for (i = 0; i < VIM; i++) {
      for (j = 0; j < VIM; j++) {
        for (k = 0; k < VIM; k++) {
          for (l = 0; l < VIM; l++) {
            for (m = 0; m < VIM; m++) {
              d_exp_s_ds[i][j][l][m] += R[i][k] * eig_values[k] * d_R_dS[j][k][l][m];
            }
          }
        }
      }
    }
    for (i = 0; i < VIM; i++) {
      for (j = 0; j < VIM; j++) {
        for (k = 0; k < VIM; k++) {
          for (l = 0; l < VIM; l++) {
            for (m = 0; m < VIM; m++) {
              d_exp_s_ds[i][j][l][m] += d_R_dS[i][k][l][m] * eig_values[k] * R[j][k];
            }
          }
        }
      }
    }
    for (i = 0; i < VIM; i++) {
      for (j = 0; j < VIM; j++) {
        for (k = 0; k < VIM; k++) {
          for (l = 0; l < VIM; l++) {
            for (m = 0; m < VIM; m++) {
              d_exp_s_ds[i][j][l][m] += R[i][k] * eig_values[k] * d_eig_dS[k][l][m] * R[j][k];
            }
          }
        }
      }
    }
  }

} // End analytical_exp_s

void compute_d_exp_s_ds(dbl s[DIM][DIM], // s - stress
                        dbl exp_s[DIM][DIM],
                        dbl d_exp_s_ds[DIM][DIM][DIM][DIM]) {
  double s_p[DIM][DIM], s_n[DIM][DIM];
  double exp_s_p[DIM][DIM], exp_s_n[DIM][DIM];
  double eig_values[DIM];
  double R1[DIM][DIM];
  int i, j, p, q;
  double ds, ds_den, fd = FD_FACTOR;

  memset(d_exp_s_ds, 0, sizeof(double) * DIM * DIM * DIM * DIM);

  for (i = 0; i < VIM; i++) {
    for (j = 0; j < VIM; j++) {
      s_p[i][j] = s_n[i][j] = s[i][j];
    }
  }

  for (i = 0; i < VIM; i++) {
    for (j = i; j < VIM; j++) {

      // perturb s
      ds = fd * s[i][j];
      ds = (fabs(ds) < fd ? fd : ds);
      s_p[i][j] += ds;
      s_n[i][j] -= ds;
      if (i != j) {
        s_p[j][i] = s_p[i][j];
        s_n[j][i] = s_n[i][j];
      }

      ds_den = 0.5 / ds;

      // find exp_s at perturbed value
      compute_exp_s(s_p, exp_s_p, eig_values, R1);
      compute_exp_s(s_n, exp_s_n, eig_values, R1);

      // approximate derivative
      for (p = 0; p < VIM; p++) {
        for (q = 0; q < VIM; q++) {
          d_exp_s_ds[p][q][i][j] = ds_den * (exp_s_p[p][q] - exp_s_n[p][q]);
          if (i != j)
            d_exp_s_ds[p][q][j][i] = d_exp_s_ds[p][q][i][j];
        }
      }
      s_p[i][j] = s_n[i][j] = s[i][j];
      if (i != j) {
        s_p[j][i] = s_n[j][i] = s[j][i];
      }
    }
  }
}
/*****************************************************************************/
void compute_saramito_model_terms(dbl *sCoeff,
                                  SARAMITO_DEPENDENCE_STRUCT *d_sCoeff,
                                  const dbl stress[DIM][DIM],
                                  const struct Generalized_Newtonian *gn_local,
                                  const int normalizeCoeff) {
  /* start by computing the norm of sigma_d (the deviatoric stress) which is written as
   * sqrt( trace(sigma_d*sigma_d)/2 )
   *
   * see the following wikipedia page for the mathematical details:
   * https://en.wikipedia.org/wiki/Cauchy_stress_tensor#Invariants_of_the_stress_deviator_tensor
   */

  const dbl yieldStress = gn_local->tau_y;
  const dbl yieldExpon = gn_local->fexp;
  const dbl m = 1. / gn_local->nexp;

  dbl traceOverVIM = 0;
  for (int i = 0; i < VIM; i++) {
    traceOverVIM += stress[i][i];
  }

  traceOverVIM /= VIM;

  // square of the deviatoric sress norm
  dbl normOfStressDSqr = 0;
  for (int i = 0; i < VIM; i++) {
    normOfStressDSqr += pow(stress[i][i] - traceOverVIM, 2) / 2.;

    for (int j = i + 1; j < VIM; j++) {
      normOfStressDSqr += pow(stress[i][j], 2);
    }
  }

  const dbl normOfStressD = sqrt(normOfStressDSqr) + 1e-16;

  const dbl sc = 1. - yieldStress / normOfStressD;

  //  phi  = Saramito coefficient with nexp = 1.
  // This is useful for computing sensitivities
  dbl phi;
  dbl expYSC = 1;
  if (yieldExpon > 0) {
    expYSC = exp(sc / yieldExpon);
    phi = log(1. + expYSC) * yieldExpon;
  } else {
    phi = fmax(0, sc);
  }
  const dbl beta = m == 1 ? 1 : pow(normOfStressD, m - 1);

  *sCoeff = m == 1 ? beta * phi : beta * pow(phi, m);

  // take care of indeterminate behavior for normOfStressD == 0
  if (normOfStressD == 0) {
    if (yieldStress <= 0) {
      *sCoeff = beta;
    } else {
      *sCoeff = 0;
    }
  }

  if (normalizeCoeff) {
    *sCoeff /= beta;
  }

  if (d_sCoeff != NULL) {
    // if normStress_d < yieldStress, set sensitivities to zero and return
    if ((*sCoeff) == 0 || yieldStress <= 0) {
      for (int i = 0; i < VIM; i++) {
        for (int j = 0; j < VIM; j++) {
          d_sCoeff->s[i][j] = 0.;
        }
      }

      d_sCoeff->tau_y = 0.;

      // otherwise, sensitivities need to be calculated
    } else {

      d_sCoeff->tau_y = m == 1 ? -1. / normOfStressD
                               : -m * pow(normOfStressD - yieldStress, m - 1) / normOfStressD;

      dbl d_sCoeff_d_normOfStressD = m == 1 ? yieldStress / normOfStressDSqr
                                            : pow(normOfStressD - yieldStress, m - 1) *
                                                  ((m - 1) * normOfStressD + yieldStress) /
                                                  normOfStressDSqr;
      if (yieldExpon > 0) {
        const dbl expYSCDerivativeTerm = expYSC / (1 + expYSC);
        d_sCoeff->tau_y *= expYSCDerivativeTerm;
        d_sCoeff_d_normOfStressD *= expYSCDerivativeTerm;
      }

      // fill d_sCoeff->s[i][j] with d(normOfStressDSqr)/d(stress[i][j])
      for (int i = 0; i < VIM; i++) {
        d_sCoeff->s[i][i] = -traceOverVIM;

        for (int j = i; j < VIM; j++) {
          if (i == j) {
            d_sCoeff->s[i][i] += stress[i][i];
          } else {
            d_sCoeff->s[i][j] = 2. * stress[i][j];
          }
        }
      }

      /* use the chain rule to computute sCoeff sensitivies to stress components
       * d(sCoeff)/d(stress)  =   d(sCoeff)/d(normOfStressD)
       *                        * d(normOfStressD)/d(normOfStressDSqr)
       *                        * d(normOfStressDSqr)/d(stress)
       *
       *                      =   d(sCoeff)/d(normOfStressD)
       *                        * 0.5/normOfStressD
       *                        * d(normOfStressDSqr)/d(stress)
       *
       */

      dbl d_normOfStressD_d_Stress = 0.5 / normOfStressD * d_sCoeff_d_normOfStressD;

      for (int i = 0; i < VIM; i++) {
        for (int j = i; j < VIM; j++) {
          d_sCoeff->s[i][j] *= d_normOfStressD_d_Stress;
          d_sCoeff->s[j][i] = d_sCoeff->s[i][j];
        }
      }
    }
  } // d_sCoeff != NULL
  // for debugging purposes.
  /*
    for(int i=0; i<VIM; i++){
    for(int j=0; j<VIM; j++){
    printf("\nstress[%d][%d] = %E", i, j, stress[i][j]);
    }
    }
    printf("\n");
    for(int i=0; i<VIM; i++){
    for(int j=0; j<VIM; j++){
    printf("\nd(sCoeff)/d(stress[%d][%d)] = %E", i, j, d_sCoeff->s[i][j]);
    }
    }
    printf("\n");
    printf("yield stress = %E", yieldStress);
    printf("\n");
    printf("|stress_d|**2 = %E", normOfStressDSqr);
    printf("\n-------------------------------------------------------------");
  */
}

void compute_a_dot_b(dbl b[DIM][DIM],
                     dbl G[DIM][DIM],
                     dbl a_dot_b[DIM][DIM],
                     dbl d_a_dot_b_db[DIM][DIM][DIM][DIM],
                     dbl d_a_dot_b_dG[DIM][DIM][DIM][DIM]) {

  if (VIM == 2) {

    dbl a12 = ((b[0][1] * G[0][0] - b[0][0] * G[0][1]) + (b[1][1] * G[1][0] - b[1][0] * G[1][1])) /
              (b[0][0] + b[1][1] + 1e-16);

    dbl a[DIM][DIM] = {{0., a12, 0.}, {-a12, 0., 0.}, {0., 0., 0.}};

    tensor_dot(a, b, a_dot_b, VIM);

    if (af->Assemble_Jacobian) {

      if (pd->v[pg->imtrx][R_GRADIENT11]) {
        for (int i = 0; i < VIM; i++) {
          for (int j = 0; j < VIM; j++) {
            dbl da12 =
                ((b[0][1] * delta(i, 0) * delta(j, 0) - b[0][0] * delta(i, 0) * delta(j, 1)) +
                 (b[1][1] * delta(i, 1) * delta(j, 0) - b[1][0] * delta(i, 1) * delta(j, 1))) /
                (b[0][0] + b[1][1] + 1e-16);

            dbl da[DIM][DIM] = {{0., da12, 0.}, {-da12, 0., 0.}, {0., 0., 0.}};

            tensor_dot(da, b, d_a_dot_b_dG[i][j], VIM);
          }
        }
      }

      for (int i = 0; i < VIM; i++) {
        for (int j = 0; j < VIM; j++) {

          dbl da12 =
              ((-delta(i, 0) * delta(j, 0) * G[0][1] +
                (delta(i, 0) * delta(j, 1) + delta(i, 1) * delta(j, 0)) * (G[0][0] - G[1][1]) +
                delta(i, 1) * delta(j, 1) * G[1][0]) -
               (1 - delta(i, 0) * delta(j, 1) - delta(i, 1) * delta(j, 0)) * a12) /
              (b[0][0] + b[1][1] + 1e-16);

          dbl da[DIM][DIM] = {{0., da12, 0.}, {-da12, 0., 0.}, {0., 0., 0.}};
          dbl db[DIM][DIM];
          for (int p = 0; p < VIM; p++) {
            for (int q = 0; q < VIM; q++) {
              db[p][q] = delta(i, p) * delta(j, q) | delta(i, q) * delta(j, p);
            }
          }
          dbl a_dot_db[DIM][DIM];
          dbl da_dot_b[DIM][DIM];

          tensor_dot(da, b, da_dot_b, VIM);
          tensor_dot(a, db, a_dot_db, VIM);
          for (int p = 0; p < VIM; p++) {
            for (int q = 0; q < VIM; q++) {
              d_a_dot_b_db[i][j][p][q] = da_dot_b[p][q] + a_dot_db[p][q];
            }
          }
        }
      }
    }

  } else { // VIM = 3
    dbl D = -b[0][1] * (b[0][1] * (b[0][0] + b[1][1]) + b[0][2] * b[1][2]) -
            b[0][2] * (b[0][1] * b[1][2] + b[0][2] * (b[0][0] + b[2][2])) +
            (b[1][1] + b[2][2]) * (-b[1][2] * b[1][2] + (b[0][0] + b[1][1]) * (b[0][0] + b[2][2])) +
            1e-16;
    dbl invD = 1.0 / D;

    dbl a12 = invD * (-pow(b[0][1], 2) + (b[0][0] + b[2][2]) * (b[1][1] + b[2][2])) *
                  (G[0][0] * b[0][1] - G[0][1] * b[0][0] + G[1][0] * b[1][1] - G[1][1] * b[0][1] +
                   G[2][0] * b[1][2] - G[2][1] * b[0][2]) -
              (b[0][1] * b[0][2] + b[1][2] * (b[1][1] + b[2][2])) *
                  (G[0][0] * b[0][2] - G[0][2] * b[0][0] + G[1][0] * b[1][2] - G[1][2] * b[0][1] +
                   G[2][0] * b[2][2] - G[2][2] * b[0][2]) +
              (b[0][1] * b[1][2] + b[0][2] * (b[0][0] + b[2][2])) *
                  (G[0][1] * b[0][2] - G[0][2] * b[0][1] + G[1][1] * b[1][2] - G[1][2] * b[1][1] +
                   G[2][1] * b[2][2] - G[2][2] * b[1][2]);

    dbl a13 =
        invD *
        ((-pow(b[0][2], 2) + (b[0][0] + b[1][1]) * (b[1][1] + b[2][2])) *
             (G[0][0] * b[0][2] - G[0][2] * b[0][0] + G[1][0] * b[1][2] - G[1][2] * b[0][1] +
              G[2][0] * b[2][2] - G[2][2] * b[0][2]) +
         (-b[0][1] * b[0][2] - b[1][2] * (b[1][1] + b[2][2])) *
             (G[0][0] * b[0][1] - G[0][1] * b[0][0] + G[1][0] * b[1][1] - G[1][1] * b[0][1] +
              G[2][0] * b[1][2] - G[2][1] * b[0][2]) -
         (b[0][1] * (b[0][0] + b[1][1]) + b[0][2] * b[1][2]) *
             (G[0][1] * b[0][2] - G[0][2] * b[0][1] + G[1][1] * b[1][2] - G[1][2] * b[1][1] +
              G[2][1] * b[2][2] - G[2][2] * b[1][2])) /
        (-b[0][1] * (b[0][1] * (b[0][0] + b[1][1]) + b[0][2] * b[1][2]) -
         b[0][2] * (b[0][1] * b[1][2] + b[0][2] * (b[0][0] + b[2][2])) +
         (b[1][1] + b[2][2]) * (-pow(b[1][2], 2) + (b[0][0] + b[1][1]) * (b[0][0] + b[2][2])));

    dbl a23 = invD * (-pow(b[1][2], 2) + (b[0][0] + b[1][1]) * (b[0][0] + b[2][2])) *
                  (G[0][1] * b[0][2] - G[0][2] * b[0][1] + G[1][1] * b[1][2] - G[1][2] * b[1][1] +
                   G[2][1] * b[2][2] - G[2][2] * b[1][2]) +
              (b[0][1] * b[1][2] + b[0][2] * (b[0][0] + b[2][2])) *
                  (G[0][0] * b[0][1] - G[0][1] * b[0][0] + G[1][0] * b[1][1] - G[1][1] * b[0][1] +
                   G[2][0] * b[1][2] - G[2][1] * b[0][2]) -
              (b[0][1] * (b[0][0] + b[1][1]) + b[0][2] * b[1][2]) *
                  (G[0][0] * b[0][2] - G[0][2] * b[0][0] + G[1][0] * b[1][2] - G[1][2] * b[0][1] +
                   G[2][0] * b[2][2] - G[2][2] * b[0][2]);

    dbl a[DIM][DIM] = {
        {0.0, a12, a13},
        {-a12, 0.0, a23},
        {-a13, -a23, 0.0},
    };

    tensor_dot(a, b, a_dot_b, VIM);

    if (af->Assemble_Jacobian) {

      dbl d_a12_dG[DIM][DIM];
      dbl d_a13_dG[DIM][DIM];
      dbl d_a23_dG[DIM][DIM];
      d_a12_dG[0][0] = (-b[0][1] * (pow(b[0][1], 2) - (b[0][0] + b[2][2]) * (b[1][1] + b[2][2])) -
                        b[0][2] * (b[0][1] * b[0][2] + b[1][2] * (b[1][1] + b[2][2]))) *
                       invD;
      d_a13_dG[0][0] = (-b[0][1] * (b[0][1] * b[0][2] + b[1][2] * (b[1][1] + b[2][2])) -
                        b[0][2] * (pow(b[0][2], 2) - (b[0][0] + b[1][1]) * (b[1][1] + b[2][2]))) *
                       invD;
      d_a23_dG[0][0] = (b[0][1] * (b[0][1] * b[1][2] + b[0][2] * (b[0][0] + b[2][2])) -
                        b[0][2] * (b[0][1] * (b[0][0] + b[1][1]) + b[0][2] * b[1][2])) *
                       invD;
      d_a12_dG[0][1] = (b[0][0] * (pow(b[0][1], 2) - (b[0][0] + b[2][2]) * (b[1][1] + b[2][2])) +
                        b[0][2] * (b[0][1] * b[1][2] + b[0][2] * (b[0][0] + b[2][2]))) *
                       invD;
      d_a13_dG[0][1] = (b[0][0] * (b[0][1] * b[0][2] + b[1][2] * (b[1][1] + b[2][2])) -
                        b[0][2] * (b[0][1] * (b[0][0] + b[1][1]) + b[0][2] * b[1][2])) *
                       invD;
      d_a23_dG[0][1] = (-b[0][0] * b[0][1] * b[1][2] + b[0][0] * b[0][2] * b[1][1] +
                        b[0][2] * b[1][1] * b[2][2] - b[0][2] * pow(b[1][2], 2)) *
                       invD;
      d_a12_dG[1][0] = (-b[1][1] * (pow(b[0][1], 2) - (b[0][0] + b[2][2]) * (b[1][1] + b[2][2])) -
                        b[1][2] * (b[0][1] * b[0][2] + b[1][2] * (b[1][1] + b[2][2]))) *
                       invD;
      d_a13_dG[1][0] = (b[0][0] * b[1][1] * b[1][2] + b[0][0] * b[1][2] * b[2][2] -
                        b[0][1] * b[0][2] * b[1][1] - pow(b[0][2], 2) * b[1][2]) *
                       invD;
      d_a23_dG[1][0] = (b[1][1] * (b[0][1] * b[1][2] + b[0][2] * (b[0][0] + b[2][2])) -
                        b[1][2] * (b[0][1] * (b[0][0] + b[1][1]) + b[0][2] * b[1][2])) *
                       invD;
      d_a12_dG[1][1] = (b[0][1] * (pow(b[0][1], 2) - (b[0][0] + b[2][2]) * (b[1][1] + b[2][2])) +
                        b[1][2] * (b[0][1] * b[1][2] + b[0][2] * (b[0][0] + b[2][2]))) *
                       invD;
      d_a13_dG[1][1] = (b[0][1] * (b[0][1] * b[0][2] + b[1][2] * (b[1][1] + b[2][2])) -
                        b[1][2] * (b[0][1] * (b[0][0] + b[1][1]) + b[0][2] * b[1][2])) *
                       invD;
      d_a23_dG[1][1] = (-b[0][1] * (b[0][1] * b[1][2] + b[0][2] * (b[0][0] + b[2][2])) -
                        b[1][2] * (pow(b[1][2], 2) - (b[0][0] + b[1][1]) * (b[0][0] + b[2][2]))) *
                       invD;
      d_a12_dG[2][2] = (b[0][2] * (b[0][1] * b[0][2] + b[1][2] * (b[1][1] + b[2][2])) -
                        b[1][2] * (b[0][1] * b[1][2] + b[0][2] * (b[0][0] + b[2][2]))) *
                       invD;
      d_a13_dG[2][2] = (b[0][2] * (pow(b[0][2], 2) - (b[0][0] + b[1][1]) * (b[1][1] + b[2][2])) +
                        b[1][2] * (b[0][1] * (b[0][0] + b[1][1]) + b[0][2] * b[1][2])) *
                       invD;
      d_a23_dG[2][2] = (b[0][2] * (b[0][1] * (b[0][0] + b[1][1]) + b[0][2] * b[1][2]) +
                        b[1][2] * (pow(b[1][2], 2) - (b[0][0] + b[1][1]) * (b[0][0] + b[2][2]))) *
                       invD;

      if (pd->v[pg->imtrx][R_GRADIENT11]) {
        for (int i = 0; i < VIM; i++) {
          for (int j = 0; j < VIM; j++) {
            dbl da12 = d_a12_dG[i][j];
            dbl da13 = d_a13_dG[i][j];
            dbl da23 = d_a23_dG[i][j];

            dbl da[DIM][DIM] = {{0., da12, da13}, {-da12, 0., da23}, {-da13, -da23, 0.}};

            tensor_dot(da, b, d_a_dot_b_dG[i][j], VIM);
          }
        }
      }

      dbl d_D_db[DIM][DIM];

      d_D_db[0][0] = -pow(b[0][1], 2) - pow(b[0][2], 2) +
                     (b[1][1] + b[2][2]) * (2 * b[0][0] + b[1][1] + b[2][2]);
      d_D_db[0][1] = -2 * b[0][1] * (b[0][0] + b[1][1]) - 2 * b[0][2] * b[1][2];
      d_D_db[0][2] = -2 * b[0][1] * b[1][2] - 2 * b[0][2] * (b[0][0] + b[2][2]);
      d_D_db[1][0] = -2 * b[0][1] * (b[0][0] + b[1][1]) - 2 * b[0][2] * b[1][2];
      d_D_db[1][1] = -pow(b[0][1], 2) - pow(b[1][2], 2) +
                     (b[0][0] + b[1][1]) * (b[0][0] + b[2][2]) +
                     (b[0][0] + b[2][2]) * (b[1][1] + b[2][2]);
      d_D_db[1][2] = -2 * b[0][1] * b[0][2] - 2 * b[1][2] * (b[1][1] + b[2][2]);
      d_D_db[2][0] = -2 * b[0][1] * b[1][2] - 2 * b[0][2] * (b[0][0] + b[2][2]);
      d_D_db[2][1] = -2 * b[0][1] * b[0][2] - 2 * b[1][2] * (b[1][1] + b[2][2]);
      d_D_db[2][2] = -pow(b[0][2], 2) - pow(b[1][2], 2) +
                     (b[0][0] + b[1][1]) * (b[0][0] + b[2][2]) +
                     (b[0][0] + b[1][1]) * (b[1][1] + b[2][2]);

      dbl d_Da12_db[DIM][DIM];
      dbl d_Da13_db[DIM][DIM];
      dbl d_Da23_db[DIM][DIM];

      d_Da12_db[0][0] =
          G[0][1] * (pow(b[0][1], 2) - (b[0][0] + b[2][2]) * (b[1][1] + b[2][2])) +
          G[0][2] * (b[0][1] * b[0][2] + b[1][2] * (b[1][1] + b[2][2])) +
          b[0][2] * (G[0][1] * b[0][2] - G[0][2] * b[0][1] + G[1][1] * b[1][2] - G[1][2] * b[1][1] +
                     G[2][1] * b[2][2] - G[2][2] * b[1][2]) +
          (b[1][1] + b[2][2]) * (G[0][0] * b[0][1] - G[0][1] * b[0][0] + G[1][0] * b[1][1] -
                                 G[1][1] * b[0][1] + G[2][0] * b[1][2] - G[2][1] * b[0][2]);
      d_Da13_db[0][0] =
          G[0][1] * (b[0][1] * b[0][2] + b[1][2] * (b[1][1] + b[2][2])) +
          G[0][2] * (pow(b[0][2], 2) - (b[0][0] + b[1][1]) * (b[1][1] + b[2][2])) -
          b[0][1] * (G[0][1] * b[0][2] - G[0][2] * b[0][1] + G[1][1] * b[1][2] - G[1][2] * b[1][1] +
                     G[2][1] * b[2][2] - G[2][2] * b[1][2]) +
          (b[1][1] + b[2][2]) * (G[0][0] * b[0][2] - G[0][2] * b[0][0] + G[1][0] * b[1][2] -
                                 G[1][2] * b[0][1] + G[2][0] * b[2][2] - G[2][2] * b[0][2]);
      d_Da23_db[0][0] = -G[0][1] * (b[0][1] * b[1][2] + b[0][2] * (b[0][0] + b[2][2])) +
                        G[0][2] * (b[0][1] * (b[0][0] + b[1][1]) + b[0][2] * b[1][2]) -
                        b[0][1] * (G[0][0] * b[0][2] - G[0][2] * b[0][0] + G[1][0] * b[1][2] -
                                   G[1][2] * b[0][1] + G[2][0] * b[2][2] - G[2][2] * b[0][2]) +
                        b[0][2] * (G[0][0] * b[0][1] - G[0][1] * b[0][0] + G[1][0] * b[1][1] -
                                   G[1][1] * b[0][1] + G[2][0] * b[1][2] - G[2][1] * b[0][2]) +
                        (2 * b[0][0] + b[1][1] + b[2][2]) *
                            (G[0][1] * b[0][2] - G[0][2] * b[0][1] + G[1][1] * b[1][2] -
                             G[1][2] * b[1][1] + G[2][1] * b[2][2] - G[2][2] * b[1][2]);
      d_Da12_db[0][1] =
          -G[0][2] * (b[0][1] * b[1][2] + b[0][2] * (b[0][0] + b[2][2])) +
          G[1][2] * (b[0][1] * b[0][2] + b[1][2] * (b[1][1] + b[2][2])) -
          2 * b[0][1] *
              (G[0][0] * b[0][1] - G[0][1] * b[0][0] + G[1][0] * b[1][1] - G[1][1] * b[0][1] +
               G[2][0] * b[1][2] - G[2][1] * b[0][2]) -
          b[0][2] * (G[0][0] * b[0][2] - G[0][2] * b[0][0] + G[1][0] * b[1][2] - G[1][2] * b[0][1] +
                     G[2][0] * b[2][2] - G[2][2] * b[0][2]) +
          b[1][2] * (G[0][1] * b[0][2] - G[0][2] * b[0][1] + G[1][1] * b[1][2] - G[1][2] * b[1][1] +
                     G[2][1] * b[2][2] - G[2][2] * b[1][2]) -
          (G[0][0] - G[1][1]) * (pow(b[0][1], 2) - (b[0][0] + b[2][2]) * (b[1][1] + b[2][2]));
      d_Da13_db[0][1] =
          G[0][2] * (b[0][1] * (b[0][0] + b[1][1]) + b[0][2] * b[1][2]) +
          G[1][2] * (pow(b[0][2], 2) - (b[0][0] + b[1][1]) * (b[1][1] + b[2][2])) -
          b[0][2] * (G[0][0] * b[0][1] - G[0][1] * b[0][0] + G[1][0] * b[1][1] - G[1][1] * b[0][1] +
                     G[2][0] * b[1][2] - G[2][1] * b[0][2]) -
          (G[0][0] - G[1][1]) * (b[0][1] * b[0][2] + b[1][2] * (b[1][1] + b[2][2])) -
          (b[0][0] + b[1][1]) * (G[0][1] * b[0][2] - G[0][2] * b[0][1] + G[1][1] * b[1][2] -
                                 G[1][2] * b[1][1] + G[2][1] * b[2][2] - G[2][2] * b[1][2]);
      d_Da23_db[0][1] =
          G[0][2] * (pow(b[1][2], 2) - (b[0][0] + b[1][1]) * (b[0][0] + b[2][2])) +
          G[1][2] * (b[0][1] * (b[0][0] + b[1][1]) + b[0][2] * b[1][2]) +
          b[1][2] * (G[0][0] * b[0][1] - G[0][1] * b[0][0] + G[1][0] * b[1][1] - G[1][1] * b[0][1] +
                     G[2][0] * b[1][2] - G[2][1] * b[0][2]) +
          (G[0][0] - G[1][1]) * (b[0][1] * b[1][2] + b[0][2] * (b[0][0] + b[2][2])) -
          (b[0][0] + b[1][1]) * (G[0][0] * b[0][2] - G[0][2] * b[0][0] + G[1][0] * b[1][2] -
                                 G[1][2] * b[0][1] + G[2][0] * b[2][2] - G[2][2] * b[0][2]);
      d_Da12_db[0][2] =
          G[0][1] * (b[0][1] * b[1][2] + b[0][2] * (b[0][0] + b[2][2])) +
          G[2][1] * (pow(b[0][1], 2) - (b[0][0] + b[2][2]) * (b[1][1] + b[2][2])) -
          b[0][1] * (G[0][0] * b[0][2] - G[0][2] * b[0][0] + G[1][0] * b[1][2] - G[1][2] * b[0][1] +
                     G[2][0] * b[2][2] - G[2][2] * b[0][2]) -
          (G[0][0] - G[2][2]) * (b[0][1] * b[0][2] + b[1][2] * (b[1][1] + b[2][2])) +
          (b[0][0] + b[2][2]) * (G[0][1] * b[0][2] - G[0][2] * b[0][1] + G[1][1] * b[1][2] -
                                 G[1][2] * b[1][1] + G[2][1] * b[2][2] - G[2][2] * b[1][2]);
      d_Da13_db[0][2] =
          -G[0][1] * (b[0][1] * (b[0][0] + b[1][1]) + b[0][2] * b[1][2]) +
          G[2][1] * (b[0][1] * b[0][2] + b[1][2] * (b[1][1] + b[2][2])) -
          b[0][1] * (G[0][0] * b[0][1] - G[0][1] * b[0][0] + G[1][0] * b[1][1] - G[1][1] * b[0][1] +
                     G[2][0] * b[1][2] - G[2][1] * b[0][2]) -
          2 * b[0][2] *
              (G[0][0] * b[0][2] - G[0][2] * b[0][0] + G[1][0] * b[1][2] - G[1][2] * b[0][1] +
               G[2][0] * b[2][2] - G[2][2] * b[0][2]) -
          b[1][2] * (G[0][1] * b[0][2] - G[0][2] * b[0][1] + G[1][1] * b[1][2] - G[1][2] * b[1][1] +
                     G[2][1] * b[2][2] - G[2][2] * b[1][2]) -
          (G[0][0] - G[2][2]) * (pow(b[0][2], 2) - (b[0][0] + b[1][1]) * (b[1][1] + b[2][2]));
      d_Da23_db[0][2] =
          -G[0][1] * (pow(b[1][2], 2) - (b[0][0] + b[1][1]) * (b[0][0] + b[2][2])) -
          G[2][1] * (b[0][1] * b[1][2] + b[0][2] * (b[0][0] + b[2][2])) -
          b[1][2] * (G[0][0] * b[0][2] - G[0][2] * b[0][0] + G[1][0] * b[1][2] - G[1][2] * b[0][1] +
                     G[2][0] * b[2][2] - G[2][2] * b[0][2]) -
          (G[0][0] - G[2][2]) * (b[0][1] * (b[0][0] + b[1][1]) + b[0][2] * b[1][2]) +
          (b[0][0] + b[2][2]) * (G[0][0] * b[0][1] - G[0][1] * b[0][0] + G[1][0] * b[1][1] -
                                 G[1][1] * b[0][1] + G[2][0] * b[1][2] - G[2][1] * b[0][2]);
      d_Da12_db[1][0] =
          -G[0][2] * (b[0][1] * b[1][2] + b[0][2] * (b[0][0] + b[2][2])) +
          G[1][2] * (b[0][1] * b[0][2] + b[1][2] * (b[1][1] + b[2][2])) -
          2 * b[0][1] *
              (G[0][0] * b[0][1] - G[0][1] * b[0][0] + G[1][0] * b[1][1] - G[1][1] * b[0][1] +
               G[2][0] * b[1][2] - G[2][1] * b[0][2]) -
          b[0][2] * (G[0][0] * b[0][2] - G[0][2] * b[0][0] + G[1][0] * b[1][2] - G[1][2] * b[0][1] +
                     G[2][0] * b[2][2] - G[2][2] * b[0][2]) +
          b[1][2] * (G[0][1] * b[0][2] - G[0][2] * b[0][1] + G[1][1] * b[1][2] - G[1][2] * b[1][1] +
                     G[2][1] * b[2][2] - G[2][2] * b[1][2]) -
          (G[0][0] - G[1][1]) * (pow(b[0][1], 2) - (b[0][0] + b[2][2]) * (b[1][1] + b[2][2]));
      d_Da13_db[1][0] =
          G[0][2] * (b[0][1] * (b[0][0] + b[1][1]) + b[0][2] * b[1][2]) +
          G[1][2] * (pow(b[0][2], 2) - (b[0][0] + b[1][1]) * (b[1][1] + b[2][2])) -
          b[0][2] * (G[0][0] * b[0][1] - G[0][1] * b[0][0] + G[1][0] * b[1][1] - G[1][1] * b[0][1] +
                     G[2][0] * b[1][2] - G[2][1] * b[0][2]) -
          (G[0][0] - G[1][1]) * (b[0][1] * b[0][2] + b[1][2] * (b[1][1] + b[2][2])) -
          (b[0][0] + b[1][1]) * (G[0][1] * b[0][2] - G[0][2] * b[0][1] + G[1][1] * b[1][2] -
                                 G[1][2] * b[1][1] + G[2][1] * b[2][2] - G[2][2] * b[1][2]);
      d_Da23_db[1][0] =
          G[0][2] * (pow(b[1][2], 2) - (b[0][0] + b[1][1]) * (b[0][0] + b[2][2])) +
          G[1][2] * (b[0][1] * (b[0][0] + b[1][1]) + b[0][2] * b[1][2]) +
          b[1][2] * (G[0][0] * b[0][1] - G[0][1] * b[0][0] + G[1][0] * b[1][1] - G[1][1] * b[0][1] +
                     G[2][0] * b[1][2] - G[2][1] * b[0][2]) +
          (G[0][0] - G[1][1]) * (b[0][1] * b[1][2] + b[0][2] * (b[0][0] + b[2][2])) -
          (b[0][0] + b[1][1]) * (G[0][0] * b[0][2] - G[0][2] * b[0][0] + G[1][0] * b[1][2] -
                                 G[1][2] * b[0][1] + G[2][0] * b[2][2] - G[2][2] * b[0][2]);
      d_Da12_db[1][1] =
          -G[1][0] * (pow(b[0][1], 2) - (b[0][0] + b[2][2]) * (b[1][1] + b[2][2])) -
          G[1][2] * (b[0][1] * b[1][2] + b[0][2] * (b[0][0] + b[2][2])) -
          b[1][2] * (G[0][0] * b[0][2] - G[0][2] * b[0][0] + G[1][0] * b[1][2] - G[1][2] * b[0][1] +
                     G[2][0] * b[2][2] - G[2][2] * b[0][2]) +
          (b[0][0] + b[2][2]) * (G[0][0] * b[0][1] - G[0][1] * b[0][0] + G[1][0] * b[1][1] -
                                 G[1][1] * b[0][1] + G[2][0] * b[1][2] - G[2][1] * b[0][2]);
      d_Da13_db[1][1] = -G[1][0] * (b[0][1] * b[0][2] + b[1][2] * (b[1][1] + b[2][2])) +
                        G[1][2] * (b[0][1] * (b[0][0] + b[1][1]) + b[0][2] * b[1][2]) -
                        b[0][1] * (G[0][1] * b[0][2] - G[0][2] * b[0][1] + G[1][1] * b[1][2] -
                                   G[1][2] * b[1][1] + G[2][1] * b[2][2] - G[2][2] * b[1][2]) -
                        b[1][2] * (G[0][0] * b[0][1] - G[0][1] * b[0][0] + G[1][0] * b[1][1] -
                                   G[1][1] * b[0][1] + G[2][0] * b[1][2] - G[2][1] * b[0][2]) +
                        (b[0][0] + 2 * b[1][1] + b[2][2]) *
                            (G[0][0] * b[0][2] - G[0][2] * b[0][0] + G[1][0] * b[1][2] -
                             G[1][2] * b[0][1] + G[2][0] * b[2][2] - G[2][2] * b[0][2]);
      d_Da23_db[1][1] =
          G[1][0] * (b[0][1] * b[1][2] + b[0][2] * (b[0][0] + b[2][2])) +
          G[1][2] * (pow(b[1][2], 2) - (b[0][0] + b[1][1]) * (b[0][0] + b[2][2])) -
          b[0][1] * (G[0][0] * b[0][2] - G[0][2] * b[0][0] + G[1][0] * b[1][2] - G[1][2] * b[0][1] +
                     G[2][0] * b[2][2] - G[2][2] * b[0][2]) +
          (b[0][0] + b[2][2]) * (G[0][1] * b[0][2] - G[0][2] * b[0][1] + G[1][1] * b[1][2] -
                                 G[1][2] * b[1][1] + G[2][1] * b[2][2] - G[2][2] * b[1][2]);
      d_Da12_db[1][2] =
          -G[1][0] * (b[0][1] * b[0][2] + b[1][2] * (b[1][1] + b[2][2])) -
          G[2][0] * (pow(b[0][1], 2) - (b[0][0] + b[2][2]) * (b[1][1] + b[2][2])) +
          b[0][1] * (G[0][1] * b[0][2] - G[0][2] * b[0][1] + G[1][1] * b[1][2] - G[1][2] * b[1][1] +
                     G[2][1] * b[2][2] - G[2][2] * b[1][2]) +
          (G[1][1] - G[2][2]) * (b[0][1] * b[1][2] + b[0][2] * (b[0][0] + b[2][2])) -
          (b[1][1] + b[2][2]) * (G[0][0] * b[0][2] - G[0][2] * b[0][0] + G[1][0] * b[1][2] -
                                 G[1][2] * b[0][1] + G[2][0] * b[2][2] - G[2][2] * b[0][2]);
      d_Da13_db[1][2] =
          -G[1][0] * (pow(b[0][2], 2) - (b[0][0] + b[1][1]) * (b[1][1] + b[2][2])) -
          G[2][0] * (b[0][1] * b[0][2] + b[1][2] * (b[1][1] + b[2][2])) -
          b[0][2] * (G[0][1] * b[0][2] - G[0][2] * b[0][1] + G[1][1] * b[1][2] - G[1][2] * b[1][1] +
                     G[2][1] * b[2][2] - G[2][2] * b[1][2]) -
          (G[1][1] - G[2][2]) * (b[0][1] * (b[0][0] + b[1][1]) + b[0][2] * b[1][2]) -
          (b[1][1] + b[2][2]) * (G[0][0] * b[0][1] - G[0][1] * b[0][0] + G[1][0] * b[1][1] -
                                 G[1][1] * b[0][1] + G[2][0] * b[1][2] - G[2][1] * b[0][2]);
      d_Da23_db[1][2] =
          -G[1][0] * (b[0][1] * (b[0][0] + b[1][1]) + b[0][2] * b[1][2]) +
          G[2][0] * (b[0][1] * b[1][2] + b[0][2] * (b[0][0] + b[2][2])) +
          b[0][1] * (G[0][0] * b[0][1] - G[0][1] * b[0][0] + G[1][0] * b[1][1] - G[1][1] * b[0][1] +
                     G[2][0] * b[1][2] - G[2][1] * b[0][2]) -
          b[0][2] * (G[0][0] * b[0][2] - G[0][2] * b[0][0] + G[1][0] * b[1][2] - G[1][2] * b[0][1] +
                     G[2][0] * b[2][2] - G[2][2] * b[0][2]) -
          2 * b[1][2] *
              (G[0][1] * b[0][2] - G[0][2] * b[0][1] + G[1][1] * b[1][2] - G[1][2] * b[1][1] +
               G[2][1] * b[2][2] - G[2][2] * b[1][2]) -
          (G[1][1] - G[2][2]) * (pow(b[1][2], 2) - (b[0][0] + b[1][1]) * (b[0][0] + b[2][2]));
      d_Da12_db[2][0] =
          G[0][1] * (b[0][1] * b[1][2] + b[0][2] * (b[0][0] + b[2][2])) +
          G[2][1] * (pow(b[0][1], 2) - (b[0][0] + b[2][2]) * (b[1][1] + b[2][2])) -
          b[0][1] * (G[0][0] * b[0][2] - G[0][2] * b[0][0] + G[1][0] * b[1][2] - G[1][2] * b[0][1] +
                     G[2][0] * b[2][2] - G[2][2] * b[0][2]) -
          (G[0][0] - G[2][2]) * (b[0][1] * b[0][2] + b[1][2] * (b[1][1] + b[2][2])) +
          (b[0][0] + b[2][2]) * (G[0][1] * b[0][2] - G[0][2] * b[0][1] + G[1][1] * b[1][2] -
                                 G[1][2] * b[1][1] + G[2][1] * b[2][2] - G[2][2] * b[1][2]);
      d_Da13_db[2][0] =
          -G[0][1] * (b[0][1] * (b[0][0] + b[1][1]) + b[0][2] * b[1][2]) +
          G[2][1] * (b[0][1] * b[0][2] + b[1][2] * (b[1][1] + b[2][2])) -
          b[0][1] * (G[0][0] * b[0][1] - G[0][1] * b[0][0] + G[1][0] * b[1][1] - G[1][1] * b[0][1] +
                     G[2][0] * b[1][2] - G[2][1] * b[0][2]) -
          2 * b[0][2] *
              (G[0][0] * b[0][2] - G[0][2] * b[0][0] + G[1][0] * b[1][2] - G[1][2] * b[0][1] +
               G[2][0] * b[2][2] - G[2][2] * b[0][2]) -
          b[1][2] * (G[0][1] * b[0][2] - G[0][2] * b[0][1] + G[1][1] * b[1][2] - G[1][2] * b[1][1] +
                     G[2][1] * b[2][2] - G[2][2] * b[1][2]) -
          (G[0][0] - G[2][2]) * (pow(b[0][2], 2) - (b[0][0] + b[1][1]) * (b[1][1] + b[2][2]));

      d_Da23_db[2][0] =
          -G[0][1] * (pow(b[1][2], 2) - (b[0][0] + b[1][1]) * (b[0][0] + b[2][2])) -
          G[2][1] * (b[0][1] * b[1][2] + b[0][2] * (b[0][0] + b[2][2])) -
          b[1][2] * (G[0][0] * b[0][2] - G[0][2] * b[0][0] + G[1][0] * b[1][2] - G[1][2] * b[0][1] +
                     G[2][0] * b[2][2] - G[2][2] * b[0][2]) -
          (G[0][0] - G[2][2]) * (b[0][1] * (b[0][0] + b[1][1]) + b[0][2] * b[1][2]) +
          (b[0][0] + b[2][2]) * (G[0][0] * b[0][1] - G[0][1] * b[0][0] + G[1][0] * b[1][1] -
                                 G[1][1] * b[0][1] + G[2][0] * b[1][2] - G[2][1] * b[0][2]);
      d_Da12_db[2][1] =
          -G[1][0] * (b[0][1] * b[0][2] + b[1][2] * (b[1][1] + b[2][2])) -
          G[2][0] * (pow(b[0][1], 2) - (b[0][0] + b[2][2]) * (b[1][1] + b[2][2])) +
          b[0][1] * (G[0][1] * b[0][2] - G[0][2] * b[0][1] + G[1][1] * b[1][2] - G[1][2] * b[1][1] +
                     G[2][1] * b[2][2] - G[2][2] * b[1][2]) +
          (G[1][1] - G[2][2]) * (b[0][1] * b[1][2] + b[0][2] * (b[0][0] + b[2][2])) -
          (b[1][1] + b[2][2]) * (G[0][0] * b[0][2] - G[0][2] * b[0][0] + G[1][0] * b[1][2] -
                                 G[1][2] * b[0][1] + G[2][0] * b[2][2] - G[2][2] * b[0][2]);
      d_Da13_db[2][1] =
          -G[1][0] * (pow(b[0][2], 2) - (b[0][0] + b[1][1]) * (b[1][1] + b[2][2])) -
          G[2][0] * (b[0][1] * b[0][2] + b[1][2] * (b[1][1] + b[2][2])) -
          b[0][2] * (G[0][1] * b[0][2] - G[0][2] * b[0][1] + G[1][1] * b[1][2] - G[1][2] * b[1][1] +
                     G[2][1] * b[2][2] - G[2][2] * b[1][2]) -
          (G[1][1] - G[2][2]) * (b[0][1] * (b[0][0] + b[1][1]) + b[0][2] * b[1][2]) -
          (b[1][1] + b[2][2]) * (G[0][0] * b[0][1] - G[0][1] * b[0][0] + G[1][0] * b[1][1] -
                                 G[1][1] * b[0][1] + G[2][0] * b[1][2] - G[2][1] * b[0][2]);
      d_Da23_db[2][1] =
          -G[1][0] * (b[0][1] * (b[0][0] + b[1][1]) + b[0][2] * b[1][2]) +
          G[2][0] * (b[0][1] * b[1][2] + b[0][2] * (b[0][0] + b[2][2])) +
          b[0][1] * (G[0][0] * b[0][1] - G[0][1] * b[0][0] + G[1][0] * b[1][1] - G[1][1] * b[0][1] +
                     G[2][0] * b[1][2] - G[2][1] * b[0][2]) -
          b[0][2] * (G[0][0] * b[0][2] - G[0][2] * b[0][0] + G[1][0] * b[1][2] - G[1][2] * b[0][1] +
                     G[2][0] * b[2][2] - G[2][2] * b[0][2]) -
          2 * b[1][2] *
              (G[0][1] * b[0][2] - G[0][2] * b[0][1] + G[1][1] * b[1][2] - G[1][2] * b[1][1] +
               G[2][1] * b[2][2] - G[2][2] * b[1][2]) -
          (G[1][1] - G[2][2]) * (pow(b[1][2], 2) - (b[0][0] + b[1][1]) * (b[0][0] + b[2][2]));
      d_Da12_db[2][2] = -G[2][0] * (b[0][1] * b[0][2] + b[1][2] * (b[1][1] + b[2][2])) +
                        G[2][1] * (b[0][1] * b[1][2] + b[0][2] * (b[0][0] + b[2][2])) +
                        b[0][2] * (G[0][1] * b[0][2] - G[0][2] * b[0][1] + G[1][1] * b[1][2] -
                                   G[1][2] * b[1][1] + G[2][1] * b[2][2] - G[2][2] * b[1][2]) -
                        b[1][2] * (G[0][0] * b[0][2] - G[0][2] * b[0][0] + G[1][0] * b[1][2] -
                                   G[1][2] * b[0][1] + G[2][0] * b[2][2] - G[2][2] * b[0][2]) +
                        (b[0][0] + b[1][1] + 2 * b[2][2]) *
                            (G[0][0] * b[0][1] - G[0][1] * b[0][0] + G[1][0] * b[1][1] -
                             G[1][1] * b[0][1] + G[2][0] * b[1][2] - G[2][1] * b[0][2]);
      d_Da13_db[2][2] =
          -G[2][0] * (pow(b[0][2], 2) - (b[0][0] + b[1][1]) * (b[1][1] + b[2][2])) -
          G[2][1] * (b[0][1] * (b[0][0] + b[1][1]) + b[0][2] * b[1][2]) -
          b[1][2] * (G[0][0] * b[0][1] - G[0][1] * b[0][0] + G[1][0] * b[1][1] - G[1][1] * b[0][1] +
                     G[2][0] * b[1][2] - G[2][1] * b[0][2]) +
          (b[0][0] + b[1][1]) * (G[0][0] * b[0][2] - G[0][2] * b[0][0] + G[1][0] * b[1][2] -
                                 G[1][2] * b[0][1] + G[2][0] * b[2][2] - G[2][2] * b[0][2]);
      d_Da23_db[2][2] =
          -G[2][0] * (b[0][1] * (b[0][0] + b[1][1]) + b[0][2] * b[1][2]) -
          G[2][1] * (pow(b[1][2], 2) - (b[0][0] + b[1][1]) * (b[0][0] + b[2][2])) +
          b[0][2] * (G[0][0] * b[0][1] - G[0][1] * b[0][0] + G[1][0] * b[1][1] - G[1][1] * b[0][1] +
                     G[2][0] * b[1][2] - G[2][1] * b[0][2]) +
          (b[0][0] + b[1][1]) * (G[0][1] * b[0][2] - G[0][2] * b[0][1] + G[1][1] * b[1][2] -
                                 G[1][2] * b[1][1] + G[2][1] * b[2][2] - G[2][2] * b[1][2]);

      for (int i = 0; i < VIM; i++) {
        for (int j = 0; j < VIM; j++) {

          dbl da12 = -a12 * invD * d_D_db[i][j] + invD * d_Da12_db[i][j];
          dbl da13 = -a13 * invD * d_D_db[i][j] + invD * d_Da13_db[i][j];
          dbl da23 = -a23 * invD * d_D_db[i][j] + invD * d_Da23_db[i][j];

          dbl da[DIM][DIM] = {{0., da12, da13}, {-da12, 0., da23}, {-da13, -da23, 0.}};

          dbl db[DIM][DIM];
          for (int p = 0; p < VIM; p++) {
            for (int q = 0; q < VIM; q++) {
              db[p][q] = delta(i, p) * delta(j, q) | delta(i, q) * delta(j, p);
            }
          }
          dbl a_dot_db[DIM][DIM];
          dbl da_dot_b[DIM][DIM];

          tensor_dot(da, b, da_dot_b, VIM);
          tensor_dot(a, db, a_dot_db, VIM);
          for (int p = 0; p < VIM; p++) {
            for (int q = 0; q < VIM; q++) {
              d_a_dot_b_db[i][j][p][q] = da_dot_b[p][q] + a_dot_db[p][q];
            }
          }
        }
      }
    }
  }
}

int sqrt_conf_source(int mode,
                     dbl b[DIM][DIM],
                     dbl source_term[DIM][DIM],
                     dbl d_source_term_db[DIM][DIM][DIM][DIM]) {
  dbl binv[DIM][DIM];
  dbl d_binv_db[DIM][DIM][DIM][DIM];
  if (VIM == 2) {
    dbl det = b[0][0] * b[1][1] - b[0][1] * b[0][1] + 1e-16;
    binv[0][0] = b[1][1] / det;
    binv[0][1] = -b[0][1] / det;
    binv[1][0] = -b[0][1] / det;
    binv[1][1] = b[0][0] / det;

    for (int p = 0; p < VIM; p++) {
      for (int q = 0; q < VIM; q++) {
        dbl ddet = delta(p, 0) * delta(q, 0) * b[1][1] + b[0][0] * delta(p, 1) * delta(q, 1) -
                   2.0 * (delta(p, 0) * delta(q, 1) + delta(p, 1) * delta(q, 0)) * b[0][1];
        d_binv_db[0][0][p][q] = (det * delta(p, 1) * delta(q, 1) - ddet * b[1][1]) / (det * det);
        d_binv_db[0][1][p][q] = (-det * delta(p, 0) * delta(q, 1) + ddet * b[0][1]) / (det * det);
        d_binv_db[1][0][p][q] = (-det * delta(p, 0) * delta(q, 1) + ddet * b[1][0]) / (det * det);
        d_binv_db[1][1][p][q] = (det * delta(p, 0) * delta(q, 0) - ddet * b[0][0]) / (det * det);
      }
    }
  } else if (VIM == 3) {
    dbl det = b[0][0] * (b[1][1] * b[2][2] - b[1][2] * b[2][1]) -
              b[0][1] * (b[1][0] * b[2][2] - b[2][0] * b[1][2]) +
              b[0][2] * (b[1][0] * b[2][1] - b[2][0] * b[1][1]) + 1e-16;

    binv[0][0] = (b[1][1] * b[2][2] - b[2][1] * b[1][2]) / (det);

    binv[0][1] = -(b[0][1] * b[2][2] - b[2][1] * b[0][2]) / (det);

    binv[0][2] = (b[0][1] * b[1][2] - b[1][1] * b[0][2]) / (det);

    binv[1][0] = -(b[1][0] * b[2][2] - b[2][0] * b[1][2]) / (det);

    binv[1][1] = (b[0][0] * b[2][2] - b[2][0] * b[0][2]) / (det);

    binv[1][2] = -(b[0][0] * b[1][2] - b[1][0] * b[0][2]) / (det);

    binv[2][0] = (b[1][0] * b[2][1] - b[1][1] * b[2][0]) / (det);

    binv[2][1] = -(b[0][0] * b[2][1] - b[2][0] * b[0][1]) / (det);

    binv[2][2] = (b[0][0] * b[1][1] - b[1][0] * b[0][1]) / (det);

    for (int p = 0; p < VIM; p++) {
      for (int q = 0; q < VIM; q++) {
        dbl db[DIM][DIM] = {{0.}};
        db[p][q] = 1.0;
        db[q][p] = 1.0;
        dbl ddet = db[0][0] * (b[1][1] * b[2][2] - b[1][2] * b[2][1]) +
                   b[0][0] * (db[1][1] * b[2][2] - db[1][2] * b[2][1] + b[1][1] * db[2][2] -
                              b[1][2] * db[2][1]) -
                   db[0][1] * (b[1][0] * b[2][2] - b[2][0] * b[1][2]) -
                   b[0][1] * (db[1][0] * b[2][2] - db[2][0] * b[1][2] + b[1][0] * db[2][2] -
                              b[2][0] * db[1][2]) +
                   db[0][2] * (b[1][0] * b[2][1] - b[2][0] * b[1][1]) +
                   b[0][2] * (db[1][0] * b[2][1] - db[2][0] * b[1][1] + b[1][0] * db[2][1] -
                              b[2][0] * db[1][1]);

        d_binv_db[0][0][p][q] =
            (db[1][1] * b[2][2] - db[2][1] * b[1][2] + b[1][1] * db[2][2] - b[2][1] * db[1][2]) /
            (det);
        d_binv_db[0][0][p][q] += (b[1][1] * b[2][2] - b[2][1] * b[1][2]) * -ddet / (det * det);

        d_binv_db[0][1][p][q] =
            -(db[0][1] * b[2][2] - db[2][1] * b[0][2] + b[0][1] * db[2][2] - b[2][1] * db[0][2]) /
            (det);
        d_binv_db[0][1][p][q] += -(b[0][1] * b[2][2] - b[2][1] * b[0][2]) * -ddet / (det * det);

        d_binv_db[0][2][p][q] =
            (db[0][1] * b[1][2] - db[1][1] * b[0][2] + b[0][1] * db[1][2] - b[1][1] * db[0][2]) /
            (det);
        d_binv_db[0][2][p][q] += (b[0][1] * b[1][2] - b[1][1] * b[0][2]) * -ddet / (det * det);

        d_binv_db[1][0][p][q] =
            -(db[1][0] * b[2][2] - db[2][0] * b[1][2] + b[1][0] * db[2][2] - b[2][0] * db[1][2]) /
            (det);
        d_binv_db[1][0][p][q] += -(b[1][0] * b[2][2] - b[2][0] * b[1][2]) * -ddet / (det * det);

        d_binv_db[1][1][p][q] =
            (db[0][0] * b[2][2] - db[2][0] * b[0][2] + b[0][0] * db[2][2] - b[2][0] * db[0][2]) /
            (det);
        d_binv_db[1][1][p][q] += (b[0][0] * b[2][2] - b[2][0] * b[0][2]) * -ddet / (det * det);

        d_binv_db[1][2][p][q] =
            -(db[0][0] * b[1][2] - db[1][0] * b[0][2] + b[0][0] * db[1][2] - b[1][0] * db[0][2]) /
            (det);
        d_binv_db[1][2][p][q] += -(b[0][0] * b[1][2] - b[1][0] * b[0][2]) * -ddet / (det * det);

        d_binv_db[2][0][p][q] =
            (db[1][0] * b[2][1] - db[1][1] * b[2][0] + b[1][0] * db[2][1] - b[1][1] * db[2][0]) /
            (det);
        d_binv_db[2][0][p][q] += (b[1][0] * b[2][1] - b[1][1] * b[2][0]) * -ddet / (det * det);

        d_binv_db[2][1][p][q] =
            -(db[0][0] * b[2][1] - db[2][0] * b[0][1] + b[0][0] * db[2][1] - b[2][0] * db[0][1]) /
            (det);
        d_binv_db[2][1][p][q] += -(b[0][0] * b[2][1] - b[2][0] * b[0][1]) * -ddet / (det * det);

        d_binv_db[2][2][p][q] =
            (db[0][0] * b[1][1] - db[1][0] * b[0][1] + b[0][0] * db[1][1] - b[1][0] * db[0][1]) /
            (det);
        d_binv_db[2][2][p][q] += (b[0][0] * b[1][1] - b[1][0] * b[0][1]) * -ddet / (det * det);
      }
    }
  } else {
    GOMA_EH(GOMA_ERROR, "Unknown VIM = %d for SQRT conformation tensor", VIM);
  }

  switch (vn->ConstitutiveEquation) {
  case WHITE_METZNER:
  case OLDROYDB: {
    for (int ii = 0; ii < VIM; ii++) {
      for (int jj = 0; jj < VIM; jj++) {
        source_term[ii][jj] = -0.5 * (binv[ii][jj] - b[ii][jj]);
      }
    }
    if (af->Assemble_Jacobian) {
      for (int ii = 0; ii < VIM; ii++) {
        for (int jj = 0; jj < VIM; jj++) {
          for (int p = 0; p < VIM; p++) {
            for (int q = 0; q < VIM; q++) {
              d_source_term_db[ii][jj][p][q] =
                  -0.5 * (d_binv_db[ii][jj][p][q] - delta(ii, p) * delta(jj, q));
            }
          }
        }
      }
    }
  } break;
  case PTT: {

    dbl d_trace_db[DIM][DIM] = {{0.0}};

    dbl trace = 0;
    for (int i = 0; i < VIM; i++) {
      for (int j = 0; j < VIM; j++) {
        trace += b[i][j] * b[i][j];
      }
    }

    if (af->Assemble_Jacobian) {
      for (int p = 0; p < VIM; p++) {
        for (int q = 0; q < VIM; q++) {
          d_trace_db[p][q] = 0.0;
          for (int i = 0; i < VIM; i++) {
            for (int j = 0; j < VIM; j++) {
              d_trace_db[p][q] +=
                  2.0 * b[i][j] * (delta(p, i) * delta(q, j) | delta(p, j) * delta(q, i));
            }
          }
        }
      }
    }

    dbl Z = 1.0;
    dbl dZ_dtrace = 0;

    // PTT exponent
    eps = ve[mode]->eps;

    if (vn->ptt_type == PTT_LINEAR) {
      Z = 1 + eps * (trace - (double)VIM);
      dZ_dtrace = eps;
    } else if (vn->ptt_type == PTT_EXPONENTIAL) {
      const double exp_max = 700;
      double inner = eps * (trace - (double)VIM);
      if ((inner > exp_max) || (inner < -exp_max)) {
        GOMA_WH_MANY(GOMA_ERROR, "Exponential overflow in PTT_EXPONENTIAL");
        return GOMA_ERROR;
      }
      Z = exp(eps * (trace - (double)VIM));
      dZ_dtrace = eps * Z;
    } else {
      GOMA_EH(GOMA_ERROR, "Unrecognized PTT Form %d", vn->ptt_type);
    }

    for (int ii = 0; ii < VIM; ii++) {
      for (int jj = 0; jj < VIM; jj++) {
        source_term[ii][jj] = -0.5 * Z * (binv[ii][jj] - b[ii][jj]);
      }
    }
    if (af->Assemble_Jacobian) {
      for (int ii = 0; ii < VIM; ii++) {
        for (int jj = 0; jj < VIM; jj++) {
          for (int p = 0; p < VIM; p++) {
            for (int q = 0; q < VIM; q++) {
              d_source_term_db[ii][jj][p][q] =
                  -0.5 * Z * (d_binv_db[ii][jj][p][q] - delta(ii, p) * delta(jj, q)) -
                  0.5 * dZ_dtrace * d_trace_db[p][q] * (binv[ii][jj] - b[ii][jj]);
            }
          }
        }
      }
    }
  } break;
  default:
    GOMA_EH(GOMA_ERROR, "Unknown Constitutive equation form for SQRT_CONF");
    break;
  }

  return GOMA_SUCCESS;
}

int assemble_stress_sqrt_conf(dbl tt, /* parameter to vary time integration from
                                       * explicit (tt = 1) to implicit (tt = 0) */
                              dbl dt, /* current time step size */
                              PG_DATA *pg_data) {
  int dim, p, q, r, w;

  int eqn, var;
  int peqn, pvar;
  int evss_gradv = 0;

  int i, j, status, mode;
  dbl v[DIM];      /* Velocity field. */
  dbl x_dot[DIM];  /* current position field derivative wrt time. */
  dbl h3;          /* Volume element (scale factors). */
  dbl dh3dmesh_pj; /* Sensitivity to (p,j) mesh dof. */

  dbl grad_v[DIM][DIM];
  dbl gamma[DIM][DIM]; /* Shear-rate tensor based on velocity */
  dbl det_J;           /* determinant of element Jacobian */

  dbl d_det_J_dmesh_pj; /* for specific (p,j) mesh dof */

  dbl mass; /* For terms and their derivatives */
  dbl mass_a, mass_b;
  dbl advection;
  dbl advection_a, advection_b, advection_c, advection_d;
  dbl diffusion;
  dbl source;
  dbl source_a = 0, source_b = 0, source_c = 0;
  int err;
  dbl alpha = 0;  /* This is the Geisekus mobility parameter */
  dbl lambda = 0; /* polymer relaxation constant */
  double xi;
  double d_xi_dF[MDE];
  dbl eps = 0; /* This is the PTT elongation parameter */
  double d_eps_dF[MDE];
  /*
   *
   * Note how carefully we avoid refering to d(phi[i])/dx[j] and refer instead
   * to the j-th component of grad_phi[j][i] so that this vector can be loaded
   * up with components that may be different in non Cartesian coordinate
   * systems.
   *
   * We will, however, insist on *orthogonal* coordinate systems, even if we
   * might permit them to be curvilinear.
   *
   * Assume all components of velocity are interpolated with the same kind
   * of basis function.
   */

  /*
   * Petrov-Galerkin weighting functions for i-th and ab-th stress residuals
   * and some of their derivatives...
   */

  dbl wt_func;

  /*
   * Interpolation functions for variables and some of their derivatives.
   */

  dbl phi_j;
  dbl wt;

  /* Variables for stress */

  int R_s[MAX_MODES][DIM][DIM];
  int v_s[MAX_MODES][DIM][DIM];
  int v_g[DIM][DIM];

  dbl b[DIM][DIM];     /* stress tensor */
  dbl b_dot[DIM][DIM]; /* stress tensor from last time step */
  dbl grad_b[DIM][DIM][DIM];
  dbl d_grad_b_dmesh[DIM][DIM][DIM][DIM]
                    [MDE]; /* derivative of grad of stress tensor for mode ve_mode */

  dbl g[DIM][DIM]; /* velocity gradient tensor */

  /* dot product tensors */

  dbl b_dot_g[DIM][DIM];

  /* polymer viscosity and derivatives */
  POLYMER_TIME_CONST_DEPENDENCE_STRUCT d_lam_struct;
  POLYMER_TIME_CONST_DEPENDENCE_STRUCT *d_lam = &d_lam_struct;

  const bool saramitoEnabled =
      (vn->ConstitutiveEquation == SARAMITO_OLDROYDB || vn->ConstitutiveEquation == SARAMITO_PTT ||
       vn->ConstitutiveEquation == SARAMITO_GIESEKUS);

  if (saramitoEnabled) {
    GOMA_EH(GOMA_ERROR, "Saramito not available for SQRT_CONF");
  }

  /*  shift function */
  dbl at = 0.0;
  dbl d_at_dT[MDE];
  dbl wlf_denom;

  /* advective terms are precalculated */
  dbl v_dot_del_b[DIM][DIM];
  dbl x_dot_del_b[DIM][DIM];

  dbl d_xdotdels_dm;

  dbl d_vdotdels_dm;

  /* SUPG variables */
  dbl supg = 0;

  if (vn->evssModel == EVSS_GRADV) {
    evss_gradv = 1;
  }

  status = 0;

  eqn = R_STRESS11;

  /*
   * Bail out fast if there's nothing to do...
   */

  if (!pd->e[pg->imtrx][eqn]) {
    return (status);
  }

  /*
   * Unpack variables from structures for local convenience...
   */

  dim = pd->Num_Dim;

  wt = fv->wt;

  det_J = bf[eqn]->detJ; /* Really, ought to be mesh eqn. */

  h3 = fv->h3; /* Differential volume element (scales). */

  /* load eqn and variable number in tensor form */
  (void)stress_eqn_pointer(v_s);
  (void)stress_eqn_pointer(R_s);

  v_g[0][0] = VELOCITY_GRADIENT11;
  v_g[0][1] = VELOCITY_GRADIENT12;
  v_g[1][0] = VELOCITY_GRADIENT21;
  v_g[1][1] = VELOCITY_GRADIENT22;
  v_g[0][2] = VELOCITY_GRADIENT13;
  v_g[1][2] = VELOCITY_GRADIENT23;
  v_g[2][0] = VELOCITY_GRADIENT31;
  v_g[2][1] = VELOCITY_GRADIENT32;
  v_g[2][2] = VELOCITY_GRADIENT33;

  /*
   * Field variables...
   */
  for (int a = 0; a < WIM; a++) {
    v[a] = fv->v[a];

    /* note, these are zero for steady calculations */
    x_dot[a] = 0.0;
    if (pd->TimeIntegration != STEADY && pd->gv[MESH_DISPLACEMENT1 + a]) {
      x_dot[a] = fv_dot->x[a];
    }
  }

  /*
   * In Cartesian coordinates, this velocity gradient tensor will
   * have components that are...
   *
   * 			grad_v[a][b] = d v_b
   *				       -----
   *				       d x_a
   */

  for (int a = 0; a < VIM; a++) {
    for (int b = 0; b < VIM; b++) {
      grad_v[a][b] = fv->grad_v[a][b];
    }
  }

  /* load up shearrate tensor based on velocity */
  for (int a = 0; a < VIM; a++) {
    for (int b = 0; b < VIM; b++) {
      gamma[a][b] = grad_v[a][b] + grad_v[b][a];
    }
  }

  for (int a = 0; a < VIM; a++) {
    for (int b = 0; b < VIM; b++) {
      if (evss_gradv) {
        g[a][b] = fv->grad_v[a][b];
      } else {
        g[a][b] = fv->G[a][b];
      }
    }
  }

  if (vn->wt_funcModel == GALERKIN) {
    supg = 0.;
  } else if (vn->wt_funcModel == SUPG) {
    supg = vn->wt_func;
  }

  SUPG_terms supg_terms;
  if (supg != 0.) {
    supg_tau(&supg_terms, dim, 1e-14, pg_data, dt, TRUE, eqn);
  }
  /* end Petrov-Galerkin addition */
  dbl yzbeta_factor = 0.0;
  dbl beta[2] = {1.0, 2.0};
  if (vn->shockcaptureModel == SC_YZBETA) {
    yzbeta_factor = vn->shockcapture;
  } else if (vn->shockcaptureModel != SC_NONE) {
    GOMA_EH(GOMA_ERROR, "Unknown shock capture model, only YZBETA supported for SQRT_CONF");
  }

  /*  shift factor  */
  if (pd->gv[TEMPERATURE]) {
    if (vn->shiftModel == CONSTANT) {
      at = vn->shift[0];
      for (j = 0; j < ei[pg->imtrx]->dof[TEMPERATURE]; j++) {
        d_at_dT[j] = 0.;
      }
    } else if (vn->shiftModel == MODIFIED_WLF) {
      wlf_denom = vn->shift[1] + fv->T - mp->reference[TEMPERATURE];
      if (wlf_denom != 0.) {
        at = exp(vn->shift[0] * (mp->reference[TEMPERATURE] - fv->T) / wlf_denom);
        for (j = 0; j < ei[pg->imtrx]->dof[TEMPERATURE]; j++) {
          d_at_dT[j] =
              -at * vn->shift[0] * vn->shift[1] / (wlf_denom * wlf_denom) * bf[TEMPERATURE]->phi[j];
        }
      } else {
        at = 1.;
      }
      for (j = 0; j < ei[pg->imtrx]->dof[TEMPERATURE]; j++) {
        d_at_dT[j] = 0.;
      }
    }
  } else {
    at = 1.;
  }

  /* Begin loop over modes */
  for (mode = 0; mode < vn->modes; mode++) {

    load_modal_pointers(mode, tt, dt, b, b_dot, grad_b, d_grad_b_dmesh);

    /* precalculate advective terms of form (v dot del tensor)*/

    /*
     * Stress tensor...(Note "anti-BSL" sign convention on deviatoric stress)
     */
    for (int ii = 0; ii < VIM; ii++) {
      for (int jj = 0; jj < VIM; jj++) {
        v_dot_del_b[ii][jj] = 0.;
        x_dot_del_b[ii][jj] = 0.;
        for (q = 0; q < WIM; q++) {
          v_dot_del_b[ii][jj] += v[q] * grad_b[q][ii][jj];
          x_dot_del_b[ii][jj] += x_dot[q] * grad_b[q][ii][jj];
        }
      }
    }

    if (saramitoEnabled == TRUE) {
      GOMA_EH(GOMA_ERROR, "Saramito not enabled sqrt");
    }

    double d_alpha_dF[MDE];
    /* get Geisekus mobility parameter */
    if (ve[mode]->alphaModel == CONSTANT) {
      alpha = ve[mode]->alpha;
    } else if (ls != NULL && ve[mode]->alphaModel == VE_LEVEL_SET) {
      double pos_alpha = ve[mode]->pos_ls.alpha;
      double neg_alpha = ve[mode]->alpha;
      double width = ls->Length_Scale;
      err = level_set_property(neg_alpha, pos_alpha, width, &alpha, d_alpha_dF);
      GOMA_EH(err, "level_set_property() failed for mobility parameter.");
    } else {
      GOMA_EH(GOMA_ERROR, "Unknown mobility parameter model");
    }

    /* get time constant */
    lambda = polymer_time_const(ve[mode]->time_const_st, gamma, d_lam);

    xi = 0;
    if (ve[mode]->xiModel == CONSTANT) {
      xi = ve[mode]->xi;
    } else if (ls != NULL && ve[mode]->xiModel == VE_LEVEL_SET) {
      double pos_xi = ve[mode]->pos_ls.xi;
      double neg_xi = ve[mode]->xi;
      double width = ls->Length_Scale;
      err = level_set_property(neg_xi, pos_xi, width, &xi, d_xi_dF);
      GOMA_EH(err, "level_set_property() failed for ptt xi parameter.");
    } else {
      GOMA_EH(GOMA_ERROR, "Unknown PTT Xi parameter model");
    }

    if (DOUBLE_NONZERO(xi)) {
      GOMA_EH(GOMA_ERROR, "PTT Xi parameter currently required to be 0 for SQRT_CONF");
    }

    if (ve[mode]->epsModel == CONSTANT) {
      eps = ve[mode]->eps;
    } else if (ls != NULL && ve[mode]->epsModel == VE_LEVEL_SET) {
      double pos_eps = ve[mode]->pos_ls.eps;
      double neg_eps = ve[mode]->eps;
      double width = ls->Length_Scale;
      err = level_set_property(neg_eps, pos_eps, width, &eps, d_eps_dF);
      GOMA_EH(err, "level_set_property() failed for ptt epsilon parameter.");
    } else {
      GOMA_EH(GOMA_ERROR, "Unknown PTT Epsilon parameter model");
    }

    (void)tensor_dot(b, g, b_dot_g, VIM);

    dbl a_dot_b[DIM][DIM];
    dbl d_a_dot_b_db[DIM][DIM][DIM][DIM];
    dbl d_a_dot_b_dG[DIM][DIM][DIM][DIM];

    compute_a_dot_b(b, g, a_dot_b, d_a_dot_b_db, d_a_dot_b_dG);

    dbl source_term[DIM][DIM];
    dbl d_source_term_db[DIM][DIM][DIM][DIM];
    goma_error err = sqrt_conf_source(mode, b, source_term, d_source_term_db);
    if (err) {
      return 1;
    }
    /*
     * Residuals_________________________________________________________________
     */

    if (af->Assemble_Residual) {
      /*
       * Assemble each component "ab" of the polymer stress equation...
       */
      for (int ii = 0; ii < VIM; ii++) {
        for (int jj = 0; jj < VIM; jj++) {

          if (ii <= jj) /* since the stress tensor is symmetric, only assemble the upper half */
          {
            eqn = R_s[mode][ii][jj];

            /*
             * In the element, there will be contributions to this many equations
             * based on the number of degrees of freedom...
             */

            for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
              wt_func = bf[eqn]->phi[i];
              /* add Petrov-Galerkin terms as necessary */
              if (supg != 0.) {
                for (w = 0; w < dim; w++) {
                  wt_func += supg * supg_terms.supg_tau * v[w] * bf[eqn]->grad_phi[i][w];
                }
              }

              mass = 0.;

              if (pd->TimeIntegration != STEADY) {
                if (pd->e[pg->imtrx][eqn] & T_MASS) {
                  mass = b_dot[ii][jj];
                  mass *= wt_func * at * lambda * det_J * wt;
                  mass *= h3;
                  mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
                }
              }

              advection = 0.;
              if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
                if (DOUBLE_NONZERO(lambda)) {

                  advection += v_dot_del_b[ii][jj] - x_dot_del_b[ii][jj];
                  advection -= b_dot_g[ii][jj];
                  advection -= a_dot_b[ii][jj];
                  advection *= wt_func * at * lambda * det_J * wt * h3;
                  advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
                }
              }

              diffusion = 0.;
              if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {
                if (vn->shockcaptureModel == SC_YZBETA) {
                  dbl Z = b_dot[ii][jj];
                  Z += 1e-16 + v_dot_del_b[ii][jj] - x_dot_del_b[ii][jj];
                  Z -= b_dot_g[ii][jj];
                  Z -= a_dot_b[ii][jj];
                  Z *= at * lambda;
                  Z += source_term[ii][jj];

                  dbl Y_inv = 1.0;
                  dbl hdc = 0;
                  dbl js = 0;
                  for (int k = 0; k < ei[pg->imtrx]->dof[eqn]; k++) {
                    for (int p = 0; p < VIM; p++) {
                      js += fabs(grad_b[p][ii][jj] * bf[eqn]->grad_phi[k][p]);
                    }
                  }
                  hdc = 1 / (js + 1e-16);

                  dbl inner = 0;
                  for (int p = 0; p < VIM; p++) {
                    inner += Y_inv * grad_b[p][ii][jj] * grad_b[p][ii][jj];
                  }

                  dbl kdc = 0;
                  for (int ib = 0; ib < 2; ib++) {
                    dbl bt = beta[ib];
                    kdc += fabs(Y_inv * Z) * pow(hdc, bt) * pow(fabs(b[ii][jj]), 1 - bt) *
                           pow(inner, bt / 2 - 1);
                    // kdc += pow(inner, bt / 2 - 1);
                  }
                  kdc *= 0.5;
                  for (int r = 0; r < VIM; r++) {
                    diffusion += kdc * grad_b[r][ii][jj] * bf[eqn]->grad_phi[i][r];
                  }

                  diffusion *= yzbeta_factor * det_J * wt * h3;
                  diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
                }
              }

              /*
               * Source term...
               */

              source = 0.;
              if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
                // consider whether saramitoCoeff should multiply here
                source = source_term[ii][jj];
                source *= wt_func * det_J * h3 * wt;

                source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
              }

              /*
               * Add contributions to this residual (globally into Resid, and
               * locally into an accumulator)
               */

              lec->R[LEC_R_INDEX(upd->ep[pg->imtrx][eqn], i)] +=
                  mass + advection + diffusion + source;
            }
          }
        }
      }
    }

    /*
     * Jacobian terms...
     */

    if (af->Assemble_Jacobian) {
      dbl R_source, R_advection; /* Places to put the raw residual portions
                                    instead of constantly recalcing them */
      for (int ii = 0; ii < VIM; ii++) {
        for (int jj = 0; jj < VIM; jj++) {
          if (ii <= jj) /* since the stress tensor is symmetric, only assemble the upper half */
          {
            eqn = R_s[mode][ii][jj];
            peqn = upd->ep[pg->imtrx][eqn];

            R_advection = v_dot_del_b[ii][jj] - x_dot_del_b[ii][jj];
            R_advection -= b_dot_g[ii][jj];
            R_advection -= a_dot_b[ii][jj];

            R_source = source_term[ii][jj];

            for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {

              wt_func = bf[eqn]->phi[i];
              /* add Petrov-Galerkin terms as necessary */
              if (supg != 0.) {
                for (w = 0; w < dim; w++) {
                  wt_func += supg * supg_terms.supg_tau * v[w] * bf[eqn]->grad_phi[i][w];
                }
              }

              /*
               * Set up some preliminaries that are needed for the (a,i)
               * equation for bunches of (b,j) column variables...
               */

              /*
               * J_S_T
               */

              var = TEMPERATURE;
              if (pd->v[pg->imtrx][var]) {
                pvar = upd->vp[pg->imtrx][var];
                for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                  phi_j = bf[var]->phi[j];

                  mass = 0.;

                  if (pd->TimeIntegration != STEADY) {
                    if (pd->e[pg->imtrx][eqn] & T_MASS) {
                      mass = b_dot[ii][jj];
                      mass *= wt_func * d_at_dT[j] * lambda * det_J * wt;
                      mass *= h3;
                      mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
                    }
                  }

                  advection = 0.;
                  if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
                    if (DOUBLE_NONZERO(lambda)) {

                      advection += R_advection;

                      advection *= wt_func * d_at_dT[j] * lambda * det_J * wt * h3;
                      advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
                    }
                  }

                  source = 0.;
                  if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
                    source = 0;
                    source *= wt_func * det_J * wt * h3;
                    source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
                  }

                  lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += mass + advection + source;
                }
              }

              /*
               * J_S_v
               */
              for (p = 0; p < WIM; p++) {
                var = VELOCITY1 + p;
                if (pd->v[pg->imtrx][var]) {
                  pvar = upd->vp[pg->imtrx][var];
                  for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                    phi_j = bf[var]->phi[j];

                    mass = 0.;

                    if (pd->TimeIntegration != STEADY) {
                      if (pd->e[pg->imtrx][eqn] & T_MASS) {
                        dbl mass_a = 0;
                        if (supg != 0.) {
                          mass_a = supg * supg_terms.supg_tau * phi_j * bf[eqn]->grad_phi[i][p];

                          for (w = 0; w < dim; w++) {
                            mass_a += supg * supg_terms.d_supg_tau_dv[p][j] * v[w] *
                                      bf[eqn]->grad_phi[i][w];
                          }

                          mass_a *= lambda * b_dot[ii][jj];
                        }
                        dbl mass_b = d_lam->v[p][j] * b_dot[ii][jj];
                        mass = mass_a + mass_b;
                        mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)] * at * det_J * wt * h3;
                      }
                    }

                    advection = 0.;

                    if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
                      if (DOUBLE_NONZERO(lambda)) {
                        advection_a = phi_j * (grad_b[p][ii][jj]);

                        advection_a *= wt_func;

                        advection_b = 0.;
                        /* Petrov-Galerkin term */
                        if (supg != 0.) {

                          advection_b =
                              supg * supg_terms.supg_tau * phi_j * bf[eqn]->grad_phi[i][p];
                          for (w = 0; w < dim; w++) {
                            advection_b += supg * supg_terms.d_supg_tau_dv[p][j] * v[w] *
                                           bf[eqn]->grad_phi[i][w];
                          }

                          advection_b *= R_advection;
                        }

                        advection_c = v_dot_del_b[ii][jj] - x_dot_del_b[ii][jj];
                        advection_c -= b_dot_g[ii][jj];
                        advection_c -= a_dot_b[ii][jj];
                        advection_c *= wt_func;
                        advection =
                            lambda * (advection_a + advection_b) + d_lam->v[p][j] * advection_c;
                        advection *= at * det_J * wt * h3;
                        advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
                      }
                    }

                    diffusion = 0.;
                    if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {
                      if (vn->shockcaptureModel == SC_YZBETA) {
                        dbl Z = b_dot[ii][jj];
                        Z += 1e-16 + v_dot_del_b[ii][jj] - x_dot_del_b[ii][jj];
                        Z -= b_dot_g[ii][jj];
                        Z -= a_dot_b[ii][jj];
                        Z *= at * lambda;
                        Z += source_term[ii][jj];

                        dbl dZ = phi_j * (grad_b[p][ii][jj]);
                        dZ *= at * lambda;
                        dbl Y_inv = 1.0;
                        dbl hdc = 0;
                        dbl js = 0;
                        for (int k = 0; k < ei[pg->imtrx]->dof[eqn]; k++) {
                          for (int r = 0; r < VIM; r++) {
                            js += fabs(grad_b[r][ii][jj] * bf[eqn]->grad_phi[k][r]);
                          }
                        }
                        hdc = 1 / (js + 1e-16);

                        dbl inner = 0;
                        for (int p = 0; p < VIM; p++) {
                          inner += Y_inv * grad_b[p][ii][jj] * grad_b[p][ii][jj];
                        }

                        dbl dkdc = 0;
                        for (int ib = 0; ib < 2; ib++) {
                          dbl bt = beta[ib];
                          dkdc += Y_inv * dZ * Y_inv * Z / fabs(Y_inv * Z) * pow(hdc, bt) *
                                  pow(fabs(b[ii][jj]), 1 - bt) * pow(inner, bt / 2 - 1);
                        }
                        dkdc *= 0.5;

                        for (int r = 0; r < VIM; r++) {
                          diffusion += dkdc * grad_b[r][ii][jj] * bf[eqn]->grad_phi[i][r];
                        }

                        diffusion *= yzbeta_factor * det_J * wt * h3;
                        diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
                      }
                    }

                    source = 0.;

                    if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
                      source_c = 0;

                      source_a = 0.;

                      source_b = 0.;
                      if (supg != 0.) {
                        source_b = supg * supg_terms.supg_tau * phi_j * bf[eqn]->grad_phi[i][p];

                        for (w = 0; w < dim; w++) {
                          source_b += supg * supg_terms.d_supg_tau_dv[p][j] * v[w] *
                                      bf[eqn]->grad_phi[i][w];
                        }

                        source_b *= R_source;
                      }

                      source = source_a + source_b + source_c;
                      source *= det_J * wt * h3;
                      source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
                    }

                    lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += mass + advection + diffusion + source;
                  }
                }
              }

              /*
               * J_S_c
               */
              var = MASS_FRACTION;
              if (pd->v[pg->imtrx][var]) {
                for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                  phi_j = bf[var]->phi[j];

                  for (w = 0; w < pd->Num_Species_Eqn; w++) {

                    source = 0.;

                    if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
                      source = 0;
                      source *= wt_func * det_J * wt * h3;
                      source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
                    }

                    if (w > 1) {
                      GOMA_EH(GOMA_ERROR, "Need more arrays for each species.");
                    }

                    lec->J[LEC_J_INDEX(peqn, MAX_PROB_VAR + w, i, j)] += source;
                  }
                }
              }

              /*
               * J_S_P
               */
              var = PRESSURE;
              if (pd->v[pg->imtrx][var]) {
                pvar = upd->vp[pg->imtrx][var];
                for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                  phi_j = bf[var]->phi[j];

                  source = 0.;
                  if (pd->e[pg->imtrx][eqn] & T_SOURCE) {

                    source *= wt_func * det_J * wt * h3;
                    source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
                  }

                  lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += source;
                }
              }

              /*
               * J_S_d
               */
              for (p = 0; p < dim; p++) {
                var = MESH_DISPLACEMENT1 + p;
                if (pd->v[pg->imtrx][var]) {
                  pvar = upd->vp[pg->imtrx][var];
                  for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                    phi_j = bf[var]->phi[j];
                    d_det_J_dmesh_pj = bf[eqn]->d_det_J_dm[p][j];
                    dh3dmesh_pj = fv->dh3dq[p] * bf[var]->phi[j];

                    mass = 0.;
                    mass_a = 0.;
                    mass_b = 0.;
                    if (pd->TimeIntegration != STEADY) {
                      if (pd->e[pg->imtrx][eqn] & T_MASS) {
                        mass_a = b_dot[ii][jj];
                        mass_a *= wt_func * (d_det_J_dmesh_pj * h3 + det_J * dh3dmesh_pj);

                        if (supg != 0.) {
                          for (w = 0; w < dim; w++) {
                            mass_b += supg * (supg_terms.supg_tau * v[w] *
                                                  bf[eqn]->d_grad_phi_dmesh[i][w][p][j] +
                                              supg_terms.d_supg_tau_dX[p][j] * v[w] *
                                                  bf[eqn]->grad_phi[i][w]);
                          }
                          mass_b *= b_dot[ii][jj] * h3 * det_J;
                        }

                        mass = mass_a + mass_b;
                        mass *= at * lambda * wt * pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
                      }
                    }

                    advection = 0.;

                    if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
                      if (DOUBLE_NONZERO(lambda)) {
                        /*
                         * Four parts:
                         *    advection_a =
                         *    	Int ( ea.(v-xdot).d(Vv)/dmesh h3 |Jv| )
                         *
                         *    advection_b =
                         *  (i)	Int ( ea.(v-xdot).Vv h3 d(|Jv|)/dmesh )
                         *  (ii)  Int ( ea.(v-xdot).d(Vv)/dmesh h3 |Jv| )
                         *  (iii) Int ( ea.(v-xdot).Vv dh3/dmesh |Jv|   )
                         *
                         * For unsteady problems, we have an
                         * additional term
                         *
                         *    advection_c =
                         *    	Int ( ea.d(v-xdot)/dmesh.Vv h3 |Jv| )
                         */

                        advection_a = R_advection;

                        advection_a *= wt_func * (d_det_J_dmesh_pj * h3 + det_J * dh3dmesh_pj);

                        d_vdotdels_dm = 0.;
                        for (q = 0; q < WIM; q++) {
                          d_vdotdels_dm += (v[q] - x_dot[q]) * d_grad_b_dmesh[q][ii][jj][p][j];
                        }

                        advection_b = d_vdotdels_dm;
                        advection_b *= wt_func * det_J * h3;

                        advection_c = 0.;
                        if (pd->TimeIntegration != STEADY) {
                          if (pd->e[pg->imtrx][eqn] & T_MASS) {
                            d_xdotdels_dm = (1. + 2. * tt) * phi_j / dt * grad_b[p][ii][jj];

                            advection_c -= d_xdotdels_dm;

                            advection_c *= wt_func * h3 * det_J;
                          }
                        }

                        advection_d = 0.;
                        if (supg != 0.) {
                          for (w = 0; w < dim; w++) {
                            advection_d += supg * (supg_terms.supg_tau * v[w] *
                                                       bf[eqn]->d_grad_phi_dmesh[i][w][p][j] +
                                                   supg_terms.d_supg_tau_dX[p][j] * v[w] *
                                                       bf[eqn]->grad_phi[i][w]);
                          }

                          advection_d *= (R_advection)*det_J * h3;
                        }

                        advection = advection_a + advection_b + advection_c + advection_d;

                        advection *= wt * at * lambda * pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
                      }
                    }

                    diffusion = 0.;
                    if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {
                      if (vn->shockcaptureModel == SC_YZBETA) {
                        dbl Z = b_dot[ii][jj];
                        Z += 1e-16 + v_dot_del_b[ii][jj] - x_dot_del_b[ii][jj];
                        Z -= b_dot_g[ii][jj];
                        Z -= a_dot_b[ii][jj];
                        Z *= at * lambda;
                        Z += source_term[ii][jj];
                        dbl dZ = 0;
                        d_vdotdels_dm = 0.;
                        for (q = 0; q < WIM; q++) {
                          d_vdotdels_dm += (v[q] - x_dot[q]) * d_grad_b_dmesh[q][ii][jj][p][j];
                        }

                        dZ = d_vdotdels_dm;

                        if (pd->TimeIntegration != STEADY) {
                          if (pd->e[pg->imtrx][eqn] & T_MASS) {
                            d_xdotdels_dm = (1. + 2. * tt) * phi_j / dt * grad_b[p][ii][jj];

                            dZ -= d_xdotdels_dm;
                          }
                        }
                        dZ *= at * lambda;

                        dbl Y_inv = 1.0;
                        dbl hdc = 0;
                        dbl djs = 0;
                        for (int k = 0; k < ei[pg->imtrx]->dof[eqn]; k++) {
                          for (int r = 0; r < VIM; r++) {
                            djs += delta(p, ii) * delta(q, jj) * grad_b[r][ii][jj] *
                                   bf[eqn]->grad_phi[k][r] * bf[var]->grad_phi[j][r] *
                                   bf[eqn]->grad_phi[k][r] /
                                   fabs(grad_b[r][ii][jj] * bf[eqn]->grad_phi[k][r] + 1e-16);
                          }
                        }
                        dbl js = 0;
                        for (int k = 0; k < ei[pg->imtrx]->dof[eqn]; k++) {
                          for (int r = 0; r < VIM; r++) {
                            js += fabs(grad_b[r][ii][jj] * bf[eqn]->grad_phi[k][r]);
                          }
                        }
                        dbl dhdc = -1 * djs / ((js + 1e-16) * (js + 1e-16));
                        hdc = 1 / (js + 1e-16);

                        dbl inner = 0;
                        for (int r = 0; r < VIM; r++) {
                          inner += Y_inv * grad_b[r][ii][jj] * grad_b[r][ii][jj];
                        }
                        dbl d_inner = 0;
                        for (int r = 0; r < VIM; r++) {
                          d_inner +=
                              Y_inv * 2.0 * grad_b[r][ii][jj] * d_grad_b_dmesh[r][ii][jj][p][j];
                        }

                        dbl kdc = 0;
                        dbl dkdc = 0;
                        for (int ib = 0; ib < 2; ib++) {
                          dbl bt = beta[ib];
                          kdc += fabs(Y_inv * Z) * pow(hdc, bt) * pow(fabs(b[ii][jj]), 1 - bt) *
                                 pow(inner, bt / 2 - 1);
                          dkdc += ((Y_inv * dZ * Y_inv * Z) / fabs(Y_inv * Z)) * pow(hdc, bt) *
                                  pow(fabs(b[ii][jj]), 1 - bt) * pow(inner, bt / 2 - 1);
                          dkdc += fabs(Y_inv * Z) * bt * dhdc * pow(hdc, bt - 1) *
                                  pow(fabs(b[ii][jj]), 1 - bt) * pow(inner, bt / 2 - 1);
                          if (DOUBLE_NONZERO(b[ii][jj])) {
                            dkdc += fabs(Y_inv * Z) * pow(hdc, bt) * (1 - bt) *
                                    (b[ii][jj] / fabs(b[ii][jj])) * bf[var]->phi[j] * delta(p, ii) *
                                    delta(q, jj) * pow(fabs(b[ii][jj]), -bt) *
                                    pow(inner, bt / 2 - 1);
                          }
                          dkdc += fabs(Y_inv * Z) * pow(hdc, bt) * pow(fabs(b[ii][jj]), 1 - bt) *
                                  d_inner * (bt / 2 - 1) * pow(inner, bt / 2 - 2);
                        }
                        kdc *= 0.5;
                        dkdc *= 0.5;

                        dbl diffusion_a = 0;
                        dbl diffusion_b = 0;
                        for (int r = 0; r < VIM; r++) {
                          diffusion_a += dkdc * grad_b[r][ii][jj] * bf[eqn]->grad_phi[i][r];
                          diffusion_a +=
                              kdc * d_grad_b_dmesh[r][ii][jj][p][j] * bf[eqn]->grad_phi[i][r];
                          diffusion_b += kdc * grad_b[r][ii][jj] * bf[eqn]->grad_phi[i][r];
                        }
                        diffusion_a *= yzbeta_factor * det_J * wt * h3;
                        diffusion_b *=
                            yzbeta_factor * (d_det_J_dmesh_pj * h3 + det_J * dh3dmesh_pj) * wt;
                        diffusion = diffusion_a + diffusion_b;
                        diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
                      }
                    }

                    /*
                     * Source term...
                     */

                    source = 0.;

                    if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
                      source_a = R_source;

                      source_a *= wt_func * (d_det_J_dmesh_pj * h3 + det_J * dh3dmesh_pj);

                      source_c = 0.;
                      if (supg != 0.) {
                        for (w = 0; w < dim; w++) {
                          source_c +=
                              supg *
                              (supg_terms.supg_tau * v[w] * bf[eqn]->d_grad_phi_dmesh[i][w][p][j] +
                               supg_terms.d_supg_tau_dX[p][j] * v[w] * bf[eqn]->grad_phi[i][w]);
                        }
                        source_c *= R_source * det_J * h3;
                      }

                      source = source_a + source_c;

                      source *= wt * pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
                    }

                    lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += mass + advection + diffusion + source;
                  }
                }
              }

              /*
               * J_S_G
               */
              if (evss_gradv == 0) {
                for (p = 0; p < VIM; p++) {
                  for (q = 0; q < VIM; q++) {
                    var = v_g[p][q];

                    if (pd->v[pg->imtrx][var]) {
                      pvar = upd->vp[pg->imtrx][var];
                      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                        phi_j = bf[var]->phi[j];
                        advection = 0.;
                        advection_a = 0.;
                        if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
                          if (DOUBLE_NONZERO(lambda)) {
                            for (int k = 0; k < VIM; k++) {
                              advection += -b[ii][k] * delta(p, k) * delta(jj, q);
                            }
                            advection += -d_a_dot_b_dG[p][q][ii][jj];
                            advection *= phi_j * h3 * det_J;

                            advection *= wt_func * wt * at * lambda *
                                         pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
                          }
                        }

                        /*
                         * Diffusion...
                         */

                        diffusion = 0.;

                        if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {
                          if (vn->shockcaptureModel == SC_YZBETA) {
                            dbl Z = b_dot[ii][jj];
                            Z += 1e-16 + v_dot_del_b[ii][jj] - x_dot_del_b[ii][jj];
                            Z -= b_dot_g[ii][jj];
                            Z -= a_dot_b[ii][jj];
                            Z *= at * lambda;
                            Z += source_term[ii][jj];
                            dbl dZ = 0;
                            for (int k = 0; k < VIM; k++) {
                              dZ += -b[ii][k] * delta(p, k) * delta(jj, q);
                            }
                            dZ += -d_a_dot_b_dG[p][q][ii][jj];
                            dZ *= bf[var]->phi[j];
                            dZ *= at * lambda;
                            dbl Y_inv = 1.0;
                            dbl hdc = 0;
                            dbl js = 0;
                            for (int k = 0; k < ei[pg->imtrx]->dof[eqn]; k++) {
                              for (int r = 0; r < VIM; r++) {
                                js += fabs(grad_b[r][ii][jj] * bf[eqn]->grad_phi[k][r]);
                              }
                            }
                            hdc = 1 / (js + 1e-16);

                            dbl inner = 0;
                            for (int p = 0; p < VIM; p++) {
                              inner += Y_inv * grad_b[p][ii][jj] * grad_b[p][ii][jj];
                            }

                            dbl dkdc = 0;
                            for (int ib = 0; ib < 2; ib++) {
                              dbl bt = beta[ib];
                              dkdc += Y_inv * dZ * Y_inv * Z / fabs(Y_inv * Z) * pow(hdc, bt) *
                                      pow(fabs(b[ii][jj]), 1 - bt) * pow(inner, bt / 2 - 1);
                            }
                            dkdc *= 0.5;
                            // dkdc = 0;
                            for (int r = 0; r < VIM; r++) {
                              diffusion += dkdc * grad_b[r][ii][jj] * bf[eqn]->grad_phi[i][r];
                            }

                            diffusion *= yzbeta_factor * det_J * wt * h3;
                            diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
                          }
                        }

                        /*
                         * Source term...
                         */

                        source = 0.;

                        lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += advection + diffusion + source;
                      }
                    }
                  }
                }
              }

              /*
               * J_S_F
               */
              var = FILL;
              if (pd->v[pg->imtrx][var]) {
                pvar = upd->vp[pg->imtrx][var];
                for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                  phi_j = bf[var]->phi[j];

                  mass = 0.;

                  if (pd->TimeIntegration != STEADY) {
                    if (pd->e[pg->imtrx][eqn] & T_MASS) {

                      mass = b_dot[ii][jj];
                      mass *= d_lam->F[j];
                      mass *= wt_func * at * det_J * wt;
                      mass *= h3;
                      mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
                    }
                  }

                  advection = 0.;

                  if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
                    if (d_lam->F[j] != 0.) {

                      advection += R_advection;
                      advection *= d_lam->F[j];
                      advection *= wt_func * at * det_J * wt * h3;
                      advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
                    }
                  }

                  diffusion = 0.;

                  if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {
                    /* add SU term in here when appropriate */

                    diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
                  }

                  source = 0.;

                  if (pd->e[pg->imtrx][eqn] & T_SOURCE) {

                    source *= wt_func * det_J * h3 * wt;

                    source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
                  }

                  lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += mass + advection + diffusion + source;
                }
              }

              /*
               * J_S_S
               */
              for (p = 0; p < VIM; p++) {
                for (q = 0; q < VIM; q++) {
                  if (q >= p) {
                    var = v_s[mode][p][q];

                    if (pd->v[pg->imtrx][var]) {
                      pvar = upd->vp[pg->imtrx][var];
                      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                        phi_j = bf[var]->phi[j];
                        mass = 0.;
                        if (pd->TimeIntegration != STEADY) {
                          if (pd->e[pg->imtrx][eqn] & T_MASS) {
                            mass = (1. + 2. * tt) * phi_j / dt * (double)delta(ii, p) *
                                   (double)delta(jj, q);
                            mass *= h3 * det_J;
                            mass *=
                                wt_func * at * lambda * wt * pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
                          }
                        }

                        advection = 0.;

                        if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
                          if (DOUBLE_NONZERO(lambda)) {
                            if ((ii == p) && (jj == q)) {
                              for (r = 0; r < WIM; r++) {
                                advection += (v[r] - x_dot[r]) * bf[var]->grad_phi[j][r];
                              }
                            }

                            for (int k = 0; k < VIM; k++) {
                              advection -=
                                  phi_j *
                                  (delta(ii, q) * delta(k, p) | delta(ii, p) * delta(k, q)) *
                                  g[k][jj];
                            }
                            advection -= phi_j * d_a_dot_b_db[p][q][ii][jj];

                            advection *= h3 * det_J;

                            advection *= wt_func * wt * at * lambda *
                                         pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
                          }
                        }

                        /*
                         * Diffusion...
                         */

                        diffusion = 0.;

                        if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {
                          if (vn->shockcaptureModel == SC_YZBETA) {
                            dbl Z = b_dot[ii][jj];
                            Z += 1e-16 + v_dot_del_b[ii][jj] - x_dot_del_b[ii][jj];
                            Z -= b_dot_g[ii][jj];
                            Z -= a_dot_b[ii][jj];
                            Z *= at * lambda;
                            Z += source_term[ii][jj];
                            dbl dZ = 0;
                            if (pd->TimeIntegration != STEADY) {
                              if (pd->e[pg->imtrx][eqn] & T_MASS) {
                                dZ = (1. + 2. * tt) * phi_j / dt * (double)delta(ii, p) *
                                     (double)delta(jj, q);
                              }
                            }
                            if ((ii == p) && (jj == q)) {
                              for (r = 0; r < WIM; r++) {
                                dZ += (v[r] - x_dot[r]) * bf[var]->grad_phi[j][r];
                              }
                            }
                            for (int k = 0; k < VIM; k++) {
                              dZ -= phi_j *
                                    (delta(ii, q) * delta(k, p) | delta(ii, p) * delta(k, q)) *
                                    g[k][jj];
                            }
                            dZ -= phi_j * d_a_dot_b_db[p][q][ii][jj];
                            dZ *= at * lambda;
                            dZ += d_source_term_db[ii][jj][p][q] * bf[var]->phi[j];

                            dbl Y_inv = 1.0;
                            dbl hdc = 0;
                            dbl djs = 0;
                            for (int k = 0; k < ei[pg->imtrx]->dof[eqn]; k++) {
                              for (int r = 0; r < VIM; r++) {
                                djs += delta(p, ii) * delta(q, jj) * grad_b[r][ii][jj] *
                                       bf[eqn]->grad_phi[k][r] * bf[var]->grad_phi[j][r] *
                                       bf[eqn]->grad_phi[k][r] /
                                       fabs(grad_b[r][ii][jj] * bf[eqn]->grad_phi[k][r] + 1e-16);
                              }
                            }
                            dbl js = 0;
                            for (int k = 0; k < ei[pg->imtrx]->dof[eqn]; k++) {
                              for (int r = 0; r < VIM; r++) {
                                js += fabs(grad_b[r][ii][jj] * bf[eqn]->grad_phi[k][r]);
                              }
                            }
                            dbl dhdc = -1 * djs / ((js + 1e-16) * (js + 1e-16));
                            hdc = 1 / (js + 1e-16);
                            dbl inner = 0;
                            for (int r = 0; r < VIM; r++) {
                              inner += Y_inv * grad_b[r][ii][jj] * grad_b[r][ii][jj];
                            }
                            dbl d_inner = 0;
                            for (int r = 0; r < VIM; r++) {
                              d_inner += Y_inv * 2.0 * delta(p, ii) * delta(q, jj) *
                                         grad_b[r][ii][jj] * bf[var]->grad_phi[j][r];
                            }

                            dbl kdc = 0;
                            dbl dkdc = 0;
                            for (int ib = 0; ib < 2; ib++) {
                              dbl bt = beta[ib];
                              kdc += fabs(Y_inv * Z) * pow(hdc, bt) * pow(fabs(b[ii][jj]), 1 - bt) *
                                     pow(inner, bt / 2 - 1);
                              dkdc += ((Y_inv * dZ * Y_inv * Z) / fabs(Y_inv * Z)) * pow(hdc, bt) *
                                      pow(fabs(b[ii][jj]), 1 - bt) * pow(inner, bt / 2 - 1);
                              dkdc += fabs(Y_inv * Z) * bt * dhdc * pow(hdc, bt - 1) *
                                      pow(fabs(b[ii][jj]), 1 - bt) * pow(inner, bt / 2 - 1);
                              if (DOUBLE_NONZERO(b[ii][jj])) {
                                dkdc += fabs(Y_inv * Z) * pow(hdc, bt) * (1 - bt) *
                                        (b[ii][jj] / fabs(b[ii][jj])) * bf[var]->phi[j] *
                                        delta(p, ii) * delta(q, jj) * pow(fabs(b[ii][jj]), -bt) *
                                        pow(inner, bt / 2 - 1);
                              }
                              dkdc += fabs(Y_inv * Z) * pow(hdc, bt) *
                                      pow(fabs(b[ii][jj]), 1 - bt) * d_inner * (bt / 2 - 1) *
                                      pow(inner, bt / 2 - 2);
                            }
                            kdc *= 0.5;
                            dkdc *= 0.5;

                            for (int r = 0; r < VIM; r++) {
                              diffusion += dkdc * grad_b[r][ii][jj] * bf[eqn]->grad_phi[i][r];
                              diffusion += kdc * delta(p, ii) * delta(q, jj) *
                                           bf[var]->grad_phi[j][r] * bf[eqn]->grad_phi[i][r];
                            }
                          }

                          diffusion *= yzbeta_factor * det_J * wt * h3;
                          diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
                        }

                        /*
                         * Source term...
                         */

                        source = 0.;

                        if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
                          source = d_source_term_db[ii][jj][p][q];
                          source *= phi_j * det_J * h3 * wt_func * wt *
                                    pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
                        }

                        lec->J[LEC_J_INDEX(peqn, pvar, i, j)] +=
                            mass + advection + diffusion + source;
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    } /* End Assemble Jacobian */
  }   /* End loop over modes */

  return (status);
}

int conf_source(int mode,
                dbl c[DIM][DIM],
                dbl source_term[DIM][DIM],
                dbl d_source_term_dc[DIM][DIM][DIM][DIM]) {
  bool compute_derivatives = d_source_term_dc != NULL && af->Assemble_Jacobian;
  switch (vn->ConstitutiveEquation) {
  case OLDROYDB: {
    for (int ii = 0; ii < VIM; ii++) {
      for (int jj = 0; jj < VIM; jj++) {
        source_term[ii][jj] = -(c[ii][jj] - delta(ii, jj));
      }
    }
    if (compute_derivatives) {
      for (int ii = 0; ii < VIM; ii++) {
        for (int jj = 0; jj < VIM; jj++) {
          for (int p = 0; p < VIM; p++) {
            for (int q = 0; q < VIM; q++) {
              d_source_term_dc[ii][jj][p][q] = -delta(ii, p) * delta(jj, q);
            }
          }
        }
      }
    }
  } break;
  case PTT: {

    dbl d_trace_dc[DIM][DIM] = {{0.0}};

    dbl trace = 0;
    for (int i = 0; i < VIM; i++) {
      trace += c[i][i];
    }

    if (compute_derivatives) {
      for (int p = 0; p < VIM; p++) {
        for (int q = 0; q < VIM; q++) {
          d_trace_dc[p][q] = 0.0;
          for (int i = 0; i < VIM; i++) {
            for (int j = 0; j < VIM; j++) {
              d_trace_dc[p][q] += 2.0 * c[i][i] * (delta(p, i) * delta(q, j));
            }
          }
        }
      }
    }

    dbl Z = 1.0;
    dbl dZ_dtrace = 0;

    // PTT exponent
    eps = ve[mode]->eps;

    if (vn->ptt_type == PTT_LINEAR) {
      Z = 1 + eps * (trace - (double)VIM);
      dZ_dtrace = eps;
    } else if (vn->ptt_type == PTT_EXPONENTIAL) {
      Z = exp(eps * (trace - (double)VIM));
      dZ_dtrace = eps * Z;
    } else {
      GOMA_EH(GOMA_ERROR, "Unrecognized PTT Form %d", vn->ptt_type);
    }

    for (int ii = 0; ii < VIM; ii++) {
      for (int jj = 0; jj < VIM; jj++) {
        source_term[ii][jj] = -Z * (c[ii][jj] - delta(ii, jj));
      }
    }
    if (compute_derivatives) {
      for (int ii = 0; ii < VIM; ii++) {
        for (int jj = 0; jj < VIM; jj++) {
          for (int p = 0; p < VIM; p++) {
            for (int q = 0; q < VIM; q++) {
              d_source_term_dc[ii][jj][p][q] =
                  -Z * (delta(ii, p) * delta(jj, q)) -
                  dZ_dtrace * d_trace_dc[p][q] * (c[ii][jj] - delta(ii, jj));
            }
          }
        }
      }
    }
  } break;
  case ROLIE_POLY:
  case ROLIE_POLY_FE: {

    dbl d_trace_dc[DIM][DIM] = {{0.0}};

    dbl trace = 1e-16;
    for (int i = 0; i < VIM; i++) {
      trace += c[i][i];
    }

    dbl tau_D = ve[mode]->time_const_st->lambda0;
    dbl tau_R = ve[mode]->stretch_time;
    dbl beta = ve[mode]->CCR_coefficient;
    dbl n = ve[mode]->polymer_exponent;
    dbl lambda_max = ve[mode]->maximum_stretch_ratio;

    dbl lambda_s = sqrt(trace / 3);

    dbl k = 1.0;
    if (vn->ConstitutiveEquation == ROLIE_POLY_FE) {
      k = ((3 - lambda_s * lambda_s / (lambda_max * lambda_max)) *
           (1 - 1 / (lambda_max * lambda_max))) /
          (1 -
           lambda_s * lambda_s / (lambda_max * lambda_max) * (3 - 1 / (lambda_max * lambda_max)));
    }

    for (int ii = 0; ii < VIM; ii++) {
      for (int jj = 0; jj < VIM; jj++) {
        source_term[ii][jj] = -(c[ii][jj] - delta(ii, jj));
        source_term[ii][jj] += -(2.0 * tau_D / tau_R) * k * (1 - sqrt(3 / (trace))) *
                               (c[ii][jj] + beta * pow(trace / 3, n) * (c[ii][jj] - delta(ii, jj)));
      }
    }

    if (compute_derivatives) {
      for (int p = 0; p < VIM; p++) {
        for (int q = 0; q < VIM; q++) {
          d_trace_dc[p][q] = 0.0;
          for (int i = 0; i < VIM; i++) {
            d_trace_dc[p][q] += delta(i, p) * delta(i, q);
          }
        }
      }
      dbl d_k[DIM][DIM] = {{0.}};
      if (vn->ConstitutiveEquation == ROLIE_POLY_FE) {
        dbl d_k_dlamda_s = 4 * lambda_s * lambda_max * lambda_max * (lambda_max * lambda_max - 1) /
                           (pow(lambda_s * lambda_s - lambda_max * lambda_max, 2.0) *
                            (3 * lambda_max * lambda_max - 1));

        for (int p = 0; p < VIM; p++) {
          for (int q = 0; q < VIM; q++) {
            d_k[p][q] = d_k_dlamda_s * d_trace_dc[p][q];
          }
        }
      }

      for (int ii = 0; ii < VIM; ii++) {
        for (int jj = 0; jj < VIM; jj++) {
          for (int p = 0; p < VIM; p++) {
            for (int q = 0; q < VIM; q++) {
              d_source_term_dc[ii][jj][p][q] = -(delta(ii, p) * delta(jj, q));
              d_source_term_dc[ii][jj][p][q] +=
                  -(2.0 * tau_D / tau_R) * d_k[p][q] * (1 - sqrt(3 / (trace))) *
                      (c[ii][jj] + beta * pow(trace / 3, n) * (c[ii][jj] - delta(ii, jj))) -
                  (2.0 * tau_D / tau_R) * k * d_trace_dc[p][q] * sqrt(3) / 2.0 *
                      pow(1 / (trace), 3.0 / 2.0) *
                      (c[ii][jj] + beta * pow(trace / 3, n) * (c[ii][jj] - delta(ii, jj))) -
                  (2.0 * tau_D / tau_R) * k * (1 - sqrt(3 / (trace))) *
                      ((delta(ii, p) * delta(jj, q)) +
                       d_trace_dc[p][q] * (n / 3) * beta * pow(trace / 3, n - 1) *
                           (c[ii][jj] - delta(ii, jj)) +
                       beta * pow(trace / 3, n) * ((delta(ii, p) * delta(jj, q))));
            }
          }
        }
      }
    }
  } break;
  default:
    GOMA_EH(GOMA_ERROR, "Unknown Constitutive equation form for SQRT_CONF");
    break;
  }

  return GOMA_SUCCESS;
}

int assemble_stress_conf(dbl tt, /* parameter to vary time integration from
                                  * explicit (tt = 1) to implicit (tt = 0) */
                         dbl dt, /* current time step size */
                         PG_DATA *pg_data) {
  int dim, p, q, r, a, b, w, k;

  int eqn, var;
  int peqn, pvar;
  int evss_gradv = 0;

  int i, j, status, mode;
  dbl v[DIM];      /* Velocity field. */
  dbl x_dot[DIM];  /* current position field derivative wrt time. */
  dbl h3;          /* Volume element (scale factors). */
  dbl dh3dmesh_pj; /* Sensitivity to (p,j) mesh dof. */

  dbl grad_v[DIM][DIM];
  dbl gamma[DIM][DIM]; /* Shear-rate tensor based on velocity */
  dbl det_J;           /* determinant of element Jacobian */

  dbl d_det_J_dmesh_pj; /* for specific (p,j) mesh dof */

  dbl mass; /* For terms and their derivatives */
  dbl mass_a, mass_b;
  dbl advection;
  dbl advection_a, advection_b, advection_c, advection_d;
  dbl diffusion;
  dbl source;
  dbl source1;
  dbl source_a = 0, source_b = 0, source_c = 0;
  int err;
  dbl alpha = 0;  /* This is the Geisekus mobility parameter */
  dbl lambda = 0; /* polymer relaxation constant */
  double xi;
  double d_xi_dF[MDE];
  dbl ucwt, lcwt; /* Upper convected derviative weight, Lower convected derivative weight */
  dbl eps = 0;    /* This is the PTT elongation parameter */
  double d_eps_dF[MDE];
  /*
   *
   * Note how carefully we avoid refering to d(phi[i])/dx[j] and refer instead
   * to the j-th component of grad_phi[j][i] so that this vector can be loaded
   * up with components that may be different in non Cartesian coordinate
   * systems.
   *
   * We will, however, insist on *orthogonal* coordinate systems, even if we
   * might permit them to be curvilinear.
   *
   * Assume all components of velocity are interpolated with the same kind
   * of basis function.
   */

  /*
   * Petrov-Galerkin weighting functions for i-th and ab-th stress residuals
   * and some of their derivatives...
   */

  dbl wt_func;

  /*
   * Interpolation functions for variables and some of their derivatives.
   */

  dbl phi_j;
  dbl wt;

  /* Variables for stress */

  int R_s[MAX_MODES][DIM][DIM];
  int v_s[MAX_MODES][DIM][DIM];
  int v_g[DIM][DIM];

  dbl s[DIM][DIM];     /* stress tensor */
  dbl s_dot[DIM][DIM]; /* stress tensor from last time step */
  dbl grad_s[DIM][DIM][DIM];
  dbl d_grad_s_dmesh[DIM][DIM][DIM][DIM]
                    [MDE]; /* derivative of grad of stress tensor for mode ve_mode */

  dbl g[DIM][DIM];  /* velocity gradient tensor */
  dbl gt[DIM][DIM]; /* transpose of velocity gradient tensor */

  /* dot product tensors */

  dbl s_dot_s[DIM][DIM];
  dbl s_dot_g[DIM][DIM];
  dbl g_dot_s[DIM][DIM];
  dbl s_dot_gt[DIM][DIM];
  dbl gt_dot_s[DIM][DIM];

  /* polymer viscosity and derivatives */
  dbl mup;
  VISCOSITY_DEPENDENCE_STRUCT d_mup_struct;
  VISCOSITY_DEPENDENCE_STRUCT *d_mup = &d_mup_struct;
  POLYMER_TIME_CONST_DEPENDENCE_STRUCT d_lam_struct;
  POLYMER_TIME_CONST_DEPENDENCE_STRUCT *d_lam = &d_lam_struct;

  SARAMITO_DEPENDENCE_STRUCT d_saramito_struct;
  SARAMITO_DEPENDENCE_STRUCT *d_saramito = &d_saramito_struct;

  // todo: will want to parse necessary parameters... for now hard code
  const bool saramitoEnabled =
      (vn->ConstitutiveEquation == SARAMITO_OLDROYDB || vn->ConstitutiveEquation == SARAMITO_PTT ||
       vn->ConstitutiveEquation == SARAMITO_GIESEKUS);

  dbl saramitoCoeff = 1.;

  dbl d_mup_dmesh_pj;

  /*  shift function */
  dbl at = 0.0;
  dbl d_at_dT[MDE];
  dbl wlf_denom;

  /* constitutive equation parameters */
  dbl Z = 1.0; /* This is the factor appearing in front of the stress tensor in PTT */

  /* advective terms are precalculated */
  dbl v_dot_del_s[DIM][DIM];
  dbl x_dot_del_s[DIM][DIM];

  dbl d_xdotdels_dm;

  dbl d_vdotdels_dm;

  dbl trace = 0.0; /* trace of the stress tensor */

  /* SUPG variables */
  dbl supg = 0;

  if (vn->evssModel == EVSS_GRADV) {
    evss_gradv = 1;
  }

  status = 0;

  eqn = R_STRESS11;

  /*
   * Bail out fast if there's nothing to do...
   */

  if (!pd->e[pg->imtrx][eqn]) {
    return (status);
  }

  /*
   * Unpack variables from structures for local convenience...
   */

  dim = pd->Num_Dim;

  wt = fv->wt;

  det_J = bf[eqn]->detJ; /* Really, ought to be mesh eqn. */

  h3 = fv->h3; /* Differential volume element (scales). */

  /* load eqn and variable number in tensor form */
  (void)stress_eqn_pointer(v_s);
  (void)stress_eqn_pointer(R_s);

  v_g[0][0] = VELOCITY_GRADIENT11;
  v_g[0][1] = VELOCITY_GRADIENT12;
  v_g[1][0] = VELOCITY_GRADIENT21;
  v_g[1][1] = VELOCITY_GRADIENT22;
  v_g[0][2] = VELOCITY_GRADIENT13;
  v_g[1][2] = VELOCITY_GRADIENT23;
  v_g[2][0] = VELOCITY_GRADIENT31;
  v_g[2][1] = VELOCITY_GRADIENT32;
  v_g[2][2] = VELOCITY_GRADIENT33;

  /*
   * Field variables...
   */
  for (a = 0; a < WIM; a++) {
    v[a] = fv->v[a];

    /* note, these are zero for steady calculations */
    x_dot[a] = 0.0;
    if (pd->TimeIntegration != STEADY && pd->gv[MESH_DISPLACEMENT1 + a]) {
      x_dot[a] = fv_dot->x[a];
    }
  }

  /*
   * In Cartesian coordinates, this velocity gradient tensor will
   * have components that are...
   *
   * 			grad_v[a][b] = d v_b
   *				       -----
   *				       d x_a
   */

  for (a = 0; a < VIM; a++) {
    for (b = 0; b < VIM; b++) {
      grad_v[a][b] = fv->grad_v[a][b];
    }
  }

  /* load up shearrate tensor based on velocity */
  for (a = 0; a < VIM; a++) {
    for (b = 0; b < VIM; b++) {
      gamma[a][b] = grad_v[a][b] + grad_v[b][a];
    }
  }

  for (a = 0; a < VIM; a++) {
    for (b = 0; b < VIM; b++) {
      if (evss_gradv) {
        g[a][b] = fv->grad_v[a][b];
        gt[a][b] = fv->grad_v[b][a];
      } else {
        g[a][b] = fv->G[a][b];
        gt[b][a] = g[a][b];
      }
    }
  }

  if (vn->wt_funcModel == GALERKIN) {
    supg = 0.;
  } else if (vn->wt_funcModel == SUPG) {
    supg = vn->wt_func;
  }

  SUPG_terms supg_terms;
  if (supg != 0.) {
    supg_tau(&supg_terms, dim, 0.0, pg_data, dt, TRUE, eqn);
  }
  /* end Petrov-Galerkin addition */

  /*  shift factor  */
  if (pd->gv[TEMPERATURE]) {
    if (vn->shiftModel == CONSTANT) {
      at = vn->shift[0];
      for (j = 0; j < ei[pg->imtrx]->dof[TEMPERATURE]; j++) {
        d_at_dT[j] = 0.;
      }
    } else if (vn->shiftModel == MODIFIED_WLF) {
      wlf_denom = vn->shift[1] + fv->T - mp->reference[TEMPERATURE];
      if (wlf_denom != 0.) {
        at = exp(vn->shift[0] * (mp->reference[TEMPERATURE] - fv->T) / wlf_denom);
        for (j = 0; j < ei[pg->imtrx]->dof[TEMPERATURE]; j++) {
          d_at_dT[j] =
              -at * vn->shift[0] * vn->shift[1] / (wlf_denom * wlf_denom) * bf[TEMPERATURE]->phi[j];
        }
      } else {
        at = 1.;
      }
      for (j = 0; j < ei[pg->imtrx]->dof[TEMPERATURE]; j++) {
        d_at_dT[j] = 0.;
      }
    }
  } else {
    at = 1.;
  }

  /* Begin loop over modes */
  for (mode = 0; mode < vn->modes; mode++) {

    load_modal_pointers(mode, tt, dt, s, s_dot, grad_s, d_grad_s_dmesh);

    /* precalculate advective terms of form (v dot del tensor)*/

    trace = 0.0;

    for (a = 0; a < VIM; a++) {
      trace += s[a][a];
      for (b = 0; b < VIM; b++) {
        v_dot_del_s[a][b] = 0.;
        x_dot_del_s[a][b] = 0.;
        for (q = 0; q < WIM; q++) {
          v_dot_del_s[a][b] += v[q] * grad_s[q][a][b];
          x_dot_del_s[a][b] += x_dot[q] * grad_s[q][a][b];
        }
      }
    }

    /*
     * Stress tensor...(Note "anti-BSL" sign convention on deviatoric stress)
     */

    /* get polymer viscosity */
    mup = viscosity(ve[mode]->gn, gamma, d_mup);

    if (saramitoEnabled == TRUE) {
      compute_saramito_model_terms(&saramitoCoeff, d_saramito, s, ve[mode]->gn, FALSE);
    } else {
      saramitoCoeff = 1.;
      d_saramito->tau_y = 0;

      for (int i = 0; i < VIM; ++i) {
        for (int j = 0; j < VIM; ++j) {
          d_saramito->s[i][j] = 0;
        }
      }
    }

    double d_alpha_dF[MDE];
    /* get Geisekus mobility parameter */
    if (ve[mode]->alphaModel == CONSTANT) {
      alpha = ve[mode]->alpha;
    } else if (ls != NULL && ve[mode]->alphaModel == VE_LEVEL_SET) {
      double pos_alpha = ve[mode]->pos_ls.alpha;
      double neg_alpha = ve[mode]->alpha;
      double width = ls->Length_Scale;
      err = level_set_property(neg_alpha, pos_alpha, width, &alpha, d_alpha_dF);
      GOMA_EH(err, "level_set_property() failed for mobility parameter.");
    } else {
      GOMA_EH(GOMA_ERROR, "Unknown mobility parameter model");
    }

    /* get time constant */
    lambda = polymer_time_const(ve[mode]->time_const_st, gamma, d_lam);

    xi = 0;
    if (ve[mode]->xiModel == CONSTANT) {
      xi = ve[mode]->xi;
    } else if (ls != NULL && ve[mode]->xiModel == VE_LEVEL_SET) {
      double pos_xi = ve[mode]->pos_ls.xi;
      double neg_xi = ve[mode]->xi;
      double width = ls->Length_Scale;
      err = level_set_property(neg_xi, pos_xi, width, &xi, d_xi_dF);
      GOMA_EH(err, "level_set_property() failed for ptt xi parameter.");
    } else {
      GOMA_EH(GOMA_ERROR, "Unknown PTT Xi parameter model");
    }

    ucwt = 1.0 - xi / 2.0;
    lcwt = xi / 2.0;

    if (ve[mode]->epsModel == CONSTANT) {
      eps = ve[mode]->eps;
    } else if (ls != NULL && ve[mode]->epsModel == VE_LEVEL_SET) {
      double pos_eps = ve[mode]->pos_ls.eps;
      double neg_eps = ve[mode]->eps;
      double width = ls->Length_Scale;
      err = level_set_property(neg_eps, pos_eps, width, &eps, d_eps_dF);
      GOMA_EH(err, "level_set_property() failed for ptt epsilon parameter.");
    } else {
      GOMA_EH(GOMA_ERROR, "Unknown PTT Epsilon parameter model");
    }

    /* get tensor dot products for future use */

    if (DOUBLE_NONZERO(alpha))
      (void)tensor_dot(s, s, s_dot_s, VIM);

    if (ucwt != 0.) {
      (void)tensor_dot(s, g, s_dot_g, VIM);
      (void)tensor_dot(gt, s, gt_dot_s, VIM);
    }

    if (lcwt != 0.) {
      (void)tensor_dot(s, gt, s_dot_gt, VIM);
      (void)tensor_dot(g, s, g_dot_s, VIM);
    }
    /*
     * Residuals_________________________________________________________________
     */

    dbl source_term[DIM][DIM];
    dbl d_source_term_dc[DIM][DIM][DIM][DIM];
    conf_source(mode, s, source_term, d_source_term_dc);
    if (af->Assemble_Residual) {
      /*
       * Assemble each component "ab" of the polymer stress equation...
       */
      for (a = 0; a < VIM; a++) {
        for (b = 0; b < VIM; b++) {

          if (a <= b) /* since the stress tensor is symmetric, only assemble the upper half */
          {
            eqn = R_s[mode][a][b];

            /*
             * In the element, there will be contributions to this many equations
             * based on the number of degrees of freedom...
             */

            for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
              wt_func = bf[eqn]->phi[i];
              /* add Petrov-Galerkin terms as necessary */
              if (supg != 0.) {
                for (w = 0; w < dim; w++) {
                  wt_func += supg * supg_terms.supg_tau * v[w] * bf[eqn]->grad_phi[i][w];
                }
              }

              mass = 0.;

              if (pd->TimeIntegration != STEADY) {
                if (pd->e[pg->imtrx][eqn] & T_MASS) {
                  mass = s_dot[a][b];
                  mass *= wt_func * at * det_J * wt;
                  mass *= h3;
                  mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
                }
              }

              advection = 0.;
              if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
                if (DOUBLE_NONZERO(lambda)) {

                  advection += v_dot_del_s[a][b] - x_dot_del_s[a][b];
                  if (ucwt != 0.)
                    advection -= ucwt * (gt_dot_s[a][b] + s_dot_g[a][b]);
                  if (lcwt != 0.)
                    advection += lcwt * (s_dot_gt[a][b] + g_dot_s[a][b]);

                  advection *= wt_func * at * det_J * wt * h3;
                  advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
                }
              }

              diffusion = 0.;
              if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {
                diffusion *= det_J * wt * h3;
                diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
              }

              /*
               * Source term...
               */

              source = 0.;
              if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
                // consider whether saramitoCoeff should multiply here
                source -= source_term[a][b];
                source *= (1 / lambda) * wt_func * det_J * h3 * wt;

                source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
              }

              /*
               * Add contributions to this residual (globally into Resid, and
               * locally into an accumulator)
               */

              lec->R[LEC_R_INDEX(upd->ep[pg->imtrx][eqn], i)] +=
                  mass + advection + diffusion + source;
            }
          }
        }
      }
    }

    /*
     * Jacobian terms...
     */

    if (af->Assemble_Jacobian) {
      dbl R_source, R_advection; /* Places to put the raw residual portions
                                    instead of constantly recalcing them */
      for (a = 0; a < VIM; a++) {
        for (b = 0; b < VIM; b++) {
          if (a <= b) /* since the stress tensor is symmetric, only assemble the upper half */
          {
            eqn = R_s[mode][a][b];
            peqn = upd->ep[pg->imtrx][eqn];

            R_advection = v_dot_del_s[a][b] - x_dot_del_s[a][b];
            if (ucwt != 0.)
              R_advection -= ucwt * (gt_dot_s[a][b] + s_dot_g[a][b]);
            if (lcwt != 0.)
              R_advection += lcwt * (s_dot_gt[a][b] + g_dot_s[a][b]);

            R_source = Z * s[a][b];

            if (DOUBLE_NONZERO(alpha))
              R_source += alpha * lambda * (s_dot_s[a][b] / mup);
            R_source *= saramitoCoeff;
            R_source += -at * mup * (g[a][b] + gt[a][b]);

            for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {

              wt_func = bf[eqn]->phi[i];
              /* add Petrov-Galerkin terms as necessary */
              if (supg != 0.) {
                for (w = 0; w < dim; w++) {
                  wt_func += supg * supg_terms.supg_tau * v[w] * bf[eqn]->grad_phi[i][w];
                }
              }

              /*
               * Set up some preliminaries that are needed for the (a,i)
               * equation for bunches of (b,j) column variables...
               */

              /*
               * J_S_T
               */

              var = TEMPERATURE;
              if (pd->v[pg->imtrx][var]) {
                pvar = upd->vp[pg->imtrx][var];
                for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                  phi_j = bf[var]->phi[j];

                  mass = 0.;

                  if (pd->TimeIntegration != STEADY) {
                    if (pd->e[pg->imtrx][eqn] & T_MASS) {
                      mass = s_dot[a][b];
                      mass *= wt_func * d_at_dT[j] * lambda * det_J * wt;
                      mass *= h3;
                      mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
                    }
                  }

                  advection = 0.;
                  if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
                    if (DOUBLE_NONZERO(lambda)) {

                      advection += v_dot_del_s[a][b] - x_dot_del_s[a][b];
                      if (ucwt != 0.)
                        advection -= ucwt * (gt_dot_s[a][b] + s_dot_g[a][b]);
                      if (lcwt != 0.)
                        advection += lcwt * (s_dot_gt[a][b] + g_dot_s[a][b]);

                      advection *= wt_func * d_at_dT[j] * lambda * det_J * wt * h3;
                      advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
                    }
                  }

                  source = 0.;
                  source1 = 0.;
                  if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
                    source = -(g[a][b] + gt[a][b]) * (at * d_mup->T[j] + mup * d_at_dT[j]);

                    if (DOUBLE_NONZERO(alpha)) {
                      source1 -= s_dot_s[a][b] / (mup * mup) * d_mup->T[j];
                      source1 *= lambda * alpha * saramitoCoeff;
                      source += source1;
                    }
                    source *= wt_func * det_J * wt * h3;
                    source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
                  }

                  lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += mass + advection + source;
                }
              }

              /*
               * J_S_v
               */
              for (p = 0; p < WIM; p++) {
                var = VELOCITY1 + p;
                if (pd->v[pg->imtrx][var]) {
                  pvar = upd->vp[pg->imtrx][var];
                  for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                    phi_j = bf[var]->phi[j];

                    mass = 0.;

                    if (pd->TimeIntegration != STEADY) {
                      if (pd->e[pg->imtrx][eqn] & T_MASS) {
                        if (supg != 0.) {
                          mass = supg * supg_terms.supg_tau * phi_j * bf[eqn]->grad_phi[i][p];

                          for (w = 0; w < dim; w++) {
                            mass += supg * supg_terms.d_supg_tau_dv[p][j] * v[w] *
                                    bf[eqn]->grad_phi[i][w];
                          }

                          mass *= s_dot[a][b];
                        }

                        mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)] * at * det_J * wt * h3;
                      }
                    }

                    advection = 0.;

                    if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
                      if (DOUBLE_NONZERO(lambda)) {
                        advection_a = phi_j * (grad_s[p][a][b]);

                        advection_a *= wt_func;

                        advection_b = 0.;
                        /* Petrov-Galerkin term */
                        if (supg != 0.) {

                          advection_b =
                              supg * supg_terms.supg_tau * phi_j * bf[eqn]->grad_phi[i][p];
                          for (w = 0; w < dim; w++) {
                            advection_b += supg * supg_terms.d_supg_tau_dv[p][j] * v[w] *
                                           bf[eqn]->grad_phi[i][w];
                          }

                          advection_b *= R_advection;
                        }

                        advection_c = 0.;
                        if (evss_gradv) {
                          if (pd->CoordinateSystem != CYLINDRICAL) {
                            if (ucwt != 0) {
                              for (k = 0; k < VIM; k++) {
                                advection_c -=
                                    ucwt * (bf[VELOCITY1 + a]->grad_phi_e[j][p][k][a] * s[k][b] +
                                            bf[VELOCITY1 + b]->grad_phi_e[j][p][k][b] * s[a][k]);
                              }
                            }
                            if (lcwt != 0.) {
                              for (k = 0; k < VIM; k++) {
                                advection_c +=
                                    lcwt * (bf[VELOCITY1 + b]->grad_phi_e[j][p][b][k] * s[a][k] +
                                            bf[VELOCITY1 + a]->grad_phi_e[j][p][a][k] * s[k][b]);
                              }
                            }
                          } else {
                            if (ucwt != 0) {
                              for (k = 0; k < VIM; k++) {
                                advection_c -=
                                    ucwt * (bf[VELOCITY1]->grad_phi_e[j][p][k][a] * s[k][b] +
                                            bf[VELOCITY1]->grad_phi_e[j][p][k][b] * s[a][k]);
                              }
                            }
                            if (lcwt != 0.) {
                              for (k = 0; k < VIM; k++) {
                                advection_c +=
                                    lcwt * (bf[VELOCITY1]->grad_phi_e[j][p][b][k] * s[a][k] +
                                            bf[VELOCITY1]->grad_phi_e[j][p][a][k] * s[k][b]);
                              }
                            }
                          }
                          advection_c *= wt_func;
                        }

                        advection = advection_a + advection_b + advection_c;
                        advection *= at * det_J * wt * h3;
                        advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
                      }
                    }

                    diffusion = 0.;
                    if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {
                      diffusion *= det_J * wt * h3;
                      diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
                    }

                    source = 0.;

                    if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
                      source = 0;
                      source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
                    }

                    lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += mass + advection + diffusion + source;
                  }
                }
              }

              /*
               * J_S_c
               */
              var = MASS_FRACTION;
              if (pd->v[pg->imtrx][var]) {
                for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                  phi_j = bf[var]->phi[j];

                  for (w = 0; w < pd->Num_Species_Eqn; w++) {

                    source = 0.;

                    if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
                      source_a = -at * d_mup->C[w][j] * (g[a][b] + gt[a][b]);

                      source_b = 0.;
                      if (DOUBLE_NONZERO(alpha)) {
                        source_b -= s_dot_s[a][b] / (mup * mup);
                        source_b *= alpha * lambda * saramitoCoeff * d_mup->C[w][j];
                      }
                      source = source_a + source_b;
                      source *= wt_func * det_J * wt * h3;
                      source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
                    }

                    if (w > 1) {
                      GOMA_EH(GOMA_ERROR, "Need more arrays for each species.");
                    }

                    lec->J[LEC_J_INDEX(peqn, MAX_PROB_VAR + w, i, j)] += source;
                  }
                }
              }

              /*
               * J_S_P
               */
              var = PRESSURE;
              if (pd->v[pg->imtrx][var]) {
                pvar = upd->vp[pg->imtrx][var];
                for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                  phi_j = bf[var]->phi[j];

                  source = 0.;
                  if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
                    source_a += -at * d_mup->P[j] * (g[a][b] + gt[a][b]);

                    source_b = 0.;
                    if (DOUBLE_NONZERO(alpha)) {
                      source_b -= (s_dot_s[a][b] / (mup * mup));
                      source_b *= d_mup->P[j] * alpha * lambda * saramitoCoeff;
                    }
                    source = source_a + source_b;
                    source *= wt_func * det_J * wt * h3;
                    source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
                  }

                  lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += source;
                }
              }

              /*
               * J_S_d
               */
              for (p = 0; p < dim; p++) {
                var = MESH_DISPLACEMENT1 + p;
                if (pd->v[pg->imtrx][var]) {
                  pvar = upd->vp[pg->imtrx][var];
                  for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                    phi_j = bf[var]->phi[j];
                    d_det_J_dmesh_pj = bf[eqn]->d_det_J_dm[p][j];
                    dh3dmesh_pj = fv->dh3dq[p] * bf[var]->phi[j];
                    d_mup_dmesh_pj = d_mup->X[p][j];

                    mass = 0.;
                    mass_a = 0.;
                    mass_b = 0.;
                    if (pd->TimeIntegration != STEADY) {
                      if (pd->e[pg->imtrx][eqn] & T_MASS) {
                        mass_a = s_dot[a][b];
                        mass_a *= wt_func * (d_det_J_dmesh_pj * h3 + det_J * dh3dmesh_pj);

                        if (supg != 0.) {
                          for (w = 0; w < dim; w++) {
                            mass_b += supg * (supg_terms.supg_tau * v[w] *
                                                  bf[eqn]->d_grad_phi_dmesh[i][w][p][j] +
                                              supg_terms.d_supg_tau_dX[p][j] * v[w] *
                                                  bf[eqn]->grad_phi[i][w]);
                          }
                          mass_b *= s_dot[a][b] * h3 * det_J;
                        }

                        mass = mass_a + mass_b;
                        mass *= at * lambda * wt * pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
                      }
                    }

                    advection = 0.;

                    if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
                      if (DOUBLE_NONZERO(lambda)) {
                        /*
                         * Four parts:
                         *    advection_a =
                         *    	Int ( ea.(v-xdot).d(Vv)/dmesh h3 |Jv| )
                         *
                         *    advection_b =
                         *  (i)	Int ( ea.(v-xdot).Vv h3 d(|Jv|)/dmesh )
                         *  (ii)  Int ( ea.(v-xdot).d(Vv)/dmesh h3 |Jv| )
                         *  (iii) Int ( ea.(v-xdot).Vv dh3/dmesh |Jv|   )
                         *
                         * For unsteady problems, we have an
                         * additional term
                         *
                         *    advection_c =
                         *    	Int ( ea.d(v-xdot)/dmesh.Vv h3 |Jv| )
                         */

                        advection_a = R_advection;

                        advection_a *= wt_func * (d_det_J_dmesh_pj * h3 + det_J * dh3dmesh_pj);

                        d_vdotdels_dm = 0.;
                        for (q = 0; q < WIM; q++) {
                          d_vdotdels_dm += (v[q] - x_dot[q]) * d_grad_s_dmesh[q][a][b][p][j];
                        }

                        advection_b = d_vdotdels_dm;
                        advection_b *= wt_func * det_J * h3;

                        advection_c = 0.;
                        if (pd->TimeIntegration != STEADY) {
                          if (pd->e[pg->imtrx][eqn] & T_MASS) {
                            d_xdotdels_dm = (1. + 2. * tt) * phi_j / dt * grad_s[p][a][b];

                            advection_c -= d_xdotdels_dm;

                            advection_c *= wt_func * h3 * det_J;
                          }
                        }

                        advection_d = 0.;
                        if (supg != 0.) {
                          for (w = 0; w < dim; w++) {
                            advection_d += supg * (supg_terms.supg_tau * v[w] *
                                                       bf[eqn]->d_grad_phi_dmesh[i][w][p][j] +
                                                   supg_terms.d_supg_tau_dX[p][j] * v[w] *
                                                       bf[eqn]->grad_phi[i][w]);
                          }

                          advection_d *= (R_advection)*det_J * h3;
                        }

                        advection = advection_a + advection_b + advection_c + advection_d;

                        advection *= wt * at * lambda * pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
                      }
                    }

                    diffusion = 0.;
                    if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {
                      diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
                    }

                    /*
                     * Source term...
                     */

                    source = 0.;

                    if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
                      source_a = R_source;
                      source_b = -at * (g[a][b] + gt[a][b]);

                      if (DOUBLE_NONZERO(alpha)) {
                        source_b += -s_dot_s[a][b] / (mup * mup) * alpha * lambda * saramitoCoeff;
                      }

                      source_a *= wt_func * (d_det_J_dmesh_pj * h3 + det_J * dh3dmesh_pj);

                      source_b *= wt_func * det_J * h3 * d_mup_dmesh_pj;

                      source_c = 0.;
                      if (supg != 0.) {
                        for (w = 0; w < dim; w++) {
                          source_c +=
                              supg *
                              (supg_terms.supg_tau * v[w] * bf[eqn]->d_grad_phi_dmesh[i][w][p][j] +
                               supg_terms.d_supg_tau_dX[p][j] * v[w] * bf[eqn]->grad_phi[i][w]);
                        }
                        source_c *= R_source * det_J * h3;
                      }

                      source = source_a + source_b + source_c;

                      source *= wt * pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
                    }

                    lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += mass + advection + diffusion + source;
                  }
                }
              }

              /*
               * J_S_G
               */
              if (evss_gradv == 0) {
                for (p = 0; p < VIM; p++) {
                  for (q = 0; q < VIM; q++) {
                    var = v_g[p][q];

                    if (pd->v[pg->imtrx][var]) {
                      pvar = upd->vp[pg->imtrx][var];
                      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                        phi_j = bf[var]->phi[j];
                        advection = 0.;
                        advection_a = 0.;
                        if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
                          if (DOUBLE_NONZERO(lambda)) {

                            advection -= ucwt * (s[p][b] * (double)delta(a, q) +
                                                 s[a][p] * (double)delta(b, q));
                            advection += lcwt * (s[a][q] * (double)delta(p, b) +
                                                 s[q][b] * (double)delta(a, p));

                            advection *= phi_j * h3 * det_J;

                            advection *=
                                wt_func * wt * at * pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
                          }
                        }

                        /*
                         * Diffusion...
                         */

                        diffusion = 0.;

                        if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {
                          diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
                        }

                        /*
                         * Source term...
                         */

                        source = 0.;

                        lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += advection + diffusion + source;
                      }
                    }
                  }
                }
              }

              /*
               * J_S_F
               */
              var = FILL;
              if (pd->v[pg->imtrx][var]) {
                pvar = upd->vp[pg->imtrx][var];
                for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                  phi_j = bf[var]->phi[j];

                  mass = 0.;

                  if (pd->TimeIntegration != STEADY) {
                    if (pd->e[pg->imtrx][eqn] & T_MASS) {

                      mass = s_dot[a][b];
                      mass *= d_lam->F[j];
                      mass *= wt_func * at * det_J * wt;
                      mass *= h3;
                      mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
                    }
                  }

                  advection = 0.;

                  if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
                    if (d_lam->F[j] != 0.) {

                      advection += v_dot_del_s[a][b] - x_dot_del_s[a][b];
                      if (ucwt != 0.)
                        advection -= ucwt * (gt_dot_s[a][b] + s_dot_g[a][b]);
                      if (lcwt != 0.)
                        advection += lcwt * (s_dot_gt[a][b] + g_dot_s[a][b]);

                      advection *= d_lam->F[j];
                      advection *= wt_func * at * det_J * wt * h3;
                      advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
                    }
                  }

                  diffusion = 0.;

                  if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {
                    /* add SU term in here when appropriate */

                    diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
                  }

                  source = 0.;

                  if (pd->e[pg->imtrx][eqn] & T_SOURCE) {

                    double invmup = 1 / mup;
                    // PTT
                    if (eps != 0) {
                      // product rule + exponential
                      source += Z *
                                ((lambda * trace * d_eps_dF[j] * invmup) +
                                 (d_lam->F[j] * trace * eps * invmup) -
                                 (lambda * trace * eps * d_mup->F[j] * invmup * invmup)) *
                                s[a][b];
                    }

                    source += -at * d_mup->F[j] * (g[a][b] + gt[a][b]);

                    // Giesekus
                    if (alpha != 0.) {
                      source += s_dot_s[a][b] *
                                (-alpha * lambda * d_mup->F[j] * invmup * invmup +
                                 d_alpha_dF[j] * lambda * invmup + alpha * d_lam->F[j] * invmup);
                    }

                    source *= wt_func * det_J * h3 * wt;

                    source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
                  }

                  lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += mass + advection + diffusion + source;
                }
              }

              /*
               * J_S_S
               */
              for (p = 0; p < VIM; p++) {
                for (q = 0; q < VIM; q++) {
                  var = v_s[mode][p][q];

                  if (pd->v[pg->imtrx][var]) {
                    pvar = upd->vp[pg->imtrx][var];
                    for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                      phi_j = bf[var]->phi[j];
                      mass = 0.;
                      if (pd->TimeIntegration != STEADY) {
                        if (pd->e[pg->imtrx][eqn] & T_MASS) {
                          mass = (1. + 2. * tt) * phi_j / dt * (double)delta(a, p) *
                                 (double)delta(b, q);
                          mass *= h3 * det_J;
                          mass *= wt_func * at * wt * pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
                        }
                      }

                      advection = 0.;

                      if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
                        if (DOUBLE_NONZERO(lambda)) {
                          if ((a == p) && (b == q)) {
                            for (r = 0; r < WIM; r++) {
                              advection += (v[r] - x_dot[r]) * bf[var]->grad_phi[j][r];
                            }
                          }
                          advection -=
                              phi_j * ucwt *
                              (gt[a][p] * (double)delta(b, q) + g[q][b] * (double)delta(a, p));
                          advection +=
                              phi_j * lcwt *
                              (gt[q][b] * (double)delta(p, a) + g[a][p] * (double)delta(q, b));

                          advection *= h3 * det_J;

                          advection *=
                              wt_func * wt * at * pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
                        }
                      }

                      /*
                       * Diffusion...
                       */

                      diffusion = 0.;

                      if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {
                        diffusion *= det_J * wt * h3;
                        diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
                      }

                      /*
                       * Source term...
                       */

                      source = 0.;

                      if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
                        source -= d_source_term_dc[a][b][p][q];
                        source *= (1 / lambda) * phi_j * det_J * h3 * wt_func * wt *
                                  pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
                      }

                      lec->J[LEC_J_INDEX(peqn, pvar, i, j)] +=
                          mass + advection + diffusion + source;
                    }
                  }
                }
              }
            }
          }
        }
      }
    } /* End Assemble Jacobian */
  }   /* End loop over modes */

  return (status);
}
/* assemble_stress_vesolid -- assemble terms (Residual &| Jacobian) for solid stress eqns
 *
 * in:
 * 	ei -- pointer to Element Indeces	structure
 *	pd -- pointer to Problem Description	structure
 *	af -- pointer to Action Flag		structure
 *	bf -- pointer to Basis Function		structure
 *	fv -- pointer to Field Variable		structure
 *	cr -- pointer to Constitutive Relation	structure
 *	md -- pointer to Mesh Derivative	structure
 *	me -- pointer to Material Entity	structure
 * 	ija -- vector of pointers into the a matrix
 * 	a  -- global Jacobian matrix
 * 	R  -- global residual vector
 *
 * out:
 *	a   -- gets loaded up with proper contribution
 *	lec -- gets loaded up with local contributions to resid, Jacobian
 * 	r  -- residual RHS vector
 *
 * Created:	Wed Dec  8 14:03:06 MST 1993 pasacki@sandia.gov
 *
 * Revised:	Sun Feb 27 06:53:12 MST 1994 pasacki@sandia.gov
 *
 *
 * Note: currently we do a "double load" into the addresses in the global
 *       "a" matrix and resid vector held in "esp", as well as into the
 *       local accumulators in "lec".
 *
 */

int assemble_stress_vesolid(const double tt,    /* parameter to vary time integration from
                                  explicit (tt = 1) to implicit (tt = 0) */
                            const double dt,    /* current time step size */
                            const int ielem,    /* current element number */
                            const int ip,       /* current integration point */
                            const int ip_total) /* total gauss integration points */
{
  int dim, p, q, r, a, b, w = -1;

  int mode; /*counter for viscoelastic modes */
  int transient_run = pd->TimeIntegration != STEADY;
  int mass_on;
  int advection_on = 0;
  int source_on = 0;
  int diffusion_on = 0;

  int err, eqn, var;
  int peqn, pvar;

  int i, j, status, imtrx;

  double mass_etm, advection_etm, diffusion_etm, source_etm;

  double vconv[MAX_PDIM];     /*Calculated convection velocity */
  double vconv_old[MAX_PDIM]; /*Calculated convection velocity at previous time*/
  CONVECTION_VELOCITY_DEPENDENCE_STRUCT d_vconv_struct;
  CONVECTION_VELOCITY_DEPENDENCE_STRUCT *d_vconv = &d_vconv_struct;
  double v_dot_del_s[DIM][DIM];
  double x_dot_del_s[DIM][DIM];

  double x_dot[DIM];  /* current position field derivative wrt time. */
  double h3;          /* Volume element (scale factors). */
  double dh3dmesh_pj; /* Sensitivity to (p,j) mesh dof. */

  double det_J; /* determinant of element Jacobian */

  double d_det_J_dmesh_pj; /* for specific (p,j) mesh dof */

  double mass; /* For terms and their derivatives */
  double advection;
  double diffusion;
  double source;
  double d_area;

  /*
   * Petrov-Galerkin weighting functions for i-th and ab-th stress residuals
   * and some of their derivatives...
   */

  dbl wt_func;

  /*
   * Interpolation functions for variables and some of their derivatives.
   */

  dbl phi_j;
  dbl wt;

  /* Variables for stress */

  int R_s[MAX_MODES][DIM][DIM];
  int v_s[MAX_MODES][DIM][DIM];

  double s[DIM][DIM];     /* stress tensor */
  double s_dot[DIM][DIM]; /* stress tensor from last time step */
  double grad_s[DIM][DIM][DIM];
  double d_grad_s_dmesh[DIM][DIM][DIM][DIM]
                       [MDE]; /* derivative of grad of stress tensor for mode ve_mode */

  double TT[DIM][DIM];                    /* Mesh stress tensor... */
  double dTT_dx[DIM][DIM][DIM][MDE];      /* Sensitivity of stress tensor...
                          to nodal displacements */
  double dTT_dp[DIM][DIM][MDE];           /* Sensitivity of stress tensor...
                                  to nodal pressure*/
  double dTT_dc[DIM][DIM][MAX_CONC][MDE]; /* Sensitivity of stress tensor...
                          to nodal concentration */
  double dTT_dp_liq[DIM][DIM][MDE];       /* Sensitivity of stress tensor...
                                           to nodal porous liquid pressure*/
  double dTT_dp_gas[DIM][DIM][MDE];       /* Sensitivity of stress tensor...
                                           to nodal porous gas pressure*/
  double dTT_dporosity[DIM][DIM][MDE];    /* Sensitivity of stress tensor...
                                        to nodal porosity*/
  double dTT_dsink_mass[DIM][DIM][MDE];   /* Sensitivity of stress tensor...
                                          to nodal sink_mass*/
  double dTT_dT[DIM][DIM][MDE];           /* Sensitivity of stress tensor...
                                      to temperature*/
  double dTT_dmax_strain[DIM][DIM][MDE];  /* Sensitivity of stress tensor...
                             to max_strain*/
  double dTT_dcur_strain[DIM][DIM][MDE];  /* Sensitivity of stress tensor...
                             to cur_strain*/
  double lame_mu = elc->lame_mu;
  double lame_lambda = elc->lame_lambda;

  /* constitutive equation parameters */
  double lambda = 0; /* polymer relaxation constant */

  status = 0;

  eqn = R_STRESS11;

  /*
   * Bail out fast if there's nothing to do...
   */

  if (!pd->e[pg->imtrx][eqn]) {
    return (status);
  }

  /*
   * Unpack variables from structures for local convenience...
   */

  dim = pd->Num_Dim;
  wt = fv->wt;
  if (pd->e[pg->imtrx][R_MESH1]) {
    det_J = bf[R_MESH1]->detJ;
  } else {
    det_J = bf[eqn]->detJ; /* Really, ought to be mesh eqn. */
  }
  h3 = fv->h3; /* Differential volume element (scales). */
  d_area = det_J * wt * h3;
  mass_on = pd->e[pg->imtrx][eqn] & T_MASS;
  advection_on = pd->e[pg->imtrx][eqn] & T_ADVECTION;
  diffusion_on = pd->e[pg->imtrx][eqn] & T_DIFFUSION;
  source_on = pd->e[pg->imtrx][eqn] & T_SOURCE;

  mass_etm = pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
  advection_etm = pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
  diffusion_etm = pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
  source_etm = pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
  /*
   * Get the deformation gradients and tensors if needed
   */
  err = belly_flop(lame_mu);
  GOMA_EH(err, "error in belly flop");
  if (err == 2)
    return (err);

  /*
   * Field variables...
   */
  /* Calculate inertia of mesh if required. PRS side note: no sensitivity
   * with respect to porous media variables here for single compont pore liquids.
   * eventually there may be sensitivity for diffusion in gas phase to p_gas */
  if (cr->MeshMotion == LAGRANGIAN || cr->MeshMotion == DYNAMIC_LAGRANGIAN) {
    err = get_convection_velocity(vconv, vconv_old, d_vconv, dt, tt);
  } else /* No inertia in an Arbitrary Mesh */
  {
    memset(vconv, 0, sizeof(double) * MAX_PDIM);
    for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
      if (pd->v[imtrx][MESH_DISPLACEMENT1])
        memset(d_vconv->X, 0, DIM * DIM * MDE * sizeof(dbl));

      if (pd->v[imtrx][VELOCITY1] || pd->v[imtrx][POR_LIQ_PRES])
        memset(d_vconv->v, 0, DIM * DIM * MDE * sizeof(dbl));

      if (pd->v[imtrx][MASS_FRACTION] || pd->v[imtrx][POR_LIQ_PRES])
        memset(d_vconv->C, 0, DIM * MAX_CONC * MDE * sizeof(dbl));

      if (pd->v[pg->imtrx][TEMPERATURE])
        memset(d_vconv->T, 0, DIM * MDE * sizeof(dbl));
    }
  }

  for (a = 0; a < WIM; a++) {
    /* note, these are zero for steady calculations */
    x_dot[a] = 0.0;
    for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
      if (transient_run && pd->v[imtrx][MESH_DISPLACEMENT1 + a]) {
        x_dot[a] = fv_dot->x[a];
      }
    }
  }

  /*
   * Total mesh stress tensor...
   */
  /* initialize some arrays */
  memset(TT, 0, sizeof(double) * DIM * DIM);
  if (af->Assemble_Jacobian) {
    memset(dTT_dx, 0, sizeof(double) * DIM * DIM * DIM * MDE);
    memset(dTT_dp, 0, sizeof(double) * DIM * DIM * MDE);
    memset(dTT_dc, 0, sizeof(double) * DIM * DIM * MAX_CONC * MDE);
    memset(dTT_dp_liq, 0, sizeof(double) * DIM * DIM * MDE);
    memset(dTT_dp_gas, 0, sizeof(double) * DIM * DIM * MDE);
    memset(dTT_dporosity, 0, sizeof(double) * DIM * DIM * MDE);
    memset(dTT_dsink_mass, 0, sizeof(double) * DIM * DIM * MDE);
    memset(dTT_dT, 0, sizeof(double) * DIM * DIM * MDE);
    memset(dTT_dmax_strain, 0, sizeof(double) * DIM * DIM * MDE);
    memset(dTT_dcur_strain, 0, sizeof(double) * DIM * DIM * MDE);
  }

  err = mesh_stress_tensor(TT, dTT_dx, dTT_dp, dTT_dc, dTT_dp_liq, dTT_dp_gas, dTT_dporosity,
                           dTT_dsink_mass, dTT_dT, dTT_dmax_strain, dTT_dcur_strain, lame_mu,
                           lame_lambda, dt, ielem, ip, ip_total);

  /* get tensor dot products for future use */

  (void)stress_eqn_pointer(v_s);
  (void)stress_eqn_pointer(R_s);

  /* Begin loop over modes */
  if (vn->modes != 1)
    GOMA_EH(GOMA_ERROR, "VE solid only set up for 1 mode at present!\n");
  for (mode = 0; mode < vn->modes; mode++) {

    load_modal_pointers(mode, tt, dt, s, s_dot, grad_s, d_grad_s_dmesh);

    /* precalculate advective terms of form (v dot del tensor)*/
    memset(v_dot_del_s, 0, sizeof(double) * DIM * DIM);
    memset(x_dot_del_s, 0, sizeof(double) * DIM * DIM);
    for (a = 0; a < VIM; a++) {
      for (b = 0; b < VIM; b++) {
        for (q = 0; q < WIM; q++) {
          v_dot_del_s[a][b] += vconv[q] * grad_s[q][a][b];
          x_dot_del_s[a][b] += x_dot[q] * grad_s[q][a][b];
        }
      }
    }

    /*
     * Stress tensor...(Note "anti-BSL" sign convention on deviatoric stress)
     */

    /* get time constant */
    if (elc->solid_retard_model == CONSTANT) {
      lambda = elc->solid_retardation;
    }

    /*
     * Residuals_________________________________________________________________
     */

    if (af->Assemble_Residual) {
      /*
       * Assemble each component "ab" of the polymer stress equation...
       */
      for (a = 0; a < VIM; a++) {
        for (b = 0; b < VIM; b++) {

          if (a <= b) /* since the stress tensor is symmetric, only assemble the upper half */
          {
            eqn = R_s[mode][a][b];
            peqn = upd->ep[pg->imtrx][eqn];

            /*
             * In the element, there will be contributions to this many equations
             * based on the number of degrees of freedom...
             */

            for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
              wt_func = bf[eqn]->phi[i];

              mass = 0.;
              if (transient_run && mass_on) {
                mass = s_dot[a][b] * wt_func * lambda * mass_etm;
              }

              advection = 0.;
              if (advection_on && DOUBLE_NONZERO(lambda)) {
                advection =
                    (v_dot_del_s[a][b] - x_dot_del_s[a][b]) * wt_func * lambda * advection_etm;
              }

              diffusion = 0.;
              if (diffusion_on) {
                /* add SU term in here when appropriate */
                diffusion = -TT[a][b] * wt_func * diffusion_etm;
              }

              /*
               * Source term...
               */

              source = 0.;
              if (source_on) {
                source = s[a][b] * wt_func * source_etm;
              }

              /*
               * Add contributions to this residual (globally into Resid, and
               * locally into an accumulator)
               */

              lec->R[LEC_R_INDEX(peqn, i)] += (mass + advection + diffusion + source) * d_area;
            }
          }
        }
      }
    }

    /*
     * Jacobian terms...
     */

    if (af->Assemble_Jacobian) {
      for (a = 0; a < VIM; a++) {
        for (b = 0; b < VIM; b++) {
          if (a <= b) /* since the stress tensor is symmetric, only assemble the upper half */
          {

            eqn = R_s[mode][a][b];
            peqn = upd->ep[pg->imtrx][eqn];

            for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {

              wt_func = bf[eqn]->phi[i];

              /*
               * Set up some preliminaries that are needed for the (a,i)
               * equation for bunches of (b,j) column variables...
               */

              /*
               * J_S_T
               */

              var = TEMPERATURE;
              if (pd->v[pg->imtrx][var]) {
                pvar = upd->vp[pg->imtrx][var];
                for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                  phi_j = bf[var]->phi[j];

                  mass = 0.;

                  advection = 0.;
                  if (advection_on && DOUBLE_NONZERO(lambda)) {
                    for (q = 0; q < WIM; q++) {
                      advection += d_vconv->T[q][j] * grad_s[q][a][b];
                    }
                    advection *= wt_func * lambda * advection_etm;
                  }

                  diffusion = 0.;
                  if (diffusion_on) {
                    diffusion = -dTT_dT[a][b][j] * wt_func * diffusion_etm;
                  }

                  source = 0.;

                  lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += (advection + diffusion) * d_area;
                }
              }

              /*
               * J_S_c
               */
              var = MASS_FRACTION;
              if (pd->v[pg->imtrx][var]) {
                pvar = MAX_PROB_VAR + w;
                for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                  phi_j = bf[var]->phi[j];

                  for (w = 0; w < pd->Num_Species_Eqn; w++) {

                    mass = 0.;

                    advection = 0.;
                    if (advection_on && DOUBLE_NONZERO(lambda)) {
                      for (q = 0; q < WIM; q++) {
                        advection += d_vconv->C[q][w][j] * grad_s[q][a][b];
                      }
                      advection *= wt_func * lambda * advection_etm;
                    }

                    diffusion = 0.;
                    if (diffusion_on) {
                      diffusion = -dTT_dc[a][b][w][j] * wt_func * diffusion_etm;
                    }

                    source = 0.;

                    lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += (advection + diffusion) * d_area;
                  }
                }
              }

              /*
               * J_S_P
               */
              var = PRESSURE;
              if (pd->v[pg->imtrx][var]) {
                pvar = upd->vp[pg->imtrx][var];
                for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                  phi_j = bf[var]->phi[j];

                  mass = 0.;

                  advection = 0.;

                  diffusion = 0.;
                  if (diffusion_on) {
                    diffusion = -dTT_dp[a][b][j] * wt_func * diffusion_etm;
                  }

                  source = 0.;

                  lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += diffusion * d_area;
                }
              }

              /*
               * J_S_p_liq
               */
              var = POR_LIQ_PRES;
              if (pd->v[pg->imtrx][var] &&
                  (cr->MeshMotion == LAGRANGIAN || cr->MeshMotion == DYNAMIC_LAGRANGIAN)) {
                pvar = upd->vp[pg->imtrx][var];
                for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                  phi_j = bf[var]->phi[j];

                  mass = 0.;

                  advection = 0.;

                  diffusion = 0.;
                  if (diffusion_on) {
                    diffusion = -dTT_dp_liq[a][b][j] * wt_func * diffusion_etm;
                  }

                  source = 0.;

                  lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += diffusion * d_area;
                }
              }

              /*
               * J_S_p_gas
               */
              var = POR_GAS_PRES;
              if (pd->v[pg->imtrx][var] &&
                  (cr->MeshMotion == LAGRANGIAN || cr->MeshMotion == DYNAMIC_LAGRANGIAN)) {
                pvar = upd->vp[pg->imtrx][var];
                for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                  phi_j = bf[var]->phi[j];

                  mass = 0.;

                  advection = 0.;

                  diffusion = 0.;
                  if (diffusion_on) {
                    diffusion = -dTT_dp_gas[a][b][j] * wt_func * diffusion_etm;
                  }

                  source = 0.;

                  lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += diffusion * d_area;
                }
              }

              /*
               * J_S_porosity
               */
              var = POR_POROSITY;
              if (pd->v[pg->imtrx][var] &&
                  (cr->MeshMotion == LAGRANGIAN || cr->MeshMotion == DYNAMIC_LAGRANGIAN)) {
                pvar = upd->vp[pg->imtrx][var];
                for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                  phi_j = bf[var]->phi[j];

                  mass = 0.;

                  advection = 0.;

                  diffusion = 0.;
                  if (diffusion_on) {
                    diffusion = -dTT_dporosity[a][b][j] * wt_func * diffusion_etm;
                  }

                  source = 0.;

                  lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += diffusion * d_area;
                }
              }

              /*
               * J_S_sink_mass
               */
              var = POR_SINK_MASS;
              if (pd->v[pg->imtrx][var] &&
                  (cr->MeshMotion == LAGRANGIAN || cr->MeshMotion == DYNAMIC_LAGRANGIAN)) {
                pvar = upd->vp[pg->imtrx][var];
                for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                  phi_j = bf[var]->phi[j];

                  mass = 0.;

                  advection = 0.;

                  diffusion = 0.;
                  if (diffusion_on) {
                    diffusion = -dTT_dsink_mass[a][b][j] * wt_func * diffusion_etm;
                  }

                  source = 0.;

                  lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += diffusion * d_area;
                }
              }

              /*
               * J_S_max_strain
               */
              var = MAX_STRAIN;
              if (pd->v[pg->imtrx][var] &&
                  (cr->MeshMotion == LAGRANGIAN || cr->MeshMotion == DYNAMIC_LAGRANGIAN)) {
                pvar = upd->vp[pg->imtrx][var];
                for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                  phi_j = bf[var]->phi[j];

                  mass = 0.;

                  advection = 0.;

                  diffusion = 0.;
                  if (diffusion_on) {
                    diffusion = -dTT_dmax_strain[a][b][j] * wt_func * diffusion_etm;
                  }

                  source = 0.;

                  lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += diffusion * d_area;
                }
              }

              /*
               * J_S_cur_strain
               */
              var = CUR_STRAIN;
              if (pd->v[pg->imtrx][var] &&
                  (cr->MeshMotion == LAGRANGIAN || cr->MeshMotion == DYNAMIC_LAGRANGIAN)) {
                pvar = upd->vp[pg->imtrx][var];
                for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                  phi_j = bf[var]->phi[j];

                  mass = 0.;

                  advection = 0.;

                  diffusion = 0.;
                  if (diffusion_on) {
                    diffusion = -dTT_dcur_strain[a][b][j] * wt_func * diffusion_etm;
                  }

                  source = 0.;

                  lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += diffusion * d_area;
                }
              }

              /*
               * J_S_d
               */
              for (p = 0; p < dim; p++) {
                var = MESH_DISPLACEMENT1 + p;
                if (pd->v[pg->imtrx][var]) {
                  pvar = upd->vp[pg->imtrx][var];
                  for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                    phi_j = bf[var]->phi[j];
                    d_det_J_dmesh_pj = bf[var]->d_det_J_dm[p][j];
                    dh3dmesh_pj = fv->dh3dq[p] * bf[var]->phi[j];

                    mass = 0.;
                    if (transient_run && mass_on) {
                      mass = s_dot[a][b] * wt_func * lambda * mass_etm;
                      mass *= wt * (d_det_J_dmesh_pj * h3 + det_J * dh3dmesh_pj);
                    }

                    advection = 0.;
                    if (advection_on && DOUBLE_NONZERO(lambda)) {
                      for (q = 0; q < WIM; q++) {
                        advection += d_vconv->X[q][p][j] * grad_s[q][a][b];
                        advection += (vconv[q] - x_dot[q]) * d_grad_s_dmesh[q][a][b][p][j];
                      }
                      advection *= wt_func * d_area * advection_etm;
                      advection += (v_dot_del_s[a][b] - x_dot_del_s[a][b]) * wt_func * lambda *
                                   advection_etm * wt *
                                   (d_det_J_dmesh_pj * h3 + det_J * dh3dmesh_pj);
                    }

                    diffusion = 0.;
                    if (diffusion_on) {
                      diffusion = -TT[a][b] * wt_func * diffusion_etm;
                      diffusion *= wt * (d_det_J_dmesh_pj * h3 + det_J * dh3dmesh_pj);
                      diffusion += -dTT_dx[a][b][p][j] * wt_func * d_area * diffusion_etm;
                    }

                    source = 0.;
                    if (source_on) {
                      source = s[a][b] * wt_func * source_etm;
                      source *= wt * (d_det_J_dmesh_pj * h3 + det_J * dh3dmesh_pj);
                    }

                    lec->J[LEC_J_INDEX(peqn, pvar, i, j)] +=
                        (mass + advection + diffusion + source);
                  }
                }
              }

              /*
               * J_S_S
               */
              var = v_s[mode][a][b];

              if (pd->v[pg->imtrx][var]) {
                pvar = upd->vp[pg->imtrx][var];
                for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                  phi_j = bf[var]->phi[j];

                  mass = 0.;
                  if (transient_run && mass_on) {
                    mass = (1. + 2. * tt) * phi_j / dt;
                    mass *= d_area * wt_func * lambda * mass_etm;
                  }

                  advection = 0.;
                  if (advection_on && DOUBLE_NONZERO(lambda)) {
                    for (r = 0; r < WIM; r++) {
                      advection += (vconv[r] - x_dot[r]) * bf[var]->grad_phi[j][r];
                    }
                    advection *= d_area * wt_func * lambda * advection_etm;
                  }

                  diffusion = 0.;

                  source = 0.;
                  if (source_on) {
                    source = phi_j * d_area * wt_func * source_etm;
                  }

                  lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += mass + advection + diffusion + source;
                }
              }
            }
          }
        }
      }
    }
  }

  return (status);
}
/*****************************************************************************/
/* END OF FILE mm_fill_stress.c */
/*****************************************************************************/
