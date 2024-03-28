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

#include "mm_fill_energy.h"
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
#include "mm_fill_ls_capillary_bcs.h"
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
/* assemble_energy -- assemble terms (Residual &| Jacobian) for energy eqns
 *
 * in:
 * 	ei -- pointer to Element Indeces	structure
 *	pd -- pointer to Problem Description	structure
 *	af -- pointer to Action Flag		structure
 *	bf -- pointer to Basis Function		structure
 *	fv -- pointer to Field Variable		structure
 *  fv_old -- pointer to old Diet Field Variable	structure
 *  fv_dot -- pointer to dot Diet Field Variable	structure
 *	cr -- pointer to Constitutive Relation	structure
 *	md -- pointer to Mesh Derivative	structure
 *	me -- pointer to Material Entity	structure
 *
 * out:
 *	a   -- gets loaded up with proper contribution
 *	lec -- gets loaded up with local contributions to resid, Jacobian
 * 	r   -- residual RHS vector
 *
 * Created:	Thu Mar  3 07:48:01 MST 1994 pasacki@sandia.gov
 *
 * Revised:	9/24/94 by RRR
 *
 * Revised:     10/98 by KSC to enable the calculation of source term due to electrode kinetics as
 *in thermal batteries.
 *
 *
 * Note: currently we do a "double load" into the addresses in the global
 *       "a" matrix and resid vector held in "esp", as well as into the
 *       local accumulators in "lec".
 *
 */
int assemble_energy(
    double time, /* present time value */
    double tt, /* parameter to vary time integration from explicit (tt = 1) to implicit (tt = 0) */
    double dt, /* current time step size */
    const PG_DATA *pg_data) {
  int eqn, var, peqn, pvar, dim, p, a, b, qq, w, i, j, status;

  dbl T_dot; /* Temperature derivative wrt time. */

  dbl q[DIM];                             /* Heat flux vector. */
  HEAT_FLUX_DEPENDENCE_STRUCT d_q_struct; /* Heat flux dependence. */
  HEAT_FLUX_DEPENDENCE_STRUCT *d_q = &d_q_struct;

  dbl rho; /* Density (no variations allowed
              here, for now) */
  // CONDUCTIVITY_DEPENDENCE_STRUCT d_k_struct; /* Thermal conductivity dependence. */
  // CONDUCTIVITY_DEPENDENCE_STRUCT *d_k = &d_k_struct;

  dbl Cp;                                      /* Heat capacity. */
  HEAT_CAPACITY_DEPENDENCE_STRUCT d_Cp_struct; /* Heat capacity dependence. */
  HEAT_CAPACITY_DEPENDENCE_STRUCT *d_Cp = &d_Cp_struct;

  dbl h;                                    /* Heat source. */
  HEAT_SOURCE_DEPENDENCE_STRUCT d_h_struct; /* Heat source dependence. */
  HEAT_SOURCE_DEPENDENCE_STRUCT *d_h = &d_h_struct;

  int v_s[MAX_MODES][DIM][DIM];
  int mode;

  dbl mass; /* For terms and their derivatives */

  dbl advection; /* For terms and their derivatives */
  dbl advection_a;
  dbl advection_b;
  dbl advection_c;
  dbl advection_d;
  dbl advection_e;
  dbl advection_f;

  dbl diffusion;
  dbl diff_a, diff_b, diff_c, diff_d;
  dbl source;

  /*
   * Galerkin weighting functions for i-th energy residuals
   * and some of their derivatives...
   */

  dbl phi_i;
  dbl grad_phi_i[DIM];

  /*
   * Petrov-Galerkin weighting functions for i-th residuals
   * and some of their derivatives...
   */

  dbl wt_func;

  /* SUPG variables */
  dbl h_elem = 0, h_elem_inv = 0, h_elem_deriv = 0, h_elem_inv_deriv = 0.;
  dbl supg, d_wt_func;

  /*
   * Interpolation functions for variables and some of their derivatives.
   */

  dbl phi_j;
  dbl grad_phi_j[DIM];

  dbl h3;          /* Volume element (scale factors). */
  dbl dh3dmesh_bj; /* Sensitivity to (b,j) mesh dof. */

  dbl det_J;

  dbl d_det_J_dmeshbj;        /* for specified (b,j) mesh dof */
  dbl dgrad_phi_i_dmesh[DIM]; /* ditto.  */
  dbl wt;

  dbl *grad_T = fv->grad_T;

  double vconv[MAX_PDIM];     /*Calculated convection velocity */
  double vconv_old[MAX_PDIM]; /*Calculated convection velocity at previous time*/
  CONVECTION_VELOCITY_DEPENDENCE_STRUCT d_vconv_struct;
  CONVECTION_VELOCITY_DEPENDENCE_STRUCT *d_vconv = &d_vconv_struct;

  /* density derivatives */
  DENSITY_DEPENDENCE_STRUCT d_rho_struct; /* density dependence */
  DENSITY_DEPENDENCE_STRUCT *d_rho = &d_rho_struct;

  int err;

  const double *hsquared = pg_data->hsquared;
  const double *vcent = pg_data->v_avg; /* Average element velocity, which is the
          centroid velocity for Q2 and the average
          of the vertices for Q1. It comes from
          the routine "element_velocity." */

  /* initialize grad_phi_j */
  for (i = 0; i < DIM; i++) {
    grad_phi_j[i] = 0;
  }

  /*   static char yo[] = "assemble_energy";*/

  status = 0;

  /*
   * Unpack variables from structures for local convenience...
   */

  dim = pd->Num_Dim;

  eqn = R_ENERGY;

  /*
   * Bail out fast if there's nothing to do...
   */

  if (!pd->e[pg->imtrx][eqn]) {
    return (status);
  }

  wt = fv->wt; /* Gauss point weight. */

  h3 = fv->h3; /* Differential volume element. */

  det_J = bf[eqn]->detJ; /* Really, ought to be mesh eqn. */

  supg = 0.;

  if (mp->Ewt_funcModel == GALERKIN) {
    supg = 0.;
  } else if (mp->Ewt_funcModel == SUPG) {
    if (!pd->gv[R_MOMENTUM1])
      GOMA_EH(GOMA_ERROR, " must have momentum equation velocity field for energy equation "
                          "upwinding. You may want to turn it off");
    supg = mp->Ewt_func;
  }

  (void)stress_eqn_pointer(v_s);

  if (supg != 0.) {
    h_elem = 0.;
    for (p = 0; p < dim; p++) {
      h_elem += vcent[p] * vcent[p] / hsquared[p];
    }
    h_elem = sqrt(h_elem) / 2.;
    if (h_elem == 0.) {
      h_elem_inv = 0.;
    } else {
      h_elem_inv = 1. / h_elem;
    }
  }

  /* end Petrov-Galerkin addition */

  /*
   * Material property constants, etc. Any variations at this Gauss point
   * are calculated with user-defined subroutine options.
   *
   * Properties to consider here are rho, Cp, k, and h.  For now we will
   * take rho as constant.   Cp, h, and k we will allow to vary with temperature,
   * spatial coordinates, and species concentration.
   */

  rho = density(d_rho, time);

  //  conductivity( d_k, time );
  Cp = heat_capacity(d_Cp, time);

  h = heat_source(d_h, time, tt, dt);

  if (pd->TimeIntegration != STEADY) {
    T_dot = fv_dot->T;
  } else {
    T_dot = 0.0;
  }

  heat_flux(q, d_q, time);

  /* get the convection velocity (it's different for arbitrary and
     lagrangian meshes) */
  if (cr->MeshMotion == ARBITRARY || cr->MeshMotion == LAGRANGIAN ||
      cr->MeshMotion == DYNAMIC_LAGRANGIAN) {
    err = get_convection_velocity(vconv, vconv_old, d_vconv, dt, tt);
    GOMA_EH(err, "Error in calculating effective convection velocity");
  } else if (cr->MeshMotion == TOTAL_ALE) {
    err = get_convection_velocity_rs(vconv, vconv_old, d_vconv, dt, tt);
    GOMA_EH(err, "Error in calculating effective convection velocity_rs");
  }

  /*
   * Residuals___________________________________________________________
   */

  if (af->Assemble_Residual) {
    eqn = R_ENERGY;
    peqn = upd->ep[pg->imtrx][eqn];
    var = TEMPERATURE;
    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {

#if 1
      /* this is an optimization for xfem */
      if (xfem != NULL) {
        int xfem_active, extended_dof, base_interp, base_dof;
        xfem_dof_state(i, pd->i[pg->imtrx][eqn], ei[pg->imtrx]->ielem_shape, &xfem_active,
                       &extended_dof, &base_interp, &base_dof);
        if (extended_dof && !xfem_active)
          continue;
      }
#endif
      phi_i = bf[eqn]->phi[i];

      mass = 0.;
      if (pd->TimeIntegration != STEADY) {
        if (pd->e[pg->imtrx][eqn] & T_MASS) {
          mass = T_dot;
          mass *= -phi_i * rho * Cp * det_J * wt;
          mass *= h3;
          mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
        }
      }

      /* only use Petrov Galerkin on advective term - if required */
      wt_func = bf[eqn]->phi[i];
      /* add Petrov-Galerkin terms as necessary */
      if (supg != 0.) {
        for (p = 0; p < dim; p++) {
          wt_func += supg * h_elem_inv * vconv[p] * bf[eqn]->grad_phi[i][p];
        }
      }

      advection = 0.;
      if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {

        for (p = 0; p < VIM; p++) {
          advection += vconv[p] * grad_T[p];
        }

        advection *= -wt_func * rho * Cp * det_J * wt;
        advection *= h3;
        advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
      }

      diffusion = 0.;
      if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {
        for (p = 0; p < VIM; p++) {
          grad_phi_i[p] = bf[eqn]->grad_phi[i][p];
        }

        for (p = 0; p < VIM; p++) {
          diffusion += grad_phi_i[p] * q[p];
        }
        diffusion *= det_J * wt;
        diffusion *= h3;
        diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
      }

      dbl divergence = 0;
      if (mp->Energy_Div_Term) {
        divergence = fv->div_v * fv->T;
        divergence *= -wt_func * rho * Cp * det_J * wt;
        divergence *= h3;
      }

      source = 0.;
      if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
        source += phi_i * h * det_J * wt;
        source *= h3;
        source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
      }

      lec->R[LEC_R_INDEX(peqn, i)] += mass + advection + diffusion + source + divergence;
    }
  }

  /*
   * Jacobian terms...
   */

  if (af->Assemble_Jacobian) {
    eqn = R_ENERGY;
    peqn = upd->ep[pg->imtrx][eqn];
    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
#if 1
      /* this is an optimization for xfem */
      if (xfem != NULL) {
        int xfem_active, extended_dof, base_interp, base_dof;
        xfem_dof_state(i, pd->i[pg->imtrx][eqn], ei[pg->imtrx]->ielem_shape, &xfem_active,
                       &extended_dof, &base_interp, &base_dof);
        if (extended_dof && !xfem_active)
          continue;
      }
#endif
      phi_i = bf[eqn]->phi[i];

      wt_func = bf[eqn]->phi[i];
      /* add Petrov-Galerkin terms as necessary */
      if (supg != 0.) {
        for (p = 0; p < dim; p++) {
          wt_func += supg * h_elem_inv * vconv[p] * bf[eqn]->grad_phi[i][p];
        }
      }

      /*
       * Set up some preliminaries that are needed for the (a,i)
       * equation for bunches of (b,j) column variables...
       */

      for (p = 0; p < VIM; p++) {
        grad_phi_i[p] = bf[eqn]->grad_phi[i][p];
      }

      /*
       * J_e_T
       */
      var = TEMPERATURE;
      if (pd->v[pg->imtrx][var]) {
        pvar = upd->vp[pg->imtrx][var];
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          phi_j = bf[var]->phi[j];

          for (p = 0; p < VIM; p++) {
            grad_phi_j[p] = bf[var]->grad_phi[j][p];
          }

          mass = 0.;
          if (pd->TimeIntegration != STEADY) {
            if (pd->e[pg->imtrx][eqn] & T_MASS) {
              mass = rho * d_Cp->T[j] * T_dot + d_rho->T[j] * Cp * T_dot +
                     rho * Cp * (1 + 2. * tt) * phi_j / dt;
              mass *= -phi_i * det_J * wt;
              mass *= h3;
              mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
            }
          }

          advection = 0.;
          if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
            for (p = 0; p < VIM; p++) {
              advection += rho * d_Cp->T[j] * vconv[p] * grad_T[p] +
                           d_rho->T[j] * Cp * vconv[p] * grad_T[p] +
                           rho * Cp * vconv[p] * grad_phi_j[p] +
                           rho * Cp * d_vconv->T[p][j] * grad_T[p];
            }
            advection *= -wt_func * det_J * wt;
            advection *= h3;
            advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
          }

          diffusion = 0.;
          if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {
            for (p = 0; p < VIM; p++) {
              diffusion += d_q->T[p][j] * grad_phi_i[p];
            }
            diffusion *= det_J * wt;
            diffusion *= h3;
            diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
          }

          dbl divergence = 0;
          if (mp->Energy_Div_Term) {
            divergence += rho * d_Cp->T[j] * fv->div_v * fv->T +
                          d_rho->T[j] * Cp * fv->div_v * fv->T + rho * Cp * fv->div_v * phi_j;
            divergence *= -wt_func * det_J * wt;
            divergence *= h3;
          }

          source = 0.;
          if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
            source += phi_i * d_h->T[j] * det_J * wt;
            source *= h3;
            source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
          }

          lec->J[LEC_J_INDEX(peqn, pvar, i, j)] +=
              mass + advection + diffusion + source + divergence;
        }
      }
      /*
       * J_e_MOM
       */
      for (b = 0; b < MAX_MOMENTS; b++) {
        var = MOMENT0 + b;
        if (pd->v[pg->imtrx][var]) {
          pvar = upd->vp[pg->imtrx][var];
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];

            for (p = 0; p < VIM; p++) {
              grad_phi_j[p] = bf[var]->grad_phi[j][p];
            }

            mass = 0.;
            if (pd->TimeIntegration != STEADY) {
              if (pd->e[pg->imtrx][eqn] & T_MASS) {
                mass = d_rho->moment[b][j] * Cp * T_dot;
                mass *= -phi_i * det_J * wt;
                mass *= h3;
                mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
              }
            }

            advection = 0.;
            if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
              for (p = 0; p < VIM; p++) {
                advection += d_rho->moment[b][j] * Cp * vconv[p] * grad_T[p];
              }
              advection *= -wt_func * det_J * wt;
              advection *= h3;
              advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
            }

            diffusion = 0.;
            if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {
              for (p = 0; p < VIM; p++) {
                diffusion += d_q->moment[b][p][j] * grad_phi_i[p];
              }
              diffusion *= det_J * wt;
              diffusion *= h3;
              diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
            }

            source = 0.;
            if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
              // source += phi_i * d_h->moment[b][j] * det_J * wt;
              source *= h3;
              source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
            }

            lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += mass + advection + diffusion + source;
          }
        }
      }

      var = DENSITY_EQN;
      if (pd->v[pg->imtrx][var]) {
        pvar = upd->vp[pg->imtrx][var];
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          phi_j = bf[var]->phi[j];

          for (p = 0; p < VIM; p++) {
            grad_phi_j[p] = bf[var]->grad_phi[j][p];
          }

          mass = 0.;
          if (pd->TimeIntegration != STEADY) {
            if (pd->e[pg->imtrx][eqn] & T_MASS) {
              mass = d_rho->rho[j] * Cp * T_dot;
              mass *= -phi_i * det_J * wt;
              mass *= h3;
              mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
            }
          }

          advection = 0.;
          if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
            for (p = 0; p < VIM; p++) {
              advection += d_rho->rho[j] * Cp * vconv[p] * grad_T[p];
            }
            advection *= -wt_func * det_J * wt;
            advection *= h3;
            advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
          }

          diffusion = 0.;
          if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {
            for (p = 0; p < VIM; p++) {
              // diffusion += d_q->rho[p][j] * grad_phi_i[p];
            }
            diffusion *= det_J * wt;
            diffusion *= h3;
            diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
          }

          source = 0.;
          if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
            // source += phi_i * d_h->rho[j] * det_J * wt;
            source *= h3;
            source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
          }

          lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += mass + advection + diffusion + source;
        }
      }
      /*
       * J_e_V
       */
      var = VOLTAGE;

      if (pd->v[pg->imtrx][var]) {
        pvar = upd->vp[pg->imtrx][var];
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          phi_j = bf[var]->phi[j];

          for (p = 0; p < dim; p++) {
            grad_phi_j[p] = bf[var]->grad_phi[j][p];
          }

          mass = 0.;

          advection = 0.;

          diffusion = 0.;

          source = 0.;
          if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
            source += phi_i * d_h->V[j] * det_J * wt;
            source *= h3;
            source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
          }

          lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += mass + advection + diffusion + source;
        }
      }

      /*
       * J_e_v
       */
      for (b = 0; b < VIM; b++) {
        var = VELOCITY1 + b;
        if (pd->v[pg->imtrx][var]) {
          pvar = upd->vp[pg->imtrx][var];
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];

            mass = 0.;

            advection = 0.;
            advection_a = 0.;
            if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
              advection_a += wt_func * rho * Cp * d_vconv->v[b][b][j] * grad_T[b];
              advection_a *= -det_J * wt;
              advection_a *= h3;
              advection_a *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
              if (supg != 0.) {
                h_elem_deriv = 0.;
                if (hsquared[b] != 0.) {
                  h_elem_deriv = vcent[b] * pg_data->dv_dnode[b][j] * h_elem_inv / 4. / hsquared[b];
                }
                if (h_elem != 0.)
                  h_elem_inv_deriv = -h_elem_deriv / h_elem / h_elem;
              }
              advection_b = 0.;
              if (supg != 0.) {
                d_wt_func = supg * h_elem_inv * d_vconv->v[b][b][j] * grad_phi_i[b] +
                            supg * h_elem_inv_deriv * vconv[b] * grad_phi_i[b];

                for (p = 0; p < dim; p++) {
                  advection_b += rho * Cp * vconv[p] * grad_T[p];
                }

                advection_b *= d_wt_func;
                advection_b *= -det_J * wt;
                advection_b *= h3;
                advection_b *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
              }
              advection = advection_a + advection_b;
            }
            dbl divergence = 0;
            if (mp->Energy_Div_Term) {
              dbl div_phi_j_e_b = 0.;
              for (p = 0; p < VIM; p++) {
                div_phi_j_e_b += bf[var]->grad_phi_e[j][b][p][p];
              }
              divergence += rho * Cp * fv->T * div_phi_j_e_b;
              divergence *= -wt_func * det_J * wt;
              divergence *= h3;
            }

            diffusion = 0.;
            source = 0.;

            if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
              source += phi_i * d_h->v[b][j] * det_J * wt;
              source *= h3;
              source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
            }

            lec->J[LEC_J_INDEX(peqn, pvar, i, j)] +=
                mass + advection + diffusion + source + divergence;
          }
        }
      }

      /*
       * J_e_d_rs
       */
      for (b = 0; b < dim; b++) {
        var = SOLID_DISPLACEMENT1 + b;
        if (pd->v[pg->imtrx][var]) {
          pvar = upd->vp[pg->imtrx][var];
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];

            dh3dmesh_bj = fv->dh3dq[b] * bf[var]->phi[j];

            d_det_J_dmeshbj = bf[eqn]->d_det_J_dm[b][j];

            if (supg != 0.) {
              h_elem_deriv = 0.;
              h_elem_inv_deriv = 0.;
              for (qq = 0; qq < dim; qq++) {
                if (pg_data->hhv[qq][b] != 0.) {
                  h_elem_deriv -= vcent[qq] * vcent[qq] * pg_data->dhv_dxnode[qq][j] *
                                  pg_data->hhv[qq][b] * h_elem_inv / 4. / hsquared[qq] /
                                  hsquared[qq];
                }
              }
              if (h_elem != 0.)
                h_elem_inv_deriv = -h_elem_deriv / h_elem / h_elem;
              // h_elem_inv_deriv = 0.; /* PRS: NOT SURE WHY THIS IS NOT RIGHT, SO SET TO ZERO */
            }

            mass = 0.;

            advection = 0.;

            if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
              /*
               * one parts:
               *	Here, by "v" we mean "vconv" which is
               *    vconv = v_real_solid - vmesh,   vmesh ~ xdot ~ dd/dt
               *
               *	(a)	Int( Cp * dvconv/dreal_solid . grad(T) h3 detJ )
               *
               */

              advection_a = 0.;
              if (pd->e[pg->imtrx][eqn] & T_MASS) {
                for (p = 0; p < dim; p++) {
                  advection_a += d_vconv->X[p][b][j] * grad_T[p];
                }
                advection_a *= -wt_func * rho * Cp * h3 * det_J * wt;
              }

              advection = advection_a;

              advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
            }

            diffusion = 0.;

            source = 0.;

            lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += mass + advection + diffusion + source;
          }
        }
      }
      /*
       * J_e_d
       */
      for (b = 0; b < dim; b++) {
        var = MESH_DISPLACEMENT1 + b;
        if (pd->v[pg->imtrx][var]) {
          pvar = upd->vp[pg->imtrx][var];
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];

            dh3dmesh_bj = fv->dh3dq[b] * bf[var]->phi[j];

            d_det_J_dmeshbj = bf[eqn]->d_det_J_dm[b][j];

            if (supg != 0.) {
              h_elem_deriv = 0.;
              h_elem_inv_deriv = 0.;
              for (qq = 0; qq < dim; qq++) {
                if (pg_data->hhv[qq][b] != 0.) {
                  h_elem_deriv -= vcent[qq] * vcent[qq] * pg_data->dhv_dxnode[qq][j] *
                                  pg_data->hhv[qq][b] * h_elem_inv / 4. / hsquared[qq] /
                                  hsquared[qq];
                }
              }
              if (h_elem != 0.)
                h_elem_inv_deriv = -h_elem_deriv / h_elem / h_elem;
              // h_elem_inv_deriv = 0.; /* PRS: NOT SURE WHY THIS IS NOT RIGHT, SO SET TO ZERO */
            }

            mass = 0.;
            if (pd->TimeIntegration != STEADY) {
              if (pd->e[pg->imtrx][eqn] & T_MASS) {
                mass = T_dot;
                mass *= -phi_i * rho *
                        (Cp * h3 * d_det_J_dmeshbj + Cp * dh3dmesh_bj * det_J +
                         d_Cp->X[b][j] * h3 * det_J) *
                        wt;
                mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
              }
            }

            advection = 0.;

            if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
              /*
               * Four parts:
               *	Here, by "v" we mean "vconv" which is
               *    vconv = v - vmesh,   vmesh ~ xdot ~ dd/dt
               *
               *	d/dmesh [ Int(... vconv.grad(T) h3 detJ ) ]
               *
               *	(a)	Int( Cp * vconv.d(grad(T))/dmesh h3 detJ )
               *	(b)	Int( Cp * vconv.grad(T) h3 ddetJ/dmesh )
               *	(c)	Int( Cp * dvconv/dmesh . grad(T) h3 detJ )
               *	(d)	Int( Cp * vconv.grad(T) dh3/dmesh detJ )
               *	(e)	Int( dCp/dmesh * vconv.grad(T) h3 detJ )
               *
               */

              advection_a = 0.;
              for (p = 0; p < dim; p++) {
                advection_a += vconv[p] * fv->d_grad_T_dmesh[p][b][j];
              }
              advection_a *= -wt_func * rho * Cp * h3 * det_J * wt;

              advection_b = 0.;
              for (p = 0; p < dim; p++) {
                advection_b += vconv[p] * grad_T[p];
              }
              advection_b *= -wt_func * rho * Cp * h3 * d_det_J_dmeshbj * wt;

              advection_c = 0.;
              if (pd->TimeIntegration != STEADY) {
                if (pd->e[pg->imtrx][eqn] & T_MASS) {
                  for (p = 0; p < dim; p++) {
                    advection_c += d_vconv->X[p][b][j] * grad_T[p];
                  }
                  advection_c *= -wt_func * rho * Cp * h3 * det_J * wt;
                }
              }

              advection_d = 0.;
              for (p = 0; p < dim; p++) {
                advection_d += vconv[p] * grad_T[p];
              }
              advection_d *= -wt_func * rho * Cp * dh3dmesh_bj * det_J * wt;

              advection_e = 0.;
              for (p = 0; p < dim; p++) {
                advection_e += vconv[p] * grad_T[p];
              }

              advection_e *= -wt_func * rho * d_Cp->X[b][j] * h3 * det_J * wt;

              advection_f = 0.;
              if (supg != 0.) {
                d_wt_func = 0.;
                for (p = 0; p < dim; p++) {
                  d_wt_func +=
                      supg * (h_elem_inv * fv->v[p] * bf[eqn]->d_grad_phi_dmesh[i][p][b][j] +
                              h_elem_inv_deriv * fv->v[p] * grad_phi_i[p]);

                  advection_f += vconv[p] * grad_T[p];
                }
                advection_f *= -d_wt_func * h3 * det_J * wt;
              }

              advection =
                  advection_a + advection_b + advection_c + advection_d + advection_e + advection_f;

              advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
            }

            /*
             * multiple parts:
             * 	diff_a = Int(...d(grad_phi_i)/dmesh.q h3 |Jv|)
             *	diff_b = Int(...grad_phi_i.d(q)/dmesh h3 |Jv|)
             *	diff_c = Int(...grad_phi_i.q h3 d(|Jv|)/dmesh)
             *	diff_d = Int(...grad_phi_i.q dh3/dmesh |Jv|  )
             */
            diffusion = 0.;
            if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {
              diff_a = 0.;
              for (p = 0; p < dim; p++) {
                dgrad_phi_i_dmesh[p] = bf[eqn]->d_grad_phi_dmesh[i][p][b][j];

                diff_a += dgrad_phi_i_dmesh[p] * q[p];
              }
              diff_a *= det_J * h3 * wt;

              diff_b = 0.;
              for (p = 0; p < VIM; p++) {
                diff_b += d_q->X[p][b][j] * grad_phi_i[p];
              }
              diff_b *= det_J * h3 * wt;

              diff_c = 0.;
              for (p = 0; p < dim; p++) {
                diff_c += grad_phi_i[p] * q[p];
              }
              diff_c *= d_det_J_dmeshbj * h3 * wt;

              diff_d = 0.;
              for (p = 0; p < dim; p++) {
                diff_d += grad_phi_i[p] * q[p];
              }
              diff_d *= det_J * dh3dmesh_bj * wt;

              diffusion = diff_a + diff_b + diff_c + diff_d;

              diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
            }

            source = 0.;

            if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
              source =
                  phi_i *
                  (h * d_det_J_dmeshbj * h3 + h * det_J * dh3dmesh_bj + d_h->X[b][j] * det_J * h3) *
                  wt;

              source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
            }

            lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += mass + advection + diffusion + source;
          }
        }
      }

      /*
       * J_e_c
       */
      var = MASS_FRACTION;
      if (pd->e[pg->imtrx][eqn] && pd->v[pg->imtrx][var]) {
        for (w = 0; w < pd->Num_Species_Eqn; w++) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];

            mass = 0.;
            if (pd->TimeIntegration != STEADY) {
              if (pd->e[pg->imtrx][eqn] & T_MASS) {
                mass = T_dot * (rho * d_Cp->C[w][j] + d_rho->C[w][j] * Cp);
                mass *= -phi_i * det_J * wt;
                mass *= h3;
                mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
              }
            }

            advection = 0.;
            if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
              for (p = 0; p < dim; p++) {
                advection += rho * Cp * d_vconv->C[p][w][j] * grad_T[p];
                advection += (rho * d_Cp->C[w][j] + d_rho->C[w][j] * Cp) * vconv[p] * grad_T[p];
              }
              advection *= -wt_func * h3 * det_J * wt;
              advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
            }

            diffusion = 0.;
            if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {
              for (p = 0; p < dim; p++) {
                diffusion += grad_phi_i[p] * d_q->C[p][w][j];
              }
              diffusion *= det_J * wt;
              diffusion *= h3;
              diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
            }

            source = 0.;
            if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
              source += phi_i * d_h->C[w][j] * det_J * wt;
              source *= h3;
              source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
            }

            lec->J[LEC_J_INDEX(peqn, MAX_PROB_VAR + w, i, j)] +=
                advection + mass + diffusion + source;
          }
        }
      }

      /*
       * J_e_F
       */
      var = FILL;
      if (pd->e[pg->imtrx][eqn] && pd->v[pg->imtrx][var]) {
        pvar = upd->vp[pg->imtrx][var];
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          phi_j = bf[var]->phi[j];

          mass = 0.;
          if (pd->TimeIntegration != STEADY) {
            if (pd->e[pg->imtrx][eqn] & T_MASS) {
              mass = T_dot * (rho * d_Cp->F[j] + d_rho->F[j] * Cp);
              mass *= -phi_i * det_J * wt;
              mass *= h3;
              mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
            }
          }

          advection = 0.;
          if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
            for (p = 0; p < dim; p++) {
              advection += (rho * d_Cp->F[j] + d_rho->F[j] * Cp) * vconv[p] * grad_T[p];
            }
            advection *= -wt_func * h3 * det_J * wt;
            advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
          }

          diffusion = 0.;
          if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {
            for (p = 0; p < dim; p++) {
              diffusion += grad_phi_i[p] * d_q->F[p][j];
            }
            diffusion *= det_J * wt;
            diffusion *= h3;
            diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
          }

          source = 0.;
          if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
            source += phi_i * d_h->F[j] * det_J * wt;
            source *= h3;
            source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
          }

          lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += advection + mass + diffusion + source;
        }
      }

      /**    add ve stress terms  **/
      /*
       * J_e_S
       */
      for (mode = 0; mode < vn->modes; mode++) {
        for (a = 0; a < VIM; a++) {
          for (b = 0; b < VIM; b++) {
            var = v_s[mode][a][b];
            if (pd->v[pg->imtrx][var]) {
              pvar = upd->vp[pg->imtrx][var];
              for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                source = 0.;
                if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
                  source += phi_i * d_h->S[mode][a][b][j] * det_J * wt;
                  source *= h3;
                  source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
                }
                lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += source;
              }
            }
          }
        }
      }
      /*
       * J_e_apr and api
       */
      var = ACOUS_PREAL;

      if (pd->v[pg->imtrx][var]) {
        pvar = upd->vp[pg->imtrx][var];
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          mass = 0.;
          advection = 0.;
          diffusion = 0.;
          source = 0.;
          if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
            source += phi_i * d_h->APR[j] * det_J * wt;
            source *= h3;
            source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
          }

          lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += mass + advection + diffusion + source;
        }
      }
      var = ACOUS_PIMAG;

      if (pd->v[pg->imtrx][var]) {
        pvar = upd->vp[pg->imtrx][var];
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          mass = 0.;
          advection = 0.;
          diffusion = 0.;
          source = 0.;
          if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
            source += phi_i * d_h->API[j] * det_J * wt;
            source *= h3;
            source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
          }

          lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += mass + advection + diffusion + source;
        }
      }

      var = LIGHT_INTP;
      if (pd->v[pg->imtrx][var]) {
        pvar = upd->vp[pg->imtrx][var];
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          source = 0.;
          if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
            source += phi_i * d_h->INT[j] * det_J * wt;
            source *= h3;
            source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
          }

          lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += source;
        }
      }

      var = LIGHT_INTM;
      if (pd->v[pg->imtrx][var]) {
        pvar = upd->vp[pg->imtrx][LIGHT_INTM];
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          source = 0.;
          if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
            source += phi_i * d_h->INT[j] * det_J * wt;
            source *= h3;
            source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
          }
          lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += source;
        }
      }

      var = LIGHT_INTD;
      if (pd->v[pg->imtrx][var]) {
        pvar = upd->vp[pg->imtrx][var];
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          source = 0.;
          if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
            source += phi_i * d_h->INT[j] * det_J * wt;
            source *= h3;
            source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
          }
          lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += source;
        }
      }
    }
  }

  return (status);
} /********************************************************************************/

double conductivity(CONDUCTIVITY_DEPENDENCE_STRUCT *d_k, dbl time)

/**************************************************************************
 *
 * conductivity
 *
 *   Calculate the thermal conductivity and its derivatives wrt to nodal dofs
 *   at the local gauss point
 *
 * Output
 * -----
 *    dbl d_k->T[j]    -> derivative of conductivity wrt the jth
 *                          Temperature unknown in an element
 *    dbl d_k->C[k][j] -> derivative of conductivity wrt the jth
 *                          species unknown of ktype, k, in the element.
 *    dbl d_k->X[a][j] -> derivative of conductivity wrt the jth
 *                          mesh displacement in the ath direction.
 *    dbl d_k->F[j]    -> derivative of conductivity wrt the jth
 *                          FILL unknown in an element
 *
 *  Return
 * --------
 *    Actual value of the conductivity
 ***************************************************************************/

{
  int dim = pd->Num_Dim;
  int var, i, j, a, w, var_offset;
  double k, Tref, tmp, temp;
  struct Level_Set_Data *ls_old;

  int i_therm_cond;

  if (pd->gv[TEMPERATURE]) {
    temp = fv->T;
  } else {
    temp = upd->Process_Temperature;
  }

  k = 0.;

  if (d_k != NULL) {
    memset(d_k->T, 0, sizeof(double) * MDE);
    memset(d_k->X, 0, sizeof(double) * DIM * MDE);
    memset(d_k->C, 0, sizeof(double) * MAX_CONC * MDE);
    memset(d_k->F, 0, sizeof(double) * MDE);
  }

  if (mp->ConductivityModel == USER) {

    usr_thermal_conductivity(mp->u_thermal_conductivity, time);

    k = mp->thermal_conductivity;

    var = TEMPERATURE;
    if (pd->v[pg->imtrx][TEMPERATURE] && d_k != NULL) {

      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_k->T[j] = mp->d_thermal_conductivity[var] * bf[var]->phi[j];
      }
    }

    if (pd->v[pg->imtrx][MESH_DISPLACEMENT1] && d_k != NULL) {
      for (a = 0; a < dim; a++) {
        var = MESH_DISPLACEMENT1 + a;
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_k->X[a][j] = mp->d_thermal_conductivity[var] * bf[var]->phi[j];
        }
      }
    }

    if (pd->v[pg->imtrx][MASS_FRACTION] && d_k != NULL) {
      for (w = 0; w < pd->Num_Species_Eqn; w++) {
        var = MASS_FRACTION;
        var_offset = MAX_VARIABLE_TYPES + w;
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_k->C[w][j] = mp->d_thermal_conductivity[var_offset] * bf[var]->phi[j];
        }
      }
    }
  } else if (mp->ConductivityModel == CONSTANT) {
    k = mp->thermal_conductivity;
  } else if (mp->ConductivityModel == THERMAL_HEAT) {
    Tref = mp->u_thermal_conductivity[4];
    tmp = temp - Tref;
    mp->thermal_conductivity =
        mp->u_thermal_conductivity[0] +
        tmp * (mp->u_thermal_conductivity[1] +
               tmp * (mp->u_thermal_conductivity[2] + tmp * mp->u_thermal_conductivity[3]));
    k = mp->thermal_conductivity;
    var = TEMPERATURE;
    if (pd->v[pg->imtrx][var] && d_k != NULL) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_k->T[j] =
            (mp->u_thermal_conductivity[1] + tmp * (2. * mp->u_thermal_conductivity[2] +
                                                    tmp * 3. * mp->u_thermal_conductivity[3])) *
            bf[var]->phi[j];
      }
    }
  } else if (mp->ConductivityModel == FOAM_PBE) {

    k = foam_pbe_conductivity(d_k, time);
  } else if (mp->ConductivityModel == FOAM_PMDI_10) {
    double rho;
    DENSITY_DEPENDENCE_STRUCT d_rho_struct;
    DENSITY_DEPENDENCE_STRUCT *d_rho = &d_rho_struct;

    rho = density(d_rho, time);
    if (mp->len_u_thermal_conductivity < 2) {
      GOMA_EH(GOMA_ERROR, "Expected at least 2 constants for thermal conductivity FOAM_PMDI_10");
      return 0;
    }
    if (mp->DensityModel != DENSITY_FOAM_PMDI_10) {
      GOMA_EH(GOMA_ERROR, "FOAM_PMDI_10 Thermal conductivity requires FOAM_PMDI_10 density");
      return 0;
    }

    double k_liq = mp->u_thermal_conductivity[0];
    double k_gas = mp->u_thermal_conductivity[1];

    double rho_liq = mp->u_density[1];

    mp->thermal_conductivity = (2.0 / 3.0) * (rho / rho_liq) * k_liq + (1 - rho / rho_liq) * k_gas;

    k = mp->thermal_conductivity;

    var = TEMPERATURE;
    if (pd->v[pg->imtrx][var] && d_k != NULL) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_k->T[j] = (2.0 / 3.0) * (d_rho->T[j] / rho_liq) * k_liq - (d_rho->T[j] / rho_liq) * k_gas;
      }
    }

    int w;
    var = MASS_FRACTION;
    if (pd->v[pg->imtrx][var] && d_k != NULL) {
      for (w = 0; w < pd->Num_Species_Eqn; w++) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_k->C[w][j] =
              (2.0 / 3.0) * (d_rho->C[w][j] / rho_liq) * k_liq - (d_rho->C[w][j] / rho_liq) * k_gas;
        }
      }
    }
  } else if (mp->ConductivityModel == TABLE) {
    struct Data_Table *table_local;
    table_local = MP_Tables[mp->thermal_conductivity_tableid];
    apply_table_mp(&mp->thermal_conductivity, table_local);
    k = mp->thermal_conductivity;

    if (d_k != NULL) {
      for (i = 0; i < table_local->columns - 1; i++) {
        var = table_local->t_index[i];
        /* currently only set up to vary w.r.t. temperature */
        switch (var) {
        case TEMPERATURE:
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            d_k->T[j] = table_local->slope[i] * bf[var]->phi[j];
          }
          break;
        default:
          GOMA_EH(GOMA_ERROR, "Variable function not yet implemented in material property table");
        }
      }
    }
  } else if (mp->ConductivityModel == LEVEL_SET) {
    ls_transport_property(mp->u_thermal_conductivity[0], mp->u_thermal_conductivity[1],
                          mp->u_thermal_conductivity[2], &mp->thermal_conductivity,
                          &mp->d_thermal_conductivity[FILL]);

    k = mp->thermal_conductivity;

    var = FILL;
    if (pd->v[pg->imtrx][var] && d_k != NULL) {

      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_k->F[j] = mp->d_thermal_conductivity[var] * bf[var]->phi[j];
      }
    }

  } else if (mp->ConductivityModel == EXTERNAL_FIELD) {
    i_therm_cond = mp->thermal_cond_external_field;
    GOMA_EH(i_therm_cond, "Thermal cond. external field not found!");
    k = fv->external_field[i_therm_cond] * mp->u_thermal_conductivity[0];
  } else {
    GOMA_EH(GOMA_ERROR, "Unrecognized thermal conductivity model");
  }

  if (ls != NULL && mp->ConductivityModel != LEVEL_SET && mp->ConductivityModel != LS_QUADRATIC &&
      mp->mp2nd != NULL &&
      mp->mp2nd->ThermalConductivityModel ==
          CONSTANT) /* Only Newtonian constitutive equation allowed for 2nd phase */
  {
    /* kludge for solidification tracking with phase function 0 */
    if (pfd != NULL && pd->e[pg->imtrx][R_EXT_VELOCITY]) {
      ls_old = ls;
      ls = pfd->ls[0];
      k = ls_modulate_thermalconductivity(k, mp->mp2nd->thermalconductivity_phase[0],
                                          ls->Length_Scale,
                                          (double)mp->mp2nd->thermalconductivitymask[0],
                                          (double)mp->mp2nd->thermalconductivitymask[1], d_k);
      ls = ls_old;
    }
    k = ls_modulate_thermalconductivity(k, mp->mp2nd->thermalconductivity, ls->Length_Scale,
                                        (double)mp->mp2nd->thermalconductivitymask[0],
                                        (double)mp->mp2nd->thermalconductivitymask[1], d_k);
  }

  return (k);
}
double heat_capacity(HEAT_CAPACITY_DEPENDENCE_STRUCT *d_Cp, dbl time)
/**************************************************************************
 *
 * heat capacity
 *
 *   Calculate the thermal heat capacity and its derivatives wrt to nodal dofs
 *   at the local gauss point
 *
 * Output
 * -----
 *    dbl d_Cp->T[j]    -> derivative of heat capacity wrt the jth
 *                          Temperature unknown in an element
 *    dbl d_Cp->V[j]    -> derivative of heat capacity wrt the jth
 *                          volage unknown in an element
 *    dbl d_Cp->v[a][j] -> derivative of heat capacity wrt the jth
 *                          velocity  in the ath direction.
 *    dbl d_Cp->C[k][j] -> derivative of heat capacity wrt the jth
 *                          species unknown of ktype, k, in the element.
 *    dbl d_Cp->X[a][j] -> derivative of heat capacity wrt the jth
 *                          mesh displacement in the ath direction.
 *    dbl d_Cp->F[j]    -> derivative of heat capacity wrt the jth
 *                          FILL unknown in an element
 *
 *  Return
 * --------
 *    Actual value of the heat capacity
 ***************************************************************************/

{
  int dim = pd->Num_Dim;
  int var, i, j, a, w, var_offset;
  double Cp, T_offset, temp;
  struct Level_Set_Data *ls_old;

  if (pd->gv[TEMPERATURE]) {
    temp = fv->T;
  } else {
    temp = upd->Process_Temperature;
  }

  Cp = 0.;

  if (d_Cp != NULL) {
    memset(d_Cp->T, 0, sizeof(double) * MDE);
    memset(d_Cp->V, 0, sizeof(double) * MDE);
    memset(d_Cp->X, 0, sizeof(double) * DIM * MDE);
    memset(d_Cp->v, 0, sizeof(double) * DIM * MDE);
    memset(d_Cp->C, 0, sizeof(double) * MAX_CONC * MDE);
    memset(d_Cp->F, 0, sizeof(double) * MDE);
  }

  if (mp->HeatCapacityModel == USER) {
    usr_heat_capacity(mp->u_heat_capacity, time);
    Cp = mp->heat_capacity;

    var = TEMPERATURE;
    if (d_Cp != NULL && pd->e[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_Cp->T[j] = mp->d_heat_capacity[var] * bf[var]->phi[j];
      }
    }

    if (d_Cp != NULL && pd->v[pg->imtrx][MESH_DISPLACEMENT1]) {
      for (a = 0; a < dim; a++) {
        var = MESH_DISPLACEMENT1 + a;
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_Cp->X[a][j] = mp->d_heat_capacity[var] * bf[var]->phi[j];
        }
      }
    }

    if (d_Cp != NULL && pd->v[pg->imtrx][FILL]) {
      var = FILL;

      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_Cp->F[j] = mp->d_heat_capacity[var] * bf[var]->phi[j];
      }
    }

    if (d_Cp != NULL && pd->v[pg->imtrx][MASS_FRACTION]) {
      for (w = 0; w < pd->Num_Species_Eqn; w++) {
        var = MASS_FRACTION;
        var_offset = MAX_VARIABLE_TYPES + w;
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_Cp->C[w][j] = mp->d_heat_capacity[var_offset] * bf[var]->phi[j];
        }
      }
    }

  } else if (mp->HeatCapacityModel == CONSTANT) {
    Cp = mp->heat_capacity;
  } else if (mp->HeatCapacityModel == THERMAL_HEAT) {
    T_offset = mp->u_heat_capacity[4];
    mp->heat_capacity = mp->u_heat_capacity[0] + mp->u_heat_capacity[1] * (temp + T_offset) +
                        mp->u_heat_capacity[2] * SQUARE(temp + T_offset) +
                        mp->u_heat_capacity[3] / SQUARE(temp + T_offset);
    Cp = mp->heat_capacity;
    var = TEMPERATURE;
    if (d_Cp != NULL && pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_Cp->T[j] = (mp->u_heat_capacity[1] + 2. * (temp + T_offset) * mp->u_heat_capacity[2] -
                      2. * mp->u_heat_capacity[3] / (SQUARE(temp + T_offset) * (temp + T_offset))) *
                     bf[var]->phi[j];
      }
    }
  } else if (mp->HeatCapacityModel == ENTHALPY) {
    Cp = enthalpy_heat_capacity_model(d_Cp);
  } else if (mp->HeatCapacityModel == FOAM_PMDI_10) {
    Cp = foam_pmdi_10_heat_cap(d_Cp, time);
  } else if (mp->HeatCapacityModel == LEVEL_SET) {
    ls_transport_property(mp->u_heat_capacity[0], mp->u_heat_capacity[1], mp->u_heat_capacity[2],
                          &mp->heat_capacity, &mp->d_heat_capacity[FILL]);

    Cp = mp->heat_capacity;

    var = FILL;
    if (d_Cp != NULL && pd->v[pg->imtrx][var]) {

      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_Cp->F[j] = mp->d_heat_capacity[var] * bf[var]->phi[j];
      }
    }

  } else if (mp->HeatCapacityModel == TABLE) {
    struct Data_Table *table_local;
    table_local = MP_Tables[mp->heat_capacity_tableid];
    apply_table_mp(&mp->heat_capacity, table_local);
    Cp = mp->heat_capacity;

    if (d_Cp != NULL) {
      for (i = 0; i < table_local->columns - 1; i++) {
        var = table_local->t_index[i];
        /* currently only set up to vary w.r.t. temperature */
        switch (var) {
        case TEMPERATURE:
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            d_Cp->T[j] = table_local->slope[i] * bf[var]->phi[j];
          }
          break;
        default:
          GOMA_EH(GOMA_ERROR, "Variable function not yet implemented in material property table");
        }
      }
    }
  } else {
    GOMA_EH(GOMA_ERROR, "Unrecognized heat capacity model");
  }

  if (ls != NULL && mp->HeatCapacityModel != LEVEL_SET && mp->HeatCapacityModel != LS_QUADRATIC &&
      mp->mp2nd != NULL &&
      mp->mp2nd->HeatCapacityModel ==
          CONSTANT) /* Only Newtonian constitutive equation allowed for 2nd phase */
  {
    /* kludge for solidification tracking with phase function 0 */
    if (pfd != NULL && pd->e[pg->imtrx][R_EXT_VELOCITY]) {
      ls_old = ls;
      ls = pfd->ls[0];
      Cp = ls_modulate_heatcapacity(Cp, mp->mp2nd->heatcapacity_phase[0], ls->Length_Scale,
                                    (double)mp->mp2nd->heatcapacitymask[0],
                                    (double)mp->mp2nd->heatcapacitymask[1], d_Cp);
      ls = ls_old;
    }
    Cp = ls_modulate_heatcapacity(Cp, mp->mp2nd->heatcapacity, ls->Length_Scale,
                                  (double)mp->mp2nd->heatcapacitymask[0],
                                  (double)mp->mp2nd->heatcapacitymask[1], d_Cp);
  }

  return (Cp);
}
double ls_modulate_thermalconductivity(double k1,
                                       double k2,
                                       double width,
                                       double pm_minus,
                                       double pm_plus,
                                       CONDUCTIVITY_DEPENDENCE_STRUCT *d_k) {
  double factor;
  int i, a, w, var;

  if (d_k == NULL) {
    k1 = ls_modulate_property(k1, k2, width, pm_minus, pm_plus, NULL, &factor);
    return (k1);
  }

  k1 = ls_modulate_property(k1, k2, width, pm_minus, pm_plus, d_k->F, &factor);

  if (pd->v[pg->imtrx][var = TEMPERATURE]) {
    for (i = 0; i < ei[pg->imtrx]->dof[var]; i++) {
      d_k->T[i] *= factor;
    }
  }

  if (pd->v[pg->imtrx][var = MASS_FRACTION]) {
    for (w = 0; w < pd->Num_Species_Eqn; w++) {
      for (i = 0; i < ei[pg->imtrx]->dof[var]; i++) {
        d_k->C[w][i] *= factor;
      }
    }
  }

  if (pd->v[pg->imtrx][var = MESH_DISPLACEMENT1]) {
    for (a = 0; a < pd->Num_Dim; a++) {
      for (i = 0; i < ei[pg->imtrx]->dof[var]; i++) {
        d_k->X[a][i] *= factor;
      }
    }
  }
  return (k1);
}
double ls_modulate_heatcapacity(double Cp1,
                                double Cp2,
                                double width,
                                double pm_minus,
                                double pm_plus,
                                HEAT_CAPACITY_DEPENDENCE_STRUCT *d_Cp) {
  double factor;
  int i, a, w, var;

  if (d_Cp == NULL) {
    Cp1 = ls_modulate_property(Cp1, Cp2, width, pm_minus, pm_plus, NULL, &factor);
    return (Cp1);
  }

  Cp1 = ls_modulate_property(Cp1, Cp2, width, pm_minus, pm_plus, d_Cp->F, &factor);

  if (pd->v[pg->imtrx][var = TEMPERATURE]) {
    for (i = 0; i < ei[pg->imtrx]->dof[var]; i++) {
      d_Cp->T[i] *= factor;
    }
  }

  if (pd->v[pg->imtrx][var = MASS_FRACTION]) {
    for (w = 0; w < pd->Num_Species_Eqn; w++) {
      for (i = 0; i < ei[pg->imtrx]->dof[var]; i++) {
        d_Cp->C[w][i] *= factor;
      }
    }
  }

  if (pd->v[pg->imtrx][var = VELOCITY1]) {
    for (a = 0; a < pd->Num_Dim; a++) {
      for (i = 0; i < ei[pg->imtrx]->dof[var]; i++) {
        d_Cp->v[a][i] *= factor;
      }
    }
  }

  if (pd->v[pg->imtrx][var = MESH_DISPLACEMENT1]) {
    for (a = 0; a < pd->Num_Dim; a++) {
      for (i = 0; i < ei[pg->imtrx]->dof[var]; i++) {
        d_Cp->X[a][i] *= factor;
      }
    }
  }

  if (pd->v[pg->imtrx][var = VOLTAGE]) {
    for (i = 0; i < ei[pg->imtrx]->dof[var]; i++) {
      d_Cp->V[i] *= factor;
    }
  }

  return (Cp1);
}
void heat_flux(double q[DIM], HEAT_FLUX_DEPENDENCE_STRUCT *d_q, double time) {
  dbl k; /* Thermal conductivity. */
  CONDUCTIVITY_DEPENDENCE_STRUCT d_k_struct;
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_k = &d_k_struct;

  dbl dq_gradT[DIM][DIM]; /*  Heat flux sensitivities  */
  dbl dq_dX[DIM][DIM];    /*  Heat flux sensitivities wrt position  */
  dbl grad_T[DIM];        /* Temperature gradient. */

  int b, j = -1, p, a, w;
  int var;

  if (d_q == NULL)
    d_k = NULL;

  k = conductivity(d_k, time);

  for (p = 0; p < VIM; p++) {
    grad_T[p] = fv->grad_T[p];
  }

  if (cr->HeatFluxModel == CR_HF_FOURIER_0) {
    for (p = 0; p < VIM; p++) {
      q[p] = -k * grad_T[p];
    }

    var = TEMPERATURE;
    if (d_q != NULL && pd->v[pg->imtrx][var]) {
      for (p = 0; p < VIM; p++) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_q->T[p][j] = -k * bf[var]->grad_phi[j][p] - d_k->T[j] * grad_T[p];
        }
      }
    }

    var = MASS_FRACTION;
    if (d_q != NULL && pd->v[pg->imtrx][var]) {
      for (w = 0; w < pd->Num_Species_Eqn; w++) {
        for (p = 0; p < VIM; p++) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            d_q->C[p][w][j] = -d_k->C[w][j] * grad_T[p];
          }
        }
      }
    }

    var = FILL;
    if (d_q != NULL && pd->v[pg->imtrx][var]) {
      for (p = 0; p < VIM; p++) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_q->F[p][j] = -d_k->F[j] * grad_T[p];
        }
      }
    }

    if (d_q != NULL && pd->v[pg->imtrx][MESH_DISPLACEMENT1]) {
      for (p = 0; p < VIM; p++) {
        for (b = 0; b < WIM; b++) {
          var = MESH_DISPLACEMENT1 + b;
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            d_q->X[p][b][j] = -k * fv->d_grad_T_dmesh[p][b][j] - d_k->X[b][j] * grad_T[p];
          }
        }
      }
    }
  } else if (cr->HeatFluxModel == CR_HF_USER) {

    usr_heat_flux(grad_T, q, dq_gradT, dq_dX, time);

    var = TEMPERATURE;
    if (d_q != NULL && pd->v[pg->imtrx][var]) {
      for (p = 0; p < VIM; p++) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_q->T[p][j] = 0.0;
          for (a = 0; a < VIM; a++) {
            d_q->T[p][j] += dq_gradT[p][a] * bf[var]->grad_phi[j][a];
          }
        }
      }
    }
    for (b = 0; b < WIM; b++) {
      var = MESH_DISPLACEMENT1 + b;
      if (d_q != NULL && pd->v[pg->imtrx][var]) {
        for (p = 0; p < VIM; p++) {
          d_q->X[p][b][j] = dq_dX[p][b] * bf[var]->phi[j];
          for (a = 0; a < VIM; a++) {
            d_q->X[p][b][j] += dq_gradT[p][a] * fv->d_grad_T_dmesh[a][b][j];
          }
        }
      }
    }
  } else {
    GOMA_EH(-1, "Unimplemented thermal constitutive relation.");
  }
}
double heat_source(HEAT_SOURCE_DEPENDENCE_STRUCT *d_h,
                   double time, /* current time */
                   double tt,   /* parameter to vary time integration from
                                   explicit (tt = 1) to implicit (tt = 0) */
                   double dt)   /* current time step size */
{
  double h = 0.;
  double h_acous = 0.;
  int j, w, a, var_offset;
  int var;

  const double F = 96487.0;        /* Faraday's constant in units of C/euiv.; KSC: 2/17/99 */
  const double R = 8.314;          /* Universal gas constant in units of J/mole K */
  dbl TT;                          /* dummy variable */
  dbl ai0, ai0_anode, ai0_cathode; /* exchange current density; KSC: 2/17/99 */
  int mn;

  struct Level_Set_Data *ls_old = ls;

  if (MAX_CONC < 3) {
    GOMA_EH(GOMA_ERROR, "heat_source expects MAX_CONC >= 3");
    return 0;
  }

  /*
   * Unpack variables from structures for local convenience...
   */

  /* initialize Heat Source sensitivities */
  if (d_h != NULL) {
    memset(d_h->T, 0, sizeof(double) * MDE);
    memset(d_h->F, 0, sizeof(double) * MDE);
    memset(d_h->V, 0, sizeof(double) * MDE);
    memset(d_h->v, 0, sizeof(double) * DIM * MDE);
    memset(d_h->X, 0, sizeof(double) * DIM * MDE);
    memset(d_h->C, 0, sizeof(double) * MAX_CONC * MDE);
    memset(d_h->S, 0, sizeof(double) * MAX_MODES * DIM * DIM * MDE);
    memset(d_h->APR, 0, sizeof(double) * MDE);
    memset(d_h->API, 0, sizeof(double) * MDE);
    memset(d_h->INT, 0, sizeof(double) * MDE);
  }

  if (mp->HeatSourceModel == USER) {
    usr_heat_source(mp->u_heat_source, time);
    h = mp->heat_source;

    var = TEMPERATURE;
    if (d_h != NULL && pd->e[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_h->T[j] = mp->d_heat_source[var] * bf[var]->phi[j];
      }
    }

    var = VOLTAGE;
    if (d_h != NULL && pd->e[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_h->V[j] = mp->d_heat_source[var] * bf[var]->phi[j];
      }
    }

    if (d_h != NULL && pd->v[pg->imtrx][VELOCITY1]) {
      for (a = 0; a < WIM; a++) {
        var = VELOCITY1 + a;
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_h->v[a][j] = mp->d_heat_source[var] * bf[var]->phi[j];
        }
      }
    }

    if (d_h != NULL && pd->v[pg->imtrx][MESH_DISPLACEMENT1]) {
      for (a = 0; a < WIM; a++) {
        var = MESH_DISPLACEMENT1 + a;
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_h->X[a][j] = mp->d_heat_source[var] * bf[var]->phi[j];
        }
      }
    }

    if (d_h != NULL && pd->v[pg->imtrx][MASS_FRACTION]) {
      for (w = 0; w < pd->Num_Species_Eqn; w++) {
        var = MASS_FRACTION;
        var_offset = MAX_VARIABLE_TYPES + w;
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_h->C[w][j] = mp->d_heat_source[var_offset] * bf[var]->phi[j];
        }
      }
    }

  } else if (mp->HeatSourceModel == PHOTO_CURING) {
    double intensity;
    double k_prop, k_inh = 0, free_rad, d_free_rad_dI;
    double *param, dhdC[MAX_CONC], dhdT, dhdI, Conc[MAX_CONC];
    int model_bit, num_mon, O2_spec = -1, rad_spec = -1, init_spec = 0;
    double k_propX = 1, k_propT = 0, k_propX_num = 0, k_propX_den = 0;
    double intensity_cgs = 2.998e+10 * 8.85e-12 / 200.0;
    double dbl_small = 1.0e-15, Xconv_denom = 0, sum_init = 0;
    double Xconv = 0.0, Xconv_init = 0.0, dXdC[MAX_CONC] = {0.0}, sum_mon = 0;

    param = mp->u_heat_source;
    model_bit = ((int)param[0]);
    h = 0;
    dhdT = 0;
    dhdI = 0;
    d_free_rad_dI = 0;
    for (a = 0; a < MAX_CONC; a++)
      dhdC[a] = 0.;

    intensity = 0.;
    if (pd->e[pg->imtrx][R_LIGHT_INTP]) {
      intensity += fv->poynt[0];
      if (pd->gv[R_LIGHT_INTM]) {
        intensity += fv->poynt[1];
      }
      if (pd->gv[R_LIGHT_INTD]) {
        intensity += fv->poynt[2];
      }
      intensity *= mp->u_species_source[init_spec][1];
      intensity = MAX(intensity, 0.0);
    } else if (pd->gv[R_ACOUS_PREAL]) {
      intensity =
          mp->u_species_source[init_spec][1] * intensity_cgs * (SQUARE(fv->apr) + SQUARE(fv->api));
    } else {
      GOMA_WH(-1, "No Intensity field found in PHOTO_CURING\n");
    }

    /* insure concentrations are positive  */
    for (j = 0; j < pd->Num_Species_Eqn; j++) {
      Conc[j] = MAX(dbl_small, fv->c[j]);
    }

    /**  heat source from momomer heat of reaction     **/
    num_mon = model_bit >> 2;

    if ((model_bit & 1) && (model_bit & 2)) {
      O2_spec = init_spec + num_mon + 2;
      rad_spec = O2_spec + 1;
      free_rad = Conc[rad_spec];
    } else if (model_bit & 1) {
      O2_spec = init_spec + num_mon + 2;
      k_inh = mp->u_species_source[O2_spec][1] *
              exp(-mp->u_species_source[O2_spec][2] *
                  (1. / fv->T - 1. / mp->u_species_source[O2_spec][3]));
      free_rad = sqrt(SQUARE(k_inh * Conc[O2_spec]) / 4. +
                      mp->u_species_source[init_spec + 1][2] * intensity * Conc[init_spec]) -
                 k_inh * Conc[O2_spec] / 2.;
      d_free_rad_dI += 0.5 * Conc[init_spec] * mp->u_species_source[init_spec + 1][2] /
                       sqrt(SQUARE(k_inh * Conc[O2_spec]) / 4. +
                            mp->u_species_source[init_spec + 1][2] * intensity * Conc[init_spec]);
    } else if (model_bit & 2) {
      rad_spec = init_spec + num_mon + 2;
      free_rad = Conc[rad_spec];
    } else {
      free_rad = sqrt(mp->u_species_source[init_spec + 1][2] * intensity * Conc[init_spec]);
      if (free_rad > 0) {
        d_free_rad_dI += 0.5 / free_rad * mp->u_species_source[init_spec + 1][2] * Conc[init_spec];
      }
    }

    switch (mp->Species_Var_Type) {
    case SPECIES_DENSITY:
      for (w = init_spec + 2; w < init_spec + 2 + num_mon; w++) {
        Xconv += fv->external_field[w] / mp->molecular_weight[w] / mp->specific_volume[w];
        sum_mon += Conc[w] / mp->molecular_weight[w];
        dXdC[w] = -1.0 / mp->molecular_weight[w];
        Xconv_init +=
            mp->u_reference_concn[w][1] / mp->molecular_weight[w] / mp->specific_volume[w];
        sum_init += mp->u_reference_concn[w][0] / mp->molecular_weight[w];
      }
      break;
    case SPECIES_CONCENTRATION:
      for (w = init_spec + 2; w < init_spec + 2 + num_mon; w++) {
        Xconv += fv->external_field[w] / mp->specific_volume[w];
        sum_mon += Conc[w];
        dXdC[w] = -1.0;
        Xconv_init += mp->u_reference_concn[w][1] / mp->specific_volume[w];
        sum_init += mp->u_reference_concn[w][0];
      }
      break;
    default:
      GOMA_EH(GOMA_ERROR, "invalid Species Type for PHOTO_CURING\n");
    }
    Xconv *= mp->specific_volume[pd->Num_Species_Eqn];
    Xconv_init *= mp->specific_volume[pd->Num_Species_Eqn];
    Xconv_denom = Xconv + sum_mon;
    Xconv /= Xconv_denom;
    Xconv = MAX(dbl_small, Xconv);
    Xconv = MIN(1.0 - dbl_small, Xconv);
    Xconv_init /= (Xconv_init + sum_init);
    for (w = init_spec + 2; w < init_spec + 2 + num_mon; w++) {
      dXdC[w] *= Xconv / Xconv_denom;
    }
    if (Xconv <= dbl_small || Xconv >= (1.0 - dbl_small)) {
      memset(dXdC, 0, sizeof(double) * MAX_CONC);
    }

    for (w = init_spec + 2; w < init_spec + num_mon + 2; w++) {
      k_prop = mp->u_species_source[w][1] *
               exp(-mp->u_species_source[w][2] * (1. / fv->T - 1. / mp->u_species_source[w][3]));
      k_propX_num = (1.0 - mp->u_species_source[w][4]) * (1. - Xconv) +
                    mp->u_species_source[w][4] * (1.0 - Xconv_init);
      k_propX_den = k_propX_num - (1.0 - mp->u_species_source[w][4]) * (1. - Xconv) *
                                      log((1. - Xconv) / (1.0 - Xconv_init));
      k_propX = SQUARE(k_propX_num) / k_propX_den;
      k_propT = k_prop * k_propX;

      h += k_propT * Conc[w] * free_rad * param[w] * mp->molecular_weight[w];
      dhdC[w] = dXdC[w] * (mp->u_species_source[w][4] - 1.0) * k_propX *
                (2. / k_propX_num + log((1. - Xconv) / (1.0 - Xconv_init)) / k_propX_den);
      dhdC[w] *= k_prop * Conc[w];
      dhdC[w] += k_propT;
      dhdC[w] *= free_rad * param[w] * mp->molecular_weight[w];

      dhdT += k_propX * k_prop * mp->u_species_source[w][2] / SQUARE(fv->T) * Conc[w] * free_rad *
              param[w] * mp->molecular_weight[w];
      dhdI += k_propT * d_free_rad_dI * param[w] * mp->molecular_weight[w];

      if (model_bit & 2) {
        dhdC[rad_spec] += k_propT * Conc[w] * param[w] * mp->molecular_weight[w];
      } else if (model_bit & 1) {
        dhdC[O2_spec] +=
            k_propT * Conc[w] * param[w] * mp->molecular_weight[w] *
            (SQUARE(k_inh / 2.) * Conc[O2_spec] /
                 sqrt(SQUARE(k_inh * Conc[O2_spec]) / 4. +
                      mp->u_species_source[init_spec + 1][2] * intensity * Conc[init_spec]) -
             k_inh / 2.);
        dhdC[init_spec] +=
            k_propT * Conc[w] * param[w] * mp->molecular_weight[w] * 0.5 *
            mp->u_species_source[init_spec + 1][2] * intensity /
            sqrt(SQUARE(k_inh * Conc[O2_spec]) / 4. +
                 mp->u_species_source[init_spec + 1][2] * intensity * Conc[init_spec]);

      } else {
        dhdC[init_spec] +=
            k_propT * Conc[w] * param[w] * mp->molecular_weight[w] * 0.5 *
            sqrt(mp->u_species_source[init_spec + 1][2] * intensity / Conc[init_spec]);
      }
    }

    /**  add heat generation from light absorption  **/

    h += param[1] * intensity * mp->light_absorption;
    dhdI = param[1] * mp->light_absorption;
    for (w = 0; w < pd->Num_Species_Eqn; w++) {
      dhdC[w] += param[1] * intensity * mp->d_light_absorption[MAX_VARIABLE_TYPES + w];
    }

    var = TEMPERATURE;
    if (d_h != NULL && pd->e[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_h->T[j] = dhdT * bf[var]->phi[j];
      }
    }
    var = LIGHT_INTP;
    if (d_h != NULL && pd->e[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_h->INT[j] = dhdI * bf[var]->phi[j] * mp->u_species_source[init_spec][1];
      }
    }
    var = ACOUS_PREAL;
    if (d_h != NULL && pd->e[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_h->APR[j] = dhdI * bf[var]->phi[j] * mp->u_species_source[init_spec][1] * intensity_cgs *
                      2.0 * fv->apr;
      }
    }
    var = ACOUS_PIMAG;
    if (d_h != NULL && pd->e[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_h->API[j] = dhdI * bf[var]->phi[j] * mp->u_species_source[init_spec][1] * intensity_cgs *
                      2.0 * fv->api;
      }
    }

    if (d_h != NULL && pd->v[pg->imtrx][MASS_FRACTION]) {
      for (w = 0; w < pd->Num_Species_Eqn; w++) {
        var = MASS_FRACTION;
        var_offset = MAX_VARIABLE_TYPES + w;
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_h->C[w][j] = dhdC[w] * bf[var]->phi[j];
        }
      }
    }

  } else if (mp->HeatSourceModel == DROP_EVAP) {
    struct Species_Conservation_Terms s_terms;
    int err, w1;
    zero_structure(&s_terms, sizeof(struct Species_Conservation_Terms), 1);

    err = get_continuous_species_terms(&s_terms, time, tt, dt, NULL);
    GOMA_EH(err, "problem in getting the species terms");

    for (w = 0; w < pd->Num_Species_Eqn; w++) {
      if (mp->SpeciesSourceModel[w] == DROP_EVAP) {
        h -= mp->latent_heat_vap[w] * s_terms.MassSource[w];

        var = TEMPERATURE;
        if (d_h != NULL && pd->e[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            d_h->T[j] -= mp->latent_heat_vap[w] * s_terms.d_MassSource_dT[w][j];
          }
        }
        if (d_h != NULL && pd->v[pg->imtrx][MASS_FRACTION]) {
          for (w1 = 0; w1 < pd->Num_Species_Eqn; w1++) {
            var = MASS_FRACTION;
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              d_h->C[w][j] -= mp->latent_heat_vap[w] * s_terms.d_MassSource_dc[w][w1][j];
            }
          }
        }
        var = RESTIME;
        if (d_h != NULL && pd->e[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            d_h->rst[j] -= mp->latent_heat_vap[w] * s_terms.d_MassSource_drst[w][j];
          }
        }
        var = PRESSURE;
        if (d_h != NULL && pd->e[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            d_h->P[j] -= mp->latent_heat_vap[w] * s_terms.d_MassSource_dP[w][j];
          }
        }
      }
    }
  } else if (mp->HeatSourceModel == CONSTANT) {
    h = mp->heat_source;
  } else if (mp->HeatSourceModel == ELECTRODE_KINETICS) /* added by KSC/GHE: 10/20/98 */
  {
    electrolyte_temperature(time, dt,
                            0); /* calculate electrolyte temperature at the present time */
    TT = mp->electrolyte_temperature;
    h = 0.0;
    for (w = 0; w < pd->Num_Species_Eqn; w++) {
      electrode_species_source(w, time, dt);
      h -= F * mp->charge_number[w] * mp->species_source[w];
    }

    ai0_anode = mp->u_reaction_rate[0];
    ai0_cathode = mp->u_reaction_rate[2];
    mn = ei[pg->imtrx]->mn;

    if (mn == 0) /* KSC: 2/17/99 */
    {
      ai0 = ai0_anode; /* exchange current density for the anode (one rxn only) */
    } else if (mn == 2) {
      ai0 = ai0_cathode; /* exchange current density for the cathode (one rxn only) */
    } else {
      ai0 = 0.0; /* no electrochemical rxn for the separator */
    }

    var = TEMPERATURE;
    if (d_h != NULL && pd->e[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_h->T[j] = -F * mp->charge_number[w] * (ai0 / R / TT) * bf[var]->phi[j];
      }
    }

    var = VOLTAGE;
    if (d_h != NULL && pd->e[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_h->V[j] = -F * mp->charge_number[w] * (-ai0 / R / TT) * bf[var]->phi[j];
      }
    }
  } else if (mp->HeatSourceModel == BUTLER_VOLMER) /* added by KSC: 4/28/06 */
  {
    h = butler_volmer_heat_source(d_h, mp->u_heat_source);
  } else if (mp->HeatSourceModel == JOULE) {
    h = joule_heat_source(d_h, tt);
  } else if (mp->HeatSourceModel == VISC_DISS) {
    h = visc_diss_heat_source(d_h, mp->u_heat_source);
  } else if (mp->HeatSourceModel == VISC_ACOUSTIC) {
    h = visc_diss_heat_source(d_h, mp->u_heat_source);
    h_acous = visc_diss_acoustic_source(d_h, mp->u_heat_source, mp->len_u_heat_source);
    h += h_acous;
  } else if (mp->HeatSourceModel == EM_DISS) {
    h = em_diss_heat_source(d_h, mp->u_heat_source, mp->len_u_heat_source);
  } else if (mp->HeatSourceModel == EM_VECTOR_DISS) {
    h = em_diss_e_curlcurl_source(d_h, mp->u_heat_source, mp->len_u_heat_source);
  } else if (mp->HeatSourceModel == EPOXY) {
    h = epoxy_heat_source(d_h, tt, dt);
  } else if (mp->HeatSourceModel == VARIABLE_DENSITY) {
    h = vary_rho_heat_source(d_h, tt, dt);
  } else if (mp->HeatSourceModel == HS_FOAM) {
    h = foam_heat_source(d_h, tt, dt);
  } else if (mp->HeatSourceModel == HS_FOAM_PBE) {
    h = foam_pbe_heat_source(d_h, tt, dt);
  } else if (mp->HeatSourceModel == HS_FOAM_PMDI_10) {
    h = foam_pmdi_10_heat_source(d_h, time, tt, dt);
  } else if (mp->HeatSourceModel == USER_GEN) {
    if (d_h == NULL) {
      dbl dhdT[MDE];
      dbl dhdX[DIM][MDE];
      dbl dhdv[DIM][MDE];
      dbl dhdC[MAX_CONC][MDE];
      dbl dhdVolt[MDE];
      usr_heat_source_gen(&h, dhdT, dhdX, dhdv, dhdC, dhdVolt, mp->u_heat_source, time);
    } else {
      usr_heat_source_gen(&h, d_h->T, d_h->X, d_h->v, d_h->C, d_h->V, mp->u_heat_source, time);
    }
  } else {
    GOMA_EH(GOMA_ERROR, "Unrecognized heat source model");
  }

  if (ls != NULL && mp->mp2nd != NULL && mp->mp2nd->HeatSourceModel == CONSTANT &&
      (pd->e[pg->imtrx][R_ENERGY] & T_SOURCE)) {
    /* kludge for solidification tracking with phase function 0 */
    if (pfd != NULL && pd->e[pg->imtrx][R_EXT_VELOCITY]) {
      ls_old = ls;
      ls = pfd->ls[0];
      ls_modulate_heatsource(&h, mp->mp2nd->heatsource_phase[0], ls->Length_Scale,
                             (double)mp->mp2nd->heatsourcemask[0],
                             (double)mp->mp2nd->heatsourcemask[1], d_h);
      ls = ls_old;
    }
    ls_modulate_heatsource(&h, mp->mp2nd->heatsource, ls->Length_Scale,
                           (double)mp->mp2nd->heatsourcemask[0],
                           (double)mp->mp2nd->heatsourcemask[1], d_h);
  }

  return (h);
}
int ls_modulate_heatsource(double *f,
                           double f2,
                           double width,
                           double pm_minus,
                           double pm_plus,
                           HEAT_SOURCE_DEPENDENCE_STRUCT *df) {
  int i, b, var;
  int dim = pd->Num_Dim;
  double factor;
  double f1 = *f;

  if (df == NULL) {
    *f = ls_modulate_property(f1, f2, width, pm_minus, pm_plus, NULL, &factor);

    return (0);
  } else {
    *f = ls_modulate_property(f1, f2, width, pm_minus, pm_plus, df->F, &factor);

    if (pd->v[pg->imtrx][var = TEMPERATURE]) {
      for (i = 0; i < ei[pg->imtrx]->dof[var]; i++) {
        df->T[i] *= factor;
      }
    }

    if (pd->v[pg->imtrx][var = MESH_DISPLACEMENT1]) {
      for (b = 0; b < dim; b++) {
        for (i = 0; i < ei[pg->imtrx]->dof[var + b]; i++) {
          df->X[b][i] *= factor;
        }
      }
    }

    if (pd->v[pg->imtrx][var = VELOCITY1]) {
      for (b = 0; b < dim; b++) {
        for (i = 0; i < ei[pg->imtrx]->dof[var]; i++) {
          df->v[b][i] *= factor;
        }
      }
    }

    if (pd->v[pg->imtrx][var = MASS_FRACTION]) {
      for (b = 0; b < pd->Num_Species; b++) {
        for (i = 0; i < ei[pg->imtrx]->dof[var]; i++) {
          df->C[b][i] *= factor;
        }
      }
    }

    if (pd->v[pg->imtrx][var = VOLTAGE]) {
      for (i = 0; i < ei[pg->imtrx]->dof[var]; i++) {
        df->V[i] *= factor;
      }
    }

    if (pd->v[pg->imtrx][var = POLYMER_STRESS11]) {
      GOMA_WH(-1, "LS modulation of heat source sensitivity wrt to polymer stress dofs not "
                  "implemented.");
    }
  }
  return (0);
}
int assemble_ls_latent_heat_source(double iso_therm,
                                   double latent_heat,
                                   double dt, /* current value of the time step  */
                                   double tt,
                                   double time,
                                   int bc_input_id,
                                   struct Boundary_Condition *BC_Types) {
  int i, j, ii, ledof;
  int eqn, peqn, var, pvar, b, w;

  struct Basis_Functions *bfm;
  double wt, det_J, h3, phi_i, phi_j, flux;
  double source;

  double sign = 1.;
  double vnorm;
  NORMAL_VELOCITY_DEPENDENCE_STRUCT d_vnorm_struct;
  NORMAL_VELOCITY_DEPENDENCE_STRUCT *d_vnorm = &d_vnorm_struct;

  /* struct Boundary_Condition *fluxbc; */

  eqn = R_ENERGY;
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

  /*
  fluxbc = BC_Types + bc_input_id;
  compute_leak_velocity_heat(&vnorm, d_vnorm, tt, dt, NULL, fluxbc);  */

  vnorm = fv->ext_v;
  flux = latent_heat * vnorm;
  memset(d_vnorm->v, 0, sizeof(dbl) * DIM * MDE);
  memset(d_vnorm->T, 0, sizeof(dbl) * MDE);
  memset(d_vnorm->C, 0, sizeof(dbl) * MAX_CONC * MDE);
  memset(d_vnorm->F, 0, sizeof(dbl) * MDE);
  memset(d_vnorm->X, 0, sizeof(dbl) * MDE * DIM);

  /*
  if ( fv->c[wspec] > 1. )
    fprintf(stderr,"flux=%g, Y_w=%g, mass_flux=%g, extv=%g,
  vnorm=%g\n",flux,fv->c[wspec],-mp->mass_flux[wspec],fv->ext_v,vnorm);
  */

  /* finite difference calculation of path dependencies for
     subelement integration
   */
  if (ls->CalcSurfDependencies) {
    peqn = TEMPERATURE;
    bfm = bf[eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][i];

      if (ei[pg->imtrx]->active_interp_ledof[ledof]) {
        ii = ei[pg->imtrx]->lvdof_to_row_lvdof[eqn][i];

        phi_i = bfm->phi[i];

        source = phi_i * flux * lsi->delta;

        source *= det_J * wt;

        source *= h3;

        /*source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];*/

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
    return (0);
  }

  /*
   * Wesiduals
   * ________________________________________________________________________________
   */
  if (af->Assemble_Residual) {
    peqn = TEMPERATURE;
    bfm = bf[eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][i];

      if (ei[pg->imtrx]->active_interp_ledof[ledof]) {
        ii = ei[pg->imtrx]->lvdof_to_row_lvdof[eqn][i];

        phi_i = bfm->phi[i];

        source = phi_i * flux * lsi->delta;

        source *= det_J * wt;

        source *= h3;

        /*source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];*/

        lec->R[LEC_R_INDEX(peqn, ii)] += source;
      }
    }
  }

  if (af->Assemble_Jacobian) {
    peqn = TEMPERATURE;
    bfm = bf[eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][i];

      if (ei[pg->imtrx]->active_interp_ledof[ledof]) {
        ii = ei[pg->imtrx]->lvdof_to_row_lvdof[eqn][i];

        phi_i = bfm->phi[i];

        /*
         * J_W_T
         */
        var = TEMPERATURE;

        if (pd->v[pg->imtrx][var]) {
          pvar = upd->vp[pg->imtrx][var];

          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];

            source = latent_heat * d_vnorm->T[j] * sign;

            source *= phi_i * lsi->delta;

            source *= det_J * wt * h3;
            /*source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];*/

            lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
          }
        }

        /*
         * J_W_vext
         */
        var = EXT_VELOCITY;

        if (pd->v[pg->imtrx][var]) {
          pvar = upd->vp[pg->imtrx][var];

          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];

            source = latent_heat * phi_j * sign;

            source *= phi_i * lsi->delta;

            source *= det_J * wt * h3;
            /*source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];*/

            lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
          }
        }

        /*
         * J_W_w
         */
        var = MASS_FRACTION;

        if (pd->v[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];

            for (w = 0; w < pd->Num_Species_Eqn; w++) {
              pvar = MAX_PROB_VAR + w;

              source = latent_heat * d_vnorm->C[w][j] * sign;

              source *= phi_i * lsi->delta;

              source *= det_J * wt * h3;
              /*source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];*/

              lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
            }
          }
        }

        /*
         * J_W_v
         */
        for (b = 0; b < WIM; b++) {
          var = VELOCITY1 + b;
          if (pd->v[pg->imtrx][var]) {
            pvar = upd->vp[pg->imtrx][var];
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              phi_j = bf[var]->phi[j];

              source = latent_heat * d_vnorm->v[b][j] * sign;

              source *= phi_i * lsi->delta;

              source *= det_J * wt * h3;
              /*source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];*/

              lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
            }
          }
        }

        /*
         * J_W_F
         */
        var = LS;

        if (pd->v[pg->imtrx][var]) {
          pvar = upd->vp[pg->imtrx][var];

          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];

            source = latent_heat * d_vnorm->F[j] * sign;

            source *= phi_i * lsi->delta;

            source += flux * phi_i * lsi->d_delta_dF[j];

            source *= det_J * wt * h3;
            /*source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];*/

            lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
          }
        }
      }
    }
  }
  return (1);
} /*
   * double visc_diss_acoustic_source (dh, param)
   *
   * ------------------------------------------------------------------------------
   * This routine is responsible for filling up the following forces and sensitivities
   * at the current gauss point:
   *     intput:
   *
   *     output:  h             - heat source
   *              d_h->T[j]    - derivative wrt temperature at node j.
   *              d_h->V[j]    - derivative wrt voltage at node j.
   *              d_h->C[i][j] - derivative wrt mass frac species i at node j
   *              d_h->v[i][j] - derivative wrt velocity component i at node j
   *              d_h->X[0][j] - derivative wrt mesh displacement components i at node j
   *              d_h->S[m][i][j][k] - derivative wrt stress mode m of component ij at node k
   *
   *   NB: The user need only supply f, dfdT, dfdC, etc....mp struct is loaded up for you
   * ---------------------------------------------------------------------------
   */

double visc_diss_acoustic_source(HEAT_SOURCE_DEPENDENCE_STRUCT *d_h,
                                 dbl *param,
                                 int num_const) /* General multipliers */
{
  /* Local Variables */
  int var, err;

  int w, j;
  double omega, visc_first, visc_second, R_gas;
  double h, temp1, temp3, ap_square;
  double R; /* Acoustic impedance. */
  CONDUCTIVITY_DEPENDENCE_STRUCT d_R_struct;
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_R = &d_R_struct;

  double k; /* Acoustic wave number. */
  CONDUCTIVITY_DEPENDENCE_STRUCT d_k_struct;
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_k = &d_k_struct;

  double alpha; /* Acoustic Absorption */
  CONDUCTIVITY_DEPENDENCE_STRUCT d_alpha_struct;
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_alpha = &d_alpha_struct;

  double visc_cmb; /* Combined viscosity term  */
  VISCOSITY_DEPENDENCE_STRUCT d_visc_cmb_struct;
  VISCOSITY_DEPENDENCE_STRUCT *d_visc_cmb = &d_visc_cmb_struct;

  dbl gamma[DIM][DIM];
  dbl mu;
  VISCOSITY_DEPENDENCE_STRUCT d_mu_struct; /* viscosity dependence */
  VISCOSITY_DEPENDENCE_STRUCT *d_mu = &d_mu_struct;

  double time = tran->time_value;

  /* Begin Execution */

  /**********************************************************/
  /* Source constant * viscosity* gammadot .. grad_v */

  /* get mu and grad_mu */
  omega = upd->Acoustic_Frequency;
  R = acoustic_impedance(d_R, time);
  k = wave_number(d_k, time);
  alpha = acoustic_absorption(d_alpha, time);

  /* modulate viscosity term according to the LS field
          (i.e. 4/3 shear viscosity plus bulk viscosity)
  */
  visc_first = param[1];
  visc_second = param[2];
  memset(gamma, 0, sizeof(dbl) * DIM * DIM);
  /*  There are some options for the appropriate shear rate which
          to evaluate the viscosity - could choose zero-shear-rate,
          oscillation frequency or rate based on grad(v).
      some of the options require extra derivatives be calculated.
          I'll use oscillation frequency for now.
  */
  /* First evaluate bulk viscosity at rest conditions
     Then evaluate shear viscosity at frequency.
     The visc_first parameter is a bulk viscosity multiplier  */

  visc_cmb = viscosity(gn, gamma, d_visc_cmb);
  gamma[0][0] = omega;
  mu = viscosity(gn, gamma, d_mu);
  visc_cmb = 4. * mu / 3. + visc_first * visc_cmb;
  for (j = 0; j < ei[pg->imtrx]->dof[TEMPERATURE]; j++) {
    d_visc_cmb->T[j] *= visc_first * d_mu->T[j];
    d_visc_cmb->T[j] += (4. / 3.) * d_mu->T[j];
  }
  {
    d_visc_cmb->F[j] *= visc_first * d_mu->F[j];
    d_visc_cmb->F[j] += (4. / 3.) * d_mu->F[j];
  }
  {
    for (w = 0; w < pd->Num_Species_Eqn; w++) {
      d_visc_cmb->C[w][j] *= visc_first * d_mu->C[w][j];
      d_visc_cmb->C[w][j] += (4. / 3.) * d_mu->C[w][j];
    }
  }
  if (ls != NULL) {
    err = ls_modulate_viscosity(
        &visc_cmb, visc_second, ls->Length_Scale, (double)mp->mp2nd->viscositymask[0],
        (double)mp->mp2nd->viscositymask[1], d_visc_cmb, mp->mp2nd->ViscosityModel);
    GOMA_EH(err, "ls_modulate_viscosity");
    /*  optionally set impedance to true value input on card	*/
    if (num_const == 4) {
      R_gas = param[3];
      R = ls_modulate_thermalconductivity(mp->acoustic_impedance, R_gas, ls->Length_Scale,
                                          (double)mp->mp2nd->acousticimpedancemask[0],
                                          (double)mp->mp2nd->acousticimpedancemask[1], d_R);
    }
  }

  ap_square = SQUARE(fv->apr) + SQUARE(fv->api);
  temp1 = SQUARE(k) / SQUARE(R);
  temp3 = (1. + 4. * SQUARE(alpha));
  h = ap_square * temp1 * visc_cmb * temp3;
  h *= param[0];
#if 0
fprintf(stderr,"visc %g %g %g %g %g %g %g\n",mu,visc_cmb,gamma[0][0], visc_first, temp1, temp3, h);
#endif

  /* Now do sensitivies */
  if (af->Assemble_Jacobian) {

    var = TEMPERATURE;
    if (d_h != NULL && pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_h->T[j] += param[0] * ap_square *
                     (visc_cmb * (temp1 * 4. * alpha * d_alpha->T[j] +
                                  temp3 * 2. * k * (R * d_k->T[j] - k * d_R->T[j]) / (R * R * R)) +
                      d_visc_cmb->T[j] * temp1 * temp3);
      }
    }

    var = MASS_FRACTION;
    if (d_h != NULL && pd->v[pg->imtrx][var]) {
      for (w = 0; w < pd->Num_Species_Eqn; w++) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_h->C[w][j] +=
              param[0] * ap_square *
              (visc_cmb * (temp1 * 4. * alpha * d_alpha->C[w][j] +
                           temp3 * 2. * k * (R * d_k->C[w][j] - k * d_R->C[w][j]) / (R * R * R)) +
               d_visc_cmb->C[w][j] * temp1 * temp3);
        }
      }
    }

    var = FILL;
    if (d_h != NULL && pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_h->F[j] += param[0] * ap_square *
                     (visc_cmb * (temp1 * 4. * alpha * d_alpha->F[j] +
                                  temp3 * 2. * k * (R * d_k->F[j] - k * d_R->F[j]) / (R * R * R)) +
                      temp1 * temp3 * d_visc_cmb->F[j]);
      }
    }

    var = ACOUS_PREAL;
    if (d_h != NULL && pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_h->APR[j] += param[0] * temp1 * visc_cmb * temp3 * 2. * fv->apr * bf[var]->phi[j];
      }
    }
    var = ACOUS_PIMAG;
    if (d_h != NULL && pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_h->API[j] += param[0] * temp1 * visc_cmb * temp3 * 2. * fv->api * bf[var]->phi[j];
      }
    }

  } /* end of if Assemble Jacobian  */

  return (h);
} /*
   * double em_diss_acoustic_source (dh, param)
   *
   * ------------------------------------------------------------------------------
   * This routine is responsible for filling up the following forces and sensitivities
   * at the current gauss point:
   *     intput:
   *
   *     output:  h             - heat source
   *              d_h->T[j]    - derivative wrt temperature at node j.
   *              d_h->V[j]    - derivative wrt voltage at node j.
   *              d_h->C[i][j] - derivative wrt mass frac species i at node j
   *              d_h->v[i][j] - derivative wrt velocity component i at node j
   *              d_h->X[0][j] - derivative wrt mesh displacement components i at node j
   *              d_h->S[m][i][j][k] - derivative wrt stress mode m of component ij at node k
   *
   *   NB: The user need only supply f, dfdT, dfdC, etc....mp struct is loaded up for you
   * ---------------------------------------------------------------------------
   */

double em_diss_heat_source(HEAT_SOURCE_DEPENDENCE_STRUCT *d_h,
                           dbl *param,
                           int num_const) /* General multipliers */
{
  /* Local Variables */
  int var;

  int w, j;
  double h, ap_square;

  double R; /* Acoustic impedance. */
  CONDUCTIVITY_DEPENDENCE_STRUCT d_R_struct;
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_R = &d_R_struct;

  double k; /* Acoustic wave number. */
  CONDUCTIVITY_DEPENDENCE_STRUCT d_k_struct;
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_k = &d_k_struct;

  double alpha; /* Acoustic Absorption */
  CONDUCTIVITY_DEPENDENCE_STRUCT d_alpha_struct;
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_alpha = &d_alpha_struct;

  double time = tran->time_value;

  /* Begin Execution */

  /**********************************************************/
  /* Source constant * viscosity* gammadot .. grad_v */

  /* get mu and grad_mu */
  R = acoustic_impedance(d_R, time);
  k = wave_number(d_k, time);
  alpha = acoustic_absorption(d_alpha, time);

  ap_square = SQUARE(fv->apr) + SQUARE(fv->api);
  h = ap_square * alpha * k / R;
  GOMA_ASSERT_ALWAYS(num_const > 0);
  h *= param[0];

  /* Now do sensitivies */
  if (af->Assemble_Jacobian) {

    var = TEMPERATURE;
    if (d_h != NULL && pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_h->T[j] += param[0] * ap_square *
                     (R * (alpha * d_k->T[j] + k * d_alpha->T[j]) - alpha * k * d_R->T[j]) /
                     SQUARE(R);
      }
    }

    var = MASS_FRACTION;
    if (d_h != NULL && pd->v[pg->imtrx][var]) {
      for (w = 0; w < pd->Num_Species_Eqn; w++) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_h->C[w][j] +=
              param[0] * ap_square *
              (R * (alpha * d_k->C[w][j] + k * d_alpha->C[w][j]) - alpha * k * d_R->C[w][j]) /
              SQUARE(R);
        }
      }
    }

    var = FILL;
    if (d_h != NULL && pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_h->F[j] += param[0] * ap_square *
                     (R * (alpha * d_k->F[j] + k * d_alpha->F[j]) - alpha * k * d_R->F[j]) /
                     SQUARE(R);
      }
    }

    var = ACOUS_PREAL;
    if (d_h != NULL && pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_h->APR[j] += param[0] * alpha * k / R * 2. * fv->apr * bf[var]->phi[j];
      }
    }
    var = ACOUS_PIMAG;
    if (d_h != NULL && pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_h->API[j] += param[0] * alpha * k / R * 2. * fv->api * bf[var]->phi[j];
      }
    }

  } /* end of if Assemble Jacobian  */

  return (h);
} /*
   * double em_diss_e_curlcurl_source (dh, param)
   *
   * ------------------------------------------------------------------------------
   * This routine is responsible for filling up the following forces and sensitivities
   * at the current gauss point:
   *     intput:
   *
   *     output:  h             - heat source
   *              d_h->T[j]    - derivative wrt temperature at node j.
   *              d_h->V[j]    - derivative wrt voltage at node j.
   *              d_h->C[i][j] - derivative wrt mass frac species i at node j
   *              d_h->v[i][j] - derivative wrt velocity component i at node j
   *              d_h->X[0][j] - derivative wrt mesh displacement components i at node j
   *              d_h->S[m][i][j][k] - derivative wrt stress mode m of component ij at node k
   *
   *   NB: The user need only supply f, dfdT, dfdC, etc....mp struct is loaded up for you
   * ---------------------------------------------------------------------------
   */

double em_diss_e_curlcurl_source(HEAT_SOURCE_DEPENDENCE_STRUCT *d_h,
                                 dbl *param,
                                 int num_const) /* General multipliers */
{
  /* Local Variables */
  int Rvar, Ivar;

  int j;
  double h;

  double time = tran->time_value;

  /* Begin Execution */

  /**********************************************************/
  /* Source = 1/2 * Re(E cross[conj{curl[E]}])
            = 1/2 * omega*epsilon_0*epsilon_c'' * (mag(E))^2 */

  double omega, epsilon_0, epsilon_cpp, h_factor, mag_E_sq;

  dbl n; /* Refractive index. */
  CONDUCTIVITY_DEPENDENCE_STRUCT d_n_struct;
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_n = &d_n_struct;

  dbl k; /* Extinction coefficient */
  CONDUCTIVITY_DEPENDENCE_STRUCT d_k_struct;
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_k = &d_k_struct;

  omega = upd->Acoustic_Frequency;

  epsilon_0 = mp->permittivity;

  n = refractive_index(d_n, time);
  k = extinction_index(d_k, time);

  epsilon_cpp = 2 * n * k;
  h_factor = omega * epsilon_0 * epsilon_cpp / 2.;

  mag_E_sq = 0;
  for (int p = 0; p < pd->Num_Dim; p++) {
    mag_E_sq += SQUARE(fv->em_er[p]) + SQUARE(fv->em_ei[p]);
  }

  h = h_factor * mag_E_sq;
  GOMA_ASSERT_ALWAYS(num_const > 0);
  h *= param[0];

  /* Now sensitivies */
  if (af->Assemble_Jacobian) {

    for (int p = 0; p < pd->Num_Dim; p++) {
      Rvar = EM_E1_REAL + p;
      Ivar = EM_E1_IMAG + p;

      if (d_h != NULL && pd->v[pg->imtrx][Rvar]) {
        for (j = 0; j < ei[pg->imtrx]->dof[Rvar]; j++) {
          d_h->EM_ER[p][j] += param[0] * h_factor * 2.0 * fv->em_er[p] * bf[Rvar]->phi[j];
        }
      }

      if (d_h != NULL && pd->v[pg->imtrx][Ivar]) {
        for (j = 0; j < ei[pg->imtrx]->dof[Ivar]; j++) {
          d_h->EM_EI[p][j] += param[0] * h_factor * 2.0 * fv->em_ei[p] * bf[Ivar]->phi[j];
        }
      }
    }

  } /* end of if Assemble Jacobian  */

  return (h);
}