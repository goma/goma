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
#include "load_field_variables.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* GOMA include files */
#define GOMA_MM_FILL_TERMS_C
#include "mm_fill_terms.h"

#include "ac_particles.h"
#include "az_aztec.h"
#include "bc_colloc.h"
#include "bc_contact.h"
#include "density.h"
#include "el_elm.h"
#include "el_elm_info.h"
#include "el_geom.h"
#include "exo_struct.h"
#include "mm_as.h"
#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "mm_dil_viscosity.h"
#include "mm_eh.h"
#include "mm_fill_aux.h"
#include "mm_fill_common.h"
#include "mm_fill_continuity.h"
#include "mm_fill_energy.h"
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

/*  _______________________________________________________________________  */

/* assemble_mesh -- assemble terms (Residual &| Jacobian) for mesh stress eqns
 *
 * in:
 * 	ei -- pointer to Element Indeces		structure
 *	pd -- pointer to Problem Description		structure
 *	af -- pointer to Action Flag			structure
 *	bf -- pointer to Basis Function			structure
 *	fv -- pointer to Field Variable			structure
 *	cr -- pointer to Constitutive Relation		structure
 *	mp -- pointer to Material Property		structure
 *	esp-- pointer to Element Stiffness Pointers	structure
 *	lec-- pointer to Local Element Contribution	structure
 *
 * out:
 *	a   -- gets loaded up with proper contribution
 *	lec -- gets loaded up with local contributions to resid, Jacobian
 * 	r   -- residual RHS vector
 *
 * Created:	Thu Mar  3 07:48:01 MST 1994 pasacki@sandia.gov
 *
 * Revised:	Sat Mar 19 16:07:51 MST 1994 pasacki@sandia.gov
 *
 *
 * Note: currently we do a "double load" into the addresses in the global
 *       "a" matrix and resid vector held in "esp", as well as into the
 *       local accumulators in "lec".
 *
 */

int assemble_mesh(double time,
                  double tt,
                  double dt,
                  int ielem,    /* current element number */
                  int ip,       /* current integration point */
                  int ip_total) /* total gauss integration points */

{
  int eqn, peqn, var, pvar;
  const int dim = pd->Num_Dim;
  int p, q, a, b, imtrx;

  int w;

  int i, j;
  int status, err;

  dbl det_J;
  dbl d_det_J_dmeshbj;
  dbl h3;          /* Volume element (scale factors). */
  dbl dh3dmesh_bj; /* Sensitivity to (b,j) mesh dof. */

  dbl VV[DIM][DIM]; /* Inertia Tensor... */

  double vconv[MAX_PDIM];     /*Calculated convection velocity */
  double vconv_old[MAX_PDIM]; /*Calculated convection velocity at previous time*/
  CONVECTION_VELOCITY_DEPENDENCE_STRUCT d_vconv_struct;
  CONVECTION_VELOCITY_DEPENDENCE_STRUCT *d_vconv = &d_vconv_struct;

  dbl TT[DIM][DIM]; /* Mesh stress tensor... */

  dbl dTT_dx[DIM][DIM][DIM][MDE];      /* Sensitivity of stress tensor...
                          to nodal displacements */
  dbl dTT_dp[DIM][DIM][MDE];           /* Sensitivity of stress tensor...
                                  to nodal pressure*/
  dbl dTT_dc[DIM][DIM][MAX_CONC][MDE]; /* Sensitivity of stress tensor...
                          to nodal concentration */
  dbl dTT_dp_liq[DIM][DIM][MDE];       /* Sensitivity of stress tensor...
                                           to nodal porous liquid pressure*/
  dbl dTT_dp_gas[DIM][DIM][MDE];       /* Sensitivity of stress tensor...
                                           to nodal porous gas pressure*/
  dbl dTT_dporosity[DIM][DIM][MDE];    /* Sensitivity of stress tensor...
                                        to nodal porosity*/
  dbl dTT_dsink_mass[DIM][DIM][MDE];   /* Sensitivity of stress tensor...
                                          to nodal sink_mass*/
  dbl dTT_dT[DIM][DIM][MDE];           /* Sensitivity of stress tensor...
                                      to temperature*/
  dbl dTT_dmax_strain[DIM][DIM][MDE];  /* Sensitivity of stress tensor...
                             to max_strain*/
  dbl dTT_dcur_strain[DIM][DIM][MDE];  /* Sensitivity of stress tensor...
                             to cur_strain*/

  dbl mu;
  dbl lambda, rho;

  dbl g[DIM]; /* Mesh body force. */

  dbl diff_a, diff_b, diff_c; /* Temporary variables hold partially */
  /* constructed diffusion terms... */
  dbl diffusion = 0., advection = 0., advect_a = 0., advect_b = 0., advect_c = 0.;
  dbl source = 0., mass = 0.;
  dbl wt;

  /*
   * Galerkin weighting functions for i-th and a-th mesh stress residuals
   * and some of their derivatives...
   */

  dbl phi_i;
  dbl(*grad_phi_i_e_a)[DIM];

  /*
   * Interpolation functions for variables and some of their derivatives.
   */

  dbl phi_j;

  /*
   * Mesh derivatives...
   */
  dbl(*dgradphi_i_e_a_dmesh)[DIM][DIM][MDE]; /* for specific (i, a, b, j) !!  */

  /* density derivatives */
  DENSITY_DEPENDENCE_STRUCT d_rho_struct; /* density dependence */
  DENSITY_DEPENDENCE_STRUCT *d_rho = &d_rho_struct;

  struct Basis_Functions *bfm;

  /*finally we have some time-derivative quantities for the Newmark scheme */
  dbl x_c[DIM];           /*base coordinates of deformed mesh */
  dbl x_old[DIM];         /*old value of base coordinates on deformed mesh */
  dbl x_dot_old[DIM];     /*Old value of mesh velocity */
  dbl x_dbl_dot_old[DIM]; /*old value of acceleration */
  dbl x_dbl_dot[DIM];     /*current value of mesh acceleration */
  dbl newmark_beta = 0.0; /*Newmark scheme beta value. */

  int transient_run = pd->TimeIntegration != STEADY;
  int mass_on;
  int advection_on = 0;
  int source_on = 0;
  int diffusion_on = 0;

  dbl mass_etm, advection_etm, diffusion_etm, source_etm;

  dbl d_area;

  dbl *_J;

  status = 0;

  /*
   * Unpack variables from structures for local convenience...
   */

  eqn = R_MESH1; /* Well, yes, there really are 3, */
                 /* but for now just use the first */
                 /* to tell whether there is anything */
                 /* to do...*/
  wt = fv->wt;
  h3 = fv->h3; /* Differential volume element (scales). */

  /*
   * Bail out fast if there's nothing to do...
   * Also bail out if we are in an a shell equation. In that case Num_Dim
   *      is equal to the problem dimension, but ielem_dim is one less than the
   *      problem dimension.
   */
  if (!pd->e[pg->imtrx][eqn] || (ei[pg->imtrx]->ielem_dim < pd->Num_Dim)) {
    return (status);
  }

  det_J = bf[eqn]->detJ;
  d_area = det_J * wt * h3;

  /*
   * Material property constants, etc. Any variations for this
   * Gauss point were evaluated in mm_fill.
   */
  rho = density(d_rho, time);
  mu = elc->lame_mu;
  lambda = elc->lame_lambda;

  if (mp->MeshSourceModel == CONSTANT) {
    for (a = 0; a < dim; a++) {
      g[a] = mp->mesh_source[a];
    }
  } else {
    GOMA_EH(GOMA_ERROR, "Unrecognized MeshSourceModel");
  }

  /*
   * Get the deformation gradients and tensors if needed
   */
  err = belly_flop(mu);
  GOMA_EH(err, "error in belly flop");
  if (err == 2)
    return (err);

  /*
   * If using volume change metric for element quality,
   * sum local contribution here.
   */
  if (nEQM > 0 && eqm->do_vol && af->Assemble_Residual) {
    eqm->vol_sum += fv->volume_change;
    if (fv->volume_change < eqm->vol_low) {
      eqm->vol_low = fv->volume_change;
    }
    eqm->vol_count++;
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
                           dTT_dsink_mass, dTT_dT, dTT_dmax_strain, dTT_dcur_strain, mu, lambda, dt,
                           ielem, ip, ip_total);

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
  for (a = 0; a < dim; a++) {
    for (b = 0; b < dim; b++) {
      VV[a][b] = vconv[a] * vconv[b];
    }
  }

  if (tran->solid_inertia) {
    newmark_beta = tran->newmark_beta;
    for (a = 0; a < dim; a++) {
      x_c[a] = fv->x[a];
      x_old[a] = fv_old->x[a];
      x_dot_old[a] = fv_dot_old->x[a];
      x_dbl_dot_old[a] = fv_dot_dot_old->x[a];

      /*generalized Newmark equations here */
      /*Notice Gama is not needed here, unless we have a damping term */
      x_dbl_dot[a] = (x_c[a] - x_old[a]) / newmark_beta / dt / dt -
                     x_dot_old[a] * (1.0 / newmark_beta / dt) -
                     x_dbl_dot_old[a] * (1 - 2. * newmark_beta) / 2. / newmark_beta;
    }
  } else {
    for (a = 0; a < dim; a++) {
      x_dbl_dot[a] = 0.;
    }
  }

  /*
   * Residuals_________________________________________________________________
   */

  if (af->Assemble_Residual) {
    /*
     * Assemble each component "a" of the momentum equation...
     */
    for (a = 0; a < dim; a++) {
      eqn = R_MESH1 + a;
      peqn = upd->ep[pg->imtrx][eqn];
      bfm = bf[eqn];

      mass_on = pd->e[pg->imtrx][eqn] & T_MASS;
      advection_on = pd->e[pg->imtrx][eqn] & T_ADVECTION;
      diffusion_on = pd->e[pg->imtrx][eqn] & T_DIFFUSION;
      source_on = pd->e[pg->imtrx][eqn] & T_SOURCE;

      mass_etm = pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
      advection_etm = pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
      diffusion_etm = pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
      source_etm = pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

      for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
        phi_i = bfm->phi[i];

        grad_phi_i_e_a = bfm->grad_phi_e[i][a];

        mass = 0.;
        if (transient_run && pd_glob[ei[pg->imtrx]->mn]->MeshMotion == DYNAMIC_LAGRANGIAN) {
          if (pd->e[pg->imtrx][eqn] & T_MASS) {
            mass = -x_dbl_dot[a];
            mass *= phi_i * rho * d_area;
            mass *= mass_etm;
          }
        }

        diffusion = 0.;
        if (diffusion_on) {
          if (cr->MeshMotion == ARBITRARY) {
            /*
             * use pseudo cartesian arbitrary mesh motion
             */
            for (p = 0; p < dim; p++) {
              diffusion += bfm->d_phi[i][p] * TT[a][p];
            }
            diffusion *= -det_J * wt;
          } else {
            for (q = 0; q < VIM; q++) {
              for (p = 0; p < VIM; p++) {
                diffusion += grad_phi_i_e_a[p][q] * TT[q][p];
              }
            }
            diffusion *= -d_area;
          }

          diffusion *= diffusion_etm;
        }

        /* add inertia of moving mesh */
        advection = 0.;
        if (advection_on) {
          GOMA_WH(-1, "Warning: mesh advection unverified.\n");
          for (p = 0; p < dim; p++) {
            for (q = 0; q < dim; q++) {

              advection += grad_phi_i_e_a[p][q] * VV[q][p];
            }
          }
          advection *= rho * d_area;
          /* Advection term only applies in Lagrangian mesh motion */
          advection *= advection_etm;
        }

        source = 0.;
        if (source_on) {
          source = phi_i * g[a] * d_area;
          /* Source term only applies in Lagrangian mesh motion */

          source *= source_etm;
        }

        /*
         * porous term removed for mesh equation
         *  - the additional effects due  to porosity are entered
         *    into the consitutive equation for stress
         */

        lec->R[LEC_R_INDEX(peqn, i)] += mass + diffusion + source + advection;
      }
    }
  }

  /*
   * Jacobian terms...
   */

  if (af->Assemble_Jacobian) {
    for (a = 0; a < dim; a++) {
      eqn = R_MESH1 + a;
      peqn = upd->ep[pg->imtrx][eqn];
      bfm = bf[eqn];

      mass_on = pd->e[pg->imtrx][eqn] & T_MASS;
      advection_on = pd->e[pg->imtrx][eqn] & T_ADVECTION;
      diffusion_on = pd->e[pg->imtrx][eqn] & T_DIFFUSION;
      source_on = pd->e[pg->imtrx][eqn] & T_SOURCE;

      mass_etm = pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
      advection_etm = pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
      diffusion_etm = pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
      source_etm = pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

      for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
        phi_i = bfm->phi[i];
        grad_phi_i_e_a = bfm->grad_phi_e[i][a];
        dgradphi_i_e_a_dmesh = bfm->d_grad_phi_e_dmesh[i][a];

        /*
         * Set up some preliminaries that are needed for the (a,i)
         * equation for bunches of (b,j) column variables...
         */

        /*
         * J_d_d
         */
        for (b = 0; b < dim; b++) {
          var = MESH_DISPLACEMENT1 + b;
          if (pd->v[pg->imtrx][var]) {
            pvar = upd->vp[pg->imtrx][var];
            _J = &(lec->J[LEC_J_INDEX(peqn, pvar, i, 0)]);

            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              phi_j = bf[var]->phi[j];

              d_det_J_dmeshbj = bfm->d_det_J_dm[b][j];

              dh3dmesh_bj = fv->dh3dq[b] * phi_j;

              mass = 0.;

              if (transient_run && pd_glob[ei[pg->imtrx]->mn]->MeshMotion == DYNAMIC_LAGRANGIAN) {
                if (mass_on) {
                  mass = -(double)delta(a, b) * phi_j / newmark_beta / dt / dt;
                  mass *= phi_i * rho * d_area;
                  mass -= x_dbl_dot[a] * phi_i * rho * wt *
                          (d_det_J_dmeshbj * h3 + dh3dmesh_bj * det_J);
                  mass *= mass_etm;
                }
              }

              diffusion = 0.;
              /* Three parts:
               *   diff_a = Int ( d(grad(phi_i e_a))/dmesh : Pi |Jv| )
               *   diff_b = Int ( grad(phi_i e_a) : d(Pi)/dmesh |Jv| )
               *   diff_c = Int ( grad(phi_i e_a) : Pi d(|Jv|)/dmesh )
               */
              if (diffusion_on) {
                if (cr->MeshMotion == ARBITRARY) {
                  /*
                   * use pseudo cartesian arbitrary mesh motion
                   */
                  diff_a = 0.;
                  diff_b = 0.;
                  diff_c = 0.;

                  for (p = 0; p < dim; p++) {
                    diff_a += bfm->d_d_phi_dmesh[i][p][b][j] * TT[a][p];

                    diff_b += bfm->d_phi[i][p] * dTT_dx[a][p][b][j];

                    diff_c += bfm->d_phi[i][p] * TT[a][p];
                  }

                  diff_a *= -det_J * wt;
                  diff_b *= -det_J * wt;
                  diff_c *= -d_det_J_dmeshbj * wt;

                } else {
                  /* LAGRANGIAN MESH - Solve in curvilinear coordinates */

                  diff_a = 0.;
                  diff_c = 0.;
                  diff_b = 0.;

                  for (p = 0; p < VIM; p++) {
                    for (q = 0; q < VIM; q++) {
                      diff_a += dgradphi_i_e_a_dmesh[p][q][b][j] * TT[q][p];

                      diff_b += grad_phi_i_e_a[p][q] * dTT_dx[q][p][b][j];

                      diff_c += grad_phi_i_e_a[p][q] * TT[q][p];
                    }
                  }
                  diff_a *= -d_area;

                  diff_b *= -d_area;

                  diff_c *= -wt * (d_det_J_dmeshbj * h3 + det_J * dh3dmesh_bj);
                }

                diffusion = diff_a + diff_b + diff_c;

                diffusion *= diffusion_etm;
              }

              /* add inertia of moving mesh */
              advect_a = 0.;
              advect_b = 0.;
              advect_c = 0.;
              if (advection_on) {
                /*
                for ( p=0; p<dim; p++)
                  {
                    for ( q=0; q<dim; q++)
                      {
                        dgradphi_i_e_a_dmeshbj[p][q] =
                          bfm->
                          d_grad_phi_e_dmesh[i][a] [p][q] [b][j];
                      }
                  }
                  */
                for (p = 0; p < dim; p++) {
                  for (q = 0; q < dim; q++) {

                    advect_a += grad_phi_i_e_a[p][q] * VV[q][p];
                    advect_b += grad_phi_i_e_a[p][q] *
                                (d_vconv->X[q][b][j] * vconv[p] + vconv[q] * d_vconv->X[p][b][j]);
                    advect_c += dgradphi_i_e_a_dmesh[p][q][b][j] * VV[q][p];
                  }
                }

                advect_a *= rho * (d_det_J_dmeshbj * h3 + det_J * dh3dmesh_bj) * wt;
                advect_b *= rho * det_J * wt;
                advect_c *= rho * det_J * wt;

                advection = advect_a + advect_b + advect_c;
                advection *= advection_etm;
              }

              source = 0.;

              if (source_on) {
                source = phi_i * g[a] * (d_det_J_dmeshbj * h3 + det_J * dh3dmesh_bj) * wt;
                source *= source_etm;
              }

              _J[j] += mass + diffusion + source + advection;
              /*  lec->J[LEC_J_INDEX(peqn,pvar,i,j)] += mass + diffusion + source + advection;*/
            }
          }
        }

        /*
         * J_d_P
         */
        /* For Lagrangian Mesh, add in the sensitivity to pressure  */
        if (pd->MeshMotion == LAGRANGIAN || pd->MeshMotion == DYNAMIC_LAGRANGIAN) {
          var = PRESSURE;
          if (pd->v[pg->imtrx][var]) {
            pvar = upd->vp[pg->imtrx][var];
            _J = &(lec->J[LEC_J_INDEX(peqn, pvar, i, 0)]);

            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              diffusion = 0.;
              phi_j = bf[var]->phi[j];
              for (p = 0; p < VIM; p++) {
                for (q = 0; q < VIM; q++) {
                  diffusion += grad_phi_i_e_a[p][q] * dTT_dp[q][p][j];
                }
              }
              diffusion *= -d_area;
              diffusion *= diffusion_etm;

              _J[j] += diffusion;
              /* lec->J[LEC_J_INDEX(peqn,pvar,i,j)] += diffusion;*/
            }
          }
        }
        /*
         * J_d_c
         */
        var = MASS_FRACTION;
        if (pd->v[pg->imtrx][var]) {

          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            for (w = 0; w < pd->Num_Species_Eqn; w++) {
              _J = &(lec->J[LEC_J_INDEX(peqn, MAX_PROB_VAR + w, i, 0)]);

              /* add inertia of moving mesh */
              advection = 0.;
              /* For Lagrangian Mesh with advection, add in the
                 sensitivity to concentration  */
              if ((cr->MeshMotion == LAGRANGIAN || cr->MeshMotion == DYNAMIC_LAGRANGIAN) &&
                  (pd->e[pg->imtrx][eqn] & T_ADVECTION)) {
                advect_a = 0.;
                advect_b = 0.;

                for (p = 0; p < dim; p++) {
                  for (q = 0; q < dim; q++) {

                    advect_a += grad_phi_i_e_a[p][q] *
                                (d_vconv->C[q][w][j] * vconv[p] + vconv[q] * d_vconv->C[p][w][j]);
                  }
                }
                advect_a *= rho * det_J * wt * h3;
                advect_a *= pd->etm[pg->imtrx][eqn][LOG2_ADVECTION];

                for (p = 0; p < dim; p++) {
                  for (q = 0; q < dim; q++) {
                    advect_b += grad_phi_i_e_a[p][q] * VV[q][p];
                  }
                }
                advect_b *= d_rho->C[w][j] * det_J * wt * h3;
                advect_b *= pd->etm[pg->imtrx][eqn][LOG2_ADVECTION];
                advection = advect_a + advect_b;
              }

              diffusion = 0.;
              if ((cr->MeshMotion == LAGRANGIAN || cr->MeshMotion == DYNAMIC_LAGRANGIAN) &&
                  (pd->e[pg->imtrx][eqn] & T_DIFFUSION)) {
                phi_j = bf[var]->phi[j];
                for (p = 0; p < VIM; p++) {
                  for (q = 0; q < VIM; q++) {
                    grad_phi_i_e_a[p][q] = bfm->grad_phi_e[i][a][p][q];

                    diffusion += grad_phi_i_e_a[p][q] * dTT_dc[q][p][w][j];
                  }
                }
                diffusion *= -det_J * wt * h3;
                diffusion *= pd->etm[pg->imtrx][eqn][LOG2_DIFFUSION];
              }

              _J[j] += advection + diffusion;
              /*lec->J[LEC_J_INDEX(peqn,MAX_PROB_VAR + w,i,j)] += advection + diffusion;*/
            }
          }
        }

        /*
         * J_d_p_liq
         */

        if (pd->v[pg->imtrx][POR_LIQ_PRES] &&
            (cr->MeshMotion == LAGRANGIAN || cr->MeshMotion == DYNAMIC_LAGRANGIAN)) {
          var = POR_LIQ_PRES;
          pvar = upd->vp[pg->imtrx][var];
          _J = &(lec->J[LEC_J_INDEX(peqn, pvar, i, 0)]);

          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            advection = 0.;
            /* For Lagrangian Mesh with advection, add in the
               sensitivity to concentration. N.B. this is zero
               for now. The only way this will be nonzero is if
               we add in diffusion-induced velocity in each porous
               phase. This may happen in multicomponent cases, but not
               now. PRS 5/7/01 */
            if (advection_on) {
              advect_a = 0.;
              advect_b = 0.;

              /*  for ( p=0; p<dim; p++)
                  {
                    for ( q=0; q<dim; q++)
                      {

                        advect_a += grad_phi_i_e_a[p][q] * (
                          dvconv_d_p_gas[q][w][j] * vconv[p] +
                          vconv[q] * dvconv_d_p_gas[p][w][j] );
                      }
                  }*/
              /* advect_a *= rho * det_J * wt * h3;
                 advect_a *= pd->etm[pg->imtrx][eqn][LOG2_ADVECTION]; */

              /*  for ( p=0; p<dim; p++)
                 {
                   for ( q=0; q<dim; q++)
                     {
                       advect_b += grad_phi_i_e_a[p][q] * VV[q][p];
                     }
                 }*/
              /* advect_b *= d_rho->C[w][j] * det_J * wt * h3;
                 advect_b *= pd->etm[pg->imtrx][eqn][LOG2_ADVECTION];
                 advection = advect_a + advect_b;  */
            }

            diffusion = 0.;
            if (diffusion_on) {
              phi_j = bf[var]->phi[j];
              for (p = 0; p < VIM; p++) {
                for (q = 0; q < VIM; q++) {
                  grad_phi_i_e_a[p][q] = bfm->grad_phi_e[i][a][p][q];

                  diffusion += grad_phi_i_e_a[p][q] * dTT_dp_liq[q][p][j];
                }
              }
              diffusion *= -d_area;
              diffusion *= diffusion_etm;
            }

            _J[j] += advection + diffusion;
            /* lec->J[LEC_J_INDEX(peqn,pvar,i,j)] += advection + diffusion;*/
          }
        }
        /*
         * J_d_p_gas
         */
        if (pd->v[pg->imtrx][POR_GAS_PRES] &&
            (cr->MeshMotion == LAGRANGIAN || cr->MeshMotion == DYNAMIC_LAGRANGIAN)) {
          var = POR_GAS_PRES;
          pvar = upd->vp[pg->imtrx][var];

          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            advection = 0.;
            /* For Lagrangian Mesh with advection, add in the
               sensitivity to concentration. N.B. this is zero
               for now. The only way this will be nonzero is if
               we add in diffusion-induced velocity in each porous
               phase. This may happen in multicomponent cases, but not
               now. PRS 5/7/01 */
            if (advection_on) {
              advect_a = 0.;
              advect_b = 0.;

              /*  for ( p=0; p<dim; p++)
                  {
                    for ( q=0; q<dim; q++)
                      {

                        advect_a += grad_phi_i_e_a[p][q] * (
                          dvconv_d_p_gas[q][w][j] * vconv[p] +
                          vconv[q] * dvconv_d_p_gas[p][w][j] );
                      }
                  }*/
              /* advect_a *= rho * det_J * wt * h3;
                 advect_a *= pd->etm[pg->imtrx][eqn][LOG2_ADVECTION]; */

              /*  for ( p=0; p<dim; p++)
                  {
                    for ( q=0; q<dim; q++)
                      {
                         advect_b += grad_phi_i_e_a[p][q] * VV[q][p];
                      }
                  }*/
              /* advect_b *= d_rho->C[w][j] * det_J * wt * h3;
                 advect_b *= advection_etm;
                 advection = advect_a + advect_b;  */
            }

            diffusion = 0.;
            if (diffusion_on) {
              phi_j = bf[var]->phi[j];
              for (p = 0; p < VIM; p++) {
                for (q = 0; q < VIM; q++) {
                  grad_phi_i_e_a[p][q] = bfm->grad_phi_e[i][a][p][q];

                  diffusion += grad_phi_i_e_a[p][q] * dTT_dp_gas[q][p][j];
                }
              }
              diffusion *= -det_J * wt * h3;
              diffusion *= diffusion_etm;
            }

            lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += advection + diffusion;
          }
        }
        /*
         * J_d_porosity
         */
        if (pd->v[pg->imtrx][POR_POROSITY]) {
          var = POR_POROSITY;
          pvar = upd->vp[pg->imtrx][var];
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            diffusion = 0.;
            if ((cr->MeshMotion == LAGRANGIAN || cr->MeshMotion == DYNAMIC_LAGRANGIAN) &&
                (pd->e[pg->imtrx][eqn] & T_DIFFUSION)) {
              phi_j = bf[var]->phi[j];
              for (p = 0; p < VIM; p++) {
                for (q = 0; q < VIM; q++) {
                  grad_phi_i_e_a[p][q] = bfm->grad_phi_e[i][a][p][q];

                  diffusion += grad_phi_i_e_a[p][q] * dTT_dporosity[q][p][j];
                }
              }
              diffusion *= -d_area;
              diffusion *= pd->etm[pg->imtrx][eqn][LOG2_DIFFUSION];
            }

            lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += diffusion;
          }
        }
        /*
         * J_d_sink_mass
         */
        if (pd->v[pg->imtrx][POR_SINK_MASS]) {
          var = POR_SINK_MASS;
          pvar = upd->vp[pg->imtrx][var];
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            diffusion = 0.;
            if ((cr->MeshMotion == LAGRANGIAN || cr->MeshMotion == DYNAMIC_LAGRANGIAN) &&
                (pd->e[pg->imtrx][eqn] & T_DIFFUSION)) {
              phi_j = bf[var]->phi[j];
              for (p = 0; p < VIM; p++) {
                for (q = 0; q < VIM; q++) {
                  grad_phi_i_e_a[p][q] = bfm->grad_phi_e[i][a][p][q];

                  diffusion += grad_phi_i_e_a[p][q] * dTT_dsink_mass[q][p][j];
                }
              }
              diffusion *= -d_area;
              diffusion *= pd->etm[pg->imtrx][eqn][LOG2_DIFFUSION];
            }

            lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += diffusion;
          }
        }
        /*
         * J_d_T
         */
        if ((cr->MeshMotion == LAGRANGIAN || cr->MeshMotion == DYNAMIC_LAGRANGIAN) &&
            pd->e[pg->imtrx][eqn]) {
          var = TEMPERATURE;
          if (pd->v[pg->imtrx][var]) {
            pvar = upd->vp[pg->imtrx][var];

            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              diffusion = 0.;
              phi_j = bf[var]->phi[j];
              for (p = 0; p < VIM; p++) {
                for (q = 0; q < VIM; q++) {
                  diffusion += grad_phi_i_e_a[p][q] * dTT_dT[q][p][j];
                }
              }
              diffusion *= -d_area;
              diffusion *= diffusion_etm;

              /* add inertia of moving mesh */
              advection = 0.;
              if (advection_on) {
                /*
                for ( p=0; p<dim; p++)
                  {
                    for ( q=0; q<dim; q++)
                      {
                        grad_phi_i_e_a[p][q] =
                          bfm->grad_phi_e[i][a] [p][q];

                      }
                  }
                  */

                for (p = 0; p < dim; p++) {
                  for (q = 0; q < dim; q++) {

                    advect_a += grad_phi_i_e_a[p][q] *
                                (d_vconv->T[q][j] * vconv[p] + vconv[q] * d_vconv->T[p][j]);
                  }
                }
                advect_a *= rho * d_area;
                advect_a *= advection_etm;

                for (p = 0; p < dim; p++) {
                  for (q = 0; q < dim; q++) {
                    advect_b += grad_phi_i_e_a[p][q] * VV[q][p];
                  }
                }
                advect_b *= d_rho->T[j] * d_area;
                advect_b *= advection_etm;

                advection = advect_a + advect_b;
              }

              lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += advection + diffusion;
            }
          }
        }
        /*
         * J_d_max_strain
         */
        if ((cr->MeshMotion == LAGRANGIAN || cr->MeshMotion == DYNAMIC_LAGRANGIAN) &&
            pd->e[pg->imtrx][eqn]) {
          var = MAX_STRAIN;
          if (pd->v[pg->imtrx][var]) {
            pvar = upd->vp[pg->imtrx][var];

            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {

              diffusion = 0.0;
              for (q = 0; q < VIM; q++) {
                for (p = 0; p < VIM; p++) {
                  diffusion += grad_phi_i_e_a[p][q] * dTT_dmax_strain[q][p][j];
                }
              }
              diffusion *= -d_area;
              diffusion *= diffusion_etm;

              lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += diffusion;
            }
          }
        }
        /*
         * J_d_cur_strain
         */
        if ((cr->MeshMotion == LAGRANGIAN || cr->MeshMotion == DYNAMIC_LAGRANGIAN) &&
            pd->e[pg->imtrx][eqn]) {
          var = CUR_STRAIN;
          if (pd->v[pg->imtrx][var]) {
            pvar = upd->vp[pg->imtrx][var];

            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {

              diffusion = 0.0;
              for (q = 0; q < VIM; q++) {
                for (p = 0; p < VIM; p++) {
                  diffusion += grad_phi_i_e_a[p][q] * dTT_dcur_strain[q][p][j];
                }
              }
              diffusion *= -d_area;
              diffusion *= diffusion_etm;

              lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += diffusion;
            }
          }
        }

      } /* end of loop over equations i  */
    }   /* end of loop over equation directions a */
  }     /* end of if jacobian */

  return (status);
}

/*  _______________________________________________________________________  */

/* end of assemble_energy */

/*  _______________________________________________________________________  */

/******************************************************************************/
/* assemble_volume
 *
 * in:
 * 	ei -- pointer to Element Indeces		structure
 *	pd -- pointer to Problem Description		structure
 *	af -- pointer to Action Flag			structure
 *	bf -- pointer to Basis Function			structure
 *	fv -- pointer to Field Variable			structure
 *	cr -- pointer to Constitutive Relation		structure
 *	mp -- pointer to Material Property		structure
 *	esp-- pointer to Element Stiffness Pointers	structure
 *	lec-- pointer to Local Element Contribution	structure
 *
 * out:
 *	a   -- gets loaded up with proper contribution
 *	lec -- gets loaded up with local contributions to resid, Jacobian
 * 	r   -- residual RHS vector
 *
 * Created:	Amy Sun 1/22/99
 *
 * Note: currently we do a "double load" into the addresses in the global
 *       "a" matrix and resid vector held in "esp", as well as into the
 *       local accumulators in "lec".
 *
 */

int assemble_volume(bool owner) {
  int eqn;
  int var;
  int dim;
  int b;

  int j;
  int status;
  int iAC = 0;
  int vj;
  int VC_mode;
  int spec_id;
  int ktype, gnn, nvdof, ledof;

  dbl det_J;
  dbl h3; /* Volume element (scale factors). */
  dbl wt;
  dbl func;
  dbl phi_j;
  dbl d_det_J_dmeshbj;
  dbl dh3dmesh_bj; /* Sensitivity to (b,j) mesh dof. */
  dbl d_detJh3_dmeshbj;
  dbl d_func_dmeshbj; /* Sensitivity of func to (b,j)mesh dof. */

  double P_amb = 1.013E+06; /* cgs units */
  double R_gas = 82.07 * 1.013E+06;
  double moles_gas;
  status = 0;
  func = 1.;

  /*
   * Unpack variables from structures for local convenience...
   */
  eqn = pd->ShapeVar;
  dim = pd->Num_Dim;
  h3 = fv->h3; /* Differential volume element (scales). */
  det_J = bf[eqn]->detJ;

  iAC = pd->VolumeIntegral;

  VC_mode = augc[iAC].VOLID;
  spec_id = augc[iAC].COMPID;

  if (pd_glob[0]->CoordinateSystem != CYLINDRICAL && pd_glob[0]->CoordinateSystem != SWIRLING) {
    wt = fv->wt;
  } else {
    wt = 2.0 * M_PIE * fv->wt;
  }

  /*
   * Load up the function the volume multiplies to
   * 1= pure volume calc for the element
   * 2= rho * volume
   * 3= mass fraction of species w * volume
   * Add as needed
   */
  if (augc[iAC].len_AC > 0) {
    R_gas = augc[iAC].DataFlt[0];
    P_amb = augc[iAC].DataFlt[1];
  }

  switch (VC_mode) {
  case 1: /* pure volume calc */
  case 11:
    break;
  case 2: /* density*volume */
  case 12:
    func = mp->density;
    break;
  case 3: /* massfraction*volume */
  case 13:
    if (pd->Num_Species_Eqn < 1) {
      fprintf(stderr, "Need at least one species to use this option");
    } else {
      func = fv->c[spec_id];
    }
    break;
  case 4:
  case 14:
    if (pd->v[pg->imtrx][VELOCITY1]) {
      func = fv->v[0];
    } else {
      GOMA_EH(GOMA_ERROR, " must have momentum equation on for this AC");
    }
    break;
  case 5: /* ideal gas version */
  case 15:
    moles_gas = 0.0;
    for (j = 0; j < pd->Num_Species; j++) {
      moles_gas += fv->c[spec_id] / mp->molecular_weight[spec_id];
    }
    func = fv->P + P_amb - R_gas * fv->T * moles_gas;
    break;

  default:
    GOMA_EH(GOMA_ERROR, "assemble_volume() unknown integral calculation\n");
    break;
  }

  /* total integral */

  if (owner) {
    augc[iAC].evol += func * det_J * wt * h3;
  }

  /* Calculate sensitivities so that we don't need to compute it numerically
     in mm_sol_nonlinear.c  */

  for (b = 0; b < dim; b++) {
    var = MESH_DISPLACEMENT1 + b;
    if (pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        eqn = R_MESH1 + b;

        phi_j = bf[var]->phi[j];

        d_det_J_dmeshbj = bf[eqn]->d_det_J_dm[b][j];

        dh3dmesh_bj = fv->dh3dq[b] * phi_j;

        d_detJh3_dmeshbj = d_det_J_dmeshbj * h3 + det_J * dh3dmesh_bj;

        /*assume no dependency of func on meshbj, may change in the future*/

        d_func_dmeshbj = 0.;

        /*grab Index_Solution for location of x. see mm_fill_ptrs.c */

        vj = ei[pg->imtrx]->gun_list[var][j];

        augc[iAC].d_evol_dx[vj] += wt * (func * d_detJh3_dmeshbj + d_func_dmeshbj * det_J * h3);
      }
    }
  }

  if (VC_mode == 3) {
    var = MASS_FRACTION;
    eqn = R_MASS;
    if (pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        phi_j = bf[var]->phi[j];
        ktype = spec_id;
        gnn = ei[pg->imtrx]->gnn_list[eqn][j];
        nvdof = ei[pg->imtrx]->Baby_Dolphin[eqn][j];
        ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][j];
        vj = Index_Solution(gnn, eqn, ktype, nvdof, ei[pg->imtrx]->matID_ledof[ledof], pg->imtrx);
        augc[iAC].d_evol_dx[vj] += phi_j * det_J * wt * h3;
      }
    }
  }
  if (VC_mode == 4) {
    var = VELOCITY1;
    eqn = R_MOMENTUM1;
    if (pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        phi_j = bf[var]->phi[j];
        ktype = 0;
        gnn = ei[pg->imtrx]->gnn_list[eqn][j];
        nvdof = ei[pg->imtrx]->Baby_Dolphin[eqn][j];
        ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][j];
        vj = Index_Solution(gnn, eqn, ktype, nvdof, ei[pg->imtrx]->matID_ledof[ledof], pg->imtrx);
        augc[iAC].d_evol_dx[vj] += phi_j * det_J * wt * h3;
      }
    }
  }
  if (VC_mode == 5) {
    var = PRESSURE;
    eqn = R_PRESSURE;
    if (pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        phi_j = bf[var]->phi[j];
        ktype = 0;
        gnn = ei[pg->imtrx]->gnn_list[eqn][j];
        nvdof = ei[pg->imtrx]->Baby_Dolphin[eqn][j];
        ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][j];
        vj = Index_Solution(gnn, eqn, ktype, nvdof, ei[pg->imtrx]->matID_ledof[ledof], pg->imtrx);
        augc[iAC].d_evol_dx[vj] += phi_j * det_J * wt * h3;
      }
    }
    var = MASS_FRACTION;
    eqn = R_MASS;
    if (pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        phi_j = bf[var]->phi[j];
        for (ktype = 0; ktype < pd->Num_Species; ktype++) {
          gnn = ei[pg->imtrx]->gnn_list[eqn][j];
          nvdof = ei[pg->imtrx]->Baby_Dolphin[eqn][j];
          ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][j];
          vj = Index_Solution(gnn, eqn, ktype, nvdof, ei[pg->imtrx]->matID_ledof[ledof], pg->imtrx);
          augc[iAC].d_evol_dx[vj] -=
              R_gas * fv->T / mp->molecular_weight[ktype] * phi_j * det_J * wt * h3;
        }
      }
    }
    var = TEMPERATURE;
    eqn = R_ENERGY;
    if (pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        phi_j = bf[var]->phi[j];
        ktype = spec_id;
        gnn = ei[pg->imtrx]->gnn_list[eqn][j];
        nvdof = ei[pg->imtrx]->Baby_Dolphin[eqn][j];
        ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][j];
        vj = Index_Solution(gnn, eqn, ktype, nvdof, ei[pg->imtrx]->matID_ledof[ledof], pg->imtrx);
        augc[iAC].d_evol_dx[vj] -=
            R_gas * (fv->c[spec_id] / mp->molecular_weight[spec_id]) * phi_j * det_J * wt * h3;
      }
    }
  }

  return (status);
}

/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/
/* assemble_curvature--
 *
 *   Assemble the residual and jacobian contributions of the
 *   curvature of the level set function:
 *           H = div( grad_F/|grad_F|)
 *
 *
 *
 * in:
 * 	ei -- pointer to Element Indices	structure
 *	pd -- pointer to Problem Description	structure
 *	af -- pointer to Action Flag		structure
 *	bf -- pointer to Basis Function		structure
 *	fv -- pointer to Field Variable		structure
 *	cr -- pointer to Constitutive Relation	structure
 *	md -- pointer to Mesh Derivative	structure
 *	me -- pointer to Material Entity	structure
 *
 * out:
 *	lec -- gets loaded up with local contributions to resid, Jacobian
 *
 * Created:	Feb 12, 2002 tabaer@sandia.gov
 *
 * Revised:
 *
 *
 *
 */

int assemble_curvature(void) /*  time step size      */
{
  int i, j, p, a;
  int peqn, pvar;
  int var;

  int dim;
  int eqn;
  int status = 0;

  double det_J;
  double h3; /* Volume element (scale factors). */
  double wt_func;
  double wt;

  double source, diffusion = 0.0;

  dim = pd->Num_Dim;

  /*
   * Bail out fast if there's nothing to do...
   */

  if (!pd->e[pg->imtrx][eqn = R_CURVATURE]) {
    return (status);
  }

  /* only relevant if level set field exists */

  if (ls == NULL)
    GOMA_EH(GOMA_ERROR, " Curvature equation can only be used in conjunction with level set.");

  peqn = upd->ep[pg->imtrx][eqn];

  wt = fv->wt;

  det_J = bf[eqn]->detJ;

  h3 = fv->h3;

  load_lsi(ls->Length_Scale);
  load_lsi_derivs();

  /*
   * Residuals_________________________________________________________________
   */

  if (af->Assemble_Residual) {
    if (pd->e[pg->imtrx][eqn]) {
      for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {

        wt_func = bf[eqn]->phi[i];

        diffusion = 0.;
        /* Now granted this really isn't a diffusion term
         * but it does involve a divergence operator integrated by parts
         */
        if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {

          for (a = 0; a < dim; a++)
            diffusion += bf[eqn]->grad_phi[i][a] * lsi->normal[a];

          diffusion *= det_J * wt * h3;
          diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
        }

        source = 0.0;

        if (pd->e[pg->imtrx][eqn] & T_SOURCE) {

          source += fv->H;

          source *= wt_func * det_J * h3 * wt;
          source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
        }

        lec->R[LEC_R_INDEX(peqn, i)] += source + diffusion;
      }
    }
  }
  /*
   * Jacobian terms_________________________________________________________________
   */

  if (af->Assemble_Jacobian) {

    if (pd->e[pg->imtrx][eqn]) {
      for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
        wt_func = bf[eqn]->phi[i];

        /*
         * J_H_H
         */

        var = CURVATURE;

        if (pd->v[pg->imtrx][var]) {
          pvar = upd->vp[pg->imtrx][var];

          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            source = 0.0;

            if (pd->e[pg->imtrx][eqn] & T_SOURCE) {

              source += bf[var]->phi[j];
              source *= wt_func * det_J * wt * h3;
              source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
            }

            lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += source;
          }
        }
        /*
         * J_H_mesh
         */

        for (a = 0; a < dim; a++) {
          var = MESH_DISPLACEMENT1 + a;

          if (pd->v[pg->imtrx][var]) {
            pvar = upd->vp[pg->imtrx][var];

            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              source = 0.0;

              if (pd->e[pg->imtrx][eqn] & T_SOURCE) {

                source += h3 * bf[var]->d_det_J_dm[a][j] + det_J * fv->dh3dmesh[a][j];
                source *= wt_func * fv->H * wt;
                source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
              }

              diffusion = 0.0;

              if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {
                double diff1 = 0.0, diff2 = 0.0, diff3 = 0.0;

                for (p = 0; p < dim; p++)
                  diff1 += bf[eqn]->d_grad_phi_dmesh[i][p][a][j] * lsi->normal[p];

                diff1 *= wt * h3 * det_J;

                for (p = 0; p < dim; p++)
                  diff2 += bf[eqn]->grad_phi[i][p] * lsi->d_normal_dmesh[p][a][j];

                diff2 *= wt * h3 * det_J;

                for (p = 0; p < dim; p++)
                  diff3 += bf[eqn]->grad_phi[p][i] * lsi->normal[p];

                diff3 *= h3 * bf[var]->d_det_J_dm[a][j] + det_J * fv->dh3dmesh[a][j];
                diff3 *= wt;

                diffusion += diff1 + diff2 + diff3;

                diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
              }

              lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += source + diffusion;
            }
          }
        }
        /*
         * J_H_F
         */

        var = LS;

        if (pd->v[pg->imtrx][var]) {

          pvar = upd->vp[pg->imtrx][var];

          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {

            if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {

              diffusion = 0.0;

              for (p = 0; p < dim; p++)
                diffusion += bf[eqn]->grad_phi[i][p] * lsi->d_normal_dF[p][j];

              diffusion *= det_J * wt * h3;
              diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
            }

            lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += diffusion;
          }
        }
      }
    }
  }
  return (0);
}

/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/
/* assemble_div_normals --
 *
 *   Assemble the residual and jacobian contributions of the
 *   curvature of the level set function via divegence of the level set normal vector field,
 *videlicit:
 *
 *           H = div( n )
 *
 *
 *
 * in:
 * 	ei -- pointer to Element Indices	structure
 *	pd -- pointer to Problem Description	structure
 *	af -- pointer to Action Flag		structure
 *	bf -- pointer to Basis Function		structure
 *	fv -- pointer to Field Variable		structure
 *	cr -- pointer to Constitutive Relation	structure
 *	md -- pointer to Mesh Derivative	structure
 *	me -- pointer to Material Entity	structure
 *
 * out:
 *	lec -- gets loaded up with local contributions to resid, Jacobian
 *
 * Created:	Aug 27, 2003 tabaer@sandia.gov
 *
 * Revised:
 *
 *
 *
 */

int assemble_div_normals(void) /*  time step size      */
{
  int i, j, p, a;
  int peqn, pvar;
  int var;

  int dim;
  int eqn;
  int status = 0;

  double det_J;
  double h3; /* Volume element (scale factors). */
  double wt_func;
  double wt;

  double source, diffusion;

  dim = pd->Num_Dim;

  /*
   * Bail out fast if there's nothing to do...
   */

  if (!pd->e[pg->imtrx][eqn = R_CURVATURE]) {
    return (status);
  }

  /* only relevant if level set field exists */

  if (ls == NULL)
    GOMA_EH(GOMA_ERROR, " Curvature equation can only be used in conjunction with level set.");

  peqn = upd->ep[pg->imtrx][eqn];

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
   * Residuals_________________________________________________________________
   */

  if (af->Assemble_Residual) {
    if (pd->e[pg->imtrx][eqn]) {
      for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
        wt_func = bf[eqn]->phi[i];

        diffusion = 0.;

        /* Now granted this really isn't a diffusion term
         * but it does involve a divergence operator
         */
        if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {
          diffusion = fv->div_n;
          diffusion *= wt_func * det_J * wt * h3;
          diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
        }

        source = 0.0;
        if (pd->e[pg->imtrx][eqn] & T_SOURCE) {

          source += fv->H;
          source *= wt_func * det_J * h3 * wt;
          source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
        }

        lec->R[LEC_R_INDEX(peqn, i)] += source + diffusion;
      }
    }
  }

  /*
   * Jacobian terms_________________________________________________________________
   */

  if (af->Assemble_Jacobian) {
    if (pd->e[pg->imtrx][eqn]) {
      for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
        wt_func = bf[eqn]->phi[i];

        /*
         * J_H_H
         */

        var = CURVATURE;

        if (pd->v[pg->imtrx][var]) {
          pvar = upd->vp[pg->imtrx][var];

          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            source = 0.0;

            if (pd->e[pg->imtrx][eqn] & T_SOURCE) {

              source += bf[var]->phi[j];
              source *= wt_func * det_J * wt * h3;
              source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
            }

            lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += source;
          }
        }

        /*
         * J_H_mesh
         */

        /*
         * J_H_n
         */

        for (a = 0; a < dim; a++) {
          var = NORMAL1 + a;

          if (pd->v[pg->imtrx][var]) {
            pvar = upd->vp[pg->imtrx][var];

            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              diffusion = 0.0;

              if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {

                for (p = 0; p < dim; p++) {
                  diffusion += bf[var]->grad_phi_e[j][a][p][p];
                }

                diffusion *= wt_func * wt * det_J * h3;
                diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
              }

              lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += diffusion;
            }
          }
        }
      }
    }
  }
  return (0);
}

/******************************************************************************/
/* assemble_LSvelocity
 *
 * in:
 * 	ei -- pointer to Element Indeces		structure
 *	pd -- pointer to Problem Description		structure
 *	af -- pointer to Action Flag			structure
 *	bf -- pointer to Basis Function			structure
 *	fv -- pointer to Field Variable			structure
 *	cr -- pointer to Constitutive Relation		structure
 *	mp -- pointer to Material Property		structure
 *	esp-- pointer to Element Stiffness Pointers	structure
 *
 * out:
 *	augc.lsvel      -- contribution to velocity of phase SIGN
 *	augc.d_lsvel_dx -- contribution to derivative of velocity wrt. v, H & x
 *	augc.d_lsvol_dx -- contribution to derivative of volume wrt. v, H & x
 *
 * Created:	Anne Grillet  2/4/2002
 *
 *
 */

int assemble_LSvelocity(bool owner, int ielem) {
  int eqn;
  int var;
  int dim;
  int vel_dim;
  int b;

  int j;
  int status = 0;
  int iAC = 0;
  int vj;

  dbl det_J;
  dbl h3; /* Volume element (scale factors). */
  dbl wt;
  dbl phi_j;
  dbl d_det_J_dmeshbj;
  dbl dh3dmesh_bj; /* Sensitivity to (b,j) mesh dof. */
  dbl d_detJh3_dmeshbj;
  dbl d_H_dmeshbj; /* Sensitivity of H to (b,j) mesh dof. */
  dbl d_H_dFj;     /* Sensitivity of H to j fill function dof. */

  dbl alpha; /* Level set length scale */
  dbl H;     /* Level set heaviside function  */

  /*
   * Unpack variables from structures for local convenience...
   */
  eqn = pd->ShapeVar;
  dim = pd->Num_Dim;
  h3 = fv->h3; /* Differential volume element (scales). */
  det_J = bf[eqn]->detJ;

  iAC = pd->LSVelocityIntegral;

  alpha = ls->Length_Scale;

  /* Added 3/15/01 */
  if (pd_glob[0]->CoordinateSystem != CYLINDRICAL && pd_glob[0]->CoordinateSystem != SWIRLING) {
    wt = fv->wt;
  } else {
    wt = 2.0 * M_PIE * fv->wt;
  }

  vel_dim = (augc[iAC].DIR == I_NEG_VX
                 ? 0
                 : (augc[iAC].DIR == I_NEG_VY ? 1 : (augc[iAC].DIR == I_NEG_VZ ? 2 : -1)));

  if (vel_dim == -1) {
    vel_dim = (augc[iAC].DIR == I_POS_VX
                   ? 0
                   : (augc[iAC].DIR == I_POS_VY ? 1 : (augc[iAC].DIR == I_POS_VZ ? 2 : -1)));
  }

  if (vel_dim == -1)
    GOMA_EH(GOMA_ERROR, "assemble_LSVelocity(): Couldn't identify velocity component and phase.\n");

  load_lsi(alpha);
  load_lsi_derivs();

  H = augc[iAC].LSPHASE == I_POS_FILL ? lsi->H : (1.0 - lsi->H);

  /* total integral */

  if (owner) {
    augc[iAC].lsvel += wt * det_J * H * h3 * fv->v[vel_dim];
    augc[iAC].lsvol += wt * det_J * H * h3;
  }

  /* Calculate sensitivities so that we don't need to compute it numerically
     in mm_sol_nonlinear.c  */

  /* sensitivity to velocity                                          */

  for (b = 0; b < dim; b++) {
    var = vel_dim;
    if (pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        vj = ei[pg->imtrx]->gun_list[var][j];

        augc[iAC].d_lsvel_dx[vj] += wt * det_J * h3 * H * bf[var]->phi[j];
        augc[iAC].d_lsvol_dx[vj] += 0.0;
      }
    }
  }

  /* sensitivity to mesh                                  */

  for (b = 0; b < dim; b++) {
    var = MESH_DISPLACEMENT1 + b;
    if (pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {

        eqn = R_MESH1 + b;

        phi_j = bf[var]->phi[j];

        d_det_J_dmeshbj = bf[eqn]->d_det_J_dm[b][j];

        dh3dmesh_bj = fv->dh3dq[b] * phi_j;

        d_detJh3_dmeshbj = d_det_J_dmeshbj * h3 + det_J * dh3dmesh_bj;

        /* now, H doesn't depend on mesh (load_lsi_derivs), but leave this in
           just in case that changes  */

        d_H_dmeshbj = lsi->d_H_dmesh[b][j];

        /*grab Index_Solution for location of x. see mm_fill_ptrs.c */

        vj = ei[pg->imtrx]->gun_list[var][j];

        augc[iAC].d_lsvel_dx[vj] +=
            wt * fv->v[vel_dim] * (H * d_detJh3_dmeshbj + d_H_dmeshbj * det_J * h3);
        augc[iAC].d_lsvol_dx[vj] += wt * (H * d_detJh3_dmeshbj + d_H_dmeshbj * det_J * h3);
      }
    }
  }

  /* sensitivity to FILL                                  */

  var = FILL;
  if (pd->v[pg->imtrx][var]) {
    for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {

      d_H_dFj = lsi->d_H_dF[j];

      /*grab Index_Solution for location of x. see mm_fill_ptrs.c */

      vj = ei[pg->imtrx]->gun_list[var][j];

      augc[iAC].d_lsvel_dx[vj] += wt * fv->v[vel_dim] * d_H_dFj * det_J * h3;
      augc[iAC].d_lsvol_dx[vj] += wt * d_H_dFj * det_J * h3;
    }
  }

  return (status);
}

/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/
/* assemble_normals--
 *
 *   Assemble the residual and jacobian contributions of the
 *     of the level set normal vector projection equation:
 *
 *           n = grad(phi);
 *
 *
 *
 * in:
 * 	ei -- pointer to Element Indices	structure
 *	pd -- pointer to Problem Description	structure
 *	af -- pointer to Action Flag		structure
 *	bf -- pointer to Basis Function		structure
 *	fv -- pointer to Field Variable		structure
 *	cr -- pointer to Constitutive Relation	structure
 *	md -- pointer to Mesh Derivative	structure
 *	me -- pointer to Material Entity	structure
 *
 * out:
 *	lec -- gets loaded up with local contributions to resid, Jacobian
 *
 * Created:	Aug 24, 2003 tabaer@sandia.gov
 *
 * Revised:
 *
 *
 *
 */
int assemble_normals(void) {

  int i, j, p, a, b;
  int ii;
  int peqn, pvar;
  int var;

  int dim;
  int eqn;
  int status = 0;
  struct Basis_Functions *bfm;

  double det_J;
  double h3; /* Volume element (scale factors). */
  double wt_func;
  double phi_i, phi_j;
  double wt;

  double normal[DIM];
  double P[DIM][DIM];

  double mag_grad_F;

  double advection = 0.0, advection1 = 0.0, source = 0.0;

  dim = pd->Num_Dim;

  /*
   * Bail out fast if there's nothing to do...
   */

  if (!pd->e[pg->imtrx][eqn = R_NORMAL1]) {
    return (status);
  }

  /* only relevant if level set field exists */

  if (ls == NULL) {
    GOMA_WH(-1, " Normals to level set can only be used in conjunction with level set.");
    return (status);
  }

  peqn = upd->ep[pg->imtrx][eqn];

  wt = fv->wt;

  det_J = bf[eqn]->detJ;

  h3 = fv->h3;

  for (a = 0; a < dim; a++)
    normal[a] = fv->grad_F[a];

  mag_grad_F = normalize_really_simple_vector(normal, dim);

  for (a = 0; a < dim; a++) {
    for (b = 0; b < dim; b++) {
      P[a][b] = (double)delta(a, b) - normal[a] * normal[b];
    }
  }

  /*
   * Residuals_________________________________________________________________
   */

  if (af->Assemble_Residual) {
    for (a = 0; a < dim; a++) {
      eqn = R_NORMAL1 + a;
      peqn = upd->ep[pg->imtrx][eqn];
      bfm = bf[eqn];

      if (pd->e[pg->imtrx][eqn]) {
        for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {

          ii = ei[pg->imtrx]->lvdof_to_row_lvdof[eqn][i];

          wt_func = bfm->phi[i];

          source = 0.0;

          /*
           * Not really a source per se, but trying to follow the GOMA paradigm, such as it is.
           */

          if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
            source = fv->n[a];
            source *= wt * wt_func * h3 * det_J;
            source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
          }

          /*
           * Ditto
           */

          advection = 0.0;

          if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
            advection = -normal[a];
            advection *= wt * wt_func * h3 * det_J;
            advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
          }

          lec->R[LEC_R_INDEX(peqn, ii)] += source + advection;
        }
      }
    }
  }
  /*
   * Jacobian terms...
   */

  if (af->Assemble_Jacobian) {
    for (a = 0; a < dim; a++) {

      int ledof;

      eqn = R_NORMAL1 + a;
      peqn = upd->ep[pg->imtrx][eqn];
      bfm = bf[eqn];

      for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
        ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][i];
        if (ei[pg->imtrx]->active_interp_ledof[ledof]) {
          ii = ei[pg->imtrx]->lvdof_to_row_lvdof[eqn][i];

          phi_i = bfm->phi[i];

          /*
           *  J_n_n
           */

          for (b = 0; b < dim; b++) {
            var = NORMAL1 + b;
            pvar = upd->vp[pg->imtrx][var];

            if (eqn == var) {
              for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                phi_j = bf[var]->phi[j];

                if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
                  source = 0.0;
                  source = phi_j;

                  source *= phi_i * wt * det_J * h3;
                  source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
                }

                lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
              }
            }
          }

          /*
           *  J_n_F
           */
          var = LS;
          pvar = upd->vp[pg->imtrx][var];

          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];

            advection = 0.0;

            if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
              for (b = 0; b < dim; b++) {
                advection += P[a][b] * bf[var]->grad_phi[j][b] / mag_grad_F;
              }

              advection *= -phi_i * wt * h3 * det_J;
              advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
            }
            lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += advection;
          }

          /*
           * J_n_d
           */

          for (b = 0; b < dim; b++) {
            var = R_MESH1 + b;
            pvar = upd->vp[pg->imtrx][var];

            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              advection = 0.0;

              if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {

                for (p = 0; p < dim; p++) {
                  advection += P[a][p] * fv->d_grad_F_dmesh[p][b][j] / mag_grad_F;
                }

                advection *= -phi_i * wt * det_J * h3;
                advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];

                advection1 =
                    -normal[a] * (bf[var]->d_det_J_dm[b][j] * h3 + fv->dh3dmesh[b][j] * det_J);
                advection1 *= phi_i * wt;
                advection1 *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
              }
              lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += advection + advection1;
            }
          }
        }
      }
    }
  }
  return (status);
}

int assemble_ls_momentum_source(void) {

  int i, j, p, q, a;
  int ii;
  int peqn, pvar;
  int var;

  int dim;
  int eqn;
  int status = 0;
  struct Basis_Functions *bfm;
  double(*grad_phi_i_e_a)[DIM] = NULL;

  double det_J;
  double h3; /* Volume element (scale factors). */
  double phi_i, phi_j;
  double wt;

  double csf[DIM][DIM];
  double d_csf_dF[DIM][DIM][MDE];

  double source = 0.0;

  dim = pd->Num_Dim;

  /*
   * Bail out fast if there's nothing to do...
   */

  if (!pd->e[pg->imtrx][eqn = R_LAGR_MULT1]) {
    return (status);
  }

  /* only relevant if level set field exists */

  if (ls == NULL) {
    GOMA_WH(-1, " Normals to level set can only be used in conjunction with level set.");
    return (status);
  }

  peqn = upd->ep[pg->imtrx][eqn];

  wt = fv->wt;

  det_J = bf[eqn]->detJ;

  h3 = fv->h3;

  memset(csf, 0, sizeof(double) * DIM * DIM);
  memset(d_csf_dF, 0, sizeof(double) * DIM * DIM * MDE);

  /* Fetch the level set interface functions. */
  load_lsi(ls->Length_Scale);

  /* Calculate the CSF tensor. */
  for (p = 0; p < VIM; p++) {
    for (q = 0; q < VIM; q++) {
      csf[p][q] = mp->surface_tension * ((double)delta(p, q) - lsi->normal[p] * lsi->normal[q]);
    }
  }

  load_lsi_derivs();

  /* Calculate the derivatives. */
  var = FILL;
  for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
    for (p = 0; p < VIM; p++) {
      for (q = 0; q < VIM; q++) {
        d_csf_dF[p][q][j] = -mp->surface_tension * lsi->d_normal_dF[p][j] * lsi->normal[q] -
                            mp->surface_tension * lsi->normal[p] * lsi->d_normal_dF[q][j];
      }
    }
  }

  /*
   * Residuals_________________________________________________________________
   */

  if (af->Assemble_Residual) {
    for (a = 0; a < dim; a++) {
      eqn = R_LAGR_MULT1 + a;
      peqn = upd->ep[pg->imtrx][eqn];
      bfm = bf[eqn];

      if (pd->e[pg->imtrx][eqn]) {
        for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {

          ii = ei[pg->imtrx]->lvdof_to_row_lvdof[eqn][i];

          phi_i = bfm->phi[i];
          grad_phi_i_e_a = bfm->grad_phi_e[i][a];

          if (pd->e[pg->imtrx][eqn]) {
            source = fv->lm[a] * phi_i;

            for (p = 0; p < VIM; p++) {
              for (q = 0; q < VIM; q++) {
                source += grad_phi_i_e_a[p][q] * csf[q][p];
              }
            }
            source *= wt * h3 * det_J;

            lec->R[LEC_R_INDEX(peqn, ii)] += source;
          }
        }
      }
    }
  }
  /*
   * Jacobian terms...
   */

  if (af->Assemble_Jacobian) {
    for (a = 0; a < dim; a++) {

      int ledof;

      eqn = R_LAGR_MULT1 + a;
      peqn = upd->ep[pg->imtrx][eqn];
      bfm = bf[eqn];

      for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
        ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][i];
        if (ei[pg->imtrx]->active_interp_ledof[ledof]) {
          ii = ei[pg->imtrx]->lvdof_to_row_lvdof[eqn][i];

          phi_i = bfm->phi[i];
          grad_phi_i_e_a = bfm->grad_phi_e[i][a];

          /*
           *  J_n_n
           */
          var = eqn;
          pvar = upd->vp[pg->imtrx][var];

          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];

            if (pd->e[pg->imtrx][eqn]) {
              source = phi_j * phi_i;

              source *= wt * det_J * h3;
            }

            lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
          }

          /*
           *  J_n_F
           */
          var = LS;
          pvar = upd->vp[pg->imtrx][var];

          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];

            source = 0.0;

            if (pd->e[pg->imtrx][eqn]) {
              for (p = 0; p < VIM; p++) {
                for (q = 0; q < VIM; q++) {
                  source += grad_phi_i_e_a[p][q] * d_csf_dF[q][p][j];
                }
              }
              source *= wt * det_J * h3;
            }
            lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
          }
        }
      }
    }
  }

  return (status);
}

int apply_ls_momentum_source(void) {
  int i, j, a, ii, ledof;
  int eqn, peqn, var, pvar;

  struct Basis_Functions *bfm;

  double wt, det_J, h3, phi_i, phi_j;

  double source;

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

          source = fv->lm[a] * lsi->delta;

          source *= phi_i * det_J * wt * h3;

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
   * Wesiduals ________________________________________________________________________________
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

          source = fv->lm[a] * lsi->delta;

          source *= phi_i * det_J * wt * h3;

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

          /* J_m_n
           */
          var = LAGR_MULT1 + a;
          if (pd->v[pg->imtrx][var]) {
            pvar = upd->vp[pg->imtrx][var];

            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              phi_j = bf[var]->phi[j];

              source = phi_j * lsi->delta;

              source *= phi_i * det_J * wt * h3;

              lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
            }
          }
          /* J_m_F
           */
          var = LS;
          if (pd->v[pg->imtrx][var]) {
            pvar = upd->vp[pg->imtrx][var];

            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              phi_j = bf[var]->phi[j];

              source = fv->lm[a] * lsi->d_delta_dF[j];

              source *= phi_i * det_J * wt * h3;

              lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
            }
          }
        }
      }
    }
  }

  return (1);
}

/********************************************************************************/

double acoustic_impedance(CONDUCTIVITY_DEPENDENCE_STRUCT *d_R, dbl time)

/**************************************************************************
 *
 * acoustic impedance
 *
 *   Calculate the acoustic impedance and its derivatives wrt to nodal dofs
 *   at the local gauss point
 *
 * Output
 * -----
 *    dbl d_R->T[j]    -> derivative of conductivity wrt the jth
 *                          Temperature unknown in an element
 *    dbl d_R->C[k][j] -> derivative of conductivity wrt the jth
 *                          species unknown of ktype, k, in the element.
 *    dbl d_R->X[a][j] -> derivative of conductivity wrt the jth
 *                          mesh displacement in the ath direction.
 *    dbl d_R->F[j]    -> derivative of conductivity wrt the jth
 *                          FILL unknown in an element
 *
 *  Return
 * --------
 *    Actual value of the acoustic impedance
 ***************************************************************************/

{
  int dim = pd->Num_Dim;
  int var, i, j, a, w, var_offset;
  double R;
  struct Level_Set_Data *ls_old;
  R = 0.;

  if (d_R != NULL) {
    memset(d_R->T, 0, sizeof(double) * MDE);
    memset(d_R->X, 0, sizeof(double) * DIM * MDE);
    memset(d_R->C, 0, sizeof(double) * MAX_CONC * MDE);
    memset(d_R->F, 0, sizeof(double) * MDE);
  }

  if (mp->Acoustic_ImpedanceModel == USER) {

    GOMA_EH(GOMA_ERROR, "user acoustic impedance code not ready yet");
    /*
      usr_acoustic_impedance(mp->u_acoustic_impedance,time);
    */

    R = mp->acoustic_impedance;

    var = TEMPERATURE;
    if (pd->v[pg->imtrx][TEMPERATURE] && d_R != NULL) {

      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_R->T[j] = mp->d_acoustic_impedance[var] * bf[var]->phi[j];
      }
    }

    if (pd->v[pg->imtrx][MESH_DISPLACEMENT1] && d_R != NULL) {
      for (a = 0; a < dim; a++) {
        var = MESH_DISPLACEMENT1 + a;
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_R->X[a][j] = mp->d_acoustic_impedance[var] * bf[var]->phi[j];
        }
      }
    }

    if (pd->v[pg->imtrx][MASS_FRACTION] && d_R != NULL) {
      for (w = 0; w < pd->Num_Species_Eqn; w++) {
        var = MASS_FRACTION;
        var_offset = MAX_VARIABLE_TYPES + w;
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_R->C[w][j] = mp->d_acoustic_impedance[var_offset] * bf[var]->phi[j];
        }
      }
    }
  } else if (mp->Acoustic_ImpedanceModel == CONSTANT) {
    R = mp->acoustic_impedance;
  } else if (mp->Acoustic_ImpedanceModel == TABLE) {
    struct Data_Table *table_local;
    table_local = MP_Tables[mp->acoustic_impedance_tableid];
    apply_table_mp(&mp->acoustic_impedance, table_local);
    R = mp->acoustic_impedance;

    if (d_R != NULL) {
      for (i = 0; i < table_local->columns - 1; i++) {
        var = table_local->t_index[i];
        /* currently only set up to vary w.r.t. temperature */
        switch (var) {
        case TEMPERATURE:
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            d_R->T[j] = table_local->slope[i] * bf[var]->phi[j];
          }
          break;
        default:
          GOMA_EH(GOMA_ERROR, "Variable function not yet implemented in material property table");
        }
      }
    }
  } else if (mp->Acoustic_ImpedanceModel == LEVEL_SET) {
    ls_transport_property(mp->u_acoustic_impedance[0], mp->u_acoustic_impedance[1],
                          mp->u_acoustic_impedance[2], &mp->acoustic_impedance,
                          &mp->d_acoustic_impedance[FILL]);

    R = mp->acoustic_impedance;

    var = FILL;
    if (pd->v[pg->imtrx][var] && d_R != NULL) {

      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_R->F[j] = mp->d_acoustic_impedance[var] * bf[var]->phi[j];
      }
    }

  } else {
    GOMA_EH(GOMA_ERROR, "Unrecognized acoustic impedance model");
  }

  if (ls != NULL && mp->Acoustic_ImpedanceModel != LEVEL_SET &&
      mp->Acoustic_ImpedanceModel != LS_QUADRATIC && mp->mp2nd != NULL &&
      mp->mp2nd->AcousticImpedanceModel ==
          CONSTANT) /* Only Newtonian constitutive equation allowed for 2nd phase */
  {
    /* kludge for solidification tracking with phase function 0 */
    if (pfd != NULL && pd->e[pg->imtrx][R_EXT_VELOCITY]) {
      ls_old = ls;
      ls = pfd->ls[0];
      R = ls_modulate_thermalconductivity(R, mp->mp2nd->acousticimpedance_phase[0],
                                          ls->Length_Scale,
                                          (double)mp->mp2nd->acousticimpedancemask[0],
                                          (double)mp->mp2nd->acousticimpedancemask[1], d_R);
      ls = ls_old;
    }
    R = ls_modulate_thermalconductivity(R, mp->mp2nd->acousticimpedance, ls->Length_Scale,
                                        (double)mp->mp2nd->acousticimpedancemask[0],
                                        (double)mp->mp2nd->acousticimpedancemask[1], d_R);
  }

  return (R);
}

double wave_number(CONDUCTIVITY_DEPENDENCE_STRUCT *d_k, dbl time) {
  int var, j, i;
  double k;
  struct Level_Set_Data *ls_old;

  k = 0.;

  if (d_k != NULL) {
    memset(d_k->T, 0, sizeof(double) * MDE);
    memset(d_k->X, 0, sizeof(double) * DIM * MDE);
    memset(d_k->C, 0, sizeof(double) * MAX_CONC * MDE);
    memset(d_k->F, 0, sizeof(double) * MDE);
  }

  if (mp->wave_numberModel == CONSTANT) {
    k = mp->wave_number;
  } else if (mp->wave_numberModel == TABLE) {
    struct Data_Table *table_local;
    table_local = MP_Tables[mp->wave_number_tableid];
    apply_table_mp(&mp->wave_number, table_local);
    k = mp->wave_number;

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
  } else if (mp->wave_numberModel == LEVEL_SET) {
    ls_transport_property(mp->u_wave_number[0], mp->u_wave_number[1], mp->u_wave_number[2],
                          &mp->wave_number, &mp->d_wave_number[FILL]);

    k = mp->wave_number;

    var = FILL;
    if (pd->v[pg->imtrx][var] && d_k != NULL) {

      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_k->F[j] = mp->d_wave_number[var] * bf[var]->phi[j];
      }
    }

  } else {
    GOMA_EH(GOMA_ERROR, "Unrecognized acoustic wave number model");
  }

  if (ls != NULL && mp->wave_numberModel != LEVEL_SET && mp->wave_numberModel != LS_QUADRATIC &&
      mp->mp2nd != NULL &&
      mp->mp2nd->wavenumberModel ==
          CONSTANT) /* Only Newtonian constitutive equation allowed for 2nd phase */
  {
    /* kludge for solidification tracking with phase function 0 */
    if (pfd != NULL && pd->e[pg->imtrx][R_EXT_VELOCITY]) {
      ls_old = ls;
      ls = pfd->ls[0];
      k = ls_modulate_thermalconductivity(k, mp->mp2nd->wavenumber_phase[0], ls->Length_Scale,
                                          (double)mp->mp2nd->wavenumbermask[0],
                                          (double)mp->mp2nd->wavenumbermask[1], d_k);
      ls = ls_old;
    }
    k = ls_modulate_thermalconductivity(k, mp->mp2nd->wavenumber, ls->Length_Scale,
                                        (double)mp->mp2nd->wavenumbermask[0],
                                        (double)mp->mp2nd->wavenumbermask[1], d_k);
  }
  return (k);
}

/*********************************************************************************/
double acoustic_absorption(CONDUCTIVITY_DEPENDENCE_STRUCT *d_alpha, dbl time)

/**************************************************************************
 *
 * acoustic absorption
 *
 *   Calculate the acoustic absorption and its derivatives wrt to nodal dofs
 *   at the local gauss point
 *
 * Output
 * -----
 *    dbl d_R->T[j]    -> derivative of conductivity wrt the jth
 *                          Temperature unknown in an element
 *    dbl d_R->C[k][j] -> derivative of conductivity wrt the jth
 *                          species unknown of ktype, k, in the element.
 *    dbl d_R->X[a][j] -> derivative of conductivity wrt the jth
 *                          mesh displacement in the ath direction.
 *    dbl d_R->F[j]    -> derivative of conductivity wrt the jth
 *                          FILL unknown in an element
 *
 *  Return
 * --------
 *    Actual value of the acoustic absorption
 ***************************************************************************/

{
  int dim = pd->Num_Dim;
  int var, i, j, a, w, var_offset;
  double alpha;
  struct Level_Set_Data *ls_old;
  alpha = 0.;

  if (d_alpha != NULL) {
    memset(d_alpha->T, 0, sizeof(double) * MDE);
    memset(d_alpha->X, 0, sizeof(double) * DIM * MDE);
    memset(d_alpha->C, 0, sizeof(double) * MAX_CONC * MDE);
    memset(d_alpha->F, 0, sizeof(double) * MDE);
  }

  if (mp->Acoustic_AbsorptionModel == USER) {

    mp->acoustic_absorption = mp->u_acoustic_absorption[0];
    mp->d_acoustic_absorption[TEMPERATURE] = 0;
    for (a = 0; a < DIM; a++) {
      mp->d_acoustic_absorption[MESH_DISPLACEMENT1 + a] = 0;
    }
    for (a = 0; a < pd->Num_Species_Eqn; a++) {
      mp->d_acoustic_absorption[MAX_VARIABLE_TYPES + a] = mp->u_acoustic_absorption[a + 1];
      mp->acoustic_absorption += mp->u_acoustic_absorption[a + 1] * fv->c[a];
    }
    alpha = mp->acoustic_absorption;

    var = TEMPERATURE;
    if (pd->v[pg->imtrx][TEMPERATURE] && d_alpha != NULL) {

      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_alpha->T[j] = mp->d_acoustic_absorption[var] * bf[var]->phi[j];
      }
    }

    if (pd->v[pg->imtrx][MESH_DISPLACEMENT1] && d_alpha != NULL) {
      for (a = 0; a < dim; a++) {
        var = MESH_DISPLACEMENT1 + a;
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_alpha->X[a][j] = mp->d_acoustic_absorption[var] * bf[var]->phi[j];
        }
      }
    }

    if (pd->v[pg->imtrx][MASS_FRACTION] && d_alpha != NULL) {
      for (w = 0; w < pd->Num_Species_Eqn; w++) {
        var = MASS_FRACTION;
        var_offset = MAX_VARIABLE_TYPES + w;
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_alpha->C[w][j] = mp->d_acoustic_absorption[var_offset] * bf[var]->phi[j];
        }
      }
    }
  } else if (mp->Acoustic_AbsorptionModel == CONSTANT) {
    alpha = mp->acoustic_absorption;
  } else if (mp->Acoustic_AbsorptionModel == TABLE) {
    struct Data_Table *table_local;
    table_local = MP_Tables[mp->acoustic_absorption_tableid];
    apply_table_mp(&mp->acoustic_absorption, table_local);
    alpha = mp->acoustic_absorption;

    if (d_alpha != NULL) {
      for (i = 0; i < table_local->columns - 1; i++) {
        var = table_local->t_index[i];
        /* currently only set up to vary w.r.t. temperature */
        switch (var) {
        case TEMPERATURE:
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            d_alpha->T[j] = table_local->slope[i] * bf[var]->phi[j];
          }
          break;
        default:
          GOMA_EH(GOMA_ERROR, "Variable function not yet implemented in material property table");
        }
      }
    }
  } else if (mp->Acoustic_AbsorptionModel == LEVEL_SET) {
    ls_transport_property(mp->u_acoustic_absorption[0], mp->u_acoustic_absorption[1],
                          mp->u_acoustic_absorption[2], &mp->acoustic_absorption,
                          &mp->d_acoustic_absorption[FILL]);

    alpha = mp->acoustic_absorption;

    var = FILL;
    if (pd->v[pg->imtrx][var] && d_alpha != NULL) {

      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_alpha->F[j] = mp->d_acoustic_absorption[var] * bf[var]->phi[j];
      }
    }

  } else {
    GOMA_EH(GOMA_ERROR, "Unrecognized acoustic absorption model");
  }

  if (ls != NULL && mp->Acoustic_AbsorptionModel != LEVEL_SET &&
      mp->Acoustic_AbsorptionModel != LS_QUADRATIC && mp->mp2nd != NULL &&
      mp->mp2nd->AcousticAbsorptionModel ==
          CONSTANT) /* Only Newtonian constitutive equation allowed for 2nd phase */
  {
    /* kludge for solidification tracking with phase function 0 */
    if (pfd != NULL && pd->e[pg->imtrx][R_EXT_VELOCITY]) {
      ls_old = ls;
      ls = pfd->ls[0];
      alpha = ls_modulate_thermalconductivity(
          alpha, mp->mp2nd->acousticabsorption_phase[0], ls->Length_Scale,
          (double)mp->mp2nd->acousticabsorptionmask[0],
          (double)mp->mp2nd->acousticabsorptionmask[1], d_alpha);
      ls = ls_old;
    }
    alpha = ls_modulate_thermalconductivity(alpha, mp->mp2nd->acousticabsorption, ls->Length_Scale,
                                            (double)mp->mp2nd->acousticabsorptionmask[0],
                                            (double)mp->mp2nd->acousticabsorptionmask[1], d_alpha);
  }

  return (alpha);
}

/*********************************************************************************/
double light_absorption(CONDUCTIVITY_DEPENDENCE_STRUCT *d_alpha, dbl time)

/**************************************************************************
 *
 * light absorption
 *
 *   Calculate the light absorption and its derivatives wrt to nodal dofs
 *   at the local gauss point
 *
 * Output
 * -----
 *    dbl d_R->T[j]    -> derivative of conductivity wrt the jth
 *                          Temperature unknown in an element
 *    dbl d_R->C[k][j] -> derivative of conductivity wrt the jth
 *                          species unknown of ktype, k, in the element.
 *    dbl d_R->X[a][j] -> derivative of conductivity wrt the jth
 *                          mesh displacement in the ath direction.
 *    dbl d_R->F[j]    -> derivative of conductivity wrt the jth
 *                          FILL unknown in an element
 *
 *  Return
 * --------
 *    Actual value of the light absorption
 ***************************************************************************/

{
  int dim = pd->Num_Dim;
  int var, i, j, a, w, var_offset;
  double alpha;
  struct Level_Set_Data *ls_old;
  alpha = 0.;

  if (d_alpha != NULL) {
    memset(d_alpha->T, 0, sizeof(double) * MDE);
    memset(d_alpha->X, 0, sizeof(double) * DIM * MDE);
    memset(d_alpha->C, 0, sizeof(double) * MAX_CONC * MDE);
    memset(d_alpha->F, 0, sizeof(double) * MDE);
  }

  if (mp->Light_AbsorptionModel == USER) {

    mp->light_absorption = mp->u_light_absorption[0];
    mp->d_light_absorption[TEMPERATURE] = 0;
    for (a = 0; a < DIM; a++) {
      mp->d_light_absorption[MESH_DISPLACEMENT1 + a] = 0;
    }
    for (a = 0; a < pd->Num_Species_Eqn; a++) {
      mp->d_light_absorption[MAX_VARIABLE_TYPES + a] = mp->u_light_absorption[a + 1];
      mp->light_absorption += mp->u_light_absorption[a + 1] * fv->c[a];
    }
    alpha = mp->light_absorption;

    var = TEMPERATURE;
    if (pd->v[pg->imtrx][TEMPERATURE] && d_alpha != NULL) {

      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_alpha->T[j] = mp->d_light_absorption[var] * bf[var]->phi[j];
      }
    }

    if (pd->v[pg->imtrx][MESH_DISPLACEMENT1] && d_alpha != NULL) {
      for (a = 0; a < dim; a++) {
        var = MESH_DISPLACEMENT1 + a;
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_alpha->X[a][j] = mp->d_light_absorption[var] * bf[var]->phi[j];
        }
      }
    }

    if (pd->v[pg->imtrx][MASS_FRACTION] && d_alpha != NULL) {
      for (w = 0; w < pd->Num_Species_Eqn; w++) {
        var = MASS_FRACTION;
        var_offset = MAX_VARIABLE_TYPES + w + 1;
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_alpha->C[w][j] = mp->d_light_absorption[var_offset] * bf[var]->phi[j];
        }
      }
    }
  } else if (mp->Light_AbsorptionModel == CONSTANT) {
    alpha = mp->light_absorption;
  } else if (mp->Light_AbsorptionModel == TABLE) {
    struct Data_Table *table_local;
    table_local = MP_Tables[mp->light_absorption_tableid];
    apply_table_mp(&mp->light_absorption, table_local);
    alpha = mp->light_absorption;

    if (d_alpha != NULL) {
      for (i = 0; i < table_local->columns - 1; i++) {
        var = table_local->t_index[i];
        /* currently only set up to vary w.r.t. temperature */
        switch (var) {
        case TEMPERATURE:
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            d_alpha->T[j] = table_local->slope[i] * bf[var]->phi[j];
          }
          break;
        default:
          GOMA_EH(GOMA_ERROR, "Variable function not yet implemented in material property table");
        }
      }
    }
  } else if (mp->Light_AbsorptionModel == LEVEL_SET) {
    ls_transport_property(mp->u_light_absorption[0], mp->u_light_absorption[1],
                          mp->u_light_absorption[2], &mp->light_absorption,
                          &mp->d_light_absorption[FILL]);

    alpha = mp->light_absorption;

    var = FILL;
    if (pd->v[pg->imtrx][var] && d_alpha != NULL) {

      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_alpha->F[j] = mp->d_light_absorption[var] * bf[var]->phi[j];
      }
    }

  } else {
    GOMA_EH(GOMA_ERROR, "Unrecognized light absorption model");
  }

  if (ls != NULL && mp->Light_AbsorptionModel != LEVEL_SET &&
      mp->Light_AbsorptionModel != LS_QUADRATIC && mp->mp2nd != NULL &&
      mp->mp2nd->LightAbsorptionModel ==
          CONSTANT) /* Only Newtonian constitutive equation allowed for 2nd phase */
  {
    /* kludge for solidification tracking with phase function 0 */
    if (pfd != NULL && pd->e[pg->imtrx][R_EXT_VELOCITY]) {
      ls_old = ls;
      ls = pfd->ls[0];
      alpha = ls_modulate_thermalconductivity(alpha, mp->mp2nd->lightabsorption_phase[0],
                                              ls->Length_Scale,
                                              (double)mp->mp2nd->lightabsorptionmask[0],
                                              (double)mp->mp2nd->lightabsorptionmask[1], d_alpha);
      ls = ls_old;
    }
    alpha = ls_modulate_thermalconductivity(alpha, mp->mp2nd->lightabsorption, ls->Length_Scale,
                                            (double)mp->mp2nd->lightabsorptionmask[0],
                                            (double)mp->mp2nd->lightabsorptionmask[1], d_alpha);
  }

  return (alpha);
}

/*********************************************************************************/
double refractive_index(CONDUCTIVITY_DEPENDENCE_STRUCT *d_n, dbl time)

/**************************************************************************
 *
 * refractive index
 *
 *   Calculate the refractive index and its derivatives wrt to nodal dofs
 *   at the local gauss point
 *
 * Output
 * -----
 *    dbl d_R->T[j]    -> derivative of ref. index wrt the jth
 *                          Temperature unknown in an element
 *    dbl d_R->C[k][j] -> derivative of ref. index wrt the jth
 *                          species unknown of ktype, k, in the element.
 *    dbl d_R->X[a][j] -> derivative of ref. index wrt the jth
 *                          mesh displacement in the ath direction.
 *    dbl d_R->F[j]    -> derivative of ref. index wrt the jth
 *                          FILL unknown in an element
 *
 *  Return
 * --------
 *    Actual value of the refractive index
 ***************************************************************************/

{
  int dim = pd->Num_Dim;
  int var, i, j, a, w, var_offset;
  double n;
  struct Level_Set_Data *ls_old;
  n = 0.;

  if (d_n != NULL) {
    memset(d_n->T, 0, sizeof(double) * MDE);
    memset(d_n->X, 0, sizeof(double) * DIM * MDE);
    memset(d_n->C, 0, sizeof(double) * MAX_CONC * MDE);
    memset(d_n->F, 0, sizeof(double) * MDE);
  }

  if (mp->Refractive_IndexModel == USER) {

    mp->refractive_index = mp->u_refractive_index[0];
    mp->d_refractive_index[TEMPERATURE] = 0;
    for (a = 0; a < DIM; a++) {
      mp->d_refractive_index[MESH_DISPLACEMENT1 + a] = 0;
    }
    for (a = 0; a < pd->Num_Species_Eqn; a++) {
      mp->d_refractive_index[MAX_VARIABLE_TYPES + a] = mp->u_refractive_index[a + 1];
      mp->refractive_index += mp->u_refractive_index[a + 1] * fv->c[a];
    }
    n = mp->refractive_index;

    var = TEMPERATURE;
    if (pd->v[pg->imtrx][TEMPERATURE] && d_n != NULL) {

      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_n->T[j] = mp->d_refractive_index[var] * bf[var]->phi[j];
      }
    }

    if (pd->v[pg->imtrx][MESH_DISPLACEMENT1] && d_n != NULL) {
      for (a = 0; a < dim; a++) {
        var = MESH_DISPLACEMENT1 + a;
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_n->X[a][j] = mp->d_refractive_index[var] * bf[var]->phi[j];
        }
      }
    }

    if (pd->v[pg->imtrx][MASS_FRACTION] && d_n != NULL) {
      for (w = 0; w < pd->Num_Species_Eqn; w++) {
        var = MASS_FRACTION;
        var_offset = MAX_VARIABLE_TYPES + w;
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_n->C[w][j] = mp->d_refractive_index[var_offset] * bf[var]->phi[j];
        }
      }
    }
  } else if (mp->Refractive_IndexModel == CONSTANT) {
    n = mp->refractive_index;
  } else if (mp->Refractive_IndexModel == TABLE) {
    struct Data_Table *table_local;
    table_local = MP_Tables[mp->refractive_index_tableid];
    apply_table_mp(&mp->refractive_index, table_local);
    n = mp->refractive_index;

    if (d_n != NULL) {
      for (i = 0; i < table_local->columns - 1; i++) {
        var = table_local->t_index[i];
        /* currently only set up to vary w.r.t. temperature */
        switch (var) {
        case TEMPERATURE:
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            d_n->T[j] = table_local->slope[i] * bf[var]->phi[j];
          }
          break;
        default:
          GOMA_EH(GOMA_ERROR, "Variable function not yet implemented in material property table");
        }
      }
    }
  } else if (mp->Refractive_IndexModel == LEVEL_SET) {
    ls_transport_property(mp->u_refractive_index[0], mp->u_refractive_index[1],
                          mp->u_refractive_index[2], &mp->refractive_index,
                          &mp->d_refractive_index[FILL]);

    n = mp->refractive_index;

    var = FILL;
    if (pd->v[pg->imtrx][var] && d_n != NULL) {

      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_n->F[j] = mp->d_refractive_index[var] * bf[var]->phi[j];
      }
    }

  } else {
    GOMA_EH(GOMA_ERROR, "Unrecognized refractive index model");
  }

  if (ls != NULL && mp->Refractive_IndexModel != LEVEL_SET &&
      mp->Refractive_IndexModel != LS_QUADRATIC && mp->mp2nd != NULL &&
      mp->mp2nd->RefractiveIndexModel ==
          CONSTANT) /* Only Newtonian constitutive equation allowed for 2nd phase */
  {
    /* kludge for solidification tracking with phase function 0 */
    if (pfd != NULL && pd->e[pg->imtrx][R_EXT_VELOCITY]) {
      ls_old = ls;
      ls = pfd->ls[0];
      n = ls_modulate_thermalconductivity(n, mp->mp2nd->refractiveindex_phase[0], ls->Length_Scale,
                                          (double)mp->mp2nd->refractiveindexmask[0],
                                          (double)mp->mp2nd->refractiveindexmask[1], d_n);
      ls = ls_old;
    }
    n = ls_modulate_thermalconductivity(n, mp->mp2nd->refractiveindex, ls->Length_Scale,
                                        (double)mp->mp2nd->refractiveindexmask[0],
                                        (double)mp->mp2nd->refractiveindexmask[1], d_n);
  }

  return (n);
}
/*********************************************************************************/
/*********************************************************************************/
double extinction_index(CONDUCTIVITY_DEPENDENCE_STRUCT *d_k, dbl time)

/**************************************************************************
 *
 * extinction index
 *
 *   Calculate the extinction index and its derivatives wrt to nodal dofs
 *   at the local gauss point
 *
 * Output
 * -----
 *    dbl d_k->T[j]    -> derivative of ext. index wrt the jth
 *                          Temperature unknown in an element
 *    dbl d_k->C[k][j] -> derivative of ext. index wrt the jth
 *                          species unknown of ktype, k, in the element.
 *    dbl d_k->X[a][j] -> derivative of ext. index wrt the jth
 *                          mesh displacement in the ath direction.
 *    dbl d_k->F[j]    -> derivative of ext. index wrt the jth
 *                          FILL unknown in an element
 *
 *  Return
 * --------
 *    Actual value of the extinction index
 ***************************************************************************/

{
  int dim = pd->Num_Dim;
  int var, i, j, a, w, var_offset;
  double k = 0;
  struct Level_Set_Data *ls_old;

  if (d_k != NULL) {
    memset(d_k->T, 0, sizeof(double) * MDE);
    memset(d_k->X, 0, sizeof(double) * DIM * MDE);
    memset(d_k->C, 0, sizeof(double) * MAX_CONC * MDE);
    memset(d_k->F, 0, sizeof(double) * MDE);
  }

  if (mp->Extinction_IndexModel == USER) {

    mp->extinction_index = mp->u_extinction_index[0];
    mp->d_extinction_index[TEMPERATURE] = 0;
    for (a = 0; a < DIM; a++) {
      mp->d_extinction_index[MESH_DISPLACEMENT1 + a] = 0;
    }
    for (a = 0; a < pd->Num_Species_Eqn; a++) {
      mp->d_extinction_index[MAX_VARIABLE_TYPES + a] = mp->u_extinction_index[a + 1];
      mp->extinction_index += mp->u_extinction_index[a + 1] * fv->c[a];
    }
    k = mp->extinction_index;

    var = TEMPERATURE;
    if (pd->v[pg->imtrx][TEMPERATURE] && d_k != NULL) {

      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_k->T[j] = mp->d_extinction_index[var] * bf[var]->phi[j];
      }
    }

    if (pd->v[pg->imtrx][MESH_DISPLACEMENT1] && d_k != NULL) {
      for (a = 0; a < dim; a++) {
        var = MESH_DISPLACEMENT1 + a;
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_k->X[a][j] = mp->d_extinction_index[var] * bf[var]->phi[j];
        }
      }
    }

    if (pd->v[pg->imtrx][MASS_FRACTION] && d_k != NULL) {
      for (w = 0; w < pd->Num_Species_Eqn; w++) {
        var = MASS_FRACTION;
        var_offset = MAX_VARIABLE_TYPES + w;
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_k->C[w][j] = mp->d_extinction_index[var_offset] * bf[var]->phi[j];
        }
      }
    }
  } else if (mp->Extinction_IndexModel == CONSTANT) {
    k = mp->extinction_index;
  } else if (mp->Extinction_IndexModel == TABLE) {
    struct Data_Table *table_local;
    table_local = MP_Tables[mp->extinction_index_tableid];
    apply_table_mp(&mp->extinction_index, table_local);
    k = mp->extinction_index;

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
  } else if (mp->Extinction_IndexModel == LEVEL_SET) {
    ls_transport_property(mp->u_extinction_index[0], mp->u_extinction_index[1],
                          mp->u_extinction_index[2], &mp->extinction_index,
                          &mp->d_extinction_index[FILL]);

    k = mp->extinction_index;

    var = FILL;
    if (pd->v[pg->imtrx][var] && d_k != NULL) {

      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_k->F[j] = mp->d_extinction_index[var] * bf[var]->phi[j];
      }
    }

  } else {
    GOMA_EH(GOMA_ERROR, "Unrecognized extinction index model");
  }

  if (ls != NULL && mp->Extinction_IndexModel != LEVEL_SET &&
      mp->Extinction_IndexModel != LS_QUADRATIC && mp->mp2nd != NULL &&
      mp->mp2nd->ExtinctionIndexModel ==
          CONSTANT) /* Only Newtonian constitutive equation allowed for 2nd phase */
  {
    /* kludge for solidification tracking with phase function 0 */
    if (pfd != NULL && pd->e[pg->imtrx][R_EXT_VELOCITY]) {
      ls_old = ls;
      ls = pfd->ls[0];
      k = ls_modulate_thermalconductivity(k, mp->mp2nd->extinctionindex_phase[0], ls->Length_Scale,
                                          (double)mp->mp2nd->extinctionindexmask[0],
                                          (double)mp->mp2nd->extinctionindexmask[1], d_k);
      ls = ls_old;
    }
    k = ls_modulate_thermalconductivity(k, mp->mp2nd->extinctionindex, ls->Length_Scale,
                                        (double)mp->mp2nd->extinctionindexmask[0],
                                        (double)mp->mp2nd->extinctionindexmask[1], d_k);
  }

  return (k);
}
/*********************************************************************************/
/*********************************************************************************/

int ls_modulate_momentumsource(double f1[DIM],
                               double f2[DIM],
                               double width,
                               double pm_minus,
                               double pm_plus,
                               MOMENTUM_SOURCE_DEPENDENCE_STRUCT *df)

{
  int i, a, b, var;
  int dim = pd->Num_Dim;
  double factor;

  for (a = 0; a < dim; a++) {

    if (df == NULL) {
      f1[a] = ls_modulate_property(f1[a], f2[a], width, pm_minus, pm_plus, NULL, &factor);

      return (0);
    } else {
      f1[a] = ls_modulate_property(f1[a], f2[a], width, pm_minus, pm_plus, df->F[a], &factor);

      if (pd->v[pg->imtrx][var = TEMPERATURE]) {
        for (i = 0; i < ei[pg->imtrx]->dof[var]; i++) {
          df->T[a][i] *= factor;
        }
      }

      if (pd->v[pg->imtrx][var = MESH_DISPLACEMENT1]) {
        for (b = 0; b < dim; b++) {
          for (i = 0; i < ei[pg->imtrx]->dof[var + b]; i++) {
            df->X[a][b][i] *= factor;
          }
        }
      }

      if (pd->v[pg->imtrx][var = VELOCITY1]) {
        for (b = 0; b < dim; b++) {
          for (i = 0; i < ei[pg->imtrx]->dof[var]; i++) {
            df->v[a][b][i] *= factor;
          }
        }
      }

      if (pd->v[pg->imtrx][var = MASS_FRACTION]) {
        for (b = 0; b < pd->Num_Species; b++) {
          for (i = 0; i < ei[pg->imtrx]->dof[var]; i++) {
            df->C[a][b][i] *= factor;
          }
        }
      }

      if (pd->v[pg->imtrx][var = EFIELD1]) {
        for (b = 0; b < dim; b++) {
          for (i = 0; i < ei[pg->imtrx]->dof[var + b]; i++) {
            df->E[a][b][i] *= factor;
          }
        }
      }
    }
  }
  return (0);
}

void apply_table_mp(double *func, struct Data_Table *table) {
  int i;
  double interp_val, var1[1], slope, temp;

  if (pd->gv[TEMPERATURE]) {
    temp = fv->T;
  } else {
    temp = upd->Process_Temperature;
  }

  for (i = 0; i < table->columns - 1; i++) {
    if (strcmp(table->t_name[i], "TEMPERATURE") == 0) {
      var1[i] = temp;
    } else if (strcmp(table->t_name[i], "MASS_FRACTION") == 0) {
      var1[i] = fv->c[table->species_eq];
    } else if (strcmp(table->t_name[i], "WAVELENGTH") == 0) {
      const double c0 = 1.0 / (sqrt(upd->Free_Space_Permittivity * upd->Free_Space_Permeability));
      dbl freq = upd->EM_Frequency;
      dbl lambda0 = c0 / freq;
      var1[i] = lambda0;
    } else if (strcmp(table->t_name[i], "FAUX_PLASTIC") == 0) {
      GOMA_EH(GOMA_ERROR, "Oops, I shouldn't be using this call for FAUX_PLASTIC.");
    } else {
      GOMA_EH(GOMA_ERROR, "Material Table Model Error-Unknown Function Column ");
    }
  }

  interp_val = interpolate_table(table, var1, &slope, NULL);
  *func = interp_val;
}
/***************************************************************************/
/***************************************************************************/
/***************************************************************************/

double interpolate_table(struct Data_Table *table, double x[], double *sloper, double dfunc_dx[])
/*
 *      A general routine that uses data supplied in a Data_Table
        structure to compute the ordinate of the table that corresponds
        to the abscissa supplied in x.

        Author: Thomas A. Baer, Org 9111
        Date  : July 16, 1998
        Revised: raroach October 12, 1999

    Parameters:
        table = pointer to Data_Table structure
        x     = array of abscissa(s) where ordinate is to be evaluated
        slope = d(ordinate)/dx

    returns:
        value of ordinate at x and the slope there too.
 */
{
  int i, N, iinter, istartx, istarty;
  double func = 0, y1, y2, y3, y4, tt, uu;
  double *t, *t2, *t3, *f;
  double cee, phi[3], xleft, xright;
  int ngrid1, ngrid2, ngrid3;

  N = table->tablelength - 1;
  t = table->t;
  t2 = table->t2;
  t3 = table->t3;
  f = table->f;
  uu = 0.0;

  switch (table->interp_method) {
  case LINEAR: /* This is knucklehead linear interpolation scheme */

    /*  maybe someday you would prefer a more logical fem-based interpolation */

    if (0) {
      for (i = 0; i < N; i++) {
        cee = (x[0] - t[i]) / (t[i + 1] - t[i]);
        if ((cee >= 0.0 && cee <= 1.0) || (cee < 0.0 && i == 0) || (cee > 1.0 && i == N - 1)) {
          phi[0] = -cee + 1.;
          phi[1] = cee;
          table->slope[0] = f[i] * phi[0] + f[i + 1] * phi[1];
          table->slope[1] = f[N + 1 + i] * phi[0] + f[N + 2 + i] * phi[1];
          break;
        }
      }
      table->slope[2] = 0.0;
    } else {
      if (x[0] < t[0]) {
        table->slope[0] = (f[1] - f[0]) / (t[1] - t[0]);
        func = f[0] + (table->slope[0]) * (x[0] - t[0]);
      }

      for (i = 0; x[0] >= t[i] && i < N; i++) {
        if (x[0] >= t[i] && x[0] < t[i + 1]) {
          table->slope[0] = (f[i + 1] - f[i]) / (t[i + 1] - t[i]);
          func = f[i] + (table->slope[0]) * (x[0] - t[i]);
        }
      }
      if (x[0] >= t[N]) {
        table->slope[0] = (f[N] - f[N - 1]) / (t[N] - t[N - 1]);
        func = f[N] + (table->slope[0]) * (x[0] - t[N]);
      }
      table->slope[1] = 0.0;
      *sloper = table->slope[0];
    }
    break;

  case QUADRATIC: /* quadratic lagrangian interpolation scheme */

    if (table->columns == 3) {
      for (i = 0; i < N; i += 2) {
        cee = (x[0] - t[i]) / (t[i + 2] - t[i]);
        if ((cee >= 0.0 && cee <= 1.0) || (cee < 0.0 && i == 0) || (cee > 1.0 && i == N - 2)) {
          phi[0] = 2. * cee * cee - 3. * cee + 1.;
          phi[1] = -4. * cee * cee + 4. * cee;
          phi[2] = 2. * cee * cee - cee;
          table->slope[0] = f[i] * phi[0] + f[i + 1] * phi[1] + f[i + 2] * phi[2];
          table->slope[1] = f[N + 1 + i] * phi[0] + f[N + 2 + i] * phi[1] + f[N + 3 + i] * phi[2];
          break;
        }
      }
      table->slope[2] = 0.0;
    } else {
      for (i = 0; i < N; i += 2) {
        cee = (x[0] - t[i]) / (t[i + 2] - t[i]);
        if ((cee >= 0.0 && cee <= 1.0) || (cee < 0.0 && i == 0) || (cee > 1.0 && i == N - 2)) {
          phi[0] = 2. * cee * cee - 3. * cee + 1.;
          phi[1] = -4. * cee * cee + 4. * cee;
          phi[2] = 2. * cee * cee - cee;
          func = f[i] * phi[0] + f[i + 1] * phi[1] + f[i + 2] * phi[2];
          phi[0] = 4. * cee - 3.;
          phi[1] = -8. * cee + 4.;
          phi[2] = 4. * cee - 1.;
          table->slope[0] =
              (f[i] * phi[0] + f[i + 1] * phi[1] + f[i + 2] * phi[2]) / (t[i + 2] - t[i]);
          break;
        }
      }
      table->slope[1] = 0.0;
      *sloper = table->slope[0];
    }
    break;

  case QUAD_GP: /* quadratic lagrangian interpolation scheme */

    if (table->columns == 3) {
      for (i = 0; i < N; i += 3) {
        xleft =
            (5. + sqrt(15.)) / 6. * t[i] - 2. / 3. * t[i + 1] + (5. - sqrt(15.)) / 6. * t[i + 2];
        xright =
            (5. - sqrt(15.)) / 6. * t[i] - 2. / 3. * t[i + 1] + (5. + sqrt(15.)) / 6. * t[i + 2];
        cee = (x[0] - xleft) / (xright - xleft);
        if ((cee >= 0.0 && cee <= 1.0) || (cee < 0.0 && i == 0) || (cee > 1.0 && i == N - 3)) {
          phi[0] = (20. * cee * cee - 2. * (sqrt(15.) + 10.) * cee + sqrt(15.) + 5.) / 6.;
          phi[1] = (-10. * cee * cee + 10. * cee - 1.) * 2. / 3.;
          phi[2] = (20. * cee * cee + 2. * (sqrt(15.) - 10.) * cee - sqrt(15.) + 5.) / 6.;
          table->slope[0] = f[i] * phi[0] + f[i + 1] * phi[1] + f[i + 2] * phi[2];
          table->slope[1] = f[N + 1 + i] * phi[0] + f[N + 2 + i] * phi[1] + f[N + 3 + i] * phi[2];
          break;
        }
      }
      table->slope[2] = 0.0;
    } else {
      for (i = 0; i < N; i += 3) {
        xleft =
            (5. + sqrt(15.)) / 6. * t[i] - 2. / 3. * t[i + 1] + (5. - sqrt(15.)) / 6. * t[i + 2];
        xright =
            (5. - sqrt(15.)) / 6. * t[i] - 2. / 3. * t[i + 1] + (5. + sqrt(15.)) / 6. * t[i + 2];
        cee = (x[0] - xleft) / (xright - xleft);
        if ((cee >= 0.0 && cee <= 1.0) || (cee < 0.0 && i == 0) || (cee > 1.0 && i == N - 3)) {
          phi[0] = (20. * cee * cee - 2. * (sqrt(15.) + 10.) * cee + sqrt(15.) + 5.) / 6.;
          phi[1] = (-10. * cee * cee + 10. * cee - 1.) * 2. / 3.;
          phi[2] = (20. * cee * cee + 2. * (sqrt(15.) - 10.) * cee - sqrt(15.) + 5.) / 6.;
          func = f[i] * phi[0] + f[i + 1] * phi[1] + f[i + 2] * phi[2];
          phi[0] = (20. * cee - sqrt(15.) + 10.) / 3.;
          phi[1] = (-40. * cee + 20.) / 3.;
          phi[2] = (20. * cee + sqrt(15.) - 10.) / 3.;
          table->slope[0] =
              (f[i] * phi[0] + f[i + 1] * phi[1] + f[i + 2] * phi[2]) / (xright - xleft);
          break;
        }
      }
      table->slope[1] = 0.0;
      *sloper = table->slope[0];
    }
    break;

  case BIQUADRATIC: /* biquadratic lagrangian interpolation scheme */

    ngrid1 = table->tablelength / table->ngrid;
    if (table->columns == 5) {
      table->slope[0] = quad_isomap_invert(x[0], x[1], 0, t, t2, NULL, f, ngrid1, table->ngrid, 1,
                                           2, 2, dfunc_dx);
      table->slope[1] = quad_isomap_invert(x[0], x[1], 0, t, t2, NULL, &f[N + 1], ngrid1,
                                           table->ngrid, 1, 2, 2, dfunc_dx);
      table->slope[2] = quad_isomap_invert(x[0], x[1], 0, t, t2, NULL, &f[2 * N + 2], ngrid1,
                                           table->ngrid, 1, 2, 2, dfunc_dx);
    } else {
      func = quad_isomap_invert(x[0], x[1], 0, t, t2, NULL, f, ngrid1, table->ngrid, 1, 2, 2,
                                dfunc_dx);
    }
    break;

  case BILINEAR: /* BILINEAR Interpolation Scheme */
    /*Find Interval of Different Values of Abscissa #1*/
    iinter = 0;
    for (i = 0; i < N; i++) {
      if (table->t[i] != table->t[i + 1] && iinter == 0) {
        iinter = i + 1;
      }
    }
    if (iinter == 1) {
      fprintf(stderr, " MP Interpolate Error - Need more than 1 point per set");
      GOMA_EH(GOMA_ERROR, "Table interpolation not implemented");
    }

    istartx = iinter;

    for (i = iinter; t[i] <= x[0] && i < N - iinter + 1; i = i + iinter) {
      istartx = i + iinter;
    }

    istarty = istartx;

    for (i = istartx + 1; t2[i] <= x[1] && i < istartx + iinter - 1; i++) {
      istarty = i;
    }

    y1 = f[istarty];
    y2 = f[istarty + 1];
    y3 = f[istarty + 1 - iinter];
    y4 = f[istarty - iinter];

    tt = (x[1] - t2[istarty]) / (t2[istarty + 1] - t2[istarty]);
    uu = (x[0] - t[istarty]) / (t[istarty + 1 - iinter] - t[istarty]);

    func = (1. - tt) * (1. - uu) * y1 + tt * (1. - uu) * y2 + tt * uu * y3 + (1. - tt) * uu * y4;

    table->slope[1] = (f[istarty + 1] - f[istarty]) / (t2[istarty + 1] - t2[istarty]);
    table->slope[0] =
        (f[istarty + 1] - f[istarty + 1 - iinter]) / (t[istarty + 1] - t[istarty + 1 - iinter]);
    *sloper = table->slope[0];
    if (dfunc_dx != NULL) {
      dfunc_dx[0] = table->slope[0];
      dfunc_dx[1] = table->slope[1];
    }
    break;

  case TRILINEAR: /* trilinear lagrangian interpolation scheme */

    ngrid1 = table->ngrid;
    ngrid2 = table->ngrid2 / table->ngrid;
    ngrid3 = table->tablelength / table->ngrid2;
    func =
        quad_isomap_invert(x[0], x[1], x[2], t, t2, t3, f, ngrid1, ngrid2, ngrid3, 1, 3, dfunc_dx);
    break;

  case TRIQUADRATIC: /* triquadratic lagrangian interpolation scheme */

    ngrid1 = table->ngrid;
    ngrid2 = table->ngrid2 / table->ngrid;
    ngrid3 = table->tablelength / table->ngrid2;
    func =
        quad_isomap_invert(x[0], x[1], x[2], t, t2, t3, f, ngrid1, ngrid2, ngrid3, 2, 3, dfunc_dx);
    break;

  default:
    GOMA_EH(GOMA_ERROR, "Table interpolation order not implemented");
  }

  return (func);
}

/***************************************************************************/

double
table_distance_search(struct Data_Table *table, double x[], double *sloper, double dfunc_dx[])
/*
 *      A routine to find the shortest distance from a point
 *      to a table represented curve

        Author: Robert B. Secor
        Date  : December 5, 2014
        Revised:

    Parameters:
        table = pointer to Data_Table structure
        x     = array of point coordinates
        slope = d(ordinate)/dx

    returns:
        value of distance
 */
{
  int i, N, iinter, istartx, istarty, basis, ordinate, i_min = 0, elem, iter;
  double func = 0, y1, y2, y3, y4, tt, uu;
  double *t, *t2, *t3, *f;
  double cee, phi[3], xleft, xright, delta, epstol = 0.0001, update, point;
  double phic[3] = {0, 0, 0};
  int ngrid1, ngrid2, ngrid3;
  double dist, dist_min = BIG_PENALTY;

  N = table->tablelength - 1;
  t = table->t;
  t2 = table->t2;
  t3 = table->t3;
  f = table->f;
  if (table->t_index[0] == MESH_POSITION1) {
    basis = 0;
    ordinate = 1;
  } else if (table->t_index[0] == MESH_POSITION2) {
    basis = 1;
    ordinate = 0;
  } else {
    basis = 0;
    ordinate = 1;
  }
  uu = 0.0;

  switch (table->interp_method) {
  case LINEAR: /* This is knucklehead linear interpolation scheme */

    /*  maybe someday you would prefer a more logical fem-based interpolation */

    if (0) {
      for (i = 0; i < N; i++) {
        cee = (x[0] - t[i]) / (t[i + 1] - t[i]);
        if ((cee >= 0.0 && cee <= 1.0) || (cee < 0.0 && i == 0) || (cee > 1.0 && i == N - 1)) {
          phi[0] = -cee + 1.;
          phi[1] = cee;
          table->slope[0] = f[i] * phi[0] + f[i + 1] * phi[1];
          table->slope[1] = f[N + 1 + i] * phi[0] + f[N + 2 + i] * phi[1];
          break;
        }
      }
      table->slope[2] = 0.0;
    } else {
      if (x[0] < t[0]) {
        table->slope[0] = (f[1] - f[0]) / (t[1] - t[0]);
        func = f[0] + (table->slope[0]) * (x[0] - t[0]);
      }

      for (i = 0; x[0] >= t[i] && i < N; i++) {
        if (x[0] >= t[i] && x[0] < t[i + 1]) {
          table->slope[0] = (f[i + 1] - f[i]) / (t[i + 1] - t[i]);
          func = f[i] + (table->slope[0]) * (x[0] - t[i]);
        }
      }
      if (x[0] >= t[N]) {
        table->slope[0] = (f[N] - f[N - 1]) / (t[N] - t[N - 1]);
        func = f[N] + (table->slope[0]) * (x[0] - t[N]);
      }
      table->slope[1] = 0.0;
      *sloper = table->slope[0];
    }
    break;

  case QUADRATIC: /* quadratic lagrangian interpolation scheme */

    if (table->columns == 3) {
      for (i = 0; i < N; i += 2) {
        cee = (x[0] - t[i]) / (t[i + 2] - t[i]);
        if ((cee >= 0.0 && cee <= 1.0) || (cee < 0.0 && i == 0) || (cee > 1.0 && i == N - 2)) {
          phi[0] = 2. * cee * cee - 3. * cee + 1.;
          phi[1] = -4. * cee * cee + 4. * cee;
          phi[2] = 2. * cee * cee - cee;
          table->slope[0] = f[i] * phi[0] + f[i + 1] * phi[1] + f[i + 2] * phi[2];
          table->slope[1] = f[N + 1 + i] * phi[0] + f[N + 2 + i] * phi[1] + f[N + 3 + i] * phi[2];
          break;
        }
      }
      table->slope[2] = 0.0;
    } else {
      /* search through points for minimum distance*/
      dist_min = dist = sqrt(SQUARE(x[basis] - t[0]) + SQUARE(x[ordinate] - f[0]));
      for (i = 1; i < N; i++) {
        dist = sqrt(SQUARE(x[basis] - t[i]) + SQUARE(x[ordinate] - f[i]));
        if (dist < dist_min) {
          i_min = i;
          dist_min = dist;
        }
      }
      elem = i_min / 2 - 1;
      cee = 0.;
      for (iter = 0; iter < 10; iter++) {
        phi[0] = 2. * cee * cee - 3. * cee + 1.;
        phi[1] = -4. * cee * cee + 4. * cee;
        phi[2] = 2. * cee * cee - cee;
        func = f[2 * elem] * phi[0] + f[2 * elem + 1] * phi[1] + f[2 * elem + 2] * phi[2];
        phic[0] = 4. * cee - 3.;
        phic[1] = -8. * cee + 4.;
        phic[2] = 4. * cee - 1.;
        delta = t[2 * elem + 2] - t[2 * elem];
        *sloper = (f[2 * elem] * phic[0] + f[2 * elem + 1] * phic[1] + f[2 * elem + 2] * phic[2]);
        point = t[2 * elem] + cee * delta;
        update = ((point - x[basis]) * delta + (func - x[ordinate]) * (*sloper)) /
                 (SQUARE(delta) + SQUARE(*sloper));
        cee -= update;
        if (fabs(update) < epstol)
          break;
        if (cee < 0.0 && elem > 0) {
          cee += 1.0;
          elem -= 1;
        }
        if (cee > 1.0 && elem < ((N - 1) / 2 - 1)) {
          cee -= 1.0;
          elem += 1;
        }
      }
      dist_min = sqrt(SQUARE(x[basis] - point) + SQUARE(x[ordinate] - func));
      table->slope[0] = dfunc_dx[basis] = (point - x[basis]) / dist_min;
      dfunc_dx[ordinate] = (func - x[ordinate]) / dist_min;
      table->slope[1] = point;
      table->slope[2] = func;
      *sloper = table->slope[0];
    }
    break;

  case QUAD_GP: /* quadratic lagrangian interpolation scheme */

    if (table->columns == 3) {
      for (i = 0; i < N; i += 3) {
        xleft =
            (5. + sqrt(15.)) / 6. * t[i] - 2. / 3. * t[i + 1] + (5. - sqrt(15.)) / 6. * t[i + 2];
        xright =
            (5. - sqrt(15.)) / 6. * t[i] - 2. / 3. * t[i + 1] + (5. + sqrt(15.)) / 6. * t[i + 2];
        cee = (x[0] - xleft) / (xright - xleft);
        if ((cee >= 0.0 && cee <= 1.0) || (cee < 0.0 && i == 0) || (cee > 1.0 && i == N - 3)) {
          phi[0] = (20. * cee * cee - 2. * (sqrt(15.) + 10.) * cee + sqrt(15.) + 5.) / 6.;
          phi[1] = (-10. * cee * cee + 10. * cee - 1.) * 2. / 3.;
          phi[2] = (20. * cee * cee + 2. * (sqrt(15.) - 10.) * cee - sqrt(15.) + 5.) / 6.;
          table->slope[0] = f[i] * phi[0] + f[i + 1] * phi[1] + f[i + 2] * phi[2];
          table->slope[1] = f[N + 1 + i] * phi[0] + f[N + 2 + i] * phi[1] + f[N + 3 + i] * phi[2];
          break;
        }
      }
      table->slope[2] = 0.0;
    } else {
      for (i = 0; i < N; i += 3) {
        xleft =
            (5. + sqrt(15.)) / 6. * t[i] - 2. / 3. * t[i + 1] + (5. - sqrt(15.)) / 6. * t[i + 2];
        xright =
            (5. - sqrt(15.)) / 6. * t[i] - 2. / 3. * t[i + 1] + (5. + sqrt(15.)) / 6. * t[i + 2];
        cee = (x[0] - xleft) / (xright - xleft);
        if ((cee >= 0.0 && cee <= 1.0) || (cee < 0.0 && i == 0) || (cee > 1.0 && i == N - 3)) {
          phi[0] = (20. * cee * cee - 2. * (sqrt(15.) + 10.) * cee + sqrt(15.) + 5.) / 6.;
          phi[1] = (-10. * cee * cee + 10. * cee - 1.) * 2. / 3.;
          phi[2] = (20. * cee * cee + 2. * (sqrt(15.) - 10.) * cee - sqrt(15.) + 5.) / 6.;
          func = f[i] * phi[0] + f[i + 1] * phi[1] + f[i + 2] * phi[2];
          phi[0] = (20. * cee - sqrt(15.) + 10.) / 3.;
          phi[1] = (-40. * cee + 20.) / 3.;
          phi[2] = (20. * cee + sqrt(15.) - 10.) / 3.;
          table->slope[0] =
              (f[i] * phi[0] + f[i + 1] * phi[1] + f[i + 2] * phi[2]) / (xright - xleft);
          break;
        }
      }
      table->slope[1] = 0.0;
      *sloper = table->slope[0];
    }
    break;

  case BIQUADRATIC: /* biquadratic lagrangian interpolation scheme */

    ngrid1 = table->tablelength / table->ngrid;
    if (table->columns == 5) {
      table->slope[0] = quad_isomap_invert(x[0], x[1], 0, t, t2, NULL, f, ngrid1, table->ngrid, 1,
                                           2, 2, dfunc_dx);
      table->slope[1] = quad_isomap_invert(x[0], x[1], 0, t, t2, NULL, &f[N + 1], ngrid1,
                                           table->ngrid, 1, 2, 2, dfunc_dx);
      table->slope[2] = quad_isomap_invert(x[0], x[1], 0, t, t2, NULL, &f[2 * N + 2], ngrid1,
                                           table->ngrid, 1, 2, 2, dfunc_dx);
    } else {
      func = quad_isomap_invert(x[0], x[1], 0, t, t2, NULL, f, ngrid1, table->ngrid, 1, 2, 2,
                                dfunc_dx);
    }
    break;

  case BILINEAR: /* BILINEAR Interpolation Scheme */
    /*Find Interval of Different Values of Abscissa #1*/
    iinter = 0;
    for (i = 0; i < N; i++) {
      if (table->t[i] != table->t[i + 1] && iinter == 0) {
        iinter = i + 1;
      }
    }
    if (iinter == 1) {
      fprintf(stderr, " MP Interpolate Error - Need more than 1 point per set");
      GOMA_EH(GOMA_ERROR, "Table interpolation not implemented");
    }

    istartx = iinter;

    for (i = iinter; t[i] <= x[0] && i < N - iinter + 1; i = i + iinter) {
      istartx = i + iinter;
    }

    istarty = istartx;

    for (i = istartx + 1; t2[i] <= x[1] && i < istartx + iinter - 1; i++) {
      istarty = i;
    }

    y1 = f[istarty];
    y2 = f[istarty + 1];
    y3 = f[istarty + 1 - iinter];
    y4 = f[istarty - iinter];

    tt = (x[1] - t2[istarty]) / (t2[istarty + 1] - t2[istarty]);
    uu = (x[0] - t[istarty]) / (t[istarty + 1 - iinter] - t[istarty]);

    func = (1. - tt) * (1. - uu) * y1 + tt * (1. - uu) * y2 + tt * uu * y3 + (1. - tt) * uu * y4;

    table->slope[1] = (f[istarty + 1] - f[istarty]) / (t2[istarty + 1] - t2[istarty]);
    table->slope[0] =
        (f[istarty + 1] - f[istarty + 1 - iinter]) / (t[istarty + 1] - t[istarty + 1 - iinter]);
    *sloper = table->slope[0];
    if (dfunc_dx != NULL) {
      dfunc_dx[0] = table->slope[0];
      dfunc_dx[1] = table->slope[1];
    }
    break;

  case TRILINEAR: /* trilinear lagrangian interpolation scheme */

    ngrid1 = table->ngrid;
    ngrid2 = table->ngrid2 / table->ngrid;
    ngrid3 = table->tablelength / table->ngrid2;
    func =
        quad_isomap_invert(x[0], x[1], x[2], t, t2, t3, f, ngrid1, ngrid2, ngrid3, 1, 3, dfunc_dx);
    break;

  case TRIQUADRATIC: /* triquadratic lagrangian interpolation scheme */

    ngrid1 = table->ngrid;
    ngrid2 = table->ngrid2 / table->ngrid;
    ngrid3 = table->tablelength / table->ngrid2;
    func =
        quad_isomap_invert(x[0], x[1], x[2], t, t2, t3, f, ngrid1, ngrid2, ngrid3, 2, 3, dfunc_dx);
    break;

  default:
    GOMA_EH(GOMA_ERROR, "Table interpolation order not implemented");
  }

  return (dist_min);
}
/*******************************************************************************
 *
 * continuous_surface_tension(): This function computes a tensor (csf) which is
 *                               a continuous, volumetric representation of the
 *           capillary stress that arises at a fluid/fluid interface (viz., the
 *           traction condition).  This is only used in the level set
 *           formulation.
 *
 * Input
 * -----
 *   Nothing.
 *
 * Ouput
 * -----
 *   csf	  = The "Continuous Surface Force" contribution to the fluid
 *                  stress tensor.
 *   d_csf_dF	  = Jacobian derivative of csf w.r.t. FILL.  If this is NULL,
 *                  then no Jacobian terms are computed.
 *
 *   NB: It is assumed that csf and d_csf_dF are zeroed.
 *
 * Returns
 * -------
 *   0 = Success.
 *  -1 = Failure.
 *
 ******************************************************************************/
int continuous_surface_tension(double st,
                               double csf[DIM][DIM],
                               double d_csf_dF[DIM][DIM][MDE],
                               double d_csf_dX[DIM][DIM][DIM][MDE]) {
  int a, b, c, j;
  int status = 0;
  int do_deriv;

  int var;

  var = ls->var;

  /*  double st = mp->surface_tension; */

  /* See if we're supposed to compute derivatives. */
  do_deriv = d_csf_dF != NULL;

  /* Fetch the level set interface functions. */
  load_lsi(ls->Length_Scale);

  /* If we're not near the zero level set, then csf is a zero tensor. */
  if (!lsi->near)
    return (status);

  /* Calculate the CSF tensor. */
  for (a = 0; a < VIM; a++) {
    for (b = 0; b < VIM; b++) {
      csf[a][b] = st * lsi->delta * ((double)delta(a, b) - lsi->normal[a] * lsi->normal[b]);
    }
  }

  /* Bail out if we're done. */
  if (!do_deriv)
    return (status);

  load_lsi_derivs();

  /* Calculate the derivatives. */
  for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
#ifdef DO_NO_UNROLL
    for (a = 0; a < VIM; a++) {
      for (b = 0; b < VIM; b++) {
        d_csf_dF[a][b][j] =
            st * lsi->d_delta_dF[j] * ((double)delta(a, b) - lsi->normal[a] * lsi->normal[b]) -
            st * lsi->delta * lsi->d_normal_dF[a][j] * lsi->normal[b] -
            st * lsi->delta * lsi->normal[a] * lsi->d_normal_dF[b][j];
      }
    }
#else
    /*(0,0) */
    d_csf_dF[0][0][j] =
        st * lsi->d_delta_dF[j] * ((double)delta(0, 0) - lsi->normal[0] * lsi->normal[0]) -
        st * lsi->delta * lsi->d_normal_dF[0][j] * lsi->normal[0] -
        st * lsi->delta * lsi->normal[0] * lsi->d_normal_dF[0][j];

    /* (1,1) */
    d_csf_dF[1][1][j] =
        st * lsi->d_delta_dF[j] * ((double)delta(1, 1) - lsi->normal[1] * lsi->normal[1]) -
        st * lsi->delta * lsi->d_normal_dF[1][j] * lsi->normal[1] -
        st * lsi->delta * lsi->normal[1] * lsi->d_normal_dF[1][j];
    /* (0,1) */
    d_csf_dF[0][1][j] =
        st * lsi->d_delta_dF[j] * ((double)delta(0, 1) - lsi->normal[0] * lsi->normal[1]) -
        st * lsi->delta * lsi->d_normal_dF[0][j] * lsi->normal[1] -
        st * lsi->delta * lsi->normal[0] * lsi->d_normal_dF[1][j];
    /* (1,0) */
    d_csf_dF[1][0][j] =
        st * lsi->d_delta_dF[j] * ((double)delta(1, 0) - lsi->normal[1] * lsi->normal[0]) -
        st * lsi->delta * lsi->d_normal_dF[1][j] * lsi->normal[0] -
        st * lsi->delta * lsi->normal[1] * lsi->d_normal_dF[0][j];
    if (VIM == 3) {
      /* (2,2) */
      d_csf_dF[2][2][j] =
          st * lsi->d_delta_dF[j] * ((double)delta(2, 2) - lsi->normal[2] * lsi->normal[2]) -
          st * lsi->delta * lsi->d_normal_dF[2][j] * lsi->normal[2] -
          st * lsi->delta * lsi->normal[2] * lsi->d_normal_dF[2][j];
      /* (2,1) */
      d_csf_dF[2][1][j] =
          st * lsi->d_delta_dF[j] * ((double)delta(2, 1) - lsi->normal[2] * lsi->normal[1]) -
          st * lsi->delta * lsi->d_normal_dF[2][j] * lsi->normal[1] -
          st * lsi->delta * lsi->normal[2] * lsi->d_normal_dF[1][j];
      /* (2,0) */
      d_csf_dF[2][0][j] =
          st * lsi->d_delta_dF[j] * ((double)delta(2, 0) - lsi->normal[2] * lsi->normal[0]) -
          st * lsi->delta * lsi->d_normal_dF[2][j] * lsi->normal[0] -
          st * lsi->delta * lsi->normal[2] * lsi->d_normal_dF[0][j];
      /* (1,2) */
      d_csf_dF[1][2][j] =
          st * lsi->d_delta_dF[j] * ((double)delta(1, 2) - lsi->normal[1] * lsi->normal[2]) -
          st * lsi->delta * lsi->d_normal_dF[1][j] * lsi->normal[2] -
          st * lsi->delta * lsi->normal[1] * lsi->d_normal_dF[2][j];
      /* (0,2)  */
      d_csf_dF[0][2][j] =
          st * lsi->d_delta_dF[j] * ((double)delta(0, 2) - lsi->normal[0] * lsi->normal[2]) -
          st * lsi->delta * lsi->d_normal_dF[0][j] * lsi->normal[2] -
          st * lsi->delta * lsi->normal[0] * lsi->d_normal_dF[2][j];
    }

#endif
  }

  if (pd->v[pg->imtrx][MESH_DISPLACEMENT1]) {
    for (a = 0; a < VIM; a++) {
      for (b = 0; b < VIM; b++) {
        for (c = 0; c < VIM; c++) {
          var = MESH_DISPLACEMENT1 + c;

          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            d_csf_dX[a][b][c][j] = st * lsi->delta *
                                   (-lsi->d_normal_dmesh[a][c][j] * lsi->normal[b] -
                                    lsi->normal[a] * lsi->d_normal_dmesh[b][c][j]);
          }
        }
      }
    }
  }
  return (status);
}

/*      routine for inverting isoparametric map
        finds cee,eta coordinates given the x,y coords
        and interpolates function

        Points on a rectangular mapped grid are assumed.
        It should work for either linear or quadratic interpolation

        Try to include 3D bricks also - 9/15/2004
*/
double quad_isomap_invert(const double x1,
                          const double y1,
                          const double z1,
                          const double x[],
                          const double y[],
                          const double z[],
                          const double p[],
                          const int ngrid1,
                          const int ngrid2,
                          const int ngrid3,
                          const int elem_order,
                          const int dim,
                          double dfunc_dx[]) {
  int i, j, k, l, iter;
  double pt[5] = {0.5, 1.0, 0.0, 0.5, 1.0};
  double phi[27], phic[27], phie[27], phig[27];
  double coord[3][27], xpt[3], pp[27];
  double jac[3][3] = {{0.0}};
  double detjt, detjti = 0.0;
  static double xi[3] = {0.5, 0.5, 0.5};
  static int nell = 0;
  int itp[27];
  int nell_xi[3], ne_xi[3];
  double dxi[3], eps, pvalue, pc, pe, pg;
  double epstol = 1.0e-4;
  double loc_tol = 0.0001;
  int max_its = 20;
  int elem_nodes;
  double bfi, bfj, dfi, dfj, bfl, dfl;

  /* Initialize dxi */
  dxi[0] = 0;
  dxi[1] = 0;
  dxi[2] = 0;

  elem_nodes = pow(elem_order + 1, dim);
  switch (dim) {
  case 2:
    ne_xi[0] = (ngrid1 - 1) / elem_order;
    ne_xi[1] = (ngrid2 - 1) / elem_order;
    nell_xi[0] = nell / ne_xi[1];
    nell_xi[1] = nell % ne_xi[1];
    break;
  case 3:
    ne_xi[0] = (ngrid1 - 1) / elem_order;
    ne_xi[1] = (ngrid2 - 1) / elem_order;
    ne_xi[2] = (ngrid3 - 1) / elem_order;
    nell_xi[2] = nell / (ne_xi[0] * ne_xi[1]);
    nell_xi[0] = nell % ne_xi[0];
    nell_xi[1] = (nell - nell_xi[2] * ne_xi[0] * ne_xi[1]) / ne_xi[0];
    break;
  }

  /*      iterate to find isoparametric coordinates  */

  iter = 0;
  eps = 10. * epstol;
  while (iter < max_its && eps > epstol) {
    iter++;

    /**     evaluate basis functions  **/

    switch (elem_order) {
    case 2:
      switch (dim) {
      case 2:
        for (i = 0; i < 3; i++) {
          k = 3 * i;
          bfi = (xi[0] - pt[i]) * (xi[0] - pt[i + 1]) /
                ((pt[i + 2] - pt[i]) * (pt[i + 2] - pt[i + 1]));
          dfi = (2. * xi[0] - pt[i] - pt[i + 1]) / ((pt[i + 2] - pt[i]) * (pt[i + 2] - pt[i + 1]));
          for (j = 0; j < 3; j++) {
            bfj = (xi[1] - pt[j]) * (xi[1] - pt[j + 1]) /
                  ((pt[j + 2] - pt[j]) * (pt[j + 2] - pt[j + 1]));
            dfj =
                (2. * xi[1] - pt[j] - pt[j + 1]) / ((pt[j + 2] - pt[j]) * (pt[j + 2] - pt[j + 1]));
            phi[k + j] = bfi * bfj;
            phic[k + j] = dfi * bfj;
            phie[k + j] = bfi * dfj;
          }
        }
        break;
      case 3:
        for (i = 0; i < 3; i++) {
          k = 3 * i;
          bfi = (xi[0] - pt[i]) * (xi[0] - pt[i + 1]) /
                ((pt[i + 2] - pt[i]) * (pt[i + 2] - pt[i + 1]));
          dfi = (2. * xi[0] - pt[i] - pt[i + 1]) / ((pt[i + 2] - pt[i]) * (pt[i + 2] - pt[i + 1]));
          for (j = 0; j < 3; j++) {
            bfj = (xi[1] - pt[j]) * (xi[1] - pt[j + 1]) /
                  ((pt[j + 2] - pt[j]) * (pt[j + 2] - pt[j + 1]));
            dfj =
                (2. * xi[1] - pt[j] - pt[j + 1]) / ((pt[j + 2] - pt[j]) * (pt[j + 2] - pt[j + 1]));
            for (l = 0; l < 3; l++) {
              bfl = (xi[1] - pt[l]) * (xi[1] - pt[l + 1]) /
                    ((pt[l + 2] - pt[l]) * (pt[l + 2] - pt[l + 1]));
              dfl = (2. * xi[1] - pt[l] - pt[l + 1]) /
                    ((pt[l + 2] - pt[l]) * (pt[l + 2] - pt[l + 1]));
              phi[k + j + 9 * l] = bfi * bfj * bfl;
              phic[k + j + 9 * l] = dfi * bfj * bfl;
              phie[k + j + 9 * l] = bfi * dfj * bfl;
              phig[k + j + 9 * l] = bfi * bfj * dfl;
            }
          }
        }
        break;
      }
      break;

    case 1:
      switch (dim) {
      case 2:
        phi[0] = (1. - xi[0]) * (1. - xi[1]);
        phi[1] = (1. - xi[0]) * xi[1];
        phi[2] = xi[0] * (1. - xi[1]);
        phi[3] = xi[0] * xi[1];
        phic[0] = xi[1] - 1.;
        phic[1] = -xi[1];
        phic[2] = 1. - xi[1];
        phic[3] = xi[1];
        phie[0] = xi[0] - 1.;
        phie[1] = 1. - xi[0];
        phie[2] = -xi[0];
        phie[3] = xi[0];
        break;
      case 3:
        phi[0] = (1. - xi[0]) * (1. - xi[1]) * (1. - xi[2]);
        phi[1] = (1. - xi[0]) * xi[1] * (1. - xi[2]);
        phi[2] = xi[0] * (1. - xi[1]) * (1. - xi[2]);
        phi[3] = xi[0] * xi[1] * (1. - xi[2]);
        phi[4] = (1. - xi[0]) * (1. - xi[1]) * xi[2];
        phi[5] = (1. - xi[0]) * xi[1] * xi[2];
        phi[6] = xi[0] * (1. - xi[1]) * xi[2];
        phi[7] = xi[0] * xi[1] * xi[2];
        phic[0] = -(1. - xi[1]) * (1. - xi[2]);
        phic[1] = -xi[1] * (1. - xi[2]);
        phic[2] = (1. - xi[1]) * (1. - xi[2]);
        phic[3] = xi[1] * (1. - xi[2]);
        phic[4] = (xi[1] - 1.) * xi[2];
        phic[5] = -xi[1] * xi[2];
        phic[6] = (1. - xi[1]) * xi[2];
        phic[7] = xi[1] * xi[2];
        phie[0] = -(1. - xi[0]) * (1. - xi[2]);
        phie[1] = (1. - xi[0]) * (1. - xi[2]);
        phie[2] = -xi[0] * (1. - xi[2]);
        phie[3] = xi[0] * (1. - xi[2]);
        phie[4] = -(1. - xi[0]) * xi[2];
        phie[5] = (1. - xi[0]) * xi[2];
        phie[6] = -xi[0] * xi[2];
        phie[7] = xi[0] * xi[2];
        phig[0] = -(1. - xi[0]) * (1. - xi[1]);
        phig[1] = -(1. - xi[0]) * xi[1];
        phig[2] = -xi[0] * (1. - xi[1]);
        phig[3] = -xi[0] * xi[1];
        phig[4] = (1. - xi[0]) * (1. - xi[1]);
        phig[5] = (1. - xi[0]) * xi[1];
        phig[6] = xi[0] * (1. - xi[1]);
        phig[7] = xi[0] * xi[1];
        break;
      }
      break;

    default:
      fprintf(stderr, "\n Unknown element order - quad_isomap\n");
      GOMA_EH(GOMA_ERROR, "Fatal Error");
    }

    /**  gather local node numbers  **/

    switch (dim) {
    case 2:
      itp[0] = elem_order * (elem_order * ne_xi[1] + 1) * nell_xi[0] + elem_order * nell_xi[1];
      itp[1] = itp[0] + 1;
      if (elem_order == 2)
        itp[2] = itp[1] + 1;
      for (i = 1; i <= elem_order; i++) {
        for (j = 0; j <= elem_order; j++) {
          itp[(elem_order + 1) * i + j] = itp[j] + (elem_order * ne_xi[1] + 1) * i;
        }
      }
      break;
    case 3:
      itp[0] = nell_xi[2] * elem_order * (elem_order * ne_xi[0] + 1) * (elem_order * ne_xi[1] + 1) +
               elem_order * (elem_order * ne_xi[0] + 1) * nell_xi[1] + elem_order * nell_xi[0];
      itp[1] = itp[0] + elem_order * (elem_order * ne_xi[0] + 1);
      if (elem_order == 2)
        itp[2] = itp[1] + elem_order * (elem_order * ne_xi[0] + 1);
      for (i = 0; i <= elem_order; i++) {
        for (j = 0; j <= elem_order; j++) {
          for (k = 0; k <= elem_order; k++) {
            itp[(elem_order + 1) * i + j + SQUARE(elem_order + 1) * k] =
                itp[j] + i + k * (elem_order * ne_xi[0] + 1) * (elem_order * ne_xi[1] + 1);
          }
        }
      }
      break;
    }

    /**   gather global variables into local arrays   */

    switch (dim) {
    case 3:
      for (i = 0; i < elem_nodes; i++) {
        coord[2][i] = z[itp[i]];
      }
      /* fall through */
    case 2:
      for (i = 0; i < elem_nodes; i++) {
        coord[0][i] = x[itp[i]];
        coord[1][i] = y[itp[i]];
      }
      break;
    }

    memset(xpt, 0, sizeof(dbl) * 3);
    memset(jac, 0, sizeof(dbl) * 9);

    switch (dim) {
    case 3:
      for (i = 0; i < elem_nodes; i++) {
        jac[0][2] += coord[0][i] * phig[i];
        jac[1][2] += coord[1][i] * phig[i];
        jac[2][0] += coord[2][i] * phic[i];
        jac[2][1] += coord[2][i] * phie[i];
        jac[2][2] += coord[2][i] * phig[i];
        xpt[2] += coord[2][i] * phi[i];
      }
      /* fall through */
    case 2:
      for (i = 0; i < elem_nodes; i++) {
        jac[0][0] += coord[0][i] * phic[i];
        jac[0][1] += coord[0][i] * phie[i];
        jac[1][0] += coord[1][i] * phic[i];
        jac[1][1] += coord[1][i] * phie[i];
        xpt[0] += coord[0][i] * phi[i];
        xpt[1] += coord[1][i] * phi[i];
      }
      break;
    }

    switch (dim) {
    case 2:
      detjt = jac[0][0] * jac[1][1] - jac[0][1] * jac[1][0];
      detjti = 1. / detjt;
      dxi[0] = (-jac[1][1] * (xpt[0] - x1) + jac[0][1] * (xpt[1] - y1)) * detjti;
      dxi[1] = (jac[1][0] * (xpt[0] - x1) - jac[0][0] * (xpt[1] - y1)) * detjti;
      break;
    case 3:
      detjt = jac[0][0] * (jac[1][1] * jac[2][2] - jac[2][1] * jac[1][2]) -
              jac[1][0] * (jac[0][1] * jac[2][2] - jac[2][1] * jac[0][2]) +
              jac[2][0] * (jac[1][2] * jac[0][1] - jac[1][1] * jac[0][2]);
      detjti = 1. / detjt;
      dxi[0] = detjti * (-(xpt[0] - x1) * (jac[2][2] * jac[1][1] - jac[2][1] * jac[1][2]) -
                         (xpt[1] - y1) * (jac[2][1] * jac[0][2] - jac[2][2] * jac[0][1]) -
                         (xpt[2] - z1) * (jac[1][2] * jac[0][1] - jac[1][1] * jac[0][2]));
      dxi[1] = detjti * (-(xpt[0] - x1) * (jac[2][0] * jac[1][2] - jac[2][2] * jac[1][0]) -
                         (xpt[1] - y1) * (jac[2][2] * jac[0][0] - jac[2][0] * jac[0][2]) -
                         (xpt[2] - z1) * (jac[1][0] * jac[0][2] - jac[1][2] * jac[0][0]));
      dxi[2] = detjti * (-(xpt[0] - x1) * (jac[2][1] * jac[1][0] - jac[2][0] * jac[1][1]) -
                         (xpt[1] - y1) * (jac[2][0] * jac[0][1] - jac[2][1] * jac[0][0]) -
                         (xpt[2] - z1) * (jac[1][1] * jac[0][0] - jac[1][0] * jac[0][1]));
      break;
    }

    /*newton updates  */

    eps = 0.0;
    for (i = 0; i < dim; i++) {
      xi[i] += dxi[i];
      eps += dxi[i] * dxi[i];
    }
    eps = sqrt(eps) / ((double)dim);

    /*check to see if crossed element boundaries   */

    for (i = 0; i < dim; i++) {
      double fnum;
      int dell;
      if (xi[i] < -loc_tol && nell_xi[i] > 0) {
        fnum = floor(xi[i]);
        dell = MAX(0, nell_xi[i] + fnum - 1) - nell_xi[i];
        if (iter > max_its / 2)
          dell = -1;
        nell_xi[i] += dell;
        xi[i] -= ((double)dell);
        if (dell != 0)
          eps = 10.0 * epstol;
      } else if (xi[i] > 1.0 + loc_tol && nell_xi[i] < ne_xi[i] - 1) {
        fnum = floor(xi[i]);
        dell = MIN(ne_xi[i] - 1, nell_xi[i] + fnum) - nell_xi[i];
        if (iter > max_its / 2)
          dell = 1;
        nell_xi[i] += dell;
        xi[i] -= ((double)dell);
        if (dell != 0)
          eps = 10.0 * epstol;
      }
    }
  }

  /*check to see if iteration was successful      */

  if (eps > epstol) {
    fprintf(stderr, " pquad iteration failed %d %g %g %g\n", iter, xi[0], xi[1], xi[2]);
    fprintf(stderr, " point %g %g %g\n", x1, y1, z1);
    GOMA_EH(GOMA_ERROR, "Fatal Error");
  }

  /*evaluate function             */

  switch (elem_order) {
  case 2:
    switch (dim) {
    case 2:
      for (i = 0; i < 3; i++) {
        k = 3 * i;
        bfi =
            (xi[0] - pt[i]) * (xi[0] - pt[i + 1]) / ((pt[i + 2] - pt[i]) * (pt[i + 2] - pt[i + 1]));
        for (j = 0; j < 3; j++) {
          bfj = (xi[1] - pt[j]) * (xi[1] - pt[j + 1]) /
                ((pt[j + 2] - pt[j]) * (pt[j + 2] - pt[j + 1]));
          phi[k + j] = bfi * bfj;
        }
      }
      break;
    case 3:
      for (i = 0; i < 3; i++) {
        k = 3 * i;
        bfi =
            (xi[0] - pt[i]) * (xi[0] - pt[i + 1]) / ((pt[i + 2] - pt[i]) * (pt[i + 2] - pt[i + 1]));
        for (j = 0; j < 3; j++) {
          bfj = (xi[1] - pt[j]) * (xi[1] - pt[j + 1]) /
                ((pt[j + 2] - pt[j]) * (pt[j + 2] - pt[j + 1]));
          for (l = 0; l < 3; l++) {
            bfl = (xi[1] - pt[l]) * (xi[1] - pt[l + 1]) /
                  ((pt[l + 2] - pt[l]) * (pt[l + 2] - pt[l + 1]));
            phi[k + j + 9 * l] = bfi * bfj * bfl;
          }
        }
      }
      break;
    }
    break;

  case 1:
    switch (dim) {
    case 2:
      phi[0] = (1. - xi[0]) * (1. - xi[1]);
      phi[1] = (1. - xi[0]) * xi[1];
      phi[2] = xi[0] * (1. - xi[1]);
      phi[3] = xi[0] * xi[1];
      if (dfunc_dx != NULL) {
        phic[0] = xi[1] - 1.;
        phic[1] = -xi[1];
        phic[2] = 1. - xi[1];
        phic[3] = xi[1];
        phie[0] = xi[0] - 1.;
        phie[1] = 1. - xi[0];
        phie[2] = -xi[0];
        phie[3] = xi[0];
      }
      break;
    case 3:
      phi[0] = (1. - xi[0]) * (1. - xi[1]) * (1. - xi[2]);
      phi[1] = (1. - xi[0]) * xi[1] * (1. - xi[2]);
      phi[2] = xi[0] * (1. - xi[1]) * (1. - xi[2]);
      phi[3] = xi[0] * xi[1] * (1. - xi[2]);
      phi[4] = (1. - xi[0]) * (1. - xi[1]) * xi[2];
      phi[5] = (1. - xi[0]) * xi[1] * xi[2];
      phi[6] = xi[0] * (1. - xi[1]) * xi[2];
      phi[7] = xi[0] * xi[1] * xi[2];
      if (dfunc_dx != NULL) {
        phic[0] = -(1. - xi[1]) * (1. - xi[2]);
        phic[1] = -xi[1] * (1. - xi[2]);
        phic[2] = (1. - xi[1]) * (1. - xi[2]);
        phic[3] = xi[1] * (1. - xi[2]);
        phic[4] = (xi[1] - 1.) * xi[2];
        phic[5] = -xi[1] * xi[2];
        phic[6] = (1. - xi[1]) * xi[2];
        phic[7] = xi[1] * xi[2];
        phie[0] = -(1. - xi[0]) * (1. - xi[2]);
        phie[1] = (1. - xi[0]) * (1. - xi[2]);
        phie[2] = -xi[0] * (1. - xi[2]);
        phie[3] = xi[0] * (1. - xi[2]);
        phie[4] = -(1. - xi[0]) * xi[2];
        phie[5] = (1. - xi[0]) * xi[2];
        phie[6] = -xi[0] * xi[2];
        phie[7] = xi[0] * xi[2];
        phig[0] = -(1. - xi[0]) * (1. - xi[1]);
        phig[1] = -(1. - xi[0]) * xi[1];
        phig[2] = -xi[0] * (1. - xi[1]);
        phig[3] = -xi[0] * xi[1];
        phig[4] = (1. - xi[0]) * (1. - xi[1]);
        phig[5] = (1. - xi[0]) * xi[1];
        phig[6] = xi[0] * (1. - xi[1]);
        phig[7] = xi[0] * xi[1];
      }
      break;
    }
    break;

  default:
    fprintf(stderr, "\n Unknown element order - quad_isomap\n");
    GOMA_EH(GOMA_ERROR, "Fatal Error");
  }

  for (i = 0; i < elem_nodes; i++) {
    pp[i] = p[itp[i]];
  }
  pvalue = 0.0;
  pc = 0.;
  pe = 0.;
  pg = 0.;
  for (i = 0; i < elem_nodes; i++) {
    pvalue += pp[i] * phi[i];
    if (dfunc_dx != NULL) {
      switch (dim) {
      case 3:
        pg += pp[i] * phig[i];
        /* fall through */
      case 2:
        pc += pp[i] * phic[i];
        pe += pp[i] * phie[i];
        break;
      }
    }
  }
  if (dfunc_dx != NULL) {
    switch (dim) {
    case 2:
      dfunc_dx[0] = (jac[1][1] * pc - jac[0][1] * pe) * detjti;
      dfunc_dx[1] = (-jac[1][0] * pc + jac[0][0] * pe) * detjti;
      break;
    case 3:
      dfunc_dx[0] = detjti * (pc * (jac[2][2] * jac[1][1] - jac[2][1] * jac[1][2]) +
                              pe * (jac[2][1] * jac[0][2] - jac[2][2] * jac[0][1]) +
                              pg * (jac[1][2] * jac[0][1] - jac[1][1] * jac[0][2]));
      dfunc_dx[1] = detjti * (pc * (jac[2][0] * jac[1][2] - jac[2][2] * jac[1][0]) +
                              pe * (jac[2][2] * jac[0][0] - jac[2][0] * jac[0][2]) +
                              pg * (jac[1][0] * jac[0][2] - jac[1][2] * jac[0][0]));
      dfunc_dx[2] = detjti * (pc * (jac[2][1] * jac[1][0] - jac[2][0] * jac[1][1]) +
                              pe * (jac[2][0] * jac[0][1] - jac[2][1] * jac[0][0]) +
                              pg * (jac[1][1] * jac[0][0] - jac[1][0] * jac[0][1]));
      break;
    }
  }
  return (pvalue);
}
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/

static double scalar_fv_fill_altmatrl(double **base, int lvdesc, int num_dofs, int var_type)

/**************************************************************************
 *   scalar_fv_fill_altmatr:
 *
 *      Kernal operation for calculating the value of a
 *       scalar field variable at the quadrature point from
 *      the basis functions.
 *      Note for the optimal results, this function should
 *      be inlined during the compile.
 **************************************************************************/
{
  int ln, idof, lvdof, lvdof_active;
  double *phi_ptr = bf[var_type]->phi, sum = 0.0;
  for (idof = 0; idof < num_dofs; idof++) {
    /*
     *  Find the local variable dof for the degree of freedom
     *  that actually belongs to another material than the one
     *  we are currently in.
     */
    lvdof = ei[pg->imtrx]->Lvdesc_to_lvdof[lvdesc][idof];
    /*
     *  Find the local node number of this degree of freedom
     */
    ln = ei[pg->imtrx]->dof_list[var_type][lvdof];
    /*
     *  Find the local variable type degree of freedom
     *  active in this element for this local node. This
     *  degree of freedom belongs to this material, and thus
     *  has a nonzero basis function associated with it.
     */
    lvdof_active = ei[pg->imtrx]->ln_to_dof[var_type][ln];
    /*
     *  Add it to the sum to get the value of the variable.
     */
    sum += *(base[lvdof]) * phi_ptr[lvdof_active];
  }
  return sum;
}
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/
double scalar_fv_fill_adjmatrl(double **base, int lvdesc, int num_dofs, int var_type)

/**************************************************************************
 *   scalar_fv_fill_adjmatr:
 *
 *      Identical to scalar_fv_fill_altmatrl but with phi_ptr array indexed
 *      by local node number.  This probably the way it should be done for interpolation
 *      orders that assign only one dof per node.  This mod was necessary to get it
 *      to work for look across BC's
 *      tabaer March 2004
 *
 **************************************************************************/
{
  int ln, idof, lvdof;
  double *phi_ptr = bf[var_type]->phi, sum = 0.0;
  for (idof = 0; idof < num_dofs; idof++) {
    /*
     *  Find the local variable dof for the degree of freedom
     *  that actually belongs to another material than the one
     *  we are currently in.
     */
    lvdof = ei[pg->imtrx]->Lvdesc_to_lvdof[lvdesc][idof];
    /*
     *  Find the local node number of this degree of freedom
     */
    ln = ei[pg->imtrx]->dof_list[var_type][lvdof];

    /*
     *  Add it to the sum to get the value of the variable.
     */
    sum += *(base[lvdof]) * phi_ptr[ln];
  }
  return sum;
}
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/

void load_matrl_statevector(MATRL_PROP_STRUCT *mp_local)

/**************************************************************************
 * load_matrl_statevector():
 *
 *  load up values of all relevant state variables at the
 *  current surface quad pt for materials that are not associated with the
 *  current element.
 *
 * input:
 * ------
 *	We assume that the state vector for the material of the
 *      current element has already been filled in. We use its values as
 *      defaults.
 *
 * output:
 * -------
 *       mp_local->StateVector[] is filled in as well, for
 *                 pertinent entries that make up the specification of 2q
 *                 the state of the material.
 *
 * HKM -> NOTE , probably not complete in terms of covering all
 *        needed state variable entries.
 *************************************************************************/
{
  int lvdesc, var_type, num_dofs, k, do_last_species = FALSE;
  double *sv = mp_local->StateVector;
  double *sv_mp = mp->StateVector, sp_sum = 0.0, rho;
  VARIABLE_DESCRIPTION_STRUCT *vd;
  int matID = mp_local->MatID;

  sv[TEMPERATURE] = sv_mp[TEMPERATURE];
  sv[FILL] = sv_mp[FILL];
  sv[VOLTAGE] = sv_mp[VOLTAGE];
  sv[PRESSURE] = sv_mp[PRESSURE];
  for (k = 0; k < mp->Num_Species; k++) {
    sv[SPECIES_UNK_0 + k] = sv_mp[SPECIES_UNK_0 + k];
  }

  for (lvdesc = 0; lvdesc < ei[pg->imtrx]->Num_Lvdesc; lvdesc++) {
    if (ei[pg->imtrx]->Lvdesc_to_MatID[lvdesc] == matID) {

      var_type = ei[pg->imtrx]->Lvdesc_to_Var_Type[lvdesc];
      num_dofs = ei[pg->imtrx]->Lvdesc_Numdof[lvdesc];

      if (var_type == MASS_FRACTION) {
        vd = ei[pg->imtrx]->Lvdesc_vd_ptr[lvdesc];
        k = vd->Subvar_Index;
        sv[SPECIES_UNK_0 + k] = scalar_fv_fill_altmatrl(esp->c[k], lvdesc, num_dofs, var_type);
        sp_sum += sv[SPECIES_UNK_0 + k];
        mp_local->StateVector_speciesVT = upd->Species_Var_Type;
        if (mp_local->Num_Species > mp_local->Num_Species_Eqn) {
          do_last_species = TRUE;
        }
      } else if (var_type == TEMPERATURE) {
        sv[var_type] = scalar_fv_fill_altmatrl(esp->T, lvdesc, num_dofs, var_type);
      } else if (var_type >= VELOCITY1 && var_type <= VELOCITY3) {
        k = var_type - VELOCITY1;
        sv[var_type] = scalar_fv_fill_altmatrl(esp->v[k], lvdesc, num_dofs, var_type);
      } else if (var_type == PRESSURE) {
        sv[var_type] = scalar_fv_fill_altmatrl(esp->P, lvdesc, num_dofs, var_type);
      } else if (var_type == FILL) {
        sv[var_type] = scalar_fv_fill_altmatrl(esp->F, lvdesc, num_dofs, var_type);
      } else if (var_type == LIGHT_INTP) {
        sv[var_type] = scalar_fv_fill_altmatrl(esp->poynt[0], lvdesc, num_dofs, var_type);
      } else if (var_type == LIGHT_INTM) {
        sv[var_type] = scalar_fv_fill_altmatrl(esp->poynt[1], lvdesc, num_dofs, var_type);
      } else {
        GOMA_EH(GOMA_ERROR, "Unimplemented");
      }
    }
  }
  if (do_last_species) {
    var_type = SPECIES_UNK_0 + mp_local->Num_Species - 1;
    switch (mp_local->Species_Var_Type) {
    case SPECIES_UNDEFINED_FORM:
    case SPECIES_MOLE_FRACTION:
    case SPECIES_MASS_FRACTION:
    case SPECIES_VOL_FRACTION:
      sv[var_type] = 1.0 - sp_sum;
      break;
    case SPECIES_CONCENTRATION:
      switch (mp_local->DensityModel) {
      case DENSITY_CONSTANT_LAST_CONC:
        sv[var_type] = mp_local->u_density[0];
        break;
      default:
        rho = calc_concentration(mp_local, FALSE, NULL);
        sv[var_type] = rho - sp_sum;
        break;
      }
      break;
    case SPECIES_DENSITY:
      /* Note, this won't work for time dependent densities - RRR */
      rho = calc_density(mp_local, FALSE, NULL, 0.0);
      sv[var_type] = rho - sp_sum;
      break;
    case SPECIES_CAP_PRESSURE:
      break;
    default:
      break;
    }
  }
}
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/

/*******************************************************************************/

/*
 * apply_distributed_sources
 *
 *    This function applies source terms to the conservation equation as narrow-banded
 *    source terms centered on the zero level set curve.
 *
 */

int apply_distributed_sources(int elem,
                              double width,
                              double x[],
                              Exo_DB *exo,
                              double dt, /* current time step size                       */
                              double theta,
                              double time,
                              const PG_DATA *pg_data,
                              int oAC,     /* Flag indicating calling function */
                              double *gAC, /* Augmenting condition arrays */
                              double **bAC,
                              double **cAC) {

  int ip, ip_total, ielem_type = ei[pg->imtrx]->ielem_type, err;
  struct LS_Embedded_BC *bcref;
  double xi[DIM];
  int ipass, num_elem_passes = 1;
  double ls_F[MDE], ad_wt[10]; /* adaptive integration weights  */
  int i;

  if (pd->e[pg->imtrx][R_PHASE1] && !pd->e[pg->imtrx][R_FILL]) {
    ls = pfd->ls[0];
  }

  if (xfem != NULL) {
    num_elem_passes = 2;
  } else {
    ls->Elem_Sign = 0;
  }

  if (ls->Integration_Depth > 0) {
    ip_total = Subgrid_Int.ip_total;
  } else {
    ip_total = elem_info(NQUAD, ielem_type);
  }
  if (ls->AdaptIntegration) {
    for (i = 0; i < ei[pg->imtrx]->num_local_nodes; i++) {
      ls_F[i] = *esp->F[i];
    }
    i = adaptive_weight(ad_wt, ip_total, ei[pg->imtrx]->ielem_dim, ls_F, ls->Length_Scale, 3,
                        ielem_type);
    GOMA_WH(i, "problem with adaptive weight routine");
  }

  for (ip = 0; ip < ip_total; ip++) {

    if (ls->Integration_Depth > 0) {
      xi[0] = Subgrid_Int.s[ip][0];
      xi[1] = Subgrid_Int.s[ip][1];
      xi[2] = Subgrid_Int.s[ip][2];
      fv->wt = Subgrid_Int.wt[ip];
    } else {
      find_stu(ip, ielem_type, &(xi[0]), &(xi[1]), &(xi[2]));
      if (ls->AdaptIntegration) {
        fv->wt = ad_wt[ip_total - 1 - ip];
      } else {
        fv->wt = Gq_weight(ip, ielem_type);
      }
    }

    for (ipass = 0; ipass < num_elem_passes; ipass++) {
      if (num_elem_passes == 2)
        ls->Elem_Sign = -1 + 2 * ipass;

      err = load_basis_functions(xi, bfd);
      GOMA_EH(err, "problem from load_basis_functions");

      err = beer_belly();
      GOMA_EH(err, "beer_belly");

      err = load_fv();
      GOMA_EH(err, "load_fv");

      err = load_bf_grad();
      GOMA_EH(err, "load_bf_grad");

      if (pd->e[pg->imtrx][R_MESH1]) {
        err = load_bf_mesh_derivs();
        GOMA_EH(err, "load_bf_mesh_derivs");
      }

      err = load_fv_grads();
      GOMA_EH(err, "load_fv_grads");

      if (pd->e[pg->imtrx][R_MESH1]) {
        err = load_fv_mesh_derivs(1);
        GOMA_EH(err, "load_fv_mesh_derivs");
      }

      load_lsi(ls->Length_Scale);
      load_lsi_derivs();

      /* Now where are the material properties in all of this.  Well, here we propagate the
       * unfortunate practice of computing the material properties in the assembly routines
       * themselves.
       */

      bcref = ls->embedded_bc;

      while (bcref != NULL) {
        BOUNDARY_CONDITION_STRUCT *bc;

        bc = BC_Types + bcref->bc_input_id;

        if (num_elem_passes == 2) {
          if ((bc->BC_ID != 0 && ls->Elem_Sign != bc->BC_ID) || (bc->BC_ID == 0 && ipass > 0)) {
            bcref = bcref->next;
            continue;
          }
        }

        switch (bc->BC_Name) {
        case LS_STRESS_JUMP_BC:
          assemble_ls_stress_jump(bc->BC_Data_Float[0], bc->BC_Data_Float[1], bc->BC_Data_Int[0]);
          break;
        case LS_CAPILLARY_BC:
          assemble_csf_tensor();
          break;
        case LS_CAP_HYSING_BC:
          assemble_cap_hysing(dt, bc->BC_Data_Float[0]);
          break;
        case LS_CAP_DENNER_DIFF_BC:
          if (pd->gv[R_NORMAL1]) {
            assemble_cap_denner_diffusion_n(dt, bc->BC_Data_Float[0]);
          } else {
            assemble_cap_denner_diffusion(dt, bc->BC_Data_Float[0]);
          }
          break;
        case LS_CAP_CURVE_BC:
          if (pd->e[pg->imtrx][R_NORMAL1])
            assemble_curvature_with_normals_source();
          else
            assemble_curvature_source();
          break;
        case LS_CAP_DIV_N_BC:
          assemble_div_n_source();
          break;
        case LS_CAP_DIV_S_N_BC:
          assemble_div_s_n_source();
          break;
        case LS_Q_BC:
          assemble_q_source(bc->BC_Data_Float[0]);
          break;
        case LS_QLASER_BC:
          assemble_qlaser_source(bc->u_BC, time);
          break;
        case LS_QVAPOR_BC:
          assemble_qvapor_source(bc->u_BC);
          break;
        case LS_QRAD_BC:
          assemble_qrad_source(bc->BC_Data_Float[0], bc->BC_Data_Float[1], bc->BC_Data_Float[2],
                               bc->BC_Data_Float[3]);
          break;
        case LS_T_BC:
          assemble_t_source(bc->BC_Data_Float[0], time);
          break;
        case LS_LATENT_HEAT_BC:
          assemble_ls_latent_heat_source(bc->BC_Data_Float[0], bc->BC_Data_Float[1], dt, theta,
                                         time, bcref->bc_input_id, BC_Types);
          break;
        case LS_YFLUX_BC:
          assemble_ls_yflux_source(bc->BC_Data_Int[0], bc->BC_Data_Float[0], bc->BC_Data_Float[1],
                                   dt, theta, time, bcref->bc_input_id, BC_Types);
          break;
        case LS_CONT_T_BC:
          assemble_cont_t_source(xi);
          break;
        case LS_CONT_VEL_BC:
          assemble_cont_vel_source(xi, exo);
          break;
        case LS_FLOW_PRESSURE_BC:
          assemble_p_source(bc->BC_Data_Float[0], bc->BC_Data_Int[0]);
          break;
        case LS_ACOUSTIC_SOURCE_BC:
          assemble_ars_source(bc->BC_Data_Float[0], bc->BC_Data_Float[1]);
          break;
        case LS_RECOIL_PRESSURE_BC:
          assemble_precoil_source(bc->BC_Data_Float);
          break;
        case LS_U_BC:
        case LS_V_BC:
        case LS_W_BC:
          assemble_uvw_source(bc->desc->equation, bc->BC_Data_Float[0]);
          break;
        case BAAIJENS_FLUID_SOLID_BC:
          assemble_LM_source(xi, oAC, gAC, bAC, cAC, x, exo);
          break;
        case LS_EXTV_FLUID_SIC_BC:
          assemble_interface_extension_velocity_sic(bc->BC_Data_Int[0]);
          break;
        case LS_EXTV_KINEMATIC_BC:
          assemble_extv_kinematic(theta, dt, time, bcref->bc_input_id, BC_Types);
          break;
        case LS_EIK_KINEMATIC_BC:
          assemble_eik_kinematic(theta, dt, time, bcref->bc_input_id, BC_Types);
          break;
        case PF_CAPILLARY_BC:
          assemble_pf_capillary(&(bc->BC_Data_Float[0]));

          break;
        default:
          break;
        }
        bcref = bcref->next;
      }
      /* equation path dependence terms */
      if (!ls->Ignore_F_deps) {
        if (pd->e[pg->imtrx][R_FILL] && ipass == 0) {
          if (tran->Fill_Equation == FILL_EQN_EIKONAL) {
            assemble_fill_path_dependence();
          }
        }
        if (pd->e[pg->imtrx][R_MOMENTUM1]) {
          err = assemble_momentum_path_dependence(time, theta, dt, pg_data);
          GOMA_EH(err, "assemble_momentum_path_dependence");
        }
        if (pd->e[pg->imtrx][R_PRESSURE]) {
#ifndef DARWIN_HACK
          err = assemble_continuity_path_dependence(time, theta, dt, pg_data);
          GOMA_EH(err, "assemble_continuity_path_dependence");
#endif
        }
        if (pd->e[pg->imtrx][R_ENERGY]) {
          assemble_energy_path_dependence(time, theta, dt, pg_data);
        }
      }
    }
  }

  /* leave the place as tidy as it was before you came */
  ls->Elem_Sign = 0;

  return (2);
}

int assemble_pf_capillary(double *pf_surf_tens) {
  int i, j, a, p, q, ii, ledof;
  int eqn, peqn, var, pvar;

  int pf;

  struct Basis_Functions *bfm;
  dbl(*grad_phi_i_e_a)[DIM] = NULL;

  double wt, det_J, h3;

  double csf[DIM][DIM];
  double d_csf_dF[DIM][DIM][MDE];
  double d_csf_dX[DIM][DIM][DIM][MDE];

  double source;
  struct Level_Set_Data *ls_save = ls;

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

  for (pf = 0; pf < pfd->num_phase_funcs; pf++) {
    memset(csf, 0, sizeof(double) * DIM * DIM);
    memset(d_csf_dF, 0, sizeof(double) * DIM * DIM * MDE);

    ls = pfd->ls[pf];

    continuous_surface_tension(mp->surface_tension * pf_surf_tens[pf], csf, d_csf_dF, d_csf_dX);

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

            source = 0.;

            grad_phi_i_e_a = bfm->grad_phi_e[i][a];

            for (p = 0; p < VIM; p++) {
              for (q = 0; q < VIM; q++) {
                source += grad_phi_i_e_a[p][q] * csf[q][p];
              }
            }

            source *= -det_J * wt;

            source *= h3;

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

            grad_phi_i_e_a = bfm->grad_phi_e[i][a];

            /* J_m_pF
             */
            var = PHASE1 + pf;
            if (pd->v[pg->imtrx][var]) {
              pvar = upd->vp[pg->imtrx][var];

              for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                source = 0.;

                for (p = 0; p < VIM; p++) {
                  for (q = 0; q < VIM; q++) {
                    source += grad_phi_i_e_a[p][q] * d_csf_dF[q][p][j];
                  }
                }

                source *= -det_J * wt * h3;
                source *= pd->etm[pg->imtrx][eqn][LOG2_SOURCE];

                lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
              }
            }
          }
        }
      }
    }
  }

  ls = ls_save;

  return (1);
}

int assemble_ls_stress_jump(double viscosity_scale, double stress_scale, int heaviside_type) {
  if (!pd->e[pg->imtrx][R_MOMENTUM1]) {
    return (0);
  }

  double wt = fv->wt;
  double h3 = fv->h3;
  double det_J = bf[R_MOMENTUM1]->detJ;
  if (ls->on_sharp_surf) /* sharp interface */
  {
    det_J = fv->sdet;
  }
  double d_area = wt * det_J * h3;

  // hard code jump for test
  double pos_mup = ve[0]->gn->pos_ls_mup;
  double neg_mup = ve[0]->gn->mu0;
  double ls_viscosity_jump = fabs(mp->viscosity - mp->mp2nd->viscosity);
  double Heaviside = 1;

  switch (heaviside_type) {
  case -1:
    Heaviside = 1 - lsi->H;
    break;
  case 1:
    Heaviside = lsi->H;
    break;
  default:
    Heaviside = 1;
    break;
  }

  load_lsi(ls->Length_Scale);
  switch (ls->ghost_stress) {
  case LS_POSITIVE:
    ls_viscosity_jump += neg_mup;
    break;
  case LS_NEGATIVE:
    ls_viscosity_jump += pos_mup;
    break;
  case LS_OFF:
  default:
    ls_viscosity_jump += fabs(pos_mup - neg_mup);
    break;
  }
  /*
   * Residuals ____________________________________________________________________________
   */
  if (af->Assemble_Residual) {

    for (int a = 0; a < WIM; a++) {
      int eqn = R_MOMENTUM1 + a;
      int peqn = upd->ep[pg->imtrx][eqn];
      for (int i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {

        int ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][i];

        if (ei[pg->imtrx]->active_interp_ledof[ledof]) {

          int ii = ei[pg->imtrx]->lvdof_to_row_lvdof[eqn][i];

          double n_gv_n = 0;
          double n_tau_n = 0;
          for (int p = 0; p < VIM; p++) {
            for (int q = 0; q < VIM; q++) {
              n_gv_n += lsi->normal[p] * (fv->grad_v[p][q] + fv->grad_v[q][p]) * lsi->normal[q];
              n_tau_n += lsi->normal[p] * fv->S[0][p][q] * lsi->normal[q];
            }
          }
          double source =
              lsi->delta * lsi->normal[a] *
              (viscosity_scale * ls_viscosity_jump * n_gv_n + Heaviside * stress_scale * n_tau_n) *
              bf[eqn]->phi[i];
          source *= d_area;
          lec->R[LEC_R_INDEX(peqn, ii)] += source;
        }
      }
    }
  }

  /*
   * Yacobian terms...
   */

  if (af->Assemble_Jacobian) {
    for (int a = 0; a < WIM; a++) {
      int eqn = R_MOMENTUM1 + a;
      int peqn = upd->ep[pg->imtrx][eqn];

      for (int i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {

        int ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][i];

        if (ei[pg->imtrx]->active_interp_ledof[ledof]) {

          int ii = ei[pg->imtrx]->lvdof_to_row_lvdof[eqn][i];

          double n_tau_n = 0;
          double n_gv_n = 0;
          for (int p = 0; p < VIM; p++) {
            for (int q = 0; q < VIM; q++) {
              n_gv_n += lsi->normal[q] * (fv->grad_v[p][q] + fv->grad_v[q][p]) * lsi->normal[p];
              n_tau_n += lsi->normal[q] * fv->S[0][q][p] * lsi->normal[p];
            }
          }
          for (int b = 0; b < WIM; b++) {
            int var = R_MOMENTUM1 + b;
            if (pd->v[pg->imtrx][var]) {
              int pvar = upd->vp[pg->imtrx][var];
              for (int j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                double n_gv_n_dv = 0;
                for (int p = 0; p < VIM; p++) {
                  for (int q = 0; q < VIM; q++) {
                    n_gv_n_dv +=
                        lsi->normal[p] *
                        (bf[var]->grad_phi_e[j][b][p][q] + bf[var]->grad_phi_e[j][b][q][p]) *
                        lsi->normal[q];
                  }
                }
                double source = lsi->normal[a] * lsi->delta * viscosity_scale * ls_viscosity_jump *
                                n_gv_n_dv * bf[eqn]->phi[i];
                source *= d_area;

                lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
              }
            }
          }

          /* J_m_F
           */
          int var = LS;
          if (pd->v[pg->imtrx][var]) {
            int pvar = upd->vp[pg->imtrx][var];

            for (int j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              double n_tau_n_dF = 0;
              double n_gv_n_dF = 0;
              for (int p = 0; p < VIM; p++) {
                for (int q = 0; q < VIM; q++) {
                  n_tau_n_dF += lsi->d_normal_dF[q][j] * fv->S[0][q][p] * lsi->normal[p];
                  n_tau_n_dF += lsi->normal[q] * fv->S[0][q][p] * lsi->d_normal_dF[p][j];
                  n_gv_n_dF += lsi->normal[q] * fv->grad_v[q][p] * lsi->d_normal_dF[p][j];
                }
              }
              double source =
                  lsi->d_delta_dF[j] * lsi->normal[a] *
                  (viscosity_scale * ls_viscosity_jump * n_gv_n + stress_scale * n_tau_n) *
                  bf[eqn]->phi[i];
              source += lsi->delta * lsi->d_normal_dF[a][j] *
                        (viscosity_scale * ls_viscosity_jump * n_gv_n + stress_scale * n_tau_n) *
                        bf[eqn]->phi[i];
              source += lsi->delta * lsi->normal[a] *
                        (viscosity_scale * n_gv_n_dF + stress_scale * n_tau_n_dF) * bf[eqn]->phi[i];
              source *= d_area;

              lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
            }
          }
        }
      }
    }
  }

  return (1);
}

int assemble_curvature_with_normals_source(void) {
  int i, j, a, b, ii, ledof;
  int eqn, peqn, var, pvar;
  int dim;

  struct Basis_Functions *bfm;

  double wt, det_J, h3, phi_i, phi_j;

  double source;

  /* Bail out with an error if CURVATURE equation not define
   */

  if (!pd->e[pg->imtrx][R_CURVATURE]) {
    GOMA_EH(GOMA_ERROR,
            "Error: Level set curvature equation needs to be activated to use LS_CAP_CURVE\n");
  }

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

  dim = pd->Num_Dim;

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

          source = mp->surface_tension * fv->H * fv->n[a] * lsi->delta;

          source *= wt * phi_i * det_J * h3;
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
           * J_m_H
           */
          var = CURVATURE;

          if (pd->v[pg->imtrx][var]) {
            pvar = upd->vp[pg->imtrx][var];

            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              phi_j = bf[var]->phi[j];

              source = mp->surface_tension * phi_j * fv->n[a] * lsi->delta;

              source *= phi_i * wt * det_J * h3;
              source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

              lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
            }
          }

          /*
           * J_m_n
           */

          var = NORMAL1 + a;

          if (pd->v[pg->imtrx][var]) {
            pvar = upd->vp[pg->imtrx][var];

            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              phi_j = bf[var]->phi[j];

              source = mp->surface_tension * fv->H * phi_j * lsi->delta;
              source *= phi_i * wt * det_J * h3;
              source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

              lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
            }
          }
          /*
           * J_m_F
           */

          var = LS;

          if (pd->v[pg->imtrx][var]) {
            pvar = upd->vp[pg->imtrx][var];

            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              source = mp->surface_tension * fv->H * fv->n[a] * lsi->d_delta_dF[j];
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
                source = mp->surface_tension * fv->H * fv->n[a] * lsi->delta;

                source *= (bf[var]->d_det_J_dm[b][j] * h3 + fv->dh3dmesh[b][j] * det_J);
                source *= wt * phi_i;
                source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
              }
            }
          }
        }
      }
    }
  }
  return (1);
}

int assemble_curvature_source(void) {
  int i, j, a, b, ii, ledof;
  int eqn, peqn, var, pvar;
  int dim;

  struct Basis_Functions *bfm;

  double wt, det_J, h3, phi_i, phi_j;

  double source;
  double f[DIM], dfdX[DIM][DIM][MDE], dfdF[DIM][MDE], dfdH[DIM][MDE];

  /* Bail out with an error if CURVATURE equation not define
   */

  if (!pd->e[pg->imtrx][R_CURVATURE]) {
    GOMA_EH(GOMA_ERROR,
            "Error: Level set curvature equation needs to be activated to use LS_CAP_CURVE\n");
  }

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

  dim = pd->Num_Dim;

  memset(dfdH, 0, DIM * MDE * sizeof(double));
  memset(dfdX, 0, DIM * DIM * MDE * sizeof(double));
  memset(dfdF, 0, DIM * MDE * sizeof(double));
  memset(f, 0, DIM * sizeof(double));

  curvature_momentum_source(f, dfdX, dfdF, dfdH);

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

          source = phi_i * f[a];

          source *= det_J * wt;

          source *= h3;

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

          source = phi_i * f[a];

          source *= det_J * wt;

          source *= h3;

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

          /* J_m_F
           */
          var = LS;
          if (pd->v[pg->imtrx][var]) {
            pvar = upd->vp[pg->imtrx][var];

            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              source = phi_i * dfdF[a][j];

              source *= det_J * wt * h3;
              source *= pd->etm[pg->imtrx][eqn][LOG2_SOURCE];

              lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
            }
          }
          /*
           * J_m_H
           */
          var = CURVATURE;

          if (pd->v[pg->imtrx][var]) {
            pvar = upd->vp[pg->imtrx][var];

            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              phi_j = bf[var]->phi[j];
              source = phi_i * dfdH[a][j];

              source *= det_J * wt * h3;
              source *= pd->etm[pg->imtrx][eqn][LOG2_SOURCE];

              lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
            }
          }

          /* J_m_X
           */

          for (b = 0; b < dim; b++) {
            var = R_MESH1 + b;

            if (pd->v[pg->imtrx][var]) {
              pvar = upd->vp[pg->imtrx][var];

              for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                phi_j = bf[var]->phi[j];
                source = phi_i * dfdX[a][b][j];

                source *= det_J * wt * h3;
                source *= pd->etm[pg->imtrx][eqn][LOG2_SOURCE];

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

int assemble_q_source(double flux) {
  int i, j, ii, ledof;
  int eqn, peqn, var, pvar;

  struct Basis_Functions *bfm;
  double wt, det_J, h3, phi_i, phi_j;

  double source;

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

  /* finite difference calculation of path dependencies for
     subelement integration
   */
  if (ls->CalcSurfDependencies) {
    peqn = upd->ep[pg->imtrx][eqn];
    bfm = bf[eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][i];

      if (ei[pg->imtrx]->active_interp_ledof[ledof]) {
        ii = ei[pg->imtrx]->lvdof_to_row_lvdof[eqn][i];

        phi_i = bfm->phi[i];

        source = phi_i * flux * lsi->delta;

        source *= det_J * wt * h3;
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
    return (0);
  }

  /*
   * Wesiduals
   * ________________________________________________________________________________
   */
  if (af->Assemble_Residual) {
    peqn = upd->ep[pg->imtrx][eqn];
    bfm = bf[eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][i];

      if (ei[pg->imtrx]->active_interp_ledof[ledof]) {
        ii = ei[pg->imtrx]->lvdof_to_row_lvdof[eqn][i];

        phi_i = bfm->phi[i];

        source = phi_i * flux * lsi->delta;

        source *= det_J * wt * h3;
        source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

        lec->R[LEC_R_INDEX(peqn, ii)] += source;
      }
    }
  }
  /*
   * Yacobian terms...
   */

  if (af->Assemble_Jacobian) {
    peqn = upd->ep[pg->imtrx][eqn];
    bfm = bf[eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][i];

      if (ei[pg->imtrx]->active_interp_ledof[ledof]) {
        ii = ei[pg->imtrx]->lvdof_to_row_lvdof[eqn][i];

        phi_i = bfm->phi[i];

        /* diffuse interface version of path dependence integral */
        /*
         * J_e_F
         */
        var = FILL;

        if (pd->v[pg->imtrx][var]) {
          pvar = upd->vp[pg->imtrx][var];

          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            source = phi_i * flux * lsi->d_delta_dF[j];

            source *= det_J * wt * h3;
            source *= pd->etm[pg->imtrx][eqn][LOG2_SOURCE];

            lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
          }
        }
      }
    }
  }

  return (1);
}

#if 1
int assemble_t_source(double T, double time) {
  int i, j, ii, ledof, p;
  int eqn, peqn, var, pvar;
  int xfem_active, extended_dof, base_interp, base_dof;
  int incomplete, other_side;
  int apply_NOBC[MDE], apply_SIC[MDE];

  struct Basis_Functions *bfm;
  double wt, det_J, h3, phi_i, phi_j;
  double sign = ls->Elem_Sign;

  double F_i;

  double q[DIM];
  HEAT_FLUX_DEPENDENCE_STRUCT d_q_struct;
  HEAT_FLUX_DEPENDENCE_STRUCT *d_q = &d_q_struct;

  double source;

  eqn = R_ENERGY;
  if (!pd->e[pg->imtrx][eqn]) {
    return (0);
  }

  /* For LS "Dirichlet" conditions, die if attempted on
     non-discontinuous interpolations */
  if (pd->i[pg->imtrx][eqn] != I_Q1_XV && pd->i[pg->imtrx][eqn] != I_Q2_XV &&
      pd->i[pg->imtrx][eqn] != I_Q1_GP && pd->i[pg->imtrx][eqn] != I_Q2_GP &&
      pd->i[pg->imtrx][eqn] != I_Q1_GN && pd->i[pg->imtrx][eqn] != I_Q2_GN &&
      pd->i[pg->imtrx][eqn] != I_Q1_G && pd->i[pg->imtrx][eqn] != I_Q2_G) {
    GOMA_EH(GOMA_ERROR,
            "LS_T requires discontinuous temperature enrichment (Q1_XV, Q2_XV, Q1_G, Q2_G)");
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
   * non-xfem dofs will receive nobc and xfem
   * dofs will receive SIC
   */

  /* DRN: This was previously set up to add the penalty term everywhere.
     For the test cases I ran, it seemed to work great.  But when this was
     tried on the momentum eqns (uvw_source), it didn't work.  So the mixed
     nobc/SIC was invented.  So this now applied here too, even though it
     may be overkill.
   */

  heat_flux(q, d_q, time);

  for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
    ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][i];

    if (ei[pg->imtrx]->active_interp_ledof[ledof]) {
      xfem_dof_state(i, pd->i[pg->imtrx][eqn], ei[pg->imtrx]->ielem_shape, &xfem_active,
                     &extended_dof, &base_interp, &base_dof);
      incomplete = dof_incomplete(base_dof, ei[pg->imtrx]->ielem_type, base_interp,
                                  ei[pg->imtrx]->ielem_shape);
      F_i = dof_distance(eqn, i);
      other_side = sign_change(F_i, (double)ls->Elem_Sign);

      if (extended_dof) {
        if (other_side) {
          if (incomplete) {
            /* extended_dof&&other_side&&incomplete -> NOBC + SIC */
            apply_NOBC[i] = TRUE;
            apply_SIC[i] = TRUE;
          } else {
            /* extended_dof&&other_side&&!incomplete -> NOBC */
            apply_NOBC[i] = TRUE;
            apply_SIC[i] = FALSE;
          }
        } else {
          /* extended_dof&&!other_side -> do nothing */
          apply_NOBC[i] = FALSE;
          apply_SIC[i] = FALSE;
        }
      } else {
        if (other_side) {
          /* !extended_dof&&other_side -> do nothing */
          apply_NOBC[i] = FALSE;
          apply_SIC[i] = FALSE;
        } else {
          /* !extended_dof&&!other_side -> NOBC */
          apply_NOBC[i] = TRUE;
          apply_SIC[i] = FALSE;
        }
      }
    }
  }

#if 0
{
   double nobc_flux = 0., sic_flux = 0.;
   double r = sqrt(fv->x[0]*fv->x[0] + fv->x[1]*fv->x[1]);
   double angle = atan2(fv->x[1], fv->x[0]);
   int side = 0;
   for ( p=0; p<VIM; p++) nobc_flux += q[p] * lsi->normal[p] * sign;
   sic_flux = (T - fv->T) * BIG_PENALTY;
   if ( r > 0.3 ) side = 1;
   if ( ls->Elem_Sign == -1 ) 
   fprintf(stderr,"FLUX: side, angle, T, NOBC, SIC: %d %g %g %g %g\n",side,angle,fv->T,nobc_flux,sic_flux);
}
#endif

  /* finite difference calculation of path dependencies for
     subelement integration
   */
  if (ls->CalcSurfDependencies) {
    peqn = upd->ep[pg->imtrx][eqn];
    bfm = bf[eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][i];

      if (ei[pg->imtrx]->active_interp_ledof[ledof]) {
        ii = ei[pg->imtrx]->lvdof_to_row_lvdof[eqn][i];

        phi_i = bfm->phi[i];

        if (apply_NOBC[i]) {
          source = 0.;
          for (p = 0; p < VIM; p++)
            source += q[p] * lsi->normal[p];
          source *= phi_i * lsi->delta * sign;
          source *= det_J * wt * h3;
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
        if (apply_SIC[i]) {
          source = phi_i * lsi->delta * (T - fv->T);
          source *= det_J * wt * h3 * BIG_PENALTY;
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
    peqn = upd->ep[pg->imtrx][eqn];
    bfm = bf[eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][i];

      if (ei[pg->imtrx]->active_interp_ledof[ledof]) {
        ii = ei[pg->imtrx]->lvdof_to_row_lvdof[eqn][i];

        phi_i = bfm->phi[i];

        if (apply_NOBC[i]) {
          source = 0.;
          for (p = 0; p < VIM; p++)
            source += q[p] * lsi->normal[p];
          source *= phi_i * lsi->delta * sign;
          source *= det_J * wt * h3;
          source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

          lec->R[LEC_R_INDEX(peqn, ii)] += source;
        }
        if (apply_SIC[i]) {
          source = phi_i * lsi->delta * (T - fv->T);
          source *= det_J * wt * h3 * BIG_PENALTY;
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
    eqn = R_ENERGY;
    peqn = upd->ep[pg->imtrx][eqn];
    bfm = bf[eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {

      ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][i];

      if (ei[pg->imtrx]->active_interp_ledof[ledof]) {

        ii = ei[pg->imtrx]->lvdof_to_row_lvdof[eqn][i];

        phi_i = bfm->phi[i];

        if (apply_NOBC[i]) {
          /*
           * J_e_T
           */
          var = TEMPERATURE;

          if (pd->v[pg->imtrx][var]) {
            pvar = upd->vp[pg->imtrx][var];

            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              phi_j = bf[var]->phi[j];

              source = 0.;
              for (p = 0; p < VIM; p++)
                source += d_q->T[p][j] * lsi->normal[p];
              source *= phi_i * lsi->delta * sign;
              source *= det_J * wt * h3;
              source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

              lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
            }
          }
          /*
           * J_e_F
           */
          var = FILL;

          if (pd->v[pg->imtrx][var]) {
            pvar = upd->vp[pg->imtrx][var];

            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              phi_j = bf[var]->phi[j];

              source = 0.;
              for (p = 0; p < VIM; p++)
                source += d_q->F[p][j] * lsi->normal[p] * lsi->delta +
                          q[p] * lsi->d_normal_dF[p][j] * lsi->delta +
                          q[p] * lsi->normal[p] * lsi->d_delta_dF[j];

              source *= phi_i * sign;
              source *= det_J * wt * h3;
              source *= pd->etm[pg->imtrx][eqn][LOG2_SOURCE];

              lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
            }
          }
          /* NOTE: dependencies for X and C (and anything else that modifies q)
             need to be added here */
        }
        if (apply_SIC[i]) {
          /*
           * J_T_T
           */
          var = TEMPERATURE;

          if (pd->v[pg->imtrx][var]) {
            pvar = upd->vp[pg->imtrx][var];

            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              phi_j = bf[var]->phi[j];

              source = -phi_i * lsi->delta * phi_j;

              source *= det_J * wt * h3 * BIG_PENALTY;
              source *= pd->etm[pg->imtrx][eqn][LOG2_SOURCE];

              lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
            }
          }
        }
      }
    }
  }
  return (1);
}
#endif

int assemble_qlaser_source(const double p[], double time) {
  int i, j, ii, ledof;
  int eqn, peqn, var, pvar;

  struct Basis_Functions *bfm;
  double wt, det_J, h3, phi_i, phi_j, flux, qlaser;
  double absorp, absorp_base;
  double normal[3] = {0., 0., 0.};

  double source, d_laser_dx[MAX_PDIM];

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

  absorp_base = p[2];
  absorp = absorp_base;

  for (i = 0; i < pd->Num_Dim; i++)
    normal[i] = lsi->normal[i] * ls->Elem_Sign;

  memset(d_laser_dx, 0, sizeof(double) * MAX_PDIM);

  /* Lets CALCULATE LASER FLUX  RARR */
  /* WARNING! sending p as x until we find what method can be used for use_pth and level
   * sets */
  qlaser = calculate_laser_flux(p, time, p, normal, 0, 1, d_laser_dx);

  /* NO Derivatives for d_laser_dx used!!!
   * MUST FIX THIS LATER!!
   */

  flux = absorp * qlaser;
  /* END CALCULATE LASER FLUX */

  /* finite difference calculation of path dependencies for
     subelement integration
   */
  if (ls->CalcSurfDependencies) {
    peqn = upd->ep[pg->imtrx][eqn];
    bfm = bf[eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][i];

      if (ei[pg->imtrx]->active_interp_ledof[ledof]) {
        ii = ei[pg->imtrx]->lvdof_to_row_lvdof[eqn][i];

        phi_i = bfm->phi[i];

        source = phi_i * flux * lsi->delta;

        source *= det_J * wt;

        source *= h3;

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
    return (0);
  }

  /*
   * Wesiduals
   * ________________________________________________________________________________
   */
  if (af->Assemble_Residual) {
    peqn = upd->ep[pg->imtrx][eqn];
    bfm = bf[eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][i];

      if (ei[pg->imtrx]->active_interp_ledof[ledof]) {
        ii = ei[pg->imtrx]->lvdof_to_row_lvdof[eqn][i];

        phi_i = bfm->phi[i];

        source = phi_i * flux * lsi->delta;

        source *= det_J * wt;

        source *= h3;

        source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

        lec->R[LEC_R_INDEX(peqn, ii)] += source;
      }
    }
  }

  if (af->Assemble_Jacobian) {
    peqn = upd->ep[pg->imtrx][eqn];
    bfm = bf[eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][i];

      if (ei[pg->imtrx]->active_interp_ledof[ledof]) {
        ii = ei[pg->imtrx]->lvdof_to_row_lvdof[eqn][i];

        phi_i = bfm->phi[i];

        /*
         * J_T_T
         */
        var = TEMPERATURE;

        if (pd->v[pg->imtrx][var]) {
          pvar = upd->vp[pg->imtrx][var];

          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];

            source = phi_i * lsi->delta * (0.0) * phi_j;

            source *= det_J * wt * h3;
            source *= pd->etm[pg->imtrx][eqn][LOG2_SOURCE];

            lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
          }
        }
        /*
         * J_e_F
         */
        var = FILL;

        if (pd->v[pg->imtrx][var]) {
          pvar = upd->vp[pg->imtrx][var];

          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            source = phi_i * flux * lsi->d_delta_dF[j];

            source *= det_J * wt * h3;
            source *= pd->etm[pg->imtrx][eqn][LOG2_SOURCE];

            lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
          }
        }
      }
    }
  }
  return (1);
}

int assemble_qvapor_source(const double p[]) {
  int i, j, ii, ledof;
  int eqn, peqn, var, pvar;

  struct Basis_Functions *bfm;
  double wt, det_J, h3, phi_i, phi_j, flux;
  double qvaporloss, d_evap_loss = 0.0;
  double source;

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

  /***************************** EXECUTION BEGINS *******************************/
  /* Lets CALCULATE LASER VAPOR HEAT LOSS  RARR */
  qvaporloss = calculate_vapor_cool(p, &d_evap_loss, 0);

  flux = -qvaporloss;
  /* END CALCULATE Q VAPOR */

  /* finite difference calculation of path dependencies for
     subelement integration
   */
  if (ls->CalcSurfDependencies) {
    peqn = upd->ep[pg->imtrx][eqn];
    bfm = bf[eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][i];

      if (ei[pg->imtrx]->active_interp_ledof[ledof]) {
        ii = ei[pg->imtrx]->lvdof_to_row_lvdof[eqn][i];

        phi_i = bfm->phi[i];

        source = phi_i * flux * lsi->delta;

        source *= det_J * wt;

        source *= h3;

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
    return (0);
  }

  /*
   * Wesiduals
   * ________________________________________________________________________________
   */
  if (af->Assemble_Residual) {
    peqn = upd->ep[pg->imtrx][eqn];
    bfm = bf[eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][i];

      if (ei[pg->imtrx]->active_interp_ledof[ledof]) {
        ii = ei[pg->imtrx]->lvdof_to_row_lvdof[eqn][i];

        phi_i = bfm->phi[i];

        source = phi_i * flux * lsi->delta;

        source *= det_J * wt;

        source *= h3;

        source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

        lec->R[LEC_R_INDEX(peqn, ii)] += source;
      }
    }
  }

  if (af->Assemble_Jacobian) {
    peqn = upd->ep[pg->imtrx][eqn];
    bfm = bf[eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][i];

      if (ei[pg->imtrx]->active_interp_ledof[ledof]) {
        ii = ei[pg->imtrx]->lvdof_to_row_lvdof[eqn][i];

        phi_i = bfm->phi[i];

        /*
         * J_T_T
         */
        var = TEMPERATURE;

        if (pd->v[pg->imtrx][var]) {
          pvar = upd->vp[pg->imtrx][var];

          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];

            source = phi_i * lsi->delta * (-d_evap_loss) * phi_j;

            source *= det_J * wt * h3;
            source *= pd->etm[pg->imtrx][eqn][LOG2_SOURCE];

            lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
          }
        }
        /*
         * J_e_F
         */
        var = FILL;

        if (pd->v[pg->imtrx][var]) {
          pvar = upd->vp[pg->imtrx][var];

          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            source = phi_i * flux * lsi->d_delta_dF[j];

            source *= det_J * wt * h3;
            source *= pd->etm[pg->imtrx][eqn][LOG2_SOURCE];

            lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
          }
        }
      }
    }
  }
  return (1);
}

int assemble_qrad_source(double htc, double Tref, double emiss, double sigma) {
  int i, j, ii, ledof;
  int eqn, peqn, var, pvar;

  struct Basis_Functions *bfm;
  double wt, det_J, h3, phi_i, phi_j, flux;
  double d_qrad_dT;
  double source;

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

  flux = htc * (Tref - fv->T) + emiss * sigma * (pow(Tref, 4.0) - pow(fv->T, 4.0));
  d_qrad_dT = -htc - 4. * emiss * sigma * pow(fv->T, 3.0);

  /* finite difference calculation of path dependencies for
     subelement integration
   */
  if (ls->CalcSurfDependencies) {
    peqn = upd->ep[pg->imtrx][eqn];
    bfm = bf[eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][i];

      if (ei[pg->imtrx]->active_interp_ledof[ledof]) {
        ii = ei[pg->imtrx]->lvdof_to_row_lvdof[eqn][i];

        phi_i = bfm->phi[i];

        source = phi_i * flux * lsi->delta;

        source *= det_J * wt;

        source *= h3;

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
    return (0);
  }

  /*
   * Wesiduals
   * ________________________________________________________________________________
   */
  if (af->Assemble_Residual) {
    peqn = upd->ep[pg->imtrx][eqn];
    bfm = bf[eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][i];

      if (ei[pg->imtrx]->active_interp_ledof[ledof]) {
        ii = ei[pg->imtrx]->lvdof_to_row_lvdof[eqn][i];

        phi_i = bfm->phi[i];

        source = phi_i * flux * lsi->delta;

        source *= det_J * wt;

        source *= h3;

        source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

        lec->R[LEC_R_INDEX(peqn, ii)] += source;
      }
    }
  }

  if (af->Assemble_Jacobian) {
    peqn = upd->ep[pg->imtrx][eqn];
    bfm = bf[eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][i];

      if (ei[pg->imtrx]->active_interp_ledof[ledof]) {
        ii = ei[pg->imtrx]->lvdof_to_row_lvdof[eqn][i];

        phi_i = bfm->phi[i];

        /*
         * J_T_T
         */
        var = TEMPERATURE;

        if (pd->v[pg->imtrx][var]) {
          pvar = upd->vp[pg->imtrx][var];

          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];

            source = phi_i * lsi->delta * d_qrad_dT * phi_j;

            source *= det_J * wt * h3;
            source *= pd->etm[pg->imtrx][eqn][LOG2_SOURCE];

            lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
          }
        }
        /*
         * J_e_F
         */
        var = FILL;

        if (pd->v[pg->imtrx][var]) {
          pvar = upd->vp[pg->imtrx][var];

          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            source = phi_i * flux * lsi->d_delta_dF[j];

            source *= det_J * wt * h3;
            source *= pd->etm[pg->imtrx][eqn][LOG2_SOURCE];

            lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
          }
        }
      }
    }
  }
  return (1);
}

int assemble_cont_t_source(double *xi) {
  int i, j, p, ii, ledof;
  int eqn, peqn, var, pvar;

  double phi_diff[MDE];
  double T_diff;
  double grad_T_diff[DIM];

  struct Basis_Functions *bfm;
  double wt, det_J, h3, phi_i, phi_j;
  int xfem_active, extended_dof, base_interp, base_dof;

  double source;

  eqn = R_ENERGY;
  if (!pd->e[pg->imtrx][eqn]) {
    return (0);
  }

  /* For LS matching conditions, just warn and exit if attempted on
     non-discontinuous interpolations */
  if (pd->i[pg->imtrx][eqn] != I_Q1_XV && pd->i[pg->imtrx][eqn] != I_Q2_XV &&
      pd->i[pg->imtrx][eqn] != I_Q1_G && pd->i[pg->imtrx][eqn] != I_Q2_G) {
    GOMA_WH(-1, "Warning: attempting to apply LS_CONT_T without discontinuous enrichment.\n");
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
   * non-xfem dofs do not need a bc and xfem
   * dofs will receive SIC
   */

  /* DRN: to get this to work for _G based enrichment
     will require making the "natural" dofs act like
     the continuous dofs in _XG type enrichment.
     I think that this could be accomplished via
     a collocated bc, but haven't implemented it yet.
   */
  eqn = R_ENERGY;
  if (pd->i[pg->imtrx][eqn] == I_Q1_G || pd->i[pg->imtrx][eqn] == I_Q2_G)
    GOMA_EH(GOMA_ERROR, "LS_CONT_T_BC not yet implemented for _G type enrichment.");

  /* need difference in T between other side and this side */
  var = TEMPERATURE;
  xfem_var_diff(var, &T_diff, phi_diff, grad_T_diff);

  /* finite difference calculation of path dependencies for
     subelement integration
   */
  if (ls->CalcSurfDependencies) {
    /* this is currently handled below using the grad_T_diff */
    return (0);
  }
  /*
   * Wesiduals
   * ________________________________________________________________________________
   */
  if (af->Assemble_Residual) {
    peqn = upd->ep[pg->imtrx][eqn];
    bfm = bf[eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][i];

      if (ei[pg->imtrx]->active_interp_ledof[ledof]) {
        xfem_dof_state(i, pd->i[pg->imtrx][eqn], ei[pg->imtrx]->ielem_shape, &xfem_active,
                       &extended_dof, &base_interp, &base_dof);

        if (extended_dof) {
          ii = ei[pg->imtrx]->lvdof_to_row_lvdof[eqn][i];
          phi_i = bfm->phi[i];

          source = phi_i * lsi->delta * T_diff;

          source *= det_J * wt * h3 * BIG_PENALTY;
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
    eqn = R_ENERGY;
    peqn = upd->ep[pg->imtrx][eqn];
    bfm = bf[eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {

      ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][i];

      if (ei[pg->imtrx]->active_interp_ledof[ledof]) {

        xfem_dof_state(i, pd->i[pg->imtrx][eqn], ei[pg->imtrx]->ielem_shape, &xfem_active,
                       &extended_dof, &base_interp, &base_dof);

        if (extended_dof) {
          ii = ei[pg->imtrx]->lvdof_to_row_lvdof[eqn][i];
          phi_i = bfm->phi[i];

          /*
           * J_T_T
           */
          var = TEMPERATURE;

          if (pd->v[pg->imtrx][var]) {
            pvar = upd->vp[pg->imtrx][var];

            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              phi_j = bf[var]->phi[j];

              source = phi_i * lsi->delta * phi_diff[j];

              source *= det_J * wt * h3 * BIG_PENALTY;
              source *= pd->etm[pg->imtrx][eqn][LOG2_SOURCE];

              lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
            }
          }
          /* path dependence integral */
          if (ls->on_sharp_surf && !ls->AdaptIntegration) {
            /*
             * J_e_F
             */
            var = FILL;

            if (pd->v[pg->imtrx][var]) {
              pvar = upd->vp[pg->imtrx][var];

              for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                phi_j = bf[var]->phi[j];

                source = 0.;
                for (p = 0; p < VIM; p++) {
                  source -= phi_i * phi_j * lsi->normal[p] * grad_T_diff[p];
                }

                for (p = 0; p < VIM; p++) {
                  source -= phi_j * lsi->normal[p] * bfm->grad_phi[i][p] * T_diff;
                }

                source /= lsi->gfmag;
                source *= det_J * wt * h3 * BIG_PENALTY;
                source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

                lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
              }
            }
          } else {
            /* diffuse interface version of path dependence integral */
            /*
             * J_e_F
             */
            var = FILL;
            if (pd->v[pg->imtrx][var]) {
              pvar = upd->vp[pg->imtrx][var];

              for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                source = phi_i * T_diff;
                source *= lsi->d_delta_dF[j];
                source *= det_J * wt * h3 * BIG_PENALTY;
                source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

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

int assemble_ls_yflux_source(int wspec,              /* species number of this boundary condition */
                             double mass_tran_coeff, /* Mass transfer coefficient       */
                             double Y_c,             /* bath concentration 	                     */
                             double dt,              /* current value of the time step            */
                             double tt,
                             double time,
                             int bc_input_id,
                             struct Boundary_Condition *BC_Types) {
  int i, j, ii, ledof;
  int eqn, peqn, var, pvar, b, w;

  struct Basis_Functions *bfm;
  double wt, det_J, h3, phi_i, phi_j, flux;
  double source;

  double Y_w; /* local concentration of current species */
  double sign = 1.;
  double vnorm;
  NORMAL_VELOCITY_DEPENDENCE_STRUCT d_vnorm_struct;
  NORMAL_VELOCITY_DEPENDENCE_STRUCT *d_vnorm = &d_vnorm_struct;
  struct Boundary_Condition *fluxbc;

  eqn = R_MASS;
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

  Y_w = fv->c[wspec];
  /*if ( Y_w > 1. ) Y_w = 1.;*/
  mass_flux_surf_mtc(mp->mass_flux, mp->d_mass_flux, fv->T, fv->c, wspec, mass_tran_coeff, Y_c);

  fluxbc = BC_Types + bc_input_id;
  compute_leak_velocity(&vnorm, d_vnorm, tt, dt, NULL, fluxbc);

  flux = -mp->mass_flux[wspec];

  flux += Y_w * sign * vnorm;

  /*
  if ( fv->c[wspec] > 1. )
    fprintf(stderr,"flux=%g, Y_w=%g, mass_flux=%g, extv=%g,
  vnorm=%g\n",flux,fv->c[wspec],-mp->mass_flux[wspec],fv->ext_v,vnorm);
  */

  /* finite difference calculation of path dependencies for
     subelement integration
   */
  if (ls->CalcSurfDependencies) {
    peqn = MAX_PROB_VAR + wspec;
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
    peqn = MAX_PROB_VAR + wspec;
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
    peqn = MAX_PROB_VAR + wspec;
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

            source = -mp->d_mass_flux[wspec][TEMPERATURE] * phi_j;
            source += Y_w * d_vnorm->T[j] * sign;

            source *= phi_i * lsi->delta;

            source *= det_J * wt * h3;
            /*source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];*/

            lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
          }
        }

        /*
         * J_W_V
         */
        var = VOLTAGE;

        if (pd->v[pg->imtrx][var]) {
          pvar = upd->vp[pg->imtrx][var];

          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];

            source = -mp->d_mass_flux[wspec][VOLTAGE] * phi_j;
            source += Y_w * d_vnorm->V[j] * sign;

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

              source = -mp->d_mass_flux[wspec][MAX_VARIABLE_TYPES + w] * phi_j;
              if (w == wspec)
                source += phi_j * vnorm * sign;
              source += Y_w * d_vnorm->C[w][j] * sign;

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

              source = Y_w * d_vnorm->v[b][j] * sign;

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

            source = Y_w * d_vnorm->F[j] * sign;

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
}

int assemble_cont_vel_source(double *xi, Exo_DB *exo) {
  int a, i, j, ii, ledof;
  int p;
  int eqn, peqn, var, pvar;
  double source;

  double phi_diff[3][MDE];
  double v_diff[3];
  double grad_v_diff[3][3];

  struct Basis_Functions *bfm;
  double wt, det_J, h3, phi_i, phi_j;
  int xfem_active, extended_dof, base_interp, base_dof;

  eqn = R_MOMENTUM1;
  if (!pd->e[pg->imtrx][eqn]) {
    return (0);
  }

  /* For LS matching conditions, just warn and exit if attempted on
     non-discontinuous interpolations */
  if (pd->i[pg->imtrx][eqn] != I_Q1_XV && pd->i[pg->imtrx][eqn] != I_Q2_XV &&
      pd->i[pg->imtrx][eqn] != I_Q1_G && pd->i[pg->imtrx][eqn] != I_Q2_G) {
    GOMA_WH(-1, "Warning: attempting to apply LS_CONT_VEL without discontinuous enrichment.\n");
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
   * non-xfem dofs do not need a bc and xfem
   * dofs will receive SIC
   */

  /* DRN: to get this to work for _G based enrichment
     will require making the "natural" dofs act like
     the continuous dofs in _XG type enrichment.
     I think that this could be accomplished via
     a collocated bc, but haven't implemented it yet.
   */
  eqn = R_MOMENTUM1;
  if (pd->i[pg->imtrx][eqn] == I_Q1_G || pd->i[pg->imtrx][eqn] == I_Q2_G)
    GOMA_EH(GOMA_ERROR, "LS_CONT_VEL_BC not yet implemented for _G type enrichment.");

  for (a = 0; a < WIM; a++) {
    var = VELOCITY1 + a;

    /* need difference in vel between other side and this side */
    /* note difference in convection for indices of grad_v_diff from grad_v */
    xfem_var_diff(var, &(v_diff[a]), &(phi_diff[a][0]), &(grad_v_diff[a][0]));
  }

  /* finite difference calculation of path dependencies for
     subelement integration
   */
  if (ls->CalcSurfDependencies) {
    /* this is currently handled below using the grad_v_diff */
    return (0);
  }

  /*
   * Wesiduals
   * ________________________________________________________________________________
   */
  for (a = 0; a < WIM; a++) {
    eqn = R_MOMENTUM1 + a;

    if (af->Assemble_Residual) {
      peqn = upd->ep[pg->imtrx][eqn];
      bfm = bf[eqn];

      for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
        ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][i];

        if (ei[pg->imtrx]->active_interp_ledof[ledof]) {
          xfem_dof_state(i, pd->i[pg->imtrx][eqn], ei[pg->imtrx]->ielem_shape, &xfem_active,
                         &extended_dof, &base_interp, &base_dof);

          if (extended_dof) {
            ii = ei[pg->imtrx]->lvdof_to_row_lvdof[eqn][i];
            phi_i = bfm->phi[i];

            source = phi_i * lsi->delta * v_diff[a];
            source *= det_J * wt * h3 * BIG_PENALTY;
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
      peqn = upd->ep[pg->imtrx][eqn];
      bfm = bf[eqn];

      for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
        ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][i];

        if (ei[pg->imtrx]->active_interp_ledof[ledof]) {
          xfem_dof_state(i, pd->i[pg->imtrx][eqn], ei[pg->imtrx]->ielem_shape, &xfem_active,
                         &extended_dof, &base_interp, &base_dof);

          if (extended_dof) {
            ii = ei[pg->imtrx]->lvdof_to_row_lvdof[eqn][i];
            phi_i = bfm->phi[i];

            /*
             * J_m_v
             */
            var = eqn;

            if (pd->v[pg->imtrx][var]) {
              pvar = upd->vp[pg->imtrx][var];

              for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                source = phi_i * lsi->delta * phi_diff[a][j];

                source *= det_J * wt * h3 * BIG_PENALTY;
                source *= pd->etm[pg->imtrx][eqn][LOG2_SOURCE];

                lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
              }
            }

            /* path dependence integral */
            if (ls->on_sharp_surf && !ls->AdaptIntegration) {
              /*
               * J_m_F
               */
              var = FILL;

              if (pd->v[pg->imtrx][var]) {
                pvar = upd->vp[pg->imtrx][var];

                for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                  phi_j = bf[var]->phi[j];

                  source = 0.;
                  for (p = 0; p < VIM; p++) {
                    source -= phi_i * phi_j * lsi->normal[p] * grad_v_diff[a][p];
                  }

                  for (p = 0; p < VIM; p++) {
                    source -= phi_j * lsi->normal[p] * bfm->grad_phi[i][p] * v_diff[a];
                  }

                  source *= lsi->gfmaginv;
                  source *= det_J * wt * h3 * BIG_PENALTY;
                  source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

                  lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
                }
              }
            } else {
              /* diffuse interface version of path dependence integral */
              /*
               * J_m_F
               */
              var = FILL;
              if (pd->v[pg->imtrx][var]) {
                pvar = upd->vp[pg->imtrx][var];

                for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                  source = phi_i * v_diff[a];
                  source *= lsi->d_delta_dF[j];
                  source *= det_J * wt * h3 * BIG_PENALTY;
                  source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

                  lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
                }
              }
            }
          }
        }
      }
    }
  }

  return (1);
}

int assemble_extv_kinematic(dbl tt, /* parameter to vary time integration from
                                     * explicit (tt = 1) to implicit (tt = 0)    */
                            dbl dt, /* current value of the time step            */
                            dbl time,
                            int bc_input_id,
                            struct Boundary_Condition *BC_Types) {
  int a, i, j, ii, ledof, w;
  int eqn, peqn, var, pvar;
  int dim;

  struct Basis_Functions *bfm;
  double wt, det_J, h3, phi_i, phi_j;

  double source;
  double vnorm = 0., coeff = 1.;
  double sign = 1.;
  double tau = 0.;
  NORMAL_VELOCITY_DEPENDENCE_STRUCT d_vnorm_struct;
  NORMAL_VELOCITY_DEPENDENCE_STRUCT *d_vnorm = &d_vnorm_struct;
  NORMAL_VELOCITY_DEPENDENCE_STRUCT d_coeff_struct;
  NORMAL_VELOCITY_DEPENDENCE_STRUCT *d_coeff = &d_coeff_struct;

  struct Boundary_Condition *bc;
  struct Boundary_Condition *fluxbc;

  eqn = R_EXT_VELOCITY;
  if (!pd->e[pg->imtrx][eqn]) {
    return (0);
  }

  dim = pd->Num_Dim;

  wt = fv->wt;
  h3 = fv->h3;

  if (ls->on_sharp_surf) /* sharp interface */
  {
    det_J = fv->sdet;
  } else /* diffuse interface */
  {
    det_J = bf[eqn]->detJ;
  }

  /*  initialize derivatives upfront	*/
  memset(d_vnorm->v, 0, sizeof(dbl) * DIM * MDE);
  memset(d_vnorm->T, 0, sizeof(dbl) * MDE);
  memset(d_vnorm->C, 0, sizeof(dbl) * MAX_CONC * MDE);
  memset(d_vnorm->V, 0, sizeof(dbl) * MDE);
  memset(d_vnorm->F, 0, sizeof(dbl) * MDE);
  memset(d_vnorm->X, 0, sizeof(dbl) * MDE * DIM);
  memset(d_coeff->v, 0, sizeof(dbl) * DIM * MDE);
  memset(d_coeff->T, 0, sizeof(dbl) * MDE);
  memset(d_coeff->C, 0, sizeof(dbl) * MAX_CONC * MDE);
  memset(d_coeff->V, 0, sizeof(dbl) * MDE);
  memset(d_coeff->F, 0, sizeof(dbl) * MDE);
  memset(d_coeff->X, 0, sizeof(dbl) * MDE * DIM);

  /* precompute vnorm and its dependencies */
  switch (BC_Types[bc_input_id].BC_Name) {
  case -1:
#if 0
  case LS_EXTV_KINEMATIC_USER:
#endif
  {
    double r = sqrt(fv->x[0] * fv->x[0] + fv->x[1] * fv->x[1]);
    double theta = atan2(fv->x[1], fv->x[0]);
    double r0 = 3.;
    vnorm = ((r - r0) + 1.) * (2. + sin(4. * theta));
  } break;

  case LS_EXTV_KINEMATIC_BC: /* fluid velocity */
  {
    /* sign is correct on this bc */
    sign = 1.;

    vnorm = 0.;
    for (a = 0; a < VIM; a++)
      vnorm += fv->v[a] * lsi->normal[a];

    for (a = 0; a < VIM; a++) {
      var = VELOCITY1 + a;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          phi_j = bf[var]->phi[j];
          d_vnorm->v[a][j] = lsi->normal[a] * phi_j;
        }
      }
    }

    var = FILL;
    if (pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_vnorm->F[j] = 0.;
        for (a = 0; a < WIM; a++) {
          d_vnorm->F[j] += fv->v[a] * lsi->d_normal_dF[a][j];
        }
      }
    }
  } break;

  case LS_EXTV_KIN_LEAK_BC: {
    /* sign needs to be corrected for this bc */
    sign = 1. * ls->Elem_Sign;

    bc = BC_Types + bc_input_id;
    if (bc->BC_Data_Int[1] == -1)
      fluxbc = NULL;
    else
      fluxbc = BC_Types + bc->BC_Data_Int[1];
    compute_leak_velocity(&vnorm, d_vnorm, tt, dt, bc, fluxbc);

    /* DRN:  This isn't too pretty.  Sorry.
       Basically, we need to know if the fluid velocity is to be included
       in the extension velocity or if it will be handled separately in
       the level set advection equation.
       If tran->Fill_Equation == FILL_EQN_EXT_V, we will assume you want
       the extension velocity to INCLUDE the fluid velocity.
       Otherwise, we will EXCLUDE the fluid velocity and assume that
       the level set advection equation will take this into account
    */
    if (pd->v[pg->imtrx][VELOCITY1] && tran->Fill_Equation == FILL_EQN_EXT_V) {
      for (a = 0; a < VIM; a++)
        vnorm += sign * fv->v[a] * lsi->normal[a];

      for (a = 0; a < VIM; a++) {
        var = VELOCITY1 + a;
        if (pd->v[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];
            d_vnorm->v[a][j] = sign * lsi->normal[a] * phi_j;
          }
        }
      }

      var = FILL;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_vnorm->F[j] = 0.;
          for (a = 0; a < WIM; a++) {
            d_vnorm->F[j] += sign * fv->v[a] * lsi->d_normal_dF[a][j];
          }
        }
      }
    }
  } break;

  case LS_EXTV_LATENT_BC: {
    /* sign needs to be corrected for this bc */

    bc = BC_Types + bc_input_id;
    if (bc->BC_Data_Int[1] == -1)
      fluxbc = NULL;
    else
      fluxbc = BC_Types + bc->BC_Data_Int[1];

    tau = bc->BC_Data_Float[0];
    vnorm = bc->BC_Data_Float[1] * (fv->T - fluxbc->BC_Data_Float[0]);

    /* determine sign based on normal temperature gradient	*/
    coeff = 0.0;
    for (a = 0; a < dim; a++)
      coeff += fv->grad_T[a] * lsi->normal[a];
    coeff *= tran->delta_t_avg;

    sign = -1.;
    var = TEMPERATURE;
    if (pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_vnorm->T[j] += bc->BC_Data_Float[1] * bf[var]->phi[j];
        for (a = 0; a < dim; a++)
          d_coeff->T[j] += bf[var]->grad_phi[j][a] * lsi->normal[a];
        d_coeff->T[j] *= tran->delta_t_avg;
      }
    }

    var = ls->var;
    if (pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        for (a = 0; a < dim; a++) {
          d_coeff->F[j] += fv->grad_T[a] * lsi->d_normal_dF[a][j];
        }
        d_coeff->F[j] *= tran->delta_t_avg;
      }
    }

    /* DRN:  This isn't too pretty.  Sorry.
       Basically, we need to know if the fluid velocity is to be included
       in the extension velocity or if it will be handled separately in
       the level set advection equation.
       If tran->Fill_Equation == FILL_EQN_EXT_V, we will assume you want
       the extension velocity to INCLUDE the fluid velocity.
       Otherwise, we will EXCLUDE the fluid velocity and assume that
       the level set advection equation will take this into account
    */
    if (pd->v[pg->imtrx][VELOCITY1] && tran->Fill_Equation == FILL_EQN_EXT_V) {
      for (a = 0; a < VIM; a++)
        vnorm += fv->v[a] * lsi->normal[a];

      for (a = 0; a < VIM; a++) {
        var = VELOCITY1 + a;
        if (pd->v[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];
            d_vnorm->v[a][j] = lsi->normal[a] * phi_j;
          }
        }
      }

      var = ls->var;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_vnorm->F[j] = 0.;
          for (a = 0; a < WIM; a++) {
            d_vnorm->F[j] += fv->v[a] * lsi->d_normal_dF[a][j];
          }
        }
      }
    }
  } break;

  default:
    sprintf(Err_Msg, "BC %s not found", BC_Types[bc_input_id].desc->name1);
    GOMA_EH(GOMA_ERROR, Err_Msg);
    break;
  }
#if 0
  DPRINTF(stderr,"kinematic extv at (%g,%g) vnorm=%g, fdot=%g\n",fv->x[0],fv->x[1],vnorm,fv_dot->F);
#endif

#define USE_GFMAG 1
  /* finite difference calculation of path dependencies for
     subelement integration
   */
  if (ls->CalcSurfDependencies) {
    peqn = upd->ep[pg->imtrx][eqn];
    bfm = bf[eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][i];

      if (ei[pg->imtrx]->active_interp_ledof[ledof]) {
        ii = ei[pg->imtrx]->lvdof_to_row_lvdof[eqn][i];

        phi_i = bfm->phi[i];

#if USE_GFMAG
        source = 2. * lsi->delta * (-coeff * (tau * fv_dot->ext_v + fv->ext_v) + vnorm * sign) *
                 lsi->gfmag;
#else
        source = 2. * lsi->delta * (-coeff * (tau * fv_dot->ext_v + fv->ext_v) + vnorm * sign);
#endif

        source *= phi_i * det_J * wt * h3;
        /* no such term multiplier yet
          source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
        */

        /* J_m_F
         */
        var = ls->var;
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
    peqn = upd->ep[pg->imtrx][eqn];
    bfm = bf[eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][i];

      if (ei[pg->imtrx]->active_interp_ledof[ledof]) {
        ii = ei[pg->imtrx]->lvdof_to_row_lvdof[eqn][i];

        phi_i = bfm->phi[i];

#if USE_GFMAG
        source = 2. * lsi->delta * (-coeff * (tau * fv_dot->ext_v + fv->ext_v) + vnorm * sign) *
                 lsi->gfmag;
#else
        source = 2. * lsi->delta * (-coeff * (tau * fv_dot->ext_v + fv->ext_v) + vnorm * sign);
#endif

        source *= phi_i * det_J * wt * h3;
        /* no such term multiplier yet
          source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
        */
        lec->R[LEC_R_INDEX(peqn, ii)] += source;
      }
    }
  }

  /*
   * Yacobian terms...
   */

  if (af->Assemble_Jacobian) {
    peqn = upd->ep[pg->imtrx][eqn];
    bfm = bf[eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {

      ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][i];

      if (ei[pg->imtrx]->active_interp_ledof[ledof]) {

        ii = ei[pg->imtrx]->lvdof_to_row_lvdof[eqn][i];

        phi_i = bfm->phi[i];

        /*
         * J_ext_v_ext_v
         */
        /*              var = eqn;  */
        var = EXT_VELOCITY;
        if (pd->v[pg->imtrx][var]) {
          pvar = upd->vp[pg->imtrx][var];

          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];

#if USE_GFMAG
            source = 2. * lsi->delta * lsi->gfmag *
                     (-coeff * (tau * (1. + 2. * tt) * phi_j / dt + phi_j));
#else
            source = 2. * lsi->delta * (-coeff * (tau * (1. + 2. * tt) * phi_j / dt + phi_j));
#endif

            source *= phi_i * det_J * wt * h3;
            /* no such term multiplier yet
              source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
            */

            lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
          }
        }

        /*
         * J_ext_v_v
         */
        for (a = 0; a < VIM; a++) {
          var = VELOCITY1 + a;

          if (pd->v[pg->imtrx][var]) {
            pvar = upd->vp[pg->imtrx][var];

            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              phi_j = bf[var]->phi[j];

#if USE_GFMAG
              source =
                  2. * lsi->delta * lsi->gfmag *
                  (-d_coeff->v[a][j] * (tau * fv_dot->ext_v + fv->ext_v) + sign * d_vnorm->v[a][j]);
#else
              source =
                  2. * lsi->delta *
                  (-d_coeff->v[a][j] * (tau * fv_dot->ext_v + fv->ext_v) + sign * d_vnorm->v[a][j]);
#endif

              source *= phi_i * det_J * wt * h3;
              /* no such term multiplier yet
                source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
              */

              lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
            }
          }
        }

        /*
         * J_ext_v_F
         */
        var = ls->var;
        if (pd->v[pg->imtrx][var]) {
          pvar = upd->vp[pg->imtrx][var];

          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];

#if USE_GFMAG
            source =
                2. * lsi->delta * lsi->gfmag *
                    (-d_coeff->F[j] * (tau * fv_dot->ext_v + fv->ext_v) + d_vnorm->F[j] * sign) +
                2. * (-coeff * (tau * fv_dot->ext_v + fv->ext_v) + vnorm * sign) *
                    (lsi->delta * lsi->d_gfmag_dF[j] + lsi->d_delta_dF[j] * lsi->gfmag);
#else
            source =
                2. * lsi->delta *
                    (-d_coeff->F[j] * (tau * fv_dot->ext_v + fv->ext_v) + d_vnorm->F[j] * sign) +
                2. * (-coeff * (tau * fv_dot->ext_v + fv->ext_v) + vnorm * sign) *
                    lsi->d_delta_dF[j];
#endif

            source *= phi_i * det_J * wt * h3;
            /* no such term multiplier yet
              source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
            */

            lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
          }
        }

        /*
         * J_ext_v_T
         */
        var = TEMPERATURE;
        if (pd->v[pg->imtrx][var]) {
          pvar = upd->vp[pg->imtrx][var];

          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];

#if USE_GFMAG
            source = 2. * lsi->delta * lsi->gfmag *
                     (-d_coeff->T[j] * (tau * fv_dot->ext_v + fv->ext_v) + d_vnorm->T[j] * sign);
#else
            source = 2. * lsi->delta *
                     (-d_coeff->T[j] * (tau * fv_dot->ext_v + fv->ext_v) + d_vnorm->T[j] * sign);
#endif
            source *= phi_i * det_J * wt * h3;
            /* no such term multiplier yet
              source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
            */

            lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
          }
        }

        /*
         * J_ext_v_V
         */
        var = VOLTAGE;
        if (pd->v[pg->imtrx][var]) {
          pvar = upd->vp[pg->imtrx][var];

          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];

#if USE_GFMAG
            source = 2. * lsi->delta * lsi->gfmag *
                     (-d_coeff->V[j] * (tau * fv_dot->ext_v + fv->ext_v) + d_vnorm->V[j] * sign);
#else
            source = 2. * lsi->delta *
                     (-d_coeff->V[j] * (tau * fv_dot->ext_v + fv->ext_v) + d_vnorm->V[j] * sign);
#endif
            source *= phi_i * det_J * wt * h3;
            /* no such term multiplier yet
              source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
            */

            lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
          }
        }

        /*
         * J_ext_v_C
         */
        var = MASS_FRACTION;
        if (pd->v[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            for (w = 0; w < pd->Num_Species_Eqn; w++) {
              phi_j = bf[var]->phi[j];
#if USE_GFMAG
              source =
                  2. * lsi->delta * lsi->gfmag *
                  (-d_coeff->C[w][j] * (tau * fv_dot->ext_v + fv->ext_v) + d_vnorm->C[w][j] * sign);
#else
              source =
                  2. * lsi->delta *
                  (-d_coeff->C[w][j] * (tau * fv_dot->ext_v + fv->ext_v) + d_vnorm->C[w][j] * sign);
#endif
              source *= phi_i * det_J * wt * h3;
              /* no such term multiplier yet
                 source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
               */

              lec->J[LEC_J_INDEX(peqn, MAX_PROB_VAR + w, ii, j)] += source;
            }
          }
        }
      }
    }
  }

  return (1);
}

int assemble_interface_extension_velocity_sic(int ext_vel_sign) {
  int a, i, j, ii, ledof;
  int eqn, peqn, var, pvar;

  struct Basis_Functions *bfm;
  double wt, det_J, h3, phi_i, phi_j;

  double source;
  double F_i;

  eqn = R_EXT_VELOCITY;
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
   * Wesiduals
   * ________________________________________________________________________________
   */
  if (af->Assemble_Residual) {
    peqn = upd->ep[pg->imtrx][eqn];
    bfm = bf[eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][i];

      F_i = dof_distance(eqn, i);

      if (ei[pg->imtrx]->active_interp_ledof[ledof] && !sign_change(F_i, (double)ext_vel_sign)) {
        ii = ei[pg->imtrx]->lvdof_to_row_lvdof[eqn][i];

        phi_i = bfm->phi[i];

        source = phi_i * lsi->delta * fv->ext_v;

        for (a = 0; a < VIM; a++) {
          source -= phi_i * lsi->delta * fv->v[a] * lsi->normal[a];
        }

        source *= det_J * wt * h3 * BIG_PENALTY;
        /* no such term multiplier yet
          source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
        */

        lec->R[LEC_R_INDEX(peqn, ii)] += source;
      }
    }
  }

  /*
   * Yacobian terms...
   */

  if (af->Assemble_Jacobian) {
    peqn = upd->ep[pg->imtrx][eqn];
    bfm = bf[eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {

      ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][i];

      F_i = dof_distance(eqn, i);
      if (ei[pg->imtrx]->active_interp_ledof[ledof] && !sign_change(F_i, (double)ext_vel_sign)) {

        ii = ei[pg->imtrx]->lvdof_to_row_lvdof[eqn][i];

        phi_i = bfm->phi[i];

        /*
         * J_ext_v_ext_v
         */
        var = eqn;

        if (pd->v[pg->imtrx][var]) {
          pvar = upd->vp[pg->imtrx][var];

          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];

            source = phi_i * lsi->delta * phi_j;

            source *= det_J * wt * h3 * BIG_PENALTY;
            /* no such term multiplier yet
              source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
            */

            lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
          }
        }

        /*
         * J_ext_v_v
         */
        for (a = 0; a < VIM; a++) {
          var = VELOCITY1 + a;

          if (pd->v[pg->imtrx][var]) {
            pvar = upd->vp[pg->imtrx][var];

            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              phi_j = bf[var]->phi[j];

              source = -phi_i * lsi->delta * phi_j * lsi->normal[a];

              source *= det_J * wt * h3 * BIG_PENALTY;
              /* no such term multiplier yet
                source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
              */

              lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
            }
          }
        }

        /*
         * J_ext_v_F
         */
        var = FILL;

        if (pd->v[pg->imtrx][var]) {
          pvar = upd->vp[pg->imtrx][var];

          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];

            source = phi_i * lsi->d_delta_dF[j] * fv->ext_v;

            for (a = 0; a < VIM; a++) {
              source -= phi_i * lsi->d_delta_dF[j] * fv->v[a] * lsi->normal[a];
            }

            for (a = 0; a < VIM; a++) {
              source -= phi_i * lsi->delta * fv->v[a] * lsi->d_normal_dF[a][j];
            }

            source *= det_J * wt * h3 * BIG_PENALTY;
            /* no such term multiplier yet
              source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
            */

            lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
          }
        }
      }
    }
  }

  return (1);
}

int assemble_eik_kinematic(dbl tt, /* parameter to vary time integration from
                                    * explicit (tt = 1) to implicit (tt = 0)    */
                           dbl dt, /* current value of the time step            */
                           dbl time,
                           int bc_input_id,
                           struct Boundary_Condition *BC_Types) {
  int a, i, j, ii, ledof, w;
  int eqn, peqn, var, pvar;

  struct Basis_Functions *bfm;
  double wt, det_J, h3, phi_i, phi_j;
  double vnorm;
  double sign = 1.;
  NORMAL_VELOCITY_DEPENDENCE_STRUCT d_vnorm_struct;
  NORMAL_VELOCITY_DEPENDENCE_STRUCT *d_vnorm = &d_vnorm_struct;
  double mass = 0.;
  double advection = 0.;
  double source = 0.;
  double penalty = 1.;

  struct Boundary_Condition *bc;
  struct Boundary_Condition *fluxbc;

  eqn = R_FILL;
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

  /* precompute vnorm and its dependencies */
  switch (BC_Types[bc_input_id].BC_Name) {
  case -1:
#if 0
  case LS_EIK_KINEMATIC_USER:
#endif
  {
    double r = sqrt(fv->x[0] * fv->x[0] + fv->x[1] * fv->x[1]);
    double theta = atan2(fv->x[1], fv->x[0]);
    double r0 = 3.;
    vnorm = ((r - r0) + 1.) * (2. + sin(4. * theta));

    memset(d_vnorm->v, 0, sizeof(dbl) * DIM * MDE);
    memset(d_vnorm->T, 0, sizeof(dbl) * MDE);
    memset(d_vnorm->C, 0, sizeof(dbl) * MAX_CONC * MDE);
    memset(d_vnorm->V, 0, sizeof(dbl) * MDE);
    memset(d_vnorm->F, 0, sizeof(dbl) * MDE);
  } break;

  case LS_EIK_KINEMATIC_BC: /* fluid velocity */
  {
    /* sign is correct on this bc */
    sign = 1.;

    vnorm = 0.;
#if 0
      for ( a=0; a < VIM; a++ ) vnorm += fv->v[a] * lsi->normal[a];
#else
    for (a = 0; a < VIM; a++)
      vnorm += fv_old->v[a] * lsi->normal[a];
#endif

    memset(d_vnorm->T, 0, sizeof(dbl) * MDE);
    memset(d_vnorm->C, 0, sizeof(dbl) * MAX_CONC * MDE);
    memset(d_vnorm->V, 0, sizeof(dbl) * MDE);

#if 0
      for ( a=0; a < WIM; a++ )
        {
          var = VELOCITY1 + a;
          if ( pd->v[pg->imtrx][var] )
            {
              for( j=0; j<ei[pg->imtrx]->dof[var]; j++)
                {
	          phi_j = bf[var]->phi[j];
		  d_vnorm->v[a][j] = lsi->normal[a] * phi_j;
	        }
	    }
        }
      var = FILL;
      if ( pd->v[pg->imtrx][var] )
        {
          for( j=0; j<ei[pg->imtrx]->dof[var]; j++)
            {
	      d_vnorm->F[j] = 0.;
              for ( a=0; a < WIM; a++ ) 
                {
                  d_vnorm->F[j] += fv->v[a] * lsi->d_normal_dF[a][j];
                }
            }
	}
#else
    memset(d_vnorm->v, 0, sizeof(dbl) * DIM * MDE);

    var = FILL;
    if (pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_vnorm->F[j] = 0.;
        for (a = 0; a < WIM; a++) {
          d_vnorm->F[j] += fv_old->v[a] * lsi->d_normal_dF[a][j];
        }
      }
    }
#endif
  } break;

  case LS_EIK_KIN_LEAK_BC: {
    /* sign needs to be corrected for this bc */
    sign = 1. * ls->Elem_Sign;

    bc = BC_Types + bc_input_id;
    if (bc->BC_Data_Int[1] == -1)
      fluxbc = NULL;
    else
      fluxbc = BC_Types + bc->BC_Data_Int[1];
    compute_leak_velocity(&vnorm, d_vnorm, tt, dt, bc, fluxbc);

    /* DRN:  This isn't too pretty.  Sorry.
       Basically, we need to know if the fluid velocity is to be included
       in the extension velocity or if it will be handled separately in
       the level set advection equation.
       If tran->Fill_Equation == FILL_EQN_EXT_V, we will assume you want
       the extension velocity to INCLUDE the fluid velocity.
       Otherwise, we will EXCLUDE the fluid velocity and assume that
       the level set advection equation will take this into account
    */
    if (pd->v[pg->imtrx][VELOCITY1] && tran->Fill_Equation == FILL_EQN_EXT_V) {
      for (a = 0; a < VIM; a++)
        vnorm += sign * fv->v[a] * lsi->normal[a];

      for (a = 0; a < VIM; a++) {
        var = VELOCITY1 + a;
        if (pd->v[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];
            d_vnorm->v[a][j] = sign * lsi->normal[a] * phi_j;
          }
        }
      }

      var = FILL;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_vnorm->F[j] = 0.;
          for (a = 0; a < WIM; a++) {
            d_vnorm->F[j] += sign * fv->v[a] * lsi->d_normal_dF[a][j];
          }
        }
      }
    }
  } break;

  default:
    sprintf(Err_Msg, "BC %s not found", BC_Types[bc_input_id].desc->name1);
    GOMA_EH(GOMA_ERROR, Err_Msg);
    break;
  }
#if 0
  DPRINTF(stderr,"kinematic extv at (%g,%g) vnorm=%g, fdot=%g\n",fv->x[0],fv->x[1],vnorm,fv_dot->F);
#endif

  if (tt != 0.)
    GOMA_EH(GOMA_ERROR, "LS_EIK_KINEMATIC currently requires backward Euler");

  /* finite difference calculation of path dependencies for
     subelement integration
   */
  if (ls->CalcSurfDependencies) {
    peqn = upd->ep[pg->imtrx][eqn];
    bfm = bf[eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {

      ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][i];

      if (ei[pg->imtrx]->active_interp_ledof[ledof]) {
        ii = ei[pg->imtrx]->lvdof_to_row_lvdof[eqn][i];

        phi_i = bfm->phi[i];

#if 0
	       mass = lsi->delta * 2. * (-fv_old->F);
	       mass *= phi_i * det_J * wt * h3 * penalty;
	       mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
#endif
#if 1
        mass = lsi->delta * 2. * (fv->F - fv_old->F);
        mass *= phi_i * det_J * wt * h3 * penalty;
        mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
#endif
#if 0
	       advection = lsi->delta * 2. * (dt * sign*vnorm);
	       advection *= phi_i * det_J * wt * h3 * penalty;
	       advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
#endif
#if 1
        advection = lsi->delta * 2. * (dt * sign * vnorm * lsi->gfmag);
        advection *= phi_i * det_J * wt * h3 * penalty;
        advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
#endif

        /* J_m_F
         */
        var = LS;
        pvar = upd->vp[pg->imtrx][var];

        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          phi_j = bf[var]->phi[j];
          lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += (mass + advection + source) * phi_j;
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
    peqn = upd->ep[pg->imtrx][eqn];
    bfm = bf[eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {

      ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][i];

      if (ei[pg->imtrx]->active_interp_ledof[ledof]) {
        ii = ei[pg->imtrx]->lvdof_to_row_lvdof[eqn][i];

        phi_i = bfm->phi[i];

#if 0
	       mass = lsi->delta * 2. * (-fv_old->F);
	       mass *= phi_i * det_J * wt * h3 * penalty;
	       mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
#endif
#if 1
        mass = lsi->delta * 2. * (fv->F - fv_old->F);
        mass *= phi_i * det_J * wt * h3 * penalty;
        mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
#endif
#if 0
	       advection = lsi->delta * 2. * (dt * sign*vnorm);
	       advection *= phi_i * det_J * wt * h3 * penalty;
	       advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
#endif
#if 1
        advection = lsi->delta * 2. * (dt * sign * vnorm * lsi->gfmag);
        advection *= phi_i * det_J * wt * h3 * penalty;
        advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
#endif

        lec->R[LEC_R_INDEX(peqn, ii)] += mass + advection + source;
      }
    }

#if 0
{
  fprintf(stderr,"checkkinematic %g %g %g %g %g %g\n",fv->x[0],fv->x[1],fv_dot->F + vnorm * lsi->gfmag,
          -(fv->F - fv_old->F)/dt, vnorm, lsi->gfmag);
}
#endif
#if 0
{
  double theta = atan2( fv->x[1], fv->x[0] );
  double r = sqrt(fv->x[0] * fv->x[0] + fv->x[1] * fv->x[1]);
  double rsoln = 2.+exp(time*(2. + sin(4.*theta)));
  rsoln = 2.-exp(time*(2. - sin(4.*theta))) + exp(time) + exp(3.*time);
  fprintf(stderr,"checkkinematic %g %g %g %g %g %g\n",theta,r,fv_dot->F + vnorm * lsi->gfmag, 
          -(fv->F - fv_old->F)/dt, vnorm, lsi->gfmag);
}
#endif
  }

  /*
   * Yacobian terms...
   */

  if (af->Assemble_Jacobian) {
    peqn = upd->ep[pg->imtrx][eqn];
    bfm = bf[eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {

      ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][i];

      if (ei[pg->imtrx]->active_interp_ledof[ledof]) {

        ii = ei[pg->imtrx]->lvdof_to_row_lvdof[eqn][i];

        phi_i = bfm->phi[i];

        /*
         * J_ls_v
         */
        for (a = 0; a < VIM; a++) {
          var = VELOCITY1 + a;

          if (pd->v[pg->imtrx][var]) {
            pvar = upd->vp[pg->imtrx][var];

            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              phi_j = bf[var]->phi[j];

#if 0
                          advection = lsi->delta * 2. * dt * sign*d_vnorm->v[a][j];
                          advection *= phi_i * det_J * wt * h3 * penalty;
                          advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
#endif
#if 1
              advection = lsi->delta * 2. * dt * lsi->gfmag * sign * d_vnorm->v[a][j];
              advection *= phi_i * det_J * wt * h3 * penalty;
              advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
#endif

              lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += advection;
            }
          }
        }

        /*
         * J_ls_T
         */
        var = TEMPERATURE;
        if (pd->v[pg->imtrx][var]) {
          pvar = upd->vp[pg->imtrx][var];

          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];

#if 0
                      advection = lsi->delta * 2. * dt * sign*d_vnorm->T[j];
                      advection *= phi_i * det_J * wt * h3 * penalty;
                      advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
#endif
#if 1
            advection = lsi->delta * 2. * dt * lsi->gfmag * sign * d_vnorm->T[j];
            advection *= phi_i * det_J * wt * h3 * penalty;
            advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
#endif

            lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += advection;
          }
        }

        /*
         * J_ls_V
         */
        var = VOLTAGE;
        if (pd->v[pg->imtrx][var]) {
          pvar = upd->vp[pg->imtrx][var];

          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];

#if 0
                      advection = lsi->delta * 2. * dt * sign*d_vnorm->V[j];
                      advection *= phi_i * det_J * wt * h3 * penalty;
                      advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
#endif
#if 1
            advection = lsi->delta * 2. * dt * lsi->gfmag * sign * d_vnorm->V[j];
            advection *= phi_i * det_J * wt * h3 * penalty;
            advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
#endif

            lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += advection;
          }
        }

        /*
         * J_ls_C
         */
        var = MASS_FRACTION;
        if (pd->v[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            for (w = 0; w < pd->Num_Species_Eqn; w++) {
              phi_j = bf[var]->phi[j];

#if 0
                          advection = lsi->delta * 2. * dt * sign*d_vnorm->C[w][j];
                          advection *= phi_i * det_J * wt * h3 * penalty;
                          advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
#endif
#if 1
              advection = lsi->delta * 2. * dt * lsi->gfmag * sign * d_vnorm->C[w][j];
              advection *= phi_i * det_J * wt * h3 * penalty;
              advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
#endif

              lec->J[LEC_J_INDEX(peqn, MAX_PROB_VAR + w, ii, j)] += advection;
            }
          }
        }

        /*
         * J_ls_F
         */
        var = FILL;

        if (pd->v[pg->imtrx][var]) {
          pvar = upd->vp[pg->imtrx][var];

          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];

#if 0
		      mass = (lsi->d_delta_dF[j] * 2.) * (-fv_old->F);
		      mass *= phi_i * det_J * wt * h3 * penalty;
		      mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
#endif
#if 1
            mass = lsi->d_delta_dF[j] * 2. * (fv->F - fv_old->F) + lsi->delta * 2. * phi_j;
            mass *= phi_i * det_J * wt * h3 * penalty;
            mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
#endif
#if 0
		      advection = 2. * dt *
		                  (lsi->d_delta_dF[j] * vnorm + 
		                   lsi->delta * sign*d_vnorm->F[j]);
		      advection *= phi_i * det_J * wt * h3 * penalty;
		      advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
#endif
#if 1
            advection =
                2. * dt *
                (vnorm * lsi->gfmag * lsi->d_delta_dF[j] + lsi->delta * vnorm * lsi->d_gfmag_dF[j] +
                 lsi->gfmag * lsi->delta * sign * d_vnorm->F[j]);
            advection *= phi_i * det_J * wt * h3 * penalty;
            advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
#endif

            lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += mass + advection + source;
          }
        }
      }
    }
  }

  return (1);
}

int assemble_p_source(double pressure, const int bcflag) {
  int i, j, a, b, p, q, ii, ledof, w;
  int eqn, peqn, var, pvar;
  int dim;
  int err, elem_sign;

  struct Basis_Functions *bfm;
  double wt, det_J, h3, phi_i, phi_j;

  double source;

  double TT[MAX_PDIM][MAX_PDIM]; /**  solid stresses  **/
  double dTT_drs[MAX_PDIM][MAX_PDIM][DIM][MDE];
  double dTT_dx[MAX_PDIM][MAX_PDIM][MAX_PDIM][MDE];
  double dTT_dp[MAX_PDIM][MAX_PDIM][MDE];
  double dTT_dc[MAX_PDIM][MAX_PDIM][MAX_CONC][MDE];
  double dTT_dp_liq[DIM][DIM][MDE];
  double dTT_dp_gas[DIM][DIM][MDE];
  double dTT_dporosity[DIM][DIM][MDE];
  double dTT_dsink_mass[DIM][DIM][MDE];
  double dTT_dT[MAX_PDIM][MAX_PDIM][MDE];
  double dTT_dmax_strain[DIM][DIM][MDE];
  double dTT_dcur_strain[DIM][DIM][MDE];
  double elast_modulus;

  double Pi[DIM][DIM]; /** liquid stresses  **/
  STRESS_DEPENDENCE_STRUCT d_Pi_struct;
  STRESS_DEPENDENCE_STRUCT *d_Pi = &d_Pi_struct;

  double force[MAX_PDIM], d_force[MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE];

  eqn = R_MOMENTUM1;
  if (!pd->e[pg->imtrx][eqn]) {
    return (0);
  }

  wt = fv->wt;
  h3 = fv->h3;

  if (ls->on_sharp_surf) /* sharp interface */
  {
    det_J = fv->sdet;
    elem_sign = ls->Elem_Sign;
  } else /* diffuse interface */
  {
    det_J = bf[eqn]->detJ;
    elem_sign = 1.;
  }

  dim = pd->Num_Dim;

  memset(force, 0, MAX_PDIM * sizeof(double));
  memset(d_force, 0, MAX_PDIM * (MAX_VARIABLE_TYPES + MAX_CONC) * MDE * sizeof(double));

  if (bcflag == -1) /**   add constant, isotropic stress  **/
  {
    for (a = 0; a < WIM; a++) {
      force[a] = -pressure * lsi->normal[a];
    }
    if (af->Assemble_Jacobian) {
      var = ls->var;
      if (pd->v[pg->imtrx][var]) {
        pvar = upd->vp[pg->imtrx][var];
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          for (a = 0; a < WIM; a++) {
            d_force[a][var][j] -= pressure * lsi->d_normal_dF[a][j];
          }
        }
      }
    }

  } else if (bcflag == 0) /**  subtract off viscous stress  **/
  {
    if (!af->Assemble_Jacobian)
      d_Pi = NULL;
    /* compute stress tensor and its derivatives */
    fluid_stress(Pi, d_Pi);

    for (a = 0; a < WIM; a++) {
      for (i = 0; i < WIM; i++) {
        force[a] = -pressure * Pi[a][i] * lsi->normal[i];
      }
    }

    if (af->Assemble_Jacobian) {
      var = TEMPERATURE;
      if (pd->v[pg->imtrx][var]) {
        for (p = 0; p < pd->Num_Dim; p++) {
          for (q = 0; q < pd->Num_Dim; q++) {
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              d_force[p][var][j] -= pressure * lsi->normal[q] * d_Pi->T[p][q][j];
            }
          }
        }
      }
      var = BOND_EVOLUTION;
      if (pd->v[pg->imtrx][var]) {

        for (p = 0; p < pd->Num_Dim; p++) {
          for (q = 0; q < pd->Num_Dim; q++) {
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              d_force[p][var][j] -= pressure * lsi->normal[q] * d_Pi->nn[p][q][j];
            }
          }
        }
      }

      var = RESTIME;
      if (pd->v[pg->imtrx][var]) {

        for (p = 0; p < pd->Num_Dim; p++) {
          for (q = 0; q < pd->Num_Dim; q++) {
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              d_force[p][var][j] -= pressure * lsi->normal[q] * d_Pi->degrade[p][q][j];
            }
          }
        }
      }

      var = ls->var;
      if (pd->v[pg->imtrx][var]) {
        for (p = 0; p < pd->Num_Dim; p++) {
          for (q = 0; q < pd->Num_Dim; q++) {
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              d_force[p][var][j] -= pressure * (lsi->normal[q] * d_Pi->F[p][q][j] +
                                                lsi->d_normal_dF[p][j] * Pi[p][q]);
            }
          }
        }
      }
      var = MASS_FRACTION;
      if (pd->v[pg->imtrx][var]) {
        for (p = 0; p < pd->Num_Dim; p++) {
          for (q = 0; q < pd->Num_Dim; q++) {
            for (w = 0; w < pd->Num_Species_Eqn; w++) {
              for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                d_force[p][MAX_VARIABLE_TYPES + w][j] -=
                    pressure * lsi->normal[q] * d_Pi->C[p][q][w][j];
              }
            }
          }
        }
      }
      var = PRESSURE;
      if (pd->v[pg->imtrx][var]) {
        for (p = 0; p < pd->Num_Dim; p++) {
          for (q = 0; q < pd->Num_Dim; q++) {
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              d_force[p][var][j] -= pressure * lsi->normal[q] * d_Pi->P[p][q][j];
            }
          }
        }
      }
      if (pd->v[pg->imtrx][VELOCITY1]) {
        for (p = 0; p < pd->Num_Dim; p++) {
          for (q = 0; q < pd->Num_Dim; q++) {
            for (b = 0; b < VIM; b++) {
              var = VELOCITY1 + b;
              for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                d_force[p][var][j] -= pressure * lsi->normal[q] * d_Pi->v[p][q][b][j];
              }
            }
          }
        }
      }
      if (pd->v[pg->imtrx][VORT_DIR1]) {
        for (p = 0; p < pd->Num_Dim; p++) {
          for (q = 0; q < pd->Num_Dim; q++) {
            for (b = 0; b < VIM; b++) {
              var = VORT_DIR1 + b;
              for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                d_force[p][var][j] -= pressure * lsi->normal[q] * d_Pi->vd[p][q][b][j];
              }
            }
          }
        }
      }
      if (pd->v[pg->imtrx][MESH_DISPLACEMENT1]) {
        for (p = 0; p < pd->Num_Dim; p++) {
          for (q = 0; q < pd->Num_Dim; q++) {
            for (b = 0; b < VIM; b++) {
              var = MESH_DISPLACEMENT1 + b;
              for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                d_force[p][var][j] -= pressure * lsi->normal[q] * d_Pi->X[p][q][b][j];
              }
            }
          }
        }
      }
    } /*  if Jacobian  */
  } else if (bcflag == 1 && pd->e[pg->imtrx][R_MESH1] && cr->MeshMotion != ARBITRARY) {
    /* initialize some arrays */
    memset(TT, 0, sizeof(double) * DIM * DIM);
    if (af->Assemble_Jacobian) {
      memset(dTT_dx, 0, sizeof(double) * DIM * DIM * DIM * MDE);
      memset(dTT_dp, 0, sizeof(double) * DIM * DIM * MDE);
      memset(dTT_drs, 0, sizeof(double) * DIM * DIM * MDE);
      memset(dTT_dc, 0, sizeof(double) * DIM * DIM * MAX_CONC * MDE);
      memset(dTT_dp_liq, 0, sizeof(double) * DIM * DIM * MDE);
      memset(dTT_dp_gas, 0, sizeof(double) * DIM * DIM * MDE);
      memset(dTT_dporosity, 0, sizeof(double) * DIM * DIM * MDE);
      memset(dTT_dsink_mass, 0, sizeof(double) * DIM * DIM * MDE);
      memset(dTT_dT, 0, sizeof(double) * DIM * DIM * MDE);
      memset(dTT_dmax_strain, 0, sizeof(double) * DIM * DIM * MDE);
      memset(dTT_dcur_strain, 0, sizeof(double) * DIM * DIM * MDE);
    }

    err = belly_flop(elc->lame_mu);
    GOMA_EH(err, "error in belly flop");
    if (err == 2)
      exit(-1);
    /*
     * Total mesh stress tensor...
     */
    err = mesh_stress_tensor(TT, dTT_dx, dTT_dp, dTT_dc, dTT_dp_liq, dTT_dp_gas, dTT_dporosity,
                             dTT_dsink_mass, dTT_dmax_strain, dTT_dcur_strain, dTT_dT, elc->lame_mu,
                             elc->lame_lambda, tran->delta_t, ei[pg->imtrx]->ielem, 0, 1);
    /* For LINEAR ELASTICITY */
    if (cr->MeshFluxModel == LINEAR) {
      if (dim == 2) {
        TT[2][2] = 1.;
        TT[1][2] = 0.;
        TT[0][2] = 0.;
      }
    }

    /*  For Hookian Elasticity and shrinkage */
    else {
      if (dim == 2) {
        elast_modulus = elc->lame_mu;
        if (cr->MeshMotion == ARBITRARY) {
          TT[2][2] = (1. - fv->volume_change) * elast_modulus;
        } else {
          if (cr->MeshFluxModel == NONLINEAR || cr->MeshFluxModel == HOOKEAN_PSTRAIN ||
              cr->MeshFluxModel == INCOMP_PSTRAIN)
            TT[2][2] = (1. - pow(fv->volume_change, 2. / 3.)) * elast_modulus - fv->P;
          /*              if (cr->MeshFluxModel == INCOMP_PSTRESS) */
          else
            TT[2][2] = 0.;
        }
        TT[1][2] = 0.;
        TT[0][2] = 0.;
      }
    }
    for (a = 0; a < WIM; a++) {
      for (i = 0; i < WIM; i++) {
        force[a] = -pressure * TT[a][i] * lsi->normal[i];
      }
    }

    /*  Jacobian contributions            */
    if (af->Assemble_Jacobian) {
      for (b = 0; b < dim; b++) {
        var = MESH_DISPLACEMENT1 + b;
        if (pd->v[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            for (p = 0; p < dim; p++) {
              for (q = 0; q < dim; q++) {
                d_force[p][var][j] -= pressure * lsi->normal[q] * dTT_dx[p][q][b][j];
              }
            }
          }
        }
      }
      /* For Lagrangian Mesh, add in the sensitivity to pressure  */
      if (pd->MeshMotion == LAGRANGIAN || pd->MeshMotion == DYNAMIC_LAGRANGIAN) {
        var = PRESSURE;
        if (pd->v[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            for (p = 0; p < dim; p++) {
              for (q = 0; q < dim; q++) {
                d_force[p][var][j] -= pressure * lsi->normal[q] * dTT_dp[p][q][j];
              }
            }
          }
        }
      }
      var = MASS_FRACTION;
      if (pd->v[pg->imtrx][var]) {
        for (p = 0; p < pd->Num_Dim; p++) {
          for (q = 0; q < pd->Num_Dim; q++) {
            for (w = 0; w < pd->Num_Species_Eqn; w++) {
              for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                d_force[p][MAX_VARIABLE_TYPES + w][j] -=
                    pressure * lsi->normal[q] * dTT_dc[p][q][w][j];
              }
            }
          }
        }
      }
      if (pd->v[pg->imtrx][POR_LIQ_PRES] &&
          (cr->MeshMotion == LAGRANGIAN || cr->MeshMotion == DYNAMIC_LAGRANGIAN)) {
        var = POR_LIQ_PRES;
        if (pd->v[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            for (p = 0; p < dim; p++) {
              for (q = 0; q < dim; q++) {
                d_force[p][var][j] -= pressure * lsi->normal[q] * dTT_dp_liq[p][q][j];
              }
            }
          }
        }
      }
      if (pd->v[pg->imtrx][POR_GAS_PRES] &&
          (cr->MeshMotion == LAGRANGIAN || cr->MeshMotion == DYNAMIC_LAGRANGIAN)) {
        var = POR_GAS_PRES;
        if (pd->v[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            for (p = 0; p < dim; p++) {
              for (q = 0; q < dim; q++) {
                d_force[p][var][j] -= pressure * lsi->normal[q] * dTT_dp_gas[p][q][j];
              }
            }
          }
        }
      }
      var = POR_POROSITY;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          for (p = 0; p < dim; p++) {
            for (q = 0; q < dim; q++) {
              d_force[p][var][j] -= pressure * lsi->normal[q] * dTT_dporosity[p][q][j];
            }
          }
        }
      }

      var = POR_SINK_MASS;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          for (p = 0; p < dim; p++) {
            for (q = 0; q < dim; q++) {
              d_force[p][var][j] -= pressure * lsi->normal[q] * dTT_dsink_mass[p][q][j];
            }
          }
        }
      }
      if ((cr->MeshMotion == LAGRANGIAN || cr->MeshMotion == DYNAMIC_LAGRANGIAN) &&
          pd->e[pg->imtrx][eqn]) {
        var = TEMPERATURE;
        if (pd->v[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            for (p = 0; p < dim; p++) {
              for (q = 0; q < dim; q++) {
                d_force[p][var][j] -= pressure * lsi->normal[q] * dTT_dT[p][q][j];
              }
            }
          }
        }
      }

      if ((cr->MeshMotion == LAGRANGIAN || cr->MeshMotion == DYNAMIC_LAGRANGIAN) &&
          pd->e[pg->imtrx][eqn]) {
        var = MAX_STRAIN;
        if (pd->v[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            for (p = 0; p < dim; p++) {
              for (q = 0; q < dim; q++) {
                d_force[p][var][j] += pressure * lsi->normal[q] * dTT_dmax_strain[p][q][j];
              }
            }
          }
        }
      }

      var = ls->var;
      if (pd->v[pg->imtrx][var]) {
        for (p = 0; p < pd->Num_Dim; p++) {
          for (q = 0; q < pd->Num_Dim; q++) {
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              d_force[p][var][j] -= pressure * (lsi->d_normal_dF[p][j] * TT[p][q]);
            }
          }
        }
      }
    } /*  if Jacobian  */

  } /* end of STRESS_TENSOR */
  else {
    GOMA_EH(GOMA_ERROR, "Invalid LS_FLOW_PRESSURE boundary condition flag (-1|0|1)\n");
  }

#if 0
fprintf(stderr,"pf %g %g %g %d %g %g\n",fv->x[0],fv->x[1], lsi->delta, elem_sign, force[0],force[1]);
#endif

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

          source = phi_i * force[a] * lsi->delta * elem_sign;

          source *= det_J * wt * h3;
          source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

          /* J_m_F
           */
          var = ls->var;
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
          source = phi_i * force[a] * lsi->delta * elem_sign;
          source *= det_J * wt * h3;
          source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
          lec->R[LEC_R_INDEX(peqn, ii)] += source;
        }
      }
    }
  }

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
           * J_m_F
           */
          var = ls->var;
          if (pd->v[pg->imtrx][var]) {
            pvar = upd->vp[pg->imtrx][var];

            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              phi_j = bf[var]->phi[j];

              source = phi_i * elem_sign *
                       (lsi->delta * d_force[a][var][j] + lsi->d_delta_dF[j] * force[a]);

              source *= det_J * wt * h3;
              source *= pd->etm[pg->imtrx][eqn][LOG2_SOURCE];

              lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
            }
          }
          /* extra dependencies for bcflag = 0 or 1	*/
          if (bcflag == 0 || bcflag == 1) {
            /*
             * J_m_Pressure if bc_Int == 0
             */

            var = PRESSURE;
            if (pd->v[pg->imtrx][var]) {
              pvar = upd->vp[pg->imtrx][var];

              for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                phi_j = bf[var]->phi[j];

                source = phi_i * d_force[a][var][j] * elem_sign * lsi->delta;

                source *= det_J * wt * h3;
                source *= pd->etm[pg->imtrx][eqn][LOG2_SOURCE];

                lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
              }
            }

            var = TEMPERATURE;
            if (pd->v[pg->imtrx][var]) {
              pvar = upd->vp[pg->imtrx][var];

              for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                phi_j = bf[var]->phi[j];

                source = phi_i * d_force[a][var][j] * elem_sign * lsi->delta;

                source *= det_J * wt * h3;
                source *= pd->etm[pg->imtrx][eqn][LOG2_SOURCE];

                lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
              }
            }

            var = BOND_EVOLUTION;
            if (pd->v[pg->imtrx][var]) {
              pvar = upd->vp[pg->imtrx][var];

              for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                phi_j = bf[var]->phi[j];

                source = phi_i * d_force[a][var][j] * elem_sign * lsi->delta;

                source *= det_J * wt * h3;
                source *= pd->etm[pg->imtrx][eqn][LOG2_SOURCE];

                lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
              }
            }

            var = MASS_FRACTION;
            if (pd->v[pg->imtrx][var]) {
              for (w = 0; w < pd->Num_Species_Eqn; w++) {
                pvar = MAX_PROB_VAR + w;
                for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                  phi_j = bf[var]->phi[j];

                  source = phi_i * d_force[a][MAX_VARIABLE_TYPES + w][j] * elem_sign * lsi->delta;

                  source *= det_J * wt * h3;
                  source *= pd->etm[pg->imtrx][eqn][LOG2_SOURCE];

                  lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
                }
              }
            }

            var = POR_LIQ_PRES;
            if (pd->v[pg->imtrx][var]) {
              pvar = upd->vp[pg->imtrx][var];

              for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                phi_j = bf[var]->phi[j];

                source = phi_i * d_force[a][var][j] * elem_sign * lsi->delta;

                source *= det_J * wt * h3;
                source *= pd->etm[pg->imtrx][eqn][LOG2_SOURCE];

                lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
              }
            }

            var = POR_GAS_PRES;
            if (pd->v[pg->imtrx][var]) {
              pvar = upd->vp[pg->imtrx][var];

              for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                phi_j = bf[var]->phi[j];

                source = phi_i * d_force[a][var][j] * elem_sign * lsi->delta;

                source *= det_J * wt * h3;
                source *= pd->etm[pg->imtrx][eqn][LOG2_SOURCE];

                lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
              }
            }

            var = POR_POROSITY;
            if (pd->v[pg->imtrx][var]) {
              pvar = upd->vp[pg->imtrx][var];

              for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                phi_j = bf[var]->phi[j];

                source = phi_i * d_force[a][var][j] * elem_sign * lsi->delta;

                source *= det_J * wt * h3;
                source *= pd->etm[pg->imtrx][eqn][LOG2_SOURCE];

                lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
              }
            }

            var = POR_SINK_MASS;
            if (pd->v[pg->imtrx][var]) {
              pvar = upd->vp[pg->imtrx][var];

              for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                phi_j = bf[var]->phi[j];

                source = phi_i * d_force[a][var][j] * elem_sign * lsi->delta;

                source *= det_J * wt * h3;
                source *= pd->etm[pg->imtrx][eqn][LOG2_SOURCE];

                lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
              }
            }

            for (b = 0; b < dim; b++) {
              var = MESH_DISPLACEMENT1 + b;

              if (pd->v[pg->imtrx][var]) {
                pvar = upd->vp[pg->imtrx][var];

                for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                  phi_j = bf[var]->phi[j];

                  source = phi_i * d_force[a][var][j] * elem_sign * lsi->delta;

                  source *= det_J * wt * h3;
                  source *= pd->etm[pg->imtrx][eqn][LOG2_SOURCE];

                  lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
                }
              }
            }
            for (b = 0; b < dim; b++) {
              var = VELOCITY1 + b;

              if (pd->v[pg->imtrx][var]) {
                pvar = upd->vp[pg->imtrx][var];

                for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                  phi_j = bf[var]->phi[j];

                  source = phi_i * d_force[a][var][j] * elem_sign * lsi->delta;

                  source *= det_J * wt * h3;
                  source *= pd->etm[pg->imtrx][eqn][LOG2_SOURCE];

                  lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
                }
              }
            }

            for (b = 0; b < dim; b++) {
              var = VORT_DIR1 + b;

              if (pd->v[pg->imtrx][var]) {
                pvar = upd->vp[pg->imtrx][var];

                for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                  phi_j = bf[var]->phi[j];

                  source = phi_i * d_force[a][var][j] * elem_sign * lsi->delta;

                  source *= det_J * wt * h3;
                  source *= pd->etm[pg->imtrx][eqn][LOG2_SOURCE];

                  lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
                }
              }
            }
          }
        }
      }
    }
  }
  return (1);
}

int assemble_precoil_source(const double p[]) {
  int i, j, a, ii, ledof;
  int eqn, peqn, var, pvar;

  struct Basis_Functions *bfm;
  double wt, det_J, h3, phi_i, phi_j;
  double pressure;
  double source;

  double dprtemp, theta;
  double pabl_c0 = 0.0, pabl_c1 = 0.0, pabl_c2 = 0.0, pabl_c3 = 0.0;
  double T_boil, T_scale, P_scale;

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

  /* Mike Kanouff's curve-fit for ablation pressure in Pascals               */
  /*                                              assume iron if T_boil>2000 */
  /*                                                      ice if T_boil<2000 */
  T_boil = mp->melting_point_solidus; /* T_boil stored as T_solidus in mat-file */
  P_scale = p[5];
  T_scale = p[6];

  if (T_boil > 2000.0 * T_scale) {
    theta = fv->T - T_boil;
    if (theta > 0.0 && theta <= 170.0 * T_scale) {
      pabl_c0 = 0.0; /* iron coefficients */
      pabl_c1 = 1.8272e-4 * 1.0133e5 * (1.0 / T_scale);
      pabl_c2 = -1.9436e-6 * 1.0133e5 * (1.0 / T_scale) * (1.0 / T_scale);
      pabl_c3 = 1.5732e-8 * 1.0133e5 * (1.0 / T_scale) * (1.0 / T_scale) * (1.0 / T_scale);
    }
    if (theta > 170.0 * T_scale) {
      pabl_c0 = 0.0; /* iron coefficients */
      pabl_c1 = -5.7333e-4 * 1.0133e5 * (1.0 / T_scale);
      pabl_c2 = 4.5500e-6 * 1.0133e5 * (1.0 / T_scale) * (1.0 / T_scale);
      pabl_c3 = 2.3022e-9 * 1.0133e5 * (1.0 / T_scale) * (1.0 / T_scale) * (1.0 / T_scale);
    }
    /*pabl_c0 =  0.0;          */ /* MKS coefficients for iron */
    /*pabl_c1 =  3.723086e+02*(1.0/T_scale);*/
    /*pabl_c2 =  -6.328050e-02*(1.0/T_scale)*(1.0/T_scale);*/
    /*pabl_c3 =  5.559470e-04*(1.0/T_scale)*(1.0/T_scale)*(1.0/T_scale);*/
  } else {
    pabl_c0 = 0.0; /* MKS coefficients for ice  */
    pabl_c1 = 3.294180e+03 * (1.0 / T_scale);
    pabl_c2 = -7.726940e+00 * (1.0 / T_scale) * (1.0 / T_scale);
    pabl_c3 = 5.480973e-01 * (1.0 / T_scale) * (1.0 / T_scale) * (1.0 / T_scale);
  }

  /* Calculate ablation pressure */
  theta = fv->T - T_boil;

  if (theta < 0.) {
    theta = 0.;
    pressure = 0.;
    dprtemp = 0;
  } else {
    pressure = P_scale *
               (pabl_c0 + pabl_c1 * theta + pabl_c2 * pow(theta, 2.0) + pabl_c3 * pow(theta, 3.0));
    dprtemp = P_scale * (pabl_c1 + 2. * pabl_c2 * theta + 3. * pabl_c3 * pow(theta, 2.0));
  }
  /* END Kanouff Curve fir for pressure */

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

          source = phi_i * lsi->normal[a] * pressure * lsi->delta * ls->Elem_Sign;

          source *= det_J * wt * h3;
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

          source = phi_i * lsi->normal[a] * pressure * lsi->delta * ls->Elem_Sign;

          source *= det_J * wt * h3;
          source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

          lec->R[LEC_R_INDEX(peqn, ii)] += source;
        }
      }
    }
  }

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
           * J_m_T
           */
          var = TEMPERATURE;

          if (pd->v[pg->imtrx][var]) {
            pvar = upd->vp[pg->imtrx][var];

            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              phi_j = bf[var]->phi[j];

              source = phi_i * lsi->normal[a] * lsi->delta * ls->Elem_Sign * dprtemp * phi_j;

              source *= det_J * wt * h3;
              source *= pd->etm[pg->imtrx][eqn][LOG2_SOURCE];

              lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
            }
          }
          /*
           * J_m_F
           */
          var = FILL;

          if (pd->v[pg->imtrx][var]) {
            pvar = upd->vp[pg->imtrx][var];

            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              phi_j = bf[var]->phi[j];

              source = phi_i * pressure * ls->Elem_Sign *
                       (lsi->delta * lsi->d_normal_dF[a][j] + lsi->d_delta_dF[j] * lsi->normal[a]);

              source *= det_J * wt * h3;
              source *= pd->etm[pg->imtrx][eqn][LOG2_SOURCE];

              lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
            }
          }
        }
      }
    }
  }
  return (1);
}

#if 0
int
assemble_uvw_source ( int eqn, double val )
{
  int i,j, a, ii,ledof;
  int peqn, var, pvar;
  int b, p;

  struct Basis_Functions *bfm;
  double wt, det_J, h3, phi_i, phi_j;

  double source;

  int xfem_active, extended_dof, base_interp, base_dof;
  double F_i;
  int sign = 1;

  /*
   * Variables for stress tensor and derivative
   */

  dbl Pi[DIM][DIM];
  STRESS_DEPENDENCE_STRUCT d_Pi_struct;
  STRESS_DEPENDENCE_STRUCT *d_Pi = &d_Pi_struct;

  if ( ! pd->e[pg->imtrx][eqn] )
    {
      return(0);
    }

  /* For LS "Dirichlet" conditions, die if attempted on
     non-discontinuous interpolations */
  if ( pd->i[pg->imtrx][eqn] != I_Q1_XV &&
       pd->i[pg->imtrx][eqn] != I_Q2_XV &&
       pd->i[pg->imtrx][eqn] != I_Q1_GP  &&
       pd->i[pg->imtrx][eqn] != I_Q2_GP  &&
       pd->i[pg->imtrx][eqn] != I_Q1_GN  &&
       pd->i[pg->imtrx][eqn] != I_Q2_GN  &&
       pd->i[pg->imtrx][eqn] != I_Q1_G  &&
       pd->i[pg->imtrx][eqn] != I_Q2_G )
    {
      GOMA_EH(GOMA_ERROR, "LS_UVW requires discontinuous velocity enrichment (Q1_XV, Q2_XV, Q1_G, Q2_G)");
    }

  a = eqn - R_MOMENTUM1;

  wt = fv->wt;
  h3 = fv->h3;

  if ( ls->on_sharp_surf ) /* sharp interface */
	{det_J = fv->sdet;}
  else              /* diffuse interface */
    	{det_J = bf[eqn]->detJ; }

  /*
   * non-xfem dofs will receive nobc and xfem
   * dofs will receive SIC
   */

  /* compute stress tensor and its derivatives */
  fluid_stress( Pi, d_Pi );

  /*
   * Wesiduals ________________________________________________________________________________
   */
  if ( af->Assemble_Residual )
    {
      peqn = upd->ep[pg->imtrx][eqn];
      bfm = bf[eqn];

      for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++)
        {
          ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][i];

          if (ei[pg->imtrx]->active_interp_ledof[ledof])
            {
              ii = ei[pg->imtrx]->lvdof_to_row_lvdof[eqn][i];

              phi_i = bfm->phi[i];

              xfem_dof_state( i, pd->i[pg->imtrx][eqn], ei[pg->imtrx]->ielem_shape,
                              &xfem_active, &extended_dof, &base_interp, &base_dof );
	      F_i = dof_distance( eqn, i );

	      if ( !extended_dof && !sign_change( F_i, (double) ls->Elem_Sign ) )
                {
                  source = 0.;
                  for ( p=0; p<VIM; p++) source += Pi[p][a] * lsi->normal[p];
                  source *= phi_i * lsi->delta * ls->Elem_Sign;
                  source *= -det_J * wt * h3;
                  source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

                  lec->R[LEC_R_INDEX(peqn,ii)] += source;
                }
	      if ( extended_dof && sign_change( F_i, (double) ls->Elem_Sign ) )
                {
                  source = phi_i * lsi->delta * (fv->v[a] - val);
                  source *= det_J * wt * h3 * BIG_PENALTY * sign;
                  source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

                  lec->R[LEC_R_INDEX(peqn,ii)] += source;
                }
            }
        }
    }

  /*
   * Yacobian terms...
   */

  if( af->Assemble_Jacobian )
    {
      peqn = upd->ep[pg->imtrx][eqn];
      bfm = bf[eqn];

      for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++)
        {
          ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][i];

          if (ei[pg->imtrx]->active_interp_ledof[ledof])
            {
              ii = ei[pg->imtrx]->lvdof_to_row_lvdof[eqn][i];

              phi_i = bfm->phi[i];

              xfem_dof_state( i, pd->i[pg->imtrx][eqn], ei[pg->imtrx]->ielem_shape,
                              &xfem_active, &extended_dof, &base_interp, &base_dof );
	      F_i = dof_distance( eqn, i );

	      if ( !extended_dof && !sign_change( F_i, (double) ls->Elem_Sign ) )
                {
                  /*
                   * J_m_m
                   */
                  for ( b=0; b<WIM; b++)
                    {
                      var = VELOCITY1 + b;

                      if ( pd->v[pg->imtrx][var] )
                        {
                          pvar = upd->vp[pg->imtrx][var];

                          for( j=0; j<ei[pg->imtrx]->dof[var]; j++)
                            {
                              phi_j = bf[var]->phi[j];

                              source = 0.;
                              for ( p=0; p<VIM; p++) source += d_Pi->v[b][p][a][j] * lsi->normal[p];
                              source *= phi_i * lsi->delta * ls->Elem_Sign;
                              source *= -det_J * wt * h3;
                              source *= pd->etm[pg->imtrx][eqn][LOG2_SOURCE];

                              lec->J[LEC_J_INDEX(peqn,pvar,ii,j)] += source;
                            }
                        }
                    }
                  /*
                   * J_m_P
                   */
                  var = PRESSURE;

                  if ( pd->v[pg->imtrx][var] )
                    {
                      pvar = upd->vp[pg->imtrx][var];

                      for( j=0; j<ei[pg->imtrx]->dof[var]; j++)
                        {
                          phi_j = bf[var]->phi[j];

                          source = 0.;
                          for ( p=0; p<VIM; p++) source += d_Pi->P[p][a][j] * lsi->normal[p];
                          source *= phi_i * lsi->delta * ls->Elem_Sign;
                          source *= -det_J * wt * h3;
                          source *= pd->etm[pg->imtrx][eqn][LOG2_SOURCE];

                          lec->J[LEC_J_INDEX(peqn,pvar,ii,j)] += source;
                        }
                    }
                  /*
                   * J_m_F
                   */
                  var = FILL;

                  if ( pd->v[pg->imtrx][var] )
                    {
                      pvar = upd->vp[pg->imtrx][var];

                      for( j=0; j<ei[pg->imtrx]->dof[var]; j++)
                        {
                          phi_j = bf[var]->phi[j];

                          source = 0.;
                          for ( p=0; p<VIM; p++) source += d_Pi->F[p][a][j] * lsi->normal[p] * lsi->delta +
                                                           Pi[p][a] * lsi->d_normal_dF[p][j] * lsi->delta +
                                                           Pi[p][a] * lsi->normal[p] * lsi->d_delta_dF[j];
                          
                          source *= phi_i * ls->Elem_Sign;
                          source *= -det_J * wt * h3;
                          source *= pd->etm[pg->imtrx][eqn][LOG2_SOURCE];

                          lec->J[LEC_J_INDEX(peqn,pvar,ii,j)] += source;
                        }
                    }
                }

	      if ( extended_dof && sign_change( F_i, (double) ls->Elem_Sign ) )
                {
                  /*
                   * J_m_m
                   */
                  var = eqn;

                  if ( pd->v[pg->imtrx][var] )
                    {
                      pvar = upd->vp[pg->imtrx][var];

                      for( j=0; j<ei[pg->imtrx]->dof[var]; j++)
                        {
                          phi_j = bf[var]->phi[j];

                          source = phi_i * phi_j * lsi->delta;
                          source *= det_J * wt * h3 * BIG_PENALTY * sign;
                          source *= pd->etm[pg->imtrx][eqn][LOG2_SOURCE];

                          lec->J[LEC_J_INDEX(peqn,pvar,ii,j)] += source;
                        }
                    }
                  /*
                   * J_m_F
                   */
                  var = FILL;
                  if( pd->v[pg->imtrx][var])
                    {
                      pvar = upd->vp[pg->imtrx][var];
	      	      
	      	      /* path dependence integral */
                      if ( ls->on_sharp_surf && !ls->AdaptIntegration)
                	{
                	  for( j=0; j<ei[pg->imtrx]->dof[var]; j++)
                	    {
                	      phi_j = bf[var]->phi[j];

                	      source = 0.;
/*
                	      for ( p=0; p < VIM; p++ )
                		{
                		  source -= phi_i * phi_j * lsi->normal[p] * fv->grad_v[p][a];
                		}
*/

                	      for ( p=0; p < VIM; p++ )
                		{
                		  source -= phi_j * lsi->normal[p] * bfm->grad_phi[i][p] * 
	      	 			   ( fv->v[a] - val );
                		}

                	      source *= lsi->gfmaginv;
                	      source *= det_J * wt * h3 * BIG_PENALTY * sign;
                	      source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

                              lec->J[LEC_J_INDEX(peqn,pvar,ii,j)] += source;
                	    }
                	}
                      else
                	{
                	  /* diffuse interface version of path dependence integral */
                	  for( j=0; j<ei[pg->imtrx]->dof[var]; j++)
                	    {
                	      source = phi_i * ( fv->v[a] - val );
                	      source *= lsi->d_delta_dF[j];
                	      source *= det_J * wt * h3 * BIG_PENALTY;
                	      source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

                              lec->J[LEC_J_INDEX(peqn,pvar,ii,j)] += source;
                	    }
                	}
	      	    }
                }
            }
        }
    }
  return ( 1 );
}
#endif
#if 1
int assemble_uvw_source(int eqn, double val) {
  int i, j, a, ii, ledof;
  int peqn, var, pvar;
  int b, p;
  int xfem_active, extended_dof, base_interp, base_dof;
  int incomplete, other_side;
  int apply_NOBC[MDE], apply_SIC[MDE];

  struct Basis_Functions *bfm;
  double wt, det_J, h3, phi_i, phi_j;

  double source;
  double F_i;
  int sign = 1;

  /*
   * Variables for stress tensor and derivative
   */

  dbl Pi[DIM][DIM];
  STRESS_DEPENDENCE_STRUCT d_Pi_struct;
  STRESS_DEPENDENCE_STRUCT *d_Pi = &d_Pi_struct;

  if (!pd->e[pg->imtrx][eqn]) {
    return (0);
  }

  /* For LS "Dirichlet" conditions, die if attempted on
     non-discontinuous interpolations */
  if (pd->i[pg->imtrx][eqn] != I_Q1_XV && pd->i[pg->imtrx][eqn] != I_Q2_XV &&
      pd->i[pg->imtrx][eqn] != I_Q1_XG && pd->i[pg->imtrx][eqn] != I_Q2_XG &&
      pd->i[pg->imtrx][eqn] != I_Q1_GP && pd->i[pg->imtrx][eqn] != I_Q2_GP &&
      pd->i[pg->imtrx][eqn] != I_Q1_GN && pd->i[pg->imtrx][eqn] != I_Q2_GN &&
      pd->i[pg->imtrx][eqn] != I_Q1_G && pd->i[pg->imtrx][eqn] != I_Q2_G) {
    GOMA_EH(GOMA_ERROR, "LS_UVW requires discontinuous velocity enrichment");
  }

  a = eqn - R_MOMENTUM1;

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
   * non-xfem dofs will receive nobc and xfem
   * dofs will receive SIC
   */

  if (pd->i[pg->imtrx][eqn] == I_Q1_XG || pd->i[pg->imtrx][eqn] == I_Q2_XG) {
    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      apply_NOBC[i] = FALSE;
      apply_SIC[i] = TRUE;
    }
  } else {
    /* compute stress tensor and its derivatives */
    fluid_stress(Pi, d_Pi);

    /* determine if NOBC or SIC or both are to be applied */
    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][i];

      if (ei[pg->imtrx]->active_interp_ledof[ledof]) {
        xfem_dof_state(i, pd->i[pg->imtrx][eqn], ei[pg->imtrx]->ielem_shape, &xfem_active,
                       &extended_dof, &base_interp, &base_dof);
        incomplete = dof_incomplete(base_dof, ei[pg->imtrx]->ielem_type, base_interp,
                                    ei[pg->imtrx]->ielem_shape);
        F_i = dof_distance(eqn, i);
        other_side = sign_change(F_i, (double)ls->Elem_Sign);

        if (extended_dof) {
          if (other_side) {
            if (incomplete) {
              /* extended_dof&&other_side&&incomplete -> NOBC + SIC */
              apply_NOBC[i] = TRUE;
              apply_SIC[i] = TRUE;
            } else {
              /* extended_dof&&other_side&&!incomplete -> NOBC */
              apply_NOBC[i] = TRUE;
              apply_SIC[i] = FALSE;
            }
          } else {
            /* extended_dof&&!other_side -> do nothing */
            apply_NOBC[i] = FALSE;
            apply_SIC[i] = FALSE;
          }
        } else {
          if (other_side) {
            /* !extended_dof&&other_side -> do nothing */
            apply_NOBC[i] = FALSE;
            apply_SIC[i] = FALSE;
          } else {
            /* !extended_dof&&!other_side -> NOBC */
            apply_NOBC[i] = TRUE;
            apply_SIC[i] = FALSE;
          }
        }
      }
    }
  }

  /* finite difference calculation of path dependencies for
     subelement integration
   */
  if (ls->CalcSurfDependencies) {
    peqn = upd->ep[pg->imtrx][eqn];
    bfm = bf[eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][i];

      if (ei[pg->imtrx]->active_interp_ledof[ledof]) {
        ii = ei[pg->imtrx]->lvdof_to_row_lvdof[eqn][i];

        phi_i = bfm->phi[i];

        if (apply_NOBC[i]) {
          source = 0.;
          for (p = 0; p < VIM; p++)
            source += Pi[p][a] * lsi->normal[p];
          source *= phi_i * lsi->delta * ls->Elem_Sign;
          source *= -det_J * wt * h3;
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
        if (apply_SIC[i]) {
          source = phi_i * lsi->delta * (fv->v[a] - val);
          source *= det_J * wt * h3 * BIG_PENALTY * sign;
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
    peqn = upd->ep[pg->imtrx][eqn];
    bfm = bf[eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][i];

      if (ei[pg->imtrx]->active_interp_ledof[ledof]) {
        ii = ei[pg->imtrx]->lvdof_to_row_lvdof[eqn][i];

        phi_i = bfm->phi[i];

        if (apply_NOBC[i]) {
          source = 0.;
          for (p = 0; p < VIM; p++)
            source += Pi[p][a] * lsi->normal[p];
          source *= phi_i * lsi->delta * ls->Elem_Sign;
          source *= -det_J * wt * h3;
          source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

          lec->R[LEC_R_INDEX(peqn, ii)] += source;
        }
        if (apply_SIC[i]) {
          source = phi_i * lsi->delta * (fv->v[a] - val);
          source *= det_J * wt * h3 * BIG_PENALTY * sign;
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
    peqn = upd->ep[pg->imtrx][eqn];
    bfm = bf[eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][i];

      if (ei[pg->imtrx]->active_interp_ledof[ledof]) {
        ii = ei[pg->imtrx]->lvdof_to_row_lvdof[eqn][i];

        phi_i = bfm->phi[i];

        xfem_dof_state(i, pd->i[pg->imtrx][eqn], ei[pg->imtrx]->ielem_shape, &xfem_active,
                       &extended_dof, &base_interp, &base_dof);
        F_i = dof_distance(eqn, i);

        if (apply_NOBC[i]) {
          /*
           * J_m_v
           */
          for (b = 0; b < WIM; b++) {
            var = VELOCITY1 + b;

            if (pd->v[pg->imtrx][var]) {
              pvar = upd->vp[pg->imtrx][var];

              for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                phi_j = bf[var]->phi[j];

                source = 0.;
                for (p = 0; p < VIM; p++)
                  source += d_Pi->v[p][a][b][j] * lsi->normal[p];
                source *= phi_i * lsi->delta * ls->Elem_Sign;
                source *= -det_J * wt * h3;
                source *= pd->etm[pg->imtrx][eqn][LOG2_SOURCE];

                lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
              }
            }
          }
          /*
           * J_m_vd
           */
          for (b = 0; b < DIM; b++) {
            var = VORT_DIR1 + b;

            if (pd->v[pg->imtrx][var]) {
              pvar = upd->vp[pg->imtrx][var];

              for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                phi_j = bf[var]->phi[j];

                source = 0.;
                for (p = 0; p < VIM; p++)
                  source += d_Pi->vd[p][a][b][j] * lsi->normal[p];
                source *= phi_i * lsi->delta * ls->Elem_Sign;
                source *= -det_J * wt * h3;
                source *= pd->etm[pg->imtrx][eqn][LOG2_SOURCE];

                lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
              }
            }
          }
          /*
           * J_m_P
           */
          var = PRESSURE;

          if (pd->v[pg->imtrx][var]) {
            pvar = upd->vp[pg->imtrx][var];

            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              phi_j = bf[var]->phi[j];

              source = 0.;
              for (p = 0; p < VIM; p++)
                source += d_Pi->P[p][a][j] * lsi->normal[p];
              source *= phi_i * lsi->delta * ls->Elem_Sign;
              source *= -det_J * wt * h3;
              source *= pd->etm[pg->imtrx][eqn][LOG2_SOURCE];

              lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
            }
          }
          /*
           * J_m_F
           */
          var = FILL;

          if (pd->v[pg->imtrx][var]) {
            pvar = upd->vp[pg->imtrx][var];

            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              phi_j = bf[var]->phi[j];

              source = 0.;
              for (p = 0; p < VIM; p++)
                source += d_Pi->F[p][a][j] * lsi->normal[p] * lsi->delta +
                          Pi[p][a] * lsi->d_normal_dF[p][j] * lsi->delta +
                          Pi[p][a] * lsi->normal[p] * lsi->d_delta_dF[j];

              source *= phi_i * ls->Elem_Sign;
              source *= -det_J * wt * h3;
              source *= pd->etm[pg->imtrx][eqn][LOG2_SOURCE];

              lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
            }
          }
        }

        if (apply_SIC[i]) {
          /*
           * J_m_m
           */
          var = eqn;

          if (pd->v[pg->imtrx][var]) {
            pvar = upd->vp[pg->imtrx][var];

            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              phi_j = bf[var]->phi[j];

              source = phi_i * phi_j * lsi->delta;
              source *= det_J * wt * h3 * BIG_PENALTY * sign;
              source *= pd->etm[pg->imtrx][eqn][LOG2_SOURCE];

              lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
            }
          }
          /*
           * J_m_F
           */
          var = FILL;
          if (pd->v[pg->imtrx][var]) {
            pvar = upd->vp[pg->imtrx][var];

            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              source = phi_i * (fv->v[a] - val);
              source *= lsi->d_delta_dF[j];
              source *= det_J * wt * h3 * BIG_PENALTY;
              source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

              lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
            }
          }
        }
      }
    }
  }
  return (1);
}
#endif

int assemble_extension_velocity_path_dependence(void) {
  int a, i, j, ii, ledof;
  int eqn, peqn, var, pvar;
  int dim;

  struct Basis_Functions *bfm;
  double wt, det_J, h3, phi_i;

  double source = 0.;
  double resid;
  double *grad_F;

  eqn = R_EXT_VELOCITY;
  if (!pd->e[pg->imtrx][eqn]) {
    return (0);
  }

  dim = pd->Num_Dim;

  wt = fv->wt;
  h3 = fv->h3;
  if (ls->on_sharp_surf) /* sharp interface */
  {
    det_J = fv->sdet;
  } else /* diffuse interface */
  {
    det_J = bf[eqn]->detJ;
  }
  if (pfd == NULL) {
    grad_F = fv->grad_F;
  } else {
    grad_F = fv->grad_pF[ls->var - PHASE1];
  }

#define GRADF_GRADEXTV  1
#define NORMAL_GRADEXTV 0
#if GRADF_GRADEXTV
  resid = 0.;
  for (a = 0; a < dim; a++)
    resid += grad_F[a] * fv->grad_ext_v[a];
#endif
#if NORMAL_GRADEXTV
  resid = 0.;
  for (a = 0; a < VIM; a++)
    resid += lsi->normal[a] * fv->grad_ext_v[a];
#endif
  /*
   * Yacobian terms...
   */

  if (af->Assemble_Jacobian) {
    peqn = upd->ep[pg->imtrx][eqn];
    bfm = bf[eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {

      ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][i];

      if (ei[pg->imtrx]->active_interp_ledof[ledof]) {

        ii = ei[pg->imtrx]->lvdof_to_row_lvdof[eqn][i];

        phi_i = bfm->phi[i];

        /* path dependence integral from volume equation */

        /*
         * J_ls_F
         */
        var = ls->var;

        if (pd->v[pg->imtrx][var]) {
          pvar = upd->vp[pg->imtrx][var];

          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {

            source = resid * (2. * lsi->d_H_dF[j]) * phi_i;
            source *= det_J * wt * h3;

            lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
          }
        }
      }
    }
  }

  return (1);
}

int assemble_fill_path_dependence(void) {
  int i, j, ii, ledof;
  int eqn, peqn, var, pvar;

  struct Basis_Functions *bfm;
  double wt, det_J, h3, phi_i;

  double source = 0.;

  eqn = R_FILL;
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
   * Yacobian terms...
   */

  if (af->Assemble_Jacobian) {
    peqn = upd->ep[pg->imtrx][eqn];
    bfm = bf[eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {

      ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][i];

      if (ei[pg->imtrx]->active_interp_ledof[ledof]) {

        ii = ei[pg->imtrx]->lvdof_to_row_lvdof[eqn][i];

        phi_i = bfm->phi[i];

        /* path dependence integral from volume equation */

        /*
         * J_ls_F
         */
        var = FILL;

        if (pd->v[pg->imtrx][var]) {
          pvar = upd->vp[pg->imtrx][var];

          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            switch (tran->Fill_Weight_Fcn) {
            case FILL_WEIGHT_G:
            case FILL_WEIGHT_GLS:

              source = (lsi->gfmag - 1.) * (2. * lsi->d_H_dF[j]) * phi_i;

              break;
            default:
              GOMA_EH(GOMA_ERROR, "Unknown Fill_Weight_Fcn");
            }

            source *= det_J * wt * h3;
            source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

            lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
          }
        }
      }
    }
  }

  return (1);
}

int assemble_energy_path_dependence(
    double time,            /* present time value */
    double tt,              /* parameter to vary time integration from
                               explicit (tt = 1) to implicit (tt = 0) */
    double dt,              /* current time step size */
    const PG_DATA *pg_data) /* average element and velocity infor for SUPG and PSPG */
{
  int eqn, var, peqn, pvar, dim, p, i, j, status;

  dbl T_dot; /* Temperature derivative wrt time. */

  dbl q[DIM]; /* Heat flux vector. */
  dbl rho;    /* Density */
  dbl Cp;     /* Heat capacity. */
  dbl h;      /* Heat source. */

  dbl mass;      /* For terms and their derivatives */
  dbl advection; /* For terms and their derivatives */
  dbl diffusion;
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
  dbl h_elem = 0;
  dbl supg;
  const dbl *vcent, *hsquared;

  /*
   * Interpolation functions for variables and some of their derivatives.
   */

  dbl h3; /* Volume element (scale factors). */
  dbl det_J;
  dbl wt;

  dbl *grad_T = fv->grad_T;

  double vconv[MAX_PDIM];     /*Calculated convection velocity */
  double vconv_old[MAX_PDIM]; /*Calculated convection velocity at previous time*/

  int err;

  int sign = ls->Elem_Sign;
  double energy_residual;

  /*   static char yo[] = "assemble_energy";*/

  status = 0;

  eqn = R_ENERGY;

  /*
   * Bail out fast if there's nothing to do...
   */

  if (!pd->e[pg->imtrx][eqn]) {
    return (status);
  }

  dim = pd->Num_Dim;

  wt = fv->wt;
  h3 = fv->h3;

  if (ls->on_sharp_surf) /* sharp interface */
  {
    det_J = fv->sdet;
  } else /* diffuse interface */
  {
    det_J = bf[eqn]->detJ;
  }

  supg = 0.0;

  if (mp->Ewt_funcModel == GALERKIN) {
    supg = 0.;
  } else if (mp->Ewt_funcModel == SUPG) {
    supg = mp->Ewt_func;
  }

  if (supg != 0.0) {
    h_elem = 0.;
    vcent = pg_data->v_avg;
    hsquared = pg_data->hsquared;

    for (p = 0; p < dim; p++) {
      h_elem += vcent[p] * vcent[p] * hsquared[p];
    }
    h_elem = sqrt(h_elem) / 2.;
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

  rho = density(NULL, time);

  /* CHECK FOR REMOVAL */
  conductivity(NULL, time);

  Cp = heat_capacity(NULL, time);

  h = heat_source(NULL, time, tt, dt);

  if (pd->TimeIntegration != STEADY) {
    T_dot = fv_dot->T;
  } else {
    T_dot = 0.0;
  }

  heat_flux(q, NULL, time);

  /* get the convection velocity (it's different for arbitrary and
     lagrangian meshes) */
  if (cr->MeshMotion == ARBITRARY || cr->MeshMotion == LAGRANGIAN ||
      cr->MeshMotion == DYNAMIC_LAGRANGIAN) {
    err = get_convection_velocity(vconv, vconv_old, NULL, dt, tt);
    GOMA_EH(err, "Error in calculating effective convection velocity");
  } else if (cr->MeshMotion == TOTAL_ALE) {
    err = get_convection_velocity_rs(vconv, vconv_old, NULL, dt, tt);
    GOMA_EH(err, "Error in calculating effective convection velocity_rs");
  }

  /*
   * Residuals___________________________________________________________
   */

  if (af->Assemble_Jacobian) {
    eqn = R_ENERGY;
    peqn = upd->ep[pg->imtrx][eqn];
    var = TEMPERATURE;
    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
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
          wt_func += supg * h_elem * vconv[p] * bf[eqn]->grad_phi[i][p];
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

      source = 0.;
      if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
        source += phi_i * h * det_J * wt;
        source *= h3;
        source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
      }

      energy_residual = mass + advection + diffusion + source;

      var = FILL;
      pvar = upd->vp[pg->imtrx][var];
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += lsi->d_H_dF[j] * energy_residual * sign;
      }
    }
  }
  return (status);
}

int assemble_momentum_path_dependence(dbl time, /* currentt time step */
                                      dbl tt,   /* parameter to vary time integration from
                                                   explicit (tt = 1) to implicit (tt = 0) */
                                      dbl dt,   /* current time step size */
                                      const PG_DATA *pg_data) {
  int i, j, a, p;
  int ledof, eqn, var, ii, peqn, pvar;
  int status;
  struct Basis_Functions *bfm;

  dbl zero[3] = {0.0, 0.0, 0.0}; /* A zero array, for convenience. */
  dbl *v_dot;                    /* time derivative of velocity field. */
  dbl *x_dot;                    /* current position field derivative wrt time. */

  dbl h3;    /* Volume element (scale factors). */
  dbl det_J; /* determinant of element Jacobian */

  /* field variables */
  dbl *grad_v[DIM];
  dbl *v = fv->v;

  dbl rho; /* Density. */

  dbl f[DIM]; /* Body force. */

  dbl mass; /* For terms and their derivatives */
  dbl advection;
  dbl porous;
  dbl diffusion;
  dbl source;

  /*
   * Galerkin weighting functions for i-th and a-th momentum residuals
   * and some of their derivatives...
   */
  dbl phi_i;
  dbl(*grad_phi_i_e_a)[DIM] = NULL;
  /*
   * Interpolation functions for variables and some of their derivatives.
   */

  dbl Pi[DIM][DIM];

  dbl wt;

  dbl d_area;

  /* coefficient variables for the Brinkman Equation
     in flows through porous media: KSC on 5/10/95 */
  dbl por;       /* porosity of porous media */
  dbl por2;      /* square of porosity */
  dbl per = 0.0; /* permeability of porous media */
  /* derivative of permeability wrt concentration */
  dbl d_per_dc[MAX_CONC][MDE];

  dbl vis; /* flowing-liquid viscosity */
  /* Flowing-liquid viscosity sensitivities */
  VISCOSITY_DEPENDENCE_STRUCT d_flow_mu_struct; /* density dependence */
  VISCOSITY_DEPENDENCE_STRUCT *d_flow_mu = &d_flow_mu_struct;
  dbl sc;    /* inertial coefficient */
  dbl speed; /* magnitude of the velocity vector */

  int v_s[MAX_MODES][DIM][DIM];

  int sign;
  double momentum_residual;

  /* the residual is essentially calculated from
     R = R+ * H[F] + R- * (1-H[F]), so
     dR/dFj = delta[F] * Nj * (R+ - R-)

     So, here we are to assemble for one side of the interface,
     dR/Fj += delta[F] * Nj * sign * R[sign]
   */

  /* Variables used for the modified fluid momentum equations when looking
   * at the particle momentum model.  Note that pd->MomentumFluxModel ==
   * SUSPENSION_PM when the Buyevich particle momentum equations are active.
   */
  int particle_momentum_on; /* Boolean for particle momentum eq.'s active */
  int species = -1;         /* Species number of particle phase */
  double p_vol_frac = 0;    /* Particle volume fraction (phi) */
  double ompvf = 1;         /* 1 - p_vol_frac "One Minus Particle */
  /*    Volume Fraction" */
  int mass_on;
  int advection_on;
  int diffusion_on;
  int source_on;
  int porous_brinkman_on;
  int transient_run = (pd->TimeIntegration != STEADY);

  int *pde = pd->e[pg->imtrx];
  dbl mass_etm;
  dbl advection_etm;
  dbl diffusion_etm;
  dbl porous_brinkman_etm;
  dbl source_etm;

  dbl h_elem_avg;

  /*Continuity stabilization*/
  dbl cont_gls, continuity_stabilization;

  status = 0;

  /*
   * Unpack variables from structures for local convenience...
   */

  eqn = R_MOMENTUM1;
  var = FILL;

  /*
   * Bail out fast if there's nothing to do...
   */

  if (!pd->e[pg->imtrx][eqn] || !pd->v[pg->imtrx][var]) {
    return (status);
  }

  sign = ls->Elem_Sign;

  wt = fv->wt;
  h3 = fv->h3;

  if (ls->on_sharp_surf) /* sharp interface */
  {
    det_J = fv->sdet;
  } else /* diffuse interface */
  {
    det_J = bf[eqn]->detJ;
  }

  d_area = wt * h3 * det_J;

  if (pd->v[pg->imtrx][POLYMER_STRESS11]) {
    (void)stress_eqn_pointer(v_s);
  }

  /* Set up variables for particle/fluid momentum coupling.
   */
  if (pd->e[pg->imtrx][R_PMOMENTUM1]) {
    particle_momentum_on = 1;
    /* This is the species number of the particle phase. */
    species = (int)mp->u_density[0];
    p_vol_frac = fv->c[species];
    ompvf = 1.0 - p_vol_frac;
    /* Uncomment this to check for when the particle volume fraction
     * becomes non-physical.  Beware, however, that the intermediate
     * solutions may, indeed, become negative while converging to a
     * physical solution.
    if(p_vol_frac<0.0 || p_vol_frac>1.0)
      {
        if(fabs(p_vol_frac)<1e-14)
          p_vol_frac=0.0;
        else
          {
            printf("assemble_momentum: p_vol_frac=%g, exiting\n",p_vol_frac);
            exit(0);
          }
      }
    */
  } else
    particle_momentum_on = 0;

  /*
   * Material property constants, etc. Any variations for this
   * Gauss point were evaluated in load_material_properties().
   */
  h_elem_avg = pg_data->h_elem_avg;

  /*** Density ***/

  rho = density(NULL, time);

  if (pd->e[pg->imtrx][eqn] & T_POROUS_BRINK) {
    if (mp->PorousMediaType != POROUS_BRINKMAN)
      GOMA_WH(-1, "Set Porous term multiplier in continuous medium");

    /* Load variable FlowingLiquid_viscosity */
    vis = flowing_liquid_viscosity(d_flow_mu);

    if (mp->PermeabilityModel == SOLIDIFICATION) {
      /* This is a permeability model that slows down
       * the flow with increasing solid fraction.
       * It should be useful for solidification problems.
       */
      per = solidification_permeability(h_elem_avg, d_per_dc);
    } else if (mp->PermeabilityModel == CONSTANT) {
      per = mp->permeability;
    } else {
      GOMA_EH(GOMA_ERROR, "Unrecognizable Permeability model");
    }

    /* Load up remaining parameters for the Brinkman Equation. */
    por = mp->porosity;
    sc = mp->Inertia_coefficient;
  } else {
    por = 1.;
    per = 1.;
    vis = mp->viscosity;
    sc = 0.;
  }

  eqn = R_MOMENTUM1;
  /*
   * Field variables...
   */

  if (transient_run && pd->v[pg->imtrx][MESH_DISPLACEMENT1])
    x_dot = fv_dot->x;
  else
    x_dot = zero;

  if (transient_run)
    v_dot = fv_dot->v;
  else
    v_dot = zero;

  /* for porous media stuff */
  speed = 0.0;
  for (a = 0; a < WIM; a++) {
    speed += v[a] * v[a];
  }
  speed = sqrt(speed);

  for (a = 0; a < VIM; a++)
    grad_v[a] = fv->grad_v[a];

  /*
   * Stress tensor, but don't need dependencies
   */
  fluid_stress(Pi, NULL);

  (void)momentum_source_term(f, NULL, time);

  // Call continuity stabilization if desired
  if (Cont_GLS) {
    calc_cont_gls(&cont_gls, NULL, time, pg_data);
  }

  if (af->Assemble_Jacobian) {
    /*
     * Assemble each component "a" of the momentum equation...
     */
    for (a = 0; a < WIM; a++) {
      eqn = R_MOMENTUM1 + a;
      peqn = upd->ep[pg->imtrx][eqn];
      bfm = bf[eqn];

      mass_on = pde[eqn] & T_MASS;
      advection_on = pde[eqn] & T_ADVECTION;
      diffusion_on = pde[eqn] & T_DIFFUSION;
      source_on = pde[eqn] & T_SOURCE;
      porous_brinkman_on = pde[eqn] & T_POROUS_BRINK;

      mass_etm = pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
      advection_etm = pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
      diffusion_etm = pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
      porous_brinkman_etm = pd->etm[pg->imtrx][eqn][(LOG2_POROUS_BRINK)];
      source_etm = pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
      /*
       * In the element, there will be contributions to this many equations
       * based on the number of degrees of freedom...
       */

      for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
        ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][i];
        if (ei[pg->imtrx]->active_interp_ledof[ledof]) {
          /*
           *  Here is where we figure out whether the row is to placed in
           *  the normal spot (e.g., ii = i), or whether a boundary condition
           *  require that the volumetric contribution be stuck in another
           *  ldof pertaining to the same variable type.
           */
          ii = ei[pg->imtrx]->lvdof_to_row_lvdof[eqn][i];

          phi_i = bfm->phi[i];

          mass = 0.;

          grad_phi_i_e_a = bfm->grad_phi_e[i][a];

          if (transient_run) {
            if (mass_on) {
              mass = v_dot[a] * rho;
              mass *= -phi_i * d_area;
              mass *= mass_etm;
            }

            if (porous_brinkman_on) {
              mass /= por;
            }

            if (particle_momentum_on) {
              mass *= ompvf;
            }
          }

          advection = 0.;
          if (advection_on) {
#ifdef DO_NO_UNROLL
            for (p = 0; p < WIM; p++) {
              advection += (v[p] - x_dot[p]) * grad_v[p][a];
            }

#else
            advection += (v[0] - x_dot[0]) * grad_v[0][a];
            advection += (v[1] - x_dot[1]) * grad_v[1][a];
            if (WIM == 3)
              advection += (v[2] - x_dot[2]) * grad_v[2][a];

#endif
            advection *= rho;
            advection *= -phi_i * d_area;
            advection *= advection_etm;

            if (porous_brinkman_on) {
              por2 = por * por;
              advection /= por2;
            }

            if (particle_momentum_on) {
              advection *= ompvf;
            }
          }

          porous = 0.;
          if (porous_brinkman_on) {
            porous = v[a] * (rho * sc * speed / sqrt(per) + vis / per);
            porous *= -phi_i * d_area;
            porous *= porous_brinkman_etm;
          }

          diffusion = 0.;
          if (diffusion_on) {
#ifdef DO_NO_UNROLL
            for (p = 0; p < VIM; p++) {
              for (int q = 0; q < VIM; q++) {
                diffusion += grad_phi_i_e_a[p][q] * Pi[q][p];
              }
            }
#else
            diffusion += grad_phi_i_e_a[0][0] * Pi[0][0];
            diffusion += grad_phi_i_e_a[1][1] * Pi[1][1];
            diffusion += grad_phi_i_e_a[0][1] * Pi[1][0];
            diffusion += grad_phi_i_e_a[1][0] * Pi[0][1];

            if (VIM == 3) {
              diffusion += grad_phi_i_e_a[2][2] * Pi[2][2];
              diffusion += grad_phi_i_e_a[2][1] * Pi[1][2];
              diffusion += grad_phi_i_e_a[2][0] * Pi[0][2];
              diffusion += grad_phi_i_e_a[1][2] * Pi[2][1];
              diffusion += grad_phi_i_e_a[0][2] * Pi[2][0];
            }

#endif

            diffusion *= -d_area;
            diffusion *= diffusion_etm;
          }

          source = 0.0;
          if (source_on) {
            source += f[a];

            source *= phi_i * d_area;
            source *= source_etm;
          }

          /* MMH For massful particles. */
          if (Particle_Dynamics &&
              (Particle_Model == SWIMMER_EXPLICIT || Particle_Model == SWIMMER_IMPLICIT)) {
            if (a == pd->Num_Dim - 1)
              /* These data values should hold the entire
               * source term. */
              source = element_particle_info[ei[pg->imtrx]->ielem].source_term[i];
          }

          // Continuity residual
          continuity_stabilization = 0.0;
          if (Cont_GLS) {
            for (p = 0; p < VIM; p++) {
              continuity_stabilization += grad_phi_i_e_a[p][p];
            }
            continuity_stabilization *= cont_gls * d_area;
          }

          momentum_residual =
              mass + advection + porous + diffusion + source + continuity_stabilization;

          var = FILL;
          pvar = upd->vp[pg->imtrx][var];
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {

            lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += lsi->d_H_dF[j] * momentum_residual * sign;
          }

        } /*end if (active_dofs) */
      }   /* end of for (i=0,ei[pg->imtrx]->dofs...) */
    }
  }
  return (status);

} /* end of function assemble_momentum_path_dependence                */
/**********************************************************************/

int assemble_continuity_path_dependence(dbl time_value,
                                        dbl tt, /* parameter to vary time integration fromexplicit
                                                   (tt = 1) to implicit (tt = 0)    */
                                        dbl dt, /* current time step size                    */
                                        const PG_DATA *pg_data) {
  int dim;
  int p, a;
  int dofs;

  int eqn, var;
  int peqn, pvar;
  int w;

  int i;
  int j;

  int status, err;

  dbl time = 0.0; /*  RSL 6/6/02  */

  dbl *v = fv->v;        /* Velocity field. */
  dbl div_v = fv->div_v; /* Divergence of v. */

  dbl epsilon = 0, derivative, sum; /*  RSL 7/24/00  */

  dbl advection;
  dbl source;
  dbl pressure_stabilization;

  dbl volsolvent = 0; /* volume fraction of solvent                */
  dbl initial_volsolvent = 0;

  dbl det_J;
  dbl h3;
  dbl wt;

  /*
   * Galerkin weighting functions...
   */

  dbl phi_i;

  dbl(*grad_phi)[DIM]; /* weight-function for PSPG term */

  /*
   * Variables for Pressure Stabilization Petrov-Galerkin...
   */
  int meqn;

  dbl pspg[DIM];

  dbl mass;

  dbl rho = 0;
  DENSITY_DEPENDENCE_STRUCT d_rho_struct; /* density dependence */
  DENSITY_DEPENDENCE_STRUCT *d_rho = &d_rho_struct;

  struct Species_Conservation_Terms s_terms;
  dbl rhos = 0, rhof = 0;
  dbl h_flux = 0;
  int w0 = 0;

  /* For particle momentum model.
   */
  int species;              /* species number for particle phase,  */
  dbl ompvf;                /* 1 - partical volume fraction */
  int particle_momentum_on; /* boolean. */

  /* Foaming model TAB */
  double dFVS_dv[DIM][MDE];
  double dFVS_dT[MDE];
  double dFVS_dx[DIM][MDE];
  double dFVS_dC[MAX_CONC][MDE];
  double dFVS_dF[MDE];

  int sign;
  double continuity;

  status = 0;

  /*
   * Unpack variables from structures for local convenience...
   */

  eqn = R_PRESSURE;
  peqn = upd->ep[pg->imtrx][eqn];

  /*
   * Bail out fast if there's nothing to do...
   */

  if (!pd->e[pg->imtrx][eqn] || !pd->v[pg->imtrx][FILL]) {
    return (status);
  }

  sign = ls->Elem_Sign;

  dim = pd->Num_Dim;

  wt = fv->wt;
  h3 = fv->h3; /* Differential volume element (scales). */

  if (ls->on_sharp_surf) /* sharp interface */
  {
    det_J = fv->sdet;
  } else {
    det_J = bf[eqn]->detJ;
  }

  grad_phi = bf[eqn]->grad_phi;

  /*
   * Get the deformation gradients and tensors if needed
   */

  if ((cr->MeshMotion == LAGRANGIAN || cr->MeshMotion == DYNAMIC_LAGRANGIAN) &&
      pd->e[pg->imtrx][R_MESH1]) {
    err = belly_flop(elc->lame_mu);
    GOMA_EH(err, "error in belly flop");
    if (err == 2)
      return (err);
  }

  if ((cr->MeshMotion == TOTAL_ALE && !pd->v[pg->imtrx][VELOCITY1]) && pd->e[pg->imtrx][R_SOLID1]) {
    err = belly_flop_rs(elc_rs->lame_mu);
    GOMA_EH(err, "error in belly flop for real solid");
    if (err == 2)
      return (err);
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

  if (PSPG) {
    calc_pspg(pspg, NULL, time_value, tt, dt, pg_data);
  }

  if ((cr->MeshMotion == LAGRANGIAN || cr->MeshMotion == DYNAMIC_LAGRANGIAN ||
       (cr->MeshMotion == TOTAL_ALE && !pd->v[pg->imtrx][VELOCITY1])) &&
      (mp->PorousMediaType == CONTINUOUS)) {
    initial_volsolvent = elc->Strss_fr_sol_vol_frac;
    volsolvent = 0.;
    for (w = 0; w < pd->Num_Species_Eqn; w++)
      volsolvent += fv->c[w];
    if (particle_momentum_on)
      volsolvent -= fv->c[species];
  }

  if (mp->MomentumSourceModel == SUSPENSION_PM || mp->SpeciesSourceModel[0] == ELECTRODE_KINETICS) {
    err = get_continuous_species_terms(&s_terms, 0.0, tt, dt, pg_data->hsquared);
    GOMA_EH(err, "problem in getting the species terms");
  }

  if (mp->SpeciesSourceModel[0] == ION_REACTIONS) {
    zero_structure(&s_terms, sizeof(struct Species_Conservation_Terms), 1);
    err = get_continuous_species_terms(&s_terms, time, tt, dt, pg_data->hsquared);
    GOMA_EH(err, "problem in getting the species terms");
  }

  if ((cr->MassFluxModel == HYDRODYNAMIC) && (mp->DensityModel == SUSPENSION) &&
      (mp->MomentumSourceModel == SUSPENSION)) {
    /*
     * Compute hydrodynamic/sedimentation flux and sensitivities.
     */

    w0 = (int)mp->u_density[0]; /* This is the species number that is transported
                                   HYDRODYNAMICally  */

    hydro_flux(&s_terms, w0, tt, dt, pg_data->hsquared);

    rhof = mp->u_density[1];
    rhos = mp->u_density[2];
  }

  if (mp->SpeciesSourceModel[0] == ELECTRODE_KINETICS ||
      mp->SpeciesSourceModel[0] == ION_REACTIONS) {
    rho = density(d_rho, time_value);
  }

  if (af->Assemble_Jacobian) {
    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      phi_i = bf[eqn]->phi[i];

      /* This terms shows up because the fluid phase is not assumed to be
       * incompressible.
       */
      mass = 0.0;

      if (particle_momentum_on) {
        if (pd->TimeIntegration != STEADY) {
          mass = -s_terms.Y_dot[species];
          mass *= phi_i * det_J * h3 * wt;
        }
      }

      if (mp->SpeciesSourceModel[0] == ELECTRODE_KINETICS ||
          mp->SpeciesSourceModel[0] == ION_REACTIONS) {
        mass = 0.0;
        if (pd->TimeIntegration != STEADY) {
          if (mp->PorosityModel == CONSTANT) {
            epsilon = mp->porosity;
          } else if (mp->PorosityModel == THERMAL_BATTERY) {
            epsilon = mp->u_porosity[0];
          } else {
            GOMA_EH(GOMA_ERROR, "invalid porosity model");
          }
          var = MASS_FRACTION;
          for (w = 0; w < pd->Num_Species - 1; w++) {
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              if (bf[var]->phi[j] > 0.0)
                break;
            }
            derivative = d_rho->C[w][j] / bf[var]->phi[j];
            mass += derivative * s_terms.Y_dot[w];
          }
          mass *= epsilon / rho;
          mass *= phi_i * det_J * h3 * wt;
        }
      }

      advection = 0.;

      if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
        if (pd->v[pg->imtrx][VELOCITY1]) /* then must be solving fluid mechanics in this material */
        {

          /*
           * Standard incompressibility constraint means we have
           * a solenoidal velocity field...if we don't, then
           * we might be in serious trouble...
           *
           * HKM -> This is where a variable density implementation must go
           */

          advection = div_v;

          /* We get a more complicated advection term because the
           * particle phase is not incompressible.
           */
          if (particle_momentum_on)
            advection *= ompvf;

          advection *= phi_i * h3 * det_J * wt;
          advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
        } else if (cr->MeshMotion == LAGRANGIAN || cr->MeshMotion == DYNAMIC_LAGRANGIAN ||
                   cr->MeshMotion == TOTAL_ALE)
        /* use divergence of displacement for linear elasticity */
        {
          advection = fv->volume_change;

          if (particle_momentum_on)
            advection *= ompvf;

          advection *= phi_i * h3 * det_J * wt;
          advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
        }
        if (mp->SpeciesSourceModel[0] == ELECTRODE_KINETICS ||
            mp->SpeciesSourceModel[0] == ION_REACTIONS) {
          advection = div_v;
          var = MASS_FRACTION;
          for (w = 0; w < pd->Num_Species - 1; w++) {
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              if (bf[var]->phi[j] > 0.0)
                break;
            }
            derivative = d_rho->C[w][j] / bf[var]->phi[j];
            sum = 0.;
            for (p = 0; p < dim; p++) {
              sum += s_terms.conv_flux[w][p];
            }
            advection += derivative * sum / rho;
          }
          advection *= phi_i * h3 * det_J * wt;
          advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
        }
      }

      source = 0.;

      if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
        /*
         * Maybe you want to penalize the equation to give a
         * nonzero diagonal entry...
         */

        if (pd->v[pg->imtrx][VELOCITY1]) {
          if (mp->DensityModel == CONSTANT || mp->DensityModel == DENSITY_FILL ||
              mp->DensityModel == DENSITY_LEVEL_SET || mp->DensityModel == DENSITY_IDEAL_GAS) {
            /* These density models assume constant local density.  That is,
               the density may be different in different parts of the flow
               but in a sufficiently small region the density behavior is flat
            */
            /* DRN (07/13/05):
               This was previously:
               source     =  P;
               But this messes with level set problems that have a density source
               over part of the domain and constant in other regions.  If someone was
               counting on this behavior as a form of a penalty method to give a non-zero
               diagonal entry, we should implement a new density model that accomplishes
               this. I really don't know what you want for DENSITY_IDEAL_GAS, though?!?
            */
            source = 0.;
          } else if (mp->DensityModel == DENSITY_FOAM || mp->DensityModel == DENSITY_FOAM_CONC ||
                     mp->DensityModel == DENSITY_FOAM_TIME ||
                     mp->DensityModel == DENSITY_FOAM_TIME_TEMP) {
            /*  These density models locally permit a time and spatially varying
             *  density.  Consequently, neither the gradient nor the time derivative
             *  of the density term in the continuity equation is zero. These terms
             *  must then be included here as a "source" term
             */
            source =
                FoamVolumeSource(time_value, dt, tt, dFVS_dv, dFVS_dT, dFVS_dx, dFVS_dC, dFVS_dF);
          } else if (mp->DensityModel == REACTIVE_FOAM) {
            /* These density models locally permit a time and spatially varying
               density.  Consequently, the Lagrangian derivative of the density
               terms in the continuity equation are not zero and are
               included here as a source term
            */
            source = REFVolumeSource(time_value, dt, tt, dFVS_dv, dFVS_dT, dFVS_dx, dFVS_dC);
          } else if (mp->DensityModel == SUSPENSION || mp->DensityModel == SUSPENSION_PM) {
            /* Although the suspension density models meet the definition of
               a locally variable density model, the Lagrangian derivative
               of their densities can be represented as a divergence of
               mass flux.  This term must be integrated by parts and so is
               included separately later on and is not include as part of the "source"
               terms */
            source = 0.0;
          }

          /* To include, or not to include, that is the question
           * when considering the particle momentum coupled eq.'s...
           */

          /* We get this source term because the fluid phase is not
           * incompressible.
           */
          if (particle_momentum_on) {
            source = 0.0;
            for (a = 0; a < WIM; a++) {
              /* Cannot use s_terms.conv_flux[a] here because that
               * is defined in terms of the particle phase
               * velocities.
               */
              source -= fv->grad_c[species][a] * v[a];
            }
          }
          source *= phi_i * h3 * det_J * wt;
          source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
        }

        if ((cr->MeshMotion == LAGRANGIAN || cr->MeshMotion == DYNAMIC_LAGRANGIAN ||
             (cr->MeshMotion == TOTAL_ALE && !pd->v[pg->imtrx][VELOCITY1])))
        /* add swelling as a source of volume */
        {
          if (mp->PorousMediaType == CONTINUOUS) {
            source = -(1. - initial_volsolvent) / (1. - volsolvent);
            source *= phi_i * h3 * det_J * wt;
            source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
          }
        }
        if (mp->SpeciesSourceModel[0] == ELECTRODE_KINETICS ||
            mp->SpeciesSourceModel[0] == ION_REACTIONS) {
          source = 0.;
          for (j = 0; j < pd->Num_Species; j++) {
            source -= s_terms.MassSource[j] * mp->molecular_weight[j];
          }
          source /= rho;
          source *= phi_i * h3 * det_J * wt;
          source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
        }
      }

      /* add Pressure-Stabilized Petrov-Galerkin term
       * if desired.
       */

      pressure_stabilization = 0.;

      if (PSPG) {
        for (a = 0; a < WIM; a++) {
          meqn = R_MOMENTUM1 + a;
          if (pd->e[pg->imtrx][meqn]) {
            pressure_stabilization += grad_phi[i][a] * pspg[a];
          }
        }
        pressure_stabilization *= h3 * det_J * wt;
      }

      h_flux = 0.0;

      if ((cr->MassFluxModel == HYDRODYNAMIC) && (mp->MomentumSourceModel == SUSPENSION)) {
        /* add divergence of particle phase flux as source term */

        /* The particle flux terms has been integrated by parts.
         * No boundary integrals are included in this formulation
         * so it is tacitly assumed that the particle phase
         * relative mass flux over all boundaries is zero
         */

        for (p = 0; p < dim; p++) {
          h_flux += grad_phi[i][p] * s_terms.diff_flux[w0][p];
        }

        h_flux *= (rhos - rhof) / rhof;
        h_flux *= h3 * det_J * wt;
        h_flux *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
        /*  h_flux = 0.0; */
      }

      continuity = mass + advection + source + pressure_stabilization + h_flux;

      var = FILL;
      pvar = upd->vp[pg->imtrx][var];
      dofs = ei[pg->imtrx]->dof[var];
      j = 0;
      while (j < dofs) {

        lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += lsi->d_H_dF[j] * continuity * sign;

        j++;
      }
    }
  }

  return (status);

} /* end of function assemble_continuity_path_dependence              */
/**********************************************************************/

int assemble_LM_source(double *xi,
                       int oAC,     /* Flag indicating calling function */
                       double *gAC, /* Augmenting condition arrays */
                       double **bAC,
                       double **cAC,
                       double x[],
                       Exo_DB *exo) {
  int i, j, a, ii, ledof, v;
  int eqn, peqn;
  int pass;
  int ac_lm = Do_Overlap;
  int id_side, nu, iAC = 0, ioffset = 0;
  int dof_l[DIM];
  double phi_l[DIM][MDE];
  int ielem = ei[pg->imtrx]->ielem;
  int iconn_ptr = exo->elem_ptr[ielem];
  double lagrange_mult[3] = {0.0, 0.0, 0.0};

  struct Basis_Functions *bfm;

  double wt, det_J, h3, phi_i;

  double source;

  load_lsi(ls->Length_Scale);
  if (lsi->delta == 0.)
    return (0);

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

  pass = ((Do_Overlap && oAC >= 0) ? 2 : 1);

  ielem = ei[pg->imtrx]->ielem;
  if ((ls->Evolution == LS_EVOLVE_SLAVE) && (ls->init_surf_list->start->type == LS_SURF_SS)) {
    struct LS_Surf *ss_surf;
    struct LS_Surf_Closest_Point *cp;

    ss_surf = closest_surf(ls->init_surf_list, x, exo, fv->x);
    cp = ss_surf->closest_point;

    /* use kitchen sink approach for now at solids location */
    if (cp->elem == -1)
      GOMA_EH(GOMA_ERROR, "Invalid element at closest_point");

    /* Associate the correct AC number with this point */
    id_side = cp->elem_side;
    ioffset = first_overlap_ac(cp->elem, id_side);
    if (ioffset == -1)
      GOMA_EH(GOMA_ERROR, "Bad AC index");

    setup_shop_at_point(cp->elem, cp->xi, exo);

    for (a = 0; a < ei[pg->imtrx]->ielem_dim; a++) {
      if (ac_lm) {
        lagrange_mult[a] = augc[ioffset + a].lm_value;
      } else {
        lagrange_mult[a] = fv->lm[a];
      }

      /* Grab stuff for cross term evaluation */
      if (pass == 2) {
        v = LAGR_MULT1 + a;
        if (pd->v[pg->imtrx][v] && !ac_lm) {
          dof_l[a] = ei[pg->imtrx]->dof[v];
          for (j = 0; j < dof_l[a]; j++) {
            phi_l[a][j] = bf[v]->phi[j];
          }
        } else {
          dof_l[a] = 1;
          phi_l[a][0] = 1.0;
        }
      }
    }

    /* Yep, that's all I needed from closest_point.
     * Now go back to the original element
     */
    setup_shop_at_point(ielem, xi, exo);
  }

  /*
   * Wesiduals
   * ________________________________________________________________________________
   */
  for (a = 0; a < WIM; a++) {
    eqn = R_MOMENTUM1 + a;
    peqn = upd->ep[pg->imtrx][eqn];
    bfm = bf[eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {

      ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][i];

      if (ei[pg->imtrx]->active_interp_ledof[ledof]) {

        ii = ei[pg->imtrx]->lvdof_to_row_lvdof[eqn][i];

        phi_i = bfm->phi[i];

        /* Fluid momentum residuals */
        /*lec->R[LEC_R_INDEX(peqn,i)] += wt *(lagrange_mult[a]) *
          det_J * 0.5 * length; */
        source = -phi_i * (lagrange_mult[a]) * lsi->delta;

        source *= det_J * wt;
        source *= h3;
        source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

        if (pass == 1)
          lec->R[LEC_R_INDEX(peqn, ii)] += source;

        /* Sensitivities to solid Lagrange multiplier
           (cross-mesh term) */
        if (Do_Overlap && pass == 2) {
          v = LAGR_MULT1 + a;
          iAC = ioffset + a;
          nu = lookup_active_dof(eqn, i, iconn_ptr);

          for (j = 0; j < dof_l[a]; j++) {
            if (nu >= 0) {
              source = -phi_i * phi_l[a][j] * lsi->delta;

              source *= det_J * wt;
              source *= h3;
              source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

              bAC[iAC][nu] += source;
            }
          }
        }
      }
    }
  }

  return (1);
}

int assemble_max_strain(void) {
  /*****************************************************************************
   * assemble_max_strain ()
   *
   * Author:  Scott A Roberts, sarober@sandia.gov, 1514
   *
   * Date:    April 12, 2012
   *
   * Purpose: Calculates the maximum von Mises strain that the material has
   *          experienced over the length of the simulation.  This was primarily
   *          implemented for the FAUX_PLASTICITY model for the modulus in
   *          solid mechanics.
   *****************************************************************************/

  // Define variables
  int eqn, peqn, var, pvar;
  int a, b, i, j, k;
  dbl phi_i, phi_j;
  dbl mass, source;
  dbl vmE, d_vmE_dstrain[DIM][DIM];

  // Initialize output status function
  int status = 0;

  // Bail out of function if there is nothing to do
  eqn = R_MAX_STRAIN;
  if (!pd->e[pg->imtrx][eqn])
    return (status);

  // Unpack FEM variables from global structures
  dbl wt = fv->wt;           // Gauss point weight
  dbl h3 = fv->h3;           // Differential volume element
  dbl det_J = bf[eqn]->detJ; // Jacobian of transformation

  // FEM weights
  dbl dA = det_J * wt * h3;
  dbl etm_mass = pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
  dbl etm_source = pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

  // Calculate strain
  /*
   * N.B.:  This strain is calculated in the same way the von Mises
   * stress is calculated, sqrt(3*II(T_dev)).  However, for this application,
   * we want to use a slightly different version, sqrt(4*II(E_dev)/3).
   * Therefore, we modify the returned values afterwards.
   */
  vmE = calc_tensor_invariant(fv->strain, d_vmE_dstrain, 4);
  vmE *= 2.0 / 3.0;
  for (a = 0; a < DIM; a++) {
    for (b = 0; b < DIM; b++) {
      d_vmE_dstrain[a][b] *= 2.0 / 3.0;
    }
  }

  // Prepare tests
  int use_new = vmE >= fv_old->max_strain;
  int mass_lump = 0;

  /* --- Assemble residuals --------------------------------------------------*/
  eqn = R_MAX_STRAIN;
  if (af->Assemble_Residual) {
    peqn = upd->ep[pg->imtrx][eqn];

    // Loop over DOF (i)
    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {

      // Load basis functions
      phi_i = bf[eqn]->phi[i];

      // Assemble mass term
      mass = 0.0;
      if (pd->e[pg->imtrx][eqn] & T_MASS) {
        if (mass_lump) {
          mass -= *esp->max_strain[i] * phi_i;
        } else {
          mass -= fv->max_strain * phi_i;
        }
      }
      mass *= dA * etm_mass;

      // Assemble source term
      source = 0.0;
      if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
        if (use_new) {
          source += vmE * phi_i;
        } else {
          source += fv_old->max_strain * phi_i;
        }
      }
      source *= dA * etm_source;

      // Assemble full residual
      lec->R[LEC_R_INDEX(peqn, i)] += mass + source;

    } // End of loop over DOF (i)

  } // End of residual assembly of R_MAX_STRAIN

  /* --- Assemble Jacobian --------------------------------------------------*/
  eqn = R_MAX_STRAIN;
  if (af->Assemble_Jacobian) {
    peqn = upd->ep[pg->imtrx][eqn];

    // Loop over DOF (i)
    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {

      // Load basis functions
      phi_i = bf[eqn]->phi[i];

      // Assemble sensitivities for MAX_STRAIN
      var = MAX_STRAIN;
      if (pd->v[pg->imtrx][var]) {
        pvar = upd->vp[pg->imtrx][var];

        // Loop over DOF (j)
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {

          // Load basis functions
          phi_j = bf[eqn]->phi[j];

          // Assemble mass terms
          mass = 0.0;
          if (pd->e[pg->imtrx][eqn] & T_MASS) {
            if (mass_lump) {
              mass -= phi_i * delta(i, j);
            } else {
              mass -= phi_i * phi_j;
            }
          }
          mass *= dA * etm_mass;

          // Assemble source term
          source = 0.0;

          // Assemble full Jacobian
          lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += mass + source;

        } // End of loop over DOF (j)

      } // End of MAX_STRAIN sensitivities

      // Assemble sensitivities for MESH_DISPLACEMENT
      for (k = 0; k < DIM; k++) {
        var = MESH_DISPLACEMENT1 + k;
        if (pd->v[pg->imtrx][var]) {
          pvar = upd->vp[pg->imtrx][var];

          // Loop over DOF (j)
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {

            // Load basis functions
            phi_j = bf[eqn]->phi[j];

            // Assemble mass terms
            mass = 0.0;
            if (pd->e[pg->imtrx][eqn] & T_MASS) {
              if (mass_lump) {
                mass -= *esp->max_strain[i] * phi_i * fv->dh3dmesh[k][j] * det_J;
                mass -= *esp->max_strain[i] * phi_i * h3 * bf[eqn]->d_det_J_dm[k][j];
              } else {
                mass -= fv->max_strain * phi_i * fv->dh3dmesh[k][j] * det_J;
                mass -= fv->max_strain * phi_i * h3 * bf[eqn]->d_det_J_dm[k][j];
              }
            }
            mass *= wt * etm_mass;

            // Assemble source term
            source = 0.0;
            if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
              if (use_new) {
                for (a = 0; a < DIM; a++) {
                  for (b = 0; b < DIM; b++) {
                    source += d_vmE_dstrain[a][b] * fv->d_strain_dx[a][b][k][j] * h3 * det_J;
                  }
                }
                source += vmE * fv->dh3dmesh[k][j] * det_J;
                source += vmE * h3 * bf[eqn]->d_det_J_dm[k][j];
              } else {
                source += fv_old->max_strain * fv->dh3dmesh[k][j] * det_J;
                source += fv_old->max_strain * h3 * bf[eqn]->d_det_J_dm[k][j];
              }
            }
            source *= phi_i * wt * etm_source;

            // Assemble full Jacobian
            lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += mass + source;

          } // End of loop over DOF (j)

        } // End of MESH_DISPLACEMENTi sensitivities

      } // End loop over MESH_DISPLACEMENT vector components

    } // End of loop over DOF (i)

  } // End of Jacobian assembly

  return (status);
} // End of assemble_max_strain()

int assemble_cur_strain(void) {
  /*****************************************************************************
   * assemble_cur_strain ()
   *
   * Author:  Scott A Roberts, sarober@sandia.gov, 1514
   *
   * Date:    April 23, 2012
   *
   * Purpose: Calculates the current von Mises strain.  This was primarily
   *          implemented for the FAUX_PLASTICITY model for the modulus in
   *          solid mechanics.
   *****************************************************************************/

  // Define variables
  int eqn, peqn, var, pvar;
  int a, b, i, j, k;
  dbl phi_i, phi_j;
  dbl mass, source;
  dbl vmE, d_vmE_dstrain[DIM][DIM];

  // Initialize output status function
  int status = 0;

  // Bail out of function if there is nothing to do
  eqn = R_CUR_STRAIN;
  if (!pd->e[pg->imtrx][eqn])
    return (status);

  // Unpack FEM variables from global structures
  dbl wt = fv->wt;           // Gauss point weight
  dbl h3 = fv->h3;           // Differential volume element
  dbl det_J = bf[eqn]->detJ; // Jacobian of transformation

  // FEM weights
  dbl dA = det_J * wt * h3;
  dbl etm_mass = pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
  dbl etm_source = pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

  // Calculate strain
  /*
   * N.B.:  This strain is calculated in the same way the von Mises
   * stress is calculated, sqrt(3*II(T_dev)).  However, for this application,
   * we want to use a slightly different version, sqrt(4*II(E_dev)/3).
   * Therefore, we modify the returned values afterwards.
   */
  vmE = calc_tensor_invariant(fv->strain, d_vmE_dstrain, 4);
  vmE *= 2.0 / 3.0;
  for (a = 0; a < DIM; a++) {
    for (b = 0; b < DIM; b++) {
      d_vmE_dstrain[a][b] *= 2.0 / 3.0;
    }
  }

  // Prepare tests
  int mass_lump = 0;

  /* --- Assemble residuals --------------------------------------------------*/
  eqn = R_CUR_STRAIN;
  if (af->Assemble_Residual) {
    peqn = upd->ep[pg->imtrx][eqn];

    // Loop over DOF (i)
    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {

      // Load basis functions
      phi_i = bf[eqn]->phi[i];

      // Assemble mass term
      mass = 0.0;
      if (pd->e[pg->imtrx][eqn] & T_MASS) {
        if (mass_lump) {
          mass -= *esp->cur_strain[i] * phi_i;
        } else {
          mass -= fv->cur_strain * phi_i;
        }
      }
      mass *= dA * etm_mass;

      // Assemble source term
      source = 0.0;
      if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
        source += vmE * phi_i;
      }
      source *= dA * etm_source;

      // Assemble full residual
      lec->R[LEC_R_INDEX(peqn, i)] += mass + source;

    } // End of loop over DOF (i)

  } // End of residual assembly of R_CUR_STRAIN

  /* --- Assemble Jacobian --------------------------------------------------*/
  eqn = R_CUR_STRAIN;
  if (af->Assemble_Jacobian) {
    peqn = upd->ep[pg->imtrx][eqn];

    // Loop over DOF (i)
    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {

      // Load basis functions
      phi_i = bf[eqn]->phi[i];

      // Assemble sensitivities for CUR_STRAIN
      var = CUR_STRAIN;
      if (pd->v[pg->imtrx][var]) {
        pvar = upd->vp[pg->imtrx][var];

        // Loop over DOF (j)
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {

          // Load basis functions
          phi_j = bf[eqn]->phi[j];

          // Assemble mass terms
          mass = 0.0;
          if (pd->e[pg->imtrx][eqn] & T_MASS) {
            if (mass_lump) {
              mass -= phi_i * delta(i, j);
            } else {
              mass -= phi_i * phi_j;
            }
          }
          mass *= dA * etm_mass;

          // Assemble source term
          source = 0.0;

          // Assemble full Jacobian
          lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += mass + source;

        } // End of loop over DOF (j)

      } // End of CUR_STRAIN sensitivities

      // Assemble sensitivities for MESH_DISPLACEMENT
      for (k = 0; k < DIM; k++) {
        var = MESH_DISPLACEMENT1 + k;
        if (pd->v[pg->imtrx][var]) {
          pvar = upd->vp[pg->imtrx][var];

          // Loop over DOF (j)
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {

            // Load basis functions
            phi_j = bf[eqn]->phi[j];

            // Assemble mass terms
            mass = 0.0;
            if (pd->e[pg->imtrx][eqn] & T_MASS) {
              if (mass_lump) {
                mass -= *esp->cur_strain[i] * phi_i * fv->dh3dmesh[k][j] * det_J;
                mass -= *esp->cur_strain[i] * phi_i * h3 * bf[eqn]->d_det_J_dm[k][j];
              } else {
                mass -= fv->cur_strain * phi_i * fv->dh3dmesh[k][j] * det_J;
                mass -= fv->cur_strain * phi_i * h3 * bf[eqn]->d_det_J_dm[k][j];
              }
            }
            mass *= wt * etm_mass;

            // Assemble source term
            source = 0.0;
            if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
              for (a = 0; a < DIM; a++) {
                for (b = 0; b < DIM; b++) {
                  source += d_vmE_dstrain[a][b] * fv->d_strain_dx[a][b][k][j] * h3 * det_J;
                }
              }
              source += vmE * fv->dh3dmesh[k][j] * det_J;
              source += vmE * h3 * bf[eqn]->d_det_J_dm[k][j];
            }
            source *= phi_i * wt * etm_source;

            // Assemble full Jacobian
            lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += mass + source;

          } // End of loop over DOF (j)

        } // End of MESH_DISPLACEMENTi sensitivities

      } // End loop over MESH_DISPLACEMENT vector components

    } // End of loop over DOF (i)

  } // End of Jacobian assembly

  return (status);
} // End of assemble_cur_strain()

/*  _______________________________________________________________________  */

/* assemble_acoustic -- assemble terms (Residual &| Jacobian) for acoustic harmonic
 *				wave equations
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
 * Created:	Thu June 28, 2005 - Robert Secor
 *
 * Note: currently we do a "double load" into the addresses in the global
 *       "a" matrix and resid vector held in "esp", as well as into the
 *       local accumulators in "lec".
 *
 */

int assemble_acoustic(double time, /* present time value */
                      double tt,   /* parameter to vary time integration from
                                    * explicit (tt = 1) to implicit (tt = 0) */
                      double dt,   /* current time step size */
                      const PG_DATA *pg_data,
                      const int ac_eqn, /* acoustic eqn id and var id	*/
                      const int ac_var) {
  int eqn, var, peqn, pvar, dim, p, b, w, i, j, status;
  int conj_var = 0; /* identity of conjugate variable  */

  dbl P = 0, P_conj = 0, sign_conj = 0; /* acoustic pressure	*/
  dbl q[DIM];
  ACOUSTIC_FLUX_DEPENDENCE_STRUCT d_q_struct;
  ACOUSTIC_FLUX_DEPENDENCE_STRUCT *d_q = &d_q_struct;

  dbl R; /* Acoustic impedance. */
  CONDUCTIVITY_DEPENDENCE_STRUCT d_R_struct;
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_R = &d_R_struct;

  dbl k; /* Acoustic wave number. */
  CONDUCTIVITY_DEPENDENCE_STRUCT d_k_struct;
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_k = &d_k_struct;
  dbl ksqr_sign; /* Sign of wavenumber squared */

  dbl alpha; /* Acoustic Absorption */
  CONDUCTIVITY_DEPENDENCE_STRUCT d_alpha_struct;
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_alpha = &d_alpha_struct;

  dbl mass; /* For terms and their derivatives */

  dbl advection; /* For terms and their derivatives */

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
   * Interpolation functions for variables and some of their derivatives.
   */

  dbl phi_j;

  dbl h3;          /* Volume element (scale factors). */
  dbl dh3dmesh_bj; /* Sensitivity to (b,j) mesh dof. */

  dbl det_J;

  dbl d_det_J_dmeshbj;        /* for specified (b,j) mesh dof */
  dbl dgrad_phi_i_dmesh[DIM]; /* ditto.  */
  dbl wt;

  /* initialize grad_phi_i */
  for (i = 0; i < DIM; i++) {
    grad_phi_i[i] = 0;
  }

  /*   static char yo[] = "assemble_acoustic";*/
  status = 0;
  /*
   * Unpack variables from structures for local convenience...
   */
  dim = pd->Num_Dim;
  eqn = ac_eqn;
  /*
   * Bail out fast if there's nothing to do...
   */
  if (!pd->e[pg->imtrx][eqn]) {
    return (status);
  }

  wt = fv->wt;           /* Gauss point weight. */
  h3 = fv->h3;           /* Differential volume element. */
  det_J = bf[eqn]->detJ; /* Really, ought to be mesh eqn. */

  /*
   * Material property constants, etc. Any variations at this Gauss point
   * are calculated with user-defined subroutine options.
   *
   * Properties to consider here are R and k. R and k we will allow to vary
   * with temperature, spatial coordinates, and species concentration.
   */

  R = acoustic_impedance(d_R, time);

  k = wave_number(d_k, time);

  alpha = acoustic_absorption(d_alpha, time);

  ksqr_sign = mp->acoustic_ksquared_sign;

  acoustic_flux(q, d_q, time, ac_eqn, ac_var);

  if (ac_eqn == R_ACOUS_PREAL) {
    P = fv->apr;
    P_conj = fv->api;
    conj_var = ACOUS_PIMAG;
    sign_conj = 1.;
  } else if (eqn == R_ACOUS_PIMAG) {
    P = fv->api;
    P_conj = -fv->apr;
    conj_var = ACOUS_PREAL;
    sign_conj = -1.;
  } else {
    GOMA_EH(GOMA_ERROR, "Invalid Acoustic eqn");
  }
  /*
   * Residuals___________________________________________________________
   */

  if (af->Assemble_Residual) {
    eqn = ac_eqn;
    peqn = upd->ep[pg->imtrx][eqn];
    var = ac_var;
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

      advection = 0.;
      if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
        advection += phi_i * (k / R) * 2 * alpha * P_conj * det_J * wt;
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

      source = 0.;
      if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
        source -= phi_i * (ksqr_sign * k / R) * P * det_J * wt;
        source *= h3;
        source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
      }

      lec->R[LEC_R_INDEX(peqn, i)] += mass + advection + diffusion + source;
    }
  }

  /*
   * Jacobian terms...
   */

  if (af->Assemble_Jacobian) {
    eqn = ac_eqn;
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

      /*
       * Set up some preliminaries that are needed for the (a,i)
       * equation for bunches of (b,j) column variables...
       */

      for (p = 0; p < VIM; p++) {
        grad_phi_i[p] = bf[eqn]->grad_phi[i][p];
      }

      /*
       * J_e_ap
       */
      var = ac_var;
      if (pd->v[pg->imtrx][var]) {
        pvar = upd->vp[pg->imtrx][var];
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          phi_j = bf[var]->phi[j];

          mass = 0.;

          advection = 0.;

          diffusion = 0.;
          if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {
            for (p = 0; p < VIM; p++) {
              diffusion += d_q->P[p][j] * grad_phi_i[p];
            }
            diffusion *= det_J * wt;
            diffusion *= h3;
            diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
          }

          source = 0.;
          if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
            source -= phi_i * (ksqr_sign * k / R) * phi_j * det_J * wt;
            source *= h3;
            source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
          }

          lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += mass + advection + diffusion + source;
        }
      }
      /*
       *  Conjugate pressure variable
       */
      var = conj_var;
      if (pd->v[pg->imtrx][var]) {
        pvar = upd->vp[pg->imtrx][var];
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          phi_j = bf[var]->phi[j];

          advection = 0.;
          if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
            advection += phi_i * (k / R) * 2 * alpha * sign_conj * phi_j * det_J * wt;
            advection *= h3;
            advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
          }
          lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += advection;
        }
      }
      /*
       * J_e_T
       */
      var = TEMPERATURE;
      if (pd->v[pg->imtrx][var]) {
        pvar = upd->vp[pg->imtrx][var];
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          phi_j = bf[var]->phi[j];

          mass = 0.;

          advection = 0.;
          if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
            advection += phi_i * 2 * P_conj *
                         (R * (alpha * d_k->T[j] + k * d_alpha->T[j]) - k * alpha * d_R->T[j]) /
                         (R * R) * det_J * wt;
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

          source = 0.;
          if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
            source -=
                phi_i * ksqr_sign * (R * d_k->T[j] - k * d_R->T[j]) / (R * R) * P * det_J * wt;
            source *= h3;
            source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
          }

          lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += mass + advection + diffusion + source;
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

            mass = 0.;

            advection = 0.;
            if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
              advection +=
                  phi_i * 2 * P_conj *
                  ((k * alpha / R) * (d_det_J_dmeshbj * h3 + det_J * dh3dmesh_bj) +
                   (R * (alpha * d_k->X[b][j] + k * d_alpha->X[b][j]) - k * alpha * d_R->X[b][j]) /
                       (R * R) * det_J * h3) *
                  wt;
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
              source -= phi_i * P *
                        ((ksqr_sign * k / R) * (d_det_J_dmeshbj * h3 + det_J * dh3dmesh_bj) +
                         ksqr_sign * (R * d_k->X[b][j] - k * d_R->X[b][j]) / (R * R) * det_J * h3) *
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

            advection = 0.;
            if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
              advection +=
                  phi_i * 2 * P_conj *
                  (R * (alpha * d_k->C[w][j] + k * d_alpha->C[w][j]) - k * alpha * d_R->C[w][j]) /
                  (R * R) * det_J * wt;
              advection *= h3;
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
              source -= phi_i * ksqr_sign * (R * d_k->C[w][j] - k * d_R->C[w][j]) / (R * R) * P *
                        det_J * wt;
              source *= h3;
              source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
            }

            lec->J[LEC_J_INDEX(peqn, MAX_PROB_VAR + w, i, j)] +=
                advection + mass + diffusion + source;
          }
        }
      }
    }
  }

  return (status);
} /* end of assemble_acoustic */

/* assemble_acoustic_reynolds_stress -- assemble terms for acoustic reynolds stress
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
 * Created:	Thu Aug. 18, 2005 - Robert Secor
 *
 * Note: currently we do a "double load" into the addresses in the global
 *       "a" matrix and resid vector held in "esp", as well as into the
 *       local accumulators in "lec".
 *
 */

int assemble_acoustic_reynolds_stress(double time, /* present time value */
                                      double tt,   /* parameter to vary time integration from
                                                    * explicit (tt = 1) to implicit (tt = 0) */
                                      double dt,   /* current time step size */
                                      const PG_DATA *pg_data) {
  int eqn, var, peqn, pvar, dim, p, b, w, i, j, status;

  dbl R; /* Acoustic impedance. */
  CONDUCTIVITY_DEPENDENCE_STRUCT d_R_struct;
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_R = &d_R_struct;

  dbl k; /* Acoustic wave number. */
  CONDUCTIVITY_DEPENDENCE_STRUCT d_k_struct;
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_k = &d_k_struct;

  dbl omega, acous_pgrad = 0;

  dbl mass; /* For terms and their derivatives */

  dbl advection; /* For terms and their derivatives */

  dbl source;

  /*
   * Galerkin weighting functions for i-th energy residuals
   * and some of their derivatives...
   */

  dbl phi_i;

  /*
   * Interpolation functions for variables and some of their derivatives.
   */

  dbl phi_j;
  dbl grad_phi_j[DIM];

  dbl h3;          /* Volume element (scale factors). */
  dbl dh3dmesh_bj; /* Sensitivity to (b,j) mesh dof. */

  dbl det_J;

  dbl d_det_J_dmeshbj; /* for specified (b,j) mesh dof */
  dbl wt;

  /*   static char yo[] = "assemble_acoustic";*/
  status = 0;
  /*
   * Unpack variables from structures for local convenience...
   */
  dim = pd->Num_Dim;
  eqn = R_ACOUS_REYN_STRESS;
  /*
   * Bail out fast if there's nothing to do...
   */
  if (!pd->e[pg->imtrx][eqn]) {
    return (status);
  }

  wt = fv->wt;           /* Gauss point weight. */
  h3 = fv->h3;           /* Differential volume element. */
  det_J = bf[eqn]->detJ; /* Really, ought to be mesh eqn. */

  /*
   * Material property constants, etc. Any variations at this Gauss point
   * are calculated with user-defined subroutine options.
   *
   * Properties to consider here are R and k. R and k we will allow to vary
   * with temperature, spatial coordinates, and species concentration.
   */

  R = acoustic_impedance(d_R, time);

  k = wave_number(d_k, time);
  omega = upd->Acoustic_Frequency;

  for (p = 0; p < dim; p++) {
    acous_pgrad += fv->grad_api[p] * fv->grad_api[p];
    acous_pgrad += fv->grad_apr[p] * fv->grad_apr[p];
  }

  /*
   * Residuals___________________________________________________________
   */

  if (af->Assemble_Residual) {
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

      mass = 0.;
      if (pd->e[pg->imtrx][eqn] & T_MASS) {
        mass += phi_i * fv->ars * det_J * wt;
        mass *= h3;
        mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
      }

      advection = 0.;
      if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
        advection -= phi_i * (1. / (4. * k * R * omega)) * (acous_pgrad)*det_J * wt;
        advection *= h3;
        advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
      }

      source = 0.;
      if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
        source -=
            phi_i * k / (4. * R * omega) * (fv->apr * fv->apr + fv->api * fv->api) * det_J * wt;
        source *= h3;
        source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
      }

      lec->R[LEC_R_INDEX(peqn, i)] += mass + advection + source;
    }
  }

  /*
   * Jacobian terms...
   */

  if (af->Assemble_Jacobian) {
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

      /*
       * Set up some preliminaries that are needed for the (a,i)
       * equation for bunches of (b,j) column variables...
       */

      /*
       * J_e_ars
       */
      var = ACOUS_REYN_STRESS;
      if (pd->v[pg->imtrx][var]) {
        pvar = upd->vp[pg->imtrx][var];
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          phi_j = bf[var]->phi[j];

          mass = 0.;
          if (pd->e[pg->imtrx][eqn] & T_MASS) {
            mass += phi_i * phi_j * det_J * wt;
            mass *= h3;
            mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
          }

          advection = 0.;

          source = 0.;

          lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += mass + advection + source;
        }
      }
      /*
       * J_e_ap
       */
      var = ACOUS_PREAL;
      if (pd->v[pg->imtrx][var]) {
        pvar = upd->vp[pg->imtrx][var];
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          phi_j = bf[var]->phi[j];

          for (p = 0; p < VIM; p++) {
            grad_phi_j[p] = bf[var]->grad_phi[j][p];
          }

          mass = 0.;

          advection = 0.;
          if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
            for (p = 0; p < VIM; p++) {
              advection += 2. * fv->grad_apr[p] * grad_phi_j[p];
            }
            advection *= -phi_i * (1. / (4. * k * R * omega)) * det_J * wt;
            advection *= h3;
            advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
          }

          source = 0.;
          if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
            source -= phi_i * k / (4. * R * omega) * (2. * fv->apr * phi_j) * det_J * wt;
            source *= h3;
            source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
          }

          lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += mass + advection + source;
        }
      }
      /*
       *  Conjugate pressure variable
       */
      var = ACOUS_PIMAG;
      if (pd->v[pg->imtrx][var]) {
        pvar = upd->vp[pg->imtrx][var];
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          phi_j = bf[var]->phi[j];

          for (p = 0; p < VIM; p++) {
            grad_phi_j[p] = bf[var]->grad_phi[j][p];
          }

          advection = 0.;
          if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
            for (p = 0; p < VIM; p++) {
              advection += 2. * fv->grad_api[p] * grad_phi_j[p];
            }
            advection *= -phi_i * (1. / (4. * k * R * omega)) * det_J * wt;
            advection *= h3;
            advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
          }

          source = 0.;
          if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
            source -= phi_i * k / (4. * R * omega) * (2. * fv->api * phi_j) * det_J * wt;
            source *= h3;
            source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
          }

          lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += advection + source;
        }
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

          advection = 0.;
          if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
            advection -= phi_i / (4. * omega) * (-1. / (k * k * R * R)) *
                         (k * d_R->T[j] + R * d_k->T[j]) * (acous_pgrad)*det_J * wt;
            advection *= h3;
            advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
          }

          source = 0.;
          if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
            source -= phi_i / (4. * omega) * (R * d_k->T[j] - k * d_R->T[j]) / (R * R) *
                      (fv->apr * fv->apr + fv->api * fv->api) * det_J * wt;
            source *= h3;
            source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
          }

          lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += mass + advection + source;
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

            mass = 0.;
            if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
              mass += phi_i * fv->ars * (d_det_J_dmeshbj * h3 + det_J * dh3dmesh_bj) * wt;
              mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
            }
            advection = 0.;
            if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
              for (p = 0; p < VIM; p++) {
                advection += 2. * (fv->grad_apr[p] * fv->d_grad_apr_dmesh[p][b][j] +
                                   fv->grad_api[p] * fv->d_grad_api_dmesh[p][b][j]);
              }
              advection *= -phi_i * (1. / (4. * k * R * omega)) * det_J * h3 * wt;
              advection -= phi_i * (1. / (4. * k * R * omega)) * (acous_pgrad) *
                           (d_det_J_dmeshbj * h3 + det_J * dh3dmesh_bj) * wt;
              advection -= phi_i / (4. * omega) * (-1. / (k * R * k * R)) *
                           (k * d_R->X[b][j] + R * d_k->X[b][j]) * (acous_pgrad)*det_J * h3 * wt;
              advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
            }

            source = 0.;

            if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
              source -= phi_i / (4. * omega) *
                        ((k / R) * (fv->apr * fv->apr + fv->api * fv->api) *
                             (d_det_J_dmeshbj * h3 + det_J * dh3dmesh_bj) +
                         det_J * h3 * (R * d_k->X[b][j] - k * d_R->X[b][j]) / (R * R)) *
                        wt;
              source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
            }

            lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += mass + advection + source;
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

            advection = 0.;
            if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
              advection -= phi_i / (4. * omega) * (-1. / (k * k * R * R)) *
                           (k * d_R->C[w][j] + R * d_k->C[w][j]) * acous_pgrad * det_J * wt;
              advection *= h3;
              advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
            }

            source = 0.;
            if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
              source -= phi_i / (4. * omega) * (R * d_k->C[w][j] - k * d_R->C[w][j]) / (R * R) *
                        (fv->apr * fv->apr + fv->api * fv->api) * det_J * wt;
              source *= h3;
              source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
            }

            lec->J[LEC_J_INDEX(peqn, MAX_PROB_VAR + w, i, j)] += advection + mass + source;
          }
        }
      }
    }
  }

  return (status);
} /* end of assemble_acoustic_reynolds_stress */

void acoustic_flux(double q[DIM],
                   ACOUSTIC_FLUX_DEPENDENCE_STRUCT *d_q,
                   double time,
                   const int eqn,
                   const int ac_var) {
  dbl R;                                     /* Acoustic impedance. */
  CONDUCTIVITY_DEPENDENCE_STRUCT d_R_struct; /* Acoustic impedance dependence. */
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_R = &d_R_struct;

  dbl k;                                     /* Acoustic wave number. */
  CONDUCTIVITY_DEPENDENCE_STRUCT d_k_struct; /* Acoustic wave number. */
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_k = &d_k_struct;

  dbl kR_inv; /* k*R product inverse	*/

  dbl grad_P[DIM]; /* Acoustic pressure. */

  int b, j, p, w;
  int var;

  if (d_q == NULL) {
    d_R = NULL;
    d_k = NULL;
  }

  R = acoustic_impedance(d_R, time);
  k = wave_number(d_k, time);
  kR_inv = 1. / (k * R);

  if (eqn == R_ACOUS_PREAL) {
    for (p = 0; p < VIM; p++) {
      grad_P[p] = fv->grad_apr[p];
    }
  } else if (eqn == R_ACOUS_PIMAG) {
    for (p = 0; p < VIM; p++) {
      grad_P[p] = fv->grad_api[p];
    }
  } else if (eqn == R_LIGHT_INTP) {
    for (p = 0; p < VIM; p++) {
      grad_P[p] = fv->grad_poynt[0][p];
    }
  } else if (eqn == R_LIGHT_INTM) {
    for (p = 0; p < VIM; p++) {
      grad_P[p] = fv->grad_poynt[1][p];
    }
  } else if (eqn == R_LIGHT_INTD) {
    for (p = 0; p < VIM; p++) {
      grad_P[p] = fv->grad_poynt[2][p];
    }
  } else {
    GOMA_EH(GOMA_ERROR, "Invalid Acoustic eqn");
  }

  for (p = 0; p < VIM; p++) {
    q[p] = kR_inv * grad_P[p];
  }

  var = ac_var;
  if (d_q != NULL && pd->v[pg->imtrx][var]) {
    for (p = 0; p < VIM; p++) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_q->P[p][j] = kR_inv * bf[var]->grad_phi[j][p];
      }
    }
  }

  var = TEMPERATURE;
  if (d_q != NULL && pd->v[pg->imtrx][var]) {
    for (p = 0; p < VIM; p++) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_q->T[p][j] = -kR_inv * kR_inv * (k * d_R->T[j] + R * d_k->T[j]) * grad_P[p];
      }
    }
  }

  var = MASS_FRACTION;
  if (d_q != NULL && pd->v[pg->imtrx][var]) {
    for (w = 0; w < pd->Num_Species_Eqn; w++) {
      for (p = 0; p < VIM; p++) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_q->C[p][w][j] = -kR_inv * kR_inv * (k * d_R->C[w][j] + R * d_k->C[w][j]) * grad_P[p];
        }
      }
    }
  }

  var = FILL;
  if (d_q != NULL && pd->v[pg->imtrx][var]) {
    for (p = 0; p < VIM; p++) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_q->F[p][j] = -kR_inv * kR_inv * (k * d_R->F[j] + R * d_k->F[j]) * grad_P[p];
      }
    }
  }

  if (d_q != NULL && pd->v[pg->imtrx][MESH_DISPLACEMENT1]) {
    for (p = 0; p < VIM; p++) {
      for (b = 0; b < WIM; b++) {
        var = MESH_DISPLACEMENT1 + b;
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          if (eqn == R_ACOUS_PREAL) {
            d_q->X[p][b][j] = kR_inv * fv->d_grad_apr_dmesh[p][b][j] -
                              kR_inv * kR_inv * (k * d_R->X[b][j] + R * d_k->X[b][j]) * grad_P[p];
          } else {
            d_q->X[p][b][j] = kR_inv * fv->d_grad_api_dmesh[p][b][j] -
                              kR_inv * kR_inv * (k * d_R->X[b][j] + R * d_k->X[b][j]) * grad_P[p];
          }
        }
      }
    }
  }
  return;
}

int assemble_ars_source(double ars_jump, double grad_jump) {
  int i, j, a, p, ii, ledof;
  int eqn, peqn, var, pvar;
  int dim;

  struct Basis_Functions *bfm;
  double wt, det_J, h3, phi_i, phi_j;

  double source;
  double omega, force, temp, acous_pgrad = 0, force1, temp1;

  double R; /* Acoustic impedance. */
  CONDUCTIVITY_DEPENDENCE_STRUCT d_R_struct;
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_R = &d_R_struct;

  double k; /* Acoustic wave number. */
  CONDUCTIVITY_DEPENDENCE_STRUCT d_k_struct;
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_k = &d_k_struct;

  double time = tran->time_value;

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

  dim = pd->Num_Dim;
  /* local force caused by Acoustic Reynolds Stress jump at interface	*/

  omega = upd->Acoustic_Frequency;
  temp = ars_jump / (4. * omega);
  force = (SQUARE(fv->apr) + SQUARE(fv->api)) * temp;

  R = acoustic_impedance(d_R, time);

  k = wave_number(d_k, time);

  for (p = 0; p < dim; p++) {
    acous_pgrad += fv->grad_api[p] * fv->grad_api[p];
    acous_pgrad += fv->grad_apr[p] * fv->grad_apr[p];
  }
  temp1 = grad_jump / (4. * omega) / SQUARE(k * R);
  force1 = acous_pgrad * temp1;
  force += force1;

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

          source = phi_i * lsi->normal[a] * force * lsi->delta;

          source *= det_J * wt * h3;
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

          source = phi_i * lsi->normal[a] * force * lsi->delta;

          source *= det_J * wt * h3;
          source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

          lec->R[LEC_R_INDEX(peqn, ii)] += source;
        }
      }
    }
  }

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
           * J_m_F
           */
          var = FILL;

          if (pd->v[pg->imtrx][var]) {
            pvar = upd->vp[pg->imtrx][var];

            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              phi_j = bf[var]->phi[j];

              source = phi_i * force *
                       (lsi->delta * lsi->d_normal_dF[a][j] + lsi->d_delta_dF[j] * lsi->normal[a]);

              source *= det_J * wt * h3;
              source *= pd->etm[pg->imtrx][eqn][LOG2_SOURCE];

              lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
            }
          }
          /*
           * Acoustic Pressure - Real
           */
          var = ACOUS_PREAL;

          if (pd->v[pg->imtrx][var]) {
            pvar = upd->vp[pg->imtrx][var];

            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              phi_j = bf[var]->phi[j];
              source = 2 * fv->apr * phi_j * temp;
              for (p = 0; p < dim; p++) {
                source += 2. * fv->grad_apr[p] * bf[var]->grad_phi[j][p] * temp1;
              }

              source *= phi_i * lsi->delta * lsi->normal[a];

              source *= det_J * wt * h3;
              source *= pd->etm[pg->imtrx][eqn][LOG2_SOURCE];

              lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
            }
          }
          /*
           * Acoustic Pressure - Imag
           */
          var = ACOUS_PIMAG;

          if (pd->v[pg->imtrx][var]) {
            pvar = upd->vp[pg->imtrx][var];

            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              phi_j = bf[var]->phi[j];
              source = 2 * fv->api * phi_j * temp;
              for (p = 0; p < dim; p++) {
                source += 2. * fv->grad_api[p] * bf[var]->grad_phi[j][p] * temp1;
              }

              source *= phi_i * lsi->delta * lsi->normal[a];

              source *= det_J * wt * h3;
              source *= pd->etm[pg->imtrx][eqn][LOG2_SOURCE];

              lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
            }
          }
        }
      }
    }
  }
  return (1);
}
/*
 * Viscous Dissipation Heating Source due to Acoustic Waves Model
 */

/*  _______________________________________________________________________  */
/*  _______________________________________________________________________  */
// end of em_diss_e_curlcurl_source
/*  _______________________________________________________________________  */

/* assemble_poynting -- assemble terms (Residual &| Jacobian) for Beer's law
 *            type light intensity absorption equations
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
 * Created:	Thu February 16, 2011 - Robert Secor
 *
 * Note: currently we do a "double load" into the addresses in the global
 *       "a" matrix and resid vector held in "esp", as well as into the
 *       local accumulators in "lec".
 *
 */

int assemble_poynting(double time, /* present time value */
                      double tt,   /* parameter to vary time integration from
                                    * explicit (tt = 1) to implicit (tt = 0) */
                      double dt,   /* current time step size */
                      const PG_DATA *pg_data,
                      const int py_eqn, /* eqn id and var id	*/
                      const int py_var) {
  int err, eqn, var, peqn, pvar, dim, p, b, w, w1, i, j, status, light_eqn = 0;
  int petrov = 0, Beers_Law = 0, explicit_deriv = 0;

  dbl P = 0;                 /* Light Intensity	*/
  dbl grad_P = 0, Psign = 0; /* grad intensity */

  dbl alpha; /* Acoustic Absorption */
  CONDUCTIVITY_DEPENDENCE_STRUCT d_alpha_struct;
  CONDUCTIVITY_DEPENDENCE_STRUCT *d_alpha = &d_alpha_struct;

  dbl diffusion = 0;
  //  dbl diff_a;
  dbl diff_b, diff_c, diff_d;

  dbl advection = 0, source = 0, advection_b = 0;
  /*
   * Galerkin weighting functions for i-th energy residuals
   * and some of their derivatives...
   */

  dbl phi_i, wt_func;

  /*
   * Interpolation functions for variables and some of their derivatives.
   */

  dbl phi_j;
  dbl grad_phi_j[DIM], grad_phi_i[DIM];

  dbl h3;          /* Volume element (scale factors). */
  dbl dh3dmesh_bj; /* Sensitivity to (b,j) mesh dof. */

  dbl det_J;

  dbl d_det_J_dmeshbj; /* for specified (b,j) mesh dof */
  dbl wt;

  const double *hsquared = pg_data->hsquared;
  const double *vcent = pg_data->v_avg; /* Average element velocity, which is the
          centroid velocity for Q2 and the average
          of the vertices for Q1. It comes from
          the routine "element_velocity." */

  /* SUPG variables */
  dbl h_elem = 0, h_elem_inv = 0, h_elem_inv_deriv = 0;
  dbl d_wt_func;

  double vconv[MAX_PDIM];     /*Calculated convection velocity */
  double vconv_old[MAX_PDIM]; /*Calculated convection velocity at previous time*/
  CONVECTION_VELOCITY_DEPENDENCE_STRUCT d_vconv_struct;
  CONVECTION_VELOCITY_DEPENDENCE_STRUCT *d_vconv = &d_vconv_struct;
  struct Species_Conservation_Terms s_terms;
  double d_drop_source[MAX_CONC];
  /*
   *    Radiative transfer equation variables - connect to input file someday
   */
  double svect[DIM] = {0., -1., 0.};
  double v_grad[DIM] = {0., 0., 0.};
  double mucos = 1.0;
  double diff_const = 1. / LITTLE_PENALTY;
  double time_source = 0., d_time_source = 0.;

  zero_structure(&s_terms, sizeof(struct Species_Conservation_Terms), 1);

  /*   static char yo[] = "assemble_poynting";*/
  status = 0;
  /*
   * Unpack variables from structures for local convenience...
   */
  dim = pd->Num_Dim;
  eqn = py_eqn;

  /*
   * Bail out fast if there's nothing to do...
   */
  if (!pd->e[pg->imtrx][eqn]) {
    return (status);
  }

  wt = fv->wt;           /* Gauss point weight. */
  h3 = fv->h3;           /* Differential volume element. */
  det_J = bf[eqn]->detJ; /* Really, ought to be mesh eqn. */

  /*
   * Material property constants, etc. Any variations at this Gauss point
   * are calculated with user-defined subroutine options.
   *
   * Properties to consider here are R and k. R and k we will allow to vary
   * with temperature, spatial coordinates, and species concentration.
   */

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

  /* get the convection velocity (it's different for arbitrary and
     lagrangian meshes) */
  if (cr->MeshMotion == ARBITRARY || cr->MeshMotion == LAGRANGIAN ||
      cr->MeshMotion == DYNAMIC_LAGRANGIAN) {
    err = get_convection_velocity(vconv, vconv_old, d_vconv, dt, tt);
    GOMA_EH(err, "Error in calculating effective convection velocity");
  }
  /* end Petrov-Galerkin addition */

  alpha = light_absorption(d_alpha, time);

  switch (py_eqn) {
  case R_LIGHT_INTP:
    light_eqn = 0;
    Psign = 1.;
    break;
  case R_LIGHT_INTM:
    light_eqn = 1;
    Psign = -1.;
    break;
  case R_LIGHT_INTD:
    light_eqn = 2;
    Psign = 0.;
    break;
  case R_RESTIME:
    petrov = (mp->Rst_func_supg > 0 ? 1 : 0);
    break;
  default:
    GOMA_EH(GOMA_ERROR, "light intensity equation");
    break;
  }
  if (py_eqn >= R_LIGHT_INTP && py_eqn <= R_LIGHT_INTD) {
    Beers_Law = 1;
    P = fv->poynt[light_eqn];
    grad_P = 0.;
    for (i = 0; i < dim; i++) {
      v_grad[i] = fv->grad_poynt[light_eqn][i];
      grad_P += svect[i] * v_grad[i];
    }
  } else {
    Beers_Law = 0;
    P = fv->restime;
    diff_const = mp->Rst_diffusion;
    for (i = 0; i < dim; i++) {
      v_grad[i] = fv->grad_restime[i];
    }
    time_source = 0.;
    d_time_source = 0.;
    switch (mp->Rst_funcModel) {
    case CONSTANT:
      time_source = mp->Rst_func;
      d_time_source = 0.;
      break;
    case LINEAR_TIMETEMP:
      if (pd->gv[R_ENERGY] && (fv->T > upd->Process_Temperature)) {
        time_source = mp->Rst_func * (fv->T - upd->Process_Temperature);
        d_time_source = mp->Rst_func;
      }
      break;
    case EXPONENTIAL_TIMETEMP:
      if (pd->gv[R_ENERGY] && (fv->T > upd->Process_Temperature)) {
        time_source = exp(mp->Rst_func * (fv->T - upd->Process_Temperature));
        d_time_source = mp->Rst_func * time_source;
      }
      break;
    case DROP_EVAP:
      if (pd->gv[R_MASS]) {
        double init_radius = 0.010, num_density = 10., denom = 1;
        err = get_continuous_species_terms(&s_terms, time, tt, dt, hsquared);
        GOMA_EH(err, "problem in getting the species terms");

        for (w = 0; w < pd->Num_Species_Eqn; w++) {
          if (mp->SpeciesSourceModel[w] == DROP_EVAP) {
            init_radius = mp->u_species_source[w][1];
            num_density = mp->u_species_source[w][2];
            P = MAX(DBL_SMALL, fv->restime);
/**  Droplet radius or volume formulation **/
#if 1
            denom = MAX(DBL_SMALL, num_density * 4 * M_PIE * CUBE(init_radius) * SQUARE(P));
#else
            denom = MAX(DBL_SMALL, num_density * (4. / 3.) * M_PIE * CUBE(init_radius));
#endif
          }
          time_source -= mp->molar_volume[w] * s_terms.MassSource[w] / denom;
#if 1
          if (P > DBL_SMALL) {
            d_time_source += mp->molar_volume[w] * s_terms.MassSource[w] / denom * 2. / P;
          }
#endif
          d_drop_source[w] = mp->Rst_func * mp->molar_volume[w] / denom;
        }
        time_source *= mp->Rst_func;
        d_time_source *= mp->Rst_func;
      }
      explicit_deriv = 1;
      break;
    default:
      GOMA_EH(GOMA_ERROR, "Residence Time Weight Function");
      break;
    }
  }
  /*
   * Residuals___________________________________________________________
   */

  if (af->Assemble_Residual) {
    eqn = py_eqn;
    peqn = upd->ep[pg->imtrx][eqn];
    var = py_var;

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {

      /* this is an optimization for xfem */
      if (xfem != NULL) {
        int xfem_active, extended_dof, base_interp, base_dof;
        xfem_dof_state(i, pd->i[pg->imtrx][eqn], ei[pg->imtrx]->ielem_shape, &xfem_active,
                       &extended_dof, &base_interp, &base_dof);
        if (extended_dof && !xfem_active)
          continue;
      }
      phi_i = bf[eqn]->phi[i];
      wt_func = (1. - mp->Rst_func_supg) * phi_i;
      /* add Petrov-Galerkin terms as necessary */
      if (petrov) {
        for (p = 0; p < dim; p++) {
          wt_func += mp->Rst_func_supg * h_elem_inv * vconv[p] * bf[eqn]->grad_phi[i][p];
        }
      }

      advection = 0.;
      if ((pd->e[pg->imtrx][eqn] & T_ADVECTION) && !Beers_Law) {
        source = -time_source;
        for (p = 0; p < dim; p++) {
          grad_phi_i[p] = bf[var]->grad_phi[i][p];
          advection += wt_func * vconv[p] * v_grad[p];
          advection += diff_const * grad_phi_i[p] * v_grad[p];
        }
        advection += wt_func * source;

        advection *= det_J * wt;
        advection *= h3;
        advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
      }
      diffusion = 0.;
      if ((pd->e[pg->imtrx][eqn] & T_DIFFUSION) && Beers_Law) {

        diffusion += phi_i * (mucos * grad_P + Psign * alpha * P);
        diffusion *= det_J * wt;
        diffusion *= h3;
        diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
      }

      lec->R[LEC_R_INDEX(peqn, i)] += diffusion + advection;
    }
  }

  /*
   * Jacobian terms...
   */

  if (af->Assemble_Jacobian) {
    eqn = py_eqn;
    peqn = upd->ep[pg->imtrx][eqn];
    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      /* this is an optimization for xfem */
      if (xfem != NULL) {
        int xfem_active, extended_dof, base_interp, base_dof;
        xfem_dof_state(i, pd->i[pg->imtrx][eqn], ei[pg->imtrx]->ielem_shape, &xfem_active,
                       &extended_dof, &base_interp, &base_dof);
        if (extended_dof && !xfem_active)
          continue;
      }
      phi_i = bf[eqn]->phi[i];
      wt_func = (1. - mp->Rst_func_supg) * phi_i;
      if (petrov) {
        for (p = 0; p < dim; p++) {
          wt_func += h_elem_inv * vconv[p] * bf[eqn]->grad_phi[i][p];
        }
      }

      /*
       * Set up some preliminaries that are needed for the (a,i)
       * equation for bunches of (b,j) column variables...
       */
      for (p = 0; p < dim; p++) {
        grad_phi_i[p] = bf[eqn]->grad_phi[i][p];
      }

      /*
       * J_e_ap
       */
      var = py_var;
      if (pd->v[pg->imtrx][var]) {
        pvar = upd->vp[pg->imtrx][var];
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          phi_j = bf[var]->phi[j];

          for (p = 0; p < VIM; p++) {
            grad_phi_j[p] = bf[var]->grad_phi[j][p];
          }

          advection = 0.;
          if ((pd->e[pg->imtrx][eqn] & T_ADVECTION) && !Beers_Law) {
            for (p = 0; p < dim; p++) {
              advection += wt_func * vconv[p] * grad_phi_j[p];
              advection += diff_const * grad_phi_i[p] * grad_phi_j[p];
            }
            if (explicit_deriv) {
              for (w = 0; w < pd->Num_Species_Eqn; w++) {
                advection += wt_func * d_drop_source[w] * s_terms.d_MassSource_drst[w][j];
              }
              advection -= wt_func * d_time_source * phi_j;
            }

            advection *= det_J * wt;
            advection *= h3;
            advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
          }
          diffusion = 0.;
          if ((pd->e[pg->imtrx][eqn] & T_DIFFUSION) && Beers_Law) {
            for (p = 0; p < VIM; p++) {
              diffusion += grad_phi_j[p] * svect[p];
            }
            diffusion *= phi_i * mucos;
            diffusion += phi_i * Psign * alpha * phi_j;
            diffusion *= h3 * det_J * wt;
            diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
          }

          lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += diffusion + advection;
        }
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

          advection = diffusion = 0;
          if ((pd->e[pg->imtrx][eqn] & T_ADVECTION) && !Beers_Law) {
            /*	             advection += (diff_const/fv->T)*grad_phi_i[p]*v_grad[p];  */
            if (explicit_deriv) {
              for (w = 0; w < pd->Num_Species_Eqn; w++) {
                advection += wt_func * d_drop_source[w] * s_terms.d_MassSource_dT[w][j];
              }
            } else {
              advection += -wt_func * d_time_source * phi_j;
            }

            advection *= det_J * wt;
            advection *= h3;
            advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
          }
          if ((pd->e[pg->imtrx][eqn] & T_DIFFUSION) && Beers_Law) {
            diffusion = phi_i * d_alpha->T[j] * P;
            diffusion *= det_J * wt;
            diffusion *= h3;
            diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
          }

          lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += diffusion + advection;
        }
      }

      /*
       * J_e_Pressure
       */
      var = PRESSURE;
      if (pd->v[pg->imtrx][var]) {
        pvar = upd->vp[pg->imtrx][var];
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          phi_j = bf[var]->phi[j];

          for (p = 0; p < VIM; p++) {
            grad_phi_j[p] = bf[var]->grad_phi[j][p];
          }

          advection = diffusion = 0;
          if ((pd->e[pg->imtrx][eqn] & T_ADVECTION) && !Beers_Law) {
            if (explicit_deriv) {
              advection = 0;
              for (w = 0; w < pd->Num_Species_Eqn; w++) {
                advection += wt_func * d_drop_source[w] * s_terms.d_MassSource_dP[w][j];
              }
            } else {
              advection = -wt_func * d_time_source * phi_j;
            }

            advection *= det_J * wt;
            advection *= h3;
            advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
          }
          if ((pd->e[pg->imtrx][eqn] & T_DIFFUSION) && Beers_Law) {
            diffusion = phi_i * d_alpha->T[j] * P;
            diffusion *= det_J * wt;
            diffusion *= h3;
            diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
          }

          lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += diffusion + advection;
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

            advection = diffusion = 0;

            if ((pd->e[pg->imtrx][eqn] & T_ADVECTION) && !Beers_Law) {
              for (p = 0; p < dim; p++) {
                advection += wt_func * vconv[p] * fv->d_grad_restime_dmesh[p][b][j];
                advection += wt_func * d_vconv->X[p][b][j] * v_grad[p];
                advection += diff_const * grad_phi_i[p] * fv->d_grad_restime_dmesh[p][b][j];
              }
              advection *= det_J * wt;
              advection *= h3;
              advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
            }
            /*
             * multiple parts:
             * 	diff_a = Int(...d(grad_phi_i)/dmesh.q h3 |Jv|)
             *	diff_b = Int(...grad_phi_i.d(q)/dmesh h3 |Jv|)
             *	diff_c = Int(...grad_phi_i.q h3 d(|Jv|)/dmesh)
             *	diff_d = Int(...grad_phi_i.q dh3/dmesh |Jv|  )
             */
            if ((pd->e[pg->imtrx][eqn] & T_DIFFUSION) && Beers_Law) {
              diff_b = 0.;
              for (p = 0; p < VIM; p++) {
                diff_b += svect[p] * fv->d_grad_poynt_dmesh[light_eqn][p][b][j];
              }
              diff_b *= mucos;
              diff_b += Psign * d_alpha->X[b][j] * P;
              diff_b *= phi_i * det_J * h3 * wt;

              diff_c = phi_i * (mucos * grad_P + Psign * alpha * P);
              diff_c *= d_det_J_dmeshbj * h3 * wt;

              diff_d = phi_i * (mucos * grad_P + Psign * alpha * P);
              diff_d *= det_J * dh3dmesh_bj * wt;

              diffusion = diff_b + diff_c + diff_d;

              diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
            }

            lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += diffusion + advection;
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

            advection = diffusion = 0;
            if ((pd->e[pg->imtrx][eqn] & T_DIFFUSION) && Beers_Law) {
              diffusion = phi_i * Psign * d_alpha->C[w][j] * P;
              diffusion *= det_J * wt;
              diffusion *= h3;
              diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];

              lec->J[LEC_J_INDEX(peqn, MAX_PROB_VAR + w, i, j)] += diffusion;
            }
            if ((pd->e[pg->imtrx][eqn] & T_DIFFUSION) && !Beers_Law) {
              advection = 0;
              for (w1 = 0; w1 < pd->Num_Species_Eqn; w1++) {
                advection += wt_func * d_drop_source[w] * s_terms.d_MassSource_dc[w][w1][j];
              }

              advection *= det_J * wt;
              advection *= h3;
              advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
              lec->J[LEC_J_INDEX(peqn, MAX_PROB_VAR + w, i, j)] += advection;
            }
          }
        }
      }
      /*
       * J_e_v
       */
      for (b = 0; b < dim; b++) {
        var = VELOCITY1 + b;
        if (pd->v[pg->imtrx][var]) {
          pvar = upd->vp[pg->imtrx][var];
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];

            advection = 0;
            advection_b = 0;
            if ((pd->e[pg->imtrx][eqn] & T_ADVECTION) && !Beers_Law) {
              for (p = 0; p < dim; p++) {
                advection += wt_func * phi_j * v_grad[p];
              }

              advection *= det_J * wt;
              advection *= h3;
              advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];

              d_wt_func = h_elem_inv * d_vconv->v[b][b][j] * grad_phi_i[b] +
                          h_elem_inv_deriv * vconv[b] * grad_phi_i[b];
              d_wt_func *= mp->Rst_func_supg;

              for (p = 0; p < dim; p++) {
                advection_b += vconv[p] * v_grad[p];
              }

              advection_b *= d_wt_func * det_J * wt;
              advection_b *= h3;
              advection_b *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];

              lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += advection + advection_b;
            }
          }
        }
      }
    }
  }

  return (status);
} /* end of assemble_poynting */

void restime_nobc_surf(

    double func[MAX_PDIM], double d_func[MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE])
/******************************************************************************
 *
 *  Function which calculates the surface integral for the "no bc" restime boundary
 *condition
 *
 ******************************************************************************/
{

  /* Local variables */

  int j, b, p;
  int var, dim;

  double v_grad[DIM], diff_const;

  /***************************** EXECUTION BEGINS *******************************/

  if (af->Assemble_LSA_Mass_Matrix)
    return;

  dim = pd->Num_Dim;

  diff_const = mp->Rst_diffusion;

  for (j = 0; j < dim; j++) {
    v_grad[j] = fv->grad_restime[j];
  }

  if (af->Assemble_Jacobian) {
    var = RESTIME;
    if (pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        for (p = 0; p < dim; p++) {
          d_func[0][var][j] += fv->snormal[p] * diff_const * bf[var]->grad_phi[j][p];
          /*		  d_func[0][var][j] +=
           * fv->snormal[p]*diff_const*(-1./fv->restime)*v_grad[p];  */
        }
      }
    }

    if (pd->v[pg->imtrx][MESH_DISPLACEMENT1]) {
      for (b = 0; b < dim; b++) {
        var = MESH_DISPLACEMENT1 + b;
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          for (p = 0; p < dim; p++) {
            d_func[0][var][j] += diff_const * (fv->snormal[p] * fv->d_grad_restime_dmesh[p][b][j] +
                                               v_grad[p] * fv->dsnormal_dx[p][b][j]);
          }
        }
      }
    }
  }

  /* Calculate the residual contribution	     			     */
  for (p = 0; p < dim; p++) {
    *func += fv->snormal[p] * diff_const * v_grad[p];
  }

  return;
} /* END of routine restime_nobc_surf                                               */
/*****************************************************************************/
