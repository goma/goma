#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* GOMA include files */
#define GOMA_MM_FILL_MOMENTUM_C
#include "ac_particles.h"
#include "ad_turbulence.h"
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
#include "polymer_time_const.h"
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

/* assemble_momentum -- assemble terms (Residual &| Jacobian) for momentum eqns
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

int assemble_momentum(dbl time,       /* current time */
                      dbl tt,         /* parameter to vary time integration from
                                         explicit (tt = 1) to implicit (tt = 0) */
                      dbl dt,         /* current time step size */
                      dbl h_elem_avg, /* average global element size for PSPG*/
                      const PG_DATA *pg_data,
                      double xi[DIM], /* Local stu coordinates */
                      const Exo_DB *exo) {
#ifdef DEBUG_MOMENTUM_JAC
  int adx;
#endif

  int dim;
  int i, j, jk, p, q, a, b, c;

  int ledof, eqn, var, ii, peqn, pvar, w;

  int *pde = pd->e[pg->imtrx];
  int *pdv = pd->v[pg->imtrx];

  int status;
  struct Basis_Functions *bfm;

  dbl zero[3] = {0.0, 0.0, 0.0}; /* A zero array, for convenience. */

  /*dbl v_dot[DIM];*/ /* time derivative of velocity field. */
  dbl *v_dot;         /* time derivative of velocity field. */

  /*dbl x_dot[DIM];*/ /* current position field derivative wrt time. */
  dbl *x_dot;         /* current position field derivative wrt time. */

  dbl h3;          /* Volume element (scale factors). */
  dbl dh3dmesh_bj; /* Sensitivity to (b,j) mesh dof. */

  /* field variables */
  /*dbl grad_v[DIM][DIM];*/ /* Gradient of v. */
  dbl *grad_v[DIM];
  dbl *v = fv->v;

  dbl rho; /* Density. */

  dbl f[DIM];                                  /* Body force. */
  MOMENTUM_SOURCE_DEPENDENCE_STRUCT df_struct; /* Body force dependence */
  MOMENTUM_SOURCE_DEPENDENCE_STRUCT *df = &df_struct;

  dbl det_J; /* determinant of element Jacobian */

  dbl d_det_J_dmesh_bj; /* for specific (b,j) mesh dof */

  dbl d_area;

  dbl mass; /* For terms and their derivatives */
  dbl advection;
  dbl advection_a, advection_b, advection_c;
  dbl porous;
  dbl diffusion;
  dbl diff_a, diff_b, diff_c;
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
  dbl phi_i;
  dbl(*grad_phi_i_e_a)[DIM] = NULL;
  /*
   * Interpolation functions for variables and some of their derivatives.
   */

  dbl phi_j;
  dbl(*d_grad_phi_i_e_a_dmesh)[DIM][DIM][MDE] = NULL;

  double *phi_i_vector, *phi_j_vector;

  dbl Pi[DIM][DIM];
  STRESS_DEPENDENCE_STRUCT d_Pi_struct;
  STRESS_DEPENDENCE_STRUCT *d_Pi = &d_Pi_struct;

  dbl wt;

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

  /* density derivatives */
  DENSITY_DEPENDENCE_STRUCT d_rho_struct; /* density dependence */
  DENSITY_DEPENDENCE_STRUCT *d_rho = &d_rho_struct;

  int v_g[DIM][DIM];
  int v_s[MAX_MODES][DIM][DIM];
  int mode;

  /* Variables used for the modified fluid momentum equations when looking
   * at the particle momentum model.  Note that pd->MomentumFluxModel ==
   * SUSPENSION_PM when the Buyevich particle momentum equations are active.
   */
  int particle_momentum_on; /* Boolean for particle momentum eq.'s active */
  int species = -1;         /* Species number of particle phase */
  double p_vol_frac = 0;    /* Particle volume fraction (phi) */
  double ompvf = 1;         /* 1 - p_vol_frac "One Minus Particle */
  /*    Volume Fraction" */
  int transient_run = pd->TimeIntegration != STEADY;
  int mass_on;
  int advection_on = 0;
  int source_on = 0;
  int diffusion_on = 0;
  int porous_brinkman_on = 0;
  int mesh_disp_on = 0;

  dbl mass_etm, advection_etm, diffusion_etm, source_etm, porous_brinkman_etm;

  double *J = NULL;
  double *R;

  /*
   * Petrov-Galerkin weighting functions for i-th residuals
   * and some of their derivatives...
   */

  dbl wt_func;

  /* SUPG variables */
  momentum_tau_terms supg_terms;

  // Continuity stabilization
  dbl continuity_stabilization;
  dbl cont_gls;
  CONT_GLS_DEPENDENCE_STRUCT d_cont_gls_struct;
  CONT_GLS_DEPENDENCE_STRUCT *d_cont_gls = &d_cont_gls_struct;

  int *n_dof = NULL;
  int dof_map[MDE];

  status = 0;

  /*
   * Unpack variables from structures for local convenience...
   */

  eqn = R_MOMENTUM1;

  /*
   * Bail out fast if there's nothing to do...
   */
  if (!pd->e[pg->imtrx][eqn]) {
    return (status);
  }

  dim = pd->Num_Dim;

  wt = fv->wt;

  det_J = bf[eqn]->detJ; /* Really, ought to be mesh eqn. */

  h3 = fv->h3; /* Differential volume element (scales). */

  d_area = det_J * wt * h3;

  dbl supg = 0.;

  if (mp->Mwt_funcModel == GALERKIN) {
    supg = 0.;
  } else if (mp->Mwt_funcModel == SUPG || mp->Mwt_funcModel == SUPG_GP ||
             mp->Mwt_funcModel == SUPG_SHAKIB) {
    supg = mp->Mwt_func;
  }

  /*** Density ***/
  rho = density(d_rho, time);

  if (supg != 0.) {
    tau_momentum_shakib(&supg_terms, dim, dt, FALSE);
  }
  /* end Petrov-Galerkin addition */

  if (pd->gv[POLYMER_STRESS11]) {
    (void)stress_eqn_pointer(v_s);
  }

  if (pd->gv[VELOCITY_GRADIENT11]) {
    v_g[0][0] = VELOCITY_GRADIENT11;
    v_g[0][1] = VELOCITY_GRADIENT12;
    v_g[1][0] = VELOCITY_GRADIENT21;
    v_g[1][1] = VELOCITY_GRADIENT22;
    v_g[0][2] = VELOCITY_GRADIENT13;
    v_g[1][2] = VELOCITY_GRADIENT23;
    v_g[2][0] = VELOCITY_GRADIENT31;
    v_g[2][1] = VELOCITY_GRADIENT32;
    v_g[2][2] = VELOCITY_GRADIENT33;
  }

  /* Set up variables for particle/fluid momentum coupling.
   */
  particle_momentum_on = 0;
  if (pd->gv[R_PMOMENTUM1]) {
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
  }

  /*
   * Material property constants, etc. Any variations for this
   * Gauss point were evaluated in load_material_properties().
   */

  if (pd->e[pg->imtrx][eqn] & T_POROUS_BRINK) {
    /* Load up remaining parameters for the Brinkman Equation. */
    por = mp->porosity;
    vis = mp->FlowingLiquid_viscosity;
    /* Load variable FlowingLiquid_viscosity */
    sc = mp->Inertia_coefficient;

    /*
     * Special lubrication case
     */
    if (mp->FSIModel > 0) {

      /* Store proper Gauss weights */
      wt = fv->wt;
      h3 = fv->h3;

      /* Prepare geometry */
      n_dof = (int *)array_alloc(1, MAX_VARIABLE_TYPES, sizeof(int));
      lubrication_shell_initialize(n_dof, dof_map, -1, xi, exo, 0);

      /* Adjust FEM weights */
      dim = pd->Num_Dim;
      fv->wt = wt;
      det_J = fv->sdet;

      /*
       * Calculate the velocity from the velocity calculator, of course!
       * But then why do I even need this equation?
       * I DON'T!!!
       */
      calculate_lub_q_v(R_LUBP, time, dt, xi, exo);

      fv->wt = wt; /*load_neighbor_var_data screws fv->wt up */
    }

    if (mp->PorousMediaType != POROUS_BRINKMAN)
      GOMA_WH(-1, "Set Porous term multiplier in continuous medium");

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
  } else {
    por = 1.;
    per = 1.;
    vis = mp->viscosity;
    sc = 0.;
  }

  /*   eqn = R_MOMENTUM1; */
  /*
   * Field variables...
   */
  if (pd->gv[MESH_DISPLACEMENT1]) {
    mesh_disp_on = 1;
  }

  x_dot = zero;
  if (transient_run && mesh_disp_on)
    x_dot = fv_dot->x;

  if (transient_run)
    v_dot = fv_dot->v;
  else
    v_dot = zero;

    /* for porous media stuff */
#ifdef DO_NO_UNROLL
  speed = 0.0;
  for (a = 0; a < WIM; a++) {
    speed += v[a] * v[a];
  }
  speed = sqrt(speed);

  for (a = 0; a < VIM; a++)
    grad_v[a] = fv->grad_v[a];
#else
  speed = 0.0;
  speed += v[0] * v[0] + v[1] * v[1];
  if (WIM == 3)
    speed += v[2] * v[2];

  speed = sqrt(speed);

  grad_v[0] = fv->grad_v[0];
  grad_v[1] = fv->grad_v[1];
  if (VIM == 3)
    grad_v[2] = fv->grad_v[2];
#endif

  /*
   * Calculate the momentum stress tensor at the current gauss point
   */
  if (af->Assemble_Jacobian) {
    fluid_stress(Pi, d_Pi);
  } else {
    fluid_stress(Pi, NULL);
  }

  (void)momentum_source_term(f, df, time);

  // Call continuity stabilization if desired
  if (Cont_GLS) {
    calc_cont_gls(&cont_gls, d_cont_gls, time, pg_data);
  }

  double pspg[3] = {0.0, 0.0, 0.0};
  PSPG_DEPENDENCE_STRUCT d_pspg;

  if (upd->PSPG_advection_correction) {
    calc_pspg(pspg, &d_pspg, time, tt, dt, pg_data);
  }

  /*
   * Residuals_________________________________________________________________
   */

  if (af->Assemble_Residual) {
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

      R = &(lec->R[LEC_R_INDEX(peqn, 0)]);
      phi_i_vector = bfm->phi;

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

          phi_i = phi_i_vector[i];
          /* only use Petrov Galerkin on advective term - if required */
          wt_func = phi_i;
          /* add Petrov-Galerkin terms as necessary */
          if (supg != 0.) {
            for (p = 0; p < dim; p++) {
              wt_func += supg * supg_terms.tau * v[p] * bfm->grad_phi[i][p];
            }
          }
          grad_phi_i_e_a = bfm->grad_phi_e[i][a];

          /* this is an optimization for xfem */
          if (xfem != NULL) {
            int xfem_active, extended_dof, base_interp, base_dof;
            xfem_dof_state(i, pd->i[pg->imtrx][eqn], ei[pg->imtrx]->ielem_shape, &xfem_active,
                           &extended_dof, &base_interp, &base_dof);
            if (extended_dof && !xfem_active)
              continue;
          }

          mass = 0.;
          if (transient_run) {
            if (mass_on) {
              mass = v_dot[a] * rho;
              mass *= -wt_func * d_area;
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

            if (upd->PSPG_advection_correction) {
              advection -= pspg[0] * grad_v[0][a];
              advection -= pspg[1] * grad_v[1][a];
              if (WIM == 3)
                advection -= pspg[2] * grad_v[2][a];
            }

            advection *= rho;
            advection *= -wt_func * d_area;
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
            if (vis != 0.) {
              porous = v[a] * (rho * sc * speed / sqrt(per) + vis / per);
              porous *= -phi_i * d_area;
              porous *= porous_brinkman_etm;
            } else if (mp->FSIModel > 0) {
              porous = v[a] - LubAux->v_avg[a];
              porous *= phi_i * wt * fv->sdet * h3;
              porous *= porous_brinkman_etm;
            } else if (vis == 0. && mp->viscosity == 0.) {
              GOMA_EH(GOMA_ERROR,
                      "cannot have both flowing liquid viscosity and mp->viscosity equal to zero");
            }
          }

          diffusion = 0.;
          if (diffusion_on) {
#ifdef DO_NO_UNROLL
            for (p = 0; p < VIM; p++) {
              for (q = 0; q < VIM; q++) {
                diffusion += grad_phi_i_e_a[p][q] * Pi[q][p];
              }
            }
#else
            diffusion += grad_phi_i_e_a[0][0] * Pi[0][0];
            diffusion += grad_phi_i_e_a[1][1] * Pi[1][1];
            diffusion += grad_phi_i_e_a[1][0] * Pi[0][1];
            diffusion += grad_phi_i_e_a[0][1] * Pi[1][0];
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

          /*
           * Source term...
           */

          source = 0.0;
          if (source_on) {
            source += f[a];
            source *= wt_func * d_area;
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

          /*
           * Add contributions to this residual (globally into Resid, and
           * locally into an accumulator)
           */

          /*lec->R[LEC_R_INDEX(peqn,ii)] += mass + advection + porous + diffusion + source;*/
          R[ii] += mass + advection + porous + diffusion + source + continuity_stabilization;

#ifdef DEBUG_MOMENTUM_RES
          printf("R_m[%d][%d] += %10f %10f %10f %10f %10f\n", a, i, mass, advection, porous,
                 diffusion, source);
#endif /* DEBUG_MOMENTUM_RES */

        } /*end if (active_dofs) */
      }   /* end of for (i=0,ei[pg->imtrx]->dofs...) */
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

      phi_i_vector = bfm->phi;

      for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
        ii = ei[pg->imtrx]->lvdof_to_row_lvdof[eqn][i];

        ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][i];
        if (ei[pg->imtrx]->active_interp_ledof[ledof]) {
          /*
           *  Here is where we figure out whether the row is to placed in
           *  the normal spot (e.g., ii = i), or whether a boundary condition
           *  require that the volumetric contribution be stuck in another
           *  ldof pertaining to the same variable type.
           */
          ii = ei[pg->imtrx]->lvdof_to_row_lvdof[eqn][i];

          phi_i = phi_i_vector[i];

          /* Assign pointers into the bf structure */

          grad_phi_i_e_a = bfm->grad_phi_e[i][a];

          /* this is an optimization for xfem */
          if (xfem != NULL) {
            int xfem_active, extended_dof, base_interp, base_dof;
            xfem_dof_state(i, pd->i[pg->imtrx][eqn], ei[pg->imtrx]->ielem_shape, &xfem_active,
                           &extended_dof, &base_interp, &base_dof);
            if (extended_dof && !xfem_active)
              continue;
          }

          d_grad_phi_i_e_a_dmesh = bfm->d_grad_phi_e_dmesh[i][a];

#ifdef DEBUG_MOMENTUM_JAC
          if (ei[pg->imtrx]->ielem == 0) {
            printf("\nASSEMBLE_MOMENTUM, a = %d, dof = %d\n", a, i);
            printf("\tphi_i = %g\n", phi_i);
          }
#endif /* DEBUG_MOMENTUM_JAC */

          wt_func = phi_i;
          /* add Petrov-Galerkin terms as necessary */
          if (supg != 0.) {
            for (p = 0; p < dim; p++) {
              wt_func += supg * supg_terms.tau * v[p] * bfm->grad_phi[i][p];
            }
          }

          /*
           * J_m_T
           */
          if (pdv[TEMPERATURE]) {
            var = TEMPERATURE;
            pvar = upd->vp[pg->imtrx][var];
            phi_j_vector = bf[var]->phi;
            J = &(lec->J[LEC_J_INDEX(peqn, pvar, ii, 0)]);

            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              phi_j = phi_j_vector[j];

              dbl mass_supg = 0;
              dbl advection_supg = 0;
              dbl source_supg = 0;
              if (supg != 0.) {
                dbl d_wt_func = 0;
                for (p = 0; p < dim; p++) {
                  d_wt_func += supg * supg_terms.d_tau_dT[j] * v[p] * bfm->grad_phi[i][p];
                }
                if (transient_run) {
                  if (mass_on) {
                    mass_supg = v_dot[a] * rho;
                    mass_supg *= -d_wt_func * d_area;
                    mass_supg *= mass_etm;
                  }

                  if (porous_brinkman_on) {
                    mass_supg /= por;
                  }

                  if (particle_momentum_on) {
                    mass_supg *= ompvf;
                  }
                }

                if (advection_on) {
#ifdef DO_NO_UNROLL
                  for (p = 0; p < WIM; p++) {
                    advection_supg += (v[p] - x_dot[p]) * grad_v[p][a];
                  }
#else
                  advection_supg += (v[0] - x_dot[0]) * grad_v[0][a];
                  advection_supg += (v[1] - x_dot[1]) * grad_v[1][a];
                  if (WIM == 3)
                    advection_supg += (v[2] - x_dot[2]) * grad_v[2][a];
#endif

                  if (upd->PSPG_advection_correction) {
                    advection_supg -= pspg[0] * grad_v[0][a];
                    advection_supg -= pspg[1] * grad_v[1][a];
                    if (WIM == 3)
                      advection_supg -= pspg[2] * grad_v[2][a];
                  }

                  advection_supg *= rho;
                  advection_supg *= -d_wt_func * d_area;
                  advection_supg *= advection_etm;

                  if (porous_brinkman_on) {
                    por2 = por * por;
                    advection_supg /= por2;
                  }

                  if (particle_momentum_on) {
                    advection_supg *= ompvf;
                  }
                }

                /*
                 * Source term...
                 */

                if (source_on) {
                  source_supg += f[a];
                  source_supg *= d_wt_func * d_area;
                  source_supg *= source_etm;
                }
              }

              mass = 0.;

              if (transient_run) {
                if (mass_on) {
                  mass = d_rho->T[j] * v_dot[a];
                  mass *= -wt_func * d_area;
                  mass *= mass_etm;
                }

                /* if porous flow is considered. KSC on 5/10/95 */
                if (porous_brinkman_on) {
                  mass /= por;
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

                advection *= -wt_func * d_rho->T[j] * d_area;
                advection *= advection_etm;

                if (porous_brinkman_on) {
                  por2 = por * por;
                  advection /= por2;
                }
              }

              porous = 0.;
              if (porous_brinkman_on) {
                if (vis != 0.) {
                  porous = v[a] * (d_rho->T[j] * sc * speed / sqrt(per));
                  porous += v[a] * (d_flow_mu->T[j] / per);
                  porous *= -phi_i * d_area;
                  porous *= porous_brinkman_etm;
                } else {
                  porous = 0.;
                }
              }

              diffusion = 0.;
              if (diffusion_on) {
#ifdef DO_NO_UNROLL
                for (p = 0; p < VIM; p++) {
                  for (q = 0; q < VIM; q++) {
                    diffusion += grad_phi_i_e_a[p][q] * d_Pi->T[q][p][j];
                  }
                }
#else
                diffusion += grad_phi_i_e_a[0][0] * d_Pi->T[0][0][j];
                diffusion += grad_phi_i_e_a[1][1] * d_Pi->T[1][1][j];
                diffusion += grad_phi_i_e_a[1][0] * d_Pi->T[0][1][j];
                diffusion += grad_phi_i_e_a[0][1] * d_Pi->T[1][0][j];
                if (VIM == 3) {
                  diffusion += grad_phi_i_e_a[2][0] * d_Pi->T[0][2][j];
                  diffusion += grad_phi_i_e_a[2][1] * d_Pi->T[1][2][j];
                  diffusion += grad_phi_i_e_a[2][2] * d_Pi->T[2][2][j];
                  diffusion += grad_phi_i_e_a[0][2] * d_Pi->T[2][0][j];
                  diffusion += grad_phi_i_e_a[1][2] * d_Pi->T[2][1][j];
                }
#endif
                diffusion *= -d_area;
                diffusion *= diffusion_etm;
              }

              source = 0.;
              if (source_on) {
                source = wt_func * df->T[a][j] * d_area;
                source *= source_etm;
              }

              if (particle_momentum_on) {
                mass *= (1.0 - p_vol_frac);
                advection *= (1.0 - p_vol_frac);
              }

              /*lec->J[LEC_J_INDEX(peqn,pvar,ii,j)] += mass + advection + porous + diffusion +
               * source;*/
              J[j] += mass + advection + porous + diffusion + source;
              if (supg != 0) {
                J[j] += mass_supg + advection_supg + source_supg;
              }
            }
          }

          for (b = 0; b < MAX_MOMENTS; b++) {
            var = MOMENT0 + b;
            if (pdv[var]) {
              pvar = upd->vp[pg->imtrx][var];
              phi_j_vector = bf[var]->phi;
              J = &(lec->J[LEC_J_INDEX(peqn, pvar, ii, 0)]);

              for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                phi_j = phi_j_vector[j];

                mass = 0.;

                if (transient_run) {
                  if (mass_on) {
                    mass = d_rho->moment[b][j] * v_dot[a];
                    mass *= -wt_func * d_area;
                    mass *= mass_etm;
                  }

                  /* if porous flow is considered. KSC on 5/10/95 */
                  if (porous_brinkman_on) {
                    mass /= por;
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

                  advection *= -wt_func * d_rho->moment[b][j] * d_area;
                  advection *= advection_etm;

                  if (porous_brinkman_on) {
                    por2 = por * por;
                    advection /= por2;
                  }
                }

                porous = 0.;
                if (porous_brinkman_on) {
                  if (vis != 0.) {
                    porous = v[a] * (d_rho->moment[b][j] * sc * speed / sqrt(per));
                    porous *= -phi_i * d_area;
                    porous *= porous_brinkman_etm;
                  } else {
                    porous = 0.;
                  }
                }

                if (particle_momentum_on) {
                  mass *= (1.0 - p_vol_frac);
                  advection *= (1.0 - p_vol_frac);
                }

                /*lec->J[peqn][pvar][ii][j] += mass + advection + porous + diffusion + source;*/
                J[j] += mass + advection + porous;
              }
            }
          }
          /*
           * J_m_nn
           */

          if (pdv[BOND_EVOLUTION]) {
            var = BOND_EVOLUTION;
            pvar = upd->vp[pg->imtrx][var];
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              phi_j = bf[var]->phi[j];

              mass = 0;
              advection = 0;
              source = 0;
              if (supg != 0.) {
                dbl d_wt_func = 0;
                for (p = 0; p < dim; p++) {
                  d_wt_func += supg * supg_terms.d_tau_dnn[j] * v[p] * bfm->grad_phi[i][p];
                }
                if (transient_run) {
                  if (mass_on) {
                    mass = v_dot[a] * rho;
                    mass *= -d_wt_func * d_area;
                    mass *= mass_etm;
                  }

                  if (porous_brinkman_on) {
                    mass /= por;
                  }

                  if (particle_momentum_on) {
                    mass *= ompvf;
                  }
                }

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

                  if (upd->PSPG_advection_correction) {
                    advection -= pspg[0] * grad_v[0][a];
                    advection -= pspg[1] * grad_v[1][a];
                    if (WIM == 3)
                      advection -= pspg[2] * grad_v[2][a];
                  }

                  advection *= rho;
                  advection *= -d_wt_func * d_area;
                  advection *= advection_etm;

                  if (porous_brinkman_on) {
                    por2 = por * por;
                    advection /= por2;
                  }

                  if (particle_momentum_on) {
                    advection *= ompvf;
                  }
                }

                /*
                 * Source term...
                 */

                if (source_on) {
                  source += f[a];
                  source *= d_wt_func * d_area;
                  source *= source_etm;
                }
              }
              diffusion = 0.;
              if (diffusion_on) {
                for (p = 0; p < VIM; p++) {
                  for (q = 0; q < VIM; q++) {
                    diffusion += grad_phi_i_e_a[p][q] * d_Pi->nn[q][p][j];
                  }
                }
                diffusion *= -d_area;
                diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
              }

              lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += mass + advection + diffusion + source;
            }
          }

          /*
           * J_m_degrade
           */

          if (pdv[RESTIME]) {
            var = RESTIME;
            pvar = upd->vp[pg->imtrx][var];
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              phi_j = bf[var]->phi[j];

              diffusion = 0.;
              if (diffusion_on) {
                for (p = 0; p < VIM; p++) {
                  for (q = 0; q < VIM; q++) {
                    diffusion += grad_phi_i_e_a[p][q] * d_Pi->degrade[q][p][j];
                  }
                }
                diffusion *= -d_area;
                diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
              }

              lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += diffusion;
            }
          }

          /*
           * J_m_EDDY_NU
           */

          if (pdv[EDDY_NU]) {
            var = EDDY_NU;
            pvar = upd->vp[pg->imtrx][var];
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              phi_j = bf[var]->phi[j];

              mass = 0.;
              advection = 0.;
              source = 0.0;
              if (supg != 0.) {
                dbl d_wt_func = 0;
                for (p = 0; p < dim; p++) {
                  d_wt_func += supg * supg_terms.d_tau_dEDDY_NU[j] * v[p] * bfm->grad_phi[i][p];
                }
                if (transient_run) {
                  if (mass_on) {
                    mass = v_dot[a] * rho;
                    mass *= -d_wt_func * d_area;
                    mass *= mass_etm;
                  }

                  if (porous_brinkman_on) {
                    mass /= por;
                  }

                  if (particle_momentum_on) {
                    mass *= ompvf;
                  }
                }

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

                  if (upd->PSPG_advection_correction) {
                    advection -= pspg[0] * grad_v[0][a];
                    advection -= pspg[1] * grad_v[1][a];
                    if (WIM == 3)
                      advection -= pspg[2] * grad_v[2][a];
                  }

                  advection *= rho;
                  advection *= -d_wt_func * d_area;
                  advection *= advection_etm;

                  if (porous_brinkman_on) {
                    por2 = por * por;
                    advection /= por2;
                  }

                  if (particle_momentum_on) {
                    advection *= ompvf;
                  }
                }

                /*
                 * Source term...
                 */

                if (source_on) {
                  source += f[a];
                  source *= d_wt_func * d_area;
                  source *= source_etm;
                }
              }

              diffusion = 0.;
              if (diffusion_on) {
                for (p = 0; p < VIM; p++) {
                  for (q = 0; q < VIM; q++) {
                    diffusion += grad_phi_i_e_a[p][q] * d_Pi->eddy_nu[q][p][j];
                  }
                }
                diffusion *= -d_area;
                diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
              }

              lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += mass + advection + diffusion + source;
            }
          }

          if (pdv[TURB_K]) {
            var = TURB_K;
            pvar = upd->vp[pg->imtrx][var];
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              phi_j = bf[var]->phi[j];

              mass = 0.;
              advection = 0.;
              source = 0.0;
              if (supg != 0.) {
                dbl d_wt_func = 0;
                for (p = 0; p < dim; p++) {
                  d_wt_func += supg * supg_terms.d_tau_dturb_k[j] * v[p] * bfm->grad_phi[i][p];
                }
                if (transient_run) {
                  if (mass_on) {
                    mass = v_dot[a] * rho;
                    mass *= -d_wt_func * d_area;
                    mass *= mass_etm;
                  }

                  if (porous_brinkman_on) {
                    mass /= por;
                  }

                  if (particle_momentum_on) {
                    mass *= ompvf;
                  }
                }

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

                  if (upd->PSPG_advection_correction) {
                    advection -= pspg[0] * grad_v[0][a];
                    advection -= pspg[1] * grad_v[1][a];
                    if (WIM == 3)
                      advection -= pspg[2] * grad_v[2][a];
                  }

                  advection *= rho;
                  advection *= -d_wt_func * d_area;
                  advection *= advection_etm;

                  if (porous_brinkman_on) {
                    por2 = por * por;
                    advection /= por2;
                  }

                  if (particle_momentum_on) {
                    advection *= ompvf;
                  }
                }

                /*
                 * Source term...
                 */

                if (source_on) {
                  source += f[a];
                  source *= d_wt_func * d_area;
                  source *= source_etm;
                }
              }

              diffusion = 0.;
              if (diffusion_on) {
                for (p = 0; p < VIM; p++) {
                  for (q = 0; q < VIM; q++) {
                    diffusion += grad_phi_i_e_a[p][q] * d_Pi->turb_k[q][p][j];
                  }
                }
                diffusion *= -d_area;
                diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
              }

              lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += mass + advection + diffusion + source;
            }
          }
          if (pdv[TURB_OMEGA]) {
            var = TURB_OMEGA;
            pvar = upd->vp[pg->imtrx][var];
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              phi_j = bf[var]->phi[j];

              mass = 0.;
              advection = 0.;
              source = 0.0;
              if (supg != 0.) {
                dbl d_wt_func = 0;
                for (p = 0; p < dim; p++) {
                  d_wt_func += supg * supg_terms.d_tau_dturb_omega[j] * v[p] * bfm->grad_phi[i][p];
                }
                if (transient_run) {
                  if (mass_on) {
                    mass = v_dot[a] * rho;
                    mass *= -d_wt_func * d_area;
                    mass *= mass_etm;
                  }

                  if (porous_brinkman_on) {
                    mass /= por;
                  }

                  if (particle_momentum_on) {
                    mass *= ompvf;
                  }
                }

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

                  if (upd->PSPG_advection_correction) {
                    advection -= pspg[0] * grad_v[0][a];
                    advection -= pspg[1] * grad_v[1][a];
                    if (WIM == 3)
                      advection -= pspg[2] * grad_v[2][a];
                  }

                  advection *= rho;
                  advection *= -d_wt_func * d_area;
                  advection *= advection_etm;

                  if (porous_brinkman_on) {
                    por2 = por * por;
                    advection /= por2;
                  }

                  if (particle_momentum_on) {
                    advection *= ompvf;
                  }
                }

                /*
                 * Source term...
                 */

                if (source_on) {
                  source += f[a];
                  source *= d_wt_func * d_area;
                  source *= source_etm;
                }
              }

              diffusion = 0.;
              if (diffusion_on) {
                for (p = 0; p < VIM; p++) {
                  for (q = 0; q < VIM; q++) {
                    diffusion += grad_phi_i_e_a[p][q] * d_Pi->turb_omega[q][p][j];
                  }
                }
                diffusion *= -d_area;
                diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
              }

              lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += mass + advection + diffusion + source;
            }
          }

          /*
           * J_m_F
           */

          if (pdv[FILL]) {
            var = FILL;
            pvar = upd->vp[pg->imtrx][var];

            J = &(lec->J[LEC_J_INDEX(peqn, pvar, ii, 0)]);

            phi_j_vector = bf[var]->phi;

            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              phi_j = phi_j_vector[j];

              dbl mass_supg = 0;
              dbl advection_supg = 0;
              dbl source_supg = 0;
              if (supg != 0.) {
                dbl d_wt_func = 0;
                for (p = 0; p < dim; p++) {
                  d_wt_func += supg * supg_terms.d_tau_dF[j] * v[p] * bfm->grad_phi[i][p];
                }
                if (transient_run) {
                  if (mass_on) {
                    mass_supg = v_dot[a] * rho;
                    mass_supg *= -d_wt_func * d_area;
                    mass_supg *= mass_etm;
                  }

                  if (porous_brinkman_on) {
                    mass_supg /= por;
                  }

                  if (particle_momentum_on) {
                    mass_supg *= ompvf;
                  }
                }

                if (advection_on) {
#ifdef DO_NO_UNROLL
                  for (p = 0; p < WIM; p++) {
                    advection_supg += (v[p] - x_dot[p]) * grad_v[p][a];
                  }
#else
                  advection_supg += (v[0] - x_dot[0]) * grad_v[0][a];
                  advection_supg += (v[1] - x_dot[1]) * grad_v[1][a];
                  if (WIM == 3)
                    advection_supg += (v[2] - x_dot[2]) * grad_v[2][a];
#endif

                  if (upd->PSPG_advection_correction) {
                    advection_supg -= pspg[0] * grad_v[0][a];
                    advection_supg -= pspg[1] * grad_v[1][a];
                    if (WIM == 3)
                      advection_supg -= pspg[2] * grad_v[2][a];
                  }

                  advection_supg *= rho;
                  advection_supg *= -d_wt_func * d_area;
                  advection_supg *= advection_etm;

                  if (porous_brinkman_on) {
                    por2 = por * por;
                    advection_supg /= por2;
                  }

                  if (particle_momentum_on) {
                    advection_supg *= ompvf;
                  }
                }

                /*
                 * Source term...
                 */

                if (source_on) {
                  source_supg += f[a];
                  source_supg *= d_wt_func * d_area;
                  source_supg *= source_etm;
                }
              }
              mass = 0.;
              if (transient_run) {
                if (mass_on) {
                  mass = d_rho->F[j] * v_dot[a];
                  mass *= -wt_func * d_area;
                  mass *= mass_etm;
                }

                /* if porous flow is considered. KSC on 5/10/95 */
                if (porous_brinkman_on) {
                  if (vis != 0.) {
                    mass /= por;
                  } else if (mp->viscosity == 0.0) {
                    GOMA_EH(GOMA_ERROR, "Error in Porous Brinkman specs.  See PRS ");
                  }
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
                advection *= -wt_func * d_rho->F[j] * d_area;
                advection *= advection_etm;

                /* if porous flow is considered. KSC on 5/10/95 */
                if (porous_brinkman_on) {

                  if (vis != 0.) {
                    por2 = por * por;
                    advection /= por2;
                  }
                }
              }

              porous = 0.;
              if (porous_brinkman_on) {
                if (vis != 0) {
                  porous = v[a] * (d_rho->F[j] * sc * speed / sqrt(per));
                  porous += v[a] * (d_flow_mu->F[j] / per);
                  porous *= -phi_i * d_area;
                  porous *= porous_brinkman_etm;
                } else if (mp->FSIModel > 0) {
                  porous = 0.0;
                  porous -= LubAux->dv_avg_df[a][j] * phi_i;
                  porous *= wt * fv->sdet * h3;
                  porous *= porous_brinkman_etm;
                } else {
                  porous = 0.;
                }
              }

              diffusion = 0.;
              if (diffusion_on) {
#ifdef DO_NO_UNROLL
                for (p = 0; p < VIM; p++) {
                  for (q = 0; q < VIM; q++) {
                    diffusion += grad_phi_i_e_a[p][q] * d_Pi->F[q][p][j];
                  }
                }
#else
                diffusion += grad_phi_i_e_a[0][0] * d_Pi->F[0][0][j];
                diffusion += grad_phi_i_e_a[1][1] * d_Pi->F[1][1][j];
                diffusion += grad_phi_i_e_a[1][0] * d_Pi->F[0][1][j];
                diffusion += grad_phi_i_e_a[0][1] * d_Pi->F[1][0][j];

                if (VIM == 3) {
                  diffusion += grad_phi_i_e_a[2][0] * d_Pi->F[0][2][j];
                  diffusion += grad_phi_i_e_a[2][1] * d_Pi->F[1][2][j];
                  diffusion += grad_phi_i_e_a[2][2] * d_Pi->F[2][2][j];
                  diffusion += grad_phi_i_e_a[0][2] * d_Pi->F[2][0][j];
                  diffusion += grad_phi_i_e_a[1][2] * d_Pi->F[2][1][j];
                }
#endif
                diffusion *= -d_area;
                diffusion *= diffusion_etm;
              }

              source = 0.;

              /* Stay away from evil Hessians. */
              if (source_on) {
                source = phi_i * df->F[a][j] * d_area;
                source *= source_etm;
              }

              if (particle_momentum_on) {
                mass *= (1.0 - p_vol_frac);
                advection *= (1.0 - p_vol_frac);
              }
              J[j] += mass + advection + porous + diffusion + source;
              if (supg != 0) {
                J[j] += mass_supg + advection_supg + source_supg;
              }
              /*lec->J[LEC_J_INDEX(peqn,pvar,ii,j)] +=mass + advection + porous + diffusion +
               * source; */
            }
          }

          /*
           * J_m_v
           */
          for (b = 0; b < WIM; b++) {
            var = VELOCITY1 + b;
            if (pdv[var]) {
              pvar = upd->vp[pg->imtrx][var];

              J = &(lec->J[LEC_J_INDEX(peqn, pvar, ii, 0)]);

              phi_j_vector = bf[var]->phi;

              for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {

                phi_j = phi_j_vector[j];
                double d_wt_func = 0;
                if (supg != 0.) {
                  d_wt_func = supg * supg_terms.tau * phi_j * bfm->grad_phi[i][b];

                  for (p = 0; p < dim; p++) {
                    d_wt_func += supg * supg_terms.d_tau_dv[b][j] * v[p] * bfm->grad_phi[i][p];
                  }
                }

                mass = 0.;
                if (transient_run) {
                  if (mass_on && (a == b)) {
                    /*mass = (1.+2.*tt) * phi_j/dt * (double)delta(a,b); */
                    mass = (1. + 2. * tt) * phi_j / dt;
                    mass *= -wt_func * rho * d_area;
                    if (supg != 0) {
                      double mass_b = v_dot[a];
                      mass_b *= -d_wt_func * rho * d_area;
                      mass += mass_b;
                    }
                    mass *= mass_etm;
                  }

                  if (porous_brinkman_on) {
                    if (vis != 0.) {
                      mass /= por;
                    } else if (mp->viscosity == 0) {
                      GOMA_EH(GOMA_ERROR, "incorrect viscosity settings on porous_brinkman");
                    }
                  }

                  if (particle_momentum_on) {
                    mass *= ompvf;
                  }
                }

                porous = 0.0;
                if (porous_brinkman_on) {
                  if (vis != 0) {
                    porous = ((rho * sc / sqrt(per)) * (2. * v[b]) * v[a] +
                              (rho * sc * speed / sqrt(per) + vis / per) * (double)delta(a, b));
                    porous *= -phi_i * phi_j * d_area;
                    porous *= porous_brinkman_etm;
                  } else if (mp->viscosity != 0) {
                    porous = delta(a, b);
                    porous *= phi_i * phi_j * wt * fv->sdet * h3;
                    // porous *= phi_i*phi_j* wt * h3;
                    porous *= porous_brinkman_etm;
                  } else {
                    GOMA_EH(GOMA_ERROR, "incorrect viscosity settings on porous_brinkman");
                  }
                }

                advection = 0.;
                advection_a = 0.;
                if (advection_on) {
                  advection_a += phi_j * grad_v[b][a];

#ifdef DO_NO_UNROLL
                  for (p = 0; p < WIM; p++) {
                    advection_a += (v[p] - x_dot[p]) * bf[var]->grad_phi_e[j][b][p][a];
                  }
#else
                  advection_a += (v[0] - x_dot[0]) * bf[var]->grad_phi_e[j][b][0][a];
                  advection_a += (v[1] - x_dot[1]) * bf[var]->grad_phi_e[j][b][1][a];
                  if (WIM == 3)
                    advection_a += (v[2] - x_dot[2]) * bf[var]->grad_phi_e[j][b][2][a];
#endif

                  if (upd->PSPG_advection_correction) {
                    advection_a -= pspg[0] * bf[var]->grad_phi_e[j][b][0][a];
                    advection_a -= pspg[1] * bf[var]->grad_phi_e[j][b][1][a];
                    if (WIM == 3)
                      advection_a -= pspg[2] * bf[var]->grad_phi_e[j][b][2][a];
                    advection_a -= d_pspg.v[0][b][j] * grad_v[0][a];
                    advection_a -= d_pspg.v[1][b][j] * grad_v[1][a];
                    if (WIM == 3)
                      advection_a -= d_pspg.v[2][b][j] * grad_v[2][a];
                  }

                  advection_a *= rho * -wt_func * d_area;
                  advection_b = 0.;
                  if (supg != 0.) {
                    advection_b += (v[0] - x_dot[0]) * grad_v[0][a];
                    advection_b += (v[1] - x_dot[1]) * grad_v[1][a];
                    if (WIM == 3)
                      advection_b += (v[2] - x_dot[2]) * grad_v[2][a];
                    if (upd->PSPG_advection_correction) {
                      advection_b -= pspg[0] * grad_v[0][a];
                      advection_b -= pspg[1] * grad_v[1][a];
                      if (WIM == 3)
                        advection_b -= pspg[2] * grad_v[2][a];
                    }
                    advection_b *= -d_wt_func * rho * d_area;
                  }
                  advection = advection_a + advection_b;
                  advection *= advection_etm;

                  /* if porous flow is considered. KSC on 5/10/95 */
                  if (porous_brinkman_on) {
                    if (vis != 0) {
                      por2 = por * por;
                      advection /= por2;
                    }
                  }

                  if (particle_momentum_on) {
                    advection *= ompvf;
                  }
                }

                diffusion = 0.;
                if (diffusion_on) {
#ifdef DO_NO_UNROLL
                  for (p = 0; p < VIM; p++) {
                    for (q = 0; q < VIM; q++) {
                      diffusion += grad_phi_i_e_a[p][q] * d_Pi->v[q][p][b][j];
                    }
                  }
#else
                  diffusion += grad_phi_i_e_a[0][0] * d_Pi->v[0][0][b][j];
                  diffusion += grad_phi_i_e_a[1][1] * d_Pi->v[1][1][b][j];
                  diffusion += grad_phi_i_e_a[1][0] * d_Pi->v[0][1][b][j];
                  diffusion += grad_phi_i_e_a[0][1] * d_Pi->v[1][0][b][j];

                  if (VIM == 3) {
                    diffusion += grad_phi_i_e_a[2][0] * d_Pi->v[0][2][b][j];
                    diffusion += grad_phi_i_e_a[2][1] * d_Pi->v[1][2][b][j];
                    diffusion += grad_phi_i_e_a[0][2] * d_Pi->v[2][0][b][j];
                    diffusion += grad_phi_i_e_a[1][2] * d_Pi->v[2][1][b][j];
                    diffusion += grad_phi_i_e_a[2][2] * d_Pi->v[2][2][b][j];
                  }
#endif

                  diffusion *= -d_area;
                  diffusion *= diffusion_etm;
                }

                source = 0.;
                if (source_on) {
                  source = wt_func * df->v[a][b][j] * d_area;
                  source += d_wt_func * f[a] * d_area;
                  source *= source_etm;
                }

                continuity_stabilization = 0.0;
                if (Cont_GLS) {
                  for (p = 0; p < VIM; p++) {
                    continuity_stabilization += grad_phi_i_e_a[p][p];
                  }
                  continuity_stabilization *= d_cont_gls->v[b][j] * d_area;
                }

                J[j] += mass + advection + porous + diffusion + source + continuity_stabilization;

                /*lec->J[LEC_J_INDEX(peqn,pvar,ii,j)] +=  mass + advection + porous + diffusion +
                 * source; */
              }
            }
          }

          /*
           * J_m_vd
           */
          if (pdv[VORT_DIR1]) {
            for (b = 0; b < DIM; b++) {
              var = VORT_DIR1 + b;
              if (pdv[var]) {
                pvar = upd->vp[pg->imtrx][var];
                for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {

                  phi_j = bf[var]->phi[j];

                  diffusion = 0.;
                  if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {

                    for (p = 0; p < VIM; p++) {
                      for (q = 0; q < VIM; q++) {
                        diffusion += grad_phi_i_e_a[p][q] * d_Pi->vd[q][p][b][j];
                      }
                    }
                    diffusion *= -det_J * wt;
                    diffusion *= h3;
                    diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
                  }

                  lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += diffusion;
                }
              }
            }
          }

          /*
           * J_m_E (E_field from EHD source)
           */
          if (pdv[EFIELD1]) {
            for (b = 0; b < WIM; b++) {
              var = EFIELD1 + b;
              if (pdv[var]) {
                pvar = upd->vp[pg->imtrx][var];
                for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {

                  phi_j = bf[var]->phi[j];

                  source = 0.;
                  if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
                    source = phi_i * df->E[a][b][j] * det_J * h3 * wt;
                    source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
                  }

                  lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
                }
              }
            }
          }

          /*
           * J_m_lubp (lubrication pressure field for porous-brinkman)
           */
          if (pdv[LUBP]) {
            var = LUBP;
            if (pdv[var]) {
              pvar = upd->vp[pg->imtrx][var];
              for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {

                porous = 0.;
                if (porous_brinkman_on) {
                  /*This needs no insulation as the code shouldn't get here
                   *for traditional Brinkman form */

                  /* Need a few more basis functions */
                  dbl grad_phi_j[DIM], grad_II_phi_j[DIM], d_grad_II_phi_j_dmesh[DIM][DIM][MDE];
                  ShellBF(var, j, &phi_j, grad_phi_j, grad_II_phi_j, d_grad_II_phi_j_dmesh,
                          n_dof[MESH_DISPLACEMENT1], dof_map);

                  /* Assemble */
                  for (b = 0; b < DIM; b++) {
                    porous -= LubAux->dv_dgradp[a][b] * grad_II_phi_j[b];
                  }
                  porous -= LubAux->dv_avg_dp2[a] * phi_j;
                  porous *= phi_i * fv->sdet * wt * h3;
                  // porous *= phi_i * wt * h3;
                  porous *= porous_brinkman_etm;
                }

                lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += porous;
              }
            }
          }

          /*
           * J_m_d_delta_h (Delta_H equation for shell for porous-brinkman)
           */
          if (pdv[SHELL_DELTAH]) {
            var = SHELL_DELTAH;
            if (pdv[var]) {

              pvar = upd->vp[pg->imtrx][var];
              for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {

                phi_j = bf[var]->phi[j];

                porous = 0.;
                if (porous_brinkman_on) {
                  porous -= LubAux->dv_avg_ddh[a] * phi_j;
                  porous *= phi_i * fv->sdet * wt * h3;
                  porous *= porous_brinkman_etm;
                }

                lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += porous;
              }
            }
          }

          /*
           * J_m_d_rs (realsolid equation for shell)
           */
          if (upd->vp[pg->imtrx][SOLID_DISPLACEMENT1]) {
            for (b = 0; b < dim; b++) {
              var = SOLID_DISPLACEMENT1 + b;
              if (pdv[var]) {
                if (mp->FSIModel == FSI_REALSOLID_CONTINUUM) {
                  for (j = 0; j < ei[pg->imtrx]->dof[SOLID_DISPLACEMENT1]; j++) {
                    jk = dof_map[j];

                    phi_j = bf[var]->phi[j];

                    porous = 0.;
                    if (porous_brinkman_on) {
                      porous -= LubAux->dv_avg_drs[a][b][j] * fv->sdet * h3;
                    }
                    porous *= phi_i * wt;

                    J[jk] += porous;
                  }
                }
              }
            }
          }

          /*
           * J_m_d_curv (lubrication curvature equation for shell)
           */
          if (pdv[SHELL_LUB_CURV]) {
            var = SHELL_LUB_CURV;
            if (pdv[var]) {

              pvar = upd->vp[pg->imtrx][var];
              for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {

                phi_j = bf[var]->phi[j];

                porous = 0.;
                if (porous_brinkman_on) {
                  porous -= LubAux->dv_avg_dk[a] * phi_j;
                  porous *= phi_i * fv->sdet * wt * h3;
                  porous *= porous_brinkman_etm;
                }

                lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += porous;
              }
            }
          }

          /*
           * J_m_sh_fp (film pressure field for porous-brinkman)
           */
          if (pdv[SHELL_FILMP]) {
            var = SHELL_FILMP;
            if (pdv[var]) {
              pvar = upd->vp[pg->imtrx][var];
              for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {

                porous = 0.;
                if (porous_brinkman_on) {
                  /*This needs no insulation as the code shouldn't get here
                   *for traditional Brinkman form */

                  /* Need a few more basis functions */
                  dbl grad_phi_j[DIM], grad_II_phi_j[DIM], d_grad_II_phi_j_dmesh[DIM][DIM][MDE];
                  ShellBF(var, j, &phi_j, grad_phi_j, grad_II_phi_j, d_grad_II_phi_j_dmesh,
                          n_dof[MESH_DISPLACEMENT1], dof_map);

                  /* Assemble */
                  for (b = 0; b < DIM; b++) {
                    porous -= LubAux->dv_dgradp[a][b] * grad_II_phi_j[b];
                  }
                  porous *= phi_i * fv->sdet * wt * h3;
                  porous *= porous_brinkman_etm;
                }

                lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += porous;
              }
            }
          }

          /*
           * J_m_sh_fh (film thickness for porous-brinkman)
           */
          if (pdv[SHELL_FILMH]) {
            var = SHELL_FILMH;
            if (pdv[var]) {
              pvar = upd->vp[pg->imtrx][var];
              for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                porous = 0.;
                if (porous_brinkman_on) {
                  /*This needs no insulation as the code shouldn't get here
                   *for traditional Brinkman form */

                  dbl grad_phi_j[DIM], grad_II_phi_j[DIM], d_grad_II_phi_j_dmesh[DIM][DIM][MDE];
                  ShellBF(var, j, &phi_j, grad_phi_j, grad_II_phi_j, d_grad_II_phi_j_dmesh,
                          n_dof[MESH_DISPLACEMENT1], dof_map);

                  /* Assemble */
                  porous -= LubAux->dv_avg_dh1[a][j] * grad_II_phi_j[a];
                  porous -= LubAux->dv_avg_dh2[a][j] * phi_j;
                  porous *= phi_i * fv->sdet * wt * h3;
                  porous *= porous_brinkman_etm;
                }

                lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += porous;
              }
            }
          }

          /*
           * J_m_c
           */

          if (pdv[MASS_FRACTION]) {
            var = MASS_FRACTION;
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              phi_j = bf[var]->phi[j];

              for (w = 0; w < pd->Num_Species_Eqn; w++) {
                mass = 0.;
                if (transient_run) {
                  if (pd->e[pg->imtrx][eqn] & T_MASS) {
                    /* For particle momentum stuff, d_rho->* = 0, since the density per phase
                     * is actually constant.  However, the mass term was multiplied by (1-phi),
                     * so that we get a J_m_c entry now.
                     */
                    if (particle_momentum_on && w == species)
                      mass += -rho * v_dot[a];
                    else
                      mass += d_rho->C[w][j] * v_dot[a];
                    mass *= -phi_i * det_J * wt * h3;
                    mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
                  }

                  if (pd->e[pg->imtrx][eqn] &
                      T_POROUS_BRINK) /* if porous flow is considered. KSC on 5/10/95 */
                  {
                    mass /= por;
                  }
                }

                advection = 0.;
                if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
                  for (p = 0; p < WIM; p++) {
                    advection += (v[p] - x_dot[p]) * grad_v[p][a];
                  }
                  if (particle_momentum_on && w == species)
                    advection *= -phi_j * (-wt_func * rho * det_J * wt * h3);
                  else
                    advection *= -wt_func * d_rho->C[w][j] * det_J * wt * h3;
                  advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];

                  if (pd->e[pg->imtrx][eqn] &
                      T_POROUS_BRINK) /* if porous flow is considered. KSC on 5/10/95 */
                  {
                    por2 = por * por;
                    advection /= por2;
                  }
                }

                porous = 0.;
                if (pd->e[pg->imtrx][eqn] & T_POROUS_BRINK) {
                  if (vis != 0.) {
                    porous = v[a] * (d_rho->C[w][j] * sc * speed / sqrt(per) -
                                     0.5 * rho * sc * speed / (per * sqrt(per)) * d_per_dc[w][j] -
                                     vis / (per * per) * d_per_dc[w][j] + d_flow_mu->C[w][j] / per);
                    porous *= -phi_i * det_J * wt * h3;
                    porous *= pd->etm[pg->imtrx][eqn][(LOG2_POROUS_BRINK)];
                  }
                }

                diffusion = 0.;
                if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {
                  for (p = 0; p < VIM; p++) {
                    for (q = 0; q < VIM; q++) {
                      diffusion += grad_phi_i_e_a[p][q] * d_Pi->C[q][p][w][j];
                    }
                  }
                  diffusion *= -det_J * wt;
                  diffusion *= h3;
                  diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
                }

                source = 0.;
                if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
                  /* df->C was calculated in mm_std_models.c */
                  source = phi_i * df->C[a][w][j] * det_J * h3 * wt;
                  source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
                }

                lec->J[LEC_J_INDEX(peqn, MAX_PROB_VAR + w, ii, j)] +=
                    mass + advection + porous + diffusion + source;
              }
            }
          }

          /*
           * J_m_P
           */

          if (pdv[PRESSURE]) {
            var = PRESSURE;
            pvar = upd->vp[pg->imtrx][var];

            J = &(lec->J[LEC_J_INDEX(peqn, pvar, ii, 0)]);

            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {

              /*  phi_j = bf[var]->phi[j]; */

              advection = 0;

              if (upd->PSPG_advection_correction) {
                advection -= d_pspg.P[0][j] * grad_v[0][a];
                advection -= d_pspg.P[1][j] * grad_v[1][a];
                if (WIM == 3)
                  advection -= d_pspg.P[2][j] * grad_v[2][a];
              }
              advection *= rho;
              advection *= -wt_func * d_area;
              advection *= advection_etm;

              diffusion = 0.;

              if (diffusion_on) {
#ifdef DO_NO_UNROLL
                for (p = 0; p < VIM; p++) {
                  for (q = 0; q < VIM; q++) {

                    diffusion -= grad_phi_i_e_a[p][q] * d_Pi->P[q][p][j];

#ifdef DEBUG_MOMENTUM_JAC
                    if (ei[pg->imtrx]->ielem == 0)
                      printf("P diffusion, %d, %d  -= %g %g\n", p, q, grad_phi_i_e_a[p][q],
                             d_Pi->P[j][q][p]);
#endif /* DEBUG_MOMENTUM_JAC */
                  }
                }
#else
                diffusion -= grad_phi_i_e_a[0][0] * d_Pi->P[0][0][j];
                diffusion -= grad_phi_i_e_a[1][1] * d_Pi->P[1][1][j];
                diffusion -= grad_phi_i_e_a[1][0] * d_Pi->P[0][1][j];
                diffusion -= grad_phi_i_e_a[0][1] * d_Pi->P[1][0][j];

                if (VIM == 3) {
                  diffusion -= grad_phi_i_e_a[2][2] * d_Pi->P[2][2][j];
                  diffusion -= grad_phi_i_e_a[2][1] * d_Pi->P[1][2][j];
                  diffusion -= grad_phi_i_e_a[2][0] * d_Pi->P[0][2][j];
                  diffusion -= grad_phi_i_e_a[0][2] * d_Pi->P[2][0][j];
                  diffusion -= grad_phi_i_e_a[1][2] * d_Pi->P[2][1][j];
                }
#endif
                diffusion *= d_area;
                diffusion *= diffusion_etm;
              }

#ifdef DEBUG_MOMENTUM_JAC
              if (ei[pg->imtrx]->ielem == 0)
                fprintf(stdout, "\tJ_m_P[%d][%d][%d] += %10f %10f %10f %10f\n", a, i, j, mass,
                        advection, porous, diffusion);
#endif /* DEBUG_MOMENTUM_JAC */

              /*lec->J[LEC_J_INDEX(peqn,pvar,ii,j)] += diffusion ;  */
              J[j] += advection + diffusion;
            }
          }

          /*
           * J_m_S
           */
          if (pdv[POLYMER_STRESS11]) {
            for (mode = 0; mode < vn->modes; mode++) {
              for (b = 0; b < VIM; b++) {
                for (c = 0; c < VIM; c++) {
                  var = v_s[mode][b][c];

                  if (pdv[var]) {

                    pvar = upd->vp[pg->imtrx][var];

                    for (j = 0; j < ei[pg->imtrx]->dof[var]; j++)

                    {
                      phi_j = bf[var]->phi[j];

                      diffusion = 0.;

                      if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {
                        for (p = 0; p < VIM; p++) {
                          for (q = 0; q < VIM; q++) {
                            diffusion -= grad_phi_i_e_a[p][q] * d_Pi->S[p][q][mode][b][c][j];
                          }
                        }
                        diffusion *= det_J * wt * h3;
                        diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
                      }

                      lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += diffusion;
                    }
                  }
                }
              }
            }
          }

          /*
           * J_m_G
           */
          if (gn->ConstitutiveEquation == BINGHAM_MIXED ||
              (pd->gv[POLYMER_STRESS11] &&
               (vn->evssModel == EVSS_F || vn->evssModel == LOG_CONF ||
                vn->evssModel == SQRT_CONF || vn->evssModel == CONF ||
                vn->evssModel == EVSS_GRADV || vn->evssModel == LOG_CONF_GRADV))) {
            for (b = 0; b < VIM; b++) {
              for (c = 0; c < VIM; c++) {
                var = v_g[b][c];

                if (pdv[var]) {

                  pvar = upd->vp[pg->imtrx][var];

                  for (j = 0; j < ei[pg->imtrx]->dof[var]; j++)

                  {
                    phi_j = bf[var]->phi[j];

                    diffusion = 0.;

                    if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {
                      for (p = 0; p < VIM; p++) {
                        for (q = 0; q < VIM; q++) {
                          diffusion += grad_phi_i_e_a[p][q] * d_Pi->g[q][p][b][c][j];
                        }
                      }
                      diffusion *= det_J * wt * h3;
                      diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
                    }

                    lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += diffusion;
                  }
                }
              }
            }
          }

          /*
           * J_m_d
           */
          if (pdv[MESH_DISPLACEMENT1]) {
            for (b = 0; b < dim; b++) {
              var = MESH_DISPLACEMENT1 + b;
              if (pdv[var]) {
                pvar = upd->vp[pg->imtrx][var];

                J = &(lec->J[LEC_J_INDEX(peqn, pvar, ii, 0)]);

                for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                  phi_j = bf[var]->phi[j];

                  d_det_J_dmesh_bj = bfm->d_det_J_dm[b][j];

                  dh3dmesh_bj = fv->dh3dq[b] * phi_j;

                  mass = 0.;

                  if (transient_run) {
                    if (mass_on) {
                      mass = v_dot[a];
                      mass *= -phi_i * rho * (d_det_J_dmesh_bj * h3 + det_J * dh3dmesh_bj) * wt;
                      mass *= mass_etm;
                    }

                    /* if porous flow is considered. KSC on 5/10/95 */
                    if (porous_brinkman_on) {
                      mass /= por;
                    }
                  }

                  porous = 0.;
                  if (porous_brinkman_on) {
                    if (vis != 0.) {
                      porous += v[a] * (rho * sc * speed / sqrt(per) + vis / per);
                      porous *= -phi_i * wt * (d_det_J_dmesh_bj * h3 + det_J * dh3dmesh_bj);
                      porous *= porous_brinkman_etm;
                    }

                    /* Case for sheet-only  */
                    else if (mp->FSIModel == FSI_SHELL_ONLY) {
                      jk = dof_map[j];
                      porous += (v[a] - LubAux->v_avg[a]) * fv->dsurfdet_dx[b][jk] * h3;
                      porous += (v[a] - LubAux->v_avg[a]) * fv->sdet * dh3dmesh_bj;
                      porous -= LubAux->dv_avg_dx[a][b][j] * fv->sdet * h3;
                      porous *= phi_i * wt;
                    }

                    else if (vis == 0. && mp->viscosity == 0.) {
                      GOMA_EH(GOMA_ERROR, "cannot have both flowing liquid viscosity and "
                                          "mp->viscosity equal to zero");
                    }
                  }

                  advection = 0.;

                  if (advection_on) {
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
#ifdef DO_NO_UNROLL
                    advection_a = 0.;
                    for (p = 0; p < WIM; p++) {
                      advection_a += (v[p] - x_dot[p]) * fv->d_grad_v_dmesh[p][a][b][j];
                    }
                    advection_a *= -wt_func * rho * d_area;

                    advection_b = 0.;
                    for (p = 0; p < WIM; p++) {
                      advection_b += (v[p] - x_dot[p]) * grad_v[p][a];
                    }
#else
                    advection_a = 0.;
                    advection_a += (v[0] - x_dot[0]) * fv->d_grad_v_dmesh[0][a][b][j];
                    advection_a += (v[1] - x_dot[1]) * fv->d_grad_v_dmesh[1][a][b][j];
                    if (WIM == 3)
                      advection_a += (v[2] - x_dot[2]) * fv->d_grad_v_dmesh[2][a][b][j];

                    advection_a *= -wt_func * rho * d_area;

                    advection_b = 0.;
                    advection_b += (v[0] - x_dot[0]) * grad_v[0][a];
                    advection_b += (v[1] - x_dot[1]) * grad_v[1][a];
                    if (WIM == 3)
                      advection_b += (v[2] - x_dot[2]) * grad_v[2][a];
#endif
                    advection_b *=
                        -wt_func * rho * wt * (d_det_J_dmesh_bj * h3 + det_J * dh3dmesh_bj);

                    advection_c = 0.;
                    if (transient_run) {
                      if (mass_on) {
#ifdef DO_NO_UNROLL
                        for (p = 0; p < WIM; p++) {
                          advection_c +=
                              (-(1. + 2. * tt) * phi_j / dt * (double)delta(p, b)) * grad_v[p][a];
                        }
#else
                        advection_c +=
                            (-(1. + 2. * tt) * phi_j / dt * (double)delta(0, b)) * grad_v[0][a];
                        advection_c +=
                            (-(1. + 2. * tt) * phi_j / dt * (double)delta(1, b)) * grad_v[1][a];
                        if (WIM == 3)
                          advection_c +=
                              (-(1. + 2. * tt) * phi_j / dt * (double)delta(2, b)) * grad_v[2][a];
#endif
                        advection_c *= -wt_func * rho * d_area;
                      }
                    }

                    advection = advection_a + advection_b + advection_c;

                    advection *= advection_etm;
                    /* if porous flow is considered. KSC on 5/10/95 */
                    if (porous_brinkman_on) {
                      por2 = por * por;
                      advection /= por2;
                    }
                  }

                  /*
                   * Diffusion...
                   */

                  diffusion = 0.;

                  if (diffusion_on) {

                    /* Three parts:
                     *   diff_a =
                     *   Int ( d(grad(phi_i e_a))/dmesh : Pi h3 |Jv|)
                     *
                     *   diff_b =
                     *   Int ( grad(phi_i e_a) : d(Pi)/dmesh h3 |Jv|)
                     *
                     *   diff_c =
                     *   Int ( grad(phi_i e_a) : Pi d(h3|Jv|)/dmesh )
                     */

                    diff_a = 0.;
                    diff_b = 0.;
                    diff_c = 0.;
#ifdef DO_NO_UNROLL
                    for (p = 0; p < VIM; p++) {
                      for (q = 0; q < VIM; q++) {
                        diff_a += d_grad_phi_i_e_a_dmesh[p][q][b][j] * Pi[q][p];

                        diff_b += grad_phi_i_e_a[p][q] * d_Pi->X[q][p][b][j];

                        diff_c += grad_phi_i_e_a[p][q] * Pi[q][p];
                      }
                    }
#else
                    diff_a += d_grad_phi_i_e_a_dmesh[0][0][b][j] * Pi[0][0];
                    diff_a += d_grad_phi_i_e_a_dmesh[1][1][b][j] * Pi[1][1];
                    diff_a += d_grad_phi_i_e_a_dmesh[0][1][b][j] * Pi[1][0];
                    diff_a += d_grad_phi_i_e_a_dmesh[1][0][b][j] * Pi[0][1];

                    diff_b += grad_phi_i_e_a[0][0] * d_Pi->X[0][0][b][j];
                    diff_b += grad_phi_i_e_a[1][1] * d_Pi->X[1][1][b][j];
                    diff_b += grad_phi_i_e_a[0][1] * d_Pi->X[1][0][b][j];
                    diff_b += grad_phi_i_e_a[1][0] * d_Pi->X[0][1][b][j];

                    diff_c += grad_phi_i_e_a[0][0] * Pi[0][0];
                    diff_c += grad_phi_i_e_a[1][1] * Pi[1][1];
                    diff_c += grad_phi_i_e_a[0][1] * Pi[1][0];
                    diff_c += grad_phi_i_e_a[1][0] * Pi[0][1];

                    if (VIM == 3) {
                      diff_a += d_grad_phi_i_e_a_dmesh[2][2][b][j] * Pi[2][2];
                      diff_a += d_grad_phi_i_e_a_dmesh[2][1][b][j] * Pi[1][2];
                      diff_a += d_grad_phi_i_e_a_dmesh[2][0][b][j] * Pi[0][2];
                      diff_a += d_grad_phi_i_e_a_dmesh[1][2][b][j] * Pi[2][1];
                      diff_a += d_grad_phi_i_e_a_dmesh[0][2][b][j] * Pi[2][0];

                      diff_b += grad_phi_i_e_a[2][2] * d_Pi->X[2][2][b][j];
                      diff_b += grad_phi_i_e_a[2][1] * d_Pi->X[1][2][b][j];
                      diff_b += grad_phi_i_e_a[2][0] * d_Pi->X[0][2][b][j];
                      diff_b += grad_phi_i_e_a[1][2] * d_Pi->X[2][1][b][j];
                      diff_b += grad_phi_i_e_a[0][2] * d_Pi->X[2][0][b][j];

                      diff_c += grad_phi_i_e_a[2][2] * Pi[2][2];
                      diff_c += grad_phi_i_e_a[2][1] * Pi[1][2];
                      diff_c += grad_phi_i_e_a[2][0] * Pi[0][2];
                      diff_c += grad_phi_i_e_a[1][2] * Pi[2][1];
                      diff_c += grad_phi_i_e_a[0][2] * Pi[2][0];
                    }

#endif
                    diff_a *= -d_area;
                    diff_b *= -d_area;
                    diff_c *= -wt * (d_det_J_dmesh_bj * h3 +

                                     det_J * dh3dmesh_bj);
                    diffusion = diff_a + diff_b + diff_c;

                    diffusion *= diffusion_etm;
                  }

                  /*
                   * Source term...
                   */

                  source = 0.;
                  if (source_on) {
                    source += phi_i * wt *
                              (f[a] * d_det_J_dmesh_bj * h3 + f[a] * det_J * dh3dmesh_bj +
                               df->X[a][b][j] * det_J * h3);
                    source *= source_etm;
                  }

                  if (particle_momentum_on) {
                    mass *= (1.0 - p_vol_frac);
                    advection *= (1.0 - p_vol_frac);
                  }

                  continuity_stabilization = 0.0;
                  if (Cont_GLS) {
                    for (p = 0; p < VIM; p++) {
                      continuity_stabilization +=
                          d_grad_phi_i_e_a_dmesh[p][p][b][j] * cont_gls * d_area;
                      continuity_stabilization +=
                          grad_phi_i_e_a[p][p] * d_cont_gls->X[b][j] * d_area;
                      continuity_stabilization += grad_phi_i_e_a[p][p] * cont_gls *
                                                  (d_det_J_dmesh_bj * h3 + det_J * dh3dmesh_bj) *
                                                  wt;
                    }
                  }

                  J[j] += mass + advection + porous + diffusion + source + continuity_stabilization;
                }

                /* Special Case for shell lubrication velocities with bounding continuum.  Note
                 * j-loop is over bulk element */
                if (mp->FSIModel == FSI_MESH_CONTINUUM || mp->FSIModel == FSI_REALSOLID_CONTINUUM) {
                  for (j = 0; j < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; j++) {
                    jk = dof_map[j];

                    phi_j = bf[var]->phi[j];
                    dh3dmesh_bj = fv->dh3dq[b] * phi_j;

                    porous = 0.;
                    if (porous_brinkman_on) {
                      porous += (v[a] - LubAux->v_avg[a]) * fv->dsurfdet_dx[b][jk] * h3;
                      porous += (v[a] - LubAux->v_avg[a]) * fv->sdet * dh3dmesh_bj;
                      porous -= LubAux->dv_avg_dx[a][b][j] * fv->sdet * h3;
                    }
                    porous *= phi_i * wt;

                    J[jk] += porous;
                  }
                }
              }
            }
          }
        } /* end of if(active_dofs) */
      }   /* end of for(i=ei[pg->imtrx]->dof*/
    }
  }
  safe_free((void *)n_dof);
  return (status);
}

void ve_polymer_stress(double gamma[DIM][DIM],
                       double stress[DIM][DIM],
                       STRESS_DEPENDENCE_STRUCT *d_stress) {
  for (int i = 0; i < VIM; i++) {
    for (int j = 0; j < VIM; j++) {
      stress[i][j] = 0;
    }
  }
  switch (vn->evssModel) {
  case LOG_CONF:
  case LOG_CONF_GRADV:
  case LOG_CONF_TRANSIENT_GRADV:
  case LOG_CONF_TRANSIENT: {
    dbl exp_s[DIM][DIM];
    dbl R1[DIM][DIM];
    dbl eig_values[DIM];

    for (int mode = 0; mode < vn->modes; mode++) {
      dbl mup = viscosity(ve[mode]->gn, gamma, NULL);
      dbl lambda = polymer_time_const(ve[mode]->time_const_st, gamma, NULL);
      if (vn->evssModel == LOG_CONF_TRANSIENT || vn->evssModel == LOG_CONF_TRANSIENT_GRADV) {
#ifdef ANALEIG_PLEASE
        analytical_exp_s(fv->S[mode], exp_s[mode], eig_values, R1, d_exp_s_ds[mode]);
#else

        if (pg->imtrx == upd->matrix_index[POLYMER_STRESS11]) {
          compute_exp_s(fv_old->S[mode], exp_s, eig_values, R1);
        } else {
          compute_exp_s(fv->S[mode], exp_s, eig_values, R1);
        }
#endif
      } else {
        compute_exp_s(fv->S[mode], exp_s, eig_values, R1);
      }

      for (int i = 0; i < VIM; i++) {
        for (int j = 0; j < VIM; j++) {
          stress[i][j] += mup / lambda * (exp_s[i][j] - (double)delta(i, j));
        }
      }
    }
  } break;
  case SQRT_CONF: {
    for (int mode = 0; mode < vn->modes; mode++) {
      /* get polymer viscosity */
      dbl mup = viscosity(ve[mode]->gn, gamma, NULL);
      dbl lambda = polymer_time_const(ve[mode]->time_const_st, gamma, NULL);

      dbl bdotb[DIM][DIM];
      dbl b[DIM][DIM];
      for (int ii = 0; ii < VIM; ii++) {
        for (int jj = 0; jj < VIM; jj++) {
          if (ii <= jj) {
            b[ii][jj] = fv->S[mode][ii][jj];
            b[jj][ii] = b[ii][jj];
          }
        }
      }

      tensor_dot(b, b, bdotb, VIM);

      for (int ii = 0; ii < VIM; ii++) {
        for (int jj = 0; jj < VIM; jj++) {
          stress[ii][jj] += -(mup / lambda) * (delta(ii, jj) - bdotb[ii][jj]);
        }
      }
    }

  } break;
  case CONF: {
    for (int mode = 0; mode < vn->modes; mode++) {
      dbl k = 1.0;
      switch (vn->ConstitutiveEquation) {
      case ROLIE_POLY_FE: {

        dbl trace = 0;
        for (int i = 0; i < VIM; i++) {
          for (int j = 0; j < VIM; j++) {
            trace += fv->S[mode][i][j] * fv->S[mode][i][j];
          }
        }
        dbl lambda_max = ve[mode]->maximum_stretch_ratio;

        dbl lambda_s = sqrt(trace / 3);

        if (vn->ConstitutiveEquation == ROLIE_POLY_FE) {
          k = ((3 - lambda_s * lambda_s / (lambda_max * lambda_max)) *
               (1 - 1 / (lambda_max * lambda_max))) /
              (1 - lambda_s * lambda_s / (lambda_max * lambda_max) *
                       (3 - 1 / (lambda_max * lambda_max)));
        }
      } break;
      default:
        break;
      }
      dbl mup = viscosity(ve[mode]->gn, gamma, NULL);
      dbl lambda = polymer_time_const(ve[mode]->time_const_st, gamma, NULL);
      for (int i = 0; i < VIM; i++) {
        for (int j = 0; j < VIM; j++) {
          stress[i][j] += mup * k / lambda * (fv->S[mode][i][j] - (double)delta(i, j));
        }
      }
    }
  } break;
  default: // Regular stress formulations
  {
    for (int mode = 0; mode < vn->modes; mode++) {
      for (int i = 0; i < VIM; i++) {
        for (int j = 0; j < VIM; j++) {
          stress[i][j] += fv->S[mode][i][j];
        }
      }
    }
  } break;
  }
  if (d_stress != NULL) {
    int v_s[MAX_MODES][DIM][DIM];
    (void)stress_eqn_pointer(v_s);
    switch (vn->evssModel) {
    case LOG_CONF:
    case LOG_CONF_GRADV: {
      dbl exp_s[DIM][DIM];
      dbl d_exp_s_ds[DIM][DIM][DIM][DIM];

      for (int mode = 0; mode < vn->modes; mode++) {
        VISCOSITY_DEPENDENCE_STRUCT d_mup_struct;
        VISCOSITY_DEPENDENCE_STRUCT *d_mup = &d_mup_struct;
        POLYMER_TIME_CONST_DEPENDENCE_STRUCT d_lam_struct;
        POLYMER_TIME_CONST_DEPENDENCE_STRUCT *d_lam = &d_lam_struct;
        dbl mup = viscosity(ve[mode]->gn, gamma, d_mup);
        dbl lambda = polymer_time_const(ve[mode]->time_const_st, gamma, d_lam);
        compute_d_exp_s_ds(fv->S[mode], exp_s, d_exp_s_ds);

        for (int p = 0; p < VIM; p++) {
          for (int q = 0; q < VIM; q++) {
            for (int b = 0; b < VIM; b++) {
              for (int c = 0; c < VIM; c++) {
                int var = v_s[mode][b][c];
                if (c >= b) {
                  for (int j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                    d_stress->S[p][q][mode][b][c][j] =
                        mup / lambda * d_exp_s_ds[p][q][b][c] * bf[var]->phi[j];
                  }
                } else {
                  for (int j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                    d_stress->S[p][q][mode][b][c][j] = 0;
                  }
                }
              }
            }
          }
        }

        for (int p = 0; p < VIM; p++) {
          for (int q = 0; q < VIM; q++) {
            for (int b = 0; b < WIM; b++) {
              int var = VELOCITY1 + b;
              for (int j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                d_stress->v[p][q][b][j] += d_mup->v[b][j] / lambda * exp_s[p][q] -
                                           mup * d_lam->v[b][j] / (lambda * lambda) * exp_s[p][q];
              }
            }
          }
        }
      }
    } break;
    case LOG_CONF_TRANSIENT_GRADV:
    case LOG_CONF_TRANSIENT: {
      for (int mode = 0; mode < vn->modes; mode++) {
        for (int p = 0; p < VIM; p++) {
          for (int q = 0; q < VIM; q++) {
            for (int r = 0; r < VIM; r++) {
              for (int c = 0; c < VIM; c++) {
                int var = v_s[mode][r][c];
                for (int j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                  d_stress->S[p][q][mode][r][c][j] = 0;
                }
              }
            }
          }
        }
      }

    } break;
    case SQRT_CONF: {
      VISCOSITY_DEPENDENCE_STRUCT d_mup_struct;
      VISCOSITY_DEPENDENCE_STRUCT *d_mup = &d_mup_struct;
      POLYMER_TIME_CONST_DEPENDENCE_STRUCT d_lam_struct;
      POLYMER_TIME_CONST_DEPENDENCE_STRUCT *d_lam = &d_lam_struct;
      for (int mode = 0; mode < vn->modes; mode++) {
        dbl mup = viscosity(ve[mode]->gn, gamma, d_mup);
        dbl lambda = polymer_time_const(ve[mode]->time_const_st, gamma, d_lam);

        dbl b[DIM][DIM];
        for (int ii = 0; ii < VIM; ii++) {
          for (int jj = 0; jj < VIM; jj++) {
            if (ii <= jj) {
              b[ii][jj] = fv->S[mode][ii][jj];
              b[jj][ii] = b[ii][jj];
            }
          }
        }
        dbl bdotb[DIM][DIM];
        tensor_dot(b, b, bdotb, VIM);

        for (int p = 0; p < VIM; p++) {
          for (int q = 0; q < VIM; q++) {
            for (int r = 0; r < VIM; r++) {
              for (int c = 0; c < VIM; c++) {
                int var = v_s[mode][r][c];
                dbl db[DIM][DIM] = {{0.}};
                db[r][c] = 1.0;
                db[c][r] = 1.0;

                dbl dbdotb[DIM][DIM];
                dbl bdotdb[DIM][DIM];

                tensor_dot(b, db, bdotdb, VIM);
                tensor_dot(db, b, dbdotb, VIM);

                dbl conf[DIM][DIM];
                for (int ii = 0; ii < VIM; ii++) {
                  for (int jj = 0; jj < VIM; jj++) {
                    conf[ii][jj] = (mup / lambda) * (b[c][jj] * (double)delta(ii, r) +
                                                     b[ii][r] * (double)delta(jj, c));
                  }
                }
                for (int j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                  d_stress->S[p][q][mode][r][c][j] = conf[p][q] * bf[var]->phi[j];
                }
              }
            }
          }
        }

        for (int p = 0; p < VIM; p++) {
          for (int q = 0; q < VIM; q++) {
            for (int b = 0; b < WIM; b++) {
              int var = VELOCITY1 + b;
              for (int j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                d_stress->v[p][q][b][j] +=
                    -(d_mup->v[b][j] / lambda) * (delta(p, q) - bdotb[p][q]) +
                    (d_lam->v[b][j] * mup / (lambda * lambda)) * (delta(p, q) - bdotb[p][q]);
              }
            }
          }
        }
      }
    } break;
    case CONF: {
      for (int mode = 0; mode < vn->modes; mode++) {
        dbl k = 1.0;
        dbl d_k[DIM][DIM] = {{0.}};
        switch (vn->ConstitutiveEquation) {
        case ROLIE_POLY_FE: {

          dbl trace = 0;
          for (int i = 0; i < VIM; i++) {
            for (int j = 0; j < VIM; j++) {
              trace += fv->S[mode][i][j] * fv->S[mode][i][j];
            }
          }
          dbl lambda_max = ve[mode]->maximum_stretch_ratio;
          dbl d_trace_dc[DIM][DIM] = {{0.0}};
          for (int p = 0; p < VIM; p++) {
            for (int q = 0; q < VIM; q++) {
              d_trace_dc[p][q] = 0.0;
              for (int i = 0; i < VIM; i++) {
                for (int j = 0; j < VIM; j++) {
                  d_trace_dc[p][q] += 2.0 * fv->S[mode][i][j] *
                                      (delta(p, i) * delta(q, j) | delta(p, j) * delta(q, i));
                }
              }
            }
          }

          dbl lambda_s = sqrt(trace / 3);

          k = ((3 - lambda_s * lambda_s / (lambda_max * lambda_max)) *
               (1 - 1 / (lambda_max * lambda_max))) /
              (1 - lambda_s * lambda_s / (lambda_max * lambda_max) *
                       (3 - 1 / (lambda_max * lambda_max)));
          dbl d_k_dlamda_s = 4 * lambda_s * lambda_max * lambda_max *
                             (lambda_max * lambda_max - 1) /
                             (pow(lambda_s * lambda_s - lambda_max * lambda_max, 2.0) *
                              (3 * lambda_max * lambda_max - 1));

          for (int p = 0; p < VIM; p++) {
            for (int q = 0; q < VIM; q++) {
              d_k[p][q] = d_k_dlamda_s * d_trace_dc[p][q];
            }
          }
        } break;
        default:
          break;
        }
        VISCOSITY_DEPENDENCE_STRUCT d_mup_struct;
        VISCOSITY_DEPENDENCE_STRUCT *d_mup = &d_mup_struct;
        POLYMER_TIME_CONST_DEPENDENCE_STRUCT d_lam_struct;
        POLYMER_TIME_CONST_DEPENDENCE_STRUCT *d_lam = &d_lam_struct;
        dbl mup = viscosity(ve[mode]->gn, gamma, d_mup);
        dbl lambda = polymer_time_const(ve[mode]->time_const_st, gamma, d_lam);
        for (int a = 0; a < VIM; a++) {
          for (int b = 0; b < VIM; b++) {
            for (int p = 0; p < VIM; p++) {
              for (int q = 0; q < VIM; q++) {
                int var = v_s[mode][p][q];
                for (int j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                  d_stress->S[a][b][mode][p][q][j] +=
                      mup * k / lambda * delta(a, p) * delta(b, q) * bf[var]->phi[j] +
                      (mup * d_k[p][q] / lambda) * bf[var]->phi[j] * (fv->S[mode][a][b]);
                }
              }
            }
          }
        }
        for (int p = 0; p < VIM; p++) {
          for (int q = 0; q < VIM; q++) {
            for (int b = 0; b < WIM; b++) {
              int var = VELOCITY1 + b;
              for (int j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                d_stress->v[p][q][b][j] +=
                    -(d_mup->v[b][j] * k / lambda) * (delta(p, q) - fv->S[mode][p][q]) +
                    (d_lam->v[b][j] * k * mup / (lambda * lambda)) *
                        (delta(p, q) - fv->S[mode][p][q]);
              }
            }
          }
        }
      }
    } break;
    default: // Regular stress formulations
    {
      for (int mode = 0; mode < vn->modes; mode++) {
        for (int i = 0; i < VIM; i++) {
          for (int j = 0; j < VIM; j++) {
            for (int p = 0; p < VIM; p++) {
              for (int q = 0; q < VIM; q++) {
                int var = v_s[mode][p][q];
                for (int k = 0; k < ei[pg->imtrx]->dof[var]; k++) {
                  d_stress->S[i][j][mode][p][q][k] = delta(i, p) * delta(j, q) * bf[var]->phi[k];
                }
              }
            }
          }
        }
      }
    } break;
    }
  }
}

/*
 * Calculate the total stress tensor for a fluid at a single gauss point
 *  This includes the diagonal pressure contribution
 *
 *  Pi = stress tensor
 *  d_Pi = dependence of the stress tensor on the independent variables
 */
void fluid_stress(double Pi[DIM][DIM], STRESS_DEPENDENCE_STRUCT *d_Pi) {

  /*
   * Variables for vicosity and derivative
   */
  dbl mu = 0.0;
  VISCOSITY_DEPENDENCE_STRUCT d_mu_struct; /* viscosity dependence */
  VISCOSITY_DEPENDENCE_STRUCT *d_mu = &d_mu_struct;

  /* polymer viscosity and derivatives */

  dbl mup;
  VISCOSITY_DEPENDENCE_STRUCT d_mup_struct; /* viscosity dependence */
  VISCOSITY_DEPENDENCE_STRUCT *d_mup = &d_mup_struct;

  /*  shift function */
  dbl at = 0.0;
  dbl d_at_dT[MDE];
  dbl wlf_denom, temp;
  dbl mu_over_mu_num = 0.0;

  /* solvent viscosity and derivatives */

  dbl mus = 0.0;
  VISCOSITY_DEPENDENCE_STRUCT d_mus_struct; /* viscosity dependence */
  VISCOSITY_DEPENDENCE_STRUCT *d_mus = &d_mus_struct;

  /* numerical "adaptive" viscosity and derivatives */

  dbl mu_num;
  dbl d_mun_dS[MAX_MODES][DIM][DIM][MDE];
  dbl d_mun_dG[DIM][DIM][MDE];

  /* Dilational viscosity */
  dbl kappa = 0.0;
  DILVISCOSITY_DEPENDENCE_STRUCT d_dilMu_struct;
  DILVISCOSITY_DEPENDENCE_STRUCT *d_dilMu = &d_dilMu_struct;
  int kappaWipesMu = 1;
  dbl dilmuMult = 1.0;

  /* particle stress for suspension balance model*/
  dbl tau_p[DIM][DIM];
  dbl d_tau_p_dv[DIM][DIM][DIM][MDE];
  dbl d_tau_p_dvd[DIM][DIM][DIM][MDE];
  dbl d_tau_p_dy[DIM][DIM][MAX_CONC][MDE];
  dbl d_tau_p_dmesh[DIM][DIM][DIM][MDE];
  dbl d_tau_p_dp[DIM][DIM][MDE];
  int w0; /* suspension species number */

  dbl gamma[DIM][DIM];      /* shrearrate tensor based on velocity */
  dbl s[DIM][DIM];          /* polymer stress tensor */
  dbl gamma_cont[DIM][DIM]; /* shearrate tensor based on continuous gradient of velocity */
  dbl P;

  dbl evss_f; /* flag to add in the continuous velocity
                 gradient tensor for Fortin's formulation */
  int v_s[MAX_MODES][DIM][DIM];
  int v_g[DIM][DIM];
  int mode; /* index for modal viscoelastic counter */

  //! Flag for doing dilational viscosity contributions.
  int do_dilational_visc = 0;
  int dim;

  int a, b, p, q, j, w, c, var;

  dbl(*grad_phi_e)[DIM][DIM][DIM] = NULL;

  int eqn = R_MOMENTUM1;

  if (pd->gv[TEMPERATURE]) {
    temp = fv->T;
  } else {
    temp = upd->Process_Temperature;
  }

  dim = pd->Num_Dim;

  if (pd->gv[POLYMER_STRESS11]) {
    (void)stress_eqn_pointer(v_s);
  }
  if (pd->gv[VELOCITY_GRADIENT11]) {
    v_g[0][0] = VELOCITY_GRADIENT11;
    v_g[0][1] = VELOCITY_GRADIENT12;
    v_g[1][0] = VELOCITY_GRADIENT21;
    v_g[1][1] = VELOCITY_GRADIENT22;
    v_g[0][2] = VELOCITY_GRADIENT13;
    v_g[1][2] = VELOCITY_GRADIENT23;
    v_g[2][0] = VELOCITY_GRADIENT31;
    v_g[2][1] = VELOCITY_GRADIENT32;
    v_g[2][2] = VELOCITY_GRADIENT33;
  }

  /*
   * Field variables...
   */

  P = fv->P;

  /* if d_Pi == NULL, then the dependencies aren't needed,
     so we won't need the viscosity dependencies, either
  */
  if (d_Pi == NULL) {
    d_mu = NULL;
    d_mus = NULL;
    d_mup = NULL;
    d_dilMu = NULL;
  }

  /*
   * In Cartesian coordinates, this velocity gradient tensor will
   * have components that are...
   *
   * 			grad_v[a][b] = d v_b
   *				       -----
   *				       d x_a
   */

  memset(tau_p, 0, sizeof(double) * DIM * DIM);
  memset(d_tau_p_dv, 0, sizeof(double) * DIM * DIM * DIM * MDE);
  memset(d_tau_p_dvd, 0, sizeof(double) * DIM * DIM * DIM * MDE);
  memset(d_tau_p_dy, 0, sizeof(double) * DIM * DIM * MAX_CONC * MDE);
  memset(d_tau_p_dmesh, 0, sizeof(double) * DIM * DIM * DIM * MDE);
  memset(d_tau_p_dp, 0, sizeof(double) * DIM * DIM * MDE);

  /*   if( cr->MassFluxModel == DM_SUSPENSION_BALANCE|| mp->QTensorDiffType[w0] == CONSTANT)
   */
  if (cr->MassFluxModel == DM_SUSPENSION_BALANCE || cr->MassFluxModel == HYDRODYNAMIC_QTENSOR) {
    w0 = gn->sus_species_no;
    particle_stress(tau_p, d_tau_p_dv, d_tau_p_dvd, d_tau_p_dy, d_tau_p_dmesh, d_tau_p_dp, w0);
  }

  evss_f = 0;
  if (pd->gv[POLYMER_STRESS11] && is_evss_f_model(vn->evssModel)) {
    evss_f = 1.0;
  }

  double Heaviside = 1;
  if (ls != NULL && ls->ghost_stress) {
    load_lsi(ls->Length_Scale);
    switch (ls->ghost_stress) {
    case LS_OFF:
      Heaviside = 1;
      break;
    case LS_POSITIVE:
      Heaviside = lsi->H;
      break;
    case LS_NEGATIVE:
      Heaviside = 1 - lsi->H;
      break;
    default:
      GOMA_EH(GOMA_ERROR, "Unknown Level Set Ghost Stress value");
      break;
    }
  }

  if (evss_f) {
    for (a = 0; a < VIM; a++) {
      for (b = 0; b < VIM; b++) {
        gamma_cont[a][b] = fv->G[a][b] + fv->G[b][a];
      }
    }
  } else {
    memset(gamma_cont, 0, sizeof(double) * DIM * DIM);
  }

  /* load up shear rate tensor based on velocity */
  for (a = 0; a < VIM; a++) {
    for (b = 0; b < VIM; b++) {
      gamma[a][b] = fv->grad_v[a][b] + fv->grad_v[b][a];
    }
  }
  if (do_dilational_visc) {
    // use previously calculated div_v.
  }

  mu = viscosity(gn, gamma, d_mu);
  if (pd->gv[POLYMER_STRESS11]) {
    mus = viscosity(gn, gamma, d_mus);

    /* initialize the derivative wrt to stress and velocity gradient */

    memset(d_mun_dS, 0, sizeof(double) * MAX_MODES * DIM * DIM * MDE);
    memset(d_mun_dG, 0, sizeof(double) * DIM * DIM * MDE);

    /* This is the adaptive viscosity from Sun et al., 1999.
     * The term multiplies the continuous and discontinuous
     * shear-rate, so it should cancel out and not affect the
     * solution, other than increasing the stability of the
     * algorithm in areas of high shear and stress.
     */

    mu_num = 1;
    if (DOUBLE_NONZERO(vn->eps)) {
      for (int mode = 0; mode < vn->modes; mode++) {
        for (a = 0; a < VIM; a++) {
          for (b = 0; b < VIM; b++) {
            s[a][b] += fv->S[mode][a][b];
          }
        }
      }

      mu_num = numerical_viscosity(s, gamma_cont, d_mun_dS, d_mun_dG);
    }

    mu_over_mu_num = mus;
    mu = mu_num * mus;

    /* first add the solvent viscosity to the total viscosity
     * including all the derivatives.
     */
    var = VELOCITY1;
    if (d_Pi != NULL && pd->v[pg->imtrx][var]) {
      for (a = 0; a < WIM; a++) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_mu->v[a][j] = mu_num * d_mus->v[a][j];
        }
      }
    }

    var = MESH_DISPLACEMENT1;
    if (d_Pi != NULL && pd->v[pg->imtrx][var]) {
      for (a = 0; a < dim; a++) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_mu->X[a][j] = mu_num * d_mus->X[a][j];
        }
      }
    }

    var = TEMPERATURE;
    if (d_Pi != NULL && pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_mu->T[j] = mu_num * d_mus->T[j];
      }
    }

    var = BOND_EVOLUTION;
    if (d_Pi != NULL && pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_mu->nn[j] = mu_num * d_mus->nn[j];
      }
    }

    var = RESTIME;
    if (d_Pi != NULL && pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_mu->degrade[j] = mu_num * d_mus->degrade[j];
      }
    }

    var = FILL;
    if (d_Pi != NULL && pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_mu->F[j] = mu_num * d_mus->F[j];
      }
    }

    var = EDDY_NU;
    if (d_Pi != NULL && pd->v[pg->imtrx][var]) {
      for (p = 0; p < VIM; p++) {
        for (q = 0; q < VIM; q++) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            d_mu->eddy_nu[j] = mu_num * d_mu->eddy_nu[j];
          }
        }
      }
    }

    var = TURB_K;
    if (d_Pi != NULL && pd->v[pg->imtrx][var]) {
      for (p = 0; p < VIM; p++) {
        for (q = 0; q < VIM; q++) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            d_mu->turb_k[j] = mu_num * d_mu->turb_k[j];
          }
        }
      }
    }

    var = TURB_OMEGA;
    if (d_Pi != NULL && pd->v[pg->imtrx][var]) {
      for (p = 0; p < VIM; p++) {
        for (q = 0; q < VIM; q++) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            d_mu->turb_omega[j] = mu_num * d_mu->turb_omega[j];
          }
        }
      }
    }

    var = PRESSURE;
    if (d_Pi != NULL && pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_mu->P[j] = mu_num * d_mus->P[j];
      }
    }

    var = MASS_FRACTION;
    if (d_Pi != NULL && pd->v[pg->imtrx][var]) {
      for (w = 0; w < pd->Num_Species_Eqn; w++) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_mu->C[w][j] = mu_num * d_mus->C[w][j];
        }
      }
    }

    /*  shift factor  */
    if (pd->e[pg->imtrx][TEMPERATURE]) {
      if (vn->shiftModel == CONSTANT) {
        at = vn->shift[0];
        for (j = 0; j < ei[pg->imtrx]->dof[TEMPERATURE]; j++) {
          d_at_dT[j] = 0.;
        }
      } else if (vn->shiftModel == MODIFIED_WLF) {
        wlf_denom = vn->shift[1] + temp - mp->reference[TEMPERATURE];
        if (wlf_denom != 0.) {
          at = exp(vn->shift[0] * (mp->reference[TEMPERATURE] - temp) / wlf_denom);
          for (j = 0; j < ei[pg->imtrx]->dof[TEMPERATURE]; j++) {
            d_at_dT[j] = -at * vn->shift[0] * vn->shift[1] / (wlf_denom * wlf_denom) *
                         bf[TEMPERATURE]->phi[j];
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

    for (mode = 0; mode < vn->modes; mode++) {
      /* get polymer viscosity */
      mup = viscosity(ve[mode]->gn, gamma, d_mup);

      mu_over_mu_num += Heaviside * at * mup;
      mu += Heaviside * mu_num * at * mup;

      var = VELOCITY1;
      if (d_Pi != NULL && pd->v[pg->imtrx][var]) {
        for (a = 0; a < WIM; a++) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            d_mu->v[a][j] += mu_num * at * d_mup->v[a][j];
          }
        }
      }

      var = MESH_DISPLACEMENT1;
      if (d_Pi != NULL && pd->v[pg->imtrx][var]) {
        for (a = 0; a < dim; a++) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            d_mu->X[a][j] += mu_num * at * d_mup->X[a][j];
          }
        }
      }

      var = TEMPERATURE;
      if (d_Pi != NULL && pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_mu->T[j] += mu_num * (at * d_mup->T[j] + mup * d_at_dT[j]);
        }
      }

      var = BOND_EVOLUTION;
      if (d_Pi != NULL && pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_mu->nn[j] += mu_num * at * d_mup->nn[j];
        }
      }

      var = RESTIME;
      if (d_Pi != NULL && pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_mu->degrade[j] += mu_num * at * d_mup->degrade[j];
        }
      }

      var = FILL;
      if (d_Pi != NULL && pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_mu->F[j] += mu_num * at * d_mup->F[j];
        }
      }

      var = PRESSURE;
      if (d_Pi != NULL && pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_mu->P[j] += mu_num * at * d_mup->P[j];
        }
      }

      var = MASS_FRACTION;
      if (d_Pi != NULL && pd->v[pg->imtrx][var]) {
        for (w = 0; w < pd->Num_Species_Eqn; w++) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            d_mu->C[w][j] += mu_num * at * d_mup->C[w][j];
          }
        }
      }
    } // for mode
  }   // if POLYMER_STRESS

  /*
   * Calculate the dilational viscosity, if necessary
   */
  if (mp->DilationalViscosityModel != DILVISCM_KAPPAWIPESMU) {
    kappa = dil_viscosity(gn, mu, d_mu, d_dilMu);
    dilmuMult = mp->dilationalViscosityMultiplier;
    // kappa = 0.0;
    kappaWipesMu = 0;
  }

  /*
   * Viscoelastic Stress contributions
   */
  dbl polymer_stress[DIM][DIM];
  STRESS_DEPENDENCE_STRUCT d_polymer_stress_struct;
  STRESS_DEPENDENCE_STRUCT *d_polymer_stress = NULL;
  if (pd->gv[POLYMER_STRESS11]) {
    if (d_Pi != NULL) {
      d_polymer_stress = &d_polymer_stress_struct;
      memset(d_polymer_stress, 0, sizeof(STRESS_DEPENDENCE_STRUCT));
    }
    ve_polymer_stress(gamma, polymer_stress, d_polymer_stress);
  }

  /* xxxxx HKM xxxx  -> operational point */
  /*
   *  tau_p[a][b] is the particle stress contribution, so
   *              it's zero for non-particle cases
   *  gamma[a][b] is the usual gradV + gradV_t term
   *  delta(a,b) is a delta function
   */
  for (a = 0; a < VIM; a++) {
    for (b = 0; b < VIM; b++) {
      Pi[a][b] = -P * (double)delta(a, b) + mu * gamma[a][b] - tau_p[a][b];
    }

    // Add in the diagonal contribution
    if (!kappaWipesMu) {
      Pi[a][a] -= dilmuMult * (mu / 3.0 - 0.5 * kappa) * gamma[a][a];
    }
  }

  if (pd->gv[POLYMER_STRESS11]) {
    for (a = 0; a < VIM; a++) {
      for (b = 0; b < VIM; b++) {
        // TODO : derivative may be missing here...
        Pi[a][b] += -evss_f * (mu - mus) * gamma_cont[a][b] + Heaviside * polymer_stress[a][b];
      }
    }
  }
  if (gn->ConstitutiveEquation == BINGHAM_MIXED) {
    dbl tau_y = gn->tau_y;
    if (pd->gv[VELOCITY_GRADIENT11]) {
      for (a = 0; a < VIM; a++) {
        for (b = 0; b < VIM; b++) {
          Pi[a][b] += tau_y * fv->G[a][b];
        }
      }
    } else {
      GOMA_EH(-1, "BINGHAM_MIXED but no VELOCITY_GRADIENT equations");
    }
  }

  /*
   * OK, FIND THE JACOBIAN
   */
  var = TEMPERATURE;
  if (d_Pi != NULL && pd->v[pg->imtrx][var]) {
    for (p = 0; p < VIM; p++) {
      for (q = 0; q < VIM; q++) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_Pi->T[p][q][j] = d_mu->T[j] * gamma[p][q];
        }
      }
    }
    if (!kappaWipesMu) {
      for (p = 0; p < VIM; p++) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_Pi->T[p][p][j] -= dilmuMult * (d_mu->T[j] / 3.0 - 0.5 * d_dilMu->T[j]) * gamma[p][p];
        }
      }
    }
    if (pd->gv[POLYMER_STRESS11]) {
      for (p = 0; p < VIM; p++) {
        for (q = 0; q < VIM; q++) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            d_Pi->T[p][q][j] -= evss_f * (d_mu->T[j] - d_mus->T[j]) * gamma_cont[p][q];
          }
        }
      }
    }
  }

  var = BOND_EVOLUTION;
  if (d_Pi != NULL && pd->v[pg->imtrx][var]) {
    for (p = 0; p < VIM; p++) {
      for (q = 0; q < VIM; q++) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_Pi->nn[p][q][j] = d_mu->nn[j] * gamma[p][q];
        }
      }
    }
    if (!kappaWipesMu) {
      for (p = 0; p < VIM; p++) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_Pi->nn[p][p][j] -= dilmuMult * (d_mu->nn[j] / 3.0 - 0.5 * d_dilMu->nn[j]) * gamma[p][p];
        }
      }
    }
    if (pd->gv[POLYMER_STRESS11]) {
      for (p = 0; p < VIM; p++) {
        for (q = 0; q < VIM; q++) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            d_Pi->nn[p][q][j] -= evss_f * (d_mu->nn[j] - d_mus->nn[j]) * gamma_cont[p][q];
          }
        }
      }
    }
  }

  var = RESTIME;
  if (d_Pi != NULL && pd->v[pg->imtrx][var]) {
    for (p = 0; p < VIM; p++) {
      for (q = 0; q < VIM; q++) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_Pi->degrade[p][q][j] = d_mu->degrade[j] * gamma[p][q];
        }
      }
    }
    if (!kappaWipesMu) {
      for (p = 0; p < VIM; p++) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_Pi->degrade[p][p][j] -=
              (d_mu->degrade[j] / 3.0 - 0.5 * d_dilMu->degrade[j]) * gamma[p][p];
        }
      }
    }
    if (pd->v[pg->imtrx][POLYMER_STRESS11]) {
      for (p = 0; p < VIM; p++) {
        for (q = 0; q < VIM; q++) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            d_Pi->degrade[p][q][j] -=
                evss_f * (d_mu->degrade[j] - d_mus->degrade[j]) * gamma_cont[p][q];
          }
        }
      }
    }
  }

  var = FILL;
  if (d_Pi != NULL && pd->v[pg->imtrx][var]) {
    for (p = 0; p < VIM; p++) {
      for (q = 0; q < VIM; q++) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_Pi->F[p][q][j] = d_mu->F[j] * gamma[p][q];
        }
      }
    }
    if (!kappaWipesMu) {
      for (p = 0; p < VIM; p++) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_Pi->F[p][p][j] -= dilmuMult * (d_mu->F[j] / 3.0 - 0.5 * d_dilMu->F[j]) * gamma[p][p];
        }
      }
    }
    if (pd->gv[POLYMER_STRESS11]) {
      for (p = 0; p < VIM; p++) {
        for (q = 0; q < VIM; q++) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            d_Pi->F[p][q][j] -= evss_f * (d_mu->F[j] - d_mus->F[j]) * gamma_cont[p][q];
          }
        }
      }
    }
  }

  if (d_Pi != NULL && pd->v[pg->imtrx][PHASE1]) {
    for (p = 0; p < VIM; p++) {
      for (q = 0; q < VIM; q++) {
        for (a = 0; a < pfd->num_phase_funcs; a++) {
          var = PHASE1 + a;
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            d_Pi->pf[p][q][a][j] = d_mu->pf[a][j] * gamma[p][q];
          }
        }
      }
    }
    if (!kappaWipesMu) {
      for (p = 0; p < VIM; p++) {
        for (a = 0; a < pfd->num_phase_funcs; a++) {
          var = PHASE1 + a;
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            d_Pi->pf[p][p][a][j] -=
                dilmuMult * (d_mu->pf[a][j] / 3.0 - 0.5 * d_dilMu->pf[a][j]) * gamma[p][p];
          }
        }
      }
    }
    if (pd->gv[POLYMER_STRESS11]) {
      for (p = 0; p < VIM; p++) {
        for (q = 0; q < VIM; q++) {
          for (a = 0; a < pfd->num_phase_funcs; a++) {
            var = PHASE1 + a;
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              d_Pi->pf[p][q][a][j] -=
                  evss_f * (d_mu->pf[a][j] - d_mus->pf[a][j]) * gamma_cont[p][q];
            }
          }
        }
      }
    }
  }

  if (d_Pi != NULL && pd->v[pg->imtrx][VELOCITY1]) {
    /* Damn... It is unfortunate that grad_phi_e is 1) assumed to be
     * the same for all velocity components (this is bad for 3d
     * stability of a 2d flow), and 2) there is no bf[] structure
     * for the theta-velocity in cylindrical coordinates (I believe
     * there is for swirling).  That's why this needed to be
     * split... */
    if (pd->CoordinateSystem != CYLINDRICAL) {
      for (p = 0; p < VIM; p++) {
        for (q = 0; q < VIM; q++) {
          for (b = 0; b < WIM; b++) {
            for (j = 0; j < ei[pg->imtrx]->dof[VELOCITY1]; j++) {
              /* grad_phi_e cannot be the same for all
               * velocities for 3d stab of 2d flow!!
               * Compare with the old way in the CYLINDRICAL
               * chunk below... */
              d_Pi->v[p][q][b][j] = mu * (bf[VELOCITY1 + q]->grad_phi_e[j][b][p][q] +
                                          bf[VELOCITY1 + p]->grad_phi_e[j][b][q][p]) +
                                    d_mu->v[b][j] * gamma[p][q] - d_tau_p_dv[p][q][b][j];
            }
          }
        }

        if (!kappaWipesMu) {
          for (b = 0; b < WIM; b++) {
            for (j = 0; j < ei[pg->imtrx]->dof[VELOCITY1]; j++) {
              d_Pi->v[p][p][b][j] -=
                  dilmuMult *
                  ((2.0 * mu / 3.0 - kappa) * (bf[VELOCITY1 + p]->grad_phi_e[j][b][p][p]) +
                   (d_mu->v[b][j] / 3.0 - 0.5 * d_dilMu->v[b][j]) * gamma[p][p]);
            }
          }
        }
      }

    } else {
      /* For CYLINDRICAL, we can assume that all of the velocity
       * components share the same interpolating basis function.
       * This was not allowable for 3d stability of 2d flow and
       * the PROJECTED_CARTESIAN coordinate system.  In fact, the
       * behavior above is "more" correct than this, so I
       * defaulted everyrthing else to it.
       */
      grad_phi_e = bf[eqn]->grad_phi_e;
      for (p = 0; p < VIM; p++) {
        for (q = 0; q < VIM; q++) {
          for (b = 0; b < WIM; b++) {
            for (j = 0; j < ei[pg->imtrx]->dof[VELOCITY1]; j++) {
              d_Pi->v[p][q][b][j] = mu * (grad_phi_e[j][b][p][q] + grad_phi_e[j][b][q][p]) +
                                    d_mu->v[b][j] * gamma[p][q] - d_tau_p_dv[p][q][b][j];
            }
          }
        }
        if (!kappaWipesMu) {
          for (b = 0; b < WIM; b++) {
            for (j = 0; j < ei[pg->imtrx]->dof[VELOCITY1]; j++) {
              d_Pi->v[p][p][b][j] -=
                  dilmuMult * ((2.0 * mu / 3.0 - kappa) * (grad_phi_e[j][b][p][p]) +
                               (d_mu->v[b][j] / 3.0 - 0.5 * d_dilMu->v[b][j]) * gamma[p][p]);
            }
          }
        }
      }
    }
    if (pd->v[pg->imtrx][POLYMER_STRESS11]) {
      for (p = 0; p < VIM; p++) {
        for (q = 0; q < VIM; q++) {
          for (b = 0; b < WIM; b++) {
            for (j = 0; j < ei[pg->imtrx]->dof[VELOCITY1]; j++) {
              d_Pi->v[p][q][b][j] -= evss_f * (d_mu->v[b][j] - d_mus->v[b][j]) * gamma_cont[p][q];
            }
          }
        }
      }
    }
  }

  // Vorticity direction dependence for qtensor
  if (d_Pi != NULL && pd->v[pg->imtrx][VORT_DIR1]) {
    memset(d_Pi->vd, 0, DIM * DIM * DIM * MDE * sizeof(double));
    for (p = 0; p < VIM; p++) {
      for (q = 0; q < VIM; q++) {
        for (b = 0; b < dim; b++) {
          for (j = 0; j < ei[pg->imtrx]->dof[VORT_DIR1]; j++) {
            d_Pi->vd[p][q][b][j] = -d_tau_p_dvd[p][q][b][j];
          }
        }
      }
    }
  }

  // Mesh Dependence
  if (d_Pi != NULL && pd->v[pg->imtrx][MESH_DISPLACEMENT1]) {
    for (p = 0; p < VIM; p++) {
      for (q = 0; q < VIM; q++) {
        for (b = 0; b < dim; b++) {
          for (j = 0; j < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; j++) {
            d_Pi->X[p][q][b][j] =
                mu * (fv->d_grad_v_dmesh[p][q][b][j] + fv->d_grad_v_dmesh[q][p][b][j]) +
                +d_mu->X[b][j] * gamma[p][q] - d_tau_p_dmesh[p][q][b][j];
          }
        }
      }
      if (!kappaWipesMu) {
        for (b = 0; b < dim; b++) {
          for (j = 0; j < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; j++) {
            d_Pi->X[p][p][b][j] -=
                dilmuMult * ((2.0 * mu / 3.0 - kappa) * (fv->d_grad_v_dmesh[p][p][b][j]) +
                             (d_mu->X[b][j] / 3.0 - 0.5 * d_dilMu->X[b][j]) * gamma[p][p]);
          }
        }
      }
    }
    if (pd->gv[POLYMER_STRESS11]) {
      for (p = 0; p < VIM; p++) {
        for (q = 0; q < VIM; q++) {
          for (b = 0; b < dim; b++) {
            for (j = 0; j < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; j++) {
              d_Pi->X[p][q][b][j] -= evss_f * (d_mu->X[b][j] - d_mus->X[b][j]) * gamma_cont[p][q];
            }
          }
        }
      }
    }
  }

  if (d_Pi != NULL && pd->v[pg->imtrx][POLYMER_STRESS11]) {
    for (mode = 0; mode < vn->modes; mode++) {
      for (p = 0; p < VIM; p++) {
        for (q = 0; q < VIM; q++) {
          for (b = 0; b < VIM; b++) {
            for (c = 0; c < VIM; c++) {
              var = v_s[mode][b][c];
              for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                d_Pi->S[p][q][mode][b][c][j] = d_polymer_stress->S[p][q][mode][b][c][j] +
                                               mu_over_mu_num * d_mun_dS[mode][b][c][j] *
                                                   (gamma[p][q] - evss_f * gamma_cont[p][q]);
              }
            }
          }
        }
      }
    }
  }

  if (d_Pi != NULL && pd->v[pg->imtrx][VELOCITY_GRADIENT11] && pd->gv[POLYMER_STRESS11]) {
    for (p = 0; p < VIM; p++) {
      for (q = 0; q < VIM; q++) {
        for (b = 0; b < VIM; b++) {
          for (c = 0; c < VIM; c++) {
            var = v_g[b][c];
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              d_Pi->g[p][q][b][c][j] =
                  mu_over_mu_num * d_mun_dG[b][c][j] * (gamma[p][q] - evss_f * gamma_cont[p][q]) +
                  evss_f * ((mu - mus) * bf[var]->phi[j] *
                            ((double)delta(c, p) * (double)delta(b, q) +
                             (double)delta(b, p) * (double)delta(c, q)));
            }
          }
        }
      }
    }
  }
  if (d_Pi != NULL && pd->gv[VELOCITY_GRADIENT11] && gn->ConstitutiveEquation == BINGHAM_MIXED) {
    dbl tau_y = gn->tau_y;
    for (p = 0; p < VIM; p++) {
      for (q = 0; q < VIM; q++) {
        if (p <= q) {
          for (b = 0; b < VIM; b++) {
            for (c = 0; c < VIM; c++) {
              var = v_g[b][c];
              for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                d_Pi->g[p][q][b][c][j] = tau_y * bf[var]->phi[j] * (delta(p, b) * delta(q, c));
                if (q > p) {
                  d_Pi->g[q][p][b][c][j] = tau_y * bf[var]->phi[j] * (delta(p, b) * delta(q, c));
                }
              }
            }
          }
        }
      }
    }
  }

  var = MASS_FRACTION;
  if (d_mu != NULL && pd->v[pg->imtrx][var]) {
    for (p = 0; p < VIM; p++) {
      for (q = 0; q < VIM; q++) {
        for (w = 0; w < pd->Num_Species_Eqn; w++) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            d_Pi->C[p][q][w][j] = d_mu->C[w][j] * gamma[p][q] - d_tau_p_dy[p][q][w][j];
          }
        }
      }
      if (!kappaWipesMu) {
        for (w = 0; w < pd->Num_Species_Eqn; w++) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            d_Pi->C[p][p][w][j] -=
                dilmuMult * (d_mu->C[w][j] / 3.0 - 0.5 * d_dilMu->C[w][j]) * gamma[p][p];
          }
        }
      }
    }

    if (pd->gv[POLYMER_STRESS11]) {
      for (p = 0; p < VIM; p++) {
        for (q = 0; q < VIM; q++) {
          for (w = 0; w < pd->Num_Species_Eqn; w++) {
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              d_Pi->C[p][q][w][j] -= evss_f * (d_mu->C[w][j] - d_mus->C[w][j]) * gamma_cont[p][q];
            }
          }
        }
      }
    }
  }

  var = EDDY_NU;
  if (d_Pi != NULL && pd->v[pg->imtrx][var]) {
    for (p = 0; p < VIM; p++) {
      for (q = 0; q < VIM; q++) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_Pi->eddy_nu[p][q][j] = d_mu->eddy_nu[j] * gamma[p][q];
        }
      }
    }
  }

  var = TURB_K;
  if (d_Pi != NULL && pd->v[pg->imtrx][var]) {
    for (p = 0; p < VIM; p++) {
      for (q = 0; q < VIM; q++) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_Pi->turb_k[p][q][j] = d_mu->turb_k[j] * gamma[p][q];
        }
      }
    }
  }
  var = TURB_OMEGA;
  if (d_Pi != NULL && pd->v[pg->imtrx][var]) {
    for (p = 0; p < VIM; p++) {
      for (q = 0; q < VIM; q++) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_Pi->turb_omega[p][q][j] = d_mu->turb_omega[j] * gamma[p][q];
        }
      }
    }
  }

  var = PRESSURE;
  if (d_Pi != NULL && pd->v[pg->imtrx][var]) {
    for (p = 0; p < VIM; p++) {
      for (q = 0; q < VIM; q++) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_Pi->P[p][q][j] = -(double)delta(p, q) * bf[var]->phi[j] + d_mu->P[j] * gamma[p][q] -
                             d_tau_p_dp[p][q][j];
        }
      }
      if (!kappaWipesMu) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_Pi->P[p][p][j] -= dilmuMult * (d_mu->P[j] / 3.0 - 0.5 * d_dilMu->P[j]) * gamma[p][p];
        }
      }
    }
    if (pd->gv[POLYMER_STRESS11]) {
      for (p = 0; p < VIM; p++) {
        for (q = 0; q < VIM; q++) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            d_Pi->P[p][q][j] -= evss_f * (d_mu->P[j] - d_mus->P[j]) * gamma_cont[p][q];
          }
        }
      }
    }
  }
} /******************************************************************************
   * momentum_source_term(): Computes the body force term for the momentum balance.
   *
   * Input
   * -----
   *   f    == Body force
   *   df->T == Derivative w.r.t. temperature
   *   df->X == Derivative w.r.t. mesh displacements
   *   df->C == Derivative w.r.t. concentration
   *   df->v == Derivative w.r.t. velocity
   *   df->F == Derivative w.r.t. FILL
   *   df->E == Derivative w.r.t. electric field
   *
   ******************************************************************************/
int momentum_source_term(dbl f[DIM], /* Body force. */
                         MOMENTUM_SOURCE_DEPENDENCE_STRUCT *df,
                         dbl time) {
  int j, a, b, w;
  int eqn, var, var_offset;
  int err;
  const int dim = pd->Num_Dim;
  int siz;
  int status = 0;
  double *phi;
  struct Level_Set_Data *ls_old;

  /* initialize everything to zero */
  siz = sizeof(double) * DIM;
  memset(f, 0, siz);

  if (df != NULL) {
    siz = sizeof(double) * DIM * MDE;
    memset(df->T, 0, siz);
    memset(df->F, 0, siz);
    memset(df->ars, 0, siz);

    siz = sizeof(double) * DIM * DIM * MDE;
    memset(df->X, 0, siz);
    memset(df->v, 0, siz);
    memset(df->E, 0, siz);

    siz = sizeof(double) * DIM * MAX_CONC * MDE;
    memset(df->C, 0, siz);
  }

  /****Momentum Source Model******/
  if (mp->MomentumSourceModel == USER) {
    err = usr_momentum_source(mp->u_momentum_source);

    for (a = 0; a < dim; a++) {
      eqn = R_MOMENTUM1 + a;
      if (pd->e[upd->matrix_index[eqn]][eqn] & T_SOURCE) {
        f[a] = mp->momentum_source[a];
        var = TEMPERATURE;
        if (df != NULL && pd->v[pg->imtrx][var]) {
          phi = bf[var]->phi;
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            df->T[a][j] = mp->d_momentum_source[a][var] * phi[j];
          }
        }
        if (df != NULL && pd->v[pg->imtrx][FILL]) {
          var = FILL;
          phi = bf[var]->phi;
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            df->F[a][j] = mp->d_momentum_source[a][var] * phi[j];
          }
        }
        if (df != NULL && pd->v[pg->imtrx][MESH_DISPLACEMENT1]) {
          for (b = 0; b < dim; b++) {
            var = MESH_DISPLACEMENT1 + b;
            phi = bf[var]->phi;
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              df->X[a][b][j] = mp->d_momentum_source[a][var] * phi[j];
            }
          }
        }

        if (df != NULL && pd->v[pg->imtrx][MASS_FRACTION]) {
          for (w = 0; w < pd->Num_Species_Eqn; w++) {
            var = MASS_FRACTION;
            var_offset = MAX_VARIABLE_TYPES + w;
            phi = bf[var]->phi;
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              df->C[a][w][j] = mp->d_momentum_source[a][var_offset] * phi[j];
            }
          }
        }

        if (df != NULL && pd->v[pg->imtrx][VELOCITY1]) {
          for (b = 0; b < DIM; b++) {
            var = VELOCITY1 + b;
            phi = bf[var]->phi;
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              df->v[a][b][j] = mp->d_momentum_source[a][var] * phi[j];
            }
          }
        }
      }
    }
  } else if (mp->MomentumSourceModel == CONSTANT) {
    int force_dim = dim;
    if (pd->CoordinateSystem == CARTESIAN_2pt5D) {
      force_dim = 3;
    }
    for (a = 0; a < force_dim; a++) {
      eqn = R_MOMENTUM1 + a;
      if (pd->e[upd->matrix_index[eqn]][eqn] & T_SOURCE) {
        f[a] = mp->momentum_source[a];
      }
    }
  } else if (mp->MomentumSourceModel == VARIABLE_DENSITY) {
    if (mp->DensityModel == SOLVENT_POLYMER) {
      DENSITY_DEPENDENCE_STRUCT d_rho_struct; /* density dependence */
      DENSITY_DEPENDENCE_STRUCT *d_rho = &d_rho_struct;
      double rho = density(d_rho, time);
      for (a = 0; a < dim; a++) {
        eqn = R_MOMENTUM1 + a;
        if (pd->e[upd->matrix_index[eqn]][eqn] & T_SOURCE) {
          f[a] = rho * mp->momentum_source[a];
        }
        if (df != NULL && pd->v[pg->imtrx][MASS_FRACTION]) {
          for (w = 0; w < pd->Num_Species_Eqn; w++) {
            var = MASS_FRACTION;
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              df->C[a][w][j] = d_rho->C[w][j] * mp->momentum_source[a];
            }
          }
        }
      }
    } else if (mp->DensityModel == DENSITY_FOAM_PMDI_10 || mp->DensityModel == DENSITY_FOAM_PBE ||
               mp->DensityModel == DENSITY_FOAM || mp->DensityModel == DENSITY_FOAM_TIME ||
               mp->DensityModel == DENSITY_FOAM_TIME_TEMP ||
               mp->DensityModel == DENSITY_MOMENT_BASED) {
      DENSITY_DEPENDENCE_STRUCT d_rho_struct; /* density dependence */
      DENSITY_DEPENDENCE_STRUCT *d_rho = &d_rho_struct;
      double rho = density(d_rho, time);
      for (a = 0; a < dim; a++) {
        eqn = R_MOMENTUM1 + a;
        if (pd->e[upd->matrix_index[eqn]][eqn] & T_SOURCE) {
          f[a] = rho * mp->momentum_source[a];
        }
        if (df != NULL && pd->v[pg->imtrx][MASS_FRACTION]) {
          for (w = 0; w < pd->Num_Species_Eqn; w++) {
            var = MASS_FRACTION;
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              df->C[a][w][j] = d_rho->C[w][j] * mp->momentum_source[a];
            }
          }
        }
        if (df != NULL && pd->v[pg->imtrx][TEMPERATURE]) {
          var = TEMPERATURE;
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            df->T[a][j] = d_rho->T[j] * mp->momentum_source[a];
          }
        }
        if (df != NULL && pd->v[pg->imtrx][FILL]) {
          var = FILL;
          phi = bf[var]->phi;
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            df->F[a][j] = d_rho->F[j] * mp->momentum_source[a];
          }
        }
      }
    } else {
      GOMA_EH(GOMA_ERROR, "Unknown density model for variable density");
    }

  } else if (mp->MomentumSourceModel == VARIABLE_DENSITY_NO_GAS) {
    if (mp->DensityModel == DENSITY_FOAM_PMDI_10 || mp->DensityModel == DENSITY_FOAM_PBE ||
        mp->DensityModel == DENSITY_FOAM || mp->DensityModel == DENSITY_FOAM_TIME ||
        mp->DensityModel == DENSITY_FOAM_TIME_TEMP || mp->DensityModel == DENSITY_MOMENT_BASED) {
      DENSITY_DEPENDENCE_STRUCT d_rho_struct; /* density dependence */
      DENSITY_DEPENDENCE_STRUCT *d_rho = &d_rho_struct;
      double rho = density(d_rho, time);
      load_lsi(ls->Length_Scale);
      double Heaviside;
      if (mp->mp2nd->densitymask[0] == 0) {
        Heaviside = 1 - lsi->H;
      } else {
        Heaviside = lsi->H;
      }
      for (a = 0; a < dim; a++) {
        eqn = R_MOMENTUM1 + a;
        if (pd->e[upd->matrix_index[eqn]][eqn] & T_SOURCE) {
          f[a] = Heaviside * rho * mp->momentum_source[a];
        }
        if (df != NULL && pd->v[pg->imtrx][MASS_FRACTION]) {
          for (w = 0; w < pd->Num_Species_Eqn; w++) {
            var = MASS_FRACTION;
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              df->C[a][w][j] = Heaviside * d_rho->C[w][j] * mp->momentum_source[a];
            }
          }
        }
        if (df != NULL && pd->v[pg->imtrx][TEMPERATURE]) {
          var = TEMPERATURE;
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            df->T[a][j] = Heaviside * d_rho->T[j] * mp->momentum_source[a];
          }
        }
        if (df != NULL && pd->v[pg->imtrx][FILL]) {
          var = FILL;
          phi = bf[var]->phi;
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            df->F[a][j] = Heaviside * d_rho->F[j] * mp->momentum_source[a];
          }
        }
      }
    } else {
      GOMA_EH(GOMA_ERROR, "Unknown density model for variable density no gas");
    }

  } else if (mp->MomentumSourceModel == SUSPENSION_PM) {
    err = suspension_pm_fluid_momentum_source(f, df);
    GOMA_EH(err, "Problems in suspension_pm_fluid_momentum_source");
  } else if (mp->MomentumSourceModel == SUSPEND || mp->MomentumSourceModel == SUSPENSION) {
    err = suspend_momentum_source(f, df);
    GOMA_EH(err, "Problems in suspend_momentum_source");
  } else if (mp->MomentumSourceModel == BOUSS) {
    err = bouss_momentum_source(f, df, 0, TRUE);
    GOMA_EH(err, "Problems in bouss_momentum_source");
  } else if (mp->MomentumSourceModel == BOUSS_JXB) {
    err = bouss_momentum_source(f, df, 1, FALSE);
    GOMA_EH(err, "Problems in bouss_momentum_source");
  } else if (mp->MomentumSourceModel == BOUSSINESQ) {
    err = bouss_momentum_source(f, df, 0, FALSE);
    GOMA_EH(err, "Problems in bouss_momentum_source");
  } else if (mp->MomentumSourceModel == EHD_POLARIZATION) {
    err = EHD_POLARIZATION_source(f, df);
    GOMA_EH(err, "Problems in EHD_POLARIZATION force routine");
  } else if (mp->MomentumSourceModel == GRAV_VIBRATIONAL) {
    err = gravity_vibrational_source(f, df, time);
    GOMA_EH(err, "Problems in GRAV_VIBRATIONAL force routine");
  } else if (mp->MomentumSourceModel == FILL_SRC) {
    err = fill_momentum_source(f);
  } else if (mp->MomentumSourceModel == LEVEL_SET && pd->gv[FILL]) {
    DENSITY_DEPENDENCE_STRUCT d_rho_struct; /* density dependence */
    DENSITY_DEPENDENCE_STRUCT *d_rho = &d_rho_struct;
    double rho = density(d_rho, time);

    for (a = 0; a < pd->Num_Dim; a++) {
      f[a] = mp->momentum_source[a] * rho;
      if (df != NULL) {
        for (j = 0; j < ei[pg->imtrx]->dof[FILL]; j++)
          df->F[a][j] = mp->momentum_source[a] * d_rho->F[j];
      }
    }
  } else if (mp->MomentumSourceModel == LEVEL_SET && pd->gv[PHASE1]) {
    DENSITY_DEPENDENCE_STRUCT d_rho_struct; /* density dependence */
    DENSITY_DEPENDENCE_STRUCT *d_rho = &d_rho_struct;
    double rho = density(d_rho, time);

    for (a = 0; a < pd->Num_Dim; a++) {
      f[a] = mp->momentum_source[a] * rho;

      if (df != NULL) {
        for (b = 0; b < pfd->num_phase_funcs; b++) {
          for (j = 0; j < ei[pg->imtrx]->dof[PHASE1]; j++) {
            df->pf[a][b][j] = mp->momentum_source[a] * d_rho->pf[b][j];
          }
        }
      }
    }
  } else if (mp->MomentumSourceModel == ACOUSTIC) {
    for (a = 0; a < dim; a++) {
      eqn = R_MOMENTUM1 + a;
      if (pd->e[upd->matrix_index[eqn]][eqn] & T_SOURCE) {
        /*  Graviational piece	*/
        f[a] = mp->momentum_source[a];
        /*  Acoustic Reynolds Stress piece	*/
        f[a] += mp->u_momentum_source[0] * fv->grad_ars[a];
      }
      if (df != NULL && pd->v[pg->imtrx][ACOUS_REYN_STRESS]) {
        var = ACOUS_REYN_STRESS;
        phi = bf[var]->phi;
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          df->ars[a][j] = mp->u_momentum_source[0] * bf[var]->grad_phi[j][a];
        }
      }
    }
  } else {
    GOMA_EH(GOMA_ERROR, "No such Navier-Stokes Source Model");
  }

  if (ls != NULL && mp->mp2nd != NULL && mp->MomentumSourceModel != LEVEL_SET &&
      mp->mp2nd->MomentumSourceModel == CONSTANT && (pd->e[pg->imtrx][R_MOMENTUM1] & T_SOURCE)) {
    /* kludge for solidification tracking with phase function 0 */
    if (pfd != NULL && pd->e[pg->imtrx][R_EXT_VELOCITY]) {
      ls_old = ls;
      ls = pfd->ls[0];
      ls_modulate_momentumsource(f, mp->mp2nd->momentumsource_phase[0], ls->Length_Scale,
                                 (double)mp->mp2nd->momentumsourcemask[0],
                                 (double)mp->mp2nd->momentumsourcemask[1], df);
      ls = ls_old;
    }
    ls_modulate_momentumsource(f, mp->mp2nd->momentumsource, ls->Length_Scale,
                               (double)mp->mp2nd->momentumsourcemask[0],
                               (double)mp->mp2nd->momentumsourcemask[1], df);
  }

  return (status);
}
