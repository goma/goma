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

#include "mm_fill_continuity.h"
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

/* assemble_continuity -- assemble Residual &| Jacobian for continuity eqns
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
 *	a  -- gets loaded up with proper contribution
 * 	r  -- residual RHS vector
 *
 * Created:	Wed Mar  2 09:27:30 MST 1994 pasacki@sandia.gov
 *
 * Revised:	Sun Mar 20 13:24:50 MST 1994 pasacki@sandia.gov
 */
int assemble_continuity(dbl time_value, /* current time */
                        dbl tt,         /* parameter to vary time integration from
                                           explicit (tt = 1) to implicit (tt = 0)    */
                        dbl dt,         /* current time step size                    */
                        const PG_DATA *pg_data) {
  int dim;
  int p, q, a, b;

  int eqn, var;
  int peqn, pvar;
  int w;

  int i, j;
  int status, err;

  dbl time = 0.0; /*  RSL 6/6/02  */

  dbl *v = fv->v;        /* Velocity field. */
  dbl div_v = fv->div_v; /* Divergence of v. */

  dbl epsilon = 0.0, derivative, sum; /*  RSL 7/24/00  */
  dbl sum1, sum2;                     /*  RSL 8/15/00  */
  dbl sum_a, sum_b;                   /*  RSL 9/28/01  */
  int jj;                             /*  RSL 7/25/00  */

  dbl advection;
  dbl source;
  dbl pressure_stabilization;

  dbl volsolvent = 0;         /* volume fraction of solvent                */
  dbl initial_volsolvent = 0; /* initial solvent volume fraction
                               * (in stress-free state) input as source
                               * constant from input file                  */

  dbl det_J;
  dbl h3;
  dbl wt;
  dbl d_area;

  dbl d_h3detJ_dmesh_bj; /* for specific (b,j) mesh dof */

  /*
   * Galerkin weighting functions...
   */

  dbl phi_i;
  dbl phi_j;
  dbl div_phi_j_e_b;
  dbl(*grad_phi)[DIM]; /* weight-function for PSPG term */

  dbl div_v_dmesh; /* for specific (b,j) mesh dof */

  /*
   * Variables for Pressure Stabilization Petrov-Galerkin...
   */
  int meqn;
  int v_s[MAX_MODES][DIM][DIM], v_g[DIM][DIM];
  int mode;

  int *pdv = pd->v[pg->imtrx];

  dbl pspg[DIM];
  PSPG_DEPENDENCE_STRUCT d_pspg_struct;
  PSPG_DEPENDENCE_STRUCT *d_pspg = &d_pspg_struct;

  dbl mass, mass_a;
  dbl source_a;
  dbl sourceBase = 0.0;

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
  double dFVS_drho[MDE];
  double dFVS_dMOM[MAX_MOMENTS][MDE];

  int transient_run = pd->TimeIntegration != STEADY;
  int advection_on = 0;
  int source_on = 0;
  int ion_reactions_on = 0, electrode_kinetics_on = 0;
  int lagrangian_mesh_motion = 0, total_ale_on = 0;
  int hydromassflux_on = 0, suspensionsource_on = 0;
  int foam_volume_source_on = 0;
  int total_ale_and_velo_off = 0;

  dbl advection_etm, source_etm;

  double *J;

  status = 0;

  /*
   * Unpack variables from structures for local convenience...
   */

  eqn = R_PRESSURE;
  peqn = upd->ep[pg->imtrx][eqn];

  /*
   * Bail out fast if there's nothing to do...
   */

  if (!pd->e[pg->imtrx][eqn]) {
    return (status);
  }

  dim = pd->Num_Dim;

  if (pd->gv[POLYMER_STRESS11]) {
    err = stress_eqn_pointer(v_s);

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

  wt = fv->wt;
  det_J = bf[eqn]->detJ; /* Really, ought to be mesh eqn. */
  h3 = fv->h3;           /* Differential volume element (scales). */

  d_area = wt * det_J * h3;

  grad_phi = bf[eqn]->grad_phi;

  /*
   * Get the deformation gradients and tensors if needed
   */

  lagrangian_mesh_motion = (cr->MeshMotion == LAGRANGIAN || cr->MeshMotion == DYNAMIC_LAGRANGIAN);
  electrode_kinetics_on = (mp->SpeciesSourceModel[0] == ELECTRODE_KINETICS);
  ion_reactions_on = (mp->SpeciesSourceModel[0] == ION_REACTIONS);
  total_ale_on = (cr->MeshMotion == TOTAL_ALE);
  hydromassflux_on = (cr->MassFluxModel == HYDRODYNAMIC);
  suspensionsource_on = (mp->MomentumSourceModel == SUSPENSION);

  if (lagrangian_mesh_motion && pd->gv[R_MESH1]) {
    err = belly_flop(elc->lame_mu);
    GOMA_EH(err, "error in belly flop");
    if (err == 2)
      return (err);
  }

  if (total_ale_on && !pd->gv[VELOCITY1]) {
    total_ale_and_velo_off = 1;
  }

  if (total_ale_and_velo_off && pd->gv[R_SOLID1]) {
    err = belly_flop_rs(elc_rs->lame_mu);
    GOMA_EH(err, "error in belly flop for real solid");
    if (err == 2)
      return (err);
  }

  particle_momentum_on = 0;
  species = -1;
  ompvf = 1.0;

  if (pd->gv[R_PMOMENTUM1]) {
    particle_momentum_on = 1;
    species = (int)mp->u_density[0];
    ompvf = 1.0 - fv->c[species];
  }

  if (PSPG) {
    calc_pspg(pspg, d_pspg, time_value, tt, dt, pg_data);
  }

  if ((lagrangian_mesh_motion || total_ale_and_velo_off) && (mp->PorousMediaType == CONTINUOUS)) {
    initial_volsolvent = elc->Strss_fr_sol_vol_frac;
    volsolvent = 0.;
    for (w = 0; w < pd->Num_Species_Eqn; w++)
      volsolvent += fv->c[w];
    if (particle_momentum_on)
      volsolvent -= fv->c[species];
  }

  if (electrode_kinetics_on || ion_reactions_on) {
    if (mp->PorosityModel == CONSTANT) {
      epsilon = mp->porosity;
    } else if (mp->PorosityModel == THERMAL_BATTERY) {
      epsilon = mp->u_porosity[0];
    } else {
      GOMA_EH(GOMA_ERROR, "invalid porosity model");
    }
  }

  if (mp->MomentumSourceModel == SUSPENSION_PM || electrode_kinetics_on) /*  RSL 7/25/00  */
  {
    err = get_continuous_species_terms(&s_terms, 0.0, tt, dt, pg_data->hsquared);
    GOMA_EH(err, "problem in getting the species terms");
  }

  if (ion_reactions_on) /*  RSL 3/19/01 and 6/6/02  */
  {
    zero_structure(&s_terms, sizeof(struct Species_Conservation_Terms), 1);
    err = get_continuous_species_terms(&s_terms, time, tt, dt, pg_data->hsquared);
    GOMA_EH(err, "problem in getting the species terms");
  }

  if ((hydromassflux_on) && (mp->DensityModel == SUSPENSION) && (suspensionsource_on)) {
    /*
     * Compute hydrodynamic/sedimentation flux and sensitivities.
     */

    w0 =
        (int)mp->u_density[0]; /* This is the species number that is transported HYDRODYNAMICally */

    hydro_flux(&s_terms, w0, tt, dt, pg_data->hsquared);

    rhof = mp->u_density[1];
    rhos = mp->u_density[2];
  }
  rho = density(d_rho, time_value);

  advection_on = pd->e[pg->imtrx][eqn] & T_ADVECTION;
  source_on = pd->e[pg->imtrx][eqn] & T_SOURCE;

  advection_etm = pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
  source_etm = pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];

  dbl ls_disable_pspg = 1;
  //  if (ls != NULL && (fabs(fv->F) < ls->Length_Scale)) {
  //    ls_disable_pspg = 0;
  //  }

  if (af->Assemble_Residual) {
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
       *  Mass Terms: drhodt terms (usually though problem dependent)
       */
      mass = 0.0;
      if (particle_momentum_on) {
        if (transient_run) {
          mass = -s_terms.Y_dot[species];
          mass *= phi_i * d_area;
        }
      }

      if (electrode_kinetics_on || ion_reactions_on) {
        mass = 0.0;
        if (transient_run) {
          var = MASS_FRACTION;
          if (pd->gv[var]) {
            for (w = 0; w < pd->Num_Species - 1; w++) {
              derivative = 0.0;
              for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                if (bf[var]->phi[j] > 0.0)
                  break;
              }
              derivative = d_rho->C[w][j] / bf[var]->phi[j];
              mass += derivative * s_terms.Y_dot[w];
            }
          }
          mass *= epsilon / rho;
          mass *= phi_i * d_area;
        }
      }

      /*
       *  Advection:
       *    This term refers to the standard del dot v .
       *
       *    int (phi_i div_v d_omega)
       *
       *   Note density is not multiplied into this term normally
       */
      advection = 0.0;
      if (advection_on) {
        if (pd->gv[VELOCITY1]) /* then must be solving fluid mechanics in this material */
        {

          /*
           * Standard incompressibility constraint means we have
           * a solenoidal velocity field
           */

          advection = div_v;

          /* We get a more complicated advection term because the
           * particle phase is not incompressible.
           */
          if (particle_momentum_on)
            advection *= ompvf;

          advection *= phi_i * d_area;
          advection *= advection_etm;
        } else if (lagrangian_mesh_motion || total_ale_on)
        /* use divergence of displacement for linear elasticity */
        {
          advection = fv->volume_change;

          if (particle_momentum_on)
            advection *= ompvf;

          advection *= phi_i * h3 * det_J * wt;
          advection *= advection_etm;
        }

        if (electrode_kinetics_on || ion_reactions_on) {
          advection = div_v;
          var = MASS_FRACTION;
          if (pd->gv[var]) {
            for (w = 0; w < pd->Num_Species - 1; w++) {
              derivative = 0.0;
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
          }
          advection *= phi_i * d_area;
          advection *= advection_etm;
        }
      }

      source = 0.0;
      sourceBase = 0.0;
      if (source_on) {
        if (pd->gv[VELOCITY1]) {
          /* DRN (07/13/05):
             This was previously:
             source     =  P;
             But this messes with level set problems that have a density source
             over part of the domain and constant in other regions.  If someone was
             counting on this behavior as a form of a penalty method to give a non-zero
             diagonal entry, we should implement a new density model that accomplishes this.
             I really don't know what you want for DENSITY_IDEAL_GAS, though?!?

             source     =  0.;
             }*/
          if (mp->DensityModel == DENSITY_FOAM || mp->DensityModel == DENSITY_FOAM_CONC ||
              mp->DensityModel == DENSITY_FOAM_TIME || mp->DensityModel == DENSITY_FOAM_TIME_TEMP ||
              mp->DensityModel == DENSITY_MOMENT_BASED ||
              mp->DensityModel == DENSITY_FOAM_PMDI_10) {
            /* These density models locally permit a time and spatially varying
               density.  Consequently, the Lagrangian derivative of the density
               terms in the continuity equation are not zero and are
               included here as a source term
            */
            source =
                FoamVolumeSource(time_value, dt, tt, dFVS_dv, dFVS_dT, dFVS_dx, dFVS_dC, dFVS_dF);
            sourceBase = source;
            foam_volume_source_on = 1;
          } else if (mp->DensityModel == REACTIVE_FOAM) {
            /* These density models locally permit a time and spatially varying
               density.  Consequently, the Lagrangian derivative of the density
               terms in the continuity equation are not zero and are
               included here as a source term
            */
            source = REFVolumeSource(time_value, dt, tt, dFVS_dv, dFVS_dT, dFVS_dx, dFVS_dC);
            sourceBase = source;
            foam_volume_source_on = 1;
          } else if (mp->DensityModel == DENSITY_FOAM_PBE) {
            memset(dFVS_dv, 0, sizeof(double) * DIM * MDE);
            memset(dFVS_dT, 0, sizeof(double) * MDE);
            memset(dFVS_dx, 0, sizeof(double) * DIM * MDE);
            memset(dFVS_dMOM, 0, sizeof(double) * MAX_MOMENTS * MDE);
            memset(dFVS_dF, 0, sizeof(double) * MDE);
            memset(dFVS_drho, 0, sizeof(double) * MDE);
            memset(dFVS_dC, 0, sizeof(double) * MAX_CONC * MDE);
            source =
                PBEVolumeSource(time_value, dt, tt, dFVS_dv, dFVS_dT, dFVS_dx, dFVS_dC, dFVS_dMOM);
            foam_volume_source_on = 1;
          } else if (mp->DensityModel == DENSITY_FOAM_PBE_EQN) {
            memset(dFVS_dv, 0, sizeof(double) * DIM * MDE);
            memset(dFVS_dT, 0, sizeof(double) * MDE);
            memset(dFVS_dx, 0, sizeof(double) * DIM * MDE);
            memset(dFVS_dMOM, 0, sizeof(double) * MAX_MOMENTS * MDE);
            memset(dFVS_dF, 0, sizeof(double) * MDE);
            memset(dFVS_drho, 0, sizeof(double) * MDE);
            memset(dFVS_dC, 0, sizeof(double) * MAX_CONC * MDE);
            source = PBEVolumeSource_rhoeqn(time_value, dt, tt, dFVS_drho);
            foam_volume_source_on = 1;
          }

          /*
            else if
            ( mp->DensityModel == SUSPENSION ||
            mp->DensityModel == SUSPENSION_PM )
            {
          */
          /* Although the suspension density models meet the definition of
             a locally variable density model, the Lagrangian derivative
             of their densities can be represented as a divergence of
             mass flux.  This term must be integrated by parts and so is
             included separately later on and is not include as part of the "source"
             terms
             source = 0.0;
             }  */

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
            sourceBase = source;
          }
          source *= phi_i * d_area;
          source *= source_etm;
        }

        if ((lagrangian_mesh_motion || total_ale_and_velo_off))
        /* add swelling as a source of volume */
        {
          if (mp->PorousMediaType == CONTINUOUS) {
            source = -(1. - initial_volsolvent) / (1. - volsolvent);
            sourceBase = source;
            source *= phi_i * d_area;
            source *= source_etm;
          }
        }
        if (electrode_kinetics_on || ion_reactions_on) {
          source = 0.0;
          for (j = 0; j < pd->Num_Species; j++) {
            source -= s_terms.MassSource[j] * mp->molecular_weight[j];
          }
          source /= rho;
          sourceBase = source;
          source *= phi_i * d_area;
          source *= source_etm;
        }
      }
      /* add Pressure-Stabilized Petrov-Galerkin term
       * if desired.
       */

      pressure_stabilization = 0.0;
      if (PSPG) {
        for (a = 0; a < WIM; a++) {
          meqn = R_MOMENTUM1 + a;
          if (pd->gv[meqn]) {
            pressure_stabilization += grad_phi[i][a] * pspg[a];
          }
        }
        pressure_stabilization *= d_area * ls_disable_pspg;
      }

      h_flux = 0.0;
      if ((hydromassflux_on) && (suspensionsource_on)) {
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

#ifdef DEBUG_CONTINUITY_RES
      printf("R_c[%d] += %10f + %10f + %10f + %10f + %10f\n", i, mass, advection, source,
             pressure_stabilization, h_flux);
#endif
      /*
       *  Add up the individual contributions and sum them into the local element
       *  contribution for the total continuity equation for the ith local unknown
       */
      lec->R[LEC_R_INDEX(peqn, i)] += mass + advection + source + pressure_stabilization + h_flux;
    }
  }

  if (af->Assemble_Jacobian) {
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
       * J_c_v NOTE that this is applied whenever velocity is a variable
       */
      for (b = 0; b < WIM; b++) {
        var = VELOCITY1 + b;
        if (pdv[var]) {
          pvar = upd->vp[pg->imtrx][var];

          J = &(lec->J[LEC_J_INDEX(peqn, pvar, i, 0)]);

          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];

            h_flux = 0.;
            if ((hydromassflux_on) && (suspensionsource_on)) {
              for (p = 0; p < dim; p++) {
                h_flux += grad_phi[i][p] * s_terms.d_diff_flux_dv[w0][p][b][j];
              }
              h_flux *= (rhos - rhof) / rhof * d_area * advection_etm;
            }

            advection = 0.;

            if (advection_on) {
              div_phi_j_e_b = 0.;
#ifdef DO_NO_UNROLL
              for (p = 0; p < VIM; p++) {
                div_phi_j_e_b += bf[var]->grad_phi_e[j][b][p][p];
              }
#else
              div_phi_j_e_b += bf[var]->grad_phi_e[j][b][0][0] + bf[var]->grad_phi_e[j][b][1][1];
              if (VIM == 3)
                div_phi_j_e_b += bf[var]->grad_phi_e[j][b][2][2];
#endif

              if (electrode_kinetics_on || ion_reactions_on) {
                sum = 0.;
                for (jj = 0; jj < pd->Num_Species - 1; jj++) {
                  derivative = 0.0;
                  for (q = 0; q < ei[pg->imtrx]->dof[MASS_FRACTION]; q++) {
                    if (bf[MASS_FRACTION]->phi[q] > 0.0)
                      break;
                  }
                  derivative = d_rho->C[jj][q] / bf[MASS_FRACTION]->phi[q];
                  sum += derivative * s_terms.d_conv_flux_dv[jj][b][b][j];
                }
                div_phi_j_e_b += sum / rho;
              }

              advection = phi_i * div_phi_j_e_b * d_area;

              if (particle_momentum_on)
                advection *= ompvf;

              advection *= advection_etm;
            }

            source = 0.;

            if (source_on) {

              if (foam_volume_source_on) {
                source = dFVS_dv[b][j];
              }

              if (particle_momentum_on) {
                /* From residual calculation.
                   source -= fv->grad_c[species][a] * v[a];
                */
                source = -s_terms.grad_Y[species][b] * phi_j;
              }

              source *= phi_i * d_area;
              source *= source_etm;
            }

            /* add Pressure-Stabilized Petrov-Galerkin term
             * if desired.
             */
            pressure_stabilization = 0.;
            if (PSPG) {
              for (a = 0; a < WIM; a++) {
                meqn = R_MOMENTUM1 + a;
                if (pd->e[pg->imtrx][meqn]) {
                  pressure_stabilization += grad_phi[i][a] * d_pspg->v[a][b][j];
                }
              }
              pressure_stabilization *= d_area * ls_disable_pspg;
            }

            J[j] += advection + source + pressure_stabilization + h_flux;
            /* lec->J[LEC_J_INDEX(peqn,pvar,i,j)] += advection + source + pressure_stabilization +
             * h_flux; */
#ifdef DEBUG_CONTINUITY_JAC
            printf("J_c_v[%d] [%d][%d] += %10f + %10f + %10f + %10f\n", i, b, j, advection, source,
                   pressure_stabilization, h_flux);
#endif
          }
        }
      }

      var = EDDY_NU;
      if (PSPG && pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          pvar = upd->vp[pg->imtrx][var];

          phi_j = bf[var]->phi[j];

          /* add Pressure-Stabilized Petrov-Galerkin term
           * if desired.
           */
          pressure_stabilization = 0.;

          for (a = 0; a < WIM; a++) {
            meqn = R_MOMENTUM1 + a;
            if (pd->e[pg->imtrx][meqn]) {
              pressure_stabilization += grad_phi[i][a] * d_pspg->eddy_nu[a][j];
            }
          }
          pressure_stabilization *= d_area * ls_disable_pspg;

          lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += pressure_stabilization;
        }
      }

      /*
       * J_c_T This term comes from the temperature dependency of the momentume source
       * which comes from the Pressure-Stabilized Petrov-Galerkin term
       */
      var = TEMPERATURE;
      if (PSPG && pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          pvar = upd->vp[pg->imtrx][var];

          phi_j = bf[var]->phi[j];

          /* add Pressure-Stabilized Petrov-Galerkin term
           * if desired.
           */
          pressure_stabilization = 0.;

          for (a = 0; a < WIM; a++) {
            meqn = R_MOMENTUM1 + a;
            if (pd->e[pg->imtrx][meqn]) {
              pressure_stabilization += grad_phi[i][a] * d_pspg->T[a][j];
            }
          }
          pressure_stabilization *= d_area * ls_disable_pspg;

          lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += pressure_stabilization;
        }
      }

      if (source_on) {
        if (foam_volume_source_on) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            pvar = upd->vp[pg->imtrx][var];

            lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += d_area * phi_i * dFVS_dT[j] * source_etm;
          }
        }
      }

      if ((hydromassflux_on) && pdv[var]) {
        if (suspensionsource_on) {
          pvar = upd->vp[pg->imtrx][var];

          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            h_flux = 0.;

            for (p = 0; p < dim; p++) {
              h_flux += grad_phi[i][p] * s_terms.d_diff_flux_dT[w0][p][j];
            }

            h_flux *= d_area * (rhos - rhof) / rhof;

            /*  h_flux = 0.0; */

            lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += h_flux;
          }
        }
      }

      /*
       *  J_c_F this is primarily the foam volume source terms
       */

      var = FILL;
      if (source_on && pdv[var] && ls != NULL) {
        pvar = upd->vp[pg->imtrx][var];

        J = &(lec->J[LEC_J_INDEX(peqn, pvar, i, 0)]);

        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          source = 0.0;

          if (foam_volume_source_on) {
            source = dFVS_dF[j];
            source *= phi_i * d_area;
          }

          J[j] += source;

          /*lec->J[LEC_J_INDEX(peqn,pvar,i,j)] += source;*/
        }
      }

      /*
       * J_c_P here species act as a volume source in continuous lagrangian mesh motion
       */
      var = PRESSURE;
      if (pdv[var]) {
        pvar = upd->vp[pg->imtrx][var];

        J = &(lec->J[LEC_J_INDEX(peqn, pvar, i, 0)]);

        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          phi_j = bf[var]->phi[j];

          advection = 0.;

          if (advection_on && (lagrangian_mesh_motion || total_ale_and_velo_off)) {
            /*Need to compute this for total ALE.  Not done yet */
            advection = fv->d_volume_change_dp[j];

            advection *= phi_i * d_area;

            advection *= advection_etm;
          }

          source = 0.;

          /* add Pressure-Stabilized Petrov-Galerkin term
           * if desired.
           */
          pressure_stabilization = 0.;
          if (PSPG) {
            for (a = 0; a < WIM; a++) {
              meqn = R_MOMENTUM1 + a;
              if (pd->e[pg->imtrx][meqn] & T_DIFFUSION) {
                pressure_stabilization += grad_phi[i][a] * d_pspg->P[a][j];
              }
            }
            pressure_stabilization *= d_area * ls_disable_pspg;
          }
          J[j] += advection + source + pressure_stabilization;
          /*lec->J[LEC_J_INDEX(peqn,pvar,i,j)] += advection  + source + pressure_stabilization; */
#ifdef DEBUG_CONTINUITY_JAC
          printf("J_c_P[%d] [%d] += %10f %10f %10f\n", i, j, advection, source,
                 pressure_stabilization);
#endif
        }
      }

      /*
       * J_c_S this term is only present for PSPG
       */
      var = POLYMER_STRESS11;
      if (PSPG && pdv[var]) {
        for (mode = 0; mode < vn->modes; mode++) {
          for (p = 0; p < VIM; p++) {
            for (q = 0; q < VIM; q++) {
              var = v_s[mode][p][q];
              pvar = upd->vp[pg->imtrx][var];
              if (pd->v[pg->imtrx][var]) {
                for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                  pressure_stabilization = 0.;

                  for (a = 0; a < WIM; a++) {
                    meqn = R_MOMENTUM1 + a;
                    if (pd->e[pg->imtrx][meqn]) {
                      pressure_stabilization += grad_phi[i][a] * d_pspg->S[a][mode][p][q][j];
                    }
                  }

                  pressure_stabilization *= h3 * det_J * wt * ls_disable_pspg;

                  lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += pressure_stabilization;
                }
              }
            }
          }
        }
      }

      /*
       * J_c_G this term is only present for PSPG
       */
      var = VELOCITY_GRADIENT11;
      if (PSPG && pdv[var]) {
        for (p = 0; p < VIM; p++) {
          for (q = 0; q < VIM; q++) {
            var = v_g[p][q];
            pvar = upd->vp[pg->imtrx][var];
            {
              for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                pressure_stabilization = 0.;

                for (a = 0; a < WIM; a++) {
                  meqn = R_MOMENTUM1 + a;
                  if (pd->e[pg->imtrx][meqn]) {
                    pressure_stabilization += grad_phi[i][a] * d_pspg->g[a][p][q][j];
                  }
                }
                pressure_stabilization *= h3 * det_J * wt * ls_disable_pspg;

                lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += pressure_stabilization;
              }
            }
          }
        }
      }

      /*
       * J_c_SH this term is only present for HYDRODYNAMIC mass flux and SUSPENSION
       *        momentum source
       */

      var = SHEAR_RATE;

      if (hydromassflux_on && pdv[var]) {
        if (suspensionsource_on) {
          if (pd->v[pg->imtrx][var]) {
            pvar = upd->vp[pg->imtrx][var];

            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              h_flux = 0.0;
              for (a = 0; a < dim; a++) {
                h_flux += grad_phi[i][a] * s_terms.d_diff_flux_dSH[w0][a][j];
              }

              h_flux *= h3 * det_J * wt * (rhos - rhof) / rhof;

              /*  h_flux = 0.0; */

              lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += h_flux;
            }
          }
        }
      }

      /*
       * J_c_d
       */

      for (b = 0; b < dim; b++) {
        var = MESH_DISPLACEMENT1 + b;
        if (pdv[var]) {
          pvar = upd->vp[pg->imtrx][var];
          J = &(lec->J[LEC_J_INDEX(peqn, pvar, i, 0)]);
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];
            /* derivative of |J| with extra term for axisymmetry e.g.
               d/dmesh [ r|J| ] */
            d_h3detJ_dmesh_bj =
                (h3 * bf[eqn]->d_det_J_dm[b][j] + det_J * fv->dh3dq[b] * bf[var]->phi[j]);

            mass = 0.0;
            if (electrode_kinetics_on || ion_reactions_on) {
              if (transient_run) {
                for (w = 0; w < pd->Num_Species - 1; w++) {
                  derivative = 0.0;
                  for (q = 0; q < ei[pg->imtrx]->dof[MASS_FRACTION]; q++) {
                    if (bf[MASS_FRACTION]->phi[q] > 0.0)
                      break;
                  }
                  derivative = d_rho->C[w][q] / bf[MASS_FRACTION]->phi[q];
                  mass += derivative * s_terms.Y_dot[w];
                }
                mass *= epsilon / rho;
                mass *= phi_i * d_h3detJ_dmesh_bj * wt;
              }
            }

            advection = 0.0;
            if (advection_on) {
              if (pdv[VELOCITY1]) {
                h_flux = 0.0;
                if ((hydromassflux_on) && (suspensionsource_on)) {
                  for (p = 0; p < dim; p++) {
                    h_flux += grad_phi[i][p] * s_terms.diff_flux[w0][p] * d_h3detJ_dmesh_bj +
                              grad_phi[i][p] * s_terms.d_diff_flux_dmesh[w0][p][b][j] * det_J * h3 +
                              bf[eqn]->d_grad_phi_dmesh[i][p][b][j] * s_terms.diff_flux[w0][p] *
                                  det_J * h3;
                  }
                  h_flux *= (rhos - rhof) / rhof * wt * pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
                }

                div_v_dmesh = fv->d_div_v_dmesh[b][j];

                advection += div_v_dmesh * det_J * h3 + div_v * (d_h3detJ_dmesh_bj);

                if (electrode_kinetics_on || ion_reactions_on) /*  RSL  9/28/01  */
                {
                  sum_a = 0.;
                  sum_b = 0.;
                  for (w = 0; w < pd->Num_Species - 1; w++) {
                    derivative = 0.0;
                    for (q = 0; q < ei[pg->imtrx]->dof[MASS_FRACTION]; q++) {
                      if (bf[MASS_FRACTION]->phi[q] > 0.0)
                        break;
                    }
                    derivative = d_rho->C[w][q] / bf[MASS_FRACTION]->phi[q];
                    sum1 = 0.;
                    sum2 = 0.;
                    for (p = 0; p < dim; p++) {
                      sum1 += s_terms.conv_flux[w][p];
                      sum2 += s_terms.d_conv_flux_dmesh[w][p][b][j];
                    }
                    sum_a += derivative * sum1;
                    sum_b += derivative * sum2;
                  }
                  sum_a /= rho;
                  sum_b /= rho;
                  advection += sum_b * det_J * h3 + sum_a * d_h3detJ_dmesh_bj;
                }

                advection *= phi_i * wt;
              } else if (lagrangian_mesh_motion || total_ale_and_velo_off) {
                advection += fv->volume_change * (d_h3detJ_dmesh_bj);

                advection += fv->d_volume_change_dx[b][j] * h3 * det_J;

                advection *= phi_i * wt;
              }
              advection *= advection_etm;
            }

            source = 0.0;
            if (source_on) {
              if (mp->DensityModel == DENSITY_FOAM || mp->DensityModel == DENSITY_FOAM_CONC ||
                  mp->DensityModel == DENSITY_FOAM_TIME || mp->DensityModel == DENSITY_FOAM_PBE ||
                  mp->DensityModel == DENSITY_FOAM_TIME_TEMP ||
                  mp->DensityModel == DENSITY_FOAM_PMDI_10) {
                source = sourceBase * d_h3detJ_dmesh_bj * wt + dFVS_dx[b][j] * d_area;
                source *= phi_i * source_etm;
              } else if (mp->DensityModel == REACTIVE_FOAM) {
                source = sourceBase * d_h3detJ_dmesh_bj * wt + dFVS_dx[b][j] * d_area;
                source *= phi_i * source_etm;
              }
            }
            if (lagrangian_mesh_motion || total_ale_and_velo_off) {
              /* add swelling as a source of volume */
              if (mp->PorousMediaType == CONTINUOUS) {
                source = -phi_i * (d_h3detJ_dmesh_bj)*wt * (1. - initial_volsolvent) /
                         (1. - volsolvent) * source_etm;
              }
            }

            if (electrode_kinetics_on || ion_reactions_on) {
              sum1 = 0.;
              sum2 = 0.;
              for (q = 0; q < pd->Num_Species; q++) {
                sum1 -= s_terms.MassSource[q] * mp->molecular_weight[q];
                sum2 -= s_terms.d_MassSource_dmesh[q][b][j] * mp->molecular_weight[q];
              }
              sum1 /= rho;
              sum2 /= rho;
              source += sum2 * det_J * h3 + sum1 * d_h3detJ_dmesh_bj;
              source *= phi_i * wt;
              source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
            }

            /* add Pressure-Stabilized Petrov-Galerkin term
             * if desired.
             */
            pressure_stabilization = 0.0;
            if (PSPG) {
              for (a = 0; a < WIM; a++) {
                meqn = R_MOMENTUM1 + a;
                if (pd->e[pg->imtrx][meqn]) {
                  pressure_stabilization +=
                      grad_phi[i][a] * d_pspg->X[a][b][j] * h3 * det_J * wt +
                      grad_phi[i][a] * pspg[a] * wt * d_h3detJ_dmesh_bj +
                      bf[eqn]->d_grad_phi_dmesh[i][a][b][j] * pspg[a] * wt * h3 * det_J;
                }
              }
              pressure_stabilization *= ls_disable_pspg;
            }

            /*lec->J[LEC_J_INDEX(peqn,pvar,i,j)] += advection  + source + pressure_stabilization +
             * h_flux + mass ; */
            J[j] += advection + source + pressure_stabilization + h_flux + mass;
          }
        }
      }

      /*
       * J_c_d_rs
       */

      for (b = 0; b < dim; b++) {
        var = SOLID_DISPLACEMENT1 + b;
        if (pdv[var]) {
          pvar = upd->vp[pg->imtrx][var];

          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];

            advection = 0.;

            if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
              if (total_ale_and_velo_off) {
                advection += fv->d_volume_change_drs[b][j] * h3 * det_J;

                advection *= phi_i * wt;
              }

              advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
            }

            source = 0.;

            lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += advection + source;
          }
        }
      }

      /*
       * J_c_c
       */

      var = MASS_FRACTION;
      if (pdv[var]) {
        pvar = upd->vp[pg->imtrx][var];

        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          phi_j = bf[var]->phi[j];

          for (w = 0; w < pd->Num_Species_Eqn; w++) {
            /* add swelling as a source of volume */
            source = 0.;
            if ((mp->PorousMediaType == CONTINUOUS) &&
                (cr->MeshMotion == LAGRANGIAN || cr->MeshMotion == DYNAMIC_LAGRANGIAN)) {
              source = -phi_j * phi_i * h3 * det_J * wt * (1. - initial_volsolvent) /
                       (1. - volsolvent) / (1. - volsolvent) * pd->etm[pg->imtrx][eqn][LOG2_SOURCE];
            }
            if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
              /* Foaming volume source term */

              if (mp->DensityModel == REACTIVE_FOAM || mp->DensityModel == DENSITY_FOAM ||
                  mp->DensityModel == DENSITY_FOAM_PBE) {
                source += dFVS_dC[w][j];
                source *= phi_i * h3 * det_J * wt * pd->etm[pg->imtrx][eqn][LOG2_SOURCE];
              }
            }

            /* add Pressure-Stabilized Petrov-Galerkin term
             * if desired.
             */
            pressure_stabilization = 0.;
            if (PSPG) {
              for (a = 0; a < WIM; a++) {
                meqn = R_MOMENTUM1 + a;
                if (pd->e[pg->imtrx][meqn] & T_DIFFUSION) {
                  pressure_stabilization += grad_phi[i][a] * d_pspg->C[a][w][j];
                }
              }
              pressure_stabilization *= h3 * det_J * wt * ls_disable_pspg;
            }

            /* The fluid phase is not incompressible in the
             * SUSPENION_PM model.
             */
            mass_a = 0.0;
            advection = 0.0;
            if (particle_momentum_on && w == species) {
              mass_a = -(1.0 + 2.0 * tt) * phi_j / dt;
              mass_a *= phi_i * h3 * det_J * wt;

              if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
                advection = -phi_j * div_v;
                advection *= phi_i * det_J * h3 * wt;
              }

              if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
                source_a = 0.0;
                for (a = 0; a < WIM; a++)
                  source_a -= grad_phi[j][a] * v[a];
                source_a *= phi_i * det_J * h3 * wt;

                source += source_a;
              }
            }

            if (mp->SpeciesSourceModel[0] == ELECTRODE_KINETICS ||
                mp->SpeciesSourceModel[0] == ION_REACTIONS) /*  RSL 3/19/01  */
            {
              sum = 0.;
              for (jj = 0; jj < pd->Num_Species - 1; jj++) {
                derivative = 0.0;
                for (q = 0; q < ei[pg->imtrx]->dof[MASS_FRACTION]; q++) {
                  if (bf[MASS_FRACTION]->phi[q] > 0.0)
                    break;
                }
                derivative = d_rho->C[jj][q] / bf[MASS_FRACTION]->phi[q];
                sum += derivative *
                       (s_terms.d_Y_dot_dc[jj][w][j] - d_rho->C[w][j] * s_terms.Y_dot[jj] / rho);
              }
              mass_a = sum * epsilon / rho;
              mass_a *= phi_i * h3 * det_J * wt;
            }

            if (mp->SpeciesSourceModel[0] == ELECTRODE_KINETICS ||
                mp->SpeciesSourceModel[0] == ION_REACTIONS) /*  RSL 3/19/01  */
            {
              sum = 0.;
              for (jj = 0; jj < pd->Num_Species - 1; jj++) {
                derivative = 0.0;
                for (q = 0; q < ei[pg->imtrx]->dof[MASS_FRACTION]; q++) {
                  if (bf[MASS_FRACTION]->phi[q] > 0.0)
                    break;
                }
                derivative = d_rho->C[jj][q] / bf[MASS_FRACTION]->phi[q];
                sum1 = 0.;
                sum2 = 0.;
                for (p = 0; p < VIM; p++) {
                  sum1 += s_terms.d_conv_flux_dc[jj][p][w][j];
                  sum2 += s_terms.conv_flux[jj][p];
                }
                sum += derivative * (sum1 - d_rho->C[w][j] * sum2 / rho);
              }
              advection = sum / rho;
              advection *= phi_i * h3 * det_J * wt;
            }

            if (mp->SpeciesSourceModel[0] == ELECTRODE_KINETICS ||
                mp->SpeciesSourceModel[0] == ION_REACTIONS) /*  RSL 3/19/01  */
            {
              source_a = 0.0;

              for (jj = 0; jj < pd->Num_Species; jj++) {
                source_a +=
                    mp->molecular_weight[jj] * (s_terms.d_MassSource_dc[jj][w][j] -
                                                d_rho->C[w][j] * s_terms.MassSource[jj] / rho);
              }
              source_a /= (-rho);
              source_a *= phi_i * h3 * det_J * wt;

              source += source_a;
            }

            lec->J[LEC_J_INDEX(peqn, MAX_PROB_VAR + w, i, j)] +=
                mass_a + advection + source + pressure_stabilization;
          }
        }
      }

      if (cr->MassFluxModel == HYDRODYNAMIC && pd->v[pg->imtrx][var]) {
        if (mp->MomentumSourceModel == SUSPENSION) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];

            h_flux = 0.0;

            for (a = 0; a < dim; a++) {
              h_flux += grad_phi[i][a] * s_terms.d_diff_flux_dc[w0][a][w0][j];
            }

            h_flux *= h3 * det_J * wt * (rhos - rhof) / rhof;
            /*  h_flux = 0.0; */

            lec->J[LEC_J_INDEX(peqn, MAX_PROB_VAR + w0, i, j)] += h_flux;
          }
        }
      }

      /*
       * J_c_MOM
       */
      for (b = 0; b < WIM; b++) {
        var = MOMENT0 + b;
        if (pdv[var]) {
          pvar = upd->vp[pg->imtrx][var];

          J = &(lec->J[LEC_J_INDEX(peqn, pvar, i, 0)]);

          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];

            source = 0.;

            if (source_on) {

              if (foam_volume_source_on) {
                source = dFVS_dMOM[b][j];
              }

              source *= phi_i * d_area;
              source *= source_etm;
            }

            J[j] += source;
          }
        }
      }

      var = DENSITY_EQN;
      if (pdv[var]) {
        pvar = upd->vp[pg->imtrx][var];

        J = &(lec->J[LEC_J_INDEX(peqn, pvar, i, 0)]);

        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          phi_j = bf[var]->phi[j];

          source = 0.;

          if (source_on) {

            if (foam_volume_source_on) {
              source = dFVS_drho[j];
            }

            source *= phi_i * d_area;
            source *= source_etm;
          }

          J[j] += source;
        }
      }
    }
  }

  return (status);
} /*
   * These density models locally permit a time and spatially varying
   *  density.  Consequently, neither the gradient nor the time derivative
   *  of the density term in the continuity equation is zero. These terms
   *  must then be included here as a "source" term via the expresion:
   *
   *   0 = del dot V + FoamVolumeSource
   *
   *  Therefore, this implies that
   *
   *     FoamVolumeSource = 1/rho * drho/dt + 1/rho * (v dot del rho)
   *
   */
double FoamVolumeSource(double time,
                        double dt,
                        double tt,
                        double dFVS_dv[DIM][MDE],
                        double dFVS_dT[MDE],
                        double dFVS_dx[DIM][MDE],
                        double dFVS_dC[MAX_CONC][MDE],
                        double dFVS_dF[MDE]) {
  double DT_Dt, vol, rho, rho2, T, Press;
  double deriv_x, deriv_T;
  double source = 0.0, source_a, source_b, source_c;
  int var, j, a, dim, w, species;
  double d_rho_dt, drho_T;

  /* density derivatives */
  /*   DENSITY_DEPENDENCE_STRUCT d_rho_struct;  /\* density dependence *\/ */
  /*   DENSITY_DEPENDENCE_STRUCT *d_rho = &d_rho_struct; */

  dim = pd->Num_Dim;

  memset(dFVS_dv, 0, sizeof(double) * DIM * MDE);
  memset(dFVS_dT, 0, sizeof(double) * MDE);
  memset(dFVS_dx, 0, sizeof(double) * DIM * MDE);
  memset(dFVS_dF, 0, sizeof(double) * MDE);
  memset(dFVS_dC, 0, sizeof(double) * MDE * MAX_CONC);

  if (mp->DensityModel == DENSITY_FOAM) {
    dbl x0, MW, Rgas, rho_epoxy, rho_fluor, Dx_Dt, drho_x;
    species = (int)mp->u_density[0]; /* species number fluorinert */
    x0 = mp->u_density[1];           /* Initial fluorinert mass fraction */
    Rgas = mp->u_density[2];         /* Gas constant in appropriate units */
    MW = mp->u_density[3];
    rho_epoxy = mp->u_density[4]; /* Density of epoxy resin */
    rho_fluor = mp->u_density[5]; /* Density of liquid fluorinert */

    if (fv->c[species] > 0.)
      vol = fv->c[species];
    else
      vol = 0.;
    if (vol > x0)
      vol = x0;

    T = fv->T;
    Press = upd->Pressure_Datum;

    rho = 1. / ((x0 - vol) * Rgas * T / (Press * MW) + (1.0 - x0) / rho_epoxy + vol / rho_fluor);
    rho2 = rho * rho;
    /*   rho = density(d_rho); */

    Dx_Dt = fv_dot->c[species];
    for (a = 0; a < dim; a++) {
      Dx_Dt += fv->v[a] * fv->grad_c[species][a];
    }

    DT_Dt = fv_dot->T;
    for (a = 0; a < dim; a++) {
      DT_Dt += fv->v[a] * fv->grad_T[a];
    }
    deriv_x = (-Rgas * T / (Press * MW) + 1. / rho_fluor);
    deriv_T = (x0 - vol) * Rgas / (Press * MW);
    drho_x = -rho2 * deriv_x;
    drho_T = -rho2 * deriv_T;
    if ((vol > 0.) && (vol < x0))
      source = -rho * (deriv_x * Dx_Dt + deriv_T * DT_Dt);
    else
      source = 0.;

    if (af->Assemble_Jacobian && ((vol > 0.) && (vol < x0))) {

      var = TEMPERATURE;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          source_a = bf[var]->phi[j] * (1 + 2. * tt) / dt;
          for (a = 0; a < dim; a++) {
            source_a += fv->v[a] * bf[var]->grad_phi[j][a];
          }
          source_a *= -rho * deriv_T; /* derivative of convective terms */
          source_b = -drho_T * (deriv_x * Dx_Dt + deriv_T * DT_Dt) *
                     bf[var]->phi[j]; /* derivative of density */
          source_c = -rho * Rgas / (Press * MW) * Dx_Dt *
                     bf[var]->phi[j]; /* derivative of deriv_x wrt T */
          dFVS_dT[j] = source_a + source_b + source_c;
        }
      }

      var = VELOCITY1;
      if (pd->v[pg->imtrx][var]) {
        for (a = 0; a < dim; a++) {
          var = VELOCITY1 + a;
          for (j = 0; pd->v[pg->imtrx][var] && j < ei[pg->imtrx]->dof[var]; j++) {
            dFVS_dv[a][j] = -rho * (deriv_x * bf[var]->phi[j] * fv->grad_c[species][a] +
                                    deriv_T * bf[var]->phi[j] * fv->grad_T[a]);
          }
        }
      }

      var = MASS_FRACTION;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          source_a = bf[var]->phi[j] * (1 + 2. * tt) / dt;
          for (a = 0; a < dim; a++) {
            source_a += fv->v[a] * bf[var]->grad_phi[j][a];
          }
          source_a *= -rho * deriv_x; /* derivative of convective terms */
          source_b = -drho_x * (deriv_x * Dx_Dt + deriv_T * DT_Dt) *
                     bf[var]->phi[j]; /* derivative of density */
          source_c = -rho * Rgas / (Press * MW) * DT_Dt *
                     bf[var]->phi[j]; /* derivative of deriv_T wrt x */
          if ((vol > 0.) && (vol < x0))
            dFVS_dC[species][j] = source_a + source_b + source_c;
        }
      }
    }
  } else if (mp->DensityModel == DENSITY_FOAM_CONC) {
    double Rgas, MW_f, MW_a, rho_epoxy, rho_fluor, T;
    int species_l, species_v, species_a;
    double Dcl_Dt, Dcv_Dt, Dca_Dt, DT_Dt;
    /* density derivatives */
    double drho_T, rho_v_inv, d_rho_v_inv_dT, rho_a_inv, d_rho_a_inv_dT;
    double drho_c_v, drho_c_a, drho_c_l, d2rho_T_c_v, d2rho_T_c_a;
    double source_a_v, source_a_a, source_a_l, source_b_v, source_b_a, source_b_l, source_c_v,
        source_c_a;

    species_l = (int)mp->u_density[0]; /* species number fluorinert liquid */
    species_v = (int)mp->u_density[1]; /* species number fluorinert vapor */
    species_a = (int)mp->u_density[2]; /* species number air vapor */
    Rgas = mp->u_density[3];           /* Gas constant in appropriate units */
    MW_f = mp->u_density[4];           /* molecular weight fluorinert */
    MW_a = mp->u_density[5];           /* molecular weight air */
    rho_epoxy = mp->u_density[6];      /* Density of epoxy resin */
    rho_fluor = mp->u_density[7];      /* Density of liquid fluorinert */

    T = fv->T;
    Press = upd->Pressure_Datum;
    rho_v_inv = Rgas * T / (Press * MW_f);
    d_rho_v_inv_dT = Rgas / (Press * MW_f);
    rho_a_inv = Rgas * T / (Press * MW_a);
    d_rho_a_inv_dT = Rgas / (Press * MW_a);

    rho = rho_epoxy + fv->c[species_v] * (1. - rho_epoxy * rho_v_inv) +
          fv->c[species_a] * (1. - rho_epoxy * rho_a_inv) +
          fv->c[species_l] * (1. - rho_epoxy / rho_fluor);

    rho2 = rho * rho;

    drho_c_v = 1. - rho_epoxy * rho_v_inv;
    drho_c_a = 1. - rho_epoxy * rho_a_inv;
    drho_c_l = 1. - rho_epoxy / rho_fluor;

    drho_T = -fv->c[species_v] * rho_epoxy * d_rho_v_inv_dT -
             fv->c[species_a] * rho_epoxy * d_rho_a_inv_dT;

    d2rho_T_c_v = -rho_epoxy * d_rho_v_inv_dT;
    d2rho_T_c_a = -rho_epoxy * d_rho_a_inv_dT;

    Dcl_Dt = fv_dot->c[species_l];
    Dcv_Dt = fv_dot->c[species_v];
    Dca_Dt = fv_dot->c[species_a];
    for (a = 0; a < dim; a++) {
      Dcl_Dt += fv->v[a] * fv->grad_c[species_l][a];
      Dcv_Dt += fv->v[a] * fv->grad_c[species_v][a];
      Dca_Dt += fv->v[a] * fv->grad_c[species_a][a];
    }

    DT_Dt = fv_dot->T;
    for (a = 0; a < dim; a++) {
      DT_Dt += fv->v[a] * fv->grad_T[a];
    }

    source = (drho_c_v * Dcv_Dt + drho_c_a * Dca_Dt + drho_c_l * Dcl_Dt + drho_T * DT_Dt) / rho;

    if (af->Assemble_Jacobian) {
      var = TEMPERATURE;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          source_a = bf[var]->phi[j] * (1 + 2. * tt) / dt;
          for (a = 0; a < dim; a++) {
            source_a += fv->v[a] * bf[var]->grad_phi[j][a];
          }
          source_a *= drho_T / rho;                            /* derivative of convective terms */
          source_b = -drho_T * source * bf[var]->phi[j] / rho; /* derivative of 1/density */
          source_c = +(d2rho_T_c_v * Dcv_Dt + d2rho_T_c_a * Dca_Dt) * bf[var]->phi[j] /
                     rho; /* derivative of drho_c_v & drho_c_a wrt T */
          dFVS_dT[j] = source_a + source_b + source_c;
        }
      }

      var = VELOCITY1;
      if (pd->v[pg->imtrx][var]) {
        for (a = 0; a < dim; a++) {
          var = VELOCITY1 + a;
          for (j = 0; pd->v[pg->imtrx][var] && j < ei[pg->imtrx]->dof[var]; j++) {
            dFVS_dv[a][j] =
                +bf[var]->phi[j] *
                (drho_c_v * fv->grad_c[species_v][a] + drho_c_a * fv->grad_c[species_a][a] +
                 drho_c_l * fv->grad_c[species_l][a] + drho_T * fv->grad_T[a]) /
                rho;
          }
        }
      }

      var = MASS_FRACTION;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          source_a = bf[var]->phi[j] * (1 + 2. * tt) / dt;
          for (a = 0; a < dim; a++) {
            source_a += fv->v[a] * bf[var]->grad_phi[j][a];
          }
          source_a_v = source_a * drho_c_v / rho; /* derivative of conc_v convective terms */
          source_a_a = source_a * drho_c_a / rho; /* derivative of conc_a convective terms */
          source_a_l = source_a * drho_c_l / rho; /* derivative of conc_l convective terms */

          source_b_v = -drho_c_v * bf[var]->phi[j] / rho * source; /* derivative of density */
          source_b_a = -drho_c_a * bf[var]->phi[j] / rho * source; /* derivative of density */
          source_b_l = -drho_c_l * bf[var]->phi[j] / rho * source; /* derivative of density */

          source_c_v = d2rho_T_c_v * DT_Dt * bf[var]->phi[j] /
                       rho; /* cross-deriv: derivative of drho_T wrt x */
          source_c_a = d2rho_T_c_a * DT_Dt * bf[var]->phi[j] /
                       rho; /* cross-deriv: derivative of drho_T wrt x */

          dFVS_dC[species_v][j] = source_a_v + source_b_v + source_c_v;
          dFVS_dC[species_a][j] = source_a_a + source_b_a + source_c_a;
          dFVS_dC[species_l][j] = source_a_l + source_b_l;
        }
      }
    }
  } else if (mp->DensityModel == DENSITY_FOAM_TIME) {
    double rho_init, rho_final, aexp, time_delay, realtime;
    rho_init = mp->u_density[0];   /* Initial density */
    rho_final = mp->u_density[1];  /* final density */
    aexp = mp->u_density[2];       /* Arhennius constant for time exponent*/
    time_delay = mp->u_density[3]; /* time delay before foaming starts */

    if (time > time_delay) {
      realtime = time - time_delay;
      rho = rho_final + (rho_init - rho_final) * exp(-aexp * realtime);
      d_rho_dt = -aexp * (rho_init - rho_final) * exp(-aexp * realtime);
    } else {
      rho = rho_init;
      d_rho_dt = 0.0;
    }

    source = d_rho_dt / rho;
  } else if (mp->DensityModel == DENSITY_FOAM_TIME_TEMP) {
    double T;
    double rho_init = mp->u_density[0];   /* Initial density (units of gm cm-3 */
    double rho_final = mp->u_density[1];  /* Final density (units of gm cm-3 */
    double cexp = mp->u_density[2];       /* cexp in units of Kelvin sec */
    double coffset = mp->u_density[3];    /* offset in units of seconds */
    double time_delay = mp->u_density[4]; /* time delay in units of seconds */
    double drhoDT;
    double d_source_dT = 0.0;
    double d2_rho_dt_dT = 0.0;
    T = fv->T;
    if (time > time_delay) {
      double realtime = time - time_delay;
      double delRho = (rho_init - rho_final);
      double cdenom = cexp - coffset * T;
      double expT = exp(-realtime * T / cdenom);
      rho = rho_final + delRho * expT;
      drhoDT = delRho * expT * (-realtime * cexp / (cdenom * cdenom));
      d_rho_dt = -delRho * expT * T / cdenom;
      d2_rho_dt_dT = -delRho * expT *
                     (-realtime * cexp * T / (cdenom * cdenom * cdenom) + cexp / (cdenom * cdenom));
    } else {
      rho = rho_init;
      drhoDT = 0.0;
      d_rho_dt = 0.0;
    }
    source = d_rho_dt / rho;
    d_source_dT = d2_rho_dt_dT / rho - d_rho_dt / rho / rho * drhoDT;
    var = TEMPERATURE;
    if (pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        dFVS_dT[j] = d_source_dT * bf[var]->phi[j];
      }
    }

  } else if (mp->DensityModel == DENSITY_FOAM_PMDI_10) {
    if (pd->gv[MOMENT1]) {
      int wCO2Liq;
      int wCO2Gas;
      int wH2O;
      int w;
      DENSITY_DEPENDENCE_STRUCT d_rho_struct;
      DENSITY_DEPENDENCE_STRUCT *d_rho = &d_rho_struct;

      rho = density(d_rho, time);

      wCO2Liq = -1;
      wCO2Gas = -1;
      wH2O = -1;
      for (w = 0; w < pd->Num_Species; w++) {
        switch (mp->SpeciesSourceModel[w]) {
        case FOAM_PMDI_10_CO2_LIQ:
          wCO2Liq = w;
          break;
        case FOAM_PMDI_10_CO2_GAS:
          wCO2Gas = w;
          break;
        case FOAM_PMDI_10_H2O:
          wH2O = w;
          break;
        default:
          break;
        }
      }

      if (wCO2Liq == -1) {
        GOMA_EH(GOMA_ERROR, "Expected a Species Source of FOAM_PMDI_10_CO2_LIQ");
      } else if (wCO2Gas == -1) {
        GOMA_EH(GOMA_ERROR, "Expected a Species Source of FOAM_PMDI_10_CO2_GAS");
      } else if (wH2O == -1) {
        GOMA_EH(GOMA_ERROR, "Expected a Species Source of FOAM_PMDI_10_H2O");
      }

      double M_CO2 = mp->u_density[0];
      double rho_liq = mp->u_density[1];
      double ref_press = mp->u_density[2];
      double Rgas_const = mp->u_density[3];
      double rho_gas = 0;

      if (fv->T > 0) {
        rho_gas = (ref_press * M_CO2 / (Rgas_const * fv->T));
      }

      double nu = 0;

      double nu_dot = 0;

      double grad_nu[DIM];

      nu = fv->moment[1];

      nu_dot = fv_dot->moment[1];
      for (a = 0; a < dim; a++) {
        grad_nu[a] = fv->grad_moment[1][a];
      }

      double inv1 = 1 / (1 + nu);
      double inv2 = inv1 * inv1;

      // double volF = nu * inv1;
      // double d_volF_dC = (d_nu_dC) * inv2;
      // double d_volF_dT = (d_nu_dT) * inv2;

      double volF_dot = (nu_dot)*inv2;

      double grad_volF[DIM];
      for (a = 0; a < dim; a++) {
        grad_volF[a] = (grad_nu[a]) * inv2;
      }

      // double rho = rho_gas * volF + rho_liq * (1 - volF);
      double rho_dot = rho_gas * volF_dot - rho_liq * volF_dot;

      double grad_rho[DIM] = {0.0};
      double d_grad_rho_dT[DIM][MDE] = {{0.0}};

      for (a = 0; a < dim; a++) {
        grad_rho[a] = rho_gas * grad_volF[a] + rho_liq * (grad_volF[a]);
        if (af->Assemble_Jacobian) {
          var = TEMPERATURE;
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            if (fv->T > 0) {
              d_grad_rho_dT[a][j] += -(rho_gas / fv->T) * bf[var]->phi[j] * grad_volF[a];
            }
          }
        }
      }

      double inv_rho = 1 / rho;

      source = rho_dot;
      for (a = 0; a < dim; a++) {
        source += fv->v[a] * grad_rho[a];
      }
      source *= inv_rho;

      if (af->Assemble_Jacobian) {
        var = TEMPERATURE;
        if (pd->v[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            source_a = 0;
            source_b = 0;
            for (a = 0; a < dim; a++) {
              source_b += fv->v[a] * d_grad_rho_dT[a][j];
            }
            source_c = inv_rho * d_rho->T[j] * source;

            dFVS_dT[j] = inv_rho * (source_a + source_b) + source_c;
          }
        }

        var = VELOCITY1;
        if (pd->v[pg->imtrx][var]) {
          for (a = 0; a < dim; a++) {
            var = VELOCITY1 + a;
            for (j = 0; pd->v[pg->imtrx][var] && j < ei[pg->imtrx]->dof[var]; j++) {
              dFVS_dv[a][j] = inv_rho * (bf[var]->phi[j] * grad_rho[a]);
            }
          }
        }

        var = MASS_FRACTION;
        if (pd->v[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            source_a = 0;
            source_b = 0;
            for (a = 0; a < dim; a++) {
              source_b += fv->v[a] * 0;
            }
            source_c = inv_rho * d_rho->C[wCO2Liq][j] * source;

            dFVS_dC[wCO2Liq][j] = inv_rho * (source_a + source_b) + source_c;
          }
        }
        var = FILL;
        if (pd->v[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            source_c = inv_rho * d_rho->F[j] * source;

            dFVS_dF[j] = source_c;
          }
        }
      }
    } else {
      int wCO2;
      int wH2O;
      int w;
      DENSITY_DEPENDENCE_STRUCT d_rho_struct;
      DENSITY_DEPENDENCE_STRUCT *d_rho = &d_rho_struct;

      wCO2 = -1;
      wH2O = -1;
      for (w = 0; w < pd->Num_Species; w++) {
        switch (mp->SpeciesSourceModel[w]) {
        case FOAM_PMDI_10_CO2:
          wCO2 = w;
          break;
        case FOAM_PMDI_10_H2O:
          wH2O = w;
          break;
        default:
          break;
        }
      }

      if (wCO2 == -1) {
        GOMA_EH(GOMA_ERROR, "Expected a Species Source of FOAM_PMDI_10_CO2");
      } else if (wH2O == -1) {
        GOMA_EH(GOMA_ERROR, "Expected a Species Source of FOAM_PMDI_10_H2O");
      }

      double M_CO2 = mp->u_density[0];
      double rho_liq = mp->u_density[1];
      double ref_press = mp->u_density[2];
      double Rgas_const = mp->u_density[3];
      double rho_gas = 0;

      if (fv->T > 0) {
        rho_gas = (ref_press * M_CO2 / (Rgas_const * fv->T));
      }

      double nu = 0;
      double d_nu_dC = 0;
      double d_nu_dT = 0;

      double nu_dot = 0;
      double d_nu_dot_dC = 0;
      double d_nu_dot_dT = 0;

      double grad_nu[DIM];
      double d_grad_nu_dC[DIM][MDE];
      double d_grad_nu_dT[DIM];

      if (fv->T > 0) {
        nu = M_CO2 * fv->c[wCO2] / rho_gas;
        d_nu_dC = M_CO2 / rho_gas;
        d_nu_dT = M_CO2 * fv->c[wCO2] * Rgas_const / (ref_press * M_CO2);

        nu_dot = M_CO2 * fv_dot->c[wCO2] / rho_gas;
        d_nu_dot_dC = M_CO2 * (1 + 2 * tt) / (dt * rho_gas);
        d_nu_dot_dT = M_CO2 * fv_dot->c[wCO2] * Rgas_const / (ref_press * M_CO2);

        for (a = 0; a < dim; a++) {
          grad_nu[a] = M_CO2 * fv->grad_c[wCO2][a] / rho_gas;
          if (af->Assemble_Jacobian) {
            d_grad_nu_dT[a] = M_CO2 * fv->grad_c[wCO2][a] * Rgas_const / (ref_press * M_CO2);
            for (j = 0; j < ei[upd->matrix_index[MASS_FRACTION]]->dof[MASS_FRACTION]; j++) {
              d_grad_nu_dC[a][j] = M_CO2 * bf[MASS_FRACTION]->grad_phi[j][a] / rho_gas;
            }
          }
        }
      } else {
        nu = 0;
        d_nu_dC = 0;
        d_nu_dT = 0;

        nu_dot = 0;
        d_nu_dot_dC = 0;
        d_nu_dot_dT = 0;

        for (a = 0; a < dim; a++) {
          grad_nu[a] = 0;
          if (af->Assemble_Jacobian) {
            d_grad_nu_dT[a] = 0;
            for (j = 0; j < ei[upd->matrix_index[MASS_FRACTION]]->dof[MASS_FRACTION]; j++) {
              d_grad_nu_dC[a][j] = 0;
            }
          }
        }
      }

      double inv1 = 1 / (1 + nu);
      double inv2 = inv1 * inv1;
      double inv3 = inv1 * inv2;

      // double volF = nu * inv1;
      // double d_volF_dC = (d_nu_dC) * inv2;
      // double d_volF_dT = (d_nu_dT) * inv2;

      double volF_dot = (nu_dot)*inv2;
      double d_volF_dot_dC = (d_nu_dot_dC * (1 + nu) - nu_dot * 2 * d_nu_dC) * inv3;
      double d_volF_dot_dT = (d_nu_dot_dT * (1 + nu) - nu_dot * 2 * d_nu_dT) * inv3;

      double grad_volF[DIM];
      double d_grad_volF_dC[DIM][MDE];
      double d_grad_volF_dT[DIM];
      for (a = 0; a < dim; a++) {
        grad_volF[a] = (grad_nu[a]) * inv2;
        if (af->Assemble_Jacobian) {
          d_grad_volF_dT[a] = (d_grad_nu_dT[a] * (1 + nu) - 2 * grad_nu[a] * d_nu_dT) * inv3;
          for (j = 0; j < ei[upd->matrix_index[MASS_FRACTION]]->dof[MASS_FRACTION]; j++) {
            d_grad_volF_dC[a][j] = (d_grad_nu_dC[a][j] * (1 + nu) -
                                    2 * grad_nu[a] * d_nu_dC * bf[MASS_FRACTION]->phi[j]) *
                                   inv3;
          }
        }
      }

      // double rho = rho_gas * volF + rho_liq * (1 - volF);
      double rho_dot = rho_gas * volF_dot - rho_liq * volF_dot;
      double d_rho_dot_dT = rho_gas * d_volF_dot_dT - rho_liq * d_volF_dot_dT;
      if (fv->T > 0) {
        d_rho_dot_dT =
            -(rho_gas / fv->T) * volF_dot + rho_gas * d_volF_dot_dT - rho_liq * d_volF_dot_dT;
      }
      double d_rho_dot_dC = rho_gas * d_volF_dot_dC - rho_liq * d_volF_dot_dC;

      double grad_rho[DIM];
      double d_grad_rho_dT[DIM][MDE];
      double d_grad_rho_dC[DIM][MDE];

      for (a = 0; a < dim; a++) {
        grad_rho[a] = rho_gas * grad_volF[a] + rho_liq * (grad_volF[a]);
        if (af->Assemble_Jacobian) {
          var = TEMPERATURE;
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            d_grad_rho_dT[a][j] = rho_gas * d_grad_volF_dT[a] * bf[var]->phi[j] +
                                  rho_liq * (d_grad_volF_dT[a] * bf[var]->phi[j]);
            if (fv->T > 0) {
              d_grad_rho_dT[a][j] += -(rho_gas / fv->T) * bf[var]->phi[j] * grad_volF[a];
            }
          }

          var = MASS_FRACTION;
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            d_grad_rho_dC[a][j] = rho_gas * d_grad_volF_dC[a][j] + rho_liq * (d_grad_volF_dC[a][j]);
          }
        }
      }

      double volF = mp->volumeFractionGas;
      double rho = rho_gas * volF + rho_liq * (1 - volF);
      double inv_rho = 1 / rho;
      var = MASS_FRACTION;
      if (volF > 0. && d_rho != NULL) {
        if (pd->v[pg->imtrx][var]) {
          for (w = 0; w < pd->Num_Species; w++) {
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              d_rho->C[w][j] = (rho_gas * mp->d_volumeFractionGas[MAX_VARIABLE_TYPES + w] -
                                rho_liq * mp->d_volumeFractionGas[MAX_VARIABLE_TYPES + w]) *
                               bf[var]->phi[j];
            }
          }
        }
      }

      var = TEMPERATURE;
      if (d_rho != NULL) {
        if (pd->v[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            if (fv->T > 0) {
              d_rho->T[j] = (-rho_gas / fv->T * volF + rho_gas * mp->d_volumeFractionGas[var] -
                             rho_liq * mp->d_volumeFractionGas[var]) *
                            bf[var]->phi[j];
            } else {
              d_rho->T[j] = (rho_gas * mp->d_volumeFractionGas[var] -
                             rho_liq * mp->d_volumeFractionGas[var]) *
                            bf[var]->phi[j];
            }
          }
        }
      }
      source = rho_dot;
      for (a = 0; a < dim; a++) {
        source += fv->v[a] * grad_rho[a];
      }
      source *= inv_rho;

      if (af->Assemble_Jacobian) {
        var = TEMPERATURE;
        if (pd->v[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            source_a = d_rho_dot_dT * bf[var]->phi[j];
            source_b = 0;
            for (a = 0; a < dim; a++) {
              source_b += fv->v[a] * d_grad_rho_dT[a][j];
            }
            source_c = inv_rho * d_rho->T[j] * source;

            dFVS_dT[j] = inv_rho * (source_a + source_b) + source_c;
          }
        }

        var = VELOCITY1;
        if (pd->v[pg->imtrx][var]) {
          for (a = 0; a < dim; a++) {
            var = VELOCITY1 + a;
            for (j = 0; pd->v[pg->imtrx][var] && j < ei[pg->imtrx]->dof[var]; j++) {
              dFVS_dv[a][j] = inv_rho * (bf[var]->phi[j] * grad_rho[a]);
            }
          }
        }

        var = MASS_FRACTION;
        if (pd->v[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            source_a = d_rho_dot_dC * bf[var]->phi[j];
            source_b = 0;
            for (a = 0; a < dim; a++) {
              source_b += fv->v[a] * d_grad_rho_dC[a][j];
            }
            source_c = inv_rho * d_rho->C[wCO2][j] * source;

            dFVS_dC[wCO2][j] = inv_rho * (source_a + source_b) + source_c;
          }
        }

        var = FILL;
        if (pd->v[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            source_c = inv_rho * d_rho->F[j] * source;

            dFVS_dF[j] = source_c;
          }
        }
      }
    }
  } else if (mp->DensityModel == DENSITY_MOMENT_BASED) {
    if (pd->gv[MOMENT1]) {
      DENSITY_DEPENDENCE_STRUCT d_rho_struct;
      DENSITY_DEPENDENCE_STRUCT *d_rho = &d_rho_struct;

      rho = density(d_rho, time);

      double rho_gas = mp->u_density[0];
      double rho_liq = mp->u_density[1];

      double nu = 0;

      double nu_dot = 0;

      double grad_nu[DIM];

      nu = fv->moment[1];

      nu_dot = fv_dot->moment[1];
      for (a = 0; a < dim; a++) {
        grad_nu[a] = fv->grad_moment[1][a];
      }

      double inv1 = 1 / (1 + nu);
      double inv2 = inv1 * inv1;

      // double volF = nu * inv1;
      // double d_volF_dC = (d_nu_dC) * inv2;
      // double d_volF_dT = (d_nu_dT) * inv2;

      double volF_dot = (nu_dot)*inv2;

      double grad_volF[DIM];
      for (a = 0; a < dim; a++) {
        grad_volF[a] = (grad_nu[a]) * inv2;
      }

      // double rho = rho_gas * volF + rho_liq * (1 - volF);
      double rho_dot = rho_gas * volF_dot - rho_liq * volF_dot;

      double grad_rho[DIM] = {0.0};
      double d_grad_rho_dT[DIM][MDE] = {{0.0}};

      for (a = 0; a < dim; a++) {
        grad_rho[a] = rho_gas * grad_volF[a] + rho_liq * (grad_volF[a]);
        if (af->Assemble_Jacobian) {
          var = TEMPERATURE;
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            if (fv->T > 0) {
              d_grad_rho_dT[a][j] += -(rho_gas / fv->T) * bf[var]->phi[j] * grad_volF[a];
            }
          }
        }
      }

      double inv_rho = 1 / rho;

      source = rho_dot;
      for (a = 0; a < dim; a++) {
        source += fv->v[a] * grad_rho[a];
      }
      source *= inv_rho;

      if (af->Assemble_Jacobian) {
        var = TEMPERATURE;
        if (pd->v[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            source_a = 0;
            source_b = 0;
            for (a = 0; a < dim; a++) {
              source_b += fv->v[a] * d_grad_rho_dT[a][j];
            }
            source_c = inv_rho * d_rho->T[j] * source;

            dFVS_dT[j] = inv_rho * (source_a + source_b) + source_c;
          }
        }

        var = VELOCITY1;
        if (pd->v[pg->imtrx][var]) {
          for (a = 0; a < dim; a++) {
            var = VELOCITY1 + a;
            for (j = 0; pd->v[pg->imtrx][var] && j < ei[pg->imtrx]->dof[var]; j++) {
              dFVS_dv[a][j] = inv_rho * (bf[var]->phi[j] * grad_rho[a]);
            }
          }
        }

        var = MASS_FRACTION;
        if (pd->v[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            dFVS_dC[0][j] = 0;
          }
        }
        var = FILL;
        if (pd->v[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            source_c = inv_rho * d_rho->F[j] * source;

            dFVS_dF[j] = source_c;
          }
        }
      }
    }
  }

  if (ls != NULL && mp->mp2nd != NULL && (mp->DensityModel != LEVEL_SET) &&
      (mp->mp2nd->DensityModel == CONSTANT)) {
    double factor;

    source = ls_modulate_property(source, 0.0, ls->Length_Scale, (double)mp->mp2nd->densitymask[0],
                                  (double)mp->mp2nd->densitymask[1], dFVS_dF, &factor);

    var = TEMPERATURE;

    for (j = 0; pd->v[pg->imtrx][var] && j < ei[pg->imtrx]->dof[var]; j++) {
      dFVS_dT[j] *= factor;
    }

    var = VELOCITY1;
    if (pd->v[pg->imtrx][var]) {
      for (a = 0; a < dim; a++) {
        var = VELOCITY1 + a;
        for (j = 0; pd->v[pg->imtrx][var] && j < ei[pg->imtrx]->dof[var]; j++) {
          dFVS_dv[a][j] *= factor;
        }
      }
    }

    var = MASS_FRACTION;
    if (pd->v[pg->imtrx][var]) {
      for (w = 0; w < pd->Num_Species_Eqn; w++) {

        /* source = ls_modulate_property( source,
           0.0,
           ls->Length_Scale,
           (double) mp->mp2nd->densitymask[0],
           (double) mp->mp2nd->densitymask[1],
           dFVS_dC[w],
           &factor );*/
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {

          dFVS_dC[w][j] *= factor;
        }
      }
    }
  }

  return (source);
} /*******************************************************************************/
double REFVolumeSource(double time,
                       double dt,
                       double tt,
                       double dFVS_dv[DIM][MDE],
                       double dFVS_dT[MDE],
                       double dFVS_dx[DIM][MDE],
                       double dFVS_dC[MAX_CONC][MDE]) {
  int a, b, j, w, w1, dim = pd->Num_Dim, var;
  double x_dot[DIM] = {0., 0., 0.};
  double source = 0.;
  double phi_j;
  double rho = 0.;
  double sv[MAX_CONC], sv_p = 0.;
  double dtmp[MAX_CONC], d2tmp[MAX_CONC][MAX_CONC];

  memset(dFVS_dv, 0, sizeof(double) * DIM * MDE);
  memset(dFVS_dT, 0, sizeof(double) * MDE);
  memset(dFVS_dx, 0, sizeof(double) * DIM * MDE);
  memset(dtmp, 0, sizeof(double) * MAX_CONC);
  memset(sv, 0, sizeof(double) * MAX_CONC);
  memset(d2tmp, 0, sizeof(double) * MAX_CONC * MAX_CONC);

  if (mp->SpeciesSourceModel[0] == FOAM) {
    foam_species_source(mp->u_species_source[0]);

  } else {
    GOMA_EH(GOMA_ERROR, "Must specify FOAM species source in the material's file");
  }

  if (pd->TimeIntegration != STEADY) {
    if (pd->v[pg->imtrx][MESH_DISPLACEMENT1]) {
      for (a = 0; a < dim; a++) {
        x_dot[a] = fv_dot->x[a];
      }
    }
  }
  rho = density(NULL, time);

  /* must include ref reaction source */

  sv_p = mp->specific_volume[pd->Num_Species_Eqn];

  for (w = 0; w < pd->Num_Species_Eqn; w++) {
    sv[w] = mp->specific_volume[w];
    dtmp[w] = -rho * (sv[w] - sv_p);
    source += dtmp[w] * mp->species_source[w];

    /* must account for spatial variation of density */

    for (a = 0; a < dim; a++) {
      source += (fv->v[a] - x_dot[a]) * dtmp[w] * fv->grad_c[w][a];
    }
  }

  for (w = 0; w < pd->Num_Species_Eqn; w++) {
    for (w1 = 0; w1 < pd->Num_Species_Eqn; w1++) {
      d2tmp[w][w1] = -rho * (sv[w] - sv_p) * dtmp[w1];
    }
  }

  if (af->Assemble_Jacobian) {
    var = VELOCITY1;

    for (j = 0; pd->v[pg->imtrx][var] && j < ei[pg->imtrx]->dof[var]; j++) {
      phi_j = bf[var]->phi[j];

      for (a = 0; a < dim; a++) {
        dFVS_dv[a][j] = 0.;
        for (w = 0; w < pd->Num_Species_Eqn; w++) {
          dFVS_dv[a][j] += phi_j * dtmp[w] * fv->grad_c[w][a];
        }
      }
    }

    var = TEMPERATURE;

    for (j = 0; pd->v[pg->imtrx][var] && j < ei[pg->imtrx]->dof[var]; j++) {
      dFVS_dT[j] = 0.;
      phi_j = bf[var]->phi[j];

      for (w = 0; w < pd->Num_Species_Eqn; w++) {
        dFVS_dT[j] -= (sv[w] - sv_p) * (rho * mp->d_species_source[MAX_VARIABLE_TYPES + w] +
                                        mp->species_source[w] * 3. * rho * rho);
      }
      dFVS_dT[j] *= phi_j;
    }

    for (b = 0; b < dim; b++) {
      var = MESH_DISPLACEMENT1 + b;

      for (j = 0; pd->v[pg->imtrx][var] && j < ei[pg->imtrx]->dof[var]; j++) {
        dFVS_dx[b][j] = 0.0;
        phi_j = bf[var]->phi[j];

        for (w = 0; w < pd->Num_Species_Eqn; w++) {
          for (a = 0; a < dim; a++) {
            dFVS_dx[b][j] += (fv->v[a] - x_dot[a]) * dtmp[w] * fv->d_grad_c_dmesh[a][w][b][j];
          }
          if (TimeIntegration != STEADY)
            dFVS_dx[b][j] += -(1.0 + 2 * tt) * phi_j * dtmp[w] * fv->grad_c[w][b] / dt;
        }
      }
    }

    var = MASS_FRACTION;

    for (w = 0; w < pd->Num_Species_Eqn; w++) {
      for (j = 0; pd->v[pg->imtrx][var] && j < ei[pg->imtrx]->dof[var]; j++) {
        dFVS_dC[w][j] = 0.0;
        phi_j = bf[var]->phi[j];

        for (w1 = 0; w1 < pd->Num_Species_Eqn; w1++) {
          for (a = 0; a < dim; a++) {
            dFVS_dC[w][j] += (fv->v[a] - x_dot[a]) * (d2tmp[w1][w] * fv->grad_c[w1][a]);

            if (w1 == w)
              dFVS_dC[w][j] += (fv->v[a] - x_dot[a]) * dtmp[w] * bf[var]->grad_phi[j][a];
          }

          /* dFVS_dC[w][j] -= (sv[w1]-sv_p)*
            (rho*mp->Jac_Species_Source[w1*pd->Num_Species+w] +
            +dtmp[w]*mp->species_source[w1]); */

          dFVS_dC[w][j] += (dtmp[w1] * mp->Jac_Species_Source[w1 * pd->Num_Species + w] +
                            +d2tmp[w1][w] * mp->species_source[w1]);
        }
      }
    }
  }
  return (source);
}
int assemble_projection_time_stabilization(Exo_DB *exo, double time, double tt, double dt) {
  double phi_i, phi_j;
  int var, pvar, ip, i, j, a, b, ipass, num_passes;
  double xi[DIM];
  /* Variables for vicosity and derivative */
  double gamma[DIM][DIM]; /* shrearrate tensor based on velocity */
  double mu;
  VISCOSITY_DEPENDENCE_STRUCT d_mu_struct; /* viscosity dependence */
  VISCOSITY_DEPENDENCE_STRUCT *d_mu = &d_mu_struct;

  int ip_total = elem_info(NQUAD, ei[pg->imtrx]->ielem_type);
  int eqn = R_PRESSURE;
  int peqn = upd->ep[pg->imtrx][eqn];
  double wt, det_J, h3, d_vol;

  if (ls == NULL) {
    num_passes = 1;
  } else {
    if (ls->elem_overlap_state)
      num_passes = 2;
    else
      num_passes = 1;
  }

  for (ipass = 0; ipass < num_passes; ipass++) {
    if (ls != NULL) {
      if (num_passes == 2)
        ls->Elem_Sign = -1 + 2 * ipass;
      else
        ls->Elem_Sign = 0;
    }

    double P_dot_avg = 0.;
    double phi_avg[MDE];
    double vol = 0.;
    double ls_scale_factor = 1.0; /* this is to turn off PSPP near the interface */

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      phi_avg[i] = 0.;
    }

    for (ip = 0; ip < ip_total; ip++) {
      find_stu(ip, ei[pg->imtrx]->ielem_type, &xi[0], &xi[1], &xi[2]); /* find quadrature point */
      wt = Gq_weight(ip, ei[pg->imtrx]->ielem_type); /* find quadrature weights for */

      setup_shop_at_point(ei[pg->imtrx]->ielem, xi, exo);

      computeCommonMaterialProps_gp(time);

      fv->wt = wt;
      h3 = fv->h3;
      det_J = bf[eqn]->detJ;
      d_vol = wt * det_J * h3;

      P_dot_avg += fv_dot->P * d_vol;

      for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
        phi_i = bf[eqn]->phi[i];
        phi_avg[i] += phi_i * d_vol;
      }
      vol += d_vol;
    }

    double inv_vol = 1. / vol;
    P_dot_avg *= inv_vol;
    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      phi_avg[i] *= inv_vol;
    }

    for (ip = 0; ip < ip_total; ip++) {
      find_stu(ip, ei[pg->imtrx]->ielem_type, &xi[0], &xi[1], &xi[2]); /* find quadrature point */
      wt = Gq_weight(ip, ei[pg->imtrx]->ielem_type); /* find quadrature weights for */

      setup_shop_at_point(ei[pg->imtrx]->ielem, xi, exo);
      fv->wt = wt;
      h3 = fv->h3;
      det_J = bf[eqn]->detJ;
      d_vol = wt * det_J * h3;

      computeCommonMaterialProps_gp(time);

      if (ls != NULL && ls->PSPP_filter) {
        load_lsi(ls->Length_Scale);

        ls_scale_factor =
            ls->on_sharp_surf
                ? 1.0
                : 1.0 -
                      lsi->delta / lsi->delta_max; /*Punting in the case of subelementintegration */
        if (af->Assemble_Jacobian)
          load_lsi_derivs();
      }

      /* load up shearrate tensor based on velocity */
      for (a = 0; a < VIM; a++) {
        for (b = 0; b < VIM; b++) {
          gamma[a][b] = fv->grad_v[a][b] + fv->grad_v[b][a];
        }
      }
      mu = viscosity(gn, gamma, d_mu);

      if (af->Assemble_Residual) {
        for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
          phi_i = bf[eqn]->phi[i];
          lec->R[LEC_R_INDEX(peqn, i)] += ls_scale_factor * PS_scaling * dt *
                                          (fv_dot->P - P_dot_avg) * (phi_i - phi_avg[i]) / mu *
                                          d_vol;
        }
      }
      if (af->Assemble_Jacobian) {

        for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
          phi_i = bf[eqn]->phi[i];

          var = PRESSURE;
          pvar = upd->vp[pg->imtrx][var];
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];
            lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += ls_scale_factor * dt * PS_scaling *
                                                     (((1 + 2. * tt) / dt) * (phi_j - phi_avg[j])) *
                                                     (phi_i - phi_avg[i]) / mu * d_vol;
          }

          if (pd->v[pg->imtrx][var = LS]) {
            pvar = upd->vp[pg->imtrx][var];
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              phi_j = bf[var]->phi[j];

              if (ls != NULL && ls->PSPP_filter) {

                lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += -(lsi->d_delta_dF[j] / lsi->delta_max) *
                                                         dt * PS_scaling * (fv_dot->P - P_dot_avg) *
                                                         (phi_i - phi_avg[i]) / mu * d_vol;
              }
            }
          }
        }
      }
    }
  }
  return (1);
}
int assemble_projection_stabilization(Exo_DB *exo, double time)
/* This routine applies the Dohrmann-Bochev Polynomial Projection Pressure Stabilization
 * to the continuity equation to help aid in iterative solutions the the Navier-Stokes
 * equations This routine also includes an expedient to ignore such terms near a level-set
 * interface. Not sure how those work, though (PRS 2/27/2007).    this routine is basic, and
 * doesn't necessarily apply for graded meshes.  See assemble_PPPS_generalized routine for
 * that.
 */
{
  double phi_i, phi_j;
  int var, pvar, ip, i, j, a, b, ipass, num_passes;
  double xi[DIM];
  /* Variables for vicosity and derivative */
  double gamma[DIM][DIM]; /* shrearrate tensor based on velocity */
  double mu;
  VISCOSITY_DEPENDENCE_STRUCT d_mu_struct; /* viscosity dependence */
  VISCOSITY_DEPENDENCE_STRUCT *d_mu = &d_mu_struct;

  int ip_total = elem_info(NQUAD, ei[pg->imtrx]->ielem_type);
  int eqn = R_PRESSURE;
  int peqn = upd->ep[pg->imtrx][eqn];
  double wt, det_J, h3, d_vol;

  if (ls == NULL) {
    num_passes = 1;
  } else {
    if (ls->elem_overlap_state)
      num_passes = 2;
    else
      num_passes = 1;
  }

  for (ipass = 0; ipass < num_passes; ipass++) {
    if (ls != NULL) {
      if (num_passes == 2)
        ls->Elem_Sign = -1 + 2 * ipass;
      else
        ls->Elem_Sign = 0;
    }

    double P_avg = 0.;
    double phi_avg[MDE];
    double vol = 0.;
    double ls_scale_factor = 1.0; /* this is to turn off PSPP near the interface */

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      phi_avg[i] = 0.;
    }

    for (ip = 0; ip < ip_total; ip++) {
      find_stu(ip, ei[pg->imtrx]->ielem_type, &xi[0], &xi[1], &xi[2]); /* find quadrature point */
      wt = Gq_weight(ip, ei[pg->imtrx]->ielem_type); /* find quadrature weights for */

      setup_shop_at_point(ei[pg->imtrx]->ielem, xi, exo);

      computeCommonMaterialProps_gp(time);

      fv->wt = wt;
      h3 = fv->h3;
      det_J = bf[eqn]->detJ;
      d_vol = wt * det_J * h3;

      P_avg += fv->P * d_vol;

      for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
        phi_i = bf[eqn]->phi[i];
        phi_avg[i] += phi_i * d_vol;
      }
      vol += d_vol;
    }

    double inv_vol = 1. / vol;
    P_avg *= inv_vol;
    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      phi_avg[i] *= inv_vol;
    }

    for (ip = 0; ip < ip_total; ip++) {
      find_stu(ip, ei[pg->imtrx]->ielem_type, &xi[0], &xi[1], &xi[2]); /* find quadrature point */
      wt = Gq_weight(ip, ei[pg->imtrx]->ielem_type); /* find quadrature weights for */

      setup_shop_at_point(ei[pg->imtrx]->ielem, xi, exo);
      fv->wt = wt;
      h3 = fv->h3;
      det_J = bf[eqn]->detJ;
      d_vol = wt * det_J * h3;
      computeCommonMaterialProps_gp(time);

      if (ls != NULL && ls->PSPP_filter) {
        load_lsi(ls->Length_Scale);

        ls_scale_factor =
            ls->on_sharp_surf
                ? 1.0
                : 1.0 -
                      lsi->delta / lsi->delta_max; /*Punting in the case of subelementintegration */
        if (af->Assemble_Jacobian)
          load_lsi_derivs();
      }

      /* load up shearrate tensor based on velocity */
      for (a = 0; a < VIM; a++) {
        for (b = 0; b < VIM; b++) {
          gamma[a][b] = fv->grad_v[a][b] + fv->grad_v[b][a];
        }
      }
      mu = viscosity(gn, gamma, d_mu);

      if (af->Assemble_Residual) {
        for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
          phi_i = bf[eqn]->phi[i];
          lec->R[LEC_R_INDEX(peqn, i)] +=
              ls_scale_factor * PS_scaling * (fv->P - P_avg) * (phi_i - phi_avg[i]) / mu * d_vol;
        }
      }
      if (af->Assemble_Jacobian) {

        for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
          phi_i = bf[eqn]->phi[i];

          var = PRESSURE;
          pvar = upd->vp[pg->imtrx][var];
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];
            lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += ls_scale_factor * PS_scaling *
                                                     (phi_j - phi_avg[j]) * (phi_i - phi_avg[i]) /
                                                     mu * d_vol;
          }

          if (pd->v[pg->imtrx][var = LS]) {
            pvar = upd->vp[pg->imtrx][var];
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              phi_j = bf[var]->phi[j];

              if (ls != NULL && ls->PSPP_filter) {

                lec->J[LEC_J_INDEX(peqn, pvar, i, j)] -= ls_scale_factor * PS_scaling *
                                                         (fv->P - P_avg) * (phi_i - phi_avg[i]) *
                                                         d_mu->F[j] / (mu * mu) * d_vol;

                lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += -(lsi->d_delta_dF[j] / lsi->delta_max) *
                                                         PS_scaling * (fv->P - P_avg) *
                                                         (phi_i - phi_avg[i]) / mu * d_vol;
              }
            }
          }
        }
      }
    }
  }
  return (1);
}
int assemble_PPPS_generalized(Exo_DB *exo)

/* This routine applies the Dohrmann-Bochev Polynomial Projection Pressure Stabilization
 * to the continuity equation to help aid in iterative solutions the the Navier-Stokes
 * equations Unlike its couterpart assemble_projection_stabilization, this routine adds the
 * expedients necessary to stabilized elements that have significant aspect ratios.   Viz.
 * they most solve an eigenproblem at the element level to do so.
 */
{
  int var, pvar, ip, i, j, a, b, m, n, k;
  double xi[DIM];
  /* Variables for vicosity and derivative */
  double gamma[DIM][DIM]; /* shrearrate tensor based on velocity */
  double mu, visc_e = 0.0, d_visc_e_dv[DIM][MDE], d_visc_e_dx[DIM][MDE];
  VISCOSITY_DEPENDENCE_STRUCT d_mu_struct; /* viscosity dependence */
  VISCOSITY_DEPENDENCE_STRUCT *d_mu = &d_mu_struct;

  int ip_total = elem_info(NQUAD, ei[pg->imtrx]->ielem_type);
  int eqn = R_PRESSURE;
  int peqn = upd->ep[pg->imtrx][eqn];
  int INFO, ipass, npass;
  double wt, det_J, h3, d_vol, sum;

  double vol = 0.;
  double xi_dV = 0.;
  double wm_dV = 0.;
  double EeT_We_Ee = 0.;

  /* Declare and initialize some element-level matrices needed for this stabilization (cf.
   * Dorhmann/Bochev notes/paper) IJNMF 46:183-201. 2004.
   */

  double x_bar[DIM], Vij[DIM][DIM], em[MDE], eim[DIM][MDE], eim_tilde[DIM][MDE];
  double Vij_tilde[DIM][DIM], Ee_tilde[DIM + 1][MDE], Ee_hat[DIM + 1][MDE];
  double Me[MDE][MDE], Ce[MDE][MDE], W[DIM], eval[DIM], evec[DIM][DIM];
  double De_hat[DIM + 1][DIM + 1], We[DIM + 1][DIM + 1], S[DIM + 1][DIM + 1], WeEe[DIM + 1][MDE];
  double WORK[100];
  double *A;
  int LWORK = 100;
  double scaling_over_visc_e, my_multiplier = 1.0;
  double R_new[MDE], R_old[MDE]; /*old and new resid vect pieces for the pressure equation */

  asdv(&A, ei[pg->imtrx]->ielem_dim * ei[pg->imtrx]->ielem_dim);

  npass = 1;

  if (pd->v[pg->imtrx][MESH_DISPLACEMENT1] && af->Assemble_Jacobian) {
    npass = 2;
    GOMA_EH(GOMA_ERROR, "need to work out numerical jacobian for PSPP and moving mesh");
  }

  for (ipass = 0; ipass < npass; ipass++) {

    memset(Me, 0, sizeof(double) * MDE * MDE);
    memset(eim, 0, sizeof(double) * DIM * MDE);
    memset(em, 0, sizeof(double) * MDE);
    memset(x_bar, 0, sizeof(double) * DIM);
    memset(Vij, 0, sizeof(double) * DIM * DIM);
    memset(d_visc_e_dv, 0, sizeof(double) * DIM * MDE);
    memset(d_visc_e_dx, 0, sizeof(double) * DIM * MDE);
    memset(R_new, 0, sizeof(double) * MDE);
    memset(R_old, 0, sizeof(double) * MDE);
    // phi_avg[MDE];
    vol = 0.;
    xi_dV = 0.;
    wm_dV = 0.;
    EeT_We_Ee = 0.;
    visc_e = 0.;

    for (ip = 0; ip < ip_total; ip++) {

      find_stu(ip, ei[pg->imtrx]->ielem_type, &xi[0], &xi[1], &xi[2]); /* find quadrature point */

      wt = Gq_weight(ip, ei[pg->imtrx]->ielem_type); /* find quadrature weights for */

      setup_shop_at_point(ei[pg->imtrx]->ielem, xi, exo);

      fv->wt = wt;
      h3 = fv->h3;
      det_J = bf[eqn]->detJ;
      d_vol = wt * det_J * h3;
      vol += d_vol;

      /* load up shearrate tensor based on velocity */
      for (a = 0; a < VIM; a++) {
        for (b = 0; b < VIM; b++) {
          gamma[a][b] = fv->grad_v[a][b] + fv->grad_v[b][a];
        }
      }
      mu = viscosity(gn, gamma, d_mu);
      if (gn->ConstitutiveEquation != NEWTONIAN)
        GOMA_EH(GOMA_ERROR, "PSPP not available for non_Newtonian yet.   Sorry.");

      visc_e += mu * d_vol;

      for (i = 0; i < ei[pg->imtrx]->ielem_dim; i++) {
        xi_dV = fv->x[i] * d_vol;
        x_bar[i] += xi_dV;
        for (j = 0; j < ei[pg->imtrx]->ielem_dim; j++) {
          Vij[i][j] += fv->x[j] * xi_dV;
        }
        for (m = 0; m < ei[pg->imtrx]->dof[eqn]; m++) {
          eim[i][m] += bf[eqn]->phi[m] * xi_dV;
        }
      }

      for (m = 0; m < ei[pg->imtrx]->dof[eqn]; m++) {
        wm_dV = bf[eqn]->phi[m] * d_vol;
        em[m] += wm_dV;
        for (n = 0; n < ei[pg->imtrx]->dof[eqn]; n++) {
          Me[m][n] += bf[eqn]->phi[n] * wm_dV;
        }
      }

      /* Pick up some Jacobian pieces for later down in the routine */

      if (af->Assemble_Jacobian && ipass != 1) {
        var = VELOCITY1;
        if (pd->v[pg->imtrx][var]) {
          for (a = 0; a < ei[pg->imtrx]->ielem_dim; a++) {
            var = VELOCITY1 + a;

            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              d_visc_e_dv[a][j] += d_vol * (*d_mu->v[a]);
            }
          }
        }

        var = MESH_DISPLACEMENT1;
        if (pd->v[pg->imtrx][var]) {
          /* We are doing all of mesh derivatives numerically */
        }
      }
    } /* End of integration loop */

    visc_e /= vol;
    if (af->Assemble_Jacobian) {
      var = VELOCITY1;
      if (pd->v[pg->imtrx][var]) {
        for (a = 0; a < ei[pg->imtrx]->ielem_dim; a++) {
          var = VELOCITY1 + a;

          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            d_visc_e_dv[a][j] /= vol;
          }
        }
      }
    }

    for (i = 0; i < ei[pg->imtrx]->ielem_dim; i++) {
      x_bar[i] /= vol;
    }

    for (i = 0; i < ei[pg->imtrx]->ielem_dim; i++) {
      for (m = 0; m < ei[pg->imtrx]->dof[eqn]; m++) {
        eim_tilde[i][m] = eim[i][m] - x_bar[i] * em[m];
      }
      for (j = 0; j < ei[pg->imtrx]->ielem_dim; j++) {
        Vij_tilde[i][j] = Vij[i][j] - vol * x_bar[i] * x_bar[j];
      }
    }
    for (m = 0; m < ei[pg->imtrx]->dof[eqn]; m++) {
      Ee_tilde[0][m] = em[m]; /* First row of Ee_tilde is special */
      for (i = 0; i < ei[pg->imtrx]->ielem_dim; i++) {
        Ee_tilde[i + 1][m] = eim_tilde[i][m];
      }
    }

    /* Next we need to solve the eigen problem */

    for (i = 0; i < ei[pg->imtrx]->ielem_dim; i++) {
      W[i] = 0;
      for (j = 0; j < ei[pg->imtrx]->ielem_dim; j++) {
        A[j * ei[pg->imtrx]->ielem_dim + i] = Vij_tilde[i][j];
      }
    }

    dsyev_("V", "U", &(ei[pg->imtrx]->ielem_dim), A, &(ei[pg->imtrx]->ielem_dim), W, WORK, &LWORK,
           &INFO, 1, 1);

    if (INFO > 0)
      GOMA_EH(GOMA_ERROR, "dsyev falied to converge");
    if (INFO < 0)
      GOMA_EH(GOMA_ERROR, "an argument of dsyev had an illegal value");

    for (j = 0; j < ei[pg->imtrx]->ielem_dim; j++) {
      eval[j] = W[j];
      for (i = 0; i < ei[pg->imtrx]->ielem_dim; i++) {
        evec[j][i] = A[j * ei[pg->imtrx]->ielem_dim + i];
      }
    }

    /*user De_hat as scratch for calculating Ee_hat */

    memset(De_hat, 0, sizeof(double) * (DIM + 1) * (DIM + 1));
    De_hat[0][0] = 1.;
    for (i = 0; i < ei[pg->imtrx]->ielem_dim; i++) {
      for (j = 0; j < ei[pg->imtrx]->ielem_dim; j++) {
        De_hat[i + 1][j + 1] = evec[i][j];
      }
    }

    for (i = 0; i < ei[pg->imtrx]->ielem_dim + 1; i++) {
      for (j = 0; j < ei[pg->imtrx]->dof[eqn]; j++) {
        sum = 0.;
        for (k = 0; k < ei[pg->imtrx]->ielem_dim + 1; k++) {
          sum += De_hat[i][k] * Ee_tilde[k][j];
        }
        Ee_hat[i][j] = sum;
      }
    }

    /* Compute the real De_hat */
    memset(De_hat, 0, sizeof(double) * (DIM + 1) * (DIM + 1));
    De_hat[0][0] = vol;
    for (i = 0; i < ei[pg->imtrx]->ielem_dim; i++) {
      De_hat[i + 1][i + 1] = eval[i];
    }

    memset(S, 0, sizeof(double) * (DIM + 1) * (DIM + 1));
    S[0][0] = 1.;
    for (i = 1; i < ei[pg->imtrx]->ielem_dim; i++) {
      S[i + 1][i + 1] = 1.0 - sqrt(eval[0] / eval[i]);
    }

    memset(We, 0, sizeof(double) * (DIM + 1) * (DIM + 1));
    We[0][0] = 1. / vol;
    for (i = 1; i < ei[pg->imtrx]->ielem_dim + 1; i++) {
      We[i][i] = S[i][i] * (2. - S[i][i]) / De_hat[i][i];
    }

    memset(WeEe, 0, sizeof(double) * (DIM + 1) * MDE);
    for (i = 0; i < ei[pg->imtrx]->ielem_dim + 1; i++) {
      for (j = 0; j < ei[pg->imtrx]->dof[eqn]; j++) {
        sum = 0.;
        for (k = 0; k < ei[pg->imtrx]->ielem_dim + 1; k++) {
          sum += We[i][k] * Ee_hat[k][j];
        }
        WeEe[i][j] = sum;
      }
    }

    scaling_over_visc_e = my_multiplier * PS_scaling / visc_e;
    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      for (j = 0; j < ei[pg->imtrx]->dof[eqn]; j++) {
        EeT_We_Ee = 0.;
        for (k = 0; k < ei[pg->imtrx]->ielem_dim + 1; k++) {
          EeT_We_Ee += Ee_hat[k][i] * WeEe[k][j];
        }
        Ce[i][j] = scaling_over_visc_e * (Me[i][j] - EeT_We_Ee);
      }
    }

    /* Now compute and add stabilizing term to residual equation for continuity*/
    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      for (j = 0; j < ei[pg->imtrx]->dof[eqn]; j++) {
        lec->R[LEC_R_INDEX(peqn, i)] += Ce[i][j] * (*esp->P[j]);
      }
    }

    if (pd->v[pg->imtrx][MESH_DISPLACEMENT1] && af->Assemble_Jacobian) {
      for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
        for (j = 0; j < ei[pg->imtrx]->dof[eqn]; j++) {
          if (ipass == 0) {
            R_old[i] += Ce[i][j] * (*esp->P[j]);
          } else {
            R_new[i] += Ce[i][j] * (*esp->P[j]);
          }
        }
      }
    }

  } /* End ipass loop */

  /*Now compute and jacobian contributions */

  if (af->Assemble_Jacobian) {
    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      var = PRESSURE;
      pvar = upd->vp[pg->imtrx][var];
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += Ce[i][j];
      }
    }

    var = MESH_DISPLACEMENT1;

    if (pd->v[pg->imtrx][var]) {
      for (a = 0; a < ei[pg->imtrx]->ielem_dim; a++) {
        var += a;
        /* first the pieces for viscosity */
        for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
          pvar = upd->vp[pg->imtrx][var];
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            /* lec->J[LEC_J_INDEX(peqn,pvar,i,j)] += numerical piece; */
            GOMA_EH(GOMA_ERROR, "need to finish off numerical jacobian for PPSP mesh displ");
          }
        }
      }
    }
    var = VELOCITY1;
    /* need this piece too for d(visc_e)/dv */

    for (a = 0; a < ei[pg->imtrx]->ielem_dim; a++) {
      var += a;
      pvar = upd->vp[pg->imtrx][var];
      for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {

        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          EeT_We_Ee = 0.;
          for (k = 0; k < ei[pg->imtrx]->ielem_dim + 1; k++) {
            EeT_We_Ee += Ee_hat[k][i] * WeEe[k][j];
          }

          lec->J[LEC_J_INDEX(peqn, pvar, i, j)] -= my_multiplier * PS_scaling * d_visc_e_dv[a][j] *
                                                   (Me[i][j] - EeT_We_Ee) / visc_e / visc_e;
        }
      }
    }
  }

  return (1);
}