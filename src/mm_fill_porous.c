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

/* mm_fill_porous -- auxiliary routines helpful in matrix & residual
 *                   assembly for porous media equations
 */

/* Standard include files */
#include "load_field_variables.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

/* GOMA include files */
#define GOMA_MM_FILL_POROUS_C
#include "bc_colloc.h"
#include "el_elm.h"
#include "el_elm_info.h"
#include "exo_struct.h"
#include "mm_as.h"
#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "mm_eh.h"
#include "mm_elem_block_structs.h"
#include "mm_fill_aux.h"
#include "mm_fill_ls.h"
#include "mm_fill_solid.h"
#include "mm_fill_species.h"
#include "mm_fill_terms.h"
#include "mm_fill_util.h"
#include "mm_mp.h"
#include "mm_mp_const.h"
#include "mm_mp_structs.h"
#include "mm_numjac.h"
#include "mm_qtensor_model.h"
#include "mm_shell_util.h"
#include "mm_std_models_shell.h"
#include "rf_allo.h"
#include "rf_bc_const.h"
#include "rf_element_storage_struct.h"
#include "rf_fem.h"
#include "rf_fem_const.h"
#include "rf_mp.h"
#include "std.h"
#include "user_mp.h"

#include "density.h"
#include "mm_fill_porous.h"

/*
 * Global variables defined here. Declared frequently via rf_bc.h
 */

/*********** R O U T I N E S   I N   T H I S   F I L E *************************
 *
 *						-All routines in this
 *				NOTE:		 file are called only by
 *						 mm_fill.c: matrix_fill
 *						 except possibly for static
 *						 functions.
 *
 *       NAME			TYPE			CALLED BY
 *  -----------------------------------------------------------------
 *  assemble_porous_transport        int
 *  load_porous_properties           int
 *  get_porous_part_sat_terms        int
 *  get_porous_fully_sat_terms       int
 *  porous_mass_flux_surf_bc         void
 *  porous_press_flux_surf_bc        void
 *  porous_convection_bc             void
 *  porous_kinematic_bc              void
 *  porous_normal_velocity_bc        void
 *  put_gas_flux_in_pores            void
 *  porous_vapor_equil_bc            void
 *  load_permeability                double
 *  load_saturation                  double
 *  load_cap_pres                    double
 *  load_enthalpy                    void
 *  load_gas_conc                    void
 *  load_gas_conc_flat               void
 *  load_bulk_density                void
 *  load_liq_perm                    void
 *  load_gas_perm                    void
 *  load_gas_diff                    void
 *  load_mass_flux                   void
 *  porous_pressure                  void
 *  porous_pressure_lub              void
 *  sat_darcy_continuous_bc          void
 *
 *******************************************************************************/
void calc_grad_Ywg(int, double *);
/******************************************************************************/
/******************************************************************************/

int assemble_porous_transport(double time, /* present time valuel; KSC           */
                              double tt,   /* parameter to vary time integration from
                                            * explicit (tt = 1) to implicit (tt = 0) */
                              double dt)

/* assemble_porous_transport --
 *
 * assemble terms (Residual &| Jacobian) for porous media eqns
 *
 * input
 * -------
 * 	ei -- pointer to Element Indeces	structure
 *	pd -- pointer to Problem Description	structure
 *	af -- pointer to Action Flag		structure
 *	bf -- pointer to Basis Function		structure
 *	fv -- pointer to Field Variable		structure
 *  fv_old -- pointer to old Field Variable	structure
 *  fv_dot -- pointer to dot Field Variable	structure
 *	cr -- pointer to Constitutive Relation	structure
 *	md -- pointer to Mesh Derivative	structure
 *	me -- pointer to Material Entity	structure
 *
 *  Output
 * --------
 *	a   -- gets loaded up with proper contribution
 *	lec -- gets loaded up with local contributions to resid, Jacobian
 * 	r   -- residual RHS vector
 *
 * Created:	Tue May  2 17:28:00 MST 2000 susubia@sandia.gov
 *
 */
{
  int var, ii, peqn, pvar, eqn, ledof;
  const int dim = pd->Num_Dim;
  int w, w1, p, b, i, j, status = 0;
  const int i_pl = 0; /* counters for porous equation numbers */
  struct Porous_Media_Terms pm_terms;

  dbl mass, advection;
  dbl advection_a, advection_b, advection_c, advection_d, advection_supg;
  dbl diffusion, diff_a, diff_b, diff_c, source;

  /*
   * Galerkin weighting functions for i-th energy residuals
   * and some of their derivatives...
   */
  dbl phi_i, *grad_phi_i;

  /*
   * Petrov-Galerkin weighting functions for i-th residuals
   * and some of their derivatives...
   */
  dbl wt_func;

  /* SUPG variables */
  dbl supg = 0.0;

  /*
   * Interpolation functions for variables and some of their derivatives.
   */
  dbl phi_j;
  dbl h3; /* Product of all 3 scale factors used for */
  /* volume integrals in this coordinate system*/
  dbl dh3dmesh_bj; /* Mesh derivative of same. */
  dbl det_J = 0.0;
  dbl d_det_J_dmeshbj;        /* for specified (b,j) mesh dof */
  dbl dgrad_phi_i_dmesh[DIM]; /* ditto.  */
  dbl wt, wt_total = 0.0;
  int err;

  /*
   * Bail out fast if there's nothing to do...
   */
  if (!pd->e[pg->imtrx][R_POR_LIQ_PRES]) {
    return (status);
  }

  wt = fv->wt; /* Gauss point weight. */
  h3 = fv->h3; /* Differential volume element. */

  /*
   * Store the weighting factor for SUPG in a temporary variable
   */
  if (mp->Porous_wt_funcModel == GALERKIN) {
    supg = 0.0;
  } else if (mp->Porous_wt_funcModel == SUPG) {
    supg = mp->Porous_wt_func;
  }

  /************************************************************************/
  /*                       START OF POROUS ASSEMBLE                       */
  /************************************************************************/

  if (mp->PorousMediaType == POROUS_SATURATED) {
    err = get_porous_fully_sat_terms(&pm_terms, tt, dt);
    GOMA_EH(err, "problem in getting the saturated porous darcy  terms");
    if (neg_elem_volume)
      return (err);
  } else if ((mp->PorousMediaType == POROUS_UNSATURATED ||
              mp->PorousMediaType == POROUS_TWO_PHASE) &&
             !efv->ev_porous_decouple) {
    err = get_porous_part_sat_terms(&pm_terms, tt, dt);
    GOMA_EH(err, "problem in getting the partially-saturated porous  terms");
    if (neg_elem_volume)
      return (err);
  } else if ((mp->PorousMediaType == POROUS_UNSATURATED ||
              mp->PorousMediaType == POROUS_TWO_PHASE) &&
             efv->ev_porous_decouple) {
    err = get_porous_part_sat_terms_decoupled(&pm_terms, tt, dt);
    GOMA_EH(err, "problem in getting the partially-saturated porous  terms, decoupled");
    if (neg_elem_volume)
      return (err);
  }

  /*
   * Residuals_________________________________________________________________
   */

  if (af->Assemble_Residual) {

    /*
     *   START porous equations. This is still a loop and is
     *   somewhat awkward so that we can keep redundant lines of code down.
     *   We will loop to the MAX_PMV number, skipping equations that are
     *   not active.  The test for activity is pd->e[pg->imtrx][var] were
     *   var is one of these types:
     *      R_POR_LIQ_PRES, R_POR_GAS_PRES, R_POR_ENERGY, AND R_POR_POROSITY
     */
    for (w = 0; w < MAX_PMV; w++) {
      eqn = R_POR_LIQ_PRES + w;
      if (pd->e[pg->imtrx][eqn]) {
        /*
         * Store the row location in the local element resid and Jacobina
         */
        peqn = upd->ep[pg->imtrx][eqn];

        /*
         * Store the determinant of the basis function at the local point
         */
        det_J = bf[eqn]->detJ;

        /*
         * Store the common integral weighting factor
         */
        wt_total = wt * det_J * h3;

        /*
         *  Loop over the degrees of freedom in the element corresponding
         *  to liquid-phase pressure unknowns. These are at different nodes in the
         *  element, usually. However, there can be more than one set
         *  of degrees of freedom at each node. Note, this
         *  step doesn't depend upon the porous media equation number, so
         *  we might think about exchanging the order of the loops!
         */
        for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
          ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][i];
          if (ei[pg->imtrx]->active_interp_ledof[ledof]) {
            /*
             *  Here is where we figure out whether the row is to placed in
             *  the normal spot (e.g., ii = i), or whether a boundary condition
             *  require that the volumetric contribution be stuck in another
             *  ldof pertaining to the same variable type.
             *  ->Note, there are no such bc's for porous flow equations.
             */
            ii = i; /*default behavor for all elements */

            /*
             * Store the local values of the basis function and
             * gradient of the basis functions
             */
            phi_i = bf[eqn]->phi[i];
            grad_phi_i = bf[eqn]->grad_phi[i];
            wt_func = phi_i;
            /*
             *  add Petrov-Galerkin terms as necessary
             *  NOTE:  We do not want to make this adjustment for
             *  a porosity equation
             */
            if (supg != 0.0 && eqn == R_POR_LIQ_PRES) {
              wt_func += supg * pm_terms.pi_supg[i];
            }

            mass = 0.0;
            if (pd->TimeIntegration != STEADY) {
              if (pd->e[pg->imtrx][eqn] & T_MASS) {
                if (mp->Porous_Mass_Lump) {
                  mass = -(pmv_ml->Inventory_Solvent_dot[i][w] * wt_func * wt_total *
                           pd->etm[pg->imtrx][eqn][LOG2_MASS]);
                } else {
                  mass = -(pm_terms.Inventory_solvent_dot[w] * wt_func * wt_total *
                           pd->etm[pg->imtrx][eqn][LOG2_MASS]);
                }
              }
            }

            /*
             *   Advection is velocity times gradient of the porous media
             *   unknown variable, possibly multiplied by a density or total
             *   concentration
             */
            advection = 0.0;
            if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
              advection_a = 0.0;
              for (p = 0; p < dim; p++) {
                advection_a += pm_terms.conv_flux[w][p];
              }
              advection_a *= -wt_func;

              advection_b = 0.0;

              advection = advection_a + advection_b;
              advection *= wt_total * pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
            }

            /*
             * Extra advection term from SUPG implementation
             */
            advection_supg = 0.0;
            if (supg > 0.0 && eqn == R_POR_LIQ_PRES) {
              advection_supg = -supg * pm_terms.pi_supg[i] * pm_terms.conv_flux_supg[w] * wt_total *
                               pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
            }

            /*
             * the diffusion term contains all the fluxes that are in
             * a divergence in the mass conservation equation and are
             * integrated by parts in finite elements
             */
            diffusion = 0.0;
            if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {
              for (p = 0; p < VIM; p++) {
                diffusion += grad_phi_i[p] * pm_terms.diff_flux[w][p];
              }
              diffusion *= wt_total * pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
            }
            /*
             * HKM -> Note the addition of a species molecular weight
             *        term is currently done in the source term
             *        calculation section of the code
             */
            source = 0.0;
            if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
              source = (pm_terms.MassSource[w] * wt_func * wt_total *
                        pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)]);
            }

            /*
             *  Sum up all of the individual contributions and store it
             *  in the local element residual vector.
             */
            lec->R[LEC_R_INDEX(peqn, ii)] += mass + advection + advection_supg + diffusion + source;

          } /* if active_dofs */
        }   /* end of loop over equations */
      }     /* If(pd->e[pg->imtrx][eqn]) */
    }       /* loop over PMV equation-types with w */
  }         /* end of assemble residuals */

  /*
   * Jacobian terms...
   */

  if (af->Assemble_Jacobian) {

    /*
     *         START loop over the rows corresponding to different
     *	    porous media equations
     *             w = row porous media index
     *             i = node (i.e., dof) where porous media
     *                 equation is located.
     */

    for (w = 0; w < MAX_PMV; w++) {
      eqn = R_POR_LIQ_PRES + w;
      peqn = upd->ep[pg->imtrx][eqn];
      if (pd->e[pg->imtrx][eqn]) {
        for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
          ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][i];
          if (ei[pg->imtrx]->active_interp_ledof[ledof]) {
            ii = i;

            /*
             * Store the local values of the basis function and
             * gradient of the basis functions
             */
            phi_i = bf[eqn]->phi[i];
            grad_phi_i = bf[eqn]->grad_phi[i];

            /*
             *  add Petrov-Galerkin terms as necessary
             */
            wt_func = phi_i;
            if (supg != 0.0 && w != 2) {
              wt_func += supg * pm_terms.pi_supg[i];
            }

            /*
             * J_pm_pm derivative of residual pieces w.r.t. the porous media
             *       unknowns
             */
            for (w1 = 0; w1 < MAX_PMV; w1++) {
              var = POR_LIQ_PRES + w1;
              if (pd->v[pg->imtrx][var]) {
                pvar = upd->vp[pg->imtrx][var];

                /*
                 * Add terms for mass lumping of the time derivative
                 * corresponding to block diagonal contributions.
                 */
                if (mp->Porous_Mass_Lump) {
                  mass = -pmv_ml->d_Inventory_Solvent_dot_dpmv[i][w][w1] * wt_func * wt_total *
                         pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
                  lec->J[LEC_J_INDEX(peqn, pvar, ii, ii)] += mass;
                }

                /*
                 * Loop over column degrees of freedom
                 */
                for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                  phi_j = bf[var]->phi[j];

                  mass = 0.0;
                  if (pd->TimeIntegration != STEADY) {
                    /*
                     * Include the dependence of the capacitance on the
                     * Liquid and Gas pressure unknowns. The liquid will
                     * generally be compressible.
                     */
                    if (pd->e[pg->imtrx][eqn] & T_MASS) {

                      if (mp->Porous_Mass_Lump) {
                        /*
                         * SUPG Section
                         */
                        mass = 0.0;
                        if (var == POR_LIQ_PRES && supg != 0.0 && eqn == POR_LIQ_PRES) {
                          /*
                           * Include the dependence of the SUPG basis function
                           * on Pl
                           */
                          mass = -supg * pm_terms.d_pi_supg_dpmv[i][i_pl][j] *
                                 pmv_ml->Inventory_Solvent_dot[i][w] * wt_total *
                                 pd->etm[pg->imtrx][eqn][LOG2_MASS];
                        }
                      } else {
                        mass = -pm_terms.d_Inventory_solvent_dot_dpmv[w][w1][j] * wt_func *
                               wt_total * pd->etm[pg->imtrx][eqn][(LOG2_MASS)];

                        /*
                         * SUPG Section
                         */
                        if (var == POR_LIQ_PRES && supg != 0.0 && eqn == POR_LIQ_PRES) {
                          /*
                           * Include the dependence of the SUPG basis function
                           * on Pl
                           */
                          mass -= supg * pm_terms.d_pi_supg_dpmv[i][i_pl][j] *
                                  pm_terms.Inventory_solvent_dot[w] * wt_total *
                                  pd->etm[pg->imtrx][eqn][LOG2_MASS];
                        }
                      }
                    }
                  }

                  advection = 0.0;
                  if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
                    advection_a = 0.0;
                    for (p = 0; p < dim; p++) {
                      advection_a -= pm_terms.d_conv_flux_dpmv[w][p][w1][j];
                    }
                    advection_a *= wt_func;

                    advection_b = 0.0;

                    advection = advection_a + advection_b;

                    /*
                     * SUPG Section
                     */
                    if (var == POR_LIQ_PRES && supg != 0.0 && eqn == POR_LIQ_PRES) {
                      /*
                       * Include the dependence of the SUPG basis function
                       * on Pl
                       */
                      for (p = 0, advection_c = 0.0; p < dim; p++) {
                        advection_c += pm_terms.conv_flux[w][p];
                      }
                      advection -= supg * pm_terms.d_pi_supg_dpmv[i][i_pl][j] * advection_c;
                    }
                    advection *= wt_total * pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
                  }

                  /*
                   * Extra advection term from SUPG implementation
                   */
                  advection_supg = 0.0;
                  if (supg > 0.0) {

                    if (var == POR_LIQ_PRES && supg != 0.0 && eqn == POR_LIQ_PRES) {
                      advection_supg =
                          (-pm_terms.d_pi_supg_dpmv[i][i_pl][j] * pm_terms.conv_flux_supg[w] -
                           pm_terms.pi_supg[i] * pm_terms.d_conv_flux_supg_dpmv[i_pl][i_pl][j]);
                    }
                    advection_supg *= supg * wt_total * pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
                  }

                  diffusion = 0.0;
                  if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {
                    for (p = 0; p < VIM; p++) {
                      diffusion += grad_phi_i[p] * pm_terms.d_diff_flux_dpmv[w][p][w1][j];
                    }
                    diffusion *= wt_total * pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
                  }

                  source = 0.0;
                  if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
                    source = (pm_terms.d_MassSource_dpmv[w][w1][j] * wt_func * wt_total *
                              pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)]);
                    if (var == POR_LIQ_PRES && supg != 0.0 && eqn == POR_LIQ_PRES) {
                      source += (pm_terms.MassSource[w] * pm_terms.d_pi_supg_dpmv[i][i_pl][j] *
                                 wt_total * pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)]);
                    }
                  }

                  lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] +=
                      mass + advection + advection_supg + diffusion + source;
                }
              }
            }

            /*
             * J_pm_d  sensitivity of porous medium equations
             *         w.r.t. mesh displacement
             */

            for (b = 0; b < dim; b++) {
              var = MESH_DISPLACEMENT1 + b;
              if (pd->v[pg->imtrx][var]) {
                pvar = upd->vp[pg->imtrx][var];
                for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                  phi_j = bf[var]->phi[j];
                  d_det_J_dmeshbj = bf[var]->d_det_J_dm[b][j];

                  dh3dmesh_bj = fv->dh3dq[b] * phi_j;

                  mass = 0.0;
                  if (pd->TimeIntegration != STEADY) {
                    if (pd->e[pg->imtrx][eqn] & T_MASS) {
                      mass = -(pm_terms.Inventory_solvent_dot[w] * wt_func * wt *
                               pd->etm[pg->imtrx][eqn][(LOG2_MASS)] *
                               (h3 * d_det_J_dmeshbj + dh3dmesh_bj * det_J));
                    }
                  }

                  advection = 0.0;
                  if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
                    /*
                     * Two parts:
                     *
                     *	Int( (v- xdot).d(Vc)/dmesh h3 |Jv| )
                     *	Int( v.Vc d(h3|Jv|)/dmesh )
                     */
                    advection_a = 0.0;
                    for (p = 0; p < dim; p++) {
                      advection_a += pm_terms.d_conv_flux_dmesh[w][p][b][j];
                    }
                    advection_a *= -wt_func * h3 * det_J * wt;

                    advection_b = 0.0;
                    for (p = 0; p < VIM; p++) {
                      advection_b += pm_terms.conv_flux[w][p];
                    }
                    advection_b *= -wt_func;

                    advection_b *= (h3 * d_det_J_dmeshbj + dh3dmesh_bj * det_J) * wt;

                    advection_d = 0.0;

                    advection = advection_a + advection_b + advection_d;
                    advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
                  }

                  diffusion = 0.0;
                  if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {
                    /*
                     * Three parts:
                     *  diff_a = Int(d(grad_phi_i)/dmesh.q h3 |Jv|)
                     *  diff_b = Int(grad_phi_i.d(q)/dmesh h3 |Jv|)
                     *  diff_c = Int(grad_phi_i.q d(h3 |Jv|)/dmesh)
                     */

                    diff_a = 0.0;
                    for (p = 0; p < VIM; p++) {
                      dgrad_phi_i_dmesh[p] = bf[eqn]->d_grad_phi_dmesh[i][p][b][j];
                      diff_a += dgrad_phi_i_dmesh[p] * pm_terms.diff_flux[w][p];
                    }
                    diff_a *= h3 * det_J * wt;

                    diff_b = 0.0;
                    for (p = 0; p < VIM; p++) {
                      diff_b += bf[eqn]->grad_phi[i][p] * pm_terms.d_diff_flux_dmesh[w][p][b][j];
                    }
                    diff_b *= h3 * det_J * wt;

                    diff_c = 0.0;
                    for (p = 0; p < VIM; p++) {
                      diff_c += bf[eqn]->grad_phi[i][p] * pm_terms.diff_flux[w][p];
                    }
                    diff_c *= (h3 * d_det_J_dmeshbj + dh3dmesh_bj * det_J) * wt;

                    diffusion = diff_a + diff_b + diff_c;
                    diffusion *= pd->etm[pg->imtrx][eqn][LOG2_DIFFUSION];
                  }

                  source = 0.0;
                  if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
                    source = pm_terms.MassSource[w];
                    source *= (h3 * d_det_J_dmeshbj + dh3dmesh_bj * det_J) * wt * phi_i;
                    source += pm_terms.d_MassSource_dmesh[w][b][j] * det_J * h3 * wt * phi_i;
                    source *= pd->etm[pg->imtrx][eqn][LOG2_SOURCE];
                  }

                  lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += mass + advection + diffusion + source;
                }
              }
            }

            /*
             * J_pm_T  sensitivity of porous medium equations w.r.t. temperature
             *
             *         Only compute this if the energy equation is NOT solved;
             *         if it is solved, the sensitivity is included in the above.
             *
             */
            var = TEMPERATURE;
            if (pd->e[pg->imtrx][eqn] && pd->v[pg->imtrx][var] && !pd->e[pg->imtrx][R_POR_ENERGY]) {
              pvar = upd->vp[pg->imtrx][var];
              for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                phi_j = bf[var]->phi[j];

                advection = 0.0;
                if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
                  for (p = 0; p < VIM; p++) {
                    advection += pm_terms.d_conv_flux_dT[w][p][j];
                  }
                  advection *= -wt_func;

                  advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)] * h3 * det_J * wt;
                }

                diffusion = 0.0;
                if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {
                  for (p = 0; p < VIM; p++) {
                    diffusion += grad_phi_i[p] * pm_terms.d_diff_flux_dT[w][p][j];
                  }
                  diffusion *= h3 * det_J * wt * pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
                }

                source = 0.0;
                if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
                  source = pm_terms.d_MassSource_dT[w][j] * det_J * h3 * wt * phi_i;
                  source *= pd->etm[pg->imtrx][eqn][LOG2_SOURCE];
                }

                lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += advection + diffusion + source;
              }
            }
            /*
             * J_pm_sink_mass  sensitivity of porous medium equations w.r.t.
             *                 absorbed liquid mass
             *
             *         Only compute this if the sink mass equation is solved;
             *
             */
            var = POR_SINK_MASS;
            if (pd->v[pg->imtrx][var]) {
              source = 0.0;
              pvar = upd->vp[pg->imtrx][var];
              if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
                for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                  source = pm_terms.d_MassSource_dSM[w][j] * det_J * h3 * wt * phi_i;
                  source *= pd->etm[pg->imtrx][eqn][LOG2_SOURCE];
                  lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += source;
                }
              }

              if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {

                for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                  advection = 0.0;
                  for (p = 0; p < VIM; p++) {
                    advection += pm_terms.d_conv_flux_dSM[w][p][j];
                  }
                  advection *= -wt_func;

                  advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)] * h3 * det_J * wt;
                  lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += advection;
                }
              }

              if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {

                for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                  diffusion = 0.0;
                  for (p = 0; p < VIM; p++) {
                    diffusion += pm_terms.d_diff_flux_dSM[w][p][j];
                  }
                  diffusion *= -wt_func;

                  diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)] * h3 * det_J * wt;
                  lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] += diffusion;
                }
              }
            }
          }
        }
      }
    }
  }

  return (status);
} /* end of assemble_porous_transport */
/****************************************************************************/
/******************************************************************************/

#ifdef HAVE_AGM_MODEL
#include "AGM_model.c"
#else
int assemble_pore_sink_mass(double time, /* present time valuel; KSC           */
                            double tt,   /* parameter to vary time integration from
                                          * explicit (tt = 1) to implicit (tt = 0) */
                            double dt)

/* assemble_pore_sink_mass --
 *
 * assemble terms (Residual &| Jacobian) for porous mass-sink equation
 *
 * input
 * -------
 * 	ei -- pointer to Element Indeces	structure
 *	pd -- pointer to Problem Description	structure
 *	af -- pointer to Action Flag		structure
 *	bf -- pointer to Basis Function		structure
 *	fv -- pointer to Field Variable		structure
 *  fv_old -- pointer to old Field Variable	structure
 *  fv_dot -- pointer to dot Field Variable	structure
 *
 *  Output
 * --------
 *	a   -- gets loaded up with proper contribution
 *	lec -- gets loaded up with local contributions to resid, Jacobian
 * 	r   -- residual RHS vector
 *
 * Created:	29 August 2005 by PRS
 *
 */
{

  int var, peqn, pvar, eqn;
  const int dim = pd->Num_Dim;
  int b, i, j, status = 0;

  double MassSource = 0.0;
  double d_MassSource[MAX_VARIABLE_TYPES + MAX_CONC][MDE];

  dbl mass, source;

  /*
   * Galerkin weighting functions for i-th energy residuals
   * and some of their derivatives...
   */
  dbl phi_i;

  /*
   * Interpolation functions for variables and some of their derivatives.
   */
  dbl phi_j;
  dbl h3;       /* Product of all 3 scale factors used for */
                /* volume integrals in this coordinate system*/
  dbl dh3dmesh; /* Mesh derivative of same. */
  dbl det_J;
  dbl d_det_J_dmesh; /* for specified (b,j) mesh dof */
  dbl wt, wt_total;

  dbl sink_dot[MDE];

  /*
   * Bail out fast if there's nothing to do...
   */
  if (!pd->e[pg->imtrx][R_POR_SINK_MASS]) {
    return (status);
  }

  /* Get the mass term - mass lumping */
  for (i = 0; i < ei[pg->imtrx]->dof[R_POR_SINK_MASS]; i++) {
    sink_dot[i] = *esp_dot->sink_mass[i];
  }

  wt = fv->wt;
  h3 = fv->h3;

  /* Get the MassSouece and their Jacobian sensitivities */
  MassSource = por_mass_source_model(d_MassSource);

  /*
   * Residuals_________________________________________________________________
   */

  if (af->Assemble_Residual) {
    eqn = R_POR_SINK_MASS;
    if (pd->e[pg->imtrx][eqn]) {
      peqn = upd->ep[pg->imtrx][eqn];

      det_J = bf[eqn]->detJ;

      wt_total = wt * det_J * h3;

      for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
        phi_i = bf[eqn]->phi[i];

        mass = 0.0;
        if (pd->TimeIntegration != STEADY) {
          if (pd->e[pg->imtrx][eqn] & T_MASS) {
            mass += (sink_dot[i] * phi_i * wt_total * pd->etm[pg->imtrx][eqn][LOG2_MASS]);
          }
        }

        source = 0.0;
        if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
          source += phi_i * mp->density * MassSource;
          source *= wt_total * pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
        }

        lec->R[LEC_R_INDEX(peqn, i)] += mass + source;
      }
    }
  }

  /*
   * Jacobian terms...
   */

  if (af->Assemble_Jacobian) {
    eqn = R_POR_SINK_MASS;
    peqn = upd->ep[pg->imtrx][eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      phi_i = bf[eqn]->phi[i];

      det_J = bf[eqn]->detJ;

      wt_total = wt * det_J * h3;

      /*
       * J_sink_mass_d_sink_mass
       *  derivative of residual pieces w.r.t. the sink mass
       */

      var = POR_SINK_MASS;
      if (pd->v[pg->imtrx][var]) {
        pvar = upd->vp[pg->imtrx][var];

        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          phi_j = bf[var]->phi[j];

          mass = 0.0;

          if (pd->TimeIntegration != STEADY) {
            if (pd->e[pg->imtrx][eqn] & T_MASS) {
              mass += delta(i, j) * (1 + 2. * tt) / dt;
              mass *= phi_i * wt_total * pd->etm[pg->imtrx][eqn][LOG2_MASS];
            }
          }

          source = 0.0;
          if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
            source += phi_i * mp->density * d_MassSource[var][j];
            source *= wt_total * pd->etm[pg->imtrx][eqn][LOG2_SOURCE];
          }

          lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += mass + source;
        }
      }

      /*
       * J_sink_mass_d_por_liq_pres
       *  derivative of residual pieces w.r.t. the pore liquid pressure
       */

      var = POR_LIQ_PRES;
      if (pd->v[pg->imtrx][var]) {
        pvar = upd->vp[pg->imtrx][var];
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {

          source = 0.0;
          if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
            source += phi_i * mp->density * d_MassSource[var][j];
            source *= wt_total * pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
          }

          lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += source;
        }
      }

      /*
       * J_sink_mass_d_shell_press_open
       *  derivative of residual pieces w.r.t. the shell porous open pressure
       */

      var = SHELL_PRESS_OPEN;
      if (pd->v[pg->imtrx][var]) {
        pvar = upd->vp[pg->imtrx][var];
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {

          source = 0.0;
          if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
            source += phi_i * mp->density * d_MassSource[var][j];
            source *= wt_total * pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
          }

          lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += source;
        }
      }

      var = SHELL_SAT_1;
      if (pd->v[pg->imtrx][var]) {
        pvar = upd->vp[pg->imtrx][var];
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {

          source = 0.0;
          if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
            source += phi_i * mp->density * d_MassSource[var][j];
            source *= wt_total * pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
          }

          lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += source;
        }
      }

      /*
       * J_sink_mass_d_mesh
       *  derivative of residual pieces w.r.t. mesh displacement
       */

      var = MESH_DISPLACEMENT1;

      if (pd->v[pg->imtrx][var]) {
        for (b = 0; b < dim; b++) {
          var = MESH_DISPLACEMENT1 + b;
          pvar = upd->vp[pg->imtrx][var];

          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];

            d_det_J_dmesh = bf[var]->d_det_J_dm[b][j];
            dh3dmesh = fv->dh3dq[b] * phi_j;

            mass = 0.0;
            if (pd->TimeIntegration != STEADY) {
              if (pd->e[pg->imtrx][eqn] & T_MASS) {
                mass += sink_dot[i] * phi_i * pd->etm[pg->imtrx][eqn][LOG2_MASS] *
                        (h3 * d_det_J_dmesh + dh3dmesh * det_J) * wt;
              }
            }

            source = 0.0;
            if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
              source += phi_i * MassSource * (h3 * d_det_J_dmesh + dh3dmesh * det_J);
              source += phi_i * d_MassSource[var][j] * h3 * det_J;
              source *= wt * pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
            }

            lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += mass + source;
          }
        }
      }

    } /* End of loop over i-th equation */
  }   /* End of if Assemble_Jacobian */

  //   GOMA_EH(GOMA_ERROR,"if you got here you need to contact P. R. Schunk.  Pore-Sink-mass");
  return status;
}
#endif
/***************************************************************************/
/***************************************************************************/

/***************************************************************************/
/***************************************************************************/
int load_porous_properties(void)

/***********************************************************************
 * load_porous_properties --
 *
 *  Calculate properties and state variables in a porous
 *  media from the unknowns used in the problem. This routine is called
 *  at a gauss point level, before the routine,
 *  assemble_porous_transport(), is called.
 *
 *  input:
 *  -----
 *
 *  assume that load_fv and load_fv_grads and load_fv_mesh_derivs have
 *  all been called already, so unknowns and their sensitivies are
 *  known at this gauss point
 *
 *
 *  output:
 *  -------
 *
 *   Calculates the porosity, permeability, saturation, capillary
 *   pressure, etc. along with gradients and time derivatives
 *   of all the relevant quantities
 *
 *      mp->porosity
 *      mp->permeability
 *
 **********************************************************************/
{
  int a, b, w, w1, i_por_ev;
  double cap_pres = 1e12, p_gas_star = 0., saturation, porosity = 1e12;
  double d_cap_pres[4], n_pow;
  const int i_pl = 0, i_pg = 1, i_pore = 2, i_pe = 3;
  int siz;
  int MTV; /* maximum number of variables types plus species unknowns */

  static int is_initialized = FALSE; /* This static flag is set TRUE after the first pass through.
                                      * This is to reduce the tremendous overhead associated with
                                      * unnecessarily reinializing the mp, pmt, pmv structures
                                      */

  if (mp->PorousMediaType == CONTINUOUS)
    return 0; /* shouldn't get here */
  // if (mp->PorousMediaType == POROUS_SHELL) return 0;

  MTV = MAX_VARIABLE_TYPES + MAX_CONC;

  /*
   * There are several types of media classifications for which GOMA is designed:
   *  CONTINUOUS: a liquid, solid or gas phase which can be treated as a continuum
   *  POROUS_SATURATED: a porous medium which is saturated with a solvent or gas -
   *      in this type of medium we solve for the liquid pressure which drives flow
   *      through the pores
   *  POROUS_UNSATURATED: a porous medium which has both liquid and gas in the pores, but
   *      the gas is at a uniform pressure so we only solve for the capillary pressure
   *      of the liquid (there may or may not be evaporation of the liquid)
   *  POROUS_SHELL_UNSATURATED: Same as POROUS_UNSATURATED, but on a shell
   *  POROUS_TWO_PHASE: a porous medium which has both liquid and gas in the pores, and
   *      we need to solve for both the liquid pressure and gas pressure with local
   *      equilibrium between the phases
   *  POROUS_BRINKMAN: a saturated porous medium in which shear stresses may be important -
   *      we solve an augmented form of the Navier-Stokes equations with extra parameters
   *      account for the porous nature of the medium.
   *
   * This file directs the calculation of the properties needed for each type of porous medium.
   */

  /*
   * Collect and Store the porosity into the materials property structure
   */
  if (!is_initialized) {
    siz = MTV * sizeof(double);
    memset(mp->d_porosity, 0, siz);
  }

  if (mp->PorosityModel == CONSTANT) {
    porosity = mp->porosity;
    if (pd->TimeIntegration == TRANSIENT) {
      mp_old->porosity = mp->porosity;
      if (pd->v[pg->imtrx][POR_POROSITY])
        mp_old->d_porosity[POR_POROSITY] = 0.0;
    }
  } else if (mp->PorosityModel == EXTERNAL_FIELD) {
    i_por_ev = mp->porosity_external_field_index;
    GOMA_EH(i_por_ev, "Porosity external field not found!");
    porosity = fv->external_field[i_por_ev];
    // mp->porosity = fv->porosity;
    if (pd->TimeIntegration == TRANSIENT) {
      mp_old->porosity = fv->porosity;
    }
  } else if (mp->PorosityModel == DEFORM &&
             (cr->MeshMotion == LAGRANGIAN || cr->MeshMotion == DYNAMIC_LAGRANGIAN)) {
    if (mp->PorousMediaType == POROUS_BRINKMAN) {
      GOMA_EH(GOMA_ERROR, "Variable porosity not allowed with brinkman");
    }
    if (!pd->v[pg->imtrx][POR_POROSITY])
      GOMA_EH(GOMA_ERROR, "Can't solve for DEFORM w/o solving for porosity");
    /* two ways to calculate porosity */
    porosity = mp->porosity = fv->porosity;
    if (pd->TimeIntegration == TRANSIENT) {
      mp_old->porosity = fv_old->porosity;
      mp_old->d_porosity[POR_POROSITY] = 1.0;
    }
    mp->d_porosity[POR_POROSITY] = 1.;
  } else {
    GOMA_EH(GOMA_ERROR, "No models for porosity");
  }

  /*
   * Calculate permeability
   */
  if (!is_initialized) {
    siz = MTV * sizeof(double);
    memset(mp->d_permeability, 0, siz);
  }

  /*  This is a Brinkman permeability model. */
  /*  It is currently set in mm_fill_terms.c by
   *  the routine solidification_permeability in
   *  mm_std_models.c
   */
  if (mp->PermeabilityModel == CONSTANT || mp->PermeabilityModel == K_TENSOR ||
      mp->PermeabilityModel == SOLIDIFICATION) {
    /* do nothing */
  } else if (mp->PermeabilityModel == ORTHOTROPIC || mp->PermeabilityModel == KC_TENSOR ||
             mp->PermeabilityModel == SM_TENSOR) {
    load_permeability_tensor();
  } else if (mp->PermeabilityModel == EXTERNAL_FIELD) {
    // For now we will load this into tensor form.
    i_por_ev = mp->perm_external_field_index;
    GOMA_EH(i_por_ev, "Permeability external field not found!");
    mp->permeability = mp->u_permeability[0] * fv->external_field[i_por_ev];
    if (pd->TimeIntegration == TRANSIENT) {
      mp_old->permeability = fv->external_field[i_por_ev];
    }
    load_permeability_tensor();
  } else {
    if (mp->PorousMediaType == POROUS_BRINKMAN) {
      GOMA_WH(GOMA_ERROR, "Variable permeability not allowed with brinkman");
    }
    (void)load_permeability();
  }

  /*
   * Early return for simple porosity models
   */
  if (mp->PorousMediaType == POROUS_SATURATED) {
    mp->saturation = 1.0;
    pmv->cap_pres = 0.0;
    pmv->liq_Xvol_solvents[i_pl] = 1.0;
    pmv->bulk_density[i_pl] = mp->density * mp->porosity;

    /*
     * Calculate the Darcy convection velocity of the liquid phase
     */
    if (mp->PermeabilityModel != K_TENSOR && mp->PermeabilityModel != ORTHOTROPIC &&
        mp->PermeabilityModel != KC_TENSOR && mp->PermeabilityModel != SM_TENSOR) {

      if (cr->PorousFluxModel == DARCY_FICKIAN) {
        for (a = 0; a < VIM; a++) {
          pmv->liq_darcy_velocity[a] =
              -mp->permeability / mp->viscosity * (fv->grad_p_liq[a] - mp->momentum_source[a]);
        }
      } else if (cr->PorousFluxModel == POWERLAW_DARCY_FICKIAN) {
        n_pow = 1. / 0.5;
        for (a = 0; a < VIM; a++) {
          pmv->liq_darcy_velocity[a] = -mp->permeability / mp->viscosity *
                                       (pow(fv->grad_p_liq[a], n_pow) - mp->momentum_source[a]);
        }
      }
    } else {
      if (cr->PorousFluxModel == DARCY_FICKIAN) {
        for (a = 0; a < VIM; a++) {
          pmv->liq_darcy_velocity[a] = 0.0;
          for (b = 0; b < VIM; b++) {
            pmv->liq_darcy_velocity[a] -= mp->perm_tensor[a][b] / mp->viscosity *
                                          (fv->grad_p_liq[b] - mp->momentum_source[b]);
          }
        }
      } else {
        GOMA_EH(GOMA_ERROR, "No other PorousFluxModel allowed for tensor perm");
      }
    }
    return 0;
  }
  if (mp->PorousMediaType == POROUS_BRINKMAN) {
    mp->saturation = 1.0;
    return 0;
  }

  /* MEDIUM IS PARTIALLY-SATURATED
   *  need to calculate the necessary properties for unsaturated
   *  or two-phase calculations in porous medium:
   *
   *  UNSATURATED: medium is filled with liquid and air, but air
   *               is assumed to be at atmospheric pressure everywhere
   *               so we just solve for capillary pressure
   *               which is atmospheric pressure minus the liquid pressure
   *  TWO_PHASE:   medium is filled with liquid and air, and we solve
   *               for both liquid pressure and air pressure
   *
   *
   *  Calculate capillary pressure and it's sensitivity to liquid
   *  and gas pressures
   *
   *  Very important here that p_gas_star, which comes from the Porous
   *  Gas Constants card is set consistently with all other pressures. It
   *  may heavily affect the saturationif not.  (PRS: 6/5/01)
   */
  if (mp->PorousMediaType == POROUS_UNSATURATED) {
    p_gas_star = mp->u_porous_gas_constants[3];
    cap_pres = pmv->cap_pres = p_gas_star - fv->p_liq;
    d_cap_pres[i_pl] = -1.;
    if (pd->TimeIntegration == TRANSIENT) {
      pmv_old->cap_pres = p_gas_star - fv_old->p_liq;
    }
  } else if (mp->PorousMediaType == POROUS_SHELL_UNSATURATED) {
    dbl p_liq_open = 0, p_liq_open_old = 0;
    p_gas_star = mp->u_porous_gas_constants[3];
    if (pd->v[pg->imtrx][SHELL_PRESS_OPEN]) {
      p_liq_open = fv->sh_p_open;
      p_liq_open_old = fv_old->sh_p_open;
    } else if (pd->v[pg->imtrx][SHELL_PRESS_OPEN_2]) {
      p_liq_open = fv->sh_p_open_2;
      p_liq_open_old = fv_old->sh_p_open_2;
    }

    cap_pres = pmv->cap_pres = p_gas_star - p_liq_open;
    d_cap_pres[i_pl] = -1.;
    if (pd->TimeIntegration == TRANSIENT) {
      pmv_old->cap_pres = p_gas_star - p_liq_open_old;
    }
  } else if (mp->PorousMediaType == POROUS_TWO_PHASE) {

    cap_pres = pmv->cap_pres = fv->p_gas - fv->p_liq;
    d_cap_pres[i_pl] = -1.;
    d_cap_pres[i_pg] = 1.;
    if (pd->TimeIntegration == TRANSIENT) {
      pmv_old->cap_pres = fv_old->p_gas - fv_old->p_liq;
    }

  } else
    GOMA_EH(GOMA_ERROR, "Illegal medium type");

  /*
   * this is a big mess - for each property of the porous media, we need
   * to calculate it's sensitivity to all the unknowns, and the second
   * derivative too - because we need to use the first derivative in the
   * chain rule to get gradients of the property - then we need the second
   * derivative to get sensitivities
   *
   * Thanks to Sam, Harry, Allen, and me, this is just a little better now. (5/01)
   */

  /* initialize the saturation derivatives */
  if (!is_initialized) {
    siz = MTV * sizeof(double);
    memset(mp->d_saturation, 0, siz);
    memset(mp_old->d_saturation, 0, siz);
    siz *= MTV;
    memset(mp->d_d_saturation, 0, siz);
  }

  /* Calculate saturation and it's sensitivity to liquid and gas
   * pressures and porosity
   *   - may also be sensitive to liquid conc. of other species, temperature
   *     or liquid density
   */
  if (mp->SaturationModel == CONSTANT) {
    saturation = mp->saturation;
    if (pd->TimeIntegration == TRANSIENT) {
      mp_old->saturation = mp->saturation;
    }
  } else {
    saturation = load_saturation(porosity, cap_pres, d_cap_pres);
  }

  if (pd->e[pg->imtrx][R_POR_ENERGY]) {

    if (!is_initialized) {
      siz = 3 * MTV * sizeof(double);
      memset(pmv->d_enthalpy, 0, siz);
      memset(pmv_old->d_enthalpy, 0, siz);

      siz *= MTV;
      memset(pmv->d_d_enthalpy, 0, siz);
    }

    /* Calculate enthalpy and it's sensitivity to liquid and gas
     * pressures and temperature
     *   - may also be sensitive to liquid conc. of other species
     */
    load_enthalpy(saturation, p_gas_star);
  }

  /*
   * initialize the liquid & gas phase concentration derivatives
   */
  if (!is_initialized) {
    /* this section really should be masked by some if( pd->v[pg->imtrx][ var ] ).  But which one ?.
     * I cannot tell */
    siz = MAX_PMV * sizeof(double);
    memset(pmv->liq_Xvol_solvents, 0, siz);
    memset(pmv_old->liq_Xvol_solvents, 0, siz);
    memset(pmv->gas_density_solvents, 0, siz);
    memset(pmv_old->gas_density_solvents, 0, siz);

    siz *= MTV;
    memset(pmv->d_liq_Xvol_solvents, 0, siz);
    memset(pmv_old->d_liq_Xvol_solvents, 0, siz);
    memset(pmv->d_gas_density_solvents, 0, siz);
    memset(pmv_old->d_gas_density_solvents, 0, siz);

    siz *= MTV;
    memset(pmv->d_d_gas_density_solvents, 0, siz);
    memset(pmv->d_d_liq_Xvol_solvents, 0, siz);

    siz = MAX_PMV * sizeof(double);
    memset(pmv->d_rhog, 0, siz);
    memset(pmv->d_Ywg, 0, siz);
    memset(pmv->d_Yag, 0, siz);

    siz *= MAX_PMV;
    memset(pmv->d_drhog, 0, siz);
    memset(pmv->d_dYwg, 0, siz);
    memset(pmv->d_dYag, 0, siz);
  }

  /* calculate the gas phase concentrations from local equilibrium assumption - the Kelvin
   *   equation, assumption of ideal gas behavior, and any activity effects
   *   - sensitive to pressures, temperature, and concentration
   *
   * the gas phase contributions to the energy equation (if solving it) are also calculated
   *
   */

  if (mp->PorousVaporPressureModel[i_pl] == FLAT) {
    /* load_gas_conc_flat(porosity, cap_pres, d_cap_pres); */
    load_gas_conc_EOS(porosity, cap_pres, d_cap_pres);
  } else if (mp->PorousVaporPressureModel[i_pl] == NON_VOLATILE) {
    /* do nothing because the gas phase concentration of solvent is zero */
  } else {
    load_gas_conc(porosity, cap_pres, d_cap_pres);
  }

  /* calculate liquid volume fraction of all species - currently pure solvent */
  /* initialize the liquid volume fraction derivatives */
  /* gas is insoluble in liquid */
  if (mp->PorousMediaType == POROUS_TWO_PHASE) {
    pmv->liq_Xvol_solvents[i_pg] = 0.0;
    for (w = POR_LIQ_PRES; w < POR_LIQ_PRES + MAX_PMV; w++)
      pmv->d_liq_Xvol_solvents[i_pg][w] = 0.0;
    if (pd->TimeIntegration == TRANSIENT) {
      pmv_old->liq_Xvol_solvents[i_pg] = 0.0;
    }
  }

  /* solid is insoluble in liquid */
  if (pd->e[pg->imtrx][R_POR_POROSITY]) {
    pmv->liq_Xvol_solvents[i_pore] = 0.0;
    for (w = POR_LIQ_PRES; w < POR_LIQ_PRES + MAX_PMV; w++)
      pmv->d_liq_Xvol_solvents[i_pore][w] = 0.0;
    if (pd->TimeIntegration == TRANSIENT) {
      pmv_old->liq_Xvol_solvents[i_pore] = 0.0;
    }
  }

  /*
   * What's left over is solvent.
   * Currently the volume fraction of solvent in the solvent
   * phase is equal to 1.0. Later, when we make the solvent
   * phase multicomponent, this might change.
   */
  pmv->liq_Xvol_solvents[i_pl] = 1.0;
  for (w = POR_LIQ_PRES; w < POR_LIQ_PRES + MAX_PMV; w++)
    pmv->d_liq_Xvol_solvents[i_pl][w] = 0.0;
  if (pd->TimeIntegration == TRANSIENT) {
    pmv_old->liq_Xvol_solvents[i_pl] = 1.0;
  }

  /*
   * Since porous equation assembly loops over components,
   * adding in contributions from the various phases, we
   * exploit this in assembling the porous energy equation.
   * Thus, the liq_Xvol_solvent for the energy equation
   * corresponds to the enthalpy of the liquid.
   *
   */
  if (pd->e[pg->imtrx][R_POR_ENERGY]) {
    pmv->liq_Xvol_solvents[i_pe] = pmv->enthalpy[0];
    for (w = POR_LIQ_PRES; w < POR_LIQ_PRES + MAX_PMV; w++) {
      pmv->d_liq_Xvol_solvents[i_pe][w] = pmv->d_enthalpy[0][w];
      for (w1 = POR_LIQ_PRES; w1 < POR_LIQ_PRES + MAX_PMV; w1++)
        pmv->d_d_liq_Xvol_solvents[i_pe][w][w1] = pmv->d_d_enthalpy[0][w][w1];
    }
    if (pd->TimeIntegration == TRANSIENT) {
      pmv_old->liq_Xvol_solvents[i_pe] = pmv_old->enthalpy[0];
      for (w = POR_LIQ_PRES; w < POR_LIQ_PRES + MAX_PMV; w++)
        pmv_old->d_liq_Xvol_solvents[i_pe][w] = pmv_old->d_enthalpy[0][w];
    }
  }

  /* initialize the bulk mass concentration derivatives */
  if (!is_initialized) {
    siz = MAX_PMV * sizeof(double);
    memset(pmv->bulk_density, 0, siz);
    memset(pmv_old->bulk_density, 0, siz);

    siz *= MTV;
    memset(pmv->d_bulk_density, 0, siz);
    memset(pmv_old->d_bulk_density, 0, siz);

    siz *= MTV;
    memset(pmv->d_d_bulk_density, 0, siz);
  }

  /*
   *   Calculate bulk mass concentrations of all species.
   *    - The bulk concentration includes the mass of the species
   *      in both gas and liquid phases
   */
  load_bulk_density(porosity, cap_pres, saturation, d_cap_pres);

  /* initialize the relative permeabilities derivatives */

  if (!is_initialized) {
    siz = MTV * sizeof(double);
    memset(mp->d_rel_liq_perm, 0, siz);
    memset(mp->d_rel_gas_perm, 0, siz);
    memset(mp_old->d_rel_liq_perm, 0, siz);
    memset(mp_old->d_rel_gas_perm, 0, siz);
  }

  /*
   *  Calculate relative permeabilities in the liquid and gas phases
   *  These are stored in the material property structure
   *
   * NOTE - actually calculate krel/viscosity
   *   HKM -> Let's change that!
   */
  if (mp->RelLiqPermModel != CONSTANT) {
    load_liq_perm(porosity, cap_pres, saturation, d_cap_pres);
  }
  if (mp->RelGasPermModel != CONSTANT) {
    load_gas_perm(porosity, cap_pres, saturation, d_cap_pres);
  }

  /*
   * Gas-phase Diffusivity model for liquid-phase solvent
   */
  if (mp->PorousDiffusivityModel[i_pl] != CONSTANT) {
    load_gas_diff(porosity, cap_pres, saturation, d_cap_pres, i_pl);
  }

  /*
   * Get the darcy fluxes for liquid and gas phase transport
   * in a porous phase
   */

  /*if (pd->e[pg->imtrx][R_POR_ENERGY]) {*/
  load_MandE_flux(porosity, cap_pres, saturation, d_cap_pres);
  /*}
 else {
  load_mass_flux(porosity, cap_pres, saturation, d_cap_pres);
 }*/

  /*  printf("\t rel_mass_flux[i_pl]=%g \n",pmv->rel_mass_flux[i_pl][0]);
  printf("\t rel_mass_flux[i_pl]=%g \n",pmv->rel_mass_flux[i_pl][1]);
  printf("\t rel_mass_flux[i_pl]=%g \n",pmv->rel_mass_flux[i_pl][2]);
  for (j=0; j<8; j++){
  printf("\t j=%d d_rel_mass_flux_dpmv[i_pl][0][i_pl]=%g
  \n",j,pmv->d_rel_mass_flux_dpmv[i_pl][0][i_pl][j]); printf("\t j=%d
  d_rel_mass_flux_dpmv[i_pl][1][i_pl]=%g \n",j,pmv->d_rel_mass_flux_dpmv[i_pl][1][i_pl][j]);
  printf("\t j=%d d_rel_mass_flux_dpmv[i_pl][2][i_pl]=%g
  \n",j,pmv->d_rel_mass_flux_dpmv[i_pl][2][i_pl][j]);
  }*/

  /*
   * Load the supg terms if requested.
   */
  if (mp->Porous_wt_funcModel == SUPG) {
    h_elem_siz(Stab->hsquared, Stab->hhv, Stab->dhv_dxnode, pd->e[pg->imtrx][R_MESH1]);
    Stab->Grid_Peclet_Number[POR_LIQ_PRES] = get_supg_terms_porous(Stab->hsquared, Stab->hhv);
  }
  is_initialized = TRUE;
  return 0;
} /*    end load_porous_properties  */
/**************************************************************************/
/**************************************************************************/
/**************************************************************************/
/*ARGSUSED*/
static void load_porous_properties_nodes(int lnn)

/**********************************************************************
 *
 * load_porous_properties_nodes()
 *
 *  This is called at a local node in order to calculate porous
 *  properties at a node.
 *
 * Input
 * -------
 *  It may be assumed that the basis functions are currently evaluated
 *  prior to this routine being called. However, the derivatives
 *  of the basis functions may not be calculated, before this routine
 *  is called.
 *
 *  lnn = Local node number
 *
 * Output
 * --------
 *  mp->porosity and derivatives wrt dependent variables
 *  mp->saturation and derivatives wrt dependent variables
 *  mp->cap_pres
 *  mp->gas
 *
 *  pmv->bulk_density[][]
 *
 **********************************************************************/
{
  double p_gas_star;
  const int i_pl = 0, i_pg = 1, i_pore = 2, i_pe = 3;
  double d_cap_pres[4];
  int w;

  /*
   *  Calculate porosity
   *
   */

  if (mp->PorosityModel == CONSTANT) {
    if (pd->TimeIntegration == TRANSIENT) {
      mp_old->porosity = mp->porosity;
      if (pd->v[pg->imtrx][POR_POROSITY])
        mp_old->d_porosity[POR_POROSITY] = 0.0;
      else
        mp_old->d_porosity[POR_POROSITY] = 1.0;
    }
  } else if (mp->PorosityModel == EXTERNAL_FIELD) {
    mp->porosity = fv->porosity;
    if (pd->TimeIntegration == TRANSIENT) {
      mp_old->porosity = mp->porosity;
    }
  } else if (mp->PorosityModel == DEFORM &&
             (cr->MeshMotion == LAGRANGIAN || cr->MeshMotion == DYNAMIC_LAGRANGIAN)) {
    mp->porosity = fv->porosity;
    mp->d_porosity[POR_POROSITY] = 1.0;
    if (pd->TimeIntegration == TRANSIENT) {
      mp_old->porosity = fv_old->porosity;
      mp_old->d_porosity[POR_POROSITY] = 1.0;
    }
  } else {
    GOMA_EH(GOMA_ERROR, "Unknown model for porosity or incorrect MeshMotion option");
  }

  if (mp->PorousMediaType == POROUS_SATURATED) {
    mp->saturation = 1.0;
    pmv->cap_pres = 0.0;
    pmv->liq_Xvol_solvents[i_pl] = 1.0;
    pmv->bulk_density[i_pl] = mp->density * mp->porosity;
  } else {
    /* MEDIUM IS PARTIALLY-SATURATED
     *  need to calculate the necessary properties for unsaturated
     *  or two-phase calculations in porous medium:
     *
     *  UNSATURATED: medium is filled with liquid and air, but air
     *               is assumed to be at atmospheric pressure everywhere
     *               so we just solve for capillary pressure
     *               which is atmospheric pressure minus the liquid pressure
     *  TWO_PHASE:   medium is filled with liquid and air, and we solve
     *               for both liquid pressure and air pressure
     *
     *
     *  Calculate capillary pressure and it's sensitivity to liquid
     *  and gas pressures
     *
     */
    if (mp->PorousMediaType == POROUS_UNSATURATED) {
      p_gas_star = mp->u_porous_gas_constants[3];
      pmv->cap_pres = p_gas_star - fv->p_liq;
      d_cap_pres[i_pl] = -1.;
      d_cap_pres[i_pg] = 1.;
      if (pmv->cap_pres < 0.0) {
        pmv->cap_pres = 0.0;
        d_cap_pres[i_pl] = 0.0;
      }
      if (pd->TimeIntegration == TRANSIENT) {
        pmv_old->cap_pres = p_gas_star - fv_old->p_liq;
        if (pmv_old->cap_pres < 0.0) {
          pmv_old->cap_pres = 0.0;
        }
      }
    } else if (mp->PorousMediaType == POROUS_TWO_PHASE) {
      p_gas_star = fv->p_gas;
      pmv->cap_pres = p_gas_star - fv->p_liq;
      d_cap_pres[i_pl] = -1.;
      d_cap_pres[i_pg] = 1.;
      d_cap_pres[i_pe] = 0.;

      if (pd->TimeIntegration == TRANSIENT) {
        pmv_old->cap_pres = fv_old->p_gas - fv_old->p_liq;
      }
    }

    /* calculate liquid volume fraction of all species - currently pure solvent */
    /* initialize the liquid volume fraction derivatives */
    /* gas is insoluble in liquid */

    if (mp->PorousMediaType == POROUS_TWO_PHASE) {
      pmv->liq_Xvol_solvents[i_pg] = 0.0;
      for (w = POR_LIQ_PRES; w < POR_LIQ_PRES + MAX_PMV; w++)
        pmv->d_liq_Xvol_solvents[i_pg][w] = 0.0;
      if (pd->TimeIntegration == TRANSIENT) {
        pmv_old->liq_Xvol_solvents[i_pg] = 0.0;
      }
    }

    /* solid is insoluble in liquid */
    if (pd->e[pg->imtrx][R_POR_POROSITY]) {
      pmv->liq_Xvol_solvents[i_pore] = 0.0;
      for (w = POR_LIQ_PRES; w < POR_LIQ_PRES + MAX_PMV; w++)
        pmv->d_liq_Xvol_solvents[i_pore][w] = 0.0;
      if (pd->TimeIntegration == TRANSIENT) {
        pmv_old->liq_Xvol_solvents[i_pore] = 0.0;
      }
    }

    /*
     * What's left over is solvent.
     * Currently the volume fraction of solvent in the solvent
     * phase is equal to 1.0. Later, when we make the solvent
     * phase multicomponent, this might change.
     */
    pmv->liq_Xvol_solvents[i_pl] = 1.0;
    for (w = POR_LIQ_PRES; w < POR_LIQ_PRES + MAX_PMV; w++)
      pmv->d_liq_Xvol_solvents[i_pl][w] = 0.0;
    if (pd->TimeIntegration == TRANSIENT) {
      pmv_old->liq_Xvol_solvents[i_pl] = 1.0;
    }

    /* Calculate saturation and it's sensitivity to liquid presssure, gas
     * pressure, and porosity.
     */
    MMH_ip = lnn; /*This is necessary for hysteresis, and we are arbitrarily
                    connecting the local node number to the gauss point number.
                    Gauss-points better be equal to nodes per element, or there
                    are no guarantees what will happen... */

    load_saturation(mp->porosity, pmv->cap_pres, d_cap_pres);
  }
  /*
   *   Calculate bulk mass concentrations of all species.
   *    - The bulk concentration includes the mass of the species
   *      in both gas and liquid phases
   */
  load_bulk_density(mp->porosity, pmv->cap_pres, mp->saturation, d_cap_pres);
}
/**************************************************************************/
/**************************************************************************/
/**************************************************************************/

double get_supg_terms_porous(double hsquare[DIM], double hh[DIM][DIM])

/*********************************************************************
 *
 * get_supg_terms_porous():
 *
 *    This subroutine calculates components of the supg expressions
 * This includes the local peclet number.
 *
 *  Input
 * -------
 *  Element size info along local element coodinate directions [lecd]:
 *
 *   hsquare[lecd] = Length of each local element coordinate
 *   hh[lecd][dim] = length of each local element coordinate
 *                   along each of the physical dimension axes.
 *
 *  The following constituitive relations are assumed to have
 *  been calculated at the gauss point before calling this routine:
 *
 *   mp->density
 *   mp->permeability
 *   mp->d_rel_liq_perm[POR_LIQ_PRES]
 *   mp->rel_liq_perm
 *   mp->viscosity
 *
 *   fv->grad_p_liq[a]
 *
 *  Output
 * -------
 *  This routine fills in the arrays in pmv pertaining to the
 *  SUPG stabilization.
 *
 *  pmv->U_supg[]       = The effective "velocities" to use in the
 *                        supg expressions.
 *  pmv->U_supg_squared = U_supg dot U_supg
 *  pmv->hveloc_lcd[lcd]= U_supg[] dot hh[][]
 *  Stab->h_veloc_elem  = ||  pmv->hveloc_lcd[lcd] ||
 *  pmv->zeta           = coth(peclet) - 1.0 / peclet
 *  pmv->k_art_diff     = h_veloc_elem * pmv->zeta / 2.0
 *  pmv->U_supg_hnorm[] = u_supg[a] * pmv->k_art_diff / u_supg_squared
 *                        (has units of length)
 *
 *  Return
 * ---------
 *  This function returns the grid peclet number
 *
 *********************************************************************/
{
  int a;
  int dim = pd->Num_Dim;
  double lambda_liq, tmp, u_supg_squared = 0.0, peclet;
  double h_veloc_elem = 0.0;
  double *h_veloc_lcd = pmv->h_veloc_lcd;
  double *u_supg = pmv->U_supg;
  double *u_supg_hnorm = pmv->U_supg_hnorm;

  /* Note:  we are going to take the permeability used here for the supg term as the
   * [1][1] component of the perm_tensor, if a perm_tensor model is used.   This
   * means ORTHOTROPIC, KC_TENSOR, or SM_TENSOR models.
   */

  if (mp->PermeabilityModel == ORTHOTROPIC || mp->PermeabilityModel == SM_TENSOR ||
      mp->PermeabilityModel == KC_TENSOR) {
    mp->permeability = mp->perm_tensor[0][0];
  }
  tmp = mp->density * mp->permeability * mp->d_rel_liq_perm[POR_LIQ_PRES];
  lambda_liq = mp->density * mp->permeability * mp->rel_liq_perm;
  for (a = 0; a < dim; a++) {
    u_supg[a] = -tmp * (fv->grad_p_liq[a]);
    u_supg_squared += u_supg[a] * u_supg[a];
  }
  pmv->U_supg_squared = u_supg_squared;

  /*
   * Calculate the h_veloc_elem = ||h dot u ||
   *
   *   In multiple dimensions, we use:
   *
   *          h_veloc_elem = sqrt(sum_lcd((h_lcd * u_lcd)**2)
   *  where
   *      lcd = local coordinate direction, i.e., kqsi and nu.
   */
  switch (dim) {
  case 2:
    h_veloc_lcd[0] = u_supg[0] * hh[0][0] + u_supg[1] * hh[0][1];
    h_veloc_lcd[1] = u_supg[0] * hh[1][0] + u_supg[1] * hh[1][1];
    h_veloc_elem = (h_veloc_lcd[0] * h_veloc_lcd[0] + h_veloc_lcd[1] * h_veloc_lcd[1]);
    break;
  case 3:
    h_veloc_lcd[0] = (u_supg[0] * hh[0][0] + u_supg[1] * hh[0][1] + u_supg[2] * hh[0][2]);
    h_veloc_lcd[1] = (u_supg[0] * hh[1][0] + u_supg[1] * hh[1][1] + u_supg[2] * hh[1][2]);
    h_veloc_lcd[2] = (u_supg[0] * hh[2][0] + u_supg[1] * hh[2][1] + u_supg[2] * hh[2][2]);
    h_veloc_elem = (h_veloc_lcd[0] * h_veloc_lcd[0] + h_veloc_lcd[1] * h_veloc_lcd[1] +
                    h_veloc_lcd[2] * h_veloc_lcd[2]);
    break;
  }
  h_veloc_elem = sqrt(h_veloc_elem);
  Stab->h_veloc_elem = h_veloc_elem;
  if (lambda_liq > 0.0) {
    peclet = h_veloc_elem / (lambda_liq * 2.0);
  } else {
    peclet = 1.0E10;
  }
  if (peclet == 0.0) {
    pmv->zeta = 0.0;
    pmv->k_art_diff = 0.0;
    for (a = 0; a < dim; a++) {
      u_supg_hnorm[a] = 0.0;
    }
  } else {
    pmv->zeta = 1.0 / tanh(peclet) - 1.0 / peclet;
    pmv->k_art_diff = h_veloc_elem * pmv->zeta * 0.5;
    if (u_supg_squared > 0.0) {
      for (a = 0; a < dim; a++) {
        u_supg_hnorm[a] = u_supg[a] * pmv->k_art_diff / u_supg_squared;
      }
    }
  }
  return peclet;
}
/**************************************************************************/
/**************************************************************************/
/**************************************************************************/

void load_nodal_porous_properties(double tt, double dt)

/*********************************************************************
 *
 * load_nodal_porous_properties():
 *
 *   This routine calculates the mass lumped quantities at node
 *   points
 *
 *  Input
 * -------
 *
 *  The following constituitive relations are assumed to have
 *  been calculated at the gauss point before calling this routine:
 *
 *   mp->density
 *   mp->permeability
 *   mp->d_rel_liq_perm[POR_LIQ_PRES]
 *   mp->rel_liq_perm
 *   mp->viscosity
 *
 *   fv->grad_p_liq[a]
 *
 *  Output
 * -------
 *
 *
 *
 *********************************************************************/
{
  int eqn, idof, lnn, i_lvdesc, lvd, err, w, w1;
  double p_liq, p_liq_old, xi[3];
  double p_gas, p_gas_old, p_porosity, p_porosity_old, p_T, p_T_old;
  int *lvdesc_to_lnn, *lvdesc_to_idof;
  const int i_pl = 0, i_pg = 1, i_pe = 3;

  eqn = POR_LIQ_PRES;
  i_lvdesc = ei[pg->imtrx]->Lvdesc_First_Var_Type[eqn];
  lvdesc_to_lnn = ei[pg->imtrx]->Lvdesc_to_Lnn[i_lvdesc];
  lvdesc_to_idof = ei[pg->imtrx]->Lvdesc_to_lvdof[i_lvdesc];
  for (lvd = 0; lvd < ei[pg->imtrx]->Lvdesc_Numdof[i_lvdesc]; lvd++) {
    idof = lvdesc_to_idof[lvd];
    p_liq = *(esp->p_liq[idof]);
    p_liq_old = *(esp_old->p_liq[idof]);
    fv->p_liq = p_liq;
    fv_old->p_liq = p_liq_old;
    if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
      p_gas = *(esp->p_gas[idof]);
      p_gas_old = *(esp_old->p_gas[idof]);
      fv->p_gas = p_gas;
      fv_old->p_gas = p_gas_old;
    }
    if (pd->e[pg->imtrx][R_POR_POROSITY]) {
      p_porosity = *(esp->porosity[idof]);
      p_porosity_old = *(esp_old->porosity[idof]);
      fv->porosity = p_porosity;
      fv_old->porosity = p_porosity_old;
    }
    if (pd->e[pg->imtrx][R_POR_ENERGY]) {
      p_T = *(esp->T[idof]);
      p_T_old = *(esp_old->T[idof]);
      fv->T = p_T;
      fv_old->T = p_T_old;
    }

    /*
     * Translate the degree of freedom into the local node number
     */
    lnn = lvdesc_to_lnn[lvd];

    /*
     * Find the correct local element coordinates, xi[], at the
     * the current local node number, lnn.
     */
    find_nodal_stu(lnn, ei[pg->imtrx]->ielem_type, xi, xi + 1, xi + 2);

    /*
     * Load up basis function information for each basis function
     * type needed by the current material at the current
     * location, xi[], in the element. This is done in terms
     * of local element coordinates.
     */
    err = load_basis_functions(xi, bfd);
    err = beer_belly();
    err = load_fv();

    /*
     * Note: we do not need derivatives of the field variables to
     *       calculate the capacitance terms in this routine.
     *       Therefore, we will skip these calculations. Also, we
     *       do not need
     */
    if (pd->e[pg->imtrx][R_MESH1]) {
      err = load_bf_grad();
      err = load_bf_mesh_derivs();
      GOMA_EH(err, "load_bf_mesh_derivs");
    }

    /*
     * Call constitutive relations with gauss point defined at the node
     */
    load_porous_properties_nodes(lnn);

    /*
     * Calculate Inventory_Solvent_dot[MDE][MAX_PMV]
     * and d_Inventory_solvent_dot_dpmv[MDE][MAX_PMV][MAX_PMV];
     */

    /*
     * Store calculated quantities in the pmv_ml structure
     */
    for (w = 0; w < MAX_PMV; w++) {
      pmv_ml->Bulk_Density[idof][w] = pmv->bulk_density[w];
      pmv_ml->Bulk_Density_old[idof][w] = pmv_old->bulk_density[w];
      pmv_ml->d_Bulk_Density[idof][w][POR_LIQ_PRES] = pmv->d_bulk_density[w][POR_LIQ_PRES];
      pmv_ml->d_Bulk_Density[idof][w][POR_GAS_PRES] = pmv->d_bulk_density[w][POR_GAS_PRES];
      pmv_ml->d_Bulk_Density[idof][w][POR_POROSITY] = pmv->d_bulk_density[w][POR_POROSITY];
      pmv_ml->d_Bulk_Density[idof][w][POR_TEMP] = pmv->d_bulk_density[w][POR_TEMP];

      pmv_ml->d_Bulk_Density_old[idof][w][POR_LIQ_PRES] = pmv_old->d_bulk_density[w][POR_LIQ_PRES];
      pmv_ml->d_Bulk_Density_old[idof][w][POR_GAS_PRES] = pmv_old->d_bulk_density[w][POR_GAS_PRES];
      pmv_ml->d_Bulk_Density_old[idof][w][POR_POROSITY] = pmv_old->d_bulk_density[w][POR_POROSITY];
      pmv_ml->d_Bulk_Density_old[idof][w][POR_TEMP] = pmv_old->d_bulk_density[w][POR_TEMP];
    }

    pmv_ml->Inventory_Solvent_old[idof][i_pl] = pmv_old->bulk_density[i_pl];
    pmv_ml->Inventory_Solvent_old[idof][i_pg] = pmv_old->bulk_density[i_pg];
    pmv_ml->Inventory_Solvent_old[idof][i_pe] = pmv_old->bulk_density[i_pe];

    pmv_ml->Inventory_Solvent[idof][i_pl] = pmv->bulk_density[i_pl];
    pmv_ml->Inventory_Solvent[idof][i_pg] = pmv->bulk_density[i_pg];
    pmv_ml->Inventory_Solvent[idof][i_pe] = pmv->bulk_density[i_pe];
    /*
     * Calculate Inventory_Solvent_dot[MDE][MAX_PMV]
     *
     */
    if (tt > 0.0) {
      pmv_ml->Inventory_Solvent_dot_old[idof][i_pl] =
          fv_dot_old->p_liq * pmv_ml->d_Bulk_Density_old[idof][i_pl][POR_LIQ_PRES] +
          fv_dot_old->p_gas * pmv_ml->d_Bulk_Density_old[idof][i_pl][POR_GAS_PRES] +
          fv_dot_old->porosity * pmv_ml->d_Bulk_Density_old[idof][i_pl][POR_POROSITY] +
          fv_dot_old->T * pmv_ml->d_Bulk_Density_old[idof][i_pl][POR_TEMP];

      pmv_ml->Inventory_Solvent_dot[idof][i_pl] =
          (1.0 + 2.0 * tt) *
              (pmv_ml->Inventory_Solvent[idof][i_pl] - pmv_ml->Inventory_Solvent_old[idof][i_pl]) /
              dt -
          2.0 * tt * pmv_ml->Inventory_Solvent_dot_old[idof][i_pl];

      pmv_ml->Inventory_Solvent_dot_old[idof][i_pg] =
          fv_dot_old->p_liq * pmv_ml->d_Bulk_Density_old[idof][i_pg][POR_LIQ_PRES] +
          fv_dot_old->p_gas * pmv_ml->d_Bulk_Density_old[idof][i_pg][POR_GAS_PRES] +
          fv_dot_old->porosity * pmv_ml->d_Bulk_Density_old[idof][i_pg][POR_POROSITY] +
          fv_dot_old->T * pmv_ml->d_Bulk_Density_old[idof][i_pg][POR_TEMP];

      pmv_ml->Inventory_Solvent_dot[idof][i_pg] =
          (1.0 + 2.0 * tt) *
              (pmv_ml->Inventory_Solvent[idof][i_pg] - pmv_ml->Inventory_Solvent_old[idof][i_pg]) /
              dt -
          2.0 * tt * pmv_ml->Inventory_Solvent_dot_old[idof][i_pg];

      pmv_ml->Inventory_Solvent_dot_old[idof][i_pe] =
          fv_dot_old->p_liq * pmv_ml->d_Bulk_Density_old[idof][i_pe][POR_LIQ_PRES] +
          fv_dot_old->p_gas * pmv_ml->d_Bulk_Density_old[idof][i_pe][POR_GAS_PRES] +
          fv_dot_old->porosity * pmv_ml->d_Bulk_Density_old[idof][i_pe][POR_POROSITY] +
          fv_dot_old->T * pmv_ml->d_Bulk_Density_old[idof][i_pe][POR_TEMP];

      pmv_ml->Inventory_Solvent_dot[idof][i_pe] =
          (1.0 + 2.0 * tt) *
              (pmv_ml->Inventory_Solvent[idof][i_pe] - pmv_ml->Inventory_Solvent_old[idof][i_pe]) /
              dt -
          2.0 * tt * pmv_ml->Inventory_Solvent_dot_old[idof][i_pe];

    } else {
      pmv_ml->Inventory_Solvent_dot[idof][i_pl] =
          (pmv_ml->Inventory_Solvent[idof][i_pl] - pmv_ml->Inventory_Solvent_old[idof][i_pl]) / dt;

      pmv_ml->Inventory_Solvent_dot[idof][i_pg] =
          (pmv_ml->Inventory_Solvent[idof][i_pg] - pmv_ml->Inventory_Solvent_old[idof][i_pg]) / dt;

      pmv_ml->Inventory_Solvent_dot[idof][i_pe] =
          (pmv_ml->Inventory_Solvent[idof][i_pe] - pmv_ml->Inventory_Solvent_old[idof][i_pe]) / dt;
    }

    /*
     * Store d_Inventory_solvent_dot_dpmv[MDE][MAX_PMV][MAX_PMV]
     * calculated in other places
     */
    if (af->Assemble_Jacobian) {
      for (w1 = 0; w1 < MAX_PMV; w1++) {
        pmv_ml->d_Inventory_Solvent_dot_dpmv[idof][i_pl][w1] =
            (1 + 2. * tt) * pmv_ml->d_Bulk_Density[idof][i_pl][POR_LIQ_PRES + w1] / dt;

        if (pd->e[pg->imtrx][R_POR_GAS_PRES])
          pmv_ml->d_Inventory_Solvent_dot_dpmv[idof][i_pg][w1] +=
              (1 + 2. * tt) * pmv_ml->d_Bulk_Density[idof][i_pg][POR_LIQ_PRES + w1] / dt;

        if (pd->e[pg->imtrx][R_POR_ENERGY])
          pmv_ml->d_Inventory_Solvent_dot_dpmv[idof][i_pe][w1] +=
              (1 + 2. * tt) * pmv_ml->d_Bulk_Density[idof][i_pe][POR_LIQ_PRES + w1] / dt;
      }
    }
  }
} /* END load_nodal_porous_properties() */
/**************************************************************************/
/**************************************************************************/
/**************************************************************************/

void load_nodal_shell_porous_properties(double tt, double dt, int eqn)

/*********************************************************************
 *
 * load_nodal_shell_porous_properties():
 *
 *   This routine calculates the mass lumped quantities at node
 *   points - This is shell version of load_nodal_porous_properties
 *
 *  Input
 * -------
 *
 *
 *  Output
 * -------
 *
 *  pmv_ml->Inventory_Solvent
 *  pmv_ml->Inventory_Solvent_old
 *  pmv_ml->Inventory_Solvent_dot
 *
 *
 *********************************************************************/
{
  int idof, lnn, i_lvdesc, lvd, err;
  int i_ext_field = 0;
  double xi[3];
  double phi, H, Patm, cap_pres, d_cap_pres[2];
  int *lvdesc_to_lnn, *lvdesc_to_idof;
  int i_pl = 0; /*Piggyback pmv_ml structure, i_pl = 0 --> shell_press_open */
                /*                            i_pl = 1 --> shell_press_open_2*/

  if (eqn == R_SHELL_SAT_OPEN_2)
    i_pl = 1;

  i_lvdesc = ei[pg->imtrx]->Lvdesc_First_Var_Type[eqn];
  lvdesc_to_lnn = ei[pg->imtrx]->Lvdesc_to_Lnn[i_lvdesc];
  lvdesc_to_idof = ei[pg->imtrx]->Lvdesc_to_lvdof[i_lvdesc];
  for (lvd = 0; lvd < ei[pg->imtrx]->Lvdesc_Numdof[i_lvdesc]; lvd++) {
    idof = lvdesc_to_idof[lvd];

    /*
     * Translate the degree of freedom into the local node number
     */
    lnn = lvdesc_to_lnn[lvd];

    /*
     * Find the correct local element coordinates, xi[], at the
     * the current local node number, lnn.
     */
    find_nodal_stu(lnn, ei[pg->imtrx]->ielem_type, xi, xi + 1, xi + 2);

    /*
     * Load up basis function information for each basis function
     * type needed by the current material at the current
     * location, xi[], in the element. This is done in terms
     * of local element coordinates.
     */
    err = load_basis_functions(xi, bfd);
    GOMA_EH(err, "load_basis_functions");
    err = beer_belly();
    GOMA_EH(err, "beer_belly");
    err = load_fv();
    GOMA_EH(err, "load_fv");

    /*
     * Evaluate constitutive relations with gauss point defined at the node
     */

    /* First, get porosity */
    if (mp->PorosityModel == CONSTANT) {
      if (pd->TimeIntegration == TRANSIENT) {
        mp_old->porosity = mp->porosity;
      }
    } else if (mp->PorosityModel == EXTERNAL_FIELD) {
      i_ext_field = mp->porosity_external_field_index;
      mp->porosity = fv->external_field[i_ext_field];
      if (pd->TimeIntegration == TRANSIENT) {
        mp_old->porosity = mp->porosity;
      }
    } else {
      GOMA_EH(GOMA_ERROR, "Only CONSTANT and EXTERNAL_FIELD porosity models are supported in shell "
                          "porous open equation");
    }

    phi = mp->porosity; /* Porosity */

    H = porous_shell_closed_height_model(); /* Pore height (vertical */
    Patm = mp->PorousShellPatm;             /* Gas pressure - always constant */
    cap_pres = Patm - fv->sh_p_open;
    if (eqn == R_SHELL_SAT_OPEN_2)
      cap_pres = Patm - fv->sh_p_open_2;
    if (pd->TimeIntegration == TRANSIENT) {
      pmv_old->cap_pres = Patm - fv_old->sh_p_open;
      if (eqn == R_SHELL_SAT_OPEN_2)
        pmv_old->cap_pres = Patm - fv_old->sh_p_open_2;
    }
    load_saturation(phi, cap_pres, d_cap_pres);

    /*
     * Calculate Inventory_Solvent_dot[MDE][MAX_PMV]
     * and d_Inventory_solvent_dot_dpmv[MDE][MAX_PMV][MAX_PMV];
     */

    /*
     * Store calculated quantities in the pmv_ml structure
     */

    pmv_ml->Inventory_Solvent_old[idof][i_pl] = H * mp->porosity * mp_old->saturation;
    pmv_ml->Inventory_Solvent[idof][i_pl] = H * mp->porosity * mp->saturation;

    /*
     * Calculate Inventory_Solvent_dot[MDE][MAX_PMV]
     *
     */
    if (tt > 0.0) {
      pmv_ml->Inventory_Solvent_dot_old[idof][i_pl] =
          fv_dot_old->sh_p_open * H * mp->porosity * mp_old->d_saturation[SHELL_PRESS_OPEN];
      if (eqn == R_SHELL_SAT_OPEN_2)
        pmv_ml->Inventory_Solvent_dot_old[idof][i_pl] = fv_dot_old->sh_p_open_2 * H *
                                                        mp_old->porosity *
                                                        mp_old->d_saturation[SHELL_PRESS_OPEN_2];

      pmv_ml->Inventory_Solvent_dot[idof][i_pl] =
          (1.0 + 2.0 * tt) *
              (pmv_ml->Inventory_Solvent[idof][i_pl] - pmv_ml->Inventory_Solvent_old[idof][i_pl]) /
              dt -
          2.0 * tt * pmv_ml->Inventory_Solvent_dot_old[idof][i_pl];

    } else {
      pmv_ml->Inventory_Solvent_dot[idof][i_pl] =
          (pmv_ml->Inventory_Solvent[idof][i_pl] - pmv_ml->Inventory_Solvent_old[idof][i_pl]) / dt;
    }

    /*
     * Store d_Inventory_solvent_dot_dpmv[MDE][MAX_PMV][MAX_PMV]
     * calculated in other places
     *
     * Right now it only has R_SHELL_SAT_OPEN as eqn and SHELL_PRESS_OPEN as var
     * Or R_SHELL_SAT_OPEN_2 as eqn and SHELL_PRESS_OPEN_2 as var
     */
    if (af->Assemble_Jacobian) {
      pmv_ml->d_Inventory_Solvent_dot_dpmv[idof][i_pl][i_pl] =
          (1 + 2. * tt) * H * mp->porosity * mp->d_saturation[SHELL_PRESS_OPEN] / dt;

      if (eqn == R_SHELL_SAT_OPEN_2)
        pmv_ml->d_Inventory_Solvent_dot_dpmv[idof][i_pl][i_pl] =
            (1 + 2. * tt) * H * mp->porosity * mp->d_saturation[SHELL_PRESS_OPEN_2] / dt;
    }
  }
} /* END load_nodal_shell_porous_properties() */
/**************************************************************************/
/**************************************************************************/
/**************************************************************************/

int get_porous_part_sat_terms(struct Porous_Media_Terms *pmt,
                              double tt, /* time integration scheme param */
                              double dt) /* current time step size */

/*********************************************************************
 *  get_porous_part_sat_terms():
 *
 *
 *   Calculate the capacity, flux, convection, and source terms
 *   for a porous medium
 *   This routine calculates all of the terms that will go into the
 *   fill routines at the volumetric gauss point level.
 *
 *   Input
 * ---------
 *   pmv  - All terms in this structure are considered to have already
 *          been evaluated at the guass point.
 *
 *   Output
 * --------
 *   pmt  - This routine is responsible for filling out all terms
 *          in this structure.
 **********************************************************************/
{
  int w, w1, i, j, a, b, err, var, eqn;
  const int i_pl = 0, i_pg = 1, i_pore = 2, i_pe = 3;
  double value, delta_val, cap_pres, p_gas_star = 0.0;
  int dim = pd->Num_Dim;
  double *grad_phi_i, dlamdPl, tmp = 0.0, *grad_phi_j;
  double phi_j, d_cap_pres[3];

  double MassSource = 0.0;
  double d_MassSource[MAX_VARIABLE_TYPES + MAX_CONC][MDE];

  /* Storage for call to get_convection_velocity()
   *  -> lagrangian motion of solid
   */
  double vconv[MAX_PDIM];     /*Calculated convection velocity */
  double vconv_old[MAX_PDIM]; /*Calculated convection velocity at previous time*/
  CONVECTION_VELOCITY_DEPENDENCE_STRUCT d_vconv_struct;
  CONVECTION_VELOCITY_DEPENDENCE_STRUCT *d_vconv = &d_vconv_struct;

  /*
   * Intialize all terms in the pmt structure to zero or NULL.
   *  HKM -> note in the future, we may elect to zero only
   *         terms that are actually used.
   */

  if (1) {
    zero_structure(pmt, sizeof(struct Porous_Media_Terms), 1);
  }
  /*
   * Get the convection velocity for Lagrangian motion
   * (it's different for arbitrary and lagrangian meshes)
   */
  if ((pd->MeshMotion == LAGRANGIAN || cr->MeshMotion == DYNAMIC_LAGRANGIAN)) {
    err = belly_flop(elc->lame_mu);
    GOMA_EH(err, "error in belly flop");
    if (neg_elem_volume)
      return (err);
  }

  /*
   * Get the stress-free-state velocity piece for arbitrary and
   * Lagrangian motion.
   */
  err = get_convection_velocity(vconv, vconv_old, d_vconv, dt, tt);
  GOMA_EH(err, "Error in calculating effective convection velocity");

  /* Get mass source term if PORE_SINK_MASS is on */

  if (pd->e[pg->imtrx][R_POR_SINK_MASS]) {
    MassSource = por_mass_source_model(d_MassSource);
  }

  /***********************************************************************
   * CAPACITY TERM -
   *     mass density (gm/cm**3) of liquid and gas solvents
   *     in all phases present at the gauss point. also, energy (J/cm**3)
   *     in all phases present, if the energy eqn is to be solved.
   ***********************************************************************/

  pmt->Inventory_solvent[i_pl] = pmv->bulk_density[i_pl]; /*rho*S*phi for PART_SAT */
  pmt->Inventory_solvent[i_pg] = pmv->bulk_density[i_pg];
  pmt->Inventory_solvent[i_pe] = pmv->bulk_density[i_pe];

  if (pd->TimeIntegration != STEADY) {
    /*
     * Liquid solvent component first, e.g. water
     */
    pmt->Inventory_solvent_old[i_pl] = pmv_old->bulk_density[i_pl];
    /*
     *  convert Inventory_solvent_dot_old by multiplying by
     *  capacitance matrix with chain rule
     */

    pmt->Inventory_solvent_dot_old[i_pl] =
        fv_dot_old->p_liq * pmv_old->d_bulk_density[i_pl][POR_LIQ_PRES] +
        fv_dot_old->p_gas * pmv_old->d_bulk_density[i_pl][POR_GAS_PRES] +
        fv_dot_old->porosity * pmv_old->d_bulk_density[i_pl][POR_POROSITY] +
        fv_dot_old->T * pmv_old->d_bulk_density[i_pl][POR_TEMP];

    pmt->Inventory_solvent_dot[i_pl] =
        (1.0 + 2.0 * tt) * (pmt->Inventory_solvent[i_pl] - pmt->Inventory_solvent_old[i_pl]) / dt -
        2.0 * tt * pmt->Inventory_solvent_dot_old[i_pl];

    /*
     * Gas solvent component, e.g. air
     */

    pmt->Inventory_solvent_old[i_pg] = pmv_old->bulk_density[i_pg];

    /* convert Inventory_solvent_dot_old by multiplying by
     * capacitance matrix with chain rule
     */
    if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
      pmt->Inventory_solvent_dot_old[i_pg] =
          fv_dot_old->p_liq * pmv_old->d_bulk_density[i_pg][POR_LIQ_PRES] +
          fv_dot_old->p_gas * pmv_old->d_bulk_density[i_pg][POR_GAS_PRES] +
          fv_dot_old->porosity * pmv_old->d_bulk_density[i_pg][POR_POROSITY] +
          fv_dot_old->T * pmv_old->d_bulk_density[i_pg][POR_TEMP];

      /*
       * Formulate the current derivative of the total amount of water density
       * at the current gauss point
       */
      pmt->Inventory_solvent_dot[i_pg] =
          (1 + 2. * tt) * (pmt->Inventory_solvent[i_pg] - pmt->Inventory_solvent_old[i_pg]) / dt -
          2. * tt * pmt->Inventory_solvent_dot_old[i_pg];
    }

    /*
     * Energy
     */
    if (pd->e[pg->imtrx][R_POR_ENERGY]) {
      pmt->Inventory_solvent_old[i_pe] = pmv_old->bulk_density[i_pe];
      /*
       *  convert Inventory_solvent_dot_old by multiplying by
       *  capacitance matrix with chain rule
       */
      pmt->Inventory_solvent_dot_old[i_pe] =
          fv_dot_old->p_liq * pmv_old->d_bulk_density[i_pe][POR_LIQ_PRES] +
          fv_dot_old->p_gas * pmv_old->d_bulk_density[i_pe][POR_GAS_PRES] +
          fv_dot_old->porosity * pmv_old->d_bulk_density[i_pe][POR_POROSITY] +
          fv_dot_old->T * pmv_old->d_bulk_density[i_pe][POR_TEMP];

      pmt->Inventory_solvent_dot[i_pe] =
          (1.0 + 2.0 * tt) * (pmt->Inventory_solvent[i_pe] - pmt->Inventory_solvent_old[i_pe]) /
              dt -
          2.0 * tt * pmt->Inventory_solvent_dot_old[i_pe];
    }

  } else {
    pmt->Inventory_solvent_dot[i_pl] = 0.0;
    pmt->Inventory_solvent_dot[i_pg] = 0.0;
    pmt->Inventory_solvent_dot[i_pe] = 0.0;
  }

  /*
   * No capacitance term for solid phase.
   */
  pmt->Inventory_solvent[i_pore] = 0.0;
  pmt->Inventory_solvent_dot[i_pore] = 0.0;

  /*
   * DIFFUSIVE FLUX TERM - this is the term that is integrated by parts in the
   *   species conservation residual equation.  In a porous medium it contains
   *   five parts of the species flux:
   *    1) convection of the bulk concentration of each species with a moving solid
   *    2) convection of each species in gas phase relative to solid (from Darcy)
   *    3) convection of each species in liquid phase relative to solid (from Darcy)
   *    4) diffusion in gas phase (relative to gas Darcy Flux)
   *    5) diffusion in liquid phase (relative to liquid Darcy Flux)
   *         - currently this is neglected
   */

  if (cr->PorousFluxModel == DARCY_FICKIAN || cr->PorousFluxModel == POWERLAW_DARCY_FICKIAN) {
    /*
     *  Here, vconv is the velocity of the solid phase
     *  rel_mass_flux[][] contains both the Darcy flux and the diffusive
     *  flux components.
     */
    for (a = 0; a < VIM; a++) {
      pmt->diff_flux[i_pl][a] = (pmv->rel_mass_flux[i_pl][a] + vconv[a] * pmv->bulk_density[i_pl]);
      if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
        pmt->diff_flux[i_pg][a] =
            (pmv->rel_mass_flux[i_pg][a] + vconv[a] * pmv->bulk_density[i_pg]);
      }

      if (pd->e[pg->imtrx][R_POR_ENERGY])
        pmt->diff_flux[i_pe][a] =
            (pmv->rel_mass_flux[i_pe][a] + vconv[a] * pmv->bulk_density[i_pe]);
    }
  } else {
    GOMA_EH(GOMA_ERROR, "Unimplemented flux constitutive relation in porous media. "
                        "Check Porous Diffusion Constitutive Equation Card");
  }

  /*
   * CONVECTIVE FLUX TERM
   *
   *  This term stems from the fact that Goma actually tracks the material derivative
   *  of the inventory, instead of the ordinary time derivative.
   */

  /* memset(pmt->conv_flux, 0, sizeof(double)*MAX_PMV*DIM); */
  if (pd->v[pg->imtrx][MESH_DISPLACEMENT1]) {
    for (a = 0; a < VIM; a++) {
      pmt->conv_flux[i_pl][a] = pmv->bulk_density[i_pl] * fv->grad_d_dot[a][a];

      if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
        pmt->conv_flux[i_pg][a] = pmv->bulk_density[i_pg] * fv->grad_d_dot[a][a];
      }
      if (pd->e[pg->imtrx][R_POR_ENERGY])
        pmt->conv_flux[i_pe][a] = pmv->bulk_density[i_pe] * fv->grad_d_dot[a][a];
    }
  }

  /*
   *  Need to calculate the dilation of the medium and relate it
   *  to the porosity. This is really a Kludge for the porosity equation
   *  - replace first mass-flux piece with
   *    the relationship between porosity and deformation
   */
  if (pd->e[pg->imtrx][R_POR_POROSITY]) {
    if (pd->e[pg->imtrx][R_POR_SINK_MASS]) {
      pmt->conv_flux[i_pore][0] =
          (fv->volume_change -
           (mp->u_porous_sink_constants[3] / mp->u_porous_sink_constants[4] +
            mp->u_porous_sink_constants[2] / mp->u_porous_sink_constants[5] *
                (1 + fv->sink_mass * mp->u_porous_sink_constants[5] / mp->density)) /
               (1.0 - fv->porosity));

    } else {
      pmt->conv_flux[i_pore][0] =
          (fv->volume_change - (1.0 - mp->u_porosity[0]) / (1.0 - fv->porosity));
    }
  }

  /*
   * SOURCE TERM - leave reaction rates set to zero for now
   */

  pmt->MassSource[i_pl] = 0.0;
  pmt->MassSource[i_pore] = 0.0;
  pmt->MassSource[i_pe] = 0.0;

  if (pd->e[pg->imtrx][R_POR_SINK_MASS]) {
    pmt->MassSource[i_pl] = MassSource;
  }

  /*** Heat Source - constant only for now ****/

  if (mp->HeatSourceModel == USER) {
    GOMA_EH(GOMA_ERROR, "User heat source model not defined");
  } else if (mp->HeatSourceModel == CONSTANT) {
    pmt->MassSource[i_pe] = mp->heat_source;
  } else {
    GOMA_EH(GOMA_ERROR, "Unrecognized heat source model");
  }

  /*
   * SUPG TERMS
   *
   *   pi_supg[i] = contribution to the basis function from SUPG.
   *   conv_flux_supg[i_pl] = New term representing convective upwinding
   *                          contribution.
   */
  if (mp->Porous_wt_funcModel == SUPG) {
    eqn = POR_LIQ_PRES;
    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {

      grad_phi_i = bf[eqn]->grad_phi[i];

      pmt->pi_supg[i] = 0.0;

      for (a = 0; a < dim; a++) {
        pmt->pi_supg[i] += pmv->U_supg_hnorm[a] * grad_phi_i[a];
      }
    }
    dlamdPl = (mp->density * mp->permeability * mp->d_rel_liq_perm[POR_LIQ_PRES]);
    pmt->conv_flux_supg[i_pl] = 0.0;
    for (a = 0; a < dim; a++) {
      pmt->conv_flux_supg[i_pl] -= fv->grad_p_liq[a] * fv->grad_p_liq[a];
    }
    pmt->conv_flux_supg[i_pl] *= dlamdPl;
  }

  /********************************************************************************
   * NOW, CALCULATE SENSITIVITIES for the Jacobian, if needed
   ********************************************************************************/

  if (af->Assemble_Jacobian) {

    /*
     * sensitivity of CAPACITY TERM - concentration (volume fraction)
     */
    var = POR_LIQ_PRES;

    for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
      for (w = 0; w < MAX_PMV; w++) {
        for (w1 = 0; w1 < MAX_PMV; w1++) {
          pmt->d_Inventory_solvent_dot_dpmv[w][w1][j] = 0.;
        }
      }
    }

    if (pd->TimeIntegration != STEADY) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        for (w1 = 0; w1 < MAX_PMV; w1++) {
          pmt->d_Inventory_solvent_dot_dpmv[i_pl][w1][j] =
              (1 + 2. * tt) * pmv->d_bulk_density[i_pl][POR_LIQ_PRES + w1] * bf[var]->phi[j] / dt;

          if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
            pmt->d_Inventory_solvent_dot_dpmv[i_pg][w1][j] =
                (1 + 2. * tt) * pmv->d_bulk_density[i_pg][POR_LIQ_PRES + w1] * bf[var]->phi[j] / dt;

            if (pd->e[pg->imtrx][R_POR_ENERGY])
              pmt->d_Inventory_solvent_dot_dpmv[i_pe][w1][j] =
                  (1 + 2. * tt) * pmv->d_bulk_density[i_pe][POR_LIQ_PRES + w1] * bf[var]->phi[j] /
                  dt;
          }
        }
      }
    }

    /*
     * sensitivity of DIFFUSIVE FLUX TERM - concentration, temp, displacement
     */
    var = POR_LIQ_PRES;

    for (w = 0; w < MAX_PMV; w++) {
      for (a = 0; a < VIM; a++) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          for (w1 = 0; w1 < MAX_PMV; w1++) {
            pmt->d_diff_flux_dpmv[w][a][w1][j] = 0.;
          }
        }
      }
    }

    if (cr->PorousFluxModel == DARCY_FICKIAN || cr->PorousFluxModel == POWERLAW_DARCY_FICKIAN) {
      for (a = 0; a < VIM; a++) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          for (w1 = 0; w1 < MAX_PMV; w1++) {
            pmt->d_diff_flux_dpmv[i_pl][a][w1][j] = pmv->d_rel_mass_flux_dpmv[i_pl][a][w1][j];

            if (pd->e[pg->imtrx][R_POR_GAS_PRES])
              pmt->d_diff_flux_dpmv[i_pg][a][w1][j] = pmv->d_rel_mass_flux_dpmv[i_pg][a][w1][j];

            if (pd->e[pg->imtrx][R_POR_ENERGY])
              pmt->d_diff_flux_dpmv[i_pe][a][w1][j] = pmv->d_rel_mass_flux_dpmv[i_pe][a][w1][j];
          }

          if (pd->e[pg->imtrx][R_POR_SINK_MASS])
            pmt->d_diff_flux_dSM[i_pl][a][j] = pmv->d_rel_mass_flux_dSM[i_pl][a][j];
        }
      }
    } else {
      GOMA_EH(GOMA_ERROR, "Unimplemented flux constitutive relation in porous media. "
                          "Check Porous Diffusion Constitutive Equation Card");
    }

    for (b = 0; b < pd->Num_Dim; b++) {
      var = MESH_DISPLACEMENT1 + b;
      if (pd->v[pg->imtrx][var]) {

        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          for (a = 0; a < VIM; a++) {
            for (w = 0; w < MAX_PMV; w++) {
              pmt->d_diff_flux_dmesh[w][a][b][j] = 0.0;
            }
          }
          if (cr->PorousFluxModel == DARCY || cr->PorousFluxModel == DARCY_FICKIAN ||
              cr->PorousFluxModel == POWERLAW_DARCY_FICKIAN) {
            for (a = 0; a < VIM; a++) {
              pmt->d_diff_flux_dmesh[i_pl][a][b][j] +=
                  pmv->d_rel_mass_flux_dmesh[i_pl][a][b][j] +
                  d_vconv->X[a][b][j] * pmv->bulk_density[i_pl];

              if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
                pmt->d_diff_flux_dmesh[i_pg][a][b][j] +=
                    pmv->d_rel_mass_flux_dmesh[i_pg][a][b][j] +
                    d_vconv->X[a][b][j] * pmv->bulk_density[i_pg];
              }

              if (pd->e[pg->imtrx][R_POR_ENERGY])
                pmt->d_diff_flux_dmesh[i_pe][a][b][j] +=
                    pmv->d_rel_mass_flux_dmesh[i_pe][a][b][j] +
                    d_vconv->X[a][b][j] * pmv->bulk_density[i_pe];
            }
          }
        }
      }
    }

    var = TEMPERATURE;
    if (pd->v[pg->imtrx][var] && !pd->e[pg->imtrx][R_POR_ENERGY]) {
      if (cr->PorousFluxModel == DARCY_FICKIAN || cr->PorousFluxModel == POWERLAW_DARCY_FICKIAN) {
        for (a = 0; a < VIM; a++) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            pmt->d_diff_flux_dT[i_pl][a][j] = 0.;
            pmt->d_diff_flux_dT[i_pg][a][j] = 0.;
          }
        }
      }
      /* NOTE - the gas concentrations, etc should also be functions of temperature */
      for (a = 0; a < VIM; a++) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          pmt->d_diff_flux_dT[i_pl][a][j] = -pmv->d_rel_mass_flux_dT[i_pl][a][j];
          pmt->d_diff_flux_dT[i_pg][a][j] = -pmv->d_rel_mass_flux_dT[i_pg][a][j];
        }
      }
    }

    /*
     * CONVECTIVE FLUX TERM
     */

    if (pd->v[pg->imtrx][MESH_DISPLACEMENT1]) {
      for (a = 0; a < VIM; a++) {
        var = POR_LIQ_PRES;
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          pmt->d_conv_flux_dpmv[i_pl][a][i_pl][j] =
              fv->grad_d_dot[a][a] * pmv->d_bulk_density[i_pl][POR_LIQ_PRES] * bf[var]->phi[j];

          if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
            pmt->d_conv_flux_dpmv[i_pg][a][i_pl][j] =
                fv->grad_d_dot[a][a] * pmv->d_bulk_density[i_pg][POR_LIQ_PRES] * bf[var]->phi[j];

            if (pd->e[pg->imtrx][R_POR_ENERGY])
              pmt->d_conv_flux_dpmv[i_pe][a][i_pl][j] +=
                  fv->grad_d_dot[a][a] * pmv->d_bulk_density[i_pe][POR_LIQ_PRES] * bf[var]->phi[j];
          }
        }
        var = POR_GAS_PRES;
        if (pd->e[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            pmt->d_conv_flux_dpmv[i_pl][a][i_pg][j] =
                fv->grad_d_dot[a][a] * pmv->d_bulk_density[i_pl][POR_GAS_PRES] * bf[var]->phi[j];

            pmt->d_conv_flux_dpmv[i_pg][a][i_pg][j] =
                fv->grad_d_dot[a][a] * pmv->d_bulk_density[i_pg][POR_GAS_PRES] * bf[var]->phi[j];

            pmt->d_conv_flux_dpmv[i_pe][a][i_pg][j] =
                fv->grad_d_dot[a][a] * pmv->d_bulk_density[i_pe][POR_GAS_PRES] * bf[var]->phi[j];
          }
        }

        var = POR_POROSITY;
        if (pd->e[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            pmt->d_conv_flux_dpmv[i_pl][a][i_pore][j] =
                fv->grad_d_dot[a][a] * pmv->d_bulk_density[i_pl][POR_POROSITY] * bf[var]->phi[j];

            if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
              pmt->d_conv_flux_dpmv[i_pg][a][i_pore][j] =
                  fv->grad_d_dot[a][a] * pmv->d_bulk_density[i_pg][POR_POROSITY] * bf[var]->phi[j];
            }

            if (pd->e[pg->imtrx][R_POR_ENERGY])
              pmt->d_conv_flux_dpmv[i_pe][a][i_pore][j] =
                  fv->grad_d_dot[a][a] * pmv->d_bulk_density[i_pe][POR_POROSITY] * bf[var]->phi[j];
          }
        }

        var = POR_TEMP;
        if (pd->e[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            pmt->d_conv_flux_dpmv[i_pl][a][i_pe][j] =
                fv->grad_d_dot[a][a] * pmv->d_bulk_density[i_pl][POR_TEMP] * bf[var]->phi[j];

            pmt->d_conv_flux_dpmv[i_pg][a][i_pe][j] =
                fv->grad_d_dot[a][a] * pmv->d_bulk_density[i_pg][POR_TEMP] * bf[var]->phi[j];

            pmt->d_conv_flux_dpmv[i_pe][a][i_pe][j] =
                fv->grad_d_dot[a][a] * pmv->d_bulk_density[i_pe][POR_TEMP] * bf[var]->phi[j];
          }
        }
      }
    }
    for (b = 0; b < pd->Num_Dim; b++) {
      var = MESH_DISPLACEMENT1 + b;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          for (a = 0; a < VIM; a++) {
            for (w = 0; w < MAX_PMV; w++) {
              pmt->d_conv_flux_dmesh[w][a][b][j] = 0.0;
            }
          }
          for (a = 0; a < VIM; a++) {
            pmt->d_conv_flux_dmesh[i_pl][a][b][j] =
                pmv->bulk_density[i_pl] * fv->d_grad_d_dot_dmesh[a][a][b][j];

            if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
              pmt->d_conv_flux_dmesh[i_pg][a][b][j] =
                  pmv->bulk_density[i_pg] * fv->d_grad_d_dot_dmesh[a][a][b][j];
            }
          }
        }
      }
    }
    /* Recall from above that until we resolve this term, the only contribution
     * here goes to the solid-phase solvent balance
     */
    if (pd->e[pg->imtrx][R_POR_POROSITY]) {
      var = POR_POROSITY;
      if (pd->e[pg->imtrx][R_POR_SINK_MASS]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          pmt->d_conv_flux_dpmv[i_pore][0][i_pore][j] =
              -(mp->u_porous_sink_constants[3] / mp->u_porous_sink_constants[4] +
                mp->u_porous_sink_constants[2] / mp->u_porous_sink_constants[5] *
                    (1 + fv->sink_mass * mp->u_porous_sink_constants[5] / mp->density)) /
              (1. - fv->porosity) / (1. - fv->porosity) * bf[var]->phi[j];
        }
      } else {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          pmt->d_conv_flux_dpmv[i_pore][0][i_pore][j] = -(1. - mp->u_porosity[0]) /
                                                        (1. - fv->porosity) / (1. - fv->porosity) *
                                                        bf[var]->phi[j];
        }
      }

      var = POR_SINK_MASS;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          pmt->d_conv_flux_dSM[i_pore][0][j] =
              -(mp->u_porous_sink_constants[2] / mp->u_porous_sink_constants[5] * bf[var]->phi[j] *
                mp->u_porous_sink_constants[5] / mp->density) /
              (1.0 - fv->porosity);
        }
      }
      for (b = 0; b < pd->Num_Dim; b++) {
        var = MESH_DISPLACEMENT1 + b;
        if (pd->v[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            /* Need to calculate the dilation of the medium and relate
             * it to the porosity
             * this is really a Kludge for the porosity equation
             * - replace first mass-flux piece with
             * the relationship between porosity and deformation */
            pmt->d_conv_flux_dmesh[i_pore][0][b][j] = fv->d_volume_change_dx[b][j];
          }
        }
      }
    }

    /*
     * SOURCE TERM - set reaction rates to zero for now; constant heat source assumed
     */
    if (pd->e[pg->imtrx][R_POR_SINK_MASS]) {

      var = POR_SINK_MASS;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          pmt->d_MassSource_dSM[i_pl][j] = d_MassSource[var][j];
        }
      }

      var = POR_LIQ_PRES;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          pmt->d_MassSource_dpmv[i_pl][i_pl][j] = d_MassSource[var][j];
        }
      }

      var = POR_POROSITY;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          pmt->d_MassSource_dpmv[i_pl][i_pore][j] = d_MassSource[var][j];
        }
      }

      var = MESH_DISPLACEMENT1;
      if (pd->v[pg->imtrx][var]) {
        for (b = 0; b < pd->Num_Dim; b++) {
          var = MESH_DISPLACEMENT1 + b;
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            pmt->d_MassSource_dmesh[i_pl][b][j] = d_MassSource[var][j];
          }
        }
      }
    }

    /*
     * SUPG TERMS -> these are calculated numerically, due to
     *               their complexity
     */
    if (mp->Porous_wt_funcModel == SUPG) {
      eqn = var = POR_LIQ_PRES;

      /*
       * Install the base state first
       */
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        switch (dim) {
        case 2:
          pmv->d_U_supg_hnorm_dmde[0][j] = pmv->U_supg_hnorm[0];
          pmv->d_U_supg_hnorm_dmde[1][j] = pmv->U_supg_hnorm[1];
          break;
        case 3:
          pmv->d_U_supg_hnorm_dmde[0][j] = pmv->U_supg_hnorm[0];
          pmv->d_U_supg_hnorm_dmde[1][j] = pmv->U_supg_hnorm[1];
          pmv->d_U_supg_hnorm_dmde[2][j] = pmv->U_supg_hnorm[2];
          break;
        }
      }
      /*
       * Turn off calculation of Jacobian terms - we don't know them
       */
      AF_assemble_Residual_Only();
      /*
       * loop over the columns
       */
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        phi_j = bf[var]->phi[j];
        grad_phi_j = bf[var]->grad_phi[j];

        /*
         * Obtain a delta
         */
        value = *(esp->p_liq[j]);
        delta_val = calc_numerical_delta(value);
        /*
         * Calculate changes in the value of interpolated
         * base variables and their derivatives at the gauss point
         */
        fv->p_liq += delta_val * phi_j;
        for (a = 0; a < dim; a++) {
          fv->grad_p_liq[a] += delta_val * grad_phi_j[a];
        }

        /*
         * Calculate changes in transport properties
         * -> duplicate the functionality of load_porous_properties()
         *    for properties that would change as a function of changes
         *    in P_liq
         */
        p_gas_star = mp->u_porous_gas_constants[3];
        cap_pres = pmv->cap_pres = p_gas_star - fv->p_liq;
        d_cap_pres[i_pl] = -1.;
        mp->saturation = load_saturation(mp->porosity, cap_pres, d_cap_pres);
        if (mp->RelLiqPermModel != CONSTANT) {
          load_liq_perm(mp->porosity, cap_pres, mp->saturation, d_cap_pres);
        }
        load_mass_flux(mp->porosity, cap_pres, mp->saturation, d_cap_pres);
        /*
         * Calculate U_supg
         */
        Stab->Grid_Peclet_Number[POR_LIQ_PRES] = get_supg_terms_porous(Stab->hsquared, Stab->hhv);
        dlamdPl = (mp->density * mp->permeability * mp->d_rel_liq_perm[POR_LIQ_PRES]);

        /*
         * Calculate the numerical derivative
         */
        switch (dim) {
        case 2:
          pmv->d_U_supg_hnorm_dmde[0][j] =
              (pmv->U_supg_hnorm[0] - pmv->d_U_supg_hnorm_dmde[0][j]) / delta_val;
          pmv->d_U_supg_hnorm_dmde[1][j] =
              (pmv->U_supg_hnorm[1] - pmv->d_U_supg_hnorm_dmde[1][j]) / delta_val;
          tmp = (fv->grad_p_liq[0] * fv->grad_p_liq[0] + fv->grad_p_liq[1] * fv->grad_p_liq[1]);
          break;
        case 3:
          pmv->d_U_supg_hnorm_dmde[0][j] =
              (pmv->U_supg_hnorm[0] - pmv->d_U_supg_hnorm_dmde[0][j]) / delta_val;
          pmv->d_U_supg_hnorm_dmde[1][j] =
              (pmv->U_supg_hnorm[1] - pmv->d_U_supg_hnorm_dmde[1][j]) / delta_val;
          pmv->d_U_supg_hnorm_dmde[2][j] =
              (pmv->U_supg_hnorm[2] - pmv->d_U_supg_hnorm_dmde[2][j]) / delta_val;
          tmp = (fv->grad_p_liq[0] * fv->grad_p_liq[0] + fv->grad_p_liq[1] * fv->grad_p_liq[1] +
                 fv->grad_p_liq[2] * fv->grad_p_liq[2]);
          break;
        }
        tmp *= -dlamdPl;
        pmt->d_conv_flux_supg_dpmv[i_pl][i_pl][j] = (tmp - pmt->conv_flux_supg[i_pl]) / delta_val;

        /*
         * Restore the base state
         */
        fv->p_liq -= delta_val * phi_j;
        for (a = 0; a < dim; a++) {
          fv->grad_p_liq[a] -= delta_val * grad_phi_j[a];
        }
      }
      /*
       * Further restore the base state
       */
      cap_pres = pmv->cap_pres = p_gas_star - fv->p_liq;
      mp->saturation = load_saturation(mp->porosity, cap_pres, d_cap_pres);
      if (mp->RelLiqPermModel != CONSTANT) {
        load_liq_perm(mp->porosity, cap_pres, mp->saturation, d_cap_pres);
      }
      load_mass_flux(mp->porosity, cap_pres, mp->saturation, d_cap_pres);
      /*
       * Calculate U_supg
       */
      Stab->Grid_Peclet_Number[POR_LIQ_PRES] = get_supg_terms_porous(Stab->hsquared, Stab->hhv);
      /*
       * Restore calculation of Jacobian terms
       */
      AF_restore_Jacobian_Flag();
      /*
       * Complete the current calculation
       */
      for (i = 0; i < ei[pg->imtrx]->dof[var]; i++) {
        grad_phi_i = bf[var]->grad_phi[i];
        switch (dim) {
        case 2:
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            pmt->d_pi_supg_dpmv[i][i_pl][j] = pmv->d_U_supg_hnorm_dmde[0][j] * grad_phi_i[0] +
                                              pmv->d_U_supg_hnorm_dmde[1][j] * grad_phi_i[1];
          }
          break;
        case 3:
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            pmt->d_pi_supg_dpmv[i][i_pl][j] = pmv->d_U_supg_hnorm_dmde[0][j] * grad_phi_i[0] +
                                              pmv->d_U_supg_hnorm_dmde[1][j] * grad_phi_i[1] +
                                              pmv->d_U_supg_hnorm_dmde[2][j] * grad_phi_i[2];
          }
          break;
        }
      }
    } /* end of SUPG */
  }   /* end of if Jacobian */
  return 0;
}
/* end of get_porous_part_sat_terms */
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

/**************************************************************************/
/**************************************************************************/

int get_porous_part_sat_terms_decoupled(struct Porous_Media_Terms *pmt,
                                        double tt, /* time integration scheme param */
                                        double dt) /* current time step size */

/*********************************************************************
 *  get_porous_part_sat_terms_decoupled():
 *
 *
 *   Calculate the capacity, flux, convection, and source terms
 *   for a porous medium WITH THE ASSUMPTION THAT A SOLID-MECHANICS CODE IS
 *   PROVIDING ANY EFFECTS OF DEFORMATION OR POROELASTICITY.   ORIGINAL
 *   IMPLEMENTATION BY PRSCHUN FOR JAS/GOMA COUPLING
 *
 *   Input
 * ---------
 *   pmv  - All terms in this structure are considered to have already
 *          been evaluated at the guass point.
 *
 *   Output
 * --------
 *   pmt  - This routine is responsible for filling out all terms
 *          in this structure.
 **********************************************************************/
{
  int w, w1, i, j, a, var, eqn;
  const int i_pl = 0, i_pg = 1, i_pe = 3;
  double value, delta_val, cap_pres, p_gas_star = 0.0;
  int dim = pd->Num_Dim;
  double *grad_phi_i, dlamdPl, tmp = 0.0, *grad_phi_j = 0;
  double phi_j, d_cap_pres[3];

  /*
   * Intialize all terms in the pmt structure to zero or NULL.
   *  HKM -> note in the future, we may elect to zero only
   *         terms that are actually used.
   */
  zero_structure(pmt, sizeof(struct Porous_Media_Terms), 1);

  /* THIS ROUTINE NEEDS TO BE UPDATED FOR TWO-PHASE OPTION. YOU NEED TO
   * RE-DERIVE THE EQUATIONS ALA' TIME DERIVATIVES FOR DECOUPLED APPROACH AND
   * CHANGE ACCORDINGLY
   */
  if (pd->e[pg->imtrx][R_POR_GAS_PRES])
    GOMA_EH(GOMA_ERROR, "Must update get_porous_part_sat_terms_decouple for two phase");

  /***********************************************************************
   * CAPACITY TERM -
   *     mass density (gm/cm**3) of liquid and gas solvents
   *     in all phases present at the gauss point. also, energy (J/cm**3)
   *     in all phases present, if the energy eqn is to be solved.
   ***********************************************************************/

  pmt->Inventory_solvent[i_pl] = pmv->bulk_density[i_pl];
  pmt->Inventory_solvent[i_pg] = pmv->bulk_density[i_pg];
  pmt->Inventory_solvent[i_pe] = pmv->bulk_density[i_pe];

  if (pd->TimeIntegration != STEADY) {

    /*
     * Porosity is extracted from its external field earlier,
     * so just use as is.
     */
    /* OLD pmt->Inventory_solvent_dot[i_pl] = pmv->bulk_density[i_pl] *
       mp->porous_compressibility * fv_dot->p_liq / fv->porosity/(1. - fv->porosity); */

    pmt->Inventory_solvent_dot[i_pl] =
        pmv->bulk_density[i_pl] * mp->porous_compressibility * fv_dot->p_liq;

    pmt->Inventory_solvent_dot[i_pl] += fv->porosity *
                                        (mp->density - pmv->gas_density_solvents[i_pl]) *
                                        mp->d_saturation[POR_LIQ_PRES] * fv_dot->p_liq;

    /*
     * Gas solvent component, e.g. air
     */
    if (pd->e[pg->imtrx][R_POR_GAS_PRES])
      GOMA_EH(GOMA_ERROR, "Need to implement gas phase transport in decoupled version");
    if (pd->e[pg->imtrx][R_POR_POROSITY])
      GOMA_EH(GOMA_ERROR, "Need to implement solid transport in decoupled version");

  } else {
    pmt->Inventory_solvent_dot[i_pl] = 0.0;
    pmt->Inventory_solvent_dot[i_pg] = 0.0;
  }

  /*
   * CONVECTIVE FLUX TERM
   *
   *  Basically all lumped together with the capacitance term.
   */
  for (a = 0; a < VIM; a++) {
    /* old pmt->conv_flux[i_pl][a] =  0.; */
  }
  pmt->conv_flux[i_pl][0] =
      pmv_old->bulk_density[i_pl] * fv->external_field[efv->ev_dpordt_index] / (1. - fv->porosity);
  /*
   * DIFFUSIVE FLUX TERM - this is the term that is integrated by parts in the
   *   species conservation residual equation.  In a porous medium it contains
   *   five parts of the species flux:
   *    1) convection of the bulk concentration of each species with a moving solid
   *    2) convection of each species in gas phase relative to solid (from Darcy)
   *    3) convection of each species in liquid phase relative to solid (from Darcy)
   *    4) diffusion in gas phase (relative to gas Darcy Flux)
   *    5) diffusion in liquid phase (relative to liquid Darcy Flux)
   *         - currently this is neglected
   */

  if (cr->PorousFluxModel == DARCY_FICKIAN || cr->PorousFluxModel == POWERLAW_DARCY_FICKIAN) {
    /*
     *  Here, vconv is the velocity of the solid phase
     *  rel_mass_flux[][] contains both the Darcy flux and the diffusive
     *  flux components.
     */
    for (a = 0; a < VIM; a++) {
      pmt->diff_flux[i_pl][a] = (pmv->rel_mass_flux[i_pl][a]);
    }
  } else {
    GOMA_EH(GOMA_ERROR, "Unimplemented flux constitutive relation in porous media. "
                        "Check Porous Diffusion Constitutive Equation Card");
  }

  /*
   * SUPG TERMS
   *
   *   pi_supg[i] = contribution to the basis function from SUPG.
   *   conv_flux_supg[i_pl] = New term representing convective upwinding
   *                          contribution.
   */
  if (mp->Porous_wt_funcModel == SUPG) {
    eqn = POR_LIQ_PRES;
    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      grad_phi_i = bf[eqn]->grad_phi[i];
      for (a = 0; a < dim; a++) {
        pmt->pi_supg[i] += pmv->U_supg_hnorm[a] * grad_phi_i[a];
      }
    }
    dlamdPl = (mp->density * mp->permeability * mp->d_rel_liq_perm[POR_LIQ_PRES]);
    pmt->conv_flux_supg[i_pl] = 0.0;
    for (a = 0; a < dim; a++) {
      pmt->conv_flux_supg[i_pl] -= fv->grad_p_liq[a] * fv->grad_p_liq[a];
    }
    pmt->conv_flux_supg[i_pl] *= dlamdPl;
  }

  /********************************************************************************
   * NOW, CALCULATE SENSITIVITIES for the Jacobian, if needed
   ********************************************************************************/

  if (af->Assemble_Jacobian) {

    /*
     * sensitivity of CAPACITY TERM - concentration (volume fraction)
     */
    var = POR_LIQ_PRES;
    for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
      for (w = 0; w < MAX_PMV; w++) {
        for (w1 = 0; w1 < MAX_PMV; w1++) {
          pmt->d_Inventory_solvent_dot_dpmv[w][w1][j] = 0.;
        }
      }
    }
    if (pd->TimeIntegration != STEADY) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        pmt->d_Inventory_solvent_dot_dpmv[i_pl][i_pl][j] +=
            pmv->bulk_density[i_pl] * mp->porous_compressibility * (1. + 2. * tt) *
                bf[var]->phi[j] / dt +
            pmv->d_bulk_density[i_pl][POR_LIQ_PRES] * bf[var]->phi[j] * mp->porous_compressibility *
                fv_dot->p_liq;

        pmt->d_Inventory_solvent_dot_dpmv[i_pl][i_pl][j] +=
            fv->porosity * (mp->density - pmv->gas_density_solvents[i_pl]) *
                (mp->d_d_saturation[POR_LIQ_PRES][POR_LIQ_PRES] * bf[var]->phi[j] * fv_dot->p_liq +
                 (1. + 2. * tt) * mp->d_saturation[POR_LIQ_PRES] * bf[var]->phi[j] / dt) -
            fv->porosity * (pmv->d_gas_density_solvents[i_pl][i_pl] * bf[var]->phi[j]) *
                mp->d_saturation[POR_LIQ_PRES] * fv_dot->p_liq;
        ;
      }
    }

    /*
     * sensitivity of convective flux terms
     */

    /*
     * sensitivity of DIFFUSIVE FLUX TERM - concentration, temp, displacement
     */
    var = POR_LIQ_PRES;
    for (w = 0; w < MAX_PMV; w++) {
      for (a = 0; a < VIM; a++) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          for (w1 = 0; w1 < MAX_PMV; w1++) {
            pmt->d_diff_flux_dpmv[w][a][w1][j] = 0.;
          }
        }
      }
    }
    if (cr->PorousFluxModel == DARCY_FICKIAN || cr->PorousFluxModel == POWERLAW_DARCY_FICKIAN) {
      for (a = 0; a < VIM; a++) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          for (w1 = 0; w1 < MAX_PMV; w1++) {
            pmt->d_diff_flux_dpmv[i_pl][a][w1][j] +=
                pmv->d_rel_mass_flux_dpmv[i_pl][a][w1][j] +
                bf[var]->phi[j] * pmv->gas_darcy_velocity[a] *
                    pmv->d_gas_density_solvents[i_pl][POR_LIQ_PRES + w1];
          }

          if (pd->e[pg->imtrx][R_POR_SINK_MASS])
            pmt->d_diff_flux_dSM[i_pl][a][j] += pmv->d_rel_mass_flux_dSM[i_pl][a][j];
        }
      }
    }

    var = TEMPERATURE;
    if (pd->v[pg->imtrx][var]) {
      if (cr->PorousFluxModel == DARCY_FICKIAN || cr->PorousFluxModel == POWERLAW_DARCY_FICKIAN) {
        for (a = 0; a < VIM; a++) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            pmt->d_diff_flux_dT[i_pl][a][j] = 0.;
            pmt->d_diff_flux_dT[i_pg][a][j] = 0.;
          }
        }
      }
      /* NOTE - the gas concentrations, etc should also be functions of temperature */
      for (a = 0; a < VIM; a++) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          pmt->d_diff_flux_dT[i_pl][a][j] -= pmv->d_rel_mass_flux_dT[i_pl][a][j];
          pmt->d_diff_flux_dT[i_pg][a][j] -= pmv->d_rel_mass_flux_dT[i_pg][a][j];
        }
      }
    }

    /*
     * CONVECTIVE FLUX TERM
     */

    /*
     * SOURCE TERM - set reaction rates to zero for now
     */

    /*
     * SUPG TERMS -> these are calculated numerically, due to
     *               their complexity
     */
    if (mp->Porous_wt_funcModel == SUPG) {
      eqn = var = POR_LIQ_PRES;

      /*
       * Install the base state first
       */
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        switch (dim) {
        case 2:
          pmv->d_U_supg_hnorm_dmde[0][j] = pmv->U_supg_hnorm[0];
          pmv->d_U_supg_hnorm_dmde[1][j] = pmv->U_supg_hnorm[1];
          break;
        case 3:
          pmv->d_U_supg_hnorm_dmde[0][j] = pmv->U_supg_hnorm[0];
          pmv->d_U_supg_hnorm_dmde[1][j] = pmv->U_supg_hnorm[1];
          pmv->d_U_supg_hnorm_dmde[2][j] = pmv->U_supg_hnorm[2];
          break;
        }
      }
      /*
       * Turn off calculation of Jacobian terms - we don't know them
       */
      AF_assemble_Residual_Only();
      /*
       * loop over the columns
       */
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        phi_j = bf[var]->phi[j];
        grad_phi_j = bf[var]->grad_phi[j];

        /*
         * Obtain a delta
         */
        value = *(esp->p_liq[j]);
        delta_val = calc_numerical_delta(value);
        /*
         * Calculate changes in the value of interpolated
         * base variables and their derivatives at the gauss point
         */
        fv->p_liq += delta_val * phi_j;
        for (a = 0; a < dim; a++) {
          fv->grad_p_liq[a] += delta_val * grad_phi_j[a];
        }

        /*
         * Calculate changes in transport properties
         * -> duplicate the functionality of load_porous_properties()
         *    for properties that would change as a function of changes
         *    in P_liq
         */
        p_gas_star = mp->u_porous_gas_constants[3];
        cap_pres = pmv->cap_pres = p_gas_star - fv->p_liq;
        d_cap_pres[i_pl] = -1.;
        mp->saturation = load_saturation(mp->porosity, cap_pres, d_cap_pres);
        if (mp->RelLiqPermModel != CONSTANT) {
          load_liq_perm(mp->porosity, cap_pres, mp->saturation, d_cap_pres);
        }
        load_mass_flux(mp->porosity, cap_pres, mp->saturation, d_cap_pres);
        /*
         * Calculate U_supg
         */
        Stab->Grid_Peclet_Number[POR_LIQ_PRES] = get_supg_terms_porous(Stab->hsquared, Stab->hhv);
        dlamdPl = (mp->density * mp->permeability * mp->d_rel_liq_perm[POR_LIQ_PRES]);

        /*
         * Calculate the numerical derivative
         */
        switch (dim) {
        case 2:
          pmv->d_U_supg_hnorm_dmde[0][j] =
              (pmv->U_supg_hnorm[0] - pmv->d_U_supg_hnorm_dmde[0][j]) / delta_val;
          pmv->d_U_supg_hnorm_dmde[1][j] =
              (pmv->U_supg_hnorm[1] - pmv->d_U_supg_hnorm_dmde[1][j]) / delta_val;
          tmp = (fv->grad_p_liq[0] * fv->grad_p_liq[0] + fv->grad_p_liq[1] * fv->grad_p_liq[1]);
          break;
        case 3:
          pmv->d_U_supg_hnorm_dmde[0][j] =
              (pmv->U_supg_hnorm[0] - pmv->d_U_supg_hnorm_dmde[0][j]) / delta_val;
          pmv->d_U_supg_hnorm_dmde[1][j] =
              (pmv->U_supg_hnorm[1] - pmv->d_U_supg_hnorm_dmde[1][j]) / delta_val;
          pmv->d_U_supg_hnorm_dmde[2][j] =
              (pmv->U_supg_hnorm[2] - pmv->d_U_supg_hnorm_dmde[2][j]) / delta_val;
          tmp = (fv->grad_p_liq[0] * fv->grad_p_liq[0] + fv->grad_p_liq[1] * fv->grad_p_liq[1] +
                 fv->grad_p_liq[2] * fv->grad_p_liq[2]);
          break;
        }
        tmp *= -dlamdPl;
        pmt->d_conv_flux_supg_dpmv[i_pl][i_pl][j] = (tmp - pmt->conv_flux_supg[i_pl]) / delta_val;

        /*
         * Restore the base state
         */
        fv->p_liq -= delta_val * phi_j;
        for (a = 0; a < dim; a++) {
          fv->grad_p_liq[a] -= delta_val * grad_phi_j[a];
        }
      }
      /*
       * Further restore the base state
       */
      cap_pres = pmv->cap_pres = p_gas_star - fv->p_liq;
      mp->saturation = load_saturation(mp->porosity, cap_pres, d_cap_pres);
      if (mp->RelLiqPermModel != CONSTANT) {
        load_liq_perm(mp->porosity, cap_pres, mp->saturation, d_cap_pres);
      }
      load_mass_flux(mp->porosity, cap_pres, mp->saturation, d_cap_pres);
      /*
       * Calculate U_supg
       */
      Stab->Grid_Peclet_Number[POR_LIQ_PRES] = get_supg_terms_porous(Stab->hsquared, Stab->hhv);
      /*
       * Restore calculation of Jacobian terms
       */
      AF_restore_Jacobian_Flag();
      /*
       * Complete the current calculation
       */
      for (i = 0; i < ei[pg->imtrx]->dof[var]; i++) {
        grad_phi_i = bf[var]->grad_phi[i];
        switch (dim) {
        case 2:
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            pmt->d_pi_supg_dpmv[i][i_pl][j] = pmv->d_U_supg_hnorm_dmde[0][j] * grad_phi_i[0] +
                                              pmv->d_U_supg_hnorm_dmde[1][j] * grad_phi_i[1];
          }
          break;
        case 3:
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            pmt->d_pi_supg_dpmv[i][i_pl][j] = pmv->d_U_supg_hnorm_dmde[0][j] * grad_phi_i[0] +
                                              pmv->d_U_supg_hnorm_dmde[1][j] * grad_phi_i[1] +
                                              pmv->d_U_supg_hnorm_dmde[2][j] * grad_phi_i[2];
          }
          break;
        }
      }
    }
  } /* end of if Jacobian */
  return 0;
}
/* end of get_porous_part_sat_terms_decoupled */
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

int get_porous_fully_sat_terms(struct Porous_Media_Terms *pmt, double tt, double dt)

/************************************************************************
 * get_porous_fully_sat_terms()
 * ----------------------------------------------------------------------
 * ----- Calculate the capacity, flux, convection, and source terms
 *       for a saturated porous medium.
 *            -> This routine fills up all of the terms in the
 *               Porous_Media_Terms structure.
 *
 * Input
 * ---------
 * tt,	  parm to vary time integration from
 *		  explicit (tt = 1) to
 *		  implicit (tt = 0)
 * dt     current time step size
 * ----------------------------------------------------------------------
 *  Written by: Richard Cairncross
 *              5/10/95
 *  Revised by: Randy Schunk and Sam Subia - March-May 2001
 ************************************************************************/
{
  int w, j, a, b, p, err, var;
  const int i_pl = 0, i_pore = 2, i_pe = 3;
  double *phi_ptr;
  double n_pow = 1. / 0.5;

  /*
   * storage for call to get_convection_velocity()
   */
  double vconv[MAX_PDIM];     /*Calculated convection velocity */
  double vconv_old[MAX_PDIM]; /*Calculated convection velocity at previous time*/
  CONVECTION_VELOCITY_DEPENDENCE_STRUCT d_vconv_struct;
  CONVECTION_VELOCITY_DEPENDENCE_STRUCT *d_vconv = &d_vconv_struct;

  zero_structure(pmt, sizeof(struct Porous_Media_Terms), 1);
#ifdef DEBUG_HKM
  /*
   * If we are not solving the por_liq_pres equation, doop!
   */
  if (!pd->v[pg->imtrx][POR_LIQ_PRES]) {
    GOMA_EH(GOMA_ERROR, "We shouldn't be here if POR_LIQ_PRESS eqn isn't defined");
  }
#endif
  if (efv->ev_porous_decouple) {
    /*compute porosity using linear pressure expansion model, incremental
     *beyond what comes in from mechanics code porosity
     */

    /*
     * Porosity is extracted from its external field earlier,
     * so just use as is.
     */

    /* pmt->Inventory_solvent_dot[i_pl] =mp->density *
       mp->porous_compressibility * fv_dot->p_liq / (1. - fv->porosity); */

    pmt->Inventory_solvent_dot[i_pl] = mp->density * mp->porous_compressibility * fv_dot->p_liq;

    /* pmt->conv_flux[i_pl][0] = mp->density*fv->porosity*(fv->porosity - 0.947)/(1 -
     * fv->porosity)/5.e-9; */

    pmt->conv_flux[i_pl][0] =
        mp->density * fv->porosity * fv->external_field[efv->ev_dpordt_index] / (1. - fv->porosity);

    for (a = 0; a < VIM; a++) {
      /* pmt->conv_flux[i_pl][a] = 0.; */
    }
  }

  else /* if(!efv->ev_porous_decouple) */
  {
    /*
     * get the convection velocity (it's different for arbitrary and
     * lagrangian meshes)
     */
    if ((pd->MeshMotion == LAGRANGIAN || cr->MeshMotion == DYNAMIC_LAGRANGIAN)) {
      err = belly_flop(elc->lame_mu);
      GOMA_EH(err, "error in belly flop");
      if (neg_elem_volume)
        return (err);
    }

    /* Calculate inertia of mesh if required. PRS side note: no sensitivity
     * with respect to porous media variables here for single compont pore liquids.
     * eventually there may be sensitivity for diffusion in gas phase to p_gas */
    /* Further note:  This call JUST gets the advective lagrangian velocity */

    err = get_convection_velocity(vconv, vconv_old, d_vconv, dt, tt);
    GOMA_EH(err, "Error in calculating effective convection velocity");

    /*
     * CAPACITY TERM -
     *    -> no capacity term for the solid phase
     *
     */
    for (w = 0; w < MAX_PMV; w++) {
      pmt->Inventory_solvent[w] = 0.0;
      pmt->Inventory_solvent_dot[w] = 0.0;
    }

    /*
     * HKM -> Need to update this section for the variable density specification
     */
    if (pd->v[pg->imtrx][POR_POROSITY]) {
      if (pd->TimeIntegration != STEADY) {

        pmt->Inventory_solvent[i_pl] = mp->porosity * mp->density;
        pmt->Inventory_solvent_old[i_pl] = mp_old->porosity * mp->density;
        pmt->Inventory_solvent_dot_old[i_pl] = 0.0;

        /* convert Inventory_dot_old by multiplying by capacitance matrix with chain rule*/
        /* This computes  d(inventory_liq_solvent)/dt = dp_liq/dt* d(porosity*rho)/d(p_liq) +...*/

        pmt->Inventory_solvent_dot_old[i_pl] +=
            fv_dot->p_liq * mp_old->d_porosity[POR_LIQ_PRES] * mp_old->density +
            fv_dot->porosity * mp_old->d_porosity[POR_POROSITY] * mp->density;

        pmt->Inventory_solvent_dot[i_pl] =
            ((1 + 2. * tt) * (pmt->Inventory_solvent[i_pl] - pmt->Inventory_solvent_old[i_pl]) /
                 dt -
             2. * tt * pmt->Inventory_solvent_dot_old[i_pl]);
      } else {
        pmt->Inventory_solvent_dot[i_pl] = 0.0;
      }
    } else {
      pmt->Inventory_solvent[i_pl] = 0.0;
      pmt->Inventory_solvent_dot[i_pl] = 0.0;
    }

    /*
     * CONVECTIVE FLUX TERM
     * Note that div_d_dot = sum_a(grad_d_dot[a][a]).  For convenience, that sum
     * is done in assemble_porous.
     */
    if (pd->v[pg->imtrx][POR_LIQ_PRES] && pd->TimeIntegration != STEADY) {
      for (a = 0; a < VIM; a++) {
        pmt->conv_flux[i_pl][a] = mp->density * mp->porosity * fv->grad_d_dot[a][a];
      }
    }
  } /*End of if(!efv->ev_porous_decouple) */

  /*
   * DIFFUSIVE FLUX TERM - this is the term that is integrated by
   *   parts in the species conservation residual equation.
   *   In a saturated porous medium it contains three parts of
   *   the species flux
   *
   *    1) convection of the bulk concentration of each species
   *       with the solid
   *    2) convection of liquid relative to solid (from Darcy)
   *    3) diffusion in liquid phase (relative to liquid convection)
   *       - currently this is neglected
   */

  if (cr->PorousFluxModel == DARCY_FICKIAN || cr->PorousFluxModel == POWERLAW_DARCY_FICKIAN) {
    /*
     * add convection fluxes into the total diffusion flux term
     * here vconv from get_convection_velocity is
     *   (Imposed Velocity of Solid relative to control volume mesh, v_s)
     *  +(Mesh velocity)
     *
     * or
     *
     *   V_sfs*F + x_dot
     *   (where V_sfs comes from the Advective Lagrangian Velocity Card)
     */
    for (a = 0; a < VIM; a++) {
      pmt->diff_flux[i_pl][a] =
          mp->density * (pmv->liq_darcy_velocity[a] + vconv[a] * mp->porosity);
    }
  } else {
    GOMA_EH(GOMA_ERROR, "Unimplemented flux constitutive relation in saturated media. Check Porous "
                        "Diffusion Constitutive Equation card");
  }

  if (pd->v[pg->imtrx][POR_POROSITY]) {
    /*
     *  Solid phase conservation equation:
     *
     *       det(F) = (1-phi_0) / (1-phi)
     *
     *  Need to calculate the dilation of the medium and
     *  relate it to the porosity.
     *  This is really a Kludge for the porosity equation
     *   - replace first mass-flux piece with the relationship
     *     between porosity and deformation
     **/
    for (a = 0; a < VIM; a++) {
      pmt->conv_flux[i_pore][a] = 0.0;
    }

    if (pd->e[pg->imtrx][R_POR_SINK_MASS]) {
      pmt->conv_flux[i_pore][0] =
          (fv->volume_change -
           (mp->u_porous_sink_constants[3] / mp->u_porous_sink_constants[4] +
            mp->u_porous_sink_constants[2] / mp->u_porous_sink_constants[5] *
                (1 + fv->sink_mass * mp->u_porous_sink_constants[5] / mp->density)) /
               (1.0 - fv->porosity));
    } else {
      pmt->conv_flux[i_pore][0] =
          (fv->volume_change - (1.0 - mp->u_porosity[0]) / (1.0 - fv->porosity));
    }
  } else {
    for (a = 0; a < VIM; a++) {
      pmt->conv_flux[i_pore][a] = 0.0;
    }
  }

  /*
   * SOURCE TERM - set reaction rates to zero for now
   */
  pmt->MassSource[i_pl] = 0.0;
  pmt->MassSource[i_pore] = 0.0;
  pmt->MassSource[i_pe] = 0.0;

  if (pd->e[pg->imtrx][R_POR_SINK_MASS] && mp->PorousSinkConstantsModel == LINEAR) {
    pmt->MassSource[i_pl] -= mp->u_porous_sink_constants[0] * mp->u_porous_sink_constants[2] *
                             (mp->u_porous_sink_constants[1] - fv->sink_mass) /
                             mp->u_porous_sink_constants[1] / fv->volume_change;

    /*  I think this term was mistakingly added in the source document due to their
     *  lack of understanding of Reynolds transport theorem
     */
    /*  for(a=0; a < VIM; a++)
      {
        pmt->MassSource[i_pl] += fv_dot->x[a] * mp->density *
          (fv->grad_porosity[a]);

      }
    */
  } else {
    GOMA_WH(0, "Porous sink equation not being solved, or constants not set");
  }

  /*** Heat Source - constant only for now ****/

  if (mp->HeatSourceModel == USER) {
    GOMA_EH(GOMA_ERROR, "User heat source model not defined");
  } else if (mp->HeatSourceModel == CONSTANT) {
    pmt->MassSource[i_pe] = mp->heat_source;
  } else {
    GOMA_EH(GOMA_ERROR, "Unrecognized heat source model");
  }

  /*
   * NOW, CALCULATE SENSITIVITIES for the Jacobian, if needed
   *
   *  HKM -> Note, all elements of the pmt structure were zeroed before
   *         calling this routine. Therefore, we may assume that
   *         the Jacobian terms are all zero before this section.
   */
  if (af->Assemble_Jacobian) {

    /*Add on sensitivity to pore pressure for decoupled approach */

    if (efv->ev_porous_decouple) {
      /*
       * dependence of liq darcy flux on the liquid pressure gradient
       */
      var = POR_LIQ_PRES;
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        phi_ptr = bf[var]->phi;

        pmt->d_Inventory_solvent_dot_dpmv[i_pl][i_pl][j] =
            mp->density * mp->porous_compressibility * (1. + 2. * tt) * phi_ptr[j] / dt;
      }
    }

    /*
     * sensitivity of CAPACITY TERM - concentration (volume fraction)
     *
     *   dependence of capacity term on the change in porosity:
     */
    var = POR_POROSITY;
    if (pd->TimeIntegration != STEADY && pd->v[pg->imtrx][POR_POROSITY]) {
      phi_ptr = bf[var]->phi;
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        pmt->d_Inventory_solvent_dot_dpmv[i_pl][i_pore][j] =
            (1 + 2. * tt) * phi_ptr[j] * mp->density / dt;
      }
    }

    /*
     * sensitivity of DIFFUSIVE FLUX TERM
     *
     *  - dependence of vconv on the porosity
     */
    if (pd->v[pg->imtrx][POR_POROSITY]) {
      var = POR_POROSITY;
      phi_ptr = bf[var]->phi;
      for (a = 0; a < VIM; a++) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          pmt->d_diff_flux_dpmv[i_pl][a][i_pore][j] =
              mp->density * vconv[a] * mp->d_porosity[POR_POROSITY] * phi_ptr[j];
        }
      }
    }

    /*
     * - Various dependencies on liquid flux term
     */

    if (mp->PermeabilityModel != K_TENSOR && mp->PermeabilityModel != ORTHOTROPIC &&
        mp->PermeabilityModel != SM_TENSOR && mp->PermeabilityModel != KC_TENSOR) {
      for (a = 0; a < VIM; a++) {
        var = POR_POROSITY;
        if (pd->v[pg->imtrx][var]) {
          phi_ptr = bf[var]->phi;
          /*
           * Dependence of the permeability on the porosity
           */
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            pmt->d_diff_flux_dpmv[i_pl][a][i_pore][j] += mp->density * pmv->liq_darcy_velocity[a] *
                                                         mp->d_permeability[POR_POROSITY] /
                                                         mp->permeability * bf[var]->phi[j];
          }
        }

        /*
         * dependence of liq darcy flux on the liquid pressure gradient
         */
        if (cr->PorousFluxModel == DARCY_FICKIAN) {
          var = POR_LIQ_PRES;
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            pmt->d_diff_flux_dpmv[i_pl][a][i_pl][j] -=
                mp->density * mp->permeability / mp->viscosity * bf[var]->grad_phi[j][a];
          }
        } else if (cr->PorousFluxModel == POWERLAW_DARCY_FICKIAN) {
          var = POR_LIQ_PRES;
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            pmt->d_diff_flux_dpmv[i_pl][a][i_pl][j] -=
                mp->density * mp->permeability / mp->viscosity * n_pow *
                pow(fv->grad_p_liq[a], n_pow - 1) * bf[var]->grad_phi[j][a];
          }
        }

        /*
         * dependence of liq darcy flux on the sink mass
         */
        var = POR_SINK_MASS;
        if (pd->v[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            pmt->d_diff_flux_dSM[i_pl][a][j] += mp->density * pmv->liq_darcy_velocity[a] *
                                                mp->d_permeability[POR_SINK_MASS] /
                                                mp->permeability * bf[var]->phi[j];
          }
        }
      }
    } else {
      /*
       *  Permeability is a tensor, and we assume it is a constant
       */
      if (pd->v[pg->imtrx][POR_POROSITY] &&
          (mp->PermeabilityModel == KOZENY_CARMAN || mp->PermeabilityModel == SINK_MASS_PERM)) {
        GOMA_EH(GOMA_ERROR, "Permeability tensor can only have constant components for now");
      }

      var = POR_LIQ_PRES;
      for (a = 0; a < VIM; a++) {
        for (b = 0; b < VIM; b++) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            /*
             * dependence of liq darcy flux on the liquid pressure gradient
             */
            pmt->d_diff_flux_dpmv[i_pl][a][i_pl][j] -=
                mp->density * mp->perm_tensor[a][b] / mp->viscosity * bf[var]->grad_phi[j][b];
          }
        }
      }

      if (pd->v[pg->imtrx][POR_POROSITY]) {
        var = POR_POROSITY;
        for (a = 0; a < VIM; a++) {
          for (b = 0; b < VIM; b++) {
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              /*
               * dependence of liq darcy flux on the liquid pressure gradient
               */
              pmt->d_diff_flux_dpmv[i_pl][a][i_pore][j] -=
                  mp->density * mp->d_perm_tensor[a][b][POR_POROSITY] * bf[var]->phi[j] /
                  mp->viscosity * fv->grad_p_liq[b];
            }
          }
        }
      }
    } /* else tensor */

    for (b = 0; b < pd->Num_Dim; b++) {
      var = MESH_DISPLACEMENT1 + b;
      phi_ptr = bf[var]->phi;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          for (a = 0; a < VIM; a++) {
            /*
             * Dependence of vconv on mesh position
             */
            pmt->d_diff_flux_dmesh[i_pl][a][b][j] +=
                mp->density * d_vconv->X[a][b][j] * mp->porosity;
            /*
             * Dependence of basis function gradients on mesh position
             */
            if (mp->PermeabilityModel != K_TENSOR && mp->PermeabilityModel != ORTHOTROPIC &&
                mp->PermeabilityModel != SM_TENSOR && mp->PermeabilityModel != KC_TENSOR) {
              pmt->d_diff_flux_dmesh[i_pl][a][b][j] -=
                  mp->density * mp->permeability / mp->viscosity * fv->d_grad_p_liq_dmesh[a][b][j];
            } else {
              for (p = 0; p < VIM; p++) {
                pmt->d_diff_flux_dmesh[i_pl][a][b][j] -= mp->density * mp->perm_tensor[a][p] /
                                                         mp->viscosity *
                                                         fv->d_grad_p_liq_dmesh[p][b][j];
                if (mp->PermeabilityModel == ORTHOTROPIC || mp->PermeabilityModel == SM_TENSOR ||
                    mp->PermeabilityModel == KC_TENSOR) {
                  pmt->d_diff_flux_dmesh[i_pl][a][b][j] -= mp->density *
                                                           mp->d_perm_tensor_dx[a][p][b][j] /
                                                           mp->viscosity * fv->grad_p_liq[p];
                }
              }
            }
          }
        }
      }
    }

    /*
     * CONVECTIVE FLUX TERM
     */

    /*
     * sensitivity of the dilation equation on the porosity
     */
    var = POR_POROSITY;
    if (pd->v[pg->imtrx][var]) {
      if (pd->e[pg->imtrx][R_POR_SINK_MASS]) {

        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          pmt->d_conv_flux_dpmv[i_pore][0][i_pore][j] -=
              (mp->u_porous_sink_constants[3] / mp->u_porous_sink_constants[4] +
               mp->u_porous_sink_constants[2] / mp->u_porous_sink_constants[5] *
                   (1 + fv->sink_mass * mp->u_porous_sink_constants[5] / mp->density)) /
              (1. - fv->porosity) / (1. - fv->porosity) * bf[var]->phi[j];
        }
      } else {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          pmt->d_conv_flux_dpmv[i_pore][0][i_pore][j] -= (1. - mp->u_porosity[0]) /
                                                         (1. - fv->porosity) / (1. - fv->porosity) *
                                                         bf[var]->phi[j];
        }
      }
    }

    var = POR_SINK_MASS;
    if (pd->v[pg->imtrx][var]) {

      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++)

      {
        pmt->d_conv_flux_dSM[i_pore][0][j] -=
            (mp->u_porous_sink_constants[2] / mp->u_porous_sink_constants[5] * bf[var]->phi[j] *
             mp->u_porous_sink_constants[5] / mp->density) /
            (1.0 - fv->porosity);
      }
    }

    /* Dependence of convective term of p_liq equation wrt porosity */
    var = POR_POROSITY;
    if (pd->v[pg->imtrx][var] && pd->TimeIntegration != STEADY) {
      for (a = 0; a < VIM; a++) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          pmt->d_conv_flux_dpmv[i_pl][a][i_pore][j] +=
              mp->density * fv->grad_d_dot[a][a] * bf[var]->phi[j];
        }
      }
    }

    /*
     * sensitivities of p_liq equation convective term
     * wrt mesh position
     */

    if (pd->v[pg->imtrx][POR_POROSITY] && pd->TimeIntegration != STEADY) {
      for (b = 0; b < pd->Num_Dim; b++) {
        var = MESH_DISPLACEMENT1 + b;
        phi_ptr = bf[var]->phi;
        if (pd->v[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            for (a = 0; a < VIM; a++) {
              /*
               * Dependence of rho*porosity*div_d_dot on the mesh position and
               * note that d_grad_dot_dmesh trace is incomplete. See load_fv_mesh_derivs
               */
              pmt->d_conv_flux_dmesh[i_pl][a][b][j] +=
                  mp->density * mp->porosity * fv->d_grad_d_dot_dmesh[a][a][b][j];
            }
          }
        }
      }
    }

    /*Dependence of porosity equation wrt mesh */

    if (pd->v[pg->imtrx][POR_POROSITY]) {
      for (b = 0; b < pd->Num_Dim; b++) {
        var = MESH_DISPLACEMENT1 + b;
        phi_ptr = bf[var]->phi;
        if (pd->v[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            /*
             * Dependence of dilation on the mesh positions
             */
            pmt->d_conv_flux_dmesh[i_pore][0][b][j] = fv->d_volume_change_dx[b][j];
          }
        }
      }
    }

    /* Source TERM Sensitivities */

    /*
     * SOURCE TERM - set reaction rates to zero for now; constant heat source assumed
     */
    if (pd->e[pg->imtrx][R_POR_SINK_MASS] && mp->PorousSinkConstantsModel == LINEAR) {
      var = POR_SINK_MASS;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          pmt->d_MassSource_dSM[i_pl][j] += mp->u_porous_sink_constants[0] *
                                            mp->u_porous_sink_constants[2] * bf[var]->phi[j] /
                                            mp->u_porous_sink_constants[1] / fv->volume_change;
        }
      }

      for (b = 0; b < pd->Num_Dim; b++) {
        var = R_MESH1 + b;
        if (pd->v[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            pmt->d_MassSource_dmesh[i_pl][b][j] +=
                fv->d_volume_change_dx[b][j] * mp->u_porous_sink_constants[0] *
                mp->u_porous_sink_constants[2] * (mp->u_porous_sink_constants[1] - fv->sink_mass) /
                mp->u_porous_sink_constants[1] / fv->volume_change / fv->volume_change;

            /*  I think this term was mistakingly added in the source document due to their
             *  lack of understanding of Reynolds transport theorem
             */
            /* pmt->d_MassSource_dmesh[i_pl][b][j] += ((1.+2.*tt) * bf[var]->phi[j]/dt) *
              mp->density *
              (fv->grad_porosity[b]);
            */

            /*
              for ( a=0; a < pd->Num_Dim; a++)
              {
                pmt->d_MassSource_dmesh[i_pl][b][j] += fv_dot->x[a] * mp->density *
                  (fv->d_grad_porosity_dmesh[a][b][j]);
              }
            */
          }
        }
      }

      var = POR_LIQ_PRES;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          for (a = 0; a < VIM; a++) {
            pmt->d_MassSource_dpmv[i_pl][i_pl][j] =
                fv_dot->x[a] * mp->density *
                (mp->porosity * mp->d_saturation[POR_LIQ_PRES] * bf[var]->grad_phi[j][a]);
          }
          /*PRS Caveate: I think an additional sensitivity to of sat(p_liq) needs to be in here */
        }
      }
      var = POR_POROSITY;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          for (a = 0; a < VIM; a++) {
            pmt->d_MassSource_dpmv[i_pl][i_pore][j] +=
                fv_dot->x[a] * mp->density * (bf[var]->grad_phi[j][a]);
          }
        }
      }
    }

  } /* end of if Jacobian */
  return 0;
}
/* end of get_porous_fully_sat_terms */
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void porous_mass_flux_surf_bc(double func[DIM],
                              double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                              int wspec, /* species number of this boundary condition */
                              double mass_tran_coeff,
                              /* Mass transfer coefficient in units consistent with
                               * gas phase concentration driving force */
                              double Y_c, /* bath concentration 	            */
                              double mass_tran_coeff1,
                              /* Mass transfer coefficient for forced liquid
                                 exciting a partially saturated domain */
                              double p_0, /* sink pressure for liquid extraction */
                              dbl dt,     /* current value of the time step          */
                              dbl tt)     /* parameter to vary time integration from */
/***********************************************************************
 *
 *  Function which calculates the surface integral for convective
 *  mass transfer from a porous medium.
 *
 ***********************************************************************/
{
  int j_id, w1, dim, a, p;
  int var;
  double phi_j;
  int err;
  int first_porous_var = POR_LAST - MAX_POROUS_NUM + 1;
  double vconv[MAX_PDIM];     /*Calculated convection velocity */
  double vconv_old[MAX_PDIM]; /*Calculated convection velocity at previous time*/
  CONVECTION_VELOCITY_DEPENDENCE_STRUCT d_vconv_struct;
  CONVECTION_VELOCITY_DEPENDENCE_STRUCT *d_vconv = &d_vconv_struct;

  if (af->Assemble_LSA_Mass_Matrix)
    return;

  dim = pd->Num_Dim;

  /* get deformation gradient, etc. if lagrangian mesh with inertia */
  if ((pd->MeshMotion == LAGRANGIAN || cr->MeshMotion == DYNAMIC_LAGRANGIAN)) {
    err = belly_flop(elc->lame_mu);
    GOMA_EH(err, "error in belly flop");
    if (neg_elem_volume)
      return;
    err = get_convection_velocity(vconv, vconv_old, d_vconv, dt, tt);
    GOMA_EH(err, "Error in calculating effective convection velocity");
  }

  /* get a mass transfer coefficient from gas phase calculations */

  if (af->Assemble_Residual) {
    *func -= mp->porosity * mass_tran_coeff * (pmv->gas_density_solvents[0] - Y_c);
    for (p = 0; p < dim; p++) {
      *func -= fv->snormal[p] * vconv[p] * pmv->bulk_density[0];
    }
    if (fv->p_liq >= p_0) {
      *func -= mp->porosity * mass_tran_coeff1 * (fv->p_liq - p_0);
    }
  }

  if (af->Assemble_Jacobian) {
    /* sum the contributions to the global stiffness matrix
     * for Porous Media variables
     */

    /*
     * J_pmv_pmv
     */
    var = POR_LIQ_PRES;
    for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++) {
      phi_j = bf[var]->phi[j_id];

      /* WARNING d_gas_density_solvents goes over first_porous_var to last_porous_var */
      for (w1 = 0; w1 < pd->Num_Porous_Eqn; w1++) {
        d_func[0][first_porous_var + w1][j_id] -=
            (mp->porosity * pmv->d_gas_density_solvents[0][first_porous_var + w1] +
             mp->d_porosity[first_porous_var + w1] * (pmv->gas_density_solvents[0] - Y_c)) *
            mass_tran_coeff * phi_j;
        for (p = 0; p < dim; p++) {
          d_func[0][first_porous_var + w1][j_id] -=
              fv->snormal[p] * (d_vconv->v[p][w1][j_id] * pmv->bulk_density[0] +
                                vconv[p] * pmv->d_bulk_density[0][first_porous_var + w1] * phi_j);
        }
      }
      if (fv->p_liq >= p_0) {
        d_func[0][first_porous_var][j_id] -= mp->porosity * mass_tran_coeff1 * phi_j;
      }
    }

    /*
     * J_pmv_T
     */
    var = TEMPERATURE;
    if (pd->v[pg->imtrx][var]) {
      for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++) {
        phi_j = bf[var]->phi[j_id];

        d_func[0][var][j_id] -=
            mp->porosity * mass_tran_coeff * pmv->d_gas_density_solvents[0][var] * phi_j;
        for (p = 0; p < dim; p++) {
          d_func[0][var][j_id] -= fv->snormal[p] * (d_vconv->T[p][j_id] * pmv->bulk_density[0] +
                                                    vconv[p] * pmv->d_bulk_density[0][var] * phi_j);
        }
      }
    }

    /*
     * J_pmv_d
     */

    for (a = 0; a < dim; a++) {
      var = MESH_DISPLACEMENT1 + a;
      if (pd->v[pg->imtrx][var]) {
        for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++) {
          phi_j = bf[var]->phi[j_id];

          for (p = 0; p < dim; p++) {
            d_func[0][var][j_id] -=
                pmv->bulk_density[0] *
                (fv->dsnormal_dx[p][a][j_id] * vconv[p] + fv->snormal[p] * d_vconv->X[p][a][j_id]);
          }
        }
      }
    }

  } /* End of if Assemble_Jacobian */

} /* end of routine porous_mass_flux_surf_bc                                                */
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void porous_convection_bc(double func[DIM],
                          double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                          int wspec, /* species number of this boundary condition    */
                          dbl dt,    /* current value of the time step               */
                          dbl tt)    /* parameter to vary time integration from */
/*******************************************************************************
 *
 *  Function which calculates the surface integral for convective
 *  mass transfer from a porous medium.
 *
 *******************************************************************************/
{
  int j_id, w1, dim, a, p, i_pore, var;
  double phi_j;
  int err;
  int first_porous_var = POR_LAST - MAX_POROUS_NUM + 1;
  double vconv[MAX_PDIM];     /*Calculated convection velocity */
  double vconv_old[MAX_PDIM]; /*Calculated convection velocity at previous time*/
  CONVECTION_VELOCITY_DEPENDENCE_STRUCT d_vconv_struct;
  CONVECTION_VELOCITY_DEPENDENCE_STRUCT *d_vconv = &d_vconv_struct;
  /***************************** EXECUTION BEGINS *******************************/
  if (af->Assemble_LSA_Mass_Matrix)
    return;

  dim = pd->Num_Dim;
  i_pore = 2;
  /* get deformation gradient, etc. if lagrangian mesh with inertia */
  if (pd->MeshMotion == LAGRANGIAN && pd->MeshInertia == 1) {
    err = belly_flop(elc->lame_mu);
    GOMA_EH(err, "error in belly flop");
    if (neg_elem_volume)
      return;
    err = get_convection_velocity(vconv, vconv_old, d_vconv, dt, tt);
    GOMA_EH(err, "Error in calculating effective convection velocity");
  }

  if (af->Assemble_Residual) {
    /*N.B. PRS:  May need an x_dot on here.  Check Rich's work */
    if (mp->PorousMediaType == POROUS_SATURATED) {
      for (p = 0; p < dim; p++) {
        *func += fv->snormal[p] * vconv[p] * mp->density;
      }
    } else {
      for (p = 0; p < dim; p++) {
        *func -= fv->snormal[p] * vconv[p] * pmv->bulk_density[0];
      }
    }
  }

  if (af->Assemble_Jacobian) {
    /* sum the contributions to the global stiffness matrix
     * for Porous Media variables
     */

    /*
     * J_s_c
     */
    var = POR_LIQ_PRES;
    if (mp->PorousMediaType == POROUS_SATURATED) {
      if (i_pore != -1) {
        for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++) {
          phi_j = bf[var]->phi[j_id];

          for (p = 0; p < dim; p++) {
            d_func[0][first_porous_var + i_pore][j_id] +=
                fv->snormal[p] * mp->density * (d_vconv->v[p][i_pore][j_id]);
          }
        }
      }
    } else {
      for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++) {
        phi_j = bf[var]->phi[j_id];

        /* WARNING: bulk_density goes from first_porous_var to last_porous_var */

        for (w1 = 0; w1 < pd->Num_Porous_Eqn; w1++) {
          for (p = 0; p < dim; p++) {
            d_func[0][first_porous_var + w1][j_id] -=
                fv->snormal[p] * (d_vconv->v[p][w1][j_id] * pmv->bulk_density[0] +
                                  vconv[p] * pmv->d_bulk_density[0][first_porous_var + w1] * phi_j);
          }
        }
      }
    }

    /*
     * J_s_T
     */
    var = TEMPERATURE;
    if (pd->v[pg->imtrx][var]) {
      if (mp->PorousMediaType != POROUS_SATURATED) {
        for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++) {
          phi_j = bf[var]->phi[j_id];

          for (p = 0; p < dim; p++) {
            d_func[0][var][j_id] -=
                fv->snormal[p] * (d_vconv->T[p][j_id] * pmv->bulk_density[0] +
                                  vconv[p] * pmv->d_bulk_density[0][var] * phi_j);
          }
        }
      }
    }

    /*
     * J_s_d
     */

    for (a = 0; a < dim; a++) {
      var = MESH_DISPLACEMENT1 + a;
      if (pd->v[pg->imtrx][var]) {
        if (mp->PorousMediaType == POROUS_SATURATED) {
          for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++) {
            phi_j = bf[var]->phi[j_id];

            for (p = 0; p < dim; p++) {
              d_func[0][var][j_id] += mp->density * (fv->dsnormal_dx[p][a][j_id] * vconv[p] +
                                                     fv->snormal[p] * d_vconv->X[p][a][j_id]);
            }
          }
        } else {
          for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++) {
            phi_j = bf[var]->phi[j_id];

            for (p = 0; p < dim; p++) {
              d_func[0][var][j_id] -=
                  pmv->bulk_density[0] * (fv->dsnormal_dx[p][a][j_id] * vconv[p] +
                                          fv->snormal[p] * d_vconv->X[p][a][j_id]);
            }
          }
        }
      }
    }

  } /* End of if Assemble_Jacobian */

} /* end of routine porous_convection_bc                                        */
/********************************************************************************/
/********************************************************************************/
/********************************************************************************/

void porous_normal_velocity_bc(double func[],
                               double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                               double x_dot[DIM], /* mesh velocity vector        */
                               dbl tt,            /* parameter to vary time integration from
                                                   * explicit (tt = 1) to implicit (tt = 0) */
                               dbl dt,            /* current value of the time step         */
                               int bc_input_id,
                               struct Boundary_Condition *BC_Types,
                               int eb_mat_solid, /* elem block id of porous phase  */
                               int eb_mat_fluid, /* elem block id of gas phase     */
                               int i_pl,         /*HARDWIRED to liquid-solvent
                                                   vapor in gas for now until
                                                   multiple species is allowed. */
                               double dens_vap)  /* density of pure solvent vapor */

/*******************************************************************************
 *
 *  Function which evaluates the kinematic boundary condition
 * for adjacent porous medium and fluid phases
 *			 Author: R. A. Cairncross (3/21/96)
 *******************************************************************************/
{
  int j_id, var, jvar, a, w, dim = pd->Num_Dim;
  double phi_j;
  if (af->Assemble_LSA_Mass_Matrix)
    return;

  /***************************** SOLID SIDE *******************************/
  /*
   *  If current material is the solid phase, calculate normal flux
   *  of solvent through the porous medium
   */
  if (Current_EB_ptr->Elem_Blk_Id == eb_mat_solid) {
    /* add convection fluxes into the total diffusion flux term
     * here vconv is the velocity of the solid phase */
    for (a = 0; a < VIM; a++) {
      *func = pmv->rel_mass_flux[i_pl][a] / dens_vap * fv->snormal[a];
    }

    if (af->Assemble_Jacobian) {

      /* sum the contributions to the global stiffness matrix */

      var = TEMPERATURE;
      for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++) {
        if (pd->v[pg->imtrx][var]) {
          phi_j = bf[var]->phi[j_id];
          /*     d( )/dx        */
          for (a = 0; a < VIM; a++) {
            d_func[0][var][j_id] -=
                pmv->d_rel_mass_flux_dT[i_pl][a][j_id] / dens_vap * fv->snormal[a];
          }
        }
      }

      var = POR_LIQ_PRES;
      for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++) {
        if (pd->v[pg->imtrx][var]) {
          phi_j = bf[var]->phi[j_id];
          for (a = 0; a < VIM; a++) {
            for (w = 0; w < pd->Num_Porous_Eqn; w++) {
              /*     d( )/dx        */
              d_func[0][POR_LIQ_PRES + w][j_id] -=
                  pmv->d_rel_mass_flux_dpmv[i_pl][a][w][j_id] / dens_vap * fv->snormal[a];
            }
          }
        }
      }

      for (jvar = 0; jvar < dim; jvar++) {
        var = MESH_DISPLACEMENT1 + jvar;
        for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++) {
          if (pd->v[pg->imtrx][var]) {
            phi_j = bf[var]->phi[j_id];
            for (a = 0; a < VIM; a++) {
              /*     d( )/dx        */
              d_func[0][var][j_id] -=
                  pmv->d_rel_mass_flux_dpmv[i_pl][a][jvar][j_id] / dens_vap * fv->snormal[a] +
                  pmv->rel_mass_flux[i_pl][a] / dens_vap * fv->dsnormal_dx[a][jvar][j_id];
            }
          }
        }
      }
    }
  }

  /***************************** FLUID SIDE *******************************/
  /*
   *  If current material is the fluid phase, calculate normal fluid velocity
   *  in liquid at the interface
   */
  else if (Current_EB_ptr->Elem_Blk_Id == eb_mat_fluid) {
    /* Calculate the residual contribution from the velocity */
    *func = 0;
    for (a = 0; a < VIM; a++) {
      *func += fv->v[a] * fv->snormal[a];
    }

    if (af->Assemble_Jacobian) {

      /* sum the contributions to the global stiffness matrix */

      for (jvar = 0; jvar < dim; jvar++) {
        var = VELOCITY1 + jvar;
        for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++) {
          if (pd->v[pg->imtrx][var]) {
            phi_j = bf[var]->phi[j_id];
            d_func[0][var][j_id] += phi_j * fv->snormal[jvar];
          }
        }
      }

      for (jvar = 0; jvar < dim; jvar++) {
        var = MESH_DISPLACEMENT1 + jvar;
        for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++) {
          if (pd->v[pg->imtrx][var]) {
            for (a = 0; a < dim; a++) {
              phi_j = bf[var]->phi[j_id];
              d_func[0][var][j_id] += fv->v[a] * fv->dsnormal_dx[a][jvar][j_id];
            }
          }
        }
      }
    }

  } else
    GOMA_EH(GOMA_ERROR, "No Slip called with incorrect block id");
  /* note, we may not want to quit at this point */
  return;
} /*  end of porous_normal_velocity_bc */
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void put_gas_flux_in_pores(double func[],
                           double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                           double x_dot[DIM], /* mesh velocity vector        */
                           dbl tt,            /* parameter to vary time integration from
                                               * explicit (tt = 1) to implicit (tt = 0) */
                           dbl dt,            /* current value of the time step         */
                           int bc_input_id,
                           struct Boundary_Condition *BC_Types,
                           int eb_mat_solid, /* elem block id of porous phase  */
                           int eb_mat_fluid, /* elem block id of gas phase     */
                           int wspec,
                           double dens_vap, /* density of pure solvent vapor */
                           double vapor_recoil)

/******************************************************************************
 *
 *  Function which equates the normal flux of solvent
 * for adjacent porous medium and fluid phases
 *			 Author: R. A. Cairncross (3/22/96)
 ******************************************************************************/
{
  int j_id, var, jvar, a, w, dim = pd->Num_Dim;

  if (af->Assemble_LSA_Mass_Matrix)
    return;

  /***************************** SOLID SIDE *******************************/
  /*
   *  If current material is the solid phase, calculate normal flux
   *  of solvent through the porous medium
   */
  if (Current_EB_ptr->Elem_Blk_Id == eb_mat_solid) {
    /* add convection fluxes into the total diffusion flux term
     * here vconv is the velocity of the solid phase */
    for (a = 0; a < VIM; a++) {
      *func = pmv->rel_mass_flux[wspec][a] / dens_vap * fv->snormal[a];
    }

    if (af->Assemble_Jacobian) {

      /* sum the contributions to the global stiffness matrix */
      var = POR_LIQ_PRES;
      for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++) {
        if (pd->v[pg->imtrx][var]) {
          for (a = 0; a < VIM; a++) {
            for (w = 0; w < pd->Num_Porous_Eqn; w++) {
              d_func[0][POR_LIQ_PRES + w][j_id] +=
                  pmv->d_rel_mass_flux_dpmv[wspec][a][w][j_id] / dens_vap * fv->snormal[a];
            }
          }
        }
      }

      for (jvar = 0; jvar < dim; jvar++) {
        var = MESH_DISPLACEMENT1 + jvar;
        for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++) {
          if (pd->v[pg->imtrx][var]) {
            for (a = 0; a < VIM; a++) {
              /*     d( )/dx        */
              d_func[0][var][j_id] +=
                  pmv->d_rel_mass_flux_dmesh[wspec][a][jvar][j_id] / dens_vap * fv->snormal[a] +
                  pmv->rel_mass_flux[wspec][a] / dens_vap * fv->dsnormal_dx[a][jvar][j_id];
            }
          }
        }
      }
    }
  } /* end of mat solid */

  /***************************** FLUID SIDE *******************************/
  /*
   *  If current material is the fluid phase, calculate normal fluid velocity
   *  in liquid at the interface
   */
  else if (Current_EB_ptr->Elem_Blk_Id == eb_mat_fluid) {
    /* Calculate the residual contribution from the velocity */
    *func = 0;
    for (a = 0; a < VIM; a++) {
      *func +=
          fv->snormal[a] * (-mp->thermal_conductivity * fv->grad_T[a]) / (1 - vapor_recoil * fv->T);
    }

    if (af->Assemble_Jacobian) {

      /* sum the contributions to the global stiffness matrix */
      for (jvar = 0; jvar < dim; jvar++) {
        var = MESH_DISPLACEMENT1 + jvar;
        for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++) {
          if (pd->v[pg->imtrx][var]) {
            for (a = 0; a < dim; a++) {
              d_func[0][var][j_id] -= (fv->dsnormal_dx[a][jvar][j_id] * fv->grad_T[a] +
                                       fv->snormal[a] * fv->d_grad_T_dmesh[a][jvar][j_id]) *
                                      mp->thermal_conductivity / (1 - vapor_recoil * fv->T);
            }
          }
        }
      }

      var = TEMPERATURE;
      for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++) {
        if (pd->v[pg->imtrx][var]) {
          for (a = 0; a < dim; a++) {
            d_func[0][var][j_id] -=
                mp->thermal_conductivity * fv->snormal[a] *
                (fv->grad_T[a] / (1 - vapor_recoil * fv->T) * bf[var]->phi[j_id] +
                 bf[var]->grad_phi[j_id][a]) /
                (1 - vapor_recoil * fv->T);
          }
        }
      }
    }

  } else
    GOMA_EH(GOMA_ERROR, "No Slip called with incorrect block id");
  /* note, we may not want to quit at this point */
  return;
} /*  end of  put_gas_flux_in_pores  */
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void porous_vapor_equil_bc(double func[],
                           double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                           double x_dot[DIM], /* mesh velocity vector        */
                           dbl tt,            /* parameter to vary time integration from
                                               * explicit (tt = 1) to implicit (tt = 0) */
                           dbl dt,            /* current value of the time step         */
                           int bc_input_id,
                           struct Boundary_Condition *BC_Types,
                           int eb_mat_solid, /* elem block id of porous phase  */
                           int eb_mat_fluid, /* elem block id of gas phase     */
                           int wspec,
                           double amb_pres) /* ambient pressure              */

/***************************************************************************
 *
 *  Function which evaluates vapor-liquid equilibrium
 * for adjacent porous medium and fluid phases using different variables
 * in each phase
 *			 Author: R. A. Cairncross (3/21/96)
 ***************************************************************************/
{
  int j_id, var, w;
  double phi_j;
  if (af->Assemble_LSA_Mass_Matrix)
    return;

  /***************************** SOLID SIDE *******************************/
  if (Current_EB_ptr->Elem_Blk_Id == eb_mat_solid) {
    /*
     *  calculate volume fraction of solvent in the vapor phase on
     *  porous side of interface
     */
    *func = mp->porous_vapor_pressure[wspec] / amb_pres;

    if (af->Assemble_Jacobian) {

      var = POR_LIQ_PRES;
      for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++) {
        if (pd->v[pg->imtrx][var]) {
          phi_j = bf[var]->phi[j_id];
          for (w = 0; w < pd->Num_Porous_Eqn; w++) {
            d_func[0][POR_LIQ_PRES + w][j_id] +=
                mp->d_porous_vapor_pressure[wspec][POR_LIQ_PRES + w] / amb_pres * phi_j;
          }
        }
      }
    }
  }
  /***************************** FLUID SIDE *******************************/
  else if (Current_EB_ptr->Elem_Blk_Id == eb_mat_fluid) {
    /*
     * HKM-> I believe this should be the mole fraction of air in the
     *       the continuous phase. However, don't have time to fix
     *       this at the moment.
     */
    *func = -fv->T;
    GOMA_EH(GOMA_ERROR, "porous_vapor_equil_bc is broken. Fix it first");

    if (af->Assemble_Jacobian) {
      var = TEMPERATURE;
      for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++) {
        if (pd->v[pg->imtrx][var]) {
          phi_j = bf[var]->phi[j_id];
          d_func[0][var][j_id] -= phi_j;
        }
      }
    }
  } else
    GOMA_EH(GOMA_ERROR, "porous_vapor_equil_bc");
  return;
} /*  end of   porous_vapor_equil_bc  */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

double load_permeability(void)

/*************************************************************************
 *
 * load_permeability()
 *
 *  input:
 *  -----
 *    mp->porosity
 *
 *  output:
 *  -----
 *    mp->permeability
 *    mp->d_permeability[]
 *
 *        also, for psd models:
 *    pmv->r_pore
 *    pmv->d_r_pore[]
 *    pmv->d_d_r_pore[][]
 *************************************************************************/
{
  int i, j;
  double factor;

  if (mp->PermeabilityModel == CONSTANT || mp->PermeabilityModel == K_TENSOR) {

    /*
     * do nothing
     */
  } else if (mp->PermeabilityModel == KOZENY_CARMAN) {
    /*FOR KOZENY_CARMAN - Porosity dependence on permeability with
     *mp->u_permeability[0] Kozeny-Carman const (=5 for particulate beds)
     *mp->u_permeability[1] surface area per unity volume
     */
    mp->permeability = pow(mp->porosity, 3.0) / pow(1. - mp->porosity, 2.0) /
                       mp->u_permeability[0] / pow(mp->u_permeability[1], 2.0);
    /* Permeablity only depends on porosity currently */
    if (pd->v[pg->imtrx][POR_POROSITY]) {
      mp->d_permeability[POR_POROSITY] =
          (pow(1. - mp->porosity, 2.0) * mp->u_permeability[0] * pow(mp->u_permeability[1], 2.0) *
               3.0 * pow(mp->porosity, 2.0) +
           pow(mp->porosity, 3.0) * 2.0 * mp->u_permeability[0] * pow(mp->u_permeability[1], 2.0) *
               (1.0 - mp->porosity)) /
          (pow(pow(1. - mp->porosity, 2.0) * mp->u_permeability[0] *
                   pow(mp->u_permeability[1], 2.0),
               2.0));
    }
  } else if (mp->PermeabilityModel == SINK_MASS_PERM) {
    /*FOR particle SINK_MASS_PERM - Porosity dependence on permeability with
     *mp->u_permeability[0] k_max (maximum saturated permeability or conductivity)
     *mp->u_permeability[1] k_rmin is the relative minumum permeability
     *mp->u_permeability[2] krscale is the relative scaled permeability
     *mp->u_permeability[3] k_exp is the powerlaw exponent

     * NB: For some reason this model NaNs real easily. Must initialize sink_mass
     * to nonzero value for one thing.   Need to investigate.
     */
    mp->permeability = mp->u_permeability[0] *
                       (mp->u_permeability[1] + (1. - mp->u_permeability[1]) /
                                                    (1 + pow(fv->sink_mass * mp->u_permeability[2],
                                                             mp->u_permeability[3])));

    /* Permeablity only depends on porosity currently */
    if (pd->v[pg->imtrx][POR_SINK_MASS]) {
      mp->d_permeability[POR_SINK_MASS] =
          mp->u_permeability[0] * (mp->u_permeability[1] - 1.) * mp->u_permeability[3] *
          pow(fv->sink_mass * mp->u_permeability[2], mp->u_permeability[3] - 1.0) *
          mp->u_permeability[2] /
          (1 + pow(fv->sink_mass * mp->u_permeability[2], mp->u_permeability[3])) /
          (1 + pow(fv->sink_mass * mp->u_permeability[2], mp->u_permeability[3]));
    }
  } else if (mp->PermeabilityModel == USER) {
    usr_permeability(mp->u_permeability);
  } else if (mp->PermeabilityModel == PSD_VOL || mp->PermeabilityModel == PSD_WEXP ||
             mp->PermeabilityModel == PSD_SEXP) {
    /* FOR PSD_VOL - pore-size distribution function with equal
     * volume in all pore sizes
     *  mp->u_permeability[0] is the porosity in the undeformed state
     *  mp->u_permeability[1] is the pore size in the undeformed state
     *  mp->u_permeability[2] is the fractional size of smallest pore/largest pore
     *  mp->u_permeability[3] is the tortuosity and geometric fudge factor
     */
    /* calculate the current maximum (or mean) pore radius as a function
     * of the porosity */
    factor = mp->u_permeability[1] * pow(1. - mp->u_permeability[0], 1. / 3.) /
             pow(mp->u_permeability[0], 0.5);

    pmv->r_pore = factor * pow(mp->porosity, 0.5) / pow(1. - mp->porosity, 1. / 3.);

    for (i = 0; i < MAX_VARIABLE_TYPES + MAX_CONC; i++) {
      pmv->d_r_pore[i] = 0.;
      for (j = 0; j < MAX_VARIABLE_TYPES + MAX_CONC; j++) {
        pmv->d_d_r_pore[i][j] = 0.;
      }
    }

    /* Permeablity only depends on porosity currently */
    if (pd->v[pg->imtrx][POR_POROSITY]) {
      pmv->d_r_pore[POR_POROSITY] =
          factor * mp->d_porosity[POR_POROSITY] *
          (0.5 * pow(mp->porosity, -0.5) / pow(1. - mp->porosity, 1. / 3.) +
           1. / 3. * pow(mp->porosity, 0.5) / pow(1. - mp->porosity, 4. / 3.));
      pmv->d_d_r_pore[POR_POROSITY][POR_POROSITY] =
          factor * mp->d_porosity[POR_POROSITY] *
          (-0.25 * pow(mp->porosity, -1.5) / pow(1. - mp->porosity, 1. / 3.) +
           2. / 6. * pow(mp->porosity, -0.5) / pow(1. - mp->porosity, 4. / 3.) +
           4. / 9. * pow(mp->porosity, 0.5) / pow(1. - mp->porosity, 7. / 3.));
    }

    if (pd->TimeIntegration == TRANSIENT) {
      pmv_old->r_pore = factor * pow(mp_old->porosity, 0.5) / pow(1. - mp_old->porosity, 1. / 3.);

      for (i = 0; i < MAX_VARIABLE_TYPES + MAX_CONC; i++) {
        pmv_old->d_r_pore[i] = 0.;
      }

      if (pd->v[pg->imtrx][POR_POROSITY]) {
        pmv_old->d_r_pore[POR_POROSITY] =
            factor * mp_old->d_porosity[POR_POROSITY] *
            (0.5 * pow(mp_old->porosity, -0.5) / pow(1. - mp_old->porosity, 1. / 3.) +
             1. / 3. * pow(mp_old->porosity, 0.5) / pow(1. - mp_old->porosity, 4. / 3.));
      }
    }

    /* use pore size to calculate the permeability */
    if (mp->PermeabilityModel == PSD_VOL) {
      factor = mp->u_permeability[3] / 60. * (1. - pow(mp->u_permeability[2], 3.0)) /
               (1. - mp->u_permeability[2]);
    } else if (mp->PermeabilityModel == PSD_WEXP) {
      factor = mp->u_permeability[3] * 5. / 24.;
    } else if (mp->PermeabilityModel == PSD_SEXP) {
      factor = mp->u_permeability[3] * 32. / 9. / 24. / M_PIE;
    }
    mp->permeability = factor * mp->porosity * pmv->r_pore * pmv->r_pore;

    /* Permeablity only depends on porosity currently */
    if (pd->v[pg->imtrx][POR_POROSITY]) {
      mp->d_permeability[POR_POROSITY] =
          factor * (mp->d_porosity[POR_POROSITY] * pmv->r_pore * pmv->r_pore +
                    2. * mp->porosity * pmv->r_pore * pmv->d_r_pore[POR_POROSITY]);
    }
  } else {
    GOMA_EH(GOMA_ERROR, "Unknown model for permeability");
  }

  return mp->permeability;
}
/* end of load_permeability */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void load_permeability_tensor(void)

/*************************************************************************
 *
 * load_permeability_tensor()
 *
 *  author:  P. R. Schunk.  April Fool's Day, 2006.
 *
 *  input:
 *  -----
 *    Possibly: mp->porosity, sink_mass
 *    Definitely:  K_11, K_22, Basis Vector components for rest state.
 *
 *  output:
 *  -----
 *    mp->perm_tensor[a][b]
 *    mp->d_perm_tensor[a][b][variable][j]
 *
 *************************************************************************/
{
  int aa, bb, cc, dd;
  int var, j;
  double K11 = 0.0, K22 = 0.0, K33 = 0.0, a[DIM][DIM];
  double a1_dot_F[DIM], a2_dot_F[DIM], a3_dot_F[DIM];
  double d_a1_dot_F_dx[DIM][DIM][MDE], d_a2_dot_F_dx[DIM][DIM][MDE], d_a3_dot_F_dx[DIM][DIM][MDE];

  memset(mp->perm_tensor, 0, sizeof(double) * DIM * DIM);
  memset(mp->d_perm_tensor, 0,
         sizeof(double) * DIM * DIM * (MAX_VARIABLE_TYPES + MAX_CONC + MAX_PMV));
  memset(a1_dot_F, 0, sizeof(double) * DIM);
  memset(a2_dot_F, 0, sizeof(double) * DIM);
  memset(a3_dot_F, 0, sizeof(double) * DIM);
  memset(d_a1_dot_F_dx, 0, sizeof(double) * DIM * DIM * MDE);
  memset(d_a2_dot_F_dx, 0, sizeof(double) * DIM * DIM * MDE);
  memset(d_a3_dot_F_dx, 0, sizeof(double) * DIM * DIM * MDE);

  /*FOR ORTHOTROPIC and other models
   *
   *                  Here we compute the local permeability tensor as
   *                    K = a_1.F.F_T.a_1 * K11 + a_2.F.F_T.a_2 * K22 + a_3.F.F_T.a_3 * K33
   *
   *                    K11, K22, K33 are the directional permeabilities in the product computed
   *                                  standard permeability models (imbedded here)
   *                    a_1 is the undeformed, initial-state vector parallel to the K11 Direction
   *                    a_2 is the undeformed, initial-state vector parallel to the K22 Direction
   *                    a_3 is the undeformed, initial-state vector parallel to the K33 Direction
   *                    F is the deformation gradient tensor.
   */

  if (mp->PermeabilityModel == ORTHOTROPIC) {

    /* preload some things */
    K11 = mp->u_permeability[0];
    K22 = mp->u_permeability[1];
    K33 = mp->u_permeability[2];
    /*  N.B.   there are 3 base vectors at the rest state that define the
     * orientation of the material relative to the orthotropic directions
     * If we have sheet laying on the x-y plane, then this matrix of vectors
     * is the Idemfactor.  This will most likely be the case.
     */

    a[0][0] = mp->u_permeability[3];
    a[0][1] = mp->u_permeability[4];
    a[0][2] = mp->u_permeability[5];

    a[1][0] = mp->u_permeability[6];
    a[1][1] = mp->u_permeability[7];
    a[1][2] = mp->u_permeability[8];

    a[2][0] = mp->u_permeability[9];
    a[2][1] = mp->u_permeability[10];
    a[2][2] = mp->u_permeability[11];

  } else if (mp->PermeabilityModel == SM_TENSOR) {
    GOMA_EH(GOMA_ERROR, "SM_TENSOR perm model is implemented but not yet tested");
    /*FOR particle SM_TENSOR (SinkMass_Tensor) - Porosity dependence on permeability with
     *mp->u_permeability[0] k_max in K11 direction (maximum saturated permeability or conductivity)
     *mp->u_permeability[1] k_rmin is the relative minumum permeability
     *mp->u_permeability[2] krscale is the relative scaled permeability
     *mp->u_permeability[3] k_exp is the powerlaw exponent
     *mp->u_permeability[4] K22 factor of K11/K_max  (0<=K22factor<=K11)
     *mp->u_permeability[5] K33 factor of K11/K_max

     * NB: For some reason this model NaNs real easily. Must initialize sink_mass
     * to nonzero value for one thing.   Need to investigate.
     */
    K11 = mp->u_permeability[0] *
          (mp->u_permeability[1] +
           (1. - mp->u_permeability[1]) /
               (1 + pow(fv->sink_mass * mp->u_permeability[2], mp->u_permeability[3])));
    K22 = K11 * mp->u_permeability[4];
    K33 = K11 * mp->u_permeability[5];
    /*  N.B.   there are 3 base vectors at the rest state that define the
     * orientation of the material relative to the orthotropic directions
     * If we have sheet laying on the x-y plane, then this matrix of vectors
     * is the Idemfactor.  This will most likely be the case.
     */

    a[0][0] = mp->u_permeability[6];
    a[0][1] = mp->u_permeability[7];
    a[0][2] = mp->u_permeability[8];

    a[1][0] = mp->u_permeability[9];
    a[1][1] = mp->u_permeability[10];
    a[1][2] = mp->u_permeability[11];

    a[2][0] = mp->u_permeability[12];
    a[2][1] = mp->u_permeability[13];
    a[2][2] = mp->u_permeability[14];

    /* Permeablity only depends on porosity currently */
    if (pd->v[pg->imtrx][POR_SINK_MASS]) {
      mp->d_perm_tensor[0][0][POR_SINK_MASS] =
          mp->u_permeability[0] * (mp->u_permeability[1] - 1.) * mp->u_permeability[3] *
          pow(fv->sink_mass * mp->u_permeability[2], mp->u_permeability[3] - 1.0) *
          mp->u_permeability[2] /
          (1 + pow(fv->sink_mass * mp->u_permeability[2], mp->u_permeability[3])) /
          (1 + pow(fv->sink_mass * mp->u_permeability[2], mp->u_permeability[3]));

      mp->d_perm_tensor[1][1][POR_SINK_MASS] =
          mp->d_perm_tensor[0][0][POR_SINK_MASS] * mp->u_permeability[4];
      mp->d_perm_tensor[2][2][POR_SINK_MASS] =
          mp->d_perm_tensor[0][0][POR_SINK_MASS] * mp->u_permeability[5];
    }

    /*  N.B.   PLEASE ADD SENSITIVIES TO REST OF CODE */
  } else if (mp->PermeabilityModel == KC_TENSOR) {
    /*FOR KOZENY_CARMAN - Porosity dependence on permeability with
     *mp->u_permeability[0] Kozeny-Carman const (=5 for particulate beds)
     *mp->u_permeability[1] surface area per unity volume
     *mp->u_permeability[2] K22 factor of K11/K_max
     *mp->u_permeability[3] K33 factor of K11/K_max
     */
    K11 = pow(mp->porosity, 3.0) / pow(1. - mp->porosity, 2.0) / mp->u_permeability[0] /
          pow(mp->u_permeability[1], 2.0);

    K22 = K11 * mp->u_permeability[2];

    if (pd->Num_Dim == 3)
      K33 = K11 * mp->u_permeability[3];

    a[0][0] = mp->u_permeability[4];
    a[0][1] = mp->u_permeability[5];
    a[0][2] = mp->u_permeability[6];

    a[1][0] = mp->u_permeability[7];
    a[1][1] = mp->u_permeability[8];
    a[1][2] = mp->u_permeability[9];

    a[2][0] = mp->u_permeability[10];
    a[2][1] = mp->u_permeability[11];
    a[2][2] = mp->u_permeability[12];

    /* Permeablity only depends on porosity currently */
    if (pd->v[pg->imtrx][POR_POROSITY]) {
      mp->d_perm_tensor[0][0][POR_POROSITY] =
          (pow(1. - mp->porosity, 2.0) * mp->u_permeability[0] * pow(mp->u_permeability[1], 2.0) *
               3.0 * pow(mp->porosity, 2.0) +
           pow(mp->porosity, 3.0) * 2.0 * mp->u_permeability[0] * pow(mp->u_permeability[1], 2.0) *
               (1.0 - mp->porosity)) /
          (pow(pow(1. - mp->porosity, 2.0) * mp->u_permeability[0] *
                   pow(mp->u_permeability[1], 2.0),
               2.0));

      mp->d_perm_tensor[1][1][POR_POROSITY] =
          mp->d_perm_tensor[0][0][POR_POROSITY] * mp->u_permeability[2];

      if (pd->Num_Dim == 3)
        mp->d_perm_tensor[2][2][POR_POROSITY] =
            mp->d_perm_tensor[0][0][POR_POROSITY] * mp->u_permeability[3];
    }
  } else if (mp->PermeabilityModel == EXTERNAL_FIELD) {

    /* preload some things */
    K11 = mp->permeability;
    K22 = mp->permeability;
    K33 = mp->permeability;
    /*  N.B.   there are 3 base vectors at the rest state that define the
     * orientation of the material relative to the orthotropic directions
     * If we have sheet laying on the x-y plane, then this matrix of vectors
     * is the Idemfactor.  This will most likely be the case.
     * PRS NOTE: for this field we are making it ISOTROPIC until need arises.
     */

    a[0][0] = 1.;
    a[0][1] = 0.;
    a[0][2] = 0.;

    a[1][0] = 0.;
    a[1][1] = 1.;
    a[1][2] = 0.;

    a[2][0] = 0.;
    a[2][1] = 0.;
    a[2][2] = 1.;

  } else {
    GOMA_EH(GOMA_ERROR, "Unknown model for tensor permeability");
  }

  /* If fixed mesh and not JAS decoupled */
  if (!pd->e[pg->imtrx][R_MESH1] && !efv->ev_porous_decouple) {
    for (aa = 0; aa < DIM; aa++) {
      for (bb = 0; bb < DIM; bb++) {
        mp->perm_tensor[aa][bb] =
            a[0][aa] * a[0][bb] * K11 + a[1][aa] * a[1][bb] * K22 + a[2][aa] * a[2][bb] * K33;
      }
    }
  } else if (pd->e[pg->imtrx][R_MESH1] && !efv->ev_porous_decouple) {
    for (aa = 0; aa < VIM; aa++) {
      for (bb = 0; bb < VIM; bb++) {
        for (cc = 0; cc < VIM; cc++) {
          a1_dot_F[cc] += a[0][aa] * fv->deform_grad[bb][cc] * delta(aa, bb);
          a2_dot_F[cc] += a[1][aa] * fv->deform_grad[bb][cc] * delta(aa, bb);
          a3_dot_F[cc] += a[2][aa] * fv->deform_grad[bb][cc] * delta(aa, bb);

          if (af->Assemble_Jacobian) {

            for (dd = 0; dd < VIM; dd++) {
              var = MESH_DISPLACEMENT1 + dd;
              for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                d_a1_dot_F_dx[cc][dd][j] +=
                    a[0][aa] * fv->d_deform_grad_dx[bb][cc][dd][j] * delta(aa, bb);
                d_a2_dot_F_dx[cc][dd][j] +=
                    a[1][aa] * fv->d_deform_grad_dx[bb][cc][dd][j] * delta(aa, bb);
                d_a3_dot_F_dx[cc][dd][j] +=
                    a[2][aa] * fv->d_deform_grad_dx[bb][cc][dd][j] * delta(aa, bb);
              }
            }
          }
        }
      }
    }
    for (aa = 0; aa < DIM; aa++) {
      for (bb = 0; bb < DIM; bb++) {
        mp->perm_tensor[aa][bb] = a1_dot_F[aa] * a1_dot_F[bb] * K11 +
                                  a2_dot_F[aa] * a2_dot_F[bb] * K22 +
                                  a3_dot_F[aa] * a3_dot_F[bb] * K33;
      }
    }

    if (af->Assemble_Jacobian) {
      for (aa = 0; aa < DIM; aa++) {
        for (bb = 0; bb < DIM; bb++) {
          for (dd = 0; dd < VIM; dd++) {
            var = MESH_DISPLACEMENT1 + dd;
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              mp->d_perm_tensor_dx[aa][bb][dd][j] = d_a1_dot_F_dx[aa][dd][j] * a1_dot_F[bb] * K11 +
                                                    a1_dot_F[aa] * d_a1_dot_F_dx[bb][dd][j] * K11 +
                                                    d_a2_dot_F_dx[aa][dd][j] * a2_dot_F[bb] * K22 +
                                                    a2_dot_F[aa] * d_a2_dot_F_dx[bb][dd][j] * K22 +
                                                    d_a3_dot_F_dx[aa][dd][j] * a3_dot_F[bb] * K33 +
                                                    a3_dot_F[aa] * d_a3_dot_F_dx[bb][dd][j] * K33;
            }
          }
        }
      }
    }

  }

  else if (efv->ev_porous_decouple) {

    GOMA_WH(GOMA_ERROR, "BIG WARNING! Any significant deformation will cast doubt on these results "
                        "as we have no Def Grad Tens for OrthoPerms");

    for (aa = 0; aa < DIM; aa++) {
      for (bb = 0; bb < DIM; bb++) {
        mp->perm_tensor[aa][bb] =
            a[0][aa] * a[0][bb] * K11 + a[1][aa] * a[1][bb] * K22 + a[2][aa] * a[2][bb] * K33;
      }
    }

  } else {
    GOMA_EH(GOMA_ERROR, "Not sure where to go from here dude: load_permeability_tensor");
  }
  /* Permeablity only depends on mesh motion,  currently */
  if (pd->v[pg->imtrx][MESH_DISPLACEMENT1]) {
    /* mp->d_perm_tensor[a][b][0][j] = 0.; */
  }

  return;
}
/* end of load_permeability_tensor */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void load_shell_permeability_tensor(void)

/*************************************************************************
 *
 * load_shell_permeability_tensor()
 *
 *  author:  S. A. Roberts. January 27, 2011
 *
 *  input:
 *  -----
 *    Material properties necessary for the calculation
 *
 *  output:
 *  -----
 *    mp->perm_tensor[a][b]
 *    mp->d_perm_tensor[a][b][variable][j]
 *
 *  note:
 *  -----
 *    As this routine is specific for shells, we will originally compile
 *    the tensor based on normal-tangent directions.  We will then need
 *    the shell normal to rotate it to the proper coordinates for the
 *    shell location.
 *
 *************************************************************************/
{

  /* Initialize */
  memset(mp->perm_tensor, 0, sizeof(double) * DIM * DIM);
  memset(mp->d_perm_tensor, 0,
         sizeof(double) * DIM * DIM * (MAX_VARIABLE_TYPES + MAX_CONC + MAX_PMV));
  memset(mp->d_perm_tensor_dx, 0, sizeof(double) * DIM * DIM * DIM * MDE);

  /* Global internal variables */
  int a, b, i, j, k;
  double NTperm[DIM][DIM] = {{0.0}};
  double ROT[DIM][DIM] = {{0.0}};
  double d_ROT_dx[DIM][DIM][DIM][MDE] = {{{{0.0}}}};

  /* === SHELL_CYLINDER_SQUARE ===
   * Represents a bundle of cylindrical tubes in a square packing.  The
   * cylinders are aligned in the direction normal to the shell.  Required
   * material properties are the cylinder radius and density.
   */
  if (mp->PermeabilityModel == SHELL_CYLINDER_SQUARE) {

    /* Load material properties */
    double r = 0.01; // Cylinder radius
    double e = 0.2;  // Cylinder density

    /* Parse properties */
    if (e > 0.3)
      GOMA_WH(GOMA_ERROR, "Your cylinder density is out of range.  Model only applicable for e < "
                          "0.3.  SHELL_CYLINDER_SQUARE");

    /* Calculate primary component to permeability */
    double k;
    k = log(1.0 / e) - 1.47633597;
    k += (2 * e - 0.79589781 * pow(e, 2)) / (1 + 0.48919241 * e - 1.60486942 * pow(e, 2));
    k *= M_PIE * pow(r, 2) / (8.0 * e);

    /* Compile permeability tensor in normal-tangent coordinates */
    NTperm[0][0] = k;
    NTperm[1][1] = k;
    NTperm[2][2] = 2 * k;

  } else {
    GOMA_EH(GOMA_ERROR, "Woah, pick another permeability model.  load_shell_permeability_tensor()");
  }

  /* Calculate rotation matrix */
  for (i = 0; i < DIM; i++) {
    ROT[0][i] = fv->stangent[0][i];
    ROT[1][i] = fv->stangent[1][i];
    ROT[2][i] = -fv->snormal[i];
    for (a = 0; a < DIM; a++) {
      for (b = 0; b < MDE; b++) {
        d_ROT_dx[0][i][a][b] = fv->dstangent_dx[0][i][a][b];
        d_ROT_dx[1][i][a][b] = fv->dstangent_dx[1][i][a][b];
        d_ROT_dx[2][i][a][b] = fv->dsnormal_dx[i][a][b];
      }
    }
  }

  /* Rotate normal-tangent permeability tensor */
  for (i = 0; i < DIM; i++) {
    for (j = 0; j < DIM; j++) {
      for (k = 0; k < DIM; k++) {
        mp->perm_tensor[i][j] += ROT[i][k] * NTperm[k][j];
        for (a = 0; a < DIM; a++) {
          for (b = 0; b < MDE; b++) {
            mp->d_perm_tensor_dx[i][j][a][b] = d_ROT_dx[i][k][a][b];
          }
        }
      }
    }
  }

  return;
}
/* end of load_shell_permeability_tensor */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

double load_saturation(double porosity, double cap_pres, double d_cap_pres[2])

/*************************************************************************
 * load_saturation -- calculate saturation in a porous
 *                    media from the unknowns used in the problem
 *
 *  input:   Assume that load_fv and load_fv_grads and load_fv_mesh_derivs
 *  -----    have all been called already, so unknowns and their
 *           sensitivies are known at this gauss point.
 *
 *  output:  calculates the saturation and its first and second
 *  -----    derivatives with respect to all the problem unknowns
 *           For transient calculations, the saturation at the
 *           old time step must be recalculated as well, in order
 *           to be used in the capacitance term.
 *
 *      mp->saturation
 *      mp->d_saturation
 *      mp->d_d_saturation
 *	    mp_old->saturation
 ************************************************************************/
{
  int i, j;
  double suction, saturation = 1.0e12, expon2;
  double suction_old = 1.0e12;
  double r_pore;
  const int i_pl = 0, i_pg = 1;
  double rad, d_sat_d_rad = 1e12, d_d_sat_d_rad = 1e12, frac;
  /* new variables for table option */
  double interp_val, var1[3], varold[3], slope;
  int var;
  /* new variables for tanh option */
  double con_a, con_b, con_c, con_d, cap_pres_clip;

  /*local variables for tanh_hyst option */
  double alpha_drain = 0.0, beta_drain = 0.0, s_min, sat_switch, pc_switch, pc_switch_clip;
  double alpha_wet = 0.0, beta_wet = 0.0, s_max;
  int ip, mat_ielem;
  extern int PRS_mat_ielem;

  /*
   *  Find which model is used for saturation and calculate:
   *    1)  The saturation, mp->saturation, mp_old->saturation
   *    2)  The sensitivity of saturation to all concentrations, porosity
   *        and temperature
   *          i.e., the first derivative of saturation w.r.t. each variable
   *           put this in mp->d_saturation[var]
   *    3)  The second derivatives of saturation w.r.t. each variable,
   *        including cross-terms, put this in mp->d_d_saturation[var][var]
   *        The second derivative is needed for sensitivity of fluxes
   *        which depend on the first derivative of saturation
   */
  if (mp->SaturationModel == VAN_GENUCHTEN) {

    /*
     * FOR VAN_GENUCHTEN EQUATION
     *  mp->u_saturation[0] is the irreduceable water saturation
     *  mp->u_saturation[1] is the irreduceable air saturation
     *  mp->u_saturation[2] is the exponent, beta
     *  mp->u_saturation[3] is the suction factor alpha/liquid density/gravity
     */
    expon2 = -(mp->u_saturation[2] - 1.) / mp->u_saturation[2];
    if (cap_pres <= 0.0) {
      saturation = mp->saturation = 1.0;
      mp->d_saturation[POR_LIQ_PRES] = 0.0;
      mp->d_d_saturation[POR_LIQ_PRES][POR_LIQ_PRES] = 0.0;
    } else {
      suction = cap_pres * mp->u_saturation[3];
      saturation = mp->saturation =
          mp->u_saturation[0] + (1.0 - mp->u_saturation[0] - mp->u_saturation[1]) *
                                    pow((1.0 + pow(suction, mp->u_saturation[2])), expon2);

      if (saturation > 1.0) {
        saturation = mp->saturation = 1.0;
      }
      if (saturation < 0.0) {
        saturation = mp->saturation = 0.0;
      }

      /*
       *  Calculate the first derivative of the saturation
       *  wrt to the dependent variables in the problem
       *  -> Note we need these terms even for a residual
       *     calculation
       *
       *   *** NOTE ---- this model does not currently have any
       *       sensitivity w.r.t. porosity or T
       *
       * ****Further NOTE----PRS added a negative here due to the redefinition
       * of capillary pressure to p_gas-p_liq for UNSATURATED case, from p_liq
       */
      mp->d_saturation[POR_LIQ_PRES] = -(1. - mp->u_saturation[0] - mp->u_saturation[1]) * expon2 *
                                       pow((1. + pow(suction, mp->u_saturation[2])), expon2 - 1.0) *
                                       mp->u_saturation[2] * mp->u_saturation[3] *
                                       pow(suction, mp->u_saturation[2] - 1);
      if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
        mp->d_saturation[POR_GAS_PRES] = -mp->d_saturation[POR_LIQ_PRES];
      }

      if (af->Assemble_Jacobian) {
        /*
         *  Calculate the second derivative of the dependence of saturation
         *  on the dependent variables in the problem.
         */
        /*n.g. this negative sign here is in question. Check out
          as I think it should be positive due to double derivative..*/
        mp->d_d_saturation[POR_LIQ_PRES][POR_LIQ_PRES] =
            -(1. - mp->u_saturation[0] - mp->u_saturation[1]) * expon2 * mp->u_saturation[2] *
            pow(mp->u_saturation[3], 2.0) *
            (pow((1.0 + pow(suction, mp->u_saturation[2])), expon2 - 2.0) *
                 pow(suction, 2.0 * mp->u_saturation[2] - 2.0) * (expon2 - 1.0) *
                 mp->u_saturation[2] +
             pow((1. + pow(suction, mp->u_saturation[2])), expon2 - 1.0) *
                 (mp->u_saturation[2] - 1.) * pow(suction, mp->u_saturation[2] - 2.0));
        if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
          mp->d_d_saturation[POR_GAS_PRES][POR_LIQ_PRES] =
              -mp->d_d_saturation[POR_GAS_PRES][POR_GAS_PRES];
          mp->d_d_saturation[POR_LIQ_PRES][POR_GAS_PRES] =
              mp->d_d_saturation[POR_GAS_PRES][POR_LIQ_PRES];
          mp->d_d_saturation[POR_GAS_PRES][POR_GAS_PRES] =
              mp->d_d_saturation[POR_LIQ_PRES][POR_LIQ_PRES];
        }
      }
    }

    /*
     * Now Calculate the old saturation value if we are doing a transient
     * calculation
     */
    if (pd->TimeIntegration == TRANSIENT) {
      if (pmv_old->cap_pres <= 0.0) {
        mp_old->saturation = 1.0;
        mp_old->d_saturation[POR_LIQ_PRES] = 0.0;
        mp_old->d_saturation[POR_GAS_PRES] = 0.0;
      } else {
        suction_old = pmv_old->cap_pres * mp->u_saturation[3];

        mp_old->saturation =
            mp->u_saturation[0] + (1. - mp->u_saturation[0] - mp->u_saturation[1]) *
                                      pow((1. + pow(suction_old, mp->u_saturation[2])), expon2);

        mp_old->d_saturation[POR_LIQ_PRES] =
            -(1.0 - mp->u_saturation[0] - mp->u_saturation[1]) * expon2 *
            pow((1.0 + pow(suction_old, mp->u_saturation[2])), expon2 - 1.0) * mp->u_saturation[2] *
            mp->u_saturation[3] * pow(suction_old, mp->u_saturation[2] - 1.0);
        if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
          mp_old->d_saturation[POR_GAS_PRES] = -mp_old->d_saturation[POR_LIQ_PRES];
        }
        if (pd->e[pg->imtrx][R_SHELL_SAT_OPEN]) {
          mp_old->d_saturation[SHELL_PRESS_OPEN] =
              -(1.0 - mp->u_saturation[0] - mp->u_saturation[1]) * expon2 *
              pow((1.0 + pow(suction_old, mp->u_saturation[2])), expon2 - 1.0) *
              mp->u_saturation[2] * mp->u_saturation[3] *
              pow(suction_old, mp->u_saturation[2] - 1.0);
        }
        if (pd->e[pg->imtrx][R_SHELL_SAT_OPEN_2]) {
          mp_old->d_saturation[SHELL_PRESS_OPEN_2] =
              -(1.0 - mp->u_saturation[0] - mp->u_saturation[1]) * expon2 *
              pow((1.0 + pow(suction_old, mp->u_saturation[2])), expon2 - 1.0) *
              mp->u_saturation[2] * mp->u_saturation[3] *
              pow(suction_old, mp->u_saturation[2] - 1.0);
        }
      }
    }

  } else if (mp->SaturationModel == PSD_VOL || mp->SaturationModel == PSD_WEXP ||
             mp->SaturationModel == PSD_SEXP) {
    /*
     * FOR PSD_ relationships
     *  mp->u_saturation[0] is the liquid-vapor surface tension
     *  mp->u_saturation[1] is the contact angle of menisci on the pore walls
     *  mp->u_saturation[2] is the Hammaker constant for adsorption
     *  mp->u_saturation[3] is the mean pore-size in undeformed state
     */

    /* first calculate the radius of the largest pores filled with
     * fluid (capillary pressure sets this radius)
     * put in check so that negative capillary pressures just result in a
     * saturated media
     */
    if (pmv->cap_pres > 1.0e-20) {
      pmv->r_cap = 2. * mp->u_saturation[0] * cos(mp->u_saturation[1]) / pmv->cap_pres;
    } else {
      pmv->r_cap = 2.0 * mp->u_saturation[0] * cos(mp->u_saturation[1]) / 1.0e-20;
    }

    for (i = 0; i < MAX_VARIABLE_TYPES + MAX_CONC; i++) {
      pmv->d_r_cap[i] = 0.0;
      pmv_old->d_r_cap[i] = 0.0;
      for (j = 0; j < MAX_VARIABLE_TYPES + MAX_CONC; j++) {
        pmv->d_d_r_cap[i][j] = 0.0;
      }
    }

    if (pmv->cap_pres > 1.0e-20) {
      pmv->d_r_cap[POR_LIQ_PRES] = -pmv->r_cap / pmv->cap_pres * d_cap_pres[i_pl];
      pmv->d_d_r_cap[POR_LIQ_PRES][POR_LIQ_PRES] =
          -2.0 * pmv->d_r_cap[POR_LIQ_PRES] / pmv->cap_pres * d_cap_pres[i_pl];

      if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
        pmv->d_r_cap[POR_GAS_PRES] = -pmv->r_cap / pmv->cap_pres * d_cap_pres[i_pg];
        pmv->d_d_r_cap[POR_GAS_PRES][POR_GAS_PRES] =
            -2. * pmv->d_r_cap[POR_GAS_PRES] / pmv->cap_pres * d_cap_pres[i_pg];
        pmv->d_d_r_cap[POR_GAS_PRES][POR_LIQ_PRES] =
            -2. * pmv->d_r_cap[POR_GAS_PRES] / pmv->cap_pres * d_cap_pres[i_pl];
        pmv->d_d_r_cap[POR_LIQ_PRES][POR_GAS_PRES] = pmv->d_d_r_cap[POR_GAS_PRES][POR_LIQ_PRES];
      }
    } else {
      pmv->d_r_cap[POR_LIQ_PRES] = 0.0;
      pmv->d_d_r_cap[POR_LIQ_PRES][POR_LIQ_PRES] = 0.0;
      if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
        pmv->d_r_cap[POR_GAS_PRES] = 0.0;
        pmv->d_d_r_cap[POR_GAS_PRES][POR_GAS_PRES] = 0.0;
        pmv->d_d_r_cap[POR_GAS_PRES][POR_LIQ_PRES] = 0.0;
        pmv->d_d_r_cap[POR_LIQ_PRES][POR_GAS_PRES] = 0.0;
      }
    }

    if (pmv->cap_pres > 1e-20 && mp->u_saturation[2] != 0) {
      /* check for adsorbed film (assume additive thickness) */
      if (mp->u_saturation[2] < 0) {
        GOMA_WH(GOMA_ERROR, "adsorption allowed with negative Hamaker constant, disjoining");
        pmv->r_cap += pow(-mp->u_saturation[2] / pmv->cap_pres, 1. / 3.);
      }
      if (mp->u_saturation[2] > 0) {
        GOMA_WH(GOMA_ERROR, "adsorption allowed with positive Hamaker constant, conjoining");
        pmv->r_cap -= pow(mp->u_saturation[2] / pmv->cap_pres, 1. / 3.);
        if (pmv->r_cap < 0)
          pmv->r_cap = 0;
      }

      /* check for adsorbed film (assume additive thickness) */
      if (mp->u_saturation[2] < 0) {
        pmv->d_r_cap[POR_GAS_PRES] += 1. / 3. * pow(-mp->u_saturation[2], 1. / 3.) /
                                      pow(pmv->cap_pres, 2. / 3.) * d_cap_pres[i_pl];
        if (pd->e[pg->imtrx][R_POR_GAS_PRES])
          pmv->d_r_cap[POR_GAS_PRES] += 1. / 3. * pow(-mp->u_saturation[2], 1. / 3.) /
                                        pow(pmv->cap_pres, 2. / 3.) * d_cap_pres[i_pg];
      }
      if (mp->u_saturation[2] > 0) {
        pmv->d_r_cap[POR_LIQ_PRES] -= 1. / 3. * pow(mp->u_saturation[2], 1. / 3.) /
                                      pow(pmv->cap_pres, 2. / 3.) * d_cap_pres[i_pl];
        if (pd->e[pg->imtrx][R_POR_GAS_PRES])
          pmv->d_r_cap[POR_GAS_PRES] -= 1. / 3. * pow(mp->u_saturation[2], 1. / 3.) /
                                        pow(pmv->cap_pres, 2. / 3.) * d_cap_pres[i_pg];
        if (pmv->r_cap == 0) {
          pmv->d_r_cap[POR_LIQ_PRES] = 0;
          if (pd->e[pg->imtrx][R_POR_GAS_PRES])
            pmv->d_r_cap[POR_GAS_PRES] = 0;
        }
      }

      /* check for adsorbed film (assume additive thickness) */
      if (mp->u_saturation[2] < 0) {
        pmv->d_d_r_cap[POR_LIQ_PRES][POR_LIQ_PRES] -=
            2. / 9. * pow(-mp->u_saturation[2], 1. / 3.) / pow(pmv->cap_pres, 5. / 3.);
        if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
          pmv->d_d_r_cap[POR_GAS_PRES][POR_GAS_PRES] -=
              2. / 9. * pow(-mp->u_saturation[2], 1. / 3.) / pow(pmv->cap_pres, 5. / 3.);
          pmv->d_d_r_cap[POR_GAS_PRES][POR_LIQ_PRES] -=
              2. / 9. * pow(-mp->u_saturation[2], 1. / 3.) / pow(pmv->cap_pres, 5. / 3.) *
              d_cap_pres[i_pg] * d_cap_pres[i_pl];
          pmv->d_d_r_cap[POR_LIQ_PRES][POR_GAS_PRES] = pmv->d_d_r_cap[POR_GAS_PRES][POR_LIQ_PRES];
        }
      }
      if (mp->u_saturation[2] > 0) {
        pmv->d_d_r_cap[POR_LIQ_PRES][POR_LIQ_PRES] -=
            2. / 9. * pow(mp->u_saturation[2], 1. / 3.) / pow(pmv->cap_pres, 5. / 3.);
        if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
          pmv->d_d_r_cap[POR_GAS_PRES][POR_GAS_PRES] -=
              2. / 9. * pow(mp->u_saturation[2], 1. / 3.) / pow(pmv->cap_pres, 5. / 3.);
          pmv->d_d_r_cap[POR_GAS_PRES][POR_LIQ_PRES] -=
              2. / 9. * pow(mp->u_saturation[2], 1. / 3.) / pow(pmv->cap_pres, 5. / 3.) *
              d_cap_pres[i_pg] * d_cap_pres[i_pl];
          pmv->d_d_r_cap[POR_LIQ_PRES][POR_GAS_PRES] = pmv->d_d_r_cap[POR_GAS_PRES][POR_LIQ_PRES];
        }
        if (pmv->r_cap == 0) {
          pmv->d_d_r_cap[POR_LIQ_PRES][POR_LIQ_PRES] = 0.;
          if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
            pmv->d_d_r_cap[POR_GAS_PRES][POR_GAS_PRES] = 0.;
            pmv->d_d_r_cap[POR_GAS_PRES][POR_LIQ_PRES] = 0.;
            pmv->d_d_r_cap[POR_LIQ_PRES][POR_GAS_PRES] = 0.;
          }
        }
      }
    }

    /* for transient, calculate old values */
    if (pd->TimeIntegration == TRANSIENT) {
      if (pmv_old->cap_pres > 1e-20) {
        pmv_old->r_cap = 2. * mp->u_saturation[0] * cos(mp->u_saturation[1]) / pmv_old->cap_pres;
      } else {
        pmv_old->r_cap = 2. * mp->u_saturation[0] * cos(mp->u_saturation[1]) / 1e-20;
      }
      if (pmv_old->cap_pres > 1e-20) {
        pmv_old->d_r_cap[POR_LIQ_PRES] = -pmv_old->r_cap / pmv_old->cap_pres * d_cap_pres[i_pl];
        if (pd->e[pg->imtrx][R_POR_GAS_PRES])
          pmv_old->d_r_cap[POR_GAS_PRES] = -pmv_old->r_cap / pmv_old->cap_pres * d_cap_pres[i_pg];
      } else {
        pmv_old->d_r_cap[POR_LIQ_PRES] = 0;
        if (pd->e[pg->imtrx][R_POR_GAS_PRES])
          pmv_old->d_r_cap[POR_GAS_PRES] = 0;
      }
    }

    if (mp->SaturationModel == PSD_VOL) {
      /* now find saturation, for which there are three regimes:
         1) r_cap > r_pore  => Saturated - all pores full
         2) alpha r_pore < r_cap < r_pore  => Partially saturated - some pores full
         3) r_cap > alpha r_pore => No Saturation - no pores full
         */
      r_pore = mp->u_saturation[3];
      if (pd->e[pg->imtrx][R_POR_POROSITY])
        r_pore = pmv->r_pore;
      if (pmv->r_cap > r_pore) {
        saturation = mp->saturation = 1.0;
      } else if (pmv->r_cap < mp->u_permeability[2] * r_pore) {
        saturation = mp->saturation = 0.0;
      } else {
        saturation = mp->saturation =
            ((pmv->r_cap / r_pore) - mp->u_permeability[2]) / (1. - mp->u_permeability[2]);

        mp->d_saturation[POR_LIQ_PRES] =
            pmv->d_r_cap[POR_LIQ_PRES] / r_pore / (1. - mp->u_permeability[2]);
        if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
          mp->d_saturation[POR_GAS_PRES] =
              pmv->d_r_cap[POR_GAS_PRES] / r_pore / (1. - mp->u_permeability[2]);
        }
        if (pd->e[pg->imtrx][R_POR_POROSITY]) {
          mp->d_saturation[POR_POROSITY] = -pmv->r_cap / r_pore / r_pore *
                                           pmv->d_r_pore[POR_POROSITY] /
                                           (1. - mp->u_permeability[2]);
        }

        if (af->Assemble_Jacobian) {
          mp->d_d_saturation[POR_LIQ_PRES][POR_LIQ_PRES] =
              pmv->d_d_r_cap[POR_LIQ_PRES][POR_LIQ_PRES] / r_pore / (1. - mp->u_permeability[2]);

          if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
            mp->d_d_saturation[POR_GAS_PRES][POR_GAS_PRES] =
                pmv->d_d_r_cap[POR_GAS_PRES][POR_GAS_PRES] / r_pore / (1. - mp->u_permeability[2]);
            mp->d_d_saturation[POR_GAS_PRES][POR_LIQ_PRES] =
                pmv->d_d_r_cap[POR_GAS_PRES][POR_LIQ_PRES] / r_pore / (1. - mp->u_permeability[2]);
            mp->d_d_saturation[POR_LIQ_PRES][POR_GAS_PRES] =
                pmv->d_d_r_cap[POR_LIQ_PRES][POR_GAS_PRES] / r_pore / (1. - mp->u_permeability[2]);
          }

          if (pd->e[pg->imtrx][R_POR_POROSITY]) {
            mp->d_d_saturation[POR_POROSITY][POR_POROSITY] =
                -(2 * mp->d_saturation[POR_POROSITY] / r_pore * pmv->d_r_pore[POR_POROSITY] +
                  pmv->r_cap / r_pore / r_pore * pmv->d_d_r_pore[POR_POROSITY][POR_POROSITY]) /
                (1. - mp->u_permeability[2]);
            mp->d_d_saturation[POR_POROSITY][POR_LIQ_PRES] = -pmv->d_r_cap[POR_LIQ_PRES] / r_pore /
                                                             r_pore * pmv->d_r_pore[POR_POROSITY] /
                                                             (1. - mp->u_permeability[2]);
            mp->d_d_saturation[POR_LIQ_PRES][POR_POROSITY] =
                mp->d_d_saturation[POR_POROSITY][POR_LIQ_PRES];

            if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
              mp->d_d_saturation[POR_POROSITY][POR_GAS_PRES] =
                  -pmv->d_r_cap[POR_GAS_PRES] / r_pore / r_pore * pmv->d_r_pore[POR_POROSITY] /
                  (1. - mp->u_permeability[2]);
              mp->d_d_saturation[POR_GAS_PRES][POR_POROSITY] =
                  mp->d_d_saturation[POR_POROSITY][POR_GAS_PRES];
            }
          }
        }
      }

      if (pd->TimeIntegration == TRANSIENT) {
        r_pore = mp->u_saturation[3];
        if (pd->e[pg->imtrx][R_POR_POROSITY])
          r_pore = pmv_old->r_pore;
        if (pmv_old->r_cap > r_pore) {
          mp_old->saturation = 1.0;
        } else if (pmv_old->r_cap < mp->u_permeability[2] * r_pore) {
          mp_old->saturation = 0.0;
        } else {
          mp_old->saturation =
              ((pmv_old->r_cap / r_pore) - mp->u_permeability[2]) / (1.0 - mp->u_permeability[2]);

          mp_old->d_saturation[POR_LIQ_PRES] =
              pmv_old->d_r_cap[POR_LIQ_PRES] / r_pore / (1.0 - mp->u_permeability[2]);
          if (pd->e[pg->imtrx][R_POR_GAS_PRES])
            mp_old->d_saturation[POR_GAS_PRES] =
                pmv_old->d_r_cap[POR_GAS_PRES] / r_pore / (1.0 - mp->u_permeability[2]);
          if (pd->e[pg->imtrx][R_POR_POROSITY])
            mp_old->d_saturation[POR_POROSITY] = -pmv_old->r_cap / r_pore / r_pore *
                                                 pmv_old->d_r_pore[POR_POROSITY] /
                                                 (1.0 - mp->u_permeability[2]);
        }
      }
    } else if (mp->SaturationModel == PSD_WEXP || mp->SaturationModel == PSD_SEXP) {
      r_pore = mp->u_saturation[3];
      if (pd->e[pg->imtrx][R_POR_POROSITY])
        r_pore = pmv->r_pore;
      rad = pmv->r_cap / r_pore;

      /* WEIGHTED EXPONENTIAL PORE_SIZE DISTRIBUTION */
      if (mp->SaturationModel == PSD_WEXP) {
        saturation = mp->saturation =
            1 - exp(-2. * rad) * (1. + 2. * rad + 2. * rad * rad + 4. / 3. * pow(rad, 3.0));

        d_sat_d_rad =
            (2. * exp(-2. * rad) * (1. + 2. * rad + 2 * rad * rad + 4. / 3. * pow(rad, 3.0)) -
             exp(-2. * rad) * (2. + 4 * rad + 4 * rad * rad));

        d_d_sat_d_rad = -2. * d_sat_d_rad + 2. * exp(-2. * rad) * (2. + 4. * rad + 4. * rad * rad) -
                        exp(-2. * rad) * (4. + 8. * rad);
      }

      /* SQUARED WEIGHTED EXPONENTIAL PORE_SIZE DISTRIBUTION - a narrower distribution*/
      if (mp->SaturationModel == PSD_SEXP) {
        frac = 9. * M_PIE / 16.;
        saturation = mp->saturation = 1 - exp(-frac * rad * rad) * (1. + frac * rad * rad);

        d_sat_d_rad = (2 * frac * rad * exp(-frac * rad * rad) * (1. + frac * rad * rad) -
                       exp(-frac * rad * rad) * (2 * frac * rad));

        d_d_sat_d_rad = -2 * frac * rad * d_sat_d_rad - exp(-frac * rad * rad) * (2 * frac) +
                        2 * frac * rad * exp(-frac * rad * rad) * (2 * frac * rad) +
                        2 * frac * exp(-frac * rad * rad) * (1. + frac * rad * rad);
      }

      mp->d_saturation[POR_LIQ_PRES] = pmv->d_r_cap[POR_LIQ_PRES] / r_pore * d_sat_d_rad;
      if (pd->e[pg->imtrx][R_POR_GAS_PRES])
        mp->d_saturation[POR_GAS_PRES] = -mp->d_saturation[POR_LIQ_PRES];
      if (pd->e[pg->imtrx][R_POR_POROSITY])
        mp->d_saturation[POR_POROSITY] = -rad / r_pore * pmv->d_r_pore[POR_POROSITY] * d_sat_d_rad;

      if (af->Assemble_Jacobian) {
        mp->d_d_saturation[POR_LIQ_PRES][POR_LIQ_PRES] =
            pmv->d_d_r_cap[POR_LIQ_PRES][POR_LIQ_PRES] / r_pore * d_sat_d_rad +
            pmv->d_r_cap[POR_LIQ_PRES] / r_pore * d_d_sat_d_rad * pmv->d_r_cap[POR_LIQ_PRES] /
                r_pore;

        if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
          mp->d_d_saturation[POR_GAS_PRES][POR_GAS_PRES] =
              mp->d_d_saturation[POR_LIQ_PRES][POR_LIQ_PRES];

          mp->d_d_saturation[POR_GAS_PRES][POR_LIQ_PRES] =
              pmv->d_d_r_cap[POR_GAS_PRES][POR_LIQ_PRES] / r_pore * d_sat_d_rad +
              pmv->d_r_cap[POR_GAS_PRES] / r_pore * d_d_sat_d_rad * pmv->d_r_cap[POR_LIQ_PRES] /
                  r_pore;

          mp->d_d_saturation[POR_LIQ_PRES][POR_GAS_PRES] =
              pmv->d_d_r_cap[POR_LIQ_PRES][POR_GAS_PRES] / r_pore * d_sat_d_rad +
              pmv->d_r_cap[POR_LIQ_PRES] / r_pore * d_d_sat_d_rad * pmv->d_r_cap[POR_GAS_PRES];
        }

        if (pd->e[pg->imtrx][R_POR_POROSITY]) {
          mp->d_d_saturation[POR_POROSITY][POR_POROSITY] =
              -2 * mp->d_saturation[POR_POROSITY] / r_pore * pmv->d_r_pore[POR_POROSITY] -
              rad / r_pore * d_sat_d_rad * pmv->d_d_r_pore[POR_POROSITY][POR_POROSITY] +
              rad / r_pore * rad / r_pore * d_d_sat_d_rad * pmv->d_r_pore[POR_POROSITY] *
                  pmv->d_r_pore[POR_POROSITY];
          mp->d_d_saturation[POR_POROSITY][POR_LIQ_PRES] =
              -pmv->d_r_cap[POR_LIQ_PRES] / r_pore / r_pore *
              (pmv->d_r_pore[POR_POROSITY] * d_sat_d_rad +
               rad * pmv->d_r_pore[POR_POROSITY] * d_d_sat_d_rad);
          mp->d_d_saturation[POR_LIQ_PRES][POR_POROSITY] =
              mp->d_d_saturation[POR_POROSITY][POR_LIQ_PRES];

          if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
            mp->d_d_saturation[POR_POROSITY][POR_GAS_PRES] =
                -pmv->d_r_cap[POR_GAS_PRES] / r_pore / r_pore *
                (pmv->d_r_pore[POR_POROSITY] * d_sat_d_rad +
                 rad * pmv->d_r_pore[POR_POROSITY] * d_d_sat_d_rad);
            mp->d_d_saturation[POR_GAS_PRES][POR_POROSITY] =
                mp->d_d_saturation[POR_POROSITY][POR_GAS_PRES];
          }
        }
      }

      if (pd->TimeIntegration == TRANSIENT) {
        r_pore = mp->u_saturation[3];
        if (pd->e[pg->imtrx][R_POR_POROSITY])
          r_pore = pmv_old->r_pore;
        rad = pmv_old->r_cap / r_pore;
        if (mp->SaturationModel == PSD_WEXP) {
          mp_old->saturation =
              1 - exp(-2. * rad) * (1. + 2. * rad + 2. * rad * rad + 4. / 3. * pow(rad, 3.0));

          d_sat_d_rad =
              (2. * exp(-2. * rad) * (1. + 2. * rad + 2. * rad * rad + 4 / 3. * pow(rad, 3.0)) -
               exp(-2. * rad) * (2. + 4. * rad + 4. * rad * rad));

          mp_old->d_saturation[POR_LIQ_PRES] =
              pmv_old->d_r_cap[POR_LIQ_PRES] / r_pore * d_sat_d_rad;
          if (pd->e[pg->imtrx][R_POR_GAS_PRES])
            mp_old->d_saturation[POR_GAS_PRES] = -mp_old->d_saturation[POR_LIQ_PRES];
          if (pd->e[pg->imtrx][R_POR_POROSITY])
            mp_old->d_saturation[POR_POROSITY] =
                -rad / r_pore * pmv_old->d_r_pore[POR_POROSITY] * d_sat_d_rad;
        }
        if (mp->SaturationModel == PSD_SEXP) {
          frac = 9. * M_PIE / 16.;
          mp_old->saturation = 1 - exp(-frac * rad * rad) * (1. + frac * rad * rad);

          d_sat_d_rad = (2 * frac * rad * exp(-frac * rad * rad) * (1. + frac * rad * rad) -
                         exp(-frac * rad * rad) * (2 * frac * rad));

          mp_old->d_saturation[POR_LIQ_PRES] =
              pmv_old->d_r_cap[POR_LIQ_PRES] / r_pore * d_sat_d_rad;
          if (pd->e[pg->imtrx][R_POR_GAS_PRES])
            mp_old->d_saturation[POR_GAS_PRES] = -mp_old->d_saturation[POR_LIQ_PRES];
          if (pd->e[pg->imtrx][R_POR_POROSITY])
            mp_old->d_saturation[POR_POROSITY] =
                -rad / r_pore * pmv_old->d_r_pore[POR_POROSITY] * d_sat_d_rad;
        }
      }
    }
  }
  /**********************************************************************
   *                   TANH MODEL FOR SATURATION
   **********************************************************************/
  else if (mp->SaturationModel == TANH) {
    /*
     * FOR TANH EQUATION
     *  mp->u_saturation[0] is the irreduceable water saturation
     *  mp->u_saturation[1] is the irreduceable air saturation
     *  mp->u_saturation[2] is variable 1
     *  mp->u_saturation[3] is variable 2
     */
    /* PRS Old 7/21/04:  con_a=0.5+mp->u_saturation[0]/2.0-(mp->u_saturation[1])/2.0; */
    /* PRS Old 7/21/04:  con_b=0.5-mp->u_saturation[0]/2.0-(mp->u_saturation[1])/2.0; */
    con_c = mp->u_saturation[2];
    con_d = mp->u_saturation[3];
    //    con_a=mp->u_saturation[0];
    //    con_b=mp->u_saturation[1];
    con_a = (mp->u_saturation[0] + tanh(con_c)) / (1. + tanh(con_c));
    con_b = (1. - mp->u_saturation[0]) / (1. + tanh(con_c));

    if (cap_pres <= 0.0) {
      cap_pres_clip = 1e-5;
    } else {
      cap_pres_clip = cap_pres;
    }
    saturation = mp->saturation = con_a - con_b * tanh(con_c - con_d / cap_pres_clip);
    /*
     *  Calculate the first derivative of the saturation
     *  wrt to the dependent variables in the problem
     *  -> Note we need these terms even for a residual
     *     calculation
     *
     *   *** NOTE ---- this model does not currently have any
     *       sensitivity w.r.t. porosity
     *
     * ****Further NOTE----PRS added a negative here due to the redefinition
     * of capillary pressure to p_gas-p_liq for UNSATURATED case, from p_liq
     * SO for this derivative it is positive
     */
    mp->d_saturation[POR_LIQ_PRES] = con_b * con_d *
                                     (1.0 - pow(tanh(con_c - con_d / cap_pres_clip), 2.0)) /
                                     cap_pres_clip / cap_pres_clip;
    if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
      mp->d_saturation[POR_GAS_PRES] = -mp->d_saturation[POR_LIQ_PRES];
    }

    if (pd->e[pg->imtrx][R_SHELL_SAT_OPEN]) {
      /*Ignore the above sensitivities and use this one because you cannot have both*/
      mp->d_saturation[SHELL_PRESS_OPEN] = con_b * con_d *
                                           (1.0 - pow(tanh(con_c - con_d / cap_pres_clip), 2.0)) /
                                           cap_pres_clip / cap_pres_clip;
    }

    if (pd->e[pg->imtrx][R_SHELL_SAT_OPEN_2]) {
      /*Ignore the above sensitivities and use this one because you cannot have both*/
      mp->d_saturation[SHELL_PRESS_OPEN_2] = con_b * con_d *
                                             (1.0 - pow(tanh(con_c - con_d / cap_pres_clip), 2.0)) /
                                             cap_pres_clip / cap_pres_clip;
    }

    if (af->Assemble_Jacobian) {
      /*
       *  Calculate the second derivative of the dependence of saturation
       *  on the dependent variables in the problem.
       */
      /*n.b. this negative sign here is in question. Check out
        as I think it should be positive due to double derivative..*/
      mp->d_d_saturation[POR_LIQ_PRES][POR_LIQ_PRES] =
          -(2.0 * con_b * con_d * con_d / pow(cap_pres_clip, 4.0) *
                tanh(con_c - con_d / cap_pres_clip) *
                (1.0 - pow(tanh(con_c - con_d / cap_pres_clip), 2.0)) +
            2.0 * con_b * con_d / pow(cap_pres_clip, 3.0) *
                (1.0 - pow(tanh(con_c - con_d / cap_pres_clip), 2.0)));
      if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
        mp->d_d_saturation[POR_GAS_PRES][POR_LIQ_PRES] =
            -mp->d_d_saturation[POR_GAS_PRES][POR_GAS_PRES];
        mp->d_d_saturation[POR_LIQ_PRES][POR_GAS_PRES] =
            mp->d_d_saturation[POR_GAS_PRES][POR_LIQ_PRES];
        mp->d_d_saturation[POR_GAS_PRES][POR_GAS_PRES] =
            mp->d_d_saturation[POR_LIQ_PRES][POR_LIQ_PRES];
      }

      if (pd->e[pg->imtrx][R_SHELL_SAT_OPEN]) {
        mp->d_d_saturation[SHELL_PRESS_OPEN][SHELL_PRESS_OPEN] =
            -(2.0 * con_b * con_d * con_d / pow(cap_pres_clip, 4.0) *
                  tanh(con_c - con_d / cap_pres_clip) *
                  (1.0 - pow(tanh(con_c - con_d / cap_pres_clip), 2.0)) +
              2.0 * con_b * con_d / pow(cap_pres_clip, 3.0) *
                  (1.0 - pow(tanh(con_c - con_d / cap_pres_clip), 2.0)));
      }

      if (pd->e[pg->imtrx][R_SHELL_SAT_OPEN_2]) {
        mp->d_d_saturation[SHELL_PRESS_OPEN_2][SHELL_PRESS_OPEN_2] =
            -(2.0 * con_b * con_d * con_d / pow(cap_pres_clip, 4.0) *
                  tanh(con_c - con_d / cap_pres_clip) *
                  (1.0 - pow(tanh(con_c - con_d / cap_pres_clip), 2.0)) +
              2.0 * con_b * con_d / pow(cap_pres_clip, 3.0) *
                  (1.0 - pow(tanh(con_c - con_d / cap_pres_clip), 2.0)));
      }
    }

    /*
     * Now Calculate the old saturation value if we are doing a transient
     * calculation
     */
    if (pd->TimeIntegration == TRANSIENT) {
      if (pmv_old->cap_pres <= 0.0) {
        cap_pres_clip = 1e-5;
      } else {
        cap_pres_clip = pmv_old->cap_pres;
      }
      mp_old->saturation = con_a - con_b * tanh(con_c - con_d / cap_pres_clip);
      mp_old->d_saturation[POR_LIQ_PRES] = con_b * con_d / cap_pres_clip / cap_pres_clip *
                                           (1 - pow(tanh(con_c - con_d / cap_pres_clip), 2.0));
      if (pd->e[pg->imtrx][R_SHELL_SAT_OPEN]) {
        mp_old->d_saturation[SHELL_PRESS_OPEN] =
            con_b * con_d / cap_pres_clip / cap_pres_clip *
            (1 - pow(tanh(con_c - con_d / cap_pres_clip), 2.0));
      }
      if (pd->e[pg->imtrx][R_SHELL_SAT_OPEN_2]) {
        mp_old->d_saturation[SHELL_PRESS_OPEN_2] =
            con_b * con_d / cap_pres_clip / cap_pres_clip *
            (1 - pow(tanh(con_c - con_d / cap_pres_clip), 2.0));
      }
      if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
        mp_old->d_saturation[POR_GAS_PRES] = -mp->d_saturation[POR_LIQ_PRES];
      }
    }
  } else if (mp->SaturationModel == TANH_EXTERNAL) {
    /* To scale these for with a 0-1 external field:
     * const = fv->external_field[efv->ev_SAT_index] * (mp->u_saturation[0] - mp->u_saturation[4]) -
     * mp->u_saturation[4])
     */
    /*
     * FOR TANH EQUATION
     *  mp->u_saturation[0] is the irreduceable water saturation
     *  mp->u_saturation[1] is the irreduceable air saturation
     *  mp->u_saturation[2] is variable 1
     *  mp->u_saturation[3] is variable 2
     */
    /* PRS Old 7/21/04:  con_a=0.5+mp->u_saturation[0]/2.0-(mp->u_saturation[1])/2.0; */
    /* PRS Old 7/21/04:  con_b=0.5-mp->u_saturation[0]/2.0-(mp->u_saturation[1])/2.0; */
    int i_sat_ev = mp->SAT_external_field_index;
    GOMA_EH(i_sat_ev, "Saturation external field not found!");
    con_c = fv->external_field[i_sat_ev] * (mp->u_saturation[2] - mp->u_saturation[6]) +
            mp->u_saturation[6];
    con_d = fv->external_field[i_sat_ev] * (mp->u_saturation[3] - mp->u_saturation[7]) +
            mp->u_saturation[7];
    con_a = (mp->u_saturation[0] + tanh(con_c)) / (1. + tanh(con_c));
    con_b = (1. - mp->u_saturation[0]) / (1. + tanh(con_c));

    if (cap_pres <= 0.0) {
      cap_pres_clip = 1e-5;
    } else {
      cap_pres_clip = cap_pres;
    }
    saturation = mp->saturation = con_a - con_b * tanh(con_c - con_d / cap_pres_clip);
    /*
     *  Calculate the first derivative of the saturation
     *  wrt to the dependent variables in the problem
     *  -> Note we need these terms even for a residual
     *     calculation
     *
     *   *** NOTE ---- this model does not currently have any
     *       sensitivity w.r.t. porosity
     *
     * ****Further NOTE----PRS added a negative here due to the redefinition
     * of capillary pressure to p_gas-p_liq for UNSATURATED case, from p_liq
     * SO for this derivative it is positive
     */
    mp->d_saturation[POR_LIQ_PRES] = con_b * con_d *
                                     (1.0 - pow(tanh(con_c - con_d / cap_pres_clip), 2.0)) /
                                     cap_pres_clip / cap_pres_clip;
    if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
      mp->d_saturation[POR_GAS_PRES] = -mp->d_saturation[POR_LIQ_PRES];
    }

    if (pd->e[pg->imtrx][R_SHELL_SAT_OPEN]) {
      /*Ignore the above sensitivities and use this one because you cannot have both*/
      mp->d_saturation[SHELL_PRESS_OPEN] = con_b * con_d *
                                           (1.0 - pow(tanh(con_c - con_d / cap_pres_clip), 2.0)) /
                                           cap_pres_clip / cap_pres_clip;
    }

    if (pd->e[pg->imtrx][R_SHELL_SAT_OPEN_2]) {
      /*Ignore the above sensitivities and use this one because you cannot have both*/
      mp->d_saturation[SHELL_PRESS_OPEN_2] = con_b * con_d *
                                             (1.0 - pow(tanh(con_c - con_d / cap_pres_clip), 2.0)) /
                                             cap_pres_clip / cap_pres_clip;
    }

    if (af->Assemble_Jacobian) {
      /*
       *  Calculate the second derivative of the dependence of saturation
       *  on the dependent variables in the problem.
       */
      /*n.b. this negative sign here is in question. Check out
        as I think it should be positive due to double derivative..*/
      mp->d_d_saturation[POR_LIQ_PRES][POR_LIQ_PRES] =
          -(2.0 * con_b * con_d * con_d / pow(cap_pres_clip, 4.0) *
                tanh(con_c - con_d / cap_pres_clip) *
                (1.0 - pow(tanh(con_c - con_d / cap_pres_clip), 2.0)) +
            2.0 * con_b * con_d / pow(cap_pres_clip, 3.0) *
                (1.0 - pow(tanh(con_c - con_d / cap_pres_clip), 2.0)));
      if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
        mp->d_d_saturation[POR_GAS_PRES][POR_LIQ_PRES] =
            -mp->d_d_saturation[POR_GAS_PRES][POR_GAS_PRES];
        mp->d_d_saturation[POR_LIQ_PRES][POR_GAS_PRES] =
            mp->d_d_saturation[POR_GAS_PRES][POR_LIQ_PRES];
        mp->d_d_saturation[POR_GAS_PRES][POR_GAS_PRES] =
            mp->d_d_saturation[POR_LIQ_PRES][POR_LIQ_PRES];
      }

      if (pd->e[pg->imtrx][R_SHELL_SAT_OPEN]) {
        mp->d_d_saturation[SHELL_PRESS_OPEN][SHELL_PRESS_OPEN] =
            -(2.0 * con_b * con_d * con_d / pow(cap_pres_clip, 4.0) *
                  tanh(con_c - con_d / cap_pres_clip) *
                  (1.0 - pow(tanh(con_c - con_d / cap_pres_clip), 2.0)) +
              2.0 * con_b * con_d / pow(cap_pres_clip, 3.0) *
                  (1.0 - pow(tanh(con_c - con_d / cap_pres_clip), 2.0)));
      }

      if (pd->e[pg->imtrx][R_SHELL_SAT_OPEN_2]) {
        mp->d_d_saturation[SHELL_PRESS_OPEN_2][SHELL_PRESS_OPEN_2] =
            -(2.0 * con_b * con_d * con_d / pow(cap_pres_clip, 4.0) *
                  tanh(con_c - con_d / cap_pres_clip) *
                  (1.0 - pow(tanh(con_c - con_d / cap_pres_clip), 2.0)) +
              2.0 * con_b * con_d / pow(cap_pres_clip, 3.0) *
                  (1.0 - pow(tanh(con_c - con_d / cap_pres_clip), 2.0)));
      }
    }

    /*
     * Now Calculate the old saturation value if we are doing a transient
     * calculation
     */
    if (pd->TimeIntegration == TRANSIENT) {
      if (pmv_old->cap_pres <= 0.0) {
        cap_pres_clip = 1e-5;
      } else {
        cap_pres_clip = pmv_old->cap_pres;
      }
      mp_old->saturation = con_a - con_b * tanh(con_c - con_d / cap_pres_clip);
      mp_old->d_saturation[POR_LIQ_PRES] = con_b * con_d / cap_pres_clip / cap_pres_clip *
                                           (1 - pow(tanh(con_c - con_d / cap_pres_clip), 2.0));
      if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
        mp_old->d_saturation[POR_GAS_PRES] = -mp->d_saturation[POR_LIQ_PRES];
      }
      if (pd->e[pg->imtrx][R_SHELL_SAT_OPEN]) {
        mp_old->d_saturation[SHELL_PRESS_OPEN] =
            con_b * con_d / cap_pres_clip / cap_pres_clip *
            (1 - pow(tanh(con_c - con_d / cap_pres_clip), 2.0));
      }
      if (pd->e[pg->imtrx][R_SHELL_SAT_OPEN_2]) {
        mp_old->d_saturation[SHELL_PRESS_OPEN_2] =
            con_b * con_d / cap_pres_clip / cap_pres_clip *
            (1 - pow(tanh(con_c - con_d / cap_pres_clip), 2.0));
      }
    }
  }
  /**********************************************************************
   *                   TANH_HYST MODEL FOR SATURATION
   **********************************************************************/
  else if (mp->SaturationModel == TANH_HYST) {
    /*
     * FOR TANH_HYST EQUATION
     *  For Initial or bounding Imbibe/wetting Curve
     *  mp->u_saturation[0] is the maximum water saturation (S_min)
     *  mp->u_saturation[1] is the irreduceable air saturation
     *  mp->u_saturation[2] is variable 1 (or beta_w in PRSs model)
     *  mp->u_saturation[3] is variable 2 (or alpha_w in PRSs model)
     *
     *  For DRAIN Curve
     *  mp->u_saturation[4] is the irreduceable water saturation
     *  mp->u_saturation[5] is the irreduceable air saturation
     *  mp->u_saturation[6] is variable 1 (or beta_d in PRSs model)
     *  mp->u_saturation[7] is variable 2 (or alpha_d in PRSs model)
     *
     *  Other input
     *  mp->u_saturation[8] Initial saturation curve for the material, viz.
     *                      if 1.0 then on the draining curve, and 0.0 then
     *                      on the wetting curve.
     *  mp->u_saturation[9] is Liquid_inventory_rate threshold for curve switching
     */

    /*
     *  Before we do anything, let us make sure we have the correct
     *  conditions for hysteretic type saturation functions.  If we are
     *  at the beginning of a time step, and the previous time step was
     *  successful, then test to see wether we are imbibing or draining
     *  and re-adjust the appropriate curve parameters accordingly.
     */

    ip = MMH_ip; /*someone should go change the name of this now oft-used
                  *expedient for gauss point datastructures. Does Matt
                  *really want his name attached to this?
                  */
    mat_ielem = PRS_mat_ielem;
    /*Hmm.  Now I am going to attached my name to another gem.
     *This is necessitated by the fact that I am too lazy to
     *bring Exo_DB struct all the way down here in the
     *dungeon. */

    if (cap_pres <= 0.0) {
      cap_pres_clip = 1e-5;
    } else {
      cap_pres_clip = cap_pres;
    }

    /*Load up con_a and con_b from element storage array */

    if (Element_Blocks[ei[pg->imtrx]->elem_blk_index].ElemStorage[mat_ielem].sat_curve_type[ip] ==
        1.0) {
      alpha_drain = mp->u_saturation[7];
      beta_drain = mp->u_saturation[6];
      s_min = mp->u_saturation[4];

      sat_switch =
          Element_Blocks[ei[pg->imtrx]->elem_blk_index].ElemStorage[mat_ielem].Sat_QP_tn[ip];
      pc_switch = Element_Blocks[ei[pg->imtrx]->elem_blk_index].ElemStorage[mat_ielem].p_cap_QP[ip];

      if (pc_switch <= 0.0) {
        pc_switch_clip = 1e-5;
      } else {
        pc_switch_clip = pc_switch;
      }
      con_b = (sat_switch - s_min) / (-1 + tanh(beta_drain - alpha_drain / pc_switch_clip));
      con_a = s_min - con_b;

      saturation = mp->saturation = con_a + con_b * tanh(beta_drain - alpha_drain / cap_pres_clip);

      if (saturation > 1.0) {
        /*	    fprintf(stderr,"sat way greater than unity, %lf %lf %lf %lf %lf", saturation,
           cap_pres_clip, s_max, con_a, con_b); GOMA_EH(GOMA_ERROR," Stopped! in tanh_hyst"); */
      }
      if (saturation > 1.0)
        saturation = mp->saturation = 1.0;
      if (saturation < 0.0) {
        fprintf(stderr, "Warning: Negative saturation detected on draining curve.\n");
        fprintf(stderr, "Proc: %d,  Element: %d, IP: %d, Saturation: %6.4g \n", ProcID,
                ei[pg->imtrx]->ielem + 1, ip, saturation);
      }
    } else /*on wetting curve */
    {
      alpha_wet = mp->u_saturation[3];
      beta_wet = mp->u_saturation[2];
      s_max = mp->u_saturation[0];

      sat_switch =
          Element_Blocks[ei[pg->imtrx]->elem_blk_index].ElemStorage[mat_ielem].Sat_QP_tn[ip];
      pc_switch = Element_Blocks[ei[pg->imtrx]->elem_blk_index].ElemStorage[mat_ielem].p_cap_QP[ip];

      if (pc_switch <= 0.0) {
        pc_switch_clip = 1e-5;
      } else {
        pc_switch_clip = pc_switch;
      }

      con_b = (sat_switch - s_max) / (tanh(beta_wet - alpha_wet / pc_switch_clip) + 1.0);
      con_a = 1. + con_b;

      saturation = mp->saturation = con_a + con_b * tanh(beta_wet - alpha_wet / cap_pres_clip);

      if (saturation > 1.0) {
        saturation = mp->saturation = 1.0;
        /*  fprintf(stderr,"sat way greater than unity, %lf %lf %lf %lf %lf", saturation,
           cap_pres_clip, s_max, con_a, con_b); GOMA_EH(GOMA_ERROR," Stopped! in tanh_hyst"); */
      }
      if (saturation < 0.0) {
        fprintf(stderr, "Warning: Negative saturation detected on draining curve.\n");
        fprintf(stderr, "Proc: %d,  Element: %d, IP: %d, Saturation: %6.4g \n", ProcID,
                ei[pg->imtrx]->ielem + 1, ip, saturation);
      }
    }

    /*
     *  Calculate the first derivative of the saturation
     *  wrt to the dependent variables in the problem
     *  -> Note we need these terms even for a residual
     *     calculation
     *
     */

    if (Element_Blocks[ei[pg->imtrx]->elem_blk_index].ElemStorage[mat_ielem].sat_curve_type[ip] ==
        1.0) {
      mp->d_saturation[POR_LIQ_PRES] =
          -con_b * alpha_drain / cap_pres_clip / cap_pres_clip *
          (1 - pow(tanh(beta_drain - alpha_drain / cap_pres_clip), 2.0));
      if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
        mp->d_saturation[POR_GAS_PRES] = -mp->d_saturation[POR_LIQ_PRES];
      }
    } else {
      mp->d_saturation[POR_LIQ_PRES] = -con_b * alpha_wet / cap_pres_clip / cap_pres_clip *
                                       (1 - pow(tanh(beta_wet - alpha_wet / cap_pres_clip), 2.0));
      if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
        mp->d_saturation[POR_GAS_PRES] = -mp->d_saturation[POR_LIQ_PRES];
      }
    }

    if (pd->e[pg->imtrx][R_SHELL_SAT_OPEN]) {
      if (Element_Blocks[ei[pg->imtrx]->elem_blk_index].ElemStorage[mat_ielem].sat_curve_type[ip] ==
          1.0) {
        mp->d_saturation[SHELL_PRESS_OPEN] =
            -con_b * alpha_drain / cap_pres_clip / cap_pres_clip *
            (1 - pow(tanh(beta_drain - alpha_drain / cap_pres_clip), 2.0));
      } else {
        mp->d_saturation[SHELL_PRESS_OPEN] =
            -con_b * alpha_wet / cap_pres_clip / cap_pres_clip *
            (1 - pow(tanh(beta_wet - alpha_wet / cap_pres_clip), 2.0));
      }
    }

    if (pd->e[pg->imtrx][R_SHELL_SAT_OPEN_2]) {
      if (Element_Blocks[ei[pg->imtrx]->elem_blk_index].ElemStorage[mat_ielem].sat_curve_type[ip] ==
          1.0) {
        mp->d_saturation[SHELL_PRESS_OPEN_2] =
            -con_b * alpha_drain / cap_pres_clip / cap_pres_clip *
            (1 - pow(tanh(beta_drain - alpha_drain / cap_pres_clip), 2.0));
      } else {
        mp->d_saturation[SHELL_PRESS_OPEN_2] =
            -con_b * alpha_wet / cap_pres_clip / cap_pres_clip *
            (1 - pow(tanh(beta_wet - alpha_wet / cap_pres_clip), 2.0));
      }
    }

    if (af->Assemble_Jacobian) {
      /*
       *  Calculate the second derivative of the dependence of saturation
       *  on the dependent variables in the problem.
       */
      /*n.g. this negative sign here is in question. Check out
        as I think it should be positive due to double derivative..*/
      if (Element_Blocks[ei[pg->imtrx]->elem_blk_index].ElemStorage[mat_ielem].sat_curve_type[ip] ==
          1.0) {
        mp->d_d_saturation[POR_LIQ_PRES][POR_LIQ_PRES] =
            (2.0 * con_b * alpha_drain * alpha_drain / pow(cap_pres_clip, 4.0) *
                 tanh(beta_drain - alpha_drain / cap_pres_clip) *
                 (1.0 - pow(tanh(beta_drain - alpha_drain / cap_pres_clip), 2.0)) +
             2.0 * con_b * alpha_drain / pow(cap_pres_clip, 3.0) *
                 (1.0 - pow(tanh(beta_drain - alpha_drain / cap_pres_clip), 2.0)));
      } else {
        mp->d_d_saturation[POR_LIQ_PRES][POR_LIQ_PRES] =
            +(2.0 * con_b * alpha_wet * alpha_wet / pow(cap_pres_clip, 4.0) *
                  tanh(beta_wet - alpha_wet / cap_pres_clip) *
                  (1.0 - pow(tanh(beta_wet - alpha_wet / cap_pres_clip), 2.0)) +
              2.0 * con_b * alpha_wet / pow(cap_pres_clip, 3.0) *
                  (1.0 - pow(tanh(beta_wet - alpha_wet / cap_pres_clip), 2.0)));
      }
      if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
        mp->d_d_saturation[POR_GAS_PRES][POR_LIQ_PRES] =
            -mp->d_d_saturation[POR_GAS_PRES][POR_GAS_PRES];
        mp->d_d_saturation[POR_LIQ_PRES][POR_GAS_PRES] =
            mp->d_d_saturation[POR_GAS_PRES][POR_LIQ_PRES];
        mp->d_d_saturation[POR_GAS_PRES][POR_GAS_PRES] =
            mp->d_d_saturation[POR_LIQ_PRES][POR_LIQ_PRES];
      }

      if (pd->e[pg->imtrx][R_SHELL_SAT_OPEN]) {
        if (Element_Blocks[ei[pg->imtrx]->elem_blk_index]
                .ElemStorage[mat_ielem]
                .sat_curve_type[ip] == 1.0) {
          mp->d_d_saturation[SHELL_PRESS_OPEN][SHELL_PRESS_OPEN] =
              (2.0 * con_b * alpha_drain * alpha_drain / pow(cap_pres_clip, 4.0) *
                   tanh(beta_drain - alpha_drain / cap_pres_clip) *
                   (1.0 - pow(tanh(beta_drain - alpha_drain / cap_pres_clip), 2.0)) +
               2.0 * con_b * alpha_drain / pow(cap_pres_clip, 3.0) *
                   (1.0 - pow(tanh(beta_drain - alpha_drain / cap_pres_clip), 2.0)));
        } else {
          mp->d_d_saturation[SHELL_PRESS_OPEN][SHELL_PRESS_OPEN] =
              +(2.0 * con_b * alpha_wet * alpha_wet / pow(cap_pres_clip, 4.0) *
                    tanh(beta_wet - alpha_wet / cap_pres_clip) *
                    (1.0 - pow(tanh(beta_wet - alpha_wet / cap_pres_clip), 2.0)) +
                2.0 * con_b * alpha_wet / pow(cap_pres_clip, 3.0) *
                    (1.0 - pow(tanh(beta_wet - alpha_wet / cap_pres_clip), 2.0)));
        }
      }

      if (pd->e[pg->imtrx][R_SHELL_SAT_OPEN_2]) {
        if (Element_Blocks[ei[pg->imtrx]->elem_blk_index]
                .ElemStorage[mat_ielem]
                .sat_curve_type[ip] == 1.0) {
          mp->d_d_saturation[SHELL_PRESS_OPEN_2][SHELL_PRESS_OPEN_2] =
              (2.0 * con_b * alpha_drain * alpha_drain / pow(cap_pres_clip, 4.0) *
                   tanh(beta_drain - alpha_drain / cap_pres_clip) *
                   (1.0 - pow(tanh(beta_drain - alpha_drain / cap_pres_clip), 2.0)) +
               2.0 * con_b * alpha_drain / pow(cap_pres_clip, 3.0) *
                   (1.0 - pow(tanh(beta_drain - alpha_drain / cap_pres_clip), 2.0)));
        } else {
          mp->d_d_saturation[SHELL_PRESS_OPEN_2][SHELL_PRESS_OPEN_2] =
              +(2.0 * con_b * alpha_wet * alpha_wet / pow(cap_pres_clip, 4.0) *
                    tanh(beta_wet - alpha_wet / cap_pres_clip) *
                    (1.0 - pow(tanh(beta_wet - alpha_wet / cap_pres_clip), 2.0)) +
                2.0 * con_b * alpha_wet / pow(cap_pres_clip, 3.0) *
                    (1.0 - pow(tanh(beta_wet - alpha_wet / cap_pres_clip), 2.0)));
        }
      }
    }

    /*
     * Now Calculate the old saturation value if we are doing a transient
     * calculation
     */
    if (pd->TimeIntegration == TRANSIENT) {
      if (pmv_old->cap_pres <= 0.0) {
        cap_pres_clip = 1e-5;
      } else {
        cap_pres_clip = pmv_old->cap_pres;
      }

      if (Element_Blocks[ei[pg->imtrx]->elem_blk_index].ElemStorage[mat_ielem].sat_curve_type[ip] ==
          1.0) {
        mp_old->saturation = con_a + con_b * tanh(beta_drain - alpha_drain / cap_pres_clip);
        mp_old->d_saturation[POR_LIQ_PRES] =
            con_b * alpha_drain / cap_pres_clip / cap_pres_clip *
            (1 - pow(tanh(beta_drain - alpha_drain / cap_pres_clip), 2.0));

        if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
          mp_old->d_saturation[POR_GAS_PRES] = -mp->d_saturation[POR_LIQ_PRES];
        }

        if (pd->e[pg->imtrx][R_SHELL_SAT_OPEN]) {
          mp_old->d_saturation[SHELL_PRESS_OPEN] = mp_old->d_saturation[POR_LIQ_PRES];
        }
        if (pd->e[pg->imtrx][R_SHELL_SAT_OPEN_2]) {
          mp_old->d_saturation[SHELL_PRESS_OPEN_2] = mp_old->d_saturation[POR_LIQ_PRES];
        }

      } else {
        mp_old->saturation = con_a + con_b * tanh(beta_wet - alpha_wet / cap_pres_clip);
        mp_old->d_saturation[POR_LIQ_PRES] =
            -con_b * alpha_wet / cap_pres_clip / cap_pres_clip *
            (1 - pow(tanh(beta_wet - alpha_wet / cap_pres_clip), 2.0));
        if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
          mp_old->d_saturation[POR_GAS_PRES] = -mp->d_saturation[POR_LIQ_PRES];
        }
        if (pd->e[pg->imtrx][R_SHELL_SAT_OPEN]) {
          mp_old->d_saturation[SHELL_PRESS_OPEN] = mp_old->d_saturation[POR_LIQ_PRES];
        }
        if (pd->e[pg->imtrx][R_SHELL_SAT_OPEN_2]) {
          mp_old->d_saturation[SHELL_PRESS_OPEN_2] = mp_old->d_saturation[POR_LIQ_PRES];
        }
      }
    }
  }
  /**********************************************************************
   *                   TABLE MODEL FOR SATURATION
   **********************************************************************/
  else if (mp->SaturationModel == TABLE) {
    struct Data_Table *table_local;
    /* next section is essentially apply_table_mp */
    table_local = MP_Tables[mp->saturation_tableid];
    for (i = 0; i < table_local->columns - 1; i++) {
      if (strcmp(table_local->t_name[i], "CAP_PRES") == 0) {
        if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
          var1[i] = cap_pres;
          varold[i] = pmv_old->cap_pres;
        } else {
          var1[i] = cap_pres;
          varold[i] = pmv_old->cap_pres;
        }
      } else {
        GOMA_EH(GOMA_ERROR, "Material Table Model Error-Unknown Function Column ");
      }
    }

    /*******************************************************************/
    /* Call appropriate interpolation schemes */
    interp_val = interpolate_table(table_local, var1, &slope, NULL);
    /*******************************************************************/
    saturation = mp->saturation = interp_val;

    /* Now compute sensitivities for each column */
    for (i = 0; i < table_local->columns - 1; i++) {
      var = table_local->t_index[i];

      switch (var) {
      case CAP_PRES:
        /*PRS: added sign change for POROUS_UNSAT case in rewrite */
        mp->d_saturation[POR_LIQ_PRES] = -table_local->slope[i];
        if (pd->e[pg->imtrx][R_POR_GAS_PRES])
          mp->d_saturation[POR_GAS_PRES] = -mp->d_saturation[POR_LIQ_PRES];

        if (pd->TimeIntegration == TRANSIENT) {
          /* ******************************************************************/
          /* Call appropriate interpolation schemes */
          interp_val = interpolate_table(table_local, varold, &slope, NULL);
          /* ******************************************************************/
          mp_old->saturation = interp_val;

          /* *** NOTE ---- this model does not currently have any
             sensitivity w.r.t. porosity */
          /*PRS: added sign change for POROUS_UNSAT case in rewrite */
          mp_old->d_saturation[POR_LIQ_PRES] = -table_local->slope[i];
          if (pd->e[pg->imtrx][R_POR_GAS_PRES])
            mp_old->d_saturation[POR_GAS_PRES] = -mp_old->d_saturation[POR_LIQ_PRES];
        }

        if (af->Assemble_Jacobian) {
          mp->d_d_saturation[POR_LIQ_PRES][POR_LIQ_PRES] = 0.0;
          if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
            mp->d_d_saturation[POR_GAS_PRES][POR_LIQ_PRES] =
                -mp->d_d_saturation[POR_GAS_PRES][POR_GAS_PRES];
            mp->d_d_saturation[POR_LIQ_PRES][POR_GAS_PRES] =
                mp->d_d_saturation[POR_GAS_PRES][POR_LIQ_PRES];
            mp->d_d_saturation[POR_GAS_PRES][POR_GAS_PRES] =
                mp->d_d_saturation[POR_LIQ_PRES][POR_LIQ_PRES];
          }
        }
        break;
      default:
        GOMA_EH(GOMA_ERROR, "Variable function not yet implemented in material property table");
      }
    }
  } else if (mp->SaturationModel == SHELL_TANH) {
    dbl dSdP, dSdP_P;
    saturation = mp->saturation = shell_saturation_pressure_curve(-cap_pres, &dSdP, &dSdP_P);
    if (pd->v[pg->imtrx][SHELL_PRESS_OPEN]) {
      mp->d_saturation[SHELL_PRESS_OPEN] = dSdP;
      mp->d_d_saturation[SHELL_PRESS_OPEN][SHELL_PRESS_OPEN] = dSdP_P;
      GOMA_WH(GOMA_ERROR,
              "SHELL_TANH MODEL not yet fitted for pressure gradients in both gas and liquid");
      // Do nothing because the saturation model is called from assemble_porous_shell_open
    }
    if (pd->v[pg->imtrx][SHELL_PRESS_OPEN_2]) {
      mp->d_saturation[SHELL_PRESS_OPEN_2] = dSdP;
      mp->d_d_saturation[SHELL_PRESS_OPEN_2][SHELL_PRESS_OPEN_2] = dSdP_P;
      GOMA_WH(GOMA_ERROR,
              "SHELL_TANH MODEL not yet fitted for pressure gradients in both gas and liquid");
      // Do nothing because the saturation model is called from assemble_porous_shell_open
    }

    /*
     * Now Calculate the old saturation value if we are doing a transient
     * calculation
     */
    if (pd->TimeIntegration == TRANSIENT) {
      cap_pres = pmv_old->cap_pres;
      mp_old->saturation = shell_saturation_pressure_curve(-cap_pres, &dSdP, &dSdP_P);

      if (pd->e[pg->imtrx][R_SHELL_SAT_OPEN]) {
        mp_old->d_saturation[SHELL_PRESS_OPEN] = dSdP;
      }
      if (pd->e[pg->imtrx][R_SHELL_SAT_OPEN_2]) {
        mp_old->d_saturation[SHELL_PRESS_OPEN_2] = dSdP;
      }
    }
  } else {
    GOMA_EH(GOMA_ERROR, "No models for saturation");
  }
  return saturation;
} /* end   load_saturation  */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

double load_cap_pres(int ipore, int ilnode, int ignode, double saturation)

/*************************************************************************
 * load_cap_pres -- calculate capillary pressure  in a porous
 *                    media from the unknowns used in the problem
 *
 *  input:   Assume that load_fv and load_fv_grads and load_fv_mesh_derivs
 *  -----    have all been called already, so unknowns and their
 *           sensitivies are known at this gauss point.
 *
 *  output:  calculates the saturation and its first and second
 *  -----    derivatives with respect to all the problem unknowns
 *           For transient calculations, the saturation at the
 *           old time step must be recalculated as well, in order
 *           to be used in the capacitance term.
 *
 *      mp->cap_pres
 *      mp->d_cap_pres
 *      mp->d_d_cap_pres
 *	    mp_old->cap_pres
 ************************************************************************/
{
  double con_a, con_b, con_c, con_d;
  double sat_min, sat_max;
  double sat_norm, sat_clip, d_sat_norm_d_saturation, d_sat_clip_d_saturation;
  double cap_pres = 0.0;
  double mexp, nexp, n_inv, m_inv;
  double brack_in, d_brack_in_d_sat_norm, d_d_brack_in_d_sat_norm;

  int i_ext_field;
  double val_ext_field;
  double sat_min_1, sat_max_1, sat_min_2, sat_max_2;
  double con_c_1, nexp_1, con_c_2, nexp_2;

  int draining_curve, draining_curve_old;
  int switch_now;
  double cap_pres_switch, sat_switch, sat_norm_switch;
  double brack_in_inv, sat_min_wet, sat_max_dry;

  /*
   *  Find which model is used for saturation and calculate:
   *    1)  The capillary pressure, mp->cap_pres,
   *    2)  The sensitivity of saturation to all concentrations, porosity
   *        and temperature
   *          i.e., the first derivative of saturation w.r.t. each variable
   *           put this in mp->d_cap_pres[var]
   *    3)  The second derivatives of saturation w.r.t. each variable,
   *        including cross-terms, put this in mp->d_d_cap_pres[var][var]
   *        The second derivative is needed for sensitivity of fluxes
   *        which depend on the first derivative of capillary pressure
   */

  /**********************************************************************
   *                   ATANH MODEL FOR CAPILLARY PRESSURE
   **********************************************************************/
  if (mp->PorousShellCapPresModel[ipore] == ATANH) {
    /*
     * FOR ATANH EQUATION
     *  mp->u_saturation[0] is the irreduceable water saturation
     *  mp->u_saturation[1] is the irreduceable air saturation
     *  mp->u_saturation[2] is shift factor
     *  mp->u_saturation[3] is multiplier in tanh function
     */
    //    sat_min = mp->u_cap_pres[0];
    //    sat_max = mp->u_cap_pres[1];
    //    con_c = mp->u_cap_pres[2];
    //    con_d = mp->u_cap_pres[3];

    /* Normalized saturation - should be between -1 to +1 */
    //    sat_norm = (2.0 * saturation - sat_max - sat_min) / (sat_max - sat_min);
    //    d_sat_norm_d_saturation = 2.0/(sat_max - sat_min);

    con_a = mp->u_PorousShellCapPres[ipore][0];
    con_b = mp->u_PorousShellCapPres[ipore][1];
    con_c = mp->u_PorousShellCapPres[ipore][2];
    con_d = mp->u_PorousShellCapPres[ipore][3];

    sat_norm = (con_a - saturation) / con_b;
    d_sat_norm_d_saturation = -1.0 / con_b;

    /* Clip if necessary */
    if (sat_norm <= -0.995) {
      sat_clip = -0.995;
      d_sat_clip_d_saturation = 0.0;
    } else if (sat_norm >= 0.995) {
      sat_clip = 0.995;
      d_sat_clip_d_saturation = 0.0;
    } else {
      sat_clip = sat_norm;
      d_sat_clip_d_saturation = d_sat_norm_d_saturation;
    }
    //    cap_pres = mp->cap_pres = con_c - con_d * atanh(sat_clip);
    cap_pres = mp->cap_pres = con_d / (con_c - atanh(sat_clip));

    /*
     *  Calculate the first derivative of the capillary pressure
     *  wrt to the dependent variables in the problem
     *  -> Note we need these terms even for a residual
     *     calculation
     *
     *   *** NOTE ---- this model as it is only works for porous shell
     *                 saturation formulations
     *
     */
    //    mp->d_cap_pres[SHELL_SAT_1] = -con_d/(1.0 - sat_clip * sat_clip) *
    //    d_sat_clip_d_saturation;

    switch (ipore) {
    case 0:
      mp->d_cap_pres[SHELL_SAT_1] = con_d /
                                    ((con_c - atanh(sat_clip)) * (con_c - atanh(sat_clip))) /
                                    (1.0 - (sat_clip * sat_clip)) * d_sat_clip_d_saturation;
      break;

    case 1:
      mp->d_cap_pres[SHELL_SAT_2] = con_d /
                                    ((con_c - atanh(sat_clip)) * (con_c - atanh(sat_clip))) /
                                    (1.0 - (sat_clip * sat_clip)) * d_sat_clip_d_saturation;
      break;

    case 2:
      mp->d_cap_pres[SHELL_SAT_3] = con_d /
                                    ((con_c - atanh(sat_clip)) * (con_c - atanh(sat_clip))) /
                                    (1.0 - (sat_clip * sat_clip)) * d_sat_clip_d_saturation;
      break;

    default:
      GOMA_EH(GOMA_ERROR, "Not valid porous shell index");
      break;
    }
  } else if (mp->PorousShellCapPresModel[ipore] == SINH) {
    /*
     * FOR SINH EQUATION
     *  mp->u_saturation[0] is the irreducable water saturation
     *  mp->u_saturation[1] is the irreduceable air saturation
     *  mp->u_saturation[2] is shift factor
     *  mp->u_saturation[3] is multiplier in tanh function
     */

    sat_min = mp->u_PorousShellCapPres[ipore][0];
    sat_max = mp->u_PorousShellCapPres[ipore][1];
    con_c = mp->u_PorousShellCapPres[ipore][2];
    con_d = mp->u_PorousShellCapPres[ipore][3];

    /* Normalized saturation - should be between -6 to +6 */
    sat_norm = (12.0 * saturation - 6.0 * (sat_max + sat_min)) / (sat_max - sat_min);
    d_sat_norm_d_saturation = 12.0 / (sat_max - sat_min);

    /* Clip if necessary */
    //    if (saturation <= 0.0) {
    //      sat_clip = -6.0;
    //      d_sat_clip_d_saturation = 0.0;
    //    }else if (saturation >= 0.995) {
    //      sat_clip = 6.0;
    //      d_sat_clip_d_saturation = 0.0;
    //    } else {
    sat_clip = sat_norm;
    d_sat_clip_d_saturation = d_sat_norm_d_saturation;
    //    }
    cap_pres = con_c - con_d * sinh(sat_clip);

    /*
     *  Calculate the first derivative of the capillary pressure
     *  wrt to the dependent variables in the problem
     *  -> Note we need these terms even for a residual
     *     calculation
     *
     *   *** NOTE ---- this model as it is only works for porous shell
     *                 saturation formulations
     *
     */

    switch (ipore) {
    case 0:
      mp->d_cap_pres[SHELL_SAT_1] = -con_d * cosh(sat_clip) * d_sat_clip_d_saturation;
      break;

    case 1:
      mp->d_cap_pres[SHELL_SAT_2] = -con_d * cosh(sat_clip) * d_sat_clip_d_saturation;
      break;

    case 2:
      mp->d_cap_pres[SHELL_SAT_3] = -con_d * cosh(sat_clip) * d_sat_clip_d_saturation;
      break;

    default:
      GOMA_EH(GOMA_ERROR, "Not valid porous shell index");
      break;
    }
  } else if (mp->PorousShellCapPresModel[ipore] == VAN_GENUCHTEN) {
    /*
     * FOR VAN_GENUCHTEN EQUATION
     *  mp->u_saturation[0] is the irreducable water saturation
     *  mp->u_saturation[1] is the irreduceable air saturation
     *  mp->u_saturation[2] is entry pressure
     *  mp->u_saturation[3] is exponent n
     */

    sat_min = mp->u_PorousShellCapPres[ipore][0];
    sat_max = mp->u_PorousShellCapPres[ipore][1];
    con_c = mp->u_PorousShellCapPres[ipore][2];
    nexp = mp->u_PorousShellCapPres[ipore][3];

    n_inv = 1.0 / nexp;
    mexp = 1.0 - n_inv;
    m_inv = -1.0 / mexp;

    /* Normalized saturation - should be between 0 to 1 */
    sat_norm = (saturation - sat_min) / (sat_max - sat_min);
    d_sat_norm_d_saturation = 1.0 / (sat_max - sat_min);

    brack_in = pow(sat_norm, m_inv) - 1.0;
    d_brack_in_d_sat_norm = m_inv * pow(sat_norm, (m_inv - 1.0));
    d_d_brack_in_d_sat_norm = m_inv * (m_inv - 1.0) * pow(sat_norm, (m_inv - 2.0));

    cap_pres = con_c * pow(brack_in, n_inv);

    /*
     *  Calculate the first derivative of the capillary pressure
     *  wrt to the dependent variables in the problem
     *  -> Note we need these terms even for a residual
     *     calculation
     *
     *   *** NOTE ---- this model as it is only works for porous shell
     *                 saturation formulations
     *
     */

    switch (ipore) {
    case 0:
      mp->d_cap_pres[SHELL_SAT_1] = n_inv * con_c * pow(brack_in, (n_inv - 1.0)) *
                                    d_brack_in_d_sat_norm * d_sat_norm_d_saturation;
      mp->d_d_cap_pres[SHELL_SAT_1][SHELL_SAT_1] =
          n_inv * (n_inv - 1.0) * con_c * pow(brack_in, (n_inv - 2.0)) * d_brack_in_d_sat_norm *
              d_brack_in_d_sat_norm * d_sat_norm_d_saturation * d_sat_norm_d_saturation +
          n_inv * con_c * pow(brack_in, (n_inv - 1.0)) * d_d_brack_in_d_sat_norm *
              d_sat_norm_d_saturation * d_sat_norm_d_saturation;
      break;

    case 1:
      mp->d_cap_pres[SHELL_SAT_2] = n_inv * con_c * pow(brack_in, (n_inv - 1.0)) *
                                    d_brack_in_d_sat_norm * d_sat_norm_d_saturation;
      mp->d_d_cap_pres[SHELL_SAT_2][SHELL_SAT_2] =
          n_inv * (n_inv - 1.0) * con_c * pow(brack_in, (n_inv - 2.0)) * d_brack_in_d_sat_norm *
              d_brack_in_d_sat_norm * d_sat_norm_d_saturation * d_sat_norm_d_saturation +
          n_inv * con_c * pow(brack_in, (n_inv - 1.0)) * d_d_brack_in_d_sat_norm *
              d_sat_norm_d_saturation * d_sat_norm_d_saturation;
      break;

    case 2:
      mp->d_cap_pres[SHELL_SAT_3] = n_inv * con_c * pow(brack_in, (n_inv - 1.0)) *
                                    d_brack_in_d_sat_norm * d_sat_norm_d_saturation;
      mp->d_d_cap_pres[SHELL_SAT_3][SHELL_SAT_3] =
          n_inv * (n_inv - 1.0) * con_c * pow(brack_in, (n_inv - 2.0)) * d_brack_in_d_sat_norm *
              d_brack_in_d_sat_norm * d_sat_norm_d_saturation * d_sat_norm_d_saturation +
          n_inv * con_c * pow(brack_in, (n_inv - 1.0)) * d_d_brack_in_d_sat_norm *
              d_sat_norm_d_saturation * d_sat_norm_d_saturation;

      break;

    default:
      GOMA_EH(GOMA_ERROR, "Not valid porous shell index");
      break;
    }

  } else if (mp->PorousShellCapPresModel[ipore] == VAN_GENUCHTEN_EXTERNAL) {
    /*
     * FOR VAN_GENUCHTEN EXTERNAL EQUATION
     *  mp->u_saturation[0] is the irreducable water saturation for curve 1
     *  mp->u_saturation[1] is the irreduceable air saturation for curve 1
     *  mp->u_saturation[2] is entry pressure for curve 1
     *  mp->u_saturation[3] is exponent n for curve 1
     *
     *  mp->u_saturation[4] is the irreducable water saturation for curve 2
     *  mp->u_saturation[5] is the irreduceable air saturation for curve 2
     *  mp->u_saturation[6] is entry pressure for curve 2
     *  mp->u_saturation[7] is exponent n for curve 2
     *
     */

    i_ext_field = mp->por_shell_cap_pres_ext_field_index[ipore];
    /* Here I assume the external field value ranges from 0 to 1 */
    val_ext_field = *evp->external_field[i_ext_field][ilnode];

    sat_min_1 = mp->u_PorousShellCapPres[ipore][0];
    sat_max_1 = mp->u_PorousShellCapPres[ipore][1];
    con_c_1 = mp->u_PorousShellCapPres[ipore][2];
    nexp_1 = mp->u_PorousShellCapPres[ipore][3];

    sat_min_2 = mp->u_PorousShellCapPres[ipore][4];
    sat_max_2 = mp->u_PorousShellCapPres[ipore][5];
    con_c_2 = mp->u_PorousShellCapPres[ipore][6];
    nexp_2 = mp->u_PorousShellCapPres[ipore][7];

    /* Interpolate the fitting parameters with external field value */
    sat_min = sat_min_1 + val_ext_field * (sat_min_2 - sat_min_1);
    sat_max = sat_max_1 + val_ext_field * (sat_max_2 - sat_max_1);
    con_c = con_c_1 + val_ext_field * (con_c_2 - con_c_1);
    nexp = nexp_1 + val_ext_field * (nexp_2 - nexp_1);

    n_inv = 1.0 / nexp;
    mexp = 1.0 - n_inv;
    m_inv = -1.0 / mexp;

    /* Normalized saturation - should be between 0 to 1 */
    sat_norm = (saturation - sat_min) / (sat_max - sat_min);
    d_sat_norm_d_saturation = 1.0 / (sat_max - sat_min);

    brack_in = pow(sat_norm, m_inv) - 1.0;
    d_brack_in_d_sat_norm = m_inv * pow(sat_norm, (m_inv - 1.0));

    cap_pres = con_c * pow(brack_in, n_inv);

    /*
     *  Calculate the first derivative of the capillary pressure
     *  wrt to the dependent variables in the problem
     *  -> Note we need these terms even for a residual
     *     calculation
     *
     *   *** NOTE ---- this model as it is only works for porous shell
     *                 saturation formulations
     *
     */

    switch (ipore) {
    case 0:
      mp->d_cap_pres[SHELL_SAT_1] = n_inv * con_c * pow(brack_in, (n_inv - 1.0)) *
                                    d_brack_in_d_sat_norm * d_sat_norm_d_saturation;
      break;

    case 1:
      mp->d_cap_pres[SHELL_SAT_2] = n_inv * con_c * pow(brack_in, (n_inv - 1.0)) *
                                    d_brack_in_d_sat_norm * d_sat_norm_d_saturation;
      break;

    case 2:
      mp->d_cap_pres[SHELL_SAT_3] = n_inv * con_c * pow(brack_in, (n_inv - 1.0)) *
                                    d_brack_in_d_sat_norm * d_sat_norm_d_saturation;
      break;

    default:
      GOMA_EH(GOMA_ERROR, "Not valid porous shell index");
      break;
    }
  } else if (mp->PorousShellCapPresModel[ipore] == VAN_GENUCHTEN_HYST) {
    /*
     * FOR VAN_GENUCHTEN_HYST EQUATION
     *  mp->u_saturation[0] is the irreducable water saturation for wetting curve
     *  mp->u_saturation[1] is the irreduceable air saturation for wetting curve
     *  mp->u_saturation[2] is entry pressure for wetting curve
     *  mp->u_saturation[3] is exponent n for wetting curve
     *
     *  mp->u_saturation[4] is the irreducable water saturation for draining curve
     *  mp->u_saturation[5] is the irreduceable air saturation for draining curve
     *  mp->u_saturation[6] is entry pressure for draining curve
     *  mp->u_saturation[7] is exponent n for draining curve
     *
     *  Other input
     *  mp->u_saturation[8] Initial saturation curve for the material, viz.
     *                      if 1.0 then on the draining curve, and 0.0 then
     *                      on the wetting curve.
     *  mp->u_saturation[9] is Liquid_inventory_rate threshold for curve switching
     */

    /* Quality check of global node number ignode*/
    if (ignode < 0) {
      if (ipore == 0) {
        ignode = ei[pg->imtrx]->gnn_list[SHELL_SAT_1][ilnode];
      } else if (ipore == 1) {
        ignode = ei[pg->imtrx]->gnn_list[SHELL_SAT_2][ilnode];
      } else if (ipore == 2) {
        ignode = ei[pg->imtrx]->gnn_list[SHELL_SAT_3][ilnode];
      } else {
        GOMA_EH(GOMA_ERROR, "Something's wrong here!");
      }
    }

    /* First, see if we are switching */
    switch_now = pmv_hyst->curve_switch[ipore][ignode];

    /* Get the minimum and maximum saturation for the scanning curves */
    sat_min_wet = pmv_hyst->sat_min_imbibe[ipore][ignode];
    sat_max_dry = pmv_hyst->sat_max_drain[ipore][ignode];

    /* Find out whether on draining or wetting curve */
    draining_curve = pmv_hyst->curve_type[ipore][ignode];
    draining_curve_old = pmv_hyst->curve_type_old[ipore][ignode];

    /* If we are not switching */
    if (switch_now == 0) {
      if (draining_curve == 1) /* Stay on main draining curve */
      {
        sat_min = mp->u_PorousShellCapPres[ipore][4];
        sat_max = sat_max_dry;
        con_c = mp->u_PorousShellCapPres[ipore][6];
        nexp = mp->u_PorousShellCapPres[ipore][7];

        n_inv = 1.0 / nexp;
        mexp = 1.0 - n_inv;
        m_inv = -1.0 / mexp;

        /* Normalized saturation - should be between 0 to 1 */
        sat_norm = (saturation - sat_min) / (sat_max - sat_min);
        d_sat_norm_d_saturation = 1.0 / (sat_max - sat_min);

        brack_in = pow(sat_norm, m_inv) - 1.0;
        d_brack_in_d_sat_norm = m_inv * pow(sat_norm, (m_inv - 1.0));

        /* Evaluate capillary pressure and its sensitivities */
        cap_pres = con_c * pow(brack_in, n_inv);

        switch (ipore) {
        case 0:
          mp->d_cap_pres[SHELL_SAT_1] = n_inv * con_c * pow(brack_in, (n_inv - 1.0)) *
                                        d_brack_in_d_sat_norm * d_sat_norm_d_saturation;
          break;

        case 1:
          mp->d_cap_pres[SHELL_SAT_2] = n_inv * con_c * pow(brack_in, (n_inv - 1.0)) *
                                        d_brack_in_d_sat_norm * d_sat_norm_d_saturation;
          break;

        case 2:
          mp->d_cap_pres[SHELL_SAT_3] = n_inv * con_c * pow(brack_in, (n_inv - 1.0)) *
                                        d_brack_in_d_sat_norm * d_sat_norm_d_saturation;
          break;

        default:
          GOMA_EH(GOMA_ERROR, "Not valid porous shell index");
          break;
        }

      } else if (draining_curve == 0) /* Stay on main wetting curve */
      {
        sat_min = sat_min_wet;
        sat_max = mp->u_PorousShellCapPres[ipore][1];
        con_c = mp->u_PorousShellCapPres[ipore][2];
        nexp = mp->u_PorousShellCapPres[ipore][3];

        n_inv = 1.0 / nexp;
        mexp = 1.0 - n_inv;
        m_inv = -1.0 / mexp;

        /* Normalized saturation - should be between 0 to 1 */
        sat_norm = (saturation - sat_min) / (sat_max - sat_min);
        d_sat_norm_d_saturation = 1.0 / (sat_max - sat_min);

        brack_in = pow(sat_norm, m_inv) - 1.0;
        d_brack_in_d_sat_norm = m_inv * pow(sat_norm, (m_inv - 1.0));

        /* Evaluate capillary pressure and its sensitivities */
        cap_pres = con_c * pow(brack_in, n_inv);
        switch (ipore) {
        case 0:
          mp->d_cap_pres[SHELL_SAT_1] = n_inv * con_c * pow(brack_in, (n_inv - 1.0)) *
                                        d_brack_in_d_sat_norm * d_sat_norm_d_saturation;
          break;

        case 1:
          mp->d_cap_pres[SHELL_SAT_2] = n_inv * con_c * pow(brack_in, (n_inv - 1.0)) *
                                        d_brack_in_d_sat_norm * d_sat_norm_d_saturation;
          break;

        case 2:
          mp->d_cap_pres[SHELL_SAT_3] = n_inv * con_c * pow(brack_in, (n_inv - 1.0)) *
                                        d_brack_in_d_sat_norm * d_sat_norm_d_saturation;
          break;

        default:
          GOMA_EH(GOMA_ERROR, "Not valid porous shell index");
          break;
        }
      } else {
        GOMA_EH(GOMA_ERROR, "Either we are on draining curve or wetting curve!");
      }
    } else if (switch_now == 1) /* We are switching */
    {
      if ((draining_curve == 1) && (draining_curve_old == 0)) /* Wetting going to draining */
      {
        /* Get all parameters from draining curve */
        sat_min = mp->u_PorousShellCapPres[ipore][4];
        sat_max = mp->u_PorousShellCapPres[ipore][5];
        con_c = mp->u_PorousShellCapPres[ipore][6];
        nexp = mp->u_PorousShellCapPres[ipore][7];

        n_inv = 1.0 / nexp;
        mexp = 1.0 - n_inv;
        m_inv = -1.0 / mexp;

        /* Get the saturation and capillary pressure at switching or reversal point */
        sat_switch = pmv_hyst->sat_switch[ipore][ignode];
        cap_pres_switch = pmv_hyst->cap_pres_switch[ipore][ignode];

        /* Calculate normalized saturation at reversal point */
        brack_in_inv = pow((cap_pres_switch / con_c), nexp) + 1.0;
        sat_norm_switch = pow(brack_in_inv, (-mexp));

        /* Quality check */
        if ((sat_norm_switch > 1.0) || (sat_norm_switch < 0.0)) {
          GOMA_EH(GOMA_ERROR, "Invalid value of normalized saturation at switching point");
        }

        /* Calculate sat_max for the scanning draining curve */
        sat_max_dry = (1.0 / sat_norm_switch) * (sat_switch - sat_min * (1.0 - sat_norm_switch));

        /* More quality check */
        if (sat_max_dry > sat_max) {
          GOMA_EH(GOMA_ERROR, "Invalid value of maximum saturation at the scanning drying curve");
        }
        /* Update sat_max_dry in pmv_hyst structure*/
        pmv_hyst->sat_max_drain[ipore][ignode] = sat_max_dry;

        /* Calculate normalized saturation based on newly calculated sat_max_dry */
        sat_norm = (saturation - sat_min) / (sat_max_dry - sat_min);
        d_sat_norm_d_saturation = 1.0 / (sat_max_dry - sat_min);

        brack_in = pow(sat_norm, m_inv) - 1.0;
        d_brack_in_d_sat_norm = m_inv * pow(sat_norm, (m_inv - 1.0));

        /* Evaluate capillary pressure and its sensitivities */
        cap_pres = con_c * pow(brack_in, n_inv);

        switch (ipore) {
        case 0:
          mp->d_cap_pres[SHELL_SAT_1] = n_inv * con_c * pow(brack_in, (n_inv - 1.0)) *
                                        d_brack_in_d_sat_norm * d_sat_norm_d_saturation;
          break;

        case 1:
          mp->d_cap_pres[SHELL_SAT_2] = n_inv * con_c * pow(brack_in, (n_inv - 1.0)) *
                                        d_brack_in_d_sat_norm * d_sat_norm_d_saturation;
          break;

        case 2:
          mp->d_cap_pres[SHELL_SAT_3] = n_inv * con_c * pow(brack_in, (n_inv - 1.0)) *
                                        d_brack_in_d_sat_norm * d_sat_norm_d_saturation;
          break;

        default:
          GOMA_EH(GOMA_ERROR, "Not valid porous shell index");
          break;
        }
      } else if ((draining_curve == 0) && (draining_curve_old == 1)) /* draining going to wetting */
      {
        /* Get all parameters from wetting curve */
        sat_min = mp->u_PorousShellCapPres[ipore][0];
        sat_max = mp->u_PorousShellCapPres[ipore][1];
        con_c = mp->u_PorousShellCapPres[ipore][2];
        nexp = mp->u_PorousShellCapPres[ipore][3];

        n_inv = 1.0 / nexp;
        mexp = 1.0 - n_inv;
        m_inv = -1.0 / mexp;

        /* Get the saturation and capillary pressure at switching or reversal point */
        sat_switch = pmv_hyst->sat_switch[ipore][ignode];
        cap_pres_switch = pmv_hyst->cap_pres_switch[ipore][ignode];

        /* Calculate normalized saturation at reversal point */
        brack_in_inv = pow((cap_pres_switch / con_c), nexp) + 1.0;
        sat_norm_switch = pow(brack_in_inv, (-mexp));

        /* Quality check */
        if ((sat_norm_switch > 1.0) || (sat_norm_switch < 0.0)) {
          GOMA_EH(GOMA_ERROR, "Invalid value of normalized saturation at switching point");
        }

        /* Calculate sat_min for the scanning wetting curve */
        sat_min_wet = (sat_switch - sat_max * sat_norm_switch) / (1.0 - sat_norm_switch);

        /* More quality check */
        if (sat_min_wet < sat_min) {
          GOMA_EH(GOMA_ERROR, "Invalid value of minimum saturation at the scanning wetting curve");
        }
        /* Update sat_min_wet in pmv_hyst structure*/
        pmv_hyst->sat_min_imbibe[ipore][ignode] = sat_min_wet;

        /* Calculate normalized saturation based on newly calculated sat_min_wet */
        sat_norm = (saturation - sat_min_wet) / (sat_max - sat_min_wet);
        d_sat_norm_d_saturation = 1.0 / (sat_max - sat_min_wet);

        brack_in = pow(sat_norm, m_inv) - 1.0;
        d_brack_in_d_sat_norm = m_inv * pow(sat_norm, (m_inv - 1.0));

        /* Evaluate capillary pressure and its sensitivities */
        cap_pres = con_c * pow(brack_in, n_inv);
        switch (ipore) {
        case 0:
          mp->d_cap_pres[SHELL_SAT_1] = n_inv * con_c * pow(brack_in, (n_inv - 1.0)) *
                                        d_brack_in_d_sat_norm * d_sat_norm_d_saturation;
          break;

        case 1:
          mp->d_cap_pres[SHELL_SAT_2] = n_inv * con_c * pow(brack_in, (n_inv - 1.0)) *
                                        d_brack_in_d_sat_norm * d_sat_norm_d_saturation;
          break;

        case 2:
          mp->d_cap_pres[SHELL_SAT_3] = n_inv * con_c * pow(brack_in, (n_inv - 1.0)) *
                                        d_brack_in_d_sat_norm * d_sat_norm_d_saturation;
          break;

        default:
          GOMA_EH(GOMA_ERROR, "Not valid porous shell index");
          break;
        }
      }
    }
  } else if (mp->PorousShellCapPresModel[ipore] == VAN_GENUCHTEN_HYST_EXT) {
    /*
     * FOR VAN_GENUCHTEN_HYST_EXT EQUATION
     *  mp->u_saturation[0] is the irreducable water saturation for wetting curve 1
     *  mp->u_saturation[1] is the irreduceable air saturation for wetting curve 1
     *  mp->u_saturation[2] is entry pressure for wetting curve 1
     *  mp->u_saturation[3] is exponent n for wetting curve 1
     *
     *  mp->u_saturation[4] is the irreducable water saturation for draining curve 1
     *  mp->u_saturation[5] is the irreduceable air saturation for draining curve 1
     *  mp->u_saturation[6] is entry pressure for draining curve 1
     *  mp->u_saturation[7] is exponent n for draining curve 1
     *
     *  mp->u_saturation[8] is the irreducable water saturation for wetting curve 2
     *  mp->u_saturation[9] is the irreduceable air saturation for wetting curve 2
     *  mp->u_saturation[10] is entry pressure for wetting curve 2
     *  mp->u_saturation[11] is exponent n for wetting curve 2
     *
     *  mp->u_saturation[12] is the irreducable water saturation for draining curve 2
     *  mp->u_saturation[13] is the irreduceable air saturation for draining curve 2
     *  mp->u_saturation[14] is entry pressure for draining curve 2
     *  mp->u_saturation[15] is exponent n for draining curve 2

     *  Other input
     *  mp->u_saturation[16] Initial saturation curve for the material, viz.
     *                       if 1.0 then on the draining curve, and 0.0 then
     *                       on the wetting curve.
     *  mp->u_saturation[17] is Liquid_inventory_rate threshold for curve switching
     */

    /* Quality check of global node number ignode*/
    if (ignode < 0) {
      if (ipore == 0) {
        ignode = ei[pg->imtrx]->gnn_list[SHELL_SAT_1][ilnode];
      } else if (ipore == 1) {
        ignode = ei[pg->imtrx]->gnn_list[SHELL_SAT_2][ilnode];
      } else if (ipore == 2) {
        ignode = ei[pg->imtrx]->gnn_list[SHELL_SAT_3][ilnode];
      } else {
        GOMA_EH(GOMA_ERROR, "Something's wrong here!");
      }
    }

    /* First, see if we are switching */
    switch_now = pmv_hyst->curve_switch[ipore][ignode];

    i_ext_field = mp->por_shell_cap_pres_ext_field_index[ipore];
    /* Here I assume the external field value ranges from 0 to 1 */
    val_ext_field = *evp->external_field[i_ext_field][ilnode];

    /* Get the minimum and maximum saturation for the scanning curves */
    sat_min_wet = pmv_hyst->sat_min_imbibe[ipore][ignode];
    sat_max_dry = pmv_hyst->sat_max_drain[ipore][ignode];

    /* Find out whether on draining or wetting curve */
    draining_curve = pmv_hyst->curve_type[ipore][ignode];
    draining_curve_old = pmv_hyst->curve_type_old[ipore][ignode];

    /* If we are not switching */
    if (switch_now == 0) {
      if (draining_curve == 1) /* Stay on main draining curve */
      {

        /* Get parameters from draining curve 1*/
        sat_min_1 = mp->u_PorousShellCapPres[ipore][4];
        con_c_1 = mp->u_PorousShellCapPres[ipore][6];
        nexp_1 = mp->u_PorousShellCapPres[ipore][7];

        /* Get parameters from draining curve 2*/
        sat_min_2 = mp->u_PorousShellCapPres[ipore][12];
        con_c_2 = mp->u_PorousShellCapPres[ipore][14];
        nexp_2 = mp->u_PorousShellCapPres[ipore][15];

        /* Interpolate the fitting parameters with external field value */
        sat_min = sat_min_1 + val_ext_field * (sat_min_2 - sat_min_1);
        con_c = con_c_1 + val_ext_field * (con_c_2 - con_c_1);
        nexp = nexp_1 + val_ext_field * (nexp_2 - nexp_1);

        sat_max = sat_max_dry;

        n_inv = 1.0 / nexp;
        mexp = 1.0 - n_inv;
        m_inv = -1.0 / mexp;

        /* Normalized saturation - should be between 0 to 1 */
        sat_norm = (saturation - sat_min) / (sat_max - sat_min);
        d_sat_norm_d_saturation = 1.0 / (sat_max - sat_min);

        brack_in = pow(sat_norm, m_inv) - 1.0;
        d_brack_in_d_sat_norm = m_inv * pow(sat_norm, (m_inv - 1.0));

        /* Evaluate capillary pressure and its sensitivities */
        cap_pres = con_c * pow(brack_in, n_inv);
        mp->d_cap_pres[SHELL_SAT_1] = n_inv * con_c * pow(brack_in, (n_inv - 1.0)) *
                                      d_brack_in_d_sat_norm * d_sat_norm_d_saturation;
      } else if (draining_curve == 0) /* Stay on main wetting curve */
      {
        sat_min = sat_min_wet;

        /* Get parameters from wetting curve 1*/
        sat_max_1 = mp->u_PorousShellCapPres[ipore][1];
        con_c_1 = mp->u_PorousShellCapPres[ipore][2];
        nexp_1 = mp->u_PorousShellCapPres[ipore][3];

        /* Get parameters from wetting curve 2*/
        sat_max_2 = mp->u_PorousShellCapPres[ipore][9];
        con_c_2 = mp->u_PorousShellCapPres[ipore][10];
        nexp_2 = mp->u_PorousShellCapPres[ipore][11];

        /* Interpolate the fitting parameters with external field value */
        sat_max = sat_max_1 + val_ext_field * (sat_max_2 - sat_max_1);
        con_c = con_c_1 + val_ext_field * (con_c_2 - con_c_1);
        nexp = nexp_1 + val_ext_field * (nexp_2 - nexp_1);

        n_inv = 1.0 / nexp;
        mexp = 1.0 - n_inv;
        m_inv = -1.0 / mexp;

        /* Normalized saturation - should be between 0 to 1 */
        sat_norm = (saturation - sat_min) / (sat_max - sat_min);
        d_sat_norm_d_saturation = 1.0 / (sat_max - sat_min);

        brack_in = pow(sat_norm, m_inv) - 1.0;
        d_brack_in_d_sat_norm = m_inv * pow(sat_norm, (m_inv - 1.0));

        /* Evaluate capillary pressure and its sensitivities */
        cap_pres = con_c * pow(brack_in, n_inv);
        mp->d_cap_pres[SHELL_SAT_1] = n_inv * con_c * pow(brack_in, (n_inv - 1.0)) *
                                      d_brack_in_d_sat_norm * d_sat_norm_d_saturation;
      } else {
        GOMA_EH(GOMA_ERROR, "Either we are on draining curve or wetting curve!");
      }
    } else if (switch_now == 1) /* We are switching */
    {
      if ((draining_curve == 1) && (draining_curve_old == 0)) /* Wetting going to draining */
      {
        /* Get all parameters from draining curve 1 */
        sat_min_1 = mp->u_PorousShellCapPres[ipore][4];
        sat_max_1 = mp->u_PorousShellCapPres[ipore][5];
        con_c_1 = mp->u_PorousShellCapPres[ipore][6];
        nexp_1 = mp->u_PorousShellCapPres[ipore][7];

        /* Get all parameters from draining curve 2 */
        sat_min_2 = mp->u_PorousShellCapPres[ipore][12];
        sat_max_2 = mp->u_PorousShellCapPres[ipore][13];
        con_c_2 = mp->u_PorousShellCapPres[ipore][14];
        nexp_2 = mp->u_PorousShellCapPres[ipore][15];

        /* Interpolate the fitting parameters with external field value */
        sat_min = sat_min_1 + val_ext_field * (sat_min_2 - sat_min_1);
        sat_max = sat_max_1 + val_ext_field * (sat_max_2 - sat_max_1);
        con_c = con_c_1 + val_ext_field * (con_c_2 - con_c_1);
        nexp = nexp_1 + val_ext_field * (nexp_2 - nexp_1);

        n_inv = 1.0 / nexp;
        mexp = 1.0 - n_inv;
        m_inv = -1.0 / mexp;

        /* Get the saturation and capillary pressure at switching or reversal point */
        sat_switch = pmv_hyst->sat_switch[ipore][ignode];
        cap_pres_switch = pmv_hyst->cap_pres_switch[ipore][ignode];

        /* Calculate normalized saturation at reversal point */
        brack_in_inv = pow((cap_pres_switch / con_c), nexp) + 1.0;
        sat_norm_switch = pow(brack_in_inv, (-mexp));

        /* Quality check */
        if ((sat_norm_switch > 1.0) || (sat_norm_switch < 0.0)) {
          GOMA_EH(GOMA_ERROR, "Invalid value of normalized saturation at switching point");
        }

        /* Calculate sat_max for the scanning draining curve */
        sat_max_dry = (1.0 / sat_norm_switch) * (sat_switch - sat_min * (1.0 - sat_norm_switch));

        /* More quality check */
        if (sat_max_dry > sat_max) {
          GOMA_EH(GOMA_ERROR, "Invalid value of maximum saturation at the scanning drying curve");
        }
        /* Update sat_max_dry in element storage structure*/
        pmv_hyst->sat_max_drain[ipore][ignode] = sat_max_dry;

        /* Calculate normalized saturation based on newly calculated sat_max_dry */
        sat_norm = (saturation - sat_min) / (sat_max_dry - sat_min);
        d_sat_norm_d_saturation = 1.0 / (sat_max_dry - sat_min);

        brack_in = pow(sat_norm, m_inv) - 1.0;
        d_brack_in_d_sat_norm = m_inv * pow(sat_norm, (m_inv - 1.0));

        /* Evaluate capillary pressure and its sensitivities */
        cap_pres = con_c * pow(brack_in, n_inv);
        mp->d_cap_pres[SHELL_SAT_1] = n_inv * con_c * pow(brack_in, (n_inv - 1.0)) *
                                      d_brack_in_d_sat_norm * d_sat_norm_d_saturation;
      } else if ((draining_curve == 0) && (draining_curve_old == 1)) /* draining going to wetting */
      {
        /* Get all parameters from wetting curve 1 */
        sat_min_1 = mp->u_PorousShellCapPres[ipore][0];
        sat_max_1 = mp->u_PorousShellCapPres[ipore][1];
        con_c_1 = mp->u_PorousShellCapPres[ipore][2];
        nexp_1 = mp->u_PorousShellCapPres[ipore][3];

        /* Get all parameters from wetting curve 2 */
        sat_min_2 = mp->u_PorousShellCapPres[ipore][8];
        sat_max_2 = mp->u_PorousShellCapPres[ipore][9];
        con_c_2 = mp->u_PorousShellCapPres[ipore][10];
        nexp_2 = mp->u_PorousShellCapPres[ipore][11];

        /* Interpolate the fitting parameters with external field value */
        sat_min = sat_min_1 + val_ext_field * (sat_min_2 - sat_min_1);
        sat_max = sat_max_1 + val_ext_field * (sat_max_2 - sat_max_1);
        con_c = con_c_1 + val_ext_field * (con_c_2 - con_c_1);
        nexp = nexp_1 + val_ext_field * (nexp_2 - nexp_1);

        n_inv = 1.0 / nexp;
        mexp = 1.0 - n_inv;
        m_inv = -1.0 / mexp;

        /* Get the saturation and capillary pressure at switching or reversal point */
        sat_switch = pmv_hyst->sat_switch[ipore][ignode];
        cap_pres_switch = pmv_hyst->cap_pres_switch[ipore][ignode];

        /* Calculate normalized saturation at reversal point */
        brack_in_inv = pow((cap_pres_switch / con_c), nexp) + 1.0;
        sat_norm_switch = pow(brack_in_inv, (-mexp));

        /* Quality check */
        if ((sat_norm_switch > 1.0) || (sat_norm_switch < 0.0)) {
          GOMA_EH(GOMA_ERROR, "Invalid value of normalized saturation at switching point");
        }

        /* Calculate sat_min for the scanning wetting curve */
        sat_min_wet = (sat_switch - sat_max * sat_norm_switch) / (1.0 - sat_norm_switch);

        /* More quality check */
        if (sat_min_wet < sat_min) {
          GOMA_EH(GOMA_ERROR, "Invalid value of minimum saturation at the scanning wetting curve");
        }
        /* Update sat_min_wet in element storage structure*/
        pmv_hyst->sat_min_imbibe[ipore][ignode] = sat_min_wet;

        /* Calculate normalized saturation based on newly calculated sat_min_wet */
        sat_norm = (saturation - sat_min_wet) / (sat_max - sat_min_wet);
        d_sat_norm_d_saturation = 1.0 / (sat_max - sat_min_wet);

        brack_in = pow(sat_norm, m_inv) - 1.0;
        d_brack_in_d_sat_norm = m_inv * pow(sat_norm, (m_inv - 1.0));

        /* Evaluate capillary pressure and its sensitivities */
        cap_pres = con_c * pow(brack_in, n_inv);
        mp->d_cap_pres[SHELL_SAT_1] = n_inv * con_c * pow(brack_in, (n_inv - 1.0)) *
                                      d_brack_in_d_sat_norm * d_sat_norm_d_saturation;
      }
    } else {
      GOMA_EH(GOMA_ERROR, "To switch or not to switch; that is the question");
    }
  } else {
    GOMA_EH(GOMA_ERROR, "No models for capillary pressure");
  }
  return cap_pres;
} /* end   load_cap_pres  */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void load_gas_conc(double porosity, double cap_pres, double d_cap_pres[2])

/************************************************************************
 * load_gas_conc -- calculate gas phase concentrations in a porous
 *  media from the unknowns used in the problem.  USES KELVIN EQUATION.
 *  see load_gas_conc_flat for flat interface vapor pressure
 *
 *  input:
 *  -----
 *
 *
 *  output:
 *  -----
 *   Calculates the gas phase concentrations and its first
 *   and second derivatives with respect to all the problem unknowns
 *      pmv->gas_density_solvents[w]
 *      mp->porous_vapor_pressure[w]
 ************************************************************************/
{
  int w, w1, w2;
  double m_rta, m_rtw, ma_mw, pc_rhort = 0.0;
  double m_rta_old = 0.0, m_rtw_old = 0.0, pc_rhort_old;
  double rho_l, rho_sat, Pg, Psat, T;
  double drhosat_dT = 0.0, d_drhosat_dT = 0.0, dPsat_dT, d_dPsat_dT;
  double Pg_old = 0.0, rhosat_old = 0.0, Psat_old = 0.0, T_old = 0.0;
  double drhosat_dT_old = 0.0, dPsat_dT_old, d_dPsat_dT_old;
  double term_1 = 0.0, term_2 = 0.0, sum_1 = 0.0;
  double term_1o, term_2o;
  const int i_pl = 0, i_pg = 1, i_pe = 3;
  int first_porous_var = POR_LAST - MAX_POROUS_NUM + 1;

  for (w = 0; w < MAX_PMV; w++) {
    pmv->gas_density_solvents[w] = 0.;
    mp->porous_vapor_pressure[w] = 0.;
  }
  /*
   *  Find which model is used for gas conc. and calculate:
   *    1)  the saturation, pmv->gas_density_solvents
   *    2)  the sensitivity of gas conc. to all concentrations, porosity and temperature
   *          i.e. the first derivative of gas conc. w.r.t. each variable
   *           put this in pmv->d_gas_density_solvents[var]
   *    3)  the second derivatives of gas conc. w.r.t. each variable, including
   *          cross-terms, put this in pmv->d_d_gas_density_solvents[var][var]
   *          second derivative is needed for sensitivity of fluxes which depend
   *          on the first derivative of gas conc.
   */
  if (mp->PorousVaporPressureModel[i_pl] == KELVIN &&
      ((mp->PorousMediaType == POROUS_TWO_PHASE && mp->PorousGasConstantsModel) ||
       mp->PorousMediaType == POROUS_UNSATURATED ||
       mp->PorousMediaType == POROUS_SHELL_UNSATURATED)) {
    /* currently assume isothermal unless the energy equation is being solved.
     * if it is being solved, replace the temperatures below with the
     * temperature field variable, and compute sensitivities.
     *
     * FOR KELVIN EQUATION
     *  mp->u_porous_vapor_pressure[i_pl][0] is the flat interface vapor pressure
     *  mp->u_porous_vapor_pressure[i_pl][1] is the liquid density
     *  mp->u_porous_vapor_pressure[i_pl][2] is the water molecular weight
     *  mp->u_porous_vapor_pressure[i_pl][3] is the gas law constant
     *  mp->u_porous_vapor_pressure[i_pl][4] is the temperature
     *
     * FOR IDEAL_GAS EQUATION --- for TWO_PHASE FLOW and UNSATURATED Flow
     *  mp->u_porous_gas_constants[0] is the air molecular weight
     *  mp->u_porous_gas_constants[1] is the gas law constant
     *  mp->u_porous_gas_constants[2] is the temperature
     *  mp->u_porous_gas_constants[3] is the ambient base pressure of gas phase
     */

    /* Note - liquid density is a constant taken from the input file
     *        and is ASSUMED CONSTANT
     */
    rho_l = mp->u_porous_vapor_pressure[i_pl][1];

    /* If solving the gas/energy equation, replace constants in the above model
     * with values consistent with current solution variables
     */
    if (pd->e[pg->imtrx][R_POR_GAS_PRES])
      Pg = fv->p_gas;
    else
      Pg = mp->u_porous_gas_constants[3];
    if (pd->e[pg->imtrx][R_POR_ENERGY]) {
      T = fv->T;
      Psat = P_water_sat_EOS(T, &dPsat_dT, &d_dPsat_dT);
      T = T + 273.15;
    } else {
      T = mp->u_porous_vapor_pressure[i_pl][4];
      Psat = mp->u_porous_vapor_pressure[i_pl][0];
    }

    if (pd->TimeIntegration == TRANSIENT) {
      if (pd->e[pg->imtrx][R_POR_GAS_PRES])
        Pg_old = fv_old->p_gas;
      else
        Pg_old = Pg;
      if (pd->e[pg->imtrx][R_POR_ENERGY]) {
        T_old = fv_old->T;
        Psat_old = P_water_sat_EOS(T_old, &dPsat_dT_old, &d_dPsat_dT_old);
        T_old = T_old + 273.15;
      } else {
        T_old = T;
        Psat_old = Psat;
        dPsat_dT_old = 0.0;
        d_dPsat_dT_old = 0.0;
      }
    }

    /* precalculate some parameters */

    m_rtw = mp->u_porous_vapor_pressure[i_pl][2] / mp->u_porous_vapor_pressure[i_pl][3] / T;

    if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
      m_rta = mp->u_porous_gas_constants[0] / mp->u_porous_gas_constants[1] / T;
      ma_mw = mp->u_porous_gas_constants[0] / mp->u_porous_vapor_pressure[i_pl][2];
    } else {
      m_rta = mp->u_porous_gas_constants[0] / mp->u_porous_vapor_pressure[i_pl][3] / T;
      ma_mw = mp->u_porous_gas_constants[0] / mp->u_porous_vapor_pressure[i_pl][2];
    }

    rho_sat = m_rtw * Psat;

    if (pd->e[pg->imtrx][R_POR_ENERGY]) {
      drhosat_dT = m_rtw * (dPsat_dT - Psat / T);
      d_drhosat_dT = m_rtw * (d_dPsat_dT - 2.0 / T * dPsat_dT + 2.0 / T / T * Psat);
    }

    /* KELVIN for water conc in gas */

    pmv->gas_density_solvents[i_pl] = rho_sat * exp(-cap_pres * m_rtw / rho_l);
    mp->porous_vapor_pressure[i_pl] = pmv->gas_density_solvents[i_pl] / m_rtw;

    /* IDEAL for air conc in gas ->
     *      air pressure = total pressure - water pressure
     */
    if (pd->e[pg->imtrx][R_POR_GAS_PRES])
      pmv->gas_density_solvents[i_pg] = Pg * m_rta - pmv->gas_density_solvents[i_pl] * ma_mw;

    if (pd->e[pg->imtrx][R_POR_ENERGY])
      pmv->gas_density_solvents[i_pe] = pmv->gas_density_solvents[i_pl] * pmv->enthalpy[1] +
                                        pmv->gas_density_solvents[i_pg] * pmv->enthalpy[2] - Pg;

    if (pd->TimeIntegration == TRANSIENT) {
      m_rtw_old =
          mp->u_porous_vapor_pressure[i_pl][2] / mp->u_porous_vapor_pressure[i_pl][3] / T_old;
      m_rta_old = mp->u_porous_gas_constants[0] / mp->u_porous_gas_constants[1] / T_old;
      rhosat_old = m_rtw_old * Psat_old;
      drhosat_dT_old = m_rtw_old * (dPsat_dT_old - Psat_old / T_old);

      pmv_old->gas_density_solvents[i_pl] =
          rhosat_old * exp(-pmv_old->cap_pres * m_rtw_old / rho_l);
      mp_old->porous_vapor_pressure[i_pl] = pmv_old->gas_density_solvents[i_pl] / m_rtw_old;

      if (pd->e[pg->imtrx][R_POR_GAS_PRES])
        pmv_old->gas_density_solvents[i_pg] =
            Pg_old * m_rta_old - pmv_old->gas_density_solvents[i_pl] * ma_mw;

      if (pd->e[pg->imtrx][R_POR_ENERGY])
        pmv_old->gas_density_solvents[i_pe] =
            pmv_old->gas_density_solvents[i_pl] * pmv_old->enthalpy[1] +
            pmv_old->gas_density_solvents[i_pg] * pmv_old->enthalpy[2] - Pg_old;
    }

    /* *** NOTE --- this model of vapor-liquid equilibrium has no dependence
     *              upon porosity, saturation, or pore-size
     */

    /* sensitivities */
    pmv->d_gas_density_solvents[i_pl][POR_LIQ_PRES] =
        -pmv->gas_density_solvents[i_pl] * m_rtw / rho_l * d_cap_pres[i_pl];
    mp->d_porous_vapor_pressure[i_pl][POR_LIQ_PRES] =
        pmv->d_gas_density_solvents[i_pl][POR_LIQ_PRES] / m_rtw;

    if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
      pmv->d_gas_density_solvents[i_pl][POR_GAS_PRES] =
          -pmv->d_gas_density_solvents[i_pl][POR_LIQ_PRES];
      mp->d_porous_vapor_pressure[i_pl][POR_GAS_PRES] =
          pmv->d_gas_density_solvents[i_pl][POR_GAS_PRES] / m_rtw;

      pmv->d_gas_density_solvents[i_pg][POR_LIQ_PRES] =
          -pmv->d_gas_density_solvents[i_pl][POR_LIQ_PRES] * ma_mw;
      pmv->d_gas_density_solvents[i_pg][POR_GAS_PRES] =
          m_rta - pmv->d_gas_density_solvents[i_pl][POR_GAS_PRES] * ma_mw;
    }

    if (pd->e[pg->imtrx][R_POR_ENERGY]) {

      /*  Pc / (rho_l * R * T) */
      pc_rhort = pmv->cap_pres / m_rtw / rho_l;

      /* Add energy terms to water eqn senstivities */
      term_1 = pmv->gas_density_solvents[i_pl] * pc_rhort / T;
      term_2 = drhosat_dT * pmv->gas_density_solvents[i_pl] / rho_sat;
      sum_1 = -1.0 / T + (term_1 + term_2) / pmv->gas_density_solvents[i_pl];
      pmv->d_gas_density_solvents[i_pl][POR_TEMP] = term_1 + term_2;
      pmv->d_gas_density_solvents[i_pe][POR_LIQ_PRES] =
          pmv->enthalpy[1] * pmv->d_gas_density_solvents[i_pl][POR_LIQ_PRES] +
          pmv->enthalpy[2] * pmv->d_gas_density_solvents[i_pg][POR_LIQ_PRES];

      /* Add energy terms to air eqn sensitivies */
      if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
        pmv->d_gas_density_solvents[i_pg][POR_TEMP] =
            -Pg * m_rta / T - pmv->d_gas_density_solvents[i_pl][POR_TEMP] * ma_mw;
        pmv->d_gas_density_solvents[i_pe][POR_GAS_PRES] =
            pmv->gas_density_solvents[i_pl] * pmv->d_enthalpy[1][POR_GAS_PRES] +
            pmv->enthalpy[1] * pmv->d_gas_density_solvents[i_pl][POR_GAS_PRES] +
            pmv->gas_density_solvents[i_pg] * pmv->d_enthalpy[2][POR_GAS_PRES] +
            pmv->enthalpy[2] * pmv->d_gas_density_solvents[i_pg][POR_GAS_PRES] - 1.0;
      }

      pmv->d_gas_density_solvents[i_pe][POR_TEMP] =
          pmv->gas_density_solvents[i_pl] * pmv->d_enthalpy[1][POR_TEMP] +
          pmv->enthalpy[1] * pmv->d_gas_density_solvents[i_pl][POR_TEMP] +
          pmv->gas_density_solvents[i_pg] * pmv->d_enthalpy[2][POR_TEMP] +
          pmv->enthalpy[2] * pmv->d_gas_density_solvents[i_pg][POR_TEMP];
    }

    if (pd->TimeIntegration == TRANSIENT) {
      pmv_old->d_gas_density_solvents[i_pl][POR_LIQ_PRES] =
          -pmv_old->gas_density_solvents[i_pl] * m_rtw_old / rho_l * d_cap_pres[i_pl];
      mp_old->d_porous_vapor_pressure[i_pl][POR_LIQ_PRES] =
          pmv_old->d_gas_density_solvents[i_pl][POR_LIQ_PRES] / m_rtw_old;

      if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
        pmv_old->d_gas_density_solvents[i_pl][POR_GAS_PRES] =
            -pmv_old->d_gas_density_solvents[i_pl][POR_LIQ_PRES];

        mp_old->d_porous_vapor_pressure[i_pl][POR_GAS_PRES] =
            pmv_old->d_gas_density_solvents[i_pl][POR_GAS_PRES] / m_rtw_old;

        pmv_old->d_gas_density_solvents[i_pg][POR_GAS_PRES] =
            m_rta_old - pmv_old->d_gas_density_solvents[i_pl][POR_GAS_PRES] * ma_mw;

        pmv_old->d_gas_density_solvents[i_pg][POR_LIQ_PRES] =
            -pmv_old->d_gas_density_solvents[i_pl][POR_LIQ_PRES] * ma_mw;
      }

      if (pd->e[pg->imtrx][R_POR_ENERGY]) {

        /*  Pc / (rho_l * R * T) */
        pc_rhort_old = pmv_old->cap_pres / m_rtw_old / rho_l;

        /* Add energy terms to water eqn senstivities */
        term_1o = pmv_old->gas_density_solvents[i_pl] * pc_rhort_old / T_old;
        term_2o = drhosat_dT_old * pmv_old->gas_density_solvents[i_pl] / rhosat_old;
        pmv_old->d_gas_density_solvents[i_pl][POR_TEMP] = term_1o + term_2o;
        pmv_old->d_gas_density_solvents[i_pe][POR_LIQ_PRES] =
            pmv_old->enthalpy[1] * pmv_old->d_gas_density_solvents[i_pl][POR_LIQ_PRES] +
            pmv_old->enthalpy[2] * pmv_old->d_gas_density_solvents[i_pg][POR_LIQ_PRES];

        /* Add energy terms to air eqn sensitivies */
        if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
          pmv_old->d_gas_density_solvents[i_pg][POR_TEMP] =
              -Pg_old * m_rta_old / T_old - pmv_old->d_gas_density_solvents[i_pl][POR_TEMP] * ma_mw;
          pmv_old->d_gas_density_solvents[i_pe][POR_GAS_PRES] =
              pmv_old->gas_density_solvents[i_pl] * pmv_old->d_enthalpy[1][POR_GAS_PRES] +
              pmv_old->enthalpy[1] * pmv_old->d_gas_density_solvents[i_pl][POR_GAS_PRES] +
              pmv_old->gas_density_solvents[i_pg] * pmv_old->d_enthalpy[2][POR_GAS_PRES] +
              pmv_old->enthalpy[2] * pmv_old->d_gas_density_solvents[i_pg][POR_GAS_PRES] - 1.0;
        }

        pmv_old->d_gas_density_solvents[i_pe][POR_TEMP] =
            pmv_old->gas_density_solvents[i_pl] * pmv_old->d_enthalpy[1][POR_TEMP] +
            pmv_old->enthalpy[1] * pmv_old->d_gas_density_solvents[i_pl][POR_TEMP] +
            pmv_old->gas_density_solvents[i_pg] * pmv_old->d_enthalpy[2][POR_TEMP] +
            pmv_old->enthalpy[2] * pmv_old->d_gas_density_solvents[i_pg][POR_TEMP];
      }
    }

    /* second derivatives for full Newton of gradients */

    if (af->Assemble_Jacobian) {
      pmv->d_d_gas_density_solvents[i_pl][POR_LIQ_PRES][POR_LIQ_PRES] =
          pmv->gas_density_solvents[i_pl] * m_rtw * m_rtw / rho_l / rho_l * d_cap_pres[i_pl] *
          d_cap_pres[i_pl];
      /* replace above with the next statement */
      pmv->d_d_gas_density_solvents[i_pl][POR_LIQ_PRES][POR_LIQ_PRES] =
          pmv->d_gas_density_solvents[i_pl][POR_LIQ_PRES] *
          pmv->d_gas_density_solvents[i_pl][POR_LIQ_PRES] / pmv->gas_density_solvents[i_pl];

      if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
        pmv->d_d_gas_density_solvents[i_pg][POR_LIQ_PRES][POR_LIQ_PRES] =
            -pmv->d_d_gas_density_solvents[i_pl][POR_LIQ_PRES][POR_LIQ_PRES] * ma_mw;

        pmv->d_d_gas_density_solvents[i_pl][POR_LIQ_PRES][POR_GAS_PRES] =
            pmv->gas_density_solvents[i_pl] * m_rtw * m_rtw / rho_l / rho_l * d_cap_pres[i_pl] *
            d_cap_pres[i_pg];
        /* replace above with the next statement */
        pmv->d_d_gas_density_solvents[i_pl][POR_LIQ_PRES][POR_GAS_PRES] =
            pmv->d_gas_density_solvents[i_pl][POR_GAS_PRES] *
            pmv->d_gas_density_solvents[i_pl][POR_LIQ_PRES] / pmv->gas_density_solvents[i_pl];

        pmv->d_d_gas_density_solvents[i_pg][POR_LIQ_PRES][POR_GAS_PRES] =
            -pmv->d_d_gas_density_solvents[i_pl][POR_LIQ_PRES][POR_GAS_PRES] * ma_mw;

        pmv->d_d_gas_density_solvents[i_pl][POR_GAS_PRES][POR_LIQ_PRES] =
            pmv->gas_density_solvents[i_pl] * m_rtw * m_rtw / rho_l / rho_l * d_cap_pres[i_pg] *
            d_cap_pres[i_pl];
        /* replace above with the next statement */
        pmv->d_d_gas_density_solvents[i_pl][POR_GAS_PRES][POR_LIQ_PRES] =
            pmv->d_gas_density_solvents[i_pl][POR_LIQ_PRES] *
            pmv->d_gas_density_solvents[i_pl][POR_GAS_PRES] / pmv->gas_density_solvents[i_pl];

        pmv->d_d_gas_density_solvents[i_pg][POR_GAS_PRES][POR_LIQ_PRES] =
            -pmv->d_d_gas_density_solvents[i_pl][POR_GAS_PRES][POR_LIQ_PRES] * ma_mw;

        pmv->d_d_gas_density_solvents[i_pl][POR_GAS_PRES][POR_GAS_PRES] =
            pmv->gas_density_solvents[i_pl] * m_rtw * m_rtw / rho_l / rho_l * d_cap_pres[i_pg] *
            d_cap_pres[i_pg];
        /* replace above with the next statement */
        pmv->d_d_gas_density_solvents[i_pl][POR_GAS_PRES][POR_GAS_PRES] =
            pmv->d_gas_density_solvents[i_pl][POR_GAS_PRES] *
            pmv->d_gas_density_solvents[i_pl][POR_GAS_PRES] / pmv->gas_density_solvents[i_pl];

        pmv->d_d_gas_density_solvents[i_pg][POR_GAS_PRES][POR_GAS_PRES] =
            -pmv->d_d_gas_density_solvents[i_pl][POR_GAS_PRES][POR_GAS_PRES] * ma_mw;
      }

      if (pd->e[pg->imtrx][R_POR_ENERGY]) {
        /* Add energy terms to water eqn senstivities */
        pmv->d_d_gas_density_solvents[i_pl][POR_LIQ_PRES][POR_TEMP] =
            pmv->d_gas_density_solvents[i_pl][POR_LIQ_PRES] *
            (pmv->d_gas_density_solvents[i_pl][POR_TEMP] / pmv->gas_density_solvents[i_pl] -
             1.0 / T);

        pmv->d_d_gas_density_solvents[i_pl][POR_GAS_PRES][POR_TEMP] =
            pmv->d_gas_density_solvents[i_pl][POR_GAS_PRES] *
            (pmv->d_gas_density_solvents[i_pl][POR_TEMP] / pmv->gas_density_solvents[i_pl] -
             1.0 / T);

        pmv->d_d_gas_density_solvents[i_pl][POR_TEMP][POR_LIQ_PRES] =
            pmv->d_gas_density_solvents[i_pl][POR_LIQ_PRES] * sum_1;

        pmv->d_d_gas_density_solvents[i_pl][POR_TEMP][POR_GAS_PRES] =
            pmv->d_gas_density_solvents[i_pl][POR_GAS_PRES] * sum_1;

        pmv->d_d_gas_density_solvents[i_pl][POR_TEMP][POR_TEMP] =
            term_1 * (-2.0 / T + 1.0 / pmv->gas_density_solvents[i_pl] *
                                     (pmv->d_gas_density_solvents[i_pl][POR_TEMP] + term_2)) +
            d_drhosat_dT * exp(-pc_rhort);

        /* Add energy terms to air eqn sensitivies */
        if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
          pmv->d_d_gas_density_solvents[i_pg][POR_LIQ_PRES][POR_TEMP] =
              -pmv->d_d_gas_density_solvents[i_pl][POR_LIQ_PRES][POR_TEMP] * ma_mw;
          pmv->d_d_gas_density_solvents[i_pg][POR_GAS_PRES][POR_TEMP] =
              -pmv->d_d_gas_density_solvents[i_pl][POR_GAS_PRES][POR_TEMP] * ma_mw -
              1.0 / m_rta / T;

          pmv->d_d_gas_density_solvents[i_pg][POR_TEMP][POR_LIQ_PRES] =
              -pmv->d_d_gas_density_solvents[i_pl][POR_TEMP][POR_LIQ_PRES] * ma_mw;
          pmv->d_d_gas_density_solvents[i_pg][POR_TEMP][POR_GAS_PRES] =
              -pmv->d_d_gas_density_solvents[i_pl][POR_TEMP][POR_GAS_PRES] * ma_mw -
              1.0 / m_rta / T;
          pmv->d_d_gas_density_solvents[i_pg][POR_TEMP][POR_TEMP] =
              -pmv->d_d_gas_density_solvents[i_pl][POR_TEMP][POR_TEMP] * ma_mw +
              2.0 * Pg / m_rta / T / T;
        }

        /* Generalize the energy equation sensitivities to allow looping; only the
         * appropriate contributing terms are nonzero.
         */

        for (w1 = first_porous_var; w1 < POR_LAST + 1; w1++) {
          for (w2 = first_porous_var; w2 < POR_LAST + 1; w2++)
            pmv->d_d_gas_density_solvents[i_pe][w1][w2] =
                pmv->d_d_gas_density_solvents[i_pl][w1][w2] * pmv->enthalpy[1] +
                pmv->d_d_gas_density_solvents[i_pg][w1][w2] * pmv->enthalpy[2] +
                pmv->d_gas_density_solvents[i_pl][w1] * pmv->d_enthalpy[1][w2] +
                pmv->d_gas_density_solvents[i_pg][w1] * pmv->d_enthalpy[2][w2] +
                pmv->d_gas_density_solvents[i_pl][w2] * pmv->d_enthalpy[1][w1] +
                pmv->d_gas_density_solvents[i_pg][w2] * pmv->d_enthalpy[2][w1] +
                pmv->gas_density_solvents[i_pl] * pmv->d_d_enthalpy[1][w1][w2] +
                pmv->gas_density_solvents[i_pg] * pmv->d_d_enthalpy[2][w1][w2];
        }
      }
    }
  } else {
    GOMA_EH(GOMA_ERROR, "No models for gas phase conc.");
  }
  return;
} /*  end of   load_gas_conc    */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*
 * load_gas_conc_flat -- calculate gas phase concentrations in a porous
 *  media from the unknowns used in the problem.  USES curvature INDEPENDENT model.
 *  see load_gas_conc for KELVIN equivalent
 *
 *  input:   assume that load_fv and load_fv_grads and load_fv_mesh_derivs have all
 *  -----    been called already, so unknowns and their sensitivies are known at this
 *           gauss point
 *
 *  output:  calculates the gas phase concentrations and its first and second derivatives with
 *  -----    respect to all the problem unknowns
 *
 */

void load_gas_conc_flat(double porosity, double cap_pres, double d_cap_pres[2])

/*************************************************************************/
/*
 * load_gas_conc_flat -- calculate gas phase concentrations in a porous
 *  media from the unknowns used in the problem.  USES curvature INDEPENDENT model.
 *  see load_gas_conc for KELVIN equivalent
 *
 *  input:   assume that load_fv and load_fv_grads and load_fv_mesh_derivs have all
 *  -----    been called already, so unknowns and their sensitivies are known at this
 *           gauss point
 *
 *  output:  calculates the gas phase concentrations and its first and second derivatives with
 *  -----    respect to all the problem unknowns
 */
{
  int w, w1, w2;
  double m_rta, m_rtw, ma_mw;
  double m_rta_old = 0.0, m_rtw_old = 0.0;
  double rho_sat, Pg, Psat, T;
  double TC;
  double drho_wgsat_dP, drho_wgsat_dT, d_drho_wgsat_dP[2], d_drho_wgsat_dT[2];
  double drhosat_dT = 0.0, d_drhosat_dT = 0.0, dPsat_dT = 0.0, d_dPsat_dT = 0.0;
  double Pg_old = 0.0, rhosat_old = 0.0, Psat_old = 0.0, T_old = 0.0;
  double drhosat_dT_old = 0.0, dPsat_dT_old, d_dPsat_dT_old;
  const int i_pl = 0, i_pg = 1, i_pe = 3;
  int first_porous_var = POR_LAST - MAX_POROUS_NUM + 1;

  for (w = 0; w < MAX_PMV; w++) {
    pmv->gas_density_solvents[w] = 0.;
    mp->porous_vapor_pressure[w] = 0.;
  }
  /*
   *  Find which model is used for gas conc. and calculate:
   *    1)  the saturation, pmv->gas_density_solvents
   *    2)  the sensitivity of gas conc. to all concentrations, porosity and temperature
   *          i.e. the first derivative of gas conc. w.r.t. each variable
   *           put this in pmv->d_gas_density_solvents[var]
   *    3)  the second derivatives of gas conc. w.r.t. each variable, including
   *          cross-terms, put this in pmv->d_d_gas_density_solvents[var][var]
   *          second derivative is needed for sensitivity of fluxes which depend
   *          on the first derivative of gas conc.
   */
  if (mp->PorousVaporPressureModel[i_pl] == FLAT &&
      ((mp->PorousMediaType == POROUS_TWO_PHASE && mp->PorousGasConstantsModel) ||
       mp->PorousMediaType == POROUS_UNSATURATED)) {
    /* currently assume isothermal
     * FOR FLAT EQUATION
     *  mp->u_porous_vapor_pressure[i_pl][0] is the flat interface vapor pressure
     *  mp->u_porous_vapor_pressure[i_pl][1] is the liquid density
     *  mp->u_porous_vapor_pressure[i_pl][2] is the water molecular weight
     *  mp->u_porous_vapor_pressure[i_pl][3] is the gas law constant
     *  mp->u_porous_vapor_pressure[i_pl][4] is the temperature
     *
     * FOR IDEAL_GAS EQUATION --- for TWO_PHASE FLOW or UNSATURATED FLOW
     *  mp->u_porous_gas_constants[0] is the air molecular weight
     *  mp->u_porous_gas_constants[1] is the gas law constant
     *  mp->u_porous_gas_constants[2] is the temperature
     *  mp->u_porous_gas_constants[3] is the ambient base pressure of gas phase
     */

    /* If solving the gas/energy equation, replace constants in the above model
     * with values consistent with current solution variables
     */
    if (pd->e[pg->imtrx][R_POR_GAS_PRES])
      Pg = fv->p_gas;
    else
      Pg = mp->u_porous_gas_constants[3];
    if (pd->e[pg->imtrx][R_POR_ENERGY]) {
      T = fv->T;
      Psat = P_water_sat_EOS(T, &dPsat_dT, &d_dPsat_dT);
      T = T + 273.15;
    } else {
      T = mp->u_porous_vapor_pressure[i_pl][4];
      Psat = mp->u_porous_vapor_pressure[i_pl][0];
    }

    if (pd->TimeIntegration == TRANSIENT) {
      if (pd->e[pg->imtrx][R_POR_GAS_PRES])
        Pg_old = fv_old->p_gas;
      else
        Pg_old = Pg;
      if (pd->e[pg->imtrx][R_POR_ENERGY]) {
        T_old = fv_old->T;
        Psat_old = P_water_sat_EOS(T_old, &dPsat_dT_old, &d_dPsat_dT_old);
        T_old = T_old + 273.15;
      } else {
        T_old = T;
        Psat_old = Psat;
        dPsat_dT_old = 0.0;
        d_dPsat_dT_old = 0.0;
      }
    }

    /* precalculate some parameters */

    m_rtw = mp->u_porous_vapor_pressure[i_pl][2] / mp->u_porous_vapor_pressure[i_pl][3] / T;

    if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
      m_rta = mp->u_porous_gas_constants[0] / mp->u_porous_gas_constants[1] / T;
      ma_mw = mp->u_porous_gas_constants[0] / mp->u_porous_vapor_pressure[i_pl][2];
    } else {
      m_rta = mp->u_porous_gas_constants[0] / mp->u_porous_vapor_pressure[i_pl][3] / T;
      ma_mw = mp->u_porous_gas_constants[0] / mp->u_porous_vapor_pressure[i_pl][2];
    }

    rho_sat = m_rtw * Psat;
    TC = T - 273.15;
    rho_sat_water_vap_EOS(Psat, TC, &drho_wgsat_dP, &drho_wgsat_dT, d_drho_wgsat_dP,
                          d_drho_wgsat_dT);

    if (pd->e[pg->imtrx][R_POR_ENERGY]) {
      drhosat_dT = m_rtw * (dPsat_dT - Psat / T);
      d_drhosat_dT = m_rtw * (d_dPsat_dT - 2.0 / T * dPsat_dT + 2.0 / T / T * Psat);
    }

    /* Curvature independent model for VP and water conc in gas
     * Note here how the mass concenration is weighted by the saturation level
     * I.e. them more menisci that abound, the lower saturation and thelower the conc
     */
    pmv->gas_density_solvents[i_pl] = rho_sat * mp->saturation;
    mp->porous_vapor_pressure[i_pl] = pmv->gas_density_solvents[i_pl] / m_rtw / mp->saturation;

    /* IDEAL for air conc in gas ->
     *      air pressure = total pressure - water pressure
     */

    if (pd->e[pg->imtrx][R_POR_GAS_PRES])
      pmv->gas_density_solvents[i_pg] = Pg * m_rta - pmv->gas_density_solvents[i_pl] * ma_mw;

    if (pd->e[pg->imtrx][R_POR_ENERGY])
      pmv->gas_density_solvents[i_pe] = pmv->gas_density_solvents[i_pl] * pmv->enthalpy[1] +
                                        pmv->gas_density_solvents[i_pg] * pmv->enthalpy[2] - Pg;

    if (pd->TimeIntegration == TRANSIENT) {
      m_rtw_old =
          mp->u_porous_vapor_pressure[i_pl][2] / mp->u_porous_vapor_pressure[i_pl][3] / T_old;
      m_rta_old = mp->u_porous_gas_constants[0] / mp->u_porous_gas_constants[1] / T_old;
      rhosat_old = m_rtw_old * Psat_old;
      drhosat_dT_old = m_rtw_old * (dPsat_dT_old - Psat_old / T_old);

      pmv_old->gas_density_solvents[i_pl] = rhosat_old * mp_old->saturation;
      mp_old->porous_vapor_pressure[i_pl] =
          pmv_old->gas_density_solvents[i_pl] / (m_rtw_old * mp_old->saturation);

      if (pd->e[pg->imtrx][R_POR_GAS_PRES])
        pmv_old->gas_density_solvents[i_pg] =
            Pg_old * m_rta_old - pmv_old->gas_density_solvents[i_pl] * ma_mw;

      if (pd->e[pg->imtrx][R_POR_ENERGY])
        pmv_old->gas_density_solvents[i_pe] =
            pmv_old->gas_density_solvents[i_pl] * pmv_old->enthalpy[1] +
            pmv_old->gas_density_solvents[i_pg] * pmv_old->enthalpy[2] - Pg_old;
    }

    /* *** NOTE --- this model of vapor-liquid equilibrium has no dependence
     *              upon porosity, saturation, or pore-size
     */

    /* sensitivities */

    pmv->d_gas_density_solvents[i_pl][POR_LIQ_PRES] = rho_sat * mp->d_saturation[POR_LIQ_PRES];
    mp->d_porous_vapor_pressure[i_pl][POR_LIQ_PRES] =
        pmv->d_gas_density_solvents[i_pl][POR_LIQ_PRES] / m_rtw;

    if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
      pmv->d_gas_density_solvents[i_pl][POR_GAS_PRES] =
          -pmv->d_gas_density_solvents[i_pl][POR_LIQ_PRES];
      mp->d_porous_vapor_pressure[i_pl][POR_GAS_PRES] =
          pmv->d_gas_density_solvents[i_pl][POR_GAS_PRES] / m_rtw;

      pmv->d_gas_density_solvents[i_pg][POR_LIQ_PRES] =
          -pmv->d_gas_density_solvents[i_pl][POR_LIQ_PRES] * ma_mw;
      pmv->d_gas_density_solvents[i_pg][POR_GAS_PRES] =
          m_rta - pmv->d_gas_density_solvents[i_pl][POR_GAS_PRES] * ma_mw;
    }

    if (pd->e[pg->imtrx][R_POR_ENERGY]) {

      /* Add energy terms to water eqn senstivities */
      pmv->d_gas_density_solvents[i_pl][POR_TEMP] = drhosat_dT * mp->saturation;
      pmv->d_gas_density_solvents[i_pe][POR_LIQ_PRES] =
          pmv->enthalpy[1] * pmv->d_gas_density_solvents[i_pl][POR_LIQ_PRES] +
          pmv->enthalpy[2] * pmv->d_gas_density_solvents[i_pg][POR_LIQ_PRES];

      /* Add energy terms to air eqn sensitivies */
      if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
        pmv->d_gas_density_solvents[i_pg][POR_TEMP] =
            -Pg * m_rta / T - pmv->d_gas_density_solvents[i_pl][POR_TEMP] * ma_mw;
        pmv->d_gas_density_solvents[i_pe][POR_GAS_PRES] =
            pmv->gas_density_solvents[i_pl] * pmv->d_enthalpy[1][POR_GAS_PRES] +
            pmv->enthalpy[1] * pmv->d_gas_density_solvents[i_pl][POR_GAS_PRES] +
            pmv->gas_density_solvents[i_pg] * pmv->d_enthalpy[2][POR_GAS_PRES] +
            pmv->enthalpy[2] * pmv->d_gas_density_solvents[i_pg][POR_GAS_PRES] - 1.0;
      }

      pmv->d_gas_density_solvents[i_pe][POR_TEMP] =
          pmv->gas_density_solvents[i_pl] * pmv->d_enthalpy[1][POR_TEMP] +
          pmv->enthalpy[1] * pmv->d_gas_density_solvents[i_pl][POR_TEMP] +
          pmv->gas_density_solvents[i_pg] * pmv->d_enthalpy[2][POR_TEMP] +
          pmv->enthalpy[2] * pmv->d_gas_density_solvents[i_pg][POR_TEMP];
    }

    if (pd->TimeIntegration == TRANSIENT) {
      pmv_old->d_gas_density_solvents[i_pl][POR_LIQ_PRES] =
          rhosat_old * mp_old->d_saturation[POR_LIQ_PRES];
      mp_old->d_porous_vapor_pressure[i_pl][POR_LIQ_PRES] =
          pmv_old->d_gas_density_solvents[i_pl][POR_LIQ_PRES] / m_rtw_old;

      if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
        pmv_old->d_gas_density_solvents[i_pl][POR_GAS_PRES] =
            -pmv_old->d_gas_density_solvents[i_pl][POR_LIQ_PRES];

        mp_old->d_porous_vapor_pressure[i_pl][POR_GAS_PRES] =
            pmv_old->d_gas_density_solvents[i_pl][POR_GAS_PRES] / m_rtw_old;

        pmv_old->d_gas_density_solvents[i_pg][POR_GAS_PRES] =
            m_rta_old - pmv_old->d_gas_density_solvents[i_pl][POR_GAS_PRES] * ma_mw;

        pmv_old->d_gas_density_solvents[i_pg][POR_LIQ_PRES] =
            -pmv_old->d_gas_density_solvents[i_pl][POR_LIQ_PRES] * ma_mw;
      }

      if (pd->e[pg->imtrx][R_POR_ENERGY]) {

        /* Add energy terms to water eqn senstivities */
        pmv_old->d_gas_density_solvents[i_pl][POR_TEMP] = drhosat_dT_old * mp_old->saturation;
        pmv_old->d_gas_density_solvents[i_pe][POR_LIQ_PRES] =
            pmv_old->enthalpy[1] * pmv_old->d_gas_density_solvents[i_pl][POR_LIQ_PRES] +
            pmv_old->enthalpy[2] * pmv_old->d_gas_density_solvents[i_pg][POR_LIQ_PRES];

        /* Add energy terms to air eqn sensitivies */
        if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
          pmv_old->d_gas_density_solvents[i_pg][POR_TEMP] =
              -Pg_old * m_rta_old / T_old - pmv_old->d_gas_density_solvents[i_pl][POR_TEMP] * ma_mw;
          pmv_old->d_gas_density_solvents[i_pe][POR_GAS_PRES] =
              pmv_old->gas_density_solvents[i_pl] * pmv_old->d_enthalpy[1][POR_GAS_PRES] +
              pmv_old->enthalpy[1] * pmv_old->d_gas_density_solvents[i_pl][POR_GAS_PRES] +
              pmv_old->gas_density_solvents[i_pg] * pmv_old->d_enthalpy[2][POR_GAS_PRES] +
              pmv_old->enthalpy[2] * pmv_old->d_gas_density_solvents[i_pg][POR_GAS_PRES] - 1.0;
        }
        pmv_old->d_gas_density_solvents[i_pe][POR_TEMP] =
            pmv_old->gas_density_solvents[i_pl] * pmv_old->d_enthalpy[1][POR_TEMP] +
            pmv_old->enthalpy[1] * pmv_old->d_gas_density_solvents[i_pl][POR_TEMP] +
            pmv_old->gas_density_solvents[i_pg] * pmv_old->d_enthalpy[2][POR_TEMP] +
            pmv_old->enthalpy[2] * pmv_old->d_gas_density_solvents[i_pg][POR_TEMP];
      }
    }

    /* second derivatives for full Newton of gradients */

    if (af->Assemble_Jacobian) {
      pmv->d_d_gas_density_solvents[i_pl][POR_LIQ_PRES][POR_LIQ_PRES] =
          rho_sat * mp->d_d_saturation[POR_LIQ_PRES][POR_LIQ_PRES];

      if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
        pmv->d_d_gas_density_solvents[i_pg][POR_LIQ_PRES][POR_LIQ_PRES] =
            -pmv->d_d_gas_density_solvents[i_pl][POR_LIQ_PRES][POR_LIQ_PRES] * ma_mw;

        pmv->d_d_gas_density_solvents[i_pl][POR_LIQ_PRES][POR_GAS_PRES] =
            rho_sat * mp->d_d_saturation[POR_LIQ_PRES][POR_GAS_PRES];

        pmv->d_d_gas_density_solvents[i_pg][POR_LIQ_PRES][POR_GAS_PRES] =
            -pmv->d_d_gas_density_solvents[i_pl][POR_LIQ_PRES][POR_GAS_PRES] * ma_mw;

        pmv->d_d_gas_density_solvents[i_pg][POR_GAS_PRES][POR_LIQ_PRES] =
            -pmv->d_d_gas_density_solvents[i_pl][POR_GAS_PRES][POR_LIQ_PRES] * ma_mw;

        pmv->d_d_gas_density_solvents[i_pl][POR_GAS_PRES][POR_GAS_PRES] =
            rho_sat * mp->d_d_saturation[POR_GAS_PRES][POR_GAS_PRES];

        pmv->d_d_gas_density_solvents[i_pg][POR_GAS_PRES][POR_GAS_PRES] =
            -pmv->d_d_gas_density_solvents[i_pl][POR_GAS_PRES][POR_GAS_PRES] * ma_mw;
      }

      if (pd->e[pg->imtrx][R_POR_ENERGY]) {
        /* Add energy terms to water eqn senstivities */
        pmv->d_d_gas_density_solvents[i_pl][POR_LIQ_PRES][POR_TEMP] =
            rho_sat * mp->d_d_saturation[POR_LIQ_PRES][POR_TEMP] +
            drhosat_dT * mp->d_saturation[POR_LIQ_PRES];

        pmv->d_d_gas_density_solvents[i_pl][POR_TEMP][POR_LIQ_PRES] =
            drhosat_dT * mp->d_saturation[POR_LIQ_PRES];

        pmv->d_d_gas_density_solvents[i_pl][POR_TEMP][POR_TEMP] =
            drhosat_dT * mp->d_saturation[POR_TEMP] + d_drhosat_dT * mp->saturation;

        /* Add energy terms to air eqn sensitivies */
        if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {

          pmv->d_d_gas_density_solvents[i_pl][POR_GAS_PRES][POR_TEMP] =
              rho_sat * mp->d_d_saturation[POR_GAS_PRES][POR_TEMP] +
              drhosat_dT * mp->d_saturation[POR_GAS_PRES];

          pmv->d_d_gas_density_solvents[i_pl][POR_TEMP][POR_GAS_PRES] =
              drhosat_dT * mp->d_saturation[POR_GAS_PRES];

          pmv->d_d_gas_density_solvents[i_pg][POR_LIQ_PRES][POR_TEMP] =
              -pmv->d_d_gas_density_solvents[i_pl][POR_LIQ_PRES][POR_TEMP] * ma_mw;
          pmv->d_d_gas_density_solvents[i_pg][POR_GAS_PRES][POR_TEMP] =
              -pmv->d_d_gas_density_solvents[i_pl][POR_GAS_PRES][POR_TEMP] * ma_mw -
              1.0 / m_rta / T;

          pmv->d_d_gas_density_solvents[i_pg][POR_TEMP][POR_LIQ_PRES] =
              -pmv->d_d_gas_density_solvents[i_pl][POR_TEMP][POR_LIQ_PRES] * ma_mw;
          pmv->d_d_gas_density_solvents[i_pg][POR_TEMP][POR_GAS_PRES] =
              -pmv->d_d_gas_density_solvents[i_pl][POR_TEMP][POR_GAS_PRES] * ma_mw -
              1.0 / m_rta / T;
          pmv->d_d_gas_density_solvents[i_pg][POR_TEMP][POR_TEMP] =
              -pmv->d_d_gas_density_solvents[i_pl][POR_TEMP][POR_TEMP] * ma_mw +
              2.0 * Pg / m_rta / T / T;
        }

        /* Generalize the energy equation sensitivities to allow looping; only the
         * appropriate contributing terms are nonzero.
         */

        for (w1 = first_porous_var; w1 < POR_LAST + 1; w1++) {
          for (w2 = first_porous_var; w2 < POR_LAST + 1; w2++)
            pmv->d_d_gas_density_solvents[i_pe][w1][w2] =
                pmv->d_d_gas_density_solvents[i_pl][w1][w2] * pmv->enthalpy[1] +
                pmv->d_d_gas_density_solvents[i_pg][w1][w2] * pmv->enthalpy[2] +
                pmv->d_gas_density_solvents[i_pl][w1] * pmv->d_enthalpy[1][w2] +
                pmv->d_gas_density_solvents[i_pg][w1] * pmv->d_enthalpy[2][w2] +
                pmv->d_gas_density_solvents[i_pl][w2] * pmv->d_enthalpy[1][w1] +
                pmv->d_gas_density_solvents[i_pg][w2] * pmv->d_enthalpy[2][w1] +
                pmv->gas_density_solvents[i_pl] * pmv->d_d_enthalpy[1][w1][w2] +
                pmv->gas_density_solvents[i_pg] * pmv->d_d_enthalpy[2][w1][w2];
        }
      }
    }
  } else {
    GOMA_EH(GOMA_ERROR, "No models for gas phase conc.");
  }
  return;
} /*  end of   load_gas_conc_flat    */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void load_gas_conc_EOS(double porosity, double cap_pres, double d_cap_pres[2])

/*************************************************************************/
/*
 * load_gas_conc_EOS -- calculate gas phase concentrations in a porous
 *  media from the unknowns used in the problem.  USES curvature INDEPENDENT model.
 *  see load_gas_conc for KELVIN equivalent
 *
 *  input:   assume that load_fv and load_fv_grads and load_fv_mesh_derivs have all
 *  -----    been called already, so unknowns and their sensitivies are known at this
 *           gauss point
 *
 *  output:  calculates the gas phase concentrations and its first and second derivatives with
 *  -----    respect to all the problem unknowns
 */
{
  int w, w1, w2;
  double Pg, Pvsat, T;
  double rho_gas, drho_gas_dPg, drho_gas_dT;
  double d_drho_gas_dPg[2], d_drho_gas_dT[2];
  double rho_wg_sat, drho_wgsat_dT, d_drho_wgsat_dT;
  double Ywg, dYwg_dPg, dYwg_dT, d_dYwg_dPg[2], d_dYwg_dT[2];
  double Yag, dYag_dPg, dYag_dT, d_dYag_dPg[2], d_dYag_dT[2];
  double dPvsat_dT;

  const int i_pl = 0, i_pg = 1, i_pe = 3;
  int i, first_porous_var = POR_LAST - MAX_POROUS_NUM + 1;

  for (w = 0; w < MAX_PMV; w++) {
    pmv->gas_density_solvents[w] = 0.;
    mp->porous_vapor_pressure[w] = 0.;
  }
  /*
   *  Find which model is used for gas conc. and calculate:
   *    1)  the saturation, pmv->gas_density_solvents
   *    2)  the sensitivity of gas conc. to all concentrations, porosity and temperature
   *          i.e. the first derivative of gas conc. w.r.t. each variable
   *           put this in pmv->d_gas_density_solvents[var]
   *    3)  the second derivatives of gas conc. w.r.t. each variable, including
   *          cross-terms, put this in pmv->d_d_gas_density_solvents[var][var]
   *          second derivative is needed for sensitivity of fluxes which depend
   *          on the first derivative of gas conc.
   */
  if (mp->PorousVaporPressureModel[i_pl] == FLAT &&
      ((mp->PorousMediaType == POROUS_TWO_PHASE && mp->PorousGasConstantsModel) ||
       mp->PorousMediaType == POROUS_UNSATURATED)) {
    /* currently assume isothermal
     * FOR FLAT EQUATION
     *  mp->u_porous_vapor_pressure[i_pl][0] is the flat interface vapor pressure
     *  mp->u_porous_vapor_pressure[i_pl][1] is the liquid density
     *  mp->u_porous_vapor_pressure[i_pl][2] is the water molecular weight
     *  mp->u_porous_vapor_pressure[i_pl][3] is the gas law constant
     *  mp->u_porous_vapor_pressure[i_pl][4] is the temperature
     *
     * FOR IDEAL_GAS EQUATION --- for TWO_PHASE FLOW or UNSATURATED FLOW
     *  mp->u_porous_gas_constants[0] is the air molecular weight
     *  mp->u_porous_gas_constants[1] is the gas law constant
     *  mp->u_porous_gas_constants[2] is the temperature
     *  mp->u_porous_gas_constants[3] is the ambient base pressure of gas phase
     */

    /* If solving the gas/energy equation, replace constants in the above model
     * with values consistent with current solution variables
     */
    if (pd->e[pg->imtrx][R_POR_GAS_PRES])
      Pg = fv->p_gas;
    else
      Pg = mp->u_porous_gas_constants[3];
    if (Pg == 0.)
      return;
    if (pd->e[pg->imtrx][R_POR_ENERGY])
      T = fv->T;
    else
      T = mp->u_porous_vapor_pressure[i_pl][4] - 273.15;

    /* precalculate some parameters */

    rho_gas = calc_rho_gas(Pg, T, &Pvsat, &dPvsat_dT, &drho_gas_dPg, &drho_gas_dT, d_drho_gas_dPg,
                           d_drho_gas_dT, &rho_wg_sat, &drho_wgsat_dT, &d_drho_wgsat_dT);

    Ywg = calc_Ywg(rho_wg_sat, drho_wgsat_dT, d_drho_wgsat_dT, rho_gas, drho_gas_dPg, drho_gas_dT,
                   d_drho_gas_dPg, d_drho_gas_dT, &dYwg_dPg, &dYwg_dT, d_dYwg_dPg, d_dYwg_dT);

    /* Curvature independent model for VP and water conc in gas */
    /* PRS: I don't understand this!!!! Why set and then reset? Also occurs for i_pg below*/
    pmv->gas_density_solvents[i_pl] = Ywg * rho_gas;
    pmv->gas_density_solvents[i_pl] = Ywg * rho_gas * mp->saturation;
    /*      printf("\t Pg=%g, T=%g, Pvsat=%g, rho_gas=%g, sat=%g, Ywg=%g
       \n",Pg,T,Pvsat,rho_gas,mp->saturation,Ywg); printf("\t pmv->gas_density_solvents[i_pl]=%g
       \n", pmv->gas_density_solvents[i_pl]);   */
    mp->porous_vapor_pressure[i_pl] = Pvsat;

    /* For now, assume binary so Yag = 1 - Ywg */
    Yag = 1.0 - Ywg;
    dYag_dPg = -dYwg_dPg;
    dYag_dT = -dYwg_dT;
    for (i = 0; i < 2; i++) {
      d_dYag_dPg[i] = -d_dYwg_dPg[i];
      d_dYag_dT[i] = -d_dYwg_dT[i];
    }
    pmv->gas_density_solvents[i_pg] = Yag * rho_gas;
    pmv->gas_density_solvents[i_pg] = Yag * rho_gas * mp->saturation;

    /* Store non-zero gas density & mass fraction derivatives for use in
     * calculations of diffusive flux */

    pmv->rhog = rho_gas;
    pmv->d_rhog[i_pg] = drho_gas_dPg;
    pmv->d_rhog[i_pe] = drho_gas_dT;
    pmv->d_drhog[i_pg][i_pg] = d_drho_gas_dPg[0];
    pmv->d_drhog[i_pg][i_pe] = d_drho_gas_dPg[1];
    pmv->d_drhog[i_pe][i_pg] = d_drho_gas_dT[0];
    pmv->d_drhog[i_pe][i_pe] = d_drho_gas_dT[1];

    pmv->d_Ywg[i_pg] = dYwg_dPg;
    pmv->d_Ywg[i_pe] = dYwg_dT;
    pmv->d_dYwg[i_pg][i_pg] = d_dYwg_dPg[0];
    pmv->d_dYwg[i_pg][i_pe] = d_dYwg_dPg[1];
    pmv->d_dYwg[i_pe][i_pg] = d_dYwg_dT[0];
    pmv->d_dYwg[i_pe][i_pe] = d_dYwg_dT[1];

    pmv->d_Yag[i_pg] = dYag_dPg;
    pmv->d_Yag[i_pe] = dYag_dT;
    pmv->d_dYag[i_pg][i_pg] = d_dYag_dPg[0];
    pmv->d_dYag[i_pg][i_pe] = d_dYag_dPg[1];
    pmv->d_dYag[i_pe][i_pg] = d_dYag_dT[0];
    pmv->d_dYag[i_pe][i_pe] = d_dYag_dT[1];

    if (pd->e[pg->imtrx][R_POR_ENERGY])
      pmv->gas_density_solvents[i_pe] = pmv->gas_density_solvents[i_pl] * pmv->enthalpy[1] +
                                        pmv->gas_density_solvents[i_pg] * pmv->enthalpy[2];

    /* *** NOTE --- this model of vapor-liquid equilibrium has no dependence
     *              upon porosity, saturation, or pore-size
     */

    /* sensitivities */

    /* old way ...
    pmv->d_gas_density_solvents[i_pl][POR_LIQ_PRES] = 0.0;
    pmv->d_gas_density_solvents[i_pl][POR_LIQ_PRES] = Ywg * rho_gas *
    mp->d_saturation[POR_LIQ_PRES]; mp->d_porous_vapor_pressure[i_pl][POR_LIQ_PRES] = 0.0;


        pmv->d_gas_density_solvents[i_pl][POR_GAS_PRES] =
                                       dYwg_dPg*rho_gas + Ywg*drho_gas_dPg;
        mp->d_porous_vapor_pressure[i_pl][POR_GAS_PRES] = 0.0;

        pmv->d_gas_density_solvents[i_pg][POR_LIQ_PRES] = 0.0;
        pmv->d_gas_density_solvents[i_pg][POR_GAS_PRES] =
                                       dYag_dPg*rho_gas + Yag*drho_gas_dPg;

    if (pd->e[pg->imtrx][R_POR_ENERGY])
      {


        pmv->d_gas_density_solvents[i_pl][POR_TEMP] =
                                        dYwg_dT*rho_gas + Ywg*drho_gas_dT;
        pmv->d_gas_density_solvents[i_pe][POR_LIQ_PRES] = 0.0;


            pmv->d_gas_density_solvents[i_pg][POR_TEMP] =
                                        dYag_dT*rho_gas + Yag*drho_gas_dT;
            pmv->d_gas_density_solvents[i_pe][POR_GAS_PRES] =
               pmv->enthalpy[1]*pmv->d_gas_density_solvents[i_pl][POR_GAS_PRES] +
               pmv->gas_density_solvents[i_pl]*pmv->d_enthalpy[1][POR_GAS_PRES] +
               pmv->enthalpy[2]*pmv->d_gas_density_solvents[i_pg][POR_GAS_PRES] +
               pmv->gas_density_solvents[i_pg]*pmv->d_enthalpy[2][POR_GAS_PRES];

        pmv->d_gas_density_solvents[i_pe][POR_TEMP] =
             pmv->enthalpy[1]*pmv->d_gas_density_solvents[i_pl][POR_TEMP] +
             pmv->gas_density_solvents[i_pl]*pmv->d_enthalpy[1][POR_TEMP] +
             pmv->enthalpy[2]*pmv->d_gas_density_solvents[i_pg][POR_TEMP] +
             pmv->gas_density_solvents[i_pg]*pmv->d_enthalpy[2][POR_TEMP];
      }
      .... end of old way */

    /* new way */
    /* Water equation */

    for (w1 = 0; w1 < MAX_PMV; w1++)
      pmv->d_gas_density_solvents[i_pl][POR_LIQ_PRES + w1] =
          pmv->d_Ywg[w1] * rho_gas * mp->saturation + Ywg * pmv->d_rhog[w1] * mp->saturation +
          Ywg * rho_gas * mp->d_saturation[POR_LIQ_PRES + w1];

    /* Air equation */

    for (w1 = 0; w1 < MAX_PMV; w1++)
      pmv->d_gas_density_solvents[i_pg][POR_LIQ_PRES + w1] =
          pmv->d_Yag[w1] * rho_gas * mp->saturation + Yag * pmv->d_rhog[w1] * mp->saturation +
          Yag * rho_gas * mp->d_saturation[POR_LIQ_PRES + w1];

    /* Energy equation */
    if (pd->e[pg->imtrx][R_POR_ENERGY]) {
      for (w1 = 0; w1 < MAX_PMV; w1++)
        pmv->d_gas_density_solvents[i_pe][POR_LIQ_PRES + w1] =
            pmv->enthalpy[1] * pmv->d_gas_density_solvents[i_pl][POR_LIQ_PRES + w1] +
            pmv->gas_density_solvents[i_pl] * pmv->d_enthalpy[1][POR_LIQ_PRES + w1] +
            pmv->enthalpy[2] * pmv->d_gas_density_solvents[i_pg][POR_LIQ_PRES + w1] +
            pmv->gas_density_solvents[i_pg] * pmv->d_enthalpy[2][POR_LIQ_PRES + w1];
    }

    /* second derivatives for full Newton of gradients */

    if (af->Assemble_Jacobian) {
      /* old way ...
      pmv->d_d_gas_density_solvents[i_pl][POR_LIQ_PRES][POR_LIQ_PRES] =
                 Ywg * rho_gas * mp->d_d_saturation[POR_LIQ_PRES][POR_LIQ_PRES];
          pmv->d_d_gas_density_solvents[i_pg][POR_LIQ_PRES][POR_LIQ_PRES] = 0.0;
          pmv->d_d_gas_density_solvents[i_pl][POR_LIQ_PRES][POR_GAS_PRES] = 0.0;
          pmv->d_d_gas_density_solvents[i_pg][POR_LIQ_PRES][POR_GAS_PRES] = 0.0;
          pmv->d_d_gas_density_solvents[i_pg][POR_GAS_PRES][POR_LIQ_PRES] = 0.0;

          pmv->d_d_gas_density_solvents[i_pl][POR_GAS_PRES][POR_GAS_PRES] =
               d_dYwg_dPg[0]*rho_gas + dYwg_dPg*drho_gas_dPg +
               dYwg_dPg*drho_gas_dPg + Ywg*d_drho_gas_dPg[0];

          pmv->d_d_gas_density_solvents[i_pg][POR_GAS_PRES][POR_GAS_PRES] =
               d_dYag_dPg[0]*rho_gas + dYag_dPg*drho_gas_dPg +
               dYag_dPg*drho_gas_dPg + Yag*d_drho_gas_dPg[0];


      if (pd->e[pg->imtrx][R_POR_ENERGY])
        {
          pmv->d_d_gas_density_solvents[i_pl][POR_LIQ_PRES][POR_TEMP] = 0.0;
          pmv->d_d_gas_density_solvents[i_pl][POR_TEMP][POR_LIQ_PRES] = 0.0;

          pmv->d_d_gas_density_solvents[i_pl][POR_TEMP][POR_TEMP] =
               d_dYwg_dT[1]*rho_gas + dYwg_dT*drho_gas_dT +
               dYwg_dT*drho_gas_dT  + Ywg*d_drho_gas_dT[1];


              pmv->d_d_gas_density_solvents[i_pl][POR_GAS_PRES][POR_TEMP] =
                   d_dYwg_dPg[1]*rho_gas + dYwg_dPg*drho_gas_dT +
                   dYwg_dT*drho_gas_dPg  + Ywg*d_drho_gas_dPg[1];

              pmv->d_d_gas_density_solvents[i_pl][POR_TEMP][POR_GAS_PRES] =
                   d_dYwg_dT[0]*rho_gas  + dYwg_dT*drho_gas_dPg +
                   dYwg_dPg*drho_gas_dT  + Ywg*d_drho_gas_dT[0];

              pmv->d_d_gas_density_solvents[i_pg][POR_LIQ_PRES][POR_TEMP] = 0.0;

              pmv->d_d_gas_density_solvents[i_pg][POR_GAS_PRES][POR_TEMP] =
                   d_dYag_dPg[1]*rho_gas + dYag_dPg*drho_gas_dT +
                   dYag_dT*drho_gas_dPg  + Yag*d_drho_gas_dPg[1];


              pmv->d_d_gas_density_solvents[i_pg][POR_TEMP][POR_LIQ_PRES] = 0.0;

              pmv->d_d_gas_density_solvents[i_pg][POR_TEMP][POR_GAS_PRES] =
                   d_dYag_dT[0]*rho_gas  + dYag_dT*drho_gas_dPg +
                   dYag_dPg*drho_gas_dT  + Yag*d_drho_gas_dT[0];

              pmv->d_d_gas_density_solvents[i_pg][POR_TEMP][POR_TEMP] =
                   d_dYag_dT[1]*rho_gas + dYag_dT*drho_gas_dT +
                   dYag_dT*drho_gas_dT  + Yag*d_drho_gas_dT[1];



      for (w1 = first_porous_var; w1 < POR_LAST+1; w1++) {
        for (w2 = first_porous_var; w2 < POR_LAST+1; w2++)
         pmv->d_d_gas_density_solvents[i_pe][w1][w2] =
           pmv->d_d_gas_density_solvents[i_pl][w1][w2]*pmv->enthalpy[1] +
           pmv->d_d_gas_density_solvents[i_pg][w1][w2]*pmv->enthalpy[2] +
           pmv->d_gas_density_solvents[i_pl][w1] * pmv->d_enthalpy[1][w2] +
           pmv->d_gas_density_solvents[i_pg][w1] * pmv->d_enthalpy[2][w2] +
           pmv->d_gas_density_solvents[i_pl][w2] * pmv->d_enthalpy[1][w1] +
           pmv->d_gas_density_solvents[i_pg][w2] * pmv->d_enthalpy[2][w1] +
           pmv->gas_density_solvents[i_pl] * pmv->d_d_enthalpy[1][w1][w2] +
           pmv->gas_density_solvents[i_pg] * pmv->d_d_enthalpy[2][w1][w2];
      }

  }
 .... end of old way */
      /* new way */

      /* Water equation */

      for (w1 = 0; w1 < MAX_PMV; w1++) {
        for (w2 = 0; w2 < MAX_PMV; w2++)
          pmv->d_d_gas_density_solvents[i_pl][POR_LIQ_PRES + w1][POR_LIQ_PRES + w2] =
              pmv->d_dYwg[w1][w2] * rho_gas * mp->saturation +
              pmv->d_Ywg[w1] * pmv->d_rhog[w2] * mp->saturation +
              Ywg * rho_gas * mp->d_saturation[POR_LIQ_PRES + w2] +
              pmv->d_Ywg[w2] * pmv->d_rhog[w1] * mp->saturation +
              Ywg * pmv->d_drhog[w1][w2] * mp->saturation +
              Ywg * pmv->d_rhog[w1] * mp->d_saturation[POR_LIQ_PRES + w2] +
              pmv->d_Ywg[w2] * rho_gas * mp->d_saturation[POR_LIQ_PRES + w1] +
              Ywg * pmv->d_rhog[w2] * mp->d_saturation[POR_LIQ_PRES + w1] +
              Ywg * rho_gas * mp->d_d_saturation[POR_LIQ_PRES + w1][POR_LIQ_PRES + w2];
      }

      /* Air equation */

      for (w1 = 0; w1 < MAX_PMV; w1++) {
        for (w2 = 0; w2 < MAX_PMV; w2++)
          pmv->d_d_gas_density_solvents[i_pl][POR_LIQ_PRES + w1][POR_LIQ_PRES + w2] =
              pmv->d_dYag[w1][w2] * rho_gas * mp->saturation +
              pmv->d_Yag[w1] * pmv->d_rhog[w2] * mp->saturation +
              Yag * rho_gas * mp->d_saturation[POR_LIQ_PRES + w2] +
              pmv->d_Yag[w2] * pmv->d_rhog[w1] * mp->saturation +
              Yag * pmv->d_drhog[w1][w2] * mp->saturation +
              Yag * pmv->d_rhog[w1] * mp->d_saturation[POR_LIQ_PRES + w2] +
              pmv->d_Yag[w2] * rho_gas * mp->d_saturation[POR_LIQ_PRES + w1] +
              Yag * pmv->d_rhog[w2] * mp->d_saturation[POR_LIQ_PRES + w1] +
              Yag * rho_gas * mp->d_d_saturation[POR_LIQ_PRES + w1][POR_LIQ_PRES + w2];
      }

      /* Energy equation */

      if (pd->e[pg->imtrx][R_POR_ENERGY]) {
        for (w1 = first_porous_var; w1 < POR_LAST + 1; w1++) {
          for (w2 = first_porous_var; w2 < POR_LAST + 1; w2++)
            pmv->d_d_gas_density_solvents[i_pe][w1][w2] =
                pmv->d_d_gas_density_solvents[i_pl][w1][w2] * pmv->enthalpy[1] +
                pmv->d_d_gas_density_solvents[i_pg][w1][w2] * pmv->enthalpy[2] +
                pmv->d_gas_density_solvents[i_pl][w1] * pmv->d_enthalpy[1][w2] +
                pmv->d_gas_density_solvents[i_pg][w1] * pmv->d_enthalpy[2][w2] +
                pmv->d_gas_density_solvents[i_pl][w2] * pmv->d_enthalpy[1][w1] +
                pmv->d_gas_density_solvents[i_pg][w2] * pmv->d_enthalpy[2][w1] +
                pmv->gas_density_solvents[i_pl] * pmv->d_d_enthalpy[1][w1][w2] +
                pmv->gas_density_solvents[i_pg] * pmv->d_d_enthalpy[2][w1][w2];
        }
      }

    } /* if Assemble Jacobian */

    if (pd->TimeIntegration == TRANSIENT) {

      if (pd->e[pg->imtrx][R_POR_GAS_PRES])
        Pg = fv_old->p_gas;
      else
        Pg = mp->u_porous_gas_constants[3];
      if (pd->e[pg->imtrx][R_POR_ENERGY])
        T = fv_old->T;
      else
        T = mp->u_porous_vapor_pressure[i_pl][4] - 273.15;

      /* precalculate some parameters */

      rho_gas = calc_rho_gas(Pg, T, &Pvsat, &dPvsat_dT, &drho_gas_dPg, &drho_gas_dT, d_drho_gas_dPg,
                             d_drho_gas_dT, &rho_wg_sat, &drho_wgsat_dT, &d_drho_wgsat_dT);
      Ywg = calc_Ywg(rho_wg_sat, drho_wgsat_dT, d_drho_wgsat_dT, rho_gas, drho_gas_dPg, drho_gas_dT,
                     d_drho_gas_dPg, d_drho_gas_dT, &dYwg_dPg, &dYwg_dT, d_dYwg_dPg, d_dYwg_dT);

      /* Curvature independent model for VP and water conc in gas */
      pmv_old->gas_density_solvents[i_pl] = Ywg * rho_gas;
      pmv_old->gas_density_solvents[i_pl] = Ywg * rho_gas * mp_old->saturation;
      mp_old->porous_vapor_pressure[i_pl] = Pvsat;

      /* For now, assume binary so Yag = 1 - Ywg */
      Yag = 1.0 - Ywg;
      dYag_dPg = -dYwg_dPg;
      dYag_dT = -dYag_dT;
      for (i = 0; i < 2; i++) {
        d_dYag_dPg[i] = -d_dYwg_dPg[i];
        d_dYag_dT[i] = -d_dYwg_dT[i];
      }

      pmv_old->gas_density_solvents[i_pg] = Yag * rho_gas * mp_old->saturation;

      /* Store non-zero gas density & mass fraction derivatives for use in
       * calculations of diffusive flux */

      pmv_old->rhog = rho_gas;
      pmv_old->d_rhog[i_pg] = drho_gas_dPg;
      pmv_old->d_rhog[i_pe] = drho_gas_dT;
      pmv_old->d_drhog[i_pg][i_pg] = d_drho_gas_dPg[0];
      pmv_old->d_drhog[i_pg][i_pe] = d_drho_gas_dPg[1];
      pmv_old->d_drhog[i_pe][i_pg] = d_drho_gas_dT[0];
      pmv_old->d_drhog[i_pe][i_pe] = d_drho_gas_dT[1];

      pmv_old->d_Ywg[i_pg] = dYwg_dPg;
      pmv_old->d_Ywg[i_pe] = dYwg_dT;
      pmv_old->d_dYwg[i_pg][i_pg] = d_dYwg_dPg[0];
      pmv_old->d_dYwg[i_pg][i_pe] = d_dYwg_dPg[1];
      pmv_old->d_dYwg[i_pe][i_pg] = d_dYwg_dT[0];
      pmv_old->d_dYwg[i_pe][i_pe] = d_dYwg_dT[1];

      pmv_old->d_Yag[i_pg] = dYag_dPg;
      pmv_old->d_Yag[i_pe] = dYag_dT;
      pmv_old->d_dYag[i_pg][i_pg] = d_dYag_dPg[0];
      pmv_old->d_dYag[i_pg][i_pe] = d_dYag_dPg[1];
      pmv_old->d_dYag[i_pe][i_pg] = d_dYag_dT[0];
      pmv_old->d_dYag[i_pe][i_pe] = d_dYag_dT[1];

      if (pd->e[pg->imtrx][R_POR_ENERGY])
        pmv_old->gas_density_solvents[i_pe] =
            pmv_old->gas_density_solvents[i_pl] * pmv_old->enthalpy[1] +
            pmv_old->gas_density_solvents[i_pg] * pmv_old->enthalpy[2];

      /* *** NOTE --- this model of vapor-liquid equilibrium has no dependence
       *              upon porosity, saturation, or pore-size
       */

      /* sensitivities */
      /* old way ...
      pmv_old->d_gas_density_solvents[i_pl][POR_LIQ_PRES] = 0.0;
      pmv_old->d_gas_density_solvents[i_pl][POR_LIQ_PRES] =
                                         Ywg * rho_gas * mp_old->d_saturation[POR_LIQ_PRES];
      mp_old->d_porous_vapor_pressure[i_pl][POR_LIQ_PRES] = 0.0;

          pmv_old->d_gas_density_solvents[i_pl][POR_GAS_PRES] =
                                         dYwg_dPg*rho_gas + Ywg*drho_gas_dPg;
          mp_old->d_porous_vapor_pressure[i_pl][POR_GAS_PRES] = 0.0;

          pmv_old->d_gas_density_solvents[i_pg][POR_LIQ_PRES] = 0.0;
          pmv_old->d_gas_density_solvents[i_pg][POR_GAS_PRES] =
                                         dYag_dPg*rho_gas + Yag*drho_gas_dPg;

      if (pd->e[pg->imtrx][R_POR_ENERGY])
        {


           pmv_old->d_gas_density_solvents[i_pl][POR_TEMP] =
                                          dYwg_dT*rho_gas + Ywg*drho_gas_dT;
          pmv_old->d_gas_density_solvents[i_pe][POR_LIQ_PRES] = 0.0;


              pmv_old->d_gas_density_solvents[i_pg][POR_TEMP] =
                                          dYag_dT*rho_gas + Yag*drho_gas_dT;
              pmv_old->d_gas_density_solvents[i_pe][POR_GAS_PRES] =
                 pmv_old->enthalpy[1]*pmv_old->d_gas_density_solvents[i_pl][POR_GAS_PRES] +
                 pmv_old->gas_density_solvents[i_pl]*pmv_old->d_enthalpy[1][POR_GAS_PRES] +
                 pmv_old->enthalpy[2]*pmv_old->d_gas_density_solvents[i_pg][POR_GAS_PRES] +
                 pmv_old->gas_density_solvents[i_pg]*pmv_old->d_enthalpy[2][POR_GAS_PRES];


          pmv_old->d_gas_density_solvents[i_pe][POR_TEMP] =
               pmv_old->enthalpy[1]*pmv_old->d_gas_density_solvents[i_pl][POR_TEMP] +
               pmv_old->gas_density_solvents[i_pl]*pmv_old->d_enthalpy[1][POR_TEMP] +
               pmv_old->enthalpy[2]*pmv_old->d_gas_density_solvents[i_pg][POR_TEMP] +
               pmv_old->gas_density_solvents[i_pg]*pmv_old->d_enthalpy[2][POR_TEMP];
        }
        ... end of old way */
      /* new way */
      /* Water equation */

      for (w1 = 0; w1 < MAX_PMV; w1++)
        pmv_old->d_gas_density_solvents[i_pl][POR_LIQ_PRES + w1] =
            pmv_old->d_Ywg[w1] * rho_gas * mp_old->saturation +
            Ywg * pmv_old->d_rhog[w1] * mp_old->saturation +
            Ywg * rho_gas * mp_old->d_saturation[POR_LIQ_PRES + w1];

      /* Air equation */

      for (w1 = 0; w1 < MAX_PMV; w1++)
        pmv_old->d_gas_density_solvents[i_pg][POR_LIQ_PRES + w1] =
            pmv_old->d_Yag[w1] * rho_gas * mp_old->saturation +
            Yag * pmv_old->d_rhog[w1] * mp_old->saturation +
            Yag * rho_gas * mp_old->d_saturation[POR_LIQ_PRES + w1];

      /* Energy equation */
      if (pd->e[pg->imtrx][R_POR_ENERGY]) {
        for (w1 = 0; w1 < MAX_PMV; w1++)
          pmv_old->d_gas_density_solvents[i_pe][POR_LIQ_PRES + w1] =
              pmv_old->enthalpy[1] * pmv_old->d_gas_density_solvents[i_pl][POR_LIQ_PRES + w1] +
              pmv_old->gas_density_solvents[i_pl] * pmv_old->d_enthalpy[1][POR_LIQ_PRES + w1] +
              pmv_old->enthalpy[2] * pmv_old->d_gas_density_solvents[i_pg][POR_LIQ_PRES + w1] +
              pmv_old->gas_density_solvents[i_pg] * pmv_old->d_enthalpy[2][POR_LIQ_PRES + w1];
      }

    } /* if TRANSIENT */

  } else {
    GOMA_EH(GOMA_ERROR, "No models for gas phase conc.");
  }
  return;
} /*  end of   load_gas_conc_EOS    */
/*****************************************************************************/

void load_bulk_density(double porosity, double cap_pres, double saturation, double d_cap_pres[2])

/************************************************************************
 * load_bulk_density -- calculate bulk phase densities in a porous
 *                      media from the unknowns used in the problem
 *
 *  input:
 *  -----
 *
 *
 *  output:
 *  -----
 *
 *   calculates the bulk phase concentrations and its first and
 *   second derivatives with  respect to all the problem unknowns
 ************************************************************************/
{
  int w, w1, w2;
  int eqn, var, var1, var2;
  int i_pl = 0, i_pg = 1, i_pe = 3;
  /* the first calculation really is to get the current liquid density */
  dbl liquid_density = 0.0, d_liquid_density_dpliq = 0.0, liquid_density_old = 0.0,
      d_liquid_density_dpliq_old = 0.0;
  dbl rhoCp = 0.0;
  /*
   *  Find which model is used for bulk conc. and calculate:
   *    1)  the saturation, pmv->bulk_density_solvents
   *    2)  the sensitivity of bulk conc. to all concentrations, porosity and temperature
   *          i.e. the first derivative of bulk conc. w.r.t. each variable
   *           put this in pmv->d_bulk_density_solvents[var]
   *    3)  the second derivatives of bulk conc. w.r.t. each variable, including
   *          cross-terms, put this in pmv->d_d_bulk_density_solvents[var][var]
   *          second derivative is needed for sensitivity of fluxes which depend
   *          on the first derivative of bulk conc.
   */

  /*
   * The first calculation really is to get the current liquid density.
   * We assume a linear compressibility relation here.
   */
  liquid_density = mp->density * (1 + mp->PorousLiqCompress * (fv->p_liq - mp->PorousLiqRefPress));
  d_liquid_density_dpliq = mp->density * mp->PorousLiqCompress;
  d_liquid_density_dpliq_old = d_liquid_density_dpliq;

  /* N.B. this first calculation should just make pmv->bulk_density for the solid
   * phase equal to zero.  Also note that this assumes liquid is single component or that
   * "solvents" are cumulative
   *
   * If the energy equation is being solved, there is a contribution from the
   * solid phase.
   */

  if (pd->TimeIntegration == TRANSIENT) {
    liquid_density_old =
        mp->density * (1 + mp->PorousLiqCompress * (fv_old->p_liq - mp->PorousLiqRefPress));
  }

  for (w = 0; w < MAX_PMV; w++) {
    /*short circuit if the equation is not active */
    if (pd->e[pg->imtrx][R_POR_LIQ_PRES + w]) {
      pmv->bulk_density[w] = pmv->gas_density_solvents[w] * porosity * (1. - saturation) +
                             liquid_density * pmv->liq_Xvol_solvents[w] * porosity * saturation;

      if (pd->TimeIntegration == TRANSIENT) {
        pmv_old->bulk_density[w] =
            pmv_old->gas_density_solvents[w] * mp_old->porosity * (1. - mp_old->saturation) +
            liquid_density_old * mp_old->porosity * mp_old->saturation *
                pmv_old->liq_Xvol_solvents[w];
      }
    }
  }

  if (pd->e[pg->imtrx][R_SHELL_SAT_OPEN] || pd->e[pg->imtrx][R_SHELL_SAT_OPEN_2]) {
    pmv->bulk_density[0] = liquid_density * pmv->liq_Xvol_solvents[0] * porosity * saturation;
  }

  if (pd->e[pg->imtrx][R_POR_ENERGY]) {
    /* matrix density and specific heat are constants */
    rhoCp = mp->matrix_density * mp->specific_heat;
    pmv->bulk_density[i_pe] += (1.0 - porosity) * rhoCp * fv->T;
    if (pd->TimeIntegration == TRANSIENT)
      pmv_old->bulk_density[i_pe] += (1.0 - mp_old->porosity) * rhoCp * fv_old->T;
  }

  /*
   * Calculate the first order derivatives of the inventory terms
   * wrt the unknowns in the problem
   */
  for (w = 0; w < MAX_PMV; w++) {

    eqn = R_POR_LIQ_PRES + w;
    /*
     * First, let's calculate the dependence on the liquid phase
     * pressure
     */

    var = POR_LIQ_PRES;

    if (pd->e[pg->imtrx][eqn] && pd->v[pg->imtrx][var]) {
      pmv->d_bulk_density[w][POR_LIQ_PRES] =
          pmv->d_gas_density_solvents[w][POR_LIQ_PRES] * porosity * (1. - saturation);
      pmv->d_bulk_density[w][POR_LIQ_PRES] +=
          mp->d_saturation[POR_LIQ_PRES] * porosity *
          (liquid_density * pmv->liq_Xvol_solvents[w] - pmv->gas_density_solvents[w]);
      pmv->d_bulk_density[w][POR_LIQ_PRES] +=
          d_liquid_density_dpliq * porosity * saturation * pmv->liq_Xvol_solvents[w];
      pmv->d_bulk_density[w][POR_LIQ_PRES] +=
          mp->d_porosity[POR_LIQ_PRES] * (pmv->gas_density_solvents[w] * (1. - saturation) +
                                          liquid_density * saturation * pmv->liq_Xvol_solvents[w]);

      /* variation of mass fraction of component w in the liquid phase:
       *  d_liq_Xvol_solvents currently zero for all but energy eqn
       */
      pmv->d_bulk_density[w][POR_LIQ_PRES] +=
          pmv->d_liq_Xvol_solvents[w][POR_LIQ_PRES] * porosity * liquid_density * saturation;
    }

    /*
     * Calculate the dependence on the Gas phase pressure
     */

    var = POR_GAS_PRES;

    if (pd->e[pg->imtrx][eqn] && pd->v[pg->imtrx][var]) {
      pmv->d_bulk_density[w][POR_GAS_PRES] =
          pmv->d_gas_density_solvents[w][POR_GAS_PRES] * porosity * (1. - saturation);
      pmv->d_bulk_density[w][POR_GAS_PRES] +=
          mp->d_saturation[POR_GAS_PRES] * porosity *
          (liquid_density * pmv->liq_Xvol_solvents[w] - pmv->gas_density_solvents[w]);

      pmv->d_bulk_density[w][POR_GAS_PRES] +=
          mp->d_porosity[POR_GAS_PRES] * (pmv->gas_density_solvents[w] * (1. - saturation) +
                                          liquid_density * saturation * pmv->liq_Xvol_solvents[w]);

      /* variation of mass fraction of component w in the liquid phase:
       *  d_liq_Xvol_solvents currently zero for all but energy eqn
       */
      pmv->d_bulk_density[w][POR_GAS_PRES] +=
          pmv->d_liq_Xvol_solvents[w][POR_GAS_PRES] * porosity * liquid_density * saturation;
    }

    /*
     * Calculate the dependence on the porosity
     */
    var = POR_POROSITY;

    if (pd->e[pg->imtrx][eqn] && pd->v[pg->imtrx][var]) {

      pmv->d_bulk_density[w][POR_POROSITY] =
          pmv->d_gas_density_solvents[w][POR_POROSITY] * porosity * (1. - saturation);
      pmv->d_bulk_density[w][POR_POROSITY] -=
          pmv->gas_density_solvents[w] * porosity * mp->d_saturation[POR_POROSITY];
      pmv->d_bulk_density[w][POR_POROSITY] +=
          liquid_density * porosity * mp->d_saturation[POR_POROSITY] * pmv->liq_Xvol_solvents[w];
      pmv->d_bulk_density[w][POR_POROSITY] +=
          mp->d_porosity[POR_POROSITY] * (pmv->gas_density_solvents[w] * (1. - saturation) +
                                          liquid_density * saturation * pmv->liq_Xvol_solvents[w]);

      /* variation of mass fraction of component w in the liquid phase:
       *  d_liq_Xvol_solvents currently zero for all but energy eqn
       */
      pmv->d_bulk_density[w][POR_POROSITY] +=
          pmv->d_liq_Xvol_solvents[w][POR_POROSITY] * porosity * liquid_density * saturation;
    }
    /*
     * Calculate the dependence on the temperature
     */

    var = POR_TEMP;

    if (pd->e[pg->imtrx][eqn] && pd->v[pg->imtrx][var]) {
      pmv->d_bulk_density[w][POR_TEMP] =
          pmv->d_gas_density_solvents[w][POR_TEMP] * porosity * (1. - saturation);
      /* saturation not a function of temperature, d_saturation[POR_TEMP] = 0 */
      pmv->d_bulk_density[w][POR_TEMP] +=
          mp->d_saturation[POR_TEMP] * porosity *
          (liquid_density * pmv->liq_Xvol_solvents[w] - pmv->gas_density_solvents[w]);

      /* porosity not a function of temperature, d_porosity[POR_TEMP] = 0 */
      pmv->d_bulk_density[w][POR_TEMP] +=
          mp->d_porosity[POR_TEMP] * (pmv->gas_density_solvents[w] * (1. - saturation) +
                                      liquid_density * saturation * pmv->liq_Xvol_solvents[w]);

      /* variation of mass fraction of component w in the liquid phase:
       *  d_liq_Xvol_solvents currently zero for all but energy eqn
       */
      pmv->d_bulk_density[w][POR_TEMP] +=
          pmv->d_liq_Xvol_solvents[w][POR_TEMP] * porosity * liquid_density * saturation;
    }
  }
  if (pd->e[pg->imtrx][R_POR_ENERGY]) {
    for (w1 = 0; w1 < MAX_PMV; w1++) {
      pmv->d_bulk_density[i_pe][POR_LIQ_PRES + w1] += -mp->d_porosity[w1] * rhoCp * fv->T;
    }

    pmv->d_bulk_density[i_pe][POR_TEMP] += (1.0 - porosity) * rhoCp;
  }

  if (pd->TimeIntegration == TRANSIENT) {
    /* N.B. PRS: I am trying these loops just a little differently then the ones above,
     * by not unrolling them.   Note that this assumes we are always looping over all
     * porous equations, with those inactive just getting zeros.
     */
    for (w = 0; w < MAX_PMV; w++) {
      for (w1 = 0; w1 < MAX_PMV; w1++) {
        eqn = R_POR_LIQ_PRES + w;
        var = POR_LIQ_PRES + w1;
        if (pd->e[pg->imtrx][eqn] && pd->v[pg->imtrx][var]) {
          pmv_old->d_bulk_density[w][var] = pmv_old->d_gas_density_solvents[w][var] *
                                            mp_old->porosity * (1. - mp_old->saturation);
          pmv_old->d_bulk_density[w][var] -=
              pmv_old->gas_density_solvents[w] * mp_old->porosity * mp_old->d_saturation[var];
          pmv_old->d_bulk_density[w][var] += liquid_density_old * mp_old->porosity *
                                             mp_old->d_saturation[var] *
                                             pmv_old->liq_Xvol_solvents[w];
          pmv_old->d_bulk_density[w][var] += d_liquid_density_dpliq_old * mp_old->porosity *
                                             mp_old->saturation * pmv_old->liq_Xvol_solvents[w];
          pmv_old->d_bulk_density[w][var] +=
              mp_old->d_porosity[var] *
              (pmv_old->gas_density_solvents[w] * (1. - mp_old->saturation) +
               liquid_density_old * mp_old->saturation * pmv_old->liq_Xvol_solvents[w]);
          /* variation of mass fraction of component w in the liquid phase:
           *  d_liq_Xvol_solvents currently zero for all but energy eqn
           */
          pmv_old->d_bulk_density[w][var] += pmv_old->d_liq_Xvol_solvents[w][var] *
                                             mp_old->porosity * liquid_density_old *
                                             mp_old->saturation;
        }
      }
    }
    /* Additional "solid" terms for energy equation */

    if (pd->e[pg->imtrx][R_POR_ENERGY]) {
      for (w1 = 0; w1 < MAX_PMV; w1++) {
        pmv_old->d_bulk_density[i_pe][POR_LIQ_PRES + w1] +=
            -mp_old->d_porosity[w1] * rhoCp * fv_old->T;
      }
      pmv_old->d_bulk_density[i_pe][POR_TEMP] += (1.0 - mp_old->porosity) * rhoCp;
    }
  } /* if Transient */

  if (af->Assemble_Jacobian) {
    /* assume that gas (species i_pg) doesn't dissolve in liquid, so no sensitivity yet
     * to liquid_volume_fraction, unless energy equation where the liquid_volume_fraction
     * "holds" the liquid enthalpy; so include senstivity for looping convenience and
     * possible future dependencies; noncontributing derivatives are zero.
     *
     * Also, the liquid density is only a function of liquid pressure (currently), but
     * fill a temporary vector for looping convenience and possible future dependencies.
     * Use the outdated d_density already present in mm_mp_structs.h */
    for (w = 0; w < MAX_PMV; w++)
      mp->d_density[w] = 0.0;
    mp->d_density[POR_LIQ_PRES] = d_liquid_density_dpliq;

    for (w = 0; w < MAX_PMV; w++) {
      for (w1 = 0; w1 < MAX_PMV; w1++) {
        for (w2 = 0; w2 < MAX_PMV; w2++) {

          eqn = R_POR_LIQ_PRES + w;
          var1 = POR_LIQ_PRES + w1;
          var2 = POR_LIQ_PRES + w2;

          if (pd->e[pg->imtrx][eqn] && pd->v[pg->imtrx][var1] && pd->v[pg->imtrx][var2]) {
            pmv->d_d_bulk_density[w][POR_LIQ_PRES + w1][POR_LIQ_PRES + w2] =
                porosity *
                (pmv->d_d_gas_density_solvents[w][POR_LIQ_PRES + w1][POR_LIQ_PRES + w2] *
                     (1. - saturation) -
                 pmv->d_gas_density_solvents[w][POR_LIQ_PRES + w1] *
                     mp->d_saturation[POR_LIQ_PRES + w2] -
                 pmv->d_gas_density_solvents[w][POR_LIQ_PRES + w2] *
                     mp->d_saturation[POR_LIQ_PRES + w1] -
                 pmv->gas_density_solvents[w] *
                     mp->d_d_saturation[POR_LIQ_PRES + w1][POR_LIQ_PRES + w2] +
                 liquid_density *
                     (mp->d_d_saturation[POR_LIQ_PRES + w1][POR_LIQ_PRES + w2] *
                          pmv->liq_Xvol_solvents[w] +
                      mp->d_saturation[POR_LIQ_PRES + w1] *
                          pmv->d_liq_Xvol_solvents[w][POR_LIQ_PRES + w2] +
                      mp->d_saturation[POR_LIQ_PRES + w2] *
                          pmv->d_liq_Xvol_solvents[w][POR_LIQ_PRES + w1] +
                      saturation *
                          pmv->d_d_liq_Xvol_solvents[w][POR_LIQ_PRES + w1][POR_LIQ_PRES + w2]) +
                 mp->d_density[POR_LIQ_PRES + w2] *
                     (mp->d_saturation[POR_LIQ_PRES + w1] * pmv->liq_Xvol_solvents[w] +
                      saturation * pmv->d_liq_Xvol_solvents[w][POR_LIQ_PRES + w1]) +
                 mp->d_density[POR_LIQ_PRES + w1] *
                     (mp->d_saturation[POR_LIQ_PRES + w2] * pmv->liq_Xvol_solvents[w] +
                      saturation * pmv->d_liq_Xvol_solvents[w][POR_LIQ_PRES + w2]));
            pmv->d_d_bulk_density[w][POR_LIQ_PRES + w1][POR_LIQ_PRES + w2] +=
                mp->d_porosity[POR_LIQ_PRES + w2] *
                (pmv->d_gas_density_solvents[w][POR_LIQ_PRES + w1] * (1. - saturation) -
                 pmv->gas_density_solvents[w] * mp->d_saturation[POR_LIQ_PRES + w1] +
                 liquid_density * mp->d_saturation[POR_LIQ_PRES + w1] * pmv->liq_Xvol_solvents[w] +
                 liquid_density * saturation * pmv->d_liq_Xvol_solvents[w][POR_LIQ_PRES + w1] +
                 saturation * pmv->liq_Xvol_solvents[w] * mp->d_density[w1]);
            pmv->d_d_bulk_density[w][POR_LIQ_PRES + w1][POR_LIQ_PRES + w2] +=
                mp->d_porosity[POR_LIQ_PRES + w1] *
                (pmv->d_gas_density_solvents[w][POR_LIQ_PRES + w2] * (1. - saturation) -
                 pmv->gas_density_solvents[w] * mp->d_saturation[POR_LIQ_PRES + w2] +
                 liquid_density * mp->d_saturation[POR_LIQ_PRES + w2] * pmv->liq_Xvol_solvents[w] +
                 liquid_density * saturation * pmv->d_liq_Xvol_solvents[w][POR_LIQ_PRES + w2] +
                 saturation * pmv->liq_Xvol_solvents[w] * mp->d_density[w2]);
          }
        }
      }
    }
    if (pd->e[pg->imtrx][R_POR_ENERGY]) {
      pmv->d_d_bulk_density[i_pe][i_pl][POR_TEMP] -= mp->d_porosity[POR_LIQ_PRES] * rhoCp;
      pmv->d_d_bulk_density[i_pe][i_pg][POR_TEMP] -= mp->d_porosity[POR_GAS_PRES] * rhoCp;
      pmv->d_d_bulk_density[i_pe][i_pe][POR_TEMP] -= mp->d_porosity[POR_TEMP] * rhoCp;
    }
  }

  return;
} /*  end  load_bulk_density */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void load_liq_perm(double porosity, double cap_pres, double saturation, double d_cap_pres[2])

/*************************************************************************
 * load_liq_perm --
 *
 *    This routine calculate the relative permeability of the liquid
 *    phase in a porous media.
 *    (HKM - actually it return k_rel / viscocity. However, this is
 *           slated to be simplified in the near future)
 *
 *  input:
 *  -----
 *
 *  Some porous properties are assumed to have been already calculated
 *  and stored within the mp structures. In particular the saturation
 *  and its derivatives are assumed to have been already calculated.
 *
 *  output:
 *  -----
 *
 *  Calculates the liquid permeability and its first derivatives
 *  with respect to all the problem unknowns. The first derivatives
 *  are calculated even if only a residual is needed for the case of
 *  supg weighting
 *
 *    mp->rel_liq_perm
 *    mp->d_rel_liq_perm[POR_LIQ_PRES]
 *    mp->d_rel_liq_perm[POR_GAS_PRES]
 *    mp->d_rel_liq_perm[POR_POROSITY]
 *    mp->d_rel_liq_perm[SHELL_PRESS_OPEN] OR mp->d_rel_liq_perm[SHELL_PRESS_OPEN_2]
 **************************************************************************/
{
  double s_eff, d_s_eff, rel_liq_perm, a1, factor, r_pore, rad, factor2;
  double d_rel_d_rad, expon2, sat_min, sat_max, viscosity, lambda;
  int i_rel_perm_ev;
  double scale;
  /*
   * Set these two derivatives to zero -> they are almost always zero
   */
  mp->d_rel_liq_perm[POR_GAS_PRES] = 0.0;
  mp->d_rel_liq_perm[POR_POROSITY] = 0.0;

  /*
   * Loop over the types of porosity models
   */
  if (mp->RelLiqPermModel == VAN_GENUCHTEN) {
    /*
     *
     * FOR VAN_GENUCHTEN EQUATION
     *  mp->u_rel_liq_perm[0] is the irreduceable water saturation
     *  mp->u_rel_liq_perm[1] is the irreduceable air saturation
     *  mp->u_rel_liq_perm[2] is the exponent, 1 - 1/beta
     *  mp->u_rel_liq_perm[3] is the liquid viscosity
     *
     *  Store some temporary variables
     */
    sat_min = mp->u_rel_liq_perm[0];
    sat_max = 1.0 - mp->u_rel_liq_perm[1];
    s_eff = (saturation - sat_min) / (sat_max - sat_min);
    viscosity = mp->u_rel_liq_perm[3];
    lambda = mp->u_rel_liq_perm[2];

    /*
     *  Clip the relative permeability to zero if the effective saturation
     *  is equal to or less than zero. -> there can be no transport
     *  in a liquid phase if there is no continguous pathway in that phase.
     */
    if (s_eff < 0.0) {
      mp->rel_liq_perm = 0.0;
      mp->d_rel_liq_perm[POR_LIQ_PRES] = 0.0;
      if (pd->e[pg->imtrx][R_SHELL_SAT_OPEN]) {
        mp->d_rel_liq_perm[SHELL_PRESS_OPEN] = 0.0;
      } else if (pd->e[pg->imtrx][R_SHELL_SAT_OPEN_2]) {
        mp->d_rel_liq_perm[SHELL_PRESS_OPEN_2] = 0.0;
      }

      if (pd->e[pg->imtrx][R_SHELL_SAT_1]) {
        mp->d_rel_liq_perm[SHELL_SAT_1] = 0.0;
      }
    }
    /*
     *  Clip the relative permeability at one -> it can never be
     *  greater than one.  Actually, note that if somehow s_eff is
     *  very close to 1.0 and fails this test, then you are dividing
     *  by zero as factor=1.0 below.    Now and then GOMA aborts due
     *  to this.
     */
    else if (s_eff >= 0.99999) {
      mp->rel_liq_perm = 1.0 / viscosity;
      mp->d_rel_liq_perm[POR_LIQ_PRES] = 0.0;
      if (pd->e[pg->imtrx][R_SHELL_SAT_OPEN]) {
        mp->d_rel_liq_perm[SHELL_PRESS_OPEN] = 0.0;
      } else if (pd->e[pg->imtrx][R_SHELL_SAT_OPEN_2]) {
        mp->d_rel_liq_perm[SHELL_PRESS_OPEN_2] = 0.0;
      }

      if (pd->e[pg->imtrx][R_SHELL_SAT_1]) {
        mp->d_rel_liq_perm[SHELL_SAT_1] = 0.0;
      }
    } else {
      expon2 = 1.0 / lambda;
      factor = pow(s_eff, expon2);
      factor2 = pow(1.0 - factor, lambda);
      a1 = 1.0 - factor2;
      rel_liq_perm = mp->rel_liq_perm = sqrt(s_eff) * a1 * a1 / viscosity;

      if (af->Assemble_Jacobian || mp->Porous_wt_funcModel == SUPG) {
        /*
         * Set these two derivatives to zero -> they are almost always zero
         */
        mp->d_rel_liq_perm[POR_GAS_PRES] = 0.0;
        mp->d_rel_liq_perm[POR_POROSITY] = 0.0;
        d_s_eff = 1.0 / (sat_max - sat_min);
        if (a1 == 0.0) {
          mp->d_rel_liq_perm[POR_LIQ_PRES] = 0.0;
          if (pd->e[pg->imtrx][R_SHELL_SAT_OPEN]) {
            mp->d_rel_liq_perm[SHELL_PRESS_OPEN] = 0.0;
          } else if (pd->e[pg->imtrx][R_SHELL_SAT_OPEN_2]) {
            mp->d_rel_liq_perm[SHELL_PRESS_OPEN_2] = 0.0;
          }

        } else {
          mp->d_rel_liq_perm[POR_LIQ_PRES] =
              d_s_eff * rel_liq_perm * mp->d_saturation[POR_LIQ_PRES] *
              (0.5 / s_eff + 2. * factor2 * factor / ((1.0 - factor) * s_eff * a1));
          if (pd->e[pg->imtrx][R_SHELL_SAT_OPEN]) {
            mp->d_rel_liq_perm[SHELL_PRESS_OPEN] =
                d_s_eff * rel_liq_perm * mp->d_saturation[SHELL_PRESS_OPEN] *
                (0.5 / s_eff + 2. * factor2 * factor / ((1.0 - factor) * s_eff * a1));
          } else if (pd->e[pg->imtrx][R_SHELL_SAT_OPEN_2]) {
            mp->d_rel_liq_perm[SHELL_PRESS_OPEN_2] =
                d_s_eff * rel_liq_perm * mp->d_saturation[SHELL_PRESS_OPEN_2] *
                (0.5 / s_eff + 2. * factor2 * factor / ((1.0 - factor) * s_eff * a1));
          }

          if (pd->e[pg->imtrx][R_SHELL_SAT_1]) {
            mp->d_rel_liq_perm[SHELL_SAT_1] =
                d_s_eff * rel_liq_perm *
                (0.5 / s_eff + 2. * factor2 * factor / ((1.0 - factor) * s_eff * a1));
          }
        }
        if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
          if (a1 == 0.0) {
            mp->d_rel_liq_perm[POR_LIQ_PRES] = 0.0;
          } else {
            mp->d_rel_liq_perm[POR_GAS_PRES] =
                d_s_eff * rel_liq_perm * mp->d_saturation[POR_GAS_PRES] *
                (0.5 / s_eff + 2. * factor2 * factor / ((1.0 - factor) * s_eff * a1));
          }
        }
      }
    }
  } else if (mp->RelLiqPermModel == PSD_VOL) {
    /*
     * FOR PSD_VOL EQUATION -
     *
     * all parameters are loaded from saturation and permeability
     * except liquid viscosity, which is parameter 0
     *
     *  mp->u_rel_liq_perm[0] is the liquid viscosity
     *
     * from saturation, there are three regimes:
     *  1) r_cap > r_pore  => Saturated - all pores full
     *  2) alpha r_pore < r_cap < r_pore  => Partially saturated - some pores full
     *  3) r_cap > alpha r_pore => No Saturation - no pores full
     */
    viscosity = mp->u_rel_liq_perm[0];
    if (pd->e[pg->imtrx][R_POR_POROSITY])
      r_pore = pmv->r_pore;
    else
      r_pore = mp->u_saturation[3];
    rad = pmv->r_cap / r_pore;

    if (pmv->r_cap > r_pore) {
      mp->rel_liq_perm = 1.0 / viscosity;
      mp->d_rel_liq_perm[POR_LIQ_PRES] = 0.0;
    } else if (pmv->r_cap < mp->u_permeability[2] * r_pore) {
      mp->rel_liq_perm = 0.0;
      mp->d_rel_liq_perm[POR_LIQ_PRES] = 0.0;
    } else {
      factor = 1.0 / (1. - CUBE(mp->u_permeability[2])) / viscosity;
      mp->rel_liq_perm = (CUBE(rad) - CUBE(mp->u_permeability[2])) * factor;

      if (af->Assemble_Jacobian) {
        mp->d_rel_liq_perm[POR_LIQ_PRES] =
            3. * pow(pmv->r_cap, 2.0) * pmv->d_r_cap[POR_LIQ_PRES] / pow(r_pore, 3.0) * factor;

        if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
          mp->d_rel_liq_perm[POR_GAS_PRES] =
              3.0 * pow(pmv->r_cap, 2.0) * pmv->d_r_cap[POR_GAS_PRES] / CUBE(r_pore) * factor;
        }
        if (pd->e[pg->imtrx][R_POR_POROSITY]) {
          mp->d_rel_liq_perm[POR_POROSITY] =
              -3. * CUBE(pmv->r_cap) / pow(r_pore, 4.0) * pmv->d_r_pore[POR_POROSITY] * factor;
        }
      }
    }
  } else if (mp->RelLiqPermModel == PSD_WEXP) {
    r_pore = mp->u_saturation[3];
    if (pd->e[pg->imtrx][R_POR_POROSITY])
      r_pore = pmv->r_pore;
    rad = pmv->r_cap / r_pore;

    mp->rel_liq_perm =
        ((1.0 - exp(-2.0 * rad) * (1 + 2. * rad + 2. * rad * rad + 4. / 3. * pow(rad, 3.0) +
                                   2. / 3. * pow(rad, 4.0) + 4. / 15. * pow(rad, 5.0))) /
         mp->u_rel_liq_perm[0]);

    d_rel_d_rad = ((2. * exp(-2. * rad) *
                        (1 + 2. * rad + 2. * rad * rad + 4. / 3. * pow(rad, 3.0) +
                         2. / 3. * pow(rad, 4.0) + 4. / 15. * pow(rad, 5.0)) -
                    exp(-2. * rad) * (2. + 4. * rad + 4. * rad * rad + 8. / 3. * pow(rad, 3.0) +
                                      4. / 3. * pow(rad, 4.0))) /
                   mp->u_rel_liq_perm[0]);

    if (af->Assemble_Jacobian) {
      mp->d_rel_liq_perm[POR_LIQ_PRES] = d_rel_d_rad * pmv->d_r_cap[POR_LIQ_PRES] / r_pore;

      if (pd->e[pg->imtrx][R_POR_GAS_PRES])
        mp->d_rel_liq_perm[POR_GAS_PRES] = d_rel_d_rad * pmv->d_r_cap[POR_GAS_PRES] / r_pore;
      if (pd->e[pg->imtrx][R_POR_POROSITY])
        mp->d_rel_liq_perm[POR_POROSITY] =
            -d_rel_d_rad * pmv->d_r_pore[POR_POROSITY] * rad / r_pore;
    }
  } else if (mp->RelLiqPermModel == PSD_SEXP) {
    r_pore = mp->u_saturation[3];
    if (pd->e[pg->imtrx][R_POR_POROSITY])
      r_pore = pmv->r_pore;
    factor = 9 / M_PIE / 16.;
    rad = pmv->r_cap / r_pore;
    mp->rel_liq_perm = (1. - exp(-factor * rad * rad) *
                                 (1 + factor * factor / 2. * pow(rad, 4.0) + factor * rad * rad)) /
                       mp->u_rel_liq_perm[0];
    d_rel_d_rad = (2 * factor * rad * exp(-factor * rad * rad) *
                       (1 + factor * factor / 2. * pow(rad, 4.0) + factor * rad * rad) -
                   exp(-factor * rad * rad) *
                       (4. * factor * factor / 2. * pow(rad, 3.0) + 2. * factor * rad)) /
                  mp->u_rel_liq_perm[0];

    if (af->Assemble_Jacobian) {
      mp->d_rel_liq_perm[POR_LIQ_PRES] = d_rel_d_rad * pmv->d_r_cap[POR_LIQ_PRES] / r_pore;

      if (pd->e[pg->imtrx][R_POR_GAS_PRES])
        mp->d_rel_liq_perm[POR_GAS_PRES] = d_rel_d_rad * pmv->d_r_cap[POR_GAS_PRES] / r_pore;
      if (pd->e[pg->imtrx][R_POR_POROSITY])
        mp->d_rel_liq_perm[POR_POROSITY] =
            -d_rel_d_rad * pmv->d_r_pore[POR_POROSITY] * rad / r_pore;
    }
  } else if (mp->RelLiqPermModel == EXTERNAL_FIELD) {
    /*
     *
     * FOR EXTERNAL FIELD
     *  mp->u_rel_liq_perm[0] is the scaling factor for the read-in external field variable
     */

    i_rel_perm_ev = mp->rel_liq_perm_external_field_index;
    scale = mp->u_rel_liq_perm[0];

    mp->rel_liq_perm = scale * fv->external_field[i_rel_perm_ev];
    mp->d_rel_liq_perm[POR_LIQ_PRES] = 0.0;
    if (pd->e[pg->imtrx][R_SHELL_SAT_OPEN]) {
      mp->d_rel_liq_perm[SHELL_PRESS_OPEN] = 0.0;
    } else if (pd->e[pg->imtrx][R_SHELL_SAT_OPEN_2]) {
      mp->d_rel_liq_perm[SHELL_PRESS_OPEN_2] = 0.0;
    }

    if (pd->e[pg->imtrx][R_SHELL_SAT_1]) {
      mp->d_rel_liq_perm[SHELL_SAT_1] = 0.0;
    }

  } else if (mp->RelLiqPermModel == VAN_GENUCHTEN_EXTERNAL) {
    /*
     *
     * FOR VAN_GENUCHTEN_EXTERNAL
     *  mp->u_rel_liq_perm[0] is the irreduceable water saturation
     *  mp->u_rel_liq_perm[1] is the irreduceable air saturation
     *  mp->u_rel_liq_perm[2] is the exponent, 1 - 1/beta for external field value of 0
     *  mp->u_rel_liq_perm[3] is the liquid viscosity
     *  mp->u_rel_liq_perm[4] is the exponent, 1 - 1/beta for external field value of 1
     */

    i_rel_perm_ev = mp->rel_liq_perm_external_field_index;

    sat_min = mp->u_rel_liq_perm[0];
    sat_max = 1.0 - mp->u_rel_liq_perm[1];
    s_eff = (saturation - sat_min) / (sat_max - sat_min);
    viscosity = mp->u_rel_liq_perm[3];

    /* Here I assume that efv is bounded between 0 and 1 */
    lambda = fv->external_field[i_rel_perm_ev] * (mp->u_rel_liq_perm[4] - mp->u_rel_liq_perm[2]) +
             mp->u_rel_liq_perm[2];

    /*
     *  Clip the relative permeability to zero if the effective saturation
     *  is equal to or less than zero. -> there can be no transport
     *  in a liquid phase if there is no continguous pathway in that phase.
     */
    if (s_eff < 0.0) {
      mp->rel_liq_perm = 0.0;
      mp->d_rel_liq_perm[POR_LIQ_PRES] = 0.0;
      if (pd->e[pg->imtrx][R_SHELL_SAT_OPEN]) {
        mp->d_rel_liq_perm[SHELL_PRESS_OPEN] = 0.0;
      } else if (pd->e[pg->imtrx][R_SHELL_SAT_OPEN_2]) {
        mp->d_rel_liq_perm[SHELL_PRESS_OPEN_2] = 0.0;
      }

      if (pd->e[pg->imtrx][R_SHELL_SAT_1]) {
        mp->d_rel_liq_perm[SHELL_SAT_1] = 0.0;
      }
    }
    /*
     *  Clip the relative permeability at one -> it can never be
     *  greater than one.  Actually, note that if somehow s_eff is
     *  very close to 1.0 and fails this test, then you are dividing
     *  by zero as factor=1.0 below.    Now and then GOMA aborts due
     *  to this.
     */
    else if (s_eff >= 0.99999) {
      mp->rel_liq_perm = 1.0 / viscosity;
      mp->d_rel_liq_perm[POR_LIQ_PRES] = 0.0;
      if (pd->e[pg->imtrx][R_SHELL_SAT_OPEN]) {
        mp->d_rel_liq_perm[SHELL_PRESS_OPEN] = 0.0;
      } else if (pd->e[pg->imtrx][R_SHELL_SAT_OPEN_2]) {
        mp->d_rel_liq_perm[SHELL_PRESS_OPEN_2] = 0.0;
      }
      if (pd->e[pg->imtrx][R_SHELL_SAT_1]) {
        mp->d_rel_liq_perm[SHELL_SAT_1] = 0.0;
      }
    } else {
      expon2 = 1.0 / lambda;
      factor = pow(s_eff, expon2);
      factor2 = pow(1.0 - factor, lambda);
      a1 = 1.0 - factor2;
      rel_liq_perm = mp->rel_liq_perm = sqrt(s_eff) * a1 * a1 / viscosity;

      if (af->Assemble_Jacobian || mp->Porous_wt_funcModel == SUPG) {
        /*
         * Set these two derivatives to zero -> they are almost always zero
         */
        mp->d_rel_liq_perm[POR_GAS_PRES] = 0.0;
        mp->d_rel_liq_perm[POR_POROSITY] = 0.0;
        d_s_eff = 1.0 / (sat_max - sat_min);
        if (a1 == 0.0) {
          mp->d_rel_liq_perm[POR_LIQ_PRES] = 0.0;
          if (pd->e[pg->imtrx][R_SHELL_SAT_OPEN]) {
            mp->d_rel_liq_perm[SHELL_PRESS_OPEN] = 0.0;
          } else if (pd->e[pg->imtrx][R_SHELL_SAT_OPEN_2]) {
            mp->d_rel_liq_perm[SHELL_PRESS_OPEN_2] = 0.0;
          }
        } else {
          mp->d_rel_liq_perm[POR_LIQ_PRES] =
              d_s_eff * rel_liq_perm * mp->d_saturation[POR_LIQ_PRES] *
              (0.5 / s_eff + 2. * factor2 * factor / ((1.0 - factor) * s_eff * a1));
          if (pd->e[pg->imtrx][R_SHELL_SAT_OPEN]) {
            mp->d_rel_liq_perm[SHELL_PRESS_OPEN] =
                d_s_eff * rel_liq_perm * mp->d_saturation[SHELL_PRESS_OPEN] *
                (0.5 / s_eff + 2. * factor2 * factor / ((1.0 - factor) * s_eff * a1));
          } else if (pd->e[pg->imtrx][R_SHELL_SAT_OPEN_2]) {
            mp->d_rel_liq_perm[SHELL_PRESS_OPEN_2] =
                d_s_eff * rel_liq_perm * mp->d_saturation[SHELL_PRESS_OPEN_2] *
                (0.5 / s_eff + 2. * factor2 * factor / ((1.0 - factor) * s_eff * a1));
          }

          if (pd->e[pg->imtrx][R_SHELL_SAT_1]) {
            mp->d_rel_liq_perm[SHELL_SAT_1] =
                d_s_eff * rel_liq_perm *
                (0.5 / s_eff + 2. * factor2 * factor / ((1.0 - factor) * s_eff * a1));
          }
        }
        if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
          if (a1 == 0.0) {
            mp->d_rel_liq_perm[POR_LIQ_PRES] = 0.0;
          } else {
            mp->d_rel_liq_perm[POR_GAS_PRES] =
                d_s_eff * rel_liq_perm * mp->d_saturation[POR_GAS_PRES] *
                (0.5 / s_eff + 2. * factor2 * factor / ((1.0 - factor) * s_eff * a1));
          }
        }
      }
    }
  } else {
    GOMA_EH(GOMA_ERROR, "No models for liquid relative permeability");
  }
  return;
} /* end  load_liq_perm  */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void load_gas_perm(double porosity, double cap_pres, double saturation, double d_cap_pres[2])

/************************************************************************
 * load_gas_perm --
 *
 *
 *  input:
 *  -----
 *
 *
 *  output:
 *  -----
 *   mp->rel_gas_perm
 *   mp->d_rel_gas_perm[]
 *************************************************************************/
{
  double liq_visc;

  if (mp->RelGasPermModel == SUM_TO_ONE) {
    /*
     * FOR SUM_TO_ONE EQUATION
     *  mp->u_rel_gas_perm[0] is the gas viscosity
     */
    if (mp->RelLiqPermModel == VAN_GENUCHTEN) {
      liq_visc = mp->u_rel_liq_perm[3];
    } else if (mp->RelLiqPermModel == PSD_VOL || mp->RelLiqPermModel == PSD_WEXP ||
               mp->RelLiqPermModel == PSD_SEXP) {
      liq_visc = mp->u_rel_liq_perm[0];
    } else {
      liq_visc = mp->viscosity;
      GOMA_WH(GOMA_ERROR, "Need liq visc in calc_rel_gas_perm");
    }

    mp->rel_gas_perm = (1. - mp->rel_liq_perm * liq_visc) / mp->u_rel_gas_perm[0];

    if (mp->rel_gas_perm < 0.0) {
      mp->rel_gas_perm = 0.0;
      mp->d_rel_gas_perm[POR_LIQ_PRES] = 0.0;
      mp->d_rel_gas_perm[POR_GAS_PRES] = 0.0;
      mp->d_rel_gas_perm[POR_POROSITY] = 0.0;
    } else {
      if (af->Assemble_Jacobian) {
        mp->d_rel_gas_perm[POR_LIQ_PRES] =
            -mp->d_rel_liq_perm[POR_LIQ_PRES] * liq_visc / mp->u_rel_gas_perm[0];

        if (pd->e[pg->imtrx][R_POR_GAS_PRES])
          mp->d_rel_gas_perm[POR_GAS_PRES] =
              -mp->d_rel_liq_perm[POR_GAS_PRES] * liq_visc / mp->u_rel_gas_perm[0];

        if (pd->e[pg->imtrx][R_POR_POROSITY])
          mp->d_rel_gas_perm[POR_POROSITY] =
              -mp->d_rel_liq_perm[POR_POROSITY] * liq_visc / mp->u_rel_gas_perm[0];
      }
    }
  } else {
    GOMA_EH(GOMA_ERROR, "No models for gas permeability");
  }
  return;
} /* end  load_gas_perm  */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void load_gas_diff(
    double porosity, double cap_pres, double saturation, double d_cap_pres[2], int spec)

/************************************************************************
 * load_gas_diff --
 *
 *       Calculate the gas phase diffusion coefficient for a
 *       single species, spec,  in a porous media from the unknowns
 *       used in the problem.
 *
 * The binary diffusion coefficient is represented by the following
 * expression:
 *
 *       D = D_o_binary * Porosity / tau * (1 - Sat) *
 *                        T_star**beta / P_star
 *
 *   where   T_star = T / T_ref
 *           P_star = P / P_ref
 *           tau = tortuosity
 *           D_o_binary = binary diffusion coefficient of species, spec,
 *                        with respect to the solvent gas species at
 *                        (T_star, P_star). For the solvent itself,
 *                        this is the self diffusion coefficient.
 *
 *  input:
 *  -----
 *
 *  FOR gas diffusion in POROUS medium:
 *  mp->u_porous_diffusivity[spec][0] is the binary diffusion
 *                                    coefficient in free space
 *  mp->u_porous_diffusivity[spec][1] is the tortuosity
 *  mp->u_porous_diffusivity[spec][2] is the reference gas pressure
 *  mp->u_porous_diffusivity[spec][3] is the reference temperature
 *  mp->u_porous_diffusivity[spec][4] is the exponent on the
 *                                    temperature dependence
 *
 *
 *  output:
 *  -----
 *       mp->porous_diffusivity[spec]
 *       mp->d_porous_diffusivity[spec][...]
 ************************************************************************/
{
  int w;
  double factor, T_star, P_star = 0.0, factorT = 0.0;

  if (mp->PorousDiffusivityModel[spec] == POROUS) {

    factor = (mp->u_porous_diffusivity[spec][0] / mp->u_porous_diffusivity[spec][1]);
    mp->porous_diffusivity[spec] = factor * mp->porosity * (1 - mp->saturation);

    if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
      P_star = fv->p_gas / mp->u_porous_diffusivity[spec][2];
      mp->porous_diffusivity[spec] /= P_star;
    }
    if (pd->v[pg->imtrx][TEMPERATURE]) {
      T_star = fv->T / mp->u_porous_diffusivity[spec][3];
      factorT = pow(T_star, mp->u_porous_diffusivity[spec][4]);
      mp->porous_diffusivity[spec] *= factorT;
    }

    if (af->Assemble_Jacobian) {
      /*
       * First add in the dependence on the porosity and saturation
       */
      for (w = 0; w < MAX_PMV; w++) {
        mp->d_porous_diffusivity[spec][POR_LIQ_PRES + w] =
            (factor * mp->d_porosity[POR_LIQ_PRES + w] * (1 - mp->saturation) -
             factor * mp->porosity * mp->d_saturation[POR_LIQ_PRES + w]);
        if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
          mp->d_porous_diffusivity[spec][POR_LIQ_PRES + w] /= P_star;
        }
        if (pd->v[pg->imtrx][TEMPERATURE]) {
          mp->d_porous_diffusivity[spec][POR_LIQ_PRES + w] *= factorT;
        }
      }
      if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
        mp->d_porous_diffusivity[spec][POR_GAS_PRES] -= mp->porous_diffusivity[spec] / fv->p_gas;
      }
      if (pd->v[pg->imtrx][TEMPERATURE]) {
        mp->d_porous_diffusivity[spec][TEMPERATURE] =
            mp->porous_diffusivity[spec] * mp->u_porous_diffusivity[spec][4] / fv->T;
      }
    }
  } else {
    GOMA_EH(GOMA_ERROR, "No models for gas phase diffusion in a porous medium");
  }
  return;
} /*  end  load_gas_diff  */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void load_MandE_flux(double porosity, double cap_pres, double saturation, double d_cap_pres[2])

/*************************************************************************
 *
 * load_MandE_flux --
 *
 *   Calculate mass flux of the solvent species, and energy flux, in a porous
 *   media from the unknowns used in the problem
 *   This routine only calculates the darcy fluxes and fickian fluxes
 *
 *  input:
 *  -----
 *  Constituitive properties and derivaties of indendent variables
 *  are assumed to be known before this call is invoked.
 *
 *  output:
 *  -----
 *
 *  pmv->liq_darcy_velocity[]  Darcy velocity of liquid phase
 *  pmv->gas_darcy_velocity[]  Darcy velocity of gas phase
 *  pmv->rel_mass_flux[i_pl][p]
 *                             Mass flux of liquid solvent due to phase
 *                             Darcy velocities and diffusion within
 *                             each phase
 *  pmv->rel_mass_flux[i_pg][p]
 *                             Mass flux of gas solvent due to phase
 *                             Darcy velocities and diffusion within
 *                             each phase
 *  pmv->d_rel_mass_flux[i_pl][p][j_pmv][j]
 *  pmv->d_rel_mass_flux[i_pg][p][j_pmv][j]
 *                             Dependence of the above quantity
 *                             on the j_pmv'th porous media variable's
 *                             j'th degree of freedom.
 ************************************************************************/
{
  int a, w, w1, w2, j, p, b, var, WIM;
  const int i_pl = 0, i_pg = 1, i_pore = 2, i_pe = 3;
  double t1, t2, t3, t4;
  double difflux[MAX_PMV][DIM], difflux_term;
  double grad_Ywg[DIM], d_grad_Ywg[DIM][MAX_PMV][MDE];

  if (pd->CoordinateSystem == CARTESIAN || pd->CoordinateSystem == CYLINDRICAL ||
      pd->CoordinateSystem == SWIRLING || pd->CoordinateSystem == PROJECTED_CARTESIAN ||
      pd->CoordinateSystem == CARTESIAN_2pt5D) {
    WIM = pd->Num_Dim;
  } else {
    WIM = VIM;
  }

  /*
   * DIFFUSIVE FLUX TERM - this is the term that is integrated by parts in the
   *   species conservation residual equation.  In a porous medium it contains
   *   five parts of the species flux:
   *    1) convection of the bulk concentration of each species with the solid
   *    2) convection of each species in gas phase relative to solid (from Darcy)
   *    3) convection of each species in liquid phase relative to solid (from Darcy)
   *    4) diffusion in gas phase (relative to gas convection)
   *    5) diffusion in liquid phase (relative to liquid convection) - currently this is neglected
   */

  for (a = 0; a < VIM; a++) {
    pmv->rel_mass_flux[i_pl][a] = 0.0;
    pmv->rel_mass_flux[i_pg][a] = 0.0;
    pmv->rel_mass_flux[i_pe][a] = 0.0;
  }

  if (cr->PorousFluxModel == DARCY_FICKIAN) { /*PRS (12/19/01): add other models here */

    /* Calculate the darcy velocity of the liquid and gas phases,
     * as well as their derivatives */

    calc_darcy_velocity();

    /* add convective contributions to total fluxes (Note: energy eqn is not
     * "mass" flux but uses same structure for consistency
     */

    for (a = 0; a < WIM; a++) {
      pmv->rel_mass_flux[i_pl][a] =
          pmv->liq_darcy_velocity[a] * pmv->liq_Xvol_solvents[i_pl] * mp->density +
          pmv->gas_darcy_velocity[a] * pmv->gas_density_solvents[i_pl];

      if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
        pmv->rel_mass_flux[i_pg][a] =
            pmv->liq_darcy_velocity[a] * pmv->liq_Xvol_solvents[i_pg] * mp->density +
            pmv->gas_darcy_velocity[a] * pmv->gas_density_solvents[i_pg];
      }
      if (pd->e[pg->imtrx][R_POR_ENERGY]) {
        pmv->rel_mass_flux[i_pe][a] =
            pmv->liq_darcy_velocity[a] * pmv->liq_Xvol_solvents[i_pl] * mp->density *
                pmv->enthalpy[0] +
            pmv->gas_darcy_velocity[a] * pmv->gas_density_solvents[i_pl] * pmv->enthalpy[1] +
            pmv->gas_darcy_velocity[a] * pmv->gas_density_solvents[i_pg] * pmv->enthalpy[2];
      }
    }
  } /* PorousFluxModel */

  /* add in diffusive flux in the gas phase only !!
   * the diffusive flux is proportional to the gradient in gas phase
   * concentration times the diffusion coefficient.
   * Here we get the gradient in gas phase concentration by the
   * chain rule - grad(cgas) = sum over variables (grad(var) * d(cgas)/d(var))
   * where the variables summed over are pressure of liq and gas, porosity, and temp
   *
   * For the energy equation, multiply the diffusive flux of water in the gas times
   * the enthalpy of water in the gas (pmv->enthalpy[1]), and the diffusive flux of
   * air in the gas times the enthalpy of air in the gas (pmv->enthalpy[2])
   */

  if (cr->PorousFluxModel == DARCY_FICKIAN) {

    /* In the binary system currently considered (water & air), the diffusive flux
     * of the air is equal and opposite to the diffusive flux of vapor (ie, water
     * in the gas phase - exploit this fact here and do not loop
     *
     *  for (w = 0; w < MAX_PMV; w++){
     *    if (w != i_pore && w != i_pe) {
     */
    w = 0;
    var = POR_LIQ_PRES;

    for (a = 0; a < WIM; a++) {
      /*
       * convert pressure, porosity and temperature gradients into gradients
       * of gas phase concentration
       */
      grad_Ywg[a] = pmv->d_Ywg[i_pl] * fv->grad_p_liq[a] + pmv->d_Ywg[i_pg] * fv->grad_p_gas[a] +
                    pmv->d_Ywg[i_pore] * fv->grad_porosity[a] + pmv->d_Ywg[i_pe] * fv->grad_T[a];

      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        for (w2 = 0; w2 < MAX_PMV; w2++) {
          d_grad_Ywg[a][w2][j] = pmv->d_Ywg[w2] * bf[var]->grad_phi[j][a] +
                                 bf[var]->phi[j] * (pmv->d_dYwg[i_pl][w2] * fv->grad_p_liq[a] +
                                                    pmv->d_dYwg[i_pg][w2] * fv->grad_p_gas[a] +
                                                    pmv->d_dYwg[i_pore][w2] * fv->grad_porosity[a] +
                                                    pmv->d_dYwg[i_pe][w2] * fv->grad_T[a]);
        }
      }

      difflux[w][a] = -mp->porous_diffusivity[w] * pmv->rhog * grad_Ywg[a];
      difflux[w + 1][a] = -difflux[w][a];
      pmv->rel_mass_flux[i_pl][a] += difflux[w][a];
      pmv->rel_mass_flux[i_pg][a] += difflux[w + 1][a];
      if (pd->e[pg->imtrx][R_POR_ENERGY])
        pmv->rel_mass_flux[i_pe][a] += difflux[i_pl][a] * (pmv->enthalpy[1] - pmv->enthalpy[2]);
    }

  } else {
    GOMA_EH(GOMA_ERROR, "Unimplemented flux constitutive relation in continuous media. Note: "
                        "POWERLAW_DARCY_FICKIAN not available for POROUS_TWO_PHASE");
  }

  /* Add conduction in the solid */

  if (pd->e[pg->imtrx][R_POR_ENERGY]) {
    for (a = 0; a < WIM; a++)
      pmv->rel_mass_flux[i_pe][a] += -mp->thermal_conductivity * fv->grad_T[a];
  }

  if (af->Assemble_Jacobian) {
    var = POR_LIQ_PRES;

    for (w = 0; w < MAX_PMV; w++) {
      for (a = 0; a < WIM; a++)

      {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          for (w1 = 0; w1 < MAX_PMV; w1++) {
            pmv->d_rel_mass_flux_dpmv[w][a][w1][j] = 0.;
          }
        }
      }
    }

    if (cr->PorousFluxModel == DARCY_FICKIAN) {
      for (w = 0; w < MAX_PMV; w++) {

        if (w != i_pore && w != i_pe) { /* no term for solid phase */

          for (a = 0; a < WIM; a++) {
            for (w1 = 0; w1 < MAX_PMV; w1++) {
              var = POR_LIQ_PRES + w1;
              if (pd->v[pg->imtrx][var]) {
                for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {

                  t1 = bf[var]->phi[j] * pmv->gas_darcy_velocity[a] *
                       pmv->d_gas_density_solvents[w][POR_LIQ_PRES + w1];

                  t2 = bf[var]->phi[j] * pmv->liq_darcy_velocity[a] * mp->density *
                       pmv->d_liq_Xvol_solvents[w][POR_LIQ_PRES + w1];

                  t3 =
                      pmv->d_liq_darcy_velocity[a][w1][j] * pmv->liq_Xvol_solvents[w] * mp->density;

                  t4 = pmv->d_gas_darcy_velocity[a][w1][j] * pmv->gas_density_solvents[w];

                  pmv->d_rel_mass_flux_dpmv[w][a][w1][j] += t1 + t2 + t3 + t4;

                  if (pd->e[pg->imtrx][R_POR_ENERGY]) {
                    pmv->d_rel_mass_flux_dpmv[i_pe][a][w1][j] +=
                        (t2 + t3) * pmv->enthalpy[0] + (t1 + t4) * pmv->enthalpy[w + 1] +
                        bf[var]->phi[j] *
                            (pmv->liq_darcy_velocity[a] * pmv->liq_Xvol_solvents[w] * mp->density *
                                 pmv->d_enthalpy[0][POR_LIQ_PRES + w1] +
                             pmv->gas_darcy_velocity[a] * pmv->gas_density_solvents[w] *
                                 pmv->d_enthalpy[w + 1][POR_LIQ_PRES + w1]);
                  }
                } /* for (j=0 ... */
              }
            } /* for (w1=0 ... */
          }   /* for (a=0 ... */
        }     /* if(w != ... */
      }       /* for (w=0 ... */
    }         /* PorousFluxModel */

    var = POR_SINK_MASS;
    if (pd->v[pg->imtrx][var]) {
      for (a = 0; a < WIM; a++) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          for (w1 = 0; w1 < MAX_PMV; w1++) {
            pmv->d_rel_mass_flux_dSM[w1][a][j] = 0.;
          }
        }
      }

      if (cr->PorousFluxModel == DARCY_FICKIAN) {
        for (w = 0; w < MAX_PMV; w++) {
          for (a = 0; a < WIM; a++) {
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              /*for now we know it only liquid phase depends on sink_mass */

              t3 = pmv->d_liq_darcy_velocity_dSM[a][j] * pmv->liq_Xvol_solvents[w] * mp->density;

              pmv->d_rel_mass_flux_dSM[i_pl][a][j] += t3;

            } /* for (j=0 ... */
          }   /* for (a=0 ... */
        }     /* for (w = 0 ... */
      }       /* PorousFluxModel */
    }

    if (cr->PorousFluxModel == DARCY_FICKIAN) {
      /* In the binary system currently considered (water & air), the diffusive flux
       * of the air is equal and opposite to the diffusive flux of vapor (ie, water
       * in the gas phase - exploit this fact here and do not loop
       *
       *  for (w = 0; w < MAX_PMV; w++){
       *    if (w != i_pore && w != i_pe) {
       */

      w = 0; /* Water equation */

      /* Derivative with respect to Pliq */

      for (a = 0; a < WIM; a++) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          difflux_term =
              -bf[var]->phi[j] * ((mp->porous_diffusivity[w] * pmv->d_rhog[i_pl] +
                                   mp->d_porous_diffusivity[w][POR_LIQ_PRES] * pmv->rhog) *
                                  grad_Ywg[a]) -
              mp->porous_diffusivity[w] * pmv->rhog * d_grad_Ywg[a][i_pl][j];
          pmv->d_rel_mass_flux_dpmv[w][a][i_pl][j] += difflux_term;
          pmv->d_rel_mass_flux_dpmv[w + 1][a][i_pl][j] += -difflux_term;
          if (pd->e[pg->imtrx][R_POR_ENERGY])
            pmv->d_rel_mass_flux_dpmv[i_pe][a][i_pl][j] +=
                difflux_term * (pmv->enthalpy[1] - pmv->enthalpy[2]) +
                difflux[i_pl][a] * bf[var]->phi[j] *
                    (pmv->d_enthalpy[1][POR_LIQ_PRES] - pmv->d_enthalpy[2][POR_LIQ_PRES]);

          /* Derivative with respect to Pgas */
          if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
            difflux_term =
                -bf[var]->phi[j] * ((mp->porous_diffusivity[w] * pmv->d_rhog[i_pg] +
                                     mp->d_porous_diffusivity[w][POR_GAS_PRES] * pmv->rhog) *
                                    grad_Ywg[a]) -
                mp->porous_diffusivity[w] * pmv->rhog * d_grad_Ywg[a][i_pg][j];
            pmv->d_rel_mass_flux_dpmv[w][a][i_pg][j] += difflux_term;
            pmv->d_rel_mass_flux_dpmv[w + 1][a][i_pg][j] += -difflux_term;
            if (pd->e[pg->imtrx][R_POR_ENERGY])
              pmv->d_rel_mass_flux_dpmv[i_pe][a][i_pg][j] +=
                  difflux_term * (pmv->enthalpy[1] - pmv->enthalpy[2]) +
                  difflux[i_pl][a] * bf[var]->phi[j] *
                      (pmv->d_enthalpy[1][POR_GAS_PRES] - pmv->d_enthalpy[2][POR_GAS_PRES]);
          }

          if (pd->e[pg->imtrx][R_POR_POROSITY]) {
            difflux_term = pmv->d_Ywg[i_pore] * pmv->rhog * mp->porous_diffusivity[w] *
                           bf[var]->grad_phi[j][a];
            pmv->d_rel_mass_flux_dpmv[w][a][i_pl][j] -= difflux_term;
            pmv->d_rel_mass_flux_dpmv[w + 1][a][i_pl][j] += difflux_term;
            if (pd->e[pg->imtrx][R_POR_ENERGY])
              pmv->d_rel_mass_flux_dpmv[i_pe][a][i_pl][j] -=
                  difflux_term * (pmv->enthalpy[1] - pmv->enthalpy[2]);
            for (w2 = 0; w2 < MAX_PMV; w2++) {
              difflux_term =
                  bf[var]->phi[j] * (mp->porous_diffusivity[w] * pmv->rhog * fv->grad_porosity[a] *
                                         pmv->d_dYwg[i_pore][w2] +
                                     (mp->d_porous_diffusivity[w][POR_LIQ_PRES + w2] * pmv->rhog +
                                      mp->porous_diffusivity[w] * pmv->d_rhog[w2]) *
                                         pmv->d_Ywg[i_pore] * fv->grad_porosity[a]);
              pmv->d_rel_mass_flux_dpmv[w][a][w2][j] -= difflux_term;
              pmv->d_rel_mass_flux_dpmv[w + 1][a][w2][j] += difflux_term;
              if (pd->e[pg->imtrx][R_POR_ENERGY])
                pmv->d_rel_mass_flux_dpmv[i_pe][a][w2][j] -=
                    difflux_term * (pmv->enthalpy[1] - pmv->enthalpy[2]) +
                    difflux[i_pore][a] * bf[var]->phi[j] *
                        (pmv->d_enthalpy[1][POR_LIQ_PRES + w2] -
                         pmv->d_enthalpy[2][POR_LIQ_PRES + w2]);
            }
          }
          /* Derivative with respect to T */

          if (pd->e[pg->imtrx][R_POR_ENERGY]) {
            difflux_term = -bf[var]->phi[j] * ((mp->porous_diffusivity[w] * pmv->d_rhog[i_pe] +
                                                mp->d_porous_diffusivity[w][POR_TEMP] * pmv->rhog) *
                                               grad_Ywg[a]) -
                           mp->porous_diffusivity[w] * pmv->rhog * d_grad_Ywg[a][i_pe][j];
            pmv->d_rel_mass_flux_dpmv[w][a][i_pe][j] += difflux_term;
            pmv->d_rel_mass_flux_dpmv[w + 1][a][i_pe][j] += -difflux_term;

            pmv->d_rel_mass_flux_dpmv[i_pe][a][i_pe][j] +=
                difflux_term * (pmv->enthalpy[1] - pmv->enthalpy[2]) +
                difflux[i_pl][a] * bf[var]->phi[j] *
                    (pmv->d_enthalpy[1][POR_TEMP] - pmv->d_enthalpy[2][POR_TEMP]);
          }

        } /* j<ei[pg->imtrx]->dof[var] */
      }   /* a<VIM */
    }     /* DARCY_FICKIAN */

    /* Additional terms for energy equation - conduction in solid
     * Assume thermal conductivity a function of temperature only
     */
    if (pd->e[pg->imtrx][R_POR_ENERGY]) {
      var = POR_TEMP;
      for (a = 0; a < WIM; a++) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          pmv->d_rel_mass_flux_dpmv[i_pe][a][i_pe][j] +=
              -mp->thermal_conductivity * bf[var]->grad_phi[j][a] -
              mp->d_thermal_conductivity[var] * bf[var]->phi[j] * fv->grad_T[a];
        }
      }
    }

    for (b = 0; b < pd->Num_Dim; b++) {
      var = MESH_DISPLACEMENT1 + b;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          for (a = 0; a < WIM; a++) {
            for (w = 0; w < MAX_PMV; w++) {
              pmv->d_rel_mass_flux_dmesh[w][a][b][j] = 0.;
            }
          }
          for (a = 0; a < WIM; a++) {
            if (cr->PorousFluxModel == DARCY_FICKIAN) {
              for (w = 0; w < MAX_PMV; w++) {
                if (w != i_pore) { /* no diffusion term for solid phase */
                  if (mp->PermeabilityModel != K_TENSOR && mp->PermeabilityModel != ORTHOTROPIC &&
                      mp->PermeabilityModel != KC_TENSOR && mp->PermeabilityModel != SM_TENSOR) {
                    pmv->d_rel_mass_flux_dmesh[w][a][b][j] +=
                        -mp->permeability * mp->rel_liq_perm * pmv->liq_Xvol_solvents[w] *
                        fv->d_grad_p_liq_dmesh[a][b][j] * mp->density;
                    if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
                      pmv->d_rel_mass_flux_dmesh[w][a][b][j] +=
                          -mp->permeability * mp->rel_gas_perm * pmv->gas_density_solvents[w] *
                          fv->d_grad_p_gas_dmesh[a][b][j];
                    }
                  } else {
                    for (p = 0; p < WIM; p++) {
                      pmv->d_rel_mass_flux_dmesh[w][a][b][j] +=
                          -mp->perm_tensor[a][p] * mp->rel_liq_perm * pmv->liq_Xvol_solvents[w] *
                          fv->d_grad_p_liq_dmesh[p][b][j] * mp->density;
                      if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
                        pmv->d_rel_mass_flux_dmesh[w][a][b][j] +=
                            -mp->perm_tensor[a][p] * mp->rel_gas_perm *
                            pmv->gas_density_solvents[w] * fv->d_grad_p_gas_dmesh[p][b][j];
                      }
                      if (mp->PermeabilityModel == ORTHOTROPIC ||
                          mp->PermeabilityModel == KC_TENSOR ||
                          mp->PermeabilityModel == SM_TENSOR) {
                        pmv->d_rel_mass_flux_dmesh[w][a][b][j] +=
                            -mp->d_perm_tensor_dx[a][p][b][j] * mp->rel_liq_perm *
                            pmv->liq_Xvol_solvents[w] * fv->grad_p_liq[p] * mp->density;
                        if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
                          pmv->d_rel_mass_flux_dmesh[w][a][b][j] +=
                              -mp->d_perm_tensor_dx[a][p][b][j] * mp->rel_gas_perm *
                              pmv->gas_density_solvents[w] * fv->grad_p_gas[p];
                        }
                      }
                    }
                  }
                } /* w!=i_pore */
              }   /* w<MAX_PMV */

              /* Additional terms for energy equation */
              if (pd->e[pg->imtrx][R_POR_GAS_PRES])
                pmv->d_rel_mass_flux_dmesh[i_pe][a][b][j] +=
                    -mp->thermal_conductivity * fv->d_grad_T_dmesh[a][b][j];

            } /* DARCY_FICKIAN */
            if (cr->PorousFluxModel == DARCY_FICKIAN) {
              for (w = 0; w < MAX_PMV; w++) {
                if (w != i_pore) { /* no diffusion term for solid phase */
                  pmv->d_rel_mass_flux_dmesh[w][a][b][j] -=
                      mp->porous_diffusivity[w] * (pmv->d_gas_density_solvents[w][POR_LIQ_PRES] *
                                                       fv->d_grad_p_liq_dmesh[a][b][j] +
                                                   pmv->d_gas_density_solvents[w][POR_GAS_PRES] *
                                                       fv->d_grad_p_gas_dmesh[a][b][j] +
                                                   pmv->d_gas_density_solvents[w][POR_POROSITY] *
                                                       fv->d_grad_porosity_dmesh[a][b][j]);
                }
              }
            } /* DARCY_FICKIAN */
          }   /* a<WIM */
        }     /* j<ei[pg->imtrx]->dof[var] */
      }       /*if ( pd->v[pg->imtrx][var] ) */
    }         /*  b<pd->Num_Dim */

    var = TEMPERATURE;
    if (pd->v[pg->imtrx][var] && !pd->e[pg->imtrx][R_POR_ENERGY]) {
      if (cr->PorousFluxModel == DARCY_FICKIAN) {
        /* NOTE - the gas concentrations, etc should also be functions of temperature */
        for (w = 0; w < MAX_PMV; w++) {
          if (w != i_pore) { /* no diffusion term for solid phase */
            for (a = 0; a < WIM; a++) {
              for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                pmv->d_rel_mass_flux_dT[w][a][j] = 0.;
                pmv->d_rel_mass_flux_dT[w][a][j] -=
                    mp->d_porous_diffusivity[w][var] * bf[var]->phi[j] *
                    (pmv->d_gas_density_solvents[w][POR_LIQ_PRES] * fv->grad_p_liq[a] +
                     pmv->d_gas_density_solvents[w][POR_GAS_PRES] * fv->grad_p_gas[a] +
                     pmv->d_gas_density_solvents[w][POR_POROSITY] * fv->grad_porosity[a]);
              }
            }
          }
        }
      }
    }
  }
  return;
} /* end  load_MandE_flux  */
/*******************************************************************************/

void calc_grad_Ywg(int var, double *grad_Ywg) {
  int dofs, p, i;
  int i_pl = 0;
  double Pg, T, Ywg;
  double Pvsat, dPvsat_dT, rho_gas, drho_gas_dPg, drho_gas_dT;
  double d_drho_gas_dPg[2], d_drho_gas_dT[2];
  double rho_wg_sat, drho_wgsat_dT, d_drho_wgsat_dT;

  dofs = ei[pg->imtrx]->dof[var];
  for (p = 0; p < VIM; p++) {
    grad_Ywg[p] = 0.;
    for (i = 0; i < dofs; i++) {

      if (pd->e[pg->imtrx][R_POR_GAS_PRES])
        Pg = *esp->p_gas[i];
      else
        Pg = mp->u_porous_gas_constants[3];
      /* perturb temperature */
      if (pd->e[pg->imtrx][R_POR_ENERGY])
        T = *esp->T[i] * 1.01;
      else
        T = mp->u_porous_vapor_pressure[i_pl][4] - 273.15;

      rho_gas = calc_rho_gas(Pg, T, &Pvsat, &dPvsat_dT, &drho_gas_dPg, &drho_gas_dT, d_drho_gas_dPg,
                             d_drho_gas_dT, &rho_wg_sat, &drho_wgsat_dT, &d_drho_wgsat_dT);

      Ywg = rho_wg_sat / rho_gas;
      /* test case */
      Ywg = 1.0e-3 * T;

      grad_Ywg[p] += Ywg * bf[var]->grad_phi[i][p];
    }
  }

} /* end calc_grad_Ywg */
/******************************************************************************/
void calc_darcy_velocity() {
  int a, b, j, w1, var, WIM;
  const int i_pl = 0, i_pg = 1, i_pore = 2;

  /*initialize some important things */
  memset(pmv->liq_darcy_velocity, 0, sizeof(double) * DIM);
  memset(pmv->gas_darcy_velocity, 0, sizeof(double) * DIM);
  memset(pmv->d_liq_darcy_velocity, 0, sizeof(double) * MDE * DIM * MAX_PMV);
  memset(pmv->d_gas_darcy_velocity, 0, sizeof(double) * MDE * DIM * MAX_PMV);

  if (pd->CoordinateSystem == CARTESIAN || pd->CoordinateSystem == CYLINDRICAL ||
      pd->CoordinateSystem == SWIRLING || pd->CoordinateSystem == PROJECTED_CARTESIAN ||
      pd->CoordinateSystem == CARTESIAN_2pt5D) {
    WIM = pd->Num_Dim;
  } else {
    WIM = VIM;
  }

  /*
   * Calculate convection velocity in liquid and gas phases
   */
  if (mp->PermeabilityModel != K_TENSOR && mp->PermeabilityModel != ORTHOTROPIC &&
      mp->PermeabilityModel != SM_TENSOR && mp->PermeabilityModel != KC_TENSOR) {

    for (a = 0; a < WIM; a++) {
      pmv->liq_darcy_velocity[a] =
          -mp->permeability * mp->rel_liq_perm * (fv->grad_p_liq[a] - mp->momentum_source[a]);

      if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
        /* neglect gravitational effect in gas, for now */
        pmv->gas_darcy_velocity[a] = -mp->permeability * mp->rel_gas_perm * fv->grad_p_gas[a];
      } else {
        pmv->gas_darcy_velocity[a] = 0.0;
      }
    }

  } else {

    for (a = 0; a < WIM; a++) {
      for (b = 0; b < WIM; b++) {
        pmv->liq_darcy_velocity[a] -=
            mp->perm_tensor[a][b] * mp->rel_liq_perm * (fv->grad_p_liq[b] - mp->momentum_source[b]);

        if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
          /* neglect gravitational effect in gas, for now */
          pmv->gas_darcy_velocity[a] -=
              mp->perm_tensor[a][b] * mp->rel_gas_perm * fv->grad_p_gas[b];
        } else {
          pmv->gas_darcy_velocity[a] = 0.0;
        }
      } /* b < WIM */
    }   /* a < WIM */
  }     /* PermeabilityModel */

  if (af->Assemble_Jacobian) {

    if (mp->PermeabilityModel != K_TENSOR && mp->PermeabilityModel != ORTHOTROPIC &&
        mp->PermeabilityModel != SM_TENSOR && mp->PermeabilityModel != KC_TENSOR) {
      for (a = 0; a < WIM; a++) {
        var = POR_LIQ_PRES;
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          for (w1 = 0; w1 < MAX_PMV; w1++) {

            pmv->d_liq_darcy_velocity[a][w1][j] -=
                bf[var]->phi[j] * (fv->grad_p_liq[a] - mp->momentum_source[a]) *
                (mp->permeability * mp->d_rel_liq_perm[POR_LIQ_PRES + w1] +
                 mp->d_permeability[POR_LIQ_PRES + w1] * mp->rel_liq_perm);

            if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
              pmv->d_gas_darcy_velocity[a][w1][j] -=
                  bf[var]->phi[j] * fv->grad_p_gas[a] *
                  (mp->permeability * mp->d_rel_gas_perm[POR_LIQ_PRES + w1] +
                   mp->d_permeability[POR_LIQ_PRES + w1] * mp->rel_gas_perm);
            }
          }

          /* sensitivity to liquid pressure gradient */
          pmv->d_liq_darcy_velocity[a][i_pl][j] -=
              mp->permeability * mp->rel_liq_perm * bf[var]->grad_phi[j][a];

          /* sensitivity to gas pressure gradient and permeability */
          if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
            pmv->d_gas_darcy_velocity[a][i_pg][j] -=
                mp->permeability * mp->rel_gas_perm * bf[var]->grad_phi[j][a];
          }
        }

        /* sensitivity to POR_SINK_MASS */
        var = POR_SINK_MASS;
        if (pd->e[pg->imtrx][R_POR_SINK_MASS]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {

            pmv->d_liq_darcy_velocity_dSM[a][j] -=
                bf[var]->phi[j] * (fv->grad_p_liq[a] - mp->momentum_source[a]) *
                (mp->d_permeability[POR_SINK_MASS] * mp->rel_liq_perm);
          }
        }
      }
    } else { /* Permeability is a tensor, and it is constant */

      var = POR_LIQ_PRES;
      for (a = 0; a < WIM; a++) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          for (w1 = 0; w1 < MAX_PMV; w1++) {

            for (b = 0; b < WIM; b++) {

              pmv->d_liq_darcy_velocity[a][w1][j] -=
                  bf[var]->phi[j] * (mp->perm_tensor[a][b] * mp->d_rel_liq_perm[POR_LIQ_PRES + w1] *
                                     (fv->grad_p_liq[b] - mp->momentum_source[b]));

              if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
                pmv->d_gas_darcy_velocity[a][w1][j] -= bf[var]->phi[j] * mp->perm_tensor[a][b] *
                                                       mp->d_rel_gas_perm[POR_LIQ_PRES + w1] *
                                                       fv->grad_p_gas[b];
              }
            }
          }

          /* sensitivity to liquid pressure gradient */

          for (b = 0; b < WIM; b++) {

            pmv->d_liq_darcy_velocity[a][i_pl][j] -=
                mp->perm_tensor[a][b] * mp->rel_liq_perm * bf[var]->grad_phi[j][b];

            /* sensitivity to gas pressure gradient and permeability */
            if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
              pmv->d_gas_darcy_velocity[a][i_pg][j] -=
                  mp->perm_tensor[a][b] * mp->rel_gas_perm * bf[var]->grad_phi[j][b];
            }
          } /* for (b=0 ...) */
        }   /* for (J=0 ...) */
      }     /* for (a=0 ...) */

      var = POR_POROSITY;
      if (pd->v[pg->imtrx][var] && mp->PermeabilityModel == KC_TENSOR) {
        for (a = 0; a < WIM; a++) {
          for (b = 0; b < WIM; b++) {
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              pmv->d_liq_darcy_velocity[a][i_pore][j] -= bf[var]->phi[j] *
                                                         mp->d_perm_tensor[a][b][POR_POROSITY] *
                                                         mp->rel_liq_perm * fv->grad_p_liq[b];
            }
          }
        }
      }

      /* sensitivity to POR_SINK_MASS */
      var = POR_SINK_MASS;
      if (pd->e[pg->imtrx][R_POR_SINK_MASS]) {
        for (a = 0; a < WIM; a++) {
          for (b = 0; b < WIM; b++) {
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {

              pmv->d_liq_darcy_velocity_dSM[a][j] -=
                  bf[var]->phi[j] * (fv->grad_p_liq[a] - mp->momentum_source[a]) *
                  (mp->d_perm_tensor[a][b][POR_SINK_MASS] * mp->rel_liq_perm);
            }
          }
        }
      }

    } /*end if (PERM TENSOR) */
  }   /* Assemble_Jacobian */
  return;
}

/*****************************************************************************/

/*****************************************************************************/

void load_mass_flux(double porosity, double cap_pres, double saturation, double d_cap_pres[2])

/*************************************************************************
 *
 * load_mass_flux --
 *
 *   Calculate mass flux of the solvent species in a porous
 *   media from the unknowns used in the problem
 *   This routine only calculates the darcy fluxes and fickian fluxes
 *
 *  input:
 *  -----
 *  Constituitive properties and derivaties of indendent variables
 *  are assumed to be known before this call is invoked.
 *
 *  output:
 *  -----
 *
 *  pmv->liq_darcy_velocity[]  Darcy velocity of liquid phase
 *  pmv->gas_darcy_velocity[]  Darcy velocity of gas phase
 *  pmv->rel_mass_flux[i_pl][p]
 *                             Mass flux of liquid solvent due to phase
 *                             Darcy velocities and diffusion within
 *                             each phase
 *  pmv->rel_mass_flux[i_pg][p]
 *                             Mass flux of gas solvent due to phase
 *                             Darcy velocities and diffusion within
 *                             each phase
 *  pmv->d_rel_mass_flux[i_pl][p][j_pmv][j]
 *  pmv->d_rel_mass_flux[i_pg][p][j_pmv][j]
 *                             Dependence of the above quantity
 *                             on the j_pmv'th porous media variable's
 *                             j'th degree of freedom.
 ************************************************************************/
{
  int a, w, w1, w2, j, p, b, var;
  const int i_pl = 0, i_pg = 1, i_pore = 2, i_pe = 3;
  double difflux[MAX_PMV][DIM], difflux_term, n_pow;

  /*
   * DIFFUSIVE FLUX TERM - this is the term that is integrated by parts in the
   *   species conservation residual equation.  In a porous medium it contains
   *   five parts of the species flux:
   *    1) convection of the bulk concentration of each species with the solid
   *    2) convection of each species in gas phase relative to solid (from Darcy)
   *    3) convection of each species in liquid phase relative to solid (from Darcy)
   *    4) diffusion in gas phase (relative to gas convection)
   *    5) diffusion in liquid phase (relative to liquid convection) - currently this is neglected
   */

  memset(pmv->liq_darcy_velocity, 0, sizeof(double) * DIM);
  memset(pmv->gas_darcy_velocity, 0, sizeof(double) * DIM);
  memset(pmv->d_liq_darcy_velocity, 0, sizeof(double) * MDE * DIM * MAX_PMV);
  memset(pmv->d_gas_darcy_velocity, 0, sizeof(double) * MDE * DIM * MAX_PMV);

  if (cr->PorousFluxModel == DARCY_FICKIAN) { /*PRS (12/19/01): add other models here */
    /*
     * Calculate convection velocity in liquid and gas phases
     */
    if (mp->PermeabilityModel != K_TENSOR && mp->PermeabilityModel != ORTHOTROPIC &&
        mp->PermeabilityModel != SM_TENSOR && mp->PermeabilityModel != KC_TENSOR) {
      for (a = 0; a < VIM; a++) {
        pmv->liq_darcy_velocity[a] =
            -mp->permeability * mp->rel_liq_perm * (fv->grad_p_liq[a] - mp->momentum_source[a]);
        if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
          pmv->gas_darcy_velocity[a] = -mp->permeability * mp->rel_gas_perm * fv->grad_p_gas[a];
        } else {
          pmv->gas_darcy_velocity[a] = 0.0;
        }
        /* neglect gravitational effect in gas, for now */
      }
    } else {
      for (a = 0; a < VIM; a++) {
        pmv->liq_darcy_velocity[a] = 0.0;
        pmv->gas_darcy_velocity[a] = 0.0;
        for (b = 0; b < VIM; b++) {
          pmv->liq_darcy_velocity[a] -= mp->perm_tensor[a][b] * mp->rel_liq_perm *
                                        (fv->grad_p_liq[b] - mp->momentum_source[b]);
          if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
            pmv->gas_darcy_velocity[a] -=
                mp->perm_tensor[a][b] * mp->rel_gas_perm * fv->grad_p_gas[b];
          }
          /* neglect gravitational effect in gas, for now */
        }
      }
    }
    /* add convection fluxes into the total diffusion flux term
     * here vconv is the velocity of the solid phase */
    for (a = 0; a < VIM; a++) {
      pmv->rel_mass_flux[i_pl][a] = 0.0;
      pmv->rel_mass_flux[i_pg][a] = 0.0;
      pmv->rel_mass_flux[i_pe][a] = 0.0;
    }

    for (a = 0; a < VIM; a++) {
      pmv->rel_mass_flux[i_pl][a] =
          pmv->gas_darcy_velocity[a] * pmv->gas_density_solvents[i_pl] +
          pmv->liq_darcy_velocity[a] * pmv->liq_Xvol_solvents[i_pl] * mp->density;
      if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
        pmv->rel_mass_flux[i_pg][a] =
            pmv->gas_darcy_velocity[a] * pmv->gas_density_solvents[i_pg] +
            pmv->liq_darcy_velocity[a] * pmv->liq_Xvol_solvents[i_pg] * mp->density;
      }
      /* Recall that in order to maintain a common structure with the above eqns,
       * the heat equation uses pmv->gas_density_solvents[i_pe] to hold the product
       * of gas density and gas enthalpy, and pmv->liq_Xvol_solvents[i_pe] holds
       * the enthalpy of the liquid.
       */
      if (pd->e[pg->imtrx][R_POR_ENERGY]) {
        pmv->rel_mass_flux[i_pe][a] =
            -mp->thermal_conductivity * fv->grad_T[a] +
            pmv->liq_darcy_velocity[a] * pmv->liq_Xvol_solvents[i_pe] * mp->density *
                pmv->enthalpy[0] +
            pmv->gas_darcy_velocity[a] * (pmv->gas_density_solvents[i_pl] * pmv->enthalpy[1] +
                                          pmv->gas_density_solvents[i_pg] * pmv->enthalpy[2]);
      }
    }
  } else if (cr->PorousFluxModel == POWERLAW_DARCY_FICKIAN) {
    /*Darcy powerlaw model goes here*/
    /*
     * Calculate convection velocity in liquid and gas phases
     */

    n_pow = 1. / 1.;

    if (mp->PermeabilityModel != K_TENSOR && mp->PermeabilityModel != ORTHOTROPIC &&
        mp->PermeabilityModel != SM_TENSOR && mp->PermeabilityModel != KC_TENSOR) {
      for (a = 0; a < VIM; a++) {
        pmv->liq_darcy_velocity[a] = -mp->permeability * mp->rel_liq_perm *
                                     (pow(fv->grad_p_liq[a], n_pow) - mp->momentum_source[a]);
        if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
          pmv->gas_darcy_velocity[a] = -mp->permeability * mp->rel_gas_perm * fv->grad_p_gas[a];
        } else {
          pmv->gas_darcy_velocity[a] = 0.0;
        }
        /* neglect gravitational effect in gas, for now */
      }
    } else {
      for (a = 0; a < VIM; a++) {
        pmv->liq_darcy_velocity[a] = 0.0;
        pmv->gas_darcy_velocity[a] = 0.0;
        for (b = 0; b < VIM; b++) {
          pmv->liq_darcy_velocity[a] -= mp->perm_tensor[a][b] * mp->rel_liq_perm *
                                        (pow(fv->grad_p_liq[b], n_pow) - mp->momentum_source[b]);
          if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
            pmv->gas_darcy_velocity[a] -=
                mp->perm_tensor[a][b] * mp->rel_gas_perm * fv->grad_p_gas[b];
          }
          /* neglect gravitational effect in gas, for now */
        }
      }
    }
    /* add convection fluxes into the total diffusion flux term
     * here vconv is the velocity of the solid phase */
    for (a = 0; a < VIM; a++) {
      pmv->rel_mass_flux[i_pl][a] = 0.0;
      pmv->rel_mass_flux[i_pg][a] = 0.0;
    }

    for (a = 0; a < VIM; a++) {
      pmv->rel_mass_flux[i_pl][a] =
          pmv->gas_darcy_velocity[a] * pmv->gas_density_solvents[i_pl] +
          pmv->liq_darcy_velocity[a] * pmv->liq_Xvol_solvents[i_pl] * mp->density;
      if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
        pmv->rel_mass_flux[i_pg][a] =
            pmv->gas_darcy_velocity[a] * pmv->gas_density_solvents[i_pg] +
            pmv->liq_darcy_velocity[a] * pmv->liq_Xvol_solvents[i_pg] * mp->density;
      }
    }

    /* add in diffusive flux in the gas phase only !!
     * the diffusive flux is proportional to the gradient in gas phase
     * concentration times the diffusion coefficient.
     * Here we get the gradient in gas phase concentration by the
     * chain rule - grad(cgas) = sum over variables (grad(var) * d(cgas)/d(var))
     * where the variables summed over are pressure of liq and gas, porosity, and temp
     */
    for (w = 0; w < MAX_PMV; w++) {
      if (w != i_pore) { /* no capacity term for solid phase */
        for (a = 0; a < VIM; a++) {
          /*
           * convert pressure and porosity gradients into gradients
           * of gas phase concentration
           */
          pmv->rel_mass_flux[w][a] -=
              mp->porous_diffusivity[w] *
              (pmv->d_gas_density_solvents[w][POR_LIQ_PRES] * fv->grad_p_liq[a] +
               pmv->d_gas_density_solvents[w][POR_GAS_PRES] * fv->grad_p_gas[a] +
               pmv->d_gas_density_solvents[w][POR_POROSITY] * fv->grad_porosity[a]);
        }
      }
    }
  } else {
    GOMA_EH(GOMA_ERROR, "Unimplemented flux constitutive relation in continuous media.");
  }

  /* add in diffusive flux in the gas phase only !!
   * the diffusive flux is proportional to the gradient in gas phase
   * concentration times the diffusion coefficient.
   * Here we get the gradient in gas phase concentration by the
   * chain rule - grad(cgas) = sum over variables (grad(var) * d(cgas)/d(var))
   * where the variables summed over are pressure of liq and gas, porosity, and temp
   *
   * For the energy equation, multiply the diffusive flux of water in the gas times
   * the enthalpy of water in the gas (pmv->enthalpy[1]), and the diffusive flux of
   * air in the gas times the enthalpy of air in the gas (pmv->enthalpy[2])
   */

  if (cr->PorousFluxModel == DARCY_FICKIAN) {

    /* In the binary system currently considered (water & air), the diffusive flux
     * of the air is equal and opposite to the diffusive flux of vapor (ie, water
     * in the gas phase - exploit this fact here and do not loop
     *
     *  for (w = 0; w < MAX_PMV; w++){
     *    if (w != i_pore && w != i_pe) {
     */
    w = 0;
    for (a = 0; a < VIM; a++) {
      /*
       * convert pressure and porosity gradients into gradients
       * of gas phase concentration
       */
      difflux[w][a] = -mp->porous_diffusivity[w] *
                      (pmv->d_gas_density_solvents[w][POR_LIQ_PRES] * fv->grad_p_liq[a] +
                       pmv->d_gas_density_solvents[w][POR_GAS_PRES] * fv->grad_p_gas[a] +
                       pmv->d_gas_density_solvents[w][POR_POROSITY] * fv->grad_porosity[a] +
                       pmv->d_gas_density_solvents[w][POR_TEMP] * fv->grad_T[a]);
      difflux[w + 1][a] = -difflux[w][a];
      pmv->rel_mass_flux[w][a] += difflux[w][a];
      if (pd->e[pg->imtrx][R_POR_GAS_PRES])
        pmv->rel_mass_flux[i_pg][a] += difflux[w + 1][a];
      if (pd->e[pg->imtrx][R_POR_ENERGY])
        pmv->rel_mass_flux[i_pe][a] += difflux[w][a] * (pmv->enthalpy[1] - pmv->enthalpy[2]);
    }

  } else {
    GOMA_EH(GOMA_ERROR, "Unimplemented flux constitutive relation in continuous media.");
  }

  if (af->Assemble_Jacobian) {
    var = POR_LIQ_PRES;
    for (w = 0; w < MAX_PMV; w++) {
      for (a = 0; a < VIM; a++) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          for (w1 = 0; w1 < MAX_PMV; w1++) {
            pmv->d_rel_mass_flux_dpmv[w][a][w1][j] = 0.;
          }
        }
      }
    }
    if (cr->PorousFluxModel == DARCY_FICKIAN) {
      for (w = 0; w < MAX_PMV; w++) {
        if (w != i_pore) { /* no diffusion term for solid phase */

          for (a = 0; a < VIM; a++) {
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              for (w1 = 0; w1 < MAX_PMV; w1++) {
                pmv->d_rel_mass_flux_dpmv[w][a][w1][j] +=
                    bf[var]->phi[j] * pmv->gas_darcy_velocity[a] *
                    pmv->d_gas_density_solvents[w][POR_LIQ_PRES + w1];
                pmv->d_rel_mass_flux_dpmv[w][a][w1][j] +=
                    bf[var]->phi[j] * pmv->liq_darcy_velocity[a] * mp->density *
                    pmv->d_liq_Xvol_solvents[w][POR_LIQ_PRES + w1];
              }
            }
          }

          if (mp->PermeabilityModel != K_TENSOR && mp->PermeabilityModel != ORTHOTROPIC &&
              mp->PermeabilityModel != SM_TENSOR && mp->PermeabilityModel != KC_TENSOR) {
            for (a = 0; a < VIM; a++) {
              for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                for (w1 = 0; w1 < MAX_PMV; w1++) {
                  pmv->d_rel_mass_flux_dpmv[w][a][w1][j] -=
                      bf[var]->phi[j] * pmv->liq_Xvol_solvents[w] * mp->density *
                      (fv->grad_p_liq[a] - mp->momentum_source[a]) *
                      (mp->permeability * mp->d_rel_liq_perm[POR_LIQ_PRES + w1] +
                       mp->d_permeability[POR_LIQ_PRES + w1] * mp->rel_liq_perm);

                  if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
                    pmv->d_rel_mass_flux_dpmv[w][a][w1][j] -=
                        bf[var]->phi[j] * pmv->gas_density_solvents[w] * fv->grad_p_gas[a] *
                        (mp->permeability * mp->d_rel_gas_perm[POR_LIQ_PRES + w1] +
                         mp->d_permeability[POR_LIQ_PRES + w1] * mp->rel_gas_perm);
                  }
                }
                /* sensitivity to liquid pressure gradient */
                pmv->d_rel_mass_flux_dpmv[w][a][i_pl][j] -= mp->permeability * mp->rel_liq_perm *
                                                            bf[var]->grad_phi[j][a] *
                                                            pmv->liq_Xvol_solvents[w] * mp->density;

                /* sensitivity to gas pressure gradient and permeability */
                if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
                  pmv->d_rel_mass_flux_dpmv[w][a][i_pg][j] -=
                      mp->permeability * pmv->gas_density_solvents[w] * mp->rel_gas_perm *
                      bf[var]->grad_phi[j][a];
                }
              }
            }
          } else { /* Permeability is a tensor, and it is constant */
            for (a = 0; a < VIM; a++) {
              for (b = 0; b < VIM; b++) {
                for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                  for (w1 = 0; w1 < MAX_PMV; w1++) {
                    pmv->d_rel_mass_flux_dpmv[w][a][w1][j] -=
                        bf[var]->phi[j] *
                        (mp->perm_tensor[a][b] * mp->d_rel_liq_perm[POR_LIQ_PRES + w1] *
                         pmv->liq_Xvol_solvents[w] * mp->density *
                         (fv->grad_p_liq[b] - mp->momentum_source[b]));

                    pmv->d_rel_mass_flux_dpmv[w][a][w1][j] -=
                        bf[var]->phi[j] *
                        (mp->d_perm_tensor[a][b][POR_LIQ_PRES + w1] * mp->rel_liq_perm *
                         pmv->liq_Xvol_solvents[w] * mp->density *
                         (fv->grad_p_liq[b] - mp->momentum_source[b]));

                    if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
                      pmv->d_rel_mass_flux_dpmv[w][a][w1][j] -=
                          bf[var]->phi[j] * mp->perm_tensor[a][b] * pmv->gas_density_solvents[w] *
                          mp->d_rel_gas_perm[POR_LIQ_PRES + w1] * fv->grad_p_gas[b];

                      pmv->d_rel_mass_flux_dpmv[w][a][w1][j] -=
                          bf[var]->phi[j] *
                          (mp->d_perm_tensor[a][b][POR_LIQ_PRES + w1] * mp->rel_gas_perm *
                           pmv->liq_Xvol_solvents[w] * mp->density * (fv->grad_p_gas[b]));
                    }
                  }
                  /* sensitivity to liquid pressure gradient */
                  pmv->d_rel_mass_flux_dpmv[w][a][i_pl][j] -=
                      mp->perm_tensor[a][b] * mp->rel_liq_perm * bf[var]->grad_phi[j][b] *
                      pmv->liq_Xvol_solvents[w] * mp->density;

                  /* sensitivity to gas pressure gradient and permeability */
                  if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
                    pmv->d_rel_mass_flux_dpmv[w][a][i_pg][j] -=
                        mp->perm_tensor[a][b] * pmv->gas_density_solvents[w] * mp->rel_gas_perm *
                        bf[var]->grad_phi[j][b];
                  }
                }
              }
            }
          }
        }
      }
    } else if (cr->PorousFluxModel == POWERLAW_DARCY_FICKIAN) {
      for (w = 0; w < MAX_PMV; w++) {
        if (w != i_pore) { /* no diffusion term for solid phase */

          for (a = 0; a < VIM; a++) {
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              for (w1 = 0; w1 < MAX_PMV; w1++) {
                pmv->d_rel_mass_flux_dpmv[w][a][w1][j] +=
                    bf[var]->phi[j] * pmv->gas_darcy_velocity[a] *
                    pmv->d_gas_density_solvents[w][POR_LIQ_PRES + w1];
              }
            }
          }

          if (mp->PermeabilityModel != K_TENSOR && mp->PermeabilityModel != ORTHOTROPIC &&
              mp->PermeabilityModel != SM_TENSOR && mp->PermeabilityModel != KC_TENSOR) {
            for (a = 0; a < VIM; a++) {
              for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                for (w1 = 0; w1 < MAX_PMV; w1++) {
                  pmv->d_rel_mass_flux_dpmv[w][a][w1][j] -=
                      bf[var]->phi[j] * pmv->liq_Xvol_solvents[w] * mp->density *
                      (pow(fv->grad_p_liq[a], n_pow) - mp->momentum_source[a]) *
                      (mp->permeability * mp->d_rel_liq_perm[POR_LIQ_PRES + w1] +
                       mp->d_permeability[POR_LIQ_PRES + w1] * mp->rel_liq_perm);

                  if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
                    pmv->d_rel_mass_flux_dpmv[w][a][w1][j] -=
                        bf[var]->phi[j] * pmv->gas_density_solvents[w] * fv->grad_p_gas[a] *
                        (mp->permeability * mp->d_rel_gas_perm[POR_LIQ_PRES + w1] +
                         mp->d_permeability[POR_LIQ_PRES + w1] * mp->rel_gas_perm);
                  }
                }
                /* sensitivity to liquid pressure gradient */
                pmv->d_rel_mass_flux_dpmv[w][a][i_pl][j] -=
                    mp->permeability * mp->rel_liq_perm * n_pow *
                    pow(fv->grad_p_liq[a], n_pow - 1) * bf[var]->grad_phi[j][a] *
                    pmv->liq_Xvol_solvents[w] * mp->density;

                /* sensitivity to gas pressure gradient and permeability */
                if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
                  pmv->d_rel_mass_flux_dpmv[w][a][i_pg][j] -=
                      mp->permeability * pmv->gas_density_solvents[w] * mp->rel_gas_perm *
                      bf[var]->grad_phi[j][a];
                }
              }
            }
          } else { /* Permeability is a tensor, and it is constant */
            for (a = 0; a < VIM; a++) {
              for (b = 0; b < VIM; b++) {
                for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                  for (w1 = 0; w1 < MAX_PMV; w1++) {
                    pmv->d_rel_mass_flux_dpmv[w][a][w1][j] -=
                        bf[var]->phi[j] *
                        (mp->perm_tensor[a][b] * mp->d_rel_liq_perm[POR_LIQ_PRES + w1] *
                         pmv->liq_Xvol_solvents[w] * mp->density *
                         (pow(fv->grad_p_liq[b], n_pow) - mp->momentum_source[b]));

                    pmv->d_rel_mass_flux_dpmv[w][a][w1][j] -=
                        bf[var]->phi[j] *
                        (mp->d_perm_tensor[a][b][POR_LIQ_PRES + w1] * mp->rel_liq_perm *
                         pmv->liq_Xvol_solvents[w] * mp->density *
                         (pow(fv->grad_p_liq[b], n_pow) - mp->momentum_source[b]));

                    if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
                      pmv->d_rel_mass_flux_dpmv[w][a][w1][j] -=
                          bf[var]->phi[j] * mp->perm_tensor[a][b] * pmv->gas_density_solvents[w] *
                          mp->d_rel_gas_perm[POR_LIQ_PRES + w1] * fv->grad_p_gas[b];

                      pmv->d_rel_mass_flux_dpmv[w][a][w1][j] -=
                          bf[var]->phi[j] * mp->d_perm_tensor[a][b][POR_LIQ_PRES + w1] *
                          pmv->gas_density_solvents[w] * mp->rel_gas_perm * fv->grad_p_gas[b];
                    }
                  }
                  /* sensitivity to liquid pressure gradient */
                  pmv->d_rel_mass_flux_dpmv[w][a][i_pl][j] -=
                      mp->perm_tensor[a][b] * mp->rel_liq_perm * n_pow *
                      pow(fv->grad_p_liq[a], n_pow - 1) * bf[var]->grad_phi[j][b] *
                      pmv->liq_Xvol_solvents[w] * mp->density;

                  /* sensitivity to gas pressure gradient and permeability */
                  if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
                    pmv->d_rel_mass_flux_dpmv[w][a][i_pg][j] -=
                        mp->perm_tensor[a][b] * pmv->gas_density_solvents[w] * mp->rel_gas_perm *
                        bf[var]->grad_phi[j][b];
                  }
                }
              }
            }
          }
        }
      }
    } else {
      GOMA_EH(GOMA_ERROR, "unimplemented darcy flux model");
    }

    if (cr->PorousFluxModel == DARCY_FICKIAN) {
      for (w = 0; w < MAX_PMV; w++) {
        if (w != i_pore) { /* no diffusion term for solid phase */
          for (a = 0; a < VIM; a++) {
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              difflux_term = pmv->d_gas_density_solvents[w][POR_LIQ_PRES] *
                             mp->porous_diffusivity[w] * bf[var]->grad_phi[j][a];
              pmv->d_rel_mass_flux_dpmv[w][a][i_pl][j] -= difflux_term;
              if (pd->e[pg->imtrx][R_POR_ENERGY])
                pmv->d_rel_mass_flux_dpmv[i_pe][a][i_pl][j] -= difflux_term * pmv->enthalpy[1];
              for (w2 = 0; w2 < MAX_PMV; w2++) {
                difflux_term =
                    bf[var]->phi[j] *
                    (mp->porous_diffusivity[w] * fv->grad_p_liq[a] *
                         pmv->d_d_gas_density_solvents[w][POR_LIQ_PRES][POR_LIQ_PRES + w2] +
                     mp->d_porous_diffusivity[w][POR_LIQ_PRES + w2] *
                         pmv->d_gas_density_solvents[w][POR_LIQ_PRES] * fv->grad_p_liq[a]);
                pmv->d_rel_mass_flux_dpmv[w][a][w2][j] -= difflux_term;
                if (pd->e[pg->imtrx][R_POR_ENERGY])
                  pmv->d_rel_mass_flux_dpmv[i_pe][a][w2][j] -=
                      difflux_term * pmv->enthalpy[1] +
                      difflux[i_pl][a] * pmv->d_enthalpy[1][POR_LIQ_PRES + w2];
              }

              if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
                difflux_term = pmv->d_gas_density_solvents[w][POR_GAS_PRES] *
                               mp->porous_diffusivity[w] * bf[var]->grad_phi[j][a];
                pmv->d_rel_mass_flux_dpmv[w][a][i_pg][j] -= difflux_term;
                if (pd->e[pg->imtrx][R_POR_ENERGY])
                  pmv->d_rel_mass_flux_dpmv[i_pe][a][i_pg][j] -= difflux_term * pmv->enthalpy[2];
                for (w2 = 0; w2 < MAX_PMV; w2++) {
                  difflux_term =
                      bf[var]->phi[j] *
                      (mp->porous_diffusivity[w] * fv->grad_p_gas[a] *
                           pmv->d_d_gas_density_solvents[w][POR_GAS_PRES][POR_LIQ_PRES + w2] +
                       mp->d_porous_diffusivity[w][POR_LIQ_PRES + w2] *
                           pmv->d_gas_density_solvents[w][POR_GAS_PRES] * fv->grad_p_gas[a]);
                  pmv->d_rel_mass_flux_dpmv[w][a][w2][j] -= difflux_term;
                  if (pd->e[pg->imtrx][R_POR_ENERGY])
                    pmv->d_rel_mass_flux_dpmv[i_pe][a][w2][j] -=
                        difflux_term * pmv->enthalpy[2] +
                        difflux[i_pg][a] * pmv->d_enthalpy[2][POR_LIQ_PRES + w2];
                }
              }

              if (pd->e[pg->imtrx][R_POR_POROSITY]) {
                pmv->d_rel_mass_flux_dpmv[w][a][i_pore][j] -=
                    pmv->d_gas_density_solvents[w][POR_POROSITY] * mp->porous_diffusivity[w] *
                    bf[var]->grad_phi[j][a];
                for (w2 = 0; w2 < MAX_PMV; w2++) {
                  pmv->d_rel_mass_flux_dpmv[w][a][w2][j] -=
                      bf[var]->phi[j] *
                      (mp->porous_diffusivity[w] * fv->grad_porosity[a] *
                           pmv->d_d_gas_density_solvents[w][POR_POROSITY][POR_LIQ_PRES + w2] +
                       mp->d_porous_diffusivity[w][POR_LIQ_PRES + w2] *
                           pmv->d_gas_density_solvents[w][POR_POROSITY] * fv->grad_porosity[a]);
                }
              }
              /* Additional terms for energy equation */
              if (pd->e[pg->imtrx][R_POR_ENERGY]) {
                pmv->d_rel_mass_flux_dpmv[i_pe][a][i_pore][j] -=
                    mp->thermal_conductivity * bf[var]->grad_phi[j][a];
                for (w2 = 0; w2 < MAX_PMV; w2++) {
                  pmv->d_rel_mass_flux_dpmv[i_pe][a][i_pore][j] -=
                      mp->d_thermal_conductivity[POR_LIQ_PRES + w2] * fv->grad_T[a];
                }
              }

            } /* j<ei[pg->imtrx]->dof[var] */
          }   /* a<VIM */
        }     /* w!=i_pore */
      }       /* w<MAX_PMV */

    } /* DARCY_FICKIAN */

    else if (cr->PorousFluxModel == POWERLAW_DARCY_FICKIAN) {
      for (w = 0; w < MAX_PMV; w++) {
        if (w != i_pore) { /* no diffusion term for solid phase */
          for (a = 0; a < VIM; a++) {
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              pmv->d_rel_mass_flux_dpmv[w][a][i_pl][j] -=
                  pmv->d_gas_density_solvents[w][POR_LIQ_PRES] * mp->porous_diffusivity[w] *
                  bf[var]->grad_phi[j][a];
              for (w2 = 0; w2 < MAX_PMV; w2++) {
                pmv->d_rel_mass_flux_dpmv[w][a][w2][j] -=
                    bf[var]->phi[j] *
                    (mp->porous_diffusivity[w] * pow(fv->grad_p_liq[a], n_pow) *
                         pmv->d_d_gas_density_solvents[w][POR_LIQ_PRES][POR_LIQ_PRES + w2] +
                     mp->d_porous_diffusivity[w][POR_LIQ_PRES + w2] *
                         pmv->d_gas_density_solvents[w][POR_LIQ_PRES] *
                         pow(fv->grad_p_liq[a], n_pow));
              }

              if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
                pmv->d_rel_mass_flux_dpmv[w][a][i_pg][j] -=
                    pmv->d_gas_density_solvents[w][POR_GAS_PRES] * mp->porous_diffusivity[w] *
                    bf[var]->grad_phi[j][a];
                for (w2 = 0; w2 < MAX_PMV; w2++) {
                  pmv->d_rel_mass_flux_dpmv[w][a][w2][j] -=
                      bf[var]->phi[j] *
                      (mp->porous_diffusivity[w] * fv->grad_p_gas[a] *
                           pmv->d_d_gas_density_solvents[w][POR_GAS_PRES][POR_LIQ_PRES + w2] +
                       mp->d_porous_diffusivity[w][POR_LIQ_PRES + w2] *
                           pmv->d_gas_density_solvents[w][POR_GAS_PRES] * fv->grad_p_gas[a]);
                }
              }

              if (pd->e[pg->imtrx][R_POR_POROSITY]) {
                pmv->d_rel_mass_flux_dpmv[w][a][i_pore][j] -=
                    pmv->d_gas_density_solvents[w][POR_POROSITY] * mp->porous_diffusivity[w] *
                    bf[var]->grad_phi[j][a];
                for (w2 = 0; w2 < MAX_PMV; w2++) {
                  pmv->d_rel_mass_flux_dpmv[w][a][w2][j] -=
                      bf[var]->phi[j] *
                      (mp->porous_diffusivity[w] * fv->grad_porosity[a] *
                           pmv->d_d_gas_density_solvents[w][POR_POROSITY][POR_LIQ_PRES + w2] +
                       mp->d_porous_diffusivity[w][POR_LIQ_PRES + w2] *
                           pmv->d_gas_density_solvents[w][POR_POROSITY] * fv->grad_porosity[a]);
                }
              }
            }
          }
        }
      }
    } else {
      GOMA_EH(GOMA_ERROR, "Unimplemented flux constitutive relation in porous media. "
                          "Check Porous Diffusion Constitutive Equation Card");
    }

    for (b = 0; b < pd->Num_Dim; b++) {
      var = MESH_DISPLACEMENT1 + b;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          for (a = 0; a < VIM; a++) {
            for (w = 0; w < MAX_PMV; w++) {
              pmv->d_rel_mass_flux_dmesh[w][a][b][j] = 0.;
            }
          }
          for (a = 0; a < VIM; a++) {
            if (cr->PorousFluxModel == DARCY_FICKIAN) {
              for (w = 0; w < MAX_PMV; w++) {
                if (w != i_pore) { /* no diffusion term for solid phase */
                  if (mp->PermeabilityModel != K_TENSOR && mp->PermeabilityModel != ORTHOTROPIC &&
                      mp->PermeabilityModel != SM_TENSOR && mp->PermeabilityModel != KC_TENSOR) {
                    pmv->d_rel_mass_flux_dmesh[w][a][b][j] +=
                        -mp->permeability * mp->rel_liq_perm * pmv->liq_Xvol_solvents[w] *
                        fv->d_grad_p_liq_dmesh[a][b][j] * mp->density;
                    if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
                      pmv->d_rel_mass_flux_dmesh[w][a][b][j] +=
                          -mp->permeability * mp->rel_gas_perm * pmv->gas_density_solvents[w] *
                          fv->d_grad_p_gas_dmesh[a][b][j];
                    }
                  } else {
                    for (p = 0; p < VIM; p++) {
                      pmv->d_rel_mass_flux_dmesh[w][a][b][j] +=
                          -mp->perm_tensor[a][p] * mp->rel_liq_perm * pmv->liq_Xvol_solvents[w] *
                          fv->d_grad_p_liq_dmesh[p][b][j] * mp->density;
                      if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
                        pmv->d_rel_mass_flux_dmesh[w][a][b][j] +=
                            -mp->perm_tensor[a][p] * mp->rel_gas_perm *
                            pmv->gas_density_solvents[w] * fv->d_grad_p_gas_dmesh[p][b][j];
                      }
                      if (mp->PermeabilityModel == ORTHOTROPIC ||
                          mp->PermeabilityModel == SM_TENSOR ||
                          mp->PermeabilityModel == KC_TENSOR) {
                        pmv->d_rel_mass_flux_dmesh[w][a][b][j] +=
                            -mp->d_perm_tensor_dx[a][p][b][j] * mp->rel_liq_perm *
                            pmv->liq_Xvol_solvents[w] * fv->grad_p_liq[p] * mp->density;
                        if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
                          pmv->d_rel_mass_flux_dmesh[w][a][b][j] +=
                              -mp->d_perm_tensor_dx[a][p][b][j] * mp->rel_gas_perm *
                              pmv->gas_density_solvents[w] * fv->grad_p_gas[p];
                        }
                      }
                    }
                  }
                } /* w!=i_pore */
              }   /* w<MAX_PMV */

              /* Additional terms for energy equation */
              if (pd->e[pg->imtrx][R_POR_GAS_PRES])
                pmv->d_rel_mass_flux_dmesh[i_pe][a][b][j] +=
                    -mp->thermal_conductivity * fv->d_grad_T_dmesh[a][b][j];

            } /* DARCY_FICKIAN */

            else if (cr->PorousFluxModel == POWERLAW_DARCY_FICKIAN) {
              for (w = 0; w < MAX_PMV; w++) {
                if (w != i_pore) { /* no diffusion term for solid phase */
                  if (mp->PermeabilityModel != K_TENSOR && mp->PermeabilityModel != ORTHOTROPIC &&
                      mp->PermeabilityModel != SM_TENSOR && mp->PermeabilityModel != KC_TENSOR) {
                    pmv->d_rel_mass_flux_dmesh[w][a][b][j] +=
                        -mp->permeability * mp->rel_liq_perm * pmv->liq_Xvol_solvents[w] * n_pow *
                        pow(fv->grad_p_liq[a], n_pow - 1) * fv->d_grad_p_liq_dmesh[a][b][j] *
                        mp->density;
                    if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
                      pmv->d_rel_mass_flux_dmesh[w][a][b][j] +=
                          -mp->permeability * mp->rel_gas_perm * pmv->gas_density_solvents[w] *
                          fv->d_grad_p_gas_dmesh[a][b][j];
                    }
                  } else {
                    for (p = 0; p < VIM; p++) {
                      pmv->d_rel_mass_flux_dmesh[w][a][b][j] +=
                          -mp->perm_tensor[a][p] * mp->rel_liq_perm * pmv->liq_Xvol_solvents[w] *
                          n_pow * pow(fv->grad_p_liq[a], n_pow - 1) *
                          fv->d_grad_p_liq_dmesh[p][b][j] * mp->density;
                      if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
                        pmv->d_rel_mass_flux_dmesh[w][a][b][j] +=
                            -mp->perm_tensor[a][p] * mp->rel_gas_perm *
                            pmv->gas_density_solvents[w] * fv->d_grad_p_gas_dmesh[p][b][j];
                      }
                      if (mp->PermeabilityModel == ORTHOTROPIC ||
                          mp->PermeabilityModel == SM_TENSOR ||
                          mp->PermeabilityModel == KC_TENSOR) {
                        pmv->d_rel_mass_flux_dmesh[w][a][b][j] +=
                            -mp->d_perm_tensor_dx[a][p][b][j] * mp->rel_liq_perm *
                            pmv->liq_Xvol_solvents[w] * n_pow * pow(fv->grad_p_liq[a], n_pow - 1) *
                            fv->grad_p_liq[p] * mp->density;
                        if (pd->e[pg->imtrx][R_POR_GAS_PRES]) {
                          pmv->d_rel_mass_flux_dmesh[w][a][b][j] +=
                              -mp->d_perm_tensor_dx[a][p][b][j] * mp->rel_gas_perm *
                              pmv->gas_density_solvents[w] * fv->grad_p_gas[p];
                        }
                      }
                    }
                  }
                }
              }
            } else {
              GOMA_EH(GOMA_ERROR, "Unimplemented Porous Flux model");
            }

            if (cr->PorousFluxModel == DARCY_FICKIAN) {
              for (w = 0; w < MAX_PMV; w++) {
                if (w != i_pore) { /* no diffusion term for solid phase */
                  pmv->d_rel_mass_flux_dmesh[w][a][b][j] -=
                      mp->porous_diffusivity[w] * (pmv->d_gas_density_solvents[w][POR_LIQ_PRES] *
                                                       fv->d_grad_p_liq_dmesh[a][b][j] +
                                                   pmv->d_gas_density_solvents[w][POR_GAS_PRES] *
                                                       fv->d_grad_p_gas_dmesh[a][b][j] +
                                                   pmv->d_gas_density_solvents[w][POR_POROSITY] *
                                                       fv->d_grad_porosity_dmesh[a][b][j]);
                }
              }
            } /* DARCY_FICKIAN */
          }   /* a<VIM */
        }     /* j<ei[pg->imtrx]->dof[var] */
      }       /*if ( pd->v[pg->imtrx][var] ) */
    }         /*  b<pd->Num_Dim */

    var = TEMPERATURE;
    if (pd->v[pg->imtrx][var] && !pd->e[pg->imtrx][R_POR_ENERGY]) {
      if (cr->PorousFluxModel == DARCY_FICKIAN) {
        /* NOTE - the gas concentrations, etc should also be functions of temperature */
        for (w = 0; w < MAX_PMV; w++) {
          if (w != i_pore) { /* no diffusion term for solid phase */
            for (a = 0; a < VIM; a++) {
              for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                pmv->d_rel_mass_flux_dT[w][a][j] = 0.;
                pmv->d_rel_mass_flux_dT[w][a][j] -=
                    mp->d_porous_diffusivity[w][var] * bf[var]->phi[j] *
                    (pmv->d_gas_density_solvents[w][POR_LIQ_PRES] * fv->grad_p_liq[a] +
                     pmv->d_gas_density_solvents[w][POR_GAS_PRES] * fv->grad_p_gas[a] +
                     pmv->d_gas_density_solvents[w][POR_POROSITY] * fv->grad_porosity[a]);
              }
            }
          }
        }
      }
    }
  }

  return;
} /* end  load_mass_flux  */

/*******************************************************************************/
/*******************************************************************************/

void porous_pressure(double *func, double d_func[], int eb_mat_solid, int eb_mat_fluid)

/***************************************************************************
 *
 *  Function which equates the pressure of a continuous fluid medium
 *  with the pressure in a porous saturated medium.
 *
 *
 *  This function calculates the expression:
 *
 *          func = P_continuum - P_porous
 *
 *  It does this by adding in P_continuum from the continuum side and
 *  P_porous from the porous media side of the interface.
 *
 *			 Author: p. r. schunk (9/28/98)
 **************************************************************************/
{
  if (af->Assemble_LSA_Mass_Matrix)
    return;

  /***************************** SOLID SIDE *******************************/
  /*
   *  If current material is the porous solid phase, add on the pressure
   */
  if (Current_EB_ptr->Elem_Blk_Id == eb_mat_solid) {
    if (mp->PorousMediaType == POROUS_SATURATED || mp->PorousMediaType == POROUS_UNSATURATED ||
        mp->PorousMediaType == POROUS_TWO_PHASE) {
      *func = -fv->p_liq;
      if (af->Assemble_Jacobian) {
        d_func[POR_LIQ_PRES] = -1.0;
      }
    }
  }

  /***************************** FLUID SIDE *******************************/
  /*
   *  If current material is the fluid phase, add on the fluid pressure
   */
  else if (Current_EB_ptr->Elem_Blk_Id == eb_mat_fluid) {
    *func = fv->P;
    if (af->Assemble_Jacobian) {
      d_func[PRESSURE] = 1.0;
    }
  } else {
    GOMA_EH(GOMA_ERROR, "Penetration called with incorrect block id");
  }
  return;
} /* end  porous_pressure  */
/******************************************************************************/
/******************************************************************************/
/*******************************************************************************/

void porous_pressure_lub(double *func,
                         double d_func[],
                         const int id_side,
                         double xi[DIM],
                         const Exo_DB *exo,
                         const double scale)

/***************************************************************************
 *
 *  Function which equates the pressure of a lubrication layer
 *  with the pressure in a porous medium.
 *
 *
 *  This function calculates the expression:
 *
 *          func = P_lub - P_porous
 *
 *  It does this by adding in P_lub from the lubrication side and
 *  P_porous from the porous media side of the interface.
 *
 *			 Author: p. r. schunk (8/09/11)
 **************************************************************************/
{
  int *n_dof = NULL;

  if (af->Assemble_LSA_Mass_Matrix)
    return;

  /*
   * Prepare geometry
   */
  int dof_map[MDE];
  n_dof = (int *)array_alloc(1, MAX_VARIABLE_TYPES, sizeof(int));
  lubrication_shell_initialize(n_dof, dof_map, id_side, xi, exo, 1);

  if (mp->PorousMediaType == POROUS_SATURATED || mp->PorousMediaType == POROUS_UNSATURATED ||
      mp->PorousMediaType == POROUS_TWO_PHASE) {
    *func = -fv->p_liq;
    if (af->Assemble_Jacobian) {
      d_func[POR_LIQ_PRES] = -1.0;
    }
  }
  *func += fv->lubp;
  if (af->Assemble_Jacobian) {
    d_func[LUBP] = 1.0;
  }

  return;
} /* end  porous_pressure_lub  */

/*******************************************************************************/
/******************************************************************************/

void sat_darcy_continuous_bc(double func[],
                             double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                             const double tt,
                             const double dt,
                             const double time_value,
                             const int eb_mat_solid,
                             const int eb_mat_fluid,
                             const double length_scale, /* length scale for "attached" */
                             const double vnormal)      /* attachment velocity */

/******************************************************************************
 *
 * sat_darcy_continuous_bc() - kinematic BC for adjacent porous media & fluid
 *
 * Function which evaluates the kinematic boundary condition
 * for adjacent porous medium and fluid phases
 *			 Author: P. R. Schunk (sometime in 1999)
 *
 * Input
 * -------
 *
 *   tt - parameter to vary time integration  method.
 *   dt - current time step size
 *   eb_mat_solid - element block id of porous  phase
 *   eb_mat_fluid - element block id of gas phase
 *   vnormal --normal attachment velocity for level-set applications
 *   length_scale --length scale over which a meniscus is declared "attached" to
 *                  a porous fluid surface, viz. if within, darcy flux is added.
 *

 *****************************************************************************/
{
  int j_id, var, jvar, a, b, p, i_pore, dim;
  double phi_j, factor = 1.0e12, d_factor[3] = {0., 0., 0.};

  double vconv[MAX_PDIM];     /*Calculated convection velocity */
  double vconv_old[MAX_PDIM]; /*Calculated convection velocity at previous time*/
  CONVECTION_VELOCITY_DEPENDENCE_STRUCT d_vconv_struct;
  CONVECTION_VELOCITY_DEPENDENCE_STRUCT *d_vconv = &d_vconv_struct;

  double x_dot[DIM]; /* mesh velocity vector   */
  int err;
  double attached, near;

  if (af->Assemble_LSA_Mass_Matrix)
    return;

  dim = pd->Num_Dim;
  i_pore = 2;

  /* When applying various terms from solid and liquid side, establish the booleans
   * "attached" and "near" by computing distances and level-set function values.
   * We are distinguishing three situations (note we can only evaluate in the fluid phase)
   *  if(!attached && !near)
   *   n.(v-v_s) = 0
   *  else if (!attached && near)
   *   n.(v-v_s) = v_normal
   *  else if (attached)
   *   n.(v-v_s) = darcy_velocity
   *  else (if not level set active)
   *   n.(v-v_s) = darcy_velocity
   *
   * Note that from the fluid side we can only establish nearness.  From the solid side we need
   * to establish attachement/contact so we know whether to apply the darcy velocity
   */

  attached = 1.0;
  near = 0.; /*default case */

  /***************************** SOLID SIDE *******************************/
  /*
   *  If current material is the solid phase, calculate normal flux
   *  of solvent through the porous medium
   */
  if (Current_EB_ptr->Elem_Blk_Id == eb_mat_solid) {
    /* add convection fluxes into the total diffusion flux term
     * here vconv is the velocity of the solid phase */

    if (mp->PorousMediaType == CONTINUOUS)
      GOMA_EH(GOMA_ERROR, "DARCY_CONTINUOUS requires a porous phase: Are your mat IDs switched?");
    if (mp->PorousMediaType == POROUS_UNSATURATED) {
      factor = mp->rel_liq_perm;
      d_factor[0] = mp->d_rel_liq_perm[POR_LIQ_PRES];
    }
    if (mp->PorousMediaType == POROUS_TWO_PHASE) {
      factor = mp->rel_liq_perm;
      d_factor[0] = mp->d_rel_liq_perm[POR_LIQ_PRES];
    }
    if (mp->PorousMediaType == POROUS_SATURATED)
      factor = 1 / mp->viscosity;

      /* establish attachment if this is a level set problem */
      /* Now establish attachment and wetted */
#if 0
    mn_liq = map_mat_index ( eb_mat_fluid );

    if( pd_glob[mn_liq]->v[pg->imtrx][LS] )
      {
	load_lsi_adjmatr( 4.*length_scale );
	attached = 1.0 - lsi->H;
      }
#endif

    /*
     * Get the lagrangian velocity of the grid
     */
    if ((pd->MeshMotion == LAGRANGIAN || cr->MeshMotion == DYNAMIC_LAGRANGIAN) &&
        pd->MeshInertia == 1) {
      err = belly_flop(elc->lame_mu);
      GOMA_EH(err, "error in belly flop");
      if (neg_elem_volume)
        return;
    }
    err = get_convection_velocity(vconv, vconv_old, d_vconv, dt, tt);
    GOMA_EH(err, "Error in calculating effective convection velocity");

    if (mp->PermeabilityModel != K_TENSOR && mp->PermeabilityModel != ORTHOTROPIC &&
        mp->PermeabilityModel != SM_TENSOR && mp->PermeabilityModel != KC_TENSOR) {
      for (a = 0; a < VIM; a++) {
        pmv->liq_darcy_velocity[a] =
            -factor * mp->permeability * (fv->grad_p_liq[a] - mp->momentum_source[a]);
      }
    } else {
      for (a = 0; a < VIM; a++) {
        pmv->liq_darcy_velocity[a] = 0.0;
        for (b = 0; b < VIM; b++) {
          pmv->liq_darcy_velocity[a] +=
              -factor * mp->perm_tensor[a][b] * (fv->grad_p_liq[b] - mp->momentum_source[b]);
        }
      }
    }
    *func = 0.0;
    for (a = 0; a < VIM; a++) {
      *func += attached * fv->snormal[a] * mp->density *
               (vconv[a] * mp->porosity + pmv->liq_darcy_velocity[a]);
    }

    if (af->Assemble_Jacobian) {

      /*
       * J_s_pmv
       */
      var = POR_POROSITY;
      if (pd->v[pg->imtrx][var]) {
        for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++) {
          phi_j = bf[var]->phi[j_id];
          for (p = 0; p < dim; p++) {
            d_func[0][POR_POROSITY][j_id] +=
                attached * fv->snormal[p] * mp->density *
                (mp->porosity * d_vconv->v[p][i_pore][j_id] + /*isn't this zero? */
                 mp->d_porosity[POR_POROSITY] * vconv[p] * phi_j +
                 pmv->liq_darcy_velocity[p] * mp->d_permeability[POR_POROSITY] * phi_j /
                     mp->permeability);
          }
        }
      }
      var = POR_LIQ_PRES;
      for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++) {
        if (mp->PermeabilityModel != K_TENSOR && mp->PermeabilityModel != ORTHOTROPIC &&
            mp->PermeabilityModel != SM_TENSOR && mp->PermeabilityModel != KC_TENSOR) {
          for (a = 0; a < VIM; a++) {
            d_func[0][var][j_id] -= attached * fv->snormal[a] * mp->density * factor *
                                    mp->permeability * bf[var]->grad_phi[j_id][a];
            d_func[0][var][j_id] -= attached * fv->snormal[a] * mp->density * d_factor[0] *
                                    bf[var]->phi[j_id] * mp->permeability * fv->grad_p_liq[a];
          }
        } else {
          for (a = 0; a < VIM; a++) {
            for (b = 0; b < VIM; b++) {
              d_func[0][var][j_id] -= attached * fv->snormal[a] * mp->density * factor *
                                      mp->perm_tensor[a][b] * bf[var]->grad_phi[j_id][b];
              d_func[0][var][j_id] -= attached * fv->snormal[a] * mp->density * d_factor[0] *
                                      bf[var]->phi[j_id] * mp->perm_tensor[a][b] *
                                      fv->grad_p_liq[b];
            }
          }
        }
      }

      var = POR_POROSITY;
      for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++) {
        if (mp->PermeabilityModel != K_TENSOR && mp->PermeabilityModel != ORTHOTROPIC &&
            mp->PermeabilityModel != SM_TENSOR && mp->PermeabilityModel != KC_TENSOR) {
          for (a = 0; a < VIM; a++) {
            d_func[0][var][j_id] -= attached * fv->snormal[a] * mp->density * factor *
                                    mp->d_permeability[POR_POROSITY] * fv->grad_p_liq[a];
          }
        } else {
          for (a = 0; a < VIM; a++) {
            for (b = 0; b < VIM; b++) {
              d_func[0][var][j_id] -= attached * fv->snormal[a] * mp->density * factor *
                                      mp->d_perm_tensor[a][b][POR_POROSITY] * fv->grad_p_liq[b];
            }
          }
        }
      }

      /*
       * J_s_d
       */

      for (a = 0; a < dim; a++) {
        var = MESH_DISPLACEMENT1 + a;
        if (pd->v[pg->imtrx][var]) {
          for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++) {

            for (p = 0; p < dim; p++) {
              d_func[0][var][j_id] += attached * fv->dsnormal_dx[p][a][j_id] * mp->density *
                                      (vconv[p] * mp->porosity + pmv->liq_darcy_velocity[p]);
              d_func[0][var][j_id] +=
                  attached * fv->snormal[p] * mp->density * (mp->porosity * d_vconv->X[p][a][j_id]);
            }
            if (mp->PermeabilityModel != K_TENSOR && mp->PermeabilityModel != ORTHOTROPIC &&
                mp->PermeabilityModel != SM_TENSOR && mp->PermeabilityModel != KC_TENSOR) {
              for (b = 0; b < VIM; b++) {
                d_func[0][var][j_id] -= attached * fv->snormal[b] * mp->density * factor *
                                        mp->permeability * fv->d_grad_p_liq_dmesh[b][a][j_id];
              }
            } else {
              for (b = 0; b < VIM; b++) {
                for (p = 0; p < VIM; p++) {
                  d_func[0][var][j_id] -= attached * fv->snormal[b] * mp->density * factor *
                                          mp->perm_tensor[b][p] *
                                          fv->d_grad_p_liq_dmesh[p][a][j_id];
                  if (mp->PermeabilityModel == ORTHOTROPIC || mp->PermeabilityModel == KC_TENSOR ||
                      mp->PermeabilityModel == SM_TENSOR) {
                    d_func[0][var][j_id] -= attached * fv->snormal[b] * mp->density * factor *
                                            mp->d_perm_tensor_dx[b][p][a][j_id] * fv->grad_p_liq[p];
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  /***************************** FLUID SIDE *******************************/
  /*
   *  If current material is the fluid phase, calculate normal fluid velocity
   *  in liquid at the interface
   */
  else if (Current_EB_ptr->Elem_Blk_Id == eb_mat_fluid) {

    double rho_liq, d_rho_dF[MDE];
    DENSITY_DEPENDENCE_STRUCT d_rho;

    memset(d_rho_dF, 0, MDE * sizeof(double));

    rho_liq = density(&d_rho, time_value);

    /* first establish nearness */
    /* length_scale = dim == 2 ? 2.0*fv->sdet : 2.0*pow( fv->sdet, 0.5) ; */
    if (ls != NULL && vnormal != 0.0) {

      load_lsi(length_scale);

      if (lsi->near)
        near = 1.0;
    } else {
      near = 0.0;
    }

    for (a = 0; a < VIM; a++)
      x_dot[a] = 0.0;
    if (pd->TimeIntegration != STEADY) {
      if (pd->v[pg->imtrx][MESH_DISPLACEMENT1]) {
        for (a = 0; a < VIM; a++) {
          x_dot[a] = fv_dot->x[a];
        }
      }
    }

    /* Calculate the residual contribution from the velocity */
    *func = 0;
    for (a = 0; a < VIM; a++) {
      *func += fv->snormal[a] * rho_liq * (fv->v[a] - x_dot[a]);
    }

    if (ls != NULL && near == 1.0) {
      *func += vnormal * rho_liq;
    }

    if (af->Assemble_Jacobian) {

      for (jvar = 0; jvar < dim; jvar++) {
        var = VELOCITY1 + jvar;
        if (pd->v[pg->imtrx][var]) {
          for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++) {
            phi_j = bf[var]->phi[j_id];
            d_func[0][var][j_id] += fv->snormal[jvar] * rho_liq * phi_j;
          }
        }
      }

      var = LS;
      if (pd->v[pg->imtrx][var]) {
        for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++) {
          for (a = 0; a < dim; a++) {
            d_func[0][var][j_id] += fv->snormal[a] * d_rho.F[j_id] * (fv->v[a] - x_dot[a]);
          }
        }
      }

      for (jvar = 0; jvar < dim; jvar++) {
        var = MESH_DISPLACEMENT1 + jvar;
        if (pd->v[pg->imtrx][var]) {
          for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++) {
            for (a = 0; a < dim; a++) {
              d_func[0][var][j_id] +=
                  fv->dsnormal_dx[a][jvar][j_id] * rho_liq * (fv->v[a] - x_dot[a]);
            }
          }
        }
      }

      if (pd->TimeIntegration != STEADY) {
        for (jvar = 0; jvar < dim; jvar++) {
          var = MESH_DISPLACEMENT1 + jvar;
          if (pd->v[pg->imtrx][var]) {
            for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++) {
              phi_j = bf[var]->phi[j_id];
              for (a = 0; a < dim; a++) {
                d_func[0][var][j_id] -= fv->snormal[a] * rho_liq * (1.0 + 2.0 * tt) * (phi_j / dt) *
                                        (double)delta(a, jvar);
              }
            }
          }
        }
      }
    }
  } else {
    GOMA_EH(GOMA_ERROR, "No Slip called with incorrect block id");
  }
  return;
} /* end  sat_darcy_continuous_bc  */
/******************************************************************************/

void por_liq_flux_fill(double *func,
                       double d_func[MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                       double tt,
                       double dt,
                       const int eb_mat_solid,
                       const int eb_mat_fluid,
                       double perm,
                       double pc1,
                       double tau,
                       double lsls,
                       int couple_pressure) {
  double H;
  int var, j;

  /***************************** SOLID SIDE *******************************/
  /*
   *  If current material is the solid phase, set pliq
   */
  if (Current_EB_ptr->Elem_Blk_Id == eb_mat_solid) {

    if (mp->PorousMediaType == CONTINUOUS)
      GOMA_EH(GOMA_ERROR, "POR_LIQ_FLUX_FILL requires a porous phase");

    /* Check if we're in the mushy zone. */
    /* Add contributions from level set side of boundary to flux */

    load_lsi_adjmatr(lsls);

    H = lsi->H;

    *func = -perm * (1.0 - H) * (fv->p_liq - pc1);
    *func -= tau * fv_dot->p_liq;

    var = POR_LIQ_PRES;
    if (pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_func[0][var][j] = -perm * (1.0 - H) * bf[var]->phi[j];
        d_func[0][var][j] -= tau * (1.0 + 2.0 * tt) * (bf[var]->phi[j] / dt);
      }
    }

    var = LS;
    for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
      d_func[0][var][j] = perm * lsi->d_H_dF[j] * (fv->p_liq - pc1);
    }

  }

  /***************************** FLUID SIDE *******************************/
  /*
   *  If current material is the fluid phase, calculate level set and set
   *   p_liq value at the interface
   */
  else if (Current_EB_ptr->Elem_Blk_Id == eb_mat_fluid) {

    if (FALSE) {
      load_lsi(lsls);

      H = lsi->H;

      *func = -perm * (1.0 - H) * (-fv->P);

      var = LS;
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_func[0][var][j] = perm * lsi->d_H_dF[j] * (fv->P);
      }

      var = PRESSURE;
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_func[0][var][j] = -perm * (1.0 - H) * (bf[var]->phi[j]);
      }
    }

  } else {
    GOMA_EH(GOMA_ERROR, "POROUS_LIQ_PRESSURE_FILL called with incorrect block id");
  }

  return;
}

void porous_liq_fill(double *func,
                     double d_func[],
                     const int eb_mat_solid,
                     const int eb_mat_fluid,
                     double pc1,
                     double pc2,
                     double lsls,
                     MATRL_PROP_STRUCT *mp2)

{
  double H;
  int var;

  /***************************** SOLID SIDE *******************************/
  /*
   *  If current material is the solid phase, set pliq
   */
  if (Current_EB_ptr->Elem_Blk_Id == eb_mat_solid) {

    if (mp->PorousMediaType == CONTINUOUS)
      GOMA_EH(GOMA_ERROR, "POROUS_LIQ_PRESSURE requires a porous phase");

    /* Check if we're in the mushy zone. */
    /* Add contributions from level set side of boundary to flux */

    load_lsi_adjmatr(lsls);

    *func += (1.0 - lsi->H) * fv->p_liq;

    if (af->Assemble_Jacobian) {
      var = POR_LIQ_PRES;
      if (pd->v[pg->imtrx][var]) {
        d_func[var] += (1.0 - lsi->H);
      }

      var = FILL;

      if (upd->ep[pg->imtrx][var]) {
        d_func[var] += -fv->p_liq * lsi->dH;
      }
    }
  }

  /***************************** FLUID SIDE *******************************/
  /*
   *  If current material is the fluid phase, calculate level set and set
   *   p_liq value at the interface
   */
  else if (Current_EB_ptr->Elem_Blk_Id == eb_mat_fluid) {

    /* Fetch the level set interfacial functions. */
    load_lsi(lsls);
    load_lsi_derivs();
    /**func = (1.0-lsi->H)*(fv->p_liq - (pc1 + (pc2 - pc1) * lsi->H));*/
    H = lsi->H;

    *func += -(1.0 - H) * (pc1 + (pc2 - pc1) * H);

    var = FILL;
    if (pd->v[pg->imtrx][var]) {
      d_func[var] = -(1.0 - H) * ((pc2 - pc1) * lsi->dH) + lsi->dH * (pc1 + (pc2 - pc1) * H);
    }
  } else {
    GOMA_EH(GOMA_ERROR, "POROUS_LIQ_PRESSURE_FILL called with incorrect block id");
  }
  return;
} /* end  porous_liq_pressure_fill_bc  */
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

double interpolate_table_sat(struct Data_Table *table, double x[DIM])

{
  int i, N, Np1, iinter, istartx, istarty, iad;
  double func = 0, y1, y2, y3, y4, tt, uu;
  double *t, *t2, *f;

  N = table->tablelength - 1;
  Np1 = table->tablelength;
  t = table->t;
  t2 = table->t2;
  f = table->f;
  uu = 0.0;

  switch (table->interp_method) { /* switch */
  case LINEAR:                    /* This is knucklehead linear interpolation scheme */
    /* check if absorb or desorb */
    if (x[1] > 0) {
      iad = 1;
    } else {
      iad = 0;
    }

    if (x[0] < t[0]) {
      table->slope[0] = (f[1 + Np1 * iad] - f[0 + Np1 * iad]) / (t[1] - t[0]);
      func = f[0 + Np1 * iad] + (table->slope[0]) * (x[0] - t[0]);
    }

    for (i = 0; x[0] >= t[i] && i < N; i++) {
      if (x[0] >= t[i] && x[0] < t[i + 1]) {
        table->slope[0] = (f[i + 1 + Np1 * iad] - f[i + Np1 * iad]) / (t[i + 1] - t[i]);
        func = f[i + Np1 * iad] + (table->slope[0]) * (x[0] - t[i]);
      }
    }
    if (x[0] >= t[N]) {
      table->slope[0] = (f[N + Np1 * iad] - f[N - 1 + Np1 * iad]) / (t[N] - t[N - 1]);
      func = f[N + Np1 * iad] + (table->slope[0]) * (x[0] - t[N]);
    }
    table->slope[1] = 0.0;
    break;
  case BILINEAR: /* BILINEAR Interpolation Scheme */
    /* check if absorb or desorb */
    if (x[2] > 0) {
      iad = 1;
    } else {
      iad = 0;
    }

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

    for (i = 0; t[i] <= x[0] && i < N + 1; i = i + iinter) {
      istartx = i;
    }

    istarty = istartx;

    for (i = istartx + 1; t2[i] <= x[1] && i < istartx + iinter - 1; i++) {
      istarty = i;
    }

    y1 = f[istarty + Np1 * iad];
    y2 = f[istarty + 1 + Np1 * iad];
    y3 = f[istarty + 1 + iinter + Np1 * iad];
    y4 = f[istarty + iinter + Np1 * iad];

    tt = (x[1] - t2[istarty]) / (t2[istarty + 1] - t2[istarty]);
    uu = (x[0] - t[istarty]) / (t[istarty + 1 + iinter] - t[istarty]);

    func = (1. - tt) * (1. - uu) * y1 + tt * (1. - uu) * y2 + tt * uu * y3 + (1. - tt) * uu * y4;

    table->slope[1] =
        (1.0 - uu) * ((f[istarty + 1 + Np1 * iad] - f[istarty + Np1 * iad]) /
                      (t2[istarty + 1] - t2[istarty])) +
        uu * ((f[istarty + 1 + iinter + Np1 * iad] - f[istarty + iinter + Np1 * iad]) /
              (t2[istarty + 1 + iinter] - t2[istarty + iinter]));
    /* slope0 not needed */
    table->slope[0] = (f[istarty + 1 + Np1 * iad] - f[istarty + 1 + iinter + Np1 * iad]) /
                      (t[istarty + 1] - t[istarty + 1 + iinter]);
    break;
  default:
    GOMA_EH(GOMA_ERROR, "Table interpolation order for Saturation not implemented");
  }

  return (func);
}
/******************************************************************************/
/******************************************************************************/
int evaluate_sat_hyst_criterion(int ip, int ielem, struct Porous_Media_Terms *pmt, dbl tt, dbl dt) {
  const int i_pl = 0;
  dbl liq_inv_dot = 0.0, cap_press = 0.0;

  if (pmt != NULL) {
    pmt->Inventory_solvent[i_pl] = pmv->bulk_density[i_pl];
  }

  if (pd->TimeIntegration != STEADY) {
    /*
     * Liquid solvent component first, e.g. water
     */
    /*  pmt->Inventory_solvent_old[i_pl] = pmv_old->bulk_density[i_pl];

    pmt->Inventory_solvent_dot[i_pl] =
        (1.0 + 2.0 * tt) *
        (pmt->Inventory_solvent[i_pl] - pmt->Inventory_solvent_old[i_pl])/dt
        - 2.0 * tt  * pmt->Inventory_solvent_dot_old[i_pl];

    liq_inv_dot =  pmt->Inventory_solvent_dot[i_pl];
    for (a = 0; a < DIM; a++)
      {
        liq_inv_dot += pmt->conv_flux[0][a];
      }
    */

    /* OK, we are going to equate liq_inv_dot to fv_dot->p_liq for
     * sport, to see if it works better.   Here if d(P_liq)/dt is > 0
     * means that p_cap is decreasing with time and the medium is saturating
     * which meants we should switch to the imbibing curve if we aren't there
     * already.   Vice Versa if it is less than zero
     */

    if (pd->v[pg->imtrx][POR_LIQ_PRES]) {
      liq_inv_dot = fv_dot->p_liq;
      cap_press = fv->p_gas - fv->p_liq;
    } else if (pd->v[pg->imtrx][SHELL_PRESS_OPEN]) {
      liq_inv_dot = fv_dot->sh_p_open;
      cap_press = -fv->sh_p_open;
      // PRS NOTE: need to add sh_p_open_2 for second layer here
    }

    /*If the accumulation is positive, above a certain
     *tolerance level, and was positive previously, then we will
     *remain on the same wetting curve with same previous curve
     *parameters.   Likewise for negative accumulation with the
     *drying curve.    We will switch if the accumulation rate
     *changes sign and the magnitude is above a certan tolerance
     *level.
     */

    if (liq_inv_dot > 0.0 &&
        Element_Blocks[ei[pg->imtrx]->elem_blk_index].ElemStorage[ielem].sat_curve_type[ip] ==
            0.0) {
      /* We were on a wetting curve, and will remain so */
    } else if (liq_inv_dot <= 0.0 && Element_Blocks[ei[pg->imtrx]->elem_blk_index]
                                             .ElemStorage[ielem]
                                             .sat_curve_type[ip] == 1.0) {
      /* We were on a drying/draining curve, and will remain so */
    } else if (liq_inv_dot > 0.0 && Element_Blocks[ei[pg->imtrx]->elem_blk_index]
                                            .ElemStorage[ielem]
                                            .sat_curve_type[ip] == 1.0) {
      /* We were on a drying/draining curve but now may potentially switch to a wetting curve */

      if ((pd->v[pg->imtrx][POR_LIQ_PRES]) || (pd->v[pg->imtrx][SHELL_PRESS_OPEN])) {
        if (fabs(liq_inv_dot) > mp->u_saturation[9] && mp->saturation <= 0.9999) {
          Element_Blocks[ei[pg->imtrx]->elem_blk_index].ElemStorage[ielem].sat_curve_type[ip] = 0.0;
          Element_Blocks[ei[pg->imtrx]->elem_blk_index].ElemStorage[ielem].Sat_QP_tn[ip] =
              mp->saturation;
          Element_Blocks[ei[pg->imtrx]->elem_blk_index].ElemStorage[ielem].p_cap_QP[ip] = cap_press;
        }
      }
    } else if (liq_inv_dot <= 0.0 && Element_Blocks[ei[pg->imtrx]->elem_blk_index]
                                             .ElemStorage[ielem]
                                             .sat_curve_type[ip] == 0.0) {
      /* We were on a wetting curve but now may potentially switch to a drying curve */

      if ((pd->v[pg->imtrx][POR_LIQ_PRES]) || (pd->v[pg->imtrx][SHELL_PRESS_OPEN])) {

        if (fabs(liq_inv_dot) > mp->u_saturation[9]) {
          Element_Blocks[ei[pg->imtrx]->elem_blk_index].ElemStorage[ielem].sat_curve_type[ip] = 1.0;
          Element_Blocks[ei[pg->imtrx]->elem_blk_index].ElemStorage[ielem].Sat_QP_tn[ip] =
              mp->saturation;
          Element_Blocks[ei[pg->imtrx]->elem_blk_index].ElemStorage[ielem].p_cap_QP[ip] = cap_press;
        }
      }
    }
  }

  return (1);
}

/*****************************************************************************/
double por_mass_source_model(double d_MassSource[MAX_VARIABLE_TYPES + MAX_CONC][MDE])
/******************************************************************************
 *
 *  A function which computes mass source term in porous media transport and their
 *  Jacobian sensitivities
 *
 *  Kristianto Tjiptowidjojo (June 2015)
 *
 *
 ******************************************************************************/
{
  int b, j, var;
  double phi_j;
  double MassSource = 0.0;

  double sat_min;
  double width, alpha;
  double sat_center, sat_normalized;
  double Hside = 0.0;
  double dHside_dS = 0.0;
  double tau = 0.0;
  double sink_mass_max = 0.0, nexp = 0.0;
  double sink_mass_clip = fv->sink_mass;

  double saturation = 0.0;

  int dim = pd->Num_Dim;

  /* Get saturation based on the formulation used*/

  if (pd->v[pg->imtrx][SHELL_SAT_1]) /* Saturation formulation */
  {
    saturation = fv->sh_sat_1;
  } else if ((pd->v[pg->imtrx][SHELL_PRESS_OPEN]) || (pd->v[pg->imtrx][POR_LIQ_PRES])) {
    saturation = mp->saturation;
  }

  /*** Calculate MassSource based on the selected models ***/

  switch (mp->PorousSinkConstantsModel) {

  case LINEAR:

    if (saturation >= mp->u_porous_sink_constants[6]) {
      tau = mp->u_porous_sink_constants[0];
    } else {
      tau = 0.;
    }

    sink_mass_max = mp->u_porous_sink_constants[1];

    MassSource = -tau * mp->u_porous_sink_constants[2] * (sink_mass_max - fv->sink_mass) *
                 saturation / sink_mass_max / fv->volume_change;

    /* Again, I don't think this term belongs, contrary to the paper from which
       it came.   Disappears with application of the Reynolds Transport theorem */

    /*
    for(b=0; b < VIM; b++)
       {
        MassSource += fv_dot->x[b] * mp->density *
                      (mp->porosity * mp->d_saturation[POR_LIQ_PRES]*fv->grad_p_liq[b] +
                       mp->saturation * fv->grad_porosity[b]);
       }
    */

    break;

  case POWER_LAW:

    /* Evaluate heaviside function based on minimum saturation */
    sat_min = mp->u_porous_sink_constants[3];
    width = mp->u_porous_sink_constants[4];
    alpha = 0.5 * width;
    sat_center = sat_min - alpha;
    sat_normalized = saturation - sat_center;

    if (saturation >= sat_min) {
      Hside = 1.0;
      dHside_dS = 0.0;
    } else if (saturation <= (sat_min - width)) {
      Hside = 0.0;
      dHside_dS = 0.0;
    } else {
      Hside = 0.5 * (1. + sat_normalized / alpha + sin(M_PIE * sat_normalized / alpha) / M_PIE);
      dHside_dS = 0.5 * (1.0 / alpha + cos(M_PIE * sat_normalized / alpha) / alpha);
    }

    tau = mp->u_porous_sink_constants[0];
    sink_mass_max = mp->u_porous_sink_constants[1];
    nexp = mp->u_porous_sink_constants[2];

    if (fv->sink_mass < sink_mass_max) {
      sink_mass_clip = fv->sink_mass;
    } else {
      sink_mass_clip = sink_mass_max;
    }

    MassSource = -tau * pow(((sink_mass_max - sink_mass_clip) / sink_mass_max), nexp) * saturation /
                 mp->density;

    MassSource *= Hside;

    break;

  default:

    GOMA_EH(GOMA_ERROR, "No valid porous sink models were found");
  }

  /****** Calculate their Jacobian sensitivities */

  memset(d_MassSource, 0.0, sizeof(double) * (MAX_VARIABLE_TYPES + MAX_CONC) * MDE);

  switch (mp->PorousSinkConstantsModel) {

  case LINEAR:

    var = POR_SINK_MASS;
    if (pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        phi_j = bf[var]->phi[j];

        d_MassSource[var][j] = tau * mp->u_porous_sink_constants[2] * phi_j * saturation /
                               sink_mass_max / fv->volume_change;
      }
    }

    var = POR_LIQ_PRES;
    if (pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        phi_j = bf[var]->phi[j];

        d_MassSource[var][j] = -tau * mp->u_porous_sink_constants[2] *
                               (sink_mass_max - fv->sink_mass) * mp->d_saturation[POR_LIQ_PRES] /
                               sink_mass_max / fv->volume_change * phi_j;

        /*
        for(a=0; a < VIM; a++)
           {
            d_MassSource[var][j] +=  fv_dot->x[a] * mp->density *
                                     (mp->porosity *
        mp->d_saturation[POR_LIQ_PRES]*bf[var]->grad_phi[j][a]);
           }
        */
        /*PRS Caveate: I think an additional sensitivity to of sat(p_liq) needs to be in here */
      }
    }

    var = MESH_DISPLACEMENT1;
    if (pd->v[pg->imtrx][var]) {
      for (b = 0; b < dim; b++) {
        var = MESH_DISPLACEMENT1 + b;

        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_MassSource[var][j] = fv->d_volume_change_dx[b][j] * tau *
                                 mp->u_porous_sink_constants[2] * (sink_mass_max - fv->sink_mass) *
                                 saturation / sink_mass_max / fv->volume_change / fv->volume_change;

          /*see comment on residual equation above about this term not belonging */

          /*
          d_MassSource[var][j] = ((1.+2.*tt) * bf[var]->phi[j]/dt) *
                                 mp->density * (mp->porosity *
          mp->d_saturation[POR_LIQ_PRES]*fv->grad_p_liq[b] + mp->saturation * fv->grad_porosity[b]);

          d_MassSource[var][j] = fv_dot->x[a] * mp->density *
                                 (mp->porosity *
          mp->d_saturation[POR_LIQ_PRES]*fv->d_grad_p_liq_dmesh[b][b][j] + mp->saturation *
          fv->d_grad_porosity_dmesh[b][b][j]);
          */
        }
      }
    }

    break;

  case POWER_LAW:

    var = POR_SINK_MASS;
    if (pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        phi_j = bf[var]->phi[j];
        d_MassSource[var][j] =
            Hside * tau * nexp *
            pow(((sink_mass_max - sink_mass_clip) / sink_mass_max), (nexp - 1.0)) *
            (phi_j / sink_mass_max) * saturation / mp->density;
      }
    }

    var = POR_LIQ_PRES;
    if (pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        phi_j = bf[var]->phi[j];

        d_MassSource[var][j] = -Hside * tau *
                                   pow(((sink_mass_max - sink_mass_clip) / sink_mass_max), nexp) *
                                   mp->d_saturation[var] / mp->density * phi_j -
                               dHside_dS * mp->d_saturation[var] * phi_j * tau *
                                   pow(((sink_mass_max - sink_mass_clip) / sink_mass_max), nexp) *
                                   saturation / mp->density;
      }
    }

    var = SHELL_PRESS_OPEN;
    if (pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        phi_j = bf[var]->phi[j];

        d_MassSource[var][j] = -Hside * tau *
                                   pow(((sink_mass_max - sink_mass_clip) / sink_mass_max), nexp) *
                                   mp->d_saturation[var] / mp->density * phi_j -
                               dHside_dS * mp->d_saturation[var] * phi_j * tau *
                                   pow(((sink_mass_max - sink_mass_clip) / sink_mass_max), nexp) *
                                   saturation / mp->density;
      }
    }

    var = SHELL_SAT_1;
    if (pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        phi_j = bf[var]->phi[j];

        d_MassSource[var][j] = -Hside * tau *
                                   pow(((sink_mass_max - sink_mass_clip) / sink_mass_max), nexp) *
                                   1.0 / mp->density * phi_j -
                               dHside_dS * phi_j * tau *
                                   pow(((sink_mass_max - sink_mass_clip) / sink_mass_max), nexp) *
                                   saturation / mp->density;
      }
    }

    break;
  }

  return (MassSource);
}
