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

/* mm_fill_pthings.c
 *
 * MMH
 * I cp'ed mm_fill_terms.c and pruned down to assemble_momentum.
 * Modified that for particle momentum.
 */

/* #define DEBUG_MOMENTUM */
/* #define DEBUG */

#define Almost_ONE 0.9999

/* Standard include files */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* GOMA include files */

#include "mm_fill_pthings.h"

#include "density.h"
#include "el_elm.h"
#include "mm_as.h"
#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "mm_eh.h"
#include "mm_fill_rs.h"
#include "mm_fill_solid.h"
#include "mm_fill_terms.h"
#include "mm_mp.h"
#include "mm_mp_const.h"
#include "mm_mp_structs.h"
#include "mm_qtensor_model.h"
#include "mm_viscosity.h"
#include "rf_fem.h"
#include "rf_fem_const.h"
#include "rf_solver.h"
#include "std.h"

#define GOMA_MM_FILL_PTHINGS_C

/*
 *********************************************************************
 * A few static function definitions.
 */

/* Assuming x != 1.0 */

static double N_func(double x) {
  double omx;
  if (x >= 1.0)
    x = Almost_ONE;
  omx = 1.0 - x;
  omx *= omx;
  omx *= (1.0 - x);
  return ((1.0 + x + x * x - x * x * x) / omx);
}

static double d_N_func_d_phi(double x) {
  double omx;
  if (x >= 1.0)
    x = Almost_ONE;
  omx = 1.0 - x;
  omx *= omx;
  omx *= omx;
  return ((4.0 + 4.0 * x - 2.0 * x * x) / omx);
}

static double Enskog(double x) {
  /*
  if(x==0)
    return 1;
  else
    {
      omx=(1.0-x);
      omx*=omx;
      omx*=(1.0-x);
      return( 0.5*(2.0-x)/omx );
    }
    */
  if (x >= 1.0)
    x = Almost_ONE;
  if (x > 0.0)
    return ((N_func(x) - Almost_ONE) / (4.0 * x));
  else
    return 1.0;
}

static double d_Enskog_d_phi(double x) {
  if (x >= 1.0)
    x = Almost_ONE;
  if (x > 0.0)
    return ((x * d_N_func_d_phi(x) - N_func(x) + 1) / (4.0 * x * x));
  else
    return (20.0 / 8.0);
}
/**************************************************************************/

/* assemble_pmomentum -- assemble terms (Residual &| Jacobian) for momentum eqns
 *
 * MMH
 * Again, I have just cp'ed this from mm_fill_terms.c and made
 * some modifications.  In particular, I have left lots of things in
 * that are currently not used.  If we want to put them in later,
 * they are already set up to go.  There are variables that are accessed
 * here that do not belong here.  For example, fv->S corresponds to the
 * fluid velocities, not the particle velocities...
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
 */

int assemble_pmomentum(dbl time_value, /* current time for density model */
                       dbl tt,         /* parameter to vary time integration from
                                          explicit (tt = 1) to implicit (tt = 0) */
                       dbl dt)         /* current time step size */
{
  int dim, p, q, a, b, eqn, var, ii, peqn, pvar, w, ledof;
  int i, j, m, status;
  struct Basis_Functions *bfm;

  dbl pv[DIM];     /* Velocity field. */
  dbl pv_dot[DIM]; /* time derivative of velocity field. */
  dbl x_dot[DIM];  /* current position field derivative wrt time. */

  dbl h3;          /* Volume element (scale factors). */
  dbl dh3dmesh_bj; /* Sensitivity to (b,j) mesh dof. */

  dbl grad_pv[DIM][DIM]; /* Gradient of pv. */
  dbl gamma[DIM][DIM];   /* shrearrate tensor based on velocity */

  dbl Pi[DIM][DIM];     /* Total stress tensor (multiplied by coeff). */
  dbl Pi_raw[DIM][DIM]; /* Rate of strain tensor (no coeff). */

  dbl rho; /* Density. */

  dbl coeff;                    /* Takes place of viscosity...  Coefficient
                                 * of the particle stress tensor.
                                 */
  dbl f[DIM];                   /* Body force. */
  dbl dfdT[DIM][MDE];           /* For temperature dependence. */
  dbl dfdX[DIM][DIM][MDE];      /* For spatial dependence. */
  dbl dfdv[DIM][DIM][MDE];      /* For velocity dependence. */
  dbl dfdC[DIM][MAX_CONC][MDE]; /* For concentration dependence. */

  dbl det_J; /* determinant of element Jacobian */

  dbl d_det_J_dmesh_bj; /* for specific (b,j) mesh dof */

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
  dbl(*grad_phi_e)[DIM][DIM][DIM] = NULL;
  dbl(*grad_phi_i_e_a)[DIM];
  /*
   * Interpolation functions for variables and some of their derivatives.
   */

  dbl phi_j;
  dbl(*d_grad_phi_i_e_a_dmesh)[DIM][DIM][MDE];

  /* MMH
   * For the fluid phase (assemble_momentum), Pi = mu * gamma.
   * For the particle phase, Pi = coeff * Pi_raw.
   */
  dbl d_gamma_dvbj[DIM][DIM][DIM][MDE];
  dbl d_gamma_doubledot_dvbj[DIM][MDE];
  dbl d_Pi_pv[DIM][DIM][DIM][MDE];
  dbl d_Pi_mesh[DIM][DIM][DIM][MDE];
  dbl d_Pi_C[DIM][DIM][MAX_CONC][MDE];
  dbl d_Pi_T[DIM][DIM][MDE];
  dbl hold1;

  dbl wt;

  /*
   * Variables for vicosity and derivative
   */

  dbl d_coeff_d_pv[DIM][MDE];
  dbl d_coeff_d_mesh[DIM][MDE];
  dbl d_coeff_d_C[MAX_CONC][MDE];

  /* coefficient variables for the Brinkman Equation
     in flows through porous media: KSC on 5/10/95 */
  dbl por;   /* porosity of porous media */
  dbl por2;  /* square of porosity */
  dbl per;   /* permeability of porous media */
  dbl vis;   /* flowing-liquid viscosity */
  dbl sc;    /* inertial coefficient */
  dbl speed; /* magnitude of the velocity vector */

  /* density derivatives */
  DENSITY_DEPENDENCE_STRUCT d_rho_struct; /* density dependence */
  DENSITY_DEPENDENCE_STRUCT *d_rho = &d_rho_struct;

  /* set porous-flow parameters depending on which zone we are in. KSC on 5/10/95 */
  /* dbl evss_f; */ /* flag to add in the continuous velocity
                  gradient tensor for Fortin's formulation */

  /* MMH
   * Variables used for the particle momentum equations when looking
   * at the particle momentum model.
   */
  int species = -1;         /* Species number of particle phase */
  double rho_p = 1e12;      /* Fluid and solid densities */
  double p_vol_frac = 1e12; /* Particle volume fraction (phi) */
  double mul1;              /* Used for the strain tensor */
  double EpEp[DIM][DIM];    /* For tensor . tensor */
  double Epinv;             /* Ep 2nd invariant, 1/2*Ep:Ep */
  double gammadot;          /* Magnitude of fluid shear. */
  double Ensval;            /* Enskog(p_vol_frac) */
  double Nval;              /* N_func(p_vol_frac) */

  double particle_radius;
  double tmp;

  const double mul2 = 1;
  const double mul3 = M_PIE * M_PIE / 4 - 1;
  const double mul4 = 0.5;

  memset(d_Pi_T, 0, sizeof(double) * DIM * DIM * MDE);

  status = 0;

  /*
   * Unpack variables from structures for local convenience...
   */

  eqn = R_PMOMENTUM1;

  /*
   * Bail out fast if there's nothing to do...
   */

  if (!pd->e[pg->imtrx][eqn]) {
    return (status);
  }

  dim = pd->Num_Dim;

  /* MMH
   * I don't need any of this, but might want to put in something analogous
   * later.
   */
  /*
  pv_s[0][0] = POLYMER_STRESS11;
  pv_s[0][1] = POLYMER_STRESS12;
  pv_s[0][2] = POLYMER_STRESS13;
  pv_s[1][0] = POLYMER_STRESS12;
  pv_s[1][1] = POLYMER_STRESS22;
  pv_s[1][2] = POLYMER_STRESS23;
  pv_s[2][0] = POLYMER_STRESS13;
  pv_s[2][1] = POLYMER_STRESS23;
  pv_s[2][2] = POLYMER_STRESS33;

  pv_g[0][0] = VELOCITY_GRADIENT11;
  pv_g[0][1] = VELOCITY_GRADIENT12;
  pv_g[1][0] = VELOCITY_GRADIENT21;
  pv_g[1][1] = VELOCITY_GRADIENT22;
  pv_g[0][2] = VELOCITY_GRADIENT13;
  pv_g[1][2] = VELOCITY_GRADIENT23;
  pv_g[2][0] = VELOCITY_GRADIENT31;
  pv_g[2][1] = VELOCITY_GRADIENT32;
  pv_g[2][2] = VELOCITY_GRADIENT33;
  */

  wt = fv->wt;

  det_J = bf[eqn]->detJ; /* Really, ought to be mesh eqn. */

  h3 = fv->h3; /* Differential volume element (scales). */

  /*
   * Material property constants, etc. Any variations for this
   * Gauss point were evaluated in load_material_properties().
   */

  if (pd->e[pg->imtrx][R_PMOMENTUM1]) {
    /* This is the species number of the particle phase. */
    species = (int)mp->u_density[0];
    rho_p = mp->u_density[2];
    p_vol_frac = fv->c[species];
    if (p_vol_frac < 0.0 || p_vol_frac > 1.0) {
      printf("assemble_pmomentum: p_vol_frac=%g, exiting\n", p_vol_frac);
      exit(0);
    }
  }

  /*** Density ***/
  /*  rho  = mp->density; */
  /* This sets the d_rho_->* stuff */
  rho = density(d_rho, time_value);
  /* MMH
   * We have to fix rho b/c it defaults to the fluid density.
   * Currently there aren't any fancy dependencies...
   */
  rho = rho_p;

  if (pd->e[pg->imtrx][eqn] & T_POROUS_BRINK) {
    if (mp->PorousMediaType != POROUS_BRINKMAN)
      GOMA_WH(-1, "Set Porous term multiplier in continuous medium");
    /* Short-hand notation for the four parameters in the Brinkman Equation. */
    por = mp->porosity;
    per = mp->permeability;
    vis = mp->FlowingLiquid_viscosity;
    sc = mp->Inertia_coefficient;

  } else {
    por = 1.;
    per = 1.;
    vis = mp->viscosity;
    sc = 0.;
  }

  pmomentum_source_term(f, dfdT, dfdX, dfdC, dfdv);

  eqn = R_PMOMENTUM1;

  /*
   * Field variables...
   */
  for (a = 0; a < WIM; a++) {
    pv[a] = fv->pv[a];
    if (pd->TimeIntegration != STEADY && pd->v[pg->imtrx][MESH_DISPLACEMENT1 + a]) {
      x_dot[a] = fv_dot->x[a];
    } else {
      x_dot[a] = 0.;
    }
    if (pd->TimeIntegration != STEADY) {
      pv_dot[a] = fv_dot->pv[a];
    } else {
      pv_dot[a] = 0.;
    }
  }

  /* for porous media stuff */
  speed = 0.0;
  for (a = 0; a < WIM; a++) {
    speed += pv[a] * pv[a];
  }
  speed = sqrt(speed);

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
      grad_pv[a][b] = fv->grad_pv[a][b];
    }
  }

  /*
  if ( pd->v[pg->imtrx][POLYMER_STRESS11] && (vn->evssModel == EVSS_F) )
    {
      evss_f = 1.;
    }
  else
    {
      evss_f = 0.;
    }
  */
  /* load up shearrate tensor based on velocity */
  for (a = 0; a < VIM; a++) {
    for (b = 0; b < VIM; b++) {
      gamma[a][b] = grad_pv[a][b] + grad_pv[b][a];
    }
  }

  /*
   * Stress tensor...(Note "anti-BSL" sign convention on deviatoric stress)
   */

  if (cr->MomentumFluxModel == CR_MF_NEWTON_0) {
    /* *** section commented-out by MMH: clarified by dal for AIX version ***
    temp = viscosity(gn, &mu, gamma, &d_mu_dgd, d_mu_dpv, d_mu_dmesh,
                     d_mu_dT, d_mu_dp, d_mu_dC);

    if ( pd->v[pg->imtrx][POLYMER_STRESS11] )
      {
        for ( mode=0; mode<vn->modes; mode++)
          {
             COMMENT: get polymer viscosity
            temp = viscosity(ve[mode]->gn, &mup, gamma, &d_mup_dgd, d_mup_dpv,
                             d_mup_dmesh, d_mup_dT, d_mup_dp, d_mup_dC);

            var = PVELOCITY1;
            if (pd->v[pg->imtrx][var] )
              {
                for ( a=0; a<WIM; a++)
                  {
                    for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
                      {
                        d_mu_dpv[a][j] =  d_mup_dpv[a][j] +  d_mu_dpv[a][j];
                      }
                  }
              }

        var = TEMPERATURE;
        if (pd->v[pg->imtrx][var] )
          {
            for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
              {
                d_mu_dT[j] = d_mu_dT[j] + d_mup_dT[j];
              }
          }
        var = PRESSURE;
        if (pd->v[pg->imtrx][var] )
          {
            for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
              {
                d_mu_dp[j] = d_mu_dp[j] + d_mup_dp[j];
              }
          }

            var = MESH_DISPLACEMENT1;
            if (pd->v[pg->imtrx][var] )
              {
                for ( a=0; a<dim; a++)
                  {
                    for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
                      {
                        d_mu_dmesh[a][j] = d_mup_dmesh[a][j] + d_mu_dmesh[a][j];
                      }
                  }
              }

            var = TEMPERATURE;
            if (pd->v[pg->imtrx][var] )
              {
                for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
                  {
                    d_mu_dT[j] = d_mu_dT[j] + d_mup_dT[j];
                  }
              }
            var = PRESSURE;
            if (pd->v[pg->imtrx][var] )
              {
                for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
                  {
                    d_mu_dp[j] = d_mu_dp[j] + d_mup_dp[j];
                  }
              }

            var = MASS_FRACTION;
            if (pd->v[pg->imtrx][var] )
              {
                for ( w=0; w<pd->Num_Spec; w++)
                  {
                    for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
                      {
                        d_mu_dC[w][j] = d_mup_dC[w][j]+d_mu_dC[w][j];
                      }
                  }
              }

            mu = mu + mup;
          }
      }
      **** end of MMH commented-out section *** */

    /* MMH
     * I don't pick up any of the d_mu_d_* things b/c there is no call
     * to viscosity here.  I just kept it all in this function, instead
     * of making up a new viscosity_* call.  There is no real "viscosity",
     * as in the fluid phase, but there is a coefficient of a tensor,
     * so we can pretend that it is mu...
     */

    memset(d_coeff_d_pv, 0, sizeof(double) * DIM * MDE);
    memset(d_coeff_d_mesh, 0, sizeof(double) * DIM * MDE);
    memset(d_coeff_d_C, 0, sizeof(double) * MAX_CONC * MDE);

    Ensval = Enskog(p_vol_frac);
    Nval = N_func(p_vol_frac);

    gammadot = 0.0;
    for (a = 0; a < VIM; a++)
      for (b = 0; b < VIM; b++) {
        tmp = (fv->grad_v[a][b] + fv->grad_v[b][a]);
        gammadot += tmp * tmp;
      }
    gammadot = sqrt(0.5 * gammadot);

    particle_radius = 0.10; /* in cm's */
    mul1 = p_vol_frac * gammadot * particle_radius * Ensval;
    mul1 *= mul1;
    mul1 *= rho * p_vol_frac * Nval;

    coeff = mul1 * mul2;

    for (j = 0; j < ei[pg->imtrx]->dof[MASS_FRACTION]; j++) {
      d_coeff_d_C[species][j] =
          mul2 * rho * gammadot * gammadot * particle_radius * particle_radius * p_vol_frac *
          p_vol_frac * Ensval *
          (3.0 * Ensval * Nval + 2.0 * p_vol_frac * Nval * d_Enskog_d_phi(p_vol_frac) +
           p_vol_frac * Ensval * d_N_func_d_phi(p_vol_frac));
      d_coeff_d_C[species][j] *= bf[MASS_FRACTION]->phi[j];
    }

    eqn = R_PMOMENTUM1;

    /* MMH
     * Calculating the meat of the YAB model.  This is the
     * stress tensor for the particle momentum equations.
     */

    Epinv = 0.0;
    for (a = 0; a < VIM; a++)
      for (b = 0; b < VIM; b++) {
        Epinv += gamma[a][b] * gamma[b][a];
        EpEp[a][b] = 0.0;
        for (p = 0; p < VIM; p++)
          EpEp[a][b] += gamma[a][p] * gamma[p][b];
      }

    for (a = 0; a < VIM; a++)
      for (b = 0; b < VIM; b++) {
        Pi_raw[a][b] = mul3 * EpEp[a][b] + mul4 * Epinv * (double)delta(a, b);
        Pi[a][b] = coeff * Pi_raw[a][b];
      }

    grad_phi_e = bf[eqn]->grad_phi_e;

    /*
     * d_Pi_d_pv
     */
    for (p = 0; p < VIM; p++)
      for (q = 0; q < VIM; q++)
        for (b = 0; b < WIM; b++)
          for (j = 0; j < ei[pg->imtrx]->dof[PVELOCITY1]; j++)
            d_gamma_dvbj[p][q][b][j] = grad_phi_e[j][b][p][q] + grad_phi_e[j][b][q][p];
    for (b = 0; b < WIM; b++)
      for (j = 0; j < ei[pg->imtrx]->dof[PVELOCITY1]; j++) {
        d_gamma_doubledot_dvbj[b][j] = 0.0;
        for (p = 0; p < VIM; p++)
          for (q = 0; q < VIM; q++)
            d_gamma_doubledot_dvbj[b][j] +=
                gamma[p][q] * d_gamma_dvbj[q][p][b][j] + d_gamma_dvbj[p][q][b][j] * gamma[q][p];
      }

    for (p = 0; p < VIM; p++) {
      for (q = 0; q < VIM; q++) {
        for (b = 0; b < WIM; b++) {
          for (j = 0; j < ei[pg->imtrx]->dof[PVELOCITY1]; j++) {
            /* MMH Ugh! */
            hold1 = 0.0;
            for (m = 0; m < VIM; m++)
              hold1 +=
                  (gamma[p][m] * d_gamma_dvbj[m][q][b][j] + d_gamma_dvbj[p][m][b][j] * gamma[m][q]);

            d_Pi_pv[p][q][b][j] =
                mul3 * hold1 + mul4 * d_gamma_doubledot_dvbj[b][j] * (double)delta(p, q);

            d_Pi_pv[p][q][b][j] *= coeff;
            d_Pi_pv[p][q][b][j] += d_coeff_d_pv[b][j] * Pi_raw[p][q];
          }
        }
      }
    }

    /*
     * MMH:
     * This looks like it would have to be changed just like d_Pi_dpv
     * above.  Something to remember if moving meshes are included.
     *
     * d_Pi_d_mesh
     */
    for (p = 0; p < VIM; p++) {
      for (q = 0; q < VIM; q++) {
        for (b = 0; b < dim; b++) {
          for (j = 0; j < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; j++) {
            d_Pi_mesh[p][q][b][j] =
                coeff * (fv->d_grad_pv_dmesh[p][q][b][j] + fv->d_grad_pv_dmesh[q][p][b][j]) +
                d_coeff_d_mesh[b][j] * Pi_raw[p][q];
          }
        }
      }
    }

    /*
     * d_Pi_d_C
     */
    for (p = 0; p < VIM; p++) {
      for (q = 0; q < VIM; q++) {
        for (j = 0; j < ei[pg->imtrx]->dof[MASS_FRACTION]; j++) {
          for (w = 0; w < pd->Num_Species_Eqn; w++) {
            d_Pi_C[p][q][w][j] = d_coeff_d_C[w][j] * Pi_raw[p][q];
          }
        }
      }
    }
  } else {
    GOMA_EH(-1, "Unimplemented momentum constitutive relation.");
  }

  /*
   * Residuals_________________________________________________________________
   */

  if (af->Assemble_Residual) {
    /*
     * Assemble each component "a" of the momentum equation...
     */
    for (a = 0; a < WIM; a++) {
      eqn = R_PMOMENTUM1 + a;
      peqn = upd->ep[pg->imtrx][eqn];
      bfm = bf[eqn];

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

          mass = 0.0;
          grad_phi_i_e_a = grad_phi_e[i][a];
          if (pd->TimeIntegration != STEADY) {
            if (pd->e[pg->imtrx][eqn] & T_MASS) {
              mass = pv_dot[a] * rho;
              mass *= p_vol_frac;
              mass *= -phi_i * det_J * wt;
              mass *= h3;
              mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
            }
            /* if porous flow is considered. KSC on 5/10/95 */
            if (pd->e[pg->imtrx][eqn] & T_POROUS_BRINK) {
              mass /= por;
            }
          }

          /* Residual */
          advection = 0.0;
          if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
            for (p = 0; p < WIM; p++) {
              advection += (pv[p] - x_dot[p]) * grad_pv[p][a];
            }

            advection *= rho;

            advection *= p_vol_frac;

            advection *= -phi_i * det_J * wt * h3;

            advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];

            /* if porous flow is considered. KSC on 5/10/95 */
            if (pd->e[pg->imtrx][eqn] & T_POROUS_BRINK) {
              por2 = por * por;
              advection /= por2;
            }
          }

          porous = 0.;
          if (pd->e[pg->imtrx][eqn] & T_POROUS_BRINK) {
            porous = pv[a] * (rho * sc * speed / sqrt(per) + vis / per);
            porous *= -phi_i * det_J * wt;
            porous *= h3;
            porous *= pd->etm[pg->imtrx][eqn][(LOG2_POROUS_BRINK)];
          }

          /* Residual */
          diffusion = 0.0;
          if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {

            for (p = 0; p < VIM; p++) {
              for (q = 0; q < VIM; q++) {
                diffusion += grad_phi_i_e_a[p][q] * Pi[q][p];
              }
            }
            diffusion *= -det_J * wt;
            diffusion *= h3;
            diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
          }

          /*
           * Source term...
           */

          source = 0.0;
          if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
            source += f[a];
            source *= phi_i * det_J * wt;
            source *= h3;
            source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
          }

          /*
           * Add contributions to this residual (globally into Resid, and
           * locally into an accumulator)
           */

          lec->R[LEC_R_INDEX(peqn, ii)] += mass + advection + porous + diffusion + source;

        } /* end of if (active_dof) */
      }   /* end of for i=0,ei[pg->imtrx]->dof statement */
    }
  }

  /*
   * Jacobian terms...
   */

  if (af->Assemble_Jacobian) {
    for (a = 0; a < WIM; a++) {
      eqn = R_PMOMENTUM1 + a;
      peqn = upd->ep[pg->imtrx][eqn];
      bfm = bf[eqn];

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

          /* Assign pointers into the bf structure */

          grad_phi_i_e_a = grad_phi_e[i][a];

          d_grad_phi_i_e_a_dmesh = bfm->d_grad_phi_e_dmesh[i][a];

          /*
           * J_pm_T
           */

          /* MMH Don't worry about this one yet ... */
          var = TEMPERATURE;
          if (pd->v[pg->imtrx][var]) {
            pvar = upd->vp[pg->imtrx][var];
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              phi_j = bf[var]->phi[j];

              mass = 0.;

              if (pd->TimeIntegration != STEADY) {
                if (pd->e[pg->imtrx][eqn] & T_MASS) {
                  mass = d_rho->T[j] * pv_dot[a];
                  mass *= -phi_i * det_J * wt * h3;
                  mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
                }
                /* if porous flow is considered. KSC on 5/10/95 */
                if (pd->e[pg->imtrx][eqn] & T_POROUS_BRINK) {
                  mass /= por;
                }
                mass *= p_vol_frac;
              }

              /* This is Temperature */
              advection = 0.;
              if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
                for (p = 0; p < WIM; p++) {
                  advection += (pv[p] - x_dot[p]) * grad_pv[p][a];
                }
                if (d_rho->T[j] != 0.0) {
                  printf("d_rho->T[%d]=%g\n", j, d_rho->T[j]);
                  exit(10101);
                }
                advection *= -phi_i * d_rho->T[j] * det_J * wt;
                advection *= h3;
                advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];

                /* if porous flow is considered. KSC on 5/10/95 */
                if (pd->e[pg->imtrx][eqn] & T_POROUS_BRINK) {
                  por2 = por * por;
                  advection /= por2;
                }
                advection *= p_vol_frac;
              }

              porous = 0.;
              if (pd->e[pg->imtrx][eqn] & T_POROUS_BRINK) {
                porous = pv[a] * (d_rho->T[j] * sc * speed / sqrt(per));
                porous *= -phi_i * det_J * wt;
                porous *= h3;
                porous *= pd->etm[pg->imtrx][eqn][(LOG2_POROUS_BRINK)];
              }

              diffusion = 0.;
              if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {
                for (p = 0; p < VIM; p++) {
                  for (q = 0; q < VIM; q++) {
                    diffusion += grad_phi_i_e_a[p][q] * d_Pi_T[q][p][j];
                  }
                }
                diffusion *= -det_J * wt;
                diffusion *= h3;
                diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
              }

              source = 0.0;
              if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
                source = phi_i * dfdT[a][j] * det_J * h3 * wt;
                source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
              }

              lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] +=
                  mass + advection + porous + diffusion + source;
            }
          }

          /*
           * J_pm_pv
           */
          for (b = 0; b < WIM; b++) {
            var = PVELOCITY1 + b;
            if (pd->v[pg->imtrx][var]) {
              pvar = upd->vp[pg->imtrx][var];
              for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                phi_j = bf[var]->phi[j];

                mass = 0.0;
                if (pd->TimeIntegration != STEADY) {
                  if (pd->e[pg->imtrx][eqn] & T_MASS) {
                    mass = (1. + 2. * tt) * phi_j / dt * (double)delta(a, b);
                    mass *= rho;
                    mass *= p_vol_frac;
                    mass *= -phi_i * det_J * h3 * wt;
                    mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
                  }

                  /* if porous flow is considered.
                   * KSC on 5/10/95
                   */
                  if (pd->e[pg->imtrx][eqn] & T_POROUS_BRINK) {
                    mass /= por;
                  }
                }

                porous = 0.;
                if (pd->e[pg->imtrx][eqn] & T_POROUS_BRINK) {
                  porous = ((rho * sc / sqrt(per)) * (2. * pv[b]) * pv[a] +
                            (rho * sc * speed / sqrt(per) + vis / per) * (double)delta(a, b));
                  porous *= -phi_i * phi_j * det_J * wt;
                  porous *= h3;
                  porous *= pd->etm[pg->imtrx][eqn][(LOG2_POROUS_BRINK)];
                }

                /* This is J_pm_pv */
                advection = 0.0;
                if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
                  advection += phi_j * grad_pv[b][a];
                  for (p = 0; p < WIM; p++) {
                    advection += (pv[p] - x_dot[p]) * grad_phi_e[j][b][p][a];
                  }
                  advection *= rho;
                  advection *= p_vol_frac;
                  advection *= -phi_i * det_J * wt;
                  advection *= h3;
                  advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];

                  /* if porous flow is considered. KSC on 5/10/95 */
                  if (pd->e[pg->imtrx][eqn] & T_POROUS_BRINK) {
                    por2 = por * por;
                    advection /= por2;
                  }
                }

                /* J_pm_pv */
                diffusion = 0.;
                if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {
                  for (p = 0; p < VIM; p++) {
                    for (q = 0; q < VIM; q++) {
                      diffusion += grad_phi_i_e_a[p][q] * d_Pi_pv[q][p][b][j];
                    }
                  }
                  diffusion *= -det_J * wt;
                  diffusion *= h3;
                  diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
                }

                source = 0.0;
                if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
                  source = phi_i * dfdv[a][b][j] * det_J * h3 * wt;
                  source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
                }

                lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] +=
                    mass + advection + porous + diffusion + source;
              }
            }
          }

          /*
           * J_pm_c
           */
          var = MASS_FRACTION;
          if (pd->v[pg->imtrx][var]) {
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {

              phi_j = bf[var]->phi[j];

              for (w = 0; w < pd->Num_Species_Eqn; w++) {

                mass = 0.0;
                if (pd->TimeIntegration != STEADY) {
                  if (pd->e[pg->imtrx][eqn] & T_MASS) {
                    /*mass = d_rho->C[w][j] * pv_dot[a];*/
                    if (w == species)
                      mass = rho * pv_dot[a];
                    else
                      mass = d_rho->C[w][j] * pv_dot[a];
                    mass *= -phi_i * det_J * wt * h3;
                    mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
                  }
                  /* if porous flow is considered. KSC on 5/10/95 */
                  if (pd->e[pg->imtrx][eqn] & T_POROUS_BRINK) {
                    mass /= por;
                  }
                }

                /* This is species */
                advection = 0.0;
                if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
                  for (p = 0; p < WIM; p++) {
                    advection += (pv[p] - x_dot[p]) * grad_pv[p][a];
                  }

                  if (w == species)
                    advection *= rho;
                  else
                    advection *= d_rho->C[w][j];

                  advection *= -phi_i * det_J * wt * h3;

                  advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];

                  /* if porous flow is considered. KSC on 5/10/95 */
                  if (pd->e[pg->imtrx][eqn] & T_POROUS_BRINK) {
                    por2 = por * por;
                    advection /= por2;
                  }
                }

                porous = 0.;
                if (pd->e[pg->imtrx][eqn] & T_POROUS_BRINK) {
                  porous = pv[a] * (d_rho->C[w][j] * sc * speed / sqrt(per));
                  porous *= -phi_i * det_J * wt * h3;
                  porous *= pd->etm[pg->imtrx][eqn][(LOG2_POROUS_BRINK)];
                }

                diffusion = 0.;
                if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {
                  for (p = 0; p < VIM; p++) {
                    for (q = 0; q < VIM; q++) {
                      diffusion += grad_phi_i_e_a[p][q] * d_Pi_C[q][p][w][j];
                    }
                  }
                  diffusion *= -det_J * wt;
                  diffusion *= h3;
                  diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
                }

                source = 0.0;
                if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
                  source = phi_i * dfdC[a][w][j] * det_J * h3 * wt;
                  source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
                }

                if (w > 1) {
                  GOMA_EH(GOMA_ERROR, "Need more arrays for each species.");
                }

                lec->J[LEC_J_INDEX(peqn, MAX_PROB_VAR + w, ii, j)] +=
                    mass + advection + porous + diffusion + source;
              }
            }
          }

          /*
           * Pressure isn't used for the particle momentum equations.
           * J_pm_P
           */
          /*
            var = PRESSURE;
            if ( pd->v[pg->imtrx][var] )
            {
            pvar = upd->vp[pg->imtrx][var];

            for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
            {

            begin MMH commented-out section; dal clarify for AIX ***

            phi_j = bf[var]->phi[j];

            mass    = 0.;

            porous    = 0.;

            advection = 0.;

            diffusion = 0.;

            if ( pd->e[pg->imtrx][eqn] & T_DIFFUSION )
            {
            for ( p=0; p<VIM; p++)
            {
            for ( q=0; q<VIM; q++)
            {

            diffusion -=  grad_phi_i_e_a[p][q] * d_Pi_P[q][p][j];

            }
            }

            diffusion *= det_J * h3 * wt;
            diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
            }

            source    = 0.;

            lec->J[LEC_J_INDEX(peqn,pvar,ii,j)] +=
            mass + advection + porous + diffusion + source;
            }
            }
            *** end of MMH section  *** */

          /* MMH Shear rate is currently not used as a variable for
           * SUSPENSION_PM flow.
           * J_pm_S
           */

          /* begin MMH commented-out section; dal clarify for AIX ***
             for ( b=0; b<VIM; b++)
             {
             for ( c=0; c<VIM; c++)
             {
             var = pv_s[b][c];

             if ( pd->v[pg->imtrx][var] )
             {
             pvar = upd->vp[pg->imtrx][var];

             for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
             {

             phi_j = bf[var]->phi[j];


             mass = 0.;

             porous    = 0.;

             advection = 0.;

             diffusion = 0.;

             if ( pd->e[pg->imtrx][eqn] & T_DIFFUSION )
             {

             start comment
             for ( p=0; p<VIM; p++)
             {
             for ( q=0; q<VIM; q++)
             {
             diffusion =
             grad_phi_i_e_a[p][q] *  (double)delta(c,p) * (double)delta(b,q);
             }
             } *
             diffusion = -grad_phi_i_e_a[c][b];
             diffusion *= phi_j * det_J * wt *h3;
             diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
             }

             source    = 0.;

             mass *= p_vol_frac;
             advection *= p_vol_frac;

             lec->J[LEC_J_INDEX(peqn,pvar,ii,j)] +=
             mass + advection + porous + diffusion + source;
             }
             }
             }
             }
             *** end of MMH section  *** */

          /* MMH Not used for SUSPENSION_PM flow.
           * J_pm_G
           */
          /*
            if ( pd->v[pg->imtrx][POLYMER_STRESS11] && (vn->evssModel == EVSS_F) )
            {
            for ( b=0; b<VIM; b++)
            {
            for ( c=0; c<VIM; c++)
            {
            var = pv_g[b][c];

            if ( pd->v[pg->imtrx][var] )
            {
            pvar = upd->vp[pg->imtrx][var];

            for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
            {

            phi_j = bf[var]->phi[j];

            diffusion = 0.;

            if ( pd->e[pg->imtrx][eqn] & T_DIFFUSION )
            {
            for ( p=0; p<VIM; p++)
            {
            for ( q=0; q<VIM; q++)
            {
            diffusion +=
            grad_phi_i_e_a[p][q] *
            evss_f * mu * ((double)delta(c,p) * (double)delta(b,q) + (double)delta(b,p) *
            (double)delta(c,q));
            }
            }
            diffusion *= phi_j * det_J * wt *h3;
            diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
            }

            mass *= p_vol_frac;
            advection *= p_vol_frac;

            lec->J[LEC_J_INDEX(peqn,pvar,ii,j)] +=
            diffusion;
            }
            }
            }
            }
            }
            */
          /*
           * J_pm_d
           */
          for (b = 0; b < dim; b++) {
            var = MESH_DISPLACEMENT1 + b;
            if (pd->v[pg->imtrx][var]) {
              pvar = upd->vp[pg->imtrx][var];
              for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                phi_j = bf[var]->phi[j];

                d_det_J_dmesh_bj = bfm->d_det_J_dm[b][j];

                dh3dmesh_bj = fv->dh3dq[b] * phi_j;

                /* d_mu_dmesh_bj  = d_mu_dmesh [b][j]; */

                mass = 0.;

                if (pd->TimeIntegration != STEADY) {
                  if (pd->e[pg->imtrx][eqn] & T_MASS) {
                    mass = pv_dot[a];
                    mass *= -phi_i * rho * (d_det_J_dmesh_bj * h3 + det_J * dh3dmesh_bj) * wt;
                    mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
                  }

                  if (pd->e[pg->imtrx][eqn] &
                      T_POROUS_BRINK) /* if porous flow is considered. KSC on 5/10/95 */
                  {
                    mass /= por;
                  }
                }

                porous = 0.;
                if (pd->e[pg->imtrx][eqn] & T_POROUS_BRINK) {
                  for (p = 0; p < WIM; p++) {
                    porous +=
                        pv[p] * (rho * sc * speed / sqrt(per) + vis / per) * (double)delta(p, a);
                    porous *= -phi_i * wt * (d_det_J_dmesh_bj * h3 + det_J * dh3dmesh_bj);
                    porous *= pd->etm[pg->imtrx][eqn][(LOG2_POROUS_BRINK)];
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

                  advection_a = 0.;
                  for (p = 0; p < WIM; p++) {
                    advection_a += (pv[p] - x_dot[p]) * fv->d_grad_pv_dmesh[p][a][b][j];
                  }
                  advection_a *= -phi_i * rho * h3 * det_J * wt;

                  advection_b = 0.;
                  for (p = 0; p < WIM; p++) {
                    advection_b += (pv[p] - x_dot[p]) * grad_pv[p][a];
                  }
                  advection_b *= -phi_i * rho * wt * (d_det_J_dmesh_bj * h3 + det_J * dh3dmesh_bj);

                  advection_c = 0.;
                  if (pd->TimeIntegration != STEADY) {
                    if (pd->e[pg->imtrx][eqn] & T_MASS) {

                      for (p = 0; p < WIM; p++) {
                        advection_c +=
                            (-(1. + 2. * tt) * phi_j / dt * (double)delta(p, b)) * grad_pv[p][a];
                      }
                      advection_c *= -phi_i * rho * wt * h3 * det_J;
                    }
                  }

                  advection = advection_a + advection_b + advection_c;

                  advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
                  if (pd->e[pg->imtrx][eqn] &
                      T_POROUS_BRINK) /* if porous flow is considered. KSC on 5/10/95 */
                  {
                    por2 = por * por;
                    advection /= por2;
                  }
                }

                /*
                 * Diffusion...
                 */

                diffusion = 0.;

                if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {

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

                  for (p = 0; p < VIM; p++) {
                    for (q = 0; q < VIM; q++) {
                      diff_a += d_grad_phi_i_e_a_dmesh[p][q][b][j] * Pi[q][p];

                      diff_b += grad_phi_i_e_a[p][q] * d_Pi_mesh[q][p][b][j];

                      diff_c += grad_phi_i_e_a[p][q] * Pi[q][p];
                    }
                  }
                  diff_a *= -det_J * h3 * wt;
                  diff_b *= -det_J * h3 * wt;
                  diff_c *= -wt * (d_det_J_dmesh_bj * h3 +

                                   det_J * dh3dmesh_bj);
                  diffusion = diff_a + diff_b + diff_c;

                  diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
                }

                /*
                 * Source term...
                 */

                source = 0.;

                if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
                  source += phi_i * wt *
                            (f[a] * d_det_J_dmesh_bj * h3 + f[a] * det_J * dh3dmesh_bj +
                             dfdX[a][b][j] * det_J * h3);

                  source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
                }

                /* MMH */
                mass *= p_vol_frac;
                advection *= p_vol_frac;

                lec->J[LEC_J_INDEX(peqn, pvar, ii, j)] +=
                    mass + advection + porous + diffusion + source;
              }
            }
          }
        } /* end of if(active_dof) */
      }   /* end_of for(i=0,ei[pg->imtrx]->dof[eqn])*/
    }
  }

  return (status);
}

/* MMH_assemble_continuity -- assemble Residual &| Jacobian for continuity eqns
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

/* MMH
 * I want to modify this to include the particle volume fraction for
 * the particle phase after I take care of the particle momentum
 * equations.
 */

int MMH_assemble_continuity(double time_value,       /* current value of time */
                            double tt,               /* parameter to vary time integration from
                                                        explicit (tt = 1) to implicit (tt = 0) */
                            double dt,               /* current time step size */
                            double h_elem_avg,       /* average global element size for PSPG,
                                                        taken to be constant wrt to Jacobian entries */
                            double hsquared[],       /* (DIM) element size information for PSPG */
                            double hh[][DIM],        /* (DIM)(DIM) these aren't currently used, but
                                                        left in just in case they are needed  later */
                            double dh_dxnode[][MDE], /* (DIM)(MDE) */
                            double U_norm, /* global velocity norm for PSPG calculations */
                            double mu_avg) /* element viscosity for PSPG calculations */

{
  int dim;
  int p, q, a, b;

  int eqn, var;
  int peqn, pvar;
  int w;

  int i, j;
  int status, err;

  dbl v[DIM]; /* Velocity field. */
  dbl div_pv; /* Divergence of v. */

  dbl P; /* Pressure. */

  dbl mu = 0.0; /* Viscosity. */

  dbl advection;
  dbl source;
  dbl pressure_stabilization;

  dbl volsolvent = 1e12; /* volume fraction of solvent */

  /*
   * Initial solvent volume fraction (in stress-free state) input as source
   * constant from input file
   */

  dbl initial_volsolvent = elc->Strss_fr_sol_vol_frac;

  dbl det_J;
  dbl h3;

  dbl d_h3detJ_dmesh_bj; /* for specific (b,j) mesh dof */

  /*
   * Galerkin weighting functions...
   */

  dbl phi_i;

  /*
   * Interpolation functions...
   */

  dbl phi_j;
  dbl div_phi_j_e_b;

  dbl div_pv_dmesh = 0.; /* for specific (b,j) mesh dof */

  dbl wt;

  /*
   * Variables for Pressure Stabilization Petrov-Galerkin...
   */
  int meqn;
  int var1;
  int r;
  int pv_s[DIM][DIM], pv_g[DIM][DIM];

  dbl mass;
  dbl diffusion;
  dbl source_a;
  dbl advection_a;
  dbl stress;
  dbl pressure;
  dbl velocity_gradient;
  dbl stabilization_a;
  dbl stabilization_b;
  dbl momentum_residual[DIM] = {0.0, 0.0, 0.0}; /* momentum residual for PSPG */
  dbl x_dot[DIM] = {0.0, 0.0, 0.0};
  dbl pv_dot[DIM] = {0.0, 0.0, 0.0};
  dbl grad_P[DIM];
  dbl grad_pv[DIM][DIM];
  dbl div_s[DIM];
  dbl div_G[DIM] = {0.0, 0.0, 0.0};
  dbl grad_phi[MDE][DIM]; /* weight-function for PSPG term */

  /* variables for Brinkman porous flow */
  dbl por = 0.0, por2 = 0.0, per = 0.0, vis, sc = 0.0, speed = 0.0;
  dbl porous;

  dbl h_elem, h_elem_inv;
  dbl rho = 0.;
  dbl Re;
  dbl tau_pspg = 0.;
  dbl d_tau_pspg_dm[DIM][MDE];

  dbl f[DIM];                   /* Body force. */
  dbl dfdT[DIM][MDE];           /* For temperature dependence. */
  dbl dfdX[DIM][DIM][MDE];      /* For spatial dependence. */
  dbl dfdv[DIM][DIM][MDE];      /* For velocity dependence. */
  dbl dfdC[DIM][MAX_CONC][MDE]; /* For concentration dependence. */

  /*
   * Variables for vicosity and derivative
   */
  VISCOSITY_DEPENDENCE_STRUCT d_mu_struct; /* viscosity dependence */
  VISCOSITY_DEPENDENCE_STRUCT *d_mu = &d_mu_struct;

  dbl gamma[DIM][DIM]; /* shrearrate tensor based on velocity */

  /*
   * Species diffusive flux and sensitivity terms
   */
  /* density derivatives */
  DENSITY_DEPENDENCE_STRUCT d_rho_struct; /* density dependence */
  DENSITY_DEPENDENCE_STRUCT *d_rho = &d_rho_struct;

  struct Species_Conservation_Terms s_terms;
  dbl rhos = 0;
  dbl rhof = 0;
  dbl h_flux = 0;

  int w0 = -1;

  status = 0;

  memset(d_mu, 0, sizeof(VISCOSITY_DEPENDENCE_STRUCT));
  memset(d_rho, 0, sizeof(DENSITY_DEPENDENCE_STRUCT));
  memset(grad_pv, 0, sizeof(double) * DIM * DIM);
  memset(dfdT, 0, sizeof(double) * DIM * MDE);
  memset(dfdX, 0, sizeof(double) * DIM * DIM * MDE);
  memset(dfdv, 0, sizeof(double) * DIM * DIM * MDE);
  memset(dfdC, 0, sizeof(double) * DIM * MAX_CONC * MDE);

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

  pv_s[0][0] = POLYMER_STRESS11;
  pv_s[0][1] = POLYMER_STRESS12;
  pv_s[0][2] = POLYMER_STRESS13;
  pv_s[1][0] = POLYMER_STRESS12;
  pv_s[1][1] = POLYMER_STRESS22;
  pv_s[1][2] = POLYMER_STRESS23;
  pv_s[2][0] = POLYMER_STRESS13;
  pv_s[2][1] = POLYMER_STRESS23;
  pv_s[2][2] = POLYMER_STRESS33;

  pv_g[0][0] = VELOCITY_GRADIENT11;
  pv_g[0][1] = VELOCITY_GRADIENT12;
  pv_g[1][0] = VELOCITY_GRADIENT21;
  pv_g[1][1] = VELOCITY_GRADIENT22;
  pv_g[0][2] = VELOCITY_GRADIENT13;
  pv_g[1][2] = VELOCITY_GRADIENT23;
  pv_g[2][0] = VELOCITY_GRADIENT31;
  pv_g[2][1] = VELOCITY_GRADIENT32;
  pv_g[2][2] = VELOCITY_GRADIENT33;

  wt = fv->wt;

  det_J = bf[eqn]->detJ; /* Really, ought to be mesh eqn. */

  h3 = fv->h3; /* Differential volume element (scales). */

  for (a = 0; a < dim; a++) {
    v[a] = fv->v[a];
  }

  P = fv->P;
  div_pv = fv->div_pv;

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

  if (cr->MeshMotion == TOTAL_ALE && pd->e[pg->imtrx][R_SOLID1]) {
    err = belly_flop_rs(elc_rs->lame_mu);
    GOMA_EH(err, "error in belly flop for real solid");
    if (err == 2)
      return (err);
  }

  if (PSPG && 0) {
    /*      h_elem = 0.;
            for ( p=0; p<dim; p++)
            {
            h_elem += h[p];
            }
            h_elem = sqrt(h_elem)/2.; */

    /* use global average for element size */

    h_elem = h_elem_avg;
    h_elem_inv = 1. / h_elem;

    /*** Density ***/
    rho = density(d_rho, time_value);

    /* Now calculate the element Reynolds number based on a global
       norm of the velocity */

    Re = rho * U_norm * h_elem / (2.0 * mu_avg);

    if (Re <= 3.0) {
      tau_pspg = PS_scaling * h_elem * h_elem / (12.0 * mu_avg);
    } else if (Re > 3.0) {
      tau_pspg = PS_scaling * h_elem / (2.0 * rho * U_norm);
    }

    /* load up shearrate tensor based on velocity */
    for (a = 0; a < VIM; a++) {
      for (b = 0; b < VIM; b++) {
        gamma[a][b] = fv->grad_pv[a][b] + fv->grad_pv[b][a];
      }
    }

    /* get viscosity for velocity second derivative/diffusion term in PSPG stuff */
    mu = viscosity(gn, gamma, d_mu);

    /* set up mesh derivative for tau_pspg, if necessary */
    if (pd->v[pg->imtrx][MESH_DISPLACEMENT1]) {
      for (b = 0; b < dim; b++) {
        var = MESH_DISPLACEMENT1 + b;
        if (pd->v[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            d_tau_pspg_dm[b][j] = 0.;
            for (w = 0; w < dim; w++) {
              if (Re <= 3.0) {
                d_tau_pspg_dm[b][j] += PS_scaling * hh[w][b] * dh_dxnode[w][j] / (24.0 * mu_avg);
              } else {
                d_tau_pspg_dm[b][j] +=
                    PS_scaling * hh[w][b] * dh_dxnode[w][j] * h_elem_inv / (8.0 * rho * U_norm);
              }
            }
          }
        }
      }
    }

    /* get variables we will need for assembly */

    for (a = 0; a < WIM; a++) {
      v[a] = fv->v[a];
      grad_P[a] = fv->grad_P[a];
      if (pd->TimeIntegration != STEADY && pd->v[pg->imtrx][MESH_DISPLACEMENT1 + a]) {
        x_dot[a] = fv_dot->x[a];
      } else {
        x_dot[a] = 0.;
      }
      if (pd->TimeIntegration != STEADY) {
        pv_dot[a] = fv_dot->v[a];
      } else {
        pv_dot[a] = 0.;
      }
    }

    for (a = 0; a < VIM; a++) {
      for (b = 0; b < VIM; b++) {
        grad_pv[a][b] = fv->grad_pv[a][b];
      }
    }

    if (pd->v[pg->imtrx][POLYMER_STRESS11]) {
      for (p = 0; p < WIM; p++) {
        div_s[p] = fv->div_S[0][p];
      }
    } else {
      for (p = 0; p < WIM; p++) {
        div_s[p] = 0.;
      }
    }

    if (pd->v[pg->imtrx][VELOCITY_GRADIENT11]) {
      for (p = 0; p < WIM; p++) {
        div_G[p] = fv->div_G[p];
      }
    } else {
      for (p = 0; p < WIM; p++) {
        div_G[p] = 0.;
      }
    }
    if (pd->e[pg->imtrx][R_MOMENTUM1] & T_POROUS_BRINK) {
      if (mp->PorousMediaType != POROUS_BRINKMAN)
        GOMA_WH(-1, "Set Porous term multiplier in continuous medium");
      /* Short-hand notation for the four parameters in the Brinkman Equation. */
      por = mp->porosity;
      por2 = por * por;
      per = mp->permeability;
      vis = mp->FlowingLiquid_viscosity;
      sc = mp->Inertia_coefficient;
    } else {
      por = 1.;
      por2 = 1.;
      per = 1.;
      vis = mp->viscosity;
      sc = 0.;
    }

    /* for porous media stuff */
    speed = 0.0;
    for (a = 0; a < WIM; a++) {
      speed += v[a] * v[a];
    }
    speed = sqrt(speed);

    /* get momentum source term */
    err = pmomentum_source_term(f, dfdT, dfdX, dfdC, dfdv);

    for (a = 0; a < WIM; a++) {
      meqn = R_PMOMENTUM1 + a;
      momentum_residual[a] = rho * pv_dot[a] / por * pd->etm[pg->imtrx][meqn][(LOG2_MASS)] +
                             grad_P[a] * pd->etm[pg->imtrx][meqn][(LOG2_DIFFUSION)] -
                             div_s[a] * pd->etm[pg->imtrx][meqn][(LOG2_DIFFUSION)] -
                             mu * div_G[a] * pd->etm[pg->imtrx][meqn][(LOG2_DIFFUSION)] -
                             f[a] * pd->etm[pg->imtrx][meqn][(LOG2_SOURCE)] +
                             v[a] * (rho * sc * speed / sqrt(per) + vis / per) *
                                 pd->etm[pg->imtrx][meqn][(LOG2_POROUS_BRINK)];
      for (b = 0; b < WIM; b++) {
        momentum_residual[a] += rho * (v[b] - x_dot[b]) * grad_pv[b][a] / por2 *
                                pd->etm[pg->imtrx][meqn][(LOG2_ADVECTION)];
      }
    }
    for (a = 0; a < VIM; a++) {
      for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
        grad_phi[i][a] = bf[eqn]->grad_phi[i][a];
      }
    }
  } /* end if (PSPG && 0) */

  if ((cr->MeshMotion == LAGRANGIAN || cr->MeshMotion == DYNAMIC_LAGRANGIAN ||
       cr->MeshMotion == TOTAL_ALE) &&
      (!mp->PorousMediaType)) {
    initial_volsolvent = elc->Strss_fr_sol_vol_frac;
    volsolvent = 0.;
    for (w = 0; w < pd->Num_Species_Eqn; w++)
      volsolvent += fv->c[w];
  }

  if ((cr->MassFluxModel == HYDRODYNAMIC) && (mp->DensityModel == SUSPENSION) &&
      (mp->MomentumSourceModel == SUSPENSION)) {
    /*
     * Compute hydrodynamic/sedimentation flux and sensitivities.
     */

    w0 =
        (int)mp->u_density[0]; /* This is the species number that is transported HYDRODYNAMICally */

    hydro_flux(&s_terms, w0, tt, dt, hsquared);

    rhof = mp->u_density[1];
    rhos = mp->u_density[2];

    for (a = 0; a < VIM; a++) {
      for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
        grad_phi[i][a] = bf[eqn]->grad_phi[i][a];
      }
    }
  }

  if (af->Assemble_Residual) {
    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      phi_i = bf[eqn]->phi[i];

      advection = 0.;

      if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {

        if (pd->v[pg->imtrx]
                 [PVELOCITY1]) /* then must be solving fluid mechanics in this material */
        {

          /*
           * Standard incompressibility constraint means we have
           * a solenoidal velocity field...if we don't, then
           * we might be in serious trouble...
           */

          advection = div_pv;

          advection *= phi_i * h3 * det_J * wt;

          advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];

        } else if (cr->MeshMotion == LAGRANGIAN || cr->MeshMotion == DYNAMIC_LAGRANGIAN ||
                   cr->MeshMotion == TOTAL_ALE)
        /* use divergence of displacement for linear elasticity */
        {
          advection = fv->volume_change;

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

        if (pd->v[pg->imtrx][PVELOCITY1]) {
          source = P;
          source *= phi_i * h3 * det_J * wt;
          source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
        }

        if ((cr->MeshMotion == LAGRANGIAN || cr->MeshMotion == DYNAMIC_LAGRANGIAN ||
             cr->MeshMotion == TOTAL_ALE))
        /* add swelling as a source of volume */
        {
          if (!mp->PorousMediaType) {
            source = -(1. - initial_volsolvent) / (1. - volsolvent);
            source *= phi_i * h3 * det_J * wt;
            source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
          }
        }
      }
      /* add Pressure-Stabilized Petrov-Galerkin term
       * if desired.
       */
      pressure_stabilization = 0.;
      if (PSPG && 0) {
        for (a = 0; a < WIM; a++) {
          meqn = R_PMOMENTUM1 + a;
          if (pd->e[pg->imtrx][meqn]) {
            pressure_stabilization += grad_phi[i][a] * momentum_residual[a];
          }
        }
        pressure_stabilization *= tau_pspg * h3 * det_J * wt;
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

      lec->R[LEC_R_INDEX(peqn, i)] += advection + source + pressure_stabilization + h_flux;
    }
  }

  if (af->Assemble_Jacobian) {
    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      phi_i = bf[eqn]->phi[i];

      /*
       * J_c_v NOTE that this is applied whenever velocity is a variable
       */
      for (b = 0; b < WIM; b++) {
        var = PVELOCITY1 + b;
        if (pd->v[pg->imtrx][var]) {
          pvar = upd->vp[pg->imtrx][var];
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];

            h_flux = 0.;
            if ((cr->MassFluxModel == HYDRODYNAMIC) && (mp->MomentumSourceModel == SUSPENSION)) {
              for (p = 0; p < dim; p++) {
                h_flux += grad_phi[i][p] * s_terms.d_diff_flux_dv[w0][p][b][j] * det_J * h3;
              }
              h_flux *= (rhos - rhof) / rhof * wt * pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
            }

            advection = 0.;

            if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
              div_phi_j_e_b = 0.;
              for (p = 0; p < VIM; p++) {
                div_phi_j_e_b += bf[var]->grad_phi_e[j][b][p][p];
              }

              advection = phi_i * div_phi_j_e_b * h3 * det_J * wt;

              advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
            }

            source = 0.;

            /* add Pressure-Stabilized Petrov-Galerkin term
             * if desired.
             */
            pressure_stabilization = 0.;
            if (PSPG && 0) {
              for (a = 0; a < WIM; a++) {
                meqn = R_PMOMENTUM1 + a;

                mass = 0.;
                if (pd->TimeIntegration != STEADY) {
                  if (pd->e[pg->imtrx][meqn] & T_MASS) {
                    mass += (1. + 2. * tt) * phi_j / dt * (double)delta(a, b);
                    mass *= rho / por * pd->etm[pg->imtrx][meqn][(LOG2_MASS)];
                  }
                }

                diffusion = 0.;
                if (pd->e[pg->imtrx][meqn] & T_DIFFUSION) {
                  diffusion -= d_mu->v[b][j] * div_G[a];

                  diffusion *= pd->etm[pg->imtrx][meqn][(LOG2_DIFFUSION)];
                }

                advection_a = 0.;
                if (pd->e[pg->imtrx][meqn] & T_ADVECTION) {
                  advection_a += phi_j * grad_pv[b][a];
                  for (p = 0; p < WIM; p++) {
                    advection_a += (v[p] - x_dot[p]) * bf[var]->grad_phi_e[j][b][p][a];
                  }
                  advection_a *= rho / por2 * pd->etm[pg->imtrx][meqn][(LOG2_ADVECTION)];
                }

                source_a = 0.;
                if (pd->e[pg->imtrx][meqn] & T_SOURCE) {
                  source_a -= dfdv[a][b][j] * pd->etm[pg->imtrx][meqn][(LOG2_SOURCE)];
                }

                porous = 0.;
                if (pd->e[pg->imtrx][meqn] & T_POROUS_BRINK) {
                  for (p = 0; p < WIM; p++) {
                    porous +=
                        (rho * sc / sqrt(per) * (2. * v[p]) * v[a] +
                         (rho * sc * speed / sqrt(per) + vis / per) * (double)delta(a, p) * phi_i) *
                        pd->etm[pg->imtrx][meqn][(LOG2_POROUS_BRINK)];
                  }
                }

                /* MMH */
                pressure_stabilization +=
                    0 * (mass + diffusion + advection_a + source_a + porous) * grad_phi[i][a];
              }
              /* MMH */
              /*
              pressure_stabilization *= tau_pspg * h3 * det_J * wt;
              */
            }

            /* MMH */
            lec->J[LEC_J_INDEX(peqn, pvar, i, j)] +=
                0 * (advection + source + pressure_stabilization + h_flux);
          }
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
          porous = 0.;
          source = 0.;
          diffusion = 0.;
          mass = 0.;
          advection = 0.;

          for (a = 0; a < WIM; a++) {
            meqn = R_PMOMENTUM1 + a;

            if (pd->TimeIntegration != STEADY) {
              if (pd->e[pg->imtrx][meqn] & T_MASS) {
                mass = d_rho->T[j] / por * pv_dot[a] * grad_phi[i][a] *
                       pd->etm[pg->imtrx][meqn][(LOG2_MASS)];
              }
            }

            if (pd->e[pg->imtrx][meqn] & T_ADVECTION) {
              advection = 0.;
              for (p = 0; p < WIM; p++) {
                advection += (v[p] - x_dot[p]) * grad_pv[p][a];
              }
              advection *=
                  grad_phi[i][a] * d_rho->T[j] / por2 * pd->etm[pg->imtrx][meqn][(LOG2_ADVECTION)];
            }

            if (pd->e[pg->imtrx][meqn] & T_POROUS_BRINK) {
              porous = v[a] * (d_rho->T[j] * sc * speed / sqrt(per)) * grad_phi[i][a] *
                       pd->etm[pg->imtrx][meqn][(LOG2_POROUS_BRINK)];
            }

            if (pd->e[pg->imtrx][meqn] & T_SOURCE) {
              source = -grad_phi[i][a] * dfdT[a][j] * pd->etm[pg->imtrx][meqn][(LOG2_SOURCE)];
            }

            if (pd->e[pg->imtrx][meqn] & T_DIFFUSION) {
              diffusion = -grad_phi[i][a] * d_mu->T[j] * div_G[a] *
                          pd->etm[pg->imtrx][meqn][(LOG2_DIFFUSION)];
            }
            pressure_stabilization += mass + advection + porous + source + diffusion;
          }

          pressure_stabilization *= tau_pspg * h3 * det_J * wt;

          lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += pressure_stabilization;
        }
      }
      if ((cr->MassFluxModel == HYDRODYNAMIC) && pd->v[pg->imtrx][var]) {
        if (mp->MomentumSourceModel == SUSPENSION) {
          pvar = upd->vp[pg->imtrx][var];

          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            h_flux = 0.;

            for (p = 0; p < dim; p++) {
              h_flux += grad_phi[i][p] * s_terms.d_diff_flux_dT[w0][p][j];
            }

            h_flux *= h3 * det_J * wt * (rhos - rhof) / rhof;

            /*  h_flux = 0.0; */

            lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += h_flux;
          }
        }
      }

      /*
       * J_c_P here species act as a volume source in continuous lagrangian mesh motion
       */
      var = PRESSURE;
      if (pd->v[pg->imtrx][var]) {
        pvar = upd->vp[pg->imtrx][var];

        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          phi_j = bf[var]->phi[j];

          advection = 0.;

          if ((cr->MeshMotion == LAGRANGIAN || cr->MeshMotion == DYNAMIC_LAGRANGIAN ||
               cr->MeshMotion == TOTAL_ALE) &&
              pd->e[pg->imtrx][eqn] & T_ADVECTION) {
            /*Need to compute this for total ALE.  Not done yet */
            advection = fv->d_volume_change_dp[j];

            advection *= phi_i * h3 * det_J * wt;

            advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
          }

          source = 0.;
          if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
            if (pd->v[pg->imtrx][PVELOCITY1]) {
              source = phi_j * h3 * det_J * wt;

              source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
            }
          }

          /* add Pressure-Stabilized Petrov-Galerkin term
           * if desired.
           */
          pressure_stabilization = 0.;
          if (PSPG) {
            for (a = 0; a < WIM; a++) {
              meqn = R_PMOMENTUM1 + a;
              if (pd->e[pg->imtrx][meqn] & T_DIFFUSION) {
                pressure_stabilization += grad_phi[i][a] *
                                          (bf[var]->grad_phi[j][a] - d_mu->P[j] * div_G[a]) *
                                          pd->etm[pg->imtrx][meqn][(LOG2_DIFFUSION)];
              }
            }
            pressure_stabilization *= tau_pspg * h3 * det_J * wt;
          }

          lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += advection + source + pressure_stabilization;
        }
      }

      /*
       * J_c_S this term is only present for PSPG
       */
      var = POLYMER_STRESS11;
      if (PSPG && pd->v[pg->imtrx][var]) {
        pvar = upd->vp[pg->imtrx][var];

        for (p = 0; p < VIM; p++) {
          for (q = 0; q < VIM; q++) {
            var = pv_s[p][q];
            if (pd->v[pg->imtrx][var]) {
              for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                phi_j = bf[var]->phi[j];

                pressure_stabilization = -grad_phi[i][q] * bf[var]->grad_phi[j][p] *
                                         pd->etm[pg->imtrx][R_PMOMENTUM1 + q][(LOG2_DIFFUSION)];

                if (pd->CoordinateSystem != CARTESIAN) {
                  for (r = 0; r < VIM; r++) {
                    pressure_stabilization -=
                        grad_phi[i][q] * phi_j * fv->grad_e[p][r][q] *
                        pd->etm[pg->imtrx][R_PMOMENTUM1 + a][(LOG2_DIFFUSION)];
                  }
                  for (a = 0; a < WIM; a++) {
                    pressure_stabilization -=
                        grad_phi[i][a] * phi_j * fv->grad_e[q][p][a] *
                        pd->etm[pg->imtrx][R_PMOMENTUM1 + a][(LOG2_DIFFUSION)];
                  }
                }

                pressure_stabilization *= tau_pspg * h3 * det_J * wt;

                lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += pressure_stabilization;
              }
            }
          }
        }
      }

      /*
       * J_c_G this term is only present for PSPG
       */
      var = VELOCITY_GRADIENT11;
      if (PSPG && pd->v[pg->imtrx][var]) {
        pvar = upd->vp[pg->imtrx][var];

        for (p = 0; p < VIM; p++) {
          for (q = 0; q < VIM; q++) {
            var = pv_g[p][q];
            if (pd->v[pg->imtrx][var]) {
              for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                phi_j = bf[var]->phi[j];

                pressure_stabilization = -grad_phi[i][q] * bf[var]->grad_phi[j][p] *
                                         pd->etm[pg->imtrx][R_PMOMENTUM1 + q][(LOG2_DIFFUSION)];

                if (pd->CoordinateSystem != CARTESIAN) {
                  for (r = 0; r < VIM; r++) {
                    pressure_stabilization -=
                        grad_phi[i][q] * phi_j * fv->grad_e[p][r][q] *
                        pd->etm[pg->imtrx][R_PMOMENTUM1 + a][(LOG2_DIFFUSION)];
                  }
                  for (a = 0; a < WIM; a++) {
                    pressure_stabilization -=
                        grad_phi[i][a] * phi_j * fv->grad_e[q][p][a] *
                        pd->etm[pg->imtrx][R_PMOMENTUM1 + a][(LOG2_DIFFUSION)];
                  }
                }

                pressure_stabilization *= mu * tau_pspg * h3 * det_J * wt;

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

      if (cr->MassFluxModel == HYDRODYNAMIC && pd->v[pg->imtrx][var]) {
        if (mp->MomentumSourceModel == SUSPENSION) {
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
        if (pd->v[pg->imtrx][var]) {
          pvar = upd->vp[pg->imtrx][var];

          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];
            /* derivative of |J| with extra term for axisymmetry e.g.
               d/dmesh [ r|J| ] */
            d_h3detJ_dmesh_bj =
                h3 * bf[eqn]->d_det_J_dm[b][j] + det_J * fv->dh3dq[b] * bf[var]->phi[j];

            advection = 0.;

            if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
              if (pd->v[pg->imtrx][PVELOCITY1]) {
                h_flux = 0.;
                if ((cr->MassFluxModel == HYDRODYNAMIC) &&
                    (mp->MomentumSourceModel == SUSPENSION)) {
                  for (p = 0; p < dim; p++) {
                    h_flux += grad_phi[i][p] * s_terms.diff_flux[w0][p] * d_h3detJ_dmesh_bj +
                              grad_phi[i][p] * s_terms.d_diff_flux_dmesh[w0][p][b][j] * det_J * h3 +
                              bf[eqn]->d_grad_phi_dmesh[i][p][b][j] * s_terms.diff_flux[w0][p] *
                                  det_J * h3;
                  }
                  h_flux *= (rhos - rhof) / rhof * wt * pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];

                  /*  h_flux = 0.0; */
                }

                /*
                div_pv_dmesh = fv->d_div_pv_dmesh[b][j];
                */

                advection += div_pv_dmesh * det_J * h3 + div_pv * (d_h3detJ_dmesh_bj);

                advection *= phi_i * wt;
              } else if (cr->MeshMotion == LAGRANGIAN || cr->MeshMotion == DYNAMIC_LAGRANGIAN ||
                         cr->MeshMotion == TOTAL_ALE) {
                advection += fv->volume_change * (d_h3detJ_dmesh_bj);

                advection += fv->d_volume_change_dx[b][j] * h3 * det_J;

                advection *= phi_i * wt;
              }

              advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
            }

            source = 0.;
            if ((cr->MeshMotion == LAGRANGIAN || cr->MeshMotion == DYNAMIC_LAGRANGIAN ||
                 cr->MeshMotion == TOTAL_ALE)) {
              /* add swelling as a source of volume */
              if (!mp->PorousMediaType) {
                source = -phi_i * (d_h3detJ_dmesh_bj)*wt * (1. - initial_volsolvent) /
                         (1. - volsolvent) * pd->etm[pg->imtrx][eqn][LOG2_SOURCE];
              }
            }

            /* add Pressure-Stabilized Petrov-Galerkin term
             * if desired.
             */
            pressure_stabilization = 0.;
            if (PSPG) {
              for (a = 0; a < WIM; a++) {
                meqn = R_PMOMENTUM1 + a;
                if (pd->e[pg->imtrx][meqn]) {

                  advection_a = 0.;
                  if (pd->e[pg->imtrx][meqn] & T_ADVECTION) {
                    if (pd->TimeIntegration != STEADY) {
                      advection_a = -rho / por2 * (1. + 2. * tt) * phi_j / dt * grad_pv[b][a] *
                                    grad_phi[i][a] * wt * tau_pspg * h3 * det_J *
                                    pd->etm[pg->imtrx][meqn][(LOG2_ADVECTION)];
                    }

                    for (p = 0; p < WIM; p++) {
                      advection_a += rho / por2 * (v[p] - x_dot[p]) *
                                     fv->d_grad_pv_dmesh[p][a][b][j] * grad_phi[i][a] * wt *
                                     tau_pspg * h3 * det_J *
                                     pd->etm[pg->imtrx][meqn][(LOG2_ADVECTION)];
                    }
                  }

                  diffusion = 0.;
                  if (pd->e[pg->imtrx][meqn] & T_DIFFUSION) {
                    diffusion -=
                        d_mu->X[b][j] * div_G[a] * pd->etm[pg->imtrx][meqn][(LOG2_DIFFUSION)];
                  }

                  pressure = 0.;
                  var1 = PRESSURE;
                  if (pd->v[pg->imtrx][var1]) {
                    if (pd->e[pg->imtrx][meqn] & T_DIFFUSION) {
                      pressure = fv->d_grad_P_dmesh[a][b][j] * grad_phi[i][a] * wt * tau_pspg * h3 *
                                 det_J * pd->etm[pg->imtrx][meqn][(LOG2_DIFFUSION)];
                    }
                  }

                  velocity_gradient = 0.;
                  var1 = VELOCITY_GRADIENT11;
                  if (pd->v[pg->imtrx][var1]) {
                    if (pd->e[pg->imtrx][meqn] & T_DIFFUSION) {
                      velocity_gradient -= fv->d_div_G_dmesh[a][b][j] * grad_phi[i][a] * wt *
                                           tau_pspg * h3 * det_J *
                                           pd->etm[pg->imtrx][meqn][(LOG2_DIFFUSION)];
                    }
                  }

                  stress = 0.;
                  var1 = POLYMER_STRESS11;
                  if (pd->v[pg->imtrx][var1]) {
                    if (pd->e[pg->imtrx][meqn] & T_DIFFUSION) {
                      stress -= fv->d_div_S_dmesh[0][a][b][j] * grad_phi[i][a] * wt * tau_pspg *
                                h3 * det_J * pd->etm[pg->imtrx][meqn][(LOG2_DIFFUSION)];
                    }
                  }

                  source_a = 0.;
                  if (pd->e[pg->imtrx][meqn] & T_SOURCE) {
                    source_a -= dfdX[a][b][j] * grad_phi[i][a] * wt * tau_pspg * h3 * det_J *
                                pd->etm[pg->imtrx][meqn][(LOG2_SOURCE)];
                  }

                  stabilization_a =
                      momentum_residual[a] * grad_phi[i][a] * wt * tau_pspg * d_h3detJ_dmesh_bj;

                  stabilization_b = momentum_residual[a] * bf[eqn]->d_grad_phi_dmesh[i][a][b][j] *
                                    wt * tau_pspg * h3 * det_J;

                  pressure_stabilization += advection_a + source_a + diffusion + pressure + stress +
                                            velocity_gradient + stabilization_a + stabilization_b;
                }
              }
            }

            lec->J[LEC_J_INDEX(peqn, pvar, i, j)] +=
                advection + source + pressure_stabilization + h_flux;
          }
        }
      }

      /*
       * J_c_d_rs
       */

      for (b = 0; b < dim; b++) {
        var = SOLID_DISPLACEMENT1 + b;
        if (pd->v[pg->imtrx][var]) {
          pvar = upd->vp[pg->imtrx][var];

          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];

            advection = 0.;

            if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
              if (cr->MeshMotion == TOTAL_ALE) {
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
      if (pd->v[pg->imtrx][var]) {
        pvar = upd->vp[pg->imtrx][var];

        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          phi_j = bf[var]->phi[j];

          for (w = 0; w < pd->Num_Species_Eqn; w++) {
            /* add swelling as a source of volume */
            source = 0.;
            if (!mp->PorousMediaType &&
                ((cr->MeshMotion == LAGRANGIAN || cr->MeshMotion == DYNAMIC_LAGRANGIAN))) {
              source = -phi_j * phi_i * h3 * det_J * wt * (1. - initial_volsolvent) /
                       (1. - volsolvent) / (1. - volsolvent) * pd->etm[pg->imtrx][eqn][LOG2_SOURCE];
            }

            pressure_stabilization = 0.;
            if (PSPG) {
              mass = 0.;
              advection = 0.;
              porous = 0.;
              source_a = 0.;
              diffusion = 0.;
              for (a = 0; a < WIM; a++) {
                meqn = R_PMOMENTUM1 + a;

                if (pd->TimeIntegration != STEADY) {
                  if (pd->e[pg->imtrx][meqn] & T_MASS) {
                    mass = d_rho->C[w][j] / por * pv_dot[a] * grad_phi[i][a] *
                           pd->etm[pg->imtrx][meqn][(LOG2_MASS)];
                  }
                }

                if (pd->e[pg->imtrx][meqn] & T_ADVECTION) {
                  advection = 0.;
                  for (p = 0; p < WIM; p++) {
                    advection += (v[p] - x_dot[p]) * grad_pv[p][a];
                  }
                  advection *= grad_phi[i][a] * d_rho->C[w][j] / por2 *
                               pd->etm[pg->imtrx][meqn][(LOG2_ADVECTION)];
                }

                if (pd->e[pg->imtrx][meqn] & T_POROUS_BRINK) {
                  porous = v[a] * (d_rho->C[w][j] * sc * speed / sqrt(per)) * grad_phi[i][a] *
                           pd->etm[pg->imtrx][meqn][(LOG2_POROUS_BRINK)];
                }

                if (pd->e[pg->imtrx][meqn] & T_SOURCE) {
                  source_a =
                      -grad_phi[i][a] * dfdC[a][w][j] * pd->etm[pg->imtrx][meqn][(LOG2_SOURCE)];
                }

                if (pd->e[pg->imtrx][meqn] & T_DIFFUSION) {
                  diffusion = -grad_phi[i][a] * d_mu->C[w][j] * div_G[a] *
                              pd->etm[pg->imtrx][meqn][(LOG2_DIFFUSION)];
                }
                pressure_stabilization += mass + advection + porous + source_a + diffusion;
              }
              pressure_stabilization *= tau_pspg * h3 * det_J * wt;
            }

            lec->J[LEC_J_INDEX(peqn, MAX_PROB_VAR + w, i, j)] += source + pressure_stabilization;
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
    }
  }

  return (status);
} /* end of function MMH_assemble_continuity */
/**********************************************************************/

int pmomentum_source_term(dbl f[DIM],                   /* Body force */
                          dbl dfdT[DIM][MDE],           /* For temperature dependence */
                          dbl dfdX[DIM][DIM][MDE],      /* For spatial dependence */
                          dbl dfdC[DIM][MAX_CONC][MDE], /* For concentration dependence */
                          dbl dfdv[DIM][DIM][MDE])      /* For velocity dependence */
{
  int err;
  int status = 0;

  /* initialize everything to zero */

  memset(f, 0, sizeof(double) * DIM);

  memset(dfdT, 0, sizeof(double) * DIM * MDE);

  memset(dfdX, 0, sizeof(double) * DIM * DIM * MDE);

  memset(dfdv, 0, sizeof(double) * DIM * DIM * MDE);

  memset(dfdC, 0, sizeof(double) * DIM * MAX_CONC * MDE);

  if (mp->MomentumSourceModel == SUSPENSION_PM) {
    /*
     * There were a couple of "extra" arguments that the current prototype
     * does not take. Hope they weren't critical.
     *    , 0, FALSE);
     */

    err = suspension_pm_particle_momentum_source(f, dfdT, dfdX, dfdC, dfdv);

    GOMA_EH(err, "Problems in suspension_pm_particle_momentum_source");
  } else {
    GOMA_EH(GOMA_ERROR, "No such Navier-Stokes Model");
  }
  return (status);
} /* end of function pmomentum_source_term */
/**********************************************************************/
/* end of file mm_fill_pthings.c */
/**********************************************************************/
