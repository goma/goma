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
 *$Id: mm_std_models.c,v 5.31 2010-07-30 20:48:38 prschun Exp $
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* GOMA include files */

#include "el_elm.h"
#include "el_geom.h"
#include "mm_as.h"
#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "mm_eh.h"
#include "mm_fill_ls.h"
#include "mm_fill_species.h"
#include "mm_fill_stress.h"
#include "mm_fill_terms.h"
#include "mm_fill_util.h"
#include "mm_mp.h"
#include "mm_mp_const.h"
#include "mm_mp_structs.h"
#include "mm_qtensor_model.h"
#include "mm_viscosity.h"
#include "rf_bc_const.h"
#include "rf_fem.h"
#include "rf_fem_const.h"
#include "rf_io_const.h"
#include "std.h"
#include "user_mp.h"

#define GOMA_MM_STD_MODELS_C
#include "density.h"
#include "mm_std_models.h"

/*********** R O U T I N E S   I N   T H I S   F I L E ************************
 *
 *       NAME            TYPE            CALLED_BY
 *    ------------             ---------               --------------
 * hydro_flux
 * suspension_balance
 * epoxy_dea_species_source
 * epoxy_species_source
 * foam_epoxy_species_source
 * epoxy_heat_source
 * bouss_momentum_source   int          assemble_momentum
 * bouss_and_jxb           int          assemble_momentum
 * joule_heat_source       int          assemble_energy
 * visc_diss_heat_source       int          assemble_energy
 * calc_KOH_Si_etch_rate_100  double
 ******************************************************************************/
/*
 * This file contains all implemented models for material properties and
 * source terms, and may eventually include constitutive equations.  Examples
 * are the Boussinesq model for momentum source, Joule heating and viscous
 * dissipation for heat source, etc.  As user-defined subroutines are built in
 * user_mp.c that are fairly standard and widely used they should be migrated
 * to this file and altered appropriately.
 ******************************************************************************/
/*
 * Boussinesq Momentum Source Model
 */

/*
 * int bouss_momentum_source (f, df, jxb, hydrostatic)
 *
 * ----------------------------------------------------------------------------
 * This routine is responsible for filling up the following forces and sensitivities
 * at the current gauss point:
 *     input:
 *
 *     output:  f[a]          - body force in direction [a]
 *              df            - dependence of body force, where:
 *
 *              df->T[a][j]    - derivative wrt temperature at node j.
 *              df->C[a][i][j] - derivative wrt mass frac species i at node j
 *              df->V[a][i][j] - derivative wrt velocity component i at node j
 *              df->X[a][0][j] - derivative wrt mesh disp. comp i at node j
 *
 *     input:  jxb = 1  add JXB lorentz term to momentum source
 *             jxb = 0  don't add JXB lorentz term
 *             N.B. Right now J field and B field have to come in as external
 *                  nodal fields using the "External Nodal Field" cards.
 *
 *	       hydrostatic = TRUE  -- add in the hydrostatic component so that
 *				      the entire body force term appears as
 *
 *				         rho * g * ( 1 - beta * (T - T_ref) )
 *
 *			     FALSE -- subtract off the hydrostatic pressure head
 *				      so that the body force term appears as
 *
 *				         rho * g * ( -beta * ( T - T_ref ) )
 *
 *   NB: The user need only supply f, dfdT, dfdC, etc.
 *       The mp struct is loaded up for you.
 *
 *       The FALSE hydrostatic option is often useful if buoyant forces are
 *	 weak and could be overwhelmed by the hydrostatic pressure field.
 *       It is currently bound to the BOUSSINESQ card, while the
 *       TRUE option is invoked by the BOUSS card to provide backward
 *	 compatibility (some example problems were performed with a hydrostatic
 *       mode incorporated into the pressure field.)
 *
 * ----------------------------------------------------------------------------
 */

int bouss_momentum_source(dbl f[DIM], /* Body force. */
                          MOMENTUM_SOURCE_DEPENDENCE_STRUCT *df,
                          int jxb,         /* Flag for turning on and off jxb */
                          int hydrostatic) /* Boolean for including hydrostatic *
                                            * head in thermal buoyancy term */
{
  int eqn, var;
  int a, b, c;

  int w;

  dbl T, C[MAX_CONC]; /* Convenient local variables */

  dbl J[DIM], B[DIM];     /* Local versions for REAL j-field and b-field */
  dbl J_I[DIM], B_I[DIM]; /* Local versions for Imag j-field and b-field */
  dbl g[DIM];

  int i, j;

  /* Begin Execution */

  /**********************************************************/

  /***********Load up convenient local variables*************/
  /*NB This ought to be done once for all fields at gauss pt*/

  T = fv->T;
  for (i = 0; i < pd->Num_Species_Eqn; i++)
    C[i] = fv->c[i];

  /* Next, if jxb is to be computed, load up J and B vectors */
  if (jxb) {
    /* First zero out */
    for (a = 0; a < 3; a++) {
      B[a] = J[a] = 0.;
      B_I[a] = J_I[a] = 0.;
    }

    for (a = 0; a < efv->Num_external_field; a++) {
      if (!strcmp(efv->name[a], "JX_REAL")) {
        if (pd->CoordinateSystem == CYLINDRICAL || pd->CoordinateSystem == SWIRLING ||
            pd->CoordinateSystem == PROJECTED_CARTESIAN) {
          J[1] = fv->external_field[a];

        } else {
          J[0] = fv->external_field[a];
        }
      }
      if (!strcmp(efv->name[a], "JX_IMAG")) {
        if (pd->CoordinateSystem == CYLINDRICAL || pd->CoordinateSystem == SWIRLING ||
            pd->CoordinateSystem == PROJECTED_CARTESIAN) {
          J_I[1] = fv->external_field[a];
        } else {
          J_I[0] = fv->external_field[a];
        }
      }
      if (!strcmp(efv->name[a], "JY_REAL")) {
        if (pd->CoordinateSystem == CYLINDRICAL || pd->CoordinateSystem == SWIRLING ||
            pd->CoordinateSystem == PROJECTED_CARTESIAN) {
          J[0] = -fv->external_field[a] / sqrt(2.);
        } else {
          J[1] = fv->external_field[a];
        }
      }
      if (!strcmp(efv->name[a], "JY_IMAG")) {
        if (pd->CoordinateSystem == CYLINDRICAL || pd->CoordinateSystem == SWIRLING ||
            pd->CoordinateSystem == PROJECTED_CARTESIAN) {
          J_I[0] = -fv->external_field[a];
        } else {
          J_I[1] = fv->external_field[a];
        }
      }
      if (!strcmp(efv->name[a], "JZ_REAL")) {
        J[2] = fv->external_field[a];
      }
      if (!strcmp(efv->name[a], "JZ_IMAG")) {
        J_I[2] = fv->external_field[a];
      }

      if (!strcmp(efv->name[a], "BX_REAL")) {
        if (pd->CoordinateSystem == CYLINDRICAL || pd->CoordinateSystem == SWIRLING ||
            pd->CoordinateSystem == PROJECTED_CARTESIAN) {
          B[1] = fv->external_field[a];
        } else {
          B[0] = fv->external_field[a];
        }
      }
      if (!strcmp(efv->name[a], "BX_IMAG")) {
        if (pd->CoordinateSystem == CYLINDRICAL || pd->CoordinateSystem == SWIRLING ||
            pd->CoordinateSystem == PROJECTED_CARTESIAN) {
          B_I[1] = fv->external_field[a];
        } else {
          B_I[0] = fv->external_field[a];
        }
      }
      if (!strcmp(efv->name[a], "BY_REAL")) {
        if (pd->CoordinateSystem == CYLINDRICAL || pd->CoordinateSystem == SWIRLING ||
            pd->CoordinateSystem == PROJECTED_CARTESIAN) {
          B[0] = -fv->external_field[a];
        } else {
          B[1] = fv->external_field[a];
        }
      }
      if (!strcmp(efv->name[a], "BY_IMAG")) {
        if (pd->CoordinateSystem == CYLINDRICAL || pd->CoordinateSystem == SWIRLING ||
            pd->CoordinateSystem == PROJECTED_CARTESIAN) {
          B_I[0] = -fv->external_field[a];
        } else {
          B_I[1] = fv->external_field[a];
        }
      }
      if (!strcmp(efv->name[a], "BZ_REAL")) {
        if (pd->CoordinateSystem == CYLINDRICAL || pd->CoordinateSystem == SWIRLING ||
            pd->CoordinateSystem == PROJECTED_CARTESIAN) {
          B[2] = -fv->external_field[a];
        } else {
          B[2] = fv->external_field[a];
        }
      }
      if (!strcmp(efv->name[a], "BZ_IMAG")) {
        if (pd->CoordinateSystem == CYLINDRICAL || pd->CoordinateSystem == SWIRLING ||
            pd->CoordinateSystem == PROJECTED_CARTESIAN) {
          B_I[2] = -fv->external_field[a];
        } else {
          B_I[2] = fv->external_field[a];
        }
      }
    }
  }

  /*components of gravity vector for input card */
  g[0] = mp->momentum_source[0];
  g[1] = mp->momentum_source[1];
  g[2] = mp->momentum_source[2];

  /**********************************************************/
  /* Temperature piece */
  if (pd->v[pg->imtrx][TEMPERATURE]) {
    for (a = 0; a < DIM; a++) {
      eqn = R_MOMENTUM1 + a;
      if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
        if (hydrostatic) {
          /*
           * This ought to be equivalent to the previous model.
           * All nonzero initialization is here now.
           */
          f[a] += mp->density * mp->momentum_source[a] *
                  (1 - mp->Volume_Expansion * (T - mp->reference[TEMPERATURE]));
        } else {
          /*
           * This definition implies that a hydrostatic pressure field
           * has been subtracted from the pressure field. Now you're
           * solving for P_prime!
           */
          f[a] += mp->density * mp->momentum_source[a] *
                  (-mp->Volume_Expansion * (T - mp->reference[TEMPERATURE]));
        }
      }
    }
  }

  /* Species piece */
  if (pd->v[pg->imtrx][MASS_FRACTION]) {
    for (a = 0; a < DIM; a++) {
      if (hydrostatic) {
        f[a] = g[a] * mp->density;
      }

      for (w = 0; w < pd->Num_Species_Eqn; w++) {
        eqn = R_MOMENTUM1 + a;
        if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
          f[a] +=
              -g[a] * mp->density * mp->species_vol_expansion[w] * (C[w] - mp->reference_concn[w]);
        } else if (pd->gv[R_LUBP]) {
          f[a] +=
              -g[a] * mp->density * mp->species_vol_expansion[w] * (C[w] - mp->reference_concn[w]);
        }
      }
    }
  }

  if (jxb) {
    /* If appropriate, jxb piece */
    /* NOTE Cludge  to stirring force only */
    for (a = 0; a < 3; a++) {
      eqn = R_MOMENTUM1 + a;
      for (b = 0; b < 3; b++) {
        for (c = 0; c < 3; c++)

        {
          if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
            /*
             * Use PRS e_ijk alternator permute macro in std.h
             * rather than eps, with hi name clash potential...
             */
            f[a] +=
                0.5 * (mp->u_momentum_source[0] * permute(a + 1, b + 1, c + 1) * J[b] * B[c] +
                       mp->u_momentum_source[0] * permute(a + 1, b + 1, c + 1) * J_I[b] * B_I[c]);
          }
        }
      }
    }
  }

  /* Now do sensitivies */

  var = TEMPERATURE;
  if (pd->v[pg->imtrx][TEMPERATURE]) {
    for (a = 0; a < DIM; a++) {
      eqn = R_MOMENTUM1 + a;
      if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          df->T[a][j] += -g[a] * mp->density * mp->Volume_Expansion * bf[var]->phi[j];
        }
      }
    }
  }

  if (pd->v[pg->imtrx][MASS_FRACTION]) {
    var = MASS_FRACTION;
    for (a = 0; a < DIM; a++) {
      eqn = R_MOMENTUM1 + a;
      if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
        for (w = 0; w < pd->Num_Species_Eqn; w++) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            df->C[a][w][j] += -g[a] * mp->density * mp->species_vol_expansion[w] * bf[var]->phi[j];
          }
        }
      } else if (pd->e[pg->imtrx][R_LUBP]) {
        for (w = 0; w < pd->Num_Species_Eqn; w++) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            df->C[a][w][j] += -g[a] * mp->density * mp->species_vol_expansion[w] * bf[var]->phi[j];
          }
        }
      }
    }
  }

  return (0); /* failsafe */
}

/*
 * EHD Polarization Source Model
 */

/*
 * int EHD_POLARIZATION_source (f, df)
 *
 * ----------------------------------------------------------------------------
 * This routine is responsible for filling up the following forces and sensitivities
 * at the current gauss point:
 *     intput:
 *
 *     output:  f[a]          - body force in direction [a]
 *              df            - dependence of body force, where:
 *
 *              df=>E[a][1][j] - derivative wrt electric field comp 1 at node j.
 *              df->X[a][0][j] - derivative wrt mesh disp. comp i at node j
 *
 * ----------------------------------------------------------------------------
 */

int EHD_POLARIZATION_source(dbl f[DIM], /* Body force. */
                            MOMENTUM_SOURCE_DEPENDENCE_STRUCT *df) {
  int eqn, var;
  int dim;
  int a, b, p;
  dbl phi_j;
  dbl advection_a;
  dbl Efield[DIM]; /* Local versions for E-field    */
  dbl g[DIM];
  int j;
  dim = pd->Num_Dim;

  /**********************************************************/

  /***********Load up convenient local variables*************/
  /*NB This ought to be done once for all fields at gauss pt*/

  for (a = 0; a < DIM; a++)
    Efield[a] = fv->E_field[a];

  /*components of whatever (charge density?) from input card */
  g[0] = mp->momentum_source[0];

  /**********************************************************/
  for (a = 0; a < dim; a++) {
    eqn = R_MOMENTUM1 + a;

    for (b = 0; b < dim; b++) {
      if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
        f[a] += g[0] * Efield[b] * fv->grad_E_field[b][a];
      }
    }
  }
  /* Now do sensitivies */

  var = EFIELD1;
  if (pd->v[pg->imtrx][var]) {
    for (a = 0; a < dim; a++) {
      eqn = R_MOMENTUM1 + a;
      for (b = 0; b < dim; b++) {
        var = R_EFIELD1 + b;
        if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];
            df->E[a][b][j] += g[0] * phi_j * fv->grad_E_field[b][a];

            for (p = 0; p < dim; p++) {
              df->E[a][b][j] += g[0] * Efield[p] * bf[var]->grad_phi_e[j][b][p][a];
            }
          }
        }
      }
    }
  }

  var = EFIELD1;
  if (pd->v[pg->imtrx][MESH_DISPLACEMENT1]) {
    for (a = 0; a < dim; a++) {
      eqn = R_MOMENTUM1 + a;

      for (b = 0; b < dim; b++) {
        var = R_EFIELD1 + b;
        if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            advection_a = 0.;
            for (p = 0; p < dim; p++) {
              advection_a += Efield[p] * fv->d_grad_E_field_dmesh[p][a][b][j];
            }

            df->X[a][b][j] += advection_a;
          }
        }
      }
    }
  }

  return (0); /* failsafe */
}

/*
 * EHD Polarization Momentum Source Model
 */

/*
 * int gravity_vibrational_source (f, df, time)
 * ----------------------------------------------------------------------------
 * This routine is responsible for filling up the following forces and sensitivities
 * at the current gauss point for a problem with gravity plus a cosine vibrational term
 * in the gravity direction. We will assume that rho is constant for this first
 * implementation:
 *
 *     f=rho*g+rho*(omega^2*A*cos(omega*time))
 *
 *     input:
 *
 *     output:  f[a]          - body force in direction [a]
 *              df            - dependence of body force, where there are no Jacobian terms:
 * ----------------------------------------------------------------------------
 */

int gravity_vibrational_source(dbl f[DIM], /* Body force. */
                               MOMENTUM_SOURCE_DEPENDENCE_STRUCT *df,
                               dbl time) {
  int eqn;
  int dim;
  int a;
  dbl g[DIM];
  dbl omega;  /* frequency of vibration  */
  dbl omega2; /* frequency of vibration squared */
  dbl A;      /* amplitude of vibration  */
  dbl g_mag;  /* magnitude of gravity vector for direction of vibration */

  DENSITY_DEPENDENCE_STRUCT d_rho_struct; /* density dependence */
  DENSITY_DEPENDENCE_STRUCT *d_rho = &d_rho_struct;
  double rho = density(d_rho, time);

  dim = pd->Num_Dim;

  /***********Load up convenient local variables*************/
  /*NB This ought to be done once for all fields at gauss pt*/

  /*components of gravity vector from input card */
  for (a = 0; a < DIM; a++)
    g[a] = mp->momentum_source[a];
  A = mp->u_momentum_source[0];
  omega = mp->u_momentum_source[1];
  omega2 = omega * omega;

  g_mag = 0.;
  for (a = 0; a < DIM; a++) {
    g_mag += g[a] * g[a];
  }
  g_mag = sqrt(g_mag);
  if (g_mag == 0.) {
    GOMA_EH(GOMA_ERROR,
            "Problems in GRAV_VIBRATIONAL force routine - must have nonzero gravity vector ");
  }

  /**********************************************************/
  for (a = 0; a < dim; a++) {
    eqn = R_MOMENTUM1 + a;

    if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
      f[a] += rho * (g[a] + omega2 * A * cos(omega * time) * g[a] / g_mag);
    }
  }
  /* Now do sensitivies - currently for constant density, there are no sensitivities */

  return (0); /* failsafe */
}

/*
 * This is the source term for a suspension model
 * where the fluid and particles have different
 * densities
 */
int suspend_momentum_source(dbl f[DIM], /* Body force. */
                            MOMENTUM_SOURCE_DEPENDENCE_STRUCT *df) {
  int eqn, var;
  int a;

  int species;

  dbl delta_rho_g[DIM];
  dbl del_rho;
  dbl vol;

  int j;
  int status = 1;

  /* Begin Execution */

  /* difference between solid phase density and fluid phase density */
  if (mp->DensityModel == SUSPENSION) {
    species = (int)mp->u_density[0];
    del_rho = mp->u_density[2] - mp->u_density[1];
  } else {
    GOMA_EH(GOMA_ERROR, "No suspension momentum source without suspension density model.");
    return (status);
  }

  /*  if(fv->c[species] > 0.)
    {
      C[species] = fv->c[species];
    }
  else
    {
      C[species] = 0.;
    }  */

  /*  vol = C[species] -  mp->u_momentum_source[0];  */

  vol = fv->c[species] - mp->u_momentum_source[0];

  /*components of gravity vector for input card */
  delta_rho_g[0] = mp->momentum_source[0] * del_rho;
  delta_rho_g[1] = mp->momentum_source[1] * del_rho;
  delta_rho_g[2] = mp->momentum_source[2] * del_rho;

  /**********************************************************/

  /* Species piece */
  if (pd->v[pg->imtrx][MASS_FRACTION]) {
    for (a = 0; a < DIM; a++)

    {
      eqn = R_MOMENTUM1 + a;
      if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
        f[a] += delta_rho_g[a] * vol;
      }
    }
  }

  if (pd->v[pg->imtrx][MASS_FRACTION]) {
    var = MASS_FRACTION;
    for (a = 0; a < DIM; a++) {
      eqn = R_MOMENTUM1 + a;
      if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          df->C[a][species][j] += delta_rho_g[a] * bf[var]->phi[j];
        }
      }
    }
  }
  return (status);
}

/*
 * This the momentum source term for problems
 * involving filling of one fluid into another
 * The density is a function of degree of filling
 * and momentum source is product of this density
 * and the gravity vector.
 */
int fill_momentum_source(double f[DIM]) {
  /* Local Variables */

  double rho;
  double *param = mp->u_density;

  if (mp->DensityModel != FILL) {
    GOMA_EH(GOMA_ERROR, "No fill momentum source without fill density model.");
    return (0);
  }

  if (ls != NULL) {
    GOMA_EH(GOMA_ERROR, "Use LEVEL_SET instead of FILL density model.");
    return (0);
  }

  rho = param[0] * fv->F + param[1] * (1.0 - fv->F);

  /* Now fill in momentum source vector using gravity vector
   * stored in mp->momentum_source .  Note that since FILL
   * is an explicit variable if makes no contributions to the
   * Jacobian.
   */

  if (pd->v[pg->imtrx][FILL]) {
    f[0] = mp->momentum_source[0] * rho;
    f[1] = mp->momentum_source[1] * rho;
    f[2] = mp->momentum_source[2] * rho;
  }

  return (1);
}

/*
 * This is the reaction source term for the species
 * equation using an extent-of-reaction model,
 * specifically for epoxy828/DEA reaction kinetics
 */

int epoxy_dea_species_source(int species_no, /* Current species number */
                             double *param)  /* pointer to user-defined parameter list */

{
  int eqn, var, var_offset, imtrx;

  dbl T; /* temperature for rate constants */
  dbl A1, E1, A2, E2, A3;
  dbl m, n;
  dbl k1, k2;
  dbl alpha, alpha_m, alpha_m1, alpha_n, alpha_n1;

  /* Begin Execution */

  if (pd->gv[TEMPERATURE]) {
    T = fv->T;
  } else {
    T = upd->Process_Temperature;
  }
  /* extent of reaction, alpha */
  alpha = fv->c[species_no];
  /*  if(alpha <= 0.) alpha = 0.0001; */
  A1 = param[0];
  E1 = param[1];
  A2 = param[2];
  E2 = param[3];
  A3 = param[4];
  n = 1.6;
  m = 2.2;

  k1 = A1 * exp(-E1 / T);
  k2 = A2 * exp(-E2 / T);

  /* Three rate expressions for Epoxy-DEA system,
     one for T<65C, one between 65 and 90C, and one for T>90C.
  */
  if (T <= 0.) {
    k1 = A1;
    k2 = A2;
  } else if (T > 338.15 && T < 363.15) {
    k2 = A3 * (90. - (T - 273.15)) * pow((T - 273.15), -6.0);
    m = 74. * k2 * 60.;
  } else if (T >= 363.15) {
    m = 0.;
    k2 = 0.;
  }

  if (alpha > 0.0) {
    alpha_m = pow(alpha, m);
    alpha_m1 = pow(alpha, m - 1);
  } else {
    alpha_m = 0.;
    alpha_m1 = 0.;
  }

  alpha_n = pow((1. - alpha), n);
  alpha_n1 = pow((1. - alpha), n - 1);

  /**********************************************************/

  /* Species piece */
  eqn = MASS_FRACTION;
  for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
    if (pd->e[imtrx][eqn] & T_SOURCE) {
      mp->species_source[species_no] = (k1 + k2 * alpha_m) * alpha_n;

      /* Jacobian entries for source term */
      var = MASS_FRACTION;
      if (pd->v[pg->imtrx][var]) {
        var_offset = MAX_VARIABLE_TYPES + species_no;
        mp->d_species_source[var_offset] =
            (m * k2 * alpha_m1) * alpha_n - (k1 + k2 * alpha_m) * n * alpha_n1;
      }

      var = TEMPERATURE;
      if (pd->v[pg->imtrx][var]) {
        if (T <= 0.) {
          mp->d_species_source[var] = 0.;
        } else if (T > 338.15 && T < 363.15) {
          if (alpha > 0) {
            mp->d_species_source[var] =
                (k1 * E1 / (T * T) +
                 (A2 * (5. * pow((T - 273.15), -6.0) - 540. * pow((T - 273.15), -7.0))) *
                     (alpha_m * log(alpha_m) + alpha_m)) *
                alpha_n;
          } else {
            mp->d_species_source[var] = (k1 * E1 + k2 * alpha_m) * alpha_n / (T * T);
          }
        } else {
          mp->d_species_source[var] = (k1 * E1 + k2 * E2 * alpha_m) * alpha_n / (T * T);
        }
      }
    }
  }
  return 0;
}
/*****************************************************************************/
/* END of routine epoxy_dea_species_source */
/*****************************************************************************/

/*
 * This is the reaction source term for the species
 * equation using an extent-of-reaction model
 */

int epoxy_species_source(int species_no, /* Current species number */
                         double *param)  /* pointer to user-defined parameter list */

{
  /* Local Variables */
  int eqn, var, var_offset, imtrx;
  /*  int p, q, a, b, c;*/

  /*  int v,w;*/

  /*  int mdofs,vdofs;*/

  /*  dbl C[MAX_CONC]; Convenient local variables */
  dbl T; /* temperature for rate constants */
  dbl A1, E1, A2, E2;
  dbl m, n;
  dbl k1, k2;
  dbl alpha, alpha_m, alpha_m1, alpha_n, alpha_n1;

  /* Begin Execution */

  if (pd->gv[TEMPERATURE]) {
    T = fv->T;
  } else {
    T = upd->Process_Temperature;
  }

  /* extent of reaction, alpha */
  alpha = fv->c[species_no];
  /*  if(alpha <= 0.) alpha = 0.0001; */
  A1 = param[0];
  E1 = param[1];
  A2 = param[2];
  E2 = param[3];
  m = param[4];
  n = param[5];

  /* Arhenius type rate constants for extent of reaction model */
  k1 = A1 * exp(-E1 / T);
  k2 = A2 * exp(-E2 / T);

  if (alpha > 0.0) {
    alpha_m = pow(alpha, m);
    alpha_m1 = pow(alpha, m - 1);
  } else {
    alpha_m = 0.;
    alpha_m1 = 0.;
  }

  alpha_n = pow((1. - alpha), n);
  alpha_n1 = pow((1. - alpha), n - 1);

  /**********************************************************/

  /* Species piece */
  eqn = MASS_FRACTION;
  for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
    if (pd->e[imtrx][eqn] & T_SOURCE) {
      mp->species_source[species_no] = (k1 + k2 * alpha_m) * alpha_n;

      /* Jacobian entries for source term */
      var = MASS_FRACTION;
      if (pd->v[pg->imtrx][var]) {
        var_offset = MAX_VARIABLE_TYPES + species_no;
        mp->d_species_source[var_offset] =
            (m * k2 * alpha_m1) * alpha_n - (k1 + k2 * alpha_m) * n * alpha_n1;
      }

      var = TEMPERATURE;
      if (pd->v[pg->imtrx][var]) {
        mp->d_species_source[var] = (k1 * E1 + k2 * E2 * alpha_m) * alpha_n / (T * T);
      }
    }
  }
  return 0;
}
/*****************************************************************************/
/* END of routine epoxy_species_source */
/*****************************************************************************/

int bond_species_source(int species_no, /* Current species number */
                        double *param)  /* pointer to user-defined parameter list */

{
  /* Local Variables */
  int eqn, var, var_offset;
  /*  int p, q, a, b, c;*/

  /*  int v,w;*/

  /*  int mdofs,vdofs;*/

  /*  dbl C[MAX_CONC]; Convenient local variables */
  dbl n0, aexp, bexp, nn;
  dbl k1, k2;
  dbl shear;
  dbl gterm_a, gterm_b;
  dbl d_gterm_a, d_gterm_b;
  dbl offset = 0.00001;

  /* Begin Execution */

  /* structure factor, nn */
  nn = fv->c[species_no];
  shear = fv->SH;
  k1 = param[0];
  k2 = param[1];
  n0 = param[2];
  aexp = param[3];
  bexp = param[4];

  shear = fv->SH;
  if (shear >= .0) {
    gterm_a = pow(shear + offset, aexp);
    gterm_b = pow(shear + offset, bexp);
    d_gterm_a = aexp * pow(shear + offset, aexp - 1.);
    d_gterm_b = bexp * pow(shear + offset, bexp - 1.);
  } else {
    gterm_a = 0.;
    gterm_b = 0.;
    d_gterm_a = 0.;
    d_gterm_b = 0.;
  }

  /* clip negative values */
  if (nn < 1.e-5)
    nn = 1.e-5;

  /**********************************************************/

  /* Species piece */
  eqn = MASS_FRACTION;
  if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
    mp->species_source[species_no] = -k1 * nn * gterm_a + k2 * (n0 - nn) * gterm_b;

    /* Jacobian entries for source term */
    var = MASS_FRACTION;
    if (pd->v[pg->imtrx][var]) {
      var_offset = MAX_VARIABLE_TYPES + species_no;
      mp->d_species_source[var_offset] = -k1 * gterm_a - k2 * gterm_b;
    }

    var = SHEAR_RATE;
    if (pd->v[pg->imtrx][var]) {
      mp->d_species_source[var] = -k1 * nn * d_gterm_a + k2 * (n0 - nn) * d_gterm_b;
    }
  }
  return 0;
}
/*****************************************************************************/
/* END of routine bond_species_source */
/*****************************************************************************/

/*
 * foam_epoxy_species_source:
 *
 * This is the foam epoxy source term for the species
 * equation using an equilibrium model
 *
 *  The source term is
 *
 *       S = aT * (C - C0)   for T > Tb
 *       S = 0               for T < Tb
 *
 *  I added the Tb - 0.001 term to eliminate jacobian errors that occur when
 *  T = Tb is a startup condition.
 *
 *  C0 was added because that emulates boiling much better. C0 should also
 *  be a function of temperature. but, that hasn't been added yet.
 */
int foam_epoxy_species_source(int species_no, /* Current species number */
                              double *param,
                              double tt,
                              double dt)
/* param - pointer to user-defined parameter list */
/* tt, dt - time derivative parameters */
{
  int eqn, var, imtrx;
  /*  int p, q, a, b, c;*/
  double Press, rho, rho2, rho_v_inv, rho_v, d_rho_v_dT, d_rho_v_inv_dT;
  double rho_a_inv, d_rho_a_inv_dT;
  double drho_c_v, drho_c_a, drho_c_l, drho_T;
  double aT, bT, p_vap, vch, Cc, Ce, ff_c, ff_e, sigma;
  double Rc_1 = 0.0, Rc_2 = 0.0, Rc = 0.0, dRc_dc_v = 0.0, dRc_dc_a = 0.0, dRc_dc_l = 0.0,
         dRc_dT = 0.0;
  double Re_1 = 0.0, Re_2 = 0.0, Re = 0.0, dRe_dc_v = 0.0, dRe_dc_a = 0.0, dRe_dc_l = 0.0,
         dRe_dT = 0.0;
  double Rgas = 0.0, MW_f = 0.0, MW_a = 0.0, rho_epoxy = 0.0, rho_fluor = 0.0, T = 0.0;
  int species_l = -1, species_v = -1, species_a = -1;

  if (mp->DensityModel == DENSITY_FOAM_CONC) {
    species_l = (int)mp->u_density[0]; /* species number fluorinert liquid */
    species_v = (int)mp->u_density[1]; /* species number fluorinert vapor */
    species_a = (int)mp->u_density[2]; /* species number air vapor */
    Rgas = mp->u_density[3];           /* Gas constant in appropriate units */
    MW_f = mp->u_density[4];           /* molecular weight fluorinert */
    MW_a = mp->u_density[5];           /* molecular weight air */
    rho_epoxy = mp->u_density[6];      /* Density of epoxy resin */
    rho_fluor = mp->u_density[7];      /* Density of liquid fluorinert */
  } else {
    GOMA_EH(-1, "Foam_epoxy_species_source only works with foam_conc density model");
  }

  /* Begin Execution */

  if (pd->gv[TEMPERATURE]) {
    T = fv->T;
  } else {
    T = upd->Process_Temperature;
  }

  Press = upd->Pressure_Datum;
  rho_v_inv = Rgas * T / (Press * MW_f);
  rho_v = 1. / rho_v_inv;
  d_rho_v_dT = -(Press * MW_f) / (Rgas * T * T);
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

  aT = param[0];
  bT = param[1];

  /* vapor pressure is a function of temperature only */
  p_vap = aT - bT / T;
  /* characteristic velocity */
  vch = param[2];
  /* coefficient for condensation */
  Cc = param[3];
  /* coefficient for evaporation */
  Ce = param[4];

  /* grab surface tension from material properties */
  sigma = mp->surface_tension;

  /* Linear model for dx/dT */

  ff_c = Cc * vch / sigma;
  ff_e = Ce * vch / sigma;

  /* Use atmospheric pressure instead of the dynamic pressure to simplify things */
  if (Press > p_vap) {
    Rc_1 = ff_c * rho_fluor * fv->c[species_v] / rho;
    Rc_2 = pow(rho_fluor * (Press - p_vap) * 2. / 3., 0.5);
    Rc = Rc_1 * Rc_2;
    Re = 0.;

    dRc_dc_v = (ff_c * rho_fluor / rho) * (1. - fv->c[species_v] * drho_c_v / rho) * Rc_2;
    dRc_dc_a = -(ff_c * rho_fluor * fv->c[species_v] / rho2 * drho_c_a) * Rc_2;
    dRc_dc_l = -(ff_c * rho_fluor * fv->c[species_v] / rho2 * drho_c_l) * Rc_2;

    dRc_dT = (ff_c * rho_fluor * fv->c[species_v] / rho2 * drho_T) * Rc_2 -
             Rc_1 / (Rc_2 * 3.) * rho_fluor * bT;

    dRe_dc_v = 0.;
    dRe_dc_a = 0.;
    dRe_dc_l = 0.;

    dRe_dT = 0.;
  } else if (p_vap > Press) {
    Rc_1 = ff_c * rho_fluor * fv->c[species_v] / rho;
    Rc_2 = pow(rho_fluor * (Press - p_vap) * 2. / 3., 0.5);
    Rc = 0.;
    Re_1 = ff_e * rho_v * fv->c[species_l] / rho;
    Re_2 = pow(rho_fluor * (p_vap - Press) * 2. / 3., 0.5);
    Re = Re_1 * Re_2;

    dRc_dc_v = 0.;
    dRc_dc_a = 0.;
    dRc_dc_l = 0.;

    dRc_dT = 0.;

    dRe_dc_v = -(ff_e * rho_v * fv->c[species_l] / rho2 * drho_c_v) * Re_2;
    dRe_dc_a = -(ff_e * rho_v * fv->c[species_l] / rho2 * drho_c_a) * Re_2;
    dRe_dc_l = (ff_e * rho_v / rho - ff_e * rho_v * fv->c[species_l] / rho2 * drho_c_a) * Re_2;

    dRe_dT = (ff_e * d_rho_v_dT * fv->c[species_l] / rho -
              ff_e * rho_v * fv->c[species_l] / rho2 * drho_T) *
                 Re_2 +
             Rc_1 / (Rc_2 * 3.) * rho_fluor * bT;
  }

  /**********************************************************/

  /* Species piece */
  eqn = MASS_FRACTION;
  for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
    if (pd->e[imtrx][eqn] & T_SOURCE) {
      mp->species_source[species_no] = Rc - Re;

      /* Jacobian entries for source term */
      var = MASS_FRACTION;
      if (pd->v[pg->imtrx][var]) {
        mp->d_species_source[MAX_VARIABLE_TYPES + species_v] = dRc_dc_v - dRe_dc_v;
        mp->d_species_source[MAX_VARIABLE_TYPES + species_a] = dRc_dc_a - dRe_dc_a;
        mp->d_species_source[MAX_VARIABLE_TYPES + species_l] = dRc_dc_l - dRe_dc_l;
      }

      var = TEMPERATURE;
      if (pd->v[pg->imtrx][var]) {
        mp->d_species_source[var] = dRc_dT - dRe_dT;
      }
    }
  }
  return 0;
}

int foam_pmdi10_rxn_species_source(int species_no, /* Current species number */
                                   double *param,
                                   double tt,
                                   double dt)
/* param - pointer to user-defined parameter list */
/* tt, dt - time derivative parameters */
{
  int eqn, var;

  double k0 = param[0];
  double w_rxn = param[1];
  double beta = param[2];
  double C_1 = param[3];
  double C_2 = param[4];
  double m = param[5];
  double n = param[6];
  double b = param[7];
  double T_g0 = param[8];
  double T_ginf = param[9];
  double A = param[10];
  double E_norm = param[11];

  double T = fv->T;

  double source = 0;

  if (T <= 0) {
    source = 0;
    mp->d_species_source[MAX_VARIABLE_TYPES + species_no] = 0;
    mp->d_species_source[TEMPERATURE] = 0;
    return (source);
  }
  double xi = fv->c[species_no];

  double T_g = (T_g0 * (1 - xi) + A * xi * T_ginf) / (1 - xi + A * xi);
  double d_T_g_dT = 0;
  double d_T_g_dC = ((-(T_g0 + A * T_ginf) * (1 - xi + A * xi)) -
                     (T_g0 * (1 - xi) + A * xi * T_ginf) * (-1 + A)) /
                    ((1 - xi + A * xi) * (1 - xi + A * xi));

  double frac = -C_1 * (T - T_g) / (C_2 + T - T_g);
  double d_frac_dT =
      ((-C_1 * (1 - d_T_g_dT)) * (C_2 + T - T_g) - (-C_1 * (T - T_g)) * (1 - d_T_g_dT)) /
      ((C_2 + T - T_g) * (C_2 + T - T_g));
  double d_frac_dC = ((-C_1 * (-d_T_g_dC)) * (C_2 + T - T_g) - (-C_1 * (T - T_g)) * (d_T_g_dC)) /
                     ((C_2 + T - T_g) * (C_2 + T - T_g));

  double a_T = pow(10, frac);
  double d_a_T_dT = log(10) * a_T * d_frac_dT;
  double d_a_T_dC = log(10) * a_T * d_frac_dC;

  double k = (pow(1 + w_rxn * a_T, -beta)) * k0 * exp(-E_norm / T);

  double d_k_dT = w_rxn * d_a_T_dT * (-beta) * k / (1 + w_rxn * a_T) + E_norm / (fv->T * fv->T) * k;
  double d_k_dC = w_rxn * d_a_T_dC * (-beta) * k / (1 + w_rxn * a_T);

  if (xi < 0) {
    source = k * b;
  } else {
    source = k * (b + pow(xi, m)) * pow(1 - xi, n);
  }

  /**********************************************************/

  /* Species piece */
  eqn = MASS_FRACTION;
  if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
    mp->species_source[species_no] = source;

    /* Jacobian entries for source term */
    var = MASS_FRACTION;
    if (pd->v[pg->imtrx][var]) {
      if (xi < 0) {
        mp->d_species_source[MAX_VARIABLE_TYPES + species_no] = d_k_dC * b;
      } else {
        mp->d_species_source[MAX_VARIABLE_TYPES + species_no] =
            d_k_dC * (b + pow(xi, m)) * pow(1 - xi, n);
        mp->d_species_source[MAX_VARIABLE_TYPES + species_no] += n * source / (1 - xi);
        if (xi > 0) {
          mp->d_species_source[MAX_VARIABLE_TYPES + species_no] +=
              m * k * pow(1 - xi, n) * pow(xi, m) / xi;
        }
      }
    }

    var = TEMPERATURE;
    if (pd->v[pg->imtrx][var]) {
      if (xi < 0) {
        mp->d_species_source[var] = d_k_dT * (b + pow(xi, m)) * pow(1 - xi, n);
      } else {
        mp->d_species_source[var] = d_k_dT * b;
      }
    }
  }

  return 0;
}

int foam_pmdi10_h2o_species_source(int species_no, /* Current species number */
                                   double *param,
                                   double time,
                                   double tt,
                                   double dt)
/* param - pointer to user-defined parameter list */
/* tt, dt - time derivative parameters */
{
  int eqn, var;

  double CH2O = fv->c[species_no];
  double T = fv->T;
  double n = param[0];
  double t_nuc = param[1];
  double A = param[2];
  double norm_E = param[3];

  double N = 0.5 * (1 + tanh((time - t_nuc) / t_nuc));

  double source = 0;

  if (T <= 0) {
    source = 0;
    mp->d_species_source[MAX_VARIABLE_TYPES + species_no] = 0;
    mp->d_species_source[TEMPERATURE] = 0;
    return (source);
  }

  if (CH2O <= 0) {
    source = 0;
    mp->species_source[species_no] = 0;
  } else {
    source = -N * A * exp(-norm_E / T) * pow(CH2O, n);
  }
  /**********************************************************/

  /* Species piece */
  eqn = MASS_FRACTION;
  if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
    mp->species_source[species_no] = source;

    /* Jacobian entries for source term */
    var = MASS_FRACTION;
    if (pd->v[pg->imtrx][var]) {
      if (CH2O > 0) {
        mp->d_species_source[MAX_VARIABLE_TYPES + species_no] = source * n / CH2O;
      }
    }

    var = TEMPERATURE;
    if (pd->v[pg->imtrx][var]) {
      mp->d_species_source[var] = -norm_E / (T * T) * source;
    }
  }

  return source;
}

int foam_pmdi10_co2_species_source(int species_no, /* Current species number */
                                   double *param,
                                   double time,
                                   double tt,
                                   double dt) {
  int eqn, var;

  double T = fv->T;
  int wH2O = -1;
  int w;

  for (w = 0; w < pd->Num_Species; w++) {
    if (mp->SpeciesSourceModel[w] == FOAM_PMDI_10_H2O) {
      wH2O = w;
      break;
    }
  }

  if (wH2O == -1) {
    GOMA_EH(GOMA_ERROR, "Expected to find a speices with source FOAM_PMDI_10_H2O");
    return 0;
  }
  double CH2O = fv->c[wH2O];
  double n = mp->u_species_source[wH2O][0];
  double t_nuc = mp->u_species_source[wH2O][1];
  double A = mp->u_species_source[wH2O][2];
  double norm_E = mp->u_species_source[wH2O][3];

  double N = 0.5 * (1 + tanh((time - t_nuc) / t_nuc));

  double source;

  if (T <= 0) {
    source = 0;
    mp->d_species_source[MAX_VARIABLE_TYPES + species_no] = 0;
    mp->d_species_source[MAX_VARIABLE_TYPES + wH2O] = 0;
    mp->d_species_source[TEMPERATURE] = 0;
    return (source);
  }

  if (CH2O <= 0) {
    source = 0;
    mp->species_source[species_no] = 0;
  } else {
    source = N * A * exp(-norm_E / T) * pow(CH2O, n);
  }

  /**********************************************************/

  /* Species piece */
  eqn = MASS_FRACTION;
  if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
    mp->species_source[species_no] = source;

    /* Jacobian entries for source term */
    var = MASS_FRACTION;
    if (pd->v[pg->imtrx][var]) {
      if (CH2O > 0) {
        mp->d_species_source[MAX_VARIABLE_TYPES + wH2O] = source * n / CH2O;
      }
    }

    var = TEMPERATURE;
    if (pd->v[pg->imtrx][var]) {
      mp->d_species_source[var] = -norm_E / (T * T) * source;
    }
  }

  return source;
}

int foam_pmdi10_co2_liq_species_source(int species_no, /* Current species number */
                                       struct Species_Conservation_Terms *st,
                                       double *param,
                                       double time,
                                       double tt,
                                       double dt) {
  int eqn, var;

  double T = fv->T;
  int wH2O = -1;
  int w;
  struct moment_growth_rate *MGR;

  for (w = 0; w < pd->Num_Species; w++) {
    if (mp->SpeciesSourceModel[w] == FOAM_PMDI_10_H2O) {
      wH2O = w;
      break;
    }
  }

  if (!pd->gv[MOMENT1]) {
    GOMA_EH(GOMA_ERROR, "Expected to find moment equations for FOAM_PMDI_10_CO2_LIQ");
    return -1;
  }

  if (wH2O == -1) {
    GOMA_EH(GOMA_ERROR, "Expected to find a speices with source FOAM_PMDI_10_H2O");
    return -1;
  }

  MGR = calloc(sizeof(struct moment_growth_rate), 1);
  int err = get_moment_growth_rate_term(MGR);
  if (err) {
    free(MGR);
    return -1;
  }

  double Rgas_const = 8.31;
  double ref_press = 1e6;

  if (mp->DensityModel == DENSITY_FOAM_PMDI_10) {
    ref_press = mp->u_density[2];
    Rgas_const = mp->u_density[3];
  } else {
    GOMA_EH(GOMA_ERROR, "Expected DENSITY_FOAM_PMDI_10 density model");
  }

  double CH2O = fv->c[wH2O];
  double n = mp->u_species_source[wH2O][0];
  double t_nuc = mp->u_species_source[wH2O][1];
  double A = mp->u_species_source[wH2O][2];
  double norm_E = mp->u_species_source[wH2O][3];

  double N = 0.5 * (1 + tanh((time - t_nuc) / t_nuc));

  double source;

  if (T <= 0) {
    source = 0;
    mp->d_species_source[MAX_VARIABLE_TYPES + species_no] = 0;
    mp->d_species_source[MAX_VARIABLE_TYPES + wH2O] = 0;
    mp->d_species_source[TEMPERATURE] = 0;
    free(MGR);
    return (source);
  }

  double source_a;

  if (CH2O <= 0) {
    source = 0;
    source_a = 0;
    mp->species_source[species_no] = 0;
  } else {
    source_a = N * A * exp(-norm_E / T) * pow(CH2O, n);
  }
  source = source_a - MGR->G[species_no][1] * ref_press / (Rgas_const * T);

  /**********************************************************/

  /* Species piece */
  eqn = MASS_FRACTION;
  if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
    int j;
    mp->species_source[species_no] = source;
    st->MassSource[species_no] = source;

    /* Jacobian entries for source term */
    var = MASS_FRACTION;
    if (pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        if (CH2O > 0) {
          st->d_MassSource_dc[species_no][wH2O][j] = source_a * n / CH2O * bf[var]->phi[j];
        }
        st->d_MassSource_dc[species_no][species_no][j] =
            -MGR->d_G_dC[species_no][1][j] * ref_press / (Rgas_const * T);
      }
    }

    var = TEMPERATURE;
    if (pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        st->d_MassSource_dT[species_no][j] =
            -norm_E / (T * T) * source_a * bf[var]->phi[j] -
            MGR->d_G_dT[species_no][1][j] * ref_press / (Rgas_const * T);
      }
    }
  }
  free(MGR);
  return 0;
}

int foam_pmdi10_co2_gas_species_source(int species_no, /* Current species number */
                                       struct Species_Conservation_Terms *st,
                                       double *param,
                                       double time,
                                       double tt,
                                       double dt) {
  int eqn, var;

  double T = fv->T;
  int wH2O = -1;
  int wCO2Liq = -1;
  int w;
  struct moment_growth_rate *MGR;

  for (w = 0; w < pd->Num_Species; w++) {
    if (mp->SpeciesSourceModel[w] == FOAM_PMDI_10_H2O) {
      wH2O = w;
    } else if (mp->SpeciesSourceModel[w] == FOAM_PMDI_10_CO2_LIQ) {
      wCO2Liq = w;
    }
  }

  if (!pd->gv[MOMENT1]) {
    GOMA_EH(GOMA_ERROR, "Expected to find moment equations for FOAM_PMDI_10_CO2_LIQ");
  }

  if (wH2O == -1) {
    GOMA_EH(GOMA_ERROR, "Expected to find a speices with source FOAM_PMDI_10_H2O");
    return -1;
  } else if (wCO2Liq == -1) {
    GOMA_EH(GOMA_ERROR, "Expected to find a speices with source FOAM_PMDI_10_CO2_LIQ");
    return -1;
  }

  MGR = calloc(sizeof(struct moment_growth_rate), 1);
  int err = get_moment_growth_rate_term(MGR);
  if (err) {
    free(MGR);
    return -1;
  }

  double Rgas_const = 8.31;
  double ref_press = 1e6;

  if (mp->DensityModel == DENSITY_FOAM_PMDI_10) {
    ref_press = mp->u_density[2];
    Rgas_const = mp->u_density[3];
  } else {
    GOMA_EH(GOMA_ERROR, "Expected DENSITY_FOAM_PMDI_10 density model");
  }

  double source;
  if (T <= 0) {
    source = 0;
    mp->d_species_source[MAX_VARIABLE_TYPES + species_no] = 0;
    mp->d_species_source[MAX_VARIABLE_TYPES + wH2O] = 0;
    mp->d_species_source[TEMPERATURE] = 0;
    free(MGR);
    return 0;
  }

  source = MGR->G[wCO2Liq][1] * ref_press / (Rgas_const * T);

  /**********************************************************/

  /* Species piece */
  eqn = MASS_FRACTION;
  if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
    int j;
    mp->species_source[species_no] = source;
    st->MassSource[species_no] = source;

    /* Jacobian entries for source term */
    var = MASS_FRACTION;
    if (pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        st->d_MassSource_dc[species_no][wCO2Liq][j] =
            MGR->d_G_dC[wCO2Liq][1][j] * ref_press / (Rgas_const * T);
      }
    }

    var = TEMPERATURE;
    if (pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        st->d_MassSource_dT[species_no][j] =
            MGR->d_G_dT[species_no][1][j] * ref_press / (Rgas_const * T);
      }
    }
  }

  free(MGR);
  return 0;
}

/*
 * This is the source term for a suspension model
 * where the fluid and particles have different
 * densities
 */

double epoxy_heat_source(HEAT_SOURCE_DEPENDENCE_STRUCT *d_h,
                         double tt, /* parameter to vary time integration from
                                     * explicit (tt = 1) to implicit (tt = 0) */

                         double dt) /* current time step size */
{
  int eqn, var;
  int j;
  int w;
  int species_no = -1;

  double h = 0.;

  /*  dbl A1, E1, A2, E2;*/
  dbl alpha_dot;
  dbl delta_h;

  /* Begin Execution */

  /* find equation that has extent of reaction */
  for (w = 0; w < pd->Num_Species_Eqn; w++) {
    if (mp->SpeciesSourceModel[w] == EPOXY || mp->SpeciesSourceModel[w] == EPOXY_DEA) {
      /* extent of reaction, alpha */
      species_no = w;
    }
  }
  if (pd->TimeIntegration != STEADY) {
    alpha_dot = fv_dot->c[species_no];
  } else {
    alpha_dot = 0.0;
  }

  delta_h = mp->heat_source;

  /* Species piece */
  eqn = TEMPERATURE;
  if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
    h = delta_h * alpha_dot;

    /* Jacobian entries for source term */
    var = MASS_FRACTION;
    if (d_h != NULL && pd->v[pg->imtrx][var]) {
      if (mp->SpeciesSourceModel[species_no] == EPOXY ||
          (mp->SpeciesSourceModel[species_no] == EPOXY_DEA)) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_h->C[species_no][j] += delta_h * (1 + 2. * tt) / dt * bf[var]->phi[j];
        }
      }
    }
  }
  return (h);
}
/*****************************************************************************/
/* END of routine epoxy_heat_source */
/*****************************************************************************/

/*
 *
 * double butler_volmer_heat_source(dh, a):
 *
 * This routine computes the current source in the solid electrode phase
 * due to an electrochemical reaction that is described by the Butler-Volmer
 * kinetic model - this current source is due to an electrochemical
 * reaction involing species wspec, which either produces electrons
 * (as in hydrogen oxidation reaction in a PEM fuel cell) or consumes
 * electrons (as in oxygen reduction reaction in a PEM fuel cell).
 * Here, the energy transport equation is being used to model the
 * electron transport in the electrode phase such that the electrical
 * potential unknown is actually represented by the temperature unknown.
 *
 *  KSC: 4/25/2006
 *
 *  Modified: 5/9/2006
 *
 * ------------------------------------------------------------------------------
 * This routine is responsible for computing the electrode-phase
 * current source term h, and associated sensitivities at the
 * present gauss point.
 *
 * intput: a[ ]         - parameter array
 * output: h            - current source, A/cm^3
 *         d_h->T[j]    - derivative wrt electrode potential at node j.
 *         d_h->V[j]    - derivative wrt electrolyte potential at node j.
 *         d_h->C[i][j] - derivative wrt concentration of species i at node j
 * ---------------------------------------------------------------------------
 */

double butler_volmer_heat_source(HEAT_SOURCE_DEPENDENCE_STRUCT *d_h, dbl *a) {
  /* Local Variables */
  dbl h;            /* current source or sink, A/cm^3 */
  dbl dh[3];        /* sensitivity vector: dh[0] = dhdT,
                                           dh[1] = dhdV, dh[2] = dhdc;   */
  dbl dhdT;         /* derivative wrt electrode potential                */
  dbl dhdV;         /* derivative wrt electrolyte potential              */
  dbl dhdc;         /* derivative wrt concentration of species wspec     */
  int wspec = a[0]; /* index of the species involved in the reaction   */
  int w, j, var;
  dbl phi_j;

  /* Compute current source and its sensitivities by calling butler_volmer_source */
  h = butler_volmer_source(a, 1, dh);
  dhdT = dh[0];
  dhdV = dh[1];
  dhdc = dh[2];

  /* sensitivies of current source wrt PHI1, PHI2, and c  */
  if (af->Assemble_Jacobian) {
    /* J_s_T --- sensitivity wrt electrode potential */
    var = TEMPERATURE;
    if (d_h != NULL && pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        phi_j = bf[var]->phi[j];
        d_h->T[j] = dhdT * phi_j;
      }
    }

    /* J_s_V --- sensitivity wrt electrolyte potential */
    var = VOLTAGE;
    if (d_h != NULL && pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        phi_j = bf[var]->phi[j];
        d_h->V[j] = dhdV * phi_j;
      }
    }

    /* J_s_c --- sensitivity wrt species concentrations */
    var = MASS_FRACTION;
    if (d_h != NULL && pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        phi_j = bf[var]->phi[j];
        for (w = 0; w < pd->Num_Species_Eqn; w++) {
          d_h->C[w][j] = 0.0; /* no dependency other than wspec */
        }

        d_h->C[wspec][j] = dhdc * phi_j;
      }
    }

  } /* end of if Assemble Jacobian  */

  return (h); /* return the current source */
}
/*****************************************************************************/
/* END of routine butler_volmer_heat_source */
/*****************************************************************************/

/*
 *
 * double butler_volmer_source(a, key, dh):
 *
 * This routine computes the current source/sink or the species source/sink
 * due to an electrochemical reaction that is described by the Butler-Volmer
 * kinetic model - this current or species source/sink is due to an electrochemical
 * reaction involving species wspec, which either consumes species and produces
 * electrons (as in hydrogen oxidation reaction in a PEM fuel cell) or consumes
 * both species and electrons (as in oxygen reduction reaction in a PEM fuel cell).
 * Here, the energy transport equation is being used to model the electron
 * transport in the electrode phase such that the electrical potential
 * unknown is actually represented by the temperature unknown.
 *
 *  KSC: 4/25/2006;
 *
 *  Modified: 5/9/2006
 *
 * ------------------------------------------------------------------------------
 * This routine is responsible for computing the source term h, and
 * associated sensitivities at the current gauss point.
 *
 * intput: a[ ]         - parameter array
 *         key          - 1 denotes current source/sink and 2 indicates species source/sink
 * output: h            - species source (moles/cm^3) or current source (A/cm^3)
 *         dh[0]=dhdT   - derivative of h wrt electrode potential
 *         dh[1]=dhdV   - derivative of h wrt electrolyte potential
 *         dh[2]=dcdV   - derivative of h wrt concentration of species wspec
 *
 * ---------------------------------------------------------------------------
 */

double butler_volmer_source(dbl *a, int key, dbl *dh) {
  /* Local Variables */
  const double R = 8.314;  /* Universal gas constant in units of Joules/mole K */
  const double F = 96487.; /* Faraday's constant in units of Coulomb/mole      */
  dbl FRT;
  dbl s, ai0, beta, cref, alphaa, alphac, T, U0, n = 0.0, PHI1, PHI2, c;
  dbl cratio, conc, conc1, grpa, grpc;
  dbl h, dhdT, dhdV, dhdc;
  int wspec;

  /***************** Source parameters ******************************/
  wspec = a[0];  /* species number */
  s = a[1];      /* stoichiometric coefficient */
  ai0 = a[2];    /* product of interfacial area by
                    exchange current density, A/cm^3 */
  beta = a[3];   /* reaction order */
  cref = a[4];   /* reference species concentration, moles/cm^3 */
  alphaa = a[5]; /* anodic transfer coeficient */
  alphac = a[6]; /* cathodic transfer coefficient */
  T = a[7];      /* temperature, K */
  U0 = a[8];     /* open-circuit potential, V */
  if (key == 2)
    n = a[9];      /* for species source/sink, get the number
                      of electrons involved in the reaction */
  FRT = F / R / T; /* F/R/T */

  /***********Load up convenient local variables*************/
  PHI1 = fv->T;     /* electrode-phase potential   */
  PHI2 = fv->V;     /* electrolyte-phase potential */
  c = fv->c[wspec]; /* load up species concentration */
  if (c == 0)
    c = cref;
  if (c < 0.0)
    c = 1.0e-10;

  /* auxiliary terms */
  cratio = c / cref;
  conc = pow(cratio, beta);
  conc1 = beta * pow(cratio, beta - 1.0);
  grpa = alphaa * FRT * (PHI1 - PHI2 - U0);
  grpc = alphac * FRT * (PHI1 - PHI2 - U0);

  /* compute current produced or consumed by species wspec using Butler-Volmer kinetics */
  h = -s * ai0 * conc * (exp(grpa) - exp(-grpc));

  /* sensitivies of source or sink term h wrt PHI1, PHI2, and c */
  dhdT = -s * ai0 * conc * FRT * (alphaa * exp(grpa) + alphac * exp(-grpc));
  dhdV = s * ai0 * conc * FRT * (alphaa * exp(grpa) + alphac * exp(-grpc));
  dhdc = -s * ai0 * conc1 * (exp(grpa) - exp(-grpc));

  /* when computing species source/sink is desired, convert current
     source/sink to species source/sink using Faraday's law */
  if (key == 2) {
    h /= n * F;
    dhdT /= n * F;
    dhdV /= n * F;
    dhdc /= n * F;
  }

  dh[0] = dhdT;
  dh[1] = dhdV;
  dh[2] = dhdc;

  return (h);
}
/*****************************************************************************/
/* END of routine butler_volmer_source */
/*****************************************************************************/

/*
 *
 * This is the source term for a drying model where
 * the film density is changing as more solvent evaporates
 * Amy C. Sun  6/4/99
 */

double vary_rho_heat_source(HEAT_SOURCE_DEPENDENCE_STRUCT *d_h,
                            double tt, /* parameter to vary time integration from
                                        * explicit (tt = 1) to implicit (tt = 0) */
                            double dt) /* current time step size */
{
  int eqn, var;
  int w, j;

  dbl T, Cp, sv_p;
  dbl sv[MAX_CONC], d_rho_dot_dc[MAX_CONC];
  dbl rho_dot = 0;
  double h = 0.;

  /* Begin Execution */
  /* this is mass_concentration based
      7/99 ACS */
  T = fv->T;
  Cp = mp->heat_capacity;

  if (mp->DensityModel == SOLVENT_POLYMER) {
    /* added ACSun 7/99 */
    rho_dot = 0.;

    for (w = 0; w < pd->Num_Species_Eqn; w++) {
      sv[w] = mp->specific_volume[w];
    }
    sv_p = mp->specific_volume[pd->Num_Species_Eqn];

    /* mass concentration formulation   */
    if (pd->TimeIntegration != STEADY) {
      for (w = 0; w < pd->Num_Species_Eqn; w++) {
        if (fv->c[w] < 0) {
          rho_dot += (fv_dot->c[w]) * (1. - sv[w] / sv_p);
          d_rho_dot_dc[w] = 0.;
        } else {
          rho_dot += (fv_dot->c[w]) * (1. - sv[w] / sv_p);
          d_rho_dot_dc[w] = (1 + 2. * tt) * (1. - sv[w] / sv_p) / dt;
        }
      }
    }
  } else {
    GOMA_EH(GOMA_ERROR, "Must have Density Model set to SOLVENT_POLYMER");
  }

  eqn = TEMPERATURE;

  if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
    h = T * Cp * rho_dot;

    /* Jacobian entries for source term */
    var = MASS_FRACTION;
    if (d_h != NULL && pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        for (w = 0; w < pd->Num_Species_Eqn; w++) {
          d_h->C[w][j] -= T * Cp * d_rho_dot_dc[w] * bf[var]->phi[j];
        }
      }
    }
    var = TEMPERATURE;
    if (pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_h->T[j] -= Cp * rho_dot * bf[var]->phi[j];
      }
    }
  }
  return (h);
}
/******************************************************************************/
/* END of routine vary_rho_heat_source */
/******************************************************************************/

/*
 * For use with fluorinert foam expansion modelling
 */

double foam_heat_source(
    HEAT_SOURCE_DEPENDENCE_STRUCT *d_h,
    double tt, /* parameter to vary time integration from explicit (tt = 1) to implicit (tt = 0) */
    double dt) /* current time step size */
{
  double hT, Tb, a0, phi0;
  double h;

  int j, var;

  hT = mp->u_heat_source[1];
  Tb = mp->u_heat_source[2];
  a0 = mp->u_heat_source[3];
  phi0 = mp->u_heat_source[4];

  h = -hT * (fv->T - Tb) * (3.0 * phi0 / 2.0 / a0);

  var = TEMPERATURE;

  if (pd->v[pg->imtrx][var]) {
    for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
      d_h->T[j] += -bf[var]->phi[j] * (hT * 3.0 * phi0 / 2.0 / a0);
    }
  }
  return (h);
}

double foam_pmdi_10_heat_source(
    HEAT_SOURCE_DEPENDENCE_STRUCT *d_h,
    double time,
    double tt, /* parameter to vary time integration from explicit (tt = 1) to implicit (tt = 0) */
    double dt) /* current time step size */
{
  double h;
  double rho;
  DENSITY_DEPENDENCE_STRUCT d_rho_struct;
  DENSITY_DEPENDENCE_STRUCT *d_rho = &d_rho_struct;

  if (mp->DensityModel != DENSITY_FOAM_PMDI_10) {
    GOMA_EH(GOMA_ERROR, "Expected FOAM_PMDI_10 Density Model for FOAM_PMDI Heat Source");
    mp->heat_capacity = 0.0;
    return 0.0;
  }

  rho = density(d_rho, time);

  int wRXN = -1;
  int w;

  for (w = 0; w < pd->Num_Species_Eqn; w++) {
    if (mp->SpeciesSourceModel[w] == FOAM_PMDI_10_RXN) {
      wRXN = w;
      break;
    }
  }

  if (wRXN == -1) {
    GOMA_EH(GOMA_ERROR, "Expected to find a species with source type FOAM_PMDI_10_RXN");
    return 0;
  }

  int j, var;

  double Delta_H_RXN = mp->u_heat_source[0];
  double M_CO2 = mp->u_density[0];
  double ref_press = mp->u_density[2];

  double Rgas_const = mp->u_density[3];

  double rho_gas = 0;

  rho_gas = (ref_press * M_CO2 / (Rgas_const * fv->T));

  double Y = 1.0 - rho_gas / rho;

  h = Delta_H_RXN * Y * rho * fv_dot->c[wRXN];

  if (d_h != NULL) {
    var = TEMPERATURE;

    if (pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_h->T[j] = 0;
        d_h->T[j] = Delta_H_RXN * (-rho_gas / (rho * fv->T)) * rho * fv_dot->c[wRXN];
        d_h->T[j] += Delta_H_RXN * Y * d_rho->T[j] * fv_dot->c[wRXN] +
                     Delta_H_RXN * (-rho_gas / (rho * rho)) * d_rho->T[j] * fv_dot->c[wRXN];
      }
    }

    var = MASS_FRACTION;
    if (pd->v[pg->imtrx][var]) {
      for (w = 0; w < pd->Num_Species_Eqn; w++) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_h->C[w][j] = Delta_H_RXN * Y * d_rho->C[w][j] * fv_dot->c[wRXN] +
                         Delta_H_RXN * rho_gas / rho * d_rho->C[w][j] * fv_dot->c[wRXN];
          if (w == wRXN) {
            d_h->C[w][j] += Delta_H_RXN * Y * rho * bf[var]->phi[j] * (1 + 2. * tt) / dt;
          }
        }
      }
    }
  }
  return (h);
}

/*
 * Joule Heating Source Model
 */

/*
 * double joule_heat_source (d_h, time)
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
 *
 *   NB: The user need only supply f, dfdT, dfdC, etc....mp struct is loaded up for you
 * --------------------------------------------------------------------------------
 */

double joule_heat_source(HEAT_SOURCE_DEPENDENCE_STRUCT *d_h, dbl time)

{
  int var;
  int dim;
  int a = 0, b;

  int w;

  dbl k = 1e12;       /* electrical conductivity. */
  dbl dkdT[MDE];      /* Temperature derivative of electrical conductivity. */
  dbl dkdV[MDE];      /* Voltage derivative of electrical conductivity. */
  dbl dkdX[DIM][MDE]; /* Spatial derivatives of t.c. */
  int j;
  double h = 0.;

  dim = pd->Num_Dim;

  memset(dkdT, 0, sizeof(double) * MDE);
  memset(dkdV, 0, sizeof(double) * MDE);
  memset(dkdX, 0, sizeof(double) * DIM * MDE);
  if (mp->Elec_ConductivityModel == USER) {
    usr_electrical_conductivity(mp->u_electrical_conductivity, time);

    k = mp->electrical_conductivity;

    var = TEMPERATURE;
    for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
      dkdT[j] = mp->d_electrical_conductivity[var] * bf[var]->phi[j];
    }

    var = VOLTAGE;
    for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
      dkdV[j] = mp->d_electrical_conductivity[var] * bf[var]->phi[j];
    }

    if (pd->v[pg->imtrx][MESH_DISPLACEMENT1]) {
      for (a = 0; a < dim; a++) {
        var = MESH_DISPLACEMENT1 + a;
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          dkdX[a][j] = mp->d_electrical_conductivity[var] * bf[var]->phi[j];
        }
      }
    }

    if (pd->v[pg->imtrx][MASS_FRACTION]) {
      for (w = 0; w < pd->Num_Species_Eqn; w++) {
        var = MASS_FRACTION;
      }
    }
  } else if (mp->Elec_ConductivityModel == CONSTANT) {

    k = mp->electrical_conductivity;

    if (pd->v[pg->imtrx][MASS_FRACTION]) {
      for (w = 0; w < pd->Num_Species_Eqn; w++) {
        var = MASS_FRACTION;
      }
    }
  } else {
    GOMA_EH(GOMA_ERROR, "Unimplemented electrical conductivity model");
  }

  /* Load current density */
  dbl J[DIM];
  if (mp->len_u_heat_source == 0) {
    for (a = 0; a < dim; a++) {
      J[a] = fv->grad_V[a];
    }
  } else if (mp->len_u_heat_source == 1) {
#if MAX_EXTERNAL_FIELD < 3
    GOMA_EH(GOMA_ERROR, "User Joule Heating source expected MAX_EXTERNAL_FIELD >= 3");
#else
    char err_msg[MAX_CHAR_IN_INPUT];
    b = (int)mp->u_heat_source[0];
    sprintf(err_msg, "Joule heating using external fields %s %s %s.\n", efv->name[a + 0],
            efv->name[a + 1], efv->name[a + 2]);
    GOMA_WH(-1, err_msg);
    for (a = 0; a < dim; a++) {
      J[a] = fv->external_field[a + b] / k;
    }
#endif
  } else if (mp->len_u_heat_source == 2) {
#if MAX_EXTERNAL_FIELD < 3
    GOMA_EH(GOMA_ERROR, "User Joule Heating source expected MAX_EXTERNAL_FIELD >= 3");
#else
    char err_msg[MAX_CHAR_IN_INPUT];
    b = (int)mp->u_heat_source[0];
    dbl scale = mp->u_heat_source[1];
    sprintf(err_msg, "Joule heating using external fields %s %s %s, with J scaled by %e.\n",
            efv->name[a + 0], efv->name[a + 1], efv->name[a + 2], scale);
    GOMA_WH(-1, err_msg);
    for (a = 0; a < dim; a++) {
      J[a] = scale * fv->external_field[a + b] / k;
    }
#endif
  } else {
    GOMA_EH(GOMA_ERROR, "Woah, not sure you have the right inputs to the Joule heating source");
  }

  /**********************************************************/
  /* Source sigma*Grad(phi).Grad(phi) */
  h = 0.;
  for (a = 0; a < dim; a++) {
    h += k * J[a] * J[a];
  }

  /* Now do sensitivies */

  var = VOLTAGE;
  if (d_h != NULL && pd->v[pg->imtrx][var]) {
    for (a = 0; a < dim; a++) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_h->V[j] += k * 2. * fv->grad_V[a] * bf[var]->grad_phi[j][a] +
                     dkdV[j] * fv->grad_V[a] * fv->grad_V[a];
      }
    }
  }

  var = TEMPERATURE;
  if (d_h != NULL && pd->v[pg->imtrx][var]) {
    for (a = 0; a < dim; a++) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_h->T[j] += dkdT[j] * fv->grad_V[a] * fv->grad_V[a];
      }
    }
  }

  if (d_h != NULL && pd->v[pg->imtrx][MESH_DISPLACEMENT1]) {
    for (a = 0; a < dim; a++) {
      for (b = 0; b < dim; b++) {
        var = MESH_DISPLACEMENT1 + b;
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_h->X[a][j] += k * 2. * fv->grad_V[a] * fv->d_grad_V_dmesh[a][b][j] +
                          dkdX[a][j] * fv->grad_V[a] * fv->grad_V[a];
        }
      }
    }
  }

  return (h);
}
/*
 * Viscous Dissipation Heating Source Model
 */

/*
 * double visc_diss_heat_source (dh, param)
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

double visc_diss_heat_source(HEAT_SOURCE_DEPENDENCE_STRUCT *d_h,
                             dbl *param) /* General multipliers */
{
  int var;
  int p, q, a, b;

  dbl gamma_dot[DIM][DIM];
  VISCOSITY_DEPENDENCE_STRUCT d_mu_struct; /* viscosity dependence */
  VISCOSITY_DEPENDENCE_STRUCT *d_mu = &d_mu_struct;
  dbl mu;

  dbl gammadot; /* strain rate invariant */

  dbl d_gd_dv[DIM][MDE];    /* derivative of strain rate invariant
                               wrt velocity */
  dbl d_gd_dmesh[DIM][MDE]; /* derivative of strain rate invariant
                               wrt mesh */

  dbl s[DIM][DIM]; /* viscoelastic stress variable */
  int mode;
  int v_s[MAX_MODES][DIM][DIM];

  int j;
  double h = 0.;

  /* Begin Execution */

  /**********************************************************/
  /* Source constant * viscosity* gammadot .. grad_v */

  /* Compute gamma_dot[][] */

  for (a = 0; a < VIM; a++) {
    for (b = 0; b < VIM; b++) {
      gamma_dot[a][b] = fv->grad_v[a][b] + fv->grad_v[b][a];
    }
  }

  /* set stress to zero if it is not defined */
  memset(s, 0, sizeof(dbl) * DIM * DIM);
  if (pd->v[pg->imtrx][POLYMER_STRESS11]) {

    (void)stress_eqn_pointer(v_s);

    for (a = 0; a < VIM; a++) {
      for (b = 0; b < VIM; b++) {
        for (mode = 0; mode < vn->modes; mode++) {
          s[a][b] += fv->S[mode][a][b];
        }
        h += s[a][b] * gamma_dot[b][a];
      }
    }
  }

  h *= 0.5;

  /* get mu and grad_mu */

  mu = viscosity(gn, gamma_dot, d_mu);
  calc_shearrate(&gammadot, gamma_dot, d_gd_dv, d_gd_dmesh);

  h += mu * gammadot * gammadot;

  h *= param[0];

  /* Now do sensitivies */
  if (af->Assemble_Jacobian) {

    var = TEMPERATURE;
    if (d_h != NULL && pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_h->T[j] += param[0] * d_mu->T[j] * gammadot * gammadot;
      }
    }

    if (d_h != NULL && pd->v[pg->imtrx][MESH_DISPLACEMENT1]) {
      for (b = 0; b < VIM; b++) {
        var = MESH_DISPLACEMENT1 + b;
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_h->X[b][j] += param[0] * (2. * mu * gammadot * d_gd_dmesh[b][j] +
                                      gammadot * gammadot * d_mu->X[b][j]);
          if (pd->v[pg->imtrx][POLYMER_STRESS11]) {
            for (p = 0; p < VIM; p++) {
              for (q = 0; q < VIM; q++) {
                d_h->X[b][j] += param[0] * s[p][q] * fv->d_grad_v_dmesh[q][p][b][j];
              }
            }
          }
        }
      }
    }

    if (d_h != NULL && pd->v[pg->imtrx][VELOCITY1]) {
      for (b = 0; b < VIM; b++) {
        var = VELOCITY1 + b;
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_h->v[b][j] +=
              param[0] * (2. * mu * gammadot * d_gd_dv[b][j] + gammadot * gammadot * d_mu->v[b][j]);
          if (pd->v[pg->imtrx][POLYMER_STRESS11]) {
            for (p = 0; p < VIM; p++) {
              for (q = 0; q < VIM; q++) {
                d_h->v[b][j] += param[0] * s[p][q] * bf[VELOCITY1]->grad_phi_e[j][b][q][p];
              }
            }
          }
        }
      }
    }

    if (d_h != NULL && pd->v[pg->imtrx][POLYMER_STRESS11]) {
      for (mode = 0; mode < vn->modes; mode++) {
        for (a = 0; a < VIM; a++) {
          for (b = 0; b < VIM; b++) {
            var = v_s[mode][a][b];
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              d_h->S[mode][a][b][j] = 0.5 * param[0] * gamma_dot[b][a] * bf[var]->phi[j];
            }
          }
        }
      }
    }

  } /* end of if Assemble Jacobian  */

  return (h);
}

double foam_pmdi_10_heat_cap(HEAT_CAPACITY_DEPENDENCE_STRUCT *d_Cp, double time) {
  int var, j;
  int w;
  double Cp, Cp_liq, Cp_gas;
  int wCO2;
  int wH2O;
  DENSITY_DEPENDENCE_STRUCT d_rho_struct;
  DENSITY_DEPENDENCE_STRUCT *d_rho = &d_rho_struct;
  double rho;

  Cp_liq = mp->u_heat_capacity[0];
  Cp_gas = mp->u_heat_capacity[1];

  if (mp->DensityModel != DENSITY_FOAM_PMDI_10) {
    GOMA_EH(GOMA_ERROR, "Expected FOAM_PMDI_10 Density Model for FOAM_PMDI Heat Capacity");
    mp->heat_capacity = 0.0;
    return 0.0;
  }

  rho = density(d_rho, time);

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

  if (wCO2 == -1 && !pd->gv[MOMENT1]) {
    GOMA_EH(GOMA_ERROR, "Expected a Species Source of FOAM_PMDI_10_CO2 or Moment equations");
  } else if (wH2O == -1) {
    GOMA_EH(GOMA_ERROR, "Expected a Species Source of FOAM_PMDI_10_H2O");
  }

  double volF = mp->volumeFractionGas;

  double M_CO2 = mp->u_density[0];
  double rho_liq = mp->u_density[1];
  double ref_press = mp->u_density[2];
  double Rgas_const = mp->u_density[3];

  double rho_gas = 0;

  rho_gas = (ref_press * M_CO2 / (Rgas_const * fv->T));

  Cp = (Cp_liq * rho_liq * (1 - volF) + Cp_gas * rho_gas * volF) / rho;

  /* Now do sensitivies */

  if (d_Cp != NULL && af->Assemble_Jacobian) {
    if (pd->v[pg->imtrx][TEMPERATURE]) {
      var = TEMPERATURE;
      if (pd->v[pg->imtrx][var]) {
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          d_Cp->T[j] = -Cp / rho * d_rho->T[j];
          d_Cp->T[j] += (Cp_liq * rho_liq * mp->d_volumeFractionGas[var]) / rho * bf[var]->phi[j];
          d_Cp->T[j] += ((Cp_gas * rho_gas * mp->d_volumeFractionGas[var]) / rho) * bf[var]->phi[j];
        }
      }
    }

    if (pd->v[pg->imtrx][MASS_FRACTION]) {
      var = MASS_FRACTION;
      if (pd->v[pg->imtrx][var]) {
        for (w = 0; w < pd->Num_Species_Eqn; w++) {
          int wvar = MAX_VARIABLE_TYPES + w;
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            d_Cp->C[w][j] = -Cp / rho * d_rho->C[w][j];
            d_Cp->C[w][j] += ((Cp_liq * rho_liq * (mp->d_volumeFractionGas[wvar]) +
                               Cp_gas * rho_gas * mp->d_volumeFractionGas[wvar]) /
                              rho) *
                             bf[var]->phi[j];
          }
        }
      }
    }
  }

  mp->heat_capacity = Cp;
  return (Cp);
}

/*
 * Enthalpy Heat Capacity Model
 */

/*
 * double enthalpy_heat_capacity_model ( d_Cp)
 *
 * ---------------------------------------------------------------------------
 * This routine is responsible for filling up the following properties and
 * sensitivities at the current gauss point:
 *
 *     input:   <none>
 *
 *     output:  Cp             - heat source
 *              d_Cp->T[j]    - derivative wrt temperature at node j.
 *              d_Cp->V[j]    - derivative wrt voltage at node j.
 *              d_Cp->C[i][j] - derivative wrt mass frac species i at node j
 *              d_Cp->v[i][j] - derivative wrt velocity component i at node j
 *              d_Cp->X[0][j] - derivative wrt mesh displacement components i at node j
 *
 *   NB: The user need only supply f, dfdT, dfdC, etc....mp struct is loaded up
 *       for you
 * ----------------------------------------------------------------------------
 */

double enthalpy_heat_capacity_model(HEAT_CAPACITY_DEPENDENCE_STRUCT *d_Cp) {
  int var;
  int a, b;
  dbl latent_heat, T_liq, T_sol;
  dbl Cp;

  int w;

  dbl grad_H[DIM], d_grad_H_dT[DIM][MDE], norm_grad_H, norm_grad_T, Cp_base, T_ave;
  int i, j;

  /* Begin Execution */

  latent_heat = mp->latent_heat_fusion[0]; /*single component for now*/
                                           /* if this depends on conc'n */
                                           /* then add section here     */
  T_liq = mp->melting_point_liquidus;
  T_sol = mp->melting_point_solidus;
  Cp_base = mp->heat_capacity;

  /************Initialize everything for saftey**************/

  Cp = 0.;

  var = VOLTAGE;
  if (d_Cp != NULL && pd->e[pg->imtrx][var]) {
    for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
      d_Cp->V[j] = 0.;
    }
  }

  if (d_Cp != NULL && pd->v[pg->imtrx][MESH_DISPLACEMENT1]) {
    for (b = 0; b < DIM; b++) {
      var = MESH_DISPLACEMENT1 + b;
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_Cp->X[b][j] = 0.;
      }
    }
  }

  if (d_Cp != NULL && pd->v[pg->imtrx][MASS_FRACTION]) {
    for (w = 0; w < pd->Num_Species_Eqn; w++) {
      var = MASS_FRACTION;
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_Cp->C[w][j] = 0.;
      }
    }
  }

  if (d_Cp != NULL && pd->v[pg->imtrx][VELOCITY1]) {
    for (b = 0; b < DIM; b++) {
      var = VELOCITY1 + b;
      ;
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_Cp->v[b][j] = 0.;
      }
    }
  }

  for (a = 0; a < DIM; a++)
    grad_H[a] = 0.;

  /**********************************************************/
  /* This is basically a hat function using tanh to smooth out
     the corners*/

  /* compute the average temperature of all nodes in this element */
  T_ave = 0.;
  for (i = 0; i < ei[pg->imtrx]->dof[TEMPERATURE]; i++) {
    T_ave += *esp->T[i];
  }
  T_ave = T_ave / ei[pg->imtrx]->dof[TEMPERATURE];

  for (a = 0; a < ei[pg->imtrx]->ielem_dim; a++) {
    for (i = 0; i < ei[pg->imtrx]->dof[TEMPERATURE]; i++) {
      if (*esp->T[i] < (T_sol)) {
        grad_H[a] += bf[TEMPERATURE]->grad_phi[i][a] * (*esp->T[i] * Cp_base);
      } else if (*esp->T[i] > T_sol && *esp->T[i] < T_liq) {
        grad_H[a] += bf[TEMPERATURE]->grad_phi[i][a] *
                     (*esp->T[i] * Cp_base + latent_heat * (*esp->T[i]) / (T_liq - T_sol));
      } else if (*esp->T[i] > T_liq) {
        grad_H[a] += bf[TEMPERATURE]->grad_phi[i][a] * (*esp->T[i] * Cp_base);
      }
    }
  }

  /* calculate norms and final heat capacity to be sent back */

  norm_grad_H = norm_grad_T = 0.;

  for (a = 0; a < ei[pg->imtrx]->ielem_dim; a++)

  {
    norm_grad_H += grad_H[a] * grad_H[a];
    norm_grad_T += fv->grad_T[a] * fv->grad_T[a];
  }

  if (fabs(norm_grad_T) <= 1.e-5) {
    Cp = Cp_base;
    if (T_ave > T_sol && T_ave < T_liq) {
      Cp = Cp_base + latent_heat / (T_liq - T_sol);
    }
  } else {
    Cp = sqrt(norm_grad_H) / sqrt(norm_grad_T);
  }

  /* Now do sensitivies */

  if (af->Assemble_Jacobian) {
    if (fabs(norm_grad_T) > 1.e-5) {
      if (d_Cp != NULL && pd->v[pg->imtrx][TEMPERATURE]) {
        var = TEMPERATURE;
        if (pd->v[pg->imtrx][var]) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            for (a = 0; a < ei[pg->imtrx]->ielem_dim; a++) {

              if (*esp->T[j] <= T_sol) {
                d_grad_H_dT[a][j] = bf[TEMPERATURE]->grad_phi[j][a] * (Cp_base);
              } else if (*esp->T[j] > T_sol && *esp->T[j] < T_liq) {
                d_grad_H_dT[a][j] = bf[TEMPERATURE]->grad_phi[j][a] *
                                    (Cp_base + latent_heat * (1. / (T_liq - T_sol)));
              } else if (*esp->T[j] > T_liq) {
                d_grad_H_dT[a][j] = bf[TEMPERATURE]->grad_phi[j][a] * (Cp_base);
              }
            }
          }
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            for (a = 0; a < ei[pg->imtrx]->ielem_dim; a++) {
              d_Cp->T[j] +=
                  (sqrt(norm_grad_T) * 0.5 * (2. * grad_H[a] * d_grad_H_dT[a][j]) /
                       sqrt(norm_grad_H) -
                   sqrt(norm_grad_H) * 0.5 *
                       (2. * fv->grad_T[a] * bf[TEMPERATURE]->grad_phi[j][a]) / sqrt(norm_grad_T)) /
                  norm_grad_T;
            }
          }
        }
      }
    }
  }

  return (Cp);
}

/*
 * CONVECTIVE LAGRANGIAN VELOCITY ROTATIONAL MODEL
 */

/*
 * int V_mesh_sfs_model(param, v_mesh_sfs, gnn)
 *
 * ----------------------------------------------------------------------------
 * This routine is responsible for filling up the following velocity vector
 *     input:       v_mesh_sfs = param[0-3]
 *                  param[0] = rotation rate
 *                  param[1] = x-coordinate of center of rotation
 *                  param[2] = y-coordinate of center of rotation
 *                  param[3] = z-coordinate of center of rotation
 *                  gnn      = if -1 means routine is being called at a Gauss point,
 *                              otherwise it is the local dof number of an element and
 *                              being called at a nodal point.
 *
 *
 *     output:  v_mesh_sfs[p]
 *
 * ----------------------------------------------------------------------------
 */

int V_mesh_sfs_model(dbl *param,      /* User parameters (u_v_mesh_sfs) */
                     dbl *v_mesh_sfs, /* The three rotational component*/
                     const int model, /* rotational model	*/
                     int gnn)         /* -1 for gauss point, global node number
                                       * for node point*/
{
  const int dim = pd->Num_Dim;
  int a;

  double X[3], X_0[3], R, theta_;
  double origin[3], dir_angle[3];
  double axis_pt[3], rad_dir[3], v_dir[3], t;

  /* Begin Execution */

  /* First figure out radius of curvature
     at current location */

  for (a = 0; a < dim; a++) {
    if (gnn == -1) {
      X[a] = fv->x0[a]; /*we only care about
                          the stress-free st.*/
    } else {
      X[a] = Coor[a][gnn];
    }
    X_0[a] = param[1 + a]; /*the first param entry
                             is the rotation rate */
  }

  if (model == ROTATIONAL) {
    /* Compute current radial position about
       axis of rotation   */
    R = 0.;
    for (a = 0; a < dim; a++) {
      R += SQUARE(X[a] - X_0[a]);
    }
    R = sqrt(R);

    /* compute the angle                  */
    /* NB for here on this only works     */
    /* for rotation about the Z axis      */
    /*  Need to furbish for general case  */
    theta_ = acos((X[0] - X_0[0]) / R);

    /* Quadrant dependence */
    if ((X[1] - X_0[1]) > 0) {
      v_mesh_sfs[0] = param[0] * R * sin(theta_);
      v_mesh_sfs[1] = -param[0] * R * cos(theta_);
      v_mesh_sfs[2] = 0;
    } else {
      v_mesh_sfs[0] = -param[0] * R * sin(theta_);
      v_mesh_sfs[1] = -param[0] * R * cos(theta_);
      v_mesh_sfs[2] = 0;
    }
  } else if (model == ROTATIONAL_3D) {
    /*  origin and direction of rotation axis	*/
    origin[0] = X_0[0];
    origin[1] = X_0[1];
    origin[2] = X_0[2];
    dir_angle[0] = param[4];
    dir_angle[1] = param[5];
    dir_angle[2] = param[6];

    /*  find intersection of axis with normal plane - i.e., locate point on
            axis that intersects plane normal to axis that contains local point. */

    t = (dir_angle[0] * (X[0] - origin[0]) + dir_angle[1] * (X[1] - origin[1]) +
         dir_angle[2] * (X[2] - origin[2])) /
        (SQUARE(dir_angle[0]) + SQUARE(dir_angle[1]) + SQUARE(dir_angle[2]));
    axis_pt[0] = origin[0] + dir_angle[0] * t;
    axis_pt[1] = origin[1] + dir_angle[1] * t;
    axis_pt[2] = origin[2] + dir_angle[2] * t;

    /*  compute radius and radial direction	*/

    R = sqrt(SQUARE(X[0] - axis_pt[0]) + SQUARE(X[1] - axis_pt[1]) + SQUARE(X[2] - axis_pt[2]));
    rad_dir[0] = (X[0] - axis_pt[0]) / R;
    rad_dir[1] = (X[1] - axis_pt[1]) / R;
    rad_dir[2] = (X[2] - axis_pt[2]) / R;

    /* compute velocity direction as perpendicular to both axis and radial
            direction.  Positive direction is determined by right hand rule */

    v_dir[0] = dir_angle[1] * rad_dir[2] - dir_angle[2] * rad_dir[1];
    v_dir[1] = dir_angle[2] * rad_dir[0] - dir_angle[0] * rad_dir[2];
    v_dir[2] = dir_angle[0] * rad_dir[1] - dir_angle[1] * rad_dir[0];

    v_mesh_sfs[0] = param[0] * R * v_dir[0];
    v_mesh_sfs[1] = param[0] * R * v_dir[1];
    v_mesh_sfs[2] = param[0] * R * v_dir[2];
  } else if (model == OSC_LINEAR) {
    double v_magn = 0, stroke = elc_rs->u_v_mesh_sfs[3];
    double end_stroke = elc_rs->u_v_mesh_sfs[4], end_dt, time_pd, time_red;
    double interval;
    for (a = 0; a < DIM; a++) {
      v_magn += SQUARE(elc_rs->u_v_mesh_sfs[a]);
    }
    v_magn = sqrt(v_magn);
    time_pd = stroke / v_magn;
    end_dt = end_stroke / v_magn;
    time_red = modf(tran->time_value / (4 * time_pd), &interval);
    if (time_red > (time_pd - end_dt) && time_red <= (time_pd + end_dt)) {
      for (a = 0; a < DIM; a++) {
        elc_rs->v_mesh_sfs[a] = -sin((time_red - time_pd) * 2 / M_PIE) * elc_rs->u_v_mesh_sfs[a];
      }
    } else if (time_red > (time_pd + end_dt) && time_red <= (3 * time_pd - end_dt)) {
      for (a = 0; a < DIM; a++) {
        elc_rs->v_mesh_sfs[a] = -elc_rs->u_v_mesh_sfs[a];
      }
    } else if (time_red > (3 * time_pd - end_dt) && time_red <= (3 * time_pd + end_dt)) {
      for (a = 0; a < DIM; a++) {
        elc_rs->v_mesh_sfs[a] = sin((time_red - time_pd) * 2 / M_PIE) * elc_rs->u_v_mesh_sfs[a];
      }
    } else {
      for (a = 0; a < DIM; a++) {
        elc_rs->v_mesh_sfs[a] = elc_rs->u_v_mesh_sfs[a];
      }
    }

  } else {
    GOMA_EH(GOMA_ERROR, "Invalid Rotational V_mesh_sfs_model.");
  }

  return (1);
}
/**********************************************************************************/
/**********************************************************************************/
/**********************************************************************************/

int Diffusivity(void)

/********************************************************************
 * Diffusivity() sets the FICKIAN diffusivity and its derivatives
 * for all species and loads them into the mp structure:
 *
 *  This routine sets the following values:
 *
 *   mp->diffusivity[w]
 *   mp->d_diffusivity[w][var+type]
 *
 * HKM -> Note, I changed the loop count from num_species_eqn to
 *        num_species, since every nondilute diffusion models needs the
 *        diffusivity of the last species even if there isn't an equation
 *        for the last species.
 *        The units for the diffusivity are cm^2 sec-1 in almost all
 *        cases where Species_Var_Type is not SPECIES_UNDEFINED_FORM.
 *        If Species_Var_Type is SPECIES_UNDEFINED_FORM, the usual
 *        free-for-all exists.
 *        Another question is whether DiffusivityModel needs to be
 *        a vector. I would think that the diffusion model would
 *        be the same for all species. If it is, then this routine
 *        could be speeded up.
 *
 * Author: Tom Baer
 * Date  : 11/6/96
 ********************************************************************/
{
  int i, j, w, var, err = 0;
  double T = 298.0, P = 101.325;
  ;
  for (w = 0; w < pd->Num_Species; w++) {
    switch (mp->DiffusivityModel[w]) {

    case CONSTANT:
      /* do nothing...mp->diffusivity set in rd_mp_specs
         ...mp->d_diffusivity zeroed in set_mp_to_unity
      */
      break;

    case USER:
      err = usr_diffusivity(w, mp->u_diffusivity[w]);
      break;

    case POROUS:
      /* do nothing... Rich loads this diffusivity model in
         mm_fill_porous.c.  Maybe this should be
         moved over here sometime ? tab
      */
      break;

    case FREE_VOL:
      err = Free_Vol_Theory_Diffusivity(w, mp->u_diffusivity[w]);
      break;

    case GENERALIZED_FREE_VOL:
      err = Generalized_FV_Diffusivity(w);
      break;

    case TABLE: {
      struct Data_Table *table_local;
      table_local = MP_Tables[mp->diffusivity_tableid[w]];
      apply_table_mp(&mp->diffusivity[w], table_local);
      for (i = 0; i < table_local->columns - 1; i++) {
        var = table_local->t_index[i];
        switch (var) {
        case TEMPERATURE:
          mp->d_diffusivity[w][TEMPERATURE] = table_local->slope[i];
          break;
        case MASS_FRACTION:
          if (pd->v[pg->imtrx][MASS_FRACTION]) {
            for (j = 0; j < pd->Num_Species; j++) {
              mp->d_diffusivity[w][MAX_VARIABLE_TYPES + j] = table_local->slope[i];
            }
          }
          break;
        default:
          GOMA_EH(GOMA_ERROR, "Variable function not yet implemented in material property table");
          break;
        }
      }
    } break;

    case LEVEL_SET:
      err = ls_transport_property(mp->u_diffusivity[w][0], mp->u_diffusivity[w][1],
                                  mp->u_diffusivity[w][2], &mp->diffusivity[w],
                                  &mp->d_diffusivity[w][FILL]);
      break;
    case CHAPMAN_GAS:
      if (pd->gv[TEMPERATURE]) {
        T = fv->T;
      } else {
        T = upd->Process_Temperature;
      }

      if (pd->gv[PRESSURE] && fv->P > 0) {
        P = fv->P;
      } else {
        P = upd->Pressure_Datum / 10000.0;
      }
      /* First coefficient is A/sigma12/omega (0.001859 atm-cm^2-(g/mol)^1/2/(K^1.5-s))  */
      mp->diffusivity[w] = mp->u_diffusivity[w][0] *
                           sqrt(1. / mp->molecular_weight[w] + 1. / mp->u_diffusivity[w][1]) *
                           pow(T, 1.5) / P;
      if (pd->e[pg->imtrx][TEMPERATURE]) {
        mp->d_diffusivity[w][TEMPERATURE] = 1.5 * mp->diffusivity[w] / T;
      }
      if (pd->e[pg->imtrx][PRESSURE]) {
        mp->d_diffusivity[w][PRESSURE] = mp->diffusivity[w] * P * log(P);
      }

      mp->d_diffusivity[w][MAX_VARIABLE_TYPES + w] = 0.0;
      break;
    default:
      break;
    }
  }
  return err;
} /* end of Diffusivity() */

/*
 * int Free_Vol_Theory_Diffusivity(param)
 *
 *      Free Volume Theory Diffusivity Model
 * ------------------------------------------------------------------------------
 * This routine is responsible computing the diffusivity according to
 * free volume theory. See Paper by Duda et al., AIChE J., Vol 28, Pg 279.
 *
 *     input:  mp->u_diffusivity = param[0-16]
 *                  param[0] = Specific critical hole free volume of component 0
 *                             required for a jump
 *                  param[1] = Specific critical hole free volume of component 1
 *                             required for a jump
 *                  param[2] = K11/gamma = free volume parameter 11 divided by overlap
 *                              factor for free volume
 *                  param[3] =  K12/gamma = free volume parameter 12 divided by overlap
 *                              factor for free volume
 *                  param[4] =  K21 - T_gl = free volume parameter 21 minus glass transition
 *                                           temperature for component 1
 *                  param[5] =  K22 - T_2l = free volume parameter 22 minus glass transition
 *                                           temperature for component 2
 *                  param[6] =  chi = interaction parameter for flory-huggins theory
 *                  param[7] = x_si= ratio of critical molar volume of solvent jumping unit
 *                              to critical molar volume of jumping unit of polymer.
 *                  param[8] = pre-exponential factor for diffusivity (D0)
 *                  param[9]= critical energy per mole to overcome attractive forces
 *                            divided by the Natural Gas law constant (K).
 *                  param[10]= specific volume of pure component 1
 *                  param[11]= specific volume of pure component 2
 *
 *                  param[12] = free volume model number
 *                  param[12] = 0 = FV 0  D_11 = (Vrentas and Duda, 1977)
 *                  param[12] = 1 = FV 2  (Zelinsky and Hanley , 1999)
 *                  param[12] = 2 = FV 3   self
 *                  param[12] = 3 = FV 4 (Alsoy and Duda, 1999)
 *                  param[12] = 4 = FVTO (friction-based)
 *
 *                  param[13] = solvent molecular weight
 *                  param[14] = polymer molecular weight
 *                  param[15] = D_o_polymer - polymer D_o for polymer
 *						 self-diffusion coefficient
 *                  param[16] = E_div_R_polymer - E_activation/R - polymer
 *						 self-diffusion model
 *                  Note: D_o_polymer and E_div_R_polymer are used SOLELY
 *		  		in the friction-based model !!
 *
 *     output:  mp->diffusivity[w]
 *              mp->d_diffusivity[w][w1]
 *
 * NOTE: The diffusivity according to this theory does not approach to
 * self-diffusivity (of solvent component) as the mixture approaches pure solvent.
 * One can make a correction by making sure that pow((1 - vol_frac_1),2)*(1.-2.*chi*vol_frac_1)
 * does not go to zero as vol_frac_1 -> 1.  For example, use a cutoff
 *  if (vol_frac_1 > 0.7)  vol_frac_1 = 0.7;
 *
 * NOTE: This subroutine has been modified to based fv->c on species_density,
 *       and NOT MASS_FRACTION. ACSun 5/00
 *
 * Modified to include 4 additional free volume models and is based on
 * Price et al., AIChE J., Vol 49, Pg 309.  APErvin 10/03
 * ------------------------------------------------------------------------------
 */
int Free_Vol_Theory_Diffusivity(int species_no, /* current species number*/
                                double *param)  /* free volume params from mat file */
{
  int w, fv_model_number;

  double T, C[MAX_CONC], dDdT;
  double V_1s, V_2s, K11_gamma, K12_gamma, K21mTg1, K22mTg2, chi, x_si, D_o, E_div_R, V_10, V_20;
  double vol_frac_1, vol_frac_2;
  double d_vol_frac_1_dc[MAX_CONC], dC0dc[MAX_CONC];
  double V_fh_gamma;
  double exponen;
  double d_V_fh_gamma_dc[MAX_CONC], dexponen_dc[MAX_CONC];
  double c0, density_tot; /* solvent dens, total density*/

  double Q_thermo = 0.0; /* the thermodynamic factor */
  double beta;           /* a fudge factor in Q_thermo */
  double MW_1 = 0.0;     /* solvent molecular weight */
  double MW_2 = 0.0;     /* polymer molecular weight */
  double D1 = 0.0;       /* the solvent self diffusion coefficient */
  double D_o_polymer = 0.0, E_div_R_polymer = 0.0, A_friction = 0.0, ratio_D1_D2_term = 0.0;
  /* used solely for friction model */
  double exponen_D2 = 0.0, D2 = 0.0; /* used solely for friction model */
  double d_D1_dc[MAX_CONC];          /* derivative of D1, the solvent
                                        self diffusion coefficient wrt the solvent mf */
  double dQ_dc[MAX_CONC];            /* derivative of Q_thermo wrt the solvent mf */
  double dA_dT;                      /* friction factor model derivative of A wrt T */
  double d_D1_dT;                    /* dD1/dT */
  double d_D2_D1_dT = 0;             /*  d(D2/D1)/dT */
  double d_D2_dc[MAX_CONC];          /*  d(D2/D1)/dw1 */
  double dA_dc[MAX_CONC];            /* friction factor model derivative of A wrt w1 */

  double term1, term2, term3, term4; /* intermediate parameters */
  PROPERTYJAC_STRUCT *densityJac = NULL;
  propertyJac_realloc(&densityJac, mp->Num_Species + 1);

  /* load up convenient quantities */
  V_1s = param[0];
  V_2s = param[1];
  K11_gamma = param[2];
  K12_gamma = param[3];
  K21mTg1 = param[4];
  K22mTg2 = param[5];
  chi = param[6];
  x_si = param[7];
  D_o = param[8];
  E_div_R = param[9];
  V_10 = param[10];
  V_20 = param[11];
  fv_model_number = (int)param[12];

  T = fv->T;

  if (fv_model_number != 0) {
    /* beta = solvent molecular weight/ polymer molecular weight */
    MW_1 = param[13];
    MW_2 = param[14];
    beta = 1.0;
  }
  if (fv_model_number == 4) /* friction model requires these two parameters */
  {
    D_o_polymer = param[15];
    E_div_R_polymer = param[16];
  }

  /*
   * this minor kludge is added for backward compatibility with the testsuite problems
   */

  if (pd->Num_Species_Eqn == 1) {
    if (mp->specific_volume[0] == -1.0)
      mp->specific_volume[0] = V_10;
    if (mp->specific_volume[1] == -1.0)
      mp->specific_volume[1] = V_20;
  }

  /*  initialization  */
  C[0] = fv->c[0];
  C[1] = 1.0 - C[0];
  vol_frac_2 = 1.0;
  /*  multi-component kludge	*/
  c0 = 0.0;
  vol_frac_1 = 0.0;
  memset(d_vol_frac_1_dc, 0, sizeof(dbl) * MAX_CONC);
  memset(dC0dc, 0, sizeof(dbl) * MAX_CONC);
  density_tot = calc_density(mp, TRUE, densityJac, 0.0);
  /*  density_tot = 1./mp->specific_volume[pd->Num_Species_Eqn]; */
  switch (mp->Species_Var_Type) {
  case SPECIES_DENSITY:
    for (w = 0; w < pd->Num_Species_Eqn; w++) {
      if (mp->FreeVolSolvent[w]) {
        c0 += fv->c[w];
        vol_frac_1 += fv->c[w] * mp->specific_volume[w];
      }
    }
    C[0] = c0 / density_tot; /* w1 - the solvent mass fraction */
    C[1] = 1. - C[0];        /* w2 - the polymer mass fraction */
    vol_frac_2 = 1.0 - vol_frac_1;
    if (af->Assemble_Jacobian) {
      for (w = 0; w < pd->Num_Species_Eqn; w++) {
        if (mp->FreeVolSolvent[w]) {
          d_vol_frac_1_dc[w] = mp->specific_volume[w];
          dC0dc[w] =
              (density_tot - c0 * mp->d_density[MAX_VARIABLE_TYPES + w]) / SQUARE(density_tot);
        } else {
          d_vol_frac_1_dc[w] = 0.0;
          dC0dc[w] = -c0 * mp->d_density[MAX_VARIABLE_TYPES + w] / SQUARE(density_tot);
        }
      }
    }
    break;

  case SPECIES_UNDEFINED_FORM:
  case SPECIES_MASS_FRACTION:
    for (w = 0; w < pd->Num_Species_Eqn; w++) {
      if (mp->FreeVolSolvent[w]) {
        c0 += fv->c[w];
        vol_frac_1 += fv->c[w] * mp->specific_volume[w];
      }
    }
    C[0] = c0 / density_tot; /* w1 - the solvent mass fraction */
    C[1] = 1. - C[0];        /* w2 - the polymer mass fraction */
    vol_frac_2 = 1.0 - vol_frac_1;
    if (af->Assemble_Jacobian) {
      for (w = 0; w < pd->Num_Species_Eqn; w++) {
        if (mp->FreeVolSolvent[w]) {
          d_vol_frac_1_dc[w] = mp->specific_volume[w] *
                               (density_tot + fv->c[w] * mp->d_density[MAX_VARIABLE_TYPES + w]);
          dC0dc[w] = 1.0;
        } else {
          d_vol_frac_1_dc[w] =
              mp->specific_volume[w] * fv->c[w] * mp->d_density[MAX_VARIABLE_TYPES + w];
          dC0dc[w] = 0.0;
        }
      }
    }
    break;

  case SPECIES_CONCENTRATION:
    for (w = 0; w < pd->Num_Species_Eqn; w++) {
      if (mp->FreeVolSolvent[w]) {
        c0 += fv->c[w] * mp->molecular_weight[w];
        vol_frac_1 += fv->c[w] * mp->specific_volume[w] * mp->molecular_weight[w];
      }
    }
    C[0] = c0 / density_tot; /* w1 - the solvent mass fraction */
    C[1] = 1. - C[0];        /* w2 - the polymer mass fraction */
    vol_frac_2 = 1.0 - vol_frac_1;
    if (af->Assemble_Jacobian) {
      for (w = 0; w < pd->Num_Species_Eqn; w++) {
        if (mp->FreeVolSolvent[w]) {
          d_vol_frac_1_dc[w] = mp->specific_volume[w] * mp->molecular_weight[w];
          dC0dc[w] =
              (density_tot * mp->molecular_weight[w] - c0 * mp->d_density[MAX_VARIABLE_TYPES + w]) /
              SQUARE(density_tot);
        } else {
          d_vol_frac_1_dc[w] = 0.0;
          dC0dc[w] = -c0 * mp->d_density[MAX_VARIABLE_TYPES + w] / SQUARE(density_tot);
        }
      }
    }
    break;

  case SPECIES_MOLE_FRACTION:
  default:
    GOMA_EH(GOMA_ERROR, "Undefined Species formulation in FREE_VOL_DIFFUSIVITY\n");
    break;
  }

  V_fh_gamma = K11_gamma * C[0] * (K21mTg1 + T) + K12_gamma * C[1] * (K22mTg2 + T);

  D_o *= exp(-E_div_R / T);

  if (fv_model_number != 0) {

    /* Calculate D1, the solvent self-diffusion coefficient */

    exponen = exp(-(C[0] * V_1s + C[1] * x_si * V_2s) / V_fh_gamma);
    D1 = D_o * exponen;

    /* Calculate Q_thermo, the thermodynamic factor as defined
       in the Price 2003 paper*/

    Q_thermo = ((1.0 - vol_frac_1) * (1.0 - (2.0 * vol_frac_1 * chi))) +
               ((vol_frac_1 * V_10 * MW_1) / (V_20 * MW_2 * beta));
  }

  /* For the friction-based model  param[12] = 4 */
  /* we calculate D2 - the polymer self-diffusion coefficient */

  if (fv_model_number == 4) {
    D_o_polymer *= exp(-E_div_R_polymer / T);
    exponen_D2 = exp(-((C[0] * V_1s / x_si) + (C[1] * V_2s)) / V_fh_gamma);
    D2 = D_o_polymer * exponen_D2;
    ratio_D1_D2_term = 1.0 - ((D2 * V_20 * MW_2) / (D1 * V_10 * MW_1));
    ratio_D1_D2_term *= vol_frac_1;
    A_friction = 1.0 - ratio_D1_D2_term;
  }

  switch (fv_model_number) {

  case 0:
    mp->diffusivity[species_no] = D_o * pow((1 - vol_frac_1), 2.0) * (1. - 2. * chi * vol_frac_1) *
                                  exp(-(C[0] * V_1s + C[1] * x_si * V_2s) / V_fh_gamma);
    break;
  case 1:
    mp->diffusivity[species_no] = (vol_frac_2 / C[1]) * Q_thermo * D1;
    break;

  case 2:
    mp->diffusivity[species_no] = D1;
    break;

  case 3:
    mp->diffusivity[species_no] = Q_thermo * D1;
    break;

  case 4:
    mp->diffusivity[species_no] = A_friction * D1 * Q_thermo;
    break;

  default:
    GOMA_EH(GOMA_ERROR, "Invalid Free Volume model number.");
    break;

  } /* end of switch(fv_model_number) */

  if (af->Assemble_Jacobian) {
    if (pd->v[pg->imtrx][TEMPERATURE]) {

      if (fv_model_number != 4) {
        dDdT = mp->diffusivity[0] * ((C[0] * V_1s + C[1] * x_si * V_2s) / pow(V_fh_gamma, 2.0)) *
                   (K11_gamma * C[0] + K12_gamma * C[1]) +
               mp->diffusivity[0] * (E_div_R / T / T);

        mp->d_diffusivity[species_no][TEMPERATURE] = dDdT;
      }
      if (fv_model_number == 4) {
        term1 = (E_div_R_polymer - E_div_R) / (pow(T, 2));
        term2 = V_1s * C[0] * ((1.0 / x_si) - 1.0);
        term3 = V_2s * C[1] * (1.0 - x_si);
        term4 = (K11_gamma * C[0] + K12_gamma * C[1]) * (term2 + term3);

        d_D2_D1_dT = (term4 / (pow(V_fh_gamma, 2))) + term1;
        dA_dT = ((vol_frac_1 * MW_2 * V_20) / (MW_1 * V_10)) * d_D2_D1_dT;
        d_D1_dT = (D1 * (E_div_R / T / T)) +
                  (D1 * ((C[0] * V_1s + C[1] * x_si * V_2s) / pow(V_fh_gamma, 2.0)) *
                   (K11_gamma * C[0] + K12_gamma * C[1]));
        dDdT = (dA_dT * Q_thermo * D1) + (A_friction * Q_thermo * d_D1_dT);

        mp->d_diffusivity[species_no][TEMPERATURE] = dDdT;
      }
    }

    if (pd->v[pg->imtrx][MASS_FRACTION]) {
      /* FOR NOW ONLY ALLOW 2 COMPONENT SYSTEMS, I.E. 1 SPECIES EQUATION) */
#if 0
	  d_vol_frac_1_dC[0] = ((C[0]*V_10 + C[1]*V_20)*V_10 - C[0]*V_10*(V_10-V_20))
	    /pow((C[0]*V_10 + C[1]*V_20),2.0);

#endif
      exponen = exp(-(C[0] * V_1s + C[1] * x_si * V_2s) / V_fh_gamma);

      for (w = 0; w < pd->Num_Species_Eqn; w++) {

        d_V_fh_gamma_dc[w] = dC0dc[w] * (K11_gamma * (K21mTg1 + T) - K12_gamma * (K22mTg2 + T));
        dexponen_dc[w] = (-(V_fh_gamma * (V_1s - x_si * V_2s) * dC0dc[w] -
                            (C[0] * V_1s + C[1] * V_2s * x_si) * d_V_fh_gamma_dc[w]) /
                          pow(V_fh_gamma, 2.0));
        if (fv_model_number != 0) {
          d_D1_dc[w] = D_o * exponen * dexponen_dc[w];
          dQ_dc[w] = d_vol_frac_1_dc[w] * (4.0 * chi * vol_frac_1 - 2.0 * chi +
                                           (V_10 * MW_1) / (V_20 * MW_2 * beta) - 1.0);
        }

        if (fv_model_number == 4) {
          term1 = C[0] * V_1s / x_si + C[1] * V_2s;
          term2 = V_1s / x_si - V_2s;
          d_D2_dc[w] =
              D_o_polymer * exponen_D2 *
              (-(V_fh_gamma * dC0dc[w] * term2 - term1 * d_V_fh_gamma_dc[w]) / SQUARE(V_fh_gamma));
          term3 = V_20 * MW_2 / (V_10 * MW_1);
          dA_dc[w] = vol_frac_1 * term3 * (D1 * d_D2_dc[w] - D2 * d_D1_dc[w]) / SQUARE(D1) -
                     d_vol_frac_1_dc[w] * (1. - term3 * D2 / D1);
        }
      }

      switch (fv_model_number) {
      case 0:
        for (w = 0; w < pd->Num_Species_Eqn; w++) {
          mp->d_diffusivity[species_no][MAX_VARIABLE_TYPES + w] =
              2. * (1. - vol_frac_1) * d_vol_frac_1_dc[w] * (1. - 2. * chi * vol_frac_1) * exponen +
              pow(1. - vol_frac_1, 2.0) * (-2. * chi) * d_vol_frac_1_dc[w] * exponen +
              pow(1. - vol_frac_1, 2.0) * (1. - 2. * chi * vol_frac_1) * exponen * dexponen_dc[w];
          mp->d_diffusivity[species_no][MAX_VARIABLE_TYPES + w] *= D_o;
        }
        break;

      case 1:

        for (w = 0; w < pd->Num_Species_Eqn; w++) {
          mp->d_diffusivity[species_no][MAX_VARIABLE_TYPES + w] =
              (-d_vol_frac_1_dc[w] * Q_thermo * D1 + vol_frac_2 * dQ_dc[w] * D1 +
               vol_frac_2 * Q_thermo * d_D1_dc[w]) /
                  C[1] -
              (vol_frac_2 * Q_thermo * D1) * (-dC0dc[w]) / SQUARE(C[1]);
        }
        break;

      case 2:
        for (w = 0; w < pd->Num_Species_Eqn; w++) {
          mp->d_diffusivity[species_no][MAX_VARIABLE_TYPES + w] = d_D1_dc[w];
        }
        break;

      case 3:
        for (w = 0; w < pd->Num_Species_Eqn; w++) {
          mp->d_diffusivity[species_no][MAX_VARIABLE_TYPES + w] =
              dQ_dc[w] * D1 + Q_thermo * d_D1_dc[w];
        }
        break;

      case 4:

        for (w = 0; w < pd->Num_Species_Eqn; w++) {
          mp->d_diffusivity[species_no][MAX_VARIABLE_TYPES + w] =
              dA_dc[w] * Q_thermo * D1 + dQ_dc[w] * A_friction * D1 +
              A_friction * Q_thermo * d_D1_dc[w];
        }
        break;
      }
    }
  }
  propertyJac_destroy(&densityJac);
  return (0);
}
/*
 * Generalized_Diffusivity() sets the GENERALIZED FICKIAN
 * diffusivity and its derivatives for all species
 * and loads them into the mp structure.
 *
 * Author: Amy Sun
 * Date  : 4/13/00
 *
 */
int Generalized_Diffusivity(void) {
  int w, w1, err = 0; /* return err != 0 to signal error */
  for (w = 0; w < pd->Num_Species_Eqn; w++) {
    switch (mp->DiffusivityModel[w]) {
    case GENERALIZED:
      /* constant mutual diffusivity */
      for (w1 = 0; w1 < pd->Num_Species_Eqn; w1++) {
        mp->diffusivity_gen_fick[w][w1] = mp->u_diffusivity[w][w1];
      }
      break;

    case GENERALIZED_FREE_VOL:
      err = Generalized_FV_Diffusivity(w);
      break;

    default:
      break;
    } /* end of switch(DiffusivityModel) */

  } /* end of for( w ) */

  return err;

} /* end of Generalized_Diffusivity() */

/*
 * Generalized Free Volume Theory Diffusivity Model
 */
/*
 * int Generalized_FV_Diffusivity()
 *
 * ------------------------------------------------------------------------------
 * This routine is responsible computing the multicomponent diffusivity according to
 * the free volume theory.  This is generalization of ternary systems.
 * See Paper by Vrentas, Duda, and Ling. J. of App.Poly.Sci. vol30, 4499 (1985).
 *
 *     input:   mp->u_diffusivity[w][] = Free volume parameters for component w
 *              Note that this is the same order as the last routine
 *
 * The way the input is set up it is identical to the format of binary, except that
 * the original "2" designation always refers to the polymer.  Hence, constants
 * associated with pure polymer alone is repeated several times. i.e.
 *              param[w][1] = Specific critical hole free volume of polymer
 *                             required for a jump
 *              param[w][3] =  K1p/gamma for polymer
 *              param[w][5] =  K2p - T_2l for polymer
 *              param[w][11]= specific volume of polymer
 *
 * for each solvent w. Only one set needs to be assigned.
 * For solvent specific free volume parameters. The following is some guidance.
 *
 *              param[w][0] = Specific critical hole free volume of solvent w
 *                            required for a jump
 *              param[w][2] = K1w/gamma for solvent w
 *              param[w][4] = K21 - T_gl for solvent w
 *              param[w][6] =  chi = interaction parameter for flory-huggins theory
 *              param[w][7] =  ratio of critical molar volume of solvent jumping unit
 *                              to critical molar volume of jumping unit of polymer.
 *              param[w][8] = pre-exponential factor for diffusivity (D0)
 *              param[w][9]= critical energy per mole to overcome attractive forces
 *                            divided by the Natural Gas law constant (K).
 *              param[w][10]= specific volume of pure component w
 *
 *     output:  mp->diffusivity[w] <- standard placeholder for pure diffusivity
 *              mp->d_diffusivity[w][var]
 *              mp->diffusivity_gen_fick[w][w1]
 *              mp->d_diffusivity_gf[w][w1][var]
 *
 *     Author:  A.C. Sun 5/00
 *
 * --------------------------------------------------------------------------------
 */
int Generalized_FV_Diffusivity(int species_no) /* current species number*/
{
  int w, w1, mode;
  int w2 = pd->Num_Species_Eqn;

  double activity[MAX_CONC], dcoeff[MAX_CONC][MAX_CONC], d2coeff[MAX_CONC][MAX_CONC][MAX_CONC];

  double T;

  double d_Vfh_domega[MAX_CONC];
  double d_num_domega[MAX_CONC][MAX_CONC];
  double d_omega_dc[MAX_CONC][MAX_CONC];
  double numerator[MAX_CONC], omega[MAX_CONC];
  double xsi[MAX_CONC][MAX_CONC];
  double vs[MAX_CONC], K1_gamma[MAX_CONC], K2mTg[MAX_CONC], Do[MAX_CONC], E_div_R[MAX_CONC];
  double vsp, K1p_gamma, K2pmTg;
  double allomega = 0.;
  double V_fh_gamma = 0.;
  double d_Vfh_dT = 0.;
  double density_tot, rho[MAX_CONC], drho_dc[MAX_CONC];
  PROPERTYJAC_STRUCT *densityJac = NULL;
  propertyJac_realloc(&densityJac, mp->Num_Species + 1);

  /* Begin Execution */

  /* load up activity coefficient information, this routine
   * there is no need to calculate mass transfer across the
   * interface  */
  /* the mode is hardcoded to FLORY for now */

  memset(dcoeff, 0, sizeof(dbl) * MAX_CONC * MAX_CONC);
  memset(d2coeff, 0, sizeof(dbl) * MAX_CONC * MAX_CONC * MAX_CONC);

  /* ln(activity_coeff) is coded based on flory-huggins for now
     there may be interest to expand this to other models,
     such as UNIQUAC, NRTL, etc. */

  mode = FLORY;
  act_coeff(activity, dcoeff, d2coeff, fv->c, mode, species_no);

  /* load up polymer FV parameters */

  vsp = mp->u_diffusivity[0][1];
  K1p_gamma = mp->u_diffusivity[0][3];
  K2pmTg = mp->u_diffusivity[0][5];

  /* load up solvent FV parameters */
  for (w = 0; w < pd->Num_Species_Eqn; w++) {
    vs[w] = mp->u_diffusivity[w][0];
    K1_gamma[w] = mp->u_diffusivity[w][2];
    K2mTg[w] = mp->u_diffusivity[w][4];
    Do[w] = mp->u_diffusivity[w][8];
    E_div_R[w] = mp->u_diffusivity[w][9];
  }

  T = fv->T;

  /* load up binary FV parameters
   * the diagonals are unity  */
  for (w = 0; w < pd->Num_Species_Eqn; w++) {
    for (w1 = w + 1; w1 < pd->Num_Species_Eqn; w1++) {
      xsi[w][w1] = mp->u_diffusivity[w][7] / mp->u_diffusivity[w1][7];
      xsi[w1][w] = 1 / xsi[w][w1];
    }
    xsi[w][w] = 1.;
    xsi[w][w2] = mp->u_diffusivity[w][7];
  }

  /* currently, c[] is assumed to be
     species_density, this needs to be explicitly
     differentiated via HKM's modifications .
     Convert density to weight fractions */

  density_tot = calc_density(mp, TRUE, densityJac, 0.0);
  memset(rho, 0, sizeof(dbl) * MAX_CONC);
  memset(drho_dc, 0, sizeof(dbl) * MAX_CONC);
  switch (mp->Species_Var_Type) {
  case SPECIES_DENSITY:
    for (w = 0; w < pd->Num_Species_Eqn; w++) {
      rho[w] = fv->c[w];
      drho_dc[w] = 1.0;
    }
    break;
  case SPECIES_UNDEFINED_FORM:
  case SPECIES_MASS_FRACTION:
    for (w = 0; w < pd->Num_Species_Eqn; w++) {
      rho[w] = fv->c[w] * density_tot;
      drho_dc[w] = density_tot + fv->c[w] * mp->d_density[MAX_VARIABLE_TYPES + w];
    }
    break;
  case SPECIES_CONCENTRATION:
    for (w = 0; w < pd->Num_Species_Eqn; w++) {
      rho[w] = fv->c[w] * mp->molecular_weight[w];
      drho_dc[w] = mp->molecular_weight[w];
    }
    break;
  case SPECIES_MOLE_FRACTION:
  case SPECIES_VOL_FRACTION:
  default:
    GOMA_EH(GOMA_ERROR, "Undefined Species formulation in Generalized_FV_Diffusivity\n");
    break;
  }

  for (w = 0; w < pd->Num_Species_Eqn; w++) {
    Do[w] *= exp(-E_div_R[w] / T);
    omega[w] = rho[w] / density_tot;
    allomega += omega[w];
  }
  omega[w2] = 1. - allomega;

  numerator[species_no] = 0.;

  for (w = 0; w < pd->Num_Species_Eqn; w++) {
    V_fh_gamma += K1_gamma[w] * omega[w] * (K2mTg[w] + T);
    numerator[species_no] += omega[w] * vs[w] * xsi[species_no][w];
    d_Vfh_dT += K1_gamma[w] * omega[w];
  }
  V_fh_gamma += K1p_gamma * omega[w2] * (K2pmTg + T);
  numerator[species_no] += omega[w2] * vsp * xsi[species_no][w2];
  d_Vfh_dT += K1p_gamma * omega[w2];

  /* mp->diffusivity is a placeholder for calculating
   * mutual diffusivity, which is a 2-dim matrix */

  mp->diffusivity[species_no] = Do[species_no] * exp(-numerator[species_no] / V_fh_gamma);

  for (w = 0; w < pd->Num_Species_Eqn; w++) {
    mp->diffusivity_gen_fick[species_no][w] =
        mp->diffusivity[species_no] *
        (fv->c[species_no] * dcoeff[species_no][w] + delta(species_no, w));
  }

  if (af->Assemble_Jacobian) {
    if (pd->v[pg->imtrx][TEMPERATURE]) {
      mp->d_diffusivity[species_no][TEMPERATURE] =
          mp->diffusivity[species_no] * d_Vfh_dT * numerator[species_no] / pow(V_fh_gamma, 2.0);

      for (w = 0; w < pd->Num_Species_Eqn; w++) {
        mp->d_diffusivity_gf[species_no][w][TEMPERATURE] =
            mp->d_diffusivity[species_no][TEMPERATURE] *
            (fv->c[species_no] * dcoeff[species_no][w] + delta(species_no, w));
      }
    }
    mp->d_diffusivity[species_no][TEMPERATURE] = 0.;

    if (pd->v[pg->imtrx][MASS_FRACTION]) {

      for (w = 0; w < pd->Num_Species_Eqn; w++) {
        d_Vfh_domega[w] = K1_gamma[w] * (K2mTg[w] + T) - K1p_gamma * (K2pmTg + T);
        d_num_domega[species_no][w] = vs[w] * xsi[species_no][w] - vsp * xsi[species_no][w2];
        d_omega_dc[species_no][w] =
            -rho[species_no] * mp->d_density[MAX_VARIABLE_TYPES + w] / SQUARE(density_tot);
      }
      d_omega_dc[species_no][species_no] += drho_dc[species_no] / density_tot;

      for (w = 0; w < pd->Num_Species_Eqn; w++) {

        mp->d_diffusivity[species_no][MAX_VARIABLE_TYPES + w] =
            mp->diffusivity[species_no] * d_omega_dc[species_no][w] *
            (-d_num_domega[species_no][w] * V_fh_gamma + d_Vfh_domega[w] * numerator[species_no]) /
            pow(V_fh_gamma, 2.0);

        for (w1 = 0; w1 < pd->Num_Species_Eqn; w1++) {

          mp->d_diffusivity_gf[species_no][w1][MAX_VARIABLE_TYPES + w] =
              (mp->d_diffusivity[species_no][MAX_VARIABLE_TYPES + w] * dcoeff[species_no][w1] +
               mp->diffusivity[species_no] * d2coeff[species_no][w1][w]) *
              fv->c[species_no];
        }
        mp->d_diffusivity_gf[species_no][species_no][MAX_VARIABLE_TYPES + w] +=
            mp->diffusivity[species_no] * dcoeff[species_no][w] +
            mp->d_diffusivity[species_no][MAX_VARIABLE_TYPES + w];

        mp->d_diffusivity[species_no][MAX_VARIABLE_TYPES + w] = 0.;
      }
    }
  }
  mp->diffusivity[species_no] = 0.;

  return (0);

} /* End of Generalized_FV_Diffusivity */

/******************************************************************************
 *
 *  Function that computes the "diffusive" mass flux resulting from
 *  velocity gradients and viscosity gradients, i.e. hydroynamics
 *  and loads st->diff_flux[][].  Returns mass flux for species w in *st->diff_flux
 *
 *     Author: T.A. Baer, 11/96
 *
 *  Added curvature-induced migration as described by Krishnan et al. (1996)
 *     Author: A.C. Sun, 7/98
 *
 ******************************************************************************/
int hydro_flux(struct Species_Conservation_Terms *st,
               int w,     /* species number */
               double tt, /* parameter to vary time integration from
                             explicit (tt = 1) to implicit (tt = 0) */
               double dt, /* current time step size */
               const double hsquared[DIM]) {
  int a, i, j, l, p, b, var;
  int status = 1;
  int dim;
  dbl gammadot, *grad_gd, gamma_dot[DIM][DIM];
  VISCOSITY_DEPENDENCE_STRUCT d_mu_struct; /* viscosity dependence */
  VISCOSITY_DEPENDENCE_STRUCT *d_mu = &d_mu_struct;
  dbl mu, grad_mu[DIM];
  dbl grad_mu1[DIM], Dmu1;
  dbl mu0;
  dbl mu_ppt = 0.0;
  dbl Dc, Dmu, Dg;
  dbl Dd[DIM];
  dbl dDd_dy[DIM];
  dbl dDd_dgrady[DIM];
  dbl dDd_dv[DIM];

  /* For shock capturing operator */
  dbl sign_adv, sign_flux;
  dbl y_dot, adv, conv_flux;
  dbl hh_siz, h_elem;

  dbl *Y, (*grad_Y)[DIM], *dmu_dY, *d2mu_dY2;
  /*  dbl d2mu_dgd_dY; */ /* cross derivative of viscosity wrt
                      shear rate and concentration */
  dbl f, maxpack, nexp, rzexp;
  dbl df_dy, df_dmu, df_dmu0 = 0;

  dbl(*d_grad_gd_dmesh)[DIM][MDE];
  dbl d_grad_mu_dmesh[DIM][DIM][MDE];
  dbl c_term = 0.0, mu_term = 0.0, g_term = 0.0, d_term = 0.0, r_term = 0.0;

  dbl dDc_dy = 0.;
  dbl dDmu_dy = 0.;
  dbl dDc_dT = 0.0; /* These are just place holders at the moment */
  dbl dDmu_dT = 0.0;
  dbl del_rho = 0.0;
  dbl dmu0_dT = 0.;

  dbl rel_mu_denom;

  /* Initialize arrays */
  for (i = 0; i < DIM; i++) {
    Dd[i] = 0;
    dDd_dy[i] = 0;
    dDd_dgrady[i] = 0;
    dDd_dv[i] = 0;
  }

  /* Jump out if we really want to do the Q-tensor diffusive flux.
   * If/when things get more complicated, it might be worth making
   * available two separate calls in the calling routines (hydro_flux
   * and hydro_qtensor_flux), instead of just calling
   * hydro_qtensor_flux from here.  The downside of putting this call
   * here, instead of outside this routine, is that we have all this
   * extra junk declared (but at least we don't waste time executing
   * unnecessary code).  The upside is that it will be transparent to
   * all of the calling routines as to exactly which diffusive flux
   * model we are implementing. -MMH
   */
  /*   if(mp->QTensorDiffType[w] == CONSTANT || */
  /*      mp->QtensorExtensionPModel == CONSTANT) */
  /*     { */
  /*       hydro_qtensor_flux(st, w); */
  /*       return(1); */
  /*     } */

  if (cr->MassFluxModel == HYDRODYNAMIC_QTENSOR_OLD) {
    hydro_qtensor_flux(st, w);
    return (1);
  } else if (cr->MassFluxModel == HYDRODYNAMIC_QTENSOR) {
    hydro_qtensor_flux_new(st, w);
    return (1);
  }

  /*
   * Put in a warning about these species variable types until
   * the equations can be worked out
   */
  if (mp->Species_Var_Type == SPECIES_MASS_FRACTION ||
      mp->Species_Var_Type == SPECIES_MOLE_FRACTION ||
      mp->Species_Var_Type == SPECIES_VOL_FRACTION) {
    GOMA_EH(GOMA_ERROR, "Possible Conflict: Hydro Flux Expression hasn't been checked out for this "
                        "species var type");
  }

  /* Set up some convenient local variables and pointers */

  Y = fv->c;
  grad_Y = fv->grad_c;
  grad_gd = fv->grad_SH;
  gammadot = fv->SH;
  d_grad_gd_dmesh = fv->d_grad_SH_dmesh;
  maxpack = gn->maxpack;
  nexp = gn->nexp;

  dim = pd->Num_Dim;

  dmu_dY = &(mp->d_viscosity[MAX_VARIABLE_TYPES]);
  d2mu_dY2 = &(mp->d2_viscosity[MAX_VARIABLE_TYPES]);
  /*  d2mu_dgd_dY = mp->d2_viscosity[MAX_VARIABLE_TYPES + MAX_CONC]; */

  /* Compute gamma_dot[][], acceleration vector accel[],
     and the norm of the velocity vector */

  /* Compute gammadot, grad(gammadot), gamma_dot[][], d_gd_dG, and d_grad_gd_dG */

  memset(gamma_dot, 0, DIM * DIM * sizeof(dbl));

  rzexp = 0.;

  for (a = 0; a < VIM; a++) {
    for (b = 0; b < VIM; b++) {
      gamma_dot[a][b] = fv->grad_v[a][b] + fv->grad_v[b][a];
    }
  }

  /* get mu and grad_mu */

  mu = viscosity(gn, gamma_dot, d_mu);

  if (gn->ConstitutiveEquation == SUSPENSION || gn->ConstitutiveEquation == CARREAU_SUSPENSION ||
      gn->ConstitutiveEquation == POWERLAW_SUSPENSION) {
    mu0 = gn->mu0; /* viscosity of pure fluid */
    rel_mu_denom = (1.0 - Y[w] / maxpack);
    if (rel_mu_denom <= 0)
      rel_mu_denom = 0.01;

    mu_ppt = pow(rel_mu_denom, nexp);

  } else if (gn->ConstitutiveEquation == FILLED_EPOXY)
  /*The extent of cure and temperature must be accounted for in
    calculating mu0 */
  {
    rel_mu_denom = (1.0 - Y[w] / maxpack);
    if (rel_mu_denom <= 0)
      rel_mu_denom = 0.01;
    mu_ppt = pow(rel_mu_denom, nexp);

    mu0 = mu / mu_ppt;
  } else {
    maxpack = .68;
    rel_mu_denom = (1.0 - Y[w] / maxpack);
    if (rel_mu_denom <= 0)
      rel_mu_denom = 0.01;
    mu0 = mu; /* viscosity of pure fluid */
    nexp = 0.;
  }

  memset(grad_mu, 0, DIM * sizeof(dbl));
  memset(grad_mu1, 0, DIM * sizeof(dbl));

  /* Controversial as to whether the grad should
   *  depend on all variables or just concentration.
   * For now we will just do concentration.
   */

  /* Separate grad_mu into two parts, so that we can
   * use two different coefficients. The first part
   * gives the sensitivity wrt to concentration and
   * the second gives the sensitivity wrt to other
   * fields such as shear-rate etc.
   * This is to get a better match with the data.
   */

  for (a = 0; a < dim; a++) {
    grad_mu[a] += dmu_dY[w] * fv->grad_c[w][a];
  }

  for (a = 0; a < dim; a++) {
    grad_mu1[a] += d_mu->gd * fv->grad_SH[a];
  }

  /* If non-neutrally bouyant suspension, compute density difference
   *  or it defaults to zero
   */

  if (mp->DensityModel == SUSPENSION) {
    del_rho = (mp->u_density[2] - mp->u_density[1]);
  }

  /* Compute HYDRODYNAMIC diffusive fluxes */
  /* Assign diffusivity values to each term */

  if (Y[w] > 0.) {
    if (mp->GamDiffType[w] == LINEAR) {
      Dc = mp->u_gadiffusivity[w][0] * 1.4 * Y[w];
      dDc_dy = mp->u_gadiffusivity[w][0] * 1.4;
    } else if (mp->GamDiffType[w] == INGBER) {
      Dc = mp->u_gadiffusivity[w][0] * Y[w] * 1.4;
      dDc_dy = mp->u_gadiffusivity[w][0];
    } else if (mp->GamDiffType[w] == LEVEL_SET) {
      double width;
      if (ls == NULL)
        GOMA_EH(GOMA_ERROR,
                "Need to activate to Level Set Interface Tracking to use this model.\n");

      width = (mp->u_gadiffusivity[w][2] == 0.0) ? ls->Length_Scale : mp->u_gadiffusivity[w][2];

      ls_transport_property(mp->u_gadiffusivity[w][0], mp->u_gadiffusivity[w][1], width, &Dc, NULL);
      dDc_dy = 0.0;
    } else {
      Dc = mp->gam_diffusivity[w];
      dDc_dy = 0.0;
    }

    if (mp->MuDiffType[w] == LINEAR) {
      Dmu = mp->u_mdiffusivity[w][0] * 1.4 * Y[w];
      dDmu_dy = mp->u_mdiffusivity[w][0] * 1.4;
    } else if (mp->MuDiffType[w] == LEVEL_SET) {
      double width;
      if (ls == NULL)
        GOMA_EH(GOMA_ERROR,
                "Need to activate to Level Set Interface Tracking to use this model.\n");

      width = (mp->u_mdiffusivity[w][2] == 0.0) ? ls->Length_Scale : mp->u_mdiffusivity[w][2];

      ls_transport_property(mp->u_mdiffusivity[w][0], mp->u_mdiffusivity[w][1], width, &Dmu, NULL);
      dDmu_dy = 0.0;
    } else {
      Dmu = mp->mu_diffusivity[w];
      dDmu_dy = 0.0;
    }

    if (mp->GravDiffType[w] == BISECTION) {
      Dg = mp->u_gdiffusivity[w][0] * del_rho;
    } else if (mp->GravDiffType[w] == RZBISECTION) {
      Dg = mp->u_gdiffusivity[w][0] * del_rho;
      rzexp = mp->u_gdiffusivity[w][1];
    } else if (mp->GravDiffType[w] == RICHARDSON_ZAKI) {
      Dg = mp->u_gdiffusivity[w][0] * del_rho;
      rzexp = mp->u_gdiffusivity[w][1];
    } else if (mp->GravDiffType[w] == LEVEL_SET) {
      double width;
      if (ls == NULL)
        GOMA_EH(GOMA_ERROR,
                "Need to activate to Level Set Interface Tracking to use this model.\n");

      width = (mp->u_gdiffusivity[w][2] == 0.0) ? ls->Length_Scale : mp->u_gdiffusivity[w][2];

      ls_transport_property(mp->u_gdiffusivity[w][0], mp->u_gdiffusivity[w][1], width, &Dg, NULL);
      Dg *= del_rho;
    } else {
      Dg = mp->g_diffusivity[w] * del_rho;
    }
  } else {
    Dg = 0.;
    Dmu = 0.;
    Dc = 0.;
  }

  /*  anisotropic diffusion coefficient set
   *  in mat file it is direction dependent
   */

  if (mp->FickDiffType[w] == ANISOTROPIC) {
    for (a = 0; a < dim; a++) {
      Dd[a] = mp->u_fdiffusivity[w][a];
      dDd_dy[a] = 0.;
    }
  } else if (mp->FickDiffType[w] == EXP_DECAY) {
    for (a = 0; a < dim; a++) {
      if (Y[w] >= 0. && Y[w] <= maxpack) {
        Dd[a] = mp->u_fdiffusivity[w][0] * (exp(-mp->u_fdiffusivity[w][1] * Y[w]) +
                                            exp(-mp->u_fdiffusivity[w][1] * fabs(maxpack - Y[w])));
        dDd_dy[a] = -mp->u_fdiffusivity[w][0] * mp->u_fdiffusivity[w][1] *
                    (exp(-mp->u_fdiffusivity[w][1] * Y[w]) -
                     exp(-mp->u_fdiffusivity[w][1] * fabs(maxpack - Y[w])));
      } else if (Y[w] < 0.) {
        Dd[a] = mp->u_fdiffusivity[w][0];
        dDd_dy[a] = 0.;
      } else if (Y[w] > maxpack) {
        Dd[a] = mp->u_fdiffusivity[w][0];
        dDd_dy[a] = 0.;
      }
    }
  } else if (mp->FickDiffType[w] == SHOCK) {
    if (pd->TimeIntegration != STEADY)
      y_dot = fv_dot->c[w];
    else
      y_dot = 0.0;

    adv = 0.;
    for (a = 0; a < dim; a++) {
      adv += fv->v[a] * grad_Y[w][a];
    }

    if (adv < 0)
      sign_adv = -1;
    else
      sign_adv = 1;

    conv_flux = y_dot + adv;
    if (conv_flux < 0)
      sign_flux = -1;
    else
      sign_flux = 1;

    hh_siz = 0.;
    for (a = 0; a < dim; a++) {
      hh_siz += hsquared[a];
    }
    /* This is the average value of h**2 in the element */
    hh_siz = hh_siz / ((double)dim);

    /* This is the size of the element */
    h_elem = sqrt(hh_siz);

    for (a = 0; a < dim; a++) {
      Dd[a] =
          mp->u_fdiffusivity[w][0] * h_elem * sign_flux * (y_dot + adv) / (sign_adv * adv + h_elem);

      if (pd->TimeIntegration != STEADY)
        dDd_dy[a] = mp->u_fdiffusivity[w][0] * h_elem * sign_flux * (1 + 2. * tt) / dt /
                    (sign_adv * adv + h_elem);
      else
        dDd_dy[a] = 0.;

      dDd_dgrady[a] = mp->u_fdiffusivity[w][0] * h_elem * sign_flux * (fv->v[a]) /
                      (sign_adv * adv + h_elem) * (1. - (y_dot + adv) / (sign_adv * adv + h_elem));
      dDd_dv[a] = mp->u_fdiffusivity[w][0] * h_elem * sign_flux * (grad_Y[w][a]) /
                  (sign_adv * adv + h_elem) * (1. - (y_dot + adv) / (sign_adv * adv + h_elem));
    }
  }

  /* The form of hindered settling function is dependent
     upon the model_name. If CONSTANT, it
     is a monotonically decreasing function
     as Y[w] changes.  If RICHARDSON-ZAKI, it follows
     the R-Z formula.  Also, a BISECTION version for
     each formula is coded. ACSun 9/98 CAR 10/98*/

  if (Y[w] >= maxpack)
    Y[w] = maxpack;
  if (mp->GravDiffType[w] == RICHARDSON_ZAKI) {
    f = pow(1. - Y[w], rzexp);
    df_dmu = 0.;
    df_dmu0 = 0.;
    df_dy = -rzexp * f / (1. - Y[w]);
  } else {
    f = (1. - mp->reference_concn[w]) / mu_ppt;
    /*    f = (1. - Y[w])/mu_ppt;  */
    df_dmu = 0;
    df_dmu0 = 0.0;
    df_dy = f * nexp / maxpack / rel_mu_denom;
    /*     df_dy = -1./mu_ppt + f*nexp/maxpack/rel_mu_denom;*/
  }

  /* Hardwire this here - later move to input deck if it works! */
  Dmu1 = 0.;

  /* assemble residual */
  for (a = 0; a < dim; a++) {
    st->diff_flux[w][a] = -Y[w] * Dc * (Y[w] * grad_gd[a] + gammadot * grad_Y[w][a]);
    st->diff_flux[w][a] += -Y[w] * Y[w] * gammadot * (Dmu * grad_mu[a] + Dmu1 * grad_mu1[a]) / mu;

    st->diff_flux[w][a] += (Dg * f * Y[w]) * mp->momentum_source[a] / mu0;
    st->diff_flux[w][a] += -Dd[a] * grad_Y[w][a];
  }

  if (af->Assemble_Jacobian) {
    /* Compute derivative of viscosity gradient with respect to mesh displacement */
    if (pd->v[pg->imtrx][MESH_DISPLACEMENT1]) {
      memset(d_grad_mu_dmesh, 0, DIM * DIM * MDE * sizeof(dbl));

      /* d_grad_mu_dmesh */

      for (a = 0; a < dim; a++) {
        for (p = 0; p < dim; p++) {
          for (j = 0; j < ei[pg->imtrx]->dof[MESH_DISPLACEMENT1]; j++) {
            d_grad_mu_dmesh[a][p][j] = 0.0;

            var = MASS_FRACTION;
            for (l = 0; l < ei[pg->imtrx]->dof[var]; l++) {
              d_grad_mu_dmesh[a][p][j] +=
                  dmu_dY[w] * (*esp->c[w][l]) * bf[var]->d_grad_phi_dmesh[l][a][p][j];
            }

            d_grad_mu_dmesh[a][p][j] += d_mu->gd * d_grad_gd_dmesh[a][p][j];
          }
        }
      }
    }

    var = MASS_FRACTION;

    for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
      for (a = 0; a < dim && pd->v[pg->imtrx][var]; a++) {

        c_term = -2.0 * Y[w] * grad_gd[a] * bf[var]->phi[j];
        c_term += -gammadot * grad_Y[w][a] * bf[var]->phi[j];
        c_term += -Y[w] * gammadot * bf[var]->grad_phi[j][a];
        c_term *= Dc;
        c_term -= dDc_dy * bf[var]->phi[j] * Y[w] * (Y[w] * grad_gd[a] + gammadot * grad_Y[w][a]);

        mu_term =
            -2.0 * Y[w] * gammadot * (Dmu * grad_mu[a] + Dmu1 * grad_mu1[a]) * bf[var]->phi[j] / mu;
        mu_term += Y[w] * Y[w] * gammadot * d_mu->C[w][j] *
                   (Dmu * grad_mu[a] + Dmu1 * grad_mu1[a]) / mu / mu;
        mu_term +=
            -Dmu * Y[w] * Y[w] * gammadot *
            (dmu_dY[w] * bf[var]->grad_phi[j][a] + d2mu_dY2[w] * grad_Y[w][a] * bf[var]->phi[j]) /
            mu;

        /* 	      mu_term +=
         * -Dmu*Y[w]*Y[w]*gammadot*(d2mu_dgd_dY*bf[var]->phi[j]*grad_Y[w][a])/mu; */
        /* 	      mu_term +=
         * -Dmu1*Y[w]*Y[w]*gammadot*(d2mu_dgd_dY*bf[var]->phi[j]*fv->grad_SH[a])/mu; */

        mu_term -= dDmu_dy * bf[var]->phi[j] * Y[w] * Y[w] * gammadot * grad_mu[a] / mu;

        g_term = ((f + df_dy * Y[w]) * bf[var]->phi[j] + Y[w] * df_dmu * d_mu->C[w][j]);
        g_term *= Dg * mp->momentum_source[a] / mu0;

        d_term = -bf[var]->grad_phi[j][a] * Dd[a] - grad_Y[w][a] * dDd_dy[a] * bf[var]->phi[j] -
                 grad_Y[w][a] * dDd_dgrady[a] * bf[var]->grad_phi[j][a];

        st->d_diff_flux_dc[w][a][w][j] = c_term + mu_term + g_term + d_term;

        /* if filled_epoxy is used, there is a dependency of mu0 to
           the cure species */
        if (gn->ConstitutiveEquation == FILLED_EPOXY) {
          st->d_diff_flux_dc[w][a][gn->cure_species_no][j] =
              Y[w] * Y[w] * gammadot * d_mu->C[gn->cure_species_no][j] *
              (Dmu * grad_mu[a] + Dmu1 * grad_mu1[a]) / mu / mu;

          st->d_diff_flux_dc[w][a][gn->cure_species_no][j] +=
              (-(Dg * f * Y[w] * mp->momentum_source[a]) / mu0) * d_mu->C[gn->cure_species_no][j] /
              mu_ppt / mu0;
        }
      }
    }

    var = SHEAR_RATE;

    for (a = 0; a < dim && pd->v[pg->imtrx][var]; a++) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {

        c_term = Y[w] * Y[w] * bf[var]->grad_phi[j][a];
        c_term += Y[w] * grad_Y[w][a] * bf[var]->phi[j];
        c_term *= -Dc;

        mu_term = -Dmu * Y[w] * Y[w] * grad_mu[a] * bf[var]->phi[j] / mu;

        mu_term -= Dmu1 * Y[w] * Y[w] * grad_mu1[a] * bf[var]->phi[j] / mu;

        mu_term += (Y[w] * Y[w] * gammadot * d_mu->gd * bf[var]->phi[j] *
                    (Dmu1 * grad_mu1[a] + Dmu * grad_mu[a]) / mu / mu);

        mu_term -= Dmu1 * (Y[w] * Y[w] * gammadot *
                           (d_mu->gd * bf[var]->grad_phi[j][a] +
                            mp->d2_viscosity[SHEAR_RATE] * fv->grad_SH[a] * bf[var]->phi[j]) /
                           mu);

        st->d_diff_flux_dSH[w][a][j] = c_term + mu_term;
      }
    }

    var = TEMPERATURE;
    for (a = 0; a < dim && pd->v[pg->imtrx][var]; a++) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++)

      {
        c_term = dDc_dT * (Y[w] * Y[w] * grad_gd[a] + Y[w] * gammadot * grad_Y[w][a]);
        mu_term = (dDmu_dT * grad_mu[a] - Dmu * grad_mu[a] * d_mu->T[j] / mu);
        mu_term *= -Y[w] * Y[w] * gammadot / mu;
        g_term = df_dmu * d_mu->T[j] + (df_dmu0 - f / mu0) * dmu0_dT * bf[var]->phi[j];
        g_term *= (Dg * Y[w] * mp->momentum_source[a]) / mu0;

        st->d_diff_flux_dT[w][a][j] = c_term + mu_term + g_term;
      }
    }

    /* Compute velocity derivatives of normal acceleration
       vector, the velocity norm.  Did not observe any
       anomalies, but the Jacobian is slightly off when
       debug=-2 */

    var = VELOCITY1;
    memset(st->d_diff_flux_dv, 0, MAX_CONC * DIM * DIM * MDE * sizeof(dbl));

    if (pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        for (a = 0; a < dim && pd->v[pg->imtrx][var]; a++) {
          d_term = -grad_Y[w][a] * dDd_dv[a] * bf[var]->phi[j];
          st->d_diff_flux_dv[w][a][a][j] = c_term + mu_term + g_term + d_term;
        }
      }
    }

    var = MESH_DISPLACEMENT1;

    if (pd->v[pg->imtrx][var]) {
      for (a = 0; a < dim; a++) {
        for (p = 0; p < dim; p++) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++)

          {
            c_term = -Y[w] * Y[w] * d_grad_gd_dmesh[a][p][j];
            c_term += -Y[w] * gammadot * fv->d_grad_c_dmesh[a][w][p][j];
            c_term *= Dc;

            mu_term = -Y[w] * Y[w] * gammadot * d_grad_mu_dmesh[a][p][j] / mu;
            mu_term *= Dmu;

            d_term = -fv->d_grad_c_dmesh[a][w][p][j] * Dd[a];
            r_term = 0.;

            st->d_diff_flux_dmesh[w][a][p][j] = c_term + mu_term + d_term + r_term;
          }
        }
      }
    }
  }

  return (status);
}

/*****************************************************************************/
/*****************************************************************************/
/******************************************************************************
 *
 *  suspension_balance
 *
 *  Function that computes the mass flux resulting from
 *  velocity gradients and viscosity gradients for the
 *  suspension balance model (Nott & Brady with improvements
 *  from Morris and Boulay), i.e. hydroynamics
 *  and loads st->diff_flux[][].
 *
 *  Returns mass flux for species w in *st->diff_flux
 *
 *     Author: R.R. Rao, 1/2010
 *
 *
 ******************************************************************************/

int suspension_balance(struct Species_Conservation_Terms *st, int w) /* species number */
{
  int a, j, p, b, var;
  int status = 1;
  int dim;
  dbl gamma_dot[DIM][DIM];
  VISCOSITY_DEPENDENCE_STRUCT d_mu_struct; /* viscosity dependence */
  VISCOSITY_DEPENDENCE_STRUCT *d_mu = &d_mu_struct;
  dbl mu;
  dbl mu0;
  dbl Dg;
  dbl Dd[DIM];
  dbl dDd_dy[DIM];
  dbl *Y, (*grad_Y)[DIM];

  dbl f, maxpack, rzexp = 0;
  dbl df_dy, df_dmu;
  dbl M; /* hindrance function */
  dbl dM_dy, dM_dmu, gammadot = 0.;

  dbl c_term, mu_term, g_term, d_term;

  dbl div_tau_p[DIM];                     /* divergence of the particle stress*/
  dbl d_div_tau_p_dgd[DIM][MDE];          /* derivative wrt shear_rate_invariant */
  dbl d_div_tau_p_dy[DIM][MAX_CONC][MDE]; /* derivative wrt concentration */
  dbl d_div_tau_p_dv[DIM][DIM][MDE];      /* derivative wrt velocity */
  dbl d_div_tau_p_dmesh[DIM][DIM][MDE];   /* derivative wrt mesh */
  dbl d_div_tau_p_dp[DIM][MDE];           /* derivative wrt pressure */

  dbl d_gd_dv[DIM][MDE];    /* derivative of strain rate invariant
                               wrt velocity */
  dbl d_gd_dmesh[DIM][MDE]; /* derivative of strain rate invariant
                               wrt mesh */

  dbl df_dmu0 = 0.0, dmu0_dcure = 0.0, dmu0_dT = 0.0;
  dbl del_rho = 0.0;

  /* Set up some convenient local variables and pointers */
  Y = fv->c;
  grad_Y = fv->grad_c;

  dim = pd->Num_Dim;

  /* Compute gamma_dot[][] */

  /* Compute gammadot, grad(gammadot), gamma_dot[][], d_gd_dG, and d_grad_gd_dG */

  memset(gamma_dot, 0, DIM * DIM * sizeof(dbl));

  for (a = 0; a < VIM; a++) {
    for (b = 0; b < VIM; b++) {
      gamma_dot[a][b] = fv->grad_v[a][b] + fv->grad_v[b][a];
    }
  }

  calc_shearrate(&gammadot, gamma_dot, d_gd_dv, d_gd_dmesh);

  mu = viscosity(gn, gamma_dot, d_mu);

  if (gn->ConstitutiveEquation == SUSPENSION || gn->ConstitutiveEquation == CARREAU_SUSPENSION ||
      gn->ConstitutiveEquation == POWERLAW_SUSPENSION || gn->ConstitutiveEquation == FILLED_EPOXY) {
    maxpack = gn->maxpack;
    mu0 = gn->mu0; /* viscosity of pure fluid */
  } else {
    maxpack = .68;
    mu0 = mu; /* viscosity of pure fluid */
  }

  /* If non-neutrally bouyant suspension, compute density difference
   *  or it defaults to zero
   */

  if (mp->DensityModel == SUSPENSION) {
    del_rho = (mp->u_density[2] - mp->u_density[1]);
  }

  /* Compute HYDRODYNAMIC diffusive fluxes */
  /* Assign diffusivity values to each term */

  if (Y[w] > 0.) {
    if (mp->GravDiffType[w] == RICHARDSON_ZAKI) {
      Dg = mp->u_gdiffusivity[w][0];
      rzexp = mp->u_gdiffusivity[w][1];
    } else {
      Dg = mp->g_diffusivity[w];
    }
  } else {
    Dg = 0.;
  }

  /*  anisotropic diffusion coefficient set
   *  in mat file it is direction dependent
   */

  if (mp->FickDiffType[w] == ANISOTROPIC) {
    for (a = 0; a < dim; a++) {
      Dd[a] = mp->u_fdiffusivity[w][a];
      dDd_dy[a] = 0.;
    }
  } else if (mp->FickDiffType[w] == EXP_DECAY) {
    for (a = 0; a < dim; a++) {
      Dd[a] = mp->u_fdiffusivity[w][0] * exp(-mp->u_fdiffusivity[w][1] * Y[w]);
      dDd_dy[a] = -mp->u_fdiffusivity[w][1] * Dd[a];
    }
  }

  /* The form of hindered settling function is dependent
     upon the model_name. If CONSTANT, it
     is a monotonically decreasing function
     as Y[w] changes.  If RICHARDSON-ZAKI, it follows
     the R-Z formula.  */

  if (Y[w] >= maxpack)
    Y[w] = maxpack;
  if (mp->GravDiffType[w] == RICHARDSON_ZAKI) {
    if (Y[w] / maxpack < 0.95) {
      f = pow(1. - Y[w], rzexp) / mu0;
      f *= (1. - Y[w] / maxpack);
      // df_dmu =  -f/mu;
      df_dmu = 0.;
      df_dy = -rzexp * f / (1. - Y[w]);
      df_dy += -pow(1. - Y[w], rzexp) / (mu0 * maxpack);
    } else {
      f = 0.;
      df_dmu = 0.;
      df_dy = 0.;
    }
  } else {
    f = (1. - Y[w]) / mu;
    df_dmu = -f / mu;
    df_dy = -1. / mu;
  }

  memset(div_tau_p, 0, sizeof(double) * DIM);
  memset(d_div_tau_p_dgd, 0, sizeof(double) * DIM * MDE);
  memset(d_div_tau_p_dy, 0, sizeof(double) * DIM * MAX_CONC * MDE);
  memset(d_div_tau_p_dv, 0, sizeof(double) * DIM * DIM * MDE);
  memset(d_div_tau_p_dmesh, 0, sizeof(double) * DIM * DIM * MDE);
  memset(d_div_tau_p_dp, 0, sizeof(double) * DIM * MDE);

  /* This is the divergence of the particle stress  */
  divergence_particle_stress(div_tau_p, d_div_tau_p_dgd, d_div_tau_p_dy, d_div_tau_p_dv,
                             d_div_tau_p_dmesh, d_div_tau_p_dp, w);

  /* this is the hindered settling term that modifies the flux */
  M = Dg * f;
  dM_dy = Dg * df_dy;
  // dM_dmu = Dg * df_dmu;
  dM_dmu = 0.;

  /* assemble residual */
  for (a = 0; a < dim; a++) {
    st->diff_flux[w][a] = -M * div_tau_p[a];
    st->diff_flux[w][a] += M * Y[w] * mp->momentum_source[a] * del_rho;
    st->diff_flux[w][a] += -Dd[a] * grad_Y[w][a];
  }

  if (af->Assemble_Jacobian) {

    var = MASS_FRACTION;
    for (a = 0; a < dim && pd->v[pg->imtrx][var]; a++) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {

        c_term = -dM_dy * bf[var]->phi[j] * div_tau_p[a];

        c_term += -M * d_div_tau_p_dy[a][w][j];

        mu_term = -dM_dmu * d_mu->C[w][j] * div_tau_p[a];

        g_term = ((f + df_dy * Y[w]) * bf[var]->phi[j] + Y[w] * df_dmu * d_mu->C[w][j]);
        g_term *= Dg * mp->momentum_source[a] * del_rho;

        d_term = -bf[var]->grad_phi[j][a] * Dd[a] - grad_Y[w][a] * dDd_dy[a] * bf[var]->phi[j];

        st->d_diff_flux_dc[w][a][w][j] = c_term + mu_term + g_term + d_term;
      }

      /* if filled_epoxy is used, there is a dependency of viscosity on
         the cure species */

      if (gn->ConstitutiveEquation == FILLED_EPOXY) {
        for (a = 0; a < dim && pd->v[pg->imtrx][var]; a++) {

          st->d_diff_flux_dc[w][a][gn->cure_species_no][j] =
              (df_dmu * d_mu->C[gn->cure_species_no][j]) +
              (df_dmu0 - f / mu0) * dmu0_dcure * bf[var]->phi[j];

          st->d_diff_flux_dc[w][a][gn->cure_species_no][j] *=
              (Dg * Y[w] * mp->momentum_source[a] * del_rho) / mu0;
        }
      }
    }

    var = PRESSURE;
    for (a = 0; a < dim && pd->v[pg->imtrx][var]; a++) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        c_term = 0.;

        c_term += -M * d_div_tau_p_dp[a][j];

        st->d_diff_flux_dP[w][a][j] = c_term;
      }
    }

    var = SHEAR_RATE;
    for (a = 0; a < dim && pd->v[pg->imtrx][var]; a++) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {

        c_term = -M * d_div_tau_p_dgd[a][j];

        st->d_diff_flux_dSH[w][a][j] = c_term;
      }
    }

    var = TEMPERATURE;
    for (a = 0; a < dim && pd->v[pg->imtrx][var]; a++) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        c_term = -Y[w] * dM_dmu * d_mu->T[j] * div_tau_p[a];

        mu_term = 0.;

        g_term = df_dmu * d_mu->T[j] + (df_dmu0 - f / mu0) * dmu0_dT * bf[var]->phi[j];
        g_term *= (Dg * Y[w] * mp->momentum_source[a] * del_rho) / mu0;

        st->d_diff_flux_dT[w][a][j] = c_term + mu_term + g_term;
      }
    }

    var = VELOCITY1;
    memset(st->d_diff_flux_dv, 0, MAX_CONC * DIM * DIM * MDE * sizeof(dbl));

    if (pd->v[pg->imtrx][var]) {
      for (a = 0; a < VIM; a++) {
        for (p = 0; p < VIM; p++) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            c_term = -M * d_div_tau_p_dv[a][p][j];

            mu_term = 0.;

            st->d_diff_flux_dv[w][a][p][j] = c_term + mu_term;
          }
        }
      }
    }

    var = MESH_DISPLACEMENT1;
    memset(st->d_diff_flux_dmesh, 0, MAX_CONC * DIM * DIM * MDE * sizeof(dbl));

    if (pd->v[pg->imtrx][var]) {
      for (a = 0; a < dim; a++) {
        for (p = 0; p < dim; p++) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++)

          {
            c_term = -M * d_div_tau_p_dmesh[a][p][j];

            mu_term = 0.;

            d_term = -fv->d_grad_c_dmesh[a][w][p][j] * Dd[a];

            st->d_diff_flux_dmesh[w][a][p][j] = c_term + mu_term + d_term;
          }
        }
      }
    }

  } /* End Jacobian */

  return (status);
}

/******************************************************************************
 *
 *  Routine which calculates the particle stress for the
 *  suspension balance model
 *     Author: R.R. Rao, 6/11
 *
 ******************************************************************************/

int particle_stress(dbl tau_p[DIM][DIM],                     /* particle stress */
                    dbl d_tau_p_dv[DIM][DIM][DIM][MDE],      /* derivative wrt velocity */
                    dbl d_tau_p_dvd[DIM][DIM][DIM][MDE],     /* derivative wrt vorticity dir */
                    dbl d_tau_p_dy[DIM][DIM][MAX_CONC][MDE], /* derivative wrt concentration */
                    dbl d_tau_p_dmesh[DIM][DIM][DIM][MDE],   /* derivative wrt mesh */
                    dbl d_tau_p_dp[DIM][DIM][MDE],           /* derivative wrt pressure */
                    int w)                                   /* species number */
{
  /*local variables */
  int a, b, p, q, var, j, dofs;
  int status = 1;
  int dim;
  dbl gammadot, gamma_dot[DIM][DIM];

  dbl d_gd_dv[DIM][MDE];    /* derivative of strain rate invariant
                               wrt velocity */
  dbl d_gd_dmesh[DIM][MDE]; /* derivative of strain rate invariant
                               wrt mesh */
  dbl mu0;
  dbl *Y;

  dbl mu;
  dbl Kn;

  dbl pp = 0, d_pp_dy = 0;
  dbl y_norm, comp = 0, comp1;

  dbl maxpack;

  /* Q tensor info */
  dbl qtensor[DIM][DIM]; /* I - 1/2 v^t v at a gauss point */
  dbl d_qtensor_dvd[DIM][DIM][DIM][MDE];
  dbl vort_dir_local[DIM], tmp;
  int print = 0;

  dbl v1[DIM], v2[DIM], v3[DIM];

  memset(v1, 0, DIM * sizeof(double));
  memset(v2, 0, DIM * sizeof(double));
  memset(v3, 0, DIM * sizeof(double));

  Y = fv->c;
  dim = pd->Num_Dim;

  /* Compute gamma_dot[][] */

  /* Compute gammadot, grad(gammadot), gamma_dot[][], d_gd_dG, and d_grad_gd_dG */

  for (a = 0; a < VIM; a++) {
    for (b = 0; b < VIM; b++) {
      gamma_dot[a][b] = fv->grad_v[a][b] + fv->grad_v[b][a];
    }
  }

  /* This is the shear rate based on velocity */
  calc_shearrate(&gammadot, gamma_dot, d_gd_dv, d_gd_dmesh);

  /* get mu and grad_mu */

  mu = viscosity(gn, gamma_dot, NULL);

  if (gn->ConstitutiveEquation == SUSPENSION || gn->ConstitutiveEquation == CARREAU_SUSPENSION ||
      gn->ConstitutiveEquation == POWERLAW_SUSPENSION || gn->ConstitutiveEquation == FILLED_EPOXY) {
    maxpack = gn->maxpack;
    mu0 = gn->mu0; /* viscosity of pure fluid */
  } else {
    maxpack = .68;
    mu0 = mu; /* viscosity of pure fluid */
  }

  y_norm = Y[w] / maxpack;
  if ((y_norm < 0.95) && (y_norm > 0.)) {
    comp = pow((1. - y_norm), -2.);
    comp1 = 2. / maxpack * pow((1. - y_norm), -3.);
  } else if (y_norm >= 0.95) {
    comp = pow((0.05), -2.);
    comp1 = 0.;
  } else if (y_norm <= 0.) {
    comp = 1.;
    comp1 = 0.;
  }

  /* Migrate this input deck later */
  Kn = 0.75;

  /* This is the particle pressure */
  if ((y_norm > 0.0) && (y_norm < 0.95)) {
    pp = Kn * y_norm * y_norm * comp;
    d_pp_dy = 2. * Kn * y_norm / maxpack * comp + Kn * y_norm * y_norm * comp1;
  } else if (y_norm <= 0.) {
    pp = 0.;
    d_pp_dy = 0.;
  } else if (y_norm >= 0.95) {
    /* do some clipping near max packing so it doesn't go unstable */
    pp = Kn * 0.95 * 0.95 * comp;
    d_pp_dy = 0.;
  }
  memset(tau_p, 0, DIM * DIM * sizeof(dbl));

  if (cr->MassFluxModel == HYDRODYNAMIC_QTENSOR_OLD) {
    /* Get Q tensor */
    for (a = 0; a < VIM; a++) {
      vort_dir_local[a] = 0.0;
    }
    find_super_special_eigenvector(gamma_dot, vort_dir_local, v1, v2, v3, &tmp, print);

    for (a = 0; a < VIM; a++) {
      for (b = 0; b < VIM; b++) {
        qtensor[a][b] = (dbl)delta(a, b) - 0.5 * vort_dir_local[a] * vort_dir_local[b];
      }
    }
  } else if (cr->MassFluxModel == HYDRODYNAMIC_QTENSOR) {
    /* Get Q tensor */
    for (a = 0; a < VIM; a++) {
      for (b = 0; b < VIM; b++) {
        qtensor[a][b] = (dbl)delta(a, b) - 0.5 * fv->vd[a] * fv->vd[b];
      }
    }
  } else {
    /* assume a diagonal Q tensor */
    memset(qtensor, 0, DIM * DIM * sizeof(dbl));
    for (a = 0; a < VIM; a++) {
      qtensor[a][a] = mp->u_qdiffusivity[w][a];
    }
  }

  dbl cp;
  cp = 0.;
  memset(tau_p, 0, DIM * DIM * sizeof(dbl));
  for (a = 0; a < VIM; a++) {
    tau_p[a][a] +=
        cp * Y[w] * fv->P; /* suspension pressure adds an isotropic term on the diagonal */
    for (b = 0; b < VIM; b++) {
      tau_p[a][b] += mu0 * pp * gammadot * qtensor[a][b];
    }
  }

  memset(d_tau_p_dv, 0, DIM * DIM * DIM * MDE * sizeof(dbl));
  memset(d_tau_p_dvd, 0, DIM * DIM * DIM * MDE * sizeof(dbl));
  memset(d_tau_p_dy, 0, DIM * DIM * MAX_CONC * MDE * sizeof(dbl));
  memset(d_tau_p_dmesh, 0, DIM * DIM * DIM * MDE * sizeof(dbl));
  memset(d_tau_p_dp, 0, DIM * DIM * MDE * sizeof(dbl));

  if (af->Assemble_Jacobian) {

    var = VELOCITY1;
    if (pd->e[pg->imtrx][var]) {
      dofs = ei[pg->imtrx]->dof[var];
      if (gammadot != 0.0) {
        for (a = 0; a < VIM; a++) {
          for (b = 0; b < VIM; b++) {
            for (p = 0; p < VIM; p++) {
              for (j = 0; j < dofs; j++) {
                d_tau_p_dv[a][b][p][j] = mu0 * pp * d_gd_dv[p][j] * qtensor[a][b];
              }
            }
          }
        }
      }
    }

    var = MESH_DISPLACEMENT1;
    if (pd->e[pg->imtrx][var]) {
      dofs = ei[pg->imtrx]->dof[var];
      if (gammadot != 0.0) {
        for (a = 0; a < VIM; a++) {
          for (b = 0; b < VIM; b++) {
            for (p = 0; p < dim; p++) {
              for (j = 0; j < dofs; j++) {
                d_tau_p_dmesh[a][b][p][j] = mu0 * pp * d_gd_dmesh[p][j] * qtensor[a][b];
              }
            }
          }
        }
      }
    }

    var = VORT_DIR1;
    if (pd->e[pg->imtrx][var]) {

      memset(d_qtensor_dvd, 0, DIM * DIM * DIM * MDE * sizeof(dbl));
      for (p = 0; p < DIM; p++) {
        for (q = 0; q < DIM; q++) {
          for (b = 0; b < DIM; b++) {
            for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
              d_qtensor_dvd[p][q][b][j] =
                  -0.5 * bf[var]->phi[j] * (fv->vd[q] * delta(p, b) + fv->vd[p] * delta(q, b));
            }
          }
        }
      }

      dofs = ei[pg->imtrx]->dof[var];
      if (gammadot != 0.0) {
        for (a = 0; a < VIM; a++) {
          for (b = 0; b < VIM; b++) {
            for (p = 0; p < VIM; p++) {
              for (j = 0; j < dofs; j++) {
                d_tau_p_dvd[a][b][p][j] = mu0 * pp * gammadot * d_qtensor_dvd[a][b][p][j];
              }
            }
          }
        }
      }
    }

    var = MASS_FRACTION;
    if (pd->e[pg->imtrx][var]) {
      dofs = ei[pg->imtrx]->dof[var];
      for (a = 0; a < VIM; a++) {
        for (j = 0; j < dofs; j++) {
          d_tau_p_dy[a][a][w][j] += cp * bf[var]->phi[j] * fv->P;

          for (b = 0; b < VIM; b++) {
            d_tau_p_dy[a][b][w][j] += mu0 * d_pp_dy * gammadot * bf[var]->phi[j] * qtensor[a][b];
          }
        }
      }
    }

    var = PRESSURE;
    if (pd->e[pg->imtrx][var]) {
      dofs = ei[pg->imtrx]->dof[var];
      for (a = 0; a < VIM; a++) {
        for (j = 0; j < dofs; j++) {
          d_tau_p_dp[a][a][j] += cp * Y[w] * bf[var]->phi[j];
        }
      }
    }
  }

  return (status);
}

/******************************************************************************
 *
 *  Routine which calculates the divergence of the particle stress for the
 *  suspension balance model. This term is needed for the particle flux.
 *     Author: R.R. Rao, 6/11
 *  Input is the species number and output is the divergence of the
 *  particle stress including derivative wrt velocity gradient, concentration,
 *  and mesh.
 *
 ******************************************************************************/

int divergence_particle_stress(
    dbl div_tau_p[DIM],                     /* divergence of the stress*/
    dbl d_div_tau_p_dgd[DIM][MDE],          /* derivative wrt shear rate inv. */
    dbl d_div_tau_p_dy[DIM][MAX_CONC][MDE], /* derivative wrt concentration */
    dbl d_div_tau_p_dv[DIM][DIM][MDE],      /* derivative wrt velocity */
    dbl d_div_tau_p_dmesh[DIM][DIM][MDE],   /* derivative wrt mesh */
    dbl d_div_tau_p_dp[DIM][MDE],           /* derivative wrt pressure */
    int w)                                  /* species number */
{
  /*local variables */
  int a, j, p, b, var;
  int status = 1;
  int dim, dofs;

  dbl gammadot, *grad_gd, gamma_dot[DIM][DIM];
  dbl mu0;
  dbl *Y, (*grad_Y)[DIM];

  dbl mu;
  dbl grad_mu[DIM];
  dbl *dmu_dY;

  VISCOSITY_DEPENDENCE_STRUCT d_mu_struct; /* viscosity dependence */
  VISCOSITY_DEPENDENCE_STRUCT *d_mu = &d_mu_struct;

  dbl qtensor[DIM][DIM];

  dbl maxpack, maxpack2, Kn;
  dbl pp, d_pp_dy, d_pp2_dy2;
  dbl comp, comp1, comp2 = 0, y_norm;
  dbl gamma_nl;

  dbl vort_dir_local[DIM];
  int print = 0;
  dbl v1[DIM], v2[DIM], v3[DIM], tmp = 0.;
  dbl Q_prime[DIM][DIM], R[DIM][DIM];
  dbl radius_p, L_char, U_max;
  dbl v_bias[DIM], vort_bias[DIM], vy_bias[DIM];

  Y = fv->c;
  grad_Y = fv->grad_c;
  grad_gd = fv->grad_SH;
  gammadot = fv->SH;

  if (gammadot < 1.e-10) {
    gammadot = 1.e-10;
  }

  dim = pd->Num_Dim;

  /* Compute gammadot, grad(gammadot), gamma_dot[][], d_gd_dG, and d_grad_gd_dG */

  memset(gamma_dot, 0, DIM * DIM * sizeof(dbl));

  for (a = 0; a < VIM; a++) {
    for (b = 0; b < VIM; b++) {
      gamma_dot[a][b] = fv_old->grad_v[a][b] + fv_old->grad_v[b][a];
    }
  }

  /* get mu and grad_mu */
  mu = viscosity(gn, gamma_dot, d_mu);
  memset(grad_mu, 0, DIM * sizeof(dbl));

  dmu_dY = &(mp->d_viscosity[MAX_VARIABLE_TYPES]);
  for (a = 0; a < DIM; a++) {
    grad_mu[a] += dmu_dY[w] * fv->grad_c[w][a];
  }

  /* assume a diagonal Q tensor */
  memset(qtensor, 0, DIM * DIM * sizeof(dbl));
  for (a = 0; a < DIM; a++) {
    qtensor[a][a] = mp->u_qdiffusivity[w][a];
  }

  /* Solve for the eigenvalues of gamma_dot   */

  for (a = 0; a < VIM; a++) {
    vort_dir_local[a] = 0.0;
  }
  find_super_special_eigenvector(gamma_dot, vort_dir_local, v1, v2, v3, &tmp, print);

  memset(v_bias, 0, DIM * sizeof(dbl));
  memset(vort_bias, 0, DIM * sizeof(dbl));
  memset(vy_bias, 0, DIM * sizeof(dbl));

  v_bias[0] = 1.;
  vy_bias[1] = 1.;
  vort_bias[2] = 1.;

  bias_eigenvector_to(v1, v_bias);
  bias_eigenvector_to(v2, vy_bias);

  v3[0] = v2[1] * v1[2] - v2[2] * v1[1];
  v3[1] = v2[2] * v1[0] - v2[0] * v1[2];
  v3[2] = v2[0] * v1[1] - v2[1] * v1[0];

  memset(Q_prime, 0, sizeof(dbl) * DIM * DIM);

  Q_prime[0][0] = (qtensor[0][0] + qtensor[1][1]) / 2.;
  Q_prime[0][2] = (-qtensor[0][0] + qtensor[1][1]) / 2.;
  Q_prime[1][1] = qtensor[2][2];
  Q_prime[2][0] = (-qtensor[0][0] + qtensor[1][1]) / 2.;
  Q_prime[2][2] = (qtensor[0][0] + qtensor[1][1]) / 2.;

  memset(R, 0, DIM * DIM * sizeof(dbl));

  for (a = 0; a < DIM; a++) {
    R[a][0] = v2[a];
    R[a][1] = v3[a];
    R[a][2] = v1[a];
  }

  memset(qtensor, 0, DIM * DIM * sizeof(dbl));
  rotate_tensor(Q_prime, qtensor, R, 0);

  if (gn->ConstitutiveEquation == SUSPENSION || gn->ConstitutiveEquation == CARREAU_SUSPENSION ||
      gn->ConstitutiveEquation == POWERLAW_SUSPENSION || gn->ConstitutiveEquation == FILLED_EPOXY) {
    maxpack = gn->maxpack;
    mu0 = gn->mu0; /* viscosity of pure fluid */
  } else {
    maxpack = .68;
    mu0 = mu; /* viscosity of pure fluid */
  }

  y_norm = Y[w] / maxpack;
  maxpack2 = maxpack * maxpack;
  if (y_norm < 0.95) {
    comp = pow((1. - y_norm), -2.);
    comp1 = 2. / maxpack * pow((1. - y_norm), -3.);
    comp2 = 6.0 / maxpack2 * pow((1. - y_norm), -4.);
  } else {
    comp = 0.;
    comp1 = 0.;
  }

  /* Migrate this input deck later */
  Kn = 0.75;

  /* This is the particle pressure */
  pp = Kn * y_norm * y_norm * comp;

  d_pp_dy = 2. * Kn * y_norm / maxpack * comp + Kn * y_norm * y_norm * comp1;

  d_pp2_dy2 =
      2. * Kn / maxpack2 * comp + 4. * Kn * y_norm / maxpack * comp1 + Kn * y_norm * y_norm * comp2;

  if (mp->SBM_Length_enabled) {
    radius_p = mp->SBM_Lengths2[w][0];
    L_char = mp->SBM_Lengths2[w][1];
    U_max = mp->SBM_Lengths2[w][2];
  } else {
    radius_p = 0.;
    L_char = 1.;
    U_max = 1.;
  }

  gamma_nl = radius_p * U_max / (L_char * L_char);

  memset(div_tau_p, 0, DIM * sizeof(dbl));

  for (a = 0; a < WIM; a++) {
    for (b = 0; b < WIM; b++) {
      div_tau_p[a] +=
          mu0 * qtensor[a][b] * (pp * grad_gd[b] + (gammadot + gamma_nl) * d_pp_dy * grad_Y[w][b]);
    }
  }

  if (af->Assemble_Jacobian) {
    memset(d_div_tau_p_dgd, 0, DIM * MDE * sizeof(dbl));
    memset(d_div_tau_p_dv, 0, DIM * DIM * MDE * sizeof(dbl));
    memset(d_div_tau_p_dy, 0, DIM * MAX_CONC * MDE * sizeof(dbl));
    memset(d_div_tau_p_dmesh, 0, DIM * DIM * MDE * sizeof(dbl));
    memset(d_div_tau_p_dp, 0, DIM * MDE * sizeof(dbl));

    var = MASS_FRACTION;
    if (pd->v[pg->imtrx][var]) {
      dofs = ei[pg->imtrx]->dof[var];
      for (a = 0; a < WIM; a++) {
        for (j = 0; j < dofs; j++) {
          for (b = 0; b < WIM; b++) {
            d_div_tau_p_dy[a][w][j] +=
                mu0 * qtensor[a][b] *
                (d_pp_dy * bf[var]->phi[j] * grad_gd[b] +
                 (gammadot + gamma_nl) * d_pp2_dy2 * bf[var]->phi[j] * grad_Y[w][b] +
                 (gammadot + gamma_nl) * d_pp_dy * bf[var]->grad_phi[j][b]);
          }
        }
      }
    }

    var = SHEAR_RATE;
    if (pd->v[pg->imtrx][var]) {
      dofs = ei[pg->imtrx]->dof[var];
      for (a = 0; a < WIM; a++) {
        for (j = 0; j < dofs; j++) {
          for (b = 0; b < WIM; b++) {
            d_div_tau_p_dgd[a][j] +=
                mu0 * qtensor[a][b] *
                (pp * bf[var]->grad_phi[j][b] + bf[var]->phi[j] * d_pp_dy * grad_Y[w][b]);
          }
        }
      }
    }

    var = MESH_DISPLACEMENT1;
    if (pd->v[pg->imtrx][var]) {
      dofs = ei[pg->imtrx]->dof[var];
      for (p = 0; p < dim; p++) {
        for (j = 0; j < dofs; j++) {
          for (a = 0; a < WIM; a++) {
            for (b = 0; b < WIM; b++) {
              d_div_tau_p_dmesh[a][p][j] +=
                  mu0 * qtensor[a][b] *
                  (pp * fv->d_grad_SH_dmesh[b][p][j] +
                   (gammadot + gamma_nl) * d_pp_dy * fv->d_grad_c_dmesh[b][w][p][j]);
            }
          }
        }
      }
    }

  } /* end if Jacobian */

  return (status);
}

void rotate_tensor(double A[DIM][DIM], double A_prime[DIM][DIM], double R0[DIM][DIM], int dir) {

  /* Rotates a tensor from A to A_prime using the orthogonal tensor R */
  /* dir = 0 : A_prime = R * A * Rt */
  /* dir = 1 : A_prime = Rt * A * R */

  double R[DIM][DIM];
  double Rt[DIM][DIM];
  double R_tmp[DIM][DIM];
  double R_dot_A[DIM][DIM];
  int i, j, k;

  memset(R, 0, sizeof(dbl) * DIM * DIM);
  memset(R_tmp, 0, sizeof(dbl) * DIM * DIM);
  memset(A_prime, 0, sizeof(dbl) * DIM * DIM);
  memset(R_dot_A, 0, sizeof(dbl) * DIM * DIM);

  for (i = 0; i < DIM; i++) {
    for (j = 0; j < DIM; j++) {
      R[i][j] = R0[i][j];
    }
  }

  if (dir == 1) {
    for (i = 0; i < DIM; i++) {
      for (j = 0; j < DIM; j++) {
        R_tmp[i][j] = R[i][j];
        Rt[i][j] = R[i][j];
      }
    }
    for (i = 0; i < DIM; i++) {
      for (j = 0; j < DIM; j++) {
        R[i][j] = R_tmp[j][i];
      }
    }
  } else {
    for (i = 0; i < DIM; i++) {
      for (j = 0; j < DIM; j++) {
        Rt[i][j] = R[j][i];
      }
    }
  }

  for (i = 0; i < DIM; i++) {
    for (j = 0; j < DIM; j++) {
      R_dot_A[i][j] = 0;
      for (k = 0; k < DIM; k++) {
        R_dot_A[i][j] += R[i][k] * A[k][j];
      }
    }
  }

  for (i = 0; i < DIM; i++) {
    for (j = 0; j < DIM; j++) {
      A_prime[i][j] = 0;
      for (k = 0; k < DIM; k++) {
        A_prime[i][j] += R_dot_A[i][k] * Rt[k][j];
      }
    }
  }
}

/******************************************************************************
 *
 *  Routine which calculates the Vapor Pressure of saturated liquid
 *  based on Antoine's equation
 *     Author: A.C. Sun 6/98
 *
 ******************************************************************************/

int antoine_psat(double param[], double *f, double *dfdt)

{
  /*local variables */

  dbl A, B, C;    /*Antoine coefficients*/
  dbl CONV;       /*Unit conversion factor*/
  dbl TMIN, TMAX; /*Temperature range */
  dbl T;

  T = fv->T;

  CONV = param[0];
  A = param[1];
  B = param[2];
  C = param[3];
  TMIN = param[4];
  TMAX = param[5];

  /* calculate the vapor pressure based on Antoine relation */

  if (T <= TMIN || T >= TMAX) {
    *f = CONV * 0.1;
    *dfdt = 0.0;
  } else {
    *f = CONV * exp(A - B / (T + C));

    /* calculate the temperature derivative */

    *dfdt = CONV * exp(A - B / (T + C)) * (B / ((T + C) * (T + C)));
  }

  return 0;
}
/*****************************************************************************/
/* END of routine antoine_psat */
/*****************************************************************************/

/******************************************************************************
 *
 *  Routine which calculates the Vapor Pressure of saturated liquid
 *  based on Riedel's Equation
 *     Author: A.C. Sun 6/98
 *
 ******************************************************************************/

int riedel_psat(double param[], double *f, double *dfdt) {
  /*local variables */

  dbl A, B, C, D, E; /*Riedel coefficients*/
  dbl CONV;          /*Unit conversion factor*/
  dbl func, dfunc;
  dbl T, TMIN, TMAX;

  T = fv->T;

  CONV = param[0];
  A = param[1];
  B = param[2];
  C = param[3];
  D = param[4];
  E = param[5];
  TMIN = param[6];
  TMAX = param[7];

  if (T <= TMIN || T >= TMAX) {
    *f = CONV * 0.1;
    *dfdt = 0.0;
  } else {
    /* calculate the vapor pressure based on Riedel's
       relation */

    func = exp(A + (B / T) + C * log(T) + D * pow(T, E));
    *f = CONV * func;

    /* calculate the temperature derivative */

    dfunc = func * ((-B / (T * T)) + (C / T) + D * E * pow(T, E - 1));
    *dfdt = CONV * dfunc;
  }
  return 0;
}
/*****************************************************************************/
/* END of routine riedel_psat */
/*****************************************************************************/

/* MMH
 * Calculate the source term and functional dependencies
 * for the SUSPENSION_PM model that apply to the FLUID
 * momentum equation(s).
 */
int suspension_pm_fluid_momentum_source(dbl f[DIM], /* Body force. */
                                        MOMENTUM_SOURCE_DEPENDENCE_STRUCT *df) {
  int status = 1;
  int dim, a, j;
  int eqn, species;
  dbl rho_f, rho_p, delta_rho;
  dbl p_vol_frac, pvf2;

  if (mp->DensityModel == SUSPENSION_PM) {
    species = (int)mp->u_density[0];
    rho_f = mp->u_density[1];
    rho_p = mp->u_density[2];
    delta_rho = rho_p - rho_f;
    pvf2 = p_vol_frac = fv->c[species];
    pvf2 *= pvf2;
  } else {
    GOMA_EH(GOMA_ERROR,
            "MomentumSourceModel is SUSPENSION_PM, but DensityModel is not.\nThis is bad.");
    return (status);
  }

  dim = pd->Num_Dim;

  /* MMH
   * The fluid momentum buoyancy source term due to particles.
   * And the sensitivity of the source term to the species concentration.
   */
  for (a = 0; a < dim; a++) {
    eqn = R_MOMENTUM1 + a;
    if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
      f[a] = mp->momentum_source[a] * (pvf2 * rho_p + (1.0 - pvf2) * rho_f);
      /* For sensitivity of source to species.  The only
       * species that should be present here is the particle
       * phase.  */
      for (j = 0; j < ei[pg->imtrx]->dof[MASS_FRACTION]; j++)
        /* MMH here, I switched from += to = (???) */
        df->C[a][species][j] =
            2.0 * p_vol_frac * delta_rho * mp->momentum_source[a] * bf[MASS_FRACTION]->phi[j];
    }
  }
  return (status);
}

/* MMH
 * Calculate the source term and functional dependencies
 * for the SUSPENSION_PM model that apply to the PARTICLE
 * momentum equation(s).
 */
int suspension_pm_particle_momentum_source(
    dbl f[DIM],                   /* Body force. */
    dbl dfdT[DIM][MDE],           /* For temperature dependence */
    dbl dfdX[DIM][DIM][MDE],      /* For spatial dependence */
    dbl dfdC[DIM][MAX_CONC][MDE], /* For concentration dependence */
    dbl dfdv[DIM][DIM][MDE])      /* For velocity dependence */
{
  int status = 1;
  const int dim = pd->Num_Dim;
  int a, j;
  int eqn, species;
  dbl rho_f, rho_p, delta_rho;
  dbl p_vol_frac;

  if (mp->DensityModel == SUSPENSION_PM) {
    species = (int)mp->u_density[0];
    rho_f = mp->u_density[1];
    rho_p = mp->u_density[2];
    delta_rho = rho_p - rho_f;
    p_vol_frac = fv->c[species];
  } else {
    GOMA_EH(GOMA_ERROR,
            "MomentumSourceModel is SUSPENSION_PM, but DensityModel is not.\nThis is bad.");
    return (status);
  }

  /* MMH
   * The particle momentum buoyancy source term due to particles.
   * And the sensitivity of the source term to the species concentration.
   */
  for (a = 0; a < dim; a++) {
    eqn = R_PMOMENTUM1 + a;
    if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
      f[a] = mp->momentum_source[a] * (p_vol_frac * (1.0 - p_vol_frac) * delta_rho);
      /* For sensitivity of source to species.  The only
       * species that should be present here is the particle
       * phase.
       */
      for (j = 0; j < ei[pg->imtrx]->dof[MASS_FRACTION]; j++) {
        /* MMH: I changed this from += to = (??? what was += there for before?!) */
        dfdC[a][species][j] = (1.0 - 2.0 * p_vol_frac) * delta_rho * mp->momentum_source[a] *
                              bf[MASS_FRACTION]->phi[j];
      }
    }
  }
  return (status);
}

/*
 * Molten Glass Viscosity Model
 */

/*
 * int molten_glass_viscosity ( vis,dvis_dT)
 *
 * ----------------------------------------------------------------------------
 * This routine is responsible for filling up the following properties and sensitivities
 * at the current gauss point:
 *     intput:
 *
 *     output:  vis         - Flowing liquid viscosity
 *              dvis_dT[j]  - derivative wrt temperature at node j.
 *   NB: The user need only supply f, dfdT, dfdC, etc....mp struct is loaded up for you
 * ----------------------------------------------------------------------------
 */

int molten_glass_viscosity(dbl *vis,         /* Base FLOWING LIQUID VISCOITY  */
                           dbl dvis_dT[MDE], /* temperature dependence. */
                           dbl *param)       /* parameter list */
{
  int eqn, var;
  dbl AA, BB, CC;
  dbl T; /* Convenient local variables */
  int j;

  /* Begin Execution */

  AA = mp->u_FlowingLiquid_viscosity[0];
  BB = mp->u_FlowingLiquid_viscosity[1];
  CC = mp->u_FlowingLiquid_viscosity[2];

  /************Initialize everything for safety**************/

  eqn = R_MOMENTUM1;
  if (pd->e[pg->imtrx][eqn] & T_POROUS_BRINK) {
    *vis = 0.;

    if (dvis_dT != NULL) {
      var = TEMPERATURE;
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        dvis_dT[j] = 0.;
      }
    }
  }
  /**********************************************************/

  /***********Load up convenient local variables*************/
  /*NB This ought to be done once for all fields at gauss pt*/

  T = fv->T;

  /**********************************************************/
  *vis = pow(10.0, AA + BB / (T - CC));

  /* Now do sensitivities */

  if (dvis_dT != NULL) {
    var = TEMPERATURE;
    for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
      dvis_dT[j] = -BB * exp(AA + BB / (T - CC)) / pow(T - CC, 2.0);
    }
  }

  return (0);
}
/*****************************************************************************/
/* END of routine molten_glass_viscosity */
/*****************************************************************************/

/*
 * Molten Glass Viscosity Model
 */

/*
 * epoxy_flowing_liquid_viscosity ( vis, d_flow_mu)
 *
 * ----------------------------------------------------------------------------
 * This routine is responsible for filling up the following properties and sensitivities
 * at the current gauss point:
 *     intput:
 *
 *     output:  vis         - Flowing liquid viscosity
 *              d_flow_mu   - derivative wrt DOF at node j.
 * ----------------------------------------------------------------------------
 */

int epoxy_flowing_liquid_viscosity(dbl *vis, /* Base FLOWING LIQUID VISCOITY  */
                                   VISCOSITY_DEPENDENCE_STRUCT *d_flow_mu, /* its dependence. */
                                   dbl *param)                             /* parameter list */
{
  /* Local Variables */

  dbl mu0 = param[0];          /* monomer reference temperature viscosity */
  dbl alpha_g = param[1];      /* extent of reaction at the gel point */
  dbl A = param[2];            /* exponent for constitutive equation */
  dbl B = param[3];            /* exponent for constitutive equation */
  dbl Aexp = param[4];         /* exponent for thermal viscosity dependence */
  int species = (int)param[5]; /* species number for cure equation */

  dbl mu;    /* viscosity */
  dbl alpha; /* extent of reaction */
  dbl exponent;
  dbl ratio;
  dbl deriv; /* stuff for the first derivative */
  dbl T;     /* Convenient local variables */
  int j, var;

  /* Begin Execution */

  /******** Evaluate terms for viscosity *********/

  alpha = fv->c[species]; /* extent of reaction */

  if (alpha < alpha_g) {
    ratio = (alpha_g) / (alpha_g - alpha);
    exponent = A + B * alpha;
  } else /* do something special at the gel point */
  {
    ratio = 100000;
    exponent = A + B * alpha_g;
  }

  if (pd->gv[TEMPERATURE]) {
    T = fv->T;
  } else {
    T = upd->Process_Temperature;
  }

  if (T <= 0.) {
    mu = mu0 * pow(ratio, exponent);
  } else {
    mu = mu0 * exp(Aexp / T) * pow(ratio, exponent);
  }

  *vis = mu;

  /******** Evaluate viscosity senstivities *********/

  if (d_flow_mu != NULL) {

    /* d_flow_mu_dT */
    var = TEMPERATURE;
    if (pd->v[pg->imtrx][var]) {
      if (T <= 0.) {
        mp->d_FlowingLiquid_viscosity[var] = 0.;
      } else {
        mp->d_FlowingLiquid_viscosity[var] = -mu * Aexp / (T * T);
      }

      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_flow_mu->T[j] = mp->d_FlowingLiquid_viscosity[var] * bf[var]->phi[j];
      }
    }

    /* d_flow_mu_dC */
    var = MASS_FRACTION;
    if (pd->v[pg->imtrx][var]) {
      if (alpha < alpha_g) {
        deriv = exponent / (alpha_g - alpha) + B * log(ratio);
        mp->d_FlowingLiquid_viscosity[MAX_VARIABLE_TYPES + species] = mu * deriv;
      } else {
        deriv = 0.0;
        mp->d_FlowingLiquid_viscosity[MAX_VARIABLE_TYPES + species] = 0.;
      }

      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_flow_mu->C[species][j] = mp->d_FlowingLiquid_viscosity[var] * bf[var]->phi[j];
      }
    }
  }
  return (0);
}
/*****************************************************************************/
/* END of routine epoxy_flowing_liquid_viscosity */
/*****************************************************************************/

/*
 * Solidification Permeability Model
 */

/*
 * int solidification_permeability (h_elem_avg )
 *
 * ----------------------------------------------------------------------------
 * This routine is responsible for filling up the following properties and sensitivities
 * at the current gauss point:
 *     input:   h_elem_avg - average element size
 *              mp->u_permeability - permeability info from input deck
 *
 *     return:  per       - permeability that depends on solid fraction
 *     output:  d_per_dc  - derivative wrt solid fraction
 * ----------------------------------------------------------------------------
 */

dbl solidification_permeability(dbl h_elem_avg, /* average element size */
                                dbl d_per_dc[MAX_CONC][MDE]) {
  int eqn, var;
  int species;
  int j;
  dbl siz, vol; /* Convenient local variables */
  dbl per;
  dbl mu0, maxpack;

  /* Begin Execution */

  /* this is the species that indicates solidification */
  species = (int)mp->u_permeability[0];
  /* viscosity of pure fluid */
  mu0 = gn->mu0;

  if (gn->ConstitutiveEquation == SUSPENSION || gn->ConstitutiveEquation == CARREAU_SUSPENSION ||
      gn->ConstitutiveEquation == POWERLAW_SUSPENSION || gn->ConstitutiveEquation == FILLED_EPOXY) {
    maxpack = gn->maxpack;
  } else {
    maxpack = .68;
  }

  /* good default behavior */
  per = 1.;

  maxpack = 1.;

  vol = fv->c[species] / maxpack;

  siz = 0.5 * h_elem_avg * h_elem_avg;
  if (vol > 0.999)
    vol = 0.999;

  eqn = R_MOMENTUM1;
  if (pd->e[pg->imtrx][eqn] & T_POROUS_BRINK) {
    if (vol > 0.) {
      per = siz * (1. - vol) * (1. - vol) / (vol * (1.43 - vol)) / mu0;
    } else {
      per = siz * 1.e12;
    }
  }

  memset(d_per_dc, 0, MAX_CONC * MDE * sizeof(dbl));

  /* Now do sensitivies */
  var = MASS_FRACTION;
  if (pd->v[pg->imtrx][var]) {
    if (vol > 0.) {
      mp->d_permeability[MAX_VARIABLE_TYPES + species] =
          siz * (1. - vol) / (vol * (1.43 - vol)) *
          (-2. - (1. - vol) / vol + (1. - vol) / (1.43 - vol)) / mu0 / maxpack;
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        d_per_dc[species][j] = mp->d_permeability[MAX_VARIABLE_TYPES + species] * bf[var]->phi[j];
      }
    } else {
      mp->d_permeability[MAX_VARIABLE_TYPES + species] = 0.;
    }
  }

  return (per);
}
/*****************************************************************************/
/* END of solidification_permeability */
/*****************************************************************************/
/*
 * This is the reaction source term for the species
 * equation using an extent-of-reaction model
 * There are 2 reactions for REF foam.  This is based on the unit of
 * mass-surface flux, m/l^2/t
 * Due to the fact that GOMA usually solves species terms by specifying
 * number of independent variables.  In this routine, all of the source
 * terms are written and are explicit for Num_Species.
 */

int foam_species_source(double *param) {
  int eqn, imtrx;
  int j;
  dbl foam, gas, s1; /*mass fractions of gas and solid_1 */
  dbl T;             /* temperature for rate constants */
  dbl A1, E1, A2, E2, expon1, expon2, sigma1, sigma2;
  dbl k1, k2, r1, r2;
  dbl dr1_dT, dr2_dT, delT, ref_T1, ref_T2, T_actual;
  // dbl char_dist; /* characteristic distance, 1mm burnt depth */
  /*  dbl xx[] = {3.4, 3.2,3.,2.8,2.6,2.4,2.2,2,1.8,1.6,1.4,1.2,
      1.,0.8,0.6,0.4,0.2,0. };
      dbl yy[] = {.9997,.9993,.9987,.9974,.9953,.9918,.9861,.9772,
      .9641,.9452,.9192,.8849,.8413,.7881,.7257,.6554,.5793,.5}; */
  dbl rpar[] = {0., 0.}, ext[] = {0., 0.}, drpar[] = {0., 0.};
  dbl dr1_df, dr1_dg, dr1_ds1, dr2_df, dr2_dg, dr2_ds1;

  if (MAX_CONC < 3) {
    GOMA_EH(GOMA_ERROR, "foam_species_source expects MAX_CONC >= 3");
    return -1;
  }

  /* Begin Execution */

  T = fv->T;

  foam = fv->c[0];
  gas = fv->c[1];
  s1 = fv->c[2];

  for (j = 0; j < pd->Num_Species; j++) {
    if (fv->c[j] <= 1.e-10)
      fv->c[j] = 1.e-10;
  }

  ext[0] = (1. - foam);
  ext[1] = gas;

  for (j = 0; j < 2; j++) {
    rpar[j] = exp((MAX((1. - ext[j]), ext[j]) - .8415) / .1767);
    drpar[j] = (1. / .1767) * rpar[j];
    if (ext[j] < 0.5)
      drpar[j] = -drpar[j];
    if (ext[j] == 0.)
      drpar[j] = 0.;
  }
  drpar[0] = -drpar[0];

  /* s2 is the mass fraction of another fragment */
  A1 = param[0];
  expon1 = param[1];
  sigma1 = param[2];
  A2 = param[3];
  expon2 = param[4];
  sigma2 = param[5];
  ref_T1 = param[6];
  ref_T2 = param[7];

  delT = ref_T2 - ref_T1;
  T_actual = delT * T + ref_T1;

  E1 = (expon1 + rpar[0] * sigma1) / 1.987;
  E2 = (expon2 + rpar[1] * sigma2) / 1.987;

  /* Arhenius type rate constants for extent of reaction model */
  if (T <= 0.) {
    k1 = A1;
    k2 = A2;
  } else {
    k1 = exp(log(A1) - E1 / T_actual);
    k2 = exp(log(A2) - E2 / T_actual);
  }
  /* the evolution of the first reaction is only dependent on the volume
     shrinkage of the 1-D foam, thus far anyway */

  r1 = k1 * foam;
  r2 = k2 * s1;

  dr1_df = k1 * (1. - drpar[0] * sigma1 * foam / (1.987 * T_actual));
  dr1_dg = 0.;
  dr1_ds1 = 0.;

  dr2_df = 0.;
  dr2_dg = 0.;
  dr2_ds1 = k2 * (1. - drpar[1] * sigma2 / (1.987 * T_actual));

  dr1_dT = r1 * (E1 * delT / pow(T_actual, 2.0));
  dr2_dT = r2 * (E2 * delT / pow(T_actual, 2.0));

  /**********************************************************/

  /* Species piece , this is num_species+1 Jacobian */
  eqn = MASS_FRACTION;
  for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
    if (pd->e[imtrx][eqn] & T_SOURCE) {
      mp->species_source[0] = -r1;
      mp->species_source[1] = (0.3 * r1 + 0.943 * r2);
      mp->species_source[2] = (0.7 * r1 - r2);

      /* Jacobian entries for source term */

      mp->Jac_Species_Source[0] = -dr1_df;
      mp->Jac_Species_Source[1] = -dr1_dg;
      mp->Jac_Species_Source[2] = -dr1_ds1;

      mp->Jac_Species_Source[3] = 0.3 * dr1_df + 0.943 * dr2_df;
      mp->Jac_Species_Source[4] = 0.3 * dr1_dg + 0.943 * dr2_dg;
      mp->Jac_Species_Source[5] = 0.3 * dr1_ds1 + 0.943 * dr2_ds1;

      mp->Jac_Species_Source[6] = 0.7 * dr1_df - dr2_df;
      mp->Jac_Species_Source[7] = 0.7 * dr1_dg - dr2_dg;
      mp->Jac_Species_Source[8] = 0.7 * dr1_ds1 - dr2_ds1;

      mp->d_species_source[MAX_VARIABLE_TYPES + 0] = -dr1_dT;
      mp->d_species_source[MAX_VARIABLE_TYPES + 1] = (0.3 * dr1_dT + 0.943 * dr2_dT);
      mp->d_species_source[MAX_VARIABLE_TYPES + 2] = (0.7 * dr1_dT - dr2_dT);
    }
  }
  return 0;
}

/*
 * This routine computes the molar rate of electrolyte-species consumption
 * due to electrochemical reactions in porous electrode regions of a
 * thermal battery cell. The molar rate of electrolyte-species consumption
 * is evaluated using Bulter-Volmer kinetics along with Faraday's law.
 * (see SAND 2000-0207 of Chen et al. 2000 for reference).
 *
 * KSC (1/2002): This latest routine was based on an earlier version (8/99);
 *               it also incorporates sensitivity addition made by RSL (7/00)
 *               with slight modification.
 */
int electrode_species_source(int species_no, /* Current species number */
                             double time,    /* present time value */
                             double delta_t) /* present time step  */
{
  int status = 0;
  int j;

  const double R = 8.314;   /* Universal gas constant in units of Joules/mole K */
  const double F = 96487.0; /* Faraday's constant in units of Coulomb/mole */
  dbl ai0 = 0.0;            /* Product of surface area by exchange current density */
  dbl U0 = 0.0;             /* Thermodynamic open-circuity potential */
  dbl alphaa = 0.0;         /* Anodic direction transfer coefficient  */
  dbl alphac = 0.0;         /* Cathodic direction transfer coefficient */
  dbl util;                 /* Fractional utilization of an electrode */
  dbl tau;                  /* Parameter for calculating electrode fractional utilization */
  dbl E;                    /* Activation energy for exchange current density */
  dbl x[MAX_CONC];          /* Mole fraction of electrolyte species */
  dbl T;                    /* Temperature of electrolyte solution */
  dbl PHI1;                 /* Electrical potential in solid electrode phase */
  dbl PHI2;                 /* Electrical potential in liquid electrolyte phase */
  dbl La, Lc;               /* Electrode thicknesses */
  dbl ea, ec;               /* Electrode porosities */
  dbl Va, Vc;               /* Molar volume of active electrode materials */
  dbl i;                    /* Current density */
  dbl na, nc;               /* number of electrons involved in the reduction/oxidation rxns */
  int mn;                   /* Region index: mn=0 indicates anode and mn=2 cathode */
  dbl FRT;                  /* F/R/T */
  dbl eta;
  dbl sum;
  dbl x0p, T0, ai00, util0, ER;
  dbl util1, util2, util3;

  if (MAX_CONC < 3) {
    GOMA_EH(GOMA_ERROR, "electrode_species_source expects MAX_CONC >= 3");
    return -1;
  }

  /* Begin Execution */

  mn = ei[pg->imtrx]
           ->mn; /* get region index: mn=1 for anode, mn=2 for cathode, mn=1 for separator */
  PHI1 = fv->T;  /* Electrode potential solved using energy transport equation */
  PHI2 = fv->V;  /* Electrolyte potential solved using charge conservation eq. */
  electrolyte_temperature(time, delta_t, 0); /* get electrolyte temperature at present time */
  T = mp->electrolyte_temperature;
  FRT = F / R / T;

  for (j = 0; j < pd->Num_Species_Eqn; j++) {
    x[j] = fv->c[j];
  }

  if (mn == 0) {
    alphaa = mp->u_reaction_rate[0];
    alphac = mp->u_reaction_rate[1];

    if (mp->InterfacialAreaModel == CONSTANT) {
      ai0 = mp->interfacial_area;
    } else if (mp->InterfacialAreaModel == THERMAL_BATTERY) {
      GOMA_EH(GOMA_ERROR, "Non-constant interfacial-area model for anode awaits to be implemented");
    }

    if (mp->ThermodynamicPotentialModel ==
        LiSi) /* Calculate temperature-dependent thermodynamic potential for */
    {         /* Li(Si) anode; see SAND2000-0207 by Chen et al. 2000, p.18-19 */
      util1 = mp->u_thermodynamic_potential[0];
      util2 = mp->u_thermodynamic_potential[1];
      La = mp->u_thermodynamic_potential[2]; /* anode thickness */
      ea = mp->u_thermodynamic_potential[3]; /* anode porosity */
      Va = mp->u_thermodynamic_potential[4]; /* molar volume of active anode material */
      i = mp->u_thermodynamic_potential[5];  /* current density */
      na = mp->u_thermodynamic_potential[6]; /* number of electrons involved in anode rxn */
      tau = na * La * ea * F / Va / i;
      util = time / tau;

      if (util <= util1) {
        U0 = -0.187529 + 0.0000731 * T;
      } else if (util > util1 && util <= util2) {
        U0 = -0.088097 + 0.0001122 * T;
      } else if (util > util2) {
        U0 = -0.0345 + 0.0001056 * T;
      }
    } else if (mp->ThermodynamicPotentialModel == CONSTANT) /* Constant thermodynamic potential */
    {
      U0 = mp->thermodynamic_potential;
    }

    eta = PHI1 - PHI2 - U0; /* overpotential */
    mp->species_source[0] = (ai0 / F) * x[0] * (exp(alphaa * FRT * eta) - exp(-alphac * FRT * eta));
    mp->species_source[1] = 0.0; /* species 1 is being neither produced nor consumed */
    mp->species_source[2] = 0.0; /* species 2 is being neither produced nor consumed */
  } else if (mn == 2) {
    alphaa = mp->u_reaction_rate[0];
    alphac = mp->u_reaction_rate[1];

    if (mp->InterfacialAreaModel == CONSTANT) {
      ai0 = mp->interfacial_area;
    } else if (mp->InterfacialAreaModel == THERMAL_BATTERY) {
      ai00 = mp->u_interfacial_area[0];
      util0 = mp->u_interfacial_area[1];
      E = mp->u_interfacial_area[2] * 4.184; /* convert E from cal/mol-K to J/mol-k */
      T0 = mp->u_interfacial_area[3];
      Lc = mp->u_interfacial_area[4]; /* cathode thickness */
      ec = mp->u_interfacial_area[5]; /* cathode porosity */
      Vc = mp->u_interfacial_area[6]; /* molar volume of active cathode material */
      i = mp->u_interfacial_area[7];  /* current density */
      nc = mp->u_interfacial_area[8]; /* number of electrons invloved in cathode rxn */
      tau = nc * Lc * ec * F / Vc / i;
      util = time / tau;
      ER = E / R;
      ai0 = ai00 * (1.0 - util / util0) * exp(-ER * (1.0 / T - 1.0 / T0));
    } else {
      GOMA_EH(
          GOMA_ERROR,
          "Variable interfacial-area model other than THERMAL_BATTERY awaits to be implemented");
    }

    if (mp->ThermodynamicPotentialModel ==
        FeS2) /* Calculate temperature-dependent thermodynamic potential for */
    {         /* FeS2 cathode; see SAND2000-0207 by Chen et al. 2000, p.19   */
      util1 = mp->u_thermodynamic_potential[0];
      util2 = mp->u_thermodynamic_potential[1];
      util3 = mp->u_thermodynamic_potential[2];
      Lc = mp->u_thermodynamic_potential[3]; /* cathode thickness */
      ec = mp->u_thermodynamic_potential[4]; /* cathode porosity */
      Vc = mp->u_thermodynamic_potential[5]; /* molar volume of active cathode material */
      i = mp->u_thermodynamic_potential[6];  /* current density */
      nc = mp->u_thermodynamic_potential[7]; /* number of electrons involved in cathode rxn */
      tau = nc * Lc * ec * F / Vc / i;
      util = time / tau;

      if (util < util1) {
        U0 = 1.4251 + 0.0004785 * T;
      } else if (util >= util1 && util <= util2) {
        U0 = 1.208771 + 0.00065142 * T;
      } else if (util >= util2 && util <= util3) {
        x0p = 0.91658 - 9.24e-05 * (T - 273.);
        U0 = 1.208771 + 0.00065142 * T +
             (0.130129 - 0.00063812 * T) /
                 (1. - (2.2 * (4.0 * x0p - 2.0) / (2.0 * x0p - 0.8) - 3.0)) * (util - 0.434) /
                 (0.5 - 0.434) * 0.53;
      } else if (util > util3) {
        U0 = 1.43211 - 0.000147 * T;
      }
    } else if (mp->ThermodynamicPotentialModel == CONSTANT) /* Constant thermodynamic potential */
    {
      U0 = mp->thermodynamic_potential;
    }

    eta = PHI1 - PHI2 - U0; /* overpotential */
    mp->species_source[0] = (ai0 / F) * x[0] * (exp(alphaa * FRT * eta) - exp(-alphac * FRT * eta));
    mp->species_source[1] = 0.0; /* species 1 is being neither produced nor consumed */
    mp->species_source[2] = 0.0; /* species 2 is being neither produced nor consumed */
  } else /* no electrochemical reactions take place in the separator region */
  {
    mp->species_source[0] = 0.0;
    mp->species_source[1] = 0.0;
    mp->species_source[2] = 0.0;
  }

  /* Sensitivities */
  mp->d_species_source[MAX_VARIABLE_TYPES + 1] = 0.; /* wrt species # 2 */
  if (species_no == 1 || species_no == 2 || mn == 1) /* for separator region and species 2 & 3 */
  {
    mp->d_species_source[MAX_VARIABLE_TYPES + 0] = 0.; /* wrt species # 1 */
    mp->d_species_source[TEMPERATURE] = 0.;            /* wrt electrode potential */
    mp->d_species_source[VOLTAGE] = 0.;                /* wrt electrolyte potential */
  } else /* for anode and cathode regions and when species_no = 0 */
  {
    eta = PHI1 - PHI2 - U0;
    sum = alphaa * exp(alphaa * FRT * eta) + alphac * exp(-alphac * FRT * eta);
    mp->d_species_source[MAX_VARIABLE_TYPES + 0] = mp->species_source[0] / x[0]; /* wrt species 1 */
    mp->d_species_source[TEMPERATURE] = (ai0 / F) * x[0] * FRT * sum; /* wrt electrode potential */
    mp->d_species_source[VOLTAGE] = -mp->d_species_source[TEMPERATURE]; /* wrt electrolyte poten. */
  }

  return (status);
}
/*****************************************************************************/
/* END of routine electrode_species_source                                   */
/*****************************************************************************/

int ion_reaction_source(int species_no) /* current species number */

/*
 * This is the species source term that accounts for the homogeneous chemical reactions
 * involved in nickel electroplating.  Obviously, the reactions are not considered to
 * be at equilibrium here, although the rate constants can be set large enough to make
 * this effectively true.
 *
 * RSL 3/21/01
 *
 */

{

#if MAX_CONC < 7
  GOMA_EH(GOMA_ERROR, "ion_reaction_source expects MAX_CONC to be >= 7");
  return -1;
#else

  int eqn, var, j;
  int four, five;
  dbl k1, k2, k3, K1, K2, K3;
  dbl Q1 = 0.0, Q2 = 0.0, Q3 = 0.0;
  dbl dQ1dx2 = 0.0, dQ1dx3 = 0.0, dQ2dx5 = 0.0, dQ2dx1 = 0.0, dQ2dx2 = 0.0, dQ3dx4 = 0.0,
      dQ3dx0 = 0.0, dQ3dx3 = 0.0;
  dbl c, rho, M_mix, x[MAX_CONC] = {0};

  /* Begin Execution */

  four = 4;
  five = 5;

  k1 = 1.3e-06;
  k2 = 1.3e+09;
  k3 = 5.9e+06;
  K1 = 1.01e-20;
  K2 = 1.01e-05;
  K3 = 4.5e-08;

  M_mix = 0.;
  for (j = 0; j < pd->Num_Species; j++) {
    x[j] = fv->c[j];
    M_mix += x[j] * mp->molecular_weight[j];
  }
  rho = density(NULL, 0.0); /*  RSL 6/22/02  */
  c = rho / M_mix;

  if (species_no == 2 || species_no == 3) {
    Q1 = k1 * (1. - c * c * x[2] * x[3] / K1);
    dQ1dx2 = -k1 * c * c * x[3] / K1;
    dQ1dx3 = -k1 * c * c * x[2] / K1;
  }

  if (species_no == 1 || species_no == 2 || species_no == 5) {
    Q2 = k2 * c * (x[five] - c * x[1] * x[2] / K2);
    dQ2dx5 = k2 * c;
    dQ2dx1 = -k2 * c * c * x[2] / K2;
    dQ2dx2 = -k2 * c * c * x[1] / K2;
  }

  if (species_no == 0 || species_no == 3 || species_no == 4) {
    Q3 = k3 * c * (x[four] - c * x[0] * x[3] / K3);
    dQ3dx4 = k3 * c;
    dQ3dx0 = -k3 * c * c * x[3] / K3;
    dQ3dx3 = -k3 * c * c * x[0] / K3;
  }

  eqn = MASS_FRACTION;
  for (int imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
    if (pd->e[imtrx][eqn] & T_SOURCE) {
      switch (species_no) {
      case 0:
        mp->species_source[species_no] = Q3;

        var = TEMPERATURE;
        if (pd->v[pg->imtrx][var]) {
          mp->d_species_source[var] = 0.;
        }

        var = MASS_FRACTION;
        if (pd->v[pg->imtrx][var]) {
          mp->d_species_source[MAX_VARIABLE_TYPES + 0] = dQ3dx0;
          mp->d_species_source[MAX_VARIABLE_TYPES + 1] = 0.;
          mp->d_species_source[MAX_VARIABLE_TYPES + 2] = 0.;
          mp->d_species_source[MAX_VARIABLE_TYPES + 3] = dQ3dx3;
          mp->d_species_source[MAX_VARIABLE_TYPES + four] = dQ3dx4;
          mp->d_species_source[MAX_VARIABLE_TYPES + five] = 0.;
        }
        break;
      case 1:
        mp->species_source[species_no] = Q2;

        var = TEMPERATURE;
        if (pd->v[pg->imtrx][var]) {
          mp->d_species_source[var] = 0.;
        }

        var = MASS_FRACTION;
        if (pd->v[pg->imtrx][var]) {
          mp->d_species_source[MAX_VARIABLE_TYPES + 0] = 0.;
          mp->d_species_source[MAX_VARIABLE_TYPES + 1] = dQ2dx1;
          mp->d_species_source[MAX_VARIABLE_TYPES + 2] = dQ2dx2;
          mp->d_species_source[MAX_VARIABLE_TYPES + 3] = 0.;
          mp->d_species_source[MAX_VARIABLE_TYPES + four] = 0.;
          mp->d_species_source[MAX_VARIABLE_TYPES + five] = dQ2dx5;
        }
        break;
      case 2:
        mp->species_source[species_no] = Q1 + Q2;

        var = TEMPERATURE;
        if (pd->v[pg->imtrx][var]) {
          mp->d_species_source[var] = 0.;
        }

        var = MASS_FRACTION;
        if (pd->v[pg->imtrx][var]) {
          mp->d_species_source[MAX_VARIABLE_TYPES + 0] = 0.;
          mp->d_species_source[MAX_VARIABLE_TYPES + 1] = dQ2dx1;
          mp->d_species_source[MAX_VARIABLE_TYPES + 2] = dQ1dx2 + dQ2dx2;
          mp->d_species_source[MAX_VARIABLE_TYPES + 3] = dQ1dx3;
          mp->d_species_source[MAX_VARIABLE_TYPES + four] = 0.;
          mp->d_species_source[MAX_VARIABLE_TYPES + five] = dQ2dx5;
        }
        break;
      case 3:
        mp->species_source[species_no] = Q1 + Q3;

        var = TEMPERATURE;
        if (pd->v[pg->imtrx][var]) {
          mp->d_species_source[var] = 0.;
        }

        var = MASS_FRACTION;
        if (pd->v[pg->imtrx][var]) {
          mp->d_species_source[MAX_VARIABLE_TYPES + 0] = dQ3dx0;
          mp->d_species_source[MAX_VARIABLE_TYPES + 1] = 0.;
          mp->d_species_source[MAX_VARIABLE_TYPES + 2] = dQ1dx2;
          mp->d_species_source[MAX_VARIABLE_TYPES + 3] = dQ1dx3 + dQ3dx3;
          mp->d_species_source[MAX_VARIABLE_TYPES + four] = dQ3dx4;
          mp->d_species_source[MAX_VARIABLE_TYPES + five] = 0.;
        }
        break;
      case 4:
        mp->species_source[species_no] = -Q3;

        var = TEMPERATURE;
        if (pd->v[pg->imtrx][var]) {
          mp->d_species_source[var] = 0.;
        }

        var = MASS_FRACTION;
        if (pd->v[pg->imtrx][var]) {
          mp->d_species_source[MAX_VARIABLE_TYPES + 0] = -dQ3dx0;
          mp->d_species_source[MAX_VARIABLE_TYPES + 1] = 0.;
          mp->d_species_source[MAX_VARIABLE_TYPES + 2] = 0.;
          mp->d_species_source[MAX_VARIABLE_TYPES + 3] = -dQ3dx3;
          mp->d_species_source[MAX_VARIABLE_TYPES + four] = -dQ3dx4;
          mp->d_species_source[MAX_VARIABLE_TYPES + five] = 0.;
        }
        break;
      case 5:
        mp->species_source[species_no] = -Q2;

        var = TEMPERATURE;
        if (pd->v[pg->imtrx][var]) {
          mp->d_species_source[var] = 0.;
        }

        var = MASS_FRACTION;
        if (pd->v[pg->imtrx][var]) {
          mp->d_species_source[MAX_VARIABLE_TYPES + 0] = 0.;
          mp->d_species_source[MAX_VARIABLE_TYPES + 1] = -dQ2dx1;
          mp->d_species_source[MAX_VARIABLE_TYPES + 2] = -dQ2dx2;
          mp->d_species_source[MAX_VARIABLE_TYPES + 3] = 0.;
          mp->d_species_source[MAX_VARIABLE_TYPES + four] = 0.;
          mp->d_species_source[MAX_VARIABLE_TYPES + five] = -dQ2dx5;
        }
        break;
      case 6:
        mp->species_source[species_no] = 0.;

        var = TEMPERATURE;
        if (pd->v[pg->imtrx][var]) {
          mp->d_species_source[var] = 0.;
        }

        var = MASS_FRACTION;
        if (pd->v[pg->imtrx][var]) {
          mp->d_species_source[MAX_VARIABLE_TYPES + 0] = 0.;
          mp->d_species_source[MAX_VARIABLE_TYPES + 1] = 0.;
          mp->d_species_source[MAX_VARIABLE_TYPES + 2] = 0.;
          mp->d_species_source[MAX_VARIABLE_TYPES + 3] = 0.;
          mp->d_species_source[MAX_VARIABLE_TYPES + four] = 0.;
          mp->d_species_source[MAX_VARIABLE_TYPES + five] = 0.;
        }
        break;
      }
    }
  } // for imtrx
  return 0;
#endif // MAX_CONC < 7
}
/*****************************************************************************/
/* END of routine ion_reaction_source                                        */
/*****************************************************************************/

/*
   This routine calculates the electrolyte temperature as a function of time
   using a lumped energy-transport model.
   KSC: 10/24/98, based on an earlier routine written by GHE on 10/23/98;
   Revised: KSC 3/3/99 */

int electrolyte_temperature(double t,       /* present value of time */
                            double delta_t, /* present time step     */
                            int print)      /* key for printing:
                                                  print = 1: printing;
                                                  print = 0: no printing.  */
{
  /* Local Variables */

  int status = 0;
  FILE *fp;
  dbl T; /* electrolyte temperature in K - Joule heating neglected */
  dbl t0 = 0;
  dbl T0;    /* initial electrolyte temperature in K */
  dbl Ta;    /* ambient temperature in K */
  dbl A;     /* cross section area from which heat is lost to the ambient in m^2 */
  dbl h0;    /* heat transfer coefficient W/m^2/K */
  dbl m;     /* mass of the battery cell in kg */
  dbl Cp;    /* heat capacity of electrolyte in J/kg/K */
  dbl AhmCp; /* = 2Ah0/m/Cp with a unit of time */

  T0 = mp->u_solution_temperature[0];
  Ta = mp->u_solution_temperature[1];
  A = mp->u_solution_temperature[2];
  h0 = mp->u_solution_temperature[3];
  m = mp->u_solution_temperature[4];
  Cp = mp->u_solution_temperature[5];

  AhmCp = 2.0 * A * h0 / (m * Cp);

  T = Ta + (T0 - Ta) * exp(-AhmCp * t); /* neglect Joule heating */

  if (print == 1) /* when so desired, write out electrolyte temperature to a file, T_vs_t.out  */
  {
    if (t == 0.0 || t == delta_t) {
      t0 = delta_t;
    } else if (t == tran->init_time || t == (tran->init_time + delta_t)) {
      t0 = tran->init_time;
    }
#ifndef tflop
    if (t <= t0) {
      fp = fopen("T_vs_t.out", "w");
    } else {
      fp = fopen("T_vs_t.out", "a");
    }
#else
    fp = fopen("T_vs_t.out", "a");
#endif
    if (fp != NULL) {
      fprintf(fp, "%15e %15e\n", t, T);
    } else {
      printf("Unable to open output file, T_vs_t.out ");
    }
    fclose(fp);
  } else if (print == 0) {
    mp->electrolyte_temperature = T; /* put electrolyte temperature in the mp structure */
  }
  return (status);
}
/*****************************************************************************/
/* END of routine electrolyte_temperature                                    */
/*****************************************************************************/

/*  _______________________________________________________________________  */

/* assemble_suspension_temperature -- assemble terms (Residual &| Jacobian) for energy eqns
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

/******************************************************************************
 *
 *  Routine which calculates the bond evolution through an evolution equation.
 *  It is called suspension temperature for historical reasons
 *     Author: R.R. Rao, 10/02
 *
 ******************************************************************************/

int assemble_bond_evolution(double time, /* present time value */
                            double tt,   /* parameter to vary time integration from
                                            explicit (tt = 1) to implicit (tt = 0) */
                            double dt /* current time step size */) {
  /*local variables */
  int eqn, var, peqn, pvar;
  int a, b, p, i, j;
  int status = 1;
  int dim;

  dbl nn;     /* number of bonds */
  dbl nn_dot; /* bond derivative wrt time. */

  dbl gterm_a, gterm_b;
  dbl d_gterm_a, d_gterm_b;

  dbl mass; /* For terms and their derivatives */

  dbl advection, advection_a, advection_b, advection_c, advection_d;

  dbl diffusion;
  dbl diff_a, diff_b, diff_c, diff_d;
  dbl diff_coef;

  dbl source;
  dbl source1, source2;
  dbl d_source1, d_source2;
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
  dbl grad_phi_j[DIM];

  dbl h3;          /* Volume element (scale factors). */
  dbl dh3dmesh_bj; /* Sensitivity to (b,j) mesh dof. */

  dbl det_J;

  dbl d_det_J_dmeshbj;        /* for specified (b,j) mesh dof */
  dbl dgrad_phi_i_dmesh[DIM]; /* ditto.  */
  dbl wt;

  /* adjustable parameters for the bond model */
  dbl k1, k2, n0, aexp, bexp;

  dbl gamma_dot[DIM][DIM];

  dbl gammadot; /* strain rate invariant */

  dbl d_gd_dv[DIM][MDE];    /* derivative of strain rate invariant
                               wrt velocity */
  dbl d_gd_dmesh[DIM][MDE]; /* derivative of strain rate invariant
                               wrt mesh */
  dbl offset = DBL_SMALL;

  for (i = 0; i < DIM; i++) {
    grad_phi_i[i] = 0;
  }

  /* Compute gamma_dot[][] */

  /* Compute gammadot, grad(gammadot), gamma_dot[][], d_gd_dG, and d_grad_gd_dG */

  eqn = R_BOND_EVOLUTION;

  /*
   * Bail out fast if there's nothing to do...
   */
  if (!pd->e[pg->imtrx][eqn]) {
    return (status);
  }

  dim = ei[pg->imtrx]->ielem_dim;

  wt = fv->wt; /* Gauss point weight. */

  h3 = fv->h3; /* Differential volume element. */

  det_J = bf[eqn]->detJ;

  memset(gamma_dot, 0, DIM * DIM * sizeof(dbl));

  for (a = 0; a < VIM; a++) {
    for (b = 0; b < VIM; b++) {
      gamma_dot[a][b] = fv->grad_v[a][b] + fv->grad_v[b][a];
    }
  }

  calc_shearrate(&gammadot, gamma_dot, d_gd_dv, d_gd_dmesh);

  nn = fv->nn;
  if (pd->TimeIntegration != STEADY) {
    nn_dot = fv_dot->nn;
  } else {
    nn_dot = 0.0;
  }

  /* clip negative values
  if(nn < 1.e-5) nn = 0.;*/
  nn = MAX(nn, 0.);

  k2 = gn->k2;
  k1 = gn->k1;

  n0 = gn->n0;
  aexp = gn->pexp;
  bexp = gn->qexp;

  diff_coef = gn->diff;

  gterm_a = pow(gammadot + offset, aexp);
  gterm_b = pow(gammadot + offset, bexp);
  /*  gterm_b =1.;*/

  d_gterm_a = aexp * pow(gammadot + offset, aexp - 1.);
  d_gterm_b = bexp * pow(gammadot + offset, bexp - 1.);
  /*   d_gterm_b =0.; */

  /*
   * Residuals___________________________________________________________
   */

  /* equation is from Mujumdar et al., (102, J. Non-Newt. Fluid Mech., 2002)
   * and takes the form:
   * -dn/dt=k1*n*(shear-rate invariant)**a -
   *          k2*(n0-n)*(shear-rate invariant)**b
   */

  if (af->Assemble_Residual) {
    peqn = upd->ep[pg->imtrx][eqn];
    var = BOND_EVOLUTION;
    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      phi_i = bf[eqn]->phi[i];
      for (p = 0; p < VIM; p++) {
        grad_phi_i[p] = bf[var]->grad_phi[i][p];
      }

      mass = 0.;
      if (pd->TimeIntegration != STEADY) {
        if (pd->e[pg->imtrx][eqn] & T_MASS) {
          mass = nn_dot;
          mass *= phi_i * det_J * wt * h3;
          mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
        }
      }
      advection = 0.;
      if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {

        for (p = 0; p < VIM; p++) {
          advection += (fv->v[p] - fv_dot->x[p]) * fv->grad_nn[p];
        }

        advection *= phi_i * wt * det_J * h3;
        advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
      }

      diffusion = 0.;
      if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {
        for (p = 0; p < VIM; p++) {
          diffusion += grad_phi_i[p] * fv->grad_nn[p];
        }
        diffusion *= diff_coef * det_J * wt * h3;
        diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
      }

      source = 0.;
      if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
        if (DOUBLE_ZERO(nn)) {
          source = -k2 * (n0)*gterm_b;
          source *= phi_i * det_J * wt * h3;
          source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
        } else {
          source = k1 * nn * gterm_a - k2 * (n0 - nn) * gterm_b;
          source *= phi_i * det_J * wt * h3;
          source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
        }
      }

      lec->R[LEC_R_INDEX(peqn, i)] += mass + advection + diffusion + source;
    }
  }

  /*
   * Jacobian terms...
   */

  if (af->Assemble_Jacobian) {
    peqn = upd->ep[pg->imtrx][eqn];
    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {

      phi_i = bf[eqn]->phi[i];
      for (p = 0; p < VIM; p++) {
        grad_phi_i[p] = bf[eqn]->grad_phi[i][p];
      }

      /*
       * Set up some preliminaries that are needed for the (a,i)
       * equation for bunches of (b,j) column variables...
       */

      /*
       * J_nn_nn
       */
      var = BOND_EVOLUTION;
      if (pd->v[pg->imtrx][var]) {
        pvar = upd->vp[pg->imtrx][var];
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          phi_j = bf[var]->phi[j];
          for (p = 0; p < VIM; p++) {
            grad_phi_j[p] = bf[var]->grad_phi[j][p];
          }

          mass = 0.;
          if (pd->TimeIntegration != STEADY) {
            /*       printf("\t mass = %d and pd->TimeIntegration= %d \n", mass,
             * pd->TimeIntegration); */
            if (pd->e[pg->imtrx][eqn] & T_MASS) {
              mass = (1. + 2. * tt) * phi_j / dt;
              /* printf("\t mass = %d and dt= %d and tt= %d\n", mass, dt, tt); */

              mass *= phi_i * det_J * wt * h3;

              /* printf("\t det_J = %d and wt= %d and h3= %d\n", det_J, wt, h3); */
              mass *= pd->etm[pg->imtrx][eqn][(LOG2_MASS)];
            }
          }
          /*  printf("finished mass term\n"); */

          advection = 0.;
          if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
            for (p = 0; p < VIM; p++) {
              advection += (fv->v[p] - fv_dot->x[p]) * grad_phi_j[p];
            }
            advection *= phi_i * det_J * wt * h3;
            advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
          }

          diffusion = 0.;
          if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {
            for (p = 0; p < VIM; p++) {
              diffusion += grad_phi_j[p] * grad_phi_i[p];
            }
            diffusion *= diff_coef * det_J * wt * h3;
            diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
          }

          source = 0.;
          if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
            if (DOUBLE_ZERO(nn)) {
              source = k2 * gterm_b * phi_j;
              source *= phi_i * det_J * wt * h3;
              source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
            } else {
              source = (k1 * gterm_a + k2 * gterm_b) * phi_j;
              source *= phi_i * det_J * wt * h3;
              source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
            }
          }

          lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += mass + advection + diffusion + source;
        }
      }

      /*
       * J_nn_v
       */
      for (b = 0; b < VIM; b++) {
        var = VELOCITY1 + b;
        if (pd->v[pg->imtrx][var]) {
          pvar = upd->vp[pg->imtrx][var];
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];
            grad_phi_j[b] = bf[var]->grad_phi[j][b];

            mass = 0.;

            advection = 0.;
            if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
              advection += phi_i * bf[var]->grad_phi[j][b] * fv->grad_nn[b];
              advection *= det_J * wt * h3;
              advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
            }

            source = 0.;
            if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
              if (DOUBLE_ZERO(nn)) {
                source2 = k2 * (n0)*d_gterm_b;
                source -= (source2)*d_gd_dv[b][j];
                source *= phi_i * det_J * wt * h3;
                source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
              } else {
                source1 = k1 * nn * d_gterm_a;
                source2 = k2 * (n0 - nn) * d_gterm_b;
                source += (source1 - source2) * d_gd_dv[b][j];
                source *= phi_i * det_J * wt * h3;
                source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
              }
            }

            lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += mass + advection + source;
          }
        }
      }

      /*
       * J_nn_d
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
            if (pd->TimeIntegration != STEADY) {
              if (pd->e[pg->imtrx][eqn] & T_MASS) {
                mass += nn_dot;
                mass *= phi_i * wt * (h3 * d_det_J_dmeshbj + dh3dmesh_bj * det_J);
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
               *
               */

              advection_a = 0.;
              for (p = 0; p < dim; p++) {
                advection_a += (fv->v[p] - fv_dot->x[p]) * fv->d_grad_nn_dmesh[p][b][j];
              }
              advection_a *= phi_i * h3 * det_J * wt;

              advection_b = 0.;
              for (p = 0; p < dim; p++) {
                advection_b += (fv->v[p] - fv_dot->x[p]) * fv->grad_nn[p];
              }
              advection_b *= phi_i * h3 * d_det_J_dmeshbj * wt;

              advection_c = 0.;
              if (pd->TimeIntegration != STEADY) {
                if (pd->e[pg->imtrx][eqn] & T_MASS) {

                  advection_c -= bf[var]->phi[i] * (1 + 2. * tt) / dt * fv->grad_nn[b];

                  advection_c *= phi_i * h3 * det_J * wt;
                }
              }

              advection_d = 0.;
              for (p = 0; p < dim; p++) {
                advection_d += (fv->v[p] - fv_dot->x[p]) * fv->grad_nn[p];
              }
              advection_d *= phi_i * dh3dmesh_bj * det_J * wt;

              advection = advection_a + advection_b + advection_c + advection_d;

              advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
            }

            diffusion = 0.;
            if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {
              diff_a = 0.;
              for (p = 0; p < dim; p++) {
                dgrad_phi_i_dmesh[p] = bf[eqn]->d_grad_phi_dmesh[i][p][b][j];

                diff_a += dgrad_phi_i_dmesh[p] * fv->grad_nn[p];
              }
              diff_a *= det_J * h3 * wt;

              diff_b = 0.;
              for (p = 0; p < VIM; p++) {
                diff_b += fv->d_grad_nn_dmesh[p][b][j] * grad_phi_i[p];
              }
              diff_b *= det_J * h3 * wt;

              diff_c = 0.;
              for (p = 0; p < dim; p++) {
                diff_c += grad_phi_i[p] * fv->grad_nn[p];
              }
              diff_c *= d_det_J_dmeshbj * h3 * wt;

              diff_d = 0.;
              for (p = 0; p < dim; p++) {
                diff_d += grad_phi_i[p] * fv->grad_nn[p];
              }
              diff_d *= det_J * dh3dmesh_bj * wt;

              diffusion = diff_a + diff_b + diff_c + diff_d;

              diffusion *= diff_coef * pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
            }

            source = 0.;

            if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
              if (nn <= 0.) {
                source2 = k2 * (n0 - nn) * gterm_b;

                d_source2 = k2 * (n0 - nn) * d_gterm_b;

                source -= (source2) * (d_det_J_dmeshbj * h3 + det_J * dh3dmesh_bj) +
                          (d_source2)*d_gd_dmesh[b][j] * h3 * det_J;

                source *= phi_i * wt * pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
              } else {
                source1 = k1 * nn * gterm_a;
                source2 = k2 * (n0 - nn) * gterm_b;

                d_source1 = k1 * nn * d_gterm_a;
                d_source2 = k2 * (n0 - nn) * d_gterm_b;

                source += (source1 - source2) * (d_det_J_dmeshbj * h3 + det_J * dh3dmesh_bj) +
                          (d_source1 - d_source2) * d_gd_dmesh[b][j] * h3 * det_J;

                source *= phi_i * wt * pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
              }
            }

            lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += mass + advection + diffusion + source;
          }
        }
      }
    }
  }

  return (status);
} /* end of assemble_bond_evolution*/
#ifdef NOT_USED
/*****************************************************************************/
static void cal_current_density(double x[],           /* global nodal solution vector  */
                                double time,          /* present time value            */
                                double delta_t,       /* present time step             */
                                int nsID[],           /* node set number               */
                                double local_value[]) /* values of current density     */

/******************************************************************************
 *
 *  A function that outputs electrical current density at the anode and
 *  cathode current collectors for a thermal-battery cell
 *
 *  Ken S. Chen (10/2004)
 *
 *
 ******************************************************************************/
{
  double i_acc;         /* current density at anode current collector, A/cm^2 */
  double i_ccc;         /* current density at cathode current collector, A/cm^2 */
  double T;             /* electrolyte or cell temperature, K */
  double sigma_anode;   /* anode electrical conductivity, S/cm */
  double sigma_cathode; /* cathode electrical conductivity, S/cm */
  int nodew_acc;        /* node number at anode current collector */
  int nodee_acc;        /* node number adjacent to anode current collector */
  int nodew_ccc;        /* node number adjacent to cathode current collector */
  int nodee_ccc;        /* node number at cathode current collector */

  int nsw_acc, nse_acc, nsw_ccc, nse_ccc;

  double PHI1e_acc, PHI1w_acc;      /* electrode potential in acc, V */
  double PHI1e_ccc, PHI1w_ccc;      /* electrode potential in ccc, V */
  int idx_PHI1w_acc, idx_PHI1e_acc; /* acc potential indices */
  int idx_PHI1w_ccc, idx_PHI1e_ccc; /* ccc potential indices */
  double xw_acc, xe_acc;            /* locations of acc and its neighbor */
  double xw_ccc, xe_ccc;            /* locations of ccc and its neighbor */

  sigma_anode = mp_glob[0]->thermal_conductivity;   /* anode conductivity, S/cm*/
  sigma_cathode = mp_glob[2]->thermal_conductivity; /* cathode conductivity, S/cm */

  nsw_acc = nsID[0];
  nse_acc = nsID[1];
  nsw_ccc = nsID[2];
  nse_ccc = nsID[3];
  nodew_acc = nsid2nn(nsw_acc); /* node number at acc current collector */
  nodee_acc = nsid2nn(nse_acc); /* node number adjacent to acc */
  nodew_ccc = nsid2nn(nsw_ccc); /* node number adjacent to ccc */
  nodee_ccc = nsid2nn(nse_ccc); /* node number at ccc current collector */

  if (nodew_acc == -1 || nodee_acc == -1 || nodew_ccc == -1 || nodee_ccc == -1)
    return;

  idx_PHI1w_acc = Index_Solution(nodew_acc, TEMPERATURE, 0, 0, -1, pg->imtrx);
  idx_PHI1e_acc = Index_Solution(nodee_acc, TEMPERATURE, 0, 0, -1, pg->imtrx);
  idx_PHI1w_ccc = Index_Solution(nodew_ccc, TEMPERATURE, 0, 0, -1, pg->imtrx);
  idx_PHI1e_ccc = Index_Solution(nodee_ccc, TEMPERATURE, 0, 0, -1, pg->imtrx);

  if (idx_PHI1w_acc == -1 || idx_PHI1e_acc == -1 || idx_PHI1w_ccc == -1 || idx_PHI1e_ccc == -1)
    return;

  PHI1w_acc = x[idx_PHI1w_acc];
  PHI1e_acc = x[idx_PHI1e_acc];
  PHI1w_ccc = x[idx_PHI1w_ccc];
  PHI1e_ccc = x[idx_PHI1e_ccc];
  xw_acc = Coor[0][nodew_acc];
  xe_acc = Coor[0][nodee_acc];
  xw_ccc = Coor[0][nodew_ccc];
  xe_ccc = Coor[0][nodee_ccc];

  i_acc = -sigma_anode * ((PHI1e_acc - PHI1w_acc) / (xe_acc - xw_acc));
  i_ccc = -sigma_cathode * ((PHI1e_ccc - PHI1w_ccc) / (xe_ccc - xw_ccc));

  electrolyte_temperature(time, delta_t, 0); /* get cell temperature */
  T = mp->electrolyte_temperature;

  local_value[0] = i_acc; /* return anode current density to rf_solve.c */
  local_value[1] = i_ccc; /* return cathode current density to rf_solve.c */

  fprintf(stderr, "In cal_current_density, time =%g\n", time);
  fprintf(stderr, "In cal_current_density: T =%g\n", T);
  fprintf(stderr, "In cal_current_density: i_acc =%g\n", i_acc);
  fprintf(stderr, "In cal_current_density: i_ccc =%g\n", i_ccc);

} /* END of routine cal_current_density */
/*****************************************************************************/
#endif

int etching_KOH_source(int wspec,     /* Current species number */
                       double *param) /* pointer to source model parameter list */
/******************************************************************************
*
*  A function that calculate source (or sink) terms in species equation (shell only for now)
*
*  Output:
*     mp->species_source[w]-------  Source/sink rate of species number w
*     mp->d_species_source[var] --- Sesntivity of the source/sink of species w
*                                   w.r.t. variables var
*
*  Right now, it works, ONLY AND IF ONLY you use the following:
*
*  Species type = SPECIES_DENSITY
*  Concentration Units: CGS, i.e. g/cm^3
*
*  Species ordering:
*                   0: H2O - water
*                   1: KOH - potassium hydroxide
*                   2: H2 - hydrogen
*                   3: Silicon hydroxyl byproducts

*
*  Kristianto Tjiptowidjojo (10/2018)
*
******************************************************************************/
{
  double etch_rate = 0.0;
  double d_etch_rate_d_C[2] = {0.0};

  /* Bulk density of crystalline silicon (g/cm^3) */
  double rho_bulk_Si = 2.3290;

  /* Molecular weight in mole/g */
  double MW_H2O = 18.01528;
  double MW_OH = 17.008;
  double MW_Si = 28.0855;
  double MW_H2 = (2.0 * 1.00794);
  double MW_SiO2OH2 = (28.0855 + 2.0 * 15.9994 + 2.0 * 17.008);

  /* Get mass concentration of each species
     Mass concentration unit is g/cm^3 */
  double rho_H2O = fv->c[0];
  double rho_KOH = fv->c[1];

  /* Area fraction, read from external field when available */
  double a_frac = 1.0;
  int i_ext_field = -1;
  if ((efv->ev) && (mp->SpeciesSourceModel[wspec] == ETCHING_KOH_EXT)) {
    i_ext_field = mp->species_source_external_field_index;
    if (i_ext_field < 0)
      GOMA_EH(-1, "Trouble getting external field index in ETCHING_KOH_EXT");
    a_frac = fv->external_field[i_ext_field];
  }

  /* Right now it only handles KOH wet etching on plane 100 of crystalline silicon*/
  etch_rate = calc_KOH_Si_etch_rate_100(rho_H2O, rho_KOH, d_etch_rate_d_C);

  /* Export it to mp->species_source, depending on their stochiometric coefficient */
  switch (wspec) {
  case 0: /* Water */
    mp->species_source[wspec] = a_frac * 2.0 * rho_bulk_Si / MW_Si * MW_H2O * etch_rate;
    break;

  case 1: /* OH */
    mp->species_source[wspec] = a_frac * 2.0 * rho_bulk_Si / MW_Si * MW_OH * etch_rate;
    break;

  case 2: /* H2 */
    mp->species_source[wspec] = a_frac * -2.0 * rho_bulk_Si / MW_Si * MW_H2 * etch_rate;
    break;

  case 3: /* SiO2OH2 */
    mp->species_source[wspec] = a_frac * -1.0 * rho_bulk_Si / MW_Si * MW_SiO2OH2 * etch_rate;
    break;
  }

  /* Export sensitivity to mp->d_species_source */
  switch (wspec) {
  case 0: /* Water */
    mp->d_species_source[MAX_VARIABLE_TYPES + 0] =
        a_frac * 2.0 * rho_bulk_Si / MW_Si * MW_H2O * d_etch_rate_d_C[0];
    mp->d_species_source[MAX_VARIABLE_TYPES + 1] =
        a_frac * 2.0 * rho_bulk_Si / MW_Si * MW_H2O * d_etch_rate_d_C[1];
    break;

  case 1: /* KOH */
    mp->d_species_source[MAX_VARIABLE_TYPES + 0] =
        a_frac * 2.0 * rho_bulk_Si / MW_Si * MW_OH * d_etch_rate_d_C[0];
    mp->d_species_source[MAX_VARIABLE_TYPES + 1] =
        a_frac * 2.0 * rho_bulk_Si / MW_Si * MW_OH * d_etch_rate_d_C[1];
    break;

  case 2: /* H2 */
    mp->d_species_source[MAX_VARIABLE_TYPES + 0] =
        a_frac * -2.0 * rho_bulk_Si / MW_Si * MW_H2 * d_etch_rate_d_C[0];
    mp->d_species_source[MAX_VARIABLE_TYPES + 1] =
        a_frac * -2.0 * rho_bulk_Si / MW_Si * MW_H2 * d_etch_rate_d_C[1];
    break;

  case 3: /* SiO2OH2 */
    mp->d_species_source[MAX_VARIABLE_TYPES + 0] =
        a_frac * -1.0 * rho_bulk_Si / MW_Si * MW_SiO2OH2 * d_etch_rate_d_C[0];
    mp->d_species_source[MAX_VARIABLE_TYPES + 1] =
        a_frac * -1.0 * rho_bulk_Si / MW_Si * MW_SiO2OH2 * d_etch_rate_d_C[1];
    break;
  }

  mp->d_species_source[MAX_VARIABLE_TYPES + 2] = 0.0;
  mp->d_species_source[MAX_VARIABLE_TYPES + 3] = 0.0;

  return 0;

} /* END of calc_KOH_Si_etch_rate_100 */

double calc_KOH_Si_etch_rate_100(double rho_H2O,            /* Concentration of water */
                                 double rho_KOH,            /* Concentration of KOH */
                                 double d_etch_rate_d_C[2]) /* Sensitivity of etch rate w.r.t.
                                                              concentration of water, and KOH */
/******************************************************************************
 *
 *  A function that outputs KOH wet etch rate of silicon surface (100 plane for now)
 *  based on kinetic model proposed by
 *
 *   Seidel, H., et al."Anisotropic etching of crystalline silicon in alkaline solutions I.
 *                      Orientation dependence and behavior of passivation layers."
 *   Journal of the electrochemical society 137.11 (1990): 3612-3626.
 *
 *
 *  Kinetic model is listed in Equation A-1
 *
 *  etch_rate = k0 * conc_H2O^4 * conc_KOH^0.25 * exp(-Ea/Kb T)
 *
 *  UNITS:
 *
 *  Kinetic models used above required units as follow:
 *
 *  Etch rate: micron/hour
 *  Rate constant k0: (micron/hr) (mole/liter)^-4.25
 *  Species concentration: mole/liter
 *
 *  Right now, unit conversion is handled automatically, ONLY AND IF ONLY
 *  you use the following:
 *
 *  Species type = SPECIES_DENSITY
 *  Concentration Units: CGS, i.e. g/cm^3
 *  Species ordering:
 *                   0: H2O - water
 *                   1: KOH - potassium hydroxide
 *                   2: H2 - hydrogen
 *                   3: Silicon hydroxyl byproducts
 *  Temperature: From Process Temperature card in general specification
 *
 *  Kristianto Tjiptowidjojo (5/2017)
 *
 *
 ******************************************************************************/
{
  double etch_rate, d_etch_rate_d_H2O, d_etch_rate_d_KOH;

  /* Boltzmann constant in eV/K */
  double k_B = 8.6173305e-5;

  /* Activation energy in eV */
  double E_a = 0.595;

  /* Temperature in K */
  double T = upd->Process_Temperature;

  /* Rate constant in (micron/hr) (mole/liter)^-4.25  */
  double k0 = 2480.0;

  /* Molecular weight in mole/g */
  double MW_H2O = 18.01528;
  double MW_KOH = 56.1056;

  /* Mole concentration in mol/liter */
  double C_H2O = rho_H2O * 1000.0 / MW_H2O;
  double C_KOH = rho_KOH * 1000.0 / MW_KOH;

  /* Evaluate heaviside function based on minimum concentration */
  double Hside = 1.0;
  double dHside_drho_KOH = 0.0;
  double rho_KOH_min = 1.0e-6;
  double rho_KOH_max = 1.0e-4;
  double width = rho_KOH_max - rho_KOH_min;
  double alpha = 0.5 * width;
  double rho_KOH_center = rho_KOH_max - alpha;
  double rho_KOH_normalized = rho_KOH - rho_KOH_center;

  if (rho_KOH >= rho_KOH_max) {
    Hside = 1.0;
    dHside_drho_KOH = 0.0;
  } else if (rho_KOH <= rho_KOH_min) {
    Hside = 0.0;
    dHside_drho_KOH = 0.0;
  } else {
    Hside =
        0.5 * (1. + rho_KOH_normalized / alpha + sin(M_PIE * rho_KOH_normalized / alpha) / M_PIE);
    dHside_drho_KOH = 0.5 * (1.0 / alpha + cos(M_PIE * rho_KOH_normalized / alpha) / alpha);
  }

  /* Calculate etch rate (micron/hr) */
  if (rho_KOH > rho_KOH_min) {
    etch_rate = Hside * k0 * pow(C_H2O, 4.0) * pow(C_KOH, 0.25) * exp(-E_a / k_B / T);
  } else {
    etch_rate = 0.0;
  }

  /* Convert to cm/s */
  etch_rate = etch_rate / 1.0e4 / 3600.0;

  /* Calculate sensitivity of etch rate w.r.t. concentration */

  if (d_etch_rate_d_C != NULL) {
    if (rho_KOH > rho_KOH_min) {

      d_etch_rate_d_H2O = Hside * 4.0 * k0 * pow(C_H2O, 3.0) * pow(C_KOH, 0.25) *
                          exp(-E_a / k_B / T) / (1.0e4 * 3600.0) * (1000.0 / MW_H2O);

      d_etch_rate_d_KOH = Hside * 0.25 * k0 * pow(C_H2O, 4.0) / pow(C_KOH, 0.75) *
                          exp(-E_a / k_B / T) / (1.0e4 * 3600.0) * (1000.0 / MW_KOH);

      d_etch_rate_d_KOH += k0 * pow(C_H2O, 4.0) * pow(C_KOH, 0.25) * exp(-E_a / k_B / T) / 1.0e4 /
                           3600.0 * dHside_drho_KOH;
    } else {
      d_etch_rate_d_H2O = 0.0;
      d_etch_rate_d_KOH = 0.0;
    }

    /* Export the etch rate and its sensitivities */
    d_etch_rate_d_C[0] = d_etch_rate_d_H2O;
    d_etch_rate_d_C[1] = d_etch_rate_d_KOH;

  } /* End of if calculating Jacobian entries */
  return etch_rate;

} /* END of calc_KOH_Si_etch_rate_100 */
