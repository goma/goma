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

/* Standard include files */

#include <math.h>
#include <stdio.h>
#include <string.h>

/* GOMA include files */

#include "mm_fill_potential.h"

#include "density.h"
#include "el_elm.h"
#include "mm_as.h"
#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "mm_eh.h"
#include "mm_fill_ls.h"
#include "mm_fill_solid.h"
#include "mm_fill_species.h"
#include "mm_fill_terms.h"
#include "mm_mp.h"
#include "mm_mp_const.h"
#include "mm_mp_structs.h"
#include "mm_qtensor_model.h"
#include "rf_bc_const.h"
#include "rf_fem.h"
#include "rf_fem_const.h"
#include "std.h"
#include "user_mp.h"

#define GOMA_MM_FILL_POTENTIAL_C

/*********** R O U T I N E S   I N   T H I S   F I L E ***********************
 *
 *       NAME			TYPE			CALLED BY
 *  -----------------------------------------------------------------
 *  assemble_potential    	void	             matrix_fill()
 *  assemble_Enorm               int                  matrix_fill()
 *  surface_charge_surf          void                 apply_integrated_bc()
 *  current_BV_surf              void                 apply_integrated_bc()
 *  current_ORR_surf             void                 apply_integrated_bc()
 *  apply_potential_grad_bc      void                 apply_integrated_bc()
 ******************************************************************************/

/*  _______________________________________________________________________  */

/* assemble_current energy -- assemble terms (Residual & Jacobian) for scalar
 * transport current eqn
 *
 * in:
 * 	ei -- pointer to Element Indices	structure
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
 * Revised:     10/98 by KSC  (to enable calculation of electrical conductivity and source term
 *              that are functions of species mole fractions and electrical potentials for
 *              for modeling electrochemical processes such as thermal batteries).
 * Revised:     1/99 by KSC (to add the mass term to the residual and jacobian equations)
 * Revised:     4/00 by RSL (to enable computation of Jacobian terms with respect to mole
 *              fractions in the Poisson equation describing charged species transport)
 * Revised:     8/2000 by KSC (to clean up and fix minor bugs)
 * Revised:     9/2000 by KSC (to implement electrolyte conductivity and electrical current density
                due to diffusion of charged species in the Poisson equation governing electrical
                potential of dilute electrolyte in which solute species fluxes can be described
                by Fick's first law and in which electroneutrality is valid everywhere).
 * Revised:     10/2000 by KSC (to add missing jacobian terms and current source term due to
                homogeneous reactions, and to multiply k and ii2 by Faraday's constant F
                so that they have the proper units as electrolyte conductivity and
                current density, respectively, in the case of Fickian flux model).
 * Revised:     12/2001 by KSC (to cleaned up calls to the electrode_species_source routine)
 * Revised:     4/2002 by KSC (to add ELECTRONEUTRALITY_FICKIAN and ELECTRONEUTRALITY_SM
 *                             flags for the electrical conductivity models that are
 *                             associated with these names).
 * Revised:     5/2002 by KSC (to add a NET_CHARGE current source card)
 * Revised:     7/2003 by ACS (add a new current source card to set up DEBYE-HUCKLE approx.)
 * Revised:     7/2004 by KSC (to modify the electrical conductivity and
                               current source terms so that the
                               electrolyte-potential equation can be
                               solved in its original form -- see
                               Eq. 2.13 of thermal-battery SAND report).
 * Revised:     1/2006 by KSC (to add routines for computing current densities
 *                             using, respectively, linearized and Tafel
 *                             kinetic models for electrochemical rxns
 *                             such as the hydrogen oxidation and the oxygen
 *                             reduction rxns in PEM fuel cells).
 *
 *
 * Note: currently we do a "double load" into the addresses in the global
 *       "a" matrix and resid vector held in "esp", as well as into the
 *       local accumulators in "lec".
 *
 */
/*ARGSUSED*/
int assemble_potential(double time, /* present time value */
                       double tt,   /* parameter to vary time integration from
                                     * explicit (tt = 1) to implicit (tt = 0)    */
                       double dt)   /* current time step size */
{
  int eqn, var, peqn, pvar, var_offset, dim, p, a, b, w;
  int i = -1;
  int j, status;

  dbl ii[DIM];     /* current density due to electrical potential gradient */
  dbl ii2[DIM];    /* current density due to concentration
                      gradients of charged species */
  dbl grad_V[DIM]; /* voltage potential gradient. */

  dbl k = 1e12;            /* electrical conductivity. */
  dbl dkdT[MDE];           /* Temperature derivative of electrical conductivity. */
  dbl dkdV[MDE];           /* Voltage derivative of electrical conductivity. */
  dbl dkdC[MAX_CONC][MDE]; /* Concentration derivative of electrical conductivity. */
  dbl dkdX[DIM][MDE];      /* Spatial derivatives of t.c. */

  dbl h = 1e12;            /* current source. */
  dbl dhdT[MDE];           /* Temperature derivative of c.s. */
  dbl dhdC[MAX_CONC][MDE]; /* Concentration derivative of c.s. */
  dbl dhdX[DIM][MDE];      /* Spatial derivative of c.s. */
  dbl dhdv[DIM][MDE];      /* Velocity derivative for the c.s. */
  /* The following declaration appears to be incorrect, so I have replaced it -- RSL 3/28/00 */
  /* dbl dhdV[DIM];	*/ /* Voltage derivative for the c.s. */
  dbl dhdV[MDE];           /* Voltage derivative for the c.s. */

  dbl didmesh[DIM]; /* At a specific (b,j) mesh dof. */
  dbl diffusion, diff_a, diff_b, diff_c, diff_d, source;
  dbl lambda_d;

  /*
   * Galerkin weighting functions for i-th energy residuals
   * and some of their derivatives...
   */
  dbl phi_i;
  dbl grad_phi_i[DIM];

  /*
   * Interpolation functions for variables and some of their derivatives.
   */
  dbl grad_phi_j[DIM];
  dbl d_i_V_j[DIM];

  dbl h3;          /* Volume element (scale factors). */
  dbl dh3dmesh_bj; /* Sensitivity to (b,j) mesh dof. */
  dbl det_J;

  dbl d_det_J_dmeshbj;        /* for specified (b,j) mesh dof */
  dbl dgrad_phi_i_dmesh[DIM]; /* ditto.  */
  dbl wt;

  dbl A_inv[MAX_CONC * DIM]
           [MAX_CONC * DIM]; /* inverse of Stefan_Maxwell coefficient matrix; KSC: 10/28/98 */
  dbl *D = mp->diffusivity;  /* Fickian diffusion coefficients */
  dbl T = 298.15;            /* Electrolyte solution temperature default = 25 C */
  dbl c = 0.0;               /* Total molar concentration of electrolyte solution */
                             /* ERROR -> appears to be an error because this is never set!!!! */
  dbl rho = 0.0;             /* Mass density of electrolyte solution */
                             /* ERROR -> appears to be an error because this is never set!!!! */
  DENSITY_DEPENDENCE_STRUCT d_rho_struct; /* density dependence */
  DENSITY_DEPENDENCE_STRUCT *d_rho = &d_rho_struct;
  dbl *M = mp->molecular_weight; /* Ptr to species molecular weights */
  dbl M_mix = 0.0;               /* Molecular weight of electrolyte solution */
                                 /* ERROR -> appears to be an error because this is never set!!!! */
  dbl grad_x[MAX_CONC][DIM];     /* mole fraction gradient */

  dbl *z = mp->charge_number; /* ptr to the charge numbers for species */
  dbl sumdelx;                /* dummy variables */
  const double F = 96487.0;   /* Faraday's constant in units of C/equiv. */
  dbl h_term1;
  dbl dii2dC[DIM][MAX_CONC][MDE]; /* concentration derivative of ii2 -- RSL 3/28/00 */
  int n = -1, l, m, kk, w1;
  dbl sum1, sum2, sum3, sumzdrt, sumzdrv, sumzdrc;

  memset((void **)A_inv, 0, MAX_CONC * DIM * MAX_CONC * DIM * sizeof(double));
  memset((void **)d_rho, 0, sizeof(DENSITY_DEPENDENCE_STRUCT));

  /*
   * Set the number of species variables from the mp structure
   */
  int n_species = mp->Num_Species;

  /* initialize grad_phi_j, grad_phi_i */
  for (j = 0; j < DIM; j++) {
    grad_phi_i[j] = 0;
    grad_phi_j[j] = 0;
  }

  /*   static char yo[] = "assemble_current";*/

  status = 0;

  /*
   * Unpack variables from structures for local convenience...
   */
  dim = pd->Num_Dim;

  /*
   * Bail out fast if there's nothing to do...
   */
  eqn = R_POTENTIAL;
  if (!pd->e[pg->imtrx][eqn])
    return (status);

  wt = fv->wt;           /* Gauss point weight. */
  h3 = fv->h3;           /* Differential volume element. */
  det_J = bf[eqn]->detJ; /* Really, ought to be mesh eqn. */

  /*
   * Material property constants, etc. Any variations at this Gauss point
   * are calculated with user-defined subroutine options.
   *
   * Properties to consider here are rho, Cp, k, and h.  For now we will
   * take rho as constant.  Cp, h, and k we will allow to vary with temperature,
   * spatial coordinates, and species concentration.
   */

  /*** Electrical Conductivity moved to separate routine  ****/
  k = electrical_conductivity(dkdT, dkdV, dkdC, dkdX, time, dt);

  /*** Current Source ****/
  memset(dhdT, 0, sizeof(double) * MDE);
  memset(dhdV, 0, sizeof(double) * MDE);
  memset(dhdX, 0, sizeof(double) * DIM * MDE);
  memset(dhdC, 0, sizeof(double) * MAX_CONC * MDE);
  if (mp->CurrentSourceModel == USER) {
    (void)usr_current_source(mp->u_current_source);
    h = mp->current_source;

    var = TEMPERATURE;
    for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
      dhdT[j] = mp->d_current_source[var] * bf[var]->phi[j];
    }

    var = VOLTAGE;
    for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
      dhdV[j] = mp->d_current_source[var] * bf[var]->phi[j];
    }

    if (pd->v[pg->imtrx][VELOCITY1]) {
      for (a = 0; a < dim; a++) {
        var = VELOCITY1 + a;
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          dhdv[a][j] = mp->d_current_source[var] * bf[var]->phi[j];
        }
      }
    }

    if (pd->v[pg->imtrx][MESH_DISPLACEMENT1]) {
      for (a = 0; a < dim; a++) {
        var = MESH_DISPLACEMENT1 + a;
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          dhdX[a][j] = mp->d_current_source[var] * bf[var]->phi[j];
        }
      }
    }

    if (pd->v[pg->imtrx][MASS_FRACTION]) {
      for (w = 0; w < pd->Num_Species_Eqn; w++) {
        var = MASS_FRACTION;
        var_offset = MAX_VARIABLE_TYPES + w;
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          dhdC[w][j] = mp->d_current_source[var_offset] * bf[var]->phi[j];
        }
      }
    }

  } else if (mp->CurrentSourceModel ==
             ELECTRODE_KINETICS) /* KSC: 10/21/98; modified, KSC: 8/30/00 */
  {
    n_species = pd->Num_Species_Eqn + 1;
    n = n_species * VIM;

    if (mp->SolutionTemperatureModel == CONSTANT) /* constant solution temperature */
    {
      T = mp->solution_temperature;
    } else if (mp->SolutionTemperatureModel ==
               THERMAL_BATTERY) /* thermal battery temperature model */
    {
      electrolyte_temperature(time, dt, 0); /* calculate electrolyte temperature at present time */
      T = mp->electrolyte_temperature;
    } else {
      GOMA_EH(GOMA_ERROR,
              "Solution-temperature model other than THERMAL_BATTERY awaits future implementation");
    }
    /* set the solution temperature to the 298 K if it is zero - safety feature */
    if (T == 0.0)
      T = 298.0;

    h_term1 = 0.0;
    for (w = 0; w < n_species; w++) /*  upper limit changed -- RSL 8/7/00  */
    {
      electrode_species_source(w, time, dt);
      h_term1 += mp->charge_number[w] * mp->species_source[w];
    }
    /* h = h_term1; */ /* let h represent the first source term; KSC: 4/20/99 */
    h = F * h_term1;   /* h represent the first source term; KSC: 7-29-04 */

    /*  Redo calculations of dhdT and dhdV; previous versions apply only to linearized Butler-
        Volmer kinetics -- RSL 8/7/00  */
    sumzdrt = 0.;
    sumzdrv = 0.;
    for (w1 = 0; w1 < n_species; w1++) {
      electrode_species_source(w1, time, dt);
      sumzdrt += mp->charge_number[w1] * mp->d_species_source[TEMPERATURE];
      sumzdrv += mp->charge_number[w1] * mp->d_species_source[VOLTAGE];
    }
    var = TEMPERATURE;
    for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
      /* dhdT[j] = sumzdrt*bf[var]->phi[j]; */
      dhdT[j] = F * sumzdrt * bf[var]->phi[j]; /* KSC: 7-29-04 */
    }
    var = VOLTAGE;
    for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
      /* dhdV[j] = sumzdrv*bf[var]->phi[j]; */
      dhdV[j] = F * sumzdrv * bf[var]->phi[j]; /* KSC: 7-29-04 */
    }

    /* Derivatives of the scalar h with respect to mole fraction are computed here; h is the
       single sum in Equation (2.13) of the thermal battery LDRD report. -- RSL 8/8/00 */

    if (pd->v[pg->imtrx][MASS_FRACTION]) {
      for (w = 0; w < pd->Num_Species_Eqn; w++) {
        var = MASS_FRACTION;
        var_offset = MAX_VARIABLE_TYPES + w;
        sumzdrc = 0.;
        for (w1 = 0; w1 < n_species; w1++) {
          electrode_species_source(w1, time, dt);
          sumzdrc += mp->charge_number[w1] * mp->d_species_source[var_offset];
        }
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          /* dhdC[w][j] = sumzdrc*bf[var]->phi[j]; */
          dhdC[w][j] = F * sumzdrc * bf[var]->phi[j]; /* KSC: 7-29-04 */
        }

      } /* end of loop over w */

    } /* end of var = MASS_FRACTION */

  } /* end of CurrentSourceModel == ELECTRODE_KINETICS */

  else if (mp->CurrentSourceModel == CONSTANT) {
    h = mp->current_source;

    /*Sensitivities were already set to zero */
  } else if (mp->CurrentSourceModel == BUTLER_VOLMER) /* added by KSC: 04/28/06 */
  {
    dbl dh[3];
    int wspec = mp->u_current_source[0];

    /* Computing current source by calling butler_volmer_source */
    h = butler_volmer_source(mp->u_current_source, 1, dh);

    var = TEMPERATURE;
    for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
      dhdT[j] = dh[0] * bf[var]->phi[j];
    }

    var = VOLTAGE;
    for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
      dhdV[j] = dh[1] * bf[var]->phi[j];
    }

    if (pd->v[pg->imtrx][MASS_FRACTION]) {
      for (w = 0; w < pd->Num_Species_Eqn; w++) {
        var = MASS_FRACTION;
        dhdC[wspec][j] = dh[2] * bf[var]->phi[j];
      }
    }
  } else if (mp->CurrentSourceModel == FICKIAN_CHARGED ||
             mp->CurrentSourceModel == STEFAN_MAXWELL_CHARGED) /* added by KSC: 10/4/00 */
  {
    h = 0.0;
    for (w = 0; w < n_species; w++) {
      h += mp->charge_number[w] * mp->species_source[w];
    }
    h *= F;

  } else if (mp->CurrentSourceModel == NET_CHARGE) {
    h = 0.0;
    for (w = 0; w < pd->Num_Species_Eqn; w++) {
      h += mp->charge_number[w] * fv->c[w];
    }
    h *= F;

    if (pd->v[pg->imtrx][MASS_FRACTION]) {
      for (w = 0; w < pd->Num_Species_Eqn; w++) {
        var = MASS_FRACTION;
        var_offset = MAX_VARIABLE_TYPES + w;
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          dhdC[w][j] = F * mp->charge_number[w] * bf[var]->phi[j];
        }
      }
    }
  } else if (mp->CurrentSourceModel == DEBYE_HUCKEL)
  /* This is the implementation of resolving potential field in the
  region close of the EDL.  The Debye-Huckle approximation assumes
  small EDL relative to the length of flow dimension.  Also,
  this assumes a Boltzmann distribution of charges near the wall.
  This is developed for the use of electroosmotic flow. ACS 7/03 */
  {
    h = -fv->V;
    lambda_d = mp->u_current_source[0];
    h /= (lambda_d * lambda_d);

  } else {
    GOMA_EH(GOMA_ERROR, "Unrecognized current source model");
  }
  /********** End of specification of the current source model **********/

  for (a = 0; a < dim; a++) {
    grad_V[a] = fv->grad_V[a];
  }

  for (a = 0; a < dim; a++) {
    ii[a] = -k * grad_V[a];
  }

  /* Calculate ii2, current density due to concentration gradients of charged
     species (as versus ii, current density due to electrical potential gradients),
     for processes involving charged species and in which electroneutrality can be
     taken to be valid. KSC: (8/30/00, 9/13/00) and RSL 9/8/00 */

  if (cr->MassFluxModel == STEFAN_MAXWELL_CHARGED || /* Stefan-Maxwell diffusion of charged species
                                                        in concentrated electrolyte */
      cr->MassFluxModel == STEFAN_MAXWELL_VOLUME)    /* added by RSL 9/14/00 */
  {
    for (i = 0; i < n_species - 1; i++) {
      for (a = 0; a < VIM; a++) {
        grad_x[i][a] = fv->grad_c[i][a];
      }
    }
    for (a = 0; a < VIM; a++) {
      sumdelx = 0.0;
      for (i = 0; i < n_species - 1; i++) {
        sumdelx += grad_x[i][a];
      }
      grad_x[n_species - 1][a] = -sumdelx;
    }

    for (a = 0; a < VIM; a++) {
      ii2[a] = 0.0;
      for (i = 0; i < n; i += VIM) {
        l = i / VIM;
        for (j = VIM; j < n; j += VIM) {
          m = j / VIM;
          ii2[a] += z[l] * A_inv[i][j] * grad_x[m][a];
        }
      }
      ii2[a] *= F; /* KSC: 7-29-04 */
    }

    /* Derivatives of the vector ii2 with respect to mole fraction are computed
       here; ii2 is the first double sum in Equation (2.13) of the thermal battery LDRD report. --
       RSL 4/2/00 */

    if (pd->v[pg->imtrx][MASS_FRACTION]) {
      for (w = 0; w < pd->Num_Species_Eqn; w++) {
        var = MASS_FRACTION;
        var_offset = MAX_VARIABLE_TYPES + w;
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          for (a = 0; a < dim; a++) {
            dii2dC[a][w][j] = 0.0;
            sum1 = 0.0;
            for (l = 0; l < n_species; l++) {
              if (w == 0) {
                sum1 -= z[l] * A_inv[VIM * l][n - VIM];
              } else {
                sum1 += z[l] * (A_inv[VIM * l][VIM * w] - A_inv[VIM * l][n - VIM]);
              }
            }
            sum2 = 0.0;
            for (l = 0; l < n_species; l++) {
              sum3 = 0.0;
              for (kk = 1; kk < n_species; kk++) {
                sum3 += fv->giant_C_matrix[w][j][VIM * l][VIM * kk] * grad_x[kk][a];
              }
              sum2 += z[l] * sum3;
            }
            dii2dC[a][w][j] = bf[var]->grad_phi[j][a] * sum1 + sum2;
            dii2dC[a][w][j] *= F; /* KSC: 7-29-04 */

          } /* end of loop over a */

        } /* end of loop over j */

      } /* end of loop over w */

    } /* end of var = MASS_FRACTION */

  } /* end of cr->MassFluxModel == STEFAN_MAXWELL_CHARGED */
  /*
   *  Fickian diffusion of charged species in dilute electrolyte
   */
  else if (cr->MassFluxModel == FICKIAN_CHARGED) {
    if (z[n_species - 1] != 0.0) {
      for (a = 0; a < VIM; a++) {
        sumdelx = 0.0;
        for (i = 0; i < n_species - 1; i++) {
          sumdelx += fv->grad_c[i][a];
        }
        fv->grad_c[n_species - 1][a] = -sumdelx;
      }
    }

    for (a = 0; a < VIM; a++) {
      ii2[a] = 0.0;
      for (i = 0; i < n_species; i++) {
        ii2[a] -= z[i] * D[i] * fv->grad_c[i][a];
      }
      ii2[a] *= F;
    }

    /* Derivatives of the vector ii2 (actually ii2/F) with respect to mole fraction -- RSL 9/8/00 */

    /*  if (pd->v[pg->imtrx][MASS_FRACTION] )
        {
         for ( w=0; w<pd->Num_Species_Eqn; w++)
            {
             var = MASS_FRACTION;
             var_offset = MAX_VARIABLE_TYPES + w;
             for ( j=0; j<ei[pg->imtrx]->dof[var]; j++)
                {
                 for ( a=0; a<dim; a++)
                    {
                     dii2dC[a][w][j] = -c*(z[w]*D[w] - z[n_species-1]*D[n_species-1])*
                                       bf[var]->grad_phi[j][a] +
                                       ii2[a]*(d_rho->C[w][j]/rho -
                                               bf[var]->phi[j]*(M[w] - M[n_species-1])/M_mix);
                    }
                }
            }
        } */
  }

  else if (cr->MassFluxModel == FICKIAN_CHARGED_X) /* RSL 9/18/00 */
  {
    for (i = 0; i < n_species - 1; i++) {
      for (a = 0; a < VIM; a++) {
        grad_x[i][a] = fv->grad_c[i][a];
      }
    }

    for (a = 0; a < VIM; a++) {
      ii2[a] = 0.0;
      for (i = 0; i < n_species - 1; i++) {
        ii2[a] -= z[i] * D[i] * grad_x[i][a];
      }
      ii2[a] *= c;
    }

    /* Derivatives of the vector ii2 (actually ii2/F) with respect to mole fraction */

    if (pd->v[pg->imtrx][MASS_FRACTION]) {
      for (w = 0; w < pd->Num_Species_Eqn; w++) {
        var = MASS_FRACTION;
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          for (a = 0; a < dim; a++) {
            dii2dC[a][w][j] = -c * z[w] * D[w] * bf[var]->grad_phi[j][a] +
                              ii2[a] * (d_rho->C[w][j] / rho -
                                        bf[var]->phi[j] * (M[w] - M[n_species - 1]) / M_mix);
          }
        }
      }
    }
  } /* end of cr->MassFluxModel == FICKIAN_CHARGED_X */

  /*
   * Residuals_________________________________________________________________
   */

  if (af->Assemble_Residual) {
    eqn = R_POTENTIAL;
    peqn = upd->ep[pg->imtrx][eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
      phi_i = bf[eqn]->phi[i];
      for (p = 0; p < dim; p++) {
        grad_phi_i[p] = bf[eqn]->grad_phi[i][p];
      }

      diffusion = 0.0;
      if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {
        for (p = 0; p < dim; p++) {
          diffusion += grad_phi_i[p] * ii[p];
        }
        /* add the ii2 contribution in the Poisson equation
         * describing electrolyte potential, KSC: 8/30/00
         */
        if (mp->VoltageFormulation == V_CONDUCTIVITY &&
            (cr->MassFluxModel == STEFAN_MAXWELL_CHARGED ||
             cr->MassFluxModel == STEFAN_MAXWELL_VOLUME || cr->MassFluxModel == FICKIAN_CHARGED ||
             cr->MassFluxModel == FICKIAN_CHARGED_X)) {
          for (p = 0; p < dim; p++) {
            diffusion += grad_phi_i[p] * ii2[p];
          }
        }
        diffusion *= det_J * wt * h3 * pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
      }

      source = 0.0;
      if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
        source = phi_i * h;
        source *= det_J * wt * h3 * pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
      }

      lec->R[LEC_R_INDEX(peqn, i)] += diffusion + source;

    } /* end of loop over i */
  }   /* end of Assemble_Residual */

  /*
   * Jacobian terms...
   */

  if (af->Assemble_Jacobian) {
    eqn = R_POTENTIAL;
    peqn = upd->ep[pg->imtrx][eqn];

    for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {

      phi_i = bf[eqn]->phi[i];

      /*
       * Set up some preliminaries that are needed for the (a,i)
       * equation for bunches of (b,j) column variables...
       */

      for (p = 0; p < dim; p++) {
        grad_phi_i[p] = bf[eqn]->grad_phi[i][p];
      }

      /*
       * J_p_V
       */
      var = VOLTAGE;

      if (pd->v[pg->imtrx][var]) {
        pvar = upd->vp[pg->imtrx][var];
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          for (p = 0; p < dim; p++) {
            grad_phi_j[p] = bf[var]->grad_phi[j][p];
          }

          diffusion = 0.0;
          if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {
            for (p = 0; p < dim; p++) {
              d_i_V_j[p] = -k * grad_phi_j[p] - dkdV[j] * grad_V[p];
            }

            for (p = 0; p < dim; p++) {
              diffusion += grad_phi_i[p] * d_i_V_j[p];
            }
            diffusion *= det_J * wt;
            diffusion *= h3;
            diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
          }

          source = 0.0;
          if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
            source += phi_i * dhdV[j] * det_J * wt;
            source *= h3;
            source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
          }

          lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += diffusion + source;
        }
      }

      /*
       * J_p_T
       */
      var = TEMPERATURE;
      if (pd->v[pg->imtrx][var]) {
        pvar = upd->vp[pg->imtrx][var];
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          for (p = 0; p < dim; p++) {
            grad_phi_j[p] = bf[var]->grad_phi[j][p];
          }

          diffusion = 0.0;
          if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {
            for (p = 0; p < dim; p++) {
              diffusion -= grad_phi_i[p] * dkdT[j] * grad_V[p];
            }
            diffusion *= det_J * wt;
            diffusion *= h3;
            diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
          }

          source = 0.0;
          if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
            source += phi_i * dhdT[j] * det_J * wt;
            source *= h3;
            source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
          }

          lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += diffusion + source;
        }
      }

      /*
       * J_p_v
       */
      for (b = 0; b < dim; b++) {
        var = VELOCITY1 + b;
        if (pd->v[pg->imtrx][var]) {
          pvar = upd->vp[pg->imtrx][var];
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            source = 0.0;
            if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
              source += phi_i * dhdv[b][j] * det_J * wt;
              source *= h3;
              source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
            }

            lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += source;
          }
        }
      }

      /*
       * J_p_d
       */
      for (b = 0; b < dim; b++) {
        var = MESH_DISPLACEMENT1 + b;
        if (pd->v[pg->imtrx][var]) {
          pvar = upd->vp[pg->imtrx][var];
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            dh3dmesh_bj = fv->dh3dq[b] * bf[var]->phi[j];

            d_det_J_dmeshbj = bf[eqn]->d_det_J_dm[b][j];

            diffusion = 0.0;

            if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {
              /*
               * Three parts:
               * 	diff_a = Int(...d(grad_phi_i)/dmesh.q h3 |Jv|)
               *	diff_b = Int(...grad_phi_i.d(q)/dmesh h3 |Jv|)
               *	diff_c = Int(...grad_phi_i.q h3 d(|Jv|)/dmesh)
               *	diff_d = Int(...grad_phi_i.q dh3/dmesh |Jv|  )
               */

              diff_a = 0.;
              for (p = 0; p < dim; p++) {
                dgrad_phi_i_dmesh[p] = bf[eqn]->d_grad_phi_dmesh[i][p][b][j];

                diff_a += dgrad_phi_i_dmesh[p] * ii[p];

                if (mp->VoltageFormulation == V_CONDUCTIVITY &&
                    (cr->MassFluxModel ==
                         FICKIAN_CHARGED || /* add ii2 contribution, KSC: 11/10/00 */
                     cr->MassFluxModel == FICKIAN_CHARGED_X)) /*  RSL 6/24/02  */
                {
                  diff_a += dgrad_phi_i_dmesh[p] * ii2[p];
                }
              }
              diff_a *= det_J * h3 * wt;

              diff_b = 0.0;
              for (p = 0; p < dim; p++) {
                didmesh[p] = -k * fv->d_grad_V_dmesh[p][b][j] - dkdX[b][j] * grad_V[p];

                if (mp->VoltageFormulation == V_CONDUCTIVITY &&
                    cr->MassFluxModel == FICKIAN_CHARGED) /* add ii2 contribution, KSC: 11/10/00 */
                {
                  for (w = 0; w < pd->Num_Species_Eqn; w++) {
                    didmesh[p] -= F * z[w] * D[w] * fv->d_grad_c_dmesh[p][w][b][j];
                  }
                }

                if (mp->VoltageFormulation == V_CONDUCTIVITY &&
                    cr->MassFluxModel == FICKIAN_CHARGED_X) /*  RSL 6/24/02  */
                {
                  for (w = 0; w < n_species - 1; w++) {
                    didmesh[p] -= c * z[w] * D[w] * fv->d_grad_c_dmesh[p][w][b][j];
                  }
                }
              }
              for (p = 0; p < dim; p++) {
                diff_b += grad_phi_i[p] * didmesh[p];
              }
              diff_b *= det_J * h3 * wt;

              diff_c = 0.;
              for (p = 0; p < dim; p++) {
                diff_c += grad_phi_i[p] * ii[p];

                if (mp->VoltageFormulation == V_CONDUCTIVITY &&
                    (cr->MassFluxModel ==
                         FICKIAN_CHARGED || /* add ii2 contribution, KSC: 11/10/00 */
                     cr->MassFluxModel == FICKIAN_CHARGED_X)) /*  RSL 6/24/02  */
                {
                  diff_c += grad_phi_i[p] * ii2[p];
                }
              }
              diff_c *= d_det_J_dmeshbj * h3 * wt;

              diff_d = 0.0;
              for (p = 0; p < dim; p++) {
                diff_d += grad_phi_i[p] * ii[p];
                if (mp->VoltageFormulation == V_CONDUCTIVITY &&
                    (cr->MassFluxModel ==
                         FICKIAN_CHARGED || /* add ii2 contribution, KSC: 11/10/00 */
                     cr->MassFluxModel == FICKIAN_CHARGED_X)) /*  RSL 6/24/02  */
                {
                  diff_d += grad_phi_i[p] * ii2[p];
                }
              }
              diff_d *= det_J * dh3dmesh_bj * wt;

              diffusion = diff_a + diff_b + diff_c + diff_d;

              diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
            }

            source = 0.0;
            if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
              source =
                  phi_i *
                  (h * d_det_J_dmeshbj * h3 + h * det_J * dh3dmesh_bj + dhdX[b][j] * det_J * h3) *
                  wt;

              source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
            }

            lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += diffusion + source;
          }
        }
      }

      /*
       * J_p_c
       */
      var = MASS_FRACTION;
      if (pd->e[pg->imtrx][eqn] && pd->v[pg->imtrx][var]) {
        for (w = 0; w < pd->Num_Species_Eqn; w++) {
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            for (p = 0; p < dim; p++) {
              grad_phi_j[p] = bf[var]->grad_phi[j][p];
            }

            diffusion = 0.0;
            if (pd->e[pg->imtrx][eqn] & T_DIFFUSION) {
              for (p = 0; p < dim; p++) {
                diffusion -= grad_phi_i[p] * dkdC[w][j] * grad_V[p];

                if (mp->VoltageFormulation == V_CONDUCTIVITY) {
                  if (cr->MassFluxModel ==
                          STEFAN_MAXWELL_CHARGED || /* add ii2 contribution for
     Stefan-Maxwell diffusion of charged species in concentrated solutions, KSC: 9/10/00 */
                      cr->MassFluxModel == STEFAN_MAXWELL_VOLUME ||
                      cr->MassFluxModel == FICKIAN_CHARGED_X) /* RSL 9/18
/00 */
                  {
                    diffusion += grad_phi_i[p] * dii2dC[p][w][j];
                  } else if (cr->MassFluxModel ==
                             FICKIAN_CHARGED) /* add ii2 contribution for
transport of dilute charged species in electrolyte, KSC: 9/10/00; modified: 10/4/00 */
                  {
                    diffusion -= grad_phi_i[p] * F * z[w] * D[w] * grad_phi_j[p];
                  }
                } /* end of if VoltageFormulation  */
              }

              diffusion *= det_J * wt;
              diffusion *= h3;
              diffusion *= pd->etm[pg->imtrx][eqn][(LOG2_DIFFUSION)];
            }

            source = 0.0;
            if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
              source += phi_i * dhdC[w][j];
              source *= det_J * wt;
              source *= h3;
              source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
            }

            lec->J[LEC_J_INDEX(peqn, MAX_PROB_VAR + w, i, j)] += diffusion + source;
          } /* end of loop over j */
        }   /* end of loop over w */
      }     /* end of var = MASS_FRACTION */
    }       /* end of loop over i */
  }         /* end of Assemble_Jacobian */

  return (status);
} /* end of assemble_potential */

/****************************************************************************/
void surface_charge_surf(
    double *func,
    double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
    const double sigma) /* value of surface charge, 0 indicates electroneutrality */

/******************************************************************************
 * A routine for imposing the sum of (c^i z^i) = sigma condition at the
 * specified boundary for the specified species so that positively or
 * negatively charged or electrically neutral surfaces can be described.
 * This routine is the generalization of the original y2_electroneutrality_surf
 * routine (written by KSC on 8/31/00) for imposing the electroneutrality
 * constraint on the second species of a three-species electrolyte system.
 * Here, for ease of programming and also for solution accuracy, the nth
 * species is taken to be the electrically neutral solvent (e.g. water
 * in the case of electrodeposition with aqueous electrolyte bath).
 *
 * Author: Ken S. Chen (12/21/01)
 *
 *****************************************************************************/
{
  int var, i, j;
  double sum;
  double z[MAX_CONC];

  for (i = 0; i < pd->Num_Species; i++) {
    z[i] = mp->charge_number[i];
  }

  sum = 0.0;
  for (i = 0; i < pd->Num_Species - 1; i++) {
    sum += fv->c[i] * z[i];
  }

  *func = sum - sigma;

  var = MASS_FRACTION;
  if (pd->v[pg->imtrx][var]) {
    for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
      for (i = 0; i < pd->Num_Species; i++) {
        d_func[0][MAX_VARIABLE_TYPES + i][j] = z[i] * bf[var]->phi[j];
      }
    }
  }

  return;
} /* End of surface_charge_surf                                */
/*****************************************************************************/

/*****************************************************************************/
void current_BV_surf(double func[DIM],
                     double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                     int wspec,     /* species number of this boundary condition */
                     double nu,     /* stoichiometric coefficient                */
                     double k,      /* kinetic rate constant                     */
                     double beta,   /* reaction order                            */
                     double alphaa, /* anodic direction transfer coefficient     */
                     double alphac, /* cathodic direction transfer coefficient   */
                     double V,      /* electrode potential                       */
                     double U0,     /* theoretical open-circuit potential        */
                     double dt,     /* current value of the time step            */
                     double tt)     /* parameter to vary time integration        */

/******************************************************************************
 *
 *  A function that calculates the electrical current density using Faraday's law
 *  with the species molar rate given by Butler-Volmer kinetics.
 *
 *  Ken S. Chen (9/2001)
 *
 *  Revised by KSC on 12/17/2001
 *  Revised by KSC on 1/7/2002
 *
 ******************************************************************************/
{
  double T = -1.0;                           /* electrolyte solution temperature            */
  double z;                                  /* charge number of species wspec              */
  double r;                                  /* consum. or gen. molar rate of species wspec */
  double i;                                  /* current density                             */
  double d_i[MAX_VARIABLE_TYPES + MAX_CONC]; /* Jacobian for cuurent density  */

  double vconv[MAX_PDIM];     /*Calculated convection velocity */
  double vconv_old[MAX_PDIM]; /*Calculated convection velocity at previous time*/
  CONVECTION_VELOCITY_DEPENDENCE_STRUCT d_vconv_struct;
  CONVECTION_VELOCITY_DEPENDENCE_STRUCT *d_vconv = &d_vconv_struct;

  int err = 0;

  const double R = 8.314;   /* Universal gas constant in units of J/mol K  */
  const double F = 96487.0; /* Faraday's constant in units of C/equiv.     */
  double FRT;               /* product of F/R/T                            */
  double PHI;               /* electrolyte potential on electrode surface  */
  double c;                 /* species concentration on electrode surface  */
  double conc, conc1, grpa, grpc;

  int w, j, j_id, w1, dim, kdir, var, jvar;
  double phi_j;

  z = mp->charge_number[wspec];
  if (mp->SolutionTemperatureModel == CONSTANT) {
    T = mp->solution_temperature;
  } else {
    GOMA_EH(GOMA_ERROR, "Solution-temperature model not yet implemented");
  }
  if (pd->e[pg->imtrx][R_ENERGY]) /* if energy equation is active, re-set electrolyte temperature */
  {
    T = fv->T;
  }
  /* set solution temperature to 298 K if it is zero - safety feature */
  if (T == 0.0) {
    T = 298.0;
    fprintf(stderr, "Warning!: a default electrolyte temperature of 298 K is being used!");
  }

  if (nAC) {                     /* if augmenting condition card is active, then */
    V = mp->electrode_potential; /* get electrode potential from AC calculation */
  }

  c = fv->c[wspec];
  PHI = fv->V;

  FRT = F / R / T;
  conc = pow(c, beta);
  conc1 = beta * pow(c, beta - 1.0);
  grpa = alphaa * FRT * (V - PHI - U0);
  grpc = alphac * FRT * (V - PHI - U0);

  dim = pd->Num_Dim;

  /* get deformation gradient, etc. if lagrangian mesh with inertia */
  if (pd->MeshMotion == LAGRANGIAN && pd->MeshInertia == 1) {
    err = belly_flop(elc->lame_mu);
    GOMA_EH(err, "error in belly flop");
    if (neg_elem_volume)
      return;
  }

  /* get the convection velocity (it's different for arbitrary and
     lagrangian meshes) */
  if (mp->PorousMediaType == CONTINUOUS) {
    err = get_convection_velocity(vconv, vconv_old, d_vconv, dt, tt);
  }
  GOMA_EH(err, "Error in calculating effective convection velocity");

  /* compute molar rate using Butler-Volmer kinetics */
  r = nu * k * conc * (exp(grpa) - exp(-grpc));

  /* compute current density using Faraday's law */
  i = z * F * r;
  *func = -i;

  /* Calculate the residual contribution from convective flux       */
  for (kdir = 0; kdir < dim; kdir++) {
    *func += c * vconv[kdir] * fv->snormal[kdir];
  }

  /* sum the contributions to the global stiffness matrix  for Species*/

  for (w = 0; w < pd->Num_Species_Eqn; w++) /* no dependence on other species */
  {
    d_i[MAX_VARIABLE_TYPES + w] = 0.;
  }
  d_i[MAX_VARIABLE_TYPES + wspec] = z * F * nu * k * conc1 * (exp(grpa) - exp(-grpc));
  d_i[TEMPERATURE] = z * F * nu * k * conc * (FRT / T) * (-V + PHI + U0) *
                     (alphaa * exp(grpa) + alphac * exp(-grpc));
  d_i[VOLTAGE] = -z * F * nu * k * conc * FRT * (alphaa * exp(grpa) + alphac * exp(-grpc));

  /* J_s_c --- sensitivity wrt species concentrations */
  var = MASS_FRACTION;
  for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++) {
    phi_j = bf[var]->phi[j_id];

    for (w1 = 0; w1 < pd->Num_Species_Eqn; w1++) {
      d_func[0][MAX_VARIABLE_TYPES + w1][j_id] -= d_i[MAX_VARIABLE_TYPES + w1] * phi_j;
    }

    for (kdir = 0; kdir < dim; kdir++) {
      d_func[0][MAX_VARIABLE_TYPES + wspec][j_id] += phi_j * vconv[kdir] * fv->snormal[kdir];

      for (w1 = 0; w1 < pd->Num_Species_Eqn; w1++) {
        d_func[0][MAX_VARIABLE_TYPES + w1][j_id] +=
            c * d_vconv->C[kdir][w1][j_id] * fv->snormal[kdir];
      }
    }
  }

  /* J_s_T --- sensitivity wrt electrolyte solution temperature */
  var = TEMPERATURE;
  if (pd->v[pg->imtrx][var]) {
    for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
      j_id = j;
      phi_j = bf[var]->phi[j_id];

      d_func[0][var][j_id] -= d_i[TEMPERATURE] * phi_j;

      for (kdir = 0; kdir < dim; kdir++) {
        d_func[0][var][j_id] += c * d_vconv->T[kdir][j_id] * fv->snormal[kdir];
      }
    }
  }

  /* J_s_V --- sensitivity wrt electrolyte potential */
  var = VOLTAGE;
  if (pd->v[pg->imtrx][var]) {
    for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
      j_id = j;
      phi_j = bf[var]->phi[j_id];

      d_func[0][var][j_id] -= d_i[VOLTAGE] * phi_j;

      for (kdir = 0; kdir < dim; kdir++) {
        d_func[0][var][j_id] += c * d_vconv->T[kdir][j_id] * fv->snormal[kdir];
      }
    }
  }

  /* J_s_d --- sensitivity wrt to nodal mesh displacements */
  for (jvar = 0; jvar < dim; jvar++) {
    var = MESH_DISPLACEMENT1 + jvar;
    if (pd->v[pg->imtrx][var]) {
      for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++) {
        phi_j = bf[var]->phi[j_id];
        for (kdir = 0; kdir < dim; kdir++) {
          d_func[0][var][j_id] += c * vconv[kdir] * fv->dsnormal_dx[kdir][jvar][j_id] +
                                  c * d_vconv->X[kdir][jvar][j_id] * fv->snormal[kdir];
        }
      }
    }
  }

  for (jvar = 0; jvar < dim; jvar++) {
    var = VELOCITY1 + jvar;
    if (pd->v[pg->imtrx][var]) {
      for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++) {
        phi_j = bf[var]->phi[j_id];
        d_func[0][var][j_id] += c * d_vconv->v[jvar][jvar][j_id] * fv->snormal[jvar];
      }
    }
  }

  return;
} /* END of routine current_BV_surf */
/*****************************************************************************/

void current_ORR_surf(double func[],
                      double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                      int wspec,     /* species number of this boundary condition */
                      double ai0,    /* product of interfacial area by            */
                                     /* exchange current density (A/cm^3)         */
                      double H,      /* thickness of catalyst layer (cm)          */
                      double cref,   /* species ref. concentration (moles/cm^3)   */
                      double alphac, /* cathodic direction transfer coefficient   */
                      double T,      /* cell temperature (K)                      */
                      double V,      /* cell voltage (K)                          */
                      double U0,     /* Open-circuit potential for HOR (V)        */
                      double beta,   /* reaction order                            */
                      double dt,     /* current value of the time step            */
                      double tt)     /* parameter to vary time integration        */

/******************************************************************************
 *
 *  A function that calculates the electrical current density using Faraday's law
 *  with the species molar rate given by Tafel kinetics in an electrochemical
 *  reaction such as the oxygen reduction reaction in PEM fuel cells.
 *
 *  Ken S. Chen (1/2006)
 *  Modified: 5/22/2006 by KSC.
 *
 ******************************************************************************/
{
  double i;                 /* current density in units of A/cm^2          */
  double didc, didT, didV;  /* sensitivities                               */
  double c;                 /* species concentration on electrode surface  */
  double PHI;               /* electrical potential in electrolyte phase   */
  const double F = 96487.0; /* Faraday's constant in units of C/equiv.     */
  const double R = 8.314;   /* Universal gas constant in units of J/mole K */
  double FRT;               /* F/R/T                                       */
  double grp, cratio;
  int j, j_id, w1, var;
  double phi_j;

  /* compute current density */
  FRT = F / R / T;
  PHI = fv->V;
  c = fv->c[wspec];
  if (c < 0.0)
    c = 1.0e-10;
  if (c == 0)
    c = cref;
  cratio = pow(c / cref, beta);
  grp = alphac * FRT * (V - PHI - U0);

  i = ai0 * H * cratio * exp(-grp);

  didc = i * (beta / c);
  didT = i * (grp / T);
  didV = i * alphac * FRT;

  /* compute the current desnity residual */
  *func = -i;

  /* J_s_c --- sensitivity wrt species concentrations */
  var = MASS_FRACTION;
  for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++) {
    phi_j = bf[var]->phi[j_id];

    for (w1 = 0; w1 < pd->Num_Species_Eqn; w1++) {
      d_func[0][MAX_VARIABLE_TYPES + w1][j_id] = 0.0;
    }

    d_func[0][MAX_VARIABLE_TYPES][j_id] = -didc * phi_j;
  }

  /* J_s_T --- sensitivity wrt electrolyte solution temperature */
  var = TEMPERATURE;
  if (pd->v[pg->imtrx][var]) {
    for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
      j_id = j;
      phi_j = bf[var]->phi[j_id];

      d_func[0][var][j_id] = -didT * phi_j;
    }
  }

  /* J_s_V --- sensitivity wrt electrolyte potential */
  var = VOLTAGE;
  if (pd->v[pg->imtrx][var]) {
    for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
      j_id = j;
      phi_j = bf[var]->phi[j_id];

      d_func[0][var][j_id] = -didV * phi_j;
    }
  }

  return;
} /* END of routine current_ORR_surf */
/*****************************************************************************/

void current_HOR_surf(double func[],
                      double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                      int wspec,     /* species number of this boundary condition */
                      double ai0,    /* product of interfacial area by            */
                                     /* exchange current density (A/cm^3)         */
                      double H,      /* thickness of anode catalyst layer (cm)    */
                      double cref,   /* species ref. concentration (moles/cm^3)   */
                      double alphaa, /* anodic direction transfer coefficient     */
                      double alphac, /* cathodic direction transfer coefficient   */
                      double T,      /* cell temperature (K)                      */
                      double U0,     /* open-circuit potential (V)                */
                      double beta,   /* reaction order                            */
                      double V,      /* electrode potential (V)                   */
                      double dt,     /* current value of the time step            */
                      double tt)     /* parameter to vary time integration        */

/******************************************************************************
 *
 *  A function that calculates the electrical current density using Faraday's law
 *  with the species molar rate given by the linearized kinetics in the hydrogen
 *  oxidation reaction
 *
 *  Ken S. Chen (1/2006)
 *  Modified: 5/22/2006 by KSC.
 *
 ******************************************************************************/
{
  double i;                 /* current density in units of A/cm^2          */
  double didc, didT, didV;  /* sensitivities                               */
  double c;                 /* species concentration on electrode surface  */
  double PHI;               /* electrical potential in electrolyte phase   */
  const double F = 96487.0; /* Faraday's constant in units of C/equiv.     */
  const double R = 8.314;   /* Universal gas constant in units of J/mole K */
  double FRT;               /* F/R/T                                       */
  double cratio;
  int j, j_id, w1, var;
  double phi_j;

  /* compute current density */
  FRT = F / R / T;
  PHI = fv->V;
  c = fv->c[wspec];
  if (c < 0.0)
    c = 1.0e-10;
  if (c == 0)
    c = cref;
  cratio = pow(c / cref, beta);

  i = ai0 * H * cratio * (alphaa + alphac) * FRT * (V - PHI - U0);

  didc = i * (beta / c);
  didT = -i / T;
  didV = -ai0 * H * cratio * (alphaa + alphac) * FRT;

  /* compute the current desnity residual */
  *func = -i;

  /* J_s_c --- sensitivity wrt species concentrations */
  var = MASS_FRACTION;
  for (j_id = 0; j_id < ei[pg->imtrx]->dof[var]; j_id++) {
    phi_j = bf[var]->phi[j_id];

    for (w1 = 0; w1 < pd->Num_Species_Eqn; w1++) {
      d_func[0][MAX_VARIABLE_TYPES + w1][j_id] = 0.0;
    }

    d_func[0][MAX_VARIABLE_TYPES][j_id] = -didc * phi_j;
  }

  /* J_s_T --- sensitivity wrt electrolyte solution temperature */
  var = TEMPERATURE;
  if (pd->v[pg->imtrx][var]) {
    for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
      j_id = j;
      phi_j = bf[var]->phi[j_id];

      d_func[0][var][j_id] = -didT * phi_j;
    }
  }

  /* J_s_V --- sensitivity wrt electrolyte potential */
  var = VOLTAGE;
  if (pd->v[pg->imtrx][var]) {
    for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
      j_id = j;
      phi_j = bf[var]->phi[j_id];

      d_func[0][var][j_id] = -didV * phi_j;
    }
  }

  return;
} /* END of routine current_HOR_surf */
/*****************************************************************************/

int assemble_Enorm(void) {
  int dim, i, j, k;
  dbl wt, detJ, h3, wt_func;
  dbl Enorm, gradEcomp;
  dbl d_Enorm_d_V[MDE];
  dbl advection, source;

  dim = pd->Num_Dim;
  if (!pd->e[pg->imtrx][R_ENORM])
    return 0;

  if (!pd->v[pg->imtrx][VOLTAGE])
    GOMA_EH(GOMA_ERROR, "Must have VOLTAGE equation active to use ENORM.");

  if (pd->v[pg->imtrx][MESH_DISPLACEMENT1])
    GOMA_EH(GOMA_ERROR, "Sorry, assemble_Enorm has not been rigged for ALE.");

  wt = fv->wt;
  h3 = fv->h3;
  detJ = bf[R_ENORM]->detJ;

  /* Compute the value of |E|, where E = grad(V). */
  Enorm = 0.0;
  for (j = 0; j < dim; j++) {
    gradEcomp = fv->grad_V[j];
    Enorm += gradEcomp * gradEcomp;
  }
  Enorm = MAX(0.0, Enorm);
  Enorm = sqrt(Enorm);

  memset(d_Enorm_d_V, 0, MDE * sizeof(dbl));
  /* Get d(Enormsq)/d(V_k). */
  /* Kludge to avoid singularity.  Should be updated after 1st iteration with real values. */
  if (Enorm > 0.0)
    for (k = 0; k < ei[pg->imtrx]->dof[VOLTAGE]; k++) {
      for (j = 0; j < dim; j++)
        d_Enorm_d_V[k] += fv->grad_Enorm[j] * bf[VOLTAGE]->grad_phi[k][j];
      d_Enorm_d_V[k] /= Enorm;
    }

  if (af->Assemble_Residual) {
    for (i = 0; i < ei[pg->imtrx]->dof[R_ENORM]; i++) {
      wt_func = bf[R_ENORM]->phi[i];

      /* Might as well call the one term, |E|^2, the advection term... */
      advection = 0.0;
      if (pd->e[pg->imtrx][R_ENORM] & T_ADVECTION) {
        advection = -Enorm;
        advection *= wt_func * detJ * wt * h3;
        advection *= pd->etm[pg->imtrx][R_ENORM][(LOG2_ADVECTION)];
      }

      /* And of course, the glorious LHS, "ENORM", shall just be the source term... */
      source = 0.0;
      if (pd->e[pg->imtrx][ENORM] & T_SOURCE) {
        source = fv->Enorm;
        source *= wt_func * detJ * wt * h3;
        source *= pd->etm[pg->imtrx][R_ENORM][(LOG2_SOURCE)];
      }

      lec->R[LEC_R_INDEX(upd->ep[pg->imtrx][R_ENORM], i)] += advection + source;
    } /* for i = 0 ... dof[R_ENORM] */
  }   /* if(Assemble_Residual) */

  if (af->Assemble_Jacobian) {
    for (i = 0; i < ei[pg->imtrx]->dof[R_ENORM]; i++) {
      wt_func = bf[R_ENORM]->phi[i];

      /* J_Enorm_V */
      for (k = 0; k < ei[pg->imtrx]->dof[VOLTAGE]; k++) {
        /* Only the advection term for _V. */
        advection = 0.0;
        if (pd->e[pg->imtrx][R_ENORM] & T_ADVECTION) {
          advection = d_Enorm_d_V[k];
          advection *= wt_func * detJ * wt * h3;
          advection *= pd->etm[pg->imtrx][R_ENORM][(LOG2_ADVECTION)];
        }
        lec->J[LEC_J_INDEX(upd->ep[pg->imtrx][R_ENORM], upd->vp[pg->imtrx][VOLTAGE], i, k)] +=
            advection;
      } /* for k = 0 ... dof[VOLTAGE] */

      /* J_Enorm_Enorm */
      for (k = 0; k < ei[pg->imtrx]->dof[ENORM]; k++) {
        /* Only the source term for _Enorm. */
        source = 0.0;
        if (pd->e[pg->imtrx][R_ENORM] & T_SOURCE) {
          source = bf[ENORM]->phi[k];
          source *= wt_func * detJ * wt * h3;
          source *= pd->etm[pg->imtrx][R_ENORM][(LOG2_SOURCE)];
        }

        lec->J[LEC_J_INDEX(upd->ep[pg->imtrx][R_ENORM], upd->vp[pg->imtrx][ENORM], i, k)] += source;
      } /* for k = 0 ... dof[ENORM] */
    }   /* for i = 0 ... dof[R_ENORM] */
  }     /* if (Assemble_Jacobian) */
  return 0;
}

int assemble_electric_field(void) /* Least square equation Efield = grad (voltage) */
{
  int dim;
  int p, a;

  int eqn, var;
  int peqn, pvar;
  int i, j;
  int status;

  dbl h3;          /* Volume element (scale factors). */
  dbl dh3dmesh_pj; /* Sensitivity to (p,j) mesh dof. */

  dbl grad_V[DIM]; /* gradient of voltage */
  dbl Efield[DIM]; /* Electric field vector */

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

  status = 0;

  /*
   * Unpack variables from structures for local convenience...
   */

  dim = pd->Num_Dim;

  eqn = R_EFIELD1;

  /*
   * Bail out fast if there's nothing to do...
   */

  if (!pd->e[pg->imtrx][eqn]) {
    return (status);
  }

  wt = fv->wt;

  det_J = bf[eqn]->detJ; /* Really, ought to be mesh eqn. */

  h3 = fv->h3; /* Differential volume element (scales). */

  /*
   * Field variables...
   */

  /* Gradient of voltage equation */

  for (a = 0; a < VIM; a++) {
    grad_V[a] = fv->grad_V[a];
    Efield[a] = fv->E_field[a];
  }

  /*
   * Residuals_________________________________________________________________
   */

  if (af->Assemble_Residual) {
    /*
     * Assemble each component "ab" of the velocity gradient equation...
     */
    for (a = 0; a < VIM; a++) {
      eqn = R_EFIELD1 + a;
      /*
       * In the element, there will be contributions to this many equations
       * based on the number of degrees of freedom...
       */

      for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {

        wt_func = bf[eqn]->phi[i]; /* add Petrov-Galerkin terms as necessary */

        advection = 0.;

        if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
          advection -= grad_V[a];
          advection *= wt_func * det_J * wt * h3;
          advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
        }

        /*
         * Source term...
         */

        source = 0;

        if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
          source += Efield[a];
          source *= wt_func * det_J * h3 * wt;
          source *= pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
        }

        lec->R[LEC_R_INDEX(upd->ep[pg->imtrx][eqn], i)] += advection + source;
      }
    }
  }

  /*
   * Jacobian terms...
   */

  if (af->Assemble_Jacobian) {
    for (a = 0; a < VIM; a++) {
      eqn = R_EFIELD1 + a;
      peqn = upd->ep[pg->imtrx][eqn];

      for (i = 0; i < ei[pg->imtrx]->dof[eqn]; i++) {
        wt_func = bf[eqn]->phi[i]; /* add Petrov-Galerkin terms as necessary */

        /*
         * J_Efield_V
         */
        var = VOLTAGE;
        if (pd->v[pg->imtrx][var]) {
          pvar = upd->vp[pg->imtrx][var];
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];

            advection = 0.;

            if (pd->e[pg->imtrx][eqn] & T_ADVECTION) {
              advection -= bf[var]->grad_phi[j][a];
              advection *= wt_func * det_J * wt * h3;
              advection *= pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
            }

            source = 0.;

            lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += advection + source;
          }
        }
        /*
         * J_Efield_d
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

                advection_a = -grad_V[a];

                advection_a *= (d_det_J_dmesh_pj * h3 + det_J * dh3dmesh_pj);

                advection_b = -fv->d_grad_V_dmesh[a][p][j];

                advection_b *= det_J * h3;

                advection = advection_a + advection_b;

                advection *= wt_func * wt * pd->etm[pg->imtrx][eqn][(LOG2_ADVECTION)];
              }

              source = 0.;

              if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
                source += Efield[a];

                source *= d_det_J_dmesh_pj * h3 + det_J * dh3dmesh_pj;

                source *= wt_func * wt * pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
              }

              lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += advection + source;
            }
          }
        }

        /*
         * J_Efield_Efield
         */

        var = EFIELD1 + a;

        if (pd->v[pg->imtrx][var]) {
          pvar = upd->vp[pg->imtrx][var];
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            phi_j = bf[var]->phi[j];

            source = 0.;

            if (pd->e[pg->imtrx][eqn] & T_SOURCE) {
              source = phi_j * det_J * h3 * wt_func * wt * pd->etm[pg->imtrx][eqn][(LOG2_SOURCE)];
            }

            lec->J[LEC_J_INDEX(peqn, pvar, i, j)] += source;
          }
        }
      }
    }
  }

  return (status);
} /* end of assemble_electric_field */

/*******************************************************************************/
void apply_potential_grad_bc(double func[],
                             double d_func[][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                             struct Boundary_Condition *BC_Type,
                             const double time,
                             const double dt)

/*
  Apply "CURRENT" bc as a strong integrated condition

  Author : Robert B. Secor, 3M
  Date   : March 9, 2004

Parameters:
  func = pointer to double that carries back the residual value
  d_func = array of sensitivities of residual wrt to all variables.
  BC_Type = pointer to Boundary Condition structure,
            i.e. &(BC_Types[bc_input_id]
*/
{
  int a, b, j, p, q, var;
  dbl perm;                  /* electrical permittivity */
  dbl d_p_dT[MDE];           /* Temperature derivative of permittivity. */
  dbl d_p_dV[MDE];           /* Potential derivative of permittivity. */
  dbl d_p_dC[MAX_CONC][MDE]; /* Concentration derivative of permittivity. */
  dbl d_p_dX[DIM][MDE];      /* Spatial derivatives of permittivity. */
  dbl efield[MAX_PDIM];      /* electric field */
  double grad_phi_j;

  if (af->Assemble_LSA_Mass_Matrix)
    return;

  *func = BC_Type->BC_Data_Float[0];

  /* Get electrical conductivity
      if (mp->VoltageFormulation == V_PERMITTIVITY)
         {
           perm = mp->permittivity;
         }
       else
         {
           perm = mp->electrical_conductivity;
         }

  memset( d_p_dV,0, sizeof(double)*MDE);
  memset( d_p_dT,0, sizeof(double)*MDE);
  memset( d_p_dX,0, sizeof(double)*MDE*DIM);
  memset( d_p_dC,0, sizeof(double)*MAX_CONC*MDE);
*/
  perm = electrical_conductivity(d_p_dT, d_p_dV, d_p_dC, d_p_dX, time, dt);

  for (a = 0; a < VIM; a++) {
    efield[a] = -fv->grad_V[a];
  }

  if (af->Assemble_Jacobian) {

    var = VOLTAGE;
    if (pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        for (p = 0; p < VIM; p++) {
          grad_phi_j = bf[var]->grad_phi[j][p];
          d_func[0][var][j] += fv->snormal[p] * (perm * grad_phi_j - d_p_dV[j] * efield[p]);
        }
      }
    }

    var = TEMPERATURE;
    if (pd->v[pg->imtrx][var]) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        for (p = 0; p < VIM; p++) {
          grad_phi_j = bf[var]->grad_phi[j][p];
          d_func[0][var][j] += fv->snormal[p] * (-d_p_dT[j] * efield[p]);
        }
      }
    }

    if (pd->v[pg->imtrx][MESH_DISPLACEMENT1]) {
      for (b = 0; b < VIM; b++) {
        var = MESH_DISPLACEMENT1 + b;
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          for (q = 0; q < VIM; q++) {
            d_func[0][var][j] +=
                fv->snormal[q] * (perm * fv->d_grad_V_dmesh[q][b][j] - d_p_dX[b][j] * efield[q]) -
                fv->dsnormal_dx[q][b][j] * (perm * efield[q]);
          }
        }
      }
    }
  }

  /* Calculate the residual contribution                                       */

  for (a = 0; a < VIM; a++) {
    *func += -perm * fv->snormal[a] * efield[a];
  }

  return;
}
/** end of apply_potential_grad_bc */

/*******************************************************************************/
double electrical_conductivity(double dkdT[MDE],
                               double dkdV[MDE],
                               double dkdC[MAX_CONC][MDE],
                               double dkdX[DIM][MDE],
                               const double time,
                               const double dt)

{
  int dim, var, var_offset;
  int i, j, a, w;
  dbl k = 1e12; /* electrical conductivity. */
  int n_species = mp->Num_Species;
  dbl sumx, sumk; /* dummy variables */
  dbl A_inv[MAX_CONC * DIM]
           [MAX_CONC * DIM]; /* inverse of Stefan_Maxwell coefficient matrix; KSC: 10/28/98 */
  dbl *D = mp->diffusivity;  /* Fickian diffusion coefficients */
  dbl T = 298.15;            /* Electrolyte solution temperature default = 25 C */
  dbl c;                     /* Total molar concentration of electrolyte solution */
  const double R = 8.314;    /* Universal gas constant in units of J/mole K */
  const double F = 96487.0;  /* Faraday's constant in units of C/equiv. */
  dbl FRT;                   /* product of F/R/T */
  dbl rho;                   /* Mass density of electrolyte solution */
  DENSITY_DEPENDENCE_STRUCT d_rho_struct; /* density dependence */
  DENSITY_DEPENDENCE_STRUCT *d_rho = &d_rho_struct;
  dbl *M = mp->molecular_weight; /* Ptr to species molecular weights */
  dbl M_mix;                     /* Molecular weight of electrolyte solution */
  dbl x[MAX_CONC];               /* mole fraction */

  dbl *z = mp->charge_number; /* ptr to the charge numbers for species */

  int i_elec_cond; /* index of external field elec. conductivity */
  int n, l, m, mn, kk;
  dbl sum1, sum2, sum3;

  dim = pd->Num_Dim;
  /*** Electrical Conductivity ****/
  if (mp->Elec_ConductivityModel == USER) {
    (void)usr_electrical_conductivity(mp->u_electrical_conductivity, time);

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
        var_offset = MAX_VARIABLE_TYPES + w;
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          dkdC[w][j] = mp->d_electrical_conductivity[var_offset] * bf[var]->phi[j];
        }
      }
    }

  }

  else if (mp->Elec_ConductivityModel == ELECTRONEUTRALITY_FICKIAN ||
           mp->Elec_ConductivityModel == ELECTRONEUTRALITY_SM ||
           mp->Elec_ConductivityModel == ELECTRODE_KINETICS) {
    if (cr->MassFluxModel == STEFAN_MAXWELL_CHARGED || cr->MassFluxModel == STEFAN_MAXWELL_VOLUME) {
      n = n_species * VIM;

      if (mp->SolutionTemperatureModel == CONSTANT) {
        T = mp->solution_temperature;
      } else if (mp->SolutionTemperatureModel == THERMAL_BATTERY) {
        electrolyte_temperature(time, dt, 0);
        T = mp->electrolyte_temperature;
      } else {
        GOMA_EH(
            GOMA_ERROR,
            "Solution-temperature model other than THERMAL_BATTERY awaits future implementation");
      }
      /*
       * set the solution temperature to the 298 K if it is zero - safety feature
       */
      if (T == 0.0)
        T = 298.0;

      FRT = F / (R * T);
      /*
       * get density from the density() routine defined in mm_fill_terms.c
       */
      rho = density(d_rho, time);

      /*
       * Calculate the mole fraction of the species, even for the last species
       */
      sumx = 0.0;
      for (i = 0; i < n_species - 1; i++) {
        x[i] = fv->c[i];
        sumx += x[i];
      }
      x[n_species - 1] = 1.0 - sumx;
      /*
       * Calculate the molecular weight of the fluid, M_mix (gm/mole), and
       * the molar concentration of the  fluid, c (mole/cm**3).
       */
      M_mix = 0.0;
      for (i = 0; i < n_species; i++) {
        M_mix += M[i] * x[i];
      }
      c = rho / M_mix;
      /*
       * retrieve the inverse matrix from the fv structure
       */
      for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
          A_inv[i][j] = fv->SM_matrix_inv[i][j];
        }
      }

      sumk = 0.0;
      for (i = 0; i < n; i += VIM) {
        l = i / VIM;
        for (j = VIM; j < n; j += VIM) {
          m = j / VIM;
          /* sumk -= FRT*z[l]*z[m]*A_inv[i][j]*x[m]; */
          sumk -= F * FRT * z[l] * z[m] * A_inv[i][j] * x[m]; /* KSC: 7-29-04 */
        }
      }
      k = sumk; /* This is the electrical conductivity value in the PHI2 equation; KSC: 10/98  */

      var = TEMPERATURE;
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        dkdT[j] = 0.0;
      }

      var = VOLTAGE;
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        dkdV[j] = 0.0;
      }

      if (pd->v[pg->imtrx][MESH_DISPLACEMENT1]) {
        for (a = 0; a < dim; a++) {
          var = MESH_DISPLACEMENT1 + a;
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            dkdX[a][j] = 0.0;
          }
        }
      }

      /* Derivatives of electrical conductivity with respect to mole fraction
         are computed here; kappa is defined by Equation (2.14) of the thermal battery LDRD report.
         -- RSL 3/31/00 */

      if (pd->v[pg->imtrx][MASS_FRACTION]) {
        for (w = 0; w < pd->Num_Species_Eqn; w++) {
          var = MASS_FRACTION;
          var_offset = MAX_VARIABLE_TYPES + w;
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            dkdC[w][j] = 0.0;
            sum1 = 0.0;
            for (l = 0; l < n_species; l++) {
              if (w == 0) {
                sum1 -= z[l] * z[n_species - 1] * A_inv[VIM * l][n - VIM];
              } else {
                sum1 += z[l] * (z[w] * A_inv[VIM * l][VIM * w] -
                                z[n_species - 1] * A_inv[VIM * l][n - VIM]);
              }
            }
            sum2 = 0.0;
            for (l = 0; l < n_species; l++) {
              sum3 = 0.0;
              for (kk = 1; kk < n_species; kk++) {
                sum3 += z[kk] * x[kk] * fv->giant_C_matrix[w][j][VIM * l][VIM * kk];
              }
              sum2 += z[l] * sum3;
            }
            /* dkdC[w][j] = -FRT * (bf[var]->phi[j]*sum1 + sum2); */
            dkdC[w][j] = -F * FRT * (bf[var]->phi[j] * sum1 + sum2); /* KSC: 7-29-04 */

          } /* end of loop over j */

        } /* end of loop over w */

      } /* end of var = MASS_FRACTION */

    } /* end of cr->MassFluxModel == STEFAN_MAXWELL_CHARGED */

    else if (cr->MassFluxModel == FICKIAN_CHARGED) {
      if (mp->SolutionTemperatureModel == CONSTANT) {
        T = mp->solution_temperature;
      } else if (mp->SolutionTemperatureModel == THERMAL_BATTERY) {
        electrolyte_temperature(time, dt, 0);
        T = mp->electrolyte_temperature;
      } else {
        GOMA_EH(
            GOMA_ERROR,
            "Solution-temperature model other than THERMAL_BATTERY awaits future implementation");
      }
      if (T == 0.0)
        T = 298.0; /* set the solution temperature to the 298 K if it is zero - safety feature */

      FRT = F / R / T;

      /*
       * In the next block, assert an electroneutrality condition by
       * setting the concentration of the last species in the
       * mechanism by solving the equation,
       *
       *  sum (z_i * c_i) = 0
       *
       * HKM -> Note, this section of code is misplaced. It should be put in
       *        a general place, mm_fill_terms.c, so that fv->c[] may be
       *        used consistently throughout the fill.
       *        Also, it seems very restrictive, as there are many situations
       *        in solids when electroneutrality shouldn't be imposed.
       */
      if (z[n_species - 1] != 0.0) {
        sumx = 0.0;
        for (i = 0; i < n_species - 1; i++) {
          sumx += (z[i] / z[n_species - 1]) * fv->c[i];
        }
        fv->c[n_species - 1] = -sumx;
      }

      sumk = 0.0;
      for (i = 0; i < n_species; i++) {
        sumk += z[i] * z[i] * D[i] * fv->c[i];
      }
      k = F * FRT * sumk;

      var = TEMPERATURE;
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        dkdT[j] = -bf[var]->phi[j] * k / T;
      }

      var = VOLTAGE;
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        dkdV[j] = 0.0;
      }

      if (pd->v[pg->imtrx][MESH_DISPLACEMENT1]) {
        for (a = 0; a < dim; a++) {
          var = MESH_DISPLACEMENT1 + a;
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            dkdX[a][j] = 0.0;
          }
        }
      }

      if (pd->v[pg->imtrx][MASS_FRACTION]) {
        for (w = 0; w < pd->Num_Species_Eqn; w++) {
          var = MASS_FRACTION;
          var_offset = MAX_VARIABLE_TYPES + w;
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            dkdC[w][j] = F * FRT * z[w] * z[w] * D[w] * bf[var]->phi[j];
          }
        }
      }
    } /* end of if (cr->MassFluxModel == FICKIAN_CHARGED) */

    else if (cr->MassFluxModel == FICKIAN_CHARGED_X) /* RSL 9/18/00 */
    {
      mn = ei[pg->imtrx]->mn;

      if (pd_glob[mn]->e[pg->imtrx][R_ENERGY]) /* if the energy equation is being solved */
      {
        T = fv->T;
      } else /* if the energy equation is NOT being solved */
      {
        if (mp->SolutionTemperatureModel == CONSTANT) {
          T = mp->solution_temperature;
        } else if (mp->SolutionTemperatureModel == THERMAL_BATTERY) {
          electrolyte_temperature(time, dt, 0);
          T = mp->electrolyte_temperature;
        } else {
          GOMA_EH(GOMA_ERROR, "Invalid solution temperature model");
        }
      }

      FRT = F / (R * T);

      rho = density(d_rho, time); /* RSL 6/22/02 */

      sumx = 0.0;
      for (i = 0; i < n_species - 1; i++) {
        x[i] = fv->c[i];
        sumx += x[i];
      }
      x[n_species - 1] = 1.0 - sumx;

      M_mix = 0.0;
      for (i = 0; i < n_species; i++) {
        M_mix += M[i] * x[i];
      }
      c = rho / M_mix;

      sumk = 0.0;
      for (i = 0; i < n_species - 1; i++) {
        sumk += z[i] * z[i] * D[i] * x[i];
      }
      k = FRT * c * sumk;

      var = TEMPERATURE;
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        dkdT[j] = -bf[var]->phi[j] * k / T;
      }

      var = VOLTAGE;
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
        dkdV[j] = 0.0;
      }

      if (pd->v[pg->imtrx][MESH_DISPLACEMENT1]) {
        for (a = 0; a < dim; a++) {
          var = MESH_DISPLACEMENT1 + a;
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            dkdX[a][j] = 0.0;
          }
        }
      }

      /* Derivatives of electrical conductivity (actually kappa/F) with respect to mole fraction
       */

      if (pd->v[pg->imtrx][MASS_FRACTION]) {
        for (w = 0; w < pd->Num_Species_Eqn; w++) {
          var = MASS_FRACTION;
          var_offset = MAX_VARIABLE_TYPES + w;
          for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
            dkdC[w][j] =
                FRT * c * bf[var]->phi[j] * z[w] * z[w] * D[w] +
                k * (d_rho->C[w][j] / rho - bf[var]->phi[j] * (M[w] - M[n_species - 1]) / M_mix);
          }
        }
      }

    } /* end of cr->MassFluxModel == FICKIAN_CHARGED_X */

  } /* end of if (mp->Elec_ConductivityModel == ELECTRONEUTRALITY_FICKIAN) etc. */

  else if (mp->Elec_ConductivityModel == CONSTANT) {

    /*
     * Use Voltage Formulation to choose between electrical conductivity
     * (default) or permittivity (must have been specified).
     */
    if (mp->VoltageFormulation == V_PERMITTIVITY) {
      k = mp->permittivity;
    } else {
      k = mp->electrical_conductivity;
    }

    var = TEMPERATURE;
    for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
      dkdT[j] = 0.;
    }

    var = VOLTAGE;
    for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
      dkdV[j] = 0.0;
    }
    if (pd->v[pg->imtrx][MESH_DISPLACEMENT1]) {
      for (a = 0; a < dim; a++) {
        var = MESH_DISPLACEMENT1 + a;
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          dkdX[a][j] = 0.0;
        }
      }
    }
    if (pd->v[pg->imtrx][MASS_FRACTION]) {
      for (w = 0; w < pd->Num_Species_Eqn; w++) {
        var = MASS_FRACTION;
        var_offset = MAX_VARIABLE_TYPES + w;
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          dkdC[w][j] = 0.0;
        }
      }
    }

    /*Sensitivities were already set to zero */
  } else if (mp->Elec_ConductivityModel == LEVEL_SET) {
    /*
     * Use Voltage Formulation to choose between electrical conductivity
     * (default) or permittivity (must have been specified).
     */
    if (mp->VoltageFormulation == V_PERMITTIVITY) {
      level_set_property(mp->u_electrical_conductivity[0], mp->u_electrical_conductivity[1],
                         mp->u_electrical_conductivity[2], &k, NULL);
    } else {
      k = mp->electrical_conductivity;
    }

    var = TEMPERATURE;
    for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
      dkdT[j] = 0.;
    }

    var = VOLTAGE;
    for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
      dkdV[j] = 0.0;
    }
    if (pd->v[pg->imtrx][MESH_DISPLACEMENT1]) {
      for (a = 0; a < dim; a++) {
        var = MESH_DISPLACEMENT1 + a;
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          dkdX[a][j] = 0.0;
        }
      }
    }
    if (pd->v[pg->imtrx][MASS_FRACTION]) {
      for (w = 0; w < pd->Num_Species_Eqn; w++) {
        var = MASS_FRACTION;
        var_offset = MAX_VARIABLE_TYPES + w;
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          dkdC[w][j] = 0.0;
        }
      }
    }

    /*Sensitivities were already set to zero */
  } else if (mp->Elec_ConductivityModel == EXTERNAL_FIELD) {
    i_elec_cond = mp->elec_cond_external_field;
    GOMA_EH(i_elec_cond, "Electrical conductivity external field not found!");
    k = fv->external_field[i_elec_cond] * mp->u_electrical_conductivity[0];

    // Maybe all sensitivies should be set to 0 by default, and changed
    // only as needed?
    var = TEMPERATURE;
    for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
      dkdT[j] = 0.;
    }

    var = VOLTAGE;
    for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
      dkdV[j] = 0.0;
    }
    if (pd->v[pg->imtrx][MESH_DISPLACEMENT1]) {
      for (a = 0; a < dim; a++) {
        var = MESH_DISPLACEMENT1 + a;
        for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
          dkdX[a][j] = 0.0;
        }
      }
    }
  } else {
    GOMA_EH(GOMA_ERROR, "Unrecognized electrical conductivity model");
  }
  return (k);
}
/*   end of electrical conductivity evaluation		*/

/* end of mm_fill_potential.c */
