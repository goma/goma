/************************************************************************ *
* Goma - Multiphysics finite element software                             *
* Sandia National Laboratories                                            *
*                                                                         *
* Copyright (c) 2014 Sandia Corporation.                                  *
*                                                                         *
* Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,  *
* the U.S. Government retains certain rights in this software.            *
*                                                                         *
* This software is distributed under the GNU General Public License.      *
\************************************************************************/

/* This routine contains all the viscosity models along with the routine 
 * that calls them. If you add a new viscosity model, it would go in here.
 */

/* Standard include files */

#include <stdio.h>

/* GOMA include files */

#include "std.h"
#include "rf_fem_const.h"
#include "mm_mp_const.h"
#include "mm_as_structs.h"
#include "mm_as.h"
#include "mm_mp_structs.h"
#include "mm_mp.h"
#include "mm_eh.h"
#include "mm_fill_ls.h"
#include "rf_allo.h"

#define GOMA_MM_DIL_VISCOSITY_C
#include "mm_dil_viscosity.h"

static void transferMultipleOfDerivatives(const dbl ratioVisc,
					  const VISCOSITY_DEPENDENCE_STRUCT *d_mu,
					  DILVISCOSITY_DEPENDENCE_STRUCT *d_dilMu) {
  int j, a, w, var;
  int dim = ei[pg->imtrx]->ielem_dim;
  if (pd->v[pg->imtrx][TEMPERATURE]) {
    for (j = 0; j < ei[pg->imtrx]->dof[TEMPERATURE]; j++) {
      d_dilMu->T[j] = ratioVisc * d_mu->T[j]; 
    }
  }
  if (pd->v[pg->imtrx][FILL]) {
    for (j = 0; j < ei[pg->imtrx]->dof[FILL]; j++) {
      d_dilMu->F[j] = ratioVisc * d_mu->F[j]; 
    }
  }
  if (pd->v[pg->imtrx][MESH_DISPLACEMENT1]) {
    for (a = 0; a < dim; a++) {
      var = MESH_DISPLACEMENT1 + a;
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
	d_dilMu->X[a][j] = ratioVisc * d_mu->X[a][j]; 
      }
    }
  }
  if (pd->v[pg->imtrx][MASS_FRACTION]) {
    for (w = 0; w < pd->Num_Species_Eqn; w++) {
      for (j = 0; j < ei[pg->imtrx]->dof[MASS_FRACTION]; j++) {
	d_dilMu->C[w][j] = ratioVisc * d_mu->C[w][j];
      }
    }
  }
  if (pd->v[pg->imtrx][VELOCITY1]) {
    for (a = 0; a < dim; a++) {
      var = VELOCITY1 + a;
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
	d_dilMu->v[a][j] = ratioVisc * d_mu->v[a][j]; 
      }
    }
  }
  if (pd->v[pg->imtrx][PRESSURE]) {
    for (j = 0; j < ei[pg->imtrx]->dof[PRESSURE]; j++) {
      d_dilMu->P[j] = ratioVisc * d_mu->P[j]; 
    }
  }

  if (pd->v[pg->imtrx][PHASE1]) {
    for (a = 0; a < pfd->num_phase_funcs; a++) {
      var = PHASE1 + a;
      for(j = 0 ; j < ei[pg->imtrx]->dof[var] ; j++) {
	d_dilMu->pf[a][j] = ratioVisc * d_mu->pf[a][j];
      }
    }
  }

#ifdef COUPLED_FILL
  for (j = 0; j < ei[pg->imtrx]->dof[PRESSURE]; j++) {
    d_dilMu->F[j] = ratioVisc * d_mu->F[j]; 
  }
#endif

  if (pd->v[pg->imtrx][BOND_EVOLUTION]) {
    for (j = 0; j < ei[pg->imtrx]->dof[BOND_EVOLUTION]; j++) {
      d_dilMu->nn[j] = ratioVisc * d_mu->nn[j]; 
    }
  }

}


//! Transfer Gauss Point Derivatives to a viscosity dependence
//! structure 
/*!
 *  Transfer the derivatives of a quantity evaluated at a gauss point to 
 *  the  DILVISCOSITY_DEPENDENCE_STRUCT dependency. This dependency
 *  is on a per independent variable basis. 
 *
 *    mu depends on a bunch of variables evaluated at the gauss point:
 *             mu = func(indVar_i_gp)
 *    then the independent variables depend
 *          where:
 *             indVar_i_gp = sum{j] ( bf_j * indVar_i_j)
 *
 * NOTE: We have just included temperature and MF dependence here. As, that's
 *       all that is necessary for the current problem.
 */
static void transferGPDerivatives(const dbl multFac,
				  const dbl * gpDerivatives,
				  DILVISCOSITY_DEPENDENCE_STRUCT *d_dilMu) {
  int j, w;
  int var = TEMPERATURE;
  if (pd->v[pg->imtrx][var]) {
    for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
      d_dilMu->T[j] += multFac * gpDerivatives[var] * bf[var]->phi[j];
    }
  }
  var = MASS_FRACTION;
  if (pd->v[pg->imtrx][var])  {
    for (w = 0; w < pd->Num_Species_Eqn; w++) {
      for (j = 0; j < ei[pg->imtrx]->dof[var]; j++)  {
	d_dilMu->C[w][j] += multFac * gpDerivatives[MAX_VARIABLE_TYPES + w] * bf[var]->phi[j];
      }
    }
  }
}

/*******************************************************************************
 * dil_viscosity(): 
 *
 *              Calculate the viscosity and derivatives of viscosity
 *              with respect to solution unknowns at the Gauss point. Most 
 *              non-Newtonian purely viscous models depend on the shear rate
 *              invariant and may also depend on temperature and pressure.
 *
 * Input
 * -----
 *   gn_local	= Pointer to the Generalized_Newtonian structure for solvent or polymer
 *   mu		= viscosity
 *   gamma_dot	= Strain rate tensor
 *   d_mu       = dependence of viscosity, where:
 *
 *   d_mu->gd	= derivative of viscosity wrt to strain rate invariant (which one?)
 *   d_mu->v	= derivative of viscosity wrt to velocity
 *   d_mu->X	= derivative of viscosity wrt to mesh
 *   d_mu->T	= derivative of viscosity wrt to temperature
 *   d_mu->P	= derivative of viscosity wrt to pressure
 *   d_mu->C	= derivative of viscosity wrt to concentration/species
 *   d_mu->F	= derivative of viscosity wrt to FILL (Level Set / VoF)
 *   d_mu->nn	= derivative of viscosity wrt to bond concentration
 *
 *  This routine takes care of zeroing out the d_mu structure every time it is
 *  called with a nonnull d_mu.
 *
 *******************************************************************************/
double
dil_viscosity (GEN_NEWT_STRUCT *gn_local, const dbl muValue, const VISCOSITY_DEPENDENCE_STRUCT *d_mu, DILVISCOSITY_DEPENDENCE_STRUCT *d_dilMu)
{

  int  w, j, var;
  dbl kappa = 0.0;
  dbl ratioVisc = mp->dilationalViscosityRatio;
  if (d_dilMu != 0) {
    if (d_mu == 0) {
      EH(GOMA_ERROR,"shouldn't be here");
    }
  }

  /* Zero out all of the sensitivities */
  zeroStructures(d_dilMu, 1);
  
  /* this section is for all Newtonian models */

  if (gn_local->ConstitutiveEquation == NEWTONIAN) {
  
    if(mp->DilationalViscosityModel == USER ) {
      EH(GOMA_ERROR,"User Dilational Viscosity Model is Unimplemented");
    } else if (mp->DilationalViscosityModel == USER_GEN) {
      EH(GOMA_ERROR,"UserGen Dilational Viscosity Model is Unimplemented");
    } else if (mp->DilationalViscosityModel == FILL) {
      EH(GOMA_ERROR," Dilational Viscosity Model for FILL is Unimplemented");
    } else if (mp->DilationalViscosityModel == TABLE) {
      EH(GOMA_ERROR," Dilational Viscosity Model for TABLE is Unimplemented");
    } else if (mp->DilationalViscosityModel == LEVEL_SET) {
      EH(GOMA_ERROR, "Dilational Viscosity Model Unimplemented for LEVEL_SET");
    } else if (mp->DilationalViscosityModel == DILVISCM_KAPPAWIPESMU ) {
      // Don't need to fill in dependencies since the terms disappear
      kappa = 2.0 *  muValue / 3.0;
    } else if (mp->DilationalViscosityModel == DILVISCM_KAPPACONSTANT) {
      kappa = mp->dilationalViscosity;
    } else if (mp->DilationalViscosityModel == DILVISCM_KAPPAFIXEDRATIO) {

      kappa = ratioVisc * muValue;
      if (d_dilMu != 0) {
	transferMultipleOfDerivatives(ratioVisc, d_mu, d_dilMu);
      }

    } else if (mp->DilationalViscosityModel == DILVISCM_KAPPABUBBLES) {
      EH(GOMA_ERROR," Newtonian model with KappaBubbles dil visc doesn't make sense");
  
    } else {
      EH(GOMA_ERROR,"Unrecognized viscosity model for Newtonian fluid");
    }

      /*       return(status); */
  } /* end Newtonian section */

  /*
   *    CONSTANT consitutive model. We can handle DILVISCM_KAPPAWIPESMU 
   *    and DILVISCM_KAPPACONSTANT cases here
   */
  else if (gn_local->ConstitutiveEquation == CONSTANT) {
    if (mp->DilationalViscosityModel == DILVISCM_KAPPAWIPESMU) {
      kappa = 2.0 *  muValue / 3.0;
    } else if (mp->DilationalViscosityModel == DILVISCM_KAPPACONSTANT) {
      kappa = mp->dilationalViscosity;
    } else {
      EH(GOMA_ERROR, "unsupported Kappa Option");
    }
  } 

  else if (gn_local->ConstitutiveEquation == FOAM_EPOXY) {

   if (mp->DilationalViscosityModel == DILVISCM_KAPPAWIPESMU) {
      kappa = 2.0 *  muValue / 3.0;
    } else if (mp->DilationalViscosityModel == DILVISCM_KAPPACONSTANT) {
      kappa = mp->dilationalViscosity;
    } else if (mp->DilationalViscosityModel == DILVISCM_KAPPAFIXEDRATIO) {
      kappa = ratioVisc * muValue;
      if (d_dilMu != 0) {
	transferMultipleOfDerivatives(ratioVisc, d_mu, d_dilMu);
      }
    } else if (mp->DilationalViscosityModel == DILVISCM_KAPPABUBBLES) {
      /*
       * implement kappa 4/3 mu_L * (1 - volF) / volF
       *
       * Note, this is mu_L, the pure liquid phase viscosity. So, we can't
       * just use mu and d_mu, because it contains the dependence of the
       * mixture viscosity on the volume fraction. 
       */
      /*
       * Obtain a cropped value of the volume fraction of air in the foam. Call it volF,
       * and make sure it isn't zero. Here we restrict volF to being greater than
       * 1.0E-4.
       */
      double volF;
      if (mp->volumeFractionGas < 1.0E-4) {
	volF = 1.0E-4;
      } else if (mp->volumeFractionGas < 1.0) {
	volF = mp->volumeFractionGas;
      } else {
	volF = 1.0;
      }
     
  
      // Get the pure liquid phase viscosity from somewhere
      double muLValue = mp->FlowingLiquid_viscosity;
   
      double ratio = 4. / 3. * (1.0 - volF ) / volF;
      kappa = ratio * muLValue;
      
      if (d_dilMu != 0) {
	// Ok, to get the derivatives, we copy a multiple of the dependencies from the pure 
	// species viscosities into the kappa dependencies
	transferGPDerivatives(ratio, mp->d_FlowingLiquid_viscosity, d_dilMu);

	// Then, we add in the explicit (1 - volF) / volF dependency, which only
	// depends on the concentration unknowns, and we are done.
	
	var = MASS_FRACTION;
	if (pd->v[pg->imtrx][var]) {
	  double tmp = 4. * muLValue / 3. / (volF * volF);
	  double * dVolFdMF = &(mp->d_volumeFractionGas[0]) + MAX_VARIABLE_TYPES;
	  for (w = 0; w < pd->Num_Species_Eqn; w++) {
	    for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
	      d_dilMu->C[w][j] -= tmp * dVolFdMF[w] * (bf[var]->phi[j]);
	    }
	  }
	}
	
      }

    } else {
      EH(GOMA_ERROR, "unsupported Kappa Option");
    }

  }
  else if (gn_local->ConstitutiveEquation == FOAM_PMDI_10) {

    if (mp->DilationalViscosityModel == DILVISCM_KAPPABUBBLES) {
      //      double muL = mp->FlowingLiquid_viscosity;
      double volF = mp->volumeFractionGas;

      if (volF < 1e-6) {
	volF = 1e-6;
      }

      double ratio = 4. / 3. * (1.0 - volF ) / volF;
      kappa = ratio * muValue ;

      if (d_dilMu != 0) {
	// Ok, to get the derivatives, we copy a multiple of the dependencies from the pure
	// species viscosities into the kappa dependencies
	transferGPDerivatives(ratio, mp->d_FlowingLiquid_viscosity, d_dilMu);

	// Then, we add in the explicit (1 - volF) / volF dependency, which only
	// depends on the concentration unknowns, and we are done.

	var = MASS_FRACTION;
	if (pd->v[pg->imtrx][var]) {
	  double tmp = 4. * muValue/ 3. / (volF * volF);
	  double * dVolFdMF = &(mp->d_volumeFractionGas[0]) + MAX_VARIABLE_TYPES;
	  for (w = 0; w < pd->Num_Species_Eqn; w++) {
	    for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
	      d_dilMu->C[w][j] -= ratio * d_mu->C[w][j]  + tmp * dVolFdMF[w] * (bf[var]->phi[j]);
	    }
	  }
	}

	var = TEMPERATURE;
	if (pd->v[pg->imtrx][var]) {
	  double tmp = 4. * muValue/ 3. / (volF * volF);
	  for (w = 0; w < pd->Num_Species_Eqn; w++) {
	    for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
	      d_dilMu->T[j] -= ratio * d_mu->T[j]  + tmp * mp->d_volumeFractionGas[var] * (bf[var]->phi[j]);
	    }
	  }
	}

      }

      if (ls != NULL)
        {
          ls_modulate_dilviscosity(&kappa, 0, ls->Length_Scale, mp->mp2nd->viscositymask[0],
                                   mp->mp2nd->viscositymask[1], d_dilMu);
        }


    } else {
      EH(GOMA_ERROR, "unsupported Kappa Option");
    }

  }
  /*
   *     For the Other constitutive models, the default is  DILVISCM_KAPPAWIPESMU
   *     Anything else, we throw an error, because we haven't handled it
   *     yet.
   */
  else {
    if (mp->DilationalViscosityModel == DILVISCM_KAPPAWIPESMU) {
      kappa = 2.0 * muValue / 3.0;
    } else if (mp->DilationalViscosityModel == DILVISCM_KAPPAFIXEDRATIO) {
      kappa = ratioVisc * muValue;
      if (d_dilMu != 0) {
	transferMultipleOfDerivatives(ratioVisc, d_mu, d_dilMu);
      }
    } else {
      EH(GOMA_ERROR, "unsupported Kappa Option");
    }
  }

  return(kappa);
}

int
ls_modulate_dilviscosity ( double *kappa1,
                        double  kappa2,
                        double width,
                        double pm_minus,
                        double pm_plus,
                        DILVISCOSITY_DEPENDENCE_STRUCT *d_dilMu)
{
  double factor;
  int i,a, w, var;

  if ( d_dilMu == NULL )
    {
      *kappa1 = ls_modulate_property( *kappa1, kappa2, width, pm_minus, pm_plus, NULL, &factor);
      return(1);
    }

  *kappa1 = ls_modulate_property( *kappa1, kappa2, width, pm_minus, pm_plus, d_dilMu->F, &factor);

  d_dilMu->gd *= factor;

  if ( pd->v[pg->imtrx][var=TEMPERATURE ] )
    {
      for(i=0; i<ei[pg->imtrx]->dof[var]; i++)
        {
          d_dilMu->T[i] *= factor;
        }
    }

  if ( pd->v[pg->imtrx][var=MASS_FRACTION ] )
    {
      for ( w=0; w<pd->Num_Species_Eqn; w++)
        {
          for( i=0; i<ei[pg->imtrx]->dof[var]; i++)
            {
              d_dilMu->C[w][i] *= factor;
            }
        }
    }

  if( pd->v[pg->imtrx][var=VELOCITY1] )
    {
      for( a=0; a<pd->Num_Dim; a++ )
        {
          for( i=0; i<ei[pg->imtrx]->dof[var]; i++ )
            {
              d_dilMu->v[a][i] *= factor;
            }
        }
    }

  if( pd->v[pg->imtrx][var=MESH_DISPLACEMENT1] )
    {
      for( a=0; a<pd->Num_Dim; a++ )
        {
          for( i=0; i<ei[pg->imtrx]->dof[var]; i++ )
            {
              d_dilMu->X [a][i] *= factor;
            }
        }
    }

  if( pd->v[pg->imtrx][var=PRESSURE] )
    {
      for( i=0; i<ei[pg->imtrx]->dof[var]; i++ )
        {
          d_dilMu->P[i] *= factor;
        }
    }

  if ( pd->v[pg->imtrx][var=BOND_EVOLUTION] )
    {
      for( i=0 ; i<ei[pg->imtrx]->dof[var]; i++)
        {
          d_dilMu->nn[i] *= factor;
        }
    }
  return ( 1 );

}