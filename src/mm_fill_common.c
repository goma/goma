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
 
/*
 *$Id: mm_fill_common.c,v 5.2 2008-03-22 00:55:49 hkmoffa Exp $
 */

#ifdef USE_RCSID
static char rcsid[] = "$Id: mm_fill_common.c,v 5.2 2008-03-22 00:55:49 hkmoffa Exp $";
#endif

#include <math.h>

/* GOMA include files */

#include "std.h"
#include "rf_fem_const.h"
#include "mm_mp_const.h"
#include "mm_as_structs.h"
#include "mm_as.h"
#include "mm_mp_structs.h"
#include "mm_mp.h"
#include "mm_eh.h"
#include "mm_fill_common.h"

/*********** R O U T I N E S   I N   T H I S   F I L E ************************
*
*       NAME            TYPE            CALLED_BY
*    ------------             ---------               --------------
* computeCommonMaterialProps_gp()                   
* 
******************************************************************************/
/*
 * This file contains all common models
******************************************************************************/


static void calculateVolumeFractionGasPhase(const dbl time) {
  dbl rho; 
  if (mp->DensityModel == DENSITY_FOAM)
    {
      EH(-1, "Not completed");
    }
  else if (mp->DensityModel == DENSITY_FOAM_CONC)
    {
      EH(-1, "Not completed");
    }
  else if (mp->DensityModel == DENSITY_FOAM_TIME)
    {
      double rho_init, rho_final, aexp, time_delay, realtime;
      rho_init   = mp->u_density[0]; /* Initial density */
      rho_final  = mp->u_density[1]; /* final density */
      aexp       = mp->u_density[2]; /* Arhennius constant for time exponent*/
      time_delay = mp->u_density[3]; /* time delay before foaming starts */
      
      if (time > time_delay)
	{
	  realtime = time - time_delay;
	  rho = rho_final + (rho_init-rho_final) * exp(-aexp*realtime);
	}
      else
	{
	  rho = rho_init;
	}        
      double theta = rho / rho_init;
      double MolecWeightGas = 30.;
      double rhoGas = MolecWeightGas / (82.05 * 300.);
      mp->volumeFractionGas = rho_init / (rho_init - rhoGas) * (1.0 - theta);
    }
  else if (mp->DensityModel == DENSITY_FOAM_TIME_TEMP)
    {
      double T;
      double rho_init  = mp->u_density[0];  /* Initial density (units of gm cm-3 */
      double rho_final = mp->u_density[1];  /* Final density (units of gm cm-3 */
      double cexp      = mp->u_density[2];  /* cexp in units of Kelvin sec */
      double coffset   = mp->u_density[3];  /* offset in units of seconds */
      double time_delay= mp->u_density[4];  /* time delay in units of seconds */
      double drhoDT;
      T = fv->T;
      if (time > time_delay)
	{
	  double realtime = time - time_delay;
	  double delRho = (rho_init - rho_final);
	  double cdenom = cexp - coffset * T;
	  double expT = exp(-realtime * T / cdenom);
	  rho = rho_final + delRho * expT;
	  drhoDT = delRho * expT *(- realtime * cexp / (cdenom * cdenom));
	}
      else
	{
	  rho = rho_init;
	  drhoDT = 0.0;
	}
      double theta = rho / rho_init;
      double MolecWeightGas = 30.;
      // Have to make rhoGas independent of T, because rho_final is independent of T
      double rhoGas = MolecWeightGas / (82.05 * 300.);
      mp->volumeFractionGas = rho_init / (rho_init - rhoGas) * (1.0 - theta);
      mp->d_volumeFractionGas[TEMPERATURE] = - drhoDT / (rho_init - rhoGas);
    }
  else if (mp->DensityModel == DENSITY_FOAM_PMDI_10)
    {
      if (pd->gv[MOMENT1]) {
	mp->volumeFractionGas = (fv->moment[1] / (1 + fv->moment[1]));
	mp->d_volumeFractionGas[MOMENT1] = 1 / ((1 + fv->moment[1])* (1 + fv->moment[1]));
      } else {
	int wCO2;
	int wH2O;
	int w;

	for (w = 0; w < MAX_VARIABLE_TYPES+MAX_CONC; w++) {
	  mp->d_volumeFractionGas[w] = 0.0;
	}

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
	  EH(-1, "Expected a Species Source of FOAM_PMDI_10_CO2");
	} else if (wH2O == -1) {
	  EH(-1, "Expected a Species Source of FOAM_PMDI_10_H2O");
	}

	double M_CO2 = mp->u_density[0];
	//double rho_liq = mp->u_density[1];
	double ref_press = mp->u_density[2];
	double Rgas_const = mp->u_density[3];
	double rho_gas = 0;
	double nu;
	double d_nu_dC;
	double d_nu_dT;

	double phi = 0.0;
	double d_phi_dC = 0.0;
	double d_phi_dT = 0.0;


	if (fv->T > 0) {
	  rho_gas = (ref_press * M_CO2 / (Rgas_const * fv->T));
	  nu = M_CO2 * fv->c[wCO2] / rho_gas;
	  d_nu_dC = M_CO2 / rho_gas;
	  d_nu_dT = M_CO2 * fv->c[wCO2] * Rgas_const / (ref_press * M_CO2);

	  phi = nu / (1 + nu);
	  d_phi_dC = (d_nu_dC) / ((1 + nu)*(1 + nu));
	  d_phi_dT = (d_nu_dT) / ((1 + nu)*(1 + nu));
	}


	mp->volumeFractionGas = phi;
	mp->d_volumeFractionGas[TEMPERATURE] = d_phi_dT;

	mp->d_volumeFractionGas[MAX_VARIABLE_TYPES+wCO2] = d_phi_dC;
      }
    }
  else if (mp->DensityModel == DENSITY_MOMENT_BASED)
    {
    mp->volumeFractionGas = (fv->moment[1] / (1 + fv->moment[1]));
    mp->d_volumeFractionGas[MOMENT1] = 1 / ((1 + fv->moment[1])* (1 + fv->moment[1]));
    }
}

/*
 * Compute common material properties at the gauss point
 *
 * These get loaed into the material property structure.
 */
void computeCommonMaterialProps_gp(const dbl time)
{
  /*
   * Decide whether we have a two-phase flow and need to 
   * calculate the volume fraction
   */
  if (mp->DensityModel == DENSITY_FOAM ||
      mp->DensityModel == DENSITY_FOAM_CONC ||
      mp->DensityModel == DENSITY_FOAM_TIME ||
      mp->DensityModel == DENSITY_FOAM_TIME_TEMP ||
          mp->DensityModel == DENSITY_MOMENT_BASED ||
      mp->DensityModel == DENSITY_FOAM_PMDI_10) {
    calculateVolumeFractionGasPhase(time);
  }
}
