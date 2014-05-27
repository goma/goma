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

#include <stdlib.h>
#include <stdio.h>
#include <strings.h>
#include <string.h>
#include <math.h>

/* GOMA include files */

#include "std.h"
#include "rf_fem_const.h"
#include "rf_fem.h"
#include "rf_masks.h"
#include "rf_io_const.h"
#include "rf_io_structs.h"
#include "rf_io.h"
#include "rf_mp.h"
#include "el_elm.h"
#include "el_geom.h"
#include "rf_bc_const.h"
#include "rf_solver_const.h"
#include "rf_fill_const.h"
#include "rf_vars_const.h"
#include "mm_mp_const.h"
#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "mm_as.h"

#include "mm_mp_structs.h"
#include "mm_mp.h"

#include "mm_eh.h"

#include "mm_fill_species.h"
#include "mm_std_models.h"
#include "mm_fill_common.h"

#define _MM_STD_MODELS_C
#include "goma.h"

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
      mp->DensityModel == DENSITY_FOAM_TIME_TEMP) {
    calculateVolumeFractionGasPhase(time);
  }
}
