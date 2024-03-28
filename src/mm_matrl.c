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
 *$Id: mm_matrl.c,v 5.10 2008-05-02 19:07:56 hkmoffa Exp $
 */

/*************** R O U T I N E S   I N   T H E   F I L E ***********************
 *
 *    NAME				TYPE		CALLED_BY
 *--------------------------------------------------------------------
 *
 *    chemkin_mat_prop_init
 ******************************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "density.h"
#include "mm_as.h"
#include "mm_as_structs.h"
#include "mm_eh.h"
#include "mm_fill_population.h"
#include "mm_fill_terms.h"
#include "mm_input.h"
#include "mm_mp.h"
#include "mm_mp_const.h"
#include "mm_mp_structs.h"
#include "mm_qtensor_model.h"
#include "mm_species.h"
#include "rd_mesh.h"
#include "rf_allo.h"
#include "rf_bc.h"
#include "rf_bc_const.h"
#include "rf_fem.h"
#include "rf_fem_const.h"
#include "rf_solve.h"
#include "std.h"
#include "user_mp.h"

#ifdef USE_CHEMKIN
#include "ck_chemkin_const.h"
#include "cpc_defs.h"
#endif

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void print_char(char ch, int num)

/*************************************************************************
 *
 * print_char:
 *
 * Repeatedly print a single character
 *************************************************************************/
{
  int i;
  for (i = 0; i < num; i++)
    (void)printf("%c", ch);
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void StrngDefaultDatabase(char *str, int ID)

/*************************************************************************
 *
 * StrngDefaultDatabase:
 *
 * StrngDefaultDatabase: translates an integer ID into a string
 * description
 *************************************************************************/
{
  switch (ID) {
  case DB_GOMA_MAT:
    strcpy(str, "Goma Material Property Default");
    break;
  case DB_CHEMKIN_MAT:
    strcpy(str, "Chemkin Material Property Default");
    break;
  default:
    strcpy(str, "Unknown Material Database Default");
  }
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void matrl_prop_print(MATRL_PROP_STRUCT *mat_ptr, int mn)

/**************************************************************************
 *
 * matrl_prop_print:
 *
 *  Parameter list:
 *
 *  mat_ptr     = Pointer to the materials property to be printed out.
 *  mn          = Integer ID of the material
 *
 *  This utility string will print out the current definition of the
 *  all of the material properties for a given material.
 *
 *************************************************************************/
{
  int i;
  char pString[256];
  PROBLEM_DESCRIPTION_STRUCT *pd_ptr = pd_glob[mn];

  /* Code cleanup crew is not happy with the unimportant information provided by this function to
   * the terminal.  In the calling program we now make this call dependent on the Debug_Flag level
   */

  printf("\n\n");
  print_char('=', 80);
  printf("\n\n");
  printf("\tMaterial Property Printout for Material %d: %s\n", mn, (char *)mat_ptr->Material_Name);
  StrngDefaultDatabase(pString, mat_ptr->DefaultDatabase);
  printf("\tDefault Database = %s\n", pString);
  printf("\n");

  printf("Media Type = %d\n", mat_ptr->PorousMediaType);
  printf("FlowingViscosity Model = %d\n", mat_ptr->FlowingLiquidViscosityModel);

  if (pd_ptr->Num_Species > 0) {
    printf("\n\t\tDescription of Species in the Material\n");
    printf("\n  Name    Phase  Weight  Species_Var_Type\n");
    print_char('-', 80);
    printf("\n");
    for (i = 0; i < pd_ptr->Num_Species; i++) {
      printf("%s\t%d\t%g\t", mat_ptr->Species_Names[i], mat_ptr->PhaseID[i],
             mat_ptr->molecular_weight[i]);
      /*
       * Convert species var type to a string. Then print it out
       */
      species_type_int_to_str(pString, mat_ptr->Species_Var_Type);
      printf("%s", pString);
      if (mat_ptr->Num_Species_Eqn < mat_ptr->Num_Species) {
        if (i == (mat_ptr->Num_Species_Eqn)) {
          printf("\t(SUM OF FRACTIONS DIRICHLET CONDITION)");
        }
      }
      printf("\n");
    }
    print_char('-', 80);
    printf("\n");
  }
  printf("\n");
  printf("\n");
  printf("\n\n");
  print_char('=', 80);
  printf("\n\n");
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int goma_mat_prop_init(MATRL_PROP_STRUCT *mat_ptr, int mn, PROBLEM_DESCRIPTION_STRUCT *pd_ptr)

/*************************************************************************
 *
 * goma_mat_prop_init():
 *
 *  This does some initializations of the materials property
 *  database, based on the now known number of species. This function
 *  is done in lieu of calling a chemkin initialization and thus
 *  should perform functions done there as well.
 *
 *   Input
 *  ---------
 *  mat_ptr -> pointer to the current material problem structure
 *             The number of species in the material should already
 *             have been input
 *  mn         Index number for the material
 *  pd_ptr ->  Pointer to the problem description structure.
 *
 *   Return
 *    1 = Everything went ok
 *   -1 = Something went terribly amiss.
 ************************************************************************/
{
  int i;
  char pstring[256];
  /*
   *  Malloc space for default names for the Species in the material
   *  structure. We will malloc one more than is necessary, because in
   *  some instances the number of species is increased by one.
   */
  if (mat_ptr->Species_Names == NULL) {
    mat_ptr->Species_Names = alloc_VecFixedStrings(mat_ptr->Num_Species + 1, sizeof(CK_NAME_STR));
    if (mat_ptr->Species_Names != NULL) {
      for (i = 0; i < mat_ptr->Num_Species + 1; i++) {
        (void)sprintf(pstring, "Species_%d", i);
        (void)strcpy(mat_ptr->Species_Names[i], pstring);
      }
    } else {
      return -1;
    }
  }
  /*
   *  Add other defaults based on a loop over the species
   */
  if (mat_ptr->Species_Var_Type == 0) {
    mat_ptr->Species_Var_Type = pd_ptr->Species_Var_Type;
  }
  mat_ptr->StateVector_speciesVT = mat_ptr->Species_Var_Type;
  for (i = 0; i < mat_ptr->Num_Species; i++) {
    mat_ptr->PhaseID[i] = 0;
    mat_ptr->Volumetric_Dirichlet_Cond[i] = 0;
  }

  /*
   *  Malloc space for default names for the porous media phases in the material
   *  structure.
   */

  if (mat_ptr->Porous_Names == NULL) {
    mat_ptr->Porous_Names = alloc_VecFixedStrings(mat_ptr->Num_Porous_Eqn, sizeof(CK_NAME_STR));
    if (mat_ptr->Porous_Names != NULL) {
      for (i = 0; i < mat_ptr->Num_Porous_Eqn; i++) {
        (void)sprintf(pstring, "Porous_Phase_%d", i);
        (void)strcpy(mat_ptr->Porous_Names[i], pstring);
      }
    } else {
      return -1;
    }
  }

  return 1;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void reconcile_bc_to_matrl(void)

/*************************************************************************
 *
 * reconcile_bc_to_matrl():
 *
 *  This does some initializations of the materials property
 *  database, based on the now known number of species. This function
 *  is done in lieu of calling a chemkin initialization and thus
 *  should perform functions done there as well.
 *
 *   Input
 *  ---------
 *  mat_ptr -> pointer to the current material problem structure
 *             The number of species in the material should already
 *             have been input
 *  mn         Index number for the material
 *  pd_ptr ->  Pointer to the problem description structure.
 *
 *   Return
 *    1 = Everything went ok
 *   -1 = Something went terribly amiss.
 ************************************************************************/
{
  int ibc, eb_index, matID, num_species;
  static char *yo = "reconcile_bc_to_matrl: ";
  BOUNDARY_CONDITION_STRUCT *bc;
  MATRL_PROP_STRUCT *matrl_ptr;

  for (ibc = 0; ibc < Num_BC; ibc++) {
    /*
     *  Create a couple of pointers to cut down on the
     *  amount of indirect addressing
     */
    bc = BC_Types + ibc;

    switch (bc->BC_Name) {

    case VL_EQUIL_BC:
      /*
       *  Check the pressure datum
       */
      if (bc->BC_Data_Float[0] != upd->Pressure_Datum) {
        fprintf(stderr,
                "%s WARNING Global Pressure Datum %g and VL_EQUIL Pressure Datum %g Differ\n", yo,
                upd->Pressure_Datum, bc->BC_Data_Float[0]);
        if (upd->Pressure_Datum != 1.0132500000E6) {
          fprintf(stderr, "%s ERROR Global Pressure Datum set differently in two places\n", yo);
          GOMA_EH(GOMA_ERROR, "Duplicate and different pressure datums");
        } else {
          fprintf(stderr, "%s WARNING Global Pressure Datum set to %g by VL_EQUIL card\n", yo,
                  bc->BC_Data_Float[0]);
          upd->Pressure_Datum = bc->BC_Data_Float[0];
        }
        /*
         * Check the molecular weights
         */
        eb_index = bc->BC_Data_Int[1];
        matID = map_mat_index(eb_index);
        if (matID < 0) {
          GOMA_EH(GOMA_ERROR, "matID not found for eb_index");
        }
        matrl_ptr = mp_glob[matID];
        num_species = matrl_ptr->Num_Species;
        if (bc->BC_Data_Float[1] != matrl_ptr->molecular_weight[0]) {
          if (matrl_ptr->molecular_weight[0] == -1.0) {
            matrl_ptr->molecular_weight[0] = bc->BC_Data_Float[1];
          } else {
            fprintf(stderr, "%s ERROR molecular weights for species 0 differ: %g %g\n", yo,
                    bc->BC_Data_Float[1], matrl_ptr->molecular_weight[0]);
            GOMA_EH(GOMA_ERROR, "Duplicate and different molecular weights");
          }
        }
        if (bc->BC_Data_Float[2] > 0.0) {
          if (bc->BC_Data_Float[2] != matrl_ptr->molecular_weight[1]) {
            if (matrl_ptr->molecular_weight[1] == -1.0) {
              matrl_ptr->molecular_weight[1] = bc->BC_Data_Float[2];
            } else {
              fprintf(stderr, "%s ERROR molecular weights for species 0 differ: %g %g\n", yo,
                      bc->BC_Data_Float[2], matrl_ptr->molecular_weight[1]);
              GOMA_EH(GOMA_ERROR, "Duplicate and different molecular weights");
            }
          }
        }
        if (bc->BC_Data_Float[3] != matrl_ptr->molecular_weight[num_species - 1]) {
          if (matrl_ptr->molecular_weight[num_species - 1] == -1.0) {
            matrl_ptr->molecular_weight[num_species - 1] = bc->BC_Data_Float[3];
          } else {
            fprintf(stderr, "%s ERROR molecular weights for species %d differ: %g %g\n", yo,
                    num_species - 1, bc->BC_Data_Float[3],
                    matrl_ptr->molecular_weight[num_species - 1]);
            GOMA_EH(GOMA_ERROR, "Duplicate and different molecular weights");
          }
        }

        /*
         * Check the molecular weights for the second phase, the gas phase
         */
        eb_index = bc->BC_Data_Int[2];
        matID = map_mat_index(eb_index);
        matrl_ptr = mp_glob[matID];
        num_species = matrl_ptr->Num_Species;
        if (bc->BC_Data_Float[1] != matrl_ptr->molecular_weight[0]) {
          if (matrl_ptr->molecular_weight[0] == -1.0) {
            matrl_ptr->molecular_weight[0] = bc->BC_Data_Float[1];
          } else {
            fprintf(stderr, "%s ERROR molecular weights for species 0 differ: %g %g\n", yo,
                    bc->BC_Data_Float[1], matrl_ptr->molecular_weight[0]);
            GOMA_EH(GOMA_ERROR, "Duplicate and different molecular weights");
          }
        }
        if (bc->BC_Data_Float[2] > 0.0) {
          if (bc->BC_Data_Float[2] != matrl_ptr->molecular_weight[1]) {
            if (matrl_ptr->molecular_weight[1] == -1.0) {
              matrl_ptr->molecular_weight[1] = bc->BC_Data_Float[2];
            } else {
              fprintf(stderr, "%s ERROR molecular weights for species 1 differ: %g %g\n", yo,
                      bc->BC_Data_Float[2], matrl_ptr->molecular_weight[1]);
              GOMA_EH(GOMA_ERROR, "Duplicate and different molecular weights");
            }
          }
        }
        if (bc->BC_Data_Float[4] != matrl_ptr->molecular_weight[num_species - 1]) {
          if (matrl_ptr->molecular_weight[num_species - 1] == -1.0) {
            matrl_ptr->molecular_weight[num_species - 1] = bc->BC_Data_Float[4];
          } else {
            fprintf(stderr, "%s ERROR molecular weights for species %d differ: %g %g\n", yo,
                    num_species - 1, bc->BC_Data_Float[4],
                    matrl_ptr->molecular_weight[num_species - 1]);
            GOMA_EH(GOMA_ERROR, "Duplicate and different molecular weights");
          }
        }
      }
      break;
    }
  }
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

double
calc_density(MATRL_PROP_STRUCT *matrl, int doJac, PROPERTYJAC_STRUCT *densityJac, double time)

/**************************************************************************
 *
 * calc_density():
 *
 *   Calculate the density and its derivatives using the properties
 *   evaluated from the local materials state vector. Note, we strive
 *   to have a complete description of the materials models in this
 *   function.
 *
 * Output
 * -----
 *    If doJac is true, this routine also fills up the densityJac
 *    structure with derivative information.
 *
 *
 *  Return
 * --------
 *    Actual value of the density.
 ***************************************************************************/
{
  int w;
  int species, num_dep, matID, num_species, num_species_eqn;
  dbl F, vol = 0, rho = 0, rho_f, rho_s, pressureThermo, RGAS_CONST;
  double avgMolecWeight = 0.0, tmp;
  double *mw = matrl->molecular_weight;
  double *stateVector = matrl->StateVector;
  int speciesVT = matrl->StateVector_speciesVT;
  matID = matrl->MatID;
  /*
   * Figure out an initial estimate of the number of dependencies
   * and augment densityJac
   */
  if (densityJac == NULL)
    doJac = FALSE;
  if (doJac) {
    densityJac->Num_Terms = 0;
    if (upd->vp[pg->imtrx][TEMPERATURE] != -1) {
      num_dep = 1;
    } else {
      num_dep = 0;
    }
    if (upd->vp[pg->imtrx][MASS_FRACTION] != -1 || upd->vp[pg->imtrx][SPECIES_UNK_0] != -1) {
      num_dep += matrl->Num_Species;
    }
    propertyJac_realloc(&densityJac, num_dep);
  }

  num_species = matrl->Num_Species;
  num_species_eqn = matrl->Num_Species_Eqn;

  /*
   *   Branch according to the density model from the material properties
   *   structure
   */
  if (matrl->DensityModel == CONSTANT) {
    rho = matrl->density;

  } else if (matrl->DensityModel == USER) {
    (void)usr_density(matrl->u_density);
    rho = matrl->density;

  } else if (matrl->DensityModel == DENSITY_THERMEXP) {
    rho = matrl->u_density[0] /
          (1. + matrl->u_density[1] * (stateVector[TEMPERATURE] - matrl->reference[TEMPERATURE]));

  } else if (matrl->DensityModel == FILL) {
    F = stateVector[FILL];
    if (F < 0.1) {
      rho = matrl->u_density[1];
    } else {
      rho = matrl->u_density[0];
    }

  } else if (matrl->DensityModel == DENSITY_FOAM || matrl->DensityModel == DENSITY_FOAM_CONC ||
             matrl->DensityModel == DENSITY_FOAM_TIME_TEMP ||
             matrl->DensityModel == DENSITY_FOAM_TIME) {
    double param = matrl->u_density[0];
    double S = fv->c[0];

    rho = param * (1.0 - S);

    if (matrl->DensityModel == DENSITY_FOAM_TIME) {
      double rho_init, rho_final, aexp, time_delay, realtime;
      rho_init = matrl->u_density[0];   /* Initial density */
      rho_final = matrl->u_density[1];  /* final density */
      aexp = matrl->u_density[2];       /* Arhennius constant for time exponent*/
      time_delay = matrl->u_density[3]; /* time delay before foaming starts */
      if (time > time_delay) {
        realtime = time - time_delay;
        rho = rho_final + (rho_init - rho_final) * exp(-aexp * realtime);
      } else {
        rho = rho_init;
      }
    } else if (matrl->DensityModel == DENSITY_FOAM_TIME_TEMP) {
      double T;
      double rho_init = mp->u_density[0];   /* Initial density (units of gm cm-3 */
      double rho_final = mp->u_density[1];  /* Final density (units of gm cm-3 */
      double cexp = mp->u_density[2];       /* cexp in units of Kelvin sec */
      double coffset = mp->u_density[3];    /* offset in units of seconds */
      double time_delay = mp->u_density[4]; /* time delay in units of seconds */
      double drhoDT;
      T = stateVector[TEMPERATURE];
      if (time > time_delay) {
        double realtime = time - time_delay;
        double delRho = (rho_init - rho_final);
        double cdenom = cexp - coffset * T;
        double expT = exp(-realtime * T / cdenom);
        rho = rho_final + delRho * expT;
        drhoDT = delRho * expT * (-realtime * cexp / (cdenom * cdenom));
      } else {
        rho = rho_init;
        drhoDT = 0.0;
      }
      if (doJac) {
        propertyJac_addEnd(densityJac, TEMPERATURE, matID, 0, drhoDT, rho);
      }
    } else {
      GOMA_EH(GOMA_ERROR, "calc_density called");
    }

  } else if (matrl->DensityModel == DENSITY_FOAM_PMDI_10) {
    int var;
    int w;
    double volF = matrl->volumeFractionGas;

    double M_CO2 = matrl->u_density[0];
    double rho_liq = matrl->u_density[1];
    double ref_press = matrl->u_density[2];
    double Rgas_const = matrl->u_density[3];

    double rho_gas = (ref_press * M_CO2 / (Rgas_const * fv->T));

    rho = rho_gas * volF + rho_liq * (1 - volF);

    /* Now do sensitivies */

    var = MASS_FRACTION;
    if (volF > 0. && doJac) {
      if (pd->v[pg->imtrx][var]) {
        for (w = 0; w < pd->Num_Species; w++) {
          double drhodC = (rho_gas * mp->d_volumeFractionGas[MAX_VARIABLE_TYPES + w] -
                           rho_liq * mp->d_volumeFractionGas[MAX_VARIABLE_TYPES + w]);
          propertyJac_addEnd(densityJac, MASS_FRACTION, matID, w, drhodC, rho);
        }
      }
    }

    var = TEMPERATURE;
    if (volF > 0. && doJac) {
      if (pd->v[pg->imtrx][var]) {
        double drhoDT;
        drhoDT = (rho_gas / fv->T * volF + rho_gas * mp->d_volumeFractionGas[var] -
                  rho_liq * mp->d_volumeFractionGas[var]);
        propertyJac_addEnd(densityJac, TEMPERATURE, matID, 0, drhoDT, rho);
      }
    }

  } else if (matrl->DensityModel == DENSITY_MOMENT_BASED) {
    int var;
    int w;
    double volF = matrl->volumeFractionGas;

    double rho_gas = matrl->u_density[0];
    double rho_liq = matrl->u_density[1];

    rho = rho_gas * volF + rho_liq * (1 - volF);

    /* Now do sensitivies */

    var = MASS_FRACTION;
    if (volF > 0. && doJac) {
      if (pd->v[pg->imtrx][var]) {
        for (w = 0; w < pd->Num_Species; w++) {
          double drhodC = 0;
          propertyJac_addEnd(densityJac, MASS_FRACTION, matID, w, drhodC, rho);
        }
      }
    }

    var = TEMPERATURE;
    if (volF > 0. && doJac) {
      if (pd->v[pg->imtrx][var]) {
        double drhoDT;
        drhoDT = 0;
        propertyJac_addEnd(densityJac, TEMPERATURE, matID, 0, drhoDT, rho);
      }
    }

  } else if (matrl->DensityModel == LEVEL_SET) {
    double *param = matrl->u_density;
    double alpha = param[2];
    F = stateVector[FILL];
    if (F > alpha) {
      rho = param[1];
    } else if (F < -alpha) {
      rho = param[0];
    } else {
      rho = ((param[0] + param[1]) + (param[1] - param[0]) * sin(M_PIE * F / (2.0 * alpha))) / 2.0;
    }

  } else if (matrl->DensityModel == DENSITY_CONST_PHASE_FUNC) {
    double scratch_rho;

    scratch_rho = density(NULL, time);
    rho = scratch_rho;
  } else if (matrl->DensityModel == SUSPENSION) {
    species = (int)matrl->u_density[0];
    rho_f = matrl->u_density[1];
    rho_s = matrl->u_density[2];
    vol = stateVector[SPECIES_UNK_0 + species];
    if (vol < 0.) {
      rho = rho_f;
    } else {
      rho = rho_f + (rho_s - rho_f) * vol;
      if (doJac) {
        propertyJac_add1SpEnd(densityJac, matrl, species, (rho_s - rho_f), rho);
      }
    }

  } else if (matrl->DensityModel == SUSPENSION_PM) {
    /* MMH
     * For the SUSPENSION_PM model, the density is "usually" just the
     * fluid density.  It is only in two cases that the density we want is
     * actually the particle phase density.  1) In the particle phase
     * momentum equations. 2) In the particle phase advection equation
     * (like a chemical species).  In those two cases, this must be "undone".
     * This seemed the best way to do it so there was a reasonable default
     * behavior, and you didn't need a bunch of if()'s bracketing every
     * density() call...
     */

    rho = matrl->u_density[1];

    /* MMH
     * This doesn't account for other species present when calculating
     * density.  But they are currently just "collected" into one rho_f,
     * so all of the d_rho_dC's should be zero!.
     */

  } else if (matrl->DensityModel == THERMAL_BATTERY) {
    rho = matrl->u_density[0] - matrl->u_density[1] * (2.0 * stateVector[SPECIES_UNK_0]);
    if (doJac) {
      propertyJac_addEnd(densityJac, MASS_FRACTION, matID, 0, -matrl->u_density[1] * 2.0, rho);
    }

  } else if (matrl->DensityModel == DENSITY_IDEAL_GAS) {
    double mw_last, c_total;
    mw_last = mp->molecular_weight[pd->Num_Species_Eqn];
    /*
     * Specify the thermodynamic pressure here
     *   -> For the moment we ignore any coupling between
     *      the dynamic pressure field and the thermodynamic
     *      pressure field.
     */
    pressureThermo = upd->Pressure_Datum; /* This is always in cgs units, so not that useful*/

    /*
     * Specify the gas constant according to user's choice
     */
    RGAS_CONST = matrl->u_density[0]; /*(cgs:cm^3-atm/g-mole K)  (si_mm: mm^3-kPa/kg-mole-degK) */
    pressureThermo = matrl->u_density[1]; /*Pressure Datum in units of choice  */

    /*
     *  Base density calculation depends on species var type
     */
    switch (speciesVT) {
    case SPECIES_MASS_FRACTION:
      avgMolecWeight = wt_from_Yk(num_species, stateVector + SPECIES_UNK_0, mw);
      rho = pressureThermo * avgMolecWeight / (RGAS_CONST * stateVector[TEMPERATURE]);
      break;
    case SPECIES_MOLE_FRACTION:
      avgMolecWeight = wt_from_Xk(num_species, stateVector + SPECIES_UNK_0, mw);
      rho = pressureThermo * avgMolecWeight / (RGAS_CONST * stateVector[TEMPERATURE]);
      break;
    case SPECIES_CONCENTRATION:
      rho = 0.0;
      if (matrl->molecular_weight[num_species_eqn] < 0) {
        for (w = 0; w < num_species_eqn; w++) {
          rho += stateVector[SPECIES_UNK_0 + w] * mw[w];
        }
      } else {
        c_total = pressureThermo / (RGAS_CONST * stateVector[TEMPERATURE]);
        for (w = 0; w < num_species_eqn; w++) {
          rho += stateVector[SPECIES_UNK_0 + w] * (mw[w] - mw_last);
        }
        rho += c_total * mw_last;
      }
      break;
    case SPECIES_DENSITY:
      rho = 0.0;
      for (w = 0; w < num_species; w++) {
        rho += stateVector[SPECIES_UNK_0 + w];
      }
      break;
    default:
      fprintf(stderr, "Density error: species var type not handled: %d\n", speciesVT);
      exit(-1);
    }

    if (doJac) {
      if (pd->v[pg->imtrx][TEMPERATURE]) {
        propertyJac_addEnd(densityJac, TEMPERATURE, matID, 0, -rho / stateVector[TEMPERATURE], rho);
      }
      if (pd->v[pg->imtrx][MASS_FRACTION]) {
        switch (speciesVT) {
        case SPECIES_MASS_FRACTION:
          for (w = 0; w < num_species - 1; w++) {
            tmp = -rho * avgMolecWeight * (1.0 / mw[w] - 1.0 / mw[num_species - 1]);
            propertyJac_addEnd(densityJac, MASS_FRACTION, matID, w, tmp, rho);
          }
          propertyJac_addEnd(densityJac, MASS_FRACTION, matID, num_species - 1, 0.0, rho);
          break;
        case SPECIES_MOLE_FRACTION:
          for (w = 0; w < num_species - 1; w++) {
            tmp = rho / avgMolecWeight * (mw[w] - mw[num_species - 1]);
            propertyJac_addEnd(densityJac, MASS_FRACTION, matID, w, tmp, rho);
          }
          propertyJac_addEnd(densityJac, MASS_FRACTION, matID, num_species - 1, 0.0, rho);
          break;
        case SPECIES_CONCENTRATION:
          for (w = 0; w < num_species_eqn; w++) {
            tmp = mw[w] - mw_last;
            propertyJac_addEnd(densityJac, MASS_FRACTION, matID, w, tmp, rho);
          }
          /*propertyJac_addEnd(densityJac, MASS_FRACTION,
                             matID, num_species-1, 0.0, rho);  */
          break;
        default:
          fprintf(stderr, "Density error: species var type not handled: %d\n", speciesVT);
          exit(-1);
        }
      }
    }
  } else if (matrl->DensityModel == REACTIVE_FOAM)

  /* assume Amagat's law, no volume change upon
   * mixing. ACS
   */
  {
    double *param = matrl->u_density;
    double sv_p = param[0];
    double sum_sv = 0, drho_dT = 0., drho_dc = 0.;
    double *stateVector = matrl->StateVector;

    if (matrl->SpeciesSourceModel[0] == FOAM) {
      foam_species_source(matrl->u_species_source[0]);
    } else {
      GOMA_EH(GOMA_ERROR, "Must specify FOAM species source in the material's file");
    }

    for (w = 0; w < matrl->Num_Species; w++) {
      sum_sv += (matrl->specific_volume[w] - sv_p) * stateVector[SPECIES_UNK_0 + w];
    }
    sum_sv += sv_p;
    rho = 1. / sum_sv;

    drho_dc = 0.;
    drho_dT = 0.;

    if (doJac) {
      if (pd->v[pg->imtrx][TEMPERATURE]) {
        for (w = 0; w < matrl->Num_Species_Eqn; w++) {
          drho_dT -= (matrl->specific_volume[w] - sv_p) * rho * rho * 3. * matrl->species_source[w];
        }
        propertyJac_addEnd(densityJac, TEMPERATURE, matID, 0, drho_dT, rho);
      }
      if (pd->v[pg->imtrx][MASS_FRACTION]) {
        for (w = 0; w < matrl->Num_Species_Eqn; w++) {
          drho_dc = -(matrl->specific_volume[w] - sv_p) * rho * rho;
          propertyJac_addEnd(densityJac, MASS_FRACTION, matID, w, drho_dc, rho);
        }
      }
    }

  } else if (matrl->DensityModel == DENSITY_FOAM_PBE) {
    int species_BA_g;
    int species_CO2_g;
    int species_CO2_l;
    int species_BA_l;
    int err;
    int var;
    err = get_foam_pbe_indices(NULL, NULL, &species_BA_l, &species_BA_g, &species_CO2_l,
                               &species_CO2_g);
    if (err)
      return 0;

    double M_BA = mp->u_species_source[species_BA_l][0];
    double M_CO2 = mp->u_species_source[species_CO2_l][0];
    double rho_bubble = 0;
    double rho_foam = matrl->u_density[0];
    double ref_press = matrl->u_density[1];
    double Rgas_const = matrl->u_density[2];
    double d_rho_dT = 0;
    double d_rho_dM1 = 0;

    if (fv->c[species_BA_g] > PBE_FP_SMALL || fv->c[species_CO2_g] > PBE_FP_SMALL) {
      rho_bubble = (ref_press / (Rgas_const * fv->T)) *
                   (fv->c[species_CO2_g] * M_CO2 + fv->c[species_BA_g] * M_BA) /
                   (fv->c[species_CO2_g] + fv->c[species_BA_g]);
    }

    double inv_mom_frac = 1 / (1 + fv->moment[1]);
    rho = rho_bubble * (fv->moment[1] * inv_mom_frac) + rho_foam * inv_mom_frac;

    if (doJac) {
      var = TEMPERATURE;
      if (pd->v[pg->imtrx][var]) {
        d_rho_dT = (-rho) / fv->T;
        propertyJac_addEnd(densityJac, var, matID, 0, d_rho_dT, rho);
      }

      var = MOMENT1;
      if (pd->v[pg->imtrx][var]) {
        d_rho_dM1 = 2 * inv_mom_frac * inv_mom_frac;
        propertyJac_addEnd(densityJac, var, matID, 0, d_rho_dM1, rho);
      }
    }
  } else if (matrl->DensityModel == DENSITY_FOAM_PBE_EQN) {
    rho = fv->rho;

    if (doJac) {
      int var = DENSITY_EQN;
      if (pd->v[pg->imtrx][var]) {
        propertyJac_addEnd(densityJac, var, matID, 0, 1, rho);
      }
    }

  } else if (matrl->DensityModel == SOLVENT_POLYMER)

  /* assume Amagat's law, no volume change upon
   * mixing. ACS
   */
  {
    double *param = matrl->u_density;
    double sv_p = param[0];
    double sum_sv = 0, drho_dc = 0.;
    double *stateVector = matrl->StateVector;

    switch (speciesVT) {
    case SPECIES_MASS_FRACTION:
      for (w = 0; w < matrl->Num_Species; w++) {
        sum_sv += (matrl->specific_volume[w] - sv_p) * stateVector[SPECIES_UNK_0 + w];
      }
      sum_sv += sv_p;
      rho = 1. / sum_sv;

      if (doJac) {
        for (w = 0; w < matrl->Num_Species_Eqn; w++) {
          drho_dc = -(matrl->specific_volume[w] - sv_p) * rho * rho;
          propertyJac_addEnd(densityJac, MASS_FRACTION, matID, w, drho_dc, rho);
        }
      }
      break;
    case SPECIES_DENSITY:
      for (w = 0; w < matrl->Num_Species; w++) {
        rho += (1. - matrl->specific_volume[w] / sv_p) * fv->c[w];
      }
      rho += 1. / sv_p;

      if (doJac) {
        for (w = 0; w < matrl->Num_Species_Eqn; w++) {
          drho_dc = (1. - matrl->specific_volume[w] / sv_p);
          matrl->d_density[MAX_VARIABLE_TYPES + w] = drho_dc;
          propertyJac_addEnd(densityJac, MASS_FRACTION, matID, w, drho_dc, rho);
        }
      }
      break;
    case SPECIES_CONCENTRATION:
      for (w = 0; w < matrl->Num_Species; w++) {
        rho += (1. - matrl->specific_volume[w] / sv_p) * fv->c[w] * matrl->molecular_weight[w];
      }
      rho += 1. / sv_p;

      if (doJac) {
        for (w = 0; w < matrl->Num_Species_Eqn; w++) {
          drho_dc = matrl->molecular_weight[w] * (1. - matrl->specific_volume[w] / sv_p);
          matrl->d_density[MAX_VARIABLE_TYPES + w] = drho_dc;
          propertyJac_addEnd(densityJac, MASS_FRACTION, matID, w, drho_dc, rho);
        }
      }
      break;
    case SPECIES_MOLE_FRACTION:
      for (w = 0; w < matrl->Num_Species; w++) {
        sum_sv += (matrl->molar_volume[w] - sv_p) * stateVector[SPECIES_UNK_0 + w];
      }
      sum_sv += sv_p;
      rho = 1. / sum_sv;

      if (doJac) {
        for (w = 0; w < matrl->Num_Species_Eqn; w++) {
          drho_dc = -(matrl->molar_volume[w] - sv_p) * rho * rho;
          propertyJac_addEnd(densityJac, MASS_FRACTION, matID, w, drho_dc, rho);
        }
      }
      break;
    default:
      GOMA_WH(-1, "SOLVENT_POLYMER defaulting to MASS_FRACTION\n");
      for (w = 0; w < matrl->Num_Species; w++) {
        sum_sv += (matrl->specific_volume[w] - sv_p) * stateVector[SPECIES_UNK_0 + w];
      }
      sum_sv += sv_p;
      rho = 1. / sum_sv;

      if (doJac) {
        for (w = 0; w < matrl->Num_Species_Eqn; w++) {
          drho_dc = -(matrl->specific_volume[w] - sv_p) * rho * rho;
          propertyJac_addEnd(densityJac, MASS_FRACTION, matID, w, drho_dc, rho);
        }
      }
      break;
    }
  } else if (matrl->DensityModel == DENSITY_CONSTANT_LAST_CONC) {
    /*
     *
     */
    if (matrl->Species_Var_Type != SPECIES_CONCENTRATION) {
      GOMA_EH(GOMA_ERROR, "unimplemented");
    }
    if (pd->Num_Species > pd->Num_Species_Eqn) {
      w = matrl->Num_Species_Eqn;
      rho = mw[w] * fv->c[w];
    } else {
      rho = 0.0;
    }
    for (w = 0; w < matrl->Num_Species_Eqn; w++) {
      rho += mw[w] * stateVector[SPECIES_UNK_0 + w];
    }
    if (doJac) {
      for (w = 0; w < matrl->Num_Species_Eqn; w++) {
        propertyJac_addEnd(densityJac, MASS_FRACTION, matID, w, mw[w], rho);
      }
    }

  } else {
    GOMA_EH(GOMA_ERROR, "Unrecognized density model");
  }
  return (rho);
}
/**************************************************************************/
/**************************************************************************/
/**************************************************************************/

double calc_concentration(MATRL_PROP_STRUCT *matrl, int doJac, PROPERTYJAC_STRUCT *concJac)

/**********************************************************************
 *
 * calc_concentration():
 *
 *   Calculate the mixture concentration and its derivatives using
 *   the properties evaluated from the local materials state vector.
 *   Note, we strive
 *   to have a complete description of the materials models in this
 *   function.
 *
 * Output
 * -----
 *    If doJac is true, this routine also fills up the concJac
 *    structure with derivative information.
 *
 *
 *  Return
 * --------
 *    Actual value of the concentration
 *************************************************************************/
{
  int k;
  double C_mix = 0.0;
  double *C_base;
  if (matrl->DensityModel == DENSITY_CONSTANT_LAST_CONC) {
    if (matrl->StateVector_speciesVT != SPECIES_CONCENTRATION) {
      GOMA_EH(GOMA_ERROR, "unimplemented");
    }
    C_base = matrl->StateVector + SPECIES_UNK_0;
    for (k = 0; k < matrl->Num_Species_Eqn; k++) {
      C_mix += C_base[k];
    }
    C_mix += matrl->u_density[0];
  } else if (matrl->DensityModel == DENSITY_CONSTANT_CONC) {

  } else if (matrl->DensityModel == CONSTANT) {

  } else {
    GOMA_EH(GOMA_ERROR, "unimplemented");
  }

  return C_mix;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void load_properties(MATRL_PROP_STRUCT *matrl, double time)

/*************************************************************************
 *
 *  load_properties():
 *
 *  Commonly used properties are calculated here based on the state
 *  variable in the material structure.
 *************************************************************************/
{
  if (af->Assemble_Jacobian) {
    if (!matrl->DensityJac) {
      matrl->DensityJac = alloc_struct_1(PROPERTYJAC_STRUCT, 1);
    }
    matrl->density = calc_density(matrl, TRUE, matrl->DensityJac, time);
  } else {
    matrl->density = calc_density(matrl, FALSE, NULL, time);
  }
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
