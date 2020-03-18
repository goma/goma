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
 *$Id: mm_species.c,v 5.2 2007-09-18 18:53:46 prschun Exp $
 */

/*************** R O U T I N E S   I N   T H E   F I L E *********************
 *
 *    NAME				TYPE		CALLED_BY
 *--------------------------------------------------------------------
 *
 *    normalize_species_fractions         int
 *    check_consistent_fraction_vector    int
 *    wt_from_Ck                          void
 *    wt_from_Xk                          void
 *    wt_from_Yk                          void
 *    Yk_from_Xk                          void
 *    Xk_from_Yk                          void
 *    Ck_from_Xk                          void
 *    Xk_from_Ck                          void
 *    Ck_from_Dk                          void
 *    Dk_from_Ck                          void
 *    convert_species_var                 void
 ************************************************************************** **/


#include <stdio.h>
#include <string.h>
#include <math.h>

#include "std.h"
#include "mm_species.h"
#include "mm_as.h"
#include "mm_as_structs.h"
#include "mm_eh.h"
#include "mm_mp.h"
#include "mm_mp_const.h"
#include "mm_mp_structs.h"
#include "rf_fem_const.h"


/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int  normalize_species_fractions(double frac_vec[], const int num_species)

    /*************************************************************************
     *
     * normalize_species_fractions:
     *
     *   This utility routine normalizes a set of mass or mole fraction
     *   vectors so that they sum exactly to one.
     *
     * Input
     * --------
     *  frac_vec[] = Vector of input mole or mass fractions
     *  num_species = Number of species (and length of the frac_vec[] vector)
     *
     * Output
     * -------
     *  frac_vec[] = Normalized values of the mass fractions
     *
     * Return
     * -------
     *  1 = Successful completion of task
     * -1 = Input mole fractions were amiss. The fractions were reassigned to
     *      all be equal and to add to one.
     *
     * IMPLEMENTATION NOTE
     *      This routine can be used to implement cropping schemes in a
     *      unified manor.
     *************************************************************************/
{
  int i;
  double sum = 0.0;
  if (num_species <= 0) return 1;
  for (i = 0; i < num_species; i++) {
     sum += frac_vec[i];
  }
  if (sum <= 0.0) {
    for (i = 0; i < num_species; i++) {
      frac_vec[i] = 1.0 / num_species;
    }
    return -1;
  }
  sum = 1.0 / sum;
  for (i = 0; i < num_species; i++) {
      frac_vec[i] *= sum;
  }
  return 1;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int check_consistent_fraction_vector(struct Variable_Initialization *var_init,
				     int num_entries, int num_species,
				     double frac_vec[])

    /*************************************************************************
     *
     * check_consistent_fraction_vector:
     *
     *   This routine processes the initialization information for species
     *   and processes the entries for mass, mole, and volume fraction in
     *   a consistent manner. A check that the total sum of fractions is
     *   equal to one is made. The final normalized fractions are returned
     *   in frac_vec[]. The function returns the type of fractions that
     *   frac_vec[] represents.
     *
     * Input
     * ---------
     *  var_init->var   int representing the variable ID or fraction type
     *                of entry i
     *  var_init->ktype  int representing the species number for entry i
     *  var_init->init_val   double representing the value for entry i
     *  num_entries = Vector lenght of structures pointed to be
     *                Variable_Initialization
     *
     * Output
     * -------
     *
     * frac_vec[k] : On return this is either all zeroes, when no fractions
     *              are found in type_vec[] or it is the normalized
     *              fraction vectors for the species in the problem
     *
     * Return
     * ---------
     *  Value = 0: no species fractions found in the list
     *          Pos: Indicates the indentitiy of the fractions.
     *               The identity of the fraction is specified in
     *               rf_fem_const.h
     *************************************************************************/
{
  int i, frac_type = 0, kspec;
  double sum = 1.0;

  /*
   * Zero out the return vector
   */
  for (i = 0; i < num_species; i++) {
    frac_vec[i] = 0.0;
  }

  /*
   * Parse the initialization data looking for one type of
   * specification of fraction and only one type. Store the values found
   * in the correct slots in the return vector, frac_vec[].
   */
  for (i = 0; i < num_entries; i++) {
    switch (var_init->var) {
    case SPECIES_MASS_FRACTION :
	kspec = var_init->ktype;
	if (kspec < 0 || kspec >= num_species) goto bad_ktype;
	frac_vec[kspec] = var_init->init_val;
	if (frac_type > 0 && frac_type != var_init->var) goto bad_data;
	frac_type = SPECIES_MASS_FRACTION;	
	break;
    case SPECIES_MOLE_FRACTION :
	kspec = var_init->ktype;
	if (kspec < 0 || kspec >= num_species) goto bad_ktype;
	frac_vec[kspec] = var_init->init_val;
        if (frac_type > 0 && frac_type != var_init->var) goto bad_data;
	frac_type = SPECIES_MOLE_FRACTION;
	break;
    case SPECIES_VOL_FRACTION :
	kspec = var_init->ktype;
	if (kspec < 0 || kspec >= num_species) goto bad_ktype;
	frac_vec[kspec] = var_init->init_val;
        if (frac_type > 0 && frac_type != var_init->var) goto bad_data;
	frac_type = SPECIES_VOL_FRACTION;
	break;
    case MASS_FRACTION :
        if (frac_type > 0 && frac_type != var_init->var) goto bad_data;
	break;
    default :
	break;
    }
    var_init++;
  }

  /*
   * Calculate the sum of "mass" fraction equals one condition
   */
  for (i = 0; i < num_species; i++) {
    sum -= frac_vec[i];
  }
  /*
   * Decide whether the input is good enough to continue,
   * whether it needs to be fixed up or whether everything is OK.
   */
  if (fabs(sum) > DBL_SMALL) {
    if (fabs(sum) < 0.01) {
      (void) normalize_species_fractions(frac_vec, num_species);
    } else {
      printf("check_consistent_fraction_vector ERROR: species fraction vector");
      printf(" is very far from summing to one: %g\n", 1.0 - sum);
      EH(-1, "check_consistent_fraction_vector\n");
    }
  }

  /*
   * Return the type of fraction found
   */
  return frac_type;
  /*
   * Error exits
   */
 bad_data:;
  printf("check_consistent_fraction_vector ERROR!\n");
  printf("\t Mismatched species initialization types: %d %d\n", frac_type,
	 var_init->var);
  EH(-1, "check_consistent_fraction_vector\n");
  return -1;
 bad_ktype:;
  printf("check_consistent_fraction_vector ERROR!\n");
  printf("Species nunumber is out of bounds: %d, %d", kspec, num_species);
  EH(-1, "check_consistent_fraction_vector\n");
  return -1;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

double wt_from_Ck(const int num_species, const double *Ck,
		  const double *MolecWeight)

    /************************************************************************
     *
     * wt_from_Ck:
     *
     *   This small utility function calculates the average molecular
     *   weight of a phase given the concentrations and molecular weight
     *   vectors of all of the species in a phase. A continguous species
     *   vector is assumed.
     ************************************************************************/
{
  int i;
  double sum = 0.0;
  double sumc = 0.0;
  for (i = 0; i < num_species; i++) {
    sum += Ck[i] * MolecWeight[i];
    sumc += Ck[i];
  }
  if (sumc > 0.0) return sum / sumc;
  return 0.0;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

double wt_from_Xk(const int num_species, const double *Xk,
		  const double *MolecWeight)

    /************************************************************************
     *
     * wt_from_Xk:
     *
     *   This small utility function calculates the average molecular
     *   weight of a phase given the mole fraction and molecular weight
     *   vectors of all of the species in a phase. A continguous species
     *   vector is assumed.
     ************************************************************************/
{
  int i;
  double sum = 0.0;
  for (i = 0; i < num_species; i++) {
     sum += Xk[i] * MolecWeight[i];
  }
  return sum;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

double 
wt_from_Yk(const int num_species,  const double *Yk,
	   const double *molecWeight)

    /************************************************************************
     *
     * wt_from_Yk():
     *
     *   This small utility function calculates the average molecular
     *   weight of a phase given the mass fraction and molecular weight
     *   vectors of all of the species in a phase.  A continguous species
     *   vector is assumed.
     ************************************************************************/
{
  int i;
  double sum = 0.0;
  for (i = 0; i < num_species; i++) {
    sum += Yk[i] / molecWeight[i];
  }
  if (fabs(sum) < DBL_SMALL) return 0.0;    
  return 1.0 / sum;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void Yk_from_Xk(const int num_species, double *Yk, double *Xk,
		const double *MolecWeight)

    /************************************************************************
     *
     * Yk_from_Xk:
     *
     *   This small utility function calculates the mass fraction vector
     *   of a given phase, given the molefraction vector and the
     *   species molecular weights.  A continguous species
     *   vector is assumed.
     ************************************************************************/
{
  int i;
  double wt = wt_from_Xk(num_species, Xk, MolecWeight);
  if (fabs(wt) < DBL_SMALL) {
    for (i = 0; i < num_species; i++) Yk[i] = Xk[i];
  } else {
    for (i = 0; i < num_species; i++) {
      Yk[i] = Xk[i] * MolecWeight[i] / wt;
    }
  }
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void Xk_from_Yk(const int num_species, double *Xk, double *Yk,
		const double *MolecWeight)

    /************************************************************************
     *
     * Xk_from_Yk:
     *
     *   This small utility function calculates the mole fraction vector
     *   of a given phase, given the mass fraction vector and the
     *   species molecular weights.  A continguous species
     *   vector is assumed.
     ************************************************************************/
{
  int i;
  double wt = wt_from_Yk(num_species, Yk, MolecWeight);  
  if (fabs(wt) < DBL_SMALL) {
    for (i = 0; i < num_species; i++) Xk[i] = Yk[i];
  } else {
    for (i = 0; i < num_species; i++) {
      if (MolecWeight[i] < DBL_SMALL) {
        Xk[i] = Yk[i];
      } else {
        Xk[i] = Yk[i] * wt / MolecWeight[i];
      }
    }
  }
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void Ck_from_Xk(const int num_species, double *Ck, double *Xk,
		MATRL_PROP_STRUCT *matrl, const double time)

    /************************************************************************
     *
     * Ck_from_Xk:
     *
     *   This small utility function calculates the species concentration
     *   vector from the species mole fraction vector.
     *
     *   Properties, such as the temperature and the thermodynamic pressure,
     *   are assumed to be already loaded into the State Vector Field in the 
     *   materials property structure. These properties are passed down
     *   to the calc_density, which calculates the density of the material.
     ************************************************************************/
{
  int i;
  double rho, wt = wt_from_Xk(num_species, Xk, matrl->molecular_weight);
  double *sv = matrl->StateVector;
  if (sv + SPECIES_UNK_0 != Xk) {
    for (i = 0; i < num_species; i++) {
      sv[SPECIES_UNK_0 + i] = Xk[i];
    }
    matrl->StateVector_speciesVT = SPECIES_MOLE_FRACTION;
  }
  rho = calc_density(matrl, FALSE, NULL, time);
  rho /= wt;
  for (i = 0; i < num_species; i++) {
    Ck[i] = Xk[i] * rho;
  }
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void Xk_from_Ck(const int num_species, double *Xk, double *Ck)

    /************************************************************************
     *
     * Xk_from_Ck:
     *
     *   This small utility function calculates the species mole
     *   fraction vector from the species concentration vector.
     *
     ************************************************************************/
{
  int i;
  double Ctot = 0.0;
  for (i = 0; i < num_species; i++) Ctot += Ck[i];
  if (Ctot > 0.0) {
    for (i = 0; i < num_species; i++) Xk[i] = Ck[i] / Ctot;
  } else {
    for (i = 0; i < num_species; i++) Xk[i] = 1.0 / num_species;
  }
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void Yk_from_Ck(const int num_species, double *Yk, double *Ck, 
		MATRL_PROP_STRUCT *matrl)

    /************************************************************************
     *
     * Yk_from_Ck:
     *
     *   This small utility function calculates the species mass
     *   fraction vector from the species concentration vector.
     ************************************************************************/
{
  int i;
  double Ctot = 0.0;
  double *mw = matrl->molecular_weight;
  for (i = 0; i < num_species; i++) Ctot += Ck[i] * mw[i];
  if (Ctot > 0.0) {
    for (i = 0; i < num_species; i++) Yk[i] = mw[i] * Ck[i] / Ctot;
  } else {
    for (i = 0; i < num_species; i++) Yk[i] = 1.0 / num_species;
  }
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void Ck_from_Yk(const int num_species, double *Ck, double *Yk,
		MATRL_PROP_STRUCT *matrl, const double time)

    /************************************************************************
     *
     * Ck_from_Yk:
     *
     *   This small utility function calculates the species concentration
     *   vector from the species mole fraction vector.
     *
     *   Properties, such as the temperature and the thermodynamic pressure,
     *   are assumed to be already loaded into the State Vector Field in the 
     *   materials property structure. These properties are passed down
     *   to the calc_density, which calculates the density of the material.
     ************************************************************************/
{
  int i;
  double rho;
  double *mw = matrl->molecular_weight;
  double *sv = matrl->StateVector;
  if (sv + SPECIES_UNK_0 != Yk) {
    for (i = 0; i < num_species; i++) {
      sv[SPECIES_UNK_0 + i] = Yk[i];
    }
    matrl->StateVector_speciesVT = SPECIES_MASS_FRACTION;
  }
  rho = calc_density(matrl, FALSE, NULL, time);
  for (i = 0; i < num_species; i++) {
    Ck[i] = Yk[i] * rho / mw[i];
  }
}
/*****************************************************************************/
/*****************************************************************************/

void Dk_from_Ck(const int num_species, double *Dk, double *Ck, 
		MATRL_PROP_STRUCT *matrl)

    /************************************************************************
     *
     * Dk_from_Ck:
     *
     *   This small utility function calculates the species density
     *   vector from the species concentration vector.
     ************************************************************************/
{
  int i;
  double *mw = matrl->molecular_weight;
  for (i = 0; i < num_species; i++) Dk[i] = mw[i] * Ck[i];
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void Ck_from_Dk(const int num_species, double *Ck, double *Dk,
		MATRL_PROP_STRUCT *matrl)

    /************************************************************************
     *
     * Ck_from_Dk:
     *
     *   This small utility function calculates the species concentration
     *   vector from the species density vector.
     *
     *   Properties, such as the temperature and the thermodynamic pressure,
     *   are assumed to be already loaded into the State Vector Field in the 
     *   materials property structure. These properties are passed down
     *   to the calc_density, which calculates the density of the material.
     ************************************************************************/
{
  int i;
  double *mw = matrl->molecular_weight;
  double *sv = matrl->StateVector;
  if (sv + SPECIES_UNK_0 != Dk) {
    for (i = 0; i < num_species; i++) {
      sv[SPECIES_UNK_0 + i] = Dk[i];
    }
    matrl->StateVector_speciesVT = SPECIES_DENSITY;
  }
  for (i = 0; i < num_species; i++) {
    Ck[i] = Dk[i] / mw[i];
  }
}
/*****************************************************************************/
/*****************************************************************************/

int convert_species_var(int species_Var_Type, 
			struct Material_Properties *mp_local,
		        int frac_vec_speciesVT, double *frac_vec,
			double time)

    /************************************************************************
     *
     * convert_species_var:
     *
     *   This converts the various species unknown types amongst themselves.
     *
     * Input
     * -----
     * species_Var_Type: Species variable type desired
     * mat_num         : Material number of the current material
     * mp_local        : Pointer to the material properties structure
     * frac_vec_speciesVT : Current form of the species variable type
     *                 for the vector frac_vec
     * frac_vec        : Vector of species unknowns: Form is specified
     *                   by frac_vec_speciesVT
     * -> Any other properties that this function needs is obtained
     *    from the state variable contained in the materials property
     *    structure, mp_local.
     *
     * Output
     * --------
     * frac_vec        : Output value for the species fraction vector in the
     *                   form specified by species_Var_Type.
     *
     * Return
     * -------
     * -1              : Case is currently not handled.
     *  1              : successful return
     *
     *   NOTE:
     *      1) Needs to be upgraded to include for the possibility that there
     *         may be more than one phase in each material domain.
     *      2) Needs to be upgraded to include the possibility that a
     *         particular species equation doesn't contribute to the sum
     *         MF = 1 condition (i.e., conduction electrons)
     ************************************************************************/
{
  int retn = species_Var_Type;
  double *molecWeight = mp_local->molecular_weight;
  int num_species = mp_local->Num_Species;

  if (species_Var_Type == SPECIES_MASS_FRACTION) {
    switch (frac_vec_speciesVT) {
    case SPECIES_MOLE_FRACTION:
	Yk_from_Xk(num_species, frac_vec, frac_vec, molecWeight);
	break;
    case SPECIES_MASS_FRACTION:
	break;
    case SPECIES_CONCENTRATION:
	Yk_from_Ck(num_species, frac_vec, frac_vec, mp_local);
	break;
    default:
	retn = -1;
	printf("Case not covered\n");
	break;
    }
  } else if (species_Var_Type == SPECIES_MOLE_FRACTION) {
   switch (frac_vec_speciesVT) {
    case SPECIES_MASS_FRACTION:
	Xk_from_Yk(num_species, frac_vec, frac_vec, molecWeight);
	break;
    case SPECIES_MOLE_FRACTION:
	break;
   case SPECIES_CONCENTRATION:
        Xk_from_Ck(num_species, frac_vec, frac_vec);
        break;
    default:
	retn = -1;
	printf("Case not covered\n");
	break;
    }
  } else if (species_Var_Type == SPECIES_CONCENTRATION) {
   switch (frac_vec_speciesVT) {
    case SPECIES_MASS_FRACTION:
	Ck_from_Yk(num_species, frac_vec, frac_vec, mp_local, time);
	break;
    case SPECIES_MOLE_FRACTION:
	Ck_from_Xk(num_species, frac_vec, frac_vec, mp_local, time);
	break;
    case SPECIES_CONCENTRATION:
	break;
    case SPECIES_DENSITY:
	Ck_from_Dk(num_species, frac_vec, frac_vec, mp_local);
	break;
    default:
	retn = -1;
	printf("Case not covered\n");
	break;
    }
  } else if (species_Var_Type == SPECIES_DENSITY) {
   switch (frac_vec_speciesVT) {
    case SPECIES_CONCENTRATION:
	Dk_from_Ck(num_species, frac_vec, frac_vec, mp_local);
	break;
    case SPECIES_DENSITY:
	break;
    default:
	retn = -1;
	printf("Case not covered\n");
	break;
    }

  } else {
    retn = -1;
  }
  return retn;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void
deriv1_Ck_to_Yk(double *deriv_wrt_species, MATRL_PROP_STRUCT *mp_local,
		double *Ck, const double time)

    /************************************************************************
     *
     * deriv1_Ck_to_Yk():
     *
     *   This routine changes a derivative wrt species concentration to a
     *   derivative wrt mass fraction, for a single vector of dependencies
     *   to species concentrations.
     *
     *  Input
     * --------
     *     deriv_wrt_species[k] = Derivative wrt species concentration
     *     mp_local = pointer to the materials property
     *     Ck[k] = Current value of the species concentration
     *
     *  We use the formula:
     *     
     *     d c_i        1                       C_i * M_i   d rho
     *     ------ =    --  ( delta_i_j * rho -  --------- * ------ )
     *     d Y_j       M_i                         rho      d Y_k
     *
     *  To create the resulting formula:
     *     
     *     d S                d S     rho * delta_i_j      C_i       d rho
     *     ------ = Sum_1toN[ ----- ( --------------- -  --------- * ------ )]
     *     d Y_j       i      d C_i         M_i            rho       d Y_j
     ************************************************************************/
{
  int i;
  int num_species = mp_local->Num_Species, index;
  double tmp, rho, *mp_ptr, *sp_ptr;
  /*
   *  It would be prudent to move the property jac structure to the mp
   *  structure, where it doesn't have to be malloced.
   */
  PROPERTYJAC_STRUCT *densityJac = NULL;
  propertyJac_realloc(&densityJac, num_species+1);
  rho = calc_density(mp_local, TRUE, densityJac, time);
  mp_ptr = mp_local->molecular_weight;
  for (i = 0, tmp = 0.0; i < num_species; i++) {
    tmp += deriv_wrt_species[i] * Ck[i];
    deriv_wrt_species[i] *= rho / mp_ptr[i];
  }
  tmp /= rho;
  index = propertyJac_find_species_unk(densityJac);
  if (index >= 0) { 
    sp_ptr = densityJac->JacVector + index;
    for (i = 0; i < num_species; i++) {
      deriv_wrt_species[i] -= tmp * sp_ptr[i];
    }
  }
  propertyJac_destroy(&densityJac);
} 
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void
assign_species_var_type(const int mn, const int speciesVT, const int optional)

    /***********************************************************************
     *
     *  assign_species_var_Type
     *
     ***********************************************************************/
{
  if (mp_glob[mn]->Species_Var_Type != SPECIES_UNDEFINED_FORM &&
      mp_glob[mn]->Species_Var_Type != speciesVT &&
      mp_glob[mn]->Species_Var_Type != 0) {
    if (!optional) {
      sprintf(Err_Msg, "assign_species_var_type ERROR, mat %d: species var type "
	      " assigned to %d while it was already assigned to %d",
	      mn, speciesVT, mp_glob[mn]->Species_Var_Type);
      EH(-1, Err_Msg);
    } else {
    printf("assign_species_var_type WARNING, mat %d: species var type "
	   "failed opt assigned to %d while it was already assigned to %d",
	   mn, speciesVT, mp_glob[mn]->Species_Var_Type);
    }
  } else {
    mp_glob[mn]->Species_Var_Type = speciesVT;
    pd_glob[mn]->Species_Var_Type = speciesVT;
  }
  if (upd->Species_Var_Type != SPECIES_UNDEFINED_FORM &&
      upd->Species_Var_Type != speciesVT &&
      upd->Species_Var_Type != 0) {
    if (!optional) {
      sprintf(Err_Msg, "assign_species_var_type ERROR, upd: species var type "
	      " assigned to %d while it was already assigned to %d",
	      speciesVT, mp_glob[mn]->Species_Var_Type);
      EH(-1, Err_Msg);
    } else {
    printf("assign_species_var_type WARNING, upd: species var type "
	   "failed opt assigned to %d while it was already assigned to %d",
	   speciesVT, mp_glob[mn]->Species_Var_Type);
    }
  } else {
    upd->Species_Var_Type = speciesVT;
  }
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void
assign_global_species_var_type(const int speciesVT, const int optional)

    /***********************************************************************
     *
     *  assign_global_species_var_Type
     *
     ***********************************************************************/
{
  int mn;
  if (upd->Species_Var_Type != SPECIES_UNDEFINED_FORM &&
      upd->Species_Var_Type != speciesVT &&
      upd->Species_Var_Type != 0) {
    if (!optional) {
      sprintf(Err_Msg, "assign_species_var_type ERROR, upd: species var type "
	      " assigned to %d while it was already assigned to %d",
	      speciesVT, upd->Species_Var_Type);
      EH(-1, Err_Msg);
    } else {
    printf("assign_species_var_type WARNING, upd: species var type "
	   "failed opt assigned to %d while it was already assigned to %d",
	   speciesVT, upd->Species_Var_Type);
    }
  } else {
    upd->Species_Var_Type = speciesVT;
    for (mn = 0; mn < MAX_NUMBER_MATLS; mn++) {
      assign_species_var_type(mn, speciesVT, optional);
    }
  }
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void assign_species_prefix(const int species_var_name, char * retn_string)

   /**************************************************************************
    *
    * assign_species_prefix:
    *
    *  This function assigns the one to three letter prefix for the species
    *  unknown depending upon the species type.
    *
    *
    *  Input
    * --------
    *  species_var_name: Valid species type
    *                   (The valid species types are listed in rf_fem_const.h)
    * Output
    * --------
    *  retn_string 
    **************************************************************************/
{
  if (retn_string == NULL) {
   EH(-1, "assign_species_prefix: bad interface\n");
  }
  switch (species_var_name) {
  case SPECIES_MASS_FRACTION:
      (void) strcpy(retn_string, "YK_");
      break;
  case SPECIES_MOLE_FRACTION:
      (void) strcpy(retn_string, "XK_");
      break;
  case SPECIES_VOL_FRACTION:
      (void) strcpy(retn_string, "VK_");
      break;
  case SPECIES_DENSITY:
      (void) strcpy(retn_string, "DK_");
      break;
  case SPECIES_CONCENTRATION:
      (void) strcpy(retn_string, "CK_");
      break;
  case SPECIES_CAP_PRESSURE:
      (void) strcpy(retn_string, "PK_");
      break;
  case SPECIES_UNDEFINED_FORM:
      (void) strcpy(retn_string, "Y");
      break;
  default:
      (void) strcpy(retn_string, "Y");
      printf("assign_species_prefix: WARNING unknown form of species vars: %d\n",
	     species_var_name);
  } 
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
