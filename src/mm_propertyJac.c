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
 *$Id: mm_propertyJac.c,v 5.1 2007-09-18 18:53:45 prschun Exp $
 */

/*************** R O U T I N E S   I N   T H E   F I L E ***********************
 *
 *    NAME				TYPE		CALLED_BY
 *--------------------------------------------------------------------
 *
 *    chemkin_mat_prop_init   
 ******************************************************************************/


#include <stdio.h>
#include <math.h>

#include "std.h"
#include "rf_fem_const.h"
#include "mm_eh.h"
#include "rf_vars_const.h"
#include "mm_as_structs.h"
#include "mm_mp_structs.h"
#include "rf_allo.h"

/**************************************************************************/
/**************************************************************************/
/**************************************************************************/

void
propertyJac_realloc(PROPERTYJAC_STRUCT **jac_ptr, int num)

    /**********************************************************************
     *
     * propertyJac_realloc()
     *
     *  Reallocates (or potentially allocates for the very first time) 
     *  a PropertyJac Structure to a different size. All information
     *  is preserved.
     *
     **********************************************************************/
{
  PROPERTYJAC_STRUCT *jac;
  if (jac_ptr == NULL) {
    EH(GOMA_ERROR,"Interface error");
  }
  jac = *jac_ptr;
  if (jac == NULL) {
    jac = alloc_struct_1(PROPERTYJAC_STRUCT, 1);
    *jac_ptr = jac;
    jac->Species_Type = SPECIES_UNDEFINED_FORM;
  }
  if (num > jac->NUM_TERMS_MALLOC) {
    realloc_ptr_1((void ***) &(jac->Var_List), num, jac->NUM_TERMS_MALLOC);
    realloc_int_1(&(jac->idof),      num, jac->NUM_TERMS_MALLOC);
    realloc_int_1(&(jac->Var_Type),  num, jac->NUM_TERMS_MALLOC);
    realloc_int_1(&(jac->MatID),     num, jac->NUM_TERMS_MALLOC);
    realloc_dbl_1(&(jac->Var_Value), num, jac->NUM_TERMS_MALLOC);

    realloc_dbl_1(&(jac->JacVector), num, jac->NUM_TERMS_MALLOC);
  }
  jac->NUM_TERMS_MALLOC = num;
}
/**************************************************************************/
/**************************************************************************/
/**************************************************************************/

void
propertyJac_free(PROPERTYJAC_STRUCT *jac)

    /**********************************************************************
     *
     * propertyJac_free()
     *
     *  Frees the underlying memory in the a PropertyJac structure
     **********************************************************************/
{
  if (jac == NULL) return;
  if (jac->NUM_TERMS_MALLOC > 0) {
    safer_free((void **) &(jac->Var_List));
    safer_free((void **) &(jac->idof));
    safer_free((void **) &(jac->Var_Type));
    safer_free((void **) &(jac->Var_Value));
    safer_free((void **) &(jac->MatID));
    safer_free((void **) &(jac->JacVector));
  }
}
/**************************************************************************/
/**************************************************************************/
/**************************************************************************/

void
propertyJac_destroy(PROPERTYJAC_STRUCT **jac_ptr)

    /**********************************************************************
     *
     * propertyJac_destroy()
     *
     *  Frees the underlying memory in the a PropertyJac structure, and
     *  then frees the structure itself.
     **********************************************************************/
{
  propertyJac_free(*jac_ptr);
  safer_free((void **)jac_ptr);
}
/**************************************************************************/
/**************************************************************************/
/**************************************************************************/

void
propertyJac_addEnd(PROPERTYJAC_STRUCT *jac, int varType, int matID,
		   int isubvar, double jacEntry, double propValue)

    /**********************************************************************
     *
     * propertyJac_addentryEnd()
     *
     *  Adds an entry to the propertyjac structure. Note: this function
     *  is complete and slow. This is a candidate for inlining.
     **********************************************************************/
{
  int iTerm;
  VARIABLE_DESCRIPTION_STRUCT *vd;
  iTerm = jac->Num_Terms;
  if (iTerm >= jac->NUM_TERMS_MALLOC) {
    propertyJac_realloc(&jac, iTerm + 1);
  }

  vd = get_vd_ptr(varType, matID, isubvar);
  if (vd == NULL) {
    EH(GOMA_ERROR,"Error can't find variable description");
  }
  jac->Var_List[iTerm] = vd;
  jac->idof[iTerm] = isubvar; 
  jac->Var_Type[iTerm] = varType;
  jac->MatID[iTerm] = matID;
  jac->JacVector[iTerm] = jacEntry;
  if (iTerm > 0) {
    if (fabs(propValue - jac->Property_Value) > 1.0E-5) {
      EH(GOMA_ERROR,"Incompatible propertyValue search");
    }    
  }
  jac->Property_Value = propValue;
  jac->Num_Terms++;
}
/**************************************************************************/
/**************************************************************************/
/**************************************************************************/

void
propertyJac_add1SpEnd(PROPERTYJAC_STRUCT *jac,
		      MATRL_PROP_STRUCT *mat_ptr,
		      int isubvar, double jacEntry, double propValue)

    /**********************************************************************
     *
     * propertyJac_add1SpEnd()
     *
     *  Adds a vector of species entries for the current property. Only
     *  one of the species entries, isubvar, is nonzero. This satisfies
     *  the constraint that species be included as a contiguous complete
     *  vector dependency.
     **********************************************************************/
{
  int iTerm, i;
  int num_species = mat_ptr->Num_Species;
  int matID = mat_ptr->MatID;
  VARIABLE_DESCRIPTION_STRUCT *vd;
  iTerm = jac->Num_Terms;
  /*
   * Check to make sure the species havn't been added anywhere
   * else
   */
  i = propertyJac_find_species_unk(jac);
  if (i != -1) {
    EH(GOMA_ERROR,"Tried to add species twice");
  }
  if (iTerm + num_species > jac->NUM_TERMS_MALLOC) {
    propertyJac_realloc(&jac, iTerm + num_species);
  }
  for (i = 0; i < num_species; i++) {
    vd = get_vd_ptr(SPECIES_UNK_0 + i, matID, 0);
    if (vd == NULL) {
      vd = get_vd_ptr(MASS_FRACTION, matID, i);
      if (vd == NULL) {
	EH(GOMA_ERROR, "Error");
      }
    }
    jac->Var_List[iTerm + i] = vd;
    jac->idof[iTerm + i] = i; 
    jac->Var_Type[iTerm + i] = SPECIES_UNK_0;
    jac->MatID[iTerm + i] = matID;
    jac->JacVector[iTerm + i] = jacEntry;
  }
  jac->JacVector[iTerm + isubvar] = jacEntry; 
  jac->Property_Value = propValue;
  jac->Num_Terms += num_species;
}
/**************************************************************************/
/**************************************************************************/
/**************************************************************************/

void
propertyJac_searchadd(PROPERTYJAC_STRUCT *jac, int varType, int matID, 
		      int isubvar, double jacEntry, double propValue)

    /**********************************************************************
     *
     * propertyJac_searchadd()
     *
     *  Adds an entry to the propertyjac structure. Note: this function
     *  is complete and slow. This is a candidate for inlining.
     **********************************************************************/
{
  int i, iTerm, notfound = TRUE;
  VARIABLE_DESCRIPTION_STRUCT *vd;
  for (i = 0; i <jac->Num_Terms && notfound; i++) {
    if (jac->Var_Type[i] == varType) {
      if (jac->MatID[i] == matID) {
	if (jac->idof[i] == isubvar) {
	  if (fabs(propValue - jac->Property_Value) > 1.0E-5) {
            EH(GOMA_ERROR,"Incompatible propertyValue search");
	  }
	  iTerm = i;
	}
      }
    }
  } 
  if (notfound) iTerm = jac->Num_Terms;
  if (iTerm >= jac->NUM_TERMS_MALLOC) {
    propertyJac_realloc(&jac, iTerm + 1);
  }
  vd = get_vd_ptr(varType, matID, isubvar);
  jac->Var_List[iTerm] = vd;
  jac->idof[iTerm] = isubvar; 
  jac->Var_Type[iTerm] = varType;
  jac->JacVector[iTerm] = jacEntry;
  jac->Num_Terms++;
}
/**************************************************************************/
/**************************************************************************/
/**************************************************************************/

int
propertyJac_find_species_unk(PROPERTYJAC_STRUCT *jac)

    /**********************************************************************
     *
     * propertyJac_addentry()
     *
     *  This function will return the position of the first MASS_Fraction
     *  or SPECIES_UNK_0 dependency for the current propertyJac structure.
     *  If it doesn't find any in the structure it returns -1. 
     **********************************************************************/
{
  int i;
  for (i = 0; i < jac->Num_Terms; i++) {
    if (jac->Var_Type[i] == MASS_FRACTION ||
	jac->Var_Type[i] == SPECIES_UNK_0) {
      return i;
    }
  }
  return -1;
}
/**************************************************************************/
/**************************************************************************/
/**************************************************************************/
