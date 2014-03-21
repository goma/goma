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
 *$Id: mm_interface.c,v 5.3 2008-11-06 15:37:04 hkmoffa Exp $
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "std.h"
#include "rf_allo.h"
#include "rf_vars_const.h"
#include "mm_mp_const.h"
#include "mm_as_structs.h"
#include "mm_mp_structs.h"
#include "mm_mp.h"
#include "mm_species.h"
#include "rf_bc_const.h"
#include "mm_fill_jac.h"
#include "mm_interface.h"
#include "mm_eh.h"

/************************************************************************/
/************************************************************************/
/************************************************************************/

INTERFACE_SOURCE_STRUCT *
interface_source_alloc(int num_terms, int is_type, int do_Jac)

   /********************************************************************
    *
    * interface_source_alloc():
    *
    *    Allocate and zero an Interface_Source structure.
    *
    * Return:
    *   This function returns the pointer to the malloced and initialized
    *   structure.
    *********************************************************************/
{
  INTERFACE_SOURCE_STRUCT *is;
  is = alloc_struct_1(INTERFACE_SOURCE_STRUCT, 1);
  is->Num_Terms = num_terms;
  is->SpeciesVT = SPECIES_UNDEFINED_FORM;
  is->IS_Type = is_type;
  is->Var_List = (VARIABLE_DESCRIPTION_STRUCT **) alloc_ptr_1(num_terms);
  is->idof = alloc_int_1(num_terms, 0);
  is->Var_Value = alloc_dbl_1(num_terms, 0.0);
  is->SourceTerm = alloc_dbl_1(num_terms, 0.0);
  is->JacMatrix = alloc_dbl_2(num_terms, num_terms, 0.0);
  is->Do_Jac = do_Jac;
  /*
   * Possibly branch depending upon is_type
   */
  return is;
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

void
interface_source_free(INTERFACE_SOURCE_STRUCT *is)

   /********************************************************************
    *
    * interface_source_free():
    *
    *    This routine frees all of the underlying memory associated
    *  with an Interface_Source structure.
    ********************************************************************/
{
  /*
   *  Free underlying memory assocated with the state of the surface,
   *  If there is any. IS_Type tells us what to do.
   */
  if (is->StateInterface) {

  }
  if (is->StateInterfaceOld) {

  }
  /*
   * Free memory allocated directly in the structure
   */
  safer_free((void **) &(is->Var_List));
  safer_free((void **) &(is->idof));
  safer_free((void **) &(is->Var_Value));
  safer_free((void **) &(is->SourceTerm));
  safer_free((void **) &(is->JacMatrix));
  is->Num_Terms = 0;
  is->IS_Type = 0;
  is->Processed = FALSE;
  is->SpeciesVT = SPECIES_UNDEFINED_FORM;
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

void
interface_source_destroy(INTERFACE_SOURCE_STRUCT **is_hdl)

   /********************************************************************
    *
    * interface_source_destroy():
    *
    *    This routine frees all of the underlying memory associated
    *  with an Interface_Source structure. It then goes on to free
    *  the structure itself.
    *********************************************************************/
{
  interface_source_free(*is_hdl);
  safer_free((void **) is_hdl);
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

void
interface_source_zero(INTERFACE_SOURCE_STRUCT *is)

   /********************************************************************
    *
    * interface_source_zero():
    *
    *    Zero the source terms and Jacobian in the structure. The state
    * of the interface and information concerning what variables are
    * defined at the interface aren't touched. Set the Processed flag
    * to FALSE to denote that there no longer is useful source term
    * information in the structure.
    *********************************************************************/
{
  (void) memset((void *) is->SourceTerm, 0, sizeof(double) * is->Num_Terms);
  if (is->Do_Jac) {
    (void) memset((void *) is->JacMatrix[0], 0,
		  sizeof(double) * is->Num_Terms * is->Num_Terms);
  }
  is->Processed = FALSE;
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

void 
is_change1_speciesVT(INTERFACE_SOURCE_STRUCT *is, const int is_species_entry, 
                     const int is_species_start, MATRL_PROP_STRUCT *mp_local, 
		     const int speciesVT, const double time)


   /********************************************************************
    *
    * is_change1_speciesVT
    *
    *   Change the Jacobian dependent
    *   Variable, if the species variable type changes.
    *
    * Input
    *   
    *   speciesVT = Desired dependent variable for the Jacobian entry
    *   
    *********************************************************************/
{
  double *d_ptr;
  if (speciesVT == SPECIES_MASS_FRACTION) {
    switch (is->SpeciesVT) {
    case SPECIES_MOLE_FRACTION:
	/*
	 * Need to change the source term units as well
	 */
	EH(-1,"not implemented");
	break;
    case SPECIES_MASS_FRACTION:
	break;
    case SPECIES_CONCENTRATION:
        if (is->Do_Jac) {
	  d_ptr = &(is->JacMatrix[is_species_entry][is_species_start]);
	  deriv1_Ck_to_Yk(d_ptr, mp_local, is->Var_Value + is_species_start, time);
	}
	break;
    default:
	EH(-1,"not implemented");
	break;
    }
  } else if (speciesVT == SPECIES_MOLE_FRACTION) {
   switch (is->SpeciesVT) {
    case SPECIES_MASS_FRACTION:	
        /*
	 * Need to change the source term units as well
	 */
	EH(-1,"not implemented");
	break;
    case SPECIES_MOLE_FRACTION:
	break;
   case SPECIES_CONCENTRATION:
  	EH(-1,"not implemented");
        break;
    default:
	EH(-1,"not implemented");
	break;
    }
  } else if (speciesVT == SPECIES_CONCENTRATION) {
   switch (is->SpeciesVT) {
    case SPECIES_MASS_FRACTION:
	/*
	 * Need to change the source term units as well
	 */
	EH(-1,"not implemented");
	break;
    case SPECIES_MOLE_FRACTION:
	EH(-1,"not implemented");
	break;
    case SPECIES_CONCENTRATION:
	break;
    default:
	EH(-1,"not implemented");
	break;
    }

  }
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

void 
is_change1_lastspecies(INTERFACE_SOURCE_STRUCT *is,
		       const int is_species_entry, 
		       const int is_species_start,
		       MATRL_PROP_STRUCT *mp_local)

   /********************************************************************
    *
    * is_change1_lastspecies
    *
    *   This changes what gets held constant during the 
    *   partial derivative of a single entry as a function of the
    *   species variable.
    *
    *   On input it is assumed that a derivative wrt Y_k is carried
    *   out assuming [Y_j, j .ne. k, 1 < j < N ] is held constant.
    *   On output is is assumed that a derivative wrt Y_k is carried
    *   out assuming [Y_j, j .ne. k, 1 < j < N-1 and sum Y_i = 1]
    *   is held constant. What this means in practice is that a change
    *   in Y_j is conteracted by allied change in Y_N such that the
    *   sum Y_i = 1 always.
    *
    *********************************************************************/
{
  int i, nsm1;
  double *d_ptr, tmp;
  
  switch (is->SpeciesVT) {
  case SPECIES_MOLE_FRACTION:
  case SPECIES_MASS_FRACTION:
  case SPECIES_CONCENTRATION:
  case SPECIES_DENSITY:
      if (mp_local->Num_Species_Eqn < mp_local->Num_Species) {
	if (is->Do_Jac) {
	  d_ptr = &(is->JacMatrix[is_species_entry][is_species_start]);
	  nsm1 = mp_local->Num_Species - 1;
	  tmp = d_ptr[nsm1];
	  if (DOUBLE_NONZERO(tmp)) {
	    for (i = 0; i < nsm1; i++) {
	      d_ptr[i] -= tmp;
	    }
	    d_ptr[nsm1] = 0.0;
	  }
	}
      }
      break;
  default:
      EH(-1,"not implemented");
      break;
  }
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

int
match_interface_source_string(char *istring)

    /********************************************************************
     *
     *  match_interface_source_string():
     *
     *
     *
     ********************************************************************/
{
  int match = -1;
  if (!strcasecmp(istring, "IS_EQUIL_PSEUDORXN") ||
      !strcasecmp(istring, "is_equil_prxn_bc")) {
    match = IS_EQUIL_PRXN_BC;
  }
  else if (!strcasecmp(istring, "VL_EQUIL_PSEUDORXN") ||
	   !strcasecmp(istring, "vl_equil_prxn_bc")) {
    match = VL_EQUIL_PRXN_BC;
  }
  else if (!strcasecmp(istring, "SURFDOMAINCHEMKIN_SURFRXN") ||
	   !strcasecmp(istring, "SDC_SURFRXN_BC")) {
    match = SDC_SURFRXN_BC;
  }
  if (match == -1) {
    fprintf(stderr,
	    "match_interface_source_string ERORR: string %s not recognized\n",
	    istring);
    EH(-1,"input error");
  }
  return match;
}
/************************************************************************/
/************************************************************************/
/************************************************************************/
