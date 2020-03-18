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
#include <strings.h> /* strcasecmp and strncasecmp moved here for POSIX.1 */

#include "std.h"
#include "rf_allo.h"
#include "rf_vars_const.h"
#include "rf_fem_const.h"
#include "rf_fem.h"
#include "mm_mp_structs.h"
#include "mm_species.h"
#include "rf_bc_const.h"
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
  int i;
  INTERFACE_SOURCE_STRUCT *is;
  is = alloc_struct_1(INTERFACE_SOURCE_STRUCT, Num_Interface_Srcs);
  for(i=0; i<Num_Interface_Srcs ; i++)    {
    is[i].Num_Terms = num_terms;
    is[i].SpeciesVT = SPECIES_UNDEFINED_FORM;
    is[i].IS_Type = is_type;
    is[i].Var_List = (VARIABLE_DESCRIPTION_STRUCT **) alloc_ptr_1(num_terms);
    is[i].idof = alloc_int_1(num_terms, 0);
    is[i].Var_Value = alloc_dbl_1(num_terms, 0.0);
    is[i].SourceTerm = alloc_dbl_1(num_terms, 0.0);
    is[i].Processed = alloc_int_1(MAX_CONC, FALSE);
    is[i].JacMatrix = alloc_dbl_2(num_terms, num_terms, 0.0);
    is[i].Do_Jac = do_Jac;
    }
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
  int i;
  /*
   *  Free underlying memory assocated with the state of the surface,
   *  If there is any. IS_Type tells us what to do.
   */
 for(i=0 ; i<Num_Interface_Srcs ; i++)    {
  if (is[i].StateInterface) {

  }
  if (is[i].StateInterfaceOld) {

  }
  /*
   * Free memory allocated directly in the structure
   */
  safer_free((void **) &(is[i].Var_List));
  safer_free((void **) &(is[i].idof));
  safer_free((void **) &(is[i].Var_Value));
  safer_free((void **) &(is[i].SourceTerm));
  safer_free((void **) &(is[i].JacMatrix));
  safer_free((void **) &(is[i].Processed));
  is[i].Num_Terms = 0;
  is[i].IS_Type = 0;
  is[i].SpeciesVT = SPECIES_UNDEFINED_FORM;
 }
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
  (void) memset((void *) is->Processed, FALSE, sizeof(int) * MAX_CONC);
  if (is->Do_Jac) {
    (void) memset((void *) is->JacMatrix[0], 0,
		  sizeof(double) * is->Num_Terms * is->Num_Terms);
  }
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

void 
is_change1_speciesVT(INTERFACE_SOURCE_STRUCT *is, const int is_species_entry, 
                     const int is_species_start, MATRL_PROP_STRUCT *mp_local, 
		     const int speciesVT, const double time, const int intf_id)


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
    switch (is[intf_id].SpeciesVT) {
    case SPECIES_MOLE_FRACTION:
	/*
	 * Need to change the source term units as well
	 */
	EH(-1,"not implemented");
	break;
    case SPECIES_MASS_FRACTION:
	break;
    case SPECIES_CONCENTRATION:
        if (is[intf_id].Do_Jac) {
	  d_ptr = &(is[intf_id].JacMatrix[is_species_entry][is_species_start]);
	  deriv1_Ck_to_Yk(d_ptr, mp_local, is[intf_id].Var_Value + is_species_start, time);
	}
	break;
    default:
	EH(-1,"not implemented");
	break;
    }
  } else if (speciesVT == SPECIES_MOLE_FRACTION) {
   switch (is[intf_id].SpeciesVT) {
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
   switch (is[intf_id].SpeciesVT) {
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
		       MATRL_PROP_STRUCT *mp_local, const int intf_id)

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
  
  switch (is[intf_id].SpeciesVT) {
  case SPECIES_MOLE_FRACTION:
  case SPECIES_MASS_FRACTION:
  case SPECIES_CONCENTRATION:
  case SPECIES_DENSITY:
      if (mp_local->Num_Species_Eqn < mp_local->Num_Species) {
	if (is[intf_id].Do_Jac) {
	  d_ptr = &(is[intf_id].JacMatrix[is_species_entry][is_species_start]);
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
  else if (!strcasecmp(istring, "YFLUX_DISC_RXN_BC") ||
	   !strcasecmp(istring, "yflux_disc_rxn_bc")) {
    match = YFLUX_DISC_RXN_BC;
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
