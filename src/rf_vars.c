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
 *$Id: rf_vars.c,v 5.1 2007-09-18 18:53:47 prschun Exp $
 */

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <strings.h>

#ifndef DBL_MAX
#define DBL_MAX 1.0E300
#endif

#include "ac_conti.h"
#include "ac_hunt.h"
#include "ac_particles.h"
#include "ac_stability.h"
#include "ac_stability_util.h"
#include "ac_update_parameter.h"
#include "bc_colloc.h"
#include "bc_contact.h"
#include "bc_curve.h"
#include "bc_dirich.h"
#include "bc_integ.h"
#include "bc/rotate.h"
#include "bc_special.h"
#include "bc_surfacedomain.h"
#include "dp_comm.h"
#include "dp_map_comm_vec.h"
#include "dp_types.h"
#include "dp_utils.h"
#include "dp_vif.h"
#include "dpi.h"
#include "el_elm.h"
#include "el_elm_info.h"
#include "el_geom.h"
#include "el_quality.h"
#include "exo_conn.h"
#include "exo_struct.h"
#include "loca_const.h"
#include "md_timer.h"
#include "mm_as.h"
#include "mm_as_alloc.h"
#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "mm_augc_util.h"
#include "mm_bc.h"
#include "mm_chemkin.h"
#include "mm_dil_viscosity.h"
#include "mm_eh.h"
#include "mm_elem_block_structs.h"
#include "mm_fill.h"
#include "mm_fill_aux.h"
#include "mm_fill_fill.h"
#include "mm_fill_jac.h"
#include "mm_fill_ls.h"
#include "mm_fill_porous.h"
#include "mm_fill_potential.h"
#include "mm_fill_pthings.h"
#include "mm_fill_ptrs.h"
#include "mm_fill_rs.h"
#include "mm_fill_shell.h"
#include "mm_fill_solid.h"
#include "mm_fill_species.h"
#include "mm_fill_stress.h"
#include "mm_fill_terms.h"
#include "mm_fill_util.h"
#include "mm_flux.h"
#include "mm_input.h"
#include "mm_interface.h"
#include "mm_more_utils.h"
#include "mm_mp.h"
#include "mm_mp_const.h"
#include "mm_mp_structs.h"
#include "mm_ns_bc.h"
#include "mm_numjac.h"
#include "mm_post_def.h"
#include "mm_post_proc.h"
#include "mm_prob_def.h"
#include "mm_qtensor_model.h"
#include "mm_shell_bc.h"
#include "mm_shell_util.h"
#include "mm_sol_nonlinear.h"
#include "mm_species.h"
#include "mm_std_models.h"
#include "mm_unknown_map.h"
#include "mm_viscosity.h"
#include "rd_dpi.h"
#include "rd_exo.h"
#include "rd_mesh.h"
#include "rf_allo.h"
#include "rf_bc.h"
#include "rf_bc_const.h"
#include "rf_element_storage_const.h"
#include "rf_element_storage_struct.h"
#include "rf_fem.h"
#include "rf_fem_const.h"
#include "rf_fill_const.h"
#include "rf_io.h"
#include "rf_io_const.h"
#include "rf_io_structs.h"
#include "rf_masks.h"
#include "rf_mp.h"
#include "rf_node_const.h"
#include "rf_pre_proc.h"
#include "rf_shape.h"
#include "rf_solve.h"
#include "rf_solver.h"
#include "rf_solver_const.h"
#include "rf_util.h"
#include "rf_vars_const.h"
#include "sl_aux.h"
#include "sl_auxutil.h"
#include "sl_eggroll.h"
#include "sl_lu.h"
#include "sl_matrix_util.h"
#include "sl_umf.h"
#include "sl_util.h"
#include "std.h"
#include "user_ac.h"
#include "user_bc.h"
#include "user_mp.h"
#include "user_mp_gen.h"
#include "user_post.h"
#include "user_pre.h"
#include "wr_dpi.h"
#include "wr_exo.h"
#include "wr_side_data.h"
#include "wr_soln.h"


/*
 * Global variables
 */
/*
 * Var_Info Array
 *
 * Array of VARIABLE_DSCRIPTION structures which is indexed by the
 * tags at each node contained in the Variable_Tag array defined
 * above.  These entries store information on
 * the existence, type, size and ordering of variable information at nodes.
 */

VARIABLE_DESCRIPTION_STRUCT **Var_Info = NULL;

/* Global number of different Var_Info entries */

int Num_Var_Info_Records = 0;

NODAL_VARS_STRUCT **Nodal_Vars_List = NULL;
int Nodal_Vars_List_Length = 0;


/************************************************************************/
/************************************************************************/
/************************************************************************/
int
find_base_variable_type(int var_type)

    /*****************************************************************
     *
     * find_base_variable_type()
     *
     *  Returns the base variable type given the variable type
     *
     *  NOTE: see rf_fem_const.h for information.
     *****************************************************************/
{
  if (var_type < V_FIRST) return -1;
  if (var_type >= V_LAST) return -1;
  if (var_type >= VELOCITY1 && var_type <= VELOCITY3) return VELOCITY1;
  if (var_type >= VELOCITY_GRADIENT11 &&
      var_type <= VELOCITY_GRADIENT33) return VELOCITY_GRADIENT11;
  if (var_type >= POLYMER_STRESS11_1 &&
      var_type >= POLYMER_STRESS33_7) return  POLYMER_STRESS11_1;
  if (var_type >= SPECIES_UNK_0 &&
      var_type <= SPECIES_UNK_LAST) return SPECIES_UNK_0;
  if (var_type >= VOLF_PHASE_0 &&
      var_type <= VOLF_PHASE_LAST) return VOLF_PHASE_0;  
  if(var_type >= VORT_DIR1 && var_type <= VORT_DIR3)
    return VORT_DIR1;
  /*
   *  If var_type doesn't fall into any special cases, return the
   *  var_type as the base variable type
   */
  return (var_type);
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

int variable_description_comparison(VARIABLE_DESCRIPTION_STRUCT *vd1,
				    VARIABLE_DESCRIPTION_STRUCT *vd2)

    /*****************************************************************
     *
     * variable_description_comparison():
     *
     *  Compare two the information in two variable types. If they
     *  are the same then return TRUE. If they differ in the
     *  significant fields, return FALSE.
     *
     *  NOTE:
     *    Currently, the 4 significant fields are Variable_Type
     *    and MatID and Ndof and subvar_index
     *****************************************************************/
{
  if (vd1 == NULL)                              return FALSE;
  if (vd2 == NULL)                              return FALSE;
  if (vd1->Variable_Type != vd2->Variable_Type) return FALSE;
  if (vd1->MatID         != vd2->MatID)         return FALSE;
  if (vd1->Ndof          != vd2->Ndof)          return FALSE;
  if (vd1->Subvar_Index  != vd2->Subvar_Index)  return FALSE;
  return TRUE;
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

int
variable_description_init(VARIABLE_DESCRIPTION_STRUCT **vd_hdl)

    /*******************************************************************
     *
     * variable_description_init:
     *
     *  Mallocs and initializes a variable description structure
     *  to its default state
     *
     *  Input
     *   *el_hdl => If NULL will allocate a structure. If nonnull it will
     *              assume that the structure has already been malloced.
     *  Return:
     *    Success:           CPC_SUCCESS
     *    Interface Failure: CPC_PUB_BAD
     *    *el_hdl   => address of the new structure that has just been
     *                 malloced.
     ******************************************************************/
{
  if (*vd_hdl == NULL) {
    *vd_hdl = smalloc(sizeof(VARIABLE_DESCRIPTION_STRUCT));
    if (*vd_hdl == NULL) EH(GOMA_ERROR, "variable_description_init");
  }
  return variable_description_default(*vd_hdl);
}

/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

int
variable_description_default(VARIABLE_DESCRIPTION_STRUCT *vd)

   /*******************************************************************
   *
   * variable_description_default():
   *
   *  Initialize a variable description structure to its default state.
   *  This means that all pointers are set to the NULL value.
   ********************************************************************/
{
  if (vd == NULL) return 1;
  /*
   *  First set everything in the structure to zero. This has the effect of
   *  setting pointers to NULL, ints to zero and/or false, and doubles to
   *  the value of zero.
   */
  (void) memset((void *)vd, 0, sizeof(VARIABLE_DESCRIPTION_STRUCT));
  return 0;
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

int
variable_description_alloc(VARIABLE_DESCRIPTION_STRUCT *vd,
                           const int var_type, const int matID,
			   const int ndof, const int subvarIndex)
    
   /*********************************************************************
    *
    *   variable_description_alloc():
    *
    *  Allocates space in a single variable description structure
    *  for the substructures of that structure, and then fills in a
    *  couple key parameters in that structure.
    *
    *  *Var_Name
    *  *double_work
    *  *integer_work
    *********************************************************************/
{
  char tmp[132];
  int i;
  vd->Variable_Type = var_type;
  vd->Ndof          = ndof;
  vd->MatID         = matID;
  vd->Subvar_Index  = subvarIndex;
  if (ndof > 1) {
    for (i = 0; i < ndof; i++) {
      strcpy(tmp, Var_Name[var_type].name1);
      sprintf(tmp + strlen(tmp), "_ndof%d", i);
      vd->Var_Name[i] = alloc_copy_string(tmp);
    }
  } else {
    vd->Var_Name[0] = alloc_copy_string(Var_Name[var_type].name1);
  }

  /*
   * Fill in base variable type here
   */
  vd->Base_Variable_Type = find_base_variable_type(var_type);

  /* 
   *  Bounds information gets initialized here
   */
  for (i = 0; i < ndof; i++) {
    vd->Upper_Bound[i] = DBL_MAX;
    vd->Lower_Bound[i] = -DBL_MAX;
    vd->Delta_Bound[i] = 100.;
    vd->Common_Value[i] = 1.0;	
  }
  return 0;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void
variable_description_free(VARIABLE_DESCRIPTION_STRUCT *vd)

    /*************************************************************************
     *
     * variable_description_free():
     *
     *  Frees memory allocated underneath the variable_description
     *  structure, and returns the structure to its default condition.
     *  This is the reverse of the _alloc() operation.
     ************************************************************************/   
{
  int i;
  for (i = 0; i < vd->Ndof; i++) {
    safer_free((void **) &(vd->Var_Name[i]));
  }
  safer_free((void **) &(vd->double_work));
  safer_free((void **) &(vd->integer_work));  
  (void) variable_description_default(vd);
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void
variable_description_destroy(VARIABLE_DESCRIPTION_STRUCT **vd_hdl)

    /*************************************************************************
   *
   * variable_description_destroy():
   *
   * Free memory required by the variable description structure.
   * This routine basically unmallocs whatever variable_description_alloc
   * malloced previously, and then unmallocs the variable
   * description structure itself.
   *
   *  Note, this is the opposite of the _create routine.
   ***************************************************************************/
{
  variable_description_free(*vd_hdl);
  safer_free((void **) vd_hdl);
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

VARIABLE_DESCRIPTION_STRUCT *
variable_description_create(const int var_type, const int matID,
			    const int ndof, const int subvarIndex)

    /********************************************************
     *
     * variable_description_create():
     *
     *  Allocates space for a single variable description structure,
     *  and then fills it with a couple key parameters.
     *
     ***********************************************************/
{
  int retn;
  VARIABLE_DESCRIPTION_STRUCT  *vd = NULL;
  variable_description_init(&vd);
  if (vd != NULL) {
    retn = variable_description_alloc(vd, var_type, matID, ndof,
				      subvarIndex);
    if (retn != 0) {
      variable_description_destroy(&vd);
    }
  }
  return vd;
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

VARIABLE_DESCRIPTION_STRUCT  *
find_or_create_vd(const int var_type, const int ndof, const int mn,
		  const int subvarIndex, const int imtrx)
				
    /********************************************************************
     *
     * find_or_create_vd():
     *
     *  Allocates space for a single variable description structure,
     *  and then fill it with a couple key parameters.
     *
     ********************************************************************/
{
  int matID, i;
  VARIABLE_DESCRIPTION_STRUCT *vd;
  PROBLEM_DESCRIPTION_STRUCT *pd_ptr;
  
  if (ndof <= 0) return NULL;
  /*
   * Determine if matID should have a generic -1 value or the
   * index of the material in the global list of materials
   * If mn equals -1 on input, we will assume that matID
   * should equal -1 also.
   */
  if (mn >= 0) {
    pd_ptr = pd_glob[mn];
    matID = -1;
    if (pd_ptr->v[imtrx][var_type] & V_MATSPECIFIC) {
      matID = mn;
    }
  } else {
    matID = -1;
  }

  /*
   * Create a new structure holding the variable description
   */
  vd = variable_description_create(var_type, matID, ndof, subvarIndex);
  
  /*
   * Try to find a match with existing variable description
   * structures. If you find a match, then unmalloc the temporary
   * variable description structure and return the address of
   * the existing structure.
   */
  for (i = 0; i < Num_Var_Info_Records; i++) {
    if (variable_description_comparison(Var_Info[i], vd)) {
      variable_description_destroy(&vd);
      return Var_Info[i];
    }
  }

  /*
   * Then, attach the new unique vd struct to the end of the global list
   */
  realloc_ptr_1((void ***)&Var_Info,
		Num_Var_Info_Records+1, Num_Var_Info_Records);
  Var_Info[Num_Var_Info_Records] = vd;
  Num_Var_Info_Records++;

  /*
   * Fill in the list index field here -> use the index in the
   * processor list of variable types.
   */
  vd->List_Index = Num_Var_Info_Records - 1;
  
  /*
   * Perhaps do a better job with the specifying the variable name here
   */

  /*
   * Do a better job with filling in the bounds information here
   */

  /*
   * Perhaps do something with the work arrays here
   */
  
  /*
   * Return the pointer
   */
  return vd;
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

VARIABLE_DESCRIPTION_STRUCT *
Variable_Description_Match(VARIABLE_DESCRIPTION_STRUCT *vd_match)

    /********************************************************************
     *
     * Variable_Description_Match():
     *
     * Tries to find a match of an input variable description structure 
     * to the variable description structores storred in the global list
     * If it finds a match it returns the pointer to the matching 
     * structure in the global list. If it can't find a match, it returns
     * NULL.
     ********************************************************************/
{
  int i, mat_save;
  if (! vd_match) return NULL;
  for (i = 0; i < Num_Var_Info_Records; i++) {
    if (variable_description_comparison(Var_Info[i], vd_match)) {
      return Var_Info[i];
    }
  }
  if (vd_match->MatID != -1) {
    mat_save = vd_match->MatID;
    vd_match->MatID = -1;
    for (i = 0; i < Num_Var_Info_Records; i++) {
      if (variable_description_comparison(Var_Info[i], vd_match)) {
	vd_match->MatID = mat_save;
	return Var_Info[i];
      }
    }
    vd_match->MatID = mat_save;
  }
  return NULL;
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

VARIABLE_DESCRIPTION_STRUCT *
get_vd_ptr(int var_type, int matID, int idof)

    /********************************************************************
     *
     *  get_vd_ptr()
     *
     *   This is a quick lookup function for looking up a matching
     *  variable description structure. It returns the pointer to the
     *  matching variable description structure in the global list. If
     *  it can't find a match, it returns NULL.
     ********************************************************************/
{ 
  int i;
  VARIABLE_DESCRIPTION_STRUCT *possible = NULL, *vd;
  for (i = 0; i < Num_Var_Info_Records; i++) {
    vd = Var_Info[i];
    if (vd->Variable_Type != var_type) continue;
    if (vd->MatID == matID) {
      if (vd->Variable_Type == MASS_FRACTION) {
        if (vd->Subvar_Index != idof) continue;
      }
      return vd;
    } else if (vd->MatID == -1) {
      if (vd->Variable_Type == MASS_FRACTION) {
        if (vd->Subvar_Index != idof) continue;
      }
      possible = vd;
    } else {
      continue;
    }
  }
  return possible;
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

void
vdesc_augment(void)

    /****************************************************************
     *
     * vdesc_augment():
     *
     *    This routine adds likely variable descriptions of variable
     *    that are not included in the solution vector, but are of
     *    utility. 
     ****************************************************************/
{
  int i, k, have_generic = FALSE, have_MF = FALSE, 
      have_specific_MF = FALSE;
  int imtrx;
  MATRL_PROP_STRUCT *mp_local;
  VARIABLE_DESCRIPTION_STRUCT *vd;
  /*
   * First look for a generic mass fraction
   */
   vd = get_vd_ptr(MASS_FRACTION, -1, 0);
   if (vd) 
     {
      have_generic = TRUE;
      have_MF = TRUE;
     }
  /*
   * Next look for a specific mass fraction unknown
   */
  for (i = 0; i < upd->Num_Mat; i++) 
     {
      mp_local = mp_glob[i];
      vd = get_vd_ptr(MASS_FRACTION, i, 0);
      if (vd) 
        {
         have_MF = TRUE;
         have_specific_MF = TRUE;
        }
     }
  if (have_MF) 
    {
     if (have_generic) 
       {
        for (k = 0; k < upd->Max_Num_Species; k++) 
           {
            for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++)
               {
	        vd = find_or_create_vd(MASS_FRACTION, 1, -1, k, imtrx);
               }
           }
       }
     if (have_specific_MF) 
       {
        for (i = 0; i < upd->Num_Mat; i++) 
           {
	    mp_local = mp_glob[i];
	    for (k = 0; k < mp_local->Num_Species; k++) 
               {
                for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++)
                   {
	            vd = find_or_create_vd(MASS_FRACTION, 1, i, k, imtrx);
                   }
	       }
           }
       }
    }
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void 
assign_species_desc_suffix(const int species_var_name, char * retn_string)

   /***************************************************************************
    *
    * assign_species_desc_suffix:
    *
    *  This function assigns a string to the return string, associated with the
    *  type of species variable being considered.
    *
    *  Input
    * --------
    *  species_var_name: Valid species type
    *                   (The valid species types are listed in rf_fem_const.h)
    * Output
    * --------
    *  retn_string:  Has to have an input length of at least 24
    ***************************************************************************/
{
  if (retn_string == NULL) {
   EH(GOMA_ERROR, "assign_species_desc_suffix: bad interface\n");
  }
  switch (species_var_name) {
  case SPECIES_MASS_FRACTION:
      (void) strcpy(retn_string, "Mass Fraction");
      break;
  case SPECIES_MOLE_FRACTION:
      (void) strcpy(retn_string, "Mole Fraction");
      break;
  case SPECIES_VOL_FRACTION:
      (void) strcpy(retn_string, "Volume Fraction");
      break;
  case SPECIES_DENSITY:
      (void) strcpy(retn_string, "Density");
      break;
  case SPECIES_CONCENTRATION:
      (void) strcpy(retn_string, "Concentration");
      break;
  case SPECIES_CAP_PRESSURE:
      (void) strcpy(retn_string, "Capillary Pressure");
      break;
  case SPECIES_UNDEFINED_FORM:
      (void) strcpy(retn_string, "Species Unknown");
      break;
  default:
      (void) strcpy(retn_string, "Mass Fraction(default)");
      printf("assign_species_desc_suffix: WARNING unknown form of "
	     "species vars: %d\n", species_var_name);
  }
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void 
assign_species_name(const int kspec, struct Material_Properties *mp_ptr,
		    char *species_name, char *species_desc, int matID)

    /*************************************************************************
     *
     * assign_species_name
     *
     *   Given a species index and a material property structure this
     *   function returns the Exodus species name and a descriptive
     *   string.
     *
     *  Input
     * ---------
     *  kspec         Index of the species in the current material domain
     *  mp_ptr        Pointer to the Material_Properties struct for the
     *                current material domain
     *  matID         Value of the matID field in the variable description
     *                structure. Usually, this is -1. However, if it is 
     *                positive, then the variable is material specific,
     *                and this routine will add a suffix unto the names.
     *  Output
     * ---------
     *  species_name  Exodus variable name
     *  species_desc  Description of the degree of freedom.
     ************************************************************************/
{
  int sr;
  char pstring[256], species_desc_end[256];
  if (species_name == NULL) {
   EH(GOMA_ERROR, "assign_species_name: bad interface\n");
  }

  if (mp_ptr->Num_Species > kspec) {
    /*
     *  Obtain the species prefix name and species desc suffix
     */
    assign_species_prefix(mp_ptr->Species_Var_Type, species_name);
    assign_species_desc_suffix(mp_ptr->Species_Var_Type, species_desc_end);

    /*
     *   HKM -
     * If we have a real species name, let's use that!
     * Also use prefixes to denote the known types of species variables.
     * For backwards compatibility, if the species names are at their
     * defaults, then we will use the old "Y0" and "Y1" nomenclature.
     */
    sprintf(pstring, "Species_%d", kspec);
    if ((mp_ptr->Species_Names != NULL)                 &&
	(strlen(mp_ptr->Species_Names[kspec]) > 0)      &&
	(strcmp(pstring, mp_ptr->Species_Names[kspec]))
	) {
      strcat(species_name, mp_ptr->Species_Names[kspec]);
      if (species_desc != NULL) {
	sprintf(species_desc, "Species %s %s",
		mp_ptr->Species_Names[kspec], species_desc_end);
      }
    } else {
      sr = strlen(species_name);
      sprintf(species_name + sr, "%d", kspec);
      if (species_desc != NULL) {
	sprintf(species_desc, "Species %d %s", kspec, species_desc_end);
      }
    }
    if (matID >= 0) {
      if ( strlen(mp_ptr->Material_Name) != 0 ) {
	sr = strlen(species_name);
	sprintf(species_name + sr, "_%-.8s", mp_ptr->Material_Name);
	sr = strlen(species_desc);
	sprintf(species_desc + sr," for %s",  mp_ptr->Material_Name);
      } else {
	sr = strlen(species_name);
	sprintf(species_name + sr, "_mat%d", matID);
	sr = strlen(species_desc);
	sprintf(species_desc + sr," for Material %d", matID);
      }
    }
  } else {
    printf("assign_species_name warning: Called with a bad species index: %d",
	   kspec);
  }
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void 
assign_var_name(const int varType, const int idof, 
		MATRL_PROP_STRUCT *mp_ptr,
		char *species_name, char *species_desc, int matID)

    /**************************************************************************
     *
     * assign_var_name
     *
     *   Given a species index and a material property structure this
     *   function returns the Exodus species name and a descriptive
     *   string.
     *
     *  Input
     * ---------
     *  kspec         Index of the species in the current material domain
     *  mp_ptr        Pointer to the Material_Properties struct for the
     *                current material domain
     *  matID         Value of the matID field in the variable description
     *                structure. Usually, this is -1. However, if it is 
     *                positive, then the variable is material specificy, 
     *                and this routine will add a suffix unto the names.
     *  Output
     * ---------
     *  species_name  Exodus variable name
     *  species_desc  Description of the degree of freedom.
     ***********************************************************************/
{
  int sr;
  if (varType == MASS_FRACTION) {
    assign_species_name(idof, mp_ptr, species_name, species_desc, matID);
  } else {
    strcpy(species_name, Var_Name[varType].name2);
    strcpy(species_desc, Var_Name[varType].name1);

    if (matID >= 0) {
      if ( strlen(mp_ptr->Material_Name) != 0 ) {
	sr = strlen(species_name);
	sprintf(species_name + sr, "_%-.8s", mp_ptr->Material_Name);
	sr = strlen(species_desc);
	sprintf(species_desc + sr," for %s",  mp_ptr->Material_Name);
      } else {
	sr = strlen(species_name);
	sprintf(species_name + sr, "_mat%d", matID);
	sr = strlen(species_desc);
	sprintf(species_desc + sr," for Material %d", matID);
      }
    }
  }
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
