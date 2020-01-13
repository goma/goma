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
 *$Id: rf_node_vars.c,v 5.3 2009-04-24 23:42:33 hkmoffa Exp $
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>


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
#include "bc_rotate.h"
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

/***************************************************************************/
/***************************************************************************/
/***************************************************************************/

int
dof_lnode_var_type(const int n, const int Element_Type,
		   const int proc_node_num, const int var_type,
		   PROBLEM_DESCRIPTION_STRUCT *pd_ptr,
                   const int imtrx)

    /**********************************************************************
     *
     * dof_lnode_var_type:
     *
     *  This is a wrapper around dof_lnode_interp_type(). It calculates
     *  the number of degrees of freedom located at this local node number,
     *  that is actually interpolated on this element (i.e., active on this
     *  element)
     *  given the interpolation type, given the variable type and the current
     *  problem description structure. Within that structure it defines
     *  the interpolation to be used for that variable type on that
     *  element type.
     *  See the dof_lnode_interp_type() description for more information.
     *
     *  proc_node_num is needed to look up whether this node is on the
     *  edge of a domain.
     *
     *  So, if a variable is located at that node, but is not active for
     *  that element, the answer is zero.
     *
     *  The return value is equal to the number of degrees of freedom at
     *  a local node for a variable type corresponding to the current
     *  element, given the problem description structure.
     *********************************************************************/
{
  int retn, interp_type, edge;
  /*
   * Store the interpolation type for the input variable
   */
  interp_type =  pd_ptr->i[imtrx][var_type];

  /*
   * Store whether this node is on an edge of the grid or not
   */
  edge = (int) Nodes[proc_node_num]->EDGE;

  /*
   * Now, given the local node number, n, the element type, and the
   * interpolation type, return the number of degrees of freedom
   */
  retn = dof_lnode_interp_type(n, Element_Type, interp_type, edge);
  EH(retn, "dof_lnode_var_type: ERROR");  
  return retn;
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

int
num_varType_at_node(const int inode, const int varType)

    /********************************************************************
     *
     * num_varType_at_node():
     *
     *  This routine counts up the number of degrees of freedom of 
     * a particular variable type exists at the current node. Note,
     * multiple subvariable types will count extra in this routine.
     * 
     * Note: this is just a wrapper around get_nv_ndofs() that accepts
     *       the node number as input
     *******************************************************************/
{
  NODE_INFO_STRUCT *node = Nodes[inode];
  int num = get_nv_ndofs(node->Nodal_Vars_Info[pg->imtrx], varType);
  return num;
}
/************************************************************************/

int
get_nv_ndofs(NODAL_VARS_STRUCT *nv, const int varType)
    
    /********************************************************************
     *
     * get_vd_ndofs:
     *
     *   This counts up how many degrees of freedom in the nodal
     *   variable correspond to the variable type varType.
     *******************************************************************/
{
  int i, index, sum = 0;
  int num = nv->Num_Var_Desc_Per_Type[varType];
  VARIABLE_DESCRIPTION_STRUCT *vd;
  for (i = 0; i < num; i++) {
    index = nv->Var_Type_Index[varType][i];
    vd	  = nv->Var_Desc_List[index];
    sum	 += vd->Ndof;
  }
  return sum;
}
/************************************************************************/

int
get_nv_ndofs_modMF(NODAL_VARS_STRUCT *nv, const int varType)
    
    /********************************************************************
     *
     * get_vd_ndofs_modMF:
     *
     *   This counts up how many degrees of freedom in the nodal
     *   variable correspond to the variable type varType.
     *
     *   This version counts MASS_FRACTION variable types with
     *   multiple subvariable degrees of freedom as a single
     *   degree of freedom. Note, the algorithm below works now.
     *   However, it may not work later.
     *******************************************************************/
{
  int ndofs = get_nv_ndofs(nv, varType);
  if (varType == MASS_FRACTION) {
    if (upd->Max_Num_Species_Eqn > 0) {
      return ndofs / upd->Max_Num_Species_Eqn;
    }
  }
  return ndofs;
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

int
get_nv_index_varType(NODAL_VARS_STRUCT *nv, const int varType,
		     const int mn, const int subVarType)

    /********************************************************************
     *
     * get_nv_index_varType():
     *
     *
     *       Given a variable type, varType, and material index, mn, and
     * subVarType (for a MASS_FRACTION variable type only), this routine
     * will return the index of the variable description structure in the
     * nodal_Vars structure list, nv->Var_Type_list, that
     * matches. If none matches, this routine will return -1.
     ********************************************************************/
{
  int i, index, good_index = -1, matID;
  int num = nv->Num_Var_Desc_Per_Type[varType];
  VARIABLE_DESCRIPTION_STRUCT *vd;
  for (i = 0; i < num; i++) {
    index = nv->Var_Type_Index[varType][i];
    vd = nv->Var_Desc_List[index];
    matID = vd->MatID;
    if (varType == MASS_FRACTION) {
      if (vd->Subvar_Index == subVarType) {
         if (matID == mn) return index;
         if (matID == -1) good_index = index;       
      }
    } else {
      if (matID == mn) return index;
      if (matID == -1) good_index = index;
    }
  }
  return good_index;
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

int
get_nodal_unknown_offset(NODAL_VARS_STRUCT *nv, const int varType,
			 const int mn, const int subVarType,
			 VARIABLE_DESCRIPTION_STRUCT **vd_ptr)

    /********************************************************************
     *
     * get_nodal_unknown_offset()
     *
     *
     *       Given a variable type, varType, and material index, mn, and
     * subVarType (this is the species number for a MASS_FRACTION 
     * variable type  and a dof for variable descriptions which have
     * multiple degrees of freedom associated with each variable type),
     * this routine returns the offset of the variable from the beginning
     * of the solution vector at this node.
       If none matches, this routine will return -1.
     ********************************************************************/
{
  int i, index, i_match = -1;
  int num = nv->Num_Var_Desc_Per_Type[varType];
  VARIABLE_DESCRIPTION_STRUCT *vd;
  if (varType == MASS_FRACTION) {
    for (i = 0; i < num; i++) {
      index = nv->Var_Type_Index[varType][i];
      vd = nv->Var_Desc_List[index];
      if (vd->Subvar_Index == subVarType) {
	if (vd->MatID == mn || mn == -2) {
	  i_match = index;
	  break;
	} else if (vd->MatID == -1) {
	  i_match = index;
	}
      }
    }
    if (i_match == -1) {
      if (vd_ptr) *vd_ptr = NULL;
      return -1;
    }
    index = nv->Nodal_Offset[i_match];
    if (vd_ptr) *vd_ptr = nv->Var_Desc_List[index];
    return (index);
  } 
  for (i = 0; i < num; i++) {
    index = nv->Var_Type_Index[varType][i];
    vd = nv->Var_Desc_List[index];
    if (vd->MatID == mn || mn == -2) {
      i_match = index;
      break;
    } else if (vd->MatID == -1) {
      i_match = index;
    }
  }
  if (i_match == -1) {
    if (vd_ptr) *vd_ptr = NULL;
    return -1;
  }
  index = nv->Nodal_Offset[i_match] + subVarType;
  if (vd_ptr) *vd_ptr = nv->Var_Desc_List[index];
  return (index);
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

int 
get_nv_offset_idof(NODAL_VARS_STRUCT *nv, const int varType, 
		   const int idof, int subVarType,
		   VARIABLE_DESCRIPTION_STRUCT **vd_ptr)

    /**************************************************************
     *
     * get_nv_offset_idof():
     *
     *      This function returns the offset into the local solution
     * vector of the idof'th occurrence of variable type, varType.
     * MASS_FRACTION variables types also distinguish between 
     * subvartypes. It also returns the address of the variable
     * description structure pertaining to this unknown.
     **************************************************************/
{
  int i, index, offset = -1, nfound;
  int num = nv->Num_Var_Desc_Per_Type[varType];
  VARIABLE_DESCRIPTION_STRUCT *vd;
  if (varType == MASS_FRACTION) {
    for (i = 0, nfound = 0; i < num; i++) {
      index = nv->Var_Type_Index[varType][i];
      vd = nv->Var_Desc_List[index];
      if (vd->Subvar_Index == subVarType) {
	if (nfound == idof) {
	  if (vd_ptr) *vd_ptr = vd;
	  offset = nv->Nodal_Offset[index];
	  return offset;
	} else {
	  nfound++;
	}
      }
    }
    EH(-1,"NOT FOUND");
  } else {
    index = nv->Var_Type_Index[varType][idof];
    if (vd_ptr) *vd_ptr = nv->Var_Desc_List[index];    
    offset = nv->Nodal_Offset[index] + subVarType;
  }
  return offset;
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

void
get_nv_vd_from_offset(NODAL_VARS_STRUCT *nv, int offset,
		      VARIABLE_DESCRIPTION_STRUCT **vd_ptr,  int *idof)

    /*
     *
     */
{
  int i, sum = 0;
  VARIABLE_DESCRIPTION_STRUCT *vd;
  for (i = 0; i < nv->Num_Var_Desc; i++) {
    vd = nv->Var_Desc_List[i];
    if (vd->Ndof + sum > offset) {
      *vd_ptr = vd;
      *idof = offset - sum;
      return;
    }
    sum += vd->Ndof; 
  }
  EH(-1,"RAN OUT OF vds");
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

int
nodal_vars_comparison(NODAL_VARS_STRUCT *nv1, NODAL_VARS_STRUCT *nv2)

    /********************************************************************
     *
     * nodal_vars_comparison():
     *
     *  Compare the information in two nodal_vars structures. If they
     *  are the same then return TRUE. If they differ in the
     *  significant fields, return FALSE.
     *
     *  Currently, the significant fields in a Nodal_Vars
     *  structure are defined to be Num_Var_Desc, Num_Unknowns,
     *  and the information in each Variable_Description structure
     *  listed in Node_Vars_list, as defined by the the function,
     *  variable_description_comparison():
     *
     *  Return:
     *    TRUE if the two structures have the same information
     *    FALSE if they refer to different information
     *******************************************************************/
{
  int i;
  VARIABLE_DESCRIPTION_STRUCT *vd1, *vd2;
  if (nv1 == NULL)                            return FALSE;
  if (nv2 == NULL)                            return FALSE;
  if (nv1->Num_Var_Desc != nv2->Num_Var_Desc) return FALSE;
  if (nv1->Num_Unknowns != nv2->Num_Unknowns) return FALSE;
  for (i = 0; i < nv1->Num_Var_Desc; i++) {
    vd1 = nv1->Var_Desc_List[i];
    vd2 = nv2->Var_Desc_List[i];
    if (!variable_description_comparison(vd1, vd2)) {
      return FALSE;
    }
  }
  return TRUE;
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

int
MatchAdd_Node_Vars_List(NODAL_VARS_STRUCT *newNV,
		        NODAL_VARS_STRUCT **matchNV_hdl)

    /********************************************************************
     *
     * MatchAdd_To_Node_Vars_list():
     *
     *       This routines compares the contents of a new Nodal_Vars
     * struct (assumed to have been previously malloced)
     * to the global list of Nodal_Vars structs
     * kept by the processor. If the new Nodal_Vars struct, newNV,
     * contains new information then the pointer to it is added onto
     * the global list. If it is found to contain duplicate information
     * to an existing item, then no addition to the list is performed.
     *
     *  Input
     * --------
     * newNV -> Pointer to the Nodal_Vars struct
     *
     *  Output
     * --------
     * *matchNV_hdl -> pointer to the Nodal_Vars Struct
     *                that matches the input Nodal_Vars Struct.
     *                If no match was found, this is equal to newNV
     *                on output.
     * Return
     * ----------
     *  If a match was found this returns true. If not, this
     *  returns false
     *
     * NOTES:
     *   Since the list of unique Nodal_Vars structures is assumed to
     *   be very small (a handful) speed in searching through old lists
     *   is not deemed to be important. Therefore, we use a simple
     *   O(n*m) searching algorithm.
     ********************************************************************/
{
  int i;
  NODAL_VARS_STRUCT *oldNV;
  /*
   * Lets determine if we can match any existing entry in the
   * global processor list
   */
  for (i = 0; i < Nodal_Vars_List_Length; i++) {
    oldNV = Nodal_Vars_List[i];
    if (nodal_vars_comparison(oldNV, newNV)) {
      *matchNV_hdl = oldNV;
      return TRUE;
    }
  }
  
  /*
   *  If a match is not found, we need to add the new Nodal
   *  Variables structure to the end of the list.
   *  We do this by reallocating the Nodal_Vars list,
   *  Proc_Nodal_Vars_List.
   */
  realloc_ptr_1((void ***) &Nodal_Vars_List,
		Nodal_Vars_List_Length+1, Nodal_Vars_List_Length);
  Nodal_Vars_List[Nodal_Vars_List_Length] = newNV;
  newNV->list_index = Nodal_Vars_List_Length;
  Nodal_Vars_List_Length++;
  *matchNV_hdl = newNV;
  return FALSE;
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

void
add_var_to_nv_struct(VARIABLE_DESCRIPTION_STRUCT *vd,
		     NODAL_VARS_STRUCT *nv)

    /*****************************************************************
     *
     * add_var_to_nv_struct():
     *
     *   This function will add a variable description structure
     *   into a nodal variable description structure if it is
     *   determined that the variable is unique.
     *****************************************************************/
{
  int i, num_desc;
  int var_type = vd->Variable_Type;
  VARIABLE_DESCRIPTION_STRUCT *vd_old;
 
  if (nv == NULL) EH(-1,"add_var_to_nv_struct error");
  if (vd == NULL) EH(-1,"add_var_to_nv_struct error");

  /*
   * Check to see whether this variable description structure is
   * a duplicate of an existing variable description structure.
   * If it is a duplicate, stop and return.
   */
  for (i = 0; i < nv->Num_Var_Desc; i++) {
    vd_old = nv->Var_Desc_List[i];
    if (variable_description_comparison(vd, vd_old)) {
      return;
    }
  }

  /*
   * Add the new variable description to end of the list and
   * increment the number of variable descriptions.
   */
  realloc_ptr_1((void ***)&(nv->Var_Desc_List),
		nv->Num_Var_Desc + 1, nv->Num_Var_Desc);
  nv->Var_Desc_List[nv->Num_Var_Desc] = vd;

  /*
   *  Add another index entry onto nv->Var_Type_Index[var_type]
   */
  num_desc = nv->Num_Var_Desc_Per_Type[var_type];
  realloc_int_1(&(nv->Var_Type_Index[var_type]), num_desc + 1, num_desc);
  nv->Var_Type_Index[var_type][num_desc] = nv->Num_Var_Desc;

  /*
   * Increment the Num_Var_Desc_Per_Type field
   */
  nv->Num_Var_Desc_Per_Type[var_type]++;

  /*
   * Calculate the new offset variable.
   */
  realloc_int_1(&nv->Nodal_Offset, nv->Num_Var_Desc + 1,
	        nv->Num_Var_Desc);
  nv->Nodal_Offset[nv->Num_Var_Desc] = nv->Num_Unknowns;

  /*
   * Increment the number of variable description structures
   */
  nv->Num_Var_Desc++;
  
  /*
   * Increment the number of unknowns
   */
  nv->Num_Unknowns += vd->Ndof;
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

static void
rearrange_mat_increasing(int var_type, int subvarIndex,
			 NODAL_VARS_STRUCT *nv_old,
			 NODAL_VARS_STRUCT *nv_new, int *rec_used)
    
      /*****************************************************************
       *
       * rearrange_mat_increasing():
       *
       *   Search for the record with the matching variable type and
       *   subvar index that  has the lowest MatID field, that is
       *   currently unused. It then copies that variable type from
       *   the old nodal variable structure to the new.
       *****************************************************************/
{
  int j, mn_min, var_rec, i, num;
  VARIABLE_DESCRIPTION_STRUCT *vd;

  /*
   *  Loop over the number of records in the old structure that
   *  have the given variable type.
   */
  num =  nv_old->Num_Var_Desc_Per_Type[var_type];
  if (var_type == MASS_FRACTION) {
    num =  nv_old->Num_Var_Desc_Per_Type[var_type] /
	   upd->Max_Num_Species_Eqn;
  }
  
  for (j = 0; j < num; j++) {
 
    /*
     * Search for the record with the matching variable type and
     * subvar index that
     * has the lowest MatID field, that is currently unused.
     */
    mn_min = INT_MAX;
    var_rec = -1;
    for (i = 0; i < nv_old->Num_Var_Desc; i++) {
      if (! rec_used[i]) {
	vd = nv_old->Var_Desc_List[i];
	if ((var_type == vd->Variable_Type) &&
	    (subvarIndex == (int) vd->Subvar_Index) &&
	    (mn_min > vd->MatID)) {
	  var_rec = i;
	  mn_min = vd->MatID;
	}
      } 
    }
    if (var_rec == -1) {
      EH(-1,"AUGH! we shouln't be here");
    }
    rec_used[var_rec] = 1;
      
    /*
     * Now fill in fields in the nodal vars struct
     */
    add_var_to_nv_struct(nv_old->Var_Desc_List[var_rec],
			 nv_new);
  }
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

NODAL_VARS_STRUCT *
copy_and_rearrange_nv(NODAL_VARS_STRUCT *nv_old)

    /*****************************************************************
     *
     * copy_and_rearrange_nv():
     *
     *   This function copies an old nodal variable structure into
     *   a new one, rearranging the order of the variables as it
     *   goes. 
     *   This function sets the order in the solution vector. Right
     *   now, the order of the solution vector is set to its original
     *   GOMA one:
     *      Outer to Inner:
     *            var_type
     *              sub_var_type
     *                 matID
     *   This is not necessarily how it should be for optimal
     *   speed, however. In particular the species unknowns should
     *   be grouped together. Therefore, rearrangement of the
     *   solution vector  might occur in the future.
     *****************************************************************/
{
  int var_type, subvarIndex = 0;
  int *rec_used;
  NODAL_VARS_STRUCT *nv_new = alloc_struct_1(NODAL_VARS_STRUCT, 1);
  nv_new->list_index = nv_old->list_index;

  /*
   *  Note: It is permissible to have nodes with zero unknowns defined
   *        on them. For example, if the element is quadratic and the
   *        variable interpolations are all Q1, the extra nodes in each
   *        element will not have any unknowns assigned to them. This
   *        is a waste, but it is permissible.
   */
  if (nv_old->Num_Var_Desc == 0) {
    return nv_new;
  }
  
  rec_used = alloc_int_1(nv_old->Num_Var_Desc, 0);
  for (var_type = V_FIRST; var_type < V_LAST; var_type++) {
    if (var_type == MASS_FRACTION) {
      for (subvarIndex = 0; subvarIndex < upd->Max_Num_Species_Eqn;
	   subvarIndex++) {
        rearrange_mat_increasing(var_type, subvarIndex, nv_old, nv_new,
				 rec_used);
      }
    } else {
      rearrange_mat_increasing(var_type, 0, nv_old, nv_new, rec_used);
    }
  }

  safer_free((void **) &rec_used);
  return nv_new;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void
nodal_vars_free(NODAL_VARS_STRUCT *nv)

    /*************************************************************************
     *
     * nodal_vars_free():
     *
     *  Frees memory allocated underneath the nodal_vars
     *  structure, and returns the structure to its default condition.
     *  This is the reverse of the _alloc() operation. Note, there is
     *  no _alloc() procedure for this structure. It is actually populated
     *  by repeated calls to add_var_to_nv_struct(). Thus, this routine
     *  serves to undo whatever was malloced by repeated calls to
     *  add_var_to_nv_struct().
     ************************************************************************/   
{
  int v;
  for (v = 0; v < V_LAST; v++) {
    safer_free((void **) &(nv->Var_Type_Index[v]));
  }
  safer_free((void **) &(nv->Var_Desc_List));
  safer_free((void **) &(nv->Nodal_Offset));
  (void) memset((void *)nv, 0, sizeof(NODAL_VARS_STRUCT));
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void
nodal_vars_destroy(NODAL_VARS_STRUCT **nv_hdl)

    /*************************************************************************
     *
     * nodal_vars_destroy():
     *
     * Free memory required by the nodal_vars structure.
     * This routine basically unmallocs whatever nodal_vars_alloc
     * malloced previously, and then unmallocs the variable
     * description structure itself.
     *
     *  Note, this is the opposite of the _create routine.
     *************************************************************************/
{
  nodal_vars_free(*nv_hdl);
  safer_free((void **) nv_hdl);
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int
find_var_type_index_in_nv(const int var_type, NODAL_VARS_STRUCT *nv)

    /*************************************************************************
     *
     * find_var_type_index_in_nv():
     *
     *    This routine returns the index of the first variable struct in the 
     *    nodal variable structure that matches var_type input from the
     *    argument list. Note there may be more than one variable description
     *    structure with the same variable type in a single nodal_vars
     *    structure.
     *    If the var_type doesn't exist in the nodal variable structure,
     *    this function returns -1.
     *************************************************************************/
{
  if (! (nv->Var_Type_Index[var_type])) return -1;
  return (nv->Var_Type_Index[var_type][0]);
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int
find_vd_index_in_nv(const int var_type, const int matID, NODAL_VARS_STRUCT *nv)

    /*************************************************************************
     *
     * find_vd_index_in_nv():
     *
     *    This routine returns the index of the matching variable struct in the 
     *    nodal variable structure that matches var_type and matID input from
     *    the argument list. Note there may be more than one variable
     *    description structure with the same variable type in a single
     *    nodal_vars structure.
     *     Generic variables are denoted by matID = -1.
     *         This function will allow a match for a nongeneric matID argument
     *    with a generic (matID = -1) variable description structure.
     *
     *          matID_arg    vd_matID's      matches
     *       ---------------------------------------------
     *           -1            -1, 2         first
     *           -1             2  3         no_match
     *            2            -1, 2         second
     *            2            -1, 1         first
     *  
     *    If no match is found in the nodal variable structure,
     *    this function returns -1.
     *************************************************************************/
{
  int i, ifound = -1;
  VARIABLE_DESCRIPTION_STRUCT *vd;
  int *list = nv->Var_Type_Index[var_type];
  if (! list) return -1;
  for (i = 0; i < nv->Num_Var_Desc_Per_Type[var_type]; i++) {
    vd = nv->Var_Desc_List[list[i]];
    if (matID == -1) {
      if (vd->MatID == -1) return i;
    } else {
      if (vd->MatID == -1) ifound = i;
      else if (vd->MatID == matID) return i;
    }
  }
  return ifound;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void
pack_nv (struct nv_packed *nvp, NODAL_VARS_STRUCT *nv)
    
    /*************************************************************************
     *
     * find_var_type_index_in_nv():
     *
     *    Pack the information contained in nv into an nv_packed structure
     *************************************************************************/
{
  int i;
  struct vd_packed *vdp;
  VARIABLE_DESCRIPTION_STRUCT *vd;
  for (i = 0; i < nv->Num_Var_Desc; i++) {
    vd = nv->Var_Desc_List[i];
    nvp->num = nv->Num_Var_Desc;
    vdp = &(nvp->VDP_i);
    vdp->Variable_Type = vd->Variable_Type;
    vdp->Ndof = vd->Ndof;
    vdp->Subvar_Index = vd->Subvar_Index;
    vdp->MatID = vd->MatID;
    nvp++;
  }
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

NODAL_VARS_STRUCT *
unpack_nv (struct nv_packed *nvp, const int imtrx)
    
    /*************************************************************************
     *
     * unpack_nv():
     *
     *    Pack the information contained into an nv_packed structure into
     * variable description structures and a new nodal var structure
     *************************************************************************/
{
  int i, numVD;
  struct vd_packed *vdp;
  VARIABLE_DESCRIPTION_STRUCT *var_ptr;
  NODAL_VARS_STRUCT *nv = alloc_struct_1(NODAL_VARS_STRUCT, 1);
  numVD = (int) nvp->num;
  for (i = 0; i < numVD; i++) {
    vdp = &(nvp->VDP_i);

    /*
     * Find a variable structure that matches this one 
     */ 
    var_ptr = find_or_create_vd((int) vdp->Variable_Type, (int) vdp->Ndof,
				(int) vdp->MatID, (int) vdp->Subvar_Index, imtrx);
    
    /*
     * Add the variable definition structure into the nodal variable
     * structure
     */
    add_var_to_nv_struct(var_ptr, nv);

    nvp++;
  }
  /*
   * Return the new nodal variables structure
   */
  return nv;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
