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
 *$Id: mm_fill_ptrs.c,v 5.16 2010-04-07 22:27:00 prschun Exp $
 */


/* Standard include files */

#include <stdlib.h>
#include <stdio.h>

/* GOMA include files */

#include "mm_fill_ptrs.h"

#include "std.h"
#include "rf_fem_const.h"
#include "rf_fem.h"
#include "el_elm.h"
#include "el_geom.h"
#include "rf_vars_const.h"
#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "mm_as.h"
#include "mm_mp.h"
#include "mm_mp_structs.h"
#include "dpi.h"
#include "el_elm_info.h"
#include "exo_struct.h"
#include "mm_elem_block_structs.h"
#include "mm_fill_ls.h"
#include "mm_fill_stress.h"
#include "mm_mp_const.h"
#include "mm_shell_util.h"
#include "mm_unknown_map.h"
#include "rd_mesh.h"
#include "rf_node_const.h"
#include "sl_util.h"
#include "sl_util_structs.h"
#include "stdbool.h"

extern	dbl *p0;		/* Defined in mm_as_alloc.c */


#include "mm_eh.h"

#define GOMA_MM_FILL_PTRS_C

/***********************************************************************/
#ifdef DEBUG_NOT
static void print_ei_index_desc(void);
#endif
/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
/*
 * Externals for the solution vector and its derivative at the
 * current and old times. These are filled in every time an element
 * dofptrs are set up.
 */
double *x_static = NULL;
double *x_old_static = NULL;
double *xdot_static = NULL;
double *xdot_old_static = NULL;
double *x_dbl_dot_static = NULL;
double *x_dbl_dot_old_static = NULL;
double *x_pred_static = NULL;

/***********************************************************************/
/***********************************************************************/
/***********************************************************************/

static int
offset_matid_vt(const int lvdof, const int mn_new,
		NODAL_VARS_STRUCT *nv, const int mn_curr,
		const int varType, const int subvarType)

    /*
     *
     */
{
  int index_new, index_old;
  index_new = get_nv_index_varType(nv, varType, mn_new, subvarType);
  index_old = get_nv_index_varType(nv, varType, mn_curr, subvarType);
  return (lvdof + index_new - index_old);
}
/***************************************************************************/

static void 
linkParent_ei(const int elem, const int relatedElem, const int v, 
	      const Exo_DB *exo, int imtrx)
{
  int i;
  /*
   *
   *
   */
  int slot = -1;
  int mn = ei[imtrx]->mn;
  //int ebn = ei_ptr->elem_blk_index;
  bool alreadyDone = 0;
  for (i = 0; i < MAX_ELEMENT_INDICES_RELATED; i++) {
    if (ei[imtrx]->linkedEIelems[i] == relatedElem) {
      slot = i;
      alreadyDone = 1;
      break;
    }
  }
  if (slot == -1) 
    {
      for (i = 0; i < MAX_ELEMENT_INDICES_RELATED; i++) {
	if (ei[imtrx]->linkedEIelems[i] == -1) {
	  slot = i;
	  break;
	}
      }
      if (slot == -1) {
	fprintf(stderr, "Ran out of slots in \n");
	exit(-1);
      }
      ei[imtrx]->linkedEIelems[i] = relatedElem;
    }
  /*
   *  We check this condition here. But, I think it should already be
   *  redundant.
   */
  if (ei[imtrx]->owningElementForColVar[v] != relatedElem)
    {
      fprintf(stderr, "inconsistent entry in ei[imtrx]->owningElementForColVar[%d], %d, %d\n",
	      v, ei[imtrx]->owningElementForColVar[v], relatedElem);
      EH(-1, "linkedParent_ei() error\n");
    }


  /*
   *  Link the current element's ei[imtrx]->owningElement_ei_ptr[v] to a slot
   *  in the eiRelated vector of pointers. What we are saying is that
   *  the column unknowns for this variable will be handled by eiRelated[slot] 
   */
  ei[imtrx]->owningElement_ei_ptr[v] = eiRelated[slot];


  /*
   * Ok, if we haven't already done so, fill up the slave ei structure.
   */
  if (!alreadyDone) 
    {
      load_ei(relatedElem, exo, ei[imtrx]->owningElement_ei_ptr[v], imtrx);

      pd = pd_glob[mn];
      Current_EB_ptr = Element_Blocks + ei[imtrx]->elem_blk_index;
      pd	 = pd_glob[mn];      /* Problem Description */
      mp	 = mp_glob[mn];      /* Material Properties */
      cr	 = cr_glob[mn];      /* Constitutive Relations */
      elc	 = elc_glob[mn];     /* Elastic Constitutive */
      elc_rs = elc_rs_glob[mn];  /* Elastic Constitutive */
      gn	 = gn_glob[mn];      /* Generalized Newtonian */
      vn	 = vn_glob[mn];      /* Viscoelastic Nonmodal */
      evpl	 = evpl_glob[mn]; 
    }

}
/***************************************************************************/

int
load_ei(const int elem, const Exo_DB *exo, struct Element_Indices *ei_ptr_fill, int imtrx)

     /***********************************************************************
      * load_ei -- load data members of the Element_Indeces structure
      * 
      * input 
      * ------
      *   ei   -- Pointer to the Element_Indeces structure
      *           (this is a variable with global scope)
      *   elem -- current element number for the element to be 
      *           processed.
      *   exo  -- pointer to the exodus data base
      *
      *
      *  Output
      * ---------
      *   The Element_Indeces structure will be completely filled up
      *   with pertinent information for the current element at
      *   the end of this routine.
      ***********************************************************************/
{
  int mn, v, etype, nnodes, lvdof, ledof, ln, gnn, owned, enunks;
  int i, nunks, index, ivarDes, indof, good_index, gnn_Parent;
  int iLvdesc; /* Running index for the number of local variable 
		* description indeces found in the current element.
		*/

  int MFdofSVNE0, j_vdesc;
  int n_eb, n_mn, mn_retn;
  static int first_time = TRUE;
  int not_matched, iLv_found, i_vd, i_vdesc, iLv_start, lvdesc;
  VARIABLE_DESCRIPTION_STRUCT *vd0;
  NODAL_VARS_STRUCT *nv;
  VARIABLE_DESCRIPTION_STRUCT *vd;
  NODE_INFO_STRUCT *ni;
  VCRR_STRUCT **vcrr;
  bool ImAChild, ImAParent;

  int n_eb_Parent, elem_Parent, n_mn_Parent, n_ielem_type_Parent;
  int iconnect_ptr_Parent;
  int  num_local_nodes_Parent;
  struct Element_Indices *ei_ptr;

  int mode = 0;
  ei_ptr = ei[imtrx];
  if (ei_ptr_fill != 0) 
    {
      mode = 1;
      ei_ptr = ei_ptr_fill;
    }


  /*
   * Clear the list of related elements whose Element_Indices structs
   * are calculated via a linked list of pointers.
   */
  if (mode == 0) {
    for (i = 0; i < MAX_ELEMENT_INDICES_RELATED; i++) {
      ei_ptr->linkedEIelems[i] = -1;
    }
  }

  /*
   * Look up the element block index for the current element
   */

  ei_ptr->elem_blk_index = find_elemblock_index(elem, exo);

  /*
   * Store the  variable, current_EB_ptr, which points to the
   * current element block structure
   */
  ei_ptr->current_EB_ptr = Element_Blocks + ei_ptr->elem_blk_index;

  /*
   * Obtain the material index from the previously calculated
   * element block to material index array, Matilda
   */
  mn  = Matilda[ei_ptr->elem_blk_index];
  if (mn < 0) {
    fprintf(stderr, "Negative material index for eb block %d \n",
            ei_ptr->elem_blk_index);
    EH(-1, "Ignore element block logic error");
  }

  /*   
   *  Assign the Globally-Scoped pointers,  pd, 
   *  mp, cr, elc, and gn, to the appropriate material index for
   *  the current element
   */
  if (mode == 0) {
    Current_EB_ptr = Element_Blocks + ei_ptr->elem_blk_index;
    pd	 = pd_glob[mn];      /* Problem Description */
    mp	 = mp_glob[mn];      /* Material Properties */
    cr	 = cr_glob[mn];      /* Constitutive Relations */
    elc	 = elc_glob[mn];     /* Elastic Constitutive */
    elc_rs = elc_rs_glob[mn];  /* Elastic Constitutive */
    gn	 = gn_glob[mn];      /* Generalized Newtonian */
    vn	 = vn_glob[mn];      /* Viscoelastic Nonmodal */
    evpl	 = evpl_glob[mn];    /* Viscplastic Constitutive */
  } else {
    // set to zero so that it will core dump and we can see what's going on
    pd	 = 0;
    Current_EB_ptr = 0;
    mp = 0;
    cr = 0;
    elc = 0;
    elc_rs = 0;
    gn = 0;
    vn = 0;
    evpl = 0;
  }
  PROBLEM_DESCRIPTION_STRUCT    *pd_ptr = pd_glob[mn];

  /*
   *  Assign values to the scalar members of ei. 
   */
  ei_ptr->mn	         = mn;
  ei_ptr->ielem	         = elem;
  ei_ptr->ielem_type     = Elem_Type(exo, elem);
  ei_ptr->ielem_shape    = type2shape(ei_ptr->ielem_type);
  ei_ptr->elem_blk_index = find_elemblock_index(elem, exo);
  ei_ptr->elem_blk_id    = exo->eb_id[ei_ptr->elem_blk_index];
  ei_ptr->ielem_dim	 = elem_info(NDIM, ei_ptr->ielem_type);
  /*
   *   Store the pointer to the beginning of this element's
   *   connectivity list
   */
  ei_ptr->iconnect_ptr    = Proc_Connect_Ptr[elem];
  ei_ptr->num_local_nodes = elem_info(NNODES, ei_ptr->ielem_type);
  ei_ptr->num_sides	  = shape2sides(ei_ptr->ielem_shape);

  /*
   *  ei_ptr->deforming_mesh is TRUE if there are mesh equations on
   *  this block OR on a parent block (if this is a shell element)
   */
  ei_ptr->deforming_mesh = FALSE;
  if (pd_ptr->e[imtrx][R_MESH1])
    {
      ei_ptr->deforming_mesh = TRUE;
    }
  else if (is_shell_element(ei_ptr->ielem, exo))
    {
      for (i = 0; i < num_elem_friends[ei_ptr->ielem]; i++)
        {
          n_eb = find_elemblock_index(elem_friends[ei_ptr->ielem][i], exo);
          n_mn = Matilda[n_eb];
          if (pd_glob[n_mn]->e[imtrx][R_MESH1]) 
            {
              ei_ptr->deforming_mesh = TRUE;
            }
        }
    }
  
  etype  = ei_ptr->ielem_type;
  nnodes = ei_ptr->num_local_nodes;

  /*
   * Initialize...arrays
   *
   *  HKM -> Later, we need to figure out whether we can completely
   *         cut this section.
   */
  if (first_time)
    {
      first_time = FALSE;
      for (v = V_FIRST; v < V_LAST; v++)
        {
          ei_ptr->Num_Lvdesc_Per_Var_Type[v] = 0;
          ei_ptr->Lvdesc_First_Var_Type[v]   = -1;
          ei_ptr->dof[v] = 0;
          for (lvdof = 0; lvdof < MDE; lvdof++)
            {
              ei_ptr->Baby_Dolphin[v][lvdof]    = -1;
              ei_ptr->dof_list[v][lvdof]        = -1;
              ei_ptr->gun_list[v][lvdof]        = -1;
              ei_ptr->gnn_list[v][lvdof]        = -1;
              ei_ptr->ln_to_dof[v][lvdof]       = -1;
              ei_ptr->ln_to_first_dof[v][lvdof] = -1;
            }
        }
    } 
  else
    {
      for (v = V_FIRST; v < V_LAST; v++) {
        if (Num_Var_In_Type[imtrx][v])
          {
            ei_ptr->dof[v] = 0;
            for (lvdof = 0; lvdof < MDE; lvdof++) 
              {
                ei_ptr->Baby_Dolphin[v][lvdof]    = -1;
                ei_ptr->dof_list[v][lvdof]        = -1;
                ei_ptr->gun_list[v][lvdof]        = -1;
                ei_ptr->gnn_list[v][lvdof]        = -1;
                ei_ptr->ln_to_dof[v][lvdof]       = -1;
                ei_ptr->ln_to_first_dof[v][lvdof] = -1;
              }
          }
      }
    }
  
  for (v = 0; v < Num_Var_Info_Records; v++) 
    {
      ei_ptr->VDindex_to_Lvdesc[v] = -1;
    } 

  /*
   *  Set the local element dof index to zero. We will
   *  index this counter for each dof we find below to keep a running
   *  total
   */
  ledof = 0;

  /*
   *  Set the local variable description index, iLvdesc, to zero. We will
   *  Keep a running track of this index.
   */
  iLvdesc = 0;
  /*
   * Loop over variable types -> this will set how the local variable
   * description indeces will be ordered within the element.
   */
  for (v = V_FIRST; v < V_LAST; v++)
    {
      ei_ptr->owningElementForColVar[v] = -1;
      ei_ptr->owningElement_ei_ptr[v] = NULL;
   
      if (Num_Var_In_Type[imtrx][v])
        {
          iLv_start = iLvdesc;
          iLv_found = 0;
          for (ln = 0; ln < nnodes; ln++)
            {
              gnn = Proc_Elem_Connect[ei_ptr->iconnect_ptr + ln];
              ni = Nodes[gnn];
              nv = ni->Nodal_Vars_Info[imtrx];
              for (i_vd = 0; i_vd < nv->Num_Var_Desc_Per_Type[v]; i_vd++)
                {
                  index = nv->Var_Type_Index[v][i_vd];
                  vd = nv->Var_Desc_List[index];
                  not_matched = TRUE;
                  /*
                   * Try to find a match amongst the already found
                   * variable description structures
                   */
                  for (i_vdesc = 0; i_vdesc < iLv_found; i_vdesc++)
                    {
                      vd0 = ei_ptr->Lvdesc_vd_ptr[iLv_start + i_vdesc];
                      if (vd == vd0) 
                        {
                          not_matched = FALSE; break;
                        }
                    }
                  /*
                   * If we have found a new local variable description
                   * then add it to the list
                   */
                  if (not_matched) 
                    {
                      ei_ptr->Lvdesc_vd_ptr[iLvdesc] = vd;
                      /*
                       *  Make sure the first vdesc found is the one
                       *  pertaining to this element. We have problems
                       *  with the logic below for mass fraction variables
                       *  with subvar index greater than zero. However,
                       *  let's go on.
                       */
                      if (iLv_found > 0 && vd->Subvar_Index == 0) 
                        {
                          vd0 = ei_ptr->Lvdesc_vd_ptr[iLv_start];
                          if (vd0->MatID != mn)
                            {
                              if ((vd->MatID == mn) ||
                                  (vd->MatID == -1 && vd0->MatID != -1))
                                {
                                  ei_ptr->Lvdesc_vd_ptr[iLv_start] = vd;
                                  ei_ptr->Lvdesc_vd_ptr[iLvdesc] = vd0;
                                }
                            }
                        }
                      iLv_found++;
                      iLvdesc++;
                    }
                }
            } /* End of loop over the local nodes */

          /*
           * Establish a pointer to the first variable description
           * index corresponding to this variable type
           */
          if (iLv_found > 0) 
            {
              ei_ptr->Lvdesc_First_Var_Type[v] = iLv_start;
            } 
          else
            {
              ei_ptr->Lvdesc_First_Var_Type[v] = -1;
            }
          /*
           * Store the total number of variable descriptions found
           * for this variable type
           */
          ei_ptr->Num_Lvdesc_Per_Var_Type[v] = iLv_found;

          /*
           * Loop over the variable descriptions found for this
           * variable type, filling in various mappings.
           */
          for (i_vdesc = 0; i_vdesc < iLv_found; i_vdesc++) 
            {
              i = iLv_start + i_vdesc;
              vd = ei_ptr->Lvdesc_vd_ptr[i];
              /*
               *  Fill in the mapping between the variable
               *  description list index and the local variable
               *  description index.
               */
              ei_ptr->VDindex_to_Lvdesc[vd->List_Index] = i;
              /*
               * Fill in the variable_type field
               */
              ei_ptr->Lvdesc_to_Var_Type[i] = v;
              /*
               * Fill in the material Id field
               */
              ei_ptr->Lvdesc_to_MatID[i] = vd->MatID;
              /*
               * Initialize the number of degrees of freedom 
               * corresponding to this variable description
               */
              ei_ptr->Lvdesc_Numdof[i] = 0;
              /*
               * for Mass fraction unknowns we need to store the
               * the species index
               */
              ei_ptr->Lvdesc_to_MFSubvar[i] = vd->Subvar_Index;
              /*
               * Fill in the local offset array
               */
              for (ln = 0; ln < nnodes; ln++) 
                {
                  gnn = Proc_Elem_Connect[ei_ptr->iconnect_ptr + ln];
                  ni = Nodes[gnn];
                  nv = ni->Nodal_Vars_Info[imtrx];
                  ei_ptr->Lvdesc_Lnn_to_Offset[i][ln] = 
                    get_nodal_unknown_offset(nv, v, vd->MatID, 
                                             vd->Subvar_Index, NULL);
                }
            }
        }
    }
  /*
   * Fill in the total number of variable descriptions found
   * for this element
   */
  ei_ptr->Num_Lvdesc = iLvdesc;

  /*
   * Need to initialize the Lvdesc_Lnn_to_lvdof array
   * to make sure that values are -1. We may not get to
   * all locations via the logic in the loops below
   */
  for (i_vdesc = 0; i_vdesc < ei_ptr->Num_Lvdesc; i_vdesc++) {
    for (ln = 0; ln < nnodes; ln++) {
      ei_ptr->Lvdesc_Lnn_to_lvdof[i_vdesc][ln] = -1;
      ei_ptr->Lvdesc_Lnn_Numdof[i_vdesc][ln] = 0;
    }
  }

  /*
   *  Loop over the variable types defined in the problem
   */
  for (v = V_FIRST; v < V_LAST; v++)
    {
      ei_ptr->owningElementForColVar[v] = -1;
      if (Num_Var_In_Type[imtrx][v])
        {
          lvdof = 0;
          for (ln = 0; ln < nnodes; ln++)
            {
              /*
               * For this local node "ln" what is the global node number, "gnn"?
               */
              gnn = Proc_Elem_Connect[ei_ptr->iconnect_ptr + ln];

              /*
               * store pointers to the Node Info structure and the Node_Vars
               * Info structure for the current global node.
               */
              ni = Nodes[gnn];
              nv = ni->Nodal_Vars_Info[imtrx];
	
              /*
               *  Is this global node number owned by the processor?
               */
              owned = (gnn < DPI_ptr->num_owned_nodes);
	      
              /*
               * For this variable at this local node, how many dofs are 
               * there?
               * Note this can be zero...
               *
               * WARNING..you will need to account for multiple species k here...
               *
               * Warning..  node_info() isn't correct for elements next to
               *        materials which have extra equations assigned to them
               *        than the ones in the current material. They will
               *        report no degrees of freedom at the element, while
               *        in fact, there are degrees of freedom at that node
               *        corresponding to variables, which are not interpolated
               *        in the current element.
               *          The new system catches these degrees of freedom,
               *        and assigns them to be inactive.
               */
              nunks = get_nv_ndofs_modMF(nv, v);

              /*
               * For this variable at this local node, how many dofs are
               * part of a valid element interpolation?
               */
              enunks = nunks;
#ifdef DEBUG_HKM
              if (elem == 165) {
                if (v == 0) {
                  // if (gnn == 319) {
                  //printf("we are here2\n");
                  //}
                }
              }
#endif
              if (nunks) 
                {
                  /*
                   *  The return value, enunks, is equal to the number of degrees 
                   *  of freedom at a local node, ln, at global node, gnn,
                   *  for a variable type, v,
                   *  corresponding to the current
                   *  element with element type, etype, given the problem
                   *  description structure, pd. So, we just pick up unknowns
                   *  which are actively interpolated as unknowns on this
                   *  element.
                   */
                  enunks = dof_lnode_var_type(ln, etype, gnn, v, pd_ptr, imtrx);
                }
	      
              if (enunks) 
                {
#ifdef DEBUG_HKM
                  if (ei_ptr->ielem == 165) {
                    //  if ( ei_ptr->ielem != elem) {
                    // printf("we are here \n");
                    // }
                  }
#endif
                  ei_ptr->owningElementForColVar[v] = elem;
                }
              /*
               * For element variables that masquerade as nodal variables
               * assigned to local node 0 in the element (P0 and P1
               * interpolations), we need to correct the calculation above.
               * nunks will include degrees of freedom corresponding to
               * local node numbers greater than 0. These dof's are really
               * element variable dofs for neighboring elements and don't
               * belong in this element's list of dofs under any circumstance.
               * Therefore, they are taken out of the lvdof mapping for the
               * element. Note, I believe discontinuous galerkin
               * interpolations will have the same treatment in some cases.
               * For other cases, the dofs are located at the centroid
               * node of an element, so a special case is not needed.
               * However, I
               * have not included a special case for discontinuous
               * galerkin interpolations, here, yet.
               */
              if (pd_ptr->i[imtrx][v] == I_P0 || pd_ptr->i[imtrx][v] == I_P1 ||
                  pd_ptr->i[imtrx][v] == I_P0_G || pd_ptr->i[imtrx][v] == I_P1_G ||
                  pd_ptr->i[imtrx][v] == I_P0_GP || pd_ptr->i[imtrx][v] == I_P1_GP ||
                  pd_ptr->i[imtrx][v] == I_P0_GN || pd_ptr->i[imtrx][v] == I_P1_GN ||
                  pd_ptr->i[imtrx][v] == I_P0_XV || pd_ptr->i[imtrx][v] == I_P1_XV ||
                  pd_ptr->i[imtrx][v] == I_P1_XG)
                {
                  nunks = enunks;
                }
	
              /*
               *  We need the if statement below (if nunks > 0) because
               *  of the possibility that nunks and nv->Num_Var_Desc_Per_Type[v]
               *  may have become uncorreled due to the nunks = enunks line
               *  special case for the I_P0 and I_P1 interpolations.
               *  There are others ways, such as inspecting the variable
               *  description structure in the loop below, to achieve the
               *  same ends. However, the present one seems to work.
               */
              if (nunks > 0)
                {
                  /*
                   *  Increase the count of the number of degrees of freedom
                   *  for the variable in the element, ei_ptr->dof[]
                   */
                  ei_ptr->dof[v] += nunks;
                  i = 0;
                  for (ivarDes = 0; ivarDes < (int) nv->Num_Var_Desc_Per_Type[v];
                       ivarDes++) 
                    {
                      /*
                       *  Find the variable description corresponding
                       *  to this lvdof
                       */
                      index = (int) nv->Var_Type_Index[v][ivarDes];
                      vd = nv->Var_Desc_List[index];

                      /*
                       *  Find the lvdesc index corresponding to this
                       *  variable description 
                       */
                      iLv_start = ei_ptr->Lvdesc_First_Var_Type[v];
                      not_matched = TRUE;
                      for (i_vd = 0; 
                           i_vd < ei_ptr->Num_Lvdesc_Per_Var_Type[v] && not_matched; 
                           i_vd++) 
                        {
                          i_vdesc = iLv_start + i_vd;
                          vd0 = ei_ptr->Lvdesc_vd_ptr[i_vdesc];
                          if (vd == vd0)
                            {
                              not_matched = FALSE;
                            }
                        }
                      if (not_matched) 
                        {
                          EH(-1, "ERROR");
                        }

                      /*
                       * lvdof counters are not incremented for mass
                       * fraction unknowns whose subvar index is greater
                       * than 0. However, all other counters are. Thus
                       * we create a special boolean here to be used
                       * below.
                       */
                      if (vd->Variable_Type == MASS_FRACTION &&
                          vd->Subvar_Index != 0)
                        {
                          MFdofSVNE0 = TRUE;
                        } 
                      else
                        {
                          MFdofSVNE0 = FALSE;
                        }

                      /*
                       * Store a mapping between the (local variable 
                       * description number, local node number) pair
                       * to the local variable degree of freedom index.
                       * We do this above the next loop, because we
                       * want to store the value of the first lvdof
                       * corresponding to an i_vdesc at a local node.
                       */
                      if (MFdofSVNE0) 
                        {
                          for (j_vdesc = 0; j_vdesc < ei_ptr->Num_Lvdesc; j_vdesc++) 
                            {
                              vd0 = ei_ptr->Lvdesc_vd_ptr[j_vdesc];
                              if (vd0->Variable_Type == MASS_FRACTION &&
                                  vd0->Subvar_Index == 0 &&
                                  vd0->MatID == vd->MatID)
                                {
                                  break;
                                }
                            }
                          ei_ptr->Lvdesc_Lnn_to_lvdof[i_vdesc][ln] = ei_ptr->Lvdesc_Lnn_to_lvdof[j_vdesc][ln];
                        } 
                      else
                        {
                          ei_ptr->Lvdesc_Lnn_to_lvdof[i_vdesc][ln] = lvdof;
                        }

                      /*
                       * Loop over Ndof in a vd. Note, only centroid
                       * pressures currently have more than one Ndof
                       * in a variable description structure
                       */
                      for (indof = 0; indof < vd->Ndof; indof++) {
                        ei_ptr->dof_list[v][lvdof] = ln;
                        ei_ptr->gnn_list[v][lvdof] = gnn;
                        lvdesc = ei_ptr->Lvdesc_Numdof[i_vdesc];
                        ei_ptr->Lvdesc_to_Lnn[i_vdesc][lvdesc] = ln;
                        ei_ptr->Lvdesc_to_Gnn[i_vdesc][lvdesc] = gnn;
                        ei_ptr->Lvdesc_to_ledof[i_vdesc][lvdesc] = ledof;
                        ei_ptr->owned_ledof[ledof] = owned;
                        ei_ptr->ieqn_ledof[ledof] = (ni->First_Unknown[imtrx] +
                                                     nv->Nodal_Offset[index] + indof);
                        ei_ptr->Lvdesc_to_Gun[i_vdesc][lvdesc] = ei_ptr->ieqn_ledof[ledof];
                        /*
                         * lvdof indexing skips subvariables 
                         * greater than zero in its indexing scheme.
                         * It assumes an equal number of species per domain
                         * with the following indexing scheme.
                         *  base  +  SubVar*Dolphin[imtrx][I][e]
                         * to point to higher SubVar's from the base
                         * MF unknown corresponding to that phase.
                         * (note this lvdof indexing scheme breaks down
                         *  when the number of species is different in each
                         *  domain).
                         * Thus, we point to the sub variable = 0
                         * entries in the arrays below that depend on
                         * lvdof as an array component and we do not overwrite
                         * these entries with info about high subvar entries.
                         */
                        if (!MFdofSVNE0) 
                          {
                            ei_ptr->lvdof_to_lvdesc[v][lvdof] = i_vdesc;
                            ei_ptr->lvdof_to_lvdesc_dof[v][lvdof] = lvdesc;
                            ei_ptr->lvdof_to_ledof[v][lvdof] = ledof;
                            ei_ptr->gun_list[v][lvdof] = ei_ptr->ieqn_ledof[ledof];
                          }

                        if (vd->MatID == -1) 
                          ei_ptr->matID_ledof[ledof] = mn;
                        else              
                          ei_ptr->matID_ledof[ledof] = vd->MatID;


                        /*
                         * Find the local variable type dof that matches the
                         * current unknown. This would be easy, except for the case
                         * of mass fractions with subvariable types greater than 0.
                         * In this case, the local variable type dof is degenerate.
                         * One must look up the subvariable type zero case, done
                         * previously in this same loop. Thus, we find the
                         * lvdof for the subvar = 0 case for the correct matID, in
                         * the case of discontinuous variables.
                         */
                        if (MFdofSVNE0) 
                          {
                            for (j_vdesc = 0; j_vdesc < ei_ptr->Num_Lvdesc; j_vdesc++) 
                              {
                                vd0 = ei_ptr->Lvdesc_vd_ptr[j_vdesc];
                                if (vd0->Variable_Type == MASS_FRACTION &&
                                    vd0->Subvar_Index == 0 &&
                                    vd0->MatID    == vd->MatID) {
                                  break;
                                }
                              }
                            ei_ptr->Lvdesc_to_lvdof[i_vdesc][lvdesc] =
                              ei_ptr->Lvdesc_Lnn_to_lvdof[j_vdesc][ln];
                          } 
                        else
                          {
                            ei_ptr->Lvdesc_to_lvdof[i_vdesc][lvdesc] = lvdof;
                          }
                        /*
                         * Decide whether the current dof is part of an active
                         * interpolation of the variable type within this
                         * element
                         */
                        ei_ptr->active_interp_ledof[ledof] = 1;
                        if (enunks == 0) 
                          {
                            ei_ptr->active_interp_ledof[ledof] = 0;
                          } 
                        else if (nv->Num_Var_Desc_Per_Type[v] > 1) 
                          {
                            good_index = get_nv_index_varType(nv, v, mn, vd->Subvar_Index);
                            if (good_index != index)
                              {
                                ei_ptr->active_interp_ledof[ledof] = 0;
                              }
                          }
	    
                        /*
                         * If this dof isn't part of an active degree of freedom
                         * then set it's material ID index back to -1.
                         */
                        if (!ei_ptr->active_interp_ledof[ledof]) 
                          {
                            if (ei_ptr->matID_ledof[ledof] >= 0) 
                              {
                                ei_ptr->matID_ledof[ledof] = vd->MatID;
                              }
                          }

                        if (!MFdofSVNE0) 
                          {
                            ei_ptr->Baby_Dolphin[v][lvdof] = i;
	      

                            /*
                             * Construct mappings between (variable type, local node
                             * number) and the local variable dof. Because there may
                             * be more than one dof with the same variable type at
                             * a node, this is a very imperfect vehicle. A better
                             * approach would be to specify the material id of the
                             * dof along with variable Type and node number. This
                             * is readily available in the new data structures.
                             */
                            if (ei_ptr->active_interp_ledof[ledof]) 
                              {
                                ei_ptr->ln_to_dof[v][ln] = lvdof;
                              } 
                            else
                              {
                                if (ei_ptr->ln_to_dof[v][ln] == -1)
                                  {
                                    ei_ptr->ln_to_dof[v][ln] = lvdof;
                                  }
                              }
                            if (i == 0) 
                              {
                                ei_ptr->ln_to_first_dof[v][ln] = lvdof;
                              }
                          }

                        if (!MFdofSVNE0)
                          {
                            /*
                             * Look up to see if there are any volume contribution
                             * row redirection for this dof
                             */
                            if (!ei_ptr->active_interp_ledof[ledof])
                              {
                                ei_ptr->lvdof_to_row_lvdof[v][lvdof] = -1;
                              } 
                            else 
                              {
                                ei_ptr->lvdof_to_row_lvdof[v][lvdof] = lvdof;
                                if (ni->Resid_Wksp)
                                  {
                                    vcrr = ni->Resid_Wksp->Vcrr;
                                    mn_retn = vcrr_lookup(vcrr, v, ei_ptr->mn);
                                    if (mn_retn != ei_ptr->mn)
                                      {
                                        /*
                                         * search the unknowns at this node with var type v
                                         * to get the offset
                                         */
                                        ei_ptr->lvdof_to_row_lvdof[v][lvdof] = 
                                          offset_matid_vt(lvdof, mn_retn, nv, ei_ptr->mn,
                                                          v, vd->Subvar_Index);
                                      }
                                  }
                              }	
                          }  
		    
                        /*
                         * Increment counters
                         *  Note: we do not increment the counter if we are
                         *        currently processing a mass fraction unknown
                         *        and the subvariable index is greater than one.
                         *        Multiple subvariable dofs are degenerate under
                         *        the lvdof nomenclature.
                         */
                        if (!MFdofSVNE0)
                          {
                            lvdof++;
                            i++;
                          }
                        ledof++;
                        ei_ptr->Lvdesc_Numdof[i_vdesc]++;
                        ei_ptr->Lvdesc_Lnn_Numdof[i_vdesc][ln]++;
                      } /* for (indof = 0; indof < vd->Ndof; indof++) */
                      if (ei_ptr->dof[v] > MDE)
                        {
                          fprintf(stderr,"ERROR ei_ptr->dof[v]=%d is greater than MDE=%d\n",ei_ptr->dof[v],MDE);
                          EH(-1, "INCREASE MDE!");
                        }
                    } /* for (ivarDes = 0; ivarDes < nv->Num_Var_Desc_Per_Type[v] */
                }
              else
                {
                  /*
                   * We are at a point here where the variable type
                   * exists in the problem, but just not at this local node
                   */
                  ei_ptr->ln_to_dof[v][ln] = -1;
                  ei_ptr->ln_to_first_dof[v][ln] = -1;
                  iLv_start = ei_ptr->Lvdesc_First_Var_Type[v];
                  for (i_vd = 0; i_vd < ei_ptr->Num_Lvdesc_Per_Var_Type[v]; 
                       i_vd++) 
                    {
                      i_vdesc = iLv_start + i_vd;
                      ei_ptr->Lvdesc_Lnn_Numdof[i_vdesc][ln] = 0;
                    }
                }
            }
	  
        }
    }
  if (ledof > MDE * MAX_PROB_VAR) 
    {
      fprintf(stderr,"ERROR ledof is greater than MDE*MAX_PROB_VAR\n");
      EH(-1, "INCREASE MDE or MAX_PROB_VAR!"); 
    }
  /*
    #ifdef DEBUG_HKM
    if (ei_ptr->ielem == 6) {
    print_ei_index_desc();
    }
    #endif
  */
  /* 
   * Here we'll count up the dofs for the external variables, 
   * if any exist
   */
  if (efv->ev) {
    for (v = 0; v < efv->Num_external_field; v++) {
      ei_ptr->dof_ext[v] = 0;
      for (ln = 0; ln < nnodes; ln++) {

        /*
         * For this variable at this local node, how many dofs are needed?
         * (according to this element) Note this can be zero...
         *
         * WARNING..you will need to account for multiple species k here...
         */
        nunks = 1;  /* assumes only one external field variable of
                     * this type can be defined 
                     */

        EH(nunks, "problem with nun  for this external  var.");
        ei_ptr->dof_ext[v] += nunks;
      }
    }
  }

  if (num_elem_friends) 
    {
      ImAChild = is_shell_element(elem, exo);
      if (num_elem_friends[elem] == 0)
        {
          ImAChild = 0; 
        }
      ImAParent = 0;
      if (! ImAChild) 
        {
          if (num_elem_friends[elem] > 0)
            {
              ImAParent = 1;
            }
        }
      /*
       * Special section for children elements. i.e., those elements
       * such as shells that have parents.
       */
      if (ImAChild)
        {
          // for (i = 0; i < num_elem_friends[elem]; i++)
          for (i = 0; i < MIN(num_elem_friends[elem], 1); i++)
            {
              /* Get the element number of the parent */
              elem_Parent = elem_friends[elem][i];
              /* Get the element block number of the parent */
              n_eb_Parent = find_elemblock_index(elem_Parent, exo);
              /* Get the material number of the parent element */
              n_mn_Parent = Matilda[n_eb_Parent];
              PROBLEM_DESCRIPTION_STRUCT * pd_Parent = pd_glob[n_mn_Parent];
	 
              /* Get the element type of the parent */
              n_ielem_type_Parent = Elem_Type(exo, elem_Parent);	      
              /*
               *   Store the pointer to the beginning of this element's
               *   connectivity list
               */
              iconnect_ptr_Parent    = Proc_Connect_Ptr[elem_Parent];
              num_local_nodes_Parent = elem_info(NNODES, n_ielem_type_Parent);
              // ei_ptr->num_sides	      = shape2sides(ei_ptr->ielem_shape);
      
              for (v = V_FIRST; v < V_LAST; v++) 
                {
                  if (Num_Var_In_Type[imtrx][v])
                    {
                      lvdof = 0;
                      for (ln = 0; ln < num_local_nodes_Parent; ln++) 
                        {
                          /*
                           * For this local node "ln" what is the global node number, "gnn"?
                           */
                          gnn_Parent = Proc_Elem_Connect[iconnect_ptr_Parent + ln];
 
                          ni = Nodes[gnn_Parent];
                          nv = ni->Nodal_Vars_Info[imtrx];

                          /*
                           * Get the number of unknowns of the v at the node corresponding to
                           * the parent
                           */
                          nunks = get_nv_ndofs_modMF(nv, v);
                          enunks = nunks; 
                          if (nunks) 
                            {
                              /*
                               *  Modify the calculation for interpolations within the
                               *  element
                               */
                              enunks = dof_lnode_var_type(ln, etype, gnn_Parent, v, pd_Parent, imtrx);
                            }
                          /*
                           *  If the child element doesn't have a variable of that type
                           *  then we assign the column variable information to the 
                           *  parent element.
                           *
                           *  The pd->e[v] stuff is for the mesh equations.
                           */
                          if (ei_ptr->owningElementForColVar[v] == -1 || 
                              ((pd_glob[mn])->e[imtrx][v] == 0))
                            {
                              if (enunks) 
                                {
                                  ei_ptr->owningElementForColVar[v] = elem_Parent;
                                  if (mode == 0) 
                                    {
                                      linkParent_ei(elem, elem_Parent, v, exo, imtrx); 
                                      pd	 = pd_glob[mn];    
                                    }
                                } 
                            } 
                          else if (ei_ptr->owningElementForColVar[v] == elem)
                            {
                              if (enunks) 
                                {
                                  if ((v < MESH_DISPLACEMENT1) || (v > MESH_DISPLACEMENT3))
                                    {
                                      printf("Inconsistency unless mesh\n");
                                      exit(-1);
                                    }
                                } 
                            }
			  
                        }
                    }
                }
            }

        }

      /*
       * Special section for children elements. i.e., those elements
       * such as shells that have parents.
       */
      if (ImAParent)
        {
          for (i = 0; i < num_elem_friends[elem]; i++)
            {
              /* Everything below here that says Parent really means Child */
              /* Get the elment number of the parent */
              elem_Parent = elem_friends[elem][i];
              /* Get the element block number of the parent */
              n_eb_Parent = find_elemblock_index(elem_Parent, exo);
              /* Get the material number of the parent element */
              n_mn_Parent = Matilda[n_eb_Parent];

              PROBLEM_DESCRIPTION_STRUCT * pd_Parent = pd_glob[n_mn_Parent];
	 
              /* Get the element type of the parent */
              n_ielem_type_Parent = Elem_Type(exo, elem_Parent);
	      
              /*
               *   Store the pointer to the beginning of this element's
               *   connectivity list
               */
              iconnect_ptr_Parent    = Proc_Connect_Ptr[elem_Parent];
              num_local_nodes_Parent = elem_info(NNODES, n_ielem_type_Parent);
              // ei_ptr->num_sides	      = shape2sides(ei_ptr->ielem_shape);
      
              for (v = V_FIRST; v < V_LAST; v++) 
                {
                  if (Num_Var_In_Type[imtrx][v])
                    {
                      lvdof = 0;
                      for (ln = 0; ln < num_local_nodes_Parent; ln++) 
                        {
                          /*
                           * For this local node "ln" what is the global node number, "gnn"?
                           */
                          gnn = Proc_Elem_Connect[iconnect_ptr_Parent + ln];
 
                          ni = Nodes[gnn];
                          nv = ni->Nodal_Vars_Info[imtrx];

                          nunks = get_nv_ndofs_modMF(nv, v);
                          enunks = nunks; 
                          if (nunks) 
                            {
                              enunks = dof_lnode_var_type(ln, etype, gnn, v, pd_Parent, imtrx);
                            }
                          if (ei_ptr->owningElementForColVar[v] == -1)
                            {
                              if (enunks) 
                                {
                                  ei_ptr->owningElementForColVar[v] = elem_Parent;
                                  if (mode == 0) 
                                    {
                                      linkParent_ei(elem, elem_Parent, v, exo, imtrx); 
                                      pd	 = pd_glob[mn];    
                                    }
                                } 
                            } 
                          else if (ei_ptr->owningElementForColVar[v] == elem)
                            {
                              if (enunks) 
                                {
                                  if ((v < MESH_DISPLACEMENT1) || (v > MESH_DISPLACEMENT3)) 
                                    {
                                      printf("Inconsistency unless mesh\n");
                                      exit(-1);
                                    }
                                } 
                            }
			  
                        }
                    }
                }
            }

        }

      if (is_shell_element(elem, exo))
        {
          for (i = 0; i < num_elem_friends[elem]; i++)
            {
              n_eb = find_elemblock_index(elem_friends[elem][i], exo);
              n_mn = Matilda[n_eb];
              if (pd_glob[n_mn]->e[imtrx][R_MESH1]) 
                {
                  ei_ptr->deforming_mesh = TRUE;
                }
            }
        }
    }
  return (mn);
}
/**************************************************************************/

static void
load_varType_Interpolation_ptrs(const int varType, double **esp_ptr,
				double **esp_old_ptr, double **esp_dot_ptr)

     /***********************************************************************
      *
      * load_varType_Interpolation_ptrs:
      *
      *  Utility routine to fill up pointers to solution variables.
      *  (Might think about unrolling the loop in this routine in the
      *   future).
      ************************************************************************/
{ 
  int i, ie, dofs, ledof;
  int *lvdof_to_ledof_tmp = ei[pg->imtrx]->lvdof_to_ledof[varType];
  dofs = ei[pg->imtrx]->dof[varType];
  for (i = 0; i < dofs; i++) {
    ledof = lvdof_to_ledof_tmp[i];
    ie = ei[pg->imtrx]->ieqn_ledof[ledof];
    esp_ptr[i]     = x_static     + ie;
    esp_old_ptr[i] = x_old_static + ie;
    esp_dot_ptr[i] = xdot_static  + ie;
  }
}

static void
load_varType_Interpolation_ptrs_mat(int imtrx, const int varType, double **esp_ptr,
                                    double **esp_old_ptr, double **esp_dot_ptr)

     /***********************************************************************
      *
      * load_varType_Interpolation_ptrs:
      *
      *  Utility routine to fill up pointers to solution variables.
      *  (Might think about unrolling the loop in this routine in the
      *   future).
      ************************************************************************/
{ 
  int i, ie, dofs, ledof;
  int *lvdof_to_ledof_tmp = ei[imtrx]->lvdof_to_ledof[varType];
  dofs = ei[imtrx]->dof[varType];
  for (i = 0; i < dofs; i++) {
    ledof = lvdof_to_ledof_tmp[i];
    ie = ei[imtrx]->ieqn_ledof[ledof];
    esp_ptr[i]     = x_static     + ie;
    esp_old_ptr[i] = x_old_static + ie;
    esp_dot_ptr[i] = xdot_static  + ie;
  }
}
/***************************************************************************/
/***************************************************************************/
/***************************************************************************/

/* load_elem_dofptr -- load info for ea var how many dof's in elem; where to put
 *
 * requirements and input:
 * -----------------------
 *
 *	ei -- pointer to Element_Indeces structure
 *
 *	Comments:  "ei" should be already be partially filled with the
 * 	correct information about this element. Here, we simply load
 * 	up how many degrees of freedom are associated with each type of
 * 	variable (and kind) in this element. Suppose that TEMPERATURE has
 *	four degrees of freedom in this element. Then,
 *
 *		ei[pg->imtrx]->dof[TEMPERATURE] = 4
 *		ei[pg->imtrx]->dof[PRESSURE]    = 3    ( in 2D with piecewise linear)
 *
 *	and...
 *		ei[pg->imtrx]->dof_list[TEMPERATURE] = { 0, 1, 2, 3}
 *		ei[pg->imtrx]->dof_list[PRESSURE]    = { 8, 8, 8}
 *
 *	If none, then...
 *		ei[pg->imtrx]->dof[SOLUTE]	     = 0;
 *		ei[pg->imtrx]->dof_list[SOLUTE] = NULL;
 *
 *	*esp->T[i] = x[?] such that *esp->T[i] points to temperature at
 *		the local degree of freedom in the element.
 *
 * return values:
 * --------------
 *		0 -- everything went fine
 *	       -1 -- something went wrong
 */

int 
load_elem_dofptr(const int ielem,
		 const Exo_DB * exo,
                 dbl *x,
                 dbl *x_old,
                 dbl *xdot,
                 dbl *xdot_old,
		 const int early_return) /* if true we just get 
					  * some of the element information for 
					  * routines that need it like 
					  * global_h_elem_siz, but not all the 
					  * computationally expensive pointers */
{
  int eqn;			/* equation, variable name indeces */
  int gnn;			/* Global Node Number */
  int dofs=0;
  int nvdof, ledof, mode;
  int ie;			/* index into global solution */
  int i;			/* indeces for dofs */
  int b, c;			/* index for concentration */
  int p;			/* for vector eq.'s. */
  int status;
  int k;
  int R_s[MAX_MODES][DIM][DIM], R_g[DIM][DIM];
  int err;
  struct Level_Set_Data *ls_old;

#ifdef DEBUG
  int dim, eshape, etype, nnodes;
#endif /* DEBUG */

  x_static	       = x;
  x_old_static	       = x_old;
  xdot_static	       = xdot;
  xdot_old_static      = xdot_old;
  x_dbl_dot_static     = tran->xdbl_dot;
  x_dbl_dot_old_static = tran->xdbl_dot_old;

  /* load eqn and variable number in tensor form */
  (void) stress_eqn_pointer(R_s);

  R_g[0][0] = R_GRADIENT11;
  R_g[0][1] = R_GRADIENT12;
  R_g[1][0] = R_GRADIENT21;
  R_g[1][1] = R_GRADIENT22;
  R_g[0][2] = R_GRADIENT13;
  R_g[1][2] = R_GRADIENT23;
  R_g[2][0] = R_GRADIENT31;
  R_g[2][1] = R_GRADIENT32; 
  R_g[2][2] = R_GRADIENT33; 

  /*
   * Count how many dofs/node for each kind of variable...
   */

  status = 0;

  /*
   * Looking for the stuff that loaded up gun_list, ln_to_dof,
   * and all their friends ?
   * They've moved to a more suitable residence in the function load_ei.
   */
  int imtrx;
  for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
    load_ei(ielem, exo, 0, imtrx);
  }

#ifdef DEBUG
  fprintf(stderr, "P_%d: Loading elem dof ptrs in (d%d, s%d, t%d, n%d)\n",
	  ProcID, dim, eshape, etype, nnodes);
  for (i = 0; i < nnodes; i++) {
    gnn = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + i];
    fprintf(stderr, "P_%d: Node name %d is %d\n", ProcID, i, gnn);
  }

#endif  

  /*
   * Load up pointers into ...
   * 	global unknown vector,	"x"
   *	global residual vector,	"resid_vector"
   * 	global Jacobian matrix,	"a"
   *    external field vectors (if required)
   */

  /*
   * Look at the problem description to figure out which pointers make sense
   * to set up...then load them up with the right pointers into the global
   * residual and into the global A matrix...
   *
   * The equations and varibles at each node have a close correspondence so
   * we do them at the same time...
   */
  eqn = R_MESH1;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->d[0], esp_old->d[0],
				    esp_dot->d[0]);
    if (tran->solid_inertia) {
      dofs = ei[pg->imtrx]->dof[eqn];
      for ( i=0; i<dofs; i++) {
	gnn   = ei[pg->imtrx]->gnn_list[eqn][i];
	ie = Index_Solution(gnn, eqn, 0, 0, -1, pg->imtrx);
	esp_dbl_dot->d[0][i]   = tran->xdbl_dot     + ie;
      }
    }
  }

  eqn = R_MESH2;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->d[1], esp_old->d[1],
				    esp_dot->d[1]);
    if (tran->solid_inertia) {
      dofs = ei[pg->imtrx]->dof[eqn];
      for ( i=0; i<dofs; i++) {
	gnn   = ei[pg->imtrx]->gnn_list[eqn][i];
	ie = Index_Solution(gnn, eqn, 0, 0, -1, pg->imtrx);
	esp_dbl_dot->d[1][i]   = tran->xdbl_dot     + ie;
      }
    }
  }

  eqn   = R_MESH3;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->d[2], esp_old->d[2],
				    esp_dot->d[2]);
    if (tran->solid_inertia) {
      dofs = ei[pg->imtrx]->dof[eqn];
      for ( i=0; i<dofs; i++) {
	gnn   = ei[pg->imtrx]->gnn_list[eqn][i];
	ie = Index_Solution(gnn, eqn, 0, 0, -1, pg->imtrx);
	esp_dbl_dot->d[2][i]   = tran->xdbl_dot     + ie;
      }
    }
  }
  else if ((pd_glob[0]->CoordinateSystem == CYLINDRICAL ||
            pd_glob[0]->CoordinateSystem == SWIRLING ||
            pd_glob[0]->CoordinateSystem == CARTESIAN_2pt5D ||
	    pd_glob[0]->CoordinateSystem == PROJECTED_CARTESIAN) && upd->ep[pg->imtrx][R_MESH1] >= 0) {
    dofs = ei[pg->imtrx]->dof[R_MESH1];
    for (i = 0; i < dofs; i++) {
      esp->d[2][i]       = p0;
      esp_old->d[2][i]   = p0;
      esp_dot->d[2][i]   = p0;
      esp_dbl_dot->d[2][i] = p0;
    }
  }

  if (early_return) {
      return(status); /* return now with the element information and the pointers
			 into mesh displacement, if available */
  }

  eqn = R_SOLID1;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->d_rs[0], esp_old->d_rs[0],
				    esp_dot->d_rs[0]);
    if (tran->solid_inertia) {
      dofs = ei[pg->imtrx]->dof[eqn];
      for ( i=0; i<dofs; i++) {
	gnn   = ei[pg->imtrx]->gnn_list[eqn][i];
	ie = Index_Solution(gnn, eqn, 0, 0, -1, pg->imtrx);
	esp_dbl_dot->d_rs[0][i]   = tran->xdbl_dot     + ie;
      }
    }
  }

  eqn = R_SOLID2;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->d_rs[1], esp_old->d_rs[1],
				    esp_dot->d_rs[1]);
    if (tran->solid_inertia) {
      dofs = ei[pg->imtrx]->dof[eqn];
      for ( i=0; i<dofs; i++) {
	gnn   = ei[pg->imtrx]->gnn_list[eqn][i];
	ie = Index_Solution(gnn, eqn, 0, 0, -1, pg->imtrx);
	esp_dbl_dot->d_rs[1][i]   = tran->xdbl_dot     + ie;
      }
    }
  }

  eqn = R_SOLID3;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->d_rs[2], esp_old->d_rs[2],
				    esp_dot->d_rs[2]);
    if (tran->solid_inertia) {
      dofs = ei[pg->imtrx]->dof[eqn];
      for ( i=0; i<dofs; i++) {
	gnn   = ei[pg->imtrx]->gnn_list[eqn][i];
	ie = Index_Solution(gnn, eqn, 0, 0, -1, pg->imtrx);
	esp_dbl_dot->d_rs[2][i]   = tran->xdbl_dot     + ie;
      }
    }
  }
  else if((pd_glob[0]->CoordinateSystem == CYLINDRICAL ||
           pd_glob[0]->CoordinateSystem == SWIRLING ||
           pd_glob[0]->CoordinateSystem == CARTESIAN_2pt5D ||
	   pd_glob[0]->CoordinateSystem == PROJECTED_CARTESIAN) && upd->ep[pg->imtrx][R_SOLID1] >= 0) {
    dofs = ei[pg->imtrx]->dof[R_SOLID1];
    for ( i=0; i<dofs; i++) {
      esp->d_rs[2][i]      = p0;
      esp_old->d_rs[2][i]  = p0;
      esp_dot->d_rs[2][i]  = p0;
      esp_dbl_dot->d_rs[2][i] = p0;
    }
  }

  eqn = R_ENERGY;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->T, esp_old->T, esp_dot->T);
  }

  eqn = R_POTENTIAL;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->V, esp_old->V, esp_dot->V);
  }

  eqn = R_SURF_CHARGE;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->qs, esp_old->qs, esp_dot->qs);
  }

  eqn = R_FILL;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->F, esp_old->F, esp_dot->F);
  }

  eqn = R_CURVATURE;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->H, esp_old->H, esp_dot->H);
    }

  eqn = R_NORMAL1;
  if( upd->ep[pg->imtrx][eqn] >=0) {
    load_varType_Interpolation_ptrs(eqn, esp->n[0], esp_old->n[0], esp_dot->n[0]);
  }
     
  eqn = R_NORMAL2;
  if( upd->ep[pg->imtrx][eqn] >=0) {
    load_varType_Interpolation_ptrs(eqn, esp->n[1], esp_old->n[1], esp_dot->n[1]);
  }
     
  eqn = R_NORMAL3;
  if( upd->ep[pg->imtrx][eqn] >=0) {
    load_varType_Interpolation_ptrs(eqn, esp->n[2], esp_old->n[2], esp_dot->n[2]);
  }
  else if((pd_glob[0] ->CoordinateSystem == CYLINDRICAL ||
	   pd_glob[0]->CoordinateSystem == SWIRLING    ||
	   pd_glob[0]->CoordinateSystem == CARTESIAN_2pt5D    ||
	   pd_glob[0]->CoordinateSystem == PROJECTED_CARTESIAN ) &&
	  upd->ep[pg->imtrx][R_NORMAL1] >= 0) {
    dofs = ei[pg->imtrx]->dof[R_NORMAL1];
    for (i = 0; i < dofs; i++) 
      {
	esp->n[2][i]       = p0;
	esp_old->n[2][i]   = p0;
	esp_dot->n[2][i]   = p0;
      }
  }
     
  
  for(p = 0; p < DIM; p++)
    {
      eqn = R_VORT_DIR1 + p;
      if(upd->ep[pg->imtrx][eqn] >= 0)
	{
	  dofs = ei[pg->imtrx]->dof[eqn];
	  for(i = 0; i < dofs; i++)
	    {
	      ie = ei[pg->imtrx]->gun_list[eqn][i];
	      esp->vd[p][i] = x + ie;
	    }
	}
    }

  eqn = R_VORT_LAMBDA;
  if(upd->ep[pg->imtrx][eqn] >= 0)
    {
      dofs = ei[pg->imtrx]->dof[eqn];
      for(i = 0; i < dofs; i++)
	{
	  ie = ei[pg->imtrx]->gun_list[eqn][i];
	  esp->vlambda[i] = x + ie;
	}
    }

  eqn = R_BOND_EVOLUTION;
  if(pd->e[pg->imtrx][eqn])
    {
      load_varType_Interpolation_ptrs(eqn, esp->nn, esp_old->nn, esp_dot->nn);
    }

  eqn = R_SHEAR_RATE;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    dofs = ei[pg->imtrx]->dof[eqn];
    for ( i=0; i<dofs; i++) {
      ie = ei[pg->imtrx]->gun_list[eqn][i];
      esp->SH[i] = x + ie;
    }
  }

  eqn = R_ENORM;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    dofs = ei[pg->imtrx]->dof[eqn];
    for ( i=0; i<dofs; i++) {
      ie = ei[pg->imtrx]->gun_list[eqn][i];
      esp->Enorm[i] = x + ie;
      esp_old->Enorm[i] = x_old + ie;
    }
  }

  eqn = R_MOMENTUM1;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->v[0], esp_old->v[0],
				    esp_dot->v[0]);
  }
 
  eqn   = R_MOMENTUM2;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->v[1], esp_old->v[1],
				    esp_dot->v[1]);
  }

  if (*p0 != 0.) {
      EH(-1, "Hey, this zero is not zero!");
  }

  eqn = R_MOMENTUM3;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->v[2], esp_old->v[2],
				    esp_dot->v[2]);
  }
  else if((pd_glob[0] ->CoordinateSystem == CYLINDRICAL) &&
	  upd->ep[pg->imtrx][R_MOMENTUM1] >= 0) {
    dofs = ei[pg->imtrx]->dof[R_MOMENTUM1];
    for (i = 0; i < dofs; i++) {
      esp->v[2][i]       = p0;
      esp_old->v[2][i]   = p0;
      esp_dot->v[2][i]   = p0;
    }
  }


  eqn = R_EXT_VELOCITY;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->ext_v, esp_old->ext_v,
				    esp_dot->ext_v);
  }
 
  eqn = R_EFIELD1;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->E_field[0], esp_old->E_field[0],
				    esp_dot->E_field[0]);
  }
 
  eqn   = R_EFIELD2;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->E_field[1], esp_old->E_field[1],
				    esp_dot->E_field[1]);
  }

  eqn = R_EFIELD3;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->E_field[2], esp_old->E_field[2],
				    esp_dot->E_field[2]);
  }
  else if((pd_glob[0] ->CoordinateSystem == CYLINDRICAL) &&
	  upd->ep[pg->imtrx][R_EFIELD1] >= 0) {
    dofs = ei[pg->imtrx]->dof[R_EFIELD1];
    for (i = 0; i < dofs; i++) {
      esp->E_field[2][i]       = p0;
      esp_old->E_field[2][i]   = p0;
      esp_dot->E_field[2][i]   = p0;
    }
  }

  eqn = R_PMOMENTUM1;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->pv[0], esp_old->pv[0],
				    esp_dot->pv[0]);
  }

  eqn = R_PMOMENTUM2;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->pv[1], esp_old->pv[1],
				    esp_dot->pv[1]);
  }
  
  eqn   = R_PMOMENTUM3;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->pv[2], esp_old->pv[2],
				    esp_dot->pv[2]);
  }
  else if((pd_glob[0]->CoordinateSystem == CYLINDRICAL ||
           pd_glob[0]->CoordinateSystem == SWIRLING ||
           pd_glob[0]->CoordinateSystem == CARTESIAN_2pt5D ||
	   pd_glob[0]->CoordinateSystem == PROJECTED_CARTESIAN) &&
	  upd->ep[pg->imtrx][R_PMOMENTUM1] >= 0) {
    dofs = ei[pg->imtrx]->dof[R_PMOMENTUM1];
    for (i = 0; i < dofs; i++) {
      esp->pv[2][i]       = p0;
      esp_old->pv[2][i]   = p0;
      esp_dot->pv[2][i]   = p0;
    }
  }

  eqn = R_PRESSURE;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->P, esp_old->P, esp_dot->P);
  }


  eqn = R_LAGR_MULT1;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->lm[0], esp_old->lm[0],
				    esp_dot->lm[0]);
  }
 
  eqn = R_LAGR_MULT2;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->lm[1], esp_old->lm[1],
				    esp_dot->lm[1]);
  }

  eqn = R_LAGR_MULT3;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->lm[2], esp_old->lm[2],
				    esp_dot->lm[2]);
  }

  eqn = R_SHELL_CURVATURE;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->sh_K, esp_old->sh_K, esp_dot->sh_K);
  }

  eqn = R_SHELL_CURVATURE2;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->sh_K2, esp_old->sh_K2, esp_dot->sh_K2);
  }

  eqn = R_SHELL_TENSION;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->sh_tens, esp_old->sh_tens, esp_dot->sh_tens);
  }

  eqn = R_SHELL_X;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->sh_x, esp_old->sh_x, esp_dot->sh_x);
  }

  eqn = R_SHELL_Y;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->sh_y, esp_old->sh_y, esp_dot->sh_y);
  }
  
  eqn = R_SHELL_USER;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->sh_u, esp_old->sh_u, esp_dot->sh_u);
  }
  
  eqn = R_SHELL_ANGLE1;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->sh_ang[0], esp_old->sh_ang[0], esp_dot->sh_ang[0]);
    for (i = 0; i < ei[pg->imtrx]->dof[R_SHELL_ANGLE1]; i++) {
      if ( *esp->sh_ang[0][i] > M_PIE ) *esp->sh_ang[0][i] -= 2. * M_PIE;
      else if ( *esp->sh_ang[0][i] < -M_PIE ) *esp->sh_ang[0][i] += 2. * M_PIE;
    }
  }
  
  eqn = R_SHELL_ANGLE2;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->sh_ang[1], esp_old->sh_ang[1], esp_dot->sh_ang[1]);
    for (i = 0; i < ei[pg->imtrx]->dof[R_SHELL_ANGLE2]; i++) {
      if ( *esp->sh_ang[1][i] > M_PIE ) *esp->sh_ang[1][i] -= 2. * M_PIE;
      else if ( *esp->sh_ang[1][i] < -M_PIE ) *esp->sh_ang[1][i] += 2. * M_PIE;
    }
  }

 eqn = R_SHELL_SURF_DIV_V;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->div_s_v, esp_old->div_s_v, esp_dot->div_s_v);
    }

 eqn = R_SHELL_SURF_CURV;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->curv, esp_old->curv, esp_dot->curv);
  }

 eqn = R_N_DOT_CURL_V;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->n_dot_curl_s_v, esp_old->n_dot_curl_s_v, esp_dot->n_dot_curl_s_v);
  }

 eqn = R_GRAD_S_V_DOT_N1;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->grad_v_dot_n[0], esp_old->grad_v_dot_n[0], esp_dot->grad_v_dot_n[0]);
  }

 eqn = R_GRAD_S_V_DOT_N2;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->grad_v_dot_n[1], esp_old->grad_v_dot_n[1], esp_dot->grad_v_dot_n[1]);
  }

 eqn = R_GRAD_S_V_DOT_N3;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->grad_v_dot_n[2], esp_old->grad_v_dot_n[2], esp_dot->grad_v_dot_n[2]);
  }

 eqn = R_SHELL_DIFF_FLUX;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->sh_J, esp_old->sh_J, esp_dot->sh_J);
  }
 
 eqn = R_SHELL_DIFF_CURVATURE;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->sh_Kd, esp_old->sh_Kd, esp_dot->sh_Kd);
  }
 
 eqn = R_SHELL_NORMAL1;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->n[0], esp_old->n[0], esp_dot->n[0]);
  }
 
 eqn = R_SHELL_NORMAL2;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->n[1], esp_old->n[1], esp_dot->n[1]);
  }

 eqn = R_SHELL_NORMAL3;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->n[2], esp_old->n[2], esp_dot->n[2]);
  }

  for( b=0; b<MAX_PHASE_FUNC; b++) {
    eqn =  R_PHASE1+b;
    if( upd->ep[pg->imtrx][eqn]>=0 ) {
      load_varType_Interpolation_ptrs(eqn, esp->pF[b], esp_old->pF[b], esp_dot->pF[b]);
    }
  }

  eqn = R_ACOUS_PREAL;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->apr, esp_old->apr, esp_dot->apr);
  }
  eqn = R_ACOUS_PIMAG;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->api, esp_old->api, esp_dot->api);
  }
  eqn = R_ACOUS_REYN_STRESS;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->ars, esp_old->ars, esp_dot->ars);
  }
  eqn = R_SHELL_BDYVELO;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->sh_bv, esp_old->sh_bv, esp_dot->sh_bv);
  }
  eqn = R_SHELL_LUBP;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->sh_p, esp_old->sh_p, esp_dot->sh_p);
  }
  eqn = R_LUBP;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->lubp, esp_old->lubp, esp_dot->lubp);
  }
  eqn = R_LUBP_2;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->lubp_2, esp_old->lubp_2, esp_dot->lubp_2);
  }
  eqn = R_SHELL_FILMP;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->sh_fp, esp_old->sh_fp, esp_dot->sh_fp);
  }
  eqn = R_SHELL_FILMH;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->sh_fh, esp_old->sh_fh, esp_dot->sh_fh);
  }
  eqn = R_SHELL_PARTC;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->sh_pc, esp_old->sh_pc, esp_dot->sh_pc);
  }
  eqn = R_SHELL_SAT_CLOSED;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->sh_sat_closed, esp_old->sh_sat_closed, esp_dot->sh_sat_closed);
  }
  eqn = R_SHELL_SAT_OPEN;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->sh_p_open, esp_old->sh_p_open, esp_dot->sh_p_open);
  }
  eqn = R_SHELL_SAT_OPEN_2;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->sh_p_open_2, esp_old->sh_p_open_2, esp_dot->sh_p_open_2);
  }
  eqn = R_SHELL_ENERGY;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->sh_t, esp_old->sh_t, esp_dot->sh_t);
  }
  eqn = R_SHELL_DELTAH;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->sh_dh, esp_old->sh_dh, esp_dot->sh_dh);
  }
  eqn = R_SHELL_LUB_CURV;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->sh_l_curv, esp_old->sh_l_curv, esp_dot->sh_l_curv);
  }
  eqn = R_SHELL_LUB_CURV_2;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->sh_l_curv_2, esp_old->sh_l_curv_2, esp_dot->sh_l_curv_2);
  }
  eqn = R_SHELL_SAT_GASN;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->sh_sat_gasn, esp_old->sh_sat_gasn, esp_dot->sh_sat_gasn);
  }
  eqn = R_POR_SINK_MASS;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->sink_mass, esp_old->sink_mass, esp_dot->sink_mass);
  }
  eqn = R_LIGHT_INTP;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->poynt[0], esp_old->poynt[0],
				    esp_dot->poynt[0]);
  }
  eqn   = R_LIGHT_INTM;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->poynt[1], esp_old->poynt[1],
				    esp_dot->poynt[1]);
  }
  eqn = R_LIGHT_INTD;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->poynt[2], esp_old->poynt[2],
				    esp_dot->poynt[2]);
  }
  eqn = R_EM_E1_REAL;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->em_er[0], esp_old->em_er[0],
				    esp_dot->em_er[0]);
  }
  eqn = R_EM_E1_IMAG;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->em_ei[0], esp_old->em_ei[0],
				    esp_dot->em_ei[0]);
  }
  eqn = R_EM_E2_REAL;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->em_er[1], esp_old->em_er[1],
				    esp_dot->em_er[1]);
  }
  eqn = R_EM_E2_IMAG;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->em_ei[1], esp_old->em_ei[1],
				    esp_dot->em_ei[1]);
  }
  eqn = R_EM_E3_REAL;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->em_er[2], esp_old->em_er[2],
				    esp_dot->em_er[2]);
  }
  eqn = R_EM_E3_IMAG;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->em_ei[2], esp_old->em_ei[2],
				    esp_dot->em_ei[2]);
  }
  eqn = R_EM_H1_REAL;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->em_hr[0], esp_old->em_hr[0],
				    esp_dot->em_hr[0]);
  }
  eqn = R_EM_H1_IMAG;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->em_hi[0], esp_old->em_hi[0],
				    esp_dot->em_hi[0]);
  }
  eqn = R_EM_H2_REAL;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->em_hr[1], esp_old->em_hr[1],
				    esp_dot->em_hr[1]);
  }
  eqn = R_EM_H2_IMAG;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->em_hi[1], esp_old->em_hi[1],
				    esp_dot->em_hi[1]);
  }
  eqn = R_EM_H3_REAL;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->em_hr[2], esp_old->em_hr[2],
				    esp_dot->em_hr[2]);
  }
  eqn = R_EM_H3_IMAG;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->em_hi[2], esp_old->em_hi[2],
				    esp_dot->em_hi[2]);
  }

  eqn = R_RESTIME;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->restime, esp_old->restime,
				    esp_dot->restime);
  }  

  eqn = R_SHELL_SHEAR_TOP;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->sh_shear_top, esp_old->sh_shear_top, esp_dot->sh_shear_top);
  }

  eqn = R_SHELL_SHEAR_BOT;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->sh_shear_bot, esp_old->sh_shear_bot, esp_dot->sh_shear_bot);
  }

  eqn = R_SHELL_CROSS_SHEAR;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->sh_cross_shear, esp_old->sh_cross_shear, esp_dot->sh_cross_shear);
  }

  eqn = R_MAX_STRAIN;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->max_strain, esp_old->max_strain, esp_dot->max_strain);
  }

  eqn = R_CUR_STRAIN;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    load_varType_Interpolation_ptrs(eqn, esp->cur_strain, esp_old->cur_strain, esp_dot->cur_strain);
  }

 
  eqn = R_STRESS11;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    /* This should loop through all the stress variables 
     * for all the modes.
     */
    for (mode = 0; mode < vn->modes; mode++) {
      for (b = 0; b < VIM; b++) {
	for (c = 0; c < VIM; c++) {
	  if (b <= c) {
	    eqn = R_s[mode][b][c];
	    if (upd->ep[pg->imtrx][eqn] >= 0) {
              load_varType_Interpolation_ptrs(eqn, esp->S[mode][b][c],
					      esp_old->S[mode][b][c],
					      esp_dot->S[mode][b][c]);
            } else {
	      dofs = ei[pg->imtrx]->dof[R_STRESS11];
	      for (i = 0; i < dofs; i++) {
		esp->S[mode][b][c][i]       = p0;
		esp_old->S[mode][b][c][i]   = p0;
		esp_dot->S[mode][b][c][i]   = p0;
	      }
	    }
	  }
	}
      }
    }
    if((pd_glob[0]->CoordinateSystem == CYLINDRICAL ||
	pd_glob[0]->CoordinateSystem == SWIRLING ||
	pd_glob[0]->CoordinateSystem == CARTESIAN_2pt5D ||
	pd_glob[0]->CoordinateSystem == PROJECTED_CARTESIAN) ) {
      if( upd->ep[pg->imtrx][R_STRESS33] == -1 )
	  EH(-1,"Hey,the STRESS33 is needed in CYLINDRICAL VE problems!");
    }
  } 
  
  eqn = R_GRADIENT11;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    /* This should loop through all the velocity gradient 
     * components of the tensor
     */
    for (b = 0; b < VIM; b++) {
      for (c = 0; c < VIM; c++) {
	eqn = R_g[b][c];
	if (upd->ep[pg->imtrx][eqn] >= 0) {
	  load_varType_Interpolation_ptrs(eqn, esp->G[b][c],
					  esp_old->G[b][c], esp_dot->G[b][c]);
        } else {
	  dofs = ei[pg->imtrx]->dof[R_GRADIENT11];
	  for (i = 0; i < dofs; i++) {
	    esp->G[b][c][i]       = p0;
	    esp_old->G[b][c][i]   = p0;
	    esp_dot->G[b][c][i]   = p0;
	  }
	}
      }
    }
  }

  eqn = R_POR_LIQ_PRES;
  if ( upd->ep[pg->imtrx][eqn] >= 0 )
    {
      load_varType_Interpolation_ptrs(eqn, esp->p_liq, esp_old->p_liq,
				      esp_dot->p_liq);
    }

  eqn = R_POR_GAS_PRES;
  if ( upd->ep[pg->imtrx][eqn] >= 0 )
    {
      load_varType_Interpolation_ptrs(eqn, esp->p_gas, esp_old->p_gas,
				      esp_dot->p_gas);
    }

  eqn = R_POR_POROSITY;
  if ( upd->ep[pg->imtrx][eqn] >= 0 )
    {
      load_varType_Interpolation_ptrs(eqn, esp->porosity, esp_old->porosity,
				      esp_dot->porosity);
    }

  eqn = R_POR_ENERGY;
  if ( upd->ep[pg->imtrx][eqn] >= 0 )
    {
      load_varType_Interpolation_ptrs(eqn, esp->T, esp_old->T,
				      esp_dot->T);
    }

  eqn = R_POR_SATURATION;
  if ( upd->ep[pg->imtrx][eqn] >= 0)
    {
      EH(-1,"Saturation-based formulation not implemented yet");
    }

  eqn = R_TFMP_MASS;
  if ( upd->ep[pg->imtrx][eqn] >= 0)
    {
      load_varType_Interpolation_ptrs(eqn, esp->tfmp_sat, esp_old->tfmp_sat,
				      esp_dot->tfmp_sat);
    }

  eqn = R_TFMP_BOUND;
  if ( upd->ep[pg->imtrx][eqn] >= 0)
    {
      load_varType_Interpolation_ptrs(eqn, esp->tfmp_pres, esp_old->tfmp_pres,
				      esp_dot->tfmp_pres);
    }
  
  eqn  = R_MASS;
  if (upd->ep[pg->imtrx][eqn] >= 0) {
    for (k = 0; k < pd->Num_Species_Eqn; k++) {
      dofs = ei[pg->imtrx]->dof[eqn];
      for ( i=0; i<dofs; i++) {
	gnn   = ei[pg->imtrx]->gnn_list[eqn][i];
	nvdof = ei[pg->imtrx]->Baby_Dolphin[eqn][i];
	ledof = ei[pg->imtrx]->lvdof_to_ledof[eqn][i]; 
	ie = Index_Solution(gnn, R_MASS, k, nvdof, ei[pg->imtrx]->matID_ledof[ledof], pg->imtrx);
	esp->c[k][i]       = x        + ie;
	esp_old->c[k][i]   = x_old    + ie;
	esp_dot->c[k][i]   = xdot     + ie;
      }
    }
  }

  for (k = 0; k < MAX_MOMENTS; k++) {
    eqn = R_MOMENT0 + k;
    if (upd->ep[pg->imtrx][eqn] >= 0) {
      load_varType_Interpolation_ptrs(eqn, esp->moment[k], esp_old->moment[k],
				      esp_dot->moment[k]);
    }
  }
  eqn = R_DENSITY_EQN;
  if ( upd->ep[pg->imtrx][eqn] >= 0)
    {
      load_varType_Interpolation_ptrs(eqn, esp->rho, esp_old->rho,
				      esp_dot->rho);
    }

  /*
   * External field variables
   * Warning - here the distinction between nodes and dof's gets
   * kindof hazy
   */
  
  if (efv->ev) {
    for (k = 0; k < efv->Num_external_field; k++) {
	if( efv->i[k] != I_TABLE)
	  {
      	   dofs = ei[pg->imtrx]->dof_ext[k]; 
      	   for (i = 0; i < dofs; i++) {
		ie = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + i];
		evp->external_field[k][i] = efv->ext_fld_ndl_val[k] + ie;
      		}
	  }
    }
  }

  /*
   * initial displacement variables
   * Warning - here the distinction between nodes and dof's gets
   * kindof hazy
   */
  if (efv->TALE) {
    for (k = 0; k < exo->num_dim; k++) {
      dofs = ei[pg->imtrx]->dof[MESH_DISPLACEMENT1 + k];
      for (i = 0; i < dofs; i++) {
	ie = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + i];
	evp->initial_displacements[k][i] =
	    efv->init_displacement_ndl_val[k] + ie;
      }
      dofs = ei[pg->imtrx]->dof[SOLID_DISPLACEMENT1+k];
      for (i = 0; i < dofs; i++) {
	ie = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + i];
	evp->initial_displacements[k + DIM][i] =
	    efv->init_displacement_ndl_val[k+ DIM] + ie;
      }
    }
  }

  err = load_elem_dofptr_all(ielem, exo);
  EH(err, "load_elem_dofptr_all");

  /* extra setup needed for level set problems */
  if ( ls != NULL ) determine_ls_elem_overlap_state();
  ls_old = ls;
  if (pfd != NULL)
    {
      if (pfd->ls[0]->Evolution == LS_EVOLVE_SLAVE || pd->e[pg->imtrx][R_EXT_VELOCITY])
	{
	  ls = pfd->ls[0];
	  determine_ls_elem_overlap_state();
	}
    }
  ls = ls_old;

  /* extra setup needed for XFEM */
  if ( xfem != NULL ) load_xfem_for_elem( x, exo );
  


  return (status);
}

/* load_elem_dofptr -- load info for ea var how many dof's in elem; where to put
 *
 * requirements and input:
 * -----------------------
 *
 *	ei -- pointer to Element_Indeces structure
 *
 *	Comments:  "ei" should be already be partially filled with the
 * 	correct information about this element. Here, we simply load
 * 	up how many degrees of freedom are associated with each type of
 * 	variable (and kind) in this element. Suppose that TEMPERATURE has
 *	four degrees of freedom in this element. Then,
 *
 *		ei[pg->imtrx]->dof[TEMPERATURE] = 4
 *		ei[pg->imtrx]->dof[PRESSURE]    = 3    ( in 2D with piecewise linear)
 *
 *	and...
 *		ei[pg->imtrx]->dof_list[TEMPERATURE] = { 0, 1, 2, 3}
 *		ei[pg->imtrx]->dof_list[PRESSURE]    = { 8, 8, 8}
 *
 *	If none, then...
 *		ei[pg->imtrx]->dof[SOLUTE]	     = 0;
 *		ei[pg->imtrx]->dof_list[SOLUTE] = NULL;
 *
 *	*esp->T[i] = x[?] such that *esp->T[i] points to temperature at
 *		the local degree of freedom in the element.
 *
 * return values:
 * --------------
 *		0 -- everything went fine
 *	       -1 -- something went wrong
 */
/* TODO: REFACTOR OUT */
int 
load_elem_dofptr_all(const int ielem,
                     const Exo_DB * exo)
{
  int eqn;			/* equation, variable name indeces */
  int gnn;			/* Global Node Number */
  int dofs=0;
  int nvdof, ledof, mode;
  int ie;			/* index into global solution */
  int i;			/* indeces for dofs */
  int b, c;			/* index for concentration */
  int p;			/* for vector eq.'s. */
  int status;
  int k;
  int R_s[MAX_MODES][DIM][DIM], R_g[DIM][DIM];

#ifdef DEBUG
  int dim, eshape, etype, nnodes;
#endif /* DEBUG */

  if (upd->Total_Num_Matrices == 1) return 0;

  /* load eqn and variable number in tensor form */
  (void) stress_eqn_pointer(R_s);

  R_g[0][0] = R_GRADIENT11;
  R_g[0][1] = R_GRADIENT12;
  R_g[1][0] = R_GRADIENT21;
  R_g[1][1] = R_GRADIENT22;
  R_g[0][2] = R_GRADIENT13;
  R_g[1][2] = R_GRADIENT23;
  R_g[2][0] = R_GRADIENT31;
  R_g[2][1] = R_GRADIENT32; 
  R_g[2][2] = R_GRADIENT33; 

  /*
   * Count how many dofs/node for each kind of variable...
   */

  status = 0;

  if (upd->Total_Num_Matrices == 1) return 0;

  /*
   * Looking for the stuff that loaded up gun_list, ln_to_dof,
   * and all their friends ?
   * They've moved to a more suitable residence in the function load_ei.
   */
  int imtrx;
  for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
    if (pg->imtrx == imtrx) continue;
    load_ei(ielem, exo, 0, imtrx);
    x_static	       = pg->matrices[imtrx].x;
    x_old_static	       = pg->matrices[imtrx].x_old;
    xdot_static	       = pg->matrices[imtrx].xdot;
    xdot_old_static      = pg->matrices[imtrx].xdot_old;
    x_dbl_dot_static     = NULL;
    x_dbl_dot_old_static = NULL;

    double * x = x_static;
    double * x_old = x_old_static;
    double * xdot = xdot_static;

    /*
     * Load up pointers into ...
     * 	global unknown vector,	"x"
     *	global residual vector,	"resid_vector"
     * 	global Jacobian matrix,	"a"
     *    external field vectors (if required)
     */

    /*
     * Look at the problem description to figure out which pointers make sense
     * to set up...then load them up with the right pointers into the global
     * residual and into the global A matrix...
     *
     * The equations and varibles at each node have a close correspondence so
     * we do them at the same time...
     */
    eqn = R_MESH1;
    if (upd->ep[imtrx][eqn] >= 0) {
      load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->d[0], esp_old->d[0],
                                      esp_dot->d[0]);
      if (tran->solid_inertia) {
        dofs = ei[imtrx]->dof[eqn];
        for ( i=0; i<dofs; i++) {
          gnn   = ei[imtrx]->gnn_list[eqn][i];
          ie = Index_Solution(gnn, eqn, 0, 0, -1, imtrx);
          esp_dbl_dot->d[0][i]   = tran->xdbl_dot     + ie;
        }
      }
    }

    eqn = R_MESH2;
    if (upd->ep[imtrx][eqn] >= 0) {
      load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->d[1], esp_old->d[1],
                                      esp_dot->d[1]);
      if (tran->solid_inertia) {
        dofs = ei[imtrx]->dof[eqn];
        for ( i=0; i<dofs; i++) {
          gnn   = ei[imtrx]->gnn_list[eqn][i];
          ie = Index_Solution(gnn, eqn, 0, 0, -1, imtrx);
          esp_dbl_dot->d[1][i]   = tran->xdbl_dot     + ie;
        }
      }
    }

    eqn   = R_MESH3;
    if (upd->ep[imtrx][eqn] >= 0) {
      load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->d[2], esp_old->d[2],
                                      esp_dot->d[2]);
      if (tran->solid_inertia) {
        dofs = ei[imtrx]->dof[eqn];
        for ( i=0; i<dofs; i++) {
          gnn   = ei[imtrx]->gnn_list[eqn][i];
          ie = Index_Solution(gnn, eqn, 0, 0, -1, imtrx);
          esp_dbl_dot->d[2][i]   = tran->xdbl_dot     + ie;
        }
      }
    }
    else if ((pd_glob[0]->CoordinateSystem == CYLINDRICAL ||
              pd_glob[0]->CoordinateSystem == SWIRLING ||
              pd_glob[0]->CoordinateSystem == PROJECTED_CARTESIAN) && upd->ep[imtrx][R_MESH1] >= 0) {
      dofs = ei[imtrx]->dof[R_MESH1];
      for (i = 0; i < dofs; i++) {
        esp->d[2][i]       = p0;
        esp_old->d[2][i]   = p0;
        esp_dot->d[2][i]   = p0;
        esp_dbl_dot->d[2][i] = p0;
      }
    }

    eqn = R_SOLID1;
    if (upd->ep[imtrx][eqn] >= 0) {
      load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->d_rs[0], esp_old->d_rs[0],
                                      esp_dot->d_rs[0]);
      if (tran->solid_inertia) {
        dofs = ei[imtrx]->dof[eqn];
        for ( i=0; i<dofs; i++) {
          gnn   = ei[imtrx]->gnn_list[eqn][i];
          ie = Index_Solution(gnn, eqn, 0, 0, -1, imtrx);
          esp_dbl_dot->d_rs[0][i]   = tran->xdbl_dot     + ie;
        }
      }
    }

    eqn = R_SOLID2;
    if (upd->ep[imtrx][eqn] >= 0) {
      load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->d_rs[1], esp_old->d_rs[1],
                                      esp_dot->d_rs[1]);
      if (tran->solid_inertia) {
        dofs = ei[imtrx]->dof[eqn];
        for ( i=0; i<dofs; i++) {
          gnn   = ei[imtrx]->gnn_list[eqn][i];
          ie = Index_Solution(gnn, eqn, 0, 0, -1, imtrx);
          esp_dbl_dot->d_rs[1][i]   = tran->xdbl_dot     + ie;
        }
      }
    }

    eqn = R_SOLID3;
    if (upd->ep[imtrx][eqn] >= 0) {
      load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->d_rs[2], esp_old->d_rs[2],
                                      esp_dot->d_rs[2]);
      if (tran->solid_inertia) {
        dofs = ei[imtrx]->dof[eqn];
        for ( i=0; i<dofs; i++) {
          gnn   = ei[imtrx]->gnn_list[eqn][i];
          ie = Index_Solution(gnn, eqn, 0, 0, -1, imtrx);
          esp_dbl_dot->d_rs[2][i]   = tran->xdbl_dot     + ie;
        }
      }
    }
    else if((pd_glob[0]->CoordinateSystem == CYLINDRICAL ||
             pd_glob[0]->CoordinateSystem == SWIRLING ||
             pd_glob[0]->CoordinateSystem == PROJECTED_CARTESIAN) && upd->ep[imtrx][R_SOLID1] >= 0) {
      dofs = ei[imtrx]->dof[R_SOLID1];
      for ( i=0; i<dofs; i++) {
        esp->d_rs[2][i]      = p0;
        esp_old->d_rs[2][i]  = p0;
        esp_dot->d_rs[2][i]  = p0;
        esp_dbl_dot->d_rs[2][i] = p0;
      }
    }

    eqn = R_ENERGY;
    if (upd->ep[imtrx][eqn] >= 0) {
      load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->T, esp_old->T, esp_dot->T);
    }

    eqn = R_POTENTIAL;
    if (upd->ep[imtrx][eqn] >= 0) {
      load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->V, esp_old->V, esp_dot->V);
    }

    eqn = R_SURF_CHARGE;
    if (upd->ep[imtrx][eqn] >= 0) {
      load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->qs, esp_old->qs, esp_dot->qs);
    }

    eqn = R_FILL;
    if (upd->ep[imtrx][eqn] >= 0) {
      load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->F, esp_old->F, esp_dot->F);
    }

    eqn = R_CURVATURE;
    if (upd->ep[imtrx][eqn] >= 0) {
      load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->H, esp_old->H, esp_dot->H);
    }

    eqn = R_NORMAL1;
    if( upd->ep[imtrx][eqn] >=0) {
      load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->n[0], esp_old->n[0], esp_dot->n[0]);
    }
     
    eqn = R_NORMAL2;
    if( upd->ep[imtrx][eqn] >=0) {
      load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->n[1], esp_old->n[1], esp_dot->n[1]);
    }
     
    eqn = R_NORMAL3;
    if( upd->ep[imtrx][eqn] >=0) {
      load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->n[2], esp_old->n[2], esp_dot->n[2]);
    }
    else if((pd_glob[0] ->CoordinateSystem == CYLINDRICAL ||
             pd_glob[0]->CoordinateSystem == SWIRLING    ||
             pd_glob[0]->CoordinateSystem == PROJECTED_CARTESIAN ) &&
            upd->ep[imtrx][R_NORMAL1] >= 0) {
      dofs = ei[imtrx]->dof[R_NORMAL1];
      for (i = 0; i < dofs; i++) 
        {
          esp->n[2][i]       = p0;
          esp_old->n[2][i]   = p0;
          esp_dot->n[2][i]   = p0;
        }
    }
     
  
    for(p = 0; p < DIM; p++)
      {
        eqn = R_VORT_DIR1 + p;
        if(upd->ep[imtrx][eqn] >= 0)
          {
            dofs = ei[imtrx]->dof[eqn];
            for(i = 0; i < dofs; i++)
              {
                ie = ei[imtrx]->gun_list[eqn][i];
                esp->vd[p][i] = x + ie;
              }
          }
      }

    eqn = R_VORT_LAMBDA;
    if(upd->ep[imtrx][eqn] >= 0)
      {
        dofs = ei[imtrx]->dof[eqn];
        for(i = 0; i < dofs; i++)
          {
            ie = ei[imtrx]->gun_list[eqn][i];
            esp->vlambda[i] = x + ie;
          }
      }

    eqn = R_BOND_EVOLUTION;
    if(pd->e[imtrx][eqn])
      {
        load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->nn, esp_old->nn, esp_dot->nn);
      }

    eqn = R_SHEAR_RATE;
    if (upd->ep[imtrx][eqn] >= 0) {
      dofs = ei[imtrx]->dof[eqn];
      for ( i=0; i<dofs; i++) {
        ie = ei[imtrx]->gun_list[eqn][i];
        esp->SH[i] = x + ie;
      }
    }

    eqn = R_ENORM;
    if (upd->ep[imtrx][eqn] >= 0) {
      dofs = ei[imtrx]->dof[eqn];
      for ( i=0; i<dofs; i++) {
        ie = ei[imtrx]->gun_list[eqn][i];
        esp->Enorm[i] = x + ie;
        esp_old->Enorm[i] = x_old + ie;
      }
    }

    eqn = R_MOMENTUM1;
    if (upd->ep[imtrx][eqn] >= 0) {
      load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->v[0], esp_old->v[0],
                                      esp_dot->v[0]);
    }
 
    eqn   = R_MOMENTUM2;
    if (upd->ep[imtrx][eqn] >= 0) {
      load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->v[1], esp_old->v[1],
                                      esp_dot->v[1]);
    }

    if (*p0 != 0.) {
      EH(-1, "Hey, this zero is not zero!");
    }

    eqn = R_MOMENTUM3;
    if (upd->ep[imtrx][eqn] >= 0) {
      load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->v[2], esp_old->v[2],
                                      esp_dot->v[2]);
    }
    else if((pd_glob[0] ->CoordinateSystem == CYLINDRICAL) &&
            upd->ep[imtrx][R_MOMENTUM1] >= 0) {
      dofs = ei[imtrx]->dof[R_MOMENTUM1];
      for (i = 0; i < dofs; i++) {
        esp->v[2][i]       = p0;
        esp_old->v[2][i]   = p0;
        esp_dot->v[2][i]   = p0;
      }
    }


    eqn = R_EXT_VELOCITY;
    if (upd->ep[imtrx][eqn] >= 0) {
      load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->ext_v, esp_old->ext_v,
                                      esp_dot->ext_v);
    }
 
    eqn = R_EFIELD1;
    if (upd->ep[imtrx][eqn] >= 0) {
      load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->E_field[0], esp_old->E_field[0],
                                      esp_dot->E_field[0]);
    }
 
    eqn   = R_EFIELD2;
    if (upd->ep[imtrx][eqn] >= 0) {
      load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->E_field[1], esp_old->E_field[1],
                                      esp_dot->E_field[1]);
    }

    eqn = R_EFIELD3;
    if (upd->ep[imtrx][eqn] >= 0) {
      load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->E_field[2], esp_old->E_field[2],
                                      esp_dot->E_field[2]);
    }
    else if((pd_glob[0] ->CoordinateSystem == CYLINDRICAL) &&
            upd->ep[imtrx][R_EFIELD1] >= 0) {
      dofs = ei[imtrx]->dof[R_EFIELD1];
      for (i = 0; i < dofs; i++) {
        esp->E_field[2][i]       = p0;
        esp_old->E_field[2][i]   = p0;
        esp_dot->E_field[2][i]   = p0;
      }
    }

    eqn = R_PMOMENTUM1;
    if (upd->ep[imtrx][eqn] >= 0) {
      load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->pv[0], esp_old->pv[0],
                                      esp_dot->pv[0]);
    }

    eqn = R_PMOMENTUM2;
    if (upd->ep[imtrx][eqn] >= 0) {
      load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->pv[1], esp_old->pv[1],
                                      esp_dot->pv[1]);
    }
  
    eqn   = R_PMOMENTUM3;
    if (upd->ep[imtrx][eqn] >= 0) {
      load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->pv[2], esp_old->pv[2],
                                      esp_dot->pv[2]);
    }
    else if((pd_glob[0]->CoordinateSystem == CYLINDRICAL ||
             pd_glob[0]->CoordinateSystem == SWIRLING ||
             pd_glob[0]->CoordinateSystem == PROJECTED_CARTESIAN) &&
            upd->ep[imtrx][R_PMOMENTUM1] >= 0) {
      dofs = ei[imtrx]->dof[R_PMOMENTUM1];
      for (i = 0; i < dofs; i++) {
        esp->pv[2][i]       = p0;
        esp_old->pv[2][i]   = p0;
        esp_dot->pv[2][i]   = p0;
      }
    }

    eqn = R_PRESSURE;
    if (upd->ep[imtrx][eqn] >= 0) {
      load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->P, esp_old->P, esp_dot->P);
    }

    eqn = R_LAGR_MULT1;
    if (upd->ep[imtrx][eqn] >= 0) {
      load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->lm[0], esp_old->lm[0],
                                      esp_dot->lm[0]);
    }
 
    eqn = R_LAGR_MULT2;
    if (upd->ep[imtrx][eqn] >= 0) {
      load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->lm[1], esp_old->lm[1],
                                      esp_dot->lm[1]);
    }

    eqn = R_LAGR_MULT3;
    if (upd->ep[imtrx][eqn] >= 0) {
      load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->lm[2], esp_old->lm[2],
                                      esp_dot->lm[2]);
    }

    eqn = R_SHELL_CURVATURE;
    if (upd->ep[imtrx][eqn] >= 0) {
      load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->sh_K, esp_old->sh_K, esp_dot->sh_K);
    }

    eqn = R_SHELL_TENSION;
    if (upd->ep[imtrx][eqn] >= 0) {
      load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->sh_tens, esp_old->sh_tens, esp_dot->sh_tens);
    }

    eqn = R_SHELL_X;
    if (upd->ep[imtrx][eqn] >= 0) {
      load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->sh_x, esp_old->sh_x, esp_dot->sh_x);
    }

    eqn = R_SHELL_Y;
    if (upd->ep[imtrx][eqn] >= 0) {
      load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->sh_y, esp_old->sh_y, esp_dot->sh_y);
    }
  
    eqn = R_SHELL_USER;
    if (upd->ep[imtrx][eqn] >= 0) {
      load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->sh_u, esp_old->sh_u, esp_dot->sh_u);
    }
  
    eqn = R_SHELL_ANGLE1;
    if (upd->ep[imtrx][eqn] >= 0) {
      load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->sh_ang[0], esp_old->sh_ang[0], esp_dot->sh_ang[0]);
      for (i = 0; i < ei[imtrx]->dof[R_SHELL_ANGLE1]; i++) {
        if ( *esp->sh_ang[0][i] > M_PIE ) *esp->sh_ang[0][i] -= 2. * M_PIE;
        else if ( *esp->sh_ang[0][i] < -M_PIE ) *esp->sh_ang[0][i] += 2. * M_PIE;
      }
    }
  
    eqn = R_SHELL_ANGLE2;
    if (upd->ep[imtrx][eqn] >= 0) {
      load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->sh_ang[1], esp_old->sh_ang[1], esp_dot->sh_ang[1]);
      for (i = 0; i < ei[imtrx]->dof[R_SHELL_ANGLE2]; i++) {
        if ( *esp->sh_ang[1][i] > M_PIE ) *esp->sh_ang[1][i] -= 2. * M_PIE;
        else if ( *esp->sh_ang[1][i] < -M_PIE ) *esp->sh_ang[1][i] += 2. * M_PIE;
      }
    }

    eqn = R_SHELL_SURF_DIV_V;
    if (upd->ep[imtrx][eqn] >= 0) {
      load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->div_s_v, esp_old->div_s_v, esp_dot->div_s_v);
    }

    eqn = R_SHELL_SURF_CURV;
    if (upd->ep[imtrx][eqn] >= 0) {
      load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->curv, esp_old->curv, esp_dot->curv);
    }

    eqn = R_N_DOT_CURL_V;
    if (upd->ep[imtrx][eqn] >= 0) {
      load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->n_dot_curl_s_v, esp_old->n_dot_curl_s_v, esp_dot->n_dot_curl_s_v);
    }

    eqn = R_GRAD_S_V_DOT_N1;
    if (upd->ep[imtrx][eqn] >= 0) {
      load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->grad_v_dot_n[0], esp_old->grad_v_dot_n[0], esp_dot->grad_v_dot_n[0]);
    }

    eqn = R_GRAD_S_V_DOT_N2;
    if (upd->ep[imtrx][eqn] >= 0) {
      load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->grad_v_dot_n[1], esp_old->grad_v_dot_n[1], esp_dot->grad_v_dot_n[1]);
    }

    eqn = R_GRAD_S_V_DOT_N3;
    if (upd->ep[imtrx][eqn] >= 0) {
      load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->grad_v_dot_n[2], esp_old->grad_v_dot_n[2], esp_dot->grad_v_dot_n[2]);
    }

    eqn = R_SHELL_DIFF_FLUX;
    if (upd->ep[imtrx][eqn] >= 0) {
      load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->sh_J, esp_old->sh_J, esp_dot->sh_J);
    }
 
    eqn = R_SHELL_DIFF_CURVATURE;
    if (upd->ep[imtrx][eqn] >= 0) {
      load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->sh_Kd, esp_old->sh_Kd, esp_dot->sh_Kd);
    }
 
    eqn = R_SHELL_NORMAL1;
    if (upd->ep[imtrx][eqn] >= 0) {
      load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->n[0], esp_old->n[0], esp_dot->n[0]);
    }
 
    eqn = R_SHELL_NORMAL2;
    if (upd->ep[imtrx][eqn] >= 0) {
      load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->n[1], esp_old->n[1], esp_dot->n[1]);
    }

    for( b=0; b<MAX_PHASE_FUNC; b++) {
      eqn =  R_PHASE1+b;
      if( upd->ep[imtrx][eqn]>=0 ) {
        load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->pF[b], esp_old->pF[b], esp_dot->pF[b]);
      }
    }

    eqn = R_ACOUS_PREAL;
    if (upd->ep[imtrx][eqn] >= 0) {
      load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->apr, esp_old->apr, esp_dot->apr);
    }
    eqn = R_ACOUS_PIMAG;
    if (upd->ep[imtrx][eqn] >= 0) {
      load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->api, esp_old->api, esp_dot->api);
    }
    eqn = R_ACOUS_REYN_STRESS;
    if (upd->ep[imtrx][eqn] >= 0) {
      load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->ars, esp_old->ars, esp_dot->ars);
    }
    eqn = R_SHELL_BDYVELO;
    if (upd->ep[imtrx][eqn] >= 0) {
      load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->sh_bv, esp_old->sh_bv, esp_dot->sh_bv);
    }
    eqn = R_SHELL_LUBP;
    if (upd->ep[imtrx][eqn] >= 0) {
      load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->sh_p, esp_old->sh_p, esp_dot->sh_p);
    }
    eqn = R_LUBP;
    if (upd->ep[imtrx][eqn] >= 0) {
      load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->lubp, esp_old->lubp, esp_dot->lubp);
    }
    eqn = R_LUBP_2;
    if (upd->ep[imtrx][eqn] >= 0) {
      load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->lubp_2, esp_old->lubp_2, esp_dot->lubp_2);
    }
    eqn = R_SHELL_FILMP;
    if (upd->ep[imtrx][eqn] >= 0) {
      load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->sh_fp, esp_old->sh_fp, esp_dot->sh_fp);
    }
    eqn = R_SHELL_FILMH;
    if (upd->ep[imtrx][eqn] >= 0) {
      load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->sh_fh, esp_old->sh_fh, esp_dot->sh_fh);
    }
    eqn = R_SHELL_PARTC;
    if (upd->ep[imtrx][eqn] >= 0) {
      load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->sh_pc, esp_old->sh_pc, esp_dot->sh_pc);
    }
    eqn = R_SHELL_SAT_CLOSED;
    if (upd->ep[imtrx][eqn] >= 0) {
      load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->sh_sat_closed, esp_old->sh_sat_closed, esp_dot->sh_sat_closed);
    }
    eqn = R_SHELL_SAT_OPEN;
    if (upd->ep[imtrx][eqn] >= 0) {
      load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->sh_p_open, esp_old->sh_p_open, esp_dot->sh_p_open);
    }
    eqn = R_SHELL_SAT_OPEN_2;
    if (upd->ep[imtrx][eqn] >= 0) {
      load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->sh_p_open_2, esp_old->sh_p_open_2, esp_dot->sh_p_open_2);
    }
    eqn = R_SHELL_ENERGY;
    if (upd->ep[imtrx][eqn] >= 0) {
      load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->sh_t, esp_old->sh_t, esp_dot->sh_t);
    }
    eqn = R_SHELL_DELTAH;
    if (upd->ep[imtrx][eqn] >= 0) {
      load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->sh_dh, esp_old->sh_dh, esp_dot->sh_dh);
    }
    eqn = R_SHELL_LUB_CURV;
    if (upd->ep[imtrx][eqn] >= 0) {
      load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->sh_l_curv, esp_old->sh_l_curv, esp_dot->sh_l_curv);
    }
    eqn = R_SHELL_LUB_CURV_2;
    if (upd->ep[imtrx][eqn] >= 0) {
      load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->sh_l_curv_2, esp_old->sh_l_curv_2, esp_dot->sh_l_curv_2);
    }
    eqn = R_SHELL_SAT_GASN;
    if (upd->ep[imtrx][eqn] >= 0) {
      load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->sh_sat_gasn, esp_old->sh_sat_gasn, esp_dot->sh_sat_gasn);
    }
    eqn = R_POR_SINK_MASS;
    if (upd->ep[imtrx][eqn] >= 0) {
      load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->sink_mass, esp_old->sink_mass, esp_dot->sink_mass);
    }
    eqn = R_LIGHT_INTP;
    if (upd->ep[imtrx][eqn] >= 0) {
      load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->poynt[0], esp_old->poynt[0],
                                      esp_dot->poynt[0]);
    }
    eqn   = R_LIGHT_INTM;
    if (upd->ep[imtrx][eqn] >= 0) {
      load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->poynt[1], esp_old->poynt[1],
                                      esp_dot->poynt[1]);
    }
    eqn = R_LIGHT_INTD;
    if (upd->ep[imtrx][eqn] >= 0) {
      load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->poynt[2], esp_old->poynt[2],
                                      esp_dot->poynt[2]);
    }

    eqn = R_SHELL_SHEAR_TOP;
    if (upd->ep[imtrx][eqn] >= 0) {
      load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->sh_shear_top, esp_old->sh_shear_top, esp_dot->sh_shear_top);
    }

    eqn = R_SHELL_SHEAR_BOT;
    if (upd->ep[imtrx][eqn] >= 0) {
      load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->sh_shear_bot, esp_old->sh_shear_bot, esp_dot->sh_shear_bot);
    }

    eqn = R_SHELL_CROSS_SHEAR;
    if (upd->ep[imtrx][eqn] >= 0) {
      load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->sh_cross_shear, esp_old->sh_cross_shear, esp_dot->sh_cross_shear);
    }

    eqn = R_MAX_STRAIN;
    if (upd->ep[imtrx][eqn] >= 0) {
      load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->max_strain, esp_old->max_strain, esp_dot->max_strain);
    }

    eqn = R_CUR_STRAIN;
    if (upd->ep[imtrx][eqn] >= 0) {
      load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->cur_strain, esp_old->cur_strain, esp_dot->cur_strain);
    }

 
    eqn = R_STRESS11;
    if (upd->ep[imtrx][eqn] >= 0) {
      /* This should loop through all the stress variables 
       * for all the modes.
       */
      for (mode = 0; mode < vn->modes; mode++) {
        for (b = 0; b < VIM; b++) {
          for (c = 0; c < VIM; c++) {
            if (b <= c) {
              eqn = R_s[mode][b][c];
              if (upd->ep[imtrx][eqn] >= 0) {
                load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->S[mode][b][c],
                                                esp_old->S[mode][b][c],
                                                esp_dot->S[mode][b][c]);
              } else {
                dofs = ei[imtrx]->dof[R_STRESS11];
                for (i = 0; i < dofs; i++) {
                  esp->S[mode][b][c][i]       = p0;
                  esp_old->S[mode][b][c][i]   = p0;
                  esp_dot->S[mode][b][c][i]   = p0;
                }
              }
            }
          }
        }
      }
      if((pd_glob[0]->CoordinateSystem == CYLINDRICAL ||
          pd_glob[0]->CoordinateSystem == SWIRLING ||
          pd_glob[0]->CoordinateSystem == PROJECTED_CARTESIAN) ) {
        if( upd->ep[imtrx][R_STRESS33] == -1 )
	  EH(-1,"Hey,the STRESS33 is needed in CYLINDRICAL VE problems!");
      }
    } 
  
    eqn = R_GRADIENT11;
    if (upd->ep[imtrx][eqn] >= 0) {
      /* This should loop through all the velocity gradient 
       * components of the tensor
       */
      for (b = 0; b < VIM; b++) {
        for (c = 0; c < VIM; c++) {
          eqn = R_g[b][c];
          if (upd->ep[imtrx][eqn] >= 0) {
            load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->G[b][c],
                                            esp_old->G[b][c], esp_dot->G[b][c]);
          } else {
            dofs = ei[imtrx]->dof[R_GRADIENT11];
            for (i = 0; i < dofs; i++) {
              esp->G[b][c][i]       = p0;
              esp_old->G[b][c][i]   = p0;
              esp_dot->G[b][c][i]   = p0;
            }
          }
        }
      }
    }

    eqn = R_POR_LIQ_PRES;
    if ( upd->ep[imtrx][eqn] >= 0 )
      {
        load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->p_liq, esp_old->p_liq,
                                        esp_dot->p_liq);
      }

    eqn = R_POR_GAS_PRES;
    if ( upd->ep[imtrx][eqn] >= 0 )
      {
        load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->p_gas, esp_old->p_gas,
                                        esp_dot->p_gas);
      }

    eqn = R_POR_POROSITY;
    if ( upd->ep[imtrx][eqn] >= 0 )
      {
        load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->porosity, esp_old->porosity,
                                        esp_dot->porosity);
      }

    eqn = R_POR_ENERGY;
    if ( pd->e[imtrx][eqn] )
      {
        load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->T, esp_old->T,
                                        esp_dot->T);
      }

    eqn = R_POR_SATURATION;
    if ( upd->ep[imtrx][eqn] >= 0 )
      {
        EH(-1,"Saturation-based formulation not implemented yet");
      }
  
    eqn  = R_MASS;
    if (upd->ep[imtrx][eqn] >= 0) {
      for (k = 0; k < pd->Num_Species_Eqn; k++) {
        dofs = ei[imtrx]->dof[eqn];
        for ( i=0; i<dofs; i++) {
          gnn   = ei[imtrx]->gnn_list[eqn][i];
          nvdof = ei[imtrx]->Baby_Dolphin[eqn][i];
          ledof = ei[imtrx]->lvdof_to_ledof[eqn][i]; 
          ie = Index_Solution(gnn, R_MASS, k, nvdof, ei[imtrx]->matID_ledof[ledof], imtrx);
          esp->c[k][i]       = x        + ie;
          esp_old->c[k][i]   = x_old    + ie;
          esp_dot->c[k][i]   = xdot     + ie;
        }
      }
    }

    for (k = 0; k < MAX_MOMENTS; k++) {
      eqn = R_MOMENT0 + k;
      if (upd->ep[imtrx][eqn] >= 0) {
	load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->moment[k], esp_old->moment[k],
					esp_dot->moment[k]);
      }
    }
    eqn = R_DENSITY_EQN;
    if ( pd->e[imtrx][eqn] )
      {
	load_varType_Interpolation_ptrs_mat(imtrx, eqn, esp->rho, esp_old->rho,
					esp_dot->rho);
      }

  }

  x_static	       = pg->matrices[pg->imtrx].x;
  x_old_static	       = pg->matrices[pg->imtrx].x_old;
  xdot_static	       = pg->matrices[pg->imtrx].xdot;
  xdot_old_static      = pg->matrices[pg->imtrx].xdot_old;
  x_dbl_dot_static     = NULL;
  x_dbl_dot_old_static = NULL;

  return (status);
}
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
#ifdef DEBUG_NOT
static void
pvn(int var_type, VARIABLE_DESCRIPTION_STRUCT *vd)
{
  if (var_type == 4) {
    printf("4(%-d)    ", vd->Subvar_Index);   
  } else {
    printf("%-5d   ", var_type);
  }
}
#endif
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
#ifdef DEBUG_NOT
static void
print_ei_index_desc()

    /******************************************************************
     *
     * print_ei_index_desc()
     *
     * This routine will print out information about element
     * variable indexing mappings in the global structure, ei,
     * for the current element. 
     * Mostly a debugging tool, so I have ifdefed it out for normal
     * goma usage. However, it's very valuable as a debugging aid.
     *
     ******************************************************************/
{
  int i, j, ln, ledof, *numidofs;
  VARIABLE_DESCRIPTION_STRUCT *vd;
  int nnodes = ei[pg->imtrx]->num_local_nodes;
  numidofs = alloc_int_1(nnodes, 0);
  printf("\n");
  fprint_line(stdout,"+", 80);
  printf(" Element Number = %d, MatID = %d\n", ei[pg->imtrx]->ielem, ei[pg->imtrx]->mn); 
  printf(" Number of variable descriptions in the element: %d\n", 
	 ei[pg->imtrx]->Num_Lvdesc);

  /*
   *  SECTION 1
   * Printout of all member arrays using the local variable
   * description as a single index
   */ 
  printf("  Printout of Local Variables Descriptions:\n");
  printf("     ");
  printf("lvdesc  Var_Type  MatID  Numdof");
  printf("\n");
  printf("     "); fprint_line(stdout,"-", 60); printf("\n");
  for (i = 0; i < ei[pg->imtrx]->Num_Lvdesc; i++) {
    vd = ei[pg->imtrx]->Lvdesc_vd_ptr[i];
    printf("     ");
    printf("%-5d", i);
    pvn(ei[pg->imtrx]->Lvdesc_to_Var_Type[i], vd);
    printf("  %-5d", ei[pg->imtrx]->Lvdesc_to_MatID[i]);
    printf("  %-5d", ei[pg->imtrx]->Lvdesc_Numdof[i]);
    printf("\n");
  }
  printf("     "); fprint_line(stdout,"-", 60); printf("\n\n");

  /*
   *  SECTION 2
   * Printout of all member arrays using the local variable 
   * description - local vdesc dof pair as arguments. We also
   * include local variable description - local node number
   * arrays here too.
   */
  printf("     ");fprint_line(stdout,"=", 75);
  printf("     Printout of LVDESC DOFS and LNN INFORMATION:\n");
  printf("     ");fprint_line(stdout,"=", 75); printf("\n");
  for (i = 0; i < ei[pg->imtrx]->Num_Lvdesc; i++) {
    vd = ei[pg->imtrx]->Lvdesc_vd_ptr[i];
    printf("     ");fprint_line(stdout,"=", 65); printf("\n");
    printf("     "); 
    printf("LVDESC = %-d, Var = ", i);
    pvn(ei[pg->imtrx]->Lvdesc_to_Var_Type[i], vd);
    printf(", MatID = %-5d", ei[pg->imtrx]->Lvdesc_to_MatID[i]);
    printf("\n");
    printf("        ");
    printf("i_vdesc  lnn  lvdof   ledof   gnn   gun  active owned\n");
    printf("        "); fprint_line(stdout,"-", 60);
    for (j = 0; j < ei[pg->imtrx]->Lvdesc_Numdof[i]; j++) {
      printf("        ");   
      printf("  %-5d", j);
      printf("  %-5d ", ei[pg->imtrx]->Lvdesc_to_Lnn[i][j]);
      printf("  %-5d ", ei[pg->imtrx]->Lvdesc_to_lvdof[i][j]);
      ledof = ei[pg->imtrx]->Lvdesc_to_ledof[i][j];
      printf(" %-5d ", ledof);
      printf("%-6d", ei[pg->imtrx]->Lvdesc_to_Gnn[i][j]);
      printf("%-6d ", ei[pg->imtrx]->Lvdesc_to_Gun[i][j]);
      printf("%-6d ", ei[pg->imtrx]->active_interp_ledof[ledof]);
      printf("%-6d ", ei[pg->imtrx]->owned_ledof[ledof]);
      if (ei[pg->imtrx]->ieqn_ledof[ledof] != ei[pg->imtrx]->Lvdesc_to_Gun[i][j]) {
	printf("WARNING ieqn_edof = %5d", ei[pg->imtrx]->ieqn_ledof[ledof]);
      }
      printf("\n");
    }
    printf("        "); fprint_line(stdout,"-", 60);

    printf("        ");
    printf("lnn Numdof  lvdof  soln_offset\n");
    printf("        "); fprint_line(stdout,"-", 60);
    for (ln = 0; ln < nnodes; ln++) {
      printf("        ");   
      printf(" %-5d", ln);
      printf("  %-5d", ei[pg->imtrx]->Lvdesc_Lnn_Numdof[i][ln]);
      printf("  %-5d", ei[pg->imtrx]->Lvdesc_Lnn_to_lvdof[i][ln]);
      printf("    %-5d", ei[pg->imtrx]->Lvdesc_Lnn_to_Offset[i][ln]);
      printf("\n");
    }
    printf("        "); fprint_line(stdout,"-", 60);
  }

  /*
   *  SECTION 3
   * Printout of all member arrays using the local variable
   * type as a single index. Var types with no dofs are not
   * printed.
   */
  printf("     ");fprint_line(stdout,"=", 75);
  printf("     Printout of LVDOF INFORMATION:\n");
  printf("     ");fprint_line(stdout,"=", 75); printf("\n");
  printf("     ");
  printf("var_type   Num_dofs  Num_Lvdesc");
  printf("\n");
  printf("     "); fprint_line(stdout,"-", 60); printf("\n");
  for (i = 0; i <  MAX_VARIABLE_TYPES; i++) {
    if (ei[pg->imtrx]->dof[i] > 0) {
      printf("     ");
      printf("     %-5d", i);
      printf("    %-5d", ei[pg->imtrx]->dof[i]);
      printf("    %-5d", ei[pg->imtrx]->Num_Lvdesc_Per_Var_Type[i]);
      printf("\n");
    }
  }
  printf("     "); fprint_line(stdout,"-", 60); printf("\n\n");

  /*
   *  SECTION 4
   * Printout of all structures using the local variable
   * type dof arrays
   */
  printf("     "); fprint_line(stdout,"=", 75);
  printf("     Printout of LVDOF IDOF and LNN ARRAY INFORMATION:\n");
  printf("     ");fprint_line(stdout,"=", 75); printf("\n");
  for (i = 0; i < MAX_VARIABLE_TYPES; i++) {
    if (ei[pg->imtrx]->dof[i] > 0) {
      memset(numidofs, 0, sizeof(int) * nnodes);
      printf("     ");fprint_line(stdout,"=", 75); printf("\n");
      printf("     "); 
      printf("VAR_TYPE = %-d:", i);
      if (i == 4) {
	printf(" note-> only subvar=0 dofs are included below");
      }
      printf("\n");
      printf("        ");
      printf("idof  lnn  ledof lvdesc,lvdesc_dof VT_offset gnn  "
	     " gun row_lvdof active\n");
      printf("        "); fprint_line(stdout,"-", 70);
      for (j = 0; j < ei[pg->imtrx]->dof[i]; j++) {
	printf("        ");   
	printf("  %-5d", j);
	ln =  ei[pg->imtrx]->dof_list[i][j];
	numidofs[ln]++;
	printf("%-5d ", ln);
	ledof = ei[pg->imtrx]->lvdof_to_ledof[i][j];
	printf("%-5d ", ledof);
	printf(" %-5d ", ei[pg->imtrx]->lvdof_to_lvdesc[i][j]);
	printf("  %-5d ", ei[pg->imtrx]->lvdof_to_lvdesc_dof[i][j]);
	printf("     %-5d ", ei[pg->imtrx]->Baby_Dolphin[i][j]);
	printf(" %-6d", ei[pg->imtrx]->gnn_list[i][j]);
	printf(" %-6d", ei[pg->imtrx]->gun_list[i][j]);
	printf(" %-5d", ei[pg->imtrx]->lvdof_to_row_lvdof[i][j]);
	printf(" %-5d", ei[pg->imtrx]->active_interp_ledof[ledof]);
	printf("\n");
      }
      printf("        "); fprint_line(stdout,"-", 70);

      printf("        ");
      printf("lnn Numidofs  lvdof  lvdof_first \n");
      printf("        "); fprint_line(stdout,"-", 60);
      for (ln = 0; ln < nnodes; ln++) {
	printf("        ");   
	printf(" %-5d", ln);
	printf("  %-5d ", numidofs[ln]);
	printf("  %-5d ", ei[pg->imtrx]->ln_to_dof[i][ln]);
	printf("  %-5d ", ei[pg->imtrx]->ln_to_first_dof[i][ln]);
	printf("\n");
      }
      printf("        "); fprint_line(stdout,"-", 60);
    }
  }
  fprint_line(stdout,"+", 80);

  safer_free((void **) &numidofs);
}
#endif
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
