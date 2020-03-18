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
 *$Id: rf_node.c,v 5.3 2008-05-02 19:07:57 hkmoffa Exp $
 */

#include <stdio.h>
#include <stdlib.h>

#include "std.h"
#include "rf_fem.h"
#include "el_elm.h"
#include "el_geom.h"
#include "exo_struct.h"
#include "rf_vars_const.h"
#include "mm_as_structs.h"
#include "rf_node_const.h"
#include "rf_mp.h"
#include "rf_allo.h"
#include "dpi.h"
#include "el_elm_info.h"
#include "mm_as.h"
#include "mm_eh.h"
#include "rf_bc.h"
#include "rf_bc_const.h"

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

int
bin_search_max(const int list[], const int num, const int value)

    /********************************************************************
     *
     *  bin_search_max():
     *
     *     Searches a monotonically increasing list of values for the value,
     * value.
     * It returns the index of the first position found, which matches value.
     * The list is assumed to be monotonically increasing, and
     * consist of elements list[0], ..., list[n-1].
     *   If no position in list matches value, the next highest position
     * in the list, or 0 if the value is lower than all values in list []
     * is returned. If the value is greater than all elements of the
     * list, then the value of num is returned.
     * Error conditions are specified with a return value of -1.
     *
     *   example:
     * -------------
     *    list [] = {3, 5, 7, 9}, num = 4
     *
     *       value       return_value
     *      ---------   --------------
     *          1             0
     *          3             0
     *          6             2
     *          8             3
     *          9             3
     *          13            4
     *
     **********************************************************************/
{
  int top, bottom, diff, middle, val_mid, val_top;
  if (num <= 0) return (-1);
  bottom = 0;
  if (value <= list[0]) return 0;
  top  = num - 1;
  val_top  = list[top];
  if (value >=  val_top) {
    if (value > val_top) return top + 1;
    return top;
  } 
  while ( (diff = top - bottom) > 1) {
    middle = bottom + diff/2;
    val_mid = list[middle];
    if (value > val_mid) {
      bottom = middle;
    } else if (value < val_mid) {
      top = middle;
    } else {
      return middle;
    }
  }
  return top;
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

void add_to_umi_int_list(UMI_LIST_STRUCT *ls, int addition)

    /********************************************************************
     *
     *  add_to_umi_int_list
     *
     *   This code snipet will add an entry to a Unique Monotonically
     *   Increasing (umi) structure list of integers.
     *
     *
     *  Input
     * --------
     *  ls = Pointer to an UMI_List structure. it is assumed that
     *       if this structure is virgin, then all of its members
     *       are zero.
     *  addition = Number to be added to the list. Note, it will be
     *             added to the list only if it is unique. Also,
     *             it will be added in the correct slot so as to
     *             
     *
     *  It calls malloc and realloc directly.
     *
     *
     *  Note: These UMI structs are really meant for short lists.
     *        The realloc and copy operations done in this routine
     *        don't scale well to high numbers. However, that's exactly
     *        what we are doing in the current applications -> sorting
     *        of lists with less than 10 entries ))). 
     ********************************************************************/
{
  int list_loc, i;
  int *list = ls->List;
  if (list == NULL) {
    //ls->List = (int *) smalloc(sizeof(int));
    ls->List = alloc_int_1(1, INT_NOINIT);
    if (ls->List == NULL) {
      printf("Error exit from add_to_umi_int_list\n");
      exit(-1);
    }
    ls->List[0] = addition;
    ls->Length = 1;
  } else {
    list_loc = bin_search_max(ls->List, ls->Length, addition);
    if (list_loc < ls->Length) {
      if (list[list_loc] == addition) return;
    }
    ls->Length++;
    //list = realloc(list, sizeof(int)*(ls->Length));
    realloc_int_1(&list, ls->Length, ls->Length-1);
    if (list == NULL) {
      printf("Error exit from add_to_umi_int_list()\n");
      exit(-1);
    }
    ls->List = list;
    for (i = ls->Length - 1; i > list_loc; i--) {
      list[i] = list[i-1];
    }
    list[list_loc] = addition;
  }
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

void
init_nodes (Exo_DB *exo, Dpi *dpi)

    /*****************************************************************
     *
     * init_nodes()
     *
     *
     *  Initializes the information in the Nodes_Info structure, Nodes.
     *  and related structures, like Variable_Masks.
     *
     *  -> Flag all nodes on side sets with the EDGE bit
     *  -> Flag all nodes on discontinuous sides with the DISC_BNDRY bit.
     *
     *****************************************************************/
{
  int i, ss_index, side_index, node_num, total_nodes, k, ibc, where;
  NODE_INFO_STRUCT *tptr;
  int number_proc_nodes = exo->num_nodes;
  struct BC_descriptions *bc_desc;
  
  /*
   * Calculate the total number of nodes and do a consistency check
   * with the fe database.
   */
  total_nodes = Num_Internal_Nodes + Num_Border_Nodes
      + Num_External_Nodes;
  if (total_nodes != number_proc_nodes) {
    fprintf(stderr,
  "Proc %d : init_nodes ERROR total nodes, %d, differs from num_proc_nodes, %d\n",
	    ProcID, total_nodes, number_proc_nodes);
    EH(-1, "init_nodes: Shouldn't be here");
  }

  /*
   *  We initialize in two steps here fore compatibility with MPSalsa
   *  where nodes structs are actually moved around dynamically.
   */
  Nodes = (NODE_INFO_STRUCT **) alloc_ptr_1(number_proc_nodes);
  tptr = alloc_struct_1(NODE_INFO_STRUCT, number_proc_nodes);

  for (i = 0; i < number_proc_nodes; i++) 
     {
      /*
       *  Point to one Node_Info structure
       *  for the current node.
       */
      Nodes[i] = tptr + i;
      /*
       *  None of this has been made compatible with parallel processing
       */
      Nodes[i]->Proc_Node_Num = i;
      Nodes[i]->First_Unknown = (int *) malloc( (upd->Total_Num_Matrices) * sizeof(int) );
      Nodes[i]->Global_Node_Num = dpi->node_index_global[i];
      Nodes[i]->Type.Internal = FALSE;
      Nodes[i]->Type.Border  = FALSE;    
      Nodes[i]->Type.External = FALSE;    
      Nodes[i]->Proc = ProcID;
      Nodes[i]->Nodal_Vars_Info = (NODAL_VARS_STRUCT **) alloc_ptr_1(upd->Total_Num_Matrices);
      Nodes[i]->DBC = malloc(sizeof(short int *) * upd->Total_Num_Matrices);

      for (k = 0; k < upd->Total_Num_Matrices; k++) {
        Nodes[i]->DBC[k] = NULL;
      }

      if (i < Num_Internal_Nodes) 
        {
         Nodes[i]->Type.Internal = TRUE;
         Nodes[i]->Type.Owned = TRUE;       
        } 
      else if (i < Num_Internal_Nodes + Num_Border_Nodes) 
        {
         Nodes[i]->Type.Border = TRUE;
         Nodes[i]->Type.Owned = TRUE;
        } 
      else 
        {
         Nodes[i]->Type.External = TRUE;
         Nodes[i]->Type.Owned = FALSE;
         where = dpi->ptr_set_membership[i];
         Nodes[i]->Proc = dpi->set_membership[where];
        }
     }

  /*
   * Make a list of what materials each node belongs to
   */
  make_node_to_elem_matrl_map(exo);

  /*
   *      loop through all side sets to determine which has
   *      a boundary
   */
 
  for (ss_index = 0; ss_index < exo->num_side_sets; ss_index++) {
    for (side_index = 0; side_index < exo->ss_num_sides[ss_index];
	 side_index++) {
      for (k = exo->ss_node_side_index[ss_index][side_index];
	   k < exo->ss_node_side_index[ss_index][side_index+1]; k++) {
	node_num = exo->ss_node_list[ss_index][k];

	/*
	 *  Flag nodes with side sets using the EDGE bit field
	 */
	Nodes[node_num]->EDGE = 1;

	/*
	 * Search for nodes that have multivalued, Discontinuous BCs 
	 * Flag these using the DISC_BNDRY bit field
	 */
	for (ibc = 0; ibc < Num_BC; ibc++) {
	  if (BC_Types[ibc].BC_ID == exo->ss_id[ss_index]) {
            bc_desc = BC_Types[ibc].desc;
            if (bc_desc->i_apply == CROSS_PHASE_DISCONTINUOUS) {
	      Nodes[node_num]->DISC_BNDRY = 1;
	    }
	  }
	}
      }
    }
  }
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

void
free_nodes (void)

    /********************************************************************
     *
     * free_nodes():
     *
     * Operation to free the Nodes structure and all of its underlying
     * substructures. See init_nodes() for a description of the way
     * it was malloced in the first place.
     ********************************************************************/
{
  int i,k;
  NODE_INFO_STRUCT *node_ptr;
  for (i = 0; i < Num_Node; i++) {
    node_ptr = Nodes[i];
    free_umi_list(&(node_ptr->Mat_List));
    for (k = 0; k < upd->Total_Num_Matrices; k++) {
      safer_free((void **) &(node_ptr->DBC[k]));
    }
    free(node_ptr->First_Unknown);
    free(node_ptr->Nodal_Vars_Info);
    free(node_ptr->DBC);
  }
  /*
   *  free_umi_list(&(node_ptr->Element_List));
   */
  safer_free((void **) &(Nodes[0]));
  safer_free((void **) &(Nodes));
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

void
make_node_to_elem_matrl_map(Exo_DB *exo)	      
    
  /*********************************************************************
   *
   * make_node_to_elem_matrl_map():
   *
   * This function creates a node to material mapping 
   *********************************************************************/
{
  int eb_index, mn, e_start, e_end, ielem, ielem_type, num_local_nodes;
  int iconnect_ptr, i, I;
#ifdef DEBUG_NODE
  UMI_LIST_STRUCT *matrlLP;
#endif

  /*
   *  Loop over the element blocks
   */
  for (eb_index = 0; eb_index < exo->num_elem_blocks; eb_index++)  {

    /*
     *  Identify the materials index with the element block
     */
    mn = Matilda[eb_index];
    if (mn < 0) {
      continue;
    }

    /* 
     * Find the beginning and end element number for the elements
     *  in this element block. Note, they are all numbered consequatively
     *  on each processor (right ???!)
     */
    e_start = exo->eb_ptr[eb_index];
    e_end   = exo->eb_ptr[eb_index+1];

    /*
     *  Look up the element type for the element block
     */
    ielem_type = exo->eb_elem_itype[eb_index];
    
    /*
     *  based on the element type, get the number of quad points,
     *  number of local nodes, and the element dimension
     */
    num_local_nodes = elem_info(NNODES, ielem_type); 
    
    /*
     *  Loop over all elements in the element block
     */
    for (ielem = e_start; ielem < e_end; ielem++) {

      /*
       *  find ptr to the beginning of this element's connectivity list
       */
      iconnect_ptr    = Proc_Connect_Ptr[ielem];

      for ( i = 0; i < num_local_nodes; i++) {
	I = Proc_Elem_Connect[iconnect_ptr + i];
        add_to_umi_int_list(&(Nodes[I]->Mat_List), mn);
	/*
	  add_to_umi_int_list(&(Nodes[I]->Element_List), ielem);
	*/
      }
    }
  }
#ifdef DEBUG_NODE
  print_sync_start(TRUE);
  printf("Processor %d: Num nodes = %d\n", ProcID, exo->num_nodes);
  for (I = 0; I < exo->num_nodes; I++) {
    matrlLP = &(Nodes[I]->Mat_List);
    printf(" %5d  %5d : ", I, matrlLP->Length);
    for (i = 0; i < matrlLP->Length; i++) {
      printf(" %d ", matrlLP->List[i]);
    }
    printf("\n");
  }
  fflush(stdout);
  print_sync_end(TRUE);
#endif
  return;
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

int
index_umi_list(const UMI_LIST_STRUCT *umi_ptr, const int value)
    
   /*********************************************************************
    *
    * index_umi_list
    * 
    * This function returns the index of value in a umi structure
    * If the value doesn't exist, it return -1.
    *
    *  Input
    * --------
    *  umi_ptr : pointer to the umi list
    *  value   : value to look up in the list
    *
    *  Return
    * ---------
    *  Index position in the list that contains the entry, value
    *  Or -1 if the list doesn't contain thevalue.
    *********************************************************************/
{
  int list_loc;
  if (umi_ptr == NULL) return -1;
  if (umi_ptr->Length <= 0) return -1;
  list_loc = bin_search_max(umi_ptr->List, umi_ptr->Length, value);
  if (list_loc < umi_ptr->Length) {
    if (umi_ptr->List[list_loc] == value) return list_loc;
  }
  return -1;
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

int
node_matrl_index(const int nodeNum, const int mn)
    
   /*********************************************************************
    *
    * node_matrl_index():
    * 
    * This function returns the index of the material in a global nodes
    * list of materials. If the node is not in a material, it returns
    * -1. The first material for a node will have an index value of 0.
    *
    *  Input
    *  ------
    *  nodeNum = local node number on this processor
    *  mn      = Index id of the material
    *********************************************************************/
{
  UMI_LIST_STRUCT *umi_ptr = &(Nodes[nodeNum]->Mat_List);
  int index = index_umi_list(umi_ptr, mn);
  return index;
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

int
first_matID_at_node(const int nodeNum)
    
   /*********************************************************************
    *
    * first_matID_at_node:
    * 
    * This function returns the material index of the first material
    * located at a processor node number, nodeNum. This will always be
    * the lowest number material index at that node.
    *
    *  Input
    *  ------
    *  nodeNum =  node number on this processor
    *********************************************************************/
{
  UMI_LIST_STRUCT *umi_ptr = &(Nodes[nodeNum]->Mat_List);
  if (umi_ptr == NULL || umi_ptr->Length <= 0) {
    fprintf(stderr,
     "first_matID_at_node ERROR P_%d: mat list is empty at node %d\n",
	    ProcID, nodeNum);
    EH(-1,"first_matID_at_node ERROR");
    return -1;
  }
  return umi_ptr->List[0];
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

void
free_umi_list(UMI_LIST_STRUCT *umi_ptr)
    
    /*********************************************************************
    *
    * free_umi_list():
    * 
    * This function frees the underlying umi (unique monotomically
    * increasing) list structure and set the length variable to zero.
    *********************************************************************/
{
 if (umi_ptr == NULL) return;
 safer_free((void **) &(umi_ptr->List));
 umi_ptr->Length = 0;
}

/************************************************************************/
/************************************************************************/
/************************************************************************/

void
node_info_tmp_free (void)

    /********************************************************************
     *
     *  nodal_resid_tmp_free()
     *
     *     This routine frees all malloced memory that was created as
     *  temporary memory during a residual and/or Jacobian fill 
     *  operation. Freed pointers have their values set to NULL, as per
     *  our convention.
     *******************************************************************/
{
  int I;
  NODE_INFO_STRUCT *node;
  /*
   *  Loop over all of the nodes, whether owned or not, calling 
   *  a function to destroy the memory
   */
  for (I = 0; I < EXO_ptr->num_nodes; I++) {
    node = Nodes[I];
    if (node->Resid_Wksp) {
      nodal_resid_wksp_destroy(&(node->Resid_Wksp));
    }
  }
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

NODAL_RESID_WKSP_STRUCT *
nodal_resid_wksp_alloc (void)

    /********************************************************************
     *
     *  nodal_resid_wksp_alloc()
     *
     *    This function allocates a Nodal_Resid_Wksp structure and then
     *    zeroes all of its members. The pointer to the function is
     *    returned.
     *******************************************************************/
{
  NODAL_RESID_WKSP_STRUCT *ptr;
  ptr = alloc_struct_1(NODAL_RESID_WKSP_STRUCT, 1);
  return ptr;
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

void
nodal_resid_wksp_destroy(NODAL_RESID_WKSP_STRUCT **residWkspHdl)

    /********************************************************************
     *
     *  nodal_resid_wksp_destroy()
     *
     *    Given a handle to a Nodal_Resid_Wksp structure, this routine
     * will free both the underlying memory for the structure, and the
     * structure itself.
     *******************************************************************/
{
  nodal_resid_wksp_free(*residWkspHdl);
  safer_free((void **) residWkspHdl);
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

void
nodal_resid_wksp_free(NODAL_RESID_WKSP_STRUCT *residWksp)

    /********************************************************************
     *
     *  nodal_resid_wksp_free()
     *
     *     Given a pointer to a nodal_resid_wksp structure, this routine
     *  will take care of freeing the underlying memory assocated with
     *  that structure.
     *******************************************************************/
{
  VCRR_STRUCT *next;
  /*
   *  Look at each pointer entry to see if a subroutine needs to be
   *  called to clear it.
   *
   */
  for (next = *(residWksp->Vcrr); next == NULL; next++) {
    safer_free((void **) &next);
  }
  safer_free((void **) &(residWksp->Vcrr));
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

void
vcrr_add(VCRR_STRUCT ***vcrr_hdl, const int eqn, const int mn_curr,
	 const int mn_apply)

    /*******************************************************************
     *
     * vcrr_add()
     *
     *    Add an entry onto the list of Volumetric Contrinbution Row
     *    Redirect entries for the current node.
     *    We do this by first making sure that the entry doesn't exist
     *    Then, we add the entry onto the end of the list by first
     *    lengthening the size of the pointer array and then mallocing
     *    space for the entry itself. Each entry consists of the
     *    following:
     *
     *       eqn: -> Variable type
     *       mn_curr -> material id of the variable whose residual
     *                  contribution will be redirected.
     *       mn_apply -> The ending material ID of the redirection.
     *
     *
     *  In other words, the residual contribution from the continuition
     *  equation for the variable (eqn, mn_curr) 
     *  will end up in the slot for the continuity equation for the
     *  variable (eqn, mn_apply).
     *
     *******************************************************************/
{
  int num = 0, nf = TRUE, ft = 0;
  VCRR_STRUCT **vcrr = *vcrr_hdl, **v;
  VCRR_STRUCT *vcrr_ptr;
  if (vcrr) {
    num++;
    for (v = vcrr; (*v != NULL) && nf; v++) {
      vcrr_ptr = *v;
      num++;
      if (vcrr_ptr->Var_Type == eqn) {
	if (vcrr_ptr->MatID_Volume == mn_curr) {
          if (vcrr_ptr->MatID_Row_Redirect == mn_apply) {
	    nf = FALSE;
	  } else {
	    EH(-1, "incompatibility");
	  }
	}
      }
    }
  } else {
    ft = 1;
  }
  if (nf) {
    realloc_ptr_1((void ***) &vcrr, num + 1 + ft, num);
    if (num > 0) num--;
    vcrr[num] = alloc_struct_1(VCRR_STRUCT, 1);
    vcrr_ptr = vcrr[num];
    vcrr_ptr->Var_Type = eqn;
    vcrr_ptr->MatID_Volume = mn_curr;
    vcrr_ptr->MatID_Row_Redirect = mn_apply;
  }
  *vcrr_hdl = vcrr;
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

int
vcrr_lookup(VCRR_STRUCT **vcrr, const int eqn, const int mn_curr)

    /********************************************************************
     *
     * vcrr_lookup()
     *
     *   Lookup whether a particular variable identified by a
     *   variable type, material ID pair (eqn, mn_curr), has a
     *   vcrr entry. Return the mn_apply value if it does. If it 
     *   doesn't have an entry, return mn_curr, indicating tthat the
     *   volumetric contribution to the total continuity equation for
     *   that variable stays in its own slot. 
     *
     *  If mn_apply isn't equal to mn_curr, then the residual
     *  contribution from the continuition equation for the variable
     *  (eqn, mn_curr) will end up in the slot for the continuity 
     *  equation for the variable (eqn, mn_apply).
     ********************************************************************/
{
  int nf;
  VCRR_STRUCT **v;
  VCRR_STRUCT *vcrr_ptr;
  if (vcrr) {
    for (v = vcrr, nf = TRUE; (*v != NULL) && nf; v++) {
      vcrr_ptr = *v;
      if (vcrr_ptr->Var_Type == eqn) {
	if (vcrr_ptr->MatID_Volume == mn_curr) {
          return vcrr_ptr->MatID_Row_Redirect;
	}
      }
    }
  } 
  return mn_curr;
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

void
nullify_dirichlet_bcs(void) 

    /**************************************************************
     *
     * nullify_dirichlet_bcs()
     *
     * Initialize Variable_Mask to say, "No BCs have been set anywhere."
     **************************************************************/
{
  int i, v;
  int total_nodes;
  NODE_INFO_STRUCT *node;
  NODAL_VARS_STRUCT *nv;
  total_nodes = Num_Internal_Nodes + Num_Border_Nodes + Num_External_Nodes;
  for (i = 0; i < total_nodes; i++) {
    node = Nodes[i];
    if (node->DBC[pg->imtrx]) {
      nv = node->Nodal_Vars_Info[pg->imtrx];
      for (v = 0; v < nv->Num_Unknowns; v++) {
	node->DBC[pg->imtrx][v] = -1;
      }
    }
  }
}
/************************************************************************/
/************************************************************************/
/************************************************************************/
