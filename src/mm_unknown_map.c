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
 *$Id: mm_unknown_map.c,v 5.13 2010-04-07 22:27:00 prschun Exp $
 */

/*
 * Notes:
 * ------
 *	[1] Modifications made to handle variables with different degrees of
 *	    freedom per element. Added Local_Offset arrays to help find the
 *	    global unknown number for variables in these kinds of elements.
 *	    This was done so that (v,P) = (Q2,P1) would be easily implemented,
 *	    etc.
 *
 * Modified:  Wed Feb 16 06:20:52 MST 1994 pasacki@sandia.gov
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "dp_map_comm_vec.h"
#include "dp_types.h"
#include "dp_utils.h"
#include "dpi.h"
#include "el_elm.h"
#include "el_elm_info.h"
#include "el_geom.h"
#include "exo_struct.h"
#include "mm_as.h"
#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "mm_eh.h"
#include "mm_mp.h"
#include "mm_mp_structs.h"
#include "mm_unknown_map.h"
#include "rd_mesh.h"
#include "rf_allo.h"
#include "rf_bc_const.h"
#include "rf_fem.h"
#include "rf_fem_const.h"
#include "rf_io.h"
#include "rf_mp.h"
#include "rf_node_const.h"
#include "rf_util.h"
#include "rf_vars_const.h"
#include "std.h"

#ifdef PARALLEL
#include "mpi.h"
#endif

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

/*
 * Definitions. These variables are used many other places via the
 * extern declarations in rf_fem.h
 */

int *num_internal_dofs;
int *num_boundary_dofs;
int *num_external_dofs;
int *num_universe_dofs;
int num_personal_elems = 0;

/*
 * Local_Offset:
 *
 *	This array of pointers gets set to lists of integers that, for every
 *	node, indicates how much offset from the first unknown of the node
 *	you need to get to the variable of interest...
 *
 *	For example, if we're solving just the energy equation and temperature
 *	is interpolated at every node, then
 *
 *		Local_Offset[node][TEMPERATURE] = 0; (the first unknown is T)
 *	and
 *		Local_Offset[node][VELOCITY1] = -1; (undefined offset)
 *
 *	Yes, this does duplicate some functionality of First_Y, First_MeshD,
 *	etc. Also, Index_P, will not really be needed anymore, since pressure
 *	is getting lumped together with other unknowns at a node.
 *
 */
int ***Local_Offset = NULL;

/*
 * Dolphin:  The scheme above is fine except for cases that occur when
 *	     neighboring elements that share a node have different ideas
 *	     about which variables and degrees of freedom need to be solved
 *	     at the node.
 *
 *	     In order to get the true accounting of offsets for each type of
 *	     variable at a node, we need to know what every element expects
 *	     of the node.
 *
 *
 *	     Thus, poll every element and save the largest estimate for
 *	     the degrees of freedom required to represent each variable...
 *
 *	     Dolphin[matrix][node][variable] = dof
 *
 * where:
 *		node  == node number
 *
 *		var   == variable index
 *
 *		dof   == number of degrees of freedom for this type of
 *			 variable at this node.
 *
 *			 Note that variables with multiple k types must all
 *			 have the same representation in terms of dof/node.
 *			 Thus, the number of degrees of freedom for each
 *			 species concentrations must be multiplied by the
 * 			 total number of concentrations that are active, too!
 */
int ***Dolphin = NULL;

/*
 * Other Globally Defined Variables
 */
int *NumUnknowns;      /* Number of unknown variables updated by this   */
                       /* processor (internal plus boundary)            */
int *NumExtUnknowns;   /* Number of external variables which are        */
                       /* copied and stored on the local processor      */
int MaxVarPerNode = 0; /* Global Maximum number of unknowns at any      */
                       /* node on any processor                         */
int Num_Var_In_Type[MAX_NUM_MATRICES][MAX_VARIABLE_TYPES];
/*
 * First_Unknown:
 *     Index to start of the first unknown in the processor
 *     solution vector at each node at each matrix.
 *        Dimension = [upd->Total_Num_Matrices][total_num_nodes];
 */
int **First_Unknown = NULL;

/*
 * dofname -- holds strings telling the name of the variable (u1, T, P, etc)
 *            and the global node number associated with a particular gdof.
 *            This should aid in debugging, etc. Allocation and setup in
 *	      mm_unknown_map.c.
 *
 * Now, idv[matrix][dof][0] = VELOCITY1, etc.
 *      idv[dof][1] = local nodal dof (0, except pressure& conc., for example)
 *	idv[dof][2] = associated global node number (0-based)
 */
int ***idv = NULL;      /* Integer variable name, nodal dof, node. */
char ***dofname = NULL; /* Names of variables. */
char ***resname = NULL; /* Names of residual equations. */

/*
 * Definitions of variables which are declared in
 * rf_mask.h
 */
int Inter_Mask[MAX_NUM_MATRICES][MAX_VARIABLE_TYPES][MAX_VARIABLE_TYPES] = {{{0}}};
int Ignore_Deps[MAX_NUM_MATRICES][MAX_VARIABLE_TYPES][MAX_VARIABLE_TYPES];

/*
 * Prototypes declarations of static functions defined in this file.
 */
static void set_interaction_masks(Exo_DB *);

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void setup_local_nodal_vars(Exo_DB *exo, Dpi *dpi)

/********************************************************************
 *
 * Routine to produce the map of the (node,variable) -> solution
 * vector. This routine is to produce the expanded solution vector
 * and the necessary arrays of pointers into this vector.
 *
 * Output:
 *
 * Variable array mask: (see rf_masks.h)
 *
 * Variable_Mask[global_node]     --	integer structure. Which contains
 * 					information on which variables are
 *   				        defined at a given node. In this
 *					routine the unknown info is set.
 *					The structure also contains info on
 *					the boundary conditiones which are
 *					set in find_and_set_Dirichlet.
 *
 *  Equation                      --   is a structure of integers defined in
 * 					"rf_masks.h" which indicates which
 *					equations are being set up and
 *					solved for.
 *
 *  Num_Unknowns_Node[global_node] --	entry is the total number of unknowns
 *  					which are defined at this node.
 *
 *  NumUnknowns                   -- 	total number of unknowns solved for
 *					by this Proc.
 *
 *  NumExtUnknowns                --	total number of unknowns which are
 *					stored on this processor and
 *					communicated from neighboring processors
 *
 *  MaxVarPerNode 		  --	value is equal to the maximum number
 *					of variables defined at any node on
 *                                      any processor.
 *  First_Unknown[matrix][global_node]     --	index into solution vector for which
 *  					        the first unknown at local node i
 *
 * Local_Offset[matrix][global_node][variable] -- how much offset from the first
 *					unknown to the variable of interest
 *
 ******************************************************************************/
{
  int nun[MAX_NUM_MATRICES][MAX_VARIABLE_TYPES]; /* Number of unknowns per node for each */
                                                 /* variable. Note that if there are any */
                                                 /* concentration variables, that they are */
                                                 /* all represented with the same ndof/node. */

  int e, e_start, e_end, i, j, n, index, num_nodes, ebi, elementType;
  int total_nodes, var_type;
  int mn;
  int imtrx;
  NODAL_VARS_STRUCT **tmp_nodal_vars = NULL, *nv, *nv_match;
  VARIABLE_DESCRIPTION_STRUCT *var_ptr;

  /*
   * Calculate the total number of nodes on this processor
   */

  total_nodes = Num_Internal_Nodes + Num_Border_Nodes + Num_External_Nodes;

  /*
   *  Allocate memory for the variable mask array of structures and
   *  other pointer arrays into the solution vector. The Variable_Mask
   *  structure is an array of bit fields.
   */
  First_Unknown = (int **)malloc((upd->Total_Num_Matrices) * sizeof(int *));
  for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
    First_Unknown[imtrx] = (int *)malloc(total_nodes * sizeof(int));
  }
  for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
    for (i = 0; i < total_nodes; i++) {
      First_Unknown[imtrx][i] = -1;
    }
  }
  /*
   *  Allocate temporary memory for holding the Nodal vars information
   *  for this processor. It will be released at the end of the routine.
   */
  tmp_nodal_vars = (NODAL_VARS_STRUCT **)alloc_ptr_1(upd->Total_Num_Matrices);
  for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
    tmp_nodal_vars[imtrx] = alloc_struct_1(NODAL_VARS_STRUCT, total_nodes);
    for (i = 0; i < total_nodes; i++) {
      tmp_nodal_vars[imtrx][i].list_index = i;
    }
  }
  /*
   * Introducing some new arrays for your programming pleasure...
   *
   * Difference between a Dolphin and a Local_Offset?
   *
   * The Dolphin keeps track of the true count of how many dof's for each kind
   * of variable are at each node. Dolphin does *not* account for
   * multiplicity of species concentrations or any other variable.
   *
   * The Local_Offset is an *accumulation* that helps to point to any
   * particular kind of variable. Local_Offset *is* adjusted to account for
   * multiplicity of species concentrations.
   *
   */

  Dolphin = (int ***)malloc((upd->Total_Num_Matrices) * sizeof(int **));
  Local_Offset = (int ***)malloc((upd->Total_Num_Matrices) * sizeof(int **));
  for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
    Dolphin[imtrx] = (int **)malloc(total_nodes * sizeof(int *));
    Local_Offset[imtrx] = (int **)malloc(total_nodes * sizeof(int *));
    for (i = 0; i < total_nodes; i++) {
      Dolphin[imtrx][i] = (int *)malloc(MAX_VARIABLE_TYPES * sizeof(int));
      Local_Offset[imtrx][i] = (int *)malloc(MAX_VARIABLE_TYPES * sizeof(int));
    }
  }
  for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
    for (i = 0; i < total_nodes; i++) {
      for (var_type = V_FIRST; var_type < V_LAST; var_type++) {
        Dolphin[imtrx][i][var_type] = -1;
        Local_Offset[imtrx][i][var_type] = -1;
      }
    }
  }

  /*
   * During this first loop over all elements we try to get
   * the overall picture accurate in terms of how many unknowns need
   * to be represented at each node. We need to do this for the benefit
   * of nodes shared by different elements with different variables
   * represented within the neighboring elements...
   */

  /*
   * For the distributed processing version, note that this step will
   * accumulate some representation of the number and kinds of active
   * variables at external nodes that are actually owned by other processors.
   *
   * That's OK, as long as we understand that the owning processor may
   * well have a different idea of the number and kinds of variables that
   * are active at that node. This processor's idea of the number of active
   * variables at that node will be a subset of the actual list of active
   * variables. Only in the optimistic case will this processor's concept
   * of the active variables at an external node be identical with the
   * nodes true count of dofs in a global context. That true count will
   * only be met with certainty on the processor that owns that node
   * during this initial step. In subsequent steps, the owning node will
   * tell ghost nodes exactly what unknowns exist at the node.
   */

  /*
   * Loop through element blocks on this processor - there may well be
   * fewer than the number of materials globally.
   */
  for (ebi = 0; ebi < exo->num_elem_blocks; ebi++) {
    e_start = exo->eb_ptr[ebi];
    e_end = exo->eb_ptr[ebi + 1];

    /*
     * Determine the material index for this particular element block
     * (if Matilda[ebi] is negative, then the element block is not in the goma problem)
     */
    mn = Matilda[ebi];
    if (mn < 0) {
      continue;
    }
    pd = pd_glob[mn];
    mp = mp_glob[mn];

    for (e = e_start; e < e_end; e++) {

      /*
       * Store the element type for this element
       */
      elementType = Elem_Type(exo, e);
      /*
       * Find out how many local nodes are in this element type
       */
      num_nodes = elem_info(NNODES, elementType);

      /*
       * For this particular element, e, find the offset into the
       * connectivity matrix, index. Index will be used to find the value
       * of the processor node number.
       */
      index = Proc_Connect_Ptr[e];

      /*
       *  Loop over the local element nodes in the current element
       */
      for (n = 0; n < num_nodes; n++) {

        /*
         *  Find the processor node number, i, corresponding to the local
         *  node number, n, by indexing into the processor element
         *  connectivity array, Proc_Elem_Connect[].
         *  (note: i may or may not be owned by this processor)
         */
        i = Proc_Elem_Connect[index++];

        /*
         * For each kind of variable type in the problem, depending on its
         * interpolation, decide how many, if any, degrees of freedom
         * it has at this node...
         */
        for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
          for (var_type = V_FIRST; var_type < V_LAST; var_type++) {

            /*
             *  Find the number of degrees of freedom at the local node
             *  corresponding to the current element for the current
             *  variable type given the problem description corresponding
             *  to the current material index.
             */
            nun[imtrx][var_type] =
                dof_lnode_var_type(n, elementType, i, var_type, pd_glob[mn], imtrx);
            if (nun[imtrx][var_type] > 0) {
              /*
               * Find a variable structure that matches this one
               */
              var_ptr = find_or_create_vd(var_type, nun[imtrx][var_type], mn, 0, imtrx);

              /*
               * Add this variable to the temporary list of variables
               * defined at this node
               */
              add_var_to_nv_struct(var_ptr, tmp_nodal_vars[imtrx] + i);

              /*
               * special case section for MASS_FRACTION variable type
               * Add a variable description structure for each of the
               * subvariables specified for the current material.
               * -> HKM: this will go away eventually
               */
              if (var_type == MASS_FRACTION) {
                for (j = 1; j < mp->Num_Species_Eqn; j++) {
                  var_ptr = find_or_create_vd(var_type, nun[imtrx][var_type], mn, j, imtrx);
                  add_var_to_nv_struct(var_ptr, tmp_nodal_vars[imtrx] + i);
                }
              }
            }

            /*
             * Calculate the old way:
             *  the calculation depends on pd for its interpolation type
             *  Also, discontinuous variables are treated here as well.
             *  sometimes leading to a doubling of the dofs at a node.
             *  The new way also gets to a doubling of the dof at the node.
             *  However, it's done by looping through different elements
             *  and then comparing the MatID fields to determine uniqueness.
             *  This is done in the loop below. Therefore, a printout from
             *  below is just informational.
             */
          } /* Loop over variables list */
        }   /* Loop over matrices */
      }     /* Loop over local nodes in an element */
    }       /* loop over elements in element blocks */
  }         /* loop over materials */

  /*
   *  Loop though nodes this time calculating the real nodal_vars
   *  structures for each node.
   */
  for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
    for (i = 0; i < total_nodes; i++) {
      nv = copy_and_rearrange_nv(tmp_nodal_vars[imtrx] + i);

      /*
       * Decide whether this nodal variable structure is new.
       * Destroy it if it isn't. If it unique, it gets added to
       * a list of unique nodal variable structures.
       */
      if (MatchAdd_Node_Vars_List(nv, &nv_match)) {
        nodal_vars_destroy(&nv);
      }

      /*
       * Copy the pointer to the node vars structure into the node_info
       * structure for the current node.
       */
      Nodes[i]->Nodal_Vars_Info[imtrx] = nv_match;

      /*
       * Do a consistency check between the new and old way of
       * calculating the degrees of freedom at a node
       */
      for (var_type = V_FIRST; var_type < V_LAST; var_type++) {
        ebi = get_nv_ndofs_modMF(nv_match, var_type);
        /*
         * Now populate the Dolphin array with the new way
         * -> The DEBUG block above has proved that old vs new
         *    is the same.
         */
        if (ebi > 0)
          Dolphin[imtrx][i][var_type] = ebi;
      }
    }
  }
  /*
   *  Unmalloc the tmp_nodal_vars variable and all of the underlying
   *  malloced structures
   */
  for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
    for (i = 0; i < total_nodes; i++) {
      nodal_vars_free(tmp_nodal_vars[imtrx] + i);
    }
    safer_free((void **)&(tmp_nodal_vars[imtrx]));
  }
  free(tmp_nodal_vars);

  /*
   *  When in debug mode, print out a complete listing of variables at
   *  every node
   */
}
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void print_vars_at_nodes(void)

/************************************************************************
 *
 * print_vars_at_nodes();
 *
 *
 *************************************************************************/
{
  int i, j, total_nodes, maxLength;
  int imtrx;
  NODE_INFO_STRUCT *node_ptr;
  NODAL_VARS_STRUCT *nv;
  VARIABLE_DESCRIPTION_STRUCT *vd;
  UMI_LIST_STRUCT *list;
  /*
   * Calculate the total number of nodes on this processor
   */
  total_nodes = Num_Internal_Nodes + Num_Border_Nodes + Num_External_Nodes;

  maxLength = -1;
  for (i = 0; i < total_nodes; i++) {
    node_ptr = Nodes[i];
    list = &(node_ptr->Mat_List);
    maxLength = MAX(maxLength, list->Length);
  }
  maxLength = MAX(maxLength, 2);

  print_sync_start(TRUE);
  printf("\n\n  Processor %d: Table of Variables at nodes:\n\n", ProcID);
  printf("Local  Global|   MatIDs |  Node   Num    |   VARIABLES:\n");
  printf("Node#  Node# |   Present|  Type  Unknowns|             \n ");
  fprint_line(stdout, "-", 80);
  for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
    for (i = 0; i < total_nodes; i++) {
      node_ptr = Nodes[i];
      nv = node_ptr->Nodal_Vars_Info[imtrx];
      printf("%-6d %-6d|", i, node_ptr->Global_Node_Num);
      list = &(node_ptr->Mat_List);
      for (j = 0; j < list->Length; j++) {
        printf("%3d ", list->List[j]);
      }
      for (j = list->Length; j < maxLength; j++) {
        printf("    ");
      }
      printf("|");
      if (i < Num_Internal_Nodes) {
        printf(" internal ");
      } else if (i < Num_Internal_Nodes + Num_Border_Nodes) {
        printf(" border   ");
      } else {
        printf(" ext(%3d) ", node_ptr->Proc);
      }
      printf(" %-5d  |", nv->Num_Unknowns);
      for (j = 0; j < nv->Num_Var_Desc; j++) {
        vd = nv->Var_Desc_List[j];
        printf(" %d", vd->Variable_Type);
        if (vd->Variable_Type == MASS_FRACTION) {
          printf("_sub%d", vd->Subvar_Index);
        }
        if (vd->MatID != -1) {
          printf("_mat%d", vd->MatID);
        }
        if (vd->Ndof > 1) {
          printf("*%d", vd->Ndof);
        }
      }
      printf("\n");
    }
  }
  fprint_line(stdout, "-", 80);
  print_sync_end(TRUE);
}
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void setup_external_nodal_vars(Exo_DB *exo, Dpi *dpi, Comm_Ex **cx)

/************************************************************************
 *
 * setup_external_nodal_vars():
 *
 *   This routine exchanges information about the solution vector
 * from each owned node to each ghost node. We first pack the information
 * about the solution vector on each owned node on this current
 * processor into a compact form. Then, we use the normal exchange
 * node information routine to exchange this information with the
 * neighboring processors.
 *   Then, we unpack the information obtained from neighboring
 * processors into nodal variable structures. And, then we use
 * the same procedure that we used in setup_local_nodal_vars() to
 * assign nodal var structures to external nodes.
 *   At the end of this procedure, we are assured that ghost nodes will
 * have the same picture of the solution vector as owned nodes.
 *
 ************************************************************************/
{
  int i, p, node_num, var_type;
  int imtrx;
  struct nv_packed **nvp = NULL, **nvp_send = NULL, **nvp_recv = NULL;
  NODAL_VARS_STRUCT *nv, *nv_match;
  COMM_NP_STRUCT *np_ptr, **np_base = NULL;
  /*
   *  Pack information for sending
   */

  if (dpi->num_neighbors > 0) {

    nvp_send = (struct nv_packed **)alloc_ptr_1(upd->Total_Num_Matrices);
    nvp_recv = (struct nv_packed **)alloc_ptr_1(upd->Total_Num_Matrices);
    nvp = (struct nv_packed **)alloc_ptr_1(upd->Total_Num_Matrices);

    for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
      nvp_send[imtrx] =
          alloc_struct_1(struct nv_packed, ptr_node_send[dpi->num_neighbors] * MaxVarPerNode);

      nvp_recv[imtrx] =
          alloc_struct_1(struct nv_packed, ptr_node_recv[dpi->num_neighbors] * MaxVarPerNode);
      //         nvp[imtrx] = alloc_struct_1(struct nv_packed,
      //			             ptr_node_recv[dpi->num_neighbors] * MaxVarPerNode);
    }

    np_base = (COMM_NP_STRUCT **)alloc_ptr_1(upd->Total_Num_Matrices);
    for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
      np_base[imtrx] = alloc_struct_1(COMM_NP_STRUCT, dpi->num_neighbors);
    }

    for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
      for (i = 0; i < ptr_node_send[dpi->num_neighbors]; i++) {
        node_num = list_node_send[i];
        nvp[imtrx] = nvp_send[imtrx] + i * MaxVarPerNode;
        pack_nv(nvp[imtrx], Nodes[node_num]->Nodal_Vars_Info[imtrx]);
      }

      np_ptr = np_base[imtrx];
      for (p = 0; p < dpi->num_neighbors; p++) {
        np_ptr->neighbor_ProcID = cx[imtrx][p].neighbor_name;
        np_ptr->send_message_buf = (void *)(nvp_send[imtrx] + ptr_node_send[p] * MaxVarPerNode);
        np_ptr->send_message_length =
            sizeof(struct nv_packed) * cx[imtrx][p].num_nodes_send * MaxVarPerNode;
        np_ptr->recv_message_buf = (void *)(nvp_recv[imtrx] + ptr_node_recv[p] * MaxVarPerNode);
        np_ptr->recv_message_length =
            sizeof(struct nv_packed) * cx[imtrx][p].num_nodes_recv * MaxVarPerNode;
        np_ptr++;
      }
      exchange_neighbor_proc_info(dpi->num_neighbors, np_base[imtrx]);

      for (i = 0; i < dpi->num_external_nodes; i++) {
        node_num = dpi->num_internal_nodes + dpi->num_boundary_nodes + i;
        nvp[imtrx] = nvp_recv[imtrx] + i * MaxVarPerNode;

        /*
         * Unpack the nodal variable structure into a new structure
         */
        nv = unpack_nv(nvp[imtrx], imtrx);

        /*
         * Decide whether this nodal variable structure is new.
         * Destroy it if it isn't. If it unique, it gets added to
         * a list of unique nodal variable structures.
         */
        if (MatchAdd_Node_Vars_List(nv, &nv_match)) {
          nodal_vars_destroy(&nv);
        }

        /*
         * Copy the pointer to the node vars structure into the node_info
         * structure for the current node.
         */
        Nodes[node_num]->Nodal_Vars_Info[imtrx] = nv_match;

        /*
         * Fill in other node structure information
         */
        Nodes[node_num]->Proc_Node_Num = node_num;

        /*
         *  Fill in the correct results for the number of degrees of freedom
         *  for external nodes into Dolphin[][]
         */
        for (var_type = V_FIRST; var_type < V_LAST; var_type++) {
          p = get_nv_ndofs_modMF(nv_match, var_type);
          if (p > 0) {
            Dolphin[imtrx][node_num][var_type] = p;
          }
        }
      }
    }
    /*
     *  Free memory allocated in this routine
     */

    if (dpi->num_neighbors > 0) {
      for (int i = 0; i < upd->Total_Num_Matrices; i++) {
        free(np_base[i]);
        //        free(nvp_save[i]);
        free(nvp_send[i]);
        free(nvp_recv[i]);
      }
    }
    safer_free((void **)&np_base);
    safer_free((void **)&nvp_send);
    safer_free((void **)&nvp_recv);
    safer_free((void **)&nvp);
  }
  /*
   *  When in debug mode, print out a complete listing of variables at
   *  every node
   */
}
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

int find_MaxUnknownNode(void)

/************************************************************************
 *
 * find_MaxUnknownNode():
 *
 *  Find the maximum number of unknowns at any one node.
 *
 * Return
 * --------
 *    Returns the maximum number of unknowns located at any one node.
 *
 ************************************************************************/
{
  int i, retn;
  int maxUnknowns = -1, gmax = 0;
  NODAL_VARS_STRUCT *nv;

  for (i = 0; i < Nodal_Vars_List_Length; i++) {
    nv = Nodal_Vars_List[i];
    if (nv->Num_Unknowns > maxUnknowns) {
      maxUnknowns = nv->Num_Unknowns;
    }
  }

#ifdef PARALLEL
  if (Num_Proc > 1) {
    retn = MPI_Allreduce(&maxUnknowns, &gmax, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (retn != MPI_SUCCESS) {
      fprintf(stderr, "find_MaxUnknownNode MPI error, P_%d, error = %d\n", ProcID, retn);
    }
  } else {
    gmax = maxUnknowns;
  }
#else
  gmax = maxUnknowns;
#endif

  return gmax;
}
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void set_unknown_map(Exo_DB *exo, Dpi *dpi)

/************************************************************************
 *
 * set_unknown_map():
 *
 *
 *
 ************************************************************************/
{
  int e, gcount, i, var_type, total_nodes, liver, count, i1, ii;
  int index, kspec, dofkspec, num;
  int imtrx;
  NODAL_VARS_STRUCT *nv;
  char position_label[MAX_DOFNAME];
  int gdof = 0;

  /*
   * Calculate the total number of nodes on this processor
   */
  total_nodes = Num_Internal_Nodes + Num_Border_Nodes + Num_External_Nodes;

  /*
   * Again loop through all of the nodes, this time figuring out how many
   * unknowns are really there, and how to set the Local_Offsets from the
   * First Unknown...
   */

  gcount = 0;
  for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
    for (i = 0; i < total_nodes; i++) {

      /*
       * Store the pointer to the nodal variables info at the
       * current node
       */
      nv = Nodes[i]->Nodal_Vars_Info[imtrx];

      count = 0;
      for (var_type = V_FIRST; var_type < V_LAST; var_type++) {
        if (Dolphin[imtrx][i][var_type] > 0) {
          /*
           * Notice how if you want any concentration unknowns you get
           * space for *all* concentration unknowns!
           *
           * Just consider it a bonus in the Sandia tradition...
           */
          if (var_type == MASS_FRACTION) {
            Local_Offset[imtrx][i][var_type] = count;
            liver = Dolphin[imtrx][i][var_type] * upd->Max_Num_Species_Eqn;
            count += liver;
          } else {
            Local_Offset[imtrx][i][var_type] = count;
            liver = Dolphin[imtrx][i][var_type];
            count += liver;
          }
        }
      }
      if (count != nv->Num_Unknowns) {
        fprintf(stderr, "P_%d: at node %d, Num_Unknowns_Node = %d, count = %d\n", ProcID, i,
                nv->Num_Unknowns, count);
        fprintf(stderr, "nv->Num_Unknowns lost count.\n");
        GOMA_EH(-1, "Num_Unknowns_Node incompatibility");
      }
    }
  }
  /*
   * Should verify that all nodes have been assigned something other
   * than "-1" as a number of unknowns...
   */
  for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
    for (i = 0; i < total_nodes; i++) {
      nv = Nodes[i]->Nodal_Vars_Info[imtrx];
      if (nv->Num_Unknowns == -1) {
        fprintf(stderr, "P_%d: local node index (%d) at matrix %d lost track.\n", ProcID, i + 1,
                imtrx + 1);
        GOMA_EH(-1, "Num_Unknowns_Node incompatibility");
      }
    }
  }

  /*
   * Determine total number of Unknowns to be solved for in this Proc
   * and the number that will be communicated from other Procs
   */

  num_internal_dofs = (int *)malloc(upd->Total_Num_Matrices * sizeof(int));
  num_boundary_dofs = (int *)malloc(upd->Total_Num_Matrices * sizeof(int));
  num_external_dofs = (int *)malloc(upd->Total_Num_Matrices * sizeof(int));
  num_universe_dofs = (int *)malloc(upd->Total_Num_Matrices * sizeof(int));

  /*
   * How many degrees of freedom I own that I do not have to communicate to
   * any other processor.
   */
  for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
    num_internal_dofs[imtrx] = 0;
    for (i = 0; i < Num_Internal_Nodes; i++) {
      nv = Nodes[i]->Nodal_Vars_Info[imtrx];
      num_internal_dofs[imtrx] += nv->Num_Unknowns;
    }
  }
  /*
   * How many degrees of freedom I own from which I will have to communicate
   * some to neighboring processors.
   */
  for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
    num_boundary_dofs[imtrx] = 0;
    for (i = Num_Internal_Nodes; i < Num_Internal_Nodes + Num_Border_Nodes; i++) {
      nv = Nodes[i]->Nodal_Vars_Info[imtrx];
      num_boundary_dofs[imtrx] += nv->Num_Unknowns;
    }
  }

  /*
   * How many degrees of freedom at nodes I do not own that I expect to receive
   * from other processors.
   */
  for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
    num_external_dofs[imtrx] = 0;
    for (i = Num_Internal_Nodes + Num_Border_Nodes;
         i < Num_Internal_Nodes + Num_Border_Nodes + Num_External_Nodes; i++) {
      nv = Nodes[i]->Nodal_Vars_Info[imtrx];
      num_external_dofs[imtrx] += nv->Num_Unknowns;
    }
  }

  /*
   * Count up the total number of degrees of freedom.
   */
  for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
    num_universe_dofs[imtrx] =
        num_internal_dofs[imtrx] + num_boundary_dofs[imtrx] + num_external_dofs[imtrx];
  }
  /*
   * Some element based quantities...
   */
  num_personal_elems = 0;

  if (Num_Proc == 1) {
    /*
     * In serial, they're all mine!
     */
    num_personal_elems = exo->num_elems;
  } else {

    /*
     * Not in parallel, though. We must check them all.
     */

    for (e = 0; e < exo->num_elems; e++) {
      if (dpi->elem_owner[e] == ProcID) {
        num_personal_elems++;
      }
    }
#if (DEBUG_LEVEL > 1)
    log_dbg("num_personal_elems = %d", num_personal_elems);
#endif
  }

  /*
   * Setup legacy variables...
   */

  NumUnknowns = (int *)malloc(upd->Total_Num_Matrices * sizeof(int));
  NumExtUnknowns = (int *)malloc(upd->Total_Num_Matrices * sizeof(int));

  for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
    NumUnknowns[imtrx] = num_internal_dofs[imtrx] + num_boundary_dofs[imtrx];
    NumExtUnknowns[imtrx] = num_external_dofs[imtrx];
  }

  /*
   * Produce the equation - variable interaction mask array (rf_masks.h)
   * This array is used durring fill process to determine the interation
   * between variables in a specific equation.
   *
   * This needs to be updated to use the pd->e[MATRIX_NUMBER][EQUATION_NAME] structure
   * instead of the many innumerable Booleans like "FluidMechanics"...
   */
  set_interaction_masks(exo);

  /*
   * Pointers into the solution vector...
   */

  /*
   * Define pointers into the solution vector on this Proc.
   * This solution vector is localy numbered so that all internal+
   * border variables (including pressure) come first.
   *
   * After all internal+border node unknowns
   * are numbered the external unknowns updated on other processors are
   * numbered. In this section of the vector each external nodes unknowns
   * are numbered (including pressure). To produce this vector the internal+
   * border nodes are first considered followed by the external nodes.
   */

  /*
   * Define pointer to first unknown, first mass fraction, and the pressure
   * unknown at each local node. This vector is first defined only for the
   * internal+border nodes. Also determine an upper bound on number of
   * unknowns at each node.
   */

  for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
    index = 0;
    for (i = 0; i < total_nodes; i++) {
      First_Unknown[imtrx][i] = index;
      Nodes[i]->First_Unknown[imtrx] = index;
      nv = Nodes[i]->Nodal_Vars_Info[imtrx];
      index += nv->Num_Unknowns;
    }
  }
  /*
   * Assign descriptive variable dof and residual equation names for
   * every degree of freedom and constraint residual equation.
   * The node numbers used in this description will be 1-based, in
   * compliance with BLOT, etc.
   */

  /*
   *   Only allocate dofname if doing numerical jacobian checking or writing
   *   intermediate solution or Debug_Flag > 1
   */
  if (dofname == NULL) {
    if (Debug_Flag < 0 || Debug_Flag > 1 || Write_Intermediate_Solutions || (Iout == 1)) {
      dofname = (char ***)malloc(upd->Total_Num_Matrices * sizeof(int **));
      for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
        dofname[imtrx] = alloc_VecFixedStrings(NumUnknowns[imtrx] + NumExtUnknowns[imtrx], 80);
      }
    }
  }

  /*
   *  Only allocate resname if doing numerical jacobian checking
   */
  if (resname == NULL) {
    if (Debug_Flag < 0 || Iout == 1) {
      resname = (char ***)malloc(upd->Total_Num_Matrices * sizeof(int **));
      for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
        resname[imtrx] = alloc_VecFixedStrings(NumUnknowns[imtrx] + NumExtUnknowns[imtrx], 80);
      }
    }
  }

  /*
   * Allocate only if unallocated...
   */
  if (idv == NULL) {
    idv = (int ***)malloc(upd->Total_Num_Matrices * sizeof(int **));
    for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
      idv[imtrx] = alloc_int_2(NumUnknowns[imtrx] + NumExtUnknowns[imtrx], 3, -1);
    }
  }

  /* Here this is a worst-case approach, used mainly for the post processing
   * We will have to revisit this if we want to go on with parallel processing
   */
  for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
    for (i = 0; i < total_nodes; i++) {
      nv = Nodes[i]->Nodal_Vars_Info[imtrx];
      gcount = First_Unknown[imtrx][i];

      if (Num_Proc > 1) {
        gdof = gcount;
        i1 = dpi->node_index_global[i] + 1; /* convert 0-based to 1-based */
      } else {
        gdof = gcount;
        i1 = i + 1; /* convert 0-based internal node numbering */
      }             /* to 1-based node numbering of BLOT, etc. */

      /*
       * Load up the coordinates of the node into a handy label
       * to help find errant dofs...
       */
      switch (Num_Dim) {
      case 1:
        sprintf(position_label, "x=%-10.6g", Coor[0][i]);
        break;
      case 2:
        sprintf(position_label, "x=%-10.6g y=%-10.6g", Coor[0][i], Coor[1][i]);
        break;
      case 3:
        sprintf(position_label, "x=%-9.6g y=%-9.6g z=%-9.6g", Coor[0][i], Coor[1][i], Coor[2][i]);
        break;
      }

      count = 0;

      for (var_type = V_FIRST; var_type < V_LAST; var_type++) {
        if (Dolphin[imtrx][i][var_type] > 0) {
          /*
           * Notice how if there are concentration unknowns at the
           * node, then we assign a  fixed number of
           * upd->Max_Num_Species_Eqn of them. If there are no
           * concentration unknowns, then Dolphin[matrix][i][MASS_FRACTION]
           * will be zero, so no space will be assigned for them
           */

          if (var_type == MASS_FRACTION) {
            num = upd->Max_Num_Species_Eqn;
          } else {
            num = 1;
          }
          liver = Dolphin[imtrx][i][var_type] * num;
          for (ii = 0; ii < liver; ii++) {
            idv[imtrx][gcount + count + ii][0] = var_type;
            idv[imtrx][gcount + count + ii][1] = ii;
            idv[imtrx][gcount + count + ii][2] = i;
            if (dofname != NULL) {
              if (var_type == MASS_FRACTION) {
                kspec = ii / Dolphin[imtrx][i][var_type];
                dofkspec = ii - kspec * Dolphin[imtrx][i][var_type];
                sprintf(dofname[imtrx][gcount + count + ii], "dof=%-5d %s%d_%d n=%05d %s",
                        gdof + count + ii, Var_Name[var_type].name2, kspec, dofkspec, i1,
                        position_label);
              } else {
                sprintf(dofname[imtrx][gcount + count + ii], "dof=%-5d %2s_%d n=%05d",
                        gdof + count + ii, Var_Name[var_type].name2, ii, i1);
              }
            }
            if (resname != NULL) {
              if (var_type == MASS_FRACTION) {
                kspec = ii / Dolphin[imtrx][i][var_type];
                dofkspec = ii - kspec * Dolphin[imtrx][i][var_type];
                sprintf(resname[imtrx][gcount + count + ii], "dof=%-5d %s%d_%d n=%05d %s",
                        gcount + count + ii, EQ_Name[var_type].name2, kspec, dofkspec, i1,
                        position_label);
              } else {
                sprintf(resname[imtrx][gcount + count + ii], "dof=%-5d %2s_%d n=%05d",
                        gcount + count + ii, EQ_Name[var_type].name2, ii, i1);
              }
            }
          } /* for: ii */
          count += liver;
        }
      }

      if (count != nv->Num_Unknowns) {
        fprintf(stderr, "P_%d: at node %d, Num_Unknowns_Node = %d, count = %d\n", ProcID, i,
                nv->Num_Unknowns, count);
        fprintf(stderr, "nv->Num_Unknowns lost count.");
        GOMA_EH(-1, "Sorry, Charlie.");
      }
    }
  }
  /*
   *    Debug_Flag > 2 output
   *      -> print out a mapping of what variables are active at
   *         each node.
   */
  if (Debug_Flag > 2) {
    print_vars_at_nodes();
  }

  /* Now that the dust has settled, let us translate these quantities
   * into Hoodian frontal solver form if the front method has been requested.
   * While you are at it, run prefront and also load up the element sweep map
   * regardless of solver.
   */
  return;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void set_interaction_masks(Exo_DB *exo)

/*
 *	Set interaction mask for the equation level interactions of
 *	the defined variables.
 *
 *	Author: 	John Shadid (SNL, 1421)
 *	Date: 		1/13/93
 *	Revised:	Tue Feb 15 07:20:14 MST 1994 pasacki@sandia.gov
 *
 *   Example:  Inter_Mask[matrix][equation][variable]
 *
 *   _______________________V A R I A B L E S  ( U N K N O W N S )____________
 *
 *                 U    U    U    T    Y    D    D    D    S    P
 *                  1    2    3         i    1    2    3    i
 *                 --   --   --   --   --   --   --   --   --   --
 * EQUATIONS
 * ---------
 *
 * MOM_1            1    1         1    1    1    1              1
 *
 * MOM_2            1    1         1    1    1    1              1
 *
 * MOM_3                      1              1    1
 *
 * ENERGY           1    1         1    1    1    1              1
 *
 * SPECIES_i        1    1         1    1    1    1              1
 *
 * MESH_1           1    1         1    1    1    1              1
 *
 * MESH_2           1    1         1    1    1    1              1
 *
 * MESH_3                                              1
 *
 * POLYMER_STRESS                                           1
 *
 * CONTINUITY       1    1         1    1    1    1              1
 *
 *
 *        U1      U2      U3      T       Y       S       P
 *      |-----|-------|-------|------|--------|-------|------|
 *   U1 |  1  |   1   |   1   |   1  |    1   |   1   |   1  |
 *      |-----|-------|-------|------|--------|-------|------|
 *   U2 |  1  |   1   |   1   |   1  |    1   |   1   |   1  |
 *      |-----|-------|-------|------|--------|-------|------|
 *   U3 |  1  |   1   |   1   |   1  |    1   |   1   |   1  |
 *      |-----|-------|-------|------|--------|-------|------|
 *   T  |  1  |   1   |   0   |   1  |    1   |   0   |   1  |
 *      |-----|-------|-------|------|--------|-------|------|
 *   Y  |  1  |   1   |   0   |   1  |    1   |   0   |   1  |
 *      |-----|-------|-------|------|--------|-------|------|
 *   S  |  1  |   1   |   1   |   1  |    1   |   1   |   1  |
 *      |-----|-------|-------|------|--------|-------|------|
 *   P  |  1  |   1   |   1   |   1  |    1   |   1   |   1  |
 *      |-----|-------|-------|------|--------|-------|------|
 *
 *	 1 - Implies this type of variable (row) interacts with this type of
 *	     variable (column)
 *	 0 - Implies no interaction
 *
 *	Variable Types:
 *		U1, U2, U3 - Three components of velocity
 *		T          - Temperature
 *		Y          - Mass fractions   (there can be more than one)
 *		D1, D2, D3 - Three components of mesh displacement...
 *		S          - Surface unknowns (there can be more than one)
 *		P          - Pressure
 *            S11,S12,
 *            S13,S22,   - Polymer Stress Tensor
 *            S23,S33
 *            G11,G12,G13,
 *            G21,G22,G23,- Velocity Gradient Tensor
 *            G31,G32,G33
 *            V           - Voltage potential
 *            F           - Fill
 *            PU1,PU2,PU3 - Three components of particle velocity
 *            P_LIQ, P_GAS, P_POR - Porous media variables
 */
{
  int ebi; /* element block indeces... */
  int e, i, n, v, mn;
  int imtrx;
  int eqn_var_mask[MAX_NUM_MATRICES][MAX_VARIABLE_TYPES][MAX_VARIABLE_TYPES];
  /* Default: Set the matrix Inter_Mask[][][] to all 0s*/

  for (imtrx = 0; imtrx < MAX_NUM_MATRICES; imtrx++) {
    for (n = 0; n < MAX_VARIABLE_TYPES; n++) {
      for (i = 0; i < MAX_VARIABLE_TYPES; i++) {
        Inter_Mask[imtrx][n][i] = 0;
        eqn_var_mask[imtrx][n][i] = 0;
      }
    }
  }
  /* try setting up Inter_Mask based on variables actually used
   *   in the equations
   */
  for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
    for (e = 0; e < MAX_EQNS; e++) {
      switch (e) {
      case R_MOMENTUM1:
      case R_MOMENTUM2:
      case R_MOMENTUM3:
        v = VELOCITY1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = TEMPERATURE;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MASS_FRACTION;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SOLID_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SOLID_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SOLID_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = PRESSURE;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = FILL;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = PHASE1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = PHASE2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = PHASE3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = PHASE4;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = PHASE5;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;

        v = VORT_DIR1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VORT_DIR2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VORT_DIR3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;

        v = POLYMER_STRESS11;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS12;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS13;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS22;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS23;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS33;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;

        v = POLYMER_STRESS11_1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS12_1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS13_1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS22_1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS23_1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS33_1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS11_2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS12_2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS13_2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS22_2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS23_2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS33_2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS11_3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS12_3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS13_3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS22_3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS23_3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS33_3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS11_4;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS12_4;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS13_4;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS22_4;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS23_4;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS33_4;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS11_5;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS12_5;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS13_5;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS22_5;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS23_5;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS33_5;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS11_6;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS12_6;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS13_6;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS22_6;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS23_6;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS33_6;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS11_7;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS12_7;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS13_7;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS22_7;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS23_7;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS33_7;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;

        v = VELOCITY_GRADIENT11;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT12;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT13;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT21;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT22;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT23;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT31;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT32;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT33;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;

        v = VOLTAGE;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POR_LIQ_PRES;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POR_GAS_PRES;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POR_POROSITY;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POR_TEMP;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POR_SATURATION;
        v = BOND_EVOLUTION;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;

        v = LAGR_MULT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = LAGR_MULT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = LAGR_MULT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;

        v = CURVATURE;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;

        v = EFIELD1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = EFIELD2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = EFIELD3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;

        v = NORMAL1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = NORMAL2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = NORMAL3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_SURF_DIV_V;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_SURF_CURV;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_NORMAL1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_NORMAL2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_NORMAL3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;

        v = N_DOT_CURL_V;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = GRAD_S_V_DOT_N1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = GRAD_S_V_DOT_N2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = GRAD_S_V_DOT_N3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_BDYVELO;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = LUBP;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_DELTAH;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_LUB_CURV;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_LUB_CURV_2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_FILMP;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_FILMH;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;

        v = MOMENT0;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MOMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MOMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MOMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = DENSITY_EQN;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;

        v = EDDY_NU;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = TURB_K;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = TURB_OMEGA;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;

        break;
      case R_USTAR:
      case R_VSTAR:
      case R_WSTAR:
        v = USTAR;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VSTAR;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = WSTAR;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        break;
      case R_PSTAR:
        v = PSTAR;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        break;
      case R_PMOMENTUM1:
      case R_PMOMENTUM2:
      case R_PMOMENTUM3:
        v = PVELOCITY1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = PVELOCITY2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = PVELOCITY3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHEAR_RATE;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = TEMPERATURE;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MASS_FRACTION;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SOLID_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SOLID_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SOLID_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = PRESSURE;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;

        break;

      case R_EXT_VELOCITY:
        v = EXT_VELOCITY;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MASS_FRACTION;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = FILL;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;

        break;

      case R_ENERGY:
        v = VELOCITY1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = TEMPERATURE;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MASS_FRACTION;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = FILL;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;

        v = VOLTAGE;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POR_LIQ_PRES;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POR_GAS_PRES;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POR_POROSITY;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POR_TEMP;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POR_SATURATION;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;

        v = MOMENT0;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MOMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MOMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MOMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = DENSITY_EQN;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        break;

      case R_MASS:
        v = VELOCITY1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = PVELOCITY1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = PVELOCITY2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = PVELOCITY3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = TEMPERATURE;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MASS_FRACTION;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = PRESSURE;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = FILL;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;

        v = VORT_DIR1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VORT_DIR2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VORT_DIR3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;

        v = SHEAR_RATE;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VOLTAGE;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POR_LIQ_PRES;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POR_GAS_PRES;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POR_POROSITY;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POR_TEMP;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POR_SATURATION;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;

        v = VELOCITY_GRADIENT11;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT12;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT13;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT21;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT22;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT23;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT31;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT32;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT33;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;

        v = LUBP;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;

        break;

      case R_MESH1:
      case R_MESH2:
      case R_MESH3:
        v = VELOCITY1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = PVELOCITY1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = PVELOCITY2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = PVELOCITY3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = TEMPERATURE;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MASS_FRACTION;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SOLID_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SOLID_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SOLID_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_CURVATURE;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_CURVATURE2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = PRESSURE;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POR_LIQ_PRES;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POR_GAS_PRES;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POR_POROSITY;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POR_TEMP;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POR_SATURATION;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POR_SINK_MASS;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_DIFF_FLUX;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT11;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT12;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT13;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT21;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT22;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT23;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT31;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT32;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT33;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS11;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS12;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS13;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS22;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS23;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS33;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;

        v = POLYMER_STRESS11_1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS12_1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS13_1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS22_1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS23_1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS33_1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS11_2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS12_2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS13_2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS22_2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS23_2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS33_2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS11_3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS12_3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS13_3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS22_3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS23_3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS33_3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS11_4;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS12_4;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS13_4;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS22_4;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS23_4;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS33_4;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS11_5;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS12_5;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS13_5;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS22_5;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS23_5;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS33_5;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS11_6;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS12_6;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS13_6;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS22_6;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS23_6;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS33_6;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS11_7;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS12_7;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS13_7;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS22_7;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS23_7;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS33_7;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_TENSION;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_USER;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_LUBP;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = ACOUS_PREAL;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = ACOUS_PIMAG;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = ACOUS_REYN_STRESS;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_BDYVELO;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = LUBP;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_DELTAH;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MAX_STRAIN;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = CUR_STRAIN;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_NORMAL1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_NORMAL2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_NORMAL3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        break;
        v = TFMP_PRES;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = TFMP_SAT;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        break;

      case R_MASS_SURF:
        v = VELOCITY1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = PVELOCITY1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = PVELOCITY2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = PVELOCITY3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = TEMPERATURE;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MASS_FRACTION;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = PRESSURE;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        break;

      case R_SOLID1:
      case R_SOLID2:
      case R_SOLID3:
        v = VELOCITY1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = TEMPERATURE;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MASS_FRACTION;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SOLID_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SOLID_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SOLID_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = PRESSURE;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_LUBP;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = ACOUS_PREAL;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = ACOUS_PIMAG;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = ACOUS_REYN_STRESS;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_BDYVELO;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = LUBP;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_DELTAH;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MAX_STRAIN;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = CUR_STRAIN;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        break;

      case R_PRESSURE:
        v = VELOCITY1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = TEMPERATURE;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MASS_FRACTION;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SOLID_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SOLID_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SOLID_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = PRESSURE;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS11;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS12;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS13;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS22;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS23;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS33;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;

        v = POLYMER_STRESS11_1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS12_1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS13_1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS22_1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS23_1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS33_1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS11_2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS12_2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS13_2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS22_2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS23_2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS33_2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS11_3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS12_3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS13_3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS22_3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS23_3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS33_3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS11_4;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS12_4;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS13_4;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS22_4;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS23_4;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS33_4;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS11_5;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS12_5;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS13_5;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS22_5;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS23_5;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS33_5;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS11_6;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS12_6;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS13_6;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS22_6;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS23_6;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS33_6;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS11_7;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS12_7;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS13_7;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS22_7;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS23_7;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS33_7;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;

        v = VELOCITY_GRADIENT11;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT12;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT13;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT21;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT22;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT23;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT31;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT32;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT33;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHEAR_RATE;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = FILL;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;

        v = MOMENT0;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MOMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MOMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MOMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = DENSITY_EQN;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = EDDY_NU;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = TURB_K;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = TURB_OMEGA;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        break;

      case R_STRESS11:
      case R_STRESS12:
      case R_STRESS13:
      case R_STRESS22:
      case R_STRESS23:
      case R_STRESS33:
        v = VELOCITY1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = TEMPERATURE;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MASS_FRACTION;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = PRESSURE;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS11;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS12;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS13;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS22;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS23;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS33;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;

        v = VELOCITY_GRADIENT11;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT12;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT13;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT21;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT22;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT23;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT31;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT32;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT33;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;

        v = SHEAR_RATE;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = FILL;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        break;

      case R_GRADIENT11:
      case R_GRADIENT12:
      case R_GRADIENT13:
      case R_GRADIENT21:
      case R_GRADIENT22:
      case R_GRADIENT23:
      case R_GRADIENT31:
      case R_GRADIENT32:
      case R_GRADIENT33:
        v = VELOCITY1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT11;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT12;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT13;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT21;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT22;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT23;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT31;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT32;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT33;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        break;

      case R_EFIELD1:
      case R_EFIELD2:
      case R_EFIELD3:
        v = R_EFIELD1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = R_EFIELD2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = R_EFIELD3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VOLTAGE;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        break;

      case R_POTENTIAL:
        v = VELOCITY1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = TEMPERATURE;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MASS_FRACTION;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = PRESSURE;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VOLTAGE;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SURF_CHARGE;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        break;

      case R_SURF_CHARGE:
        v = VELOCITY1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VOLTAGE;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SURF_CHARGE;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        break;

      case R_SHELL_CURVATURE:
        v = SHELL_TENSION;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_CURVATURE;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = PRESSURE;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;

        v = LUBP;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_FILMP;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;

        v = SHELL_NORMAL1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_NORMAL2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_NORMAL3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = TFMP_PRES;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        break;

      case R_SHELL_CURVATURE2:
        v = SHELL_CURVATURE2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;

        v = SHELL_NORMAL1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_NORMAL2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_NORMAL3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;

        break;

      case R_SHELL_ANGLE1:
        v = SHELL_ANGLE1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_ANGLE2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;

        break;

      case R_SHELL_ANGLE2:
        v = SHELL_ANGLE1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_ANGLE2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;

        break;

      case R_SHELL_TENSION:
        v = SHELL_TENSION;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_CURVATURE;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = PRESSURE;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = LUBP;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_FILMP;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        break;

      case R_SHELL_X:
        v = SHELL_CURVATURE;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_X;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_Y;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        break;

      case R_SHELL_Y:
        v = SHELL_CURVATURE;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_X;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_Y;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        break;

      case R_SHELL_USER:
        v = SHELL_USER;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = TEMPERATURE;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MASS_FRACTION;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SOLID_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SOLID_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SOLID_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = PRESSURE;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VOLTAGE;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        break;

      case R_SHELL_SURF_DIV_V:
        v = SHELL_SURF_DIV_V;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_NORMAL1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_NORMAL2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_SURF_CURV;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;

        break;

      case R_SHELL_SURF_CURV:
        v = SHELL_SURF_CURV;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_NORMAL1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_NORMAL2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        break;

      case R_N_DOT_CURL_V:
        v = N_DOT_CURL_V;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        break;

      case R_GRAD_S_V_DOT_N1:
      case R_GRAD_S_V_DOT_N2:
      case R_GRAD_S_V_DOT_N3:
        v = GRAD_S_V_DOT_N1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = GRAD_S_V_DOT_N2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = GRAD_S_V_DOT_N3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_NORMAL1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_NORMAL2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        break;

        /*
         * EDW Note: For the following 4 shell equation entries,
         * all three mesh components are included, even though
         * the equations are not yet 3D-compatible!
         */
      case R_SHELL_DIFF_FLUX:
        v = SHELL_DIFF_FLUX;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_DIFF_CURVATURE;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        break;

      case R_SHELL_DIFF_CURVATURE:
        v = SHELL_DIFF_CURVATURE;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_NORMAL1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_NORMAL2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        break;

      case R_SHELL_NORMAL1:
        v = SHELL_NORMAL1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        break;

      case R_SHELL_NORMAL2:
        v = SHELL_NORMAL2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        break;

      case R_SHELL_NORMAL3:
        v = SHELL_NORMAL3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        break;

      case R_EDDY_NU:
        v = EDDY_NU;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;

        v = MESH_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;

        v = VELOCITY1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;

        break;
      case R_TURB_K:
        v = TURB_K;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = TURB_OMEGA;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;

        v = MESH_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;

        v = VELOCITY1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;

        break;
      case R_TURB_OMEGA:
        v = TURB_OMEGA;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = TURB_K;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;

        v = MESH_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;

        v = VELOCITY1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;

        break;

      case R_ACOUS_PREAL:
      case R_ACOUS_PIMAG:
      case R_ACOUS_REYN_STRESS:
      case R_SHELL_BDYVELO:
        v = ACOUS_PREAL;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = ACOUS_PIMAG;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = ACOUS_REYN_STRESS;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_BDYVELO;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MASS_FRACTION;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = FILL;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        break;

      case R_LIGHT_INTP:
      case R_LIGHT_INTM:
      case R_LIGHT_INTD:
        v = LIGHT_INTP;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = LIGHT_INTM;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = LIGHT_INTD;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = FILL;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        break;

      case R_RESTIME:
        v = RESTIME;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = TEMPERATURE;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MASS_FRACTION;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = FILL;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        break;

      case R_SHELL_LUBP:
        v = SHELL_LUBP;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_BDYVELO;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MASS_FRACTION;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = FILL;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        break;

      case R_LUBP:
        v = LUBP;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = FILL;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_SAT_CLOSED;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_PRESS_OPEN;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_TEMPERATURE;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_DELTAH;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_LUB_CURV;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_SAT_GASN;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_FILMP;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_PARTC;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SOLID_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SOLID_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SOLID_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_SHEAR_TOP;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_SHEAR_BOT;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_CROSS_SHEAR;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;

        v = SHELL_NORMAL1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_NORMAL2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_NORMAL3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;

        v = VELOCITY1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = PRESSURE;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MASS_FRACTION;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;

        /* Need to add the height-var here and velocity var here */
        break;

      case R_LUBP_2:
        v = LUBP_2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = FILL;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = PHASE1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_SAT_CLOSED;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_PRESS_OPEN;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_PRESS_OPEN_2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_TEMPERATURE;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_DELTAH;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_LUB_CURV_2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_SAT_GASN;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_FILMP;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_PARTC;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SOLID_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SOLID_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SOLID_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_SHEAR_TOP;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_SHEAR_BOT;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_CROSS_SHEAR;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;

        v = SHELL_NORMAL1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_NORMAL2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_NORMAL3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;

        /* Need to add the height-var here and velocity var here */
        break;

      case R_SHELL_FILMP:
        v = SHELL_FILMP;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_FILMH;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_PARTC;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = LUBP;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = PRESSURE;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_NORMAL1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_NORMAL2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_NORMAL3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;

        break;

      case R_SHELL_FILMH:
        v = SHELL_FILMH;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_FILMP;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_PARTC;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_NORMAL1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_NORMAL2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_NORMAL3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        break;

      case R_SHELL_PARTC:
        v = SHELL_PARTC;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_FILMP;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_FILMH;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = LUBP;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        break;

      case R_SHELL_SAT_CLOSED:
        v = SHELL_SAT_CLOSED;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = LUBP;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = FILL;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_LUB_CURV;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_LUB_CURV_2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_SAT_GASN;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        break;

      case R_SHELL_SAT_OPEN:
        v = SHELL_PRESS_OPEN;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = LUBP;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = LUBP_2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = FILL;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = PHASE1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_PRESS_OPEN_2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POR_SINK_MASS;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        break;

      case R_SHELL_SAT_OPEN_2:
        v = PHASE1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = LUBP_2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = FILL;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_PRESS_OPEN_2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        break;

      case R_SHELL_SAT_1:
        v = SHELL_SAT_1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_SAT_2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_SAT_3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POR_SINK_MASS;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        break;

      case R_SHELL_SAT_2:
        v = SHELL_SAT_1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_SAT_2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_SAT_3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        break;

      case R_SHELL_SAT_3:
        v = SHELL_SAT_1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_SAT_2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_SAT_3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        break;

      case R_SHELL_ENERGY:
        v = SHELL_TEMPERATURE;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_FILMH;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = LUBP;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = FILL;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_DELTAH;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SOLID_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SOLID_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SOLID_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;

        break;

      case R_SHELL_DELTAH:
        v = SHELL_DELTAH;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_TEMPERATURE;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = LUBP;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SOLID_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SOLID_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SOLID_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        break;

      case R_SHELL_LUB_CURV:
        v = SHELL_LUB_CURV;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = FILL;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        break;

      case R_SHELL_LUB_CURV_2:
        v = SHELL_LUB_CURV_2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = PHASE1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        break;

      case R_SHELL_SAT_GASN:
        v = SHELL_SAT_GASN;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_SAT_CLOSED;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = FILL;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        break;

      case R_POR_SINK_MASS:
        v = POR_SINK_MASS;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POR_POROSITY;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POR_LIQ_PRES;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POR_GAS_PRES;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_PRESS_OPEN;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_SAT_1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;

        break;

      case R_SHELL_SHEAR_TOP:
        v = LUBP;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_SHEAR_TOP;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_SHEAR_BOT;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_CROSS_SHEAR;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;

        break;

      case R_SHELL_SHEAR_BOT:
        v = LUBP;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_SHEAR_TOP;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_SHEAR_BOT;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_CROSS_SHEAR;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;

        break;

      case R_SHELL_CROSS_SHEAR:
        v = LUBP;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_SHEAR_TOP;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_SHEAR_BOT;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;

        break;

      case R_MAX_STRAIN:
        v = MAX_STRAIN;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;

        break;

      case R_CUR_STRAIN:
        v = CUR_STRAIN;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;

        break;

      case R_FILL:
        v = FILL;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SOLID_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SOLID_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SOLID_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = EXT_VELOCITY;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        /* for eikonal methods, R_FILL can depend on anything the
           various kinematic conditions can depend on
         */
        v = MASS_FRACTION;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = LUBP;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_LUB_CURV;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_DELTAH;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_NORMAL1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_NORMAL2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_NORMAL3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        break;

      case R_SHEAR_RATE:
        v = VELOCITY1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHEAR_RATE;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        break;

      case R_ENORM:
        v = VOLTAGE;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = ENORM;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        break;

      case R_BOND_EVOLUTION:
        v = VELOCITY1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHEAR_RATE;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = BOND_EVOLUTION;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        break;

      case R_CURVATURE:
        v = CURVATURE;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = FILL;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = NORMAL1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = NORMAL2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = NORMAL3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        break;

      case R_PHASE1:
      case R_PHASE2:
      case R_PHASE3:
      case R_PHASE4:
      case R_PHASE5:
        v = VELOCITY1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;

        v = MESH_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;

        v = PHASE1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = PHASE2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = PHASE3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = PHASE4;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = PHASE5;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;

        v = LUBP_2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;

        v = SHELL_LUB_CURV_2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        break;

      case R_NORMAL1:
      case R_NORMAL2:
      case R_NORMAL3:
        v = (e == R_NORMAL1 ? NORMAL1
                            : (e == R_NORMAL2 ? NORMAL2 : (e == R_NORMAL3 ? NORMAL3 : -1)));

        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = FILL;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        break;

      case R_VORT_DIR1:
      case R_VORT_DIR2:
      case R_VORT_DIR3:
      case R_VORT_LAMBDA:
        v = VELOCITY1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VORT_DIR1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VORT_DIR2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VORT_DIR3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHEAR_RATE;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MASS_FRACTION;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = PRESSURE;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        break;

      case R_LAGR_MULT1:
      case R_LAGR_MULT2:
      case R_LAGR_MULT3:
        v = VELOCITY1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = TEMPERATURE;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MASS_FRACTION;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SOLID_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SOLID_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SOLID_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = PRESSURE;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VOLTAGE;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POR_LIQ_PRES;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POR_GAS_PRES;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POR_POROSITY;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POR_TEMP;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POR_SATURATION;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = FILL;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHEAR_RATE;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = LAGR_MULT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = LAGR_MULT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = LAGR_MULT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        break;

      case R_STRESS11_1:
      case R_STRESS12_1:
      case R_STRESS13_1:
      case R_STRESS22_1:
      case R_STRESS23_1:
      case R_STRESS33_1:
        v = VELOCITY1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS11_1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS12_1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS13_1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS22_1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS23_1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS33_1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = PRESSURE;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT11;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT12;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT13;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT21;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT22;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT23;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT31;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT32;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT33;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = FILL;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        break;

      case R_STRESS11_2:
      case R_STRESS12_2:
      case R_STRESS13_2:
      case R_STRESS22_2:
      case R_STRESS23_2:
      case R_STRESS33_2:
        v = VELOCITY1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;

        v = POLYMER_STRESS11_2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS12_2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS13_2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS22_2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS23_2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS33_2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = PRESSURE;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT11;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT12;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT13;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT21;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT22;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT23;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT31;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT32;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT33;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = FILL;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        break;

      case R_STRESS11_3:
      case R_STRESS12_3:
      case R_STRESS13_3:
      case R_STRESS22_3:
      case R_STRESS23_3:
      case R_STRESS33_3:
        v = VELOCITY1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS11_3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS12_3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS13_3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS22_3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS23_3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS33_3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = PRESSURE;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT11;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT12;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT13;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT21;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT22;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT23;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT31;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT32;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT33;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = FILL;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        break;

      case R_STRESS11_4:
      case R_STRESS12_4:
      case R_STRESS13_4:
      case R_STRESS22_4:
      case R_STRESS23_4:
      case R_STRESS33_4:
        v = VELOCITY1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS11_4;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS12_4;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS13_4;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS22_4;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS23_4;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS33_4;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = PRESSURE;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT11;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT12;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT13;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT21;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT22;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT23;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT31;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT32;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT33;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = FILL;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        break;

      case R_STRESS11_5:
      case R_STRESS12_5:
      case R_STRESS13_5:
      case R_STRESS22_5:
      case R_STRESS23_5:
      case R_STRESS33_5:
        v = VELOCITY1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS11_5;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS12_5;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS13_5;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS22_5;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS23_5;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS33_5;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = PRESSURE;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT11;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT12;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT13;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT21;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT22;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT23;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT31;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT32;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT33;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = FILL;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        break;

      case R_STRESS11_6:
      case R_STRESS12_6:
      case R_STRESS13_6:
      case R_STRESS22_6:
      case R_STRESS23_6:
      case R_STRESS33_6:
        v = VELOCITY1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS11_6;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS12_6;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS13_6;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS22_6;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS23_6;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS33_6;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = PRESSURE;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT11;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT12;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT13;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT21;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT22;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT23;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT31;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT32;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT33;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = FILL;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        break;

      case R_STRESS11_7:
      case R_STRESS12_7:
      case R_STRESS13_7:
      case R_STRESS22_7:
      case R_STRESS23_7:
      case R_STRESS33_7:
        v = VELOCITY1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS11_7;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS12_7;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS13_7;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS22_7;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS23_7;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POLYMER_STRESS33_7;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = PRESSURE;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT11;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT12;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT13;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT21;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT22;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT23;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT31;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT32;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY_GRADIENT33;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = FILL;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        break;

      case R_POR_LIQ_PRES:
      case R_POR_GAS_PRES:
      case R_POR_POROSITY:
      case R_POR_SATURATION:
        v = VELOCITY1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = TEMPERATURE;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = PRESSURE;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MASS_FRACTION;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POR_LIQ_PRES;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POR_GAS_PRES;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POR_POROSITY;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POR_TEMP;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POR_SATURATION;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POR_SINK_MASS;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = LS;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        break;

      case R_MOMENT0:
      case R_MOMENT1:
      case R_MOMENT2:
      case R_MOMENT3:
        v = VELOCITY1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MOMENT0;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MOMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MOMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MOMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = TEMPERATURE;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MASS_FRACTION;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        break;
      case R_DENSITY_EQN:
        v = DENSITY_EQN;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MOMENT0;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MOMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MOMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MOMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = TEMPERATURE;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MASS_FRACTION;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        break;

      case R_TFMP_MASS:
      case R_TFMP_BOUND:
        v = TFMP_PRES;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = TFMP_SAT;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_NORMAL1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_NORMAL2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_NORMAL3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_CURVATURE;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = SHELL_TENSION;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        break;

      case R_SPECIES_UNK_0:
      case R_SPECIES_UNK_1:
      case R_SPECIES_UNK_2:
      case R_SPECIES_UNK_3:
      case R_SPECIES_UNK_4:
      case R_SPECIES_UNK_5:
      case R_SPECIES_UNK_6:
      case R_SPECIES_UNK_7:
      case R_SPECIES_UNK_8:
      case R_SPECIES_UNK_9:
      case R_SPECIES_UNK_10:
      case R_SPECIES_UNK_11:
      case R_SPECIES_UNK_12:
      case R_SPECIES_UNK_13:
      case R_SPECIES_UNK_14:
      case R_SPECIES_UNK_15:
      case R_SPECIES_UNK_16:
      case R_SPECIES_UNK_17:
      case R_SPECIES_UNK_18:
      case R_SPECIES_UNK_19:
      case R_SPECIES_UNK_20:
      case R_SPECIES_UNK_21:
      case R_SPECIES_UNK_22:
      case R_SPECIES_UNK_23:
      case R_SPECIES_UNK_24:
      case R_SPECIES_UNK_25:
      case R_SPECIES_UNK_26:
      case R_SPECIES_UNK_27:
      case R_SPECIES_UNK_28:
      case R_SPECIES_UNK_29:

        v = VELOCITY1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = VELOCITY3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT1;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT2;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = MESH_DISPLACEMENT3;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POR_LIQ_PRES;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POR_GAS_PRES;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POR_POROSITY;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POR_TEMP;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = POR_SATURATION;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        break;

      case R_EM_E1_REAL:
      case R_EM_E2_REAL:
      case R_EM_E3_REAL:
      case R_EM_E1_IMAG:
      case R_EM_E2_IMAG:
      case R_EM_E3_IMAG:
      case R_EM_H1_REAL:
      case R_EM_H2_REAL:
      case R_EM_H3_REAL:
      case R_EM_H1_IMAG:
      case R_EM_H2_IMAG:
      case R_EM_H3_IMAG:
        v = EM_E1_REAL;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = EM_E2_REAL;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = EM_E3_REAL;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = EM_E1_IMAG;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = EM_E2_IMAG;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = EM_E3_IMAG;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = EM_H1_REAL;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = EM_H2_REAL;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = EM_H3_REAL;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = EM_H1_IMAG;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = EM_H2_IMAG;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = EM_H3_IMAG;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = EM_CONT_REAL;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = EM_CONT_IMAG;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        break;
      case R_EM_CONT_REAL:
      case R_EM_CONT_IMAG:
        v = EM_E1_REAL;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = EM_E2_REAL;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = EM_E3_REAL;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = EM_E1_IMAG;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = EM_E2_IMAG;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        v = EM_E3_IMAG;
        if (Num_Var_In_Type[imtrx][v])
          eqn_var_mask[imtrx][e][v] = 1;
        break;
      }
    }

  } /* End of loop over matrices */

  /* DRN: This attempts to shrink the matrix if we are using
     explicit fill function and are choosing to ignore the
     dependencies on fill
   */
  if (ls != NULL && ls->Ignore_F_deps) {
    for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
      for (e = 0; e < MAX_EQNS; e++) {
        for (v = 0; v < MAX_VARIABLE_TYPES; v++) {
          if ((e == R_FILL && v != FILL) || (e != R_FILL && v == FILL))
            eqn_var_mask[imtrx][e][v] = 0;
        }
      }
    }
  }

  for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
    for (e = 0; e < MAX_EQNS; e++) {
      for (v = 0; v < MAX_VARIABLE_TYPES; v++) {
        if (Ignore_Deps[imtrx][e][v])
          eqn_var_mask[imtrx][e][v] = 0;
      }
    }
  }

  /* Here we are simply looking for the worst-case scenario, so
     that the boundaries are handled correctly.  The individual
     pieces of the domain will not have extraneous memory allocated
     at them for variables not defined at the nodes.  This is
     automatically taken care of, I believe.  */

  for (ebi = 0; ebi < exo->num_elem_blocks; ebi++) {
    mn = Matilda[ebi];
    if (mn < 0) {
      continue;
    }
    for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
      for (e = 0; e < MAX_EQNS; e++) {
        for (v = 0; v < MAX_VARIABLE_TYPES; v++) {
          Inter_Mask[imtrx][e][v] |= (pd_glob[mn]->e[imtrx][e] & eqn_var_mask[imtrx][e][v]) ? 1 : 0;
        }
      }
    }
  }

  if (Debug_Flag > 2) {
    printf("\n\nProc: %d\n", ProcID);
    for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
      for (n = 0; n < MAX_VARIABLE_TYPES; n++) {
        printf("\n");
        for (i = 0; i < MAX_VARIABLE_TYPES; i++) {
          printf(" %d ", Inter_Mask[imtrx][n][i]);
        }
      }
    }
    printf("\n\n");
  }
}
/**************************************************************************/
/**************************************************************************/
/**************************************************************************/

int Index_Solution(const int nodeNum,
                   const int varType,
                   const int subvarIndex,
                   const int iNdof,
                   const int matID,
                   const int imtrx)

/********************************************************************
 *
 * Index_Solution():
 *
 * Routine to return the index in the solution vector for the
 * iNunk'th variable of type iVarType at the node iGlobNum
 * corresponding to the material, matID.
 * Return -1 if undefined.
 *
 * Input
 * --------
 *  nodeNum     = processor Node Number
 *  varType     = Variable Type
 *  subvarIndex = Subvariable index -> only for MASS_FRACTION
 *                variable types is this active. Will go
 *                away.
 *  iNdof       = Local nodal degree of freedom. This is
 *                equal to zero, except for centroid
 *                pressures, and hermite cubic interpolated
 *                variables.
 *                    (HKM -> see backwards compatibility section
 *                            below)
 *  matID       = material id -> or -1 if the specification
 *                of the material doesn't matter in this
 *                instance. It would matter if the soln
 *                number at this node varied depending on
 *                the material id.  That is, if this variable has
 *                potentially more than one value at a given node.
 *              = -2: For this special input value, a match
 *                    will be returned for the first variable
 *                    description structure found that
 *                    matches the varType, subvartype pair.
 *
 * imtrx        = Matrix ID, goes from 0 < imtrx < upd->Total_Num_Matrices
 *
 * Output
 * --------
 * This function returns the processor index into the solution
 * vector for this degree of freedom.
 * If this degree of freedom doesn't exist, this function returns
 * a value of -1.
 *
 * NOTES:
 *   This function should execute as fast as possible, because
 *   it is unfortunately called during low levels of the jacobian
 *   and residual fills.
 **********************************************************************/
{
  int dofp, index, i_match = -1, i, ifound;
  /*
   * Pointer to the node info struct for this node
   */
  NODE_INFO_STRUCT *node_ptr = Nodes[nodeNum];
  NODAL_VARS_STRUCT *nv = node_ptr->Nodal_Vars_Info[imtrx];
  VARIABLE_DESCRIPTION_STRUCT *vd_match = NULL, *vd;
  /*
   * Do extra debugging of argument list when in debug mode
   */

  /*
   * Quick return for var types not present at this node
   */
  if (nv->Num_Var_Desc_Per_Type[varType] == 0)
    return -1;

  if (varType == MASS_FRACTION) {
    ifound = 0;
    if (subvarIndex >= upd->Max_Num_Species_Eqn) {
      printf("ERROR Index_Solution: subvarIndex is bad: %d\n", subvarIndex);
      GOMA_EH(GOMA_ERROR, "ERROR Index_Solution: subvarIndex is bad");
    }
    index = nv->Num_Var_Desc_Per_Type[varType] / upd->Max_Num_Species_Eqn;
    for (i = 0; i < nv->Num_Var_Desc; i++) {
      vd = nv->Var_Desc_List[i];
      if ((vd->Variable_Type == MASS_FRACTION) && (vd->Subvar_Index == subvarIndex)) {
        if (vd->MatID == matID || matID == -2) {
          vd_match = vd;
          i_match = i;
          break;
        } else if (vd->MatID == -1) {
          vd_match = vd;
          i_match = i;
        }
        ifound++;
        if (ifound == index)
          break;
      }
    }
  } else {
    ifound = 0;
    for (i = 0; i < nv->Num_Var_Desc; i++) {
      vd = nv->Var_Desc_List[i];
      if (vd->Variable_Type == varType) {
        if (vd->MatID == matID || matID == -2) {
          vd_match = vd;
          i_match = i;
          break;
        } else if (vd->MatID == -1) {
          vd_match = vd;
          i_match = i;
        }
        ifound++;
        if (ifound == (int)nv->Num_Var_Desc_Per_Type[varType])
          break;
      }
    }
  }

  if (vd_match == NULL) {
    return -1;
  }

  /*
   * Gather the local number of degrees of this variable type at this
   * node into a local variable, dofp.
   */
  dofp = vd_match->Ndof;

  /*
   * Check to see whether this variable type is defined to exist at
   * this global node number. If it doesn't exist, return immediately.
   */
  if (dofp < 1)
    return (-1);

  /*
   *  The index is determined by several arrays:
   *    First_Unknown[matrix number][nodeNum] -> index into the processor node list for
   *              the first unknown pertaining to the node number
   *    Local_offset[][][] -> offset for the particular processor node
   *              var type from the First_Unknown.
   *        -> note this can be redefined to be a short int array
   *           if space is needed.
   *    iNdof -> offset from Local_Offset! -> they are necessary
   *             defined as continguous within the solution vector.
   *
   * NOTE -> Local_Offset will break under the new treatment with
   *         mat_id. Will have to do away with it, or surplant it
   *         for cases where discontinuous vars occur.
   */

  /*
   * HKM -> Special backwards compatibility section that will go away
   *
   *  Under the old method, iNdof was used as an index into the solution
   *  vector for the case of discontinuous variables. Below we will protect
   *  against this functionality grandfathered into Index_Solution()'s
   *  API. It's actually already taken into account in the Nodal_Offsets
   *  in the new method. However, it's needed by the old method to get
   *  the right index and we check old vs. new in the DEBUG_HKM section
   *  below. Therefore, we will add a section of code to discard iNdof
   *  for almost all variable types, except pressure.
   *
   */
  if (vd_match->Ndof > 1) {
    index = node_ptr->First_Unknown[imtrx] + nv->Nodal_Offset[i_match] + iNdof;
  } else {
    index = node_ptr->First_Unknown[imtrx] + nv->Nodal_Offset[i_match];
  }

  return index;
}
/***************************************************************************/
/***************************************************************************/
/***************************************************************************/

int variable_type_nodalInterp(int var)

/***********************************************************************
 *
 * variable_type_nodalInterp()
 *
 *   This code checks loops over all materials checking to see if a
 *   variable type exists in the problem, and is either part of the
 *   solution vector or is part of an interpolation using the finite
 *   element basis function structure.  It does this (currently) by
 *   checking to see if the variable has a non-zero interpolation defined
 *   in any of the Problem_Description structures pertaining to different
 *   materials.
 *   I don't know if this is the most straightforward way to accomplish
 *   this goal. However, it is what is used currently.
 *   The constants pertaining to the types of interpolations
 *   are located in mm_as_const.h.
 *
 *   Only certain types of variable interpolations are considered, those
 *   that basically consist of continuous interpolation using one
 *   value at each node.
 *
 *   The vector, Num_Var_In_Type[][], is checked for non-null too. Again,
 *   I am not sure why both checks are needed.
 *
 *   Implementation note -> I do not know why just some of the
 *      interpolations are included in the list and not other. For
 *      example, why aren't hermite cubics in the list below?
 *
 *   Input
 *   -------
 *    var => Variable type (defined in rf_fem_const.h)
 *
 *   Return
 *   ---------
 *    TRUE -> unknowns corresponding to Variable types exist somewhere
 *            in the domain.
 *    FALSE -> Variable type doesn't exist in the problem.
 **********************************************************************/
{
  int i, post_flag, i_type;
  struct Problem_Description *pd_local;
  post_flag = 0;
  if (Num_Var_In_Type[pg->imtrx][var]) {
    for (i = 0; i < upd->Num_Mat; i++) {
      pd_local = pd_glob[i];
      i_type = pd_local->i[pg->imtrx][var];
      if (i_type == I_Q1 || i_type == I_Q2 || i_type == I_Q1_G || i_type == I_Q2_G ||
          i_type == I_Q1_GP || i_type == I_Q2_GP || i_type == I_Q1_GN || i_type == I_Q2_GN ||
          i_type == I_Q1_XV || i_type == I_Q2_XV || i_type == I_Q1_XG || i_type == I_Q2_XG ||
          i_type == I_Q1_HV || i_type == I_Q1_HG || i_type == I_Q1_HVG || i_type == I_Q2_HV ||
          i_type == I_Q2_HG || i_type == I_Q2_HVG || i_type == I_Q2_D || i_type == I_Q2_LSA ||
          i_type == I_Q2_D_LSA || i_type == I_Q1_D || i_type == I_H3 || i_type == I_S2 ||
          i_type == I_SP || i_type == I_N1)
        post_flag = 1;
    }
  }
  return post_flag;
}
/***************************************************************************/
/***************************************************************************/
/***************************************************************************/

VARIABLE_DESCRIPTION_STRUCT *Index_Solution_Inv(
    const int gindex, int *inode, int *i_Var_Desc, int *i_offset, int *idof, int imtrx)

/************************************************************************
 *
 * Index_Solution_Inv():
 *
 * Routine to return the node number and variable description for
 * a given processor unknown index.
 *
 *   Input
 *       gindex =  Processor value of the unknown number, i.e., the
 *                 index into the processor solution vector.
 *       *inode =  User supplied guess as to the node number pertaining
 *                 to gindex. Setting it to zero when in doubt
 *                 is a good choice.
 *       imtrx  =  Matrix ID: 0 <= imtrx < upd->Total_Num_Matrices
 *   Output
 *      *inode = local node number
 *      *i_offset = offset from the first degree of freedom at the node
 *                  in the solution vector
 *      vd = pointer to the corresponding variable description struct
 *      *i_Var_Desc = The index of this variable description structure
 *                    in the list of variable description structures
 *                    for this node.
 *      *idof  = ID of the degree of freedom if there are more than
 *               one degree of freedom contained in this variable
 *               description structure (usually there is not).
 *   Return
 *    nonNULL = pointer to the corresponding variable description
 *              struct
 *     NULL   = Global degree of freedom doesn't exist
 *
 *   NOTES:
 *      Any or all of the output variables, except for vd, may be set
 *      to NULL on input. In this case, nothing is returned in that
 *      particular location.
 ************************************************************************/
{
  int val_mid, val_top, diff, middle, top, bottom;
  int ivd, Inode, I_offset;
  NODAL_VARS_STRUCT *nv;
  VARIABLE_DESCRIPTION_STRUCT *vd;
  if (gindex < 0)
    return NULL;

  /*
   * Do a modified binary search to find the global node number
   * This assumes that First_Unknown is a unique, monotonically,
   * increasing function of the processor node number.
   *  Note: log(N) cost so this is well worth it.
   */
  top = DPI_ptr->num_owned_nodes + DPI_ptr->num_external_nodes - 1;
  val_top = Nodes[top]->First_Unknown[imtrx];
  if (gindex >= val_top) {
    Inode = top;
    goto found_node;
  }
  /*
   * First check to see if user supplied a good guess for inode,
   * based on a previous call to this routine, ie., we are
   * monotomically increasing gindex as we call this routine.
   */

  if (inode) {
    middle = (int)(*inode);
    if (middle >= 0 && middle < top) {
      if (gindex < Nodes[middle + 1]->First_Unknown[imtrx]) {
        if (gindex >= Nodes[middle]->First_Unknown[imtrx]) {
          Inode = middle;
          goto found_node;
        }
      } else {
        if (gindex < Nodes[middle + 2]->First_Unknown[imtrx]) {
          Inode = middle + 1;
          goto found_node;
        }
      }
    }
  }
  /*
   * Ok, that didn't work. Do an arbitrary binary search
   */
  bottom = 0;
  while ((diff = top - bottom) > 1) {
    middle = bottom + diff / 2;
    val_mid = Nodes[middle]->First_Unknown[imtrx];
    if (gindex > val_mid) {
      bottom = middle;
    } else if (gindex < val_mid) {
      top = middle;
    } else {
      Inode = middle;
      goto found_node;
    }
  }
  Inode = bottom;

found_node:;

  /*
   * Handle the problematic case of zero unknowns at the
   * node just found. The way First_Unknown field is
   * set up currently when the node has zero unknowns
   * associated with it is that it is equal to the value
   * of the next valid First_Unknown field counting
   * upwards. Therefore, we just increase Inode by 1
   * until we hit a node with unknowns.
   */

  nv = Nodes[Inode]->Nodal_Vars_Info[imtrx];
  while (nv->Num_Unknowns == 0) {
    Inode++;
    nv = Nodes[Inode]->Nodal_Vars_Info[imtrx];
  }
  /*
   * OK, we found the global node number, Inode.
   */
  if (inode)
    *inode = Inode;
  I_offset = gindex - Nodes[Inode]->First_Unknown[imtrx];
  if (i_offset)
    *i_offset = I_offset;
  bottom = 0;
  for (ivd = 0; ivd < nv->Num_Var_Desc; ivd++) {
    vd = nv->Var_Desc_List[ivd];
    if (I_offset < (bottom + vd->Ndof)) {
      if (i_Var_Desc)
        *i_Var_Desc = ivd;
      if (idof)
        *idof = I_offset - bottom;
      return vd;
    }
    bottom += vd->Ndof;
  }
  GOMA_EH(GOMA_ERROR, "Index_Solution_Inv ERROR");
  return NULL;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void dofname40(const int unknown, char *str)

/*********************************************************************
 *
 * dofname40:
 *
 *  This routine creates a descriptive string describing a variable
 *  The length of the string will not exceed 40 characters
 *
 * Input
 * ---------
 *  unknown -> processor unknown number
 *
 *  Output
 * ---------
 *  str -> return string (must have been previously allocated with
 *         space for 40 characters);
 *********************************************************************/
{
  VARIABLE_DESCRIPTION_STRUCT *vd = NULL;
  char position_label[40], *str_pos = str;
  int inode = 0, i_Var_Desc, i_offset, idof, var_type, kspec, matID, len = 0;
  if (str == NULL)
    return;
  *str = '\0';
  if (unknown < 0)
    return;

  vd = Index_Solution_Inv(unknown, &inode, &i_Var_Desc, &i_offset, &idof, pg->imtrx);
  if (vd == NULL)
    return;
  var_type = vd->Variable_Type;
  matID = vd->MatID;

  /*
   *  Add the name of the processor if we are doing a multiprocessor run
   */
  if (Num_Proc > 1) {
    sprintf(str, "P_%d, ", ProcID);
    len += strlen(str);
    str += strlen(str);
  }
  if (var_type == MASS_FRACTION) {
    kspec = vd->Subvar_Index;
    sprintf(str, "%d %s_%d n=%d", unknown, Var_Name[var_type].name2, kspec, inode + 1);
  } else {
    sprintf(str, "%d %s n=%d", unknown, Var_Name[var_type].name2, inode + 1);
  }
  len += strlen(str);
  str += strlen(str);
  if (matID != -1) {
    sprintf(str, "Mn = %d", matID);
    len += strlen(str);
    str += strlen(str);
  }
  switch (Num_Dim) {
  case 1:
    sprintf(position_label, " (x=%.3g)", Coor[0][inode]);
    break;
  case 2:
    sprintf(position_label, " (x=%.3g y=%.3g)", Coor[0][inode], Coor[1][inode]);
    break;
  case 3:
    sprintf(position_label, " (x=%.3g y=%.3g z=%.3g)", Coor[0][inode], Coor[1][inode],
            Coor[2][inode]);
    break;
  }
  if (len < 39) {
    strncat(str_pos, position_label, 39 - len);
  }
  str_pos[39] = '\0';
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
