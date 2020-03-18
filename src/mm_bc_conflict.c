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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "std.h"
#include "rf_fem_const.h"
#include "rf_fem.h"
#include "rf_io_const.h"
#include "rf_io.h"
#include "rf_mp.h"
#include "el_elm.h"
#include "rf_allo.h"
#include "rf_bc.h"
#include "mm_eh.h"
#include "exo_struct.h"
#include "dpi.h"
#include "rf_vars_const.h"
#include "mm_as_structs.h"
#include "mm_as.h"
#include "mm_as_const.h"
#include "rf_node_const.h"
#include "rd_mesh.h"
#include "mm_bc.h"
#include "dp_utils.h"
#include "el_elm_info.h"
#include "rf_bc_const.h"
#ifndef MAX_NODAL_BCS
#define MAX_NODAL_BCS  35
#endif

/*
 * These two workhorse variables help get more specific information out
 * to the generic error handler. They're declared here so they may be used
 * in every routine in this file without worrying about local declaration.
 */



/*
 * Variable DEFINITIONS...
 */

int **BC_dup_nodes = NULL;
int ****BC_dup_list = NULL;
int *BC_dup_ptr = NULL;
int **mesh_rotate_node = NULL;
int **mesh_rotate_ss = NULL;
int *num_mesh_rotate;
int **mom_rotate_node = NULL;
int **mom_rotate_ss = NULL;
int *num_mom_rotate;

/********** R O U T I N E S   D E F I N E D   I N   T H I S   F I L E **********
*
*     NAME OF FUNCTION			TYPE 		CALLED BY
* ------------------------       ---------------     -------------------

*
*******************************************************************************/

/** P R O T O   D E F I N I T I O N S   O  F   S T AT I C   F U N C T IO N S **/

/*****************************************************************************/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static int **alloc_bc_unk_list_node(int inode)
    
    /**************************************************************************
     *
     * alloc_bc_unk_list_node():
     *
     *  This function allocates a vector of pointers to ints of
     *  length MAX(nv->Num_Unknowns, 1) and then populates
     *  those pointers by allocating a vector of ints of length
     *  MAX_NODAL_BCS. The resulting array looks like an integer
     *  array of size:
     *
     *          bc_list_node[nv->Num_Unknowns][MAX_NODAL_BCS]
     *
     *  The elements of the array are initialized to the value of -1.
     *************************************************************************/
{
  int unk;
  NODE_INFO_STRUCT *node = Nodes[inode];
  NODAL_VARS_STRUCT *nv = node->Nodal_Vars_Info[pg->imtrx];
  int num = nv->Num_Unknowns;
  int **bc_unk_list_node;

  if (num <= 0 && upd->Total_Num_Matrices > 1) {
    num = 1;
  }
  bc_unk_list_node = (int **) alloc_ptr_1(num);
  for (unk = 0; unk < num; unk++) {
    bc_unk_list_node[unk] = alloc_int_1(MAX_NODAL_BCS, -1);
  }

  int sum_unknowns = 0;
  int imtrx;
  for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
    NODAL_VARS_STRUCT *nv = node->Nodal_Vars_Info[imtrx];
    sum_unknowns += nv->Num_Unknowns;
  }

  if (sum_unknowns <=  0) {
    fprintf(stderr,
            "P_%d: Warning: node %d with zero unknowns has an applied bc\n",
            ProcID, inode);
  }

  return bc_unk_list_node;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void free_bc_unk_list_node(int ***bc_list_node_ptr, int inode)
    
    /**************************************************************************
     *
     * free_bc_list_node():
     *
     *  This function frees the malloced arrays created by
     *  alloc_bc_list_node(); All previously malloced pointers are set
     *  back to the NULL value for better error checking
     *
     * Input
     *
     * bc_list_node_ptr =  pointer to the structure malloced by
     *                     the function alloc_bc_list_node().
     *************************************************************************/
{
  int unk;
  NODE_INFO_STRUCT *node = Nodes[inode];
  NODAL_VARS_STRUCT *nv = node->Nodal_Vars_Info[pg->imtrx];
  int num = nv->Num_Unknowns;
  int **bc_list_node = *bc_list_node_ptr;
  for (unk = 0; unk < num; unk++) {
    safer_free((void **) (bc_list_node + unk));
  }
  safer_free((void **) bc_list_node_ptr);
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static int
find_first_opening(const int *list, const int lengthList)
    
    /**************************************************************************
     *
     * find_first_opening():
     *
     *  This function find the first entry in an integer list that is
     *  equal to negative one. An error exit from the code is
     *  created if the end of the list is encountered with all entries
     *  not equal to -1.
     *
     * Input
     *
     * list[]     = integer list
     * lengthList = Length of the list
     *
     * Output
     *
     * return -> first index in the array corresponding to a number which
     *           is equal to -1 (the default fill value in this case).
     *************************************************************************/
{
  int i;
  for (i = 0; i < lengthList; i++) {
    if (list[i] == -1) return i;
  }
  /*
   *  create an error exit here
   */
  EH(GOMA_ERROR, "too many duplications");
  return -1;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void
delete_bc_entry(int *bc_list, const int j)

    /*
     *
     */
{
  int i;
  bc_list[j] = -1;
  for (i = j; i < MAX_NODAL_BCS-1 && (bc_list[i+1] != -1); i++) {
    bc_list[i] = bc_list[i+1];
  }
} 
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
static void
compress_bc_list(int *bc_list) 

{
  int i, j = 0;
  for (i = 0; i < MAX_NODAL_BCS; i++) {
    if (bc_list[i] >= 0) {
      bc_list[j] = bc_list[i];
      if (i != j) bc_list[i] = -1;
      j++;
    }
  }
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
static void
add_to_node_unk_bc_list(int *list, int ibc)
{
  int i;
  for (i = 0; i < MAX_NODAL_BCS; i++) {
    if (list[i] == -1) {
      list[i] = ibc;
      return;
    }
    if (list[i] == ibc) {
      return;
    }
  }
  /*
   *  create an error exit here
   */
  EH(GOMA_ERROR, "add_to_node_unk_bc_list out of room");
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int
search_bc_dup_list(const int bc_item, int *list)
{
  int i;
  for (i = 0; list[i] != -1; i++) {
    if (list[i] == bc_item) return i;
  } 
  return -1;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void 
check_for_bc_conflicts2D(Exo_DB *exo, Dpi *dpi)

    /*************************************************************************
   *
   * check_for_bc_conflicts2D():
   *
   * This routine takes care of several tasks up front for the application
   * of boundary conditions.  Among other things, 
   * It looks for multiple bc's being applied to the same equation at the 
   * same node, and sets up criteria for deciding which BC has precedence
   * Last revised by Rich Cairncross 12/05/95
   ***************************************************************************/
{
  int es, i, ibc, iss, num_bc_nodes, used_BC, varType;
  int ibc1, ibc2, inode, bct1, bct2=-1, k;
  int ins, action, idum, imove, istay, i_delete, i_dont_move;
  int eqn, idup, idup1, idup2, p, j, num_total_nodes, dups;
  int ***BC_Unk_List;
  double sum, a1, a2, a3, a4, a5, a6, a7, b1, b2, b3, b4, b5, b6, b7, max, c1;
  double  c2, det, det1, det2, det3;
  int ivar, jvar, var=-1;
  int count_BC, count_DC, count_rotate, count_weak, count_strong;
  int count_coll, count_special;
  int ss_rot, dir_rot, i_rotate, two_ss_rot;
  int divert;			/* output diag of duplications to file */
  int bickel;			/* BC list temp variable */
  char bc_divert_fn[MAX_FNL];
  FILE *bc_dup_out;
  int matIndex = -1, offset, retn_matIndex;
  int ielem, elem_block_index, offset_p;
  int offset_mom1, offset_mom2, offset_mesh1, offset_mesh2, ndofs, idof;
  VARIABLE_DESCRIPTION_STRUCT *vd_retn1, *vd_retn2;
  NODE_INFO_STRUCT *node;
  NODAL_VARS_STRUCT *nv;
  VARIABLE_DESCRIPTION_STRUCT *vd;
  struct BC_descriptions *bc_desc;

  static const char yo[] = "check_for_bc_conflicts2D";

  action   = -1;
  i_delete = -1;

  BC_dup_nodes = malloc(sizeof(int *) * (size_t) upd->Total_Num_Matrices);
  BC_dup_list = malloc(sizeof(int ***) * (size_t) upd->Total_Num_Matrices);
  BC_dup_ptr = malloc(sizeof(int) * (size_t) upd->Total_Num_Matrices);

/*****************************************************************************/
/*                              BLOCK 3                                      */
/*      Start Checking for Duplicate Boundary Conditions at the same node    */
/*****************************************************************************/

  num_total_nodes = ( dpi->num_internal_nodes + 
		      dpi->num_boundary_nodes +
		      dpi->num_external_nodes);


  /***************************************************************************
   *
   * FIRST make a list of all the boundary conditions from the input deck
   *  at each node on the current processor.
   *
   ***************************************************************************/
/*****************************************************************************/
/*      BOUNDARY CONDITIONS SPECIFIED BY NODE SETS                           */
/*  visit each node on all node sets and list the BC's applied there         */
/*****************************************************************************/

  int *matrix_used_BC = NULL;

  if (ProcID == 0 )
    {
      matrix_used_BC = calloc((size_t) Num_BC, sizeof(int));
    }

  for (pg->imtrx = 0; pg->imtrx < upd->Total_Num_Matrices; pg->imtrx++) {
  /* 
   * BC_list, what is it? A particular equation at a node can look to get
   * multiple specifications of what to do about the boundary condition.
   * All the possibilities are saved in BC_list[][][] and then combed through
   * to determine precedence order and which to keep, which to discard.
   * Lots of heuristics determine what to do - see below.
   *
   * BC_list[node_number][equation_number][potential_boundary_condition] = 
   * boundary_condition_index
   * 
   * create an array which will hold boundary condition information temporarily
   * for each node
   */
    BC_Unk_List = (int ***) alloc_ptr_1(num_total_nodes);
    num_bc_nodes = 0;

    for (ibc = 0; ibc < Num_BC; ibc++) {
      /* First check to see if this boundary condition is a node set */
      if (!strcmp(BC_Types[ibc].Set_Type, "NS")) {
        /*
         * Set a flag to indicate that a boundary condition is used
         */
        used_BC = 0;
        matIndex = -2;
        if (BC_Types[ibc].BC_EBID_Apply != -1) {
          matIndex = map_mat_index(BC_Types[ibc].BC_EBID_Apply);
        }
        /*
         * Loop over the total number of node sets defined for the
         * current processor
         */
        for (ins = 0; ins < exo->num_node_sets; ins++) {
          /* Check for a match between the ID of the current node set
           * and the node set ID specified in the input file
           *  - continue if a match is found
           */
          if (exo->ns_id[ins] == BC_Types[ibc].BC_ID) {

            /* Loop over all the global nodes in the current node set */
            for (i = 0; i < exo->ns_num_nodes[ins]; i++) {

              /* Get the ith local node in the current node set */
              inode = exo->ns_node_list[exo->ns_node_index[ins] + i];
              /*
               *  Allocate temp BC storage for this node and initialize
               *  the storage to -1
               */
              if (BC_Unk_List[inode] == NULL) {
                num_bc_nodes++;
                BC_Unk_List[inode] = alloc_bc_unk_list_node(inode);
              }
              /*
               * Fill in an entry into the list for each component of the
               * boundary condition
               */
              for (p = 0; p < BC_Types[ibc].desc->vector; p++) {
                offset = find_bc_unk_offset(BC_Types + ibc, matIndex, inode, p,
                                            &retn_matIndex, &vd_retn1);
                if (offset >= 0) {
                  add_to_node_unk_bc_list(BC_Unk_List[inode][offset], ibc);
                  used_BC = 1;
                }
              }
            } /* for (i = 0; i < exo->ns_num_nodes[ins]; i++) { */
          } /* if (exo->ns_id[ins] == BC_Types[ibc].BC_ID) { */
        } /* for (ins = 0; ins < exo->num_node_sets; ins++) { */
        used_BC = gsum_Int(used_BC);
        if (ProcID == 0)
          {
            matrix_used_BC[ibc] |= used_BC;
          }


      }  /* END if (!strcmp(BC_Types[ibc].Set_Type, "NS"))   */
    } /* for (ibc = 0; ibc < Num_BC; ibc++) */

    /***************************************************************************/
    /*      BOUNDARY CONDITIONS SPECIFIED BY SIDE SETS                         */
    /*  visit each node on all side sets and list the BC's applied there       */
    /***************************************************************************/

    for (ibc = 0; ibc < Num_BC; ibc++) {
      /* First check to see if this boundary condition is a side set */
      if (!strcmp(BC_Types[ibc].Set_Type, "SS")) {
        /*
         * Set a flag to indicate that a boundary condition is used
         */
        used_BC = 0;
        /* Loop over the total number of side sets defined on the current processor */
        for (iss = 0; iss < exo->num_side_sets; iss++) {
          /* Check for a match between the ID of the current side set and the side set
           * ID specified in the input file - continue if a match is found 
           */
          if (exo->ss_id[iss] == BC_Types[ibc].BC_ID) {

            /*
             * Step through every elem/side and look at every node on that
             * elem/side...
             */
            for (es = 0; es < exo->ss_num_sides[iss]; es++) {
              /*
               * Find the element, element_block_index, and the material
               * index corresponding to this side
               */
              ielem = exo->ss_elem_list[exo->ss_elem_index[iss] + es];
              elem_block_index = find_elemblock_index(ielem, exo);
              matIndex = Matilda[elem_block_index];
              /*
               * Loop over all of the local nodes on the side 
               */
              for (k = exo->ss_node_side_index[iss][es]; 
                   k < exo->ss_node_side_index[iss][es+1]; k++) {
                /*
                 * Find the global node index corresponding to the local
                 * node index
                 */
                inode = exo->ss_node_list[iss][k];

                /* allocate temp BC storage for this node and initialize*/
                if (BC_Unk_List[inode] == NULL) {
                  num_bc_nodes++;
                  BC_Unk_List[inode] = alloc_bc_unk_list_node(inode);
                }
           
                /*
                 * Fill in an entry into the list for each component of the
                 * boundary condition
                 */
                for (p = 0; p < BC_Types[ibc].desc->vector; p++) {
                  offset = find_bc_unk_offset(BC_Types + ibc, matIndex, inode, p,
                                              &retn_matIndex, &vd_retn1);
                  if (offset >= 0) {
                    add_to_node_unk_bc_list(BC_Unk_List[inode][offset], ibc);
                    used_BC = 1;
                  }
                }

                /*FIRST REAL SPECIAL CASE */
                /* Now, at this node and check and see if on */
                /* R_MOMENTUM1 and R_MOMENTUM2 positions there */
                /* are VELO_NORMAL and a VELO_TANGENT coming from */
                /* different side sets. Here we must make sure that */
                /* the point actually makes the distinguished dup list*/
		    
                node = Nodes[inode];
                nv = node->Nodal_Vars_Info[pg->imtrx];
                ndofs = get_nv_ndofs(nv, R_MOMENTUM1);
                for (idof = 0; idof < ndofs; idof++) {
                  if ( idof >= nv->Num_Var_Desc_Per_Type[R_MOMENTUM1] ) continue; /*check needed for XFEM */
                  offset_mom1 = get_nv_offset_idof(nv, R_MOMENTUM1, idof, 0, &vd_retn1);
                  offset_mom2 = get_nv_offset_idof(nv, R_MOMENTUM2, idof, 0, &vd_retn2);
                  if (vd_retn1->MatID != vd_retn2->MatID) {
                    EH(GOMA_ERROR,"Unforseen occurrence");
                  }
                  if (offset_mom1 < 0) break;
                  idup1 = find_first_opening(BC_Unk_List[inode][offset_mom1],
                                             MAX_NODAL_BCS);
                  idup2 = find_first_opening(BC_Unk_List[inode][offset_mom2],
                                             MAX_NODAL_BCS);
                  for (i = 0; i < idup1; i++) {
                    ibc1 = BC_Unk_List[inode][offset_mom1][i];
                    if (BC_Types[ibc1].BC_Name == VELO_NORMAL_BC ||
                        BC_Types[ibc1].BC_Name == VELO_NORMAL_LS_BC ) {
                      for (j = 0; j < idup2; j++) {
                        ibc2 = BC_Unk_List[inode][offset_mom2][j]; 
                        if (BC_Types[ibc2].BC_Name == VELO_TANGENT_BC ||
                            BC_Types[ibc2].BC_Name == VELO_TANGENT_USER_BC ||
                            BC_Types[ibc2].BC_Name == VELO_STREAMING_BC ) {
                          if(BC_Types[ibc2].BC_ID != BC_Types[ibc1].BC_ID) {
                            /* Now you have established a velo_normal and velo_tangent
                             *  coming from two separate side sets, at the same node.  We know
                             *  if these sides meet at 90 deg., the conditions are redundant
                             *  so you must do whatever you have to here to make sure that
                             *  one is dropped.  So, what we will do is make all of the current
                             *  entries to include the velo_normal, and, to make sure that it
                             *  makes the dup_list, we need two entries so we add another
                             *  to the end of the list
                             */
                            BC_Unk_List[inode][offset_mom1][idup1] = ibc1;
                          }
                        }
                      }
                    }
                  }			
                }
              }   /* end of nodes per side-set-side loop k */
            }

            /* since we've found the matching side set for this BC*/

            /* Special Case */
            /* Check for table boundary conditions which may be used by */
            /* a discontinuous polymer stress formulation on the inlet   */
 
            if(BC_Types[ibc].BC_Name == TABLE_BC)
              {
                varType=BC_Types[ibc].table->f_index;
                switch(varType)
                  {
                  case POLYMER_STRESS11:
                  case POLYMER_STRESS12:
                  case POLYMER_STRESS22:
                  case POLYMER_STRESS11_1:
                  case POLYMER_STRESS12_1:
                  case POLYMER_STRESS22_1:
                  case POLYMER_STRESS11_2:
                  case POLYMER_STRESS12_2:
                  case POLYMER_STRESS22_2:
                  case POLYMER_STRESS11_3:
                  case POLYMER_STRESS12_3:
                  case POLYMER_STRESS22_3:
                  case POLYMER_STRESS11_4:
                  case POLYMER_STRESS12_4:
                  case POLYMER_STRESS22_4:
                  case POLYMER_STRESS11_5:
                  case POLYMER_STRESS12_5:
                  case POLYMER_STRESS22_5:
                  case POLYMER_STRESS11_6:
                  case POLYMER_STRESS12_6:
                  case POLYMER_STRESS22_6:
                  case POLYMER_STRESS11_7:
                  case POLYMER_STRESS12_7:
                  case POLYMER_STRESS22_7:
                    {
                      /* used material index determined for last element on the current side set */
                      switch(pd_glob[matIndex]->i[pg->imtrx][varType])
                        {
                        case I_P0:
                        case I_P1:
                        case I_PQ1:
                        case I_PQ2:
                          {
                            fprintf(stderr,
                                    "WARNING: Boundary condition %d, %s, applied on SS %d, may be used as an inlet condition \n",
                                    ibc, BC_Types[ibc].desc->name1, BC_Types[ibc].BC_ID);
                            break;
                          }
                        default:
                          {
                            break;
                          }
                        }
                      break;
                    }
		
                  default:
                    {
                      break;
                    }
                  }
              }


          }  /* if (Proc_SS_Ids[iss] == BC_Types[ibc].BC_ID)  */
        }  /*  for (iss = 0; iss < Proc_Num_Side_Sets; iss++) */
        /*
         */



        used_BC = gsum_Int(used_BC);
        if (ProcID == 0)
          {
            matrix_used_BC[ibc] |= used_BC;
          }
      }  /* END if (!strcmp(BC_Types[ibc].Set_Type, "SS")) 		      */
    }  /* END for (ibc = 0; ibc < Num_BC; ibc++)				      */

    /****************************************************************************/
    /*      RESOLVE BOUNDARY CONDITION CONFLICTS SPECIFIED                      */
    /****************************************************************************/
    /* I know this section looks like a real mess, and I'm sorry that it'll be hard
     * to follow.  Should anyone want to make changes here, or try to understand what
     * logic (?) is built into this routine - see the BC conflict resolution flow-chart
     * in the GOMA manual documentation and try to believe that this routine follows 
     * that flow chart somewhat religiously
     */


    /*
     * When running distributed, some processors have small chunks of the overall
     * problem with absolutely no boundary conditions. This should be permissible.
     */

    if (num_bc_nodes > 0) {
      BC_dup_nodes[pg->imtrx] = alloc_int_1(num_bc_nodes, -1);
      BC_dup_list[pg->imtrx] = (int ***) alloc_ptr_1(num_bc_nodes);
    }

    /* resolve conflicts between BC's at each node */
    /* first initialize the BC duplications list*/
    BC_dup_ptr[pg->imtrx] = 0;

    /* 
     * loop through nodes and remove nodes that don't have conflicts 
     * count nodes with conflicts and build an array to store them 
     */
    for (inode = 0; inode < num_total_nodes; inode++) {
      if (BC_Unk_List[inode] != NULL) {
        node = Nodes[inode];
        nv = node->Nodal_Vars_Info[pg->imtrx];
        /* check for duplications */
        dups = 0;
        for (offset = 0; offset < nv->Num_Unknowns; offset++) {
          get_nv_vd_from_offset(nv, offset, &vd, &idof);
          eqn = vd->Variable_Type;
          /* count duplications */
          idup = find_first_opening(BC_Unk_List[inode][offset], MAX_NODAL_BCS);
          /*
           * RESOLVE DUPLICATIONS AT THIS NODE
           */
          if (idup > 1) {
            /* initialize new entry in the duplications list, if not already listed */
            if (in_list(inode, 0, BC_dup_ptr[pg->imtrx], BC_dup_nodes[pg->imtrx]) == -1) {
              /* Put the node in the list */
              BC_dup_nodes[pg->imtrx][BC_dup_ptr[pg->imtrx]] = inode;
              /* Put the pointer to the BC_Unk_List for the node
               * into another list
               */
              BC_dup_list[pg->imtrx][BC_dup_ptr[pg->imtrx]] = BC_Unk_List[inode];
              BC_dup_ptr[pg->imtrx]++;
              if (Debug_Flag > 0) {
                fprintf(stderr, "New node in duplication list %d %d\n",inode+1,BC_dup_ptr[pg->imtrx]);
              }
            }

            if (Debug_Flag) {
              idup = 0;
              fprintf(stderr, "  DUPLICATE BC at node %d equation %d, material %d: ", 
                      inode+1, eqn, vd->MatID);
              while ((ibc = BC_Unk_List[inode][offset][idup]) != -1) {
                fprintf(stderr, "%s, ", BC_Types[ibc].desc->name2);
                idup++;
              }
              fprintf(stderr, "\n");
            }
            /* increment count of # duplications at this node */
            dups++;

            /*****************************************************************************/
            /* try to resolve the differences quickly */
            /* when two BC's are applied at a point, there are three possibilities
             * for how to choose to apply the BC's:
             *  1- One of the BC's takes precedence and the other is discarded
             *     regardless of their coefficients
             *     - the first one in the input deck takes precedence
             *  DIRICHLET Conditions and UNROTATED POINTWISE COLLOCATED Conditions
             *  2- Both BC's get applied - regardless of their coefficients
             *  3- One of the BC's shifts to another equation if it's coefficients
             *     are different
             *  BUT FIRST LOOK AT THE REAL PATHOLOGICAL VELO_NORMAL/VELO_TANGENT CASES
             */
            /*****************************************************************************/
            /*FIRST REAL SPECIAL CASE */
            /* Now, at this node and check and see if on */
            /* R_MOMENTUM1 and R_MOMENTUM2 positions there */
            /* are VELO_NORMAL and a VELO_TANGENT coming from */
            /* different side sets. */
	     
            if (eqn == R_MOMENTUM1) {
              offset_mom1 = offset;
              idup1 = find_first_opening(BC_Unk_List[inode][offset_mom1],
                                         MAX_NODAL_BCS);
              /*
               * Get the corresponding offset for the v velocity for the
               * current material. vd still contains the description
               * of the current u velocity unknown.
               */
              offset_mom2 = 
		get_nodal_unknown_offset(nv, R_MOMENTUM2, vd->MatID, 0, &vd_retn2);
              idup2 = find_first_opening(BC_Unk_List[inode][offset_mom2],
                                         MAX_NODAL_BCS);
	      
              for (i = 0; i < idup1; i++) {
                ibc1 = BC_Unk_List[inode][offset_mom1][i];
                if (BC_Types[ibc1].BC_Name == VELO_NORMAL_BC ||
                    BC_Types[ibc1].BC_Name == VELO_NORMAL_LS_BC ) {
                  for (j = 0; j < idup2; j++) {
                    ibc2 = BC_Unk_List[inode][offset_mom2][j];
                    if (ibc2 != -1) {
                      if (BC_Types[ibc2].BC_Name == VELO_TANGENT_BC ||
                          BC_Types[ibc2].BC_Name == VELO_TANGENT_USER_BC ||
                          BC_Types[ibc2].BC_Name == VELO_STREAMING_BC ) {
                        /* Hmm, we got a hit. 
                         * Now make sure they are not coming from the same
                         * side set, which is something that is common and
                         * permissible.
                         */
                        if (BC_Types[ibc2].BC_ID != BC_Types[ibc1].BC_ID) {
                          /* Now you have established a velo_normal and velo_tangent
                           *  coming from two separate side sets, at the same node.  We know
                           *  if these sides meet at 90 deg., the conditions are redundant
                           *  so you must do whatever you have to here to make sure that
                           *  one is dropped.  So what we will do is make all of the current
                           *  entries to include the velo_normal, and to make sure that it
                           *  makes the dup_list we need two entries.
                           */	  
                          delete_bc_entry(BC_Unk_List[inode][offset_mom2], j);
                        }
                      }
                    }
                  }
                }
              }
	      	      
              /*SECOND REAL SPECIAL CASE */
              /* Now, at this node and check and see if on */
              /* R_MESH1 and R_MOMENTUM2 positions there */
              /* are KINEMATIC and a VELO_TANGENT coming from */
              /* different side sets  */
              offset_mesh1 =   
		get_nodal_unknown_offset(nv, R_MESH1, vd->MatID, 0, &vd_retn1);
              if (offset_mesh1 < 0) {
                idup1 = 0;
              } else {
                idup1 = find_first_opening(BC_Unk_List[inode][offset_mesh1],
                                           MAX_NODAL_BCS);
              }
              if (offset_mom2 < 0) {
                idup2 = 0;
              } else {
                idup2 = find_first_opening(BC_Unk_List[inode][offset_mom2],
                                           MAX_NODAL_BCS);
              }
              for (i = 0; i < idup1; i++) {
                ibc1 = BC_Unk_List[inode][offset_mesh1][i];
                if (BC_Types[ibc1].BC_Name == KINEMATIC_BC ||
                    BC_Types[ibc1].BC_Name == KINEMATIC_DISC_BC ||
                    BC_Types[ibc1].BC_Name == KINEMATIC_PETROV_BC ||
                    BC_Types[ibc1].BC_Name == KINEMATIC_COLLOC_BC) {
                  for (j = 0; j < idup2; j++) {
                    ibc2 = BC_Unk_List[inode][offset_mom2][j];
                    if (BC_Types[ibc2].BC_Name == VELO_TANGENT_BC ||
 		  	BC_Types[ibc2].BC_Name == VELO_TANGENT_USER_BC ||
 		  	BC_Types[ibc2].BC_Name == VELO_SLIP_FLUID_BC ||
 		  	BC_Types[ibc2].BC_Name == VELO_SLIP_ROT_FLUID_BC ||
 			BC_Types[ibc2].BC_Name == VELO_STREAMING_BC ) {
                      /* Hmm, we got a hit.  Now make sure they are not coming from the same
                         side set, which is something that is common */
                      if (BC_Types[ibc2].BC_ID != BC_Types[ibc1].BC_ID) {
                        /* Now you have established a KINEMATIC and velo_tangent
                         *  makes the dup_list, lets get rid of the velo_tangent, 
                         *  unless we have a rolling motion condition at a dynamic
                         *  contact line.  In that case
                         *  we would like to retain the velo_tangent and the kinematic
                         *  as they are not redundant.  The user signifies this case
                         *  with a -1 in the Data_int[0] slot
                         */
		      if (BC_Types[ibc2].BC_Data_Int[0] > -1)  /*this is a last defense
                                                                    to retain velo_tangent*/
                          {
                            delete_bc_entry(BC_Unk_List[inode][offset_mom2], j);
                          }
                      }
                    }
                  }
                }
              }
            }  /* End of if (eqn == R_MOMENTUM1) */
            /*****************************************************************************/
	  
            /* first give all DIRICHLET Conditions Precedence 
             * - they swamp out all other boundary conditions and shift rotated
             *   conditions to another coordinate direction
             */
            for (j = 0; j < idup;  j++) {
              ibc1 = BC_Unk_List[inode][offset][j];
              if (ibc1 != -1) {
                bct1 = BC_Types[ibc1].BC_Name;
                if (BC_Types[ibc1].desc->method == DIRICHLET) {
                  for (k = 0; k < idup; k++) {
                    /* flag indicating that the k bc should be deleted from the
                     * list */
                    i_delete = 1;
                    ibc2 = BC_Unk_List[inode][offset][k];
                    if (ibc2 != -1 && ibc1 != ibc2) {
                      bct2 = BC_Types[ibc2].BC_Name;

                      /* deal with identical conditions */
                      sum = 0.0;
                      max = 0.0;
                      if (bct1 == bct2) { /* these conditions are the same type */
                        if (bct1 != FIX_BC) {
                          sum += fabs(BC_Types[ibc1].BC_Data_Float[0]
                                      - BC_Types[ibc2].BC_Data_Float[0]); 
                          if (fabs(BC_Types[ibc1].BC_Data_Float[0]) > max) 
			    max = fabs(BC_Types[ibc1].BC_Data_Float[0]);
                          if (fabs(BC_Types[ibc2].BC_Data_Float[0]) > max) 
			    max = fabs(BC_Types[ibc2].BC_Data_Float[0]);
                        }
                        /* sum the differences between coefficients and normalize*/
                        if (max > 0. && sum/max > 1e-8) {
                          if (Debug_Flag > 0) {
                            fprintf(stderr, 
                                    "    TWO BC's %s at node %d have conflicting values, "
                                    " defaulting to 1st occurence = %g\n", 
                                    BC_Types[ibc1].desc->name2,
                                    inode, BC_Types[ibc2].BC_Data_Float[0]);
                          }
                          i_delete = 2; /* delete 2nd occurence */
                        } else {
                          if (Debug_Flag > 0) {
                            fprintf(stderr, 
                                    "    IDENTICAL BC %s at node %d equation %d\n", 
                                    BC_Types[ibc1].desc->name2,
                                    inode+1, eqn);
                          }
                        }
                        action = 1;	/* discard extra bc */
                      }

                      /* move rotated conditions */
                      if (BC_Types[ibc2].desc->rotate != NO_ROT ||
                          BC_Types[ibc2].desc->BC_Name == SDC_STEFANFLOW_BC) {
                        /*  move this condition ibc2 to the next coordinate direction, 
                         *  if available 
                         * HKM -> Note this seems to preclude that a tangential
                         *   mesh condition can move to the x coordinate direction
                         *   due to a dv dirichlet condition ?!??
                         */
                        if (BC_Types[ibc2].desc->rotate == R_MESH1) {
                          p = eqn - R_MESH1;
                          if (p + 1 < pd_glob[0]->Num_Dim) {
                            offset_mesh2 = 
			      get_nodal_unknown_offset(nv, eqn+1, vd->MatID, 0,
						       &vd_retn2);
                            idum = find_first_opening(BC_Unk_List[inode][offset_mesh2],
                                                      MAX_NODAL_BCS);
                            if (Debug_Flag > 0) {
                              fprintf(stderr, "    MOVING BC %s at node %d because of %s\n", 
                                      BC_Types[ibc2].desc->name2, 
                                      inode+1, BC_Types[ibc1].desc->name2);
                            }
                            BC_Unk_List[inode][offset_mesh2][idum] = ibc2;
                          } else {
                            if (Debug_Flag > 0) {
                              fprintf(stderr,
                                      "    REMOVING BC %s at node %d because of %s\n", 
                                      BC_Types[ibc2].desc->name2, 
                                      inode+1, BC_Types[ibc1].desc->name2);
                            }
                          }
                        } else if (BC_Types[ibc2].desc->rotate == R_MOMENTUM1) {
                          p = eqn - R_MOMENTUM1;
                          if (p+1 < pd_glob[0]->Num_Dim) {
                            offset_mom2 = 
			      get_nodal_unknown_offset(nv, eqn+1, vd->MatID, 0,
						       &vd_retn2);
                            idum = find_first_opening(BC_Unk_List[inode][offset_mom2],
                                                      MAX_NODAL_BCS);
	       
                            if (Debug_Flag > 0) {
                              fprintf(stderr,
                                      "    MOVING BC %s at node %d because of %s\n", 
                                      BC_Types[ibc2].desc->name2, 
                                      inode+1, BC_Types[ibc1].desc->name2);
                            }
                            BC_Unk_List[inode][offset_mom2][idum] = ibc2;
                          } else {
                            if (Debug_Flag > 0) {
                              fprintf(stderr,
                                      "    REMOVING BC %s at node %d because of %s\n", 
                                      BC_Types[ibc2].desc->name2, 
                                      inode+1, BC_Types[ibc1].desc->name2);
                            }
                          }
                        } else if (BC_Types[ibc2].desc->rotate == R_SOLID1 ) {
                          p = eqn - R_SOLID1;
                          if (p+1 < pd_glob[0]->Num_Dim) {
                            offset_mom2 = 
			      get_nodal_unknown_offset(nv, eqn+1, vd->MatID, 0,
						       &vd_retn2);
                            idum = find_first_opening(BC_Unk_List[inode][offset_mom2],
                                                      MAX_NODAL_BCS);
	       
                            if (Debug_Flag > 0) {
                              fprintf(stderr,
                                      "    MOVING BC %s at node %d because of %s\n", 
                                      BC_Types[ibc2].desc->name2, 
                                      inode+1, BC_Types[ibc1].desc->name2);
                            }
                            BC_Unk_List[inode][offset_mom2][idum] = ibc2;
                          } else {
                            if (Debug_Flag > 0) {
                              fprintf(stderr,
                                      "    REMOVING BC %s at node %d because of %s\n", 
                                      BC_Types[ibc2].desc->name2, 
                                      inode+1, BC_Types[ibc1].desc->name2);
                            }
                          }
                        } else {
                          WH(-1, "Missing type of rotation in BC conflict resolution");
                        }
                      }

                      /* remove this BC from list */
                      if (i_delete == 1) {
                        if (Debug_Flag > 0) {
                          fprintf(stderr,
                                  "    DISCARDING %s at node %d from dup_list because of %s\n", 
                                  BC_Types[ibc2].desc->name2, inode+1,
                                  BC_Types[ibc1].desc->name2);
                        }
                        BC_Unk_List[inode][offset][k] = -1;
                      }
                      if (i_delete == 2) {
                        if (Debug_Flag > 0) {
                          fprintf(stderr,
                                  "    DISCARDING %s at node %d from dup_list because of %s\n", 
                                  BC_Types[ibc1].desc->name2, inode+1, 
                                  BC_Types[ibc2].desc->name2);
                        }
                        BC_Unk_List[inode][offset][j] = -1;
                      }
                    } /* if ((ibc2 != ibc1) && (ibc2 != -1)) */
                  }  /* for (k=j+1; k<idup; k++) */
                } /* if (BC_Types[ibc1].desc->method == DIRICHLET) */
              } /* if (ibc1 != -1) */
            } /*  for (j=0; j<idup-1; j++)*/

            /*
             * Compress the list and find the number of bcs
             */
            compress_bc_list(BC_Unk_List[inode][offset]);
            idup = find_first_opening(BC_Unk_List[inode][offset], MAX_NODAL_BCS);
            /**********************************************************************/
            /* second give precedence to unrotated strong (integrated or 
             * collocated)  boundary conditions 
             * - they swamp out all other boundary conditions and shift rotated
             *   conditions to another coordinate direction 
             */
            for (j = 0; j < idup; j++) {
              ibc1 = BC_Unk_List[inode][offset][j];
              if (ibc1 != -1) {
                bct1 = BC_Types[ibc1].BC_Name;
                if (BC_Types[ibc1].desc->rotate == NO_ROT &&
                    (BC_Types[ibc1].desc->method == COLLOCATE_SURF  ||
                     BC_Types[ibc1].desc->method == STRONG_INT_SURF)) {
		
                  for (k = 0; k < idup; k++) {
                    ibc2 = BC_Unk_List[inode][offset][k];
                    if ((ibc2 != ibc1) && (ibc2 != -1)) {
                      bct2 = BC_Types[ibc2].BC_Name;
                      /* Flag indicating that this bc should be deleted */
                      i_delete = 1;

                      action = 1;
                      /* need to deal with GD conditions here 
                       *  - two GD conditions on the same SS and
                       *    applied to the same equation do not conflict
                       *    with each other 
                       */
                      if ((bct1 <= GD_TIME_BC && bct1 >= GD_CONST_BC) && 
                          (bct2 <= GD_TIME_BC && bct2 >= GD_CONST_BC)) {
                        if (BC_Types[ibc1].BC_ID == BC_Types[ibc2].BC_ID) {
                          action = -1;
                          if (Debug_Flag > 0) {
                            fprintf(stderr, "    MULTIPLE GD %s at node %d\n", 
                                    BC_Types[ibc2].desc->name2, inode+1);
                          }
                        }
                      }

                      /*
                       *  allow multiple FLUID_SOLID or SOLID_FLUID
                       *  conditions at one point 
                       */
                      if ((bct1 == FLUID_SOLID_BC && bct2 == FLUID_SOLID_BC) ||
                          (bct1 == SOLID_FLUID_BC && bct2 == SOLID_FLUID_BC) ||
                          (bct1 == SOLID_FLUID_RS_BC && bct2 == SOLID_FLUID_RS_BC)) {
                        action = -1;
                      }

                      /*
                       * Boundary conditions of the same time at 
                       * a node are simplified. The first one is kept, while
                       * the second one is discarded.
                       */
                      if (bct1 == bct2 && action == 1) {
                        if (Debug_Flag > 0) {
                          fprintf(stderr, 
                                  "    TWO BC's %s at node %d have same type, "
                                  " defaulting to 1st occurence = %g\n",
                                  BC_Types[ibc1].desc->name2,
                                  inode+1, BC_Types[ibc1].BC_Data_Float[0]);
                        }
                      }

                      /*
                       *  SPECIAL CASE:
                       *   Weakly integrated BC's in Sheep's Clothing --
                       *
                       *  The boundary conditions below are really weakly 
                       *  integrated boundary conditions. However, they are 
                       *  named COLLOCATE_SURF boundary conditions. They 
                       *  involve replacing the difference in the  
                       *  normal components of the stress between two phases
                       *  on opposite sides of a boundary with the volumetric
                       *  contribution of the stress on one side of the boundary.
                       *  Thus, they are really weakly integrated boundary
                       *  conditions. Because of this, other weakly integrated
                       *  boundary conditions should be allowed to be 
                       *  additively applied on top of them. 
                       */
                      if (bct1 == FLUID_SOLID_BC  || 
                          bct1 == SOLID_FLUID_BC  ||
                          bct1 == SOLID_FLUID_RS_BC) {
                        if (BC_Types[ibc2].desc->method == WEAK_INT_SURF ||
                            BC_Types[ibc2].desc->method == WEAK_INT_EDGE || 
                            BC_Types[ibc2].desc->method == WEAK_SHARP_INT ) {
                          action = -1;
                          if (Debug_Flag > 0) {
                            fprintf(stderr, 
                                    "    BC %s is a Weak BC masquerading as colloc BC"
                                    " Allow Weak BC %s to stack on top of it\n",
                                    BC_Types[ibc1].desc->name2, 
                                    BC_Types[ibc2].desc->name2);
                          }
                        }
                        /*
                         *  Also, because these are considered weakly integrated
                         *  Dirichlet conditions, other COLLOCATE_SURF and 
                         *  STRONG_INT_SURF boundary conditions will take 
                         *  precedence over these at the current node.
                         */
                        if (BC_Types[ibc2].desc->method == COLLOCATE_SURF ||
                            BC_Types[ibc2].desc->method == STRONG_INT_SURF) {
                          if (bct2 != bct1) {
                            i_delete = -1;
                            action = -1;
                          }
                        }
                      }

                      /*
                       *  move rotated conditions -> to the next coordinate
                       *  direction if available.
                       *   
                       */
                      if (action == 1 && 
                          BC_Types[ibc2].desc->rotate != NO_ROT) {
                        if (BC_Types[ibc2].desc->rotate == R_MESH1) {
                          p = eqn - R_MESH1;
                          if (p+1 < pd_glob[0]->Num_Dim) {
                            offset_mesh2 = 
			      get_nodal_unknown_offset(nv, eqn+1, vd->MatID, 0,
						       &vd_retn2);
                            idum = find_first_opening(BC_Unk_List[inode][offset_mesh2],
                                                      MAX_NODAL_BCS);
                            if (Debug_Flag > 0) {
                              fprintf(stderr, "    MOVING BC %s at node %d because of %s\n", 
                                      BC_Types[ibc2].desc->name2, 
                                      inode+1, BC_Types[ibc1].desc->name2);
                              BC_Unk_List[inode][offset_mesh2][idum] = ibc2;
                            }
                          } else {
                            if (Debug_Flag > 0) {
                              fprintf(stderr, 
                                      "    REMOVING BC %s at node %d because of %s and next coord not available\n", 
                                      BC_Types[ibc2].desc->name2, 
                                      inode+1, BC_Types[ibc1].desc->name2);
                            }
                          }
                        } else if (BC_Types[ibc2].desc->rotate == R_MOMENTUM1) {
                          p = eqn - R_MOMENTUM1;
                          if (p+1 < pd_glob[0]->Num_Dim) {
                            offset_mom2 = get_nodal_unknown_offset(nv, eqn+1, vd->MatID, 
                                                                   0, &vd_retn2);
                            idum = find_first_opening(BC_Unk_List[inode][offset_mom2],
                                                      MAX_NODAL_BCS);
                            if (Debug_Flag > 0) {
                              fprintf(stderr,
                                      "    MOVING BC %s at node %d to next coord because of %s\n", 
                                      BC_Types[ibc2].desc->name2, 
                                      inode+1, BC_Types[ibc1].desc->name2);
                            }
                            BC_Unk_List[inode][offset_mom2][idum] = ibc2;
                          } else {
                            if (Debug_Flag > 0) {
                              fprintf(stderr, 
                                      "    REMOVING BC %s at node %d because of %s and next coord not avail\n", 
                                      BC_Types[ibc2].desc->name2, 
                                      inode+1, BC_Types[ibc1].desc->name2);
                            }
                          }
                        } else {
                          WH(-1, "Missing type of rotation in BC conflict resolution");
                        }
                      }
                      /*
                       * Special Kluge
                       */
                      if (BC_Types[ibc1].desc->BC_Name ==  CONT_TANG_VEL_BC &&
                          BC_Types[ibc2].desc->BC_Name ==  SDC_STEFANFLOW_BC ) {
                        if (i_delete == 1) i_delete = -1;
                      }

                      /* remove this BC from list */
                      if (action == 1 && i_delete == 1) {
	
                        if (Debug_Flag > 0) {
                          fprintf(stderr, 
                                  "    DISCARDING %s at node %d from dup_list because of %s\n", 
                                  BC_Types[ibc2].desc->name2, 
                                  inode+1, BC_Types[ibc1].desc->name2);
                        }
                        BC_Unk_List[inode][offset][k] = -1;
                      }

                      if (i_delete == -1) {
                        if (Debug_Flag > 0) {
                          fprintf(stderr, 
                                  "    DISCARDING %s at node %d from dup_list because of %s\n", 
                                  BC_Types[ibc1].desc->name2, 
                                  inode+1, BC_Types[ibc2].desc->name2);
                        }
                        BC_Unk_List[inode][offset][j] = -1;
                      }
                    }		    
                  }	/* for (k=j+1; k<idup; k++) */
                } /* if (BC_Types[ibc1].desc->method == DIRICHLET) */
              } /* if (ibc1 != -1) */
            } /*  for (j=0; j<idup-1; j++)*/
            compress_bc_list(BC_Unk_List[inode][offset]);
            idup = find_first_opening(BC_Unk_List[inode][offset], MAX_NODAL_BCS);
            /*****************************************************************************/
            /* third give precedence to rotated strong (integrated or collocated)
             * boundary conditions 
             * - they swamp out weak boundary conditions
             * - Normal conditions have precedence over tangential */

            /* at this point, the only BC's left to conflict with are weak
             * conditions, special conditions, 
             * and other strong rotated conditions 
             * - if conflicting with a weak condition, dump the weak condition
             * - if conflicting with a special condition, leave it in there
             * - if conflicting with another rotated condition, 
             *   - check for collinearity, 
             *   - then choose normal over tangent, 
             *   - then then choose 1st in input file */

            for (j = 0; j < idup; j++) {
              ibc1 = BC_Unk_List[inode][offset][j];
              if (ibc1 != -1) {
                bct1 = BC_Types[ibc1].BC_Name;
                if (BC_Types[ibc1].desc->rotate != NO_ROT &&
                    (BC_Types[ibc1].desc->method == COLLOCATE_EDGE 
                     || BC_Types[ibc1].desc->method == SPECIAL)) {
                  /* 
                   * make sure that CONTACT ANGLE conditions end up 
                   * on the correct equation - it has 
                   * precedence over other rotated conditions
                   */
                  for (k = 0; k < idup; k++) {
                    /* flag indicating that this bc should be deleted from list */
                    i_delete = 1; i_dont_move = 0;
                    ibc2 = BC_Unk_List[inode][offset][k];
                    if ((ibc2 != ibc1) && (ibc2 != -1)) {
                      bct2 = BC_Types[ibc2].BC_Name;

                      if ( (BC_Types[ibc1].desc->BC_Name == SHEET_ENDSLOPE_BC) &&
                           (BC_Types[ibc2].desc->BC_Name == TENSION_SHEET_BC) ) 
                        {
                          i_dont_move = 1;  /* This is a kludge.  I'm hoping to do better tab*/
                          i_delete = 0;
                        }

                      /* move rotated conditions */
                      if ( ( i_dont_move != 1) && BC_Types[ibc2].desc->rotate != NO_ROT) {
                        /* move this condition ibc2 to the next coordinate 
                         * direction, if available 
                         */
                        if (BC_Types[ibc2].desc->rotate == R_MESH1) {
                          p = eqn - R_MESH1;
                          if (p+1 < pd_glob[0]->Num_Dim) {
                            offset_mesh2 = 
			      get_nodal_unknown_offset(nv, eqn+1, vd->MatID, 0,
						       &vd_retn2);
                            idum = find_first_opening(BC_Unk_List[inode][offset_mesh2],
                                                      MAX_NODAL_BCS);
                            if (Debug_Flag > 0) {
                              fprintf(stderr, 
                                      "    MOVING BC %s at node %d because of %s\n", 
                                      BC_Types[ibc2].desc->name2, 
                                      inode+1, BC_Types[ibc1].desc->name2);
                            }
                            BC_Unk_List[inode][offset_mesh2][idum] = ibc2;
                          } else {
                            if (Debug_Flag > 0) {
                              fprintf(stderr, 
                                      "    REMOVING BC %s at node %d because of %s\n", 
                                      BC_Types[ibc2].desc->name2, 
                                      inode+1, BC_Types[ibc1].desc->name2);
                            }
                          }
                        } else if (BC_Types[ibc2].desc->rotate == R_MOMENTUM1) {
                          p = eqn - R_MOMENTUM1;
                          if (p+1 < pd_glob[0]->Num_Dim) {
                            offset_mom2 = get_nodal_unknown_offset(nv, eqn+1, vd->MatID,
                                                                   0, &vd_retn2);
                            idum = find_first_opening(BC_Unk_List[inode][offset_mom2],
                                                      MAX_NODAL_BCS);
                            if (Debug_Flag > 0) {
                              fprintf(stderr, 
                                      "    MOVING BC %s at node %d because of %s\n", 
                                      BC_Types[ibc2].desc->name2, 
                                      inode+1, BC_Types[ibc1].desc->name2);
                            }
                            BC_Unk_List[inode][offset_mom2][idum] = ibc2;
                          } else {
                            if (Debug_Flag > 0) {
                              fprintf(stderr, 
                                      "    REMOVING BC %s at node %d because of %s\n", 
                                      BC_Types[ibc2].desc->name2, 
                                      inode+1, BC_Types[ibc1].desc->name2);
                            }
                          }
                        } else {
                          WH(-1, "Missing type of rotation in BC conflict resolution");
                        }
                      }

                      /* remove this BC from list */
                      if (i_delete == 1) {
                        if (Debug_Flag > 0) {
                          fprintf(stderr, "    DISCARDING %s at node %d from dup_list because of %s\n", 
                                  BC_Types[ibc2].desc->name2, 
                                  inode+1, BC_Types[ibc1].desc->name2);
                        }
                        BC_Unk_List[inode][offset][k] = -1;
                      }
                    }
                  }	/* for (k=j+1; k<idup; k++) */

                } else if (BC_Types[ibc1].desc->rotate != NO_ROT &&
                           (BC_Types[ibc1].desc->method == COLLOCATE_SURF || 
                            BC_Types[ibc1].desc->method == STRONG_INT_SURF)) {

                  for (k = 0; k < idup; k++) {
                    ibc2 = BC_Unk_List[inode][offset][k];
                    if ((ibc2 != ibc1) && (ibc2 != -1) && (ibc1 != -1)) {
                      bct2 = BC_Types[ibc2].BC_Name;
                      /* flag indicating whether this bc should be deleted from list */
                      /*
                       * HKM -> Why is this turned off ?
                       *        It is only turned on, if the second boundary condition
                       *        is a weakly integrated bc. Two strongly integrated bc's
                       *        are allowed to be applied to the same unknown.
                       */
                      i_delete = 0;

                      action = 1;
                      /* need to deal with GD conditions here - two GD conditions on the same SS and
                       * applied to the same equation do not conflict with each other
                       */
                      if ((bct1 <= GD_TIME_BC && bct1 >= GD_CONST_BC) && 
                          (bct2 <= GD_TIME_BC && bct2 >= GD_CONST_BC)) {
                        if (BC_Types[ibc1].BC_ID == BC_Types[ibc2].BC_ID) {
                          action = -1;
                          if (Debug_Flag > 0) {
                            fprintf(stderr, "    MULTIPLE GD %s at node %d\n", 
                                    BC_Types[ibc2].desc->name2, inode+1);
                          }
                        }
                      }

                      if (BC_Types[ibc2].desc->method == WEAK_INT_SURF ||
                          BC_Types[ibc2].desc->method == WEAK_SHARP_INT ) {
                        i_delete = 1;
                      } else if (BC_Types[ibc1].desc->method == COLLOCATE_EDGE || 
                                 BC_Types[ibc2].desc->method == SPECIAL) {
                        /* put additional cross check for bct1 being
                         * on a mesh equation.  E.G. Currently goes in
                         * here if bct1 is a qside 
                         */
                        if (bct2 == CA_BC ||
                            bct2 == CA_MOMENTUM_BC ||
                            bct2 == MOVING_CA_BC ||
                            bct2 == VELO_THETA_TPL_BC ||
                            bct2 == CA_OR_FIX_BC ||
                            bct2 == SHEET_ENDSLOPE_BC ) {
                          /* put additional check for equation type MESH */
                          WH(-1, "Illegal action possible for CA_BC w/o MESH");
                          i_delete = 0;
                          /* 
                           * make sure this ends up on the same equation as kinematic or capillary 
                           */
                          sprintf(Err_Msg, 
                                  "Contact angle at node %d from %s and %s?",
                                  inode+1, BC_Types[ibc1].desc->name1,
                                  BC_Types[ibc2].desc->name1);
                          WH(-1, Err_Msg);
                        } else {
                          if (Debug_Flag > 0) {
                            fprintf(stderr,
                                    "    REMOVING SPECIAL BC %s at node %d from dup_list because of %s\n",
                                    BC_Types[ibc2].desc->name2, 
                                    inode+1, BC_Types[ibc1].desc->name2);
                          }
                          i_delete = 1;
                        }
                      } else {

                        /* Only way to get here is if multiple rotated 
                         * conditions exist for this node
                         */
                        det = 0.;
                        if (bct1 == bct2 && action == 1) { /* these conditions are the same type */
                          /* Check for collinear planes, splines, geom, or GD
                           */
                          switch (bct2) {
                          case PLANE_BC:
			case ROLL_FLUID_BC:
                          case GEOM_BC:	/* aka SPLINE note that the relation here may depend on 
					 * the functional form of the geometry  */
                          case SLOPEX_BC:
                          case SLOPEY_BC:
                          case SLOPEZ_BC:
                          case SLOPE_BC:
			    a1 = BC_Types[ibc1].BC_Data_Float[0] ;
			    b1 = BC_Types[ibc1].BC_Data_Float[1] ;			
			    c1 = BC_Types[ibc1].BC_Data_Float[2] ;			
			    a2 = BC_Types[ibc2].BC_Data_Float[0] ;
			    b2 = BC_Types[ibc2].BC_Data_Float[1] ;
			    c2 = BC_Types[ibc2].BC_Data_Float[2] ;
			    /* determine if the lines are linearly independent */
			    det1 = a1*b2-a2*b1;
			    det2 = a1*c2-a2*c1;
			    det3 = c1*b2-c2*b1;
			    det = fabs(det1) + fabs(det2) + fabs(det3);
			    break;
			case FILLET_BC:
			case DOUBLE_RAD_BC:
			    a1 = BC_Types[ibc1].BC_Data_Float[0] ;
			    a2 = BC_Types[ibc1].BC_Data_Float[1] ;			
			    a3 = BC_Types[ibc1].BC_Data_Float[2] ;			
			    a4 = BC_Types[ibc1].BC_Data_Float[3] ;			
			    a5 = BC_Types[ibc1].BC_Data_Float[4] ;			
			    b1 = BC_Types[ibc2].BC_Data_Float[0] ;
			    b2 = BC_Types[ibc2].BC_Data_Float[1] ;
			    b3 = BC_Types[ibc2].BC_Data_Float[2] ;
			    b4 = BC_Types[ibc2].BC_Data_Float[3] ;
			    b5 = BC_Types[ibc2].BC_Data_Float[4] ;
			      det = fabs(a1 - b1) +  fabs(a2 - b2) +  fabs(a3 - b3)
                                       + fabs(a4-b4) + fabs(a5 - b5);
			    break;
                          case GD_CONST_BC:
			    if ( (BC_Types[ibc1].BC_Data_Int[2] == BC_Types[ibc2].BC_Data_Int[2]) &&
				 (BC_Types[ibc1].BC_Data_Int[1] == BC_Types[ibc2].BC_Data_Int[1]) ) {
			      /* the functional form of the geometry  */
			      a1 = BC_Types[ibc1].BC_Data_Float[0] ;
			      a2 = BC_Types[ibc2].BC_Data_Float[0] ;
			      /* determine if the lines are linearly independent */
			      det = a1-a2;
			    }
			    else det = 1;
			    break;
                          case GD_TIME_BC:
                          case GD_LINEAR_BC:
			    if ( (BC_Types[ibc1].BC_Data_Int[2] == BC_Types[ibc2].BC_Data_Int[2]) &&
				 (BC_Types[ibc1].BC_Data_Int[1] == BC_Types[ibc2].BC_Data_Int[1]) ) {
			      a1 = BC_Types[ibc1].BC_Data_Float[0] ;
			      b1 = BC_Types[ibc2].BC_Data_Float[0] ;		       
			      a2 = BC_Types[ibc1].BC_Data_Float[1] ;
			      b2 = BC_Types[ibc2].BC_Data_Float[1] ;
			      /* determine if the lines are linearly independent */
			      det = fabs(a1 - b1) +  fabs(a2 - b2);
			    }
			    else det = 1;
			    break;
                          case GD_CIRC_BC:
                          case GD_PARAB_BC:
			    if ( (BC_Types[ibc1].BC_Data_Int[2] == BC_Types[ibc2].BC_Data_Int[2]) &&
				 (BC_Types[ibc1].BC_Data_Int[1] == BC_Types[ibc2].BC_Data_Int[1]) ) {
			      a1 = BC_Types[ibc1].BC_Data_Float[0] ;
			      b1 = BC_Types[ibc2].BC_Data_Float[0] ;		
			      a2 = BC_Types[ibc1].BC_Data_Float[1] ;
			      b2 = BC_Types[ibc2].BC_Data_Float[1] ;
			      a3 = BC_Types[ibc1].BC_Data_Float[2] ;
			      b3 = BC_Types[ibc2].BC_Data_Float[2] ;
			      det = fabs(a1 - b1) +  fabs(a2 - b2) +  fabs(a3 - b3);
			    }
			    break;
                          case GD_PARAB_OFFSET_BC:
			    if ( (BC_Types[ibc1].BC_Data_Int[2] == BC_Types[ibc2].BC_Data_Int[2]) &&
				 (BC_Types[ibc1].BC_Data_Int[1] == BC_Types[ibc2].BC_Data_Int[1]) ) {
			      a1 = BC_Types[ibc1].BC_Data_Float[0] ;
			      b1 = BC_Types[ibc2].BC_Data_Float[0] ;
			      a2 = BC_Types[ibc1].BC_Data_Float[1] ;
			      b2 = BC_Types[ibc2].BC_Data_Float[1] ;
			      a3 = BC_Types[ibc1].BC_Data_Float[2] ;
			      b3 = BC_Types[ibc2].BC_Data_Float[2] ;
			      a4 = BC_Types[ibc1].BC_Data_Float[3] ;
			      b4 = BC_Types[ibc2].BC_Data_Float[3] ;
			      det = fabs(a1 - b1) +  fabs(a2 - b2) +  fabs(a3 - b3) + fabs(a4 - b4);
			    }
			    break;				
                          case GD_POLYN_BC:
			    if ( (BC_Types[ibc1].BC_Data_Int[2] == BC_Types[ibc2].BC_Data_Int[2]) &&
				 (BC_Types[ibc1].BC_Data_Int[1] == BC_Types[ibc2].BC_Data_Int[1]) ) {
			      a1 = BC_Types[ibc1].BC_Data_Float[0] ;
			      b1 = BC_Types[ibc2].BC_Data_Float[0] ;		
			      a2 = BC_Types[ibc1].BC_Data_Float[1] ;
			      b2 = BC_Types[ibc2].BC_Data_Float[1] ;
			      a3 = BC_Types[ibc1].BC_Data_Float[2] ;
			      b3 = BC_Types[ibc2].BC_Data_Float[2] ;
			      a4 = BC_Types[ibc1].BC_Data_Float[3] ;
			      b4 = BC_Types[ibc2].BC_Data_Float[3] ;
			      a5 = BC_Types[ibc1].BC_Data_Float[4] ;
			      b5 = BC_Types[ibc2].BC_Data_Float[4] ;
			      a6 = BC_Types[ibc1].BC_Data_Float[5] ;
			      b6 = BC_Types[ibc2].BC_Data_Float[5] ;
			      a7 = BC_Types[ibc1].BC_Data_Float[6] ;
			      b7 = BC_Types[ibc2].BC_Data_Float[6] ;
			      det = fabs(a1 - b1) +  fabs(a2 - b2) +  fabs(a3 - b3)
                                + fabs(a4 - b4) +  fabs(a5 - b5) +  fabs(a6 - b6) +  fabs(a7 - b7);
			    }
			    break;
                          case VNORM_LEAK_BC:
			    a1 = BC_Types[ibc1].BC_Data_Float[0] ;
			    b1 = BC_Types[ibc2].BC_Data_Float[0] ;
			    a2 = BC_Types[ibc1].BC_Data_Float[1] ;
			    b2 = BC_Types[ibc2].BC_Data_Float[1] ;
			    /* determine if the lines are linearly independent */
			    det = fabs(a1 - b1) +  fabs(a2 - b2);
			    break;
                          case DISTNG_BC:
                          case VELO_NORMAL_BC:
                          case VELO_NORMAL_DISC_BC:
                          case VELO_NORMAL_LS_BC:
			    /* the functional form of the geometry  */
			    a1 = BC_Types[ibc1].BC_Data_Float[0] ;
			    a2 = BC_Types[ibc2].BC_Data_Float[0] ;
			    /* determine if the lines are linearly independent */
			    det = a1-a2;
			    /* Actually, also determine if the user is trying to
			       override the code and get one velo_normal shifted to
			       the next available equation, as you would want to do
			       if the lines are at least a bit anti-parallel */
			    if( BC_Types[ibc1].BC_Data_Int[0] == 1 ||
				BC_Types[ibc2].BC_Data_Int[0] == 1) det = 1.;
			
			    break;
                          case VELO_TANGENT_BC:
                          case VELO_TANGENT_USER_BC:
                          case VELO_STREAMING_BC:
			    a1 = BC_Types[ibc1].BC_Data_Float[0] ;
			    b1 = BC_Types[ibc2].BC_Data_Float[0] ;		
			    a2 = BC_Types[ibc1].BC_Data_Float[1] ;
			    b2 = BC_Types[ibc2].BC_Data_Float[1] ;
			    a3 = BC_Types[ibc1].BC_Data_Float[2] ;
			    b3 = BC_Types[ibc2].BC_Data_Float[2] ;
			    det = fabs(a1 - b1) +  fabs(a2 - b2) +  fabs(a3 - b3);
			    det += fabs( (double) BC_Types[ibc1].BC_Data_Int[0] + 
					 (double) BC_Types[ibc2].BC_Data_Int[0] );
			    break;
                          case VELO_SLIP_BC:
                          case VELO_SLIP_ROT_BC:
			case VELO_SLIP_FLUID_BC:
			case VELO_SLIP_ROT_FLUID_BC:
			case AIR_FILM_BC:
			case AIR_FILM_ROT_BC:
			    a1 = BC_Types[ibc1].BC_Data_Float[0] ;
			    b1 = BC_Types[ibc2].BC_Data_Float[0] ;		
			    a2 = BC_Types[ibc1].BC_Data_Float[1] ;
			    b2 = BC_Types[ibc2].BC_Data_Float[1] ;
			    a3 = BC_Types[ibc1].BC_Data_Float[2] ;
			    b3 = BC_Types[ibc2].BC_Data_Float[2] ;
			    a4 = BC_Types[ibc1].BC_Data_Float[3] ;
			    b4 = BC_Types[ibc2].BC_Data_Float[3] ;
			    a5 = BC_Types[ibc1].BC_Data_Float[4] ;
			    b5 = BC_Types[ibc2].BC_Data_Float[4] ;
			    det = fabs(a1-b1) + fabs(a2-b2) + fabs(a3-b3)+fabs(a4-b4)+fabs(a5-b5);
			    det += fabs( (double) BC_Types[ibc1].BC_Data_Int[0] + 
					 (double) BC_Types[ibc2].BC_Data_Int[0] );
			    break;
			    /*
			     * Conditions which always get applied along both 
			     * intersecting boundaries at the same node
			     */
                          case KIN_LEAK_BC:
                          case KINEMATIC_BC:
                          case KINEMATIC_DISC_BC:
                          case KINEMATIC_COLLOC_BC:
                          case KINEMATIC_PETROV_BC:
			    det = 0.;	/* default to zero */
			    break;
                          default:
			    WH(-1, "Intersection between unreconcilable rotated conditions");
			    det = 1.;	/* default to nonzero value */
			    break;
                          }
		    
                          if (fabs(det) < 1.0e-10) {
                            if (Debug_Flag > 0) {
                              fprintf(stderr, "    COLINEAR SIDES %s at node %d \n", 
                                      BC_Types[ibc1].desc->name2, inode+1);
                            }
                            i_delete = 1;
                            /* colinear strong conditions don't get deleted */
                            if (BC_Types[ibc2].desc->method == STRONG_INT_SURF) {
                              i_delete = 2;
                            }
                          }
                        }

                        if (i_delete == 0 && action == 1) {
                          /* both ibc1 and ibc2 are rotated strong conditions that
                           * are not colinear 
                           *  - move ibc1 to the next available direction 
                           *  - unless ibc1 is a tangent condition and ibc2 
                           *    is a normal condition
                           *
                           *  HKM -> The algorithm here could stand an improvement.
                           *         One idea is to develop a directional cosine
                           *         approximation for each boundary condition
                           *         with respect to the axes. The boundary condition
                           *         with the largest absolute value of the
                           *         directional cosine wrt the x axis would be
                           *         assigned to the x coordinate, etc. I have
                           *         seen cases (two PLANE_BC conditions) where
                           *         this logic would prevent a zero being put
                           *         on the diagonal. Note, no rotation is 
                           *         performed when bc's are rotated bc's, so
                           *         the x, y, z directional cosine issue is the
                           *         pertinent one. 
                           */
                          imove = ibc2;
                          istay = ibc1;
                          i_delete = 1;
                          /*
                           * move integrated conditions in favor of collocated conditions 
                           */
                          if (BC_Types[ibc1].desc->method == COLLOCATE_SURF &&
                              BC_Types[ibc2].desc->method == STRONG_INT_SURF) {
                            /* do nothing */
                          }
                          else if (BC_Types[ibc2].desc->method == COLLOCATE_SURF && 
                                   BC_Types[ibc1].desc->method == STRONG_INT_SURF) {
                            imove = ibc1;
                            istay = ibc2;
                            i_delete = -1;
                          }
                          /*
                           * use trick to determine if ibc1 is tangent and ibc2 is normal 
                           * If the first boundary condition is a tangent condition, move it
                           * and keep the second boundary condition
                           */
                          else if (BC_Types[ibc1].desc->equation > BC_Types[ibc2].desc->equation) {
                            imove = ibc1;
                            istay = ibc2;
                            i_delete = -1; /* indicates that first bc should be deleted */
                          }
                          if (BC_Types[imove].desc->rotate == R_MESH1){
                            p = eqn - R_MESH1;
                            if (p+1 < pd_glob[0]->Num_Dim) {
                              offset_mesh2 = 
				get_nodal_unknown_offset(nv, eqn+1, vd->MatID, 0,
							 &vd_retn2);
                              idum = find_first_opening(BC_Unk_List[inode][offset_mesh2],
                                                        MAX_NODAL_BCS);
                              if (Debug_Flag > 0) {
                                fprintf(stderr, "    MOVING BC %s at node %d because of %s\n", 
                                        BC_Types[imove].desc->name2, 
                                        inode+1, BC_Types[istay].desc->name2);
                              }
                              BC_Unk_List[inode][offset_mesh2][idum] = imove;
                            } else {
                              if (Debug_Flag > 0) {
                                fprintf(stderr, "    REMOVING BC %s at node %d because of %s\n", 
                                        BC_Types[imove].desc->name2, 
                                        inode+1, BC_Types[istay].desc->name2);
                              }
                            }
                          } else if (BC_Types[imove].desc->rotate == R_MOMENTUM1) {
                            p = eqn - R_MOMENTUM1;
                            if (p+1 < pd_glob[0]->Num_Dim) {
                              offset_mom2 = 
				get_nodal_unknown_offset(nv, eqn+1, vd->MatID, 0,
							 &vd_retn2);
                              idum = find_first_opening(BC_Unk_List[inode][offset_mom2],
                                                        MAX_NODAL_BCS);
                              if (Debug_Flag > 0) {
                                fprintf(stderr, "    MOVING BC %s at node %d because of %s\n", 
                                        BC_Types[imove].desc->name2, 
                                        inode+1, BC_Types[istay].desc->name2);
                              }
                              BC_Unk_List[inode][offset_mom2][idum] = imove;
                            } else {
                              if (Debug_Flag > 0) {
                                fprintf(stderr, "    REMOVING BC %s at node %d because of %s\n", 
                                        BC_Types[imove].desc->name2, 
                                        inode+1, BC_Types[istay].desc->name2);
                              }
                            }
                          } else {
                            WH(-1, "Missing type of rotation in BC conflict resolution");
                          }
                        }

                        /* remove this BC from list */
                        if (action == 1 && i_delete == 1) {
                          if (Debug_Flag > 0) {
                            fprintf(stderr,
                                    "    DISCARDING %s at node %d from dup_list because of %s\n", 
                                    BC_Types[ibc2].desc->name2, 
                                    inode+1, BC_Types[ibc1].desc->name2);
                          }
                          BC_Unk_List[inode][offset][k] = -1;
                        }

                        if (action == 1 && i_delete == -1) {
                          if (Debug_Flag > 0) {
                            fprintf(stderr,
                                    "    DISCARDING %s at node %d from dup_list because of %s\n", 
                                    BC_Types[ibc1].desc->name2, 
                                    inode+1, BC_Types[ibc2].desc->name2);
                          }
                          BC_Unk_List[inode][offset][j] = -1;
                        }

                      }
                    }		    
                  }	/* for (k=j+1; k<idup; k++) */
                } /* if (BC_Types[ibc1].desc->method == DIRICHLET) */
              } /* if (ibc1 != -1) */
            } /*  for (j=0; j<idup-1; j++)*/

            /* 
             * WEAK conditions remaining in dup list apply as they are, 
             * multiple weak conditions are ok
             */

          } /* if idup>1 */
        }	/*  for (eqn=0; eqn<MAX_VARIABLE_TYPES + MAX_CONC; eqn++)  */

        /*
         * Make sure the bc list is compressed
         */
        for (offset = 0; offset < nv->Num_Unknowns; offset++) {
          compress_bc_list(BC_Unk_List[inode][offset]);
        }
      } /* if (BC_list[inode] != NULL) */
    } /*for (inode=0; inode<num_total_nodes; inode++)  */

    /***************************************************************************/
    /* CLEAN UP                                                                */
    /***************************************************************************/

    /* now print out BC_dup list for 2-D problems */
    /* also check for potentially conflicting combinations */

    /*
     * Conditionally divert out from stderr to a file if it looks like
     * it might be voluminous...
     * MP jobs sync decision to divert by doing a gsum_Int().
     */

    /* Code output cleanup crew says always divert BC dup spew to file from here on out.
     * If you care enough about it to understand it you can take two seconds and open the file 
     *
     ****old code  *****
     divert = 0;  
     if (Num_Proc > 1) {
     divert = 1;
     } else {
     if (exo->num_node_sets > DUP_THRESHHOLD_NODESETS) divert = 1;
     if (exo->num_side_sets > DUP_THRESHHOLD_SIDESETS) divert = 1;
     if (Num_BC > DUP_THRESHHOLD_NUMBCS)               divert = 1;
     }
    */

    divert = 1;

    divert = gsum_Int(divert);

    if (Debug_Flag > 0 || pd_glob[0]->Num_Dim < 3 )  {
      print_sync_start(TRUE);
      if (divert) {
        strcpy(bc_divert_fn, DUP_THRESHHOLD_FILENAME);
        if (ProcID == 0) {
          bc_dup_out = fopen(bc_divert_fn, "w");
        } else {
          bc_dup_out = fopen(bc_divert_fn, "a");	
        }
        if (BC_dup_ptr[pg->imtrx] > 0) {
          fprintf(bc_dup_out,"BC Duplication list for Processor %d:\n",
                  ProcID);
          fprintf(bc_dup_out,"-------- --------- ----------- ----\n");
        }
      } else {
        bc_dup_out = stderr;
      }
    

      for (i = 0; i < BC_dup_ptr[pg->imtrx]; i++) {
        inode = BC_dup_nodes[pg->imtrx][i];
        if (inode < (dpi->num_internal_nodes + dpi->num_boundary_nodes)) {
          node = Nodes[inode];
          nv = node->Nodal_Vars_Info[pg->imtrx];
          fprintf(bc_dup_out, "\nAt global node %d",
                  node->Global_Node_Num + 1);
          fprintf(bc_dup_out, " (proc node %d)", inode + 1);
          if (pd_glob[0]->Num_Dim < 3) {
            fprintf(bc_dup_out, " (%g, %g)\n", 
                    exo->x_coord[inode],  exo->y_coord[inode]);
          } else {  
            fprintf(bc_dup_out, " (%g, %g, %g)\n", 
                    exo->x_coord[inode],  exo->y_coord[inode], 
                    exo->z_coord[inode]);
          }
          for (offset = 0; offset < nv->Num_Unknowns; offset++) {
            get_nv_vd_from_offset(nv, offset, &vd, &idof);
            eqn = vd->Variable_Type;
            if ((bickel = BC_Unk_List[inode][offset][0]) != -1) {
              if (eqn != R_MASS) {
                fprintf(bc_dup_out, " %-15s", EQ_Name[eqn].name1);
                if (vd->MatID != -1) {
                  fprintf(bc_dup_out, " (MatID=%d)", vd->MatID);
                }
                fprintf(bc_dup_out, " gets ");
              } else {
                fprintf(bc_dup_out, " %-15s",
                        EQ_Name[MAX_VARIABLE_TYPES + vd->Subvar_Index].name1);
                if (vd->MatID != -1) {
                  fprintf(bc_dup_out, " (MatID=%d)", vd->MatID);
                }
                fprintf(bc_dup_out, " gets ");
              }
              idup = 0;
              while ((bickel = BC_Unk_List[inode][offset][idup]) != -1) {
                if (idup == 0) {
                  fprintf(bc_dup_out, "  %-15s (%d) from %s %d\n", 
                          BC_Types[bickel].desc->name1, bickel,
                          BC_Types[bickel].Set_Type,
                          BC_Types[bickel].BC_ID);
                } else {
                  fprintf(bc_dup_out, "\t\t\t%-15s (%d) from %s %d\n", 
                          BC_Types[bickel].desc->name1, bickel,
                          BC_Types[bickel].Set_Type,
                          BC_Types[bickel].BC_ID);
                }
                idup++;
              }
            }
          }
        }
      }
      fprintf(bc_dup_out, "\n");

      if (divert) {
        fflush(bc_dup_out);
        fclose(bc_dup_out);
      }
      print_sync_end(TRUE);
    }

    /************************************************************************/
    /* ROTATION - determine whether mesh or momentum equations are rotated  
     * for each node and determine wrt which SS they are rotated
     */
    /************************************************************************/
    /* make a list of momentum and mesh nodes that will have rotated conditions
     * and the SS#'s of the sides on which they should be rotated
     */
    num_mesh_rotate[pg->imtrx]  = 0;
    num_mom_rotate[pg->imtrx]   = 0;

    if( num_bc_nodes > 0 ) {
      mesh_rotate_node[pg->imtrx] = alloc_int_1(num_bc_nodes, INT_NOINIT);
      mesh_rotate_ss[pg->imtrx]   = alloc_int_1(num_bc_nodes, INT_NOINIT);
      mom_rotate_node[pg->imtrx]  = alloc_int_1(num_bc_nodes, INT_NOINIT);
      mom_rotate_ss[pg->imtrx]    = alloc_int_1(num_bc_nodes, INT_NOINIT);
    }

    /* loop through nodes and remove nodes that don't have conflicts 
     * count nodes with conflicts and build an array to store them 
     */
    for (inode = 0; inode < num_total_nodes; inode++) {
      node = Nodes[inode];
      nv =  node->Nodal_Vars_Info[pg->imtrx];
      if (BC_Unk_List[inode] != NULL) {
        for (ivar = 0; ivar < 2; ivar++) {
          /* if ivar=0 checking the MESH equations
           * if ivar=1 checking the MOMENTUM equations
           */
          /*
           * Have to deal with the possibility of multiple momentum equations at
           * this node
           */
          if (ivar == 0) eqn = R_MESH1;
          else           eqn = R_MOMENTUM1;
          ndofs = get_nv_ndofs(nv, eqn);
          for (idof = 0; idof < ndofs; idof++) {
            if ( idof >= nv->Num_Var_Desc_Per_Type[eqn] ) continue; /*check needed for XFEM*/
            offset = get_nv_offset_idof(nv, eqn, idof, 0, &vd_retn1);
            /* LOOK at MESH equations and BC's applied to those equations 
             * COUNT the occurences of:
             */
            count_BC = 0;
            count_DC = 0;
            count_rotate = 0;
            count_weak = 0;
            count_strong = 0;
            count_coll = 0;
            count_special = 0;
            ss_rot = -1;
            two_ss_rot = -1;
            dir_rot = -1;
            for (jvar = 0; jvar < pd_glob[0]->Num_Dim; jvar++) {
              var = eqn + jvar;
              offset_p = get_nv_offset_idof(nv, var, idof, 0, &vd_retn1);
              idup = 0;
              while (BC_Unk_List[inode][offset_p] != NULL && 
                     BC_Unk_List[inode][offset_p][idup] != -1) {
                count_BC++;
                ibc = BC_Unk_List[inode][offset_p][idup];
                bc_desc = BC_Types[ibc].desc;
                if (bc_desc->method == DIRICHLET)       count_DC++;
                if (bc_desc->method == WEAK_INT_SURF)   count_weak++;
                if (bc_desc->method == WEAK_SHARP_INT)   count_weak++;
                if (bc_desc->method == STRONG_INT_SURF) count_strong++;
                if (bc_desc->method == COLLOCATE_SURF)  count_coll++;
                if (bc_desc->method == SPECIAL)         count_special++;
                /* count independent rotated conditions */
                if (bc_desc->rotate != NO_ROT) {
                  if (count_rotate == 0 || jvar != dir_rot) {
                    if (ss_rot != -1 &&
                        ss_rot != BC_Types[ibc].BC_ID) {
                      two_ss_rot = 1;
                    }
                    ss_rot = BC_Types[ibc].BC_ID;
                    dir_rot = jvar;
                    count_rotate++;
                  } else {
                    /* current rotated BC must be an allowed duplicate of 
                     * the previous rotated BC because it is in the same direction, so don't 
                     * double-count it
                     */
                  }
                }
                idup++;
              }
            }
            if (Debug_Flag > 1) {
	      fprintf(stderr, "%d BC's at node %d on %d: DC=%d "
		      "WIC=%d SIC=%d PC=%d SP=%d ROT=%d\n",
		      count_BC,inode+1,ivar, count_DC,count_weak,
		      count_strong,count_coll,count_special,count_rotate);
            }
            /* evaluate if node needs rotation */
            i_rotate = 0;
            if (count_rotate > 0 &&
                count_rotate < pd_glob[0]->Num_Dim) {
              i_rotate = 1;
            }
            if (count_rotate >= pd_glob[0]->Num_Dim && ivar == 0 && 
                (pd_glob[0]->Num_Dim < 3 || Debug_Flag >= 2) ) {
	      log_msg("Don't rotate MESH at node %d because too many rotations\n",
		      inode + 1);
            }
            if (count_rotate >= pd_glob[0]->Num_Dim && ivar == 1 &&  
                (pd_glob[0]->Num_Dim < 3 || Debug_Flag >= 2) ) {
              log_msg("Don't rotate MOMENTUM at node %d because too many rotations\n",
                      inode + 1);
            }
            if (count_DC > 0) {
              i_rotate = 0;
              if (count_rotate > 0) {
                if (ivar == 0) {

                  log_msg("Don't rotate MESH at node %d because there is a Dirichlet condition\n",
                          inode + 1);
                } else {

                  log_msg("Don't rotate MOMENTUM at node %d because there is a Dirichlet condition\n",
                          inode + 1);
                }
              }
            }
            if (i_rotate && ivar == 0) {
              if (Debug_Flag > 0) {
                fprintf(stderr, "Rotate MESH for BC's at node %d\n",
                        inode+1);
              }
              if (num_mesh_rotate[pg->imtrx] == 0 ||
                  mesh_rotate_node[pg->imtrx][num_mesh_rotate[pg->imtrx]-1] != inode) {
                mesh_rotate_node[pg->imtrx][num_mesh_rotate[pg->imtrx]] = inode;
                mesh_rotate_ss[pg->imtrx][num_mesh_rotate[pg->imtrx]] = ss_rot;
                num_mesh_rotate[pg->imtrx]++;
              }
            }
            if (i_rotate && ivar == 1) {
              if (Debug_Flag > 0)
		fprintf(stderr, "Rotate MOMENTUM for BC's at node %d\n",
			inode+1);
              if (num_mom_rotate[pg->imtrx] == 0 ||
                  mom_rotate_node[pg->imtrx][num_mom_rotate[pg->imtrx]-1] != inode) {
                mom_rotate_node[pg->imtrx][num_mom_rotate[pg->imtrx]] = inode;
                mom_rotate_ss[pg->imtrx][num_mom_rotate[pg->imtrx]] = ss_rot;
                num_mom_rotate[pg->imtrx]++;
              }
            }
      
            /* process warnings for this node */
            if (count_DC > 0 && count_rotate > 0) {

              log_msg("WARNING: combining dirichlet and rotated conditions "
                      "at node %d may be dangerous\n", inode+1);
            }
            if (ivar == 1 && count_rotate > 1 && two_ss_rot == 1) {

              log_msg( "WARNING: applying 2 rotated momentum conditions on MOMENTUM at "
                       "junction node %d may be dangerous", inode+1); 
            }
          }
        } /* end of loop over ivar */
      }/* BC_Unk_List == NULL */
    }

    if ( num_bc_nodes > 0 ) {
      mesh_rotate_node[pg->imtrx] = (int *)realloc(mesh_rotate_node[pg->imtrx], (size_t) num_mesh_rotate[pg->imtrx] * sizeof(int));
      mesh_rotate_ss[pg->imtrx]   = (int *)realloc(mesh_rotate_ss[pg->imtrx], (size_t) num_mesh_rotate[pg->imtrx] * sizeof(int));
      mom_rotate_node[pg->imtrx]  = (int *)realloc(mom_rotate_node[pg->imtrx], (size_t) num_mom_rotate[pg->imtrx] * sizeof(int));
      mom_rotate_ss[pg->imtrx]    = (int *)realloc(mom_rotate_ss[pg->imtrx], (size_t) num_mom_rotate[pg->imtrx] * sizeof(int));
    }

    if( mesh_rotate_node[pg->imtrx] != NULL || mom_rotate_node[pg->imtrx] != NULL)
      {
	int ebi=0;
	
        while( !Use_2D_Rotation_Vectors && ebi < exo->num_elem_blocks )
          {
			Use_2D_Rotation_Vectors = ( exo->eb_elem_itype[ebi] == LINEAR_TRI || exo->eb_elem_itype[ebi] == QUAD_TRI);

            ebi++;	
          }
      }
	
    if (Debug_Flag > 0) {
      for (i = 0; i < num_mesh_rotate[pg->imtrx]; i++) {
        fprintf(stderr, "\tMESH node %d rotate on SS %d\n",
                mesh_rotate_node[pg->imtrx][i]+1,mesh_rotate_ss[pg->imtrx][i]);
      }
    }
    if (Debug_Flag > 0 && num_mom_rotate[pg->imtrx] > 0) {
      fprintf(stderr, "Rotate MOMENTUM for BC's at %d NODES \n",num_mom_rotate[pg->imtrx]);
    }
    if (Debug_Flag > 0) {
      for (i = 0; i < num_mom_rotate[pg->imtrx]; i++) {
        fprintf(stderr, "\tMOMENTUM node %d rotate on SS %d\n",
                mom_rotate_node[pg->imtrx][i]+1,mom_rotate_ss[pg->imtrx][i]);
      }
    }
  
    /* free memory associated with nodes that don't have duplicate BC's
     */
    for (inode = 0; inode < num_total_nodes; inode++) {
      if (BC_Unk_List[inode] != NULL) {
        if (in_list(inode, 0, BC_dup_ptr[pg->imtrx], BC_dup_nodes[pg->imtrx]) == -1) {
          free_bc_unk_list_node(BC_Unk_List + inode, inode);
        }
      }
    }
    /* remove array which held boundary condition information temporarily
     * for each node */
    safer_free((void **)&BC_Unk_List);
  }

  if (ProcID == 0)
    {
      for (int ibc = 0; ibc < Num_BC; ibc++)
        {
          if (matrix_used_BC[ibc] == 0)
            {
              WH_MANY(-1,
               "Boundary condition %d, %s, applied on NS %d, is never used\n",
                        ibc, BC_Types[ibc].desc->name1, BC_Types[ibc].BC_ID);
            }
        }
      free(matrix_used_BC);
    }


  pg->imtrx = 0;

  return;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void 
check_for_bc_conflicts3D(Exo_DB *exo, Dpi *dpi)

    /*************************************************************************
     *
     * check_for_bc_conflicts3D():
     *
     * this routine takes care of several tasks up front for the application
     * of boundary conditions.  Among other things, 
     * It looks for multiple bc's being applied to the same equation at the 
     * same node, and sets up criteria for deciding which BC has precedence
     * Last revised by Rich Cairncross 12/05/95
     *************************************************************************/
{
  int side_index, lni, len_ssid_list;
  int *ssid_list;
  int i, ibc, iss, num_bc_nodes;
  int  ibc1, ibc2, inode, bct1, bct2, k, ins;
  int ***BC_Unk_List, **NS_list, eqn, idup, p, q, j, num_total_nodes, dups;
  int num_rot_nodes, irc, eq = -1, node_ok, iptr;
#ifndef PARALLEL
  int bc_found;
#endif
  int dim, bct, bcSS, bcSS1, eprint;
  int standard_BC, ndup, j_DC, j_PC, j_SI, j_CSWI;
  int save_this_bc[MAX_SS_PER_NODE];
  int ndup1[DIM], save_this_bc1[DIM][MAX_SS_PER_NODE];
  int faceBC, vertexBC, edgeBC;
  FILE *ofbc = NULL;		/* output file stream handle for BC info */
  char ofbc_fn[MAX_FNL];	/* output file name for BC info */
  int matIndex, offset, retn_matIndex;
  int ielem, elem_block_index;
  int offset_mom2, idof;
  NODE_INFO_STRUCT *node;
  NODAL_VARS_STRUCT *nv;
  VARIABLE_DESCRIPTION_STRUCT *vd, *vd_retn2;
#ifndef PARALLEL
  int tilt_bct = -1;
  int tilt_bcSS = -1;
  int tilt_save = -1;
  int val = 0;
#endif
  strcpy(ofbc_fn, BC_3D_INFO_FILENAME); /* def in rf_bc_const.h */

  if (Unlimited_Output) {
    print_sync_start(TRUE);
    if (ProcID == 0) {
      ofbc = fopen(ofbc_fn, "w");
    } else {
      ofbc = fopen(ofbc_fn, "a");
    }
    if (ofbc == NULL) {
      sprintf(Err_Msg, 
	      "P_%d: cannot open file \"%s\" to write 3D BC info",
	      ProcID, ofbc_fn);
      EH(GOMA_ERROR, Err_Msg);
    }
    if (Num_Proc > 1) {
      fprintf(ofbc, "BC Duplication and Rotation list for Processor %d:\n",
	      ProcID);
      fprintf(ofbc, "----------- ------------ -------------- -------------\n");
    }
  }

  BC_dup_nodes = malloc(sizeof(int *) * (size_t) upd->Total_Num_Matrices);
  BC_dup_list = malloc(sizeof(int ***) * (size_t) upd->Total_Num_Matrices);
  BC_dup_ptr = malloc(sizeof(int) * (size_t) upd->Total_Num_Matrices);

  for (pg->imtrx = 0; pg->imtrx < upd->Total_Num_Matrices; pg->imtrx++) {

    /***************************************************************************/
    /*                              BLOCK 3                                    */
    /*      Start Checking for Duplicate Boundary Conditions at the same node  */
    /***************************************************************************/
    dim = pd_glob[0]->Num_Dim;

    num_total_nodes = ( dpi->num_internal_nodes + 
                        dpi->num_boundary_nodes +
                        dpi->num_external_nodes );

    /* create an array which will hold boundary condition information temporarily
     * for each node 
     */

    BC_Unk_List = (int ***) alloc_ptr_1(num_total_nodes);
    NS_list = (int **) alloc_ptr_1(num_total_nodes);
    for (i = 0; i < num_total_nodes; i++) {
      NS_list[i] = alloc_int_1(MAX_SS_PER_NODE, -1);
    }

    num_bc_nodes = 0;

    /***************************************************************************/
    /*
     * FIRST make a list of all the boundary conditions from the input deck
     *  at each node
     */
    /***************************************************************************/
    /*****************************************************************************/
    /*      BOUNDARY CONDITIONS SPECIFIED BY NODE SETS                           */
    /*  visit each node on all node sets and list the BC's applied there         */
    /*****************************************************************************/

    for (ibc = 0; ibc < Num_BC; ibc++) {

      /* First check to see if this boundary condition is a node set 	      */
      if (!strcmp(BC_Types[ibc].Set_Type, "NS")) {
        matIndex = -2;
        if (BC_Types[ibc].BC_EBID_Apply != -1) {
          matIndex = map_mat_index(BC_Types[ibc].BC_EBID_Apply);
        }
        /*
         *  Loop over the total number of node sets defined for the 
         *  current processor 
         */
        for (ins = 0; ins < exo->num_node_sets; ins++) {

          /* Check for a match between the ID of the current node set 
           * and the node set ID specified in the input file 
           * - continue if a match is found 
           */
          if (exo->ns_id[ins] == BC_Types[ibc].BC_ID) {

            /* Loop over all the global nodes in the current node set */
            for (i = 0; i < exo->ns_num_nodes[ins]; i++) {

              /* Get the ith local node in the current node set */
              inode = exo->ns_node_list[exo->ns_node_index[ins] + i];

              /*
               * Find the first available opening in the list of node set
               * id's applied at the current node. Then add the current
               *  node set id
               */
              iptr = find_first_opening(NS_list[inode], MAX_SS_PER_NODE);
              NS_list[inode][iptr] = exo->ns_id[ins];
              /*
               *  Allocate temp BC storage for this node and initialize
               *  the storage to -1
               */
              if (BC_Unk_List[inode] == NULL) {
                num_bc_nodes++;
                BC_Unk_List[inode] = alloc_bc_unk_list_node(inode);
              }

              /*
               * Fill in an entry into the list for each component of the
               * boundary condition
               */
              for (p = 0; p < BC_Types[ibc].desc->vector; p++) {
                offset = find_bc_unk_offset(BC_Types + ibc, matIndex, inode, p,
                                            &retn_matIndex, &vd);
                if (offset >= 0) {
                  add_to_node_unk_bc_list(BC_Unk_List[inode][offset], ibc);
                }
              }  
            } /* for (i = 0; i < exo->ns_num_nodes[ins]; i++) */
          } /* if (exo->ns_id[ins] == BC_Types[ibc].BC_ID) */
        } /*for (ins = 0; ins < Proc_Num_Node_Sets; ins++) */
      }  /* END if (!strcmp(BC_Types[ibc].Set_Type, "NS"))   */
    } /*for (ibc = 0; ibc < Num_BC; ibc++) */

    /*****************************************************************************/
    /*      BOUNDARY CONDITIONS SPECIFIED BY SIDE SETS                           */
    /*  visit each node on all side sets and list the BC's applied there         */
    /*  Also, make a list of all SS at each node                                 */
    /*****************************************************************************/

    SS_list = (int **) alloc_ptr_1(num_total_nodes);
    for (i = 0; i < num_total_nodes; i++) {
      SS_list[i] = alloc_int_1(MAX_SS_PER_NODE, -1);
    }

    /* Loop over the total number of side sets defined on the current
     * processor 
     */
    iptr = 0;
    for (iss = 0; iss < exo->num_side_sets; iss++) {
      /*
       * Step through every elem/side and look at every node on that
       * elem/side...
       */
      for (side_index = 0; side_index < exo->ss_num_sides[iss]; side_index++) {
        /*
         * Find the element, element_block_index, and the material
         * index corresponding to this side
         */
        ielem = exo->ss_elem_list[exo->ss_elem_index[iss] + side_index];
        elem_block_index = find_elemblock_index(ielem, exo);
        matIndex = Matilda[elem_block_index];
        /*
         * Loop over all of the local nodes on the side 
         */
        for (lni = exo->ss_node_side_index[iss][side_index];
             lni < exo->ss_node_side_index[iss][side_index+1]; lni++ ) {
          inode = exo->ss_node_list[iss][lni];

          if (in_list(exo->ss_id[iss], 0, MAX_SS_PER_NODE, SS_list[inode]) == -1) {
            iptr = 0;
            while (SS_list[inode][iptr] != -1) iptr++;
            if (iptr >= MAX_SS_PER_NODE) {
              sprintf(Err_Msg, 
                      "Node (%d) wants > %d sidesets - boost MAX_SS_PER_NODE.", 
                      inode+1, MAX_SS_PER_NODE);
              EH(GOMA_ERROR, Err_Msg);
            }
	
            SS_list[inode][iptr] = exo->ss_id[iss];
          }
        }
      }
    }

  
    for (ibc = 0; ibc < Num_BC; ibc++) {

      /* First check to see if this boundary condition is a side set      */
      if (!strcmp(BC_Types[ibc].Set_Type, "SS")) {

        /*
         * Determine what entity this boundary condition is applied to
         * - face, edge, or vertex.
         */
        faceBC = edgeBC = vertexBC = FALSE;
        if (BC_Types[ibc].BC_ID  != -1 && 
            BC_Types[ibc].BC_ID2 == -1 && 
            BC_Types[ibc].BC_ID3 == -1 ) {
          faceBC = TRUE;
        } else if (BC_Types[ibc].BC_ID  != -1 && 
                   BC_Types[ibc].BC_ID2 != -1 && 
                   BC_Types[ibc].BC_ID3 == -1 ) {
          edgeBC = TRUE;
        } else if (BC_Types[ibc].BC_ID  != -1 && 
                   BC_Types[ibc].BC_ID2 != -1 && 
                   BC_Types[ibc].BC_ID3 != -1 ) {
          vertexBC = TRUE;
        }

        /* Loop over the total number of side sets defined on the
         * current processor 
         */
        for (iss = 0; iss < exo->num_side_sets; iss++) {
          /* Check for a match between the ID of the current side set
           * and the side set ID specified in the input file 
           * - continue if a match is found 
           */
          if (exo->ss_id[iss] == BC_Types[ibc].BC_ID) {

            for (side_index = 0; side_index < exo->ss_num_sides[iss]; 
                 side_index++) {
              /*
               * Find the element, element_block_index, and the material
               * index corresponding to this side
               */
              ielem = exo->ss_elem_list[exo->ss_elem_index[iss] + side_index];
              elem_block_index = find_elemblock_index(ielem, exo);
              matIndex = Matilda[elem_block_index];
              /*
               * Loop over all of the local nodes on the side 
               */
              for (lni = exo->ss_node_side_index[iss][side_index];
                   lni < exo->ss_node_side_index[iss][side_index+1]; lni++ ) {
                inode = exo->ss_node_list[iss][lni];
	
                /* 
                 * Determine if this node belongs to those the entity to 
                 * which this bc applies
                 */
                node_ok = 0;

                /* Face conditions */
                if (faceBC) node_ok = 1;

                /* Edge conditions */
                if (edgeBC &&
                    (in_list(BC_Types[ibc].BC_ID2, 0, MAX_SS_PER_NODE, 
                             SS_list[inode]) != -1)) {
		  node_ok = 1;
                }

                /* Vertex conditons */
                if (vertexBC &&
                    (in_list(BC_Types[ibc].BC_ID2, 0, MAX_SS_PER_NODE, 
                             SS_list[inode]) != -1)   &&
                    (in_list(BC_Types[ibc].BC_ID3, 0, MAX_SS_PER_NODE, 
                             SS_list[inode]) != -1)) {
		  node_ok = 1;
                }

                if (node_ok) {
                  /*
                   *  Allocate temp BC storage for this node and initialize
                   *  the storage to -1
                   */
                  if (BC_Unk_List[inode] == NULL) {
                    num_bc_nodes++;
                    BC_Unk_List[inode] = alloc_bc_unk_list_node(inode);
                  }

                  /*
                   * Fill in an entry into the list for each component of the
                   * boundary condition
                   */
                  for (p = 0; p < BC_Types[ibc].desc->vector; p++) {
                    offset = find_bc_unk_offset(BC_Types + ibc, matIndex, inode, p,
                                                &retn_matIndex, &vd);
                    if (offset >= 0) {
                      add_to_node_unk_bc_list(BC_Unk_List[inode][offset], ibc);
                    }
                  }
                } /* if ( node_ok ) */
              } /* for ( lni=exo->ss_node_side_index[iss][side_index]; */
            } /*for ( side_index=0; side_index<exo->ss_num_sides[iss]; side_index++) */
          }  /* if (Proc_SS_Ids[iss] == BC_Types[ibc].BC_ID)  */
        }  /*  for (iss = 0; iss < Proc_Num_Side_Sets; iss++) */
      }  /* END if (!strcmp(BC_Types[ibc].Set_Type, "SS")) 		      */
    }  /* END for (ibc = 0; ibc < Num_BC; ibc++)				      */

    /***************************************************************************/
    /*
     * SECOND make a list of all the rotation conditions from the input deck
     *  at each node
     */
    /***************************************************************************/
    /* create an array which will hold boundary condition information temporarily
     * for each node 
     */
    ROT_list = (int **) alloc_ptr_1(num_total_nodes);
    num_rot_nodes = 0;
  
    /***************************************************************************/
    /*      ROTATION CONDITIONS ONLY ON SIDE SETS                              */
    /*  visit each node on all side sets and list the ROT's applied there      */
    /***************************************************************************/

    for (irc = 0; irc < Num_ROT; irc++) {

      /* Check for consistency. I.E. make sure all side sets on cards exist */

      /*
       * New consistency check - the side set ID for this ROT condition may
       * well exist on another processor and not on this processor. This should
       * be "OK".
       */

      if (Num_Proc == 1) {
        len_ssid_list = exo->num_side_sets;
        ssid_list     = exo->ss_id;
      } else {
        len_ssid_list = dpi->num_side_sets_global;
        ssid_list     = dpi->ss_id_global;
      }

      /*
       * If ss_id1 for this ROT condition ain't in the ss_id list then
       * it's invalid.
       */

      if ( -1 == in_list(ROT_Types[irc].ss_id[0], 0, len_ssid_list, ssid_list)) {
        sprintf(Err_Msg, 
                "ROT[%d] (ROT = %s %s %d ...) has ss_id_1 = %d not in mesh.",
                irc, rot_eq_type2string(ROT_Types[irc].eq_type),
                rotopology_type2string(ROT_Types[irc].type),
                ROT_Types[irc].ss_id[0], 
                ROT_Types[irc].ss_id[0]);
        fprintf(stderr, "Warning: %s\n", Err_Msg);
#if 0	
        /* Leave this barn door open to help test further... */
        EH(GOMA_ERROR, err_msg);
#endif
      }

      /* 
       * For each rotation card, Loop over the total number of nodes
       * checking to see if the rotation card affects that node.
       */
      ROT_Types[irc].node = -1;
      for (inode = 0; inode < num_total_nodes; inode++) {
        /*
         * Check for a match between the ID of the current side set and the 
         * side set ID specified in the input file - continue if a match 
         * is found 
         */
        node_ok = 0;
        if (ROT_Types[irc].type == FACE &&
            (in_list(ROT_Types[irc].ss_id[0], 0, MAX_SS_PER_NODE, 
                     SS_list[inode]) != -1)) {
          node_ok = 1;
        }
        if (ROT_Types[irc].type == CURVE &&
            (in_list(ROT_Types[irc].ss_id[0], 0, MAX_SS_PER_NODE, 
                     SS_list[inode]) != -1) &&
            (in_list(ROT_Types[irc].ss_id[1], 0, MAX_SS_PER_NODE,
                     SS_list[inode]) != -1)) {
          node_ok = 1;
        }
        if (ROT_Types[irc].type == VERTEX &&
            (in_list(ROT_Types[irc].ss_id[0], 0, MAX_SS_PER_NODE,
                     SS_list[inode]) != -1) &&
            (in_list(ROT_Types[irc].ss_id[1], 0, MAX_SS_PER_NODE,
                     SS_list[inode]) != -1) &&
            (in_list(ROT_Types[irc].ss_id[2], 0, MAX_SS_PER_NODE,
                     SS_list[inode]) != -1) ) {
          ROT_Types[irc].node = inode;
          node_ok = 1;
        }

        if (node_ok) {

          /*
           * allocate temp BC storage for this node and initialize
           */
          if (ROT_list[inode] == NULL) {
            num_rot_nodes++;
            ROT_list[inode] = alloc_int_1(NUM_VECTOR_EQUATIONS, -1);
          }
          if      (ROT_Types[irc].eq_type == R_MESH1) eq = VECT_EQ_MESH;
          else if (ROT_Types[irc].eq_type == R_MOMENTUM1) eq = VECT_EQ_MOM;
          else EH(GOMA_ERROR,"Illegal vector equation");

          if (ROT_list[inode][eq] == -1) {
            ROT_list[inode][eq] = irc;
          } else {
            /*
             * Resolve conflicts between conflicting rotation conditions.
             * Give precedence to vertex, then edge, then face then body 
             */
            if (ROT_Types[irc].type > ROT_Types[ROT_list[inode][eq]].type) {
              ROT_list[inode][eq] = irc;
            } else if (ROT_Types[irc].type == ROT_Types[ROT_list[inode][eq]].type) {
              /* 
               * If two rotations of same type - give precedence to rotation which
               * is first in input file
               */
              if (irc < ROT_list[inode][eq]) {
                ROT_list[inode][eq] = irc;
              }
            }
          }
        }  /* if (node_ok)  */
      }  /*  end of loop over nodes */
    }  /* END for (irc = 0; irc < Num_ROT; irc++)				     */

#ifndef PARALLEL
    fprintf (stderr, "Found %d nodes at which to rotate equations\n", 
             num_rot_nodes);
#endif
#ifdef PARALLEL
    DPRINTF(stdout, "P_%d: rotating equations at %d nodes.\n", ProcID,
            num_rot_nodes);
#endif

    /******************************************************************************/
    /*      RESOLVE BOUNDARY CONDITION CONFLICTS SPECIFIED                        */
    /******************************************************************************/
    /* I know this section looks like a real mess, and I'm sorry that it'll be hard
     * to follow.  Should anyone want to make changes here, or try to understand what
     * logic (?) is built into this routine - see the BC conflict resolution flow-chart
     * in the GOMA manual documentation and try to believe that this routine follows 
     * that flow chart somewhat religiously
     */
    
    BC_dup_nodes[pg->imtrx] = alloc_int_1(num_bc_nodes, -1);
    BC_dup_list[pg->imtrx] = (int ***) alloc_ptr_1(num_bc_nodes);
    BC_dup_ptr[pg->imtrx] = 0;

    /* 
     * loop through nodes and rearrange the BC list 
     * - remove nodes that don't have conflicts 
     * - count nodes with conflicts and build an array to store them
     */
    for (inode = 0; inode < num_total_nodes; inode++) {
      if (BC_Unk_List[inode] != NULL) {
        node = Nodes[inode];
        nv = node->Nodal_Vars_Info[pg->imtrx];
        /* Currently don't worry about rotation unless a BC exists at a 
         * node 
         */
        dups = 0;
        for (offset = 0; offset < nv->Num_Unknowns; offset++) {
          get_nv_vd_from_offset(nv, offset, &vd, &idof);
          eqn = vd->Variable_Type;

          /* count duplications */
          ndup = 0;
          while (BC_Unk_List[inode][offset][ndup] != -1) {
            ndup++;
          }
	
          if (ROT_list[inode] == NULL) {
            standard_BC = 1;
          } 
          else if ((eqn == R_MESH1 || eqn == R_MESH2 || eqn == R_MESH3)
                   && ROT_list[inode][VECT_EQ_MESH] == -1) {
            standard_BC = 1;
          } 
          else if ((eqn == R_MOMENTUM1 || eqn == R_MOMENTUM2 || 
                    eqn == R_MOMENTUM3)
                   && ROT_list[inode][VECT_EQ_MOM] == -1) {
            standard_BC = 1;
          }  else  {
            /* use ROT_Types to resolve BC's */
            standard_BC = 0;

            /* Resolve all three mesh or momentum BC's the first time 
             * through 
             */
            if (eqn == R_MESH1 || eqn == R_MOMENTUM1) {
              /* search BC_list to find BC's and SS or NS numbers that
               * match with ROT_list then change BC_list to include 
               * only those BC's and all WEAK BC's 
               */
              if (eqn == R_MESH1)     eq = VECT_EQ_MESH;
              if (eqn == R_MOMENTUM1) eq = VECT_EQ_MOM;
              irc = ROT_list[inode][eq];

              /* initialize new entry in the duplications list, 
               * if not already listed 
               */
              if (in_list(inode, 0, BC_dup_ptr[pg->imtrx], BC_dup_nodes[pg->imtrx]) == -1) {
                BC_dup_nodes[pg->imtrx][BC_dup_ptr[pg->imtrx]] = inode; /* put node in 
                                                   * list 
                                                   */
                BC_dup_list[pg->imtrx][BC_dup_ptr[pg->imtrx]]  = BC_Unk_List[inode]; /* ptr to
                                                                * eqn 
                                                                * BC_list
                                                                */
                BC_dup_ptr[pg->imtrx]++;
                if (Debug_Flag > 1 && Unlimited_Output) {
                  fprintf(ofbc, 
                          "New Rot node in dup list %d %d\n",
                          inode+1,BC_dup_ptr[pg->imtrx]);
                }
              }

              /* Now place the BC's in their correct spots, 
               * allowing all weak conditions 
               */
              /* count duplications */
              for (p = 0; p < dim; p++) {
                ndup1[p] = 0;
                offset_mom2 = 
		  get_nodal_unknown_offset(nv, eqn + p, vd->MatID, 0, 
					   &vd_retn2);
                while (BC_Unk_List[inode][offset_mom2][ndup1[p]] != -1) {
                  ndup1[p]++;
                }
                for (j = 0; j < ndup1[p]; j++) {
                  save_this_bc1[p][j] = 0;
                }
              }

              /* loop over dimensions in Rotation specification */
              for (p=0; p<dim; p++) {
                bct  = ROT_Types[irc].BC_Type[p];
                bcSS = ROT_Types[irc].BC_SS[p];
#ifndef PARALLEL
                bc_found = 0;
#endif

                if (bct < 0) {
                  /* don't do anything this is a rotated piece */
#ifndef PARALLEL
                  bc_found = 1;
#endif
                } else {
                  /* find BC in dup list */

                  /*
                   * But it might not be here on this processor!
                   */

                  for (q = 0; q < dim; q++) {
                    for (j = 0; j < ndup1[q]; j++) {
                      offset_mom2 = 
			get_nodal_unknown_offset(nv, eqn + q, vd->MatID, 0, 
						 &vd_retn2);
                      ibc1 = BC_Unk_List[inode][offset_mom2][j];
                      bct1 = BC_Types[ibc1].BC_Name;
                      bcSS1 = BC_Types[ibc1].BC_ID;

                      /* check to see that the BC's are 
                       * the same and that this BC hasn't 
                       * already been checked 
                       */

                      /*
                       * More detailed diagnostic why
                       * not found...
                       */

#ifndef PARALLEL
                      tilt_bct  = ( bct != bct1 );
                      tilt_bcSS = ( bcSS != bcSS1 );
                      tilt_save = ( save_this_bc1[q][j] != 0 );
                      val       = save_this_bc1[q][j];
#endif

                      if ( bct == bct1 && 
                           bcSS == bcSS1 && 
                           save_this_bc1[q][j] == 0) {
                        /* 
                         * Found the BC which is specified 
                         * by ROT condition!! 
                         */
#ifndef PARALLEL
                        bc_found = 1;
#endif
                        ROT_Types[irc].BC_desc[p] = 
			  BC_Types[ibc1].desc;
                        ROT_Types[irc].BC_id[p] = ibc1;

                        /* find a spot for this BC in the 
                         * dup_list */
                        save_this_bc1[p][ndup1[p]] = 1;
                        offset_mom2 = 
			  get_nodal_unknown_offset(nv, eqn + p, vd->MatID, 0, 
                                                   &vd_retn2);
                        BC_Unk_List[inode][offset_mom2][ndup1[p]] = ibc1;
                        ndup1[p]++;
                        /* need to deal with GD conditions 
                         * here!! */

                        if ( BC_Types[ibc1].desc->method == 
                             DIRICHLET && p != q) { 
                          EH(GOMA_ERROR, "Can't change equation type of DIRICHLET!");
                        }
                      }

                      /* retain all weak conditions */
                      if ( BC_Types[ibc1].desc->method == WEAK_INT_SURF ||
                           BC_Types[ibc1].desc->method == WEAK_SHARP_INT  ){ 
                        save_this_bc1[q][j] = 1;
		      
                        /* retain all weak conditions */
                      } 
                      else if ( BC_Types[ibc1].desc->method ==
                                WEAK_SHIFT) { 
                        save_this_bc1[q][j] = 1;

                        /* retain all weak conditions */
                      } 
                      else if ( BC_Types[ibc1].desc->method ==
                                WEAK_INT_EDGE) { 
                        save_this_bc1[q][j] = 1;
                      }
		  
                      /* retain all Dirichlet conditions */
                      /* 		    } else if (BC_Types[ibc1].desc->method == DIRICHLET) {  */
                      /* 		      save_this_bc1[q][j] = 1; */
                      /* 		    } */
                    }  
                  }
                }

                /*
                 * It's got to be "OK" for the BC not to exist on 
                 * this processor.
                 */
#ifndef PARALLEL		/* This is dangerous!!!!!!! */
                if ( bc_found != 1 ) {
                  fprintf(ofbc, 
                          "Cannot find BC to match ROT [%d], coord %d at node %d\n",
                          irc, p, inode);
                  fprintf(stderr, 
                          "Problem on ROT [%d], coord [%d], node (%d)\n",
                          irc, p, inode+1);
                  if ( tilt_bct && tilt_bcSS && tilt_save ) {
                    fprintf(stderr, 
                            "P_%d: bct = %d, bct1 = %d\n",
                            ProcID, bct, bct1);
                    fprintf(stderr, 
                            "P_%d: bcSS = %d, bcSS1 = %d\n",
                            ProcID, bcSS, bcSS1);
                    fprintf(stderr, 
                            "P_%d: save_this_bc1[q][j]=%d\n",
                            ProcID, val);
                  }
                  fprintf(stderr, "P_%d death, nss = %d\n",
                          ProcID, exo->num_side_sets);
                  EH(GOMA_ERROR,"Rotation BC not found");
                }
#endif

              } /* end of loop over rotation directions */
              for (p = 0; p < dim; p++) {
                offset_mom2 = 
		  get_nodal_unknown_offset(nv, eqn + p, vd->MatID, 0, 
					   &vd_retn2);
                /* redo list keeping only those in save_this_bc */
                k = 0;
                for (j = 0; j < ndup1[p]; j++) {
                  if (save_this_bc1[p][j]) {
                    BC_Unk_List[inode][offset_mom2][k] = 
		      BC_Unk_List[inode][offset_mom2][j];
                    k++;
                  }
                }
                for (j = k; j < ndup1[p]; j++) {
                  BC_Unk_List[inode][offset_mom2][j] = -1;
                }
              }
            }
          }
	
          if (standard_BC) {
            /* treat as in 2D, but don't shift any mesh or momentum BC's */
            /* ORDER of PRECEDENCE
             *  - allow all weak conditions in addition to:
             *  1. 1st DIRICHLET or
             *  2. 1st COLLOCATED non CSWI or 
             *  3. 1st INTEGRATED
             *  4. 1st COLLOCATED CSWI
             *
             * (here CSWI stands for Collocated Surface bc which are
             *  weak integrals in reality)
             */
            /*
             * RESOLVE DUPLICATIONS AT THIS NODE
             */
            if (ndup > 1) {

              /* initialize new entry in the duplications list,
               *  if not already listed 
               */
              if (in_list(inode, 0, BC_dup_ptr[pg->imtrx], BC_dup_nodes[pg->imtrx]) == -1) {
                BC_dup_nodes[pg->imtrx][BC_dup_ptr[pg->imtrx]] = inode;
                BC_dup_list[pg->imtrx][BC_dup_ptr[pg->imtrx]] = BC_Unk_List[inode];
                BC_dup_ptr[pg->imtrx]++;
                if (Debug_Flag > 1 && Unlimited_Output ) {
                  fprintf(ofbc, 
                          "New node in duplication list %d %d\n",inode+1,BC_dup_ptr[pg->imtrx]);
                }
              }

              idup=0;
              if (Debug_Flag > 0 && Unlimited_Output) { 
                fprintf(ofbc,
                        "  DUPLICATE BC at node %d equation %d: ", inode+1, eqn);
              }
              while (BC_Unk_List[inode][offset][idup] != -1) {
                ibc = BC_Unk_List[inode][offset][idup];
                if (Debug_Flag > 0 && Unlimited_Output) { 
                  fprintf(ofbc, "%s, ", BC_Types[ibc].desc->name2);
                }
                idup++;
              }
              if (Debug_Flag > 0 && Unlimited_Output) { 
                fprintf(ofbc, "\n");
              }
              dups++;		/* increment count of # duplications at this node */
	    
              /*********************************************************************/
              /* try to resolve the differences quickly */
              /* when two BC's are applied at a point, there are three possibilities
               * for how to choose to apply the BC's:
               *  1- One of the BC's takes precedence and the other is discarded
               *     regardless of their coefficients
               *     - the first one in the input deck takes precedence
               *  1. 1st DIRICHLET or
               *  2. 1st COLLOCATED or 
               *  3. 1st INTEGRATED
               *  2- Both BC's get applied - regardless of their coefficients
               */
              /*********************************************************************/

              /* make array which indicates if a given BC is to be retained */
              for (j = 0; j < ndup; j++) save_this_bc[j] = 0;
              j_DC = -1;
              j_PC = -1;
              j_SI = -1;
              j_CSWI = -1;
              for (j = 0; j < ndup; j++) {
                ibc1 = BC_Unk_List[inode][offset][j];
                bct1 = BC_Types[ibc1].BC_Name;

                /* Check to make sure user isn't trying to rotate without 
                 * specifying a rotation condition */
                if (BC_Types[ibc1].desc->rotate != NO_ROT) {
                  sprintf(Err_Msg, 
                          "BC [%d] %s %2s %d is a rotated condition, but lacks explicit ROT",
                          ibc1, BC_Types[ibc1].desc->name1,
                          BC_Types[ibc1].Set_Type, BC_Types[ibc1].BC_ID);
                  EH(GOMA_ERROR, Err_Msg);
                }

                /* save all WEAK conditions */
                if (BC_Types[ibc1].desc->method == WEAK_INT_SURF ||
                    BC_Types[ibc1].desc->method == WEAK_SHARP_INT ) { 
                  save_this_bc[j] = 1;

                } else if (BC_Types[ibc1].desc->method == WEAK_SHIFT) { 
                  save_this_bc[j] = 1;

                } else if (BC_Types[ibc1].desc->method == WEAK_INT_EDGE) { 
                  save_this_bc[j] = 1;

                  /* DIRICHLET has precedence */
                } else if (BC_Types[ibc1].desc->method == DIRICHLET) {
                  if (j_DC != -1) {
                    ibc2 = BC_Unk_List[inode][offset][j_DC];
                    bct2 = BC_Types[ibc2].BC_Name;
                    /* RESOLVE CONFLICT */
                    if (Debug_Flag > 0 && Unlimited_Output) {
                      ibc = BC_Unk_List[inode][offset][j_DC];
                      fprintf(ofbc, "    TWO DBC's %s and %s at node %d \n", 
                              BC_Types[ibc1].desc->name2, 
                              BC_Types[ibc].desc->name2,
                              inode+1);
                    }
                  } else {
                    save_this_bc[j] = 1;
                    j_DC = j;
                  }
                }
              } /* end of loop over dups */

              /* only check for collocated condition if no Dirichlet exist */
              if (j_DC == -1) {
                for (j = 0; j < ndup; j++) {
                  ibc1 = BC_Unk_List[inode][offset][j];
                  bct1 = BC_Types[ibc1].BC_Name;
                  if ((BC_Types[ibc1].desc->method == COLLOCATE_SURF) &&
                      (bct1 != FLUID_SOLID_BC &&
                       bct1 != SOLID_FLUID_BC &&
                       bct1 != SOLID_FLUID_RS_BC                    )) {
                    if (j_PC != -1) {
                      ibc2 = BC_Unk_List[inode][offset][j_PC];
                      bct2 = BC_Types[ibc2].BC_Name;
                      /* RESOLVE CONFLICT */
                      /* need to deal with GD conditions here 
                       * - two GD conditions on the same SS and
                       * applied to the same equation do not conflict with each other 
                       */
                      if ((bct1 <= GD_TIME_BC && bct1 >= GD_CONST_BC) && 
                          (bct2 <= GD_TIME_BC && bct2 >= GD_CONST_BC)) {
                        if (BC_Types[ibc1].BC_ID == BC_Types[ibc2].BC_ID) {
                          if (Debug_Flag > 0 && Unlimited_Output) {
                            fprintf(ofbc, 
                                    "    MULTIPLE GD %s at node %d\n", 
                                    BC_Types[ibc2].desc->name2, inode+1);
                          }
                          save_this_bc[j] = 1;
                        } else {
                          if (Debug_Flag > 0 && Unlimited_Output) {
                            fprintf(ofbc, 
                                    "    REMOVE %s in favor of %s at node %d\n", 
                                    BC_Types[ibc2].desc->name2, 
                                    BC_Types[ibc1].desc->name2, inode+1);
                          }
                        }
                      } else {
                        if (Debug_Flag > 0 && Unlimited_Output) {
                          fprintf(ofbc, 
                                  "    REMOVE %s in favor of %s at node %d\n",
                                  BC_Types[ibc2].desc->name2, 
                                  BC_Types[ibc1].desc->name2, inode+1);
                        }
                      }
		  
                    } else {
                      save_this_bc[j] = 1;
                      j_PC = j;
                    }
                  }
                }
	    
                /* only check for integrated if no collocated or Dirichlet exist */
                if (j_PC == -1) {
                  for (j = 0; j < ndup; j++) {
                    ibc1 = BC_Unk_List[inode][offset][j];
                    bct1 = BC_Types[ibc1].BC_Name;
                    if (BC_Types[ibc1].desc->method == STRONG_INT_SURF) {
                      if (j_SI != -1) {
                        ibc2 = BC_Unk_List[inode][offset][j_PC];
                        bct2 = BC_Types[ibc2].BC_Name;
                        /* RESOLVE CONFLICT */
                        if (Debug_Flag > 0 && Unlimited_Output) { 
                          fprintf(ofbc, 
                                  "    REMOVE %s in favor of %s at node %d\n", 
                                  BC_Types[ibc2].desc->name2, 
                                  BC_Types[ibc1].desc->name2, inode+1);
                        }
                      } else {
                        save_this_bc[j] = 1;
                        j_SI = j;
                      }
                    }
                  }
                  /* Only check if no Dirichlet, (other) collocated,
                   * or strongly integrated bc's exist
                   */
                  if (j_SI == -1) {
                    for (j = 0; j < ndup; j++) {
                      ibc1 = BC_Unk_List[inode][offset][j];
                      bct1 = BC_Types[ibc1].BC_Name;
                      if ((BC_Types[ibc1].desc->method == COLLOCATE_SURF) &&
                          (bct1 == FLUID_SOLID_BC ||
                           bct1 == SOLID_FLUID_BC ||
                           bct1 == SOLID_FLUID_RS_BC)) {
                        if (j_CSWI != -1) {
                          ibc2 = BC_Unk_List[inode][offset][j_CSWI];
                          bct2 = BC_Types[ibc2].BC_Name;
                          if (Debug_Flag > 0 && Unlimited_Output) { 
                            fprintf(ofbc, 
                                    "    REMOVE %s in favor of %s at node %d\n", 
                                    BC_Types[ibc1].desc->name2, 
                                    BC_Types[ibc2].desc->name2, inode+1);
                          }
                        } else {
                          save_this_bc[j] = 1;
                          j_CSWI = j;
                        }
                      }
                    }	  

                  } /* end of if SI doesn't exist */
                }	/* end of if PC doesnt exist */
              } /* end of if DC doesnt exist */
	      
              /* redo list keeping only those in save_this_bc */
              k = 0;
              for (j = 0; j < ndup; j++) {
                if (save_this_bc[j]) {
                  BC_Unk_List[inode][offset][k] = BC_Unk_List[inode][offset][j];
                  k++;
                }
              }
              for (j = k; j < ndup; j++) {
                BC_Unk_List[inode][offset][j] = -1;
              }
            } /* end of if ndup > 1 */
          } /* end of if standard BC */
        }	/* end of loop over equations */
      } /* end of if BC_list */
    } /* end of loop over nodes */

    /******************************************************************************/
    /* CLEAN UP                                                                   */
    /******************************************************************************/

    /* Now print out BC_dup list for 3-D problems */
    /* Also check for potentially conflicting combinations */
    if (pd_glob[0]->Num_Dim == 3 ) {
      if (Unlimited_Output) {
        fprintf(ofbc, 
                "\nBoundary Condition Duplication List (numbers in paren are BC_ID):\n");
      }
      for (i = 0; i < BC_dup_ptr[pg->imtrx]; i++) {
        inode = BC_dup_nodes[pg->imtrx][i];
        node = Nodes[inode];
        if (Unlimited_Output) {
	  fprintf(ofbc, "N=%d:", BC_dup_nodes[pg->imtrx][i]+1);
        }
        nv = node->Nodal_Vars_Info[pg->imtrx];
        /* check for duplications */
        dups = 0;
        for (offset = 0; offset < nv->Num_Unknowns; offset++) {
          get_nv_vd_from_offset(nv, offset, &vd, &idof);
          eqn = vd->Variable_Type;
          eprint = 1;
          if (eqn >= R_MOMENTUM1 && eqn <= R_MOMENTUM3) {
            if (ROT_list[inode] != NULL && ROT_list[inode][VECT_EQ_MOM] >= 0) {
              if (ROT_Types[ROT_list[inode][VECT_EQ_MOM]].BC_Type[eqn - R_MOMENTUM1] < 0
                  && ROT_Types[ROT_list[inode][VECT_EQ_MOM]].ROTATE) {
                if( Unlimited_Output )
		  fprintf(ofbc, " E%d=ROTATE", eqn); 
                eprint = 0;
              }
            }
          }
          if (eqn >= R_MESH1 && eqn <= R_MESH3) {
            if (ROT_list[inode] != NULL && ROT_list[inode][VECT_EQ_MESH] >= 0) {
              if (ROT_Types[ROT_list[inode][VECT_EQ_MESH]].BC_Type[eqn - R_MESH1] < 0
                  && ROT_Types[ROT_list[inode][VECT_EQ_MESH]].ROTATE) {
                if( Unlimited_Output )
		  fprintf(ofbc, " E%d=ROTATE", eqn); 
                eprint = 0;
              }
            }
          }
          if (BC_Unk_List[inode][offset][0] != -1) {
            if (eprint && Unlimited_Output ) 
	      fprintf(ofbc, " E%d=", eqn);
            idup=0;
            while (BC_Unk_List[inode][offset][idup] != -1) {
              if (Unlimited_Output) {
                fprintf(ofbc, "%s(%d", 
                        BC_Types[BC_Unk_List[inode][offset][idup]].desc->name1, 
                        BC_Types[BC_Unk_List[inode][offset][idup]].BC_ID);
              }
              if (BC_Types[BC_Unk_List[inode][offset][idup]].BC_ID2 == -1) {
                if (Unlimited_Output) {
                  fprintf(ofbc, ")");
                }
              } else {
                if (Unlimited_Output) {
                  fprintf(ofbc, ",%d)", 
                          BC_Types[BC_Unk_List[inode][offset][idup]].BC_ID2);
                }
              }
              idup++;
            }
          }
        }
        if (Unlimited_Output) fprintf(ofbc, "\n");
      }
    }
    if (Unlimited_Output) fprintf(ofbc, "\n");

    /* free memory associated with nodes that don't have duplicate BC's */
    for (inode = 0; inode < num_total_nodes; inode++) {
      if (BC_Unk_List[inode] != NULL) {
        if (in_list(inode, 0, BC_dup_ptr[pg->imtrx], BC_dup_nodes[pg->imtrx]) == -1) {
          free_bc_unk_list_node(BC_Unk_List + inode, inode);
        }
      }
    }

    /* remove array which held boundary condition information temporarily
     * for each node 
     */
    safer_free((void **) &BC_Unk_List);
    for ( i=0; i<num_total_nodes; i++) {
      safer_free((void **) (NS_list + i));
    }
    safer_free((void **) &NS_list);

  }

  if (Unlimited_Output) {
    fclose(ofbc);
    print_sync_end(TRUE);
  }

  pg->imtrx = 0;

  return;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
