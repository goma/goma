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
 *$Id: rf_setup_problem.c,v 5.4 2009-04-23 22:49:05 hkmoffa Exp $
 */

#include <stdio.h>
#include <string.h>

#include "std.h"
#include "rf_allo.h"
#include "rf_fem_const.h"
#include "rf_fem.h"
#include "rf_io.h"
#include "rf_vars_const.h"
#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "mm_as.h"
#include "rf_mp.h"
#include "exo_struct.h"
#include "dp_types.h"
#include "dpi.h"
#include "dp_map_comm_vec.h"
#include "mm_unknown_map.h"
#include "mm_bc.h"
#include "rd_mesh.h"
#include "rf_node_const.h"
#include "rf_solve.h"
#include "mm_eh.h"
#include "rf_util.h"
#include "rf_bc_const.h"
#include "rf_bc.h"
#include "rf_element_storage_const.h"
#include "mm_shell_util.h"
#include "dp_utils.h"
#include "mpi.h"

static void
associate_bc_to_matrix(void);

static void
set_bc_equation(void)
{
  int ibc;

  for (ibc = 0; ibc < Num_BC; ibc++) {
    /*
     *  Create a couple of pointers to cut down on the
     *  amount of indirect addressing
     */
    int eqn = (BC_Types[ibc].desc)->equation;
    int set_eqn = BC_Types[ibc].equation;
    if (set_eqn < 0 || set_eqn >= V_LAST) {
      if (eqn >= V_FIRST && eqn < V_LAST) {
        BC_Types[ibc].equation = eqn;
      } else {
        switch (eqn) {
        case R_MESH_NORMAL:
        case R_MESH_TANG1:
        case R_MESH_TANG2:
        case MESH_POSITION1:
          BC_Types[ibc].equation = R_MESH1;
          break;
        case MESH_POSITION2:
          BC_Types[ibc].equation = R_MESH2;
          break;
        case MESH_POSITION3:
          BC_Types[ibc].equation = R_MESH3;
          break;
        case R_MOM_NORMAL:
        case R_MOM_TANG1:
        case R_MOM_TANG2:
           BC_Types[ibc].equation = R_MOMENTUM1;
          break;
        case R_SOLID_NORMAL:
        case R_SOLID_TANG1:
        case R_SOLID_TANG2:
        case SOLID_NORM:
        case SOLID_TANG1:
        case SOLID_TANG2:
        case SOLID_POSITION1:
          BC_Types[ibc].equation = R_SOLID1;
          break;
        case SOLID_POSITION2:
          BC_Types[ibc].equation = R_SOLID2;
          break;
        case SOLID_POSITION3:
          BC_Types[ibc].equation = R_SOLID3;
          break;
        default:
          EH(-1, "Error in linking BC to equation");
          break;
        }
      }
    }
  }
}

static void
associate_bc_to_matrix(void)
{
  int ibc;
  int imtrx;
  int mn;

  if (upd->Total_Num_Matrices == 1) {
    /* Preserve legacy behavior */
    for (ibc = 0; ibc < Num_BC; ibc++) {
    
      BC_Types[ibc].matrix = 0;
      
    }
  } else {
    for (ibc = 0; ibc < Num_BC; ibc++) {
      /*
       *  Create a couple of pointers to cut down on the
       *  amount of indirect addressing
       */
      int eqn = BC_Types[ibc].equation;

      BC_Types[ibc].matrix = -1;

      if (eqn >= V_FIRST && eqn < V_LAST) {
	for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
	  for (mn = 0; mn < upd->Num_Mat; mn++) {
	    if (pd_glob[mn]->e[imtrx][eqn]) {
	      BC_Types[ibc].matrix = imtrx;
	    }
	  }
	}
      }
      if (BC_Types[ibc].matrix == -1) {
	char errstr[512];
	snprintf(errstr, 512, "Could not find matching matrix for BC #%d %s, BC will not be used", 
		 ibc, (BC_Types[ibc].desc)->name1);
	WH(-1, errstr);
      }
    }
  }
}

void
set_matrix_index_and_global_v(void)
{
  int mn;
  int i;
  int imtrx;

  /* Initialize matrix indices */
  for (mn = 0; mn < upd->Num_Mat; mn++) {
    for ( i=0; i<MAX_VARIABLE_TYPES; i++) {
      pd_glob[mn]->mi[i] = -1;
      upd->matrix_index[i] = -1;
    }
  }

  for (mn = 0; mn < upd->Num_Mat; mn++) {
    for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
      for ( i=0; i<MAX_VARIABLE_TYPES; i++) {
        if (pd_glob[mn]->v[imtrx][i]) {
          pd_glob[mn]->mi[i] = imtrx;
          upd->matrix_index[i] = imtrx;
        }
      }
    }
  }

  for (mn = 0; mn < upd->Num_Mat; mn++) {
    for ( i=0; i<MAX_VARIABLE_TYPES; i++) {
      pd_glob[mn]->gv[i] = 0;
    }
  }

  for (mn = 0; mn < upd->Num_Mat; mn++) {
    for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
      for ( i=0; i<MAX_VARIABLE_TYPES; i++) {
        if (pd_glob[mn]->v[imtrx][i]) {
          pd_glob[mn]->gv[i] = 1;
        }
      }
    }
  }
}


/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

int setup_problem(Exo_DB *exo,	/* ptr to the finite element mesh database */
	          Dpi *dpi)	/* distributed processing information */

    /********************************************************************
     *
     * setup_problem():
     *
     *      Setup_problem() determines the degrees of freedom at each
     * node and formulates the solution vector. It then determines the
     * communications pattern for exchanging that solution vector between
     * processors.
     *     Lastly, it sets up structures that help to carry out the
     * boundary condition integrals on side sets.
     *
     * NOTE:
     *   This function was formed by taking common parts out of 
     **********************************************************************/
{
  int i;
  static char *yo = "setup_problem:";

  /*
   *  Initialize nodal based structures pertaining to properties
   *  and boundary conditions
   */
  init_nodes(exo, dpi);

  /*
   *  Determine the side set or node set index that matches all of
   *  those boundary conditions. Store it in the boundary condition
   *  structure.
   */
  bc_set_index(exo);

  /*
   *  Determine what boundary conditions are on side sets that
   *  are internal-boundary side sets
   */
  bc_internal_boundary(exo);

  /*
   *  Determine the neighboring material indecises on all boundary
   *  conditions whether or not they be specified in the input deck
   */
  bc_matrl_index(exo);
  
  /*
   * Determine from the boundary conditions and problem description
   * structures what degrees of freedom are shared between materials
   * and which degrees of freedom generate multiple degrees of freedom
   * at each node
   */
  coordinate_discontinuous_variables(exo, dpi);
 
  /*
   * Reconcile the boundary conditions with material properties
   * entered into the databases.
   */
  reconcile_bc_to_matrl();

  /*
   * Determine the values of the DV_Indexing_Type field for boundary
   * conditions
   */
  determine_dvi_index();


  set_bc_equation();

  /* Link BCs to matrix by equation */
  associate_bc_to_matrix();


  /*
   * Enumerate my own degrees of freedom on this processor.
   */
  setup_local_nodal_vars(exo, dpi);

  /*
   * setup communications patterns between ghost and owned
   * nodes.
   */
  setup_nodal_comm_map(exo, dpi, cx);

  /*
   * Find the global maximum number of unknowns located at any one
   * node on any processor
   */
  MaxVarPerNode = find_MaxUnknownNode();
  
  /*
   * Exchange my idea of what materials I have at each node with
   * my surrounding processors. Make sure we are all in sync
   */
  setup_external_nodal_matrls(exo, dpi, cx[0]);

  /*
   * Exchange my idea of what degrees of freedom I have with my
   * surrounding processors. Make sure we are all in sync.
   * Owned nodes tell ghost nodes what variables are active
   * at that node.
   */
  setup_external_nodal_vars(exo, dpi, cx);

  /*
   * Finish setting the unknown map on this processor
   */
  set_unknown_map(exo, dpi);
  
  /*
   * Now determine the communications pattern that is necessary to
   * exchange the solution vector in as efficient a manner as
   * possible
   */
  log_msg("setup_dof_comm_map...");
  setup_dof_comm_map(exo, dpi, cx);

  /*
   * Output some statistics concerning the communications pattern
   */
  if (Num_Proc > 1) output_comm_stats(dpi, cx);

#ifdef COUPLED_FILL
  /*
   * I extracted this from setup_fill_comm_map because some of the
   * renormalization routines make use of them.
   */
  num_fill_unknowns      = count_vardofs(FILL, dpi->num_universe_nodes);
  internal_fill_unknowns = count_vardofs(FILL, dpi->num_internal_nodes);
  owned_fill_unknowns    = count_vardofs(FILL, 
		           (dpi->num_internal_nodes + dpi->num_boundary_nodes));
  boundary_fill_unknowns = owned_fill_unknowns - internal_fill_unknowns;
  external_fill_unknowns = num_fill_unknowns - owned_fill_unknowns;
#else /* COUPLED_FILL */
  /*
   * Set up the communications pattern for doing a solution
   * of the FILL variable type equations alone, if needed
   */
  if (Explicit_Fill) {
#ifdef DEBUG
    fprintf(stderr, "P_%d: setup_fill_comm_map() begins...\n", ProcID);
#endif /* DEBUG */
    log_msg("setup_fill_comm_map...");
    setup_fill_comm_map(exo, dpi, cx);
  }
#endif /* COUPLED_FILL */

  /*
   *  Possibly increase the number of variable descriptions to include
   *  those that are not part of the solution vector
   */
  vdesc_augment();

  /*
   *  Setup Boundary condition inter-connectivity
   *  -> Stefan flow bc's need to know what bc's contain the reaction
   *     info.
   *  -> YFLUX bc's need to be connected -> future implementation
   */
  set_up_BC_connectivity();

  /*
   *  Set up the structures necessary to carry out
   *  surface integrals
   */
#ifdef DEBUG
  fprintf(stderr, "P_%d set_up_Surf_BC() begins...\n",ProcID);
#endif
  log_msg("set_up_Surf_BC...");
  set_up_Surf_BC(First_Elem_Side_BC_Array, exo, dpi); 

  /*
   *  Set up the Edge boundary condition structures
   */
#ifdef DEBUG
  fprintf(stderr, "P_%d: set_up_Edge_BC() begins...\n",ProcID);
#endif
  log_msg("set_up_Edge_BC...");
  set_up_Edge_BC(First_Elem_Edge_BC_Array, exo, dpi); 

  /*
   * Set up "boundary" conditions on level set surfaces
   */
  set_up_Embedded_BC();
  
  /* Special 1D
   * Set up surface integral boundary conditions
   * that apply at single nodes.
   */

  setup_Point_BC(First_Elem_Side_BC_Array, exo, dpi);

  /*
   *  Print out the edge boudary condition structures
   *  if necessary
   */
  if (Debug_Flag) {
    for (pg->imtrx = 0; pg->imtrx < upd->Total_Num_Matrices; pg->imtrx++) {
      print_setup_Surf_BC(First_Elem_Side_BC_Array[pg->imtrx]);
    }
  }

  /*
   *  Malloc structures of size Num_Var_Info_Records
   */
  for (pg->imtrx = 0; pg->imtrx < upd->Total_Num_Matrices; pg->imtrx++) {
    ei[pg->imtrx]->VDindex_to_Lvdesc = alloc_int_1(Num_Var_Info_Records, -1);
  }
  for (i = 0; i < MAX_ELEMENT_INDICES_RELATED; i++) {
    eiRelated[i]->VDindex_to_Lvdesc = alloc_int_1(Num_Var_Info_Records, -1);
  }

  /*
   * Setup the storage for temporary quantities of interest that
   * are storred on a "per processor element" basis. These include
   * temporary storage of volumetric quadrature information
   */
  setup_element_storage();

  /*
   * Setup some structures for solving problems with shell elements.
   */
  init_shell_element_blocks(exo);

  /* Communicate non-shared but needed BC information */
  exchange_bc_info();

  return 0;
}
/************************************************************************/

int
free_problem(Exo_DB *exo,     /* ptr to the finite element mesh database */
	     Dpi *dpi)	      /* distributed processing information */
     /*
      *  this is the deconstructor routine to the setup_problem routine.
      *  It's designed to deallocate all the arrays and structures allocated in setup_problem
      */
{
  /*
   * Free up the First_Elem_Side_BC_Array array 
   */
    free_Surf_BC(First_Elem_Side_BC_Array, exo);
  free_Edge_BC(First_Elem_Edge_BC_Array, exo, dpi);
  return 0;
}
/************************************************************************/

static void
check_discontinuous_interp_type(PROBLEM_DESCRIPTION_STRUCT *curr_pd,
			        int var_type,
                                int imtrx)
     
    /********************************************************************
     *
     * check_discontinuous_interp_type
     *
     * -> Test whether an interpolation specified in the input deck
     *    is consistent with a discontinuous interpolation at
     *    the interface
     ********************************************************************/
{
  int *v_ptr = curr_pd->v[imtrx];
  int interp_type = curr_pd->i[imtrx][var_type];
  if (v_ptr[var_type] & V_MATSPECIFIC) {
    switch (interp_type) {
    case I_NOTHING:
    case I_Q1:
    case I_Q2:
    case I_Q2_LSA:
    case I_H3:
    case I_S2:
    case I_B3:
    case I_Q3:
    case I_Q4:
    case I_SP:
	fprintf(stderr,"check_interpolation_discontinuous ERROR");
	fprintf(stderr," var type to be discontinuous in mat %s\n,",
		curr_pd->MaterialName);
	fprintf(stderr,"\tbut incompatible interp type specified: %d\n",
		interp_type);
	EH(-1,"incompatible interp type");
	break;
	    case I_G0:
    case I_G1:
    case I_P0:
    case I_P1:
    case I_Q1_D:
    case I_Q2_D:
    case I_Q2_D_LSA:
    case I_PQ1:
    case I_PQ2:
         break;
    default:
     	fprintf(stderr,"check_interpolation_discontinuous ERROR");
	fprintf(stderr,"unknown interpolation type\n");
	EH(-1,"unknown interpolation type");
    }
  }
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

static void
turn_on_discontinuous(PROBLEM_DESCRIPTION_STRUCT *curr_pd,
		      const int var_type,
                      const int imtrx)
    
    /********************************************************************
     *
     * turn_on_discontinuous():
     *
     * -> Test whether a variable is active. If it is, then turn on
     *    the bit that denotes the interpolation is discontinuous at
     *    material boundaries.
     ********************************************************************/
{
  int *v_ptr = &(curr_pd->v[imtrx][var_type]); 
  if ((*v_ptr) & (V_SOLNVECTOR)) {
    *v_ptr |= V_MATSPECIFIC;
    check_discontinuous_interp_type(curr_pd, var_type, imtrx);
  }
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

int
coordinate_discontinuous_variables(Exo_DB *exo,	Dpi *dpi)	
    
    /********************************************************************
     *
     * coordinate_discontinuous_variables():
     *
     * -> Make sure we have the correct designations for the v field
     *    in the problem description structure for each material.
     *
     * -> Make sure that we have the same v field for all material
     *    types on all processors.
     *     
     *
     *******************************************************************/
{
  int ibc, eqn_type, ss_index, side_index, k, node_num, imat;
  int imtrx;
  int num_mat, mat_index, var_type, *ivec;
  UMI_LIST_STRUCT *curr_mat_list;
  NODE_INFO_STRUCT *node_ptr;
  PROBLEM_DESCRIPTION_STRUCT *curr_pd;

  /*
   *  Loop over the boundary conditions. If we have a cross
   *  phase discontinuous boundary condition, then we need to set the
   *  v field for the appropriate variable types on both sides of the
   *  interface to denote a discontinuous interpolation at the
   *  interface.
   */
  for (ibc = 0; ibc < Num_BC; ibc++) {
    if (BC_Types[ibc].desc->i_apply == CROSS_PHASE_DISCONTINUOUS) {
      eqn_type = BC_Types[ibc].desc->equation;
      /*
       * If we are applying a bc on the momentum equations
       * let's assign it a base equation type
       */
      if (eqn_type == R_MOMENTUM1  || eqn_type == R_MOMENTUM2 ||
	  eqn_type == R_MOMENTUM3  ||
	  eqn_type == R_MOM_NORMAL || eqn_type == R_MOM_TANG1 ||
	  eqn_type == R_MOM_TANG2   ) {
	eqn_type = R_MOMENTUM1;	
      }
      /*
       *  If we are applying a discontinuous bc on one species
       *  equation, then we must apply it to all species equations.
       */
      if (eqn_type == R_MASS ||
	  (eqn_type >= R_SPECIES_UNK_0 && eqn_type <= R_SPECIES_UNK_LAST)
	  ) {
	eqn_type =  R_SPECIES_UNK_0;
      }

      for (ss_index = 0; ss_index < exo->num_side_sets; ss_index++) {
	if (BC_Types[ibc].BC_ID == exo->ss_id[ss_index]) {
	  for (side_index = 0; side_index < exo->ss_num_sides[ss_index];
	       side_index++) {
	    for (k = exo->ss_node_side_index[ss_index][side_index];
		 k < exo->ss_node_side_index[ss_index][side_index+1]; k++) {
	      node_num = exo->ss_node_list[ss_index][k];
              node_ptr = Nodes[node_num];
	      curr_mat_list = &(node_ptr->Mat_List);
	      num_mat = curr_mat_list->Length;

	      /*
	       *  Now make sure that we have the discontinuous var turned
	       *  on
	       */
	      for (imat = 0; imat < num_mat; imat++) 
                 {
		  mat_index= (curr_mat_list->List)[imat];
		  curr_pd = pd_glob[mat_index];
                  for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++)
                     {
                      if (eqn_type == R_MOMENTUM1) 
                        {
		         turn_on_discontinuous(curr_pd, R_MOMENTUM1, imtrx);
		         turn_on_discontinuous(curr_pd, R_MOMENTUM2, imtrx);
		         turn_on_discontinuous(curr_pd, R_MOMENTUM3, imtrx);
		         turn_on_discontinuous(curr_pd, PRESSURE, imtrx);
                        }
		      else if (eqn_type == R_SPECIES_UNK_0) 
                        {
		         turn_on_discontinuous(curr_pd, R_MASS, imtrx);	
		         for (var_type = R_SPECIES_UNK_0;
		              var_type < R_SPECIES_UNK_LAST; var_type++) 
                            {
		             turn_on_discontinuous(curr_pd, var_type, imtrx);
		            }
		        } 
                      else 
                        {
		         turn_on_discontinuous(curr_pd, eqn_type, imtrx);
		        }
                     }
	         }
	    }
	  }
	}
      }
    }
  }

  /*
   *  Just to dot the eyes, make sure that v fields are uniform on
   *  distributed processor problems. We will use the MPI_BOR
   *  operation on a Reduce operation to processor zero, followed
   *  by a broadcast from zero, to accomplish this.
   */
#ifdef PARALLEL
  ivec = alloc_int_1(V_LAST, 0);
  for (imat = 0; imat < upd->Num_Mat; imat++) 
     {
      curr_pd = pd_glob[imat];
      for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++)
         {
          for (k = 0; k <  V_LAST; k++) 
             {
              ivec[k] = curr_pd->v[imtrx][k];
             }
      ReduceBcast_BOR(ivec, V_LAST);
          for (k = 0; k < V_LAST; k++) 
             {
              curr_pd->v[imtrx][k] = ivec[k];
             }
         }    
     }
  safer_free((void **) &ivec);
#endif
  
  return 0;
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

static int
find_next_max(int *bin_list, int *flag_list, int num)

    /********************************************************************
     *
     * find_next_max()
     *
     *     Find the next maximum element in list that is above zero.
     *     Set the flag_list[] for that element so that it won't be
     *     chosen again.
     *     If no elements are above zero, return -1.
     *
     ********************************************************************/
{
  int i, max_id = -1, max_num = 0;
  for (i = 0; i < num; i++) {
    if (bin_list[i] > max_num) {
      if (! flag_list[i]) {
        max_id = i;
	max_num = bin_list[i];
      }
    }
  }
  if (max_id >= 0) {
    flag_list[max_id] = 1;
  }
  return max_id;
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

static int assign_matrl_2(struct Boundary_Condition *bc_ptr, int matrl)
{
  if (bc_ptr->BC_matrl_index_1 == matrl ||
      bc_ptr->BC_matrl_index_2 == matrl) {
    return 0;
  }
  if (bc_ptr->BC_matrl_index_1 < 0) {
    bc_ptr->BC_matrl_index_1 = matrl;
    return 1;
  }
  if (bc_ptr->BC_matrl_index_2 < 0) {
    bc_ptr->BC_matrl_index_2 = matrl;
    return 1;
  }
  return -1;
}
/************************************************************************/
/************************************************************************/
/************************************************************************/
static void
bc_matrl_index_print(struct Boundary_Condition *bc_ptr, int *bin_matrl,
		     int *bin_matrl_elem, int min_node_matrl,
		     int max_node_matrl, int node_matrl_1,
		     int node_matrl_2, int node_matrl_3,
		     int node_matrl_4, int ibc)

    /*******************************************************************
     *
     * bc_matrl_index_print():
     *
     *  Print out a summary of what bc_matrl_index() determined for
     *  the current bc on this processor.
     *******************************************************************/
{
  int i;
  NODE_INFO_STRUCT *node_ptr;
  UMI_LIST_STRUCT *matrlLP;
  printf("=============================================================\n");
  printf("\tbc_matrl_index statistics for ibc = %d on Proc %d:\n",
	 ibc, ProcID);
  printf("\t  Matrl_index    bin_matrl    [bin_matrl_elem]\n");
  printf("\t===================================================\n");
  for (i = 0; i < upd->Num_Mat; i++) {
    printf("\t %6d            %6d ", i, bin_matrl[i]);
    if (!strcmp(bc_ptr->Set_Type, "SS")) {
      printf("      %6d", bin_matrl_elem[i]);
    }
    printf("\n");
  }
  printf("\t===================================================\n");
  if (max_node_matrl >= 0) {
    node_ptr = Nodes[max_node_matrl];
    matrlLP = &(node_ptr->Mat_List);
    printf("max materials at a node (%d) = %d,",
	   max_node_matrl,  matrlLP->Length);
    printf("  their id's = ");
    for (i = 0; i < matrlLP->Length; i++) printf("%d ", matrlLP->List[i]);
    printf("\n");
  
    node_ptr = Nodes[min_node_matrl];
    matrlLP = &(node_ptr->Mat_List);
    printf("min materials at a node (%d) = %d,",
	   min_node_matrl,  matrlLP->Length);
    printf("  their id's = ");
    for (i = 0; i < matrlLP->Length; i++) printf("%d ", matrlLP->List[i]);
    printf("\n");

    printf("\t\tbc->matrl_1 = %d\n", bc_ptr->BC_matrl_index_1);
    if (bc_ptr->BC_matrl_index_2 != -1)
	printf("\t\tbc->matrl_2 = %d\n", bc_ptr->BC_matrl_index_2);
    if (bc_ptr->BC_matrl_index_3 != -1)
	printf("\t\tbc->matrl_3 = %d\n", bc_ptr->BC_matrl_index_3);
    if (bc_ptr->BC_matrl_index_4 != -1)
	printf("\t\tbc->matrl_4 = %d\n", bc_ptr->BC_matrl_index_4);
  }
  printf("=============================================================\n");
  fflush(stdout);
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

void
bc_set_index(Exo_DB *exo)

    /********************************************************************
     *
     * bc_set_index()
     *
     *   Calculates the set index that matches the BC_ID input from the
     * data card. It storres the result in Set_Index for later repeated
     * usage.
     *  MPI USAGE NOTE: Goma does not necessarily supply all side set
     *                  IDs to all processors. If a processor does not
     *                  know about any sides in a side set, then the 
     *                  side set is left out of that processors 
     *                  knowledgebase altogeter. This may or may not be
     *                  the way to go. In any case, for mp problems if
     *                  the side set is not known on a processor, then
     *                  the search here will fail on that processor.
     *                  All we care about is that it successful on one
     *                  processor.
     ********************************************************************/
{
  int ibc, ss_index;
  struct Boundary_Condition *bc_ptr;
  for (ibc = 0; ibc < Num_BC; ibc++) {
    bc_ptr = BC_Types + ibc;
    bc_ptr->Set_Index = -1;
    /*
     * Side Sets
     */
    if (!strcmp(bc_ptr->Set_Type, "SS")) {
      for (ss_index = 0; ss_index < exo->num_side_sets; ss_index++) {
	if (bc_ptr->BC_ID == exo->ss_id[ss_index]) {
          bc_ptr->Set_Index = ss_index;
	  break;
	}
      }
    }
    /*
     * Node Sets
     */
    if (!strcmp(bc_ptr->Set_Type, "NS")) {
      for (ss_index = 0; ss_index < exo->num_node_sets; ss_index++) {
	if (bc_ptr->BC_ID == exo->ns_id[ss_index]) {
	  bc_ptr->Set_Index = ss_index;
	  break;
	}
      }
    }
    /*
     * Level Set Conditions
     * we just need to fake the code below into thinking we identified
     * the right surface to apply this condition to
     */
    if (!strcmp(bc_ptr->Set_Type, "LS") || 
	!strcmp(bc_ptr->Set_Type, "PF")    ) {
      bc_ptr->Set_Index = 0;
    }
    /*
     *  Find the maximum value of the set index across all of the 
     *  processors. All we are interested in is finding out if one
     *  of the processors has found the index.
     */
    ss_index = gmax_int(bc_ptr->Set_Index);
    if (ss_index == -1) {
      DPRINTF(stderr,"On P_%d, bc_set_index: ERROR Can't find ns or ss id %d for"
	      " boundary condition %d with type %s \n",
	      ProcID, bc_ptr->BC_ID, ibc, bc_ptr->Set_Type);
      EH(-1,"Error in specification of boundary condition");
    }
  }
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

void
bc_internal_boundary(Exo_DB *exo)

    /********************************************************************
     *
     * bc_internal_boundary():
     *
     *  Determine if a boundary condition is on an internal boundary
     *  specified by a side set. If it is, then flag the boolean
     *  Internal_Boundary. (Note: most of the work has already been
     *  done in rd_mesh.c )
     ********************************************************************/
{
  int ibc, ss_index;
  struct Boundary_Condition *bc_ptr;
  for (ibc = 0; ibc < Num_BC; ibc++) {
    bc_ptr = BC_Types + ibc;
    bc_ptr->Internal_Boundary = FALSE;
    if (!strcmp(bc_ptr->Set_Type, "SS")) {
      ss_index = bc_ptr->Set_Index;
      if (ss_index >= 0) {
        bc_ptr->Internal_Boundary = SS_Internal_Boundary[ss_index];
      } else {
        bc_ptr->Internal_Boundary = -1;
      }
    }
  }
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

void
bc_matrl_index(Exo_DB *exo)
    
    /********************************************************************
     *
     * bc_matrl_index():
     *
     *    Find out what materials are on each side of a boundary
     *    condition. Note, some boundary conditions require one
     *    to specify this information in the form of an element block
     *    id. However, some others do not. This procedure attempts to
     *    calculate this for all boundary conditions and then fill in
     *    the BC_matrl_index_# elements of the Boundary_Condition
     *    structure.
     *
     *    The basic algorithm involves finding some information about
     *    the side or node set on which the boundary condition is
     *    applied. Then, given the bc name and this information,
     *    a decision is made as to what the BC_matrl_index_#'s should be
     *    set to.
     *
     *    We gather the following information to make this decision:
     *        1) element block indecises input from deck
     *        2) Number of nodes in the bc set in each material
     *        3) Number of elements, whoses sides are in the side set,
     *           in each material
     *        4) Node in the bc set which contains the minimum
     *           number of materials
     *        5) Node in the bc set  which contains the maximum
     *           number of materials
     *        6) A random node in the bc set which contains a
     *           specific number of materials (1, 2, 3, 4)
     *           
     *
     * HKM NOTE:
     *    This routine is a work in progress. Calculation of EDGE's and
     *    VERTICES are not done yet. Also, mp aspects haven't been
     *    figured out yet.
     *
     *******************************************************************/    
{
  int ibc, ss_index, side_index, k, node_num, i;
  int i_apply_meth, num_matrl_needed = -1;
  int found = FALSE, min_node_matrl, max_node_matrl, min_matrl,
      max_matrl, node_matrl_1, node_matrl_2, node_matrl_3,
      node_matrl_4, matrl_first;
  int *bin_matrl, *ind_matrl, *bin_matrl_elem,
      *node_flag_1ss, node_count;
  int mat_index, success, ielem;
  NODE_INFO_STRUCT *node_ptr;
  UMI_LIST_STRUCT *matrlLP;
  struct Boundary_Condition *bc_ptr;
  static char *yo = "bc_matrl_index :";

  bin_matrl = alloc_int_1(upd->Num_Mat * 4,   INT_NOINIT);
  ind_matrl = bin_matrl + upd->Num_Mat;
  bin_matrl_elem = ind_matrl +  upd->Num_Mat;
  node_flag_1ss = alloc_int_1(exo->num_nodes, INT_NOINIT);

  for (ibc = 0; ibc < Num_BC; ibc++) {
    bc_ptr = BC_Types + ibc;
    found = FALSE;
    /*
     *  Some boundary condition specifications already require
     *  you to specify the element blocks on either side of the
     *  side set.
     */
     switch (bc_ptr->BC_Name) {
       /*
	* For these boundary conditions, the element block ID numbers
	* are in the firest two integer slots
	*/
     case POROUS_PRESSURE_BC:
     case DARCY_CONTINUOUS_BC:
     case Y_DISCONTINUOUS_BC:
     case POROUS_GAS_BC:
     case VP_EQUIL_BC:
     case VN_POROUS_BC:
     case FLUID_SOLID_BC:
     case SOLID_FLUID_BC:
     case NO_SLIP_BC:
     case FLUID_SOLID_CONTACT_BC:
     case SOLID_FLUID_CONTACT_BC:
     case T_CONTACT_RESIS_BC:
     case T_CONTACT_RESIS_2_BC:
     case LIGHTP_JUMP_BC:
     case LIGHTM_JUMP_BC:
     case LIGHTP_JUMP_2_BC:
     case LIGHTM_JUMP_2_BC:
	 bc_ptr->BC_matrl_index_1 = map_mat_index(bc_ptr->BC_Data_Int[0]);
	 bc_ptr->BC_matrl_index_2 = map_mat_index(bc_ptr->BC_Data_Int[1]);
	 break;

	 /*
	  * For this boundary condition the element block numbers
	  * are in the second and third integer slots
	  */
     case VL_EQUIL_BC:
     case YFLUX_DISC_RXN_BC:
     case DISCONTINUOUS_VELO_BC:	 
	 bc_ptr->BC_matrl_index_1 = map_mat_index(bc_ptr->BC_Data_Int[1]);
	 bc_ptr->BC_matrl_index_2 = map_mat_index(bc_ptr->BC_Data_Int[2]);
	 break;
     }
     /*
      * Initialize quantities
      */
     for (i = 0; i < 4 * upd->Num_Mat; i++) bin_matrl[i] = 0;
     for (i = 0; i < exo->num_nodes; i++) node_flag_1ss[i] = 0;
     max_node_matrl = min_node_matrl = -1;
     min_matrl = 4000000;
     max_matrl = -1;
     node_matrl_1 = node_matrl_2 = node_matrl_3 = node_matrl_4 = -1;
     /*
      *  Determine how many materials each boundary condition needs
      *  based on the type of the boundary condition
      */
     i_apply_meth = BC_Types[ibc].desc->i_apply;
     if (i_apply_meth == CROSS_PHASE_DISCONTINUOUS ||
	 i_apply_meth == CROSS_PHASE) {
       num_matrl_needed = 2;
     } else if (i_apply_meth == SINGLE_PHASE) {
       num_matrl_needed = 1;
     }
     /*
      *  Loop over the nodes in the side set looking up what
      *  materials are located at each node
      *    -> this will work for side sets. Need to do node sets
      *       as well.
      */
     if (!strcmp(bc_ptr->Set_Type, "SS")) {
       for (ss_index = 0; ss_index < exo->num_side_sets; ss_index++) {
	 /*
	  * This logic works for one side set specifications. For
	  * two side set specifications (EDGES), we will have to go
	  * with a calculation of the union of side sets.
	  */
	 if (bc_ptr->BC_ID == exo->ss_id[ss_index]) {
	   found = TRUE;
	   for (side_index = 0;
		side_index < exo->ss_num_sides[ss_index];
		side_index++) {
	     /*
	      * Locate the element number, find the material index,
	      * then bin the result.
	      */
	     ielem = exo->ss_elem_list[exo->ss_elem_index[ss_index]+side_index];
	     mat_index = find_mat_number(ielem, exo);
	     bin_matrl_elem[mat_index]++;
	   
	     for (k = exo->ss_node_side_index[ss_index][side_index];
		  k < exo->ss_node_side_index[ss_index][side_index+1];
		  k++) {
	       node_num = exo->ss_node_list[ss_index][k];
	       if (!node_flag_1ss[node_num]) {
		 node_flag_1ss[node_num] = 1;
		 node_ptr = Nodes[node_num];
		 matrlLP = &(node_ptr->Mat_List);
		 /*
		  * Bin the materials at this node for later usage.
		  */
	      
		 for (i = 0; i < matrlLP->Length; i++) {
#ifdef DEBUG_IGNORE_ELEMENT_BLOCK_CAPABILITY
                   if (matrlLP->List[i] < 0) {
                     fprintf(stderr,"Material list contains negative number\n");
                     EH(-1,"logic error in ignoring an element block");
                   }
#endif
		   bin_matrl[matrlLP->List[i]]++;
		 }
		 /*
		  *  Find the max and min number of materials for a
		  *  node in this side set
		  */
		 if (matrlLP->Length > max_matrl) {
		   max_matrl = matrlLP->Length;
		   max_node_matrl = node_num;
		 }
		 if (matrlLP->Length < min_matrl) {
		   min_matrl = matrlLP->Length;
		   min_node_matrl = node_num;
		 }
		 /*
		  * Find representative nodes with specific
		  * numbers of materials
		  */
		 if (matrlLP->Length == 1) node_matrl_1 = node_num;
		 if (matrlLP->Length == 2) node_matrl_2 = node_num;
		 if (matrlLP->Length == 3) node_matrl_3 = node_num;
		 if (matrlLP->Length == 4) node_matrl_4 = node_num;
	       }
	     }
	   }
	 }
       } /* End of side set loop */
     } /* if SS */
     /*
      * Node Sets
      */
     if (!strcmp(bc_ptr->Set_Type, "NS")) {
       for (ss_index = 0; ss_index < exo->num_node_sets; ss_index++) {
	 /*
	  * This logic works for one side set specifications. For
	  * two side set specifications (EDGES), we will have to go
	  * with a calculation of the union of side sets.
	  */
	 if (bc_ptr->BC_ID == exo->ns_id[ss_index]) {
	   found = TRUE;

	   /*
	    * Loop over the number of nodes
	    */
	   for (k = 0; k < exo->ns_num_nodes[ss_index]; k++) {
	     node_num = exo->ns_node_list[exo->ns_node_index[ss_index]+k];
	     if (!node_flag_1ss[node_num]) {
	       node_flag_1ss[node_num] = 1;
	       node_ptr = Nodes[node_num];
	       matrlLP = &(node_ptr->Mat_List);
	       /*
		* Bin the materials at this node for later usage.
		*/
	      
	       for (i = 0; i < matrlLP->Length; i++) {
		 bin_matrl[matrlLP->List[i]]++;
	       }
	       /*
		*  Find the max and min number of materials for a
		*  node in this side set
		*/
	       if (matrlLP->Length > max_matrl) {
		 max_matrl = matrlLP->Length;
		 max_node_matrl = node_num;
	       }
	       if (matrlLP->Length < min_matrl) {
		 min_matrl = matrlLP->Length;
		 min_node_matrl = node_num;
	       }
	       /*
		* Find representative nodes with specific
		* numbers of materials
		*/
	       if (matrlLP->Length == 1) node_matrl_1 = node_num;
	       if (matrlLP->Length == 2) node_matrl_2 = node_num;
	       if (matrlLP->Length == 3) node_matrl_3 = node_num;
	       if (matrlLP->Length == 4) node_matrl_4 = node_num;
	     }
	   }	   
	 }
       } /* End of side set loop */
     } /* if NS */     

     /*
      *  Ok, we have obtained statistics on the side and node sets
      *  let's make a decision
      */
     
     if (found) {
       matrl_first = find_next_max(bin_matrl, ind_matrl, upd->Num_Mat);
       success = assign_matrl_2(bc_ptr, matrl_first);
       if (success < 0) {
	 printf("%s P_%d: problem in assigning first matrl index in ibc %d, %d:\n",
		yo, ProcID, ibc, matrl_first);
	 bc_matrl_index_print(bc_ptr, bin_matrl, bin_matrl_elem,
			      min_node_matrl, max_node_matrl, 
			      node_matrl_1, node_matrl_2, node_matrl_3,
			      node_matrl_4, ibc);
       }
       matrl_first = find_next_max(bin_matrl, ind_matrl, upd->Num_Mat);
       if (matrl_first >= 0) {
	 success = assign_matrl_2(bc_ptr, matrl_first);
	 if (success < 0) {
	   printf("%s P_%d: problem in assigning second matrl index in ibc %d, %d:\n",
		  yo, ProcID, ibc, matrl_first);
	   bc_matrl_index_print(bc_ptr, bin_matrl, bin_matrl_elem,
				min_node_matrl, max_node_matrl,
				node_matrl_1, node_matrl_2, node_matrl_3,
				node_matrl_4, ibc);
	 }
       } else {
         if (num_matrl_needed > 1) {
	   printf("%s P_%d: problem in finding a needed second matrl index:\n",
		  yo, ProcID);
	   EH(-1,"bc_matrl_index ERROR");
	 }
       }
     }
     /*
      * For debug purposes, print out everything that we have found
      * out and decided about this bc on all of the processors.
      */

     /*
      * MP Fix: We may not get the same results on all processors
      *         In this case, just take the processor with the most
      *         nodes in this bc and with a valid result, and use
      *         that. Broadcast that result to all nodes. Cross your
      *         fingers.
      */
       
     node_count = 0;
     for (i = 0; i < exo->num_nodes; i++) {
       node_count += node_flag_1ss[i];
     }
#ifdef PARALLEL

     k = ProcWithMaxInt(node_count, &i);
     MPI_Bcast(&(bc_ptr->BC_matrl_index_1), 1, MPI_INT, k,
	       MPI_COMM_WORLD);
     MPI_Bcast(&(bc_ptr->BC_matrl_index_2), 1, MPI_INT, k,
	       MPI_COMM_WORLD);
     MPI_Bcast(&(bc_ptr->BC_matrl_index_3), 1, MPI_INT, k,
	       MPI_COMM_WORLD);
     MPI_Bcast(&(bc_ptr->BC_matrl_index_4), 1, MPI_INT, k,
	       MPI_COMM_WORLD);


#endif
     
  }
  safer_free((void **) &bin_matrl);
  safer_free((void **) &node_flag_1ss);
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

void
determine_dvi_index(void)

    /********************************************************************
     *
     * determine_dvi_index():
     *
     *    
     ********************************************************************/
{
  int ibc;
  struct Boundary_Condition *bc_ptr;
  struct BC_descriptions *bc_desc;
  for (ibc = 0; ibc < Num_BC; ibc++) {
    bc_ptr = BC_Types + ibc;
    bc_desc = bc_ptr->desc;

    /*
     * set the default indexing type for boundary conditions based on
     * the value of the DV_Index_Default field in the BC_description
     * field. They will be overridden for
     * boundary conditions that don't fit the defaults. Generally,
     * the following defaults are valid:
     *
     *  SINGLE_PHASE                 DVI_SINGLE_PHASE_DB
     *  CROSS_PHASE                  DVI_CROSS_PHASE_CONJUGATE
     *  CROSS_PHASE_DISCONTINUOUS    DVI_SIDTIE
     *
     *  Set the default value of DV_Indexing_MatID to the first 
     *  material index (ie., lowest valued index if there are more
     *  than one.
     */
    bc_ptr->DV_Indexing_Type  = bc_desc->DV_Index_Default;
    bc_ptr->DV_Indexing_MatID = bc_ptr->BC_matrl_index_1;

    /*
     * Override the defaults. Also check to see if there is a
     * conflict with the type of boundary condition and whether
     * or not the current side set is an internal boundary or
     * not.
     */
    switch (bc_ptr->BC_Name) {

    case CAP_REPULSE_BC:
    case CAPILLARY_BC:
    case CAPILLARY_TABLE_BC:
    case CAP_REPULSE_ROLL_BC:
    case CAP_REPULSE_USER_BC:
    case CAP_REPULSE_TABLE_BC:
    case CAPILLARY_SHEAR_VISC_BC:
    case LATENT_HEAT_BC:
    case YFLUX_USER_BC:
    case SDC_HEATRXN_BC:
	if (bc_ptr->Internal_Boundary) {
	  bc_ptr->DV_Indexing_Type = DVI_DVVSIG;
	} else {
	  bc_ptr->DV_Indexing_Type = DVI_VSIG;
	}
	break;

    case FLUID_SOLID_BC:
    case FLUID_SOLID_RS_BC:
    case SOLID_FLUID_BC:
    case SOLID_FLUID_RS_BC:
    case LATENT_HEAT_INTERNAL_BC:
    case POROUS_CONV_BC:
    case POROUS_GAS_BC:
        if (! bc_ptr->Internal_Boundary) {
          fprintf(stderr,
		  "Boundary condition %s can't be put on external domain bc\n",
		  bc_desc->name1);
	  EH(-1,"Boundary Condition incompatibility");
	}
	bc_ptr->DV_Indexing_Type = DVI_VSIG;
	break;

    case POROUS_FLUX_BC:
	if (bc_ptr->Internal_Boundary) {
          fprintf(stderr,
		  " WARNING: Boundary condition %s can't be put on internal domain bc and expect mass conservation\n",
		  bc_desc->name1);
	  fprintf(stderr,
		  "\t No attempt at matching fluxes on both sides of interfaces is made\n");
	}
	break;

    case KINEMATIC_DISC_BC:
        if (! bc_ptr->Internal_Boundary) {
          fprintf(stderr,
		  "Boundary condition %s can't be put on external domain bc\n",
		  bc_desc->name1);
	  EH(-1,"Boundary Condition incompatibility");
	}
	bc_ptr->DV_Indexing_Type = DVI_MULTI_PHASE_SINGLE;
        break;

    case KINEMATIC_SPECIES_BC:
        if (! bc_ptr->Internal_Boundary) {
          fprintf(stderr,
		  "Boundary condition %s can't be put on external domain bc\n",
		  bc_desc->name1);
	  EH(-1,"Boundary Condition incompatibility");
	}
	bc_ptr->DV_Indexing_Type = DVI_DVVSIG;
	break;

    case VELO_NORMAL_EDGE_BC:
    case VELO_NORMAL_EDGE_INT_BC:
    case VELO_TANGENT_BC:
    case VELO_STREAMING_BC:
    case VELO_TANGENT_USER_BC:
    case VELO_TANGENT_3D_BC: 
    case VELO_TANGENT_EDGE_BC:
    case VELO_TANGENT_EDGE_INT_BC:
    case VELO_SLIP_BC:
    case VELO_SLIP_FILL_BC:
    case VELO_SLIP_ROT_BC:
    case VELO_SLIP_ROT_FILL_BC:
    case VELO_SLIP_FLUID_BC:
    case VELO_SLIP_ROT_FLUID_BC:
    case AIR_FILM_BC:
    case AIR_FILM_ROT_BC:
	bc_ptr->DV_Indexing_Type = DVI_SID;
	break;

    case VL_EQUIL_PRXN_BC:
    case YFLUX_DISC_RXN_BC:
    case SDC_SURFRXN_BC:
        if (! bc_ptr->Internal_Boundary) {
          fprintf(stderr,
		  "Boundary condition %s can't be put on external domain bc\n",
		  bc_desc->name1);
	  EH(-1,"Boundary Condition incompatibility");
	}
	bc_ptr->DV_Indexing_Type = DVI_MULTI_PHASE_VD;
	break;

    case VL_EQUIL_BC:
	if (! bc_ptr->Internal_Boundary) {
          fprintf(stderr,
		  "Boundary condition %s can't be put on external domain bc\n",
		  bc_desc->name1);
	  EH(-1,"Boundary Condition incompatibility");
	}

    }

  }
}
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void
setup_external_nodal_matrls(Exo_DB *exo, Dpi *dpi, Comm_Ex *cx)

    /************************************************************************
     *
     * setup_external_nodal_matrls():
     *
     *   This routine exchanges information about the materials
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
  int i, k, p, node_num, max_Matrl;
  int *mesg_send = NULL, *mesg_recv = NULL, *istart;
  NODE_INFO_STRUCT *node;
  COMM_NP_STRUCT *np_ptr, *np_base = NULL;

  /*
   * Find out the maximum number of materials per node
   */
  max_Matrl = find_MaxMatrlPerNode();

  /*
   *  Pack information for sending
   */

  if (dpi->num_neighbors > 0) {
    mesg_send = alloc_int_1(ptr_node_send[dpi->num_neighbors] * max_Matrl, -1);
    mesg_recv = alloc_int_1(ptr_node_recv[dpi->num_neighbors] * max_Matrl, -1);
 
    for (i = 0; i < ptr_node_send[dpi->num_neighbors]; i++) {
      node_num = list_node_send[i];
      node = Nodes[node_num];
      istart = mesg_send + i*max_Matrl;
      for (k = 0; k < node->Mat_List.Length; k++) {
	istart[k] = node->Mat_List.List[k];
      }
    }
    np_base = alloc_struct_1(COMM_NP_STRUCT, dpi->num_neighbors);
  }
  
  np_ptr = np_base;
  for (p = 0; p < dpi->num_neighbors; p++) {
    np_ptr->neighbor_ProcID = cx[p].neighbor_name;
    np_ptr->send_message_buf =
	(void *) (mesg_send + ptr_node_send[p] * max_Matrl);
    np_ptr->send_message_length =
	sizeof(int) * cx[p].num_nodes_send * max_Matrl;
    np_ptr->recv_message_buf =
	(void *) (mesg_recv + ptr_node_recv[p] * max_Matrl);
    np_ptr->recv_message_length =
	sizeof(int) * cx[p].num_nodes_recv * max_Matrl;
    np_ptr++;
  }
  exchange_neighbor_proc_info(dpi->num_neighbors, np_base);
  

  for (i = 0; i < dpi->num_external_nodes; i++) {
    node_num = dpi->num_internal_nodes + dpi->num_boundary_nodes + i;
    node = Nodes[node_num];
    istart = mesg_recv + i*max_Matrl;


    /*
     * Unpack ther materials and then add them to the existing
     * list.
     */
    for (k = 0; k < max_Matrl; k++) {
      if (istart[k] != -1) {
	add_to_umi_int_list(&(node->Mat_List), istart[k]);
      } else {
	break;
      }
    }   
  }

  /*
   *  Free memory allocated in this routine
   */
  safer_free((void **) &np_base);
  safer_free((void **) &mesg_send);
  safer_free((void **) &mesg_recv);

  /*
   *  When in debug mode, print out a complete listing of variables at
   *  every node
   */

}
/***************************************************************************/
/***************************************************************************/
/***************************************************************************/

int
find_MaxMatrlPerNode(void)

     /***************************************************************
      *
      * find_MaxMatrlPerNode:
      *
      *    This function returns the maximum number of materials present
      * at any one node in the domain.
      ***************************************************************/
{
  int number_proc_nodes = EXO_ptr->num_nodes;
  int I, lmax = 0;
  for (I = 0; I < number_proc_nodes; I++) {
    lmax = MAX(lmax, Nodes[I]->Mat_List.Length);
  }
  lmax = gmax_int(lmax);
  return lmax;
}
/***************************************************************************/
/***************************************************************************/
/***************************************************************************/
