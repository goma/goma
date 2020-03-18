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
 *$Id: mm_bc.c,v 5.11 2010-07-21 16:39:26 hkmoffa Exp $
 */


#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "std.h"
#include "rf_fem_const.h"
#include "rf_fem.h"
#include "rf_io_const.h"
#include "rf_io.h"
#include "rf_mp.h"
#include "el_elm.h"
#include "el_geom.h"
#include "rf_allo.h"
#include "rf_bc.h"
#include "mm_names.h"
#include "mm_eh.h"
#include "dpi.h"
#include "rf_vars_const.h"
#include "mm_as_structs.h"
#include "mm_as.h"
#include "mm_as_const.h"
#include "bc_colloc.h"
#include "el_elm_info.h"
#include "exo_struct.h"
#include "mm_bc.h"
#include "mm_elem_block_structs.h"
#include "mm_ns_bc.h"
#include "mm_post_proc.h"
#include "mm_unknown_map.h"
#include "rd_mesh.h"
#include "rf_bc_const.h"
#include "rf_node_const.h"
#include "rf_shape.h"
#include "stdbool.h"

#define GOMA_MM_BC_C

/*
 * These two workhorse variables help get more specific information out
 * to the generic error handler. They're declared here so they may be used
 * in every routine in this file without worrying about local declaration.
 */

static Spfrtn sr;

/*
 * Variable DEFINITIONS...
 */

int **ROT_list = NULL;

int **SS_list = NULL;

int Use_2D_Rotation_Vectors=FALSE;

int PRESSURE_DATUM = FALSE; /* flag to determine if a pressure datum is set */
int pressure_datum_element = 0; /* element in which the pressure datum is set */
double pressure_datum_value = 0.0; /* value of the pressure datum */



/********** R O U T I N E S   D E F I N E D   I N   T H I S   F I L E **********
 *
 *     NAME OF FUNCTION			TYPE 		CALLED BY
 * ------------------------       ---------------     -------------------
 *   find_and_set_Dirichlet ()		void	      rf_sol_nonlinear
 *       set_nodal_Dirichlet_BC() static void	      find_and_set_Dirichlet
 *   alloc_First_Elem_BC()	 	void	      pre_process
 *   set_up_Surf_BC ()		 	void	      rf_sol_nonlinear
 *       setup_Elem_BC ()	 static void	      set_up_Surf_BC
 *       same_side ()	 	 static  int	      set_up_Surf_BC
 *   print_setup_Surf_BC ()	        void	      rf_sol_nonlinear
 *       print_Elem_Surf_BC ()    static void          print_setup_Surf_BC
 *
 *******************************************************************************/

/** P R O T O   D E F I N I T I O N S   O  F   S T AT I C   F U N C T IO N S **/

static int 
set_nodal_Dirichlet_BC ( 
			       int 		        inode,
			       int                        ibc,
			       struct Boundary_Condition *boundary_condition,
			       double		        x[],
			       double                     xdot[],
			       int
			       );

static void
setup_Elem_BC (  struct elem_side_bc_struct **elem_side_bc, 
    			struct Boundary_Condition  *bc_type, 
    			int ibc, int num_nodes_on_side,  int ielem, int *,
			Exo_DB *exo);

static struct elem_edge_bc_struct *setup_Elem_Edge_BC	/* mm_bc.c */
(struct elem_edge_bc_struct **,			/* elem_edge_bc */
	struct Boundary_Condition  *,			/* bc_type */
	int ,						/* ibc */
	int ,						/* num_nodes_on_edge */
	int ,						/* ipin */
	int ,                                           /* shared */
	int ,						/* ielem */
	int [],						/* local_edge_node_list
							 */
	Exo_DB *);					/* exo */

static int 
same_side     (  int [], int [], int );

static void 
print_Elem_Surf_BC ( int ielem, struct elem_side_bc_struct *elem_side_bc
			    );
static void elem_side_matrl_list(struct elem_side_bc_struct *);
static void vcrr_determination(struct elem_side_bc_struct *,
			       struct Boundary_Condition *, int, int);

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void
find_and_set_Dirichlet(double x[],    /* solution vector at this processor */
		       double xdot[], /* time derivative of solution vector */
		       Exo_DB *exo,   /* ptr to EXODUS II FE database */
		       Dpi *dpi)      /* ptr to distrib processing FE database */
     /*
      *      Function which finds the nodes in the node-sets and side-sets which have
      *      Dirichlet boundary conditions set.  It then fills in the values of the
      *      Dirichlet boundary condtions into the solution vector, x[].
      *     
      *      Author:          Scott Hutchinson (1421)
      *      Date:            22 January 1993
      *      Revised:         22 January 1993
      *      
      *      ---------------------------------------------------------------------
      *
      *	Migrating towards the new means of accessing EXODUS II finite element
      *	database information. [1997/08/23 10:53 MDT pasacki@sandia.gov]
      *      
      */
{
  int side_index, I, mn, ielem, e_start, e_end, eb_index, ielem_type,
    num_local_nodes, matID, elem_block_index, offset;
  int *been_there = NULL;
  int lni;			/* local node index on one side of a SS */
 
  static const char yo[] = "find_and_set_Dirichlet";

#ifdef DEBUG_BC
  extern int 			ProcID;
#endif /* DEBUG_BC */
  int apply_pressure_datum;	/* boolean */
  int i, ibc, ins, inode = -1, error_cond = FALSE;
  int num_nodes, n, iconnect_ptr, datum_set, ie, local_elem;
  NODE_INFO_STRUCT *node;
  NODAL_VARS_STRUCT *nv;
  VARIABLE_DESCRIPTION_STRUCT *vd;

  /* Loop over Num_BC, the number of boundary conditions defined in the 
     input file */

  for (ibc = 0; ibc < Num_BC; ibc++) {

    /*****************************************************************************/
    /*                             BLOCK 1                                       */
    /*         DIRICHLET BOUNDARY CONDITIONS SPECIFIED BY NODE SETS              */
    /*****************************************************************************/

    /* First check to see if this boundary condition is a node set */
    if (!strcmp(BC_Types[ibc].Set_Type, "NS")) {

      /*
       *   Next, check to see if this node set boundary condition is a Dirichlet
       *  boundary condition
       */
      if (BC_Types[ibc].desc->method == DIRICHLET ||
	  BC_Types[ibc].BC_Name == SURFTANG_BC ||
	  BC_Types[ibc].BC_Name == SURFTANG_SCALAR_BC ||
	  BC_Types[ibc].BC_Name == CAP_ENDFORCE_BC ||
	  BC_Types[ibc].BC_Name == CAP_ENDFORCE_SCALAR_BC ||
	  BC_Types[ibc].BC_Name == CA_BC ||
	  BC_Types[ibc].BC_Name == CA_MOMENTUM_BC ||
	  BC_Types[ibc].BC_Name == CA_OR_FIX_BC ||
	  BC_Types[ibc].BC_Name == MOVING_CA_BC ||
	  BC_Types[ibc].BC_Name == VELO_THETA_TPL_BC ||
	  BC_Types[ibc].BC_Name == VELO_THETA_HOFFMAN_BC ||
	  BC_Types[ibc].BC_Name == VELO_THETA_COX_BC ||
	  BC_Types[ibc].BC_Name == VELO_THETA_SHIK_BC ||
	  BC_Types[ibc].BC_Name == SH_SLOPE_X_BC ||
	  BC_Types[ibc].BC_Name == SH_SLOPE_Y_BC  ||
	  BC_Types[ibc].BC_Name == SHEET_ENDSLOPE_BC ) {

	/*
	 * Resolve what material's variables this boundary condition will
	 * be applied on. If no specification of element block ID is made,
	 * set matID = -2, so that set_nodal_Dirichlet_BC() will assign the
	 * Dirichlet condition to the first variable of a variable type
	 * found at a node. If the element block ID is specified, then
	 * assign matID to the corresponding material index, so that
	 * only the variable associated with that particular material is
	 * assigned to the Dirichlet condition.  
	 */
	if (BC_Types[ibc].BC_EBID_Apply == -1) {
	  matID = -2;
	} else {
	  matID = map_mat_index(BC_Types[ibc].BC_EBID_Apply);
	}

	/* Loop over the total number of node sets defined for the current processor */
  	for (ins = 0; ins < exo->num_node_sets; ins++) {

	  /* Check for a match between the ID of the current node set and the node set
	     ID specified in the input file - continue if a match is found */
	  if (exo->ns_id[ins] == BC_Types[ibc].BC_ID) {

            /* Loop over all the global nodes in the current node set */
	    for (i = 0; i < exo->ns_num_nodes[ins]; i++) {

              /* Get the ith local node in the current node set */
	      inode = exo->ns_node_list[ exo->ns_node_index[ins] + i ];

	      /* Set the bit field in Variable mask to indicate the presence
	       * of a dirichlet BC and also set the boundary condition in x[],
	       * for the current node, inode, and boundary condition
	       */
	      error_cond = 
		(set_nodal_Dirichlet_BC(inode, ibc, &BC_Types[ibc], x, xdot,
					matID)
		 || error_cond);

	    }  /* END for (i = 0; i < Proc_NS_Count[ins]; i++)               */
	  }  /* END if (Proc_NS_Ids[ins] == BC_Types[ibc].BC_ID)	     */
	}  /* END for (ins = 0; ins < Proc_Num_Node_Sets; ins++) 	     */
      }  /* END if (BC_Types[ibc].BC_Name == U_BC  || ....		     */
    }  /* END if (!strcmp(BC_Types[ibc].Set_Type, "NS")) 		     */

    /*****************************************************************************/
    /*                             BLOCK 2                                       */
    /*         DIRICHLET BOUNDARY CONDITIONS SPECIFIED BY SIDE SETS              */
    /*****************************************************************************/

    /* First check to see if this boundary condition is a side set 	      */
    if (!strcmp(BC_Types[ibc].Set_Type, "SS")) {

      /* Next, check to see if this boundary condition is a Dirichlet
       * boundary condition specified by a side set 
       */
      if (BC_Types[ibc].desc->method == DIRICHLET) {
	fprintf(stderr, 
		"\nWARNING: using side-sets for Dirichlet BC %s", 
		BC_Types[ibc].desc->name1);
	/* Loop over the total number of side sets defined for the current processor */
	for (ins = 0; ins < exo->num_side_sets; ins++) {
	  /*
	   * Check for a match between the ID of the current side set and the side set
	   * ID specified in the input file - continue if a match is found 
	   */
	  if (Proc_SS_Ids[ins] == BC_Types[ibc].BC_ID) {
	    /* 
	     * Loop over all the sides in the current side set 
	     */
	    for (side_index = 0; side_index < exo->ss_num_sides[ins]; side_index++) {
	      /*
	       * Discover what element and what material index this side refers
	       * to. If the boundary condition is restricted to being applied to
	       * a particular element block ID, then check this condition. If
	       * condition fails, go on to the next side set. If it passes, pass
	       * the specific matID to the  set_nodal_Dirichlet_BC() so that
	       * the specific variable corresponding to that material can be
	       * flagged for the Dirichlet condition. If the boundary condition
	       * is unspecific to Element Block ID, pass matID to 
	       * set_nodal_Dirichlet_BC() as well. We know the material that
	       * we are currently in. Therefore, we should use that information
	       * to specify the pertinent variable wherever possible.
	       */
	      ielem = exo->ss_elem_list[exo->ss_elem_index[ins]+side_index];
	      elem_block_index = find_elemblock_index(ielem, exo);
	      matID = Matilda[elem_block_index];
	      if (BC_Types[ibc].BC_EBID_Apply != -1) {
		if (Element_Blocks[elem_block_index].Elem_Blk_Id !=
		    BC_Types[ibc].BC_EBID_Apply) {
		  continue;
		}
	      }
	      /*
	       * Loop over the nodes on the current side
	       */
	      for (lni = exo->ss_node_side_index[ins][side_index];
		   lni < exo->ss_node_side_index[ins][side_index+1]; lni++) {
		inode = exo->ss_node_list[ins][lni];
		/*
		 * Set the bit field in Variable mask to indicate the presence
		 *  of a dirichlet BC and also set the boundary condition in x[],
		 * for the current node, inode, and boundary condition	
		 */
		error_cond |= set_nodal_Dirichlet_BC(inode, ibc, &BC_Types[ibc],
						     x, xdot, matID);
	      }
	    }
	  }  /* END if (Proc_SS_Ids[ins] == BC_Types[ibc].BC_ID)	      */
	}  /* END for (ins = 0; ins < Proc_Num_Side_Sets; ins++) 	      */
      } /* END if (BC_Types[ibc].BC_Name ==USIDE_BC ...			      */
    } /* END if (!strcmp(BC_Types[ibc].Set_Type, "SS")) 		      */

    /*****************************************************************************/
    /*                             BLOCK 3                                       */
    /*         DIRICHLET BOUNDARY CONDITIONS SPECIFIED BY MATERIAL NUMBERS       */
    /*****************************************************************************/
    /*
     * First check to see if this boundary condition is volumetric bc
     * based on material numbers
     */
    if (!strcmp(BC_Types[ibc].Set_Type, "MN")) {
      /*
       *  Malloc a temporary vector of ints, that will be used
       *  to ensure that we vist each node that needs a bc, once
       *  and only once.
       */
      been_there = alloc_int_1(Num_Node, 0);
      /*
       *  Check to see if this is indeed a Dirichlet condition
       */
      if (BC_Types[ibc].desc->method == DIRICHLET) {
        /*
	 * Loop over all of the element blocks in the problem
         */
	for (eb_index = 0; eb_index < exo->num_elem_blocks; eb_index++)  {
	  /*
	   *  Identify the materials id with the element block
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
	   *  based on the element type, get the number of local nodes
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
	      if (! been_there[I]) {
		been_there[I] = TRUE;
		/*
		 * Set the bit field in Variable mask to indicate the presence
		 * of a Dirichlet BC and also set the boundary condition in x[],
		 * for the current node, inode, and boundary condition
		 */
		error_cond |= 
		  set_nodal_Dirichlet_BC(inode, ibc, &BC_Types[ibc], x, 
					 xdot, mn);
	      }
	    }
	  } /*  for (ielem = e_start; ielem < e_end; ielem++) */
	} /* for (eb_index = 0; eb_index < exo->num_elem_blocks; eb_index++)  */
      } /* if (BC_Types[ibc].desc->method == DIRICHLET) */
      safer_free((void **) &been_there);
    } /* if (!strcmp(BC_Types[ibc].Set_Type, "MN")) */
  }  /* END for (ibc = 0; ibc < Num_BC; ibc++) */

  if (error_cond) {
    printf("find_and_set_Dirichlet WARNING: error condition found\n");
    printf("\t Will continue on anyway against my better judgement\n");
  }

  /*
   * Modify this so that only the processor owning the pressure datum
   * element attempts to set the datum value. Other processors should
   * not attempt to do anything here.
   */

  apply_pressure_datum = FALSE;	/* standard default to begin... */

  if ( PRESSURE_DATUM == 1 )
    {
      if ( pressure_datum_element < 0 || 
	   pressure_datum_element >= dpi->num_elems_global )
	{
	  log_err("Specified pressure datum element %d out of range.",
		  pressure_datum_element);
	}

      local_elem = in_list(pressure_datum_element, 0, dpi->num_elems, 
			   dpi->elem_index_global);

      /*
       * Still need to check whether the element is actually owned by
       * this processor. That's because our nodal dominated decomposition
       * results in some overassembly of elements that are shared by more
       * than one processor. Therefore, we'll check the unique element
       * decomposition information in our copy of the element ownership
       * assignment for all of the elements this processor traverses.
       */

      if ( local_elem > -1 )
	{
	  if ( dpi->elem_owner[local_elem] == ProcID )
	    {
	      apply_pressure_datum   = TRUE;

	      /*
	       * Reset so the index of the pressure datum element is local
	       * instead of global.
	       */

	      pressure_datum_element = local_elem;
	    }
	}
    }

  if (apply_pressure_datum) {
    iconnect_ptr = Proc_Connect_Ptr[pressure_datum_element];
    /*
     * determine the node at which to apply the pressure datum
     * it is always applied to the first pressure unknown in the designated 
     * element
     *
     * HKM -> It would be best to use a new boundary condition here
     *        so that PRESSURE DATUM could be treated the same as
     *        all of the other boundary conditions.
     */
    num_nodes = elem_info(NNODES, Elem_Type(exo, pressure_datum_element));
    datum_set = 1;
    mn = find_mat_number(pressure_datum_element, exo);
    for (n = 0; n < num_nodes; n++) {
      inode = Proc_Elem_Connect[iconnect_ptr + n];
      node = Nodes[inode];
      nv = node->Nodal_Vars_Info[pg->imtrx];
      if (Dolphin[pg->imtrx][inode][PRESSURE] > 0 && datum_set) {
	datum_set = 0;
	if (!node->DBC[pg->imtrx]) {
	  node->DBC[pg->imtrx] = alloc_short_1(nv->Num_Unknowns, -1);
	}
	offset = get_nodal_unknown_offset(nv, PRESSURE, mn, 0, &vd);

	node->DBC[pg->imtrx][offset] = 0;
	ie = Index_Solution(inode, PRESSURE, 0, 0, mn, pg->imtrx);
	x[ie]  = pressure_datum_value;
	xdot[ie]  = 0.0; 
	log_msg("Setting pressure datum");
	log_msg("        pressure datum element (local)  = %d", 
		pressure_datum_element);
	log_msg("        pressure datum element (global) = %d", 
		dpi->elem_index_global[pressure_datum_element]);
	log_msg("        pressure datum node (local)     = %d", 
		inode);
	log_msg("        pressure datum node (global)    = %d", 
		dpi->node_index_global[inode]);
	log_msg("        pressure datum dof (local)      = %d", 
		ie);
	log_msg("        pressure datum value            = %g", 
		pressure_datum_value);
      }
    }
  }
  return;
} /* END of ROUTINE find_and_set_Dirichlet ***********************************/
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static int
set_nodal_Dirichlet_BC(int inode, int ibc,
		       struct Boundary_Condition *boundary_condition,
		       double x[], double xdot[], int matID)

     /********************************************************************
      *    set_nodal_Dirichlet_BC: 
      *
      *	This function returns a 0 if the boundary condition is 
      *  set successfully.  It returns a 1, if the dirichlet condition has
      *  already been set before.  In this case, the dirichlet condition is not
      *  set again!
      *
      ********************************************************************/
{
  int eqn, w, ieqn, offset;
  int ndof=0;
  NODE_INFO_STRUCT *node = Nodes[inode];
  NODAL_VARS_STRUCT *nv = node->Nodal_Vars_Info[pg->imtrx];
  VARIABLE_DESCRIPTION_STRUCT *vd;
  /*
   * Fill in the correct bit field in Variable_Mask depending
   * on the name of the current Node set, boundary_condition->BC_Name
   */
  if (boundary_condition->desc->method == DIRICHLET 
      && boundary_condition->BC_Name != FIX_BC
      && boundary_condition->BC_Name != PERIODIC_BC) {
    
    eqn = boundary_condition->desc->equation; 
    if (eqn > V_LAST) {
      EH(GOMA_ERROR,"Bad eqn");
    }
    if (eqn != R_MASS) {
      offset = get_nodal_unknown_offset(nv, eqn, matID, 0, &vd);
      if (offset >= 0) {
	if (node->DBC[pg->imtrx] == NULL) {
	  node->DBC[pg->imtrx] = alloc_short_1(nv->Num_Unknowns, -1);
	}
	node->DBC[pg->imtrx][offset] = (short int) ibc;
	ieqn = node->First_Unknown[pg->imtrx] + offset;
        /*
         *  If BC_relax is set to default, set the boundary
         *  condition here without including it in the residual.
         *  We also need to set the time derivative, as it can't
         *  be determined via the normal procedure for this
         *  case.
         *  We will then put a 1 on the diagonal of the
         *  matrix and a zero in the residual corresponding
         *  to this unknown.
         */
	if (boundary_condition->BC_Name == DX_USER_NODE_BC ||
	    boundary_condition->BC_Name == DY_USER_NODE_BC ||
	    boundary_condition->BC_Name == DZ_USER_NODE_BC )
           {
              double epsilon, X1, X2, DX1, DX2;

              X1  = boundary_condition->BC_Data_Float[0];
              DX1 = boundary_condition->BC_Data_Float[1];
              X2  = boundary_condition->BC_Data_Float[2];
              DX2 = boundary_condition->BC_Data_Float[3];

              epsilon=(DX2-DX1)/(X2-X1);
              x[ieqn] = epsilon*(Coor[0][inode]-X1)+DX1;
	      if (boundary_condition->BC_Name == DY_USER_NODE_BC) 
                       {x[ieqn] = epsilon*(Coor[1][inode]-X1)+DX1;}
	      if (boundary_condition->BC_Name == DZ_USER_NODE_BC) 
                       {x[ieqn] = epsilon*(Coor[2][inode]-X1)+DX1;}
              xdot[ieqn] = 0.0;
           }
        else if (boundary_condition->BC_relax == -1.0 ) 
        {
	  if (boundary_condition->BC_Name == F_DIODE_BC)
	    {
	      EH(GOMA_ERROR, "F_DIODE_BC has special needs. You need to comment this line out and follow instructions");
	      /*real hack here for diode level-set.  Let liquid out, but not back in.  
	        if you want to use this comment out the EH above and go to bc_dirich.c and uncomment the if{} protection
                clearly marked for this BC */
	      if(x[ieqn] <= 0.0) 
		{
		  x[ieqn] = boundary_condition->BC_Data_Float[0];
		  xdot[ieqn] = 0.0;
		}
	    } 
	  else
	    {
	      x[ieqn] = boundary_condition->BC_Data_Float[0];
	      xdot[ieqn] = 0.0;
	    }
	}
      }
    } else {
      /*
       *  MASS_FRACTION SPECIAL SECTION
       * Obtain the species number from the first integer
       */
      w = boundary_condition->BC_Data_Int[0];
      if (w >= (upd->Max_Num_Species)) {
	EH(GOMA_ERROR, "Species number on BC Y card exceeds number of species available");
      }
      offset = get_nodal_unknown_offset(nv, eqn, matID, w, &vd);
      if (offset < 0) {
	WH(-1, "Y BC applied to material without species equation");
      } else {
	if (node->DBC[pg->imtrx] == NULL) {
	  node->DBC[pg->imtrx] = alloc_short_1(nv->Num_Unknowns, -1);
	}
	node->DBC[pg->imtrx][offset] = (short int) ibc;
	ieqn = node->First_Unknown[pg->imtrx] + offset;	   
        /*
         *  If BC_relax is set to default, set the boundary
         *  condition here without including it in the residual
         *  We also need to set the time derivative, as it can't
         *  be determined via the normal procedure for this
         *  case.
         */
	if (boundary_condition->BC_relax == -1.0) {
	  x[ieqn] = boundary_condition->BC_Data_Float[0];
	  xdot[ieqn] = 0.0;
	}
      }
    }
  } else if ((boundary_condition->BC_Name == CA_BC) || 
	     (boundary_condition->BC_Name == CA_MOMENTUM_BC) ||
	     (boundary_condition->BC_Name == MOVING_CA_BC) ||
	     (boundary_condition->BC_Name == VELO_THETA_TPL_BC) ||
	     (boundary_condition->BC_Name == VELO_THETA_HOFFMAN_BC) ||
	     (boundary_condition->BC_Name == VELO_THETA_SHIK_BC) ||
	     (boundary_condition->BC_Name == VELO_THETA_COX_BC)) {
    /*
     *  Here we need to be able to cross reference the node flagged
     *  to have a contact angle condition with the CA BC card to which
     *  it corresponds.  Since the BC_Data_Int field of the 
     *  struct boundary_condition is not used for this card, we can place
     *  the global node number in it.
     */
    boundary_condition->BC_Data_Int[0] = 1000000*inode+1000000;
    node->DBCA = 1;

  } else if (boundary_condition->BC_Name == SURFTANG_SCALAR_BC ||
	     boundary_condition->BC_Name == CAP_ENDFORCE_SCALAR_BC ) {
    /*
     *  Here we need to be able to cross reference the node flagged
     *  to have a surface tangent condition with the SURFTANG card to which
     *  it corresponds.  Since the BC_Data_Int field of the 
     *  struct boundary_condition is not used for this card, we can place
     *  the global node number in it.
     */
    boundary_condition->BC_Data_Int[0] = inode;
    node->DBSTS = 1;

  } else if (boundary_condition->BC_Name == SURFTANG_BC || 
	     boundary_condition->BC_Name == CAP_ENDFORCE_BC ) {
    /*
     *  Here we need to be able to cross reference the node flagged
     *  to have a surface tangent condition with the SURFTANG card to which
     *  it corresponds.  Since the BC_Data_Int field of the 
     *  struct boundary_condition is not used for this card, we can place
     *  the global node number in it.
     */
    boundary_condition->BC_Data_Int[0] = inode;
    node->DBST = 1;

  } else if (boundary_condition->BC_Name == CA_OR_FIX_BC) {	
    /*  This is a special condition which takes one of two paths
     *  on startup:  The first is to apply the specified contact
     *  angle if the contact line position is significantly far from
     *  a sharp location point, the location of which can be entered
     *  on the card. 
     *  The criterion for detachment is the Gibbs inequality condition
     *
     *  Here we need to be able to cross reference the node flagged
     *  to have a contact angle condition with the CA BC card to which
     *  it corresponds.  Since the BC_Data_Int field of the 
     *  struct boundary_condition is not used for this card, we can place
     *  the global node number in it.  This will be BC_Data_Int[0]
     *
     *  We also need to keep track of whether this node has been pinned previously
     *  or is still allowed to satisfy a contact angle condition.  We will use
     *  the variable ipin for this within the routines, but the space for it will be
     *  BC_Data_Int[1].  Initialize it to 0, for the pinned case, so that way if you
     *  are restarting with it pinned, then it will immediately evaluate the gibbs
     *  criterion.  If it fails the pinned test, then it will free the line.
     */  
    boundary_condition->BC_Data_Int[1] = 1;
    /*
     * Here we will not make the decision (or maybe) on the pinned versus
     * moving mode.  See the routine mm_fill.c
     *
     *    default case, free contact line:
     */
    boundary_condition->BC_Data_Int[0] = 1000000*inode+1000000;
    node->DBCA = 1;
    /* 
     * Encode original position for reference, using the last available chunks
     * of the bc input float data
     */
    ieqn = Index_Solution(inode, MESH_DISPLACEMENT1, 0, ndof, matID, pg->imtrx);
    boundary_condition->BC_Data_Float[7] =
      boundary_condition->BC_Data_Float[4] - (Coor[0][inode]);
    if (boundary_condition->BC_Data_Int[0]) {
      xdot[ieqn] = 0.;
    }
    ieqn = Index_Solution(inode, MESH_DISPLACEMENT2, 0, ndof, matID, pg->imtrx);
    boundary_condition->BC_Data_Float[8] =
      boundary_condition->BC_Data_Float[5] - (Coor[1][inode]);
    if (boundary_condition->BC_Data_Int[0]) {
      xdot[ieqn] = 0.;
    }
    boundary_condition->BC_Data_Float[9] = 0.;
    if (pd_glob[0]->Num_Dim == 3) {
      ieqn = Index_Solution(inode, MESH_DISPLACEMENT3, 0, ndof, matID, pg->imtrx);
      boundary_condition->BC_Data_Float[9] =
	boundary_condition->BC_Data_Float[6] - (Coor[2][inode]);
      if (boundary_condition->BC_Data_Int[0]) {
	xdot[ieqn] = 0.;
      }
    }
  } else if (boundary_condition->BC_Name == FIX_BC) {
    /*
     * This boundary condition is like a dirichlet condition in that
     * it replaces the equation with a one on the diagonal and sets
     * the residual to zero, but does not replace the solution vector
     * because the solution vector is considered 'fixed'
     */
    eqn = boundary_condition->BC_Data_Int[0];
    if (eqn <= V_LAST && eqn >= V_FIRST) {
      if (eqn != R_MASS) {
	offset = get_nodal_unknown_offset(nv, eqn, matID, 0, &vd);
	if (offset >= 0) {
	  if (node->DBC[pg->imtrx] == NULL) {
	    node->DBC[pg->imtrx] = alloc_short_1(nv->Num_Unknowns, -1);
	  }
	  node->DBC[pg->imtrx][offset] = (short int) ibc;
	}
      } else {
	/*
	 *  MASS_FRACTION SPECIAL SECTION
	 * Obtain the species number from the first integer
	 */
	w = boundary_condition->BC_Data_Int[1];
	if (w >= (upd->Max_Num_Species)) {
	  EH(GOMA_ERROR, "Species number on BC Y card exceeds number of species available");
	}
	offset = get_nodal_unknown_offset(nv, eqn, matID, w, &vd);
	if (offset < 0) {
	  WH(-1, "Y BC applied to material without species equation");
	} else {
	  if (node->DBC[pg->imtrx] == NULL) {
	    node->DBC[pg->imtrx] = alloc_short_1(nv->Num_Unknowns, -1);
	  }
	  node->DBC[pg->imtrx][offset] = (short int) ibc;
	}
      }
    }
  }
  else if(boundary_condition->BC_Name == PERIODIC_BC)
    {
      /* Nothing to do here.  If, one day, we want a more
       * productionized version of periodic BCs, then precomputing
       * matchings nodes, etc., could be done here. */
    }
  else if (boundary_condition->BC_Name == SH_SLOPE_X_BC ) {
    /*
     *  Here we need to be able to cross reference the node flagged
     *  to have a surface tangent condition with the SURFTANG card to which
     *  it corresponds.  Since the BC_Data_Int field of the 
     *  struct boundary_condition is not used for this card, we can place
     *  the global node number in it.
     */
    boundary_condition->BC_Data_Int[0] = inode;
    node->DBSH_SLOPE_X = 1;
  }
  else if (boundary_condition->BC_Name == SH_SLOPE_Y_BC ) {
    /*
     *  Here we need to be able to cross reference the node flagged
     *  to have a surface tangent condition with the SURFTANG card to which
     *  it corresponds.  Since the BC_Data_Int field of the 
     *  struct boundary_condition is not used for this card, we can place
     *  the global node number in it.
     */
    boundary_condition->BC_Data_Int[0] = inode;
    node->DBSH_SLOPE_Y = 1;
  }
  else if ( boundary_condition->BC_Name == SHEET_ENDSLOPE_BC ) {
    boundary_condition->BC_Data_Int[0] = inode;
    node->DBSES = 1;
  }    
  return (0);
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void 
alloc_First_Elem_BC (struct elem_side_bc_struct ****First_Elem_Side_BC_Array,
		     struct elem_edge_bc_struct ****First_Elem_Edge_BC_Array,
		     const int num_internal_elems)
     /*************************************************************************
      *
      * alloc_First_Elem_BC():
      *
      *    Allocates space for pointers to structures to hold side boundary
      *    conditions for all of the elements on the current processor:
      *    Note: the elements of these arrays are initialized to the NULL
      *          pointer
      *
      *   Input
      * --------
      *  num_internal_elems = Number of elements
      *
      *   Output
      * --------
      *  First_Elem_Side_BC_Array address for the pointer to the array
      *                           of length num_internal_elems containing
      *                           the pointers to structures
      *  First_Elem_Edge_BC_Array address for the pointer to the array
      *                           of length num_internal_elems containing
      *                           the pointers to structures
      ***********************************************************************/
{
  int imtrx;

  size_t sz_side = sizeof (struct elem_side_bc_struct ***);
  size_t sz_edge = sizeof (struct elem_edge_bc_struct ***);
  *First_Elem_Side_BC_Array = malloc(sz_side * (size_t) upd->Total_Num_Matrices);
  *First_Elem_Edge_BC_Array = malloc(sz_edge * (size_t) upd->Total_Num_Matrices);

  for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
    (*First_Elem_Side_BC_Array)[imtrx] = (struct elem_side_bc_struct **)
      alloc_ptr_1(num_internal_elems);
    (*First_Elem_Edge_BC_Array)[imtrx] = (struct elem_edge_bc_struct **)
      alloc_ptr_1(num_internal_elems);
  }
  return;
}
/**************************************************************************************************************************/

static int
searchEBForNode(int eb_id_target, Exo_DB *exo, int nodeTarget, int *ielemList, int listlength)
/*
 *  searchEBForNode:
 *
 *    Search an element block for a global node number. Return the number of elements containing
 *    that node and the global element numbers.
 *
 *      Input
 *     ------------
 *         eb_id_target   = ID of the element block. (The index is looked up in this routine)
 *         exo            = Exo_DB exodus database structure
 *         nodeTarget     = Global node number of the target node
 *         lstlength      = Length of the ielemList vector. No more will be found that this length
 *
 *      Output
 *     -------------
 *         return         = Returns the nunber of elements found
 *         ielemList[]    = vector of element numbers found that contains this node.
 *         
 */
{
  int bnoff, bn, bnn, ie, i, eid;
  /*
   * Find the target element block index
   */
  int ebindex_target = -1;
  for (i = 0; i < exo->num_elem_blocks; ++i)
    {
      eid = exo->eb_id[i];
      if (eid == eb_id_target)
	{
	  ebindex_target = i;
	  break;
	}
    }
  if (ebindex_target == -1) 
    {
      EH(GOMA_ERROR, "searchEBForNode Error: eb_id_target not found");
    }
  int nelems = exo->eb_num_elems[ebindex_target];
  int npe = exo->eb_num_nodes_per_elem[ebindex_target];
  int *bconn = exo->eb_conn[ebindex_target];
  int nfound = 0;
  // Loop over the elements in the element block
  for (ie = 0; ie < nelems; ie++) 
    {
      bnoff = npe * ie;
      for (bn = 0; bn < npe; ++bn)
	{
	  bnn = bconn[bnoff + bn];
	  if (bnn == nodeTarget) 
	    {
	      ielemList[nfound] = exo->eb_ptr[ebindex_target] + ie;
	      nfound++;
	      if (nfound >= listlength)
		{
		  return nfound;
		}
	    }
	}
    }
  return nfound;
}
/************************************************************************************************************************/

void
set_up_Surf_BC(struct elem_side_bc_struct **First_Elem_Side_BC_Array[ ],
	       Exo_DB *exo, Dpi *dpi)

/*****************************************************************
 * This routine takes care of several tasks up front for the application
 * of boundary conditions.  Among other things, 
 * It looks for multiple bc's being applied to the same equation at the 
 * same node, and sets up criteria for deciding which BC has precedence
 * Last revised by Rich Cairncross 12/05/95
 *****************************************************************/
{
  int 	    i, ibc, iss, ielem, num_nodes_on_side;
  int       ibc_id, ibcmp, ibc_compid, ibc_desc, ibc1, ibc2, inode = -1;
  int       poinbc, ins;
  int *local_node_list = 0;
  int nfound, elemList[5], ebid_target;
  int side;
  BOUNDARY_CONDITION_STRUCT *bc, *bc2;
  char err_msg[MAX_CHAR_IN_INPUT];
  static char *yo = "set_up_Surf_BC";
  
#ifdef DEBUG_BC
  extern int ProcID, Dim;
#endif

#ifdef DEBUG_BC
  printf("\nNumber of side sets on Proc %d = %d\n", ProcID, 
	 exo->num_side_sets);
#endif

  /*
   *  Loop over Num_BC, the number of boundary conditions defined in the 
   *  input file 
   */
  for (ibc = 0; ibc < Num_BC; ibc++) {
    /******************************************************************************/
    /*                              BLOCK 1                                       */
    /*      SURFACE INTEGRAL BOUNDARY CONDITIONS SPECIFIED BY NODE SETS           */
    /*  apply integrated conditions to nodesets with set_up_Point_BC()            */
    /*  special case: NS surf integral for 1d elements                            */
    /*  where the surface is the edge node                                        */
    /******************************************************************************/

    /* First check to see if this boundary condition is a node set 	      */
    if (!strcmp(BC_Types[ibc].Set_Type, "NS")) {

      /*
       *  Next, check to see if this boundary condition is a Surface Integral
       *  boundary condition specified by a node set. 
       *  We now have one such boundary condition. This is an edge condition on
       *  a shell equation in 2D problems. This must be specified via a node set
       *  because the boundary condition is specified at a single point. You can't
       *  have a side set as a single point.
       */
     
      // Isolate this case for now
      if (BC_Types[ibc].BC_Name == SH_GAMMA1_DERIV_SYMM_BC ||
          BC_Types[ibc].BC_Name == SH_GAMMA2_DERIV_SYMM_BC )
	{
	  /* 
	   *  Next, check to see if this boundary condition (not necessarily a
	   *  Surface Integral boundary condition) should be specified by  edges of surface shells.
	   *
	   */
	  if (BC_Types[ibc].desc->method == COLLOCATE_SURF ||
	      BC_Types[ibc].desc->method == STRONG_INT_SURF ||
	      BC_Types[ibc].desc->method == WEAK_INT_SURF ||
	      BC_Types[ibc].desc->method == WEAK_SHARP_INT ||	  
	      BC_Types[ibc].desc->method == CONTACT_SURF ||
	      BC_Types[ibc].desc->method == WEAK_SHELL_GRAD ||
	      BC_Types[ibc].desc->method == STRONG_SHELL_GRAD ) 
	    {
	      /*
	       * Resolve what material's variables this boundary condition will
	       * be applied on. If no specification of element block ID is made,
	       * set matID = -2, so that set_nodal_Dirichlet_BC() will assign the
	       * Dirichlet condition to the first variable of a variable type
	       * found at a node. If the element block ID is specified, then
	       * assign matID to the corresponding material index, so that
	       * only the variable associated with that particular material is
	       * assigned to the Dirichlet condition.  
	       */
	      /*
	      if (BC_Types[ibc].BC_EBID_Apply == -1) 
		{
		  matID = -2;
		}
	      else
		{
		  matID = map_mat_index(BC_Types[ibc].BC_EBID_Apply);
		}
	      */
	      /* Loop over the total number of node sets defined for the current processor */
	      for (ins = 0; ins < exo->num_node_sets; ins++) 
		{
		  /* 
		   *Check for a match between the ID of the current node set and the node set
		   *  ID specified in the input file - continue if a match is found 
		   */
		  if (exo->ns_id[ins] == BC_Types[ibc].BC_ID)
		    {
		      /* 
		       * Loop over all the global nodes in the current node set 
		       */
		      if ( exo->ns_num_nodes[ins] != 1)
			{
			  EH(GOMA_ERROR,"not 1");
			}
		      for (i = 0; i < exo->ns_num_nodes[ins]; i++)
			{
			  /*
			   *  Get the ith local node in the current node set
			   */
			  inode = exo->ns_node_list[ exo->ns_node_index[ins] + i ];
			  if (BC_Types[ibc].BC_Name == SH_GAMMA1_DERIV_SYMM_BC ||
			      BC_Types[ibc].BC_Name == SH_GAMMA2_DERIV_SYMM_BC) 
			    {
			      
			      ebid_target = BC_Types[ibc].BC_Data_Int[0];
                           
			      nfound = searchEBForNode(ebid_target, exo, inode, elemList, 5);
			      if (nfound != 1) 
				{
				  EH(GOMA_ERROR, "nfound != 1");
				}
			      num_nodes_on_side = 1;
			      local_node_list[0] = inode;
			      ielem = elemList[0];
#ifdef DEBUG_SS 
			      fprintf(stderr, "\tside = not determined yet");
			      fprintf(stderr, "\tielem = %d, nodes = ", elemList[0]);
			      for (j = 0; j < num_nodes_on_side; j++) {
				fprintf(stderr, "local_node_list[%d] = %d\n", j, inode);
			      }
#endif
                              if (BC_Types[ibc].matrix >= 0) {
                                setup_Elem_BC(&First_Elem_Side_BC_Array[BC_Types[ibc].matrix][ielem], &BC_Types[ibc],
                                              ibc, num_nodes_on_side, ielem, 
                                              local_node_list, exo);
                              }

			    }
			}  /* END for (i = 0; i < Proc_NS_Count[ins]; i++)           */
		    }  /* END if (Proc_NS_Ids[ins] == BC_Types[ibc].BC_ID)	     */
		}  /* END for (ins = 0; ins < Proc_Num_Node_Sets; ins++) 	     */
	    }
	}



    }  /* END if (!strcmp(BC_Types[ibc].Set_Type, "NS")) 		     */

    /*****************************************************************************/
    /*                              BLOCK 2                                      */
    /*      BOUNDARY CONDITIONS SPECIFIED BY SIDE SETS                           */
    /*****************************************************************************/

    /* First check to see if this boundary condition is a side set 	     */
    if (!strcmp(BC_Types[ibc].Set_Type, "SS")) {

      /* Next, check to see if this boundary condition (not necessarily a
	 Surface Integral boundary condition) should be specified by side-sets */
      if (BC_Types[ibc].desc->method == COLLOCATE_SURF ||
	  BC_Types[ibc].desc->method == STRONG_INT_SURF ||
	  BC_Types[ibc].desc->method == WEAK_INT_SURF ||
	  BC_Types[ibc].desc->method == WEAK_SHARP_INT ||	  
	  BC_Types[ibc].desc->method == CONTACT_SURF ||
          BC_Types[ibc].desc->method == WEAK_SHELL_GRAD ||
          BC_Types[ibc].desc->method == STRONG_SHELL_GRAD ) {

	/* Loop over the total number of side sets defined on the current processor */



	for (iss = 0; iss < exo->num_side_sets; iss++) {


          /* 
	   * Check for a match between the ID of the current side set
	   * and the side set ID specified in the input file 
	   * - continue if a match is found 
	   */
	  if (exo->ss_id[iss] == BC_Types[ibc].BC_ID) {
	    /*
	     * Loop over all the global elements in the current
	     * side set 
	     */
	    for (i = 0; i < exo->ss_num_sides[iss]; i++) {
	      ielem = exo->ss_elem_list[exo->ss_elem_index[iss]+i];
	      side  = exo->ss_side_list[exo->ss_elem_index[iss]+i];
	      /*
		for (i = 0; i < Proc_SS_Elem_Count[iss]; i++) {
	      */

	      /* Get the local element number of the ith element in the 
	       * current side set 
	       * ielem = Proc_SS_Elem_List[Proc_SS_Elem_Pointers[iss] + i];
	       */
	      /* Get the number of nodes on the side of the current element 
	       * NOTE - Currently, this must be the same for all elements 
	       * in the side set, a fairly big limitation	
	       */
	      num_nodes_on_side = ( exo->ss_node_side_index[iss][i+1] -
				    exo->ss_node_side_index[iss][i] );

	      if (num_nodes_on_side < 1) {
		fprintf(stderr, 
			"P_%d: SS[%d]=%d, index=%d, elem=(%d), side=%d has num_nodes_on_side=%d!\n",
			ProcID, iss, exo->ss_id[iss], i, ielem+1, side,
			num_nodes_on_side);
	      }
	      local_node_list = 
		&(exo->ss_node_list[iss][exo->ss_node_side_index[iss][i]]);
#ifdef DEBUG_SS 
	      fprintf(stderr, "\tside = %d", side);
	      fprintf(stderr, "\tielem = %d, nodes = ", ielem);
	      for (j = 0; j < num_nodes_on_side; j++) {
		fprintf(stderr, "local_node_list[%d] = %d\n",
			j, local_node_list[j]);
	      }
#endif
              if (BC_Types[ibc].matrix >= 0) {
                setup_Elem_BC(&First_Elem_Side_BC_Array[BC_Types[ibc].matrix][ielem], &BC_Types[ibc],
                              ibc, num_nodes_on_side, ielem, 
                              local_node_list, exo);
              }

	    }  /* END for (i = 0; i < Proc_SS_Elem_Count[iss]; i++)          */
	  }  /* END if (Proc_SS_Ids[iss] == BC_Types[ibc].BC_ID)	     */
	}  /* END for (iss = 0; iss < Proc_Num_Side_Sets; iss++) 	     */
      }  /* END if (BC_Types[ibc].BC_Name == TNRMLSIDE_BC || ....	     */
    }  /* END if (!strcmp(BC_Types[ibc].Set_Type, "SS")) 		     */
  }  /* END for (ibc = 0; ibc < Num_BC; ibc++)				     */

  /*****************************************************************************/
  /*                              BLOCK 4                                      */
  /*    CHECK FOR LEAKY BOUNDARY CONDITIONS AND CONNECTION TO FLUX CONDITIONS  */
  /*****************************************************************************/

  for (ibc = 0; ibc < Num_BC; ibc++) {
    bc = BC_Types + ibc;
    if (bc->BC_Name == SDC_STEFANFLOW_BC) {
      ibc_id = bc->BC_ID;
      bc->BC_Data_Int[1] = 0;
      ibc1 = 2;
      for (ibc2 = 0; ibc2 < Num_BC; ibc2++) {   
	bc2 = BC_Types + ibc2;
	if ((bc2->BC_ID == ibc_id) &&
	    (bc2->BC_Name == VL_EQUIL_PRXN_BC)) {
	  bc->BC_Data_Int[1]++;
	  bc->BC_Data_Int[ibc1] = ibc2;
	  ibc1++;
	  bc2->BC_Data_Int[3] = ibc;
	}
      }
    }
  }
  
  /* 
   * Loop over Num_BC, the number of boundary conditions defined in 
   * the input file to find any kinematic or velo-normal conditions 
   * and any heat flux conditions that might be linked to mass flux 
   * conditions by latent heat of vaporization etc.   
   */
  DPRINTF(stdout, "\n");
  for (ibc = 0; ibc < Num_BC; ibc++) 
    {
      if (BC_Types[ibc].BC_Name == KIN_LEAK_BC ||
          BC_Types[ibc].BC_Name == KIN_CHEM_BC ||
          BC_Types[ibc].BC_Name == KIN_ELECTRODEPOSITION_BC ||  /*  RSL 5/27/02  */
          BC_Types[ibc].BC_Name == VNORM_ELECTRODEPOSITION_BC ||  /*  RSL 5/30/02  */
          BC_Types[ibc].BC_Name == QRAD_BC ||
          BC_Types[ibc].BC_Name == QCONV_BC ||
          BC_Types[ibc].BC_Name == VNORM_LEAK_BC ||
          BC_Types[ibc].BC_Name == LATENT_HEAT_BC ||
	  BC_Types[ibc].BC_Name == LS_EIK_KIN_LEAK_BC ||
	  BC_Types[ibc].BC_Name == LS_EXTV_KIN_LEAK_BC ||
	  BC_Types[ibc].BC_Name == LS_EXTV_LATENT_BC )

	{
	  /*
	   * IF a LEAK condition is found, loop to find flux 
	   * Conditions that describe the rate of the leak
	   */

          ibc_id=BC_Types[ibc].BC_ID;
	  ibc_compid = BC_Types[ibc].species_eq;
          ibc_desc=BC_Types[ibc].BC_Desc_index;
	  ibc1=-1 ;
	  ibcmp = -1 ;
	  DPRINTF(stdout, "Leaky %s BC found on %2s %d\n",
		  BC_Desc[ibc_desc].name1, 
		  BC_Types[ibc].Set_Type, ibc_id);

          for (ibc2 = 0; ibc2 < Num_BC; ibc2++) 
	    {
	      /*
	       * Assume grouping for either of the 3 cases.
	       */
	      if ( ( BC_Types[ibc2].BC_ID == ibc_id ) && 
		   ( ( BC_Types[ibc2].BC_Name == YFLUX_BC ) ||
		     ( BC_Types[ibc2].BC_Name == LS_YFLUX_BC ) ||
		     ( BC_Types[ibc2].BC_Name == LS_LATENT_HEAT_BC ) ||
                     ( BC_Types[ibc2].BC_Name == YFLUX_BV_BC ) ||
                     ( BC_Types[ibc2].BC_Name == YFLUX_HOR_BC ) ||
                     ( BC_Types[ibc2].BC_Name == YFLUX_ORR_BC ) ||
                     ( BC_Types[ibc2].BC_Name == YFLUX_H2O_ANODE_BC ) ||
                     ( BC_Types[ibc2].BC_Name == YFLUX_H2O_CATHODE_BC ) ||
                     ( BC_Types[ibc2].BC_Name == YFLUX_SULFIDATION_BC ) ||
                     ( BC_Types[ibc2].BC_Name == YFLUX_BV2_BC ) ||  /*  RSL 8/8/01  */
                     ( BC_Types[ibc2].BC_Name == YFLUX_NI_BC ) ||  /*  RSL 8/8/01  */
		     ( BC_Types[ibc2].BC_Name == YFLUX_EQUIL_BC ) ||
		     ( BC_Types[ibc2].BC_Name == YFLUX_USER_BC ) ||
		     ( BC_Types[ibc2].BC_Name == YFLUX_CONST_BC ) ||
		     ( BC_Types[ibc2].BC_Name == POROUS_LIQ_FLUX_CONST_BC ) ||
		     ( BC_Types[ibc2].BC_Name == POROUS_GAS_FLUX_CONST_BC ) ||
		     ( BC_Types[ibc2].BC_Name == POROUS_FLUX_BC ) ) )
		{
		  BC_Types[ibc2].BC_Data_Int[1] = ibc1;
		  ibc1 = ibc2;
		  if ((ibc_compid != -1) && (ibc_compid == BC_Types[ibc2].species_eq)) 
		    ibcmp = ibc2;
		  DPRINTF(stdout, "Leaky BC connection found to flux %d %d\n",
			  ibc, ibc2);
		}
	    }
	  /*
	   * Set integer array for LEAK condition to BC number
	   * of last Flux condition
	   */
	  if (ibcmp != -1)
	    {
	      BC_Types[ibc].BC_Data_Int[1]=ibcmp;
	    }
	  else
	    {
	      BC_Types[ibc].BC_Data_Int[1]=ibc1;
	    }
	}
    }

  /*****************************************************************************/
  /*                              BLOCK 5                                      */
  /*CHECK FOR DYNAMIC CONTACT LINE POINBC NODESET FOR SLIP WITH VELO_TANGENT_BC*/
  /*****************************************************************************/
  /* Loop over Num_BC, the number of boundary conditions defined in the 
     input file to find any velo_tangent conditions that have a nonzero
     poinbc associated with them to give the location of the dynamic contact
     line. Replace the poinbc number in BC_Types[ibc].BC_Data_Int[0] with the
     global node number                                                        */

  for (ibc = 0; ibc < Num_BC; ibc++) {
    if (BC_Types[ibc].BC_Name == VELO_TANGENT_BC ||
 	BC_Types[ibc].BC_Name == VELO_TANGENT_USER_BC ||
	BC_Types[ibc].BC_Name == VELO_SLIP_BC ||
	BC_Types[ibc].BC_Name == VELO_SLIP_ROT_BC ||
	BC_Types[ibc].BC_Name == VELO_SLIP_FILL_BC ||
	BC_Types[ibc].BC_Name == VELO_SLIP_ROT_FILL_BC ||
	BC_Types[ibc].BC_Name == AIR_FILM_BC ||
	BC_Types[ibc].BC_Name == AIR_FILM_ROT_BC ||
	BC_Types[ibc].BC_Name == VELO_SLIP_FLUID_BC ||
	BC_Types[ibc].BC_Name == VELO_SLIP_ROT_FLUID_BC ||
	BC_Types[ibc].BC_Name == VELO_STREAMING_BC ) {
      poinbc = BC_Types[ibc].BC_Data_Int[0];
      if (poinbc != 0 && poinbc != -1) {    /*the -1 case is another new
					      velo_tangent retaining case, 
					      prs 8/98 */
        /* Set to Flag for node not found */
        inode = -2;
        /* Change -1 & DCL nset_id case to positive id  */
        if(poinbc < 0) poinbc = -poinbc;

	for (ins = 0; ins < exo->num_node_sets; ins++) {
	  if (exo->ns_id[ins] == poinbc) {
	    for (i = 0; i < exo->ns_num_nodes[ins]; i++) {
	      inode = exo->ns_node_list[exo->ns_node_index[ins] + i];
	    }
	  }	
	}
	BC_Types[ibc].BC_Data_Int[0] = inode;  
	if (Debug_Flag > 1) {
	  printf("set_up_Surf_BC: point BC= %d corresponds to node= %d\n",
	         poinbc, inode) ;
	}
	  
      }
    }
    if (BC_Types[ibc].BC_Name == VELO_TANGENT_SOLID_BC ||
	BC_Types[ibc].BC_Name == VELO_SLIP_SOLID_BC ) {
      poinbc = BC_Types[ibc].BC_Data_Int[2];
      if (poinbc != 0 && poinbc != -1) { 
	for (ins = 0; ins < exo->num_node_sets; ins++) {
	  if (exo->ns_id[ins] == poinbc) {
	    for (i = 0; i < exo->ns_num_nodes[ins]; i++) {
	      inode = exo->ns_node_list[exo->ns_node_index[ins] + i];
	    }
	  }	
	}
	BC_Types[ibc].BC_Data_Int[2] = inode;  
	if (Debug_Flag > 1) {
	  printf("set_up_Surf_BC: point BC= %d corresponds to node= %d\n",
	         poinbc, inode) ;
	}
	  
      }
    }
  }
  for (ibc = 0; ibc < Num_BC; ibc++) {
    if (BC_Types[ibc].BC_Name == ROLL_FLUID_BC)
    {
      for (ibc2 = 0; ibc2 < Num_BC; ibc2++) {
         if (BC_Types[ibc2].BC_Name == VELO_SLIP_ROT_FLUID_BC)
           {
            BC_Types[ibc].BC_Data_Int[2] = ibc2;
           }
        }
    }
  }

  /*****************************************************************************/
  /*                              BLOCK 6                                      */
  /*   CHECK FOR Q_VELO_SLIP CONDITIONS AND CORRESPONDING VELO_SLIP CONDITION  */
  /*****************************************************************************/
  /* Loop over Num_BC, the number of boundary conditions defined in 
     the input file to find any Q_VELO_SLIP conditions and identify
     the corresponding VELO_SLIP or VELO_SLIP_ROT conditions */

  DPRINTF(stdout, "\n");
  for (ibc = 0; ibc < Num_BC; ibc++) 
    {
      if ( BC_Types[ibc].BC_Name == Q_VELO_SLIP_BC )

	{
	  /*
	   * IF a Q_VELO_SLIP condition is found, loop to find 
	   * corresponding VELO_SLIP or VELO_SLIP_ROT condition
	   */

          ibc_id=BC_Types[ibc].BC_ID;
          ibc_desc=BC_Types[ibc].BC_Desc_index;
          DPRINTF(stdout, "Q_VELO_SLIP BC found on %2s %d\n", BC_Types[ibc].Set_Type, ibc_id);

          /* now search for corresponding slip condition */
          BC_Types[ibc].BC_Data_Int[0] = -1;
          for (ibc2 = 0; ibc2 < Num_BC; ibc2++) 
	    {
	      if ( ( BC_Types[ibc2].BC_ID == ibc_id ) && 
		   ( ( BC_Types[ibc2].BC_Name == VELO_SLIP_BC ) ||
		     ( BC_Types[ibc2].BC_Name == VELO_SLIP_ROT_BC ) ||
		     ( BC_Types[ibc2].BC_Name == VELO_SLIP_FLUID_BC ) ||
		     ( BC_Types[ibc2].BC_Name == VELO_SLIP_ROT_FLUID_BC ) ||
		     ( BC_Types[ibc2].BC_Name == AIR_FILM_BC ) ||
		     ( BC_Types[ibc2].BC_Name == AIR_FILM_ROT_BC ) ) ) 
		{
		  /*
		   * Set integer to BC number of matching SLIP CONDITON
		   */
                  BC_Types[ibc].BC_Data_Int[0] = ibc2;
		  DPRINTF(stdout, "%s BC found for Q_VELO_SLIP BC on %2s %d\n", 
			  BC_Types[ibc2].desc->name2, BC_Types[ibc].Set_Type, ibc_id);
                  
		}
	    }
            
	  if (BC_Types[ibc].BC_Data_Int[0] == -1)
	    {
	      sr = sprintf(err_msg, 
			   "%s: Q_VELO_SLIP specified without corresponding VELO_SLIP condition", yo);
	      EH(GOMA_ERROR, err_msg);
	    }
	}
    }


  /* Initialize rotation variables and lists */
  mesh_rotate_node = calloc((size_t) upd->Total_Num_Matrices, sizeof(int *));
  mesh_rotate_ss = calloc((size_t) upd->Total_Num_Matrices, sizeof(int *));
  num_mesh_rotate = calloc((size_t) upd->Total_Num_Matrices, sizeof(int));
  mom_rotate_node = calloc((size_t) upd->Total_Num_Matrices, sizeof(int *));
  mom_rotate_ss = calloc((size_t) upd->Total_Num_Matrices, sizeof(int *));
  num_mom_rotate = calloc((size_t) upd->Total_Num_Matrices, sizeof(int));

  if (Num_ROT == 0) check_for_bc_conflicts2D(exo, dpi);
//  if (Num_ROT == 0 && exo->num_dim == 3) setup_rotated_bc_nodes(exo, BC_Types, Num_BC);
  if (Num_ROT > 0)  check_for_bc_conflicts3D(exo, dpi);

  return;
} /* End of set_up_Surf_BC */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void
free_Surf_BC(struct elem_side_bc_struct **First_Elem_Side_BC_Array[], Exo_DB *exo)
{
  int e_start, e_end;
  int ielem;
  int imtrx;
  struct elem_side_bc_struct *current_bc;
  struct elem_side_bc_struct *next_bc;
	
  e_start = exo->eb_ptr[0];
  e_end   = exo->eb_ptr[exo->num_elem_blocks];
  for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
    for( ielem= e_start; ielem < e_end; ielem++)
      {
        if ( First_Elem_Side_BC_Array[imtrx][ ielem ] != NULL )
          {
            current_bc = First_Elem_Side_BC_Array[imtrx][ ielem ];
			
            do {
              safer_free( (void **) &current_bc->local_node_id );
              safer_free( (void **) &current_bc->local_elem_node_id );
              /*				elem_qp_storage_free( current_bc );*/
				
              next_bc = current_bc->next_side_bc;
				
              safer_free( (void **) &current_bc );
				
              current_bc = next_bc;
            }
            while ( current_bc != NULL );
			
          }
      }
  }
	
  safer_free ( (void **) &First_Elem_Side_BC_Array );

}
			
			
void
free_Edge_BC ( struct elem_edge_bc_struct **First_Elem_Edge_BC_Array[ ],
	       Exo_DB *exo, Dpi *dpi)
{
  int e_start, e_end;
  int ielem;
  int imtrx;
  struct elem_edge_bc_struct *current_bc;
  struct elem_edge_bc_struct *next_bc;
	
  e_start = exo->eb_ptr[0];
  e_end   = exo->eb_ptr[exo->num_elem_blocks];
  for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
    for( ielem= e_start; ielem < e_end; ielem++)
      {
        if ( First_Elem_Edge_BC_Array[imtrx][ ielem ] != NULL )
          {
            current_bc =  First_Elem_Edge_BC_Array[imtrx][ ielem ];
			
            do {
              safer_free( (void **) &current_bc->elem_side_bc_1->local_node_id );
              safer_free( (void **) &current_bc->elem_side_bc_1->local_elem_node_id );
              safer_free( (void **) &current_bc->elem_side_bc_2->local_node_id );
              safer_free( (void **) &current_bc->elem_side_bc_2->local_elem_node_id );
				
              safer_free ( (void **) &current_bc->elem_side_bc_1 );
              safer_free ( (void **) &current_bc->elem_side_bc_2 );
              next_bc = current_bc->next_edge_bc;
				
              safer_free( (void **) &current_bc );
				
              current_bc = next_bc;
            }
            while ( current_bc != NULL );
          }
      }
  }	
  safer_free ( (void **) &First_Elem_Edge_BC_Array );

}
	

void
setup_Point_BC (struct elem_side_bc_struct **First_Elem_Side_BC_Array[ ],
Exo_DB *exo, Dpi *dpi) {
  char err_msg[MAX_CHAR_IN_INPUT];
  //int i;
  int ibc, ielem, ns_id, ins;
  int is_ns_g, found_ns_local, found_ns_global;
  int num_nodes_on_side = 1; // This is the edge of a 1D element
  //int ipin = 0;
  int node_list[MDE]; /* Even though exodus is spec'd for up to 3 nodes
  *                      on a 1d element.
  * * * * * * * * * * * * * * * * * * */
  //struct elem_side_bc_struct *this_side_bc = NULL;

  for (ibc = 0; ibc < Num_BC; ibc++) {
    /* First check to see if this boundary condition is a nodeset         */
    if (!strcmp(BC_Types[ibc].Set_Type, "NS")) {
      /* Make sure the condition is Neumann type (necessary?) */
      //struct Boundary_Condition thisBC;
      //thisBC = BC_Types[ibc];
      if (BC_Types[ibc].desc->method == WEAK_INT_SURF ||
          BC_Types[ibc].desc->method == STRONG_INT_SURF) {

        // This is a Neumann condition that gets applied to a node set
        // What if the nodeset has more than one node?
        // Get the NS_ID for this condition
        ns_id = BC_Types[ibc].BC_ID;
        // see if I (in the mpi processors sense) have this nodeset
        ins = in_list(ns_id, 0, exo->num_node_sets,exo->ns_id);
        // record truth value of having the nodeset in the mesh
        found_ns_local = (ins > -1);

	      /*
	       * Take advantage of some limited global information that has
	       * been cached in dpi to check if some sideset complies. This
	       * beats doing the MPI_Allreduce(...Logical OR) with the associated
	       * communications overhead.
	       */

        #ifdef PARALLEL
          // check for nodeset in global nodeset list
	        is_ns_g  = in_list(ns_id, 0, dpi->num_node_sets_global,
          dpi->ns_id_global);

          if ( is_ns_g > -1) {
			      found_ns_global = TRUE;
          } else {
            found_ns_global = FALSE;
          }

        #endif
        #ifndef PARALLEL
          // if not parallel, local is global
	        found_ns_global = found_ns_local;
        #endif

        if ( !found_ns_global ) {
          sr = sprintf(err_msg, "NS_ID %d for BC (%d) not found",
          ns_id, ibc+1);
          EH(GOMA_ERROR, err_msg);
        }

        // Here we need to use setup_Elem_Edge_BC to attach this
        // BC to all its elements on this processor.
        // Shouldn't there be only 1 element attached to an 
        // edge-node BC?
        if ( found_ns_local ) {
          node_list[0] = exo->ns_node_list[ins];
          ielem = exo->node_elem_list[
                  exo->node_elem_pntr[
                  exo->ns_node_list[ins]]];

          setup_Elem_BC (&First_Elem_Side_BC_Array[BC_Types[ibc].matrix][ielem],
				    &BC_Types[ibc],
				    ibc, num_nodes_on_side,
				    ielem,
				    node_list, exo);
        }
      }
    }
  }
  return;
}

void
set_up_Edge_BC (struct elem_edge_bc_struct **First_Elem_Edge_BC_Array[ ],
		Exo_DB *exo, Dpi *dpi)
{
  int num_nodes_on_side2;	/* of the 2nd side set, when checking for
				 * edge intersections of two SS's */
  int j;
  int ss_id1;
  int ss_id2;
  int kdex;
  int k2_start;
  char err_msg[MAX_CHAR_IN_INPUT];
  /*****************************************************************
   * this routine takes care of several tasks up front for the application
   * of boundary conditions.  Among other things, 
   * It looks for multiple bc's being applied to the same equation at the 
   * same node, and sets up criteria for deciding which BC has precedence
   * Last revised by Rich Cairncross 12/05/95
   *****************************************************************/

  /* Local variables */

  int i, ibc, iss1, iss2, ielem, num_nodes_on_side, num_nodes_on_edge, ipin;
  int found_ss_primary_local;
  int found_ss_primary_global;
  int found_ss_secondary_local;
  int found_ss_secondary_global;
  int l, k, k_node;
  int node_list[MDE];
  int node_ctr = -1;
#ifdef PARALLEL
  int iss1g, iss2g;
#endif
  struct elem_edge_bc_struct *this_edge_bc = NULL;

  /***************************** EXECUTION BEGINS ******************************/

  /*****************************************************************************/
  /*                              BLOCK 1                                      */
  /*      BOUNDARY CONDITIONS SPECIFIED AS THE INTERSECTION OF TWO SS          */
  /*****************************************************************************/
  /* Loop over Num_BC, the number of boundary conditions defined in the 
     input file */

  for (ibc = 0; ibc < Num_BC; ibc++) {

    /* First check to see if this boundary condition is a side set 	      */
    if (!strcmp(BC_Types[ibc].Set_Type, "SS")) {
      
      /* Next, check to see if this boundary condition (not necessarily a
	 Surface Integral boundary condition) should be specified by side-sets */
      if (BC_Types[ibc].desc->method == COLLOCATE_EDGE ||
	  BC_Types[ibc].desc->method == WEAK_INT_EDGE ||
	  BC_Types[ibc].desc->method == STRONG_INT_EDGE ) 
	{

	  num_nodes_on_edge = -1;
	  ipin = 1;  /*Assume these nodes are pinned on startup if CA_EDGE_OR_FIX
		       BC is used.  */

	  /*
	   * This cross check needs to be generalized for parallel usage.
	   * That is, it must be OK for this processor to NOT have the
	   * requisite primary and secondary sidesets that make up each and
	   * every BC. All that is required is that *some* processor be
	   * able to successfully fulfill these requirements.
	   */


	  /* Load list of elements along this EDGE */
	  /* get list along primary SS */

	  ss_id1 = BC_Types[ibc].BC_ID;
	  iss1   = in_list(ss_id1, 0, exo->num_side_sets, exo->ss_id);

	  found_ss_primary_local  = ( iss1 > -1 );

	  /*
	   * Take advantage of some limited global information that has
	   * been cached in dpi to check if some sideset complies. This
	   * beats doing the MPI_Allreduce(...Logical OR) with the associated
	   * communications overhead.
	   */

#ifdef PARALLEL
	  iss1g  = in_list(ss_id1, 0, dpi->num_side_sets_global, 
			   dpi->ss_id_global);
	  found_ss_primary_global = ( iss1g > -1);
#endif
#ifndef PARALLEL
	  found_ss_primary_global = found_ss_primary_local;
#endif
	  if ( ! found_ss_primary_global ) /* OK, there really is a problem! */
	    {
	      sr = sprintf(err_msg, "Primary SS_ID %d for BC (%d) not found",
			   ss_id1, ibc+1);
	      EH(GOMA_ERROR, err_msg);
	    }


	  /* get list along secondary SS */
	  ss_id2 = BC_Types[ibc].BC_ID2;
	  iss2   = in_list(ss_id2, 0, exo->num_side_sets, exo->ss_id);

	  found_ss_secondary_local = ( iss2 > -1 );

#ifdef PARALLEL
	  iss2g  = in_list(ss_id2, 0, dpi->num_side_sets_global, 
			   dpi->ss_id_global);
	  found_ss_secondary_global = ( iss2g > -1);
#endif
#ifndef PARALLEL
	  found_ss_secondary_global = found_ss_secondary_local;
#endif

	  if ( ! found_ss_secondary_global ) /* OK, there really is a problem! */
	    {
	      sr = sprintf(err_msg, "Secondary SS_ID %d for BC (%d) not found",
			   ss_id2, ibc+1);
	      EH(GOMA_ERROR, err_msg);
	    }


	  /* 
	   * If any of these elements exist on this processor, then....
	   *
	   *
	   * Loop over all the global elements in the primary side set 
	   *
	   * NOTE: currently assuming both side-sets are in the same element 
	   * (i.e. no neclacing) 
	   */

	  if ( found_ss_primary_local && found_ss_secondary_local )
	    {
	      for ( i=0; i<exo->ss_num_sides[iss1]; i++)
		{
		  ielem = exo->ss_elem_list[exo->ss_elem_index[iss1]+i];
		  num_nodes_on_side = ( exo->ss_node_side_index[iss1][i+1] -
					exo->ss_node_side_index[iss1][i] );

		  /*
		   * find out if this element exists on the secondary SS and
		   * if so, then where
		   */

		  j = in_list(ielem, 0, exo->ss_num_sides[iss2], 
			      &(exo->ss_elem_list[exo->ss_elem_index[iss2]]));
	   
		  /* j is the element number of element in SS2 */
	  
		  if ( j != -1 )
		    {
		      /*
		       * Found that ielem is part of SS2, too!
		       *
		       * Now count up the nodes that are common to both 
		       * SS1 side i and SS2 side j 
		       * and build a list to contain their names.
		       */

		      node_ctr = 0;

		    
		      k2_start = exo->ss_node_side_index[iss2][j];

		      num_nodes_on_side2 = 
			( exo->ss_node_side_index[iss2][j+1] -
			  exo->ss_node_side_index[iss2][j] );

		      for (k=0; k<num_nodes_on_side; k++)
			{

			  kdex   = exo->ss_node_side_index[iss1][i] + k;
			  k_node = exo->ss_node_list[iss1][kdex];

			  l = in_list(k_node, 0, num_nodes_on_side2, 
				      &(exo->ss_node_list[iss2][k2_start]));
	      

			  if ( l != -1 )
			    {
			      node_list[node_ctr] = k_node;
			      node_ctr++;
			    }
			}

		      if (node_ctr == 0) EH(GOMA_ERROR, "No nodes on edge");

		      if (num_nodes_on_edge != -1) 
			{ 
			  /* check to make sure the number of edge nodes
			   * doesn't change 
			   */
			  if (num_nodes_on_edge != node_ctr) 
			    {
			      EH(GOMA_ERROR, "Number of edge nodes varies w/ element!");
			    }
			  num_nodes_on_edge = node_ctr;
			} 
		      else 
			{
			  num_nodes_on_edge = node_ctr;
			}
		    
                      if (BC_Types[ibc].matrix >= 0) {
                        this_edge_bc = setup_Elem_Edge_BC (&First_Elem_Edge_BC_Array[BC_Types[ibc].matrix][ielem],
                                                           &BC_Types[ibc],
                                                           ibc, num_nodes_on_edge, ipin, FALSE,
                                                           ielem, 
                                                           node_list, exo);

                        setup_Elem_BC (&(this_edge_bc->elem_side_bc_1), 
                                       &BC_Types[ibc],
                                       ibc, num_nodes_on_side, ielem, 
                                       &(exo->ss_node_list[iss1][exo->ss_node_side_index[iss1][i]]), exo);


                        setup_Elem_BC (&(this_edge_bc->elem_side_bc_2), 
                                       &BC_Types[ibc],
                                       ibc, num_nodes_on_side, ielem, 
                                       &(exo->ss_node_list[iss2][exo->ss_node_side_index[iss2][j]]), exo);	
                      }
   
		    }
		  else
		    {
		      /*
		       * Handle the case when the edge in question is shared by two elements
		       * That is, the "necklaced" cases
		       */

		      /*  First, look through all neighbor elements for the one that 
		       *  belongs to ss2.  If more than one neighbor belongs give up for 
		       *  now
		       */
		    
		      int ielem2 = -1;
		      int found;

		      for( l = exo->elem_elem_pntr[ielem], found = FALSE; 
			   found == FALSE && l < exo->elem_elem_pntr[ielem+1]; 
			   l++)
			{

			  j = in_list(exo->elem_elem_list[l], 0, exo->ss_num_sides[iss2], 
				      &(exo->ss_elem_list[exo->ss_elem_index[iss2]]));

			  /* j is the side number of the neighbor whose side is on ss2 */

			  if ( j != -1 )
			    {
			      /* 
			       * One of neighbor elements has a side on sideset two
			       *
			       * Does this side share any nodes with the side on the primary side ?
			       */

			      ielem2 =  exo->elem_elem_list[l];

			      node_ctr = 0;
			    
			      k2_start = exo->ss_node_side_index[iss2][j];

			      num_nodes_on_side2 =  ( exo->ss_node_side_index[iss2][j+1] -
						      exo->ss_node_side_index[iss2][j] );
			    
			      /*
			       * Next we search for nodes on sideset 1 which also appear in sideset2 at side 
			       * indices i and j respectively
			       * 
			       * This set is the global nodes on the edge.  
			       *
			       */

			      for ( k = exo->ss_node_side_index[iss1][i]; k<exo->ss_node_side_index[iss1][i+1];  k++)
				{
				
				  k_node = exo->ss_node_list[iss1][ k ];

				  if( in_list(k_node, 0, num_nodes_on_side2, 
					      &(exo->ss_node_list[iss2][k2_start])) != -1 )
				    {
				      node_list[node_ctr++] = k_node;
				    }
				}

			    
			      found = node_ctr > 1 ? TRUE : FALSE;

 
			    } /* end of j != -1 */
			} /* end of for(l = ... */

		      if ( found == TRUE )
			{


			  if (num_nodes_on_edge != -1) 
			    { 
			      /* check to make sure the number of edge nodes
			       * doesn't change 
			       */
			      if (num_nodes_on_edge != node_ctr) 
				{
				  EH(GOMA_ERROR, "Number of edge nodes varies w/ element!");
				}
			      num_nodes_on_edge = node_ctr;
			    } 
			  else 
			    {
			      num_nodes_on_edge = node_ctr;
			    }

			  /*
			   * Set up elem_edge_bc structure for this edge */

                          if (BC_Types[ibc].matrix >= 0) {			
                            this_edge_bc = setup_Elem_Edge_BC (&First_Elem_Edge_BC_Array[BC_Types[ibc].matrix][ielem],
                                                               &BC_Types[ibc],
                                                               ibc, num_nodes_on_edge, ipin, TRUE,
                                                               ielem, 
                                                               node_list, exo);

                            /* create elem_side_bc's for SS1 and SS2 */


                            setup_Elem_BC (&(this_edge_bc->elem_side_bc_1), 
                                           &BC_Types[ibc],
                                           ibc, num_nodes_on_side, ielem, 
                                           &(exo->ss_node_list[iss1][exo->ss_node_side_index[iss1][i]]), exo);


                            setup_Elem_BC (&(this_edge_bc->elem_side_bc_2), 
                                           &BC_Types[ibc],
                                           ibc, num_nodes_on_side, ielem2, 
                                           &(exo->ss_node_list[iss2][exo->ss_node_side_index[iss2][j]]), exo);
                          }
			} /* end of if (found ) */	


		    } /* end of if ( j!= -1)	*/

		    
		}  /* END for (i = 0; i < Proc_SS_Elem_Count[iss]; i++)        */

	    } /* END if ( found_ss_primary_local && found_ss_secondary_local ) */

	}  /* END if (BC_Types[ibc].BC_Name == TNRMLSIDE_BC || ....      */

    } /* END if (!strcmp(BC_Types[ibc].Set_Type, "SS")) 	      */

  }  /* END for (ibc = 0; ibc < Num_BC; ibc++)			      */

  return;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void
set_up_Embedded_BC (void)

     /*****************************************************************
      * this routine creates a list of bc's applied on level
      * set surfaces
      *****************************************************************/
{
  int 	                ibc;
  struct LS_Embedded_BC *bc;
  /* int                   found_capillary=FALSE; */
  int                   pf;

  /*****************************************************************************/
  /*      "BOUNDARY" CONDITIONS SPECIFIED BY LEVEL SET                         */
  /*****************************************************************************/

  /*
   *  Loop over Num_BC, the number of boundary conditions defined in the 
   *  input file and make list of those applied on level set surfaces
   */
  
  if (ls == NULL && pfd == NULL) return;
  
  if ( ls != NULL ) ls->embedded_bc = NULL;
  if ( pfd != NULL )
    {
      for(pf=0; pf<pfd->num_phase_funcs; pf++)
        {
          pfd->ls[pf]->embedded_bc = NULL;
        }
    }
  
  for (ibc = 0; ibc < Num_BC; ibc++) {
  
    if (!strcmp(BC_Types[ibc].Set_Type, "LS")) {
      
      bc = (struct LS_Embedded_BC *) smalloc( sizeof( struct LS_Embedded_BC ) );
      
      if ( ls == NULL)
        {
          EH(GOMA_ERROR,"Attempt to apply bc on LS but level set is not active, did you intend PF?\n");
        }
        
      bc->next = ls->embedded_bc;
      ls->embedded_bc = bc;
      
      bc->bc_input_id = ibc;
      
      /* if ( BC_Types[ibc].BC_Name == LS_CAPILLARY_BC ) found_capillary = TRUE; */

    } else if (!strcmp(BC_Types[ibc].Set_Type, "PF")) {
      
      bc = (struct LS_Embedded_BC *) smalloc( sizeof( struct LS_Embedded_BC ) );
      pf = abs(BC_Types[ibc].BC_ID) - 1;
      
      if ( pfd == NULL )
        {
          EH(GOMA_ERROR,"Attempt to apply bc on PF but phase functions are not active.\n");
        }
        
      if ( pf >= pfd->num_phase_funcs )
        {
          EH(GOMA_ERROR,"Attempt to apply bc to invalid phase function.\n");
        }
        
      bc->next = pfd->ls[pf]->embedded_bc;
      pfd->ls[pf]->embedded_bc = bc;
      
      bc->bc_input_id = ibc;
      
      /* if ( BC_Types[ibc].BC_Name == PF_CAPILLARY_BC ) found_capillary = TRUE; */

    }
  }  /* END for (ibc = 0; ibc < Num_BC; ibc++)			     */

  /*   if (!found_capillary) { */
  /*     WH(-1,"No LS_CAPILLARY condition specified on the level set surface.\n"); */
  /*     WH(-1,"To add surface tension add 'BC = LS_CAPILLARY LS 0' in the input file.\n"); */
  /*   } */
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int
find_bc_unk_offset(struct Boundary_Condition *bc, int curr_mat,
		   int inode, int p, int *retn_matIndex,
		   VARIABLE_DESCRIPTION_STRUCT **vd_ptr)

     /*********************************************************
      *
      * find_bc_unk_offset():
      *
      *   Finds the offset (and therefore the unknown id )
      *   for applying a boundary condition at a node
      *
      *  Input
      *    bc -> boundary condition structure
      *    curr_mat -> current material ID of the material that
      *                will receive the boundary conditon
      *                If this field is -2, then the  
      *                boundary condition is applied to the
      *                unknown in the first material containing
      *                the requesite variable type.
      *    inode -> current node number
      *    p     -> Current value of the vector, if this is a
      *             vector boundary condition
      *  
      * Return
      *
      *  Offset of the unknown in the solution vector.
      *  -1, if this boundary condition is not to be applied
      *  at this node.
      *  
      * Output            
      *
      *  retn_matIndex = Returns the material index of the
      *                  the unknown on which this boundary
      *                  condition will be applied. (can not 
      *                  be equal to -1, even if the underlying
      *                  variable type has a matID value of -1.
      *  vd_ptr -> pointer to the variable descriptions struct
      *            corresponding to the unknown on which this
      *            boundary condition is applied
      *********************************************************/
{
  int matIndex = -1, ndofs_ieqn, offset, imat;
  struct BC_descriptions *bc_desc = bc->desc;
  int ieqn = bc_desc->equation + p;
  NODE_INFO_STRUCT *node = Nodes[inode];
  NODAL_VARS_STRUCT *nv = node->Nodal_Vars_Info[pg->imtrx];
  PROBLEM_DESCRIPTION_STRUCT *pd_curr = 0;
  int bc_name = bc_desc->BC_Name;
  /*
   * Not at an intersection:
   *    For boundary conditions which refer to rotated
   *    coordinates, assign an equation number to them
   *    base upon a (normal, tang1, tang2) = ( x, y, z)
   *    mapping.
   */
  if ((ieqn >= R_MESH_NORMAL) && (ieqn <= R_MESH_TANG2)) {
    ieqn = ieqn - R_MESH_NORMAL + R_MESH1;
  }
  if ((ieqn >= R_SOLID_NORMAL) && (ieqn <= R_SOLID_TANG2)) {
    ieqn = ieqn - R_SOLID_NORMAL + R_SOLID1;
  }
  if ((ieqn >= R_MOM_NORMAL) && (ieqn <= R_MOM_TANG2)) {
    ieqn = ieqn - R_MOM_NORMAL + R_MOMENTUM1;
  }
  /*
   * If we haven't assigned ieqn by this point for rotated bc's
   * we have made an error
   */
  /*
   * If the equation has zero degrees of freedom at this node
   * return. Otherwise, get the number of variable description
   * structures with this variable type. Don't count more
   * than one species unknown per material at this point.
   */
  if ((ndofs_ieqn = get_nv_ndofs_modMF(nv, ieqn)) <= 0) {
    return -1;
  }
  /*
   *  Nominally, you apply the boundary condition on the equation
   *  for the degree of freedom corresponding to the current
   *  material of the element.
   *  (if there is more than one dof for the specified varType at
   *   the current node).
   *  An input material index of -2 means the first material at the node
   *  while -1 matid is not allowed as input to this function.
   *  Thus, below if the input material id is -2, we find the
   *  first material at the current node that contains the
   *  ieqn variable type. In any case, we set matIndex and
   *  pd_curr appropriately duing this section.
   */
  if (curr_mat == -2) {
    for (imat = 0; imat <  node->Mat_List.Length; imat++) {
      matIndex = node->Mat_List.List[imat];
      pd_curr = pd_glob[matIndex];
      if (pd_curr->e[pg->imtrx][ieqn]) break;
    }
  } else {
    matIndex = curr_mat;
    pd_curr = pd_glob[matIndex];
  }
 
  /*
   * Now, let's handle the case where the current boundary
   * condition is applied on an equation type for which there
   * isn't a valid interpolation within the current material.
   *
   * If this condition should be applied only from the side of
   * the interface where the variable is defined, return -1.
   * For SINGLE_PHASE boundary conditions, this is not an error.
   * It is a feature. We will just not collect a contribution
   * from this side of the boundary.
   *
   * Note, there are some single phase boundary
   * conditions, such as CAPILLARY, which collect contributions
   * from both sides of an internal boundary. However, 
   * the corresponding equation, in this case velocity, 
   * exists on both sides of the boundary.
   */
  if (!pd_curr->e[pg->imtrx][ieqn]) {
    /*
     *  Return here if we are on part of the interface contained
     *  in a material that doesn't have a valid volumetric
     *  interpolation for the pertinent variable. Note, at
     *  this point we have determined that the variable exists
     *  at the node. However, for SINGLE_PHASE boundary 
     *  conditions, we wish to discard these cases.
     */
    if (bc_desc->i_apply == SINGLE_PHASE) {
      return -1;
    }
    /*
     * If the boundary condition is a cross-phase boundary
     * condition and the equation isn't in the material,
     * then we look for the equation in the material
     * across the interface from the current material.
     * We do this by looking up what materials are on
     * the other side of the interface. This was previously
     * calculated and storred in the BC structure.
     */
    else {
      matIndex = bc->BC_matrl_index_1;
      if (matIndex == curr_mat) {
	matIndex = bc->BC_matrl_index_2;
      }
      /*
       * Just to dot the eyes, make sure that the selected matIndex
       * is present at this node
       */
      if (in_list(matIndex, 0, node->Mat_List.Length, 
		  node->Mat_List.List) == -1) {
	EH(GOMA_ERROR,"TROUBLE");
      }
    }
  }

  /*
   * Next handle DIVTIE boundary conditions where the
   * boundary condition is applied as a strong boundary condition
   * on the material with the high material index
   */
  if (ndofs_ieqn > 1) {
 
    /*
     *  For TIE boundary conditions involving interfaces
     *  with discontinuous variables, choose to replace the
     *  normal continuity equation in the material having
     *  the larger material index with the strongly integrated
     *  TIE Dirichlet condition.
     */
    if (bc_name == CONT_TANG_VEL_BC ||
	bc_name == CONT_NORM_VEL_BC ||
	bc_name == DISCONTINUOUS_VELO_BC ||
	bc_name == VL_EQUIL_BC ||
	bc_name == VL_POLY_BC ||
	bc_name == SDC_STEFANFLOW_BC) {
      matIndex = MAX(bc->BC_matrl_index_1, bc->BC_matrl_index_2);
    }
    /*
     * For DVVSIG (Discontinuous Variable Volumetric Surface
     * Integral Galerkin) boundary conditions applied to
     * interfaces with discontinuous variables, choose to
     * add the surface integral boundary condition to
     * the dofs from the material with lowest number material id,
     * irrespective of what material we are currently in.
     */
    else if (bc_name == KINEMATIC_SPECIES_BC ||
	     bc_name == CAPILLARY_BC ||
	     bc_name == CAP_REPULSE_BC ||
	     bc_name == CAP_REPULSE_ROLL_BC ||
	     bc_name == CAP_REPULSE_USER_BC ||
	     bc_name == CAP_REPULSE_TABLE_BC ||
	     bc_name == CAPILLARY_TABLE_BC ||
	     bc_name == CAP_RECOIL_PRESS_BC) {
      matIndex = MIN(bc->BC_matrl_index_1, bc->BC_matrl_index_2);
    }
  }
  *retn_matIndex = matIndex;
  /*
   * OK, now that we know the matIndex of the equation
   * let's find the offset
   */

  if (ieqn == MASS_FRACTION) {
    offset = get_nodal_unknown_offset(nv, ieqn, matIndex, bc->species_eq,
				      vd_ptr);
  } else {
    offset = get_nodal_unknown_offset(nv, ieqn, matIndex, 0, vd_ptr);
  }
  if (offset < 0) {
    return -1;
  }
  return offset;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void
setup_Elem_BC(struct elem_side_bc_struct **elem_side_bc,
	      struct Boundary_Condition  *bc_type,
	      int ibc, int num_nodes_on_side, int ielem,
	      int local_ss_node_list[], Exo_DB *exo)

     /*********************************************************************
      *
      * setup_Elem_BC()
      *
      *  This routine adds a boundary condition to an element side_bc_struct
      *  It fills in all pertinent fields in the  elem_side_bc_struct
      *  structure. It is called recursively, because elem_side_bc_struct
      *  is formulated as a linked list.
      *********************************************************************/
{
  int 		id_local_elem_coord[num_nodes_on_side];
  struct elem_side_bc_struct *side = *elem_side_bc;
  char err_msg[MAX_CHAR_IN_INPUT];
  /* Check to see if we are on a new side of the element that doesn't have
   * a boundary condition specified on it */

  if (side == NULL) {
    side = alloc_struct_1(struct elem_side_bc_struct, 1);
    side->ielem = ielem;
    *elem_side_bc = side;
    side->id_side = find_id_side_BC(ielem, num_nodes_on_side,
				    local_ss_node_list, ibc, id_local_elem_coord, exo);
    side->num_nodes_on_side = num_nodes_on_side;
    side->BC_applied = bc_type->BC_Name;
    side->local_node_id = alloc_int_1(num_nodes_on_side, INT_NOINIT);
    side->local_elem_node_id = alloc_int_1(num_nodes_on_side, INT_NOINIT);
    for (int i = 0; i < num_nodes_on_side ; i++) {
      side->local_node_id[i] = local_ss_node_list[i];
      side->local_elem_node_id[i] = id_local_elem_coord[i];
    }
    side->BC_input_id = alloc_int_1(2, -1);
    side->BC_input_id[0] = ibc;
    side->Num_BC = 1;
    side->next_side_bc = NULL;
    /*
     * Determine how many materials comprise this side.
     */
    elem_side_matrl_list(side);
    vcrr_determination(side, bc_type, ielem, ibc);
  } else {
    if (same_side(local_ss_node_list, side->local_node_id, 
		  num_nodes_on_side)) {
      /* More than one surface integral BC on the side of the element -
       * Add to the existing structure 
       */
      side->BC_applied += bc_type->BC_Name;
      realloc_int_1(&(side->BC_input_id), side->Num_BC + 2,
		    side->Num_BC + 1);
      side->BC_input_id[side->Num_BC] = ibc;
      side->Num_BC++;
      side->BC_input_id[side->Num_BC] = -1;
      if (side->Num_BC > MAX_BC_PER_SIDE) {
	sprintf(err_msg, 
		"MAX_BC_PER_SIDE(%d) exceeded - elem [%d], BC %d",
		MAX_BC_PER_SIDE, ielem, ibc);
	EH(GOMA_ERROR, err_msg);
      }
      vcrr_determination(side, bc_type, ielem, ibc);
    } else {
      /*
       * Go to the next link for the current element in a
       * recursive loop
       */  
      setup_Elem_BC(&(side->next_side_bc), bc_type, ibc, num_nodes_on_side, 
		    ielem, local_ss_node_list, exo);
    }
  }
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static struct elem_edge_bc_struct *
setup_Elem_Edge_BC (struct elem_edge_bc_struct **elem_edge_bc,
		    struct Boundary_Condition  *bc_type,
		    int ibc, 
		    int num_nodes_on_edge,
                    int ipin,
		    int shared,
		    int ielem, 
		    int local_edge_node_list[],
		    Exo_DB *exo)
{
  int 		i;
  int 		id_local_elem_coord[MAX_NODES_PER_SIDE];
  int           param_dir;
  struct elem_edge_bc_struct *this_edge_bc = *elem_edge_bc;
  /* Check to see if we are on a new edge of the element that doesn't have
     a boundary condition specified on it */

  if (*elem_edge_bc == NULL) 
    {

      *elem_edge_bc = alloc_struct_1(struct elem_edge_bc_struct, 1);

      /* Initialize the structure */
      (*elem_edge_bc)->ielem = ielem;
      (*elem_edge_bc)->id_edge = find_id_edge(ielem, num_nodes_on_edge,
					      local_edge_node_list,
					      id_local_elem_coord, 
					      &param_dir, exo);
      (*elem_edge_bc)->num_nodes_on_edge = num_nodes_on_edge;
      (*elem_edge_bc)->ipin = ipin;
      (*elem_edge_bc)->shared = shared;
      (*elem_edge_bc)->BC_applied = bc_type->BC_Name;
      for (i = 0; i < num_nodes_on_edge ; i++) {
	(*elem_edge_bc)->local_node_id[i] = local_edge_node_list[i];
	(*elem_edge_bc)->edge_elem_node_id[i] = id_local_elem_coord[i];
      }
      (*elem_edge_bc)->BC_input_id[0] = ibc;
      for (i = 1; i < MAX_BC_PER_SIDE ; i++)
	(*elem_edge_bc)->BC_input_id[i] = -1;
      (*elem_edge_bc)->next_edge_bc = NULL;
	
      (*elem_edge_bc)->elem_side_bc_1 = NULL;
      (*elem_edge_bc)->elem_side_bc_2 = NULL;
      this_edge_bc = *elem_edge_bc;
    } 
  else 
    {

      /* Go to the next link for the current element in a recursive loop */

      this_edge_bc = setup_Elem_Edge_BC(&((*elem_edge_bc)->next_edge_bc), 
					bc_type, ibc, num_nodes_on_edge, 
					ipin, shared, ielem,
					local_edge_node_list, exo);
    }

  return this_edge_bc;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static int
same_side (int ilist1[], int ilist2[], int list_length)

     /* 
      *  BOOLEAN function that returns 0 if the two integer lists don't 
      *  contain the same values and 1 if they do.
      *
      *  Currently, it is assumed that the two lists are ordered in the same way.
      */
{
  int i;
  for (i = 0; i< list_length; i++) {
    if (ilist1[i] != (int) ilist2[i]) return (0);
  }
  return (1);
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void
initialize_Boundary_Condition (struct Boundary_Condition *bc_ptr)

     /************************************************************************
      *
      * initialize_Boundary_Condition():
      *
      *    Sets the initial conditions for the Boundary_COndition structure
      *  before being filled in rd_bc_specs().
      *
      *************************************************************************/
{
  if (bc_ptr == NULL) {
    EH(GOMA_ERROR, "initialize_Boundary_Condition FATAL ERROR");
  }
  
  /*
   * set all ints to zero and pointers to NULL
   */
  (void) memset(bc_ptr, 0, sizeof(struct Boundary_Condition));

  bc_ptr->BC_Data_Int[0] = -1;

  /*
   *  The default is to assume no relaxation with this dirichlet condition
   */
  bc_ptr->BC_relax = -1.0;

  /*
   * Initially there are no user constants for this bc.
   * If, through a read_constants() call you allocate space for
   * the "double *u_BC" member, then update len_u_BC to assist
   * transporting this list across processors.
   */
  bc_ptr->len_u_BC = 0;
  bc_ptr->max_DFlt = 0;

  /*
   * Also, the majority of boundary conditions may rely upon
   * one of Rich's statically declared BC_Descriptions. However,
   * some new fancy BCs seem to want to make their own description.
   * These new dynamically allocated versions contain information that
   * may be required on other processors where the BC will actually
   * be applied.
   *  Higher values of index_dad mean this BC
   * points to a dynamically-allocated 
   * BC_Description... 
   */
  bc_ptr->index_dad = -1;

  /*
   * Set some quantities to -1 to indicate that they are unspecificed
   */
  bc_ptr->BC_ID = -1;
  bc_ptr->BC_ID2 = -1; 
  bc_ptr->BC_ID3 = -1;
  bc_ptr->BC_matrl_index_1 = -1;
  bc_ptr->BC_matrl_index_2 = -1;
  bc_ptr->BC_matrl_index_3 = -1;
  bc_ptr->BC_matrl_index_4 = -1;
  bc_ptr->BC_EBID_Apply = -1;
  bc_ptr->species_eq = -1;
  bc_ptr->equation = -1;
  bc_ptr->matrix = 0;
  
  /*
   *  Since the table structs are dynamically allocated. This index is used
   *  to record  which of the table structs goes  with this BC.
   *  Used when reconstructing the BC on another processor.
   *  The value of negative one indicates that there is no corresponding
   *  table to go with this boundary condition. The Table pointer
   *  has also been set to NULL.
   */
  bc_ptr->table_index = -1;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int
find_id_side(const int ielem,			/* element index number */
	     const int num_nodes_on_side,	/* guess what?          (in) */
	     const int local_ss_node_list[],	/* nodes, ordered	(in) */
	     int id_local_elem_coord[],	/* Local node numbers of the side (out) */
	     const Exo_DB *exo)		/* ptr to FE db         (in) */
{

  /* TAB certifies that this function conforms to the exo/patran side numbering convention
   * 11/9/98.  The return value is exo/patran side number ! */

  int ielem_type;
  int num_local_nodes;
  int ielem_dim;
  int iconnect_ptr;
  int i;
  double sum;
  
  char err_msg[MAX_CHAR_ERR_MSG];  

  static char *yo = "find_id_side";

  /*----------------------------Start Execution----------------------------*/

  /* Find out what type of element this is */   

  ielem_type = Elem_Type(exo, ielem);
    
  /* Find out how many local basis functions there are */ 

  num_local_nodes = elem_info(NNODES, ielem_type);
    
  /* Find out the physical dimension of the element */  

  ielem_dim = elem_info(NDIM, ielem_type);

  /* find the pointer the beginning of this element's connectivity list */

  iconnect_ptr = exo->elem_node_pntr[ielem];

  /* Find the local element coordinates - note, this algorithm keeps the 
   * ordering in local_ss_node_list - might be necessary 
   */

  for (i = 0; i < num_nodes_on_side; i++)
    {
      if ( ( id_local_elem_coord[i] = 
	     in_list (local_ss_node_list[i], 0, num_local_nodes, 
		      &exo->elem_node_list[iconnect_ptr]) ) == -1)
	{
	  sr = sprintf(err_msg, "%s: lost SS node [%d] = %d  %d  %d",
		       yo, i, local_ss_node_list[i], num_nodes_on_side, ielem );
 	  EH(GOMA_ERROR, err_msg); 
	}
    }

  switch (ielem_dim) 
    {
    case 3:
      for (i = 0, sum = 0.0; i < num_nodes_on_side; i++)
	sum += shape( 0.0, 0.0, 1.0, ielem_type, PSI, id_local_elem_coord[i]);
      if (sum > 0.999) return (6);
      for (i = 0, sum = 0.0; i < num_nodes_on_side; i++)
	sum += shape( 0.0, 0.0,-1.0, ielem_type, PSI, id_local_elem_coord[i]);
      if (sum > 0.999) return (5);
      /* fall through */
    case 2:
      /* newly added for triangles */
      for (i = 0, sum = 0.0; i < num_nodes_on_side; i++)
	sum += shape( 0.5, 0.0, 0.0, ielem_type, PSI, id_local_elem_coord[i]);
      if (sum > 0.999) return (3);
	   
      for (i = 0, sum = 0.0; i < num_nodes_on_side; i++)
	sum += shape( 0.5, 0.5, 0.0, ielem_type, PSI, id_local_elem_coord[i]);
      if (sum > 0.999) return (1);

      for (i = 0, sum = 0.0; i < num_nodes_on_side; i++)
	sum += shape( 0.0, 0.5, 0.0, ielem_type, PSI, id_local_elem_coord[i]);
      if (sum > 0.999) return (2);

      for (i = 0, sum = 0.0; i < num_nodes_on_side; i++)
	sum += shape( -1.0, 0.0, 0.0, ielem_type, PSI, id_local_elem_coord[i]);
      if (sum > 0.999) return (4);
      for (i = 0, sum = 0.0; i < num_nodes_on_side; i++)
	sum += shape( 0.0,1.0, 0.0, ielem_type, PSI, id_local_elem_coord[i]);
      if (sum > 0.999) return (3);

      for (i = 0, sum = 0.0; i < num_nodes_on_side; i++)
	sum += shape( 1.0, 0.0, 0.0, ielem_type, PSI, id_local_elem_coord[i]);
      if (sum > 0.999) return (2);
      for (i = 0, sum = 0.0; i < num_nodes_on_side; i++)
	sum += shape(0.0, -1.0, 0.0, ielem_type, PSI, id_local_elem_coord[i]);
      if (sum > 0.999) return (1);
      break;
    case 1:
      //  Side 1 is next to node number 0 at s = -1 on the left side.
      //  Side 2 is next to node number 1 at s = +1 on the right side. 
      for (i = 0, sum = 0.0; i < num_nodes_on_side; i++)
	sum += shape( 1.0, 0.0, 0.0, ielem_type, PSI, id_local_elem_coord[i]);
      if (sum > 0.999) return (2);
      for (i = 0, sum = 0.0; i < num_nodes_on_side; i++)
	sum += shape(-1.0, 0.0, 0.0, ielem_type, PSI, id_local_elem_coord[i]);
      if (sum > 0.999) return (1);
      break;
    } /* END switch ielem_dim */

  /* An error condition has occurred, if here */

  sr = sprintf(err_msg, 
	       "%s: problem for elem (%d), dimension %d, %d nodes on side.", 
	       yo, ielem+1, ielem_dim, num_nodes_on_side);

  EH(GOMA_ERROR, err_msg);

  return (-1);
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int
find_id_side_BC(const int ielem,		/* element index number */
		const int num_nodes_on_side,	/* guess what?          (in) */
		const int local_ss_node_list[],	/* nodes, ordered	(in) */
		const int ibc,                  /* BC index             (in) */
		int id_local_elem_coord[],	/* Local node numbers of the side (out) */
		const Exo_DB *exo)		/* ptr to FE db         (in) */
{
  int sideid, ielem_type;

  /* Run the usual routine */
  sideid = find_id_side(ielem, num_nodes_on_side,
			local_ss_node_list, id_local_elem_coord, exo);

  /* Find out what type of element this is */   
  ielem_type = Elem_Type(exo, ielem);

  /* If we're working with tetrahedral elements, re-work the side id */
  if ( ielem_type == LINEAR_TET ) {
    bool nodes[4] = {false,false,false,false};
    for (int i = 0; i < 4; i++) {
      int node_id = exo->elem_node_list[exo->elem_node_pntr[ielem]+i];
      for (int j = 0; j < num_nodes_on_side; j++) {
        if (node_id == local_ss_node_list[j]) {
          nodes[i] = true;
        }
      }
    }
    if (nodes[0] && nodes[1] && nodes[3]) {
      return 1;
    } else if (nodes[1] && nodes[2] && nodes[3]) {
      return 2;
    } else if (nodes[0] && nodes[2] && nodes[3]) {
      return 3;
    } else if (nodes[0] && nodes[1] && nodes[2]) {
      return 4;
    } else {
      EH(GOMA_ERROR, "Unknown tet layout");
    }
  }

  return(sideid);    
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int 
find_id_side_SS(const int ielem,    /* ielem - element index number     (in) */
		const int iss,      /* iss - sideset index number       (in) */
		const Exo_DB *exo)  /* exo - pointer to FE db           (in) */
/* 
 * Finds the side id based on a side set number and element number 
 * Author:  Scott A. Roberts, 1514   November 2, 2011
 */
{
  int i;
  int nsides = exo->ss_num_sides[iss];   // Number of sides in sideset
  int eindx = exo->ss_elem_index[iss];   // Index to begnning of sideset
  for ( i = 0; i < nsides; i++) {
    if ( ielem == exo->ss_elem_list[eindx+i] ) {
      return(exo->ss_side_list[eindx+i]);
    }
  }
  return(-1);
}  /* End of find_id_side_SS  */

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/


/* print_setup_Surf_BC() - dump out bc data
 *
 */

void
print_setup_Surf_BC(struct elem_side_bc_struct *First_Elem_Side_BC_Array[ ])
{
  int			ielem;

  if (ProcID == 0) {
    (void) printf(
		  "\t\t PRINTOUT OF SURFACE BOUNDARY INTEGRAL SET-UP\n");
  }
  (void) printf ("\n\nSurface Integrals on Processor %d:\n", ProcID);
  (void) printf("Global_Elem_Num ID_side Num_Nodes_On_Side BC_applied  ");
  (void) printf(" BC_input_id   |  Global_node_id\'s ");
  (void) printf("       |  Local_Elem_node_ID\'s\n");

  for (ielem = 0; ielem < Num_Internal_Elems; ielem++) {
    if (First_Elem_Side_BC_Array[ielem] != NULL)
      print_Elem_Surf_BC (ielem, First_Elem_Side_BC_Array[ielem]);
  }

  return;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void
print_Elem_Surf_BC (int ielem,
		    struct elem_side_bc_struct *elem_side_bc)
{
  int i, ibc = 0;

  while (elem_side_bc->BC_input_id[ibc] != -1) {
    (void) printf (
		   "   %5d       %5d        %5d         %9d     %5d        |    ", 
		   ielem,
		   (int) elem_side_bc->id_side,
		   (int) elem_side_bc->num_nodes_on_side,
		   elem_side_bc->BC_applied,
		   (int) elem_side_bc->BC_input_id[ibc]);
    for (i = 0; i < (int) elem_side_bc->num_nodes_on_side ; i++)
      (void) printf (" %d  ",  (int)elem_side_bc->local_node_id[i]);
    (void) printf ("   |  ");
    for (i = 0; i < (int) elem_side_bc->num_nodes_on_side ; i++)
      (void) printf (" %2d ", (int)elem_side_bc->local_elem_node_id[i]);
    (void) printf ("\n");
    ibc++;
  }
  if (elem_side_bc->next_side_bc != NULL)
    print_Elem_Surf_BC (ielem, elem_side_bc->next_side_bc);
}
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/

struct BC_descriptions *
alloc_BC_description (struct BC_descriptions *old_ptr)

     /********************************************************************
      *
      * alloc_BD_description():
      *
      * Function to allocate space for a new BC_description structure 
      * and then initialize it to the values of a previously created
      * BC_description structure, old_ptr.
      *
      * Return
      * --------
      *   This funtion returns the address of the newly malloced
      *   structure.
      ********************************************************************/
{
  struct BC_descriptions *tmp_ptr;
  tmp_ptr = alloc_struct_1(struct BC_descriptions, 1);
  (void) memcpy((void *) tmp_ptr, (void *) old_ptr,
		sizeof(struct BC_descriptions));
  return tmp_ptr;
}
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*
 * Convert the integer flag for rotation condition equation type into a string
 * for more informative diagnostics.
 *
 * Created: 1999/01/28 08:06 MST pasacki@sandia.gov
 */

char *
rot_eq_type2string(const int eq_type)
{
  static char *values[] = { "MOM", "MESH", "(bad)" };

  if ( eq_type == R_MESH1 )
    {
      return(values[1]);
    }
  else if ( eq_type == R_MOMENTUM1 )
    {
      return(values[0]);
    }

  /*
   * Ought not to get here.
   */

  return(values[2]);
}
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*
 * Convert the integer flag for rotation condition topology type into a string
 * for more informative diagnostics.
 *
 * Created: 1999/01/28 09:01 MST pasacki@sandia.gov
 */

char *
rotopology_type2string(const int type)
{
  static char *values[] = { "SURFACE", "EDGE", "VERTEX", "VOLUME", "(bad)" };

  if ( type == FACE )
    {
      return(values[0]);
    }
  else if ( type == CURVE )
    {
      return(values[1]);
    }
  else if ( type == VERTEX )
    {
      return(values[2]);
    }
  else if ( type == BODY )
    {
      return(values[3]);
    }
  
  return(values[4]);
}
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/

static void
elem_side_matrl_list(struct elem_side_bc_struct *side)

     /*********************************************************************
      *
      *  elem_side_matrl_list()
      *
      *  This routine develops the list of materials that are adjacent
      *  to an element side. It does this by looking at all of the nodes
      *  that comprise a side. Then, a material is deemed part of the
      *  side if it is included in all nodes that are part of the side.
      *********************************************************************/
{
  int i, j, *matrl_list, num_nodes, *node_list;
  UMI_LIST_STRUCT *mlist;
  NODE_INFO_STRUCT *node;
  num_nodes = side->num_nodes_on_side;
  node_list = side->local_node_id;
  matrl_list = alloc_int_1(upd->Num_Mat, 0);
  for (i = 0; i < num_nodes; i++) {
    node = Nodes[node_list[i]];
    mlist = &(node->Mat_List);
    for (j = 0; j < mlist->Length; j++) {
#ifdef DEBUG_IGNORE_ELEMENT_BLOCK_CAPABILITY
      if (mlist->List[j] < 0) {
        fprintf(stderr,"node list contains negative mn number\n");
        EH(GOMA_ERROR, "logic error with ignored element block");
      }
#endif
      matrl_list[mlist->List[j]]++;
    }
  }
  for (i = 0, j = 0; i < upd->Num_Mat; i++) {
    if (matrl_list[i] == num_nodes) {
      side->MatID_List[j] = i;
      j++;
    }
  }
  side ->Num_MatID = j;
  safer_free((void **) &matrl_list);
}
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/

static void
vcrr_determination(struct elem_side_bc_struct *side,
		   struct Boundary_Condition *bc, int ielem, int ibc)

     /*********************************************************************
      *
      * vcrr_determination()
      *
      ********************************************************************/
{
  int i, node_num, lnn, mn, p, ebn, mn_low, nf, index_eq, matID_apply;
  int eqn;
  struct BC_descriptions *bc_desc = bc->desc;
  NODE_INFO_STRUCT *node_info;
  VARIABLE_DESCRIPTION_STRUCT *vd;
  NODAL_RESID_WKSP_STRUCT *wksp;


  /*
   * First determine if this boundary condition requires the addition
   * of an entry into this list
   */
  if (bc->DV_Indexing_Type == DVI_DVVSIG) {
    if (side->Num_MatID > 1) {
      /*
       * Find the current material
       */
      ebn = find_elemblock_index(ielem, EXO_ptr);
      mn  = Matilda[ebn];
      mn_low = MIN(bc->BC_matrl_index_1, bc->BC_matrl_index_2);
      if (mn_low != mn) {
        for (i = 0, nf = TRUE; i < side->Num_MatID && nf ; i++) {
          if (side->MatID_List[i] == mn_low) {
            nf = FALSE;
	  }
	}
	if (nf) {
	  EH(GOMA_ERROR,"Problems with material mappings");
	}
        /*
	 * If we are here, then we need to add an entry for
	 * every node on the surface
	 */
	for (i = 0; i < side->num_nodes_on_side; i++) {
	  node_num =  side->local_node_id[i];
	  lnn = side->local_elem_node_id[i];
	  node_info = Nodes[node_num];
          if (!(node_info->Resid_Wksp)) {
            node_info->Resid_Wksp = nodal_resid_wksp_alloc();
	  }
	  wksp = node_info->Resid_Wksp;
	  for (p = 0; p < bc_desc->vector; p++) {
	    index_eq = bc_eqn_index(lnn, node_num, ibc, mn,
				    p, &eqn, &matID_apply, &vd);
            if (index_eq >= 0) {
	      vcrr_add(&(wksp->Vcrr), eqn, mn, mn_low);
	    }
	  }
	}
      }
    }
  }
}
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/

void
set_up_BC_connectivity(void)

     /************************************************************
      *
      * set_up_BC_connectivity()
      *
      *  This routine consists of a loop over the boundary conditions.
      *  In this loop boundary conditions get information that they
      *  need from other boundary conditions and other information
      *  in the problem set.
      *
      *  SDC_KIN_SF_BC:
      *      These stefan flow boundary conditions get the bc_id 
      *  number for the associated bc that handles the determination
      *  of the surface reaction rate on the surface.
      ************************************************************/
{
  int ibc, jbc, ifound = TRUE;
  BOUNDARY_CONDITION_STRUCT *bc, *bc2;
  struct BC_descriptions *bc_desc, *bc_desc2;
  for (ibc = 0; ibc < Num_BC; ibc++) {
    bc = BC_Types + ibc;
    bc_desc = bc->desc;
    switch (bc_desc->BC_Name) {
    case SDC_KIN_SF_BC:
      for (jbc = 0, ifound = FALSE; jbc < Num_BC && !ifound; jbc++) {
	bc2 = BC_Types + jbc;
	bc_desc2 = bc2->desc;
	if (bc2->BC_ID == bc->BC_ID) {
	  if (bc_desc2->BC_Name == bc->BC_Data_Int[1]) {
	    ifound = TRUE;
	    bc->BC_Data_Int[2] = jbc;
	  }
	}
      }
      if (!ifound) {
	EH(GOMA_ERROR,"set_up_BC_connectivity failure");
      }
      break;
    default:
      break;
    }
  }
}

/*************************************************************************/
/*************************************************************************/
/*************************************************************************/

/**
 * Exchange needed boundary condition information that is only
 * present on one processor but needed on all or some of the other
 * processors.
 *
 * Runs in serial as well so calculations that are shared can be moved
 * to this exchange
 *
 * Called from setup_problem(), only called before solve
 * Returns 0 if success
 *
 * Return -1 if error
 */
int exchange_bc_info(void)
{
  int ibc;
  int error = 0;

  /* loop over boundary conditions */
  for (ibc = 0; ibc < Num_BC; ibc++) {
    /* check if BC needs special exchange information */
    switch (BC_Types[ibc].BC_Name) {
    case VELO_SLIP_BC:
    case VELO_SLIP_ROT_BC:
    case VELO_SLIP_FLUID_BC:
    case VELO_SLIP_ROT_FLUID_BC:
    case ROLL_FLUID_BC:
    case AIR_FILM_BC:
    case AIR_FILM_ROT_BC:
      exchange_fvelo_slip_bc_info(ibc);
      break;
    default:
      break;
    }

  } /* end loop over BC_Types */

  return error;
}
/************************************************************************/
/*			END of mm_bc.c			                */
/************************************************************************/
