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
 * These routines setup the MPI datatypes used to communicate information
 * between this processor and its neighbors. Both node based data (eg fill eqn)
 * and full degree of freedom based information is considered. The saved
 * data types are in the Comm_Ex structures. These structures have handy
 * pointer variables to lists of nodes and dofs (in the local processor 
 * numbering scheme) that are aliases to appropriate spots in the big lists
 * defined globally in rf_fem.h
 *
 * The PARALLEL C preprocessor variable is used to protect inadvertent access
 * to MPI routines when they are not linked in. The MPI datatypes embedded
 * in the Comm_Ex structure are replaced with innocuous int's and otherwise
 * the dpi data is set up to work for the serial case as much as possible.
 *
 * Created: 1997/08/16 14:16 MDT pasacki@sandia.gov
 *
 * Revised: 
 */


#include <stdlib.h>
#include <stdio.h>

#include "std.h"
#include "el_elm_info.h"
#include "rf_allo.h"
#include "rf_fem_const.h"
#include "rf_mp.h"
#include "rf_vars_const.h"
#include "mm_as_structs.h"
#include "mm_as.h"
#include "mm_eh.h"
#include "exo_struct.h"
#include "dpi.h"
#include "dp_types.h"
#include "dp_map_comm_vec.h"
#include "dp_utils.h"
#include "el_elm_info.h"
#include "rf_node_const.h"
#include "rf_util.h"

#ifdef PARALLEL
#include "mpi.h"
#endif

#define GOMA_DP_MAP_COMM_VEC_C

/*
 * This is the single place where these variables are defined. Usually they
 * are declared extern for global use in rf_fem.h
 */

int *node_index_fill = NULL;
int *fill_node_list = NULL;
int num_fill_unknowns = 0;
int owned_fill_unknowns = 0;
int internal_fill_unknowns = 0;
int boundary_fill_unknowns = 0;
int external_fill_unknowns = 0;

/*
 *      more singly defined pointers
 *
 *  Pointers to array lists being transported, 
 *  and the lists of integers to be transported,
 *  likewise declared extern for global use in rf_fem.h
 */

int *ptr_node_recv = NULL;
int *ptr_fill_node_recv = NULL;
int **ptr_dof_send = NULL;
int **list_dof_send = NULL;
int *ptr_node_send = NULL;
int *list_node_send = NULL;
int *ptr_fill_node_send = NULL;
int *list_fill_node_send = NULL;

static void build_node_recv_indeces /* dp_map_comm_vec.c */
(Exo_DB *,		/* exo */
       Dpi *);			/* dpi */

/*********************************************************************************/
/*********************************************************************************/
/*********************************************************************************/

void
setup_nodal_comm_map(Exo_DB *exo, Dpi *dpi, Comm_Ex **cx)

    /****************************************************************************
     *
     * setup_nodal_comm_map():
     *
     *  This routine sets up the communications pattern between nodes needed
     *  for the owning node to ghost node communications. It sets up a list
     *  of nodes that this processor needs information from neighboring nodes,
     *  and it sets up a list of nodes that this processor must send to 
     *  its neighbors processors.
     ****************************************************************************/
{
  int i, p, gnn, lnn, local_node_num;
  int *list_gnode_recv = NULL;
  static char yo[] = "setup_nodal_comm_map ERROR:";
  COMM_NP_STRUCT *np_base=NULL, *np_ptr;
  int imtrx;
  /*
   * Fill in the global vectors ptr_node_recv[] and
   * list_node_recv[]
   */
  build_node_recv_indeces(exo, dpi);
  
  /*
   * For more convenience, set up the Comm_Ex structures.
   *    p is the index of the neighboring processor, 
   *    neighbor[p] is its processor number
   */
  for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
    for (p = 0; p < dpi->num_neighbors; p++) {
      /*
       *  Processor number of the neighboring processor
       *  (0 to Nproc-1)
       */
      cx[imtrx][p].neighbor_name      = dpi->neighbor[p];
      /*
       * num_nodes_recv: The number of nodes to receive from this processor
       */
      cx[imtrx][p].num_nodes_recv     = (ptr_node_recv[p+1] - ptr_node_recv[p]);
    }
  }

  /*
   * Set up a neighboring processor communications structure to
   * facilitate communications between neighboring processors
   */
#ifdef PARALLEL
  if (dpi->num_neighbors > 0) {
    np_base = alloc_struct_1(COMM_NP_STRUCT, dpi->num_neighbors);
  }
#endif

  for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
  /*
   * Tell each of my neighbors how many nodes that I want from them.
   * Receive from my neighbors how many of my owned nodes I have
   * to send information to them, cx[p].num_nodes_send.
   */
    np_ptr = np_base;
    for (p = 0; p < dpi->num_neighbors; p++) {
      np_ptr->neighbor_ProcID = cx[0][p].neighbor_name;
      np_ptr->send_message_buf = (void *) &(cx[imtrx][p].num_nodes_recv);
      np_ptr->send_message_length = sizeof(int);
      np_ptr->recv_message_buf = (void *) &(cx[imtrx][p].num_nodes_send);
      np_ptr->recv_message_length = sizeof(int);
      np_ptr++;
    }
    exchange_neighbor_proc_info(dpi->num_neighbors, np_base);
  }
  /*
   *  Create space for storage of the node numbers that I must
   *  send to other processors, int list_node_send[].
   *  Also, create a list of pointers into that vector to
   *  denote the start of each neighboring processor's section,
   *  int ptr_node_send[]
   */
  ptr_node_send = alloc_int_1(dpi->num_neighbors + 1, 0);
  for (p = 0; p < dpi->num_neighbors; p++) {
    ptr_node_send[p+1] = ptr_node_send[p] + cx[0][p].num_nodes_send;
  }
  if (ptr_node_send[dpi->num_neighbors] > 0) {
    list_node_send =  alloc_int_1(ptr_node_send[dpi->num_neighbors],
                                  -1);
  }
  /*
   * Fill in the address of the location in list_node_send for each
   * processor's send information into the Comm_ex structure,
   * cx[p].local_nodeces_send
   */
  for (p = 0; p < dpi->num_neighbors; p++) {
    cx[0][p].local_nodeces_send = list_node_send + ptr_node_send[p];
  }

  /*
   *  Create a list of global node numbers that I need from the
   *  my neighboring processors. 
   *  Thus, we can use ptr_node_recv[p] to form individual messages for
   *  each processor.
   */
  if (dpi->num_external_nodes > 0) {
    list_gnode_recv = alloc_int_1(dpi->num_external_nodes, -1);
    for (i =0; i < dpi->num_external_nodes; i++) {
      local_node_num = (dpi->num_internal_nodes +
			dpi->num_boundary_nodes + i);
      list_gnode_recv[i] =  dpi->node_index_global[local_node_num];
    }
  }

  /*
   * Exchange global node number information with the surrounding
   * processors. The nodes that I need to provide information about
   * are now listed in list_node_send[i] after this operation. I must
   * still convert from global node number to processor node number
   * however.
   */
  np_ptr = np_base;
  for (p = 0; p < dpi->num_neighbors; p++) {
    np_ptr->send_message_buf = (void *)
	                       (list_gnode_recv + ptr_node_recv[p]);
    np_ptr->send_message_length = sizeof(int) * cx[0][p].num_nodes_recv;
    np_ptr->recv_message_buf = (void *) cx[0][p].local_nodeces_send;
    np_ptr->recv_message_length = sizeof(int) * cx[0][p].num_nodes_send;
    np_ptr++;
  }
  exchange_neighbor_proc_info(dpi->num_neighbors, np_base);
  safer_free((void **) &list_gnode_recv);

  /*
   *  Map back from global node number to processor node number
   *   -> Note, this is an o(N^2) search algorithm. It may be very
   *            slow as the number of nodes increase.
   */
  for (i = 0; i < ptr_node_send[dpi->num_neighbors]; i++) {
    gnn = list_node_send[i];
    lnn = in_list(gnn, dpi->num_internal_nodes, exo->num_nodes,
		  dpi->node_index_global);
    if (lnn == -1) {
      fprintf(stderr,"%s %d: couldn't find gnn = %d\n", yo, ProcID, gnn);
      EH(GOMA_ERROR,"failed dp mapping");
    }
    list_node_send[i] = lnn;
  }
  
  /*
   * Free memory in the routine
   */
  safer_free((void **) &(np_base));
}
/*********************************************************************************/
/*********************************************************************************/
/*********************************************************************************/

void
exchange_neighbor_proc_info(int num_neighbors, COMM_NP_STRUCT *np_ptr)

    /****************************************************************************
     *
     * exchange_neighbor_proc_info():
     *
     *  This routine performs a canonical exchange operation between
     *  neighboring processors. It uses nonblocking communications and
     *  MPI_BYTE datatypes. Therefore, it is fairly fast. It expects a
     *  vector of COMM_GN_STRUCT structures of length num_neighbors on
     *  each processor as input.
     *
     *  Because it uses MPI_BYTE, it will fail in heterogeneous dp environments.
     *  However, so will aztec, so this is no loss in current generality.
     ****************************************************************************/
{
  static int mtype = 115;
  char *yo = "exchange_neighbor_proc_info ERROR: ";
#ifdef PARALLEL
  int p, retn;
  COMM_NP_STRUCT *np_base = np_ptr;

  /*
   * Post receives for all messages
   */
  for (p = 0; p < num_neighbors; p++) {
    retn = MPI_Irecv((np_ptr->recv_message_buf), np_ptr->recv_message_length,
	             MPI_BYTE, np_ptr->neighbor_ProcID, mtype,
		     MPI_COMM_WORLD, &(np_ptr->recv_request));
    if (retn != MPI_SUCCESS) {
      fprintf(stderr,"%s Proc %d: Irecv to %d failed post: %d\n", yo, ProcID,
	      np_ptr->neighbor_ProcID, retn);
      EH(GOMA_ERROR,"MPI failure");
    }
    np_ptr++;
  }

  /*
   * Send out the data
   */
  np_ptr = np_base;
  for (p = 0; p < num_neighbors; p++) {
    retn = MPI_Isend(np_ptr->send_message_buf, np_ptr->send_message_length,
	             MPI_BYTE, np_ptr->neighbor_ProcID, mtype,
		     MPI_COMM_WORLD, &(np_ptr->send_request));
    if (retn != MPI_SUCCESS) {
      fprintf(stderr,"%s Proc %d: Isend to %d failed post: %d\n", yo, ProcID,
	      np_ptr->neighbor_ProcID, retn);
      EH(GOMA_ERROR,"MPI failure");
    }
    np_ptr++;    
  }

 /*
  *  Wait until all messages are sent and received.
  */
  np_ptr = np_base;
  for (p = 0; p < num_neighbors; p++) {
    retn = MPI_Wait(&(np_ptr->recv_request), &(np_ptr->recv_status));
    if (retn != MPI_SUCCESS) {
      fprintf(stderr,"%s Proc %d: Irecv to %d failed: %d\n", yo, ProcID,
	      np_ptr->neighbor_ProcID, retn);
      EH(GOMA_ERROR,"MPI failure");
    }
    retn = MPI_Wait(&(np_ptr->send_request), &(np_ptr->send_status));
    if (retn != MPI_SUCCESS) {
      fprintf(stderr,"%s Proc %d: Isend to %d failed: %d\n", yo, ProcID,
	      np_ptr->neighbor_ProcID, retn);
      EH(GOMA_ERROR,"MPI failure");
    }    
    np_ptr++;  
  }
#else
  if (num_neighbors > 0) {
    fprintf(stderr,"%s this processor has neighbors but PARALLEL ifdef not on\n",
	    yo); 
    EH(GOMA_ERROR,"MPI failure");
  }
#endif
  mtype++;
  if (mtype > 199) mtype = 115;
}
/*********************************************************************************/
/*********************************************************************************/
/*********************************************************************************/

void
setup_dof_comm_map(Exo_DB *exo, Dpi *dpi, Comm_Ex **cx)

    /*****************************************************************************
     *
     * setup_dof_comm_map():
     *
     *    Set up the communications pattern for exchanging the solution vector
     * from owned nodes to ghost nodes.
     *****************************************************************************/
{
  int i, p, owner, index_owner, dofs_i_want, local_node_number;
  int imtrx;
  int index, base;

  /*
   * Gather the frequency profile of requested dofs for each external
   * processor. Look through every external node.
   */
  for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
    for (p = 0; p < dpi->num_neighbors; p++) {
      cx[imtrx][p].num_dofs_recv = 0;
    }
  }
  for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++)
     {
      for (i = 0; i < dpi->num_external_nodes; i++) 
         {
          local_node_number = dpi->num_internal_nodes +
	                      dpi->num_boundary_nodes + i;
          owner = dpi->node_owner[local_node_number];
          index_owner = in_list(owner, 0, dpi->num_neighbors, dpi->neighbor);
          dofs_i_want = Nodes[local_node_number]->Nodal_Vars_Info[imtrx]->Num_Unknowns;
          cx[imtrx][index_owner].num_dofs_recv += dofs_i_want;
         }
     }

  /*
   *  Count up the number of degrees of freedom that need to be sent
   *  from this processor to the p'th processor. 
   *
   */
  for (p = 0; p < dpi->num_neighbors; p++) 
     {
      for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++)
         {
           cx[imtrx][p].num_dofs_send = 0;
           for (i = ptr_node_send[p]; i < ptr_node_send[p+1]; i++) 
             {
               local_node_number = list_node_send[i];
               cx[imtrx][p].num_dofs_send +=
                 Nodes[local_node_number]->Nodal_Vars_Info[imtrx]->Num_Unknowns;
             }
         }
     }
  /*
   *  Set up the index array for pointing into list_dof_send, ptr_dof_send.
   *  The pth entry points to the first entry in list_dof_send that
   *  gets sent to the pth processor.
   *
   *    ptr_dof_send[] gets allocated here, and is never freed.
   */
  ptr_dof_send = malloc(sizeof(int *) * upd->Total_Num_Matrices);
  for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
    ptr_dof_send[imtrx] = alloc_int_1(dpi->num_neighbors + 1, 0);
    for (p = 0; p < dpi->num_neighbors; p++) {
      ptr_dof_send[imtrx][p+1]  = ptr_dof_send[imtrx][p]  + cx[imtrx][p].num_dofs_send;
    }
  }

  /*
   *  Set up the index array, list_dof_send, that points to degrees of
   *  freedom that get sent from this processor to neighboring
   *  processors. list_dof_send[i] is the offset from the base of the
   *  solution vector to the i'th dof that needs to be sent to
   *  a neighboring processor.
   *
   *     list_dof_send[] is allocated here, and is never freed.
   */
  list_dof_send = malloc(sizeof(int *) * upd->Total_Num_Matrices);
  for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
    if (ptr_dof_send[imtrx][dpi->num_neighbors] > 0) {
      list_dof_send[imtrx] = alloc_int_1(ptr_dof_send[imtrx][dpi->num_neighbors], -1);
    }
  }

  for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++)
     {
       index = 0;
      for (i = 0; i < ptr_node_send[dpi->num_neighbors]; i++) 
         {
          local_node_number = list_node_send[i];
          base = Nodes[local_node_number]->First_Unknown[imtrx];
          for (p = 0; p < Nodes[local_node_number]->Nodal_Vars_Info[imtrx]->Num_Unknowns; p++) 
             {
              list_dof_send[imtrx][index] = base + p;
              index++;
             }
         }
     }
}
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/


static void 
build_node_recv_indeces(Exo_DB *exo, Dpi *dpi)

    /*****************************************************************
     *
     * build_node_recv_indeces() -- for ea neighbor proc, what this proc wants
     *
     * This is a simple list based on the concept of nodes. It is used to
     * created indeces on a per processor basis that can be used to communicate
     * things like fill equation information.
     *
     * Using the mesh and distributed processing information, this builds the
     * integer vectors
     *
     *		ptr_node_recv
     *
     * that are defined in rf_fem.h that decribe this processor's communication
     * with all of its neighbors. Then, the individual pointers in the Comm_Ex
     * structures are pointed into this list.
     ********************************************************************/
{
  int i, local_node_number, index_owner, owner,  prev_owner, num_owner_changes;
  NODE_INFO_STRUCT *node_ptr;

  ptr_node_recv = alloc_int_1(dpi->num_neighbors+1, 0);

  /*
   * Gather the frequency profile for each processor among the external
   * nodes...
   */
  prev_owner = -1;
  num_owner_changes = 0;
  for (i = 0; i < dpi->num_external_nodes; i++) {
    local_node_number = dpi->num_internal_nodes + dpi->num_boundary_nodes + i;
    node_ptr = Nodes[local_node_number];
    owner = dpi->node_owner[local_node_number];
    node_ptr->Proc = owner;
    index_owner = in_list(owner, 0, dpi->num_neighbors, dpi->neighbor);
    ptr_node_recv[index_owner+1]++;
    
    if (prev_owner != owner) num_owner_changes++;
    prev_owner = owner;
  }

  /*
   *  Put in a check to see if the external nodes on this processor
   *  are numbered so that they are contiguous wrt the owning neighboring
   *  processor. Aztec's communications pattern is contingent on this
   *  requirement.
   *   We do this check by counting the number of ownership changes as
   *  we raster through the external nodes. 
   */
  if (num_owner_changes != dpi->num_neighbors) {
    fprintf(stderr,
	    "Proc %d: num_owner_changes %d not equal to dpi->num_neighbors %d\n",
	    ProcID, num_owner_changes, dpi->num_neighbors);
    EH(GOMA_ERROR, "External node ownership inconsistency!");
  }

  /*
   * Convert frequency into pointer list.
   */
  for (i = 0; i < dpi->num_neighbors; i++) {
    ptr_node_recv[i+1] += ptr_node_recv[i];
  }
  if (ptr_node_recv[dpi->num_neighbors] != dpi->num_external_nodes) {
    EH(GOMA_ERROR, "Mismatch in external node ownership.");
  }
  
  return;
}
/*********************************************************************************/
/*********************************************************************************/
/*********************************************************************************/

void
output_comm_stats(Dpi *dpi, Comm_Ex **cx)

    /*******************************************************************
     *
     * output_comm_stats:
     *
     *   Output some global statistics concerning the dof communications
     *   pattern.
     *******************************************************************/
{
  int max_mesg = 0, min_mesg = 1000000, p;
  int gmax_neighbor, gmax_neighbor_proc;
  double avg_mesg = 0.0, gavg_mesg, gavg_neighbor;
  int gmin_neighbor, gmin_neighbor_proc;
  int gmax_mesg_proc;
  int gmin_mesg_proc;
  int imtrx;
  
  gmax_neighbor = gmaxloc_int(dpi->num_neighbors, ProcID,
			      &gmax_neighbor_proc);
  gmin_neighbor = gminloc_int(dpi->num_neighbors, ProcID,
			      &gmin_neighbor_proc);
  gavg_neighbor = gavg_double((double) dpi->num_neighbors);

  for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
    for (p = 0; p < dpi->num_neighbors; p++) {
      max_mesg = MAX(max_mesg, cx[imtrx][p].num_dofs_send);
      min_mesg = MIN(min_mesg, cx[imtrx][p].num_dofs_send);
      avg_mesg += cx[imtrx][p].num_dofs_send;
    }
  }

  gmaxloc_int(max_mesg, ProcID, &gmax_mesg_proc);
  gminloc_int(min_mesg, ProcID, &gmin_mesg_proc);
  if (dpi->num_neighbors > 0) {
    avg_mesg /= dpi->num_neighbors;
  }
  gavg_mesg = gavg_double(avg_mesg);
 
  if (ProcID == 0) {
    printf("\n\n----------------- Communications Stats ------------\n");
    printf("\n\t\tNumber of Neighbors\n");
    printf("\t\tMax: %5d (Proc = %d)\t\t Avg: %6.3f\t\tMin: %d (Proc = %d)\n",
	   gmax_neighbor, gmax_neighbor_proc, gavg_neighbor,
	   gmin_neighbor, gmin_neighbor_proc);
    printf("\n\t\tSize of Messages\n");
    printf("\t\tMax: %5d (Proc = %d)\t\t Avg: %6.3f\t\tMin: %d (Proc = %d)\n",
	   max_mesg, gmax_mesg_proc, gavg_mesg, min_mesg, gmin_mesg_proc);
    printf("--------------------------------------------------\n");
    fflush(stdout);
  }
}
/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
