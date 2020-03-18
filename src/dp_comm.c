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
 
#include "dp_comm.h"
#include "dp_map_comm_vec.h"
#include "dp_types.h"
#include "dpi.h"
#include "rf_allo.h"
#include "rf_fem.h"
#ifdef USE_RCSID
static char rcsid[] = "$Id: dp_comm.c,v 5.1 2007-09-18 18:53:41 prschun Exp $";
#endif

/* System Include files */

/* User include files */
/*
#include "std.h"
#include "el_elm.h"

#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "mm_as.h"
#include "rf_masks.h"
#include "rf_bc_const.h"

#include "mm_eh.h"

#include "exo_struct.h"
#include "dpi.h"
#include "dp_types.h"
*/

#define GOMA_DP_COMM_C

/********************************************************************/
/********************************************************************/
/********************************************************************/

void 
exchange_dof(Comm_Ex *cx,  Dpi *dpi,  double *x, int imtrx)

    /************************************************************
     *
     *  exchange_dof():
     *
     *  send/recv appropriate pieces of a dof-based double array
     ************************************************************/
{
  COMM_NP_STRUCT *np_base, *np_ptr;
  double *ptr_send_list, *ptr_recv_list;
  register double *ptrd;
  register int *ptr_int, i;
  int p;
  int num_neighbors = dpi->num_neighbors;
  int total_num_send_unknowns;

  if (num_neighbors == 0) return;

#ifdef PARALLEL
  total_num_send_unknowns = ptr_dof_send[imtrx][dpi->num_neighbors];
  np_base = alloc_struct_1(COMM_NP_STRUCT, dpi->num_neighbors);
  ptrd = (double *) alloc_dbl_1(total_num_send_unknowns, DBL_NOINIT);    
  ptr_send_list = ptrd;

  /*
   * gather up the list of send unknowns
   */
  ptr_int = list_dof_send[imtrx];
  for (i = total_num_send_unknowns; i > 0; i--) {
    *ptrd++ = x[*ptr_int++];
  }

  /*
   * store base address for the start of the external degrees of freedom
   * in this vector
   */
  ptr_recv_list = x + num_internal_dofs[imtrx] + num_boundary_dofs[imtrx];
  
  np_ptr = np_base;
  for (p = 0; p < dpi->num_neighbors; p++) {
    np_ptr->neighbor_ProcID = cx[p].neighbor_name;
    np_ptr->send_message_buf = (void *)
	                       (ptr_send_list + ptr_dof_send[imtrx][p]);
    np_ptr->send_message_length = sizeof(double) * cx[p].num_dofs_send;
    np_ptr->recv_message_buf = (void *) ptr_recv_list;
    np_ptr->recv_message_length = sizeof(double) * cx[p].num_dofs_recv;
    ptr_recv_list += cx[p].num_dofs_recv;
    np_ptr++;
  }
  exchange_neighbor_proc_info(dpi->num_neighbors, np_base);
  safer_free((void **) &np_base);
  safer_free((void **) &ptr_send_list);
#endif /* PARALLEL */
}
/********************************************************************/
/********************************************************************/
/********************************************************************/
/*    
{
#ifdef PARALLEL
  int p;

  for ( p=0; p<d->num_neighbors; p++)
    {
      MPI_Irecv(a,			
		1,			
		cx[p].mpidt_d_dof_recv,	
		cx[p].neighbor_name,	
		555,			
  		MPI_COMM_WORLD,
		Request + Num_Requests*p + 2 );   

      MPI_Isend(a,			
		1,				
		cx[p].mpidt_d_dof_send,	
		cx[p].neighbor_name,
		555,
		MPI_COMM_WORLD,
		( Request + Num_Requests*p + 3 ) );
    }

  for ( p=0; p<d->num_neighbors; p++)
    {
      MPI_Wait( Request + Num_Requests*p + 2, Status  + Num_Requests*p + 2 );
      MPI_Wait( Request + Num_Requests*p + 3, Status  + Num_Requests*p + 3 );
    }
#endif

  return;
}
*/
/********************************************************************/
/********************************************************************/
/********************************************************************/

void 
exchange_node(Comm_Ex *cx,  Dpi *dpi,  double *x)

    /************************************************************
     *
     *  exchange_dof():
     *
     *  send/recv appropriate pieces of a node-based double array
     ************************************************************/
{
  COMM_NP_STRUCT *np_base, *np_ptr;
  double *ptr_send_list, *ptr_recv_list;
  register double *ptrd;
  register int *ptr_int, i;
  int p;
  int num_neighbors = dpi->num_neighbors;
  int total_num_send_unknowns;

  if (num_neighbors == 0) return;

#ifdef PARALLEL
  total_num_send_unknowns = ptr_node_send[dpi->num_neighbors];
  np_base = alloc_struct_1(COMM_NP_STRUCT, dpi->num_neighbors);
  ptrd = (double *) alloc_dbl_1(total_num_send_unknowns, DBL_NOINIT);    
  ptr_send_list = ptrd;

  /*
   * gather up the list of send unknowns
   */
  ptr_int = list_node_send;
  for (i = total_num_send_unknowns; i > 0; i--) {
    *ptrd++ = x[*ptr_int++];
  }

  /*
   * store base address for the start of the entries corresponding
   * to external nodes in this vector
   */
  ptr_recv_list = x + dpi->num_internal_nodes +
                      dpi->num_boundary_nodes;
  
  np_ptr = np_base;
  for (p = 0; p < dpi->num_neighbors; p++) {
    np_ptr->neighbor_ProcID = cx[p].neighbor_name;
    np_ptr->send_message_buf = (void *)
	                       (ptr_send_list + ptr_node_send[p]);
    np_ptr->send_message_length = sizeof(double) * cx[p].num_nodes_send;
    np_ptr->recv_message_buf = (void *) ptr_recv_list;
    np_ptr->recv_message_length = sizeof(double) * cx[p].num_nodes_recv;
    ptr_recv_list += cx[p].num_nodes_recv;
    np_ptr++;
  }
  exchange_neighbor_proc_info(dpi->num_neighbors, np_base);
  safer_free((void **) &np_base);
  safer_free((void **) &ptr_send_list);
#endif
}
/********************************************************************/
/********************************************************************/
/********************************************************************/
