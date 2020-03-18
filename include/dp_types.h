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
 * 
 * A data structure is used to contain information about what is to be sent
 * and received from a given neighboring processor. This includes both nodal
 * data (eg, fill equation) and the full degree of freedom (dof) information.
 *
 * Once the MPI indexed datatypes have been committed, they can be used with
 * both the x[] and the R[] vectors since the indexing values should be the
 * same.
 *
 * Created: 1997/08/13 10:40 MDT pasacki@sandia.gov
 * Recreated: 1999/04/27 10:42 MDT pasacki@sandia.gov
 */

/*
 *$Id: dp_types.h,v 5.1 2007-09-18 18:53:41 prschun Exp $
 */

#ifndef GOMA_DP_TYPES_H
#define GOMA_DP_TYPES_H

#include <mpi.h>

/*
 * Used during the initial allocation of the data description structure.
 * This limit is not hard; checking during data member addition is done
 * with reallocation as required so that no limit exists on the number of
 * members in the datatype.
 */

#ifndef SZ_DDD_INIT_MAX_MEMBERS
#define SZ_DDD_INIT_MAX_MEMBERS 100
#endif

/*
 * dalloc() -- C preprocessor macro conditionally allocates space for doubles
 *
 *	Only works if the number of elements is greater than zero.
 *
 *	Used in conjunction with registering double vectors for MPI transport
 *	cross processor (see crdv() below).
 *
 */

#define dalloc(size, ptr) if (size > 0) ptr= alloc_dbl_1(size, 0.0)

/*
 * crdv() -- C preprocessor macro conditionally register double vector 
 *
 *	     This registration is done prior to MPI transport if size > 0.
 *
 * Presumably, any nontrivial length vectors of user defined data were
 * allocated receiving space on other processors.
 *
 * Assume: 
 *		(1) "n" refers to some "DDD"
 *
 *		(2) Data type is MPI_DOUBLE.
 */

#define crdv(len, addr) if((len)>0 && (addr) != NULL) ddd_add_member(n, (addr), (len), MPI_DOUBLE); if ( (len)>0 && (addr) == NULL ) { printf("P_%d: crdv ERROR: member %d NULL address but nonzero length=%d %s line %d!\n", ProcID, n->num_members, (len), __FILE__, __LINE__); fflush(stdout); }




/*
 * All of the necessary pieces to define an MPI derived datatype.
 */

struct Derived_Datatype_Description
{
  int num_members;
  int max_members;
  int *block_count;
  MPI_Datatype *data_type;
  MPI_Aint     *address;
  MPI_Datatype  new_type;

  MPI_Aint extent;              /* extent of new derived data type */
  MPI_Aint lb;
  /*  int count;                        count of new derived data type */
  int size;                     /* size of new derived data type */
};

typedef struct Derived_Datatype_Description *DDD;


struct Communication_Exchange
{
  int neighbor_name;

  int num_nodes_recv;
  int num_nodes_send;
  int *local_nodeces_send;

  /*
   * Once these are registered with MPI_blah(), they can be used for
   * offsets from both the x_fill[] and the R_fill[] vectors! The integer
   * versions are used to verify global node names between processors. The
   * dbl versions are used to communcate actual data.
   */
/*
  MPI_Datatype mpidt_i_node_send;
  MPI_Datatype mpidt_d_node_send;
  MPI_Datatype mpidt_i_node_recv;
  MPI_Datatype mpidt_d_node_recv;
*/
  /* likewise for fill_nodes */

  int num_fill_nodes_recv;
  int num_fill_nodes_send;
  int *local_fill_nodeces_send;
/*
  MPI_Datatype mpidt_i_fill_node_send;
  MPI_Datatype mpidt_d_fill_node_send;
  MPI_Datatype mpidt_i_fill_node_recv;
  MPI_Datatype mpidt_d_fill_node_recv;
*/
  /* yet again for dofs */

  int num_dofs_recv;
  int num_dofs_send;
  int *local_dofdeces_send;

  /*
   * Once these are registered with MPI_blah(), they can be used for
   * offsets from both the x[] and the R[] vectors! The int versions are used
   * to communicate global dof names between processors to verify what is
   * being sent and received. The dbl types are used to communicate the actual
   * values.
   */
/*
  MPI_Datatype mpidt_i_dof_send;
  MPI_Datatype mpidt_d_dof_send;
  MPI_Datatype mpidt_i_dof_recv;
  MPI_Datatype mpidt_d_dof_recv;
*/  
};

typedef struct Communication_Exchange Comm_Ex;

/*
 *  This structure below is needed for the kernal operation of
 *  communications between neighboring processors.
 */

struct Comm_Neighbor_Proc {
    int neighbor_ProcID;         /* Integer containing the ProcID of the neighbor
				    with which to communicate */
    void *send_message_buf;      /* send message buffer (contiguous in space) -
				    This is sent to neighbor_ProcID */
    int send_message_length;     /* Length of the send message buffer in bytes */
    void *recv_message_buf;      /* receive message buffer (contiguous in space)
				    This is received from neighbor_ProcID */
    int recv_message_length;     /* Length of the receive message buffer */
    MPI_Request send_request;      /* Request object for send */
    MPI_Request recv_request;    /* request object for receive */
    MPI_Status  send_status;     /* status object for send operation */
    MPI_Status  recv_status;     /* status object for receive operation */
};
typedef struct Comm_Neighbor_Proc COMM_NP_STRUCT;


#endif /* GOMA_DP_TYPES_H */
