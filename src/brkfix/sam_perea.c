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

/* assess_weights() -- for a goma-like FE problem, assess vertex and edge wts
 *
 * Notes:
 *	[1] The interactions are used to assess weights of interactions. The
 *          total number of terms collected for each equation/dof associated
 *          with a finite element node are used to compute vertex weights.
 *
 *	Assess weights used to build the overall problem graph.
 *
 *      The graph verteces will correspond to nodes of the
 *      finite element mesh. Graph edges correspond to nonzero
 *	interaction between equations for dofs at a given node
 *	and variable dofs at a given node. The interaction is "thru"
 *      an element.
 *
 *	Two kinds of weights will be accumulated in a survey
 *	of the node-node interactions. 
 *
 *	First, the "or" weights represent the logical OR of
 *	all potential interactions between eqn_node & var_node.
 *
 *	These weights represent the potential communication 
 *	cost and be expressed as graph edge weights. The actual
 *      communication costs will be somewhat less, since
 *      variables communicated for one equation may be used
 *      for another. A different examination of the or matrix
 *	can describe the number of nonzero matrix entries in
 *      the overall sparse matrix.
 *
 *	Second, the  "add" weights represent the arithmetic sum
 *	of all the interactions between an eqn_node & var_node.
 *	These weights represent the potential computation cost
 *	and will be expressed as graph vertex weights. This
 *	add matrix can be used to determine how many individual
 *      element level contributions will be made to the
 *      sparse matrix.
 *
 *	Finally, this scheme really requires directed graphs to hold
 *      the information it generates. However, we'll base our cost
 *      estimates for edge cutting on the largest of the two costs.
 *      Therefore, scan through and symmetrize the edge weights so that
 *
 *		edge_weight( i depends j ) = edge_weight ( j depends i )
 *	
 *			= max( original values)
 *
 *
 * Still todo: subparametric special elements, extra costs associated
 *             with doing boundary conditions, better handling of lazy
 *             nodes with absolutely no eqns or vars...
 *
 * Created: 1997/05/09 07:39 MDT pasacki@sandia.gov
 *
 * Revised:
 */

#define _SAM_PEREA_C

#include <config.h>

#include <stdio.h>

#ifdef STDC_HEADERS
#include <stdlib.h>
#endif

#include "map_names.h"
#include "std.h"
#include "eh.h"
#include "aalloc.h"
#include "exo_struct.h"
#include "nodesc.h"
#include "brkfix_types.h"
#include "sam_perea.h"

extern int in_list		/* utils.c */
PROTO((int ,			/* val    - what integer value to seek */
       int *,			/* start  - where to begin looking */
       int ));			/* length - how far to search from start */

extern int fence_post		/* utils.c */
PROTO((int ,			/* val    - integer whose category we seek */
       int *,			/* array  - where to look  */
       int ));			/* length - how far to search in array */



void
assess_weights(Exo_DB *x, 
	       Bevm ***mult,
	       int ***evd,
	       int *ebl,
	       int *np,
	       int *el,
	       int *ep,
	       int *nl,
	       int *node_kind,
	       Node_Description **pnd,
	       int *num_basic_eqnvars,
	       int *eqn_node_names,
	       int *var_node_names,
	       int *nnz_contribute,
	       int *nat_contribute,
	       int *ccs_contribute)
{
  int col;
  int column_max;
  int *common_elements;

  int e;
  int eb_index;
  int elem;
  int eqn_node;
  int evid;
  int ewt;

  int i;
  int ieb;
  int ii;
  int inn;
  int iv;

  int j;
  int jj;

  int l;

  int m;
  int map_e_index[MAX_EQNVARS];
  int map_v_index[MAX_EQNVARS];

  int n;
  int nbev;
  /*  int ne;*/
  int neb;
  int *neighbor_nodes;
  int nelems_this_node;
  int nn;
  /*  int nns;*/
  int node;
  /*  int nss;*/

  int *nn_eids;
  int *nn_ewts;
  int *nn_vids;
  int *nn_vwts;

  int **nn_or;
  int **nn_add;
  
  int num_connecting_elems;
  int num_neighbor_nodes;
  int num_assembled_terms;
  int num_comm_chunks;
  int num_matrix_nonzeroes;

  int max_common_elements;

  int row;

  int var_node;
  int vwt;

  int where;

  char  err_msg[MAX_CHAR_ERR_MSG];

  Spfrtn sr=0;

  Node_Description *end;
  Node_Description *vnd;

  /*
   * Convenience variables...
   */

  /*
   * ne  = x->num_elems;
   */

  nn  = x->num_nodes;
  neb = x->num_elem_blocks;

  /*
   *  nns = x->num_node_sets;
   *  nss = x->num_side_sets;
   */

  /*
   * Allocate space for the little matrices used for each node-node
   * interaction. For a given eqn-node to var-node interaction, we'll
   * need to have the following information compiled:
   *
   *	the columns of the interaction matrix need to be identified with
   *    the variable ID's at the var-node, as well as the corresponding
   *    weights for each of these variables. The weights are just
   *    products of all of the three kinds of multiplicities for an eqn:
   *    (i) vector/tensor/etc., (ii) concentration, (iii) nodal dof
   *
   *	the rows of the interaction matrix need to be identified with the
   *    equation ID's at the eqn-node, as well as their weights.
   *
   *	Finally, the little AND and OR matrix entries may be computed
   *    by considering the equation-variable dependency matrix for
   *    each of the elements through which the eqn-node depends upon the
   *    var-node.
   */

  nn_eids = (int *) smalloc(MAX_EQNVARS * SZ_INT);
  nn_ewts = (int *) smalloc(MAX_EQNVARS * SZ_INT);

  nn_vids = (int *) smalloc(MAX_EQNVARS * SZ_INT);
  nn_vwts = (int *) smalloc(MAX_EQNVARS * SZ_INT);

  nn_or    = (int **) smalloc(MAX_EQNVARS * sizeof(int *));
  nn_or[0] = (int *) smalloc(MAX_EQNVARS * MAX_EQNVARS * sizeof(int));
  for ( i=1; i<MAX_EQNVARS; i++)
    {
      nn_or[i] = nn_or[i-1] + MAX_EQNVARS;
    }

  nn_add    = (int **) smalloc(MAX_EQNVARS * sizeof(int *));
  nn_add[0] = (int *) smalloc(MAX_EQNVARS * MAX_EQNVARS * sizeof(int));
  for ( i=1; i<MAX_EQNVARS; i++)
    {
      nn_add[i] = nn_add[i-1] + MAX_EQNVARS;
    }


  /*
   * The maximum number of elements connecting any given eqn node with
   * any given var node will likely be the maximum number of elements
   * to which any *single* node belongs. (The worse case scenario will
   * be if the eqn node and the var node are the same.)
   */

  max_common_elements = -1;
  
  for ( n=0; n<nn; n++)
    {
      nelems_this_node = np[n+1] - np[n];

      if ( nelems_this_node > max_common_elements )
	{
	  max_common_elements = nelems_this_node;
	}
    }

  neighbor_nodes = (int *) smalloc(MAX_NEIGHBOR_NODES* SZ_INT);

#ifdef DEBUG
  fprintf(stderr, "max common elements = %d\n", max_common_elements);
#endif

  common_elements = (int *) smalloc(max_common_elements * SZ_INT);

  inn = 0;			/* 0 < inn < length_node_node */

  for ( eqn_node=0; eqn_node<nn; eqn_node++)
    {
      /*
      if ( inn >= length_node_node )
	{
	  EH(-1, "Abnormal lack of space!");
	}
	*/
      /*
       * 0. First, figure out how many *other* nodes interact with this one.
       *
       * 1. Determine interaction based solely on the mesh topology. Every
       *    node connected to the given node by an element is a candidate
       *    for interaction.
       *
       * 2. Find out how many equations(dofs) associated at this given node
       *    there are.
       *
       * 3. Find out how many variables(dofs) are associate with each of the 
       *    interacting nodes. Don't forget that a node interacts with itself.
       *    We'll remove the self loop so it doesn't become a graph edge, but
       *    properly account for the computational load of computing the
       *    self-interaction terms.
       *
       * 4. Cross check the potential interactions using the appropriate
       *    Equation Variable Dependency description that was read in for
       *    each element block into the "evd" array.
       *
       * 5. Do not forget that two nodes may interact through more
       *    than one element and that, therefore, the interaction matrices
       *    can differ betweent the eqns and vars of the two nodes, or the
       *    node and itself. The interactions need to be "or"-ed together
       *    and added together.
       */

      /*
       * Initialize the equation info. It will stay the same for
       * every var_node of the interaction.
       */

      for ( i=0; i<MAX_EQNVARS; i++)
	{
	  nn_eids[i] = UNDEFINED_EQNVARID;
	  nn_ewts[i] = 0;
	}

      end = pnd[node_kind[eqn_node]];

      for ( i=0; i<end->num_basic_eqnvars; i++)
	{
	  nn_eids[i] = end->eqnvar_ids[i];
	  nn_ewts[i] = ( end->eqnvar_wts[i][0] * end->eqnvar_wts[i][1] *
			 end->eqnvar_wts[i][2] );
	}

#ifdef DEBUG
      fprintf(stderr, "enode = (%d) has active eids = ", eqn_node+1);
      for ( i=0; i<end->num_basic_eqnvars; i++)
	{
	  fprintf(stderr, "%d (x%d) ", nn_eids[i], nn_ewts[i]);
	}
      fprintf(stderr, "\n\n");
#endif

      /*
       * Initialize, then load up names of every connecting var_node to 
       * this eqn_node.
       */

      num_neighbor_nodes = 0;

      for ( i=0; i<MAX_NEIGHBOR_NODES; i++)
	{
	  neighbor_nodes[i] = -1;
	}

      /*
       * Examine every element containing this eqn node.
       */

      for ( l=np[eqn_node]; l<np[eqn_node+1]; l++)
	{
	  elem = el[l];
	  
	  /*
	   * This element may belong to an element block where the
	   * interaction is much weaker than might be supposed from 
	   * looking at a large impressive list of active equations at
	   * the eqn_node or the large impressive list of active variables
	   * at the var_node. Fortunately, insofar as this connecting element
	   * is concerned, we need only check the interactions for the
	   * element block concerned.
	   */

	  /*
	   * Look at every node that this element contains.
	   */

	  for ( m=ep[elem]; m<ep[elem+1]; m++ )
	    {
	      var_node = nl[m];
	      
	      /*
	       * If this variable node is not already in the list of 
	       * neighbor nodes, add it.
	       */

	      BULL(var_node, neighbor_nodes, num_neighbor_nodes);
	      
	      if ( num_neighbor_nodes > MAX_NEIGHBOR_NODES )
		{
		  sr = sprintf(err_msg, "@ n = %d too many neighbors.", 
			       eqn_node);
		  EH(-1, err_msg);
		}
	      
	    }
	}
#ifdef DEBUG
      fprintf(stderr, "node (%d) has %d neighbor_nodes\n",
	      eqn_node+1, num_neighbor_nodes);
      fprintf(stderr, "and they are:");
      for ( i=0; i<num_neighbor_nodes; i++)
	{
	  fprintf(stderr, " %d", neighbor_nodes[i]+1);
	}
      fprintf(stderr, "\n");
#endif

      for ( iv=0; iv<num_neighbor_nodes; iv++)
	{
	  var_node = neighbor_nodes[iv];

	  /*
	   * Gather together a list of all elements that connect
	   * this var_node with the eqn_node. First, initialize to
	   * the empty list.
	   */

	  num_connecting_elems = 0;

	  for ( i=0; i<max_common_elements; i++)
	    {
	      common_elements[i] = -1;
	    }
	  
	  /*
	   * Look through all the elements this eqn_node touches.
	   * If any of those elements contain var_node, then this
	   * element is common.
	   */

	  for ( l=np[eqn_node]; l<np[eqn_node+1]; l++)
	    {
	      elem = el[l];
	      
	      for ( m=ep[elem]; m<ep[elem+1]; m++)
		{
		  node = nl[m];
		  
		  if ( node == var_node )
		    {
		      /*
		       * Assume this will only happen once per element.
		       */
		      common_elements[num_connecting_elems] = elem;
		      num_connecting_elems++;
		    }
		}
	    }

#ifdef DEBUG
	  fprintf(stderr, 
		  "en=(%d),vn=(%d) thru %d elems: ",
		  eqn_node+1, var_node+1, num_connecting_elems);
	  for ( i=0; i<num_connecting_elems; i++)
	    {
	      fprintf(stderr, "%d ", common_elements[i]+1);
	    }
	  fprintf(stderr, "\n");
#endif

	  /*
	   * At the var_node, find eqnvar_IDs and weights.
	   */
	  
	  for ( i=0; i<MAX_EQNVARS; i++)
	    {
	      nn_vids[i] = UNDEFINED_EQNVARID;
	      nn_vwts[i] = 0;
	    }

	  vnd = pnd[node_kind[var_node]];

	  for ( i=0; i<vnd->num_basic_eqnvars; i++)
	    {
	      nn_vids[i] = vnd->eqnvar_ids[i];
	      nn_vwts[i] = ( vnd->eqnvar_wts[i][0] * vnd->eqnvar_wts[i][1] *
			     vnd->eqnvar_wts[i][2] );
	    }

#ifdef DEBUG
	  fprintf(stderr, "vnode = (%d) has active vids = ", var_node+1);
	  for ( i=0; i<vnd->num_basic_eqnvars; i++)
	    {
	      fprintf(stderr, "%d (x%d) ", nn_vids[i], nn_vwts[i]);
	    }
	  fprintf(stderr, "\n\n");
#endif

	  /*
	   * Now, construct the ADD and OR matrices for this particular
	   * eqn_node -> var_node interaction using the Equation Variable
	   * Dependencies for each element participating in the interaction.
	   *
	   * Initially, empty the matrix.
	   */

	  for ( row=0; row<MAX_EQNVARS; row++)
	    {
	      for ( col=0; col<MAX_EQNVARS; col++)
		{
		  nn_or[row][col]  = 0;
		  nn_add[row][col] = 0;
		}
	    }

	  for ( e=0; e<num_connecting_elems; e++)
	    {
	      elem     = common_elements[e];
	      eb_index = fence_post(elem, ebl, neb+1);

	      nbev     = num_basic_eqnvars[eb_index];

	      /*
	       * For this element block index, what are the indeces
	       * corresponding to the active eqns and vars?
	       *
	       * There might be fewer eqnvars known in this eb_index than
	       * are known at either node. 
	       *
	       * Construct a map from the ev_indeces in the
	       * elem block into the ev_indeces of the nodes weight matrices.
	       */
	      
	      /*	      for ( ieb=0; ieb<end->num_basic_eqnvars; ieb++)*/

	      for ( ieb=0; ieb<nbev; ieb++)
		{
		  evid   = mult[eb_index][ieb]->eqnvar_id;

		  where  = in_list(evid, nn_eids, end->num_basic_eqnvars);

		  /*
		   * If this evid is not at this particular node, OK.
		   * Just make the map = -1 and check later
		   */

		  map_e_index[ieb] = where;

		  where  = in_list(evid, nn_vids, vnd->num_basic_eqnvars);

		  map_v_index[ieb] = where;
		}

	      /*
	       * Combine this elements idea of interaction strength
	       * into the OR and ADD matrices built to express the
	       * eqn_node - var_node interaction through all connecting
	       * elements.
	       */

	      for ( i=0; i<nbev; i++)
		{
		  ii = map_e_index[i];
		  if ( ii != -1 )
		    {
		      ewt = nn_ewts[ii];
		      for ( j=0; j<nbev; j++)
			{
			  jj  = map_v_index[j];
			  if ( jj != -1 )
			    {
			      vwt = nn_vwts[jj];
			      if ( evd[eb_index][i][j] != 0 )
				{
				  nn_or[ii][jj]  = MAX( nn_or[ii][jj], vwt);
				  nn_add[ii][jj] += ewt * vwt;
				}
			    }
			}
		    }
		}

	    }

#ifdef DEBUG	  
	  /*
	   * Dump the interaction matrices for this eqn_node-var_node
	   * interaction.
	   */
	  fprintf(stderr, "(%d)-(%d) ", eqn_node+1, var_node+1);

	  fprintf(stderr, "eqn_ids @ (%d): ", eqn_node+1);
	  for ( i=0; i<end->num_basic_eqnvars; i++)
	    {
	      fprintf(stderr, "%d ", nn_eids[i]);
	    }

	  fprintf(stderr, "var_ids @ (%d): ", var_node+1);
	  for ( i=0; i<vnd->num_basic_eqnvars; i++)
	    {
	      fprintf(stderr, "%d ", nn_vids[i]);
	    }
	  fprintf(stderr, "\n");	  

	  fprintf(stderr, "OR\n");
	  for ( i=0; i<end->num_basic_eqnvars; i++)
	    {
	      for ( j=0; j<vnd->num_basic_eqnvars; j++)
		{
		  fprintf(stderr, "%8d", nn_or[i][j]);
		}
	      fprintf(stderr, "\n");
	    }

	  fprintf(stderr, "ADD\n");
	  for ( i=0; i<end->num_basic_eqnvars; i++)
	    {
	      for ( j=0; j<vnd->num_basic_eqnvars; j++)
		{
		  fprintf(stderr, "%8d", nn_add[i][j]);
		}
	      fprintf(stderr, "\n");
	    }

	  fprintf(stderr, "\n");
#endif

	  /*
	   * All dependencies have been accumulated. Distill into scalars
	   * for the eqn_node/var_node interaction strength.
	   */

	  num_matrix_nonzeroes = 0;
	  num_assembled_terms  = 0;

	  for ( i=0; i<end->num_basic_eqnvars; i++)
	    {
	      for ( j=0; j<vnd->num_basic_eqnvars; j++)
		{
		  num_assembled_terms  += nn_add[i][j];
		  num_matrix_nonzeroes += nn_or[i][j];
		}
	    }


	  num_comm_chunks      = 0;	  
	  for ( j=0; j<vnd->num_basic_eqnvars; j++)
	    {
	      column_max = 0;
	      for ( i=0; i<end->num_basic_eqnvars; i++)
		{
		  if ( nn_or[i][j] > column_max )
		    {
		      column_max = nn_or[i][j];
		    }
		}

	      num_comm_chunks += column_max;
	    }

	  eqn_node_names[inn] = eqn_node;
	  var_node_names[inn] = var_node;

	  nnz_contribute[inn] = num_matrix_nonzeroes;
	  nat_contribute[inn] = num_assembled_terms;
	  ccs_contribute[inn] = num_comm_chunks;

	  inn++;
	  
	}
    }

  free(nn_eids);
  free(nn_ewts);

  free(nn_vids);
  free(nn_vwts);

  free(nn_or[0]);
  free(nn_add[0]);

  free(nn_or);
  free(nn_add);

  free(common_elements);
  free(neighbor_nodes);

  if ( sr < 0 ) exit(2);

  return;
}
