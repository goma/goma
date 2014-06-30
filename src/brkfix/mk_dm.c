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

/* make_goma_dofmap.c -- make goma dofmap
 *
 * Given some rudimentary descriptions of the problem and the mesh, concoct the
 *
 * Notes: 
 *		[1] Make up a goma-like degree of freedom map.
 *
 *		[2] Every dof must be associated with a node.
 *
 *		[3] Nodes are traversed in global order.
 *
 *		[4] Degrees of freedom begin with the variable with the
 *	            lowest integer identifier and move upward.
 *
 *		[5] Variables with multiple degrees of freedom at a
 *                  node are lumped together.
 *
 *		[6] Since there are likely to be only a handful of different
 *		    kinds of nodes, we employ a catalog of node types. Each
 *		    node type has different numbers and kinds of variables
 *		    associated with it.
 * 
 * Created: 1997/05/08 10:48 MDT pasacki@sandia.gov
 * 
 * Revised: 
 */

#define _MK_DM_C

#include <config.h>

#include <stdio.h>

#ifdef STDC_HEADERS
#include <stdlib.h>
#endif

#include "map_names.h"
#include "std.h"
#include "aalloc.h"
#include "eh.h"
#include "exo_struct.h"
#include "brkfix_types.h"
#include "nodesc.h"
#include "mk_dm.h"

/*
 * Function prototypes for functions defined in other files...
 */

extern int in_list		/* utils.c */
PROTO((int ,			/* val    - what integer value to seek */
       int *,			/* start  - where to begin looking */
       int ));			/* length - how far to search from start */


extern int fence_post		/* utils.c */
PROTO((int ,			/* val    - integer whose category we seek */
       int *,			/* array  - where to look  */
       int ));			/* length - how far to search in array */



/*
 * Function prototypes for functions defined in this file...
 */

void 
make_goma_dofmap(Exo_DB *x,
		 Bevm ***mult,
		 int ***evd,
		 int ***Lucky,
		 int *num_basic_eqnvars,
		 int *node_kind, 
		 int *node_dof0, 
		 Node_Description **pnd,
		 int *nkn)
{
  int current_dof;
  int num_kinds_nodes;
  int d;
  int e;
  int eb_index;
  int *ev_ids;
  int evid;
  int **ev_wts;
  int found_matching_node_type;

  int i;
  int index_min;
  int j;
  int l;
  int loc;
  int look;
  int m;
  int m2;
  int matching_node_type_index = 0;
  int min;
  int n;
  /*  int ne;*/
  int neb;
  int nn;
  /*  int nns;*/
  /*  int nss;*/
  int num_abevs;		/* number of active basic eqnvars this node */

  int where;

  int *ep  = x->elem_node_pntr;
  int *nl  = x->elem_node_list;
  int *ebl = x->eb_ptr;
  int *np  = x->node_elem_pntr;
  int *el  = x->node_elem_list;

  char  err_msg[MAX_CHAR_ERR_MSG];

  Node_Description *tnd;		/* temporary pointer */

  Spfrtn sr=0;

  /*
   * Convenience variables...
   */

  /*
   * ne  = x->num_elems;
   */

  nn  = x->num_nodes;
  neb = x->num_elem_blocks;

  /*
   * nns = x->num_node_sets;
   * nss = x->num_side_sets;
   */

  num_kinds_nodes = 0;

  /*
   * Round One; Determine how many total degrees of freedom there are in the
   *            problem. For each node, examine every surrounding element,
   *            determining how many degrees of freedom contribute from
   *            that element's corresponding local node number.
   */


  /*
   * Local information for the node under scrutiny. How many distinct
   * basic eqn/vars are active? What are their names? For each active
   * variable, what are the weights associated with the multipliers?
   */


  ev_ids = (int *) smalloc(MAX_EQNVARS * SZ_INT);

  INIT_IVEC(ev_ids, MAX_EQNVARS, UNDEFINED_EQNVARID);

  ev_wts    = (int **) smalloc(MAX_EQNVARS * sizeof(int *));
  ev_wts[0] = (int *) smalloc(MAX_EQNVARS * 3 * sizeof(int));

  for ( i=1; i<MAX_EQNVARS; i++)
    {
      ev_wts[i] = ev_wts[i-1] + 3;
    }

  for ( i=0; i<MAX_EQNVARS; i++)
    {
      for ( j=0; j<3; j++)
	{
	  ev_wts[i][j] = 0;
	}
    }

  /*
   * Main loop over each node in the entire mesh.
   */

  current_dof  = 0;

  node_dof0[0] = 0;		/* beginning of first nodes dof list */

  for ( n=0; n<nn; n++)
    {
#ifdef DEBUG
      fprintf(stderr, "\n\nglobal node [%d]\n\n", n);
#endif
  
      /*
       * Clean out lists of valid variables and their dof contributions
       * at the current node. These little arrays are catalogs of the
       * active variables at the current node that list their IDs as well as
       * their 3 different kinds of weights.
       */
      for ( i=0; i<MAX_EQNVARS; i++)
	{
	  ev_ids[i] = UNDEFINED_EQNVARID;
	  for ( j=0; j<3; j++)
	    {
	      ev_wts[i][j] = 0;
	    }
	}

      /*
       * Initialize counter of number of active basic eqnvars at this node.
       */

      num_abevs = 0;
      
      /*
       * Loop over each element that contains this node.
       */

      for (l=np[n]; l<np[n+1]; l++)
	{
	  e = el[l];

	  /*
	   * Which *local* node number in the element corresponds to the node
	   * whose global name is "n"?
	   */

	  loc = 0;

	  while ( nl[ep[e]+loc] != n && (ep[e]+loc) < ep[e+1] )
	    {
	      loc++;
	    }
#ifdef DEBUG
	  fprintf(stderr, "node (%d) is aka elem (%d), local node (%d)\n",
		  n+1, e+1, loc+1);
#endif
	  if ( ep[e]+loc >= ep[e+1] )
	    {
	      sr = sprintf(err_msg, 
			   "Difficulty finding local node in n=%d, e=%d",
			   n, e);
	      EH(-1, err_msg);
	    }

	  eb_index = fence_post(e, ebl, neb+1);

	  /*
	   * For this element, in this element block, at this local node,
	   * what kinds of variables contribute how much to the dof load?
	   */

	  /*
	   * Look through catalog variables. If the variables at this
	   * local node in this element are already in the catalog, fine.
	   * If they are not in the catalog, add them.
	   */

	  for ( j=0; j<num_basic_eqnvars[eb_index]; j++)
	    {
	      if ( Lucky[eb_index][loc][j] != 0 )
		{
		  /*
		   * Then, in this element block at this local node, the
		   * eqnvar with index j is active!
		   */
		  
		  evid = mult[eb_index][j]->eqnvar_id;
		  where = in_list(evid, ev_ids, MAX_EQNVARS);
		}
	      else
		{
		  /*
		   * The eqnvar with index j is not active at this local
		   * node.
		   */
		  evid = UNDEFINED_EQNVARID;
		  where = -1;
		}

	      if ( where == -1 && evid != UNDEFINED_EQNVARID )
		{
		  ev_ids[num_abevs]    = evid;
		  ev_wts[num_abevs][0] = mult[eb_index][j]->vect_mult;
		  ev_wts[num_abevs][1] = mult[eb_index][j]->conc_mult;
		  ev_wts[num_abevs][2] = mult[eb_index][j]->ndof_mult;
		  num_abevs++;
#ifdef DEBUG
		  fprintf(stderr, 
			  "\tevid not in this nodes catalog. adding...\n");
		  fprintf(stderr, 
			  "\t\teqnvar ID, wts = %d ( %d %d %d )\n", evid,
			  mult[eb_index][j]->vect_mult,
			  mult[eb_index][j]->conc_mult,
			  mult[eb_index][j]->ndof_mult);
#endif
		}
	    }
	}

      /*
       * Now that each element has been researched for the appropriate
       * contribution, re-order the ev_ids so that they appear in ascending
       * numerical order (sort weights, too.)
       */
      if ( num_abevs > MAX_EQNVARS )
	{
	  EH(-1, "Too many active eqnvars this node");
	}
      
      for ( m=0; m<num_abevs; m++)
	{
	  index_min = m;
	  min       = ev_ids[m];
	  for ( m2=m; m2<num_abevs; m2++)
	    {
	      if ( ev_ids[m2] < min )
		{
		  min       = ev_ids[m2];
		  index_min = m2;
		}
	    }
	  if ( index_min != m )
	    {
	      ISWAP(ev_ids[m], ev_ids[index_min]);
	      for ( j=0; j<3; j++ )
		{
		  ISWAP(ev_wts[m][j], ev_wts[index_min][j]);
		}
	    }
	}
#ifdef DEBUG
      fprintf(stderr, "Sorted list of basic eqnvars at node (%d)\n", n+1);
      for ( i=0; i<num_abevs; i++)
	{
	  fprintf(stderr, "eqnvar[%d] ID=%d wts = %d %d %d\n",
		  i, ev_ids[i], ev_wts[i][0], ev_wts[i][1], ev_wts[i][2]);
	}
#endif
      /*
       * A complete sorted list of the active eqnvars for this node has been
       * obtained. Now,
       *
       * Look through the catalog of existing Prototype Node Descriptions
       * to see if this one is like an earlier one. If it is, then good.
       * If it's not, then make up a new Prototype Node Description like
       * this one!
       */

      found_matching_node_type = FALSE;
      d = 0;

#ifdef DEBUG
      fprintf(stderr, "\nCurrently, num_kinds_nodes = %d\n", 
	      num_kinds_nodes);
#endif
      
      while ( !found_matching_node_type && d<num_kinds_nodes )
	{
#ifdef DEBUG
	  fprintf(stderr, "while checking for equiv old nodedescs\n");
#endif

	  if ( pnd[d]->num_basic_eqnvars == num_abevs )	/* maybe... */
	    {
	      look=0;
	      while ( ( pnd[d]->eqnvar_ids[look] == ev_ids[look] ) &&
		      ( look < num_abevs ) )
		{
		  look++;
		}
	      if ( look == num_abevs )
		{ 
		  found_matching_node_type = TRUE;
		  matching_node_type_index = d;
		}
	    }
	  d++;
	}

      if ( found_matching_node_type )
	{
	  /*
	   * Great, use an old node description for *this* node description!
	   */
	  node_kind[n] = matching_node_type_index;
	  tnd = pnd[matching_node_type_index];		

#ifdef DEBUG
	  fprintf(stderr, "Hey, old description %d fits this node!\n",
		  matching_node_type_index);
#endif

	}
      else
	{
	  /*
	   * Hmmm. Need to allocate and setup a new node description.
	   */

	  node_kind[n] = num_kinds_nodes;

#ifdef DEBUG
	  fprintf(stderr, "\ncreating new node description %d\n",
		  num_kinds_nodes);
#endif
	  tnd = pnd[num_kinds_nodes];

	  /*
	   * Initialize structure before tailoring it to describe
	   * this new node type.
	   */

	  tnd->num_basic_eqnvars = -1;
	  for ( i=0; i<MAX_EQNVARS; i++)
	    {
	      tnd->eqnvar_ids[i]     = UNDEFINED_EQNVARID;
	      tnd->eqnvar_wts[i][0]  = 0;
	      tnd->eqnvar_wts[i][1]  = 0;
	      tnd->eqnvar_wts[i][2]  = 0;
	    }


	  tnd->num_basic_eqnvars = num_abevs;

	  for ( i=0; i<num_abevs; i++ )
	    {
	      tnd->eqnvar_ids[i] = ev_ids[i];
	      for ( j=0; j<3; j++ )
		{
		  tnd->eqnvar_wts[i][j] = ev_wts[i][j];
		}
	    }

	  num_kinds_nodes++;	/* there is a new one now! */

#ifdef DEBUG
	  fprintf(stderr, "\n\tNew node description:\n" );
	  fprintf(stderr, "\t\ttnd->num_basic_eqnvars = %d\n",
		  tnd->num_basic_eqnvars);
	  fprintf(stderr, "\t\ttnd->eqnvar_ids[] = ");
	  for ( i=0; i<num_abevs; i++)
	    {
	      fprintf(stderr, "%5d", tnd->eqnvar_ids[i]);
	    }
	  fprintf(stderr, "\n");
	  fprintf(stderr, "\t\ttnd->eqnvar_wts[][] = ");
	  for ( i=0; i<num_abevs; i++)
	    {
	      fprintf(stderr, "\t\t\t%d wt=(", tnd->eqnvar_ids[i]);
	      for ( j=0; j<3; j++ )
		{
		  fprintf(stderr, "%d", tnd->eqnvar_wts[i][j]);
		}
	      fprintf(stderr, ")\n");
	    }
#endif
	  if ( num_kinds_nodes >= MAX_NODE_KINDS )
	    {
	      sr = sprintf(err_msg, 
			   "@ node (%d) num_kinds_nodes >= MAX_NODE_KINDS", 
			   n+1);
	      EH(-1, err_msg);
	    }
	}

      /*
       * With this node's repertoire of basic eqnvars established, we can
       * determine the initial dof of the next node.
       */

      
      for ( i=0; i<tnd->num_basic_eqnvars; i++)
	{
	  current_dof += ( tnd->eqnvar_wts[i][0] * 
			   tnd->eqnvar_wts[i][1] *
			   tnd->eqnvar_wts[i][2] );
	}


      node_dof0[n+1] = current_dof;
    } 

#ifdef DEBUG
  fprintf(stderr, "Number of node kinds = %d\n", num_kinds_nodes);
  fprintf(stderr, "Dump of node_dof0 pointers:\n");
  for ( i=0; i<nn; i++)
    {
      fprintf(stderr, "node (%d) dofs: %d <= dof < %d\n", i+1, node_dof0[i], 
	      node_dof0[i+1]);
    }

  for ( i=0; i<nn; i++)
    {
      fprintf(stderr, "node (%d) is kind %d\n", i+1, node_kind[i]);
    }

#endif
  *nkn = num_kinds_nodes;

  free(ev_ids);
  free(ev_wts[0]);
  free(ev_wts);
  
  if ( sr < 0 ) exit(2);

  return;
}
