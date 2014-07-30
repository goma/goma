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

/* emuck - element level mucking around
 *
 * $Header: /projects/seataf/CVS/ACCESS/analysis/goma/brkfix/emuck.c,v 1.2 2007-04-01 20:23:32 tabaer Exp $
 */

#define _EMUCK_C

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>

#include "map_names.h"
#include "std.h"
#include "eh.h"
#include "aalloc.h"
#include "exo_struct.h"
#include "exo_conn.h"
#include "utils.h"
#include "dpi.h"
#include "emuck.h"

/* assign_elem_ownership() - assign set proc based on nodal decomp
 *
 * This will result in the assignment of each and every element in the
 * monolith to exactly one of the set/procs that is available.
 *
 * Rigorously this might be part of the graph partition, but for expediency,
 * we'll distribute elements according to which processor owns the centroid
 * node, assuming there is a centroid node.
 *
 * If there is not a centroid node, then assign according to which processor
 * owns the most nodes in the element. Finally, break ties so as to equi-
 * distribute the number of owned elements between processors.
 *
 * Notes:
 *       1. The int ** arrays get allocated here. You're responsible for
 *          eventually freeing them.
 */


void
assign_elem_ownership
(const Exo_DB *exo,		/* monolith FE db w/ connectivity       (in) */
 const int num_sets,		/* how many pieces                      (in) */
 const int *nodal_own_pntr,	/* pointers into list                   (in) */
 const int *nodal_own_list,	/* procsets "owning" nodes              (in) */
 int **pntr_element_owner,	/* which proc owns this elem           (out) */
 int **pntr_element_owner_dist, /* num elems owned by each proc        (out) */
 int **pntr_element_bnd_stat)	/* boolean, !boundary=internal         (out) */
{
  int e;
  int err;

  int f;

  int i;

  int *element_owner;
  int *element_owner_dist;
  int *element_bnd_stat;

  int max_dex;
  int max_val;

  int min_dex;
  int min_val;

  int n;
  int neighbor_elem;
  int node;
  int node_owner=-1;
  int num_owners;
  int num_unassigned_elems;
  int num_viable_owners;

  int previous_node_owner;

  int this_owner;

  int *viable_owners;

  if ( num_sets < 2 )
    {
      fprintf(stderr, "Number of sets = %d\n", num_sets);
    }

  /*
   * Allocate space to show who owns each piece - initialize.
   */
  
  element_owner            = (int *) smalloc(exo->num_elems*sizeof(int));
  *pntr_element_owner      = element_owner;


  element_owner_dist       = (int *) smalloc(num_sets*sizeof(int));
  *pntr_element_owner_dist = element_owner_dist;


  element_bnd_stat         = (int *) smalloc(exo->num_elems*sizeof(int));
  *pntr_element_bnd_stat   = element_bnd_stat;

  /*
   * Initialize: the owner is "-1" (nobody)...
   */

  for ( e=0; e<exo->num_elems; e++)
    {
      element_owner[e] = -1;
    }

  /*
   * Initialize: each set/proc now owns exactly 0 elements...
   */

  for ( i=0; i<num_sets; i++)
    {
      element_owner_dist[i] = 0;
    }

  /*
   * Phase I - elements with every node owned by the same set/proc
   *           are assigned to that same set/proc. [i.e., private elems
   *           have clearly assigned ownership.]
   */

  num_unassigned_elems = exo->num_elems;

#ifdef DEBUG  
  D_P(info, 1, "Starting Phase I element assignment\n");

  fprintf(stderr, "Starting Phase I (homogeneous node owner) attempt.\n");
#endif

  for ( e=0; e<exo->num_elems; e++)
    {
      num_owners=0;
      previous_node_owner=-1;
      for ( n=exo->elem_node_pntr[e]; n<exo->elem_node_pntr[e+1]; n++)
	{
	  node       = exo->elem_node_list[n];
	  node_owner = nodal_own_list[nodal_own_pntr[node]];
	  if ( node_owner != previous_node_owner )
	    {
	      num_owners++;
	    }
	  previous_node_owner = node_owner;
	}
      if ( num_owners == 1 )
	{
	  element_owner[e] = node_owner;
	  element_owner_dist[node_owner]++;
	  num_unassigned_elems--;
#ifdef DEBUG  
	  fprintf(stderr, "Element [%d] has nodes owned by 1 owner (%d).\n", 
		  e, node_owner);
#endif
	}
    }

  /*
   * Phase II - for each of the unassigned elements, look for a centroid
   *            node. If there is one, then it determines who owns this
   *            element.
   *
   *            Centroid nodes may be distinguished by the fact that 
   *            "there can be only ONE element" in their node-elem 
   *            connectivity. These are the Highlander nodes...
   */

#ifdef DEBUG  
  fprintf(stderr, "Phase I assigned %d/%d elements (%g%%).\n", 
	  exo->num_elems - num_unassigned_elems, exo->num_elems,
	  1e2*(double)(exo->num_elems-num_unassigned_elems)/
	  (double)exo->num_elems);

  for ( i=0; i<exo->num_elems; i++)
    {
      fprintf(stderr, "Elem [%d] owned by set/proc %d\n", i, element_owner[i]);
    }

  for ( i=0; i<num_sets; i++)
    {
      fprintf(stderr, "Number of elements owned by set/proc %d = %d\n",
	      i, element_owner_dist[i]);
    }

  fprintf(stderr, "Starting Phase II (centroid owner) attempt.\n");
#endif

  if ( num_unassigned_elems > 0 )
    {
      for ( e=0; e<exo->num_elems; e++)
	{
	  if ( element_owner[e] == -1 )
	    {
	      if ( ( exo->elem_node_pntr[e+1] - exo->elem_node_pntr[e] ) == 1 )
		{
		  node       = exo->elem_node_list[exo->elem_node_pntr[e]];
		  node_owner = nodal_own_list[nodal_own_pntr[node]];
		  element_owner[e]   = node_owner;
		  element_owner_dist[node_owner]++;
		  num_unassigned_elems--;
		}
	    }
	}
    }

#ifdef DEBUG  
  fprintf(stderr, "Phase I, II assigned %d/%d elements (%g%%).\n", 
	  exo->num_elems - num_unassigned_elems, exo->num_elems,
	  1e2*(double)(exo->num_elems-num_unassigned_elems)/
	  (double)exo->num_elems);
  fprintf(stderr, "Starting Phase III (Robin Hood) attempt.\n");
#endif

  /*
   * Phase III - For any elements that are *still* unassigned, look
   *             at who owns the nodes in the element. Whichever of
   *             these owners has the fewest elements wins this element
   *             This is kind of like the draft pick model.
   */

  viable_owners = (int *) smalloc(MAX_NEIGHBOR_NODES*sizeof(int));
  for ( e=0; e<exo->num_elems; e++)
    {

      if ( element_owner[e] == -1 )
	{

	  /*
	   * Clear list of candidates from the last element election...
	   */

	  num_viable_owners=0;
	  for ( i=0; i<MAX_NEIGHBOR_NODES; i++)
	    {
	      viable_owners[i]=-1;
	    }
	  for ( n=exo->elem_node_pntr[e]; n<exo->elem_node_pntr[e+1]; n++)
	    {
	      node       = exo->elem_node_list[n];
	      node_owner = nodal_own_list[nodal_own_pntr[node]];
	      BULL(node_owner, viable_owners, num_viable_owners);
	    }

	  /*
	   * Pick the owner with the fewest elements so far...
	   */

	  min_val = GIGA;
	  min_dex = -1;

	  for ( i=0; i<num_viable_owners; i++)
	    {
	      if ( element_owner_dist[viable_owners[i]] < min_val )
		{
		  min_val = element_owner_dist[viable_owners[i]];
		  min_dex = i;
		}
	    }

#ifdef DEBUG
	  fprintf(stderr, 
		  "Phase III - Assigning element [%d] to set/proc %d\n",
		  e, viable_owners[min_dex]);
#endif

	  element_owner[e]   = viable_owners[min_dex];
	  element_owner_dist[viable_owners[min_dex]]++;
	  num_unassigned_elems--;
	}
    }

  if ( num_unassigned_elems != 0 )
    {
      EH(-1, "Something is wrong. Three phases didn't kill all elems.");
    }

  err = get_max_val_index(element_owner_dist, num_sets, &max_val, &max_dex);
  EH(err, "Problem finding maximum value and index.");

  err = get_min_val_index(element_owner_dist, num_sets, &min_val, &min_dex);
  EH(err, "Problem finding minimum value and index.");

#ifdef DEBUG
  fprintf(stderr, "Statistics on element distribution:\n");
  fprintf(stderr, "Minimum number of elements per setproc = %d (%d)\n",
	  min_val, min_dex);
  fprintf(stderr, "Maximum number of elements per setproc = %d (%d)\n",
	  max_val, max_dex);
#endif
  /*
   * Wander through the assigned elements. Examine the ownership of the
   * facing neighbors. If different from this element, then this is a boundary
   * element. Apologies to earlier users of that nomenclature.
   */

  for ( e=0; e<exo->num_elems; e++)
    {
      this_owner = element_owner[e];
      element_bnd_stat[e] = FALSE;
      for ( f=exo->elem_elem_pntr[e]; f<exo->elem_elem_pntr[e+1]; f++)
	{
	  neighbor_elem = exo->elem_elem_list[f];
	  if ( neighbor_elem > -1 )
	    {
	      if ( this_owner != element_owner[neighbor_elem] )
		{
		  element_bnd_stat[e] = TRUE;
		}
	    }
	}
    }

  safe_free(viable_owners);
  return;
}


/*
 * build_elem_elem_xtra() - determine more about how elements connect
 *
 * Notes: Requires an Exo_DB database that has already had the elem-elem
 *        connectivity built.
 *
 *        If the element connects to nothing (-1), then the face and twist
 *        of the connecting element are set to -1.
 *
 * Created: 1999/08/04 08:06 MDT pasacki@sandia.gov
 *
 * Revised:
 */

void
build_elem_elem_xtra(Exo_DB *exo) /* monolith FE db w/ connectivity     (in) */
{
  int length;
  int e;
  int f;
  int face;
  int hi;
  int i;
  int lo;
  int n1;
  int n2;
  int neighbor;
  int snl_this[MAX_NODES_PER_SIDE]; /* Side Node List */
  int snl_other[MAX_NODES_PER_SIDE]; /* Side Node List */
  int twist;
  int where;

  if ( ! exo->elem_elem_conn_exists )
    {
      EH(-1, "Cannot build xtra w/o basic info.");
    }

  /*
   * Allocate some space...
   */

  length = exo->elem_elem_pntr[exo->num_elems];

  exo->elem_elem_face = (int *) smalloc(length*sizeof(int));
  exo->elem_elem_twst = (int *) smalloc(length*sizeof(int));

  /*
   * Initialize...
   */

  for ( i=0; i<length; i++)
    {
      exo->elem_elem_face[i] = -1;
      exo->elem_elem_twst[i] = -1;
    }

  /*
   * The names of the neighbor's face can be built by find the neighbor
   * and searching for which of its faces connects to the originating
   * element.
   */

  for ( e=0; e<exo->num_elems; e++)
    {
      
      /*
       * Clean out the side node lists...
       */

      for ( i=0; i<MAX_NODES_PER_SIDE; i++)
	{
	  snl_this[i]  = -1;
	  snl_other[i] = -1;
	}

      /*
       * Look at each face of this element and the neighbor, if any, that
       * is there.
       */

      for ( f=exo->elem_elem_pntr[e],face=0; 
	    f<exo->elem_elem_pntr[e+1]; 
	    f++,face++)
	{
	  neighbor = exo->elem_elem_list[f];

	  if ( neighbor > -1 )
	    {

	      /*
	       * Now, go the neighbor and find out the index of this element
	       * in its list!
	       */
	      lo = exo->elem_elem_pntr[neighbor];
	      hi = exo->elem_elem_pntr[neighbor+1];

	      where = in_list(e, exo->elem_elem_list + lo, hi-lo);

	      EH(where, "Didn't find my recipricol element!");

	      /*
	       * Record where originating element was found!
	       */

	      exo->elem_elem_face[f] = where;

	      /*
	       * Now load up the nodes on this element, this face and compare
	       * to other element, other face collection of same nodes. How
	       * much cyclic twist of the OTHER element is required to get
	       * the low nodes coincident?
	       */

	      n1 = build_side_node_list(e, face, exo, snl_this);

	      n2 = build_side_node_list(neighbor, where, exo, snl_other);

	      if ( n1 != n2 )
		{
		  EH(-1, "Difft number of nodes on facing elems!?!");
		}

	      twist = in_list(snl_this[0], snl_other, n2);

	      EH(twist, "Low node not found in facing neighbor!");

	      exo->elem_elem_twst[f] = twist;
	    }
	}
    }

  return;
}


/*
 * build_elem_elem_dpi() -- construct elem elem distributed proc info
 *
 *
 * Notes:
 *
 *	Uses the connectivity information from the monolith db (mexo)
 * and the ownership assignement (owner) and the polylithic piece to
 * allocate and fill selected arrays in d.
 *
 * Created: 1999/08/10 14:42 MDT pasacki@sandia.gov
 */

void 
build_elem_elem_dpi(Exo_DB *monolith, 
		    int *elem_owner_list, 
		    Exo_DB *piece, 
		    Dpi *dpi)
{
  int econ_global;
  int elem;
  int elem_global;
  int hi;
  int i;
  int index;
  int len;
  int lo;
  int lo_global;
  int ne;


  if ( ! monolith->elem_elem_conn_exists )
    {
      EH(-1, "Need monolithic connectivity.");
    }

  if ( ! piece->elem_elem_conn_exists )
    {
      EH(-1, "Need decomposed e-e connectivity.");
    }

  ne  = piece->num_elems;

  len = piece->elem_elem_pntr[ne];

  dpi->len_elem_elem_list = len;

  /*
   * Allocate and initialize...
   */

  dpi->elem_owner = (int *) smalloc(ne*sizeof(int));

  for ( i=0; i<ne; i++)
    {
      dpi->elem_owner[i] = UNASSIGNED;
    }

  dpi->elem_elem_list_global = (int *) smalloc(len*sizeof(int));
  dpi->elem_elem_twst_global = (int *) smalloc(len*sizeof(int));
  dpi->elem_elem_face_global = (int *) smalloc(len*sizeof(int));
  dpi->elem_elem_proc_global = (int *) smalloc(len*sizeof(int));

  for ( i=0; i<len; i++)
    {
      dpi->elem_elem_list_global[i] = UNASSIGNED;
      dpi->elem_elem_twst_global[i] = UNASSIGNED;
      dpi->elem_elem_face_global[i] = UNASSIGNED;
      dpi->elem_elem_proc_global[i] = UNASSIGNED;
    }

  /*
   * The little piece has to know who owns all of the elements that it
   * traverses. Mostly this should be itself, but occassinally it will
   * be other set/procs. This array gives a local glimpse of that division
   * or assignment.
   *
   * Note that elem_owner[i] gives the owner of this element, while
   * down below, elem_elem_proc_global[] gives the owners of the connecting
   * elements (ordered globally).
   */

  for ( i=0; i<ne; i++)
    {
      index = dpi->elem_index_global[i];
      if ( index < 0 || index > monolith->num_elems - 1 )
	{
	  fprintf(stderr, "Bad index attempt %d\n", index);
	}
      dpi->elem_owner[i] = elem_owner_list[dpi->elem_index_global[i]];
    }

  /*
   * How all the global elements fit together. The names, faces, orientations 
   * and owning processors of elements adjoining a given element are shoveled
   * into a chunk for this processor...
   *
   * Note!!! The information is from the global perspective! Thus, the local
   * face ordering may well be different from the local face ordering, etc.
   */

  for ( elem=0; elem<ne; elem++)
    {
      elem_global = dpi->elem_index_global[elem];
      lo_global   = monolith->elem_elem_pntr[elem_global];

      lo = piece->elem_elem_pntr[elem];
      hi = piece->elem_elem_pntr[elem+1];
      for ( i=0; i<hi-lo; i++)
	{
	  dpi->elem_elem_list_global[lo+i] = 
	    monolith->elem_elem_list[lo_global+i];

	  dpi->elem_elem_twst_global[lo+i] = 
	    monolith->elem_elem_twst[lo_global+i];

	  dpi->elem_elem_face_global[lo+i] = 
	    monolith->elem_elem_face[lo_global+i];

	  econ_global = monolith->elem_elem_list[lo_global+i];

	  dpi->elem_elem_proc_global[lo+i] = ( econ_global == -1 ) ?
	    -1 : elem_owner_list[econ_global];
	}
    }

  /*
   * Verification...everything got assigned?
   */

  for ( i=0; i<ne; i++)
    {
      if ( dpi->elem_owner[i] < 0 )
	{
	  EH(-1, "This proc thinks an element is unclaimed!");
	}
    }

  for ( i=0; i<len; i++)
    {
      if ( dpi->elem_elem_list_global[i] == UNASSIGNED )
	{
	  EH(-1, "Unassigned piece elem_elem_list_global[] !");
	}
    }


  return;
}
