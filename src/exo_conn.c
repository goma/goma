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
 
/* exo_conn.c -- build additional connectivity information for EXODUS struct
 *
 * Description: Build additional connectivity maps for an EXODUS II mesh.
 *              These include:
 *
 *			(o) total element->node connectivity (not just
 *			    per element block basis)
 *
 *			(o) total node->element connectivity
 *
 *			(o) element->element connectivity
 *
 *              Notes:
 *
 *		[0] The PATRAN/EXODUS II convention of face naming is 
 *                  followed. The only deviation is that for convenience
 *                  in C, the sequencing is typically 0,...,n-1 instead of 
 *		    the FORTRAN convention of 1,...,n.
 *
 *		[1] Since GOMA sometimes uses a different convention we
 *		    inherited from rf_salsa, be careful about side names.
 *		    The "id_side" that is sometimes used in boundary condition
 *		    procedures is NOT the same as the side id from the
 *		    PATRAN/EXODUS II manual!!!
 *
 *		[2] Everything is predicated on building pointers and lists.
 *
 *              [3] Assume we're only interested initially in elements that
 *                  are connected by faces that are just one dimension less.
 *		    That is, for 3D elements, faces of finite area.
 *		    For 2D, edges of finite length.
 *
 *		[4] If no other element faces the given element, then denote
 *		    the element ID with OUTER_SPACE. Also, reserve the 
 *                  possibility that the neighbor element is off-processor.
 *		    For now, serially, this will be the same as outer space,
 *		    but we can identify them easily as, say, elem name < 0.
 *
 *
 *   "Dedicated to Bob Cochran, who inspired me to use pointers and lists
 *    to handle just about any kind of connectivity."
 *
 *
 *
 * Created: 1998/04/07 11:49 MDT pasacki@sandia.gov
 * 
 * Revised: 1998/04/14 08:05 MDT pasacki@sandia.gov
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "std.h"
#include "mm_eh.h"
#include "el_elm.h"
#include "rf_allo.h"
#include "rf_mp.h"
#include "exo_struct.h"
#include "dpi.h"
#include "el_elm_info.h"
#include "exo_conn.h"
#include "rd_exo.h"
#include "rd_mesh.h"
#include "rf_fem.h"
#include "rf_solver.h"
#include "rf_solver_const.h"

#define GOMA_EXO_CONN_C

/*
 * Prototype declarations of static functions defined in this file.
 */

static int int_intersect
(int *,			/* a - first integer list	        (in) */
       int *,			/* b - second integer list              (in) */
       const int ,		/* len_a - length of 1st integer list   (in) */
       const int ,		/* len_b - length of 2nd  integer list  (in) */
       int *,			/* ia - indeces intersections, 1st list (out)*/
       int *);			/* ib - indeces intersections, 2nd list (out)*/

static int sides2nodes
(const int ,		/* face - face number 0,1,...,num_faces-1    */
       const int ,		/* shape - one of LINE_SEGMENT, etc.         */
       int *);			/* local_indeces - gets filled w/ right ones */

/*
 * From the PATRAN/EXODUS convention, generate a list of the GOMA id_sides
 * that are consistent with their definitions in el_elm_info.c
 *
 * This is the one place where these arrays are defined. They are usually
 * accessible through their extern declarations in el_elm.h
 *
 * Usage: The idea is that the ordering of element-element connectivity
 *        is based on element sides or faces as indicated by PATRAN.
 *        Given the 0-based side index, we can find the "id_side" used
 *        in GOMA from this area. Those id_sides are used frequently to
 *        help determine which way is "out" of an element...
 *
 * 1D: 1:(x=-1), 2:(x=1)
 *
 * 2D: 1:(x=-1), 2:(x=1), 3:(y=-1), 4:(y=1)
 *
 * 3D: 1:(x=-1), 2:(x=1), 3:(y=-1), 4:(y=1), 5:(z=-1), 6:(z=1)
 */

int ID_Side_Quad[4] = { 3, 2, 4, 1};

int ID_Side_Hex[6] = { 3, 2, 4, 1, 5, 6};

Spfrtn sr;

/*
 * Prototype declarations of functions defined in this file.
 */
static int get_num_faces	/* exo_conn.c */
(char *);		/* elem_type */

#if FALSE /* ................just for demo and debugging...................*/
static void demo_elem_node_conn
(Exo_DB *);		/* exo - pntr to EXODUS II FE database */

extern void demo_node_elem_conn
(Exo_DB *);		/* exo - pntr to EXODUS II FE database */

static void demo_elem_elem_conn
(Exo_DB *);		/* exo - pntr to EXODUS II FE database */
#endif

/*
 * The element node connectivity is already built in rd_mesh.c, so
 * here we just assign the proper pointers to existing arrays that already
 * exist. If they don't, then we go ahead and build them and make the
 * arrays reference the same data.
 */

void
build_elem_node(Exo_DB *exo)
{
  int ebi;
  int elem;
  int i;
  int j;
  int k;
  int length;
  int npe;
  int node;
  
  
  if ( exo->elem_ptr != NULL && exo->node_list != NULL )
    {
      exo->elem_node_conn_exists = TRUE;
      exo->elem_node_pntr        = exo->elem_ptr;
      exo->elem_node_list        = exo->node_list;
      return;
    }

  /*
   * Connectivity -- consolidate from a per element block description into
   *		    a total description for all of the elements that this
   *		    processor sees from all element blocks.
   */
 
  length = 0;
  for ( i=0; i<exo->num_elem_blocks; i++)
    {
      length += (exo->eb_num_elems[i])*(exo->eb_num_nodes_per_elem[i]);
    }

  exo->node_list = (int *) smalloc(length*sizeof(int));
  exo->elem_ptr  = (int *) smalloc((exo->num_elems+1)*sizeof(int));

  exo->elem_node_conn_exists = TRUE;
  exo->elem_node_pntr        = exo->elem_ptr;
  exo->elem_node_list        = exo->node_list;

  /*
   * Load the per element block connectivities into the big connectivity.
   */

  elem = 0;
  node = 0;
  exo->elem_ptr[0] = 0;

  for ( ebi=0; ebi<exo->num_elem_blocks; ebi++)
    {
      npe = exo->eb_num_nodes_per_elem[ebi]; 
      k = 0;
      for ( i=0; i<exo->eb_num_elems[ebi]; i++ )
	{
	  for ( j=0; j<npe; j++ )
	    {
	      exo->node_list[node] = exo->eb_conn[ebi][k];
	      k++;
	      node++;
	    }
	  exo->elem_ptr[elem+1] = exo->elem_ptr[elem] + npe;
	  elem++;
	}
    }  

  return;
}

/*
 * Build the node -> node connectivity given the elem -> node connectivity
 * and the node -> elem connectivity.
 *
 * Created: 1999/01/27 08:38 MST pasacki@sandia.gov
 *
 * Revised:
 */

void
build_node_node(Exo_DB *exo)
{
  int curr_list_size;
  int dude;
  int e;
  int elem;
  int end;
  int i;
  int length;
  int max;
  int n;
  int node;
  int npe;
  int start;
  int this_node;
  int total_list_size;
  int n_elem;
  
  int *list;
  /*
   * Don't even attempt to do this without adequate preparation.
   */

  if ( ! exo->node_elem_conn_exists )
    {
      EH(-1, "node_elem conn must exist before node_node can be built.");
    }

  if ( ! exo->elem_node_conn_exists )
    {
      EH(-1, "elem_node conn must exist before node_node can be built.");
    }


  /*
   * Intial memory allocation. The length of the node list will be less than
   * this number by probably about a factor of 2. We will realloc more
   * precisely later. Repeating - this is an overestimate.
   */

  length = 0;

  max    = -1;

  for ( node=0; node<exo->num_nodes; node++)
    {
      this_node = 0;
      for ( e=exo->node_elem_pntr[node]; e<exo->node_elem_pntr[node+1]; e++)
	{
	  elem       = exo->node_elem_list[e];
	  npe        = exo->elem_node_pntr[elem+1] - exo->elem_node_pntr[elem];
	  length    += npe;
	  this_node += npe;
	}
      if ( this_node > max )
	{
	  max = this_node;
	}
    }

  exo->node_node_conn_exists = TRUE;
  exo->node_node_pntr        = (int *) smalloc((exo->num_nodes+1)*sizeof(int));
  exo->node_node_list        = (int *) smalloc(length*sizeof(int));
  exo->centroid_list         = (int *) smalloc((exo->num_nodes)*sizeof(int));

  /*
   * To build unique lists of nodes we'll need some little buffer arrays, too.
   */

  list = (int *) smalloc(max*sizeof(int));

  /*
   * Loop through each node and build a list of all the nodes to which it
   * is connected. Flatten the list by extracting duplicate entries. Finally,
   * sort it and attach it to the global list.
   */

  exo->node_node_pntr[0] = 0;

  total_list_size = 0;

  for ( node=0; node<exo->num_nodes; node++)
    {

      exo->centroid_list[node] = -1;

      /* Clear out any old garbage... */

      for ( i=0; i<max; i++)
	{
	  list[i] = -1;
	}

      curr_list_size = 0;

      for ( e=exo->node_elem_pntr[node]; e<exo->node_elem_pntr[node+1]; e++)
	{
	  elem = exo->node_elem_list[e];

	  start = exo->elem_node_pntr[elem];
	  end   = exo->elem_node_pntr[elem+1];
	  for ( n=start; n<end; n++)
	    {
	      dude = exo->elem_node_list[n];

	      /*
	       * If this dude is not in the list, then add it and make the list
	       * suitably larger.
	       */

	      if ( -1 == in_list(dude, 0, curr_list_size, list) )
		{
		  list[curr_list_size] = dude;
		  curr_list_size++;
		}
	    }
	}

      /* 
       * add centroid nodes of surrounding elements if node is a centroid node
       *  This feature is used in discontinuous Galerkin upwinding....
       */

      if (  ( exo->node_elem_pntr[node+1] - exo->node_elem_pntr[node] ) == 1
	    && ( Use_DG || TRUE ) )  /* Disable Use_DG switch */
	{
	  /* 
	   * This node belongs to only one element 
	   */

	  elem = exo->node_elem_list[ exo->node_elem_pntr[node] ];

	  if (  ( n = centroid_node( Elem_Type( exo, elem ) ) ) != -1 )
	    {

	      if ( node  == exo->elem_node_list[ exo->elem_node_pntr[elem] + n ] )
		{
		  /* This is a centroid node */
		  
		  exo->centroid_list[node] = elem;

		  for( e = exo->elem_elem_pntr[elem]; e < exo->elem_elem_pntr[elem+1]; e++ )
		    {
		      if ( (n_elem = exo->elem_elem_list[e]) != -1 )
			{
			  
			  n = centroid_node ( Elem_Type( exo , n_elem ) );
			  
			  EH(n, "No centroid node exists for element type ");
		      
			  dude = exo->elem_node_list[ exo->elem_node_pntr[n_elem] + n ];
			  
			  if ( -1 == in_list(dude, 0, curr_list_size, list) )
			    {
			      list[curr_list_size++] = dude;
			    }
			}
		    }
		}
	    }
	}


      /*
       * Sort the list before appending to the big concatenated list...
       */

      integer_sort(curr_list_size, list);

      /*
       * No, we're assuming we're never in danger of overrunning the buffer
       * since we were so profligate at the beginning...
       */

      for ( i=0; i<curr_list_size; i++, total_list_size++)
	{
	  exo->node_node_list[total_list_size] = list[i];
	}

      exo->node_node_pntr[node+1] = total_list_size;
    }

  exo->node_node_list = (int *) realloc(exo->node_node_list,
					total_list_size*sizeof(int));


  safe_free(list);

  return;
}


void
build_node_elem(Exo_DB *exo)
{
  int e;
  int i;
  int index;
  int length;
  int n;
  int node;
  char err_msg[MAX_CHAR_ERR_MSG];

  /*
   * If the regular element->node connectivity hasn't been built, then
   * this exercise won't work.
   */

  if ( ! exo->elem_node_conn_exists )
    {
      EH(-1, "Build elem->node before node->elem.");
      return;
    }

  /*
   * Proposition: The length of the elem_node_list[] and the node_elem_list[] 
   * arrays should be the same.
   */

  length = exo->elem_node_pntr[exo->num_elems];

  exo->node_elem_pntr = (int *) smalloc((exo->num_nodes+1)*sizeof(int));
  exo->node_elem_list = (int *) smalloc(length*sizeof(int));

  /*
   * Initialize, then count occurrences of each node in the 
   * element->node connectivity list.
   */

  for ( i=0; i<=exo->num_nodes; i++)
    {
      exo->node_elem_pntr[i] = 0;
    }

  for ( i=0 ; i<length; i++ )
    {
      node = exo->elem_node_list[i];
      exo->node_elem_pntr[node+1]++;	
    }

  /*
   * Transform occurence profile into pointer list.
   */

  for ( i=0; i<exo->num_nodes; i++)
    {
      exo->node_elem_pntr[i+1] += exo->node_elem_pntr[i];
    }

  if ( exo->node_elem_pntr[exo->num_nodes] != length )
    {
      EH(-1, "Inconsistency during connectivity inversion.");
    }

  for ( i=0; i<length; i++)
    {
      exo->node_elem_list[i] = -1;
    }

  /*
   * Traverse the elem->node list. For each element, look at its collection
   * of nodes. Use the node's name to place this element into our new element
   * list.
   */

  for ( e=0; e<exo->num_elems; e++)
    {
      for ( n=exo->elem_node_pntr[e]; n<exo->elem_node_pntr[e+1]; n++)
	{
	  node = exo->elem_node_list[n];

	  /*
	   * Start looking in el[] for the next available slot to inject the 
	   * element name. 
	   */

	  index = exo->node_elem_pntr[node];

	  /*
	   * If this position is filled, increment to see if we can inject
	   * later.
	   */

	  while ( exo->node_elem_list[index] != -1 && 
		  index < exo->node_elem_pntr[node+1] )
	    {
	      index++;
	    }

	  if ( index >= exo->node_elem_pntr[node+1] )
	    {
	      sr = sprintf(err_msg, "No free spot elem [%d], node [%d].",
			  e, node);
	      EH(-1, err_msg);
	    }
	  exo->node_elem_list[index] = e;
	}
    }

  /*
   * Verify the element list is completely filled with meaningful element
   * numbers...
   */

  for ( i=0; i<length; i++)
    {
      if ( exo->node_elem_list[i] == -1 )
	{
	  EH(-1, "Node->elem connectivity has holes!");
	}
    }

  exo->node_elem_conn_exists = TRUE;

  return;
}

void
build_elem_elem(Exo_DB *exo)
{
  int ce;
  int count;
  int e;
  int ebi;
  int elem;
  int ename;
  int face;
  //int his_dim, her_dim;
  int i, j;
  int index;
  int len_curr;
  int len_prev;
  int len_intr;
  int length, length_new, num_faces;
  int iel;
  int ioffset;
  int n;
  int neighbor_name = -1;
  int node;
  int num_elem_sides;
  int num_nodes;
  int snl[MAX_NODES_PER_SIDE];	/* Side Node List - NOT Saturday Night Live! */
  char err_msg[MAX_CHAR_ERR_MSG];

  /*
   * Integer arrays used to find intersection sets of node->element lists.
   */

  int prev_set[MAX_EPN];	/* list of elements attached to previous node*/
  int curr_set[MAX_EPN];	/* list of elements attached to "this" node */


  int ip[MAX_EPN];		/* indeces of hits for prev_set[] */
  int ic[MAX_EPN];		/* indeces of hits for curr_set[] */

  /*
   * If the element->node and node->element connectivities have not been
   * built, then we won't be able to do this task.
   */

  if ( ! exo->elem_node_conn_exists ||
       ! exo->node_elem_conn_exists )
    {
      EH(-1, "Build elem->node before node->elem.");
      return;
    }
 /*
  * The number of elements connected via conventional faces may be deduced
  * from the number of elements and their type.
  */

 exo->elem_elem_pntr = (int *) smalloc((exo->num_elems+1)*sizeof(int));

 length = 0;

 for ( i=0; i<exo->num_elem_blocks; i++)
   {
     length += exo->eb_num_elems[i] * get_num_faces(exo->eb_elem_type[i]);
   }



 exo->elem_elem_list = (int *) smalloc(length*sizeof(int));

 /*
  * Initialize...
  */

 for ( i=0; i<length; i++)
   {
     exo->elem_elem_list[i] = UNASSIGNED_YET;
   }

 /*

 elem = 0;
 for ( ebi=0; ebi<exo->num_elem_blocks; ebi++)
   {
     num_elem_sides = get_num_faces(exo->eb_elem_type[ebi]);     
     for ( e=0; e<exo->eb_num_elems[ebi]; e++)
       {
	 exo->elem_elem_pntr[elem] = count;
	 elem++;
	 count += num_elem_sides;
       }
   }
 */

 /*
  * Walk through the elements, block by block.
  */

 count = 0;
 elem  = 0;
 for ( ebi=0; ebi<exo->num_elem_blocks; ebi++)
   {
     num_elem_sides = get_num_faces(exo->eb_elem_type[ebi]);

     for ( e=0; e<exo->eb_num_elems[ebi]; e++,elem++)
       {
	 exo->elem_elem_pntr[elem] = count;
	 count += num_elem_sides;
	 
	 /*
	  * Look at each side of the element, collecting a unique
	  * list of integers corresponding to the minimum number of nodes
	  * needed to identify an entire side. 
	  *
	  * Typically, the same number of nodes as space dimensions are
	  * needed, with exceptions being the various "sides" of shells,
	  * beams and trusses...
	  */

	 for ( face=0; face<num_elem_sides; face++)
	   {

	     /*
	      * Given the element and the face construct the
	      * list of node numbers that determine that side.
	      */

	     /*
	      * Later, we might not need *all* the nodes on a side,
	      * particularly for high order elements. It may suffice
	      * to check only as many nodes as space dimensions that
	      * the element lives in...
	      */

	     num_nodes = build_side_node_list(elem, face, exo, snl);

	     /*
	      * Cross check: for each node in the side there is a list
	      * of elements connected to it. Beginning with all the
	      * elements connected to the first node (except for this given
	      * element), cross check with all the elements connected with
	      * the 2nd node to build an intersection set of one element.
	      */

	     for ( i=0; i<MAX_EPN; i++)
	       {
		 prev_set[i] = -1;
		 curr_set[i] = -1;
	       }
	     len_prev  = 0;
	     len_curr  = 0;
	     len_intr  = 0;

	     for ( n=0; n<num_nodes; n++)
	       {
		 /*
		  * Copy this node's element list into a clean "curr_set" array
		  * that will be intersected with any previously gathered
		  * lists of elements that qualify as promiscuously in
		  * contact with nodes...
		  */

		 node = snl[n];

		 for ( i=0; i<MAX_EPN; i++)
		   {
		     curr_set[i] = -1;
		   }

		 len_curr  = 0;

		 for ( ce=exo->node_elem_pntr[node]; 
		       ce<exo->node_elem_pntr[node+1]; ce++)
		   {
		     ename = exo->node_elem_list[ce];
		     /*
		      * Go ahead and accumulate the self element name
		      * just as a consistency check....
		      */

		     /*
		     if ( ename != e )
		       {
		       }
		     */

		     /*
		      * PKN: The current Goma use of ->elem_elem...
		      * is such that this connectivity should list
		      * connections like QUAD-BAR or HEX-SHELL.
		      * So, I'll add this dimension matching conditional
		      */

		     /* PRS (Summer 2012): Need to change this for shell stacks which 
			have the same dim*/

		     /* We need however to consider a special case (as of 8/30/2012
		      * this is a SHELL-on-SHELL stack. Viz. two materials, each a shell material
		      * which share not a side but a face.   Since faces of shells are sides
		      * in patran speak, we need some special logic.   We need to avoid adding
		      * the friend shell element (neighboring material) to the current shell element
		      * even though each material has the same number of sides.   
		      * Here goes   (BTW, I cannot find max-nodes-per-element anywhere!!!!)
		      */
		    
		     int shell_on_shell = 0; int flippy_flop = 0; 
		     int nbr_ebid; int nbr_num_elem_sides;
		     nbr_ebid = fence_post(ename, exo->eb_ptr,
                                           exo->num_elem_blocks+1);
                     EH(nbr_ebid, "Bad element block ID!");
                     nbr_num_elem_sides = get_num_faces(exo->eb_elem_type[nbr_ebid]);
		     
		     shell_on_shell = 0; 
		     flippy_flop = 0;

		     if (exo->eb_id[ebi] < 100 && exo->eb_id[nbr_ebid] >= 100) flippy_flop=1;
		     if (exo->eb_id[ebi] >= 100 && exo->eb_id[nbr_ebid] < 100) flippy_flop=1;

		     if ((nbr_ebid != ebi) && 
			 (strstr(exo->eb_elem_type[nbr_ebid], "SHELL")) &&
			 (strstr(exo->eb_elem_type[ebi],      "SHELL")) &&
			 flippy_flop)  shell_on_shell = 1;

		     // his_dim = elem_info(NDIM, exo->eb_elem_itype[ebi]);
		     // her_dim = elem_info(NDIM, exo->eb_elem_itype[exo->elem_eb[ename]]);
		     // if( his_dim == her_dim )
		     if (nbr_num_elem_sides == num_elem_sides && !shell_on_shell)
		       {
			 curr_set[len_curr] = ename;
			 len_curr++;
		       }
		   }
		 
		 /*
		  * The first node is special - we'll just compare
		  * it with itself by making the "previous set" just the
		  * same as the current set...
		  */
		 

		 if ( n == 0 )
		   {
		     for ( i=0; i<MAX_EPN; i++)
		       {
			 prev_set[i] = curr_set[i];
		       }
		     len_prev = len_curr;
		   }


		 /*
		  * First, clean the intersection list and the list of
		  * hit indeces in the previous and current lists.
		  *
		  * Then find the intersection of the previous and current 
		  * sets of elements attached to the previous and current
		  * nodes...
		  */

		 for ( i=0; i<MAX_EPN; i++)
		   {
		     ip[i]       = -1;
		     ic[i]       = -1;
		   }
		 len_intr = 0;

		 len_intr = int_intersect(prev_set, curr_set, len_prev,
					  len_curr, ip, ic);

		 /*
		  * Now, let's make the intersection set the next previous
		  * set of elements, a standard for comparison. We should
		  * eventually boil down to either one or zero elements
		  * that qualify...
		  */
		 
		 for ( i=0; i<MAX_EPN; i++)
		   {
		     prev_set[i] = -1;
		   }

		 for ( i=0; i<len_intr; i++)
		   {
		      prev_set[i] = curr_set[ic[i]];
		   }

		 len_prev = len_intr;
	       }


	     /*
	      * Now consider the different cases.
	      */

	     if ( len_intr == 2 )
	       {
		 /*
		  * The boiled list contains self and one other element.
		  */

		 if ( prev_set[0] == elem )
		   {
		     neighbor_name = prev_set[1];
		   }
		 else
		   {
		     neighbor_name = prev_set[0];
		     if ( prev_set[1] != elem )
		       {
			 sr = sprintf(err_msg, 
				      "2 elems ( %d %d ) 1 should be %d!",
				      prev_set[0], prev_set[1], elem);
			 EH(-1, err_msg);
		       }
		   }
	       }
	     else if ( len_intr == 1 &&  prev_set[0] == elem )
	       {
		 /*
		  * The boiled list has one member, this element.
		  * 
		  * The face must connect either to outer space or to
		  * another processor.
		  */
		 if ( Num_Proc == 1 )
		   {
		     neighbor_name = -1;
		   }
		 else
		   {
		     neighbor_name = -1;

		     /*
		      * I am going to punt for now. Later, revisit this
		      * condition and insert code to check for neighbor
		      * processors containing all the same face nodes.
		      *
		      * EH(-1, "Not done yet..."); 
		      *
		      */

		     /*
		      * Check if ALL the nodes on this face belong
		      * to another processors list of nodes. I.e., the
		      * node must all be in the external node list of
		      * and belong to the same external processor.
		      */
		   }
	       }

	     /*
	      * Pathological cases that normally should not occur....
	      */

	     else if ( len_intr == 0 )
	       {
		 sr = sprintf(err_msg, "Elem %d, face %d should self contain!",
			      elem, face);
		 EH(-1, err_msg);
	       }
	     else if ( len_intr == 1 && prev_set[0] != elem )
	       {
		 sr = sprintf(err_msg, 
			      "Elem %d, face %d only connects with elem %d ?",
			      elem, face, prev_set[0]);
		 EH(-1, err_msg);
	       }
	     else
	       {
		 sr = sprintf(err_msg, 
   	         "Unknown elem-elem connection elem %d, face %d, len_intr=%d",
			      elem, face, len_intr);
		 WH(-1, err_msg);
	       }

	     /*
	      * Now we know how to assign the neighbor name for this face
	      * of the element.
	      */

	     index = exo->elem_elem_pntr[elem] + face;
	     exo->elem_elem_list[index] = neighbor_name;
	     
	   } /* end face loop this elem */

       } /* end elem loop this elemblock */

   } /* end elem block loop */

 exo->elem_elem_pntr[exo->num_elems] = count; /* last fencepost */

 exo->elem_elem_conn_exists = TRUE;

 if (Linear_Solver == FRONT)
   {
     /*
      * Now that we have elem_elem_pntr and elem_elem_list for our parallel
      * world, we are going to use them also for optimal element bandwidth
      * reduction ordering.    We will use METIS, but METIS requires the CSR
      * format, which is compressed, viz. we need to remove the -1s.  Here
      * we go
      */
     
     /* First check for the assumption that all blocks have same number of
	element faces.  Stop if they don't and issue an error to the next
	aspiring developer */
     for ( i=0; i<exo->num_elem_blocks; i++)
       {
	 if(get_num_faces(exo->eb_elem_type[0]) != 
	    get_num_faces(exo->eb_elem_type[i])     )
	   {
	     EH(-1,"Stop! We cannot reorder these elements with METIS with elemement type changes");
	   }
       }
     /* Now begin */
     exo->elem_elem_xadj = (int *) smalloc((exo->num_elems+1)*sizeof(int));
     /*initialize */
     for(e=0; e<exo->num_elems+1 ; e++) 
       { 
      	 exo->elem_elem_xadj[e] = exo->elem_elem_pntr[e]; 
       } 
     
     /* Recompute length of adjacency list by removing external edges */
     
     length_new = 0; 
     for (i = 0; i < length; i++) { 
       if(exo->elem_elem_list[i] != -1) length_new++; 
     } 
     exo->elem_elem_adjncy = alloc_int_1(length_new, -1);
     
     /* Now convert */
     ioffset=0; 
     for(iel = 0; iel < exo->num_elems; iel++) 
       { 
	 /* Big assumption here that all blocks have the same   */
	 /* element type.  Can be furbished later since this is  */
	 /* just for the frontal solver   */
	 
	 num_faces = get_num_faces(exo->eb_elem_type[0]); 
      	 for (i= iel*num_faces; i < (iel+1)*num_faces; i++) 
      	   { 
      	     j = i - ioffset; 
      	     if(exo->elem_elem_list[i] == -1) 
      	       { 
      		 ioffset++; 
      		 for(e=iel +1; e <exo->num_elems+1; e++)exo->elem_elem_xadj[e]--; 
      	       } 
      	     else 
      	       { 
      		 exo->elem_elem_adjncy[j] = exo->elem_elem_list[i]; 
      	       } 
      	   } 
       } 
     /* convert to Fortran style */
     for(e=0; e<exo->num_elems+1 ; e++)  exo->elem_elem_xadj[e]++;
     for ( i=0; i<length_new; i++)   exo->elem_elem_adjncy[i]++;

   } /* End FRONTAL_SOLVER if  */

 /*
  * Verification that every element/face has assigned something besides
  * the initial default value of "unassigned".
  *
  * For your convenience - FORTRAN 1-based numbering.
  */


#if FALSE
 demo_elem_elem_conn(exo);
#endif

 return;
}



int
build_side_node_list(int elem,
		     int face,
		     Exo_DB *exo,
		     int *snl)
{
  int element_type;
  int i;
  int nodes_this_side;
  int shape;
  
  int local_nodeces[MAX_NODES_PER_SIDE];

  /*
   * Assume a canonical ordering to the local nodes in an element according
   * to the PATRAN convention. Faces, too, are conventionally numbered, so
   * all that is needed is the kind of element we have.
   */

  element_type = Elem_Type(exo, elem);
  shape        = type2shape(element_type);

  /*
   * Count up the number of nodes and provide their local 0-based
   * indeces or offsets so their global names may be more easily retrieved.
   */

  nodes_this_side = sides2nodes(face, shape, local_nodeces);
  EH(nodes_this_side, "Problem counting nodes on an element face.");

  for ( i=0; i<nodes_this_side; i++)
    {
      snl[i] = exo->elem_node_list[exo->elem_node_pntr[elem]
				  +local_nodeces[i]];
    }

  return(nodes_this_side);
}

/*
 * sides2nodes - Given element shape and a face, determine local node indeces
 *
 * Notes: Uses the PATRAN convention for face numbering and for local node
 *        numbering for the more common variety of elements, except that 
 *        a C language bias towards 0-based indexing and naming prevails.
 *
 *        Unusual elements like shells, bars, beams, trusses are not
 *        handled.
 *
 *	  Error checking is not strong, the assumption being that the calling
 *        routine is giving good face numbers and a sufficiently large array
 *        to hold the local node indeces.
 *
 *	  Note that higher order nodes are invisible - only those nodes
 *        required to represent the fundamental element shape are put into
 *        the list. Fortunately, PATRAN uses the lowest numbered nodes
 *        to outline the basic element shape so we are safe in ignoring the
 *        higher numbered nodes for quadratic and other high order elements.
 *
 *	  On the faces of 3D elements, looking at the face from the outside,
 *	  order the nodes so that we traverse them in a counter clockwise
 *	  fashion.
 *
 *
 * Created: 1998/04/09 16:33 MDT pasacki@sandia.gov
 *
 * Revised:
 */

static int
sides2nodes(const int face,	/* assume face number 0,1,...,num_faces-1 */
	    const int shape,	/* one of LINE_SEGMENT, etc. */
	    int *local_indeces)	/* get filled with right ones */
{
  int i;
  int nodes_this_side=-1;

  for ( i=0; i<MAX_NODES_PER_SIDE; i++)
    {
      local_indeces[i] = -1;
    }

  switch(shape)
    {
    case LINE_SEGMENT:
      nodes_this_side  = 1;
      local_indeces[0] = face;
      break;

    case TRIANGLE:		/* nice regular counterclockwise node/face...*/
    case TRISHELL:
      nodes_this_side  = 2;
      local_indeces[0] = face;
      local_indeces[1] = (face+1)%3;
      break;

    case QUADRILATERAL:
      nodes_this_side  = 2;
      local_indeces[0] = face;
      local_indeces[1] = (face+1)%4;
      break;

    case SHELL:
      nodes_this_side  = 2;
      local_indeces[0] = face;
      local_indeces[1] = (face+1)%4;
      break;

    case TETRAHEDRON:
      nodes_this_side  = 3;
      switch ( face )
	{
	case 0:
	  local_indeces[0] = 0;
	  local_indeces[1] = 1;
	  local_indeces[2] = 3;
	  break;

	case 1:
	  local_indeces[0] = 1;
	  local_indeces[1] = 2;
	  local_indeces[2] = 3;
	  break;

	case 2:
	  local_indeces[0] = 2;
	  local_indeces[1] = 0;
	  local_indeces[2] = 3;
	  break;

	case 3:
	  local_indeces[0] = 0;
	  local_indeces[1] = 2;
	  local_indeces[2] = 1;
	  break;

	default:
	  EH(-1, "Bad face for TETRAHEDRON.");
	  break;
	}
      break;

    case PRISM:			/* aka "WEDGE" */
      switch ( face )
	{
	case 0:
	  nodes_this_side  = 4;
	  local_indeces[0] = 0;
	  local_indeces[1] = 1;
	  local_indeces[2] = 4;
	  local_indeces[3] = 3;
	  break;

	case 1:
	  nodes_this_side  = 4;
	  local_indeces[0] = 1;
	  local_indeces[1] = 2;
	  local_indeces[2] = 5;
	  local_indeces[3] = 4;
	  break;

	case 2:
	  nodes_this_side  = 4;
	  local_indeces[0] = 2;
	  local_indeces[1] = 0;
	  local_indeces[2] = 3;
	  local_indeces[3] = 5;
	  break;

	case 3:
	  nodes_this_side  = 3;
	  local_indeces[0] = 0;
	  local_indeces[1] = 2;
	  local_indeces[2] = 1;
	  break;

	case 4:
	  nodes_this_side  = 3;
	  local_indeces[0] = 3;
	  local_indeces[1] = 4;
	  local_indeces[2] = 5;
	  break;

	default:
	  EH(-1, "Bad face for PRISM/WEDGE.");
	  break;
	}
      break;

    case HEXAHEDRON:
      nodes_this_side = 4;
      switch ( face )
	{
	case 0:
	case 1:
	case 2:
	case 3:
	  local_indeces[0] = 0 + face;
	  local_indeces[1] = (1 + face)%4;
	  local_indeces[2] = (4 + (1 + face)%4);
	  local_indeces[3] = 4 + face;
	  break;

	case 4:			/* bottom */
	case 5:			/* top */
	  local_indeces[0] = (face-1);
	  local_indeces[1] = (face-1) +(face-5)   + (face-4)  ;
	  local_indeces[2] = (face-1) +(face-5)*2 + (face-4)*2;
	  local_indeces[3] = (face-1) +(face-5)*3 + (face-4)*3;
	  break;

	default:
	  EH(-1, "Bad face for HEXAHEDRON.");
	  break;
	}
      break;

    default:
      EH(-1, "Bad shape.");
      break;
    }

  return(nodes_this_side);
}

static int 
get_num_faces(char *elem_type)
{
  int val=-1;

  if ( elem_type == NULL )
    {
      return(-1);
    }

  /*
   * Check each possibility, the more likely ones first...
   */

  if ( strstr(elem_type, "QUAD") != NULL )
    {
      val = 4;
    }
  else if ( strstr(elem_type, "HEX") != NULL )
    {
      val = 6;
    }
  else if ( strstr(elem_type, "TRI") != NULL )
    {
      val = 3;
    }
  else if ( strstr(elem_type, "TET") != NULL )
    {
      val = 4;
    }
  else if ( strstr(elem_type, "WEDGE") != NULL )
    {
      val = 5;
    }
  else if ( strstr(elem_type, "SHELL") != NULL ) /* in 3D */
    {
      val = 4;
    }
  else if ( strstr(elem_type, "TRUSS") != NULL || 
       strstr(elem_type, "BAR")   != NULL || 
       strstr(elem_type, "BEAM")  != NULL) /* in 2D */
    {
      val = 2;
    }

  return(val);
}

int 
get_exterior_faces( int elem,   
		    int *ex_faces, 
		    const Exo_DB *exo,
		    const Dpi *dpi)
{
  int index;
  int face;
  int num_ext_faces = 0;
  int *elem_elem_list;

  elem_elem_list = exo->elem_elem_list ;

#ifdef PARALLEL
  

  if(Num_Proc > 1 )
    elem_elem_list = dpi->elem_elem_list_global;  /* This global list has the 
						 * element connectivity of the original 
						 * monolith.  Including the actual exterior
						 * faces 
						 */
#endif

  for( index = exo->elem_elem_pntr[elem], face=0; 
       index < exo->elem_elem_pntr[elem+1]; index++, face++)
    {
      if (elem_elem_list[index] == -1 )
	{
	  ex_faces[num_ext_faces] = face;
	  
	  num_ext_faces++;
	}
    }
return( num_ext_faces);
}


#if 0 /* ................................. unused........................ */
static int 
get_elem_dimension(char *elem_type)
{
  int val=-1;

  if ( elem_type == NULL )
    {
      return(-1);
    }

  /*
   * Check each possibility, the more likely ones first...
   */

  if ( strstr(elem_type, "QUAD") != NULL )
    {
      val = 2;
    }

  if ( strstr(elem_type, "HEX") != NULL )
    {
      val = 3;
    }

  if ( strstr(elem_type, "TRI") != NULL )
    {
      val = 2;
    }

  if ( strstr(elem_type, "TET") != NULL )
    {
      val = 3;
    }

  if ( strstr(elem_type, "WEDGE") != NULL )
    {
      val = 3;
    }

  if ( strstr(elem_type, "SHELL") != NULL ) /* in 3D */
    {
      val = 2;
    }

  if ( strstr(elem_type, "TRUSS") != NULL || 
       strstr(elem_type, "BAR")   != NULL ||
       strstr(elem_type, "BEAM")  != NULL) /* in 2D */
    {
      val = 2;
    }

  return(val);
}
#endif /* ................................. unused........................ */

/* int_intersect - provide the number and the indeces of the intersections
 *                 of 2 lists
 *
 *		   Assume both lists are full of unique integers, that none
 *		   appear more than once. If so, you're responsible for what
 *		   happens.
 *
 *		   Adapted from Scott Hutchinson's find_inter.
 *
 * Revised: 1998/04/10 14:51 MDT pasacki@sandia.gov
 */
static int
int_intersect(int *a,		/* first integer list			(in) */
	      int *b,		/* second integer list			(in) */
	      const int len_a,	/* length of first integer list		(in) */
	      const int len_b,	/* length of second integer list	(in) */
	      int *ia,		/* indeces of intersections, first list (out)*/
	      int *ib)		/* indeces of intersections, second list(out)*/
{
  int i;
  int j;
  int num_hit=0;

  for ( i=0; i<len_a; i++)
    {
      for ( j=0; j<len_b; j++)
	{
	  if ( a[i] == b[j] )
	    {
	      ia[num_hit] = i;
	      ib[num_hit] = j;
	      num_hit++;
	    }
	}
    }

  return(num_hit);
}

#if FALSE
static void
demo_elem_node_conn(Exo_DB *exo)
{
  int e;
  int n;
  int node_name;
  int start;
  int end;
  
  if ( ! exo->elem_node_conn_exists )
    {
      EH(-1, "Attempt to access undeveloped elem->node connectivity.");
    }

  for ( e=0; e<exo->num_elems; e++)
    {
      fprintf(stdout, "Element [%d] has nodes: ", e);
      start = exo->elem_node_pntr[e];
      end   = exo->elem_node_pntr[e+1];

      for ( n=start; n<end; n++)
	{
	  node_name = exo->elem_node_list[n];
	  fprintf(stdout, "%d ", node_name);
	}
      fprintf(stdout, "\n");
    }
  return;
}
#endif

#if 0				/* Change to 1 to have a routine that shows
				 * how connectivity works in my life. */
void
demo_node_elem_conn(Exo_DB *exo)
{
  int e;
  int elem_name;
  int n;
  int start;
  int end;
  
  if ( ! exo->node_elem_conn_exists )
    {
      EH(-1, "Attempt to access undeveloped node->elem connectivity.");
    }

  for ( n=0; n<exo->num_nodes; n++)
    {
      fprintf(stdout, "Node [%d] has elements: ", n+1);
      start = exo->node_elem_pntr[n];
      end   = exo->node_elem_pntr[n+1];
      for ( e=start; e<end; e++)
	{
	  elem_name = exo->node_elem_list[e];
	  fprintf(stdout, "%d ", elem_name+1);
	}
      fprintf(stdout, "\n");
    }
  return;
}
#endif

#if FALSE
static void
demo_elem_elem_conn(Exo_DB *exo)
{
  int e;
  int elem_name;
  int en;
  int start;
  int end;
  
  if ( ! exo->elem_elem_conn_exists )
    {
      EH(-1, "Attempt to access undeveloped elem->elem connectivity.");
    }

  for ( e=0; e<exo->num_elems; e++)
    {
      fprintf(stdout, "Element [%d] has elements: ", e);
      start = exo->elem_elem_pntr[e];
      end   = exo->elem_elem_pntr[e+1];
      for ( en=start; en<end; en++)
	{
	  elem_name = exo->elem_elem_list[en];
	  fprintf(stdout, "%d ", elem_name);
	}
      fprintf(stdout, "\n");
    }
  return;
}
#endif

/*
 * Build the node -> node connectivity given the elem -> node connectivity
 * and the node -> elem connectivity.
 *
 * Created: 1999/01/27 08:38 MST pasacki@sandia.gov
 *
 * Revised:
 */

void
brk_build_node_node(Exo_DB *exo)
{
  int curr_list_size;
  int dude;
  int e;
  int elem;
  int end;
  int i;
  int length;
  int max;
  int n;
  int node;
  int npe;
  int start;
  int this_node;
  int total_list_size;
  
  int *list;
  /*
   * Don't even attempt to do this without adequate preparation.
   */

  if ( ! exo->node_elem_conn_exists )
    {
      EH(-1, "node_elem conn must exist before node_node can be built.");
    }

  if ( ! exo->elem_node_conn_exists )
    {
      EH(-1, "elem_node conn must exist before node_node can be built.");
    }


  /*
   * Intial memory allocation. The length of the node list will be less than
   * this number by probably about a factor of 2. We will realloc more
   * precisely later. Repeating - this is an overestimate.
   */

  length = 0;

  max    = -1;

  for ( node=0; node<exo->num_nodes; node++)
    {
      this_node = 0;
      for ( e=exo->node_elem_pntr[node]; e<exo->node_elem_pntr[node+1]; e++)
	{
	  elem       = exo->node_elem_list[e];
	  npe        = exo->elem_node_pntr[elem+1] - exo->elem_node_pntr[elem];
	  length    += npe;
	  this_node += npe;
	}
      if ( this_node > max )
	{
	  max = this_node;
	}
    }

  exo->node_node_conn_exists = TRUE;
  exo->node_node_pntr        = (int *) smalloc((exo->num_nodes+1)*sizeof(int));
  exo->node_node_list        = (int *) smalloc(length*sizeof(int));

  /*
   * To build unique lists of nodes we'll need some little buffer arrays, too.
   */

  list = (int *) smalloc(max*sizeof(int));

  /*
   * Loop through each node and build a list of all the nodes to which it
   * is connected. Flatten the list by extracting duplicate entries. Finally,
   * sort it and attach it to the global list.
   */

  exo->node_node_pntr[0] = 0;

  total_list_size = 0;

  for ( node=0; node<exo->num_nodes; node++)
    {

      /* Clear out any old garbage... */

      for ( i=0; i<max; i++)
	{
	  list[i] = -1;
	}

      curr_list_size = 0;

      for ( e=exo->node_elem_pntr[node]; e<exo->node_elem_pntr[node+1]; e++)
	{
	  elem = exo->node_elem_list[e];

	  start = exo->elem_node_pntr[elem];
	  end   = exo->elem_node_pntr[elem+1];
	  for ( n=start; n<end; n++)
	    {
	      dude = exo->elem_node_list[n];

	      /*
	       * If this dude is not in the list, then add it and make the list
	       * suitably larger.
	       */

	      if ( -1 == in_list(dude, 0, curr_list_size, list) )
		{
		  list[curr_list_size] = dude;
		  curr_list_size++;
		}
	    }
	}

      /*
       * Sort the list before appending to the big concatenated list...
       */

      integer_sort(curr_list_size, list);

      /*
       * No, we're assuming we're never in danger of overrunning the buffer
       * since we were so profligate at the beginning...
       */

      for ( i=0; i<curr_list_size; i++, total_list_size++)
	{
	  exo->node_node_list[total_list_size] = list[i];
	}

      exo->node_node_pntr[node+1] = total_list_size;
    }

  exo->node_node_list = (int *) realloc(exo->node_node_list,
					total_list_size*sizeof(int));


  safe_free(list);

  return;
}
