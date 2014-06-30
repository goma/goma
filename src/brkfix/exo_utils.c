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

/* exo_utils.c -- high level routines for dealing with Exo_DB structures
 *
 * Sometimes the EXODUS II API just is too low level. The routines in this
 * file, along with those in rd_exo.c and wr_exo.c, are meant to make life
 * easier when dealing with EXODUS II data.
 *
 * Created: 1998/12/03 13:01 MST pasacki@sandia.gov
 *
 * Revised:
 */

#define _EXO_UTILS_C

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef STDC_HEADERS
#include <stdlib.h>
#endif

#include <stdio.h>
#include <string.h>

#include "map_names.h"
#include "std.h"
#include "aalloc.h"
#include "eh.h"
#include "exo_struct.h"
#include "utils.h"
#include "exo_utils.h"

/*
 * zero_base() -- push down the element names and node names by one in an
 *                EXODUS II data base. This makes C language zero-based
 *                arrays work more naturally.
 */

void
zero_base(Exo_DB *E)
{
  int eb;
  int i,j;
  int l;
  int length_conn;

  /*
   * 1. Node numbers are named in the connectivity lists for each 
   *    element block. Decrement them.
   */

  for ( eb=0; eb<E->num_elem_blocks; eb++)
    {
      length_conn =  E->eb_num_elems[eb] * E->eb_num_nodes_per_elem[eb];
      for ( l=0; l<length_conn; l++)
	{
	  (E->eb_conn[eb][l])--;
	}
    }

  /*
   * 2. Node numbers are named in the node lists for each node set.
   *    Decrement them.
   */

  for ( l=0; l<E->ns_node_len; l++)
    {
      (E->ns_node_list[l])--;
    }

  /*
   * 3. Element numbers are named in the element list for the side sets.
   *    Decrement these.
   */

  for ( l=0; l<E->ss_elem_len; l++)
    {
      (E->ss_elem_list[l])--;
    }

  /*
   * 4. Node numbers in the special node list for the sidesets.
   *    Decrement these (if they exist).
   */

  if ( E->ss_node_list_exists )
    {
      for ( i=0; i<E->num_side_sets; i++)
	{
	  
	  /*
	   * Instead of doing a side at a time, do the whole list whose
	   * length has been counted up and stored in the last element of this
	   * index.
	   */
	  
	  for ( j=0; j<E->ss_node_side_index[i][ E->ss_num_sides[i] ]; j++)
	    {
	      (E->ss_node_list[i][j])--;
	    }
	}
    }

  /*
   * 5. If a node map exists, the names of the nodes might be too high.
   *    Decrement them. In some cases for very general node maps this
   *    step might be superfluous and even undesirable. Here, the idea is
   *    that the node map probably corresponds to an inherent contiguous
   *    node numbering scheme on a larger mesh. 
   */

  if ( E->node_map_exists )
    {
      for ( i=0; i<E->num_nodes; i++)
	{
	  (E->node_map[i])--;
	}
    }

  /*
   * 6. Likewise, if an element map exists, then it might well correspond to
   *    a contiguous integer sequence on a larger mesh. Those numbers could
   *    either be 1-base or 0-based arrays.
   */

  if ( E->elem_map_exists )
    {
      for ( i=0; i<E->num_elems; i++)
	{
	  (E->elem_map[i])--;
	}
    }

  return;
}

/*
 * one_base() -- push up the element names and node names by one in an
 *               EXODUS II data base. This makes C language zero-based
 *               arrays work more naturally.
 */

void
one_base(Exo_DB *E)
{
  int eb;
  int i,j,l;
  int length_conn;

  /*
   * 1. Node numbers are named in the connectivity lists for each 
   *    element block. Increment them.
   */

  for ( eb=0; eb<E->num_elem_blocks; eb++)
    {
      length_conn =  E->eb_num_elems[eb] * E->eb_num_nodes_per_elem[eb];
      for ( l=0; l<length_conn; l++)
	{
	  (E->eb_conn[eb][l])++;
	}
    }

  /*
   * 2. Node numbers are named in the node lists for each node set.
   *    Increment them.
   */

  for ( l=0; l<E->ns_node_len; l++)
    {
      (E->ns_node_list[l])++;
    }

  /*
   * 3. Element numbers are named in the element list for the side sets.
   *    Increment these.
   */

  for ( l=0; l<E->ss_elem_len; l++)
    {
      (E->ss_elem_list[l])++;
    }

  /*
   * 4. Node numbers in the special node list for the sidesets.
   *    Increment these (if they exist).
   */

  if ( E->ss_node_list_exists )
    {
      for ( i=0; i<E->num_side_sets; i++)
	{
	  
	  /*
	   * Instead of doing a side at a time, do the whole list, whose
	   * length has been counted up and stored in the last element of this
	   * index.
	   */
	  
	  for ( j=0; j<E->ss_node_side_index[i][ E->ss_num_sides[i] ]; j++)
	    {
	      (E->ss_node_list[i][j])++;
	    }
	}
    }



  /*
   * 5. If a node map exists, the names of the nodes might be too high.
   *    Increment them. In some cases for very general node maps this
   *    step might be superfluous and even undesirable. Here, the idea is
   *    that the node map probably corresponds to an inherent contiguous
   *    node numbering scheme on a larger mesh. 
   */

  if ( E->node_map_exists )
    {
      for ( i=0; i<E->num_nodes; i++)
	{
	  (E->node_map[i])++;
	}
    }

  /*
   * 6. Likewise, if an element map exists, then it might well correspond to
   *    a contiguous integer sequence on a larger mesh. Those numbers could
   *    either be 1-base or 0-based arrays.
   */

  if ( E->elem_map_exists )
    {
      for ( i=0; i<E->num_elems; i++)
	{
	  (E->elem_map[i])++;
	}
    }


}

#if 0
/*
 * type2shape() -- convert from general element types into basic shapes
 *
 * The element types tell both the interpolation and the basic shape. Sometimes
 * it's nice to be able to separate the two concepts. Hence, this code to
 * convert the legacy variables into basic shapes.
 */

Element_shape
type2shape(const Element_type et)
{
  Element_shape shape=-1;

  if ( et == BILINEAR_QUAD   ||
       et == C_BILINEAR_QUAD ||
       et == S_BIQUAD_QUAD   ||
       et == S_BIQUAD_QUAD   ||
       et == BIQUAD_QUAD     ||
       et == P0_QUAD         ||
       et == P1_QUAD         ||
       et == BILINEAR_SHELL  ||
       et == BIQUAD_SHELL    )
    {
      shape = QUADRILATERAL;
    }
  else if ( et == TRILINEAR_HEX   ||
	    et == C_TRILINEAR_HEX ||
	    et == S_TRIQUAD_HEX   ||
	    et == TRIQUAD_HEX     ||
	    et == P0_HEX          ||
	    et == P1_HEX )
    {
      shape = HEXAHEDRON;
    }
  else
    {
      EH(-1, "Element type not classifiable into a shape.");
    }

  return(shape);
}
#endif

int
shape2sides(const Element_shape es)
{
  int num_sides=-1;

  switch (es)
    {
    case QUADRILATERAL:
      num_sides = 4;
      break;

    case SHELL:
      num_sides = 4;
      break;

    case HEXAHEDRON:
      num_sides = 6;
      break;

    case TRIANGLE:
      num_sides = 3;
      break;

    case TETRAHEDRON:
      num_sides = 4;
      break;

    case LINE_SEGMENT:
      num_sides = 2;
      break;

    case PRISM:
      num_sides = 5;
      break;
      
    case PYRAMID:
      num_sides = 5;
      break;

    default:
      EH(-1, "How many sides does this new shape have?");
      break;
    }

  return(num_sides);
}

/* get_element_shape() -- given an EXODUS db and element index, return shape
 *
 * Created: 1999/08/16 08:53 MDT pasacki@sandia.gov
 */

Element_shape
get_element_shape(const Exo_DB *exo,
		  const int element)
{
  int found=FALSE;
  int ebi=0;
  Element_shape es=UNDEFINED_ELEMENT_SHAPE;
  char *string;

  if ( element < 0 )
    {
      EH(-1, "Bad element type lookup -- negative element.");
    }

  if ( element > exo->num_elems - 1 )
    {
      EH(-1, "Bad element type lookup -- element number too high.");
    }

  /*
   * Which element block index is this in?
   */

  while ( ! found && ebi < exo->num_elem_blocks )
    {
      found |= ( element >= exo->eb_ptr[ebi] && element < exo->eb_ptr[ebi+1] );
      ebi++;
    }
  
  if ( ebi > exo->num_elem_blocks )
    {
      fprintf(stderr, "eb_ptr[%d] = %d\n", 0, exo->eb_ptr[0]);
      fprintf(stderr, "eb_ptr[%d] = %d\n", 1, exo->eb_ptr[1]);
      fprintf(stderr, "exo->num_elem_blocks = %d\n", exo->num_elem_blocks);
      fprintf(stderr, "exo->num_elems       = %d\n", exo->num_elems);
      fprintf(stderr, "element              = %d\n", element);
      fprintf(stderr, "ebi                  = %d\n", ebi);
      EH(-1, "Could not find element in any element block.");
    }
  
  ebi--;

  string = exo->eb_elem_type[ebi];

  /*
   * Interpret the EXODUS string element type identifier as an element shape...
   */

  if ( 0 == strncmp(string, "QUAD", 4) )
    {
      es = QUADRILATERAL;
    }

  if ( 0 == strncmp(string, "SHELL", 4) )
    {
      es = SHELL;
    }

  if ( 0 == strncmp(string, "HEX", 3) )
    {
      es = HEXAHEDRON;
    }

  if ( 0 == strncmp(string, "TRI", 3) )
    {
      es = TRIANGLE;
    }

  if ( 0 == strncmp(string, "TET", 3) )
    {
      es = TETRAHEDRON;
    }

  if ( 0 == strncmp(string, "BAR", 3) )
    {
      es = LINE_SEGMENT;
    }

  return (es);
}




/*
 * find_element_block()
 * Returns the element block that a given element belongs to
 *
 * Scott A Roberts, 2010-08-24
 */
int 
find_element_block ( Exo_DB *exo,
		     const int e ) {
  int i, j, ei;
  ei = -1; 
  for ( i=0; i<exo->num_elem_blocks; i++) { 
    for ( j=0; j<exo->eb_num_elems[i]; j++) {  
      ei++;
      if ( ei == e ) return i;
    }
  }
  return -1;
}

/*
 * find_element_friends()
 * Returns a list of elements who share a face between HEX-SHELL
 * 
 * Scott A Roberts, 2010-08-24
 */

int 
find_element_friends( Exo_DB *exo,
		      int e,
		      int *friend_list ) {

  int i, num_friends;

  /* Find block of this element */
  int eb = find_element_block(exo,e);
  
  /* How many nodes must be in common to be a friend? */
  int common_node_req;
  if ( !strcmp( exo->eb_elem_type[eb], "SHELL4") || !strcmp(exo->eb_elem_type[eb], "shell4") || 
       !strcmp( exo->eb_elem_type[eb], "HEX8")   || !strcmp(exo->eb_elem_type[eb], "hex8")      ) {
    common_node_req = 4;
  } else {
    return 0;      
  }

  /* Generate list of elements touching this element
     and how many nodes are in common */
  int nip, ni, eip, ei, dup;
  int num_neigh = 0;
  int neigh_elem_id[50];
  int neigh_elem_ct[50];
  for ( nip = exo->elem_node_pntr[e]; nip < exo->elem_node_pntr[e+1]; nip++) {  // Loop through nodes of this element
    ni = exo->elem_node_list[nip];
    for ( eip = exo->node_elem_pntr[ni]; eip < exo->node_elem_pntr[ni+1]; eip++) { // Loop through elements belonging to that node
      ei = exo->node_elem_list[eip];
      dup = 0;
      for ( i = 0; i < num_neigh; i++) {
	if ( neigh_elem_id[i] == ei ) {
	  dup = 1;
	  neigh_elem_ct[i]++;
	}
      }
      if ( dup == 0 ) {
	neigh_elem_id[num_neigh] = ei;
	neigh_elem_ct[num_neigh] = 1;
	num_neigh++;
      }
    }
  }

  /* Loop through neighbors to see if any are friends */
  int nei, neib, fi;
  num_friends = 0;
  for ( nei = 0; nei < num_neigh; nei++) {
    if ( neigh_elem_ct[nei] == common_node_req ) {
      neib = find_element_block(exo,neigh_elem_id[nei]);
      if ( neib != eb ) {
	dup = 0;
	for ( fi = 0; fi < num_friends; fi++) {
	  if (friend_list[fi] == ei) dup = 1;
	}
	if (dup == 0) {
	  friend_list[num_friends] = neigh_elem_id[nei];
	  num_friends++;
	}
      }
    }
  }

  return num_friends;
}


/*
 * find_local_node_number()
 * Returns the local node numbers for a given global node and element
 * 
 * Scott A Roberts, 2010-08-24
 */
int find_local_node_number( Exo_DB *exo, int e, int n ) {
  int i, no, eb, ns, nn;
  eb = find_element_block(exo,e);
  if (strcmp( exo->eb_elem_type[eb], "HEX8")) EH(-1,"Only use find_local_node_number for HEX8");
  ns = exo->elem_node_pntr[e];
  nn = exo->elem_node_pntr[e+1]-ns;
  for ( i = 0; i < nn; i++) {
    no = exo->elem_node_list[ns+i];
    if (no == n) return i;
  }  
  return -1;
}

/*
 * find_edge_connected_nodes()
 * Returns the local node numbers of nodes that are connected
 * to a given local node by an edge.  Only works for HEX8.
 * 
 * Scott A Roberts, 2010-08-24
 */
void find_edge_connected_nodes(int ln, int *lcn) {
  lcn[0] = (ln+1)%4 + (ln-ln%4);
  lcn[1] = (ln+3)%4 + (ln-ln%4);
  lcn[2] = (ln+4)%8;
  return;
}

/* 
 * is_node_in_element()
 * Determines if a given global node is present in an element
 *
 * Scott A Roberts, 2010-08-24
 */
int is_node_in_element(Exo_DB *exo, int n, int e) {
  int i;
  for ( i=exo->elem_node_pntr[e]; i<exo->elem_node_pntr[e+1]; i++) {
    if ( n == exo->elem_node_list[i] ) return 1;
  }
  return 0;
}
