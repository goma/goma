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
 * utils.c -- utilities useful for building goma-like dofmaps, connectivities...
 *
 *
 * Created: 1997/04/10 13:25 MDT pasacki@sandia.gov
 *
 * Revised: 1997/04/10 13:26 MDT pasacki@sandia.gov
 */

#define _UTILS_C

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
#include "exo_utils.h"
#include "utils.h"
#include "string.h"
/*
 * Function declarations of static functions defined here.
 */

static int intcompare PROTO((const void *, const void *));
static int proc_ident PROTO((const void *, const void *));

/*
 * Count up the total number of node-node interactions in a mesh. Include
 * self interactions. This length would be the length of an overall
 * node-node connectivity list.
 */

static char err_msg[MAX_CHAR_ERR_MSG];
static Spfrtn sr=0;

/*
 * Make this variable this visible...
 */

static int keep_len_assignment;
static int *assignment_copy;

/*
 * Definitions of the functions.
 */




int
count_node_node_interactions(int num_nodes,
			     int *node_ptr,
			     int *elem_list,
			     int *elem_ptr,
			     int *node_list)
{
  int elem;
  int eqn_node;
  int i;
  int l;
  int length;
  int m;
  int next_free_spot;
  int var_node;
  int where;

  int *neighbors;

  neighbors = (int *) smalloc( MAX_NEIGHBOR_NODES * sizeof(int));

  length = 0;

  for ( eqn_node=0; eqn_node<num_nodes; eqn_node++)
    {
      next_free_spot = 0;
      
      for ( i=0; i<MAX_NEIGHBOR_NODES; i++)
	{
	  neighbors[i] = -1;
	}

      /*
       * Examine every element containing this eqn node.
       */

      for ( l=node_ptr[eqn_node]; l<node_ptr[eqn_node+1]; l++)
	{
	  elem = elem_list[l];
	  
	  /*
	   * Look at every node that this element contains.
	   */

	  for ( m=elem_ptr[elem]; m<elem_ptr[elem+1]; m++ )
	    {
	      var_node = node_list[m];
	      
	      /*
	       * If this variable node is not in our list, then add it.
	       */

	      where = in_list(var_node, neighbors, next_free_spot);

	      if ( where == -1 )
		{
		  neighbors[next_free_spot] = var_node;
		  next_free_spot++;
		}

	      if ( next_free_spot > MAX_NEIGHBOR_NODES )
		{
		  sr = sprintf(err_msg, "@ node=%d. Increase max neighbors.",
			       eqn_node);
		  EH(-1, err_msg);
		  EH(sr, err_msg);
		}
	    }
	}
#ifdef DEBUG
      fprintf(stderr, "node %d contributes %d\n", eqn_node, next_free_spot);
#endif
      length += next_free_spot;
    }

  free(neighbors);

  return(length);
}



/*
 * in_list() -- return location of 1st appearance of given integer in 
 *              integer list; otherwise return -1
 *
 * Created: 1997/03/19 10:33 MST pasacki@sandia.gov
 */

int 
in_list(int val,		/* what integer value to seek */
	int *start,		/* where to begin looking */
	int length)		/* how many to consider from the start */

{
  int i;
  for ( i=0; i<length; i++, start++)
    {
      if ( *start == val )
	{
	  return(i);
	}
    }
  return(-1);
}

/*
 * findex_mono() -- determine location of appearance of given integer in 
 *                  strictly monotone integer list.
 *
 * Return values:	int n, where n is the index, i.e., array[n] = val.
 *			       -1 if the value is not in the list or something
 *			       else is awry.
 *
 * Notes:  This should generally be faster than in_list for large lists if
 *	   you can guarantee monotonicity, with no more than O(log2(N))
 *         operation count. Bisection is used.
 *
 *	   The integer list must be strictly increasing, otherwise the
 *         results are not guaranteed.
 *
 *	   
 *
 * Created: 1997/05/01 15:58 MDT pasacki@sandia.gov
 *
 * Revised: 
 */

int 
findex_mono(int val,	/* what integer value to seek */
	    int *start,	/* where to begin looking */
	    int length)	/* how far to look from start */
{
  int rtn;
  int half;

#ifdef DEBUG
  fprintf(stderr, "iol: val=%d, start[?]=%d, start[%d]=%d\n",
	  val, start[0], length-1, start[length-1]);
#endif

  if ( length == 0 ) return(-1);

  if ( length == 1 )
    {
      if ( *start == val )
	{
	  return(0);
	}
      else
	{
	  return(-1);
	}
    }
  else if ( length > 1 )
    {
      if ( val < start[0] )
	{
	  return(-1);
	}
      else if ( val > start[length-1] )
	{
	  return(-1);
	}
      else if ( val == start[0] )
	{
	  return(0);
	}
      else if ( val == start[length-1] )
	{
	  return(length-1);
	}
      else
	{
	  half = length/2;
	  if ( val < start[half] )
	    {
	      rtn = findex_mono(val, start, half);
	      return(rtn);
	    }
	  else if ( val > start[half] )
	    {
	      rtn = findex_mono(val, &start[half], length-half);
	      if ( rtn != -1 )
		{
		  return(half+rtn);
		}
	      else
		{
		  return(-1);
		}
	    }
	  else
	    {
	      return(half);
	    }
	}
    }  
  else
    {
      fprintf(stderr, "Bad length.\n");
    }

  return(-1);
}




/*
 * fence_post() -- return the index in an integer array where the 
 *                 provided integer value is bounded between fenceposts as
 *
 *			array[index] <= val < array[index+1]
 *
 *		   
 * Assumptions:
 *		[1] The array is monotonically increasing with index.
 *
 *		[2] If val < array[0], then -1 is returned.
 *
 *		[3] If val >= array[length-1], then -1 is returned.
 *
 * (This routine was created to quickly find that element block to which a
 *  given element belongs. 
 *
 * Created: 1997/04/03 07:28 MST pasacki@sandia.gov
 *
 * Revised: 1997/04/21 08:43 MDT pasacki@sandia.gov
 */

int
fence_post(int val,		/* the integer we seek */
	   int *array,
	   int length)
{
  int i;
  int index;
  int first_val, last_val;
  int found;
  
  double frac;

  first_val = array[0];

  last_val = array[length-1];

  if ( first_val == last_val ) return(-1);

  if ( val < first_val )  return(-1);

  if ( val >= last_val )  return(-1);

  if ( val == first_val ) return(0);

  /*
   * Verify monotonicity. Turn off for efficiency later.
   */

  for ( i=1; i<length; i++)
    {
      if ( array[i-1] > array[i] )
	{
	  sr = sprintf(err_msg, "Non monotone map where a[%d]=%d is > a[%d]=%d",
		       i-1, array[i-1], i, array[i]);
	  EH(-1, err_msg);
	}
    }

  /*
   * Linear approximation to first guess.
   */

  frac = ((double)(val - first_val))/((double)(last_val - first_val));

  index = (int)( (double)length * frac);

  index = MIN(index, length-2);

  /*
   * From this starting point, look up or down accordingly.
   */

  found = FALSE;

  while ( ! found )
    {
      if ( array[index+1] <= val ) 
	{
	  index++;
	}
      else if ( array[index] > val )
	{
	  index--;
	}
      else
	{
	  found = TRUE;
	}
    }

#ifdef DEBUG
  fprintf(stdout, "fencepost: val= %d in (a[%d]=%d, a[%d]=%d)\n",
	  val, index, array[index], index+1, array[index+1]);
#endif
  return(index);
}





/*
 * gcf() -- return the greatest common factor of 2 integers
 *
 * This particular implementation was adapted from a usenet posting by 
 * Dave Hax, Prairie Village, KS [comp.lang.c++,1997/04/18]
 */

int 
gcf ( int x,
      int y )
{
  return( (y==0)?(x):gcf(y,x%y));
}


/*
 * get_node_index() -- return the local index into the local proc node list.
 *
 * Notes:
 *		[1] Assume proc_list is divided into 3 subsequences: 
 *		    internal nodes, boundary nodes, and external nodes.
 *
 *		[2] Assume that the global node numbers within each subsequence
 *		    are ordered monotonically increasing. This means that we
 *                  can use the faster findex_mono() routine.
 *
 *		[3] If not found, then return -1.
 *
 * Created: 1997/05/13 08:41 MDT pasacki@sandia.gov
 *
 * Revised:
 */

int
get_node_index( int global_node,
		int *node_list,
		int num_internal_nodes,
		int num_boundary_nodes,
		int num_external_nodes)
{
  int n;

  n = findex_mono(global_node, node_list, num_internal_nodes);
		  
  if ( n != -1 )
    {
      return(n);
    }
  else
    {
      n = findex_mono(global_node, node_list+num_internal_nodes, 
		      num_boundary_nodes);
      if ( n != -1 )
	{
	  return(n+num_internal_nodes);
	}
      else
	{

	  /*
	   * Now, with the contiguous indeces for each processor, the
	   * external_nodes are not guaranteed to be monotone. Use
	   * a less efficient but more rigorous search.
	   */

	  n = in_list(global_node,  
		      node_list+num_internal_nodes + num_boundary_nodes, 
		      num_external_nodes);
	  /*
	  n = findex_mono(global_node, ( node_list+num_internal_nodes +
					 num_boundary_nodes ), 
			  num_external_nodes);
	  */
	  if ( n != -1 )
	    {
	      return(n + num_internal_nodes + num_boundary_nodes);
	    }
	  else
	    {
	      return(-1);
	    }
	}
    }
  /*  return(-1); */
}


/*
 * get_internal_boundary_index() -- return local index into proc node list.
 *
 * Notes:
 *		[1] Assume proc_list is divided into 3 subsequences: 
 *		    internal nodes, boundary nodes, and external nodes.
 *
 *		[2] Assume that the global node numbers within each subsequence
 *		    are ordered monotonically increasing. This means that we
 *                  can use the faster findex_mono() routine.
 *
 *		[3] If not found, then return -1.
 *
 *		[4] Only look over the internal and boundary nodes.
 *
 * Created: 1997/05/14 07:21 MDT pasacki@sandia.gov
 *
 * Revised:
 */

int
get_internal_boundary_index( int global_node,
			     int *node_list,
			     int num_internal_nodes,
			     int num_boundary_nodes)
{
  int n;

  n = findex_mono(global_node, node_list, num_internal_nodes);
		  
  if ( n != -1 )
    {
      return(n);
    }
  else
    {
      n = findex_mono(global_node, node_list+num_internal_nodes, 
		      num_boundary_nodes);
      if ( n != -1 )
	{
	  return(n+num_internal_nodes);
	}
      else
	{
	  return(-1);
	}
    }
  /*  return(-1); */
}

/*
 * proc_sort() -- sort external nodes according to assigment of processors
 *
 * Notes:
 *		[1] Using the system qsort() means that not only will the
 *		    resorted list be contiguous, but that processors with
 *		    smaller integer names will occur first.
 *
 * Created: 1998/02/21 11:22 MST pasacki@sandia.gov
 *
 * Revised:
 */

void 
proc_sort(int *node_list,
	  int len,
	  int len_assignment,
	  int *assignment)
{
  keep_len_assignment = len_assignment;
  assignment_copy = assignment;

  qsort(node_list, len, sizeof(int), proc_ident);

  return;
}


/* proc_ident() -- comparison function used by qsort. 1st by owner, 2nd by num.
 *
 * Created: 1998/02/21 11:34 MST pasacki@sandia.gov
 *
 * Revised:
 */

static int 
proc_ident(const void *arg1, 
	   const void *arg2)
{
  int *a;
  int *b;
  
  short proc_a;
  short proc_b;

  a = (int *)arg1;
  b = (int *)arg2;

  if ( *a > keep_len_assignment-1 )
    {
      fprintf(stderr, "proc_ident called with arg1 = %d (too big, > %d-1)\n",
	      *a, keep_len_assignment);
    }

  if ( *a < 0 )
    {
      fprintf(stderr, "proc_ident called with arg1 = %d (too low, < 0)\n", *a);
    }

  if ( *b > keep_len_assignment-1 )
    {
      fprintf(stderr, "proc_ident called with arg2 = %d (too big, > %d-1)\n",
	      *b, keep_len_assignment);
    }

  if ( *b < 0 )
    {
      fprintf(stderr, "proc_ident called with arg2 = %d (too low, < 0)\n", *b);
    }

  proc_a = assignment_copy[*a];
  proc_b = assignment_copy[*b];
  
  /*
   * Primarily, sort according to the integer name of the owning processor.
   * Once that's done, then according to the global node number.
   */

  if ( proc_a < proc_b ) 
    {
      return(-1);
    }
  else if ( proc_a > proc_b )
    {
      return(1);
    }
  else
    {
      if ( *a < *b )
	{
	  return(-1);
	}
      else if ( *a > *b )
	{
	  return(1);
	}
      else
	{
	  return(0);		/* should not normally happen! */
	}
    }
}

void
isort(const int length,
      int *array)
{
  if ( length < 0 )
    {
      EH(-1, "Negative length array to sort.");
    }

  if ( length == 0 )
    {
      return;
    }

  qsort((void *)array, length, sizeof(int), intcompare);

  return;
}

static int
intcompare(const void *arg1,
	   const void *arg2)
{
  int *i1 = (int *)arg1;
  int *i2 = (int *)arg2;

  if ( *i1 < *i2 ) return -1;
  if ( *i1 > *i2 ) return  1;
  return 0;
}



/* get_filename_num_procs -- determine the number of processors from filenames
 *
 * When reconstituting a monolith from many pieces, the individual pieces
 * are named "basename_1ofmany.exoII", typically. Use the shell and pattern
 * matching built into it to figure out the number of processors that
 * are involved. If more than one possibility exists, then return "-1".
 * Likewise, if any error occurs, return a -1.
 *
 * My poor understanding of RE's limit the robustness of this routine
 * to reasonable filenames. It is possible to pass it filenames for which
 * the sed commands will not suffice to extract the number of processors.
 *
 * Created: 1998/10/08 10:14 MDT pasacki@sandia.gov
 *
 * Revised:
 */
int
get_filename_num_procs(const char *basename)
{
  char string_system_command[MAX_SYSTEM_COMMAND_LENGTH];
  char fixXXXXXX[] = "/tmp/fileXXXXXX";
  FILE *s;
  int val=-1;
  strcpy(fixXXXXXX,"./fileXXXXXX");

  if( mkstemp( fixXXXXXX ) == -1 ) {
    fprintf( stderr,
        "get_filename_num_procs: temporary file could not be opened\n");
    return val;
  }

  sprintf(string_system_command, 
          "ls -1t %s.exoII.*.1 | head -1 | sed -e 's/^.*\.exoII\.//' -e 's/\.1//'  > %s",
	  basename, fixXXXXXX );

  if ( -1 == system(string_system_command) )
    {
      sr = sprintf(err_msg, "Trouble passing system\n(\n    %s\n)\n", 
		   string_system_command);
      EH(-1, err_msg);
    }

  s = fopen( fixXXXXXX, "r");
  
  if ( 1 != fscanf(s, "%d", &val) )
    {
      val = -1;
    }

  /*
   *  Ok, delete the temporary file
   */
  sprintf(string_system_command, "/bin/rm -f %s", fixXXXXXX );
  system(string_system_command);


  return val;
}

int 
get_min_val_index(const int *array,	/* in */
		  const int length,	/* in */
		  int *minimum_value, /* out */
		  int *minimum_index) /* out */
{
  int i;
  int mi;
  int mv;

  if ( length < 1 ) return -1;

  mi = 0;
  mv = array[0];

  for ( i=1; i<length; i++)
    {
      if ( array[i] < mv )
	{
	  mi = i;
	  mv = array[i];
	}
    }

  *minimum_index = mi;
  *minimum_value = mv;

  return 0;
}

int 
get_max_val_index(const int *array,	/* in */
		  const int length,	/* in */
		  int *maximum_value, /* out */
		  int *maximum_index) /* out */
{
  int i;
  int mi;
  int mv;

  if ( length < 1 ) return -1;

  mi = 0;
  mv = array[0];

  for ( i=1; i<length; i++)
    {
      if ( array[i] > mv )
	{
	  mi = i;
	  mv = array[i];
	}
    }

  *maximum_index = mi;
  *maximum_value = mv;

  return 0;
}


/*
 * is_shell_type()
 * Detects if the element type passed to the function is a
 * type of shell element that we know how to deal with here.
 * 
 * Scott A Roberts, 2010-08-24
 */
int 
is_shell_type(char *elem_type) {
  if (!strcmp( elem_type, "SHELL4") ||
      !strcmp( elem_type, "shell4") ) {
    return 1;
  } else {
    return 0;
  }
}

/* is_shell_element()
 * Finds element number in Exodus database and then
 * runs is_shell_type
 */
int 
is_shell_element(Exo_DB *exo, int e) {
  int eb = find_element_block(exo,e);
  return is_shell_type(exo->eb_elem_type[eb]);
}


