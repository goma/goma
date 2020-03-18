
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
 * rd_mesh.c
 *
 * This file contains routines for reading in the finite element mesh for
 * goma. This is done both for serial and distributed versions. Some checking
 * of the consistency is done for sidesets, nodesets, and element blocks.
 *
 * The distributed processing information is read from the same file(s), if
 * we are doing distributed processing. The "DPI" was written in netcdf format
 * to augment the basic EXODUS II information by the "brk" program. This 
 * routine makes use of similar routines as the problem decomposer, brk, in
 * reading the information.
 *
 * This routine supplants el_exoII_io.c of the old way.
 *
 */

#ifndef lint
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "exodusII.h"
#include "std.h"
#include "exo_struct.h"
#include "dpi.h"
#include "el_elm.h"    /* Must be after exodusII.h */
#include "rf_fem.h"
#include "rf_mp.h"
#include "rf_io_const.h"
#include "rf_io.h"
#include "rf_allo.h"
#include "rf_bc_const.h"
#include "rf_bc.h"
#include "mm_as_structs.h"
#include "mm_as.h"
#include "mm_eh.h"
#include "dp_utils.h"
#include "el_elm_info.h"
#include "exo_conn.h"
#include "mm_elem_block_structs.h"
#include "mm_mp.h"
#include "mm_mp_structs.h"
#include "rd_dpi.h"
#include "rd_exo.h"
#include "rd_mesh.h"

#define GOMA_RD_MESH_C

/*
 * Variables defined here that are typically declared extern in el_geom.h
 *
 * These are largely legacy variables and have been deprecated in favor of
 * the exo->.... usage.
 */

int Num_Dim = 0;
int Num_Node = 0;
int Num_Elem = 0;
int Max_NP_Elem = 0;

int Num_Internal_Elems = 0;

int *GNodes = NULL;

int *GElems = NULL;

int *Proc_Connect_Ptr = NULL;
int Proc_Num_Elem_Blk;
int *Proc_Num_Elem_In_Blk = NULL;
int *Proc_Elem_Blk_Ids = NULL;
int *Proc_Elem_Blk_Types = NULL;
int *Proc_Nodes_Per_Elem = NULL;
int *Proc_Num_Attr = NULL;
int *Proc_Elem_Connect = NULL;
int Proc_Num_Node_Sets = 0;
int Proc_NS_List_Length = 0;
int *Proc_NS_Ids = NULL;
int *Proc_NS_Count = NULL;
int *Proc_NS_Pointers = NULL;
int *Proc_NS_List = NULL;
double *Proc_NS_Dist_Fact = NULL;
int Proc_Num_Side_Sets = 0;
int Proc_SS_Elem_List_Length = 0;
int Proc_SS_Node_List_Length = 0;
int *Proc_SS_Ids = NULL;
int *Proc_SS_Elem_Count = NULL;
int *Proc_SS_Node_Count = NULL;
int *Proc_SS_Elem_Pointers = NULL;
int *Proc_SS_Node_Pointers = NULL;
int *Proc_SS_Elem_List = NULL;
double *Proc_SS_Dist_Fact = NULL;

int *ss_to_blks[MAX_MAT_PER_SS+1] = {NULL};

double **Coor = NULL;

int Num_Internal_Nodes = 0;
int Num_Border_Nodes = 0;
int Num_External_Nodes = 0;

int *Matilda = NULL;
/*
 * Definitions for rf_bc_const.h
 */
int *SS_Internal_Boundary = NULL;

#ifdef HAVE_BRK
extern int _brk_
( int,
		 char **,
		 char ** );
#endif

/**************************************************************************/
/**************************************************************************/
/***************************************************************************/

/*
 * read_mesh_exoII() -- read in the finite element mesh
 *
 * Description:
 *
 * This routine reads in the finite element mesh used for solving the problem.
 * The EXODUS II API library is used extensively. It is built on top of
 * the netCDF data format specification. This routine has been completely
 * rewritten - any resemblence to Scott Hutchinson's routine of the same name
 * is purely coincidental.
 *
 * For distributed processing problems, we assume that multiple ".exoII" files
 * exist. These files contain some global problem information that is the same
 * on every processor. That global information is read just once on processor
 * zero and broadcast everywhere. Then, on the other processors, only
 * the local description of the problem is read.
 *
 * Created: 1997/07/09 11:12 MDT pasacki@sandia.gov
 */

int 
read_mesh_exoII(Exo_DB *exo,
		Dpi    *dpi)
{
  static char yo[] = "read_mesh_exoII";

  int error;
  int i;
  int len;
  int max;
  int *arr;

  multiname(ExoFile, ProcID, Num_Proc);
  error = rd_exo(exo, ExoFile, 0, ( EXODB_ACTION_RD_INIT+
				    EXODB_ACTION_RD_MESH+
				    EXODB_ACTION_RD_RES0 ) );
  check_parallel_error("Error in reading exodus file");
  /*
   *    if an error was encountered return to
   *  main continuing with a systematic shutdown
   *       rd_exo -->  rd_mesh  -->  main
   */

  if ( error == -1 ) return(error);

  /*
   * EXODUS II node names and element names are 1-based. Transform them
   * to 0-based names for convenience in goma applications.
   */

  zero_base(exo);

  if ( Num_Proc == 1 )
    {
      /*
       * Construct a likely imposter for the distributed processing information
       * that will suffice in most cases when we are running in serial mode.
       */
      uni_dpi(dpi, exo);		
    }
  else
    {
      rd_dpi(dpi, ExoFile);	/* local extra info for distributed
				   processing kept here, too.
				   some of this is stored in EXODUS names, 
				   like the element number map, but some
				   information has its own netcdf name.
				   */
      check_parallel_error("Error in reading Distributed Processing Information");
    }
	

  // SS_Internal_Boundary uses the dpi values
  SS_Internal_Boundary = alloc_int_1(exo->num_side_sets, INT_NOINIT);
  for (int ss_index = 0; ss_index < exo->num_side_sets; ss_index++)
    {
      int global_ss_index = dpi->ss_index_global[ss_index];
      SS_Internal_Boundary[ss_index] = dpi->ss_internal_global[global_ss_index];
    }

  setup_old_dpi(exo, dpi);

  setup_old_exo(exo, dpi, Num_Proc);

  /*
   * Duane mentions that blindly compiling with MDE too low (like 12)
   * and running a 27 node HEX element causes mysterious NaN results.
   * Instead, the code should stop with a meaningful message. Right on,
   * Duane!
   */

  len = dpi->num_elem_blocks_global;
  arr = dpi->eb_num_nodes_per_elem_global;

  max = -1;
  for ( i=0; i<len; i++)
    {
      if ( arr[i] > max )
	{
	  max = arr[i];
	}
    }

  if ( max > MDE )
    {
      log_msg("The mesh has elements with %d nodes.", max);
      log_err("Edit \"rf_fem_const.h\" to set MDE to %d; rebuld GOMA.", max);
    }


  check_sidesets(exo, BC_Types, Num_BC, dpi);

  check_nodesets(exo, BC_Types, Num_BC, dpi);

  check_elemblocks(exo, upd->Num_Mat, pd_glob, dpi);

  /*
   *  Allocate and set up the element block to Material ID mapping array,
   *  Matilda. Also, set up an equivalent entry in the Element Block
   *  structure.
   */
  Matilda = alloc_int_1(exo->num_elem_blocks, 0);
  setup_matilda(exo, Matilda);

  /*
   * Build up nice connectivity information. Not just the big list of
   * elem->node connectivity, but also node->elem and elem->elem and even
   * our trusty node-node.
   */

  build_elem_node(exo);

  build_node_elem(exo);

  build_elem_elem(exo);

  build_node_node(exo);

  return 0;
}


/*
 * Finally, determine whether a particular side set is internal or
 * external.
 *
 * Ground Rules:
 *	[1] A side set is external if it is not internal.
 *	[2] A side set is internal if any two element/side pairs share nodes.
 *   [3] It suffices to look for duplicate nodes only among the
 *	    element/side pairs in the same side set. I.e., an "internal" sideset
 *       is not supposed to have two side sets associated with it.
 *   [4] Note that side sets that are partially internal and partially
 *	    external are not accounted for in this scheme. Unless you rewrite
 *	    the code, you'll need to break those hybrids into pieces.
 */

int * find_ss_internal_boundary(Exo_DB *e)
{
  char err_msg[MAX_CHAR_ERR_MSG];
  int *ss_is_internal = alloc_int_1(e->num_side_sets, -1);
  int *first_side_node_list = alloc_int_1(MAX_NODES_PER_SIDE, -1);
  int *other_side_node_list = alloc_int_1(MAX_NODES_PER_SIDE, -1);

  for (int ss_index = 0; ss_index < e->num_side_sets; ss_index++)
    {
      /*
         * It suffices to check the first element/side pair. The nodes here
         * are cross-checked with the nodes in subsequent element/side pairs
         * in this same sideset.
         */
      int side  = 0;
      int start    = e->ss_node_side_index[ss_index][side];
      int end    = e->ss_node_side_index[ss_index][side+1];
      for (int i = 0; i < (end-start); i++)
        {
          first_side_node_list[i] = e->ss_node_list[ss_index][start+i];
        }

      /*
         * Sort the node numbers into ascending order.
         */
      if ((end-start) < 1)
        {
          EH(GOMA_ERROR, "Bad side node index listing!");
        }
      integer_sort((end-start), first_side_node_list);

      /*
         * Now look at the 2nd through last elem/sides nodegroups for any match,
         * but only if there are at least 2 sides in this sideset.
         *
         *	"Just one side?"
         *
         *	"You're external, buddy!
         */

      int num_sides   = e->ss_num_sides[ss_index];
      if (num_sides > 1)
        {
          side        = 1;
          int match_found = FALSE;
          do
            {
              int start    = e->ss_node_side_index[ss_index][side];
              int end    = e->ss_node_side_index[ss_index][side+1];
              for (int i = 0; i < (end-start); i++) {
                  other_side_node_list[i] = e->ss_node_list[ss_index][start+i];
                }
              if ((end-start) < 1) {
                  sprintf(err_msg,
                          "SS ID %d (%d sides), side_index[%d]=%d, side_index[%d]=%d",
                          e->ss_id[ss_index], e->ss_num_sides[ss_index],
                          side, start, side+1, end);
                  EH(GOMA_ERROR, err_msg);
                }
              integer_sort((end-start), other_side_node_list);
              int equal_vectors = TRUE;
              for (int i = 0; i < (end-start); i++)
                {
                  equal_vectors &= (other_side_node_list[i] == first_side_node_list[i]);
                }
              match_found = equal_vectors;
              side++;
            } while (side<num_sides && !match_found);

          if (match_found)
            {
              /*
               * Set this indicator to the SS ID, but any quantity not
               * equal to "0" would do just as well.
               */
              ss_is_internal[ss_index] = e->ss_id[ss_index];
            }
        }
    }
  free(other_side_node_list);
  free(first_side_node_list);
  return ss_is_internal;
}

/*
 * The new fangled read routines for EXODUS II data have filled a data
 * structure. This data structure is mined and the values of the traditional
 * named global variables are assigned appropriately. In some cases, like
 * arrays, the assignment of pointers saves space. In other cases, new
 * arrays are allocated.
 */

void
setup_old_exo(Exo_DB *e, Dpi *dpi, int num_proc)
{
  /* int blk_start; */
  int cur;			/* tmp counter for how many ebs touch a SS */
  /* int eb;			 element block index */
  int ebi;
  int elem;
  int i;
  int j;
  int k;
  int length;
  int lo;
  int mmps;			/* max EBs (~mats) per side set found */
  /*   int nb;			 number of blocks */
  int node;
  int nodes_1st_side;
  int nodes_this_side;
  int npe;
  int num_sides;
  int ss_index;			
  int ss_index_max;		/* index for SS touching the most EBs */
  /* int start; */
  int sum;

  int *ebl;			/* element block list */
  int *ebp;			/* ptrs based on element blocks */
  int *ssl;			/* lists of sidesets that EBs touch */
  int *ssp;			/* side set pointers */

  char Title[MAX_LINE_LENGTH+1];/* EXODUS II title                              */
  char err_msg[MAX_CHAR_ERR_MSG];
  strcpy(Title,          e->title);
  
  CPU_word_size        = e->comp_wordsize;
  IO_word_size         = e->io_wordsize;

  Num_Dim              = e->num_dim;
  Num_Node             = e->num_nodes;

  Num_Elem             = e->num_elems;

  Proc_Num_Side_Sets   = e->num_side_sets;

  Num_Internal_Elems   = e->num_elems;

  /*
   * Node point coordinates.
   */

  Coor = (dbl **) smalloc(Num_Dim * sizeof(dbl *));

  if ( Num_Dim > 0 )
    {
      Coor[0] = e->x_coord;
    }

  if ( Num_Dim > 1 )
    {
      Coor[1] = e->y_coord;
    }

  if ( Num_Dim > 2 )
    {
      Coor[2] = e->z_coord;
    }

  /*
   * Element Blocks.
   */

  Proc_Num_Elem_Blk    = e->num_elem_blocks;
  Proc_Nodes_Per_Elem  = e->eb_num_nodes_per_elem;
  Proc_Elem_Blk_Ids    = e->eb_id;
  Proc_Num_Elem_In_Blk = e->eb_num_elems;

  e->eb_elem_itype     = (int *) smalloc(e->num_elem_blocks * sizeof(int));

  Proc_Elem_Blk_Types  = e->eb_elem_itype;

  /*
   * Another legacy variable....
   */

  Max_NP_Elem          = -1;
  for ( i=0; i<e->num_elem_blocks; i++)
    {
      if ( e->eb_num_nodes_per_elem[i] > Max_NP_Elem )
	{
	  Max_NP_Elem = e->eb_num_nodes_per_elem[i];
	}
    }

  /*
   * Fill this array with integers for internal use...
   */

  for ( i=0; i<e->num_elem_blocks; i++)
    {
      e->eb_elem_itype[i] = get_type(e->eb_elem_type[i], 
				     e->eb_num_nodes_per_elem[i],
				     e->eb_num_attr[i]);
    }

  /*
   * Handy auxiliary pointer variable - not one from EXODUS. For each element
   * block index, this gives the lo and hi element numbers.
   */

   if ( e->eb_ptr[Proc_Num_Elem_Blk] != Num_Internal_Elems )
     {
       EH(GOMA_ERROR, "Inconsistent element count.");
     }

   /*
    * Connectivity -- consolidate from a per element block description into
    *		      a total description for all of the elements that this
    *		      processor sees from all element blocks.
    */

   length = 0;
   for ( i=0; i<e->num_elem_blocks; i++)
     {
       length += (e->eb_num_elems[i])*(e->eb_num_nodes_per_elem[i]);
     }

   e->node_list = (int *) smalloc( length * sizeof(int));
   e->elem_ptr  = (int *) smalloc( (e->num_elems+1) * sizeof(int));

   Proc_Connect_Ptr  = e->elem_ptr;
   Proc_Elem_Connect = e->node_list;

   /*
    * Load the per element block connectivities into the big connectivity.
    */

   elem = 0;
   node = 0;
   e->elem_ptr[0] = 0;

   for ( ebi=0; ebi<e->num_elem_blocks; ebi++)
     {
       npe = e->eb_num_nodes_per_elem[ebi]; 
       k = 0;
       for ( i=0; i<e->eb_num_elems[ebi]; i++ )
	 {
	   for ( j=0; j<npe; j++ )
	     {
	       e->node_list[node] = e->eb_conn[ebi][k];
	       k++;
	       node++;
	     }
	   e->elem_ptr[elem+1] = e->elem_ptr[elem] + npe;
	   elem++;
	 }
     }  

   /*
    * Attributes. Usually there are none. If you need them, go ahead and
    * use those double values via
    *
    *		e->eb_attr[eb_index][element]
    * 
    */

   Proc_Num_Attr = e->eb_num_attr;

   
   /*
    * Node Sets.
    */

   /*
    * Scalar values...
    */

   Proc_Num_Node_Sets   = e->num_node_sets;
   Proc_NS_List_Length  = e->ns_node_len;

   /*
    * These are pointers to ints...
    */

   Proc_NS_Ids          = e->ns_id;
   Proc_NS_Count        = e->ns_num_nodes;
   Proc_NS_Pointers     = e->ns_node_index;

   Proc_NS_List         = e->ns_node_list;
   Proc_NS_Dist_Fact    = e->ns_distfact_list;

   /*
    * Side Sets.
    */

   Proc_Num_Side_Sets       = e->num_side_sets;
   Proc_SS_Elem_List_Length = e->ss_elem_len;
   Proc_SS_Node_List_Length = e->ss_node_len;

   Proc_SS_Ids              = e->ss_id;
   Proc_SS_Elem_Count       = e->ss_num_sides;

   /*
    * This variable expects a single count of nodes for each side set.
    * Real life is more complicated - a side set can span element blocks
    * containing different element types with different numbers of nodes
    * per side on them. Quite heterogeneous, eh?
    *
    * Nevertheless, as a sop to the old usage, create a facsimile for the
    * old beast, doing some checking to avoid getting eaten.
    */

   Proc_SS_Node_Count       = (int *) smalloc(Proc_Num_Side_Sets*sizeof(int));

   for ( i=0; i<e->num_side_sets; i++)
     {
       nodes_1st_side = e->ss_node_cnt_list[i][0];
       for ( j=0; j<e->ss_num_sides[i]; j++)
	 {
	   nodes_this_side = e->ss_node_cnt_list[i][j];
	   if ( nodes_this_side != nodes_1st_side )
	     {
	       sprintf(err_msg, 
		    "Whoa! SS %d has sides with varying numbers of nodes.", 
		       e->ss_id[i]);
	       EH(GOMA_ERROR, err_msg);
	     }

	 }
       /*       Proc_SS_Node_Count[i] = nodes_this_side; */
     }

   /*
    * I like this better....
    */

   for ( i=0; i<e->num_side_sets; i++)
     {
       num_sides             = e->ss_num_sides[i];
       Proc_SS_Node_Count[i] = e->ss_node_side_index[i][num_sides];
     }



   Proc_SS_Elem_Pointers = e->ss_elem_index;

   /*
    * Pointer into the node list, where each SS's list of nodes begins.
    */

   if ( e->num_side_sets > 0 )	/* Thanks, Polly! */
     {
       Proc_SS_Node_Pointers    =  ( (int *) 
				     smalloc(Proc_Num_Side_Sets*sizeof(int)));
       Proc_SS_Node_Pointers[0] = 0;
     }

   for ( i=1; i<e->num_side_sets; i++)
     {
       sum = 0;
       for ( j=0; j<e->ss_num_sides[i]; j++)
	 {
	   sum += e->ss_node_cnt_list[i][j];
	 }
       Proc_SS_Node_Pointers[i] = Proc_SS_Node_Pointers[i-1] + sum;
     }

   Proc_SS_Elem_List = e->ss_elem_list;
   Proc_SS_Dist_Fact = e->ss_distfact_list;

   /*
    * The ss_node_list is doubly indexed. The Proc_SS_Node_List is a 1D array
    * that will have the same info. Now it is obsolete.
    *
    * Old fashion usage:
    * --- ------- -----
    *
    *     for ( i=0; i<Proc_SS_Node_Count[ss_index]; i++)
    *        {
    *           node = Proc_SS_Node_List[Proc_SS_Node_Pointers[ss_index] + i];
    *        }
    *
    * New fashion usage:
    * --- ------- -----
    *
    *     for ( side=0; side<exo->ss_num_sides[ss_index]; side++)
    *        {
    *           for ( l=exo->ss_node_side_index[ss_index][side];
    *                 l<exo->ss_node_side_index[ss_index][side+1]; l++ )
    *              {
    *                 node = exo->ss_node_list[ss_index][l];
    *              }
    *        }
    */

    /* Proc_SS_Node_List = e->ss_node_list; */

   /*
    * Finally, setup the famous PRS ss_to_blks connectivity that is
    * used extensively in the BC sections...
    *
    * First, find the sideset to element block connectivity.
    *
    * Second, after checking, allocate and fill Randy's array.
    */

   sseb_conn(e, &ssp, &ebl, &ebp, &ssl);

   /*
    * Verify the most highly connected side set does not touch more than
    * MAX_MAT_PER_SS element blocks...
    */

   mmps = -1;
   ss_index_max = 0;
   for ( i=0; i<e->num_side_sets; i++)
     {
       cur = ssp[i+1] - ssp[i];
       if ( cur > mmps )
	 {
	   ss_index_max = i;
	   mmps = cur;
	 }
     }

   if ( mmps > MAX_MAT_PER_SS )
     {
       sprintf(err_msg, 
	  "SS %d hits lots of EBs. Set MAX_MAX_PER_SS to %d in rf_bc_const.h",
	       e->ss_id[ss_index_max], mmps);
       EH(GOMA_ERROR, err_msg);
     }


   for ( i=0; i<MAX_MAT_PER_SS+1; i++)
     {
       ss_to_blks[i] = (int *) smalloc(e->num_side_sets * sizeof(int *));
       for ( j=0; j< e->num_side_sets; j++)
	 {
	   ss_to_blks[i][j] = -1;
	 }
     }

   if (num_proc > 1) {
    for ( ss_index=0; ss_index<e->num_side_sets; ss_index++) 
     {
        ss_to_blks[0][ss_index] = e->ss_id[ss_index];

        int global_ss_index = dpi->ss_index_global[ss_index];
        int start = dpi->ss_block_index_global[global_ss_index];
        int end = dpi->ss_block_index_global[global_ss_index +1];
        for (int bidx = start; bidx < end; bidx++) {
          // expects 1-indexed blocks
          ss_to_blks[bidx - start + 1][ss_index] = dpi->ss_block_list_global[bidx] + 1;
        }
     }
   } else {
       
    for ( ss_index=0; ss_index<e->num_side_sets; ss_index++)
        {
        ss_to_blks[0][ss_index] = e->ss_id[ss_index];

        lo = ssp[ss_index];

        for ( j=ssp[ss_index]; j<ssp[ss_index+1]; j++)
            {
            ss_to_blks[j-lo+1][ss_index] = e->eb_id[ebl[j]];
            }
        }
   }

   /*
    * Tell us about this connectivity...
    */


   /*
    * The node_map and elem_map from EXODUS are no longer appropriated by
    * our parallel processing needs. They are being returned to whatever
    * other uses you may have for them.
    */

   /*
    * Reality checks...
    */

   if (Num_Node/Num_Proc < 1) {
     sprintf(err_msg, "Whoa! Problem with %d nodes on %d processors.",
	     Num_Node, Num_Proc);
     EH(GOMA_ERROR, err_msg);
   }

   if (Num_Elem/Num_Proc < 1) {
     sprintf(err_msg, "Whoa! Problem with %d elems on %d processors.",
	     Num_Elem, Num_Proc);
     EH(GOMA_ERROR, err_msg);
   }

  /*
   * Free malloced memory from this routine
   */
  safer_free((void **) &ssp);
  safer_free((void **) &ebl);

  return;
}
/***************************************************************************/
/***************************************************************************/
/***************************************************************************/
/*
 * check_sidesets() -- cross check SS BC's with SSs defined in overall mesh.
 *
 * Purpose:
 *	    [1]	Issue warnings if any of the sideset boundary conditions that
 *		were specified in the input deck are for sideset ID's that are
 *		not contained in the mesh.
 *
 *	    [2] Issue warnings if any of the sidesets identified in the mesh
 *		are not being used for the application of sideset boundary
 *		conditions.
 * 
 * Notes:
 *
 *		This formalizes earlier code buried in the el_exoII_io.c, 
 *		extending the ideas to distributed computing using MPI to
 *		accomplish reduction.
 *
 *		Always exit if BC applies to a nonexistent SS ID.
 *	
 *		Do not always exit if an extra SS exists in the mesh.
 *
 *		The number of boundary conditions is the same quantity on
 *		every processor and is the total number of boundary conditions
 *		for the entire problem.
 *
 *		OTOH, the number of sidesets is *different* on each processor
 *		and is some subset of the total number of sidesets for the
 *		entire problem. Fortunately, the global number of sidesets
 *		and their IDs are known to every processor.
 *
 * Created: 1997/07/30 09:52 MDT pasacki@sandia.gov
 *
 * Revised:
 */

void 
check_sidesets(Exo_DB *e,	   /* EXODUS II FE db has all mesh info (in) */
	       struct Boundary_Condition bct[],		     /* BC info (in) */
	       int nbc,		       /* number of boundary conditions (in) */
	       Dpi *d)		/* distributed processing info          (in) */
{
  int found;
  int i;
  int j;
  int nss;
  int *ssid;
  int count_unused_ss;
  int *unused_ss;

  char err_msg[MAX_CHAR_IN_INPUT]; /* needed for many routines in here. */

  /*
   * The global side set IDs are known so no communication is needed. That is,
   * even if this processor does not have possession of that particular side
   * set, as long as it is in the global list we are OK.
   */

#ifndef PARALLEL
  nss  = e->num_side_sets;
  ssid = e->ss_id;
#endif

#ifdef PARALLEL
  nss  = d->num_side_sets_global;
  ssid = d->ss_id_global;
#endif

  for ( i=0; i<nbc; i++)
    {
      if ( strcmp(bct[i].Set_Type, "SS") == 0 )
	{
	  if ( in_list(bct[i].BC_ID, 0, nss, ssid) == -1)
	    {
	      sprintf(err_msg, "BC %s on SS %d is not in mesh!",
		      bct[i].desc->name1, bct[i].BC_ID);
	        EH(GOMA_ERROR, err_msg);  	/* assume this is a fatal problem */
	    }
	}
    }

  /*
   * Does every sideset in the overall mesh have a corresponding BC?
   */

  count_unused_ss = 0;
  unused_ss = (int *) smalloc(nss * sizeof(int));
  
  for ( i=0; i<nss; i++)
    {
      found = FALSE;
      for ( j=0; j<nbc; j++)
	{
	  found |= ( ( strcmp(bct[j].Set_Type, "SS") == 0 ) &&
		     ( ssid[i] == bct[j].BC_ID ) );
	  
	}
      if ( ! found )
	{
	  unused_ss[count_unused_ss] = i;
	  count_unused_ss++;
	}
    }

  
  if ( count_unused_ss > 0 && Debug_Flag > 2 )
    {
      DPRINTF(stderr, "\nUnused side sets:");
      for ( i=0; i<count_unused_ss; i++)
	{
	  DPRINTF(stderr, " %d", ssid[unused_ss[i]]);
	}
      DPRINTF(stderr, "\n");
    }

 free(unused_ss);

  return;
}

/*
 * check_nodesets() -- verify correspondence between BC NS and NS in mesh
 *
 * Purpose:	Cross check to insure that every boundary condition applied
 *		to a node set has a corresponding node set in the overall mesh.
 *
 *		Cross check to verify that every node set in the mesh has
 *		some boundary condition that applies to it. (not critical).
 *
 * Notes:	See check_sidesets above.
 *
 * Created: 1997/07/31 08:20 MDT pasacki@sandia.gov
 *	
 * Revised:
 */

void
check_nodesets(Exo_DB *e,	   /* EXODUS II FE db has all mesh info (in) */
	       struct Boundary_Condition bct[],		     /* BC info (in) */
	       int nbc,		       /* number of boundary conditions (in) */
	       Dpi *d)		/* distributed processing info          (in) */
{
  int found;
  int i;
  int j;
  int nns;
  int *nsid;
  int count_unused_ns;
  int *unused_ns;

  char err_msg[MAX_CHAR_IN_INPUT]; /* needed for many routines in here. */

  /*
   * The global node set IDs are known so no communication is needed. That is,
   * even if this processor does not have possession of that particular node
   * set, as long as it is in the global list we are OK.
   */

#ifndef PARALLEL
  nns  = e->num_node_sets;
  nsid = e->ns_id;
#endif

#ifdef PARALLEL
  nns  = d->num_node_sets_global;
  nsid = d->ns_id_global;
#endif

  for ( i=0; i<nbc; i++)
    {
      if ( strcmp(bct[i].Set_Type, "NS") == 0 )
	{
	  if ( in_list(bct[i].BC_ID, 0, nns, nsid) == -1)
	    {
	       sprintf(err_msg, "BC %s on NS %d is not in mesh!",
		       bct[i].desc->name1, bct[i].BC_ID);
	        EH(GOMA_ERROR, err_msg);  	/* assume this is a fatal problem */
	    }
	}
    }

  /*
   * Does every nodeset in the overall mesh have a corresponding BC?
   */

  count_unused_ns = 0;
  unused_ns = (int *) smalloc(nns * sizeof(int));

  for ( i=0; i<nns; i++)
    {
      found = FALSE;
      for ( j=0; j<nbc; j++)
	{
	  found |= ( ( strcmp(bct[j].Set_Type, "NS") == 0 ) &&
		     ( nsid[i] == bct[j].BC_ID ) );
	  
	}
      if ( ! found )
	{
	  unused_ns[count_unused_ns] = i;
	  count_unused_ns++;
	}
    }
  
  if ( count_unused_ns > 0 && Debug_Flag > 2)
  {
      DPRINTF(stderr, "\nUnused node sets:");
      for ( i=0; i<count_unused_ns; i++)
	  {
		  DPRINTF(stderr, " %d", nsid[unused_ns[i]]);
	  }
      DPRINTF(stderr, "\n");
  }

  free(unused_ns);

  return;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void
check_elemblocks(Exo_DB *e,	   /* EXODUS II FE db has all mesh info (in) */
		 int nmat,	   /* number of materials               (in) */
		 struct Problem_Description *p[], /* for MAT info    (in) */
		 Dpi *d)	   /* distributed processing info       (in) */

    /*************************************************************************
     *
     *  check_elemblocks() -- verify correspondence between MATs and EBs
     *
     * Purpose:	Insure that every MAT specified in the input file
     *		corresponds to an element block that really exists in the
     *		overall mesh.
     *
     *		Check to see if every element block in the mesh has a
     *		corresponding MAT associated with it. (not critical)
     *
     * Created: 1997/07/31 09:06 MDT pasacki@sandia.gov
     *
     * Revised:
     ************************************************************************/
{
  int found, i, ebid_mat, j, neb;
  int m;			/* material index */
  int *ebid;
  MATRL_PROP_STRUCT *mp_ptr;

  char err_msg[MAX_CHAR_IN_INPUT]; /* needed for many routines in here. */

  /*
   * The global node set IDs are known so no communication is needed. That is,
   * even if this processor does not know about a particular element block
   * in the e structure, it does know about that element block in the dpi
   * structure, d.
   */

#ifndef PARALLEL
  neb  = e->num_elem_blocks;
  ebid = e->eb_id;
#else
  neb  = d->num_elem_blocks_global;
  ebid = d->eb_id_global;
#endif

  /*
   * Does every element block in the list of element blocks
   * maintained by materials have a real element block associated with it?
   */

  for (i = 0; i < nmat; i++) {
    mp_ptr = mp_glob[i];
    for (j = 0; j < mp_ptr->Num_Matrl_Elem_Blk; j++) {
      ebid_mat = mp_ptr->Matrl_Elem_Blk_Ids[j];
      /*
       * Check to see if ebid_mat, obtained from the input
       * file, is indeed in the list
       * of element blocks obtained from the exodus file
       */
      if (in_list(ebid_mat, 0, neb, ebid) == -1) {
	sprintf(err_msg, "EB ID %d for MAT %s is not in mesh!",
		ebid_mat, p[i]->MaterialName);
	EH(GOMA_ERROR, err_msg);	/* assume this is a fatal problem */
      }
    }
  }

  /*
   * Does every EB in the overall mesh have a corresponding MAT?
   */

  for ( i=0; i<neb; i++) {
    found = FALSE;
    for (m = 0; (m < nmat) && !found; m++) {
      found = eb_in_matrl(ebid[i], m);
    }
    if (!found) {
      DPRINTF(stderr, "EB %d with ID %d in the mesh is unused.\n",
	      i, ebid[i]);
    }
  }

  return;
}
/***************************************************************************/
/***************************************************************************/
/***************************************************************************/

void
setup_old_dpi(Exo_DB *e,
          Dpi    *d)
{

  Num_Internal_Nodes = d->num_internal_nodes;

  Num_Border_Nodes   = d->num_boundary_nodes;

  Num_External_Nodes = d->num_external_nodes;

  return;
}
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

int
ebID_to_ebIndex(const int ebid)

    /*********************************************************************
     *
     * ebID_to_ebIndex():
     *
     *     Provides a mapping between the element block ID and the
     *     element block index.
     *
     *  Input
     * ----------
     *  ebid = Element block id
     *
     *  Output
     * ----------
     *  return = element block index in the list of element blocks.
     *********************************************************************/
{
  int eb;
  for (eb = 0; eb < EXO_ptr->num_elem_blocks; eb++) {
    if (ebid == EXO_ptr->eb_id[eb]) {
      return eb;
    }
  }
  return -1;
}
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void
setup_matilda(Exo_DB *e,  int *matilda)

    /********************************************************************
     *
     * setup_matilda() --
     *
     *  Set up element block index to material index connectivity.
     *
     * Setup the connectivity between local element block indeces and
     * the material index. The matilda array pointed to by the last
     * argument will be filled  with the appropriate indeces to permit
     * easy determination of the material index from the element block
     * index.
     *
     * Additionally, a duplicate entry in the Element Block Structure
     * is filled in.
     *
     *  Input
     * -------
     *  e - Processors exodus database structure
     *  matilda[] = previous malloced integer array,
     *              length = number of element blocks defined on proc.
     *
     *******************************************************************/
{
  int eb;			/* element block index */
  int ebid;			/* element block identifier */
  int m, found;
  ELEM_BLK_STRUCT *eb_ptr = Element_Blocks;
  
  for (eb = 0; eb < e->num_elem_blocks; eb++) {
    ebid  = e->eb_id[eb];
    found = FALSE;
    for (m = 0; m < upd->Num_Mat && (!found); m++) {
      if ((found = eb_in_matrl(ebid, m))) {
	break;
      }
    }
    if (!found) {
      /*
       *  It's currently not an error to have an element block not
       *  be part of a material
       */
#ifdef DEBUG_IGNORE_ELEMENT_BLOCK_CAPABILITY
      fprintf(stderr," filling in -1 into matilda array for missing eb block");
      matilda[eb] = -1;
      eb_ptr->MatlProp_ptr = NULL;
#else
      EH(GOMA_ERROR, "Trouble with matilda.");
#endif
    } else {
      matilda[eb] = m;
      eb_ptr->MatlProp_ptr = mp_glob[m];
    }
    eb_ptr++;
  }
  return;
}
/**************************************************************************/
/**************************************************************************/
/**************************************************************************/

int
eb_in_matrl(const int ebid, const int mn)

    /*********************************************************************
     *
     * eb_in_matrl()
     *
     *     This functions determines if a particular element block,
     *  denoted by its element block id, ebid, is in a particular
     *  material, denoted by the Material Property index, mn.      
     *  If it is, the function returns TRUE. If not, the function returns
     *  FALSE
     *
     *  Input
     * -------
     *  ebid = element block id
     *  mn   = Material index number (0 <= mn < Num_Matrl)
     *
     *  Return
     * ---------
     *  The return is boolean, i.e., either true of false.
     *********************************************************************/
{
  int pos;
  MATRL_PROP_STRUCT *mp_ptr = mp_glob[mn];
  pos = in_list(ebid, 0,  mp_ptr->Num_Matrl_Elem_Blk,
		mp_ptr->Matrl_Elem_Blk_Ids);
  if (pos < 0) return FALSE;      
  return TRUE;
}
/**************************************************************************/
/**************************************************************************/
/**************************************************************************/

int 
find_elemblock_index(const int element,
		     const Exo_DB *exo)

    /*********************************************************************
     *
     * find_elemblock_index()
     *
     *  Find element block index given the processor element number.
     *
     * Created: 1997/08/20 05:53 MDT pasacki@sandia.gov
     *
     *  Input
     * ----------
     *   element = Local element number
     *   exo     = pointer to the exodus database structure
     *
     *  Return
     * ----------
     *   return  = The element block index, ranging from 0 to
     *              exo->num_elem_blocks - 1.
     *             An error condition is indicated by setting
     *             the return value to -1.
     ********************************************************************/    
{
  int eb    = -1;
  int found = FALSE;
  if (element < 0 || element > exo->num_elems ) {
    sprintf(Err_Msg, "element %d out of range %d <= elem < %d",
            element, 0, exo->num_elems);
    EH(GOMA_ERROR, Err_Msg);
  }
  while (eb < exo->num_elem_blocks && !found ) {
    eb++;
    found = (element >= exo->eb_ptr[eb] &&
             element <  exo->eb_ptr[eb+1] );
  }
  if (! found) eb = -1;
  return (eb);
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

int
find_mat_number(const int ielem, const Exo_DB *exo)

    /********************************************************************
     *
     * find_mat_number()
     *
     *     Given an element number, this routine returns the material
     * index.
     *
     *  Input
     * ---------
     *  ielem  = Element number (processor specific)
     *  exo    = Processor exodus structure
     *
     *  Return
     *  --------
     *  material index number
     *
     *  If an error occurs, the program exits.
     ********************************************************************/
{
  int ebi;			/* element block index */
  ebi = find_elemblock_index(ielem, exo);
  EH(ebi,
     "find_mat_number: Can not find matl - unknown element block index");
  return(Matilda[ebi]);
} /* END of find_mat_num */
/*********************************************************************/
/*********************************************************************/
/*********************************************************************/

int
map_mat_index(const int ebid)

    /*************************************************************
     *
     * map_mat_index()
     *
     *   This routine will return the material index number given
     * an element block ID number. Failure to find a match 
     * produces a program error exit.
     *
     *  Input
     * --------
     *  ebid = element block ID value (must be one or greater)
     *
     *  Return
     * --------
     *  return = material index number
     *************************************************************/
{
  int m, found = FALSE;
  for (m = 0; m < upd->Num_Mat && !found; m++) {
    if (eb_in_matrl(ebid, m)) {
      found = TRUE;
      break;
    }
  }
  if (! found ) {
      if(ebid)   {
    fprintf(stderr,
	    "P_%d: Couldn't find Element block %d in any material, proceeding anyway with trepidation\n",
	    ProcID, ebid);
            }
    //EH(GOMA_ERROR, "Trouble in map_mat_index");
    return -1;
  }
  return m;
}
/*********************************************************************/
/*********************************************************************/
/*********************************************************************/
/* 
 * sseb_conn -- side set / element block connectivity builder
 *
 * Build auxiliary arrays to be used in boundary condition applications that
 * help to quickly determine answers to Randy's tough questions:
 *
 *	(i)  Which element blocks (materials) does a given sideset touch?
 *
 *	(ii) Is the sideset "internal" or "external"?
 *
 * Strategy: Build 4 generic sideset/elementblock connectivity arrays.
 *	     Use the generic connectivities to build backward-compatible 
 *	     data structures used later by goma.
 *
 *	[i]   ssp -- side set pointers (by index) into the element block list
 *	[ii]  ebl -- list of element block indeces
 *
 *	[iii] ebp -- element block pointers (by index) into the side set list
 *	[iv]  ssl -- list of side set indeces
 *
 *	NOTE!!! ebp and ssl are not yet constructed!!!!
 *
 * Assume that an element block pointer list has been developed.
 *
 * Created: 1997/07/28 10:38 MDT pasacki@sandia.gov
 *
 * Revised: 
 */

void 
sseb_conn(Exo_DB *e,		/* see exo_struct.h for full def         (in) */
	  int **side_set_pointers, /* ptrs into eb_list                 (out) */
	  int **element_block_list, /* list of ebs for ss's             (out) */
	  int **element_block_pointers,	/* ptrs into ss_list            (out) */
	  int **side_set_list)	/* lists of ss's for eb's               (out) */
{
  int begin;			/* where to start in concatenate SS elem list */
  int current_list_size;	/* of minibuffer list */
  int i;
  int j;
  /* int k; */
  int eb_index;
  int elem;
  int *list;			/* mini buffer of ebl entries temporary */
  int len_ebl;			/* current length of ebl list */
  int max_ebl;			/* keep track of needed vector length for ebl */
  int max_list;			/* maximum length of minilist */
  int *ssp;
  int *ebl;
  
  ssp = (int *) smalloc((e->num_side_sets+1)*sizeof(int));

  ssp[0]  = 0;

  max_ebl = LIST_CHUNK_SIZE;

  len_ebl = 0;

  ebl = (int *) smalloc(max_ebl*sizeof(int));

  /*
   * Search each side set. Search every element in the side set and record
   * the element block index in a growing list of distinct eb indeces.
   */

  max_list = LIST_CHUNK_SIZE;

  list     = (int *) smalloc(max_list*sizeof(int));

  for ( i=0; i<e->num_side_sets; i++)
    {
      begin = e->ss_elem_index[i];
      current_list_size = 0;
      for ( j=0; j<e->ss_num_sides[i]; j++)
	{
	  elem     = e->ss_elem_list[begin+j];
	  eb_index = fence_post(elem, e->eb_ptr, e->num_elem_blocks+1);
	  if ( eb_index == -1 )
	    {
	      EH(GOMA_ERROR, "Could not locate element in element block collection.");
	    }

	  
	  build_list(eb_index, &list, &current_list_size, &max_list);
	}

      /*
       * Now copy this side set's little list into a growing ebl, checking
       * first to see that it's big enough.
       */
      
      if ( len_ebl + current_list_size >= max_ebl )
	{
	  max_ebl += LIST_CHUNK_SIZE;
	  ebl = (int *) realloc(ebl, max_ebl*sizeof(int));
	}

      for ( j=0; j<current_list_size; j++)
	{
	  ebl[len_ebl+j] = list[j];
	}

      len_ebl += current_list_size;

      ssp[i+1] = len_ebl;
    }

  *side_set_pointers      = ssp;
  *element_block_list     = ebl;

  free(list);

  return;
}


/*
 * build_list -- add a prospective member (int) to a list if not already there.
 */

void 
build_list(int prospective_member, 
	   int **incoming_list, 
	   int *current_size,
	   int *current_max_size)
{
  int len;
  int max;
  int *list;

  len = *current_size;
  max = *current_max_size;
  list = *incoming_list;

  /*
   * If the current list has zero size, then this element is definitely
   * new.
   */

#if 0
  if ( len == 0 )
    {
      if ( len+1 >= max )
	{
	  max += LIST_CHUNK_SIZE;
	  list = (int *)realloc(list, max*sizeof(int)); 
	}
      
      list[len] = prospective_member;

      len++;
      
    }
#endif

  if ( in_list(prospective_member, 0, len, list) == -1 )
    {

      /*
       * It's not there. Before we add it, make sure the list has enough
       * space to hold it.
       */

      
      if ( len+1 >= max )
	{
	  max += LIST_CHUNK_SIZE;
	  list = (int *)realloc(list, max*sizeof(int)); 
	}
      
      list[len] = prospective_member;

      len++;
    }

  /*
   * Reciprocate from local vars to global vars...
   */

  *current_size     = len;
  *current_max_size = max;
  *incoming_list    = list;

  return;
}

/* multiname() -- translate filename string to distributed processing version
 *
 * 
 * Description:
 *
 * Many data file names will be unique to a given processor. Construct that
 * name for this processor. The names will be translate like this
 *
 * Hmmm, pathological cases will always rear their ugly heads.
 * 
 * Old names				New names
 * --- -----				--- -----
 *
 *  file.suffix				file_2of47.suffix
 *
 *  file				file_2of47
 *
 *  .suffix				_2of47.suffix
 *
 *  file.				file_2of47.
 *
 *  file.a.suffix                       file.a_2of47.suffix
 *
 *  .a.b				.a_2of47.suffix
 *
 *  a.b.				a.b_2of47.
 *
 *  a..b				a._2of47.b
 *
 *
 * Note: The input string is assumed to have sufficient space allocated to
 *       contain the revised name.
 *
 *	 Processors, while named beginning at zero, are incremented so that
 *       the string reads 1ofn, 2ofn, ...., nofn and NOT
 *	 0ofn, 1ofn, ..., n-1ofn.
 *
 *       The key rule is to search for the last occurrence of a "." if one
 *       exists. The algorithm is SUPPOSED to be robust enough to handle
 *       the pathological cases.
 *
 * Created: 1997/07/09 13:18 MDT pasacki@sandia.gov
 * 
 * Revised: 1999/09/16 16:30 MDT pasacki@sandia.gov
 */

void
multiname(char *in_name, 
	  const int processor_name, 
	  const int number_processors)
{
  char proc_string[MAX_FNL];
  char err_msg[MAX_CHAR_IN_INPUT];

  int i;

  /*
   * Don't name something as "1of1" here.
   */

  if ( number_processors == 1 ) return;

  if ( strlen( in_name ) < 1 ) return; /* Zero length names can't be multiname */

  if ( processor_name < 0 ) 
    {
      sprintf(err_msg, "processor_name = %d ( less than zero).", 
              processor_name);
      EH(GOMA_ERROR, err_msg);
    }
  else if ( processor_name > number_processors - 1 )
    {
      sprintf(err_msg, "processor_name = %d ( too high ).", 
	      processor_name);
      EH(GOMA_ERROR, err_msg);
    }

  if ( number_processors < 1 )
    {
      sprintf(err_msg, "number_processors = %d ( less than one ).", 
              number_processors);
      EH(GOMA_ERROR, err_msg);
    }
  
  for ( i=0; i<MAX_FNL; i++)
    {
      proc_string[i] = '\0';
    }

  sprintf(proc_string, ".%d.%d", number_processors, processor_name);

  strcat(in_name, proc_string);
  
  return;
}

/* get_suffix() -- copy everything including last dot into suffix; return len
 *
 * Thus "file.abc" copies ".abc" into the suffix area, which is assumed to
 * be sufficiently large.
 *
 * If the input string has no ".", then the suffix is untouched and zero is
 * returned.
 *
 * Created: 1999/09/16 16:44 MDT pasacki@sandia.gov
 */
int 
get_suffix(char *suffix_string,
	   const char *input_string)
{
  int i;
  int len;

  /*
   * Starting from the end of the string, look back until we find a
   * period "." character.
   */

  len = strlen(input_string);

  if ( len < 1 ) exit(-1);

  i = len;

  while ( i>-1 && input_string[i] != '.' ) i--;

  if ( i > -1 )			/* found a dot */
    {
      strcpy(suffix_string, &(input_string[i]));
      return (len-i+1);
    }
  else				/* did not find a dot */
    {
      suffix_string[0] = '\0';
      return 0;
    }

  /*
   * Should not be here.
   *  return -1;
   */
}

/* get_prefix() -- copy everything up to any last dot into prefix
 *
 * Thus "file.abc" copies "file" into the suffix area, which is assumed to
 * be sufficiently large.
 *
 * If the input string has no ".", then the prefix gets the entire input
 * string and the strlen is returned.
 *
 * Created: 1999/09/16 16:44 MDT pasacki@sandia.gov
 */

int 
get_prefix(char *prefix_string,
	   const char *input_string)
{
  int i;
  int len;

  /*
   * Starting from the end of the string, look back until we find a
   * period "." character.
   */

  len = strlen(input_string);

  if ( len < 1 ) exit(-1);

  /*
   * Clean out the prefix string just for safety...
   */
  for ( i=0; i<len; i++)
    {
      prefix_string[i] = '\0';
    }

  /*
   * Search backwards for first occurrence of a dot...
   */

  i = len;
  while ( i>-1 && input_string[i] != '.' ) i--;

  if ( i > -1 )			/* found a dot */
    {
      strncpy(prefix_string, &(input_string[0]), i);
      return (i);
    }
  else				/* did not find a dot */
    {
      strncpy(prefix_string, input_string, len-1);
      return len;
    }

  /*
   * Should not be here.
   *  return -1;
   */
}



/*
 * strip_suffix() -- chop off trailing characters and the period in a string
 *
 * Description:
 *	The input string is examined from back to front until a period is
 * found or the beginning of the string is reached. A new string is allocated
 * and filled with everything but the suffix. Thus "a.b" -> "a"
 *
 * Created: 1997/07/09 16:10 MDT pasacki@sandia.gov
 */

void
strip_suffix(char *result, 
	     char *in)
{
  int i,j;
  int e;

  e = strlen(in);

  if ( e < 1 ) exit(-1);

  /*
   * Starting from the end of the string, look back until we find a
   * period "." character.
   */

  i = e;

  while ( i>0 && *(in+i) != '.' ) i--;

  for ( j=0; j<i; j++)
    {
      result[j] = *(in+j);
    }
  result[i] = '\0';

  /*
   * No suffix found? Then the whole string is the basename.
   */

  if ( i == 0 )
    {
      strcpy(result, in);
    }

  return;
}

#if 0				/* kill if not missed since 1999/09/20 */
/*
 * get_suffix() -- extract the tail substring of a string past the last period
 *
 * Created: 1997/07/09 16:14 MDT pasacki@sandia.gov
 */

void
get_suffix(char *result, 
	   char *in)
{
  int i;
  int b;
  int e;

  e = strlen(in);

  b = e;

  while ( b>0 && *(in+b) != '.' ) b--;

  if ( b == 0 )
    {
      *result = '\0';
    }
  else
    {
      for ( i=b; i<e; i++)
	{
	  result[i-b] = in[i+1];
	}
    }

  return;
}

#endif

/* Elem_Type() -- return integer element type
 *
 * Return an integer identification of the element type instead of just the
 * character type.
 *
 * Input the element number on this processor.
 *
 * Notes:
 *	    [1] This function replaces an array of a similar name that
 *		stored the same information.
 *
 *
 * Created: 1997/07/23 15:19 MDT pasacki@sandia.gov
 *
 * Revised:
 */

int 
Elem_Type(const Exo_DB *exo,
	  const int element)
{
  int eb_index;
  int type;

  type = -1;			/* default */


  /*
   * Which element block index is this in?
   */

  eb_index = fence_post(element, exo->eb_ptr, (exo->num_elem_blocks)+1);
  if (eb_index < 0)
    {
      EH(GOMA_ERROR, "Fence post does not include this element.");
    }

  if (eb_index > exo->num_elem_blocks - 1)
    {
      EH(GOMA_ERROR, "Index too high.");
    }

  type = exo->eb_elem_itype[eb_index];

  return(type);
}

/************************************************************************/
/************************************************************************/
/************************************************************************/
/*
 * integer_sort() -- sort a list of integers
 *
 *
 * This makes use of the qsort "quick sort" that is part of the standard C
 * library.
 *
 * Created: 1997/09/11 07:58 MDT pasacki@sandia.gov
 */

static int
integer_compare(const void *arg1, const void *arg2)
{
  const int *a = (const int *)(arg1);
  const int *b = (const int *)(arg2);
  if (*a > *b) return  1;
  if (*a < *b) return -1;
  return(0);
}

void 
integer_sort(int length, int *array)
{
  if (length < 1) {
      EH(GOMA_ERROR, "Negative or zero length array to sort?");
  }
  qsort(array, length, sizeof(int), integer_compare);
  return;
}
/************************************************************************/
/************************************************************************/
/************************************************************************/
