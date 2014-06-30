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

/* wr_graph_file -- write out an ASCII text file of the goma problem graph
 *
 * Notes:
 *	[1] The graph file format is described in the Chaco 2.0 manual and in
 *          the Metis 2.0 release notes.
 *
 *
 *
 *
 *
 *
 *
 * Created: 1997/05/09 10:35 MDT pasacki@sandia.gov
 *
 * Revised:
 */

#define _WR_GRAPH_FILE_C

#include <stdio.h>
#include <stdlib.h>

#include "map_names.h"
#include "std.h"
#include "eh.h"
#include "wr_graph_file.h"

void 
wr_graph_file(char *out_graph_file_name,	/* gfn - graph file name */
	      char *ifn,			/* ifn - input file name */
	      char *efn,			/* efn - .exoII file name */
	      int ne,				/* number of elements */
	      int nn,				/* number of fe nodes */
	      int total_dofs,			/* num dofs whole problem */
	      int tnnz,				/* total number of nonzeroes */
	      int tnat,				/* total num assembled terms */
	      int one_neighbor_per_line,	/* boolean for graphfile fmt */
	      int *pnn,				/* ptr to node-node stuff */
	      int dragon[],			/* mins, maxes, and scales */
	      int *var_node_names,		/* to find neighbors */
	      int *nat_contribute,		/* for vertex weights */
	      int *ccs_contribute)		/* for edge weights */
{

  int edge_weight;

  int format_flag;

  int j;

  int max_eweight;
  int max_vweight;
  int min_eweight;
  int min_vweight;

  int n;
  int neighbor;
  int num_edges;
  int num_verteces;

  int scale_eweight;
  int scale_vweight;

  int vertex_name;
  int vertex_weight;

  char  err_msg[MAX_CHAR_ERR_MSG];

  FILE *gfs;

  Spfrtn sr=0;

  /*
   * Unpack these from the dragon[] array...
   */

  num_verteces  = dragon[0];
  num_edges     = dragon[1];

  min_vweight   = dragon[2];
  scale_vweight = dragon[3];
  max_vweight   = dragon[4];

  min_eweight   = dragon[5];
  scale_eweight = dragon[6];
  max_eweight   = dragon[7];

  gfs = fopen(out_graph_file_name, "w");

  if ( gfs == NULL )
    {
      sr = sprintf(err_msg, "Trouble opening graph file \"%s\"\n", 
		   out_graph_file_name);
      EH(-1, err_msg);
    }

  format_flag   = 0;
  format_flag  += 1;		/* read edge weights */
  format_flag  += 10;		/* read vertex weights */

  if ( one_neighbor_per_line )
    {
      format_flag  += 100;	/* vertex numbers will be specified */
    }

  fprintf(gfs, "%% graph file for monolith problem\n");
  fprintf(gfs, "%%\n");
  fprintf(gfs, "%% Name of input file               = %s\n", ifn);
  fprintf(gfs, "%% Name of EXODUS II FE db file     = %s\n", efn);
  fprintf(gfs, "%%\n");
  fprintf(gfs, "%% Number of elements               = %d\n", ne);
  fprintf(gfs, "%% Number of nodes                  = %d\n", nn);
  fprintf(gfs, "%% Number of degrees of freedom     = %d\n", total_dofs);
  fprintf(gfs, "%% Number of nonzero matrix entries = %d\n", tnnz);
  fprintf(gfs, "%% Number of assembled terms        = %d\n", tnat);
  fprintf(gfs, "%%\n");
  fprintf(gfs, "%% Minimum vertex weight            = %d\n", min_vweight);
  fprintf(gfs, "%% Maximum vertex weight            = %d\n", max_vweight);
  fprintf(gfs, "%% Scale   vertex weights         1 = %d\n", scale_vweight);
  fprintf(gfs, "%%\n");
  fprintf(gfs, "%% Minimum edge weight              = %d\n", min_eweight);
  fprintf(gfs, "%% Maximum edge weight              = %d\n", max_eweight);
  fprintf(gfs, "%% Scale   edge weights           1 = %d\n", scale_eweight);
  fprintf(gfs, "%%\n");
  fprintf(gfs, "%% Format:\n");

  if ( one_neighbor_per_line )
    {
      fprintf(gfs, "%% vertex_name vertex_wt neighbor_name edge_wt\n");
    }
  else
    {
      fprintf(gfs, "%% vertex_wt neighbor_name1 edge_wt1 ...\n");
    }
  fprintf(gfs, "%%\n");
  fprintf(gfs, "%% number_verteces number_edges format\n");
  fprintf(gfs, "%%\n");
  
  fprintf(gfs, "%11d %15d %8d\n", num_verteces, num_edges, format_flag);

  /*
   * Step through node by node and write graph file line of vertex weights,
   * neighbors, and edge weights.
   */
  
  if ( one_neighbor_per_line )
    {
      for ( n=0; n<nn; n++)
	{
	  vertex_name   = n+1;
	  vertex_weight = 0;
	  
	  /*
	   * Look over every contribution to assembled terms for this
	   * node to total up the vertex weight.
	   */
	  
	  for ( j=pnn[n]; j<pnn[n+1]; j++)
	    {
	      vertex_weight += nat_contribute[j];
	    }
	  
	  /*
	   * Revisit each node-node interaction, writing out the
	   * vertex name and weight each time a distinct neighbor and 
	   * edge weight is listed.
	   */
	  
	  for ( j=pnn[n]; j<pnn[n+1]; j++)
	    {
	      neighbor    = var_node_names[j] + 1;
	      edge_weight = ccs_contribute[j];
	      if ( edge_weight != 0 && neighbor != vertex_name )
		{
		  fprintf(gfs, "%3d %3d %3d %3d\n", vertex_name, 
			  (vertex_weight/scale_vweight),
			  neighbor, (edge_weight/scale_eweight));
		}
	    }
	}
    }
  else
    {
      for ( n=0; n<nn; n++)
	{
	  vertex_name   = n + 1;
	  vertex_weight = 0;
	  
	  for ( j=pnn[n]; j<pnn[n+1]; j++)
	    {
	      vertex_weight += nat_contribute[j];
	    }
	  
	  fprintf(gfs, "%3d ", vertex_weight/scale_vweight);
	  
	  for ( j=pnn[n]; j<pnn[n+1]; j++)
	    {
	      neighbor = var_node_names[j] + 1;
	      edge_weight = ccs_contribute[j];
	      if ( edge_weight != 0 && neighbor != n+1 )
		{
		  fprintf(gfs, "%3d %3d ", neighbor, edge_weight/scale_eweight);
		}
	    }
	  
	  fprintf(gfs, "\n");
	}
    }

  fclose(gfs);

  if ( sr < 0 ) exit(2);

  return;
}
