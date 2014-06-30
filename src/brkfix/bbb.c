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

/* bbb -- build_big_bones, construct monolith skeleton from polylith dpi
 *
 * For symmetry with the other chunks, lets only allocate the space and 
 * fill in what is known unambiguously about the global problem.
 *
 * The largest bulk of information will be incorporated in a 2nd phase
 * as we traverse the polyliths. This may be done more than once, depending
 * on the number of timeplanes of data that must be reconstituted.
 *
 * While the first polylith is used in this routine for some purposes, that
 * part of the information that we also want from OTHER polyliths is deferred
 * so that the 1st polylith may be treated just like any other during the
 * 2nd phase.
 *
 *
 * Created: 1998/08/07 09:29 MDT pasacki@sandia.gov
 *
 * Revised: 1998/08/11 06:49 MDT pasacki@sandia.gov
 */

#define _BBB_C

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>

#include <sys/types.h>
#include <sys/stat.h>

#ifdef HAVE_FCNTL_H
#include <fcntl.h>
#endif

#ifdef STDC_HEADERS
#include <stdlib.h>
#endif

#include <string.h>

#include "std.h"		/* useful general stuff */
#include "eh.h"			/* error handling */
#include "aalloc.h"		/* multi-dim array allocation */
#include "exo_struct.h"		/* some definitions for EXODUS II */
#include "dpi.h"		/* distributed processing information */
#include "nodesc.h"		/* node descriptions */
#include "bbb.h"

static const int sc  = sizeof(char);
static const int si  = sizeof(int);
static const int sd  = sizeof(dbl);
static const int spc = sizeof(char *);
static const int spi = sizeof(int *);
static const int spd = sizeof(dbl *);

static Spfrtn sr=0;		/* sprintf() return type */
static char err_msg[MAX_CHAR_ERR_MSG];

/*
 * Function definitions.
 */

void
build_big_bones(Exo_DB *p,	/* EXODUS info from representative polylith */
		Dpi    *d,	/* distributed processing info from polylith */
		Exo_DB *m)	/* ptr to monolithic EXODUS II database */
{
  int i;
  int j;
  int len;
  int len_conn_i;		/* Length of connectivity data for a block */

  /*
   * Check for reasonable state of polylith...
   */

  if ( p->state == EXODB_STATE_GRND )
    {
      EH(-1, "Polylith only in ground state.");
    }

  m->state = EXODB_STATE_GRND;

  /*
   * All the ex_open() type of information...
   */

  m->path  = (char *) smalloc(FILENAME_MAX_ACK*sizeof(char));

  for ( i=0; i<FILENAME_MAX_ACK; i++)
    {
      m->path[i] = '\0';
    }

  /*
   * Output monolith filename may well be overriden, but fill in a likely
   * placeholder name until such time.
   */

  mononame(p->path, m->path);

  m->mode          = EX_WRITE;
  m->comp_wordsize = p->comp_wordsize;
  m->io_wordsize   = 0;		/* for now...anyway */
  m->version       = p->version;

  /*
   * Basic ex_get_init() quantitites...
   */

  m->title = (char *) smalloc(MAX_LINE_LENGTH*sc);

  strcpy(m->title, p->title);

  m->num_dim         = p->num_dim;
  m->num_nodes       = d->num_nodes_global;
  m->num_elems       = d->num_elems_global;
  m->num_elem_blocks = d->num_elem_blocks_global;
  m->num_node_sets   = d->num_node_sets_global;
  m->num_side_sets   = d->num_side_sets_global;
  
  /*
   * Miscellaneous information from ex_inquire()...
   */

  /*
   * QA Records...
   */

  m->num_qa_rec      = p->num_qa_rec;

  if ( m->num_qa_rec > 0 )
    {
      m->qa_record = (QA_Record *) smalloc(m->num_qa_rec*sizeof(QA_Record));
      for ( i=0; i<m->num_qa_rec; i++)
	{
	  for ( j=0; j<4; j++)
	    {
	      m->qa_record[i][j] = (char *) smalloc(MAX_STR_LENGTH*sc);
	      strcpy(m->qa_record[i][j], p->qa_record[i][j]);
	    }
	}
    }

  /*
   * Info Records...
   */

  m->num_info = p->num_info;

  if ( m->num_info > 0 )
    {
      m->info = (INFO_Record *) smalloc(m->num_info*sizeof(INFO_Record));
      for ( i=0; i<m->num_info; i++)
	{
	  m->info[i] = (char *) smalloc(MAX_LINE_LENGTH * sc);
	  strcpy(m->info[i], p->info[i]);
	}
    }

  m->api_version     = p->api_version;
  m->db_version      = p->db_version;

  m->ns_node_len     = d->ns_node_len_global;
  m->ns_distfact_len = d->ns_distfact_len_global;
  m->ss_elem_len     = d->ss_elem_len_global;
  m->ss_distfact_len = d->ss_distfact_len_global;

  /*
   * m->ss_node_len  = d->ss_node_len_global;
   */

  /*
   * Protective measures against rampant replications of the EXODUS II
   * "ID" property virus in old versions of GOMA...
   */

  m->eb_num_props    = MIN(p->eb_num_props, d->num_props_eb);
  m->ns_num_props    = MIN(p->ns_num_props, d->num_props_ns);
  m->ss_num_props    = MIN(p->ss_num_props, d->num_props_ss);
  
  m->num_times       = p->num_times;

  /*
   * Make space for spatial coordinates...
   */

  if ( m->num_dim > 0 )
    {
      m->x_coord = (dbl *) smalloc(m->num_nodes * sd);
      for ( i=0; i<m->num_nodes; i++)
	{
	  m->x_coord[i] = -9e12;
	}
    }
      
  if ( m->num_dim > 1 )
    {
      m->y_coord = (dbl *) smalloc(m->num_nodes * sd);
      for ( i=0; i<m->num_nodes; i++)
	{
	  m->y_coord[i] = -9e12;
	}
    }
      
  if ( m->num_dim > 2 )
    {
      m->z_coord = (dbl *) smalloc(m->num_nodes * sd);
      for ( i=0; i<m->num_nodes; i++)
	{
	  m->z_coord[i] = -9e12;
	}
    }
      
  if ( m->num_dim > 0 )
    {
      m->coord_names = (char **) smalloc(m->num_dim * spc);
      for ( i=0; i<m->num_dim; i++ )
	{
	  m->coord_names[i] = (char *) smalloc(MAX_STR_LENGTH*sc);
	  strcpy(m->coord_names[i], p->coord_names[i]);
	}
    }

  /*
   * Node map...PASS!
   */

  /*
   * Element map...PASS!
   */

  /*
   * Element order map...PASS!
   */

  /*
   * Element blocks...
   */

  if ( m->num_elem_blocks > 0 )
    {
      m->eb_id                 = (int *) smalloc(m->num_elem_blocks*si);
      m->eb_elem_type          = (char **) smalloc(m->num_elem_blocks*spc);
      m->eb_num_elems          = (int *) smalloc(m->num_elem_blocks*si);
      m->eb_num_nodes_per_elem = (int *) smalloc(m->num_elem_blocks*si);
      m->eb_num_attr           = (int *) smalloc(m->num_elem_blocks*si);
      m->eb_conn               = (int **) smalloc(m->num_elem_blocks*spi);
      m->eb_attr               = (dbl **) smalloc((m->num_elem_blocks)*spd);
      m->eb_ptr                = (int *) smalloc((m->num_elem_blocks+1)*si);
      m->eb_ptr[0]             = 0;

      for ( i = 0; i < m->num_elem_blocks; i++)
	{
	  m->eb_id[i]                 = d->eb_id_global[i];
	  m->eb_elem_type[i]          = smalloc((MAX_STR_LENGTH+1)*sc);
	  strcpy(m->eb_elem_type[i], d->eb_elem_type_global[i]);
	  m->eb_num_elems[i]          = d->eb_num_elems_global[i];
	  m->eb_num_nodes_per_elem[i] = d->eb_num_nodes_per_elem_global[i];
	  m->eb_num_attr[i]           = d->eb_num_attr_global[i];
	  m->eb_ptr[i+1]              = m->eb_ptr[i] + m->eb_num_elems[i];

	  len_conn_i                  = ( m->eb_num_elems[i] * 
					  m->eb_num_nodes_per_elem[i] );

	  if ( len_conn_i > 0 )
	    {
	      m->eb_conn[i] = (int *) smalloc(len_conn_i * si);
	      for ( j=0; j<len_conn_i; j++)
		{
		  m->eb_conn[i][j] = -1; /* Initialize as "undefined". */
		}
	    }

	  len = MAX(1, m->eb_num_attr[i] * m->eb_num_elems[i]);
	  if ( len > 0 )
	    {
	      m->eb_attr[i] = (dbl *) smalloc(len * sd);
	      for ( j=0; j<len; j++)
		{
		  m->eb_attr[i][j] = -9.0e12; /* Initialize as "undefined". */
		}
	    }
	}
    }

  /*
   * Node sets...
   */

  if ( m->num_node_sets > 0 )
    {
      m->ns_id             = (int *) smalloc(m->num_node_sets* si);
      m->ns_num_nodes      = (int *) smalloc(m->num_node_sets* si);
      m->ns_num_distfacts  = (int *) smalloc(m->num_node_sets* si);
      m->ns_node_index     = (int *) smalloc(m->num_node_sets* si);
      m->ns_distfact_index = (int *) smalloc(m->num_node_sets* si);
      
      if ( m->ns_node_len > 0 )
	{
	  m->ns_node_list = (int *) smalloc(m->ns_node_len* si);
	  for ( j=0; j<m->ns_node_len; j++)
	    {
	      m->ns_node_list[j] = -1; /* Initialize as "undefined". */
	    }
	}
      
      if ( m->ns_distfact_len > 0 )
	{
	  m->ns_distfact_list = (dbl *) smalloc(m->ns_distfact_len* sd);
	  for ( j=0; j<m->ns_distfact_len; j++)
	    {
	      m->ns_distfact_list[j] = -9e12; /* Initialize as "undefined". */
	    }
	}

      for ( i=0; i<m->num_node_sets; i++)
	{
	  m->ns_id[i]             = d->ns_id_global[i];
	  m->ns_num_nodes[i]      = d->ns_num_nodes_global[i];
	  m->ns_num_distfacts[i]  = d->ns_num_distfacts_global[i];
	  m->ns_node_index[i]     = d->ns_node_index_global[i];
	  m->ns_distfact_index[i] = d->ns_distfact_index_global[i];
	}

      /*
       * The m->ns_..._list[] arrays will be populated as the individual
       * polyliths are traversed...
       */
    }

  /*
   * Side sets...
   */

  if ( m->num_side_sets > 0 ) 
    {
      m->ss_id             = (int *) smalloc(m->num_side_sets* si);
      m->ss_num_sides      = (int *) smalloc(m->num_side_sets* si);
      m->ss_num_distfacts  = (int *) smalloc(m->num_side_sets* si);
      m->ss_elem_index     = (int *) smalloc(m->num_side_sets* si);
      m->ss_distfact_index = (int *) smalloc(m->num_side_sets* si);
      
      /*
      m->ss_node_cnt_list  = (int **) smalloc(m->num_side_sets* spi);
      m->ss_node_list      = (int **) smalloc(m->num_side_sets* spi);
      m->ss_node_side_index  = (int **) smalloc(m->num_side_sets* spi);
      */
	  
      if ( m->ss_elem_len > 0 )
	{
	  m->ss_elem_list = (int *) smalloc(m->ss_elem_len* si);
	  m->ss_side_list = (int *) smalloc(m->ss_elem_len* si);

	  for ( j=0; j<m->ss_elem_len; j++)
	    {
	      m->ss_elem_list[j] = -1; /* Initialize as "undefined". */
	      m->ss_side_list[j] = -1; /* Initialize as "undefined". */
	    }
	}

      if ( m->ss_distfact_len > 0 )
	{
	  m->ss_distfact_list = (dbl *) smalloc(m->ss_distfact_len* sd);
	  for ( j=0; j<m->ss_distfact_len; j++)
	    {
	      m->ss_distfact_list[j] = -9e12; /* Initialize as "undefined". */
	    }
	}	    
      
      for ( i=0; i<m->num_side_sets; i++)
	{
	  m->ss_id[i]             = d->ss_id_global[i];
	  m->ss_num_sides[i]      = d->ss_num_sides_global[i];
	  m->ss_num_distfacts[i]  = d->ss_num_distfacts_global[i];
	  m->ss_elem_index[i]     = d->ss_elem_index_global[i];
	  m->ss_distfact_index[i] = d->ss_distfact_index_global[i];
	}

      /*
       * Forget about the ss_node_cnt_list[], etc. for now. They're
       * not required for ex_put_concat_side_set(), so don't bother
       * reconstructing them, at least for now...
       */

      m->ss_node_list_exists = FALSE;
    }

  /*
   * Properties of node sets...
   */

  if ( m->ns_num_props > 1 ) 
    {

      m->ns_prop_name = (char **) smalloc(m->ns_num_props* spc);
      for ( i=0; i<m->ns_num_props; i++)
	{
	  m->ns_prop_name[i] = (char *) smalloc(MAX_STR_LENGTH* sc);
	  strcpy(m->ns_prop_name[i], p->ns_prop_name[i]);
	}
	  
      m->ns_prop = (int **) smalloc(m->ns_num_props* spi);
      for ( i=0; i<m->ns_num_props; i++)
	{
	  m->ns_prop[i] = (int *)smalloc(m->num_node_sets* si);
	  for ( j=0; j<m->num_node_sets; j++)
	    {
	      m->ns_prop[i][j] = d->ns_prop_global[i][j];
	    }
	}
    }
    

  /*
   * Properties of side sets...
   */

  if ( m->ss_num_props > 1 ) 
    {

      m->ss_prop_name = (char **) smalloc(m->ss_num_props* spc);
      for ( i=0; i<m->ss_num_props; i++)
	{
	  m->ss_prop_name[i] = (char *) smalloc(MAX_STR_LENGTH* sc);
	  strcpy(m->ss_prop_name[i], p->ss_prop_name[i]);
	}
      
      m->ss_prop = (int **) smalloc(m->ss_num_props* spi);
      for ( i=0; i<m->ss_num_props; i++)
	{
	  m->ss_prop[i] = (int *)smalloc(m->num_side_sets* si);
	  for ( j=0; j<m->num_side_sets; j++)
	    {
	      m->ss_prop[i][j] = d->ss_prop_global[i][j];
	    }
	}
    }


  /*
   * Properties of element blocks...
   */
  
  if ( m->eb_num_props > 1 ) 
    {
      
      m->eb_prop_name = (char **) smalloc(m->eb_num_props* spc);
      for ( i=0; i<m->eb_num_props; i++)
	{
	  m->eb_prop_name[i] = (char *) smalloc(MAX_STR_LENGTH* sc);
	  strcpy(m->eb_prop_name[i], p->eb_prop_name[i]);
	}
      
      m->eb_prop = (int **) smalloc(m->eb_num_props* spi);
      for ( i=0; i<m->eb_num_props; i++)
	{
	  m->eb_prop[i] = (int *)smalloc(m->num_elem_blocks* si);
	  for ( j=0; j<m->num_elem_blocks; j++)
	    {
	      m->eb_prop[i][j] = d->eb_prop_global[i][j];
	    }
	}
    }
      
  /*
   * Results data...how many of each kind.
   */

  m->num_glob_vars = p->num_glob_vars;
  m->num_elem_vars = p->num_elem_vars;
  m->num_node_vars = p->num_node_vars;

  /*
   * Results data...their names...
   */

  if ( m->num_glob_vars > 0 )
    {
      m->glob_var_names = (char **) smalloc(m->num_glob_vars* spc);
      for ( i=0; i<m->num_glob_vars; i++)
	{
	  m->glob_var_names[i] = (char *) smalloc(MAX_STR_LENGTH* sc);
	  strcpy(m->glob_var_names[i], p->glob_var_names[i]);
	}
      
    }

  if ( m->num_elem_vars > 0 )
    {
      m->elem_var_names = (char **) smalloc(m->num_elem_vars* spc);
      for ( i=0; i<m->num_elem_vars; i++)
	{
	  m->elem_var_names[i] = (char *) smalloc(MAX_STR_LENGTH* sc);
	  strcpy(m->elem_var_names[i], p->elem_var_names[i]);
	}
    }

  if ( m->num_node_vars > 0 )
    {
      m->node_var_names = (char **) smalloc(m->num_node_vars* spc);
      for ( i=0; i<m->num_node_vars; i++)
	{
	  m->node_var_names[i] = (char *) smalloc(MAX_STR_LENGTH* sc);
	  strcpy(m->node_var_names[i], p->node_var_names[i]);
	}
    }

  /*
   * Time values...polyliths are assumed to simulate synchronously...
   */
  
  if ( m->num_times > 0 )
    {
      m->time_vals = (dbl *) smalloc(m->num_times* sizeof(dbl));
      for ( i=0; i<m->num_times; i++)
	{
	  m->time_vals[i] = p->time_vals[i];
	}
    }

  /*
   * Element variable truth table.
   *
   * Implicit in the definition of the element variable truth table
   * is an indexing of the element blocks that are known, whether they
   * are the full set of element blocks in the global problem or whether
   * they are some subset in which this processor participates.
   *
   * For now, take the easy route out for distributed problems. That is,
   * assume that each processor will be able to know what the global
   * problem's element_var_tab will be. Store that information in
   * the Dpi...
   */

  if ( m->num_elem_vars > 0 )
    {
      m->elem_var_tab = (int *) smalloc(m->num_elem_vars *
					m->num_elem_blocks * si);
      for ( i=0; i<(m->num_elem_vars*m->num_elem_blocks); i++)
	{
	  m->elem_var_tab[i] = d->elem_var_tab_global[i];
	}
    }

  m->state |= EXODB_STATE_MESH;
  m->state |= EXODB_STATE_RES0;

  return;
}

/*
 * build_global_conn() - contribute to monolithic connectivity from a polylith
 *
 *	- Assumes polylith mesh, monolith mesh and Dpi are in good shape.
 *
 * Created: 1998/08/11 08:13 MDT pasacki@sandia.gov
 *
 * Revised:
 */

void
build_global_conn(Exo_DB *p,	/* EXODUS info from representative polylith */
		  Dpi    *d,	/* distributed processing info from polylith */
		  Exo_DB *m)	/* ptr to monolithic EXODUS II database */
{
  int e;			/* element index */
  int ebi;			/* element block index */
  int ebi_global;		/* element block index in global problem */
  int elem_global;
  int elem_global_this_block;

  int l;

  int nd;			/* node name (local) */
  int nd_global;		/* node name (global) */

  int node;
  int node_global;
  int npe;			/* number of Nodes Per Element this block */

  /*
   * A few checks to throw out nonsense...
   */

  if ( p->state <= EXODB_STATE_GRND )
    {
      EH(-1, "Polylith needs to be constructed.");
    }

  if ( m->state <= EXODB_STATE_GRND )
    {
      EH(-1, "Basic monolith needs to be constructed.");
    }

  if ( d == NULL )
    {
      EH(-1, "Need to pass in a reasonable set of distributed proc info.");
    }

  /*
   * Now comes less serious cases that ought not to happen too frequently...
   */

  if ( p->num_elems == 0 )
    {
      return;
    }

  if ( m->num_elems == 0 )
    {
      return;
    }

  for ( ebi=0; ebi<p->num_elem_blocks; ebi++)
    {
      ebi_global = d->eb_index_global[ebi];
      npe        = p->eb_num_nodes_per_elem[ebi]; /* local & global ! */
      for ( e=0; e<p->eb_num_elems[ebi]; e++)
	{
	  /*
	   * The global element number is not enough. The connectivity
	   * arrays are based on the number of elements in a block, so
	   * subtract off the number of elements in all of the preceding
	   * blocks to get the element index within this block from the
	   * global perspective.
	   */

	  elem_global            = d->elem_index_global[p->eb_ptr[ebi]+e];
	  elem_global_this_block = elem_global - m->eb_ptr[ebi_global];
	  nd_global              = npe*(elem_global_this_block);
	  nd                     = npe*e;
	  for ( l=0; l<npe; l++)
	    {
	      node                                = p->eb_conn[ebi][nd+l];
	      node_global                         = d->node_index_global[node];
	      m->eb_conn[ebi_global][nd_global+l] = node_global;
	    }
	}
    }

  return;
}

/*
 * build_global_attr() - contribute to monolithic eb attributes from a polylith
 *
 *	- Assumes polylith mesh, monolith mesh and Dpi are in good shape.
 *
 * Created: 1998/08/12 14:40 MDT pasacki@sandia.gov
 *
 * Revised:
 */

void
build_global_attr(Exo_DB *p,	/* EXODUS info from representative polylith */
		  Dpi    *d,	/* distributed processing info from polylith */
		  Exo_DB *m)	/* ptr to monolithic EXODUS II database */
{
  int e;
  int ebi;
  int ebi_global;
  int elem_global;		/* 1 ... num_elems in global problem */
  int elem_global_this_block;	/* since beginning of this block */
  int l;
  int nattr;			/* num element block attributes this block */
  int na;			/* attribute index locally */
  int na_global;		/* number attribute index globally  */

  /*
   * A few checks to throw out nonsense...
   */

  if ( p->state <= EXODB_STATE_GRND )
    {
      EH(-1, "Polylith needs to be constructed.");
    }

  if ( m->state <= EXODB_STATE_GRND )
    {
      EH(-1, "Basic monolith needs to be constructed.");
    }

  if ( d == NULL )
    {
      EH(-1, "Need to pass in a reasonable set of distributed proc info.");
    }

  /*
   * Now comes less serious cases that ought not to happen too frequently...
   */

  if ( p->num_elems == 0 )
    {
      return;
    }

  if ( m->num_elems == 0 )
    {
      return;
    }

  for ( ebi=0; ebi<p->num_elem_blocks; ebi++)
    {
      ebi_global = d->eb_index_global[ebi];
      nattr      = p->eb_num_attr[ebi];		/* same for local & global */
      for ( e=0; e<p->eb_num_elems[ebi]; e++)
	{

	  /*
	   * The global element number is not enough. The attributes
	   * arrays are based on the number of elements in a block, so
	   * subtract off the number of elements in all of the preceding
	   * blocks to get the element index within this block from the
	   * global perspective.
	   */

	  elem_global            = d->elem_index_global[e];
	  elem_global_this_block = elem_global - m->eb_ptr[ebi_global];

	  na_global              = nattr * elem_global_this_block;
	  na                     = nattr * e;
	  for ( l=0; l<nattr; l++)
	    {
	      m->eb_attr[ebi_global][l] = 
		p->eb_attr[ebi][na+l];
	    }
	}
    }

  return;
}

/*
 * build_global_coords() - contribute to monolith coordinates from a polylith
 *
 *	- Assumes polylith mesh, monolith mesh and Dpi are in good shape.
 *
 * Created: 1998/08/12 08:04 MDT pasacki@sandia.gov
 *
 * Revised:
 */

void
build_global_coords(Exo_DB *p,	/* EXODUS info from representative polylith */
		    Dpi    *d,	/* distributed processing info from polylith */
		    Exo_DB *m)	/* ptr to monolithic EXODUS II database */
{
  int dim;
  int n;
  int n_global;			/* node number in global problem */

  /*
   * A few checks to throw out nonsense...
   */

  if ( p->state <= EXODB_STATE_GRND )
    {
      EH(-1, "Polylith needs to be constructed.");
    }

  if ( m->state <= EXODB_STATE_GRND )
    {
      EH(-1, "Basic monolith needs to be constructed.");
    }

  if ( d == NULL )
    {
      EH(-1, "Need to pass in a reasonable set of distributed proc info.");
    }

  if ( m->num_dim != p->num_dim )
    {
      EH(-1, "Monolith and polylith with mismatched spatial dimensions.");
    }

  if ( m->num_dim > 0 )
    {
      if ( m->x_coord == NULL )
	{
	  EH(-1, "Memory for coordinates must be allocated.");
	}
      if ( p->x_coord == NULL )
	{
	  EH(-1, "Transcribing coordinates from a null polylith!");
	}
    }

  if ( m->num_dim > 1 )
    {
      if ( m->y_coord == NULL )
	{
	  EH(-1, "Memory for coordinates must be allocated.");
	}
      if ( p->y_coord == NULL )
	{
	  EH(-1, "Transcribing coordinates from a null polylith!");
	}
    }

  if ( m->num_dim > 2 )
    {
      if ( m->z_coord == NULL )
	{
	  EH(-1, "Memory for coordinates must be allocated.");
	}
      if ( p->z_coord == NULL )
	{
	  EH(-1, "Transcribing coordinates from a null polylith!");
	}
    }


  dim = m->num_dim;

  for ( n=0; n<p->num_nodes; n++)
    {
      n_global = d->node_index_global[n];
      if ( dim > 0 )
	{
	  m->x_coord[n_global] = p->x_coord[n];
	}
      if ( dim > 1 )
	{
	  m->y_coord[n_global] = p->y_coord[n];	  
	}
      if ( dim > 2 )
	{
	  m->z_coord[n_global] = p->z_coord[n];	  
	}
    }

  return;
}


/*
 * build_global_ns() - contribute to monolith node sets from a polylith
 *
 *	- Assumes polylith mesh, monolith mesh and Dpi are in good shape.
 *
 * Created: 1998/08/12 16:03 MDT pasacki@sandia.gov
 *
 * Revised:
 */

void
build_global_ns(Exo_DB *p,	/* EXODUS info from representative polylith */
		Dpi    *d,	/* distributed processing info from polylith */
		Exo_DB *m)	/* ptr to monolithic EXODUS II database */
{
  int i;

  /*
   * A few checks to throw out nonsense...
   */

  if ( p->state <= EXODB_STATE_GRND )
    {
      EH(-1, "Polylith needs to be constructed.");
    }

  if ( m->state <= EXODB_STATE_GRND )
    {
      EH(-1, "Basic monolith needs to be constructed.");
    }

  if ( d == NULL )
    {
      EH(-1, "Need to pass in a reasonable set of distributed proc info.");
    }

  if ( m->num_node_sets == 0 )
    {
      return;
    }

  if ( p->num_node_sets == 0 )
    {
      return;
    }

  for ( i=0; i<p->ns_node_len; i++)
    {
      m->ns_node_list[ d->ns_node_list_index_global[i] ] = 
	d->node_index_global[ p->ns_node_list[i] ];
    }

  for ( i=0; i<p->ns_distfact_len; i++)
    {
      m->ns_distfact_list[ d->ns_distfact_list_index_global[i] ] =
	p->ns_distfact_list[i];
    }

  return;
}


/*
 * build_global_ss() - contribute to monolith side sets from a polylith
 *
 *	- Assumes polylith mesh, monolith mesh and Dpi are in good shape.
 *
 * Created: 1998/08/12 17:15 MDT pasacki@sandia.gov
 *
 * Revised:
 */

void
build_global_ss(Exo_DB *p,	/* EXODUS info from representative polylith */
		Dpi    *d,	/* distributed processing info from polylith */
		Exo_DB *m)	/* ptr to monolithic EXODUS II database */
{
  int i;

  /*
   * A few checks to throw out nonsense...
   */

  if ( p->state <= EXODB_STATE_GRND )
    {
      EH(-1, "Polylith needs to be constructed.");
    }

  if ( m->state <= EXODB_STATE_GRND )
    {
      EH(-1, "Basic monolith needs to be constructed.");
    }

  if ( d == NULL )
    {
      EH(-1, "Need to pass in a reasonable set of distributed proc info.");
    }

  if ( m->num_side_sets == 0 )
    {
      return;
    }

  if ( p->num_side_sets == 0 )
    {
      return;
    }

  for ( i=0; i<p->ss_elem_len; i++)
    {
      m->ss_elem_list[ d->ss_elem_list_index_global[i] ] = 
	d->elem_index_global[ p->ss_elem_list[i] ];
      m->ss_side_list[ d->ss_elem_list_index_global[i] ] = 
	p->ss_side_list[i];	/* Sides are the same. */
    }

  for ( i=0; i<p->ss_distfact_len; i++)
    {
      m->ss_distfact_list[ d->ss_distfact_list_index_global[i] ] =
	p->ss_distfact_list[i];	/* Same, local or global. */
    }

  return;
}

/*
 * build_global_res() - contribute to monolith results from a polylith
 *
 *	- Assumes polylith mesh, monolith mesh and Dpi are in good shape.
 *
 *	- Any global variables are directly overwritten - the last one
 *        is what you get.
 *
 *	- Also, transcribe the results descriptions in ev_time_indeces[], etc.
 *
 * Created: 1998/08/14 12:21 MDT pasacki@sandia.gov
 *
 * Revised:
 */

void
build_global_res(Exo_DB *p,	/* EXODUS info from representative polylith */
		 Dpi    *d,	/* distributed processing info from polylith */
		 Exo_DB *m)	/* ptr to monolithic EXODUS II database */
{
  int b;			/* element block index - local */
  int b_global;			/* element block index - global */
  int e;
  int eg;
  int elem_local;
  int elem_global;
  int index;
  int index_global;
  int n;
  int node_global;
  int t;			/* time plane counter */
  int v;			/* variable counter */
  int m_var;
  char *m_name = 0;
  char *p_name = 0;

  /*
   * A few checks to throw out nonsense...
   */

  if ( p->state <= EXODB_STATE_GRND )
    {
      EH(-1, "Polylith needs to be constructed.");
    }

  if ( m->state <= EXODB_STATE_GRND )
    {
      EH(-1, "Basic monolith needs to be constructed.");
    }

  if ( d == NULL )
    {
      EH(-1, "Need to pass in a reasonable set of distributed proc info.");
    }

  /*
   * Global variables...
   */

  if ( p->num_glob_vars > 0 )
    {
      for ( t = 0; t < p->num_gv_time_indeces; t++)
	{
	  m->gv_time_indeces[t] = p->gv_time_indeces[t];

	  for (m_var = 0; m_var < m->num_glob_vars; m_var++) 
	    {
	      m_name = (m->glob_var_names)[m_var];
	      for ( v = 0; v < p->num_glob_vars; v++)
		{
		  p_name = (p->glob_var_names)[v];
		  if (!strcmp(m_name, p_name)) 
		    {
		      m->gv[t][m_var] = p->gv[t][v];
		    }
		}

	    }
	}
    }

  /*
   * Nodal variables...
   */

  if ( p->num_node_vars > 0 )
    {
      for ( t = 0; t < p->num_nv_time_indeces; t++)
	{
	  m->nv_time_indeces[t] = p->nv_time_indeces[t];

          for (m_var = 0; m_var < m->num_nv_indeces; m_var++) 
	    {
	      m_name = (m->node_var_names)[m_var];
	      for ( v = 0; v < p->num_nv_indeces; v++)
		{
		  p_name = (p->node_var_names)[v];
		  if (!strcmp(m_name, p_name)) 
		    {
		      
		      //m->nv_indeces[m_var] = p->nv_indeces[v];
		      m->nv_indeces[m_var] = m_var + 1;

		      /*
		       * Should really only look at internal and boundary
		       * nodes and reject the external nodes' results.
		       */
		      for ( n = 0; n < (d->num_internal_nodes+d->num_boundary_nodes); n++)
			{
			  node_global = d->node_index_global[n];
			  m->nv[t][m_var][node_global] = p->nv[t][v][n];
#ifdef DEBUG
			  fprintf(stderr, 
				  "Monolith nv[time=%d][var=%d][node=%d] = %g (poly[%d][%d][node=%d])\n",
				  t, m_var, node_global, p->nv[t][v][n], t, v, n);
#endif
			}
		    }
		}
	    }
	}
    }

  /*
   * Element variables...
   */

  if ( p->num_elem_vars > 0 )
    {
      for ( t=0; t<p->num_ev_time_indeces; t++)
	{
	  m->ev_time_indeces[t] = p->ev_time_indeces[t];
	  
	  for ( b = 0; b < p->num_elem_blocks; b++)
	    {
	      b_global = d->eb_index_global[b];
	      for ( v = 0; v < p->num_elem_vars; v++)
		{
		  index        = b * p->num_elem_vars + v;
		  index_global = b_global * p->num_elem_vars + v;

		  if ( p->elem_var_tab[index] != 0 )
		    {
		      if ( m->elem_var_tab[index_global] == 0 )
			{
			  sr = sprintf(err_msg, 
				       "Inconsistency in element variable truth tables EBID %d",
				       p->eb_id[b]);
			  EH(-1, err_msg);
			  EH(sr, err_msg);
			}

		      for ( e = 0; e < p->eb_num_elems[b]; e++)
			{
			  elem_local  = p->eb_ptr[b] + e;
			  elem_global = d->elem_index_global[elem_local];
			  eg          = elem_global - m->eb_ptr[b_global];
			  m->ev[t][index_global][eg] = p->ev[t][index][e];
			}
		    }
		}
	    }
	}
    }

  return;
}

/* mononame() -- translate filename string from parallel to monolithic version
 *
 * Synopsis:  "a_b.c" -> "a.c"
 *
 *
 * Many data file names will be unique to a given processor. Construct that
 * name for the monolith. The names will be translate like this
 *
 * in string  = "basename_43of2187.suffix"
 *
 * out string = "basename.suffix"
 *
 *
 * Notes: 1. The input string is assumed to have sufficient space allocated to
 *           contain the revised name.
 *
 *	  2. The input string is assumed to have the underscore and period
 *	     characters as flags as to what will be transformed. If not,
 *	     results will be unpredictable.
 *
 * Created: 1998/08/07 10:18 MDT pasacki@sandia.gov
 *
 * Revised:
 */

void
mononame(char *in,
	 char *out)
{
  int i;

  char in_sans_suffix[FILENAME_MAX_ACK];
  char suffix[FILENAME_MAX_ACK];
  char err_msg[1024];

  char *p;
  char *q;

  Spfrtn sr=0;

  /*
   * Initialize...
   */

  for ( i=0; i<FILENAME_MAX_ACK; i++)
    {
      in_sans_suffix[i] = '\0';
      suffix[i]         = '\0';
    }

  /*
   * Look backwards from the end for the first period you find. Save the suffix
   * and the leading chunk separately.
   */

  p = strrchr(in, '.');
  strncpy(in_sans_suffix, in, (p-in));
  p = strrchr(in_sans_suffix, '.');
  strncpy(out, in_sans_suffix, (p-in_sans_suffix));

  return;
}


