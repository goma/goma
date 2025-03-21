/************************************************************************ *
* Goma - Multiphysics finite element software                             *
* Sandia National Laboratories                                            *
*                                                                         *
* Copyright (c) 2022 Goma Developers, National Technology & Engineering   *
*               Solutions of Sandia, LLC (NTESS)                          *
*                                                                         *
* Under the terms of Contract DE-NA0003525, the U.S. Government retains   *
* certain rights in this software.                                        *
*                                                                         *
* This software is distributed under the GNU General Public License.      *
* See LICENSE file.                                                       *
\************************************************************************/

/* wr_exo -- open, write EXODUS II finite element mesh, close.
 *
 * Input:
 *	x -- a pointer to a structure containing the information from
 *	     an EXODUS IIv2 finite element database. A detailed description
 *	     of this structure may be found in "exo_struct.h"
 *
 *	verbosity -- a integer that determines the level of output
 *		from rd_exo(). If zero, then no output is produced, larger
 *		values produce more detailed reporting of rd_exo()'s progress.
 *
 * Return values:
 *	integer -- A value of zero is returned upon normal completion of the
 *		   routine. A value of -1 is returned if abnormal conditions
 *		   were encountered.
 *
 * Notes:
 *
 *	1. wr_exo() does inquiries before attempting to write stuff out.
 *		    items need to exist before they can be written.
 *
 *	2. Assume that various memory allocation has already been done.
 *
 *	3. Based on code for the more general wr_exo() routine, but
 *	   chopped down to do just the mesh.
 *
 * Created: 1997/08/04 10:02 MDT pasacki@sandia.gov
 *
 * Revised:
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/utsname.h>
#include <time.h>
#include <unistd.h> /* for getuid() */

static int has_been_called = 0;

#include "el_elm.h"
#include "el_elm_info.h"
#include "exo_struct.h"
#include "exodusII.h"
#include "md_timer.h"
#include "mm_as.h"
#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "mm_eh.h"
#include "mm_elem_block_structs.h"
#include "mm_fill_util.h"
#include "mm_mp.h"
#include "mm_mp_const.h"
#include "mm_mp_structs.h"
#include "mm_post_def.h"
#include "rf_allo.h"
#include "rf_bc_const.h"
#include "rf_fem.h"
#include "rf_fem_const.h"
#include "rf_io_const.h"
#include "rf_io_structs.h" /* for Results_Description */
#include "rf_mp.h"         /* are we serial or parallel? */
#include "std.h"
#include "wr_exo.h"

#define GOMA_WR_EXO_C

extern char **Argv; /* global shadow of argv, def'd in main.c */

static Spfrtn sr; /* sprintf() return type, whatever it is. */

/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
int wr_mesh_exo(Exo_DB *exo,    /* def'd in exo_struct.h */
                char *filename, /* where to write */
                int verbosity)  /* how much to tell while writing */
{
  int i;
  int status = 0;
  struct Exodus_Base *exo_base = exo->base_mesh;

  /*
   * This is a sad and pathetic little hack intended only for short term
   * use.  If its been longer than a year since 9/3/99 it should be replaced
   * with better code.  TABAER */

  dbl dummy = 0;

  if (verbosity > 0) {
    fprintf(stderr, "wr_mesh_exo() begins.\n");
  }

  /*
   * Mesh data is so fundamental that we'll create the file with clobber,
   * obliterating any existing file of the same name. That is, preserving
   * other data in an EXODUS II file while writing onto it new mesh information
   * is deemed too extraordinary. If mesh information is written, it causes
   * all information in the file to be superseded.
   */

  exo->io_wordsize = 8;

  exo->cmode = EX_CLOBBER;
  exo->exoid = ex_create(filename, exo->cmode, &exo->comp_wordsize, &exo->io_wordsize);
  GOMA_EH(exo->exoid, "ex_create");

  if (verbosity > 1) {
    fprintf(stderr, "ex_open/create() rtn = %d\n", exo->exoid);
  }

  if (verbosity > 2) {
    fprintf(stderr, "\texo_base->path    = \"%s\"\n", exo->path);
    fprintf(stderr, "\texo_base->mode    = %d\n", exo->mode);
    fprintf(stderr, "\texo_base->comp_ws = %d\n", exo->comp_wordsize);
    fprintf(stderr, "\texo_base->io_ws   = %d\n", exo->io_wordsize);
    fprintf(stderr, "\texo_base->version = %g\n", exo->version);
  }

  if (verbosity > 1) {
    fprintf(stderr, "ex_put_init() call...\n");
  }
  status = ex_put_init(exo->exoid, exo->title, exo_base->num_dim, exo_base->num_nodes,
                       exo_base->num_elems, exo_base->num_elem_blocks, exo_base->num_node_sets,
                       exo_base->num_side_sets);

  GOMA_EH(status, "ex_put_init");

  if (verbosity > 0) {
    fprintf(stderr, "\texo_base->title           = \"%s\"\n", exo->title);
    fprintf(stderr, "\texo_base->num_nodes       = %d\n", exo_base->num_nodes);
    fprintf(stderr, "\texo_base->num_elems       = %d\n", exo_base->num_elems);
    fprintf(stderr, "\texo_base->num_elem_blocks = %d\n", exo_base->num_elem_blocks);
    fprintf(stderr, "\texo_base->num_node_sets   = %d\n", exo_base->num_node_sets);
    fprintf(stderr, "\texo_base->num_side_sets   = %d\n", exo_base->num_side_sets);
  }

  if (verbosity > 1) {
    fprintf(stderr, "\texo_base->num_qa_rec      = %d\n", exo->num_qa_rec);
    fprintf(stderr, "\texo_base->num_info        = %d\n", exo->num_info);
  }

  if (exo->num_qa_rec > 0) {
    status = ex_put_qa(exo->exoid, exo->num_qa_rec, exo->qa_record);
    GOMA_EH(status, "ex_put_qa");
  }

  if (exo->num_info > 0) {
    status = ex_put_info(exo->exoid, exo->num_info, exo->info);
    GOMA_EH(status, "ex_put_info");
  }

  if (verbosity > 0) {
    fprintf(stderr, "ex_put_coord()...\n");
  }

  if (exo_base->num_dim < 3) {
    exo_base->z_coord = &dummy;
  }

  if (exo_base->num_dim < 2) {
    exo_base->y_coord = &dummy;
  }

  if (exo_base->num_dim < 1) {
    exo_base->x_coord = &dummy;
  }

  status = ex_put_coord(exo->exoid, exo_base->x_coord, exo_base->y_coord, exo_base->z_coord);
  GOMA_EH(status, "ex_put_coord");

  status = ex_put_coord_names(exo->exoid, exo->coord_names);
  GOMA_EH(status, "ex_get_coord_names");

  /*
   * ELEMENT BLOCKS...
   */
  if (exo_base->num_elem_blocks > 0) {
    for (i = 0; i < exo_base->num_elem_blocks; i++) {
      if (verbosity > 0) {
        fprintf(stderr, "ex_put_elem_block()...\n");
      }
      status = ex_put_block(exo->exoid, EX_ELEM_BLOCK, exo_base->eb_id[i],
                            exo_base->eb_elem_type[i], exo_base->eb_num_elems[i],
                            exo_base->eb_num_nodes_per_elem[i], 0, 0, exo_base->eb_num_attr[i]);
      GOMA_EH(status, "ex_put_blocks elem");

      if ((exo_base->eb_num_elems[i] * exo_base->eb_num_nodes_per_elem[i]) > 0) {
        status =
            ex_put_conn(exo->exoid, EX_ELEM_BLOCK, exo_base->eb_id[i], exo_base->eb_conn[i], 0, 0);
        GOMA_EH(status, "ex_put_conn elem");
      }
    }
  }

  status = ex_put_id_map(exo->exoid, EX_NODE_MAP, exo->base_mesh->node_map);
  GOMA_EH(status, "ex_put_id_map EX_NODE_MAP");

  status = ex_put_id_map(exo->exoid, EX_ELEM_MAP, exo->base_mesh->elem_map);
  GOMA_EH(status, "ex_put_id_map EX_ELEM_MAP");

  /*
   * NODE SETS...
   */

  if (exo_base->num_node_sets > 0) {
    if (verbosity > 0) {
      fprintf(stderr, "ex_put_concat_sets() node sets...\n");
    }

    ex_set_specs ns_specs;

    ns_specs.sets_ids = exo_base->ns_id;
    ns_specs.num_entries_per_set = exo_base->ns_num_nodes;
    ns_specs.num_dist_per_set = exo_base->ns_num_distfacts;
    ns_specs.sets_entry_index = exo_base->ns_node_index;
    ns_specs.sets_dist_index = exo_base->ns_distfact_index;
    ns_specs.sets_entry_list = exo_base->ns_node_list;
    ns_specs.sets_extra_list = NULL;
    ns_specs.sets_dist_fact = exo_base->ns_distfact_list;

    status = ex_put_concat_sets(exo->exoid, EX_NODE_SET, &ns_specs);
    GOMA_EH(status, "ex_put_concat_sets node_sets");
  }

  /*
   * SIDE SETS...
   */

  if (exo_base->num_side_sets > 0) {
    if (verbosity > 0) {
      fprintf(stderr, "ex_put_concat_sets() side sets...\n");
    }

    ex_set_specs ss_specs;
    ss_specs.sets_ids = exo_base->ss_id;
    ss_specs.num_entries_per_set = exo_base->ss_num_sides;
    ss_specs.num_dist_per_set = exo_base->ss_num_distfacts;
    ss_specs.sets_entry_index = exo_base->ss_elem_index;
    ss_specs.sets_dist_index = exo_base->ss_distfact_index;
    ss_specs.sets_entry_list = exo_base->ss_elem_list;
    ss_specs.sets_extra_list = exo_base->ss_side_list;
    ss_specs.sets_dist_fact = exo_base->ss_distfact_list;

    status = ex_put_concat_sets(exo->exoid, EX_SIDE_SET, &ss_specs);

    GOMA_EH(status, "ex_put_concat_sets side_sets");
  }

  /*
   * PROPERTIES...
   */

  /*
   * EXODUS II will write out the default table of one property called
   * the "ID" for NS, SS, and EBs. Unless you actually have more properties
   * you want dumped, then we'll not write these out.
   */

  /*
   * Well, the damage s done. Very old EXODUS II data sets have spuriously
   * compounded ID properties already.
   */

  /*
   * Node sets...
   */

  if (exo_base->ns_num_props > 1) {

    if (verbosity > 0) {
      fprintf(stderr, "ex_put_prop_names(nodesets)...\n");
    }
    status = ex_put_prop_names(exo->exoid, EX_NODE_SET, exo_base->ns_num_props - 1,
                               &(exo_base->ns_prop_name[1]));
    GOMA_EH(status, "ex_put_prop_names(EX_NODE_SET)");

    /*
     * the following loop begins at 1 so as avoid writing
     * the first "ID" node set property table
     * This automatically added by ex_put_prop_array
     * as the first property table written to all exodus files
     * Consequently, if we were to write the "ID" table out
     * here it would continually be replicated as the file
     * is repeatedly rewritten
     */

    for (i = 1; i < exo_base->ns_num_props; i++) {
      if (strcmp(exo_base->ns_prop_name[i], "ID") != 0) {
        status = ex_put_prop_array(exo->exoid, EX_NODE_SET, exo_base->ns_prop_name[i],
                                   exo_base->ns_prop[i]);
        GOMA_EH(status, "ex_put_prop_array(EX_NODE_SET)");
      }
    }
  }

  /*
   * Side sets...
   */

  if (exo_base->ss_num_props > 1) {

    /*
     * Only write these out if the second property is not the same ole
     * "ID" like the first one...
     */

    if (verbosity > 0) {
      fprintf(stderr, "ex_put_prop_names(sidesets)...\n");
    }
    status = ex_put_prop_names(exo->exoid, EX_SIDE_SET, exo_base->ss_num_props - 1,
                               &(exo_base->ss_prop_name[1]));
    GOMA_EH(status, "ex_get_prop_names(EX_SIDE_SET)");

    for (i = 1; i < exo_base->ss_num_props; i++) {
      if (strcmp(exo_base->ss_prop_name[i], "ID") != 0) {

        status = ex_put_prop_array(exo->exoid, EX_SIDE_SET, exo_base->ss_prop_name[i],
                                   exo_base->ss_prop[i]);
        GOMA_EH(status, "ex_put_prop_array(EX_SIDE_SET)");
      }
    }
  }
  /*
   * Element blocks...
   */

  if (exo_base->eb_num_props > 1) {

    if (verbosity > 0) {
      fprintf(stderr, "ex_put_prop_names(elemblocks)...\n");
    }

    status = ex_put_prop_names(exo->exoid, EX_ELEM_BLOCK, exo_base->eb_num_props - 1,
                               &(exo_base->eb_prop_name[1]));
    GOMA_EH(status, "ex_put_prop_names(EX_ELEM_BLOCK)");

    for (i = 1; i < exo_base->eb_num_props; i++) {
      if (strcmp(exo_base->ss_prop_name[i], "ID") != 0) {

        status = ex_put_prop_array(exo->exoid, EX_ELEM_BLOCK, exo_base->eb_prop_name[i],
                                   exo_base->eb_prop[i]);
        GOMA_EH(status, "ex_put_prop_array(EX_ELEM_BLOCK)");
      }
    }
  }

  status = ex_close(exo->exoid);
  GOMA_EH(status, "ex_close()");

  return (status);
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

static int file_is_readable(const char *filename)

/***********************************************************************
 *
 * file_is_readable()
 *
 *  Attempts to determine if a file exists and is readable by this
 *  user
 *
 *  Input
 * -------
 *    filename = Relative pathname to the file.
 *
 *  Return
 * --------
 *    True (1) = File exists and is readable
 *    false(0) = File doesn't exist or is not readable.
 *
 *  (note: fopen will set the global variable, errno, as well. Therefore,
 *         there is more specific information if needed)
 ***********************************************************************/
{
  FILE *fp;
  if (filename == NULL)
    return 0;
  fp = fopen(filename, "r");
  if (fp == NULL)
    return 0;
  (void)fclose(fp);
  return 1;
}
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

/* wr_result_prelim_exo() -- open, setup EXODUS II db for results, close
 *
 * One time initialization that writes out how many nodal variables to expect
 * in the results and what their names are. This information comes in via
 * the Results_Description.
 *
 * This routine is based on the setup_exoII_results() of long ago, but modified
 * to use the new exo_struct.
 *
 * The filename written to is "exo->path", so set it to what you need and
 * remember to set it back for the next person.
 *
 * in:
 * ---
 *	rd		pointer to a "Results_Description" structure that
 * 			describes the results that are about to be dumped.
 *			for more information on this structure, see the
 *			file "rf_io.h"
 *
 *	exo		ptr to Exodus II finite element database that has
 *			been read in. Basically, expect that mesh data is
 *			known at this point, even if the results have not
 *			yet been computed.
 * out:
 * ----
 *	(none)
 *
 * Revised: 1997/08/26 14:06 MDT pasacki@sandia.gov
 */

void wr_result_prelim_exo(struct Results_Description *rd,
                          Exo_DB *exo,
                          char *filename,
                          double ***gvec_elem) {
  int i, error;
  int filename_exists; /* boolean */

  int num_vars;             /* number of var_type variables to be written */
  char *var_names[MAX_NNV]; /* NOTE: this array must be sized to be the
                               greatest value of MAX_NNV, MAX_NEV, MAX_NGV,
                               etc. It is reused for each of these variable
                               types. Currently (8/13/98), MAX_NNV is the
                               largest value @ 100 of them all - RRL   */

  char *gvar_names[MAX_NGV];
  /* array containing num_vars variable names */

  /*
   * We don't support history variables.
   */

  if (rd->nhv > 0) {
    GOMA_EH(GOMA_ERROR, "Not prepared to write history variables.");
  }

  if (filename == NULL) {
    GOMA_EH(GOMA_ERROR, "No file specified to write EXODUS II info.");
  }

  /*
   *  Figure out whether the file exists and is readable by this
   *  user.
   */
  filename_exists = file_is_readable(filename);

  /*
   * If this is an entirely new file, then we'll take the opportunity to
   * write out 8 byte floating point data. If it's an old file, then
   * we'll query it by passing a zero and respect it's wishes.
   *
   * The I/O wordsize is really immaterial to the memory resident structure
   * and is often obsolete when reading and writing memory contents to and
   * from different EXODUS II files. Same applies for the filename.
   */
  if (filename_exists) {
    exo->io_wordsize = 0; /* query the file */
  } else {
    exo->io_wordsize = 8;
  }

  exo->cmode = EX_WRITE;

  if (filename_exists) {

    exo->exoid =
        ex_open(filename, exo->cmode, &exo->comp_wordsize, &exo->io_wordsize, &exo->version);
    GOMA_EH(exo->exoid, "ex_open");
  } else {
    exo->exoid = ex_create(filename, exo->cmode, &exo->comp_wordsize, &exo->io_wordsize);
    GOMA_EH(exo->exoid, "ex_create");
  }

  /*
   * Analysis Results
   * -------- -------
   */

  /* --------------------- Global Variables -------------------------- */
  if (rd->ngv > 0) {
    num_vars = rd->ngv;
    error = ex_put_variable_param(exo->exoid, EX_GLOBAL, num_vars);
    GOMA_EH(error, "ex_put_variable_param global");
    for (i = 0; i < rd->ngv; i++) {
      gvar_names[i] = rd->gvname[i];
    }
    error = ex_put_variable_names(exo->exoid, EX_GLOBAL, num_vars, gvar_names);
    GOMA_EH(error, "ex_put_variable_names global");
  }

  /* -------------------- Element Variables -------------------------- */
  if (rd->nev > 0) {
    num_vars = rd->nev;
    error = ex_put_variable_param(exo->exoid, EX_ELEM_BLOCK, num_vars);
    GOMA_EH(error, "ex_put_variable_param elem block");
    for (i = 0; i < rd->nev; i++) {
      var_names[i] = rd->evname[i];
    }
    error = ex_put_variable_names(exo->exoid, EX_ELEM_BLOCK, num_vars, var_names);
    GOMA_EH(error, "ex_put_variable_names elem block");

    /* Create truth table at this time - saves mucho cycles later
       Also malloc the gvec_elem final dim. Easier to do right
       when the truth table is built. */
    create_truth_table(rd, exo, gvec_elem);
  }

  /* -------------------- Nodal Variables -------------------------- */
  if (rd->nnv > 0) {
    num_vars = rd->nnv;
    error = ex_put_variable_param(exo->exoid, EX_NODAL, num_vars);
    GOMA_EH(error, "ex_put_variable_param EX_NODAL");
    for (i = 0; i < rd->nnv; i++) {
      var_names[i] = rd->nvname[i];
    }
    error = ex_put_variable_names(exo->exoid, EX_NODAL, num_vars, var_names);
    GOMA_EH(error, "ex_put_variable_names nodal");
  }
  error = ex_close(exo->exoid);
  GOMA_EH(error, "ex_close");

  return;
}
/* End of wr_result_prelim_exo() ------------------------------------------- */
/***********************************************************************/
/***********************************************************************/
/***********************************************************************/

void wr_result_prelim_exo_segregated(struct Results_Description **rd,
                                     Exo_DB *exo,
                                     char *filename,
                                     double ****gvec_elem) {
  int i, error;
  int filename_exists; /* boolean */

  int num_vars;             /* number of var_type variables to be written */
  char *var_names[MAX_NNV]; /* NOTE: this array must be sized to be the
                               greatest value of MAX_NNV, MAX_NEV, MAX_NGV,
                               etc. It is reused for each of these variable
                               types. Currently (8/13/98), MAX_NNV is the
                               largest value @ 100 of them all - RRL   */

  char *gvar_names[MAX_NGV];
  /* array containing num_vars variable names */
  int total_nev = 0;
  int total_nnv = 0;

  for (pg->imtrx = 0; pg->imtrx < upd->Total_Num_Matrices; pg->imtrx++) {
    total_nev += rd[pg->imtrx]->nev;
    total_nnv += rd[pg->imtrx]->nnv;
  }

  /*
   * We don't support history variables.
   */

  if (rd[0]->nhv > 0) {
    GOMA_EH(GOMA_ERROR, "Not prepared to write history variables.");
  }

  if (filename == NULL) {
    GOMA_EH(GOMA_ERROR, "No file specified to write EXODUS II info.");
  }
  /*
   *  Figure out whether the file exists and is readable by this
   *  user.
   */
  filename_exists = file_is_readable(filename);

  /*
   * If this is an entirely new file, then we'll take the opportunity to
   * write out 8 byte floating point data. If it's an old file, then
   * we'll query it by passing a zero and respect it's wishes.
   *
   * The I/O wordsize is really immaterial to the memory resident structure
   * and is often obsolete when reading and writing memory contents to and
   * from different EXODUS II files. Same applies for the filename.
   */
  if (filename_exists) {
    exo->io_wordsize = 0; /* query the file */
  } else {
    exo->io_wordsize = 8;
  }

  exo->cmode = EX_WRITE;

  if (filename_exists) {

    exo->exoid =
        ex_open(filename, exo->cmode, &exo->comp_wordsize, &exo->io_wordsize, &exo->version);
    GOMA_EH(exo->exoid, "ex_open");
  } else {
    exo->exoid = ex_create(filename, exo->cmode, &exo->comp_wordsize, &exo->io_wordsize);
    GOMA_EH(exo->exoid, "ex_create");
  }

  /*
   * Analysis Results
   * -------- -------
   */

  /* --------------------- Global Variables -------------------------- */
  if (rd[0]->ngv > 0) {
    num_vars = rd[0]->ngv;
    error = ex_put_variable_param(exo->exoid, EX_GLOBAL, num_vars);
    GOMA_EH(error, "ex_put_variable_param global");
    for (i = 0; i < rd[0]->ngv; i++) {
      gvar_names[i] = rd[0]->gvname[i];
    }
    error = ex_put_variable_names(exo->exoid, EX_GLOBAL, num_vars, gvar_names);
    GOMA_EH(error, "ex_put_variable_names global");
  }

  /* -------------------- Element Variables -------------------------- */
  if (total_nev > 0) {
    error = ex_put_variable_param(exo->exoid, EX_ELEM_BLOCK, total_nev);
    GOMA_EH(error, "ex_put_variable_param elem block");
    pg->imtrx = 0;
    int count = 0;
    i = 0;
    while (count < total_nev) {
      while (i >= rd[pg->imtrx]->nev && pg->imtrx < upd->Total_Num_Matrices) {
        i = 0;
        pg->imtrx++;
      }

      if (pg->imtrx >= upd->Total_Num_Matrices) {
        GOMA_EH(GOMA_ERROR, "Error counting element variables");
        return;
      }

      var_names[count] = rd[pg->imtrx]->evname[i];
      i++;
      count++;
    }
    error = ex_put_variable_names(exo->exoid, EX_ELEM_BLOCK, total_nev, var_names);
    GOMA_EH(error, "ex_put_variable_names elem block");

    /* Create truth table at this time - saves mucho cycles later
       Also malloc the gvec_elem final dim. Easier to do right
       when the truth table is built. */

    create_truth_table_segregated(rd, exo, gvec_elem);
  }

  /* -------------------- Nodal Variables -------------------------- */
  if (total_nnv > 0) {
    error = ex_put_variable_param(exo->exoid, EX_NODAL, total_nnv);
    GOMA_EH(error, "ex_put_variable_param EX_NODAL");
    pg->imtrx = 0;
    int count = 0;
    i = 0;
    while (count < total_nnv) {
      while (i >= rd[pg->imtrx]->nnv && pg->imtrx < upd->Total_Num_Matrices) {
        i = 0;
        pg->imtrx++;
      }

      if (pg->imtrx >= upd->Total_Num_Matrices) {
        GOMA_EH(GOMA_ERROR, "Error counting nodal variables");
        return;
      }

      var_names[count] = rd[pg->imtrx]->nvname[i];
      i++;
      count++;
    }
    error = ex_put_variable_names(exo->exoid, EX_NODAL, total_nnv, var_names);
    GOMA_EH(error, "ex_put_variable_names nodal");
  }

  error = ex_close(exo->exoid);
  GOMA_EH(error, "ex_close");

  return;
}

void wr_nodal_result_exo(Exo_DB *exo,
                         char *filename,
                         double vector[],
                         int variable_index,
                         int time_step,
                         double time_value)

/*****************************************************************
 * write_nodal_result_exo()
 *     -- open/write/close EXODUS II db for 1 nodal var at one
 *        time step.
 *
 * The output EXODUS II database contains the original model
 * information with some minor QA and info additions, with new
 * nodal value solution data written.
 *
 ******************************************************************/
{
  char err_msg[MAX_CHAR_IN_INPUT];
  int error;
  exo->cmode = EX_WRITE;
  exo->io_wordsize = 0; /* query */
  exo->exoid = ex_open(filename, exo->cmode, &exo->comp_wordsize, &exo->io_wordsize, &exo->version);
  if (exo->exoid < 0) {
    sr = sprintf(err_msg, "ex_open() = %d on \"%s\" failure @ step %d, time = %g", exo->exoid,
                 filename, time_step, time_value);
    GOMA_EH(GOMA_ERROR, err_msg);
  }
  error = ex_put_time(exo->exoid, time_step, &time_value);
  GOMA_EH(error, "ex_put_time");
  dbl *base_vector = malloc(sizeof(dbl) * exo->base_mesh->num_nodes);
  // copy and transform vector to base_vector
  for (int i = 0; i < exo->num_nodes; i++) {
    int index = exo->ghost_node_to_base[i];
    if (index >= 0) {
      base_vector[index] = vector[i];
    }
  }
  error = ex_put_var(exo->exoid, time_step, EX_NODAL, variable_index, 1, exo->base_mesh->num_nodes,
                     base_vector);
  free(base_vector);
  GOMA_EH(error, "ex_put_var nodal");
  error = ex_close(exo->exoid);
  GOMA_EH(error, "ex_close");
  return;
}
/***********************************************************************/
/***********************************************************************/
/***********************************************************************/

void wr_elem_result_exo(Exo_DB *exo,
                        const char *filename,
                        double ***vector,
                        const int variable_index,
                        const int time_step,
                        const double time_value,
                        struct Results_Description *rd) {
  int error, i;
  double local_time_value = time_value;
  /* static char *yo = "wr_elem_result_exo"; */

  /*
   * This file must already exist.
   */

  exo->cmode = EX_WRITE;

  exo->io_wordsize = 0; /* query */
  exo->exoid = ex_open(filename, exo->cmode, &exo->comp_wordsize, &exo->io_wordsize, &exo->version);
  GOMA_EH(exo->exoid, "ex_open");

  error = ex_put_time(exo->exoid, time_step, &local_time_value);
  GOMA_EH(error, "ex_put_time");

  /* If the truth table has NOT been set up, this will be really slow... */

  for (i = 0; i < exo->num_elem_blocks; i++) {
    if (exo->elem_var_tab_exists == TRUE) {
      /* Only write out vals if this variable exists for the block */
      if (exo->elem_var_tab[i * rd->nev + variable_index] == 1 &&
          exo->base_mesh->eb_num_elems[i] > 0) {
        dbl *base_vector = malloc(sizeof(dbl) * exo->base_mesh->eb_num_elems[i]);
        for (int j = 0; j < exo->eb_num_elems[i]; j++) {
          int index = exo->eb_ghost_elem_to_base[i][j];
          if (index >= 0) {
            base_vector[index] = vector[i][variable_index][j];
          }
        }

        error = ex_put_var(exo->exoid, time_step, EX_ELEM_BLOCK, variable_index + 1, exo->eb_id[i],
                           exo->base_mesh->eb_num_elems[i], base_vector);
        free(base_vector);
        GOMA_EH(error, "ex_put_var elem");
      }
    } else {
      if (exo->base_mesh->eb_num_elems[i] > 0) {
        /* write it anyway (not really recommended from a performance viewpoint) */
        GOMA_WH(GOMA_ERROR,
                "Writing exodus element variable without truth table, contact developers");
        error = ex_put_var(exo->exoid, time_step, EX_ELEM_BLOCK, variable_index + 1, /* Convert to 1
                                                                             based for exodus */
                           exo->eb_id[i], exo->eb_num_elems[i], vector[i][variable_index]);
        GOMA_EH(error, "ex_put_var elem");
      }
    }
  }

  error = ex_close(exo->exoid);
  GOMA_EH(error, "ex_close");
  return;
}

void wr_global_result_exo(
    Exo_DB *exo, const char *filename, const int time_step, const int ngv, double u[]) {
  /*****************************************************************
   * write_global_result_exo()
   *     -- open/write/close EXODUS II db for all global values
   *
   * The output EXODUS II database contains the original model
   * information with some minor QA and info additions, with new
   * global data written.
   *
   ******************************************************************/
  int error;

  /*
   * This capability is deactivated for parallel processing.
   * brkfix doesn't support global variables, when this
   * changes this restriction should be removed. TAB 3/2002 */

  if (u == NULL)
    return; /* Do nothing if this is NULL */

  exo->cmode = EX_WRITE;
  exo->io_wordsize = 0; /* query */
  exo->exoid = ex_open(filename, exo->cmode, &exo->comp_wordsize, &exo->io_wordsize, &exo->version);
  if (exo->exoid < 0) {
    GOMA_EH(GOMA_ERROR, "wr_nodal_result_exo: could not open the output file");
  }

  error = ex_put_var(exo->exoid, time_step, EX_GLOBAL, 1, 0, ngv, u);

  GOMA_EH(error, "ex_put_var glob_vars");

  error = ex_close(exo->exoid);

  return;
}
/* End of write_exoII_results ----------------------------------------- */

/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
/*
 * create_truth_table() -- Test each proposed element against each element
 * block for existance and populate the truth table. This saves file
 * rewrites later (a netcdf phenom) when each variable is written sequentially
 * later.
 *
 * Created: 1998/08/19 08:17 MDT rrlober@sandia.gov
 */
void create_truth_table(struct Results_Description *rd, Exo_DB *exo, double ***gvec_elem) {
  char err_msg[MAX_CHAR_IN_INPUT], if_ev;
  int i, j, eb_indx, ev_indx, mat_num, error, check, iii;
  int imtrx;
  int tev, found_match, ip_total;
  ELEM_BLK_STRUCT *eb_ptr;

  static const char yo[] = "create_truth_table";

  tev = 0;
  i = 0;
  free(exo->elem_var_tab);
  exo->elem_var_tab = alloc_int_1((exo->num_elem_blocks * rd->nev), 0);
  exo->truth_table_existance_key = (int *)smalloc((V_LAST - V_FIRST) * sizeof(int));

  for (i = 0; i < V_LAST - V_FIRST; i++) {
    exo->truth_table_existance_key[i] = 0;
  }

  /* This first cycle is a test to detect which potential elem variables
     (treated as nodal in goma) exist in which block via tests below on
     whether the variables have been requested at all by the user and if
     requested, are of the appropriate integration order for conversion
     to an elem var.

     This is necessary since the array of the truth table cycles through
     the element var index fastest, and if a given element var is defined
     for one block but not another, the block in which it is undefined
     will not know the difference between a defined variable that is
     of the wrong interpolation order for this block, and a variable
     that is not defined at all for the problem. This first cycle scopes
     for these cases and sets up a temp array of all possible elem vars
     model wide that must be treated (1 or 0) in the truth table. RRL */

  if_ev = FALSE;

  for (eb_indx = 0; eb_indx < exo->num_elem_blocks; eb_indx++) {
    /* First test for all the potential elem vars from primary nodal vars
       for this block */
    if (exo->eb_num_elems[eb_indx] > 0) {
      mat_num = Matilda[eb_indx];
      if (mat_num < 0) {
        continue;
      }
      for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
        for (j = V_FIRST; j < V_LAST; j++) {
          if (pd_glob[mat_num]->v[imtrx][j] != V_NOTHING) {
            if (pd_glob[mat_num]->i[imtrx][j] == I_P0) {
              if (Num_Var_In_Type[imtrx][j] > 1) {
                fprintf(stderr,
                        "%s: Too many components in variable type for element variable %s (%s)\n",
                        yo, Exo_Var_Names[j].name2, Exo_Var_Names[j].name1);
                exit(-1);
              }
              if (exo->truth_table_existance_key[j - V_FIRST] == 0) {
                /* We just found a candidate for an element variable */
                tev += Num_Var_In_Type[imtrx][j];
                exo->truth_table_existance_key[j - V_FIRST] = 1;
              }
            }
            if (pd_glob[mat_num]->i[imtrx][j] == I_P1) {
              if (exo->truth_table_existance_key[j - V_FIRST] == 0) {
                /* We just found a candidate for an element variable */
                int dofs = getdofs(type2shape(exo->eb_elem_itype[eb_indx]), I_P1);
                tev += dofs;
                exo->truth_table_existance_key[j - V_FIRST] = dofs;
              }
            }
          }
        }
      }
      /* Now pick up all the post processing variables for this block
         - yes, for now they must
         each be listed separately and painfully */
      for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
        if (ERROR_ZZ_VEL != -1 && Num_Var_In_Type[imtrx][R_MOMENTUM1]) {
          tev++;
          if (ERROR_ZZ_VEL_ELSIZE != -1) {
            tev++;
          }
        }
        if (ERROR_ZZ_Q != -1 && Num_Var_In_Type[imtrx][R_ENERGY]) {
          tev++;
          if (ERROR_ZZ_Q_ELSIZE != -1) {
            tev++;
          }
        }
      }
      check = 0;
      for (i = 0; i < upd->Num_Mat; i++) {
        if (pd_glob[i]->MeshMotion == LAGRANGIAN || pd_glob[i]->MeshMotion == DYNAMIC_LAGRANGIAN)
          check = 1;
      }

      for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
        if (ERROR_ZZ_P != -1 && (Num_Var_In_Type[imtrx][R_MOMENTUM1] || check)) {
          tev++;
          if (ERROR_ZZ_P_ELSIZE != -1) {
            tev++;
          }
        }
      }
      /* Finally pick up all of the element-level-storage continuation
       * variables, e.g. for saturation hysteresis function
       */
      mat_num = Matilda[eb_indx];
      mp = mp_glob[mat_num];
      eb_ptr = Element_Blocks + eb_indx;
      ip_total = elem_info(NQUAD, eb_ptr->Elem_Type);
      if ((mp->PorousMediaType == POROUS_UNSATURATED ||
           mp->PorousMediaType == POROUS_SHELL_UNSATURATED ||
           mp->PorousMediaType == POROUS_TWO_PHASE) &&
          mp->SaturationModel == TANH_HYST && !if_ev) {
        for (j = 0; j < ip_total; j++) {
          if (SAT_CURVE_TYPE != -1)
            tev++; /*For Sat curve type */
          if (CAP_PRESS_SWITCH != -1)
            tev++; /*For saturation switch */
          if (SAT_QP_SWITCH != -1)
            tev++; /*for cap press switch point*/
        }
        if_ev = TRUE;
      }
    }
  }

  /* Sanity check */
  if (tev != rd->nev) {
    sr = sprintf(err_msg, "%s: Elem var count mismatch: tev(%d)<>rd->nev(%d)!?", yo, tev, rd->nev);
    GOMA_EH(GOMA_ERROR, err_msg);
    /*
    fprintf(stderr,
            "%s: Disagreement over number of element variables\n",
            yo );
    exit (-1);
    */
  }

  /* Now do the real loop and populate the truth table */
  i = 0;
  for (eb_indx = 0; eb_indx < exo->num_elem_blocks; eb_indx++) {
    if (exo->eb_num_elems[eb_indx] > 0) {
      /* First test for all the potential elem vars from primary nodal vars
         for this block */
      mat_num = Matilda[eb_indx];
      ev_indx = 0;
      for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
        for (j = V_FIRST; j < V_LAST; j++) {
          found_match = FALSE;
          if (pd_glob[mat_num]->v[imtrx][j] != V_NOTHING) {
            if (pd_glob[mat_num]->i[imtrx][j] == I_P0) {
              if (Num_Var_In_Type[imtrx][j] > 1) {
                fprintf(stderr,
                        "%s: Too many components in variable type for element variable %s (%s)\n",
                        yo, Exo_Var_Names[j].name2, Exo_Var_Names[j].name1);
                exit(-1);
              }
              /* We just found a candidate for an element variable */
              exo->elem_var_tab[i++] = 1;
              found_match = TRUE;
              ev_indx++;
              /* malloc the entry for this block by number of elems for this block
                 but - only if the variable exists for this block! (by the truth table) */
              if (has_been_called == 0) {
                /* NOTE: this final array dim is only to be malloc'd once; when a user
                   is annealing the mesh, anneal mesh calls wr_result_prelim_exo again,
                   and hence create_truth_table, which would realloc this dim of gvec_elem.
                   this test will prevent that. - RRL */
                asdv(&gvec_elem[eb_indx][ev_indx - 1], exo->eb_num_elems[eb_indx]);
              }
            }
            if (pd_glob[mat_num]->i[imtrx][j] == I_P1) {
              int dof = getdofs(type2shape(exo->eb_elem_itype[eb_indx]), I_P1);
              /* We just found a candidate for an element variable */
              for (int k = 0; k < dof; k++) {
                exo->elem_var_tab[i++] = 1;
                found_match = TRUE;
                ev_indx++;
                /* malloc the entry for this block by number of elems for this block
                   but - only if the variable exists for this block! (by the truth table) */

                if (has_been_called == 0) {
                  /* NOTE: this final array dim is only to be malloc'd once; when a user
                     is annealing the mesh, anneal mesh calls wr_result_prelim_exo again,
                     and hence create_truth_table, which would realloc this dim of gvec_elem.
                     this test will prevent that. - RRL */
                  asdv(&gvec_elem[eb_indx][ev_indx - 1], exo->eb_num_elems[eb_indx]);
                }
              }
            }
          }
          if (found_match == FALSE && exo->truth_table_existance_key[j - V_FIRST] >= 1) {
            for (int k = 0; k < exo->truth_table_existance_key[j - V_FIRST]; k++) {
              exo->elem_var_tab[i++] = 0;
              ev_indx++;
            }
          }
        }
      }
      /* Now pick up all the post processing variables for this block
         - yes, for now they must
         each be listed separately and painfully */

      for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
        if (ERROR_ZZ_VEL != -1 && Num_Var_In_Type[imtrx][R_MOMENTUM1]) {
          exo->elem_var_tab[i++] = 1;
          ev_indx++;
          /* malloc the entry for this block by number of elems for this block
             but - only if the variable exists for this block! (by the truth table) */
          if (has_been_called == 0) {
            /* NOTE: this final array dim is only to be malloc'd once; when a user
               is annealing the mesh, anneal mesh calls wr_result_prelim_exo again,
               and hence create_truth_table, which would realloc this dim of gvec_elem.
               this test will prevent that. - RRL */
            asdv(&gvec_elem[eb_indx][ev_indx - 1], exo->eb_num_elems[eb_indx]);
          }
          if (ERROR_ZZ_VEL_ELSIZE != -1) {
            exo->elem_var_tab[i++] = 1;
            ev_indx++;
            /* malloc the entry for this block by number of elems for this block
               but - only if the variable exists for this block! (by the truth table) */
            if (has_been_called == 0) {
              /* NOTE: this final array dim is only to be malloc'd once; when a user
                 is annealing the mesh, anneal mesh calls wr_result_prelim_exo again,
                 and hence create_truth_table, which would realloc this dim of gvec_elem.
                 this test will prevent that. - RRL */
              asdv(&gvec_elem[eb_indx][ev_indx - 1], exo->eb_num_elems[eb_indx]);
            }
          }
        }

        if (ERROR_ZZ_Q != -1 && Num_Var_In_Type[imtrx][R_ENERGY]) {
          exo->elem_var_tab[i++] = 1;
          ev_indx++;
          /* malloc the entry for this block by number of elems for this block
             but - only if the variable exists for this block! (by the truth table) */
          if (has_been_called == 0) {
            /* NOTE: this final array dim is only to be malloc'd once; when a user
               is annealing the mesh, anneal mesh calls wr_result_prelim_exo again,
               and hence create_truth_table, which would realloc this dim of gvec_elem.
               this test will prevent that. - RRL */
            asdv(&gvec_elem[eb_indx][ev_indx - 1], exo->eb_num_elems[eb_indx]);
          }
          if (ERROR_ZZ_Q_ELSIZE != -1) {
            exo->elem_var_tab[i++] = 1;
            ev_indx++;
            /* malloc the entry for this block by number of elems for this block
               but - only if the variable exists for this block! (by the truth table) */
            if (has_been_called == 0) {
              /* NOTE: this final array dim is only to be malloc'd once; when a user
                 is annealing the mesh, anneal mesh calls wr_result_prelim_exo again,
                 and hence create_truth_table, which would realloc this dim of gvec_elem.
                 this test will prevent that. - RRL */
              asdv(&gvec_elem[eb_indx][ev_indx - 1], exo->eb_num_elems[eb_indx]);
            }
          }
        }
      }
      check = 0;
      for (iii = 0; iii < upd->Num_Mat; iii++) {
        if (pd_glob[iii]->MeshMotion == LAGRANGIAN ||
            pd_glob[iii]->MeshMotion == DYNAMIC_LAGRANGIAN)
          check = 1;
      }

      for (imtrx = 0; imtrx < upd->Total_Num_Matrices; imtrx++) {
        if (ERROR_ZZ_P != -1 && (Num_Var_In_Type[imtrx][R_MOMENTUM1] || check)) {
          exo->elem_var_tab[i++] = 1;
          ev_indx++;
          /* malloc the entry for this block by number of elems for this block
             but - only if the variable exists for this block! (by the truth table) */
          if (has_been_called == 0) {
            /* NOTE: this final array dim is only to be malloc'd once; when a user
               is annealing the mesh, anneal mesh calls wr_result_prelim_exo again,
               and hence create_truth_table, which would realloc this dim of gvec_elem.
               this test will prevent that. - RRL */
            asdv(&gvec_elem[eb_indx][ev_indx - 1], exo->eb_num_elems[eb_indx]);
          }
          if (ERROR_ZZ_P_ELSIZE != -1) {
            exo->elem_var_tab[i++] = 1;
            ev_indx++;
            /* malloc the entry for this block by number of elems for this block
               but - only if the variable exists for this block! (by the truth table) */
            if (has_been_called == 0) {
              /* NOTE: this final array dim is only to be malloc'd once; when a user
                 is annealing the mesh, anneal mesh calls wr_result_prelim_exo again,
                 and hence create_truth_table, which would realloc this dim of gvec_elem.
                 this test will prevent that. - RRL */
              asdv(&gvec_elem[eb_indx][ev_indx - 1], exo->eb_num_elems[eb_indx]);
              has_been_called++;
            }
          }
        }
      }
      /*Now finally the saturation hysteresis variables */
      if (SAT_CURVE_TYPE != -1 || CAP_PRESS_SWITCH != -1 || SAT_QP_SWITCH != -1) {
        eb_ptr = Element_Blocks + eb_indx;
        ip_total = elem_info(NQUAD, eb_ptr->Elem_Type);
        for (j = 0; j < ip_total; j++) {
          /*Note that we will set these for all 3 var types because you
           *will never see them individually.
           */
          exo->elem_var_tab[i++] = 1;
          ev_indx++;
          if (has_been_called == 0) {
            asdv(&gvec_elem[eb_indx][ev_indx - 1], exo->eb_num_elems[eb_indx]);
          }
          exo->elem_var_tab[i++] = 1;
          ev_indx++;
          if (has_been_called == 0) {
            asdv(&gvec_elem[eb_indx][ev_indx - 1], exo->eb_num_elems[eb_indx]);
          }
          exo->elem_var_tab[i++] = 1;
          ev_indx++;
          if (has_been_called == 0) {
            asdv(&gvec_elem[eb_indx][ev_indx - 1], exo->eb_num_elems[eb_indx]);
          }
        }
      }
    } else {
      i += rd->nev;
    }
  }

  /* write out table */
  error = ex_put_truth_table(exo->exoid, EX_ELEM_BLOCK, exo->num_elem_blocks, rd->nev,
                             exo->elem_var_tab);
  GOMA_EH(error, "ex_put_truth_table EX_ELEM_BLOCK");
  /* Now set truth table exists flag */
  exo->elem_var_tab_exists = TRUE;
}

void create_truth_table_segregated(struct Results_Description **rd,
                                   Exo_DB *exo,
                                   double ****gvec_elem) {
  char err_msg[MAX_CHAR_IN_INPUT], if_ev;
  int i, j, eb_indx, ev_indx, mat_num, error, check, iii;
  int tev, found_match, ip_total;
  int total_nev;
  ELEM_BLK_STRUCT *eb_ptr;

  static const char yo[] = "create_truth_table";

  total_nev = 0;
  for (pg->imtrx = 0; pg->imtrx < upd->Total_Num_Matrices; pg->imtrx++) {
    total_nev += rd[pg->imtrx]->nev;
  }
  tev = 0;
  i = 0;
  free(exo->elem_var_tab);
  exo->elem_var_tab = (int *)smalloc((exo->num_elem_blocks * total_nev) * sizeof(int));
  exo->truth_table_existance_key = (int *)smalloc((V_LAST - V_FIRST) * sizeof(int));

  for (i = 0; i < V_LAST - V_FIRST; i++) {
    exo->truth_table_existance_key[i] = 0;
  }

  /* This first cycle is a test to detect which potential elem variables
     (treated as nodal in goma) exist in which block via tests below on
     whether the variables have been requested at all by the user and if
     requested, are of the appropriate integration order for conversion
     to an elem var.

     This is necessary since the array of the truth table cycles through
     the element var index fastest, and if a given element var is defined
     for one block but not another, the block in which it is undefined
     will not know the difference between a defined variable that is
     of the wrong interpolation order for this block, and a variable
     that is not defined at all for the problem. This first cycle scopes
     for these cases and sets up a temp array of all possible elem vars
     model wide that must be treated (1 or 0) in the truth table. RRL */

  if_ev = FALSE;

  for (eb_indx = 0; eb_indx < exo->num_elem_blocks; eb_indx++) {
    /* First test for all the potential elem vars from primary nodal vars
       for this block */
    mat_num = Matilda[eb_indx];
    if (mat_num < 0) {
      continue;
    }
    for (pg->imtrx = 0; pg->imtrx < upd->Total_Num_Matrices; pg->imtrx++) {
      for (j = V_FIRST; j < V_LAST; j++) {
        if (pd_glob[mat_num]->v[pg->imtrx][j] != V_NOTHING) {
          if (pd_glob[mat_num]->i[pg->imtrx][j] == I_P0) {
            if (Num_Var_In_Type[pg->imtrx][j] > 1) {
              fprintf(stderr,
                      "%s: Too many components in variable type for element variable %s (%s)\n", yo,
                      Exo_Var_Names[j].name2, Exo_Var_Names[j].name1);
              exit(-1);
            }
            if (exo->truth_table_existance_key[j - V_FIRST] == 0) {
              /* We just found a candidate for an element variable */
              tev += Num_Var_In_Type[pg->imtrx][j];
              exo->truth_table_existance_key[j - V_FIRST] = 1;
            }
          } else if (pd_glob[mat_num]->i[pg->imtrx][j] == I_P1) {
            if (exo->truth_table_existance_key[j - V_FIRST] == 0) {
              /* We just found a candidate for an element variable */
              int dofs = getdofs(type2shape(exo->eb_elem_itype[eb_indx]), I_P1);
              tev += dofs;
              exo->truth_table_existance_key[j - V_FIRST] = dofs;
            }
          }
        }
      }
    }
    /* Now pick up all the post processing variables for this block
       - yes, for now they must
       each be listed separately and painfully */
    for (pg->imtrx = 0; pg->imtrx < upd->Total_Num_Matrices; pg->imtrx++) {
      if (ERROR_ZZ_VEL != -1 && Num_Var_In_Type[pg->imtrx][R_MOMENTUM1]) {
        tev++;
        if (ERROR_ZZ_VEL_ELSIZE != -1) {
          tev++;
        }
      }
      if (ERROR_ZZ_Q != -1 && Num_Var_In_Type[pg->imtrx][R_ENERGY]) {
        tev++;
        if (ERROR_ZZ_Q_ELSIZE != -1) {
          tev++;
        }
      }
    }
    check = 0;
    for (i = 0; i < upd->Num_Mat; i++) {
      if (pd_glob[i]->MeshMotion == LAGRANGIAN || pd_glob[i]->MeshMotion == DYNAMIC_LAGRANGIAN)
        check = 1;
    }

    for (pg->imtrx = 0; pg->imtrx < upd->Total_Num_Matrices; pg->imtrx++) {
      if (ERROR_ZZ_P != -1 && (Num_Var_In_Type[pg->imtrx][R_MOMENTUM1] || check)) {
        tev++;
        if (ERROR_ZZ_P_ELSIZE != -1) {
          tev++;
        }
      }
    }
    /* Finally pick up all of the element-level-storage continuation
     * variables, e.g. for saturation hysteresis function
     */
    mat_num = Matilda[eb_indx];
    mp = mp_glob[mat_num];
    eb_ptr = Element_Blocks + eb_indx;
    ip_total = elem_info(NQUAD, eb_ptr->Elem_Type);
    if ((mp->PorousMediaType == POROUS_UNSATURATED ||
         mp->PorousMediaType == POROUS_SHELL_UNSATURATED ||
         mp->PorousMediaType == POROUS_TWO_PHASE) &&
        mp->SaturationModel == TANH_HYST && !if_ev) {
      for (j = 0; j < ip_total; j++) {
        if (SAT_CURVE_TYPE != -1)
          tev++; /*For Sat curve type */
        if (CAP_PRESS_SWITCH != -1)
          tev++; /*For saturation switch */
        if (SAT_QP_SWITCH != -1)
          tev++; /*for cap press switch point*/
      }
      if_ev = TRUE;
    }
  }

  /* Sanity check */
  if (tev != total_nev) {
    sr =
        sprintf(err_msg, "%s: Elem var count mismatch: tev(%d)<>rd->nev(%d)!?", yo, tev, total_nev);
    GOMA_EH(GOMA_ERROR, err_msg);
    /*
    fprintf(stderr,
            "%s: Disagreement over number of element variables\n",
            yo );
    exit (-1);
    */
  }

  /* Now do the real loop and populate the truth table */
  i = 0;
  for (eb_indx = 0; eb_indx < exo->num_elem_blocks; eb_indx++) {
    /* First test for all the potential elem vars from primary nodal vars
       for this block */
    mat_num = Matilda[eb_indx];
    ev_indx = 0;
    for (pg->imtrx = 0; pg->imtrx < upd->Total_Num_Matrices; pg->imtrx++) {
      for (j = V_FIRST; j < V_LAST; j++) {
        found_match = FALSE;
        if (pd_glob[mat_num]->v[pg->imtrx][j] != V_NOTHING) {
          if (pd_glob[mat_num]->i[pg->imtrx][j] == I_P0) {
            if (Num_Var_In_Type[pg->imtrx][j] > 1) {
              fprintf(stderr,
                      "%s: Too many components in variable type for element variable %s (%s)\n", yo,
                      Exo_Var_Names[j].name2, Exo_Var_Names[j].name1);
              exit(-1);
            }
            /* We just found a candidate for an element variable */
            exo->elem_var_tab[i++] = 1;
            found_match = TRUE;
            ev_indx++;
            /* malloc the entry for this block by number of elems for this block
               but - only if the variable exists for this block! (by the truth table) */
            if (has_been_called == 0) {
              /* NOTE: this final array dim is only to be malloc'd once; when a user
                 is annealing the mesh, anneal mesh calls wr_result_prelim_exo again,
                 and hence create_truth_table, which would realloc this dim of gvec_elem.
                 this test will prevent that. - RRL */
              asdv(&gvec_elem[pg->imtrx][eb_indx][ev_indx - 1], exo->eb_num_elems[eb_indx]);
            }
          } else if (pd_glob[mat_num]->i[pg->imtrx][j] == I_P1) {
            int dof = getdofs(type2shape(exo->eb_elem_itype[eb_indx]), I_P1);
            /* We just found a candidate for an element variable */
            for (int k = 0; k < dof; k++) {
              exo->elem_var_tab[i++] = 1;
              found_match = TRUE;
              ev_indx++;
              /* malloc the entry for this block by number of elems for this block
                 but - only if the variable exists for this block! (by the truth table) */

              if (has_been_called == 0) {
                /* NOTE: this final array dim is only to be malloc'd once; when a user
                   is annealing the mesh, anneal mesh calls wr_result_prelim_exo again,
                   and hence create_truth_table, which would realloc this dim of gvec_elem.
                   this test will prevent that. - RRL */
                asdv(&gvec_elem[pg->imtrx][eb_indx][ev_indx - 1], exo->eb_num_elems[eb_indx]);
              }
            }
          }
        }
        if (found_match == FALSE && exo->truth_table_existance_key[j - V_FIRST] == 1) {
          exo->elem_var_tab[i++] = 0;
          ev_indx++;
        }
      }
    }
    /* Now pick up all the post processing variables for this block
       - yes, for now they must
       each be listed separately and painfully */

    for (pg->imtrx = 0; pg->imtrx < upd->Total_Num_Matrices; pg->imtrx++) {
      if (ERROR_ZZ_VEL != -1 && Num_Var_In_Type[pg->imtrx][R_MOMENTUM1]) {
        exo->elem_var_tab[i++] = 1;
        ev_indx++;
        /* malloc the entry for this block by number of elems for this block
           but - only if the variable exists for this block! (by the truth table) */
        if (has_been_called == 0) {
          /* NOTE: this final array dim is only to be malloc'd once; when a user
             is annealing the mesh, anneal mesh calls wr_result_prelim_exo again,
             and hence create_truth_table, which would realloc this dim of gvec_elem.
             this test will prevent that. - RRL */
          asdv(&gvec_elem[pg->imtrx][eb_indx][ev_indx - 1], exo->eb_num_elems[eb_indx]);
        }
        if (ERROR_ZZ_VEL_ELSIZE != -1) {
          exo->elem_var_tab[i++] = 1;
          ev_indx++;
          /* malloc the entry for this block by number of elems for this block
             but - only if the variable exists for this block! (by the truth table) */
          if (has_been_called == 0) {
            /* NOTE: this final array dim is only to be malloc'd once; when a user
               is annealing the mesh, anneal mesh calls wr_result_prelim_exo again,
               and hence create_truth_table, which would realloc this dim of gvec_elem.
               this test will prevent that. - RRL */
            asdv(&gvec_elem[pg->imtrx][eb_indx][ev_indx - 1], exo->eb_num_elems[eb_indx]);
          }
        }
      }

      if (ERROR_ZZ_Q != -1 && Num_Var_In_Type[pg->imtrx][R_ENERGY]) {
        exo->elem_var_tab[i++] = 1;
        ev_indx++;
        /* malloc the entry for this block by number of elems for this block
           but - only if the variable exists for this block! (by the truth table) */
        if (has_been_called == 0) {
          /* NOTE: this final array dim is only to be malloc'd once; when a user
             is annealing the mesh, anneal mesh calls wr_result_prelim_exo again,
             and hence create_truth_table, which would realloc this dim of gvec_elem.
             this test will prevent that. - RRL */
          asdv(&gvec_elem[pg->imtrx][eb_indx][ev_indx - 1], exo->eb_num_elems[eb_indx]);
        }
        if (ERROR_ZZ_Q_ELSIZE != -1) {
          exo->elem_var_tab[i++] = 1;
          ev_indx++;
          /* malloc the entry for this block by number of elems for this block
             but - only if the variable exists for this block! (by the truth table) */
          if (has_been_called == 0) {
            /* NOTE: this final array dim is only to be malloc'd once; when a user
               is annealing the mesh, anneal mesh calls wr_result_prelim_exo again,
               and hence create_truth_table, which would realloc this dim of gvec_elem.
               this test will prevent that. - RRL */
            asdv(&gvec_elem[pg->imtrx][eb_indx][ev_indx - 1], exo->eb_num_elems[eb_indx]);
          }
        }
      }
    }
    check = 0;
    for (iii = 0; iii < upd->Num_Mat; iii++) {
      if (pd_glob[iii]->MeshMotion == LAGRANGIAN || pd_glob[iii]->MeshMotion == DYNAMIC_LAGRANGIAN)
        check = 1;
    }

    for (pg->imtrx = 0; pg->imtrx < upd->Total_Num_Matrices; pg->imtrx++) {
      if (ERROR_ZZ_P != -1 && (Num_Var_In_Type[pg->imtrx][R_MOMENTUM1] || check)) {
        exo->elem_var_tab[i++] = 1;
        ev_indx++;
        /* malloc the entry for this block by number of elems for this block
           but - only if the variable exists for this block! (by the truth table) */
        if (has_been_called == 0) {
          /* NOTE: this final array dim is only to be malloc'd once; when a user
             is annealing the mesh, anneal mesh calls wr_result_prelim_exo again,
             and hence create_truth_table, which would realloc this dim of gvec_elem.
             this test will prevent that. - RRL */
          asdv(&gvec_elem[pg->imtrx][eb_indx][ev_indx - 1], exo->eb_num_elems[eb_indx]);
        }
        if (ERROR_ZZ_P_ELSIZE != -1) {
          exo->elem_var_tab[i++] = 1;
          ev_indx++;
          /* malloc the entry for this block by number of elems for this block
             but - only if the variable exists for this block! (by the truth table) */
          if (has_been_called == 0) {
            /* NOTE: this final array dim is only to be malloc'd once; when a user
               is annealing the mesh, anneal mesh calls wr_result_prelim_exo again,
               and hence create_truth_table, which would realloc this dim of gvec_elem.
               this test will prevent that. - RRL */
            asdv(&gvec_elem[pg->imtrx][eb_indx][ev_indx - 1], exo->eb_num_elems[eb_indx]);
            has_been_called++;
          }
        }
      }
    }
    /*Now finally the saturation hysteresis variables */
    if (SAT_CURVE_TYPE != -1 || CAP_PRESS_SWITCH != -1 || SAT_QP_SWITCH != -1) {
      eb_ptr = Element_Blocks + eb_indx;
      ip_total = elem_info(NQUAD, eb_ptr->Elem_Type);
      for (j = 0; j < ip_total; j++) {
        /*Note that we will set these for all 3 var types because you
         *will never see them individually.
         */
        exo->elem_var_tab[i++] = 1;
        ev_indx++;
        if (has_been_called == 0) {
          asdv(&gvec_elem[pg->imtrx][eb_indx][ev_indx - 1], exo->eb_num_elems[eb_indx]);
        }
        exo->elem_var_tab[i++] = 1;
        ev_indx++;
        if (has_been_called == 0) {
          asdv(&gvec_elem[pg->imtrx][eb_indx][ev_indx - 1], exo->eb_num_elems[eb_indx]);
        }
        exo->elem_var_tab[i++] = 1;
        ev_indx++;
        if (has_been_called == 0) {
          asdv(&gvec_elem[pg->imtrx][eb_indx][ev_indx - 1], exo->eb_num_elems[eb_indx]);
        }
      }
    }
  }

  /* write out table */
  error = ex_put_truth_table(exo->exoid, EX_ELEM_BLOCK, exo->num_elem_blocks, total_nev,
                             exo->elem_var_tab);
  GOMA_EH(error, "ex_put_truth_table EX_ELEM_BLOCK");

  /* Now set truth table exists flag */
  exo->elem_var_tab_exists = TRUE;
}
/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
/*
 * add_qa_stamp() -- if space exists, add a goma QA record to the EXODUS II db.
 *
 * Notes:    [1] A new QA record is created that is slightly larger that the
 *		 old one. The old information is transcribed and a new record
 *		 is added that reflects this analysis.
 *
 *	     [2] New records are not added if there are at least 4 QA records
 *		 already accumulated. The threshhold varies according to the
 *		 C preprocessor symbol MAX_QA.
 *
 *		c.f. SAND83-0905
 *
 * Created: 1997/08/04 08:17 MDT pasacki@sandia.gov
 *
 * Revised: 2000/02/07 12:40 MST pasacki@sandia.gov
 */

void add_qa_stamp(Exo_DB *exo) {
  int i;
  int j;
  int k;
  int n;

  /*  char *new[MAX_QA][4];*/

  QA_Record *Q;

  n = exo->num_qa_rec;

  if (n > MAX_QA - 1) {
    return;
  }

  exo->num_qa_rec++;

  Q = (QA_Record *)smalloc((exo->num_qa_rec) * sizeof(QA_Record));

  for (i = 0; i < exo->num_qa_rec; i++) {
    for (j = 0; j < 4; j++) {
      Q[i][j] = (char *)smalloc(LEN_QA_RECORD * sizeof(char));

      /*
       * Initialize to null terminators...
       */

      for (k = 0; k < LEN_QA_RECORD; k++) {
        Q[i][j][k] = '\0';
      }
    }
  }

  /*
   * Copy any old stuff into the new area.
   */

  for (i = 0; i < exo->num_qa_rec - 1; i++) {
    for (j = 0; j < 4; j++) {
      strcpy(Q[i][j], exo->qa_record[i][j]);
    }
  }

  /*
   * Preload new record with all terminating nulls.
   */

  for (j = 0; j < 4; j++) {
    for (k = 0; k < MAX_STR_LENGTH; k++) {
      Q[n][j][k] = '\0';
    }
  }

  /*
   * Concoct a new record.
   */

  strcpy(Q[n][0], "GOMA");
  strcpy(Q[n][1], GOMA_VERSION); /* def'd in std.h for now */
  get_date(Q[n][2]);
  get_time(Q[n][3]);

  /*
   * Free the old record members.
   */

  for (i = 0; i < n; i++) {
    for (j = 0; j < 4; j++) {
      safer_free((void **)&(exo->qa_record[i][j]));
    }
  }

  if (n > 0) {
    safer_free((void **)&(exo->qa_record));
  }

  /*
   * Assign new pointer.
   */

  exo->qa_record = Q;

  return;
}
/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
/*
 * add_info_stamp -- add an information record to EXODUS II db if space exists.
 *
 * Created: 1997/08/04 08:22 MDT pasacki@sandia.gov
 */

void add_info_stamp(Exo_DB *exo) {
  int i;
  int k;
  int n;
  char **a;
  char buf[MAX_LINE_LENGTH + 1];
  time_t now, then;
#ifdef NO_LEAKY_GETPWUID
  struct passwd *pwe;
#endif
  struct utsname utsname;
  INFO_Record *I;

  n = exo->num_info;

  if (n + 10 > MAX_INFO) {
    return;
  }

  exo->num_info += 10;

  /*
  buf[0] = '\0';
  ni = (char **) calloc( exo->num_info, sizeof(char *));
  */

  /*
   * Allocate space for the new info record.
   */

  I = (INFO_Record *)smalloc(exo->num_info * sizeof(INFO_Record));

  for (i = 0; i < exo->num_info; i++) {
    I[i] = (char *)smalloc((MAX_LINE_LENGTH + 1) * sizeof(char));
  }

  /*
   * Transcribe any old records...
   */

  for (i = 0; i < n; i++) {
    strcpy(I[i], exo->info[i]);
  }

  /*
   * Initialize the new records to terminating nulls...
   */

  for (i = n; i < exo->num_info; i++) {
    for (k = 0; k < MAX_LINE_LENGTH; k++) {
      I[i][k] = '\0';
    }
  }

  /*
   * Fill in the new records with information about this run.
   */

  strcpy(I[n], "____");

  /*
   * -9 -- the command line issued for this simulation
   */

  /*
  for ( i=0; i<Argc; i++)
    {
      strcat(buf, Argv[i]);
      strcat(buf, " ");
    }
  */

  for (i = 0; i < MAX_LINE_LENGTH; i++) {
    buf[i] = '\0';
  }

  a = Argv;

  if (a != NULL) {
    i = 0;
    /* MMH: This needs a +2: +1 for the space, and +1 for the
     * null.
     */
    while (*a != NULL && ((strlen(buf) + strlen(*a) + 2) < MAX_LINE_LENGTH)) {
      strcat(buf, *a);
      strcat(buf, " ");
      a++;
    }
  }

  strcpy(I[n + 1], buf);

  /*
   * -8 -- the date and time of the simulation
   */

  now = time(&then);
  strftime(buf, MAX_LINE_LENGTH, "%C", localtime(&now));
  strcpy(I[n + 2], buf);

  /*
   * -7 -- current working directory
   */

  if (ProcID < 8) /* too much I/O overhead for many procs */
  {
    char *cwderr = getcwd(buf, MAX_LINE_LENGTH + 1);
    if (cwderr == NULL) {
      strcpy(buf, ".");
    }
  } else {
    strcpy(buf, ".");
  }

  strncpy(I[n + 3], buf, MAX_LINE_LENGTH + 1);

  /*
   * -6 -- the name of the user
   */

  sprintf(buf, "uid %d", (int)(getuid()));

#ifdef NO_LEAKY_GETPWUID
  pwe = getpwuid(getuid());
  strcpy(buf, pwe->pw_name);
#endif

  strcpy(I[n + 4], buf);

  /*
   * -5 through -1 -- the POSIX system information
   */

  uname(&utsname);

  /* Lets be on the safe side put one char less than MAX_LINE_LENGTH into the info buffers */

  strncpy(I[n + 5], utsname.sysname, MAX_LINE_LENGTH - 1);
  strncpy(I[n + 6], utsname.nodename, MAX_LINE_LENGTH - 1);
  strncpy(I[n + 7], utsname.release, MAX_LINE_LENGTH - 1);
  strncpy(I[n + 8], utsname.version, MAX_LINE_LENGTH - 1);
  strncpy(I[n + 9], utsname.machine, MAX_LINE_LENGTH - 1);

  /*
   * Free the old beast and assign the new one.
   */
  for (i = 0; i < n; i++) {
    safer_free((void **)&(exo->info[i]));
  }
  if (n > 0) {
    safer_free((void **)&(exo->info));
  }
  exo->info = I;
  return;
}

void wr_result_exo(
    Exo_DB *exo, char *filename, int verbosity, int write_node_vars, int write_elem_vars) {
  int i;
  int index;
  int j;
  int k;
  int status;
  int time_index;
  char err_msg[MAX_CHAR_IN_INPUT];
  /*
   * This file should already exist.
   */

  exo->cmode = EX_WRITE;

#ifdef DEBUG
  fprintf(stderr, "%s: begins\n", yo);
#endif

  exo->io_wordsize = 0; /* i.e., query */
  exo->comp_wordsize = sizeof(dbl);
  exo->exoid = ex_open(filename, exo->cmode, &exo->comp_wordsize, &exo->io_wordsize, &exo->version);

#ifdef DEBUG
  fprintf(stderr, "\t\tfilename    = \"%s\"\n", filename);
  fprintf(stderr, "\t\tcomp_ws     = %d\n", exo->comp_wordsize);
  fprintf(stderr, "\t\tio_wordsize = %d\n", exo->io_wordsize);
#endif

  /*
   * Element variable truth table and values at ONE TIME ONLY.
   */

  if (exo->num_elem_vars > 0 && write_elem_vars) {

#ifdef DEBUG
    fprintf(stderr, "\t\tneb         = %d\n", exo->num_elem_blocks);
    fprintf(stderr, "\t\tnev         = %d\n", exo->num_elem_vars);
    fprintf(stderr, "\t\tevt:        =   \n");
    for (i = 0; i < exo->num_elem_blocks; i++) {
      for (j = 0; j < exo->num_elem_vars; j++) {
        fprintf(stderr, "block index %d, elem var index %d is %d\n", i, j,
                exo->elem_var_tab[i * (exo->num_elem_vars) + j]);
      }
    }
#endif

    /*
     * This has already been done.
     *
    status = ex_put_elem_var_tab(exo->exoid,
                                 exo->num_elem_blocks,
                                 exo->num_elem_vars,
                                 exo->elem_var_tab);
    GOMA_EH(status, "ex_put_elem_var_tab");
    */

    for (i = 0; i < exo->num_ev_time_indeces; i++) {
      time_index = exo->ev_time_indeces[i];

      for (j = 0; j < exo->num_elem_blocks; j++) {
        for (k = 0; k < exo->num_elem_vars; k++) {
          index = j * exo->num_elem_vars + k;

          if (exo->elem_var_tab == NULL || exo->elem_var_tab[index] != 0) {
            status = ex_put_var(exo->exoid, time_index, EX_ELEM_BLOCK, k + 1, exo->eb_id[j],
                                exo->eb_num_elems[j], &(exo->ev[i][index][0]));
            if (status < 0) {
              sprintf(err_msg, "ex_put_var() elem bad rtn: time %d, elemvar %d, EB ID %d",
                      time_index, k + 1, exo->eb_id[j]);
              GOMA_EH(GOMA_ERROR, err_msg);
            }
          }
        }
      }
    }
  }

  /*
   * Put nodal variable values at last time step...
   */

  if (exo->num_node_vars > 0 && write_node_vars) {
    for (i = 0; i < exo->num_nv_time_indeces; i++) {
      time_index = exo->nv_time_indeces[i];

      for (j = 0; j < exo->num_nv_indeces; j++) {
        status = ex_put_var(exo->exoid, time_index, EX_NODAL, exo->nv_indeces[j], 1, exo->num_nodes,
                            &(exo->nv[i][j][0]));
        GOMA_EH(status, "ex_put_var nodal");
      }
    }
  }

  status = ex_close(exo->exoid);
  GOMA_EH(status, "ex_close()");

  return;
}

/* wr_resetup_exo() -- open/write/close EXODUS II db for results names
 *
 * Created: 1998/01/26 14:06 MST pasacki@sandia.gov
 *
 * Revised:
 */

void wr_resetup_exo(Exo_DB *exo, char *filename, int verbosity) {
  int error;
  int i;
  int status;

  /*
   * This file must already exist.
   */

  exo->cmode = EX_WRITE;

#ifdef DEBUG
  fprintf(stderr, "%s: begins\n", yo);
#endif

  exo->io_wordsize = 0; /* i.e., query */
  exo->comp_wordsize = sizeof(dbl);
  exo->exoid = ex_open(filename, exo->cmode, &exo->comp_wordsize, &exo->io_wordsize, &exo->version);

#ifdef DEBUG
  fprintf(stderr, "\t\tfilename    = \"%s\"\n", filename);
  fprintf(stderr, "\t\tcomp_ws     = %d\n", exo->comp_wordsize);
  fprintf(stderr, "\t\tio_wordsize = %d\n", exo->io_wordsize);
#endif

  /*
   * Results setup...
   */

  if (exo->num_glob_vars > 0) {
    status = ex_put_variable_param(exo->exoid, EX_GLOBAL, exo->num_glob_vars);
    GOMA_EH(status, "ex_put_variable_param global");
    status = ex_put_variable_names(exo->exoid, EX_GLOBAL, exo->num_glob_vars, exo->glob_var_names);
    GOMA_EH(status, "ex_put_variable_names global");
  }

  if (exo->num_elem_vars > 0) {
    status = ex_put_variable_param(exo->exoid, EX_ELEM_BLOCK, exo->num_elem_vars);
    GOMA_EH(status, "ex_put_variable_param elem block");
    status =
        ex_put_variable_names(exo->exoid, EX_ELEM_BLOCK, exo->num_elem_vars, exo->elem_var_names);
    GOMA_EH(status, "ex_put_variable_names elem block");
    if (exo->elem_var_tab != NULL) {
      status = ex_put_truth_table(exo->exoid, EX_ELEM_BLOCK, exo->num_elem_blocks,
                                  exo->num_elem_vars, exo->elem_var_tab);
      GOMA_EH(status, "ex_put_truth_table elem block");
    }
  }

  if (exo->num_node_vars > 0) {
    status = ex_put_variable_param(exo->exoid, EX_NODAL, exo->num_node_vars);
    GOMA_EH(status, "ex_put_variable_param nodal");
    status = ex_put_variable_names(exo->exoid, EX_NODAL, exo->num_node_vars, exo->node_var_names);
    GOMA_EH(status, "ex_put_variable_names nodal");
  }

  if (exo->num_times > 0) {
    for (i = 0; i < exo->num_times; i++) {
      status = ex_put_time(exo->exoid, i + 1, &(exo->time_vals[i]));
      GOMA_EH(status, "ex_put_times");
    }
  }

  error = ex_close(exo->exoid);
  if (error != 0)
    exit(2);
  return;
}