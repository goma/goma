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

/* rd_input() -- read the brk input file, checking against .exoII data, allocate
 *
 * Notes:
 *		[1] Space for the special data structures are allocated
 *                  and filled, checking against the EXODUS II file as
 *                  appropriate.
 *
 * Created: 1997/05/08 09:11 MDT pasacki@sandia.gov
 *
 * Revised:
 */

#define _RD_IN_C

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <string.h>

#ifdef STDC_HEADERS
#include <stdlib.h>
#endif

#include "map_names.h"
#include "std.h"
#include "aalloc.h"
#include "eh.h"
#include "exo_struct.h"
#include "brkfix_types.h"
#include "rd_in.h"
#include "utils.h"

const char delimiters[]=" 	\n"; /* whitespace delimiters (space, tab) */

void 
rd_input(char *in_file_name,	/* name of the input file */
	 Exo_DB *mono,		/* monolithic database */
	 Bevm ****p_mult,	/* basic eqnvar multiplicity */
	 int ****p_evd,		/* eqnvar dependencies */
	 int ****p_Lucky,		/* local node/dof existences */
	 int **p_num_basic_eqnvars) /* for ea element block */
{
  /*
   * Since space is allocated for these for a higher level routine, declare
   * these convenience variables to correspond to the pointers that were
   * passed who will receive addresses of space.
   */

  Bevm ***mult;
  int ***evd;
  int ***Lucky;
  int *num_basic_eqnvars;

  int beqn;
  int bvar;

  int done;

  int eb;
  int *eb_checklist;
  int ebid;
  int eb_index;
  int err;
  int ev;
  int evid;

  int i;
  int index;
  int ival;

  int j;

  int k;

  int lazy_node;

  int vmul;
  int cmul;
  int nmul;

  int nbev;

  /*  int ne;*/
  int neb;
  /*  int nn;*/
  /*  int nns;*/
  int npe;
  /*  int nss;*/

  char  line[LINE_BUFFER_LENGTH];
  char  err_msg[MAX_CHAR_ERR_MSG];
  char  string[LINE_BUFFER_LENGTH];
  char *ptr;
  char *rtn;
  char *tmp;

  FILE *in_strm;

  Spfrtn sr=0;

  in_strm = fopen(in_file_name, "r");

  if ( in_strm == NULL )
    {
      EH(-1, "Problem opening input file.");
    }

  /*
   * Convenience variables...
   */

  /*
   * ne  = mono->num_elems;
   * nn  = mono->num_nodes;
   */

  neb = mono->num_elem_blocks;

  /*
   * nns = mono->num_node_sets;
   * nss = mono->num_side_sets;
   */

  /*
   * Allocate gross space.
   */

  *p_mult              = (Bevm ***) smalloc(neb*sizeof(Bevm **));
  mult                 = *p_mult;

  *p_evd               = (int ***) smalloc(neb*sizeof(int **));
  evd                  = *p_evd;

  *p_Lucky             = (int ***) smalloc(neb*sizeof(int **));
  Lucky                = *p_Lucky;


  *p_num_basic_eqnvars = (int *) smalloc(neb*SZ_INT);
  num_basic_eqnvars    = *p_num_basic_eqnvars;

  eb_checklist         = (int *) smalloc(neb*SZ_INT);

  INIT_IVEC(num_basic_eqnvars, -1, neb);

  INIT_IVEC(eb_checklist, -1, neb);

  /*
   * In the input file, the active equations and variables are described
   * on the basis of each element block.
   */

  for ( eb=0; eb<neb; eb++)
    {
      rtn = fgets(line, LINE_BUFFER_LENGTH, in_strm);
      err = sscanf(line, "%d", &ebid);
#ifdef DEBUG
	  fprintf(stderr, "Reading input info for eb_id = %d\n", ebid);
	  fprintf(stderr, "eb_index=%d, eb_ID=%d\n", eb, ebid);
#endif

      if ( err != 1 )
	{
	  sr = sprintf(err_msg, "Trouble reading EB ID at eb=%d", eb);
	  EH(-1, err_msg);
	}
      
      /*
       * Verify this ebid exists in the fe db and has not been encountered
       * previously.
       */
      
      if ( !IUL(ebid, mono->eb_id, neb) )
	{
	  sr = sprintf(err_msg, "EB_ID = %d not in %s", ebid, mono->path);
	  EH(-1, err_msg);
	}

      for ( i=0; i<eb; i++ )
	{
	  if ( eb_checklist[i] != -1 )
	    {
	      if ( ebid == eb_checklist[i] )
		{
		  sr = sprintf(err_msg, "EB_ID = %d found twice.", ebid);
		  EH(-1, err_msg);
		}
	    }
	}

      eb_checklist[eb] = ebid;

      rtn = fgets(line, LINE_BUFFER_LENGTH, in_strm);
      err = sscanf(line, "%d", &nbev);

#ifdef DEBUG
      fprintf(stdout, "Number of basic eqnvars this eb = %d\n", nbev);
#endif

      if ( err != 1 )
	{
	  EH(-1, "Trouble reading number of basic eqnvars.");
	}

      if ( nbev < 1 )
	{
	  sr = sprintf(err_msg, "Input file says %d basic eqnvars.", nbev);
	  EH(-1, err_msg);
	}

      num_basic_eqnvars[eb] = nbev;

      /*
       * Get ready to read in the multiplicities for the basic equation/
       * variables that are present in this element block. The different kinds
       * of multiplicities exist for each of the basic equation/variables.
       *
       * Basic equations differ somewhat from equations in the goma sense.
       * Those equations are listed multiple times for vectors, for example,
       * where they are listed but once here.
       */

      mult[eb] = (Bevm **) smalloc(nbev*sizeof(Bevm *));

#ifdef DEBUG
      fprintf(stdout, "Starting loop over basic eqnvars in ebid=%d\n",
	      ebid);
#endif

      for ( ev=0; ev<nbev; ev++)
	{
	  /*
	   * Allocate and initialize...
	   */

	  mult[eb][ev] = (Bevm *) smalloc(sizeof(Bevm));

	  mult[eb][ev]->eqnvar_id = UNDEFINED_EQNVARID;
	  mult[eb][ev]->vect_mult = 0;
	  mult[eb][ev]->conc_mult = 0;
	  mult[eb][ev]->ndof_mult = 0;

	  /*
	   * Note! The eqnvar_ids are arbitrary integers that should uniquely
	   * identify each basic variable equation kind for the entire
	   * problem. The implication is that variables of the same eqnvar_id
	   * are continuous across interblock boundaries.
	   *
	   * Also, the ordering scheme used to order the unknowns is based
	   * on the numerical order implied by the choice of eqnvar_id.
	   */

	  rtn = fgets(line, LINE_BUFFER_LENGTH, in_strm);

	  err = sscanf(line, "%d %d %d %d", &evid, &vmul, &cmul, &nmul);

#ifdef DEBUG
	  fprintf(stdout, "eqnvar_id = %d, multiplicities = %d %d %d\n",
		  evid, vmul, cmul, nmul);
#endif

	  if ( err != 4 )
	    {
	      EH(-1, "Trouble reading group of eqnvarid's with multipliers.");
	    }

	  /*
	   * Insure this variable id has not appeared before in this eb's list
	   */
	  for ( k=0; k<ev; k++)
	    {
	      if ( mult[eb][k]->eqnvar_id == evid )
		{
		  EH(-1, "Multiple occurence of basic variable in this eb.");
		}
	    }

	  mult[eb][ev]->eqnvar_id = evid;

	  if ( vmul > 0 )
	    {
	      mult[eb][ev]->vect_mult = vmul;
	    }
	  else
	    {
	      EH(-1, "Vector multiplier too low.");
	    }

	  if ( cmul > 0 )
	    {
	      mult[eb][ev]->conc_mult = cmul;
	    }
	  else
	    {
	      EH(-1, "Concentration multiplier too low.");
	    }

	  if ( nmul > 0 )
	    {
	      mult[eb][ev]->ndof_mult = nmul;
	    }
	  else
	    {
	      EH(-1, "Nodal dof multiplier too low.");
	    }
	}	  

      /*
       * The equation/variable dependency matrix for this element block
       */

      evd[eb] = (int **) smalloc(nbev * sizeof(int *));
      evd[eb][0] = (int *) smalloc(nbev * nbev * sizeof(int));

      for ( i=1; i<nbev; i++)
	{
	  evd[eb][i] = evd[eb][i-1] + nbev;
	}




      for ( beqn=0; beqn<nbev; beqn++)
	{
	  rtn = fgets(line, LINE_BUFFER_LENGTH, in_strm);
	  
	  tmp = line;
	  ptr = strtok(tmp, delimiters);

	  for ( bvar=0; bvar<nbev; bvar++)
	    {
	      if ( ptr == NULL )
		{
		  fprintf(stderr, "Expecting some interactions here.\n");
		}

	      err = sscanf(ptr, "%d", &(evd[eb][beqn][bvar]));
#ifdef DEBUG
	      fprintf(stderr, "Interaction (eqn,var)=(%d,%d) is %d\n",
		      beqn, bvar, evd[eb][beqn][bvar]);
#endif
	      if ( err != 1 )
		{
		  EH(-1, "Problem reading a eqn/var dependency.");
		}

	      ptr = strtok(NULL, delimiters);
	    }
	}

      /*
       * Find out what kind of highest order conglomerate element is
       * present in this element block. Not every variable is defined
       * at every node, but at least one of some kind of variable needs
       * to be defined at a given node.
       */

      rtn = fgets(line, LINE_BUFFER_LENGTH, in_strm);
      err = sscanf(line, "%d", &npe);

      if ( err != 1 )
	{
	  EH(-1, "Trouble determining nodes per element.");
	}
      
      if ( npe < 1 )
	{
	  EH(-1, "Too few nodes per element in this block.");
	}

#ifdef DEBUG
      fprintf(stderr, "\tNumber of nodes/elem this eb = %d\n", npe);
#endif

      /*
       * Verify that the number of nodes per element for this listed block
       * agrees with the number of nodes per element as indicated by
       * the EXODUS II db information...
       */

      eb_index = in_list(ebid, mono->eb_id, mono->num_elem_blocks);
      EH(eb_index, "Did not find specified element block ID in monolith.");

      if ( npe != mono->eb_num_nodes_per_elem[eb_index] )
	{
	  fprintf(stderr, 
		  "Inconsistent estimates of nodes per element for ebid %d.\n",
		  ebid);
	  fprintf(stderr, "Input file says %d nodes per element.\n", npe);
	  fprintf(stderr, "EXODUS file says %d nodes per element.\n", 
		  mono->eb_num_nodes_per_elem[eb_index]);
	  EH(-1, "Giving up.");
	}


      /*
       * Allocate space for and initialize the local node dof existence
       * matrix for this element block.
       */

      Lucky[eb]    = (int **) smalloc(npe * sizeof(int *));
      Lucky[eb][0] = (int *) smalloc(npe * nbev * sizeof(int));

      for ( i=1; i<npe; i++)
	{
	  Lucky[eb][i] = Lucky[eb][i-1] + nbev;
	}

      for ( i=0; i<npe; i++)
	{
	  for ( j=0; j<nbev; j++ )
	    {
	      Lucky[eb][i][j] = 0;
	    }
	}

      for ( i=0; i<npe; i++)
	{
	  lazy_node = FALSE;	/* give it the benefit of the doubt to begin */
	  rtn = fgets(line, LINE_BUFFER_LENGTH, in_strm);

	  if ( rtn == NULL )
	    {
	      fprintf(stderr, 
		      "Got squat reading the %d line of %d, EB ID %d\n", i+1,
		      npe, ebid);
	      EH(-1, "Giving up.");
	    }

#ifdef DEBUG
	  fprintf(stderr, "\t[%d] line=\"%s\"\n", i, line);
#endif

	  tmp = line;
	  ptr = strtok(tmp, delimiters);

	  done = (ptr == NULL);	/* i.e., no more args left in line */
	  
	  while ( !done )
	    {
#ifdef DEBUG
	      fprintf(stderr, "\tptr = \"%s\"\n", ptr);
#endif
	      /*
	       * Check for no active variables, specified with a "none"
	       */
	      err = sscanf(ptr, "%s", string);
	      
	      if ( err == 1 )
		{
		  lazy_node = ( strcmp(string, "none") == 0 );
		}

	      if ( ! lazy_node )
		{
		  err = sscanf(ptr, "%d", &ival);
		}
	      else
		{
		  err = 0;	/* facsimile of how many ints are here... */
		}

	      if ( err == 1 )	/* there is an eqnvar at this node */
		{
#ifdef DEBUG
		  fprintf(stderr, "\teqnvar_id = %d\n", ival);
#endif
		  /* 
		   * Find correct index for this eqnvarid in the mult list.
		   */
		  index = -1;
		  for ( k=0; k<nbev; k++)
		    {
		      if ( mult[eb][k]->eqnvar_id == ival )
			{
			  index = k;
			}
		    }
		  if ( index == -1 )
		    {
		      sr = sprintf(err_msg, "eqnvar_id=%d not in mult", ival);
		      EH(-1, err_msg);
		    }
		  else
		    {
		      Lucky[eb][i][index]++;
		    }
		  ptr = strtok(NULL, delimiters);
		}
	      else
		{
		  done = TRUE;
		}
	      done |= ( ptr == NULL );
	    }
	}
#ifdef DEBUG
      /*
       * Dump out local node dof existence profile for this eb.
       */
      for ( i=0; i<npe; i++)
	{
	  for ( j=0; j<nbev; j++ )
	    {
	      fprintf(stderr, "Lucky[eb=%d][ln=%d][eqnvar=%d] = %d\n",
		      eb, i, j, Lucky[eb][i][j]);
	    }
	}
#endif
    }

  fclose(in_strm);

  free(eb_checklist);

  if ( sr < 0 ) {
    fprintf( stdout,"rd_in: error reading the brk input file\n");
    exit(2);
  }

  return;
}

