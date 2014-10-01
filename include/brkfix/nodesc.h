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

/* nodesc.h -- describes the variables & multiplicities at a given node
 * 
 *
 * Each node can be associated with one of a handful of prototypes.
 * Each prototype is extensively described by the number and names of
 * the active eqnvars at the node, as well as the cumulative weights
 * that apply to each eqnvar.
 *
 * The structure is for convenience and normal usage; the integer vector
 * packs the same information so that I/O is more easily achieved of the
 * same information.

 * Created: 1997/05/18 08:12 MDT pasacki@sandia.gov
 *
 * Revised:
 */

#ifndef _NODESC_H
#define _NODESC_H


#define MAX_EQNVARS			16

#define MAX_NODE_KINDS			32

#define LEN_NODE_DESCRIPTION (4*MAX_EQNVARS+1)

struct Node_Description
{
  int num_basic_eqnvars;
  int eqnvar_ids[MAX_EQNVARS];		/* put them in ascending order! */
  int eqnvar_wts[MAX_EQNVARS][3];	/* [0]=vect, [1]=conc, [2]=nodaldof */
};

typedef struct Node_Description Node_Description;

#define SZ_ND (sizeof(struct Node_Description))

#endif
