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
 
#ifndef lint
static char *cvs_masksh_id =
  "$Id";
#endif

/*
 * Var_Info Array
 *
 * Array of VAR_STRUCT which is indexed by the tags at each node contained in
 * the Variable_Tag array defined above.  These entries store information on
 * the existence, type, size and ordering of variable information at nodes.
 */

VAR_STRUCT *Var_Info = NULL;

/* Global number of different Var_Info entries */

u_short Num_Var_Tags;

/*********************end of rf_mask.h****************************************/


