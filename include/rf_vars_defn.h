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
 *$Id: rf_vars_defn.h,v 5.2 2008-10-02 15:36:29 hkmoffa Exp $
 */

#ifndef _RF_VARS_DEFN_H
#define _RF_VARS_DEFN_H

/*
 * Var_Info Array
 *
 * Array of VARIABLE_DSCRIPTION structures which is indexed by the
 * tags at each node contained in the Variable_Tag array defined
 * above.  These entries store information on
 * the existence, type, size and ordering of variable information at nodes.
 */

VARIABLE_DESCRIPTION_STRUCT **Var_Info = NULL;

/* Global number of different Var_Info entries */

int Num_Var_Info_Records = 0;

NODAL_VARS_STRUCT **Nodal_Vars_List = NULL;
int Nodal_Vars_List_Length = 0;

/*********************end of rf_vars_defn.h *********************************/

#endif
