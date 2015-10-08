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
 *$Id: rf_node.h,v 5.1 2007-09-18 18:53:47 prschun Exp $
 */

#ifndef _RF_NODE_H
#define _RF_NODE_H

#include "rf_node_const.h"

/*
 *  Global Pointer to the Vector of Pointers to Node Structures
 *  for all of the nodes on the current processor
 *  Index of the array is the matrix number and processor (local) node number.
 *
 */
NODE_INFO_STRUCT **Nodes = NULL;

#endif


/*****************************************************************************/
/*                       STRUCTURE DEFINITIONS                               */
/*****************************************************************************/

