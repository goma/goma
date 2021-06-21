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
 *$Id: mm_elem_block.h,v 5.1 2007-09-18 18:53:42 prschun Exp $
 */

#ifndef GOMA_MM_ELEM_BLOCK_H
#define GOMA_MM_ELEM_BLOCK_H

#include "mm_elem_block_structs.h"

/*****************************************************************************/
/*                       STRUCTURE DEFINITIONS                               */
/*****************************************************************************/

/*
 * Element_Block_Struct: element block information.
 */


ELEM_BLK_STRUCT *Element_Blocks = NULL;       /* Pointer to array of global
                                                 element block information. */

ELEM_BLK_STRUCT *Current_EB_ptr = NULL; /* Pointer to the current element block
					   structure. This is calculated in the
					   fill */

#endif
