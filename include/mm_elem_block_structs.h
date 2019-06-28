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
 *$Id: mm_elem_block_structs.h,v 5.1 2007-09-18 18:53:42 prschun Exp $
 */

#ifndef _MM_ELEM_BLOCK_STRUCTS_H
#define _MM_ELEM_BLOCK_STRUCTS_H

#include "std.h"
#include "rf_element_storage_struct.h"

/*****************************************************************************/
/*                       STRUCTURE DEFINITIONS                               */
/*****************************************************************************/


/*
 * Element_Block_Struct: element block information.
 *
 *     Information about each element block known to the current processor
 *     is storred in this structure.
 */

struct Element_Block_Struct {
  int Elem_Blk_Num;                 /* The index of this element block struct
                                     * in the Element_Blocks array.          */
  int Elem_Blk_Id;                  /* The ExodusII element block ID number. */
  int Elem_Type;                    /* Integer type of the element in this
				     * element block. The types are defined
				     * in el_elm.h */
  int IP_total;                     /* Total number of volumetric quadrature 
				     * points to be used for volumetric 
				     * integrations */
  int Num_Elems_In_Block;           /* Number of elements in the element 
				     * block known to the current proccessor.
				     * Thus, this number may vary between
				     * processors */
  int Num_Nodes_Per_Elem;           /* Number of nodes per element for
				     * elements in this element block.   */
  int Num_Attr_Per_Elem;            /* Number of attributes per element in
                                     * the element block.                    */
  MATRL_PROP_STRUCT *MatlProp_ptr;  /* Pointer to the material property 
                                     * structure applicable to this element
				     * block                                 */
  ELEMENT_STORAGE_STRUCT *ElemStorage;
};
typedef struct Element_Block_Struct ELEM_BLK_STRUCT;

/*
 *  This is the global list of element blocks defined for this problem
 */
extern ELEM_BLK_STRUCT *Element_Blocks;

/*
 *  This is a pointer to the current element block structure that pertains
 *  to the current element being processed
 */
extern ELEM_BLK_STRUCT *Current_EB_ptr;

#endif
