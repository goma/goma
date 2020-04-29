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
 *$Id: rf_pre_proc.c,v 5.1 2007-09-18 18:53:47 prschun Exp $
 */

#include "el_elm_info.h"
#include "el_geom.h"
#include "exo_struct.h"
#include "mm_bc.h"
#include "rf_bc.h"
#include "rf_fem.h"
#include "rf_fem_const.h"
#include "rf_pre_proc.h"

#define GOMA_RF_PRE_PROC_C

static void build_elem_type_list
(Exo_DB *);		/* exo */

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void 
pre_process(Exo_DB *exo)

/* pre_process() -- counts up unique element types, some bc initialization
 *      
 *      Author:          Scott Hutchinson (1421)
 *      Date:            8 January 1993
 *
 *      Revised: 1997/08/28 16:54 MDT pasacki@sandia.gov
 */
{

  build_elem_type_list(exo);

  alloc_First_Elem_BC (&First_Elem_Side_BC_Array, 
		       &First_Elem_Edge_BC_Array, Num_Internal_Elems);

  return;
} /* END of routine pre_process  */
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/ 
/* build_elem_type_list() -- build Unique_Element_Types list for the problem
 *
 * Description:
 *	Build up the integer array that describes how many different basis
 *	functions we will need to consider for this problem. This makes use
 *	of the different element shapes as well as the degree of the
 *	interpolation selected for the different variables. This information
 *	is used later so that only as many basis functions as absolutely
 *	necessary are built up during assembly.
 *		
 *	The affected global variables Num_Element_Types and Unique_Element_Types
 *	are external, defined in "rf_bc.h".
 *		
 * Created: 1997/08/04 14:20 MDT pasacki@sandia.gov
 *
 * Revised: 1997/08/28 16:51 MDT pasacki@sandia.gov
 */

static void 
build_elem_type_list(Exo_DB *exo)
{
  int i;
  int ielem_type;

  Num_Element_Types = 0;
  for ( i=0; i<MAX_ELEMENT_TYPES; i++)
    {
      Unique_Element_Types[i] = -1;
    }

  /* loop through all the element blocks on the processor, setting the element
     type for each element in each block */

  for (i = 0; i < exo->num_elem_blocks; i++)
    {
      ielem_type = exo->eb_elem_itype[i];

      /*
       * If this element type has not been encountered before, then
       * add it to our list of unique element types...
       */
      
      if ( in_list( ielem_type, 0, Num_Element_Types, 
		    Unique_Element_Types) == -1 )
	{
	  Unique_Element_Types[Num_Element_Types] = ielem_type;
	  Num_Element_Types++;
	}
    }

  return;

} /* end of build_elem_type_list() */
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/
