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
 *$Id: rf_masks.h,v 5.1 2007-09-18 18:53:46 prschun Exp $
 */

#ifndef _RF_MASKS_H
#define _RF_MASKS_H

#ifndef MAX_TOTAL_BCS_POWER
#define MAX_TOTAL_BCS_POWER  8
#endif

/*
 *  Inter_Mask:
 * ---------------------------------------------------------------------------
 *
 *       mask array which tells which unknowns interact with which unknowns...
 *
 *  Example:  Inter_Mask[equation][variable]
 *
        U1      U2      U3      T       Y    D1      D2      D3       S       P
      |-----|-------|-------|------|-------|-----|-------|--------|-------|------|
   U1 |  1  |   1   |   0   |   1  |    1  |  1  |   1   |   0    |   0   |   1  |
      |-----|-------|-------|------|-------|-----|-------|--------|-------|------|
   U2 |  1  |   1   |   0   |   1  |    1  |  1  |   1   |   0    |   0   |   1  |
      |-----|-------|-------|------|-------|-----|-------|--------|-------|------|
   U3 |  0  |   0   |   0   |   0  |    0  |  0  |   0   |   0    |   0   |   0  |
      |-----|-------|-------|------|-------|-----|-------|--------|-------|------|
   T  |  1  |   1   |   0   |   1  |    1  |  1  |   1   |   0    |   0   |   1  |
      |-----|-------|-------|------|-------|-----|-------|--------|-------|------|
   Y  |  1  |   1   |   0   |   1  |    1  |  0  |   0   |   0    |   0   |   1  |
      |-----|-------|-------|------|-------|-----|-------|--------|-------|------|
   D1 |  0  |   0   |   0   |   0  |    0  |  1  |   1   |   0    |   0   |   0  |
      |-----|-------|-------|------|-------|-----|-------|--------|-------|------|
   D2 |  0  |   0   |   0   |   0  |    0  |  1  |   1   |   0    |   0   |   0  |
      |-----|-------|-------|------|-------|-----|-------|--------|-------|------|
   D3 |  0  |   0   |   0   |   0  |    0  |  0  |   0   |   0    |   0   |   0  |
      |-----|-------|-------|------|-------|-----|-------|--------|-------|------|
   S  |  0  |   0   |   0   |   0  |    0  |  0  |   0   |   0    |   1   |   0  |
      |-----|-------|-------|------|-------|-----|-------|--------|-------|------|
   P  |  1  |   1   |   0   |   1  |    1  |  1  |   1   |   0    |   0   |   1  |
      |-----|-------|-------|------|-------|-----|-------|--------|-------|------|

	 1 - Implies this type of variable (row) interacts with this type of
	     variable (column)
	 0 - Implies no interaction
	 
	 Note: A row with all zeros implies that this equation(s) is not being
	       solved.


	Variable Types:
		U1, U2, U3           - Three components of velocity
		T                    - Temperature
		Y                    - Mass fractions   (there can be more than one)
		D1, D2, D3           - Three components of mesh displacement
		D1_RS, D2_RS, D3_RS  - Three components of mesh displacement
		S                    - Surface unknowns (there can be more than one)
		P                    - Pressure
      */

/*
 * Inter_Mask is Defined in mm_unknown_map.c
 */

extern int Inter_Mask[MAX_NUM_MATRICES][MAX_VARIABLE_TYPES][MAX_VARIABLE_TYPES];
extern int Ignore_Deps[MAX_NUM_MATRICES][MAX_VARIABLE_TYPES][MAX_VARIABLE_TYPES];

#endif
