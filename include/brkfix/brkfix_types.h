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

#ifndef _BRKFIX_TYPES_H
#define _BRKFIX_TYPES_H

/*
 * Multiplicity of each Basic equation and variable is more than one for
 * any/none/all of 3 reasons:
 *
 *		    (i) It is a vector or tensor equation with
 *			2,3,6, or 9 degrees of freedom instead of just 1.
 *
 *		   (ii) It is a concentration like variable with multiplicity
 *			equal to the number of species in the problem.
 *
 *		  (iii) It is a special interpolation, like piecewise linear
 *			discontinuous pressure that has 3 degrees of freedom
 *			associated with a single node when it appears.
 */

typedef struct Basic_Equation_Variable_Multiplicity
{
  int element_block_id;
  int eqnvar_id;
  int vect_mult;
  int conc_mult;
  int ndof_mult;
} Bevm;

#endif
