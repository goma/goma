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
 
#ifndef GOMA_BC_SURFACEDOMAIN_H
#define GOMA_BC_SURFACEDOMAIN_H

#include "el_elm.h"
#include "mm_eh.h"
#include "rf_fem_const.h"
#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_BC_SURFACEDOMAIN_C
#define EXTERN
#
#endif

#ifndef GOMA_BC_SURFACEDOMAIN_C
#define EXTERN extern
#endif

EXTERN void mass_flux_sd_bc
(double [],		/* func                                      */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
       int ,			/* wspec - species number this BC            */
       double ,			/* mass_tran_coeff - (cgs?? MKS units)       */
       double ,			/* Y_c - bath concentration 	             */
       double ,			/* dt - current value of the time step       */
       double );		/* tt - parameter varies time integration    *
				 * from explicit to implicit                 */

#endif /* GOMA_BC_SURFACEDOMAIN_H */
