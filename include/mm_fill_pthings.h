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
 
#ifndef GOMA_MM_FILL_PTHINGS_H
#define GOMA_MM_FILL_PTHINGS_H

#include "el_elm.h"
#include "mm_fill_potential.h"
#include "std.h"
#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_MM_FILL_PTHINGS_C
#define EXTERN
#
#endif

#ifndef GOMA_MM_FILL_PTHINGS_C
#define EXTERN extern
#endif

EXTERN int assemble_pmomentum
(	double ,		/* time - present time value                 */
	dbl ,			/* tt - parameter to vary time integration   *
				 * from explicit (tt = 1) to                 *
				 * implicit (tt = 0)                         */
	dbl );			/* dt - current time step size               */

EXTERN int MMH_assemble_continuity
(double,                  /* time - present time value                 */
       double,			/* parameter to vary time integration from 
				 * explicit (tt = 1) to implicit (tt = 0)    */
       double,			/* current time step size */
       double,			/* average global element size for PSPG,
				 * taken to be constant wrt to Jacobian 
				 * entries                                   */
       double [],		/* h - element size information for PSPG     */
       double [][DIM],		/* hh - (DIM)(DIM) currently unused, but     *
				 * remain just in case they are needed later */
       double [][MDE],		/* dh_dxnode - (DIM)(MDE)                    */
       double ,			/* U_norm - global velocity norm for PSPG    *
				 * calculations                              */
       double );		/* mu_avg - element viscosity for PSPG calcs */

EXTERN int pmomentum_source_term
( dbl [DIM],                   /* Body force */
        dbl [DIM][MDE],              /* For temperature dependence */
        dbl [DIM][DIM][MDE],         /* For spatial dependence */
        dbl [DIM][MAX_CONC][MDE],    /* For concentration */
        dbl [DIM][DIM][MDE] );      /* For velocity dependence */

#endif /* GOMA_MM_FILL_PTHINGS_H */
