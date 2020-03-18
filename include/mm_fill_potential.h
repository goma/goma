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
 * mm_fill_potential.h -- prototype declarations for mm_fill_potential.c
 */

#ifndef GOMA_MM_FILL_POTENTIAL_H
#define GOMA_MM_FILL_POTENTIAL_H

#include "el_elm.h"
#include "mm_fill_porous.h"
#include "rf_fem_const.h"

struct Boundary_Condition;
#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_MM_FILL_POTENTIAL_C
#define EXTERN /* do nothing */
#endif

#ifndef GOMA_MM_FILL_POTENTIAL_C
#define EXTERN extern
#endif

EXTERN int assemble_Enorm
(void);

EXTERN int assemble_potential
(double ,			/* time - present time value                 */
       double ,			/* tt - parameter to vary time integration 
				 * from explicit (tt = 1) to 
				 * implicit (tt = 0)                         */
       double );		/* dt - current time step size               */

EXTERN void current_BV_surf	/* mm_fill_potential.c                       */
(double [],		/* func                                      */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
       int ,			/* wspec - species number of this BC         */
       double ,			/* nu - stoichiometric coefficient           */
       double ,			/* k - kinetic rate constant                 */
       double ,			/* beta - reaction order    	             */
       double ,			/* alphaa - anodic transfer coefficient      */
       double ,			/* alphac - cathodic transfer coefficient    */
       double ,			/* V - electrode potential                   */
       double ,			/* U0 - electrolyte open-circuit potential   */
       double ,			/* dt - current value of the time step       */
       double );      		/* tt - parameter to vary time integration   */


EXTERN void current_ORR_surf    /* mm_fill_potential.c                       */
(double [],               /* func                                      */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
       int ,          /* wspec - species number of this boundary condition   */
       double ,       /* ai0 - product of interfacial area by                */
                      /* exchange current density (A/cm^3)                   */
       double ,       /* H - thickness of catalyst layer (cm)                */
       double ,       /* cref - species ref. concentration (moles/cm^3)      */
       double ,       /* alphac - cathodic direction transfer coefficient    */
       double ,       /* T - cell temperature (K)                            */
       double ,       /* V - call voltage (V)                                */
       double ,       /* U0 - open-circuit potential (V)                     */
       double ,       /* beta - reaction order                               */
       double ,       /* dt - current value of the time step                 */
       double );     /* tt - parameter to vary time integration             */

EXTERN void current_HOR_surf    /* mm_fill_potential.c                       */
(double [],               /* func                                      */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
       int ,          /* wspec - species number of this boundary condition   */
       double ,       /* ai0 - product of interfacial area by                */
                      /* exchange current density (A/cm^3)                   */
       double ,       /* H - thickness of catalyst layer (cm)                */
       double ,       /* cref - species ref. concentration (moles/cm^3)      */
       double ,       /* alphaa - anodic direction transfer coefficient      */
       double ,       /* alphac - cathodic direction transfer coefficient    */
       double ,       /* T - cell temperature (K)                            */
       double ,       /* U0 - open-circuit potential for HOR (V)             */
       double ,       /* beta - reaction order                               */
       double ,       /* V - electrode potential (V)                         */
       double ,       /* dt - current value of the time step                 */
       double );     /* tt - parameter to vary time integration             */

EXTERN void surface_charge_surf
(double *,
       double [DIM][MAX_VARIABLE_TYPES+MAX_CONC][MDE],
       const double);      /* sigma - value of surface charge, 0 indicates 
			      electroneutrality */

EXTERN int assemble_electric_field
(void);      /* Least square equation of Efield = grad voltage */

EXTERN void apply_potential_grad_bc
(double [],               /* func                                      */
       double [][MAX_VARIABLE_TYPES+MAX_CONC][MDE], /* d_func                */
       struct Boundary_Condition *,/* BC_Type                             */
       const double,		/* time_value	*/
       const double );                  /* time_step */

EXTERN double electrical_conductivity
(double [MDE],               /* dkdT                         */
       double [MDE],			/*dkdV			*/
       double [MAX_CONC][MDE],		/*dkdC			*/
       double [DIM][MDE],		/*dkdX			*/
       const double,		/* time	*/
       const double );                  /* time_step */

#endif /* GOMA_MM_FILL_POTENTIAL_H */
