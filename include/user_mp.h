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
 

#ifndef GOMA_USER_MP_H
#define GOMA_USER_MP_H

#include "el_elm.h"
#include "rf_fem_const.h"
#include "std.h"
#include "user_ac.h"

struct Elastic_Constitutive;
#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_USER_MP_C
#define EXTERN
#
#endif

#ifndef GOMA_USER_MP_C
#define EXTERN extern
#endif

EXTERN int usr_thermal_conductivity
(dbl *, dbl);			/* param - user-defined parameter list       */

EXTERN int usr_electrical_conductivity
(dbl *, dbl);			/* param - user-defined parameter list       */

EXTERN int usr_density
(dbl *);			/* param - user-defined parameter list       */

EXTERN int usr_heat_capacity
(dbl *, dbl);			/* param - user-defined parameter list       */

EXTERN int usr_heat_source
(dbl *, dbl);			/* param - user-defined parameter list       */

EXTERN int usr_species_source
(int ,			/* species_no - Current species number       */
       dbl *);			/* param - ptr to user-defined parm list     */

EXTERN int usr_current_source
(dbl *);			/* param - ptr to user-defined parm list     */

EXTERN int usr_viscosity
(dbl *);			/* param - ptr to user-defined parm list     */

EXTERN int usr_surface_tension
(dbl *);			/* param - ptr to user-defined parm list     */

EXTERN int usr_momentum_source
(dbl *);			/* param - ptr to user-defined parm list     */

EXTERN int usr_lame_mu
(struct Elastic_Constitutive *, dbl *); /* param - ptr to user-defined parm list     */

EXTERN int usr_lame_lambda
(struct Elastic_Constitutive *, dbl *); /* param - ptr to user-defined parm list     */

EXTERN int usr_expansion
(dbl *,			/* param - ptr to user-defined parm list     */
	double *,
	double [MAX_VARIABLE_TYPES+MAX_CONC]
	);	

EXTERN int usr_diffusivity
(int ,			/* species_no - of diffusivity etc. needed   */
       dbl *);			/* param - ptr to user-defined parm list     */

EXTERN int usr_FlowingLiquidViscosity
(dbl *);			/* param - ptr to user-defined parm list     */

#if defined SECOR_HEAT_FLUX 
EXTERN double  usr_heat_flux
( const double [],        /*   temperature gradient       */
        double [],                      /*   heat flux vector       */
        double [][DIM],                 /*   flux sens wrt gradT    */
        double [][DIM],              /*   flow sens  wrt coordinates */
	const double,			/* time */
        const double,                   /* gap  */
        const double [DIM],             /* dgap_dX */
        const double [DIM],             /* Velocity bottom */
        const double [DIM],             /* Velocity top  */
        double [][DIM],                 /*  dq_dVb  */
        double [][DIM] );              /*  dq_dVt  */
#else
EXTERN double  usr_heat_flux
( const double [],        /*   temperature gradient       */
        double [],                      /*   heat flux vector       */
        double [][DIM],                 /*   flux sens wrt gradT    */
        double [][DIM],              /*   flow sens  wrt coordinates */
        const double                   /* time */
     );
#endif

EXTERN int usr_permeability
(dbl *);			/* param - ptr to user-defined parm list     */

EXTERN int usr_yield_stress
(dbl *, dbl);		/* param - ptr to user-defined parm list, time*/

#endif /* GOMA_USER_MP_H */
