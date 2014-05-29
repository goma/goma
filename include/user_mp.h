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
 

#ifndef _USER_MP_H
#define _USER_MP_H

#ifdef EXTERN
#undef EXTERN
#endif

#ifdef _USER_MP_C
#define EXTERN
#
#endif

#ifndef _USER_MP_C
#define EXTERN extern
#endif

EXTERN int usr_thermal_conductivity
PROTO((dbl *, dbl));			/* param - user-defined parameter list       */

EXTERN int usr_electrical_conductivity
PROTO((dbl *, dbl));			/* param - user-defined parameter list       */

EXTERN int usr_density
PROTO((dbl *));			/* param - user-defined parameter list       */

EXTERN int usr_heat_capacity
PROTO((dbl *, dbl));			/* param - user-defined parameter list       */

EXTERN int usr_heat_source
PROTO((dbl *, dbl));			/* param - user-defined parameter list       */

EXTERN int usr_species_source
PROTO((int ,			/* species_no - Current species number       */
       dbl *));			/* param - ptr to user-defined parm list     */

EXTERN int usr_current_source
PROTO((dbl *));			/* param - ptr to user-defined parm list     */

EXTERN int usr_viscosity
PROTO((dbl *));			/* param - ptr to user-defined parm list     */

EXTERN int usr_surface_tension
PROTO((dbl *));			/* param - ptr to user-defined parm list     */

EXTERN int usr_momentum_source
PROTO((dbl *));			/* param - ptr to user-defined parm list     */

EXTERN int usr_lame_mu
PROTO((struct Elastic_Constitutive *, dbl *)); /* param - ptr to user-defined parm list     */

EXTERN int usr_lame_lambda
PROTO((struct Elastic_Constitutive *, dbl *)); /* param - ptr to user-defined parm list     */

EXTERN int usr_expansion
PROTO((dbl *,			/* param - ptr to user-defined parm list     */
	double *,
	double [MAX_VARIABLE_TYPES+MAX_CONC]
	));	

EXTERN int usr_diffusivity
PROTO((int ,			/* species_no - of diffusivity etc. needed   */
       dbl *));			/* param - ptr to user-defined parm list     */

EXTERN int usr_FlowingLiquidViscosity
PROTO((dbl *));			/* param - ptr to user-defined parm list     */

#if defined SECOR_HEAT_FLUX 
EXTERN double  usr_heat_flux
PROTO(( const double [],        /*   temperature gradient       */
        double [],                      /*   heat flux vector       */
        double [][DIM],                 /*   flux sens wrt gradT    */
        double [][DIM],              /*   flow sens  wrt coordinates */
	const double,			/* time */
        const double,                   /* gap  */
        const double [DIM],             /* dgap_dX */
        const double [DIM],             /* Velocity bottom */
        const double [DIM],             /* Velocity top  */
        double [][DIM],                 /*  dq_dVb  */
        double [][DIM] ));              /*  dq_dVt  */
#else
EXTERN double  usr_heat_flux
PROTO(( const double [],        /*   temperature gradient       */
        double [],                      /*   heat flux vector       */
        double [][DIM],                 /*   flux sens wrt gradT    */
        double [][DIM],              /*   flow sens  wrt coordinates */
        const double                   /* time */
     ));
#endif

EXTERN int usr_permeability
PROTO((dbl *));			/* param - ptr to user-defined parm list     */

EXTERN int usr_yield_stress
PROTO((dbl *, dbl));		/* param - ptr to user-defined parm list, time*/

#endif /* _USER_MP_H */
