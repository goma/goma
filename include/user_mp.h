/************************************************************************ *
* Goma - Multiphysics finite element software                             *
* Sandia National Laboratories                                            *
*                                                                         *
* Copyright (c) 2022 Goma Developers, National Technology & Engineering   *
*               Solutions of Sandia, LLC (NTESS)                          *
*                                                                         *
* Under the terms of Contract DE-NA0003525, the U.S. Government retains   *
* certain rights in this software.                                        *
*                                                                         *
* This software is distributed under the GNU General Public License.      *
* See LICENSE file.                                                       *
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

EXTERN int usr_thermal_conductivity(dbl *, dbl); /* param - user-defined parameter list       */

EXTERN int usr_electrical_conductivity(dbl *, dbl); /* param - user-defined parameter list       */

EXTERN int usr_density(dbl *); /* param - user-defined parameter list       */

EXTERN int usr_heat_capacity(dbl *, dbl); /* param - user-defined parameter list       */

EXTERN int usr_heat_source(dbl *, dbl); /* param - user-defined parameter list       */

EXTERN int usr_species_source(int,    /* species_no - Current species number       */
                              dbl *); /* param - ptr to user-defined parm list     */

EXTERN int usr_current_source(dbl *); /* param - ptr to user-defined parm list     */

EXTERN int usr_viscosity(dbl *); /* param - ptr to user-defined parm list     */

EXTERN int usr_surface_tension(dbl *); /* param - ptr to user-defined parm list     */

EXTERN int usr_momentum_source(dbl *); /* param - ptr to user-defined parm list     */

EXTERN int usr_lame_mu(struct Elastic_Constitutive *,
                       dbl *); /* param - ptr to user-defined parm list     */

EXTERN int usr_lame_lambda(struct Elastic_Constitutive *,
                           dbl *); /* param - ptr to user-defined parm list     */

EXTERN int usr_expansion(dbl *, /* param - ptr to user-defined parm list     */
                         double *,
                         double[MAX_VARIABLE_TYPES + MAX_CONC],
                         const int);

EXTERN int usr_diffusivity(int,    /* species_no - of diffusivity etc. needed   */
                           dbl *); /* param - ptr to user-defined parm list     */

EXTERN int usr_FlowingLiquidViscosity(dbl *); /* param - ptr to user-defined parm list     */

EXTERN int usr_solid_viscosity(dbl *, /* param - ptr to user-defined parm list     */
                               double *,
                               double[MAX_VARIABLE_TYPES + MAX_CONC],
                               const int);

EXTERN int usr_solid_dil_viscosity(dbl *, /* param - ptr to user-defined parm list     */
                                   double *,
                                   double[MAX_VARIABLE_TYPES + MAX_CONC],
                                   const int);

EXTERN double usr_heat_flux(const double[],   /*   temperature gradient       */
                            double[],         /*   heat flux vector       */
                            double[DIM][DIM], /*   flux sens wrt gradT    */
                            double[DIM][DIM], /*   flow sens  wrt coordinates */
                            const double      /* time */
);

EXTERN int usr_permeability(dbl *); /* param - ptr to user-defined parm list     */

EXTERN int usr_yield_stress(dbl *, dbl); /* param - ptr to user-defined parm list, time*/

#endif /* GOMA_USER_MP_H */
