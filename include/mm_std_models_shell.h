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
 
#ifndef GOMA_MM_STD_MODELS_SHELL_H
#define GOMA_MM_STD_MODELS_SHELL_H

#include "el_elm.h"
#include "mm_shell_bc.h"
#include "std.h"
#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_MM_STD_MODELS_SHELL_C
#define EXTERN
#endif

#ifndef GOMA_MM_STD_MODELS_SHELL_C
#define EXTERN extern
#endif

EXTERN double height_function_model /* mm_std_models_shell.c                */
(dbl *,
       dbl *,
       dbl *,
       dbl *,
       dbl [DIM],
       dbl [DIM],
       dbl *,
       dbl *,
       dbl, dbl);

EXTERN double velocity_function_model /* mm_std_models_shell.c                */
(dbl [DIM],
       dbl [DIM],
       dbl , dbl );
EXTERN double film_evaporation_model /* mm_std_models_shell.c                */
(dbl,
       dbl *,
       dbl ,
       dbl * );

EXTERN double disjoining_pressure_model /* mm_std_models_shell.c                */
(dbl,
       dbl [DIM],
       dbl [DIM],
       dbl [DIM][MDE],
       dbl [DIM][MDE] );

EXTERN double diffusion_coefficient_model /* mm_std_models_shell.c                */
(dbl,
       dbl *);

EXTERN int load_lubrication_momentum_source      /* mm_std_models_shell.c          */
(dbl, dbl);

EXTERN double porous_shell_closed_porosity_model /* mm_std_models_shell.c          */
(void);

EXTERN double porous_shell_closed_radius_model   /* mm_std_models_shell.c          */
(void);

EXTERN double porous_shell_closed_height_model   /* mm_std_models_shell.c          */
(void);

EXTERN double porous_shell_cross_perm_model     /* mm_std_models_shell.c          */
(void);

EXTERN void dynamic_contact_angle_model
(
       double *,
       double *,
       double,
       double *,
       double *
       );

EXTERN double rolling_pressure
(
       double,
       double*,
       double
       );

#endif /* GOMA_MM_STD_MODELS_H */
