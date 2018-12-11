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
 
#ifndef _MM_STD_MODELS_H
#define _MM_STD_MODELS_H

#ifdef EXTERN
#undef EXTERN
#endif

#ifdef _MM_STD_MODELS_C
#define EXTERN
#endif

#ifndef _MM_STD_MODELS_C
#define EXTERN extern
#endif

EXTERN int bouss_momentum_source /* mm_std_models.c                          */
PROTO((dbl [DIM],		/* f - Body force.                           */
       MOMENTUM_SOURCE_DEPENDENCE_STRUCT *,
       int ,			/* jxb - on/off Flag for Lorentz force       */
       int ));			/* hydrostatic - Boolean for including       *
				 * hydrostatic head in thermal buoyancy term */
EXTERN int jxb_mesh_source /* mm_std_models.c                          */
PROTO((dbl [DIM]));		/* f - Body force.                           */

EXTERN int EHD_POLARIZATION_source /*mm_std_models.c                         */
PROTO((dbl [DIM],               /* f - Body force                            */
       MOMENTUM_SOURCE_DEPENDENCE_STRUCT * ));

EXTERN int gravity_vibrational_source /* mm_std_models.c                     */
PROTO((dbl f[DIM],              /* Body force in direction [a]               */
       MOMENTUM_SOURCE_DEPENDENCE_STRUCT *df, /* df                          */
       dbl time));              /* time                                      */

EXTERN int suspend_momentum_source /* mm_std_models.c                        */
PROTO((dbl [DIM],		/* f - Body force.                           */
       MOMENTUM_SOURCE_DEPENDENCE_STRUCT * ));

EXTERN int fill_momentum_source	/* mm_std_models.c                           */
PROTO((double [DIM]));		/* f                                         */

EXTERN int epoxy_dea_species_source /* mm_std_models.c                       */
PROTO((int ,			/* species_no - Current species number       */
       double *));		/* param - pointer to user-defined parm list */

EXTERN int epoxy_species_source /* mm_std_models.c                           */
PROTO((int ,			/* species_no - Current species number       */
       double *));		/* param - user-defined parameter list       */

EXTERN double epoxy_heat_source	/* mm_std_models.c                           */
PROTO((HEAT_SOURCE_DEPENDENCE_STRUCT *,
       double ,			/* tt - parameter to vary time integration   *
				 * from explicit (tt = 1) to                 *
				 *      implicit (tt = 0)                    */
       double ));		/* dt - current time step size               */

EXTERN double joule_heat_source	/* mm_std_models.c                           */
PROTO((HEAT_SOURCE_DEPENDENCE_STRUCT *, dbl));

EXTERN double butler_volmer_heat_source       /* mm_std_models.c             */
PROTO((HEAT_SOURCE_DEPENDENCE_STRUCT *, dbl *));

EXTERN double butler_volmer_source            /* mm_std_models.c             */
PROTO((dbl *, int, dbl *));

EXTERN double visc_diss_heat_source /* mm_std_models.c                       */
PROTO((HEAT_SOURCE_DEPENDENCE_STRUCT *,
       dbl *));			/* param - General multipliers               */

EXTERN double enthalpy_heat_capacity_model /* mm_std_models.c                */
PROTO((HEAT_CAPACITY_DEPENDENCE_STRUCT *));

EXTERN int V_mesh_sfs_model	/* mm_std_models.c                           */
PROTO((dbl *,			/* param - Body force.                       */
       dbl *,			/* v_mesh_sfs - three rotational component   */
       const int,		/* v_mesh_sfs_model - rotational	     */
       int ));			/* gnn - -1 for gauss point, global node     *
				 * number for node point                     */
EXTERN int Diffusivity		/* mm_std_models.c                           */
PROTO((void ));		
      
EXTERN int Free_Vol_Theory_Diffusivity /* mm_std_models.c                    */
PROTO((int ,			/* species_no - current species number       */
       double *));		/* param - free volume params from mat file  */
   
EXTERN int hydro_flux		/* mm_std_models.c                           */
PROTO((struct Species_Conservation_Terms *, /* st                            */
       int ,			/* w - species number                        */
       dbl , 			/* time step parameter                       */
       dbl ,                    /* time step                                 */
       const dbl [DIM]));       /* element size                              */

EXTERN int suspension_balance	/* mm_std_models.c                           */
PROTO((struct Species_Conservation_Terms *, /* st                            */
       int ));			/* w - species number                        */

EXTERN int particle_stress      /* mm_std_models.c                           */
PROTO((dbl [DIM][DIM],          /* particle stress */
       dbl [DIM][DIM][DIM][MDE],/* derivative wrt velocity  */
       dbl [DIM][DIM][DIM][MDE],/* derivative wrt vorticity direction */
       dbl [DIM][DIM][MAX_CONC][MDE],/* derivative wrt concentration */
       dbl [DIM][DIM][DIM][MDE],/* derivative wrt mesh */
       dbl [DIM][DIM][MDE],     /* derivative wrt pressure */
       int ));                  /* species number                            */

EXTERN int divergence_particle_stress /* mm_std_models.c                     */
PROTO((dbl [DIM],               /* divergence of the stress*/
       dbl [DIM][MDE],          /* derivative wrt shear rate invariant */
       dbl [DIM][MAX_CONC][MDE],/* derivative wrt concentration */
       dbl [DIM][DIM][MDE],     /* derivative wrt velocity */
       dbl [DIM][DIM][MDE],     /* derivative wrt mesh */
       dbl [DIM][DIM][MDE],     /* derivative wrt vorticity direction */
       dbl [DIM][MDE],          /* derivative wrt pressure */
       int ));	                /* species number */

EXTERN int antoine_psat		/* mm_std_models.c                           */
PROTO((int ,			/* species_no                                */
       double [],		/* param                                     */
       double *,		/* f                                         */
       double *));		/* dfdt                                      */

EXTERN int riedel_psat		/* mm_std_models.c                           */
PROTO((int ,			/* species_no                                */
       double [],		/* param                                     */
       double *,		/* f                                         */
       double *));		/* dfdt                                      */

EXTERN int suspension_pm_fluid_momentum_source /* mm_std_models.c            */
PROTO((dbl [DIM],		/* f - Body force.                           */
       MOMENTUM_SOURCE_DEPENDENCE_STRUCT *df ));

EXTERN int suspension_pm_particle_momentum_source /* mm_std_models.c         */
PROTO((dbl [DIM],		/* f - Body force.                           */
       dbl [DIM][MDE],		/* dfdT - For temperature dependence.        */
       dbl [DIM][DIM][MDE],	/* dfdX - For spatial dependence.            */
       dbl [DIM][MAX_CONC][MDE], /* dfdC - For concentration dependence.     */
       dbl [DIM][DIM][MDE]));	/* dfdv - For velocity dependence.           */

EXTERN int molten_glass_viscosity /* mm_std_models.c                         */
PROTO((dbl *,			/* vis - Base FLOWING LIQUID VISCOSITY       */
       dbl [MDE],		/* dvis_dT - temperature dependence.         */
       dbl *));			/* param  - parameter list                   */

EXTERN int epoxy_flowing_liquid_viscosity /* mm_std_models.c                     */
PROTO((dbl *,			          /* vis - Base FLOWING LIQUID VISCOSITY */
       VISCOSITY_DEPENDENCE_STRUCT *,     /* vis sensitivity.                    */
       dbl *));			          /* param  - parameter list             */

EXTERN int electrode_species_source /* mm_std_models.c                       */
PROTO((int ,			/* species_no - Current species number       */
       double ,			/* time - present time value                 */
       double ));		/* delta_t - present time step               */

EXTERN int ion_reaction_source /*  mm_std_models.c -- RSL 3/20/01  */
PROTO((int));                  /*  current species number          */

EXTERN int electrolyte_temperature /* mm_std_models.c                        */
PROTO((double ,			/* t - present value of time                 */
       double ,			/* delta_t - present time step               */
       int ));			/* print - key for printing:                 *
				 * print = 1: printing;                      *
				 * print = 0: no printing.                   */

EXTERN int assemble_suspension_temperature /* mm_std_models.c */
PROTO((dbl,                     /* time */
       dbl,                     /* theta */
       dbl));                   /* delta_t */

EXTERN double vary_rho_heat_source /* mm_std_models.c */
PROTO((HEAT_SOURCE_DEPENDENCE_STRUCT *,
       double ,			/* tt - parameter varies time integration from 
				 * explicit (tt = 1) to implicit (tt = 0) */
       double ));		/* dt - current time step size */
	   
EXTERN double foam_heat_source
PROTO((HEAT_SOURCE_DEPENDENCE_STRUCT *,
	   double ,	/* parameter to vary time integration from explicit (tt = 1) to implicit (tt = 0) */
	   double ));	/* current time step size */
						

EXTERN int Generalized_Diffusivity /* mm_std_models.c */
PROTO(( void ));

EXTERN int Generalized_FV_Diffusivity /* mm_std_models.c */
PROTO((int));			/* species number */

EXTERN int foam_species_source /* mm_std_models.c                       */
PROTO(( double *));		/* param - pointer to user-defined parm list */

EXTERN int foam_epoxy_species_source /* mm_std_models.c                       */
PROTO(( int, 
	double *,		/* param - pointer to user-defined parm list */ 
	double,                 /* tt */
	double));               /* dt */

EXTERN int assemble_bond_evolution /* mm_std_models.c */
PROTO((dbl,                     /* Current time */
       double ,			/* tt - parameter varies time integration from 
				 * explicit (tt = 1) to implicit (tt = 0) */
       double ));		/* dt - current time step size */


#endif /* _MM_STD_MODELS_H */
