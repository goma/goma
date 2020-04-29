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
 
#ifndef GOMA_MM_VISCOSITY_H
#define GOMA_MM_VISCOSITY_H

#include "el_elm.h"
#include "mm_as_structs.h"
#include "mm_mp_structs.h"
#include "mm_qtensor_model.h"
#include "std.h"
#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_MM_VISCOSITY_C
#define EXTERN
#
#endif

#ifndef GOMA_MM_VISCOSITY_C
#define EXTERN extern
#endif

EXTERN double viscosity		/* mm_viscosity.c                            */
(GEN_NEWT_STRUCT *,               /* gn_local                          */
       dbl [DIM][DIM],	        	/* gamma_dot - strain rate tensor    */
       VISCOSITY_DEPENDENCE_STRUCT *); /* d_mu - viscosity dependence       */

EXTERN double power_law_viscosity	/* mm_viscosity.c                    */
(GEN_NEWT_STRUCT *,               /* gn_local                          */
       dbl [DIM][DIM],	        	/* gamma_dot - strain rate tensor    */
       VISCOSITY_DEPENDENCE_STRUCT *); /* d_mu - viscosity dependence       */

EXTERN double herschel_buckley_viscosity/* mm_viscosity.c                    */
(GEN_NEWT_STRUCT *,               /* gn_local                          */
       dbl [DIM][DIM],	        	/* gamma_dot - strain rate tensor    */
       VISCOSITY_DEPENDENCE_STRUCT *); /* d_mu - viscosity dependence       */

EXTERN double carreau_viscosity	/* mm_viscosity.c                            */
(GEN_NEWT_STRUCT *,               /* gn_local                          */
       dbl [DIM][DIM],	        	/* gamma_dot - strain rate tensor    */
       VISCOSITY_DEPENDENCE_STRUCT *); /* d_mu - viscosity dependence       */

EXTERN double bingham_viscosity	/* mm_viscosity.c                            */
(GEN_NEWT_STRUCT *,               /* gn_local                          */
       dbl [DIM][DIM],	        	/* gamma_dot - strain rate tensor    */
       VISCOSITY_DEPENDENCE_STRUCT *); /* d_mu - viscosity dependence       */

EXTERN double bingham_wlf_viscosity	/* mm_viscosity.c                    */
(GEN_NEWT_STRUCT *,               /* gn_local                          */
       dbl [DIM][DIM],	        	/* gamma_dot - strain rate tensor    */
       VISCOSITY_DEPENDENCE_STRUCT *); /* d_mu - viscosity dependence       */

EXTERN double carreau_wlf_viscosity /* mm_viscosity.c                        */
(GEN_NEWT_STRUCT *,               /* gn_local                          */
       dbl [DIM][DIM],	        	/* gamma_dot - strain rate tensor    */
       VISCOSITY_DEPENDENCE_STRUCT *); /* d_mu - viscosity dependence       */

EXTERN int fill_viscosity	/* mm_viscosity.c                            */
(dbl *);			/* param - ptr to user-defined parm list     */

EXTERN int level_set_viscosity	/* mm_viscosity.c                            */
(dbl *);			/* param - ptr to user-defined parm list     */

EXTERN int suspension_viscosity	/* mm_viscosity.c                            */
(int ,			/* species - for solid volume fraction track */
       dbl ,			/* mu0 - carrier fluid viscosity             */
       dbl ,			/* maxpack - maximum solid volume fraction   */       
       dbl ,			/* nexp - exponent for constitutive equation */
       dbl );                /* Concentration field, either in shell or continuum */

EXTERN double carreau_suspension_viscosity	/* mm_viscosity.c                    */
(GEN_NEWT_STRUCT *,               /* gn_local                          */
       dbl [DIM][DIM],	        	/* gamma_dot - strain rate tensor    */
       VISCOSITY_DEPENDENCE_STRUCT *); /* d_mu - viscosity dependence       */

EXTERN double powerlaw_suspension_viscosity /* mm_viscosity.c                   */
(GEN_NEWT_STRUCT *,               /* gn_local                          */
       dbl [DIM][DIM],	        	/* gamma_dot - strain rate tensor    */
       VISCOSITY_DEPENDENCE_STRUCT *); /* d_mu - viscosity dependence       */

EXTERN int epoxy_viscosity	/* mm_viscosity.c                            */
(int ,			/* species - species number for cure eqn     */
       dbl ,			/* mu0 - monomer reference temp viscosity    */
       dbl ,			/* alpha_g - extent of reaction at gel point */
       dbl ,			/* A - exponent for constitutive equation    */
       dbl ,			/* B - exponent for constitutive equation    */
       dbl );			/* Aexp - exponent for thermal viscosity dep */

EXTERN int sylgard_viscosity	/* mm_viscosity.c                            */
(int ,			/* species - species number for cure eqn     */
       dbl ,			/* mu0 - monomer reference temp viscosity    */
       dbl ,			/* alpha_g - extent of reaction at gel point */
       dbl ,			/* A - exponent for constitutive equation    */
       dbl );			/* Aexp - exponent for thermal viscosity dep */

EXTERN int filled_epoxy_viscosity /* mm_viscosity.c                          */
(int ,			/* species_sus - solid volume fraction       */
       int ,			/* species_cur - species num, extent of rxn  */
       dbl ,			/* mu0 - monomer reference viscosity         */
       dbl ,			/* maxpack - maximum solid volume fraction   */
       dbl ,			/* nexp - exponent for constitutive equation */
       dbl ,			/* alpha_g - extent of rxn, gel point        */
       dbl ,			/* A - exponent, cure constitutive equation  */
       dbl ,			/* B - exponent, cure constitutive equation  */
       dbl ,			/* T_g0 - gelation T unreacted filled_epoxy  */
       dbl );			/* Atexp - exponent, thermal viscosity func  */

EXTERN int foam_epoxy_viscosity /* mm_viscosity.c                          */
(int ,			/* species_sus - solid volume fraction       */
       int ,			/* species_cur - species num, extent of rxn  */
       dbl ,			/* mu0 - monomer reference viscosity         */
       dbl ,			/* alpha_g - extent of rxn, gel point        */
       dbl );			/* Aexp - exponent, thermal viscosity func  */

EXTERN int
foam_pmdi10_viscosity(int species,    /* species number for cure equation */
		      dbl mu0,        /* monomer reference temperature viscosity */
		      dbl alpha_g,    /* extent of reaction at the gel point */
		      dbl A,          /* exponent for constitutive equation */
		      dbl B,          /* exponent for constitutive equation */
		      dbl norm_E);     /* Normalized activation energy */

EXTERN int thermal_viscosity	/* mm_viscosity.c                            */
(dbl ,			/* mu0 - ref temperature fluid viscosity     */
       dbl );			/* Aexp - exponent for constitutive equation */

EXTERN int bond_viscosity
(dbl ,          /* reference zero shear rate fluid viscosity */
       dbl ,          /* reference high shear rate fluid viscosity */
       dbl );        /* exponent for constitutive equation */

EXTERN int cure_viscosity	/* mm_viscosity.c                            */
(int ,			/* species - num, solid volume fraction      */
       dbl ,			/* mu0 - carrier fluid viscosity             */
       dbl ,			/* alpha_g - extent of reaction at gel point */
       dbl ,			/* A - exponent for constitutive equation    */
       dbl );			/* B - exponent for constitutive equation    */

EXTERN double carreau_wlf_conc_viscosity /* mm_viscosity.c                   */
(GEN_NEWT_STRUCT *,               /* gn_local                          */
       dbl [DIM][DIM],	        	/* gamma_dot - strain rate tensor    */
       VISCOSITY_DEPENDENCE_STRUCT *,   /* d_mu - viscosity dependence       */
       const int ); 			/* const. model      */


EXTERN int ls_modulate_viscosity
( double *,                      /* Primary phase viscosity value */
	  double  ,                      /* Second phase viscosity */
	  double width,                  /* Transition width */
	  double pm_minus,               /* Negative phase viscosity mask */
	  double pm_plus,                /* Positive phase viscosity mask */
        VISCOSITY_DEPENDENCE_STRUCT *, /* d_mu - viscosity dependence */
        const int model);             /* Viscosity Model - CONSTANT or RATIO*/
	  
EXTERN void copy_pF_to_F 
( int  );

EXTERN double flowing_liquid_viscosity		/* mm_viscosity.c                            */
(VISCOSITY_DEPENDENCE_STRUCT *); /* d_flow_vis - flowing liquid viscosity sensitivity   */

#endif /* GOMA_MM_VISCOSITY_H */
