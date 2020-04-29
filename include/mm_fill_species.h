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
 *$Id: mm_fill_species.h,v 5.5 2008-11-06 15:53:44 hkmoffa Exp $
 */

#ifndef GOMA_MM_FILL_SPECIES_H
#define GOMA_MM_FILL_SPECIES_H

#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_MM_FILL_SPECIES_C
#define EXTERN /* do nothing */
#endif

#ifndef GOMA_MM_FILL_SPECIES_C
#define EXTERN extern
#endif

#include "mm_as_structs.h"
#include "el_elm.h"
#include "mm_fill_population.h"
#include "rf_bc_const.h"
#include "rf_fem_const.h"
#include "std.h"

EXTERN int assemble_mass_transport /* mm_fill_species.c                      */
    (double time, /* time - present time valuel; KSC           */
     double tt,   /* tt - parm to vary time integration from
                   * explicit (tt = 1) to implicit (tt = 0)    */
     double dt,   /* dt - current time step size               */
     PG_DATA *pg_data);

EXTERN int assemble_mass_transport_path_dependence
(double ,                 /* time - present time valuel; KSC           */
       double ,                 /* tt - parameter to vary time integration from
				 * explicit (tt = 1) to implicit (tt = 0)    */
       double ,                 /* dt - current time step size               */
       const dbl [DIM],         /* h - element sizes, not scale factors.     */
       const dbl [DIM][DIM],    /* hh - */
       const dbl [DIM][MDE],    /* dh_dxnode -                               */
       const dbl [DIM],         /* vcent - average element velocity, which is
				 * the centroid velocity for Q2 and
				 * the average of the vertices for Q1.
				 * From routine "element_velocity."          */
       const dbl [DIM][MDE]);  /* dvc_dnode                                 */


EXTERN void mass_flux_surf_mtc  /* mm_fill_species.c                         */
(dbl [MAX_CONC],		/* mass_flux                                 */
       dbl [MAX_CONC][MAX_VARIABLE_TYPES+MAX_CONC], /* d_mass_flux           */
       double ,			/* Tsurf_quad - temperature at the surface   */
       double [MAX_CONC],	/* Ysurf_quad - concentration at surface 
				 * quadrature point                          */
       int ,			/* wspec -                                   */
       double ,			/* mass_tran_coeff - Heat transfer coefficient 
				 * (cgs?? MKS units)                         */
       double );		/* Y_c - bath temperature (Kelvin)	     */



EXTERN void mass_flux_surf_BV
(dbl mass_flux[MAX_CONC],
       dbl d_mass_flux[MAX_CONC][MAX_VARIABLE_TYPES+MAX_CONC],
       int wspec,                    /* species number                     */
       double nu,                    /* stoichiometric coefficient         */  
       double k,                     /* kinetic rate constant              */  
       double beta,                  /* reaction order                     */
       double alphaa,                /* anodic direction transfer coef.    */
       double alphac,                /* cathodic direction transfer coef.  */
       double V,                     /* electrode potential                */
       double U0,                    /* open circuit electrolyte potential */
       double T);	             /* electrolyte solution temperature   */

EXTERN void mass_flux_surf_HOR
(dbl mass_flux[MAX_CONC],
       dbl d_mass_flux[MAX_CONC][MAX_VARIABLE_TYPES+MAX_CONC],
       int wspec,      /* species number                          */
       double ai0,     /* product of interfacial area by          */
                       /* exchange current density (A/cm^3)       */
       double H,       /* thickness of catalyst layer (cm)        */
       double cref,    /* ref. concentration (moles/cm^3)         */
       double alphaa,  /* anodic direction transfer coefficient   */
       double alphac,  /* cathodic direction transfer coefficient */
       double T,       /* cell temperature (K)                    */
       double U0,      /* open-circuit potential (V)              */
       double beta,    /* reaction order                          */
       double n,       /* number of electrons involved in rxn     */
       double V);      /* electrode potential                     */

EXTERN void mass_flux_surf_ORR
(dbl mass_flux[MAX_CONC],
       dbl d_mass_flux[MAX_CONC][MAX_VARIABLE_TYPES+MAX_CONC],
       int wspec,      /* species number                            */
       double ai0,     /* product of interfacial area by            */
       /* exchange current density (A/cm^3)         */
       double H,       /* thickness of catalyst layer (cm)          */
       double cref,    /* species ref. concentration (moles/cm^3)   */
       double alphac,  /* cathodic direction transfer coefficient   */
       double T,       /* cell temperature (K)                      */
       double V,       /* cel voltage (V)                           */
       double U0,      /* open-circuit potential (V)                */
       double beta,    /* reaction order                            */
       double n);     /* number of electrons involved in rxn       */

EXTERN void mass_flux_surf_H2O_ANODE
(dbl mass_flux[MAX_CONC],
       dbl d_mass_flux[MAX_CONC][MAX_VARIABLE_TYPES+MAX_CONC],
       int wspec,      /* species number                          */
       double ai0,     /* product of interfacial area by          */
       /* exchange current density (A/cm^3)       */
       double Ha,      /* thickness of anode catalyst layer (cm)  */
       double cH2ref,  /* ref. concentration for H2 (moles/cm^3)  */
       double alphaa,  /* anodic direction transfer coefficient   */
       double alphac,  /* cathodic direction transfer coefficient */
       double T,       /* cell temperature (K)                    */
       double U0a,     /* Open-circuit potential for HOR (V)      */
       double nd);    /* electro-osmatic drag coefficient        */

EXTERN void mass_flux_surf_H2O_CATHODE
(dbl mass_flux[MAX_CONC],
       dbl d_mass_flux[MAX_CONC][MAX_VARIABLE_TYPES+MAX_CONC],
       int wspec,      /* species number                           */
       double ai0,     /* product of interfacial area by           */
       /* exchange current density (A/cm^3)        */
       double Hc,      /* thickness of cathode catalyst layer (cm) */
       double cO2ref,  /* ref. concentration for O2 (moles/cm^3)   */
       double alphac,  /* cathodic direction transfer coefficient  */
       double T,       /* cell temperature (K)                     */
       double V,       /* cel voltage (V)                          */
       double U0c,     /* Open-circuit potential for ORR (V)       */
       double nd);    /* electro-osmatic drag coefficient         */

EXTERN void mass_flux_surf_SULFIDATION 
(dbl mass_flux[MAX_CONC],
       dbl d_mass_flux[MAX_CONC][MAX_VARIABLE_TYPES+MAX_CONC],
       int mode,      /* key word for sulfidation kinetic model */
       int wspec,     /* species number             */
       double nu,     /* stoichiometric coefficient */  
       double k1,     /* forward rate constant      */  
       double E1,     /* forward activation energy  */
       double kn1,    /* backward rate constant     */  
       double En1,    /* backward activation energy */
       double T,      /* Temperature                */
       double c_H2S,  /* bulk concentration of H2S  */
       double c_O2); /* bulk concentration of O2   */

EXTERN void mass_flux_surf_BV2    /* mm_fill_species.c - RSL 5/28/02         */
(dbl ,              /* time                                            */
       dbl [MAX_CONC],    /* mass_flux                                       */
       dbl [MAX_CONC][MAX_VARIABLE_TYPES+MAX_CONC], /* d_mass_flux           */
       int ,              /* wspec -- species number                         */
       int ,              /* flag -- output option                           */
       double ,           /* k -- rate constant                              */
       double ,           /* PHI_E -- electrode potential                    */
       double ,           /* alphaa -- anodic transfer coefficient           */
       double ,           /* alphac -- cathodic transfer coefficient         */
       double ,           /* U0 -- standard open circuit potential           */
       double );         /* T -- electrolyte temperature                    */

EXTERN void mass_flux_surf_NI    /* mm_fill_species.c - RSL 5/28/02    */
(dbl [MAX_CONC],           /* mass_flux                          */
       dbl [MAX_CONC][MAX_VARIABLE_TYPES+MAX_CONC], /* d_mass_flux     */
       dbl ,              /* time                                      */
       int ,              /* wspec -- species number                   */
       int ,              /* flag -- output option                     */
       double ,           /* PHI_E -- electrode potential              */
       double );         /* T -- electrolyte temperature              */

EXTERN void mass_flux_surf_etch    /* mm_fill_species.c    */
(dbl [DIM],         /* func                                      */
       dbl [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func        */
       const int,         /* wspec -- species number                   */
       const int,         /* etching plane                             */
       const double,      /* time                                      */
       const double,      /* dt, current value of the time step        */
       const double );   /* tt, parameter to vary time integration    */

EXTERN void mass_flux_surf_const  /* mm_fill_species.c - RSL 6/28/00   */
(dbl [MAX_CONC],    /* mass_flux                                 */
       dbl [MAX_CONC][MAX_VARIABLE_TYPES+MAX_CONC], /* d_mass_flux     */
       double ,           /* Tsurf_quad - temperature at the surface   */
       double [MAX_CONC], /* Ysurf_quad - concentration at surface     *
                           * quadrature point                          */
       int ,              /* wspec                                     */
       double );         /* const_mass_flux                           */

EXTERN void mass_flux_surf_user_bc
(double [DIM],	        /*  func                                     */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
       const int ,		/* wspec - species number                    */
       const double [],		/* p - parameterize species transfer model   */
       const dbl );		/* time  */

EXTERN void mass_flux_alloy_surf
(double [DIM],	        /*  func                                     */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
       const int ,		/* wspec - species number                    */
       const double [],		/* p - parameterize species transfer model   */
       const dbl );		/* time  */


EXTERN void yflux_disc_rxn_bc   /* mm_fill_species.c                         */
(double [],		/* func                                      */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
       int ,			/* wspec                                     */
       int ,			/* i_mat_a                                   */
       int ,			/* i_mat_b                                   */
       double ,			/* kf, forward rate constant                 */
       double ,			/* kr, forward rate constant                 */
       double ,			/* dt, current value of the time step        */
       double );		/* tt, parameter to vary time integration    */

EXTERN void raoults_law         /* mm_fill_species.c                         */
(double [],		/* func                                      */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
       int ,			/* wspec                                     */
       int ,			/* i_mat_liq                                 */
       int ,			/* i_mat_gas                                 */
       double ,			/* amb_pres - ambient pressure               */
       double ,			/* M1                                        */
       double ,			/* M2                                        */
       double ,			/* M3                                        */
       double );		/* M4 molecular weights of 4 species. M1 must
				 * be the first volatile species, M2 is the
                                 * second volatile species, M3 is the condensed
                                 * phase in the liquid, and M4 is the insoluble
                                 * gas phase, e.g. air                       */

extern void raoults_law_new(double [], 
			    double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
			    BOUNDARY_CONDITION_STRUCT *);

EXTERN void flory_huggins	/* mm_fill_species.c                         */
(double [],		/* func                                      */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
       int ,			/* wspec                                     */
       int ,			/* i_mat_gas                                 */
       int ,			/* i_mat_liq - mat id for gas and liquid     */
       int ,			/* mode VOLUME or MASS based formulation     */
       double );		/* amb_pres - ambient pressure               */

EXTERN void kinematic_species_bc /* mm_fill_species.c                        */
(double [DIM],		/* func                                      */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
       int ,			/* wspec current species no.                 */
       double ,			/* vnormal - normal velocity                 */
       double [MAX_PDIM],	/* x_dot - mesh velocity vector              */
       dbl ,			/* tt - parameter to vary time integration 
				 * from explicit (tt = 1) to 
				 * implicit (tt = 0)                         */
       dbl );			/* dt - current value of the time step       */

EXTERN void mass_flux_surf_bc	/* mm_fill_species.c                         */
(double [],		/* func                                      */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
       int ,			/* wspec - species number of this BC         */
       double ,			/* mass_tran_coeff -                         */
       double ,			/* Y_c - bath concentration 	             */
       dbl ,			/* dt - current value of the time step       */
       dbl );			/* tt - parameter to vary time integration   */

EXTERN void mass_flux_BV_surf_bc	/* mm_fill_species.c                 */
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

EXTERN void mass_flux_HOR_surf_bc
(double func[],
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
       int ,          /* wspec - species number of this boundary condition   */
       double ,       /* ai0 - product of interfacial area by                */
                      /* exchange current density (A/cm^3)                   */
       double ,       /* H - thickness of catalyst layer (cm)                */
       double ,       /* cref - ref. concentration (moles/cm^3)              */
       double ,       /* alphaa - anodic direction transfer coefficient      */
       double ,       /* alphac - cathodic direction transfer coefficient    */
       double ,       /* T - cell temperature (K)                            */
       double ,       /* U0 - open-circuit potential (V)                     */
       double ,       /* beta - reaction order                               */
       double ,       /* n - number of electrons involved in rxn             */
       double ,       /* V - cell voltage (V)                                */
       double ,       /* dt - current value of the time step                 */
       double );     /* tt - parameter to vary time integration             */

EXTERN void mass_flux_ORR_surf_bc
(double func[],
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
       int ,          /* wspec - species number of this boundary condition   */
       double ,       /* ai0 - product of interfacial area by                */
                      /* exchange current density (A/cm^3)                   */
       double ,       /* H - thickness of catalyst layer (cm)                */
       double ,       /* cref - species ref. concentration (moles/cm^3)      */
       double ,       /* alphac - cathodic direction transfer coefficient    */
       double ,       /* T - cell temperature (K)                            */
       double ,       /* V - cell voltage (V)                                */
       double ,       /* U0 - open-circuit potential                         */
       double ,       /* beta - reaction order                               */
       double ,       /* n - number of electrons involved in rxn             */
       double ,       /* dt - current value of the time step                 */
       double );     /* tt - parameter to vary time integration             */

EXTERN void mass_flux_H2O_ANODE_surf_bc
(double func[],
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
       int ,          /* wspec - species number of this boundary condition   */
       double ,       /* ai0 - product of interfacial area by                */
                      /* exchange current density (A/cm^3)                   */
       double ,       /* Ha - thickness of anode catalyst layer (cm)         */
       double ,       /* cH2ref - ref. concentration for H2 (moles/cm^3)     */
       double ,       /* alphaa - anodic direction transfer coefficient      */
       double ,       /* alphac - cathodic direction transfer coefficient    */
       double ,       /* T - cell temperature (K)                            */
       double ,       /* U0 - Open-circuit potential for HOR (V)             */
       double ,       /* nd - electro-osmatic drag coefficient               */
       double ,       /* dt - current value of the time step                 */
       double );     /* tt - parameter to vary time integration             */

EXTERN void mass_flux_H2O_CATHODE_surf_bc
(double func[],
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
       int ,          /* wspec - species number of this boundary condition   */
       double ,       /* ai0 - product of interfacial area by                */
                      /* exchange current density (A/cm^3)                   */
       double ,       /* Hc - thickness of cathode catalyst layer (cm)       */
       double ,       /* cO2ref - ref. concentration for O2 (moles/cm^3)     */
       double ,       /* alphac - cathodic direction transfer coefficient    */
       double ,       /* T - cell temperature (K)                            */
       double ,       /* V - call voltage (V)                                */
       double ,       /* U0 - Open-circuit potential for ORR (V)             */
       double ,       /* nd - electro-osmatic drag coefficient               */
       double ,       /* dt - current value of the time step                 */
       double );     /* tt - parameter to vary time integration             */

EXTERN void mass_flux_SULFIDATION_surf_bc	/* mm_fill_species.c         */
(double [],		/* func                                      */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
       int ,			/* key word for sulfidation kinetic model    */
       int ,			/* wspec - species number of this BC         */
       double ,			/* nu - stoichiometric coefficient           */
       double ,			/* k1 - forward kinetic rate constant	     */
       double ,			/* E1 - forward activation energy	     */
       double ,			/* kn1 - backward kinetic rate constant	     */
       double ,			/* En1 - backward activation energy	     */
       double ,			/* T - temperature			     */
       double ,			/* c_H2S - bulk concentration of H2S	     */
       double ,			/* c_O2 - bulk concentration of O2	     */
       double ,			/* M_solid - molecular weight of Cu2S	     */
       double ,			/* rho_solid - theor. open-circuit potential */
       double ,			/* dt - current value of the time step       */
       double );      		/* tt - parameter to vary time integration   */

EXTERN void const_mass_flux_surf_bc	/* mm_fill_species.c - RSL 6/28/00 */
(double [],         /* func                                      */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func     */
       int ,              /* wspec - species number of this BC         */
       double ,           /* const_mass_flux - specified flux          */
       dbl ,              /* time  - current value of the time         */
       dbl ,              /* dt - current value of the time step       */
       dbl );            /* tt - parameter to vary time integration   */

EXTERN void mass_flux_BV2_surf_bc       /* mm_fill_species.c - RSL 5/28/02 */
(double [],         /* func                                      */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func     */
       int ,              /* wspec -- species number of this BC        */
       double ,           /* k -- rate constant                        */
       double ,           /* PHI_E -- electrode potential              */
       double ,           /* alphaa -- anodic transfer coefficient     */
       double ,           /* alphac -- cathodic transfer coefficient   */
       double ,           /* U0 -- standard open circuit potential     */
       double ,           /* T -- electrolyte temperature              */
       dbl ,              /* time - current value of the time          */
       dbl ,              /* dt - current value of the time step       */
       dbl );            /* tt - parameter to vary time integration   */

EXTERN void mass_flux_NI_surf_bc        /* mm_fill_species.c - RSL 5/28/02 */
(double [],         /* func                                      */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func     */
       int ,              /* wspec -- species number of this BC        */
       double ,           /* PHI_E -- electrode potential              */
       double ,           /* T -- electrolyte temperature              */
       dbl ,              /* time - current value of the time          */
       dbl ,              /* dt - current value of the time step       */
       dbl );            /* tt - parameter to vary time integration   */

EXTERN void current_BV2_surf_bc /* mm_fill_species.c - RSL 5/28/02 */
(double [],         /* func                                      */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func     */
       int ,              /* wspec -- species number of this BC        */
       double ,           /* k -- rate constant                        */
       double ,           /* PHI_E -- electrode potential              */
       double ,           /* alphaa -- anodic transfer coefficient     */
       double ,           /* alphac -- cathodic transfer coefficient   */
       double ,           /* U0 -- standard open circuit potential     */
       double ,           /* T -- electrolyte temperature              */
       dbl ,              /* time - current value of the time          */
       dbl ,              /* dt - current value of the time step       */
       dbl );            /* tt - parameter to vary time integration   */

EXTERN void current_NI_surf_bc  /* mm_fill_species.c - RSL 5/28/02 */
(double [],         /* func                                      */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func     */
       int ,              /* wspec -- dummy variable                   */
       double ,           /* PHI_E -- electrode potential              */
       double ,           /* T -- electrolyte temperature              */
       dbl ,              /* time - current value of the time          */
       dbl ,              /* dt - current value of the time step       */
       dbl );            /* tt - parameter to vary time integration   */

EXTERN void mass_flux_equil_mtc	/* mm_fill_species.c                         */
(dbl [MAX_CONC],		/* mass_flux                                 */
       dbl [MAX_CONC][MAX_VARIABLE_TYPES+MAX_CONC], /* d_mass_flux           */
       double [MAX_CONC],	/* activity                                  */
       double [MAX_CONC][MAX_CONC], /* dact_dC                               */
       double [MAX_CONC],	/* y_mass - conc at boundary                 */
       int ,			/* mode - model to which the VLE is based    */
       double ,			/* amb_pres - ambient pressure               */
       int ,			/* wspec - species no.                       */
       double ,			/* mass_tran_coeff - MASS transfer coeff     */
       double [MAX_VARIABLE_TYPES+MAX_CONC],	/* d_mtc     */
       double );		/* Y_c - bulk concentration 	             */

EXTERN void mtc_chilton_coburn	/* mm_fill_species.c                         */
(double *,		/* mtc                                 */
       double [MAX_VARIABLE_TYPES+MAX_CONC], /* d_mtc           */
       int ,			/* wspec - species no.                       */
       double ,			/* htc     */
       double ,			/* T_gas               */
       double ,			/* amb_pres - ambient pressure               */
       double );		/* diffusivity_gas	             */
 
EXTERN void get_equil_surf_bc	/* mm_fill_species.c                         */
(double [],		/* func                                      */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
       int ,			/* mode - model on which the VLE is based    */
       int ,			/* wspec - species number of this BC         */
       double ,			/* amb_pres                                  */
       double ,			/* mass_tran_coeff                           */
       double ,			/* Y_c - bath concentration 	             */
       double ,			/* T_gas - Chilton-Coburn 	             */
       double ,			/* diff_gas - Chilton-Coburn 	             */
       dbl ,			/* dt - current value of the time step       */
       dbl );			/* tt - parameter to vary time integration   */

EXTERN void act_coeff           /* mm_fill_species.c                         */
(dbl [MAX_CONC],          /* ln_gamma                                  */
       dbl [MAX_CONC][MAX_CONC], /* dln_gamma_dC                             */
       dbl [MAX_CONC][MAX_CONC][MAX_CONC], /* d2ln_gamma_dC2                 */
       double [MAX_CONC],        /* concentration                            */
       int ,                     /* which activity model                     */
       int );                    /* species number                           */

EXTERN void sus_mass_flux_surf_bc /* mm_fill_species.c                       */
(double [],		/* func                                      */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
       const int ,		/* sus_species - species number of this 
				 * boundary condition                        */
       const dbl ,		/* time - current value of the time          */
       const dbl ,		/* dt - current value of the time step       */
       const dbl ,		/* tt - parameter to vary time integration 
				 * from BE(0) to CN(1/2) to FE(1)            */
       const dbl []);          /* element size                              */

EXTERN void kin_bc_leak		/* mm_fill_species.c                         */
(double [],		/* func                                      */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
       double [DIM],		/* x_dot - mesh velocity                     */
       dbl ,			/* tt - parameter to vary time integration 
				 * from explicit (tt = 1) to 
				 * implicit (tt = 0)                         */
       dbl ,			/* dt - current value of the time step       */
       int ,			/* bc_input_id                               */
       struct Boundary_Condition *); /* BC_Types                            */

EXTERN void compute_leak_velocity /* mm_fill_species.c                       */
(double *,		/* vnorm                                     */
       NORMAL_VELOCITY_DEPENDENCE_STRUCT *,
       dbl ,			/* tt - parameter to vary time integration 
				 * from explicit (tt = 1) to 
				 * implicit (tt = 0)                         */
       dbl ,			/* dt - current value of the time step       */
       struct Boundary_Condition *,   /* bc                                  */
       struct Boundary_Condition *); /* fluxbc                              */

EXTERN void compute_leak_velocity_heat /* mm_fill_species.c                       */
(double *,		/* vnorm                                     */
       NORMAL_VELOCITY_DEPENDENCE_STRUCT *,
       dbl ,			/* tt - parameter to vary time integration 
				 * from explicit (tt = 1) to 
				 * implicit (tt = 0)                         */
       dbl ,			/* dt - current value of the time step       */
       struct Boundary_Condition *,   /* bc                                  */
       struct Boundary_Condition *); /* fluxbc                              */

EXTERN void compute_leak_energy /* mm_fill_species.c                       */
(double *,		/* enorm                                     */
       NORMAL_ENERGY_DEPENDENCE_STRUCT *,
       dbl ,			/* tt - parameter to vary time integration 
				 * from explicit (tt = 1) to 
				 * implicit (tt = 0)                         */
       dbl ,			/* dt - current value of the time step       */
       struct Boundary_Condition *,   /* bc                                  */
       struct Boundary_Condition *); /* fluxbc                              */

EXTERN void kin_bc_electrodeposition /* mm_fill_species.c - RSL 5/28/02   */
(double [],         /* func                                         */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func        */
       dbl ,              /* time - current value of the time             */
       dbl ,              /* tt - parameter to vary time integration      */
       dbl ,              /* dt - current value of the time step          */
       int ,              /* bc_input_id                                  */
       struct Boundary_Condition *); /* BC_Types                         */

EXTERN void vnorm_bc_electrodeposition /* mm_fill_species.c - RSL 5/30/02 */
(double [],         /* func                                         */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func        */
       dbl ,              /* tt - parameter to vary time integration      */
       dbl ,              /* time - current time                          */
       dbl ,              /* dt - current value of the time step          */
       int ,              /* bc_input_id                                  */
       struct Boundary_Condition *); /* BC_Types                         */

EXTERN void lat_heat_bc		/* mm_fill_species.c                         */
(double [],		/* func                                      */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
       double [DIM],		/* x_dot - mesh velocity                     */
       dbl ,			/* tt - parameter to vary time integration 
				 * from explicit (tt = 1) to 
				 * implicit (tt = 0)                         */
       dbl ,			/* dt - current value of the time step       */
       int ,			/* bc_input_id                               */
       struct Boundary_Condition *); /* BC_Types                            */

EXTERN void lat_heat_internal_bc /* mm_fill_species.c                        */
(double [],		/* func                                      */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
       const double ,		/* tt - parameter to vary time integration 
				 * from explicit (tt = 1) to 
				 * implicit (tt = 0)                         */
       const double ,		/* dt - current value of the time step       */
       const int ,		/* bc_input_id                               */
       const struct Boundary_Condition *); /* BC_Types                      */

EXTERN int get_convection_velocity /* cf mm_fill_species.c */
(double [DIM],		/* vconv - Calculated convection velocity */
       double [DIM],		/* vconv_old - Calculated convection velocity
				 * at previous time*/
       CONVECTION_VELOCITY_DEPENDENCE_STRUCT *,
       double ,			/* dt */
       double );		/* tt */

EXTERN int get_continuous_species_terms	/* mm_fill_species.c                 */
(struct Species_Conservation_Terms *, /* st                            */
       double ,			/* time - present time value; KSC: 10/98     */
       double ,                 /* tt - parameter to vary time integration
				 * from explicit (tt = 1) to
				 * implicit (tt = 0)                         */
       double , 		/* dt - current time step size               */
       const double [DIM]);    /* element size                              */

EXTERN void ludcmp		/* mm_fill_species.c                         */
(double [MAX_CONC*DIM][MAX_CONC*DIM], /* a                             */
       int ,			/* n                                         */
       int [MAX_CONC*DIM],	/* indx                                      */
       double );		/* d                                         */

EXTERN void lubksb		/* mm_fill_species.c                         */
(double [MAX_CONC*DIM][MAX_CONC*DIM], /* a                             */
       int ,			/* n                                         */
       int [MAX_CONC*DIM],	/* indx                                      */
       double [MAX_CONC*DIM]);	/* b                                         */

EXTERN int Stefan_Maxwell_diff_flux /* mm_fill_species.c                     */
(struct Species_Conservation_Terms *, /* st                            */
       double ,			/* time                                      */
       double );		/* dt                                        */

EXTERN int fickian_charged_flux /* mm_fill_species.c                         */
(struct Species_Conservation_Terms *, /* st                            */
       int   );                /* w                                         */

EXTERN int fickian_charged_flux_x /* mm_fill_species.c -- RSL 9/18/00        */
(struct Species_Conservation_Terms *, /* st                            */
       double ,                 /* time                                      */
       double );               /* dt                                        */

EXTERN int fickian_flux		/* mm_fill_species.c                         */
(struct Species_Conservation_Terms *, /* st                            */
       int );			/* w                                         */

EXTERN int generalized_fickian_flux /* mm_fill_species.c                     */
(struct Species_Conservation_Terms *, /* st                            */
       int );			/* w                                         */

EXTERN int assemble_invariant	/* mm_fill_species.c                         */
(double ,			/* tt - parm to vary time integration from 
				 * explicit (tt = 1) to implicit (tt = 0)    */
       double );		/* dt - time step size                       */

EXTERN int get_particle_convection_velocity /* mm_fill_species.c             */
(double [DIM],		/* pvconv                                    */
       double [DIM],		/* pvconv_old - previous time                */
       CONVECTION_VELOCITY_DEPENDENCE_STRUCT *,
       double ,			/* dt                                        */
       double );		/* tt   */
	   
EXTERN int ls_modulate_speciessource
( int ,
	    struct Species_Conservation_Terms * );
		
EXTERN void fickian_charged_gradient_bc
( double [],
	double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
	int ,	/* species number of this boundary condition */
	double );/* Value for normal surface gradient component       */
		
#endif /* GOMA_MM_FILL_SPECIES_H */
