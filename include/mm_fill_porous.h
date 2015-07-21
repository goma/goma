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
 * mm_fill_porous.h -- prototype declarations for mm_fill_porous.c
 */

#ifndef _MM_FILL_POROUS_H
#define _MM_FILL_POROUS_H

#ifdef EXTERN
#undef EXTERN
#endif

#ifdef _MM_FILL_POROUS_C
#define EXTERN /* do nothing */
#endif

#ifndef _MM_FILL_POROUS_C
#define EXTERN extern
#endif

EXTERN int assemble_porous_transport /* mm_fill_porous.c                     */
PROTO((double ,                 /* time - present time valuel; KSC           */
       double ,                 /* tt - parm to vary time integration from
                                 * explicit (tt = 1) to implicit (tt = 0)    */
       double));                /* dt - current time step size               */
 

EXTERN int load_porous_properties
PROTO((void));

extern double
get_supg_terms_porous(double [DIM], double [DIM][DIM]);

extern void load_nodal_porous_properties(double, double);
extern void load_nodal_shell_porous_properties(double, double, int);

EXTERN int get_porous_part_sat_terms
PROTO((struct Porous_Media_Terms *, /* pm                            */
       double ,			/* tt - time integration scheme param        */
       double ));		/* dt - current time step size               */

EXTERN int get_porous_part_sat_terms_decoupled
PROTO((struct Porous_Media_Terms *, /* pm                            */
       double ,			/* tt - time integration scheme param        */
       double ));		/* dt - current time step size               */

EXTERN int get_porous_fully_sat_terms
PROTO((struct Porous_Media_Terms *, /* pm                            */
       double ,			/* tt - parm to vary time integration from 
				 * explicit (tt = 1) to implicit (tt = 0)    */
       double ));		/* dt - current time step size               */

EXTERN void porous_mass_flux_surf_bc
PROTO((double [DIM],		/* func                                      */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
       int ,			/* wspec - species number of this BC         */
       double ,			/* mass_tran_coeff - units consistent with 
				 * gas phase concentration driving force     */
       double ,			/* Y_c - bath concentration 	             */
       double ,			/* mass_tran_coeff1 - units consistent with 
				 * gas phase concentration driving force
				 * pressure driving force from liquid phase*/
       double ,			/* sink pressure 	             */
       dbl ,			/* dt - current value of the time step       */
       dbl ));                  /* tt - parm to vary time integration        */

EXTERN void porous_convection_bc
PROTO((double [DIM],		/* func                                      */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
       int ,			/* wspec -species number of this BC          */
       dbl ,			/* dt - current value of the time step       */
       dbl ));                  /* tt - parm to vary time integration        */

EXTERN void porous_kinematic_bc
PROTO((double [DIM],		/* func                                      */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
       int ,			/* wspec - species number of this BC         */
       dbl ,			/* dt - current value of the time step       */
       dbl ,			/* tt - parm to vary time integration        */
       dbl ));			/* vflux                                     */

EXTERN void porous_normal_velocity_bc 
PROTO((double [],		/* func                                      */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
       double [DIM],		/* x_dot - mesh velocity vector              */
       dbl ,			/* tt - parm to vary time integration from 
				 * explicit (tt = 1) to implicit (tt = 0)    */
       dbl ,			/* dt - current value of the time step       */
       int ,			/* bc_input_id                               */
       struct Boundary_Condition *, /* BC_Types                              */
       int ,			/* i_mat_solid - mat block id porous phase   */
       int ,			/* i_mat_fluid - mat block id gas phase      */
       int ,			/* wspec                                     */
       double ));		/* dens_vap - density of pure solvent vapor  */

EXTERN void put_gas_flux_in_pores
PROTO((double [],		/* func                                      */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
       double [DIM],		/* x_dot - mesh velocity vector              */
       dbl ,			/* tt - parm to vary time integration from 
				 * explicit (tt = 1) to implicit (tt = 0)    */
       dbl ,			/* dt - current value of the time step       */
       int ,			/* bc_input_id                               */
       struct Boundary_Condition *, /* BC_Types                              */
       int ,			/* i_mat_solid - mat block id porous phase   */
       int ,			/* i_mat_fluid - mat block id gas phase      */
       int ,			/* wspec                                     */
       double ,			/* dens_vap - density of pure solvent vapor  */
       double ));		/* vapor_recoil                              */

EXTERN void porous_vapor_equil_bc
PROTO((double [],		/* func                                      */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
       double [DIM],		/* x_dot - mesh velocity vector              */
       dbl ,			/* tt - parm to vary time integration from 
				 * explicit (tt = 1) to implicit (tt = 0)    */
       dbl ,			/* dt - current value of the time step       */
       int ,			/* bc_input_id                               */
       struct Boundary_Condition *, /* BC_Types                              */
       int ,			/* i_mat_solid - mat block id porous phase   */
       int ,			/* i_mat_fluid - mat block id gas phase      */
       int ,			/* wspec                                     */
       double ));		/* amb_pres - ambient pressure, perhaps?     */

EXTERN double load_permeability
PROTO((void));	      

EXTERN void load_permeability_tensor
PROTO((void));      

EXTERN void load_shell_permeability_tensor
PROTO((void));      

EXTERN double load_saturation
PROTO((double ,			/* porosity                                  */
       double ,			/* cap_pres                                  */ 
       double [2]));            /* d_cap_pres                         */

EXTERN void load_enthalpy
PROTO((double ,                /* saturation                                */
       double));               /* pressure                                  */

EXTERN double calc_rho_gas
PROTO((  double , 
         double ,
         double *,
         double *,
         double *,
         double *,
         double *,
         double *,
         double *,
         double *,
         double *));

EXTERN double calc_Ywg
PROTO((  double ,
         double ,
         double , 
         double ,
         double ,
         double ,
         double *,
         double *,
         double *,
         double *,
         double *,
         double *));

EXTERN double rho_sat_water_vap_EOS
PROTO(( double ,
        double ,
        double *,		/* drho_wg_dP                                    */
        double *,		/* drho_wg_dT                                    */
        double *,               /* d_drhowg_dP                                  */
        double *));		/* d_drhowg_dT                                  */

EXTERN double rho_air_EOS
PROTO((   double , 
          double ,
          double ,
          double ,
          double *,
          double *,
          double *,
          double *));

EXTERN double eval_rho_poly 
PROTO(( double *,
        double ,
        double ));

EXTERN double P_water_sat_EOS
PROTO((double ,			/* temperature                               */
       double *,		/* dPsat_dT                                  */
       double *));		/* d_dPsat_dT                                */

EXTERN double h_air_EOS
PROTO((double ,			/* temperature                               */
       double *,		/* dha_dT                                    */
       double *));		/* d_dha_dT                                  */

EXTERN double h_water_compressed_EOS
PROTO((double ,			/* pressure                                  */
       double ,                 /* temperature                               */
       double *,		/* dhl_dP                                    */
       double *,		/* dhl_dT                                    */
       double *,                /* d_dhl_dP                                  */
       double *));		/* d_dhl_dT                                  */

EXTERN double h_water_superheat_EOS
PROTO((double ,			/* pressure                                  */
       double ,                 /* temperature                               */
       double *,		/* dhwg_dP                                    */
       double *,		/* dhwg_dT                                    */
       double *,                /* d_dhwg_dP                                  */
       double *));		/* d_dhwg_dT                                  */

EXTERN double eval_poly
PROTO((double *,		/* coef[]                                    */
       double ,		        /* pressure                                  */
       double ));		/* temperature                               */

EXTERN double eval_poly_dP
PROTO((double *,		/* coef[]                                    */
       double ,		        /* pressure                                  */
       double ,		        /* temperature                               */
       double *));		/* d_dP_poly                                 */

EXTERN double eval_poly_dT
PROTO((double *,		/* coef[]                                    */
       double ,		        /* pressure                                  */
       double ,		        /* temperature                               */
       double *));		/* d_dT_poly                                 */

EXTERN double eval_pv_poly
PROTO((int ,			/* order                                     */
       double *,		/* coef[]                                    */
       double ));		/* temperature                               */

EXTERN double eval_pv_poly_dT
PROTO((int ,			/* order                                     */
       double *,		/* coef[]                                    */
       double ,  		/* temperature                               */
       double * ));             /* df_d_dT                                   */

EXTERN void load_gas_conc
PROTO((double ,			/* porosity                                  */
       double ,			/* cap_pres                                  */
       double [2]));		/* d_cap_pres                                */

EXTERN void load_gas_conc_flat
PROTO((double ,			/* porosity                                  */
       double ,			/* cap_pres                                  */
       double [2]));		/* d_cap_pres                                */

EXTERN void load_gas_conc_EOS
PROTO((double ,			/* porosity                                  */
       double ,			/* cap_pres                                  */
       double [2]));		/* d_cap_pres                                */

EXTERN void load_bulk_density
PROTO((double ,			/* porosity                                  */
       double ,			/* cap_pres                                  */
       double ,			/* saturation                                */
       double [2]));		/* d_cap_pres                                */

EXTERN void load_liq_perm
PROTO((double ,			/* porosity                                  */
       double ,			/* cap_pres                                  */
       double ,			/* saturation                                */
       double [2]));		/* d_cap_pres                                */

EXTERN void load_gas_perm
PROTO((double ,			/* porosity                                  */
       double ,			/* cap_pres                                  */
       double ,			/* saturation                                */
       double [2]));		/* d_cap_pres                                */

EXTERN void load_gas_diff
PROTO((double ,			/* porosity                                  */
       double ,			/* cap_pres                                  */
       double ,			/* saturation                                */
       double [2],		/* d_cap_pres                                */
       int ));			/* species number                            */

EXTERN void load_mass_flux
PROTO((double ,			/* porosity                                  */
       double ,			/* cap_pres                                  */
       double ,			/* saturation                                */
       double [2]));		/* d_cap_pres                                */

EXTERN void load_MandE_flux
PROTO((double ,			/* porosity                                  */
       double ,			/* cap_pres                                  */
       double ,			/* saturation                                */
       double [2]));		/* d_cap_pres                                */

EXTERN void calc_darcy_velocity
PROTO((void ));

EXTERN void load_gas_conc
PROTO((double ,			/* porosity                                  */
       double ,			/* cap_pres                                  */
       double [2]));		/* d_cap_pres                                */

EXTERN void load_gas_conc_flat
PROTO((double ,			/* porosity                                  */
       double ,			/* cap_pres                                  */
       double [2]));		/* d_cap_pres                                */

EXTERN void load_bulk_density
PROTO((double ,			/* porosity                                  */
       double ,			/* cap_pres                                  */
       double ,			/* saturation                                */
       double [2]));		/* d_cap_pres                                */

EXTERN void load_liq_perm
PROTO((double ,			/* porosity                                  */
       double ,			/* cap_pres                                  */
       double ,			/* saturation                                */
       double [2]));		/* d_cap_pres                                */

EXTERN void load_gas_perm
PROTO((double ,			/* porosity                                  */
       double ,			/* cap_pres                                  */
       double ,			/* saturation                                */
       double [2]));		/* d_cap_pres                                */

EXTERN void load_gas_diff
PROTO((double ,			/* porosity                                  */
       double ,			/* cap_pres                                  */
       double ,			/* saturation                                */
       double [2],		/* d_cap_pres                                */
       int ));			/* species number                            */

EXTERN void load_mass_flux
PROTO((double ,			/* porosity                                  */
       double ,			/* cap_pres                                  */
       double ,			/* saturation                                */
       double [2]));		/* d_cap_pres                                */

EXTERN void porous_pressure
PROTO((double *,		/* func                                      */
       double [],		/* d_func                                    */
       int ,			/* i_mat_solid                               */
       int ));			/* i_mat_fluid                               */

EXTERN void porous_pressure_lub
PROTO((double *,		/* func                                      */
       double [],		/* d_func                                    */
       const int ,	       	/* id_side                                   */
       double [DIM],            /* isoparametric coordinates                 */
       const Exo_DB *,          /* Exo struct                                */
       const double));		/* scale factor                             */

EXTERN void sat_darcy_continuous_bc
PROTO((double [],		/* func                                      */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
       const double ,		/* tt - parm vary time integration method.   */
       const double ,		/* dt - current time step size               */
       const double ,		/* time  - current time                      */
       const int ,		/* i_mat_solid - mat block id porous phase   */
       const int ,		/* i_mat_fluid - mat block id gas phase      */
       const double,            /* length scale for level-set                */
       const double));          /* v_attached for level-set case             */


EXTERN int evaluate_sat_hyst_criterion
PROTO((int,		        /* ip                                     */
       int,                     /* ielem                                  */
       struct Porous_Media_Terms *, /* pm                                 */
       const double ,		/* tt - parm vary time integration method.   */
       const double ));		/* dt - current time step size               */

EXTERN void porous_liq_fill
PROTO((double *,	       /* func                                      */
       double [],              /* d_func           */
       const int, 	       /* i_mat_solid - mat block id porous phase   */
       const int,   	       /* i_mat_fluid - mat block id gas phase      */
       double,    	       /* cap press fluid                           */ 
       double,    	       /* cap press fluid                           */ 
       double,    	       /* level set length scale                    */ 
       MATRL_PROP_STRUCT *));  /* material property struct                  */ 
	   
void por_liq_flux_fill 
PROTO(( double *,
	double [MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
	double , 
	double ,
	const int, 
	const int ,
	double ,
	double ,
	double ,
	double ,
	int      ));

double por_mass_source_model
PROTO(( double d_MassSource[MAX_VARIABLE_TYPES + MAX_CONC][MDE] ));

#endif /* _MM_FILL_POROUS_H */
