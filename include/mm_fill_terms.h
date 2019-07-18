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
 

#ifndef _MM_FILL_TERMS_H
#define _MM_FILL_TERMS_H

#ifdef EXTERN
#undef EXTERN
#endif

#ifdef _MM_FILL_TERMS_C
#define EXTERN /* do nothing */
#endif

#ifndef _MM_FILL_TERMS_C
#define EXTERN extern
#endif


EXTERN int assemble_mesh	/* mm_fill_terms.c                           */
PROTO((double ,			/* time                                      */
       double ,			/* tt                                        */
       double ,			/* dt                                        */
       int ,			/* ielem - current element number            */
       int ,			/* ip - current integration point            */
       int ));			/* ip_total - total gauss integration points */

EXTERN int assemble_energy	/* mm_fill_terms.c                           */
PROTO((	double ,		/* time - present time value                 */
	double ,		/* tt - parameter to vary time integration 
					* from explicit (tt = 1) to 
					* implicit (tt = 0)                   */
	double ,			/* dt - current time step size        */
	const PG_DATA *));	/* dvc_dnode                                 */

EXTERN int assemble_momentum	/* mm_fill_terms.c                           */
PROTO((	double ,		/* time - present time value                 */
	dbl ,			/* tt - parm to vary time integration from 
				 * explicit (tt = 1) to implicit (tt = 0)    */
	dbl ,			/* dt - current time step size               */
	dbl ,                    /* h_elem_avg - average global element size  */
	const PG_DATA * , 	/* dvc_dnode                                 */
	dbl xi[DIM],            /* Local stu coords */
	const Exo_DB *exo));      /* ExodusII database struct pointer */

EXTERN int assemble_continuity	/* mm_fill_terms.c                           */
PROTO((dbl,                     /* time_value */
       dbl ,			/* tt - to vary time integration from 
				   explicit (tt = 1) to implicit (tt = 0)    */
       dbl ,			/* dt - current time step size               */
       const PG_DATA *));


EXTERN int calc_pspg		/* mm_fill_terms.c                           */
PROTO((    dbl [DIM],
	   PSPG_DEPENDENCE_STRUCT *,
	   dbl ,                /* current time                              */
	   dbl ,		/* parameter to vary time integration from
				   explicit (tt = 1) to implicit (tt = 0)    */
	   dbl ,		/* current time step size                    */
	   const PG_DATA * ));

EXTERN int calc_cont_gls		/* mm_fill_terms.c                           */
PROTO((    dbl *,
	   CONT_GLS_DEPENDENCE_STRUCT *,
	   dbl ,                /* current time                              */
	   const PG_DATA * ));

EXTERN int assemble_momentum_path_dependence
PROTO(( double ,                /* time */
	double ,                /* tt, parameter to vary time integration from
                                 * explicit (tt = 1) to implicit (tt = 0)    */
        double ,                /* dt, current time step size                */
	const PG_DATA * ));     /* PG data needed for continuity stabilization */

EXTERN int assemble_continuity_path_dependence
PROTO((dbl,                     /* time_value */
       dbl ,			/* tt - to vary time integration from 
				   explicit (tt = 1) to implicit (tt = 0)    */
       dbl ,			/* dt - current time step size               */
       const PG_DATA *));	/* d(v_avg)dvj for PSPG calcs */

EXTERN int assemble_extension_velocity_path_dependence
PROTO(( ));			/* no args */

#if 0 
EXTERN int assemble_continuity_path_dependence
PROTO((double ,                 /* time_value                                */
       double ,                 /* tt, parameter to vary time integration from *
                                 * explicit (tt = 1) to implicit (tt = 0)    */
       double ,                 /* dt, current time step size                */
       double ,                 /* h_elem_avg,  average global element size  *
                                 * for PSPG, taken constant wrt to Jacobian  */
       double [],               /* hsquared[DIM], element size info for PSPG */
       double [][DIM],		/* hh[DIM][DIM], not currently used, but     *
                                 * left in just in case they're needed later */
       double [][MDE],		/* dh_dxnode[DIM][MDE]                       */
       double ,		        /* U_norm, global velocity norm for PSPG calcs */
       double ,		        /* mu_avg, element viscosity for PSPG calcs  */
       double ,		        /* rho_avg,element density for PSPG calcs    */
       double [],		/* v_avg[DIM], element velocity for PSPG calcs */
       double [][MDE]));        /* dv_dnode[DIM][MDE],deriv.velocity wrt nodal variables */

#endif

EXTERN int assemble_volume	/* mm_fill_terms.c                           */
PROTO((bool));

EXTERN int assemble_curvature
PROTO((void )); /* mm_fill_terms.c                         */


EXTERN int assemble_normals     /* mm_fill_terms.c                         */
PROTO ((void ));

EXTERN int assemble_ls_momentum_source
PROTO(( void ));

EXTERN int apply_ls_momentum_source
PROTO(( void ));

EXTERN int assemble_div_normals   /* mm_fill_terms.c                         */
PROTO ((void ));

EXTERN int assemble_LSvelocity	/* mm_fill_terms.c                           */
PROTO((bool, int));

EXTERN int assemble_acoustic	/* mm_fill_terms.c                           */
PROTO((	double ,		/* time - present time value         */
	double ,		/* tt - parameter to vary time integration
			        	* from explicit (tt = 1) to 
					* implicit (tt = 0)                   */
	double ,		/* dt - current time step size               */
	const PG_DATA *,	/* dvc_dnode                                 */
	const int ,		/*  acoustic eqn id and var id		     */
	const int ));	

EXTERN int assemble_poynting	/* mm_fill_terms.c                           */
PROTO((	double ,		/* time - present time value         */
	double ,		/* tt - parameter to vary time integration
			        	* from explicit (tt = 1) to 
					* implicit (tt = 0)                   */
	double ,		/* dt - current time step size               */
	const PG_DATA *,	/* dvc_dnode                                 */
	const int ,		/*  Light intensity eqn id and var id		     */
	const int ));	

EXTERN int assemble_emwave	/* mm_fill_terms.c                           */
PROTO((	double ,		/* time - present time value         */
	double ,		/* tt - parameter to vary time integration
			        	* from explicit (tt = 1) to 
					* implicit (tt = 0)                   */
	double ,		/* dt - current time step size               */
	const PG_DATA *,	/* dvc_dnode                                 */
	const int ,		/*  Light intensity eqn id and var id		     */
	const int ,		/*  Light intensity eqn id and var id		     */
	const int ));	

EXTERN int assemble_acoustic_reynolds_stress	/* mm_fill_terms.c */
PROTO(( double,					/* time */
        double,					/* tt */
        double,					/* dt */
        const PG_DATA *));			/* dvc_dnode */

EXTERN int assemble_pore_sink_mass	/* mm_fill_terms.c                           */
PROTO((	double ,		/* time - present time value         */
	double ,		/* tt - parameter to vary time integration
			        	* from explicit (tt = 1) to 
					* implicit (tt = 0)                   */
	double ));		/* dt - current time step size               */


EXTERN int assemble_acoustic_energy	/* mm_fill_terms.c                 */
PROTO((	double ,		/* time - present time value         */
	double ,		/* tt - parameter to vary time integration
			        	* from explicit (tt = 1) to 
					* implicit (tt = 0)                 */
	double ,		/* dt - current time step size               */
	const PG_DATA *));	/* dvc_dnode                                 */

EXTERN int load_fv		/* mm_fill_terms.c                           */
PROTO((void ));

EXTERN int load_fv_grads	/* mm_fill_terms.c                           */
PROTO((void ));

EXTERN int load_fv_mesh_derivs	/* mm_fill_terms.c                           */
PROTO((int ));                  /* okToZero - Turns on zeroing in the function
                                   This is usually turned on except when accumulating */

EXTERN double density		     /* mm_fill_terms.c                           */
PROTO((DENSITY_DEPENDENCE_STRUCT *,  /* density dependence */
       double   ));                  /* time */

EXTERN double conductivity		/* mm_fill_terms.c             */
PROTO((CONDUCTIVITY_DEPENDENCE_STRUCT *,
       dbl      ));             /* time */

EXTERN double heat_capacity		/* mm_fill_terms.c                  */
PROTO((HEAT_CAPACITY_DEPENDENCE_STRUCT *,
       dbl      ));             /* time */

EXTERN double heat_source		/* mm_fill_terms.c                  */
PROTO((HEAT_SOURCE_DEPENDENCE_STRUCT *,
       double ,			/* time - present time value                 */
       double ,			/* tt - parameter to vary time integration
				 * from explicit (tt = 1) to
				 * implicit (tt = 0)                         */
       double ));		/* dt - current time step size               */

EXTERN double acoustic_impedance		/* mm_fill_terms.c             */
PROTO((CONDUCTIVITY_DEPENDENCE_STRUCT *,
       dbl      ));             /* time */

EXTERN double wave_number		/* mm_fill_terms.c             */
PROTO((CONDUCTIVITY_DEPENDENCE_STRUCT *,
       dbl      ));             /* time */

EXTERN double acoustic_absorption		/* mm_fill_terms.c             */
PROTO((CONDUCTIVITY_DEPENDENCE_STRUCT *,
       dbl      ));             /* time */

EXTERN double refractive_index		/* mm_fill_terms.c             */
PROTO((CONDUCTIVITY_DEPENDENCE_STRUCT *,
       dbl      ));             /* time */

EXTERN double light_absorption		/* mm_fill_terms.c             */
PROTO((CONDUCTIVITY_DEPENDENCE_STRUCT *,
       dbl      ));             /* time */

EXTERN double extinction_index		/* mm_fill_terms.c             */
PROTO((CONDUCTIVITY_DEPENDENCE_STRUCT *,
       dbl      ));             /* time */

EXTERN int momentum_source_term	/* mm_fill_terms.c                           */
PROTO((dbl [DIM],		/* f - Body force.                           */
       MOMENTUM_SOURCE_DEPENDENCE_STRUCT *,
       double ));

EXTERN int ls_modulate_momentumsource 
PROTO (( double [DIM],
	 double [DIM],
	 double ,
	 double ,
	 double ,
	 MOMENTUM_SOURCE_DEPENDENCE_STRUCT * ));


EXTERN void apply_table_mp
PROTO(( double *func, 
	struct Data_Table *table));

EXTERN dbl solidification_permeability
PROTO(( dbl, 
	dbl [MAX_CONC][MDE]));    /* continuous surface tension */
        
EXTERN int continuous_surface_tension
PROTO(( double, 
		double [DIM][DIM],        /* continuous surface tension                 */
		double [DIM][DIM][MDE], /* derivative w.r.t. FILL                     */
		double [DIM][DIM][DIM][MDE] )); /* d with respect to mesh */

EXTERN double quad_isomap_invert
PROTO(( const double,	/*  coordinate1  */
        const double,	/*  coordinate2  */
        const double,	/*  coordinate3  */
        const double [],	/* grid points of coordinate1 */
        const double [],	/* grid points of coordinate2 */
        const double [],	/* grid points of coordinate3 */
        const double [],	/* function values at grid points */
	const int,		/* # of grid points in direction 1 */
        const int,		/* # of grid points in direction 2 */
        const int,		/* # of grid points in direction 3 */
	const int,		/* element order(2=biquadratic, 1=bilinear) */
        const int,              /* element dimension */
        double [] ));           /* gradient array */

extern void load_matrl_statevector
PROTO(( MATRL_PROP_STRUCT *));

EXTERN double FoamVolumeSource
PROTO (( double ,
	 double ,
	 double ,
	 double [DIM][MDE],
	 double [MDE],
	 double [DIM][MDE],
	 double [MAX_CONC][MDE],
	 double [MDE] ));

EXTERN double REFVolumeSource
PROTO (( double ,
	 double ,
	 double  ,
	 double [DIM][MDE],
	 double [MDE],
	 double [DIM][MDE],
	 double [MAX_CONC][MDE] ));

EXTERN int assemble_projection_stabilization 
PROTO (( Exo_DB * ));

EXTERN int assemble_PPPS_generalized
PROTO (( Exo_DB * ));
         
EXTERN int apply_distributed_sources 
PROTO (( int,
	 double,
         double [],
         Exo_DB *,
         double,
         double,
         double,
	 const PG_DATA *,
         int,
         double *,
         double **,
         double ** ));

EXTERN int assemble_csf_tensor 
     PROTO (( void ));
     
EXTERN int assemble_curvature_with_normals_source
PROTO (( void  ));

EXTERN int assemble_curvature_source
     PROTO (( void ));

EXTERN int assemble_div_n_source
     PROTO (( void  ));

EXTERN int assemble_div_s_n_source
     PROTO (( void ));

EXTERN int
assemble_cap_hysing(double dt, double scale);

EXTERN int
assemble_cap_denner_diffusion(double dt, double scale);

EXTERN int
assemble_cap_denner_diffusion_n(double dt, double scale);

EXTERN void grad_vector_fv_fill
PROTO (( double ***,
		 double ( *)[DIM][DIM][DIM],
		 int,
		 double ( * ) [DIM] ));

EXTERN void grad_scalar_fv_fill 
PROTO (( double **,
		 double (*)[DIM], 
		 int ,
		 double * ));		

EXTERN void scalar_fv_fill
PROTO(( double **, 
	double **, 
	double **,
	double *, 
	int ,
	double *, 
	double *, 
	double * ));

EXTERN double scalar_fv_fill_adjmatrl
PROTO ((double **,
	int , 
	int ,
	int ));

EXTERN int assemble_q_source
PROTO (( double ));

EXTERN int assemble_qlaser_source
PROTO (( const double [],
         double ));

EXTERN int assemble_qvapor_source
PROTO (( const double [] ));

EXTERN int assemble_qrad_source
PROTO ((double,
        double,
        double,
        double ));

EXTERN int assemble_t_source
PROTO (( double,
         double ));

EXTERN int assemble_cont_t_source
PROTO (( double * ));

EXTERN int assemble_ls_yflux_source
PROTO(( int,
        double,
        double,
        double,
        double,
        double,
	int,
        struct Boundary_Condition * ));

EXTERN int assemble_ls_latent_heat_source
PROTO(( double,
        double,
        double,
        double,
        double,
	int,
        struct Boundary_Condition * ));

EXTERN int assemble_ars_source
PROTO (( double,
	 double ));

EXTERN int assemble_cont_vel_source
PROTO (( double *,
         Exo_DB * ));

EXTERN int assemble_extv_kinematic
PROTO(( dbl,
        dbl,
        dbl,
        int,
        struct Boundary_Condition * ));

EXTERN int assemble_p_source
PROTO (( double, const int ));

EXTERN int assemble_precoil_source
PROTO (( const double [] ));

EXTERN int assemble_uvw_source
PROTO (( int,
         double ));

EXTERN int assemble_extension_velocity_path_dependence
PROTO((void));


EXTERN int assemble_LM_source
PROTO (( double *,
	 int ,
	 double *,
	 double **,
	 double **,
	 double [],
	 Exo_DB * ));

EXTERN int assemble_interface_extension_velocity_sic
PROTO (( int ));

EXTERN int assemble_eik_kinematic
PROTO((dbl,
       dbl,
       dbl,
       int,
       struct Boundary_Condition * ));

EXTERN int assemble_fill_path_dependence
PROTO (( void ));

EXTERN int assemble_energy_path_dependence
PROTO((double ,			/* time - present time value                 */
       double ,			/* tt - parameter to vary time integration 
				 * from explicit (tt = 1) to 
				 * implicit (tt = 0)                         */
       double ,
       const PG_DATA *));	/* dvc_dnode                                 */

EXTERN int assemble_continuity_path_dependence 
PROTO((dbl ,
	   dbl ,	/* parameter to vary time integration from explicit (tt = 1) to implicit (tt = 0)    */
	   dbl ,	/* current time step size                    */
       const PG_DATA *)); /* deriv. velocity wrt nodal variables   */
			  

EXTERN double ls_modulate_thermalconductivity
PROTO (( double ,
	 double ,
	 double ,
	 double ,
	 double ,
	 CONDUCTIVITY_DEPENDENCE_STRUCT * ));

EXTERN double ls_modulate_heatcapacity
PROTO (( double ,
	 double ,
	 double ,
	 double ,
	 double ,
	 HEAT_CAPACITY_DEPENDENCE_STRUCT * ));
	 
	 
EXTERN int ls_modulate_heatsource 
PROTO(( double *,
		double ,
		double ,
		double ,
		double ,
		HEAT_SOURCE_DEPENDENCE_STRUCT * ));
	   

EXTERN void fluid_stress
PROTO (( double [DIM][DIM],                      /* Pi[DIM][DIM] */
         STRESS_DEPENDENCE_STRUCT * ));          /* d_Pi         */

EXTERN void fluid_stress_conf
PROTO (( double [DIM][DIM],                      /* Pi[DIM][DIM] */
         STRESS_DEPENDENCE_STRUCT * ));          /* d_Pi         */

EXTERN void heat_flux
PROTO (( double [DIM],                      /* q[DIM] */
         HEAT_FLUX_DEPENDENCE_STRUCT *,          /* dq     */
         double ));                         /* time   */
		 
EXTERN void acoustic_flux
PROTO (( double [DIM],                      /* q[DIM] */
         ACOUSTIC_FLUX_DEPENDENCE_STRUCT *,          /* dq     */
         double,              			/* time   */
	 const int,		/* acoustic eqn id and var id	*/
	 const int ));           
		 
EXTERN int assemble_pf_capillary
PROTO (( double * ));

EXTERN int assemble_qvapor_source 
PROTO (( const double [])) ;

EXTERN int assemble_qrad_source 
PROTO ((double,
		double, 
		double,
		double));

EXTERN double visc_diss_acoustic_source 
PROTO((HEAT_SOURCE_DEPENDENCE_STRUCT *,
       dbl *,			/* param - General multipliers   */
       int));			/* number of parameters   */

EXTERN double em_diss_heat_source 
PROTO((HEAT_SOURCE_DEPENDENCE_STRUCT *,
       dbl *,			/* param - General multipliers   */
       int));			/* number of parameters   */

EXTERN int assemble_max_strain
PROTO(( void ));

EXTERN int assemble_cur_strain
PROTO(( void ));

#endif /* _MM_FILL_TERMS_H */
