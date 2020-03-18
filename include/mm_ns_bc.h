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
 

#ifndef GOMA_MM_NS_BC_H
#define GOMA_MM_NS_BC_H

#include "el_elm.h"
#include "exo_struct.h"
#include "mm_as_structs.h"
#include "mm_fill_jac.h"
#include "mm_mp_const.h"
#include "mm_numjac.h"
#include "rf_bc_const.h"
#include "rf_fem_const.h"
#include "std.h"
#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_MM_NS_BC_C
#define EXTERN /* do nothing */
#endif

#ifndef GOMA_MM_NS_BC_C
#define EXTERN extern
#endif


EXTERN void ls_attach_bc
(double [DIM],
       double [DIM][MAX_VARIABLE_TYPES+MAX_CONC][MDE],
       const double );

EXTERN void fvelo_normal_bc
(double [DIM],		/* func                                      */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
       const double ,		/* vnormal - normal velocity                 */
       const int,               /* contact flag */
       const double [MAX_PDIM], /* x_dot - Bad name, says Phil! 
				 * -mesh velocity vector                     */
       const double ,		/* tt - parameter to vary time integration 
				 * from explicit (tt = 1) to 
				 * implicit (tt = 0)                         */
       const double ,		/* dt - current value of the time step       */
       const int ,               /* bc id */
       const double ,		/* interface zone half-width           */
       const double ,		/* interface zone shift                */
       const double );		/* gas leak angle (degrees)            */

EXTERN void fmesh_etch_bc
(double *,            /* func                                      */
       double [MAX_VARIABLE_TYPES + MAX_CONC], /* d_func           */
       const int,               /* Etch plane                 */
       const int,               /* Local node ID                 */
       const double [MAX_PDIM], /* Mesh velocity */
       const double ,           /* tt - parameter to vary time integration 
                                 * from explicit (tt = 1) to 
                                 * implicit (tt = 0)                         */
       const double );         /* dt - current value of the time step       */


EXTERN void fvelo_tangential_ls_bc
(double [DIM],		/* func                                      */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
       const double ,		/* vtangent                 */
       const double [MAX_PDIM], /* x_dot - Bad name, says Phil! 
				 * -mesh velocity vector                     */
       const double ,		/* tt - parameter to vary time integration 
				 * from explicit (tt = 1) to 
				 * implicit (tt = 0)                         */
       const double ,		/* dt - current value of the time step       */
       const double ,		/* interface zone half-width           */
       const double ,		/* interface zone shift                */
       const double );		/* gas leak angle (degrees)            */

EXTERN double sdc_stefan_flow
(JACOBIAN_VAR_DESC_STRUCT *,    /*  *func_jac        */
        BOUNDARY_CONDITION_STRUCT *,   /*  *bc              */
       int,                           /*  ip               */
       ELEM_SIDE_BC_STRUCT *,         /*  *elem_side_bc    */
       const double [MAX_PDIM],       /*  x_dot[MAX_PDIM]  */
       const double,                  /*  time               */
       const double,                  /*  tt               */
       const double,               /*  dt               */
       const int );               /*  interface_id               */

EXTERN double sdc_stefan_volume_flow
(JACOBIAN_VAR_DESC_STRUCT *,    /*  *func_jac        */
       BOUNDARY_CONDITION_STRUCT *,   /*  *bc              */
       int,                           /*  ip               */
       ELEM_SIDE_BC_STRUCT *,         /*  *elem_side_bc    */
       const double [MAX_PDIM],       /*  x_dot[MAX_PDIM]  */
       const double,                  /*  tt               */
       const double );               /*  dt               */

EXTERN double mass_flux_surface
(JACOBIAN_VAR_DESC_STRUCT *,    /*  *func_jac        */
       const double [MAX_PDIM],       /*  x_dot[MAX_PDIM]  */
       const double,                  /*  tt               */
       const double );               /*  dt               */

extern double vol_flux_surface
(JACOBIAN_VAR_DESC_STRUCT *,    /*  *func_jac        */
       const double [MAX_PDIM],       /*  x_dot[MAX_PDIM]  */
       const double,                  /*  tt               */
       const double );               /*  dt               */

EXTERN void fvelo_normal_edge_bc
(double [DIM],		/* func                                      */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
       const double ,		/* vnormal - normal velocity (speed)         */
       const double [MAX_PDIM],	/* x_dot - mesh velocity vector              */
       const double ,		/* tt - parameter to vary time integration 
				 * from explicit (tt=1) to implicit (tt=0)   */
       const double );		/* dt - current value of the time step       */

EXTERN void fvelo_normal_disc_bc
(double [DIM],		/* func                                      */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
       const double ,		/* vnormal - normal velocity (speed)         */
       const double [MAX_PDIM],	/* x_dot - mesh velocity vector              */
       const double ,		/* tt - parameter to vary time integration 
				 * from explicit (tt=1) to implicit (tt=0)   */
       const double );		/* dt - current value of the time step       */

EXTERN void fvelo_tangent_edge_bc
(double [DIM],		/* func                                      */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
       double *,		/* pvv - velocity of substrate               */
       const double [MAX_PDIM],	/* x_dot - mesh velocity vector              */
       const double ,		/* tt - parameter to vary time integration 
				 * from explicit (tt=1) to implicit (tt=0)   */
       const double );		/* dt - current value of the time step       */

EXTERN void fvelo_tangential_bc
(double [],		/* func                                      */
       double [MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func      */
       double [],		/* x - Soln vector                           */
       const int ,		/* dcl_node - global node number of the 
				 * dynamic contact line (point in 2D)        */
       const double ,		/* vtangent - specified tangential speed 
				 * from input deck                           */
       const double ,		/* beta - slip coefficient that scales the 
				 * xdot slip term                            */
       const double ,		/* alpha - coefficient that scales the 
				 * exponent  of the position dependent 
				 * slip velocity (1/length, from input deck) */
       const double [MAX_PDIM],	/* xsurf - coordinates of surface Gauss point,
				 * i.e. current position                     */
       const double [MAX_PDIM],	/* x_dot - mesh velocity vector              */
       const double ,		/* tt - parameter to vary time integration 
				 * from explicit (tt=1) to implicit (tt=0)   */
       const double ,		/* dt - current value of the time step size  */
       const int ,		/* bc identifier			     */
       const int ,
       double xi[DIM],
       const Exo_DB *exo,
       const double ,                   /* time_value  */
       const double[] ,                 /* parameter_list*/
       const int );                    /* number of parameters  */

EXTERN void fvelo_tangent_3d
(double [],             /* func */
       double [MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],         /* dfunc */
       const double [],       /* free surface velocity x_dot */
       const double ,         /* specified tangential speed */
       const double ,         /* x-component of surface tangent */
       const double ,         /* y-component of surface tangent */
       const double ,         /* z-component of surface tangent */ 
       const double ,	      /* tt - parameter to vary time integration 
				 * method from BE(0) to CN(1/2) to FE(1)     */
       const double );		/* dt - current value of the time step size  */

EXTERN void fzero_velo_tangent_3d
(double [],             /* func */
       double [MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],         /* dfunc */
     const int);


EXTERN void fvelo_tangential_solid_bc
(double [],		/* func                                      */
       double [MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func      */
       double [],               /* solution vector x                         */
       double [],               /* time derivative of solution vector x      */
       double [],               /* time derivative of solution vector x_rs                         */
       const int ,              /* bc type */
       const int ,		/* i_mat_solid                               */
       const int ,		/* i_mat_fluid                               */
       dbl ,                    /* slip coefficient                          */
       const int,               /* node id for DCL                           */
       dbl,               /* position parameter for slip model         */
       const double [MAX_PDIM], /* coordinates for surface gauss point       */
       const dbl ,		/* tt - parameter to vary time integration 
				 * method from BE(0) to CN(1/2) to FE(1)     */
       const dbl );		/* dt - current value of the time step size  */

EXTERN void fvelo_normal_solid_bc
(double [],		/* func                                      */
       double [MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func      */
       double [],               /* solution vector x                         */
       double [],               /* time derivative of solution vector x      */
       double [],               /* time derivative of solution vector x_rs                         */
       const int ,              /* bc type */
       const int ,		/* i_mat_solid                               */
       const int ,		/* i_mat_fluid                               */
        const dbl ,		/* tt - parameter to vary time integration 
				 * method from BE(0) to CN(1/2) to FE(1)     */
       const dbl );		/* dt - current value of the time step size  */

void
fvelo_slip_bc(double func[MAX_PDIM],
	      double d_func[MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
	      double x[],
	      const int type,    /* whether rotational or not */
	      const int max_float,    /* max float number */
              double bc_float[MAX_BC_FLOAT_DATA],
	      const int dcl_node,/*   node id for DCL  */
	      const double xsurf[MAX_PDIM], /* coordinates of surface Gauss  *
					     * point, i.e. current position  */
	      const double tt,   /* parameter in time stepping alg           */
	      const double dt);   /* current time step value                  */

void
fvelo_slip_power_bc(double func[MAX_PDIM],
		    double d_func[MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
		    const int type,    /* whether rotational or not */
		    const int max_float,    /* Max float number from input file */
		    double bc_float[MAX_BC_FLOAT_DATA],
		    const double tt,   /* parameter in time stepping alg           */
		    const double dt);

int
exchange_fvelo_slip_bc_info(int ibc /* Index into BC_Types for VELO_SLIP_BC */);

EXTERN void
fvelo_slip_ls_heaviside(double func[MAX_PDIM],
			double d_func[MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
			double width,
			double beta_negative,
			double beta_positive,
			const double vsx,      /* velocity components of solid  */
			const double vsy,	/* surface on which slip condition   */
			const double vsz,	/* is applied           */
			const double tt,
                        const double dt);

void
fvelo_airfilm_bc(double func[MAX_PDIM],
	      double d_func[MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
	      double x[],
	      const int type,    /* whether rotational or not */
              double bc_float[MAX_BC_FLOAT_DATA],
	      const int dcl_node,/*   node id for DCL  */
	      const double xsurf[MAX_PDIM], /* coordinates of surface Gauss  *
					     * point, i.e. current position  */
	      const double tt,   /* parameter in time stepping alg           */
	      const double dt);   /* current time step value                  */

EXTERN void fvelo_slip_level
( double [MAX_PDIM],	/* func                                      */
	double [MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func      */
        const int ,		/* type - whether rotational or not          */
	double,                 /* width   */
	double,                 /* beta    */
	const double ,		/* vsx                                       */
	const double ,		/* vsy                                       */
	const double ,		/* vsz - velocity components of solid surface*/
	const double ,        	/* beta for outside of interface region      */
	const double ,		/* gas phase factor                          */
	const double ,          /* contact fraction*/
	const double ,          /* time stabilization factor*/
        const double , 	        /* tt - parameter to vary time integration */
        const double );	/* dt - current value of the time step size  */	
				   
EXTERN void fvelo_slip_electrokinetic_bc
( double [MAX_PDIM],	/* func                                      */
	double [MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func      */
	double ,	/* permitivity       */
	double );	/* zeta_potential    */

EXTERN void fvelo_electrokinetic_3d
(double [],             /* func */
       double [MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],         /* dfunc */
       const double [],       /* free surface velocity x_dot */
       const double ,         /* permitivity       */
       const double ,	      /* zeta_potential    */
       const double ,         /* x-component of surface tangent */
       const double ,         /* y-component of surface tangent */
       const double ,         /* z-component of surface tangent */ 
       const double ,	      /* tt - parameter to vary time integration */
       const double );		/* dt - current value of the time step size  */	


EXTERN void load_surface_tension
(double [][MDE]);	/* dsigma_dx - dimensions [DIM][MDE] */

EXTERN void elec_surf_stress
(double [MDE][DIM],	/* cfunc                                     */
       double [MDE][DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_cfunc     */
       const int ,		/* id_side - ID of the side of the element   */
       struct elem_side_bc_struct *, /* elem_side_bc                         */
       const int ,		/* iconnect_ptr                              */
       const double,            /* Electric Permittivity                     */
       const double,          /* Equation term multiplier                  */
       const int);          /* BC type                 */

EXTERN void fn_dot_T
(double [MDE][DIM],	/* cfunc                                     */
       double [MDE][DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_cfunc     */
       const int ,		/* id_side - ID of the side of the element   */
       const double ,		/* sigma - surface tension                   */
       const double ,		/* pb - applied pressure                     */
       struct elem_side_bc_struct *, /* elem_side_bc                         */
       const int ,		/* iconnect_ptr                              */
       double [DIM][MDE]);	/* dsigma_dx                               */

EXTERN void apply_repulsion
(double cfunc[MDE][DIM],
       double d_cfunc[MDE][DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
       const double ,		/* pr - coefficient for repulsion force to 
				 * ensure no penetration of the solid boundary
				 * by the free surface                       */
       const double ,		/* "A" in Ax + By + Cz + D = 0 equation      */
       const double ,		/* "B" in Ax + By + Cz + D = 0 equation      */
       const double ,		/* "C" in Ax + By + Cz + D = 0 equation      */
       const double ,		/* "D" in Ax + By + Cz + D = 0 equation      */
       struct elem_side_bc_struct *, /* elem_side_bc */
       const int );		/* iconnect_ptr */

EXTERN void apply_repulsion_roll
(double cfunc[MDE][DIM],
       double d_cfunc[MDE][DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
       double [],		/* Solution vector */
       const double ,		/* roll radius      */
       const double [3],		/* axis origin      */
       const double [3],		/* direction angles      */
       const double ,		/* omega - roll rotation rate    */
       const double ,		/* repulsion length scale      */
       const double ,		/* repulsion exponent     */
       const double ,		/* repulsion coefficient     */
       const double ,		/* gas viscosity   */
       const double ,		/* exclusion scale    */
       const int ,		/* DCL node id    */
       struct elem_side_bc_struct *, /* elem_side_bc */
       const int );		/* iconnect_ptr */

EXTERN void apply_repulsion_user
(double cfunc[MDE][DIM],
       double d_cfunc[MDE][DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
       const double ,		/* roll radius      */
       const double [3],		/* axis origin      */
       const double [3],		/* direction angles      */
       const double ,		/* repulsion length scale      */
       const double ,		/* repulsion exponent     */
       const double ,		/* repulsion coefficient     */
       const double ,		/* inverse slip coefficient    */
       const double ,		/* omega - roll rotation rate    */
       struct elem_side_bc_struct *, /* elem_side_bc */
       const int );		/* iconnect_ptr */

EXTERN void apply_repulsion_table
(double cfunc[MDE][DIM],
       double d_cfunc[MDE][DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
       double [],		/* Solution vector */
       const double ,		/* repulsion length scale      */
       const double ,		/* repulsion exponent     */
       const double ,		/* repulsion coefficient     */
       const double ,		/* inverse slip coefficient    */
       const double ,		/* exclusion scale    */
       const double [3],		/* wall velocity      */
       const int ,		/* DCL node id    */
       struct elem_side_bc_struct *, /* elem_side_bc */
       const int );		/* iconnect_ptr */


EXTERN void apply_vapor_recoil
(double cfunc[MDE][DIM],
       double d_cfunc[MDE][DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
       const double ,		/* Boiling Temp                              */
       const double ,		/* Melting temp                              */
       const double ,		/* Reference temp                            */
       const double ,		/* Pressure scale                            */
       const double ,    /* Temperature scale */
       struct elem_side_bc_struct *, /* elem_side_bc */
       const int );		/* iconnect_ptr */

EXTERN void flow_n_dot_T_hydro
(double [DIM],		/* func                                      */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
       const double ,		/* a for pressure variation                  */
       const double ,		/* b for pressure variation                  */
       const double ,		/* c for pressure variation                  */
       const double );		/* d - pressure variation                    */

EXTERN void flow_n_dot_T_var_density
(double [DIM],            /* func                                      */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],/* d_func            */
       const double,            /* a - reference pressure                  */
       const double);          /* time - current time                       */

#if 0
/* deprecated March 2002 by TAB */
EXTERN void hydrostatic_n_dot_T
(double *,		/* func                                      */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE]); /* d_func         */
#endif

EXTERN void flow_n_dot_T_nobc
(double [DIM],		/* func                                      */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
       const double ,		/* pdatum - pressure datum from input card   */
       const int );		/* iflag - 1 to use pdatum, otherwise use P  */

EXTERN void flow_n_dot_T_gradv
(double [DIM],		/* func                                      */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
       const double ,		/* pdatum - pressure datum from input card   */
       const int );		/* iflag - -1 to use pdatum, otherwise use P */

EXTERN void flow_n_dot_T_segregated
(double [DIM],		                            // func
 double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE]); // d_func

EXTERN void
stress_no_v_dot_gradS(double func[MAX_MODES][6],
                      double d_func[MAX_MODES][6][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                      const double dt,
                      const double tt
			   );

EXTERN void
stress_no_v_dot_gradS_logc(double func[MAX_MODES][6],
                      double d_func[MAX_MODES][6][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                      const double dt,
                      const double tt
			   );

EXTERN void flow_n_dot_T_gradv_sic(double [DIM],                  /* func  */
				   double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],  /* d_func   */
				   const double,     /* pdatum - pressure datum from input card   */
				   const int);       /* iflag - -1 to use pdatum, otherwise use P */

EXTERN void press_poisson_segregated
(double * ,                                        // func 
 double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], // d_func
 double ,                                          // current time
 double);                                          // time step

EXTERN void
PSPG_consistency_bc (double *func,
		     double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
		     const dbl x_dot[DIM],
		     const dbl time, /* current time  */
		     const dbl dt, /* time step size */
		     const dbl tt, /* time step parameter */
		     const PG_DATA *pg_data);	/* U_norm - global velocity norm             */


EXTERN void fapply_CA
(double *,		/* func                                      */
       double [MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func      */
       double [MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func_ss   */
       const double [MAX_PDIM],	/* fsnormal - free surface normal components */
       double [MAX_PDIM][MAX_PDIM][MDE] ,			 /* dfsnormal_dx - free surface 
					  * normal derivatives 
					  * ([i][j][k]) component i wrt 
					  * displacement j at node k         */
       const double [MAX_PDIM],	/* ssnormal - solid surface normal components*/
       double [MAX_PDIM][MAX_PDIM][MDE],			 /* dssnormal_dx - solid surface 
					  * normal derivatives 
					  * ([i][j][k]) component i wrt 
					  * displacement j at node k         */
       const double  );	/* contact_angle - Static or dynamic contact 
				 * angle                                     */

EXTERN void fapply_var_CA
(double *,		/* func                                      */
       double [MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func      */
       double [MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func_ss   */
       const double [MAX_PDIM],	/* fsnormal - free surface normal components */
       double [MAX_PDIM][MAX_PDIM][MDE], /* dfsnormal_dx - free surface 
						* normal derivatives 
						* ([i][j][k]) component i wrt 
						* displacement j at node k   */
       const double [MAX_PDIM],	/* ssnormal - solid surface normal components*/
       double [MAX_PDIM][MAX_PDIM][MDE], /* dssnormal_dx - solid surface 
						* normal derivatives 
						* ([i][j][k]) component i wrt 
						* displacement j at node k   */
       const double [MAX_PDIM],	/* clnormal - contact line normal components */

       double [MAX_PDIM][MAX_PDIM][MDE], /* dclnormal_dx - contact line
						* normal vector derivatives
						* ([i][j][k]) component i wrt
						* displacement j at node k   */
       const dbl [5],		/* BC_Data_Float - static contact angle, 
				 * response slope, components of web velocity
				 * vector                                    */
       const double [MAX_PDIM],	/* xdot - Current mesh velocity vector       */
       const double ,		/* tt - parameter varies to select the
				 * time integration scheme from
				 * BE(0) to CN(1/2) to FE(1)                 */
       const double);		/* dt - current time step size               */

EXTERN void fapply_var_CA_user
(double *,		/* func                                      */
       double [MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func      */
       double [MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func_ss   */
       const double [MAX_PDIM],	/* fsnormal - free surface normal components */
       double [MAX_PDIM][MAX_PDIM][MDE], /* dfsnormal_dx - free surface 
						* normal derivatives 
						* ([i][j][k]) component i wrt 
						* displacement j at node k   */
       const double [MAX_PDIM],	/* ssnormal - solid surface normal components*/
       double [MAX_PDIM][MAX_PDIM][MDE], /* dssnormal_dx - solid surface 
						* normal derivatives 
						* ([i][j][k]) component i wrt 
						* displacement j at node k   */
       const double [MAX_PDIM],	/* clnormal - contact line normal components */

       double [MAX_PDIM][MAX_PDIM][MDE], /* dclnormal_dx - contact line
						* normal vector derivatives
						* ([i][j][k]) component i wrt
						* displacement j at node k   */
       const int ,		/* num_user_const - number of user constants 
				 * in the user-defined model                 */
       const double *,		/* user_const - ptr to array of the 
				 * user constants                            */
       const double [MAX_PDIM],	/* xdot - Current mesh velocity vector       */
       const double ,		/* tt - parameter varies to select the
				 * time integration scheme from
				 * BE(0) to CN(1/2) to FE(1)                 */
       const double);		/* dt - current time step size               */

EXTERN int evaluate_gibbs_criterion
(const double [MAX_PDIM], /* fsnormal - free surface normal vector     */
       const double [MAX_PDIM], /* ssnormal - solid surface normal vector    */
       int *ipin,		/* Flag indicates pinned or not        (out) */
       const double ,		/* contact_angle - Static or dynamic contact 
				 * angle                                     */
       const double ,		/* x_pos - x-coordinate of sharp edge        */
       const double ,		/* y_pos - y-coordinate of sharp edge        */
       const double ,		/* z_pos - z-coordinate of sharp edge        */
       const double ,		/* sign_origx - original relative sign on the
				 * x-position of contact line with specified 
				 * sharp edge                                */
       const double ,		/* sign_origy - original relative sign on the
				 * y-position of contact line with specified
				 * sharp edge                                */
       const double);		/* sign_origz - original relative sign on the
				 * z-position of contact line with specified 
				 * sharp edge                                */


EXTERN void fapply_moving_CA
(double *,		/* func                                      */
       double [MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func      */
       double [MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func_ss   */
       const double [MAX_PDIM],	/* fsnormal - free surface normal components */
       double [MAX_PDIM][MAX_PDIM][MDE], /* dfsnormal_dx - free surface 
					  * normal derivatives 
					  * ([i][j][k]) component i wrt 
					  * displacement j at node k   */
       const double [MAX_PDIM],	/* ssnormal - solid surface normal components*/
       double [MAX_PDIM][MAX_PDIM][MDE], /* dssnormal_dx - solid surface 
						* normal derivatives 
						* ([i][j][k]) component i wrt 
						* displacement j at node k   */
       const double ,		/* stat_ca - Static or dynamic contact angle */
       const double ,           /* advancing_ca - Static or dynamic contact 
				 * angle                                     */
       const double ,		/* receding_ca - Static or dynamic contact 
				 * angle                                     */
       const double ,		/* scaling - Static or dynamic contact angle */
       const double ,		/* vwx - x-Velocity of wall                  */
       const double ,		/* vwy - y-Velocity of wall                  */
       const double ,		/* vwz - z-Velocity of wall                  */

       const double [MAX_PDIM],	/* xdot - Current mesh velocity vector       */
       const double ,		/* dt - current time step size               */
       const double);		/* tt - parameter varies to select the
				 * time integration scheme from
				 * BE(0) to CN(1/2) to FE(1)                 */

EXTERN void fapply_moving_CA_sinh
(double *,		/* func */
       double [MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func */
       double [MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func_ss */
       const double [MAX_PDIM],	/* fsnormal - free surface normal components */
       double [MAX_PDIM][MAX_PDIM][MDE], /* dfsnormal_dx - free surface 
						* normal derivatives 
						* ([i][j][k])    
						* component i wrt            *
						* displacement j at node k   */
       const double [MAX_PDIM],	/* ssnormal - solid surface normal vector    */
       double [MAX_PDIM][MAX_PDIM][MDE], /* dssnormal_dx - solid surface
						* normal vector sensitivity
						* derivatives ([i][j][k])    *
						* component i wrt            *
						* displacement j at node k   */
       const double ,		/* equilibrium_contact_angle */
       const double ,		/* velocity_pre_exponential */
       const double ,		/* modified_surface_tension */
       const double ,		/* relaxation_time */
       const double ,		/* old_tpl_velocity */
       const double [MAX_PDIM],	/* x_dot - mesh velocity vector              */
       const dbl , /* dt - current value of the time step            */
       const dbl ,     /* tt - explicit (tt = 1) -- implicit (tt = 0)    */
       const double ,		/* time - current time */
       const double ,		/* wall velocity */
       const double ,		/* theta_max */
       const double ,		/* dewet parameter */
       const double ,		/* dcl_shearrate */
       const int ,		/* BC identifier */
       double [MAX_PDIM][MDE],          /* wall velo derivs     */
       const int );            /* local_node_number    */




EXTERN void apply_ST
(const int ,		/* irow_index - Row index for Elemental 
				 * stiffness matrix                          */
       const int ,		/* I - Global node number                    */
       const int ,		/* iconnect_ptr - Pointer for element into 
				 * connectivity array                        */
       const int ,		/* ielem_dim - Physical dimension of element,
				 * i.e., 1, 2, 3                             */
       int [],			/* ija - column pointer array                */
       double [],		/* a - nonzero elements of matrix            */
       double [],		/* resid_vector - Residual vector	     */
       double [],		/* x - Vector containing current solution    */
       const int ,		/* num_nodes_on_side - Number of nodes on 
				 * the side of the element                   */
       const int [MAX_NODES_PER_SIDE], /* local_elem_node_id - Vector of 
					* local elem node numbers on the 
					* element side                       */
       const int ,		/* id_side - ID of the side of the element   */
       const double ,		/* tx - x-component surface tangent outflow  */
       const double ,		/* ty - y-component surface tangent outflow  */
       const double ,		/* tz - z-component surface tangent outflow  */
       const double ,		/* sigma - surface tension                   */
       double ,			/* rcoord - r coord in axisym problems: 
				 * used in r dr dz                           */
       const int ,              /* ufixed - u-flag for fixing x-component of 
				 * velocity                                  */
       const int ,              /* vfixed - v-flag for fixing y-component of 
				 * velocity                                  */
       const int);		/* wfixed - w-flag for fixing z-component of 
				 * velocity                                  */

EXTERN void fapply_ST
(double [MAX_PDIM],	/* func                                      */
       double [MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func      */
       const double ,		/* tx - x-component surface tangent outflow  */
       const double ,		/* ty - y-component surface tangent outflow  */
       const double ,		/* tz - z-component surface tangent outflow  */
       const double ,		/* sigma - surface tension                   */
       const double ,           /* rcoord - r coord in axisym problems: used 
				 * in r dr dz                                */
       const int );		/* id_side */

EXTERN void apply_ST_scalar
(double [MAX_PDIM],	/* func                                      */
       double [MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func      */
       const int ,		/* id_side HKM -> useless                    */
       const double, 		/* sigma - surface tension                   */
       const int );		/* User added sign -> handles unresolved     *
                                 * sign issues for this card                 */

EXTERN void apply_ST_3D
(double [MAX_PDIM],	/* func                                      */
       double [MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func      */
       const double ,		/* tx - x-component surface tangent outflow  */
       const double ,		/* ty - y-component surface tangent outflow  */
       const double ,		/* tz - z-component surface tangent outflow  */
       const double );		/* sigma - surface tension                   */

EXTERN void apply_ST_scalar_3D
(double [MAX_PDIM],	/* func                                      */
       double [MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func      */
       const double );		/* sigma - surface tension                   */

EXTERN void apply_CA_FILL
(double [MAX_PDIM],	/* func                                      */
       double [MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func      */
       const double );         /* contact angle                             */

EXTERN void apply_sharp_ca
(double [MAX_PDIM],	/* func                                      */
       double [MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func      */
       const double );         /* contact angle                             */

EXTERN void apply_wetting_velocity
(	double [MAX_PDIM],
		double [MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], 
		int ,
		int ,
		double ,
		double [MAX_PDIM],
		double, 
		double ,
		double  );


void apply_linear_wetting_sic
( double [MAX_PDIM],
	double [MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], 
	const double ,
	const double ,
	int ,
	int ,
	double ,
	double [MAX_PDIM],
	double ,
	double ,
	double ,
	double );
		
EXTERN void apply_wetting_tension
(double [MAX_PDIM],	/* func                                      */
       double [MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func      */
       const double );         /* wetting tension                           */

EXTERN void velocity_profile
(const int ,		/* irow_index - Elemental stiffness matrix 
				 * row index                                 */
       const int ,		/* I - Global node number                    */
       const int ,		/* iconnect_ptr - pointer for this element 
				 * into the elem-node connectivity           */
       const int ,		/* ielem_dim - physical dimension of the 
				 * element, ie., 1, 2, 3                     */
       int [],			/* ija - column pointer array                */
       double [],		/* a - nonzero matrix elements               */
       double [],		/* resid_vector - Residual vector	     */
       double [],		/* x - Vector containing current solution    */
       const int ,		/* velo_condition - integer denoting which 
				 * condition is being applied                */
       const int ,		/* num_nodes_on_side - Number of nodes on 
				 * the side of the element                   */
       const int [MAX_NODES_PER_SIDE], /* local_elem_node_id - local 
					* element's node numbers located on 
					* the side of the element 	     */
       const double ,		/* a1 - parameters from input deck card      */
       const double ,		/* a2 - parameters from input deck card      */
       const double );		/* time - at which bc's are evaluated        */

EXTERN void ftmelt_bc
(double [DIM],		/* func                                      */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
       const double ,		/* a1 - function parameters from data card   */
       const double [MAX_PDIM],	/* x_dot -  mesh velocity vector             */
       const dbl ,		/* tt - parameter to vary time integration 
				 * from explicit (tt = 1) to 
				 * implicit (tt = 0)                         */
       const dbl );		/* dt - current value of the time step size  */


EXTERN void continuous_tangent_velocity
(double [DIM],		/* func                                      */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
       const int );		/* ielem_dim - physical dimension of the 
				 * element, ie., 1, 2, 3                     */

EXTERN void continuous_normal_velocity
(double [DIM],		/* func                                      */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
       const int );		/* ielem_dim - physical dimension of the 
				 * element, ie., 1, 2, 3                     */

EXTERN void discontinuous_velocity
(double [DIM],		 /* func                                      */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
       double [MAX_PDIM],	 /* x_dot - mesh velocity vector              */
       const int ,		 /* mode - Evporation or dissolution          */
       const int ,		 /* i_mat_liquid - the blk id for the liquid  */
       const int ,		 /* i_mat_gas - the blk id for the gas phase  */
       const dbl ,		 /* tt - parameter to vary time integration 
				  * from explicit (tt = 1) to 
				  * implicit (tt = 0)                         */
       const dbl );		 /* dt - current value of the time step size  */

EXTERN void fnormal_stress_bc
(double [DIM],		  /* func                                      */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
       const double ,		  /* stress_normal -             normal stress */
       const dbl );		  /* relax - relaxation parameters             */

EXTERN void qside_directional
(double [DIM],              /* func */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /*   d_func   */
       const double ,             /* value of qx  */
       const double ,             /* value of qy  */
       const double );           /* value of qz  */
EXTERN void qside_contact_resis
(double [DIM],              /* func */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /*   d_func   */
       const int ,                /* value id_block_1  */
       const int ,                /* value id_block_2  */
       const double );           /* value Rinv  */

EXTERN void qside_light_jump
(double [DIM],              /* func */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /*   d_func   */
       const double,	/* Time                                      */
       const int,	/* bc_type */     
       const int ,                /* value id_block_1  */
       const int                 /* value id_block_2  */
        );           /* value Rinv  */

EXTERN void qside_ls
(double [DIM],              /* func */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /*  d_func */
       int ,                      /* no_LS_block_id */
       int ,                      /* LS_block_id */
       double ,                   /* Q_neg_side  */
       double );                 /* Q_pos_side  */

EXTERN void q_velo_slip_bc	/* mm_fill_terms.c                           */
(double [MAX_PDIM],	/* func                                      */
       double [MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func      */
       int,                     /* corresponding slip bc id                  */
       double [MAX_PDIM],       /* soln vector                               */
       const double [MAX_PDIM], /*gauss point coordinates                  */
       const double,
       const double);
               
EXTERN void qnobc_surf		/* mm_fill_terms.c                           */
(double [DIM],		/* func                                      */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
       const double);          /* Time                                      */
 
EXTERN void potential_nobc_surf	
(double [DIM],		/* func                                      */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
       const double);          /* Time                                      */

EXTERN void qlaser_surf		/* mm_fill_terms.c                           */
(double [DIM],		/* func                                      */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
       const double [],		/* p- parameters for the bc                  */
       const double [],		/* goma solution vector          	     */
       const double);          /* Time                                      */

EXTERN void q_vapor		/* mm_fill_terms.c                           */
(double [DIM],		/* func                                      */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
       const double [],		/* p- parameters for the bc                  */
       const double [],		/* goma solution vector          	     */
       const double);          /* Time                                      */

EXTERN double calculate_laser_flux
( const double [],
        double,
        const double [],
	double [],
	int,
	int,
	double []);

EXTERN double calculate_vapor_cool
( const double [],
        double *,
        double);

EXTERN void qrad_surf		/* mm_fill_terms.c                           */
(double [DIM],		/* func                                      */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
       double ,			/* heat_tran_coeff - (cgs units)             */
       double ,			/* T_c - bath temperature (Kelvin)	     */
       double ,			/* emissivity                                */
       double );		/* Boltzmann's constant                      */

EXTERN void qrad_surf_repulse		/* mm_ns_bc.c                           */
(double [DIM],		/* func                                      */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
       const double ,			/* heat_tran_coeff - (cgs units)             */
       const double ,			/* T_c - bath temperature (Kelvin)	     */
       const double ,			/* emissivity                                */
       const double ,		/* Boltzmann's constant                      */
       const double ,		/* roll radius      */
       const double [3],		/* axis origin      */
       const double [3],		/* direction angles      */
       const double ,		/* repulsion length scale      */
       const double ,		/* repulsion exponent     */
       const double ,		/* repulsion coefficient     */
       const double );         /* Roll Temperature    */

EXTERN void apply_sharp_wetting_velocity
(	double [MAX_PDIM],
		double [MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], 
		const int,
		const double ,
		const double ,
		const double ,
		const double ,
		const double ,
		const double  );


void apply_blake_wetting_velocity
(	double [MAX_PDIM],
		double [MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], 
		const int,
		int ,
		int ,
		double ,
		double ,
		double ,
		double ,
		double  );
                
void apply_blake_wetting_velocity_sic
(	double [MAX_PDIM],
		double [MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], 
                const double,
                const double,
		      int,
		const int ,
		const int ,
		double ,
		const double ,
		const double ,
		const double ,
		const double ,
		      double [MAX_PDIM],
		const double,
                const int ,
		const double ,
		const double ,
		const double );

EXTERN void acoustic_plane_transmission	
(double [DIM],		/* func                                      */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
       const double,	/* Time                                      */
       const int,	/* bc_type */     
       const double,    /*  boundary impedance	*/
       const double,    /*  boundary absorption	*/
       const double,    /*  boundary incident real*/
       const double,    /*  boundary incident imaginary	*/
       const int );	/* element block id */     

EXTERN void light_transmission	
(double [DIM],		/* func                                      */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
       const double,	/* Time                                      */
       const int,	/* bc_type */     
       const double,    /*  boundary impedance	*/
       const double,    /*  boundary absorption	*/
       const double,    /*  boundary incident */
       const int );	/* element block id */     

EXTERN void acoustic_nobc_surf	
(double [DIM],		/* func                                      */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
       const double,	/* Time                                      */
       const int );     /* bc_type */     

void sheet_tension 
( double [MDE][DIM],
	double [MDE][DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
	const int ,	/* ID of the side of the element             */
	const double ,
	struct elem_side_bc_struct *,
	const int ,
        double [DIM],
        const Exo_DB * );


void apply_SES 
( double [MAX_PDIM],
	 double [MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
	 struct elem_side_bc_struct *,
	 double ,
	 int ,
	 double ,
	 double ,
	 int, 
	 int ,
	 int );


void shear_to_shell 
( double [MDE][DIM],
		double [MDE][DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
		const int ,	/* ID of the side of the element */
		double ,
		struct elem_side_bc_struct *,
		const int ,
		double [DIM],
		const Exo_DB *);
		
void apply_hysteresis_wetting_sic 
(	double *, 			
		double [MAX_PDIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], 
		double ,
		double ,
		double * );

void fgamma1_deriv_bc
(double [DIM],		/* func                                      */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
       const double );		/* vnormal - normal velocity                 */
  
void fgamma2_deriv_bc
(double [DIM],		/* func                                      */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
       const double );		/* vnormal - normal velocity                 */

void dvzdr_zero_deriv_bc
(double [DIM],		/* func                                      */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
       const double nwall[DIM],
       const double );		/* vnormal - normal velocity                 */
  
  
void
ls_wall_angle_bc(double func[DIM],
                double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
                const double angle);
#endif /* GOMA_MM_NS_BC_H */
