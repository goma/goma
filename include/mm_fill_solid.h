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
 

#ifndef GOMA_MM_FILL_SOLID_H
#define GOMA_MM_FILL_SOLID_H

#include "el_elm.h"
#include "exo_struct.h"
#include "mm_fill_rs.h"
#include "rf_fem_const.h"
#include "std.h"

struct Boundary_Condition;
struct Elastic_Constitutive;
#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_MM_FILL_SOLID_C
#define EXTERN /* do nothing */
#endif

#ifndef GOMA_MM_FILL_SOLID_C
#define EXTERN extern
#endif

EXTERN int belly_flop		/* mm_fill_solid.c                           */
(dbl );			/* mu - elastic modulus (plane stress case)  */

EXTERN void invert_tensor	/* mm_fill_solid.c                           */
(double [DIM][DIM],	/* A - tensor to be inverted                 */
       double [DIM][DIM],	/* B - inverted tensor                       */
       int ,			/* dim - dimensions of tensor                */
       double [DIM][DIM][DIM][MDE], /* dA - sensitivities of tensor to be 
				     * inverted                              */
       double [DIM][DIM][DIM][MDE], /* dB - sensitivities of inverted 
				     * tensor                                */
       int ,			/* dof - number of dofs of variable for 
				 * sensitivities                             */
       int );			/* sense - flag to calculate sensitivities   */

EXTERN void slope_n_dot_n0_bc	/* mm_fill_solid.c                           */
(double [DIM],		/* func                                      */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
       double ,			/* slopex                                    */
       double ,			/* slopey                                    */
       double );		/* slopez                                    */

EXTERN void force_n_dot_f_bc	/* mm_fill_solid.c                           */
(double [DIM],		/* func                                      */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
       const double ,		/* forcex                                    */
       const double ,		/* forcey                                    */
       const double ,		/* forcez                                    */
       const int ,		/* sic_flag                                  */
       const double ,		/* delta_t         */
       const double ,		/* theta         */
       const int ,		/* ip         */
       const int ,		/* ip_total         */
       const double );		/* time_value  */

EXTERN void rep_force_n_dot_f_bc /* mm_fill_solid.c                          */
(double [DIM],		/* func                                      */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
       const double ,		/* pr - coeff for repulsion force to ensure
				 * no penetration of the solid boundary by the
				 * free surface                              */
       const double ,		/* ap - a coefficient in plane equation      */
       const double ,		/* bp - b coefficient in plane equation      */
       const double ,		/* cp - c coefficient in plane equation      */
       const double ,		/* dp - d coefficient in plane equation      */
       const double ,		/* repexp - repulsive force exponent      */
       const double ,		/* friction - friction coefficient      */
       const int );		/* BC id #      */

EXTERN void rep_force_roll_n_dot_f_bc /* mm_fill_solid.c                          */
(double [DIM],		/* func                                      */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
       const double ,		/* pr - coeff for repulsion force to ensure
				 * no penetration of the solid boundary by the
				 * free surface                              */
       const double [3],		/* ap - a coefficient in plane equation      */
       const double [3],		/* bp - b coefficient in plane equation      */
       const double ,		/* roll_rad - radius of cylindrical indentor */
       const double ,		/* repexp - repulsive force exponent      */
       const double ,		/* friction - friction coefficient      */
       const int );		/* BC id #      */

EXTERN void norm_force_n_dot_f_bc /* mm_fill_solid.c                         */
(double [DIM],		/* func                                      */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
       double ,			/* forcex                                    */
       double ,			/* forcey                                    */
       double );		/* forcez                                    */

EXTERN void put_liquid_stress_in_solid /* mm_fill_solid.c                    */
(int ,			/* id - local element node number for the 
				 * current node whose residual contribution
				 * is being sought                           */
       int ,			/* I - Global node number                    */
       int ,			/* ielem_dim - physical dimension of element,
				 * ie., 1, 2, 3                              */
       double [],		/* resid_vector - Residual vector NO DUH!    */
       int ,			/* i_mat_solid - mat block id's of solid     */
       int ,			/* i_mat_fluid - mat block id's of liquid    */
       int [],			/* local_node_list_fs - MDE list to keep track
				 * of nodes at which liquid contributions have
				 * been transfered to solid (fluid-solid 
				 * boundaries)                               */
       double );		/* scale - term scaling                      */

EXTERN void put_liquid_stress_in_solid_ALE /* mm_fill_solid.c                */
(int ,			/* id - local element node number for the 
				 * current node whose residual contribution
				 * is being sought                           */
       int ,			/* I - Global node number                    */
       int ,			/* ielem_dim - physical dimension of element,
				 * ie., 1, 2, 3                              */
       double [],		/* resid_vector - Residual vector NO DUH!    */
       int ,			/* i_mat_solid - mat block id's of solid     */
       int ,			/* i_mat_fluid - mat block id's of liquid    */
       int [],			/* local_node_list_fs - MDE list to keep track
				 * of nodes at which liquid contributions have
				 * been transfered to solid (fluid-solid 
				 * boundaries)                               */
       double );		/* scale - term scaling                      */

EXTERN void put_fluid_stress_on_shell /* mm_fill_shell.c                    */
(int ,			/* id - local element node number for the 
				 * current node whose residual contribution
				 * is being sought                           */
       int ,                    /* Local shell element node num associated   */
                                /* with id */
       int ,			/* I - Global node number                    */
       int ,			/* ielem_dim - physical dimension of element,
				 * ie., 1, 2, 3                              */
       double [],		/* resid_vector - Residual vector NO DUH!    */
       int [],			/* local_node_list_fs - MDE list to keep track
				 * of nodes at which liquid contributions have
				 * been transfered to solid (fluid-solid 
				 * boundaries)                               */
       double );		/* scale - term scaling                      */
	   
EXTERN void  put_shear_stress_on_shell
(int ,		/* local element node number for the 
					 * current node whose residual contribution
					 * is being sought                        */
		int ,		/* local shell element node number corresponding to id */
		int ,		/* Global node number                      */
		int ,		/* physical dimension of the elem  */
		int [],		/* MDE list to keep track
					 * of nodes at which 
					 * solid contributions 
					 * have been transfered
					 * to liquid (fluid-solid
					 * boundaries)          */
		double ); /* Scale factor, nondimension       */

EXTERN void penetration		/* mm_fill_solid.c                           */
(double [],		/* func                                      */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
       double [DIM],		/* x_dot - mesh velocity vector              */
       dbl ,			/* tt - parameter to vary time integration 
				 * from explicit (tt = 1) to 
				 * implicit (tt = 0)                         */
       dbl ,			/* dt - current value of the time step       */
       int ,			/* bc_input_id                               */
       struct Boundary_Condition *, /* BC_Types -                            */
       int ,			/* i_mat_solid                               */
       int );			/* i_mat_fluid                               */

EXTERN void no_slip		/* mm_fill_solid.c                           */
(double [],		/* func                                      */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
       double [DIM],		/* x_dot - mesh velocity vector              */
       double [DIM],		/* x_rs_dot - mesh velocity vector              */
       dbl ,			/* tt - parameter to vary time integration 
				 * from explicit (tt = 1) to 
				 * implicit (tt = 0)                         */
       dbl ,			/* dt - current value of the time step       */
       int ,			/* bc_input_id                               */
       struct Boundary_Condition *, /* BC_Types -                            */
       int ,			/* i_mat_solid                               */
       int );			/* i_mat_fluid                               */

EXTERN int mesh_stress_tensor	/* mm_fill_solid.c                           */
(dbl [DIM][DIM],		/* TT                                        */
       dbl [DIM][DIM][DIM][MDE], /* dTT_dx                                   */
       dbl [DIM][DIM][MDE],	/* dTT_dp                                    */
       dbl [DIM][DIM][MAX_CONC][MDE], /* dTT_dc                              */
       dbl [DIM][DIM][MDE],      /* dTT_dp_liq                               */
       dbl [DIM][DIM][MDE],      /* dTT_dp_gas                               */
       dbl [DIM][DIM][MDE],      /* dTT_dporosity                            */
       dbl [DIM][DIM][MDE],      /* dTT_dsink_mass                           */
       dbl [DIM][DIM][MDE],      /* dTT_dT                                   */
       dbl [DIM][DIM][MDE],      /* dTT_dmax_strain                          */
       dbl [DIM][DIM][MDE],      /* dTT_dcur_strain                          */
       dbl ,			/* mu                                        */
       dbl ,			/* lambda                                    */
       dbl ,			/* delta_t                                   */
       int ,			/* ielem - current element number            */
       int ,			/* ip - current integration point            */
       int );			/* ip_total - num gauss points in element    */

EXTERN int get_evp_stress_tensor /* mm_fill_solid.c                          */
(dbl [DIM][DIM],		/* TT                                        */
       dbl [DIM][DIM][DIM][MDE], /* dTT_dx                                   */
       dbl [DIM][DIM][MDE],	/* dTT_dp                                    */
       dbl [DIM][DIM][MAX_CONC][MDE], /* dTT_dc                              */
       dbl ,			/* mu                                        */
       dbl ,			/* lambda                                    */
       dbl ,			/* delta_t                                   */
       int ,			/* ielem - current element number            */
       int ,			/* ip - current integration point            */
       int );			/* ip_total - num gauss points in element    */

EXTERN int get_F_vp		/* mm_fill_solid.c                           */
(dbl [DIM][DIM],		/* F_vp                                      */
       dbl [DIM][DIM],		/* F_vp_old                                  */
       dbl [DIM][DIM],		/* TT                                        */
       dbl [DIM][DIM][DIM][MDE],/* dTT_dx				     */
       dbl [DIM][DIM][MAX_CONC][MDE],/* dTT_dc				     */
       dbl [DIM][DIM][DIM][MDE],/* dF_vp_dx				     */
       dbl [DIM][DIM][MAX_CONC][MDE],/* dF_vp_dc			     */
       dbl [MAX_CONC][MDE],	/* d_plastic_mu_dc			     */
       dbl [MAX_CONC][MDE],	/* d_yield_dc				     */
       dbl );			/* delta_t                                   */

EXTERN int load_elastic_properties /* mm_fill_solid.c                        */
(struct Elastic_Constitutive *,
       double *,	       	/* mu                                        */
       double *,		/* lambda                                    */
       double * ,		/* thermexp                                  */
       double [MAX_CONC] ,	/* speciesexp                                */
       double [DIM][MDE],	/* d_mu_dx                                   */
       double [DIM][MDE],	/* d_lambda_dx                               */
       double [MAX_VARIABLE_TYPES+MAX_CONC],	/* d_thermexp_dx             */
       double [MAX_CONC][MAX_VARIABLE_TYPES+MAX_CONC]);/* d_speciesexp_dx   */

EXTERN int load_plastic_properties /* mm_fill_solid.c                        */
(double [MAX_CONC][MDE],	/* d_plastic_mu_dc                           */
       double [MAX_CONC][MDE]); /* d_yield_dc                               */

EXTERN  int check_for_neg_elem_volume
( int,                            /* element block to check            */
        double [],                      /* x                                 */
        double [],			/* resid_vector                      */
        double [],			/* x_old                             */
        double [],			/* x_older                           */
        double [],			/* xdot                              */
        double [],			/* xdot_old                          */
        double [],			/* x_update                          */
        double *,			/* delta_t                           */
        double *,			/* theta                             */
        double *,			/* time_value                        */
        Exo_DB * );			/* exo                               */

EXTERN void friction_n_dot_f_bc /* mm_fill_solid.c                         */
(double [DIM],		/* func                                      */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
       const double ,			/* friction coefficient             */
       const int ,			/*  block id */
       const double ,			/* delta_t         */
       const double ,			/* theta         */
       const int ,			/* ip         */
       const int ,			/* ip_total         */
       const int ,			/* BC_Name     */
       const double ,			/* time_value  */
       const double[] ,			/* parameter_list*/
       const int );			/* number of parameters  */
 
#endif /* GOMA_MM_FILL_SOLID_H */
