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
 

#ifndef _USER_BC_H
#define _USER_BC_H

EXTERN dbl velo_vary_fnc
PROTO((const int ,		/* velo_condition                            */
       const dbl ,		/* x1                                        */
       const dbl ,		/* x2                                        */
       const dbl ,		/* x3                                        */
       const dbl [],		/* p                                         */
       const dbl ));		/* time                                      */

EXTERN dbl dvelo_vary_fnc_d1
PROTO((const int ,		/* velo_condition                            */
       const dbl ,		/* x1                                        */
       const dbl ,		/* x2                                        */
       const dbl ,		/* x3                                        */
       const dbl [],		/* p                                         */
       const dbl ));		/* time                                      */

EXTERN dbl dvelo_vary_fnc_d2
PROTO((const int ,		/* velo_condition                            */
       const dbl ,		/* x1                                        */
       const dbl ,		/* x2                                        */
       const dbl ,		/* x3                                        */
       const dbl [],		/* p                                         */
       const dbl ));		/* time                                      */

EXTERN dbl dvelo_vary_fnc_d3
PROTO((const int ,		/* velo_condition                            */
       const dbl ,		/* x1                                        */
       const dbl ,		/* x2                                        */
       const dbl ,		/* x3                                        */
       const dbl [],		/* p                                         */
       const dbl ));		/* time                                      */

EXTERN dbl fnc
PROTO((const dbl ,		/* x1                                        */
       const dbl ,		/* x2                                        */
       const dbl ,		/* x3                                        */
       const dbl [],		/* p                                         */
       const dbl ));		/* time                                      */

EXTERN dbl dfncd1
PROTO((const dbl ,		/* x1                                        */
       const dbl ,		/* x2                                        */
       const dbl ,		/* x3                                        */
       const dbl [],		/* p                                         */
       const dbl ));		/* time                                      */

EXTERN dbl dfncd2
PROTO((const dbl ,		/* x1                                        */
       const dbl ,		/* x2                                        */
       const dbl ,		/* x3                                        */
       const dbl [],		/* p                                         */
       const dbl ));		/* time                                      */

EXTERN dbl dfncd3
PROTO((const dbl ,		/* x1                                        */
       const dbl ,		/* x2                                        */
       const dbl ,		/* x3                                        */
       const dbl [],		/* p                                         */
       const dbl ));		/* time                                      */

EXTERN void quser_surf
PROTO((double [DIM],		/* func                                      */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
       double [],		/* p - parameterize heat transfer model      */
       const dbl ));		/* time                                      */

EXTERN void tuser
PROTO((double *,		/* func                                      */
       double [],		/* d_func - [MAX_VARIABLE_TYPES + MAX_CONC]  */
       const double [],		/* u_bc - parameterize temperature eqn model */
       const double ));		/* time */

EXTERN void yuser_surf
PROTO((double *,		/* func                                      */
       double [DIM][MAX_VARIABLE_TYPES+MAX_CONC][MDE], /* d_func             */
       const int ,		/* species                                   */
       const double [],		/* u_bc - to parameterize species eqn model  */
       const double ));		/* time */

EXTERN void y2_electroneutrality_surf
PROTO((double *,		/* func                                      */
       double [DIM][MAX_VARIABLE_TYPES+MAX_CONC][MDE], /* d_func             */
       const int ,		/* species                                   */
       const double [],		/* u_bc - to parameterize species eqn model  */
       const double ));		/* time */

EXTERN void uuser_surf
PROTO((double [DIM],		/* func                                      */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
       double [],		/* u_bc - parameterize u velocity model      */
       const dbl ));		/* time                                      */

EXTERN void vuser_surf
PROTO((double [DIM],		/* func                                      */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
       double [],		/* u_bc - parameterize u velocity model      */
       const dbl ));		/* time                                      */

EXTERN void wuser_surf
PROTO((double [DIM],		/* func                                      */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
       double [],		/* u_bc - parameterize u velocity model      */
       const dbl ));		/* time  */

EXTERN void uuser_colloc_surf
PROTO((double *,                /* func                                      */
       double [],               /* d_func           */
       const double [],         /* u_bc - parameterize u velocity model      */
       const int ,              /* Node ID */
       const dbl ));            /* time  */

EXTERN void vuser_colloc_surf
PROTO((double *,                /* func                                      */
       double [],               /* d_func           */
       const double [],         /* u_bc - parameterize v velocity model      */
       const int ,              /* Node ID */
       const dbl ));            /* time  */

EXTERN void wuser_colloc_surf
PROTO((double *,                /* func                                      */
       double [],               /* d_func           */
       const double [],         /* u_bc - parameterize w velocity model      */
       const int ,              /* Node ID */
       const dbl ));            /* time  */

EXTERN void dx_user_surf
PROTO((double *,		/* func                                      */
       double [],               /* d_func           */
       const double [],		/* u_bc - parameterize u velocity model      */
       const dbl ));		/* time  */       
EXTERN void dy_user_surf
PROTO((double *,		/* func                                      */
       double [],               /* d_func           */
       const double [],		/* u_bc - parameterize u velocity model      */
       const dbl ));		/* time  */       
EXTERN void dz_user_surf
PROTO((double *,	          	/* func                                      */
       double [],               /* d_func           */
       const double [],		/* u_bc - parameterize u velocity model      */
       const dbl ));		/* time  */       

EXTERN void p_liq_user_surf
PROTO((double *,	          	/* func                                      */
       double [],               /* d_func           */
       const double [],		/* u_bc - parameterize u velocity model      */
       const dbl ));		/* time  */       

EXTERN void shell_p_open_user_surf 
PROTO((double *func, 
       double d_func[],
       const double u_bc[],
       const double time));

EXTERN void fn_dot_T_user
PROTO((double [DIM],		/* func                                      */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
       const double [],		/* u_bc                                      */
       const dbl ));		/* time                                      */

EXTERN void flow_n_dot_T_user
PROTO((double [DIM],		/* func                                      */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
       const double [],		/* u_BC - Parameters from input deck         */
       const dbl ));		/* time                                      */

EXTERN double var_CA_user
PROTO((double ,			/* Ca_local */
       int ,			/* num */
       const double *,		/* a - parameter list from user */
       double *));		/* d_cos_CA_Ca_local */

EXTERN int user_gibbs_criterion
PROTO((const double [MAX_PDIM], /* fsnormal - Vector of free surface normal
				 * components                                */
       const double [MAX_PDIM], /* ssnormal - Vector of solid surface normal 
				 * components                                */
       const int ,		/* imodel - Flag which tracks which model    */
       int * ,			/* ipin  - Flag for pinned or not            */
       const double []));	/* p - User defined parameter list, or model
				 * spec. list                                */

EXTERN void force_user_surf	/* user_bc.c                                 */
PROTO((double [DIM],		/* func                                      */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
       double [],		/* p - user parameter list                   */
       dbl ));			/* time                                      */

EXTERN void volt_user_surf	/* user_bc.c                                 */
PROTO((double [DIM],		/* func                                      */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
       double [],		/* p - user parameter list                   */
       const dbl ));		/* time                                      */

EXTERN void current_user_surf	/* user_bc.c                                 */
PROTO((double [DIM],		/* func                                      */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func           */
       double [],		/* p - user parameter list                   */
       const dbl ));		/* time                                      */

extern void mass_flux_user_surf
PROTO((double [MAX_CONC],
       double [MAX_CONC][MAX_VARIABLE_TYPES+MAX_CONC],
       const int, const double [], const double));

#endif /* _USER_BC_H */
