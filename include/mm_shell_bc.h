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
 

#ifndef GOMA_MM_SHELL_BC_H
#define GOMA_MM_SHELL_BC_H

#include "el_elm.h"
#include "exo_struct.h"
#include "mm_eh.h"
#include "rf_fem_const.h"
#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_MM_SHELL_BC_C
#define EXTERN /* do nothing */
#endif

#ifndef GOMA_MM_SHELL_BC_C
#define EXTERN extern
#endif

EXTERN void shell_n_dot_flow_bc_confined
(double func[DIM],
       double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
       const double flowrate, /* imposed flow rate */
       const double time,     /* current time */
       const double dt,       /* current time step size */
       double xi[DIM],        /* Local stu coords */
       const Exo_DB *exo);   /* ExodusII database struct pointer */
       
EXTERN void lub_static_pressure
(double func[DIM],
       double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
       const double P_atm,    /* imposed atmospheric pressure */
       const double time,     /* current time */
       const double dt,       /* current time step size */
       double xi[DIM],        /* Local stu coords */
       const Exo_DB *exo);   /* ExodusII database struct pointer */

EXTERN void shell_n_dot_flow_bc_film
(double func[DIM],
       double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
       const double flowrate, /* imposed flow rate */
       const double time,     /* current time */
       const double dt,       /* current time step size */
       double xi[DIM],        /* Local stu coords */
       const Exo_DB *exo);   /* ExodusII database struct pointer */


EXTERN void shell_n_dot_gradp_bc
(double func[DIM],
       double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
       const double time,     /* current time */
       const double dt,       /* current time step size */
       double xi[DIM],        /* Local stu coords */
       const Exo_DB *exo);   /* ExodusII database struct pointer */


EXTERN void shell_n_dot_gradh_bc
(double [DIM],            /* func */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func */
       const double );          

EXTERN void shell_n_dot_pflux_bc
(double func[DIM],
       double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
       const double flux,     /* imposed particles flux */
       const double time,     /* current time */
       const double dt);     /* current time step size */       

EXTERN void shell_n_dot_gas_velo_bc_tfmp
(double func[DIM],
       double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
       const double flux,     /* imposed particles flux */
       const double time,     /* current time */
       const double dt,       /* current time step size */
       double xi[DIM],        /* Local stu coords */
       const Exo_DB *exo);   /* ExodusII database struct pointer */

EXTERN void shell_lubrication_outflow
(double func[DIM],
       double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
       const double time,     /* current time */
       const double dt,       /* current time step size */
       double xi[DIM],        /* Local stu coords */
       const Exo_DB *exo);   /* ExodusII database struct pointer */

EXTERN void apply_sdet
(double func[DIM],
       double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
       double xi[DIM],        /* Local stu coords */
       const Exo_DB *exo);   /* ExodusII database struct pointer */

EXTERN void apply_sh_weak
(double func[DIM],
       double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
       double xi[DIM],        /* Local stu coords */
       const Exo_DB *exo,   /* ExodusII database struct pointer */
       const double dy_ds);  /* applied dy_ds at boundary */

EXTERN void shell_n_dot_liq_velo_bc_tfmp
(double func[DIM],
       double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
       const double flux,     /* imposed particles flux */
       const double time,     /* current time */
       const double dt,       /* current time step size */
       double xi[DIM],        /* Local stu coords */
       const Exo_DB *exo);   /* ExodusII database struct pointer */

EXTERN void shell_num_diff_bc_tfmp
(double func[DIM],
       double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
       const double time,     /* current time */
       const double dt,       /* current time step size */
       double xi[DIM],        /* Local stu coords */
       const Exo_DB *exo);   /* ExodusII database struct pointer */

EXTERN void shell_tfmp_avg_plate_velo_liq
(double func[DIM],
       double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
       const double time,     /* current time */
       const double dt,       /* current time step size */
       double xi[DIM],        /* Local stu coords */
       const Exo_DB *exo);   /* ExodusII database struct pointer */

EXTERN void shell_tfmp_n_dot_grad_s
(double func[DIM],
       double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
       const double time,     /* current time */
       const double dt,       /* current time step size */
       double xi[DIM],        /* Local stu coords */
       const Exo_DB *exo);   /* ExodusII database struct pointer */

EXTERN void apply_shell_traction_bc
(double func[DIM],
       double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
       const int BC_name,    /* BC name identifier */
       const double tx,      /* traction in x-direction */
       const double ty,      /* traction in y-direction */
       const double tz);    /* traction in z-direction */

EXTERN void match_lubrication_film_pressure
(double func[DIM],
       double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
       const int eb_mat_lubp,     /* element block for lubp equation */
       const int eb_mat_filmp);  /* element block for filmp equation */


EXTERN void put_lub_flux_in_film /* mm_shell_bc.c                    */
(int ,			/* id - local element node number for the 
				 * current node whose residual contribution
				 * is being sought                           */
       int ,			/* I - Global node number                    */
       int ,			/* ielem_dim - physical dimension of element,
				 * ie., 1, 2, 3                              */
       double [],		/* resid_vector - Residual vector NO DUH!    */
       int ,			/* i_mat_lubp     */
       int ,			/* i_mat_filmp    */
       int []);		/* local_node_list_fs - MDE list to keep track
				 * of nodes at which liquid contributions have
				 * been transfered to solid (fluid-solid 
				 * boundaries)                               */
EXTERN void 
shell_tfmp_avg_plate_velo_gas(double func[DIM],
			     double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
			     const double time,     /* current time */
			     const double delta_t,       /* current time step size */
			     double xi[DIM],        /* Local stu coordinates */
			     const Exo_DB *exo); 

#endif
