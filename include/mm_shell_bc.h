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
 

#ifndef _MM_SHELL_BC_H
#define _MM_SHELL_BC_H

#ifdef EXTERN
#undef EXTERN
#endif

#ifdef _MM_SHELL_BC_C
#define EXTERN /* do nothing */
#endif

#ifndef _MM_SHELL_BC_C
#define EXTERN extern
#endif

EXTERN void shell_n_dot_flow_bc_confined
PROTO((double func[DIM],
       double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
       const double flowrate, /* imposed flow rate */
       const double time,     /* current time */
       const double dt,       /* current time step size */
       double xi[DIM],        /* Local stu coords */
       const Exo_DB *exo));   /* ExodusII database struct pointer */
       

EXTERN void shell_n_dot_flow_bc_film
PROTO((double func[DIM],
       double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
       const double flowrate, /* imposed flow rate */
       const double time,     /* current time */
       const double dt,       /* current time step size */
       double xi[DIM],        /* Local stu coords */
       const Exo_DB *exo));   /* ExodusII database struct pointer */


EXTERN void shell_n_dot_gradp_bc
PROTO((double func[DIM],
       double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
       const double time,     /* current time */
       const double dt,       /* current time step size */
       double xi[DIM],        /* Local stu coords */
       const Exo_DB *exo));   /* ExodusII database struct pointer */


EXTERN void shell_n_dot_gradh_bc
PROTO((double [DIM],            /* func */
       double [DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE], /* d_func */
       const double ));          

EXTERN void shell_n_dot_pflux_bc
PROTO((double func[DIM],
       double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
       const double flux,     /* imposed particles flux */
       const double time,     /* current time */
       const double dt));     /* current time step size */       

EXTERN void match_lubrication_film_pressure
PROTO((double func[DIM],
       double d_func[DIM][MAX_VARIABLE_TYPES + MAX_CONC][MDE],
       const int eb_mat_lubp,     /* element block for lubp equation */
       const int eb_mat_filmp));  /* element block for filmp equation */


EXTERN void put_lub_flux_in_film /* mm_shell_bc.c                    */
PROTO((int ,			/* id - local element node number for the 
				 * current node whose residual contribution
				 * is being sought                           */
       int ,			/* I - Global node number                    */
       int ,			/* ielem_dim - physical dimension of element,
				 * ie., 1, 2, 3                              */
       double [],		/* resid_vector - Residual vector NO DUH!    */
       int ,			/* i_mat_lubp     */
       int ,			/* i_mat_filmp    */
       int []));		/* local_node_list_fs - MDE list to keep track
				 * of nodes at which liquid contributions have
				 * been transfered to solid (fluid-solid 
				 * boundaries)                               */

#endif
