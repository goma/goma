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
 
#ifndef GOMA_BC_SPECIAL_H
#define GOMA_BC_SPECIAL_H

#include "el_elm.h"
#include "exo_struct.h"
#include "rf_bc_const.h"
#include "sl_util_structs.h"

#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_BC_SPECIAL_C
#define EXTERN /* do nothing */
#endif

#ifndef GOMA_BC_SPECIAL_C
#define EXTERN extern
#endif

EXTERN int apply_special_bc
(struct Aztec_Linear_Solver_System *,
       double [],		/* x - soln vector                           */
       double [],		/* resid_vector - RHS                        */
       double [],		/* x_old - soln vector previous timestep     */
       double [],		/* x_older - soln vector preprevious timestep*/
       double [],		/* xdot - dx/dt                              */
       double [],		/* xdot_old - dxdt at previous timestep      */
       double ,			/* delta_t - current time step size          */
       double ,			/* theta - parameter to vary time 
				 * integration algorithm from
				 * explicit (theta = 1) to
				 * -- implicit (theta = 0)                   */
       struct elem_side_bc_struct *[], /* first_elem_side_BC_array - array of
					* pointers to the first surface 
					* integral defined for each element. 
					* It has a length equal to the total 
					* number of elements defined on the 
					* current processor                  */
       int ,			/* ielem - element number                    */
       int ,			/* ip_total - number of gauss points         */
       int ,			/* ielem_type - element type flag            */
       int ,			/* num_local_nodes -                         */
       int ,			/* ielem_dim -                               */
       int ,			/* iconnect_ptr                              */
       struct elem_side_bc_struct *, /* elem_side_bc - Pointer to an element
				      *  side boundary condition structure   */
       int ,			/* num_total_nodes                           */
       int ,			/* bc_application                            */
       int [],                  /* CA_id CA condition id array */
       int [],                   /* CA_fselem free surface element array */
       int [],                   /* CA_sselem solid surface element array */
       Exo_DB *,		/* exo - ptr to FE EXODUS II database        */
       double );		/* time_value */

EXTERN int apply_shell_grad_bc
(double [],               /* x - Soln vector                           */
       double [],               /* resid_vector -                            */
       const double ,           /* delta_t - current time step size          */
       const double ,           /* theta - parameter (0 to 1) to vary time 
                                 * integration (implicit=0, explicit=1)      */
       const double ,           /* h_elem_avg - global average element size  */
       const double [DIM],      /* h - average element size                  */
       const double ,           /* mu_avg - average element viscosity        */
       const double ,           /* U_norm - global velocity norm             */
       const int ,              /* ielem - element number                    */
       const int ,              /* ielem_type - element type                 */
       const int ,              /* num_local_nodes -                         */
       const int ,              /* ielem_dim -                               */
       const int ,              /* iconnect_ptr                              */
       struct elem_side_bc_struct *, /* elem_side_bc - Pointer to an element 
                                      * side boundary condition structure    */
       const int ,              /* num_total_nodes                           */
       const int ,              /* bc_application - flag indicating whether 
                                 * to integrate strong or weak BC's          */
       const double ,           /* time_value                                */
       const Exo_DB *);        /* exo - ptr to FE database                  */

EXTERN int apply_sharp_integrated_bc 
(double [],           /* Solution vector for the current processor */
	double [],     /* Residual vector for the current processor */
	const double , /* current time */
	const double , /* current time step size */
	const double ,   /* parameter (0 to 1) to vary time integration
						 *  ( implicit - 0 to explicit - 1)          */
	const double[DIM], /* hsquared */
        const int ,       /* element number */
	const int ,  /* element type */
	const int ,
	const int ,
	const int ,
	ELEM_SIDE_BC_STRUCT *,
	                    /* Pointer to an element side boundary condition                                 * structure */
	const int, 
	const Exo_DB *);
        
EXTERN void assemble_sharp_integrated_bc
(double [],           /* Solution vector for the current processor */
	double [],/* Residual vector for the current processor */
	const double , /* current time  */
	const double , /* current time step size */
	const double ,   /* parameter (0 to 1) to vary time integration
						 *  ( implicit - 0 to explicit - 1)          */
        const int ,       /* element number */
	const int ,  /* element type */
	const int ,
	const int ,
	const int ,
	ELEM_SIDE_BC_STRUCT *,
	                    /* Pointer to an element side boundary condition                                 * structure */
	const int, 
	const Exo_DB *,
        double [DIM],
        double,
        int);

#endif /* GOMA_BC_SPECIAL_H */
