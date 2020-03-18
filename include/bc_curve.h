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
 
#ifndef GOMA_BC_CURVE_H
#define GOMA_BC_CURVE_H


#include "bc_colloc.h"
#include "exo_struct.h"

struct elem_edge_bc_struct;
#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_BC_CURVE_C
#define EXTERN
#
#endif

#ifndef GOMA_BC_CURVE_C
#define EXTERN extern
#endif


EXTERN int apply_integrated_curve_bc
(double [],		/* x - soln vector for the current processor */
       double [],		/* resid_vector - for the current processor  */
       const double ,		/* delta_t - current time step size          */
       const double ,		/* theta - parameter to vary time integration 
				 * from explicit (theta=1) to 
				 * implicit (theta=0)                        */
       const int ,		/* ielem - element number                    */
       const int ,		/* ielem_type - element type                 */
       const int ,		/* num_local_nodes                           */
       const int ,		/* ielem_dim                                 */
       const int ,		/* iconnect_ptr                              */
       struct elem_edge_bc_struct *, /* elem_edge_bc - Pointer to an element 
				      * side boundary condition structure    */
       const int ,		/* num_total_nodes */
       const int ,		/* bc_application - flag indicating whether 
				 * to integrate strong or weak BC's          */
       const Exo_DB *);	/* exo - ptr to FE database                  */

EXTERN int apply_point_colloc_edge_bc
(double [],		/* x - Soln vector                           */
       double [],		/* x_old - Soln vector, previous timestep    */
       double [],		/* x_older - Soln vector, previous2 timestep */
       double [],		/* xdot - current solution                   */
       double [],		/* xdot_old - previous timestep              */
       double [],		/* resid_vector - Residual vector            */
       const double ,		/* delta_t - current time step size          */
       const double ,		/* theta - parameter to vary time integration:
				 * explicit (theta = 1) -- 
				 * implicit (theta = 0)                      */
       const int ,		/* ielem - element number                    */
       const int ,		/* ip_total - total number of gauss points   */
       const int ,		/* ielem_type - element type                 */
       const int ,		/* num_local_nodes                           */
       const int ,		/* ielem_dim                                 */
       const int ,		/* iconnect_ptr                              */
       struct elem_edge_bc_struct *, /* elem_edge_bc - Ptr to an element
					    *  side BC structure             */
       const int ,		/* num_total_nodes                           */
       int [],			/* local_node_list_fs - dimensioned [MDE]; 
				 * list to keep track of nodes at which solid 
				 * contributions have been transfered to 
				 * liquid (fluid-solid boundaries)           */
       const double);		/* time_value                                */

#endif /* GOMA_BC_CURVE_H */
