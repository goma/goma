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
 
#ifndef _BC_INTEG_H
#define _BC_INTEG_H


#ifdef EXTERN
#undef EXTERN
#endif

#ifdef _BC_INTEG_C
#define EXTERN /* do nothing */
#endif

#ifndef _BC_INTEG_C
#define EXTERN extern
#endif

EXTERN int apply_integrated_bc
PROTO((double [],		/* x - Soln vector                           */
       double [],		/* resid_vector -                            */
       const double ,		/* delta_t - current time step size          */
       const double ,		/* theta - parameter (0 to 1) to vary time 
				 * integration (implicit=0, explicit=1)      */
       const PG_DATA *,
       const int ,		/* ielem - element number                    */
       const int ,		/* ielem_type - element type                 */
       const int ,		/* num_local_nodes -                         */
       const int ,		/* ielem_dim -                               */
       const int ,		/* iconnect_ptr                              */
       struct elem_side_bc_struct *, /* elem_side_bc - Pointer to an element 
				      * side boundary condition structure    */
       const int ,		/* num_total_nodes                           */
       const int ,		/* bc_application - flag indicating whether 
				 * to integrate strong or weak BC's          */
       const double ,		/* time_value                                */
	   SGRID *, 
       const Exo_DB *));	/* exo - ptr to FE database                  */

EXTERN void apply_table_wic_bc
PROTO((double [],               /* func                                      */
       double [][MAX_VARIABLE_TYPES+MAX_CONC][MDE], /* d_func                */
       struct Boundary_Condition *,/* BC_Type                             */
       double ));                  /* time_value */


#ifdef STATIC

/*
 * Prototype declarations of static functions in bc_colloc.c...
 */


#endif

#endif /* _BC_INTEG_H */
