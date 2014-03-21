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
 * mm_fill_fill.h -- prototype declarations for mm_fill_fill.c
 */


/*
 *$Id: mm_fill_fill.h,v 5.3 2008-11-04 20:12:05 hkmoffa Exp $
 */

#ifndef _MM_FILL_FILL_H
#define _MM_FILL_FILL_H

#ifdef EXTERN
#undef EXTERN
#endif

#ifdef _MM_FILL_FILL_C
#define EXTERN /* do nothing */
#endif

#ifndef _MM_FILL_FILL_C
#define EXTERN extern
#endif

EXTERN int integrate_explicit_eqn
PROTO((struct Aztec_Linear_Solver_System *, /* ams - cf "sl_util_structs.h"  */
       double [],		/* rf - residual for fill equation only      */
       double [],		/* xf - vector with fill at nodes only       */
       double [],		/* xf_old - vector with fill at nodes only   */
       double [],		/* xfdot - vector with fill at nodes only    */
       double [],		/* xfdot_old - vector with fill at nodes only*/
       double [],		/* x - Solution vector for                   */
       double [],		/* x_old - Solution vector at last timestep  */
       double [],		/* xdot - time derivative of soln vector     */
       dbl ,			/* delta_t - time step size                  */
       dbl ,   			/* theta - time integration parameter        */
       dbl *,			/* time2 -                                   */
       int  ,                   /* equation type to integrate */
       int [],			/* node_to_fill - this is a map from the 
				 * node number to the unknown number for the
				 * fill equation                             */
       Exo_DB *,		/* exo - ptr to exodus file                  */
       Dpi *,  			/* dpi - distributed processing info         */
       Comm_Ex * ));            /* cx - comm structure                       */


#ifndef COUPLED_FILL
EXTERN int assemble_fill
PROTO((double [],               /* afill - Jacobian matrix for fill equation */
       int [],                  /* ijaf - pointer to nonzeros in Jacobian    */
       double [],               /* rf - rhs vector                           */
       double ,                 /* dt - current time step size               */
       double ,                 /* tt - parameter to vary time integration from
                                 * explicit (tt = 1) to implicit (tt = 0) */
       int [] ,               /* node_to_fill -  */
       const int ));               /* equation -  */

EXTERN int assemble_fill_ext_v
PROTO((double [],               /* afill - Jacobian matrix for fill equation */
       int [],                  /* ijaf - pointer to nonzeros in Jacobian    */
       double [],               /* rf - rhs vector                           */
       double ,                 /* dt - current time step size               */
       double ,                 /* tt - parameter to vary time integration from
                                 * explicit (tt = 1) to implicit (tt = 0) */
       int [],                  /* node_to_fill -  */
       double []));             /* element size measure for GLS */  

EXTERN int assemble_fill_fake
PROTO((double ,			/* tt - parameter varies time integration from 
				 *      explicit (tt = 1) to 
				 *      implicit (tt = 0)                    */
       double ));               /* dt - current time step size               */

EXTERN int assemble_fill_gradf
PROTO((double [],               /* afill[], Jacobian matrix for fill equation      */
       int [],                  /* ijaf[],  pointer to nonzeros in Jacobian matrix */
       double [],               /* rf[],    rhs vector                             */
       double ,                 /* dt, current time step size                      */
       double ,                 /* tt, parameter to vary time integration from     *
			         * explicit (tt = 1) to implicit (tt = 0)          */
       int [],                  /* node_to_fill[],                                 */   
       double []));             /* hsquared[DIM]), element size                    */

#endif /* ifndef COUPLED_FILL */

EXTERN int assemble_fill_fake
PROTO((double ,			/* tt - parameter varies time integration from 
				 *      explicit (tt = 1) to 
				 *      implicit (tt = 0)                    */
       double ));               /* dt - current time step size               */



#ifdef COUPLED_FILL
EXTERN int assemble_fill
PROTO((double ,                      /* tt                   */
       double ,                      /* dt                   */
       double [],                    /* hsquared[DIM]        */
       double [][DIM],               /* hh[DIM][DIM]         */
       double [][MDE],               /* dh_dxnode[DIM][MDE]  */
       const int,                    /* equation -  */
       double [DIM],                 /* xi - coordinates */
       Exo_DB * const,               /* Exodus database */
       double                        /* Time */
       ));

EXTERN int assemble_fill_ext_v
PROTO((double ,			/* tt - parameter varies time integration from 
				 *      explicit (tt = 1) to 
				 *      implicit (tt = 0)                    */
       double,                  /* dt - current time step size               */
       double [],               /* hsquared[DIM]           next three parameters */
       double [][DIM],          /* hh[DIM][DIM]            are element size  */
       double [][MDE]));        /* dh_dxnode[DIM][MDE]     info for GLS */

EXTERN int assemble_fill_gradf
PROTO((double ,                      /* tt                   */
       double ,                      /* dt                   */
       double [],                    /* hsquared[DIM]        */
       double [][DIM],               /* hh[DIM][DIM]         */
       double [][MDE]));             /* dh_dxnode[DIM][MDE]  */

#endif /* COUPLED_FILL */

EXTERN void map_direction
PROTO((double *, 
       double * ));

EXTERN int fill_matrix
PROTO((double [],		/* afill - matrix for fill variables only */
       int [],			/* ijaf - pointer to nonzeros in fill matrix */
       double [],		/* rf - residual for fill equation only */
       double [],		/* xf - vector with fill at nodes only     */
       double [],		/* x - Solution vector for the current proc */
       double [],		/* x_old - Solution vector at last time step */
       double [],		/* xdot - time derivative of solution vector */
       dbl,			/* delta_t - time step size */
       dbl,			/* theta - time integration parameter */
       int ,                    /* equation type being integrated */
       int [],			/* node_to_fill - mapping */
       struct elem_side_bc_struct *[],
       Exo_DB *,	        /* exo - ptr to finite element database */
       Dpi *));			/* distributed processing  info */



EXTERN int assemble_surface
PROTO((Exo_DB *,		/* exo - ptr to basic exodus ii mesh info */
       double [],		/* x - global vector containing all unknowns */
       double [],		/* afill - Jacobian matrix fill equation */
       int [],			/* ijaf - ptr to nonzeros Jacobian matrix */
       double [],		/* rf - rhs vector   */
       double ,			/* delta_t - current time step size */
       double ,			/* theta - parameter to vary time integration 
				 * explicit (theta = 1) implicit (theta = 0) */
       int [],			/* node_to_fill  */
       int ,			/* ielem_type - element type  */
       int ,			/* ielem_type_fill - element type for 
				 * fill function */
       int ,			/* id_side - id number of current side 
				 * according to EXODUS convention  */
       double,			/* F_inlet */
       int ,			/* neighbor - element neighboring this side */
       int ,			/* ielem - current element */
       int ));			/* num_local_nodes - number of nodes per 
				 * element */
EXTERN int get_side_info
PROTO((const int ,              /* ielem_type                                */
       const int ,              /* id_side - 1-based EXODUS/PATRAN side num  */
       int *,			/* nodes_per_side                            */
       int []));		/* local_elem_node_id                        */

EXTERN int assemble_surface_species
PROTO((Exo_DB *,		/* exo - ptr to basic exodus ii mesh info    */
       double [],		/* x                                         */
       double ,			/* delta_t - current time step size          */
       double ,			/* theta - parameter to vary time integration 
				 * from explicit (tt = 1) 
				 * to implicit (tt = 0)                      */
       int ,			/* ielem_type - element type                 */
       int ,			/* ielem_type_fill - element type fill func  */
       int ,			/* id_side - id number of current side 
				 * according to EXODUS convention            */
       int ,			/* neighbor - element neighboring this side  */
       int ,			/* ielem - current element                   */
       int ));			/* num_local_nodes - number nodes per elem   */

EXTERN int elem_on_ss
PROTO((Exo_DB *,		/* exo                                       */
       int,			/* ss_id                                     */
       int ));			/* elem                                      */


#ifdef LS_CODE

EXTERN void correct_level_set
PROTO((struct Aztec_Linear_Solver_System *ams,
       double [],
       double [],
       double [],
       double [],
       double [],
       int    [],
       int      ,
       int      ,
       double   ,
       double   ,
       int      ,
       int      ,
       Exo_DB  *,
       Dpi     *,
       Comm_Ex * ));


EXTERN int apply_ls_inlet_bc 
PROTO((double [],		/* afill - Jacobian matrix fill equation */
       int [],			/* ijaf - ptr to nonzeros Jacobian matrix */
       double [],		/* x - global vector containing all unknowns */
       double [],		/* rf - rhs vector   */
       int [],
       struct elem_side_bc_struct *,
       Exo_DB * ));


EXTERN int apply_strong_fill_ca_bc 
PROTO((double [],		/* afill - Jacobian matrix fill equation */
       int [],			/* ijaf - ptr to nonzeros Jacobian matrix */
       double [],		/* x - global vector containing all unknowns */
       double [],		/* rf - rhs vector   */
       const double ,		/* delta_t - current time step size */
       const double ,		/* theta - parameter to vary time integration 
				 * explicit (theta = 1) implicit (theta = 0) */
       int [],			/* node_to_fill  */
       const int ,		/* ielem - current element */
       int ,			/* ielem_type - element type  */
       const int ,		/* num_local_nodes - number of nodes per element */
       const int ,		/* ielem_dim */
       const int ,		/* iconnect_ptr */
       struct elem_side_bc_struct *,
       const int ,		/* num_total_nodes */
       const double,			/* contact angle */
       const Exo_DB *));	/* exo - ptr to basic exodus ii mesh info */

EXTERN int elem_on_isosurface
PROTO(( int,
        double [],
        const Exo_DB *,
        int,
        double ));

EXTERN int current_elem_on_isosurface
PROTO(( int,
        double ));

EXTERN int huygens_renormalization
PROTO(( double *,
        const int ,
        Exo_DB *,
	Comm_Ex *,
        Dpi    *,
	int     ,
	int     ,
        double  ,
	int    ));

EXTERN void cgm_based_initialization 
PROTO (( double *,
         int));

EXTERN void surf_based_initialization 
PROTO (( double *,
         double *,
         double *,
	 Exo_DB  *,
         int,
	 struct LS_Surf_List *,
         double,
         double,
         double ));

EXTERN void free_surf_list
PROTO (( struct LS_Surf_List **));

EXTERN struct LS_Surf * closest_surf
PROTO (( struct LS_Surf_List *,
         double *,
         Exo_DB  *,
         double * ));

EXTERN void find_surf_closest_point 
PROTO (( struct LS_Surf *,
         double *,
         Exo_DB  *,
         double * ));

EXTERN void find_facets
PROTO((  struct LS_Surf_List *,
	 int,
         double,
	 Exo_DB * ));

EXTERN struct LS_Surf * create_surf
PROTO (( int ));

EXTERN struct LS_Surf * create_surf_point
PROTO (( double *,
         int,
         double *,
         int ));

EXTERN struct LS_Surf * create_surf_facet_line
PROTO (( struct LS_Surf *,
         struct LS_Surf * ));
         
EXTERN void append_surf
PROTO (( struct LS_Surf_List *,
         struct LS_Surf * ));

EXTERN void create_subsurfs
PROTO (( struct LS_Surf_List *,
         double *,
         Exo_DB * ));

EXTERN void free_subsurfs
PROTO (( struct LS_Surf_List * ));

EXTERN void assemble_Global_surf_list
PROTO(( struct LS_Surf_List * ));

                     
EXTERN int level_set_interface
PROTO(( const double,          /*  F                         */
	const double [DIM],    /*  grad_F[DIM]               */
	const double,          /*  width                     */
	const int,             /*  do_deriv                  */
	int *,                 /* *near                      */
	double *,              /* *H                         */
	double *,              /* *d_H_dF                    */
	double [DIM],          /*  d_H_dgradF[DIM]           */
	double *,              /* *delta                     */
	double *,              /* *d_delta_dF                */
	double [DIM],          /*  d_delta_dgradF[DIM]       */
	double [DIM],          /*  normal[DIM]               */
	double [DIM],          /*  d_normal_dF[DIM]          */
	double [DIM][DIM]));   /*  d_normal_dgradF[DIM][DIM] */
	
EXTERN int level_set_property
PROTO(( const double,          /*  p0           */
	const double,          /*  p1           */
	const double,          /*  width        */
	double *,              /* *pp           */
	double [MDE]));        /*  d_pp_dF[MDE] */
                     
EXTERN int ls_transport_property
PROTO(( const double,          /*  p0           */
	const double,          /*  p1           */
	const double,          /*  width        */
	double *,              /* *pp           */
	double * ));           /*  d_pp_dF */

EXTERN void zero_lsi
PROTO(( void));

EXTERN int load_lsi
PROTO(( const double ));       /* width */

EXTERN int load_lsi_shell_second
PROTO(( const double ));       /* width */

EXTERN int load_lsi_derivs
PROTO(( void));


EXTERN double find_LS_global_flux
PROTO(( const Exo_DB *,
	const Dpi *,
	const double *,
	      double *,
	double [],
	int ));

EXTERN double find_LS_vel
PROTO(( const Exo_DB *,       /* exo */
	const Dpi *,          /* dpi */
	const double *,       /* params[] */
	const int,            /* chosen_vel */
	double [],            /* x[] (solution vector) */
	int ));               /* num_total_unknowns */

EXTERN void semi_lagrange_step
PROTO(( int ,               /* num_total_nodes */ 
	int ,               /* num_total_unknowns */
	int ,               /* num_fill_unknowns */
	double [],          /* x  */
	double [],          /* xf */
	double [],          /*   xf_old   */
	double [],          /*    xfdot   */
	double [],          /*   xfdot_old   */
	int [],             /*  node_to_fill  */
	double ,            /*  delta_t   */
	double  ,           /*  theta     */
	Exo_DB *,           /*  exo       */
	Dpi    *,           /*  dpi      */
	Comm_Ex *));        /*  cx      */
 
#endif


EXTERN int boundary_curvature
PROTO(( double [DIM],
	double [DIM][MAX_VARIABLE_TYPES+MAX_CONC][MDE],
	double * ));

EXTERN int curvature_momentum_source
PROTO(( double [DIM],         
	double [DIM][DIM][MDE],
	double [DIM][MDE],
	double [DIM][MDE] ));

#ifdef CONST_CURVE
EXTERN int curvature_momentum_source_const
PROTO(( double [DIM],         
	double [DIM][DIM][MDE],
	double [DIM][MDE],
	double [DIM][MDE] ));
#endif


/* Phase function kernel */

EXTERN int assemble_phase_function
PROTO(( double ,         
	double ,
	double ,
	double [DIM],
	Exo_DB *));
	   
	   
EXTERN int assemble_pf_constraint
PROTO(( double, 
		double *,
		double  ,
		double *,
		double * ));


#endif /* _MM_FILL_FILL_H */
