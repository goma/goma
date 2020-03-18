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

#ifndef GOMA_MM_FILL_FILL_H
#define GOMA_MM_FILL_FILL_H

#include "dp_types.h"
#include "dpi.h"
#include "el_elm.h"
#include "exo_struct.h"
#include "mm_as_structs.h"
#include "mm_fill_aux.h"
#include "rf_fem_const.h"
#include "std.h"

struct Aztec_Linear_Solver_System;
struct elem_side_bc_struct;
#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_MM_FILL_FILL_C
#define EXTERN /* do nothing */
#endif

#ifndef GOMA_MM_FILL_FILL_C
#define EXTERN extern
#endif

struct LS_Mass_Lumped_Penalty {
  dbl penalty[MDE];
  dbl penalty_old[MDE];
  dbl d_penalty[MDE][MDE];
  dbl inv_grad_F[MDE];
};

EXTERN int integrate_explicit_eqn
(struct Aztec_Linear_Solver_System *, /* ams - cf "sl_util_structs.h"  */
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
       Comm_Ex * );            /* cx - comm structure                       */


#ifndef COUPLED_FILL
int 
assemble_fill(double tt, 
	      double dt, 
	      dbl hsquared[DIM], 
	      dbl hh[DIM][DIM], 
	      dbl dh_dxnode[DIM][MDE],
	      const int applied_eqn,
	      double xi[DIM],
	      Exo_DB *exo,
	      double time
	      );

EXTERN int assemble_fill_ext_v
(double [],               /* afill - Jacobian matrix for fill equation */
       int [],                  /* ijaf - pointer to nonzeros in Jacobian    */
       double [],               /* rf - rhs vector                           */
       double ,                 /* dt - current time step size               */
       double ,                 /* tt - parameter to vary time integration from
                                 * explicit (tt = 1) to implicit (tt = 0) */
       int [],                  /* node_to_fill -  */
       double []);             /* element size measure for GLS */  

EXTERN int assemble_fill_fake
(double ,			/* tt - parameter varies time integration from 
				 *      explicit (tt = 1) to 
				 *      implicit (tt = 0)                    */
       double );               /* dt - current time step size               */

EXTERN int assemble_fill_gradf
(double [],               /* afill[], Jacobian matrix for fill equation      */
       int [],                  /* ijaf[],  pointer to nonzeros in Jacobian matrix */
       double [],               /* rf[],    rhs vector                             */
       double ,                 /* dt, current time step size                      */
       double ,                 /* tt, parameter to vary time integration from     *
			         * explicit (tt = 1) to implicit (tt = 0)          */
       int [],                  /* node_to_fill[],                                 */   
       double []);             /* hsquared[DIM]), element size                    */

#endif /* ifndef COUPLED_FILL */

EXTERN int assemble_fill_fake
(double ,			/* tt - parameter varies time integration from 
				 *      explicit (tt = 1) to 
				 *      implicit (tt = 0)                    */
       double );               /* dt - current time step size               */



#ifdef COUPLED_FILL
int 
assemble_fill(double tt, 
	      double dt, 
	      dbl hsquared[DIM], 
	      dbl hh[DIM][DIM], 
	      dbl dh_dxnode[DIM][MDE],
	      const int applied_eqn,
	      double xi[DIM],
	      Exo_DB *exo,
	      double time
	      );

EXTERN int assemble_fill_ext_v
(double ,			/* tt - parameter varies time integration from 
				 *      explicit (tt = 1) to 
				 *      implicit (tt = 0)                    */
       double,                  /* dt - current time step size               */
       double [],               /* hsquared[DIM]           next three parameters */
       double [][DIM],          /* hh[DIM][DIM]            are element size  */
       double [][MDE]);        /* dh_dxnode[DIM][MDE]     info for GLS */

EXTERN int assemble_fill_gradf
(double ,                      /* tt                   */
       double ,                      /* dt                   */
       double [],                    /* hsquared[DIM]        */
       double [][DIM],               /* hh[DIM][DIM]         */
       double [][MDE]);             /* dh_dxnode[DIM][MDE]  */

#endif /* COUPLED_FILL */

EXTERN void map_direction
(double *, 
       double * );

EXTERN int fill_matrix
(double [],		/* afill - matrix for fill variables only */
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
       Dpi *);			/* distributed processing  info */



EXTERN int assemble_surface
(Exo_DB *,		/* exo - ptr to basic exodus ii mesh info */
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
       int );			/* num_local_nodes - number of nodes per 
				 * element */
EXTERN int get_side_info
(const int ,              /* ielem_type                                */
       const int ,              /* id_side - 1-based EXODUS/PATRAN side num  */
       int *,			/* nodes_per_side                            */
       int []);		/* local_elem_node_id                        */

int
assemble_surface_species (Exo_DB *exo,	/* ptr to basic exodus ii mesh information */
			  double x[],
			  double delta_t, /* current time step size */
			  double theta,	/* parameter to vary time integration from
					 * explicit (theta = 1) to implicit (theta = 0) */
			  int ielem_type, /* element type  */
			  int ielem_type_fill, /* element type for fill function */
			  int id_side,	/* id number of current side according to 
					 * EXODUS convention  */
			  int neighbor,	/* element neighboring this side */
			  int ielem,	/* current element */
			  int num_local_nodes);   /* number of nodes per element */

EXTERN int elem_on_ss
(Exo_DB *,		/* exo                                       */
       int,			/* ss_id                                     */
       int );			/* elem                                      */


#ifdef LS_CODE

EXTERN void correct_level_set
(struct Aztec_Linear_Solver_System *ams,
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
       Comm_Ex * );


EXTERN int apply_ls_inlet_bc 
(double [],		/* afill - Jacobian matrix fill equation */
       int [],			/* ijaf - ptr to nonzeros Jacobian matrix */
       double [],		/* x - global vector containing all unknowns */
       double [],		/* rf - rhs vector   */
       int [],
       struct elem_side_bc_struct *,
       Exo_DB * );


EXTERN int apply_strong_fill_ca_bc 
(double [],		/* afill - Jacobian matrix fill equation */
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
       const Exo_DB *);	/* exo - ptr to basic exodus ii mesh info */

EXTERN int elem_on_isosurface
( int,
        double [],
        const Exo_DB *,
        int,
        double );

EXTERN int current_elem_on_isosurface
( int,
        double );

EXTERN int huygens_renormalization
( double *,
        const int ,
        Exo_DB *,
	Comm_Ex *,
        Dpi    *,
	int     ,
	int     ,
        double  ,
	int    );

EXTERN void surf_based_initialization 
( double *,
         double *,
         double *,
	 Exo_DB  *,
         int,
	 struct LS_Surf_List *,
         double,
         double,
         double );

EXTERN void free_surf_list
( struct LS_Surf_List **);

EXTERN struct LS_Surf * closest_surf
( struct LS_Surf_List *,
         double *,
         Exo_DB  *,
         double * );

EXTERN void find_surf_closest_point 
( struct LS_Surf *,
         double *,
         Exo_DB  *,
         double * );

EXTERN void find_facets
(  struct LS_Surf_List *,
	 int,
         double,
	 Exo_DB * );

EXTERN struct LS_Surf * create_surf
( int );

EXTERN struct LS_Surf * create_surf_point
( double *,
         int,
         double *,
         int );

EXTERN struct LS_Surf * create_surf_facet_line
( struct LS_Surf *,
         struct LS_Surf * );
         
EXTERN void append_surf
( struct LS_Surf_List *,
         struct LS_Surf * );

EXTERN void create_subsurfs
( struct LS_Surf_List *,
         double *,
         Exo_DB * );

EXTERN void free_subsurfs
( struct LS_Surf_List * );

EXTERN void assemble_Global_surf_list
( struct LS_Surf_List * );

                     
EXTERN int level_set_interface
( const double,          /*  F                         */
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
	double [DIM][DIM]);   /*  d_normal_dgradF[DIM][DIM] */
	
EXTERN int level_set_property
( const double,          /*  p0           */
	const double,          /*  p1           */
	const double,          /*  width        */
	double *,              /* *pp           */
	double [MDE]);        /*  d_pp_dF[MDE] */
                     
EXTERN int ls_transport_property
( const double,          /*  p0           */
	const double,          /*  p1           */
	const double,          /*  width        */
	double *,              /* *pp           */
	double * );           /*  d_pp_dF */

EXTERN void zero_lsi
( void);

EXTERN int load_lsi
( const double );       /* width */

EXTERN int load_lsi_shell_second
( const double );       /* width */

EXTERN int load_lsi_derivs
( void);


EXTERN double find_LS_global_flux
( const Exo_DB *,
	const Dpi *,
	const double *,
	      double *,
	double [],
	int );

EXTERN double find_LS_vel
( const Exo_DB *,       /* exo */
	const Dpi *,          /* dpi */
	const double *,       /* params[] */
	const int,            /* chosen_vel */
	double [],            /* x[] (solution vector) */
	int );               /* num_total_unknowns */

EXTERN void semi_lagrange_step
( int ,               /* num_total_nodes */ 
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
	Comm_Ex *);        /*  cx      */
 
#endif


EXTERN int boundary_curvature
( double [DIM],
	double [DIM][MAX_VARIABLE_TYPES+MAX_CONC][MDE],
	double * );

EXTERN int curvature_momentum_source
( double [DIM],         
	double [DIM][DIM][MDE],
	double [DIM][MDE],
	double [DIM][MDE] );

#ifdef CONST_CURVE
EXTERN int curvature_momentum_source_const
( double [DIM],         
	double [DIM][DIM][MDE],
	double [DIM][MDE],
	double [DIM][MDE] );
#endif


/* Phase function kernel */

EXTERN int assemble_phase_function
( double ,         
	double ,
	double ,
	double [DIM],
	Exo_DB *);
	   
	   
EXTERN int assemble_pf_constraint
( double, 
		double *,
		double  ,
		double *,
		double * );

EXTERN void assemble_ls_mass_lumped_penalty(
    struct LS_Mass_Lumped_Penalty *mass_lumped_penalty, int ip_total,
    int ielem_type, PG_DATA *pg_data);

#endif /* GOMA_MM_FILL_FILL_H */
