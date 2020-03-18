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
 *$Id: mm_fill_ls.h,v 5.5 2008-03-18 19:01:58 hkmoffa Exp $
 */

#ifndef GOMA_MM_FILL_LS_H
#define GOMA_MM_FILL_LS_H

#include "dp_types.h"
#include "dpi.h"
#include "el_elm.h"
#include "exo_struct.h"
#include "mm_as_structs.h"
#include "mm_fill_fill.h"
#include "rf_vars_const.h"
#include "std.h"

struct Aztec_Linear_Solver_System;
struct Boundary_Condition;
struct elem_side_bc_struct;
#ifdef EXTERN
#undef EXTERN
#endif

#ifdef GOMA_MM_FILL_LS_C
#define EXTERN /* do nothing */
#endif

#ifndef GOMA_MM_FILL_LS_C
#define EXTERN extern
#endif

#ifndef COUPLED_FILL
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

#endif /*COUPLED_FILL */
 

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

EXTERN void free_subsurfs
( struct LS_Surf_List * );

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

EXTERN void create_subsurfs
( struct LS_Surf_List *,
         double *,
         Exo_DB * );

EXTERN void assemble_Global_surf_list
( struct LS_Surf_List * );

EXTERN int sign_change
( double,
        double );

EXTERN int elem_on_isosurface
( int,
        double [],
        const Exo_DB *,
        int,
        double );

EXTERN int
elem_near_isosurface ( int elem,
                     double x[],
                     const Exo_DB* exo,
                     int isovar,
                     double isoval );

EXTERN int current_elem_on_isosurface
( int,
        double );


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

EXTERN void find_facets
(  struct LS_Surf_List *,
	 int,
         double,
	 Exo_DB * );
         
EXTERN struct LS_Surf * create_surf
( int );

EXTERN void append_surf
( struct LS_Surf_List *,
         struct LS_Surf * );

struct Integ_Elem_Struct {
  int ielem_type;
  int num_local_nodes;
  double (* xi)[DIM];
  double * f;
  int num_subelements;
  struct Integ_Elem_Struct ** subelements;
  int num_sides;
  int * side_ids;
  int * bc_sides;
  int sign;
  int is_conformal;
};

typedef struct Integ_Elem_Struct Integ_Elem;
         
EXTERN int find_link_intersection
( double *, 
         double *,
         int,
         double,
         Integ_Elem * );

EXTERN struct LS_Surf_List * create_surf_list
( void );

EXTERN struct LS_Surf * create_surf_point
( double *,
         int,
         double *,
         int );

EXTERN struct LS_Surf * create_surf_facet_line
( struct LS_Surf *,
         struct LS_Surf *,
         int,
         int );

EXTERN void ls_var_initialization 
( double **,
	 Exo_DB *, 
	 Dpi *, 
	 Comm_Ex ** );
         

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
                     

EXTERN int level_set_property_offset
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



EXTERN double ls_modulate_property
( double ,
	 double ,
	 double ,
	 double ,
	 double ,
	 double [MDE],
	 double * );

double
ls_modulate_property_offset ( double p1,
                       double p2,
                       double width,
                       double pm_minus,
                       double pm_plus,
                       double dpdF[MDE],
                       double *factor );

EXTERN int is_xfem_interp
( const int);		/* interp type */

struct Extended_Shape_Fcn_Basics {
   struct Level_Set_Data *ls;
   int ielem;
   int elem_state;
   int elem_var_state;
   int node_var_state[MDE];
   double xi[DIM];
   int near;
   int have_XG;
   double F;
   double H;
   double delta;
   double dF_xi[DIM];
   int have_disc;
   double F_zero;
   double bf_plus;
   double F_plus;
   double grad_bf_plus[DIM];
   double grad_F_plus[DIM];
   double *active_vol;
   double *tot_vol;
};

struct Extended_Shape_Fcn_Basics * xfem;       /* This is a global structure for the basic pieces needed for XFEM */

EXTERN void load_xfem_for_elem
( double [],
        const Exo_DB  * );

EXTERN void load_xfem_for_stu
( const double [] );


EXTERN void xfem_correct
( int ,                  /* num_total_nodes    */
        double [],             /* x[]                */
        double [],             /* xdot[]             */
        double [],             /* x_old[]            */
        double [],             /* xdot_old[]         */
        double [],             /* delta_x[]         */
        double ,               /* theta_arg          */
        double  );            /* delta_t            */

EXTERN void xfem_predict
( int ,                  /* num_total_nodes    */
        int ,                  /* numProcUnknowns    */
        double ,               /* delta_t    */
        double ,               /* delta_t_old    */
        double ,               /* delta_t_older    */
        double ,               /* theta_arg    */
        double [],             /* x[]    */
        double [],             /* x_old[]    */
        double [],             /* x_older[]    */
        double [],             /* x_oldest[]    */
        double [],             /* xdot[]    */
        double [],             /* xdot_old[]    */
        double [] );          /* xdot_older[]    */

EXTERN void xfem_var_diff
( int ,
        double *,
        double [MDE],
        double [DIM] );

EXTERN void zero_lsi
( void);

EXTERN void zero_lsi_derivs
( void);

EXTERN int load_lsi
( const double );       /* width */

EXTERN int load_lsi_offset
( const double );       /* width */

EXTERN int
load_lsi_old(const double width, struct Level_Set_Interface *lsi_old);

EXTERN int load_lsi_adjmatr
( const double );      /* width */

EXTERN int load_lsi_derivs
( void);

#ifndef COUPLED_FILL

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

#endif /*COUPLED_FILL */

EXTERN int assemble_level_project
(double [],               /* Jacobian matrix for fill equation  */
       int    [],               /* pointer to nonzeros in Jacobian matrix   */
       double [],               /* rhs vector   */
       double ,                 /* current time step size */
       double ,                 /* parameter to vary time integration from
                                 * explicit (tt = 1) to implicit (tt = 0) */
       int [] );               /* node_to_fill -  */

EXTERN int assemble_level_correct
(double [],                /* Jacobian matrix for fill equation  */
       int    [],                /* pointer to nonzeros in Jacobian matrix   */
       double [],                /* rhs vector   */
       double ,                  /* current time step size */
       double ,                  /* parameter to vary time integration from
                                  * explicit (tt = 1) to implicit (tt = 0) */
       int [] );                /* node_to_fill -  */

EXTERN int print_ls_interface( double *x,
			       Exo_DB *exo,
			       Dpi    *dpi,
			       const double time,
			       char *filenm,
			       int print_all_times );

EXTERN void print_surf_list 
( struct LS_Surf_List *,
        double );

EXTERN void append_surf_isosurf 
( struct LS_Surf_List *,
	int, 
	double );

EXTERN void print_point_list
( double *,
	Exo_DB *,
	char *, 
	double );
        
EXTERN int generate_facet_list
( double (**)[DIM],
        double (**)[DIM],
        int ** ,
        double *,
	Exo_DB * );
        
EXTERN int dof_incomplete
( int,
        int,
        int,
        int );

EXTERN void determine_ls_elem_overlap_state
( void );

EXTERN void xfem_dof_state
( const int, /* ledof */
        const int, /* interpolation type */
        const int, /* element shape */
        int *,     /* flag indicating xfem affects this dof's basis functions */
        int *,     /* flag indicating if this an extended dof */
        int *,     /* base interpolation, ie, I_Q1_XG -> I_Q1 */
        int * );  /* what dof of base_interp does this dof map to */

EXTERN int is_extended_dof
( const int ,                          /* I */
        const int ,                          /* idof */
        VARIABLE_DESCRIPTION_STRUCT *,
	const double );
                     
EXTERN void assemble_interface_extension_velocity
( dbl [],
        Exo_DB *,
        Dpi *);

EXTERN void assemble_boundary_extension_velocity
( double [],
        Exo_DB *,
        Dpi * );

EXTERN int assemble_extension_velocity
( dbl [DIM],
	dbl [DIM][DIM],
	dbl [DIM][MDE]);



NTREE *Subgrid_Tree;                          /* This is a global pointer to the subgrid integration shape function tree */
NTREE_INT Subgrid_Int;                        /* This is a global structure for the subgrid integration points and weights specific to element */

SGRID *create_search_grid
( NTREE *  );

void divide_search_grid 
( SGRID *,
	int  );

EXTERN void free_search_grid 
( SGRID ** );

EXTERN void free_shape_fcn_tree 
( NTREE * );

void map_local_coordinates
( double *,
	double * );

void print_search_grid
( SGRID * );

void find_grid_intersections
( SGRID *, struct LS_Surf_List * );

NTREE *create_shape_fcn_tree
( int );

int build_integration_grid
( SGRID *,
	int,
	double );

int find_tree_integration_pts
( NTREE *,
	 int );

void compute_tree_size
( NTREE *,
	 double *,
	 double *,
	 double *,
	 double * );

int grid_overlaps_interface
( SGRID *,
	 double );

int gather_integration_pts 
( SGRID *,
	 double (*)[DIM],
	 double *,
	 int );


EXTERN int current_elem_overlaps_interface
( double  );

EXTERN int elem_overlaps_interface
( int elem,
        double x[],
        const Exo_DB* exo,
        double width );

EXTERN void map_subelement_stu
( double *,
        Integ_Elem *,
        double * );
                    
EXTERN double subelement_detJ
( Integ_Elem *,
        double *,
	int );

EXTERN double subelement_surfdet
( Integ_Elem * ,
        double * ,
        int,
	int );
                                                                                           
EXTERN int get_subgrid_integration_pts 
( NTREE *,
		SGRID **, 
		double (**)[DIM], 
		double **, 
		double  );

EXTERN int get_surface_subgrid_integration_pts
( SGRID *,
		double [DIM],
		double (*)[DIM],
		double *,
		int );

EXTERN int print_subgrid_integration_pts
( double(*)[DIM],
        double *, 
	int);

EXTERN int print_subgrid_surface_integration_pts
( double (*)[DIM],
        double *,
        int );

struct Integ_Elem_Desc_Struct {
  int num_nodes;
  int num_sides;
  double (* x)[DIM];
  int * bc_sides;
  struct Integ_Elem_Desc_Struct * next;
};

typedef struct Integ_Elem_Desc_Struct Integ_Elem_Desc;

struct Integ_Elem_Desc_List_Struct {
  Integ_Elem_Desc * start;
  int size;
};

typedef struct Integ_Elem_Desc_List_Struct Integ_Elem_Desc_List;

EXTERN void free_subelement_descriptions
( Integ_Elem_Desc ** );

EXTERN void get_subelement_descriptions
( double x[],
        Exo_DB *exo,
	Integ_Elem_Desc_List *list );
	
EXTERN void ls_surface_extents
( double x[],
	Exo_DB *exo,
	double min[3],
	double max[3] );

EXTERN void gather_subelement_descriptions
( Integ_Elem_Desc_List *,
        Integ_Elem * );

EXTERN double Courant_Time_Step
( double [],
        double [],
        double [],
        double [],
        double [],
        double [],
        int *,
	Exo_DB *exo );
        
EXTERN void subelement_mesh_output
( double [],
	Exo_DB *exo );

EXTERN int get_facet_integration_pts
( double (**)[DIM],
        double **,
        Exo_DB * );
                                                    
EXTERN int get_subelement_integration_pts
( double (**)[DIM],
	double **,
        int **,
	double,
        int,
        int );

EXTERN void get_subelement_facets
( struct LS_Surf_List *,
	double );

void
gather_subelement_facets
( struct LS_Surf_List *,
        Integ_Elem * );

EXTERN Integ_Elem * create_integ_elements
( double );

void
free_integ_elements
( Integ_Elem * );
        
EXTERN void build_integ_element
( Integ_Elem *,
        double,
        int,
        double (*)[DIM],
        int *,
	int );

int
num_subelement_integration_pts
( Integ_Elem *,
        int,
        int );

EXTERN int gather_subelement_integration_pts
( Integ_Elem *,
        double (*)[DIM],
	double *,
        int *,
        int,
        int,
        int );

EXTERN int gather_surface_subgrid_integration_pts 
( SGRID *, 
		int , 
		double [DIM], 
		double (*)[DIM], 
		double *, 
		int  );


EXTERN void subelement_J
( Integ_Elem *,
        double *,
        double [DIM][DIM],
        int  );
              
EXTERN void iso_contour_on_side 
( double,            /* isoval */
	int,               /* dim */
	int,               /* ielem_type */
	int,               /* id_side  */
	int *,             /* ip_total */
	double (**)[DIM],  /* (**s)[DIM] */
	double ** );      /* **wt  */

EXTERN void clear_xfem_contribution
( int );


EXTERN void compute_xfem_contribution
( int );


EXTERN void check_xfem_contribution
( int,
        struct Aztec_Linear_Solver_System *,
        double [],
        double [],
	Exo_DB * );

EXTERN void resolve_ls_adc_old 
( struct Boundary_Condition *,
		 Exo_DB *,
		 double *,
		 double ,
		 int * ,
		 int );

EXTERN struct LS_Surf * resolve_ls_adc
( struct LS_Surf_List *, 
		 struct Boundary_Condition *,
		 Exo_DB *,
		 double *,
		 double ,
		 int * ,
		 int );



#endif /* ifndef GOMA_MM_FILL_LS_H */
