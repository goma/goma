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

#ifndef _MM_FILL_LS_H
#define _MM_FILL_LS_H

#ifdef EXTERN
#undef EXTERN
#endif

#ifdef _MM_FILL_LS_C
#define EXTERN /* do nothing */
#endif

#ifndef _MM_FILL_LS_C
#define EXTERN extern
#endif

#ifndef COUPLED_FILL
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

#endif /*COUPLED_FILL */
 

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

EXTERN void free_subsurfs
PROTO (( struct LS_Surf_List * ));

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

EXTERN void create_subsurfs
PROTO (( struct LS_Surf_List *,
         double *,
         Exo_DB * ));

EXTERN void assemble_Global_surf_list
PROTO(( struct LS_Surf_List * ));

EXTERN int sign_change
PROTO(( double,
        double ));

EXTERN int elem_on_isosurface
PROTO(( int,
        double [],
        const Exo_DB *,
        int,
        double ));

EXTERN int current_elem_on_isosurface
PROTO(( int,
        double ));


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

EXTERN void find_facets
PROTO((  struct LS_Surf_List *,
	 int,
         double,
	 Exo_DB * ));
         
EXTERN struct LS_Surf * create_surf
PROTO (( int ));

EXTERN void append_surf
PROTO (( struct LS_Surf_List *,
         struct LS_Surf * ));

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
PROTO (( double *, 
         double *,
         int,
         double,
         Integ_Elem * ));

EXTERN struct LS_Surf_List * create_surf_list
PROTO (( void ));

EXTERN struct LS_Surf * create_surf_point
PROTO (( double *,
         int,
         double *,
         int ));

EXTERN struct LS_Surf * create_surf_facet_line
PROTO (( struct LS_Surf *,
         struct LS_Surf *,
         int,
         int ));

EXTERN void ls_var_initialization 
PROTO (( double *, 
	 Exo_DB *, 
	 Dpi *, 
	 Comm_Ex * ));
         

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



EXTERN double ls_modulate_property
PROTO (( double ,
	 double ,
	 double ,
	 double ,
	 double ,
	 double [MDE],
	 double * ));

EXTERN int is_xfem_interp
PROTO(( const int));		/* interp type */

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
PROTO(( double [],
        const Exo_DB  * ));

EXTERN void load_xfem_for_stu
PROTO(( const double [] ));


EXTERN void xfem_correct
PROTO(( int ,                  /* num_total_nodes    */
        double [],             /* x[]                */
        double [],             /* xdot[]             */
        double [],             /* x_old[]            */
        double [],             /* xdot_old[]         */
        double [],             /* delta_x[]         */
        double ,               /* theta_arg          */
        double  ));            /* delta_t            */

EXTERN void xfem_predict
PROTO(( int ,                  /* num_total_nodes    */
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
        double [] ));          /* xdot_older[]    */

EXTERN void xfem_var_diff
PROTO(( int ,
        double *,
        double [MDE],
        double [DIM] ));

EXTERN void zero_lsi
PROTO(( void));

EXTERN void zero_lsi_derivs
PROTO(( void));

EXTERN int load_lsi
PROTO(( const double ));       /* width */

EXTERN int load_lsi_adjmatr
PROTO(( const double ));      /* width */

EXTERN int load_lsi_derivs
PROTO(( void));

#ifndef COUPLED_FILL

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

#endif /*COUPLED_FILL */

EXTERN int assemble_level_project
PROTO((double [],               /* Jacobian matrix for fill equation  */
       int    [],               /* pointer to nonzeros in Jacobian matrix   */
       double [],               /* rhs vector   */
       double ,                 /* current time step size */
       double ,                 /* parameter to vary time integration from
                                 * explicit (tt = 1) to implicit (tt = 0) */
       int [] ));               /* node_to_fill -  */

EXTERN int assemble_level_correct
PROTO((double [],                /* Jacobian matrix for fill equation  */
       int    [],                /* pointer to nonzeros in Jacobian matrix   */
       double [],                /* rhs vector   */
       double ,                  /* current time step size */
       double ,                  /* parameter to vary time integration from
                                  * explicit (tt = 1) to implicit (tt = 0) */
       int [] ));                /* node_to_fill -  */

EXTERN void print_surf_list 
PROTO(( struct LS_Surf_List *,
        double ));

EXTERN void append_surf_isosurf 
PROTO(( struct LS_Surf_List *,
	int, 
	double ));

EXTERN void print_point_list
PROTO(( double *,
	Exo_DB *,
	char *, 
	double ));
        
EXTERN int generate_facet_list
PROTO(( double (**)[DIM],
        double (**)[DIM],
        int ** ,
        double *,
	Exo_DB * ));
        
EXTERN int dof_incomplete
PROTO(( int,
        int,
        int,
        int ));

EXTERN void determine_ls_elem_overlap_state
PROTO(( void ));

EXTERN void xfem_dof_state
PROTO(( const int, /* ledof */
        const int, /* interpolation type */
        const int, /* element shape */
        int *,     /* flag indicating xfem affects this dof's basis functions */
        int *,     /* flag indicating if this an extended dof */
        int *,     /* base interpolation, ie, I_Q1_XG -> I_Q1 */
        int * ));  /* what dof of base_interp does this dof map to */

EXTERN int is_extended_dof
PROTO(( const int ,                          /* I */
        const int ,                          /* idof */
        VARIABLE_DESCRIPTION_STRUCT *,
	const double ));
                     
EXTERN void assemble_interface_extension_velocity
PROTO(( dbl [],
        Exo_DB *,
        Dpi *));

EXTERN void assemble_boundary_extension_velocity
PROTO(( double [],
        Exo_DB *,
        Dpi * ));

EXTERN int assemble_extension_velocity
PROTO(( dbl [DIM],
	dbl [DIM][DIM],
	dbl [DIM][MDE]));



NTREE *Subgrid_Tree;                          /* This is a global pointer to the subgrid integration shape function tree */
NTREE_INT Subgrid_Int;                        /* This is a global structure for the subgrid integration points and weights specific to element */

SGRID *create_search_grid
PROTO (( NTREE *  ));

void divide_search_grid 
PROTO(( SGRID *,
	int  ));

EXTERN void free_search_grid 
PROTO(( SGRID ** ));

EXTERN void free_shape_fcn_tree 
PROTO (( NTREE * ));

void map_local_coordinates
PROTO(( double *,
	double * ));

void print_search_grid
PROTO (( SGRID * ));

void find_grid_intersections
PROTO (( SGRID *, struct LS_Surf_List * ));

NTREE *create_shape_fcn_tree
PROTO(( int ));

int build_integration_grid
PROTO(( SGRID *,
	int,
	double ));

int find_tree_integration_pts
PROTO (( NTREE *,
	 int ));

void compute_tree_size
PROTO (( NTREE *,
	 double *,
	 double *,
	 double *,
	 double * ));

int grid_overlaps_interface
PROTO (( SGRID *,
	 double ));

int gather_integration_pts 
PROTO (( SGRID *,
	 double (*)[DIM],
	 double *,
	 int ));


EXTERN int current_elem_overlaps_interface
PROTO(( double  ));

EXTERN int elem_overlaps_interface
PROTO(( int elem,
        double x[],
        const Exo_DB* exo,
        double width ));

EXTERN void map_subelement_stu
PROTO(( double *,
        Integ_Elem *,
        double * ));
                    
void subelement_J
PROTO(( Integ_Elem *,
        double *,
        double [DIM][DIM],
        int ));

EXTERN double subelement_detJ
PROTO(( Integ_Elem *,
        double *,
	int ));

EXTERN double subelement_surfdet
PROTO(( Integ_Elem * ,
        double * ,
        int,
	int ));
                                                                                           
EXTERN int get_subgrid_integration_pts 
PROTO(( NTREE *,
		SGRID **, 
		double (**)[DIM], 
		double **, 
		double  ));

EXTERN int get_surface_subgrid_integration_pts
PROTO(( SGRID *,
		double [DIM],
		double (*)[DIM],
		double *,
		int ));

EXTERN int print_subgrid_integration_pts
PROTO(( double(*)[DIM],
        double *, 
	int));

EXTERN int print_subgrid_surface_integration_pts
PROTO(( double (*)[DIM],
        double *,
        int ));

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
PROTO(( Integ_Elem_Desc ** ));

EXTERN void get_subelement_descriptions
PROTO(( double x[],
        Exo_DB *exo,
	Integ_Elem_Desc_List *list ));
	
EXTERN void ls_surface_extents
PROTO(( double x[],
	Exo_DB *exo,
	double min[3],
	double max[3] ));

EXTERN void gather_subelement_descriptions
PROTO(( Integ_Elem_Desc_List *,
        Integ_Elem * ));

EXTERN double Courant_Time_Step
PROTO(( double [],
        double [],
        double [],
        double [],
        double [],
        double [],
        int *,
	Exo_DB *exo ));
        
EXTERN void subelement_mesh_output
PROTO(( double [],
	Exo_DB *exo ));

EXTERN int get_facet_integration_pts
PROTO(( double (**)[DIM],
        double **,
        Exo_DB * ));
                                                    
EXTERN int get_subelement_integration_pts
PROTO(( double (**)[DIM],
	double **,
        int **,
	double,
        int,
        int ));

EXTERN void get_subelement_facets
PROTO(( struct LS_Surf_List *,
	double ));

void
gather_subelement_facets
PROTO(( struct LS_Surf_List *,
        Integ_Elem * ));

EXTERN Integ_Elem * create_integ_elements
PROTO(( double ));

void
free_integ_elements
PROTO(( Integ_Elem * ));
        
EXTERN void build_integ_element
PROTO(( Integ_Elem *,
        double,
        int,
        double (*)[DIM],
        int *,
	int ));

int
num_subelement_integration_pts
PROTO(( Integ_Elem *,
        int,
        int ));

EXTERN int gather_subelement_integration_pts
PROTO(( Integ_Elem *,
        double (*)[DIM],
	double *,
        int *,
        int,
        int,
        int ));

EXTERN int gather_surface_subgrid_integration_pts 
PROTO(( SGRID *, 
		int , 
		double [DIM], 
		double (*)[DIM], 
		double *, 
		int  ));


EXTERN void map_subelement_stu
PROTO(( double *,
        Integ_Elem *,
        double * ));


EXTERN void subelement_J
PROTO(( Integ_Elem *,
        double *,
        double [DIM][DIM],
        int  ));
              
EXTERN void iso_contour_on_side 
PROTO(( double,            /* isoval */
	int,               /* dim */
	int,               /* ielem_type */
	int,               /* id_side  */
	int *,             /* ip_total */
	double (**)[DIM],  /* (**s)[DIM] */
	double ** ));      /* **wt  */

EXTERN void clear_xfem_contribution
PROTO(( int ));


EXTERN void compute_xfem_contribution
PROTO(( int ));


EXTERN void check_xfem_contribution
PROTO(( int,
        struct Aztec_Linear_Solver_System *,
        double [],
        double [],
	Exo_DB * ));

EXTERN void resolve_ls_adc_old 
PROTO (( struct Boundary_Condition *,
		 Exo_DB *,
		 double *,
		 double ,
		 int * ,
		 int ));

EXTERN struct LS_Surf * resolve_ls_adc
PROTO (( struct LS_Surf_List *, 
		 struct Boundary_Condition *,
		 Exo_DB *,
		 double *,
		 double ,
		 int * ,
		 int ));



#endif /* ifndef _MM_FILL_LS_H */
