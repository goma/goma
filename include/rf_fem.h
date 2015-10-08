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
 *$Id: rf_fem.h,v 5.2 2009-03-12 00:00:56 hkmoffa Exp $
 *
 *
 * rf_fem.h:
 *
 *	Include file for globally occuring parameters and flags
 *	specific to the FEM problem. 
 *
 */

#ifndef _H_RF_FEM

#include "rf_fem_const.h"	/* In case you have not already done so. */
/*   max number of Interface Sources */
#ifndef MAX_INTERFACE
#define MAX_INTERFACE  5
#endif


/* Geometric Parameters  */ 

int  CoordinateSystem;  /* Indicates type of coordinate system (see fem_const.h)*/


/* FEM Interpolation Parameters (see fem_const.h) */

int  VelocityPressure;  /* Indicates which element type is used             */
	                /* for velocity and pressure interpolation.         */
int  Velocity;		/* Indicates which type of interpolation is used    */
			/* for velocity. (set by value of VelocityPressure) */
int  Pressure;		/* Indicates which type of interpolation is used    */
			/* for pressure. (set by value of VelocityPressure) */
int  Temperature;       /* Indicates which type of interpolation is used    */
			/* for temperature.                                 */
int  MeshDisplacement;  /* Indicates which element type is used             */
                        /* for mesh displacement interpolation.             */
int  MassFraction;      /* Indicates which type of interpolation is used    */
                        /* for mass fraction and density.                   */
int  nEQM;		/* Indicates one or more element quality metrics    */
			/* are to be performed.				    */
int  Use_DG;		/* Indicates when Discontinuous Galerkin	    */
                        /* inpterpolation is in use.			    */
int Do_Overlap;         /* Indicates that Overlap AC's will be used         */

/* Parameters to select Problem type (see fem_const.h)*/

int  ProblemType,	/* Select type of problem to be solved		    */
     ProblemCoupling,   /* Select fully coupled or dilute solution method   */
     StateEq,	        /* Select equation of state                         */
     Multicomponent;	/* Select fomulation for multicomponent transport   */


/* global variable to account for extra terms in axisymmetric or swirling
   flow problems: define in setup_pd */
int VIM;

/* Boundary Condition information */

int Num_BC;		/* number of boundary conditions which are defined  */

int Num_Interface_Srcs;	/* number *_D interfaces*/
int IntSrc_BCID[MAX_INTERFACE];
/* Rotation information */
int Num_ROT;		/* number of rotations which are defined  */

/* Import & Export field counts (used in coupled mode only) */
int Num_Import_NV;	/* number of nodal vars to import */
int Num_Import_EV;	/* number of nodal vars to import */
int Num_Export_XS;	/* number of solution vars to export */
int Num_Export_XP;	/* number of post-proc vars to export */
int Export_XS_ID[MAX_EXTERNAL_FIELD]; /* ID's of solution vars to export */
int Export_XP_ID[MAX_EXTERNAL_FIELD]; /* ID's of post proc vars to export */


/*
 * How many unique kinds of basis functions do we need to set up?
 * Some will be needed as Galerkin weights, some are used to interpolate
 * variables, and the elemental Jacobian matrix transforming the global
 * integral into local coordinates will use some kind of basis function.
 * For moving mesh problems, this will be the same as the interpolation
 * function for mesh displacement unknowns...
 *
 * At any rate, to economize the memory allocation, try to count up
 * from the input file specification exactly how many different basis
 * functions will be used.
 *
 * eg.; 2D fluid flow w/ mesh displacement
 *           v - Q2 (for momentum weight and for velocity interpolation)
 *           P - P1 (for continuity weight and for pressure interpolation)
 *           d - Q2 (for mesh stress weight and for displacement interpolation)
 *
 * Total: 2 unique_kinds of basis_functions
 */

int	Num_Basis_Functions;
int	Unique_Basis_Functions[MAX_BASIS_FUNCTIONS];

/*
 * Count up unique element types read in from the EXODUS II database, where
 * each element block has an element type associated with it.
 */

int	Num_Element_Types;
int	Unique_Element_Types[MAX_ELEMENT_TYPES];


/*
 * This is all very confusing...but...there really are basis function shapes
 * and interpolations. Together, they form a distinct basis function type.
 *
 * 
 */

int	Num_Interpolations;
int	Unique_Interpolations[MAX_INTERPOLATIONS];
int	Highest_Interpolation;

int	Num_Shapes;
int	Unique_Shapes[MAX_ELEMENT_SHAPES];


/* Parameters to select time integration technique                           */
int  TimeIntegration;   /* Select time integration method                    */
#ifndef COUPLED_FILL
int  Explicit_Fill;     /* Select time integration method                    */
int  exp_subcycle;      /* Subcycling frequency for Fill equation            */
#endif /* not COUPLED_FILL */
int  Use_Level_Set;     /* Global switch to turn on level set computations   */
int  Use_Phase_Field;   /* Global switch to turn on phase-field computations   */
		
int  MaxTimeSteps;      /* Maximum number of time steps in transient soution */
double  Delta_t0;	/* Initial time step                                 */
double  Delta_t_min;    /* Minimum time step size                            */
double  Delta_t_max;    /* Maximum time step size                            */
double  TimeMax;	/* Time at which to end integration                  */
/* double  theta;  */   /* Time step parameter: theta = 0. => Backward Euler
			                        theta = 1. => Forward  Euler */

double  eps;            /* Time step error                                   */
int  print_freq;
double print_delt;
double print_delt2_time, print_delt2;

/* Parameters for continuation */
int     Continuation;
int     ContType;
int     BdyCondID;
int     DataFltID;
int     MatID;
int     MatPropID;
int     MatPropSIID;
int     MaxPathSteps;
double  Delta_s0;
double  Delta_s_min;
double  Delta_s_max;
double  PathMax;

double  print_delt2_path;
double  BegParameterValue;
double  EndParameterValue;

/* Parameters for augmenting conditions */
int     nAC;

/* Parameters for multiple continuation conditions */
int     nCC, nTC, nUC, nUTC;

/* Parameters for hunting conditions */
int     nHC;

/* Internal Flags that indicate what needs to be calculated within the
   Residual and Coefficient Matrix Function                                 */

   /* The following flags can have TRUE and FALSE values */        

int  NumRegUnk_PN ;      /* Number of regular unknowns defined at each      */
                         /* global node.  This is defined as the number of  */
                         /* unknowns which have basis functions defined at  */
                         /* all global nodes in the element.  Unknowns      */
                         /* whose basis function interpolations are really  */
                         /* subparametrizations of the element are not      */
                         /* included, here.                                 */

/* Internal Parameters that calculate information needed in the
   Residual and Coefficient Matrix Function                                 */

int  NumSolnVelocities; /* Number of velocities to be solved for            */
                        /* (0, 1, 2, or 3)                                  */ 


/*
 * New quick reference to find the material index from an element block
 * index. This is more desirable for the distributed processing case, since
 * the old loops over the global number of materials will fail. Instead, 
 * use loops over the local number of element blocks that this processor
 * sees, then find the corresponding material via this integer array.
 *
 * Old:	   for ( mn=0; mn<pd_glob[0]->Num_Mat; mn ++)
 *            {
 *		.            
 *		.            
 *		.            
 * New:
 *	   for ( ebi=0; ebi<exo->num_elem_blocks; ebi++)
 *            {
 *               mn = Matilda[ebi];
 *
 *
 *	Remember! Think element block indeces, not material numbers!
 *
 */

extern int *Matilda;		/* defined and filled in rd_mesh.c */

/* Information on the number of unknowns (variables) which are define       */


extern int *num_internal_dofs;	/* I own, nobody wants. */

extern int *num_boundary_dofs;	/* I own, other procs want. */

extern int *num_external_dofs;	/* They own, I want. */

extern int *num_universe_dofs;	/* All the dofs this processor is aware of. 
				 * This is NOT the same as the number of
				 * degrees of freedom in the global problem.
				 * That number is considerably larger. */

/*
 * For convenience in assembling all internal and boundary nodal based 
 * equations, this processor traverses some of the same elements that other
 * processors traverse. However, for many purposes a unique element assignment
 * is needed. This is the number of elements on this processor that are 
 * assigned to this processor. It will never exceed exo->num_elems, for 
 * example.
 */

extern int num_personal_elems;

/*
 * ptr_node_recv[dp->num_neighbors+1]:
 *   This is a pointer list. It contains the offset
 *   of the start of each neighbor's node list in the
 *   list of external nodes for this processor.
 *   Remember that the list of external nodes in the global
 *   list of nodes may not be ordered
 *   by the ownership of that external node. We seek here to 
 *   create an ordered list with ptr_node_recv[] and
 *   list_node_recv[]. The ordering will be wrt the list of
 *   neighboring processors given by dpi->neighbor[].
 *
 *   ptr_node_recv[0] = 0
 *   ptr_node_recv[1] = a : a is the offset into the list of 
 *                          external nodes for the 
 *                          first node owned by this
 *                          processor's second neighbor,
 *                          dpi->neighbor[1].
 *   
 *   ptr_node_recv[dpi->num_neighbors] = dpi->num_external_nodes
 *
 */
extern int *ptr_node_recv;

extern int *ptr_fill_node_recv;

/*
 * This characterizes what it is that I have that others want from me.
 */
extern int *ptr_dof_send;
extern int *list_dof_send;
extern int *ptr_node_send;
extern int *list_node_send;
extern int *ptr_fill_node_send;
extern int *list_fill_node_send;


extern int *NumUnknowns; /* Number of unknown variables updated by this   */
			/* processor                                     */
extern int *NumExtUnknowns;/* Number of external variables which are      */
			/* copied and stored on the local processor      */
extern int MaxVarPerNode;/* Global Maximum number of unknowns at any     */
			/* node on any processor                         */
extern int Num_Var_In_Type[MAX_NUM_MATRICES][MAX_VARIABLE_TYPES];
                       /* The number of variables of this type  in the   */
                       /* current problem. For species variables types,  */
                       /* there are either  zero or Max_Num_Species_Eqn  */
                       /* species variables. For other variables, there  */
                       /* are either zero or one, usually.               */

/*
 *   FILL variable type variables for each processor
 */
extern int num_fill_unknowns;      /*  Number of FILL unknowns           */
extern int owned_fill_unknowns;    /*  Number of owned FILL unknowns     */
extern int internal_fill_unknowns; /*  Number of internal FILL unknowns  */
extern int boundary_fill_unknowns; /*  Number of boundary FILL unknowns  */
extern int external_fill_unknowns; /*  Number of external FILL unknowns  */

/*
 * Local_Offset:
 *
 *	This array of pointers gets set to lists of integers that, for every 
 *	node, indicates how much offset from the first unknown of the node
 *	you need to get to the variable of interest...
 *
 *	For example, if we're solving just the energy equation and temperature
 *	is interpolated at every node, then 
 *
 *		Local_Offset[matrix][node][TEMPERATURE] = 0; (the first unknown is T)
 *	and
 *		Local_Offset[matrix][node][VELOCITY1] = -1; (undefined offset)
 *	
 *	Yes, this does duplicate some functionality of First_Y, First_MeshD,
 *	etc. Also, Index_P, will not really be needed anymore, since pressure
 *	is getting lumped together with other unknowns at a node.
 *	
 */
extern int ***Local_Offset;

/*
 * Dolphin:  The scheme above is fine except for cases that occur when
 *	     neighboring elements that share a node have different ideas
 *	     about which variables and degrees of freedom need to be solved
 *	     at the node.
 *
 *	     In order to get the true accounting of offsets for each type of
 *	     variable at a node, we need to know what every element expects
 *	     of the node.
 *
 *
 *	     Thus, poll every element and save the largest estimate for
 *	     the degrees of freedom required to represent each variable...
 *
 *	     Dolphin[matrix][node][variable] = dof
 *
 * where:
 *		node  == global node number
 *
 *		var   == variable index
 *
 *		dof   == number of degrees of freedom for this type of
 *			 variable at this node.
 *
 *			 Note that variables with multiple k types must all
 *			 have the same representation in terms of dof/node.
 *			 Thus, the number of degrees of freedom for each
 *			 species concentrations must be multiplied by the 
 * 			 total number of concentrations that are active, too!
 */
extern int 	***Dolphin;	

/* Information about element Topology */
int num_vertex;
int num_curve;
int num_surface;
int num_body;
int *body_nodes;

/*
 * dofname -- holds strings telling the name of the variable (u1, T, P, etc)
 *            and the global node number associated with a particular gdof.
 *            This should aid in debugging, etc. Allocation and setup in
 *	      mm_unknown_map.c.
 *
 * Now, idv[matrix][dof][0] = VELOCITY1, etc.
 *      idv[matrix][dof][1] = local nodal dof (0, except pressure& conc., for example)
 *	idv[matrix][dof][2] = associated global node number (0-based)
 */
extern int  ***idv;    	    /* Integer variable name, nodal dof, node. */
extern char ***dofname;	    /* Names of variables. */
extern char ***resname;	    /* Names of residual equations. */

#define _H_RF_FEM

#endif
