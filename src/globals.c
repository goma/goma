/************************************************************************ *
* Goma - Multiphysics finite element software                             *
* Sandia National Laboratories                                            *
*                                                                         *
* Copyright (c) 2022 Goma Developers, National Technology & Engineering   *
*               Solutions of Sandia, LLC (NTESS)                          *
*                                                                         *
* Under the terms of Contract DE-NA0003525, the U.S. Government retains   *
* certain rights in this software.                                        *
*                                                                         *
* This software is distributed under the GNU General Public License.      *
* See LICENSE file.                                                       *
\************************************************************************/
/* File containing the many of Goma's global variables, moved for -fno-common */

#include "ac_particles.h"
#include "rf_fem.h"
#include "rf_fem_const.h"
#include "rf_node_const.h"
#include "rf_solver.h"

int CoordinateSystem; /* Indicates type of coordinate system (see fem_const.h)*/

/* FEM Interpolation Parameters (see fem_const.h) */

int VelocityPressure; /* Indicates which element type is used             */
/* for velocity and pressure interpolation.         */
int Velocity; /* Indicates which type of interpolation is used    */
/* for velocity. (set by value of VelocityPressure) */
int Pressure; /* Indicates which type of interpolation is used    */
/* for pressure. (set by value of VelocityPressure) */
int Temperature; /* Indicates which type of interpolation is used    */
/* for temperature.                                 */
int MeshDisplacement; /* Indicates which element type is used             */
/* for mesh displacement interpolation.             */
int MassFraction; /* Indicates which type of interpolation is used    */
/* for mass fraction and density.                   */
int nEQM; /* Indicates one or more element quality metrics    */
/* are to be performed.				    */
int Use_DG; /* Indicates when Discontinuous Galerkin	    */
/* inpterpolation is in use.			    */
int Do_Overlap; /* Indicates that Overlap AC's will be used         */

/* Parameters to select Problem type (see fem_const.h)*/

int ProblemType,     /* Select type of problem to be solved		    */
    ProblemCoupling, /* Select fully coupled or dilute solution method   */
    StateEq,         /* Select equation of state                         */
    Multicomponent;  /* Select fomulation for multicomponent transport   */

/* global variable to account for extra terms in axisymmetric or swirling
   flow problems: define in setup_pd */
int VIM;
int WIM;

/* Boundary Condition information */

int Num_BC; /* number of boundary conditions which are defined  */

int Num_Interface_Srcs; /* number *_D interfaces*/
int IntSrc_BCID[MAX_INTERFACE];
/* Rotation information */
int Num_ROT; /* number of rotations which are defined  */

/* Import & Export field counts (used in coupled mode only) */
int Num_Import_NV;                    /* number of nodal vars to import */
int Num_Import_EV;                    /* number of nodal vars to import */
int Num_Export_XS;                    /* number of solution vars to export */
int Num_Export_XP;                    /* number of post-proc vars to export */
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

int Num_Basis_Functions;
int Unique_Basis_Functions[MAX_BASIS_FUNCTIONS];

/*
 * Count up unique element types read in from the EXODUS II database, where
 * each element block has an element type associated with it.
 */

int Num_Element_Types;
int Unique_Element_Types[MAX_ELEMENT_TYPES];

/*
 * This is all very confusing...but...there really are basis function shapes
 * and interpolations. Together, they form a distinct basis function type.
 *
 *
 */

int Num_Interpolations;
int Unique_Interpolations[MAX_INTERPOLATIONS];
int Highest_Interpolation;

int Num_Shapes;
int Unique_Shapes[MAX_ELEMENT_SHAPES];

/* Parameters to select time integration technique                           */
int TimeIntegration; /* Select time integration method                    */
int Use_Level_Set;   /* Global switch to turn on level set computations   */
int Use_Phase_Field; /* Global switch to turn on phase-field computations   */

/* double  theta;  */ /* Time step parameter: theta = 0. => Backward Euler
                                              theta = 1. => Forward  Euler */

double eps; /* Time step error                                   */
int print_freq;
double print_delt;
double print_delt2_time, print_delt2;

/* Parameters for continuation */
int Continuation;
int ContType;
int BdyCondID;
int DataFltID;
int MatID;
int MatPropID;
int MatPropSIID;
int MaxPathSteps;
double Delta_s0;
double Delta_s_min;
double Delta_s_max;
double PathMax;

double print_delt2_path;
double BegParameterValue;
double EndParameterValue;

/* Parameters for augmenting conditions */
int nAC;

/* Parameters for multiple continuation conditions */
int nCC, nTC, nUC, nUTC;

/* Parameters for hunting conditions */
String_line Matrix_Solver;

String_line Matrix_Format;

String_line Matrix_Scaling;

String_line Matrix_Preconditioner;

String_line Matrix_Subdomain_Solver; /* new Aztec 2.x option */

String_line Matrix_Residual_Norm_Type;

String_line Matrix_Output_Type;

String_line Matrix_Factorization_Reuse;

String_line Matrix_Graph_Fillin; /* new Aztec 2.x option */

String_line Matrix_Maximum_Iterations;

String_line Matrix_Polynomial_Order;

String_line Matrix_Factor_Overlap;

String_line Matrix_Overlap_Type;

String_line Matrix_Krylov_Subspace;

String_line Matrix_Orthogonalization;

String_line Matrix_Auxiliary_Vector;

String_line Matrix_Convergence_Tolerance;

String_line Matrix_Drop_Tolerance;

String_line Matrix_Factorization_Save; /* Aztec 2 */

String_line Matrix_ILUT_Fill_Factor; /* Aztec 2 */

String_line Matrix_RILU_Relax_Factor; /* Aztec 2 */

String_line Matrix_BILU_Threshold; /* Trilinos 1 */

String_line Matrix_Relative_Threshold; /* Trilinos 2 */

String_line Matrix_Absolute_Threshold; /* Trilinos 2 */

String_line Amesos_Package;

String_line Amesos2_Package;

String_line AztecOO_Solver;

String_line Stratimikos_File[MAX_NUM_MATRICES];

String_line Amesos2_File[MAX_NUM_MATRICES];

/*
 * A new Aztec 2.0 option. There are more and difft options and our
 * previous options probably ought to be revised to reflect the newer
 * Aztec 2.0 capability specifications.
 */

String_line Matrix_Reorder;

int Linear_Solver; /* Aztec, Sparse, MA28, UMFPACK */

int UMFPACK_IDIM; /* UMFPACK STORAGE CONSTANT */
int UMFPACK_XDIM; /* UMFPACK STORAGE CONSTANT */
int LOCA_UMF_ID;  /* UMFPACK SYSTEM ID */

int Max_Newton_Steps;  /* Maximum number of Newton steps to take.     */
int Guess_Flag;        /* Indicates the type of initial guess         */
int Conformation_Flag; /* Indicates mapping from stress to log-conformation tensor */
int Print3DBCDup;

double damp_factor;
double damp_factor1; /* Relaxation factor for Newton iteration */
/* damp_factor1 = 1.0 is full Newton */
/* damp_factor1 = 0.0 is not updating our */
/*                   solution estimate */
double damp_factor2, /* Additional damping factors for custom */
    damp_factor3,    /* schemes for automatic control with    */
    custom_tol1,     /* NORM(0,0) tolerances                  */
    custom_tol2, custom_tol3;
double var_damp[MAX_VARIABLE_TYPES]; /* variable specific damp factors */

int Newt_Jacobian_Reformation_stride; /*Stride for reformation of jacobian for
                                   modified newton scheme               */
int Time_Jacobian_Reformation_stride;
int Newton_Line_Search_Type;
int modified_newton;               /*boolean flag for modified Newton */
int save_old_A;                    /*boolean flag for saving old A matrix
                                    for resolve reasons with AZTEC.   There
                                    are at least four reasons, that you
                                    can see in sl_util.c */
double convergence_rate_tolerance; /* tolerance for jacobian reformation
                                       based on convergence rate */
double modified_newt_norm_tol;     /* tolerance for jacobian reformation
                                           based on residual norm */

double Epsilon[MAX_NUM_MATRICES][3]; /* Used for determining stopping criteria.     */
int Solver_Output_Format;            /* Bitmap for Solver Output Format     */
int Output_Variable_Stats;           /* Toggle for Variable Stats Output    */
int Output_Variable_Regression;      /* Toggle for Variable Regression    */

int NZeros; /* Number of nonzeros in this procs matrix     */

int nnz_own; /* number of nonzeroes in the matrix that
              * belong to this processor
              */

int nnz_total; /* total number of nonzeroes in the matrix
                * that this processor sees, including the
                * entries for external unknowns that will
                * not be updated by this processor  */

int GNZeros; /* Number of nonzeros in global matrix         */

int fill_zeros; /* number of nonzeros in fill matrix for this
                  processor */

int Gfill_zeros; /* number of nonzeros in fill matrix for the
                   global problem */

int PSPG;          /* 1 means pressure stabilized Petrov-Galerkin is used */
int PSPP;          /* 1 means pressure stabilized polynomial projection is used */
double PS_scaling; /* This term is a constant scaling for the PSPG or PSPP term */
int Cont_GLS;      /* 1 means continuity stabilization is used */

int Filter_Species, filter_species_material_number;
double c_min, c_max;

int Include_Visc_Sens, Visc_Sens_Copy, Visc_Sens_Factor;
/* 1 means to include the sensitivities of the
viscosity functions in the jacobian matrix.
0 indicates that these sensitivities are not to
be included.  The latter case is useful for highly
shear thinning viscosity models.  By disabling
the viscosity sensititives convergence can be
achieved albeit at a less than quadratic rate
      */
int nHC;

int PRS_mat_ielem;

struct Boundary_Condition *BC_Types;
struct Rotation_Specs *ROT_Types;

NODE_INFO_STRUCT **Nodes = NULL;

struct elem_side_bc_struct ***First_Elem_Side_BC_Array;
struct elem_edge_bc_struct ***First_Elem_Edge_BC_Array;

/* Global variables extern declared in ac_particles.h. */
int Particle_Dynamics;                /* global toggle indicating particles are present. */
enum Particle_Model_t Particle_Model; /* What flavor of particle<->continuum stuff... */
dbl Particle_Model_Data[MAX_PARTICLE_MODEL_DATA_VALUES]; /* Real values for this model. */
int Particle_Number;                                     /* number of discrete particles. */
particle_filename_s Particle_Restart_Filename;           /* restart filename */
int Particle_Output_Stride;    /* How often to output particle information. */
dbl Particle_Output_Time_Step; /* Output every these units. */
int Particle_Max_Time_Steps;   /* max number of particle time steps if steady solution. */
enum Particle_Output_Format_t Particle_Output_Format; /* What kind of output file? */
int Particle_Full_Output_Stride;                   /* > 0 => full output every that many steps. */
particle_filename_s Particle_Full_Output_Filename; /* where to put them. */
int Particle_Number_Sample_Types;                  /* How many datasets to output? */
int *Particle_Number_Samples_Existing;             /* How many are tagged for each sample type?*/
int *Particle_Number_Samples;          /* How many particles to output for dataset #n? */
int *Particle_Number_Output_Variables; /* How many output vars for each sample. */
particle_variable_s *
    *Particle_Output_Variables; /* List of variable indices to output for dataset #n */
particle_filename_s *Particle_Filename_Template; /* Template of where to put the data... */

dbl Particle_Density;         /* Density of particle in problem units */
dbl Particle_Radius;          /* Radius of particle in problem units. */
dbl Particle_Ratio;           /* Real/computational particle ratio. */
int Particle_Show_Debug_Info; /* Show particle debug info. */
enum Particle_Domain_t Particle_Creation_Domain;
enum Particle_Domain_t Particle_Move_Domain;
particle_filename_s Particle_Creation_Domain_Filename;
particle_s Particle_Creation_Domain_Name;
particle_filename_s Particle_Move_Domain_Filename;
particle_s Particle_Move_Domain_Name;
dbl Particle_Creation_Domain_Reals[MAX_DOMAIN_REAL_VALUES];
dbl Particle_Move_Domain_Reals[MAX_DOMAIN_REAL_VALUES];
dbl xi_boundary_tolerances[3] = {XI_BOUNDARY_TOLERANCE0, XI_BOUNDARY_TOLERANCE1,
                                 XI_BOUNDARY_TOLERANCE2};

int Particle_Number_PBCs; /* Number of particle-related sideset BC's. */
PBC_t *PBCs;              /* Particle boundary condition structures. */

element_particle_info_t *element_particle_info;