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

/*
 *$Id: mm_post_proc.c,v 5.15 2010-07-21 16:39:27 hkmoffa Exp $
 */

/* Standard include files */

#include "mm_fill_em.h"
#include "rf_solve.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h> /* strcasecmp and strncasecmp moved here for POSIX.1 */

#define GOMA_MM_POST_PROC_C
#include "mm_post_def.h"
#include "mm_post_proc.h"
/* GOMA include files */
#include "ac_particles.h"
#include "bc/rotate.h"
#include "bc_contact.h"
#include "dpi.h"
#include "el_elm.h"
#include "el_elm_info.h"
#include "el_geom.h"
#include "exo_struct.h"
#include "exodusII.h"
#include "mm_as.h"
#include "mm_as_alloc.h"
#include "mm_as_const.h"
#include "mm_as_structs.h"
#include "mm_bc.h"
#include "mm_elem_block_structs.h"
#include "mm_fill_aux.h"
#include "mm_fill_common.h"
#include "mm_fill_ls.h"
#include "mm_fill_population.h"
#include "mm_fill_porous.h"
#include "mm_fill_ptrs.h"
#include "mm_fill_rs.h"
#include "mm_fill_solid.h"
#include "mm_fill_species.h"
#include "mm_fill_stabilization.h"
#include "mm_fill_stress.h"
#include "mm_fill_terms.h"
#include "mm_fill_util.h"
#include "mm_input.h"
#include "mm_more_utils.h"
#include "mm_mp.h"
#include "mm_mp_const.h"
#include "mm_mp_structs.h"
#include "mm_qtensor_model.h"
#include "mm_shell_util.h"
#include "mm_std_models_shell.h"
#include "mm_unknown_map.h"
#include "mm_viscosity.h"
#include "rd_mesh.h"
#include "rf_allo.h"
#include "rf_bc.h"
#include "rf_bc_const.h"
#include "rf_element_storage_struct.h"
#include "rf_fem.h"
#include "rf_fem_const.h"
#include "rf_io.h"
#include "rf_io_const.h"
#include "rf_io_structs.h"
#include "rf_mp.h"
#include "rf_shape.h"
#include "rf_solver.h"
#include "rf_util.h"
#include "rf_vars_const.h"
#include "shell_tfmp_struct.h"
#include "shell_tfmp_util.h"
#include "sl_aux.h"
#include "sl_petsc.h"
#include "std.h"
#include "user_mp.h"
#include "user_post.h"
#include "wr_exo.h"

/*
 * Global variable definitions.
 * This is the 1 place these variables are defined. If you need them
 * elsewhere, then include them via mm_post_def.h
 */
int nn_post_fluxes;      /* Dimension of the following structure */
int nn_post_fluxes_sens; /* Dimension of the following structure */
int nn_post_data;        /* Dimension of the following structure */
int nn_post_data_sens;   /* Dimension of the following structure */
int nn_error_metrics;    /* Dimension of the following structure */
int nn_particles;        /* Dimension of the following structure */
int nn_volume;           /* number of pp_volume_int structures */
int ppvi_type;           /* Maybe there's a better way to do this, Seems like globals abound! AMC*/
int nn_global; /* Hopefully these post processing thin get refactored into non-globals at some point
                */
int nn_average;
static int *listel; /* Pointer to element order map from exoII  */

struct Post_Processing_Data **pp_data;
struct Post_Processing_Data_Sens **pp_data_sens;
struct Post_Processing_Error *pp_error_data;
struct Post_Processing_Fluxes **pp_fluxes;
struct Post_Processing_Fluxes_Sens **pp_fluxes_sens;
struct Post_Processing_Particles **pp_particles;
struct Post_Processing_Volumetric **pp_volume;
pp_Global **pp_global;
pp_Average **pp_average;

static int error_presence_key[3]; /* Truth key (dim 3) as to which error
                                   * measure elem sizes are being done */

int Num_Nodal_Post_Proc_Var = 0;
int Num_Elem_Post_Proc_Var = 0;

/* The following variables are flags for input options - i.e. 0 or 1,
 * but they are set to the post-processing variable number in load_nodal_tkn
 * to be used in mm_post_proc.c  - for options which imply more than one post-
 * processing variable to be output, the flag becomes the post-processing
 * variable number of the first one of this variable type
 *
 * Hey! All upper case should denote predefined constants from the
 * preprocessor, not *variables*! Just don't let it happen again!
 */

/*
 * Legend:
 *
 *  n>0 -- "Yes, calculate this quantity during post processing." Typically, 1.
 *
 *   -1 -- "No, do not calculate this quantity during post processing."
 */

int CAPILLARY_PRESSURE = -1; /* capillary pressure in a porous media */
int CONC_CONT = -1;          /* concentration at vertex & midside nodes*/
int CONDUCTION_VECTORS = -1; /* conduction flux vectors*/
int CURL_V = -1;             /* Steve Kempka's favorite quantity */
int DARCY_VELOCITY_GAS = -1; /* Darcy velocity vectors for gas phase
                              * flow in a partially saturated porous
                              * media */
int DARCY_VELOCITY_LIQ = -1; /* Darcy velocity vectors for flow in a
                              * saturated or unsaturated medium */
int DENSITY = -1;            /* density function at vertex and midside
                              * nodes, e.g. for particle settling etc. */
int POLYMER_VISCOSITY = -1;
int POLYMER_TIME_CONST = -1;
int PTT_XI = -1;
int PTT_EPSILON = -1;
int MOBILITY_PARAMETER = -1;

int HEAVISIDE = -1;
int DIELECTROPHORETIC_FIELD = -1;
/* Dielectrophoretic force vectors. */
int DIELECTROPHORETIC_FIELD_NORM = -1;
/* Norm of DIELECTROPHORETIC_FIELD. */
int ENORMSQ_FIELD;          /* grad(|E|^2), for dielectrophoresis. */
int ENORMSQ_FIELD_NORM;     /* |grad(|E|^2)|, for dielectrophoresis. */
int DIFFUSION_VECTORS = -1; /* Diffusive mass flux vectors (units ?) */
int DIFFUSION_VECTORS_POR_LIQ_GPHASE = -1;
/* Diffusion of the solvent liquid in the
 * gas phase for porous flow problems */
int DIFFUSION_VECTORS_POR_AIR_GPHASE = -1;
/* Diffusion of the air in the
 * gas phase for porous flow problems */
int DIV_PVELOCITY = -1;      /* check the divergence of the particle phase
                              * velocities.  */
int DIV_TOTAL = -1;          /* Divergence of the sum of the fluid and
                              * particle phases.  This should be zero. */
int DIV_VELOCITY = -1;       /* incompressibility constraint at vertex and
                              * midside nodes, e.g. del*v = 0 */
int ELECTRIC_FIELD = -1;     /* Electric field vectors: E = -grad(VOLTAGE) */
int ELECTRIC_FIELD_MAG = -1; /* Electric field magnitude: sqrt(E.E) */

int ENERGY_FLUXLINES = -1;    /* energy flux function, analogous to
                               * stream function ... */
int ERROR_ZZ_P = -1;          /* Zienkiewicz-Zhu error indicator (element
                               * quantity) based solely on pressure
                               * contributions */
int ERROR_ZZ_P_ELSIZE = -1;   /* Recommended new element size from ZZ
                               * pressure measure                          */
int ERROR_ZZ_Q = -1;          /* Zienkiewicz-Zhu error indicator (element
                               * quantity) based solely on heat flux
                               * contributions                             */
int ERROR_ZZ_Q_ELSIZE = -1;   /* Recommended new element size from ZZ heat
                               * flux measure */
int ERROR_ZZ_VEL = -1;        /* Zienkiewicz-Zhu error indicator (element
                               * quantity) based solely on velocity (shear
                               * stress) contributions                     */
int ERROR_ZZ_VEL_ELSIZE = -1; /* Recommended new element size from ZZ
                               * velocity measure */
int SAT_CURVE_TYPE = -1;
int SAT_QP_SWITCH = -1;    /* Value of sat. function at hysteretic
                            * swtch point */
int CAP_PRESS_SWITCH = -1; /* Value of cap. press function at hysteretic
                            * swtch point */
int NUM_CURVE_SWITCH = -1; /* Number of hysteretic curve switch */
int SAT_HYST_MIN = -1;     /* Minimum saturation value for scanning imbibition curve */
int SAT_HYST_MAX = -1;     /* Maximum saturation value for scanning draining curve */

int EVP_DEF_GRAD_TENSOR = -1;
int EXTERNAL_POST = -1; /* external field variables read from other
                         * files */
int FILL_CONT = -1;     /* fill at vertex & midside nodes*/
int FIRST_INVAR_STRAIN = -1;
int FLUXLINES = -1;           /* mass flux function. This is analogous to
                               * stream function but represents mass flux */
int LAGRANGE_CONVECTION = -1; /* Lagrangian convection velocity */
int MEAN_SHEAR = -1;
int MM_RESIDUALS = -1; /* stress equation residuals at vertex
                        * and midside nodes*/
int NS_RESIDUALS = -1; /* Navier-Stokes residuals at vertex
                        * and midside nodes */
int POROUS_RHO_GAS_SOLVENTS = -1;
/* gas phase concentration of each solvent
 * species in a porous media */
int POROUS_RHO_LPHASE = -1; /* liquid phase density per unit volume of
                             * material in a porous medium */
int POROUS_RHO_TOTAL_SOLVENTS = -1;
/* Total density of each solvent species in a
 * porous media. Total, here means the
 * density summed up over all phases */
int POROUS_SATURATION = -1;   /* saturation in a porous media */
int POROUS_GRIDPECLET = -1;   /* Grid Peclet number for porous media */
int POROUS_SUPGVELOCITY = -1; /* Effective velocities to use in SUPG
                               * formulations in porous media */
int POROUS_LIQUID_ACCUM_RATE = -1;
/* The rate at which liquid in a partially
 * saturated porous medium is accumulating
 * at a point */
int REL_LIQ_PERM = -1; /* Relative liquid permeability in porous media */

int PRESSURE_CONT = -1;   /* pressure at vertex & midside nodes*/
int SH_DIV_S_V_CONT = -1; /* SH_DIV_S_V at midside nodes */
int SH_CURV_CONT = -1;    /* SH_SURF_CURV at midside nodes */
int REAL_STRESS_TENSOR = -1;
int SEC_INVAR_STRAIN = -1;     /* 2nd strain invariant vertex,midside nodes*/
int STRAIN_TENSOR = -1;        /* strain tensor for mesh deformation  */
int STREAM = -1;               /* stream function*/
int STREAM_NORMAL_STRESS = -1; /* streamwise normal stress function*/
int STREAM_SHEAR_STRESS = -1;  /* streamwise shear stress function*/
int STREAM_TENSION = -1;       /* streamwise tension - (Ttt-Tnn)/R*/
int STRESS_CONT = -1;          /* stress at vertex & midside nodes*/
int STRESS_TENSOR = -1;        /* stress tensor for mesh deformation
                                * (Lagrangian pressure) */
int SURFACE_VECTORS = -1;      /* vector field of normals and tangents on
                                * surfaces, curves and vertices */
int SHELL_NORMALS = -1;        /* vector field of smoothed normals computed from fv->sh_ang */
int THIRD_INVAR_STRAIN = -1;
int TIME_DERIVATIVES = -1;           /* time derivatives */
int TOTAL_STRESS11 = -1;             /* sum over all modes for multi-mode models */
int TOTAL_STRESS12 = -1;             /* sum over all modes for multi-mode models */
int TOTAL_STRESS13 = -1;             /* sum over all modes for multi-mode models */
int TOTAL_STRESS22 = -1;             /* sum over all modes for multi-mode models */
int TOTAL_STRESS23 = -1;             /* sum over all modes for multi-mode models */
int TOTAL_STRESS33 = -1;             /* sum over all modes for multi-mode models */
int USER_POST = -1;                  /* a user defined function */
int PP_Viscosity = -1;               /* Viscosity */
int PP_FlowingLiquid_Viscosity = -1; /* Flowing Liquid Viscosity - Porous Brinkman term */
int PP_VolumeFractionGas = -1; /* Volume fraction of gas in a foam or other two phase material */
int ACOUSTIC_PRESSURE = -1;
int ACOUSTIC_PHASE_ANGLE = -1;
int ACOUSTIC_ENERGY_DENSITY = -1;
int LIGHT_INTENSITY = -1;
int PRINCIPAL_STRESS = -1;
int PRINCIPAL_REAL_STRESS = -1;
int LUB_HEIGHT = -1;
int LUB_HEIGHT_2 = -1;
int LUB_VELO_UPPER = -1;
int LUB_VELO_LOWER = -1;
int LUB_VELO_FIELD = -1;
int LUB_VELO_FIELD_2 = -1;
int LUB_FLUID_SOURCE = -1;
int DISJ_PRESS = -1;
int SH_SAT_OPEN = -1;
int SH_SAT_OPEN_2 = -1;
int SH_CAP_PRES = -1;
int SH_PORE_FLUX = -1;
int SH_STRESS_TENSOR = -1;
int SH_TANG = -1;
int PP_LAME_MU = -1;
int PP_LAME_LAMBDA = -1;
int VON_MISES_STRESS = -1;
int VON_MISES_STRAIN = -1;
int UNTRACKED_SPEC = -1;
int CONF_MAP = -1;
int RHO_DOT = -1;
int MOMENT_SOURCES = -1;
int YZBETA = -1;
int TFMP_GAS_VELO = -1;
int TFMP_LIQ_VELO = -1;
int TFMP_INV_PECLET = -1;
int TFMP_KRG = -1;
int VELO_SPEED = -1;
int GIES_CRIT = -1;
int J_FLUX = -1;
int EIG = -1;
int EIG1 = -1;
int EIG2 = -1;
int EIG3 = -1;
int GRAD_SH = -1;
int GRAD_Y = -1;
int HELICITY = -1;
int LAMB_VECTOR = -1;
int Q_FCN = -1;
int POYNTING_VECTORS = -1; /* conduction flux vectors*/
int PSPG_PP = -1;
int SARAMITO_YIELD = -1;
int STRESS_NORM = -1;
int SPECIES_SOURCES = -1; /* continuous species sources*/
int VISCOUS_STRESS = -1;
int VISCOUS_STRESS_NORM = -1;
int VISCOUS_VON_MISES_STRESS = -1;
int EM_CONTOURS = -1;
int TOTAL_EM_CONTOURS = -1;
int SCATTERED_EM_CONTOURS = -1;
int ORIENTATION_VECTORS = -1; /* orientation vectors*/
int FIRST_STRAINRATE_INVAR = -1;
int SEC_STRAINRATE_INVAR = -1;
int THIRD_STRAINRATE_INVAR = -1;

int len_u_post_proc = 0; /* size of dynamically allocated u_post_proc
                          * actually is */
double *u_post_proc = 0; /* user-provided values used in calculating
                          * user defined post processing variable */

/*
 *  Post-processing Step 1: add a new variable flag to end of mm_post_def.h
 *
 *       e.g.  extern int STREAM;
 *
 *       Note that this flag is now -1 (false) or the id number of the
 *       post-processing variable in load_nodal_tkn
 *
 *                  Step 2: add the definition in mm_post_proc.c
 *
 *       e.g.  int STREAM;
 */

/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
/*
 *  Global Variables Defined in this file:
 */
static int *listel;

/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

/*
 * Prototypes of static functions
 */

static int calc_standard_fields /* mm_post
                                   _proc.c                            */
    (double **post_proc_vect,
     double **lumped_mass,
     double delta_t,
     double theta,
     int ielem,
     const int ielem_type,
     int ip,
     int ip_total,
     RESULTS_DESCRIPTION_STRUCT *rd,
     struct Porous_Media_Terms *pmt,
     double time,
     Exo_DB *exo,
     double xi[DIM],
     const PG_DATA *pg_data);

static int calc_zz_error_vel /* mm_post_proc.c                            */
    (double[],               /* x - Soln vector for the current processor */
     double[],
     double[],
     double[],
     double[],
     int,           /* ev_indx - Variable index for zz_error     */
     double ***,    /* gvec_elem - element variable values that  *
                     * are indexed: [eb_indx][ev_indx][elem]     */
     Exo_DB *const, /* exo                                       */
     Dpi *const,    /* dpi                                       */
     int);          /* compute_elem_size                         */

static int abs_error_at_elem /* mm_post_proc.c                            */
    (int,                    /* i_elem                                    */
     double ***,             /* tau_nodal_lsp                             */
     double *,               /* gvec_elem                                 */
     double *,               /* velocity_norm                             */
     double *);              /* ielem_area                                */

static int fill_lhs_lspatch /* mm_post_proc.c                            */
    (double *,              /* i_node_coords                             */
     double *,              /* wt_gp_loc                                 */
     double *,              /* xgp_loc                                   */
     double *,              /* ygp_loc                                   */
     double *,              /* zgp_loc                                   */
     double *,              /* det_gp_loc                                */
     int,                   /* max_terms                                 */
     double **,             /* s_lhs                                     */
     int,                   /* elem                                      */
     double ****);          /* tau_gp_ptch                               */

static int calc_stream_fcn   /* mm_post_proc.c                            */
    (double[],               /* x                             soln vector */
     double[4],              /* del_stream_fcn                            */
     double[MAX_PDIM][MDE]); /* vel - array for local nodal velocity     *
                              * values which must be divergence free     */

static int correct_stream_fcn /* mm_post_proc.c                            */
    (int *,                   /* kount - counter for element connectivity  */
     int,                     /* iel - current element number              */
     double[4],               /* del_stream_fcn - elemental side increments*
                               * to stream function calculated by          *
                               * calc_stream_fcn()                         */
     double[],                /* stream_fcn_vect                           */
     int[]);                  /* listnd - count times node is accessed     */

static int look_for_post_proc /* mm_post_proc.c                            */
    (FILE *,                  /* ifp - pointer to file                     */
     char *,                  /* search_string -                           */
     int *);                  /* flag_variable - integer flag for options  */

static int midsid(double[],  /* stream_fcn_vect */
                  Exo_DB *); /* exo */

void sum_average_nodal(double **avg_count, double **avg_sum, int global_node, double time);

void post_process_average(double x[],            /* Solution vector for the current processor */
                          double x_old[],        /* Solution vector at last time step */
                          double xdot[],         /* time derivative of solution vector */
                          double xdot_old[],     /* dx/dt at previous time step */
                          double resid_vector[], /* Residual vector for the
                                                    current proc */
                          Exo_DB *exo,
                          Dpi *dpi,
                          double **post_proc_vect,
                          double time);

/*
 * Prototypes of functions defined in other files that are needed here.
 */

/*
 * Here's a RECIPE for adding new post-processing options to make it easier to
 * add new ones (note that these lines are repeated at each place where you
 * need to make changes):
 *
 * Post-processing Step
 *
 *  [1] Add a new variable flag to end of mm_post_def.h. For example
 *
 *		EXTERN int FOOBAR;
 *
 *      where the integer flag FOOBAR will be -1 if it is not to be computed,
 *      or it will be some integer that denotes its identification to the
 *      routine load_nodal_tkn().
 *  [1.5] Also add the definition in mm_post_proc.c For example
 *
 *		int FOOBAR;
 *
 *  [2] Add a new line in mm_post_proc.c routine rd_post_process_specs()
 *      to search input file for indications that your new variable should
 *      be activated. For example,
 *
 *		iread = look_for_post_proc(ifp, "Foo bar", &FOOBAR, '\n');
 *
 *  [3] Add a new line to put your variable's name into the exodus II database
 *      in mm_post_proc.c, routine load_nodal_tkn(). For the sake of argument,
 *      suppose the scalar quantity "foobar" requires the energy equation
 *      to be active.
 *
 *      if ( FOOBAR != -1 && Num_Var_In_Type[pg->imtrx][R_ENERGY])
 *         {
 *            set_nv_tkud(rd, index, 0, 0, "FB","[1]", "Foobar", FALSE);
 *            index++;
 *            FOOBAR = index_post;
 *            index_post++;
 *         }
 *
 *  [4] Make provision to send the flag to calculate FOOBAR to any other
 *      processors. For simple Boolean flags, add a line in dp_vif.c for
 *      the main body transport via Noahs Ark like:
 *
 *		ddd_add_member(n, &FOOBAR, 1, MPI_INT);
 *
 *  [5] Add algorithm to calculate your new variable in mm_post_proc.c
 *      and put it in post_proc_vect[] array.
 *
 *  [6] Use and enjoy your new post-processing variable!
 *
 *
 * Notes: All these changes can be made by editing only 4 files:
 *
 *              mm_post_def.h,
 *		mm_post_proc.c,
 *		mm_post_proc.h,
 *		dp_vif.c
 *
 *        Based on Rich Cairncross' initial recipe.
 *
 *        Revised 1997/11/08 15:21 MST pasacki@sandia.gov
 */

/*________________________________________________________________________*/

static int calc_standard_fields(double **post_proc_vect,
                                double **lumped_mass,
                                double delta_t,
                                double theta,
                                int ielem,
                                const int ielem_type,
                                int ip,
                                int ip_total,
                                RESULTS_DESCRIPTION_STRUCT *rd,
                                struct Porous_Media_Terms *pmt,
                                double time,
                                Exo_DB *exo,
                                double xi[3],
                                const PG_DATA *pg_data)

/****************************************************************************
 *
 * calc_standard_fields()
 *
 * -- calculate post-processing quantities at a gauss point.
 *    Then, project their values onto all of the local nodes in the
 *    local element. This routine is called at the element-gauss point
 *    level.
 *
 ****************************************************************************/
{
  int dim, a, b, c, eqn, var, w, w1, i, j, I, index, status, mode;
  dbl det_J; /* determinant of Jacobian: has "r" in it for
                axisymmetric case, no "r" for Cartesian*/
  dbl phi_i; /* Weighting functions for i-th residuals */
  dbl phi_j; /* Interpolation functions for variables */
  dbl wt;
  dbl n[MAX_CONC][DIM]; /* Diffusion flux vector. */
  dbl qc[DIM];          /* Conduction flux vector. */
  int err, ldof, ldof_right, ileft, iright, midside, checkPorous;
  /*  int first_porous_var = POR_LAST - MAX_POROUS_NUM + 1;*/
  /*  int SPECIES = MAX_VARIABLE_TYPES; */

  double TT[MAX_PDIM][MAX_PDIM];
  double EE[MAX_PDIM][MAX_PDIM];
  double FVP[MAX_PDIM][MAX_PDIM];
  double TVP[MAX_PDIM][MAX_PDIM];
  dbl dTT_drs[DIM][DIM][DIM][MDE];
  double dTT_dx[MAX_PDIM][MAX_PDIM][MAX_PDIM][MDE];
  double dTT_dp[MAX_PDIM][MAX_PDIM][MDE];
  double dTT_dc[MAX_PDIM][MAX_PDIM][MAX_CONC][MDE];
  dbl dTT_dp_liq[DIM][DIM][MDE];    /* Sensitivity of stress tensor...
                                        to nodal porous liquid pressure*/
  dbl dTT_dp_gas[DIM][DIM][MDE];    /* Sensitivity of stress tensor...
                                        to nodal porous gas pressure*/
  dbl dTT_dporosity[DIM][DIM][MDE]; /* Sensitivity of stress tensor...
                                     to nodal porosity*/
  dbl dTT_dsink_mass[DIM][DIM][MDE];
  double dTT_dT[MAX_PDIM][MAX_PDIM][MDE];
  double dTT_dmax_strain[MAX_PDIM][MAX_PDIM][MDE];
  double dTT_dcur_strain[MAX_PDIM][MAX_PDIM][MDE];

  /*  double elast_modulus;	*/
  double stream_grad, velo_sqrd, vdelvdx[DIM], curv;
  double vconv[MAX_PDIM];     /*Calculated convection velocity */
  double vconv_old[MAX_PDIM]; /*Calculated convection velocity at previous time*/
  CONVECTION_VELOCITY_DEPENDENCE_STRUCT d_vconv_struct;
  CONVECTION_VELOCITY_DEPENDENCE_STRUCT *d_vconv = &d_vconv_struct;

  /*
   * Variables for vicosity
   */
  dbl gamma[DIM][DIM];
  dbl mu;

  /* polymer viscosity */
  dbl mup;

  /* shift factor for polymer viscosity */

  dbl at = 0, wlf_denom;

  dbl rho; /* density variables */

  dbl hs[DIM];

  int species_particle; /* For SUSPENSION_PM, which species is the
                         * particle phase? */
  dbl p_vol_frac;       /* Particle volume fraction. */

  struct Species_Conservation_Terms s_terms;

  /* DG VE stuff */
  int v_s[MAX_MODES][DIM][DIM];
  dbl ves[DIM][DIM];

  /* Variables for dielectrophoretic force calculations. */
  dbl coeff_a, coeff_b, coeff_c, coeff_d, CM_fact, dielectrophoretic_force_coeff;
  dbl dfvector[MAX_PDIM], dfnorm; /* Enormsqvector[MAX_PDIM];*/

  /*
   * Additional variables for projection
   */
  dbl Ttt, E_E, TrE, Dnn, ts, Tnt, nv[DIM] = {0., 0., 0.};
  int i_pg, i_pl, i_pore;

  double *local_post, *local_lumped;

  status = 0;
  /* Some initialization */
  memset(dTT_dx, 0, sizeof(double) * DIM * DIM * DIM * MDE);
  memset(dTT_drs, 0, sizeof(double) * DIM * DIM * DIM * MDE);
  memset(dTT_dp, 0, sizeof(double) * DIM * DIM * MDE);
  memset(dTT_dc, 0, sizeof(double) * DIM * DIM * MAX_CONC * MDE);
  memset(dTT_dp_liq, 0, sizeof(double) * DIM * DIM * MDE);
  memset(dTT_dp_gas, 0, sizeof(double) * DIM * DIM * MDE);
  memset(dTT_dporosity, 0, sizeof(double) * DIM * DIM * MDE);
  memset(dTT_dsink_mass, 0, sizeof(double) * DIM * DIM * MDE);
  memset(dTT_dT, 0, sizeof(double) * DIM * DIM * MDE);
  memset(dTT_dmax_strain, 0, sizeof(double) * DIM * DIM * MDE);
  memset(dTT_dcur_strain, 0, sizeof(double) * DIM * DIM * MDE);
  memset(FVP, 0, sizeof(double) * DIM * DIM);

  zero_structure(&s_terms, sizeof(struct Species_Conservation_Terms), 1);

  /* load eqn and variable number in tensor form */
  err = stress_eqn_pointer(v_s);

  i_pl = 0;
  i_pg = 1;
  i_pore = 2;

  /* Porous media shell variables */
  dbl S;
  dbl d_cap_pres[2], cap_pres;
  d_cap_pres[0] = d_cap_pres[1] = 0.;
  dbl Patm;

  local_post = alloc_dbl_1(rd->TotalNVPostOutput, 0.0);
  local_lumped = alloc_dbl_1(rd->TotalNVPostOutput, 0.0);

  /*
   * Unpack variables from structures for local convenience...
   */

  dim = pd_glob[0]->Num_Dim;

  if (pd->v[pg->imtrx][R_MESH1] && ei[pg->imtrx]->ielem_dim >= dim) {
    err = belly_flop(elc->lame_mu);
    GOMA_EH(err, "error in belly flop");
    if (err == 2)
      exit(-1);
  }

  /*
   * Compute desired quantities at current gauss point and add them
   * into the post_proc_vector array
   */

  if (STREAM_NORMAL_STRESS != -1 && pd->e[pg->imtrx][R_MOMENTUM1]) {
    velo_sqrd = 0.;
    stream_grad = 0.;

    for (a = 0; a < VIM; a++) {
      for (b = 0; b < VIM; b++) {
        gamma[a][b] = fv->grad_v[a][b] + fv->grad_v[b][a];
      }
    }
    for (a = 0; a < VIM; a++) {
      velo_sqrd += fv->v[a] * fv->v[a];
      for (b = 0; b < VIM; b++) {
        /* vv:gamma */
        stream_grad += fv->v[a] * gamma[a][b] * fv->v[b];
      }
    }
    stream_grad *= mp->viscosity;
    /*
     **  viscoelastic stress tensor
     **/
    if (pd->v[pg->imtrx][POLYMER_STRESS11]) {
      dbl log_c[DIM][DIM];
      dbl exp_s[DIM][DIM];
      dbl R1[DIM][DIM];
      dbl eig_values[DIM];
      dbl mup = 0.;
      dbl lambda = 0.;
      int ve_mode, p, r;
      memset(ves, 0, sizeof(dbl) * DIM * DIM);
      if (vn->evssModel == LOG_CONF || vn->evssModel == LOG_CONF_GRADV) {
        for (ve_mode = 0; ve_mode < vn->modes; ve_mode++) {
          for (p = 0; p < VIM; p++) {
            for (r = 0; r < VIM; r++) {
              log_c[p][r] = fv->S[ve_mode][p][r];
            }
          }
#ifdef ANALEIG_PLEASE
          analytical_exp_s(log_c, exp_s, eig_values, R1, NULL);
#else
          compute_exp_s(log_c, exp_s, eig_values, R1);
#endif
          mup = viscosity(ve[ve_mode]->gn, gamma, NULL);
          if (ve[ve_mode]->time_constModel == CONSTANT) {
            lambda = ve[ve_mode]->time_const;
          }
          for (a = 0; a < VIM; a++) {
            for (b = 0; b < VIM; b++) {
              ves[a][b] += mup / lambda * (exp_s[a][b] - (double)delta(a, b));
            }
          }
        }
      } else {
        for (a = 0; a < VIM; a++) {
          for (b = 0; b < VIM; b++) {
            for (ve_mode = 0; ve_mode < vn->modes; ve_mode++) {
              ves[a][b] += fv->S[ve_mode][a][b];
            }
          }
        }
      }
      for (a = 0; a < VIM; a++) {
        for (b = 0; b < VIM; b++) {
          stream_grad += fv->v[a] * ves[a][b] * fv->v[b];
        }
      }
    } // if pd->v[pg->imtrx][POLYMER_STRESS11]
    if (DOUBLE_NONZERO(velo_sqrd)) {
      Ttt = stream_grad / velo_sqrd;
    } else {
      Ttt = 0.0;
    }

    local_post[STREAM_NORMAL_STRESS] = Ttt;
    local_lumped[STREAM_NORMAL_STRESS] = 1.;
  }
  if (STREAM_SHEAR_STRESS != -1 && pd->e[pg->imtrx][R_MOMENTUM1]) {
    nv[0] = fv->v[1];
    nv[1] = -fv->v[0];
    velo_sqrd = 0.;
    stream_grad = 0.;

    for (a = 0; a < VIM; a++) {
      for (b = 0; b < VIM; b++) {
        gamma[a][b] = fv->grad_v[a][b] + fv->grad_v[b][a];
      }
    }
    for (a = 0; a < dim; a++) {
      velo_sqrd += fv->v[a] * fv->v[a];
      for (b = 0; b < dim; b++) {
        /* vv:gamma */
        stream_grad += fv->v[a] * gamma[a][b] * nv[b];
      }
    }
    stream_grad *= mp->viscosity;
    /*
     **  viscoelastic stress tensor
     **/
    if (pd->v[pg->imtrx][POLYMER_STRESS11]) {
      dbl log_c[DIM][DIM];
      dbl exp_s[DIM][DIM];
      dbl R1[DIM][DIM];
      dbl eig_values[DIM];
      dbl mup = 0.;
      dbl lambda = 0.;
      int ve_mode, p, r;
      memset(ves, 0, sizeof(dbl) * DIM * DIM);
      if (vn->evssModel == LOG_CONF || vn->evssModel == LOG_CONF_GRADV) {
        for (ve_mode = 0; ve_mode < vn->modes; ve_mode++) {
          for (p = 0; p < VIM; p++) {
            for (r = 0; r < VIM; r++) {
              log_c[p][r] = fv->S[ve_mode][p][r];
            }
          }
#ifdef ANALEIG_PLEASE
          analytical_exp_s(log_c, exp_s, eig_values, R1, NULL);
#else
          compute_exp_s(log_c, exp_s, eig_values, R1);
#endif
          mup = viscosity(ve[ve_mode]->gn, gamma, NULL);
          if (ve[ve_mode]->time_constModel == CONSTANT) {
            lambda = ve[ve_mode]->time_const;
          }
          for (a = 0; a < VIM; a++) {
            for (b = 0; b < VIM; b++) {
              ves[a][b] += mup / lambda * (exp_s[a][b] - (double)delta(a, b));
            }
          }
        }
      } else {
        for (a = 0; a < VIM; a++) {
          for (b = 0; b < VIM; b++) {
            for (ve_mode = 0; ve_mode < vn->modes; ve_mode++) {
              ves[a][b] += fv->S[ve_mode][a][b];
            }
          }
        }
      }
      for (a = 0; a < VIM; a++) {
        for (b = 0; b < VIM; b++) {
          stream_grad += fv->v[a] * ves[a][b] * nv[b];
        }
      }
    } // if pd->v[pg->imtrx][POLYMER_STRESS11]
    if (DOUBLE_NONZERO(velo_sqrd)) {
      Tnt = stream_grad / velo_sqrd;
    } else {
      Tnt = 0.0;
    }
    local_post[STREAM_SHEAR_STRESS] = Tnt;
    local_lumped[STREAM_SHEAR_STRESS] = 1.;
  }

  if (STREAM_TENSION != -1 && pd->e[pg->imtrx][R_MOMENTUM1]) {
    velo_sqrd = 0.;
    stream_grad = 0.;
    curv = 0.;

    for (a = 0; a < VIM; a++) {
      for (b = 0; b < VIM; b++) {
        gamma[a][b] = fv->grad_v[a][b] + fv->grad_v[b][a];
      }
    }
    memset(vdelvdx, 0, sizeof(dbl) * DIM);
    for (a = 0; a < dim; a++) {
      velo_sqrd += fv->v[a] * fv->v[a];
      for (b = 0; b < dim; b++) {
        vdelvdx[b] += fv->v[a] * fv->grad_v[b][a];
      }
    }
    for (a = 0; a < dim; a++) {
      for (b = 0; b < dim; b++) {
        /* (2*vv-Iv):gamma */
        stream_grad += (2. * fv->v[a] * fv->v[b] - delta(a, b) * velo_sqrd) * gamma[a][b];
      }
    }
    for (a = 0; a < dim; a++) {
      for (b = 0; b < dim; b++) {
        for (c = 0; c < dim; c++) {
          curv += SQUARE(velo_sqrd * fv->grad_v[c][b] - fv->v[b] * vdelvdx[c]);
        }
      }
    }
    if (DOUBLE_NONZERO(velo_sqrd))
      curv = sqrt(curv / CUBE(velo_sqrd));
    stream_grad *= mp->viscosity;
    /*
     **  viscoelastic stress tensor
     **/
    if (pd->v[pg->imtrx][POLYMER_STRESS11]) {
      dbl log_c[DIM][DIM];
      dbl exp_s[DIM][DIM];
      dbl R1[DIM][DIM];
      dbl eig_values[DIM];
      dbl mup = 0.;
      dbl lambda = 0.;
      int ve_mode, p, r;
      memset(ves, 0, sizeof(dbl) * DIM * DIM);
      if (vn->evssModel == LOG_CONF || vn->evssModel == LOG_CONF_GRADV) {
        for (ve_mode = 0; ve_mode < vn->modes; ve_mode++) {
          for (p = 0; p < VIM; p++) {
            for (r = 0; r < VIM; r++) {
              log_c[p][r] = fv->S[ve_mode][p][r];
            }
          }
#ifdef ANALEIG_PLEASE
          analytical_exp_s(log_c, exp_s, eig_values, R1, NULL);
#else
          compute_exp_s(log_c, exp_s, eig_values, R1);
#endif
          mup = viscosity(ve[ve_mode]->gn, gamma, NULL);
          if (ve[ve_mode]->time_constModel == CONSTANT) {
            lambda = ve[ve_mode]->time_const;
          }
          for (a = 0; a < VIM; a++) {
            for (b = 0; b < VIM; b++) {
              ves[a][b] += mup / lambda * (exp_s[a][b] - (double)delta(a, b));
            }
          }
        }
      } else {
        for (a = 0; a < VIM; a++) {
          for (b = 0; b < VIM; b++) {
            for (ve_mode = 0; ve_mode < vn->modes; ve_mode++) {
              ves[a][b] += fv->S[ve_mode][a][b];
            }
          }
        }
      }
      for (a = 0; a < VIM; a++) {
        for (b = 0; b < VIM; b++) {
          stream_grad += (2. * fv->v[a] * fv->v[b] - delta(a, b) * velo_sqrd) * ves[a][b];
        }
      }
    } // if pd->v[pg->imtrx][POLYMER_STRESS11]
    /* Need to determine curvature to multiply here...*/

    if (DOUBLE_NONZERO(velo_sqrd)) {
      Ttt = curv * stream_grad / velo_sqrd;
    } else {
      Ttt = 0.0;
    }
    local_post[STREAM_TENSION] = Ttt;
    local_lumped[STREAM_TENSION] = 1.;
  }

  if (DIV_VELOCITY != -1 && pd->e[pg->imtrx][PRESSURE]) {
    Dnn = 0.;
    for (a = 0; a < dim; a++) {
      { Dnn += fv->grad_v[a][a]; }
    }
    /* If you want everything positive, uncomment this next line */
    /* Dnn = fabs(Dnn);  */
    local_post[DIV_VELOCITY] = Dnn;
    local_lumped[DIV_VELOCITY] = 1.;
  }

  if (DIV_PVELOCITY != -1 && pd->e[pg->imtrx][R_PMOMENTUM1]) {
    Dnn = 0.0;
    for (a = 0; a < dim; a++)
      Dnn += fv->grad_pv[a][a];
    local_post[DIV_PVELOCITY] = Dnn;
    local_lumped[DIV_PVELOCITY] = 1.0;
  }

  if (DIV_TOTAL != -1 && pd->e[pg->imtrx][R_PMOMENTUM1]) {
    species_particle = (int)mp->u_density[0];
    p_vol_frac = fv->c[species_particle];
    Dnn = 0.0;
    for (a = 0; a < dim; a++) {
      dbl tempf, tempp;

      tempf = (1.0 - p_vol_frac) * fv->grad_v[a][a];
      tempp = p_vol_frac * fv->grad_pv[a][a];
      Dnn += tempf + tempp;
    }
    local_post[DIV_TOTAL] = Dnn;
    local_lumped[DIV_TOTAL] = 1.0;
  }

  if (PP_Viscosity != -1 && pd->e[pg->imtrx][R_MOMENTUM1]) {
    for (a = 0; a < VIM; a++) {
      for (b = 0; b < VIM; b++) {
        gamma[a][b] = fv->grad_v[a][b] + fv->grad_v[b][a];
      }
    }

    mu = viscosity(gn, gamma, NULL);

    if (pd->v[pg->imtrx][POLYMER_STRESS11]) {
      /*  shift factor  */
      if (pd->e[pg->imtrx][TEMPERATURE]) {
        if (vn->shiftModel == CONSTANT) {
          at = vn->shift[0];
        } else if (vn->shiftModel == MODIFIED_WLF) {
          wlf_denom = vn->shift[1] + fv->T - mp->reference[TEMPERATURE];
          if (wlf_denom != 0.) {
            at = exp(vn->shift[0] * (mp->reference[TEMPERATURE] - fv->T) / wlf_denom);
          } else {
            at = 1.;
          }
        }
      } else {
        at = 1.;
      }
      for (mode = 0; mode < vn->modes; mode++) {
        /* get polymer viscosity */
        mup = viscosity(ve[mode]->gn, gamma, NULL);
        mu += at * mup;
      }
    }

    if (gn->ConstitutiveEquation == BINGHAM_MIXED) {

      dbl gamma[DIM][DIM];
      for (a = 0; a < VIM; a++) {
        for (b = 0; b < VIM; b++) {
          gamma[a][b] = fv->grad_v[a][b] + fv->grad_v[b][a];
        }
      }

      dbl gammadot;
      calc_shearrate(&gammadot, gamma, NULL, NULL);

      dbl Du_eps = sqrt(gammadot + gn->epsilon * gn->epsilon);

      mu += gn->tau_y / Du_eps;
    }

    local_post[PP_Viscosity] = mu;
    local_lumped[PP_Viscosity] = 1.;
  }

  if (PP_Viscosity != -1 &&
      (pd->e[pg->imtrx][R_LUBP] || pd->e[pg->imtrx][R_SHELL_FILMP] || pd->e[pg->imtrx][R_LUBP_2])) {
    if (gn->ConstitutiveEquation == POWER_LAW || gn->ConstitutiveEquation == BINGHAM ||
        gn->ConstitutiveEquation == HERSCHEL_BULKLEY || gn->ConstitutiveEquation == CARREAU ||
        gn->ConstitutiveEquation == CARREAU_WLF || gn->ConstitutiveEquation == BINGHAM_WLF) {
      int *n_dof = NULL;
      int dof_map[MDE];
      n_dof = (int *)array_alloc(1, MAX_VARIABLE_TYPES, sizeof(int));
      lubrication_shell_initialize(n_dof, dof_map, -1, xi, exo, 0);
      /* Calculate Shear-rate */
      if (pd->e[pg->imtrx][R_LUBP]) {
        calculate_lub_q_v(R_LUBP, time, delta_t, xi, exo);
      } else if (pd->e[pg->imtrx][R_LUBP_2]) {
        calculate_lub_q_v(R_LUBP_2, time, delta_t, xi, exo);
      } else {
        calculate_lub_q_v(R_SHELL_FILMP, time, delta_t, xi, exo);
      }
      local_post[PP_Viscosity] = LubAux->mu_star;
    } else {
      local_post[PP_Viscosity] = viscosity(gn, NULL, NULL);
    }
    local_lumped[PP_Viscosity] = 1.0;
  }

  if (PP_FlowingLiquid_Viscosity != -1 && pd->e[pg->imtrx][R_MOMENTUM1]) {
    if (mp->PorousMediaType == POROUS_BRINKMAN) {
      mu = flowing_liquid_viscosity(NULL);
    } else {
      mu = 0.0;
    }
    local_post[PP_FlowingLiquid_Viscosity] = mu;
    local_lumped[PP_FlowingLiquid_Viscosity] = 1.0;
  }

  if (PP_VolumeFractionGas != -1 && pd->e[pg->imtrx][R_MOMENTUM1]) {
    double volF = 0.0;
    computeCommonMaterialProps_gp(time);
    volF = mp->volumeFractionGas;
    local_post[PP_VolumeFractionGas] = volF;
    local_lumped[PP_VolumeFractionGas] = 1.0;
  }

  if (DENSITY != -1 && pd->e[pg->imtrx][R_MOMENTUM1]) {
    rho = density(NULL, time);
    local_post[DENSITY] = rho;
    local_lumped[DENSITY] = 1.;
  }

  if (HEAVISIDE != -1 && ls != NULL && pd->e[pg->imtrx][R_FILL]) {
    load_lsi(ls->Length_Scale);
    local_post[HEAVISIDE] = lsi->H;
    local_lumped[HEAVISIDE] = 1.;
  }

  if (RHO_DOT != -1 && pd->e[pg->imtrx][R_MOMENTUM1]) {
    double rho_dot = 0;
    if (mp->DensityModel == DENSITY_FOAM_PMDI_10) {
      double M_CO2 = mp->u_density[0];
      double rho_liq = mp->u_density[1];
      double ref_press = mp->u_density[2];
      double Rgas_const = mp->u_density[3];
      if (pd->gv[MOMENT1]) {
        double rho_gas = 0;

        if (fv->T > 0) {
          rho_gas = (ref_press * M_CO2 / (Rgas_const * fv->T));
        }

        double nu = fv->moment[1];
        double nu_dot = fv_dot->moment[1];

        double inv1 = 1 / (1 + nu);
        double inv2 = inv1 * inv1;

        double volF_dot = (nu_dot)*inv2;

        rho_dot = rho_gas * volF_dot - rho_liq * volF_dot;
      } else {
        int wCO2;
        int wH2O;
        int w;

        wCO2 = -1;
        wH2O = -1;
        for (w = 0; w < pd->Num_Species; w++) {
          switch (mp->SpeciesSourceModel[w]) {
          case FOAM_PMDI_10_CO2:
            wCO2 = w;
            break;
          case FOAM_PMDI_10_H2O:
            wH2O = w;
            break;
          default:
            break;
          }
        }

        if (wCO2 == -1) {
          GOMA_EH(GOMA_ERROR, "Expected a Species Source of FOAM_PMDI_10_CO2");
        } else if (wH2O == -1) {
          GOMA_EH(GOMA_ERROR, "Expected a Species Source of FOAM_PMDI_10_H2O");
        }

        double rho_gas = 0;

        if (fv->T > 0) {
          rho_gas = (ref_press * M_CO2 / (Rgas_const * fv->T));
        }

        double nu = 0;
        double nu_dot = 0;

        if (fv->T > 0) {
          nu = M_CO2 * fv->c[wCO2] / rho_gas;
          nu_dot = M_CO2 * fv_dot->c[wCO2] / rho_gas;
        } else {
          nu = 0;
          nu_dot = 0;
        }

        double inv1 = 1 / (1 + nu);
        double inv2 = inv1 * inv1;

        double volF_dot = (nu_dot)*inv2;

        rho_dot = rho_gas * volF_dot - rho_liq * volF_dot;
      }
    }
    local_post[RHO_DOT] = rho_dot;
    local_lumped[RHO_DOT] = 1.;
  }

  if (POLYMER_VISCOSITY != -1 && pd->e[pg->imtrx][R_STRESS11]) {
    mode = 0;
    double mup = viscosity(ve[mode]->gn, gamma, NULL);
    local_post[POLYMER_VISCOSITY] = mup;
    local_lumped[POLYMER_VISCOSITY] = 1.;
  }

  if (POLYMER_TIME_CONST != -1 && pd->e[pg->imtrx][R_STRESS11]) {
    mode = 0;
    double lambda = 0;
    double mup = viscosity(ve[mode]->gn, gamma, NULL);
    if (ve[mode]->time_constModel == CONSTANT) {
      lambda = ve[mode]->time_const;
    } else if (ve[mode]->time_constModel == CARREAU || ve[mode]->time_constModel == POWER_LAW) {
      lambda = mup / ve[mode]->time_const;
    } else if (ls != NULL && ve[mode]->time_constModel == VE_LEVEL_SET) {
      double pos_lambda = ve[mode]->pos_ls.time_const;
      double neg_lambda = ve[mode]->time_const;
      double width = ls->Length_Scale;
      err = level_set_property(neg_lambda, pos_lambda, width, &lambda, NULL);
      GOMA_EH(err, "level_set_property() failed for polymer time constant.");
    }
    local_post[POLYMER_TIME_CONST] = lambda;
    local_lumped[POLYMER_TIME_CONST] = 1.;
  }

  if (PTT_XI != -1 && pd->e[pg->imtrx][R_STRESS11]) {
    mode = 0;
    double xi = 0;
    double mup = viscosity(ve[mode]->gn, gamma, NULL);
    if (ve[mode]->xiModel == CONSTANT) {
      xi = ve[mode]->xi;
    } else if (ve[mode]->xiModel == CARREAU || ve[mode]->xiModel == POWER_LAW) {
      xi = mup / ve[mode]->xi;
    } else if (ls != NULL && ve[mode]->xiModel == VE_LEVEL_SET) {
      double pos_xi = ve[mode]->pos_ls.xi;
      double neg_xi = ve[mode]->xi;
      double width = ls->Length_Scale;
      err = level_set_property(neg_xi, pos_xi, width, &xi, NULL);
      GOMA_EH(err, "level_set_property() failed for polymer time constant.");
    }
    local_post[PTT_XI] = xi;
    local_lumped[PTT_XI] = 1.;
  }

  if (PTT_EPSILON != -1 && pd->e[pg->imtrx][R_STRESS11]) {
    mode = 0;
    double eps = 0;
    double mup = viscosity(ve[mode]->gn, gamma, NULL);
    if (ve[mode]->epsModel == CONSTANT) {
      eps = ve[mode]->eps;
    } else if (ve[mode]->epsModel == CARREAU || ve[mode]->epsModel == POWER_LAW) {
      eps = mup / ve[mode]->eps;
    } else if (ls != NULL && ve[mode]->epsModel == VE_LEVEL_SET) {
      double pos_eps = ve[mode]->pos_ls.eps;
      double neg_eps = ve[mode]->eps;
      double width = ls->Length_Scale;
      err = level_set_property(neg_eps, pos_eps, width, &eps, NULL);
      GOMA_EH(err, "level_set_property() failed for polymer time constant.");
    }
    local_post[PTT_EPSILON] = eps;
    local_lumped[PTT_EPSILON] = 1.;
  }

  if (MOBILITY_PARAMETER != -1 && pd->e[pg->imtrx][R_STRESS11]) {
    mode = 0;
    double alpha = 0;
    double mup = viscosity(ve[mode]->gn, gamma, NULL);
    if (ve[mode]->alphaModel == CONSTANT) {
      alpha = ve[mode]->alpha;
    } else if (ve[mode]->alphaModel == CARREAU || ve[mode]->alphaModel == POWER_LAW) {
      alpha = mup / ve[mode]->alpha;
    } else if (ls != NULL && ve[mode]->alphaModel == VE_LEVEL_SET) {
      double pos_alpha = ve[mode]->pos_ls.alpha;
      double neg_alpha = ve[mode]->alpha;
      double width = ls->Length_Scale;
      err = level_set_property(neg_alpha, pos_alpha, width, &alpha, NULL);
      GOMA_EH(err, "level_set_property() failed for polymer time constant.");
    }
    local_post[MOBILITY_PARAMETER] = alpha;
    local_lumped[MOBILITY_PARAMETER] = 1.;
  }
  if (MEAN_SHEAR != -1) {
    if (pd->e[pg->imtrx][R_MOMENTUM1]) {
      double gammadot, gamma[DIM][DIM];
      for (a = 0; a < VIM; a++) {
        for (b = 0; b < VIM; b++) {
          gamma[a][b] = fv->grad_v[a][b] + fv->grad_v[b][a];
        }
      }
      /* find second invariant of strain-rate */
      calc_shearrate(&gammadot, gamma, NULL, NULL);

      local_post[MEAN_SHEAR] = gammadot;
      local_lumped[MEAN_SHEAR] = 1.;
    }
    if (pd->e[pg->imtrx][R_LUBP] || pd->e[pg->imtrx][R_SHELL_FILMP] || pd->e[pg->imtrx][R_LUBP_2]) {
      int *n_dof = NULL;
      int dof_map[MDE];
      n_dof = (int *)array_alloc(1, MAX_VARIABLE_TYPES, sizeof(int));
      lubrication_shell_initialize(n_dof, dof_map, -1, xi, exo, 0);
      /* Calculate Shear-rate */
      if (pd->e[pg->imtrx][R_LUBP]) {
        calculate_lub_q_v(R_LUBP, time, delta_t, xi, exo);
      } else if (pd->e[pg->imtrx][R_LUBP_2]) {
        calculate_lub_q_v(R_LUBP_2, time, delta_t, xi, exo);
      } else {
        calculate_lub_q_v(R_SHELL_FILMP, time, delta_t, xi, exo);
      }
      local_post[PP_Viscosity] = LubAux->srate;
      local_lumped[MEAN_SHEAR] = 1.;
    }
  }

  if (GIES_CRIT != -1 && pd->e[pg->imtrx][R_MOMENTUM1]) {
    double gammadot, gamma[DIM][DIM];
    double vorticity, omega[DIM][DIM];
    for (a = 0; a < VIM; a++) {
      for (b = 0; b < VIM; b++) {
        gamma[a][b] = fv->grad_v[a][b] + fv->grad_v[b][a];
        omega[a][b] = fv->grad_v[a][b] - fv->grad_v[b][a];
      }
    }
    /* find second invariant of strain-rate */
    calc_shearrate(&gammadot, gamma, NULL, NULL);
    calc_shearrate(&vorticity, omega, NULL, NULL);
    local_post[GIES_CRIT] = (gammadot - vorticity) / (DBL_SMALL + gammadot + vorticity);
    local_lumped[GIES_CRIT] = 1.;
  }

  if (Q_FCN != -1 && pd->e[pg->imtrx][R_MOMENTUM1]) {
    double gammadot, del_v[DIM][DIM];
    for (a = 0; a < VIM; a++) {
      for (b = 0; b < VIM; b++) {
        del_v[a][b] = fv->grad_v[a][b];
      }
    }
    /* find second invariant of velocity gradient tensor */
    calc_shearrate(&gammadot, del_v, NULL, NULL);

    local_post[Q_FCN] = gammadot;
    local_lumped[Q_FCN] = 1.;
  }

  if (PSPG_PP != -1 && pd->e[pg->imtrx][R_MOMENTUM1]) {
    dbl pspg[DIM];
    calc_pspg(pspg, NULL, time, theta, delta_t, pg_data);
    for (int i = 0; i < VIM; i++) {
      local_post[PSPG_PP + i] = pspg[i];
      local_lumped[PSPG_PP + i] = 1.;
    }
  }

  if (VELO_SPEED != -1 && pd->e[pg->imtrx][R_MOMENTUM1]) {
    velo_sqrd = 0.;
    for (a = 0; a < VIM; a++) {
      velo_sqrd += SQUARE(fv->v[a]);
    }
    local_post[VELO_SPEED] = sqrt(velo_sqrd);
    local_lumped[VELO_SPEED] = 1.;
  }

  if (PRESSURE_CONT != -1 && (pd->v[pg->imtrx][PRESSURE] || pd->v[pg->imtrx][TFMP_PRES]) &&
      (pd->e[pg->imtrx][R_MOMENTUM1] ||
       (pd->MeshMotion == LAGRANGIAN || pd->MeshMotion == DYNAMIC_LAGRANGIAN) ||
       (pd->MeshMotion == TOTAL_ALE))) {
    if (pd->v[pg->imtrx][PRESSURE]) {
      local_post[PRESSURE_CONT] = fv->P;
    } else if (pd->v[pg->imtrx][TFMP_PRES]) {
      local_post[PRESSURE_CONT] = fv->tfmp_pres;
    }
    local_lumped[PRESSURE_CONT] = 1.;
  }

  if (SH_DIV_S_V_CONT != -1 && pd->v[pg->imtrx][SHELL_SURF_DIV_V]) {
    local_post[SH_DIV_S_V_CONT] = fv->div_s_v;
    local_lumped[SH_DIV_S_V_CONT] = 1.;
  }

  if (SH_CURV_CONT != -1 && pd->v[pg->imtrx][SHELL_SURF_CURV]) {
    local_post[SH_CURV_CONT] = fv->curv;
    local_lumped[SH_CURV_CONT] = 1.;
  }

  if (FILL_CONT != -1 && pd->v[pg->imtrx][FILL]) {
    local_post[FILL_CONT] = fv->F;
    local_lumped[FILL_CONT] = 1.;
  }

  if (CONC_CONT != -1 && pd->v[pg->imtrx][MASS_FRACTION]) {
    for (w = 0; w < pd->Num_Species_Eqn; w++) {
      local_post[CONC_CONT + w] = fv->c[w];
      local_lumped[CONC_CONT + w] = 1.;
    }
  }

  if (EM_CONTOURS != -1 && pd->v[pg->imtrx][EM_E1_REAL]) {
    index = 0;
    for (b = 0; b < DIM; b++) {

      if (pd->v[pg->imtrx][EM_E1_REAL]) {
        local_post[EM_CONTOURS + index] = fv->em_er[b];
        local_lumped[EM_CONTOURS + index] = 1.;
        index++;
      }
    }

    for (b = 0; b < DIM; b++) {

      if (pd->v[pg->imtrx][EM_E1_IMAG]) {
        local_post[EM_CONTOURS + index] = fv->em_ei[b];
        local_lumped[EM_CONTOURS + index] = 1.;
        index++;
      }
    }
  }

  if (TOTAL_EM_CONTOURS != -1 && pd->v[pg->imtrx][EM_E1_REAL]) {
    index = 0;
    const double c0 = 1.0 / (sqrt(upd->Free_Space_Permittivity * upd->Free_Space_Permeability));
    dbl x = fv->x[0];
    dbl y = fv->x[1];
    dbl z = fv->x[2];
    dbl freq = upd->EM_Frequency;
    dbl lambda0 = c0 / freq;
    dbl k0 = 2 * M_PI / lambda0;
    complex double wave[3] = {0};
    complex double curl_wave[3];
    if (mp->PermittivityModel != RADIAL_PML) {
      incident_wave(x, y, z, k0, wave, curl_wave);
    }

    for (b = 0; b < DIM; b++) {

      if (pd->v[pg->imtrx][EM_E1_REAL]) {
        local_post[TOTAL_EM_CONTOURS + index] = creal(wave[b]);
        local_lumped[TOTAL_EM_CONTOURS + index] = 1.;
        index++;
      }
    }

    for (b = 0; b < DIM; b++) {

      if (pd->v[pg->imtrx][EM_E1_IMAG]) {
        local_post[TOTAL_EM_CONTOURS + index] = cimag(wave[b]);
        local_lumped[TOTAL_EM_CONTOURS + index] = 1.;
        index++;
      }
    }
  }

  if (SCATTERED_EM_CONTOURS != -1 && pd->v[pg->imtrx][EM_E1_REAL]) {
    index = 0;
    const double c0 = 1.0 / (sqrt(upd->Free_Space_Permittivity * upd->Free_Space_Permeability));
    dbl x = fv->x[0];
    dbl y = fv->x[1];
    dbl z = fv->x[2];
    dbl freq = upd->EM_Frequency;
    dbl lambda0 = c0 / freq;
    dbl k0 = 2 * M_PI / lambda0;
    complex double wave[3] = {0};
    complex double curl_wave[3];
    if (mp->PermittivityModel != RADIAL_PML) {
      incident_wave(x, y, z, k0, wave, curl_wave);
    }

    for (b = 0; b < DIM; b++) {

      if (pd->v[pg->imtrx][EM_E1_REAL]) {
        local_post[SCATTERED_EM_CONTOURS + index] = fv->em_er[b] - creal(wave[b]);
        local_lumped[SCATTERED_EM_CONTOURS + index] = 1.;
        index++;
      }
    }

    for (b = 0; b < DIM; b++) {

      if (pd->v[pg->imtrx][EM_E1_IMAG]) {
        local_post[SCATTERED_EM_CONTOURS + index] = fv->em_ei[b] - cimag(wave[b]);
        local_lumped[SCATTERED_EM_CONTOURS + index] = 1.;
        index++;
      }
    }
  }

  if (STRESS_CONT != -1 && pd->v[pg->imtrx][POLYMER_STRESS11]) {
    index = 0;
    for (mode = 0; mode < vn->modes; mode++) {
      if (pd->v[pg->imtrx][v_s[mode][0][0]]) {
        for (a = 0; a < VIM; a++) {
          for (b = 0; b < VIM; b++) {
            /* since the stress tensor is symmetric,
               only assemble the upper half */
            if (a <= b) {
              if (pd->v[pg->imtrx][v_s[mode][a][b]]) {
                local_post[STRESS_CONT + index] = fv->S[mode][a][b];
                local_lumped[STRESS_CONT + index] = 1.;
                index++;
              }
            }
          }
        }
      }
    }

    /* Get total stress tensor for multi mode calculations */
    if (vn->modes > 1) {
      for (a = 0; a < VIM; a++) {
        for (b = 0; b < VIM; b++) {
          /* since the stress tensor is symmetric, only assemble the
             upper half */
          if (a <= b) {
            ts = 0.;
            for (mode = 0; mode < vn->modes; mode++) {
              if (pd->v[pg->imtrx][v_s[mode][a][b]]) {
                ts += fv->S[mode][a][b];
              }
            }
            local_post[STRESS_CONT + index] = ts;
            local_lumped[STRESS_CONT + index] = 1.;
            index++;
          }
        }
      }
    }
  }

  if (MOMENT_SOURCES != -1 && pd->v[pg->imtrx][MOMENT0]) {
    double msource[MAX_MOMENTS];
    MOMENT_SOURCE_DEPENDENCE_STRUCT *d_msource;
    d_msource = calloc(sizeof(MOMENT_SOURCE_DEPENDENCE_STRUCT), 1);
    moment_source(msource, d_msource);

    for (int mom = 0; mom < MAX_MOMENTS; mom++) {
      if (pd->gv[MOMENT0 + mom]) {
        local_post[MOMENT_SOURCES + mom] = msource[mom];
        local_lumped[MOMENT_SOURCES + mom] = 1.;
      }
    }
    free(d_msource);
  }

  if (YZBETA != -1 && pd->v[pg->imtrx][MASS_FRACTION]) {
    struct Petrov_Galerkin_Data pg_data;
    h_elem_siz(pg_data.hsquared, pg_data.hhv, pg_data.dhv_dxnode, pd->e[pg->imtrx][R_MESH1]);
    struct Species_Conservation_Terms s_terms;
    zero_structure(&s_terms, sizeof(struct Species_Conservation_Terms), 1);
    err = get_continuous_species_terms(&s_terms, time, theta, delta_t, pg_data.hsquared);

    dbl h_elem = 0;
    for (int a = 0; a < ei[pg->imtrx]->ielem_dim; a++) {
      h_elem += pg_data.hsquared[a];
    }

    /* This is the size of the element */
    h_elem = sqrt(h_elem / ((double)ei[pg->imtrx]->ielem_dim));

    for (int w = 0; w < pd->Num_Species_Eqn; w++) {
      dbl strong_residual = 0;
      strong_residual = fv_dot->c[w];
      for (int p = 0; p < VIM; p++) {
        strong_residual += fv->v[p] * fv->grad_c[w][p];
      }
      // strong_residual -= s_terms.MassSource[w];
      dbl h_elem = 0;
      for (int a = 0; a < ei[pg->imtrx]->ielem_dim; a++) {
        h_elem += pg_data.hsquared[a];
      }
      /* This is the size of the element */
      h_elem = sqrt(h_elem / ((double)ei[pg->imtrx]->ielem_dim));

      dbl inner = 0;
      for (int i = 0; i < dim; i++) {
        inner += fv->grad_c[w][i] * fv->grad_c[w][i];
      }

      //        dbl yzbeta = 0;

      dbl inv_sqrt_inner = (1 / sqrt(inner + 1e-12));
      dbl dc1 = fabs(strong_residual) * inv_sqrt_inner * h_elem * 0.5;
      dbl dc2 = fabs(strong_residual) * h_elem * h_elem * 0.25;
      //        yzbeta = 0.5*(dc1 + dc2);

      local_post[YZBETA + w] =
          0.5 *
          (dc1 +
           dc2); // fmin(dc1,
                 // supg_terms.supg_tau);//yzbeta;//supg_terms.supg_tau;//fmin(supg_terms.supg_tau,
                 // 0.5*(dc1 + dc2));
      local_lumped[YZBETA + w] = 1.;
    }
  }

  if (FIRST_INVAR_STRAIN != -1 && pd->e[pg->imtrx][R_MESH1]) {
    TrE = 0.;
    if (pd->CoordinateSystem == CYLINDRICAL) {
      TrE = fv->strain[0][0] + 2 * fv->strain[1][1];
    } else {
      for (a = 0; a < dim; a++) {
        TrE += fv->strain[a][a];
      }
    }
    /*     TrE = fv->div_d; */
    local_post[FIRST_INVAR_STRAIN] = TrE;
    local_lumped[FIRST_INVAR_STRAIN] = 1.;
  }

  if (SEC_INVAR_STRAIN != -1 && pd->e[pg->imtrx][R_MESH1]) {
    E_E = 0.;
    /* find second invariant of strain */
    for (a = 0; a < dim; a++) {
      for (b = 0; b < dim; b++) {
        E_E += 0.5 * (fv->strain[a][b] * fv->strain[a][b] - fv->strain[a][a] * fv->strain[b][b]);
      }
    }
    local_post[SEC_INVAR_STRAIN] = E_E;
    local_lumped[SEC_INVAR_STRAIN] = 1.;
  }

  if (THIRD_INVAR_STRAIN != -1 && pd->e[pg->imtrx][R_MESH1]) {
    /* this is actually the volume change - third invarient of the Deformation Gradient! */
    local_post[THIRD_INVAR_STRAIN] = fv->volume_change;
    local_lumped[THIRD_INVAR_STRAIN] = 1.;
  }

  if (FIRST_STRAINRATE_INVAR != -1 && pd->e[pg->imtrx][R_MESH1]) {
    dbl INV, d_INV_dT[DIM][DIM];
    INV = calc_tensor_invariant(fv_dot->strain, d_INV_dT, 1);
    local_post[FIRST_STRAINRATE_INVAR] = INV;
    local_lumped[FIRST_STRAINRATE_INVAR] = 1.;
  }

  if (SEC_STRAINRATE_INVAR != -1 && pd->e[pg->imtrx][R_MESH1]) {
    dbl INV, d_INV_dT[DIM][DIM];
    INV = calc_tensor_invariant(fv_dot->strain, d_INV_dT, 2);
    local_post[SEC_STRAINRATE_INVAR] = INV;
    local_lumped[SEC_STRAINRATE_INVAR] = 1.;
  }

  if (THIRD_STRAINRATE_INVAR != -1 && pd->e[pg->imtrx][R_MESH1]) {
    dbl INV, d_INV_dT[DIM][DIM];
    INV = calc_tensor_invariant(fv_dot->strain, d_INV_dT, 3);
    local_post[THIRD_STRAINRATE_INVAR] = INV;
    local_lumped[THIRD_STRAINRATE_INVAR] = 1.;
  }

  if (DIELECTROPHORETIC_FIELD != -1 && pd->e[pg->imtrx][R_ENORM]) {
    if (Particle_Model_Data[1] <= 0.0 || Particle_Model_Data[2] <= 0.0 ||
        Particle_Model_Data[3] <= 0.0 || Particle_Model_Data[4] <= 0.0 ||
        Particle_Model_Data[5] <= 0.0 || Particle_Model_Data[6] <= 0.0 || Particle_Density <= 0.0 ||
        Particle_Radius <= 0.0)
      GOMA_WH(-1, "Not performing Dielectrophoretic Force Vectors output becuase of inconsistent "
                  "Partile_Model_Data.");
    coeff_a = Particle_Model_Data[1] - Particle_Model_Data[2];
    coeff_b = (Particle_Model_Data[4] - Particle_Model_Data[3]) / Particle_Model_Data[5];
    coeff_c = Particle_Model_Data[1] + 2.0 * Particle_Model_Data[2];
    coeff_d = (Particle_Model_Data[3] + 2.0 * Particle_Model_Data[4]) / Particle_Model_Data[5];
    CM_fact = (coeff_a * coeff_c - coeff_b * coeff_d) / (coeff_c * coeff_c + coeff_d * coeff_d);
    dielectrophoretic_force_coeff = 4.0 * M_PIE * Particle_Radius * Particle_Radius *
                                    Particle_Radius * Particle_Model_Data[2] * CM_fact *
                                    Particle_Model_Data[6] * Particle_Model_Data[6];
    for (a = 0; a < dim; a++)
      dfvector[a] = dielectrophoretic_force_coeff * fv->Enorm * fv->grad_Enorm[a];
    for (a = 0; a < dim; a++) {
      local_post[DIELECTROPHORETIC_FIELD + a] = dfvector[a];
      local_lumped[DIELECTROPHORETIC_FIELD + a] = 1.0;
    }
  }

  if (DIELECTROPHORETIC_FIELD_NORM != -1 && pd->e[pg->imtrx][R_ENORM]) {
    if (Particle_Model_Data[1] <= 0.0 || Particle_Model_Data[2] <= 0.0 ||
        Particle_Model_Data[3] <= 0.0 || Particle_Model_Data[4] <= 0.0 ||
        Particle_Model_Data[5] <= 0.0 || Particle_Model_Data[6] <= 0.0 || Particle_Density <= 0.0 ||
        Particle_Radius <= 0.0)
      GOMA_WH(-1, "Not performing Dielectrophoretic Force Vectors output becuase of inconsistent "
                  "Partile_Model_Data.");
    coeff_a = Particle_Model_Data[1] - Particle_Model_Data[2];
    coeff_b = (Particle_Model_Data[4] - Particle_Model_Data[3]) / Particle_Model_Data[5];
    coeff_c = Particle_Model_Data[1] + 2.0 * Particle_Model_Data[2];
    coeff_d = (Particle_Model_Data[3] + 2.0 * Particle_Model_Data[4]) / Particle_Model_Data[5];
    CM_fact = (coeff_a * coeff_c - coeff_b * coeff_d) / (coeff_c * coeff_c + coeff_d * coeff_d);
    dielectrophoretic_force_coeff = 4.0 * M_PIE * Particle_Radius * Particle_Radius *
                                    Particle_Radius * Particle_Model_Data[2] * CM_fact *
                                    Particle_Model_Data[6] * Particle_Model_Data[6];
    for (a = 0; a < dim; a++)
      dfvector[a] = dielectrophoretic_force_coeff * fv->Enorm * fv->grad_Enorm[a];
    dfnorm = 0.0;
    for (a = 0; a < dim; a++)
      dfnorm += dfvector[a] * dfvector[a];
    dfnorm = MAX(0.0, dfnorm);
    dfnorm = sqrt(dfnorm);
    local_post[DIELECTROPHORETIC_FIELD_NORM] = dfnorm;
    local_lumped[DIELECTROPHORETIC_FIELD_NORM] = 1.0;
  }

  if (ENORMSQ_FIELD != -1 && pd->e[pg->imtrx][R_ENORM])
    for (a = 0; a < dim; a++) {
      local_post[ENORMSQ_FIELD + a] = 2.0 * fv->Enorm * fv->grad_Enorm[a];
      local_lumped[ENORMSQ_FIELD + a] = 1.0;
    }

  if (ENORMSQ_FIELD_NORM != -1 && pd->e[pg->imtrx][R_ENORM]) {
    dfnorm = 0.0;
    for (a = 0; a < dim; a++)
      dfnorm += (2.0 * fv->Enorm * fv->grad_Enorm[a]) * (2.0 * fv->Enorm * fv->grad_Enorm[a]);
    dfnorm = MAX(0.0, dfnorm);
    dfnorm = sqrt(dfnorm);
    local_post[ENORMSQ_FIELD_NORM] = dfnorm;
    local_lumped[ENORMSQ_FIELD_NORM] = 1.0;
  }

  if (DIFFUSION_VECTORS != -1 && pd->e[pg->imtrx][R_MASS]) {
    if (mp->PorousMediaType == CONTINUOUS) {
      if (cr->MassFluxModel == FICKIAN) {
        for (w = 0; w < pd->Num_Species_Eqn; w++) {
          for (a = 0; a < dim; a++) {
            n[w][a] = -mp->diffusivity[w] * fv->grad_c[w][a];
          }
        }
      } else if (cr->MassFluxModel == HYDRODYNAMIC || cr->MassFluxModel == HYDRODYNAMIC_QTENSOR ||
                 cr->MassFluxModel == HYDRODYNAMIC_QTENSOR_OLD) {
        /* Don't bother with the element size for the shock capturing
           term when we are post processing */
        for (a = 0; a < dim; a++)
          hs[a] = 0.;

        for (w = 0; w < pd->Num_Species_Eqn; w++) {
          hydro_flux(&s_terms, w, theta, delta_t, hs);
          for (a = 0; a < dim; a++) {
            n[w][a] = s_terms.diff_flux[w][a];
          }
        }
      } else {
        for (w = 0; w < pd->Num_Species_Eqn; w++) {
          for (a = 0; a < dim; a++) {
            n[w][a] = 0.0;
          }
        }
      }
    }
    for (w = 0; w < pd->Num_Species_Eqn; w++) {
      for (a = 0; a < dim; a++) {
        local_post[DIFFUSION_VECTORS + w * dim + a] = n[w][a];
        local_lumped[DIFFUSION_VECTORS + w * dim + a] = 1.0;
      }
    }
  }

  if (DIFFUSION_VECTORS_POR_LIQ_GPHASE != -1 && pd->e[pg->imtrx][POR_LIQ_PRES] &&
      mp->PorousMediaType != POROUS_SATURATED) {
    for (a = 0; a < VIM; a++) {
      /*
       * First do the diffusive flux of solvent in the gas phase
       * UMR
       */
      n[0][a] =
          mp->diffusivity[0] * (pmv->d_gas_density_solvents[0][POR_LIQ_PRES] * fv->grad_p_liq[a] +
                                pmv->d_gas_density_solvents[0][POR_GAS_PRES] * fv->grad_p_gas[a]);
    }
    for (a = 0; a < dim; a++) {
      local_post[DIFFUSION_VECTORS_POR_LIQ_GPHASE + a] = n[0][a];
      local_lumped[DIFFUSION_VECTORS_POR_LIQ_GPHASE + a] = 1.0;
    }
  }

  if (DIFFUSION_VECTORS_POR_AIR_GPHASE != -1 && pd->e[pg->imtrx][POR_LIQ_PRES] &&
      mp->PorousMediaType != POROUS_SATURATED) {
    for (a = 0; a < VIM; a++) {
      /*
       * First do the diffusive flux of air in the gas phase
       * UMR
       */
      n[1][a] =
          mp->diffusivity[1] * (pmv->d_gas_density_solvents[1][POR_LIQ_PRES] * fv->grad_p_liq[a] +
                                pmv->d_gas_density_solvents[1][POR_GAS_PRES] * fv->grad_p_gas[a]);
    }
    for (a = 0; a < dim; a++) {
      local_post[DIFFUSION_VECTORS_POR_AIR_GPHASE + a] = n[1][a];
      local_lumped[DIFFUSION_VECTORS_POR_AIR_GPHASE + a] = 1.0;
    }
  }

  if (CONDUCTION_VECTORS != -1 && pd->e[pg->imtrx][R_ENERGY]) {
    if (cr->HeatFluxModel == CR_HF_FOURIER_0) {
      if (mp->ConductivityModel == USER) {
        err = usr_thermal_conductivity(mp->u_thermal_conductivity, time);
      }

      for (a = 0; a < dim; a++) {
        qc[a] = -mp->thermal_conductivity * fv->grad_T[a];
      }
    } else if (cr->HeatFluxModel == CR_HF_USER) {
      double dq_gradT[DIM][DIM], dq_dX[DIM][DIM];
      usr_heat_flux(fv->grad_T, qc, dq_gradT, dq_dX, time);
      printf("untested\n");
      exit(-1);
    } else {
      for (a = 0; a < dim; a++) {
        qc[a] = 0.;
      }
    }
    for (a = 0; a < dim; a++) {
      local_post[CONDUCTION_VECTORS + a] = qc[a];
      local_lumped[CONDUCTION_VECTORS + a] = 1.;
    }
  }

  if (ORIENTATION_VECTORS != -1) {
    /* parameters are sent through USER_POST input line for 3D printing orientation ...  */
    if (pd->e[pg->imtrx][R_ENERGY] && USER_POST != -1 && len_u_post_proc > 2) {
      double q_mag = 0., sign = 0.;
      double radius = u_post_proc[1]; // overlap = u_post_proc[2];
      int dir = 0;
      if (cr->HeatFluxModel == CR_HF_FOURIER_0) {
        if (mp->ConductivityModel == USER) {
          err = usr_thermal_conductivity(mp->u_thermal_conductivity, time);
        }
        for (a = 0; a < dim; a++) {
          qc[a] = -mp->thermal_conductivity * fv->grad_T[a];
          q_mag += SQUARE(qc[a]);
        }
      } else {
        for (a = 0; a < dim; a++) {
          qc[a] = 0.;
        }
      }
      if (fv->x0[2] <= radius) {
        dir = 0;
      } else {
        // double tmp = modf((fv->x0[2] - radius) / (2 * radius - overlap), &sign);
        sign += 1.;
        dir = (int)sign % 2;
      }
      sign = 2 * (double)dir - 1.;
      for (a = 0; a < dim; a++) {
        local_post[ORIENTATION_VECTORS + a] = sign * qc[a] / sqrt(q_mag);
        local_lumped[ORIENTATION_VECTORS + a] = 1.;
      }
    } else if (pd->e[pg->imtrx][R_MOMENTUM1]) {
      double velo_sqrd = 0.;
      for (a = 0; a < VIM; a++) {
        velo_sqrd += SQUARE(fv->v[a]);
      }
      for (a = 0; a < VIM; a++) {
        local_post[ORIENTATION_VECTORS + a] = fv->v[a] / sqrt(velo_sqrd);
        local_lumped[ORIENTATION_VECTORS + a] = 1.;
      }
    } else {
      GOMA_EH(GOMA_ERROR, "Orientation field not defined!");
    }
  }
  if (SHELL_NORMALS != -1 && (pd->e[pg->imtrx][R_SHELL_ANGLE1] || pd->e[pg->imtrx][R_LUBP])) {
    double sh_n[DIM];
    for (a = 0; a < DIM; a++) {
      sh_n[a] = 0;
    }
    if (pd->e[pg->imtrx][R_SHELL_ANGLE1]) {
      if (dim == 2) {
        sh_n[0] = cos(fv->sh_ang[0]);
        sh_n[1] = sin(fv->sh_ang[0]);
      } else {
        GOMA_EH(GOMA_ERROR, "Not hard at all to implement SHELL_NORMALS for 3D, so just do it!");
      }
    } else if (pd->e[pg->imtrx][R_LUBP]) {
      int *n_dof = NULL;
      int dof_map[MDE];
      dbl wt = fv->wt;
      n_dof = (int *)array_alloc(1, MAX_VARIABLE_TYPES, sizeof(int));
      lubrication_shell_initialize(n_dof, dof_map, -1, xi, exo, 0);
      for (a = 0; a < dim; a++) {
        sh_n[a] = fv->snormal[a];
      }
      fv->wt = wt;
      safe_free((void *)n_dof);
    } else {
      GOMA_EH(GOMA_ERROR, "Not sure how I got here.");
    }
    for (a = 0; a < dim; a++) {
      local_post[SHELL_NORMALS + a] = sh_n[a];
      local_lumped[SHELL_NORMALS + a] = 1.;
    }
  }

  /* calculate mesh stress here !!*/
  if (STRESS_TENSOR != -1 && pd->e[pg->imtrx][R_MESH1]) {
    /*
     * Total mesh stress tensor...
     */
    err = mesh_stress_tensor(TT, dTT_dx, dTT_dp, dTT_dc, dTT_dp_liq, dTT_dp_gas, dTT_dporosity,
                             dTT_dsink_mass, dTT_dT, dTT_dmax_strain, dTT_dcur_strain, elc->lame_mu,
                             elc->lame_lambda, delta_t, ielem, ip, ip_total);

    /* For LINEAR ELASTICITY */
    if (cr->MeshFluxModel == LINEAR) {
      if (dim == 2) {
        TT[2][2] = 1.;
        TT[1][2] = 0.;
        TT[0][2] = 0.;
      }
    }
    /*  For Hookian Elasticity and shrinkage */
    else {
      if (dim == 2) {
        if (cr->MeshMotion == ARBITRARY) {
          TT[2][2] = (1. - fv->volume_change) * elc->lame_mu;
        } else {
          if (cr->MeshFluxModel == NONLINEAR || cr->MeshFluxModel == HOOKEAN_PSTRAIN ||
              cr->MeshFluxModel == INCOMP_PSTRAIN)
            TT[2][2] = (1. - pow(fv->volume_change, 2. / 3.)) * elc->lame_mu - fv->P;
          /*              if (cr->MeshFluxModel == INCOMP_PSTRESS) */
          else
            TT[2][2] = 0.;
        }
        TT[1][2] = 0.;
        TT[0][2] = 0.;
      }
    }

    for (j = 0; j < DIM; j++) {
      local_post[STRESS_TENSOR + j] = TT[j][j];
      local_lumped[STRESS_TENSOR + j] = 1.;
    }
    local_post[STRESS_TENSOR + 3] = TT[0][1];
    local_lumped[STRESS_TENSOR + 3] = 1.;
    if (dim == 3) /*Note the T_theta_theta term for
                    Axisymm is picked up in previous loop */
    {
      local_post[STRESS_TENSOR + 4] = TT[0][2];
      local_post[STRESS_TENSOR + 5] = TT[1][2];
      local_lumped[STRESS_TENSOR + 4] = 1.;
      local_lumped[STRESS_TENSOR + 5] = 1.;
    }

  } /* end of STRESS_TENSOR */

  /* calculate mesh strain here !!*/
  if (STRAIN_TENSOR != -1 && pd->e[pg->imtrx][R_MESH1]) {
    for (i = 0; i < dim; i++) {
      for (j = 0; j < dim; j++) {
        EE[i][j] = fv->strain[i][j];
      }
    }

    if (cr->MeshFluxModel != INCOMP_PSTRESS && dim <= 2) {
      EE[2][2] = 0.;
    } else if (dim <= 2) {
      EE[2][2] = pow((fv->deform_grad[0][0] * fv->deform_grad[1][1] -
                      fv->deform_grad[1][0] * fv->deform_grad[0][1]),
                     0.5) /
                 pow((1. - fv->P / elc->lame_mu), 3. / 4.);
    }
    for (j = 0; j < 3; j++) {
      local_post[STRAIN_TENSOR + j] = EE[j][j];
      local_lumped[STRAIN_TENSOR + j] = 1.;
    }
    local_post[STRAIN_TENSOR + 3] = EE[0][1];
    local_lumped[STRAIN_TENSOR + 3] = 1.;
    if (dim > 2) {
      local_post[STRAIN_TENSOR + 4] = EE[0][2];
      local_post[STRAIN_TENSOR + 5] = EE[1][2];
      local_lumped[STRAIN_TENSOR + 4] = 1.;
      local_lumped[STRAIN_TENSOR + 5] = 1.;
    }
  }
  /* end of STRAIN_TENSOR */

  /* calculate EVP def grad here !!*/
  if (evpl->ConstitutiveEquation == EVP_HYPER && EVP_DEF_GRAD_TENSOR != -1 &&
      pd->e[pg->imtrx][R_MESH1]) {
    /*even though I changed this to VIM, I noticed that FVPs
      are not transmitted to restart file.....PRS 6/7/2002 */
    for (i = 0; i < VIM; i++) {
      for (j = 0; j < VIM; j++) {
        FVP[i][j] = evpl_glob[0]->F_vp_glob[ielem][ip][i][j];
      }
    }

    for (j = 0; j < 3; j++) {
      local_post[EVP_DEF_GRAD_TENSOR + j] = FVP[j][j];
      local_lumped[EVP_DEF_GRAD_TENSOR + j] = 1.;
    }
    local_post[EVP_DEF_GRAD_TENSOR + 3] = FVP[0][1];
    local_post[EVP_DEF_GRAD_TENSOR + 4] = FVP[1][0];
    local_lumped[EVP_DEF_GRAD_TENSOR + 3] = 1.;
    local_lumped[EVP_DEF_GRAD_TENSOR + 4] = 1.;
    if (dim > 2) {
      local_post[EVP_DEF_GRAD_TENSOR + 5] = FVP[0][2];
      local_post[EVP_DEF_GRAD_TENSOR + 6] = FVP[2][0];
      local_post[EVP_DEF_GRAD_TENSOR + 7] = FVP[1][2];
      local_post[EVP_DEF_GRAD_TENSOR + 8] = FVP[2][1];
      local_lumped[EVP_DEF_GRAD_TENSOR + 5] = 1.;
      local_lumped[EVP_DEF_GRAD_TENSOR + 6] = 1.;
      local_lumped[EVP_DEF_GRAD_TENSOR + 7] = 1.;
      local_lumped[EVP_DEF_GRAD_TENSOR + 8] = 1.;
    }
    /* To restart the EVP calculation, you also need the stress tensor.
       We dump here if mesh-stress is requested !!*/
    if (STRESS_TENSOR != -1 && pd->e[pg->imtrx][R_MESH1]) {
      for (i = 0; i < DIM; i++) {
        for (j = 0; j < DIM; j++) {
          TVP[i][j] = evpl_glob[0]->TT_glob[ielem][ip][i][j];
        }
      }
      if (dim > 2) {
        w1 = 9;
      } else {
        w1 = 5;
      }
      for (j = 0; j < 3; j++) {
        local_post[EVP_DEF_GRAD_TENSOR + j + w1] = TVP[j][j];
        local_lumped[EVP_DEF_GRAD_TENSOR + j + w1] = 1.;
      }
      local_post[EVP_DEF_GRAD_TENSOR + 3 + w1] = TVP[0][1];
      local_post[EVP_DEF_GRAD_TENSOR + 4 + w1] = TVP[1][0];
      local_lumped[EVP_DEF_GRAD_TENSOR + 3 + w1] = 1.;
      local_lumped[EVP_DEF_GRAD_TENSOR + 4 + w1] = 1.;
      if (dim > 2) {
        local_post[EVP_DEF_GRAD_TENSOR + 5 + w1] = TVP[0][2];
        local_post[EVP_DEF_GRAD_TENSOR + 6 + w1] = TVP[2][0];
        local_post[EVP_DEF_GRAD_TENSOR + 7 + w1] = TVP[1][2];
        local_post[EVP_DEF_GRAD_TENSOR + 8 + w1] = TVP[2][1];
        local_lumped[EVP_DEF_GRAD_TENSOR + 5 + w1] = 1.;
        local_lumped[EVP_DEF_GRAD_TENSOR + 6 + w1] = 1.;
        local_lumped[EVP_DEF_GRAD_TENSOR + 7 + w1] = 1.;
        local_lumped[EVP_DEF_GRAD_TENSOR + 8 + w1] = 1.;
      }
    }
  }
  /* end of EVP_DEF_GRAD_TENSOR */

  if (LAGRANGE_CONVECTION != -1 &&
      (pd->MeshMotion == LAGRANGIAN || cr->MeshMotion == DYNAMIC_LAGRANGIAN)) {
    /* get the convection velocity (it's different for arbitrary and
       lagrangian meshes) */
    err = get_convection_velocity(vconv, vconv_old, d_vconv, delta_t, theta);
    GOMA_EH(err, "Error in calculating effective convection velocity");

    for (j = 0; j < dim; j++) {
      local_post[LAGRANGE_CONVECTION + j] = vconv[j];
      local_lumped[LAGRANGE_CONVECTION + j] = 1.;
    }
  }

  /*
   * Porous Media post-processing
   */
  checkPorous = 0;
  if (mp->PorousMediaType == POROUS_UNSATURATED || mp->PorousMediaType == POROUS_TWO_PHASE ||
      mp->PorousMediaType == POROUS_SHELL_UNSATURATED || mp->PorousMediaType == POROUS_SATURATED) {
    checkPorous = 1;
  }

  if (POROUS_SATURATION != -1 && checkPorous) {
    local_post[POROUS_SATURATION] = mp->saturation;
    local_lumped[POROUS_SATURATION] = 1.;
  } /* end of POROUS_SATURATION */

  if (POROUS_RHO_TOTAL_SOLVENTS != -1 && checkPorous) {
    w = 0;
    if (pd->v[pg->imtrx][POR_LIQ_PRES]) {
      local_post[POROUS_RHO_TOTAL_SOLVENTS] = pmv->bulk_density[i_pl];
      local_lumped[POROUS_RHO_TOTAL_SOLVENTS] = 1.;
    }
    w++;
    if (Num_Var_In_Type[pg->imtrx][R_POR_GAS_PRES]) {
      if (pd->v[pg->imtrx][POR_GAS_PRES]) {
        local_post[POROUS_RHO_TOTAL_SOLVENTS + w] = pmv->bulk_density[i_pg];
        local_lumped[POROUS_RHO_TOTAL_SOLVENTS + w] = 1.;
      }
      w++;
    }
    if (pd->v[pg->imtrx][POR_POROSITY]) {
      local_post[POROUS_RHO_TOTAL_SOLVENTS + w] = pmv->bulk_density[i_pore];
      local_lumped[POROUS_RHO_TOTAL_SOLVENTS + w] = 1.;
    }
  } /* end of POROUS_RHO_TOTAL_SOLVENTS */

  if (POROUS_RHO_GAS_SOLVENTS != -1 && checkPorous) {
    if (mp->PorousMediaType == POROUS_UNSATURATED || mp->PorousMediaType == POROUS_TWO_PHASE ||
        mp->PorousMediaType == POROUS_SHELL_UNSATURATED) {
      w = 0;
      if (Num_Var_In_Type[pg->imtrx][R_POR_LIQ_PRES]) {
        if (pd->v[pg->imtrx][POR_LIQ_PRES]) {
          local_post[POROUS_RHO_GAS_SOLVENTS] = pmv->gas_density_solvents[i_pl];
          local_lumped[POROUS_RHO_GAS_SOLVENTS] = 1.;
        }
        w++;
      }
      if (Num_Var_In_Type[pg->imtrx][R_POR_GAS_PRES]) {
        if (pd->v[pg->imtrx][POR_GAS_PRES]) {
          local_post[POROUS_RHO_GAS_SOLVENTS + w] = pmv->gas_density_solvents[i_pg];
          local_lumped[POROUS_RHO_GAS_SOLVENTS + w] = 1.;
        }
        w++;
      }
      if (Num_Var_In_Type[pg->imtrx][R_POR_POROSITY]) {
        if (pd->v[pg->imtrx][POR_POROSITY]) {
          local_post[POROUS_RHO_GAS_SOLVENTS + w] = pmv->gas_density_solvents[i_pore];
          local_lumped[POROUS_RHO_GAS_SOLVENTS + w] = 1.;
        }
        w++;
      }
    }
  } /* end of POROUS_RHO_GAS_SOLVENTS */

  if (POROUS_RHO_LPHASE != -1 && checkPorous) {
    local_post[POROUS_RHO_LPHASE] = mp->density * mp->porosity * mp->saturation;
    local_lumped[POROUS_RHO_LPHASE] = 1.;
  } /* end of POROUS_RHO_LPHASE */

  if (DARCY_VELOCITY_GAS != -1 && mp->PorousMediaType == POROUS_TWO_PHASE) {
    for (a = 0; a < dim; a++) {
      local_post[DARCY_VELOCITY_GAS + a] = pmv->gas_darcy_velocity[a];
      local_lumped[DARCY_VELOCITY_GAS + a] = 1.;
    }
  } /* end of DARCY_VELOCITY_GAS */

  if (DARCY_VELOCITY_LIQ != -1 &&
      (mp->PorousMediaType == POROUS_UNSATURATED ||
       mp->PorousMediaType == POROUS_SHELL_UNSATURATED || mp->PorousMediaType == POROUS_SATURATED ||
       mp->PorousMediaType == POROUS_TWO_PHASE)) {
    for (a = 0; a < dim; a++) {
      local_post[DARCY_VELOCITY_LIQ + a] = pmv->liq_darcy_velocity[a];
      local_lumped[DARCY_VELOCITY_LIQ + a] = 1.;
    }
  } /* end of DARCY_VELOCITY_LIQ */

  if (POROUS_LIQUID_ACCUM_RATE != -1 && (mp->PorousMediaType == POROUS_UNSATURATED ||
                                         mp->PorousMediaType == POROUS_SHELL_UNSATURATED ||
                                         mp->PorousMediaType == POROUS_TWO_PHASE)) {
    local_post[POROUS_LIQUID_ACCUM_RATE] = pmt->Inventory_solvent_dot[0];
    for (a = 0; a < dim; a++) {
      local_post[POROUS_LIQUID_ACCUM_RATE] += pmt->conv_flux[0][a];
    }
    local_lumped[POROUS_LIQUID_ACCUM_RATE] = 1.;
  } /* end of POROUS_LIQUID_ACCUM_RATE */

  if (ELECTRIC_FIELD != -1 && pd->e[pg->imtrx][R_POTENTIAL]) {
    for (a = 0; a < dim; a++) {
      local_post[ELECTRIC_FIELD + a] = -fv->grad_V[a];
      local_lumped[ELECTRIC_FIELD + a] = 1.0;
    }
  } /* end of ELECTRIC_FIELD */

  if (ACOUSTIC_PRESSURE != -1 &&
      (pd->e[pg->imtrx][R_ACOUS_PREAL] || pd->e[pg->imtrx][R_ACOUS_PIMAG])) {
    local_post[ACOUSTIC_PRESSURE] = sqrt(fv->apr * fv->apr + fv->api * fv->api);
    local_lumped[ACOUSTIC_PRESSURE] = 1.0;
  } /* end of ACOUSTIC_PRESSURE */

  if (ACOUSTIC_PHASE_ANGLE != -1 &&
      (pd->e[pg->imtrx][R_ACOUS_PREAL] || pd->e[pg->imtrx][R_ACOUS_PIMAG])) {
    local_post[ACOUSTIC_PHASE_ANGLE] = atan2(fv->api, fv->apr) * 180 / M_PIE;
    local_lumped[ACOUSTIC_PHASE_ANGLE] = 1.0;
  } /* end of ACOUSTIC_PHASE_ANGLE */

  if (ACOUSTIC_ENERGY_DENSITY != -1 &&
      (pd->e[pg->imtrx][R_ACOUS_PREAL] || pd->e[pg->imtrx][R_ACOUS_PIMAG])) {
    double acous_pgrad = 0;
    double k, R, omega;
    k = wave_number(NULL, time);
    R = acoustic_impedance(NULL, time);
    omega = upd->Acoustic_Frequency;
    for (a = 0; a < dim; a++) {
      acous_pgrad += fv->grad_api[a] * fv->grad_api[a];
      acous_pgrad += fv->grad_apr[a] * fv->grad_apr[a];
    }
    local_post[ACOUSTIC_ENERGY_DENSITY] =
        0.25 * k / (R * omega) * (fv->apr * fv->apr + fv->api * fv->api + acous_pgrad / (k * k));
    local_lumped[ACOUSTIC_ENERGY_DENSITY] = 1.0;
  } /* end of ACOUSTIC_ENERGY_DENSITY */

  if (LIGHT_INTENSITY != -1 && (pd->e[pg->imtrx][R_LIGHT_INTP] || pd->e[pg->imtrx][R_LIGHT_INTM])) {
    local_post[LIGHT_INTENSITY] = fv->poynt[0] + fv->poynt[1];
    local_lumped[LIGHT_INTENSITY] = 1.0;
  } /* end of LIGHT_INTENSITY */

  if (UNTRACKED_SPEC != -1 && pd->e[pg->imtrx][R_MASS]) {
    double density_tot = 0.;
    switch (mp->Species_Var_Type) {
    case SPECIES_CONCENTRATION:
      density_tot = calc_density(mp, FALSE, NULL, 0.0);
      local_post[UNTRACKED_SPEC] = density_tot;
      for (j = 0; j < pd->Num_Species_Eqn; j++) {
        local_post[UNTRACKED_SPEC] -= fv->c[j] * mp->molecular_weight[j];
      }
      local_post[UNTRACKED_SPEC] /= mp->molecular_weight[pd->Num_Species_Eqn];
      break;
    case SPECIES_DENSITY:
      density_tot = calc_density(mp, FALSE, NULL, 0.0);
      local_post[UNTRACKED_SPEC] = density_tot;
      for (j = 0; j < pd->Num_Species_Eqn; j++) {
        local_post[UNTRACKED_SPEC] -= fv->c[j];
      }
      break;
    case SPECIES_MASS_FRACTION:
    case SPECIES_UNDEFINED_FORM:
      local_post[UNTRACKED_SPEC] = 1.0;
      for (j = 0; j < pd->Num_Species_Eqn; j++) {
        local_post[UNTRACKED_SPEC] -= fv->c[j];
      }
      break;
    default:
      GOMA_WH(-1, "Undefined Species Type in UNTRACKED_SPEC\n");
    }
    local_lumped[UNTRACKED_SPEC] = 1.0;
  } /* end of UNTRACKED_SPEC*/

  if ((TFMP_GAS_VELO != -1 || TFMP_LIQ_VELO != -1 || TFMP_KRG != -1) &&
      pd->e[pg->imtrx][R_TFMP_MASS] && pd->e[pg->imtrx][R_TFMP_BOUND]) {

    int k;
    int *n_dof = NULL;
    int dof_map[MDE];
    double v_l[DIM], v_g[DIM];
    /*
     * Prepare geometry
     */
    n_dof = (int *)array_alloc(1, MAX_VARIABLE_TYPES, sizeof(int));
    lubrication_shell_initialize(n_dof, dof_map, -1, xi, exo, 0);

    /* Gather necessary values (S, h, Krg, gradII_P)*/

    // need pure phase viscosities
    double mu_l, mu_g;

    load_tfmp_viscosity_model(&mu_l, &mu_g);

    double S = fv->tfmp_sat;
    /* Use the height_function_model */
    double H_U, dH_U_dtime, H_L, dH_L_dtime;
    double dH_U_dX[DIM], dH_L_dX[DIM], dH_U_dp, dH_U_ddh;
    double h = height_function_model(&H_U, &dH_U_dtime, &H_L, &dH_L_dtime, dH_U_dX, dH_L_dX,
                                     &dH_U_dp, &dH_U_ddh, time, delta_t);

    double dh_dtime = dH_U_dtime - dH_L_dtime;
    double gradII_h[DIM];
    for (k = 0; k < DIM; k++) {
      gradII_h[k] = dH_U_dX[k] - dH_L_dX[k];
    }
    double dh_dmesh[DIM][MDE];
    double dh_dnormal[DIM][MDE];
    double d2h_dtime_dmesh[DIM][MDE];
    double d2h_dtime_dnormal[DIM][MDE];
    double d_gradIIh_dmesh[DIM][DIM][MDE];
    double d_gradIIh_dnormal[DIM][DIM][MDE];

    // dh_dtime is not used here
    double tt = 1.0;

    load_displacement_coupling_model(tt, delta_t, &h, &dh_dtime, gradII_h, dh_dmesh, dh_dnormal,
                                     d2h_dtime_dmesh, d2h_dtime_dnormal, d_gradIIh_dmesh,
                                     d_gradIIh_dnormal, n_dof, dof_map);
    //  rel perms
    double Krl, dKrl_dS, Krg, dKrg_dS;
    load_relative_permeability_model(S, &Krl, &dKrl_dS, &Krg, &dKrg_dS);

    /* Use the velocity function model */
    double veloU[DIM], veloL[DIM], veloAVG[DIM];

    velocity_function_model(veloU, veloL, time, delta_t);

    for (k = 0; k < DIM; k++) {
      veloAVG[k] = (veloU[k] + veloL[k]) / 2.;
    }
    veloAVG[2] = 0.0;

    double gradII_P[DIM];
    Inn(fv->grad_tfmp_pres, gradII_P);

    /* Calculate Velocity */
    for (k = 0; k < DIM; k++) {
      v_l[k] = -h * h / 12.0 / mu_l * Krl * gradII_P[k];
      v_g[k] = -h * h / 12.0 / mu_g * Krg * gradII_P[k];
    }
    if (TFMP_GAS_VELO != -1) {
      local_post[TFMP_GAS_VELO] = v_g[0] + veloAVG[0];
      local_lumped[TFMP_GAS_VELO] += 1;
      local_post[TFMP_GAS_VELO + 1] = v_g[1] + veloAVG[1];
      local_lumped[TFMP_GAS_VELO + 1] += 1;
      local_post[TFMP_GAS_VELO + 2] = v_g[2] + veloAVG[2];
      local_lumped[TFMP_GAS_VELO + 2] += 1;
    }
    if (TFMP_LIQ_VELO != -1) {
      local_post[TFMP_LIQ_VELO] = v_l[0] + veloAVG[0];
      local_lumped[TFMP_LIQ_VELO] += 1;
      local_post[TFMP_LIQ_VELO + 1] = v_l[1] + veloAVG[1];
      local_lumped[TFMP_LIQ_VELO + 1] += 1;
      local_post[TFMP_LIQ_VELO + 2] = v_l[2] + veloAVG[2];
      local_lumped[TFMP_LIQ_VELO + 2] += 1;
    }
    if (TFMP_KRG != -1) {
      local_post[TFMP_KRG] = Krg;
      local_lumped[TFMP_KRG] += 1;
      //    } else {
      //      local_post[TFMP_KRG] = 0.0;
    }

    /* Cleanup */
    safe_free((void *)n_dof);
  }
  if ((TFMP_INV_PECLET != -1) && pd->e[pg->imtrx][R_TFMP_MASS] && pd->e[pg->imtrx][R_TFMP_BOUND]) {

    int k;
    int *n_dof = NULL;
    int dof_map[MDE];

    /*
     * Prepare geometry
     */
    n_dof = (int *)array_alloc(1, MAX_VARIABLE_TYPES, sizeof(int));
    lubrication_shell_initialize(n_dof, dof_map, -1, xi, exo, 0);

    /* Gather necessary values (S, h, Krg, gradII_P)*/

    // need pure phase viscosities
    double mu_l, mu_g;

    load_tfmp_viscosity_model(&mu_l, &mu_g);

    double S = fv->tfmp_sat;
    /* Use the height_function_model */
    double H_U, dH_U_dtime, H_L, dH_L_dtime;
    double dH_U_dX[DIM], dH_L_dX[DIM], dH_U_dp, dH_U_ddh;
    double h = height_function_model(&H_U, &dH_U_dtime, &H_L, &dH_L_dtime, dH_U_dX, dH_L_dX,
                                     &dH_U_dp, &dH_U_ddh, time, delta_t);

    double dh_dtime = dH_U_dtime - dH_L_dtime;
    double gradII_h[DIM];
    for (k = 0; k < DIM; k++) {
      gradII_h[k] = dH_U_dX[k] - dH_L_dX[k];
    }

    double dh_dmesh[DIM][MDE];
    double dh_dnormal[DIM][MDE];
    double d2h_dtime_dmesh[DIM][MDE];
    double d2h_dtime_dnormal[DIM][MDE];
    double d_gradIIh_dmesh[DIM][DIM][MDE];
    double d_gradIIh_dnormal[DIM][DIM][MDE];

    // dh_dtime is not used here
    double tt = 1.0;

    load_displacement_coupling_model(tt, delta_t, &h, &dh_dtime, gradII_h, dh_dmesh, dh_dnormal,
                                     d2h_dtime_dmesh, d2h_dtime_dnormal, d_gradIIh_dmesh,
                                     d_gradIIh_dnormal, n_dof, dof_map);

    // Artificial diffusion constant
    double D, Krd, dKrd_dS;
    load_molecular_diffusion_model(S, &D, &Krd, &dKrd_dS);

    //  rel perms
    double Krl, dKrl_dS, Krg, dKrg_dS;
    load_relative_permeability_model(S, &Krl, &dKrl_dS, &Krg, &dKrg_dS);

    double gradII_P[DIM], gradII_S[DIM];

    Inn(fv->grad_tfmp_pres, gradII_P);
    Inn(fv->grad_tfmp_sat, gradII_S);

    dbl mag_gradII_P = 0.0;
    dbl mag_gradII_S = 0.0;

    /* Calculate gradient magnitudes */
    for (k = 0; k < DIM; k++) {
      mag_gradII_P += gradII_P[k] * gradII_P[k];
      mag_gradII_S += gradII_S[k] * gradII_S[k];
    }
    mag_gradII_P = sqrt(mag_gradII_P);
    mag_gradII_S = sqrt(mag_gradII_S);

    // if (TFMP_GAS_VELO != -1) {
    // if ( Krd != 0.0 && mag_gradII_S >= 1.e-10 ) {
    if (Krd != 0.0 && mag_gradII_S != 0.0) {
      local_post[TFMP_INV_PECLET] =
          (Krl * h * h * h / 12.0 / mu_l * mag_gradII_P) / (D * Krd * mag_gradII_S);
    } else {
      local_post[TFMP_INV_PECLET] = 0.0;
    }
    local_lumped[TFMP_INV_PECLET] += 1;

    //}
    /* Cleanup */
    safe_free((void *)n_dof);
  }

  /*  EXTERNAL tables	*/
  if (efv->ev) {
    for (j = 0; j < efv->Num_external_field; j++) {
      if (efv->i[j] == I_TABLE) {
        local_post[EXTERNAL_POST + j] = fv->external_field[j];
        local_lumped[EXTERNAL_POST + j] = 1.0;
      }
    }
  }

  if (ELECTRIC_FIELD_MAG != -1 && pd->e[pg->imtrx][R_POTENTIAL]) {
    for (a = 0; a < dim; a++) {
      local_post[ELECTRIC_FIELD_MAG] += fv->grad_V[a] * fv->grad_V[a];
    }
    local_post[ELECTRIC_FIELD_MAG] = sqrt(local_post[ELECTRIC_FIELD_MAG]);
    local_lumped[ELECTRIC_FIELD_MAG] = 1.0;
  } /* end of ELECTRIC_FIELD_MAG */

  if (CAPILLARY_PRESSURE != -1 && (mp->PorousMediaType == POROUS_UNSATURATED ||
                                   mp->PorousMediaType == POROUS_SHELL_UNSATURATED ||
                                   mp->PorousMediaType == POROUS_TWO_PHASE)) {
    local_post[CAPILLARY_PRESSURE] = pmv->cap_pres;
    local_lumped[CAPILLARY_PRESSURE] = 1.;
  } /* end of CAPILLARY_PRESSURE */

  if (POROUS_GRIDPECLET != -1 && checkPorous) {
    local_post[POROUS_GRIDPECLET] = Stab->Grid_Peclet_Number[POR_LIQ_PRES];
    local_lumped[POROUS_GRIDPECLET] = 1.0;
  } /* end of POROUS_GRIDPECLET */

  if (POROUS_SUPGVELOCITY != -1 && checkPorous) {
    for (a = 0; a < dim; a++) {
      local_post[POROUS_SUPGVELOCITY + a] = pmv->U_supg[a];
      local_lumped[POROUS_SUPGVELOCITY + a] = 1.0;
    }
  } /* end of POROUS_SUPGVELOCITY */

  if (CURL_V != -1 && pd->e[pg->imtrx][R_MOMENTUM1]) {
    /* MMH: Note that in the SWIRLING coordinate system, we really
     * do have a 3-vector and not just a scalar.
     *
     * Same for PROJECTED_CARTESIAN, if VELOCITY3 can become
     * non-zero. -MMH
     */
    if (pd->CoordinateSystem == SWIRLING || pd->CoordinateSystem == PROJECTED_CARTESIAN ||
        pd->CoordinateSystem == CARTESIAN_2pt5D || dim == 3) {
      for (j = 0; j < VIM; j++) {
        local_post[CURL_V + j] = fv->curl_v[j];
        local_lumped[CURL_V + j] = 1.0;
      }
    } else {
      /* We are in 2D or axisymmetric.  In either case, it is the
       * last component of the vorticity vector that contains the
       * magnitude.
       */
      local_post[CURL_V] = fv->curl_v[2];
      local_lumped[CURL_V] = 1.0;
    }
  }

  if (HELICITY != -1 && pd->e[pg->imtrx][R_MOMENTUM1]) {
    /* MMH: Note that in the SWIRLING coordinate system, we really
     * do have a 3-vector and not just a scalar.
     *
     * Same for PROJECTED_CARTESIAN, if VELOCITY3 can become
     * non-zero. -MMH
     */
    if (pd->CoordinateSystem == SWIRLING || pd->CoordinateSystem == PROJECTED_CARTESIAN ||
        pd->CoordinateSystem == CARTESIAN_2pt5D || dim == 3) {
      local_post[HELICITY] = 0.0;
      for (j = 0; j < VIM; j++) {
        local_post[HELICITY] += fv->v[j] * fv->curl_v[j];
      }
      local_lumped[HELICITY] = 1.0;
    }
  }

  if (LAMB_VECTOR != -1 && pd->e[pg->imtrx][R_MOMENTUM1]) {
    for (j = 0; j < VIM; j++) {
      local_post[LAMB_VECTOR + j] = 0;
      for (a = 0; a < VIM; a++) {
        for (b = 0; b < VIM; b++) {
          local_post[LAMB_VECTOR + j] += permute(a, b, j) * fv->curl_v[a] * fv->v[b];
        }
      }
      local_lumped[LAMB_VECTOR + j] = 1.0;
    }
  }

  /* calculate poynting vectors for EM calculations here !!  */
  if (POYNTING_VECTORS != -1 &&
      ((Num_Var_In_Type[pg->imtrx][R_ACOUS_PREAL] || Num_Var_In_Type[pg->imtrx][R_ACOUS_PIMAG]))) {
    double poynt[DIM];
    memset(poynt, 0, sizeof(double) * DIM);
    /*  Acoustic analogy -- scalar version  */
    if (pd->e[pg->imtrx][R_ACOUS_PREAL] || pd->e[pg->imtrx][R_ACOUS_PIMAG]) {
      double k, R, prefactor;
      k = wave_number(NULL, time);
      R = acoustic_impedance(NULL, time);
      prefactor = 0.5 / (k * R);

      for (a = 0; a < DIM; a++) {
        poynt[a] += prefactor * (fv->api * fv->grad_apr[a] - fv->apr * fv->grad_api[a]);
      }
    }
    for (a = 0; a < dim; a++) {
      local_post[POYNTING_VECTORS + a] = poynt[a];
      local_lumped[POYNTING_VECTORS + a] = 1.;
    }
  }

  /* calculate species sources  */
  if (SPECIES_SOURCES != -1 && pd->e[pg->imtrx][R_MASS]) {
    err = get_continuous_species_terms(&s_terms, time, theta, delta_t, hs);
    for (w = 0; w < pd->Num_Species_Eqn; w++) {
      local_post[SPECIES_SOURCES + w] = s_terms.MassSource[w];
      local_lumped[SPECIES_SOURCES + w] = 1.;
    }
  }

  if (STRESS_NORM != -1 && pd->e[pg->imtrx][POLYMER_STRESS11]) {
    for (int mode = 0; mode < vn->modes; mode++) {

      dbl traceOverVim = 0;
      for (int i = 0; i < VIM; i++) {
        traceOverVim += fv->S[mode][i][i];
      }

      traceOverVim /= VIM;

      // square of the deviatoric sress norm
      dbl normOfStressDSqr = 0;
      for (int i = 0; i < VIM; i++) {
        normOfStressDSqr += pow(fv->S[mode][i][i] - traceOverVim, 2) / 2.;

        for (int j = i + 1; j < VIM; j++) {
          normOfStressDSqr += pow(fv->S[mode][i][j], 2);
        }
      }

      dbl normOfStressD = sqrt(normOfStressDSqr);
      local_post[STRESS_NORM + mode] = normOfStressD;
      local_lumped[STRESS_NORM + mode] = 1.;
    }
  }

  if (SARAMITO_YIELD != -1 && pd->e[pg->imtrx][POLYMER_STRESS11]) {
    for (int mode = 0; mode < vn->modes; mode++) {
      dbl coeff;
      compute_saramito_model_terms(&coeff, NULL, fv->S[mode], ve[mode]->gn, TRUE);
      local_post[SARAMITO_YIELD + mode] = coeff;
      local_lumped[SARAMITO_YIELD + mode] = 1.;
    }
  }

  if (VISCOUS_STRESS != -1 && pd->e[pg->imtrx][R_MOMENTUM1]) {
    for (a = 0; a < VIM; a++) {
      for (b = 0; b < VIM; b++) {
        gamma[a][b] = fv->grad_v[a][b] + fv->grad_v[b][a];
      }
    }

    mu = viscosity(gn, gamma, NULL);
    // printf("%lf", mu);
    for (a = 0; a < VIM; a++) {
      for (b = 0; b < VIM; b++) {
        local_post[VISCOUS_STRESS + a * VIM + b] = mu * gamma[a][b];
        local_lumped[VISCOUS_STRESS + a * VIM + b] = 1.;
      }
    }
  }

  if (VISCOUS_STRESS_NORM != -1 && pd->e[pg->imtrx][R_MOMENTUM1]) {
    for (a = 0; a < VIM; a++) {
      for (b = 0; b < VIM; b++) {
        gamma[a][b] = fv->grad_v[a][b] + fv->grad_v[b][a];
      }
    }

    mu = viscosity(gn, gamma, NULL);

    dbl traceOverVim = 0;
    for (int i = 0; i < VIM; i++) {
      traceOverVim += mu * (fv->grad_v[i][i] + fv->grad_v[i][i]);
    }

    traceOverVim /= VIM;

    // square of the deviatoric sress norm
    dbl normOfStressDSqr = 0;
    for (int i = 0; i < VIM; i++) {
      normOfStressDSqr += pow(mu * (fv->grad_v[i][i] + fv->grad_v[i][i]) - traceOverVim, 2) / 2.;

      for (int j = i + 1; j < VIM; j++) {
        normOfStressDSqr += pow(mu * (fv->grad_v[i][j] + fv->grad_v[j][i]), 2);
      }
    }

    dbl normOfStressD = sqrt(normOfStressDSqr);
    local_post[VISCOUS_STRESS_NORM] = normOfStressD;
    local_lumped[VISCOUS_STRESS_NORM] = 1.;
  }

  if (VISCOUS_VON_MISES_STRESS != -1 && pd->e[pg->imtrx][R_MOMENTUM1]) {
    for (a = 0; a < VIM; a++) {
      for (b = 0; b < VIM; b++) {
        gamma[a][b] = fv->grad_v[a][b] + fv->grad_v[b][a];
      }
    }

    mu = viscosity(gn, gamma, NULL);

    dbl viscousVonMisesStress = 0.0;

    if (Num_Dim < 3) {
      dbl stress_xx = mu * (fv->grad_v[0][0] + fv->grad_v[0][0]);
      dbl stress_xy = mu * (fv->grad_v[0][1] + fv->grad_v[1][0]);
      dbl stress_yy = mu * (fv->grad_v[1][1] + fv->grad_v[1][1]);

      viscousVonMisesStress = sqrt(pow(stress_xx, 2) + pow(stress_yy, 2) + 3 * pow(stress_xy, 2) -
                                   stress_xx * stress_yy);

    } else if (Num_Dim > 2) {
      dbl stress_xx = mu * (fv->grad_v[0][0] + fv->grad_v[0][0]);
      dbl stress_yy = mu * (fv->grad_v[1][1] + fv->grad_v[1][1]);
      dbl stress_zz = mu * (fv->grad_v[2][2] + fv->grad_v[2][2]);

      dbl stress_xy = mu * (fv->grad_v[0][1] + fv->grad_v[1][0]);
      dbl stress_xz = mu * (fv->grad_v[0][2] + fv->grad_v[2][0]);
      dbl stress_yz = mu * (fv->grad_v[1][2] + fv->grad_v[2][1]);

      viscousVonMisesStress =
          sqrt(0.5 * (pow((stress_xx - stress_yy), 2) + pow((stress_yy - stress_zz), 2) +
                      pow((stress_zz - stress_xx), 2) +
                      6 * (pow(stress_xy, 2) + pow(stress_yz, 2) + pow(stress_xz, 2))));
    }

    // dbl viscousVonMisesStress1 = viscousVonMisesStress;
    local_post[VISCOUS_VON_MISES_STRESS] = viscousVonMisesStress;
    local_lumped[VISCOUS_VON_MISES_STRESS] = 1.;
  }

  /* calculate real-solid stress here !!  */
  if (REAL_STRESS_TENSOR != -1 && pd->e[pg->imtrx][R_SOLID1]) {
    mu = elc_rs->lame_mu;
    err = belly_flop_rs(mu);
    GOMA_EH(err, "error in belly flop");
    if (err == 2)
      return (err);
    /*
     * Total mesh stress tensor...
     */
    err = solid_stress_tensor(TT, dTT_dx, dTT_drs, dTT_dp, dTT_dc, dTT_dp_liq, dTT_dp_gas,
                              dTT_dporosity, dTT_dT, dTT_dmax_strain, elc_rs->lame_mu,
                              elc_rs->lame_lambda);

    if (dim == 2) {
      if (cr->RealSolidFluxModel == NONLINEAR || cr->RealSolidFluxModel == HOOKEAN_PSTRAIN ||
          cr->RealSolidFluxModel == INCOMP_PSTRAIN)
        TT[2][2] = (1. - pow(fv->volume_change, 2. / 3.)) * elc_rs->lame_mu - fv->P;
      /*              if (cr->MeshFluxModel == INCOMP_PSTRESS) */
      else
        TT[2][2] = 0.;

      TT[1][2] = 0.;
      TT[0][2] = 0.;
    }

    for (j = 0; j < DIM; j++) {
      local_post[REAL_STRESS_TENSOR + j] = TT[j][j];
      local_lumped[REAL_STRESS_TENSOR + j] = 1.;
    }
    local_post[REAL_STRESS_TENSOR + 3] = TT[0][1];
    local_lumped[REAL_STRESS_TENSOR + 3] = 1.;
    if (dim == 3) {
      local_post[REAL_STRESS_TENSOR + 4] = TT[0][2];
      local_post[REAL_STRESS_TENSOR + 5] = TT[1][2];
      local_lumped[REAL_STRESS_TENSOR + 4] = 1.;
      local_lumped[REAL_STRESS_TENSOR + 5] = 1.;
    }

  } /* end of REAL_STRESS_TENSOR */

  /* calculate principal stress differences*/
  if (PRINCIPAL_STRESS != -1 && pd->e[pg->imtrx][R_MESH1]) {
    double I_T, II_T, III_T, coeff_a, coeff_b;
    double m_par = 0, theta1, evalue1, evalue2, evalue3;
    if (cr->MeshMotion != ARBITRARY) {
      /*
       * Total mesh stress tensor...
       */
      err = mesh_stress_tensor(TT, dTT_dx, dTT_dp, dTT_dc, dTT_dp_liq, dTT_dp_gas, dTT_dporosity,
                               dTT_dsink_mass, dTT_dT, dTT_dmax_strain, dTT_dcur_strain,
                               elc->lame_mu, elc->lame_lambda, delta_t, ielem, ip, ip_total);

      /* For LINEAR ELASTICITY */
      if (cr->MeshFluxModel == LINEAR) {
        if (dim == 2) {
          TT[2][2] = 1.;
          TT[1][2] = 0.;
          TT[0][2] = 0.;
        }
      }
      /*  For Hookian Elasticity and shrinkage */
      else {
        if (dim == 2) {
          if (cr->MeshMotion == ARBITRARY) {
            TT[2][2] = (1. - fv->volume_change) * elc->lame_mu;
          } else {
            if (cr->MeshFluxModel == NONLINEAR || cr->MeshFluxModel == HOOKEAN_PSTRAIN ||
                cr->MeshFluxModel == INCOMP_PSTRAIN)
              TT[2][2] = (1. - pow(fv->volume_change, 2. / 3.)) * elc->lame_mu - fv->P;
            /*              if (cr->MeshFluxModel == INCOMP_PSTRESS) */
            else
              TT[2][2] = 0.;
          }
          TT[1][2] = 0.;
          TT[0][2] = 0.;
        }
      }
    } else {
      memset(TT, 0, sizeof(dbl) * DIM * DIM);
      memset(gamma, 0, sizeof(dbl) * DIM * DIM);
      for (a = 0; a < VIM; a++) {
        for (b = 0; b < VIM; b++) {
          gamma[a][b] = fv->grad_v[a][b] + fv->grad_v[b][a];
        }
      }
      mu = viscosity(gn, gamma, NULL);
      for (a = 0; a < VIM; a++) {
        for (b = 0; b < VIM; b++) {
          TT[a][b] = mu * gamma[a][b] - fv->P * delta(a, b);
        }
      }
      if (pd->v[pg->imtrx][POLYMER_STRESS11]) {
        for (a = 0; a < VIM; a++) {
          for (b = 0; b < VIM; b++) {
            for (mode = 0; mode < vn->modes; mode++) {
              TT[a][b] += fv->S[mode][a][b];
            }
          }
        }
      }
    }
    /*  try for Trig solution of cubic equation     */
    I_T = TT[0][0] + TT[1][1] + TT[2][2];
    II_T = TT[0][0] * TT[1][1] + TT[0][0] * TT[2][2] + TT[1][1] * TT[2][2] -
           (SQUARE(TT[0][1]) + SQUARE(TT[0][2]) + SQUARE(TT[1][2]));
    III_T = TT[0][0] * TT[1][1] * TT[2][2] + 2. * (TT[0][1] * TT[1][2] * TT[0][2]) -
            TT[0][0] * SQUARE(TT[1][2]) - TT[1][1] * SQUARE(TT[0][2]) - TT[2][2] * SQUARE(TT[0][1]);
    coeff_a = (3. * II_T - SQUARE(I_T)) / 3.;
    coeff_b = (2. * (-I_T) * SQUARE(I_T) - 9. * (-I_T) * II_T + 27. * (-III_T)) / 27.;
    if (coeff_a > 0) {
      fprintf(stderr, "trouble - imaginary roots %g %g\n", coeff_a, coeff_b);
    } else {
      m_par = 2. * sqrt(-coeff_a / 3.);
    }
    theta1 = acos(3. * coeff_b / (coeff_a * m_par)) / 3.;
    evalue1 = m_par * cos(theta1) + I_T / 3.;
    evalue2 = m_par * cos(theta1 + 2. * M_PIE / 3.) + I_T / 3.;
    evalue3 = m_par * cos(theta1 + 4. * M_PIE / 3.) + I_T / 3.;
    theta1 = evalue2;
    if (fabs(theta1) > fabs(evalue1)) {
      evalue2 = evalue1;
      evalue1 = theta1;
    }
    theta1 = evalue3;
    if (fabs(theta1) > fabs(evalue2)) {
      evalue3 = evalue2;
      evalue2 = theta1;
    }
    if (fabs(theta1) > fabs(evalue1)) {
      evalue2 = evalue1;
      evalue2 = theta1;
    }
    local_post[PRINCIPAL_STRESS] = evalue1;
    local_lumped[PRINCIPAL_STRESS] = 1.;
    local_post[PRINCIPAL_STRESS + 1] = evalue2;
    local_lumped[PRINCIPAL_STRESS + 1] = 1.;
    local_post[PRINCIPAL_STRESS + 2] = evalue3;
    local_lumped[PRINCIPAL_STRESS + 2] = 1.;
  } /* end of PRINCIPAL_STRESS */

  /* calculate principal real stress differences*/
  if (PRINCIPAL_REAL_STRESS != -1 && pd->e[pg->imtrx][R_SOLID1]) {
    double I_T, II_T, III_T, coeff_a, coeff_b;
    double m_par = 0, theta1, evalue1, evalue2, evalue3;
    /*
     * Total mesh stress tensor...
     */
    mu = elc_rs->lame_mu;
    err = belly_flop_rs(mu);
    GOMA_EH(err, "error in belly flop");
    if (err == 2)
      return (err);
    err = solid_stress_tensor(TT, dTT_dx, dTT_drs, dTT_dp, dTT_dc, dTT_dp_liq, dTT_dp_gas,
                              dTT_dporosity, dTT_dT, dTT_dmax_strain, elc_rs->lame_mu,
                              elc_rs->lame_lambda);

    if (dim == 2) {
      if (cr->RealSolidFluxModel == NONLINEAR || cr->RealSolidFluxModel == HOOKEAN_PSTRAIN ||
          cr->RealSolidFluxModel == INCOMP_PSTRAIN)
        TT[2][2] = (1. - pow(fv->volume_change, 2. / 3.)) * elc_rs->lame_mu - fv->P;
      /*              if (cr->MeshFluxModel == INCOMP_PSTRESS) */
      else
        TT[2][2] = 0.;

      TT[1][2] = 0.;
      TT[0][2] = 0.;
    }
    /*  try for Trig solution of cubic equation     */
    I_T = TT[0][0] + TT[1][1] + TT[2][2];
    II_T = TT[0][0] * TT[1][1] + TT[0][0] * TT[2][2] + TT[1][1] * TT[2][2] -
           (SQUARE(TT[0][1]) + SQUARE(TT[0][2]) + SQUARE(TT[1][2]));
    III_T = TT[0][0] * TT[1][1] * TT[2][2] + 2. * (TT[0][1] * TT[1][2] * TT[0][2]) -
            TT[0][0] * SQUARE(TT[1][2]) - TT[1][1] * SQUARE(TT[0][2]) - TT[2][2] * SQUARE(TT[0][1]);
    coeff_a = (3. * II_T - SQUARE(I_T)) / 3.;
    coeff_b = (2. * (-I_T) * SQUARE(I_T) - 9. * (-I_T) * II_T + 27. * (-III_T)) / 27.;
    if (coeff_a > 0) {
      fprintf(stderr, "trouble - imaginary roots %g %g\n", coeff_a, coeff_b);
    } else {
      m_par = 2. * sqrt(-coeff_a / 3.);
    }
    theta1 = acos(3. * coeff_b / (coeff_a * m_par)) / 3.;
    evalue1 = m_par * cos(theta1) + I_T / 3.;
    evalue2 = m_par * cos(theta1 + 2. * M_PIE / 3.) + I_T / 3.;
    evalue3 = m_par * cos(theta1 + 4. * M_PIE / 3.) + I_T / 3.;
    theta1 = evalue2;
    if (fabs(theta1) > fabs(evalue1)) {
      evalue2 = evalue1;
      evalue1 = theta1;
    }
    theta1 = evalue3;
    if (fabs(theta1) > fabs(evalue2)) {
      evalue3 = evalue2;
      evalue2 = theta1;
    }
    if (fabs(theta1) > fabs(evalue1)) {
      evalue2 = evalue1;
      evalue2 = theta1;
    }
    local_post[PRINCIPAL_REAL_STRESS] = evalue1;
    local_lumped[PRINCIPAL_REAL_STRESS] = 1.;
    local_post[PRINCIPAL_REAL_STRESS + 1] = evalue2;
    local_lumped[PRINCIPAL_REAL_STRESS + 1] = 1.;
    local_post[PRINCIPAL_REAL_STRESS + 2] = evalue3;
    local_lumped[PRINCIPAL_REAL_STRESS + 2] = 1.;
  } /* end of PRINCIPAL_REAL_STRESS */

  if (LUB_HEIGHT != -1 && (pd->e[pg->imtrx][R_LUBP] || pd->e[pg->imtrx][R_SHELL_FILMP] ||
                           pd->e[pg->imtrx][R_TFMP_MASS] || pd->e[pg->imtrx][R_TFMP_BOUND])) {
    double H_U, dH_U_dtime, H_L, dH_L_dtime;
    double dH_U_dX[DIM], dH_L_dX[DIM], dH_U_dp, dH_U_ddh;

    /* Setup lubrication */
    int *n_dof = NULL;
    int dof_map[MDE];
    dbl wt = fv->wt;
    n_dof = (int *)array_alloc(1, MAX_VARIABLE_TYPES, sizeof(int));
    lubrication_shell_initialize(n_dof, dof_map, -1, xi, exo, 0);

    if (pd->e[pg->imtrx][R_LUBP]) {
      local_post[LUB_HEIGHT] = height_function_model(&H_U, &dH_U_dtime, &H_L, &dH_L_dtime, dH_U_dX,
                                                     dH_L_dX, &dH_U_dp, &dH_U_ddh, time, delta_t);
    } else if (pd->e[pg->imtrx][R_SHELL_FILMP]) {
      local_post[LUB_HEIGHT] = fv->sh_fh;
    } else if (pd->e[pg->imtrx][R_TFMP_MASS] || pd->e[pg->imtrx][R_TFMP_BOUND]) {
      double tt = 1.0;
      GAP_STRUCT gap_v;
      GAP_STRUCT *gap = &gap_v;
      gap->time = time;
      gap->tt = tt;
      gap->delta_t = delta_t;
      gap->n_dof = n_dof;
      gap->dof_map = dof_map;
      load_gap_model(gap);

      local_post[LUB_HEIGHT] = gap->h;
    }

    switch (mp->FSIModel) {
    case FSI_MESH_CONTINUUM:
    case FSI_MESH_UNDEF:
    case FSI_SHELL_ONLY_UNDEF:
      for (a = 0; a < dim; a++) {
        local_post[LUB_HEIGHT] -= fv->snormal[a] * fv->d[a];
      }
      break;
    case FSI_SHELL_ONLY_MESH:
      if (((pd->e[pg->imtrx][R_SHELL_NORMAL1]) && (pd->e[pg->imtrx][R_SHELL_NORMAL2]) &&
           (pd->e[pg->imtrx][R_SHELL_NORMAL3])) ||
          (pd->e[pg->imtrx][R_MESH1] && pd->e[pg->imtrx][R_SHELL_NORMAL1] &&
           (pd->e[pg->imtrx][R_SHELL_NORMAL2]) && pd->e[pg->imtrx][R_LUBP])) {
        double dh_dmesh[DIM][MDE];
        double dh_dnormal[DIM][MDE];
        double d2h_dtime_dmesh[DIM][MDE];
        double d2h_dtime_dnormal[DIM][MDE];
        double d_gradIIh_dmesh[DIM][DIM][MDE];
        double d_gradIIh_dnormal[DIM][DIM][MDE];

        // dh_dtime not used here
        double tt = 1.0;
        double dh_dtime = 1.0;
        double gradII_h[DIM];

        load_displacement_coupling_model(
            tt, delta_t, &(local_post[LUB_HEIGHT]), &dh_dtime, gradII_h, dh_dmesh, dh_dnormal,
            d2h_dtime_dmesh, d2h_dtime_dnormal, d_gradIIh_dmesh, d_gradIIh_dnormal, n_dof, dof_map);
      } else if (pd->e[pg->imtrx][R_TFMP_BOUND]) {
        // do nothing:
        // the normal_dot_displacement is taken care of in load_gap_model
      } else {
        for (a = 0; a < dim; a++) {
          local_post[LUB_HEIGHT] -= fv->snormal[a] * fv->d[a];
        }
      }
      break;
    case FSI_REALSOLID_CONTINUUM:
      for (a = 0; a < dim; a++) {
        local_post[LUB_HEIGHT] -= fv->snormal[a] * fv->d_rs[a];
      }
      break;
    }
    local_lumped[LUB_HEIGHT] = 1.0;

    /* Cleanup */
    fv->wt = wt;
    safe_free((void *)n_dof);

  } /* end of LUB_HEIGHT */

  if (LUB_HEIGHT_2 != -1 && (pd->e[pg->imtrx][R_LUBP_2])) {
    double H_U_2, dH_U_2_dtime, H_L_2, dH_L_2_dtime;
    double dH_U_2_dX[DIM], dH_L_2_dX[DIM], dH_U_2_dp, dH_U_2_ddh;

    /* Setup lubrication */
    int *n_dof = NULL;
    int dof_map[MDE];
    dbl wt = fv->wt;
    n_dof = (int *)array_alloc(1, MAX_VARIABLE_TYPES, sizeof(int));
    lubrication_shell_initialize(n_dof, dof_map, -1, xi, exo, 0);

    if (pd->e[pg->imtrx][R_LUBP_2]) {
      local_post[LUB_HEIGHT_2] =
          height_function_model(&H_U_2, &dH_U_2_dtime, &H_L_2, &dH_L_2_dtime, dH_U_2_dX, dH_L_2_dX,
                                &dH_U_2_dp, &dH_U_2_ddh, time, delta_t);
    } else if (pd->e[pg->imtrx][R_SHELL_FILMP]) {
      local_post[LUB_HEIGHT] = 0.;
    }

    switch (mp->FSIModel) {
    case FSI_MESH_CONTINUUM:
    case FSI_MESH_UNDEF:
    case FSI_SHELL_ONLY_UNDEF:
      for (a = 0; a < dim; a++) {
        local_post[LUB_HEIGHT_2] -= fv->snormal[a] * fv->d[a];
      }
      break;
    case FSI_SHELL_ONLY_MESH:
      if ((pd->e[pg->imtrx][R_SHELL_NORMAL1]) && (pd->e[pg->imtrx][R_SHELL_NORMAL2]) &&
          (pd->e[pg->imtrx][R_SHELL_NORMAL3])) {
        for (a = 0; a < dim; a++) {
          local_post[LUB_HEIGHT] -= fv->n[a] * fv->d[a];
        }
      } else {
        for (a = 0; a < dim; a++) {
          local_post[LUB_HEIGHT] -= fv->snormal[a] * fv->d[a];
        }
      }
      break;
    case FSI_REALSOLID_CONTINUUM:
      for (a = 0; a < dim; a++) {
        local_post[LUB_HEIGHT_2] -= fv->snormal[a] * fv->d_rs[a];
      }
      break;
    }
    local_lumped[LUB_HEIGHT_2] = 1.0;

    /* Cleanup */
    fv->wt = wt;
    safe_free((void *)n_dof);

  } /* end of LUB_HEIGHT_2 */

  if ((LUB_VELO_UPPER != -1 || LUB_VELO_LOWER != -1) && (pd->e[pg->imtrx][R_LUBP])) {
    /* Setup lubrication */
    int *n_dof = NULL;
    int dof_map[MDE];
    dbl wt = fv->wt;
    n_dof = (int *)array_alloc(1, MAX_VARIABLE_TYPES, sizeof(int));
    lubrication_shell_initialize(n_dof, dof_map, -1, xi, exo, 0);

    /* Post values */
    double veloU[DIM], veloL[DIM];
    velocity_function_model(veloU, veloL, time, delta_t);
    if (LUB_VELO_UPPER != -1) {
      local_post[LUB_VELO_UPPER] = veloU[0];
      local_lumped[LUB_VELO_UPPER] = 1.0;
      local_post[LUB_VELO_UPPER + 1] = veloU[1];
      local_lumped[LUB_VELO_UPPER + 1] = 1.0;
      local_post[LUB_VELO_UPPER + 2] = veloU[2];
      local_lumped[LUB_VELO_UPPER + 2] = 1.0;
    }
    if (LUB_VELO_LOWER != -1) {
      local_post[LUB_VELO_LOWER] = veloL[0];
      local_lumped[LUB_VELO_LOWER] = 1.0;
      local_post[LUB_VELO_LOWER + 1] = veloL[1];
      local_lumped[LUB_VELO_LOWER + 1] = 1.0;
      local_post[LUB_VELO_LOWER + 2] = veloL[2];
      local_lumped[LUB_VELO_LOWER + 2] = 1.0;
    }

    /* Cleanup */
    fv->wt = wt;
    safe_free((void *)n_dof);

  } /* end of LUB_VELO */

  if ((LUB_VELO_FIELD != -1) && (pd->e[pg->imtrx][R_LUBP] || pd->e[pg->imtrx][R_SHELL_FILMP])) {

    /* Setup lubrication */
    int *n_dof = NULL;
    int dof_map[MDE];
    n_dof = (int *)array_alloc(1, MAX_VARIABLE_TYPES, sizeof(int));
    lubrication_shell_initialize(n_dof, dof_map, -1, xi, exo, 0);

    /* Calculate velocities */
    if (pd->e[pg->imtrx][R_LUBP]) {
      calculate_lub_q_v(R_LUBP, time, delta_t, xi, exo);
    } else {
      calculate_lub_q_v(R_SHELL_FILMP, time, delta_t, xi, exo);
    }

    /* Post velocities */
    local_post[LUB_VELO_FIELD] = LubAux->v_avg[0];
    local_lumped[LUB_VELO_FIELD] = 1.0;
    local_post[LUB_VELO_FIELD + 1] = LubAux->v_avg[1];
    local_lumped[LUB_VELO_FIELD + 1] = 1.0;
    local_post[LUB_VELO_FIELD + 2] = LubAux->v_avg[2];
    local_lumped[LUB_VELO_FIELD + 2] = 1.0;

    /* Cleanup */
    safe_free((void *)n_dof);

  } /* end of LUB_VELO_FIELD */

  if ((LUB_VELO_FIELD_2 != -1) && (pd->e[pg->imtrx][R_LUBP_2])) {

    /* Setup lubrication */
    int *n_dof = NULL;
    int dof_map[MDE];
    n_dof = (int *)array_alloc(1, MAX_VARIABLE_TYPES, sizeof(int));
    lubrication_shell_initialize(n_dof, dof_map, -1, xi, exo, 0);

    /* Calculate velocities */
    calculate_lub_q_v(R_LUBP_2, time, delta_t, xi, exo);

    /* Post velocities */
    local_post[LUB_VELO_FIELD_2] = LubAux->v_avg[0];
    local_lumped[LUB_VELO_FIELD_2] = 1.0;
    local_post[LUB_VELO_FIELD_2 + 1] = LubAux->v_avg[1];
    local_lumped[LUB_VELO_FIELD_2 + 1] = 1.0;
    local_post[LUB_VELO_FIELD_2 + 2] = LubAux->v_avg[2];
    local_lumped[LUB_VELO_FIELD_2 + 2] = 1.0;

    /* Cleanup */
    safe_free((void *)n_dof);

  } /* end of LUB_VELO_FIELD_2 */

  if ((LUB_FLUID_SOURCE != -1) &&
      ((pd->e[pg->imtrx][R_LUBP]) ||
       (pd->e[pg->imtrx][R_SHELL_FILMP] && pd->e[pg->imtrx][R_SHELL_FILMH]))) {
    /* Setup lubrication */
    int *n_dof = NULL;
    int dof_map[MDE];
    dbl wt = fv->wt;
    n_dof = (int *)array_alloc(1, MAX_VARIABLE_TYPES, sizeof(int));
    lubrication_shell_initialize(n_dof, dof_map, -1, xi, exo, 0);

    /* Post values */
    double LubSourceFlux = 0.0;
    lubrication_fluid_source(&LubSourceFlux, NULL, NULL);

    local_post[LUB_FLUID_SOURCE] = LubSourceFlux;
    local_lumped[LUB_FLUID_SOURCE] = 1.0;

    /* Cleanup */
    fv->wt = wt;
    safe_free((void *)n_dof);

  } /* end of LUB_FLUID_SOURCE */

  if ((PP_LAME_MU != -1) && (pd->e[pg->imtrx][R_MESH1])) {

    /* Define parameters */
    double mu;
    double lambda;
    double thermexp = 0;
    double speciesexp[MAX_CONC];
    double viscos = 0, dil_viscos = 0;
    double d_mu_dx[DIM][MDE];
    double d_lambda_dx[DIM][MDE];
    double d_thermexp_dx[MAX_VARIABLE_TYPES + MAX_CONC];
    double d_speciesexp_dx[MAX_CONC][MAX_VARIABLE_TYPES + MAX_CONC];
    double d_viscos_dx[MAX_VARIABLE_TYPES + MAX_CONC];
    double d_dilviscos_dx[MAX_VARIABLE_TYPES + MAX_CONC];

    /* Calculate modulus */

    if (pd->MeshMotion == TOTAL_ALE) {
      load_elastic_properties(elc_rs, &mu, &lambda, &thermexp, speciesexp, &viscos, &dil_viscos,
                              d_mu_dx, d_lambda_dx, d_thermexp_dx, d_speciesexp_dx, d_viscos_dx,
                              d_dilviscos_dx);
    } else {
      load_elastic_properties(elc, &mu, &lambda, &thermexp, speciesexp, &viscos, &dil_viscos,
                              d_mu_dx, d_lambda_dx, d_thermexp_dx, d_speciesexp_dx, d_viscos_dx,
                              d_dilviscos_dx);
    }

    /* Post velocities */
    local_post[PP_LAME_MU] = mu;
    local_lumped[PP_LAME_MU] = 1.0;

  } /* end of PP_LAME_MU */

  if ((PP_LAME_LAMBDA != -1) && (pd->e[pg->imtrx][R_MESH1])) {

    /* Define parameters */
    double mu;
    double lambda;
    double thermexp = 0;
    double speciesexp[MAX_CONC];
    double viscos = 0, dil_viscos = 0;
    double d_mu_dx[DIM][MDE];
    double d_lambda_dx[DIM][MDE];
    double d_thermexp_dx[MAX_VARIABLE_TYPES + MAX_CONC];
    double d_speciesexp_dx[MAX_CONC][MAX_VARIABLE_TYPES + MAX_CONC];
    double d_viscos_dx[MAX_VARIABLE_TYPES + MAX_CONC];
    double d_dilviscos_dx[MAX_VARIABLE_TYPES + MAX_CONC];

    /* Calculate modulus */
    if (pd->MeshMotion == TOTAL_ALE) {
      load_elastic_properties(elc_rs, &mu, &lambda, &thermexp, speciesexp, &viscos, &dil_viscos,
                              d_mu_dx, d_lambda_dx, d_thermexp_dx, d_speciesexp_dx, d_viscos_dx,
                              d_dilviscos_dx);
    } else {
      load_elastic_properties(elc, &mu, &lambda, &thermexp, speciesexp, &viscos, &dil_viscos,
                              d_mu_dx, d_lambda_dx, d_thermexp_dx, d_speciesexp_dx, d_viscos_dx,
                              d_dilviscos_dx);
    }

    /* Post velocities */
    local_post[PP_LAME_LAMBDA] = lambda;
    local_lumped[PP_LAME_LAMBDA] = 1.0;

  } /* end of PP_LAME_LAMBDA */

  if (DISJ_PRESS != -1 && (pd->e[pg->imtrx][R_SHELL_FILMP])) {

    /* Setup lubrication */
    int *n_dof = NULL;
    int dof_map[MDE];
    dbl wt = fv->wt;
    n_dof = (int *)array_alloc(1, MAX_VARIABLE_TYPES, sizeof(int));
    lubrication_shell_initialize(n_dof, dof_map, -1, xi, exo, 0);

    double DisjPress;
    double grad_DisjPress[DIM];
    double dgrad_DisjPress_dH1[DIM][MDE], dgrad_DisjPress_dH2[DIM][MDE];
    double dgrad_DisjPress_dH[DIM][MDE];

    DisjPress =
        disjoining_pressure_model(fv->sh_fh, fv->grad_sh_fh, n_dof, dof_map, grad_DisjPress,
                                  dgrad_DisjPress_dH1, dgrad_DisjPress_dH2, dgrad_DisjPress_dH);

    local_post[DISJ_PRESS] = DisjPress;
    local_lumped[DISJ_PRESS] = 1.0;

    /* Cleanup */
    fv->wt = wt;
    safe_free((void *)n_dof);

  } /* end of DISJ_PRESS */

  if ((SH_SAT_OPEN != -1) && pd->e[pg->imtrx][R_SHELL_SAT_OPEN]) {

    /* Calculate saturation */
    Patm = mp->PorousShellPatm;
    cap_pres = Patm - fv->sh_p_open;
    S = load_saturation(mp->porosity, cap_pres, d_cap_pres);

    /* Post saturation */
    local_post[SH_SAT_OPEN] = S;
    local_lumped[SH_SAT_OPEN] = 1.0;
  } /* end of SH_SAT_OPEN */

  if ((SH_SAT_OPEN_2 != -1) && pd->e[pg->imtrx][R_SHELL_SAT_OPEN_2]) {

    /* Calculate saturation */
    Patm = mp->PorousShellPatm;
    cap_pres = Patm - fv->sh_p_open_2;
    S = load_saturation(mp->porosity, cap_pres, d_cap_pres);

    /* Post saturation */
    local_post[SH_SAT_OPEN_2] = S;
    local_lumped[SH_SAT_OPEN_2] = 1.0;
  } /* end of SH_SAT_OPEN_2 */

  if ((SH_CAP_PRES != -1) && pd->e[pg->imtrx][R_SHELL_SAT_1]) {

    double shell_cap_pressure_nodal[MDE] = {0.0};
    double shell_cap_pressure = 0.0;

    for (j = 0; j < ei[pg->imtrx]->dof[SHELL_SAT_1]; j++) {
      shell_cap_pressure_nodal[j] = load_cap_pres(0, j, -1, *esp->sh_sat_1[j]);
      shell_cap_pressure += shell_cap_pressure_nodal[j] * bf[SHELL_SAT_1]->phi[j];
    }

    /* Post capillary pressure */
    local_post[SH_CAP_PRES] = shell_cap_pressure;
    local_lumped[SH_CAP_PRES] = 1.0;

    if (pd->e[pg->imtrx][R_SHELL_SAT_2]) {
      shell_cap_pressure = 0.0;
      for (j = 0; j < ei[pg->imtrx]->dof[SHELL_SAT_2]; j++) {
        shell_cap_pressure_nodal[j] = load_cap_pres(1, j, -1, *esp->sh_sat_2[j]);
        shell_cap_pressure += shell_cap_pressure_nodal[j] * bf[SHELL_SAT_2]->phi[j];
      }

      /* Post capillary pressure */
      local_post[SH_CAP_PRES + 1] = shell_cap_pressure;
      local_lumped[SH_CAP_PRES + 1] = 1.0;
    }

    if (pd->e[pg->imtrx][R_SHELL_SAT_3]) {
      shell_cap_pressure = 0.0;
      for (j = 0; j < ei[pg->imtrx]->dof[SHELL_SAT_3]; j++) {
        shell_cap_pressure_nodal[j] = load_cap_pres(2, j, -1, *esp->sh_sat_3[j]);
        shell_cap_pressure += shell_cap_pressure_nodal[j] * bf[SHELL_SAT_3]->phi[j];
      }

      /* Post capillary pressure */
      local_post[SH_CAP_PRES + 2] = shell_cap_pressure;
      local_lumped[SH_CAP_PRES + 2] = 1.0;
    }

  } /* end of SH_CAP_PRES */

  if ((SH_PORE_FLUX != -1) && pd->e[pg->imtrx][R_SHELL_SAT_1] && pd->e[pg->imtrx][R_SHELL_SAT_2]) {

    double j_1_2_nodal[MDE] = {0.0}, j_2_3_nodal[MDE] = {0.0};
    double j_1_2 = 0.0, j_2_3 = 0.0;

    /* Calculate the fluxes and their sensitivities */
    porous_shell_open_source_model(j_1_2_nodal, j_2_3_nodal, NULL, NULL);

    for (j = 0; j < ei[pg->imtrx]->dof[SHELL_SAT_3]; j++) {
      j_1_2 += j_1_2_nodal[j] * bf[SHELL_SAT_1]->phi[j];
      j_2_3 += j_2_3_nodal[j] * bf[SHELL_SAT_3]->phi[j];
    }

    /* Post flux between porous shell layers 1 and 2 */
    local_post[SH_PORE_FLUX] = j_1_2;
    local_lumped[SH_PORE_FLUX] = 1.0;

    if (pd->e[pg->imtrx][R_SHELL_SAT_3]) {
      /* Post flux between porous shell layers 2 and 3 */
      local_post[SH_PORE_FLUX + 1] = j_2_3;
      local_lumped[SH_PORE_FLUX + 1] = 1.0;
    }
  } /* end of SH_PORE_FLUX */

  if (REL_LIQ_PERM != -1 && (mp->PorousMediaType == POROUS_UNSATURATED ||
                             mp->PorousMediaType == POROUS_SHELL_UNSATURATED ||
                             mp->PorousMediaType == POROUS_TWO_PHASE)) {

    /* Continuum porous media */
    if (pd->e[pg->imtrx][R_POR_LIQ_PRES]) {
      local_post[REL_LIQ_PERM] = mp->rel_liq_perm;
      local_lumped[REL_LIQ_PERM] = 1.;
    }

    /* Shell porous media */
    else if ((pd->e[pg->imtrx][R_SHELL_SAT_OPEN]) || (pd->e[pg->imtrx][R_SHELL_SAT_OPEN_2])) {

      /* Calculate saturation */
      Patm = mp->PorousShellPatm;
      cap_pres = Patm - fv->sh_p_open;
      if (pd->e[pg->imtrx][R_SHELL_SAT_OPEN_2])
        cap_pres = Patm - fv->sh_p_open_2;
      S = load_saturation(mp->porosity, cap_pres, d_cap_pres);

      /* Then get relative permeability */
      if (mp->RelLiqPermModel != CONSTANT) {
        load_liq_perm(mp->porosity, cap_pres, mp->saturation, d_cap_pres);
      }

      local_post[REL_LIQ_PERM] = mp->rel_liq_perm;
      local_lumped[REL_LIQ_PERM] = 1.;
    }
  } /* end of REL_LIQ_PERM */

  if ((SH_STRESS_TENSOR != -1) && pd->e[pg->imtrx][R_SHELL_NORMAL1] &&
      pd->e[pg->imtrx][R_SHELL_NORMAL2] && pd->e[pg->imtrx][R_SHELL_NORMAL3] &&
      pd->e[pg->imtrx][R_MESH1]) {

    dbl TT[DIM][DIM];
    dbl dTT_dx[DIM][DIM][DIM][MDE];
    dbl dTT_dnormal[DIM][DIM][DIM][MDE];

    memset(TT, 0.0, sizeof(double) * DIM * DIM);
    memset(dTT_dx, 0.0, sizeof(double) * DIM * DIM * DIM * MDE);
    memset(dTT_dnormal, 0.0, sizeof(double) * DIM * DIM * DIM * MDE);

    shell_stress_tensor(TT, dTT_dx, dTT_dnormal);

    /* Post stresses */
    local_post[SH_STRESS_TENSOR] = TT[0][0];
    local_lumped[SH_STRESS_TENSOR] = 1.0;
    local_post[SH_STRESS_TENSOR + 1] = TT[1][1];
    local_lumped[SH_STRESS_TENSOR + 1] = 1.0;
    local_post[SH_STRESS_TENSOR + 2] = TT[0][1];
    local_lumped[SH_STRESS_TENSOR + 2] = 1.0;
  }

  if ((SH_TANG != -1) && pd->e[pg->imtrx][R_SHELL_NORMAL1] && pd->e[pg->imtrx][R_SHELL_NORMAL2] &&
      pd->e[pg->imtrx][R_SHELL_NORMAL3] && pd->e[pg->imtrx][R_MESH1]) {

    dbl t0[DIM];
    dbl t1[DIM];

    shell_tangents(t0, t1, NULL, NULL, NULL, NULL);

    /* Post tangents and curvatures */

    local_post[SH_TANG] = t0[0];
    local_lumped[SH_TANG] = 1.0;
    local_post[SH_TANG + 1] = t0[1];
    local_lumped[SH_TANG + 1] = 1.0;
    local_post[SH_TANG + 2] = t0[2];
    local_lumped[SH_TANG + 2] = 1.0;

    local_post[SH_TANG + 3] = t1[0];
    local_lumped[SH_TANG + 3] = 1.0;
    local_post[SH_TANG + 4] = t1[1];
    local_lumped[SH_TANG + 4] = 1.0;
    local_post[SH_TANG + 5] = t1[2];
    local_lumped[SH_TANG + 5] = 1.0;
  }

  if (VON_MISES_STRAIN != -1 && pd->e[pg->imtrx][R_MESH1]) {

    dbl INV, d_INV_dT[DIM][DIM];
    INV = calc_tensor_invariant(fv->strain, d_INV_dT, 4);
    local_post[VON_MISES_STRAIN] = 2.0 / 3.0 * INV;
    local_lumped[VON_MISES_STRAIN] = 1.;
  }

  if (VON_MISES_STRESS != -1 && pd->e[pg->imtrx][R_MESH1]) {

    dbl INV, d_INV_dT[DIM][DIM];

    // Calculate base stress tensor
    err = mesh_stress_tensor(TT, dTT_dx, dTT_dp, dTT_dc, dTT_dp_liq, dTT_dp_gas, dTT_dporosity,
                             dTT_dsink_mass, dTT_dT, dTT_dmax_strain, dTT_dcur_strain, elc->lame_mu,
                             elc->lame_lambda, delta_t, ielem, ip, ip_total);

    // Adjust tensor for various model cases
    if (cr->MeshFluxModel == LINEAR) {
      if (dim == 2) {
        TT[2][2] = 1.;
        TT[1][2] = 0.;
        TT[0][2] = 0.;
      }
    } else {
      if (dim == 2) {
        if (cr->MeshMotion == ARBITRARY) {
          TT[2][2] = (1. - fv->volume_change) * elc->lame_mu;
        } else {
          if (cr->MeshFluxModel == NONLINEAR || cr->MeshFluxModel == HOOKEAN_PSTRAIN ||
              cr->MeshFluxModel == INCOMP_PSTRAIN)
            TT[2][2] = (1. - pow(fv->volume_change, 2. / 3.)) * elc->lame_mu - fv->P;
          else
            TT[2][2] = 0.;
        }
        TT[1][2] = 0.;
        TT[0][2] = 0.;
      }
    }

    // Calculate invariant
    INV = calc_tensor_invariant(TT, d_INV_dT, 4);
    local_post[VON_MISES_STRESS] = INV;
    local_lumped[VON_MISES_STRESS] = 1.;
  }

  if (CONF_MAP != -1 && pd->v[pg->imtrx][POLYMER_STRESS11] &&
      (vn->evssModel == LOG_CONF || vn->evssModel == LOG_CONF_GRADV || vn->evssModel == CONF ||
       vn->evssModel == LOG_CONF_TRANSIENT || vn->evssModel == LOG_CONF_TRANSIENT_GRADV)) {
    index = 0;
    VISCOSITY_DEPENDENCE_STRUCT d_mup_struct;
    VISCOSITY_DEPENDENCE_STRUCT *d_mup = &d_mup_struct;
    d_mup = NULL;
    double lambda;
    double R1[DIM][DIM];
    double eig_values[DIM];
    dbl exp_s[DIM][DIM];
    for (mode = 0; mode < vn->modes; mode++) {
      if (vn->evssModel == CONF) {

        for (a = 0; a < VIM; a++) {
          for (b = 0; b < VIM; b++) {
            exp_s[a][b] = fv->S[mode][a][b];
          }
        }
      } else {
#ifdef ANALEIG_PLEASE
        analytical_exp_s(fv->S[mode], exp_s, eig_values, R1, NULL);
#else
        compute_exp_s(fv->S[mode], exp_s, eig_values, R1);
#endif
      }
      mup = viscosity(ve[mode]->gn, gamma, d_mup);
      // Polymer time constant
      lambda = 0.0;
      if (ve[mode]->time_constModel == CONSTANT) {
        lambda = ve[mode]->time_const;
      }
      /* Looks like these models are not working right now
       *else if(ve[mode]->time_constModel == CARREAU || ve[mode]->time_constModel == POWER_LAW)
       * {
       *   lambda = mup/ve[mode]->time_const;
       * }
       */
      if (lambda == 0.0) {
        GOMA_EH(GOMA_ERROR, "The conformation tensor needs a non-zero polymer time constant.");
      }
      if (mup == 0.0) {
        GOMA_EH(GOMA_ERROR, "The conformation tensor needs a non-zero polymer viscosity.");
      }

      if (pd->v[pg->imtrx][v_s[mode][0][0]]) {
        for (a = 0; a < VIM; a++) {
          for (b = 0; b < VIM; b++) {
            /* since the stress tensor is symmetric,
               only assemble the upper half */
            if (a <= b) {
              if (pd->v[pg->imtrx][v_s[mode][a][b]]) {
                local_post[CONF_MAP + index] = (mup / lambda) * (exp_s[a][b] - delta(a, b));
                local_lumped[CONF_MAP + index] = 1.;
                index++;
              }
            }
          } // for b
        }   // for a
      }
    } // Loop over modes
  }

  if (CONF_MAP != -1 && pd->v[pg->imtrx][POLYMER_STRESS11] && (vn->evssModel == SQRT_CONF)) {
    index = 0;
    VISCOSITY_DEPENDENCE_STRUCT d_mup_struct;
    VISCOSITY_DEPENDENCE_STRUCT *d_mup = &d_mup_struct;
    d_mup = NULL;
    double lambda;
    for (mode = 0; mode < vn->modes; mode++) {
      mup = viscosity(ve[mode]->gn, gamma, d_mup);

      dbl b_dot_b[DIM][DIM];
      dbl cb[DIM][DIM];
      for (int ii = 0; ii < VIM; ii++) {
        for (int jj = 0; jj < VIM; jj++) {
          if (ii <= jj) {
            cb[ii][jj] = fv->S[mode][ii][jj];
            cb[jj][ii] = cb[ii][jj];
          }
        }
      }

      tensor_dot(cb, cb, b_dot_b, VIM);
      // Polymer time constant
      lambda = 0.0;
      if (ve[mode]->time_constModel == CONSTANT) {
        lambda = ve[mode]->time_const;
      }
      /* Looks like these models are not working right now
       *else if(ve[mode]->time_constModel == CARREAU || ve[mode]->time_constModel == POWER_LAW)
       * {
       *   lambda = mup/ve[mode]->time_const;
       * }
       */
      if (lambda == 0.0) {
        GOMA_EH(GOMA_ERROR, "The conformation tensor needs a non-zero polymer time constant.");
      }
      if (mup == 0.0) {
        GOMA_EH(GOMA_ERROR, "The conformation tensor needs a non-zero polymer viscosity.");
      }

      if (pd->v[pg->imtrx][v_s[mode][0][0]]) {
        for (a = 0; a < VIM; a++) {
          for (b = 0; b < VIM; b++) {
            /* since the stress tensor is symmetric,
               only assemble the upper half */
            if (a <= b) {
              if (pd->v[pg->imtrx][v_s[mode][a][b]]) {
                local_post[CONF_MAP + index] = -(mup / lambda) * (delta(a, b) - b_dot_b[a][b]);
                local_lumped[CONF_MAP + index] = 1.;
                index++;
              }
            }
          } // for b
        }   // for a
      }
    } // Loop over modes
  }
  /*if (cr->MassFluxModel == DM_SUSPENSION_BALANCE ) {
    index = 0;
    int w = 0;
    for (a = 0; a < dim; a++)
      {
        hs[a] = 0.;
      }

    err = get_continuous_species_terms(&s_terms, time, theta, delta_t, hs);
    for (a=0; a < dim; a++)
      {
        local_post[J_FLUX + index] = s_terms.diff_flux[w][a];
        local_lumped[J_FLUX + index] = 1.;
        index++;
      }
      }*/

  if (GRAD_SH != -1 && pd->v[pg->imtrx][SHEAR_RATE]) {
    index = 0;
    for (a = 0; a < dim; a++) {
      local_post[GRAD_SH + index] = fv->grad_SH[a];
      local_lumped[GRAD_SH + index] = 1.;
      index++;
    }
  }

  if (GRAD_Y != -1 && pd->v[pg->imtrx][MASS_FRACTION]) {
    index = 0;
    for (a = 0; a < dim; a++) {
      local_post[GRAD_Y + index] = fv->grad_c[0][a];
      local_lumped[GRAD_Y + index] = 1.;
      index++;
    }
  }

  if (EIG1 != -1) {
    index = 0;
    dbl v1[DIM];
    dbl v2[DIM], v3[DIM];
    dbl e1 = 0., e2 = 0., e3 = 0.;
    dbl gamma_dot[DIM][DIM], eigen[3];
    dbl v_bias[DIM], vy_bias[DIM], vort_bias[DIM];

    memset(v1, 0, DIM * sizeof(dbl));
    memset(v2, 0, DIM * sizeof(dbl));
    memset(v3, 0, DIM * sizeof(dbl));
    memset(v_bias, 0, DIM * sizeof(dbl));
    memset(vort_bias, 0, DIM * sizeof(dbl));
    memset(vy_bias, 0, DIM * sizeof(dbl));
    memset(gamma_dot, 0, DIM * DIM * sizeof(dbl));

    for (a = 0; a < dim; a++) {
      for (b = 0; b < dim; b++) {
        gamma_dot[a][b] = fv->grad_v[a][b] + fv->grad_v[b][a];
      }
    }

    find_eigenvalues_eigenvectors(gamma_dot, &e1, &e2, &e3, v1, v2, v3);
    eigen[0] = e1;
    eigen[1] = e2;
    eigen[2] = e3;

    v_bias[0] = 1.;
    vy_bias[1] = 1.;
    vort_bias[2] = 1.;

    bias_eigenvector_to(v1, v_bias);
    bias_eigenvector_to(v2, vy_bias);
    bias_eigenvector_to(v3, vort_bias);

    for (a = 0; a < DIM; a++) {
      local_post[EIG + index] = eigen[a];
      local_lumped[EIG + index] = 1.;
      index++;
    }

    index = 0;
    for (a = 0; a < DIM; a++) {
      local_post[EIG1 + index] = v1[a];
      local_lumped[EIG1 + index] = 1.;
      index++;
    }

    index = 0;
    for (a = 0; a < DIM; a++) {
      local_post[EIG2 + index] = v2[a];
      local_lumped[EIG2 + index] = 1.;
      index++;
    }

    index = 0;
    for (a = 0; a < DIM; a++) {
      local_post[EIG3 + index] = v3[a];
      local_lumped[EIG3 + index] = 1.;
      index++;
    }
  }

  if (USER_POST != -1) {
    /* calculate a user-specified post-processing variable */

    err = get_continuous_species_terms(&s_terms, time, theta, delta_t, hs);
    local_post[USER_POST] = user_post(u_post_proc);
    local_lumped[USER_POST] = 1.;
  }

  /*
   *  ---- NOTE do this step above this line, order of post processing
   *       variables is not important here, anymore ---
   *
   *  Post-processing Step 4: add algorithm to calculate your new variable
   *                          in mm_post_proc
   *                          and put it in post_proc_vect array
   */

  wt = fv->wt;
  eqn = pd->ProjectionVar;
  det_J = bf[eqn]->detJ;

  if (ielem_type == BILINEAR_SHELL || ielem_type == BIQUAD_SHELL ||
      ielem_type == BILINEAR_TRISHELL || ielem_type == P0_SHELL || ielem_type == P1_SHELL) {
    int *n_dof = NULL;
    int dof_map[MDE];
    n_dof = (int *)array_alloc(1, MAX_VARIABLE_TYPES, sizeof(int));
    lubrication_shell_initialize(n_dof, dof_map, -1, xi, exo, 0);
    safe_free((void *)n_dof);
    det_J = fv->sdet;
  }

  /*
   * Put local element contributions into the global right-hand side
   *
   * HKM NOTE on projection operations:
   *
   *    One could make the argument that the procedure below is not
   * the most accurate one possible. A more accurate procedure would
   * be to form completely the local element contribution to the
   * projection operation:
   *
   *     MassMatrix_lem_i_j * S_j = integral( s_gp * bf_i)
   *
   *  where MassMatrix_lem_i_j = local element mass matrix
   *                           = int( bf_i * bf_j)
   *
   *  Then, invert that equation to obtain S_j. Then, renormalize
   * the result to get the elemental contribution:
   *
   *     S_i_elem  = S_i * sum_j(MassMatrix_lem_i_j)
   *     lumpedMM  = sum_j(MassMatrix_lem_i_j)
   *
   * And then, add up the elemental contributions.
   *
   * The method described above gives an exact solution to a sample
   * problem. The lumped Mass Matrix method carried out below
   * does not. However, in most papers on projection operations,
   * (aka Hughes et al.) the mass lumping operation undertaken
   * below is deemed to be sufficient.
   */
  for (i = 0; i < ei[pg->imtrx]->num_local_nodes; i++) {
    I = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + i];
    /*
     *  check to make sure that unknowns are defined at this node,
     *  otherwise don't add anything to this node
     */
    ldof = bf[pd->ProjectionVar]->interpolation == I_N1
               ? i
               : ei[upd->matrix_index[pd->ProjectionVar]]->ln_to_dof[eqn][i];
    if (ldof >= 0) {
      phi_i = bf[eqn]->phi[ldof];
      for (var = 0; var < rd->TotalNVPostOutput; var++) {
        post_proc_vect[var][I] += local_post[var] * phi_i * wt * det_J;
        lumped_mass[var][I] += local_lumped[var] * phi_i * wt * det_J;
      }
    } else {
      /*
       *  Special Compatibility Section for the case where the projection
       *  basis functions don't span all of the nodes in the local element.
       *
       */
      midside = 0;
      if (bf[eqn]->interpolation == I_Q1 || bf[eqn]->interpolation == I_Q1_D ||
          bf[eqn]->interpolation == I_Q1_XV || bf[eqn]->interpolation == I_Q1_G ||
          bf[eqn]->interpolation == I_Q1_GP || bf[eqn]->interpolation == I_Q1_GN ||
          bf[eqn]->interpolation == I_SP) {
        if (ielem_type == S_BIQUAD_QUAD || ielem_type == BIQUAD_QUAD) {
          midside = 1;
        }
      }
      if (!midside) {
        GOMA_WH(-1, "calc_standard_fields: Couldn't find eqn whose interp spanned all local nodes");
      }
      if (i > 3 && i < 8) {
        ileft = i - 4;
        iright = i - 3;
        if (iright == 4)
          iright = 0;
        ldof = ei[pg->imtrx]->ln_to_dof[eqn][ileft];
        ldof_right = ei[pg->imtrx]->ln_to_dof[eqn][iright];
        phi_i = bf[eqn]->phi[ldof];
        phi_j = bf[eqn]->phi[ldof_right];
        for (var = 0; var < rd->TotalNVPostOutput; var++) {
          post_proc_vect[var][I] += local_post[var] * (phi_i + phi_j) * wt * det_J;
          lumped_mass[var][I] += local_lumped[var] * (phi_i + phi_j) * wt * det_J;
        }
      }
      if (i == 8 && ielem_type == BIQUAD_QUAD) {
        for (ileft = 0; ileft < 4; ileft++) {
          ldof = ei[pg->imtrx]->ln_to_dof[eqn][ileft];
          phi_i = bf[eqn]->phi[ldof];
          for (var = 0; var < rd->TotalNVPostOutput; var++) {
            post_proc_vect[var][I] += local_post[var] * phi_i * wt * det_J;
            lumped_mass[var][I] += local_lumped[var] * phi_i * wt * det_J;
          }
        }
      }
    }
  }

  safer_free((void **)&local_post);
  safer_free((void **)&local_lumped);
  return (status);
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void post_process_average(double x[],            /* Solution vector for the current processor */
                          double x_old[],        /* Solution vector at last time step */
                          double xdot[],         /* time derivative of solution vector */
                          double xdot_old[],     /* dx/dt at previous time step */
                          double resid_vector[], /* Residual vector for the
                                                    current proc */
                          Exo_DB *exo,
                          Dpi *dpi,
                          double **post_proc_vect,
                          double time) {
  int ielem;
  int ebn;
  int err;
  int e_start = exo->eb_ptr[0];
  int e_end = exo->eb_ptr[exo->num_elem_blocks];
  double xi[DIM];

  double **avg_count;
  double **avg_sum;

  if (nn_average <= 0) {
    return;
  }

  avg_count = calloc(nn_average, sizeof(double *));
  avg_sum = calloc(nn_average, sizeof(double *));
  for (int i = 0; i < nn_average; i++) {
    avg_count[i] = calloc(dpi->num_universe_nodes, sizeof(double));
    avg_sum[i] = calloc(dpi->num_universe_nodes, sizeof(double));
  }

  for (ielem = e_start, ebn = 0; ielem < e_end; ielem++) {
    ebn = find_elemblock_index(ielem, exo);
    int mn = Matilda[ebn];
    if (mn < 0) {
      continue;
    }

    pd = pd_glob[mn];
    cr = cr_glob[mn];
    elc = elc_glob[mn];
    elc_rs = elc_rs_glob[mn];
    gn = gn_glob[mn];
    mp = mp_glob[mn];
    vn = vn_glob[mn];
    evpl = evpl_glob[mn];

    err = load_elem_dofptr(ielem, exo, x, x_old, xdot, xdot_old, 0);
    GOMA_EH(err, "load_elem_dofptr");

    err = bf_mp_init(pd);
    GOMA_EH(err, "bf_mp_init");

    fv->wt = 1.e30;

    int node;
    for (node = 0; node < ei[pg->imtrx]->num_local_nodes; node++) {

      int global_node = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + node];

      find_nodal_stu(node, ei[pg->imtrx]->ielem_type, xi, xi + 1, xi + 2);

      err = load_basis_functions(xi, bfd);
      GOMA_EH(err, "problem from load_basis_functions");

      err = beer_belly();
      GOMA_EH(err, "beer_belly");

      err = load_fv();
      GOMA_EH(err, "load_fv");

      err = load_bf_grad();
      GOMA_EH(err, "load_bf_grad");

      err = load_fv_vector();

      if (ei[pg->imtrx]->deforming_mesh &&
          (pd->e[pg->imtrx][R_MESH1] || pd->v[pg->imtrx][R_MESH1])) {
        if (!pd->e[pg->imtrx][R_MESH1]) {
          // printf(" We are here\n");
        }
        err = load_bf_mesh_derivs();
        GOMA_EH(err, "load_bf_mesh_derivs");
      }

      err = load_fv_grads();
      GOMA_EH(err, "load_fv_grads");

      if (ei[pg->imtrx]->deforming_mesh &&
          (pd->e[pg->imtrx][R_MESH1] || pd->v[pg->imtrx][R_MESH1])) {
        if (!pd->e[pg->imtrx][R_MESH1]) {
          //	printf(" We are here2\n");
        }
        err = load_fv_mesh_derivs(0);
        GOMA_EH(err, "load_fv_mesh_derivs");
      }

      if (mp->PorousMediaType != CONTINUOUS) {
        err = load_porous_properties();
        GOMA_EH(err, "load_porous_properties");
      }

      computeCommonMaterialProps_gp(time);

      sum_average_nodal(avg_count, avg_sum, global_node, time);
    }
  }

  for (int i = 0; i < nn_average; i++) {
    int node;
    for (node = 0; node < dpi->num_universe_nodes; node++) {
      if ((pp_average[i]->non_variable_type && pg->imtrx == 0) ||
          (!pp_average[i]->non_variable_type && Num_Var_In_Type[pg->imtrx][pp_average[i]->type])) {
        post_proc_vect[pp_average[i]->index_post][node] = avg_sum[i][node] / avg_count[i][node];
      }
    }
  }

  for (int i = 0; i < nn_average; i++) {
    free(avg_count[i]);
    free(avg_sum[i]);
  }

  free(avg_count);
  free(avg_sum);
}

void sum_average_nodal(double **avg_count, double **avg_sum, int global_node, double time) {
  for (int i = 0; i < nn_average; i++) {
    avg_count[i][global_node] += 1;
    if (pp_average[i]->non_variable_type == 0) {
      switch (pp_average[i]->type) {
      case PRESSURE:
        avg_sum[i][global_node] += fv->P;
        break;
      case TEMPERATURE:
        avg_sum[i][global_node] += fv->T;
        break;
      case MASS_FRACTION:
        avg_sum[i][global_node] += fv->c[pp_average[i]->species_index];
        break;
      case MOMENT0:
        avg_sum[i][global_node] += fv->moment[0];
        break;
      case MOMENT1:
        avg_sum[i][global_node] += fv->moment[1];
        break;
      case MOMENT2:
        avg_sum[i][global_node] += fv->moment[2];
        break;
      case MOMENT3:
        avg_sum[i][global_node] += fv->moment[3];
        break;
      case DENSITY_EQN: {
        double rho = density(NULL, time);
        avg_sum[i][global_node] += rho;
      } break;
      default:
        GOMA_EH(GOMA_ERROR, "Unknown nodal average type %d %s", pp_average[i]->type,
                pp_average[i]->type_name);
        break;
      }
    } else {
      switch (pp_average[i]->type) {
      case AVG_DENSITY: {
        double rho = density(NULL, time);
        avg_sum[i][global_node] += rho;
      } break;
      case AVG_HEAVISIDE: {
        load_lsi(ls->Length_Scale);
        avg_sum[i][global_node] += lsi->H;
      } break;
      case AVG_VISCOSITY: {
        double gamma[DIM][DIM];
        int a;
        int b;
        for (a = 0; a < VIM; a++) {
          for (b = 0; b < VIM; b++) {
            gamma[a][b] = fv->grad_v[a][b] + fv->grad_v[b][a];
          }
        }
        double mu = viscosity(gn, gamma, NULL);

        if (gn->ConstitutiveEquation == BINGHAM_MIXED) {

          dbl gamma[DIM][DIM];
          for (a = 0; a < VIM; a++) {
            for (b = 0; b < VIM; b++) {
              gamma[a][b] = fv->grad_v[a][b] + fv->grad_v[b][a];
            }
          }

          dbl gammadot;
          calc_shearrate(&gammadot, gamma, NULL, NULL);

          dbl Du_eps = sqrt(gammadot + gn->epsilon * gn->epsilon);

          mu += gn->tau_y / Du_eps;
        }

        avg_sum[i][global_node] += mu;
      } break;
      case AVG_SHEAR: {
        double gamma[DIM][DIM];
        int a;
        int b;
        for (a = 0; a < VIM; a++) {
          for (b = 0; b < VIM; b++) {
            gamma[a][b] = fv->grad_v[a][b] + fv->grad_v[b][a];
          }
        }
        dbl gammadot = 0;
        calc_shearrate(&gammadot, /* strain rate invariant */
                       gamma,     /* strain rate tensor */
                       NULL, NULL);

        avg_sum[i][global_node] += gammadot;
      } break;
      case AVG_EMR_X: {
        avg_sum[i][global_node] += fv->em_er[0];
      } break;
      case AVG_EMR_Y: {
        avg_sum[i][global_node] += fv->em_er[1];
      } break;
      case AVG_EMR_Z: {
        avg_sum[i][global_node] += fv->em_er[2];
      } break;
      case AVG_EMI_X: {
        avg_sum[i][global_node] += fv->em_ei[0];
      } break;
      case AVG_EMI_Y: {
        avg_sum[i][global_node] += fv->em_ei[1];
      } break;
      case AVG_EMI_Z: {
        avg_sum[i][global_node] += fv->em_ei[2];
      } break;
      case AVG_EM_INC_MAG: {
        const double c0 = 1.0 / (sqrt(upd->Free_Space_Permittivity * upd->Free_Space_Permeability));
        dbl x = fv->x[0];
        dbl y = fv->x[1];
        dbl z = fv->x[2];
        dbl freq = upd->EM_Frequency;
        dbl lambda0 = c0 / freq;
        dbl k0 = 2 * M_PI / lambda0;
        complex double wave[3] = {0};
        complex double curl_wave[3];
        if (mp->PermittivityModel != RADIAL_PML) {
          incident_wave(x, y, z, k0, wave, curl_wave);
        }
        avg_sum[i][global_node] +=
            sqrt(creal(wave[0]) * creal(wave[0]) + creal(wave[1]) * creal(wave[1]) +
                 creal(wave[2]) * creal(wave[2]) + cimag(wave[0]) * cimag(wave[0]) +
                 cimag(wave[1]) * cimag(wave[1]) + cimag(wave[2]) * cimag(wave[2]));
      } break;
      case AVG_EM_SCAT_MAG: {
        const double c0 = 1.0 / (sqrt(upd->Free_Space_Permittivity * upd->Free_Space_Permeability));
        dbl x = fv->x[0];
        dbl y = fv->x[1];
        dbl z = fv->x[2];
        dbl freq = upd->EM_Frequency;
        dbl lambda0 = c0 / freq;
        dbl k0 = 2 * M_PI / lambda0;
        complex double wave[3] = {0};
        complex double curl_wave[3];
        if (mp->PermittivityModel != RADIAL_PML) {
          incident_wave(x, y, z, k0, wave, curl_wave);
        }
        avg_sum[i][global_node] +=
            sqrt((fv->em_er[0] - creal(wave[0])) * (fv->em_er[0] - creal(wave[0])) +
                 (fv->em_er[1] - creal(wave[1])) * (fv->em_er[1] - creal(wave[1])) +
                 (fv->em_er[2] - creal(wave[2])) * (fv->em_er[2] - creal(wave[2])) +
                 (fv->em_ei[0] - cimag(wave[0])) * (fv->em_ei[0] - cimag(wave[0])) +
                 (fv->em_ei[1] - cimag(wave[1])) * (fv->em_ei[1] - cimag(wave[1])) +
                 (fv->em_ei[2] - cimag(wave[2])) * (fv->em_ei[2] - cimag(wave[2])));
      } break;
      case AVG_EM_MAG: {
        avg_sum[i][global_node] += sqrt(fv->em_er[0] * fv->em_er[0] + fv->em_er[1] * fv->em_er[1] +
                                        fv->em_er[2] * fv->em_er[2] + fv->em_ei[0] * fv->em_ei[0] +
                                        fv->em_ei[1] * fv->em_ei[1] + fv->em_ei[2] * fv->em_ei[2]);
      } break;
      case AVG_EMSCATR_X: {
        const double c0 = 1.0 / (sqrt(upd->Free_Space_Permittivity * upd->Free_Space_Permeability));
        dbl x = fv->x[0];
        dbl y = fv->x[1];
        dbl z = fv->x[2];
        dbl freq = upd->EM_Frequency;
        dbl lambda0 = c0 / freq;
        dbl k0 = 2 * M_PI / lambda0;
        complex double wave[3] = {0};
        complex double curl_wave[3];
        if (mp->PermittivityModel != RADIAL_PML) {
          incident_wave(x, y, z, k0, wave, curl_wave);
        }
        avg_sum[i][global_node] += fv->em_er[0] - creal(wave[0]);
      } break;
      case AVG_EMSCATR_Y: {
        const double c0 = 1.0 / (sqrt(upd->Free_Space_Permittivity * upd->Free_Space_Permeability));
        dbl x = fv->x[0];
        dbl y = fv->x[1];
        dbl z = fv->x[2];
        dbl freq = upd->EM_Frequency;
        dbl lambda0 = c0 / freq;
        dbl k0 = 2 * M_PI / lambda0;
        complex double wave[3] = {0};
        complex double curl_wave[3];
        if (mp->PermittivityModel != RADIAL_PML) {
          incident_wave(x, y, z, k0, wave, curl_wave);
        }
        avg_sum[i][global_node] += fv->em_er[1] - creal(wave[1]);
      } break;
      case AVG_EMSCATR_Z: {
        const double c0 = 1.0 / (sqrt(upd->Free_Space_Permittivity * upd->Free_Space_Permeability));
        dbl x = fv->x[0];
        dbl y = fv->x[1];
        dbl z = fv->x[2];
        dbl freq = upd->EM_Frequency;
        dbl lambda0 = c0 / freq;
        dbl k0 = 2 * M_PI / lambda0;
        complex double wave[3] = {0};
        complex double curl_wave[3];
        if (mp->PermittivityModel != RADIAL_PML) {
          incident_wave(x, y, z, k0, wave, curl_wave);
        }
        avg_sum[i][global_node] += fv->em_er[2] - creal(wave[2]);
      } break;
      case AVG_EMSCATI_X: {
        const double c0 = 1.0 / (sqrt(upd->Free_Space_Permittivity * upd->Free_Space_Permeability));
        dbl x = fv->x[0];
        dbl y = fv->x[1];
        dbl z = fv->x[2];
        dbl freq = upd->EM_Frequency;
        dbl lambda0 = c0 / freq;
        dbl k0 = 2 * M_PI / lambda0;
        complex double wave[3] = {0};
        complex double curl_wave[3];
        if (mp->PermittivityModel != RADIAL_PML) {
          incident_wave(x, y, z, k0, wave, curl_wave);
        }
        avg_sum[i][global_node] += fv->em_ei[0] - cimag(wave[0]);
      } break;
      case AVG_EMSCATI_Y: {
        const double c0 = 1.0 / (sqrt(upd->Free_Space_Permittivity * upd->Free_Space_Permeability));
        dbl x = fv->x[0];
        dbl y = fv->x[1];
        dbl z = fv->x[2];
        dbl freq = upd->EM_Frequency;
        dbl lambda0 = c0 / freq;
        dbl k0 = 2 * M_PI / lambda0;
        complex double wave[3] = {0};
        complex double curl_wave[3];
        if (mp->PermittivityModel != RADIAL_PML) {
          incident_wave(x, y, z, k0, wave, curl_wave);
        }
        avg_sum[i][global_node] += fv->em_ei[1] - creal(wave[1]);
      } break;
      case AVG_EMSCATI_Z: {
        const double c0 = 1.0 / (sqrt(upd->Free_Space_Permittivity * upd->Free_Space_Permeability));
        dbl x = fv->x[0];
        dbl y = fv->x[1];
        dbl z = fv->x[2];
        dbl freq = upd->EM_Frequency;
        dbl lambda0 = c0 / freq;
        dbl k0 = 2 * M_PI / lambda0;
        complex double wave[3] = {0};
        complex double curl_wave[3];
        if (mp->PermittivityModel != RADIAL_PML) {
          incident_wave(x, y, z, k0, wave, curl_wave);
        }
        avg_sum[i][global_node] += fv->em_ei[2] - cimag(wave[2]);
      } break;
      case AVG_EMINCR_X: {
        const double c0 = 1.0 / (sqrt(upd->Free_Space_Permittivity * upd->Free_Space_Permeability));
        dbl x = fv->x[0];
        dbl y = fv->x[1];
        dbl z = fv->x[2];
        dbl freq = upd->EM_Frequency;
        dbl lambda0 = c0 / freq;
        dbl k0 = 2 * M_PI / lambda0;
        complex double wave[3] = {0};
        complex double curl_wave[3];
        if (mp->PermittivityModel != RADIAL_PML) {
          incident_wave(x, y, z, k0, wave, curl_wave);
        }
        avg_sum[i][global_node] += creal(wave[0]);
      } break;
      case AVG_EMINCR_Y: {
        const double c0 = 1.0 / (sqrt(upd->Free_Space_Permittivity * upd->Free_Space_Permeability));
        dbl x = fv->x[0];
        dbl y = fv->x[1];
        dbl z = fv->x[2];
        dbl freq = upd->EM_Frequency;
        dbl lambda0 = c0 / freq;
        dbl k0 = 2 * M_PI / lambda0;
        complex double wave[3] = {0};
        complex double curl_wave[3];
        if (mp->PermittivityModel != RADIAL_PML) {
          incident_wave(x, y, z, k0, wave, curl_wave);
        }
        avg_sum[i][global_node] += creal(wave[1]);
      } break;
      case AVG_EMINCR_Z: {
        const double c0 = 1.0 / (sqrt(upd->Free_Space_Permittivity * upd->Free_Space_Permeability));
        dbl x = fv->x[0];
        dbl y = fv->x[1];
        dbl z = fv->x[2];
        dbl freq = upd->EM_Frequency;
        dbl lambda0 = c0 / freq;
        dbl k0 = 2 * M_PI / lambda0;
        complex double wave[3] = {0};
        complex double curl_wave[3];
        if (mp->PermittivityModel != RADIAL_PML) {
          incident_wave(x, y, z, k0, wave, curl_wave);
        }
        avg_sum[i][global_node] += creal(wave[2]);
      } break;
      case AVG_EMINCI_X: {
        const double c0 = 1.0 / (sqrt(upd->Free_Space_Permittivity * upd->Free_Space_Permeability));
        dbl x = fv->x[0];
        dbl y = fv->x[1];
        dbl z = fv->x[2];
        dbl freq = upd->EM_Frequency;
        dbl lambda0 = c0 / freq;
        dbl k0 = 2 * M_PI / lambda0;
        complex double wave[3] = {0};
        complex double curl_wave[3];
        if (mp->PermittivityModel != RADIAL_PML) {
          incident_wave(x, y, z, k0, wave, curl_wave);
        }
        avg_sum[i][global_node] += cimag(wave[0]);
      } break;
      case AVG_EMINCI_Y: {
        const double c0 = 1.0 / (sqrt(upd->Free_Space_Permittivity * upd->Free_Space_Permeability));
        dbl x = fv->x[0];
        dbl y = fv->x[1];
        dbl z = fv->x[2];
        dbl freq = upd->EM_Frequency;
        dbl lambda0 = c0 / freq;
        dbl k0 = 2 * M_PI / lambda0;
        complex double wave[3] = {0};
        complex double curl_wave[3];
        if (mp->PermittivityModel != RADIAL_PML) {
          incident_wave(x, y, z, k0, wave, curl_wave);
        }
        avg_sum[i][global_node] += creal(wave[1]);
      } break;
      case AVG_EMINCI_Z: {
        const double c0 = 1.0 / (sqrt(upd->Free_Space_Permittivity * upd->Free_Space_Permeability));
        dbl x = fv->x[0];
        dbl y = fv->x[1];
        dbl z = fv->x[2];
        dbl freq = upd->EM_Frequency;
        dbl lambda0 = c0 / freq;
        dbl k0 = 2 * M_PI / lambda0;
        complex double wave[3] = {0};
        complex double curl_wave[3];
        if (mp->PermittivityModel != RADIAL_PML) {
          incident_wave(x, y, z, k0, wave, curl_wave);
        }
        avg_sum[i][global_node] += cimag(wave[2]);
      } break;
      default:
        GOMA_EH(GOMA_ERROR, "Unknown nodal average non-variable type");
        break;
      }
    }
  }
}

/*****************************************************************************/

void post_process_nodal(double x[],            /* Solution vector for the current processor */
                        double **x_sens_p,     /* solution sensitivities */
                        double x_old[],        /* Solution vector at last time step */
                        double xdot[],         /* time derivative of solution vector */
                        double xdot_old[],     /* dx/dt at previous time step */
                        double resid_vector[], /* Residual vector for the
                                                  current proc */
                        int ts,
                        double *time_ptr,
                        double delta_t,
                        double theta,
                        double *x_pp,
                        Exo_DB *exo,
                        Dpi *dpi,
                        RESULTS_DESCRIPTION_STRUCT *rd,
                        char filename[],
                        int matrix_offset)

/*******************************************************************************
  Function which directs calculation of stream function & other post processing
  quantities

  Author:          P. R. Schunk
  Date:            12 September 1994
  Revised:         by many

*******************************************************************************/

{
  /* TAB certifies that this function conforms to the exo/patran side numbering
   *  convention 11/10/98. */

  int idex;             /* temporary */
  int candidate;        /* prospective new unique node name in a SS
                         * node list */
  int length_node_list; /* count up length of uniquely specified node
                         * names in a SS node list */
  int eb_index;
  int mesh_exoid, cpu_word_size, io_word_size;
  float version;
  int err;             /* temp variable to hold diagnostic flags. */
  int ip;              /* ip is the local quadrature point index      */
  int ip_total;        /* ip_total is the total number of volume
                          quadrature points in the element            */
  int j, k, l;         /* local index loop counter                  */
  int i, id;           /* Index for the local node number - row     */
  int I;               /* Indeces for global node number - row      */
  int ielem = 0, iel;  /* Index for elements                        */
  int ielem_type;      /* Element type of the current element         */
  int ielem_dim;       /* Element physical dimension                  */
  int num_local_nodes; /* Number of local basis functions in the
                          current element                             */
  int iconnect_ptr;    /* Pointer to the beginning of the connectivity
                          list for the current element                */
  int eqn = -1;        /* equation name (R_ENERGY, etc) */
  int var;             /* variable name (TEMPERATURE, etc) */
  int w;               /* species counter */

  double s, t, u;         /* Gaussian-quadrature point locations        */
  int *bloated_node_list; /* temporary used to point to raw SS node list
                           * that replicates nodes needlessly */
  int num_universe_nodes; /* alias for dpi->num_universe_nodes */

  /*
   * LEGEND:
   *
   * output:
   *    dphixdm[i][j][k][n] - deriv of (deriv of phi[i] wrt to global coord[j])
   *                          wrt mesh disp. k at node n
   *
   * T                          temperature (just like T_quad)
   *                            at the jth node (dof in elem)
   * dTdx[i]                    the i-th component of temperature gradient
   */

  /* __________________________________________________________________________*/
  /* __________________________________________________________________________*/

  double xi[DIM];                 /* Local element coordinates of Gauss point. */
  double wt;                      /* Quadrature weights
                                     units - ergs/(sec*cm*K) = g*cm/(sec^3*K)    */
  dbl *pressure_elem_vect = NULL; /* vector to hold nodal pressure*/
  /* PRS Cludge for remeshing guys */

  /* new post processing array that contains the nodal values of all the
     post processing variables */
  double **post_proc_vect;

  double **lumped_mass; /* vector to hold lumped mass matrix         */
  int *listnd = NULL;
  int *listndm[MAX_CONC]; /*  */
  int *listnde = NULL;
  int check, e_start, e_end, mn;
  double del_stream_fcn[4];  /* elemental side increments to stream function
                                calculated by routine calc_stream_fcn       */
  double vel[MAX_PDIM][MDE]; /* array for local nodal velocity values */

  int ii, kount = 0;
  int kounte = 0, kountm[MAX_CONC];

  /* side-post stuff */
  double *local_post, *local_lumped;
  double phi_i, phi_j;

  int num_side_in_set; /* returned number of sides(faces or edges) is  */
  int num_elem_in_set, num_node_in_set, num_nodes_on_side;
  /* in the side set                              */
  int num_dist_fact_in_set; /* returned number of distribution factors in   */
  int *elem_list;
  int *node_list;
  int *ss_ids;
  int id_side, iss = 0, p, q, dim, ldof, jd, iapply, ss_index;
  int mode;
  int id_local_elem_coord[MAX_NODES_PER_SIDE];

  /*  particle tracking stuff  */

  double p_time, p_step, p_norm, p_normu[2], vscale = 0.0; /*   time stepping variables  */
  int p_done, p_converged, inewton, p_dim, p_N;            /*convergence criteria, etc. */
  double p_lambda, p_force[MAX_PDIM];                      /* particle time constant, ext. force */
  double p_mass[2 * MAX_PDIM][2 * MAX_PDIM + 1], pivot;    /*  matrix info for time-stepper*/
  double p_xold[MAX_PDIM], p_x[MAX_PDIM], sum;
  double p_vel[MAX_PDIM], p_velold[MAX_PDIM];
  double f_vel[MAX_PDIM], f_velold[MAX_PDIM];
  FILE *jfp = NULL;    /*  file pointer  */
  int velo_interp = 0; /*  velocity basis functions  */

  struct Porous_Media_Terms pm_terms; /*added for POROUS_LIQUID_ACCUM_RATE*/

  /*
   * BEGINNING OF EXECUTABLE STATEMENTS
   */

  num_universe_nodes = dpi->num_universe_nodes;

  /*****************************************************************************/
  /*                               BLOCK 0                                     */
  /*                        INITIALIZATION FOR ROUTINE                         */
  /*     INITIALIZATION THAT IS DEPENDENT ON THE CURRENT PROCESSOR NUMBER      */
  /*****************************************************************************/
  /*
   * For simple post processing tasks like obtaining L2 approximations
   * of discontinuous quantities and using mass matrix lumping, use
   * the total number of nodes visible to this processor, including
   * external nodes that are rightfully owned by other processors. The
   * basic idea is that any values computed for external nodes will be
   * suitably discarded during the recombination of the distributed
   * simulation results.
   */

  /*
   * Obtain an optimized element order map from the exoII database
   * if available. If it is not available, then calculate one
   */
  if (Num_Proc == 1) {
    listel = alloc_int_1(Num_Internal_Elems, 0);
    cpu_word_size = sizeof(dbl);
    io_word_size = 0;
    mesh_exoid = ex_open(ExoFile, EX_READ, &cpu_word_size, &io_word_size, &version);
    GOMA_EH(mesh_exoid, "ex_open");
    /*
     * Turn off annoying warning about nonexistent element order maps
     * Only describe the problem if it is fatal.
     */
    ex_opts(EX_ABORT);
    err = ex_get_map(mesh_exoid, listel);
    GOMA_EH(err, "ex_get_map");
    err = ex_close(mesh_exoid);
    GOMA_EH(err, "ex_close");

    /*
     * If a valid mapping was not storred in the original exodus file,
     * then let's calculate mapping using a our own poor man's algorithm,
     * that at least satifies the bare-bones criteria needed for
     * calculating a connected stream function field variable where a
     * connected mesh exists.
     */
    for (i = 0; i < Num_Internal_Elems; i++) {
      listel[i]--;
    }
    if ((err = check_elem_order(listel, exo)) != 0) {
      err = elem_order_for_nodal_connect(listel, exo);
    }
    /*
     * Convert the mapping to 1 to Num_Internal_Elems basis
     */
    for (i = 0; i < Num_Internal_Elems; i++) {
      listel[i]++;
    }
  }

  /*
   *   If we have no postprocessing vectors to calculate then
   *   perform a quick return
   */
  if (rd->TotalNVPostOutput == 0)
    return;

  /* Allocate memory for requested function vectors */

  post_proc_vect = (double **)smalloc(rd->TotalNVPostOutput * sizeof(double *));

  for (j = 0; j < rd->TotalNVPostOutput; j++) {
    post_proc_vect[j] = (double *)smalloc(num_universe_nodes * sizeof(double));
  }

  if (Num_Proc == 1) {
    if (STREAM != -1 && Num_Var_In_Type[pg->imtrx][R_MOMENTUM1]) {
      listnd = (int *)smalloc(num_universe_nodes * sizeof(int));
    }
  }

  check = 0;
  for (i = 0; i < upd->Num_Mat; i++) {
    if (pd_glob[i]->MeshMotion == LAGRANGIAN)
      check = 1;
    if (pd_glob[i]->MeshMotion == DYNAMIC_LAGRANGIAN)
      check = 1;
  }
  if (PRESSURE_CONT != -1 && (Num_Var_In_Type[pg->imtrx][R_MOMENTUM1] || check)) {
    /* PRS Cludge for remeshing guys */
    pressure_elem_vect = (double *)smalloc(Num_Internal_Elems * sizeof(double));
  }

  if (Num_Proc == 1) {
    if (FLUXLINES != -1 && Num_Var_In_Type[pg->imtrx][R_MASS]) {
      for (w = 0; w < upd->Max_Num_Species_Eqn; w++) {
        listndm[w] = (int *)smalloc(num_universe_nodes * sizeof(int));
        kountm[w] = 0;
      }

      if (DIFFUSION_VECTORS == -1)
        GOMA_EH(GOMA_ERROR, "Need diffusion vectors for flux-function");
      if (pd->TimeIntegration != STEADY) {
        DPRINTF(stdout, "#####################################################");
        DPRINTF(stdout, "# WARNING: Fluxlines invalid in transient solutions #");
        DPRINTF(stdout, "#####################################################");
      }
    }
  }

  if (Num_Proc == 1) {
    if (ENERGY_FLUXLINES != -1 && Num_Var_In_Type[pg->imtrx][R_ENERGY]) {
      listnde = (int *)smalloc(num_universe_nodes * sizeof(int));
      if (CONDUCTION_VECTORS == -1)
        GOMA_EH(GOMA_ERROR, "Need conduction vectors for flux-function");
      if (pd->TimeIntegration != STEADY) {
        DPRINTF(stdout, "######################################################");
        DPRINTF(stdout, "# WARNING: Fluxlines invalid in transient solutions #");
        DPRINTF(stdout, "######################################################");
      }
    }
  }

  /* If L2 fitting, allocate memory for lumped mass and rhs
   * wouldn't be here if it weren't */

  lumped_mass = (double **)smalloc(rd->TotalNVPostOutput * sizeof(double *));
  for (j = 0; j < rd->TotalNVPostOutput; j++) {
    lumped_mass[j] = (double *)smalloc(num_universe_nodes * sizeof(double));
  }

  /*Initialize */
  if (Num_Proc == 1) {
    if (STREAM != -1 && Num_Var_In_Type[pg->imtrx][R_MOMENTUM1]) {
      for (I = 0; I < num_universe_nodes; I++) {
        listnd[I] = 0;
      }
    }
    if (ENERGY_FLUXLINES != -1 && Num_Var_In_Type[pg->imtrx][R_ENERGY]) {
      for (I = 0; I < num_universe_nodes; I++) {
        listnde[I] = 0;
      }
    }
    if (FLUXLINES != -1 && Num_Var_In_Type[pg->imtrx][R_MASS]) {
      for (w = 0; w < upd->Max_Num_Species_Eqn; w++) {
        for (I = 0; I < num_universe_nodes; I++) {
          listndm[w][I] = 0;
        }
      }
    }
  }

  if (rd->TotalNVPostOutput > 0) {
    for (i = 0; i < rd->TotalNVPostOutput; i++) {
      for (I = 0; I < num_universe_nodes; I++) {
        lumped_mass[i][I] = 0.;
        post_proc_vect[i][I] = 0.;
      }
    }
  }

  /*  af->Assemble_Jacobian = FAE; */

  /******************************************************************************/
  /*                                BLOCK 1                                     */
  /*          LOOP OVER THE ELEMENTS DEFINED ON THE CURRENT PROCESSOR           */
  /*          INITIALIZATION THAT IS DEPENDENT ON THE ELEMENT NUMBER            */
  /*          First calculate field variables using volume integrals            */
  /*          Then calculate stream and flux functions (need fields of          */
  /*          diffusion and conduction vectors for the flux functions           */
  /******************************************************************************/

  /*
   * Loop over all the elements, calculating the required interaction
   * coefficients
   */

  /*
   * Eventually, you may want to have an outermost loop over element blocks
   * where the "pd" and material ID's could potentially change, then loop
   * over each element of each element block.
   */

  /*
   * Note more robust distributed processing version: loop is over the
   * local element block indeces. The appropriate material index is found
   * from Matilda[].
   */

  for (eb_index = 0; eb_index < exo->num_elem_blocks; eb_index++) {
    mn = Matilda[eb_index];

    pd = pd_glob[mn];
    cr = cr_glob[mn];
    elc = elc_glob[mn];
    elc_rs = elc_rs_glob[mn];
    gn = gn_glob[mn];
    mp = mp_glob[mn];
    vn = vn_glob[mn];
    evpl = evpl_glob[mn];

    for (mode = 0; mode < vn->modes; mode++) {
      ve[mode] = ve_glob[mn][mode];
    }

    e_start = exo->eb_ptr[eb_index];
    e_end = exo->eb_ptr[eb_index + 1];

    for (iel = e_start; iel < e_end; iel++) {
      ielem = iel;

      PRS_mat_ielem = ielem - e_start; /*added for hysteretic saturation func*/

      /*
       * For each variable there are generally different degrees of
       * freedom that they and their equations contribute to.
       *
       * For this element, load up arrays that tell, for each variable,
       *
       *        (i) how many degrees of freedom they contribute towords
       *        (ii) the local node number associated with each degree of
       *             freedom
       *        (iii) pointers in the "esp" structure that tell where
       *              things are located in the global scheme...
       *                (a) nodal unknowns in this element...
       *                (b) Residual equations receiving contributions in
       *                    this element.
       *                (c) where the Jacobian entries go in the global
       *                    "a" matrix in its MSR format...
       */

      err = load_elem_dofptr(ielem, exo, x, x_old, xdot, xdot_old, 0);
      GOMA_EH(err, "load_elem_dofptr");
      err = bf_mp_init(pd);
      GOMA_EH(err, "bf_mp_init");
      ielem_type = ei[pg->imtrx]->ielem_type;
      ip_total = elem_info(NQUAD, ielem_type); /* number of
                                                * quadrature pts */

      num_local_nodes = ei[pg->imtrx]->num_local_nodes; /* number of
                                                         * local basis
                                                         * functions */

      ielem_dim = ei[pg->imtrx]->ielem_dim; /* physical
                                             * dimension of this
                                             * element */

      iconnect_ptr = ei[pg->imtrx]->iconnect_ptr; /* pointer to
                                                   * beginning of this
                                                   * element's
                                                   * connectivity list */

      struct Petrov_Galerkin_Data pg_data;

      memset(pg_data.h, 0, sizeof(double) * DIM);
      memset(pg_data.hh, 0, sizeof(double) * DIM * DIM);
      memset(pg_data.dh_dxnode, 0, sizeof(double) * MDE * DIM);
      memset(pg_data.hsquared, 0, sizeof(double) * DIM);
      memset(pg_data.hhv, 0, sizeof(double) * DIM * DIM);
      memset(pg_data.dhv_dxnode, 0, sizeof(double) * MDE * DIM);
      memset(pg_data.v_avg, 0, sizeof(double) * DIM);
      memset(pg_data.dv_dnode, 0, sizeof(double) * MDE * DIM);
      pg_data.mu_avg = 0.;
      pg_data.rho_avg = 0.;

      /* get element level constants for upwinding and
         stabilized schemes, if necessary */
      /*
       * If PSPG is turned on, then calculate the centroid viscosity
       * for use in the PSPG formulas. Note, we actually call
       * load_basis_functions here. Is this big penalty necessary or
       * can be piggyback on top of one gauss point?
       */
      int pspg_local = 0;
      if (PSPG) {
        if (PSPG == 1) {
          pspg_local = 0;
        }
        /* This is the flag for the standard local PSPG */
        else if (PSPG == 2) {
          pspg_local = 1;
        }
      }

      if ((PSPG || (mp->Mwt_funcModel == SUPG)) && pd_glob[mn]->e[pg->imtrx][R_PRESSURE] &&
          pd_glob[mn]->e[pg->imtrx][R_MOMENTUM1]) {
        xi[0] = 0.0;
        xi[1] = 0.0;
        xi[2] = 0.0;
        (void)load_basis_functions(xi, bfd);
        setup_shop_at_point(ielem, xi, exo);
        pg_data.mu_avg = element_viscosity();
        pg_data.rho_avg = density(NULL, *time_ptr);

        if (pspg_local) {
          h_elem_siz(pg_data.hsquared, pg_data.hhv, pg_data.dhv_dxnode, pd_glob[mn]->gv[R_MESH1]);
          element_velocity(pg_data.v_avg, pg_data.dv_dnode, exo);
        }
      }
      /******************************************************************************/
      /*                              BLOCK 1A                                      */
      /*                   START OF VOLUME INTEGRATION LOOP                         */
      /*                LOOP OVER THE NUMBER OF QUADRATURE POINTS                   */
      /*    Here we are doing a L2 projection onto the grid points of selected      */
      /*    quantities                                                              */
      /******************************************************************************/
      /* Loop over all the Volume Quadrature integration points */

      for (ip = 0; ip < ip_total; ip++) {

        MMH_ip = ip; /*Added for hysteretic saturation func.*/

        find_stu(ip, ielem_type, &s, &t, &u);
        xi[0] = s;
        xi[1] = t;
        xi[2] = u;

        /*
         * find quadrature weights for current ip
         */
        wt = Gq_weight(ip, ielem_type);
        fv->wt = wt;

        /*
         * Load up basis function information for ea variable...
         * Old usage: fill_shape
         */
        err = load_basis_functions(xi, bfd);
        GOMA_EH(err, "problem from load_basis_functions");

        /*
         * This has elemental Jacobian transformation and some
         * basic mesh derivatives...
         * Old usage: calc_Jac, jelly_belly
         */
        err = beer_belly();
        GOMA_EH(err, "beer_belly");

        /*
         * Load up field variable values at this Gauss point.
         */
        err = load_fv();
        GOMA_EH(err, "load_fv");

        /*
         * Here, load in the final part of the necessary basis function
         * information derivatives in the physical space coordinates...
         *
         *			grad(phi_i)
         *
         *			grad(phi_i e_a)
         * where:
         *		phi_i is the basis function at i-th dof
         *		e_a   is a unit vector in the coordinate system
         *
         * 		grad() operator depends on the coordinate system.
         */
        err = load_bf_grad();
        GOMA_EH(err, "load_bf_grad");

        err = load_fv_vector();

        /*
         * Load up physical space gradients of field variables at this
         * Gauss point.
         */
        err = load_fv_grads();
        GOMA_EH(err, "load_fv_grads");

        /*
         * Load up porous media variables and properties, if needed
         */
        if (mp->PorousMediaType == POROUS_UNSATURATED || mp->PorousMediaType == POROUS_SATURATED ||
            mp->PorousMediaType == POROUS_TWO_PHASE) {
          err = load_porous_properties();
          GOMA_EH(err, "load_porous_properties");
        }
        if (mp->PorousMediaType == POROUS_SATURATED) {
          err = get_porous_fully_sat_terms(&pm_terms, *time_ptr, delta_t);
          GOMA_EH(err, "problem in getting the saturated porous darcy  terms");
        } else if (mp->PorousMediaType == POROUS_UNSATURATED ||
                   mp->PorousMediaType == POROUS_TWO_PHASE) {
          err = get_porous_part_sat_terms(&pm_terms, *time_ptr, delta_t);
          GOMA_EH(err, "problem in getting the partially-saturated porous  terms");
        }

        computeCommonMaterialProps_gp(*time_ptr);

        /*
         * Calculate the contribution from this element of the
         * projection of the standard field variables unto the
         * node variables
         */
        err = calc_standard_fields(post_proc_vect, lumped_mass, delta_t, theta, ielem, ielem_type,
                                   ip, ip_total, rd, &pm_terms, *time_ptr, exo, xi, &pg_data);
        GOMA_EH(err, "calc_standard_fields");
      } /* END  for (ip = 0; ip < ip_total; ip++)                      */
    }   /* END  for (iel = 0; iel < num_internal_elem; iel++)            */
  }     /* END for (ieb loop) */

  /* Solve linear system for requested fields and put back in rhs vector*/

  /******************************************************************************/
  /*                                BLOCK 2                                     */
  /*                            SURFACE VARIABLES                               */
  /*     Evaluate Post-Processing Variables that are only defined on surfaces   */
  /******************************************************************************/
  if (SURFACE_VECTORS != -1 && pd->e[pg->imtrx][R_MESH1] && Num_ROT == 0) {

    /* loop over side-sets and find surface vectors for all elements on each
     * side-set */

    /*
     * We already read in this information, just refer to it if possible.
     */

    ss_ids = exo->ss_id;

    local_post = (double *)smalloc(rd->TotalNVPostOutput * sizeof(double));
    local_lumped = (double *)smalloc(rd->TotalNVPostOutput * sizeof(double));

    for (iss = 0; iss < exo->num_side_sets; iss++) {
      num_side_in_set = exo->ss_num_sides[iss];
      num_dist_fact_in_set = exo->ss_num_distfacts[iss];

      num_elem_in_set = num_side_in_set;
      num_node_in_set = num_dist_fact_in_set;

      num_nodes_on_side = num_node_in_set / num_elem_in_set; /* Be it
                                                              * hereby
                                                              * recorded
                                                              * that this
                                                              * is done
                                                              * under
                                                              * protest. */

      bloated_node_list = &(exo->ss_node_list[iss][0]);
      elem_list = &(exo->ss_elem_list[exo->ss_elem_index[iss]]);

      /*
       * Create node_list from original bloated list by adding only unique
       * node numbers to it.
       */

      node_list = (int *)smalloc(num_node_in_set * sizeof(int));

      for (j = 0; j < num_node_in_set; j++) {
        node_list[j] = -1;
      }

      length_node_list = 0;

      for (j = 0; j < num_node_in_set; j++) {
        candidate = bloated_node_list[j];
        if (in_list(candidate, 0, length_node_list, node_list) == -1) {
          node_list[length_node_list] = candidate;
          length_node_list++;
        }
      }

      node_list = (int *)realloc(node_list, length_node_list * sizeof(int));

      num_node_in_set = length_node_list;

      /*******************************************************************/
      /*                      BLOCK 2A                                   */
      /*            SURFACE INTEGRATION or NODAL LOOP                    */
      /* Evaluate Post-Processing Variables that are only defined on     */
      /* surfaces.                                                       */
      /*******************************************************************/

      for (i = 0; i < num_elem_in_set; i++) {
        err = load_elem_dofptr(elem_list[i], exo, x, x_old, xdot, xdot_old, 0);
        err = bf_mp_init(pd);
        GOMA_EH(err, "load_elem_dofptr");
        iconnect_ptr = ei[pg->imtrx]->iconnect_ptr;
        ielem_type = ei[pg->imtrx]->ielem_type;
        ip_total = elem_info(NQUAD_SURF, ielem_type);
        num_local_nodes = ei[pg->imtrx]->num_local_nodes;
        ielem_dim = ei[pg->imtrx]->ielem_dim;
        dim = ielem_dim;

        id_side =
            find_id_side(ei[pg->imtrx]->ielem, num_nodes_on_side,
                         &exo->ss_node_list[iss][num_nodes_on_side * i], id_local_elem_coord, exo);
        /*
         * Here, we will either perform surface integral over element with
         * Gaussian Quadrature, or evaluate the surface properties
         * at each nodal point (I think in a fine enough mesh, the two
         * become the same thing) */

        /* use nodal points only!! */
        for (k = 0; k < num_nodes_on_side; k++) {
          /* find where to put this nodal value in ss-list */
          /* Find the local element node number for the current node */
          id = id_local_elem_coord[k];

          I = exo->node_list[ei[pg->imtrx]->iconnect_ptr + id];
          find_nodal_stu(id, ielem_type, &xi[0], &xi[1], &xi[2]);

          err = load_basis_functions(xi, bfd);
          GOMA_EH(err, "problem from load_basis_functions");

          err = beer_belly();
          GOMA_EH(err, "beer_belly");

          /* precalculate variables at  current integration pt.*/
          err = load_fv();
          GOMA_EH(err, "load_fv");

          err = load_bf_grad();
          GOMA_EH(err, "load_bf_grad");

          err = load_fv_vector();

          err = load_bf_mesh_derivs();
          GOMA_EH(err, "load_bf_mesh_derivs");

          surface_determinant_and_normal(ei[pg->imtrx]->ielem, iconnect_ptr, num_local_nodes,
                                         ielem_dim - 1, id_side, num_nodes_on_side,
                                         id_local_elem_coord);

          /* calculate the components of the surface normal and
           * mesh displacement derivatives
           */

          if (ielem_dim != 3 && ielem_dim == pd->Num_Dim) {
            calc_surf_tangent(ei[pg->imtrx]->ielem, iconnect_ptr, num_local_nodes, ielem_dim - 1,
                              num_nodes_on_side, id_local_elem_coord);
          }

          /*
           * Put local contributions into global right-hand side
           * if it is not a right-hand side variable - it won't get
           * added in (contribution is zero)
           */
          if (SURFACE_VECTORS != -1 && pd->e[pg->imtrx][R_MESH1] && dim == 2) {

            for (p = 0; p < rd->TotalNVPostOutput; p++) {
              local_post[p] = 0.;
              local_lumped[p] = 0.;
            }

            /* Set flag to indicate if we're in the right material
             * (only one) to apply
             */
            iapply = 0;
            if ((ss_index = in_list(ss_ids[iss], 0, Proc_Num_Side_Sets, ss_to_blks[0])) == -1) {
              GOMA_EH(GOMA_ERROR, "Cannot match side SSID to ss_to_blks[].");
            }

            if (exo->eb_id[find_elemblock_index(ei[pg->imtrx]->ielem, exo)] ==
                ss_to_blks[1][ss_index]) {
              iapply = 1;
            }
            if (iapply) {
              for (p = 0; p < dim; p++) {
                idex = SURFACE_VECTORS + p;
                local_post[idex] = fv->snormal[p];
                local_lumped[idex] = 1.;
                local_post[idex + dim] = fv->stangent[0][p];
                local_lumped[idex + dim] = 1.;
                /*   not for 2D problems       local_post[idex+2*dim]   = fv->stangent[1][p];
                                               local_lumped[idex+2*dim] = 1.;  */
              }
            }
          } /* end of Surface vectors */

          /*
           * Choose the same weighting function for all
           * post-processing variables  - really want this to be
           * the highest order weighting function
           */
          eqn = pd->ProjectionVar;

          /* also convert from node number to dof number */
          ldof = ei[pg->imtrx]->ln_to_dof[eqn][id];
          if (ldof < 0) {
            GOMA_EH(GOMA_ERROR, "post_process_nodal: bad surface projection");
          }
          phi_i = bf[eqn]->phi[ldof];
          for (var = 0; var < rd->TotalNVPostOutput; var++) {
            post_proc_vect[var][I] += (local_post[var] * phi_i * fv->sdet * fv->h3);
            for (j = 0; j < num_nodes_on_side; j++) {
              /* Find the local element node number for the
               * current node
               */
              jd = id_local_elem_coord[j];
              /* also convert from node number to dof number */
              phi_j = bf[eqn]->phi[ei[pg->imtrx]->ln_to_dof[eqn][jd]];
              lumped_mass[var][I] += (local_lumped[var] * fv->sdet * phi_i * phi_j * fv->h3);
            }
          }
        } /* end of node-point loop */
      }   /* end of element loop on side-set */
    }     /* end of SS loop */

    safer_free((void **)&local_post);
    safer_free((void **)&local_lumped);
  } /* end of if surface variables */

  /* if ROTATION conditions are defined, determine these tangents separately
   * and plug them directly into the post_proc_vect */
  if (SURFACE_VECTORS != -1 && Num_ROT > 0) {
    dim = pd_glob[0]->Num_Dim;
    /*
     * NORMAL, TANGENT and OTHER Vectors required for ROTATION are calculated
     * ahead of time so we don't run into anomolous behavior due to neclaced
     * elements, junction points, . . .
     */
    calculate_all_rotation_vectors(exo, x);

    for (I = 0; I < num_universe_nodes; I++) {
      if (rotation[I][VECT_EQ_MESH] != NULL) {
        for (p = 0; p < dim; p++) {
          if (rotation[I][VECT_EQ_MESH][p]->ok) {
            for (q = 0; q < dim; q++) {
              post_proc_vect[SURFACE_VECTORS + p * dim + q][I] =
                  rotation[I][VECT_EQ_MESH][p]->vector[q];
              lumped_mass[SURFACE_VECTORS + p * dim + q][I] = 0.;
            }
          }
        }
      }
    }
  }

  /******************************************************************************/
  /*                      SMOOTHING                                             */
  /******************************************************************************/

#ifdef GOMA_ENABLE_PETSC
#if !(PETSC_USE_COMPLEX)
  if (upd->petsc_solve_post_proc) {
    petsc_solve_post_proc(post_proc_vect, rd, dpi);
  } else {
#endif
#endif
    for (ii = 0; ii < rd->TotalNVPostOutput; ii++) {
      for (I = 0; I < num_universe_nodes; I++) {
        if (fabs(lumped_mass[ii][I]) > DBL_SMALL) {
          post_proc_vect[ii][I] /= lumped_mass[ii][I];
        } else {
          post_proc_vect[ii][I] = 0.0;
        }
      }
    }
#ifdef GOMA_ENABLE_PETSC
#if !(PETSC_USE_COMPLEX)
  }
#endif
#endif

  for (ii = 0; ii < rd->TotalNVPostOutput; ii++) {
    exchange_node(cx[0], dpi, post_proc_vect[ii]);
  }

  /******************************************************************************/
  /*                                BLOCK 3                                     */
  /*                CALCULATE STREAM OR FLUX FUNCTIONS                          */
  /******************************************************************************/

  if (Num_Proc == 1) {
    if (STREAM != -1 || FLUXLINES != -1 || ENERGY_FLUXLINES != -1) {
      /*
       * Loop over all the elements, in the "optimal" order.
       */

      e_start = exo->eb_ptr[0];
      e_end = exo->eb_ptr[exo->num_elem_blocks];
      for (iel = e_start; iel < e_end; iel++) {
        ielem = listel[iel] - 1;

        /*
         * For each variable there are generally different degrees of
         * freedom that they and their equations contribute to.
         *
         * For this element, load up arrays that tell, for each variable,
         *
         *        (i) how many degrees of freedom they contribute towords
         *        (ii) the local node number associated with each degree
         *             of freedom
         *        (iii) pointers in the "esp" structure that tell where
         *              things are located in the global scheme...
         *                (a) nodal unknowns in this element...
         *                (b) Residual equations receiving contributions
         *                    in this element.
         *                (c) where the Jacobian entries go in the global
         *                    "a" matrix in its MSR format...
         */

        err = load_elem_dofptr(ielem, exo, x, x_old, xdot, xdot_old, 0);

        GOMA_EH(err, "load_elem_dofptr");

        iconnect_ptr = ei[pg->imtrx]->iconnect_ptr;

        ielem_type = ei[pg->imtrx]->ielem_type;
        ip_total = elem_info(NQUAD_SURF, ielem_type);
        num_local_nodes = ei[pg->imtrx]->num_local_nodes;

        ielem_dim = ei[pg->imtrx]->ielem_dim;

        dim = ielem_dim;

        /* Initialize velocity */
        for (i = 0; i < MAX_PDIM; i++) {
          for (j = 0; j < MDE; j++) {
            vel[i][j] = 0.;
          }
        }

        /* get convection velocity at the nodal points */

        if (STREAM != -1 && Num_Var_In_Type[pg->imtrx][VELOCITY1] &&
            ei[pg->imtrx]->ielem_dim == 2) {
          /* Go for it -- Calculate the stream function at the nodes
           * of this element
           */

          /*	if (ei[pg->imtrx]->ielem_dim == 2 && ei[pg->imtrx]->num_local_nodes == 9) { */

          if (ei[pg->imtrx]->ielem_dim == 2) {
            if (pd->e[pg->imtrx][R_MOMENTUM1]) {
              for (i = 0; i < ei[pg->imtrx]->ielem_dim; i++) {
                var = VELOCITY1 + i;
                for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                  vel[i][j] = *esp->v[i][j];
                }
              }
            }
            err = calc_stream_fcn(x, del_stream_fcn, vel);
          } else {
#if 1
            GOMA_EH(GOMA_ERROR, "No stream function in 3D ");
#endif
          }
          err = correct_stream_fcn(&kount, iel, del_stream_fcn, post_proc_vect[STREAM], listnd);

        } /* END of if(STREAM) */

        check = 0;
        for (i = 0; i < upd->Num_Mat; i++) {
          if (pd_glob[i]->MeshMotion == ARBITRARY)
            check = 1;
          if ((pd_glob[i]->MeshMotion == LAGRANGIAN ||
               pd_glob[i]->MeshMotion == DYNAMIC_LAGRANGIAN) &&
              pd_glob[i]->MeshInertia)
            check = 2;
        }

        if (FLUXLINES != -1 && Num_Var_In_Type[pg->imtrx][R_MASS]) {
          for (w = 0; w < upd->Max_Num_Species_Eqn; w++) {
            /* Go for it -- Calculate the flux function at the nodes
             *  of this element
             */

            if (ei[pg->imtrx]->ielem_dim == 2) {
              for (i = 0; i < ei[pg->imtrx]->ielem_dim; i++) {
                var = MASS_FRACTION;
                for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                  I = Proc_Elem_Connect[Proc_Connect_Ptr[iel] + j];
                  /* add in contribution from diffusion */
                  vel[i][j] =
                      post_proc_vect[DIFFUSION_VECTORS + w * ei[pg->imtrx]->ielem_dim + i][I];
                  if ((check == 1) && Num_Var_In_Type[pg->imtrx][VELOCITY1]) {
                    vel[i][j] += ((*esp->v[i][j]) * (*esp->c[w][j]));
                  } else if (check == 2)
                  /* use convection velocity for lagrangian
                     mesh motion */
                  {
                    I = Proc_Elem_Connect[Proc_Connect_Ptr[iel] + j];
                    if (TimeIntegration != STEADY) {
                      printf("\n Need to update Lagrangian convection velocities for transient \n");
                      vel[i][j] += post_proc_vect[LAGRANGE_CONVECTION + i][I] * (*esp->c[w][j]);
                    } else {
                      vel[i][j] += post_proc_vect[LAGRANGE_CONVECTION + i][I] * (*esp->c[w][j]);
                    }
                  }
                }
              }
              err = calc_stream_fcn(x, del_stream_fcn, vel);
            } else {
              GOMA_EH(GOMA_ERROR, "No flux lines in 3D ");
            }
            err = correct_stream_fcn(&kountm[w], iel, del_stream_fcn, post_proc_vect[FLUXLINES + w],
                                     listndm[w]);

          } /* END of loop over components */
        }   /* END of if(FLUXLINES) */

        if (ENERGY_FLUXLINES != -1 && Num_Var_In_Type[pg->imtrx][R_ENERGY]) {
          /* Go for it -- Calculate the stream function at the nodes of this element */
          if (ei[pg->imtrx]->ielem_dim == 2) {
            for (i = 0; i < ei[pg->imtrx]->ielem_dim; i++) {
              var = TEMPERATURE;
              for (j = 0; j < ei[pg->imtrx]->dof[var]; j++) {
                I = Proc_Elem_Connect[Proc_Connect_Ptr[iel] + j];
                vel[i][j] = post_proc_vect[CONDUCTION_VECTORS + i][I];
                if ((check == 1) && Num_Var_In_Type[pg->imtrx][VELOCITY1]) {
                  vel[i][j] += *esp->v[i][j] * (*esp->T[j]);
                } else if (check == 2)
                /* use convection velocity for lagrangian mesh motion */
                {
                  I = Proc_Elem_Connect[Proc_Connect_Ptr[iel] + j];
                  if (TimeIntegration != STEADY) {
                    printf("\n Need to update Lagrangian convection velocities for transient \n");
                    vel[i][j] += post_proc_vect[LAGRANGE_CONVECTION + i][I] * (*esp->T[j]);
                  } else {
                    vel[i][j] += post_proc_vect[LAGRANGE_CONVECTION + i][I] * (*esp->T[j]);
                  }
                }
              }
            }
            err = calc_stream_fcn(x, del_stream_fcn, vel);
          } else {
            GOMA_EH(GOMA_ERROR, "No energy flux lines in 3D ");
          }
          err = correct_stream_fcn(&kounte, iel, del_stream_fcn, post_proc_vect[ENERGY_FLUXLINES],
                                   listnde);
        } /* END of if(ENERGY_FLUXLINES) */

      } /* END of loop over elements */

    } /* END of if(STREAM or FLUXLINES) */

    /*
     * do some last minute processing for special variable types
     */

    if (STREAM != -1 && Num_Var_In_Type[pg->imtrx][VELOCITY1]) {
      /* First process global vector to recover average value at nodes*/
      for (I = 0; I < num_universe_nodes; I++) {
        ii = listnd[I];
        if (ii != 0) {
          post_proc_vect[STREAM][I] /= (double)ii;
        }
      }

      /* Compute value at midside nodes if these exist */
      if (ei[pg->imtrx]->num_local_nodes == 9) {
        err = midsid(post_proc_vect[STREAM], exo);
      }

      free(listnd);
    }

    if (FLUXLINES != -1 && Num_Var_In_Type[pg->imtrx][R_MASS]) {
      for (w = 0; w < upd->Max_Num_Species_Eqn; w++) {
        /* First process global vector to recover average value at nodes*/
        for (I = 0; I < num_universe_nodes; I++) {
          ii = listndm[w][I];
          if (ii != 0) {
            post_proc_vect[FLUXLINES + w][I] /= (double)ii;
          }
        }

        /* Then compute value at midside nodes */

        if (ei[pg->imtrx]->num_local_nodes == 9) {
          err = midsid(post_proc_vect[FLUXLINES + w], exo);
        }

        free(listndm[w]);
      }
    }

    if (ENERGY_FLUXLINES != -1 && Num_Var_In_Type[pg->imtrx][R_ENERGY]) {
      /* First process global vector to recover average value at nodes*/
      for (I = 0; I < num_universe_nodes; I++) {
        ii = listnde[I];
        if (ii != 0) {
          post_proc_vect[ENERGY_FLUXLINES][I] /= (double)ii;
        }
      }

      /* Then compute value at midside nodes */

      if (ei[pg->imtrx]->num_local_nodes == 9) {
        err = midsid(post_proc_vect[ENERGY_FLUXLINES], exo);
      }

      free(listnde);
    }

  } /* end of serial processing block for streamlines */

  if (NS_RESIDUALS != -1 && Num_Var_In_Type[pg->imtrx][R_MOMENTUM1]) {
    for (I = 0; I < num_universe_nodes; I++) {
      for (j = 0; j < ei[pg->imtrx]->ielem_dim; j++) {
        ii = Index_Solution(I, R_MOMENTUM1 + j, 0, 0, -2, pg->imtrx);
        if (ii != -1) {
          post_proc_vect[NS_RESIDUALS + j][I] = resid_vector[ii];
        }
      }
    }
  }

  if (NS_RESIDUALS != -1 && Num_Var_In_Type[pg->imtrx][R_MOMENTUM1]) {
    for (I = 0; I < num_universe_nodes; I++) {
      for (j = 0; j < ei[pg->imtrx]->ielem_dim; j++) {
        ii = Index_Solution(I, R_MOMENTUM1 + j, 0, 0, -2, pg->imtrx);
        if (ii != -1) {
          post_proc_vect[NS_RESIDUALS + j][I] = resid_vector[ii];
        }
      }
    }
  }

  if (MM_RESIDUALS != -1 && Num_Var_In_Type[pg->imtrx][R_MESH1]) {
    for (I = 0; I < num_universe_nodes; I++) {
      for (j = 0; j < ei[pg->imtrx]->ielem_dim; j++) {
        ii = Index_Solution(I, R_MESH1 + j, 0, 0, -2, pg->imtrx);
        if (ii != -1) {
          post_proc_vect[MM_RESIDUALS + j][I] = resid_vector[ii];
        }
      }
    }
  }

  if (efv->ev) {
    for (j = 0; j < efv->Num_external_field; j++) {
      if (efv->i[j] != I_TABLE) {
        for (I = 0; I < num_universe_nodes; I++) {
          post_proc_vect[EXTERNAL_POST + j][I] = efv->ext_fld_ndl_val[j][I];
        }
      }
    }
  }

  /* sum up all stress modes to compute the total stress */
  if (TOTAL_STRESS11 != -1 && Num_Var_In_Type[pg->imtrx][R_STRESS11]) {
    sum_total_stress(x, R_STRESS11, 0, post_proc_vect[TOTAL_STRESS11], exo);
  }
  if (TOTAL_STRESS12 != -1 && Num_Var_In_Type[pg->imtrx][R_STRESS12]) {
    sum_total_stress(x, R_STRESS12, 0, post_proc_vect[TOTAL_STRESS12], exo);
  }

  if (TOTAL_STRESS13 != -1 && Num_Var_In_Type[pg->imtrx][R_STRESS13]) {
    sum_total_stress(x, R_STRESS13, 0, post_proc_vect[TOTAL_STRESS13], exo);
  }

  if (TOTAL_STRESS22 != -1 && Num_Var_In_Type[pg->imtrx][R_STRESS22]) {
    sum_total_stress(x, R_STRESS22, 0, post_proc_vect[TOTAL_STRESS22], exo);
  }

  if (TOTAL_STRESS23 != -1 && Num_Var_In_Type[pg->imtrx][R_STRESS23]) {
    sum_total_stress(x, R_STRESS23, 0, post_proc_vect[TOTAL_STRESS23], exo);
  }

  if (TOTAL_STRESS33 != -1 && Num_Var_In_Type[pg->imtrx][R_STRESS33]) {
    sum_total_stress(x, R_STRESS33, 0, post_proc_vect[TOTAL_STRESS33], exo);
  }

  if (nn_average > 0) {
    post_process_average(x, x_old, xdot, xdot_old, resid_vector, exo, dpi, post_proc_vect,
                         *time_ptr);
    for (int ii = 0; ii < nn_average; ii++) {
      exchange_node(cx[0], dpi, post_proc_vect[pp_average[ii]->index_post]);
    }
  }

  /*****************************************************************************/
  /*                               BLOCK 4                                     */
  /*           NOW, write results onto exodusII database                       */
  /*****************************************************************************/
  /* ----------
   * NOW, write results onto exodusII database
   * ----------*/

  for (i = 0; i < rd->TotalNVPostOutput; i++) {

    /*
     * When we're parallel, this will write out results for external
     * nodes that are owned by other processors. When we "fix" the results
     * into a monolith we may simply ignore this processor's concept of
     * the nodal result variable at external nodes.
     */
    if (filename != NULL) {
      wr_nodal_result_exo(exo, filename, post_proc_vect[i],
                          matrix_offset + rd->TotalNVSolnOutput + i + 1, ts, *time_ptr);
    }
  }

  /*
   * Compute particle traces and output to file
   */

  if (nn_particles > 0) {
    double orig_coord[3];
    /* 1) Search Element Database and find originating element
      id for each particle */

    fprintf(stderr, "  starting particle trace computation...  \n");

    p_dim = exo->num_dim;

    for (i = 0; i < nn_particles; i++) {

      for (j = 0; j < p_dim; j++) {
        orig_coord[j] = pp_particles[i]->coord[j];
      }
      if (p_dim == 2) {
        pp_particles[i]->coord[2] = 0.;
        p_x[2] = 0.;
      }
      p_lambda = pp_particles[i]->mass * pp_particles[i]->mobility;
      p_force[0] = pp_particles[i]->mobility * pp_particles[i]->force[0];
      p_force[1] = pp_particles[i]->mobility * pp_particles[i]->force[1];
      p_force[2] = pp_particles[i]->mobility * pp_particles[i]->force[2];

      /*	if the product of the particle mass and mobility is zero or negative
            then the particle will follow the fluid pathlines.  If the product
            is positive, solve a double size set of equations with particle
            inertia and possibly external forces included.		*/

      if (p_lambda > 0.0) {
        p_N = 2 * p_dim;
      } else {
        p_N = p_dim;
      }

      if (pd->TimeIntegration == STEADY || tran->time_value == tran->init_time + tran->delta_t) {
        /*   get initial position    */

        pp_particles[i]->Current_element_id =
            find_id_elem(pp_particles[i]->coord[0], pp_particles[i]->coord[1],
                         pp_particles[i]->coord[2], x, exo, 0, exo->num_elems);

        if (pp_particles[i]->Current_element_id == -1)
          GOMA_EH(GOMA_ERROR, "Cannot locate one of your particles");
      } /* if steady	*/

      /* find basis functions associated with velocity variables */

      for (j = 0; j < Num_Basis_Functions; j++) {
        if (pd_glob[ei[pg->imtrx]->mn]->i[pg->imtrx][VELOCITY1] == bfd[j]->interpolation) {
          velo_interp = j;
        }
      }

      /*    get initial isoparametric coordinates - initial isopar. coords   */

      for (j = 0; j < MAX_PDIM; j++) {
        pp_particles[i]->xi_coord[j] = 0.;
      }

      err =
          invert_isoparametric_map(&pp_particles[i]->Current_element_id, &pp_particles[i]->coord[0],
                                   &pp_particles[i]->xi_coord[0], exo, x, &velo_interp);

      if (err == 0) {
        GOMA_WH(-1, "map inversion for particle path failed");
        goto next_particle;
      }

      for (j = 0; j < p_dim; j++) {
        if (fabs(pp_particles[i]->xi_coord[j]) > 2.0) {
          fprintf(stderr, "particle %d - initial pt outside the domain\n", i);
          fprintf(stderr, "isoparametric coords  %g %g %g \n", pp_particles[i]->xi_coord[0],
                  pp_particles[i]->xi_coord[1], pp_particles[i]->xi_coord[2]);
          goto next_particle;
        }
      }

      err = load_basis_functions(&pp_particles[i]->xi_coord[0], bfd);
      GOMA_EH(err, "problem from load_basis_functions");

      err = beer_belly();
      GOMA_EH(err, "beer_belly");

      err = load_fv();
      GOMA_EH(err, "load_fv");

      err = load_bf_grad();
      GOMA_EH(err, "load_bf_grad");

      err = load_fv_vector();

      err = load_fv_grads();
      GOMA_EH(err, "load_fv_grads");

      if (pd->TimeIntegration == STEADY || tran->time_value == tran->init_time + tran->delta_t) {
        vscale = 0.;
        for (j = 0; j < p_dim; j++) {
          p_x[j] = p_xold[j] = fv->x[j];
          p_vel[j] = p_velold[j] = fv->v[j];
          f_vel[j] = f_velold[j] = fv->v[j];
          pp_particles[i]->coord[j] = p_x[j];
          vscale += fv->v[j] * fv->v[j];
        }
        vscale = sqrt(vscale);
      } else {
        for (j = 0; j < p_dim; j++) {
          p_xold[j] = pp_particles[i]->coord[j];
          p_velold[j] = pp_particles[i]->p_velo[j];
          f_velold[j] = fv_old->v[j];
        }
      }

      /* 3) Start Time integration loop */

      if (pd->TimeIntegration == STEADY) {
        p_time = pp_particles[i]->Start_Time;
        p_step = pp_particles[i]->Delta_s / vscale;
      } else if (tran->time_value == tran->init_time + tran->delta_t) {
        p_time = tran->init_time;
        p_step = tran->delta_t;
      } else {
        p_time = tran->time_value - tran->delta_t;
        p_step = tran->delta_t;
      }

      p_done = FALSE;

      if (pd->TimeIntegration == STEADY || tran->time_value == tran->init_time + tran->delta_t) {
        /**  write file headings  **/

        if ((jfp = fopen(pp_particles[i]->filenm, "a")) != NULL) {
          err = usr_ptracking(jfp, i + 1, p_x, p_vel, p_xold, p_velold, TRUE, p_time, p_step);
          GOMA_EH(err, "problem with usr_ptracking");
        }
      } /* if steady	*/
      else {
        jfp = fopen(pp_particles[i]->filenm, "a");
      }

      while (!p_done) {

        /*  Adjust particle size for droplets  */
        if (pd->v[pg->imtrx][RESTIME] && mp->SpeciesSourceModel[0] == DROP_EVAP) {
          p_lambda = pp_particles[i]->mass * pp_particles[i]->mobility * SQUARE(fv->restime);
          p_force[0] = pp_particles[i]->mobility * pp_particles[i]->force[0] * SQUARE(fv->restime);
          p_force[1] = pp_particles[i]->mobility * pp_particles[i]->force[1] * SQUARE(fv->restime);
          p_force[2] = pp_particles[i]->mobility * pp_particles[i]->force[2] * SQUARE(fv->restime);
        }
        /*   make Euler predictor step  */

        for (j = 0; j < p_dim; j++) {
          p_x[j] = p_xold[j] + p_velold[j] * p_step;
        }
        if (p_lambda > 0.) {
          for (j = 0; j < p_dim; j++) {
            p_vel[j] = p_velold[j] + (f_velold[j] - p_velold[j] + p_force[j]) * p_step / p_lambda;
          }
        }

        p_time += p_step;

        /*  trapezoidal rule corrector step  */

        p_converged = FALSE;
        inewton = 0;

        while (!p_converged && inewton < 6) {
          err = invert_isoparametric_map(&pp_particles[i]->Current_element_id, p_x,
                                         &pp_particles[i]->xi_coord[0], exo, x, &velo_interp);

          if (err == 0) {
            GOMA_WH(-1, "map inversion for particle path failed");
            goto next_particle;
          }

          for (j = 0; j < p_dim; j++) {
            if (fabs(pp_particles[i]->xi_coord[j]) > 2.0) {
              fprintf(stderr, "particle %d has exited the domain\n", i);
              fprintf(stderr, "element = %d\n", pp_particles[i]->Current_element_id);
              fprintf(stderr, "coords  %g %g %g \n", p_x[0], p_x[1], p_x[2]);
              p_done = TRUE;
              goto next_particle;
            }
          }

          err = load_basis_functions(&pp_particles[i]->xi_coord[0], bfd);
          GOMA_EH(err, "problem from load_basis_functions");

          err = beer_belly();
          GOMA_EH(err, "beer_belly");

          err = load_fv();
          GOMA_EH(err, "load_fv");

          err = load_bf_grad();
          GOMA_EH(err, "load_bf_grad");

          err = load_fv_vector();

          err = load_fv_grads();
          GOMA_EH(err, "load_fv_grads");

          /*  assemble mass matrix and rhs */

          if (p_lambda > 0.0) {
            for (j = 0; j < p_dim; j++) {
              f_vel[j] = fv->v[j];
              p_mass[2 * j][p_N] = p_x[j] - p_xold[j] - 0.5 * p_step * (p_vel[j] + p_velold[j]);
              p_mass[2 * j + 1][p_N] =
                  p_lambda * (p_vel[j] - p_velold[j]) -
                  0.5 * p_step *
                      (f_vel[j] - p_vel[j] + f_velold[j] - p_velold[j] + 2. * p_force[j]);
              for (k = 0; k < p_dim; k++) {
                p_mass[2 * j][2 * k] = -delta(j, k);
                p_mass[2 * j][2 * k + 1] = 0.5 * p_step * delta(j, k);
                p_mass[2 * j + 1][2 * k] = 0.5 * p_step * fv->grad_v[k][j];
                p_mass[2 * j + 1][2 * k + 1] = -(0.5 * p_step + p_lambda) * delta(j, k);
              }
            }
          } else {
            for (j = 0; j < p_dim; j++) {
              f_vel[j] = fv->v[j];
              p_mass[j][p_dim] = p_x[j] - p_xold[j] - 0.5 * p_step * (f_vel[j] + f_velold[j]);
              for (k = 0; k < p_dim; k++) {
                p_mass[j][k] = 0.5 * p_step * fv->grad_v[k][j] - delta(j, k);
              }
            }
          }

          p_norm = 0;
          for (j = 0; j < p_N; j++) {
            p_norm += p_mass[j][p_N] * p_mass[j][p_N];
          }

          /*   solve linear system using straight Gauss elimination assuming
                          no pivoting required    */

          for (j = 0; j < p_N - 1; j++) {
            for (k = j + 1; k < p_N; k++) {
              pivot = p_mass[k][j] / p_mass[j][j];
              for (l = j + 1; l < p_N + 1; l++) {
                p_mass[k][l] -= p_mass[j][l] * pivot;
              }
            }
          }

          /*  backsubstitution  */

          p_mass[p_N - 1][p_N] = p_mass[p_N - 1][p_N] / p_mass[p_N - 1][p_N - 1];

          for (j = p_N - 2; j >= 0; j--) {
            sum = 0;
            for (k = j + 1; k < p_N; k++) {
              sum += p_mass[j][k] * p_mass[k][p_N];
            }
            p_mass[j][p_N] = (p_mass[j][p_N] - sum) / p_mass[j][j];
          }

          /**   update solution  **/
          p_normu[0] = 0;
          p_normu[1] = 0;
          if (p_lambda > 0.0) {
            for (j = 0; j < p_dim; j++) {
              p_x[j] += p_mass[2 * j][p_N];
              p_vel[j] += p_mass[2 * j + 1][p_N];
            }
            for (j = 0; j < p_dim; j++) {
              p_normu[0] += p_mass[2 * j][p_N] * p_mass[2 * j][p_N];
              p_normu[1] += p_mass[2 * j + 1][p_N] * p_mass[2 * j + 1][p_N];
            }
          } else {
            for (j = 0; j < p_dim; j++) {
              p_x[j] += p_mass[j][p_N];
              p_vel[j] = f_vel[j];
            }
            for (j = 0; j < p_dim; j++) {
              p_normu[0] += p_mass[j][p_N] * p_mass[j][p_N];
            }
          }

          p_converged = ((p_norm + p_normu[0] + p_normu[1]) < 1.0E-12);
          inewton++;
        } /*  corrector while loop  */

        /*   check convergence of corrector step  */

        if (!p_converged) {
          p_time -= p_step;
          p_step *= 0.5;
          if ((pd->TimeIntegration == STEADY &&
               p_step < 0.0001 * pp_particles[i]->Delta_s / vscale) ||
              (pd->TimeIntegration == TRANSIENT && p_step < tran->Delta_t_min)) {
            GOMA_WH(-1, "particle timestep below minimum");
            goto next_particle;
          }
          continue;
        } else {

          /* update variables */

          for (j = 0; j < p_dim; j++) {
            pp_particles[i]->coord[j] = p_x[j];
            pp_particles[i]->p_velo[j] = p_vel[j];
          }

          /*   call usr_ptracking for computation and output of other quantities
           *   along the particle traces
           */

          err = usr_ptracking(jfp, i + 1, p_x, p_vel, p_xold, p_velold, FALSE, p_time, p_step);
          GOMA_EH(err, "problem with usr_ptracking");

          vscale = 0.;
          for (j = 0; j < p_dim; j++) {
            p_xold[j] = p_x[j];
            p_velold[j] = p_vel[j];
            f_velold[j] = f_vel[j];
            vscale += p_vel[j] * p_vel[j];
          }
          vscale = sqrt(vscale);
          if (pd->TimeIntegration == STEADY) {
            p_step = MIN(pp_particles[i]->Delta_s / vscale, 1.2 * p_step);
          } else {
            p_step = MIN(tran->delta_t, 1.2 * p_step);
          }
        }
        p_done = ((pd->TimeIntegration == STEADY && p_time >= pp_particles[i]->End_Time) ||
                  (pd->TimeIntegration == TRANSIENT && p_time >= tran->time_value));

      } /*   while loop */

    next_particle:

      fflush(jfp);
      fclose(jfp);
      for (j = 0; j < p_dim; j++) {
        pp_particles[i]->coord[j] = orig_coord[j];
      }
    } /*  particle counter */
    fprintf(stderr, "  Done tracing %d particles...  \n", nn_particles);
  }

  /*  err = usr_print(time_ptr, delta_t, x, post_proc_vect, -1);  */

  check = 0;
  for (i = 0; i < upd->Num_Mat; i++) {
    if (pd_glob[i]->MeshMotion == LAGRANGIAN || pd_glob[i]->MeshMotion == DYNAMIC_LAGRANGIAN)
      check = 1;
  }
  if (PRESSURE_CONT != -1 && (Num_Var_In_Type[pg->imtrx][VELOCITY1] || check)) {
    /*PRS cludge for remeshing guys */
    for (I = 0; I < Num_Internal_Elems; I++) {
      pressure_elem_vect[I] = (dbl)I;
    }

    safe_free(pressure_elem_vect);
  }

  /*
   * Save porous saturation vector to x_pp.
   * Later, this will be usable for arbitrary post-processing variables.
   */
  /*
  printf(" Porous sat. ID is %d\n", POROUS_SATURATION);
  printf(" Capillary pressure ID is %d\n", CAPILLARY_PRESSURE);
  printf(" Number of PP vars from rd is %d\n", rd->TotalNVPostOutput);
  printf(" Pointer check of x_pp: %p\n", x_pp);
    if (POROUS_SATURATION != -1 && x_pp != NULL)
      {
        for (I=0; I<num_universe_nodes; I++)
          {
            x_pp[I] = post_proc_vect[POROUS_SATURATION][I];
          }
      }
  */

  if (Num_Export_XP > 0 && x_pp != NULL) {
    idex = 0;
    for (j = 0; j < Num_Export_XP; j++) {
      for (I = 0; I < exo->num_nodes; I++) {
        x_pp[idex] = post_proc_vect[Export_XP_ID[j]][I];
        idex++;
      }
    }
  }

  /*
   * write out pressure along chosen side set
   */
  /*   err = usr_print(NULL,0.,x,post_proc_vect, PRESSURE_CONT);
     err = usr_print(NULL,0.,x,post_proc_vect, STRESS_TENSOR);  */

  for (j = 0; j < rd->TotalNVPostOutput; j++) {
    safe_free(post_proc_vect[j]);
  }

  safe_free(post_proc_vect);

  safe_free(listel);

  for (j = 0; j < rd->TotalNVPostOutput; j++) {
    safe_free(lumped_mass[j]);
  }

  safe_free(lumped_mass);

  return;
  /*****************************************************************************/
} /*   END OF post_process_nodal                                             */
/*****************************************************************************/

void post_process_elem(double x[],     /* soln vector */
                       double x_old[], /* soln vector at previous time step */
                       double xdot[],  /* time derivative of soln vector */
                       double xdot_old[],
                       double resid_vector[], /* Residual vector */
                       const int tev,
                       const int tev_post,
                       double ***gvec_elem, /* Triply indexed array containing
                                             * element variable values on return.
                                             * convention:
                                             * [elemblock_index][elemvar_index]
                                             *      [element_index(inblock)] */
                       const int ts,
                       const double *time_ptr,
                       const double delta_t,
                       Exo_DB *const exo,
                       Dpi *const dpi,
                       struct Results_Description *rd)

/************************************************************************
 * Function which directs calculation of element-based post processing
 *quantities
 *
 * Author:          Randy R. Lober
 * Date:            13 August 1998
 * Revised:
 *
 ************************************************************************/
{
  int i, eb_indx, ev_indx, elem, compute_elem_size;
  int mn, ip_total, ielem, ip, ev_indx_tmp;
  ELEM_BLK_STRUCT *eb_ptr;

  if (Num_Elem_Post_Proc_Var == 0)
    return;

  /* Initialize  - NOTE ONLY initialize the members that have been malloc'd */
  i = 0;
  if (tev_post > 0) {
    for (eb_indx = 0; eb_indx < exo->num_elem_blocks; eb_indx++) {
      for (ev_indx = tev; ev_indx < tev + tev_post; ev_indx++) {
        if (exo->elem_var_tab[i++] == 1) { /* Checks to see if this ev_indx
                                              exists for this eb_indx */
          for (elem = 0; elem < exo->eb_num_elems[eb_indx]; elem++) {
            gvec_elem[eb_indx][ev_indx][elem] = 0.;
          }
        }
      }
    }
  }

  /* Compute the elem post vars here */
  ev_indx = tev;

  /* Zienkiewicz-Zhu error indicator based on velocity */

  if (ERROR_ZZ_VEL != -1 && Num_Var_In_Type[pg->imtrx][R_MOMENTUM1]) {
    compute_elem_size = 0;
    if (ERROR_ZZ_VEL_ELSIZE != -1)
      compute_elem_size = 1;
    if (calc_zz_error_vel(x, x_old, xdot, xdot_old, resid_vector,
                          ev_indx,   /* Variable index for zz_error  */
                          gvec_elem, /* elem var vals[eb_indx][ev_indx][elem] */
                          exo, dpi, compute_elem_size) != 0) {
      GOMA_EH(GOMA_ERROR, " calc_zz_error_vel failure");
    } else {
      /* If we get to here, calc_zz_error_vel worked and we need to increment ev_indx */
      if (compute_elem_size) {
        /* if the elem size was also computed, 2 ev's were added, increment by 2 */
        ev_indx += 2;
      } else {
        /* If only the error was computed, increment ev_indx by 1 */
        ev_indx++;
      }
    }
  }

  if (SAT_CURVE_TYPE != -1 && CAP_PRESS_SWITCH != -1 && SAT_QP_SWITCH != -1) {
    ev_indx_tmp = ev_indx;
    if (Num_Var_In_Type[pg->imtrx][R_POR_LIQ_PRES] ||
        Num_Var_In_Type[pg->imtrx][R_SHELL_SAT_OPEN] ||
        Num_Var_In_Type[pg->imtrx][R_SHELL_SAT_OPEN_2]) {
      for (eb_indx = 0; eb_indx < exo->num_elem_blocks; eb_indx++) {
        ev_indx = ev_indx_tmp;
        mn = Matilda[eb_indx];
        mp = mp_glob[mn];
        eb_ptr = Element_Blocks + eb_indx;
        ip_total = elem_info(NQUAD, eb_ptr->Elem_Type);
        if ((mp->PorousMediaType == POROUS_UNSATURATED ||
             mp->PorousMediaType == POROUS_SHELL_UNSATURATED ||
             mp->PorousMediaType == POROUS_TWO_PHASE) &&
            mp->SaturationModel == TANH_HYST) {
          for (ip = 0; ip < ip_total; ip++) {
            for (ielem = 0; ielem < eb_ptr->Num_Elems_In_Block; ielem++) {
              gvec_elem[eb_indx][ev_indx][ielem] += eb_ptr->ElemStorage[ielem].sat_curve_type[ip];
            }
            ev_indx++;
          }
          for (ip = 0; ip < ip_total; ip++) {
            for (ielem = 0; ielem < eb_ptr->Num_Elems_In_Block; ielem++) {
              gvec_elem[eb_indx][ev_indx][ielem] += eb_ptr->ElemStorage[ielem].Sat_QP_tn[ip];
            }
            ev_indx++;
          }
          for (ip = 0; ip < ip_total; ip++) {
            for (ielem = 0; ielem < eb_ptr->Num_Elems_In_Block; ielem++) {
              gvec_elem[eb_indx][ev_indx][ielem] += eb_ptr->ElemStorage[ielem].p_cap_QP[ip];
            }
            ev_indx++;
          }
        }
      } /*End of loop over element block */
    }
  }

  for (i = 0; i < Num_Elem_Post_Proc_Var; i++) {
    /*
     * When we're parallel, this will write out results for boundary
     * elements that are shared by other processors.
     */
    wr_elem_result_exo(exo, ExoFileOut, gvec_elem, tev + i, ts, *time_ptr, rd);
  }

  /* Release the struct array memory if it exists */
  if (nn_error_metrics > 0)
    free(pp_error_data);
}
/*****************************************************************************/

void post_process_global(double *x, /* Solution vector for the current processor */
                         Exo_DB *exo,
                         Dpi *dpi,
                         double time) {
  int i;
  for (i = 0; i < nn_global; i++) {
    switch (pp_global[i]->type) {
    case PP_GLOBAL_LS_INTERFACE_PRINT: {
      if (ls != NULL) {
        print_ls_interface(x, exo, dpi, time, pp_global[i]->filenm, FALSE);
      }
    } break;
    case PP_GLOBAL_LS_INTERFACE_PRINT_ALL_TIMES: {
      if (ls != NULL) {
        print_ls_interface(x, exo, dpi, time, pp_global[i]->filenm, TRUE);
      }
    } break;
    default:
      GOMA_EH(GOMA_ERROR, "Unknown global post process type");
      break;
    }
  }
}

/********************************************************************************/
/********************************************************************************/
/********************************************************************************/

static int calc_zz_error_vel(double x[], /* Solution vector                       */
                             double x_old[],
                             double xdot[],
                             double xdot_old[],
                             double resid_vector[],
                             int ev_indx,         /* Variable index for zz_error               */
                             double ***gvec_elem, /* evar vals[eb_indx][ev_indx][elem]  */
                             Exo_DB *const exo,
                             Dpi *const dpi,
                             int compute_elem_size)

/*
  Function which calculates the Zienkiewicz-Zhu error indicator for each
  appropriate element in the model

  Author:          R. R. Lober (9113)
  Date:            9 September 1998
  Revised

*/

{
  int status, i_node, i_elem, i_start, i_end, i, j, k, kk, elem_id, max_gp, min_gp;
  int num_elems_in_patch, i_elem_type, i_elem_gp, i_elem_dim, mat_num;
  int max_dim, valid_count, max_terms, last_interp, i_eb_indx;
  int do_the_lu_decomp, *indx, max_velocity_norm_i_elem = 0, max_velocity_err_i_elem = 0;
  int *valid_elem_mask, num_local_proc_nodes, err, global_node_num;
  int remesh_status;
  double **xgp_loc, **ygp_loc, **zgp_loc, **det_gp_loc, **s_lhs, *rhs, ***tau_lsp;
  double ****tau_gp_ptch, *i_node_coords, **wt_gp_loc, tau, xgp, ygp, zgp, det, wt;
  double max_velocity_norm, max_velocity_norm_tmp, *elem_areas = NULL, h1, h2, error_ratio;
  double expansion_rate, reduction_rate, elem_size_min, elem_size_max, target_error;
  double max_velocity_err, max_x = 0, max_y = 0, max_z = 0, total_volume, pct_over_target,
                           volume_tmp;
  double pct_over_tolerance;

  status = 0;

  /* Sanity check for mesh topology data */
  if (!exo->node_elem_conn_exists) {
    GOMA_EH(GOMA_ERROR, "Attempt to access undeveloped node->elem connectivity.");
  }

  i_node_coords = (double *)smalloc(3 * sizeof(double));

  num_local_proc_nodes = dpi->num_internal_nodes + dpi->num_external_nodes;
  /* Per conversation with PAS 9/15, in the parallel world of goma, each node
     has a unique and unambiguous owning proc, while each proc has elements that
     are shared and termed boundary elements. To ensure our data requests are
     appropriately on proc, we will use node indices as follows:
     our local view of all the nodes will be from
        node index 0 to ( dpi->num_internal_nodes + dpi->num_external_nodes )
     and this index will be used to count into the Exodus Struct arrays. */

  /* Setup memory to hold least square projected derivatives for the fluid
     shear stress tensor. For 2D, this will entail storing t11, t12, & t22
     for each node. For 3D, it will be t11 t12, t13, t22, t23, & t33 for
     each node. Access will be (for t12) tau_lsp[0][1][node]. The reason
     this looks a little complicated is because tau is symetric, and we
     don't want to waste space for unneeded components*num_nodes. */
  /* VIM is used here instead of exo->num_dim because for certain 2D models
     that were specified to be solved in cylindrical or swirling coordinate
     systems, they contain 3D tensors (thus tau_* must contain 3D comps) */
  /* MMH: PROJECTED_CARTESIAN coordinate systems are similar to
   * SWIRLING in this respect, too.
   */

  tau_lsp = (double ***)smalloc(VIM * sizeof(double **));
  for (i = 0; i < VIM; i++) {
    tau_lsp[i] = (double **)smalloc(VIM * sizeof(double *));
    for (j = 0; j < VIM; j++) {
      if (j < i) {
        tau_lsp[i][j] = tau_lsp[j][i];
      } else {
        tau_lsp[i][j] = (double *)smalloc(exo->num_nodes * sizeof(double));
        for (k = 0; k < num_local_proc_nodes; k++) {
          tau_lsp[i][j][k] = 0.;
        }
      }
    }
  }

  /* Complete node loop for LS patch construction @ each node */
  for (i_node = 0; i_node < num_local_proc_nodes; i_node++) {

#ifdef RRL_DEBUG
#ifdef DBG_1
    fprintf(stdout, "At node %d and seeing elements:\n", i_node + 1);
#endif
#endif

    if (exo->num_dim == 3) {
      i_node_coords[0] = exo->x_coord[i_node];
      i_node_coords[1] = exo->y_coord[i_node];
      i_node_coords[2] = exo->z_coord[i_node];
    } else {
      i_node_coords[0] = exo->x_coord[i_node];
      i_node_coords[1] = exo->y_coord[i_node];
      i_node_coords[2] = 0.;
    }

    i_start = exo->node_elem_pntr[i_node];
    i_end = exo->node_elem_pntr[i_node + 1];

    num_elems_in_patch = i_end - i_start;

    /* Build mask for this patch and initialize */
    valid_elem_mask = (int *)smalloc(num_elems_in_patch * sizeof(int));
    for (i = 0; i < num_elems_in_patch; i++)
      valid_elem_mask[i] = 0;

    max_gp = -1;
    min_gp = 1000;
    max_dim = -1;
    last_interp = -1;
    j = 0;
    /* Loop over elems in this patch about node i_node
      ( initial scoping loop over this patch ) */
    for (i_elem = i_start; i_elem < i_end; i_elem++) {

      /* Check elem_type for each elem (for appropriate interp)
         and check pd_glob for the block of this elem to check for
         velocity existance and for variable interp for the material.
         This data will be applied to the patch mask to indicate
         which neighbor elems can be appropriately employed in the
         least square fit for this patch. */

      elem_id = exo->node_elem_list[i_elem];
      i_elem_type = Elem_Type(exo, elem_id);     /* element type */
      i_elem_gp = elem_info(NQUAD, i_elem_type); /* number of
                                                    quadrature points */
      i_elem_dim = elem_info(NDIM, i_elem_type); /* element dimension
                                                    (of course!) */
      /* Save this dimension if first time through - subsequent passes
         test if dimension ever changes. This is a no no. */
      if (i_elem == i_start)
        max_dim = i_elem_dim;
      if (i_elem_dim != max_dim || i_elem_dim != exo->num_dim) {
        GOMA_EH(GOMA_ERROR, "Cannot mix element dimensionality for error computation");
      }

      if (i_elem_gp > max_gp) {
        max_gp = i_elem_gp;
      }

      if (i_elem_gp < min_gp) {
        min_gp = i_elem_gp;
      }

      i_eb_indx = exo->elem_eb[elem_id];
      mat_num = Matilda[i_eb_indx];
      pd = pd_glob[mat_num];

      /* Confirm momentum is being solved in this material & velocity is
         turned on for each dimension that the element occupies. */
      valid_count = 0;
      for (k = 0; k < i_elem_dim; k++) {
        if (pd->e[pg->imtrx][R_MOMENTUM1 + k] && pd->v[pg->imtrx][VELOCITY1 + k]) {
          valid_count++;
        }
      }

      /* Only include element if every dimension was valid, and if this
         ev_indx exists for this i_eb_indx */
      if (valid_count == i_elem_dim &&
          exo->elem_var_tab[i_eb_indx * exo->num_elem_vars + ev_indx] == 1) {
        valid_elem_mask[j] = 1;

        /* Save this interpolation if first time through - subsequent passes
         test if interpolation ever changes. This is a no no. */
        if (i_elem == i_start)
          last_interp = pd->i[pg->imtrx][VELOCITY1];
        if (pd->i[pg->imtrx][VELOCITY1] != last_interp) {
          GOMA_EH(GOMA_ERROR, "Cannot mix velocity interpolation levels for error computation");
        }
      }

#ifdef RRL_DEBUG
#ifdef DBG_1
      fprintf(stdout, "       %d (blk %d, type %d, gp %d)\n", elem_id + 1,
              exo->eb_id[exo->elem_eb[elem_id]], i_elem_type, i_elem_gp);
#endif
#endif

      j++;
    } /* End initial scoping loop over this patch about node i_node */

    /* Need to build local arrays to store working local gp coords */
    xgp_loc = (double **)smalloc(num_elems_in_patch * sizeof(double *));
    ygp_loc = (double **)smalloc(num_elems_in_patch * sizeof(double *));
    zgp_loc = (double **)smalloc(num_elems_in_patch * sizeof(double *));
    det_gp_loc = (double **)smalloc(num_elems_in_patch * sizeof(double *));
    wt_gp_loc = (double **)smalloc(num_elems_in_patch * sizeof(double *));

    for (k = 0; k < num_elems_in_patch; k++) {
      xgp_loc[k] = (double *)smalloc(max_gp * sizeof(double));
      ygp_loc[k] = (double *)smalloc(max_gp * sizeof(double));
      zgp_loc[k] = (double *)smalloc(max_gp * sizeof(double));
      det_gp_loc[k] = (double *)smalloc(max_gp * sizeof(double));
      wt_gp_loc[k] = (double *)smalloc(max_gp * sizeof(double));

      /* Initialize the arrays for this patch */
      for (kk = 0; kk < max_gp; kk++) {
        xgp_loc[k][kk] = 0.;
        ygp_loc[k][kk] = 0.;
        zgp_loc[k][kk] = 0.;
        det_gp_loc[k][kk] = 0.;
        wt_gp_loc[k][kk] = 0.;
      }
    }

    /* Determine size of LS system to be solved for this patch */
    if (max_dim > 2 && min_gp > 9) {
      max_terms = 10; /* 3D with quadratic elements */
    } else if (max_dim > 2 && min_gp <= 9) {
      max_terms = 4; /* 3D with linear elements */
    } else if (max_dim <= 2 && min_gp > 5) {
      max_terms = 6; /* 2D with quadratic elements */
    } else {
      max_terms = 3; /* 2D with linear elements */
    }

    s_lhs = (double **)smalloc(max_terms * sizeof(double *));
    for (k = 0; k < max_terms; k++) {
      s_lhs[k] = (double *)smalloc(max_terms * sizeof(double));
      for (kk = 0; kk < max_terms; kk++) {
        s_lhs[k][kk] = 0.;
      }
    }
    rhs = (double *)smalloc(max_terms * sizeof(double));
    for (kk = 0; kk < max_terms; kk++) {
      rhs[kk] = 0.;
    }

    /* Setup memory to hold least squares' projected derivatives for the
       fluid shear stress tensor at the gauss points for elements in this
       patch.  For 2D, this will entail storing t11, t12, & t22 for each
       gauss point. For 3D, it will be t11 t12, t13, t22, t23, & t33 for each
       gauss point. Access will be (for t12) tau_gp_ptch[0][1][elem][gp]. The
       reason this looks a little complicated is because tau is symetric, and
       we don't want to waste space for unneeded components. */
    /* VIM is used here instead of exo->num_dim because for certain 2D models
       that were specified to be solved in cylindrical or swirling coordinate
       systems, they contain 3D tensors (thus tau_* must contain 3D comps) */
    /* MMH: PROJECTED_CARTESIAN coordinate systems are similar to
     * SWIRLING in this respect, too.
     */

    tau_gp_ptch = (double ****)smalloc(VIM * sizeof(double ***));
    for (i = 0; i < VIM; i++) {
      tau_gp_ptch[i] = (double ***)smalloc(VIM * sizeof(double **));
      for (j = 0; j < VIM; j++) {
        if (j < i) {
          tau_gp_ptch[i][j] = tau_gp_ptch[j][i];
        } else {
          tau_gp_ptch[i][j] = (double **)smalloc(num_elems_in_patch * sizeof(double *));
          for (k = 0; k < num_elems_in_patch; k++) {
            tau_gp_ptch[i][j][k] = (double *)smalloc(max_gp * sizeof(double));
            for (kk = 0; kk < max_gp; kk++) {
              tau_gp_ptch[i][j][k][kk] = 0.0;
            }
          }
        }
      }
    }

    i_start = exo->node_elem_pntr[i_node];
    i_end = exo->node_elem_pntr[i_node + 1];
    /* Actual workhorse LHS loop over the elems in this patch */
    for (k = 0, i_elem = i_start; i_elem < i_end; k++, i_elem++) {
      /* Only employ element if it has been deemed worthy from above scoping loop */
      if (valid_elem_mask[k] == 1) {

        /* Extract the info for this elem of the patch */
        elem_id = exo->node_elem_list[i_elem];

        err = load_elem_dofptr(elem_id, exo, x, x_old, xdot, xdot_old, 0);
        GOMA_EH(err, "load_elem_dofptr");

        err = bf_mp_init(pd);
        GOMA_EH(err, "bf_mp_init");

        mat_num = ei[pg->imtrx]->mn;

        i_elem_type = ei[pg->imtrx]->ielem_type;   /* element type */
        i_elem_gp = elem_info(NQUAD, i_elem_type); /* number of
                                                      quadrature points */
        i_elem_dim = ei[pg->imtrx]->ielem_dim;     /* element dimension
                                                     (of course!) */

        err = fill_lhs_lspatch(i_node_coords, wt_gp_loc[k], xgp_loc[k], ygp_loc[k], zgp_loc[k],
                               det_gp_loc[k], max_terms, s_lhs, k, tau_gp_ptch);

      } /* End of treating worthy nodes with valid_elem_mask values */
    }   /* End workhorse LHS loop over this patch */

    /* Now loop over needed components in tau_lsp for node i_node,
       again over the elements in the patch this time filling the rhs
       vector and solving via lu decomp and back sub for each
       component in tau_lsp in turn.*/
    /* NOTE: We're computing the needed tau components using VIM since
       that will pick up the necessary 3D components for cylindrical
       or swirling coordinate systems. The size of the system (least squares)
       being solved will be determined by the max_terms calculation, since
       that mechanism will result in the most high fidelity LS patch fitted.
       For example, if we were to take a 2D mesh of quadratic elements with
       a cylindrical coord system, and call it 3D, the best patch fit it could
       get would be linear 3D. By solving a 2D system on quadratic elements,
       the system fidelity goes up to 6 terms for the same problem. We can
       always solve for the additional tau terms with either system, so we're
       dictating the number of tau terms via VIM, and the max_terms via the
       actual element dimension and interpolation levels. Per conversation with
       RRR, 9/21/98 */
    /* MMH: PROJECTED_CARTESIAN coordinate systems are similar to
     * SWIRLING in this respect, too.
     */
    do_the_lu_decomp = 1;
    indx = (int *)smalloc(max_terms * sizeof(int));
    for (i = 0; i < VIM; i++) {
      for (j = 0; j < VIM; j++) {
        if (j < i) {
          tau_lsp[i][j][i_node] = tau_lsp[j][i][i_node];
        } else {
          for (kk = 0; kk < max_terms; kk++) {
            rhs[kk] = 0.;
          }

          for (k = 0, i_elem = i_start; i_elem < i_end; k++, i_elem++) {
            /* Only employ element if it has been deemed worthy from
             * above scoping loop */
            if (valid_elem_mask[k] == 1) {
              elem_id = exo->node_elem_list[i_elem];
              i_elem_type = Elem_Type(exo, elem_id);     /* element type */
              i_elem_gp = elem_info(NQUAD, i_elem_type); /* number of
                                                            quadrature points */
              for (kk = 0; kk < i_elem_gp; kk++) {
                tau = tau_gp_ptch[i][j][k][kk];
                xgp = xgp_loc[k][kk];
                ygp = ygp_loc[k][kk];
                zgp = zgp_loc[k][kk];
                det = det_gp_loc[k][kk];
                wt = wt_gp_loc[k][kk];
                switch (max_terms) {
                case 3: /* 2D with linear elements    */
                  rhs[0] += 1.0 * tau * det * wt;
                  rhs[1] += xgp * tau * det * wt;
                  rhs[2] += ygp * tau * det * wt;
                  break;
                case 6: /* 2D with quadratic elements */
                  rhs[0] += 1.0 * tau * det * wt;
                  rhs[1] += xgp * tau * det * wt;
                  rhs[2] += ygp * tau * det * wt;
                  rhs[3] += xgp * xgp * tau * det * wt;
                  rhs[4] += xgp * ygp * tau * det * wt;
                  rhs[5] += ygp * ygp * tau * det * wt;
                  break;
                case 4: /* 3D with linear elements    */
                  rhs[0] += 1.0 * tau * det * wt;
                  rhs[1] += xgp * tau * det * wt;
                  rhs[2] += ygp * tau * det * wt;
                  rhs[3] += zgp * tau * det * wt;
                  break;
                case 10: /* 3D with quadratic elements */
                  rhs[0] += 1.0 * tau * det * wt;
                  rhs[1] += xgp * tau * det * wt;
                  rhs[2] += ygp * tau * det * wt;
                  rhs[3] += zgp * tau * det * wt;
                  rhs[4] += xgp * xgp * tau * det * wt;
                  rhs[5] += xgp * ygp * tau * det * wt;
                  rhs[6] += xgp * zgp * tau * det * wt;
                  rhs[7] += ygp * ygp * tau * det * wt;
                  rhs[8] += ygp * zgp * tau * det * wt;
                  rhs[9] += zgp * zgp * tau * det * wt;
                  break;
                default:
                  GOMA_EH(GOMA_ERROR, "Unsupported size in building RHS of least"
                                      " squares patch for error");
                  break;
                }
              }
            }
          }

          /* Now solve the system for this patch. The first time through perform
             the LU decomposition and the back substitution. All subsequent loop
             iterations for this patch need only back substitute. */

          if (lu_decomp_backsub_driver(s_lhs, rhs, indx, max_terms, do_the_lu_decomp) == -1) {
            GOMA_EH(GOMA_ERROR, " Error occurred in calc_zz_error_vel");
            status = -1;
            return (status);
          }

          /* Next time just do the back sub */
          if (do_the_lu_decomp == 1)
            do_the_lu_decomp = 0;

          tau_lsp[i][j][i_node] = rhs[0]; /* This works since we solved the system on a local
                                             coord basis */

#ifdef RRL_DEBUG
#ifdef DBG_2
          fprintf(stdout, "Just computed lsp tau%d%d for node %d = %6.4lf\n", i + 1, j + 1,
                  i_node + 1, tau_lsp[i][j][i_node]);
#endif
#endif
        }
      }
    }

#ifdef RRL_DEBUG
#ifdef DBG_1
    fprintf(stdout, "\n  node    t11       t12       t13       t21       t22       t23       t31   "
                    "   t32      t33\n");
    if (ei[pg->imtrx]->ielem_dim > 2) {
      fprintf(stdout,
              "  %2d   %6.4lf   %6.4lf   %6.4lf   %6.4lf   %6.4lf   %6.4lf   %6.4lf   %6.4lf   "
              "%6.4lf\n",
              i_node + 1, tau_lsp[0][0][i_node], tau_lsp[0][1][i_node], tau_lsp[0][2][i_node],
              tau_lsp[1][0][i_node], tau_lsp[1][1][i_node], tau_lsp[1][2][i_node],
              tau_lsp[2][0][i_node], tau_lsp[2][1][i_node], tau_lsp[2][2][i_node]);
    } else {
      fprintf(stdout,
              "  %2d   %6.4lf   %6.4lf   < NA >   %6.4lf   %6.4lf   < NA >   < NA >   < NA >   < "
              "NA >\n",
              i_node + 1, tau_lsp[0][0][i_node], tau_lsp[0][1][i_node], tau_lsp[1][0][i_node],
              tau_lsp[1][1][i_node]);
    }
#endif
#endif

    /* Clean up memory for this patch */
    for (k = 0; k < num_elems_in_patch; k++) {
      free(xgp_loc[k]);
      free(ygp_loc[k]);
      free(zgp_loc[k]);
      free(det_gp_loc[k]);
      free(wt_gp_loc[k]);
    }
    free(xgp_loc);
    free(ygp_loc);
    free(zgp_loc);
    free(det_gp_loc);
    free(wt_gp_loc);

    free(valid_elem_mask);

    for (k = 0; k < max_terms; k++) {
      free(s_lhs[k]);
    }
    free(s_lhs);
    free(rhs);
    free(indx);

    for (i = 0; i < VIM; i++) {
      for (j = 0; j < VIM; j++) {
        /* Only the upper tri diag members were malloc'd, so only free those */
        if (j >= i) {
          for (k = 0; k < num_elems_in_patch; k++) {
            free(tau_gp_ptch[i][j][k]);
          }
          free(tau_gp_ptch[i][j]);
        }
      }
      free(tau_gp_ptch[i]);
    }
    free(tau_gp_ptch);

#ifdef RRL_DEBUG
#ifdef DBG_1
    fprintf(stdout, "\n");
#endif
#endif

  } /* End of loop building LHS for the patch about node i_node */

  /* Loop over each element in the model for the following. */
  /* Now have tau_lsp at the nodes, and tau at the gauss points (via fv)
     To get tau_lsp back at the gauss points, we will interpolate with
     shape functions. With both quantities back at the gauss points, we will
     cycle over the gauss points, computing the gfem/lsp difference for
     each tau component at each gauss point, then using these to compute
     the elemental error. Then store this absolute error val for the
     currrent element into gvec_elem for now. */

  valid_elem_mask = (int *)smalloc(exo->num_elems * sizeof(int));
  for (i = 0; i < exo->num_elems; i++)
    valid_elem_mask[i] = 0;

  /* First do scoping loop to set up the valid_elem_mask for the entire model
     this time (valid to compute zz error based on velocity)  */
  for (i_elem = 0; i_elem < exo->num_elems; i_elem++) {
    elem_id = exo->node_elem_list[i_elem];
    i_elem_type = Elem_Type(exo, elem_id);     /* element type */
    i_elem_gp = elem_info(NQUAD, i_elem_type); /* number of quadrature points */
    i_elem_dim = elem_info(NDIM, i_elem_type); /* element dimension  */
    i_eb_indx = exo->elem_eb[i_elem];
    mat_num = Matilda[i_eb_indx];
    pd = pd_glob[mat_num];

    /* Confirm momentum is being solved in this material & velocity is
       turned on for each dimension that the element occupies. */
    valid_count = 0;
    for (k = 0; k < i_elem_dim; k++) {
      if (pd->e[pg->imtrx][R_MOMENTUM1 + k] && pd->v[pg->imtrx][VELOCITY1 + k]) {
        valid_count++;
      }
    }

    /* Only include element if every dimension was valid, and if this
       ev_indx exists for this i_eb_indx */
    if (valid_count == i_elem_dim &&
        exo->elem_var_tab[i_eb_indx * exo->num_elem_vars + ev_indx] == 1) {
      valid_elem_mask[i_elem] = 1;
    }
  } /* End scoping loop for element error integration */

  /* Now do real loop for element error integration over entire model */
  /* Set up values to find max's and do the elem size stuff */
  if (compute_elem_size != 0) {
    elem_areas = (double *)smalloc(exo->num_elems * sizeof(double));
    for (i = 0; i < exo->num_elems; i++)
      elem_areas[i] = 0;
  }
  max_velocity_norm = 0.;
  max_velocity_err = 0.;
  total_volume = 0.;
  pct_over_target = 0.;
  pct_over_tolerance = 0.;
  volume_tmp = 0.;
  remesh_status = -1;

  for (i_elem = 0; i_elem < exo->num_elems; i_elem++) {
    /* Only employ element if it has been deemed worthy from above scoping loop */

    if (valid_elem_mask[i_elem] == 1) {
      err = load_elem_dofptr(i_elem, exo, x, x_old, xdot, xdot_old, 0);
      GOMA_EH(err, "load_elem_dofptr");

      err = bf_mp_init(pd);
      GOMA_EH(err, "bf_mp_init");

      /* Extract the info for this elem of the patch */

      mat_num = ei[pg->imtrx]->mn;
      i_elem_type = ei[pg->imtrx]->ielem_type;   /* element type */
      i_elem_gp = elem_info(NQUAD, i_elem_type); /* number of quadrature points */
      i_elem_dim = ei[pg->imtrx]->ielem_dim;     /* element dimension
                                                   (of course!) */

      i_eb_indx = exo->elem_eb[i_elem];

#ifdef RRL_DEBUG
#ifdef DBG_1
      fprintf(stdout, "\n\n About to integrate zz error on the worthy element %d (blk %d)\n",
              i_elem + 1, exo->eb_id[i_eb_indx]);
#endif
#endif

      /* Now interpolate the LS projected values of tau back to the gauss points,
         extract the fem values of tau at the same gauss points,
         and compute the error (absolute indicator) for elem i_elem */
      err = abs_error_at_elem(i_elem, tau_lsp, gvec_elem[i_eb_indx][ev_indx],
                              &max_velocity_norm_tmp, &elem_areas[i_elem]);

      /* sum the total area/vol of the elems being used for error computation */
      total_volume += elem_areas[i_elem];

      if (max_velocity_norm_tmp > max_velocity_norm) {
        max_velocity_norm = max_velocity_norm_tmp;
        max_velocity_norm_i_elem = i_elem;
      }

      if (gvec_elem[i_eb_indx][ev_indx][i_elem] > max_velocity_err) {
        max_velocity_err = gvec_elem[i_eb_indx][ev_indx][i_elem];
        max_velocity_err_i_elem = i_elem;

        max_x = 0;
        max_y = 0;
        max_z = 0;

        for (j = 0; j < ei[pg->imtrx]->num_local_nodes; j++) {
          /* also convert from node number to dof number */
          global_node_num = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + j];
          max_x += exo->x_coord[global_node_num];
          max_y += exo->y_coord[global_node_num];
          if (exo->num_dim > 2)
            max_z += exo->z_coord[global_node_num];
        }
        max_x /= ei[pg->imtrx]->num_local_nodes;
        max_y /= ei[pg->imtrx]->num_local_nodes;
        if (exo->num_dim > 2)
          max_z /= ei[pg->imtrx]->num_local_nodes;
      }

    } /* End of treating worthy elems with valid_elem_mask values */
  }   /* End of real loop for element error integration over entire model */

#ifdef RRL_DEBUG
#ifdef DBG_1
  fprintf(stdout, "\n\n About to integrate zz error on the worthy element %d (blk %d)\n",
          i_elem + 1, exo->eb_id[i_eb_indx]);
#endif
#endif

  fprintf(stdout, "\n Max energy norm of the velocity vector ( element %d ) \t= %6.4f\n",
          max_velocity_norm_i_elem + 1, max_velocity_norm);
  fprintf(stdout, " Max ZZ error (velocity based)          ( element %d ) \t= %6.4f\n",
          max_velocity_err_i_elem + 1, max_velocity_err);
  fprintf(stdout, "   x centroid = %6.4f\n", max_x);
  fprintf(stdout, "   y centroid = %6.4f\n", max_y);
  fprintf(stdout, "   z centroid = %6.4f\n", max_z);

  /* Now compute recommended elem sizes if needed */
  if (compute_elem_size != 0) {
    /* First grab user supplied params */
    reduction_rate = pp_error_data[0].error_params[0];
    expansion_rate = pp_error_data[0].error_params[1];
    elem_size_min = pp_error_data[0].error_params[2];
    elem_size_max = pp_error_data[0].error_params[3];
    target_error = pp_error_data[0].error_params[4];
    pct_over_tolerance = pp_error_data[0].error_params[5];

    for (i_elem = 0; i_elem < exo->num_elems; i_elem++) {
      /* Only employ element if it has been deemed worthy from above scoping loop */
      if (valid_elem_mask[i_elem] == 1) {
        i_eb_indx = exo->elem_eb[i_elem];

        /* Determine typical element length scale - sqrt if 2D, cube root if 3D */
        if (exo->num_dim == 3) {
          h1 = pow(elem_areas[i_elem], 0.3333333333);
        } else {
          h1 = sqrt(elem_areas[i_elem]);
        }

        error_ratio = gvec_elem[i_eb_indx][ev_indx][i_elem] / target_error;

        if (error_ratio <= 1.0) {
          h2 = (pow(error_ratio, -1. / expansion_rate)) * h1;
        } else {
          h2 = (pow(error_ratio, -1. / reduction_rate)) * h1;
          /* This element needs refinement - add into the volume over calculation */
          volume_tmp += elem_areas[i_elem];
        }

        /* Now ensure that the elem size floors & ceilings are observed */
        if (h2 < elem_size_min)
          h2 = elem_size_min;
        if (h2 > elem_size_max)
          h2 = elem_size_max;

        /* Store into gvec array - note that we store into the next ev_indx place
           since this value always follows the error measure that it applies to */
        gvec_elem[i_eb_indx][ev_indx + 1][i_elem] = h2;
      }
    }

    pct_over_target = volume_tmp / (total_volume * 0.01); /* Builds % of volume over error
                                                            target */

    fprintf(stdout, "\n Target error                                        \t= %6.4f\n",
            target_error);
    fprintf(stdout, " Problem volume (area) considered for error          \t= %6.4f\n",
            total_volume);
    fprintf(stdout, " Percent of problem volume over error target         \t= %6.2f [%%]\n",
            pct_over_target);
    fprintf(stdout, " Tolerance for percent of problem volume over value  \t= %6.2f [%%]\n",
            pct_over_tolerance);

    if (pct_over_target > pct_over_tolerance) {
      /* Need to remesh the problem and run it again using sizing data */
      remesh_status = 2;
      fprintf(stdout, " Remeshing status                                    \t= %d\n",
              remesh_status);
      fprintf(stdout, "\n (Recommend remeshing model using the *ELSIZE element variable)\n");
    } else if (pct_over_target <= pct_over_tolerance && max_velocity_err > target_error) {
      remesh_status = 1;
      fprintf(stdout, " Remeshing status                                    \t= %d\n",
              remesh_status);
      fprintf(stdout,
              "\n (Max error > target error, but volume < target error within tolerance)\n");
    } else if (max_velocity_err <= target_error) {
      remesh_status = 0;
      fprintf(stdout, " Remeshing status                                    \t= %d\n",
              remesh_status);
      fprintf(stdout, "\n (Max error < target error, error reduction objective met)\n");
    } else {
      GOMA_EH(GOMA_ERROR,
              "Undetermined state encountered during error/remeshing size data calculation\n");
    }

    fprintf(stdout, "\n Current number of nodes (resolution)                \t= %d\n",
            exo->num_nodes);
    fprintf(stdout, " Current number of elements                          \t= %d\n",
            exo->num_elems);

    /* Free memory */
    free(elem_areas);
  }

  /* Now normalize the raw values of ZZ velocity based error */
  for (i_elem = 0; i_elem < exo->num_elems; i_elem++) {
    /* Only employ element if it has been deemed worthy from above scoping loop */
    if (valid_elem_mask[i_elem] == 1) {
      i_eb_indx = exo->elem_eb[i_elem];
      gvec_elem[i_eb_indx][ev_indx][i_elem] /=
          (max_velocity_norm * 0.01); /* NOTE: converting to % */
    }
  }

  /* Final memory cleanup */
  for (i = 0; i < VIM; i++) {
    for (j = 0; j < VIM; j++) {
      /* Only the upper tri diag members were malloc'd, so only free those */
      if (j >= i) {
        free(tau_lsp[i][j]);
      }
    }
    free(tau_lsp[i]);
  }
  free(tau_lsp);
  free(i_node_coords);
  free(valid_elem_mask);

  return (status);
}

static int abs_error_at_elem(int i_elem,
                             double ***tau_nodal_lsp,
                             double *gvec_elem,
                             double *velocity_norm,
                             double *ielem_area)
/******************************************************************************
  Function which interpolates the LS projected values of tau back to the gauss
  points, extracts the fem values of tau at the same gauss points, then using
  the differences between the two, computes the velocity based error norm
  (the absolute indicator value) for the element.

  Author:          R. R. Lober (9113)
  Date:            30 September 1998
  Revised

******************************************************************************/
{
  int status = 0;
  int i, j, i_elem_gp, eqn, global_node_num, a, b;
  int err; /* temp variable to hold diagnostic flags. */
  double sum_tau, phi_i, ***tau_gp_fem, ***tau_gp_lsp, ***tau_gp_ext;
  double xi[DIM]; /* Local element coordinates of Gauss point. */
  double det, wt;
  double gamma[DIM][DIM]; /* shrearrate tensor based on velocity */

  double mu, elem_error, error_norm_tmp, elem_area, velocity_norm_tmp, elem_vel_norm;

  /* polymer viscosity */
  double mup;

  i_elem_gp = elem_info(NQUAD, ei[pg->imtrx]->ielem_type); /* number of
                                                   quadrature points */
  eqn = R_MOMENTUM1; /* We depend on this eqn for the velocity based error measure */

  /* Grab memory to hold the fem & lsp values of tau, only malloc the upper tri-diag
     portion of these symetric tensors */
  tau_gp_fem = (double ***)smalloc(VIM * sizeof(double **));
  tau_gp_lsp = (double ***)smalloc(VIM * sizeof(double **));
  tau_gp_ext = (double ***)smalloc(VIM * sizeof(double **));
  for (i = 0; i < VIM; i++) {
    tau_gp_fem[i] = (double **)smalloc(VIM * sizeof(double *));
    tau_gp_lsp[i] = (double **)smalloc(VIM * sizeof(double *));
    tau_gp_ext[i] = (double **)smalloc(VIM * sizeof(double *));
    for (j = 0; j < VIM; j++) {
      if (j < i) {
        tau_gp_fem[i][j] = tau_gp_fem[j][i];
        tau_gp_lsp[i][j] = tau_gp_lsp[j][i];
        tau_gp_ext[i][j] = tau_gp_ext[j][i];
      } else {
        tau_gp_fem[i][j] = (double *)smalloc(i_elem_gp * sizeof(double));
        tau_gp_lsp[i][j] = (double *)smalloc(i_elem_gp * sizeof(double));
        tau_gp_ext[i][j] = (double *)smalloc(i_elem_gp * sizeof(double));
      }
    }
  }

  elem_error = 0.;
  elem_vel_norm = 0.;
  elem_area = 0.;

  /* Now loop over all the gauss points for this elem */
  for (i = 0; i < i_elem_gp; i++) {
    find_stu(i, ei[pg->imtrx]->ielem_type, &xi[0], &xi[1], &xi[2]); /* find quadrature point */

    fv->wt = Gq_weight(i, ei[pg->imtrx]->ielem_type); /* find quadrature weights for
                                              current ip */
    wt = fv->wt;

    err = load_basis_functions(xi, bfd);
    GOMA_EH(err, "problem from load_basis_functions");

    /* This has elemental Jacobian transformation and some
       basic mesh derivatives...
       Old usage: calc_Jac, jelly_belly */
    /* NOTE: this call also updates the value of detJ inside the bf struct
       for the point xi[] in the element pointed to by the ei struct */

    err = beer_belly();
    GOMA_EH(err, "beer_belly");

    /* Load up field variable values at this Gauss point. */
    /* Now fill fv with tau goodies, and cycle over the gauss points
       filling the tau components into tau_gp_fem */
    err = load_fv();
    GOMA_EH(err, "load_fv");

    /* Load up physical space gradients of field variables at this
       Gauss point. */
    /* NOTE: load_bf_grad MUST be called before load_fv_grads as this
       call depends on it! - RRL 10/30/98 */
    err = load_bf_grad();
    GOMA_EH(err, "load_bf_grad");

    err = load_fv_vector();

    err = load_fv_grads();
    GOMA_EH(err, "load_fv_grads");

    /* Now generate & save the tau shear stress tensor components for this
       gauss point. */
    /* In Cartesian coordinates, this velocity gradient tensor will
       have components that are...

                        grad_v[a][b] = d v_b
                                       -----
                                       d x_a                         */

    for (a = 0; a < VIM; a++) {
      for (b = 0; b < VIM; b++) {
        gamma[a][b] = fv->grad_v[a][b] + fv->grad_v[b][a];
      }
    }

    mu = viscosity(gn, gamma, NULL);

    if (pd->v[pg->imtrx][POLYMER_STRESS11]) {
      /* get polymer viscosity */
      mup = viscosity(gn, gamma, NULL);
      mu = mu + mup;
    }

    for (a = 0; a < VIM; a++) {
      for (b = 0; b < VIM; b++) {
        tau_gp_fem[a][b][i] = mu * gamma[a][b];
      }
    }

    det = bf[eqn]->detJ;

    /* Now have all fem tau components for this elem, and gauss point i.
       Now interpolate the tau_lsp values back to gauss point i as well. */
    for (a = 0; a < VIM; a++) {
      for (b = 0; b < VIM; b++) {
        /* Calculate this symetrically (that's how its being stored) */
        if (b >= a) {
          sum_tau = 0.;
          for (j = 0; j < ei[pg->imtrx]->num_local_nodes; j++) {
            /* also convert from node number to dof number */
            phi_i = bf[eqn]->phi[ei[pg->imtrx]->ln_to_dof[eqn][j]];
            global_node_num = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + j];
            sum_tau += phi_i * tau_nodal_lsp[a][b][global_node_num];
          }
          tau_gp_lsp[a][b][i] = sum_tau;
#ifdef RRL_DEBUG
#ifdef DBG_1
          fprintf(stdout,
                  "          LS \ttau%d%d at gp %d = %8.6lf (x=%6.4lf\ty=%6.4lf\tz=%6.4lf)\n",
                  a + 1, b + 1, i, tau_gp_lsp[a][b][i], sumx, sumy, sumz);
#endif
#endif
        }
      }
    }

    /* Now integrate the error over the elem - this is essentially

       -- for 2D ---------------
       = {(tau11)**2 + 2(tau12)**2 + (tau22)**2}

       -- for 3D ---------------
       = {(tau11)**2 + 2(tau12)**2 + 2(tau13)**2 + (tau22)**2 +
       2(tau23)**2 + (tau33)**2}

       Note that each off diagonal term occurs twice    */

    error_norm_tmp = 0.;
    velocity_norm_tmp = 0.;
    for (a = 0; a < VIM; a++) {
      for (b = 0; b < VIM; b++) {
        if (b > a) { /* Off axis symmetric term (upper tri only) */
          error_norm_tmp += 2 * ((tau_gp_lsp[a][b][i] - tau_gp_fem[a][b][i]) *
                                 (tau_gp_lsp[a][b][i] - tau_gp_fem[a][b][i]));
          velocity_norm_tmp += 2 * (tau_gp_fem[a][b][i] * tau_gp_fem[a][b][i]);
#ifdef RRL_DEBUG
#ifdef DBG_1
          fprintf(stdout,
                  "\n    computed off diag error component \t2*(tau%d%d)**2 at gp %d = %8.6lf\n",
                  a + 1, b + 1, i,
                  2 * ((tau_gp_lsp[a][b][i] - tau_gp_fem[a][b][i]) *
                       (tau_gp_lsp[a][b][i] - tau_gp_fem[a][b][i])));
          fprintf(stdout,
                  "    computed LS             component \t(tau_gp_lsp%d%d)  gp %d = %8.6lf\n",
                  a + 1, b + 1, i, tau_gp_lsp[a][b][i]);
          fprintf(stdout,
                  "    computed FEM            component \t(tau_gp_fem%d%d)  gp %d = %8.6lf\n",
                  a + 1, b + 1, i, tau_gp_fem[a][b][i]);

          fprintf(stdout,
                  "    computed EXACT          component \t(tau%d%d exact)   gp %d = %8.6lf\n",
                  a + 1, b + 1, i, tau_gp_ext[a][b][i]);

          fprintf(stdout, "    gp %d is at \tx=%6.4lf\ty=%6.4lf\tz=%6.4lf\n", i, sumx, sumy, sumz);
          fprintf(stdout,
                  "    computed FEM          u component \t                gp %d = %8.6lf\n", i,
                  fv->v[0]);
          fprintf(stdout,
                  "    computed EXACT        u component \t                gp %d = %8.6lf\n", i,
                  2. - 2. * (sumy * sumy));
#endif
#endif
        } else if (b == a) { /* diagonal term */
          error_norm_tmp += ((tau_gp_lsp[a][b][i] - tau_gp_fem[a][b][i]) *
                             (tau_gp_lsp[a][b][i] - tau_gp_fem[a][b][i]));
          velocity_norm_tmp += (tau_gp_fem[a][b][i] * tau_gp_fem[a][b][i]);
#ifdef RRL_DEBUG
#ifdef DBG_1
          fprintf(stdout,
                  "\n    computed     diag error component \t  (tau%d%d)**2 at gp %d = %8.6lf\n",
                  a + 1, b + 1, i,
                  ((tau_gp_lsp[a][b][i] - tau_gp_fem[a][b][i]) *
                   (tau_gp_lsp[a][b][i] - tau_gp_fem[a][b][i])));
          fprintf(stdout,
                  "    computed LS             component \t(tau_gp_lsp%d%d)  gp %d = %8.6lf\n",
                  a + 1, b + 1, i, tau_gp_lsp[a][b][i]);
          fprintf(stdout,
                  "    computed FEM            component \t(tau_gp_fem%d%d)  gp %d = %8.6lf\n",
                  a + 1, b + 1, i, tau_gp_fem[a][b][i]);
          fprintf(stdout,
                  "    computed EXACT          component \t(tau%d%d exact)   gp %d = %8.6lf\n",
                  a + 1, b + 1, i, tau_gp_ext[a][b][i]);
          fprintf(stdout, "    gp %d is at \tx=%6.4lf\ty=%6.4lf\tz=%6.4lf\n", i, sumx, sumy, sumz);
          fprintf(stdout,
                  "    computed FEM          u component \t                gp %d = %8.6lf\n", i,
                  fv->v[0]);
          fprintf(stdout,
                  "    computed EXACT        u component \t                gp %d = %8.6lf\n", i,
                  2. - 2. * (sumy * sumy));
#endif
#endif
        }
      }
    }
    elem_error += error_norm_tmp * det * wt;
    elem_vel_norm += velocity_norm_tmp * det * wt;
    elem_area += det * wt;
    ;
  }

  *ielem_area = elem_area;

  /* Now its almost revealed!! */
  gvec_elem[i_elem] = sqrt(elem_error / elem_area);
  *velocity_norm = sqrt(elem_vel_norm / elem_area);

#ifdef RRL_DEBUG
#ifdef DBG_1
  fprintf(stdout, "\n==> ZZ velocity err for element %d = %8.6lf <==\n", ei[pg->imtrx]->ielem + 1,
          gvec_elem[i_elem]);
  fprintf(stdout, "    Velocity norm for element %d = %8.6lf <==\n", ei[pg->imtrx]->ielem + 1,
          *velocity_norm);
#endif
#endif

  for (i = 0; i < VIM; i++) {
    for (j = 0; j < VIM; j++) {
      /* Only the upper tri diag members were malloc'd, so only free those */
      if (j >= i) {
        free(tau_gp_fem[i][j]);
        free(tau_gp_lsp[i][j]);
      }
    }
    free(tau_gp_fem[i]);
    free(tau_gp_lsp[i]);
  }
  free(tau_gp_fem);
  free(tau_gp_lsp);

  return (status);
}

static int fill_lhs_lspatch(double *i_node_coords,
                            double *wt_gp_loc,
                            double *xgp_loc,
                            double *ygp_loc,
                            double *zgp_loc,
                            double *det_gp_loc,
                            int max_terms,
                            double **s_lhs,
                            int elem,
                            double ****tau_gp_ptch)
/******************************************************************************
  Function which calculates the local coords for an element participating in a
  patch about a given node, then fills the contributions from that element into
  the LHS matrix of the least squares patch system being formed about the given
  node.

  Author:          R. R. Lober (9113)
  Date:            15 September 1998
  Revised

******************************************************************************/

{
  int i, j, i_elem_gp, eqn, local_i, a, b;
  int status = 0;
  int err; /* temp variable to hold diagnostic flags. */
  double sumx, sumy, sumz, phi_i;
  double xi[DIM]; /* Local element coordinates of Gauss point. */
  double xgp, ygp, zgp, det, wt;
  double gamma[DIM][DIM]; /* shrearrate tensor based on velocity */

  /* viscosity */
  double mu;

  /* polymer viscosity */
  double mup;

  i_elem_gp = elem_info(NQUAD, ei[pg->imtrx]->ielem_type); /* number of
                                                   quadrature points */
  eqn = R_MOMENTUM1; /* We depend on this eqn for the velocity based error measure */

  /* First build local coord system for the gauss points of this
     element about the node in the center of this patch. This is important
     since these least-squares systems tend to be ill-conditioned when the
     usual poloynomial basis are used (Reference D. Pelletier's report n0
     1, "Implementation of Error Analysis and Norms to Computational Fluid
     Dynamics Applications, (Bilinear Finite Elements)" for SNL, dated
     6/97. (or Randy Lober for copy)). The introduction of a local coord
     system alleviates this ill-conditioning. */

  for (i = 0; i < i_elem_gp; i++) {
    find_stu(i, ei[pg->imtrx]->ielem_type, &xi[0], &xi[1], &xi[2]); /* find quadrature point */

    wt_gp_loc[i] = Gq_weight(i, ei[pg->imtrx]->ielem_type); /* find quadrature weights for
                                                    current ip */
    fv->wt = wt_gp_loc[i];

    err = load_basis_functions(xi, bfd);
    GOMA_EH(err, "problem from load_basis_functions");

    /*
      This has elemental Jacobian transformation and some
      basic mesh derivatives...
      Old usage: calc_Jac, jelly_belly */

    /* NOTE: this call also updates the value of detJ inside the bf struct
       for the point xi[] in the element pointed to by the ei struct */

    err = beer_belly();
    GOMA_EH(err, "beer_belly");

    /*
     * Load up field variable values at this Gauss point.
     */
    /* Now fill fv with tau goodies, and cycle over the gauss points
       filling the tau components into tau_gp_ptch */
    err = load_fv();
    GOMA_EH(err, "load_fv");

    /*
     * Load up physical space gradients of field variables at this
     * Gauss point.
     */
    /* NOTE: load_bf_grad MUST be called before load_fv_grads as this
       call depends on it! - RRL 10/30/98 */
    err = load_bf_grad();
    GOMA_EH(err, "load_bf_grad");

    err = load_fv_vector();

    err = load_fv_grads();
    GOMA_EH(err, "load_fv_grads");

    /* Now generate & save the tau shear stress tensor components for this
       gauss point. */
    /*
     * In Cartesian coordinates, this velocity gradient tensor will
     * have components that are...
     *
     * 			grad_v[a][b] = d v_b
     *				       -----
     *				       d x_a
     */
    for (a = 0; a < VIM; a++) {
      for (b = 0; b < VIM; b++) {
        gamma[a][b] = fv->grad_v[a][b] + fv->grad_v[b][a];
      }
    }
    mu = viscosity(gn, gamma, NULL);

    if (pd->v[pg->imtrx][POLYMER_STRESS11]) {
      /* get polymer viscosity */
      mup = viscosity(gn, gamma, NULL);
      mu = mu + mup;
    }

    for (a = 0; a < VIM; a++) {
      for (b = 0; b < VIM; b++) {
        tau_gp_ptch[a][b][elem][i] = mu * gamma[a][b];
      }
    }

#ifdef RRL_DEBUG
#ifdef DBG_2
    if (i == 0) {
      fprintf(stdout, "\n   gp     t11      t12      t13      t21      t22      t23      t31      "
                      "t32      t33      mu\n");
    }
    if (ei[pg->imtrx]->ielem_dim > 2) {
      fprintf(stdout,
              "  %2d   %6.4lf   %6.4lf   %6.4lf   %6.4lf   %6.4lf   %6.4lf   %6.4lf   %6.4lf   "
              "%6.4lf   %6.4lf\n",
              i + 1, tau_gp_ptch[0][0][elem][i], tau_gp_ptch[0][1][elem][i],
              tau_gp_ptch[0][2][elem][i], tau_gp_ptch[1][0][elem][i], tau_gp_ptch[1][1][elem][i],
              tau_gp_ptch[1][2][elem][i], tau_gp_ptch[2][0][elem][i], tau_gp_ptch[2][1][elem][i],
              tau_gp_ptch[2][2][elem][i], mu);
    } else {
      fprintf(stdout,
              "  %2d   %6.4lf   %6.4lf   < NA >   %6.4lf   %6.4lf   < NA >   < NA >   < NA >   < "
              "NA >   %6.4lf\n",
              i + 1, tau_gp_ptch[0][0][elem][i], tau_gp_ptch[0][1][elem][i],
              tau_gp_ptch[1][0][elem][i], tau_gp_ptch[1][1][elem][i], mu);
    }
#endif
#endif

    /* Now generate & save the relative coords between this gauss point and
       this patch node */
    sumx = 0.;
    sumy = 0.;
    sumz = 0.;
    for (j = 0; j < ei[pg->imtrx]->num_local_nodes; j++) {
      /* also convert from node number to dof number */
      phi_i = bf[eqn]->phi[ei[pg->imtrx]->ln_to_dof[eqn][j]];
      local_i = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + j];
      sumx += phi_i * Coor[0][local_i];
      sumy += phi_i * Coor[1][local_i];
      if (ei[pg->imtrx]->ielem_dim > 2) {
        sumz += phi_i * Coor[2][local_i];
      }
    }
    xgp_loc[i] = sumx - i_node_coords[0];
    ygp_loc[i] = sumy - i_node_coords[1];
    if (ei[pg->imtrx]->ielem_dim > 2) {
      zgp_loc[i] = sumz - i_node_coords[2];
    }
    det_gp_loc[i] = bf[eqn]->detJ;
  }

#ifdef RRL_DEBUG
#ifdef DBG_2
  fprintf(stdout, "   gp     x        y       z        det       wt\n");
  for (i = 0; i < i_elem_gp; i++) {
    fprintf(stdout, "  %2d   %6.4lf   %6.4lf   %6.4lf   %6.4lf   %6.4lf\n", i + 1, xgp_loc[i],
            ygp_loc[i], zgp_loc[i], det_gp_loc[i], wt_gp_loc[i]);
  }
#endif
#endif

  /* Now Populate the contributions to s_lhs from this element
     for each gauss point */
  for (i = 0; i < i_elem_gp; i++) {
    switch (max_terms) {
    case 3: /* 2D with linear elements    */
#ifdef RRL_DEBUG
#ifdef DBG_2
      fprintf(stdout, "fill_lhs_lspatch using max_terms = %d (2D linear case)\n", max_terms);
#endif
#endif
      xgp = xgp_loc[i];
      ygp = ygp_loc[i];
      det = det_gp_loc[i];
      wt = wt_gp_loc[i];
      s_lhs[0][0] += 1.0 * det * wt;
      s_lhs[0][1] += xgp * det * wt;
      s_lhs[0][2] += ygp * det * wt;

      s_lhs[1][0] += xgp * det * wt;
      s_lhs[1][1] += xgp * xgp * det * wt;
      s_lhs[1][2] += xgp * ygp * det * wt;

      s_lhs[2][0] += ygp * det * wt;
      s_lhs[2][1] += xgp * ygp * det * wt;
      s_lhs[2][2] += ygp * ygp * det * wt;
      break;
    case 6: /* 2D with quadratic elements */
#ifdef RRL_DEBUG
#ifdef DBG_2
      fprintf(stdout, "fill_lhs_lspatch using max_terms = %d (2D quadratic case)\n", max_terms);
#endif
#endif
      xgp = xgp_loc[i];
      ygp = ygp_loc[i];
      det = det_gp_loc[i];
      wt = wt_gp_loc[i];
      s_lhs[0][0] += 1.0 * det * wt;
      s_lhs[0][1] += xgp * det * wt;
      s_lhs[0][2] += ygp * det * wt;
      s_lhs[0][3] += xgp * xgp * det * wt;
      s_lhs[0][4] += xgp * ygp * det * wt;
      s_lhs[0][5] += ygp * ygp * det * wt;

      s_lhs[1][0] += xgp * det * wt;
      s_lhs[1][1] += xgp * xgp * det * wt;
      s_lhs[1][2] += xgp * ygp * det * wt;
      s_lhs[1][3] += xgp * xgp * xgp * det * wt;
      s_lhs[1][4] += xgp * xgp * ygp * det * wt;
      s_lhs[1][5] += xgp * ygp * ygp * det * wt;

      s_lhs[2][0] += ygp * det * wt;
      s_lhs[2][1] += ygp * xgp * det * wt;
      s_lhs[2][2] += ygp * ygp * det * wt;
      s_lhs[2][3] += ygp * xgp * xgp * det * wt;
      s_lhs[2][4] += ygp * ygp * xgp * det * wt;
      s_lhs[2][5] += ygp * ygp * ygp * det * wt;

      s_lhs[3][0] += xgp * xgp * det * wt;
      s_lhs[3][1] += xgp * xgp * xgp * det * wt;
      s_lhs[3][2] += xgp * xgp * ygp * det * wt;
      s_lhs[3][3] += xgp * xgp * xgp * xgp * det * wt;
      s_lhs[3][4] += xgp * xgp * xgp * ygp * det * wt;
      s_lhs[3][5] += xgp * xgp * ygp * ygp * det * wt;

      s_lhs[4][0] += xgp * ygp * det * wt;
      s_lhs[4][1] += xgp * xgp * ygp * det * wt;
      s_lhs[4][2] += xgp * ygp * ygp * det * wt;
      s_lhs[4][3] += xgp * xgp * xgp * ygp * det * wt;
      s_lhs[4][4] += xgp * xgp * ygp * ygp * det * wt;
      s_lhs[4][5] += xgp * ygp * ygp * ygp * det * wt;

      s_lhs[5][0] += ygp * ygp * det * wt;
      s_lhs[5][1] += ygp * ygp * xgp * det * wt;
      s_lhs[5][2] += ygp * ygp * ygp * det * wt;
      s_lhs[5][3] += ygp * ygp * xgp * xgp * det * wt;
      s_lhs[5][4] += ygp * ygp * ygp * xgp * det * wt;
      s_lhs[5][5] += ygp * ygp * ygp * ygp * det * wt;
      break;
    case 4: /* 3D with linear elements    */
#ifdef RRL_DEBUG
#ifdef DBG_2
      fprintf(stdout, "fill_lhs_lspatch using max_terms = %d (3D linear case)\n", max_terms);
#endif
#endif
      xgp = xgp_loc[i];
      ygp = ygp_loc[i];
      zgp = zgp_loc[i];
      det = det_gp_loc[i];
      wt = wt_gp_loc[i];
      s_lhs[0][0] += 1.0 * det * wt;
      s_lhs[0][1] += xgp * det * wt;
      s_lhs[0][2] += ygp * det * wt;
      s_lhs[0][3] += zgp * det * wt;

      s_lhs[1][0] += xgp * det * wt;
      s_lhs[1][1] += xgp * xgp * det * wt;
      s_lhs[1][2] += xgp * ygp * det * wt;
      s_lhs[1][3] += xgp * zgp * det * wt;

      s_lhs[2][0] += ygp * det * wt;
      s_lhs[2][1] += ygp * xgp * det * wt;
      s_lhs[2][2] += ygp * ygp * det * wt;
      s_lhs[2][3] += ygp * zgp * det * wt;

      s_lhs[3][0] += zgp * det * wt;
      s_lhs[3][1] += zgp * xgp * det * wt;
      s_lhs[3][2] += zgp * ygp * det * wt;
      s_lhs[3][3] += zgp * zgp * det * wt;
      break;
    case 10: /* 3D with quadratic elements */
#ifdef RRL_DEBUG
#ifdef DBG_2
      fprintf(stdout, "fill_lhs_lspatch using max_terms = %d (3D quadratic case)\n", max_terms);
#endif
#endif
      xgp = xgp_loc[i];
      ygp = ygp_loc[i];
      zgp = zgp_loc[i];
      det = det_gp_loc[i];
      wt = wt_gp_loc[i];
      s_lhs[0][0] += 1.0 * det * wt;
      s_lhs[0][1] += xgp * det * wt;
      s_lhs[0][2] += ygp * det * wt;
      s_lhs[0][3] += zgp * det * wt;
      s_lhs[0][4] += xgp * xgp * det * wt;
      s_lhs[0][5] += xgp * ygp * det * wt;
      s_lhs[0][6] += xgp * zgp * det * wt;
      s_lhs[0][7] += ygp * ygp * det * wt;
      s_lhs[0][8] += ygp * zgp * det * wt;
      s_lhs[0][9] += zgp * zgp * det * wt;

      s_lhs[1][0] += xgp * det * wt;
      s_lhs[1][1] += xgp * xgp * det * wt;
      s_lhs[1][2] += xgp * ygp * det * wt;
      s_lhs[1][3] += xgp * zgp * det * wt;
      s_lhs[1][4] += xgp * xgp * xgp * det * wt;
      s_lhs[1][5] += xgp * xgp * ygp * det * wt;
      s_lhs[1][6] += xgp * xgp * zgp * det * wt;
      s_lhs[1][7] += xgp * ygp * ygp * det * wt;
      s_lhs[1][8] += xgp * ygp * zgp * det * wt;
      s_lhs[1][9] += xgp * zgp * zgp * det * wt;

      s_lhs[2][0] += ygp * det * wt;
      s_lhs[2][1] += ygp * xgp * det * wt;
      s_lhs[2][2] += ygp * ygp * det * wt;
      s_lhs[2][3] += ygp * zgp * det * wt;
      s_lhs[2][4] += ygp * xgp * xgp * det * wt;
      s_lhs[2][5] += ygp * ygp * xgp * det * wt;
      s_lhs[2][6] += ygp * xgp * zgp * det * wt;
      s_lhs[2][7] += ygp * ygp * ygp * det * wt;
      s_lhs[2][8] += ygp * ygp * zgp * det * wt;
      s_lhs[2][9] += ygp * zgp * zgp * det * wt;

      s_lhs[3][0] += zgp * det * wt;
      s_lhs[3][1] += zgp * xgp * det * wt;
      s_lhs[3][2] += zgp * ygp * det * wt;
      s_lhs[3][3] += zgp * zgp * det * wt;
      s_lhs[3][4] += zgp * xgp * xgp * det * wt;
      s_lhs[3][5] += zgp * xgp * ygp * det * wt;
      s_lhs[3][6] += zgp * zgp * xgp * det * wt;
      s_lhs[3][7] += zgp * ygp * ygp * det * wt;
      s_lhs[3][8] += zgp * zgp * ygp * det * wt;
      s_lhs[3][9] += zgp * zgp * zgp * det * wt;

      s_lhs[4][0] += xgp * xgp * det * wt;
      s_lhs[4][1] += xgp * xgp * xgp * det * wt;
      s_lhs[4][2] += xgp * xgp * ygp * det * wt;
      s_lhs[4][3] += xgp * xgp * zgp * det * wt;
      s_lhs[4][4] += xgp * xgp * xgp * xgp * det * wt;
      s_lhs[4][5] += xgp * xgp * xgp * ygp * det * wt;
      s_lhs[4][6] += xgp * xgp * xgp * zgp * det * wt;
      s_lhs[4][7] += xgp * xgp * ygp * ygp * det * wt;
      s_lhs[4][8] += xgp * xgp * ygp * zgp * det * wt;
      s_lhs[4][9] += xgp * xgp * zgp * zgp * det * wt;

      s_lhs[5][0] += xgp * ygp * det * wt;
      s_lhs[5][1] += xgp * xgp * ygp * det * wt;
      s_lhs[5][2] += xgp * ygp * ygp * det * wt;
      s_lhs[5][3] += xgp * ygp * zgp * det * wt;
      s_lhs[5][4] += xgp * xgp * xgp * ygp * det * wt;
      s_lhs[5][5] += xgp * xgp * ygp * ygp * det * wt;
      s_lhs[5][6] += xgp * xgp * ygp * zgp * det * wt;
      s_lhs[5][7] += xgp * ygp * ygp * ygp * det * wt;
      s_lhs[5][8] += xgp * ygp * ygp * zgp * det * wt;
      s_lhs[5][9] += xgp * ygp * zgp * zgp * det * wt;

      s_lhs[6][0] += xgp * zgp * det * wt;
      s_lhs[6][1] += xgp * xgp * zgp * det * wt;
      s_lhs[6][2] += xgp * ygp * zgp * det * wt;
      s_lhs[6][3] += xgp * zgp * zgp * det * wt;
      s_lhs[6][4] += xgp * xgp * xgp * zgp * det * wt;
      s_lhs[6][5] += xgp * xgp * ygp * zgp * det * wt;
      s_lhs[6][6] += xgp * xgp * zgp * zgp * det * wt;
      s_lhs[6][7] += xgp * ygp * ygp * zgp * det * wt;
      s_lhs[6][8] += xgp * ygp * zgp * zgp * det * wt;
      s_lhs[6][9] += xgp * zgp * zgp * zgp * det * wt;

      s_lhs[7][0] += ygp * ygp * det * wt;
      s_lhs[7][1] += ygp * ygp * xgp * det * wt;
      s_lhs[7][2] += ygp * ygp * ygp * det * wt;
      s_lhs[7][3] += ygp * ygp * zgp * det * wt;
      s_lhs[7][4] += ygp * ygp * xgp * xgp * det * wt;
      s_lhs[7][5] += ygp * ygp * ygp * xgp * det * wt;
      s_lhs[7][6] += ygp * ygp * xgp * zgp * det * wt;
      s_lhs[7][7] += ygp * ygp * ygp * ygp * det * wt;
      s_lhs[7][8] += ygp * ygp * ygp * zgp * det * wt;
      s_lhs[7][9] += ygp * ygp * zgp * zgp * det * wt;

      s_lhs[8][0] += ygp * zgp * det * wt;
      s_lhs[8][1] += ygp * xgp * zgp * det * wt;
      s_lhs[8][2] += ygp * ygp * zgp * det * wt;
      s_lhs[8][3] += ygp * zgp * zgp * det * wt;
      s_lhs[8][4] += ygp * xgp * xgp * zgp * det * wt;
      s_lhs[8][5] += ygp * ygp * xgp * zgp * det * wt;
      s_lhs[8][6] += ygp * xgp * zgp * zgp * det * wt;
      s_lhs[8][7] += ygp * ygp * ygp * zgp * det * wt;
      s_lhs[8][8] += ygp * ygp * zgp * zgp * det * wt;
      s_lhs[8][9] += ygp * zgp * zgp * zgp * det * wt;

      s_lhs[9][0] += zgp * zgp * det * wt;
      s_lhs[9][1] += zgp * zgp * xgp * det * wt;
      s_lhs[9][2] += zgp * zgp * ygp * det * wt;
      s_lhs[9][3] += zgp * zgp * zgp * det * wt;
      s_lhs[9][4] += zgp * zgp * xgp * xgp * det * wt;
      s_lhs[9][5] += zgp * zgp * xgp * ygp * det * wt;
      s_lhs[9][6] += zgp * zgp * zgp * xgp * det * wt;
      s_lhs[9][7] += zgp * zgp * ygp * ygp * det * wt;
      s_lhs[9][8] += zgp * zgp * zgp * ygp * det * wt;
      s_lhs[9][9] += zgp * zgp * zgp * zgp * det * wt;

      break;
    default:
      GOMA_EH(GOMA_ERROR, "Unsupported size in building LHS of least squares patch for error");
      break;
    }
  }
  return (status);
}

/*****************************************************************************/
/*   Routine calc_stream_fcn                                                 */
/*****************************************************************************/

static int calc_stream_fcn(double x[], /* soln vector */
                           double del_stream_fcn[4],
                           double vel[MAX_PDIM][MDE]) /* array for local
                                                       * nodal velocity
                                                       * values which must
                                                       * be divergence free*/

/******************************************************************************
  Function which calculates the stream function variation around an individual
  element

  Author:          P. R. Schunk (1511)
  Date:            9 September 1992
  Revised          9 January 1995 (Rich Cairncross to adapt for other vector fields)

******************************************************************************/
{
  double yy[9], xx[9];
  double fact1[3][3], fact2[2];
  int ngauss;
  int DeformingMesh;
  double gspt[3], gswt[3];
  int I, i, j, k, i1, i2, i3;
  int ileft, iright, Ileft, Iright;
  double xdiff, ydiff;
  double wx1, wx2, wx3, wy1, wy2, wy3;
  double dsx[3], dsy[3];
  double s, wt, r, aa, b, c;
  int status = 0;

  /* Initialize factors for analytical integration */

  fact1[0][0] = -0.5;
  fact1[1][0] = 0.666667;
  fact1[2][0] = -0.166667;
  fact1[0][1] = -0.666667;
  fact1[1][1] = 0.0;
  fact1[2][1] = 0.666667;
  fact1[0][2] = 0.166667;
  fact1[1][2] = -0.666667;
  fact1[2][2] = 0.5;

  fact2[0] = 0.333333;
  fact2[1] = 1.333333;

  init_vec_value(xx, 0.0, 9);
  init_vec_value(yy, 0.0, 9);

  ngauss = 3;
  gspt[0] = -0.7745966692;
  gspt[1] = 0.;
  gspt[2] = 0.7745966692;
  gswt[0] = 0.55555555556;
  gswt[1] = 0.8888888889;
  gswt[2] = 0.555555556;

  DeformingMesh = pd->e[pg->imtrx][R_MESH1];

  /*
   *  need to adapt this for subparametric mapping, because all
   *  d[][] are not defined
   */
  if (DeformingMesh) {
    for (i = 0; i < ei[pg->imtrx]->num_local_nodes; i++) {
      I = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + i];
      if (Index_Solution(I, MESH_DISPLACEMENT1, 0, 0, -2, pg->imtrx) != -1) {
        xx[i] = Coor[0][I] + *esp->d[0][ei[pg->imtrx]->ln_to_dof[R_MESH1][i]];
        yy[i] = Coor[1][I] + *esp->d[1][ei[pg->imtrx]->ln_to_dof[R_MESH2][i]];
      } else { /* elements are subparametric */
        if (i < 8) {
          /* make node lie halfway between adjacent nodes */
          ileft = i - 4;
          Ileft = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + ileft];
          iright = i - 3;
          if (iright == 4)
            iright = 0;
          Iright = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + iright];
          xx[i] = Coor[0][I] +
                  0.5 * (x[Index_Solution(Ileft, MESH_DISPLACEMENT1, 0, 0, -2, pg->imtrx)] +
                         x[Index_Solution(Iright, MESH_DISPLACEMENT1, 0, 0, -2, pg->imtrx)]);
          yy[i] = Coor[1][I] +
                  0.5 * (x[Index_Solution(Ileft, MESH_DISPLACEMENT2, 0, 0, -2, pg->imtrx)] +
                         x[Index_Solution(Iright, MESH_DISPLACEMENT2, 0, 0, -2, pg->imtrx)]);
        } else {
          /* put centroid in center */
          I = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + 8];
          xx[i] = Coor[0][I];
          yy[i] = Coor[1][I];
          for (ileft = 0; ileft < 4; ileft++) {
            Ileft = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + ileft];
            xx[i] += 0.25 * x[Index_Solution(Ileft, MESH_DISPLACEMENT1, 0, 0, -2, pg->imtrx)];
            yy[i] += 0.25 * x[Index_Solution(Ileft, MESH_DISPLACEMENT2, 0, 0, -2, pg->imtrx)];
          }
        }
      }
    }
  } else {
    for (i = 0; i < ei[pg->imtrx]->num_local_nodes; i++) {
      I = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + i];
      xx[i] = Coor[0][I];
      yy[i] = Coor[1][I];
    }
  }

  if (pd->CoordinateSystem == CARTESIAN || pd->CoordinateSystem == CARTESIAN_2pt5D) {
    if (pd_glob[0]->i[pg->imtrx][R_MOMENTUM1] == I_Q1) {
      GOMA_WH(-1, "Stream function with Q1 mapping may not be accurate ");

      for (i = 0; i < ei[pg->imtrx]->num_sides; i++) {
        i1 = i;
        i2 = i + 1;
        i3 = i + ei[pg->imtrx]->num_sides;
        if (i == ei[pg->imtrx]->num_sides - 1)
          i2 = 0;
        xdiff = (xx[i2] - xx[i1]) / 2.;
        ydiff = (yy[i2] - yy[i1]) / 2.;
        del_stream_fcn[i] =
            (fact2[0] * (vel[0][i1]) + fact2[1] * (vel[0][i3]) + fact2[0] * (vel[0][i2])) * ydiff -
            (fact2[0] * (vel[1][i1]) + fact2[1] * (vel[1][i3]) + fact2[0] * (vel[1][i2])) * xdiff;
      }
    } else {
      for (i = 0; i < ei[pg->imtrx]->num_sides; i++) {
        i1 = i;
        i2 = i + 1;
        i3 = i + ei[pg->imtrx]->num_sides;
        if (i == ei[pg->imtrx]->num_sides - 1)
          i2 = 0;
        wx1 = xx[i1] * fact1[0][0] + xx[i3] * fact1[1][0] + xx[i2] * fact1[2][0];
        wx2 = xx[i1] * fact1[0][1] + xx[i3] * fact1[1][1] + xx[i2] * fact1[2][1];
        wx3 = xx[i1] * fact1[0][2] + xx[i3] * fact1[1][2] + xx[i2] * fact1[2][2];
        wy1 = yy[i1] * fact1[0][0] + yy[i3] * fact1[1][0] + yy[i2] * fact1[2][0];
        wy2 = yy[i1] * fact1[0][1] + yy[i3] * fact1[1][1] + yy[i2] * fact1[2][1];
        wy3 = yy[i1] * fact1[0][2] + yy[i3] * fact1[1][2] + yy[i2] * fact1[2][2];

        del_stream_fcn[i] = wy1 * (vel[0][i1]) + wy2 * (vel[0][i3]) + wy3 * (vel[0][i2]) -
                            wx1 * (vel[1][i1]) - wx2 * (vel[1][i3]) - wx3 * (vel[1][i2]);
      }
    }

  } else if (pd->CoordinateSystem == CYLINDRICAL || pd->CoordinateSystem == SWIRLING) {

    if (pd_glob[0]->i[pg->imtrx][R_MOMENTUM1] == I_Q1) {
      GOMA_WH(-1, "Stream function with Q1 may not be accurate ");
      for (i = 0; i < ei[pg->imtrx]->num_sides; i++) {
        i1 = i;
        i2 = i + 1;
        i3 = i + ei[pg->imtrx]->num_sides;
        if (i == ei[pg->imtrx]->num_sides - 1)
          i2 = 0;
        xdiff = (xx[i2] - xx[i1]) / 2.;
        ydiff = (yy[i2] - yy[i1]) / 2.;
        del_stream_fcn[i] =
            (fact2[0] * (vel[0][i1]) + fact2[1] * (vel[0][i3]) + fact2[0] * (vel[0][i2])) * ydiff -
            (fact2[0] * (vel[1][i1]) + fact2[1] * (vel[1][i3]) + fact2[0] * (vel[1][i2])) * xdiff;
      }
    } else {
      for (i = 0; i < ei[pg->imtrx]->num_sides; i++) {
        i1 = i;
        i2 = i + 1;
        i3 = i + ei[pg->imtrx]->num_sides;
        if (i == ei[pg->imtrx]->num_sides - 1)
          i2 = 0;
        for (j = 0; j < 3; j++) {
          dsx[j] = 0.;
          dsy[j] = 0.;
          for (k = 0; k < ngauss; k++) {
            s = gspt[k];
            wt = gswt[k];
            r = yy[i1] * (s * s - s) * 0.5 + yy[i3] * (1. - s * s) + yy[i2] * (s * s + s) * 0.5;
            aa = (s * s - s) * 0.25 * (j - 1) * (j - 2) - (1. - s * s) * (j) * (j - 2) +
                 (s * s + s) * 0.25 * (j) * (j - 1);
            b = (2. * s - 1) * xx[i1] * 0.5 + (-2. * s) * xx[i3] + (2. * s + 1.) * xx[i2] * 0.5;
            c = (2. * s - 1.) * yy[i1] * 0.5 + (-2. * s) * yy[i3] + (2. * s + 1.) * yy[i2] * 0.5;
            dsx[j] = dsx[j] + wt * r * aa * b;
            dsy[j] = dsy[j] + wt * r * aa * c;
          }
        }
        del_stream_fcn[i] = dsx[0] * (vel[1][i1]) + dsx[1] * (vel[1][i3]) + dsx[2] * (vel[1][i2]) -
                            dsy[0] * (vel[0][i1]) - dsy[1] * (vel[0][i3]) - dsy[2] * (vel[0][i2]);
      }
    }

  } else {
    GOMA_EH(GOMA_ERROR, "Stream function routine called with incompatible coord system");
  }
  return (status);
} /* End routine stream function */
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

static int midsid(double stream_fcn_vect[], /* Soln vector for current proc */
                  Exo_DB *exo)

/******************************************************************************
  Function which fills interpolates the values of the stream function from
  the vertices to the midpoint nodes

  Author:          P. R. Schunk (1511)
  Date:            9 September 1992

******************************************************************************/

{
  /* local variables */
  int i, ii, iii, I;
  int ielem;
  double stream_fcn[9];

  for (ielem = 0; ielem < Num_Internal_Elems; ielem++) {

    ei[pg->imtrx]->iconnect_ptr = exo->elem_ptr[ielem];

    ei[pg->imtrx]->ielem_shape = type2shape(ei[pg->imtrx]->ielem_type);
    ei[pg->imtrx]->num_sides = shape2sides(ei[pg->imtrx]->ielem_shape);

    for (i = 0; i < ei[pg->imtrx]->num_sides; i++) {
      I = exo->node_list[ei[pg->imtrx]->iconnect_ptr + i];
      stream_fcn[i] = stream_fcn_vect[I];
    }
    /*
     *     QUADRILATERAL ELEMENTS (9 NODE) ONLY FOR NOW
     */
    stream_fcn[4] = 0.5 * (stream_fcn[0] + stream_fcn[1]);
    stream_fcn[5] = 0.5 * (stream_fcn[1] + stream_fcn[2]);
    stream_fcn[6] = 0.5 * (stream_fcn[2] + stream_fcn[3]);
    stream_fcn[7] = 0.5 * (stream_fcn[3] + stream_fcn[0]);
    stream_fcn[8] = 0.25 * (stream_fcn[0] + stream_fcn[1] + stream_fcn[2] + stream_fcn[3]);

    ii = ei[pg->imtrx]->num_sides;
    iii = ei[pg->imtrx]->num_sides + ei[pg->imtrx]->num_sides + 1;
    for (i = ii; i < iii; i++) {
      I = exo->node_list[ei[pg->imtrx]->iconnect_ptr + i];
      stream_fcn_vect[I] = stream_fcn[i];
    }
  }

  return (0);
}
/*________________________________________________________________________*/

static int correct_stream_fcn(int *kount, /* a counter for element connectivity ??     */
                              int iel,    /* current element number                    */
                              double del_stream_fcn[4], /* elemental side increments to *
                                                         * stream function calculated by*
                                                         * calc_stream_fcn()            */
                              double stream_fcn_vect[],
                              int listnd[]) /* count number of times node is accessed */

/******************************************************************************
  Function which corrects the values of the stream function

  Author:          P. R. Schunk (1511)
  Date:            9 September 1992
  Adapted:         R. A. Cairncross 9 January 1995

******************************************************************************/
{
  int nsideq[7];
  int i, ii, iii, iiii, I, J, nprvel, nstart;
  int nel, nnel, kel, j;
  int status = 0;
  double s, stream_fcn[9];

  if (ei[pg->imtrx]->num_sides > 4) {
    GOMA_EH(GOMA_ERROR, "correct_stream_fcn not available in 3D");
    return -1;
  }

  /* initialize other vectors */
  for (i = 0; i < 4; i++)
    nsideq[i] = i;
  for (i = 4; i < 7; i++)
    nsideq[i] = i - 4;
  for (i = 0; i < 9; i++)
    stream_fcn[i] = 0;

  /*
   *  Now correct for mass balance by looping over all the sides
   *  of the element, accumulating the stream function variable.
   *  The loop should end up being equal. If it is not, then
   *  fix it up, so that it is.
   */
  s = 0.;
  for (i = 0; i < ei[pg->imtrx]->num_sides; i++) {
    s += del_stream_fcn[i];
  }
  s = s / ei[pg->imtrx]->num_sides;
  for (i = 0; i < ei[pg->imtrx]->num_sides; i++) {
    del_stream_fcn[i] -= s;
  }

  /* Calculate the stream function at the nodal points */
  if (iel == 0) {
    stream_fcn[0] = 0.;
    for (i = 0; i < ei[pg->imtrx]->num_sides - 1; i++) {
      stream_fcn[i + 1] = stream_fcn[i] + del_stream_fcn[i];
    }
    for (i = 0; i < ei[pg->imtrx]->num_sides; i++) {
      I = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + i];
      stream_fcn_vect[I] += stream_fcn[i];
      listnd[I] += 1;
    }
  }

  if (iel != 0) {

    /* locate connecting node with known value of stream_fcn */

    nprvel = iel - 1;
    for (i = 0; i < ei[pg->imtrx]->num_sides; i++) {
      I = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + i];
      nstart = i;
      for (nel = 0; nel < nprvel + 1; nel++) {
        kel = nprvel - nel;
        nnel = listel[kel] - 1;
        /*  nnel = kel; */

        /* restriction here is that all previous elements
           must be of the same type */

        for (j = 0; j < ei[pg->imtrx]->num_sides; j++) {
          J = Proc_Elem_Connect[Proc_Connect_Ptr[nnel] + j];
          if (J == I && listnd[I] > 0)
            goto found;
        }
      }
    }
    goto nextelm;
  found:
    /* Connecting node located, calculate stream function for current element */
    stream_fcn[nstart] = stream_fcn_vect[J] / listnd[J];
    for (ii = 0; ii < ei[pg->imtrx]->num_sides - 1; ii++) {
      iii = nsideq[nstart + ii];
      iiii = nsideq[nstart + ii + 1];
      stream_fcn[iiii] = stream_fcn[iii] + del_stream_fcn[iii];
    }

    /* add to global vector */

    for (ii = 0; ii < ei[pg->imtrx]->num_sides; ii++) {
      I = Proc_Elem_Connect[ei[pg->imtrx]->iconnect_ptr + ii];
      stream_fcn_vect[I] += stream_fcn[ii];
      listnd[I] += 1;
    }
    goto nextelement;
  nextelm:
    (*kount)++;
  }
  if (*kount > 1) {
    GOMA_WH(-1, "calc_stream_fcn: POSSIBLE ERROR - Connectivity not sufficiently continuous\n");
  }
nextelement:
  *kount = *kount;

  return (status);
} /* end of correct_stream_fcn */
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
/*
 * rd_post_process_specs -- read post processing specification section of input file
 *
 * Comments:	This code was lifted out of the ever-growing read_input_file()
 *           	routine above. This makes things more modular, with division
 *           	into "sections".
 *
 *		Someday we need to comb through all these placeholder options
 *		to see if they are useful.
 *
 *
 */

void rd_post_process_specs(FILE *ifp, char *input) {
  static const char yo[] = "rd_post_process_specs";

  int iread; /* status flag from look_for_optional */
  char line[81];
  char data_line_buffer[MAX_CHAR_IN_INPUT];
  char *arguments[MAX_NUMBER_PARAMS];
  char first_string[MAX_CHAR_IN_INPUT];
  char filename[MAX_FNL];
  char second_string[MAX_CHAR_IN_INPUT];
  char ts[MAX_CHAR_IN_INPUT];
  int i, k, nargs, sz, sens_vec_ct, j;
  long save_position;
  int pbits[5];

  /*
   * Variables used to read one DATA line at a time...
   */
  int num_read_items;
  char variable_name[MAX_CHAR_IN_INPUT];
  int node_set_id;
  int species_id;
  char file_name[MAX_CHAR_IN_INPUT];
  char optional_format[MAX_CHAR_IN_INPUT];

  char echo_string[MAX_CHAR_ECHO_INPUT] = "\0";
  char echo_file[MAX_CHAR_IN_INPUT] = "\0";
  int elemBlock_id = -1;

  strcpy(echo_file, Echo_Input_File);
  /*
   * Initialize count of post proc vars to be exported
   */
  Num_Export_XP = 0;

  look_for(ifp, "Post Processing Specifications", input, '=');

  ECHO("\nPost Processing Specifications =\n", echo_file);

  iread = look_for_post_proc(ifp, "Stream Function", &STREAM);
  iread = look_for_post_proc(ifp, "Streamwise normal stress", &STREAM_NORMAL_STRESS);
  iread = look_for_post_proc(ifp, "Streamwise shear stress", &STREAM_SHEAR_STRESS);
  iread = look_for_post_proc(ifp, "Streamwise Stress Difference", &STREAM_TENSION);
  iread = look_for_post_proc(ifp, "Mean shear rate", &MEAN_SHEAR);
  iread = look_for_post_proc(ifp, "Pressure contours", &PRESSURE_CONT);
  iread = look_for_post_proc(ifp, "Shell div_s_v contours", &SH_DIV_S_V_CONT);
  iread = look_for_post_proc(ifp, "Shell curv contours", &SH_CURV_CONT);
  iread = look_for_post_proc(ifp, "Fill contours", &FILL_CONT);
  iread = look_for_post_proc(ifp, "Concentration contours", &CONC_CONT);
  iread = look_for_post_proc(ifp, "Stress contours", &STRESS_CONT);
  iread = look_for_post_proc(ifp, "EM contours", &EM_CONTOURS);
  iread = look_for_post_proc(ifp, "Total EM", &TOTAL_EM_CONTOURS);
  iread = look_for_post_proc(ifp, "Scattered EM", &SCATTERED_EM_CONTOURS);
  iread = look_for_post_proc(ifp, "First Invariant of Strain", &FIRST_INVAR_STRAIN);
  iread = look_for_post_proc(ifp, "Second Invariant of Strain", &SEC_INVAR_STRAIN);
  iread = look_for_post_proc(ifp, "Third Invariant of Strain", &THIRD_INVAR_STRAIN);
  iread = look_for_post_proc(ifp, "Velocity Divergence", &DIV_VELOCITY);
  iread = look_for_post_proc(ifp, "Particle Velocity Divergence", &DIV_PVELOCITY);
  iread = look_for_post_proc(ifp, "Total Velocity Divergence", &DIV_TOTAL);
  iread = look_for_post_proc(ifp, "Electric Field", &ELECTRIC_FIELD);
  iread = look_for_post_proc(ifp, "Electric Field Magnitude", &ELECTRIC_FIELD_MAG);
  iread = look_for_post_proc(ifp, "Viscosity", &PP_Viscosity);
  iread = look_for_post_proc(ifp, "FlowingLiquid Viscosity", &PP_FlowingLiquid_Viscosity);
  iread = look_for_post_proc(ifp, "Volume Fraction of Gas Phase", &PP_VolumeFractionGas);
  iread = look_for_post_proc(ifp, "Density", &DENSITY);
  iread = look_for_post_proc(ifp, "Density Time Derivative", &RHO_DOT);
  iread = look_for_post_proc(ifp, "Polymer Viscosity", &DENSITY);
  iread = look_for_post_proc(ifp, "Polymer Viscosity", &POLYMER_VISCOSITY);
  iread = look_for_post_proc(ifp, "Polymer Time Constant", &POLYMER_TIME_CONST);
  iread = look_for_post_proc(ifp, "Mobility Parameter", &MOBILITY_PARAMETER);
  iread = look_for_post_proc(ifp, "PTT Xi parameter", &PTT_XI);
  iread = look_for_post_proc(ifp, "PTT Epsilon parameter", &PTT_EPSILON);
  iread = look_for_post_proc(ifp, "Navier Stokes Residuals", &NS_RESIDUALS);
  iread = look_for_post_proc(ifp, "Moving Mesh Residuals", &MM_RESIDUALS);
  iread = look_for_post_proc(ifp, "Mass Diffusion Vectors", &DIFFUSION_VECTORS);
  iread = look_for_post_proc(ifp, "Dielectrophoretic Field", &DIELECTROPHORETIC_FIELD);
  iread = look_for_post_proc(ifp, "Dielectrophoretic Field Norm", &DIELECTROPHORETIC_FIELD_NORM);
  iread = look_for_post_proc(ifp, "Enormsq Field", &ENORMSQ_FIELD);
  iread = look_for_post_proc(ifp, "Enormsq Field Norm", &ENORMSQ_FIELD_NORM);
  iread = look_for_post_proc(ifp, "Diffusive Mass Flux Vectors", &DIFFUSION_VECTORS);
  iread = look_for_post_proc(ifp, "Mass Fluxlines", &FLUXLINES);
  iread = look_for_post_proc(ifp, "Energy Conduction Vectors", &CONDUCTION_VECTORS);
  iread = look_for_post_proc(ifp, "Energy Fluxlines", &ENERGY_FLUXLINES);
  iread = look_for_post_proc(ifp, "Time Derivatives", &TIME_DERIVATIVES);
  iread = look_for_post_proc(ifp, "Mesh Stress Tensor", &STRESS_TENSOR);
  iread = look_for_post_proc(ifp, "Real Solid Stress Tensor", &REAL_STRESS_TENSOR);
  iread = look_for_post_proc(ifp, "Mesh Strain Tensor", &STRAIN_TENSOR);
  iread = look_for_post_proc(ifp, "Viscoplastic Def_Grad Tensor", &EVP_DEF_GRAD_TENSOR);
  iread = look_for_post_proc(ifp, "Lagrangian Convection", &LAGRANGE_CONVECTION);
  iread = look_for_post_proc(ifp, "Normal and Tangent Vectors", &SURFACE_VECTORS);
  iread = look_for_post_proc(ifp, "Shell Normal Vectors", &SHELL_NORMALS);
  iread = look_for_post_proc(ifp, "Error ZZ velocity", &ERROR_ZZ_VEL);
  iread = look_for_post_proc(ifp, "Error ZZ heat flux", &ERROR_ZZ_Q);
  iread = look_for_post_proc(ifp, "Error ZZ pressure", &ERROR_ZZ_P);
  iread = look_for_post_proc(ifp, "Map Log-Conf Stress", &CONF_MAP);
  iread = look_for_post_proc(ifp, "Map Conf Stress", &CONF_MAP);
  iread = look_for_post_proc(ifp, "Velocity Magnitude", &VELO_SPEED);
  iread = look_for_post_proc(ifp, "Giesekus Criterion", &GIES_CRIT);
  iread = look_for_post_proc(ifp, "Particle stress flux", &J_FLUX);
  iread = look_for_post_proc(ifp, "Eigenvalues", &EIG);
  iread = look_for_post_proc(ifp, "Eigenvector1", &EIG1);
  iread = look_for_post_proc(ifp, "Eigenvector2", &EIG2);
  iread = look_for_post_proc(ifp, "Eigenvector3", &EIG3);
  iread = look_for_post_proc(ifp, "Shear gradient", &GRAD_SH);
  iread = look_for_post_proc(ifp, "Concentration gradient", &GRAD_Y);
  iread = look_for_post_proc(ifp, "Vorticity Vector", &CURL_V);
  iread = look_for_post_proc(ifp, "Helicity Value", &HELICITY);
  iread = look_for_post_proc(ifp, "Moment Sources", &MOMENT_SOURCES);
  iread = look_for_post_proc(ifp, "YZbeta Species", &YZBETA);
  iread = look_for_post_proc(ifp, "Heaviside", &HEAVISIDE);
  iread = look_for_post_proc(ifp, "Lamb Vector", &LAMB_VECTOR);
  iread = look_for_post_proc(ifp, "Q Function", &Q_FCN);
  iread = look_for_post_proc(ifp, "Poynting Vectors", &POYNTING_VECTORS);
  iread = look_for_post_proc(ifp, "PSPG Post", &PSPG_PP);
  iread = look_for_post_proc(ifp, "Saramito Yield Coeff", &SARAMITO_YIELD);
  iread = look_for_post_proc(ifp, "VE Stress Norm", &STRESS_NORM);
  iread = look_for_post_proc(ifp, "Species Sources", &SPECIES_SOURCES);
  iread = look_for_post_proc(ifp, "Viscous Stress", &VISCOUS_STRESS);
  iread = look_for_post_proc(ifp, "Viscous Stress Norm", &VISCOUS_STRESS_NORM);
  iread = look_for_post_proc(ifp, "Viscous Von Mises Stress", &VISCOUS_VON_MISES_STRESS);
  iread = look_for_post_proc(ifp, "Orientation Vectors", &ORIENTATION_VECTORS);
  iread = look_for_post_proc(ifp, "First StrainRate Invariant", &FIRST_STRAINRATE_INVAR);
  iread = look_for_post_proc(ifp, "Second StrainRate Invariant", &SEC_STRAINRATE_INVAR);
  iread = look_for_post_proc(ifp, "Third StrainRate Invariant", &THIRD_STRAINRATE_INVAR);
  iread = look_for_post_proc(ifp, "User-Defined Post Processing", &USER_POST);

  /*
   * Initialize for surety before communication to other processors.
   */

  len_u_post_proc = 0;

  if (USER_POST == 1) /* may need parameters for this option */
  {
    if (fgets(line, 81, ifp) != NULL) {
      double dummy[1] = {-1.0};
      strip(line);
      len_u_post_proc = count_parameters(line);
      if (len_u_post_proc > 0) {
        u_post_proc = (dbl *)array_alloc(1, len_u_post_proc, sizeof(dbl));

        /* parse parameters into little strings */
        tokenize_by_whsp(line, arguments, MAX_NUMBER_PARAMS);
        for (i = 0; i < len_u_post_proc; i++) {
          u_post_proc[i] = atof(arguments[i]);
        }
      } else {
        u_post_proc = dummy;
      }
    }
  }

  iread = look_for_post_proc(ifp, "Acoustic Pressure Magnitude", &ACOUSTIC_PRESSURE);
  iread = look_for_post_proc(ifp, "Acoustic Phase Angle", &ACOUSTIC_PHASE_ANGLE);
  iread = look_for_post_proc(ifp, "Acoustic Energy Density", &ACOUSTIC_ENERGY_DENSITY);
  iread = look_for_post_proc(ifp, "Light Intensity", &LIGHT_INTENSITY);
  iread = look_for_post_proc(ifp, "Principal Stress", &PRINCIPAL_STRESS);
  iread = look_for_post_proc(ifp, "Principal Real Stress", &PRINCIPAL_REAL_STRESS);
  iread = look_for_post_proc(ifp, "Lubrication Height", &LUB_HEIGHT);
  iread = look_for_post_proc(ifp, "Lubrication Height 2", &LUB_HEIGHT_2);
  iread = look_for_post_proc(ifp, "Lubrication Upper Velocity", &LUB_VELO_UPPER);
  iread = look_for_post_proc(ifp, "Lubrication Lower Velocity", &LUB_VELO_LOWER);
  iread = look_for_post_proc(ifp, "Lubrication Velocity Field", &LUB_VELO_FIELD);
  iread = look_for_post_proc(ifp, "Lubrication Velocity Field 2", &LUB_VELO_FIELD_2);
  iread = look_for_post_proc(ifp, "Lubrication Fluid Source", &LUB_FLUID_SOURCE);
  iread = look_for_post_proc(ifp, "Disjoining Pressure", &DISJ_PRESS);
  iread = look_for_post_proc(ifp, "Porous Shell Open Saturation", &SH_SAT_OPEN);
  iread = look_for_post_proc(ifp, "Porous Shell Open Saturation 2", &SH_SAT_OPEN_2);
  iread = look_for_post_proc(ifp, "Porous Shell Capillary Pressure", &SH_CAP_PRES);
  iread = look_for_post_proc(ifp, "Porous Shell Inter Layer Flux", &SH_PORE_FLUX);
  iread = look_for_post_proc(ifp, "Shell Stress Tensor", &SH_STRESS_TENSOR);
  iread = look_for_post_proc(ifp, "Shell Tangents", &SH_TANG);
  iread = look_for_post_proc(ifp, "Lame MU", &PP_LAME_MU);
  iread = look_for_post_proc(ifp, "Lame LAMBDA", &PP_LAME_LAMBDA);
  iread = look_for_post_proc(ifp, "Von Mises Strain", &VON_MISES_STRAIN);
  iread = look_for_post_proc(ifp, "Von Mises Stress", &VON_MISES_STRESS);
  iread = look_for_post_proc(ifp, "Untracked Species", &UNTRACKED_SPEC);
  iread = look_for_post_proc(ifp, "Porous Saturation", &POROUS_SATURATION);
  iread = look_for_post_proc(ifp, "Total density of solvents in porous media",
                             &POROUS_RHO_TOTAL_SOLVENTS);
  iread = look_for_post_proc(ifp, "Density of solvents in gas phase in porous media",
                             &POROUS_RHO_GAS_SOLVENTS);
  iread = look_for_post_proc(ifp, "Density of liquid phase in porous media", &POROUS_RHO_LPHASE);
  iread = look_for_post_proc(ifp, "Gas phase Darcy velocity in porous media", &DARCY_VELOCITY_GAS);
  iread =
      look_for_post_proc(ifp, "Liquid phase Darcy velocity in porous media", &DARCY_VELOCITY_LIQ);
  iread = look_for_post_proc(ifp, "Local liquid accumulation rate in porous media",
                             &POROUS_LIQUID_ACCUM_RATE);
  iread = look_for_post_proc(ifp, "Capillary pressure in porous media", &CAPILLARY_PRESSURE);
  iread = look_for_post_proc(ifp, "Grid Peclet Number in porous media", &POROUS_GRIDPECLET);
  iread = look_for_post_proc(ifp, "SUPG Velocity in porous media", &POROUS_SUPGVELOCITY);
  iread = look_for_post_proc(ifp, "Relative Liquid Permeability", &REL_LIQ_PERM);

  iread = look_for_post_proc(ifp, "TFMP_gas_velo", &TFMP_GAS_VELO);
  iread = look_for_post_proc(ifp, "TFMP_liq_velo", &TFMP_LIQ_VELO);
  iread = look_for_post_proc(ifp, "TFMP_inverse_Peclet", &TFMP_INV_PECLET);
  iread = look_for_post_proc(ifp, "TFMP_Krg", &TFMP_KRG);

  /* Report count of post-proc vars to be exported */
  /*
    fprintf(stderr, "Goma will export %d post-processing variables.\n", Num_Export_XP);
  */

  /*
   * Make theses additions in the same format as the above post processing variables
   *
   *  Post-processing Step 2: add a new call to look_for_post_proc in mm_post_proc
   *                          rd_post_proc_specs
   *                          to search input file for your new variable
   */

  /*
   *  Schedule post processing error-based element size calculation if needed
   */

  for (i = 0; i < 3; i++) {
    error_presence_key[i] = 0;
  }

  i = 0;
  nn_error_metrics = 0;
  if (look_for_optional(ifp, "Error ZZ velocity elem size", input, '=') == 1) {
    if (ERROR_ZZ_VEL == -1) {
      GOMA_EH(GOMA_ERROR, "'Error ZZ velocity elem size' card REQUIRES 'Error ZZ "
                          "velocity = yes' card - please add");
    } else {
      ERROR_ZZ_VEL_ELSIZE = 1;
      error_presence_key[i] = 1;
      nn_error_metrics++;
    }
  } else {
    ERROR_ZZ_VEL_ELSIZE = -1;
  }
  i++;
  if (look_for_optional(ifp, "Error ZZ heat flux elem size", input, '=') == 1) {
    if (ERROR_ZZ_Q == -1) {
      GOMA_EH(GOMA_ERROR,
              "'Error ZZ heat flux elem size' card REQUIRES 'Error ZZ heat flux = yes' card - "
              "please add");
    } else {
      ERROR_ZZ_Q_ELSIZE = 1;
      error_presence_key[i] = 1;
      nn_error_metrics++;
    }
  } else {
    ERROR_ZZ_Q_ELSIZE = -1;
  }
  i++;
  if (look_for_optional(ifp, "Error ZZ pressure elem size", input, '=') == 1) {
    if (ERROR_ZZ_P == -1) {
      GOMA_EH(GOMA_ERROR,
              "'Error ZZ pressure elem size' card REQUIRES 'Error ZZ pressure = yes' card - "
              "please add");
    } else {
      ERROR_ZZ_P_ELSIZE = 1;
      error_presence_key[i] = 1;
      nn_error_metrics++;
    }
  } else {
    ERROR_ZZ_P_ELSIZE = -1;
  }

  /*
   *  Allocate memory to hold the error information
   */

  if (nn_error_metrics > 0) {
    pp_error_data = (struct Post_Processing_Error *)smalloc(nn_error_metrics *
                                                            sizeof(struct Post_Processing_Error));

    i = 0;
    if (error_presence_key[i]) {
      /* Need this line to line up ifp? at the right place I think */
      look_for_optional(ifp, "Error ZZ velocity elem size", input, '=');
      if (fscanf(ifp, "%lf %lf %lf %lf %lf %lf", &pp_error_data[i].error_params[0],
                 &pp_error_data[i].error_params[1], &pp_error_data[i].error_params[2],
                 &pp_error_data[i].error_params[3], &pp_error_data[i].error_params[4],
                 &pp_error_data[i].error_params[5]) != 6) {
        GOMA_EH(GOMA_ERROR,
                " Error reading 'Error ZZ velocity elem size' card values - expecting 6 floats");
      }
#ifdef RRL_DEBUG
#ifdef DBG_0
      fprintf(stdout,
              "Read params for ZZ vel elem size (%6.4f) (%6.4f) (%6.4f) (%6.4f) (%6.4f) (%6.4f)\n",
              pp_error_data[i].error_params[0], pp_error_data[i].error_params[1],
              pp_error_data[i].error_params[2], pp_error_data[i].error_params[3],
              pp_error_data[i].error_params[4], pp_error_data[i].error_params[5]);
#endif
#endif
      SPF(echo_string, "%s = %.4g %.4g %.4g %.4g %.4g %.4g", input,
          pp_error_data[i].error_params[0], pp_error_data[i].error_params[1],
          pp_error_data[i].error_params[2], pp_error_data[i].error_params[3],
          pp_error_data[i].error_params[4], pp_error_data[i].error_params[5]);

      ECHO(echo_string, echo_file);
    }
    i++;
    if (error_presence_key[i]) {
      look_for(ifp, "Error ZZ heat flux elem size", input, '=');
      if (fscanf(ifp, "%lf %lf %lf %lf %lf %lf", &pp_error_data[i].error_params[0],
                 &pp_error_data[i].error_params[1], &pp_error_data[i].error_params[2],
                 &pp_error_data[i].error_params[3], &pp_error_data[i].error_params[4],
                 &pp_error_data[i].error_params[5]) != 6) {
        GOMA_EH(GOMA_ERROR,
                " Error reading 'Error ZZ heat flux elem size' card values - expecting 6 floats");
      }
#ifdef RRL_DEBUG
#ifdef DBG_1
      fprintf(stdout,
              "Read params for ZZ heat flux elem size (%6.4f) (%6.4f) (%6.4f) (%6.4f) (%6.4f) "
              "(%6.4f)\n",
              pp_error_data[i].error_params[0], pp_error_data[i].error_params[1],
              pp_error_data[i].error_params[2], pp_error_data[i].error_params[3],
              pp_error_data[i].error_params[4], pp_error_data[i].error_params[5]);
#endif
#endif
      SPF(echo_string, "%s = %.4g %.4g %.4g %.4g %.4g %.4g", input,
          pp_error_data[i].error_params[0], pp_error_data[i].error_params[1],
          pp_error_data[i].error_params[2], pp_error_data[i].error_params[3],
          pp_error_data[i].error_params[4], pp_error_data[i].error_params[5]);

      ECHO(echo_string, echo_file);
    }
    i++;
    if (error_presence_key[i]) {
      look_for(ifp, "Error ZZ pressure elem size", input, '=');
      if (fscanf(ifp, "%lf %lf %lf %lf %lf %lf", &pp_error_data[i].error_params[0],
                 &pp_error_data[i].error_params[1], &pp_error_data[i].error_params[2],
                 &pp_error_data[i].error_params[3], &pp_error_data[i].error_params[4],
                 &pp_error_data[i].error_params[5]) != 6) {
        GOMA_EH(GOMA_ERROR,
                " Error reading 'Error ZZ pressure elem size' card values - expecting 6 floats");
      }
#ifdef RRL_DEBUG
#ifdef DBG_1
      fprintf(
          stdout,
          "Read params for ZZ pressure elem size (%6.4f) (%6.4f) (%6.4f) (%6.4f) (%6.4f) (%6.4f)\n",
          pp_error_data[i].error_params[0], pp_error_data[i].error_params[1],
          pp_error_data[i].error_params[2], pp_error_data[i].error_params[3],
          pp_error_data[i].error_params[4], pp_error_data[i].error_params[5]);
#endif
#endif
      SPF(echo_string, "%s = %.4g %.4g %.4g %.4g %.4g %.4g", input,
          pp_error_data[i].error_params[0], pp_error_data[i].error_params[1],
          pp_error_data[i].error_params[2], pp_error_data[i].error_params[3],
          pp_error_data[i].error_params[4], pp_error_data[i].error_params[5]);

      ECHO(echo_string, echo_file);
    }
  }

  /*
   *  SCHEDULE POST-PROCESSING FLUX CALCULATIONS, IF NEEDED
   *
   */

  iread = look_for_optional(ifp, "Post Processing Fluxes", input, '=');

  /* count number of post-processing flux calculation specifications */

  if (iread == 1) {
    nn_post_fluxes = count_list(ifp, "FLUX", input, '=', "END OF FLUX");
    ECHO("\nPost Processing Fluxes =\n", echo_file);
  } else {
    nn_post_fluxes = 0;
  }

  /*
   *  Allocate memory to hold the flux information
   */

  if (nn_post_fluxes > 0) {
    sz = sizeof(struct Post_Processing_Fluxes *);
    pp_fluxes = (struct Post_Processing_Fluxes **)array_alloc(1, nn_post_fluxes, sz);

    sz = sizeof(struct Post_Processing_Fluxes);

    for (i = 0; i < nn_post_fluxes; i++) {
      pp_fluxes[i] = (struct Post_Processing_Fluxes *)array_alloc(1, 1, sz);
    }
    /*Now load up information by reading cards */

    for (i = 0; i < nn_post_fluxes; i++) {
      look_for(ifp, "FLUX", input, '=');

      /* Read FLUX  type: i.e. FORCE_X, HEAT FLUX, etc. */

      if (fscanf(ifp, "%s", ts) != 1) {
        GOMA_EH(GOMA_ERROR, "error reading Post Processing Flux input variable");
      }

      pp_fluxes[i]->flux_type = -1;
      for (k = 0; k < Num_Flux_Names; k++) {
        if (!strcmp(ts, pp_flux_names[k].name)) {
          pp_fluxes[i]->flux_type = pp_flux_names[k].Index;
          strcpy(pp_fluxes[i]->flux_type_name, ts);
          break;
        }
      }
      if (pp_fluxes[i]->flux_type == -1) {
        GOMA_EH(GOMA_ERROR, "Invalid Flux name");
      }
      SPF(echo_string, "%s = %s", "FLUX", pp_fluxes[i]->flux_type_name);

      /* read ss id */

      if (fscanf(ifp, "%d", &pp_fluxes[i]->ss_id) != 1) {
        fprintf(stderr, "%s:\tError reading flux->ss_id\n", yo);
        exit(-1);
      }

      /* read block id */

      if (fscanf(ifp, "%d", &pp_fluxes[i]->blk_id) != 1) {
        fprintf(stderr, "%s:\tError reading flux->blk_id\n", yo);
        exit(-1);
      }

      SPF(endofstring(echo_string), " %d %d", pp_fluxes[i]->ss_id, pp_fluxes[i]->blk_id);

      /*	  pp_fluxes[i]->blk_id--;  */

      /* read species id*/

      if (fscanf(ifp, "%d", &pp_fluxes[i]->species_number) != 1) {
        fprintf(stderr, "%s:\tError reading flux->species_number\n", yo);
        exit(-1);
      }

      /* read file name */
      if (fscanf(ifp, "%s", ts) != 1) {
        GOMA_EH(GOMA_ERROR, "error reading Post Processing Flux filename");
      }
      strcpy(pp_fluxes[i]->flux_filenm, ts);

      SPF(endofstring(echo_string), " %d %s", pp_fluxes[i]->species_number,
          pp_fluxes[i]->flux_filenm);

      read_string(ifp, first_string, '\n');
      nargs = sscanf(first_string, "%s %d %d %d %d %d", second_string, &pbits[0], &pbits[1],
                     &pbits[2], &pbits[3], &pbits[4]);
      if (nargs > 0 && strcmp(second_string, "profile") == 0) {
        pp_fluxes[i]->profile_flag = 32;
        SPF(endofstring(echo_string), " %s", second_string);
        sz = 1;
        for (j = 1; j < nargs; j++) {
          if (pbits[j - 1])
            pp_fluxes[i]->profile_flag += sz;
          sz *= 2;
        }
        SPF(endofstring(echo_string), " %d", pp_fluxes[i]->profile_flag);
      } else {
        pp_fluxes[i]->profile_flag = FALSE;
      }
      ECHO(echo_string, echo_file);
    }
    ECHO("\nEND OF FLUX\n", echo_file);
  }

  /*
   *  SCHEDULE POST-PROCESSING FLUX SENSITIVITY CALCULATIONS, IF NEEDED
   *
   */

  /*   initialize sensitivity vector counter */
  sens_vec_ct = 0;

  iread = look_for_optional(ifp, "Post Processing Flux Sensitivities", input, '=');

  /* count number of post-processing flux calculation specifications */
  if (iread == 1) {
    nn_post_fluxes_sens = count_list(ifp, "FLUX_SENS", input, '=', "END OF FLUX_SENS");
    ECHO("\nPost Processing Flux Sensitivities =\n", echo_file);
  } else {
    nn_post_fluxes_sens = 0;
  }

  /*
   *  Allocate memory to hold the flux information
   */

  if (nn_post_fluxes_sens > 0) {
    sz = sizeof(struct Post_Processing_Fluxes_Sens *);
    pp_fluxes_sens = (struct Post_Processing_Fluxes_Sens **)array_alloc(1, nn_post_fluxes_sens, sz);

    sz = sizeof(struct Post_Processing_Fluxes_Sens);

    for (i = 0; i < nn_post_fluxes_sens; i++) {
      pp_fluxes_sens[i] = (struct Post_Processing_Fluxes_Sens *)array_alloc(1, 1, sz);
    }
    /*Now load up information by reading cards */

    for (i = 0; i < nn_post_fluxes_sens; i++) { /* nn_post_fluxes */
      look_for(ifp, "FLUX_SENS", input, '=');
      /* Read FLUX  type: i.e. FORCE_X, HEAT FLUX, etc. */
      if (fscanf(ifp, "%s", ts) != 1) {
        GOMA_EH(GOMA_ERROR, "error reading FLUX_SENS name");
      }

      pp_fluxes_sens[i]->flux_type = -1;
      for (k = 0; k < Num_Flux_Names; k++) {
        if (!strcmp(ts, pp_flux_names[k].name)) {
          pp_fluxes_sens[i]->flux_type = pp_flux_names[k].Index;
          strcpy(pp_fluxes_sens[i]->flux_type_name, ts);
          break;
        }
      }

      if (pp_fluxes_sens[i]->flux_type == -1) {
        GOMA_EH(GOMA_ERROR, "Invalid Flux name in Flux Sensitivity");
      }

      sprintf(echo_string, "%s = %s", "FLUX", pp_fluxes_sens[i]->flux_type_name);

      /* read ss id */

      if (fscanf(ifp, "%d", &pp_fluxes_sens[i]->ss_id) != 1) {
        fprintf(stderr, "%s:\tError reading flux_sens->ss_id\n", yo);
        exit(-1);
      }

      /* read block id */

      if (fscanf(ifp, "%d", &pp_fluxes_sens[i]->blk_id) != 1) {
        fprintf(stderr, "%s:\tError reading flux_sens->blk_id\n", yo);
        exit(-1);
      }

      /*	  pp_fluxes_sens[i]->blk_id--;*/

      /* read species id*/

      if (fscanf(ifp, "%d", &pp_fluxes_sens[i]->species_number) != 1) {
        fprintf(stderr, "%s:\tError reading flux_sens->species_number\n", yo);
        exit(-1);
      }
      sprintf(endofstring(echo_string), " %d %d %d", pp_fluxes_sens[i]->ss_id,
              pp_fluxes_sens[i]->blk_id, pp_fluxes_sens[i]->species_number);

      /*  read sensitivity variable type */

      if (fscanf(ifp, "%80s", input) != 1) {
        fprintf(stderr, "%s:\tError reading flux_sens->sens_type (BC, MT, AC, UM, AN, or UF)\n",
                yo);
        exit(-1);
      }
      if (!strcmp(input, "BC")) {
        pp_fluxes_sens[i]->sens_type = 1;
      } else if (!strcmp(input, "MT")) {
        pp_fluxes_sens[i]->sens_type = 2;
      } else if (!strcmp(input, "AC")) {
        pp_fluxes_sens[i]->sens_type = 3;
      } else if (!strcmp(input, "UM")) {
        pp_fluxes_sens[i]->sens_type = 4;
      } else if (!strcmp(input, "UF")) {
        pp_fluxes_sens[i]->sens_type = 5;
      } else if (!strcmp(input, "AN")) {
        pp_fluxes_sens[i]->sens_type = 6;
      } else {
        fprintf(stderr, "%s:\tImproper set_type for flux sensitivity - %s\n", yo, input);
        exit(-1);
      }

      /*  read BC id or material id */

      if (fscanf(ifp, "%d", &pp_fluxes_sens[i]->sens_id) != 1) {
        fprintf(stderr, "%s:\tError reading flux_sens->sens_id\n", yo);
        exit(-1);
      }

      sprintf(endofstring(echo_string), " %s %d", input, pp_fluxes_sens[i]->sens_id);

      if (pp_fluxes_sens[i]->sens_type == 2 || pp_fluxes_sens[i]->sens_type == 4)
        pp_fluxes_sens[i]->sens_id--;

      /*  read data float or material property number */

      if (fscanf(ifp, "%d", &pp_fluxes_sens[i]->sens_flt) != 1) {
        fprintf(stderr, "%s:\tError reading flux_sens->sens_flt\n", yo);
        exit(-1);
      }

      sprintf(endofstring(echo_string), " %4i", pp_fluxes_sens[i]->sens_flt);

      if (pp_fluxes_sens[i]->sens_type == 4) {
        if (fscanf(ifp, "%d", &pp_fluxes_sens[i]->sens_flt2) != 1) {
          fprintf(stderr, "%s:\tError reading flux_sens->sens_flt2\n", yo);
          exit(-1);
        }
        sprintf(endofstring(echo_string), " %4i", pp_fluxes_sens[i]->sens_flt2);
      } else {
        pp_fluxes_sens[i]->sens_flt2 = -1;
      }

      /*      determine sensitivity vector number  */
      /*   search for same previous sensitivity info */

      pp_fluxes_sens[i]->vector_id = -1;

      for (j = 0; j < i; j++) {
        if (pp_fluxes_sens[i]->sens_type == pp_fluxes_sens[j]->sens_type &&
            pp_fluxes_sens[i]->sens_id == pp_fluxes_sens[j]->sens_id &&
            pp_fluxes_sens[i]->sens_flt == pp_fluxes_sens[j]->sens_flt &&
            pp_fluxes_sens[i]->sens_flt2 == pp_fluxes_sens[j]->sens_flt2) {
          pp_fluxes_sens[i]->vector_id = pp_fluxes_sens[j]->vector_id;
        }
      }

      if (pp_fluxes_sens[i]->vector_id == -1) {
        pp_fluxes_sens[i]->vector_id = sens_vec_ct;

        if (Continuation == ALC_FIRST) {
          if (cont->upType == pp_fluxes_sens[i]->sens_type) {
            int id1 = -1, id2 = -1;
            int id3 = -1;

            switch (cont->upType) {
            case 1:
            case 3:
              id1 = cont->upBCID;
              id2 = cont->upDFID;
              break;
            case 2:
              id1 = cont->upMTID;
              id2 = cont->upMPID;
              break;
            case 4:
              id1 = cont->upMTID;
              id2 = cont->upMPID;
              id3 = cont->upMDID;
              break;
            case 5:
              GOMA_EH(GOMA_ERROR, "sensitivities to UF not done");
              break;
            case 6:
              GOMA_EH(GOMA_ERROR, "sensitivities to AN not done");
              break;
            }
            if (id1 == pp_fluxes_sens[i]->sens_id && id2 == pp_fluxes_sens[i]->sens_flt &&
                id3 == pp_fluxes_sens[i]->sens_flt2) {
              cont->sensvec_id = sens_vec_ct;
            }
          }
        }

        sens_vec_ct++;
      }

      /* read file name */
      if (fscanf(ifp, "%s", ts) != 1) {
        GOMA_EH(GOMA_ERROR, "error reading Post Processing Flux Sens filename");
      }
      strcpy(pp_fluxes_sens[i]->flux_filenm, ts);

      read_string(ifp, first_string, '\n');
      nargs = sscanf(first_string, "%s", second_string);
      if (nargs == 1 && strcmp(second_string, "profile") == 0) {
        pp_fluxes_sens[i]->profile_flag = TRUE;
        SPF(endofstring(echo_string), " %s", second_string);
      } else {
        pp_fluxes_sens[i]->profile_flag = FALSE;
      }
      ECHO(echo_string, echo_file);
    }
    ECHO("\nEND OF FLUX\n", echo_file);
  }

  /*
   *  SCHEDULE POST-PROCESSING DATA CALCULATIONS, IF NEEDED
   *
   */

  iread = look_for_optional(ifp, "Post Processing Data", input, '=');

  /* count number of Post Processing Data specifications */
  if (iread == 1) {
    nn_post_data = count_list(ifp, "DATA", input, '=', "END OF DATA");

    ECHO("\nPost Processing Data =\n", echo_file);

  } else {
    nn_post_data = 0;
  }

  /*
   *  Allocate memory to hold the data information
   */

  if (nn_post_data > 0) {
    sz = sizeof(struct Post_Processing_Data *);
    pp_data = (struct Post_Processing_Data **)array_alloc(1, nn_post_data, sz);

    sz = sizeof(struct Post_Processing_Data);

    for (i = 0; i < nn_post_data; i++) {
      pp_data[i] = (struct Post_Processing_Data *)array_alloc(1, 1, sz);
    }
    /*Now load up information by reading cards */

    /*
     * DATA = {vbl_name} {nodesetid} {elemblockid} {speciesnum} {file} {format}
     *
     *         vbl_name    = VELOCITY1, DMX, dmy, ..., theta
     *         nodesetid   = 2, 5, 2001, ...
     *         elemblockid = 1, 2, 12, ...
     *         speciesnum  = 1, 2, ...
     *         file        = myfile.d
     *         format      = [t|x|y|z]
     *         first_time  = FALSE
     */

    for (i = 0; i < nn_post_data; i++) {
      char *fgetsretval;
      look_for(ifp, "DATA", input, '=');

      save_position = ftell(ifp);
      fgetsretval = fgets(data_line_buffer, MAX_CHAR_IN_INPUT, ifp);
      if (fgetsretval == NULL) {
        GOMA_EH(GOMA_ERROR, "Error reading post processing line");
      }

      fseek(ifp, save_position, SEEK_SET);

      /*
       * Clean out variables.
       */

      strncpy(variable_name, "\0", 1);
      node_set_id = -1;
      species_id = -1;
      strncpy(file_name, "\0", 1);
      for (j = 0; j < MAX_CHAR_IN_INPUT; j++)
        optional_format[j] = '\0';

      num_read_items = sscanf(data_line_buffer, "%s %d %d %d %s %s", variable_name, &node_set_id,
                              &elemBlock_id, &species_id, file_name, optional_format);

      if (num_read_items < 5) {
        log_err("Insufficient DATA \n\t\"%s\".\n", data_line_buffer);
      }

      SPF(echo_string, "%s =  %s %d %d %d %s %s", "DATA", variable_name, node_set_id, elemBlock_id,
          species_id, file_name, optional_format);

      /* Read DATA  type: i.e. VELOCITY1, VOLTAGE, etc. */
      sprintf(input, "%s", variable_name);

      /* Search for appropriate integer identifier in Var_Name struct */

      /*OK the decryption:  pp->data[i]->data_type = 0.. MAX_VAR_TYPE */
      /* 	                     =>Primitive variable request */
      /*                    pp->data[i]->data_type = -1 */
      /* 		             =>Post proc variable request */
      /* 		    pp->data[i]->data_type= -2 */
      /*                            =>Cannot recognize the request. */

      pp_data[i]->data_type = -2; /* initialize first */
      for (k = 0; k < Num_Var_Names; k++) {
        if (!strncasecmp(input, Var_Name[k].name1, strlen(input)) ||
            !strncasecmp(input, Var_Name[k].name2, strlen(input))) {
          pp_data[i]->data_type = Var_Name[k].Index;
          strcpy(pp_data[i]->data_type_name, Var_Name[k].name1);
          break;
        }
      }

      if (pp_data[i]->data_type == -2) /*not a primary variable, then try post proc */
      {
        for (k = 0; k < Num_Post_Var_Names; k++) {
          if (!strncasecmp(input, Post_Var_Name[k].name1, strlen(input))) {
            pp_data[i]->data_type = -1; /* Flag Post var */
            strcpy(pp_data[i]->data_type_name, input);
            break;
          }
        }
      }

      if (pp_data[i]->data_type == -2) {
        log_err("Invalid DATA print variable \"%s\"; (see mm_names.h)\n", input);
      }

      if (node_set_id != -1) {
        pp_data[i]->ns_id = node_set_id;
      } else {
        log_err("Problem scanning post processing DATA node set id %s\n",
                pp_data[i]->data_type_name);
      }

      pp_data[i]->mat_num = -1;
      pp_data[i]->elem_blk_id = elemBlock_id;

      if (species_id != -1) {
        pp_data[i]->species_number = species_id;
      } else {
        log_err("Problem scanning post processing DATA species index %s\n",
                pp_data[i]->data_type_name);
      }

      if (1 != sscanf(file_name, "%s", filename)) {
        log_err("Problem scanning post processing DATA output filename %s\n",
                pp_data[i]->data_type_name);
      } else {
        strcpy(pp_data[i]->data_filenm, filename);
      }

      /*
       * Do NOT complain if the optional_format specifier is missing.
       */

      strncpy(pp_data[i]->format_flag, "\0\0\0\0\0\0\0\0", 8);

      if (optional_format[0] != '\0') {
        strcpy(pp_data[i]->format_flag, optional_format);
      }

      pp_data[i]->first_time = TRUE;

      ECHO(echo_string, echo_file);
    }
    ECHO("\nEND OF DATA\n", echo_file);
  }

  /*
   *  SCHEDULE POST-PROCESSING DATA SENSITIVITY CALCULATIONS, IF NEEDED
   *
   */

  iread = look_for_optional(ifp, "Post Processing Data Sensitivities", input, '=');

  /* count number of Post Processing Data specifications */
  if (iread == 1) {
    nn_post_data_sens = count_list(ifp, "DATA_SENS", input, '=', "END OF DATA_SENS");
  } else {
    nn_post_data_sens = 0;
  }

  /*
   *  Allocate memory to hold the data information
   */

  if (nn_post_data_sens > 0) {
    sz = sizeof(struct Post_Processing_Data_Sens *);
    pp_data_sens = (struct Post_Processing_Data_Sens **)array_alloc(1, nn_post_data_sens, sz);

    sz = sizeof(struct Post_Processing_Data_Sens);

    for (i = 0; i < nn_post_data_sens; i++) {
      pp_data_sens[i] = (struct Post_Processing_Data_Sens *)array_alloc(1, 1, sz);
    }
    /*Now load up information by reading cards */

    for (i = 0; i < nn_post_data_sens; i++) {
      look_for(ifp, "DATA_SENS", input, '=');
      /* Read DATA  type: i.e. VELOCITY1, VOLTAGE, etc. */
      if (fscanf(ifp, "%s", input) != 1) {
        GOMA_EH(GOMA_ERROR, "error reading DATA_SENS name");
      }
      /* Search for appropriate integer identifier in Var_Name struct */

      /*OK the decryption:  pp->data[i]->data_type = 0.. MAX_VAR_TYPE */
      /*                             =>Primitive variable request */
      /*                    pp->data[i]->data_type = -1 */
      /*                             =>Post proc variable request */
      /*                    pp->data[i]->data_type= -2 */
      /*                            =>Cannot recognize the request. */

      pp_data_sens[i]->data_type = -2; /* initialize first */
      for (k = 0; k < Num_Var_Names; k++) {
        if (!strcmp(input, Var_Name[k].name1)) {
          pp_data_sens[i]->data_type = Var_Name[k].Index;
          strcpy(pp_data_sens[i]->data_type_name, input);
          break;
        }
      }
      if (pp_data_sens[i]->data_type == -2) /*not a primary variable, then try post proc */
      {
        for (k = 0; k < Num_Post_Var_Names; k++) {
          GOMA_EH(GOMA_ERROR, "Unrecognized print variable. ");
          if (!strcmp(input, Post_Var_Name[k].name1))
            pp_data_sens[i]->data_type = -1; /* Ok this says it is a Post var */
          strcpy(pp_data_sens[i]->data_type_name, input);
        }
      }
      if (pp_data_sens[i]->data_type == -2)
        GOMA_EH(GOMA_ERROR, "Invalid choice of print_sens variable");

      /* read ns id */

      if (fscanf(ifp, "%d", &pp_data_sens[i]->ns_id) != 1) {
        fprintf(stderr, "%s:\tError reading data_sens->ns_id\n", yo);
        exit(-1);
      }

      /* read the material index number and then decrement its value */

      if (fscanf(ifp, "%d", &pp_data_sens[i]->mat_id) != 1) {
        fprintf(stderr, "%s:\tError reading data_sens->blk_id\n", yo);
        exit(-1);
      }
      /*	  pp_data_sens[i]->mat_id--;  */

      /* read species id*/

      if (fscanf(ifp, "%d", &pp_data_sens[i]->species_number) != 1) {
        fprintf(stderr, "%s:\tError reading data_sens->species_number\n", yo);
        exit(-1);
      }

      /*  read sensitivity variable type */

      if (fscanf(ifp, "%80s", input) != 1) {
        fprintf(stderr, "%s:\tError reading data_sens->sens_type (BC, MT, AC, UM, AN, or UF)\n",
                yo);
        exit(-1);
      }
      if (!strcmp(input, "BC")) {
        pp_data_sens[i]->sens_type = 1;
      } else if (!strcmp(input, "MT")) {
        pp_data_sens[i]->sens_type = 2;
      } else if (!strcmp(input, "AC")) {
        pp_data_sens[i]->sens_type = 3;
      } else if (!strcmp(input, "UM")) {
        pp_data_sens[i]->sens_type = 4;
      } else if (!strcmp(input, "UF")) {
        pp_data_sens[i]->sens_type = 5;
      } else if (!strcmp(input, "AN")) {
        pp_data_sens[i]->sens_type = 6;
      } else {
        fprintf(stderr, "%s:\tImproper set_type for data sensitivity - %s\n", yo, input);
        exit(-1);
      }

      /*  read BC id or material id */

      if (fscanf(ifp, "%d", &pp_data_sens[i]->sens_id) != 1) {
        fprintf(stderr, "%s:\tError reading data_sens->sens_id\n", yo);
        exit(-1);
      }
      if (pp_data_sens[i]->sens_type == 2 || pp_data_sens[i]->sens_type == 4)
        pp_data_sens[i]->sens_id--;

      /*  read data float or material property number */

      if (fscanf(ifp, "%d", &pp_data_sens[i]->sens_flt) != 1) {
        fprintf(stderr, "%s:\tError reading data_sens->sens_flt\n", yo);
        exit(-1);
      }

      if (pp_data_sens[i]->sens_type == 4) {
        if (fscanf(ifp, "%d", &pp_data_sens[i]->sens_flt2) != 1) {
          fprintf(stderr, "%s:\tError reading data_sens->sens_flt2\n", yo);
          exit(-1);
        }
      } else {
        pp_data_sens[i]->sens_flt2 = -1;
      }
      /*

      DETERMINE SENSITIVITY VECTOR NUMBER
      SEARCH FOR SAME PREVIOUS SENSITIVITY INFO

      */

      pp_data_sens[i]->vector_id = -1;

      /* first search flux sensitivity info */
      for (j = 0; j < nn_post_fluxes_sens; j++) {
        if (pp_data_sens[i]->sens_type == pp_fluxes_sens[j]->sens_type &&
            pp_data_sens[i]->sens_id == pp_fluxes_sens[j]->sens_id &&
            pp_data_sens[i]->sens_flt == pp_fluxes_sens[j]->sens_flt &&
            pp_data_sens[i]->sens_flt2 == pp_fluxes_sens[j]->sens_flt2) {
          pp_data_sens[i]->vector_id = pp_fluxes_sens[j]->vector_id;
        }
      }
      for (j = 0; j < i; j++) {
        if (pp_data_sens[i]->sens_type == pp_data_sens[j]->sens_type &&
            pp_data_sens[i]->sens_id == pp_data_sens[j]->sens_id &&
            pp_data_sens[i]->sens_flt == pp_data_sens[j]->sens_flt &&
            pp_data_sens[i]->sens_flt2 == pp_data_sens[j]->sens_flt2) {
          pp_data_sens[i]->vector_id = pp_data_sens[j]->vector_id;
        }
      }
      if (pp_data_sens[i]->vector_id == -1) {
        pp_data_sens[i]->vector_id = sens_vec_ct;

        /*

        IF FIRST ORDER CONTINUATION CHECK FOR SAME SENSITIVITY PARAMETER

        */

        if (Continuation == ALC_FIRST) {
          if (cont->upType == pp_data_sens[i]->sens_type) {
            int id1 = -1, id2 = -1;
            int id3 = -1;

            switch (cont->upType) {
            case 1:
            case 3:
              id1 = cont->upBCID;
              id2 = cont->upDFID;
              break;
            case 2:
              id1 = cont->upMTID;
              id2 = cont->upMPID;
              break;
            case 4:
              id1 = cont->upMTID;
              id2 = cont->upMPID;
              id3 = cont->upMDID;
              break;
            case 5:
              GOMA_EH(GOMA_ERROR, "sensitivities to UF not done");
              break;
            case 6:
              GOMA_EH(GOMA_ERROR, "sensitivities to AN not done");
              break;
            }
            if (id1 == pp_data_sens[i]->sens_id && id2 == pp_data_sens[i]->sens_flt &&
                id3 == pp_data_sens[i]->sens_flt2) {
              cont->sensvec_id = sens_vec_ct;
            }
          }
        }
        sens_vec_ct++;
      }

      /* read file name */
      read_string(ifp, first_string, '\n');
      nargs = sscanf(first_string, "%s", second_string);
      if (nargs == 0) {
        GOMA_EH(GOMA_ERROR, "Found zero arguments for the Data Sensitivity file name");
      }
      strcpy(pp_data_sens[i]->data_filenm, second_string);
    } /*   data card count loop  */
  }   /*   if data_sens conditions */

  /*
   *  SCHEDULE POST-PROCESSING PARTICLE TRACKING CALCULATIONS, IF NEEDED
   *
   */
  iread = look_for_optional(ifp, "Post Processing Particle Traces", input, '=');

  /* count number of Post Processing Particles */
  if (iread == 1) {
    nn_particles = count_list(ifp, "PARTICLE", input, '=', "END OF PARTICLES");
    ECHO("\nPost Processing Particle Traces =\n", echo_file);
  } else {
    nn_particles = 0;
  }

  /*
   *  Allocate memory to hold the particle information
   */

  if (nn_particles > 0) {
    sz = sizeof(struct Post_Processing_Particles *);
    pp_particles = (struct Post_Processing_Particles **)array_alloc(1, nn_particles, sz);

    sz = sizeof(struct Post_Processing_Particles);

    for (i = 0; i < nn_particles; i++) {
      pp_particles[i] = (struct Post_Processing_Particles *)array_alloc(1, 1, sz);
    }
    /*Now load up information by reading cards */
    for (i = 0; i < nn_particles; i++) {
      look_for(ifp, "PARTICLE", input, '=');

      if (fscanf(ifp, "%lf %lf %lf", &pp_particles[i]->coord[0], &pp_particles[i]->coord[1],
                 &pp_particles[i]->coord[2]) != 3) {
        GOMA_EH(GOMA_ERROR,
                "error reading particle input data (check to make sure you have 3 coord values");
      }
      SPF(echo_string, "%s = %g %g %g", "PARTICLE", pp_particles[i]->coord[0],
          pp_particles[i]->coord[1], pp_particles[i]->coord[2]);

      if (fscanf(ifp, "%lf %lf %lf", &pp_particles[i]->Start_Time, &pp_particles[i]->End_Time,
                 &pp_particles[i]->Delta_s) != 3) {
        GOMA_EH(GOMA_ERROR, "error reading particle time data (Start_time, End_Time, Delta_s)");
      }
      SPF(endofstring(echo_string), " %g %g %g", pp_particles[i]->Start_Time,
          pp_particles[i]->End_Time, pp_particles[i]->Delta_s);
      if (fscanf(ifp, "%lf %lf", &pp_particles[i]->mass, &pp_particles[i]->mobility) != 2) {
        GOMA_WH(-1, "defaulting to zero mass, mobility");
        pp_particles[i]->mass = 0.0;
        pp_particles[i]->mobility = 0.0;
      }
      SPF(endofstring(echo_string), " %g %g", pp_particles[i]->mass, pp_particles[i]->mobility);

      if (fscanf(ifp, "%lf %lf %lf", &pp_particles[i]->force[0], &pp_particles[i]->force[1],
                 &pp_particles[i]->force[2]) != 3) {
        GOMA_WH(-1, "defaulting to zero external force");
        pp_particles[i]->force[0] = 0.0;
        pp_particles[i]->force[1] = 0.0;
        pp_particles[i]->force[2] = 0.0;
      }
      SPF(endofstring(echo_string), " %g %g %g", pp_particles[i]->force[0],
          pp_particles[i]->force[1], pp_particles[i]->force[2]);

      read_string(ifp, first_string, '\n');
      nargs = sscanf(first_string, "%s", second_string);
      if (nargs == 0) {
        GOMA_EH(GOMA_ERROR, "Found zero arguments for the Particle file name");
      }
      strcpy(pp_particles[i]->filenm, second_string);
      SPF(endofstring(echo_string), " %s", pp_particles[i]->filenm);

      pp_particles[i]->Current_element_id = -1;
      ECHO(echo_string, echo_file);
    }

  } /*end if iread */
  ECHO("\nEND OF PARTICLES\n", echo_file);
  /*
   *  SCHEDULE POST-PROCESSING VOLUME INTEGRALS, IF NEEDED
   *
   */
  iread = look_for_optional(ifp, "Post Processing Volumetric Integration", input, '=');

  if (iread == 1) {
    nn_volume = count_list(ifp, "VOLUME_INT", input, '=', "END OF VOLUME_INT");
    ECHO("\nPost Processing Volumetric Integration =\n", echo_file);
  } else {
    nn_volume = 0;
  }

  if (nn_volume) {
    ppvi_type = PPVI_VERBOSE; // Default to default output behavior see "Volumetric Integration
                              // Output Format" for other possible behaviors
    sz = sizeof(struct Post_Processing_Volumetric *);

    pp_volume = (struct Post_Processing_Volumetric **)array_alloc(1, nn_volume, sz);

    sz = sizeof(struct Post_Processing_Volumetric);

    for (i = 0; i < nn_volume; i++) {
      pp_volume[i] = (struct Post_Processing_Volumetric *)array_alloc(1, 1, sz);
    }

    for (i = 0; i < nn_volume; i++) {
      look_for(ifp, "VOLUME_INT", input, '=');

      if (fscanf(ifp, "%s", ts) != 1) {
        GOMA_EH(GOMA_ERROR, "error reading Post Processing volume integral card \n");
      }

      pp_volume[i]->volume_type = -1;

      for (k = 0; k < Num_Vol_Names; k++) {
        if (!strcmp(ts, pp_vol_names[k].name)) {
          pp_volume[i]->volume_type = pp_vol_names[k].Index;
          strcpy(pp_volume[i]->volume_name, ts);
        }
      }
      SPF(echo_string, "%s = %s", "VOLUME_INT", pp_volume[i]->volume_name);

      if (fscanf(ifp, "%d %d", &(pp_volume[i]->blk_id), &(pp_volume[i]->species_no)) != 2) {
        fprintf(stderr, "%s:\tError reading blk_id or species number for volume integral\n", yo);
        exit(-1);
      }

      SPF(endofstring(echo_string), " %d %d", pp_volume[i]->blk_id, pp_volume[i]->species_no);

      /* read file name TAB asks RBS why read the filename this waY ?
         I don't know, I copied it from somewhere - RBS
         Actually, Randy probably had me do it this way
         -how's that for passing the buck.

         Ah, that makes sense.  Randy wasn't a nerd in high school.
         Therefore, his code is always a bit off.  Not his fault
         really. (I kid, I kid )  */

      read_string(ifp, first_string, '\n');

      nargs = sscanf(first_string, "%s", second_string);

      if (nargs == 0) {
        GOMA_EH(GOMA_ERROR, "Found zero arguments for the Volumetric Integration file name");
      }

      strcpy(pp_volume[i]->volume_fname, second_string);

      SPF(endofstring(echo_string), " %s", pp_volume[i]->volume_fname);

      {
        int num_const, k;
        char *args[MAX_NUMBER_PARAMS];

        if ((num_const = tokenize_by_whsp(first_string, args, MAX_NUMBER_PARAMS)) > 1) {
          pp_volume[i]->params = alloc_dbl_1(num_const - 1, 0.0);

          for (k = 1; k < num_const; k++) {
            pp_volume[i]->params[k - 1] = atof(args[k]);
          }

          pp_volume[i]->num_params = num_const - 1;

          SPF_DBL_VEC(endofstring(echo_string), num_const - 1, &(pp_volume[i]->params[0]));

        } else {
          pp_volume[i]->params = NULL;
          pp_volume[i]->num_params = 0;
        }
      }
      ECHO(echo_string, echo_file);

    } /* end of i < nn_volume */
    ECHO("\nEND OF VOLUME_INT\n", echo_file);
  } /*if iread */

  iread = look_for_optional(ifp, "Volumetric Integration Output Format", input, '=');
  if (iread == 1) {
    FILE *amcfp;
    int i;
    if (fscanf(ifp, "%s", ts) != 1) {
      GOMA_EH(GOMA_ERROR, "error reading Volume Integration Output Format card \n");
    }
    if (!strcasecmp(ts, "Verbose")) {
      ppvi_type = PPVI_VERBOSE;
    } else if (!strcasecmp(ts, "CSV")) {
      ppvi_type = PPVI_CSV;

      /* in order to prepend output file with data type info */
      char ts1[MAX_CHAR_IN_INPUT];
      for (i = 0; i < nn_volume; i++) {
        amcfp = fopen(pp_volume[i]->volume_fname, "a");
        if (amcfp != NULL) {
          strcpy(ts1, pp_volume[i]->volume_name);
          fprintf(amcfp, "Time, %s\n", ts1);
          fflush(amcfp);
          fclose(amcfp);
        }
      }
    } else {
      GOMA_WH(-1, "The Volumetric Integration Output Format was not recognized");
      ppvi_type = PPVI_VERBOSE;
    }
  }

  /*
   *  SCHEDULE POST-PROCESSING GLOBAL CALCULATIONS, IF NEEDED
   */

  iread = look_for_optional(ifp, "Post Processing Global", input, '=');

  /* count number of post-processing flux calculation specifications */

  if (iread == 1) {
    nn_global = count_list(ifp, "GLOBAL", input, '=', "END OF GLOBAL");
    ECHO("\nPost Processing Global =\n", echo_file);
  } else {
    nn_global = 0;
  }

  /*
   *  Allocate memory to hold the flux information
   */

  if (nn_global > 0) {
    sz = sizeof(struct Post_Processing_Global *);
    pp_global = (struct Post_Processing_Global **)array_alloc(1, nn_global, sz);

    sz = sizeof(struct Post_Processing_Global);

    for (i = 0; i < nn_global; i++) {
      pp_global[i] = (struct Post_Processing_Global *)array_alloc(1, 1, sz);
    }

    /*Now load up information by reading cards */

    for (i = 0; i < nn_global; i++) {
      look_for(ifp, "GLOBAL", input, '=');

      /* Read GLOBAL  type */

      if (fscanf(ifp, "%s", ts) != 1) {
        GOMA_EH(GOMA_ERROR, "error reading Post Processing Global input variable");
      }

      pp_global[i]->type = -1;
      for (k = 0; k < PP_GLOBAL_COUNT; k++) {
        if (!strcmp(ts, pp_global_names[k].name)) {
          pp_global[i]->type = pp_global_names[k].Index;
          strcpy(pp_global[i]->type_name, ts);
          break;
        }
      }
      if (pp_global[i]->type == -1) {
        GOMA_EH(GOMA_ERROR, "Invalid Global name");
      }
      SPF(echo_string, "%s = %s", "GLOBAL", pp_global[i]->type_name);

      /* Custom post processing card reading for global */
      switch (pp_global[i]->type) {
        /* Read in a file name */
      case PP_GLOBAL_LS_INTERFACE_PRINT:
      case PP_GLOBAL_LS_INTERFACE_PRINT_ALL_TIMES: {

        /* read file name */
        if (fscanf(ifp, "%s", ts) != 1) {
          GOMA_EH(GOMA_ERROR, "error reading Post Processing Global filename");
        }
        strcpy(pp_global[i]->filenm, ts);
        SPF(endofstring(echo_string), " %s", pp_global[i]->filenm);
      } break;
      default:
        GOMA_EH(GOMA_ERROR, "Unsupported GLOBAL Post Processing Type");
        break;
      }

      ECHO(echo_string, echo_file);
    }
    ECHO("\nEND OF GLOBAL\n", echo_file);
  }

  iread = look_for_optional(ifp, "Post Processing Averages", input, '=');

  /* count number of post-processing flux calculation specifications */

  if (iread == 1) {
    nn_average = count_list(ifp, "AVERAGE", input, '=', "END OF AVERAGES");
    ECHO("\nPost Processing Averages =\n", echo_file);
  } else {
    nn_average = 0;
  }

  if (nn_average > 0) {
    sz = sizeof(pp_Average *);
    pp_average = (pp_Average **)array_alloc(1, nn_average, sz);

    sz = sizeof(pp_Average);

    for (i = 0; i < nn_average; i++) {
      pp_average[i] = (pp_Average *)array_alloc(1, 1, sz);
    }

    /*Now load up information by reading cards */

    for (i = 0; i < nn_average; i++) {
      look_for(ifp, "AVERAGE", input, '=');
      save_position = ftell(ifp);
      char *fgetsret = fgets(data_line_buffer, MAX_CHAR_IN_INPUT, ifp);
      if (fgetsret == NULL) {
        GOMA_EH(GOMA_ERROR, "Error reading post processing line in Averages");
      }

      fseek(ifp, save_position, SEEK_SET);

      pp_average[i]->type = -1;
      pp_average[i]->non_variable_type = 0;

      int read_average_items =
          sscanf(data_line_buffer, "%s %d", variable_name, &(pp_average[i]->species_index));

      if (read_average_items < 2) {
        GOMA_EH(GOMA_ERROR, "Error in reading Average post_processing");
      }

      for (k = 0; k < Num_Var_Names; k++) {
        int st1 = strlen(variable_name);
        int st2 = strlen(Var_Name[k].name1);
        int st3 = strlen(Var_Name[k].name2);
        if ((st1 == st2 && !strncasecmp(variable_name, Var_Name[k].name1, strlen(variable_name))) ||
            (st1 == st3 && !strncasecmp(variable_name, Var_Name[k].name2, strlen(variable_name)))) {
          pp_average[i]->type = Var_Name[k].Index;
          if (pp_average[i]->type == MASS_FRACTION) {
            int err = snprintf(pp_average[i]->type_name, MAX_VAR_NAME_LNGTH, "%s%d%s",
                               Var_Name[k].name2, pp_average[i]->species_index, "_AVG");
            if (err < 0 || err >= MAX_VAR_NAME_LNGTH) {
              GOMA_EH(GOMA_ERROR, "Error writing mass fraction variable for average");
            }
          } else {
            strcpy(pp_average[i]->type_name, Var_Name[k].name2);
            strcat(pp_average[i]->type_name, "_AVG");
          }
          break;
        }
      }

      if (pp_average[i]->type == -1) {
        if (!strncasecmp(variable_name, "DENSITY_AVG", strlen(variable_name))) {
          strcpy(pp_average[i]->type_name, "DENSITY_AVG");
          pp_average[i]->non_variable_type = 1;
          pp_average[i]->type = AVG_HEAVISIDE;
        } else if (!strncasecmp(variable_name, "HEAVISIDE", strlen(variable_name))) {
          strcpy(pp_average[i]->type_name, "HEAVISIDE_AVG");
          pp_average[i]->non_variable_type = 1;
          pp_average[i]->type = AVG_HEAVISIDE;
        } else if (!strncasecmp(variable_name, "VISCOSITY", strlen(variable_name))) {
          strcpy(pp_average[i]->type_name, "VISCOSITY_AVG");
          pp_average[i]->non_variable_type = 1;
          pp_average[i]->type = AVG_VISCOSITY;
        } else if (!strncasecmp(variable_name, "SHEARRATE", strlen(variable_name))) {
          strcpy(pp_average[i]->type_name, "SHEARRATE_AVG");
          pp_average[i]->non_variable_type = 1;
          pp_average[i]->type = AVG_SHEAR;
        } else if (!strncasecmp(variable_name, "EM_SCAT_MAG", strlen(variable_name))) {
          strcpy(pp_average[i]->type_name, "EM_SCAT_MAG");
          pp_average[i]->non_variable_type = 1;
          pp_average[i]->type = AVG_EM_SCAT_MAG;
        } else if (!strncasecmp(variable_name, "EM_MAG", strlen(variable_name))) {
          strcpy(pp_average[i]->type_name, "EM_MAG");
          pp_average[i]->non_variable_type = 1;
          pp_average[i]->type = AVG_EM_MAG;
        } else if (!strncasecmp(variable_name, "EM_INC_MAG", strlen(variable_name))) {
          strcpy(pp_average[i]->type_name, "EM_INC_MAG");
          pp_average[i]->non_variable_type = 1;
          pp_average[i]->type = AVG_EM_INC_MAG;
        } else if (!strncasecmp(variable_name, "EM", strlen(variable_name))) {
          int new_items = 6 - 1;
          pp_Average **pp_average_tmp = pp_average;
          sz = sizeof(pp_Average *);
          pp_average = (pp_Average **)array_alloc(1, nn_average + new_items, sz);
          for (int k = 0; k < nn_average; k++) {
            pp_average[k] = pp_average_tmp[k];
          }
          free(pp_average_tmp);
          sz = sizeof(pp_Average);
          for (int k = nn_average; k < (new_items + nn_average); k++) {
            pp_average[k] = (pp_Average *)array_alloc(1, 1, sz);
          }
          nn_average += new_items;
          strcpy(pp_average[i]->type_name, "PP_EMR_X");
          pp_average[i]->non_variable_type = 1;
          pp_average[i]->type = AVG_EMR_X;
          pp_average[i]->species_index = 0;
          i++;
          strcpy(pp_average[i]->type_name, "PP_EMR_Y");
          pp_average[i]->non_variable_type = 1;
          pp_average[i]->type = AVG_EMR_Y;
          i++;
          strcpy(pp_average[i]->type_name, "PP_EMR_Z");
          pp_average[i]->non_variable_type = 1;
          pp_average[i]->type = AVG_EMR_Z;
          i++;
          strcpy(pp_average[i]->type_name, "PP_EMI_X");
          pp_average[i]->non_variable_type = 1;
          pp_average[i]->type = AVG_EMI_X;
          i++;
          strcpy(pp_average[i]->type_name, "PP_EMI_Y");
          pp_average[i]->non_variable_type = 1;
          pp_average[i]->type = AVG_EMI_Y;
          i++;
          strcpy(pp_average[i]->type_name, "PP_EMI_Z");
          pp_average[i]->non_variable_type = 1;
          pp_average[i]->type = AVG_EMI_Z;
        } else if (!strncasecmp(variable_name, "EMSCAT", strlen(variable_name))) {
          int new_items = 6 - 1;
          pp_Average **pp_average_tmp = pp_average;
          sz = sizeof(pp_Average *);
          pp_average = (pp_Average **)array_alloc(1, nn_average + new_items, sz);
          for (int k = 0; k < nn_average; k++) {
            pp_average[k] = pp_average_tmp[k];
          }
          free(pp_average_tmp);
          sz = sizeof(pp_Average);
          for (int k = nn_average; k < (new_items + nn_average); k++) {
            pp_average[k] = (pp_Average *)array_alloc(1, 1, sz);
          }
          nn_average += new_items;
          strcpy(pp_average[i]->type_name, "PP_EMSCATR_X");
          pp_average[i]->non_variable_type = 1;
          pp_average[i]->type = AVG_EMSCATR_X;
          pp_average[i]->species_index = 0;
          i++;
          strcpy(pp_average[i]->type_name, "PP_EMSCATR_Y");
          pp_average[i]->non_variable_type = 1;
          pp_average[i]->type = AVG_EMSCATR_Y;
          i++;
          strcpy(pp_average[i]->type_name, "PP_EMSCATR_Z");
          pp_average[i]->non_variable_type = 1;
          pp_average[i]->type = AVG_EMSCATR_Z;
          i++;
          strcpy(pp_average[i]->type_name, "PP_EMSCATI_X");
          pp_average[i]->non_variable_type = 1;
          pp_average[i]->type = AVG_EMSCATI_X;
          i++;
          strcpy(pp_average[i]->type_name, "PP_EMSCATI_Y");
          pp_average[i]->non_variable_type = 1;
          pp_average[i]->type = AVG_EMSCATI_Y;
          i++;
          strcpy(pp_average[i]->type_name, "PP_EMSCATI_Z");
          pp_average[i]->non_variable_type = 1;
          pp_average[i]->type = AVG_EMSCATI_Z;
        } else if (!strncasecmp(variable_name, "EMINC", strlen(variable_name))) {
          int new_items = 6 - 1;
          pp_Average **pp_average_tmp = pp_average;
          sz = sizeof(pp_Average *);
          pp_average = (pp_Average **)array_alloc(1, nn_average + new_items, sz);
          for (int k = 0; k < nn_average; k++) {
            pp_average[k] = pp_average_tmp[k];
          }
          free(pp_average_tmp);
          sz = sizeof(pp_Average);
          for (int k = nn_average; k < (new_items + nn_average); k++) {
            pp_average[k] = (pp_Average *)array_alloc(1, 1, sz);
          }
          nn_average += new_items;
          strcpy(pp_average[i]->type_name, "PP_EMINCR_X");
          pp_average[i]->non_variable_type = 1;
          pp_average[i]->type = AVG_EMINCR_X;
          pp_average[i]->species_index = 0;
          i++;
          strcpy(pp_average[i]->type_name, "PP_EMINCR_Y");
          pp_average[i]->non_variable_type = 1;
          pp_average[i]->type = AVG_EMINCR_Y;
          i++;
          strcpy(pp_average[i]->type_name, "PP_EMINCR_Z");
          pp_average[i]->non_variable_type = 1;
          pp_average[i]->type = AVG_EMINCR_Z;
          i++;
          strcpy(pp_average[i]->type_name, "PP_EMINCI_X");
          pp_average[i]->non_variable_type = 1;
          pp_average[i]->type = AVG_EMINCI_X;
          i++;
          strcpy(pp_average[i]->type_name, "PP_EMINCI_Y");
          pp_average[i]->non_variable_type = 1;
          pp_average[i]->type = AVG_EMINCI_Y;
          i++;
          strcpy(pp_average[i]->type_name, "PP_EMINCI_Z");
          pp_average[i]->non_variable_type = 1;
          pp_average[i]->type = AVG_EMINCI_Z;
        } else {
          fprintf(stderr, "Error reading unknown variable type: %s\n", variable_name);
          GOMA_EH(GOMA_ERROR, "Unknown variable type for post processing");
        }
      }

      SPF(echo_string, "%s = %s", "AVERAGE", variable_name);
      SPF(endofstring(echo_string), " %d", pp_average[i]->species_index);
      ECHO(echo_string, echo_file);
    }
    ECHO("\nEND OF AVERAGE\n", echo_file);
  }
}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

static int look_for_post_proc(FILE *ifp,           /* pointer to file                           */
                              char *search_string, /* input string for which to search  */
                              int *flag_variable)  /* integer flag for post proc options */

/*
 * look_for_post_proc -- generic read post processing specifications for
 *                       each entry. The call in this code specifically
 *                       rewinds the file so that searches may occur over
 *                       the entire file.
 *
 */
{
  char input[MAX_CHAR_IN_INPUT];
  int iread;
  char echo_string[MAX_CHAR_ECHO_INPUT] = "\0";
  char echo_file[MAX_CHAR_IN_INPUT] = "\0";

  strcpy(echo_file, Echo_Input_File);

  /*
   * Rewind file and look for the optional string anywhere in the input deck
   * PP vars to be exported are indicated by flag_variable=2.
   */
  iread = look_for_optional(ifp, search_string, input, '=');
  if (iread == 1) {
    if (fscanf(ifp, "%s", input) != 1) {
      strip(input);
    }
    if (strcmp(input, "yes") == 0) {
      *flag_variable = 1;
    } else if (strcmp(input, "exp") == 0) {
      *flag_variable = 2;
      Num_Export_XP++;
    } else if (strcmp(input, "no") == 0) {
      *flag_variable = -1;
    } else {
      GOMA_EH(GOMA_ERROR, "Unknown problem option");
    }
    SPF(echo_string, "%s = %s", search_string, input);
    ECHO(echo_string, echo_file);
  }

#ifndef LIBRARY_MODE
  if (Num_Export_XP > 0) {
    GOMA_WH(-1, "Post-process exp option only available in LIBRARY_MODE!");
    Num_Export_XP = 0;
  }
#endif
  return iread;
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
/* load_nodal_tkn -- load nodal variable types, kinds, and names
 *
 * input arguments:
 *	rd		pointer to results description structure that is
 *			loaded up here.
 *
 * output arguments:
 * 	nvtype[]	nodal variable type (0=VELOCITY1, etc.)
 *	nvkind[] 	nodal variable kind (0, except for mass species)
 *	nvname[] 	nodal variable name ("VX", etc.)
 *	nvunit[] 	nodal variable unit ("kg/s", etc.)
 *	nvdesc[] 	nodal variable description
 *
 * return value:
 *	int	       -1	something went wrong in this routine
 *			0	everything went OK fine hunky-dory
 *
 * Author:	Philip A. Sackinger
 * Created:	Tue Mar 23 08:00:13 MST 1993
 * Revised:	Wed Apr  7 06:46:16 MDT 1993
 */

int load_nodal_tkn(struct Results_Description *rd, int *tnv, int *tnv_post) {
  int i = -1, k, index, index_post, index_post_export, status = 0, check, w;
  int v_s[MAX_MODES][DIM][DIM];
  int a, b, mode;
  MATRL_PROP_STRUCT *matrl = mp_glob[0];
  VARIABLE_DESCRIPTION_STRUCT *vd;

  /*
   * Temporary strings hold current species name (short & long)...
   */
  char nm[MAX_VAR_NAME_LNGTH], nm1[MAX_VAR_NAME_LNGTH], nm2[MAX_VAR_NAME_LNGTH];
  char ds[MAX_VAR_DESC_LNGTH], ds1[MAX_VAR_DESC_LNGTH], ds2[MAX_VAR_DESC_LNGTH];
  char species_name[MAX_VAR_NAME_LNGTH];
  char species_desc[MAX_VAR_DESC_LNGTH];
  int post_flag, var;
  int mn; /* Material Number ID */

  /* load eqn and variable number in tensor form */

  v_s[0][0][0] = POLYMER_STRESS11;
  v_s[0][0][1] = POLYMER_STRESS12;
  v_s[0][0][2] = POLYMER_STRESS13;
  v_s[0][1][0] = POLYMER_STRESS12;
  v_s[0][1][1] = POLYMER_STRESS22;
  v_s[0][1][2] = POLYMER_STRESS23;
  v_s[0][2][0] = POLYMER_STRESS13;
  v_s[0][2][1] = POLYMER_STRESS23;
  v_s[0][2][2] = POLYMER_STRESS33;

  v_s[1][0][0] = POLYMER_STRESS11_1;
  v_s[1][0][1] = POLYMER_STRESS12_1;
  v_s[1][0][2] = POLYMER_STRESS13_1;
  v_s[1][1][0] = POLYMER_STRESS12_1;
  v_s[1][1][1] = POLYMER_STRESS22_1;
  v_s[1][1][2] = POLYMER_STRESS23_1;
  v_s[1][2][0] = POLYMER_STRESS13_1;
  v_s[1][2][1] = POLYMER_STRESS23_1;
  v_s[1][2][2] = POLYMER_STRESS33_1;

  v_s[2][0][0] = POLYMER_STRESS11_2;
  v_s[2][0][1] = POLYMER_STRESS12_2;
  v_s[2][0][2] = POLYMER_STRESS13_2;
  v_s[2][1][0] = POLYMER_STRESS12_2;
  v_s[2][1][1] = POLYMER_STRESS22_2;
  v_s[2][1][2] = POLYMER_STRESS23_2;
  v_s[2][2][0] = POLYMER_STRESS13_2;
  v_s[2][2][1] = POLYMER_STRESS23_2;
  v_s[2][2][2] = POLYMER_STRESS33_2;

  v_s[3][0][0] = POLYMER_STRESS11_3;
  v_s[3][0][1] = POLYMER_STRESS12_3;
  v_s[3][0][2] = POLYMER_STRESS13_3;
  v_s[3][1][0] = POLYMER_STRESS12_3;
  v_s[3][1][1] = POLYMER_STRESS22_3;
  v_s[3][1][2] = POLYMER_STRESS23_3;
  v_s[3][2][0] = POLYMER_STRESS13_3;
  v_s[3][2][1] = POLYMER_STRESS23_3;
  v_s[3][2][2] = POLYMER_STRESS33_3;

  v_s[4][0][0] = POLYMER_STRESS11_4;
  v_s[4][0][1] = POLYMER_STRESS12_4;
  v_s[4][0][2] = POLYMER_STRESS13_4;
  v_s[4][1][0] = POLYMER_STRESS12_4;
  v_s[4][1][1] = POLYMER_STRESS22_4;
  v_s[4][1][2] = POLYMER_STRESS23_4;
  v_s[4][2][0] = POLYMER_STRESS13_4;
  v_s[4][2][1] = POLYMER_STRESS23_4;
  v_s[4][2][2] = POLYMER_STRESS33_4;

  v_s[5][0][0] = POLYMER_STRESS11_5;
  v_s[5][0][1] = POLYMER_STRESS12_5;
  v_s[5][0][2] = POLYMER_STRESS13_5;
  v_s[5][1][0] = POLYMER_STRESS12_5;
  v_s[5][1][1] = POLYMER_STRESS22_5;
  v_s[5][1][2] = POLYMER_STRESS23_5;
  v_s[5][2][0] = POLYMER_STRESS13_5;
  v_s[5][2][1] = POLYMER_STRESS23_5;
  v_s[5][2][2] = POLYMER_STRESS33_5;

  v_s[6][0][0] = POLYMER_STRESS11_6;
  v_s[6][0][1] = POLYMER_STRESS12_6;
  v_s[6][0][2] = POLYMER_STRESS13_6;
  v_s[6][1][0] = POLYMER_STRESS12_6;
  v_s[6][1][1] = POLYMER_STRESS22_6;
  v_s[6][1][2] = POLYMER_STRESS23_6;
  v_s[6][2][0] = POLYMER_STRESS13_6;
  v_s[6][2][1] = POLYMER_STRESS23_6;
  v_s[6][2][2] = POLYMER_STRESS33_6;

  v_s[7][0][0] = POLYMER_STRESS11_7;
  v_s[7][0][1] = POLYMER_STRESS12_7;
  v_s[7][0][2] = POLYMER_STRESS13_7;
  v_s[7][1][0] = POLYMER_STRESS12_7;
  v_s[7][1][1] = POLYMER_STRESS22_7;
  v_s[7][1][2] = POLYMER_STRESS23_7;
  v_s[7][2][0] = POLYMER_STRESS13_7;
  v_s[7][2][1] = POLYMER_STRESS23_7;
  v_s[7][2][2] = POLYMER_STRESS33_7;

  index = 0;
  index_post = 0;
  index_post_export = 0;

  /*
   *  Here, if a variable and equation is turned on in any one
   *  material, then we must make provisions in the rd structure
   * for all materials.  This is for the sake of the post-processor.
   *
   *
   * Check each variable to see if it has a valid interpolation and
   * write it as a nodal variable if it has this type in any phase.
   * We check here to see if the variable has an interpolation
   * type which is amenable to extraction via extract_nodal_variable().
   * If the variable exists in the problem but isn't continuous
   * or doesn't just have one value per node, (aka pressure, often).
   * Then, we don't include it in the list.
   */

  /*  This is to make sure that the mesh displacements are first
   *  in the exodus file. Just a little something to keep BLOT
   *  happy, so it will show our deformed meshes.
   */

  for (var = MESH_DISPLACEMENT1; var < (MESH_DISPLACEMENT3 + 1); var++) {
    if (variable_type_nodalInterp(var)) {
      set_nv_tkud(rd, index, var, 0, -1, Var_Name[var].name2, "[1]", Var_Name[var].name1, FALSE);
      index++;
    }
  }

  /*
   *  Loop over all of the variables defined in the problem,
   *   (except mesh displacement variables, which are handled up above)
   *  and define their Results_Description structure entries.
   *  These will be used to define the associated exodus variables.
   */
  for (var = V_FIRST; var < V_LAST; var++) {
    if ((var < MESH_DISPLACEMENT1) || (var > MESH_DISPLACEMENT3)) {
      if (variable_type_nodalInterp(var)) {
        if (var == MASS_FRACTION) {
          /*
           * HKM:
           *  Loop over the materials defined in the problem starting
           *  with the generic material, -1. Search for variables
           *  defined in the solution vector. If there are any,
           *  then add them to the structure that defines what the
           *  Exodus output is.
           */
          for (mn = -1; mn < upd->Num_Mat; mn++) {
            if (mn == -1) {
              for (i = upd->Num_Mat - 1; i >= 0; i--) {
                if (mp_glob[i]->Num_Species == upd->Max_Num_Species) {
                  matrl = mp_glob[i];
                }
              }
            } else {
              matrl = mp_glob[mn];
            }
            for (i = 0; i < matrl->Num_Species; i++) {
              vd = get_vd_ptr(MASS_FRACTION, mn, i);
              if (vd && (vd->MatID == mn)) {
                /*
                 * Assign the exodus species name
                 */
                assign_species_name(i, matrl, species_name, species_desc, mn);
                /*
                 * Set the values in the Results_Description structure
                 */
                set_nv_tkud(rd, index, var, i, mn, species_name, "[1]", species_desc, FALSE);
                index++;
              }
            }
          }
        } else {
          for (mn = -1; mn < upd->Num_Mat; mn++) {
            if (mn == -1) {
              for (i = upd->Num_Mat - 1; i >= 0; i--) {
                if (pd_glob[i]->i[pg->imtrx][var]) {
                  matrl = mp_glob[i];
                }
              }
            } else {
              matrl = mp_glob[mn];
            }
            vd = get_vd_ptr(var, mn, i);
            if (vd && (vd->MatID == mn)) {
              assign_var_name(var, 0, matrl, species_name, species_desc, mn);
              set_nv_tkud(rd, index, var, 0, mn, species_name, "[1]", species_desc, FALSE);
              index++;
            }
          }
        }
      }
    }
  }
  /*
   *  Set the total number of Output Nodal Vectors determined by
   *  extract_nodal_vec() at this point
   */
  rd->TotalNVSolnOutput = index;

  /*
   * Start checking for post-processing options
   *  at this point they should contain 1 = yes or -1 = no
   *  convert to -1 = no and post-processing variable number = yes
   */

  if (STREAM != -1 && Num_Var_In_Type[pg->imtrx][VELOCITY1]) {
    if (Num_Dim == 3) {
      GOMA_WH(-1, "Cant do stream function in 3D");
      STREAM = -1;
    } else {
      set_nv_tkud(rd, index, 0, 0, -2, "STREAM", "[1]", "Stream Function", FALSE);
      index++;
      if (STREAM == 2) {
        Export_XP_ID[index_post_export] = index_post;
        index_post_export++;
      }
      STREAM = index_post;
      index_post++;
      /* Put the index also in the table.  We should condense STREAM
         and the corresponding entry in Post_Var_Name struct to be the same */
      for (k = 0; k < Num_Post_Var_Names; k++) {
        if (!strcmp("STREAM", Post_Var_Name[k].name1))
          Post_Var_Name[k].Index = index_post;
      }
    }
  }

  if (STREAM_NORMAL_STRESS != -1 && Num_Var_In_Type[pg->imtrx][R_MOMENTUM1]) {
    set_nv_tkud(rd, index, 0, 0, -2, "SNS", "[1]", "Streamwise normal stress", FALSE);
    index++;
    if (STREAM_NORMAL_STRESS == 2) {
      Export_XP_ID[index_post_export] = index_post;
      index_post_export++;
    }
    STREAM_NORMAL_STRESS = index_post;
    index_post++;
  }

  if (STREAM_SHEAR_STRESS != -1 && Num_Var_In_Type[pg->imtrx][R_MOMENTUM1]) {
    set_nv_tkud(rd, index, 0, 0, -2, "SSS", "[1]", "Streamwise shear stress", FALSE);
    index++;
    if (STREAM_SHEAR_STRESS == 2) {
      Export_XP_ID[index_post_export] = index_post;
      index_post_export++;
    }
    STREAM_SHEAR_STRESS = index_post;
    index_post++;
  }

  if (STREAM_TENSION != -1 && Num_Var_In_Type[pg->imtrx][R_MOMENTUM1]) {
    set_nv_tkud(rd, index, 0, 0, -2, "SSD", "[1]", "Streamwise Stress Difference", FALSE);
    index++;
    if (STREAM_TENSION == 2) {
      Export_XP_ID[index_post_export] = index_post;
      index_post_export++;
    }
    STREAM_TENSION = index_post;
    index_post++;
  }

  if (DIV_VELOCITY != -1 && Num_Var_In_Type[pg->imtrx][PRESSURE]) {
    set_nv_tkud(rd, index, 0, 0, -2, "DIVV", "[1]", "Divergence of Velocity", FALSE);
    index++;
    if (DIV_VELOCITY == 2) {
      Export_XP_ID[index_post_export] = index_post;
      index_post_export++;
    }
    DIV_VELOCITY = index_post;
    index_post++;
  }

  if (DIV_PVELOCITY != -1 && Num_Var_In_Type[pg->imtrx][R_PMOMENTUM1]) {
    set_nv_tkud(rd, index, 0, 0, -2, "DIVPV", "[1]", "Divergence of Particle Velocity", FALSE);
    index++;
    if (DIV_PVELOCITY == 2) {
      Export_XP_ID[index_post_export] = index_post;
      index_post_export++;
    }
    DIV_PVELOCITY = index_post;
    index_post++;
  }

  if (DIV_TOTAL != -1 && Num_Var_In_Type[pg->imtrx][R_PMOMENTUM1]) {
    set_nv_tkud(rd, index, 0, 0, -2, "DIVTOTAL", "[1]", "Divergence of Total Velocity", FALSE);
    index++;
    if (DIV_TOTAL == 2) {
      Export_XP_ID[index_post_export] = index_post;
      index_post_export++;
    }
    DIV_TOTAL = index_post;
    index_post++;
  }

  if (PP_Viscosity != -1 &&
      (Num_Var_In_Type[pg->imtrx][R_MOMENTUM1] || Num_Var_In_Type[pg->imtrx][R_LUBP] ||
       Num_Var_In_Type[pg->imtrx][R_SHELL_FILMP])) {
    set_nv_tkud(rd, index, 0, 0, -2, "MU", "[1]", "Viscosity", FALSE);
    index++;
    if (PP_Viscosity == 2) {
      Export_XP_ID[index_post_export] = index_post;
      index_post_export++;
    }
    PP_Viscosity = index_post;
    index_post++;
  }

  if (PP_FlowingLiquid_Viscosity != -1 && Num_Var_In_Type[pg->imtrx][R_MOMENTUM1]) {
    set_nv_tkud(rd, index, 0, 0, -2, "FLOW_MU", "[1]", "FlowingLiquid Viscosity", FALSE);
    index++;
    if (PP_FlowingLiquid_Viscosity == 2) {
      Export_XP_ID[index_post_export] = index_post;
      index_post_export++;
    }
    PP_FlowingLiquid_Viscosity = index_post;
    index_post++;
  }

  if (PP_VolumeFractionGas != -1 && Num_Var_In_Type[pg->imtrx][R_MOMENTUM1]) {
    set_nv_tkud(rd, index, 0, 0, -2, "VOLUMEFRACTIONGAS", "[1]", "Volume Fraction of Gas Phase",
                FALSE);
    index++;
    if (PP_VolumeFractionGas == 2) {
      Export_XP_ID[index_post_export] = index_post;
      index_post_export++;
    }
    PP_VolumeFractionGas = index_post;
    index_post++;
  }

  if (DENSITY != -1 && Num_Var_In_Type[pg->imtrx][R_MOMENTUM1]) {
    set_nv_tkud(rd, index, 0, 0, -2, "RHO", "[1]", "Density", FALSE);
    index++;
    if (DENSITY == 2) {
      Export_XP_ID[index_post_export] = index_post;
      index_post_export++;
    }
    DENSITY = index_post;
    index_post++;
  }

  if (HEAVISIDE != -1 && ls != NULL && Num_Var_In_Type[pg->imtrx][R_FILL]) {
    set_nv_tkud(rd, index, 0, 0, -2, "HEAVISIDE", "[1]", "HEAVISIDE", FALSE);
    index++;
    if (HEAVISIDE == 2) {
      Export_XP_ID[index_post_export] = index_post;
      index_post_export++;
    }
    HEAVISIDE = index_post;
    index_post++;
  }

  if (RHO_DOT != -1 && Num_Var_In_Type[pg->imtrx][R_MOMENTUM1]) {
    set_nv_tkud(rd, index, 0, 0, -2, "RHO_DOT", "[1]", "RHO_DOT", FALSE);
    index++;
    if (RHO_DOT == 2) {
      Export_XP_ID[index_post_export] = index_post;
      index_post_export++;
    }
    RHO_DOT = index_post;
    index_post++;
  }

  if (MEAN_SHEAR != -1 && Num_Var_In_Type[pg->imtrx][R_MOMENTUM1]) {
    set_nv_tkud(rd, index, 0, 0, -2, "SHEAR", "[1]", "Mean shear rate", FALSE);
    index++;
    if (MEAN_SHEAR == 2) {
      Export_XP_ID[index_post_export] = index_post;
      index_post_export++;
    }
    MEAN_SHEAR = index_post;
    index_post++;
  }

  if (GIES_CRIT != -1 && Num_Var_In_Type[pg->imtrx][R_MOMENTUM1]) {
    set_nv_tkud(rd, index, 0, 0, -2, "GIES", "[1]", "Giesekus Criterion", FALSE);
    index++;
    if (GIES_CRIT == 2) {
      Export_XP_ID[index_post_export] = index_post;
      index_post_export++;
    }
    GIES_CRIT = index_post;
    index_post++;
  }

  if (Q_FCN != -1 && Num_Var_In_Type[pg->imtrx][R_MOMENTUM1]) {
    set_nv_tkud(rd, index, 0, 0, -2, "Q_FCN", "[1]", "Q Function", FALSE);
    index++;
    if (Q_FCN == 2) {
      Export_XP_ID[index_post_export] = index_post;
      index_post_export++;
    }
    Q_FCN = index_post;
    index_post++;
  }

  if (PSPG_PP != -1 && Num_Var_In_Type[pg->imtrx][R_MOMENTUM1]) {
    PSPG_PP = index_post;
    set_nv_tkud(rd, index, 0, 0, -2, "PSPG_X", "[1]", "PSPG PP X", FALSE);
    index++;
    index_post++;
    set_nv_tkud(rd, index, 0, 0, -2, "PSPG_Y", "[1]", "PSPG PP Y", FALSE);
    index++;
    index_post++;
    if (VIM > 2) {
      set_nv_tkud(rd, index, 0, 0, -2, "PSPG_Z", "[1]", "PSPG PP Z", FALSE);
      index++;
      index_post++;
    }
  }

  if (VELO_SPEED != -1 && Num_Var_In_Type[pg->imtrx][R_MOMENTUM1]) {
    set_nv_tkud(rd, index, 0, 0, -2, "VELO_SPEED", "[1]", "Velocity Magnitude", FALSE);
    index++;
    if (VELO_SPEED == 2) {
      Export_XP_ID[index_post_export] = index_post;
      index_post_export++;
    }
    VELO_SPEED = index_post;
    index_post++;
  }

  check = 0;
  for (i = 0; i < upd->Num_Mat; i++) {
    if (pd_glob[i]->MeshMotion == LAGRANGIAN || pd_glob[i]->MeshMotion == DYNAMIC_LAGRANGIAN)
      check = 1;
  }
  if (PRESSURE_CONT != -1 && (Num_Var_In_Type[pg->imtrx][R_MOMENTUM1] || check)) {
    set_nv_tkud(rd, index, 0, 0, -2, "PRESSURE", "[1]", "hydrodynamic pressure", FALSE);
    index++;
    if (PRESSURE_CONT == 2) {
      Export_XP_ID[index_post_export] = index_post;
      index_post_export++;
    }
    PRESSURE_CONT = index_post;
    index_post++;
  }

  if (SH_DIV_S_V_CONT != -1 && (Num_Var_In_Type[pg->imtrx][R_MOMENTUM1])) {
    set_nv_tkud(rd, index, 0, 0, -2, "SH_DIV_S_V", "[1]", "shell div_s_v", FALSE);
    index++;
    if (SH_DIV_S_V_CONT == 2) {
      Export_XP_ID[index_post_export] = index_post;
      index_post_export++;
    }
    SH_DIV_S_V_CONT = index_post;
    index_post++;
  }

  if (SH_CURV_CONT != -1 && (Num_Var_In_Type[pg->imtrx][R_SHELL_SURF_CURV])) {
    set_nv_tkud(rd, index, 0, 0, -2, "SH_SURF_CURV", "[1]", "shell curv", FALSE);
    index++;
    if (SH_CURV_CONT == 2) {
      Export_XP_ID[index_post_export] = index_post;
      index_post_export++;
    }
    SH_CURV_CONT = index_post;
    index_post++;
  }

  if (FILL_CONT != -1 && (Num_Var_In_Type[pg->imtrx][R_FILL])) {
    set_nv_tkud(rd, index, 0, 0, -2, "FILL", "[1]", "continuous fill", FALSE);
    index++;
    if (FILL_CONT == 2) {
      Export_XP_ID[index_post_export] = index_post;
      index_post_export++;
    }
    FILL_CONT = index_post;
    index_post++;
  }

  if (CONC_CONT != -1 && (Num_Var_In_Type[pg->imtrx][R_MASS])) {
    if (CONC_CONT == 2) {
      GOMA_EH(GOMA_ERROR, "Post-processing vectors cannot be exported yet!");
    }
    CONC_CONT = index_post;
    for (w = 0; w < upd->Max_Num_Species_Eqn; w++) {
      sprintf(species_name, "YC%d", w);
      sprintf(species_desc, "Concentration of %d", w);
      set_nv_tkud(rd, index, 0, 0, -2, species_name, "[1]", species_desc, FALSE);
      index++;
      index_post++;
    }
  }

  if (STRESS_CONT != -1 && (Num_Var_In_Type[pg->imtrx][POLYMER_STRESS11])) {
    int index_post_save = index_post;

    if (STRESS_CONT == 2) {
      GOMA_EH(GOMA_ERROR, "Post-processing vectors cannot be exported yet!");
    }

    for (mode = 0; mode < MAX_MODES; mode++) {
      for (a = 0; a < VIM; a++) {
        for (b = 0; b < VIM; b++) {
          /* stress tensor is symmetric */
          if (a <= b) {
            if (Num_Var_In_Type[pg->imtrx][v_s[mode][a][b]]) {
              sprintf(species_name, "cs%d%d_%d", a + 1, b + 1, mode);
              sprintf(species_desc, "Continuous stress %d%d_%d", a + 1, b + 1, mode);
              set_nv_tkud(rd, index, 0, 0, -2, species_name, "[1]", species_desc, FALSE);
              index++;
              if (STRESS_CONT == 2) {
                Export_XP_ID[index_post_export] = index_post;
                index_post_export++;
              }
              index_post++;
            }
          }
        }
      }
    }

    STRESS_CONT = index_post_save;

    check = 0;
    for (i = 0; i < upd->Num_Mat; i++) {
      if (vn_glob[i]->modes > 1)
        check = 1;
    }

    /* write out total stress if we have more than one mode */
    if (check) {
      for (a = 0; a < VIM; a++) {
        for (b = 0; b < VIM; b++) {
          /* stress tensor is symmetric */
          if (a <= b) {
            if (Num_Var_In_Type[pg->imtrx][v_s[0][a][b]]) {
              sprintf(species_name, "tcs%d%d", a + 1, b + 1);
              sprintf(species_desc, "Total continuous stress %d%d", a + 1, b + 1);
              set_nv_tkud(rd, index, 0, 0, -2, species_name, "[1]", species_desc, FALSE);
              index++;
              index_post++;
            }
          }
        }
      }
    }
  }

  if (EM_CONTOURS != -1 && Num_Var_In_Type[pg->imtrx][EM_E1_REAL]) {
    EM_CONTOURS = index_post;
    int dim = DIM;
    if (pd->gv[EM_E1_REAL]) {
      sprintf(species_name, "EM_REALX");
      sprintf(species_desc, "EM X Vector");
      set_nv_tkud(rd, index, 0, 0, -2, species_name, "[1]", species_desc, FALSE);
      index++;
      index_post++;

      if (dim > 1) {
        sprintf(species_name, "EM_REALY");
        sprintf(species_desc, "EM Y Vector");
        set_nv_tkud(rd, index, 0, 0, -2, species_name, "[1]", species_desc, FALSE);
        index++;
        index_post++;
      }

      if (dim > 2) {
        sprintf(species_name, "EM_REALZ");
        sprintf(species_desc, "EM Z Vector");
        set_nv_tkud(rd, index, 0, 0, -2, species_name, "[1]", species_desc, FALSE);
        index++;
        index_post++;
      }
    }
    if (pd->gv[EM_E1_IMAG]) {
      sprintf(species_name, "EM_IMAGX");
      sprintf(species_desc, "EM X Vector");
      set_nv_tkud(rd, index, 0, 0, -2, species_name, "[1]", species_desc, FALSE);
      index++;
      index_post++;

      if (dim > 1) {
        sprintf(species_name, "EM_IMAGY");
        sprintf(species_desc, "EM Y Vector");
        set_nv_tkud(rd, index, 0, 0, -2, species_name, "[1]", species_desc, FALSE);
        index++;
        index_post++;
      }

      if (dim > 2) {
        sprintf(species_name, "EM_IMAGZ");
        sprintf(species_desc, "EM Z Vector");
        set_nv_tkud(rd, index, 0, 0, -2, species_name, "[1]", species_desc, FALSE);
        index++;
        index_post++;
      }
    }
  }
  if (TOTAL_EM_CONTOURS != -1 && Num_Var_In_Type[pg->imtrx][EM_E1_REAL]) {
    TOTAL_EM_CONTOURS = index_post;
    int dim = DIM;
    if (pd->gv[EM_E1_REAL]) {
      sprintf(species_name, "EM_TOT_REALX");
      sprintf(species_desc, "EM X Vector");
      set_nv_tkud(rd, index, 0, 0, -2, species_name, "[1]", species_desc, FALSE);
      index++;
      index_post++;

      if (dim > 1) {
        sprintf(species_name, "EM_TOT_REALY");
        sprintf(species_desc, "EM Y Vector");
        set_nv_tkud(rd, index, 0, 0, -2, species_name, "[1]", species_desc, FALSE);
        index++;
        index_post++;
      }

      if (dim > 2) {
        sprintf(species_name, "EM_TOT_REALZ");
        sprintf(species_desc, "EM Z Vector");
        set_nv_tkud(rd, index, 0, 0, -2, species_name, "[1]", species_desc, FALSE);
        index++;
        index_post++;
      }
    }
    if (pd->gv[EM_E1_IMAG]) {
      sprintf(species_name, "EM_TOT_IMAGX");
      sprintf(species_desc, "EM X Vector");
      set_nv_tkud(rd, index, 0, 0, -2, species_name, "[1]", species_desc, FALSE);
      index++;
      index_post++;

      if (dim > 1) {
        sprintf(species_name, "EM_TOT_IMAGY");
        sprintf(species_desc, "EM Y Vector");
        set_nv_tkud(rd, index, 0, 0, -2, species_name, "[1]", species_desc, FALSE);
        index++;
        index_post++;
      }

      if (dim > 2) {
        sprintf(species_name, "EM_TOT_IMAGZ");
        sprintf(species_desc, "EM Z Vector");
        set_nv_tkud(rd, index, 0, 0, -2, species_name, "[1]", species_desc, FALSE);
        index++;
        index_post++;
      }
    }
  }
  if (SCATTERED_EM_CONTOURS != -1 && Num_Var_In_Type[pg->imtrx][EM_E1_REAL]) {
    SCATTERED_EM_CONTOURS = index_post;
    int dim = DIM;
    if (pd->gv[EM_E1_REAL]) {
      sprintf(species_name, "EM_SCA_REALX");
      sprintf(species_desc, "EM X Vector");
      set_nv_tkud(rd, index, 0, 0, -2, species_name, "[1]", species_desc, FALSE);
      index++;
      index_post++;

      if (dim > 1) {
        sprintf(species_name, "EM_SCA_REALY");
        sprintf(species_desc, "EM Y Vector");
        set_nv_tkud(rd, index, 0, 0, -2, species_name, "[1]", species_desc, FALSE);
        index++;
        index_post++;
      }

      if (dim > 2) {
        sprintf(species_name, "EM_SCA_REALZ");
        sprintf(species_desc, "EM Z Vector");
        set_nv_tkud(rd, index, 0, 0, -2, species_name, "[1]", species_desc, FALSE);
        index++;
        index_post++;
      }
    }
    if (pd->gv[EM_E1_IMAG]) {
      sprintf(species_name, "EM_SCA_IMAGX");
      sprintf(species_desc, "EM X Vector");
      set_nv_tkud(rd, index, 0, 0, -2, species_name, "[1]", species_desc, FALSE);
      index++;
      index_post++;

      if (dim > 1) {
        sprintf(species_name, "EM_SCA_IMAGY");
        sprintf(species_desc, "EM Y Vector");
        set_nv_tkud(rd, index, 0, 0, -2, species_name, "[1]", species_desc, FALSE);
        index++;
        index_post++;
      }

      if (dim > 2) {
        sprintf(species_name, "EM_SCA_IMAGZ");
        sprintf(species_desc, "EM Z Vector");
        set_nv_tkud(rd, index, 0, 0, -2, species_name, "[1]", species_desc, FALSE);
        index++;
        index_post++;
      }
    }
  }

  if (MOMENT_SOURCES != -1 && Num_Var_In_Type[pg->imtrx][MOMENT0]) {
    MOMENT_SOURCES = index_post;
    for (int mom = 0; mom < MAX_MOMENTS; mom++) {
      if (pd->gv[MOMENT0 + mom]) {
        sprintf(species_name, "MOMSOURCE%d", mom);
        sprintf(species_desc, "Moment Source %d", mom);
        set_nv_tkud(rd, index, 0, 0, -2, species_name, "[1]", species_desc, FALSE);
        index++;
        index_post++;
      }
    }
  }

  if (YZBETA != -1 && Num_Var_In_Type[pg->imtrx][MASS_FRACTION]) {
    YZBETA = index_post;
    for (int w = 0; w < pd->Num_Species_Eqn; w++) {
      sprintf(species_name, "YZBETA_%d", w);
      sprintf(species_desc, "YZbeta %d", w);
      set_nv_tkud(rd, index, 0, 0, -2, species_name, "[1]", species_desc, FALSE);
      index++;
      index_post++;
    }
  }

  if (FIRST_INVAR_STRAIN != -1 && Num_Var_In_Type[pg->imtrx][R_MESH1]) {
    set_nv_tkud(rd, index, 0, 0, -2, "IE", "[1]", "1st invariant of strain tensor", FALSE);
    index++;
    if (FIRST_INVAR_STRAIN == 2) {
      Export_XP_ID[index_post_export] = index_post;
      index_post_export++;
    }
    FIRST_INVAR_STRAIN = index_post;
    index_post++;
  }

  if (SEC_INVAR_STRAIN != -1 && Num_Var_In_Type[pg->imtrx][R_MESH1]) {
    set_nv_tkud(rd, index, 0, 0, -2, "IIE", "[1]", "2nd invariant of strain tensor", FALSE);
    index++;
    if (SEC_INVAR_STRAIN == 2) {
      Export_XP_ID[index_post_export] = index_post;
      index_post_export++;
    }
    SEC_INVAR_STRAIN = index_post;
    index_post++;
  }

  if (THIRD_INVAR_STRAIN != -1 && Num_Var_In_Type[pg->imtrx][R_MESH1]) {
    set_nv_tkud(rd, index, 0, 0, -2, "IIIE", "[1]", "3rd invariant of strain tensor", FALSE);
    index++;
    if (THIRD_INVAR_STRAIN == 2) {
      Export_XP_ID[index_post_export] = index_post;
      index_post_export++;
    }
    THIRD_INVAR_STRAIN = index_post;
    index_post++;
  }

  if (FIRST_STRAINRATE_INVAR != -1 && Num_Var_In_Type[pg->imtrx][R_MESH1]) {
    set_nv_tkud(rd, index, 0, 0, -2, "ISR", "[1]", "1st strain rate invariant", FALSE);
    index++;
    if (FIRST_STRAINRATE_INVAR == 2) {
      Export_XP_ID[index_post_export] = index_post;
      index_post_export++;
    }
    FIRST_STRAINRATE_INVAR = index_post;
    index_post++;
  }

  if (SEC_STRAINRATE_INVAR != -1 && Num_Var_In_Type[pg->imtrx][R_MESH1]) {
    set_nv_tkud(rd, index, 0, 0, -2, "IISR", "[1]", "2nd strain rate invariant", FALSE);
    index++;
    if (SEC_STRAINRATE_INVAR == 2) {
      Export_XP_ID[index_post_export] = index_post;
      index_post_export++;
    }
    SEC_STRAINRATE_INVAR = index_post;
    index_post++;
  }

  if (THIRD_STRAINRATE_INVAR != -1 && Num_Var_In_Type[pg->imtrx][R_MESH1]) {
    set_nv_tkud(rd, index, 0, 0, -2, "IIISR", "[1]", "3rd strain rate invariant", FALSE);
    index++;
    if (THIRD_STRAINRATE_INVAR == 2) {
      Export_XP_ID[index_post_export] = index_post;
      index_post_export++;
    }
    THIRD_STRAINRATE_INVAR = index_post;
    index_post++;
  }

  if (DIELECTROPHORETIC_FIELD != -1 && Num_Var_In_Type[pg->imtrx][R_POTENTIAL]) {
    if (DIELECTROPHORETIC_FIELD == 2) {
      GOMA_EH(GOMA_ERROR, "Post-processing vectors cannot be exported yet!");
    }
    DIELECTROPHORETIC_FIELD = index_post;
    for (i = 0; i < Num_Dim; i++) {
      sprintf(nm, "DF%d", i);
      sprintf(ds, "Dielectrophoretic force field component %d", i);
      set_nv_tkud(rd, index, 0, 0, -2, nm, "[1]", ds, FALSE);
      index++;
      if (DIELECTROPHORETIC_FIELD == 2) {
        Export_XP_ID[index_post_export] = index_post;
        index_post_export++;
      }
      index_post++;
    }
  }

  if (DIELECTROPHORETIC_FIELD_NORM != -1 && Num_Var_In_Type[pg->imtrx][R_POTENTIAL]) {
    set_nv_tkud(rd, index, 0, 0, -2, "DFN", "[1]", "Dielectrophoretic force norm", FALSE);
    index++;
    if (DIELECTROPHORETIC_FIELD_NORM == 2) {
      Export_XP_ID[index_post_export] = index_post;
      index_post_export++;
    }
    DIELECTROPHORETIC_FIELD_NORM = index_post;
    index_post++;
  }

  if (ENORMSQ_FIELD != -1 && Num_Var_In_Type[pg->imtrx][R_POTENTIAL] &&
      Num_Var_In_Type[pg->imtrx][R_ENORM]) {
    if (ENORMSQ_FIELD == 2) {
      GOMA_EH(GOMA_ERROR, "Post-processing vectors cannot be exported yet!");
    }
    ENORMSQ_FIELD = index_post;
    for (i = 0; i < Num_Dim; i++) {
      sprintf(nm, "GENS%d", i);
      sprintf(ds, "grad(|E|^2) field component %d", i);
      set_nv_tkud(rd, index, 0, 0, -2, nm, "[1]", ds, FALSE);
      index++;
      if (ENORMSQ_FIELD == 2) {
        Export_XP_ID[index_post_export] = index_post;
        index_post_export++;
      }
      index_post++;
    }
  }

  if (ENORMSQ_FIELD_NORM != -1 && Num_Var_In_Type[pg->imtrx][R_POTENTIAL] &&
      Num_Var_In_Type[pg->imtrx][R_ENORM]) {
    set_nv_tkud(rd, index, 0, 0, -2, "GENSNORM", "[1]", "|grad(|E|^2)|", FALSE);
    index++;
    if (ENORMSQ_FIELD_NORM == 2) {
      Export_XP_ID[index_post_export] = index_post;
      index_post_export++;
    }
    ENORMSQ_FIELD_NORM = index_post;
    index_post++;
  }

  if (DIFFUSION_VECTORS != -1) {
    char dir_ch[] = "XYZ";
    if (DIFFUSION_VECTORS == 2) {
      GOMA_EH(GOMA_ERROR, "Post-processing vectors cannot be exported yet!");
    }
    DIFFUSION_VECTORS = index_post;
    if (upd->Max_Num_Species_Eqn == 0)
      DIFFUSION_VECTORS = -1;
    for (w = 0; w < upd->Max_Num_Species_Eqn; w++) {
      for (i = 0; i < Num_Dim; i++) {
        sprintf(species_name, "Y%dDIF%c", w, dir_ch[i]);
        sprintf(species_desc, "Diffusion of %d in %c direction", w, dir_ch[i]);
        set_nv_tkud(rd, index, 0, 0, -2, species_name, "[1]", species_desc, FALSE);
        index++;
        index_post++;
      }
    }
    DIFFUSION_VECTORS_POR_LIQ_GPHASE = 1;
    DIFFUSION_VECTORS_POR_AIR_GPHASE = 1;
  }
  /*
   * Check for presence of gas phase diffusion of liquid
   * and air in a porous flow problem
   */
  check = 0;
  for (mn = 0; mn < upd->Num_Mat; mn++) {
    if (mp_glob[mn]->PorousMediaType == POROUS_UNSATURATED ||
        mp_glob[mn]->PorousMediaType == POROUS_TWO_PHASE) {
      check = 1;
    }
  }
  if (!check) {
    DIFFUSION_VECTORS_POR_LIQ_GPHASE = -1;
    DIFFUSION_VECTORS_POR_AIR_GPHASE = -1;
  }

  if (DIFFUSION_VECTORS_POR_LIQ_GPHASE != -1) {
    DIFFUSION_VECTORS_POR_LIQ_GPHASE = index_post;
    for (i = 0; i < Num_Dim; i++) {
      sprintf(species_name, "YDIF_liq_gphase_%d", i);
      sprintf(species_desc, "Diffusion of liq in gas phase in %d direction", i);
      set_nv_tkud(rd, index, 0, 0, -2, species_name, "[1]", species_desc, FALSE);
      index++;
      index_post++;
    }
  }
  if (DIFFUSION_VECTORS_POR_AIR_GPHASE != -1) {
    DIFFUSION_VECTORS_POR_AIR_GPHASE = index_post;
    for (i = 0; i < Num_Dim; i++) {
      sprintf(species_name, "YDIF_air_gphase_%d", i);
      sprintf(species_desc, "Diffusion of air in gas phase in %d direction", i);
      set_nv_tkud(rd, index, 0, 0, -2, species_name, "[1]", species_desc, FALSE);
      index++;
      index_post++;
    }
  }

  if (FLUXLINES != -1 && Num_Var_In_Type[pg->imtrx][R_MASS]) {
    if (FLUXLINES == 2) {
      GOMA_EH(GOMA_ERROR, "Post-processing vectors cannot be exported yet!");
    }
    if (Num_Dim == 3) {
      GOMA_WH(-1, "Cant do flux function in 3D");
      FLUXLINES = -1;
    } else if (TimeIntegration != STEADY) {
      GOMA_WH(-1, "Cant do flux function in transient");
      FLUXLINES = -1;
    } else {
      FLUXLINES = index_post;
      for (w = 0; w < upd->Max_Num_Species_Eqn; w++) {
        sprintf(species_name, "Y%dFLUX", w);
        sprintf(species_desc, "Fluxfunction of %d", w);
        set_nv_tkud(rd, index, 0, 0, -2, species_name, "[1]", species_desc, FALSE);
        index++;
        index_post++;
      }
    }
  }

  if (CONDUCTION_VECTORS != -1 && Num_Var_In_Type[pg->imtrx][R_ENERGY]) {
    if (CONDUCTION_VECTORS == 2) {
      GOMA_EH(GOMA_ERROR, "Post-processing vectors cannot be exported yet!");
    }
    CONDUCTION_VECTORS = index_post;
    /* X Component */
    sprintf(nm, "TCONDX");
    sprintf(ds, "Conduction in X direction");
    set_nv_tkud(rd, index, 0, 0, -2, nm, "[1]", ds, FALSE);
    index++;
    index_post++;
    /* Y Component */
    sprintf(nm, "TCONDY");
    sprintf(ds, "Conduction in Y direction");
    set_nv_tkud(rd, index, 0, 0, -2, nm, "[1]", ds, FALSE);
    index++;
    index_post++;
    /* Z Component */
    if (Num_Dim > 2) {
      sprintf(nm, "TCONDZ");
      sprintf(ds, "Conduction in Z direction");
      set_nv_tkud(rd, index, 0, 0, -2, nm, "[1]", ds, FALSE);
      index++;
      index_post++;
    }
  }

  if (POYNTING_VECTORS != -1 &&
      ((Num_Var_In_Type[pg->imtrx][R_ACOUS_PREAL] || Num_Var_In_Type[pg->imtrx][R_ACOUS_PIMAG]) ||
       (Num_Var_In_Type[pg->imtrx][R_EM_E1_REAL] || Num_Var_In_Type[pg->imtrx][R_EM_E2_REAL] ||
        Num_Var_In_Type[pg->imtrx][R_EM_E3_REAL]))) {
    if (POYNTING_VECTORS == 2) {
      GOMA_EH(GOMA_ERROR, "Post-processing vectors cannot be exported yet!");
    }
    POYNTING_VECTORS = index_post;
    /* X Component */
    sprintf(nm, "POYNTX");
    sprintf(ds, "Poynting in X direction");
    set_nv_tkud(rd, index, 0, 0, -2, nm, "[1]", ds, FALSE);
    index++;
    index_post++;
    /* Y Component */
    sprintf(nm, "POYNTY");
    sprintf(ds, "Poynting in Y direction");
    set_nv_tkud(rd, index, 0, 0, -2, nm, "[1]", ds, FALSE);
    index++;
    index_post++;
    /* Z Component */
    if (Num_Dim > 2) {
      sprintf(nm, "POYNTZ");
      sprintf(ds, "Poynting in Z direction");
      set_nv_tkud(rd, index, 0, 0, -2, nm, "[1]", ds, FALSE);
      index++;
      index_post++;
    }
  }

  if (ORIENTATION_VECTORS != -1 && Num_Var_In_Type[pg->imtrx][R_ENERGY]) {
    if (ORIENTATION_VECTORS == 2) {
      GOMA_EH(GOMA_ERROR, "Post-processing vectors cannot be exported yet!");
    }
    ORIENTATION_VECTORS = index_post;
    /* X Component */
    sprintf(nm, "ORX");
    sprintf(ds, "Orientation in X direction");
    set_nv_tkud(rd, index, 0, 0, -2, nm, "[1]", ds, FALSE);
    index++;
    index_post++;
    /* Y Component */
    sprintf(nm, "ORY");
    sprintf(ds, "Orientation in Y direction");
    set_nv_tkud(rd, index, 0, 0, -2, nm, "[1]", ds, FALSE);
    index++;
    index_post++;
    /* Z Component */
    if (Num_Dim > 2) {
      sprintf(nm, "ORZ");
      sprintf(ds, "Orientation in Z direction");
      set_nv_tkud(rd, index, 0, 0, -2, nm, "[1]", ds, FALSE);
      index++;
      index_post++;
    }
  }

  if (SPECIES_SOURCES != -1 && Num_Var_In_Type[pg->imtrx][R_MASS]) {
    if (SPECIES_SOURCES == 2) {
      GOMA_EH(GOMA_ERROR, "Post-processing species source vectors cannot be exported yet!");
    }
    SPECIES_SOURCES = index_post;
    if (upd->Max_Num_Species_Eqn == 0)
      SPECIES_SOURCES = -1;
    for (w = 0; w < upd->Max_Num_Species_Eqn; w++) {
      sprintf(species_name, "YSRC%d", w);
      sprintf(species_desc, "Source of %d ", w);
      set_nv_tkud(rd, index, 0, 0, -2, species_name, "[1]", species_desc, FALSE);
      index++;
      index_post++;
    }
  }

  if (ELECTRIC_FIELD != -1 && Num_Var_In_Type[pg->imtrx][R_POTENTIAL]) {
    if (ELECTRIC_FIELD == 2) {
      GOMA_EH(GOMA_ERROR, "Post-processing vectors cannot be exported yet!");
    }
    ELECTRIC_FIELD = index_post;

    /* X Component */
    sprintf(nm, "EX");
    sprintf(ds, "X-component of the electric field");
    set_nv_tkud(rd, index, 0, 0, -2, nm, "[1]", ds, FALSE);
    index++;
    index_post++;
    /* Y Component */
    sprintf(nm, "EY");
    sprintf(ds, "Y-component of the electric field");
    set_nv_tkud(rd, index, 0, 0, -2, nm, "[1]", ds, FALSE);
    index++;
    index_post++;
    /* Z Component */
    if (Num_Dim > 2) {
      sprintf(nm, "EZ");
      sprintf(ds, "Z-component of the electric field");
      set_nv_tkud(rd, index, 0, 0, -2, nm, "[1]", ds, FALSE);
      index++;
      index_post++;
    }
  }

  if (ELECTRIC_FIELD_MAG != -1 && Num_Var_In_Type[pg->imtrx][R_POTENTIAL]) {

    sprintf(nm, "EE");
    sprintf(ds, "Electric field magnitude");
    set_nv_tkud(rd, index, 0, 0, -2, nm, "[1]", ds, FALSE);
    index++;
    if (ELECTRIC_FIELD_MAG == 2) {
      Export_XP_ID[index_post_export] = index_post;
      index_post_export++;
    }
    ELECTRIC_FIELD_MAG = index_post;
    index_post++;
  }

  if (ENERGY_FLUXLINES != -1 && Num_Var_In_Type[pg->imtrx][R_ENERGY]) {
    if (Num_Dim == 3) {
      GOMA_WH(-1, "Cant do flux function in 3D");
      ENERGY_FLUXLINES = -1;
    } else if (TimeIntegration != STEADY) {
      GOMA_WH(-1, "Cant do flux function in transient");
      ENERGY_FLUXLINES = -1;
    } else {
      ENERGY_FLUXLINES = index_post;
      set_nv_tkud(rd, index, 0, 0, -2, "EFLUX", "[1]", "Energy Fluxfunction", FALSE);
      index++;
      if (ENERGY_FLUXLINES == 2) {
        Export_XP_ID[index_post_export] = index_post;
        index_post_export++;
      }
      index_post++;
    }
  }

  /* there are always 6 entries in the stress tensor (3-D symmetric) */
  if (STRESS_TENSOR != -1 && Num_Var_In_Type[pg->imtrx][R_MESH1]) {
    if (STRESS_TENSOR == 2) {
      GOMA_EH(GOMA_ERROR, "Post-processing tensors cannot be exported yet!");
    }
    STRESS_TENSOR = index_post;
    set_nv_tkud(rd, index, 0, 0, -2, "T11", "[1]", "Mesh stretch stress x", FALSE);
    index++;
    index_post++;
    set_nv_tkud(rd, index, 0, 0, -2, "T22", "[1]", "Mesh stretch stress y", FALSE);
    index++;
    index_post++;
    set_nv_tkud(rd, index, 0, 0, -2, "T33", "[1]", "Mesh stretch stress z", FALSE);
    index++;
    index_post++;
    set_nv_tkud(rd, index, 0, 0, -2, "T12", "[1]", "Mesh shear stress xy", FALSE);
    index++;
    index_post++;
    if (Num_Dim > 2) {
      set_nv_tkud(rd, index, 0, 0, -2, "T13", "[1]", "Mesh shear stress xz", FALSE);
      index++;
      index_post++;
      set_nv_tkud(rd, index, 0, 0, -2, "T23", "[1]", "Mesh shear stress yz", FALSE);
      index++;
      index_post++;
    }
  }

  /* there are always 6 entries in the stress tensor (3-D symmetric) */
  if (REAL_STRESS_TENSOR != -1 && Num_Var_In_Type[pg->imtrx][R_SOLID1]) {
    if (REAL_STRESS_TENSOR == 2) {
      GOMA_EH(GOMA_ERROR, "Post-processing tensors cannot be exported yet!");
    }
    REAL_STRESS_TENSOR = index_post;
    set_nv_tkud(rd, index, 0, 0, -2, "T11_RS", "[1]", "Real_solid stretch stress x", FALSE);
    index++;
    index_post++;
    set_nv_tkud(rd, index, 0, 0, -2, "T22_RS", "[1]", "Real_solid stretch stress y", FALSE);
    index++;
    index_post++;
    set_nv_tkud(rd, index, 0, 0, -2, "T33_RS", "[1]", "Real_solid stretch stress z", FALSE);
    index++;
    index_post++;
    set_nv_tkud(rd, index, 0, 0, -2, "T12_RS", "[1]", "Real_solid shear stress xy", FALSE);
    index++;
    index_post++;
    if (Num_Dim > 2) {
      set_nv_tkud(rd, index, 0, 0, -2, "T13_RS", "[1]", "Real_solid shear stress xz", FALSE);
      index++;
      index_post++;
      set_nv_tkud(rd, index, 0, 0, -2, "T23_RS", "[1]", "Real_solid shear stress yz", FALSE);
      index++;
      index_post++;
    }
  }

  if (STRAIN_TENSOR != -1 && Num_Var_In_Type[pg->imtrx][R_MESH1]) {
    if (STRAIN_TENSOR == 2) {
      GOMA_EH(GOMA_ERROR, "Post-processing tensors cannot be exported yet!");
    }
    STRAIN_TENSOR = index_post;
    set_nv_tkud(rd, index, 0, 0, -2, "E11", "[1]", "Mesh stetch strain x", FALSE);
    index++;
    index_post++;
    set_nv_tkud(rd, index, 0, 0, -2, "E22", "[1]", "Mesh stetch strain y", FALSE);
    index++;
    index_post++;
    set_nv_tkud(rd, index, 0, 0, -2, "E33", "[1]", "Mesh stetch strain z", FALSE);
    index++;
    index_post++;
    set_nv_tkud(rd, index, 0, 0, -2, "E12", "[1]", "Mesh shear strain xy", FALSE);
    index++;
    index_post++;
    if (Num_Dim > 2) {
      set_nv_tkud(rd, index, 0, 0, -2, "E13", "[1]", "Mesh shear strain xz", FALSE);
      index++;
      index_post++;
      set_nv_tkud(rd, index, 0, 0, -2, "E23", "[1]", "Mesh shear strain yz", FALSE);
      index++;
      index_post++;
    }
  }

  if (evpl_glob[0]->ConstitutiveEquation == EVP_HYPER && EVP_DEF_GRAD_TENSOR != -1 &&
      Num_Var_In_Type[pg->imtrx][R_MESH1]) {
    if (EVP_DEF_GRAD_TENSOR == 2) {
      GOMA_EH(GOMA_ERROR, "Post-processing tensors cannot be exported yet!");
    }
    EVP_DEF_GRAD_TENSOR = index_post;
    set_nv_tkud(rd, index, 0, 0, -2, "FVP11", "[1]", "Viscoplastic deformation gradient x", FALSE);
    index++;
    index_post++;
    set_nv_tkud(rd, index, 0, 0, -2, "FVP22", "[1]", "Viscoplastic deformation gradient y", FALSE);
    index++;
    index_post++;
    set_nv_tkud(rd, index, 0, 0, -2, "FVP33", "[1]", "Viscoplastic deformation gradient z", FALSE);
    index++;
    index_post++;
    set_nv_tkud(rd, index, 0, 0, -2, "FVP12", "[1]", "Viscoplastic deformation gradient xy", FALSE);
    index++;
    index_post++;
    set_nv_tkud(rd, index, 0, 0, -2, "FVP21", "[1]", "Viscoplastic deformation gradient yx", FALSE);
    index++;
    index_post++;
    if (Num_Dim > 2) {
      set_nv_tkud(rd, index, 0, 0, -2, "FVP13", "[1]", "Viscoplastic deformation gradient xz",
                  FALSE);
      index++;
      index_post++;
      set_nv_tkud(rd, index, 0, 0, -2, "FVP31", "[1]", "Viscoplastic deformation gradient zx",
                  FALSE);
      index++;
      index_post++;
      set_nv_tkud(rd, index, 0, 0, -2, "FVP23", "[1]", "Viscoplastic deformation gradient yz",
                  FALSE);
      index++;
      index_post++;
      set_nv_tkud(rd, index, 0, 0, -2, "FVP32", "[1]", "Viscoplastic deformation gradient zy",
                  FALSE);
      index++;
      index_post++;
    }
    if (STRESS_TENSOR != -1 && Num_Var_In_Type[pg->imtrx][R_MESH1]) {
      set_nv_tkud(rd, index, 0, 0, -2, "TVP11", "[1]", "EVP stretch stress x", FALSE);
      index++;
      index_post++;
      set_nv_tkud(rd, index, 0, 0, -2, "TVP22", "[1]", "EVP stretch stress y", FALSE);
      index++;
      index_post++;
      set_nv_tkud(rd, index, 0, 0, -2, "TVP33", "[1]", "EVP stretch stress z", FALSE);
      index++;
      index_post++;
      set_nv_tkud(rd, index, 0, 0, -2, "TVP12", "[1]", "EVP shear stress xy", FALSE);
      index++;
      index_post++;
      set_nv_tkud(rd, index, 0, 0, -2, "TVP21", "[1]", "EVP shear stress yx", FALSE);
      index++;
      index_post++;
      if (Num_Dim > 2) {
        set_nv_tkud(rd, index, 0, 0, -2, "TVP13", "[1]", "EVP shear stress xz", FALSE);
        index++;
        index_post++;
        set_nv_tkud(rd, index, 0, 0, -2, "TVP31", "[1]", "EVP shear stress zx", FALSE);
        index++;
        index_post++;

        set_nv_tkud(rd, index, 0, 0, -2, "TVP23", "[1]", "EVP shear stress yz", FALSE);
        index++;
        index_post++;
        set_nv_tkud(rd, index, 0, 0, -2, "TVP32", "[1]", "EVP shear stress zy", FALSE);
        index++;
        index_post++;
      }
    }
  }

  check = 0;
  for (i = 0; i < upd->Num_Mat; i++) {
    if (pd_glob[i]->MeshMotion == LAGRANGIAN || pd_glob[i]->MeshMotion == DYNAMIC_LAGRANGIAN)
      check = 1;
  }

  if (LAGRANGE_CONVECTION != -1 && check) {
    if (LAGRANGE_CONVECTION == 2) {
      GOMA_EH(GOMA_ERROR, "Post-processing vectors cannot be exported yet!");
    }
    LAGRANGE_CONVECTION = index_post;
    for (i = 0; i < Num_Dim; i++) {
      sprintf(nm, "VL%d", i);
      sprintf(ds, "Lagrangian Convection in %d direction", i);
      set_nv_tkud(rd, index, 0, 0, -2, nm, "[1]", ds, FALSE);
      index++;
      index_post++;
    }
  }

  if (SURFACE_VECTORS != -1 && Num_Var_In_Type[pg->imtrx][MESH_DISPLACEMENT1]) {
    if (SURFACE_VECTORS == 2) {
      GOMA_EH(GOMA_ERROR, "Post-processing vectors cannot be exported yet!");
    }
    SURFACE_VECTORS = index_post;

    for (i = 0; i < Num_Dim; i++) {
      sprintf(nm, "N%d", i);
      sprintf(ds, "Normal Vector in %d direction", i);
      if (Num_Dim == 2) {
        sprintf(nm1, "T%d", i);
        sprintf(ds1, "Tangent Vector in %d direction", i);
      } else if (Num_Dim == 3) {
        sprintf(nm1, "TA%d", i);
        sprintf(ds1, "1st Tangent Vector in %d direction", i);
        sprintf(nm2, "TB%d", i);
        sprintf(ds2, "2nd Tangent Vector in %d direction", i);
      }

      set_nv_tkud(rd, index, 0, 0, -2, nm, "[1]", ds, FALSE);
      if (Num_Dim == 2) {
        set_nv_tkud(rd, index + Num_Dim, 0, 0, -2, nm1, "[1]", ds1, FALSE);
      } else if (Num_Dim == 3) {
        set_nv_tkud(rd, index + Num_Dim, 0, 0, -2, nm1, "[1]", ds1, FALSE);
        set_nv_tkud(rd, index + Num_Dim + Num_Dim, 0, 0, -2, nm2, "[1]", ds2, FALSE);
      }
      index++;
      index_post++;
    }
    index += Num_Dim * (Num_Dim - 1);
    index_post += Num_Dim * (Num_Dim - 1);
  }

  if (SHELL_NORMALS != -1 &&
      (Num_Var_In_Type[pg->imtrx][SHELL_ANGLE1] || Num_Var_In_Type[pg->imtrx][LUBP])) {
    if (SHELL_NORMALS == 2) {
      GOMA_EH(GOMA_ERROR, "Post-processing vectors cannot be exported yet!");
    }
    SHELL_NORMALS = index_post;

    /* X Component */
    sprintf(nm, "SH_NX");
    sprintf(ds, "X-component of Shell Normal Vector");
    set_nv_tkud(rd, index, 0, 0, -2, nm, "[1]", ds, FALSE);
    index++;
    index_post++;
    /* Y Component */
    sprintf(nm, "SH_NY");
    sprintf(ds, "Y-component of Shell Normal Vector");
    set_nv_tkud(rd, index, 0, 0, -2, nm, "[1]", ds, FALSE);
    index++;
    index_post++;
    /* Z Component */
    if (Num_Dim > 2) {
      sprintf(nm, "SH_NZ");
      sprintf(ds, "Z-component of Shell Normal Vector");
      set_nv_tkud(rd, index, 0, 0, -2, nm, "[1]", ds, FALSE);
      index++;
      index_post++;
    }
  }

  if (CONF_MAP != -1 && Num_Var_In_Type[pg->imtrx][POLYMER_STRESS11]) {
    CONF_MAP = index_post;
    // Loop over any additional viscoelastic modes
    for (mode = 0; mode < MAX_MODES; mode++) {
      for (a = 0; a < VIM; a++) {
        for (b = 0; b < VIM; b++) {
          if (a <= b) {
            if (Num_Var_In_Type[pg->imtrx][v_s[mode][a][b]]) {
              sprintf(species_name, "MS%d%d_%d", a + 1, b + 1, mode);
              sprintf(species_desc, "log conf stress %d%d_%d", a + 1, b + 1, mode);
              set_nv_tkud(rd, index, 0, 0, -2, species_name, "[1]", species_desc, FALSE);
              index++;
              index_post++;
            }
          }
        }
      }
    }

    if (Num_Dim > 2) {
      GOMA_EH(GOMA_ERROR, "Log Conf Stress not implemented for 3D");
    }
  }

  if (STRESS_NORM != -1 && Num_Var_In_Type[pg->imtrx][POLYMER_STRESS11]) {
    STRESS_NORM = index_post;
    int num_modes = 0;
    for (int i = 0; i < MAX_MODES; i++) {
      if (Num_Var_In_Type[pg->imtrx][v_s[i][0][0]]) {
        num_modes++;
      }
    }
    for (int mode = 0; mode < num_modes; mode++) {
      sprintf(species_name, "STRESS_NORM_%d", mode);
      sprintf(species_desc, "stress norm mode_%d", mode);
      set_nv_tkud(rd, index, 0, 0, -2, species_name, "[1]", species_desc, FALSE);
      index++;
      index_post++;
    }
  }

  if (SARAMITO_YIELD != -1 && Num_Var_In_Type[pg->imtrx][POLYMER_STRESS11]) {
    SARAMITO_YIELD = index_post;
    int num_modes = 0;
    for (int i = 0; i < MAX_MODES; i++) {
      if (Num_Var_In_Type[pg->imtrx][v_s[i][0][0]]) {
        num_modes++;
      }
    }
    for (int mode = 0; mode < num_modes; mode++) {
      sprintf(species_name, "SARAMITO_COEFF_%d", mode);
      sprintf(species_desc, "saramito coeff mode_%d", mode);
      set_nv_tkud(rd, index, 0, 0, -2, species_name, "[1]", species_desc, FALSE);
      index++;
      index_post++;
    }
  }

  if (VISCOUS_STRESS != -1 && Num_Var_In_Type[pg->imtrx][R_MOMENTUM1]) {
    if (Num_Dim > 2) {
      VISCOUS_STRESS = index_post;
      set_nv_tkud(rd, index, 0, 0, -2, "VS11", "[1]", "Viscous stress xx", FALSE);
      index++;
      index_post++;
      set_nv_tkud(rd, index, 0, 0, -2, "VS12", "[1]", "Viscous stress xy", FALSE);

      index++;
      index_post++;

      set_nv_tkud(rd, index, 0, 0, -2, "VS13", "[1]", "Viscous stress xz", FALSE);
      index++;
      index_post++;

      set_nv_tkud(rd, index, 0, 0, -2, "VS21", "[1]", "Viscous stress yx", FALSE);

      index++;
      index_post++;

      set_nv_tkud(rd, index, 0, 0, -2, "VS22", "[1]", "Viscous stress yy", FALSE);
      index++;
      index_post++;

      set_nv_tkud(rd, index, 0, 0, -2, "VS23", "[1]", "Viscous stress yz", FALSE);
      index++;
      index_post++;

      set_nv_tkud(rd, index, 0, 0, -2, "VS31", "[1]", "Viscous stress zx", FALSE);

      index++;
      index_post++;

      set_nv_tkud(rd, index, 0, 0, -2, "VS32", "[1]", "Viscous stress zy", FALSE);
      index++;
      index_post++;

      set_nv_tkud(rd, index, 0, 0, -2, "VS33", "[1]", "Viscous stress zz", FALSE);
      index++;
      index_post++;

    } else {
      VISCOUS_STRESS = index_post;
      set_nv_tkud(rd, index, 0, 0, -2, "VS11", "[1]", "Viscous stress xx", FALSE);
      index++;
      index_post++;
      set_nv_tkud(rd, index, 0, 0, -2, "VS12", "[1]", "Viscous stress xy", FALSE);

      index++;
      index_post++;

      set_nv_tkud(rd, index, 0, 0, -2, "VS21", "[1]", "Viscous stress yx", FALSE);

      index++;
      index_post++;

      set_nv_tkud(rd, index, 0, 0, -2, "VS22", "[1]", "Viscous stress yy", FALSE);
      index++;
      index_post++;
    }
  }

  if (VISCOUS_STRESS_NORM != -1 && Num_Var_In_Type[pg->imtrx][R_MOMENTUM1]) {
    SARAMITO_YIELD = -1;
    VISCOUS_STRESS_NORM = index_post;
    set_nv_tkud(rd, index, 0, 0, -2, "VSNORM", "[1]", "Viscous stress norm", FALSE);
    index++;
    index_post++;
  }

  if (VISCOUS_VON_MISES_STRESS != -1 && Num_Var_In_Type[pg->imtrx][R_MOMENTUM1]) {
    VISCOUS_VON_MISES_STRESS = index_post;
    set_nv_tkud(rd, index, 0, 0, -2, "VVMISSTRESS", "[1]", "Viscous Von Mises stress", FALSE);
    index++;
    index_post++;
  }

  /*if (J_FLUX != -1)
  {
    J_FLUX = index_post;
    set_nv_tkud(rd, index, 0, 0, -2, "J1","[1]",
                "particle flux x", FALSE);
    index++;
    index_post++;
    set_nv_tkud(rd, index, 0, 0, -2, "J2","[1]",
                "particle flux y", FALSE);

    index++;
    index_post++;
    set_nv_tkud(rd, index, 0, 0, -2, "J3","[1]",
                "particle flux z", FALSE);
    index++;
    index_post++;
    }*/

  if (GRAD_SH != -1) {
    GRAD_SH = index_post;
    set_nv_tkud(rd, index, 0, 0, -2, "GRAD_SH0", "[1]", "Shear gradient x", FALSE);
    index++;
    index_post++;

    set_nv_tkud(rd, index, 0, 0, -2, "GRAD_SH1", "[1]", "Shear gradient y", FALSE);
    index++;
    index_post++;

    set_nv_tkud(rd, index, 0, 0, -2, "GRAD_SH2", "[1]", "Shear gradient z", FALSE);
    index++;
    index_post++;
  }

  if (GRAD_Y != -1) {
    GRAD_Y = index_post;
    set_nv_tkud(rd, index, 0, 0, -2, "GRAD_Y0", "[1]", "Concentration gradient x", FALSE);
    index++;
    index_post++;

    set_nv_tkud(rd, index, 0, 0, -2, "GRAD_Y1", "[1]", "Concentration gradient y", FALSE);
    index++;
    index_post++;

    set_nv_tkud(rd, index, 0, 0, -2, "GRAD_Y2", "[1]", "Concentration gradient z", FALSE);
    index++;
    index_post++;
  }

  if (EIG != -1) {
    EIG = index_post;
    set_nv_tkud(rd, index, 0, 0, -2, "EIGEN0", "[1]", "First eigenvalue", FALSE);
    index++;
    index_post++;
    set_nv_tkud(rd, index, 0, 0, -2, "EIGEN1", "[1]", "Second eigenvalue", FALSE);

    index++;
    index_post++;
    set_nv_tkud(rd, index, 0, 0, -2, "EIGEN2", "[1]", "Third eigenvalue", FALSE);
    index++;
    index_post++;
  }

  if (EIG1 != -1) {
    EIG1 = index_post;
    set_nv_tkud(rd, index, 0, 0, -2, "EIGEN0_0", "[1]", "First eigenvector x", FALSE);
    index++;
    index_post++;
    set_nv_tkud(rd, index, 0, 0, -2, "EIGEN0_1", "[1]", "First eigenvector y", FALSE);

    index++;
    index_post++;
    set_nv_tkud(rd, index, 0, 0, -2, "EIGEN0_2", "[1]", "First eigenvector z", FALSE);
    index++;
    index_post++;
  }

  if (EIG2 != -1) {
    EIG2 = index_post;
    set_nv_tkud(rd, index, 0, 0, -2, "EIGEN1_0", "[1]", "Second eigenvector x", FALSE);
    index++;
    index_post++;
    set_nv_tkud(rd, index, 0, 0, -2, "EIGEN1_1", "[1]", "Second eigenvector y", FALSE);

    index++;
    index_post++;
    set_nv_tkud(rd, index, 0, 0, -2, "EIGEN1_2", "[1]", "Second eigenvector z", FALSE);
    index++;
    index_post++;
  }

  if (EIG3 != -1) {
    EIG3 = index_post;
    set_nv_tkud(rd, index, 0, 0, -2, "EIGEN2_0", "[1]", "Third eigenvector x", FALSE);
    index++;
    index_post++;
    set_nv_tkud(rd, index, 0, 0, -2, "EIGEN2_1", "[1]", "Third eigenvector y", FALSE);

    index++;
    index_post++;
    set_nv_tkud(rd, index, 0, 0, -2, "EIGEN2_2", "[1]", "Third eigenvector z", FALSE);
    index++;
    index_post++;
  }

  /*
   * Porous flow post-processing setup section
   */
  check = 0;
  for (i = 0; i < upd->Num_Mat; i++) {
    if (mp_glob[i]->PorousMediaType == POROUS_UNSATURATED ||
        mp_glob[i]->PorousMediaType == POROUS_SATURATED ||
        mp_glob[i]->PorousMediaType == POROUS_SHELL_UNSATURATED ||
        mp_glob[i]->PorousMediaType == POROUS_TWO_PHASE) {
      check = 1;
    }
  }

  /*
  printf("Before porous entries, IP = %d and IPE = %d\n",
  index_post, index_post_export);
  */
  if (POROUS_SATURATION != -1 && Num_Var_In_Type[pg->imtrx][R_POR_LIQ_PRES] && check) {
    set_nv_tkud(rd, index, 0, 0, -2, "SAT", "[1]", "Saturation", FALSE);
    index++;
    if (POROUS_SATURATION == 2) {
      Export_XP_ID[index_post_export] = index_post;
      index_post_export++;
    }
    POROUS_SATURATION = index_post;
    index_post++;
  }

  if (POROUS_RHO_TOTAL_SOLVENTS != -1 && Num_Var_In_Type[pg->imtrx][R_POR_LIQ_PRES] && check) {
    if (POROUS_RHO_TOTAL_SOLVENTS == 2) {
      GOMA_EH(GOMA_ERROR, "Post-processing vectors cannot be exported yet!");
    }
    POROUS_RHO_TOTAL_SOLVENTS = index_post;
    sprintf(species_name, "Rho_Total_liq");
    sprintf(species_desc, "Total density of liquid solvent");
    set_nv_tkud(rd, index, 0, 0, -2, species_name, "[1]", species_desc, FALSE);
    index++;
    index_post++;
    sprintf(species_name, "Rho_Total_air");
    sprintf(species_desc, "Total density of gas solvent");
    set_nv_tkud(rd, index, 0, 0, -2, species_name, "[1]", species_desc, FALSE);
    index++;
    index_post++;
    if (Num_Var_In_Type[pg->imtrx][R_POR_POROSITY]) {
      sprintf(species_name, "Rho_Total_solid");
      sprintf(species_desc, "Total density of solid solvent");
      set_nv_tkud(rd, index, 0, 0, -2, species_name, "[1]", species_desc, FALSE);
      index++;
      index_post++;
    }
  }

  if (POROUS_RHO_GAS_SOLVENTS != -1 && Num_Var_In_Type[pg->imtrx][R_POR_LIQ_PRES] && check) {
    if (POROUS_RHO_GAS_SOLVENTS == 2) {
      GOMA_EH(GOMA_ERROR, "Post-processing vectors cannot be exported yet!");
    }
    POROUS_RHO_GAS_SOLVENTS = index_post;
    if (Num_Var_In_Type[pg->imtrx][R_POR_LIQ_PRES]) {
      sprintf(species_name, "RhoSolv_g_liq");
      sprintf(species_desc, "Gas Phase Density of liq solvent");
      set_nv_tkud(rd, index, 0, 0, -2, species_name, "[1]", species_desc, FALSE);
      index++;
      index_post++;
    }
    if (Num_Var_In_Type[pg->imtrx][R_POR_GAS_PRES]) {
      sprintf(species_name, "RhoSolv_g_air");
      sprintf(species_desc, "Gas Phase Density of gas solvent");
      set_nv_tkud(rd, index, 0, 0, -2, species_name, "[1]", species_desc, FALSE);
      index++;
      index_post++;
    }
    if (Num_Var_In_Type[pg->imtrx][R_POR_POROSITY]) {
      sprintf(species_name, "RhoSolv_g_solid");
      sprintf(species_desc, "Gas Phase Density of solid solvent");
      set_nv_tkud(rd, index, 0, 0, -2, species_name, "[1]", species_desc, FALSE);
      index++;
      index_post++;
    }
  }

  if (POROUS_RHO_LPHASE != -1 && Num_Var_In_Type[pg->imtrx][R_POR_LIQ_PRES] && check) {
    sprintf(species_name, "Rho_Liq_Phase");
    sprintf(species_desc, "Density of the Liquid Phase");
    set_nv_tkud(rd, index, 0, 0, -2, species_name, "[1]", species_desc, FALSE);
    index++;
    if (POROUS_RHO_LPHASE == 2) {
      Export_XP_ID[index_post_export] = index_post;
      index_post_export++;
    }
    POROUS_RHO_LPHASE = index_post;
    index_post++;
  }

  if (DARCY_VELOCITY_GAS != -1 && Num_Var_In_Type[pg->imtrx][R_POR_GAS_PRES] && check) {
    if (DARCY_VELOCITY_GAS == 2) {
      GOMA_EH(GOMA_ERROR, "Post-processing vectors cannot be exported yet!");
    }
    DARCY_VELOCITY_GAS = index_post;
    for (i = 0; i < Num_Dim; i++) {
      sprintf(nm, "Darcy_Vel_g_%d", i);
      sprintf(ds, "Gas Phase Darcy Velocity in the %d direction", i);
      set_nv_tkud(rd, index, 0, 0, -2, nm, "[1]", ds, FALSE);
      index++;
      index_post++;
    }
  }

  if (DARCY_VELOCITY_LIQ != -1 && Num_Var_In_Type[pg->imtrx][R_POR_LIQ_PRES] && check) {
    if (DARCY_VELOCITY_LIQ == 2) {
      GOMA_EH(GOMA_ERROR, "Post-processing vectors cannot be exported yet!");
    }
    DARCY_VELOCITY_LIQ = index_post;
    for (i = 0; i < Num_Dim; i++) {
      sprintf(nm, "Darcy_Vel_l_%d", i);
      sprintf(ds, "Liquid Phase Darcy Velocity in the %d direction", i);
      set_nv_tkud(rd, index, 0, 0, -2, nm, "[1]", ds, FALSE);
      index++;
      index_post++;
    }
  }
  if (POROUS_LIQUID_ACCUM_RATE != -1 && Num_Var_In_Type[pg->imtrx][R_POR_LIQ_PRES] && check) {
    set_nv_tkud(rd, index, 0, 0, -2, "Liq_Inv_Dot", "[1]", "Porous Liquid Accumlation Rate", FALSE);
    index++;
    if (POROUS_LIQUID_ACCUM_RATE == 2) {
      Export_XP_ID[index_post_export] = index_post;
      index_post_export++;
    }
    POROUS_LIQUID_ACCUM_RATE = index_post;
    index_post++;
  }

  /*
  printf(" Before cap. pressure, IP = %d and IPE = %d\n",
  index_post, index_post_export);
  */
  if (CAPILLARY_PRESSURE != -1 && Num_Var_In_Type[pg->imtrx][R_POR_LIQ_PRES] && check) {
    set_nv_tkud(rd, index, 0, 0, -2, "PC", "[1]", "Capillary Pressure", FALSE);
    index++;
    if (CAPILLARY_PRESSURE == 2) {
      Export_XP_ID[index_post_export] = index_post;
      index_post_export++;
    }
    CAPILLARY_PRESSURE = index_post;
    index_post++;
  }

  if (POROUS_GRIDPECLET != -1 && Num_Var_In_Type[pg->imtrx][R_POR_LIQ_PRES] && check) {
    set_nv_tkud(rd, index, 0, 0, -2, "Por_Grid_Peclet", "[1]", "Porous Grid Peclet Number", FALSE);
    index++;
    if (POROUS_GRIDPECLET == 2) {
      Export_XP_ID[index_post_export] = index_post;
      index_post_export++;
    }
    POROUS_GRIDPECLET = index_post;
    index_post++;
  }

  if (POROUS_SUPGVELOCITY != -1 && Num_Var_In_Type[pg->imtrx][R_POR_LIQ_PRES] && check) {
    if (POROUS_SUPGVELOCITY == 2) {
      GOMA_EH(GOMA_ERROR, "Post-processing vectors cannot be exported yet!");
    }
    POROUS_SUPGVELOCITY = index_post;
    for (i = 0; i < Num_Dim; i++) {
      sprintf(nm, "U_supg_porous_%d", i);
      sprintf(ds, "Porous SUPG Velocity in the %d direction", i);
      set_nv_tkud(rd, index, 0, 0, -2, nm, "[1]", ds, FALSE);
      index++;
      index_post++;
    }
  }

  if (REL_LIQ_PERM != -1 &&
      (Num_Var_In_Type[pg->imtrx][R_POR_LIQ_PRES] || Num_Var_In_Type[pg->imtrx][R_SHELL_SAT_OPEN] ||
       Num_Var_In_Type[pg->imtrx][R_SHELL_SAT_OPEN_2]) &&
      check) {
    set_nv_tkud(rd, index, 0, 0, -2, "Rel_liq_perm", "[1]", "Relative Liquid Permeability", FALSE);
    index++;
    if (REL_LIQ_PERM == 2) {
      Export_XP_ID[index_post_export] = index_post;
      index_post_export++;
    }
    REL_LIQ_PERM = index_post;
    index_post++;
  } else {
    REL_LIQ_PERM = -1;
  }

  /*
  printf(" After porous entries, IP = %d and IPE = %d\n",
  index_post, index_post_export);
  */

  /* MMH: Output the vorticity vector.
   */
  if (CURL_V != -1 && Num_Var_In_Type[pg->imtrx][R_MOMENTUM1]) {
    if (pd->CoordinateSystem == SWIRLING || pd->CoordinateSystem == PROJECTED_CARTESIAN ||
        pd->CoordinateSystem == CARTESIAN_2pt5D || Num_Dim == 3) {
      if (CURL_V == 2) {
        GOMA_EH(GOMA_ERROR, "Post-processing vectors cannot be exported yet!");
      }
      CURL_V = index_post;
      set_nv_tkud(rd, index, 0, 0, -2, "VORTX", "[1]", "x-component of vorticity", FALSE);
      index++;
      index_post++;
      set_nv_tkud(rd, index, 0, 0, -2, "VORTY", "[1]", "y-component of vorticity", FALSE);
      index++;
      index_post++;
      set_nv_tkud(rd, index, 0, 0, -2, "VORTZ", "[1]", "z-component of vorticity", FALSE);
      index++;
      index_post++;
    } else if (Num_Dim == 2) {
      set_nv_tkud(rd, index, 0, 0, -2, "VORT", "[1]", "Vorticity (scalar)", FALSE);
      index++;
      if (CURL_V == 2) {
        Export_XP_ID[index_post_export] = index_post;
        index_post_export++;
      }
      CURL_V = index_post;
      index_post++;
    } else {
      GOMA_WH(-1, "Why are you asking for vorticity on a 1D problem?");
      CURL_V = -1;
    }
  }

  if (POLYMER_VISCOSITY != -1 && Num_Var_In_Type[pg->imtrx][R_STRESS11]) {
    set_nv_tkud(rd, index, 0, 0, -2, "MUP", "[1]", "Polymer Viscosity", FALSE);
    index++;
    if (POLYMER_VISCOSITY == 2) {
      Export_XP_ID[index_post_export] = index_post;
      index_post_export++;
    }
    POLYMER_VISCOSITY = index_post;
    index_post++;
  } else {
    POLYMER_VISCOSITY = -1;
  }
  /* RBS: Output the Lamb vector.
   */
  if (LAMB_VECTOR != -1 && Num_Var_In_Type[pg->imtrx][R_MOMENTUM1]) {
    if (pd->CoordinateSystem == SWIRLING || pd->CoordinateSystem == PROJECTED_CARTESIAN ||
        pd->CoordinateSystem == CARTESIAN_2pt5D || Num_Dim == 3) {
      if (LAMB_VECTOR == 2) {
        GOMA_EH(GOMA_ERROR, "Post-processing vectors cannot be exported yet!");
      }
      LAMB_VECTOR = index_post;
      set_nv_tkud(rd, index, 0, 0, -2, "LAMBX", "[1]", "x-component of Lamb", FALSE);
      index++;
      index_post++;
      set_nv_tkud(rd, index, 0, 0, -2, "LAMBY", "[1]", "y-component of Lamb", FALSE);
      index++;
      index_post++;
      set_nv_tkud(rd, index, 0, 0, -2, "LAMBZ", "[1]", "z-component of Lamb", FALSE);
      index++;
      index_post++;
    } else if (Num_Dim == 2) {
      LAMB_VECTOR = index_post;
      set_nv_tkud(rd, index, 0, 0, -2, "LAMBX", "[1]", "x-component of Lamb", FALSE);
      index++;
      index_post++;
      set_nv_tkud(rd, index, 0, 0, -2, "LAMBY", "[1]", "y-component of Lamb", FALSE);
      index++;
      index_post++;
      if (LAMB_VECTOR == 2) {
        Export_XP_ID[index_post_export] = index_post;
        index_post_export++;
      }
    }
  }
  /* RBS: Helicity quantity (scalar)
   */
  if (HELICITY != -1 && Num_Var_In_Type[pg->imtrx][R_MOMENTUM1]) {
    if (pd->CoordinateSystem == SWIRLING || pd->CoordinateSystem == PROJECTED_CARTESIAN ||
        pd->CoordinateSystem == CARTESIAN_2pt5D || Num_Dim == 3) {
      if (HELICITY == 2) {
        GOMA_EH(GOMA_ERROR, "Post-processing vectors cannot be exported yet!");
      }
      HELICITY = index_post;
      set_nv_tkud(rd, index, 0, 0, -2, "HELIX", "[1]", "helicity scalar", FALSE);
      index++;
      index_post++;
    } else {
      GOMA_WH(-1, "Why are you asking for Helicity in a 1D or 2D problem?");
      HELICITY = -1;
    }
  }

  if (POLYMER_TIME_CONST != -1 && Num_Var_In_Type[pg->imtrx][R_STRESS11]) {
    set_nv_tkud(rd, index, 0, 0, -2, "LAMBDA", "[1]", "Polymer Time Constant", FALSE);
    index++;
    if (POLYMER_TIME_CONST == 2) {
      Export_XP_ID[index_post_export] = index_post;
      index_post_export++;
    }
    POLYMER_TIME_CONST = index_post;
    index_post++;
  }

  if (MOBILITY_PARAMETER != -1 && Num_Var_In_Type[pg->imtrx][R_STRESS11]) {
    set_nv_tkud(rd, index, 0, 0, -2, "ALPHA", "[1]", "Mobility Parameter", FALSE);
    index++;
    if (MOBILITY_PARAMETER == 2) {
      Export_XP_ID[index_post_export] = index_post;
      index_post_export++;
    }
    MOBILITY_PARAMETER = index_post;
    index_post++;
  }

  if (PTT_XI != -1 && Num_Var_In_Type[pg->imtrx][R_STRESS11]) {
    set_nv_tkud(rd, index, 0, 0, -2, "XI", "[1]", "PTT Xi parameter", FALSE);
    index++;
    if (PTT_XI == 2) {
      Export_XP_ID[index_post_export] = index_post;
      index_post_export++;
    }
    PTT_XI = index_post;
    index_post++;
  }

  if (PTT_EPSILON != -1 && Num_Var_In_Type[pg->imtrx][R_STRESS11]) {
    set_nv_tkud(rd, index, 0, 0, -2, "EPSILON", "[1]", "PTT Epsilon parameter", FALSE);
    index++;
    if (PTT_EPSILON == 2) {
      Export_XP_ID[index_post_export] = index_post;
      index_post_export++;
    }
    PTT_EPSILON = index_post;
    index_post++;
  }

  if (USER_POST != -1) {
    set_nv_tkud(rd, index, 0, 0, -2, "USER", "[1]", "User-Defined", FALSE);
    index++;
    if (USER_POST == 2) {
      Export_XP_ID[index_post_export] = index_post;
      index_post_export++;
    }
    USER_POST = index_post;
    index_post++;
  }

  if (ACOUSTIC_PRESSURE != -1 &&
      (Num_Var_In_Type[pg->imtrx][R_ACOUS_PREAL] && Num_Var_In_Type[pg->imtrx][R_ACOUS_PIMAG])) {
    set_nv_tkud(rd, index, 0, 0, -2, "AC_PRES", "[1]", "Acoustic Pressure Magn.", FALSE);
    index++;
    if (ACOUSTIC_PRESSURE == 2) {
      Export_XP_ID[index_post_export] = index_post;
      index_post_export++;
    }
    ACOUSTIC_PRESSURE = index_post;
    index_post++;
  }

  if (ACOUSTIC_PHASE_ANGLE != -1 &&
      (Num_Var_In_Type[pg->imtrx][R_ACOUS_PREAL] && Num_Var_In_Type[pg->imtrx][R_ACOUS_PIMAG])) {
    set_nv_tkud(rd, index, 0, 0, -2, "AC_PHASE", "[1]", "Acoustic Phase Angle", FALSE);
    index++;
    if (ACOUSTIC_PHASE_ANGLE == 2) {
      Export_XP_ID[index_post_export] = index_post;
      index_post_export++;
    }
    ACOUSTIC_PHASE_ANGLE = index_post;
    index_post++;
  }

  if (ACOUSTIC_ENERGY_DENSITY != -1 &&
      (Num_Var_In_Type[pg->imtrx][R_ACOUS_PREAL] && Num_Var_In_Type[pg->imtrx][R_ACOUS_PIMAG])) {
    set_nv_tkud(rd, index, 0, 0, -2, "AC_ED", "[1]", "Acoustic Energy Density", FALSE);
    index++;
    if (ACOUSTIC_ENERGY_DENSITY == 2) {
      Export_XP_ID[index_post_export] = index_post;
      index_post_export++;
    }
    ACOUSTIC_ENERGY_DENSITY = index_post;
    index_post++;
  }

  if (LIGHT_INTENSITY != -1 &&
      (Num_Var_In_Type[pg->imtrx][R_LIGHT_INTP] && Num_Var_In_Type[pg->imtrx][R_LIGHT_INTM])) {
    set_nv_tkud(rd, index, 0, 0, -2, "L_INT", "[1]", "Light Intensity", FALSE);
    index++;
    if (LIGHT_INTENSITY == 2) {
      Export_XP_ID[index_post_export] = index_post;
      index_post_export++;
    }
    LIGHT_INTENSITY = index_post;
    index_post++;
  }

  if (UNTRACKED_SPEC != -1 && Num_Var_In_Type[pg->imtrx][R_MASS]) {
    set_nv_tkud(rd, index, 0, 0, -2, "UNT_SPEC", "[1]", "Untracked Species", FALSE);
    index++;
    if (UNTRACKED_SPEC == 2) {
      Export_XP_ID[index_post_export] = index_post;
      index_post_export++;
    }
    UNTRACKED_SPEC = index_post;
    index_post++;
  }

  if (PRINCIPAL_STRESS != -1 && Num_Var_In_Type[pg->imtrx][R_MESH1]) {
    if (PRINCIPAL_STRESS == 2) {
      GOMA_EH(GOMA_ERROR, "Post-processing vectors cannot be exported yet!");
    }
    PRINCIPAL_STRESS = index_post;

    /* First Principal Stress */
    sprintf(nm, "PS_1");
    sprintf(ds, "Principal Stress 1");
    set_nv_tkud(rd, index, 0, 0, -2, nm, "[1]", ds, FALSE);
    index++;
    index_post++;
    /* Second Principal Stress */
    sprintf(nm, "PS_2");
    sprintf(ds, "Principal Stress 2");
    set_nv_tkud(rd, index, 0, 0, -2, nm, "[1]", ds, FALSE);
    index++;
    index_post++;
    /* Third Principal Stress */
    sprintf(nm, "PS_3");
    sprintf(ds, "Principal Stress 3");
    set_nv_tkud(rd, index, 0, 0, -2, nm, "[1]", ds, FALSE);
    index++;
    index_post++;
  }

  if (PRINCIPAL_REAL_STRESS != -1 && Num_Var_In_Type[pg->imtrx][R_SOLID1]) {
    if (PRINCIPAL_REAL_STRESS == 2) {
      GOMA_EH(GOMA_ERROR, "Post-processing vectors cannot be exported yet!");
    }
    PRINCIPAL_REAL_STRESS = index_post;

    /* First Principal Stress */
    sprintf(nm, "PS_RS_1");
    sprintf(ds, "Principal Real Stress 1");
    set_nv_tkud(rd, index, 0, 0, -2, nm, "[1]", ds, FALSE);
    index++;
    index_post++;
    /* Second Principal Stress */
    sprintf(nm, "PS_RS_2");
    sprintf(ds, "Principal Real Stress 2");
    set_nv_tkud(rd, index, 0, 0, -2, nm, "[1]", ds, FALSE);
    index++;
    index_post++;
    /* Third Principal Stress */
    sprintf(nm, "PS_RS_3");
    sprintf(ds, "Principal Real Stress 3");
    set_nv_tkud(rd, index, 0, 0, -2, nm, "[1]", ds, FALSE);
    index++;
    index_post++;
  }

  if (LUB_HEIGHT != -1 &&
      (Num_Var_In_Type[pg->imtrx][R_LUBP] || Num_Var_In_Type[pg->imtrx][R_SHELL_FILMP] ||
       Num_Var_In_Type[pg->imtrx][R_TFMP_MASS] || Num_Var_In_Type[pg->imtrx][R_TFMP_BOUND])) {
    if (LUB_HEIGHT == 2) {
      GOMA_EH(GOMA_ERROR, "Post-processing vectors cannot be exported yet!");
    }
    LUB_HEIGHT = index_post;
    sprintf(nm, "LUB_H");
    sprintf(ds, "Lubrication Height");
    set_nv_tkud(rd, index, 0, 0, -2, nm, "[1]", ds, FALSE);
    index++;
    index_post++;
  }

  if (LUB_HEIGHT_2 != -1 && (Num_Var_In_Type[pg->imtrx][R_LUBP_2])) {
    if (LUB_HEIGHT_2 == 2) {
      GOMA_EH(GOMA_ERROR, "Post-processing vectors cannot be exported yet!");
    }
    LUB_HEIGHT_2 = index_post;
    sprintf(nm, "LUB_H_2");
    sprintf(ds, "Lubrication Height 2");
    set_nv_tkud(rd, index, 0, 0, -2, nm, "[1]", ds, FALSE);
    index++;
    index_post++;
  }

  if (LUB_VELO_UPPER != -1 && Num_Var_In_Type[pg->imtrx][R_LUBP]) {
    if (LUB_VELO_UPPER == 2) {
      GOMA_EH(GOMA_ERROR, "Post-processing vectors cannot be exported yet!");
    }
    LUB_VELO_UPPER = index_post;
    sprintf(nm, "LUB_VUX");
    sprintf(ds, "Lubrication Upper Velocity x-component");
    set_nv_tkud(rd, index, 0, 0, -2, nm, "[1]", ds, FALSE);
    index++;
    index_post++;
    sprintf(nm, "LUB_VUY");
    sprintf(ds, "Lubrication Upper Velocity y-component");
    set_nv_tkud(rd, index, 0, 0, -2, nm, "[1]", ds, FALSE);
    index++;
    index_post++;
    sprintf(nm, "LUB_VUZ");
    sprintf(ds, "Lubrication Upper Velocity z-component");
    set_nv_tkud(rd, index, 0, 0, -2, nm, "[1]", ds, FALSE);
    index++;
    index_post++;
  }

  if (LUB_VELO_LOWER != -1 && Num_Var_In_Type[pg->imtrx][R_LUBP]) {
    if (LUB_VELO_LOWER == 2) {
      GOMA_EH(GOMA_ERROR, "Post-processing vectors cannot be exported yet!");
    }
    LUB_VELO_LOWER = index_post;
    sprintf(nm, "LUB_VLX");
    sprintf(ds, "Lubrication Lower Velocity x-component");
    set_nv_tkud(rd, index, 0, 0, -2, nm, "[1]", ds, FALSE);
    index++;
    index_post++;
    sprintf(nm, "LUB_VLY");
    sprintf(ds, "Lubrication Lower Velocity y-component");
    set_nv_tkud(rd, index, 0, 0, -2, nm, "[1]", ds, FALSE);
    index++;
    index_post++;
    sprintf(nm, "LUB_VLZ");
    sprintf(ds, "Lubrication Lower Velocity z-component");
    set_nv_tkud(rd, index, 0, 0, -2, nm, "[1]", ds, FALSE);
    index++;
    index_post++;
  }

  if (LUB_VELO_FIELD != -1 &&
      (Num_Var_In_Type[pg->imtrx][R_LUBP] || Num_Var_In_Type[pg->imtrx][R_SHELL_FILMP])) {
    if (LUB_VELO_FIELD == 2) {
      GOMA_EH(GOMA_ERROR, "Post-processing vectors cannot be exported yet!");
    }
    LUB_VELO_FIELD = index_post;
    sprintf(nm, "LUB_VELO_X");
    sprintf(ds, "Lubrication Velocity x-component");
    set_nv_tkud(rd, index, 0, 0, -2, nm, "[1]", ds, FALSE);
    index++;
    index_post++;
    sprintf(nm, "LUB_VELO_Y");
    sprintf(ds, "Lubrication Velocity y-component");
    set_nv_tkud(rd, index, 0, 0, -2, nm, "[1]", ds, FALSE);
    index++;
    index_post++;
    sprintf(nm, "LUB_VELO_Z");
    sprintf(ds, "Lubrication Velocity z-component");
    set_nv_tkud(rd, index, 0, 0, -2, nm, "[1]", ds, FALSE);
    index++;
    index_post++;
  }

  if (LUB_VELO_FIELD_2 != -1 && (Num_Var_In_Type[pg->imtrx][R_LUBP_2])) {
    if (LUB_VELO_FIELD_2 == 2) {
      GOMA_EH(GOMA_ERROR, "Post-processing vectors cannot be exported yet!");
    }
    LUB_VELO_FIELD_2 = index_post;
    sprintf(nm, "LUB_VELO_2_X");
    sprintf(ds, "Lubrication Velocity x-component");
    set_nv_tkud(rd, index, 0, 0, -2, nm, "[1]", ds, FALSE);
    index++;
    index_post++;
    sprintf(nm, "LUB_VELO_2_Y");
    sprintf(ds, "Lubrication Velocity y-component");
    set_nv_tkud(rd, index, 0, 0, -2, nm, "[1]", ds, FALSE);
    index++;
    index_post++;
    sprintf(nm, "LUB_VELO_2_Z");
    sprintf(ds, "Lubrication Velocity z-component");
    set_nv_tkud(rd, index, 0, 0, -2, nm, "[1]", ds, FALSE);
    index++;
    index_post++;
  }

  if (LUB_FLUID_SOURCE != -1 &&
      ((Num_Var_In_Type[pg->imtrx][R_LUBP]) ||
       (Num_Var_In_Type[pg->imtrx][R_SHELL_FILMP] && Num_Var_In_Type[pg->imtrx][R_SHELL_FILMH]))) {
    if (LUB_FLUID_SOURCE == 2) {
      GOMA_EH(GOMA_ERROR, "Post-processing vectors cannot be exported yet!");
    }
    LUB_FLUID_SOURCE = index_post;
    sprintf(nm, "LUB_SOURCE");
    sprintf(ds, "Lubrication Source");
    set_nv_tkud(rd, index, 0, 0, -2, nm, "[1]", ds, FALSE);
    index++;
    index_post++;
  }

  if (DISJ_PRESS != -1 && Num_Var_In_Type[pg->imtrx][R_SHELL_FILMP]) {
    if (DISJ_PRESS == 2) {
      GOMA_EH(GOMA_ERROR, "Post-processing vectors cannot be exported yet!");
    }
    DISJ_PRESS = index_post;
    sprintf(nm, "DISJ_PRESS");
    sprintf(ds, "Disjoining Pressure");
    set_nv_tkud(rd, index, 0, 0, -2, nm, "[1]", ds, FALSE);
    index++;
    index_post++;
  }

  if (SH_SAT_OPEN != -1 && Num_Var_In_Type[pg->imtrx][R_SHELL_SAT_OPEN]) {
    if (SH_SAT_OPEN == 2) {
      GOMA_EH(GOMA_ERROR, "Post-processing vectors cannot be exported yet!");
    }
    SH_SAT_OPEN = index_post;
    sprintf(nm, "SH_SAT_OPEN");
    sprintf(ds, "Open Porous Shell Saturation");
    set_nv_tkud(rd, index, 0, 0, -2, nm, "[1]", ds, FALSE);
    index++;
    index_post++;
  }

  if (SH_SAT_OPEN_2 != -1 && Num_Var_In_Type[pg->imtrx][R_SHELL_SAT_OPEN_2]) {
    if (SH_SAT_OPEN_2 == 2) {
      GOMA_EH(GOMA_ERROR, "Post-processing vectors cannot be exported yet!");
    }
    SH_SAT_OPEN_2 = index_post;
    sprintf(nm, "SH_SAT_OPEN_2");
    sprintf(ds, "Open Porous Shell Saturation 2");
    set_nv_tkud(rd, index, 0, 0, -2, nm, "[1]", ds, FALSE);
    index++;
    index_post++;
  }

  if (SH_CAP_PRES != -1 && Num_Var_In_Type[pg->imtrx][R_SHELL_SAT_1]) {
    if (SH_CAP_PRES == 2) {
      GOMA_EH(GOMA_ERROR, "Post-processing vectors cannot be exported yet!");
    }
    SH_CAP_PRES = index_post;
    sprintf(nm, "SH_CAP_PRES_1");
    sprintf(ds, "Porous Shell Capillary Pressure Layer 1");
    set_nv_tkud(rd, index, 0, 0, -2, nm, "[1]", ds, FALSE);
    index++;
    index_post++;
    if (Num_Var_In_Type[pg->imtrx][R_SHELL_SAT_2]) {
      sprintf(nm, "SH_CAP_PRES_2");
      sprintf(ds, "Porous Shell Capillary Pressure Layer 2");
      set_nv_tkud(rd, index, 0, 0, -2, nm, "[1]", ds, FALSE);
      index++;
      index_post++;
    }
    if (Num_Var_In_Type[pg->imtrx][R_SHELL_SAT_3]) {
      sprintf(nm, "SH_CAP_PRES_3");
      sprintf(ds, "Porous Shell Capillary Pressure Layer 3");
      set_nv_tkud(rd, index, 0, 0, -2, nm, "[1]", ds, FALSE);
      index++;
      index_post++;
    }
  }

  if (SH_PORE_FLUX != -1 && Num_Var_In_Type[pg->imtrx][R_SHELL_SAT_1] &&
      Num_Var_In_Type[pg->imtrx][R_SHELL_SAT_2]) {
    if (SH_PORE_FLUX == 2) {
      GOMA_EH(GOMA_ERROR, "Post-processing vectors cannot be exported yet!");
    }
    SH_PORE_FLUX = index_post;
    sprintf(nm, "SH_PORE_FLUX_1_2");
    sprintf(ds, "Porous Shell Flux Between Layers 1 and 2");
    set_nv_tkud(rd, index, 0, 0, -2, nm, "[1]", ds, FALSE);
    index++;
    index_post++;

    if (Num_Var_In_Type[pg->imtrx][R_SHELL_SAT_3]) {
      sprintf(nm, "SH_PORE_FLUX_2_3");
      sprintf(ds, "Porous Shell Flux Between Layers 2 and 3");
      set_nv_tkud(rd, index, 0, 0, -2, nm, "[1]", ds, FALSE);
      index++;
      index_post++;
    }
  }

  if (SH_STRESS_TENSOR != -1 &&
      (Num_Var_In_Type[pg->imtrx][R_SHELL_NORMAL1] && Num_Var_In_Type[pg->imtrx][R_SHELL_NORMAL2] &&
       Num_Var_In_Type[pg->imtrx][R_SHELL_NORMAL3] && Num_Var_In_Type[pg->imtrx][R_MESH1])) {
    if (SH_STRESS_TENSOR == 2) {
      GOMA_EH(GOMA_ERROR, "Post-processing vectors cannot be exported yet!");
    }
    SH_STRESS_TENSOR = index_post;
    sprintf(nm, "SH_S11");
    sprintf(ds, "Shell stress tensor component 11");
    set_nv_tkud(rd, index, 0, 0, -2, nm, "[1]", ds, FALSE);
    index++;
    index_post++;
    sprintf(nm, "SH_S22");
    sprintf(ds, "Shell stress tensor component 22");
    set_nv_tkud(rd, index, 0, 0, -2, nm, "[1]", ds, FALSE);
    index++;
    index_post++;
    sprintf(nm, "SH_S12");
    sprintf(ds, "Shell stress tensor component 12");
    set_nv_tkud(rd, index, 0, 0, -2, nm, "[1]", ds, FALSE);
    index++;
    index_post++;
  }

  if (SH_TANG != -1 &&
      (Num_Var_In_Type[pg->imtrx][R_SHELL_NORMAL1] && Num_Var_In_Type[pg->imtrx][R_SHELL_NORMAL2] &&
       Num_Var_In_Type[pg->imtrx][R_SHELL_NORMAL3] && Num_Var_In_Type[pg->imtrx][R_MESH1])) {
    if (SH_TANG == 2) {
      GOMA_EH(GOMA_ERROR, "Post-processing vectors cannot be exported yet!");
    }
    SH_TANG = index_post;

    sprintf(nm, "SH_T1X");
    sprintf(ds, "Shell tangent 1 x - component");
    set_nv_tkud(rd, index, 0, 0, -2, nm, "[1]", ds, FALSE);
    index++;
    index_post++;
    sprintf(nm, "SH_T1Y");
    sprintf(ds, "Shell tangent 1 y - component");
    set_nv_tkud(rd, index, 0, 0, -2, nm, "[1]", ds, FALSE);
    index++;
    index_post++;
    sprintf(nm, "SH_T1Z");
    sprintf(ds, "Shell tangent 1 z - component");
    set_nv_tkud(rd, index, 0, 0, -2, nm, "[1]", ds, FALSE);
    index++;
    index_post++;

    sprintf(nm, "SH_T2X");
    sprintf(ds, "Shell tangent 2 x - component");
    set_nv_tkud(rd, index, 0, 0, -2, nm, "[1]", ds, FALSE);
    index++;
    index_post++;
    sprintf(nm, "SH_T2Y");
    sprintf(ds, "Shell tangent 2 y - component");
    set_nv_tkud(rd, index, 0, 0, -2, nm, "[1]", ds, FALSE);
    index++;
    index_post++;
    sprintf(nm, "SH_T2Z");
    sprintf(ds, "Shell tangent 2 z - component");
    set_nv_tkud(rd, index, 0, 0, -2, nm, "[1]", ds, FALSE);
    index++;
    index_post++;
  }

  if (PP_LAME_MU != -1 && Num_Var_In_Type[pg->imtrx][R_MESH1]) {
    if (PP_LAME_MU == 2) {
      GOMA_EH(GOMA_ERROR, "Post-processing vectors cannot be exported yet!");
    }
    PP_LAME_MU = index_post;
    sprintf(nm, "LAME_MU");
    sprintf(ds, "Lame MU Modulus");
    set_nv_tkud(rd, index, 0, 0, -2, nm, "[1]", ds, FALSE);
    index++;
    index_post++;
  }

  if (PP_LAME_LAMBDA != -1 && Num_Var_In_Type[pg->imtrx][R_MESH1]) {
    if (PP_LAME_LAMBDA == 2) {
      GOMA_EH(GOMA_ERROR, "Post-processing vectors cannot be exported yet!");
    }
    PP_LAME_LAMBDA = index_post;
    sprintf(nm, "LAME_LAMBDA");
    sprintf(ds, "Lame LAMBDA Modulus");
    set_nv_tkud(rd, index, 0, 0, -2, nm, "[1]", ds, FALSE);
    index++;
    index_post++;
  }

  if (VON_MISES_STRAIN != -1 && Num_Var_In_Type[pg->imtrx][R_MESH1]) {
    set_nv_tkud(rd, index, 0, 0, -2, "VME", "[1]", "Von Mises Strain", FALSE);
    index++;
    if (VON_MISES_STRAIN == 2) {
      Export_XP_ID[index_post_export] = index_post;
      index_post_export++;
    }
    VON_MISES_STRAIN = index_post;
    index_post++;
  }

  if (VON_MISES_STRESS != -1 && Num_Var_In_Type[pg->imtrx][R_MESH1]) {
    set_nv_tkud(rd, index, 0, 0, -2, "VMT", "[1]", "Von Mises Stress", FALSE);
    index++;
    if (VON_MISES_STRESS == 2) {
      Export_XP_ID[index_post_export] = index_post;
      index_post_export++;
    }
    VON_MISES_STRESS = index_post;
    index_post++;
  }

  if (TFMP_GAS_VELO != -1 && Num_Var_In_Type[pg->imtrx][R_TFMP_MASS] &&
      Num_Var_In_Type[pg->imtrx][R_TFMP_BOUND]) {
    if (TFMP_GAS_VELO == 2) {
      GOMA_EH(GOMA_ERROR, "Post-processing vectors cannot be exported yet!");
    }
    TFMP_GAS_VELO = index_post;
    sprintf(nm, "TFMP_GAS_VX");
    sprintf(ds, "TFMP Gas Velocity x-component");
    set_nv_tkud(rd, index, 0, 0, -2, nm, "[1]", ds, FALSE);
    index++;
    index_post++;
    sprintf(nm, "TFMP_GAS_VY");
    sprintf(ds, "TFMP Gas Velocity y-component");
    set_nv_tkud(rd, index, 0, 0, -2, nm, "[1]", ds, FALSE);
    index++;
    index_post++;
    sprintf(nm, "TFMP_GAS_VZ");
    sprintf(ds, "TFMP Gas Velocity z-component");
    set_nv_tkud(rd, index, 0, 0, -2, nm, "[1]", ds, FALSE);
    index++;
    index_post++;
  } else {
    TFMP_GAS_VELO = -1;
  }

  if (TFMP_LIQ_VELO != -1 && Num_Var_In_Type[pg->imtrx][R_TFMP_MASS] &&
      Num_Var_In_Type[pg->imtrx][R_TFMP_BOUND]) {
    if (TFMP_LIQ_VELO == 2) {
      GOMA_EH(GOMA_ERROR, "Post-processing vectors cannot be exported yet!");
    }
    TFMP_LIQ_VELO = index_post;
    sprintf(nm, "TFMP_LIQ_VX");
    sprintf(ds, "TFMP Liq Velocity x-component");
    set_nv_tkud(rd, index, 0, 0, -2, nm, "[1]", ds, FALSE);
    index++;
    index_post++;
    sprintf(nm, "TFMP_LIQ_VY");
    sprintf(ds, "TFMP Liq Velocity y-component");
    set_nv_tkud(rd, index, 0, 0, -2, nm, "[1]", ds, FALSE);
    index++;
    index_post++;
    sprintf(nm, "TFMP_LIQ_VZ");
    sprintf(ds, "TFMP Liq Velocity z-component");
    set_nv_tkud(rd, index, 0, 0, -2, nm, "[1]", ds, FALSE);
    index++;
    index_post++;
  } else {
    TFMP_LIQ_VELO = -1;
  }

  if (TFMP_INV_PECLET != -1 && Num_Var_In_Type[pg->imtrx][R_TFMP_MASS] &&
      Num_Var_In_Type[pg->imtrx][R_TFMP_BOUND]) {
    if (TFMP_INV_PECLET == 2) {
      GOMA_EH(GOMA_ERROR, "Post-processing vectors cannot be exported yet!");
    }
    TFMP_INV_PECLET = index_post;
    sprintf(nm, "TFMP_INV_PECLET");
    sprintf(ds, "TFMP Inverted Local Peclet Number");
    set_nv_tkud(rd, index, 0, 0, -2, nm, "[1]", ds, FALSE);
    index++;
    index_post++;
  } else {
    TFMP_INV_PECLET = -1;
  }
  if (TFMP_KRG != -1 && Num_Var_In_Type[pg->imtrx][R_TFMP_BOUND]) {
    if (TFMP_KRG == 2) {
      GOMA_EH(GOMA_ERROR, "Post-processing vectors cannot be exported yet!");
    }
    TFMP_KRG = index_post;
    sprintf(nm, "TFMP_KRG");
    sprintf(ds, "TFMP Local Gas Permeability");
    set_nv_tkud(rd, index, 0, 0, -2, nm, "[1]", ds, FALSE);
    index++;
    index_post++;
  } else {
    TFMP_KRG = -1;
  }

  /* Add external variables if they are present */

  /*
  printf(" Before external vars, IP = %d and IPE = %d\n",
  index_post, index_post_export);
  */
  if (efv->ev) {
    EXTERNAL_POST = index_post;
    for (i = 0; i < efv->Num_external_field; i++) {
      set_nv_tkud(rd, index, 0, 0, -2, efv->name[i], "[1]", "External-fixed field", FALSE);
      index++;
      index_post++;
    }
  }

  /*
   *  Post-processing Step 3: add a new line to put your variable's name into the
   *                          exodus II database in mm_post_proc load_nodal_tkn
   */
  /*
   *   ---- NOTE: any additional post processing variables should be input before
   *              this line, so they are before the NS_RESIDUALS and MM_RESIDUALS
   *              post processing options
   */

  if (NS_RESIDUALS != -1 && Num_Var_In_Type[pg->imtrx][R_MOMENTUM1]) {
    if (NS_RESIDUALS == 2) {
      GOMA_EH(GOMA_ERROR, "Post-processing vectors cannot be exported yet!");
    }
    NS_RESIDUALS = index_post;
    set_nv_tkud(rd, index, 0, 0, -2, "RMX", "[1]", "x-component of NS eqns", FALSE);
    index++;
    index_post++;
    set_nv_tkud(rd, index, 0, 0, -2, "RMY", "[1]", "y-component of NS eqns", FALSE);
    index++;
    index_post++;
    if (Num_Dim > 2) {
      set_nv_tkud(rd, index, 0, 0, -2, "RMZ", "[1]", "z-component of NS eqns", FALSE);
      index++;
      index_post++;
    }
  }

  if (MM_RESIDUALS != -1 && Num_Var_In_Type[pg->imtrx][R_MESH1]) {
    if (MM_RESIDUALS == 2) {
      GOMA_EH(GOMA_ERROR, "Post-processing vectors cannot be exported yet!");
    }
    MM_RESIDUALS = index_post;
    set_nv_tkud(rd, index, 0, 0, -2, "RDX", "[1]", "x-component of mesh eqns", FALSE);
    index++;
    index_post++;
    set_nv_tkud(rd, index, 0, 0, -2, "RDY", "[1]", "y-component of mesh eqns", FALSE);
    index++;
    index_post++;
    if (Num_Dim > 2) {
      set_nv_tkud(rd, index, 0, 0, -2, "RDZ", "[1]", "z-component of mesh eqns", FALSE);
      index++;
      index_post++;
    }
  }

  /* write out total stress */
  var = R_STRESS11;
  if (Num_Var_In_Type[pg->imtrx][var]) {
    post_flag = 0;
    for (i = 0; i < upd->Num_Mat; i++) {
      if (pd_glob[i]->i[pg->imtrx][var] == I_Q1 || pd_glob[i]->i[pg->imtrx][var] == I_Q2 ||
          pd_glob[i]->i[pg->imtrx][var] == I_Q2_D || pd_glob[i]->i[pg->imtrx][var] == I_Q1_D ||
          pd_glob[i]->i[pg->imtrx][var] == I_SP || pd_glob[i]->i[pg->imtrx][var] == I_Q2_LSA ||
          pd_glob[i]->i[pg->imtrx][var] == I_Q2_D_LSA) {
        if (vn_glob[i]->modes > 1) {
          post_flag = 1;
        }
      }
    }
    if (post_flag) {
      TOTAL_STRESS11 = index_post;
      set_nv_tkud(rd, index, var, 0, -2, "ts11", "[1]", "total stress 11", FALSE);
      index++;
      index_post++;
      var = R_STRESS12;
      if (Num_Var_In_Type[pg->imtrx][var]) {
        TOTAL_STRESS12 = index_post;
        set_nv_tkud(rd, index, var, 0, -2, "ts12", "[1]", "total stress 12", FALSE);
        index++;
        index_post++;
      }

      var = R_STRESS22;
      if (Num_Var_In_Type[pg->imtrx][var]) {
        TOTAL_STRESS22 = index_post;
        set_nv_tkud(rd, index, var, 0, -2, "ts22", "[1]", "total stress 22", FALSE);
        index++;
        index_post++;
      }

      var = R_STRESS13;
      if (Num_Var_In_Type[pg->imtrx][var]) {
        TOTAL_STRESS13 = index_post;
        set_nv_tkud(rd, index, var, 0, -2, "ts13", "[1]", "total stress 13", FALSE);
        index++;
        index_post++;
      }

      var = R_STRESS23;
      if (Num_Var_In_Type[pg->imtrx][var]) {
        TOTAL_STRESS23 = index_post;
        set_nv_tkud(rd, index, var, 0, -2, "ts23", "[1]", "total stress 23", FALSE);
        index++;
        index_post++;
      }

      var = R_STRESS33;
      if (Num_Var_In_Type[pg->imtrx][var]) {
        TOTAL_STRESS33 = index_post;
        set_nv_tkud(rd, index, var, 0, -2, "ts33", "[1]", "total stress 33", FALSE);
        index++;
        index_post++;
      }
    }
  }

  if (nn_average > 0) {
    for (i = 0; i < nn_average; i++) {
      if ((pp_average[i]->non_variable_type && pg->imtrx == 0) ||
          (!pp_average[i]->non_variable_type && Num_Var_In_Type[pg->imtrx][pp_average[i]->type])) {
        pp_average[i]->index = index;
        pp_average[i]->index_post = index_post;

        set_nv_tkud(rd, index, 0, 0, -2, pp_average[i]->type_name, "[1]", "average value", FALSE);
        index++;
        index_post++;
      }
    }
  }
  rd->TotalNVPostOutput = index - rd->TotalNVSolnOutput;

  /************************** setup output of time derivatives if requested */

  if (TIME_DERIVATIVES != -1 && (TimeIntegration != STEADY)) {
    if (TIME_DERIVATIVES == 2) {
      GOMA_EH(GOMA_ERROR, "Post-processing time derivatives cannot be exported yet!");
    }
    for (var = V_FIRST; var < V_LAST; var++) {
      if (Num_Var_In_Type[pg->imtrx][var]) {
        if (variable_type_nodalInterp(var)) {
          if (var == MASS_FRACTION) {
            for (mn = -1; mn < upd->Num_Mat; mn++) {
              if (mn == -1) {
                for (i = upd->Num_Mat - 1; i >= 0; i--) {
                  if (mp_glob[i]->Num_Species == upd->Max_Num_Species) {
                    matrl = mp_glob[i];
                  }
                }
              } else {
                matrl = mp_glob[mn];
              }
              for (i = 0; i < matrl->Num_Species; i++) {
                vd = get_vd_ptr(MASS_FRACTION, mn, i);
                if (vd && (vd->MatID == mn)) {
                  /*
                   * Assign the exodus species name
                   */
                  assign_species_name(i, matrl, species_name, ds, mn);
                  strcat(species_name, "DOT");
                  sprintf(species_desc, "Time Derivative of ");
                  strcat(species_desc, ds);
                  /*
                   * Set the values in the Results_Description structure
                   */
                  set_nv_tkud(rd, index, var, i, mn, species_name, "[1]", species_desc, TRUE);
                  index++;
                }
              }
            }
          } else {
            for (mn = -1; mn < upd->Num_Mat; mn++) {
              if (mn == -1) {
                for (i = upd->Num_Mat - 1; i >= 0; i--) {
                  if (pd_glob[i]->i[pg->imtrx][var]) {
                    matrl = mp_glob[i];
                  }
                }
              } else {
                matrl = mp_glob[mn];
              }
              vd = get_vd_ptr(var, mn, i);
              if (vd && (vd->MatID == mn)) {
                assign_var_name(var, 0, matrl, species_name, ds, mn);
                strcat(species_name, "DOT");
                sprintf(species_desc, "Time Derivative of ");
                strcat(species_desc, ds);
                set_nv_tkud(rd, index, var, 0, mn, species_name, "[1]", species_desc, TRUE);
                index++;
              }
            }
          }
        }
      }
    }
  }

  /* end of setup of time derivatives  ****************************/

  /* put index into the results description structure */
  rd->nnv = index;
  *tnv_post = index - *tnv;

  if ((TIME_DERIVATIVES != -1 && pd_glob[0]->TimeIntegration == TRANSIENT) &&
      index_post != (*tnv_post - *tnv))
    GOMA_WH(-1, "Bad nodal post process variable count ");

  if (TIME_DERIVATIVES == -1 && index_post != *tnv_post)
    GOMA_WH(-1, "Bad nodal post process variable count ");

  return (status);
}
/**************************************************************************/
/**************************************************************************/
/**************************************************************************/

int load_elem_tkn(struct Results_Description *rd, const Exo_DB *exo, int tev, int *tev_post) {
  int i, j, index = 0, index_post = 0, ipost, status = 0, check;
  char appended_name[MAX_VAR_NAME_LNGTH];
  int *ev_var_mask;
  static const char yo[] = "load_elem_tkn";
  int ip_total;
  int eb_index, mn;
  ELEM_BLK_STRUCT *eb_ptr;

  /* Here, if a variable and equation is turned on in any one
     material, then we must make provisions in the rd structure
     for all materials.  This is for the sake of the post-processor.
  */

  ev_var_mask = (int *)smalloc((V_LAST - V_FIRST) * sizeof(int));
  for (i = 0; i < V_LAST - V_FIRST; i++) {
    ev_var_mask[i] = 0;
  }
  /* First cycle through all the primary variables that are nodal looking
     for element variable candidates (currently must be interpolated with I_P0 */
  for (i = 0; i < exo->num_elem_blocks; i++) {
    int mat_id = Matilda[i];
    if (exo->eb_num_elems[i] <= 0)
      continue;
    for (j = V_FIRST; j < V_LAST; j++) {
      if (pd_glob[mat_id]->v[pg->imtrx][j] != V_NOTHING) {
        if (pd_glob[mat_id]->i[pg->imtrx][j] == I_P0) {
          if (Num_Var_In_Type[pg->imtrx][j] > 1) {
            fprintf(stderr,
                    "%s: Too many components in variable type (%s - %s) for element variable\n", yo,
                    Exo_Var_Names[j].name2, Exo_Var_Names[j].name1);
            exit(-1);
          }
          if (ev_var_mask[j - V_FIRST] == 0) {
            /* We just found a candidate for an element variable */
            /* Append a suffix onto the var name to differentiate from its
             nodal counterpart */
            sprintf(appended_name, "%s_E", Exo_Var_Names[j].name2);
            set_ev_tkud(rd, index, j, appended_name, Var_Units[j].name2, Exo_Var_Names[j].name1,
                        FALSE);
            index++;
            ev_var_mask[j - V_FIRST] = 1; /* Only count this variable once */
          }
        }
        if (pd_glob[mat_id]->i[pg->imtrx][j] == I_P1) {
          int dof = getdofs(type2shape(exo->eb_elem_itype[i]), I_P1);
          if (ev_var_mask[j - V_FIRST] == 0) {
            /* We just found a candidate for an element variable */
            /* Append a suffix onto the var name to differentiate from its
             nodal counterpart */
            for (i = 1; i <= dof; i++) {
              sprintf(appended_name, "%s_E%d", Exo_Var_Names[j].name2, i);
              set_ev_tkud(rd, index, j, appended_name, Var_Units[j].name2, Exo_Var_Names[j].name1,
                          FALSE);
              index++;
              ev_var_mask[j - V_FIRST] = 1; /* Only count this variable once */
            }
          }
        }
      }
    }
  }

  free(ev_var_mask);

  /* Now pick up all the post processing variables - yes, for now they must
     each be listed separately and painfully */

  /* ZZ error based on the fluid shear stress energy norm */
  if (ERROR_ZZ_VEL != -1 && Num_Var_In_Type[pg->imtrx][R_MOMENTUM1]) {
    set_ev_tkud(rd, index, 0, "ERR_ZZ_V", "[1]", "Zienkiewicz-Zhu error from velocity", FALSE);
    index++;
    ERROR_ZZ_VEL = index_post;
    index_post++;

    if (ERROR_ZZ_VEL_ELSIZE != -1) {
      set_ev_tkud(rd, index, 0, "ERR_ZZ_V_ELSIZE", "[1]",
                  "Target element size from Zienkiewicz-Zhu error from velocity", FALSE);
      index++;
      ERROR_ZZ_VEL_ELSIZE = index_post;
      index_post++;
    }
  } else {
    ERROR_ZZ_VEL = -1;
  }

  /* ZZ error based on the heat flux energy norm */
  if (ERROR_ZZ_Q != -1 && Num_Var_In_Type[pg->imtrx][R_ENERGY]) {
    set_ev_tkud(rd, index, 0, "ERR_ZZ_Q", "[1]", "Zienkiewicz-Zhu error from heat flux", FALSE);
    index++;
    ERROR_ZZ_Q = index_post;
    index_post++;
    if (ERROR_ZZ_Q_ELSIZE != -1) {
      set_ev_tkud(rd, index, 0, "ERR_ZZ_Q_ELSIZE", "[1]",
                  "Target element size from Zienkiewicz-Zhu error from heat flux", FALSE);
      index++;
      ERROR_ZZ_Q_ELSIZE = index_post;
      index_post++;
    }
  } else {
    ERROR_ZZ_Q = -1;
  }

  /* ZZ error based on the pressure energy norm */
  /* RRL I had to set check back to zero so the code wouldn't stop
     in wr_exo.c with a create_truth_table variable count error.  prs */
  check = 0;
  for (i = 0; i < upd->Num_Mat; i++) {
    if (pd_glob[i]->MeshMotion == LAGRANGIAN || pd_glob[i]->MeshMotion == DYNAMIC_LAGRANGIAN)
      check = 0;
  }
  if (ERROR_ZZ_P != -1 && (Num_Var_In_Type[pg->imtrx][R_MOMENTUM1] || check)) {
    set_ev_tkud(rd, index, 0, "ERR_ZZ_P", "[1]", "Zienkiewicz-Zhu error from pressure", FALSE);
    index++;
    ERROR_ZZ_P = index_post;
    index_post++;
    if (ERROR_ZZ_P_ELSIZE != -1) {
      set_ev_tkud(rd, index, 0, "ERR_ZZ_P_ELSIZE", "[1]",
                  "Target element size from Zienkiewicz-Zhu error from pressure", FALSE);
      index++;
      ERROR_ZZ_P_ELSIZE = index_post;
      index_post++;
    }
  } else {
    ERROR_ZZ_P = -1;
  }

  /* Now, if necessary, cycle through and look for Hysteretic saturation models
   * which require more element variables (viz. corresponding to the total number
   * of Gauss points) for restart capability
   */

  /* For now assume that if one block contains the element variables
   * for saturation hysteresis, then they all do.  This is inefficient but
   * to fix you need to make tev_post a element_block dependent array.
   */
  ipost = FALSE;
  for (eb_index = 0; eb_index < exo->num_elem_blocks; eb_index++) {
    if (exo->eb_num_elems[eb_index] > 0) {
      mn = Matilda[eb_index];
      if (mn < 0) {
        continue;
      }
      mp = mp_glob[mn];
      eb_ptr = Element_Blocks + eb_index;
      ip_total = elem_info(NQUAD, eb_ptr->Elem_Type);
      if ((mp->PorousMediaType == POROUS_UNSATURATED ||
           mp->PorousMediaType == POROUS_SHELL_UNSATURATED ||
           mp->PorousMediaType == POROUS_TWO_PHASE) &&
          mp->SaturationModel == TANH_HYST && !ipost) {
        ipost = TRUE;

        SAT_CURVE_TYPE = index_post;
        for (j = 0; j < ip_total; j++) {
          /* We just found more element variables */
          /* Append a index suffix onto the var name to differentiate
             between gauss point values*/
          sprintf(appended_name, "sat_curve_type%d", j);
          set_ev_tkud(rd, index, 0, appended_name, "[1]", "Saturation hysteresis curve type",
                      FALSE);
          index++;
          index_post++;
        }

        SAT_QP_SWITCH = index_post;
        for (j = 0; j < ip_total; j++) {
          sprintf(appended_name, "sat_switch%d", j);
          set_ev_tkud(rd, index, 0, appended_name, "[1]",
                      "Value of saturation at hysteresis switch", FALSE);
          index++;
          index_post++;
        }

        CAP_PRESS_SWITCH = index_post;
        for (j = 0; j < ip_total; j++) {
          sprintf(appended_name, "pc_switch%d", j);
          set_ev_tkud(rd, index, 0, appended_name, "[1]", "Value of cap press at hysteresis switch",
                      FALSE);
          index++;
          index_post++;
        }
      }
    }
  }

  /* Finally put index into the results description structure */
  rd->nev = index;
  *tev_post = index - tev;
  Num_Elem_Post_Proc_Var = index_post;

  if ((TIME_DERIVATIVES != -1 && pd_glob[0]->TimeIntegration == TRANSIENT) &&
      index_post != (*tev_post - tev))
    GOMA_WH(-1, "Bad elem post process variable count ");

  if (TIME_DERIVATIVES == -1 && index_post != *tev_post)
    GOMA_WH(-1, "Bad elem post process variable count ");

  return (status);
}

int find_id_edge(const int ielem, /* 0-based element number */
                 const int num_nodes_on_edge,
                 const int local_edge_node_list[],
                 int id_local_elem_coord[],
                 int *param_dir, /* direction of parametric edge curve */
                 const Exo_DB *exo)

/*
 * this routine finds the edge id for 3D hex elements - it is an extension of the
 * id_side arguments - the ordering of sides and edges for a 3D element are shown
 * below
 *
 *                                    Back Face = 3
 *
 *
 *                                 7------=8=--------6
 *                                /:                /|
 *                Top Face = 6   / :               / |
 *                              5  :              6  |
 *                             /  10             /  12
 *                            /    :            /    |      Right Face = 2
 *                           4-------=7=-------5     |
 *                           |     :           |     |
 *                           |     3......=4=..|.....2
 *                           9    :           11    /
 *      Left Face = 4        |   :             |   /
 *                          ^|  1              |  2
 *                          || :/              | /  Bottom Face = 5
 *                          u|:t               |/
 *                           0-------=3=-------1
 *                            s->
 *
 *                             Front Face = 1
 *
 * The routine also determines the direction of the parameterization of the edge (i.e. edge number
 * 1 is paramterized by t (local coordinate 1)
 *
 * Written by: Richard Cairncross 8 August 1996
 *
 *  TAB certifies that this function conforms to the exo/patran side numbering convention 11/9/98.
 *
 */

{
  int ielem_type;
  int num_local_nodes;
  int ielem_dim;
  int iconnect_ptr;
  int i;
  double sum;
  double u_val;

  /*-------------------------------Start Execution-----------------------------*/

  /* Find out what type of element this is */
  ielem_type = Elem_Type(exo, ielem);

  /* Find out how many local basis functions there are */
  num_local_nodes = elem_info(NNODES, ielem_type);

  /* Find out the physical dimension of the element */
  ielem_dim = elem_info(NDIM, ielem_type);

  /* find the pointer the beginning of this element's connectivity list */
  iconnect_ptr = exo->elem_ptr[ielem];

  /* Find the local element coordinates - note, this algorithm keeps the ordering
                 in local_ss_node_list - might be necessary */

  for (i = 0; i < num_nodes_on_edge; i++) {
    id_local_elem_coord[i] =
        in_list(local_edge_node_list[i], 0, num_local_nodes, &(exo->node_list[iconnect_ptr]));
    if ((id_local_elem_coord[i]) == -1) {
      GOMA_EH(GOMA_ERROR, "find_id_edge ERROR: side set nodal map error\n");
    }
  }

  *param_dir = -1;
  u_val = 0.;
  switch (ielem_dim) {
  case 3:
    u_val = 1.;
    for (i = 0, sum = 0.0; i < num_nodes_on_edge; i++)
      sum += shape(1.0, 1.0, 0.0, ielem_type, PSI, id_local_elem_coord[i]);
    if (sum > 0.999) {
      *param_dir = 2;
      return (12);
    }
    for (i = 0, sum = 0.0; i < num_nodes_on_edge; i++)
      sum += shape(1.0, -1.0, 0.0, ielem_type, PSI, id_local_elem_coord[i]);
    if (sum > 0.999) {
      *param_dir = 2;
      return (11);
    }
    for (i = 0, sum = 0.0; i < num_nodes_on_edge; i++)
      sum += shape(-1.0, 1.0, 0.0, ielem_type, PSI, id_local_elem_coord[i]);
    if (sum > 0.999) {
      *param_dir = 2;
      return (10);
    }
    for (i = 0, sum = 0.0; i < num_nodes_on_edge; i++)
      sum += shape(-1.0, -1.0, 0.0, ielem_type, PSI, id_local_elem_coord[i]);
    if (sum > 0.999) {
      *param_dir = 2;
      return (9);
    }
    for (i = 0, sum = 0.0; i < num_nodes_on_edge; i++)
      sum += shape(0.0, 1.0, 1.0, ielem_type, PSI, id_local_elem_coord[i]);
    if (sum > 0.999) {
      *param_dir = 0;
      return (8);
    }
    for (i = 0, sum = 0.0; i < num_nodes_on_edge; i++)
      sum += shape(0.0, -1.0, 1.0, ielem_type, PSI, id_local_elem_coord[i]);
    if (sum > 0.999) {
      *param_dir = 0;
      return (7);
    }
    for (i = 0, sum = 0.0; i < num_nodes_on_edge; i++)
      sum += shape(1.0, 0.0, 1.0, ielem_type, PSI, id_local_elem_coord[i]);
    if (sum > 0.999) {
      *param_dir = 1;
      return (6);
    }
    for (i = 0, sum = 0.0; i < num_nodes_on_edge; i++)
      sum += shape(-1.0, 0.0, 1.0, ielem_type, PSI, id_local_elem_coord[i]);
    if (sum > 0.999) {
      *param_dir = 1;
      return (5);
    }
    for (i = 0, sum = 0.0; i < num_nodes_on_edge; i++)
      sum += shape(0.0, 1.0, -u_val, ielem_type, PSI, id_local_elem_coord[i]);
    if (sum > 0.999) {
      *param_dir = 0;
      return (4);
    }
    for (i = 0, sum = 0.0; i < num_nodes_on_edge; i++)
      sum += shape(0.0, -1.0, -u_val, ielem_type, PSI, id_local_elem_coord[i]);
    if (sum > 0.999) {
      *param_dir = 0;
      return (3);
    }
    for (i = 0, sum = 0.0; i < num_nodes_on_edge; i++)
      sum += shape(1.0, 0.0, -u_val, ielem_type, PSI, id_local_elem_coord[i]);
    if (sum > 0.999) {
      *param_dir = 1;
      return (2);
    }
    for (i = 0, sum = 0.0; i < num_nodes_on_edge; i++)
      sum += shape(-1.0, 0.0, -u_val, ielem_type, PSI, id_local_elem_coord[i]);
    if (sum > 0.999) {
      *param_dir = 1;
      return (1);
    }
    break;
  case 2:
    for (i = 0, sum = 0.0; i < num_nodes_on_edge; i++)
      sum += shape(1.0, 1.0, 0.0, ielem_type, PSI, id_local_elem_coord[i]);
    if (sum > 0.999) {
      *param_dir = 2;
      return (3);
    }
    for (i = 0, sum = 0.0; i < num_nodes_on_edge; i++)
      sum += shape(1.0, -1.0, 0.0, ielem_type, PSI, id_local_elem_coord[i]);
    if (sum > 0.999) {
      *param_dir = 2;
      return (2);
    }
    for (i = 0, sum = 0.0; i < num_nodes_on_edge; i++)
      sum += shape(-1.0, 1.0, 0.0, ielem_type, PSI, id_local_elem_coord[i]);
    if (sum > 0.999) {
      *param_dir = 2;
      return (4);
    }
    for (i = 0, sum = 0.0; i < num_nodes_on_edge; i++)
      sum += shape(-1.0, -1.0, 0.0, ielem_type, PSI, id_local_elem_coord[i]);
    if (sum > 0.999) {
      *param_dir = 2;
      return (1);
    }
    break;
  case 1:
    /* On a 1D element, edges are points, in ExoII structures, they are
     * certainly nodes.
     * * * */
    sum = shape(1.0, 0.0, 0.0, ielem_type, PSI, id_local_elem_coord[0]);
    if (sum > 0.999) {
      *param_dir = 0;
      return 2; // These edge numbers have to match cases in find_surf_st
    } else {
      *param_dir = 0;
      return 1; // node number
    }
    /*GOMA_EH(GOMA_ERROR,"Can't have edges in 1D");*/
    break;

  } /* END switch ielem_dim */

  /* An error condition has occurred, if here */
  GOMA_EH(GOMA_ERROR, "find_id_edge ERROR: no edge was found !\n");
  return (0);
}
/****************************************************************************/
/****************************************************************************/

int find_id_edge_TET(const int ielem, /* 0-based element number */
                     const int num_nodes_on_edge,
                     const int local_edge_node_list[],
                     int id_local_elem_coord[],
                     int *param_dir, /* direction of parametric edge curve */
                     const Exo_DB *exo)

/*
 * this routine finds the edge id for 3D tet elements - it is an extension of the
 * id_side arguments - the ordering of sides and edges for a 3D element are shown
 * below
 *
 * The routine also determines the direction of the parameterization of the edge (i.e. edge number
 * 1 is paramterized by t (local coordinate 1)
 *
 * Written by: Randy Schunk 12 Feb 2012
 *
 *  TAB certifies that this function conforms to the exo/patran side numbering convention 11/9/98.
 *
 */

{
  int ielem_type;
  int num_local_nodes;
  int ielem_dim;
  int iconnect_ptr;
  int i;
  double sum;

  /*-------------------------------Start Execution-----------------------------*/

  /* Find out what type of element this is */
  ielem_type = Elem_Type(exo, ielem);

  /* Find out how many local basis functions there are */
  num_local_nodes = elem_info(NNODES, ielem_type);

  /* Find out the physical dimension of the element */
  ielem_dim = elem_info(NDIM, ielem_type);

  /* find the pointer the beginning of this element's connectivity list */
  iconnect_ptr = exo->elem_ptr[ielem];

  /* Find the local element coordinates - note, this algorithm keeps the ordering
                 in local_ss_node_list - might be necessary */

  for (i = 0; i < num_nodes_on_edge; i++) {
    if ((id_local_elem_coord[i] = in_list(local_edge_node_list[i], 0, num_local_nodes,
                                          &(exo->node_list[iconnect_ptr]))) == -1) {
      GOMA_EH(GOMA_ERROR, "find_id_edge_TET ERROR: side set nodal map error\n");
    }
  }

  *param_dir = -1;
  switch (ielem_dim) {
  case 3:
    for (i = 0, sum = 0.0; i < num_nodes_on_edge; i++)
      // sum += shape( 0.5, 0.5, 0.0, ielem_type, PSI, id_local_elem_coord[i]);
      sum += shape(0.5, 0, 0.0, ielem_type, PSI, id_local_elem_coord[i]);
    if (sum > 0.999) {
      *param_dir = 0;
      return (1);
    }
    for (i = 0, sum = 0.0; i < num_nodes_on_edge; i++)
      // sum += shape( 0.5,0.0, 0.5, ielem_type, PSI, id_local_elem_coord[i]);
      sum += shape(0.5, 0.5, 0.0, ielem_type, PSI, id_local_elem_coord[i]);
    if (sum > 0.999) {
      *param_dir = 1;
      return (2);
    }
    for (i = 0, sum = 0.0; i < num_nodes_on_edge; i++)
      // sum += shape(0.0, 0.5, 0.5, ielem_type, PSI, id_local_elem_coord[i]);
      sum += shape(0.0, 0.5, 0.0, ielem_type, PSI, id_local_elem_coord[i]);
    if (sum > 0.999) {
      *param_dir = 1;
      return (3);
    }
    for (i = 0, sum = 0.0; i < num_nodes_on_edge; i++)
      sum += shape(0.0, 0.0, 0.5, ielem_type, PSI, id_local_elem_coord[i]);
    if (sum > 0.999) {
      *param_dir = 2;
      return (4);
    }
    for (i = 0, sum = 0.0; i < num_nodes_on_edge; i++)
      // sum += shape( 0.0, 0.5, 0.0, ielem_type, PSI, id_local_elem_coord[i]);
      sum += shape(0.5, 0.0, 0.5, ielem_type, PSI, id_local_elem_coord[i]);
    if (sum > 0.999) {
      *param_dir = 0;
      return (5);
    }
    for (i = 0, sum = 0.0; i < num_nodes_on_edge; i++)
      // sum += shape( 0.0,0.0, 0.5, ielem_type, PSI, id_local_elem_coord[i]);
      sum += shape(0.0, 0.50, 0.5, ielem_type, PSI, id_local_elem_coord[i]);
    if (sum > 0.999) {
      *param_dir = 2;
      return (6);
    }
    break;
  case 2:
    GOMA_EH(GOMA_ERROR, "find_id_edge_TET not working in 2D");
    break;
  case 1:
    GOMA_EH(GOMA_ERROR, "Can't have edges in 1D");
    break;

  } /* END switch ielem_dim */

  /* An error condition has occurred, if here */
  GOMA_EH(GOMA_ERROR, "find_id_edge_TET ERROR: no edge was found !\n");
  return (0);
}
/****************************************************************************/
/****************************************************************************/
/*ARGSUSED*/
int count_nodes_on_SS(const int ss_id,  /* SS id of Primary Side Set */
                      const int ss_id2, /* SS id of 2nd Side Set for edges*/
                      const int ss_id3, /* SS id of 3rd Side Set for vertices*/
                      const int iconnect_ptr,
                      const int ielem,
                      const int num_local_nodes,
                      int local_ss_node_list[MAX_NODES_PER_SIDE],
                      int local_elem_node_id[MAX_NODES_PER_SIDE])
/*
 * This routine counts the local nodes that are located on a SS (for faces),
 * or the number of nodes which are on two SS (edges) or three SS (vertices).
 *
 * Written by Richard Cairncross 13 August 1996
 */
{
  int J, j, num_nodes_on_side;

  num_nodes_on_side = 0;
  for (j = 0; j < num_local_nodes; j++) {
    J = Proc_Elem_Connect[iconnect_ptr + j];
    if (in_list(ss_id, 0, MAX_SS_PER_NODE, SS_list[J]) != -1) {
      if (ss_id2 == -1 || in_list(ss_id2, 0, MAX_SS_PER_NODE, SS_list[J]) != -1) {
        if (ss_id3 == -1 || in_list(ss_id3, 0, MAX_SS_PER_NODE, SS_list[J]) != -1) {
          local_ss_node_list[num_nodes_on_side] = J;
          local_elem_node_id[num_nodes_on_side] = j;
          num_nodes_on_side++;
        }
      }
    }
  }
  return (num_nodes_on_side);
}
/************************************************************************************/
/************************************************************************************/
/************************************************************************************/

int elem_order_for_nodal_connect(int *listel, const Exo_DB *exo)

/********************************************************************************
 *
 * elem_order_for_nodal_connect():
 *
 *      This routine provides an ordering of the elements such that they
 * provide connectivity with the underlying nodes. In other words, the elements
 * are ordered such that the next element always includes nodes found in
 * elements previously included in the ordering (if at all possible).
 *
 * Output:
 *  listel[] = Connectivity Ordering of the elements. (space for this must
 *             have already been allocated).
 *
 * Return:
 *  This routine returns the number of times the above connectivity principle
 *  has been broken in order to include all elements in the listing.
 *
 ********************************************************************************/
{
  int nbreaks = 0, ielem_type, num_local_nodes, elem_found, J, iel, ielem;
  int *used_elem, *used_node;
  int elem_loop_first, iconnect_ptr, ln, iorder;
  used_elem = alloc_int_1(exo->num_elems, 0);
  used_node = alloc_int_1(exo->num_nodes, 0);
  /*
   * Pick the first element to start off the whole process.
   * Mark down the element's nodes as used.
   */
  iel = 0;
  ielem_type = Elem_Type(exo, iel);
  num_local_nodes = elem_info(NNODES, ielem_type);
  iconnect_ptr = exo->elem_ptr[iel];
  for (ln = 0; ln < num_local_nodes; ln++) {
    J = Proc_Elem_Connect[iconnect_ptr + ln];
    used_node[J]++;
  }
  used_elem[iel] = 1;
  elem_loop_first = 1;
  listel[0] = iel;
  /*
   * Now do the search a total of 1 less than the total number of
   * elements defined on the processor.
   */
  for (iorder = 1; iorder < exo->num_elems; iorder++) {
    elem_found = FALSE;
    /*
     *  Loop over all of the unused elements. Hopefully, the
     *  elem_loop_first gimmic will make this search more O(N)
     *  than O(N**2). However, it is not guarranteed.
     *
     *  If we find an element that has a used node, we use that element
     *  and exit the loop.
     */
    for (ielem = elem_loop_first; ielem < exo->num_elems && !elem_found; ielem++) {
      if (!used_elem[ielem]) {
        ielem_type = Elem_Type(exo, ielem);
        num_local_nodes = elem_info(NNODES, ielem_type);
        iconnect_ptr = exo->elem_ptr[ielem];
        for (ln = 0; ln < num_local_nodes; ln++) {
          J = Proc_Elem_Connect[iconnect_ptr + ln];
          if (used_node[J]) {
            elem_found = TRUE;
            iel = ielem;
            break;
          }
        }
      }
    }
    if (!elem_found) {
      for (ielem = elem_loop_first; ielem < exo->num_elems; ielem++) {
        if (!used_elem[ielem]) {
          nbreaks++;
          iel = ielem;
          elem_found = TRUE;
          break;
        }
      }
    }
    /*
     * OK, we have found the element that we are to use.
     * Mark its nodes, add it to the list, and mark
     * the element as used.
     */
    if (elem_found) {
      ielem_type = Elem_Type(exo, iel);
      num_local_nodes = elem_info(NNODES, ielem_type);
      iconnect_ptr = exo->elem_ptr[iel];
      for (ln = 0; ln < num_local_nodes; ln++) {
        J = Proc_Elem_Connect[iconnect_ptr + ln];
        used_node[J]++;
      }
      used_elem[iel] = 1;
      if (iel == elem_loop_first)
        elem_loop_first++;
    }
    listel[iorder] = iel;
  }
  safer_free((void **)&used_elem);
  safer_free((void **)&used_node);

  return nbreaks;
}
/************************************************************************************/
/************************************************************************************/
/************************************************************************************/

int check_elem_order(const int *listel, const Exo_DB *exo)

/********************************************************************************
 *
 * check_elem_order():
 *
 *      This routine checks the ordering of the elements in an element map,
 * such that they  provide connectivity with the underlying nodes.
 * In other words, the elements
 * are ordered such that the next element always includes nodes found in
 * elements previously included in the ordering (if at all possible).
 *
 * Input:
 *  listel[] = Connectivity Ordering of the elements. (space for this must
 *             have already been allocated).
 *
 * Return:
 *  This routine returns the number of times the above connectivity principle
 *  has been broken in order to include all elements in the listing.
 *  It returns -1 if not all of the elements are in the element map.
 *
 ********************************************************************************/
{
  int nbreaks = 0, ielem_type, num_local_nodes, J, iel;
  int iconnect_ptr, ln, iorder, elem_OK;
  int *used_elem, *used_node;
  used_elem = alloc_int_1(exo->num_elems, 0);
  used_node = alloc_int_1(exo->num_nodes, 0);
  /*
   * Pick the first element to start off the whole process.
   * Mark down the element's nodes as used.
   */

  iel = listel[0];
  ielem_type = Elem_Type(exo, iel);
  num_local_nodes = elem_info(NNODES, ielem_type);
  iconnect_ptr = exo->elem_ptr[iel];
  for (ln = 0; ln < num_local_nodes; ln++) {
    J = Proc_Elem_Connect[iconnect_ptr + ln];
    used_node[J]++;
  }
  used_elem[iel] = 1;
  /*
   * Now do the search a total of 1 less than the total number of
   * elements defined on the processor.
   */
  for (iorder = 1; iorder < exo->num_elems; iorder++) {
    elem_OK = FALSE;
    iel = listel[iorder];
    ielem_type = Elem_Type(exo, iel);
    num_local_nodes = elem_info(NNODES, ielem_type);
    iconnect_ptr = exo->elem_ptr[iel];
    for (ln = 0; ln < num_local_nodes; ln++) {
      J = Proc_Elem_Connect[iconnect_ptr + ln];
      if (used_node[J])
        elem_OK = TRUE;
      used_node[J]++;
    }
    if (!elem_OK) {
      nbreaks++;
    }
    used_elem[iel] = 1;
  }
  /*
   * Check to see that all of the elements are in the map
   */
  for (iorder = 0; iorder < exo->num_elems; iorder++) {
    if (!used_elem[iorder]) {
      nbreaks = -1;
    }
  }
  safer_free((void **)&used_elem);
  safer_free((void **)&used_node);
  return nbreaks;
}
/************************************************************************************/
/************************************************************************************/
/************************************************************************************/
