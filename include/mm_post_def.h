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

/* mm_post_def.h -- definitions used in post processing calculations
 *
 * Notes: contents largely transferred into here from the old "mm_post_proc.h"
 *        so that file might be used exclusively for prototype declarations
 *        of functions defined in mm_post_proc.c
 */

#ifndef GOMA_MM_POST_DEF_H
#define GOMA_MM_POST_DEF_H

#include "el_elm.h"
#include "rf_fem_const.h"
#include "rf_io_const.h"
#include "std.h"

/*
 *  These Define parameters help us parse the force and flux calculation
 *  requests.
 */

#define FORCE_NORMAL             0
#define FORCE_TANGENT1           1
#define FORCE_TANGENT2           2
#define FORCE_X                  3
#define FORCE_Y                  4
#define FORCE_Z                  5
#define VOLUME_FLUX              6
#define SPECIES_FLUX             7
#define HEAT_FLUX                8
#define TORQUE                   9
#define AVERAGE_CONC             10
#define SURF_DISSIP              11
#define AREA                     12
#define VOL_REVOLUTION           13
#define PORE_LIQ_FLUX            14
#define CHARGED_SPECIES_FLUX     15
#define CURRENT_FICKIAN          16
#define CURRENT                  17
#define NEG_LS_FLUX              18
#define POS_LS_FLUX              19
#define N_DOT_X                  20
#define ELEC_FORCE_NORMAL        21
#define ELEC_FORCE_TANGENT1      22
#define ELEC_FORCE_TANGENT2      23
#define ELEC_FORCE_X             24
#define ELEC_FORCE_Y             25
#define ELEC_FORCE_Z             26
#define NET_SURF_CHARGE          27
#define DELTA                    28
#define ACOUSTIC_FLUX_NORMAL     29
#define ACOUSTIC_FLUX_TANGENT1   301
#define ACOUSTIC_FLUX_TANGENT2   311
#define ACOUSTIC_FLUX_X          321
#define ACOUSTIC_FLUX_Y          331
#define ACOUSTIC_FLUX_Z          341
#define ACOUSTIC_INTENSITY       351
#define LS_DCA                   361
#define SHELL_VOLUME_FLUX        30
#define FORCE_X_POS              371
#define FORCE_Y_POS              372
#define FORCE_Z_POS              373
#define FORCE_X_NEG              374
#define FORCE_Y_NEG              375
#define FORCE_Z_NEG              376
#define SPECIES_FLUX_REVOLUTION  398
#define REPULSIVE_FORCE          399
#define POYNTING_X               400
#define POYNTING_Y               401
#define POYNTING_Z               402
#define SCATTERING_CROSS_SECTION 403

#define I_VOLUME                  0
#define I_DISSIP                  1
#define I_SPECIES_MASS            2
#define I_HEAT_ENERGY             3
#define I_MOMX                    4
#define I_MOMY                    5
#define I_MOMZ                    6
#define I_TRACE                   7
#define I_POS_FILL                8
#define I_NEG_FILL                9
#define I_NEG_VX                  10
#define I_NEG_VY                  11
#define I_NEG_VZ                  12
#define I_POS_VX                  13
#define I_POS_VY                  14
#define I_POS_VZ                  15
#define I_LS_ARC_LENGTH           16
#define I_RATE_OF_DEF_II          17
#define I_POROUS_LIQUID_INV       18
#define I_SPEED                   19
#define I_MAG_GRAD_FILL_ERROR     20
#define I_NEG_CENTER_X            21
#define I_NEG_CENTER_Y            22
#define I_NEG_CENTER_Z            23
#define I_POS_CENTER_X            24
#define I_POS_CENTER_Y            25
#define I_POS_CENTER_Z            26
#define I_SURF_SPECIES            27
#define I_ELOADX                  28
#define I_ELOADY                  29
#define I_ELOADZ                  30
#define I_SURF_TEMP               31
#define I_JOULE                   32
#define I_LUB_LOAD                33
#define I_VOLUME_PLANE            34
#define I_POS_VOLPLANE            35
#define I_NEG_VOLPLANE            36
#define I_SPECIES_SOURCE          37
#define I_KINETIC_ENERGY          38
#define I_SHELL_VOLUME            39
#define I_TFMP_FORCE              40
#define I_MASS                    41
#define I_MASS_NEGATIVE_FILL      42
#define I_MASS_POSITIVE_FILL      43
#define I_VORTICITY               44
#define I_GIESEKUS                45
#define I_LAMB_MAG                46
#define I_HELICITY                47
#define I_Q_FCN                   48
#define I_POROUS_LIQUID_INV_2     49
#define I_POROUS_LIQUID_INV_3     50
#define I_EM_ABSORB_CROSS_SECTION 51
#ifdef GOMA_MM_POST_PROC_C
struct Post_Processing_Flux_Names {
  char *name; /* flux string */
  int Index;  /* identifier  */
};

typedef struct Post_Processing_Flux_Names FLUX_NAME_STRUCT;

extern FLUX_NAME_STRUCT pp_flux_names[];
extern int Num_Flux_Names;

struct Post_Processing_Flux_Names pp_flux_names[50] = {
    {"FORCE_NORMAL", FORCE_NORMAL},
    {"FORCE_TANGENT1", FORCE_TANGENT1},
    {"FORCE_TANGENT2", FORCE_TANGENT2},
    {"FORCE_X", FORCE_X},
    {"FORCE_Y", FORCE_Y},
    {"FORCE_Z", FORCE_Z},
    {"VOLUME_FLUX", VOLUME_FLUX},
    {"SHELL_VOLUME_FLUX", SHELL_VOLUME_FLUX},
    {"SPECIES_FLUX", SPECIES_FLUX},
    {"HEAT_FLUX", HEAT_FLUX},
    {"TORQUE", TORQUE},
    {"AVERAGE_CONC", AVERAGE_CONC},
    {"SURF_DISSIP", SURF_DISSIP},
    {"AREA", AREA},
    {"VOL_REVOLUTION", VOL_REVOLUTION},
    {"PORE_LIQ_FLUX", PORE_LIQ_FLUX},
    {"CHARGED_SPECIES_FLUX", CHARGED_SPECIES_FLUX},
    {"CURRENT_FICKIAN", CURRENT_FICKIAN},
    {"CURRENT", CURRENT},
    {"NEG_LS_FLUX", NEG_LS_FLUX},
    {"POS_LS_FLUX", POS_LS_FLUX},
    {"N_DOT_X", N_DOT_X},
    {"ELEC_FORCE_NORMAL", ELEC_FORCE_NORMAL},
    {"ELEC_FORCE_TANGENT1", ELEC_FORCE_TANGENT1},
    {"ELEC_FORCE_TANGENT2", ELEC_FORCE_TANGENT2},
    {"ELEC_FORCE_X", ELEC_FORCE_X},
    {"ELEC_FORCE_Y", ELEC_FORCE_Y},
    {"ELEC_FORCE_Z", ELEC_FORCE_Z},
    {"NET_CHARGE", NET_SURF_CHARGE},
    {"DELTA", DELTA},
    {"ACOUSTIC_FLUX_NORMAL", ACOUSTIC_FLUX_NORMAL},
    {"ACOUSTIC_FLUX_TANGENT1", ACOUSTIC_FLUX_TANGENT1},
    {"ACOUSTIC_FLUX_TANGENT2", ACOUSTIC_FLUX_TANGENT2},
    {"ACOUSTIC_FLUX_X", ACOUSTIC_FLUX_X},
    {"ACOUSTIC_FLUX_Y", ACOUSTIC_FLUX_Y},
    {"ACOUSTIC_FLUX_Z", ACOUSTIC_FLUX_Z},
    {"ACOUSTIC_INTENSITY", ACOUSTIC_INTENSITY},
    {"LS_DCA", LS_DCA},
    {"FORCE_X_POS", FORCE_X_POS},
    {"FORCE_Y_POS", FORCE_Y_POS},
    {"FORCE_Z_POS", FORCE_Z_POS},
    {"FORCE_X_NEG", FORCE_X_NEG},
    {"FORCE_Y_NEG", FORCE_Y_NEG},
    {"FORCE_Z_NEG", FORCE_Z_NEG},
    {"SPECIES_FLUX_REVOLUTION", SPECIES_FLUX_REVOLUTION},
    {"REPULSIVE_FORCE", REPULSIVE_FORCE},
    {"POYNTING_X", POYNTING_X},
    {"POYNTING_Y", POYNTING_Y},
    {"POYNTING_Z", POYNTING_Z},
    {"SCATTERING_CROSS_SECTION", SCATTERING_CROSS_SECTION},
};

int Num_Flux_Names = sizeof(pp_flux_names) / sizeof(struct Post_Processing_Flux_Names);

struct Post_Processing_Volume_Names {
  char *name;
  int Index;
};

typedef struct Post_Processing_Volume_Names VOL_NAME_STRUCT;

extern VOL_NAME_STRUCT pp_vol_names[];
extern int Num_Vol_Names;

VOL_NAME_STRUCT pp_vol_names[] = {{"VOLUME", I_VOLUME},
                                  {"LUB_LOAD", I_LUB_LOAD},
                                  {"DISSIPATION", I_DISSIP},
                                  {"JOULE", I_JOULE},
                                  {"SPECIES_MASS", I_SPECIES_MASS},
                                  {"HEAT_ENERGY", I_HEAT_ENERGY},
                                  {"MOMENTUMX", I_MOMX},
                                  {"MOMENTUMY", I_MOMY},
                                  {"MOMENTUMZ", I_MOMZ},
                                  {"STRESS_TRACE", I_TRACE},
                                  {"POSITIVE_FILL", I_POS_FILL},
                                  {"NEGATIVE_FILL", I_NEG_FILL},
                                  {"NEGATIVE_VX", I_NEG_VX},
                                  {"NEGATIVE_VY", I_NEG_VY},
                                  {"NEGATIVE_VZ", I_NEG_VZ},
                                  {"POSITIVE_VX", I_POS_VX},
                                  {"POSITIVE_VY", I_POS_VY},
                                  {"POSITIVE_VZ", I_POS_VZ},
                                  {"NEGATIVE_CENTER_X", I_NEG_CENTER_X},
                                  {"NEGATIVE_CENTER_Y", I_NEG_CENTER_Y},
                                  {"NEGATIVE_CENTER_Z", I_NEG_CENTER_Z},
                                  {"POSITIVE_CENTER_X", I_POS_CENTER_X},
                                  {"POSITIVE_CENTER_Y", I_POS_CENTER_Y},
                                  {"POSITIVE_CENTER_Z", I_POS_CENTER_Z},
                                  {"LS_ARC", I_LS_ARC_LENGTH},
                                  {"II_GAMMA_DOT", I_RATE_OF_DEF_II},
                                  {"POROUS_LIQ_INVENTORY", I_POROUS_LIQUID_INV},
                                  {"POROUS_LIQ_INVENTORY_2", I_POROUS_LIQUID_INV_2},
                                  {"POROUS_LIQ_INVENTORY_3", I_POROUS_LIQUID_INV_3},
                                  {"SPEED_SQUARED", I_SPEED},
                                  {"SURFACE_SPECIES", I_SURF_SPECIES},
                                  {"ELECTRIC_LOAD_X", I_ELOADX},
                                  {"ELECTRIC_LOAD_Y", I_ELOADY},
                                  {"ELECTRIC_LOAD_Z", I_ELOADZ},
                                  {"SURFACE_TEMPERATURE", I_SURF_TEMP},
                                  {"VOL_PLANE", I_VOLUME_PLANE},
                                  {"POS_PLANE_FILL", I_POS_VOLPLANE},
                                  {"NEG_PLANE_FILL", I_NEG_VOLPLANE},
                                  {"SPECIES_SOURCE", I_SPECIES_SOURCE},
                                  {"KINETIC_ENERGY", I_KINETIC_ENERGY},
                                  {"SHELL_VOLUME", I_SHELL_VOLUME},
                                  {"TFMP_FORCE", I_TFMP_FORCE},
                                  {"MASS", I_MASS},
                                  {"MASS_NEGATIVE_FILL", I_MASS_NEGATIVE_FILL},
                                  {"MASS_POSITIVE_FILL", I_MASS_POSITIVE_FILL},
                                  {"VORTICITY", I_VORTICITY},
                                  {"GIESEKUS", I_GIESEKUS},
                                  {"LAMB_MAG", I_LAMB_MAG},
                                  {"HELICITY", I_HELICITY},
                                  {"Q_FCN", I_Q_FCN},
                                  {"ABSORPTION_CROSS_SECTION", I_EM_ABSORB_CROSS_SECTION}};

int Num_Vol_Names = sizeof(pp_vol_names) / sizeof(VOL_NAME_STRUCT);

#define PP_GLOBAL_START                        0
#define PP_GLOBAL_LS_INTERFACE_PRINT           0
#define PP_GLOBAL_LS_INTERFACE_PRINT_ALL_TIMES 1
#define PP_GLOBAL_COUNT                        2 /* Increase when another is added */

struct Post_Processing_Global_Names {
  char *name; /* flux string */
  int Index;  /* identifier  */
};

typedef struct Post_Processing_Global_Names GLOBAL_NAME_STRUCT;

extern GLOBAL_NAME_STRUCT pp_global_names[];
extern int Num_Global_Names;

struct Post_Processing_Global_Names pp_global_names[PP_GLOBAL_COUNT] = {
    {"LEVEL_SET_INTERFACE_PRINT", PP_GLOBAL_LS_INTERFACE_PRINT},
    {"LEVEL_SET_INTERFACE_PRINT_ALL_TIMES", PP_GLOBAL_LS_INTERFACE_PRINT_ALL_TIMES}};

#endif

/*
 * Define structure to hold details of these requests for surface integrations
 * The idea here is to declare a pointer to this structure and dimension to
 * the number of requested fluxes, or the number of flux cards
 */

struct Post_Processing_Fluxes {
  int flux_type;                    /*flux type reqested */
  char flux_type_name[MAX_DOFNAME]; /*name of flux type requested */
  int species_number;               /*species number if relevant */
  int ss_id;                        /*corresponding ss id request */
  int blk_id;                       /*material from which flux should be calculated*/
  char flux_filenm[MAX_FNL];        /*Neutral file name for output of fluxes */
  int profile_flag;                 /*  flag to control flux profile output */
};

typedef struct Post_Processing_Fluxes pp_Fluxes;

struct Post_Processing_Fluxes_Sens {
  int flux_type;                    /*flux type reqested */
  char flux_type_name[MAX_DOFNAME]; /*name of flux type requested */
  int species_number;               /*species number if relevant */
  int ss_id;                        /*corresponding ss id request */
  int blk_id;                       /*material from which flux should be calculated*/
  int sens_type;                    /* sensitivity variable type  */
  int sens_id;                      /* Boundary condition id or Material id */
  int sens_flt;                     /* data float id or matl property id */
  int sens_flt2;                    /* user matl float id */
  int vector_id;                    /* number of sensitivity vector */
  char flux_filenm[MAX_FNL];        /*Neutral file name for output of fluxes */
  int profile_flag;                 /*  flag to control flux sens. output */
};

typedef struct Post_Processing_Fluxes_Sens pp_Fluxes_Sens;

struct Post_Processing_Data {
  int data_type;                    /*data type reqested */
  char data_type_name[MAX_DOFNAME]; /*name of data type requested */
  int species_number;               /*species number if relevant */
  int ns_id;                        /*corresponding ns id request */
  int mat_num;               /* material num from which data should be calculated (starts with 0) */
  int elem_blk_id;           /* element block id from which data should be calculated*/
  char data_filenm[MAX_FNL]; /*Neutral file name for output of data */
  char format_flag[8];       /* Select output style format. */
  int first_time;            /* boolean to output header once only */
};

typedef struct Post_Processing_Data pp_Data;

struct Post_Processing_Data_Sens {
  int data_type;                    /*data type reqested */
  char data_type_name[MAX_DOFNAME]; /*name of data type requested */
  int species_number;               /*species number if relevant */
  int ns_id;                        /*corresponding ns id request */
  int mat_id;                       /* material id from which data should be calculated*/
  int sens_type;                    /* sensitivity variable type  */
  int sens_id;                      /* Boundary condition id or Material id */
  int sens_flt;                     /* data float id or matl property id */
  int sens_flt2;                    /* user matl float id */
  int vector_id;                    /* number of sensitivity vector */
  char data_filenm[MAX_FNL];        /*Neutral file name for output of data */
};

typedef struct Post_Processing_Data_Sens pp_Data_Sens;

struct Post_Processing_Volumetric {
  int volume_type;               /* volume integral type requested */
  char volume_name[MAX_DOFNAME]; /* string name of volume integral type */
  int species_no;                /* species number if relevant */
  int blk_id;                    /* Element block id on which integral should be evaluated */
  char volume_fname[MAX_FNL];    /* output filename */
  double *params;                /* some additional double constants */
  int num_params;                /* number of additional constants	*/
};

typedef struct Post_Processing_Volumetric pp_Volume;

/*
 * Define structure to hold details of particle tracking specifications
 */

struct Post_Processing_Particles {
  double coord[DIM];      /*Initial Coordinates for each particle*/
  double xi_coord[DIM];   /*Initial Coordinates for each particle*/
  double p_velo[DIM];     /*Particle velocity components*/
  int Current_element_id; /*Initial element id for each particle */
  double Start_Time;      /*Start time */
  double End_Time;        /*End time */
  double Delta_s;         /*Point Spacing */
  char filenm[MAX_FNL];   /*  file name for data output */
  double mass;            /* mass of particle  */
  double mobility;        /* mobility of particle  */
  double force[DIM];      /* external force on particle */
};

typedef struct Post_Processing_Particles pp_Particles;

struct Post_Processing_Error {
  int error_type;      /* error type reqested */
  dbl error_params[6]; /* Five values are to be supplied:
                             [0] = elem size reduction rate
                             [1] = elem size expansion rate
                             [2] = minimum elem size
                             [3] = maximum elem size
                             [4] = target error value
                             [5] = Volume over error target tolerance [%]    */
};

typedef struct Post_Processing_Error pp_Error;

struct Post_Processing_Global {
  int type;
  char type_name[MAX_DOFNAME]; /*name of flux type requested */
  char filenm[MAX_FNL];
};

typedef struct Post_Processing_Global pp_Global;

typedef struct Post_Processing_Averages {
  int type;
  char type_name[MAX_VAR_NAME_LNGTH];
  int species_index;
  int index_post;
  int index;
  int non_variable_type;
} pp_Average;

enum AverageExtraTypes {
  AVG_DENSITY,
  AVG_HEAVISIDE,
  AVG_VISCOSITY,
  AVG_SHEAR,
  AVG_EMR_X,
  AVG_EMR_Y,
  AVG_EMR_Z,
  AVG_EMI_X,
  AVG_EMI_Y,
  AVG_EMI_Z,
  AVG_EMSCATR_X,
  AVG_EMSCATR_Y,
  AVG_EMSCATR_Z,
  AVG_EMSCATI_X,
  AVG_EMSCATI_Y,
  AVG_EMSCATI_Z,
  AVG_EMINCR_X,
  AVG_EMINCR_Y,
  AVG_EMINCR_Z,
  AVG_EMINCI_X,
  AVG_EMINCI_Y,
  AVG_EMINCI_Z,
  AVG_EM_MAG,
  AVG_EM_INC_MAG,
  AVG_EM_SCAT_MAG
};

/*
 * All of these variables are actually defined in mm_post_proc.c
 *
 * The nn_* integers contain information about the actual sizes of these
 * arrays of structures...
 */

extern pp_Fluxes **pp_fluxes;
extern pp_Fluxes_Sens **pp_fluxes_sens;
extern pp_Data **pp_data;
extern pp_Data_Sens **pp_data_sens;
extern pp_Error *pp_error_data;
extern pp_Particles **pp_particles;
extern pp_Volume **pp_volume;
extern pp_Global **pp_global;
extern pp_Average **pp_average;

extern int nn_post_fluxes;
extern int nn_post_fluxes_sens;
extern int nn_post_data;
extern int nn_post_data_sens;
extern int nn_error_metrics;
extern int nn_particles;
extern int nn_volume;
extern int ppvi_type;
extern int nn_global;
extern int nn_average;

extern int Num_Nodal_Post_Proc_Var;
extern int Num_Elem_Post_Proc_Var;

#if 0 /* these are def'd in rf_bc_const.h */
extern struct Equation_Names Exo_Var_Names[];
extern int Num_Exo_Var_Names;  
extern struct Equation_Names Var_Units[];
extern int Num_Var_Units;
#endif

/* the following variables are flags for input options - i.e. 0 or 1,
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
 * -1 -- "No, do not calculate this quantity during post processing."
 */

extern int CAPILLARY_PRESSURE; /* capillary pressure in a porous media */
extern int CONC_CONT;          /* concentration at vertex & midside nodes*/
extern int CONDUCTION_VECTORS; /* conduction flux vectors*/

extern int CURL_V;             /* Steve Kempka's favorite quantity */
extern int DARCY_VELOCITY_GAS; /* Darcy velocity vectors for gas phase
                                * flow in a partially saturated porous
                                * media */
extern int DARCY_VELOCITY_LIQ; /* Darcy velocity vectors for flow in a
                                * saturated or unsaturated medium */
extern int DENSITY;            /* density function at vertex and midside
                                * nodes, e.g. for particle settling etc. */
extern int POLYMER_VISCOSITY;
extern int POLYMER_TIME_CONST;
extern int MOBILITY_PARAMETER;
extern int PTT_XI;
extern int PTT_EPSILON;
extern int DIELECTROPHORETIC_FIELD;
/* Dielectrophoretic force vectors. */
extern int DIELECTROPHORETIC_FIELD_NORM;
/* Norm of DIELECTROPHORETIC_VECTORS. */
extern int ENORMSQ_FIELD;      /* grad(|E|^2) */
extern int ENORMSQ_FIELD_NORM; /* |grad(|E|^2)| */
extern int DIFFUSION_VECTORS;  /* diffusion flux vectors*/
extern int DIFFUSION_VECTORS_POR_LIQ_GPHASE;
/* Diffusion of the solvent liquid in the
 * gas phase for porous flow problems */
extern int DIFFUSION_VECTORS_POR_AIR_GPHASE;
/* Diffusion of the air in the
 * gas phase for porous flow problems */
extern int DIV_PVELOCITY;       /* check the divergence of the particle phase
                                 * velocities.  */
extern int DIV_TOTAL;           /* Divergence of the sum of the fluid and
                                 * particle phases.  This should be zero. */
extern int DIV_VELOCITY;        /* incompressibility constraint at vertex and
                                 * midside nodes, e.g. del*v = 0 */
extern int ELECTRIC_FIELD;      /* Electric field vectors: E = -grad(VOLTAGE) */
extern int ELECTRIC_FIELD_MAG;  /* Electric field magnitude: sqrt(E.E) */
extern int ENERGY_FLUXLINES;    /* energy flux function, analogous to
                                 * stream function ... */
extern int ERROR_ZZ_P;          /* Zienkiewicz-Zhu error indicator (element
                                 * quantity) based solely on pressure
                                 * contributions */
extern int ERROR_ZZ_P_ELSIZE;   /* Recommended new element size from ZZ
                                 * pressure measure                          */
extern int ERROR_ZZ_Q;          /* Zienkiewicz-Zhu error indicator (element
                                 * quantity) based solely on heat flux
                                 * contributions                             */
extern int ERROR_ZZ_Q_ELSIZE;   /* Recommended new element size from ZZ heat
                                 * flux measure */
extern int ERROR_ZZ_VEL;        /* Zienkiewicz-Zhu error indicator (element
                                 * quantity) based solely on velocity (shear
                                 * stress) contributions                     */
extern int ERROR_ZZ_VEL_ELSIZE; /* Recommended new element size from ZZ
                                 * velocity measure */
extern int EVP_DEF_GRAD_TENSOR;
extern int EXTERNAL_POST; /* external field variables read from other
                           * files */
extern int FILL_CONT;     /* fill at vertex & midside nodes*/
extern int FIRST_INVAR_STRAIN;
extern int FLUXLINES;           /* mass flux function. This is analogous to
                                 * stream function but represents mass flux */
extern int LAGRANGE_CONVECTION; /* Lagrangian convection velocity */
extern int MEAN_SHEAR;
extern int MM_RESIDUALS; /* stress equation residuals at vertex
                          * and midside nodes*/
extern int NS_RESIDUALS; /* Navier-Stokes residuals at vertex
                          * and midside nodes */
extern int POROUS_RHO_GAS_SOLVENTS;
/* gas phase concentration of each solvent
 * species in a porous media */
extern int POROUS_RHO_LPHASE; /* liquid phase density per unit volume of
                               * material in a porous medium */
extern int POROUS_RHO_TOTAL_SOLVENTS;
/* Total density of each solvent species in a
 * porous media. Total, here means the
 * density summed up over all phases */
extern int POROUS_SATURATION;   /* saturation in a partially-saturated
                                 * porous media */
extern int POROUS_GRIDPECLET;   /* Grid Peclet number for porous media */
extern int POROUS_SUPGVELOCITY; /* Effective velocities to use in SUPG
                                 * formulations in porous media */
extern int POROUS_LIQUID_ACCUM_RATE;
/* The rate at which liquid in a partially
 * saturated porous medium is accumulating
 * at a point */
extern int REL_LIQ_PERM;    /* Relative liquid permeability in porous media */
extern int PRESSURE_CONT;   /* pressure at vertex & midside nodes*/
extern int SH_DIV_S_V_CONT; /* SH_DIV_S_V at midside nodes */
extern int SH_CURV_CONT;    /* SH_SURF_CURV at midside nodes */
extern int REAL_STRESS_TENSOR;
extern int SEC_INVAR_STRAIN;     /* 2nd strain invariant vertex,midside nodes*/
extern int STRAIN_TENSOR;        /* strain tensor for mesh deformation  */
extern int STREAM;               /* stream function*/
extern int STREAM_NORMAL_STRESS; /* streamwise normal stress function*/
extern int STREAM_SHEAR_STRESS;  /* streamwise shear stress function*/
extern int STREAM_TENSION;       /* streamwise Stress Difference*/
extern int STRESS_CONT;          /* stress at vertex & midside nodes*/
extern int STRESS_TENSOR;        /* stress tensor for mesh deformation
                                  * (Lagrangian pressure) */
extern int SURFACE_VECTORS;      /* vector field of normals and tangents on
                                  * surfaces, curves and vertices */
extern int SHELL_NORMALS;        /* vector field of smoothed normals computed from fv->sh_ang */
extern int THIRD_INVAR_STRAIN;
extern int TIME_DERIVATIVES;           /* time derivatives */
extern int TOTAL_STRESS11;             /* sum over all modes for multi-mode models */
extern int TOTAL_STRESS12;             /* sum over all modes for multi-mode models */
extern int TOTAL_STRESS13;             /* sum over all modes for multi-mode models */
extern int TOTAL_STRESS22;             /* sum over all modes for multi-mode models */
extern int TOTAL_STRESS23;             /* sum over all modes for multi-mode models */
extern int TOTAL_STRESS33;             /* sum over all modes for multi-mode models */
extern int USER_POST;                  /* a user defined function */
extern int PP_Viscosity;               /* Value of the fluid viscosity */
extern int PP_FlowingLiquid_Viscosity; /* Value of the fluid flowing liquid viscosity (Porous
                                          Brinkman term) */
extern int PP_VolumeFractionGas;       /* Value of the volume fraction of the gas component
                                          in a foam or other two phase material */

extern int len_u_post_proc;         /* size of dynamically allocated u_post_proc
                                     * actually is */
extern double *u_post_proc;         /* user-provided values used in calculating
                                     * user defined post processing variable */
extern int SAT_CURVE_TYPE;          /*Saturation hysteresis curve type */
extern int CAP_PRESS_SWITCH;        /*Capillary pressure at hysteresis switch */
extern int SAT_QP_SWITCH;           /*Saturation at hysteresis switch*/
extern int NUM_CURVE_SWITCH;        /* Number of hysteretic curve switch */
extern int SAT_HYST_MIN;            /* Minimum saturation value for scanning imbibition curve */
extern int SAT_HYST_MAX;            /* Maximum saturation value for scanning draining curve */
extern int ACOUSTIC_PRESSURE;       /* Acoustic Pressure Magnitude	*/
extern int ACOUSTIC_PHASE_ANGLE;    /* Acoustic Pressure Phase Angle	*/
extern int ACOUSTIC_ENERGY_DENSITY; /* Acoustic Energy Density	*/
extern int LIGHT_INTENSITY;         /* Light Intensity	*/
extern int PRINCIPAL_STRESS;        /* Principal Stresses*/
extern int PRINCIPAL_REAL_STRESS;   /* Principal Real Stresses*/
extern int LUB_HEIGHT;              /* Lubrication gap*/
extern int LUB_HEIGHT_2;            /* Lubrication gap, second layer*/
extern int LUB_VELO_UPPER;          /* Lubrication upper surface velocity*/
extern int LUB_VELO_LOWER;          /* Lubrication lower surface velocity*/
extern int LUB_VELO_FIELD;          /* Velocity field calculated from lubrication */
extern int LUB_VELO_FIELD_2;        /* Velocity field calculated from lubrication, second layer */
extern int DISJ_PRESS;              /* Disjoining pressure */
extern int SH_SAT_OPEN;             /* Saturation for open porous shells */
extern int SH_SAT_OPEN_2;           /* Saturation for open porous shells 2 */
extern int SH_CAP_PRES;             /* Capillary pressure for porous shell */
extern int SH_PORE_FLUX;            /* Flux between porous shell layers */
extern int SH_STRESS_TENSOR;        /* stress tensor for structural shell */
extern int SH_TANG;                 /* Tangents vectors for structural shell */
extern int PP_LAME_MU;              /* Lame MU coefficient for solid/mesh */
extern int PP_LAME_LAMBDA;          /* Lame LAMBDA coefficient for solid/mesh */
extern int VON_MISES_STRAIN;
extern int VON_MISES_STRESS;
extern int CONF_MAP; /* Map log-conformation tensor to stress */
extern int J_FLUX;   /* Particle stress flux                  */
extern int EIG;      /* Eigenvalues of rate-of-strain tensor  */
extern int EIG1;     /* Eigenvector of rate-of-strain tensor  */
extern int EIG2;     /* Eigenvector of rate-of-strain tensor  */
extern int EIG3;     /* Eigenvector of rate-of-strain tensor  */
extern int HEAVISIDE;
extern int RHO_DOT;
extern int MOMENT_SOURCES;
extern int YZBETA;
extern int GRAD_Y;         /* Concentration gradient                  */
extern int GRAD_SH;        /* Shear gradient                */
extern int UNTRACKED_SPEC; /*Untracked Species Concentration */

extern int TFMP_GAS_VELO;
extern int TFMP_LIQ_VELO;
extern int TFMP_INV_PECLET;
extern int TFMP_KRG;
extern int VELO_SPEED;       /* i.e., velocity magnitude */
extern int GIES_CRIT;        /* Giesekus Flow Character */
extern int HELICITY;         /* v dot vorticity  */
extern int LAMB_VECTOR;      /* Lamb Vector = vorticity x v  */
extern int Q_FCN;            /* 2nd invariant of grad_v  */
extern int POYNTING_VECTORS; /* EM Poynting Vectors*/
extern int PSPG_PP;          /* 2nd invariant of grad_v  */
extern int SARAMITO_YIELD;
extern int STRESS_NORM;
extern int SPECIES_SOURCES; /* Species sources */
extern int VISCOUS_STRESS;  /* Viscous stress */
extern int VISCOUS_STRESS_NORM;
extern int VISCOUS_VON_MISES_STRESS;
extern int EM_CONTOURS;
extern int TOTAL_EM_CONTOURS;
extern int SCATTERED_EM_CONTOURS;
extern int ORIENTATION_VECTORS;
extern int FIRST_STRAINRATE_INVAR;
extern int SEC_STRAINRATE_INVAR;
extern int THIRD_STRAINRATE_INVAR;
/*
 *  Post-processing Step 1: add a new variable flag to end of mm_post_proc.h
 *
 *       e.g.  int STREAM;
 *
 *       Note that this flag is now -1 (false) or the id number of the
 *       post-processing variable in load_nodal_tkn
 */

#endif /* GOMA_MM_POST_DEF_H */
