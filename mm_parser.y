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

%{



#include <stdlib.h>
#include <stdio.h>
#include <strings.h>
#include <string.h>
#include <ctype.h>	
#include <time.h>	
#include <unistd.h>

#include "std.h"
#include "rf_io_const.h"
#include "rf_vars_const.h"
#include "rf_fem_const.h"
#include "rf_fem.h"
#include "rf_solver.h"
#include "rf_mp.h"
#include "rf_io_structs.h"
/*#include "rf_io.h"*/
#include "mm_as_structs.h"
#include "rf_bc_const.h"
#include "rf_allo.h"
#include "rf_bc.h"
#include "mm_augc_util.h"

#include "mm_as_const.h"

#include "mm_as.h"
#include "mm_as_alloc.h"

#include "mm_mp_structs.h"
#include "mm_mp.h"
#include "mm_mp_const.h"

#include "mm_eh.h"

#include "mm_post_def.h"

#include "sl_util_structs.h"
#include "sl_util.h"
#include "mm_input.h"

/*#include "goma.h"*/
/*#include "rf_io_defn.h"*/
#include "mm_species.h"
#include "rd_mesh.h"
#include "ac_stability.h"
#include "sl_eggroll_def.h"
#include "math.h"



#include "mm_parser.h"
#include "mm_cards.h"

#undef yywrap
extern char yytext[];

/***********************/
/* Parser Declarations */
/***********************/

FILE	*in_file, *parser_log, *error_log;		/* string is read from in_file, parser_log is the Goma parser log file */

struct	tm *loctime;			/* needed to print out time stamp to parser log file */
time_t	curtime;			/* needed to print out time stamp to parser log file */
char	debug_clue[100];
char	msg[100];
char 	*mystring;
char    test[100];	
float	NA = -9999999;
int	read_bc_mp = -1;
int	number_of_species_source_cards_needed[MAX_NUMBER_MATLS];
int 	Num_Var_Init_Mat[MAX_NUMBER_MATLS];
float	floating_point_constant_list_array[MAX_NUMBER_PARAMS];
float	table_array[100][5];
int     first_default = TRUE;
int	accept_post_processing_particles = TRUE;
int	accept_data_sens_cards = TRUE;
int	accept_post_processing_data_until_end_of_ppdata_card = TRUE;
int	accept_ppfsens_until_end_of_ppfsens_card = TRUE;
int	accept_hunting_conditions = FALSE;
int	accept_acs_until_end_of_ac_card;
int	accept_rot_until_end_of_rot_card;
int	accept_flux_until_end_of_flux_card;
int	accept_flux_sens_until_end_of_flux_sens_card;
int	accept_volume_int_until_end_of_volume_int_card;
int	accept_data_until_end_of_data_card = FALSE;
int	accept_ppfs_until_end_of_ppf_card = FALSE;
int	accept_hcs_until_end_of_hc_card;
int	accept_bcs_until_end_of_bc_card;
int	accept_ccs_until_end_of_cc_card;
int	accept_tps_until_end_of_tp_card;
int	accept_materials_until_end_of_material_card = FALSE;
int	accept_table_data_cards = 0;
int	accept_eqs_till_end_of_eq_card;
int	accept_material_related_cards = FALSE;
int	number_of_flory_huggins_parameter_cards_to_accept = 0;
int     n_species;
int	accept_level_set_initialization_surface_cards = FALSE;
int	floating_point_constant_list_index = -1;
int	int_array[50];
int	int_array_index = 0;	/* 12/12/00 not being used */
int     model_read;
int	j, k;
int	iBC = -1;
int	iAC = -1;
int	iCC = -1;
int	iTP = -1;
int	iHC = -1;
int	number_of_error_zz_elem_size_cards_found = 0;
int	mat_equation_count = 0;
int	number_of_acs_expected;
int	number_of_ppfs_found = 0;
int	number_of_acs_found = 0;
int	number_of_ccs_expected;
int	number_of_ccs_found = 0;
int	number_of_tps_expected;
int	number_of_tps_found = 0;
int	number_of_data_expected;
int	number_of_data_found = 0;
int	number_of_hcs_found = 0;
int	number_of_hcs_expected;
int	number_of_flux_found = 0;
int	number_of_flux_sens_found = 0;
int	number_of_rots_found = 0;
int	number_of_volume_int_found = 0;
int	number_of_table_cards_read = 0;
int	number_of_bcs_expected;
int	number_of_bcs_found = 0;
int	number_of_equations_expected[9] = {0,0,0,0,0,0,0,0,0};
int	number_of_equations_found[9] = {0,0,0,0,0,0,0,0,0};
int	number_of_materials_expected = 0;
int	number_of_materials_found = 0;
int	number_of_level_set_initialization_surface_cards_expected = 0;
int	number_of_level_set_initialization_surface_cards_found = 0;
int	warning_number = 1;
int 	allocate_GD_bc_description();
int 	allocate_fix_bc_description();
int 	allocate_table_bcs_description();
/* int 	allocate_material_property();*/
int 	card_index;
int 	columns_expected;
int 	common_BC_processing(char *, char *, int);
int 	declare_error(char *msg);
int 	declare_warning ( char * );
int 	declare_defualt_warning ( char * );
int 	error_found_in_last_line = 0;	
int 	error_number = 1;		/* used in yyerror(const char *msg) to number the parse errors that bison finds */
int 	mn = -1;				/* number of the current material file being processed */
int	file_index = 0;
int 	process_equation_and_variable_names(char * equation_name, char * variable_name);
int 	read_array_into_table( enum table_types, struct Data_Table * table, int Num_Pnts );
int 	set_equation ( char *, int, int, int, int, int, float, float, float, float, float, float );
int 	setup_GD_table();
int 	setup_bc_table(int, int, int);
int 	valid_fix_bc_variable_name ( char * variable_name);
int 	yyerror(const char *msg);
int  	create_new_material(char *, int, int, int, int);
int  	set_variable(int , int );
int     eq_index;
int     line_number = 1;
int     var_index;
int	nHC = 0;
int	rdhun = FALSE;
int	AC_rd_file = FALSE;
int	pp_index;
int	i,j, cnt;
int     number_of_constants;
int     sens_vec_ct = 0;
int	mp_table_species_number = -9999999;

struct Data_Table *   	setup_mp_table (struct Data_Table *,int, char *, char *,char *, char *,int, int);

enum	table_types table_type;

struct LS_Init_Object **Surf;
/*struct Post_Processing_Fluxes pp_fluxes[];*/
void    card_read(char *, int);
int times_card_read(char *, int);
void	summary_and_closes();
int common_CC_processing( int ,int ,int ,int , float, float);
int common_TP_processing( int ,int ,int ,int , float, float);
int common_HC_processing( int ,int ,int ,int , float, float, float, float, float);
int common_error_element_size_processing( float, float, float, float, float, float);
int common_pp_flux_processing( char *, int, int, int, int, char *, int);
void postprocess_EQ_cards ();

int nn_post_fluxes = 0;
int nn_post_fluxes_sens = 0;
int nn_post_data = 0;
int nn_post_data_sens = 0;
int nn_particles = 0;
int nn_volume= 0;

MATRL_PROP_STRUCT *mat_ptr;
extern void fallback_chemkin_generic_prop(int*, int, int*, int, MATRL_PROP_STRUCT *);

extern int parse_table_file (char * tablefile, char * tabletag, int    columns);
/*extern int parse_table_file (char * tablefile, char * tabletag, int    columns);*/
extern int yyrestart();
extern int yylex();

#define YYERROR_VERBOSE 1 	

%}

%union	{
	char 	string[100];
	int 	integer;
	char	character;
	int	integer_array[50];
	float	floating;
	}
	
/* Parser File Mode Tokens: */

%token <string> SUSPEND_
%token <string> BOUSSINESQ_
%token <string> BOUSS_JXB_
%token <string> INTEGER_
%token <string> AC_
%token <string> ADAPTIVE_
%token <string> ADVECTIVE_
%token <string> AEXP_
%token <string> ALGORITHM_ 
%token <string> AND_
%token <string> ANNEAL_MESH_ON_OUTPUT
%token <string> ARBITRARY_
%token <string> AUXILIARY_ 
%token <string> A_
%token <string> BASIS_ONCE_ 
%token <string> BASIS_RESEED_ 
%token <string> BC_ 
%token <string> BINGHAM_
%token <string> BODY_
%token <string> BOUNDARY_
%token <string> BOUNDARY_CONDITION_DATA_FLOAT_TAG_
%token <string> BOUNDARY_CONDITION_ID_
%token <string> BOUSS_
%token <string> BULK_
%token <string> BUTLER_VOLMER_IJ_
%token <string> BUTLER_VOLMER_J_
%token <string> BUTLER_VOLMER_
%token <string> B_
%token <string> CAPACITY_
%token <string> CAPILLARY_
%token <string> CAP_REPULSE_
%token <string> CARREAU_
%token <string> CARTESIAN_
%token <string> CA_
%token <string> CA_EDGE_
%token <string> CA_EDGE_CURVE_
%token <string> CA_EDGE_CURVE_INT_
%token <string> CA_EDGE_INT_
%token <string> CA_EDGE_OR_FIX_
%token <string> CA_OR_FIX_
%token <string> CHARGE_NUMBER_
%token <string> CIRCLE_
%token <string> COEFFICIENT_
%token <string> COMPRESSIBLE_
%token <string> CONCENTRATION_
%token <string> CONCENTRATION_CONTOURS_
%token <string> CONDITION_
%token <string> CONDUCTION_
%token <string> CONDUCTIVITY_
%token <string> CONSTANT_
%token <string> CONSTANT_2_
%token <string> CONSTITUTIVE_
%token <string> CONTACT_LINE_
%token <string> CONTINUATIION_PRINTING_FREQUENCY_
%token <string> CONTINUATION_
%token <string> CONTINUITY_
%token <string> CONTOURS_ 
%token <string> CONT_NORM_VEL_
%token <string> CONT_TANG_VEL_
%token <string> CONVECTION_
%token <string> CONVECTIVE_
%token <string> COORDINATE_
%token <string> CORRECTION_
%token <string> CROSS_STREAM_SHEAR_RATE_
%token <string> CR_
%token <string> CURE_
%token <string> CURRENT_
%token <string> CURVATURE_DIFFUSIVITY_
%token <string> CYLINDRICAL_
%token <string> CONTINUOUS_
%token <string> D1_
%token <string> D1_RS_
%token <string> D2_
%token <string> D2_RS_
%token <string> D3_
%token <string> D3_RS_
%token <string> DARCEY_
%token <string> DARCEY_FICKIAN_
%token <string> DARCY_CONTINUOUS_
%token <string> DARCY_CONTINUOUS_
%token <string> DATABASE_
%token <string> DATA_
%token <string> DATA_SENS_
%token <string> DATUM_
%token <string> DEBUG_ 
%token <string> DEFAULT_
%token <string> DEFAULT_MATERIAL_SPECIES_TYPE_
%token <string> DEFINED_
%token <string> DEFORMABLE_
%token <string> DEFORM_
%token <string> DELTA_S_
%token <string> DELIMITER_
%token <string> DELTA_T_ 
%token <string> DENSITY_
%token <string> DERIVATIVES_ 
%token <string> DERIVATIVE_
%token <string> DESCRIPTION_
%token <string> DIFFUSION_
%token <string> DIFFUSIVITY_
%token <string> DILATATION_
%token <string> DILATION_
%token <string> DILATATION_
%token <string> DISABLE_
%token <string> DISCONTINUOUS_
%token <string> DISCONTINUOUS_VELO_
%token <string> DISTNG_
%token <string> DISTNG_
%token <string> DIVERGENCE_
%token <string> DOMAIN_
%token <string> DROP_
%token <string> DUMPS_
%token <string> DXDISTNG_
%token <string> DX_ 
%token <string> DX_RS_ 
%token <string> DYDISTNG_
%token <string> DY_
%token <string> DY_RS_
%token <string> DZDISTNG_
%token <string> DZ_
%token <string> DZ_RS_
%token <string> ELECTRONEUTRALITY_FICKIAN_	
%token <string> ELECTRONEUTRALITY_SM_	
%token <string> EDGE_
%token <string> ELECTRICAL_
%token <string> ELECTROOSMOTIC_
%token <string> ELEMENT_
%token <string> ELEM_
%token <string> END_
%token <string> ENERGY_
%token <string> ENERGY_CONDUCTION_VECTORS_
%token <string> ENERGY_FLUXLINES_
%token <string> ENTHALPY_
%token <string> EPOXY_
%token <string> EPSILON_
%token <string> EQUALS_ 
%token <string> EQUATIONS_
%token <string> EQUATION_
%token <string> EQ_ 
%token <string> ERROR_
%token <string> ERROR_ZZ_HEAT_FLUX_
%token <string> ERROR_ZZ_PRESSURE_
%token <string> ERROR_ZZ_VELOCITY_
%token <string> EVP_
%token <string> EVSS_F_
%token <string> EVSS_G_
%token <string> EXODUS_
%token <string> EXOII_FILENAME_
%token <string> EXPANSION_
%token <string> EXPONENT_
%token <string> EXPONENTIAL_
%token <string> EXTERNAL_
%token <string> FACTORIZATION_
%token <string> FACTOR_ 
%token <string> FEM_
%token <string> FICKIAN_
%token <string> FICKIAN_DIFFUSIVITY_
%token <string> FIELD_
%token <string> FILE_ 
%token <string> FILLED_EPOXY_
%token <string> FILL_
%token <string> FILL_CA_
%token <string> FILL_CONTOURS_
%token <string> FILL_INLET_
%token <string> FILL_INLET_
%token <string> FILL_SUBCYCLE_
%token <string> FILTER_CONCENTRATION_
%token <string> FINAL_PARAMETER_VALUE_
%token <string> FIRST_
%token <string> FIX_
%token <string> FLOAT_
%token <string> FLORY_HUGGINS_
%token <string> FLOWING_
%token <string> FLOW_GRADV_
%token <string> FLOW_HYDROSTATIC_
%token <string> FLOW_PRESSURE_
%token <string> FLOW_PRESS_USER_
%token <string> FLOW_STRESSNOBC_
%token <string> FLUCTUATIONS_ 
%token <string> FLUID_SOLID_
%token <string> FLUID_SOLID_RS_
%token <string> FLUXES_
%token <string> FLUXINESS_ 
%token <string> FLUXLINES_
%token <string> FLUX_
%token <string> FLUX_SENS_
%token <string> FOOBAR_
%token <string> FORCE_
%token <string> FORCE_RS_ 
%token <string> FORCE_USER_
%token <string> FORCE_USER_RS_
%token <string> FORMULATION_
%token <string> FRAC_
%token <string> FREE_
%token <string> FREE_VOL_
%token <string> FREQUENCY_
%token <string> FUNCTION_
%token <string> FUSION_
%token <string> F_
%token <string> G11_
%token <string> G12_
%token <string> G13_
%token <string> G21_
%token <string> G22_
%token <string> G23_
%token <string> G31_
%token <string> G32_
%token <string> G33_
%token <string> GALERKIN_
%token <string> GAS_
%token <string> GD_ 
%token <string> GD_CIRC_
%token <string> GD_CONST_
%token <string> GD_LINEAR_
%token <string> GD_PARAB_
%token <string> GD_POLYN_
%token <string> GD_TABLE_
%token <string> GD_TIME_
%token <string> GEL_
%token <string> GENERAL_
%token <string> GEOMX_
%token <string> GEOMY_
%token <string> GEOMZ_
%token <string> GEOM_
%token <string> GIESEKUS_
%token <string> GRADIENT11_
%token <string> GRADIENT12_
%token <string> GRADIENT13_
%token <string> GRADIENT21_
%token <string> GRADIENT22_
%token <string> GRADIENT23_
%token <string> GRADIENT31_
%token <string> GRADIENT32_
%token <string> GRADIENT33_
%token <string> GRAVITY_BASED_DIFFUSIVITY_
%token <string> GUESS_
%token <string> HC_
%token <string> HEAT_
%token <string> HIGH_
%token <string> HOOKEAN_PSTRAIN_
%token <string> HOOKEAN_PSTRESS_
%token <string> HUNTING_SPECIFICATION_
%token <string> HYDRODYNAMIC_
%token <string> HYDROSTATIC_SYMM_
%token <string> IDIM_ 
%token <string> II_ 
%token <string> INCOMP_PSTRAIN_
%token <string> INCOMP_PSTRESS_
%token <string> INERTIA_
%token <string> INITIALIZE_
%token <string> INITIAL_ 
%token <string> INITIAL_PARAMETER_VALUE_
%token <string> INITIAL_TIME_
%token <string> INTEGRATION_
%token <string> INTERMEDIATE_
%token <string> INVARIANT_
%token <string> IN_
%token <string> ISOPARAMETRIC_
%token <string> ITERATIONS_
%token <string> JACOBIAN_
%token <string> KINEMATIC_
%token <string> KINEMATIC_COLLOC_
%token <string> KINEMATIC_DISC_
%token <string> KINEMATIC_EDGE_
%token <string> KINEMATIC_PETROV_
%token <string> KINEMATIC_SPECIES_
%token <string> KIN_CHEM_
%token <string> KIN_DISPLACEMENT_
%token <string> KIN_LEAK_
%token <string> KRYLOV_
%token <string> LAGRANGIAN_
%token <string> LAMBDA_
%token <string> LAME_
%token <string> LATENT_
%token <string> LATENT_HEAT_
%token <string> LATENT_HEAT_INTERNAL_
%token <string> LAW_
%token <string> LEVEL_
%token <string> LINEAR_
%token <string> BILINEAR_
%token <string> LINEAR__STABILITY_
%token <string> STABILITY_
%token <string> LIQUIDUS_
%token <string> LIQUID_
%token <string> LIQ_
%token <string> LOW_
%token <string> MAPPING_ 
%token <string> MASS_
%token <string> MATERIALS_  
%token <string> MATERIAL_ID_
%token <string> MATERIAL_IS_NONDILUTE_
%token <string> MATERIAL_PROPERTY_TAG_
%token <string> MATERIAL_PROPERTY_TAG_SUBINDEX_
%token <string> MATRIX_
%token <string> MAT_
%token <string> MAXIMUM_
%token <string> MAXIMUM_NUMBER_OF_PATH_STEPS_
%token <string> MAXIMUM_PATH_STEP_
%token <string> MAXIMUM_PATH_VALUE_
%token <string> MAXIMUM_TIME_
%token <string> MAXIMUM_TIME_STEP_
%token <string> MEAN_
%token <string> MECHANICAL_PROPERTIES_DELIMITER_
%token <string> MEDIA_
%token <string> MESH1_
%token <string> MESH2_
%token <string> MESH3_
%token <string> MESH_
%token <string> MINIMUM_
%token <string> MINIMUM_PATH_STEP_
%token <string> MOBILITY_
%token <string> MODEL_
%token <string> MODIFIED_
%token <string> MOLAR_
%token <string> MOLECULAR_WEIGHT_
%token <string> MOMENTUM1_
%token <string> MOMENTUM2_
%token <string> MOMENTUM3_
%token <string> MOM_ 
%token <string> MOM_SOLID1_
%token <string> MOM_SOLID2_
%token <string> MOM_SOLID3_
%token <string> MOTION_ 
%token <string> MOVING_
%token <string> MOVING_CA_
%token <string> MOVING_CA_
%token <string> MOVING_PLANE_
%token <string> MOVING_PLANE_
%token <string> MU_
%token <string> N2_ 
%token <string> N3_ 
%token <string> NAME_
%token <string> NAVIER_
%token <string> NAVIER_STOKES_
%token <string> NETWORK_
%token <string> NEWTONIAN_
%token <string> NEWTON_
%token <string> NOLDROYDB_
%token <string> NONE_
%token <string> NONLINEAR_
%token <string> NONLINEAR_PLANE_STRAIN_
%token <string> INCOMP_3D_
%token <string> NON_CONDENSABLE_
%token <string> MOLECULAR_
%token <string> WEIGHT_
%token <string> NON_VOLATILE_
%token <string> NOPOLYMER_
%token <string> NORMALIZED_ 
%token <string> NORMALIZED_CORRECTION_TOLERANCE_
%token <string> NORMAL_
%token <string> NORMAL_AND_TANGENT_VECTORS_
%token <string> NORM_
%token <string> NORM_FORCE_
%token <string> NORM_FORCE_ 
%token <string> NORM_FORCE_RS_
%token <string> NORM_FORCE_RS_ 
%token <string> NO_
%token <string> NO_SLIP_
%token <string> NO_SLIP_RS_
%token <string> NUMBER_
%token <string> N_ 
%token <string> OF_
%token <string> OLDROYDB_
%token <string> ORDER_
%token <string> ORTHOGONALIZATION_ 
%token <string> OUTPUT_
%token <string> OFF_
%token <string> OVERLAP_
%token <string> P0_
%token <string> P1_
%token <string> P1_P2_Q1_OR_Q2_
%token <string> PACKING_
%token <string> PARAMETER_
%token <string> PARAMETERS_
%token <string> PARSER_IN_INPUT_FILE_MODE PARSER_IN_MAT_FILE_MODE
%token <string> PARTIALLY_WETTING_
%token <string> PARTICAL_VELOCITY_DIVERGENCE_
%token <string> PARTICLE_
%token <string> PATH_STEP_ERROR_
%token <string> PATH_STEP_PARAMETER_
%token <string> PENETRATION_
%token <string> PERMEABILITY_
%token <string> PHASE_
%token <string> PHYSICAL_
%token <string> PLANEX_
%token <string> PLANEY_ 
%token <string> PLANEZ_
%token <string> PLANE_
%token <string> PLASTIC_
%token <string> PMOMENTUM1_
%token <string> PMOMENTUM2_
%token <string> PMOMENTUM3_
%token <string> POINT_
%token <string> POLYMER_
%token <string> POLYNOMIAL_
%token <string> POROSITY_
%token <string> POROUS_
%token <string> POROUS_PART_SAT_
%token <string> POROUS_BRINKMAN_
%token <string> POROUS_CONV_
%token <string> POROUS_FLUX_
%token <string> POROUS_GAS_
%token <string> POROUS_KIN_
%token <string> POROUS_PRESSURE_
%token <string> POROUS_SATURATED_
%token <string> POROUS_TWO_PHASE_
%token <string> POROUS_UNSATURATED_
%token <string> POST_
%token <string> POWER_
%token <string> POWER_LAW_
%token <string> PQ1_
%token <string> PRECONDITIONER_
%token <string> PRECONDITIONS_ 
%token <string> PRESSURE_
%token <string> PRESSURE_USER_
%token <string> PRINTING_
%token <string> PROBLEM_
%token <string> PROCESSING_
%token <string> PROCESSORS_ 
%token <string> PROPERTIES_
%token <string> PSD_SEXP_
%token <string> PSD_VOL_
%token <string> PSD_WEXP_
%token <string> PSEUDO_SOLID_
%token <string> PSPG_
%token <string> PSPG_
%token <string> PTT_
%token <string> PU_
%token <string> PV_
%token <string> PW_
%token <string> P_
%token <string> P_EQUIL_
%token <string> Q1_
%token <string> Q1_D_
%token <string> Q1_OR_Q2_
%token <string> Q2_
%token <string> Q2_D_
%token <string> QCONV_
%token <string> QRAD_
%token <string> QSIDE_
%token <string> QUAD_GP_
%token <string> QUSER_
%token <string> Q_ANY_
%token <string> Q_TENSOR_DIFFUSIVITY_
%token <string> Q_VELO_SLIP_
%token <string> RATE_
%token <string> RATIO_
%token <string> REACTION_RATE_
%token <string> REAL_
%token <string> REFERENCE_
%token <string> REFORM_
%token <string> REL_
%token <string> REP_FORCE_
%token <string> REP_FORCE_RS_
%token <string> RESIDUALS_ 
%token <string> RESIDUAL_
%token <string> RESULTS_
%token <string> REUSE_
%token <string> RECALC_
%token <string> CALC_
%token <string> ROTATIONAL_
%token <string> ROTATIONS_
%token <string> ROT_
%token <string> S11_
%token <string> S11_1_
%token <string> S11_2_
%token <string> S11_3_
%token <string> S11_4_
%token <string> S11_5_
%token <string> S11_6_
%token <string> S11_7_
%token <string> S12_
%token <string> S12_1_
%token <string> S12_2_
%token <string> S12_3_
%token <string> S12_4_
%token <string> S12_5_
%token <string> S12_6_
%token <string> S12_7_
%token <string> S13_
%token <string> S13_1_
%token <string> S13_2_
%token <string> S13_3_
%token <string> S13_4_
%token <string> S13_5_
%token <string> S13_6_
%token <string> S13_7_
%token <string> S21_
%token <string> S22_
%token <string> S22_1_
%token <string> S22_2_
%token <string> S22_3_
%token <string> S22_4_
%token <string> S22_5_
%token <string> S22_6_
%token <string> S22_7_
%token <string> S23_
%token <string> S23_1_
%token <string> S23_2_
%token <string> S23_3_
%token <string> S23_4_
%token <string> S23_5_
%token <string> S23_6_
%token <string> S23_7_
%token <string> S33_
%token <string> S33_1_
%token <string> S33_2_
%token <string> S33_3_
%token <string> S33_4_
%token <string> S33_5_
%token <string> S33_6_
%token <string> S33_7_
%token <string> SATURATION_
%token <string> SCALING_
%token <string> SECONDARY_FREQUENCY_
%token <string> SECOND_ 
%token <string> SECOND_FREQUENCY_TIME_
%token <string> SEED_ 
%token <string> SENSITIVITIES_
%token <string> SHEAR_
%token <string> SHEAR_RATE_
%token <string> SHEAR_RATE_DIFFUSIVITY_
%token <string> SH_
%token <string> SIZE_
%token <string> SINUSOIDAL_
%token <string> SLOPEX_
%token <string> SLOPEY_
%token <string> SLOPEZ_
%token <string> SLOPE_ 
%token <string> SOLIDUS_
%token <string> SOLID_
%token <string> SOLID_FLUID_
%token <string> SOLID_FLUID_RS_
%token <string> SOLN_
%token <string> SOLUTION_
%token <string> SOLUTION_TEMPERATURE_
%token <string> SOLVENT_
%token <string> SOLVER_
%token <string> SUBDOMAIN_
%token <string> SOLVE_
%token <string> SOURCE_
%token <string> SPECIES_
%token <string> SPECIES_BULK_
%token <string> SPECIES_SURF_
%token <string> SPECIES_TIME_INTEGRATION_
%token <string> SPECIES_UNKNOWN_
%token <string> SPECIES_WEIGHT_FUNCTION_
%token <string> SPECIFICATIONS_
%token <string> SPECIFIC_
%token <string> SPHERICAL_
%token <string> SPLINEX_
%token <string> SPLINEY_ 
%token <string> SPLINEZ_
%token <string> SPLINE_
%token <string> SP_
%token <string> SP_N_
%token <string> STABILIZATION_
%token <string> STAGES_ 
%token <string> STEPS_
%token <string> STEP_ 
%token <string> STOKES_
%token <string> STRAIN_
%token <string> STREAMWISE_
%token <string> STREAM_
%token <string> STRESS11_
%token <string> STRESS12_
%token <string> STRESS13_
%token <string> STRESS22_
%token <string> STRESS23_
%token <string> STRESS33_
%token <string> STRESS_
%token <string> STRESS_CONTOURS_
%token <string> STRIDE_
%token <string> STRING_
%token <string> STRONG_FILL_CA_
%token <string> SUBCYCLE_
%token <string> SUBPARAMETRIC_
%token <string> SUBSPACE_ 
%token <string> SUM_TO_ONE_
%token <string> SUPG_
%token <string> SURFACE_
%token <string> SURFACES_
%token <string> SURFDOMAIN_ENERGY_BAL_
%token <string> SURFDOMAIN_MASS_FRACTION_
%token <string> SURFDOMAIN_STEFAN_FLOW_
%token <string> SURFTANG_
%token <string> SURFTANG_EDGE_
%token <string> SURFTANG_SCALAR_
%token <string> SURFTANG_SCALAR_EDGE_
%token <string> SUSPENSION_
%token <string> SWIRLING_
%token <string> SYSTEM_ 
%token <string> S_
%token <string> T1_ 
%token <string> T2_ 
%token <string> TABLE_
%token <string> TABLE_WICS_
%token <string> TABLE_WICV_
%token <string> TANGENT_
%token <string> DOM_DECOMP_
%token <string> TEMPERATURE_
%token <string> TENSION_
%token <string> TENSOR_
%token <string> THERMAL_
%token <string> THIRD_ 
%token <string> TIME_
%token <string> TNRMLSIDE_
%token <string> TOLERANCE_ 
%token <string> TOTAL_ALE_
%token <string> TOTAL_VELOCITY_DIVERGENCE_
%token <string> TRACES_
%token <string> TSHRSIDE_
%token <string> TYPE_ 
%token <string> T_
%token <string> T_MELT_
%token <string> T_USER_
%token <string> VOLT_USER_
%token <string> U1_
%token <string> U2_
%token <string> U3_
%token <string> UMF_
%token <string> UMFF_
%token <string> UNREACTED_
%token <string> UNSEEDED_ 
%token <string> USER_
%token <string> USER_DEFINED_
%token <string> USER_FILL_
%token <string> USER_GEN_
%token <string> UUSER_
%token <string> UVARY_ 
%token <string> U_ 
%token <string> VAN_GENUCHTEN_
%token <string> VAPORIZATION_
%token <string> VAPOR_
%token <string> VAR_CA_EDGE_
%token <string> VAR_CA_USER_
%token <string> VECTORS_ 
%token <string> VECTOR_ 
%token <string> VELOCITY_
%token <string> VELOCITY_DIVERGENCE_
%token <string> VELO_NORMAL_
%token <string> VELO_NORMAL_DISC_
%token <string> VELO_NORMAL_EDGE_
%token <string> VELO_NORMAL_EDGE_INT_
%token <string> VELO_NORM_COLLOC_
%token <string> VELO_SLIP_
%token <string> VELO_SLIP_FILL_
%token <string> VELO_SLIP_ROT_
%token <string> VELO_SLIP_SOLID_
%token <string> VELO_TANGENT_
%token <string> VELO_TANGENT_3D_
%token <string> VELO_TANGENT_EDGE_
%token <string> VELO_TANGENT_EDGE_INT_
%token <string> VELO_TANGENT_SOLID_
%token <string> VERTEX_
%token <string> VISCOPLASTIC_DEF_GRAD_TENSOR_
%token <string> VISCOSITY_
%token <string> VISCOSITY_DIFFUSIVITY_
%token <string> VL_EQUIL_
%token <string> VL_POLY_
%token <string> VNORM_LEAK_
%token <string> VN_POROUS_
%token <string> VOLTAGE_
%token <string> VOLT_
%token <string> VOLUMETRIC_
%token <string> VOLUME_
%token <string> VOLUME_INT_
%token <string> VOL_
%token <string> VORTICITY_VECTOR_
%token <string> VP_EQUIL_
%token <string> VUSER_
%token <string> VVARY_
%token <string> V_ 
%token <string> WEIGHTING_
%token <string> WEIGHT_
%token <string> WETTING_
%token <string> WHITE_METZNER_
%token <string> WLF_
%token <string> WRITE_
%token <string> WRITE_INITIAL_SOLUTION_
%token <string> WUSER_
%token <string> WVARY_ 
%token <string> W_
%token <string> XDIM_
%token <string> XI_
%token <string> X_ 
%token <string> Y2_ELECTRONEUTRALITY_
%token <string> Y2_ELECTRONEUTRALITY_   
%token <string> YEILD_
%token <string> YES_ 
%token <string> YFLUX_
%token <string> YFLUX_BV_
%token <string> YFLUX_HOR_
%token <string> YFLUX_ORR_
%token <string> YFLUX_H2O_ANODE_
%token <string> YFLUX_H2O_CATHODE_
%token <string> YFLUX_CONST_
%token <string> YFLUX_EQUIL_
%token <string> YFLUX_SUS_
%token <string> YFLUX_USER_
%token <string> YTOTALFLUX_CONST_
%token <string> YUSER_
%token <string> Y_
%token <string> Y_DISCONTINUOUS_
%token <string> YZ_
%token <string> YX_
%token <string> Z_
%token <string> NS_
%token <string> SS_
%token <string> NC_
%token <string> SC_
%token <string> COLLOCATE_EDGE_ 	
%token <string> WEAK_INT_EDGE_
%token <string> STRONG_INT_EDGE_
%token <string> RAOULT_
%token <string> FLORY_
%token <string> LU_
%token <string> FRONT_
%token <string> MA28_
%token <string> Y12M_
%token <string> GMRES_
%token <string> CG_ 
%token <string> CGS_
%token <string> TFQMR_
%token <string> BICGSTAB_
%token <string> JACOBI_
%token <string> BJACOBI_
%token <string> LS_
%token <string> ILU_
%token <string> RILU_
%token <string> ILUT_
%token <string> ICC_
%token <string> BILU_
%token <string> ROW_SYM_
%token <string> SYM_DIAG_
%token <string> SUSPENSION_PM_
%token <string> CONSTANT_LAST_CONC_
%token <string> FOAM_
%token <string> LEVEL_SET_
%token <string> IDEAL_GAS_
%token <string> THERMAL_BATTERY_
%token <string> SYM_ROW_SUM_
%token <string> SYM_GS_
%token <string> NEUMANN_
%token <string> STORAGE_
%token <string> VBR_
%token <string> MSR_
%token <string> FORMAT_
%token <string> AUGMENTING_
%token <string> CONDITIONS_
%token <string> ZY_
%token <string> ZX_
%token <string> R_MOMENTUM1_  	 
%token <string> VX_  		
%token <string> R_MOMENTUM2_  	 
%token <string> VY_  		
%token <string> R_MOMENTUM3_  	 
%token <string> VZ_  		
%token <string> R_ENERGY_  	
%token <string> T_  		
%token <string> R_MASS_  		 
%token <string> Y_  		
%token <string> R_MESH1_  	 
%token <string> DMX_  		
%token <string> R_MESH2_  	 
%token <string> DMY_  		
%token <string> R_MESH3_  	 
%token <string> DMZ_  		
%token <string> R_MASS_SURF_  	 
%token <string> S_  		
%token <string> R_PRESSURE_  	 
%token <string> P_  		
%token <string> R_STRESS11_  	 
%token <string> S11_  		
%token <string> R_STRESS12_  	 
%token <string> S12_  		
%token <string> R_STRESS22_  	 
%token <string> S22_  		
%token <string> R_STRESS13_  	 
%token <string> S13_  		
%token <string> R_STRESS23_ 	 
%token <string> S23_  		
%token <string> R_STRESS33_  	 
%token <string> S33_  		
%token <string> R_SOLID1_  	 
%token <string> DMX_  		
%token <string> R_SOLID2_  	 
%token <string> DMY_RS_  		
%token <string> R_SOLID3_  	 
%token <string> DMZ_RS_  		
%token <string> R_GRADIENT11_  	 
%token <string> G11_  		
%token <string> R_GRADIENT12_  	 
%token <string> G12_  		
%token <string> R_GRADIENT21_  	 
%token <string> G21_  		
%token <string> R_GRADIENT22_  	 
%token <string> G22_  		
%token <string> R_GRADIENT13_  	 
%token <string> G13_  		
%token <string> R_GRADIENT23_  	 
%token <string> G23_  		
%token <string> R_GRADIENT31_  	 
%token <string> G31_  		
%token <string> R_GRADIENT32_  	 
%token <string> G32_  		
%token <string> R_GRADIENT33_  	 
%token <string> G33_  		
%token <string> R_POTENTIAL_  	 
%token <string> V_  		 
%token <string> R_FILL_  		 
%token <string> F_  		 
%token <string> R_SHEAR_RATE_  	 
%token <string> SH_  		
%token <string> R_PMOMENTUM1_  	 
%token <string> PVX_  		
%token <string> R_PMOMENTUM2_  	 
%token <string> PVY_  		
%token <string> R_PMOMENTUM3_  	 
%token <string> PVZ_  		
%token <string> R_STRESS11_1_  	 
%token <string> S11_1_  		
%token <string> R_STRESS12_1_  	 
%token <string> S12_1_  		
%token <string> R_STRESS22_1_  	 
%token <string> S22_1_  		
%token <string> R_STRESS13_1_  	 
%token <string> S13_1_  		 
%token <string> R_STRESS23_1_  	 
%token <string> S23_1_  		
%token <string> R_STRESS33_1_  	 
%token <string> S33_1_  		 
%token <string> R_STRESS11_2_  	 
%token <string> S11_2_  		
%token <string> R_STRESS12_2_  	 
%token <string> S12_2_  		
%token <string> R_STRESS22_2_  	 
%token <string> S22_2_  		 
%token <string> R_STRESS13_2_  	 
%token <string> S13_2_  		
%token <string> R_STRESS23_2_  	 
%token <string> S23_2_  		
%token <string> R_STRESS33_2_  	 
%token <string> S33_2_  		 
%token <string> R_STRESS11_3_  	 
%token <string> S11_3_  		
%token <string> R_STRESS12_3_  	 
%token <string> S12_3_  		
%token <string> R_STRESS22_3_  	 
%token <string> S22_3_  		 
%token <string> R_STRESS13_3_  	 
%token <string> S13_3_  		
%token <string> R_STRESS23_3_  	 
%token <string> S23_3_  		
%token <string> R_STRESS33_3_  	 
%token <string> S33_3_  		 
%token <string> R_STRESS11_4_  	 
%token <string> S11_4_  		
%token <string> R_STRESS12_4_  	 
%token <string> S12_4_  		
%token <string> R_STRESS22_4_  	 
%token <string> S22_4_  		 
%token <string> R_STRESS13_4_  	 
%token <string> S13_4_  		
%token <string> R_STRESS23_4_  	 
%token <string> S23_4_  		
%token <string> R_STRESS33_4_  	 
%token <string> S33_4_  		 
%token <string> R_STRESS11_5_  	 
%token <string> S11_5_  		
%token <string> R_STRESS12_5_  	 
%token <string> S12_5_  		
%token <string> R_STRESS22_5_  	 
%token <string> S22_5_  		
%token <string> R_STRESS13_5_  	 
%token <string> S13_5_  		
%token <string> R_STRESS23_5_  	 
%token <string> S23_5_  		
%token <string> R_STRESS33_5_  	 
%token <string> S33_5_  		
%token <string> R_STRESS11_6_  	 
%token <string> S11_6_  		
%token <string> R_STRESS12_6_  	 
%token <string> S12_6_  		
%token <string> R_STRESS22_6_  	 
%token <string> S22_6_  		
%token <string> R_STRESS13_6_  	 
%token <string> S13_6_  		
%token <string> R_STRESS23_6_  	 
%token <string> S23_6_  		
%token <string> R_STRESS33_6_  	 
%token <string> S33_6_  		
%token <string> R_STRESS11_7_  	 
%token <string> S11_7_  		
%token <string> R_STRESS12_7_  	 
%token <string> S12_7_  		
%token <string> R_STRESS22_7_  	 
%token <string> S22_7_  		
%token <string> R_STRESS13_7_  	 
%token <string> S13_7_  		
%token <string> R_STRESS23_7_  	 
%token <string> S23_7_  		
%token <string> R_STRESS33_7_  	 
%token <string> S33_7_  		
%token <string> R_SPECIES_0_  	 
%token <string> SP_0_  		
%token <string> R_SPECIES_1_  	 
%token <string> SP_1_  		
%token <string> R_SPECIES_2_  	 
%token <string> SP_2_  		
%token <string> R_SPECIES_3_  	 
%token <string> SP_3_  		
%token <string> R_SPECIES_4_  	 
%token <string> SP_4_  		
%token <string> R_SPECIES_5_  	 
%token <string> SP_5_  		
%token <string> R_SPECIES_6_  	 
%token <string> SP_6_  		
%token <string> R_SPECIES_7_  	 
%token <string> SP_7_  		
%token <string> R_SPECIES_8_  	 
%token <string> SP_8_  		
%token <string> R_SPECIES_9_  	 
%token <string> SP_9_  		
%token <string> R_SPECIES_10_  	 
%token <string> SP_10_  		
%token <string> R_SPECIES_11_  	 
%token <string> SP_11_  		
%token <string> R_SPECIES_12_  	 
%token <string> SP_12_  		
%token <string> R_SPECIES_13_  	 
%token <string> SP_13_  		
%token <string> R_SPECIES_14_  	 
%token <string> SP_14_  		
%token <string> R_SPECIES_15_  	 
%token <string> SP_15_  		
%token <string> R_SPECIES_16_  	 
%token <string> SP_16_  		
%token <string> R_SPECIES_17_  	 
%token <string> SP_17_  		
%token <string> R_SPECIES_18_  	 
%token <string> SP_18_  		
%token <string> R_SPECIES_19_  	 
%token <string> SP_19_  		
%token <string> R_SPECIES_20_  	 
%token <string> SP_20_  		
%token <string> R_SPECIES_21_  	 
%token <string> SP_21_  		 
%token <string> R_SPECIES_22_  	 
%token <string> SP_22_  		 
%token <string> R_SPECIES_23_  	 
%token <string> SP_23_ 		 
%token <string> R_SPECIES_24_  	 
%token <string> SP_24_  		
%token <string> R_SPECIES_25_  	 
%token <string> SP_25_  		
%token <string> R_SPECIES_26_  	 
%token <string> SP_26_  		
%token <string> R_SPECIES_27_  	 
%token <string> SP_27_  		
%token <string> R_SPECIES_28_  	 
%token <string> SP_28_  		
%token <string> R_SPECIES_29_  	 
%token <string> SP_29_  		
%token <string> R_VOLFRACPH_0_  	 
%token <string> VFP_0_  		
%token <string> R_VOLFRACPH_1_  	 
%token <string> VFP_1_  		
%token <string> R_VOLFRACPH_2_  	 
%token <string> VFP_2_  		
%token <string> R_VOLFRACPH_3_  	 
%token <string> VFP_3_  		
%token <string> R_VOLFRACPH_4_  	 
%token <string> VFP_4_  		
%token <string> R_Y0_  		 
%token <string> Y0_  		
%token <string> R_Y1_  		 
%token <string> Y1_  		
%token <string> R_Y2_  		 
%token <string> Y2_  		
%token <string> R_Y3_  		 
%token <string> Y3_  		 
%token <string> R_Y4_  		 
%token <string> Y4_  		
%token <string> R_Y5_  		 
%token <string> Y5_  		
%token <string> R_Y6_  		 
%token <string> Y6_  		
%token <string> R_Y7_  		 
%token <string> Y7_  		
%token <string> R_Y8_  		 
%token <string> Y8_  		
%token <string> R_Y9_  		 
%token <string> Y9_  		
%token <string> R_Y10_  		 
%token <string> Y10_  		 
%token <string> R_Y11_  		 
%token <string> Y11_  		 
%token <string> R_Y12_  		 
%token <string> Y12_  		 
%token <string> R_Y13_  		 
%token <string> Y13_  		 
%token <string> R_Y14_  		 
%token <string> Y14_  		 
%token <string> R_Y15_  		 
%token <string> Y15_  		 
%token <string> R_Y16_  		 
%token <string> Y16_  		 
%token <string> R_Y17_  		 
%token <string> Y17_  		 
%token <string> R_Y18_  		 
%token <string> Y18_  		 
%token <string> R_Y19_  		 
%token <string> Y19_  		 
%token <string> R_Y20_  		 
%token <string> Y20_  		 
%token <string> R_Y21_  		 
%token <string> Y21_  		 
%token <string> R_Y22_  		 
%token <string> Y22_  		 
%token <string> R_Y23_  		 
%token <string> Y23_  		 
%token <string> R_Y24_  		 
%token <string> Y24_  		 
%token <string> R_Y25_  		 
%token <string> Y25_  		 
%token <string> R_Y26_  		 
%token <string> Y26_  		 
%token <string> R_Y27_  		 
%token <string> Y27_  		 
%token <string> R_Y28_  		 
%token <string> Y28_  		 
%token <string> R_Y29_  		 
%token <string> Y29_  		 
%token <string> R_MOM_NORMAL_  	 
%token <string> DN_  		 
%token <string> R_MOM_TANG1_  	 
%token <string> DT1_  		 
%token <string> R_MOM_TANG2_  	 
%token <string> DT2_  		 
%token <string> R_MESH_NORMAL_  	 
%token <string> VN_  		 
%token <string> R_MESH_TANG1_  	 
%token <string> VT1_  		 
%token <string> R_MESH_TANG2_  	 
%token <string> VT2_  		 
%token <string> GEOM_
%token <string> GEOMX_
%token <string> GEOMY_
%token <string> GEOMZ_
%token <string> W_
%token <string> SURFCOMAIN_MASS_FRACTION_
%token <string> SURFCOMAIN_ENERGY_BAL_FRACTION_
%token <string> SURFCOMAIN_STEFAN_FLOW_
%token <string> VELOCITY1_
%token <string> VELOCITY2_
%token <string> VELOCITY3_
%token <string> MASS_FRACTION_
%token <string> MESH_DISPLACEMENT1_
%token <string> MESH_DISPLACEMENT2_
%token <string> MESH_DISPLACEMENT3_
%token <string> POLYMER_STRESS11_
%token <string> POLYMER_STRESS12_
%token <string> POLYMER_STRESS13_
%token <string> POLYMER_STRESS22_
%token <string> POLYMER_STRESS23_
%token <string> POLYMER_STRESS33_
%token <string> SOLID_DISPLACEMENT1_
%token <string> DMX_RS_
%token <string> SOLID_DISPLACEMENT2_
%token <string> SOLID_DISPLACEMENT3_
%token <string> VELOCITY_GRADIENT11_
%token <string> VELOCITY_GRADIENT12_
%token <string> VELOCITY_GRADIENT22_
%token <string> VELOCITY_GRADIENT13_
%token <string> VELOCITY_GRADIENT23_
%token <string> VELOCITY_GRADIENT31_
%token <string> VELOCITY_GRADIENT32_
%token <string> VELOCITY_GRADIENT33_
%token <string> PVELOCITY1_
%token <string> PVELOCITY2_
%token <string> PVELOCITY3_
%token <string> POLYMER_STRESS11_1_
%token <string> POLYMER_STRESS12_1_
%token <string> POLYMER_STRESS22_1_
%token <string> POLYMER_STRESS13_1_
%token <string> POLYMER_STRESS23_1_
%token <string> POLYMER_STRESS33_1_
%token <string> POLYMER_STRESS11_2_
%token <string> POLYMER_STRESS12_2_
%token <string> POLYMER_STRESS13_2_
%token <string> POLYMER_STRESS23_2_
%token <string> POLYMER_STRESS33_2_
%token <string> POLYMER_STRESS11_3_
%token <string> POLYMER_STRESS12_3_
%token <string> POLYMER_STRESS22_3_
%token <string> POLYMER_STRESS13_3_
%token <string> POLYMER_STRESS23_3_
%token <string> POLYMER_STRESS33_3_
%token <string> POLYMER_STRESS11_4_
%token <string> POLYMER_STRESS12_4_
%token <string> POLYMER_STRESS13_4_
%token <string> POLYMER_STRESS23_4_
%token <string> POLYMER_STRESS33_4_
%token <string> POLYMER_STRESS11_5_
%token <string> POLYMER_STRESS12_5_
%token <string> POLYMER_STRESS22_5_
%token <string> POLYMER_STRESS13_5_
%token <string> POLYMER_STRESS23_5_
%token <string> POLYMER_STRESS33_5_
%token <string> POLYMER_STRESS11_6_
%token <string> POLYMER_STRESS12_6_
%token <string> POLYMER_STRESS22_6_
%token <string> POLYMER_STRESS13_6_
%token <string> POLYMER_STRESS23_6_
%token <string> POLYMER_STRESS33_6_
%token <string> POLYMER_STRESS11_7_
%token <string> POLYMER_STRESS12_7_
%token <string> POLYMER_STRESS22_7_
%token <string> POLYMER_STRESS13_7_
%token <string> POLYMER_STRESS23_7_
%token <string> POLYMER_STRESS33_7_
%token <string> SPECIES_CONC_0_
%token <string> SPECIES_CONC_1_
%token <string> SPECIES_CONC_2_
%token <string> SPECIES_CONC_3_
%token <string> SPECIES_CONC_4_
%token <string> SPECIES_CONC_5_
%token <string> SPECIES_CONC_6_
%token <string> SPECIES_CONC_7_
%token <string> SPECIES_CONC_8_
%token <string> SPECIES_CONC_9_
%token <string> SPECIES_CONC_10_
%token <string> SPECIES_CONC_11_
%token <string> SPECIES_CONC_12_
%token <string> SPECIES_CONC_13_
%token <string> SPECIES_CONC_14_
%token <string> SPECIES_CONC_15_
%token <string> SPECIES_CONC_16_
%token <string> SPECIES_CONC_17_
%token <string> SPECIES_CONC_18_
%token <string> SPECIES_CONC_19_
%token <string> SPECIES_CONC_20_
%token <string> SPECIES_CONC_21_
%token <string> SPECIES_CONC_22_
%token <string> SPECIES_CONC_23_
%token <string> SPECIES_CONC_24_
%token <string> SPECIES_CONC_25_
%token <string> SPECIES_CONC_26_
%token <string> SPECIES_CONC_27_
%token <string> SPECIES_CONC_28_
%token <string> SPECIES_CONC_29_
%token <string> VOLFRACPH_0_
%token <string> VOLFRACPH_1_
%token <string> VOLFRACPH_2_
%token <string> VOLFRACPH_3_
%token <string> VOLFRACPH_4_
%token <string> MESH_POSITION1_
%token <string> MESH_POSITION2_
%token <string> MESH_POSITION3_
%token <string> VEL_NORM_
%token <string> D_VEL1_DT_
%token <string> UDOT_
%token <string> D_VEL2_DT_
%token <string> VDOT_
%token <string> D_VEL3_DT_
%token <string> WDOT_
%token <string> D_T_DT_
%token <string> TDOT_
%token <string> D_C_DT_
%token <string> TDOT_
%token <string> D_C_DT_
%token <string> CDOT_
%token <string> D_X1_DT_
%token <string> XDOT_
%token <string> D_X2_DT_
%token <string> VDOT_
%token <string> D_X3_DT_
%token <string> ZDOT_
%token <string> D_S_DT_
%token <string> SDOT_
%token <string> D_P_DT_
%token <string> PDOT_
%token <string> VELOCITY_GRADIENT21_
%token <string> VELOCITY_GRADIENT22_
%token <string> VELOCITY_GRADIENT13_
%token <string> POLYMER_STRESS22_2_
%token <string> POLYMER_STRESS22_4_
%token <string> POLYMER_STRESS23_6_
%token <string> YDOT_
%token <string> QUADRATIC_
%token <string> BIQUADRATIC_
%token <string> XY_
%token <string> XZ_
%token <string> LIQUID_VAPOR_
%token <string> SOLID_LIQUID_
%token <string> EVAPORATION_
%token <string> DISSOLUTION_
%token <string> ZERO_
%token <string> RANDOM_ 
%token <string> ONE_
%token <string> SINES_
%token <string> READ_	
%token <string> READ_EXOII_
%token <string> READ_EXOII_FILE_ 
%token <string> SPECIES_UNK_0_        
%token <string> SPECIES_UNK_1_        
%token <string> SPECIES_UNK_2_       
%token <string> SPECIES_UNK_3_       
%token <string> SPECIES_UNK_4_      
%token <string> SPECIES_UNK_5_  
%token <string> SPECIES_UNK_6_     
%token <string> SPECIES_UNK_7_       
%token <string> SPECIES_UNK_8_       
%token <string> SPECIES_UNK_9_       
%token <string> SPECIES_UNK_10_      
%token <string> SPECIES_UNK_11_    
%token <string> SPECIES_UNK_12_      
%token <string> SPECIES_UNK_13_      
%token <string> SPECIES_UNK_14_      
%token <string> SPECIES_UNK_15_      
%token <string> SPECIES_UNK_16_      
%token <string> SPECIES_UNK_17_      
%token <string> SPECIES_UNK_18_      
%token <string> SPECIES_UNK_19_      
%token <string> SPECIES_UNK_20_      
%token <string> SPECIES_UNK_21_      
%token <string> SPECIES_UNK_22_      
%token <string> SPECIES_UNK_23_      
%token <string> SPECIES_UNK_24_      
%token <string> SPECIES_UNK_25_       
%token <string> SPECIES_UNK_26_      
%token <string> SPECIES_UNK_27_       
%token <string> SPECIES_UNK_28_      
%token <string> SPECIES_UNK_29_      
%token <string> SPECIES_UNK_LAST_     
%token <string> VOLF_PHASE_0_         
%token <string> VOLF_PHASE_1_        
%token <string> VOLF_PHASE_2_        
%token <string> VOLF_PHASE_3_        
%token <string> VOLF_PHASE_4_         
%token <string> VOLF_PHASE_LAST_      
%token <string> POR_LIQ_PRES_	
%token <string> POR_GAS_PRES_	
%token <string> POR_POROSITY_		
%token <string> POR_SATURATION_
%token <string> POR_LAST_            
%token <string> VORT_DIR1_           
%token <string> VORT_DIR2_           
%token <string> VORT_DIR3_            
%token <string> VORT_LAMBDA_
%token <string> Q2_LSA_
%token <string> Q2_D_LSA_
%token <string> PQ2_
%token <string> ATM_
%token <string> TORR_
%token <string> ANNEAL_
%token <string> ON_
%token <string> STEADY_
%token <string> TRANSIENT_
%token <string> SET_
%token <string> INTERFACE_
%token <string> TRACKING_
%token <string> LENGTH_
%token <string> SCALE_
%token <string> IMPLICIT_
%token <string> EMBEDDED_
%token <string> INITIALIZATION_
%token <string> METHOD_
%token <string> PROJECT_
%token <string> EXO_READ_
%token <string> SURFACES_
%token <string> CONTROL_
%token <string> WIDTH_
%token <string> RENORMALIZATION_
%token <string> HUYGENS_
%token <string> HUYGENS_CONSTRAINED_
%token <string> SPHERE_
%token <string> PROJECTION_
%token <string> NODESET_
%token <string> SURF_
%token <string> MT_
%token <string> HUNTING_
%token <string> CONTINUATION_
%token <string> ID_
%token <string> DATA_
%token <string> THE_WORD_FLOAT_
%token <string> SUBINDEX_
%token <string> VALUE_
%token <string> FINAL_
%token <string> DELTA_S_
%token <string> PATH_
%token <string> LOCA_
%token <string> AGGRESSIVENESS_
%token <string> ALC_
%token <string> DESIRED_
%token <string> MAX_
%token <string> LIMIT_
%token <string> CC_
%token <string> TP_
%token <string> MATERIAL_
%token <string> TAG_
%token <string> PROPERTY_
%token <string> FRACTION_
%token <string> SENSITIVITY_
%token <string> TC_
%token <string> HZERO_
%token <string> HFIRST_
%token <string> PF_
%token <string> USERBC_
%token <string> USERMAT_
%token <string> VC_
%token <string> FC_
%token <string> VFLUX_
%token <string> FORCE_X_
%token <string> FORCE_Y_
%token <string> FORCE_Z_
%token <string> FORCE_NORMAL_
%token <string> FORCE_TANGENT1_
%token <string> FORCE_TANGENT2_
%token <string> VOLUME_FLUX_
%token <string> HEAT_FLUX_
%token <string> SPECIES_FLUX_
%token <string> R0_
%token <string> RHS_
%token <string> ANORM_
%token <string> NOSCALED_
%token <string> SOL_
%token <string> WEIGHTED_
%token <string> NORM_
%token <string> WARNINGS_
%token <string> ALL_
%token <string> LAST_
%token <string> GRAPH_
%token <string> FILLIN_
%token <string> DIAG_
%token <string> FULL_
%token <string> SYMMETRIC_
%token <string> STANDARD_
%token <string> RESID_
%token <string> RAND_
%token <string> REORDER_
%token <string> RCM_
%token <string> SAVE_
%token <string> RELAX_
%token <string> THRESHOLD_
%token <string> ABSOLUTE_
%token <string> RELATIVE_
%token <string> CLASSIC_
%token <string> THREE_D_
%token <string> THREE_D_FILE_
%token <string> INLINE_
%token <string> FILTER_
%token <string> EIGEN_
%token <string> EIGENSOLVER_
%token <string> MODES_
%token <string> RECORD_
%token <string> SIZE_
%token <string> STEPS_
%token <string> RECYCLE_
%token <string> SHIFTS_
%token <string> WAVE_
%token <string> NUMBERS_
%token <string> TOTAL_
%token <string> DIFFUSE_
%token <string> VISCOPLASTIC_
%token <string> DEF_
%token <string> GRAD_
%token <string> ZZ_
%token <string> SOLVENTS_
%token <string> DARCY_
%token <string> GRID_
%token <string> PECLET_
%token <string> VORTICITY_
%token <string> CROSSSTREAM_
%token <string> AVERAGE_CONC_
%token <string> SURF_DISSIP_
%token <string> AREA_
%token <string> VOL_REVOLUTION_
%token <string> TORQUE_
%token <string> PROFILE_
%token <string> SOLID_POSITION1_		     
%token <string> X_RS_ 			     
%token <string> SOLID_POSITION2_		     
%token <string> Y_RS_ 			    
%token <string> SOLID_POSITION3_		     
%token <string> Z_RS_ 	
%token <string> P_LIQ_		         
%token <string> VDX_			
%token <string> VDY_		
%token <string> VDZ_			
%token <string> VLAMBDA_         	
%token <string> P_SAT_
%token <string> P_POR_
%token <string> P_GAS_
%token <string> POR_PORSITY_
%token <string> USER_POST_
%token <string> ERROR_ZZ_VEL_
%token <string> ERROR_ZZ_Q_
%token <string> ERROR_ZZ_P_
%token <string> SURFACE_VECTORS_
%token <string> LAGRANGE_CONVECTION_
%token <string> STRAIN_TENSOR_
%token <string> REAL_STRESS_TENSOR_
%token <string> STRESS_TENSOR_
%token <string> TIME_DERIVATIVES_
%token <string> CONDUCTION_VECTORS_
%token <string> DIFFUSION_VECTORS_
%token <string> MM_RESIDUALS_
%token <string> NS_RESIDUALS_
%token <string> DIV_VELOCITY_
%token <string> THIRD_INVAR_STRAIN_
%token <string> SEC_INVAR_STRAIN_
%token <string> FIRST_INVAR_STRAIN_
%token <string> FILL_CONT_
%token <string> PRESSURE_CONT_
%token <string> MEAN_SHEAR_
%token <string> CROSS_STREAM_SHEAR_
%token <string> STREAM_NORMAL_STRESS_
%token <string> SECOND_INVAR_STRAIN_
%token <string> PARTICLES_
%token <string> DISSIPATION_
%token <string> SPECIES_MASS_
%token <string> HEAT_ENERGY_
%token <string> MOMENTUMX_
%token <string> MOMENTUMY_
%token <string> MOMENTUMZ_
%token <string> STRESS_TRACE_
%token <string> POSITIVE_FILL_
%token <string> NEGATIVE_FILL_
%token <string> PLASTICITY_
%token <string> EVP_HYPER_
%token <string> SHEAR_HARDEN_
%token <string> DENSE_POWER_LAW_
%token <string> NO_MODEL_
%token <string> CARREAU_WLF_ 
%token <string> CARREAU_SUSPENSION_ 
%token <string> POWERLAW_SUSPENSION_
%token <string> SAT_DES_
%token <string> COMP_LEVEL_
%token <string> CAP_PRES_
%token <string> NON_LINEAR_
%token <string> WLF_
%token <string> CONSTANT2_
%token <string> NS_
%token <string> EB_
%token <string> BIGX1_ 
%token <string> X1_ 
%token <string> BIGX2_ 
%token <string> X2_ 
%token <string> BIGX3_ 
%token <string> X3_ 
%token <string> CURVATURE_ 
%token <string> H_ 
%token <string> ISOSURFACE_
%token <string> GOMA_MAT_
%token <string> CHEMKIN_MAT_
%token <string> STEFAN_MAXWELL_CHARGED_
%token <string> FICKIAN_CHARGED_
%token <string> ELECTRODE_KINETICS_
%token <string> EPOXY_DEA_
%token <string> JOULE_
%token <string> VISC_DISS_
%token <string> BLANK_LINE_
%token <string> SOLIDIFICATION_
%token <string> KOZENY_CARMAN_
%token <string> TENSOR_
%token <string> FLOWINGLIQUID_
%token <string> MOLTEN_GLASS_
%token <string> DARCY_FICKIAN_
%token <string> GENERALIZED_FICKIAN_
%token <string> STEFAN_MAXWELL_
%token <string> STEFAN_MAXWELL_VOLUME_
%token <string> SUSPENSION_BALANCE_
%token <string> HYDRO_ 
%token <string> SUSP_BAL_ 
%token <string> GENERALIZED_FREE_VOL_ 
%token <string> GENERALIZED_ 
%token <string> KELVIN_ 
%token <string> FLAT_ 
%token <string> ANTOINE_ 
%token <string> RIEDEL_ 
%token <string> PROJECTED_CARTESIAN_
%token <string> REALLY_MEANT_
%token <string> EXAMPLE_TOKEN2_		/*EXAMPLE CARD TOKEN 2 DELCARATION*/
%token <string> EXAMPLE_TOKEN1_		/*EXAMPLE CARD TOKEN 1 DECLARATION*/

%type <string> navier_stokes_source_card
%type <integer> navier_stokes_source_model
%type <integer> heat_source_model
%type <integer> species_source_model
%type <string> material_variable_initialize_card
%type <string> species_source_card
%type <string> solid_body_source_card
%type <string> heat_source_card
%type <string> solution_temperature_card
%type <integer> solution_temperature_model_name
%type <string> non_condensable_molecular_weight_card
%type <string> non_volatile_molar_volume_card
%type <string> non_volatile_specific_volume_card
%type <string> mass_source_card
%type <string> current_source_card
%type <integer> current_source_model
%type <string> combined
%type <string> input_file
%type <string> input_file_card
%type <string> fem_file_specifications_delimiter_card
%type <string> fem_file_card		 
%type <string> output_file_card		  
%type <string> guess_file_card	 	  
%type <string> soln_file_card		 
%type <string> write_int_results_card	 	
%type <string> problem_desc_delimiter_card
%type <string> number_of_materials_card
%type <string> mat_file_name_card
%type <string> coord_system_card
%type <integer> coord_system_type
%type <string> element_mapping_card
%type <integer> element_mapping_type	
%type <string> mesh_motion_card
%type <integer> mesh_motion_type
%type <string> bulk_species_card
%type <string> really_meant_card
%type <string> number_of_equations_card
%type <string> mesh1_equation_card
%type <string> mesh2_equation_card
%type <string> mesh3_equation_card
%type <integer> mesh_equation_weight
%type <string> mom_solid1_equation_card
%type <string> mom_solid2_equation_card
%type <string> mom_solid3_equation_card
%type <integer> mom_solid_equation_weight
%type <string> momentum1_equation_card
%type <string> momentum2_equation_card
%type <string> momentum3_equation_card%type <string> pmomentum1_equation_card
%type <string> pmomentum2_equation_card
%type <string> pmomentum3_equation_card
%type <integer> pmomentum_equation_weight
%type <string> stress11_equation_card
%type <string> stress12_equation_card
%type <string> stress13_equation_card
%type <string> stress22_equation_card
%type <string> stress23_equation_card
%type <integer> stress_equation_weight
%type <string> gradient11_equation_card
%type <string> gradient12_equation_card
%type <string> gradient13_equation_card
%type <string> gradient21_equation_card
%type <string> gradient22_equation_card
%type <string> gradient23_equation_card
%type <string> gradient31_equation_card
%type <string> gradient32_equation_card
%type <string> gradient33_equation_card
%type <integer> gradient_equation_weight
%type <integer> momentum_equation_weight
%type <string> voltage_equation_card
%type <integer> voltage_equation_weight
%type <string> continuity_equation_card
%type <integer> continuity_equation_weight
%type <string> energy_equation_card
%type <integer> energy_equation_weight
%type <string> species_bulk_equation_card
%type <integer> species_bulk_equation_weight
%type <string> species_surf_equation_card
%type <integer> species_surf_equation_weight
%type <string> level_set_equation_card
%type <integer> level_set_equation_weight
%type <string> fill_equation_card
%type <string> level_set_implicit_embedded_surface_card
%type <integer> fill_equation_weight
%type <string> mat_file
%type <string> mat_card	 		
%type <string> density_card
%type <string> default_database_card
%type <integer> database_type
%type <string> solid_constitutive_equation_card	
%type <string> convective_lagrangian_velocity_card	
%type <string> lame_mu_card
%type <integer> lame_mu_model_name
%type <string> lame_lambda_card
%type <integer> lame_lambda_model_name	
%type <string> conductivity_card
%type <integer> conductivity_model
%type <string> heat_capacity_card
%type <integer> heat_capacity_model
%type <string> stress_free_solvent_vol_frac_card
%type <string> liquid_constitutive_card
%type <integer> liquid_consitiutive_model_name	
%type <string> mechanical_properties_viscosity_card
%type <integer> viscosity_model_name	
%type <string> low_rate_viscosity_card
%type <string> low_rate_viscosity_model_name	
%type <string> power_law_exponent_equation_card
%type <string> power_law_exponent_model_name	
%type <string> high_rate_viscosity_card
%type <string> high_rate_viscosity_model_name	
%type <string> time_constant_card
%type <string> time_constant_model_name	
%type <string> aexp_card
%type <string> aexp_model_name	
%type <string> thermal_exponent_card
%type <string> thermal_exponent_model_name	
%type <string> yield_stress_card
%type <integer> yield_stress_model_name	
%type <string> yield_exponent_card
%type <integer> yield_exponent_model_name	
%type <string> suspension_maximum_packing_card
%type <string> suspension_maximum_packing_model_name	
%type <string> suspension_species_number_card
%type <string> cure_gel_point_card
%type <string> cure_gel_point_model_name	
%type <string> cure_a_exponent_card
%type <string> cure_a_exponent_model_name	
%type <string> cure_b_exponent_card
%type <string> cure_b_exponent_model_name	
%type <string> cure_species_number_card
%type <string> cure_species_number_card
%type <string> polymer_constitutive_equation_card
%type <string> polymer_constitutive_equation_model_name	
%type <string> polymer_stress_formulation_card
%type <string> polymer_stress_formulation_model_name	
%type <string> polymer_weight_function_card
%type <string> polymer_weight_function_model_name	
%type <string> polymer_weighting_card
%type <string> polymer_viscosity_card
%type <string> polymer_viscosity_model_name	
%type <string> polymer_time_constant_card
%type <string> polymer_time_constant_model_name	
%type <string> mobility_parameter_card
%type <string> mobility_parameter_model_name	
%type <string> surface_tension_card
%type <string> surface_tension_model_name
%type <string> shear_rate_equation_card
%type <integer> shear_rate_equation_weight
%type <string> species_unknown_equation_card
%type <integer> species_unknown_equation_weight
%type <string> general_specifications_delimiter_card
%type <string> debug_card
%type <string> initial_guess_card
%type <string> initialize_card
%type <string> anneal_mesh_on_output_card
%type <string> external_field_card
%type <string> number_of_processors_card
%type <string> output_level_card
%type <string> time_integration_specifications_delimiter_card
%type <string> time_integration_card
%type <string> delta_t_card
%type <string> maximum_number_of_time_steps_card
%type <string> maximum_time_card
%type <string> maximum_time_step_card
%type <string> minimum_time_step_card
%type <string> time_step_parameter_card
%type <string> printing_frequency_card
%type <string> second_frequency_time_card
%type <string> initial_time_card
%type <string> fill_subcycle_card
%type <string> solver_specifications_delimiter_card
%type <string> preconditioner_card
%type <string> matrix_scaling_card
%type <string> matrix_subdomain_solver_card
%type <string> matrix_residual_norm_type_card
%type <string> matrix_output_type_card
%type <string> matrix_factorization_reuse_card
%type <string> matrix_factorization_overlap_card
%type <string> matrix_polynomial_order_card
%type <string> size_of_krylov_subspace_card
%type <string> orthogonalization_card
%type <string> maximum_linear_solve_iterations_card
%type <string> matrix_auxiliary_vector_card
%type <string> matrix_drop_tolerance_card
%type <string> number_of_newtonian_iterations_card
%type <string> newton_correction_factor_card
%type <string> normalized_residual_tolerance_card
%type <string> modified_newton_tolerance_card
%type <string> jacobian_reform_time_stride_card
%type <string> residual_ratio_tolerance_card
%type <string> pressure_stabilization_card
%type <string> pressure_stabilization_scaling_card
%type <string> solver_pressure_datum_card
%type <string> umf_idim_card
%type <string> umf_xdim_card
%type <string> disable_viscosity_sensitivities_card
%type <string> post_processing_specification_section_delimiter_card
%type <string> stream_function_card
%type <string> streamwise_normal_stress_card
%type <string> mean_shear_rate_card
%type <string> pressure_contours_card
%type <string> first_invariant_of_strain_card
%type <string> second_invariant_of_strain_card
%type <string> third_invariant_of_strain_card
%type <string> mesh_dilatation_card
%type <string> velocity_divergence_card
%type <string> navier_stokes_residuals_card
%type <string> moving_mesh_residuals_card
%type <string> mass_diffusion_vectors_card
%type <string> mass_fluxlines_card
%type <string> energy_conduction_vectors_card
%type <string> energy_fluxlines_card
%type <string> time_derivatives_card
%type <string> mesh_stress_tensor_card
%type <string> mesh_strain_tensor_card
%type <integer> density_model
%type <string> real_solid_stress_tensor_card
%type <string> porous_saturation_card
%type <string> lagrangian_convection_card
%type <string> normal_and_tangent_vectors_card
%type <string> post_processing_viscosity_card
%type <string> user_defined_post_processing_card
%type <string> flux_card
%type <string> end_of_flux_card
%type <string> end_of_flux_sens_card
%type <integer> yes_no_option
%type <floating> floating_point_constant_list
%type <floating> optional_floating_point_constant_list
%type <string> post_processing_fluxes_section_delimiter_card
%type <string> rot_surface_card
%type <string> rot_edge_card
%type <string> rot_eqn_type
%type <string> rot_vertex_card
%type <string> rot_volume_card
%type <string> end_of_rot_card
%type <string> rot_eq_type
%type <string> bc_name_or_rotation_string
%type <string> rot_seed_method
%type <floating> optional_float
%type <integer> optional_integer
%type <integer> integer
%type <string> bc_format_1 
%type <string> bc_format_2 	
%type <string> bc_format_4 
%type <string> bc_format_6 
%type <string> bc_format_7 
%type <string> bc_format_9 
%type <string> bc_format_11 	   
%type <string> bc_format_16 	
%type <string> bc_format_17 	   
%type <string> bc_format_18 	   
%type <string> bc_format_12 	   
%type <string> bc_format_13 	   
%type <string> bc_format_14 	   
%type <string> bc_format_19 	   
%type <string> bc_format_20 	   
%type <string> bc_format_15 	   
%type <string> bc_format_22
%type <string> bc_format_23 	   
%type <string> bc_format_25 	   
%type <string> bc_format_26 	   
%type <string> bc_format_27 	   
%type <string> bc_format_28 	   
%type <string> bc_format_29 	   
%type <string> bc_format_30 	   
%type <string> bc_format_32      
%type <string> table_bcs 
%type <string> bc_card
%type <string> bc_type
%type <integer> bc_id
%type <integer> bc_keyword_1
%type <integer> bc_keyword_2
%type <integer> bc_keyword_3
%type <string> bc_format_3
%type <string> bc_format_5
%type <string> bc_format_8
%type <string> bc_format_10
%type <string> bc_format_21
%type <string> bc_format_24
%type <string> end_table_card
%type <string> file_equals_with_name_equals_card
%type <string> file_equals_card
%type <floating> table_data_card
%type <string> optional_file_equals_string
%type <string> optional_name_equals_string
%type <integer> solution_algorithm
%type <string> valid_preconditioner
%type <string> valid_subdomain_solver
%type <string> matrix_storage_format_card
%type <string> valid_storage_format
%type <string> valid_matrix_scaling
%type <string> valid_matrix_output_type
%type <string> augmenting_conditions_specifications_section_delimiter_card
%type <string> augmenting_conditions_initial_guess_card
%type <string> number_of_augmenting_conditions_card
%type <string> augmenting_condition_card
%type <string> end_of_augmenting_conditions_card
%type <string> boundary_condition_specifications_delimiter_card
%type <string> number_of_boundary_conditions_card
%type <string> end_of_BC_card
%type <string> gd_time_table_option
%type <string> optional_LINEAR
%type <string> equation_name
%type <string> var_name
%type <string> post_var_name
%type <integer> interpolation_method
%type <integer> abcissa
%type <integer> ordinate
%type <string> end_of_equations_card
%type <string> end_of_materials_card
%type <integer> discontinuous_velo_bc_keyword
%type <integer> latent_heat_internal_bc_keyword
%type <string> file_card_option
%type <string> domain_mapping_file_card
%type <string> write_initial_solution_card
%type <string> number_of_jacobian_file_dumps_card
%type <string> optional_string
%type <integer> initial_guess_option
%type <integer> variable_name_unknown
%type <integer> external_field_interpolation_method
%type <string> pressure_datum_card
%type <integer> optional_pressure_datum_units
%type <integer> time_integraton_option
%type <string> level_set_interface_tracking_card
%type <string> level_set_length_scale_card
%type <string> level_set_control_width_card
%type <string> level_set_renormalization_method_card
%type <string> level_set_renormalization_tolerance_card
%type <string> level_set_renormalization_frequency_card
%type <string> level_set_initialization_method_card
%type <integer> level_set_remornalization_method
%type <string> level_set_initialization_surface_card
%type <string> hunting_conditions_specificatins_delimiter_card
%type <string> number_of_hunting_conditions_card
%type <string> hc_card
%type <string> end_of_hc_card 
%type <integer> hc_type
%type <integer> data_sens_type
%type <string> continuation_specificatins_delimiter_card
%type <string> continuation_specificatins_delimiter_card 
%type <string> continuation_card 
%type <string> continuation_type_card 
%type <string> boundary_condition_id_card
%type <string> material_id_card 
%type <string> material_property_tag_card 
%type <string> material_property_tag_subindex_card 
%type <string> initial_parameter_value_card 
%type <string> final_parameter_value_card 
%type <string> delta_s_card 
%type <string> maximum_number_of_path_steps_card 
%type <string> maximum_path_value_card 
%type <string> minimum_path_step_card  
%type <string> maximum_path_step_card 
%type <string> path_step_parameter_card 
%type <string> path_step_error_card 
%type <string> continuation_printing_frequency_card 
%type <string> second_frequency_card 
%type <string> loca_method_card 
%type <string> continuation_order_card 
%type <string> step_control_aggressiveness_card 
%type <string> alc_desired_solution_fraction_card 
%type <string> alc_max_parameter_sensitivity_card 
%type <string> alc_tangent_factor_exponent_card 
%type <string> alc_tangent_factor_step_limit_card 
%type <string> number_of_continuation_conditions_card 
%type <string> cc_card 
%type <string> end_of_cc_card 
%type <string> tp_continuation_type_card 
%type <string> tp_boundary_condition_id_card 
%type <string> tp_bc_data_float_tag_card 
%type <string> tp_parameter_material_id_card 
%type <string> tp_parameter_material_property_tag_card 
%type <string> tp_material_property_tag_subindex_card 
%type <string> initial_guess_of_tp_parameter_card 
%type <string> tp_parameter_final_value_card 
%type <string> number_of_tp_continuation_conditions_card 
%type <string> tc_card 
%type <string> end_of_tc_card 
%type <string> boundary_condition_data_float_tag_card
%type <integer> continuation_type
%type <integer> continuation_method
%type <integer> loca_method
%type <integer> read_option
%type <integer> ac_userbc_augmenting_condition
%type <integer> ac_usermat_augmenging_condition
%type <integer> ac_volume_augmenting_condition
%type <integer> ac_flux_augmenting_condition
%type <integer> mfid
%type <integer> cc_type
%type <string> valid_residual_norm_type
%type <string> valid_matrix_reuse_type
%type <string> valid_matrix_factorization_overlap
%type <string> valid_matrix_overlap_type
%type <string> matrix_overlap_type_card
%type <string> valid_matrix_auxiliary_vector_type
%type <string> valid_matrix_reorder_type
%type <string> matrix_reorder_card
%type <string> matrix_factorization_save_card
%type <string> matrix_rilu_relax_factor_card
%type <string> matrix_relative_threshold_card
%type <string> matrix_absolute_threshold_card
%type <string> orthogonalization_option
%type <string> normalized_corection_tolerance_card
%type <string> linear_stability_card
%type <integer> valid_linear_stability_option
%type <string> filter_concentration_card
%type <string> eigensolver_specifications_delimiter_card 
%type <string> post_processing_data_var 
%type <string> eigen_number_of_modes_card  
%type <string> eigen_record_modes_card  
%type <string> eigen_size_of_krylov_subspace_card  
%type <string> eigen_maximum_iterations_card  
%type <string> eigen_number_of_filter_steps_card  
%type <string> eigen_recycle_card  
%type <string> eigen_tolerance_card  
%type <string> eigen_initial_vector_weight_card  
%type <string> eigen_initial_shifts_card  
%type <string> eigen_wave_numbers_card  
%type <string> cross_stream_shear_rate_card 
%type <string> fill_contours_card 
%type <string> concentration_contours_card 
%type <string> stress_contours_card 
%type <string> particle_velocity_divergence_card 
%type <string> total_velocity_divergence_card 
%type <string> post_processing_viscosity_card 
%type <string> diffusive_mass_flux_vectors_card 
%type <string> real_solid_stress_tensor_card 
%type <string> viscoplastic_def_grad_tensor_card 
%type <string> error_zz_velocity_card 
%type <string> error_zz_heat_flux_card 
%type <string> error_zz_pressure_card 
%type <string> total_density_of_solvents_in_porous_media_card 
%type <string> density_of_solvents_in_gas_phase_in_porous_media_card 
%type <string> density_of_liquid_phase_in_porous_media_card 
%type <string> gas_phase_darcy_velocity_in_porous_media_card 
%type <string> liquid_phase_darcy_velocity_in_porous_media_card 
%type <string> grid_peclet_number_in_porous_media_card 
%type <string> SUPG_velocity_in_porous_media_card 
%type <string> vorticity_vector_card 
%type <string> capillary_pressure_in_porous_media_card
%type <string> error_zz_element_size_card
%type <string> zz_error_type
%type <string> valid_pp_flux_name
%type <integer> optional_profile
%type <integer> valid_flux_sensitivity_type
%type <string> post_processing_fluxes_sensitivities_section_delimiter_card 
%type <string> flux_sensitivity_card 
%type <string> end_of_flux_sens_card 
%type <string> post_processing_data_specifications_delimiter_card 
%type <string> post_processing_data_card 
%type <string> end_of_data_card 
%type <string> post_processing_data_sensitivity_delimiter_card 
%type <string> data_sens_card 
%type <string> end_of_data_sens_card 
%type <string> post_processing_particle_tracking_calculations_delimiter_card 
%type <string> particle_card 
%type <string> end_of_particles_card 
%type <string> post_processing_volumetric_integration_delimiter_card 
%type <string> volume_int_card 
%type <string> end_of_volume_int_card 
%type <string> pp_vol_names
%type <floating> float
%type <integer> solid_constitutive_equation_model
%type <string> plasticity_equation_card
%type <integer> plasticity_equation_model
%type <integer> convective_lagrangian_velocity_model
%type <string> solid_thermal_expansion_card
%type <integer> solid_thermal_expansion_model
%type <string> mp_table_independent_variable_name
%type <string> optional_mp_table_independent_variable_name
%type <integer> mp_table_interpolation_method
%type <string> reference_temperature_card
%type <string> diffusivity_card
%type <integer> diffusivity_model
%type <string> volume_expansion_card
%type <string> electrical_conductivity_card
%type <integer> electrical_conductivity_model
%type <string> media_type_card
%type <string> media_type
%type <string> porosity_card
%type <integer> porosity_model
%type <string> permeability_card
%type <string> flowingliquid_viscosity_card
%type <string> inertia_coefficient_card
%type <integer> flowingliquid_viscosity_model
%type <integer> permeability_model
%type <string> porous_diffusion_constitutive_equation_card
%type <integer> porous_diffusion_constitutive_equation_model
%type <string> diffusion_constitutive_equation_card
%type <integer> diffusion_constitutive_equation_model
%type <string> liquidus_temperature_card
%type <string> solidus_temperature_card
%type <string> plastic_viscosity_card
%type <string> latent_heat_vaporization_card
%type <string> latent_heat_fusion_card
%type <string> species_volume_expansion_card
%type <string> reference_concentration_card
%type <stripg> vapor_pressure_card
%type <string> porous_latent_heat_vaporization_card
%type <integer> plastic_viscosity_model
%type <string> EVP_yeild_stress_card
%type <integer> EVP_yeild_stress_model
%type <string> pseudo_solid_constitutive_equation_card
%type <integer> pseudo_solid_constitutive_equation_model
%type <string> pseudo_solid_lame_mu_card
%type <integer> pseudo_solid_lame_mu_model
%type <string> pseudo_solid_lame_lambda_card
%type <integer> pseudo_solid_lame_lambda_model
%type <string> thermal_wlf_constant2_card
%type <integer> thermal_wlf_constant2_model
%type <string> flory_huggins_card
%type <string> flory_huggins_parameter_card
%type <integer> vapor_pressure_model
%%

combined: PARSER_IN_INPUT_FILE_MODE input_file {} | PARSER_IN_MAT_FILE_MODE mat_file {};

/************************************* FEM File Specifications rules: **************************************/
input_file:	input_file_card 
		  {  floating_point_constant_list_index = -1; 
                     line_number++;
                     if (ProcID == 0)  fprintf(parser_log,"\n%d ",line_number ); 
                  } 
		| input_file input_file_card 
		  {  if ( error_found_in_last_line == 0)
		     {
		       floating_point_constant_list_index = -1; 
                       line_number++;
                       if (ProcID == 0)  fprintf(parser_log,"\n%d ",line_number );
                     }
                  }
		| error {}
		;

input_file_card:	/* empty */ {}
                | fem_file_specifications_delimiter_card
                | fem_file_card 
		| output_file_card 
		| guess_file_card 
		| domain_mapping_file_card
		| soln_file_card
		| write_int_results_card  
		| general_specifications_delimiter_card 
		| debug_card 
		| number_of_jacobian_file_dumps_card
		| initial_guess_card 
		| initialize_card
		| external_field_card 
		| pressure_datum_card 
		| anneal_mesh_on_output_card
		| number_of_processors_card
		| output_level_card
		| time_integration_specifications_delimiter_card {}
		| time_integration_card {}
		| delta_t_card {}
		| maximum_number_of_time_steps_card {}
		| maximum_time_card {}
		| maximum_time_step_card {}
		| minimum_time_step_card {}
		| time_step_parameter_card {}
		| printing_frequency_card {}
		| second_frequency_time_card {}
		| initial_time_card {}
		| fill_subcycle_card {}
		| time_step_error_card {}
		| level_set_interface_tracking_card {}
		| level_set_length_scale_card {}
		| level_set_control_width_card {}
		| level_set_renormalization_method_card {}
		| level_set_renormalization_tolerance_card {}
		| level_set_renormalization_frequency_card {}
		| level_set_initialization_method_card {}
		| level_set_initialization_surface_card {}	
		| continuation_specificatins_delimiter_card {}
		| continuation_card {}
		| continuation_type_card {}
		| boundary_condition_id_card {}
		| boundary_condition_data_float_tag_card {}
		| material_id_card {}
		| material_property_tag_card {}
		| material_property_tag_subindex_card {}
		| initial_parameter_value_card {}
		| final_parameter_value_card {}
		| delta_s_card {}
		| maximum_number_of_path_steps_card {}
		| maximum_path_value_card {}
		| minimum_path_step_card {}
		| maximum_path_step_card {}
		| path_step_parameter_card {}
		| path_step_error_card {}
		| continuation_printing_frequency_card {}
		| second_frequency_card {}
		| loca_method_card {}
		| continuation_order_card {}
		| step_control_aggressiveness_card {}
		| alc_desired_solution_fraction_card {}
		| alc_max_parameter_sensitivity_card {}
		| alc_tangent_factor_exponent_card {}
		| alc_tangent_factor_step_limit_card {}
		| number_of_continuation_conditions_card {}
		| cc_card {}
		| end_of_cc_card {}
		| tp_continuation_type_card {}
		| tp_boundary_condition_id_card {}
		| tp_bc_data_float_tag_card {}
		| tp_parameter_material_id_card {}
		| tp_parameter_material_property_tag_card {}
		| tp_material_property_tag_subindex_card {}
		| initial_guess_of_tp_parameter_card {}
		| tp_parameter_final_value_card {}
		| number_of_tp_continuation_conditions_card {}
		| tc_card {}
		| end_of_tc_card {}
		| hunting_conditions_specificatins_delimiter_card {}
		| number_of_hunting_conditions_card {}
		| hc_card {}		
		| end_of_hc_card {}
		| solver_specifications_delimiter_card {}
		| matrix_storage_format_card {}
		| solution_algorithm_card {}
		| preconditioner_card {}
		| matrix_fill_factor_card {}
		| matrix_scaling_card {}
		| matrix_subdomain_solver_card {}
		| matrix_residual_norm_type_card {}
		| matrix_output_type_card {}
		| matrix_factorization_reuse_card {}
		| matrix_factorization_overlap_card {}
		| matrix_polynomial_order_card {}
		| size_of_krylov_subspace_card {}
		| orthogonalization_card {}
		| maximum_linear_solve_iterations_card {}
		| matrix_auxiliary_vector_card {}
		| matrix_drop_tolerance_card {}
		| number_of_newtonian_iterations_card {}
		| newton_correction_factor_card {}
		| normalized_residual_tolerance_card {}
		| modified_newton_tolerance_card {}
		| jacobian_reform_time_stride_card {}
		| residual_ratio_tolerance_card {}
		| pressure_stabilization_card {}
		| pressure_stabilization_scaling_card {}
		| solver_pressure_datum_card {}
		| umf_idim_card {}
		| umf_xdim_card {}
		| disable_viscosity_sensitivities_card {}
		| matrix_graph_fillin_card {}
		| matrix_overlap_type_card {}
		| matrix_reorder_card {}
		| matrix_factorization_save_card {}
		| matrix_ilut_fill_factor_card {}
		| matrix_rilu_relax_factor_card {}
		| matrix_bilu_threshold_card {}
		| matrix_relative_threshold_card {}
		| matrix_absolute_threshold_card {}
		| size_of_krylov_subspace_card {}
		| linear_stability_card {}
		| filter_concentration_card {}
		| normalized_corection_tolerance_card {}
		| eigensolver_specifications_delimiter_card {}
		| eigen_number_of_modes_card {}
		| eigen_record_modes_card {}
		| eigen_size_of_krylov_subspace_card {}
		| eigen_maximum_iterations_card {}
		| eigen_number_of_filter_steps_card {}
		| eigen_recycle_card {}
		| eigen_tolerance_card {}
		| eigen_initial_vector_weight_card {}
		| eigen_initial_shifts_card {}
		| eigen_wave_numbers_card {}	
		| augmenting_conditions_specifications_section_delimiter_card {}
		| augmenting_conditions_initial_guess_card {}
		| number_of_augmenting_conditions_card {}
		| augmenting_condition_card {}
		| end_of_augmenting_conditions_card {}	
		| bc_card {}
		| boundary_condition_specifications_delimiter_card {}
		| number_of_boundary_conditions_card {}
		| end_of_BC_card {}
		| rot_surface_card {}
		| rot_edge_card	{}	
		| rot_vertex_card {}
		| end_of_rot_card {}	
		| problem_desc_delimiter_card {}
		| number_of_materials_card {}
		| mat_file_name_card {}
		| coord_system_card {}
		| element_mapping_card {} 
		| mesh_motion_card {}
		| bulk_species_card {}
		| really_meant_card {}		
		| end_of_materials_card {}
		| number_of_equations_card {}
		| mesh1_equation_card {}
		| mesh2_equation_card {}
		| mesh3_equation_card {}
		| mom_solid1_equation_card {}
		| mom_solid2_equation_card {}
		| mom_solid3_equation_card {}
		| momentum1_equation_card {}
		| momentum2_equation_card {}
		| momentum3_equation_card {}
		| pmomentum1_equation_card {}
		| pmomentum2_equation_card {}
		| pmomentum3_equation_card {}
		| voltage_equation_card {}
		| continuity_equation_card {}
		| energy_equation_card {}
		| species_bulk_equation_card {}
		| species_surf_equation_card {}
		| stress11_equation_card {}
		| stress12_equation_card {}
		| stress13_equation_card {}
		| stress22_equation_card {}
		| stress23_equation_card {}
		| gradient11_equation_card {}
		| gradient12_equation_card {}
		| gradient13_equation_card {}
		| gradient21_equation_card {}
		| gradient22_equation_card {}
		| gradient23_equation_card {}
		| gradient31_equation_card {}
		| gradient32_equation_card {} 
		| gradient33_equation_card {}
		| level_set_equation_card {}
		| fill_equation_card {}
		| level_set_implicit_embedded_surface_card {}
		| shear_rate_equation_card {}
		| species_unknown_equation_card {}
		| end_of_equations_card	 {}
		| post_processing_specification_section_delimiter_card {}
		| stream_function_card {}
		| streamwise_normal_stress_card {}
		| cross_stream_shear_rate_card  {}
		| mean_shear_rate_card {}
		| pressure_contours_card {}
		| fill_contours_card  {}
		| concentration_contours_card  {}
		| stress_contours_card  {}
		| first_invariant_of_strain_card {}
		| second_invariant_of_strain_card {}
		| third_invariant_of_strain_card {}
		| mesh_dilatation_card {}
		| velocity_divergence_card {}
		| particle_velocity_divergence_card  {}
		| total_velocity_divergence_card  {}
		| post_processing_viscosity_card  {}
		| navier_stokes_residuals_card {}
		| moving_mesh_residuals_card {}
		| mass_diffusion_vectors_card {}
		| diffusive_mass_flux_vectors_card  {}
		| mass_fluxlines_card {}
		| energy_conduction_vectors_card {}
		| energy_fluxlines_card {}
		| time_derivatives_card {}
		| mesh_stress_tensor_card {}
		| real_solid_stress_tensor_card  {}
		| mesh_strain_tensor_card {}
		| viscoplastic_def_grad_tensor_card  {}
		| lagrangian_convection_card {}
		| normal_and_tangent_vectors_card {}
		| error_zz_velocity_card  {}
		| error_zz_heat_flux_card  {}
		| error_zz_pressure_card  {}
		| user_defined_post_processing_card {}
		| porous_saturation_card {}
		| total_density_of_solvents_in_porous_media_card  {}
		| density_of_solvents_in_gas_phase_in_porous_media_card  {}
		| density_of_liquid_phase_in_porous_media_card  {}
		| gas_phase_darcy_velocity_in_porous_media_card  {}
		| liquid_phase_darcy_velocity_in_porous_media_card  {}
		| capillary_pressure_in_porous_media_card {}
		| grid_peclet_number_in_porous_media_card  {}
		| SUPG_velocity_in_porous_media_card  {}
		| vorticity_vector_card  {}
		| error_zz_element_size_card {}
		| post_processing_fluxes_section_delimiter_card {}
		| flux_card {}
		| end_of_flux_card {}
		| post_processing_fluxes_sensitivities_section_delimiter_card {}
		| flux_sensitivity_card {}
		| end_of_flux_sens_card {}
		| post_processing_data_specifications_delimiter_card {}
		| post_processing_data_card {}
		| end_of_data_card {}
		| post_processing_data_sensitivity_delimiter_card {}
		| data_sens_card {}
		| end_of_data_sens_card {}
		| post_processing_particle_tracking_calculations_delimiter_card {}
		| particle_card {}
		| end_of_particles_card {}
		| post_processing_volumetric_integration_delimiter_card {}
		| volume_int_card {}
		| end_of_volume_int_card {}	
		| table_data_card {}
		| end_table_card {}
		| error {}
		;

fem_file_specifications_delimiter_card:
		/* empty */ {} 
		| FEM_ FILE_ SPECIFICATIONS_ CR_
		{ 
		  card_read("FEM_FILE_SPECIFICATIONS_CARD",file_index);
		} 
		| error {}
		;

fem_file_card:	/*empty*/ {}
		| FEM_ FILE_ EQUALS_ STRING_ CR_
		{
		  strcpy(ExoFile,$4); 
		  card_read("FEM_FILE_CARD",file_index);
		  
		}
		| error {}		
		;
		
output_file_card: /* empty */ {} 
		| OUTPUT_ EXODUS_ II_ FILE_ EQUALS_ STRING_ CR_
		{
		  strcpy(ExoFileOut,$6);
		  card_read("OUTPUT_FILE_CARD",file_index);
		  
		}
		| error {}
		;

guess_file_card:
		/* empty */ {} 
		| GUESS_ FILE_ EQUALS_ file_card_option CR_
		{
		  strcpy(Init_GuessFile,$4);
		  card_read("GUESS_FILE_CARD",file_index);
		}
		| error {}
		;


soln_file_card:	
		/* empty */ {} 
		| SOLN_ FILE_ EQUALS_ file_card_option CR_
		{
		  strcpy(Soln_OutFile,$4);
		  card_read("SOLN_FILE_CARD",file_index);
		}
		| error {}
		;

file_card_option:
		  NONE_		{strcpy($$,'\0');}
		| NO_ 		{strcpy($$,'\0');}
		| STRING_ 	{}
		| error {}
		;
		
domain_mapping_file_card:	
		/* empty */ {} 
		| DOMAIN_ MAPPING_ FILE_ EQUALS_ STRING_ CR_
		{
		  strcpy(DomainMappingFile,$5);
		  card_read("DOMAIN_MAPPING_FILE_CARD",file_index);
		}
		| error {}
		;		
		
write_int_results_card:	/* empty */ {} 
		| WRITE_ INTERMEDIATE_ RESULTS_ EQUALS_ YES_  CR_
		{  
		  Write_Intermediate_Solutions = TRUE;
		  card_read("WRITE_INTERMEDIATE_RESULTS_CARD",file_index);
		  
		}
		| WRITE_ INTERMEDIATE_ RESULTS_ EQUALS_ NO_  CR_
		{  
		  Write_Intermediate_Solutions = FALSE;
		  card_read("WRITE_INTERMEDIATE_RESULTS_CARD",file_index);
		}
		| error {}
		;

write_initial_solution_card:	/* empty */ {} 
		| WRITE_ INITIAL_ SOLUTION_ EQUALS_ YES_  CR_
		{  
		  Write_Initial_Solution = TRUE;
		  card_read("WRITE_INITIAL_SOLUTION_CARD",file_index);
		  
		}
		| WRITE_ INITIAL_ SOLUTION_ EQUALS_ NO_  CR_
		{  
		  Write_Initial_Solution = FALSE;
		  card_read("WRITE_INITIAL_SOLUTION_CARD",file_index);
		}
		| error {}
		;


general_specifications_delimiter_card:
		/* empty */ {}
		| GENERAL_ SPECIFICATIONS_  CR_
		{
		  card_read("GENERAL_SPECIFICATIONS_DELIMITER_CARD",file_index);  
		}
		| error {}
		;

number_of_processors_card:
		/* empty */ {}
		| NUMBER_ OF_ PROCESSORS_ EQUALS_ integer  CR_
		{
		  sprintf(msg, "The Number of processors card is no longer being used.  Card being ignored."); 
		  declare_warning(msg);		   
		}
		| error {}
		;

output_level_card:
		/* empty */ {}
		| OUTPUT_ LEVEL_ EQUALS_ integer  CR_
		{
		  card_read("OUTPUT_LEVEL_CARD",file_index);
		  Iout = $4;  
		}
		| error {}
		;		
		
debug_card:
		/* empty */ {}
		| DEBUG_ EQUALS_ integer  CR_
		{
		  card_read("DEBUG_CARD",file_index);
		  Debug_Flag = $3;
		    
		}
		| error {}
		;


number_of_jacobian_file_dumps_card:
		/* empty */ {}
		| NUMBER_ OF_ JACOBIAN_ FILE_ DUMPS_ EQUALS_ integer  CR_
		{
		  card_read("NUMBER_OF_JACOBIAN_FILE DUMPS_CARD",file_index);  
		  #ifdef MATRIX_DUMP
		  Number_Jac_Dump = $7;
		  #endif
		}
		| error {}
		;
				
initial_guess_card:
		/* empty */ {}
		| INITIAL_ GUESS_ EQUALS_ initial_guess_option optional_string optional_string CR_
		{
		  switch($4)
      		  { 
                    case 0: 
                    { Guess_Flag = 0; break;}
                    case 1: 
                    { Guess_Flag = 1; break;}
                    case 2: 
                    { Guess_Flag = 2; break;}
                    case 3: 
                    { Guess_Flag = 3; break;}                                                            
                    case 4: 
                    { Guess_Flag = 4; break;}
                    case 5: 
                    { Guess_Flag = 5; break;}
                    case 6: 
                    { if (strcmp($5,"NA"))
                      {
                      Guess_Flag = 6;
                      strcpy(ExoAuxFile, $5);
                      } else {
		      sprintf(msg, " Undecipherable 2 options for Initial guess."); 
		      declare_error(msg);
		      }                    
                    } /*end case 6 */  
		  }  /* end switch */
		  card_read("INITIAL_GUESS_CARD",file_index);
		}
		| error {}
		;

initial_guess_option:
		  ZERO_		{$$=0;}
		| RANDOM_ 	{$$=1;}
		| ONE_		{$$=2;}
		| SINES_	{$$=3;}
		| READ_		{$$=4;}
		| READ_EXOII_	 {$$=5;}
		| READ_EXOII_FILE_ {$$=6;}
		| error {}
		;
		
optional_string:
		/* empty */ {strcpy($$,"NA");}
		| STRING_ {}
		| error {}
		;		
		
initialize_card:
		/* empty */ {}
		| INITIALIZE_ variable_name_unknown EQUALS_ integer float  CR_
		{
		  Var_init[Num_Var_Init].var = $2;	
		  Var_init[Num_Var_Init].ktype = $4;
	       	  Var_init[Num_Var_Init].init_val = $5;
		  Num_Var_Init++;
		  card_read("INITIALIZE_CARD",file_index);  
		}
		| error {}
		;
		
variable_name_unknown:	
  VELOCITY1_		{$$=1;} 
| VELOCITY2_		{$$=2;} 
| VELOCITY3_		{$$=3;} 
| TEMPERATURE_		{$$=4;} 
| MASS_FRACTION_	{$$=5;} 
| MESH_DISPLACEMENT1_	{$$=6;} 
| MESH_DISPLACEMENT2_	{$$=7;} 
| MESH_DISPLACEMENT3_	{$$=8;} 
| SURFACE_    		{$$=9;} 
| PRESSURE_		{$$=10;} 
| POLYMER_STRESS11_     {$$=11;} 
| POLYMER_STRESS12_     {$$=12;}  
| POLYMER_STRESS22_     {$$=12;}  
| POLYMER_STRESS13_     {$$=13;}  
| POLYMER_STRESS23_     {$$=14;}  
| POLYMER_STRESS33_     {$$=15;}  
| SOLID_DISPLACEMENT1_  {$$=16;} 
| SOLID_DISPLACEMENT2_  {$$=17;} 
| SOLID_DISPLACEMENT3_  {$$=18;} 
| VELOCITY_GRADIENT11_  {$$=19;}  
| VELOCITY_GRADIENT12_  {$$=20;}  
| VELOCITY_GRADIENT21_  {$$=21;}  
| VELOCITY_GRADIENT22_  {$$=22;}  
| VELOCITY_GRADIENT13_  {$$=23;} 
| VELOCITY_GRADIENT23_  {$$=24;} 
| VELOCITY_GRADIENT31_  {$$=25;} 
| VELOCITY_GRADIENT32_  {$$=26;} 
| VELOCITY_GRADIENT33_  {$$=27;} 
| VOLTAGE_         	{$$=28;}  
| FILL_    		{$$=29;}  
| LS_                   {$$=29;}  
| SHEAR_RATE_           {$$=30;} 
| PVELOCITY1_           {$$=31;} 
| PVELOCITY2_           {$$=32;} 
| PVELOCITY3_           {$$=33;} 
| POLYMER_STRESS11_1_   {$$=34;}  
| POLYMER_STRESS12_1_   {$$=35;}  
| POLYMER_STRESS22_1_   {$$=36;}  
| POLYMER_STRESS13_1_   {$$=37;}  
| POLYMER_STRESS23_1_   {$$=38;}  
| POLYMER_STRESS33_1_   {$$=39;}  
| POLYMER_STRESS11_2_   {$$=40;}  
| POLYMER_STRESS12_2_   {$$=41;}  
| POLYMER_STRESS22_2_   {$$=42;}  
| POLYMER_STRESS13_2_   {$$=43;}  
| POLYMER_STRESS23_2_   {$$=44;}  
| POLYMER_STRESS33_2_   {$$=45;}  
| POLYMER_STRESS11_3_   {$$=46;}  
| POLYMER_STRESS12_3_   {$$=47;}  
| POLYMER_STRESS22_3_   {$$=48;}  
| POLYMER_STRESS13_3_   {$$=49;}  
| POLYMER_STRESS23_3_   {$$=50;}  
| POLYMER_STRESS33_3_   {$$=51;}  
| POLYMER_STRESS11_4_   {$$=52;}  
| POLYMER_STRESS12_4_   {$$=53;}  
| POLYMER_STRESS22_4_   {$$=54;}  
| POLYMER_STRESS13_4_   {$$=55;}  
| POLYMER_STRESS23_4_   {$$=56;}  
| POLYMER_STRESS33_4_   {$$=57;}  
| POLYMER_STRESS11_5_   {$$=58;}  
| POLYMER_STRESS12_5_   {$$=59;}  
| POLYMER_STRESS22_5_   {$$=60;}  
| POLYMER_STRESS13_5_   {$$=61;}  
| POLYMER_STRESS23_5_   {$$=62;}  
| POLYMER_STRESS33_5_   {$$=63;}  
| POLYMER_STRESS11_6_   {$$=64;}  
| POLYMER_STRESS12_6_   {$$=65;}  
| POLYMER_STRESS22_6_   {$$=66;}  
| POLYMER_STRESS13_6_   {$$=67;}  
| POLYMER_STRESS23_6_   {$$=68;}  
| POLYMER_STRESS33_6_   {$$=69;}  
| POLYMER_STRESS11_7_   {$$=70;} 
| POLYMER_STRESS12_7_   {$$=71;}  
| POLYMER_STRESS22_7_   {$$=72;}  
| POLYMER_STRESS13_7_   {$$=73;}  
| POLYMER_STRESS23_7_   {$$=74;}  
| POLYMER_STRESS33_7_   {$$=75;}  
| SPECIES_UNK_0_        {$$=76;}
| SPECIES_UNK_1_        {$$=77;}
| SPECIES_UNK_2_        {$$=78;}
| SPECIES_UNK_3_        {$$=79;}
| SPECIES_UNK_4_        {$$=80;}
| SPECIES_UNK_5_        {$$=81;}
| SPECIES_UNK_6_        {$$=82;}
| SPECIES_UNK_7_        {$$=83;}
| SPECIES_UNK_8_        {$$=84;}
| SPECIES_UNK_9_        {$$=85;}
| SPECIES_UNK_10_       {$$=86;}
| SPECIES_UNK_11_       {$$=87;}
| SPECIES_UNK_12_       {$$=88;}
| SPECIES_UNK_13_       {$$=89;}
| SPECIES_UNK_14_       {$$=90;}
| SPECIES_UNK_15_       {$$=91;}
| SPECIES_UNK_16_       {$$=92;}
| SPECIES_UNK_17_       {$$=93;}
| SPECIES_UNK_18_       {$$=94;}
| SPECIES_UNK_19_       {$$=95;}
| SPECIES_UNK_20_       {$$=96;}
| SPECIES_UNK_21_       {$$=97;}
| SPECIES_UNK_22_       {$$=98;}
| SPECIES_UNK_23_       {$$=99;}
| SPECIES_UNK_24_       {$$=100;}
| SPECIES_UNK_25_       {$$=101;}
| SPECIES_UNK_26_       {$$=102;}
| SPECIES_UNK_27_       {$$=103;}
| SPECIES_UNK_28_       {$$=104;}
| SPECIES_UNK_29_       {$$=105;}
| SPECIES_UNK_LAST_     {$$=105;}
| VOLF_PHASE_0_         {$$=106;}
| VOLF_PHASE_1_         {$$=107;}
| VOLF_PHASE_2_         {$$=108;}
| VOLF_PHASE_3_         {$$=109;}
| VOLF_PHASE_4_         {$$=110;}
| VOLF_PHASE_LAST_      {$$=110;}
| POR_LIQ_PRES_		{$$=111;}
| POR_GAS_PRES_		{$$=112;}
| POR_POROSITY_		{$$=113;}
| POR_SATURATION_	{$$=114;}
| POR_LAST_             {$$=114;}    
| VORT_DIR1_            {$$=115;}
| VORT_DIR2_            {$$=116;}
| VORT_DIR3_            {$$=117;}
| VORT_LAMBDA_          {$$=118;}
| EXTERNAL_             {$$=0;} 
| error {}
;

external_field_card:
		/* empty */ {}
		| EXTERNAL_ FIELD_ EQUALS_ STRING_ external_field_interpolation_method STRING_  CR_
		{
                  if(Num_Var_External > MAX_EXTERNAL_FIELD) 
	          {
	            sprintf(msg,">%d external field vars. Fix MAX_EXTERNAL_FIELD (rf_fem_const.h), recompile.",MAX_EXTERNAL_FIELD ); 
		    declare_error(msg);
	          }
	          else
	          {		
 		    efv->ev = T_SOMETHING;
		    strcpy(efv->name[Num_Var_External] , $4);
		    efv->i[Num_Var_External] = $5;
		    strcpy (efv->file_nm[Num_Var_External], $6);
		    Num_Var_External++;
		    card_read("EXTERNAL_FIELD_CARD",file_index);
		  }
		}
		| error {}
		;
		
external_field_interpolation_method:
		  Q1_		{$$=I_Q1;}
		| Q2_		{$$=I_Q2;}
		| Q2_LSA_	{$$=I_Q2_LSA;}
		| Q1_D_		{$$=I_Q1_D;}
		| Q2_D_		{$$=I_Q2_D;}
		| Q2_D_LSA_	{$$=I_Q2_D_LSA;}
		| PQ1_		{$$=I_PQ1;}
		| PQ2_		{$$=I_PQ2;}
		| P0_		{$$=I_P0;}
		| P1_		{$$=I_P1;}
		| SP_		{$$=I_SP;}
		| error {}
		;

pressure_datum_card:
		/* empty */ {}
		| PRESSURE_ DATUM_ EQUALS_ float optional_pressure_datum_units CR_
		{
		  if ($5 == 0) upd->Pressure_Datum = ($4)*1.01325E6;
		  if ($5 == 1) upd->Pressure_Datum = ($4)*1.01325E6 / 760;
		  if ($5 == 2) upd->Pressure_Datum = ($4)*1.0;		  
		  card_read("PRESSURE_DATUM_CARD",file_index);
		}
		| error {}
		;
		
optional_pressure_datum_units:
		/* empty */ 	{$$=4;}
		| ATM_		{$$=0;}
		| TORR_		{$$=1;}
		| CGS_		{$$=2;}
		| error 	{}
		;		     
		  		
		
anneal_mesh_on_output_card:
		/* empty */ {}
		| ANNEAL_ MESH_ ON_ OUTPUT_ EQUALS_ yes_no_option CR_
		{
		  Anneal_Mesh = $6;
		  card_read("ANNEAL_MESH_ON_OUTPUT_CARD",file_index); 
		}
		| error {}
		;		

time_integration_specifications_delimiter_card:
		/* empty */ {}
		| TIME_ INTEGRATION_ SPECIFICATIONS_  CR_
		{
		  card_read("TIME_INTEGRATION_SPECIFICATIONS_DELIMITER_CARD",file_index);  
		}
		| error {}
		;

time_integration_card:
		/* empty */ {}
		| TIME_ INTEGRATION_ EQUALS_ time_integraton_option CR_
		{
  		  for (j=0; j<MAX_NUMBER_MATLS; j++) 
  		  {
                    pd_glob[j]->TimeIntegration =  $4;		
                  }
		  card_read("TIME_INTEGRATION_CARD",file_index);  
		}
		| error {}
		;
		
time_integraton_option:
		  STEADY_ 	{$$=STEADY;}
		| TRANSIENT_	{$$=TRANSIENT;}
		| error		{}
		;
		
delta_t_card:
		/* empty */ {}
		| DELTA_T_ EQUALS_ float  CR_
		{
		  tran->Delta_t0 = $3;
		  card_read("DELTA_T_CARD",file_index);  
		}
		| error {}
		;		
		
maximum_number_of_time_steps_card:
		/* empty */ {}
		| MAXIMUM_ NUMBER_ OF_ TIME_ STEPS_ EQUALS_ integer  CR_
		{
    		  tran->MaxTimeSteps = $7; 		  
		  card_read("MAXIMUM_NUMBER_OF_TIME_STEPS_CARD",file_index);  
		}
		| error {}
		;
		
maximum_time_card:
		/* empty */ {}
		| MAXIMUM_ TIME_ EQUALS_ float  CR_
		{
		  tran->TimeMax = $4; 
		  card_read("MAXIMUM_TIME_CARD",file_index);  
		}
		| error {}
		;
		
maximum_time_step_card:
		/* empty */ {}
		| MAXIMUM_ TIME_ STEP_ EQUALS_ float  CR_
		{
		  tran->Delta_t_max = $5;
		  card_read("MAXIMUM_TIME_STEP_CARD",file_index);  
		}
		| error {}
		;
		
minimum_time_step_card:
		/* empty */ {}
		| MINIMUM_ TIME_ STEP_ EQUALS_ float  CR_
		{
		  tran->Delta_t_min = $5;
		  card_read("MINIMUM_TIME_STEP_CARD",file_index);
		}
		| error {}
		;
		
time_step_parameter_card:
		/* empty */ {}
		| TIME_ STEP_ PARAMETER_ EQUALS_ float  CR_
		{
		  tran->theta = $5;
		  card_read("TIME_STEP_PARAMETER_CARD",file_index); 
		}
		| error {}
		;
		
time_step_error_card:
		/* empty */ {}
		| TIME_ STEP_ ERROR_ EQUALS_ float  CR_
		{ 
		  tran->eps = $5;
		  card_read("TIME_STEP_ERROR_CARD",file_index);  
		}
		| TIME_ STEP_ ERROR_ EQUALS_ float integer integer integer integer integer integer integer CR_
		{
		  tran->eps = $5;
		  tran->use_var_norm[0] = $6;
		  tran->use_var_norm[1] = $7;
		  tran->use_var_norm[2] = $8;
		  tran->use_var_norm[3] = $9;
		  tran->use_var_norm[4] = $10;
		  tran->use_var_norm[5] = $11;
		  tran->use_var_norm[6] = $12; 
		  card_read("TIME_STEP_ERROR_CARD",file_index);  
		}
		| error {}
		;
		
printing_frequency_card:
		/* empty */ {}
		| PRINTING_ FREQUENCY_ EQUALS_ integer optional_float CR_
		{	
                  tran->print_freq = $4;
                  if (tran->print_freq == 0)
                  {
                    if ( $5 == NA )
                    {
                      sprintf(msg,"Error reading Printing delta time." ); 
		      declare_error(msg);
		    }
		    else
		    {
		      tran->print_delt = $5;
                      print_delt2 = -$5;
                      print_delt2_time = TimeMax;
                    }
                  }
		  card_read("PRINTING_FREQUENCY_CARD",file_index);  
		}
		| error {}
		;
		
second_frequency_time_card:
		/* empty */ {}
		| SECOND_ FREQUENCY_ TIME_ EQUALS_ float float CR_
		{
                  if(tran->print_freq == 0)
                  {
                    tran->print_delt2_time = $5;
		    tran->print_delt2 = $6;
		  }
		  card_read("SECOND_FREQUENCY_TIME_CARD",file_index);
		}
		| error {}
		;
		
initial_time_card:
		/* empty */ {}
		| INITIAL_ TIME_ EQUALS_ float CR_
		{
		  tran->init_time = $4;
		  card_read("INITIAL_TIME_CARD",file_index); 
		}
		| error {}
		;
		
fill_subcycle_card:
		/* empty */ {}
		| FILL_ SUBCYCLE_ EQUALS_ integer  CR_
		{
		  tran->exp_subcycle = $4;
		  card_read("FILL_SUBCYCLE_CARD",file_index);  
		}
		| error {}
		;
		
level_set_interface_tracking_card:
		/* empty */ {}
		| LEVEL_ SET_ INTERFACE_ TRACKING_ EQUALS_ yes_no_option CR_
		{
		  if ($6)
		  {
		    ls = (struct Level_Set_Data *) array_alloc(1, 1, sizeof( struct Level_Set_Data ) );
	    	    ls->embedded_bc = NULL;
                    ls->init_surf_list = NULL;
                    lsi = (struct Level_Set_Interface *) array_alloc(1, 1, sizeof( struct Level_Set_Interface ) );
                    ls->Use_Level_Set = Use_Level_Set = TRUE ;
	            zero_lsi();
	          }
		  card_read("LEVEL_SET_INTERFACE_TRACKING_CARD",file_index);                                      		
		}
		| error {}
		;

level_set_length_scale_card:
		/* empty */ {}
		| LEVEL_ SET_  LENGTH_ SCALE_ EQUALS_ float CR_
		{
		  if (ls != NULL)
		  {
		    ls->Length_Scale = $6;
		  }
		  card_read("LEVEL_SET_INTERFACE_TRACKING_CARD",file_index);  
		}
		| error {}
		;		
		
level_set_initialization_method_card:
		/* empty */ {}
		| LEVEL_ SET_ INITIALIZATION_ METHOD_ EQUALS_ PROJECT_ CR_
		{
		if (ls != NULL)
		  {
                    ls->Init_Method == PROJECT;
                    card_read("LEVEL_SET_INITIALIZATION_METHOD_CARD",file_index);
                  }
                }
		| LEVEL_ SET_ INITIALIZATION_ METHOD_ EQUALS_ EXODUS_ CR_                
                {
                  if (ls != NULL)
		  {
                    ls->Init_Method == EXO_READ;
                    card_read("LEVEL_SET_INITIALIZATION_METHOD_CARD",file_index);
                  }
                }
                | LEVEL_ SET_ INITIALIZATION_ METHOD_ EQUALS_ SURFACES_ integer CR_
                {
                  if (ls != NULL)
		  {
		    ls->Init_Method = SURFACES; 
		    ls->init_surf_list = create_surf_list();
		    
		    number_of_level_set_initialization_surface_cards_expected = $7;
		    accept_level_set_initialization_surface_cards = TRUE;
                    card_read("LEVEL_SET_INITIALIZATION_METHOD_CARD",file_index);
                    
                  }
                }
		| LEVEL_ SET_ INITIALIZATION_ METHOD_ EQUALS_ NODESET_ NS_ integer EB_ integer CR_
		{
		  struct LS_Surf *surf;
                  struct LS_Surf_NS_Data *s;
		  ls->Init_Method = SURFACES;
		  ls->init_surf_list = create_surf_list();
                  surf = create_surf( LS_SURF_NS );
                  append_surf( ls->init_surf_list, surf );
                  s = (struct LS_Surf_NS_Data *) surf->data;
                  s->ns_id = $8;
                  s->PosEB_id = $10;
		  card_read("LEVEL_SET_INITIALIZATION_METHOD_CARD",file_index);
		}
		| error {}
		;		  
		  
		  
		  
level_set_initialization_surface_card:
		/* empty */ {}
		| SURF_ EQUALS_ PLANE_ float float float float CR_
		  {
		    if ( accept_level_set_initialization_surface_cards )
		    {
		     /* struct LS_Surf *surf;               NOTE: I WAS NOT ABLE TO GET THE ACTIONS OF THE SURFACE CARDS TO COMPILE.
                      struct LS_Surf_Plane_Data *s;               INSTEAD OF WASTING MORE TIME ON IT, I COMMENTED THEM OUT WITH HOPES
                      surf = create_surf( LS_SURF_PLANE );        THAT SOMEONE MORE FAMILIAR WITH THEM WOULD BE ABLE TO GET THEM TO WORK.
                      append_surf( ls->init_surf_list, surf );    JSS - 3/13/02
                      s = (struct LS_Surf_Plane_Data *) surf->data;
                      s->n[0] = $4; 
		      s->n[1] = $5; 
		      s->n[2] = $6;
                      s->d    = $7; */
		      number_of_level_set_initialization_surface_cards_found++;
		      if( number_of_level_set_initialization_surface_cards_expected == number_of_level_set_initialization_surface_cards_found)
		        accept_level_set_initialization_surface_cards = FALSE;  
		          card_read("LEVEL_SET_INITIALIZATION_SURFACE_CARD",file_index);    
		  } 
		}
		| SURF_ EQUALS_ CIRCLE_ float float float CR_
 		  {
		    if ( accept_level_set_initialization_surface_cards )
		    {
		      /*  struct LS_Surf *surf;
		      struct LS_Surf_Plane_Data *s;
		      surf = create_surf( LS_SURF_CIRCLE );
                      append_surf( ls->init_surf_list, surf );
                      s = (struct LS_Surf_Sphere_Data *) surf->data;
                      s->center[2] = 0.;  
                      s->center[0] = $4;
		      s->center[1] = $5;
		      s->r         = $6;  */
		      number_of_level_set_initialization_surface_cards_found++;
		      if( number_of_level_set_initialization_surface_cards_expected == number_of_level_set_initialization_surface_cards_found)
		        accept_level_set_initialization_surface_cards = FALSE;	       
		      card_read("LEVEL_SET_INITIALIZATION_SURFACE_CARD",file_index);             
		    }
 		  }
		| SURF_ EQUALS_ SPHERE_ float float float float CR_
		  {
		    if ( accept_level_set_initialization_surface_cards )
		    {
                      /* struct LS_Surf *surf;
   		      struct LS_Surf_Sphere_Data *s;
                      surf = create_surf( LS_SURF_SPHERE );
                      append_surf( ls->init_surf_list, surf );
                      s = (struct LS_Surf_Sphere_Data *) surf->data;
                      s->center[0] = $4;
		      s->center[1] = $5;
                      s->center[2] = 0;
		      s->r = $6; */
		      number_of_level_set_initialization_surface_cards_found++;
		      if( number_of_level_set_initialization_surface_cards_expected == number_of_level_set_initialization_surface_cards_found)
		        accept_level_set_initialization_surface_cards = FALSE;
		      card_read("LEVEL_SET_INITIALIZATION_SURFACE_CARD",file_index);
		    }		    
		  } 
		| SURF_ EQUALS_ SS_ integer CR_
		  {
		    
		    if ( accept_level_set_initialization_surface_cards )
		    {
		     /* struct LS_Surf *surf;
 	              struct LS_Surf_SS_Data *s;
                      surf = create_surf( LS_SURF_SS );
                      append_surf( ls->init_surf_list, surf );
                      s = (struct LS_Surf_SS_Data *) surf->data;
                      s->ss_id = $4;*/
		      number_of_level_set_initialization_surface_cards_found++;
		      if( number_of_level_set_initialization_surface_cards_expected == number_of_level_set_initialization_surface_cards_found)
		        accept_level_set_initialization_surface_cards = FALSE;
                      card_read("LEVEL_SET_INITIALIZATION_SURFACE_CARD",file_index); 
		    }		    
		  }
		| SURF_ EQUALS_ USER_ floating_point_constant_list CR_
		  {
		    
		    if ( accept_level_set_initialization_surface_cards )
		    {
		      /*struct LS_Surf *surf;
		      struct LS_Surf_Sphere_Data *s;
                      surf = create_surf( LS_SURF_USER );
                      append_surf( ls->init_surf_list, surf );
                      s = (struct LS_Surf_User_Data *) surf->data;
 		      s->Int_Data[0] = read_floats( &(s->Real_Data), 0);         	    */
		      number_of_level_set_initialization_surface_cards_found++;
		      if( number_of_level_set_initialization_surface_cards_expected == number_of_level_set_initialization_surface_cards_found)
		        accept_level_set_initialization_surface_cards = FALSE;
		    }		    
		  card_read("LEVEL_SET_INITIALIZATION_SURFACE_CARD",file_index);      
		  }  
		| SURF_ EQUALS_ ISOSURFACE_  var_name CR_
		  {
		    
		    if ( accept_level_set_initialization_surface_cards )
		    {                     	 
		      /*struct LS_Surf *surf;   
	              struct LS_Surf_Iso_Data *s;
                      surf = create_surf( LS_SURF_ISOSURFACE );
                      append_surf( ls->init_surf_list, surf );
                      s = (struct LS_Surf_Iso_Data *) surf->data;
                      s->isovar = $4;		  */
                      number_of_level_set_initialization_surface_cards_found++;
		      if( number_of_level_set_initialization_surface_cards_expected == number_of_level_set_initialization_surface_cards_found)
		        accept_level_set_initialization_surface_cards = FALSE;
		    }		    
		  card_read("LEVEL_SET_INITIALIZATION_SURFACE_CARD",file_index);      
		  }		  
		| error {}
		;
		
level_set_control_width_card:
		/* empty */ {}
		| LEVEL_ SET_ CONTROL_ WIDTH_ EQUALS_ float CR_
		{
		  if (ls != NULL)
		  {
		    ls->Control_Width = $6;
		    if ( ls->Control_Width > 10.0)
		    {
		      sprintf(msg, " That's an awfully large Level Set Control Width."); 
		      declare_warning(msg);		    	    
		    }
		  }
		  card_read("LEVEL_SET_CONTROL_WIDTH_CARD",file_index);  
		}
		| error {}
		;

level_set_renormalization_tolerance_card:
		/* empty */ {}
		| LEVEL_ SET_ RENORMALIZATION_ TOLERANCE_ EQUALS_ float CR_
		{
		  if (ls != NULL)
		  {
		    ls->Renorm_Tolerance = $6;
		    if ( ls->Renorm_Tolerance > 1.0)
		    {
		      sprintf(msg, " That's an awfully large Level Set Renormailzation Tolerance."); 
		      declare_warning(msg);		    	    
		    }		    
		  }
		  card_read("LEVEL_SET_RENORMALIZATION_TOLERANCE_CARD",file_index);  
		}
		| error {}
		;
		
level_set_renormalization_method_card:
		/* empty */ {}
		| LEVEL_ SET_ RENORMALIZATION_ METHOD_ EQUALS_ level_set_remornalization_method optional_float CR_
		{
		  if (ls != NULL)
		  {
		    ls->Renorm_Method = $6;
		    if ( ls->Renorm_Method == HUYGENS_C )
		    {
		      if ( $7 != NA )
		      {
		        ls->Mass_Value = $7;
		        ls->Mass_Sign = ls->Mass_Value <= 0 ? I_NEG_FILL : I_POS_FILL;
		        ls->Mass_Value = fabs( ls->Mass_Value);  
		      }
		      else
		      {
		        ls->Mass_Value = 0.0;
		        ls->Mass_Sign  = I_NEG_FILL;
		      }
		    }
		  }
		  card_read("LEVEL_SET_RENORMALIZATION_METHOD_CARD",file_index); 
		}
		| error {}
		;
		
level_set_remornalization_method:
                  CORRECTION_		{$$=CORRECT;}
                | HUYGENS_		{$$=HUYGENS;}
                | HUYGENS_CONSTRAINED_	{$$=HUYGENS_C;}
                | NONE_			{$$=FALSE;}
                | error {}
                ;		
		
level_set_renormalization_frequency_card:
		/* empty */ {}
		| LEVEL_ SET_ RENORMALIZATION_ FREQUENCY_ EQUALS_ float CR_
		{
		  if (ls != NULL)
		  {
		    ls->Renorm_Countdown = ls->Renorm_Freq = $6;
		  }
		  card_read("LEVEL_SET_RENORMALIZATION_FREQUENCY_CARD",file_index); 
		}
		| error {}
		;		
		
continuation_specificatins_delimiter_card:
		/* empty */ {}
		| CONTINUATION_ SPECIFICATIONS_  CR_
		{
		  card_read("CONTINUATION_SPECIFICATIONS_DELIMITER_CARD",file_index);  
		}
		| error {}
		;
		
		
continuation_card:
		/* empty */ {}
		| CONTINUATION_ EQUALS_  continuation_method CR_
		{
		  for (j=0; j<MAX_NUMBER_MATLS; j++) 
		  {
		    pd_glob[j]->Continuation =  $3;
		  }
		  card_read("CONTINUATION_CARD",file_index);  
		}
		| error {}
		;
		
continuation_method:
		  NONE_ 	{$$=-1;}
		| ZERO_ 	{$$=ALC_ZEROTH;}
		| FIRST_ 	{$$=ALC_FIRST;}
		| HZERO_ 	{$$=HUN_ZEROTH;}
		| HFIRST_ 	{$$=HUN_FIRST;}
		| SECOND_	{$$=ALC_SECOND;}
		| LOCA_		{$$=LOCA;}
		| error {}
		;

continuation_type_card:
		/* empty */ {}
		| CONTINUATION_ TYPE_ EQUALS_ continuation_type CR_
		{                
		  if ( pd_glob[0]->Continuation == LOCA )
		  {               
		    cont->upType = ContType;
		  }
		  if ($4 == -1 )
		  {
		    declare_warning("No continuation.");
		  }
		  card_read("CONTINUATION_TYPE_CARD",file_index);  
		}
		| error {}
		;
		
continuation_type:
		  NONE_ {$$=-1;}
		| BC_ 	{$$=1;}
		| MT_ 	{$$=2;}
		| error {}
		;

boundary_condition_id_card:
		/* empty */ {}
		| BOUNDARY_ CONDITION_ ID_ EQUALS_ integer CR_
		{
		  if ( pd_glob[0]->Continuation == LOCA )
		  {		
		    cont->upBCID = $5;
		  }
		  card_read("BOUNDARY_CONDITION_ID_CARD",file_index);  
		}
		| error {}
		;

boundary_condition_data_float_tag_card:
		/* empty */ {}
		| BOUNDARY_ CONDITION_ DATA_ THE_WORD_FLOAT_ TAG_ EQUALS_ integer CR_
		{
		  if ( pd_glob[0]->Continuation == LOCA )
		  {
		    cont->upDFID = $7;
		  }
		  card_read("BOUNDARY_CONDITION_DATA_FLOAT_TAG_CARD",file_index);  
		}
		| error {}
		;


material_id_card:
		/* empty */ {}
		| MATERIAL_ ID_ EQUALS_  integer CR_
		{
		  if ( pd_glob[0]->Continuation == LOCA )
		  {
		    cont->upMTID = ($4)-1;
		  }
		  card_read("MATERIAL_ID_CARD",file_index);  
		}
		| error {}
		;


material_property_tag_card:
		/* empty */ {}
		| MATERIAL_ PROPERTY_ TAG_ EQUALS_ integer CR_
		{
		  if ( pd_glob[0]->Continuation == LOCA )
		  {
		    cont->upMPID = $5;
		  }
		  card_read("MATERIAL_PROPERTY_TAG_CARD",file_index);  
		}
		| error {}
		;


material_property_tag_subindex_card:
		/* empty */ {}
		| MATERIAL_ PROPERTY_ TAG_ SUBINDEX_ EQUALS_  integer CR_
		{
		  if ( pd_glob[0]->Continuation == LOCA )
		  {
		    cont->upMDID = $6;
		  }
		  card_read("MATERIAL_PROPERTY_TAG_SUBINDEX_CARD",file_index);  
		}
		| error {}
		;


initial_parameter_value_card:
		/* empty */ {}
		| INITIAL_ PARAMETER_ VALUE_ EQUALS_ float CR_
		{
		  if ( pd_glob[0]->Continuation == LOCA )
		  {
		    cont->BegParameterValue = $5;
		  }
		  card_read("INITIAL_PARAMETER_VALUE_CARD",file_index);  
		}
		| error {}
		;


final_parameter_value_card:
		/* empty */ {}
		| FINAL_ PARAMETER_ VALUE_ EQUALS_ float CR_
		{
		  if ( pd_glob[0]->Continuation == LOCA )
		  {
		    cont->EndParameterValue = $5;
		  }
		  card_read("FINAL_PARAMETER_VALUE_CARD",file_index);  
		}
		| error {}
		;


delta_s_card:
		/* empty */ {}
		| DELTA_S_ EQUALS_  float CR_
		{
		  if ( pd_glob[0]->Continuation == LOCA )
		  {
		    cont->Delta_s0 = $3;
		  }
		  card_read("DELTA_S_CARD",file_index);  
		}
		| error {}
		;


maximum_number_of_path_steps_card:
		/* empty */ {}
		|MAXIMUM_ NUMBER_ OF_ PATH_ STEPS_ EQUALS_  integer CR_
		{
		  if ( pd_glob[0]->Continuation == LOCA )
		  {
		    cont->MaxPathSteps = $7;
		  }
		  card_read("MAXIMUM_NUMBER_OF_PATH_STEPS_CARD",file_index);  
		}
		| error {}
		;

maximum_path_value_card:
		/* empty */ {}
		| MAXIMUM_ PATH_ VALUE_ EQUALS_  float CR_
		{
		  sprintf(msg, "The continuation maximum path value card is no longer implemented." ); 
		  declare_warning(msg); 
		 /*
		  if ( pd_glob[0]->Continuation == LOCA )
		  {
		    cont->PathMax = $5;
		  }
		    This card is currently commented out in the old parser
		    */
		  card_read("MAXIMUM_PATH_VALUE_CARD",file_index);  
		}
		| error {}
		;

minimum_path_step_card:
		/* empty */ {}
		| MINIMUM_ PATH_ STEP_ EQUALS_  float CR_
		{
		  if ( pd_glob[0]->Continuation == LOCA )
		  {
		    cont->Delta_s_min = $5;
		  }
		  card_read("MINIMUM_PATH_STEP_CARD",file_index);  
		}
		| error {}
		;

maximum_path_step_card:
		/* empty */ {}
		| MAXIMUM_ PATH_ STEP_ EQUALS_  float CR_
		{
		  if ( pd_glob[0]->Continuation == LOCA )
		  {
		    cont->Delta_s_max = $5;
		  }
		  card_read("MAXIMUM_PATH_STEP_CARD",file_index);  
		}
		| error {}
		;

path_step_parameter_card:
		/* empty */ {}
		| PATH_ STEP_ PARAMETER_ EQUALS_ float CR_
		{
		  sprintf(msg, "The continuation path step parameter card is no longer implemented." ); 
		  declare_warning(msg); 
		  /*
		  if ( pd_glob[0]->Continuation == LOCA )
		  {
		    cont->theta = $5;
		  }
		  */
		  card_read("PATH_STEP_PARAMETER_CARD",file_index);  
		}
		| error {}
		;

path_step_error_card:
		/* empty */ {}
		| PATH_ STEP_ ERROR_ EQUALS_  float integer integer integer integer integer integer CR_
		{
	          sprintf(msg, "The continuation path step error card is no longer implemented." ); 
		  declare_warning(msg); 
		  /*
		  if ( pd_glob[0]->Continuation == LOCA )
		  {
		    cont->eps = $5;
		    cont->use_var_norm[0] = $6;
		    cont->use_var_norm[1] = $7;
		    cont->use_var_norm[2] = $8;
		    cont->use_var_norm[3] = $9;
		    cont->use_var_norm[4] = $10;
		    cont->use_var_norm[5] = $11;
		  }
		  */		
		  card_read("PATH_STEP_ERROR_CARD",file_index);  
		}
		| error {}
		;

continuation_printing_frequency_card:
		/* empty */ {}
		| CONTINUATION_ PRINTING_ FREQUENCY_ EQUALS_  integer optional_float CR_
		{
		  if ( pd_glob[0]->Continuation == LOCA )
		  {
		    if ($5 == 0)
		    {
		      if( $6 != NA )
		      {
		        cont->print_delt = $6;
			cont->print_delt2 = -(cont->print_delt);
			print_delt2_time = PathMax;
		      }
		      else
		      {
		        sprintf(msg, "Error reading second frequency path step." ); 
		        declare_error(msg); 
		      }
		    }
		    else
		    {
		      cont->print_freq = $5;
		    }
		  }
		  card_read("CONTINUATION_PRINTING_FREQUENCY_CARD",file_index);  
		}
		| error {}
		;
		
second_frequency_card:
		/* empty */ {}
		| SECOND_ FREQUENCY_ EQUALS_ float float  CR_
		{
		  if ( pd_glob[0]->Continuation == LOCA )
		  {
		    if (  cont->print_freq == 0)
		    {
	  	      cont->print_delt2 = $5;
	  	      cont->print_delt2_path = $4;		  
		    }
		  }
    	          card_read("SECOND_FREQUENCY_CARD",file_index);  		  
		}
		| error {}
		;


loca_method_card:
		/* empty */ {}
		| LOCA_ METHOD_ EQUALS_ loca_method CR_
		{
		  if ( pd_glob[0]->Continuation == LOCA )
		  {
		    switch($4)
                    { 
                      case CONTINUATION:	/* SS */
                      {
                        loca_in->Cont_Alg = CONTINUATION; /* read continuation order card */
                        break;
                      }
                      case 0: 		/* zero */
                      {
                        break;		/* nothing set? */
                      }
                      case 1: 		/* first */
                      {
                        loca_in->Cont_Order = 1;
                        break;
                      }
                      case 2: 		/* alc */
                      {
                        loca_in->Cont_Order = 2;
                        break;
                      }
                      case TP_CONTINUATION: /*tp */
                      {
                        loca_in->Cont_Alg = TP_CONTINUATION;
                        break;
                      }
                      case PF_CONTINUATION: /* pf */
                      {
                        loca_in->Cont_Alg = PF_CONTINUATION;
                        break;
                      }                                                              
                    } /* end switch */                                         
                  }
		  card_read("LOCA_METHOD_CARD",file_index);  
		}
		| error {}
		;
		
loca_method:
		  SS_ 		{$$=CONTINUATION;}
		| ZERO_		{$$=0;}
		| FIRST_ 	{$$=1;}
		| ALC_		{$$=2;}
		| TP_		{$$=TP_CONTINUATION;}
		| PF_		{$$=PF_CONTINUATION;}
		| error {}
		;
		
continuation_order_card:
		/* empty */ {}
		| CONTINUATION_ ORDER_ EQUALS_ integer CR_
		{
		  if ( pd_glob[0]->Continuation == LOCA )
		  {
		    if ( loca_in->Cont_Alg == CONTINUATION)
		    {
		      loca_in->Cont_Order = $4;
		    }
		  }
		  card_read("CONTINUATION_ORDER_CARD",file_index);  
		}
		| error {}
		;

step_control_aggressiveness_card:
		/* empty */ {}
		| STEP_ CONTROL_ AGGRESSIVENESS_ EQUALS_  float CR_
		{
		  if ( pd_glob[0]->Continuation == LOCA )
		  {
		    loca_in->StepAggr = $5;
		    if (loca_in->StepAggr < -1.0e-6)
            	    {
                      loca_in->StepAggr = 0.0;
                      if (loca_in->Cont_Order != 2) cont->Delta_s0 = 
                        (cont->EndParameterValue - cont->BegParameterValue)
                        / ((double)(cont->MaxPathSteps - 1));
                    }
                  }
		  card_read("STEP_CONTROL_AGGRESSIVENESS_CARD",file_index);  
		}
		| error {}
		;

alc_desired_solution_fraction_card:
		/* empty */ {}
		| ALC_ DESIRED_ SOLUTION_ FRACTION_ EQUALS_  float CR_
		{
		  if ( pd_glob[0]->Continuation == LOCA )
		  {
		    loca_in->DpDs2 = $6;
		  }
		  card_read("ALC_DESIRED_SOLUTION_FRACTION_CARD",file_index);  
		}
		| error {}
		;

alc_max_parameter_sensitivity_card:
		/* empty */ {}
		| ALC_ MAXIMUM_ PARAMETER_ SENSITIVITY_ EQUALS_  float CR_
		{
		  if ( pd_glob[0]->Continuation == LOCA )
		  {
		    loca_in->DpDsHi = $6;
		  }
		  card_read("ALC_MAX_PARAMETER_SENSITIVITY_CARD",file_index);  
		}
		| error {}
		;
		
alc_tangent_factor_exponent_card:
		/* empty */ {}
		| ALC_ TANGENT_ FACTOR_ EXPONENT_ EQUALS_  float CR_
		{
		  if ( pd_glob[0]->Continuation == LOCA )
		  {
		    loca_in->Texp = $6;
		  }
		  card_read("ALC_TANGENT_FACTOR_EXPONENT_CARD",file_index);  
		}
		| error {}
		;
		
alc_tangent_factor_step_limit_card:
		/* empty */ {}
		| ALC_ TANGENT_ FACTOR_ STEP_ LIMIT_ EQUALS_ float CR_
		{
		  if ( pd_glob[0]->Continuation == LOCA )
		  {
		    loca_in->MaxTS = $7;
		  }
		  card_read("ALC_TANGENT_FACTOR_STEP_LIMIT_CARD",file_index);  
		}
		| error {}
		;
		
number_of_continuation_conditions_card:
		/* empty */ {}
		| NUMBER_ OF_ CONTINUATION_ CONDITIONS_ EQUALS_  integer CR_
		{
		  if ( pd_glob[0]->Continuation == LOCA )
		  {
		    if( $6 >0  )
		    {
		      number_of_ccs_expected = $6;
		      accept_ccs_until_end_of_cc_card = FALSE;
		    } 
		    else if ($6 == -1)
		    {
		      accept_ccs_until_end_of_cc_card = TRUE;
		    }
                    else if ($6 == -2)
		    {
		      accept_ccs_until_end_of_cc_card = FALSE;
		      number_of_ccs_expected = 0;
		      accept_hunting_conditions = TRUE;
		    }		    
		    else
		    {
		      declare_error("Invalid Number of Augmenting Conditions");   
		      number_of_ccs_expected = 0;
		      accept_ccs_until_end_of_cc_card = FALSE;
		    }		  
		  }
		  card_read("NUMBER_OF_CONTINUATION_CONDITIONS_CARD",file_index);  
		}
		| error {}
		;
		
cc_card:
		/* empty */ {}
		| CC_ EQUALS_  cc_type integer integer integer float float CR_
		{
		  if ( pd_glob[0]->Continuation == LOCA )
		  {
 		    common_CC_processing( $3, $4, $5, $6, $7, $8 );		  
		  }
		  card_read("CC_CARD",file_index);  
		}
		| error {}
		;
		
cc_type:
		  BC_	{$$=1;}
		| MT_  	{$$=2;}
		| error {}
		;		
		
end_of_cc_card:
		/* empty */ {}
		| END_ OF_ CC_  CR_
		{
		  if ( pd_glob[0]->Continuation == LOCA )
		  {
		    if ( !accept_ccs_until_end_of_cc_card && (number_of_ccs_expected !=  number_of_ccs_found) )
		    { 
		      sprintf(msg, "Wrong number of continuation conditions found (%i expected, %i found).", number_of_ccs_expected, number_of_ccs_found ); 
		      declare_error(msg); 
	  	    }	  	  		
                    accept_ccs_until_end_of_cc_card = FALSE; 
		  }
		  card_read("END_OF_CC_CARD",file_index);  
		}
		| error {}
		;
		
tp_continuation_type_card:
		/* empty */ {}
		| TP_ CONTINUATION_ TYPE_ EQUALS_ continuation_type CR_
		{
		  if ( (pd_glob[0]->Continuation == LOCA) && (loca_in->Cont_Alg == TP_CONTINUATION) )
		  {
		    loca_in->TPupType = $5;
		    if ($5 == -1)
		    {
		      declare_error("TP Continuation type of BC or MT required.");
		    }
		  }
		  card_read("TP_CONTINUATION_TYPE_CARD",file_index);  
		}
		| error {}
		;
				
tp_boundary_condition_id_card:
		/* empty */ {}
		| TP_ BOUNDARY_ CONDITION_ ID_ EQUALS_  integer CR_
		{
		  if ( (pd_glob[0]->Continuation == LOCA) && (loca_in->Cont_Alg == TP_CONTINUATION) )
		  {
		    loca_in->TPupBCID = $6;
		  }
		  card_read("TP_BOUNDARY_CONDITION_ID_CARD",file_index);  
		}
		| error {}
		;
		
tp_bc_data_float_tag_card:
		/* empty */ {}
		| TP_ BC_ DATA_ THE_WORD_FLOAT_ TAG_ EQUALS_  integer CR_
		{
		  if ( (pd_glob[0]->Continuation == LOCA) && (loca_in->Cont_Alg == TP_CONTINUATION) )
		  {
		    loca_in->TPupDFID = $7;
		  }
		  card_read("TP_BC_DATA_FLOAT_TAG_CARD",file_index);  
		}
		| error {}
		;
		
tp_parameter_material_id_card:
		/* empty */ {}
		| TP_ PARAMETER_ MATERIAL_ ID_ EQUALS_  integer CR_
		{
		  if ( (pd_glob[0]->Continuation == LOCA) && (loca_in->Cont_Alg == TP_CONTINUATION) )
		  {
		    loca_in->TPupMTID = ($6)-1;
		  }
		  card_read("TP_PARAMETER_MATERIAL_ID_CARD",file_index);  
		}
		| error {}
		;

tp_parameter_material_property_tag_card:
		/* empty */ {}
		| TP_ PARAMETER_ MATERIAL_ PROPERTY_ TAG_ EQUALS_ integer CR_
		{
		  if ( (pd_glob[0]->Continuation == LOCA) && (loca_in->Cont_Alg == TP_CONTINUATION) )
		  {
		    loca_in->TPupMPID = $7;
		  }
		  card_read("TP_PARAMETER_MATERIAL_PROPERTY_TAG_CARD",file_index);  
		}
		| error {}
		;
		
tp_material_property_tag_subindex_card:
		/* empty */ {}
		| TP_ MATERIAL_ PROPERTY_ TAG_ SUBINDEX_ EQUALS_ integer CR_
		{
		  if ( (pd_glob[0]->Continuation == LOCA) && (loca_in->Cont_Alg == TP_CONTINUATION) )
		  {
		    loca_in->TPupMDID = $7;
		  }
		  card_read("TP_MATERIAL_PROPERTY_TAG_SUBINDEX_CARD",file_index);  
		}
		| error {}
		;
		
initial_guess_of_tp_parameter_card:
		/* empty */ {}
		| INITIAL_ GUESS_ OF_ TP_ PARAMETER_ EQUALS_  float CR_
		{
		  if ( (pd_glob[0]->Continuation == LOCA) && (loca_in->Cont_Alg == TP_CONTINUATION) )
		  {
		    loca_in->TPGuess = $7;
		  }
		  card_read("INITIAL_GUESS_OF_TP_PARAMETER_CARD",file_index);  
		}
		| error {}
		;
		
tp_parameter_final_value_card:
		/* empty */ {}
		| TP_ PARAMETER_ FINAL_ VALUE_ EQUALS_  float CR_
		{
		  if ( (pd_glob[0]->Continuation == LOCA) && (loca_in->Cont_Alg == TP_CONTINUATION) )
		  {
		    loca_in->TPFinal = $6;
		  }
		  card_read("TP_PARAMETER_FINAL_VALUE_CARD",file_index);  
		}
		| error {}
		;
		
number_of_tp_continuation_conditions_card:
		/* empty */ {}
		| NUMBER_ OF_ TP_ CONTINUATION_ CONDITIONS_ EQUALS_  integer CR_
		{
		  if ( (pd_glob[0]->Continuation == LOCA) && (loca_in->Cont_Alg == TP_CONTINUATION) )
		  {
                    if( $7 >0  )
		    {
		      number_of_tps_expected = $7;
		      accept_tps_until_end_of_tp_card = FALSE;
		    } 
		    else if ($7 == -1)
		    {
		      accept_tps_until_end_of_tp_card = TRUE;
		    }		    
		    else
		    {
		      declare_error("Invalid Number of TP Continuation Conditions");   
		      number_of_tps_expected = 0;
		      accept_tps_until_end_of_tp_card = FALSE;
		    }			  
		  
		  }
		  card_read("NUMBER_OF_TP_CONTINUATION_CONDITIONS_CARD",file_index);  
		}
		| error {}
		;
				
tc_card:
		/* empty */ {}
		| TC_ EQUALS_  cc_type integer integer integer float float CR_
		{
		  if ( (pd_glob[0]->Continuation == LOCA) && (loca_in->Cont_Alg == TP_CONTINUATION) )
		  {
 		    common_TP_processing( $3, $4, $5, $6, $7, $8 );		    
		  }
		  card_read("TC_CARD",file_index);  
		}
		| error {}
		;
		
end_of_tc_card:	
		/* empty */ {}
		| END_ OF_ TC_  CR_
		{
		  if ( (pd_glob[0]->Continuation == LOCA) && (loca_in->Cont_Alg == TP_CONTINUATION) )
		  {
                    if ( !accept_tps_until_end_of_tp_card && (number_of_tps_expected !=  number_of_tps_found) )
		    { 
		      sprintf(msg, "Wrong number of TP continuation conditions found (%i expected, %i found).", number_of_tps_expected, number_of_tps_found ); 
		      declare_error(msg); 
	  	    }	  	  		
                    accept_tps_until_end_of_tp_card = FALSE; 		  
		  }
		  card_read("END_OF_TC_CARD",file_index);  
		}
		| error {}
		;	
				

hunting_conditions_specificatins_delimiter_card:
		/* empty */ {}
		| HUNTING_ SPECIFICATIONS_  CR_
		{
		  card_read("HUNTING_SPECIFICATIONS_DELIMITER_CARD",file_index);  
		}
		| error {}
		;
		
number_of_hunting_conditions_card:
		/* empty */ {}
		| NUMBER_ OF_ HUNTING_ CONDITIONS_ EQUALS_ integer CR_
		{
		  if ( !( (pd_glob[0]->Continuation != HUN_ZEROTH) 
		           && (pd_glob[0]->Continuation != HUN_FIRST) 
		           && (pd_glob[0]->Continuation != LOCA) ) 
		           && accept_hunting_conditions )		
		  {
		    if( $6 >0 )
		    {
		      number_of_hcs_expected = $6;
		      accept_hcs_until_end_of_hc_card = FALSE;
		    } 
		    else if ($6 == -1)
		    {
		      accept_hcs_until_end_of_hc_card = TRUE;
		    }		    
		    else
		    {
		      declare_error("Invalid Number of Hunting Conditions");   
		      number_of_hcs_expected = 0;
		      accept_hcs_until_end_of_hc_card = FALSE;
		    }			  
		  }		
		  card_read("NUMBER_OF_HUNTING_CONDITIONS_CARD",file_index); 		
		}
		| error {}
		;
		
hc_card:
		/* empty */ {}
		| HC_ EQUALS_ hc_type integer integer integer float float float float float CR_
		{
		  if ( !( (pd_glob[0]->Continuation != HUN_ZEROTH) 
		           && (pd_glob[0]->Continuation != HUN_FIRST) 
		           && (pd_glob[0]->Continuation != LOCA) )
		           && accept_hunting_conditions )		
		  {
		    common_HC_processing( $3, $4, $5, $6, $7, $8, $9, $10, $11 );
		  }
		  else
		  {
		    declare_warning("Hunting condition card being ignored.");
		  }
		  card_read("HC_CARD",file_index); 
		}
		| error {}
		;		

hc_type:
		  AC_	{$$=3;}
		| BC_	{$$=1;}
		| MT_	{$$=2;}
		| error {}
		;

end_of_hc_card:
		/* empty */ {}
		| END_ OF_ HC_ CR_
		{
		  if ( !( (pd_glob[0]->Continuation != HUN_ZEROTH) 
		           && (pd_glob[0]->Continuation != HUN_FIRST) 
		           && (pd_glob[0]->Continuation != LOCA) )
		           && accept_hunting_conditions  )		
		  {
                    if ( !accept_hcs_until_end_of_hc_card && (number_of_hcs_expected !=  number_of_hcs_found) )
		    { 
		      sprintf(msg, "Wrong number of Hunting conditions found (%i expected, %i found).", number_of_hcs_expected, number_of_hcs_found ); 
		      declare_error(msg); 
	  	    }	  	  		
                    accept_hcs_until_end_of_hc_card = FALSE; 		  
		  }
		card_read("END_OF_HC_CARD",file_index);  		
		}
		| error {}
		;
		
		
solver_specifications_delimiter_card:
		/* empty */ {}
		| SOLVER_ SPECIFICATIONS_  CR_
		{
		  card_read("SOLVER_SPECIFICATIONS_DELIMITER_CARD",file_index);  
		}
		| error {}
		;







solution_algorithm_card:
		/* empty */ {}
		| SOLUTION_ ALGORITHM_ EQUALS_ solution_algorithm CR_
		{
		  Linear_Solver = $4;
		  card_read("SOLUTION_ALGORITHM_CARD",file_index);  
		}
		| error {}
		;

solution_algorithm:
		  LU_ 	{$$=SPARSE13a;}
		| FRONT_{$$=FRONT;}
		| UMF_ 	{$$=UMFPACK2;}
		| UMFF_ {$$=UMFPACK2F;}
		| MA28_ {$$=MA28;}
	/*	| Y12M_ {}
		| GMRES_{}  These options no longer appear in mm_input and are not implemented here.
		| CG_ 	{}
		| CGS_ 	{}
		| TFQMR_	{}
		| BICGSTAB_ 	{}
	*/	| error {}
		;

matrix_storage_format_card:
		/* empty */ {}
		| MATRIX_ STORAGE_ FORMAT_ EQUALS_ valid_storage_format  CR_
		{
		  if (Linear_Solver == FRONT )
		  {
		    strcpy(Matrix_Format,$5);
		  }
		  card_read("MATRIX_STORAGE_FORMAT_CARD",file_index);  
		}
		| error {}
		;		
		
valid_storage_format:
		  VBR_ {/*strcpy($$,"vbr");*/}
		| MSR_ {/*strcpy($$,"msr");*/}
		| error {}
		;
		
preconditioner_card:
		/* empty */ {}
		| PRECONDITIONER_ EQUALS_ valid_preconditioner CR_
		{
		  strcpy(Matrix_Preconditioner,$3);
		  card_read("PRECONDITIONER_CARD",file_index);  
		}
		| error {}
		;
		
valid_preconditioner:
		  NONE_ 	{}
		| JACOBI_ 	{}
		| NEUMANN_ 	{}
		| LS_ 		{}
		| LU_ 		{}
		| ILU_ 		{}
		| ILUT_		{}
		| BILU_ 	{}
		| SYM_GS_ 	{}
		| DOM_DECOMP_	{}
		| ICC_		{}
		| ILU_		{}
		| RILU_		{}
		| error 	{}
		;

matrix_subdomain_solver_card:
		/* empty */ {}
		| MATRIX_ SUBDOMAIN_ SOLVER_ EQUALS_ valid_subdomain_solver CR_
		{
		  strcpy(Matrix_Subdomain_Solver,$5);
		  card_read("MATRIX_SUBDOMAIN_SOLVER_CARD",file_index);  
		}
		| error {}
		;
				
valid_subdomain_solver:
		  NONE_ {}
		| ILU_ 		{}
		| ILUT_		{}
		| BILU_ 	{}
		| ICC_		{}
		| ILU_		{}
		| RILU_		{}
		| error 	{}
		;
				
matrix_fill_factor_card:
		/* empty */ {}
		| MATRIX_  FILL_ FACTOR_ EQUALS_ float  CR_
		{
		  card_read("MATRIX_FILL_FACTOR_CARD",file_index);  
		}
		| error {}
		;	
		
		
matrix_scaling_card:
		/* empty */ {}
		| MATRIX_ SCALING_ EQUALS_ valid_matrix_scaling CR_
		{
		  strcpy(Matrix_Scaling,$4);
		  card_read("MATRIX_SCALING_CARD",file_index);  
		}
		| error {}
		;		
		
valid_matrix_scaling:
		  NONE_ 	{}
		| JACOBI_ 	{}
		| BJACOBI_ 	{}
		| ROW_SYM_	{}
		| SYM_DIAG_	{}
		| SYM_ROW_SUM_	{}
		| error 	{}
		;		
		
matrix_residual_norm_type_card:
		/* empty */ {}
		| MATRIX_ RESIDUAL_ NORM_ TYPE_ EQUALS_ valid_residual_norm_type CR_
		{
		  strcpy(Matrix_Residual_Norm_Type,$6);
		  card_read("MATRIX_RESIDUAL_NORM_TYPE_CARD",file_index);  
		}
		| error {}
		;
		
valid_residual_norm_type:
		  R0_ 		{}
		| RHS_ 		{}
		| ANORM_ 	{}
		| NOSCALED_	{}
		| SOL_		{}
		| WEIGHTED_	{}
		| error 	{}
		;		
		
matrix_output_type_card:
		/* empty */ {}
		| MATRIX_ OUTPUT_ TYPE_ EQUALS_ valid_matrix_output_type CR_
		{
		  strcpy(Matrix_Output_Type,$5);
		  card_read("MATRIX_OUTPUT_TYPE_CARD",file_index);  
		}
		| error {}
		;
		
valid_matrix_output_type:
		  NONE_ 	{}
		| ALL_	 	{}
		| LAST_ 	{}
		| WARNINGS_	{}
		| error 	{}
		;		
				
matrix_factorization_reuse_card:
		/* empty */ {}
		| MATRIX_ FACTORIZATION_ REUSE_ EQUALS_ valid_matrix_reuse_type CR_
		{ 
		  strcpy(Matrix_Factorization_Reuse,$5);
		  card_read("MATRIX_FACTORIZATION_REUSE_CARD",file_index);  
		}
		| error {}
		;
		
valid_matrix_reuse_type:
		  CALC_
		| RECALC_
		| REUSE_
		| error {}
		;
		
matrix_factorization_overlap_card:
		/* empty */ {}
		| MATRIX_ FACTORIZATION_ OVERLAP_ EQUALS_ valid_matrix_factorization_overlap CR_
		{
		  strcpy(Matrix_Factor_Overlap,$5);
		  card_read("MATRIX_FACTORIZATION_OVERLAP_CARD",file_index);
		}
		| error {}
		;
		
valid_matrix_factorization_overlap:
		  NONE_
		| DIAG_
		| FULL_
		| error {}
		;

matrix_overlap_type_card:
		/* empty */ {}
		| MATRIX_ OVERLAP_ TYPE_ EQUALS_ valid_matrix_overlap_type CR_
		{
		  strcpy(Matrix_Overlap_Type,$5);
		  card_read("MATRIX_OVERLAP_TYPE_CARD",file_index);
		}
		| error {}
		;
		
valid_matrix_overlap_type:
		  STANDARD_
		| SYMMETRIC_
		| error {}
		;
	
	
matrix_reorder_card:
		/* empty */ {}
		| MATRIX_ REORDER_ EQUALS_ valid_matrix_reorder_type CR_
		{
		  strcpy(Matrix_Reorder,$4);
		  card_read("MATRIX_REORDER_CARD",file_index);
		}
		| error {}
		;
		
valid_matrix_reorder_type:
		  NONE_
		| RCM_
		| error {}
		;	

size_of_krylov_subspace_card:
		/* empty */ {}
		| SIZE_ OF_ KRYLOV_ SUBSPACE_ EQUALS_ orthogonalization_option CR_
		{
		  strcpy(Matrix_Krylov_Subspace,$6);
		  card_read("SIZE_OF_KRYLOV_SUBSPACE_CARD",file_index);
		}
		| error {}
		;
					
				
matrix_graph_fillin_card:
		/* empty */ {}
		| MATRIX_ GRAPH_ FILLIN_ EQUALS_ integer CR_
		{
		  sprintf(Matrix_Graph_Fillin,"%i",$5);
		  card_read("MATRIX_GRAPH_FILLIN_CARD",file_index);
		}
		| error {}
		;
		
matrix_factorization_save_card:
		/* empty */ {}
		| MATRIX_ FACTORIZATION_ SAVE_ EQUALS_ integer CR_
		{
		  #ifdef AZTEC_2
		    sprintf(Matrix_Factorization_Save,"%i",$5);
		  #endif
		  card_read("MATRIX_FACTORIZATION_SAVE_CARD",file_index);		  
		}
		| error {}
		;
		
		
matrix_ilut_fill_factor_card:
		/* empty */ {}
		| MATRIX_ ILUT_ FILL_ FACTOR_ EQUALS_ float CR_
		{
		  #ifdef AZTEC_2
		    strcpy(Matrix_ILUT_Fill_Factor,$5);
		  #endif
		  card_read("MATRIX_ILUT_FILL_FACTOR_CARD",file_index);		  
		}
		| error {}
		;				
		
matrix_rilu_relax_factor_card:
		/* empty */ {}
		| MATRIX_ RILU_ RELAX_ FACTOR_ EQUALS_ float CR_
		{
		  #ifdef AZTEC_2
		    sprintf(Matrix_RILU_Relax_Factor,"%f",$5);
		  #endif
		  card_read("MATRIX_RILU_RELAX_FACTOR_CARD",file_index);		  
		}
		| error {}
		;				
				
matrix_bilu_threshold_card:
		/* empty */ {}
		| MATRIX_ BILU_ THRESHOLD_ EQUALS_ float CR_
		{
		  #ifdef AZTEC_2
		    sprintf(Matrix_BILU_Threshold,"%f",$5);
		  #endif
		  card_read("MATRIX_BILU_THRESHOLD_CARD",file_index);
		}
		| error {}
		;
		
matrix_relative_threshold_card:
		/* empty */ {}
		| MATRIX_ RELATIVE_ THRESHOLD_ EQUALS_ float CR_
		{
		  #ifdef TRILINOS
		    sprintf(Matrix_Relative_Threshold,"%f",$5);
		  #endif
		  card_read("MATRIX_RELATIVE_THRESHOLD_CARD",file_index);
		}
		| error {}
		;
		
matrix_absolute_threshold_card:
		/* empty */ {}
		| MATRIX_ ABSOLUTE_ THRESHOLD_ EQUALS_ float CR_
		{
		  #ifdef TRILINOS
		    sprintf(Matrix_Absolute_Threshold,"%f",$5);
		  #endif
		  card_read("MATRIX_ABSOLUTE_THRESHOLD_CARD",file_index);
		}
		| error {}
		;
				
		
matrix_polynomial_order_card:
		/* empty */ {}
		| MATRIX_ POLYNOMIAL_ ORDER_ EQUALS_ integer  CR_
		{
		  sprintf(Matrix_Polynomial_Order,"%i",$5);
		  card_read("MATRIX_POLYNOMIAL_ORDER_CARD",file_index); 
		}
		| POLYNOMIAL_ EQUALS_ integer CR_
		{
		  sprintf(Matrix_Polynomial_Order,"%i",$3); /* ??? don't remember where this form of the card came from.  Is is valid? */		
		  card_read("MATRIX_POLYNOMIAL_ORDER_CARD",file_index); 
		}		
		| error {}
		;
		
size_of_krylov_subspace_card:
		/* empty */ {}
		| SIZE_ OF_ KRYLOV_ SUBSPACE_ EQUALS_ integer  CR_
		{
		  sprintf(Matrix_Krylov_Subspace,"%i",$6);
		  card_read("SIZE_OF_KRYLOV_SUBSPACE_CARD",file_index);  
		}
		| error {}
		;
		
orthogonalization_card:
		/* empty */ {}
		| ORTHOGONALIZATION_ EQUALS_ orthogonalization_option CR_
		{
		  strcpy(Matrix_Orthogonalization,$3);
		  card_read("ORTHOGONALIZATION_CARD",file_index);  
		}
		| error {}
		;
		
orthogonalization_option:
		  CLASSIC_ {}
		| MODIFIED_ {}
		| error {}
		;		
		
maximum_linear_solve_iterations_card:
		/* empty */ {}
		|  MAXIMUM_ LINEAR_ SOLVE_ ITERATIONS_ EQUALS_ integer CR_
		{
		  sprintf(Matrix_Maximum_Iterations,"&i",$6);
		  card_read("MAXIMUM_LINEAR_ITERATIONS_SOLVE_CARD",file_index);
		}
		| error {}
		;
		
matrix_auxiliary_vector_card:
		/* empty */ {}
		| MATRIX_ AUXILIARY_ VECTOR_ EQUALS_ valid_matrix_auxiliary_vector_type CR_
		{
		  strcpy(Matrix_Auxiliary_Vector,$5);
		  card_read("MATRIX_AUXILIARY_VECTOR_CARD",file_index); 
		}
		| error {}
		;
		
valid_matrix_auxiliary_vector_type:
		  RESID_
		| RAND_
		| error {}
		;		
		
		
matrix_drop_tolerance_card:
		/* empty */ {}
		| MATRIX_ DROP_ TOLERANCE_ EQUALS_ float  CR_
		{
		  sprintf(Matrix_Drop_Tolerance,"%f",$5);
		  card_read("MATRIX_DROP_TOLERANCE_CARD",file_index);  
		}
		| error {}
		;
		
number_of_newtonian_iterations_card:
		/* empty */ {}
		| NUMBER_ OF_ NEWTON_ ITERATIONS_ EQUALS_ integer optional_integer CR_
		{
		  Max_Newton_Steps = $6;
		  if( $7 != NA )
		  {
		    Newt_Jacobian_Reformation_stride = $7;
		    if (Newt_Jacobian_Reformation_stride > 1)
		    {
		      modified_newton = TRUE;
		    }		  
		  }
		  else
		  {
		    Newt_Jacobian_Reformation_stride = 1;
		    declare_warning("Defaulting to un_modified Newton iterations.");		  
		  }
		  card_read("NUMBER_OF_NEWTONIAN_ITERATIONS_CARD",file_index);  	
		}
		| error {}
		;
		
newton_correction_factor_card:
		/* empty */ {}
		| NEWTON_ CORRECTION_ FACTOR_ EQUALS_ float  optional_float optional_float optional_float optional_float optional_float CR_
		{
		  damp_factor1 = $5;
		  if (  ($6 != NA) & ($7!= NA) & ($8 != NA) & ($9 != NA) & ($10 != NA) )
		  {
		    /* five more floats were read */
		    custom_tol1  = $6;
		    damp_factor2 = $7;
                    custom_tol2  = $8;
                    damp_factor3 = $9;
                    custom_tol3  = $10;
 		    if( ( damp_factor1 <= 1. && damp_factor1 >= 0.) &&
                       ( damp_factor2 <= 1. && damp_factor2 >= 0.) &&
                       ( damp_factor3 <= 1. && damp_factor3 >= 0.))
                    {
                      /* fprintf(stdout,"Consistent Newton Damping factor scheme found\n"); */
                    }
                    else
                    {
                      declare_error("All damping factors must be in the range 0 <= fact <=1");
                    }
		  }
		  else
		  {
		    /* five more floats were not read */
		    damp_factor2 = damp_factor3 =  -1.;
                    custom_tol1 = custom_tol2 = custom_tol3 = -1;
                    declare_warning("Defaulting to constant Newton relaxation. Use 6 flts for dynamic relaxation.");
		  }
		  card_read("NEWTON_CORRECTION_FACTOR_CARD",file_index);  
		}
		| error {}
		;		
		
normalized_residual_tolerance_card:
		/* empty */ {}
		| NORMALIZED_ RESIDUAL_ TOLERANCE_ EQUALS_ float  CR_
		{
		  Epsilon[0] = $5;
		  card_read("NORMALIZED_RESIDUAL_TOLERANCE_CARD",file_index);  
		}
		| error {}
		;

normalized_corection_tolerance_card:
		/* empty */ {}
		| NORMALIZED_ CORRECTION_ TOLERANCE_ EQUALS_ float  CR_
		{
		  Epsilon[2] = $5;
		  card_read("NORMALIZED_CORRECTION_TOLERANCE_CARD",file_index);  
		}
		| error {}
		;
		
modified_newton_tolerance_card:
		/* empty */ {}
		| MODIFIED_ NEWTON_ TOLERANCE_ EQUALS_ float float CR_
		{
		  convergence_rate_tolerance = $5;
		  modified_newt_norm_tol = $6;
		  card_read("MODIFIED_NEWTON_TOLERANCE_CARD",file_index);  
		}
		| error {}
		;
		
jacobian_reform_time_stride_card:
		/* empty */ {}
		| JACOBIAN_ REFORM_ TIME_ STRIDE_ EQUALS_ integer  CR_
		{
		  Time_Jacobian_Reformation_stride = $6;
		  if ( Time_Jacobian_Reformation_stride > 1)
		  {
		    modified_newton = TRUE;
		  }
		  card_read("JACOBIAN_REFORM_TIME_STRIDE_CARD",file_index);  
		}
		| error {}
		;
		
residual_ratio_tolerance_card:
		/* empty */ {}
		| RESIDUAL_ RATIO_ TOLERANCE_ EQUALS_ float CR_
		{
		  sprintf(Matrix_Convergence_Tolerance,"%f",$5); 
		  Epsilon[1] = $5;
		  card_read("RESIDUAL_RATIO_TOLERANCE_CARD",file_index);
		}
		| error {}
		;
		
pressure_stabilization_card:
		/* empty */ {}
		| PRESSURE_ STABILIZATION_ EQUALS_ yes_no_option CR_
		{
		  if( $4 == 1)
		  {
		    PSPG = 1;
		  }
		  else
		  {
		    PSPG = 0;
		  }
		  card_read("PRESSURE_STABILIZATION_CARD",file_index); 
		}
		| error {}
		;
		
pressure_stabilization_scaling_card:
		/* empty */ {}
		| PRESSURE_ STABILIZATION_ SCALING_ EQUALS_ float  CR_
		{
		  PSPG_scaling = $5;
		  card_read("PRESSURE_STABILIZATION_SCALING_CARD",file_index);  
		}
		| error {}
		;
		
linear_stability_card:
		/* empty */ {}
		| LINEAR_ STABILITY_ EQUALS_ valid_linear_stability_option CR_
		{
		  switch($4)
		  {
                    case 1:              /*YES*/
                    { Linear_Stability = LSA_NORMAL; break;}
                    
                    case 2:              /*NO*/
                    { Linear_Stability = LSA_NONE; break;}
                    
                    case 3:              /*INLINE SAME AS YES*/
                    { Linear_Stability = LSA_NORMAL; break;}                                                            
                    
                    case 4:              /*FILE*/
                    { Linear_Stability = LSA_SAVE; break;}
                    
                    case 5:              /*3D*/
                    { Linear_Stability = LSA_NONE; break;}
                    
                    case 6:              /*3DFILE*/
                    { Linear_Stability = LSA_3D_OF_2D_SAVE; break;}                    		  
                    
		  } /* end switch */
		  card_read("LINEAR_STABILITY_CARD",file_index); 
		}
		| error {}
		;
		
valid_linear_stability_option:
		  YES_ 		{$$=1;}
		| NO_		{$$=2;}
		| INLINE_ 	{$$=3;}
		| FILE_		{$$=4;}
		| THREE_D_	{$$=5;}
		| THREE_D_FILE_ {$$=6;}
		| error 	{}
		;		
		
		
filter_concentration_card:
		/* empty */ {}
		| FILTER_ CONCENTRATION_ EQUALS_ integer float float CR_
		{
		  filter_species_material_number = $4;
		  c_min = $5;
		  c_max = $6;
		  Filter_Species = TRUE;
		  card_read("FILTER_CONCENTRATION_CARD",file_index);
		}
		| error {}
		;
				
solver_pressure_datum_card:
		/* empty */ {}
		| PRESSURE_ DATUM_ EQUALS_ integer float CR_
		{
		  card_read("SOLVER_PRESSURE_DATUM_CARD",file_index);  
		  pressure_datum_element = $4;
		  pressure_datum_value = $5;
		}
		| error {}
		;
		
umf_idim_card:
		/* empty */ {}
		| UMF_ IDIM_ EQUALS_ integer CR_
		{
		  UMFPACK_IDIM = $4;
		  card_read("UMF_IDIM_CARD",file_index);
		}
		| error {}
		;
		
umf_xdim_card:
		/* empty */ {}
		| UMF_ XDIM_ EQUALS_ integer CR_
		{
		  UMFPACK_XDIM = $4;
		  card_read("UMF_XDIM_CARD",file_index); 
		}
		| error {}
		;
		
disable_viscosity_sensitivities_card:
		/* empty */ {}
		| DISABLE_ VISCOSITY_ SENSITIVITIES_ EQUALS_ yes_no_option CR_
		{
		  Include_Visc_Sens = $5;
		  card_read("FILL_SUBCYCLE_CARD",file_index);  
		}
		| error {}
		;
		
eigensolver_specifications_delimiter_card: 
		/* empty */ {}
		| EIGENSOLVER_ SPECIFICATIONS_ CR_
		{		  
		  card_read("EIGENSOLVER_SPECIFICATIONS_DELIMITER_CARD",file_index); 
		}
		| error {}
		;
		
eigen_number_of_modes_card:  
		/* empty */ {}
		| EIGEN_ NUMBER_ OF_ MODES_ EQUALS_ integer CR_
		{
		  eigen->Eigen_NEV_WANT = $6;
		  card_read("EIGEN_NUMBER_OF_MODES_CARD",file_index); 
		}
		| error {}
		;
		
eigen_record_modes_card:  
		/* empty */ {}
		| EIGEN_ RECORD_ MODES_ EQUALS_ integer CR_
		{
		  eigen->Eigen_Record_Modes = $5;  
		  card_read("EIGEN_RECORD_MODES_CARD",file_index); 
		}
		| error {}
		;
		
eigen_size_of_krylov_subspace_card:  
		/* empty */ {}
		| EIGEN_ SIZE_ OF_ KRYLOV_ SUBSPACE_ EQUALS_ integer CR_
		{
		  eigen->Eigen_Krylov_Subspace = $7;  
		  card_read("EIGEN_SIZE_OF_KRYLOV_SUBSPACE_CARD",file_index); 
		}
		| error {}
		;
		
eigen_maximum_iterations_card:  
		/* empty */ {}
		| EIGEN_ MAXIMUM_ ITERATIONS_ EQUALS_ integer CR_
		{
		  eigen->Eigen_Maximum_Iterations = $5;  
		  card_read("EIGEN_MAXIMUM_ITERATIONS_CARD",file_index); 
		}
		| error {}
		;
		
eigen_number_of_filter_steps_card:  
		/* empty */ {}
		| EIGEN_ NUMBER_ OF_ FILTER_ STEPS_ EQUALS_ integer CR_
		{
		  eigen->Eigen_Filter = $7;  
		  card_read("EIGEN_NUMBER_OF_FILTER_STEPS_CARD",file_index); 
		}
		| error {}
		;
		
eigen_recycle_card:  
		/* empty */ {}
		| EIGEN_ RECYCLE_ EQUALS_ yes_no_option CR_
		{
		  eigen->Eigen_Recycle = $4;  
		  card_read("EIGEN_RECYCLE_CARD",file_index); 
		}
		| error {}
		;
		
eigen_tolerance_card:  
		/* empty */ {}
		| EIGEN_ TOLERANCE_ EQUALS_ float CR_
		{
		  eigen->Eigen_Tolerance = $4;  
		  card_read("EIGEN_TOLERANCE_CARD",file_index); 
		}
		| error {}
		;
		
eigen_initial_vector_weight_card:  
		/* empty */ {}
		| EIGEN_ INITIAL_ VECTOR_ WEIGHT_ EQUALS_ float CR_
		{
		  eigen->Eigen_IV_Wt = $6;
		  card_read("EIGEN_INITIAL_VECTOR_WEIGHT_CARD",file_index); 
		}
		| error {}
		;
		
eigen_initial_shifts_card:  
		/* empty */ {}
		| EIGEN_ INITIAL_ SHIFTS_ EQUALS_ float float float float CR_
		{
		  eigen->Eigen_Shifts[0] = $5; 
		  eigen->Eigen_Shifts[1] = $6; 
		  eigen->Eigen_Shifts[2] = $7; 
		  eigen->Eigen_Shifts[3] = $8; 		  		  		      
		  card_read("EIGEN_INITIAL_SHIFTS_CARD",file_index); 
		}
		| error {}
		;
		
eigen_wave_numbers_card: 		
		/* empty */ {}
		| EIGEN_ WAVE_ NUMBERS_ EQUALS_  floating_point_constant_list CR_
		{
		    if ( floating_point_constant_list_index != -1 )
		    {
		      LSA_wave_numbers = (dbl *)malloc(number_of_constants * sizeof(dbl)); 
		      LSA_number_wave_numbers = read_floats( &LSA_wave_numbers, 0 );
		    }
		  card_read("EIGEN_WAVE_NUMBER_CARD",file_index); 
		}
		| error {}
		;		
		
		
augmenting_conditions_specifications_section_delimiter_card:
		/* empty */ {}
		| AUGMENTING_ CONDITIONS_ SPECIFICATIONS_ CR_ 
		{
		  card_read("AUGMENTING_CONDITIONS_SPECIFICATIONS_SECTION_DELIMITER_CARD",file_index);  
		}
		| error {}
		;				
				
				
augmenting_conditions_initial_guess_card:
		/* empty */ {}
		| AUGMENTING_ CONDITIONS_ INITIAL_ GUESS_ EQUALS_ read_option CR_ 
		{
		  /* if (cards[AUGMENTING_CONDITION_CARD]->times_read[0] > 0)
		  {
		    AC_rd_file = TRUE;
		  }
		  else
		  {
		    sprintf(msg, "Augmenting conditions initial guess card must occur before the first AC = card."); 
		    declare_error(msg); 		    
		  }*/
		  card_read("AUGMENTING_CONDITIONS_INITIAL_GUESS_CARD",file_index);  
		}
		| error {}
		;

read_option:
		  STRING_ {$$=FALSE;}
		| READ_   {$$=TRUE;}
		;	

number_of_augmenting_conditions_card:
		/* empty */ {}
		| NUMBER_ OF_ AUGMENTING_ CONDITIONS_ EQUALS_ integer CR_ 
		{
		  if( $6 >0  )
		  {
		    number_of_acs_expected = $6;
		    accept_acs_until_end_of_ac_card = FALSE;
		  } 
		  else if ($6 == -1)
		  {
		    accept_acs_until_end_of_ac_card = TRUE;		    
		  }
		  else
		  {
		    declare_error("Invalid Number of Augmenting Conditions");
		    number_of_acs_expected = 0;
		    accept_acs_until_end_of_ac_card = FALSE;
		  }
		card_read("NUMBER_OF_AUGMENTING_CONDITIONS_CARD",file_index);	
  		}
		| error {}
		;	

end_of_augmenting_conditions_card:
		/* empty */ {}
		| END_ OF_ AC_ CR_ 
		{
	  	  if ( !accept_acs_until_end_of_ac_card && (number_of_acs_expected !=  number_of_acs_found) )
		  { 
		    sprintf(msg, "Wrong number of augmenting conditions found (%i expected, %i found).", number_of_acs_expected, number_of_acs_found ); 
		    declare_error(msg); 
	  	  }	  	  		
                  accept_acs_until_end_of_ac_card = FALSE; 
                  card_read("END_OF_AUGMENTING_CONDITIONS_CARD",file_index);
		}
		| error {}
		;
		
augmenting_condition_card:
		/* empty */ {}
		/* 1  2       3                              4        5        6               7               8                                     9 */
		| AC_ EQUALS_ ac_userbc_augmenting_condition integer integer optional_string optional_string optional_floating_point_constant_list CR_ 
		{
		  common_AC_processing();
		  augc[iAC].Type = $3;
		  if ( floating_point_constant_list_index != -1 )
		  {
		    augc[iAC].len_AC = read_floats( &(augc[iAC].DataFlt), 0 );
		  }
		  if ( augc[iAC].BCID == APREPRO_AC_BCID )
		  {
		    if ( strcmp($6,"NA" ) )
		    {
		      strcpy(augc[iAC].Params_File, $6);
		    }
		    else
		    {
		      sprintf(msg, "Error reading Parameter File name."); 
		      declare_error(msg); 
		    }
		    if ( strcmp($7,"NA" ) )
		    {
		      strcpy(augc[iAC].AP_param, $7);		    
		    }
		    else
		    {
		      sprintf(msg, "Error reading Parameter name."); 
		      declare_error(msg); 		    
		    }	
		    if( Num_Proc != 1 )
  		    {
		      sprintf(msg, "Aprepro parameter AC not ready for parallel."); 
		      declare_error(msg);   		    
  		    } 	  
		  }		  
		  card_read("AUGMENTING_CONDITION_CARD",file_index);  
		}		
		| AC_ EQUALS_ ac_usermat_augmenging_condition integer integer optional_floating_point_constant_list CR_ 
		{
		  if ( common_AC_processing() )
		  {
		    augc[iAC].Type = $3;
		    if ( floating_point_constant_list_index != -1 )
		    {
		      augc[iAC].len_AC = read_floats( &(augc[iAC].DataFlt), 0 );
		    }
		  }		  		
		  card_read("AUGMENTING_CONDITION_CARD",file_index);  
		}
		| AC_ EQUALS_ ac_volume_augmenting_condition integer integer integer integer integer float optional_floating_point_constant_list CR_ 
		{
		  if ( common_AC_processing() )
		  {
		    augc[iAC].Type = $3;
		    if ( floating_point_constant_list_index != -1 )
		    {
		      augc[iAC].len_AC = read_floats( &(augc[iAC].DataFlt), 0 );
		    }
		  }		  
		  card_read("AUGMENTING_CONDITION_CARD",file_index);  
		}
		| AC_ EQUALS_ ac_flux_augmenting_condition integer integer integer mfid optional_string optional_string optional_floating_point_constant_list CR_
		/* 1  2       3                            4        5        6        7                    8               9                  */
		{
		  if ( common_AC_processing() )
		  {
		    augc[iAC].Type = $3;		
		    if ( floating_point_constant_list_index != -1 )
		    {
		      augc[iAC].len_AC = read_floats( &(augc[iAC].DataFlt), 0 );
		    }		  
		    if ( augc[iAC].BCID == APREPRO_AC_BCID )
		    {
		      if ( strcmp($8,"NA" ) )
		      {
		        strcpy(augc[iAC].Params_File, $8);
		      }
		      else
		      {
		        sprintf(msg, "Error reading Parameter File name."); 
		        declare_error(msg); 
		      }
		      if ( strcmp($9,"NA" ) )  /* jjj check this card */
		      {
		        strcpy(augc[iAC].AP_param, $9);		    
		      }
		      else
		      {
		        sprintf(msg, "Error reading Parameter name."); 
		        declare_error(msg); 		    
		      }	
		      if( Num_Proc != 1 )
  		      {
		        sprintf(msg, "Aprepro parameter AC not ready for parallel."); 
		        declare_error(msg);   		    
  		      } 	  
		    }
		  }
		  card_read("AUGMENTING_CONDITION_CARD",file_index);  
		}						
		| error {}
		;			
		
mfid:
		  FORCE_X_	{$$=FORCE_X;}
		| FORCE_Y_	{$$=FORCE_Y;}
		| FORCE_Z_	{$$=FORCE_Z;}
		| FORCE_NORMAL_	{$$=FORCE_NORMAL;}
		| FORCE_TANGENT1_	{$$=FORCE_TANGENT1;}
		| FORCE_TANGENT2_	{$$=FORCE_TANGENT2;}
		| VOLUME_FLUX_		{$$=VOLUME_FLUX;}
		| HEAT_FLUX_	{$$=HEAT_FLUX;}
		| SPECIES_FLUX_	{$$=SPECIES_FLUX;}
		| error {}
		;
						
ac_userbc_augmenting_condition:
		  BC_ 		{$$=AC_USERBC;}
		| USERBC_	{$$=AC_USERBC;}
		| error {}
		;
		
ac_usermat_augmenging_condition:
		  MT_ 		{$$=AC_USERMAT;}
		| USERMAT_	{$$=AC_USERMAT;}
		| error {}
		;

ac_volume_augmenting_condition:
		  VC_ 		{$$=AC_VOLUME;}
		| VOLUME_	{$$=AC_VOLUME;}
		| error {}
		;

ac_flux_augmenting_condition:
		  FC_ 		{$$=AC_FLUX;}
		| FLUX_		{$$=AC_FLUX;}
		| error {}
		;
		
boundary_condition_specifications_delimiter_card:
		/* empty */ {}
		| BOUNDARY_ CONDITION_ SPECIFICATIONS_ CR_ 
		{
		  card_read("BOUNDARY_CONDITION_SPECIFICATIONS_DELIMITER_CARD",file_index);  
		}
		| error {}
		;

number_of_boundary_conditions_card:
		/* empty */ {}
		| NUMBER_ OF_ BC_ EQUALS_ integer CR_ 
		{
		  if( $5 >0 )
		  {
		    number_of_bcs_expected = $5;
		    accept_bcs_until_end_of_bc_card = FALSE;
		  } 
		  else if ($5 == -1)
		  {
		    accept_bcs_until_end_of_bc_card = TRUE;
		  }
		  else
		  {
		    declare_error("Invalid Number of Boundary Conditions");
		    number_of_bcs_expected = 0;
		    accept_bcs_until_end_of_bc_card = FALSE;
		  }
		card_read("NUMBER_OF_BOUNDARY_CONDITIONS_CARD",file_index);
  		}
		| error {}
		;
				
end_of_BC_card:
		/* empty */ {}
		| END_ OF_ BC_ CR_ 
		{  
	  	  if ( !accept_bcs_until_end_of_bc_card && (number_of_bcs_expected !=  number_of_bcs_found) )
		  { 
		    sprintf(msg, "Wrong number of boundary conditions found (%i expected, %i found).", number_of_bcs_expected, number_of_bcs_found ); 
		    declare_error(msg); 
	  	  }	  	  		
                  accept_bcs_until_end_of_bc_card = FALSE; 
                  card_read("END_OF_BC_CARD",file_index);
		}
		| error {}
		;
			
bc_card:
	/*empty */ {}
/********************************************************************************************************/	
/*													*/
/* These are the definitions of the various BC formats types:						*/
/*													*/
/* BC Format Number	Data Reqired on card in additon to the BC name, BC type & BC ID			*/
/*													*/
/* 1	No additonal data required									*/
/*      (HYDROSTATIC_SYMM, PSPG, Q_VELO_SLIP)								*/
/*													*/
/* 2	<float> [optional integer]									*/
/*	(QSIDE, TNRMLSIDE, TSHRSIDE, CURRENT,								*/
/*	 KINEMATIC, KINEMATIC_PETROV, KINEMATIC_COLLOC, KINEMATIC_DISC, T_MELT,				*/
/*	 VELO_NORMAL, VELO_NORM_COLLOC, VELO_NORMAL_DISC, SURFTANG_SCALAR, FLOW_STRESSNOB,		*/
/*       FLOW_GRADV, FILL_INLET, FILL_CA, FLOW_PRESSURE, STRONG_FILL_CA)				*/
/*													*/
/* 3	<second bc_id> <float> [otional integer]							*/
/*	(KINEMATIC_EDGE, SURFTANG_SCALAR_EDGE)								*/
/*													*/
/* 4	<float> [optional_float] [otional integer]							*/
/*	(DISTNG, DX, DX_RS, DXDISTNG,									*/
/*	 DY, DY_RS, DYDISTNG, DZ, DZ_RS, DZDISTNG,							*/
/*	 P,T, U, V, W, PU, PV, PW,									*/
/*	 S11, S12, S13, S22, S23, S33, 									*/
/*	 S11_1, S12_1, S13_1, S22_1, S23_1, S33_1,							*/
/*	 S11_2, S12_2, S13_2, S22_2, S23_2, S33_2,							*/
/*	 S11_3, S12_3, S13_3, S22_3, S23_3, S33_3,							*/
/*	 S11_4, S12_4, S13_4, S22_4, S23_4, S33_4,							*/
/*	 S11_5, S12_5, S13_5, S22_5, S23_5, S33_5,							*/
/*	 S11_6, S12_6, S13_6, S22_6, S23_6, S33_6,							*/
/*	 S11_7, S12_7, S13_7, S22_7, S23_7, S33_7,							*/
/*	 G11, G12, G13, G21, G22, G23, G31, G32, G33, 							*/
/*	 VOLT, F, SH, CONT_TANG_VEL, CONT_NORM_VEL)							*/
/*													*/
/* 5	<second bc_id> <float> [optional float] [otional integer]					*/
/*	(VELO_NORMAL_EDGE, VELO_NORMAL_EDGE_INT,							*/
/*	 CA_EDGE_CURVE, CA_EDGE_CURVE_INT)								*/
/*                                                                                                      */
/* 6	<float> <float>											*/
/*      (QCONV, QRAD, KIN_LEAK, VNORM_LEAK)								*/
/*													*/
/* 7	<float> <float> <float>										*/
/*	(FORCE, FORCE_RS, NORM_FORCE, NORM_FORCE_RS, 							*/
/*	 SLOPEX, SLOPEY, SLOPEZ, SLOPE, CAPILLARY)							*/
/*													*/
/* 8	<second bc_id> <float> <float> <float>								*/
/*	(VELO_TANGENT_EDGE, VELO_TANGENT_EDGE_INT)							*/
/*													*/
/* 9	<float> <float> <float> <float>									*/
/*	(PLANEX, PLANEY, PLANEZ, PLANE, 								*/
/*	 CA, SURFTANG, FLOW_HYDROSTATIC,								*/
/*	 VELO_TANGENT_3D, VELO_SLIP, VELO_SLIP_ROT, VELO_SLIP_FILL)					*/
/*													*/
/* 10	<second bc_id> <float> <float> <float> <float>							*/	
/*	(SURFTANG_EDGE, CA_EDGE, CA_EDGE_INT)								*/
/*                                                                                                      */
/* 11	<float> <float> <float> <float> <float>								*/
/*	(five floats)											*/
/*	(REP_FORCE, REP_FORCE_RS)									*/
/*                                                                                                      */
/* 12	<float> <float> <float> <float> <float>	<float> <float>						*/
/*	(seven floats)											*/
/*      (CAP_REPULSE, CA_OR_FIX, MOVING_PLANE)                                                          */
/*													*/
/* 13	<second bc_id> <float> <float> <float> <float> <float> <float> <float> <float>			*/	
/*	(second bc_id + eight floats)									*/
/*	(VAR_CA_EDGE)											*/
/*													*/
/* 14	<float> <float> <float> <float> <float> <float> <float> <float> <float> <float> 		*/
/*	(ten floats)											*/
/*	(MOVING_CA)											*/
/*													*/
/* 15	<integer>											*/
/*	(POROUS_CONV, YFLUX_SUS, KIN_DISPLACEMENT)							*/
/*													*/
/* 16	<integer> <integer> <integer> <float> <float> <float> <float> <float>				*/	
/*	(VL_EQUIL)											*/
/*													*/
/* 17	<keyword> <integer> <float> <float> <float>							*/
/*	(YFLUX_EQUIL)											*/
/*													*/
/* 18	<keyword> <integer> <integer> <integer> <float>		 					*/
/*	(VL_POLY)											*/
/*													*/
/* 19	<float> <float> <float> <float> <float>..............						*/
/*	(arbitrary number of floats)									*/
/*      (UVARY, VVARY, WVARY, SPLINEX, SPLINEY, SPLINEZ,						*/
/*	 SPLINE, T_USER, VOLT_USER, UUSER, VUSER, WUSER, QUSER,							*/
/*	 FORCE_USER, FORCE_USER_RS, PRESSURE_USER, FLOW_PRESS_USER, KIN_CHEM)                           */
/*													*/
/* 20	<integer> <float> <float> <float> <float> <float>.............					*/
/*	(integer followed by an arbitrary number of floats)						*/
/*      (YFLUX_USER, YUSER, Y2_ELECTRONEUTRALITY)                                                       */
/*													*/
/* 21	<second bc_id> <keyword> <float> <float> <float> <float> <float>.............			*/
/*	(second bc_id, an integer and an arbitrary number of floats)					*/
/*      (CA_EDGE_OR_FIX)										*/
/*                                                                                                      */
/* 22	<integer> <float> [optional float]								*/	
/*	(Y, Y_DISCONTINUOUS, YFLUX_CONST, YTOTALFLUX_CONST, POROUS_KIN, KINEMATIC_SPECIES)		*/
/*													*/
/* 23	<integer> <float> <float>									*/
/*	(YFLUX, POROUS_FLUX)										*/
/*													*/
/* 24	<integer> <float> <float> <float> <float> <float> <float> <float> <float>			*/
/*	(integer followed by eight floats)								*/
/*      (YFLUX_BV)                                                                                      */
/*													*/
/* 25	<integer> <integer> [optional float]								*/	
/*	(SOLID_FLUID, FLUID_SOLID,									*/
/*	 SOLID_FLUID_RS, FLUID_SOLID_RS,								*/
/*	 PENETRATION, POROUS_PRESSURE, DARCY_CONTINUOUS,						*/
/*	 NO_SLIP, NO_SLIP_RS, VELO_TANGENT_SOLID)							*/
/*													*/
/* 26	<integer> <integer> <float> [optional integer] [optional float]					*/
/*	(VELO_SLIP_SOLID)										*/
/*													*/
/* 27	<integer> <integer> <integer>									*/
/*	(P_EQUIL)											*/
/*													*/
/* 28	<integer> <integer> <integer> <float>								*/	
/*	(VN_POROUS, VP_EQUIL)										*/
/*													*/
/* 29	<integer> <integer> <integer> <float> <float>							*/
/*	(POROUS_GAS)											*/
/*													*/
/* 30	<integer> <float> <float> <float>								*/
/*	(VELO_TANGENT, LATENT_HEAT)									*/
/*													*/
/* 31	<ordinate> <abcissa> <interpolation method> [optional FILE = string] [optional FILE = string]  	*/
/*	(TABLE BC)											*/
/*      (TABLE_WICV, TABLE_WICS, TABLE)		        						*/
/*													*/
/* 32	<second bc_id>                                                                                  */
/*      (bc requries a second ssid)									*/
/*	(COLLOCATE_EDGE, WEAK_INT_EDGE, STRONG_INT_EDGE)						*/
/*													*/
/* 33	<var_name> <integer>										*/
/*      (FIX bc)											*/
/*                                                                                                      */
/* 34   <equation_name> <integer> <var_name> <integer> <float>						*/	
/*	(GD_CONST_ bc)											*/
/*                                                                                                      */
/* 35   <equation_name> <integer> <var_name> <integer> <float> <float>					*/
/*	(GD_LINEAR_ bc)											*/
/*                                                                                                      */
/* 36	<equation_name> <integer> <var_name> <integer> <float> <float> <float>				*/
/*	(GD_PARAB_  bc)											*/
/*                                                                                                      */
/* 37	<equation_name> <integer> <var_name> <integer> <float> <float> <float>				*/	
/*	(GD_CIRC_ bc)											*/
/*                                                                                                      */
/* 38	<equation_name> <integer> <var_name> <integer> 							*/
/*      <float> <float> <float> <float> <float> <float> <float>						*/
/*	(GD_POLYN_ bc, seven floats)									*/
/*                                                                                                      */
/* 39	<equation_name> <integer> <gd_time_table_option> 						*/
/*      <integer>											*/
/*      [optional_float]										*/
/*      [optional_float]										*/
/*      [optional_LINEAR]										*/
/*      [optional_file_equals_string]									*/
/*      [optional_name_equals_string]									*/
/*      (GD_TIME_ bc)											*/
/*													*/
/* 40	<equation_name> <integer> <var_name> <integer> 							*/
/*      <float>												*/
/*      LINEAR 												*/
/*      [optional_file_equals_string]									*/
/*      [optional_name_equals_string]									*/
/*      (GD_TABLE_ bc)											*/	
/*													*/
/* 41	<discontinuous_velo_bc_keyword> <integer> <integer>						*/
/*	(DISCONTINUOUS_VELO bc)										*/
/*													*/	
/* 42	<latent_heat_internal_bc_keyword> <integer> <integer> <float>					*/
/*	(LATENT_HEAT_INTERNAL bc)									*/
/*													*/	
/*													*/
/********************************************************************************************************/		
	| BC_ EQUALS_ bc_format_1 bc_type bc_id CR_ 
  	  {
  	    card_read("BC_CARD",file_index); /* BC Format 1 */
 	    common_BC_processing($3,$4,$5);
    	  }
	| BC_ EQUALS_ bc_format_2 bc_type bc_id float optional_integer CR_ 
  	  {
    	    card_read("BC_CARD",file_index); /* BC Format 2 */
  	    if ( common_BC_processing($3,$4,$5) )
     	    {
    	      BC_Types[iBC].BC_Data_Float[0] = $6;
    	      if( NA == $7 )
     	        BC_Types[iBC].BC_Data_Int[0] = -1;
    	      else
       	        BC_Types[iBC].BC_Data_Int[0] = $7;
       	     }
    	  } 
	| BC_ EQUALS_ bc_format_3 bc_type bc_id bc_id float optional_integer CR_ 
  	  { 
  	    card_read("BC_CARD",file_index);  /* BC Format 3 */
     	    if ( common_BC_processing($3,$4,$5) )
     	    {
     	      BC_Types[iBC].BC_ID2 = $6;
    	      BC_Types[iBC].BC_Data_Float[0] = $7;
    	      if( NA == $8 )
      	        BC_Types[iBC].BC_Data_Int[0] = -1;
    	      else
  	        BC_Types[iBC].BC_Data_Int[0] = $8;
  	    }
    	  }    	  
    	  
    	| BC_ EQUALS_ bc_format_4 bc_type bc_id float optional_float optional_integer CR_ 
  	  {
  	    card_read("BC_CARD",file_index);  /* BC Format 4 */
            if ( common_BC_processing($3,$4,$5) )
    	    {
    	      BC_Types[iBC].BC_Data_Float[0] = $6;
     	      if( $7 == NA)
      	        BC_Types[iBC].BC_Data_Int[0] = -1.0;   
    	      else
      	        BC_Types[iBC].BC_relax = $7;	  
    	      if( NA == $8 )
    	        BC_Types[iBC].BC_EBID_Apply = -1;
    	      else
    	        BC_Types[iBC].BC_EBID_Apply = $8;
    	    }
      	  }
    	| BC_ EQUALS_ bc_format_5 bc_type bc_id bc_id float optional_float optional_integer CR_ 
  	  {
  	    card_read("BC_CARD",file_index);  /* BC Format 5 */
  	    if ( common_BC_processing($3,$4,$5) )
    	    {
              BC_Types[iBC].BC_ID2 = $6;
              BC_Types[iBC].BC_Data_Float[0] = $7;
              if($8 == NA)
                BC_Types[iBC].BC_Data_Int[0] = -1.0;
    	      else
    	        BC_Types[iBC].BC_relax = $8;    	  
    	      if( NA == $9 )
    	        BC_Types[iBC].BC_EBID_Apply = -1; 
    	      else
    	        BC_Types[iBC].BC_EBID_Apply = $9;
    	    }
      	  }     	  
      	| BC_ EQUALS_ bc_format_6 bc_type bc_id float float CR_ 
  	  {
  	    card_read("BC_CARD",file_index);   	  /* BC Format 6 */
    	    if ( common_BC_processing($3,$4,$5) )
    	    {
    	      BC_Types[iBC].BC_Data_Float[0] = $6;
	      BC_Types[iBC].BC_Data_Float[1] = $7;
  	    }
  	  }
	| BC_ EQUALS_ bc_format_7 bc_type bc_id float float float CR_ 
  	  {
  	    card_read("BC_CARD",file_index);  /* BC Format 7 */
      	    if ( common_BC_processing($3,$4,$5) )
    	    {
     	      BC_Types[iBC].BC_Data_Float[0] = $6;
    	      BC_Types[iBC].BC_Data_Float[1] = $7;
     	      BC_Types[iBC].BC_Data_Float[2] = $8;
  	    }
  	  }
	| BC_ EQUALS_ bc_format_8 bc_type bc_id bc_id float float float CR_ 
  	  {
  	    card_read("BC_CARD",file_index);  /* BC Format 8 */
  	    if ( common_BC_processing($3,$4,$5) )
    	    {
    	      BC_Types[iBC].BC_ID2 = $6;
     	      BC_Types[iBC].BC_Data_Float[0] = $7;
    	      BC_Types[iBC].BC_Data_Float[1] = $8;
     	      BC_Types[iBC].BC_Data_Float[2] = $9;
     	    }
  	  }  	  
	| BC_ EQUALS_ bc_format_9 bc_type bc_id float float float float CR_ 
  	  {
  	    card_read("BC_CARD",file_index);  /* BC Format 9 */
    	    if ( common_BC_processing($3,$4,$5) )
    	    {
    	      BC_Types[iBC].BC_Data_Float[0] = $6;
    	      BC_Types[iBC].BC_Data_Float[1] = $7;
     	      BC_Types[iBC].BC_Data_Float[2] = $8;
     	      BC_Types[iBC].BC_Data_Float[3] = $9; 	  
     	    }
  	  }
  	| BC_ EQUALS_ bc_format_10 bc_type bc_id bc_id float float float float CR_ 
  	  {
  	    card_read("BC_CARD",file_index); /* BC Format 10 */
    	    if ( common_BC_processing($3,$4,$5) )
    	    {
     	      BC_Types[iBC].BC_ID2 = $6;   	    
    	      BC_Types[iBC].BC_Data_Float[0] = $7;
    	      BC_Types[iBC].BC_Data_Float[1] = $8;
     	      BC_Types[iBC].BC_Data_Float[2] = $9;
     	      BC_Types[iBC].BC_Data_Float[3] = $10; 	  
     	    }
  	  }
	| BC_ EQUALS_ bc_format_11 bc_type bc_id float float float float float CR_ 
  	  {
  	    card_read("BC_CARD",file_index); /* BC Format 11 */
    	    if ( common_BC_processing($3,$4,$5) )
    	    {
    	      BC_Types[iBC].BC_Data_Float[0] = $6;
    	      BC_Types[iBC].BC_Data_Float[1] = $7;
     	      BC_Types[iBC].BC_Data_Float[2] = $8;
     	      BC_Types[iBC].BC_Data_Float[3] = $9; 
     	      BC_Types[iBC].BC_Data_Float[4] = $10; 	  
     	    }
  	  }  	
	| BC_ EQUALS_ bc_format_12 bc_type bc_id float float float float float float float CR_ 
  	  {
  	    card_read("BC_CARD",file_index); /* BC Format 12 */
    	    if ( common_BC_processing($3,$4,$5) )
    	    {
    	      BC_Types[iBC].BC_Data_Float[0] = $6;
    	      BC_Types[iBC].BC_Data_Float[1] = $7;
     	      BC_Types[iBC].BC_Data_Float[2] = $8;
     	      BC_Types[iBC].BC_Data_Float[3] = $9; 
     	      BC_Types[iBC].BC_Data_Float[4] = $10;     	  
     	      BC_Types[iBC].BC_Data_Float[5] = $11; 
     	      BC_Types[iBC].BC_Data_Float[6] = $12;         	  
     	    }
   	  }	   
	| BC_ EQUALS_ bc_format_13 bc_type bc_id bc_id float float float float float float float float CR_ 
  	  {
  	    card_read("BC_CARD",file_index); /* BC Format 13 */
      	    if (common_BC_processing($3,$4,$5) )
    	    {
     	      BC_Types[iBC].BC_ID2 = $6;    	    
   	      BC_Types[iBC].BC_Data_Float[0] = $7;
    	      BC_Types[iBC].BC_Data_Float[1] = $8;
     	      BC_Types[iBC].BC_Data_Float[2] = $9;
     	      BC_Types[iBC].BC_Data_Float[3] = $10; 
     	      BC_Types[iBC].BC_Data_Float[4] = $11;     	  
     	      BC_Types[iBC].BC_Data_Float[5] = $12; 
     	      BC_Types[iBC].BC_Data_Float[6] = $13;       	  
     	      BC_Types[iBC].BC_Data_Float[7] = $14;    	  
     	    }
   	  }	   
	| BC_ EQUALS_ bc_format_14 bc_type bc_id float float float float float float float float float float CR_ 
  	  {
  	    card_read("BC_CARD",file_index); /* BC Format 14 */
      	    if ( common_BC_processing($3,$4,$5) )
    	    {
   	      BC_Types[iBC].BC_Data_Float[0] = $6;
    	      BC_Types[iBC].BC_Data_Float[1] = $7;
     	      BC_Types[iBC].BC_Data_Float[2] = $8;
     	      BC_Types[iBC].BC_Data_Float[3] = $9; 
     	      BC_Types[iBC].BC_Data_Float[4] = $10;     	  
     	      BC_Types[iBC].BC_Data_Float[5] = $11; 
     	      BC_Types[iBC].BC_Data_Float[6] = $12;       	  
     	      BC_Types[iBC].BC_Data_Float[7] = $13;       	  
    	      BC_Types[iBC].BC_Data_Float[8] = $14;       	  
     	      BC_Types[iBC].BC_Data_Float[9] = $15;       	  
     	    }
     	  }  	
	| BC_ EQUALS_ bc_format_15 bc_type bc_id integer CR_ 
    	  {
    	    card_read("BC_CARD",file_index); /* BC Format 15 */
    	    if (common_BC_processing($3,$4,$5) )
     	    {
     	      BC_Types[iBC].species_eq = BC_Types[iBC].BC_Data_Int[0] = $6;
     	    }
  	  }
  	| BC_ EQUALS_ bc_format_16 bc_type bc_id integer integer integer
  	  float float float float float CR_ 
  	  {
  	    card_read("BC_CARD",file_index); /* BC Format 16 */
     	    if ( common_BC_processing($3,$4,$5) )
    	    {
   	      BC_Types[iBC].BC_Data_Int[0] = $6;
    	      BC_Types[iBC].BC_Data_Int[1] = $7;
     	      BC_Types[iBC].BC_Data_Int[2] = $8;
     	      BC_Types[iBC].BC_Data_Float[0] = $9; 
     	      BC_Types[iBC].BC_Data_Float[1] = $10;     	  
     	      BC_Types[iBC].BC_Data_Float[2] = $11; 
     	      BC_Types[iBC].BC_Data_Float[3] = $12;       	  
     	      BC_Types[iBC].BC_Data_Float[4] = $13;       	  
     	      BC_Types[iBC].species_eq = BC_Types[iBC].BC_Data_Int[0];
     	    }
   	  }	
	| BC_ EQUALS_ bc_format_17 bc_type bc_id bc_keyword_1 integer float float float CR_ 
      	  {
      	    card_read("BC_CARD",file_index); /* BC Format 17 */
     	    if ( common_BC_processing($3,$4,$5) )
    	    {
    	      BC_Types[iBC].BC_Data_Int[2] = $6;
      	      BC_Types[iBC].BC_Data_Int[0] = $7;
     	      BC_Types[iBC].BC_Data_Float[0] = $8; 
     	      BC_Types[iBC].BC_Data_Float[1] = $9;     	  
     	      BC_Types[iBC].BC_Data_Float[2] = $10;
     	      BC_Types[iBC].species_eq = BC_Types[iBC].BC_Data_Int[0];
     	    }
     	  }
	| BC_ EQUALS_ bc_format_18 bc_type bc_id bc_keyword_2 integer integer integer float CR_ 
  	  {
  	    card_read("BC_CARD",file_index); /* BC Format 18 */
      	    if ( common_BC_processing($3,$4,$5) )
    	    {
    	      BC_Types[iBC].BC_Data_Int[0] = $6;
      	      BC_Types[iBC].BC_Data_Int[1] = $7;
     	      BC_Types[iBC].BC_Data_Int[2] = $8; 
     	      BC_Types[iBC].BC_Data_Int[3] = $9;     	  
     	      BC_Types[iBC].BC_Data_Float[0] = $10;
     	      BC_Types[iBC].species_eq = BC_Types[iBC].BC_Data_Int[1];     	    
     	    }
  	  }   
	| BC_ EQUALS_ bc_format_19 bc_type bc_id floating_point_constant_list CR_ 
  	  {
  	    card_read("BC_CARD",file_index);  /* BC Format 19 */
    	    if ( common_BC_processing($3,$4,$5) )
    	    {
     	      BC_Types[iBC].len_u_BC = read_floats( &( BC_Types[iBC].u_BC) , 0 );
     	      /*
     	      BC_Types[iBC].u_BC = alloc_dbl_1(floating_point_constant_list_index, 0.0);
              for(j = 0; j < floating_point_constant_list_index; j++) 
              {
                BC_Types[iBC].u_BC[j] = floating_point_constant_list_array[j];
              }
    	      BC_Types[iBC].len_u_BC = floating_point_constant_list_index;
      	      floating_point_constant_list_index = -1;*/
      	    }
  	  }	   
	| BC_ EQUALS_ bc_format_20 bc_type bc_id integer floating_point_constant_list CR_ 
  	  {
  	    card_read("BC_CARD",file_index); /* BC Format 20 */
  	    if ( common_BC_processing($3,$4,$5) )
  	    {
	      BC_Types[iBC].species_eq = BC_Types[iBC].BC_Data_Int[0] = $6;
     	      BC_Types[iBC].len_u_BC = read_floats( &( BC_Types[iBC].u_BC), 0 );	      
    	      /*
    	      BC_Types[iBC].u_BC = alloc_dbl_1(floating_point_constant_list_index, 0.0);
              for(j = 0; j < floating_point_constant_list_index; j++) 
              {
                BC_Types[iBC].u_BC[j] = floating_point_constant_list_array[j];
              }
    	      BC_Types[iBC].len_u_BC = floating_point_constant_list_index;
      	      floating_point_constant_list_index = -1;        	    */
      	    }
    	  }	   
	| BC_ EQUALS_ bc_format_21 bc_type bc_id bc_id bc_keyword_3 floating_point_constant_list CR_ 
  	  {
  	    card_read("BC_CARD",file_index); /* BC Format 21 */
    	    if ( common_BC_processing($3,$4,$5) )
    	    {
    	      BC_Types[iBC].BC_Data_Int[0] = $6;
     	      BC_Types[iBC].len_u_BC = read_floats( &( BC_Types[iBC].u_BC), 0 );
  	      /*
  	      BC_Types[iBC].u_BC = alloc_dbl_1(floating_point_constant_list_index, 0.0);
              for(j = 0; j < floating_point_constant_list_index; j++) 
              {
                BC_Types[iBC].u_BC[j] = floating_point_constant_list_array[j];
              }
    	      BC_Types[iBC].len_u_BC = floating_point_constant_list_index;
      	      floating_point_constant_list_index = -1;     	    
      	      */
      	    }
  	  }	   	   
	| BC_ EQUALS_ bc_format_22 bc_type bc_id integer float optional_float CR_ /* this looks worng jjj */
  	{
  	  card_read("BC_CARD",file_index); /* BC Format 22 */
    	  if ( common_BC_processing($3,$4,$5) )
    	  {
    	    BC_Types[iBC].species_eq = BC_Types[iBC].BC_Data_Int[0] = $6;    
    	    BC_Types[iBC].BC_Data_Float[0] = $7;
            if(!strcmp(BC_Types[iBC].Set_Type, "NS"))
	    {
	      if ($7 == NA)
	      {
	        BC_Types[iBC].BC_relax = -1.0;
	      }
	      else
	      {
	        BC_Types[iBC].BC_relax = $7;
	      }
	    } 	      
	  }
  	}
	| BC_ EQUALS_ bc_format_23 bc_type bc_id integer float float CR_ 
  	{
  	  card_read("BC_CARD",file_index); /* BC Format 23 */
          if( common_BC_processing($3,$4,$5) )
    	  {
    	    BC_Types[iBC].species_eq = BC_Types[iBC].BC_Data_Int[0] = $6;
    	    BC_Types[iBC].BC_Data_Float[0] = $7;
     	    BC_Types[iBC].BC_Data_Float[1] = $8;
     	  }
  	}	   
	| BC_ EQUALS_ bc_format_24 bc_type bc_id integer float float float float float float float float CR_ 
  	{
  	  card_read("BC_CARD",file_index); /* BC Format 24 */
    	  if ( common_BC_processing($3,$4,$5) )
    	  {
    	    BC_Types[iBC].species_eq = BC_Types[iBC].BC_Data_Int[0] = $6;
    	    BC_Types[iBC].BC_Data_Float[0] = $7;
     	    BC_Types[iBC].BC_Data_Float[1] = $8;
     	    BC_Types[iBC].BC_Data_Float[2] = $9; 
     	    BC_Types[iBC].BC_Data_Float[3] = $10;     	  
     	    BC_Types[iBC].BC_Data_Float[4] = $11; 
     	    BC_Types[iBC].BC_Data_Float[5] = $12;       	  
     	    BC_Types[iBC].BC_Data_Float[6] = $13;       	  
    	    BC_Types[iBC].BC_Data_Float[7] = $14;       	  
  	  }
  	}	   
	| BC_ EQUALS_ bc_format_25 bc_type bc_id integer integer optional_float CR_ 
  	{
  	  card_read("BC_CARD",file_index); /* BC Format 25 */
    	  if ( common_BC_processing($3,$4,$5) )
    	  {
	    BC_Types[iBC].BC_Data_Int[0] = $6;
            BC_Types[iBC].BC_Data_Int[1] = $7;    	    
            /* It's not clear that the following code belongs here.  Some of it could be moved to post processing jss */
            /* POROUS_PRESSURE and DARCY_CONTINUOUS are  assumed to operate 
	     on the first species as it equations the Liquid phase pressure
	     to the continuous  phase pressure*/
  	    if(!strcmp(BC_Types[iBC].desc->name1, "POROUS_PRESSURE") ||
	       !strcmp(BC_Types[iBC].desc->name1, "DARCY_CONTINUOUS"))
	    {
	      BC_Types[iBC].species_eq = 0;
	    }
            /* add optional scaling float for fluid-structure bcs */
   	    if(!strcmp(BC_Types[iBC].desc->name1, "VELO_TANGENT_SOLID"))
	    {
	      /* these also need to be initialized */
	      BC_Types[iBC].BC_Data_Float[0] = 1.;  /*beta*/
	      BC_Types[iBC].BC_Data_Int[2] = 0;
	      BC_Types[iBC].BC_Data_Float[1] = 0.;
	    }
  	    if(!strcmp(BC_Types[iBC].desc->name1, "SOLID_FLUID_RS") ||
	      !strcmp(BC_Types[iBC].desc->name1, "SOLID_FLUID") ||
	      !strcmp(BC_Types[iBC].desc->name1, "FLUID_SOLID") ||
	      !strcmp(BC_Types[iBC].desc->name1, "FLUID_SOLID_RS"))
	    {
	      if ( $8 != NA )
	      {
                BC_Types[iBC].BC_Data_Float[0] = $8;
              }  
              else
              {
	        /* Default scaling factor to 1.0 */
	        WH(-1,"Defaulting Solid-Fluid Scaling Factor to 1.0");
                BC_Types[iBC].BC_Data_Float[0] = 1.0;               
              }
            }
 	    if( !strcmp(BC_Types[iBC].desc->name1, "FLUID_SOLID") && Linear_Solver == FRONT)
	     EH(-1,"Mucho Problemo with FLUID_SOLID boundary condition and frontal solver. USE lu");
	    if( (!strcmp(BC_Types[iBC].desc->name1, "SOLID_FLUID") ||
	      !strcmp(BC_Types[iBC].desc->name1, "SOLID_FLUID_RS")) && Linear_Solver == FRONT)
	      WH(-1,"Pacito (small?) Problemo with SOLID_FLUID boundary condition and frontal solver. Call PRS");          
          }
        }	   
	| BC_ EQUALS_ bc_format_26 bc_type bc_id integer integer float optional_integer optional_float CR_ 
  	{
  	  card_read("BC_CARD",file_index); /* BC Format 26 */
     	  if ( common_BC_processing($3,$4,$5) )
    	  {
	    BC_Types[iBC].BC_Data_Int[0] = $6;
            BC_Types[iBC].BC_Data_Int[1] = $7;        	    
    	    BC_Types[iBC].BC_Data_Float[0] = $8;
	    if ( (NA == $9) ||  ($10 == NA) )
	    {
              BC_Types[iBC].BC_Data_Float[1] = 0.0;
              BC_Types[iBC].BC_Data_Int[2] = -1;
            }  
            else
            {
              BC_Types[iBC].BC_Data_Float[1] = $10;
              BC_Types[iBC].BC_Data_Int[2] = $9;      
            }
          }    	    
    	}    	  	   
	| BC_ EQUALS_ bc_format_27 bc_type bc_id integer integer integer CR_ 
  	{
  	  card_read("BC_CARD",file_index); /* BC Format 27 */
    	  if ( common_BC_processing($3,$4,$5) )
    	  {
	    BC_Types[iBC].BC_Data_Int[0] = $6;
            BC_Types[iBC].BC_Data_Int[1] = $7;  
            BC_Types[iBC].BC_Data_Int[2] = $8;     	    	
          }
  	}	   
	| BC_ EQUALS_ bc_format_28 bc_type bc_id integer integer integer float CR_ 
  	{
  	  card_read("BC_CARD",file_index); /* BC Format 28 */
      	  if ( common_BC_processing($3,$4,$5) )
    	  {
    	    BC_Types[iBC].BC_Data_Int[0] = $6;
            BC_Types[iBC].BC_Data_Int[1] = $7;  
            BC_Types[iBC].BC_Data_Int[2] = $8;     
    	    BC_Types[iBC].BC_Data_Float[0] = $9;               	    
    	  }
  	}	   
	| BC_ EQUALS_ bc_format_29 bc_type bc_id integer integer integer float float CR_ 
 	{
 	  card_read("BC_CARD",file_index); /* BC Format 29 */
  	  if ( common_BC_processing($3,$4,$5) )
 	  {
	    BC_Types[iBC].BC_Data_Int[0] = $6;
            BC_Types[iBC].BC_Data_Int[1] = $7;  
            BC_Types[iBC].BC_Data_Int[2] = $8;     
    	    BC_Types[iBC].BC_Data_Float[0] = $9;
     	    BC_Types[iBC].BC_Data_Float[1] = $10;     
    	  }
    	}	   
	| BC_ EQUALS_ bc_format_30 bc_type bc_id integer float float float CR_ 
  	{
  	  card_read("BC_CARD",file_index); /* BC Format 30 */
      	  if ( common_BC_processing($3,$4,$5) )
    	  {
    	    BC_Types[iBC].BC_Data_Int[0] = $6;
      	    BC_Types[iBC].BC_Data_Float[0] = $7;
     	    BC_Types[iBC].BC_Data_Float[1] = $8;  
      	    BC_Types[iBC].BC_Data_Float[2] = $9;  
      	  }
  	}	      
  	/*$1  $2      $3        $4      $5    $6       $7      $8                    $9                         $10                         $11 */	
	| BC_ EQUALS_ table_bcs bc_type bc_id ordinate abcissa interpolation_method optional_file_equals_string optional_name_equals_string CR_ 
	/* BC = TABLE      SS {M} {ordinate} {abcissa}            {interpolation}                                        {File = name} */
	/* BC = TABLE      SS  4     X           V                   linear                                     FILE = newt_inlet_rot_table */
	/* BC = TABLE_WICS SS  16    X    SOLID_DISPLACEMENT2 	   {TABLE_SCALE} 	QUAD_GP 		FILE = 	FORCE1.TABLE */
	/* bc = TABLE_WICV SS  12    ZX   MESH_DISPLACEMENT1       BIQUADRATIC                                  FILE=load.table2 */
	{
	  card_read("TABLE_BC_CARD",file_index); /* BC Format 31 */
 	  if(num_BC_Tables < MAX_BC_TABLES)
	  {	  
	    if ( common_BC_processing($3,$4,$5) )
	    {  
	      allocate_table_bcs_description();
	      setup_bc_table($6,$7,$8);
	      if ( strcmp($9,"NA") )
	      {
	        table_type = BC_TABLE;
	        parse_table_file($9,$10,columns_expected);
    	        read_array_into_table(table_type, BC_Types[iBC].table, number_of_table_cards_read);
    	      }
	    }
	    else
	    {
	      accept_table_data_cards = FALSE;
	      table_type = NO_TABLE;
	    }
	  }
	  else
	  {
	    declare_error("Maximum TABLE_BCs exceeded.  Too many BC Tables."); accept_table_data_cards = FALSE; 
	    table_type = NO_TABLE;
	  }	  
    	}
   	| BC_ EQUALS_ bc_format_32 bc_type bc_id integer CR_ 
  	{
	  card_read("BC_CARD",file_index); /* BC Format 32 */
          if( common_BC_processing($3,$4,$5) )
          {
            BC_Types[iBC].BC_ID2 = $6;         
          }
  	}
	| BC_ EQUALS_ FIX_ bc_type bc_id var_name integer CR_
	{
	  card_read("FIX_BC_CARD",file_index); /* BC Format 33 */
 	  if( common_BC_processing($3,$4,$5) )
	  { 
    	    if ( valid_fix_bc_variable_name($6) )
    	    {
    	      BC_Types[iBC].BC_Data_Int[1] = $7;
    	      BC_Types[iBC].species_eq = BC_Types[iBC].BC_Data_Int[1];
    	      BC_Types[iBC].BC_relax = -1; 
    	    }
    	  }	  
	}
	| BC_ EQUALS_ GD_CONST_ bc_type bc_id equation_name integer var_name integer float CR_
	/*1   2       3         4       5     6             7       8        9        10     11     */
	{
	  card_read("GD_CONST_BC_CARD",file_index); /* BC Format 34 */
	  if(common_BC_processing($3,$4,$5) )
	  { 
    	    if ( process_equation_and_variable_names($6, $8) )
    	    {
    	      BC_Types[iBC].BC_Data_Int[0] = eq_index;
    	      BC_Types[iBC].BC_Data_Int[1] = $7;
    	      BC_Types[iBC].BC_Data_Int[2] = var_index;
    	      BC_Types[iBC].BC_Data_Int[3] = $9;    
    	      BC_Types[iBC].BC_Data_Float[0] = $10;   
    	      allocate_GD_bc_description(); 	      	      
    	    }	  
	  }
	}
	| BC_ EQUALS_ GD_LINEAR_ bc_type bc_id equation_name integer var_name integer float float CR_
	{
	  card_read("GD_LINEAR_BC_CARD",file_index); /* BC Format 35 */
 	  if( common_BC_processing($3,$4,$5) )
	  { 
    	    if ( process_equation_and_variable_names($6, $8) )
    	    {
    	      BC_Types[iBC].BC_Data_Int[0] = eq_index;
    	      BC_Types[iBC].BC_Data_Int[1] = $7;
    	      BC_Types[iBC].BC_Data_Int[2] = var_index;
    	      BC_Types[iBC].BC_Data_Int[3] = $9;    
    	      BC_Types[iBC].BC_Data_Float[0] = $10;   
    	      BC_Types[iBC].BC_Data_Float[1] = $11;       	           	      
    	      allocate_GD_bc_description();
    	    }	  
	  }
	}
	| BC_ EQUALS_ GD_PARAB_  bc_type bc_id equation_name integer var_name integer float float float CR_
	{
	  card_read("GD_PARAB_BC_CARD",file_index); /* BC Format 36 */
	  if( common_BC_processing($3,$4,$5) )
	  { 
    	    if ( process_equation_and_variable_names($6, $8) )
    	    {
    	      BC_Types[iBC].BC_Data_Int[0] = eq_index;
    	      BC_Types[iBC].BC_Data_Int[1] = $7;
    	      BC_Types[iBC].BC_Data_Int[2] = var_index;
    	      BC_Types[iBC].BC_Data_Int[3] = $9;    
    	      BC_Types[iBC].BC_Data_Float[0] = $10;   
    	      BC_Types[iBC].BC_Data_Float[1] = $11;   
    	      BC_Types[iBC].BC_Data_Float[2] = $12;      
    	      allocate_GD_bc_description();    	        	         	      
    	    }	  
	  }
	}
	| BC_ EQUALS_ GD_CIRC_  bc_type bc_id equation_name integer var_name integer float float float CR_
	{
	  card_read("GD_CIRC_BC_CARD",file_index); /* BC Format 37 */
	  if( common_BC_processing($3,$4,$5) )
	  { 
    	    if ( process_equation_and_variable_names($6, $8) )
    	    {
    	      BC_Types[iBC].BC_Data_Int[0] = eq_index;
    	      BC_Types[iBC].BC_Data_Int[1] = $7;
    	      BC_Types[iBC].BC_Data_Int[2] = var_index;
    	      BC_Types[iBC].BC_Data_Int[3] = $9;    
    	      BC_Types[iBC].BC_Data_Float[0] = $10;   
    	      BC_Types[iBC].BC_Data_Float[1] = $11;   
    	      BC_Types[iBC].BC_Data_Float[2] = $12;      	         	      
    	      allocate_GD_bc_description();    	        
    	    }	  
	  }
	}	  
	| BC_ EQUALS_ GD_POLYN_  bc_type bc_id equation_name integer var_name integer float float float float float float float CR_
	{
	  card_read("GD_POLYN_BC_CARD",file_index); /* BC Format 38 */
	  if( common_BC_processing($3,$4,$5) )
	  {
	    if ( process_equation_and_variable_names($6, $8) )  
	    { 
    	      BC_Types[iBC].BC_Data_Int[0] = eq_index;
    	      BC_Types[iBC].BC_Data_Int[1] = $7;
    	      BC_Types[iBC].BC_Data_Int[2] = var_index;
    	      BC_Types[iBC].BC_Data_Int[3] = $9;    
    	      BC_Types[iBC].BC_Data_Float[0] = $10;   
    	      BC_Types[iBC].BC_Data_Float[1] = $11;   
    	      BC_Types[iBC].BC_Data_Float[2] = $12;  
    	      BC_Types[iBC].BC_Data_Float[3] = $13;   
    	      BC_Types[iBC].BC_Data_Float[4] = $14;   
    	      BC_Types[iBC].BC_Data_Float[5] = $15;      	      
    	      BC_Types[iBC].BC_Data_Float[6] = $16;    	      
    	      allocate_GD_bc_description();    	        
    	    }
    	  }
	}
	| BC_ EQUALS_ GD_TIME_ bc_type bc_id equation_name integer gd_time_table_option integer optional_float optional_float optional_LINEAR optional_file_equals_string optional_name_equals_string CR_
        /*$1  $2      $3       $4      $5    $6                  $7       $8                    $9        $10           $11            $12                 $13                     $14                   $15 */
	/* GD_TIME card TABLE */
	{
	  card_read("GD_TIME_BC_CARD",file_index); /* BC Format 39 */
	  if ( num_BC_Tables < MAX_BC_TABLES )
	  {
	    if(common_BC_processing($3,$4,$5))
	    {
	      if ( process_equation_and_variable_names($6, $8) ) /*????*/
	      {
	        accept_table_data_cards = TRUE;
	        table_type = BC_TABLE;
	        columns_expected = 2;
    	        BC_Types[iBC].BC_Data_Int[0] = eq_index;
    	        BC_Types[iBC].BC_Data_Int[1] = $7;
    	        BC_Types[iBC].BC_Data_Int[2] = var_index;
    	        BC_Types[iBC].BC_Data_Int[3] = $9;   	          
                if ( strcmp($8,"TABLE") ) /* option is not TABLE */
                {
                  BC_Types[iBC].BC_Data_Float[0] = $10;	/* jjj need to check this */
                  BC_Types[iBC].BC_Data_Float[1] = $11;                    
                  accept_table_data_cards = FALSE;
                  table_type = NO_TABLE;
                }
                else  /* option is TABLE */
                {
	          setup_GD_table();    
		  if ( strcmp($13,"NA") )
		  {
		    accept_table_data_cards = FALSE;
		    table_type = BC_TABLE;
		    parse_table_file($13,$14,2);
    	            read_array_into_table(table_type, BC_Types[iBC].table, number_of_table_cards_read);
    	          }
    	          allocate_GD_bc_description();    	          
    	        }
              }  	          	            
    	    } 
    	    else 
    	    {
    	      accept_table_data_cards = FALSE;
    	      table_type = NO_TABLE;
    	    }
    	  } 
    	  else 
    	  {
    	    declare_error("Maximum TABLE_BCs exceeded.  Too many BC Tables."); accept_table_data_cards = FALSE; 
    	    table_type = NO_TABLE;
    	  }
	}
	| BC_ EQUALS_ GD_TABLE_ bc_type bc_id equation_name integer var_name integer float LINEAR_ optional_file_equals_string optional_name_equals_string CR_ 
	/* GD_TABLE card */
	{
 	  card_read("GD_TABLE_BC_CARD",file_index); /* BC Format 40 */
	  if ( num_BC_Tables < MAX_BC_TABLES )
	  {
	    if(common_BC_processing($3,$4,$5))
	    {
	      if ( process_equation_and_variable_names($6, $8) )
	      {
	        accept_table_data_cards = TRUE;
	        table_type = BC_TABLE;
	        columns_expected = 2;
    	        BC_Types[iBC].BC_Data_Int[0] = eq_index;
    	        BC_Types[iBC].BC_Data_Int[1] = $7;
    	        BC_Types[iBC].BC_Data_Int[2] = var_index;
    	        BC_Types[iBC].BC_Data_Int[3] = $9;      	
    	        BC_Types[iBC].BC_Data_Float[0] = $10;
	        setup_GD_table();      	       
		if ( strcmp($12,"NA") )
		{
		  parse_table_file($12,$13,2);
    	          read_array_into_table(table_type, BC_Types[iBC].table, number_of_table_cards_read);
    	        }
    	      allocate_GD_bc_description();    	          
	      }  	          	            
    	    } 
    	    else 
    	    {
    	      accept_table_data_cards = FALSE;
    	      table_type = NO_TABLE;
    	    }
    	  } 
    	  else 
    	  {
    	    declare_error("Maximum TABLE_BCs exceeded.  Too many BC Tables."); accept_table_data_cards = FALSE;  
    	    table_type = NO_TABLE;
    	  }	  	  
    	}
    	| BC_ EQUALS_ DISCONTINUOUS_VELO_ bc_type bc_id discontinuous_velo_bc_keyword integer integer CR_ 
    	{
  	  card_read("DISCONTINUOUS_VELO_BC_CARD",file_index);           /* BC Format 41 */
	  if ( common_BC_processing($3,$4,$5) )
	  {
	    BC_Types[iBC].BC_Data_Int[0] = $6;
	    BC_Types[iBC].BC_Data_Int[1] = $7;
	    BC_Types[iBC].BC_Data_Int[2] = $8;
	  }               	  
    	}
    	| BC_ EQUALS_ LATENT_HEAT_INTERNAL_ bc_type bc_id latent_heat_internal_bc_keyword integer integer float CR_
    	{
	  card_read("LATENT_HEAT_INTERNAL_BC_CARD",file_index);          /* BC Format 42 */
	  if ( common_BC_processing($3,$4,$5) )
	  {
	    BC_Types[iBC].BC_Data_Int[0] = $6;
	    BC_Types[iBC].BC_Data_Int[1] = $7;
	    BC_Types[iBC].BC_Data_Int[2] = $8;
	    BC_Types[iBC].BC_Data_Float[0] = $9;	   	  
	  }   	  
    	}    	
 	| table_data_card{}
	| end_table_card {}
	| file_equals_card {}
	| file_equals_with_name_equals_card {}
	| error {}
  	;
	
discontinuous_velo_bc_keyword:
	  EVAPORATION_ { $$ = EVAPORATION }
	| DISSOLUTION_ { $$ = DISSOLUTION }
	| error {}
	;
	
latent_heat_internal_bc_keyword:
          LIQUID_VAPOR_ { $$ = LIQUID_VAPOR }
        | SOLID_LIQUID_ { $$ = SOLID_LIQUID }
        | error {}
        ;

gd_time_table_option:
	  TABLE_       {}
	| EXPONENTIAL_ {}
	| SINUSOIDAL_  {}
	| LINEAR_      {}
	| error {}
	;	
	
file_equals_card:
        /*empty*/ {}
	| FILE_ EQUALS_ STRING_ CR_
	{
	  if (accept_table_data_cards)
	  {
	    card_read("FILE_EQUALS_CARD",file_index);
	    parse_table_file($3,"NA",columns_expected);
	    switch (table_type)
	    {
	      case BC_TABLE:
	      {
	        read_array_into_table(table_type, BC_Types[iBC].table, number_of_table_cards_read);
	        break;
	      }
	      case CONDUCTIVITY_TABLE:
	      {
	        read_array_into_table(table_type, mat_ptr->table, number_of_table_cards_read);
	        break;
	      }	      
	    }
	    accept_table_data_cards = FALSE;
	    table_type = NO_TABLE;
	  }
	  else
	  {
	    /* declare_warning("This card is being ignored.");*/
	  }
	}
	| error {}
	;
	
optional_file_equals_string:
	/* empty */ {strcpy($$,"NA");}
	| FILE_ EQUALS_ STRING_ {strcpy($$,$3);}
	| error {}
	;
	
optional_name_equals_string:
	/* empty */ {strcpy($$,"NA"); }
	| NAME_ EQUALS_ STRING_ {strcpy($$,$3); }
	| error {}
	;	
	
file_equals_with_name_equals_card:
        /* empty */ {}
        | FILE_ EQUALS_ STRING_ NAME_ EQUALS_ STRING_ CR_
	{
	  if (accept_table_data_cards)
	  {
	    card_read("FILE_EQUALS_WITH_NAME_EQUALS_CARD",file_index);
	    parse_table_file($3,$6,columns_expected);
	    switch (table_type)
	    {
	      case BC_TABLE:
	      {
	        read_array_into_table(table_type, BC_Types[iBC].table, number_of_table_cards_read);
	        break;
	      }
	      case CONDUCTIVITY_TABLE:
	      {
	        read_array_into_table(table_type, mat_ptr->table, number_of_table_cards_read);
	        break;
	      }	      
	    }	    
  	    accept_table_data_cards = FALSE;
	    table_type = NO_TABLE;
  	  }
  	  else
  	  {
	    /* declare_warning("This card is being ignored.");  	   */
  	  }
	}
	| error {}
	;
	  	
table_data_card:
	floating_point_constant_list CR_
	{ 
	  
          if (accept_table_data_cards )
          { 
            if (columns_expected == number_of_constants)
            {		
              /*card_read("table_data_CARD",file_index);*/
              /*table_array[number_of_table_cards_read][0] = atof($1);*/
              for(k=0;k<floating_point_constant_list_index;k++)
              {
                table_array[number_of_table_cards_read][k] = floating_point_constant_list_array[k];
              }
              floating_point_constant_list_index = -1;
              number_of_table_cards_read++;
              card_read("TABLE_DATA_CARD",file_index);
            }
            else
            { 	   
              sprintf(msg,"%d table columns expected, %d table columns found.", columns_expected, number_of_constants); 
              declare_error(msg);
              floating_point_constant_list_index = -1;
            }
          }
	  else 
	  {
	    declare_warning("TABLE data card being ignored."); 
	    floating_point_constant_list_index = -1;
	  }
	}
	| error {fprintf(parser_log, "--%s--", yytext); declare_error("Invalid Table data card."); }
	;
	
end_table_card:
	END_ TABLE_ CR_ 
	{ 
	  if (accept_table_data_cards)
	  { 
	    switch (table_type)
	    {   
	      case BC_TABLE:
	      {
	        read_array_into_table(table_type, BC_Types[iBC].table, number_of_table_cards_read);
	        break;
	      }
	      case CONDUCTIVITY_TABLE:
	      {
	        read_array_into_table(table_type, mat_ptr->table, number_of_table_cards_read);
	        break;
	      }	    
	      case DIFFUSIVITY_TABLE:
	      {
	        read_array_into_table(table_type, mat_ptr->table, number_of_table_cards_read);
	        break;
	      }	  	        
	    }
	    accept_table_data_cards = FALSE;
	    table_type = NO_TABLE;	    
	  }
	  else
	  {
	    /*declare_warning("END TABLE data card ignored.");*/
	  }
	  card_read("END_TABLE_CARD",file_index);
	}
	| error {}
	;
	
/* BCs which don't require any further data */
bc_format_1:
          HYDROSTATIC_SYMM_		{}
        | PSPG_ 			{}
        | Q_VELO_SLIP_	 		{}
        | error {}
        ;
         
/* BCs which require a single floating point and one additional optional parameter*/
bc_format_2:	
  	  QSIDE_ 			{}
	| TNRMLSIDE_ 			{}
	| TSHRSIDE_ 			{}
	| CURRENT_ 			{}
	| KINEMATIC_ 			{}
	| KINEMATIC_PETROV_ 		{}
	| KINEMATIC_COLLOC_ 		{}
	| KINEMATIC_DISC_ 		{}
	| T_MELT_ 			{}
	| VELO_NORMAL_ 			{}
	| VELO_NORM_COLLOC_ 		{}
	| VELO_NORMAL_DISC_ 		{}
	| SURFTANG_SCALAR_ 		{}
        | FLOW_STRESSNOBC_ 		{}
        | FLOW_GRADV_ 			{}
	| FILL_INLET_ 			{}
        | FILL_CA_ 			{}
        | FLOW_PRESSURE_ 		{}
        | STRONG_FILL_CA_ 		{}
        | error {}
        ;
        
         
/* BCs which require a second SSID and a single floating point and one additional optional parameter*/
bc_format_3:	
 	KINEMATIC_EDGE_ 		{}
 	| SURFTANG_SCALAR_EDGE_ 	{}
        | error {}
        ;
         
        
/* BCs which require a single floating and two optional parmeters*/
bc_format_4:
          DISTNG_ 		 	{}
	| DX_ 		 		{}
	| DX_RS_ 		 	{}
	| DXDISTNG_ 			{}
	| DY_ 		 		{}
	| DY_RS_ 			{}
	| DYDISTNG_ 			{}
	| DZ_ 				{}
	| DZ_RS_ 			{} 
	| DZDISTNG_ 			{} 
	| P_ 				{}
	| T_ 				{} 
	| U_ 				{} 
	| V_ 				{} 
	| W_ 				{} 
	| PU_ 				{}
	| PV_ 				{}
	| PW_ 				{}
	| S11_ 				{}
	| S12_ 				{}
	| S13_ 				{}
	| S22_ 				{}
	| S23_ 				{}
	| S33_ 				{}
	| S11_1_ 			{}
	| S12_1_	 		{}
	| S13_1_ 			{}
	| S22_1_ 			{}
	| S23_1_ 			{}
	| S33_1_ 			{}
	| S11_2_ 			{}
	| S12_2_ 			{}
	| S13_2_ 			{}
	| S22_2_ 			{}
	| S23_2_ 			{}
	| S33_2_ 			{}
	| S11_3_ 			{}
	| S12_3_ 			{}
	| S13_3_ 			{}
	| S22_3_ 			{}
	| S23_3_ 			{}
	| S33_3_ 			{}
	| S11_4_ 			{}
	| S12_4_	 		{}
	| S13_4_ 			{}
	| S22_4_ 			{}
	| S23_4_ 			{}
	| S33_4_	 		{}
	| S11_5_ 			{}
	| S12_5_ 			{}
	| S13_5_ 			{}
	| S22_5_	 		{}
	| S23_5_ 			{}
	| S33_5_ 			{}
	| S11_6_ 			{}
	| S12_6_	 		{}
	| S13_6_ 			{}
	| S22_6_ 			{}
	| S23_6_ 			{}
	| S33_6_	 		{}
	| S11_7_ 			{}
	| S12_7_ 			{}
	| S13_7_ 			{}
	| S22_7_	 		{}
	| S23_7_ 			{}
	| S33_7_ 			{}
	| G11_ 				{}
	| G12_	 			{}
	| G13_ 				{}
	| G21_ 				{}
	| G22_ 				{}
	| G23_	 			{}
	| G31_				{}
	| G32_				{}
	| G33_ 				{}
	| VOLT_ 			{}
	| F_	 			{}
	| SH_ 			 	{}
        | CONT_TANG_VEL_ 		{}
        | CONT_NORM_VEL_ 		{}
	| error {}
	;
	
/* BCs which require a second SSID and single floating */
bc_format_5:
	VELO_NORMAL_EDGE_ 		{}
	| VELO_NORMAL_EDGE_INT_ 	{}
	| CA_EDGE_CURVE_ 		{}
	| CA_EDGE_CURVE_INT_ 		{}
	| error{}
	;
	
/* BCs which require two floating points as data input */
bc_format_6:
          QCONV_ 			{} 
	| QRAD_ 			{}
	| KIN_LEAK_ 			{} 
	| VNORM_LEAK_ 			{} 
	| error {}
	;
	
/* BCs which require three floating point as data input */
bc_format_7:
	  FORCE_ 		 	{} 
	| FORCE_RS_ 			{} 
	| NORM_FORCE_ 			{} 
	| NORM_FORCE_RS_	 	{} 
	| SLOPEX_		 	{} 
	| SLOPEY_ 			{} 
	| SLOPEZ_ 			{} 
	| SLOPE_ 			{} 
        | CAPILLARY_ 			{}
	| error {}
	;
	
/* BCs which require a second SSID and three floating point as data input */
bc_format_8:
	VELO_TANGENT_EDGE_ 		{}
	| VELO_TANGENT_EDGE_INT_	{}
	| error {}
	;
		
/* BCs which require four floating point as data input */
bc_format_9:	   
	  PLANEX_ 			{}  
	| PLANEY_ 			{} 
	| PLANEZ_ 			{}
	| PLANE_ 			{}
	| CA_ 				{}
	| SURFTANG_ 			{}
        | FLOW_HYDROSTATIC_ 		{}
	| VELO_TANGENT_3D_ 		{}
	| VELO_SLIP_ 			{}
	| VELO_SLIP_ROT_		{}
	| VELO_SLIP_FILL_ 		{}
	| error {}
	;
	
/* BCs which require a second SSID and four floating point as data input */
bc_format_10:	   
	SURFTANG_EDGE_   		{}
	| CA_EDGE_ 			{}
	| CA_EDGE_INT_ 			{}
	| error {}
	;
	
/* BCs which require five floating point as data input */
bc_format_11:	
	  REP_FORCE_	  		{}
	| REP_FORCE_RS_ 		{}
	| error {}
	;
	
/* BCs which require 3 integers and five floating point values as data input  */
bc_format_16:	
	  VL_EQUIL_ 		 	{}
	| error {}
	;
	
/* BCs which require 1 keyword, 1 integer and 3 floating point values as data input */
bc_format_17:	   /* bc_keyword_1 are the associated valid keywords for this BC */
          YFLUX_EQUIL_ 		 	{}
        | error {}
        ;

/* keywords associated with BCs which require 1 keyword, 1 integer and 3 floating point values as data input */        
bc_keyword_1: 
  	  RAOULT_ 			{$$ = RAOULT;} 
	| FLORY_ 			{$$ = FLORY;}
	| error {}
	;
       
/* BCs which require 1 keyword, 3 integers and 1 floating point values as data input */
bc_format_18:	   
	  VL_POLY_		 	{}
	| error {}
	;
	
/* keywords associated with BCs which require 1 keyword, 4 integer and 1 floating point values as data input */        
bc_keyword_2: 
  	  VOLUME_ 			{$$ = VOLUME;} 
	| MASS_ 			{$$ = MASS;}
	| error {}
	;	
	
/* BCs which require seven floating points as data input */
bc_format_12:	   
	  CAP_REPULSE_	 		{}
        | CA_OR_FIX_ 			{}
	| MOVING_PLANE_ 		{}
	| error {}
	;
	
/* BCs which require a second SSID and eight floating points as data input */
bc_format_13:	   
	  VAR_CA_EDGE_ 			{}
	| error {}
	;
		 
/* BCs which require ten floating points as data input */
bc_format_14:	   
	  MOVING_CA_	 	 	{}
	| error {}
	;
	
/* BCs which require an arbitrary number of constants as data input */
bc_format_19:	   
          UVARY_ 			{} 
	| VVARY_ 			{} 
	| WVARY_ 		 	{} 
	| SPLINEX_			{}  
	| SPLINEY_	 		{} 
	| SPLINEZ_ 			{}
	| SPLINE_ 			{}
        | T_USER_			{}
        | VOLT_USER_			{}
        | UUSER_			{}
        | VUSER_ 			{}
        | WUSER_ 			{}
        | QUSER_ 			{}
	| FORCE_USER_ 			{}
	| FORCE_USER_RS_		{}
        | PRESSURE_USER_ 		{}
        | FLOW_PRESS_USER_ 		{}
        | KIN_CHEM_ 			{}
        | VAR_CA_USER_ 			{}
        | error {}
        ;
                        
/* BCs which require one integer plus an arbitrary number of double constants */ 
bc_format_20:	   
          YFLUX_USER_			{}
        | YUSER_ 			{}
        | Y2_ELECTRONEUTRALITY_ 	{}   
        | error {}
        ;
        
/* BCs which require a second SSID and Keyword plus an arbitrary number of double constants */
bc_format_21:	   
          CA_EDGE_OR_FIX_	 	{}
        | error {}
        ;
        
        /* keywords associated with BCs which require 1 keyword, 1 integer and multiple floating point values as data input */        
bc_keyword_3: 
  	  CIRCLE_ 			{$$ = CIRCLE;} 
	| USER_				{$$ = B_USER;}
	| error 
	{ if( (error_found_in_last_line |= 0) && (ProcID == 0) )fprintf(parser_log," << Keyword should be CIRCLE or USER. >>");}
	;
                
/* BCs which require one integer as data input */
bc_format_15:	   
	  POROUS_CONV_	 		{}
        | YFLUX_SUS_	 		{}
        | KIN_DISPLACEMENT_	 	{}
        | error {}
        ;
        
/* BCs which require one integer and one floating point and one optional floating point as data input */
bc_format_22:	   
	  Y_				{}
	| Y_DISCONTINUOUS_	 	{}
        | YFLUX_CONST_ 			{}
        | YTOTALFLUX_CONST_ 		{}  
	| POROUS_KIN_ 			{}
        | KINEMATIC_SPECIES_ 		{}
        | error {}
        ;
        
/* BCs which require one integer and two floating points */
bc_format_23:	   
        YFLUX_ 			{}
	| POROUS_FLUX_ 			{}
	| error {}
	;
	
/* BCs which require one integer and eight floating points */
bc_format_24:	   
        YFLUX_BV_	 		{}
	| error {}
	;	
	
/* BCs which require two integers as data input */
bc_format_25:	   
	  SOLID_FLUID_ 			{}
	| FLUID_SOLID_ 			{}
	| SOLID_FLUID_RS_ 		{}
	| FLUID_SOLID_RS_		{}
	| PENETRATION_ 			{}
	| POROUS_PRESSURE_ 		{}
	| DARCY_CONTINUOUS_		{}
	| NO_SLIP_ 			{}
	| NO_SLIP_RS_ 			{}
	| VELO_TANGENT_SOLID_ 		{}
	| error {}
	;
	
	
/* BCs which require 2 integers, then 1 float, followed by optional int. and float */
bc_format_26:	   
	  VELO_SLIP_SOLID_	 	{}
	| error {}
	;
	
	
/* BCs which require three integers as data input */
bc_format_27:	   
	  P_EQUIL_		 	{}
	| error {}
	;
	
	
/* BCs which require three integers and one float */
bc_format_28:	   
	  VN_POROUS_		 	{}
	| VP_EQUIL_ 			{}
	| error {}
	;
	
	
/* BCs which require three integers and two floats */ 
bc_format_29:	   
	  POROUS_GAS_			{}
	| error {}
	;
	
	
/* BCs which require one integer and three floats */
bc_format_30:	   
	  VELO_TANGENT_			{}
	| LATENT_HEAT_ 			{} 
	| error {}
	;	

/* BCc requiring two (string integer)s and a multiple number of floats */	   
table_bcs:	   
	  TABLE_WICV_		{ columns_expected = 3;}
        | TABLE_WICS_ 		{ columns_expected = 2;}
	| TABLE_ 		{ columns_expected = 2;}
	| error {}
	;
	
/* BCs which require two sideset side set IDs */       
bc_format_32:       
          COLLOCATE_EDGE_ 	{}
	| WEAK_INT_EDGE_	{}
	| STRONG_INT_EDGE_	{}	
	| error {}
	;
	        
bc_type:  NS_ {}
	| NC_ {}
	| SS_ {}
	| SC_ {}
	;	

bc_id:  INTEGER_ {$$=atoi($1);}
	| error {}
	;
	
equation_name:	/* Equation_Names Name in mm_names.h */
      R_MOMENTUM1_  	{} 
    | VX_  		{}
    | R_MOMENTUM2_  	{} 
    | VY_  		{}
    | R_MOMENTUM3_  	{} 
    | VZ_  		{}
    | R_ENERGY_  	{}
    | T_  		{}
    | R_MASS_  		{} 
    | Y_  		{}
    | R_MESH1_  	{} 
    | DMX_  		{}
    | R_MESH2_  	{} 
    | DMY_  		{}
    | R_MESH3_  	{} 
    | DMZ_  		{}
    | R_MASS_SURF_  	{} 
    | S_  		{}
    | R_PRESSURE_  	{} 
    | P_  		{}
    | R_STRESS11_  	{} 
    | S11_  		{}
    | R_STRESS12_  	{} 
    | S12_  		{}
    | R_STRESS22_  	{} 
    | S22_  		{}
    | R_STRESS13_  	{} 
    | S13_  		{}
    | R_STRESS23_ 	{} 
    | S23_  		{}
    | R_STRESS33_  	{} 
    | S33_  		{}
    | R_SOLID1_  	{} 
    | DMX_  		{}
    | R_SOLID2_  	{} 
    | DMY_RS_  		{}
    | R_SOLID3_  	{} 
    | DMZ_RS_  		{}
    | R_GRADIENT11_  	{} 
    | G11_  		{}
    | R_GRADIENT12_  	{} 
    | G12_  		{}
    | R_GRADIENT21_  	{} 
    | G21_  		{}
    | R_GRADIENT22_  	{} 
    | G22_  		{}
    | R_GRADIENT13_  	{} 
    | G13_  		{}
    | R_GRADIENT23_  	{} 
    | G23_  		{}
    | R_GRADIENT31_  	{} 
    | G31_  		{}
    | R_GRADIENT32_  	{} 
    | G32_  		{}
    | R_GRADIENT33_  	{} 
    | G33_  		{}
    | R_POTENTIAL_  	{} 
    | V_  		{} 
    | R_FILL_  		{} 
    | F_  		{} 
    | R_SHEAR_RATE_  	{} 
    | SH_  		{}
    | R_PMOMENTUM1_  	{} 
    | PVX_  		{}
    | R_PMOMENTUM2_  	{} 
    | PVY_  		{}
    | R_PMOMENTUM3_  	{} 
    | PVZ_  		{}
    | R_STRESS11_1_  	{} 
    | S11_1_  		{}
    | R_STRESS12_1_  	{} 
    | S12_1_  		{}
    | R_STRESS22_1_  	{} 
    | S22_1_  		{}
    | R_STRESS13_1_  	{} 
    | S13_1_  		{} 
    | R_STRESS23_1_  	{} 
    | S23_1_  		{}
    | R_STRESS33_1_  	{} 
    | S33_1_  		{} 
    | R_STRESS11_2_  	{} 
    | S11_2_  		{}
    | R_STRESS12_2_  	{} 
    | S12_2_  		{}
    | R_STRESS22_2_  	{} 
    | S22_2_  		{} 
    | R_STRESS13_2_  	{} 
    | S13_2_  		{}
    | R_STRESS23_2_  	{} 
    | S23_2_  		{}
    | R_STRESS33_2_  	{} 
    | S33_2_  		{} 
    | R_STRESS11_3_  	{} 
    | S11_3_  		{}
    | R_STRESS12_3_  	{} 
    | S12_3_  		{}
    | R_STRESS22_3_  	{} 
    | S22_3_  		{} 
    | R_STRESS13_3_  	{} 
    | S13_3_  		{}
    | R_STRESS23_3_  	{} 
    | S23_3_  		{}
    | R_STRESS33_3_  	{} 
    | S33_3_  		{} 
    | R_STRESS11_4_  	{} 
    | S11_4_  		{}
    | R_STRESS12_4_  	{} 
    | S12_4_  		{}
    | R_STRESS22_4_  	{} 
    | S22_4_  		{} 
    | R_STRESS13_4_  	{} 
    | S13_4_  		{}
    | R_STRESS23_4_  	{} 
    | S23_4_  		{}
    | R_STRESS33_4_  	{} 
    | S33_4_  		{} 
    | R_STRESS11_5_  	{} 
    | S11_5_  		{}
    | R_STRESS12_5_  	{} 
    | S12_5_  		{}
    | R_STRESS22_5_  	{} 
    | S22_5_  		{}
    | R_STRESS13_5_  	{} 
    | S13_5_  		{}
    | R_STRESS23_5_  	{} 
    | S23_5_  		{}
    | R_STRESS33_5_  	{} 
    | S33_5_  		{}
    | R_STRESS11_6_  	{} 
    | S11_6_  		{}
    | R_STRESS12_6_  	{} 
    | S12_6_  		{}
    | R_STRESS22_6_  	{} 
    | S22_6_  		{}
    | R_STRESS13_6_  	{} 
    | S13_6_  		{}
    | R_STRESS23_6_  	{} 
    | S23_6_  		{}
    | R_STRESS33_6_  	{} 
    | S33_6_  		{}
    | R_STRESS11_7_  	{} 
    | S11_7_  		{}
    | R_STRESS12_7_  	{} 
    | S12_7_  		{}
    | R_STRESS22_7_  	{} 
    | S22_7_  		{}
    | R_STRESS13_7_  	{} 
    | S13_7_  		{}
    | R_STRESS23_7_  	{} 
    | S23_7_  		{}
    | R_STRESS33_7_  	{} 
    | S33_7_  		{}
    | R_SPECIES_0_  	{} 
    | SP_0_  		{}
    | R_SPECIES_1_  	{} 
    | SP_1_  		{}
    | R_SPECIES_2_  	{} 
    | SP_2_  		{}
    | R_SPECIES_3_  	{} 
    | SP_3_  		{}
    | R_SPECIES_4_  	{} 
    | SP_4_  		{}
    | R_SPECIES_5_  	{} 
    | SP_5_  		{}
    | R_SPECIES_6_  	{} 
    | SP_6_  		{}
    | R_SPECIES_7_  	{} 
    | SP_7_  		{}
    | R_SPECIES_8_  	{} 
    | SP_8_  		{}
    | R_SPECIES_9_  	{} 
    | SP_9_  		{}
    | R_SPECIES_10_  	{} 
    | SP_10_  		{}
    | R_SPECIES_11_  	{} 
    | SP_11_  		{}
    | R_SPECIES_12_  	{} 
    | SP_12_  		{}
    | R_SPECIES_13_  	{} 
    | SP_13_  		{}
    | R_SPECIES_14_  	{} 
    | SP_14_  		{}
    | R_SPECIES_15_  	{} 
    | SP_15_  		{}
    | R_SPECIES_16_  	{} 
    | SP_16_  		{}
    | R_SPECIES_17_  	{} 
    | SP_17_  		{}
    | R_SPECIES_18_  	{} 
    | SP_18_  		{}
    | R_SPECIES_19_  	{} 
    | SP_19_  		{}
    | R_SPECIES_20_  	{} 
    | SP_20_  		{}
    | R_SPECIES_21_  	{} 
    | SP_21_  		{} 
    | R_SPECIES_22_  	{} 
    | SP_22_  		{} 
    | R_SPECIES_23_  	{} 
    | SP_23_ 		{} 
    | R_SPECIES_24_  	{} 
    | SP_24_  		{}
    | R_SPECIES_25_  	{} 
    | SP_25_  		{}
    | R_SPECIES_26_  	{} 
    | SP_26_  		{}
    | R_SPECIES_27_  	{} 
    | SP_27_  		{}
    | R_SPECIES_28_  	{} 
    | SP_28_  		{}
    | R_SPECIES_29_  	{} 
    | SP_29_  		{}
    | R_VOLFRACPH_0_  	{} 
    | VFP_0_  		{}
    | R_VOLFRACPH_1_  	{} 
    | VFP_1_  		{}
    | R_VOLFRACPH_2_  	{} 
    | VFP_2_  		{}
    | R_VOLFRACPH_3_  	{} 
    | VFP_3_  		{}
    | R_VOLFRACPH_4_  	{} 
    | VFP_4_  		{}
    | R_Y0_  		{} 
    | Y0_  		{}
    | R_Y1_  		{} 
    | Y1_  		{}
    | R_Y2_  		{} 
    | Y2_  		{}
    | R_Y3_  		{} 
    | Y3_  		{} 
    | R_Y4_  		{} 
    | Y4_  		{}
    | R_Y5_  		{} 
    | Y5_  		{}
    | R_Y6_  		{} 
    | Y6_  		{}
    | R_Y7_  		{} 
    | Y7_  		{}
    | R_Y8_  		{} 
    | Y8_  		{}
    | R_Y9_  		{} 
    | Y9_  		{}
    | R_Y10_  		{} 
    | Y10_  		{} 
    | R_Y11_  		{} 
    | Y11_  		{} 
    | R_Y12_  		{} 
    | Y12_  		{} 
    | R_Y13_  		{} 
    | Y13_  		{} 
    | R_Y14_  		{} 
    | Y14_  		{} 
    | R_Y15_  		{} 
    | Y15_  		{} 
    | R_Y16_  		{} 
    | Y16_  		{} 
    | R_Y17_  		{} 
    | Y17_  		{} 
    | R_Y18_  		{} 
    | Y18_  		{} 
    | R_Y19_  		{} 
    | Y19_  		{} 
    | R_Y20_  		{} 
    | Y20_  		{} 
    | R_Y21_  		{} 
    | Y21_  		{} 
    | R_Y22_  		{} 
    | Y22_  		{} 
    | R_Y23_  		{} 
    | Y23_  		{} 
    | R_Y24_  		{} 
    | Y24_  		{} 
    | R_Y25_  		{} 
    | Y25_  		{} 
    | R_Y26_  		{} 
    | Y26_  		{} 
    | R_Y27_  		{} 
    | Y27_  		{} 
    | R_Y28_  		{} 
    | Y28_  		{} 
    | R_Y29_  		{} 
    | Y29_  		{} 
    | R_MOM_NORMAL_  	{} 
    | DN_  		{} 
    | R_MOM_TANG1_  	{} 
    | DT1_  		{} 
    | R_MOM_TANG2_  	{} 
    | DT2_  		{} 
    | R_MESH_NORMAL_  	{} 
    | VN_  		{} 
    | R_MESH_TANG1_  	{} 
    | VT1_  		{} 
    | R_MESH_TANG2_  	{} 
    | VT2_  		{} 
    | error {}
    ;
    
  
post_var_name: /* post processing variable names */
    STREAM_			{pp_index= -1;} 
  | STREAM_NORMAL_STRESS_	{pp_index= -1;} 
  | CROSS_STREAM_SHEAR_		{pp_index= -1;} 
  | MEAN_SHEAR_			{pp_index= -1;} 
  | PRESSURE_CONT_		{pp_index= -1;} 
  | FILL_CONT_			{pp_index= -1;} 
  | FIRST_INVAR_STRAIN_		{pp_index= -1;} 
  | SEC_INVAR_STRAIN_		{pp_index= -1;} 
  | THIRD_INVAR_STRAIN_		{pp_index= -1;} 
  | DIV_VELOCITY_		{pp_index= -1;} 
  | VISCOSITY_			{pp_index= -1;} 
  | DENSITY_			{pp_index= -1;} 
  | NS_RESIDUALS_		{pp_index= -1;} 
  | MM_RESIDUALS_		{pp_index= -1;} 
  | DIFFUSION_VECTORS_		{pp_index= -1;} 
  | FLUXLINES_			{pp_index= -1;} 
  | CONDUCTION_VECTORS_		{pp_index= -1;} 
  | ENERGY_FLUXLINES_		{pp_index= -1;} 
  | TIME_DERIVATIVES_		{pp_index= -1;} 
  | STRESS_TENSOR_		{pp_index= -1;} 
  | REAL_STRESS_TENSOR_		{pp_index= -1;} 
  | STRAIN_TENSOR_		{pp_index= -1;} 
  | LAGRANGE_CONVECTION_	{pp_index= -1;} 
  | SURFACE_VECTORS_		{pp_index= -1;} 
  | ERROR_ZZ_VEL_		{pp_index= -1;} 
  | ERROR_ZZ_Q_			{pp_index= -1;} 
  | ERROR_ZZ_P_			{pp_index= -1;} 
  | USER_POST_			{pp_index= -1;} 
  | error {}
  ;

var_name:  /* Var_Names in mm_names.h */
      VELOCITY1_		{pp_index= VELOCITY1;}       
    | VX_			{pp_index= VELOCITY1;} 
    | VELOCITY2_ 		{pp_index= VELOCITY2;}        
    | VY_ 			{pp_index= VELOCITY2;}    
    | VELOCITY3_ 		{pp_index= VELOCITY3;} 
    | VZ_ 			{pp_index= VELOCITY3;}    
    | TEMPERATURE_ 		{pp_index= TEMPERATURE;}      
    | T_ 			{pp_index= TEMPERATURE;}    
    | MASS_FRACTION_ 		{pp_index= MASS_FRACTION;}      
    | Y_ 			{pp_index= MASS_FRACTION;}   
    | MESH_DISPLACEMENT1_ 	{pp_index= MESH_DISPLACEMENT1;} 
    | DMX_ 			{pp_index= MESH_DISPLACEMENT1;}   
    | MESH_DISPLACEMENT2_ 	{pp_index= MESH_DISPLACEMENT2;} 
    | DMY_ 			{pp_index= MESH_DISPLACEMENT2;}   
    | MESH_DISPLACEMENT3_ 	{pp_index= MESH_DISPLACEMENT3;}  
    | DMZ_ 			{pp_index= MESH_DISPLACEMENT3;}   
    | SURFACE_ 			{pp_index= SURFACE;}             
    | S_ 			{pp_index= SURFACE;}  
    | PRESSURE_ 		{pp_index= PRESSURE;}           
    | P_ 			{pp_index= PRESSURE;}   
    | POLYMER_STRESS11_ 	{pp_index= POLYMER_STRESS11;}   
    | S11_ 			{pp_index= POLYMER_STRESS11;}  
    | POLYMER_STRESS12_ 	{pp_index= POLYMER_STRESS12;}  
    | S12_ 			{pp_index= POLYMER_STRESS12;}   
    | POLYMER_STRESS22_ 	{pp_index= POLYMER_STRESS22;} 
    | S22_ 			{pp_index= POLYMER_STRESS22;}  
    | POLYMER_STRESS13_ 	{pp_index= POLYMER_STRESS13;}  
    | S13_ 			{pp_index= POLYMER_STRESS13;}   
    | POLYMER_STRESS23_ 	{pp_index= POLYMER_STRESS23;}  
    | S23_ 			{pp_index= POLYMER_STRESS23;}  
    | POLYMER_STRESS33_ 	{pp_index= POLYMER_STRESS33;}  
    | S33_ 			{pp_index= POLYMER_STRESS33;}  
    | SOLID_DISPLACEMENT1_ 	{pp_index= SOLID_DISPLACEMENT1;} 
    | DMX_RS_ 			{pp_index= SOLID_DISPLACEMENT1;}
    | SOLID_DISPLACEMENT2_ 	{pp_index= SOLID_DISPLACEMENT2;}  
    | DMY_RS_ 			{pp_index= SOLID_DISPLACEMENT2;}
    | SOLID_DISPLACEMENT3_ 	{pp_index= SOLID_DISPLACEMENT3;}  
    | DMZ_RS_ 			{pp_index= SOLID_DISPLACEMENT3;}  
    | VELOCITY_GRADIENT11_ 	{pp_index= VELOCITY_GRADIENT11;} 
    | G11_ 			{pp_index= VELOCITY_GRADIENT11;} 
    | VELOCITY_GRADIENT12_ 	{pp_index= VELOCITY_GRADIENT12;}  
    | G12_ 			{pp_index= VELOCITY_GRADIENT12;}
    | VELOCITY_GRADIENT21_ 	{pp_index= VELOCITY_GRADIENT21;}  
    | G21_ 			{pp_index= VELOCITY_GRADIENT21;}
    | VELOCITY_GRADIENT22_ 	{pp_index= VELOCITY_GRADIENT22;}  
    | G22_ 			{pp_index= VELOCITY_GRADIENT22;} 
    | VELOCITY_GRADIENT13_ 	{pp_index= VELOCITY_GRADIENT13;}  
    | G13_ 			{pp_index= VELOCITY_GRADIENT13;}
    | VELOCITY_GRADIENT23_ 	{pp_index= VELOCITY_GRADIENT23;}  
    | G23_ 			{pp_index= VELOCITY_GRADIENT23;}
    | VELOCITY_GRADIENT31_ 	{pp_index= VELOCITY_GRADIENT31;} 
    | G31_ 			{pp_index= VELOCITY_GRADIENT31;} 
    | VELOCITY_GRADIENT32_ 	{pp_index= VELOCITY_GRADIENT32;}  
    | G32_ 			{pp_index= VELOCITY_GRADIENT32;}
    | VELOCITY_GRADIENT33_ 	{pp_index= VELOCITY_GRADIENT33;}  
    | G33_ 			{pp_index= VELOCITY_GRADIENT33;} 
    | VOLTAGE_ 			{pp_index= VOLTAGE;}   
    | V_ 			{pp_index= VOLTAGE;}
    | FILL_ 			{pp_index= FILL;}       
    | F_ 			{pp_index= FILL;}
    | SHEAR_RATE_ 		{pp_index= SHEAR_RATE;}   
    | SH_ 			{pp_index= SHEAR_RATE;}
    | PVELOCITY1_ 		{pp_index= PVELOCITY1;}   
    | PVX_ 			{pp_index= PVELOCITY1;}
    | PVELOCITY2_ 		{pp_index= PVELOCITY2;}     
    | PVY_ 			{pp_index= PVELOCITY2;} 
    | PVELOCITY3_ 		{pp_index= PVELOCITY3;}     
    | PVZ_ 			{pp_index= PVELOCITY3;} 
    | POLYMER_STRESS11_1_ 	{pp_index= POLYMER_STRESS11_1;}  
    | S11_1_ 			{pp_index= POLYMER_STRESS11_1;} 
    | POLYMER_STRESS12_1_ 	{pp_index= POLYMER_STRESS12_1;} 
    | S12_1_ 			{pp_index= POLYMER_STRESS12_1;} 
    | POLYMER_STRESS22_1_ 	{pp_index= POLYMER_STRESS22_1;}  
    | S22_1_ 			{pp_index= POLYMER_STRESS22_1;} 
    | POLYMER_STRESS13_1_ 	{pp_index= POLYMER_STRESS13_1;}  
    | S13_1_ 			{pp_index= POLYMER_STRESS13_1;} 
    | POLYMER_STRESS23_1_ 	{pp_index= POLYMER_STRESS23_1;}
    | S23_1_ 			{pp_index= POLYMER_STRESS23_1;}
    | POLYMER_STRESS33_1_ 	{pp_index= POLYMER_STRESS33_1;} 
    | S33_1_ 			{pp_index= POLYMER_STRESS33_1;}  
    | POLYMER_STRESS11_2_ 	{pp_index= POLYMER_STRESS11_2;}  
    | S11_2_ 			{pp_index= POLYMER_STRESS11_2;} 
    | POLYMER_STRESS12_2_ 	{pp_index= POLYMER_STRESS12_2;}  
    | S12_2_ 			{pp_index= POLYMER_STRESS12_2;}  
    | POLYMER_STRESS22_2_ 	{pp_index= POLYMER_STRESS22_2;}  
    | S22_2_ 			{pp_index= POLYMER_STRESS22_2;} 
    | POLYMER_STRESS13_2_ 	{pp_index= POLYMER_STRESS13_2;}   
    | S13_2_ 			{pp_index= POLYMER_STRESS13_2;} 
    | POLYMER_STRESS23_2_ 	{pp_index= POLYMER_STRESS23_2;}  
    | S23_2_ 			{pp_index= POLYMER_STRESS23_2;} 
    | POLYMER_STRESS33_2_ 	{pp_index= POLYMER_STRESS33_2;}   
    | S33_2_ 			{pp_index= POLYMER_STRESS33_2;} 
    | POLYMER_STRESS11_3_ 	{pp_index= POLYMER_STRESS11_3;}  
    | S11_3_ 			{pp_index= POLYMER_STRESS11_3;} 
    | POLYMER_STRESS12_3_ 	{pp_index= POLYMER_STRESS12_3;}  
    | S12_3_ 			{pp_index= POLYMER_STRESS12_3;} 
    | POLYMER_STRESS22_3_ 	{pp_index= POLYMER_STRESS22_3;}  
    | S22_3_ 			{pp_index= POLYMER_STRESS22_3;} 
    | POLYMER_STRESS13_3_ 	{pp_index= POLYMER_STRESS13_3;}  
    | S13_3_ 			{pp_index= POLYMER_STRESS13_3;}  
    | POLYMER_STRESS23_3_ 	{pp_index= POLYMER_STRESS23_3;}   
    | S23_3_ 			{pp_index= POLYMER_STRESS23_3;}  
    | POLYMER_STRESS33_3_ 	{pp_index= POLYMER_STRESS33_3;}  
    | S33_3_ 			{pp_index= POLYMER_STRESS33_3;}  
    | POLYMER_STRESS11_4_ 	{pp_index= POLYMER_STRESS11_4;}  
    | S11_4_ 			{pp_index= POLYMER_STRESS11_4;}  
    | POLYMER_STRESS12_4_ 	{pp_index= POLYMER_STRESS12_4;} 
    | S12_4_ 			{pp_index= POLYMER_STRESS12_4;}  
    | POLYMER_STRESS22_4_ 	{pp_index= POLYMER_STRESS22_4;}  
    | S22_4_ 			{pp_index= POLYMER_STRESS22_4;}  
    | POLYMER_STRESS13_4_ 	{pp_index= POLYMER_STRESS13_4;} 
    | S13_4_ 			{pp_index= POLYMER_STRESS13_4;}  
    | POLYMER_STRESS23_4_ 	{pp_index= POLYMER_STRESS23_4;} 
    | S23_4_ 			{pp_index= POLYMER_STRESS23_4;}  
    | POLYMER_STRESS33_4_ 	{pp_index= POLYMER_STRESS33_4;} 
    | S33_4_ 			{pp_index= POLYMER_STRESS33_4;}  
    | POLYMER_STRESS11_5_ 	{pp_index= POLYMER_STRESS11_5;}  
    | S11_5_ 			{pp_index= POLYMER_STRESS11_5;} 
    | POLYMER_STRESS12_5_ 	{pp_index= POLYMER_STRESS12_5;} 
    | S12_5_ 			{pp_index= POLYMER_STRESS12_5;} 
    | POLYMER_STRESS22_5_ 	{pp_index= POLYMER_STRESS22_5;}  
    | S22_5_ 			{pp_index= POLYMER_STRESS22_5;} 
    | POLYMER_STRESS13_5_ 	{pp_index= POLYMER_STRESS13_5;} 
    | S13_5_ 			{pp_index= POLYMER_STRESS13_5;} 
    | POLYMER_STRESS23_5_ 	{pp_index= POLYMER_STRESS23_5;} 
    | S23_5_ 			{pp_index= POLYMER_STRESS23_5;} 
    | POLYMER_STRESS33_5_ 	{pp_index= POLYMER_STRESS33_5;}  
    | S33_5_ 			{pp_index= POLYMER_STRESS33_5;} 
    | POLYMER_STRESS11_6_ 	{pp_index= POLYMER_STRESS11_6;}  
    | S11_6_ 			{pp_index= POLYMER_STRESS11_6;} 
    | POLYMER_STRESS12_6_ 	{pp_index= POLYMER_STRESS12_6;}  
    | S12_6_ 			{pp_index= POLYMER_STRESS12_6;} 
    | POLYMER_STRESS22_6_ 	{pp_index= POLYMER_STRESS22_6;} 
    | S22_6_ 			{pp_index= POLYMER_STRESS22_6;} 
    | POLYMER_STRESS13_6_ 	{pp_index= POLYMER_STRESS13_6;} 
    | S13_6_ 			{pp_index= POLYMER_STRESS13_6;} 
    | POLYMER_STRESS23_6_ 	{pp_index= POLYMER_STRESS23_6;} 
    | S23_6_ 			{pp_index= POLYMER_STRESS23_6;} 
    | POLYMER_STRESS33_6_ 	{pp_index= POLYMER_STRESS33_6;}  
    | S33_6_ 			{pp_index= POLYMER_STRESS33_6;} 
    | POLYMER_STRESS11_7_ 	{pp_index= POLYMER_STRESS11_7;} 
    | S11_7_ 			{pp_index= POLYMER_STRESS11_7;} 
    | POLYMER_STRESS12_7_ 	{pp_index= POLYMER_STRESS12_7;}  
    | S12_7_ 			{pp_index= POLYMER_STRESS12_7;} 
    | POLYMER_STRESS22_7_ 	{pp_index= POLYMER_STRESS22_7;} 
    | S22_7_ 			{pp_index= POLYMER_STRESS22_7;} 
    | POLYMER_STRESS13_7_ 	{pp_index= POLYMER_STRESS13_7;}
    | S13_7_ 			{pp_index= POLYMER_STRESS13_7;} 
    | POLYMER_STRESS23_7_ 	{pp_index= POLYMER_STRESS23_7;} 
    | S23_7_ 			{pp_index= POLYMER_STRESS23_7;} 
    | POLYMER_STRESS33_7_ 	{pp_index= POLYMER_STRESS33_7;} 
    | S33_7_ 			{pp_index= POLYMER_STRESS33_7;} 
    | SPECIES_CONC_0_ 		{pp_index= SPECIES_UNK_0;}    
    | SP_0_ 			{pp_index= SPECIES_UNK_0;} 
    | SPECIES_CONC_1_ 		{pp_index= SPECIES_UNK_1;}   
    | SP_1_ 			{pp_index= SPECIES_UNK_1;}  
    | SPECIES_CONC_2_ 		{pp_index= SPECIES_UNK_2;}   
    | SP_2_ 			{pp_index= SPECIES_UNK_2;} 
    | SPECIES_CONC_3_ 		{pp_index= SPECIES_UNK_3;}  
    | SP_3_ 			{pp_index= SPECIES_UNK_3;}  
    | SPECIES_CONC_4_ 		{pp_index= SPECIES_UNK_4;}   
    | SP_4_ 			{pp_index= SPECIES_UNK_4;} 
    | SPECIES_CONC_5_ 		{pp_index= SPECIES_UNK_5;}   
    | SP_5_ 			{pp_index= SPECIES_UNK_5;} 
    | SPECIES_CONC_6_ 		{pp_index= SPECIES_UNK_6;}   
    | SP_6_ 			{pp_index= SPECIES_UNK_6;} 
    | SPECIES_CONC_7_ 		{pp_index= SPECIES_UNK_7;}  
    | SP_7_ 			{pp_index= SPECIES_UNK_7;} 
    | SPECIES_CONC_8_ 		{pp_index= SPECIES_UNK_8;}  
    | SP_8_ 			{pp_index= SPECIES_UNK_8;} 
    | SPECIES_CONC_9_ 		{pp_index= SPECIES_UNK_9;}  
    | SP_9_ 			{pp_index= SPECIES_UNK_9;}  
    | SPECIES_CONC_10_ 		{pp_index= SPECIES_UNK_10;} 
    | SP_10_ 			{pp_index= SPECIES_UNK_10;} 
    | SPECIES_CONC_11_ 		{pp_index= SPECIES_UNK_11;} 
    | SP_11_ 			{pp_index= SPECIES_UNK_11;} 
    | SPECIES_CONC_12_ 		{pp_index= SPECIES_UNK_12;}  
    | SP_12_ 			{pp_index= SPECIES_UNK_12;}  
    | SPECIES_CONC_13_ 		{pp_index= SPECIES_UNK_13;} 
    | SP_13_ 			{pp_index= SPECIES_UNK_13;}  
    | SPECIES_CONC_14_ 		{pp_index= SPECIES_UNK_14;}  
    | SP_14_ 			{pp_index= SPECIES_UNK_14;} 
    | SPECIES_CONC_15_ 		{pp_index= SPECIES_UNK_15;}   
    | SP_15_ 			{pp_index= SPECIES_UNK_15;} 
    | SPECIES_CONC_16_ 		{pp_index= SPECIES_UNK_16;}  
    | SP_16_ 			{pp_index= SPECIES_UNK_16;} 
    | SPECIES_CONC_17_ 		{pp_index= SPECIES_UNK_17;}   
    | SP_17_ 			{pp_index= SPECIES_UNK_17;}   
    | SPECIES_CONC_18_ 		{pp_index= SPECIES_UNK_18;}  
    | SP_18_ 			{pp_index= SPECIES_UNK_18;} 
    | SPECIES_CONC_19_ 		{pp_index= SPECIES_UNK_19;}  
    | SP_19_ 			{pp_index= SPECIES_UNK_19;}  
    | SPECIES_CONC_20_ 		{pp_index= SPECIES_UNK_20;} 
    | SP_20_ 			{pp_index= SPECIES_UNK_20;} 
    | SPECIES_CONC_21_ 		{pp_index= SPECIES_UNK_21;}  
    | SP_21_ 			{pp_index= SPECIES_UNK_21;} 
    | SPECIES_CONC_22_ 		{pp_index= SPECIES_UNK_22;}  
    | SP_22_ 			{pp_index= SPECIES_UNK_22;}   
    | SPECIES_CONC_23_ 		{pp_index= SPECIES_UNK_23;}  
    | SP_23_ 			{pp_index= SPECIES_UNK_23;}  
    | SPECIES_CONC_24_ 		{pp_index= SPECIES_UNK_24;} 
    | SP_24_ 			{pp_index= SPECIES_UNK_24;}  
    | SPECIES_CONC_25_ 		{pp_index= SPECIES_UNK_25;} 
    | SP_25_ 			{pp_index= SPECIES_UNK_25;}  
    | SPECIES_CONC_26_ 		{pp_index= SPECIES_UNK_26;}  
    | SP_26_ 			{pp_index= SPECIES_UNK_26;}   
    | SPECIES_CONC_27_ 		{pp_index= SPECIES_UNK_27;}  
    | SP_27_ 			{pp_index= SPECIES_UNK_27;} 
    | SPECIES_CONC_28_ 		{pp_index= SPECIES_UNK_28;} 
    | SP_28_ 			{pp_index= SPECIES_UNK_28;}  
    | SPECIES_CONC_29_ 		{pp_index= SPECIES_UNK_29;} 
    | SP_29_ 			{pp_index= SPECIES_UNK_29;}  
    | VOLFRACPH_0_ 		{pp_index= VOLF_PHASE_0 ;} 
    | VFP_0_ 			{pp_index= VOLF_PHASE_0 ;}  
    | VOLFRACPH_1_ 		{pp_index= VOLF_PHASE_1 ;}  
    | VFP_1_ 			{pp_index= VOLF_PHASE_1 ;} 
    | VOLFRACPH_2_ 		{pp_index= VOLF_PHASE_2 ;}   
    | VFP_2_ 			{pp_index= VOLF_PHASE_2 ;}  
    | VOLFRACPH_3_ 		{pp_index= VOLF_PHASE_3 ;}   
    | VFP_3_ 			{pp_index= VOLF_PHASE_3 ;}  
    | VOLFRACPH_4_ 		{pp_index= VOLF_PHASE_4 ;}   
    | VFP_4_ 			{pp_index= VOLF_PHASE_4 ;} 
    | MESH_POSITION1_ 		{pp_index= MESH_POSITION1;} 
    | X_ 			{pp_index= MESH_POSITION1;} 
    | MESH_POSITION2_ 		{pp_index= MESH_POSITION2;} 
    | Y_ 			{pp_index= MESH_POSITION2;} 
    | MESH_POSITION3_ 		{pp_index= MESH_POSITION3;} 
    | Z_ 			{pp_index= MESH_POSITION3;} 
    | VEL_NORM_ 		{pp_index= VEL_NORM;}      
    | VN_ 			{pp_index= VEL_NORM;}  
    | D_VEL1_DT_ 		{pp_index= D_VEL1_DT;}  
    | UDOT_ 			{pp_index= D_VEL1_DT;} 
    | D_VEL2_DT_ 		{pp_index= D_VEL2_DT;}  
    | VDOT_ 			{pp_index= D_VEL2_DT;}  
    | D_VEL3_DT_ 		{pp_index= D_VEL3_DT;} 
    | WDOT_ 			{pp_index= D_VEL3_DT;} 
    | D_T_DT_ 			{pp_index= D_T_DT;}    
    | TDOT_ 			{pp_index= D_T_DT;} 
    | D_C_DT_ 			{pp_index= D_Y_DT;}     
    | CDOT_ 			{pp_index= D_Y_DT;} 
    | D_X1_DT_ 			{pp_index= D_X1_DT;}    
    | XDOT_ 			{pp_index= D_X1_DT;} 
    | D_X2_DT_ 			{pp_index= D_X2_DT;}   
    | YDOT_ 			{pp_index= D_X2_DT;} 
    | D_X3_DT_ 			{pp_index= D_X3_DT;}   
    | ZDOT_ 			{pp_index= D_X3_DT;} 
    | D_S_DT_ 			{pp_index= D_S_DT;}   
    | SDOT_ 			{pp_index= D_S_DT;}  
    | D_P_DT_ 			{pp_index= D_P_DT;}   
    | PDOT_ 			{pp_index= D_P_DT;}     
    | SOLID_POSITION1_		{pp_index= SOLID_POSITION1;}     
    | X_RS_ 			{pp_index= SOLID_POSITION1;}     
    | SOLID_POSITION2_		{pp_index= SOLID_POSITION2;}     
    | Y_RS_ 			{pp_index= SOLID_POSITION2;}    
    | SOLID_POSITION3_		{pp_index= SOLID_POSITION3;}     
    | Z_RS_ 			{pp_index= SOLID_POSITION3;}         
    | POR_LIQ_PRES_		{pp_index= POR_LIQ_PRES;}
    | P_LIQ_			{pp_index= POR_LIQ_PRES;}
    | POR_GAS_PRES_		{pp_index= POR_GAS_PRES;}
    | P_GAS_ 			{pp_index= POR_GAS_PRES;}
    | POR_PORSITY_ 		{pp_index= POR_POROSITY;} 
    | P_POR_			{pp_index= POR_POROSITY;} 
    | POR_SATURATION_ 		{pp_index= POR_SATURATION;}
    | P_SAT_			{pp_index= POR_SATURATION;}
    | VORT_DIR1_ 		{pp_index= VORT_DIR1;}         
    | VDX_			{pp_index= VORT_DIR1;} 
    | VORT_DIR2_		{pp_index= VORT_DIR2;}
    | VDY_			{pp_index= VORT_DIR2;}
    | VORT_DIR3_		{pp_index= VORT_DIR3;}
    | VDZ_			{pp_index= VORT_DIR3;}
    | VORT_LAMBDA_		{pp_index= VORT_LAMBDA;}
    | VLAMBDA_         		{pp_index= VORT_LAMBDA;}
    | BIGX1_			{pp_index= BIGX1;}  
    | X1_			{pp_index= BIGX1;}
    | BIGX2_			{pp_index= BIGX2;} 
    | X2_			{pp_index= BIGX2;}
    | BIGX3_			{pp_index= BIGX3;} 
    | X3_			{pp_index= BIGX3;}
    | CURVATURE_		{pp_index= CURVATURE;} 
    | H_			{pp_index= CURVATURE;}
    | MESH_POSITION1_		{pp_index= MESH_POSITION1;} 
    | X_			{pp_index= MESH_POSITION1;}
    | MESH_POSITION2_		{pp_index= MESH_POSITION2;} 
    | Y_			{pp_index= MESH_POSITION2;}
    | MESH_POSITION3_		{pp_index= MESH_POSITION3;} 
    | Z_			{pp_index= MESH_POSITION3;}
    | VEL_NORM_			{pp_index= VEL_NORM;} 
    | VN_			{pp_index= VEL_NORM;}
    | D_VEL1_DT_		{pp_index= D_VEL1_DT;} 
    | UDOT_			{pp_index= D_VEL1_DT;}
    | D_VEL2_DT_		{pp_index= D_VEL2_DT;} 
    | VDOT_			{pp_index= D_VEL2_DT;}
    | D_VEL3_DT_		{pp_index= D_VEL3_DT;}
    | WDOT_			{pp_index= D_VEL3_DT;}
    | D_T_DT_			{pp_index= D_T_DT;} 
    | TDOT_			{pp_index= D_T_DT;}
    | D_C_DT_			{pp_index= D_Y_DT;} 
    | CDOT_			{pp_index= D_Y_DT;}
    | D_X1_DT_			{pp_index= D_X1_DT;} 
    | XDOT_			{pp_index= D_X1_DT;}
    | D_X2_DT_			{pp_index= D_X2_DT;} 
    | YDOT_			{pp_index= D_X2_DT;}
    | D_X3_DT_			{pp_index= D_X3_DT;}
    | ZDOT_			{pp_index= D_X3_DT;}
    | D_S_DT_			{pp_index= D_S_DT;} 
    | SDOT_			{pp_index= D_S_DT;}
    | D_P_DT_			{pp_index= D_P_DT;} 
    | PDOT_			{pp_index= D_P_DT;}
    | SOLID_POSITION1_		{pp_index= SOLID_POSITION1;} 
    | X_RS_			{pp_index= SOLID_POSITION1;}
    | SOLID_POSITION2_		{pp_index= SOLID_POSITION2;} 
    | Y_RS_			{pp_index= SOLID_POSITION2;}
    | SOLID_POSITION3_		{pp_index= SOLID_POSITION3;} 
    | Z_RS_			{pp_index= SOLID_POSITION3;}
    | error {}
    ;
	
	
rot_surface_card:
		/*empty */ {}
		| ROT_ EQUALS_ 
			rot_eq_type 			/* $3 is eq_type*/
			SURFACE_ 			/* $4 is always "SURFACE" */
			integer 			/* $5 is ss_id_1 */
			bc_name_or_rotation_string 	/* $6 is x_eqn */
			integer 			/* $7 is ss_idx */
			bc_name_or_rotation_string 	/* $8 is y_eqn */
			integer 			/* $9 is ss_idy */
			bc_name_or_rotation_string 	/* $10 is z_eqn */
			integer  			/* $11 is ss_idz*/
			rot_seed_method 		/* $12 is seed method*/
			optional_float 			/* $13 only occurs when seed method is SEED*/
			optional_float 			/* $14 only occurs when seed method is SEED*/
			optional_float 			/* $15 only occurs when seed method is SEED*/
			CR_ 
		{}
		| error {}
		;
		
rot_edge_card:
		/*empty */ {}
		| ROT_ EQUALS_ 
			rot_eqn_type 			/* $3 eq_type */	
			EDGE_ 				/* $4 always "EDGE" */
			integer 			/* $5 ss_id_1 */
			integer 			/* $6 ss_id_2 */
			bc_name_or_rotation_string 	/* $7 x_eqn */
			integer				/* $7 ss_id_x */
			bc_name_or_rotation_string 	/* $8 y_eqn */
			integer				/* $9 ss_id_y */
			bc_name_or_rotation_string	/* $10 z_eqn */
			integer				/* $11 ss_id_z */
			rot_seed_method 		/* $12 seed_method */
			optional_float 			/* $13 only occurs when seed_method is SEED*/
			optional_float 			/* $14 only occurs when seed_method is SEED*/
			optional_float 			/* $15 only occurs when seed_method is SEED*/
			CR_
		{}
		| error {}
		;
		
rot_vertex_card:
		/*empty */ {}
		| ROT_ EQUALS_ 
			rot_eqn_type 			/* $3 eq_type */	
			VERTEX_				/* $4 always "VERTEX" */
			integer 			/* $5 ss_id_1 */
			integer 			/* $6 ss_id_2 */
			integer				/* $7 ss_id_3 */
			bc_name_or_rotation_string 	/* $8 x_eqn */
			integer				/* $9 ss_id_x */
			bc_name_or_rotation_string 	/* $10 y_eqn */
			integer				/* $11 ss_id_y */
			bc_name_or_rotation_string	/* $12 z_eqn */
			integer				/* $13 ss_id_z */
			rot_seed_method 		/* $14 seed_method */
			optional_float 			/* $15 only occurs when seed_method is SEED*/
			optional_float 			/* $16 only occurs when seed_method is SEED*/
			optional_float 			/* $17 only occurs when seed_method is SEED*/
			CR_
		{}
		| error {}
		;

rot_volume_card:
		/*empty */ {}
		| ROT_ EQUALS_ 
			rot_eqn_type 			/* $3 eq_type */	
			VOLUME_ 			/* $4 always "VOLUME" */
			integer 			/* $5 ss_id_1 */
			integer 			/* $6 ss_id_2 */
			bc_name_or_rotation_string 	/* $7 x_eqn */
			integer				/* $7 ss_id_x */
			bc_name_or_rotation_string 	/* $8 y_eqn */
			integer				/* $9 ss_id_y */
			bc_name_or_rotation_string	/* $10 z_eqn */
			integer				/* $11 ss_id_z */
			rot_seed_method 		/* $12 seed_method */
			optional_float 			/* $13 only occurs when seed_method is SEED*/
			optional_float 			/* $14 only occurs when seed_method is SEED*/
			optional_float 			/* $15 only occurs when seed_method is SEED*/
			CR_
		{}
		| error {}
		;

interpolation_method:
		  LINEAR_ 	{$$=LINEAR_INT;}
		| QUADRATIC_	{$$=QUADRATIC_INT;}
		| QUAD_GP_ 	{$$=QUAD_GP_INT;}
		| BIQUADRATIC_	{$$=BIQUADRATIC_INT;}
		| error {}
		;
		
ordinate:
		  TIME_ {$$=TIME_ORD}
		| X_ 	{$$=X_ORD;}
		| Y_ 	{$$=Y_ORD;}
		| Z_ 	{$$=Z_ORD;}
		| XY_  	{$$=XY_ORD;}
		| XZ_ 	{$$=XZ_ORD;}
		| YZ_ 	{$$=YZ_ORD;}
		| YX_ 	{$$=YX_ORD;}
		| ZX_ 	{$$=ZX_ORD;}
		| ZY_ 	{$$=ZY_ORD;}
		| error {}
		;
		
abcissa:
		  VELOCITY1_ 		{$$=VELOCITY1_ABS;}
		| U_ 			{$$=VELOCITY1_ABS;}
		| VELOCITY2_ 		{$$=VELOCITY2_ABS;}
		| V_ 			{$$=VELOCITY2_ABS;}
		| VELOCITY3_ 		{$$=VELOCITY3_ABS;}
		| W_ 			{$$=VELOCITY3_ABS;}
		| TEMPERATURE_ 		{$$=TEMPERATURE_ABS;}
		| T_ 			{$$=TEMPERATURE_ABS;}
		| MASS_FRACTION_ 	{$$=MASS_FRACTION_ABS;}
		| Y_ 			{$$=MASS_FRACTION_ABS;}
		| SPECIES_ 		{$$=SPECIES_ABS;} /* ????*/
		| MESH_DISPLACEMENT1_ 	{$$=MESH_DISPLACEMENT1_ABS;}
		| DX_ 			{$$=MESH_DISPLACEMENT1_ABS;}
		| MESH_DISPLACEMENT2_ 	{$$=MESH_DISPLACEMENT2_ABS;}
		| DY_ 			{$$=MESH_DISPLACEMENT2_ABS;}
		| MESH_DISPLACEMENT3_ 	{$$=MESH_DISPLACEMENT3_ABS;}
		| DZ_  			{$$=MESH_DISPLACEMENT3_ABS;}
		| PRESSURE_ 		{$$=PRESSURE_ABS;}
		| P_ 		{$$=PRESSURE_ABS;}
		| SHEAR_RATE_ 	{$$=SHEAR_RATE_ABS;}
		| SH_  		{$$=SHEAR_RATE_ABS;}
		| S11_ 		{$$=S11_ABS;}
		| S12_ 		{$$=S12_ABS;}
		| S22_ 		{$$=S22_ABS;}
		| S13_ 		{$$=S13_ABS;}
		| S23_ 		{$$=S23_ABS;}
		| S33_ 		{$$=S33_ABS;}
		| S11_1_ 	{$$=S11_1_ABS;}
		| S12_1_ 	{$$=S12_1_ABS;}
		| S22_1_ 	{$$=S22_1_ABS;}
		| S13_1_ 	{$$=S13_1_ABS;}
		| S23_1_ 	{$$=S23_1_ABS;}
		| S33_1_ 	{$$=S33_1_ABS;}
		| S11_2_ 	{$$=S11_2_ABS;}
		| S12_2_ 	{$$=S12_2_ABS;}
		| S22_2_ 	{$$=S22_2_ABS;}
		| S13_2_ 	{$$=S13_2_ABS;}
		| S23_2_ 	{$$=S23_2_ABS;}
		| S33_2_ 	{$$=S33_2_ABS;}
		| S11_3_ 	{$$=S11_3_ABS;}
		| S12_3_ 	{$$=S12_3_ABS;}
		| S22_3_ 	{$$=S22_3_ABS;}
		| S13_3_ 	{$$=S13_3_ABS;}
		| S23_3_ 	{$$=S23_3_ABS;}
		| S33_3_ 	{$$=S33_3_ABS;}
		| S11_4_ 	{$$=S11_4_ABS;}
		| S12_4_ 	{$$=S12_4_ABS;}
		| S22_4_ 	{$$=S22_4_ABS;}
		| S13_4_ 	{$$=S13_4_ABS;}
		| S23_4_ 	{$$=S23_4_ABS;}
		| S33_4_ 	{$$=S33_4_ABS;}
		| S11_5_ 	{$$=S11_5_ABS;}
		| S12_5_ 	{$$=S12_5_ABS;}
		| S22_5_ 	{$$=S22_5_ABS;}
		| S13_5_ 	{$$=S13_5_ABS;}
		| S23_5_ 	{$$=S23_5_ABS;}
		| S33_5_ 	{$$=S33_5_ABS;}
		| S11_6_ 	{$$=S11_6_ABS;}
		| S12_6_ 	{$$=S12_6_ABS;}
		| S22_6_ 	{$$=S22_6_ABS;}
		| S13_6_ 	{$$=S13_6_ABS;}
		| S23_6_ 	{$$=S23_6_ABS;}
		| S33_6_ 	{$$=S33_6_ABS;}
		| S11_7_ 	{$$=S11_7_ABS;}
		| S12_7_ 	{$$=S12_7_ABS;}
		| S22_7_ 	{$$=S22_7_ABS;}
		| S13_7_ 	{$$=S13_7_ABS;}
		| S23_7_ 	{$$=S23_7_ABS;}
		| S33_7_ 	{$$=S33_7_ABS;}
		| error {}
		;
		
optional_float :
		/*empty*/ {$$=NA;}
		| float {}
		;

optional_LINEAR :
		/*empty*/ {strcpy($$,"NA");}
		| LINEAR_ {}
		;


optional_integer :
		/*empty*/ {$$==NA;}
		| INTEGER_ {$$=atoi($1);}
		;
integer :
		/*empty*/ {}
		| INTEGER_ {$$=atoi($1);}
		;
rot_eq_type:
		MESH_ {}
		| MOM_ {}
		;

bc_name_or_rotation_string:
		NONE_ {}
		| N_  {}
		| NO_ {}
		| N2_ {}
		| N3_ {}
		| T_  {}
		| T1_ {}
		| T2_ {}
		| B_  {}
		| S_  {}
		| X_  {}
		| Y_  {}
		| Z_  {}
		| GD_ {}
		| error {}
		;
		
rot_seed_method:
		UNSEEDED_
		| SEED_
		| NONE_ {}
		| BASIS_RESEED_
		| BASIS_ONCE_ 
		;
rot_eqn_type:
		MESH_ {}
		| MOM_ {}
		;
		
end_of_rot_card:
		/*empty */ {}
		| error {}
		;		
		
problem_desc_delimiter_card:
		/* empty */ {}
		| PROBLEM_ DESCRIPTION_  CR_
		{
		  card_read("PROBLEM_DESC_DELIMITER_CARD",file_index);
		  
		}
		| error {}
		;

number_of_materials_card:
		/* empty */ {} 
		| NUMBER_ OF_ MATERIALS_ EQUALS_ integer  CR_ 
		{
		  if( ($5 >0) && ($5 <= MAX_NUMBER_MATLS) )
		  {
		    upd->Num_Mat = $5;
		    card_read("NUMBER_OF_MATERIALS_CARD",file_index);
		    number_of_materials_expected = $5;
		  } 
		  else if ($5 == -1)
		  {
		    card_read("NUMBER_OF_MATERIALS_CARD",file_index);
                    accept_materials_until_end_of_material_card = TRUE;
		  }
		  else
		  {
		    declare_error("Invalid Number of Materials");
		    card_read("NO_VALID_CARD_READ",mn);	/* jjj ??? */	    
		  }
  
		}
		| error {}
		;

end_of_materials_card:
		/* empty */ {}
		| END_ OF_ MAT_ CR_ 
		{

	  	  if ( !accept_materials_until_end_of_material_card && (number_of_materials_expected !=  number_of_materials_found) )
		  { 
		    sprintf(msg, "Wrong number of materials found (%i expected, %i found).", number_of_materials_expected, number_of_materials_found ); 
		    declare_error(msg); 
	  	  }	  	  		
                  accept_materials_until_end_of_material_card = accept_material_related_cards = FALSE; 
                  card_read("END_OF_MAT_CARD",file_index);
		}
		| error {}
		;
		
mat_file_name_card:
		/* empty */ {} 
		| MAT_ EQUALS_ STRING_ integer  CR_
		{  
		  if ( (number_of_materials_expected > number_of_materials_found) || (accept_materials_until_end_of_material_card) )
		  {
                    if ($4>0) 
		    {
                      create_new_material($3, 1, $4, 0, 0);
		      card_read("MATERIAL_FILE_NAME_CARD",file_index);
		    } 
		    else 
		    {
		      declare_error("Invalid Block Identification Number");
		    }		    
		  }
		  else
		  {
		    accept_material_related_cards = FALSE;
		  }
		}
       		| MAT_ EQUALS_ STRING_ integer integer  CR_
		{
		  if ( (number_of_materials_expected > number_of_materials_found) || (accept_materials_until_end_of_material_card) )
		  {
	            if ($4>0 && $5>0)
		    {
		      create_new_material($3, 2, $4, $5, 0);
		      card_read("MATERIAL_FILE_NAME_CARD",file_index);
		    } 
		    else 
		    {
  		      declare_error("Invalid Block Identification Number");
		    }		    
		  }
		  else
		  {
		    accept_material_related_cards = FALSE;
		  }		  
		}
       		| MAT_ EQUALS_ STRING_ integer integer integer  CR_
		{
		  if ( (number_of_materials_expected > number_of_materials_found) || (accept_materials_until_end_of_material_card) )
		  {

                    if ($4>0 && $5>0 && $6>0)
                    {
		      card_read("MATERIAL_FILE_NAME_CARD",file_index); 
		      create_new_material($3, 3, $4, $5, $6);
		    } 
		    else 
		    {
  		      declare_error("Invalid Block Identification Number");
		    }		    
		  }
		  else
		  {
		    accept_material_related_cards = FALSE;
		  }		  
 		}
		| error {}
		;
		
coord_system_card:
 		/* empty */ {} 
		| COORDINATE_ SYSTEM_ EQUALS_ coord_system_type  CR_	
		{ 	  
		  if ( (pd_glob[mn] != NULL) && accept_material_related_cards)
		  {
		    card_read("COORDINATE_SYSTEM_CARD",file_index);
		    CoordinateSystem = $4;
		    upd->CoordinateSystem = $4;
		    pd_glob[mn]->CoordinateSystem = $4;
		  }
		  else
		  {
		    declare_warning("Coordinate system card is being ignored.");	
		  }

		}
		| error {}
		;

coord_system_type: CYLINDRICAL_ {$$ = CYLINDRICAL;	}
		|  SWIRLING_	{$$ = SWIRLING;		}
		|  CARTESIAN_ 	{$$ = CARTESIAN;	}
		|  SPHERICAL_ 	{$$ = SPHERICAL;	}
		|  PROJECTED_CARTESIAN_ {$$ = PROJECTED_CARTESIAN; }
		|  error{}
		;	
		
element_mapping_card:
		 /* empty */ {} 
		| ELEMENT_ MAPPING_ EQUALS_ element_mapping_type  CR_  	
		{
		  if ( (pd_glob[mn] != NULL) && accept_material_related_cards)
		  {
		    card_read("ELEMENT_MAPPING_CARD",file_index);
		    pd_glob[mn]->IntegrationMap = $4;
		  }
		  else
		  {
	            declare_warning("Element mapping card is being ignored.");	
		  }
		}
		| error {}
		;
		
element_mapping_type:	
		ISOPARAMETRIC_	{$$ = ISOPARAMETRIC; 	}
		| SUBPARAMETRIC_{$$ = SUBPARAMETRIC; 	} 
		| Q1_ 		{$$ = I_Q1;		}
		| Q2_ 		{$$ = I_Q2;		}
		| SP_ 		{$$ = I_SP;		}
		;

mesh_motion_card:
 		/* empty */ {} 
		| MESH_ MOTION_ EQUALS_ mesh_motion_type   CR_	
		{
		  if  ( (pd_glob[mn] != NULL) && accept_material_related_cards)
		  {
		    card_read("MESH_MOTION_CARD",file_index);
		    pd_glob[mn]->MeshMotion = $4;
		    cr_glob[mn]->MeshMotion = $4;   
		  }
		   else
		  { declare_warning("Mesh motion card is being ignored.");	}
		}
		| error {}
		;

mesh_motion_type:
		  ARBITRARY_ 	{$$ = ARBITRARY;	}
		| LAGRANGIAN_ 	{$$ = LAGRANGIAN;	}
		| TOTAL_ALE_	{$$ = TOTAL_ALE;	}
		;
	
bulk_species_card:
 		/* empty */ {} 
		| NUMBER_ OF_ BULK_ SPECIES_ EQUALS_ integer  CR_
		{
		  card_read("BULK_SPECIES_CARD",file_index); /* more exception processing needed here */
                  if (accept_material_related_cards)
                  { 
		    if ($6 >= 0)
		    {
		      if (pd_glob[mn] != NULL)
 		      {
 		        pd_glob[mn]->Num_Species_Eqn = pd_glob[mn]->Num_Species = $6;
    		      }
    		    } else {
    		      declare_error("Invalid Number of Bulk Species");
    		    }
    		  }
    		  else
    		  { 
    		    declare_warning("Number of buld species card is being ignored.");	
    		  }
		}
		| error {}
		;
		
really_meant_card:
 		/* empty */ {} 
		| REALLY_MEANT_ CR_
		{
		  card_read("REALLY_MEANT_CARD",file_index); 
		  if (accept_material_related_cards)
                  { 

    		  }
    		  else
    		  { 
    		    declare_warning("Really meant card is being ignored.");	
    		  }
		}
		| error {}
		;		
		
number_of_equations_card:
		/* empty */ {} 
		| NUMBER_ OF_ EQ_ EQUALS_ integer  CR_
		{
                  if(accept_material_related_cards)
                  {
		    if($5>0)
		    {
		      number_of_equations_expected[mn] = $5;
		      accept_eqs_till_end_of_eq_card = FALSE;
		    } 
		    else if ($5 == -1)
		    {
		      accept_eqs_till_end_of_eq_card = TRUE; 
		    }
		    else
		    {
                      sprintf(msg, "Invalid number of equations (%i).", $5 ); 
		      declare_error(msg); 
		      number_of_equations_expected[mn] = 0;
		      accept_eqs_till_end_of_eq_card = FALSE;		      		    		    
		    }
		  }
		  else 
		  {  
		    declare_warning("Number of equations card is being ignored.");	
		  }
                  card_read("NUMBER_OF_EQUATIONS_CARD",file_index);
		}
		| error {}
		;
		

mesh1_equation_card:
		/* empty */ {} 
		| EQ_ EQUALS_ MESH1_ mesh_equation_weight D1_ mesh_equation_weight float float float float float   CR_
		{
                  /*set_equation (card pointer,               mn,    name, wt,                var, intp,   ms, adv,  bnd, dif, src,  por );*/
		    set_equation ("MESH1_EQUATION_CARD", mn, R_MESH1, $4, MESH_DISPLACEMENT1,   $6,   $7,  $8,  $9, $10,  $11,  NA  )
		}
		| error {}
		;
		
mesh2_equation_card:
		/* empty */ {} 
		| EQ_ EQUALS_ MESH2_ mesh_equation_weight D2_ mesh_equation_weight float float float float float  CR_
		{		
		/*set_equation (mn,    name, wt,                var, intp,   ms, adv,  bnd, dif, src,  por );*/
		  set_equation("MESH2_EQUATION_CARD",mn, R_MESH2, $4, MESH_DISPLACEMENT2,   $6,   $7,  $8,  $9, $10,  $11,  NA  );
		}
		| error {}
		;
		
mesh3_equation_card:
		/* empty */ {} 
		| EQ_ EQUALS_ MESH3_ mesh_equation_weight D3_ mesh_equation_weight float float float float float  CR_
		{		
                /*set_equation(cards[],mn,    name, wt,                var, intp,   ms, adv,  bnd, dif, src,  por );*/
		  set_equation("MESH3_EQUATION_CARD" ,mn, R_MESH3, $4, MESH_DISPLACEMENT3,   $6,   $7,  $8,  $9, $10,  $11,  NA  );		  
		}
		| error {}
		;				

mesh_equation_weight:
		Q1_ {$$ = I_Q1;} | Q2_ {$$ = I_Q2;} ;
 

mom_solid1_equation_card:
		/* empty */ {} 
		| EQ_ EQUALS_ MOM_SOLID1_ mom_solid_equation_weight D1_RS_ mom_solid_equation_weight float float float float float  CR_
		{		
		   /*set_equation (mn,     name, wt,                 var, intp,   ms, adv,  bnd, dif, src,  por );*/
		    set_equation("MOM_SOLID1_EQUATION_CARD" ,mn, R_SOLID1, $4, SOLID_DISPLACEMENT1,   $6,   $7,  $8,  $9, $10,  $11,  NA  );
		}
		| error {}
		;

mom_solid2_equation_card:
		/* empty */ {} 
		| EQ_ EQUALS_ MOM_SOLID2_ mom_solid_equation_weight D2_RS_ mom_solid_equation_weight float float float float float  CR_
		{		
			  
		  /*set_equation (mn,     name, wt,                 var, intp,   ms, adv,  bnd, dif, src,  por );*/
		    set_equation("MOM_SOLID2_EQUATION_CARD",mn, R_SOLID2, $4, SOLID_DISPLACEMENT2,   $6,   $7,  $8,  $9, $10,  $11,  NA  );	
		}
		| error {}
		;
		
mom_solid3_equation_card:
		/* empty */ {} 
		| EQ_ EQUALS_ MOM_SOLID3_ mom_solid_equation_weight D3_RS_ mom_solid_equation_weight float float float float float  CR_
		{		

		  /*set_equation (mn,     name, wt,                 var, intp,   ms, adv,  bnd, dif, src,  por );*/
		    set_equation("MOM_SOLID3_EQUATION_CARD",mn, R_SOLID3, $4, SOLID_DISPLACEMENT3,   $6,   $7,  $8,  $9, $10,  $11,  NA  );		  

		}
		| error {}
		;
		
mom_solid_equation_weight:
		Q1_ {$$ = I_Q1;} | Q2_ {$$ = I_Q2;};


momentum1_equation_card:
		/* empty */ {} 
		| EQ_ EQUALS_ MOMENTUM1_ momentum_equation_weight U1_ momentum_equation_weight float float float float float float  CR_
		{		
	          /*set_equation (mn,        name, wt,       var, intp,   ms, adv, bnd, dif,  src,  por );*/
		    set_equation("MOMENTUM1_EQUATION_CARD",mn, R_MOMENTUM1, $4, VELOCITY1,   $6,   $7,  $8,  $9, $10,  $11,  $12 );		  
		}
		| error {}
		;
	
momentum2_equation_card:
		/* empty */ {} 
		| EQ_ EQUALS_ MOMENTUM2_ momentum_equation_weight U2_ momentum_equation_weight float float float float float float  CR_
		{		
		  /*set_equation (mn,        name, wt,       var, intp,   ms, adv, bnd, dif,  src,  por );*/
		    set_equation("MOMENTUM2_EQUATION_CARD",mn, R_MOMENTUM2, $4, VELOCITY2,   $6,   $7,  $8,  $9, $10,  $11,  $12 );		  

		}
		| error {}
		;
		
momentum3_equation_card:
		/* empty */ {} 
		| EQ_ EQUALS_ MOMENTUM3_ momentum_equation_weight U3_ momentum_equation_weight float float float float float float  CR_
		{		
		  /*set_equation (mn,        name, wt,       var, intp,   ms, adv, bnd, dif,  src,  por );*/
		    set_equation("MOMENTUM3_EQUATION_CARD",mn, R_MOMENTUM3, $4, VELOCITY3,   $6,   $7,  $8,  $9, $10,  $11,  $12 );		  
		}
		| error {}
		;
		
momentum_equation_weight:
		Q1_ {$$ = I_Q1;} | Q2_ {$$ = I_Q2;} | Q1_D_ {$$ = I_Q1_D;} | Q2_D_ {$$ = I_Q2_D;};
		
pmomentum1_equation_card:
		/* empty */ {} 
		| EQ_ EQUALS_ PMOMENTUM1_ momentum_equation_weight U1_ momentum_equation_weight float float float float float float  CR_
		{		
		  /*set_equation (mn,        name, wt,         var, intp,   ms, adv, bnd, dif,  src,  por );*/
		    set_equation("PMOMENTUM1_EQUATION_CARD",mn, R_PMOMENTUM1, $4, PVELOCITY1,   $6,   $7,  $8,  $9, $10,  $11,  $12 );		  	  
		}
		| error {}
		;
	
pmomentum2_equation_card:
		/* empty */ {} 
		| EQ_ EQUALS_ PMOMENTUM2_ momentum_equation_weight U2_ momentum_equation_weight float float float float float float  CR_
		{		
		  /*set_equation (mn,        name, wt,         var, intp,   ms, adv, bnd, dif,  src,  por );*/
		    set_equation("PMOMENTUM2_EQUATION_CARD" ,mn, R_PMOMENTUM2, $4, PVELOCITY2,   $6,   $7,  $8,  $9, $10,  $11,  $12 );		  
		}
		| error {}
		;
		
pmomentum3_equation_card:
		/* empty */ {} 
		| EQ_ EQUALS_ PMOMENTUM3_ momentum_equation_weight U3_ momentum_equation_weight float float float float float float  CR_
		{		
		  /*set_equation (mn,        name, wt,         var, intp,   ms, adv, bnd, dif,  src,  por );*/
		    set_equation("PMOMENTUM3_EQUATION_CARD" ,mn, R_PMOMENTUM3, $4, PVELOCITY3,   $6,   $7,  $8,  $9, $10,  $11,  $12 );		  	  
		}
		| error {}
		;
		
pmomentum_equation_weight:
		Q1_ {$$ = I_Q1;} | Q2_ {$$ = I_Q2;} | Q1_D_ {$$ = I_Q1_D;} | Q2_D_ {$$ = I_Q2_D;};


stress11_equation_card:
		/* empty */ {} 
		| EQ_ EQUALS_ STRESS11_ stress_equation_weight S11_ stress_equation_weight float float float float float  CR_
		{		
		  /*set_equation (mn,       name, wt,              var, intp,   ms, adv, bnd, dif,  src,  por );*/
		    set_equation("STRESS11_EQUATION_CARD",mn, R_STRESS11, $4, POLYMER_STRESS11,   $6,   $7,  $8,  $9, $10,  $11,  NA );		  
		}
		| error {}
		;

stress12_equation_card:
		/* empty */ {} 
		| EQ_ EQUALS_ STRESS12_ stress_equation_weight S12_ stress_equation_weight float float float float float  CR_
		{		
		  /*set_equation (mn,       name, wt,              var, intp,   ms, adv, bnd, dif,  src,  por );*/
		    set_equation("STRESS12_EQUATION_CARD" ,mn, R_STRESS12, $4, POLYMER_STRESS12,   $6,   $7,  $8,  $9, $10,  $11,  NA );		  
		  
		}
		| error {}
		;


stress13_equation_card:
		/* empty */ {} 
		| EQ_ EQUALS_ STRESS13_ stress_equation_weight S13_ stress_equation_weight float float float float float  CR_
		{		
 		  /*set_equation (mn,       name, wt,              var, intp,   ms, adv, bnd, dif,  src,  por );*/
		    set_equation("STRESS13_EQUATION_CARD" ,mn, R_STRESS13, $4, POLYMER_STRESS13,   $6,   $7,  $8,  $9, $10,  $11,  NA );		  
		    
		}
		| error {}
		;


stress22_equation_card:
		/* empty */ {} 
		| EQ_ EQUALS_ STRESS22_ stress_equation_weight S22_ stress_equation_weight float float float float float  CR_
		{		
		  /*set_equation (mn,       name, wt,              var, intp,   ms, adv, bnd, dif,  src,  por );*/
		    set_equation("STRESS22_EQUATION_CARD" ,mn, R_STRESS22, $4, POLYMER_STRESS22,   $6,   $7,  $8,  $9, $10,  $11,  NA );		  
		  
		}
		| error {}
		;


stress23_equation_card:
		/* empty */ {} 
		| EQ_ EQUALS_ STRESS23_ stress_equation_weight S23_ stress_equation_weight float float float float float  CR_
		{		
		  /*set_equation (mn,       name, wt,              var, intp,   ms, adv, bnd, dif,  src,  por );*/
		    set_equation("STRESS23_EQUATION_CARD" ,mn, R_STRESS23, $4, POLYMER_STRESS23,   $6,   $7,  $8,  $9, $10,  $11,  NA );		  
		  
		}
		| error {}
		;

stress_equation_weight:
		Q1_ {$$ = I_Q1;} | Q2_ {$$ = I_Q2;};

gradient11_equation_card:
		/* empty */ {} 
		/*1    2      3           4                        5    6                        7     8       9 */
		| EQ_ EQUALS_ GRADIENT11_ gradient_equation_weight G11_ gradient_equation_weight float float   CR_
		{		
		  /*set_equation (mn,       name, wt,                   var, intp,   ms, adv, bnd, dif,  src,  por );*/
		    set_equation("GRADIENT11_EQUATION_CARD",mn, R_GRADIENT11, $4, VELOCITY_GRADIENT11,   $6,   NA,  $7,  NA,  NA,   $8,   NA );		  
		  
		}
		| error {}
		;

gradient12_equation_card:
		/* empty */ {} 
		| EQ_ EQUALS_ GRADIENT12_ gradient_equation_weight G12_ gradient_equation_weight float float   CR_
		{		
		  /*set_equation (mn,       name, wt,                   var, intp,   ms, adv, bnd, dif,  src,  por );*/
		    set_equation("GRADIENT12_EQUATION_CARD",mn, R_GRADIENT12, $4, VELOCITY_GRADIENT12,   $6,   NA,  $7,  NA,  NA,   $8,   NA );		  
		  
		}
		| error {}
		;

gradient13_equation_card:
		/* empty */ {} 
		| EQ_ EQUALS_ GRADIENT13_ gradient_equation_weight G13_ gradient_equation_weight float float   CR_
		{		
		  /*set_equation (mn,       name, wt,                   var, intp,   ms, adv, bnd, dif,  src,  por );*/
		    set_equation("GRADIENT13_EQUATION_CARD" ,mn, R_GRADIENT13, $4, VELOCITY_GRADIENT13,   $6,   NA,  $7,  NA,  NA,   $8,   NA );		  
		  
		}
		| error {}
		;

gradient21_equation_card:
		/* empty */ {} 
		| EQ_ EQUALS_ GRADIENT21_ gradient_equation_weight G21_ gradient_equation_weight float float   CR_
		{		
		  /*set_equation (mn,       name, wt,                   var, intp,   ms, adv, bnd, dif,  src,  por );*/
		    set_equation("GRADIENT21_EQUATION_CARD",mn, R_GRADIENT21, $4, VELOCITY_GRADIENT21,   $6,   NA,  $7,  NA,  NA,   $8,   NA );		  
		  
		}
		| error {}
		;

gradient22_equation_card:
		/* empty */ {} 
		| EQ_ EQUALS_ GRADIENT22_ gradient_equation_weight G22_ gradient_equation_weight float float   CR_
		{		
		  /*set_equation (mn,       name, wt,                   var, intp,   ms, adv, bnd, dif,  src,  por );*/
		    set_equation("GRADIENT22_EQUATION_CARD",mn, R_GRADIENT22, $4, VELOCITY_GRADIENT22,   $6,   NA,  $7,  NA,  NA,   $8,   NA );		  
		  
		}
		| error {}
		;

gradient23_equation_card:
		/* empty */ {} 
		| EQ_ EQUALS_ GRADIENT23_ gradient_equation_weight G23_ gradient_equation_weight float float   CR_
		{		
		  /*set_equation (mn,       name, wt,                   var, intp,   ms, adv, bnd, dif,  src,  por );*/
		    set_equation("GRADIENT23_EQUATION_CARD",mn, R_GRADIENT23, $4, VELOCITY_GRADIENT23,   $6,   NA,  $7,  NA,  NA,   $8,   NA );		  
		  
		}
		| error {}
		;

gradient31_equation_card:
		/* empty */ {} 
		| EQ_ EQUALS_ GRADIENT31_ gradient_equation_weight G31_ gradient_equation_weight float float   CR_
		{		
		  /*set_equation (mn,       name, wt,                   var, intp,   ms, adv, bnd, dif,  src,  por );*/
		    set_equation("GRADIENT31_EQUATION_CARD",mn, R_GRADIENT31, $4, VELOCITY_GRADIENT31,   $6,   NA,  $7,  NA,  NA,   $8,   NA );		  
		  
		}
		| error {}
		;

gradient32_equation_card:
		/* empty */ {} 
		| EQ_ EQUALS_ GRADIENT32_ gradient_equation_weight G32_ gradient_equation_weight float float   CR_
		{		
		  /*set_equation (mn,       name, wt,                   var, intp,   ms, adv, bnd, dif,  src,  por );*/
		    set_equation("GRADIENT32_EQUATION_CARD",mn, R_GRADIENT32, $4, VELOCITY_GRADIENT32,   $6,   NA,  $7,  NA,  NA,   NA,   NA );		  
		  
		}
		| error {}
		;

gradient33_equation_card:
		/* empty */ {} 
		| EQ_ EQUALS_ GRADIENT33_ gradient_equation_weight G33_ gradient_equation_weight float float  CR_
		{		
		  /*set_equation (mn,       name, wt,                   var, intp,   ms, adv, bnd, dif,  src,  por );*/
		    set_equation("GRADIENT33_EQUATION_CARD",mn, R_GRADIENT33, $4, VELOCITY_GRADIENT33,   $6,   NA,  $7,  NA,  NA,   $8,   NA );		  
		  
		}
		| error {}
		;

gradient_equation_weight:
		Q1_ {$$ = I_Q1;} | Q2_ {$$ = I_Q2;};

voltage_equation_card:
		/* empty */ {} 
		| EQ_ EQUALS_ VOLTAGE_ voltage_equation_weight V_ voltage_equation_weight float float  CR_
		{		
		  /*set_equation (mn,        name, wt,     var, intp,   ms, adv, bnd, dif,  src,  por );*/
		    set_equation("VOLTAGE_EQUATION_CARD",mn, R_POTENTIAL, $4, VOLTAGE,   $6,   NA,  $7,  NA,  NA,   $8,   NA );		  
		}
		| error {}
		;
	
voltage_equation_weight:
		Q1_ {$$ = I_Q1;} | Q2_ {$$ = I_Q2;};
		
continuity_equation_card:
		/* empty */ {} 
		| EQ_ EQUALS_ CONTINUITY_ continuity_equation_weight P_ continuity_equation_weight float float  CR_
		{		
		  /*set_equation (mn,        name, wt,      var, intp,   ms, adv, bnd, dif,  src,  por );*/
		    set_equation("CONTINUITY_EQUATION_CARD",mn,  R_PRESSURE, $4, PRESSURE,   $6,   NA,  $7,  NA,  NA,   $8,   NA );		  
		  
		  
		}
		| error {}
		;		

continuity_equation_weight:
		P0_ {$$=I_P0;} | P1_ {$$=I_P1;}| Q1_ {$$=I_Q1;}| Q2_ {$$=I_Q2;};

energy_equation_card:
		/* empty */ {} 
		| EQ_ EQUALS_ ENERGY_ energy_equation_weight T_ energy_equation_weight float float float float float  CR_
		{		

		  
		  /*set_equation (mn,        name, wt,         var, intp,   ms, adv, bnd, dif,  src,  por );*/
		    set_equation("ENERGY_EQUATION_CARD",mn,    R_ENERGY, $4, TEMPERATURE,   $6,   $7,  $8,  $9, $10,  $11,   NA );		  
		}
		| error {}
		;		

energy_equation_weight:
		Q1_D_ {$$=I_Q1_D;}| Q2_D_ {$$=I_Q2_D;}| Q1_ {$$=I_Q1;}| Q2_ {$$=I_Q2;};

species_bulk_equation_card:
		/* empty */ {} 
		| EQ_ EQUALS_ SPECIES_BULK_ species_bulk_equation_weight Y_ species_bulk_equation_weight float float float float float  CR_
		{		
		  /*set_equation (mn,        name, wt,            var, intp,   ms, adv, bnd, dif,  src,  por );*/
		  set_equation("SPECIES_BULK_EQUATION_CARD",mn,      R_MASS, $4,  MASS_FRACTION,   $6,   $7,  $8,  $9,  $10,   NA, NA ); 
		  
		}
		| error {}
		;		

species_bulk_equation_weight:
		Q1_D_ {$$=I_Q1_D;}| Q2_D_ {$$=I_Q2_D;}| Q1_ {$$=I_Q1;}| Q2_ {$$=I_Q2;};

species_surf_equation_card:
		/* empty */ {} 
		| EQ_ EQUALS_ SPECIES_SURF_ species_surf_equation_weight  CR_ /* jjj what's the card format????*/
		{		
		  /*set_equation (name,                            mn, wt,          var,    intp,   ms,   adv, bnd, dif, src,  por );*/
		  /*  set_equation("SPECIES_SURF_EQUATION_CARD",mn, R_MASS_SURF, $4,     FILL,   $5,   NA,  NA,  NA,  NA,   NA,   NA );		  */
		}
		| error {}
		;		

species_surf_equation_weight:
		Q1_ {$$ = I_Q1;} | Q2_ {$$ = I_Q2;};

		level_set_equation_card:
		/* empty */ {} 
		| EQ_ EQUALS_ LEVEL_ SET_ level_set_equation_weight S_ level_set_equation_weight float float float  CR_
		{		
	          if ( Use_Level_Set != TRUE ) declare_error("Level Set Interface tracking must be turned on.");
      	          /* set explicit flag for fill equation and turn on level set switch*/
                  Explicit_Fill = 1;
  		  /*set_equation (mn,        name, wt,      var, intp,   ms, adv, bnd, dif,  src,  por ); jjj need more information on this card*/
		  /*  set_equation("LEVEL_ SET_EQUATION_CARD",mn, R_LEVEL_SET, $4,       NA,   NA,   NA,  NA,  NA,  NA,   NA,   NA );		  */
		}
		| error {}
		;		

level_set_equation_weight:
		Q1_ {$$ = I_Q1;} | Q2_ {$$ = I_Q2;} /* tom question jjj*/
		;
		
fill_equation_card:
		/* empty */ {} 
		| EQ_ EQUALS_ FILL_ fill_equation_weight F_ fill_equation_weight float float float  CR_
		{		
		  Explicit_Fill = 1;
  		  /*set_equation (mn,        name, wt,      var, intp,   ms, adv, bnd, dif,  src,  por ); need weights verified tom jjj */
		  /*  set_equation("FILL_EQUATION_CARD",mn,      R_FILL, $4,     FILL,   $5,   $6,  $7,  NA,  NA,   $8,   NA );		  */
		}
		| error {}
		;		

level_set_implicit_embedded_surface_card:
		/* empty */ {}
		| LEVEL_ SET_ IMPLICIT_ EMBEDDED_ SURFACE_ EQUALS_ yes_no_option CR_
		{
		    sprintf(msg, "The Level set implicit embedded surface card is no longer supported." ); 
		    declare_error(msg);
		}
		| error {}
		;

fill_equation_weight:
		PQ1_ {$$=I_PQ1;}/* tom question jjj */
		;

shear_rate_equation_card:
		/* empty */ {} 
		| EQ_ EQUALS_ SHEAR_RATE_ shear_rate_equation_weight SH_ shear_rate_equation_weight float float float  CR_
		{		
		  /*set_equation (mn,         name, wt,      var, intp,   ms, adv, bnd, dif,  src,  por ); need weights verified tom jjj */
		  /*  set_equation("SHEAR_RATE_EQUATION_CARD",mn, R_SHEAR_RATE, $4,     FILL,   $5,   NA,  $6,  NA,  $7,   $8,   NA );*/		  
		  /* jjj */
		}
		| error {}
		;		

shear_rate_equation_weight:
		Q1_ {$$ = I_Q1;} | Q2_ {$$ = I_Q2;};

species_unknown_equation_card:
		/* empty */ {} 
		/*1    2      3                4         */
		| EQ_ EQUALS_ SPECIES_UNKNOWN_ CR_
		{		
		  /* Harry M questions */;
  		  /*set_equation (mn,            name, wt,            var, intp,   ms, adv, bnd, dif,  src,  por ); jjj more work here */
		  /*   set_equation("SPECIES_UNKNOWN_EQUATION_CARD",mn, R_SPECIES_UNK_0, $4,  MASS_FRACTION,   $6,   $7,  $8,  $9,  $10,   NA, NA );	*/
		}
		| error {}
		;		

species_unknown_equation_weight:
		Q1_D_ {$$=I_Q1_D;}| Q2_D_ {$$=I_Q2_D;}| Q1_ {$$=I_Q1;}| Q2_ {$$=I_Q2;};
		
end_of_equations_card:
		/* empty */ {}
		| END_ OF_ EQ_ CR_
		{ 
	  	  if ( !accept_eqs_till_end_of_eq_card && (number_of_equations_expected[mn] !=  number_of_equations_found[mn]) )
		  { 
		    sprintf(msg, "Wrong number of equations found for material %i (%i expected, %i found).", mn+1, number_of_equations_expected[mn], number_of_equations_found[mn] ); 
		    declare_error(msg); 
	  	  }	  	  		
                  accept_eqs_till_end_of_eq_card = FALSE; 
                  card_read("END_OF_EQUATIONS_CARD",file_index);
		}
		| error {}
		;
		

post_processing_specification_section_delimiter_card:
		/* empty */ {} 
		| POST_ PROCESSING_ SPECIFICATIONS_ CR_
		{	
		  card_read("POST_PROCESSING_SPECIFICATION_SECTION_DELIMITER_CARD",file_index);		  
		}
		| error {}
		;
			
stream_function_card:
		/* empty */ {} 
		| STREAM_ FUNCTION_ EQUALS_ yes_no_option CR_
		{		
		  STREAM = $4;
		  card_read("STREAM_FUNCTION_CARD",file_index);		  
		}
		| error {}
		;
			
streamwise_normal_stress_card:
		/* empty */ {} 
		| STREAMWISE_ NORMAL_ STRESS_ EQUALS_ yes_no_option CR_
		{
		  STREAM_NORMAL_STRESS = $5;		
		  card_read("STREAMWISE_NORMAL_STRESS_CARD",file_index);		  
		}
		| error {}
		;
			
mean_shear_rate_card:
		/* empty */ {} 
		| MEAN_ SHEAR_ RATE_ EQUALS_ yes_no_option CR_
		{		
		  MEAN_SHEAR = $5;
		  card_read("MEAN_SHEAR_RATE_CARD",file_index);		  
		}
		| error {}
		;
			
pressure_contours_card:
		/* empty */ {} 
		| PRESSURE_ CONTOURS_ EQUALS_ yes_no_option CR_
		{
		  PRESSURE_CONT = $4;		
		  card_read("PRESSURE_CONTOURS_CARD",file_index);		  
		}
		| error {}
		;	

first_invariant_of_strain_card:
		/* empty */ {} 
		| FIRST_ INVARIANT_ OF_ STRAIN_ EQUALS_ yes_no_option CR_
		{
		  FIRST_INVAR_STRAIN = $6;		
  	   	  card_read("FIRST_INVARIANT_OF_STRAIN_CARD",file_index);		  
		}
		| error {}
		;	

second_invariant_of_strain_card:
		/* empty */ {} 
		| SECOND_ INVARIANT_ OF_ STRAIN_ EQUALS_ yes_no_option CR_
		{		
		  SEC_INVAR_STRAIN = $6;
		  card_read("SECOND_INVARIANT_OF_STRAIN_CARD",file_index);		  
		}
		| error {}
		;	

third_invariant_of_strain_card:
		/* empty */ {} 
		| THIRD_ INVARIANT_ OF_ STRAIN_ EQUALS_ yes_no_option CR_
		{		
		  THIRD_INVAR_STRAIN = $6;
		  card_read("THIRD_INVARIANT_OF_STRAIN_CARD",file_index);		  
		}
		| error {}
		;	

mesh_dilatation_card:
		/* empty */ {} 
		| MESH_ DILATATION_ EQUALS_ yes_no_option CR_
		{		
                  sprintf(msg, "The Mesh dilatation card is no longer being used.  Card being ignored."); 
		  declare_warning(msg);			  
		}
		| error {}
		;

velocity_divergence_card:
		/* empty */ {} 
		| VELOCITY_ DIVERGENCE_ EQUALS_ yes_no_option CR_
		{		
		  DIV_VELOCITY = $4;
		  card_read("VELOCITY_DIVERGENCE_CARD",file_index);		  
		}
		| error {}
		;	

navier_stokes_residuals_card:
		/* empty */ {} 
		| NAVIER_ STOKES_ RESIDUALS_ EQUALS_ yes_no_option CR_
		{		
		  NS_RESIDUALS = $5;
 		  card_read("NAVIER_STOKES_RESIDUALS_CARD",file_index);		  
		}
		| error {}
		;	

moving_mesh_residuals_card:
		/* empty */ {} 
		| MOVING_ MESH_ RESIDUALS_ EQUALS_ yes_no_option CR_
		{		
		  MM_RESIDUALS = $5;
		  card_read("MOVING_MESH_RESIDUALS_CARD",file_index);		  
		}
		| error {}
		;	

mass_diffusion_vectors_card:
		/* empty */ {} 
		| MASS_ DIFFUSION_ VECTORS_ EQUALS_ yes_no_option CR_
		{		
		  DIFFUSION_VECTORS = $5;
		  card_read("MASS_DIFFUSION_VECTORS_CARD",file_index);		  
		}
		| error {}
		;	

mass_fluxlines_card:
		/* empty */ {} 
		| MASS_ FLUXLINES_ EQUALS_ yes_no_option CR_
		{		
		  FLUXLINES = $4;
		  card_read("MASS_FLUXLINES_CARD",file_index);		  
		}
		| error {}
		;	

energy_conduction_vectors_card:
		/* empty */ {} 
		| ENERGY_ CONDUCTION_ VECTORS_ EQUALS_ yes_no_option CR_
		{		
		  CONDUCTION_VECTORS = $5;
		  card_read("ENERGY_CONDUCTION_VECTORS_CARD",file_index);		  
		}
		| error {}
		;	

energy_fluxlines_card:
		/* empty */ {} 
		| ENERGY_ FLUXLINES_ EQUALS_ yes_no_option CR_
		{		
		  ENERGY_FLUXLINES = $4;
		  card_read("ENERGY_FLUXLINES_CARD",file_index);		  
		}
		| error {}
		;	

time_derivatives_card:
		/* empty */ {} 
		| TIME_ DERIVATIVES_ EQUALS_ yes_no_option CR_
		{		
		  TIME_DERIVATIVES = $4;
		  card_read("TIME_DERIVATIVES_CARD",file_index);		  
		}
		| error {}
		;	

mesh_stress_tensor_card:
		/* empty */ {} 
		| MESH_ STRESS_ TENSOR_ EQUALS_ yes_no_option CR_
		{		
		  STRESS_TENSOR = $5;
		  card_read("MESH_STRESS_TENSOR_CARD",file_index);		  
		}
		| error {}
		;	

mesh_strain_tensor_card:
		/* empty */ {} 
		| MESH_ STRAIN_ TENSOR_ EQUALS_ yes_no_option CR_
		{		
		  STRAIN_TENSOR = $5;
		  card_read("MESH_STRAIN_TENSOR_CARD",file_index);		  
		}
		| error {}
		;	

real_solid_stress_tensor_card:
		/* empty */ {} 
		| REAL_ SOLID_ STRESS_ TENSOR_ EQUALS_ yes_no_option CR_
		{		
		  REAL_STRESS_TENSOR = $6;
		  card_read("REAL_SOLID_STRESS_TENSOR_CARD",file_index);		  
		}
		| error {}
		;	

porous_saturation_card:
		/* empty */ {} 
		| POROUS_ SATURATION_ EQUALS_ yes_no_option CR_
		{		
		  POROUS_SATURATION = $4;
		  card_read("POROUS_SATURATION_CARD",file_index);		  
		}
		| error {}
		;	


lagrangian_convection_card:
		/* empty */ {} 
		| LAGRANGIAN_ CONVECTION_ EQUALS_ yes_no_option CR_
		{		
		  LAGRANGE_CONVECTION = $4;
		  card_read("LAGRANGIAN_CONVECTION_CARD",file_index);		  
		}
		| error {}
		;	

normal_and_tangent_vectors_card:
		/* empty */ {} 
		| NORMAL_ AND_ TANGENT_ EQUALS_ yes_no_option CR_
		{		
		  SURFACE_VECTORS = $5;
		  card_read("NORMAL_AND_TANGENT_VECTORS_CARD",file_index);		  
		}
		| error {}
		;	
		
fill_contours_card:
		/* empty */ {}
		| FILL_ CONTOURS_ EQUALS_ yes_no_option CR_
		{		
		  FILL_CONT = $4;
		  card_read("FILL_CONTOURS_CARD",file_index);		  
		}
		| error {}
		;	
 
concentration_contours_card:
		/* empty */ {}
		| CONCENTRATION_ CONTOURS_ EQUALS_ yes_no_option CR_
		{		
		  CONC_CONT = $4;
		  card_read("CONCENTRATION_CONTOURS_CARD",file_index);		  
		}
		| error {}
		;	
		 
stress_contours_card:
		/* empty */ {}
		| STRESS_ CONTOURS_ EQUALS_ yes_no_option CR_
		{
		  STRESS_CONT = $4;		
		  card_read("STRESS_CONTOURS_CARD",file_index);		  
		}
		| error {}
		;	
		
particle_velocity_divergence_card:
		/* empty */ {}
		| PARTICLE_ VELOCITY_ DIVERGENCE_ EQUALS_ yes_no_option CR_
		{
		  DIV_PVELOCITY = $5;		
		  card_read("PARTICLE_VELOCITY_DIVERGENCE_CARD",file_index);		  
		}
		| error {}
		;	
		
total_velocity_divergence_card:
		/* empty */ {}
		| TOTAL_ VELOCITY_ DIVERGENCE_ EQUALS_ yes_no_option CR_
		{		
		  DIV_TOTAL = $5;
		  card_read("TOTAL_VELOCITY_DIVERGENCE_CARD",file_index);		  
		}
		| error {}
		;	
		
post_processing_viscosity_card:
		/* empty */ {}
		| VISCOSITY_ EQUALS_ yes_no_option CR_
		{		
		  VISCOSITY = $3;
		  card_read("VISCOSITY__CARD",file_index);		  
		}
		| error {}
		;	
		 
diffusive_mass_flux_vectors_card:
		/* empty */ {}
		| DIFFUSE_ MASS_ FLUX_ VECTORS_ EQUALS_ yes_no_option CR_
		{
		  DIFFUSION_VECTORS = $6;		
		  card_read("DIFFUSE_MASS_FLUX_VECTORS_CARD",file_index);		  
		}
		| error {}
		;	
		
real_solid_stress_tensor_card:
		/* empty */ {}
		| REAL_ SOLID_ STRESS_ TENSOR_ EQUALS_ yes_no_option CR_
		{		
		  REAL_STRESS_TENSOR = $6;
		  card_read("REAL_SOLID_STRESS_TENSOR_CARD",file_index);		  
		}
		| error {}
		;	
		
viscoplastic_def_grad_tensor_card:
		/* empty */ {}
		| VISCOPLASTIC_ DEF_ GRAD_ TENSOR_ EQUALS_ yes_no_option CR_
		{		
		  EVP_DEF_GRAD_TENSOR = $6;
		  card_read("VISCOPLASTIC_DEF_GRAD_TENSOR_CARD",file_index);		  
		}
		| error {}
		;	
		
error_zz_velocity_card:
		/* empty */ {}
		| ERROR_ ZZ_ VELOCITY_ EQUALS_ yes_no_option CR_
		{		
		  ERROR_ZZ_VEL = $5;
		  card_read("ERROR_ZZ_VELOCITY_CARD",file_index);		  
		}
		| error {}
		;	
		 
error_zz_heat_flux_card:
		/* empty */ {}
		| ERROR_ ZZ_ HEAT_ FLUX_ EQUALS_ yes_no_option CR_
		{
		  ERROR_ZZ_Q = $6;		
		  card_read("ERROR_ZZ_HEAT_FLUX_CARD",file_index);		  
		}
		| error {}
		;	
		
error_zz_pressure_card:
		/* empty */ {}
		| ERROR_ ZZ_ PRESSURE_ EQUALS_ yes_no_option CR_
		{		
		  ERROR_ZZ_P = $5;
		  card_read("ERROR_ZZ_PRESSURE_CARD",file_index);		  
		}
		| error {}
		;	
		
total_density_of_solvents_in_porous_media_card:
		/* empty */ {}
		| TOTAL_ DENSITY_ OF_ SOLVENTS_ IN_ POROUS_ MEDIA_ EQUALS_ yes_no_option CR_
		{		
		  POROUS_RHO_TOTAL_SOLVENTS = $9;
		  card_read("TOTAL_DENSITY_OF_SOLVENTS_IN_POROUS_MEDIA_CARD",file_index);		  
		}
		| error {}
		;	
		
density_of_solvents_in_gas_phase_in_porous_media_card:
		/* empty */ {}
		| DENSITY_ OF_ SOLVENT_ IN_ GAS_ PHASE_ IN_ POROUS_ MEDIA_ EQUALS_ yes_no_option CR_
		{		
		  POROUS_RHO_GAS_SOLVENTS = $11;
		  card_read("DENSITY_OF_SOLVENTS_IN_GAS_PHASE_IN_POROUS_MEDIA_CARD",file_index);		  
		}
		| error {}
		;	
		
density_of_liquid_phase_in_porous_media_card:
		/* empty */ {}
		| DENSITY_ OF_ LIQUID_ PHASE_ IN_ POROUS_ MEDIA_ EQUALS_ yes_no_option CR_
		{		
		  POROUS_RHO_LPHASE = $9;
		  card_read("DENSITY_OF_LIQUID_PHASE_IN_POROUS_MEDIA_CARD",file_index);		  
		}
		| error {}
		;	
		
gas_phase_darcy_velocity_in_porous_media_card:
		/* empty */ {}
		| GAS_ PHASE_ DARCY_ VELOCITY_ IN_ POROUS_ MEDIA_ EQUALS_ yes_no_option CR_
		{		
		  DARCY_VELOCITY_GAS = $9;
		  card_read("GAS_PHASE_DARCY_VELOCITY_IN_POROUS_MEDIA_CARD",file_index);		  
		}
		| error {}
		;	
		
liquid_phase_darcy_velocity_in_porous_media_card:
		/* empty */ {}
		| LIQUID_ PHASE_ DARCY_ VELOCITY_ IN_ POROUS_ MEDIA_ EQUALS_ yes_no_option CR_
		{		
		  DARCY_VELOCITY_LIQ = $9;
		  card_read("LIQUID_PHASE_DARCY_VELOCITY_IN_POROUS_MEDIA_CARD",file_index);		  
		}
		| error {}
		;	
		
grid_peclet_number_in_porous_media_card:
		/* empty */ {}
		| GRID_ PECLET_ NUMBER_ IN_ POROUS_ MEDIA_ EQUALS_ yes_no_option CR_
		{		
		  POROUS_GRIDPECLET = $8;
		  card_read("GRID_PECLET_NUMBER_IN_POROUS_MEDIA_CARD",file_index);		  
		}
		| error {}
		;	
		
SUPG_velocity_in_porous_media_card:
		/* empty */ {}
		| SUPG_ VELOCITY_ IN_ POROUS_ MEDIA_ EQUALS_ yes_no_option CR_
		{		
		  POROUS_SUPGVELOCITY = $7;
		  card_read("SUPG_VELOCITY_IN_POROUS_MEDIA_CARD",file_index);		  
		}
		| error {}
		;	
		
cross_stream_shear_rate_card:
		/* empty */ {}
		| CROSSSTREAM_ SHEAR_ RATE_ EQUALS_ yes_no_option CR_
		{		
		  CROSS_STREAM_SHEAR = $5;
		  card_read("CROSSSTREAM_SHEAR_RATE_CARD",file_index);		  
		}
		| error {}
		;			
		
vorticity_vector_card:
		/* empty */ {}
		| VORTICITY_ VECTOR_ EQUALS_ yes_no_option CR_
		{		
		  CURL_V = $4;
		  card_read("VORTICITY_VECTOR_CARD",file_index);		  
		}
		| error {}
		;	
		
capillary_pressure_in_porous_media_card:
		/* empty */ {}
		| CAPILLARY_ PRESSURE_ IN_ POROUS_ MEDIA_ EQUALS_ yes_no_option CR_
		{		
		  CAPILLARY_PRESSURE = $7;
		  card_read("CAPILLARY_PRESSURE_IN_POROUS_MEDIA_CARD",file_index);		  
		}
		| error {}
		;

user_defined_post_processing_card:
		/* empty */ {} 
		| USER_DEFINED_ POST_ PROCESSING_ EQUALS_ yes_no_option optional_floating_point_constant_list CR_ 
		{
		  USER_POST = $5;
		  if ( floating_point_constant_list_index > -1)
	          {
	            u_post_proc = (dbl *)array_alloc(1, number_of_constants, sizeof(dbl));
	            read_floats( &(*u_post_proc), 0);
    		  }
		  card_read("USER_DEFINED_POST_PROCESSING_CARD",file_index);	  
		}
		| error {}
		;	

floating_point_constant_list:
		  float {floating_point_constant_list_index++; number_of_constants = floating_point_constant_list_index+1; floating_point_constant_list_array[floating_point_constant_list_index] = $1; }
		| floating_point_constant_list float {floating_point_constant_list_index++;number_of_constants = floating_point_constant_list_index+1; floating_point_constant_list_array[floating_point_constant_list_index] = $2; }
		| error { if( (error_found_in_last_line |= 0) && (ProcID == 0) ) fprintf(parser_log," << Invalid floating point constant.>>");}
		;

optional_floating_point_constant_list:
		  /* empty , a floating_point_constant_list_index of -1 indicated an empty array*/ {floating_point_constant_list_index = -1;}
		|  float {floating_point_constant_list_index++; number_of_constants = floating_point_constant_list_index+1; floating_point_constant_list_array[floating_point_constant_list_index] = $1; }
		| floating_point_constant_list float {floating_point_constant_list_index++; number_of_constants = floating_point_constant_list_index+1; floating_point_constant_list_array[floating_point_constant_list_index] = $2; }
		| error { if( (error_found_in_last_line |= 0) && (ProcID == 0) ) fprintf(parser_log," << Invalid floating point constant.>>");}
		;

error_zz_element_size_card:
		/* empty */ {} 
		| ERROR_ ZZ_ zz_error_type ELEM_ SIZE_ EQUALS_ float float float float float float CR_
		{	
		  if ( !strcmp($3, "VELOCITY") )
		  {
		    if ( ERROR_ZZ_VEL == 1 )
		    { 
		      common_error_element_size_processing( $7, $8, $9, $10, $11, $12 ); 
		    }
		    else
		    {
		      declare_error("Error ZZ velocity elem size card requires 'Error ZZ velocity = YES' card.  Please add.");
		    }
		  }
		  else if ( !strcmp($3, "HEAT") )
		  {
		    if ( ERROR_ZZ_Q == 1)
		    {
		      common_error_element_size_processing( $7, $8, $9, $10, $11, $12 );
		    }
		    else
		    {
		      declare_error("Error ZZ heat elem size card requires 'Error ZZ heat = YES' card.  Please add.");
		    }
		  }
		  else if ( !strcmp( $3, "PRESSURE") )
		  {
		    if ( ERROR_ZZ_P == 1 )
		    {
		      common_error_element_size_processing( $7, $8, $9, $10, $11, $12 );
		    } 
		    else
		    {
		      declare_error("Error ZZ pressure elem size card requires 'Error ZZ pressure = YES' card.  Please add.");
		    }
		  }	    
		  card_read("ERROR_ZZ_ELEM_SIZE_CARD",file_index);		  
		}
		| error {}
		;
		
zz_error_type:    VELOCITY_ 	{strcpy($$,"VELOCITY");}
		| HEAT_ 	{strcpy($$,"HEAT");}
		| PRESSURE_ 	{strcpy($$,"PRESSURE");}
		| error {}
		;


post_processing_fluxes_section_delimiter_card:
		/* empty */ {} 
		| POST_ PROCESSING_ FLUXES_ EQUALS_ CR_
		{		
		  accept_ppfs_until_end_of_ppf_card = TRUE;
		  card_read("POST_PROCESSING_FLUXES_SECTION_DELIMITER_CARD",file_index);		  
		}
		| error {}
		;

flux_card:
		/* empty */ {} 
		| FLUX_ EQUALS_ valid_pp_flux_name integer integer integer  STRING_ optional_profile CR_
		{		 
		  common_pp_flux_processing( $3, pp_index, $4, $5, $6, $7, $8 );
		  card_read("FLUX_CARD",file_index);		  
		}
		| error {}
		;	

valid_pp_flux_name:
	  	  FORCE_NORMAL_		{pp_index=FORCE_NORMAL;}
	  	| FORCE_TANGENT1_	{pp_index=FORCE_TANGENT1;}
	  	| FORCE_TANGENT2_	{pp_index=FORCE_TANGENT2;}
	  	| FORCE_X_		{pp_index=FORCE_X;}
	  	| FORCE_Y_		{pp_index=FORCE_Y;}
	  	| FORCE_Z_		{pp_index=FORCE_Z;}
	  	| VOLUME_FLUX_		{pp_index=VOLUME_FLUX;}
	  	| SPECIES_FLUX_		{pp_index=SPECIES_FLUX;}
	  	| HEAT_FLUX_		{pp_index=HEAT_FLUX;}
	  	| TORQUE_		{pp_index=TORQUE;}
	  	| AVERAGE_CONC_		{pp_index=AVERAGE_CONC;}
	  	| SURF_DISSIP_		{pp_index=SURF_DISSIP;}
	  	| AREA_			{pp_index=AREA;}
	  	| VOL_REVOLUTION_	{pp_index=VOL_REVOLUTION;}
	  	| error {}
	  	;

optional_profile:
		  /* empty */ {$$=FALSE;}
		| PROFILE_ {$$=TRUE;}
		| error {}
		;

end_of_flux_card:
		/* empty */ {} 
		| END_ OF_ FLUX_ CR_
		{		
		  accept_ppfs_until_end_of_ppf_card = FALSE;
		  card_read("END_OF_FLUX_CARD",file_index);		  
		}
		| error {}
		;	

post_processing_fluxes_sensitivities_section_delimiter_card:
		/* empty */ {} 
		| POST_ PROCESSING_ FLUX_ SENSITIVITIES_ CR_
		{		
		  card_read("POST_PROCESSING_FLUX_SENSITIVITIES_SECTION_DELIMITER_CARD",file_index);		  
		}
		| error {}
		;

flux_sensitivity_card:
		/* empty */ {} 
		/*                                                        $6 species_number                             $9 sens_flt       $11 profile_flag     */
		/*$1         $2      $3flux_type_name   $4ss_id  $5blk_id          $7sens_type                 $8sens_id         $10flux_filenm           $12 */   
		| FLUX_SENS_ EQUALS_ valid_pp_flux_name integer integer integer valid_flux_sensitivity_type integer integer STRING_ optional_profile CR_
		{	
		  if (	 accept_ppfsens_until_end_of_ppfsens_card )
		  {
                    if ( nn_post_fluxes_sens >= 20)
                    {
                      declare_error("Parser is limited to twenty post processing flux sensitivity cards.");
                    }
                    else
                    {
                      if ( nn_post_fluxes_sens == 0)
                      {  
                        pp_fluxes_sens = (struct Post_Processing_Fluxes_Sens **) array_alloc(1, 20, sizeof(struct Post_Processing_Fluxes_Sens *));
                        for(i = 0; i < 20; i++)      		
                        {
                          pp_fluxes_sens[i] = (struct Post_Processing_Fluxes_Sens *) array_alloc(1, 1, sizeof(struct Post_Processing_Fluxes_Sens));
                        }
                      }
                      strcpy(pp_fluxes_sens[nn_post_fluxes_sens]->flux_type_name, $3);
                      pp_fluxes_sens[nn_post_fluxes_sens]->flux_type = pp_index;  
                      pp_fluxes_sens[nn_post_fluxes_sens]->ss_id = $4;
                      pp_fluxes_sens[nn_post_fluxes_sens]->blk_id = $5;
                      pp_fluxes_sens[nn_post_fluxes_sens]->species_number = $6;
                      pp_fluxes_sens[nn_post_fluxes_sens]->sens_type = $7;
                      pp_fluxes_sens[nn_post_fluxes_sens]->sens_id = $8;
                      pp_fluxes_sens[nn_post_fluxes_sens]->sens_flt = $9;                                                   
                      strcpy(pp_fluxes_sens[nn_post_fluxes_sens]->flux_filenm, $10);
                      pp_fluxes_sens[nn_post_fluxes_sens]->profile_flag = 11;
                      
                      pp_fluxes_sens[nn_post_fluxes_sens]->vector_id = -1;
                      for(j=0;j<nn_post_fluxes_sens;j++)
                      {
                       if(pp_fluxes_sens[nn_post_fluxes_sens]->sens_type == pp_fluxes_sens[j]->sens_type &&
                          pp_fluxes_sens[nn_post_fluxes_sens]->sens_id == pp_fluxes_sens[j]->sens_id &&
                          pp_fluxes_sens[nn_post_fluxes_sens]->sens_flt == pp_fluxes_sens[j]->sens_flt)
                       {pp_fluxes_sens[nn_post_fluxes_sens]->vector_id=pp_fluxes_sens[j]->vector_id;}
                      }
                      if(pp_fluxes_sens[nn_post_fluxes_sens]->vector_id == -1)
                      {
                        pp_fluxes_sens[nn_post_fluxes_sens]->vector_id=sens_vec_ct;
                        if(Continuation == ALC_FIRST)
	                {
                          if( (cont->upType == 1 && pp_fluxes_sens[nn_post_fluxes_sens]->sens_type == 1)
		           || (cont->upType == 3 && pp_fluxes_sens[nn_post_fluxes_sens]->sens_type == 3) )
	                  {
		            if(cont->upBCID == pp_fluxes_sens[nn_post_fluxes_sens]->sens_id &&
		               cont->upDFID == pp_fluxes_sens[nn_post_fluxes_sens]->sens_flt)
		            {cont->sensvec_id = sens_vec_ct;}
	                  }
	                  if(cont->upType == 2 && pp_fluxes_sens[nn_post_fluxes_sens]->sens_type == 2)
	                  {
		           if(cont->upMTID == pp_fluxes_sens[nn_post_fluxes_sens]->sens_id &&
		              cont->upMPID == pp_fluxes_sens[nn_post_fluxes_sens]->sens_flt)
		           {cont->sensvec_id = sens_vec_ct;}
	                  }
	                }
	                sens_vec_ct++;
                      }                       	    
		    }
		    nn_post_fluxes_sens++;	
		  }
		  card_read("FLUX_SENSITIVITY_CARD",file_index);		  
		}
		| error {}
		;	
		
valid_flux_sensitivity_type:
		  BC_	{$$=1;}
		| MT_	{$$=2;}
		| AC_	{$$=3;}
		| error {}
		;		


end_of_flux_sens_card:
		/* empty */ {} 
		| END_ OF_ FLUX_SENS_ CR_
		{		
		  accept_ppfsens_until_end_of_ppfsens_card = FALSE;
		  card_read("END_OF_FLUX_SENSITIVITY_CARD",file_index);		  
		}
		| error {}
		;

post_processing_data_specifications_delimiter_card:
		/* empty */ {} 
		|  POST_ PROCESSING_ DATA_ CR_
		{
		  card_read("POST_PROCESSING_DATA_DELIMITER_CARD",file_index);		  
		}
		| error {}
		;
		 
post_processing_data_card:
		/* empty */ {} 
		|  DATA_ EQUALS_ post_processing_data_var integer integer integer STRING_ CR_
		{		
 		  if ( accept_post_processing_data_until_end_of_ppdata_card )
 		  {
                    if ( nn_post_data >= 20)
                    {
                      declare_error("Parser is limited to twenty post processing data cards.");
                    }
                    else
                    {
                      if ( nn_post_data == 0)
                      {  
                        pp_data = (struct Post_Processing_Data **) array_alloc(1, 20, sizeof(struct Post_Processing_Data *));
                        for(i = 0; i < 20; i++)      		
                        {
                          pp_data[i] = (struct Post_Processing_Data *) array_alloc(1, 1, sizeof(struct Post_Processing_Data));
                        }
                      } 		  
 		    }
 		    pp_data[nn_post_data]->data_type = pp_index,
 		    strcpy( pp_data[nn_post_data]->data_type_name, $3);
 		    pp_data[nn_post_data]->ns_id = $4;
 		    pp_data[nn_post_data]->mat_id = $5;
 		    pp_data[nn_post_data]->mat_id--;
 		    pp_data[nn_post_data]->species_number = $6; 		    
 		    strcpy(pp_data[nn_post_data]->data_filenm, $7); 
 		    nn_post_data++;
 		  } /* end of accept_post.....*/ 
		  card_read("POST_PROCESSING_DATA_CARD",file_index);		  
		}
		| error {}
		;
		  
post_processing_data_var:
		  var_name {}
		| post_var_name {}
		| error {}
		;

end_of_data_card:
		/* empty */ {} 
		|  END_ OF_ DATA_ CR_
		{		
		  accept_post_processing_data_until_end_of_ppdata_card = FALSE;
		  card_read("END_OF_DATA_CARD",file_index);		  
		}
		| error {}
		;
		  
post_processing_data_sensitivity_delimiter_card:
		/* empty */ {} 
		|  POST_ PROCESSING_ DATA_ SENSITIVITIES_ CR_
		{		
		  card_read("POST_PROCESSING_DATA_SENSITIVITIES_DELIMITER_CARD",file_index);		  
		}
		| error {}
		;
		  
data_sens_card:
		/* empty */ {} 
		|  DATA_SENS_ EQUALS_ post_processing_data_var integer integer integer valid_flux_sensitivity_type integer integer STRING_ CR_
		{		
 		  if ( accept_data_sens_cards )
 		  {
                    if ( nn_post_data_sens >= 20)
                    {
                      declare_error("Parser is limited to twenty post processing data sensitivity cards.");
                    }
                    else
                    {
                      if ( nn_post_data_sens == 0)
                      {  
                        pp_data_sens = (struct Post_Processing_Data_Sens **) array_alloc(1, 20, sizeof(struct Post_Processing_Data_Sens *));
                        for(i = 0; i < 20; i++)      		
                        {
                          pp_data_sens[i] = (struct Post_Processing_Data_Sens *) array_alloc(1, 1, sizeof(struct Post_Processing_Data_Sens));
                        }
                      } 		  
 		    }
 		    pp_data_sens[nn_post_data_sens]->data_type = pp_index,
 		    strcpy( pp_data[nn_post_data]->data_type_name, $3);
 		    pp_data_sens[nn_post_data_sens]->ns_id = $4;
 		    pp_data_sens[nn_post_data_sens]->mat_id = $5;
 		    pp_data_sens[nn_post_data_sens]->mat_id--;
 		    pp_data_sens[nn_post_data_sens]->species_number = $6;
 		    pp_data_sens[nn_post_data_sens]->sens_type = $7;
 		    pp_data_sens[nn_post_data_sens]->sens_id = $8;
 		    pp_data_sens[nn_post_data_sens]->sens_flt = $9; 		    
 		    strcpy(pp_data_sens[nn_post_data_sens]->data_filenm, $10); 
                    pp_data_sens[nn_post_data_sens]->vector_id = -1;
 		    /* first search flux sensitivity info */
                    for(j=0;j<nn_post_fluxes_sens;j++)
                    {
                      if(pp_data_sens[nn_post_data_sens]->sens_type == pp_fluxes_sens[j]->sens_type &&
                      pp_data_sens[nn_post_data_sens]->sens_id == pp_fluxes_sens[j]->sens_id &&
                      pp_data_sens[nn_post_data_sens]->sens_flt == pp_fluxes_sens[j]->sens_flt)
                      {pp_data_sens[nn_post_data_sens]->vector_id=pp_fluxes_sens[j]->vector_id;}
                    }
                    for(j=0;j<nn_post_data_sens;j++)
                    {
                      if(pp_data_sens[nn_post_data_sens]->sens_type == pp_data_sens[j]->sens_type &&
                      pp_data_sens[nn_post_data_sens]->sens_id == pp_data_sens[j]->sens_id &&
                      pp_data_sens[nn_post_data_sens]->sens_flt == pp_data_sens[j]->sens_flt)
                      {pp_data_sens[nn_post_data_sens]->vector_id=pp_data_sens[j]->vector_id;}
                    }
                    if(pp_data_sens[nn_post_data_sens]->vector_id == -1)
                    {
                      pp_data_sens[nn_post_data_sens]->vector_id=sens_vec_ct;
                      /*
                      IF FIRST ORDER CONTINUATION CHECK FOR SAME SENSITIVITY PARAMETER
                      */
                      if(Continuation == ALC_FIRST)
                      {
                        if( (cont->upType == 1 && pp_data_sens[nn_post_data_sens]->sens_type == 1)
                        || (cont->upType == 3 && pp_data_sens[nn_post_data_sens]->sens_type == 3) )
                        {
                          if(cont->upBCID == pp_data_sens[nn_post_data_sens]->sens_id &&
                          cont->upDFID == pp_data_sens[nn_post_data_sens]->sens_flt)
                            {cont->sensvec_id = sens_vec_ct;}
                        }
                        if(cont->upType == 2 && pp_data_sens[nn_post_data_sens]->sens_type == 2)
                        {
                          if(cont->upMTID == pp_data_sens[nn_post_data_sens]->sens_id &&
                          cont->upMPID == pp_fluxes_sens[nn_post_data_sens]->sens_flt)
                            {cont->sensvec_id = sens_vec_ct;}
                        }
                     }
                     sens_vec_ct++; 	
                     nn_post_data_sens++;
                   }	    
 		  } /* end of accept_post.....*/ 		
		  card_read("DATA_SENS_CARD",file_index);		  
		}
		| error {}
		;
		
data_sens_type:
		  AC_	{$$=3;}
		| BC_	{$$=1;}
		| MT_	{$$=2;}
		| error {}
		;

end_of_data_sens_card:
		/* empty */ {} 
		|  END_ OF_ DATA_SENS_ CR_
		{		
		  accept_data_sens_cards = FALSE;
		  card_read("END_OF_DATA_SENS_CARD",file_index);		  
		}
		| error {}
		;
		  
post_processing_particle_tracking_calculations_delimiter_card:
		/* empty */ {} 
		|  POST_ PROCESSING_ PARTICLE_ TRACES_ CR_
		{		
		
		  card_read("POST_PROCESSING_PARTICLE_TRACKING_CALCULATIONS_DELIMITER_CARD",file_index);		  
		}
		| error {}
		;
		  
particle_card:
		/* empty */ {} 
		|  PARTICLE_ EQUALS_ float float float float float float STRING_ CR_
		{		
		  if ( accept_post_processing_particles )
 		  {
                    if ( nn_particles >= 20)
                    {
                      declare_error("Parser is limited to twenty post processing particle cards.");
                    }
                    else
                    {
                      if ( nn_particles == 0)
                      {  
                        pp_particles = (struct Post_Processing_Particles **) array_alloc(1, 20, sizeof(struct Post_Processing_Particles *));
                        for(i = 0; i < 20; i++)      		
                        {
                          pp_particles[i] = (struct Post_Processing_Particles *) array_alloc(1, 1, sizeof(struct Post_Processing_Particles));
                        }
                      } 		  
 		    }
 		    pp_particles[nn_particles]->coord[0] = $3;
 		    pp_particles[nn_particles]->coord[1] = $4;
 		    pp_particles[nn_particles]->coord[2] = $5;
 		    pp_particles[nn_particles]->Start_Time = $6;
 		    pp_particles[nn_particles]->End_Time = $7;
 		    pp_particles[nn_particles]->Delta_s = $8;
 		    strcpy(pp_particles[nn_particles]->filenm, $9); 
 		    nn_particles++;
 		  } /* end of accept_post.....*/ 		
		  card_read("PARTICLE_CARD",file_index);		  
		}
		| error {}
		;
		  
end_of_particles_card:
		/* empty */ {} 
		|  END_ OF_ PARTICLES_ CR_
		{		
		  accept_post_processing_particles = FALSE;
		  card_read("END_OF_PARTICLES_CARD",file_index);		  
		}
		| error {}
		;
		  
post_processing_volumetric_integration_delimiter_card:
		/* empty jjj ???*/ {} 
		|  POST_ PROCESSING_ VOLUMETRIC_ INTEGRATION_ CR_
		{		
		
		  card_read("POST_PROCESSING_VOLUMETRIC_INTEGRATION_DELIMITER_CARD",file_index);		  
		}
		| error {}
		;
		  
volume_int_card:
		/* empty */ {} 
		|  VOLUME_INT_ EQUALS_ CR_	/* jjj */
		{		
		
		  card_read("VOLUME_INT_CARD",file_index);		  
		}
		| error {}
		;

pp_vol_names: 
    VOLUME_		{pp_index = I_VOLUME ;}	/* jjj not used anywhere */
  | DISSIPATION_	{pp_index = I_DISSIP ;}
  | SPECIES_MASS_	{pp_index = I_SPECIES_MASS;}
  | HEAT_ENERGY_	{pp_index = I_HEAT_ENERGY ;}
  | MOMENTUMX_		{pp_index = I_MOMX ;}
  | MOMENTUMY_		{pp_index = I_MOMY ;}
  | MOMENTUMZ_		{pp_index = I_MOMZ ;}
  | STRESS_TRACE_	{pp_index = I_TRACE ;}
  | POSITIVE_FILL_	{pp_index = I_POS_FILL;}
  | NEGATIVE_FILL_	{pp_index = I_NEG_FILL;}
  | error {}
  ;

		  
end_of_volume_int_card:
		/* empty */ {} 
		|  END_ OF_ VOLUME_INT_ CR_
		{		
		
		  card_read("END_OF_VOLUME_INT_CARD",file_index);		  
		}
		| error {}
		;
		  

yes_no_option:
		YES_ {$$ = 1;} | NO_ {$$ = 0;} | ON_ {$$ = 1;} | OFF_ {$$ = 0;}
		
float:
	  FLOAT_ {$$=atof($1)}
	| INTEGER_ 
	  {
	   /*  sprintf(msg, "Integer %s found in place of floating point number. );
	    declare_error(msg);*/
	  }
	| error {}
	;
		
/*********************************** MAT File Specifications rules:***************************************/
mat_file:	mat_card 
                  {  floating_point_constant_list_index = -1; 
                     line_number++;
                     if (ProcID == 0)  fprintf(parser_log,"\n%d ",line_number );
                  } 
		| mat_file mat_card 
		  {  floating_point_constant_list_index = -1; 
                     line_number++;
                     if (ProcID == 0)  fprintf(parser_log,"\n%d ",line_number );
                  }
		| error {}
		;

mat_card:	/* empty */ {}
		| default_database_card {}
	 	| density_card {}
	 	| solid_constitutive_equation_card {}
	 	| convective_lagrangian_velocity_card {}
	 	| lame_mu_card {}
	 	| lame_lambda_card {}
	 	| conductivity_card {}
	 	| heat_capacity_card {}
	 	| reference_temperature_card {}
	 	| diffusivity_card {}
	 	| plastic_viscosity_card {}
	 	| latent_heat_vaporization_card {}
	 	| latent_heat_fusion_card {}
	 	| species_volume_expansion_card {}
	 	| reference_concentration_card {}
	 	| vapor_pressure_card {} {}	 	
	 	| porous_latent_heat_vaporization_card {}
	 	| volume_expansion_card {}
	 	| electrical_conductivity_card {}
	 	| media_type_card {}
	 	| porosity_card {}
	 	| permeability_card {}
	 	| flowingliquid_viscosity_card {}
	 	| inertia_coefficient_card {}
	 	| porous_diffusion_constitutive_equation_card {}
	 	| diffusion_constitutive_equation_card {}
	 	| liquidus_temperature_card {}
	 	| solidus_temperature_card {}
	 	| stress_free_solvent_vol_frac_card {}
		| solid_thermal_expansion_card {}
	 	| liquid_constitutive_card {}
	 	| mechanical_properties_viscosity_card {}
	 	| low_rate_viscosity_card {}
	 	| power_law_exponent_equation_card {}	   
	 	| high_rate_viscosity_card {}
	 	| time_constant_card {}
	 	| aexp_card {}
	 	| thermal_exponent_card {}
	 	| thermal_wlf_constant2_card {}
	 	| yield_stress_card {}
	 	| yield_exponent_card {}
	 	| suspension_maximum_packing_card {}
	 	| suspension_species_number_card {}
	 	| cure_gel_point_card {}
	 	| cure_a_exponent_card {}
	 	| cure_b_exponent_card {}
	 	| cure_species_number_card {}
	 	| cure_species_number_card {}
	 	| polymer_constitutive_equation_card {}
	 	| polymer_stress_formulation_card {}
	 	| polymer_weight_function_card {}
	 	| polymer_weighting_card {}
	 	| polymer_viscosity_card {}
	 	| polymer_time_constant_card {}
	 	| mobility_parameter_card {}
	 	| surface_tension_card {}
	 	| plasticity_equation_card {}
	 	| pseudo_solid_constitutive_equation_card {}
	 	| pseudo_solid_lame_mu_card {}
	 	| pseudo_solid_lame_lambda_card {}
	 	| solution_temperature_card {}
	 	| non_condensable_molecular_weight_card {}
	 	| non_volatile_molar_volume_card {}
	 	| non_volatile_specific_volume_card {}
	 	| flory_huggins_card {}
	 	| flory_huggins_parameter_card {}
	 	| mass_source_card {}
	 	| heat_source_card {}
	 	| solid_body_source_card {}	 	
	 	| species_source_card {}
	 	| current_source_card {}
	 	| material_variable_initialize_card {}
	 	| navier_stokes_source_card {}	 	
	 	| table_data_card  {}
	 	| end_table_card  {}
	 	| file_equals_card   {}
	 	| example_card {}	/*ROLL-UP OF EXAMPLE CARD INOT MATERIAL CARDS*/
		| error  {}
		;
		
default_database_card:
		/* empty */ {}  
		| DEFAULT_ DATABASE_ EQUALS_  database_type CR_
		  {
		    if($4==DB_CHEMKIN_MAT)
		    {
                      sprintf(msg, "CHEMKIN is currently not supported by the lex/yacc parser.  Material %s defaulting to Goma default database.", pd_glob[mn]->MaterialName);
                      declare_warning(msg); 		      
		    }
                    if (goma_mat_prop_init(mp_glob[mn], mn, pd_glob[mn]) < 0) 
                    {
                      sprintf(msg,"%s for mat %d had an inconsistency, BAIL!", "goma_mat_prop_int", mn);
                      declare_error(msg);
                    }		    
		    card_read("DEFAULT_DATABASE_CARD",file_index);
		  }
		| error {}
		;
		
database_type:
		  CHEMKIN_MAT_	{$$=DB_CHEMKIN_MAT;}
		| GOMA_MAT_	{$$= DB_GOMA_MAT;}
		| error {}
		;	
				
density_card:   /* empty */ {} 
		| DENSITY_ EQUALS_ density_model optional_floating_point_constant_list CR_
                {
                  if ($3 == DENSITY_FILL && number_of_constants < 2)
                  {
                    sprintf(msg, "Material %s - expected at least 2 constants for the Density FILL model.", pd_glob[mn]->MaterialName);
                    declare_error(msg);                 
                  }
                  else if ($3 == DENSITY_SUSPENSION && number_of_constants < 3) 
                  {
                    sprintf(msg, "Material %s - expected at least 3 constants for the Density SUSPENSION model.", pd_glob[mn]->MaterialName);
                    declare_error(msg);                   
                  }
                  else if ($3 == DENSITY_SUSPENSION_PM && number_of_constants < 3) 
                  {
                    sprintf(msg, "Material %s - expected at least 3 constants for the Density SUSPENSION_PM model.", pd_glob[mn]->MaterialName);
                    declare_error(msg);                      
                  }
                  else if ($3 == DENSITY_THERMAL_BATTERY && number_of_constants < 2) 
                  {
                    sprintf(msg, "Material %s - expected at least 2 constants for the Density THERMAL_BATTERY model.", pd_glob[mn]->MaterialName);
                    declare_error(msg);                      
                  }
                  else if ($3 == DENSITY_LEVEL_SET && number_of_constants < 3) 
                  {
                    sprintf(msg, "Material %s - expected at least 3 constants for the Density LEVEL_SET model.", pd_glob[mn]->MaterialName);
                    declare_error(msg);                      
                  }
                  else if ($3 == DENSITY_FOAM && number_of_constants < 2) 
                  {
                    sprintf(msg, "Material %s - expected at least 2 constants for the Density FOAM model.", pd_glob[mn]->MaterialName);
                    declare_error(msg);                      
                  }
                  else if ($3 == DENSITY_CONSTANT_LAST_CONC && number_of_constants < 1)
                  {
                    sprintf(msg, "Material %s - expected at least 1 constant for the Density CONSTANT_LAST_CONC model.", pd_glob[mn]->MaterialName);
                    declare_error(msg);                      
                  }
                  else if ($3 == DENSITY_IDEAL_GAS)
                  {
                    if (number_of_constants > 0)
                    {
                      sprintf(msg, "Material %s - expected no constants for the Density IDEAL_GAS model.  Constants being ignored.", pd_glob[mn]->MaterialName);
                      declare_warning(msg);                          
                    }
                    mp_glob[mn]->DensityModel = $3;   
                    card_read("DENSITY_CARD",file_index);      
                  }
                  else if ( ( ($3 == USER) || ($3 == USER_GEN) ) && number_of_constants == 0)
                  {
                    sprintf(msg, "Material %s - expected at least 1 constant for the Density USER or USER_GEN model.", pd_glob[mn]->MaterialName);
                    declare_error(msg);                       
                  }
                  else if ($3 == CONSTANT)
                  {
		    switch(number_of_constants)
		    {
		      case 0:
                        sprintf(msg, "Material %s - expected 1 constant for the Density CONSTANT model.", pd_glob[mn]->MaterialName);
                        declare_error(msg); 		          
		        break;
		      case 1:
         	        mp_glob[mn]->DensityModel = CONSTANT;   
                        mp_glob[mn]->density = floating_point_constant_list_array[0]; 
                        card_read("DENSITY_CARD",file_index);   
		        break;
		      default:
                        sprintf(msg, "Material %s - expected only one constant for the Density CONSTANT model.  Subsequent constants being ignored.", pd_glob[mn]->MaterialName);
                        declare_warning(msg);   
            	        mp_glob[mn]->DensityModel = CONSTANT;   
                        mp_glob[mn]->density = floating_point_constant_list_array[0]; 
                        card_read("DENSITY_CARD",file_index);	  		        
		    } 	              
                  }
                  else
                  {
                    mp_glob[mn]->DensityModel = $3;   
                    mp_glob[mn]->len_u_density = read_floats( &(mp_glob[mn]->u_density), 0); 
                    card_read("DENSITY_CARD",file_index);                                   
                  } 
                } 
		| error {}
		;
		
density_model:
		  FILL_			{$$=DENSITY_FILL;}
                | SUSPENSION_		{$$=DENSITY_SUSPENSION;}
                | SUSPENSION_PM_	{$$=DENSITY_SUSPENSION_PM;}
                | THERMAL_BATTERY_	{$$=DENSITY_THERMAL_BATTERY;}
                | LEVEL_SET_		{$$=DENSITY_LEVEL_SET;}
                | FOAM_			{$$=DENSITY_FOAM;}
                | CONSTANT_LAST_CONC_	{$$=DENSITY_CONSTANT_LAST_CONC;}
                | IDEAL_GAS_		{$$=DENSITY_IDEAL_GAS;}
                | CONSTANT_		{$$=CONSTANT;}
                | USER_			{$$=USER;}
                | USER_GEN_		{$$=USER_GEN;}
                | error {}
                ;
		
solid_constitutive_equation_card:	
		/* empty */ {}  
		| SOLID_ CONSTITUTIVE_ EQUATION_ EQUALS_ solid_constitutive_equation_model CR_
		{
		  cr_glob[mn]->MeshFluxModel = $5;		  
		  card_read("SOLID_CONSTITUTIVE_EQUATION_CARD",file_index);
		}  
		| error {}
		;

solid_constitutive_equation_model:
		  LINEAR_	{$$ = LINEAR;}
		| NONLINEAR_	{$$ = NONLINEAR;}
		| NONLINEAR_PLANE_STRAIN_	{$$ = NONLINEAR;}
		| INCOMP_PSTRAIN_	{$$ = INCOMP_PSTRAIN;}
		| INCOMP_PSTRESS_	{$$ = INCOMP_PSTRESS;}
		| HOOKEAN_PSTRAIN_	{$$ = HOOKEAN_PSTRAIN;}
		| HOOKEAN_PSTRESS_	{$$ = HOOKEAN_PSTRESS;}
		| INCOMP_3D_		{$$ = INCOMP_3D;}
		| error {}
		;

plasticity_equation_card:	
		/* empty */ {}  
		| PLASTICITY_ EQUATION_ EQUALS_ plasticity_equation_model CR_
		{
		  evpl_glob[mn]->ConstitutiveEquation = $4;		  
		  card_read("PLASTICITY_EQUATION_CARD",file_index);
		}		
		| error {}
		;
		
plasticity_equation_model:
		  EVP_HYPER_	{$$=EVP_HYPER;}
		| NO_MODEL_	{$$=NO_MODEL;}
		| error {}
		;

convective_lagrangian_velocity_card:	
		/* empty */ {} 
		| CONVECTIVE_ LAGRANGIAN_ VELOCITY_ EQUALS_ convective_lagrangian_velocity_model optional_floating_point_constant_list CR_ 
		{
		  if (pd_glob[mn]->MeshMotion == LAGRANGIAN || pd_glob[mn]->MeshMotion == TOTAL_ALE)	
		  {
		    if ($5 == CONSTANT)
		    {
		      if ( number_of_constants >= 3 )
		      {
	  	        elc_glob[mn]->v_mesh_sfs_model = pd_glob[mn]->MeshInertia = CONSTANT;  
	  	        elc_glob[mn]->v_mesh_sfs[0] = floating_point_constant_list_array[0];
	  	        elc_glob[mn]->v_mesh_sfs[1] = floating_point_constant_list_array[1];
	  	        elc_glob[mn]->v_mesh_sfs[2] = floating_point_constant_list_array[2];
	  	        if (pd_glob[mn]->TimeIntegration == TRANSIENT)
	  	        {
	      		  sprintf(msg,"Currently can't do transient mesh motion with inertia.");
	      		  declare_error(msg);
	      	        }
	      	        if ( number_of_constants > 3)
	      	        {
	      		  sprintf(msg,"Only three floats required for CONSTANT model.  Subsequent floats are being ignored.");
	      		  declare_warning(msg);	      	        
	      	        }
	      	      }
	      	      else
	      	      {
	      		  sprintf(msg,"Expected three floats for CONSTANT model.");
	      		  declare_error(msg);	      	      
	      	      }		    
		    }
		    else if ($5 == ROTATIONAL)
		    { 
		      if (number_of_constants >= 4)
		      {
	  	        pd_glob[mn]->MeshInertia = 1;
	  	        elc_glob[mn]->v_mesh_sfs_model = ROTATIONAL;
	                elc_glob[mn]->len_u_v_mesh_sfs = read_floats( &(elc_glob[mn]->u_v_mesh_sfs), 0);	 	  	        		      
		      }
		      else
		      {
		        sprintf(msg, 
	          		"Matl %s expected at least 4 constants for Convective Lagrangian Velocity ROTATIONAL model.",
			       	pd_glob[mn]->MaterialName);
		  	declare_error(msg);
		      }
		    }
		    else /* model is NONE */
		    {	  
             	      pd_glob[mn]->MeshInertia = 0;
                      elc_glob[mn]->v_mesh_sfs_model = 0;
                      elc_glob[mn]->v_mesh_sfs[0] = 0.;
                      elc_glob[mn]->v_mesh_sfs[1] = 0.;
                      elc_glob[mn]->v_mesh_sfs[2] = 0.;
		      if (number_of_constants > 0)
		      {
	      	        sprintf(msg,"No floats required for NONE.  Floats are being ignored.");
	      	        declare_warning(msg);	  			
		      }		      
		    }
		    card_read("CONVECTIVE_LAGRANGIAN_VELOCITY_CARD",file_index);
		  }
		  else
		  {
	            elc_glob[mn]->v_mesh_sfs_model = pd_glob[mn]->MeshInertia;
	            if (pd_glob[mn]->TimeIntegration == TRANSIENT)
	            {
	      	        sprintf(msg,"Currently can't do transient mesh motion with inertia");
	      	        declare_error(msg);	/* possibly should be a warning */
	            }		  
		  }
		}				
		| error {}
		;

convective_lagrangian_velocity_model:
		  CONSTANT_	{$$=CONSTANT;}
		| ROTATIONAL_	{$$=ROTATIONAL;}
		| NONE_		{$$=0;}
		| error {}
		;
		
lame_mu_card:
		/* empty */ {} 
		| LAME_ MU_ EQUALS_ lame_mu_model_name optional_floating_point_constant_list CR_ 
		{
		  if ($4 == POWER_LAW && number_of_constants < 3 )
		  {
                      sprintf(msg, "Material %s - expected at least 3 constants for the Lame MU  POWER_LAW model.", pd_glob[mn]->MaterialName);
                      declare_error(msg);               		  
		  }
		  else if ($4 == CONTACT_LINE && number_of_constants < 4 )
		  {
                      sprintf(msg, "Material %s - expected at least 4 constants for the Lame MU CONTACT_LINE model.", pd_glob[mn]->MaterialName);
                      declare_error(msg); 		  
		  } 
		  else if ($4 == SHEAR_HARDEN && number_of_constants < 2 )
		  {
                      sprintf(msg, "Material %s - expected at least 2 constants for the Lame MU SHEAR_HARDEN model.", pd_glob[mn]->MaterialName);
                      declare_error(msg); 		  
		  }
		  else if ($4 == EXPONENTIAL && number_of_constants < 3 )
		  {
                      sprintf(msg, "Material %s - expected at least 3 constants for the Lame MU EXPONENTIAL model.", pd_glob[mn]->MaterialName);
                      declare_error(msg); 		  
		  }
		  else if ($4 == DENSE_POWER_LAW && number_of_constants < 2 )
		  {
                      sprintf(msg, "Material %s - expected at least 2 constants for the Lame MU DENSE_POWER_LAW model.", pd_glob[mn]->MaterialName);
                      declare_error(msg); 		  
		  }
                  else if ( ( ($4 == USER) || ($4 == USER_GEN) ) && number_of_constants == 0)
                  {
                    sprintf(msg, "Material %s - expected at least 1 constant for the Lame MU USER or USER_GEN model.", pd_glob[mn]->MaterialName);
                    declare_error(msg);                       
                  }		  
  		  else if($4 == CONSTANT )
		  {
		    switch(number_of_constants)
		    {
		      case 0:
                        sprintf(msg, "Material %s - expected one constant for the Lame MU CONSTANT model.", pd_glob[mn]->MaterialName);
                        declare_error(msg); 		          
		        break;
		      case 1:
		        elc_glob[mn]->len_u_mu = read_floats( &(elc_glob[mn]->u_mu), 0);
		        elc_glob[mn]->lame_mu_model = $4;
		        card_read("LAME_MU_CARD",file_index);
		        break;
		      default:
                        sprintf(msg, "Material %s - expected only one constant for the Lame MU CONSTANT model.  Subsequent constants being ignored.", pd_glob[mn]->MaterialName);
                        declare_warning(msg);   
		        elc_glob[mn]->len_u_mu = floating_point_constant_list_array[0];
		        elc_glob[mn]->lame_mu_model = $4;
		        card_read("LAME_MU_CARD",file_index);	  		        
		    }	         		  
		  }                
		  else
		  {
		    elc_glob[mn]->len_u_mu = read_floats( &(elc_glob[mn]->u_mu), 0);
		    elc_glob[mn]->lame_mu_model = $4;
		    card_read("LAME_MU_CARD",file_index);		  
		  }
		}
		| error {}
		;
		
lame_mu_model_name:
		  CONSTANT_ 		{$$=CONSTANT;}
		| POWER_LAW_ 		{$$=POWER_LAW;}
		| CONTACT_LINE_ 	{$$=CONTACT_LINE;}
		| USER_ 		{$$=USER;}
		| SHEAR_HARDEN_ 	{$$=SHEAR_HARDEN;}
		| EXPONENTIAL_ 		{$$=EXPONENTIAL;}
		| DENSE_POWER_LAW_	{$$=DENSE_POWER_LAW;}
		| USER_GEN_ 		{$$=USER_GEN;}
		| error {}
		;

lame_lambda_card:
		/* empty */ {} 
		| LAME_ LAMBDA_ EQUALS_ lame_lambda_model_name optional_floating_point_constant_list CR_ 
		{
		  if ($4 == POWER_LAW && number_of_constants < 3 )
		  {
                      sprintf(msg, "Material %s - expected at least 3 constants for the Lame LAMBDA  POWER_LAW model.", pd_glob[mn]->MaterialName);
                      declare_error(msg);               		  
		  }
		  else if ($4 == EXPONENTIAL && number_of_constants < 3 )
		  {
                      sprintf(msg, "Material %s - expected at least 3 constants for the Lame LAMBDA EXPONENTIAL model.", pd_glob[mn]->MaterialName);
                      declare_error(msg); 		  
		  }
                  else if ( ( ($4 == USER) || ($4 == USER_GEN) ) && number_of_constants == 0)
                  {
                    sprintf(msg, "Material %s - expected at least 1 constant for the Lame LAMBDA USER or USER_GEN model.", pd_glob[mn]->MaterialName);
                    declare_error(msg);                       
                  }		  
  		  else if($4 == CONSTANT )
  		  {
		    switch(number_of_constants)
		    {
		      case 0:
                        sprintf(msg, "Material %s - expected one constant for the Lame LAMBDA CONSTANT model.", pd_glob[mn]->MaterialName);
                        declare_error(msg);      		          
		        break;
		      case 1:
		        elc_glob[mn]->len_u_lambda = floating_point_constant_list_array[0];
		        elc_glob[mn]->lame_lambda_model = $4;
		        card_read("LAME_LAMBDA_CARD",file_index);	 
		        break;
		      default:
                        sprintf(msg, "Material %s - expected only one constant for the Lame LAMBDA CONSTANT model.  Subsequent constants being ignored.", pd_glob[mn]->MaterialName);
                        declare_warning(msg);   
		        elc_glob[mn]->len_u_lambda = floating_point_constant_list_array[0];
		        elc_glob[mn]->lame_lambda_model = $4;
		        card_read("LAME_LAMBDA_CARD",file_index);			        
		    }      		  
		  }                
		  else
		  {
		    elc_glob[mn]->len_u_lambda = read_floats( &(elc_glob[mn]->u_lambda), 0);
		    elc_glob[mn]->lame_lambda_model = $4;
		    card_read("LAME_LAMBDA_CARD",file_index);		  
		  }
		}		
		| error {}
		;
	
conductivity_card:
		/* empty */ {}
		| CONDUCTIVITY_ EQUALS_ conductivity_model optional_floating_point_constant_list CR_
		{
 		  if($3 == CONSTANT )
		  {
		    switch(number_of_constants)
		    {
		      case 0:
                        sprintf(msg, "Material %s - expected one constant for the Conductivity CONSTANT model.", pd_glob[mn]->MaterialName);
                        declare_error(msg);     		          
		        break;
		      case 1:
		        mat_ptr->thermal_conductivity = floating_point_constant_list_array[0];
		        mat_ptr->ConductivityModel = $3;
		        
		        card_read("CONDUCTIVITY_CARD",file_index);	 
		        break;
		      default:
                        sprintf(msg, "Material %s - expected only one constant for the Plastic Viscosity CONSTANT model.  Subsequent constants being ignored.", pd_glob[mn]->MaterialName);
                        declare_warning(msg);   
		        mat_ptr->thermal_conductivity = floating_point_constant_list_array[0];
		        mat_ptr->ConductivityModel = $3;
		        card_read("CONDUCTIVITY_CARD",file_index);		        
		    }		      		  
		  }  
		  else
		  {
		    if( &(mat_ptr->u_thermal_conductivity) == NULL) 
	              {
	                mat_ptr->u_thermal_conductivity = (dbl *)array_alloc(1, number_of_constants, sizeof(dbl));
	              }	
		    mat_ptr->ConductivityModel = $3;
                    mat_ptr->len_u_thermal_conductivity = read_floats( &(mat_ptr->u_thermal_conductivity) , 0);
                    card_read("CONDUCTIVITY_CARD",file_index);	
                  }	     		
		}
		| CONDUCTIVITY_ EQUALS_ TABLE_ integer mp_table_independent_variable_name  		/* $1 $2 $3 $4 $5 */
		                                       optional_mp_table_independent_variable_name 	/* $6 */
		                                       optional_mp_table_independent_variable_name 	/* $7 */
		                                       mp_table_interpolation_method 			/* $8 */
		                                       optional_file_equals_string 			/* $9 */
		                                       optional_name_equals_string			/* $10 */
		                                       CR_
		{
                  if( num_MP_Tables == MAX_MP_TABLES )
	          {
	            sprintf(msg, "Material %s - Maximum TABLE_MPs exceeded.", pd_glob[mn]->MaterialName);
                    declare_error(msg);
	          }
	          else
                  {
		    table_type = CONDUCTIVITY_TABLE;
		    MP_Tables[num_MP_Tables] = setup_mp_table(mat_ptr->table,$4,"Conductivity",$5,$6,$7,$8,mp_table_species_number);
                    if ( (strcmp( $9, "NA") == 0) )
                    { 
                      /* Don't read table data from another file */
                      accept_table_data_cards=TRUE;
                    }
                    else
                    {
                      /* Try to read table data from another file */ 
                      accept_table_data_cards=FALSE;
		      table_type = CONDUCTIVITY_TABLE;
		      parse_table_file($9,$10,$4);
    	              read_array_into_table(table_type, mat_ptr->table, number_of_table_cards_read);
                    }		    
		  }
		  card_read("CONDUCTIVITY_CARD",file_index);
		}
		| error {}
		;	
		
conductivity_model:
		  USER_		{$$=USER;}
		| CONSTANT_	{$$=CONSTANT;}
		| error		{}
		;


heat_capacity_card:
		/* empty */ {}
		| HEAT_ CAPACITY_ EQUALS_ conductivity_model optional_floating_point_constant_list CR_
		{
 		  if($4 == CONSTANT )
		  {
		    switch(number_of_constants)
		    {
		      case 0:
                        sprintf(msg, "Material %s - expected one constant for the Conductivity CONSTANT model.", pd_glob[mn]->MaterialName);
                        declare_error(msg);     		          
		        break;
		      case 1:
		        mat_ptr->heat_capacity = floating_point_constant_list_array[0];
		        mat_ptr->HeatCapacityModel = $4;
		        card_read("HEAT_CAPACITY_CARD",file_index);	 
		        break;
		      default:
                        sprintf(msg, "Material %s - expected only one constant for the Plastic Viscosity CONSTANT model.  Subsequent constants being ignored.", pd_glob[mn]->MaterialName);
                        declare_warning(msg);   
		        mat_ptr->heat_capacity = floating_point_constant_list_array[0];
		        mat_ptr->HeatCapacityModel = $4;
		        card_read("HEAT_CAPACITY_CARD",file_index);		        
		    }		      		  
		  }  
		  else if ($4 == ENTHALPY )
		  {
		    switch(number_of_constants)
		    {
		      case 0:
		      case 1:
                        sprintf(msg, "Material %s - expected two constants for the Heat Capacity ENTHALPY model.", pd_glob[mn]->MaterialName);
                        declare_error(msg);     		          
		        break;
		      case 2:
		        mat_ptr->HeatCapacityModel = $4;
                        mat_ptr->len_u_heat_capacity = read_floats( &(mat_ptr->heat_capacity) , 0);	        
		        card_read("HEAT_CAPACITY_CARD",file_index);	 
		        break;
		      default:
                        sprintf(msg, "Material %s - expected only two constant for the Heat Capacity ENTHALPY model.  Subsequent constants being ignored.", pd_glob[mn]->MaterialName);
                        declare_warning(msg);   
		        mat_ptr->len_u_heat_capacity = read_floats( &(mat_ptr->heat_capacity) , 0);
		        mat_ptr->HeatCapacityModel = $4;
		        card_read("HEAT_CAPACITY_CARD",file_index);		        
		    }		      		  		  
		  }
		  else	/* USER */
		  {
		    mat_ptr->HeatCapacityModel = $4; /* USER */
                    mat_ptr->len_u_heat_capacity = read_floats( &(mat_ptr->heat_capacity) , 0);
                    card_read("HEAT_CAPACITY_CARD",file_index);	
                  }	     		
		}
		| HEAT_ CAPACITY_ EQUALS_ TABLE_ integer mp_table_independent_variable_name  		/* $1 $2 $3 $4 $5 $6*/
		                                       optional_mp_table_independent_variable_name 	/* $7 */
		                                       optional_mp_table_independent_variable_name 	/* $8 */
		                                       mp_table_interpolation_method 			/* $9 */
		                                       optional_file_equals_string 			/* $10 */
		                                       optional_name_equals_string			/* $11 */
		                                       CR_
		{
                  if( num_MP_Tables == MAX_MP_TABLES )
	          {
	            sprintf(msg, "Material %s - Maximum TABLE_MPs exceeded.", pd_glob[mn]->MaterialName);
                    declare_error(msg);
	          }
	          else
                  {
		    table_type = HEAT_CAPACITY_TABLE;
		    MP_Tables[num_MP_Tables] = setup_mp_table(mat_ptr->table,$5,"Heat Capacity",$6,$7,$8,$9,mp_table_species_number);
                    if ( (strcmp( $10, "NA") == 0) )
                    { 
                      /* Don't read table data from another file */
                      accept_table_data_cards=TRUE;
                    }
                    else
                    {
                      /* Try to read table data from another file */ 
                      accept_table_data_cards=FALSE;
		      table_type = HEAT_CAPACITY_TABLE;
		      parse_table_file($10,$11,$5);
    	              read_array_into_table(table_type, mat_ptr->table, number_of_table_cards_read);
                    }		    
		  }
                  card_read("HEAT_CAPACITY_CARD",file_index);
		}
		| error {}
		;	
		
heat_capacity_model:
		  USER_		{$$=USER;}
		| CONSTANT_	{$$=CONSTANT;}
		| ENTHALPY_	{$$=ENTHALPY;}
		| error		{}
		;
		
volume_expansion_card:
		/* empty */ {} 
		| VOLUME_ EXPANSION_ EQUALS_ CONSTANT_ optional_floating_point_constant_list CR_ 
		{
		  switch(number_of_constants)
		  {
		    case 0:
                      sprintf(msg, "Material %s - expected one constant for the Volume Expansion Temperature CONSTANT model.", pd_glob[mn]->MaterialName);
                      declare_error(msg);     		          
		      break;
		    case 1:
		      mat_ptr->Volume_Expansion = floating_point_constant_list_array[0];
		      mat_ptr->VolumeExpansionModel = CONSTANT;
		      card_read("VOLUME_EXPANSION_CARD",file_index);	 
		      break;
		    default:
                      sprintf(msg, "Material %s - expected only one constant for the Volume Expansion Temperature CONSTANT model.  Subsequent constants being ignored.", pd_glob[mn]->MaterialName);
                      declare_warning(msg);   
		      mat_ptr->Volume_Expansion = floating_point_constant_list_array[0];
		      mat_ptr->VolumeExpansionModel = CONSTANT;
		      card_read("VOLUME_EXPANSION_CARD",file_index);			        
		  }    		      		  
		}
		| error {}
		;			


liquidus_temperature_card:
		/* empty */ {} 
		| LIQUIDUS_ TEMPERATURE_ EQUALS_ CONSTANT_ optional_floating_point_constant_list CR_ 
		{
		  switch(number_of_constants)
		  {
		    case 0:
                      sprintf(msg, "Material %s - expected one constant for the Liquidus Temperature CONSTANT model.", pd_glob[mn]->MaterialName);
                      declare_error(msg);     		          
		      break;
		    case 1:
		      mat_ptr->melting_point_liquidus = floating_point_constant_list_array[0];
		      mat_ptr->LiquidusModel = CONSTANT;
		      card_read("LIQUIDUS_TEMPERATURE_CARD",file_index);	 
		      break;
		    default:
                      sprintf(msg, "Material %s - expected only one constant for the Liquidus Temperature CONSTANT model.  Subsequent constants being ignored.", pd_glob[mn]->MaterialName);
                      declare_warning(msg);   
		      mat_ptr->melting_point_liquidus = floating_point_constant_list_array[0];
		      mat_ptr->LiquidusModel = CONSTANT;
		      card_read("LIQUIDUS_TEMPERATURE_CARD",file_index);			        
		  }    		      		  
		}
		| error {}
		;
		
solidus_temperature_card:
		/* empty */ {} 
		| SOLIDUS_ TEMPERATURE_ EQUALS_ CONSTANT_ optional_floating_point_constant_list CR_ 
		{
		  switch(number_of_constants)
		  {
		    case 0:
                      sprintf(msg, "Material %s - expected one constant for the Solidus Temperature CONSTANT model.", pd_glob[mn]->MaterialName);
                      declare_error(msg);     		          
		      break;
		    case 1:
		      mat_ptr->melting_point_solidus = floating_point_constant_list_array[0];
		      mat_ptr->SolidusModel = CONSTANT;
		      card_read("SOLIDUS_TEMPERATURE_CARD",file_index);	 
		      break;
		    default:
                      sprintf(msg, "Material %s - expected only one constant for the Solidus Temperature CONSTANT model.  Subsequent constants being ignored.", pd_glob[mn]->MaterialName);
                      declare_warning(msg);   
		      mat_ptr->melting_point_solidus = floating_point_constant_list_array[0];
		      mat_ptr->SolidusModel = CONSTANT;
		      card_read("SOLIDUS_TEMPERATURE_CARD",file_index);			        
		  }    		      		  
		}
		| error {}
		;		
		
lame_lambda_model_name:	
		  CONSTANT_ 	{$$=CONSTANT;}
		| POWER_LAW_ 	{$$=POWER_LAW;}
		| USER_ 	{$$=USER;}
		| USER_GEN_	{$$=USER_GEN;}
		| EXPONENTIAL_	{$$=EXPONENTIAL;}
		| error {}
		;

stress_free_solvent_vol_frac_card:
		/* empty */ {} 
		| STRESS_ FREE_ SOLVENT_ VOL_ FRAC_ EQUALS_ CONSTANT_ optional_floating_point_constant_list CR_ 
		{
                  if (number_of_constants == 0)
                  {
                    sprintf(msg, "Material %s - expected one constant for the Stress Free Solvent Vol Frac CONSTANT model.", pd_glob[mn]->MaterialName);
                    declare_error(msg);                   
                  } 
                  else if (number_of_constants == 1)
                  {
		    elc_glob[mn]->Strss_fr_sol_vol_frac = floating_point_constant_list_array[0];
		    card_read("STRESS_FREE_SOLVENT_VOL_FRAC_CARD",file_index);	        		   
                  }
                  else
                  { 
                    sprintf(msg, "Material %s - expected only one constant for the Stress Free Solvent Vol Frac CONSTANT model.  Subsequent constants being ignored.", pd_glob[mn]->MaterialName);
                    declare_warning(msg);   
		    elc_glob[mn]->Strss_fr_sol_vol_frac = floating_point_constant_list_array[0];
		    card_read("STRESS_FREE_SOLVENT_VOL_FRAC_CARD",file_index);		                  
                  }        	
		}
		| error {}
		;
		
solid_thermal_expansion_card:
		/* empty */ {} 
		| SOLID_ THERMAL_ EXPANSION_ EQUALS_ solid_thermal_expansion_model optional_floating_point_constant_list CR_ 
		{
		  if ( ( ($5 == USER) || ($5 == USER_GEN) ) && number_of_constants == 0)
                  {
                    sprintf(msg, "Material %s - expected at least 1 constant for the Solid Thermal Expansion USER or USER_GEN model.", pd_glob[mn]->MaterialName);
                    declare_error(msg);                       
                  }		  
  		  else if($5 == CONSTANT )
		  {
		    switch(number_of_constants)
		    {
		      case 0:
                        sprintf(msg, "Material %s - expected one constant for the Solid Thermal Expansion CONSTANT model.", pd_glob[mn]->MaterialName);
                        declare_error(msg);     		          
		        break;
		      case 1:
		        elc_glob[mn]->thermal_expansion = floating_point_constant_list_array[0];
		        elc_glob[mn]->thermal_expansion_model = $5;
		        card_read("SOLID_THERMAL_EXPANSION_CARD",file_index);	 
		        break;
		      default:
                        sprintf(msg, "Material %s - expected only one constant for the Solid Thermal Exapnsion CONSTANT model.  Subsequent constants being ignored.", pd_glob[mn]->MaterialName);
                        declare_warning(msg);   
		        elc_glob[mn]->thermal_expansion = floating_point_constant_list_array[0];
		        elc_glob[mn]->thermal_expansion_model = $5;
		        card_read("SOLID_THERMAL_EXPANSION_CARD",file_index);			        
		    }    		      		  
		  }                
		  else
		  {
		    if( &(elc_glob[mn]->u_thermal_expansion) == NULL) 
	            {
	              elc_glob[mn]->u_thermal_expansion = (dbl *)array_alloc(1, number_of_constants, sizeof(dbl));
	            }
		    elc_glob[mn]->len_u_thermal_expansion = read_floats( &(elc_glob[mn]->u_thermal_expansion), 0);
		    elc_glob[mn]->thermal_expansion_model = $5;
		    /* model name is not saved...possible error in original mm_input_mp.c jjj ??? */
		    card_read("SOLID_THERMAL_EXPANSION_CARD",file_index);		  
		  }	
		}
		| error {}
		;		
		
solid_thermal_expansion_model:
		  CONSTANT_	{$$=CONSTANT;}
		| USER_		{$$=USER;}
		| USER_GEN_	{$$=USER_GEN;}
		| error		{}
		;
		
electrical_conductivity_card:
		/* empty */ {} 
		| ELECTRICAL_ CONDUCTIVITY_ EQUALS_ electrical_conductivity_model optional_floating_point_constant_list CR_ 
		{
		  if ( ( ($5 == USER) || ($5 == USER_GEN) ) && number_of_constants == 0)
                  {
                    sprintf(msg, "Material %s - expected at least 1 constant for the Solid Thermal Expansion USER or USER_GEN model.", pd_glob[mn]->MaterialName);
                    declare_error(msg);                       
                  }		  
  		  else if($5 == CONSTANT )
		  {
		    switch(number_of_constants)
		    {
		      case 0:
                        sprintf(msg, "Material %s - expected one constant for the Electrical Conductivity CONSTANT model.", pd_glob[mn]->MaterialName);
                        declare_error(msg);     		          
		        break;
		      case 1:
		        mat_ptr->electrical_conductivity = floating_point_constant_list_array[0];
		        mat_ptr->Elec_ConductivityModel = $4;
		        card_read("ELECTRICAL_CONDUCTIVITY_CARD",file_index);	 
		        break;
		      default:
                        sprintf(msg, "Material %s - expected only one constant for the Electrical Conductivity CONSTANT model.  Subsequent constants being ignored.", pd_glob[mn]->MaterialName);
                        declare_warning(msg);   
		        mat_ptr->electrical_conductivity = floating_point_constant_list_array[0];
		        mat_ptr->Elec_ConductivityModel = $4;
		        card_read("ELECTRICAL_CONDUCTIVITY_CARD",file_index);			        
		    }    		      		  
		  }                
		  else if ( ($4 == ELECTRONEUTRALITY_FICKIAN) || ($4 == ELECTRONEUTRALITY_SM) || ($4 == ELECTRODE_KINETICS) )
		  {
		    mat_ptr->Elec_ConductivityModel = $4;
		    card_read("ELECTRICAL_CONDUCTIVITY_CARD",file_index);
		  }
		  else
		  {
		    if( &(mat_ptr->u_electrical_conductivity) == NULL) 
	            {
	              mat_ptr->u_electrical_conductivity = (dbl *)array_alloc(1, number_of_constants, sizeof(dbl));
	            }
		    mat_ptr->len_u_electrical_conductivity = read_floats( &(mat_ptr->u_electrical_conductivity), 0);
		    mat_ptr->Elec_ConductivityModel = $4;
		    card_read("ELECTRICAL_CONDUCTIVITY_CARD",file_index);	  
		  }	
		}
		| error {}
		;
					
electrical_conductivity_model:
		  ELECTRONEUTRALITY_FICKIAN_	{$$=ELECTRONEUTRALITY_FICKIAN;}
		| ELECTRONEUTRALITY_SM_		{$$=ELECTRONEUTRALITY_SM;}
		| ELECTRODE_KINETICS_		{$$=ELECTRODE_KINETICS;}
		| CONSTANT_			{$$=CONSTANT;}
		| USER_				{$$=USER;}
		| error				{}
		;
		
media_type_card:
		/* empty */ {} 
		| MEDIA_ TYPE_ EQUALS_ media_type CR_ 
		{
                  if (!strcmp($4, "CONTINUOUS") || !strcmp($4, "NONE")   )            
	          {
	            mat_ptr->PorousMediaType = CONTINUOUS;
	            /*consistency checks */
                    if( pd_glob[mn]->e[R_POR_LIQ_PRES])
                    {
                      declare_error("Media Type must be porous.");
                    }
                    /* create a list of the porous media equations 
                    active in this material                     */      
                    i = 0;
                    /*old  for( j=0; j< pd_glob[mn]->Num_Porous_Eqn; j++ ) */
                    for( j=0; j< MAX_POROUS_NUM; j++ ) 
                    {
                      if( pd_glob[mn]->e[R_POR_LIQ_PRES+j] ) 
                      {
                         mat_ptr->Porous_Eqn[i] = R_POR_LIQ_PRES + j;
                         i++;
                      }
                    }
                    if( i != pd_glob[mn]->Num_Porous_Eqn )
                    {
	              declare_warning("Possible duplicate porous media equations.");
                    }                    
	          }
                  else if (   !strcmp($4, "POROUS_SATURATED") )
            	  {
            	    mat_ptr->PorousMediaType = POROUS_SATURATED;
                    /*consistency checks */
                    if( !pd_glob[mn]->e[R_POR_LIQ_PRES])
	            {
	              declare_error("You cannot run a porous media simulation without selecting the porous media equations.");
	            }
	            mat_ptr->CapStress = NO_MODEL; /*so proper effective stress principle is chosen*/
	          }
                  else if (   !strcmp($4, "POROUS_TWO_PHASE") )
            	  {
            	    mat_ptr->PorousMediaType = POROUS_TWO_PHASE;
                    if( !pd_glob[mn]->e[R_POR_LIQ_PRES] ||
                          !pd_glob[mn]->e[R_POR_GAS_PRES] )
            	    {
                      declare_error("You cannot run a two-phase porous media simulation without selecting both the porous_liq and porous_gas equations.");
	             }
            	  }
                  else if (  ( !strcmp($4, "POROUS_UNSATURATED") || 
		           !strcmp($4, "POROUS_PART_SAT") ) )
            	  {
            	    mat_ptr->PorousMediaType = POROUS_UNSATURATED;
                    if( !pd_glob[mn]->e[R_POR_LIQ_PRES])
	            {
	              declare_error("You cannot run a porous media simulation without selecting the porous media equations.");
	            }
            	  }
                  else if (   (!strcmp($4, "POROUS_BRINKMAN")) )
            	  {
            	    mat_ptr->PorousMediaType = POROUS_BRINKMAN;
                    mat_ptr->i_ys = 0;
            	  }
            	  card_read("MEDIA_TYPE_CARD",file_index);
		}
		| error {}
		;
				
media_type:
                  CONTINUOUS_		{}		
		| NONE_			{}		
		| POROUS_SATURATED_	{}		
		| POROUS_TWO_PHASE_	{}
		| POROUS_PART_SAT_	{}
		| POROUS_UNSATURATED_	{}
		| POROUS_BRINKMAN_	{}
		| error			{}
		;				

porosity_card:
		/* empty */ {} 
		| POROSITY_ EQUALS_ porosity_model optional_floating_point_constant_list CR_ 
		{
                  if ( mat_ptr->PorousMediaType == CONTINUOUS )
                  {		  
		    switch(number_of_constants)
		    {
		      case 0:
                        sprintf(msg, "Material %s - expected one constant for the Porosity Card.", pd_glob[mn]->MaterialName);
                        declare_error(msg);     		          
		        break;
		      case 1:
		        card_read("POROSITY_CARD",file_index);	 
		        break;
		      default:
                        sprintf(msg, "Material %s - expected only one constant for the Porosity Card. Subsequent constants being ignored.", pd_glob[mn]->MaterialName);
                        declare_warning(msg);   
		        card_read("POROSITY_CARD",file_index);			        
		    }   /* end switch(number_of_constants) */	
		    switch($3)
		    {
		      case CONSTANT:
		        mat_ptr->PorosityModel=CONSTANT;
		        mat_ptr->porosity = floating_point_constant_list_array[0];   		          
		        break;
   	              case DEFORM:
		        mat_ptr->PorosityModel=DEFORM;
		        mat_ptr->u_porosity[0] = floating_point_constant_list_array[0];
		        mat_ptr->len_u_porosity = 1;
		        if( !pd_glob[mn]->e[R_POR_POROSITY] )
	                {
                          declare_error("You cannot calculate porosity variation without the porosity equation.");
	                }
		        break;	
		        fprintf(stdout, "Active Porous Equations (1-true, 0 false): pl = %d, pg = %d, porosity = %d\n", 
	                  pd_glob[mn]->e[R_POR_LIQ_PRES], 
	                  pd_glob[mn]->e[R_POR_GAS_PRES],
	                  pd_glob[mn]->e[R_POR_POROSITY]);
		        
		    }  /* end switch($3) */ 
                    if(pd_glob[mn]->e[R_POR_POROSITY] && strcmp(mat_ptr->PorosityModel, "DEFORM")) 
	            {
	              declare_error("You must have a DEFORM porosity model with the porosity equation.");
	            }		    	    	 
		  }
		  else
		  {
		    declare_warning("Porous Media Type is not CONTINUOUS.  Porosity Card is being ignored.");
		  }   
		}
		| error {}
		;
porosity_model:
                  CONSTANT_	{$$=CONSTANT;}		
		| DEFORM_	{$$=DEFORM;}		
		| error		{}
		;

permeability_card:
		/* empty */ {} 
		| PERMEABILITY_ EQUALS_ permeability_model optional_floating_point_constant_list CR_ 
		{
		  mat_ptr->PermeabilityModel = $4; 
		  switch($3)
		  {
		    case CONSTANT:
		      switch(number_of_constants)
		      {
		        case 0:
                          sprintf(msg, "Material %s - expected one constant for the Permeability CONSTANT model.", pd_glob[mn]->MaterialName);
                          declare_error(msg);     		          
		          break;
		        case 1:
		          card_read("PERMEABILITY_CARD",file_index);	 
		          break;
		        default:
                          sprintf(msg, "Material %s - expected only one constant for the Permeability CONSTANT model.  Subsequent constants being ignored.", pd_glob[mn]->MaterialName);
                          declare_warning(msg);   
		          card_read("PEREMABILITY_CARD",file_index);			        
		      } /* end switch(number_of_constants) */ 
		      break;
		    case SOLIDIFICATION:
		      if (number_of_constants < 1  )
		      {
		        sprintf(msg,"Matl %s expected at least 1 constant for Permeability SOLIDIFICATION model.",
		            pd_glob[mn]->MaterialName);
                        declare_error(msg);
		      }
		      else
		      {
		        mat_ptr->len_u_permeability = read_floats( &(mat_ptr->u_permeability), 0);
		        card_read("PERMEABILITY_CARD",file_index);\
		      }
		      break;
		    case PSD_VOL:
		    case PSD_WEXP:
		    case PSD_SEXP:	
		      if (number_of_constants < 4  )
		      {
		        sprintf(msg,"Matl %s expected at least 4 constants for this Permeability model.",
		            pd_glob[mn]->MaterialName);
                         declare_error(msg);		          
		      }
		      else
		      {
		        mat_ptr->len_u_permeability = read_floats( &(mat_ptr->u_permeability), 0);
		        card_read("PERMEABILITY_CARD",file_index);		        
		      }
		      break;
   	            case KOZENY_CARMAN:		      	        
                      if (number_of_constants < 2 &&  ($3==KOZENY_CARMAN) )
		      {
		          sprintf(msg,"Matl %s expected at least 2 constants for Permeability KOZENY_CARMAN model.",
			      pd_glob[mn]->MaterialName);
                          declare_error(msg);
		      }
		      else
		      {
		        mat_ptr->len_u_permeability = read_floats( &(mat_ptr->u_permeability), 0);
		        card_read("PERMEABILITY_CARD",file_index);
		      }
		      break;
		    case K_TENSOR:
                      if (number_of_constants < 2  )
		      {
		        sprintf(msg,"Matl %s expected at least 2 constants for Permeability TENSOR model.",
		           pd_glob[mn]->MaterialName);
                        declare_error(msg);
		      }
		      else
		      {
		        mat_ptr->len_u_permeability = read_floats( &(mat_ptr->u_permeability), 0);
	                mat_ptr->perm_tensor[0][0]  = mat_ptr->u_permeability[0];
	                mat_ptr->perm_tensor[1][1]  = mat_ptr->u_permeability[1];
	                mat_ptr->perm_tensor[1][0]  = mat_ptr->u_permeability[2];
	                mat_ptr->perm_tensor[0][1]  = mat_ptr->u_permeability[3];
	                mat_ptr->permeability       = mat_ptr->u_permeability[0];
	                card_read("PERMEABILITY_CARD",file_index);		          
		     }
		  }  /* end switch($3) */
		}
		| error {}
		;
				
permeability_model:
                  CONSTANT_		{$$=CONSTANT;}		
		| PSD_VOL_		{$$=PSD_VOL;}		
		| PSD_WEXP_		{$$=PSD_WEXP;}		
		| PSD_SEXP_		{$$=PSD_SEXP;}
		| SOLIDIFICATION_	{$$=SOLIDIFICATION;}
		| TENSOR_		{$$=K_TENSOR;}
		| KOZENY_CARMAN_	{$$=KOZENY_CARMAN;}
		| error			{}
		;
		
inertia_coefficient_card:
		/* empty */ {} 
		| INERTIA_ COEFFICIENT_  EQUALS_ CONSTANT_ floating_point_constant_list CR_ 
		{
	          if( mat_ptr->PorousMediaType == POROUS_BRINKMAN )
		  {
		    switch(number_of_constants)
		    {
		      case 0:
                        sprintf(msg, "Material %s - expected one constant for the Inertia Coefficient CONSTANT model.", pd_glob[mn]->MaterialName);
                        declare_error(msg);     		          
		        break;
		      case 1:
                        mat_ptr->InertiaCoefficientModel =  CONSTANT;
               	        mat_ptr->Inertia_coefficient = floating_point_constant_list_array[0];
		        card_read("INERTIA_COEFFICIENT_CARD",file_index);	 
		        break;
		      default:
                        sprintf(msg, "Material %s - expected only one constant for the Inertia Coefficient CONSTANT model.  Subsequent constants being ignored.", pd_glob[mn]->MaterialName);
                        declare_warning(msg);   
		        mat_ptr->Inertia_coefficient = floating_point_constant_list_array[0];
		        mat_ptr->InertiaCoefficientModel =  CONSTANT;
		        card_read("INERTIA_COEFFICIENT_CARD",file_index);			        
		    }    		      		  
		  }                
		  else
		  {
                    sprintf(msg, "Meida type is not POROUS_BRINKMAN.  Ignoring Inertia Coefficient card and defaulting inertia coefficient to 0.0.");
                    declare_warning(msg);  	  
                    mat_ptr->Inertia_coefficient = 0.0;
		  }	
		}
		| error {}
		;				

porous_diffusion_constitutive_equation_card:
                /* empty*/ {}
                | POROUS_ DIFFUSION_ CONSTITUTIVE_ EQUATION_ EQUALS_ porous_diffusion_constitutive_equation_model CR_
                {
                  if (pd_glob[mn]->Num_Porous_Eqn > 0)
                  {
 		    switch($6)
		    {
		      case FICKIAN:
                        sprintf(msg, "Material %s - expected one constant for the Inertia Coefficient CONSTANT model.", pd_glob[mn]->MaterialName);
                        declare_error(msg);     		          
		        break;
		      case DARCY_FICKIAN:
                        pd_glob[mn]->PorousFluxModel = $6;
	                cr_glob[mn]->PorousFluxModel = $6;
		        card_read("POROUS_DIFFUSION_CONSTITUTIVE_EQUATION_CARD",file_index);	 
		        break;
                    }
                  }
                  else
                  {
                    sprintf(msg,"The number of porous equations is zero, Porous Diffusion Constitutive Equation card being ignored.");
		    declare_warning(msg);
                  }
                }
                | error {}
                ;
                
porous_diffusion_constitutive_equation_model:
		  FICKIAN_		{$$ = FICKIAN;}
		| DARCY_FICKIAN_	{$$ = DARCY_FICKIAN;}
		| error {}
		;	

diffusion_constitutive_equation_card:
                /* empty*/ {}
                | DIFFUSION_ CONSTITUTIVE_ EQUATION_ EQUALS_ diffusion_constitutive_equation_model CR_
                {
                  if (pd_glob[mn]->Num_Porous_Eqn > 0)
                  {
 		    switch($5)
		    {
		      case HYDRODYNAMIC:
                        if( !pd_glob[mn]->e[R_SHEAR_RATE] )
                        {
                          declare_error("HYDRODYNAMIC mass flux requires shear_rate dof in EQ list.");  
                        }
                        else
                        {
                          card_read("DIFFUSION_CONSTITUTIVE_EQUATION_CARD",file_index); 
                        }   		          
		        break;
		      case DM_SUSPENSION_BALANCE:
		        if( !pd_glob[mn]->e[R_GRADIENT11] )
                        {
                          declare_error("SUSPENSION_BALANCE mass flux requires velocity gradient tensor in EQ list.");  
                        }
                        else
                        {
                          card_read("DIFFUSION_CONSTITUTIVE_EQUATION_CARD",file_index); 
                        }   		          
		        break;
		      default:
		        card_read("DIFFUSION_CONSTITUTIVE_EQUATION_CARD",file_index);		      
                    }
                    pd_glob[mn]->MassFluxModel = $5;
                    cr_glob[mn]->MassFluxModel = $5;              
                  }
                  else
                  {
                    sprintf(msg,"The number of porous equations is zero, Diffusion Constitutive Equation card being ignored.");
		    declare_warning(msg);
                  }
                }
                | error {}
                ;
                
diffusion_constitutive_equation_model:
                  FICKIAN_		{$$ = FICKIAN;}
		| DARCY_FICKIAN_	{$$ = DARCY_FICKIAN;}
		| GENERALIZED_FICKIAN_	{$$ = GENERALIZED_FICKIAN;}
		| STEFAN_MAXWELL_	{$$ = STEFAN_MAXWELL;}
		| STEFAN_MAXWELL_CHARGED_	{$$ = STEFAN_MAXWELL_CHARGED;}
		| STEFAN_MAXWELL_VOLUME_	{$$ = STEFAN_MAXWELL_VOLUME;}
		| FICKIAN_CHARGED_	{$$ = FICKIAN_CHARGED;}
		| DARCY_		{$$ = DARCY;}
		| HYDRODYNAMIC_		{$$ = HYDRODYNAMIC;}
		| SUSPENSION_BALANCE_	{$$ = DM_SUSPENSION_BALANCE;}
		| NONE_			{$$ = NON_DIFFUSING;}																				
		| error {}
		;

flowingliquid_viscosity_card:
		/* empty */ {} 
		| FLOWINGLIQUID_ VISCOSITY_ EQUALS_ flowingliquid_viscosity_model optional_floating_point_constant_list CR_ 
		{
                  if ( ($5 == USER)  && number_of_constants == 0)
                  {
                    sprintf(msg, "Material %s - expected at least 1 constant for the FlowingLiquid Viscosity USER model.", pd_glob[mn]->MaterialName);
                    declare_error(msg);                       
                  }		  
  		  else if($5 == CONSTANT )
		  {
		    switch(number_of_constants)
		    {
		      case 0:
                        sprintf(msg, "Material %s - expected one constant for the FlowingLiquid Viscosity CONSTANT model.", pd_glob[mn]->MaterialName);
                        declare_error(msg);     		          
		        break;
		      case 1:
		        mat_ptr->FlowingLiquid_viscosity = floating_point_constant_list_array[0];
		        mat_ptr->FlowingLiquidViscosityModel = $4;
		        card_read("FLOWINGLIQUID_VISCOSITY_CARD",file_index);	 
		        break;
		      default:
                        sprintf(msg, "Material %s - expected only one constant for the FlowingLiquid Viscosity CONSTANT model.  Subsequent constants being ignored.", pd_glob[mn]->MaterialName);
                        declare_warning(msg);   
		        mat_ptr->FlowingLiquid_viscosity = floating_point_constant_list_array[0];
		        mat_ptr->FlowingLiquidViscosityModel = $4;
		        card_read("FLOWINGLIQUID_VISCOSITY_CARD",file_index);			        
		    }    		      		  
		  }                
		  else if ( $4 == MOLTEN_GLASS )
		  {
		    if (number_of_constants < 3) 
		    {
	              sprintf(msg, "Matl %s needs at least 3 constants for FlowingLiquid Viscosity MOLTEN_GLASS model.", pd_glob[mn]->MaterialName);
	              declare_error(msg);
                    }
                  }
		  else
		  {
		    if( &(mat_ptr->u_FlowingLiquid_viscosity) == NULL) 
	            {
	              mat_ptr->u_FlowingLiquid_viscosity = (dbl *)array_alloc(1, number_of_constants, sizeof(dbl));
	            }
		    mat_ptr->len_u_FlowingLiquid_viscosity = read_floats( &(mat_ptr->u_FlowingLiquid_viscosity), 0);
		    mat_ptr->FlowingLiquidViscosityModel = $4;
		    card_read("FLOWINGLIQUID_VISCOSITY_CARD",file_index);	  
		  }
		}
		| error {}
		;
				
flowingliquid_viscosity_model:
                  CONSTANT_		{$$=CONSTANT;}		
		| USER_			{$$=USER;}		
		| MOLTEN_GLASS_		{$$=MOLTEN_GLASS;}		 
		| error			{}
		;			
						
reference_temperature_card:
		/* empty */ {} 
		| REFERENCE_ TEMPERATURE_ EQUALS_ CONSTANT_ optional_floating_point_constant_list CR_ 
		{
		  switch(number_of_constants)
		  {
		    case 0:
                      sprintf(msg, "Material %s - expected one constant for the Solid Reference Temperature CONSTANT model.", pd_glob[mn]->MaterialName);
                      declare_error(msg);     		          
		      break;
		    case 1:
		      elc_glob[mn]->solid_reference_temp = floating_point_constant_list_array[0];
		      elc_glob[mn]->solid_reference_temp_model = CONSTANT;
		      card_read("REFERENCE_TEMPERATURE_CARD",file_index);	 
		      break;
		    default:
                      sprintf(msg, "Material %s - expected only one constant for the Solid Reference Temperature CONSTANT model.  Subsequent constants being ignored.", pd_glob[mn]->MaterialName);
                      declare_warning(msg);   
		      elc_glob[mn]->solid_reference_temp = floating_point_constant_list_array[0];
		      elc_glob[mn]->solid_reference_temp_model = CONSTANT;
		      card_read("REFERENCE_TEMPERATURE_CARD",file_index);			        
		  }    		      		  
		}
		| error {}
		;
		
diffusivity_card:
		/* empty */ {}
		| DIFFUSIVITY_ EQUALS_ diffusivity_model integer optional_floating_point_constant_list CR_
		{
		  if (pd_glob[mn]->MassFluxModel != STEFAN_MAXWELL         &&
                      pd_glob[mn]->MassFluxModel != STEFAN_MAXWELL_CHARGED &&
                      pd_glob[mn]->MassFluxModel != STEFAN_MAXWELL_VOLUME    )
                  {  
                    fallback_chemkin_generic_prop(&model_read, $4, &(mat_ptr->DiffusivityModel[j]),FALSE, mat_ptr);
 		    if ($3==CONSTANT)
		    {
		      switch(number_of_constants)
		      {
		        case 0:
                          sprintf(msg, "Material %s - expected one constant for the Diffusivity CONSTANT model.", pd_glob[mn]->MaterialName);
                          declare_error(msg);     		          
		          break;
		        case 1:
		          mat_ptr->DiffusivityModel[$4] = CONSTANT;
			  mat_ptr->diffusivity[$4] = floating_point_constant_list_array[0];
		          card_read("DIFFUSIVITY_CARD",file_index);	 
		          break;
		        default:
                          sprintf(msg, "Material %s - expected only one constant for the Diffusivity CONSTANT model.  Subsequent constants being ignored.", pd_glob[mn]->MaterialName);
                          declare_warning(msg);   
		          mat_ptr->DiffusivityModel[$4] = CONSTANT;
			  mat_ptr->diffusivity[$4] = floating_point_constant_list_array[0];
		          card_read("DIFFUSIVITY_CARD",file_index);		        
		      }	
		    }
		    else if ( ($3 == HYDRO) || ($3 == SUSP_BAL) )
		    {
		      mat_ptr->DiffusivityModel[$4] = $3;
		      card_read("DIFFUSIVITY_CARD",file_index);
		    } 
		    else
		    {
		      switch ($3)
		      {
		        case USER:
		          if (number_of_constants < 1) 
		          {
	                    sprintf(msg, "Matl %s (conc %d) needs at least 1 constants for %s %s model.\n",
			            pd_glob[mn]->MaterialName, $4,
			            "Diffusivity", "CONSTANT");
	                    declare_error(msg);	
	                  }	        
		          break;
		        case FREE_VOL:
		          if (number_of_constants < 12) 
		          {
	                    sprintf(msg, "Matl %s (conc %d) needs at least 12 constants for %s %s model.\n",
			            pd_glob[mn]->MaterialName, $4,
			            "Diffusivity", "FREE_VOL");
	                    declare_error(msg);
	                  }
	                  if( pd_glob[mn]->Num_Species_Eqn != 1 ) 
	                  {
	                    declare_error("Binary Free volume models are for 2 components, or one tracked species.");
	                  }		        
		          break;
		        case GENERALIZED_FREE_VOL: /* LEVEL_SET is also 12 jjj ???*/
		          if (number_of_constants < 12) 
		          {
	                    sprintf(msg, "Matl %s (conc %d) needs at least 12 constants for %s %s model.",
			            pd_glob[mn]->MaterialName, $4,
			            "Diffusivity", "GENERALIZED_FREE_VOL");
	                    declare_error(msg);
	                  }
	                  if( pd_glob[mn]->Num_Species_Eqn < 1 ) 
	                  {
	                    declare_error("Generalized models are for 2 or more BULK components.");
	                  }		        
		          break;
		        case GENERALIZED:
		          if (number_of_constants <  pd_glob[mn]->Num_Species_Eqn ) 
		          {
	                    sprintf(msg,"Matl %s (conc %d) needs one constant for each i-j pair %s %s model.\n",
			      pd_glob[mn]->MaterialName, $4,
			      "Diffusivity", "GENERALIZED");
	                    declare_error(msg);
	                  }
	                  if( pd_glob[mn]->Num_Species_Eqn < 2 ) 
	                  {
	                    declare_error("Generalized diffusivity model is for 2 or more BULK components.");
	                  }
		          break;
		        case POROUS:
                          if (number_of_constants < 5) 
                          {
	                    sprintf(msg, "Matl %s (conc %d) needs at least 5 constants for %s %s model.\n",
			            pd_glob[mn]->MaterialName, $4,
			            "Diffusivity", "POROUS");
			    declare_error(msg);
			  }
		          break;
		      } /* end switch on model $3 */
  
                      if (mat_ptr->u_diffusivity[$4] == NULL) 
                      {
                        sprintf(msg,"Space for pointer vector over species needs to be malloced first.");
                        declare_error(msg);
                      }
                      else
                      {
                        mat_ptr->DiffusivityModel[$4] = $3;    		
                        mat_ptr->len_u_diffusivity[$4] = read_floats( &(mat_ptr->u_diffusivity[$4]), 0);
                        card_read("DIFFUSIVITY_CARD",file_index);		      
                      }
                    }                    		    
		  }
		  else
		  {
		    sprintf(msg, "pd_glob[mn]->MassFluxModel = %i.  Diffusivity card for species %i being ignored.",pd_glob[mn]->MassFluxModel,$4);
                    declare_warning(msg);  
		  }
		}
		| DIFFUSIVITY_ EQUALS_ TABLE_ integer  mp_table_independent_variable_name  		/* $1 $2 $3 $4 $5 */
		                                       optional_mp_table_independent_variable_name 	/* $6 */
		                                       optional_mp_table_independent_variable_name 	/* $7 */
		                                       mp_table_interpolation_method 			/* $8 */
		                                       optional_file_equals_string 			/* $9 */
		                                       optional_name_equals_string			/* $10 */
		                                       CR_
		{
                  if( num_MP_Tables == MAX_MP_TABLES )
	          {
	            sprintf(msg, "Material %s - Maximum TABLE_MPs exceeded.", pd_glob[mn]->MaterialName);
                    declare_error(msg);
	          }
	          else
                  {
		    table_type = DIFFUSIVITY_TABLE;
		    MP_Tables[num_MP_Tables] = setup_mp_table(mat_ptr->table,$4,"Diffusivity",$5,$6,$7,$8,mp_table_species_number);
                    if ( (strcmp( $9, "NA") == 0) )
                    { 
                      /* Don't read table data from another file */
                      accept_table_data_cards=TRUE;
                    }
                    else
                    {
                      /* Try to read table data from another file */ 
                      accept_table_data_cards=FALSE;
		      table_type = DIFFUSIVITY_TABLE;
		      parse_table_file($9,$10,$4);
    	              read_array_into_table(table_type, mat_ptr->table, number_of_table_cards_read);
                    }		    
		  }
		  card_read("DIFFUSIVITY_CARD",file_index);
		}
		| error {}
		;	
				
diffusivity_model:
		  POROUS_ 		{$$=POROUS;}
                | CONSTANT_ 		{$$=CONSTANT;}
                | HYDRO_ 		{$$=HYDRO;}
                | SUSPENSION_ 		{$$=SUSP_BAL;}
                | FREE_VOL_ 		{$$=FREE_VOL;}
                | GENERALIZED_FREE_VOL_ {$$=GENERALIZED_FREE_VOL;}
                | GENERALIZED_ 		{$$=GENERALIZED;}
                | LEVEL_SET_ 		{$$=LEVEL_SET;}
                | USER_		  	{$$=USER;}
		| error			{}
		;


porous_latent_heat_vaporization_card: 
                /* empty*/ {}
                | POROUS_ LATENT_ HEAT_ VAPORIZATION_ EQUALS_ CONSTANT_ integer optional_floating_point_constant_list CR_
                {
                  if( mat_ptr->PorousMediaType != POROUS_SATURATED )
                  {
                    if($7 != 0)
                    {
                      sprintf(msg, "Material %s - Porous Latent Heat Vaporization species number must be 0.", pd_glob[mn]->MaterialName);
                      declare_error(msg);                   
                    }
		    switch(number_of_constants)
		        {
		          case 0:
                            sprintf(msg, "Material %s - expected one constant for the Porous Latent Heat Vaporization CONSTANT model.", pd_glob[mn]->MaterialName);
                            declare_error(msg);     		          
		            break;
		          case 1:
		            mat_ptr->PorousLatentHeatVapModel[$7] = CONSTANT;
		            mat_ptr->porous_latent_heat_vap[$7] = floating_point_constant_list_array[0];
		            card_read("POROUS_LATENT_HEAT_VAPORIZATION_CARD",file_index);	 
		            break;
		          default:
                            sprintf(msg, "Material %s - expected only one constant for the Porous Latent Heat Vaporization CONSTANT model.  Subsequent constants being ignored.", pd_glob[mn]->MaterialName);
                            declare_warning(msg);   
		            mat_ptr->PorousLatentHeatVapModel[$7] = CONSTANT;
		            mat_ptr->porous_latent_heat_vap[$7] = floating_point_constant_list_array[0];
		            card_read("POROUS_LATENT_HEAT_VAPORIZATION_CARD",file_index);		        
		        } /* end switch number_of_constants */	
                  }
                  else
                  {
		    declare_warning("Porous Media Type is POROUS_SATURATED; Porous Latent Heat Vaporization card is being ignored.");
                  }
                }
                | error {}
                ;
                					
latent_heat_vaporization_card: 
                /* empty*/ {}
                | LATENT_ HEAT_ VAPORIZATION_ EQUALS_ CONSTANT_ integer optional_floating_point_constant_list CR_
                {
                  switch(number_of_constants)
		  {
		    case 0:
                      sprintf(msg, "Material %s - expected one constant for the Latent Heat Vaporization CONSTANT model.", pd_glob[mn]->MaterialName);
                      declare_error(msg);     		          
		      break;
		    case 1:
		      mat_ptr->LatentHeatVapModel[$6] = CONSTANT;
		      mat_ptr->latent_heat_vap[$6] = floating_point_constant_list_array[0];
		      card_read("LATENT_HEAT_VAPORIZATION_CARD",file_index);	 
		      break;
		    default:
                      sprintf(msg, "Material %s - expected only one constant for the Latent Heat Vaporization CONSTANT model.  Subsequent constants being ignored.", pd_glob[mn]->MaterialName);
                      declare_warning(msg);   
		      mat_ptr->LatentHeatVapModel[$6] = CONSTANT;
		      mat_ptr->latent_heat_vap[$6] = floating_point_constant_list_array[0];
		      card_read("LATENT_HEAT_VAPORIZATION_CARD",file_index);		        
		  } /* end switch number_of_constants */
		  fallback_chemkin_generic_prop(&model_read, $6, &(mat_ptr->LatentHeatVapModel[$6]),TRUE, mat_ptr);
                }
                | error {}
                ;
    
latent_heat_fusion_card: 
                /* empty*/ {}
                | LATENT_ HEAT_ FUSION_ EQUALS_ CONSTANT_ integer optional_floating_point_constant_list CR_
                {
                  switch(number_of_constants)
		  { 
		    case 0:
                      sprintf(msg, "Material %s - expected one constant for the Latent Heat Fusion CONSTANT model.", pd_glob[mn]->MaterialName);
                      declare_error(msg);     		          
		      break;
		    case 1:
		      mat_ptr->LatentHeatFusionModel[$6] = CONSTANT;
		      mat_ptr->latent_heat_fusion[$6] = floating_point_constant_list_array[0];
		      card_read("LATENT_HEAT_FUSION_CARD",file_index);	 
		      break;
		    default:
                      sprintf(msg, "Material %s - expected only one constant for the Latent Heat Fusion CONSTANT model.  Subsequent constants being ignored.", pd_glob[mn]->MaterialName);
                      declare_warning(msg);   
		      mat_ptr->LatentHeatFusionModel[$6] = CONSTANT;
		      mat_ptr->latent_heat_fusion[$6] = floating_point_constant_list_array[0];
		      card_read("LATENT_HEAT_FUSION_CARD",file_index);		        
		  } /* end switch number_of_constants */
		  fallback_chemkin_generic_prop(&model_read, $6, &(mat_ptr->LatentHeatFusionModel[$6]),TRUE, mat_ptr);
                }
                | error {}
                ;    

species_volume_expansion_card: 
                /* empty*/ {}
                | SPECIES_ VOLUME_ EXPANSION_ EQUALS_ CONSTANT_ integer optional_floating_point_constant_list CR_
                {
                  switch(number_of_constants)
		  { 
		    case 0:
                      sprintf(msg, "Material %s - expected one constant for the Species Volume Expansion CONSTANT model.", pd_glob[mn]->MaterialName);
                      declare_error(msg);     		          
		      break;
		    case 1:
		      mat_ptr->SpecVolExpModel[$6] = CONSTANT;
		      mat_ptr->species_vol_expansion[$6] = floating_point_constant_list_array[0];
		      card_read("SPECIES_VOLUME_EXPANSION_CARD",file_index);	 
		      break;
		    default:
                      sprintf(msg, "Material %s - expected only one constant for the Species Volume Expansion CONSTANT model.  Subsequent constants being ignored.", pd_glob[mn]->MaterialName);
                      declare_warning(msg);   
		      mat_ptr->SpecVolExpModel[$6] = CONSTANT;
		      mat_ptr->species_vol_expansion[$6] = floating_point_constant_list_array[0];
		      card_read("SPECIES_VOLUME_EXPANSION_CARD",file_index);		        
		  } /* end switch number_of_constants */
		  fallback_chemkin_generic_prop(&model_read, $6,&(mat_ptr->SpecVolExpModel[$6]), TRUE, mat_ptr);
                }
                | error {}
                ;  
                
reference_concentration_card: 
                /* empty*/ {}
                | REFERENCE_ CONCENTRATION_ EQUALS_ CONSTANT_ integer optional_floating_point_constant_list CR_
                {
                  switch(number_of_constants)
		  { 
		    case 0:
                      sprintf(msg, "Material %s - expected one constant for the Reference Concentration CONSTANT model.", pd_glob[mn]->MaterialName);
                      declare_error(msg);     		          
		      break;
		    case 1:
		      mat_ptr->SpecVolExpModel[$5] = CONSTANT;
		      mat_ptr->species_vol_expansion[$5] = floating_point_constant_list_array[0];
		      card_read("REFERENCE_CONCENTRATION_CARD",file_index);	 
		      break;
		    default:
                      sprintf(msg, "Material %s - expected only one constant for the Reference Concentration CONSTANT model.  Subsequent constants being ignored.", pd_glob[mn]->MaterialName);
                      declare_warning(msg);   
		      mat_ptr->SpecVolExpModel[$5] = CONSTANT;
		      mat_ptr->species_vol_expansion[$5] = floating_point_constant_list_array[0];
		      card_read("REFERENCE_CONCENTRATION_CARD",file_index);		        
		  } /* end switch number_of_constants */
		  fallback_chemkin_generic_prop(&model_read, $5,&(mat_ptr->SpecVolExpModel[$5]), TRUE, mat_ptr);
                }
                | error {}
                ;                  
                

vapor_pressure_card:
		/* empty */ {} 
		| VAPOR_ PRESSURE_ EQUALS_ vapor_pressure_model integer optional_floating_point_constant_list CR_ 
		{
		  fallback_chemkin_generic_prop(&model_read, $5,&(mat_ptr->VaporPressureModel[$5]), TRUE, mat_ptr);
		  if ($4==CONSTANT)
		  {
		    switch(number_of_constants)
		    {
		      case 0:
                        sprintf(msg, "Material %s - expected one constant for the Vapor Pressure CONSTANT model.", pd_glob[mn]->MaterialName);
                        declare_error(msg);     		          
		        break;
		      case 1:
		        mat_ptr->VaporPressureModel[$5] = CONSTANT;
		        mat_ptr->vapor_pressure[$5] = floating_point_constant_list_array[0];
		        card_read("VAPOR_PRESSURE_CARD",file_index);	 
		        break;
		      default:
                        sprintf(msg, "Material %s - expected only one constant for the Vapor Pressure CONSTANT model.  Subsequent constants being ignored.", pd_glob[mn]->MaterialName);
                        declare_warning(msg);   
		        mat_ptr->VaporPressureModel[$5] = CONSTANT;
		        mat_ptr->vapor_pressure[$5] = floating_point_constant_list_array[0];
		        card_read("VAPOR_PRESSURE_CARD",file_index);	
		    }	  /* end switch number_of_constants */      
		  }
		  else
		  {
		     switch($4)
		     {
                       case KELVIN:
		       case FLAT:
		         mat_ptr->VaporPressureModel[$5] = $4;
		         switch(mat_ptr->PorousMediaType)
		         {
		           case POROUS_TWO_PHASE:
		             if (number_of_constants < 5)
		             {
		               declare_error( "KELVIN & FLAT Vapor Pressure models require at least five constants for the POROUS_TWO_PHASE Media Type.");  		             
		             }
		             else
		             {
		               card_read("VAPOR_PRESSURE_CARD",file_index);
		             }
		             break;
		           case POROUS_UNSATURATED:
		             if (number_of_constants < 5)
		             {
		               declare_error( "KELVIN & FLAT Vapor Pressure models require at least seven constants for the POROUS_UNSATURATED Media Type.");  		             
		             }		           
		             else
		             {
		               card_read("VAPOR_PRESSURE_CARD",file_index);
		             }
		             break;
		          default:
		             sprintf(msg, "KELVIN & FLAT Vapor Pressure models are invalid in matl %s unless Media Type is POROUS_TWO_PHASE or POROUS_UNSATURATED.", pd_glob[mn]->MaterialName);
			     declare_error(msg);		           
		         }  /* end switch mat_ptr->PorousMediaType */
		         break;
		       case NON_VOLATILE:
		         break;
		       case ANTOINE:
		         if (number_of_constants < 6)
		         {
		           declare_error( "ANTOINE Vapor Pressure model requires at least six constants for the POROUS_UNSATURATED Media Type.");  		             
		         }		           
		         else
		         {
		           card_read("VAPOR_PRESSURE_CARD",file_index);
		         }		       
		         break;
		       case RIEDEL:
		         if (number_of_constants < 8)
		         {
		           declare_error( "RIEDEL Vapor Pressure model requires at least eight constants.");  		             
		         }		           
		         else
		         {
		           card_read("VAPOR_PRESSURE_CARD",file_index);
		         }
		         break;
		       case IDEAL_GAS:
		         if (number_of_constants < 5)
		         {
		           declare_error( "IDEAL_GAS Vapor Pressure model requires at least three constants.");  		             
		         }		           
		         else
		         {
		           card_read("VAPOR_PRESSURE_CARD",file_index);
		         }
		         break;		       
		     } /* end switch $4 which is the model */	
		     mat_ptr->len_u_vapor_pressure[$5] = read_floats( &(mat_ptr->u_vapor_pressure[$5]), 0);	    
		  }			
		}
		|
		error {}
		;
                	
vapor_pressure_model:
		  CONSTANT_	{$$  = CONSTANT;}
		| KELVIN_	{$$  = KELVIN;}
		| FLAT_		{$$  = FLAT;}
		| NON_VOLATILE_	{$$  = NON_VOLATILE;}
		| ANTOINE_	{$$  = ANTOINE;}
		| RIEDEL_  	{$$  = RIEDEL;}   
		| IDEAL_GAS_	{$$  = IDEAL_GAS;}          	
		| error {}
		;
                		               
plastic_viscosity_card:
		/* empty */ {} 
		| PLASTIC_ VISCOSITY_ EQUALS_ plastic_viscosity_model optional_floating_point_constant_list CR_ 
		{	  
  		  if($4 == CONSTANT )
		  {
		    switch(number_of_constants)
		    {
		      case 0:
                        sprintf(msg, "Material %s - expected one constant for the Plastic Viscosity CONSTANT model.", pd_glob[mn]->MaterialName);
                        declare_error(msg);     		          
		        break;
		      case 1:
		        evpl_glob[mn]->plastic_mu = floating_point_constant_list_array[0];
		        evpl_glob[mn]->plastic_mu_model = $4;
		        card_read("PLASTIC_VISCOSITY_CARD",file_index);	 
		        break;
		      default:
                        sprintf(msg, "Material %s - expected only one constant for the Plastic Viscosity CONSTANT model.  Subsequent constants being ignored.", pd_glob[mn]->MaterialName);
                        declare_warning(msg);   
		        evpl_glob[mn]->plastic_mu = floating_point_constant_list_array[0];
		        evpl_glob[mn]->plastic_mu_model = $4;
		        card_read("PLASTIC_VISCOSITY_CARD",file_index);			        
		    }    		      		  
		  }                
		  else /* model has to be LINEAR which needs at least two constants*/
		  {
		    if (number_of_constants >= 2)
		    {
		      evpl_glob[mn]->plastic_mu_model = LINEAR;
		      evpl_glob[mn]->len_u_plastic_mu = read_floats( &(evpl_glob[mn]->u_plastic_mu), 0);
		      card_read("PLASTIC_VISCOSITY_CARD",file_index);	
                    }
		    else
		    {
		      sprintf(msg, 
	                     "Matl %s expected at least 2 constants for Plastic Viscosity LINEAR model.",
			     pd_glob[mn]->MaterialName);
		      declare_error(msg);
		    }
		  }	
		}
		| error {}
		;		
		
plastic_viscosity_model:
		  CONSTANT_	{$$=CONSTANT;}
		| LINEAR_	{$$=LINEAR;}
		| error		{}
		;


		
EVP_yeild_stress_card:
		/* empty */ {} 
		| EVP_ YEILD_ STRESS_ EQUALS_ EVP_yeild_stress_model optional_floating_point_constant_list CR_ 
		{	  
  		  if($5 == CONSTANT )
		  {
		    switch(number_of_constants)
		    {
		      case 0:
                        sprintf(msg, "Material %s - expected one constant for the EVP Yeild Stress CONSTANT model.", pd_glob[mn]->MaterialName);
                        declare_error(msg);     		          
		        break;
		      case 1:
		        evpl_glob[mn]->yield = floating_point_constant_list_array[0];
		        evpl_glob[mn]->yield_model = $5;
		        card_read("EVP_YEILD_STRESS_CARD",file_index);	 
		        break;
		      default:
                        sprintf(msg, "Material %s - expected only one constant for the EVP Yeild Stress CONSTANT model.  Subsequent constants being ignored.", pd_glob[mn]->MaterialName);
                        declare_warning(msg);   
		        evpl_glob[mn]->yield = floating_point_constant_list_array[0];
		        evpl_glob[mn]->yield_model = $5;
		        card_read("EVP_YEILD_STRESS_CARD",file_index);			        
		    }    		      		  
		  }                
		  else /* model has to be LINEAR which needs at least two constants*/
		  {
		    if (number_of_constants >= 2)
		    {
		      evpl_glob[mn]->yield_model = LINEAR;
		      evpl_glob[mn]->len_u_yield = read_floats( &(evpl_glob[mn]->u_yield), 0);
		      card_read("EVP_YEILD_STRESS_CARD",file_index);	
                    }
		    else
		    {
		      sprintf(msg, 
	                     "Matl %s expected at least 2 constants for EVP Yeild Stress LINEAR model.",
			     pd_glob[mn]->MaterialName);
		      declare_error(msg);
		    }
		  }	
		}
		| error {}
		;		
		
EVP_yeild_stress_model:
		  CONSTANT_	{$$=CONSTANT;}
		| LINEAR_	{$$=LINEAR;}
		| error		{}
		;
		
pseudo_solid_constitutive_equation_card:
                /* empty*/ {}
                | PSEUDO_SOLID_ CONSTITUTIVE_ EQUATION_ EQUALS_ pseudo_solid_constitutive_equation_model CR_
                {
                  if( pd_glob[mn]->MeshMotion == TOTAL_ALE)
                  {
                    cr_glob[mn]->MeshFluxModel = $5;
                    card_read("PSEUDO_SOLID_CONSTITUTIVE_EQUATION_CARD", mn);
                  }
                  else
                  {
                    sprintf(msg, 
	                     "Mesh Motion is not TOTAL_ALE, Pseudo-solid Constitutive Equation card being ignored.");
		    declare_warning(msg);
                  }
                }
                | error {}
                ;
pseudo_solid_constitutive_equation_model:
		  LINEAR_		{$$ = LINEAR;}
		| NON_LINEAR_		{$$ =  NONLINEAR;}
		| NONLINEAR_PLANE_STRAIN_	{$$ =  NONLINEAR;}
		| INCOMP_PSTRAIN_	{$$ = INCOMP_PSTRAIN ;}
		| INCOMP_PSTRESS_	{$$ =  INCOMP_PSTRESS;}
		| HOOKEAN_PSTRAIN_	{$$ =  HOOKEAN_PSTRAIN;}
		| HOOKEAN_PSTRESS_	{$$ =  HOOKEAN_PSTRESS;}
		| INCOMP_3D_		{$$ = INCOMP_3D;}
		| EVP_HYPER_		{$$ = EVP_HYPER;}
		| error {}
		;		
								
pseudo_solid_lame_mu_card:
                  /* empty */ {}
                  | PSEUDO_SOLID_ LAME_ MU_ EQUALS_ pseudo_solid_lame_mu_model optional_floating_point_constant_list CR_
                  {
                    if( pd_glob[mn]->MeshMotion == TOTAL_ALE)
                    { 
                      if ($5 == CONTACT_LINE_ && number_of_constants < 4 )
		      {
                        sprintf(msg, "Material %s - expected at least 4 constants for the pseudo_solid Lame MU CONTACT_LINE model.", pd_glob[mn]->MaterialName);
                        declare_error(msg);               		  
		      }
                      else if ( ( ($5 == USER) || ($5 == USER_GEN) ) && number_of_constants == 0)
                      {
                        sprintf(msg, "Material %s - expected at least 1 constant for the Pseudo_solid Lame MU USER or USER_GEN model.", pd_glob[mn]->MaterialName);
                        declare_error(msg);                       
                      }		  
  		      else if($5 == CONSTANT )
  		      {
		        switch(number_of_constants)
		        {
		          case 0:
                            sprintf(msg, "Material %s - expected one constant for the Pseudo-solid Lame MU CONSTANT model.", pd_glob[mn]->MaterialName);
                            declare_error(msg);      		          
		            break;
		          case 1:		          
		            elc_glob[mn]->len_u_mu = floating_point_constant_list_array[0];
		            elc_glob[mn]->lame_mu_model = $5;
		            card_read("PSEUDO_SOLID_LAME_MU_CARD",file_index);	 
		            break;
		          default:
                            sprintf(msg, "Material %s - expected only one constant for the Pseudo-solid Lame MU CONSTANT model.  Subsequent constants being ignored.", pd_glob[mn]->MaterialName);
                            declare_warning(msg);   
		            elc_glob[mn]->len_u_mu = floating_point_constant_list_array[0];
		            elc_glob[mn]->lame_mu_model = $5;
		            card_read("PSEUDO_SOLID_LAME_MU_CARD",file_index);			        
		        }
		      }          		                 
		      else
		      {
		        elc_glob[mn]->len_u_mu = read_floats( &(elc_glob[mn]->u_mu), 0);
		        elc_glob[mn]->lame_mu_model = $5;
		        card_read("PSEUDO_SOLID_LAME_MU_CARD",file_index);		  
		      }
                    }
                    else
                    {
                      sprintf(msg, 
	                   "Mesh Motion is not TOTAL_ALE, Pseudo-solid Lame MU card being ignored.");
		      declare_warning(msg);
                    }                  
                  }
                  | error {}
                  ;
                  
pseudo_solid_lame_mu_model:
		  USER_		{$$=USER;}
		| USER_GEN_  	{$$=USER_GEN;}	
		| CONSTANT_	{$$=CONSTANT;}
		| CONTACT_LINE_	{$$=CONTACT_LINE;}
		| error {}
		;
		
pseudo_solid_lame_lambda_card:
                  /* empty */ {}
                | PSEUDO_SOLID_ LAME_ LAMBDA_ EQUALS_ pseudo_solid_lame_lambda_model CR_ /* jjj this needs a floating point list */
                {
                  if( pd_glob[mn]->MeshMotion == TOTAL_ALE)
                  {
                    if ( ( ($5 == USER) || ($5 == USER_GEN) ) && number_of_constants == 0)
                    {
                      sprintf(msg, "Material %s - expected at least 1 constant for the Lame LAMBDA USER or USER_GEN model.", pd_glob[mn]->MaterialName);
                      declare_error(msg);                       
                    }		  
  		    else if($5 == CONSTANT )
  		    {
		      switch(number_of_constants)
		      {
		        case 0:
                          sprintf(msg, "Material %s - expected one constant for the Lame LAMBDA CONSTANT model.", pd_glob[mn]->MaterialName);
                          declare_error(msg);      		          
		          break;
		        case 1:
		          elc_glob[mn]->len_u_lambda = floating_point_constant_list_array[0];
		          elc_glob[mn]->lame_lambda_model = $5;
		          card_read("PSEUDO_SOLID_LAME_LAMBDA_CARD",file_index);	 
		          break;
		        default:
                          sprintf(msg, "Material %s - expected only one constant for the Lame LAMBDA CONSTANT model.  Subsequent constants being ignored.", pd_glob[mn]->MaterialName);
                          declare_warning(msg);   
		          elc_glob[mn]->len_u_lambda = floating_point_constant_list_array[0];
		          elc_glob[mn]->lame_lambda_model = $5;
		          card_read("PSEUDO_SOLID_LAME_LAMBDA_CARD",file_index);			        
		      }      		  
		    }                
		    else
		    {
		      elc_glob[mn]->len_u_lambda = read_floats( &(elc_glob[mn]->u_lambda), 0);
		      elc_glob[mn]->lame_lambda_model = $5;
		      card_read("PSEUDO_SOLID_LAME_LAMBDA_CARD",file_index);		  
		    }                 
                  }
                  else
                  {
                    sprintf(msg, 
	                 "Mesh Motion is not TOTAL_ALE, Pseudo-solid Lame LAMBDA card being ignored.");
		    declare_warning(msg);                  
                  }
                }
                
pseudo_solid_lame_lambda_model:
		  USER_		{$$=USER;}
		| USER_GEN_  	{$$=USER_GEN;}	
		| CONSTANT_	{$$=CONSTANT;}
		| error {}
		;

		
liquid_constitutive_card:
		/* empty */ {} 
		| LIQUID_ CONSTITUTIVE_ EQUATION_ EQUALS_ liquid_consitiutive_model_name CR_ 
		{
		  pd_glob[mn]->MomentumFluxModel = $5;
 		  cr_glob[mn]->MomentumFluxModel = CR_MF_NEWTON_0;
  		  gn_glob[mn]->ConstitutiveEquation = $5;
		  card_read("LIQUID_CONSTITUTIVE_CARD",file_index);
		}
		| error {}
		;
				
liquid_consitiutive_model_name:	
		  CONSTANT_ 		{$$=CONSTANT;}
		| NEWTONIAN_		{$$=NEWTONIAN;}
		| POWERLAW_SUSPENSION_ 	{$$=POWERLAW_SUSPENSION;}
		| POWER_LAW_ 		{$$=POWER_LAW;}
		| CARREAU_ 		{$$=CARREAU;}
		| CARREAU_WLF_ 		{$$=CARREAU_WLF;}
		| CARREAU_SUSPENSION_	{$$=CARREAU_SUSPENSION;}		
		| SUSPENSION_ 		{$$=SUSPENSION;}		
		| BINGHAM_ 		{$$=BINGHAM;}
		| THERMAL_ 		{$$=THERMAL;}
		| CURE_ 		{$$=CURE;}
		| SUSPENSION_ 		{$$=SUSPENSION;}
		| EPOXY_ 		{$$=EPOXY;}
		| FILLED_EPOXY_ 	{$$=FILLED_EPOXY;}
		| error {}
		; 

mechanical_properties_viscosity_card:
		/* empty */ {} 
		|  VISCOSITY_ EQUALS_ viscosity_model_name optional_floating_point_constant_list CR_ 
		{ 
		  if(pd_glob[mn]->MomentumFluxModel == NEWTONIAN)
    		  {
  		    mp_glob[mn]->table = (struct Data_Table *) smalloc (sizeof( struct Data_Table));
		    if ( $3 == CONSTANT )
		    {
		      switch(number_of_constants)
		      {
		        case 0:
                          sprintf(msg, "Material %s - expected one constant for the viscosity CONSTANT model.", pd_glob[mn]->MaterialName);
                          declare_error(msg);  		          
		          break;
		        case 1:
		          /* mp_glob[mn]->viscosity = floating_point_constant_list_array[0]; */
		          break;
		        default:
                          sprintf(msg, "Material %s - expected only one constant for the viscosity CONSTANT model.  Subsequent constants are being ignored.", pd_glob[mn]->MaterialName);
                          declare_warning(msg);  			        
		          /* mp_glob[mn]->viscosity = floating_point_constant_list_array[0];		         */
		      }
		    }
		    else if ( ( ($3 == USER) || ($3 == USER_GEN) ) && number_of_constants == 0)
                    {
                      sprintf(msg, "Material %s - expected at least 1 constant for the Solid Thermal Expansion USER or USER_GEN model.", pd_glob[mn]->MaterialName);
                      declare_error(msg);                       
                    }	    
		    else if ( $3 == FILL && number_of_constants < 2)
		    {
                      sprintf(msg, "Material %s - expected at least 2 constants for the Solid Thermal Expansion FILL model.", pd_glob[mn]->MaterialName);
                      declare_error(msg);		    
		    }
		    else if ($3 == LEVEL_SET && number_of_constants < 3)
		    {
                      sprintf(msg, "Material %s - expected at least 3 constants for the Solid Thermal Expansion LEVEL_SET model.", pd_glob[mn]->MaterialName);
                      declare_error(msg);		    
		    }
		    else if ($3 == SUSPENSION_PM && number_of_constants < 1)
		    {
                      sprintf(msg, "Material %s - expected at least 1 constant for the Solid Thermal Expansion SUSPENSION_PM model.", pd_glob[mn]->MaterialName);
                      declare_error(msg);		    
		    }
		    /* mat_ptr->ViscosityModel = $3;
		    if( &(mat_ptr->u_viscosity) == NULL) 
	              {
	                mat_ptr->u_viscosity = (dbl *)array_alloc(1, number_of_constants, sizeof(dbl));
	              }	
		    mat_ptr->len_u_viscosity = read_floats( &(mp_glob[mn]->u_viscosity) , 0);
		    if ( mat_ptr->u_viscosity[2] == 0.0 && $3 == LEVEL_SET )  mat_ptr->u_viscosity[2] = ls->Length_Scale/2.0;*/
		    card_read("MECHANICAL_PROPERTIES_VISCOSITY_CARD",file_index);
		  }
		}
		| VISCOSITY_ EQUALS_ TABLE_ integer mp_table_independent_variable_name mp_table_interpolation_method optional_file_equals_string CR_
		{
		  if(pd_glob[mn]->MomentumFluxModel == NEWTONIAN)
    		  {
		    mp_glob[mn]->table = (struct Data_Table *) smalloc (sizeof( struct Data_Table));
		    		    
		    card_read("MECHANICAL_PROPERTIES_VISCOSITY_CARD",file_index);
		  }		
		}
		| error {}
		;
		
viscosity_model_name:	
		  CONSTANT_		{$$=CONSTANT;}
                | FILL_   		{$$=FILL;}
		| LEVEL_SET_		{$$=LEVEL_SET;}
		| SUSPENSION_PM_ 	{$$=SUSPENSION_PM;}
		| USER_			{$$=USER;}
		| USER_GEN_		{$$=USER_GEN;}
		| error {}
		;	

mp_table_independent_variable_name:
		  TEMPERATURE_			  {mp_table_species_number = NA;}
		| T_				  {mp_table_species_number = NA;}
		| MASS_FRACTION_ optional_integer {mp_table_species_number = $2;}
		| CAP_PRES_			  {mp_table_species_number = NA;}
		| SPECIES_ optional_integer 	  {mp_table_species_number = $2;}
		| error 			  {}
		;

optional_mp_table_independent_variable_name:
		  /* empty */				{strcpy($$,"NA");}
		| mp_table_independent_variable_name	{strcpy($$,$1);}
		| error 		{}
		;

mp_table_interpolation_method:
		  LINEAR_		{$$ = LINEAR;}
		| BILINEAR_		{$$ = BILINEAR;}
		| error 		{}
		;
		
optional_file_equals_string:
	/* empty */ {strcpy($$,"NA");}
	| FILE_ EQUALS_ STRING_ {strcpy($$,$3);}
	| error {}
	;

low_rate_viscosity_card:
		/* empty */ {} 
		| LOW_ RATE_ VISCOSITY_ EQUALS_ low_rate_viscosity_model_name floating_point_constant_list CR_ 
		{  if(pd_glob[mn]->MomentumFluxModel == POWER_LAW || 
                      pd_glob[mn]->MomentumFluxModel == POWERLAW_SUSPENSION || 
                      pd_glob[mn]->MomentumFluxModel == CARREAU || 
                      pd_glob[mn]->MomentumFluxModel == CARREAU_SUSPENSION || 
                      pd_glob[mn]->MomentumFluxModel == BINGHAM ||
                      pd_glob[mn]->MomentumFluxModel == CARREAU_WLF ||
                      pd_glob[mn]->MomentumFluxModel == SUSPENSION ||
                      pd_glob[mn]->MomentumFluxModel == EPOXY ||
                      pd_glob[mn]->MomentumFluxModel == FILLED_EPOXY ||
                      pd_glob[mn]->MomentumFluxModel == THERMAL ||
                      pd_glob[mn]->MomentumFluxModel == CURE)
		  {
		    switch(number_of_constants)
		    {
		      case 0:
                        sprintf(msg, "Material %s - expected one constant for the Low Rate Viscosity CONSTANT model.", pd_glob[mn]->MaterialName);
                        declare_error(msg);     		          
		        break;
		      case 1:
		        gn_glob[mn]->mu0 = floating_point_constant_list_array[0];
		        gn_glob[mn]->mu0Model = CONSTANT;
		        card_read("LOW_RATE_VISCOSITY_CARD",file_index);	 
		        break;
		      default:
                        sprintf(msg, "Material %s - expected only one constant for the Low Rate Viscosity CONSTANT model.  Subsequent constants being ignored.", pd_glob[mn]->MaterialName);
                        declare_warning(msg);   
		        gn_glob[mn]->mu0 = floating_point_constant_list_array[0];
		        gn_glob[mn]->mu0Model = CONSTANT;  			    
		        card_read("LOW_RATE_VISCOSITY_CARD",file_index);
		    }/* end of switch */
		  }
		  else
		  {
                    sprintf(msg, "Low rate viscosity card is being ingored because liquid constitutive model = %d", pd_glob[mn]->MomentumFluxModel);
                    declare_warning(msg); 		
		  }
                }
		| error {}
		;
		
low_rate_viscosity_model_name:	
		  CONSTANT_ {}
		| error {}
		;

power_law_exponent_equation_card:
		/* empty */ {} 
		| POWER_ LAW_ EXPONENT_ EQUALS_ power_law_exponent_model_name floating_point_constant_list CR_ 
		{  if(pd_glob[mn]->MomentumFluxModel == POWER_LAW || 
                      pd_glob[mn]->MomentumFluxModel == POWERLAW_SUSPENSION  || 
                      pd_glob[mn]->MomentumFluxModel == CARREAU || 
                      pd_glob[mn]->MomentumFluxModel == CARREAU_SUSPENSION || 
                      pd_glob[mn]->MomentumFluxModel == BINGHAM ||
                      pd_glob[mn]->MomentumFluxModel == CARREAU_WLF ||
                      pd_glob[mn]->MomentumFluxModel == SUSPENSION ||
                      pd_glob[mn]->MomentumFluxModel == FILLED_EPOXY)
		  {
		    switch(number_of_constants)
		    {
		      case 0:
                        sprintf(msg, "Material %s - expected one constant for the Power Law Exponent CONSTANT model.", pd_glob[mn]->MaterialName);
                        declare_error(msg);     		          
		        break;
		      case 1:
		        gn_glob[mn]->nexp = floating_point_constant_list_array[0];
		        gn_glob[mn]->nexpModel = CONSTANT;
		        card_read("POWER_LAW_EXPONENT_EQUATION_CARD",file_index);	 
		        break;
		      default:
                        sprintf(msg, "Material %s - expected only one constant for the Power Law Exponent CONSTANT model.  Subsequent constants being ignored.", pd_glob[mn]->MaterialName);
                        declare_warning(msg);   
		        gn_glob[mn]->nexp = floating_point_constant_list_array[0];
		        gn_glob[mn]->nexpModel = CONSTANT;  			    
		        card_read("POWER_LAW_EXPONENT_EQUATION_CARD",file_index);
		    }/* end of switch */
		  }
		  else
		  {
                    sprintf(msg, "Power Law Exponent card is being ingored because liquid constitutive model = %d", pd_glob[mn]->MomentumFluxModel);
                    declare_warning(msg); 		
		  }
                }		 
		| error {}
		;
		
power_law_exponent_model_name:	
		CONSTANT_ ;
		
high_rate_viscosity_card:
		/* empty */ {} 
		| HIGH_ RATE_ VISCOSITY_ EQUALS_ high_rate_viscosity_model_name floating_point_constant_list CR_ 
		{  if( pd_glob[mn]->MomentumFluxModel == CARREAU || 
                       pd_glob[mn]->MomentumFluxModel == CARREAU_SUSPENSION ||
                       pd_glob[mn]->MomentumFluxModel == CARREAU_WLF ||
                       pd_glob[mn]->MomentumFluxModel == BINGHAM)
		   {
		     switch(number_of_constants)
		     {
		       case 0:
                         sprintf(msg, "Material %s - expected one constant for the High Rate Viscosity CONSTANT model.", pd_glob[mn]->MaterialName);
                         declare_error(msg);     		          
		         break;
		       case 1:
		         gn_glob[mn]->muinf = floating_point_constant_list_array[0];
		         gn_glob[mn]->muinfModel = CONSTANT;
		         card_read("HIGH_RATE_VISCOSITY_CARD",file_index);	 
		         break;
		       default:
                         sprintf(msg, "Material %s - expected only one constant for the High Rate Viscosity CONSTANT model.  Subsequent constants being ignored.", pd_glob[mn]->MaterialName);
                         declare_warning(msg);   
		         gn_glob[mn]->muinf = floating_point_constant_list_array[0];
		         gn_glob[mn]->muinfModel = CONSTANT;  			    
		         card_read("HIGH_RATE_VISCOSITY_CARD",file_index);
		     }/* end of switch */
		   }
		   else
		   {
                     sprintf(msg, "High rate viscosity card is being ingored because liquid constitutive model = %d", pd_glob[mn]->MomentumFluxModel);
                     declare_warning(msg); 		
		   }
                }		
		| error {}
		;
		
high_rate_viscosity_model_name:	
		CONSTANT_ ;

time_constant_card:
		/* empty */ {} 
		| TIME_ CONSTANT_ EQUALS_ time_constant_model_name floating_point_constant_list CR_ 
		{  if( pd_glob[mn]->MomentumFluxModel == CARREAU || 
                       pd_glob[mn]->MomentumFluxModel == CARREAU_SUSPENSION ||
                       pd_glob[mn]->MomentumFluxModel == CARREAU_WLF ||
                       pd_glob[mn]->MomentumFluxModel == BINGHAM)
		   {
		     switch(number_of_constants)
		     {
		       case 0:
                         sprintf(msg, "Material %s - expected one constant for the Time Constant CONSTANT model.", pd_glob[mn]->MaterialName);
                         declare_error(msg);     		          
		         break;
		       case 1:
		         gn_glob[mn]->lam = floating_point_constant_list_array[0];
		         gn_glob[mn]->lamModel = CONSTANT;
		         card_read("TIME_CONSTANT_CARD",file_index);	 
		         break;
		       default:
                         sprintf(msg, "Material %s - expected only one constant for the Time Constant CONSTANT model.  Subsequent constants being ignored.", pd_glob[mn]->MaterialName);
                         declare_warning(msg);   
		         gn_glob[mn]->lam = floating_point_constant_list_array[0];
		         gn_glob[mn]->lamModel = CONSTANT;  			    
		         card_read("TIME_CONSTANT_CARD",file_index);
		     }/* end of switch */
		   }
		   else
		   {
                     sprintf(msg, "Time constant card is being ingored because liquid constitutive model = %d", pd_glob[mn]->MomentumFluxModel);
                     declare_warning(msg); 		
		   }
                }		
		| error {}
		;
		
time_constant_model_name:	
		CONSTANT_ ;

aexp_card:	/* empty */ {} 
		| AEXP_ EQUALS_ aexp_model_name floating_point_constant_list  CR_
		{  if( pd_glob[mn]->MomentumFluxModel == CARREAU || 
                       pd_glob[mn]->MomentumFluxModel == CARREAU_SUSPENSION ||
                       pd_glob[mn]->MomentumFluxModel == CARREAU_WLF ||
                       pd_glob[mn]->MomentumFluxModel == BINGHAM)
		   {
		     switch(number_of_constants)
		     {
		       case 0:
                         sprintf(msg, "Material %s - expected one constant for the Aexp CONSTANT model.", pd_glob[mn]->MaterialName);
                         declare_error(msg);     		          
		         break;
		       case 1:
		         gn_glob[mn]->aexp = floating_point_constant_list_array[0];
		         gn_glob[mn]->aexpModel = CONSTANT;
		         card_read("AEXP_CARD",file_index);	 
		         break;
		       default:
                         sprintf(msg, "Material %s - expected only one constant for the Aexp CONSTANT model.  Subsequent constants being ignored.", pd_glob[mn]->MaterialName);
                         declare_warning(msg);   
		         gn_glob[mn]->aexp = floating_point_constant_list_array[0];
		         gn_glob[mn]->aexpModel = CONSTANT;  			    
		         card_read("AEXP_CARD",file_index);
		     }/* end of switch */
		   }
		   else
		   {
                     sprintf(msg, "Aexp card is being ingored because liquid constitutive model = %d", pd_glob[mn]->MomentumFluxModel);
                     declare_warning(msg); 		
		   }
                }		
		| error {}
		;
		
aexp_model_name:	
		CONSTANT_ ;

thermal_exponent_card:
		/* empty */ {} 
		| THERMAL_ EXPONENT_ EQUALS_ thermal_exponent_model_name floating_point_constant_list CR_ 
                {  if( pd_glob[mn]->MomentumFluxModel == BINGHAM || 
                       pd_glob[mn]->MomentumFluxModel == POWERLAW_SUSPENSION ||
                       pd_glob[mn]->MomentumFluxModel == CARREAU_SUSPENSION ||
                       pd_glob[mn]->MomentumFluxModel == CARREAU_WLF ||
                       pd_glob[mn]->MomentumFluxModel == EPOXY ||
                       pd_glob[mn]->MomentumFluxModel == FILLED_EPOXY ||
                       pd_glob[mn]->MomentumFluxModel == THERMAL)
		   {
		     switch(number_of_constants)
		     {
		       case 0:
                         sprintf(msg, "Material %s - expected one constant for the Thermal Exponent CONSTANT model.", pd_glob[mn]->MaterialName);
                         declare_error(msg);     		          
		         break;
		       case 1:
		         gn_glob[mn]->atexp = floating_point_constant_list_array[0];
		         gn_glob[mn]->atexpModel = CONSTANT;
		         card_read("THERMAL_EXPONENT_CARD",file_index);	 
		         break;
		       default:
                         sprintf(msg, "Material %s - expected only one constant for the Thermal Exponent CONSTANT model.  Subsequent constants being ignored.", pd_glob[mn]->MaterialName);
                         declare_warning(msg);   
		         gn_glob[mn]->atexp = floating_point_constant_list_array[0];
		         gn_glob[mn]->atexpModel = CONSTANT;  			    
		         card_read("THERMAL_EXPONENT_CARD",file_index);
		     }/* end of switch */
		   }
		   else
		   {
                     sprintf(msg, "Thermal exponent card is being ingored because liquid constitutive model = %d", pd_glob[mn]->MomentumFluxModel);
                     declare_warning(msg); 		
		   }
                }			
		| error {}
		;
		
thermal_exponent_model_name:	
		CONSTANT_ ;


thermal_wlf_constant2_card:
		/* empty */ {} 
		| THERMAL_ WLF_ CONSTANT2_ EQUALS_ thermal_wlf_constant2_model floating_point_constant_list CR_ 
                {  if( pd_glob[mn]->MomentumFluxModel == CARREAU_WLF)
		   {
		     switch(number_of_constants)
		     {
		       case 0:
                         sprintf(msg, "Material %s - expected one constant for the Thermal Exponent CONSTANT model.", pd_glob[mn]->MaterialName);
                         declare_error(msg);     		          
		         break;
		       case 1:
		         gn_glob[mn]->wlfc2 = floating_point_constant_list_array[0];
		         gn_glob[mn]->wlfc2Model = CONSTANT;
		         card_read("THERMAL_WLF_CONSTANT2_CARD",file_index);	 
		         break;
		       default:
                         sprintf(msg, "Material %s - expected only one constant for the Thermal Exponent CONSTANT model.  Subsequent constants being ignored.", pd_glob[mn]->MaterialName);
                         declare_warning(msg);   
		         gn_glob[mn]->wlfc2 = floating_point_constant_list_array[0];
		         gn_glob[mn]->wlfc2Model = CONSTANT;  			    
		         card_read("THERMAL_WLF_CONSTANT2_CARD",file_index);
		     }/* end of switch */
		   }
		   else
		   {
                     sprintf(msg, "Thermal WLF Constant2 card is being ingored because liquid constitutive model = %d", pd_glob[mn]->MomentumFluxModel);
                     declare_warning(msg); 		
		   }
                }			
		| error {}
		;

thermal_wlf_constant2_model:
		CONSTANT_ {} 
		;

yield_stress_card:
		/* empty */ {} 
		| YEILD_ STRESS_ EQUALS_ yield_stress_model_name floating_point_constant_list CR_ 
                {  if( pd_glob[mn]->MomentumFluxModel == BINGHAM)
		   {
		     switch(number_of_constants)
		     {
		       case 0:
                         sprintf(msg, "Material %s - expected one constant for the Yeild Stress CONSTANT model.", pd_glob[mn]->MaterialName);
                         declare_error(msg);     		          
		         break;
		       case 1:
		         gn_glob[mn]->tau_y = floating_point_constant_list_array[0];
		         gn_glob[mn]->tau_yModel = $4;
		         card_read("YEILD_STRESS_CARD",file_index);	 
		         break;
		       default:
                         sprintf(msg, "Material %s - expected only one constant for the Yeild Stress CONSTANT model.  Subsequent constants being ignored.", pd_glob[mn]->MaterialName);
                         declare_warning(msg);   
		         gn_glob[mn]->tau_y = floating_point_constant_list_array[0];
		         gn_glob[mn]->tau_yModel = $4;  			    
		         card_read("YEILD_STRESS_CARD",file_index);
		     }/* end of switch */
		   }
		   else
		   {
                     sprintf(msg, "Yeild Stress card is being ingored because liquid constitutive model = %d", pd_glob[mn]->MomentumFluxModel);
                     declare_warning(msg); 		
		   }
                }			
		| error {}
		;
		
yield_stress_model_name:	
		CONSTANT_ {$$=CONSTANT;}
		;

yield_exponent_card:
		/* empty */ {} 
		| YEILD_ EXPONENT_ EQUALS_ yield_exponent_model_name floating_point_constant_list CR_ 
                {  if( pd_glob[mn]->MomentumFluxModel == BINGHAM)
		   {
		     switch(number_of_constants)
		     {
		       case 0:
                         sprintf(msg, "Material %s - expected one constant for the Yeild Stress CONSTANT model.", pd_glob[mn]->MaterialName);
                         declare_error(msg);     		          
		         break;
		       case 1:
		         gn_glob[mn]->fexp = floating_point_constant_list_array[0];
		         gn_glob[mn]->fexpModel = $4;
		         card_read("YEILD_EXPONENT_CARD",file_index);	 
		         break;
		       default:
                         sprintf(msg, "Material %s - expected only one constant for the Yeild Stress CONSTANT model.  Subsequent constants being ignored.", pd_glob[mn]->MaterialName);
                         declare_warning(msg);   
		         gn_glob[mn]->fexp = floating_point_constant_list_array[0];
		         gn_glob[mn]->fexpModel = $4;  			    
		         card_read("YEILD_EXPONENT_CARD",file_index);
		     }/* end of switch */
		   }
		   else
		   {
                     sprintf(msg, "Yeild Exponent card is being ingored because liquid constitutive model = %d", pd_glob[mn]->MomentumFluxModel);
                     declare_warning(msg); 		
		   }
                }			
		| error {}
		;
		
yield_exponent_model_name:	
		CONSTANT_ {$$=CONSTANT;}
		;

suspension_maximum_packing_card:
		/* empty */ {} 
		| SUSPENSION_ MAXIMUM_ PACKING_ EQUALS_ suspension_maximum_packing_model_name floating_point_constant_list CR_ 
                {  if( pd_glob[mn]->MomentumFluxModel == SUSPENSION || 
                       pd_glob[mn]->MomentumFluxModel == POWERLAW_SUSPENSION || 
                       pd_glob[mn]->MomentumFluxModel == CARREAU_SUSPENSION || 
                       pd_glob[mn]->MomentumFluxModel == FILLED_EPOXY)
		   {
		     switch(number_of_constants)
		     {
		       case 0:
                         sprintf(msg, "Material %s - expected one constant for the Suspension Maximum Packing CONSTANT model.", pd_glob[mn]->MaterialName);
                         declare_error(msg);     		          
		         break;
		       case 1:
		         gn_glob[mn]->maxpack = floating_point_constant_list_array[0];
		         gn_glob[mn]->maxpackModel = CONSTANT;
		         card_read("SUSPENSION_MAXIMUM_PACKING_CARD",file_index);	 
		         break;
		       default:
                         sprintf(msg, "Material %s - expected only one constant for the Suspension Maximum Packing CONSTANT model.  Subsequent constants being ignored.", pd_glob[mn]->MaterialName);
                         declare_warning(msg);   
		         gn_glob[mn]->maxpack = floating_point_constant_list_array[0];
		         gn_glob[mn]->maxpackModel = CONSTANT;  			    
		         card_read("SUSPENSION_MAXIMUM_PACKING_CARD",file_index);
		     }/* end of switch */
		   }
		   else
		   {
                     sprintf(msg, "Suspension Maximum Packing card is being ingored.");
                     declare_warning(msg); 		
		   }
                }			
		| error {}
		;
		
suspension_maximum_packing_model_name:	
		CONSTANT_ ;

suspension_species_number_card:
		/* empty */ {} 
		| SUSPENSION_ SPECIES_ NUMBER_ EQUALS_ integer CR_ 
		{
                  if( 
                      ( pd_glob[mn]->MomentumFluxModel == SUSPENSION || 
                        pd_glob[mn]->MomentumFluxModel == POWERLAW_SUSPENSION || 
                        pd_glob[mn]->MomentumFluxModel == CARREAU_SUSPENSION || 
                        pd_glob[mn]->MomentumFluxModel == FILLED_EPOXY) 
                      &&
                       (times_card_read("SUSPENSION_SPECIES_NUMBER_CARD",file_index)>0)
                    )
		  {
		    gn_glob[mn]->sus_species_no = $5;
		    card_read("SUSPENSION_SPECIES_NUMBER_CARD",file_index);
		  }
		  else
		  {
                    sprintf(msg, "Suspension Species Number card is being ingored.");
                    declare_warning(msg); 		
		  }		  
		}
		| error {}
		;
		
cure_gel_point_card:
		/* empty */ {} 
		| CURE_ GEL_ POINT_ EQUALS_ cure_gel_point_model_name floating_point_constant_list CR_ 
                {  if( pd_glob[mn]->MomentumFluxModel == CURE || 
                       pd_glob[mn]->MomentumFluxModel == EPOXY || 
                       pd_glob[mn]->MomentumFluxModel == FILLED_EPOXY)
		   {
		     switch(number_of_constants)
		     {
		       case 0:
                         sprintf(msg, "Material %s - expected one constant for the Suspension Maximum Packing CONSTANT model.", pd_glob[mn]->MaterialName);
                         declare_error(msg);     		          
		         break;
		       case 1:
		         gn_glob[mn]->maxpack = floating_point_constant_list_array[0];
		         gn_glob[mn]->maxpackModel = CONSTANT;
		         card_read("CURE_GEL_POINT_CARD",file_index);	 
		         break;
		       default:
                         sprintf(msg, "Material %s - expected only one constant for the Suspension Maximum Packing CONSTANT model.  Subsequent constants being ignored.", pd_glob[mn]->MaterialName);
                         declare_warning(msg);   
		         gn_glob[mn]->maxpack = floating_point_constant_list_array[0];
		         gn_glob[mn]->maxpackModel = CONSTANT;  			    
		         card_read("CURE_GEL_POINT_CARD",file_index);
		     }/* end of switch */
		   }
		   else
		   {
                     sprintf(msg, "Cure Gel Point card is being ingored because liquid constitutive model = %d", pd_glob[mn]->MomentumFluxModel);
                     declare_warning(msg); 		
		   }
                }			
		| error {}
		;
		
cure_gel_point_model_name:	
		CONSTANT_ ;

cure_a_exponent_card:
		/* empty */ {} 
		| CURE_ A_ EXPONENT_ EQUALS_ cure_a_exponent_model_name floating_point_constant_list CR_ 
                {  if( pd_glob[mn]->MomentumFluxModel == CURE || 
                       pd_glob[mn]->MomentumFluxModel == EPOXY || 
                       pd_glob[mn]->MomentumFluxModel == FILLED_EPOXY)
		   {
		     switch(number_of_constants)
		     {
		       case 0:
                         sprintf(msg, "Material %s - expected one constant for the Cure A Exponent CONSTANT model.", pd_glob[mn]->MaterialName);
                         declare_error(msg);     		          
		         break;
		       case 1:
		         gn_glob[mn]->cureaexp = floating_point_constant_list_array[0];
		         gn_glob[mn]->cureaexpModel = CONSTANT;
		         card_read("CURE_A_EXPONENT_CARD",file_index);	 
		         break;
		       default:
                         sprintf(msg, "Material %s - expected only one constant for the Cure A Exponent CONSTANT model.  Subsequent constants being ignored.", pd_glob[mn]->MaterialName);
                         declare_warning(msg);   
		         gn_glob[mn]->cureaexp = floating_point_constant_list_array[0];
		         gn_glob[mn]->cureaexpModel = CONSTANT;  			    
		         card_read("CURE_A_EXPONENT_CARD",file_index);
		     }/* end of switch */
		   }
		   else
		   {
                     sprintf(msg, "Cure Gel Point card is being ingored because liquid constitutive model = %d", pd_glob[mn]->MomentumFluxModel);
                     declare_warning(msg); 		
		   }
                }
		| error {}
		;
		
cure_a_exponent_model_name:	
		CONSTANT_ ;

cure_b_exponent_card:
		/* empty */ {} 
		| CURE_ B_ EXPONENT_ EQUALS_ cure_b_exponent_model_name floating_point_constant_list CR_ 
                {  if( pd_glob[mn]->MomentumFluxModel == CURE || 
                       pd_glob[mn]->MomentumFluxModel == EPOXY || 
                       pd_glob[mn]->MomentumFluxModel == FILLED_EPOXY)
		   {
		     switch(number_of_constants)
		     {
		       case 0:
                         sprintf(msg, "Material %s - expected one constant for the Cure B Exponent CONSTANT model.", pd_glob[mn]->MaterialName);
                         declare_error(msg);     		          
		         break;
		       case 1:
		         gn_glob[mn]->curebexp = floating_point_constant_list_array[0];
		         gn_glob[mn]->curebexpModel = CONSTANT;
		         card_read("CURE_B_EXPONENT_CARD",file_index);	 
		         break;
		       default:
                         sprintf(msg, "Material %s - expected only one constant for the Cure B Exponent CONSTANT model.  Subsequent constants being ignored.", pd_glob[mn]->MaterialName);
                         declare_warning(msg);   
		         gn_glob[mn]->curebexp = floating_point_constant_list_array[0];
		         gn_glob[mn]->curebexpModel = CONSTANT;  			    
		         card_read("CURE_B_EXPONENT_CARD",file_index);
		      }/* end of switch */
		   }
		   else
		   {
                     sprintf(msg, "Cure Gel Point card is being ingored because liquid constitutive model = %d", pd_glob[mn]->MomentumFluxModel);
                     declare_warning(msg); 		
		   }
                }
		| error {}
		;
		
cure_b_exponent_model_name:	
		CONSTANT_ ;

cure_species_number_card:
		/* empty */ {} 
		| CURE_ SPECIES_ NUMBER_ EQUALS_ integer  CR_
		{
                  if( pd_glob[mn]->MomentumFluxModel == CURE || 
                       pd_glob[mn]->MomentumFluxModel == EPOXY || 
                       pd_glob[mn]->MomentumFluxModel == FILLED_EPOXY
                     )
		  {
		    gn_glob[mn]->cure_species_no = $5;
		    card_read("CURE_SPECIES_NUMBER_CARD",file_index);
		  }
		  else
		  {
                    sprintf(msg, "Cure Species Number card is being ingored.");
                    declare_warning(msg); 		  
		  }
		}
		| error {}
		;
		
unreacted_gel_temperature_card:
		/* empty */ {} 
		| UNREACTED_ GEL_ TEMPERATURE_ EQUALS_ CONSTANT_ float  CR_
		{
                  if( pd_glob[mn]->MomentumFluxModel == CURE || 
                       pd_glob[mn]->MomentumFluxModel == EPOXY || 
                       pd_glob[mn]->MomentumFluxModel == FILLED_EPOXY
                     )
		  {
		    gn_glob[mn]->tgel0Model = $5;
		    gn_glob[mn]->tgel0 = $6;
		    card_read("UNREACTED_GEL_TEMPERATURE_CARD",file_index);
		  }
		  else
		  {
                    sprintf(msg, "Unreacted Gel Temperature card is being ingored.");
                    declare_warning(msg); 		  
		  }
		}
		| error {}
		;		
		
polymer_constitutive_equation_card:
		/* empty */ {} 
		| POLYMER_ CONSTITUTIVE_ EQUATION_ EQUALS_ polymer_constitutive_equation_model_name floating_point_constant_list CR_
		| POLYMER_ CONSTITUTIVE_ EQUATION_ EQUALS_ NOPOLYMER_  CR_
		{card_read("POLYMER_CONSTITUTIVE_EQUATION_CARD",file_index);}
		| error {}
		;
		
polymer_constitutive_equation_model_name:	
		NOLDROYDB_ | GIESEKUS_ | WHITE_METZNER_ ;

polymer_stress_formulation_card:
		/* empty */ {} 
		| POLYMER_ STRESS_ FORMULATION_ EQUALS_ polymer_stress_formulation_model_name floating_point_constant_list CR_ 
		{card_read("POLYMER_STRESS_FORMULATION_CARD",file_index);}
		| error {}
		;
		
polymer_stress_formulation_model_name:	
		EVSS_F_ | EVSS_G_ ;

polymer_weight_function_card:
		/* empty */ {} 
		| POLYMER_ WEIGHT_ FUNCTION_ EQUALS_ polymer_weight_function_model_name floating_point_constant_list CR_ 
		{card_read("POLYMER_WEIGHT_FUNCTION_CARD",file_index);}
		| error {}
		;
		
polymer_weight_function_model_name:	
		GALERKIN_ | SUPG_ ;

polymer_weighting_card:
		/* empty */ {} 
		| POLYMER_ WEIGHTING_ EQUALS_ float  CR_
		{card_read("POLYMER_WEIGHTING_CARD",file_index);}
		| error {}
		;
		
polymer_viscosity_card:
		/* empty */ {} 
		| POLYMER_ VISCOSITY_ EQUALS_ polymer_viscosity_model_name float CR_
		{card_read("POLYMER_VISCOSITY_CARD",file_index);}
		| error {}
		;
		
polymer_viscosity_model_name:	
		CONSTANT_ | CARREAU_ | POWER_LAW_ ;	

polymer_time_constant_card:
		/* empty */ {} 
		| POLYMER_ TIME_ CONSTANT_ EQUALS_ polymer_time_constant_model_name float CR_
		{card_read("POLYMER_TIME_CONSTANT_CARD",file_index);}
		| error {}
		;
		
polymer_time_constant_model_name:	
		CONSTANT_ | POWER_LAW_ | CARREAU_ ;

mobility_parameter_card:
		/* empty */ {} 
		| MOBILITY_ PARAMETER_ EQUALS_ mobility_parameter_model_name float CR_
		{card_read("MOBILITY_PARAMETER_CARD",file_index);}
		| error {}
		;
		
mobility_parameter_model_name:	
		CONSTANT_ ;

surface_tension_card:
		/* empty */ {}
		| SURFACE_ TENSION_ EQUALS_ surface_tension_model_name float CR_
		{
		  card_read("SURFACE_TENSION_CARD",file_index);
		}
		| error {}
		;
		
surface_tension_model_name:	
		CONSTANT_ | USER_ | DILATION_ ;
		
solution_temperature_card:	
		/* empty */ {}
		| SOLUTION_ TEMPERATURE_ EQUALS_ solution_temperature_model_name integer float CR_
		{
		  if(pd_glob[mn]->MassFluxModel == FICKIAN_CHARGED)
		  {
		    /* Note: the integer species number is not being used */
		    mat_ptr->SolutionTemperatureModel = $4;
		    card_read("SOLUTION_TEMPERATURE_CARD",file_index);
		  }
		  else
		  {
		    sprintf(msg, "The FICKIAN_CHARGED flux model is not being used.  This Solution Temperature is card being ignored."); 
		    declare_warning(msg);		  
		  }
		}
		| error {}
		;
		
solution_temperature_model_name:	
		  CONSTANT_        {$$=CONSTANT;} 
		| THERMAL_BATTERY_ {$$=THERMAL_BATTERY;}
		| error {}
		;

non_condensable_molecular_weight_card:/* empty */ {}
		| NON_CONDENSABLE_ MOLECULAR_ WEIGHT_ EQUALS_ CONSTANT_ integer float CR_
		{
		  if (read_bc_mp >= 0)
		  {
		    /* Note: The integer is not being used. */
		    mat_ptr->molecular_weight[mat_ptr->Num_Species_Eqn] = $7;		  
		    card_read("NON_CONDENSABLE_MOLECULAR_WEIGHT_CARD",file_index);
		  }
		}
		| error {}
		; 
		
non_volatile_molar_volume_card:/* empty */ {}
		| NON_VOLATILE_ MOLAR_ VOLUME_ EQUALS_ CONSTANT_ integer float CR_
		{
		  if (read_bc_mp >= 0)
		  {
    		    /* Note: The integer is not being used. */
  		    mat_ptr->molar_volume[pd_glob[mn]->Num_Species_Eqn] = $7; 		  
  		    card_read("NON_VOLATILE_MOLAR_VOLUME_CARD",file_index);
  		  }
		}
		| error {}
		; 		

non_volatile_specific_volume_card:/* empty */ {}
		| NON_VOLATILE_ SPECIFIC_ VOLUME_ EQUALS_ CONSTANT_ integer float CR_
		{
		  if (read_bc_mp >= 0)
		  {
		    /* Note: The integer is not being used. */
		    mat_ptr->molar_volume[pd_glob[mn]->Num_Species_Eqn] = $7; 		  
		    card_read("NON_VOLATILE_SPECIFIC_VOLUME_CARD",file_index);
		  }  
		}
		| error {}
		; 		
		
flory_huggins_card:/* empty */ {}
		| FLORY_HUGGINS_ PARAMETERS_ EQUALS_ CONSTANT_ integer integer float CR_
		{
		  if (read_bc_mp >= 0)
		  {
		    n_species = pd_glob[mn]->Num_Species_Eqn + 1;
		    number_of_flory_huggins_parameter_cards_to_accept = ((n_species*n_species - n_species)/2)-1; 
                    mat_ptr->flory_param[$5][$6] = $7; 
	            mat_ptr->flory_param[$6][$5] = $7*mat_ptr->molar_volume[$6]/mat_ptr->molar_volume[$5]; 		  
		    for (k = 0; k<n_species; k++)
		      {
                        mat_ptr->flory_param[k][k] = 0.;
		      }		 		  
		    card_read("FLORY_HUGGINS_CARD",file_index);
		  }
		}
		| error {}
		;

flory_huggins_parameter_card:/* empty */ {}
		| integer integer float CR_
		{
		  if (read_bc_mp >= 0)
		  {
		    if(number_of_flory_huggins_parameter_cards_to_accept>0)
		    {
                      mat_ptr->flory_param[$1][$2] = $3; 
	              mat_ptr->flory_param[$2][$1] = $3*mat_ptr->molar_volume[$2]/mat_ptr->molar_volume[$1]; 
		      card_read("FLORY_HUGGINS_PARAMETER_CARD",file_index);
		      number_of_flory_huggins_parameter_cards_to_accept--;
		    }
		    else
		    {
		      sprintf(msg, "Number of Floury_Huggins parameters has been satisified.  This card being ignored."); 
		      declare_warning(msg);
		    }
		  }
		}
		| error {}
		;


solid_body_source_card:
		/* empty */ {}
		| SOLID_ BODY_ SOURCE_ EQUALS_ CONSTANT_  float float float CR_
		{
		  mp_glob[mn]->MeshSourceModel= CONSTANT;
		  mp_glob[mn]->mesh_source[0]= $6;
		  mp_glob[mn]->mesh_source[1]= $7;
		  mp_glob[mn]->mesh_source[2]= $8;
		  card_read("SOLID_BODY_SOURCE_CARD",file_index);
		}
		| error {}
		;
	
mass_source_card:
		/* empty */ {}
		| MASS_ SOURCE_ EQUALS_ CONSTANT_ float CR_
		{
		  if ($5 == 0.0)
		  {
		    mp_glob[mn]->MassSourceModel = CONSTANT;
		    mat_ptr->mass_source = 0;
		    card_read("MASS_CARD",file_index);
		  }
		  else
		  {
		    declare_error("Invalid floating point parameter.  Only 0.0 is allowed for Mass Source card.");
		  }  
		}
		| error {}
		;	
			
heat_source_card:	
		/* empty */ {}
		| HEAT_ SOURCE_ EQUALS_ heat_source_model optional_floating_point_constant_list CR_
		{
                  mp_glob[mn]->HeatSourceModel = $4;
                  switch($4)
		  {
		    case JOULE:
		      break;
		    case VISC_DISS:
                      mp_glob[mn]->len_u_heat_source = read_floats( &(mp_glob[mn]->u_heat_source), 0); 		      
		      break;
		    case EPOXY:
		      if(number_of_constants == 0)
		      {
                        sprintf(msg, "Material %s - expected at least one constant for the Heat Source EPOXY model.", pd_glob[mn]->MaterialName);
                        declare_error(msg);             		  
		      }
		      else if (number_of_constants > 1)
		      {	      
	      	        sprintf(msg,"Material %s - expected only one constant for the Heat Source EPOXY model.  Subsequent constants are being ignored.", pd_glob[mn]->MaterialName);
	      		declare_warning(msg);	
     	                mp_glob[mn]->heat_source = floating_point_constant_list_array[0];	      		        	        
	              }
	      	      else
	      	      {
     	                mp_glob[mn]->heat_source = floating_point_constant_list_array[0];
	      	      }			      
		      break;
		    case BUTLER_VOLMER:
		      if(number_of_constants == 0)
		      {
                        sprintf(msg, "Material %s - 1 integer and expected 8 floats for the Heat Source BUTLER_VOLMER model.", pd_glob[mn]->MaterialName);
                        declare_error(msg);             		  
		      }
		      else if (number_of_constants > 9)
		      {	      
	      	        sprintf(msg,"Material %s - expected only 1 integer and 8 floats for the Heat Source BUTLER_VOLMER model.  Subsequent constants are being ignored.", pd_glob[mn]->MaterialName);
	      		declare_warning(msg);	
     	                mp_glob[mn]->heat_source = floating_point_constant_list_array[0];	      		        	        
	              }
	      	      else
	      	      {
     	                mp_glob[mn]->heat_source = floating_point_constant_list_array[0];
	      	      }			      
		    case ELECTRODE_KINETICS:
		      break;	
		    case USER:
		      if( &(mat_ptr->u_heat_source) == NULL) 
	              {
	                mat_ptr->u_heat_source = (dbl *)array_alloc(1, number_of_constants, sizeof(dbl));
	              }		    
                      mp_glob[mn]->len_u_heat_source = read_floats( &(mp_glob[mn]->u_heat_source), 0);
		      break;
		    case CONSTANT:
		      if(number_of_constants == 0)
		      {
                        sprintf(msg, "Material %s - expected at one constant for the Heat Source CONSTANT model.", pd_glob[mn]->MaterialName);
                        declare_error(msg);             		  
		      }
		      else if(number_of_constants > 1)
		      {	      
	      	        sprintf(msg,"Material %s - expected only one constant for the Heat Source CONSTANT model.  Subsequent constants are being ignored.", pd_glob[mn]->MaterialName);
	      		declare_warning(msg);	
     	                mp_glob[mn]->heat_source = floating_point_constant_list_array[0];	      		        	        
	      	      }
	              else
	              {
     	                mp_glob[mn]->heat_source = floating_point_constant_list_array[0];
	              }		      
	              break;		        	        
	            default:   
	              break;   
	          }  /* end switch($4) */			
		  card_read("HEAT_SOURCE_CARD",file_index);
		}
		| error {}
		;			
		
heat_source_model:
		  JOULE_		{$$=JOULE;}
		| VISC_DISS_		{$$=VISC_DISS;}
		| EPOXY_		{$$=EPOXY;}
		| BUTLER_VOLMER_	{$$=BUTLER_VOLMER;}
		| ELECTRODE_KINETICS_	{$$=ELECTRODE_KINETICS;}
		| CONSTANT_		{$$=CONSTANT;}
		| USER_			{$$=USER;}
		| error {}
		;		
				
navier_stokes_source_card:		
		/* empty */ {}
		| NAVIER_STOKES_ SOURCE_ EQUALS_ navier_stokes_source_model floating_point_constant_list CR_
		{
		  mat_ptr->MomentumSourceModel = $4;
                  switch($4)
		  {
		    case BOUSS:
		    case BOUSSINESQ:
		    case FILL_SRC:
		    case LEVEL_SET:
		    case SUSPENSION_PM:
		      if (number_of_constants < 3)
		      {
		        sprintf(msg, "Material %s - expected at least three constants for this Navier-Stokes model.", pd_glob[mn]->MaterialName);
                        declare_error(msg);   
		      }
		      else
		      {
			mat_ptr->momentum_source[0] = floating_point_constant_list_array[0];
			mat_ptr->momentum_source[1] = floating_point_constant_list_array[1];
          		mat_ptr->momentum_source[2] = floating_point_constant_list_array[2];		      	        
		      }
		      if (number_of_constants >3)
		      {
			sprintf(msg,"Material %s - expected only three constants for this Navier-Stokes model.  Subsequent constants are being ignored.", pd_glob[mn]->MaterialName);
	      		declare_warning(msg);		      
		      }	     
		      if( (mat_ptr->DensityModel != SUSPENSION) && ($4==SUSPENSION_PM) )
	              {
	                sprintf(msg, 
			   "For matl %s, %s = \"%s\" needs %s = \"%s\".\n",
			   pd_glob[mn]->MaterialName,
			   "Navier-Stokes Source", "SUSPENSION_PM",
			   "Density", "SUSPENSION_PM");
	                declare_error(msg);
	              }			    
		      break;
		    case BOUSS_JXB:
		    case SUSPEND:
		    case SUSPENSION:
		      if (number_of_constants < 4)
		      {
		        sprintf(msg, "Material %s - expected at least three constants for this Navier-Stokes Source model.", pd_glob[mn]->MaterialName);
                        declare_error(msg);   
		      }
		      else
		      {
	                mat_ptr->len_u_momentum_source = 1;
                        mat_ptr->momentum_source[0] = floating_point_constant_list_array[0];
		        mat_ptr->momentum_source[1] = floating_point_constant_list_array[1];
          	        mat_ptr->momentum_source[2] = floating_point_constant_list_array[2];  
		        mat_ptr->u_momentum_source  = (dbl *)array_alloc(1, 1, sizeof(dbl));          	        
                        mat_ptr->u_momentum_source[0]  = floating_point_constant_list_array[3]; 		      	        
		      }
		      if (number_of_constants >3)
		      {
			sprintf(msg,"Material %s - expected only three constants for this Navier-Stokes model.  Subsequent constants are being ignored.", pd_glob[mn]->MaterialName);
	      		declare_warning(msg);		      
		      }	  
		      if( (mat_ptr->DensityModel != SUSPENSION) && ($4==SUSPENSION) )
	              {
	                sprintf(msg, 
			   "For matl %s, %s = \"%s\" needs %s = \"%s\".\n",
			   pd_glob[mn]->MaterialName,
			   "Navier-Stokes Source", "SUSPEND",
			   "Density", "SUSPENSION");
	                declare_error(msg);
	              }		    	    
		      break;	                      
		     	      		      		      	
		    case USER: 
		      if( &(mat_ptr->u_momentum_source) == NULL) 
	              {
	                mat_ptr->u_momentum_source = (dbl *)array_alloc(1, number_of_constants, sizeof(dbl));
	              }
 		      mat_ptr->len_u_momentum_source = read_floats( &(mat_ptr->u_momentum_source), 0);		    
		      break;
		    case CONSTANT:
		      if (number_of_constants < 3)
		      {
		        sprintf(msg, "Material %s - expected at least three constants for the Navier-Stokes Source CONSTANT model.", pd_glob[mn]->MaterialName);
                        declare_error(msg);   
		      }
		      else
		      {
	                mat_ptr->momentum_source[0] = floating_point_constant_list_array[0];  
	                mat_ptr->momentum_source[1] = floating_point_constant_list_array[1];   
	                mat_ptr->momentum_source[2] = floating_point_constant_list_array[2];  		        
		      }
		      if (number_of_constants >3)
		      {
			sprintf(msg,"Material %s - expected only three constants for the Navier-Stokes Source CONSTANT model.  Subsequent constants are being ignored.", pd_glob[mn]->MaterialName);
	      		declare_warning(msg);		      
		      }	      
	              break;		        	        
	            default:   
	              break;   
	          }  /* end switch($4) */	
		  card_read("NAVIER_STOKES_SOURCE_CARD",file_index);
		}
		| error {}
		;
		
navier_stokes_source_model:
		  BOUSS_		{$$=BOUSS;}
		| FILL_			{$$=FILL_SRC;}
		| LEVEL_SET_		{$$=LEVEL_SET;}
		| SUSPEND_		{$$=SUSPEND;}
		| SUSPENSION_		{$$=SUSPENSION;}
		| SUSPENSION_PM_	{$$=SUSPENSION_PM;}
		| BOUSS_JXB_		{$$=BOUSS_JXB;}		
		| BOUSSINESQ_		{$$=BOUSSINESQ;}
		| CONSTANT_		{$$=CONSTANT;}	
		| USER_			{$$=USER;}						
		| error {}
		;				
			
species_source_card:	
		/* empty */ {}
		| SPECIES_ SOURCE_ EQUALS_ integer species_source_model floating_point_constant_list CR_
		{ /* the correct implementation of this card required CHEMKIN */
		  if(number_of_species_source_cards_needed[mn] > 0)
		  {
		    number_of_species_source_cards_needed[mn]--;
		    fallback_chemkin_generic_prop(&model_read, $4, &(mp_glob[mn]->SpeciesSourceModel[i]), FALSE, mat_ptr );
		    /* NOTE: this chemkin routine is not the one in mm_input_mp.c...it's a copy of that routine 
		       and can be found at the bottom of this file.  */
		    /*fallback_chemkin_generic_prop($5, i, mp_glob[mn]->SpeciesSourceModel[i], FALSE, mp_glob[mn]);*/
                    switch($5)
		    {
		      case EPOXY:
		        if(number_of_constants < 6)
		        {
                          sprintf(msg, "Material %s - expected at least 6 constants for the Species Source EPOXY model.", pd_glob[mn]->MaterialName);
                          declare_error(msg);             		  
		        }
		        else
		        {
		          mp_glob[mn]->SpeciesSourceModel[i] = $5;
	                  mp_glob[mn]->u_species_source[i] = (dbl *)array_alloc(1,6,sizeof(dbl)); 
	                  mp_glob[mn]->len_u_species_source[i] = 6; 
	                  mp_glob[mn]->u_species_source[i][0] = floating_point_constant_list_array[0];  /* prefactor for k1 */
	                  mp_glob[mn]->u_species_source[i][1] = floating_point_constant_list_array[1];  /* exponent for k1 */ 
	                  mp_glob[mn]->u_species_source[i][2] = floating_point_constant_list_array[2];  /* prefactor for k2 */
	                  mp_glob[mn]->u_species_source[i][3] = floating_point_constant_list_array[3];  /* exponent for k2 */
	                  mp_glob[mn]->u_species_source[i][4] = floating_point_constant_list_array[4];  /* m power law coefficient */
	                  mp_glob[mn]->u_species_source[i][5] = floating_point_constant_list_array[5];  /* n power law coefficient */
	                  if(number_of_constants > 6)
	                   {
                             sprintf(msg, "Material %s - accepected only 6 constants for the Species Source EPOXY model.  Subsequent constants being ignored.", pd_glob[mn]->MaterialName);
                             declare_warning(msg);	                  
	                   }		        		        
		        } 
		        break;
		      case EPOXY_DEA:
		        if (number_of_constants < 5 )
		        {
                          sprintf(msg, "Material %s - expected at least 5 constants for the Species Source EPOXY_DEA model.", pd_glob[mn]->MaterialName);
                          declare_error(msg);               		  
		        }
		        else   		          
		        {	    
                          mp_glob[mn]->SpeciesSourceModel[i] = $5;
                          mp_glob[mn]->u_species_source[i] = (dbl *) array_alloc(1,5,sizeof(dbl));
                          mp_glob[mn]->len_u_species_source[i] = 5;
                          mp_glob[mn]->u_species_source[i][0] = floating_point_constant_list_array[0];  /* prefactor for k1 */
                          mp_glob[mn]->u_species_source[i][1] = floating_point_constant_list_array[1];  /* exponent for k1 */ 
                          mp_glob[mn]->u_species_source[i][2] = floating_point_constant_list_array[2];  /* prefactor for k2, low T */
                          mp_glob[mn]->u_species_source[i][3] = floating_point_constant_list_array[3];  /* exponent for k2, low T */
                          mp_glob[mn]->u_species_source[i][4] = floating_point_constant_list_array[4];  /* prefactor for k2, mid T */
	                  if(number_of_constants > 6)
	                  {
                            sprintf(msg, "Material %s - accepected only 5 constants for the Species Source EPOXY_DEA model.  Subsequent constants being ignored.", pd_glob[mn]->MaterialName);
                            declare_warning(msg);
	                  }
	                }
		        break;
		      case BUTLER_VOLMER:
		        if(number_of_constants < 9)
		        {
                          sprintf(msg, "Material %s - expected at least 9 floats for the Species Source BUTLER_VOLMER model.", pd_glob[mn]->MaterialName);
                          declare_error(msg);
		        }
		        else
		        {
		          mp_glob[mn]->SpeciesSourceModel[i] = $5;
	                  mp_glob[mn]->u_species_source[i] = (dbl *)array_alloc(1,9,sizeof(dbl));
	                  mp_glob[mn]->len_u_species_source[i] = 9;
	                  mp_glob[mn]->u_species_source[i][0] = floating_point_constant_list_array[0];  /* stoichiometric coefficient */
	                  mp_glob[mn]->u_species_source[i][1] = floating_point_constant_list_array[1];  /* product of interfacial area by exchange current density, ai0  */ 
	                  mp_glob[mn]->u_species_source[i][2] = floating_point_constant_list_array[2];  /* reaction order, beta */
	                  mp_glob[mn]->u_species_source[i][3] = floating_point_constant_list_array[3];  /* reference concentration, cref */
	                  mp_glob[mn]->u_species_source[i][4] = floating_point_constant_list_array[4];  /* anodic transfer coeficient, aa */
	                  mp_glob[mn]->u_species_source[i][5] = floating_point_constant_list_array[5];  /* cathodic transfer coefficient, ac */
	                  mp_glob[mn]->u_species_source[i][6] = floating_point_constant_list_array[6];  /* temperature, T */
	                  mp_glob[mn]->u_species_source[i][7] = floating_point_constant_list_array[7];  /* open-circuit potential, V */
	                  mp_glob[mn]->u_species_source[i][8] = floating_point_constant_list_array[8];  /* number of electrons, n */
	                  if(number_of_constants > 9)
	                   {
                             sprintf(msg, "Material %s - accepected only 9 floats for the Species Source BUTLER_VOLMER model.  Subsequent constants being ignored.", pd_glob[mn]->MaterialName);
                             declare_warning(msg);
	                   }
		        }
		        break;
		      case ELECTROOSMOTIC:
		        if(number_of_constants < 11)
		        {
                          sprintf(msg, "Material %s - expected at least 11 floats for the Species Source BUTLER_VOLMER model.", pd_glob[mn]->MaterialName);
                          declare_error(msg);
		        }
		        else
		        {
		          mp_glob[mn]->SpeciesSourceModel[i] = $5;
	                  mp_glob[mn]->u_species_source[i] = (dbl *)array_alloc(1,11,sizeof(dbl));
	                  mp_glob[mn]->len_u_species_source[i] = 11;
	                  mp_glob[mn]->u_species_source[i][0] = floating_point_constant_list_array[0];   /* index of species involved in rxn */
	                  mp_glob[mn]->u_species_source[i][1] = floating_point_constant_list_array[1];   /* stoichiometric coefficient */
	                  mp_glob[mn]->u_species_source[i][2] = floating_point_constant_list_array[2];   /* product of interfacial area by exchange current density, ai0  */ 
	                  mp_glob[mn]->u_species_source[i][3] = floating_point_constant_list_array[3];   /* reaction order, beta */
	                  mp_glob[mn]->u_species_source[i][4] = floating_point_constant_list_array[4];   /* reference concentration, cref */
	                  mp_glob[mn]->u_species_source[i][5] = floating_point_constant_list_array[5];   /* anodic transfer coeficient, aa */
	                  mp_glob[mn]->u_species_source[i][6] = floating_point_constant_list_array[6];   /* cathodic transfer coefficient, ac */
	                  mp_glob[mn]->u_species_source[i][7] = floating_point_constant_list_array[7];   /* temperature, T */
	                  mp_glob[mn]->u_species_source[i][8] = floating_point_constant_list_array[8];   /* open-circuit potential, V */
	                  mp_glob[mn]->u_species_source[i][9] = floating_point_constant_list_array[9];   /* number of electrons, n */
	                  mp_glob[mn]->u_species_source[i][10] = floating_point_constant_list_array[10]; /* electroosmotic drag coeff., nd */
	                  if(number_of_constants > 11)
	                   {
                             sprintf(msg, "Material %s - accepected only 11 floats for the Species Source BUTLER_VOLMER model.  Subsequent constants being ignored.", pd_glob[mn]->MaterialName);
                             declare_warning(msg);
	                   }
		        }
		        break;
		      case ELECTRODE_KINETICS:
		        break;	
		      case USER:
		        mp_glob[mn]->SpeciesSourceModel[i] = $5;  
/*jjj special allocation for species number array here */
 		        mp_glob[mn]->len_u_species_source[i] = read_floats( &(mp_glob[mn]->u_species_source), 0);
		        break;
		      case CONSTANT:
		        mp_glob[mn]->SpeciesSourceModel[i] = $5; 		       
                        mp_glob[mn]->species_source[i] = floating_point_constant_list_array[0];
		        break;		        		        	        
		      default:   
		        break;   
		    }  /* end switch($4) */			
		  }
		  else 
		  {
                    sprintf(msg, "Material %s already has %d Species Source cards.  This Species Source card is being ignored.",pd_glob[mn]->MaterialName,pd_glob[mn]->Num_Species);
                    declare_warning(msg);		  
		  }
		  card_read("SPECIES_SOURCE_CARD",file_index);
		}
		| error {}
		;			
		
species_source_model:
		  EPOXY_		{$$=EPOXY;}
		| EPOXY_DEA_		{$$=EPOXY_DEA;}
                | BUTLER_VOLMER_	{$$=BUTLER_VOLMER;}
		| ELECTRODE_KINETICS_	{$$=ELECTRODE_KINETICS;}
                | ELECTROOSMOTIC_	{$$=ELECTROOSMOTIC;}
		| CONSTANT_		{$$=CONSTANT;}
		| USER_			{$$=USER;}
		| error {}
		;
		
current_source_card:
		/* empty */ {}
		| CURRENT_ SOURCE_ EQUALS_ current_source_model optional_floating_point_constant_list CR_
		{
                  mp_glob[mn]->CurrentSourceModel = $4;
                  switch($4)
		  {
		    case BUTLER_VOLMER:
		      if(number_of_constants == 0)
		      {
                        sprintf(msg, "Material %s - expected 1 integer and 8 floats for the Current Source BUTLER_VOLMER model.", pd_glob[mn]->MaterialName);
                        declare_error(msg);             		  
		      }
		      else if (number_of_constants > 9)
		      {	      
	      	        sprintf(msg,"Material %s - expected only 1 integer and 8 floats for the Current Source BUTLER_VOLMER model.  Subsequent constants are being ignored.", pd_glob[mn]->MaterialName);
	      		declare_warning(msg);	
     	                mp_glob[mn]->current_source = floating_point_constant_list_array[0];	      		        	        
	              }
	      	      else
	      	      {
     	                mp_glob[mn]->current_source = floating_point_constant_list_array[0];
	      	      }			      
		    case ELECTRODE_KINETICS:
		      break;	
		    case USER:
		      if( &(mat_ptr->u_current_source) == NULL) 
	              {
	                mat_ptr->u_current_source = (dbl *)array_alloc(1, number_of_constants, sizeof(dbl));
	              }		    
                      mp_glob[mn]->len_u_current_source = read_floats( &(mp_glob[mn]->u_current_source), 0);
		      break;
		    case CONSTANT:
		      if(number_of_constants == 0)
		      {
                        sprintf(msg, "Material %s - expected at one constant for the Current Source CONSTANT model.", pd_glob[mn]->MaterialName);
                        declare_error(msg);             		  
		      }
		      else if(number_of_constants > 1)
		      {	      
	      	        sprintf(msg,"Material %s - expected only one constant for the Current Source CONSTANT model.  Subsequent constants are being ignored.", pd_glob[mn]->MaterialName);
	      		declare_warning(msg);	
     	                mp_glob[mn]->current_source = floating_point_constant_list_array[0];	      		        	        
	      	      }
	              else
	              {
     	                mp_glob[mn]->current_source = floating_point_constant_list_array[0];
	              }		      
	              break;		        	        
	            default:   
	              break;   
	          }  /* end switch($4) */			
		  card_read("CURRENT_SOURCE_CARD",file_index);
		}
		| error {}
                ;

current_source_model:
		  BUTLER_VOLMER_		{$$=BUTLER_VOLMER;}
		  ELECTRODE_KINETICS_		{$$=ELECTRODE_KINETICS;}
		| FICKIAN_CHARGED_		{$$=FICKIAN_CHARGED;}
		| STEFAN_MAXWELL_CHARGED_	{$$=STEFAN_MAXWELL_CHARGED;}
		| CONSTANT_			{$$=CONSTANT;}
		| error {}
                ;
      
material_variable_initialize_card:
		/* empty */ {}
		| INITIALIZE_ EQUALS_ var_name integer float CR_
		{
		  Var_init_mat[mn][Num_Var_Init_Mat[mn]].var = pp_index;
		  Var_init_mat[mn][Num_Var_Init_Mat[mn]].ktype = $4;
		  Var_init_mat[mn][Num_Var_Init_Mat[mn]].init_val = $5;
		  Num_Var_Init_Mat[mn]++;
		  if (Num_Var_Init_Mat[mn] > 0) 
                  {
                    fprintf(stderr, "Material %d Number of reinitialized vars  = %d\n", 
	            mn, Num_Var_Init_Mat[mn]);
                  }
		  card_read("MATERIAL_VARIABLE_INITIALIZE_CARD",file_index);
		}
		| error {}
		;
			    		
example_card:					/*EXAMPLE CARD RULE AND ACTION*/
		/* empty */ {}
		| EXAMPLE_TOKEN1_ EXAMPLE_TOKEN2_  EQUALS_  integer float CR_
		{
		  /* Put the desired action here. */
		  card_read("EXAMPLE_CARD",file_index);
		}
		| error {}
		;	
		
%%


extern void initialize_Boundary_Condition (struct Boundary_Condition *);
extern struct BC_descriptions * alloc_BC_description (struct BC_descriptions *);
int number_of_characters_in_this_line;
enum mode parser_file_mode;
int i;
int num_new_BC_Desc = 0;
void defaults(int);
void consistency(int);
int  file_to_string (FILE *file, char string[], int max_number_chars_in);
void initialize_input_file_values();
void initialize_material_file_values(int);
void default_it(char *,int);


parse_input_file()
{
  if (preprocessor(file_index)) /* Perform preprocessing necessary to parse the Goma input file */
  {
    yyparse();			/* Parse the Goma input file */
    defaults(file_index); 
    consistency(file_index); 
  }  
  for(mn=0;mn<number_of_materials_found;mn++)
    {
      file_index++;
      if(preprocessor(file_index)) /* Perform preprocessing necessary to parse the Goma material file */
      {
        yyparse();		/* Parse the Goma material file */
        defaults(file_index);
        consistency(file_index); 
      }
    }
  summary_and_closes();
  if (error_number > 0) exit;
  exit;
  return(0);  
}

int preprocessor
  (
  int file_index
  )
{
  int material_index;
  if (ProcID == 0)
  {	/* Processor zero processing. */
    if( file_index == 0 )  /* input file processing */ 
    {
      /* Open the error & parser log file */
      if ( (error_log = fopen( "error.log", "w")) == NULL)
        {
          fprintf(stderr, "preprocessor: Open of error.log failed.\n");
          return(0);		
        } 
        else 
        {
          curtime = time(NULL);
          loctime = localtime(&curtime);

fprintf(error_log,"\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
          fprintf(error_log,"+                                                              \n"); 
          fprintf(error_log,"+                GOMA Parser Error Log File                    \n");
          fprintf(error_log,"+                                                              \n");
          fprintf(error_log,"+                 %s",asctime(loctime) );
          fprintf(error_log,"+                                                              \n");
          fprintf(error_log,"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n");
        }
      if ( (parser_log = fopen( "parser.log", "w")) == NULL)
        {
          fprintf(stderr, "preprocessor: Open of parser.log failed.\n");
          return(0);		
        } 
        else 
        {
          curtime = time(NULL);
          loctime = localtime(&curtime);


fprintf(parser_log,"\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
          fprintf(parser_log,"+                                                              \n"); 
          fprintf(parser_log,"+                   GOMA Parser Log File                       \n");
          fprintf(parser_log,"+                                                              \n");
          fprintf(parser_log,"+                 %s",asctime(loctime) );
          fprintf(parser_log,"+                                                              \n");
          fprintf(parser_log,"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n");
          mystring=(char *)malloc(100000*sizeof(char));

          if(!mystring)
          {
            printf("Could not allocate space for big string, mystring.\n");
            exit(-1);
          }
          /* Open the input file */
          if((in_file = fopen_aprepro( Input_File, "r")) == NULL)
          {
            sprintf(msg,"<<< ERROR %i: Open of %s file failed. >>>\n\n",error_number++,Input_File);
            fprintf( parser_log,"%s",msg);
            fprintf( error_log,"%s",msg);
            fprintf( stdout,"%s",msg);
            return(0);
          } 
          else 
          {
            file_to_string(in_file, mystring, 10000);
            fprintf(parser_log,"Parse of file \"%s\":\n",Input_File);
            fprintf(error_log,"Parsing file \"%s\"\n",Input_File); 
            fprintf(stdout,"Parsing file \"%s\"\n",Input_File);         
            fprintf(parser_log, "%d ", line_number);
            initialize_input_file_values();
            parser_file_mode = input_file; /* input or mat */
            return(1);
          } 
        } 
      }
      else
      {  /* material file processing */
        material_index = file_index-1;
        mystring=(char *)malloc(10000*sizeof(char));
        if(!mystring)
        {
          printf("Could not allocate space for big string, mystring.\n");
          exit(-1);
        }
        if((in_file = fopen_aprepro( strcat(mp_glob[material_index]->Material_Name,".mat"), "r")) == NULL)
        {
          sprintf(msg,"\n<<< ERROR %i: Open of %s file failed. >>>",error_number++, mp_glob[material_index]->Material_Name );
          fprintf( parser_log,"%s",msg);
          fprintf( error_log,"%s",msg);          
          fprintf( stdout,"%s",msg);          
          return(0);
        } 
        else 
        {
          file_to_string(in_file, mystring, 10000);
          fprintf(parser_log,"\n\nParse of file \"%s\"\n",mp_glob[material_index]->Material_Name);
          fprintf(error_log,"\nParse of file \"%s\"\n",mp_glob[material_index]->Material_Name); 
          fprintf(stdout,"\nParsing file \"%s\"\n",mp_glob[material_index]->Material_Name);          
          line_number=1;
          fprintf(parser_log, "%d ", line_number);
          initialize_material_file_values(material_index);
          parser_file_mode = mat_file; /* input or mat */
          return(1);        
        }
      } 
    } 
    else 
    { /* Non-Processor zero processing  */
      /* The details of this processing will have to be worked out with Phil, et.al.  */
      if (file_index == 0)
      { /* non-procesor 0 input file preprocessing */
        /* Recieve the string from the input file.  If you don't get it, do some error processing and exit */
        initialize_input_file_values();
        parser_file_mode = input_file; /* input or mat */
        return(1);
      } 
      else 
      { /* non-processor 0 material file preprocessing */
        /* Receive the string from processor 0 for the material file.   If you don't get it, do some error processing and exit */
        initialize_material_file_values(material_index);
        parser_file_mode = mat_file; /* input or mat */
        return(1);             
      }  
    }
}

void defaults
  (
  int file_index
  )
{

  int k = 0;
  enum file_types file_type_being_parsed;
  
  /*  The data structures, etc. involved are defined in mm_cards.h as:										*/
  /*   																		*/
  /* enum file_types {MAT, INPUT};														*/
  /* enum card_types {REQUIRED_ONE_PER_FILE, REQUIRED_MN_PER_FILE, REQUIRED_SN_PER_FILE, OPTIONAL_DEFAULTABLE, OPTIONAL_DEFAULTABLE_MN,         */
  /* OPTIONAL_NOT_DEFAULTABLE};	                                                                                                                */
  /*																		*/
  /* struct card_characteristics 														*/
  /*     char * card_name;															*/
  /*     enum file_types file_type;														*/
  /*     enum card_types card_type;     													*/
  /*     int times_read_per_file[12];														*/
  /*																		*/
  /*  struct card_characteristics cards[]	= {........etc.}										*/
  /*																		*/
  /*  int num_cards = sizeof(cards) / sizeof(struct card_characteristics);									*/

  if (file_index == 0)
  {
    file_type_being_parsed = INPUT;
  }
  else
  {
    file_type_being_parsed = MAT;
  }  
  
  while (k<num_cards)
  {
    if( file_type_being_parsed == cards[k].file_type )
    {
    
 /* fprintf (stdout,"\n %d  %s",file_index, cards[k].card_name);     */

      if ( cards[k].card_type == REQUIRED_ONE_PER_FILE && cards[k].times_read_per_file[file_index] != 1 )
      {
          sprintf(msg, "Wrong number of %s cards found, %d expected, %d found.",cards[k].card_name, 1, cards[k].times_read_per_file[file_index]);
          declare_error(msg);
      }

      if ( cards[k].card_type == REQUIRED_MN_PER_FILE && cards[k].times_read_per_file[file_index] != number_of_materials_found )
      {
          sprintf(msg, "Wrong number of %s cards found, %d expected, %d found.",cards[k].card_name, number_of_materials_found, cards[k].times_read_per_file[file_index]);
          declare_error(msg); 
      }
      if ( cards[k].card_type == REQUIRED_SN_PER_FILE && cards[k].times_read_per_file[file_index] != mp_glob[mn]->Num_Species )
      {
          sprintf(msg, "Wrong number of %s cards found, %d expected, %d found.",cards[k].card_name, mp_glob[mn]->Num_Species, cards[k].times_read_per_file[file_index]);
          declare_error(msg); 
      }
     
      if (cards[k].card_type == OPTIONAL_DEFAULTABLE && cards[k].times_read_per_file[file_index] == 0)
      { 
        default_it( cards[k].card_name, mn );
      }  
 
      if (cards[k].card_type == OPTIONAL_DEFAULTABLE_MN && cards[k].times_read_per_file[file_index] != number_of_materials_found)
      { 
        default_it( cards[k].card_name, mn );
      }            
      
    }  /* end file_type_being_parsed */    
    k++;
  } /* end while */ 
}

void default_it
  (
  char * card_name, 
  int mn
  )
{

  int i;
  
  if ( !strcmp(card_name,"AEXP_CARD") )           
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"ALC_DESIRED_SOLUTION_FRACTION_CARD") )       
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"ALC_MAX_PARAMETER_SENSITIVITY_CARD") )    
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"ALC_TANGENT_FACTOR_EXPONENT_CARD") )           
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }  
  
  if ( !strcmp(card_name,"ALC_TANGENT_FACTOR_STEP_LIMIT_CARD") )        
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }  
  
  if ( !strcmp(card_name,"ANNEAL_MESH_ON_OUTPUT_CARD") )          
  {  sprintf(msg, "Defaulting %s, Anneal_Mesh = FALSE",card_name);
     Anneal_Mesh = FALSE; 
     declare_default_warning(msg); return; }  
  
  if ( !strcmp(card_name,"AUGMENTING_CONDITIONS_INITIAL_GUESS_CARD") )      
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }  
  
  if ( !strcmp(card_name,"AUGMENTING_CONDITIONS_SPECIFICATIONS_SECTION_DELIMITER_CARD") )             
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }  
  
  if ( !strcmp(card_name,"AUGMENTING_CONDITION_CARD") )           
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }  
  
  if ( !strcmp(card_name,"BOUNDARY_CONDITION_DATA_FLOAT_TAG_CARD") )           
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }  
  
  if ( !strcmp(card_name,"BOUNDARY_CONDITION_ID_CARD") )            
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }  
  
  if ( !strcmp(card_name,"BOUNDARY_CONDITION_SPECIFICATIONS_DELIMITER_CARD") )            
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; } 
  
  if ( !strcmp(card_name,"BULK_SPECIES_CARD") )           
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }  
  
  if ( !strcmp(card_name,"CAPILLARY_PRESSURE_IN_POROUS_MEDIA_CARD") )           		  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }  
  
  if ( !strcmp(card_name,"CC_CARD") )           
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }  
  
  if ( !strcmp(card_name,"CONCENTRATION_CONTOURS_CARD") )           		  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }  
  
  if ( !strcmp(card_name,"CONTINUATION_CARD") )           
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }  
  
  if ( !strcmp(card_name,"CONTINUATION_ORDER_CARD") )           
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }  
  
  if ( !strcmp(card_name,"CONTINUATION_PRINTING_FREQUENCY_CARD") )           
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }  
  
  if ( !strcmp(card_name,"CONTINUATION_SPECIFICATIONS_DELIMITER_CARD") )            
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }  
  
  if ( !strcmp(card_name,"CONTINUATION_TYPE_CARD") )           
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }  
  
  if ( !strcmp(card_name,"CONVECTIVE_LAGRANGIAN_VELOCITY_CARD") )            
  {  pd_glob[mn]->MeshInertia = 0;
     sprintf(msg, "Defaulting %s, pd_glob[mn]->MeshInertia = NONE",card_name); 
     declare_default_warning(msg); return; }  
  
  if ( !strcmp(card_name,"COORDINATE_SYSTEM_CARD") )            
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }  
  
  if ( !strcmp(card_name,"CROSSSTREAM_SHEAR_RATE_CARD") )            		  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }  
  
  if ( !strcmp(card_name,"CURE_A_EXPONENT_CARD") )            
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }  
  
  if ( !strcmp(card_name,"CURE_B_EXPONENT_CARD") )           
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }  
  
  if ( !strcmp(card_name,"CURE_GEL_POINT_CARD") )            
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }  
  
  if ( !strcmp(card_name,"CURE_SPECIES_NUMBER_CARD") )           
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }  
  
  if ( !strcmp(card_name,"DATA_SENS_CARD") )            	  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }  
  
  if ( !strcmp(card_name,"DEBUG_CARD") )           
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }   
  
  if ( !strcmp(card_name,"DELTA_S_CARD") )             
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; } 
 
  if ( !strcmp(card_name,"DELTA_T_CARD") )              
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; } 
 
  if ( !strcmp(card_name,"DENSITY_OF_LIQUID_PHASE_IN_POROUS_MEDIA_CARD") )           		  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; } 
 
  if ( !strcmp(card_name,"DENSITY_OF_SOLVENTS_IN_GAS_PHASE_IN_POROUS_MEDIA_CARD") )           		  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; } 
 
  if ( !strcmp(card_name,"DIFFUSE_MASS_FLUX_VECTORS_CARD") )            		  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; } 
 
  if ( !strcmp(card_name,"DISCONTINUOUS_VELO_BC_CARD") )                  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; } 
 
  if ( !strcmp(card_name,"DOMAIN_MAPPING_FILE_CARD") )           
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; } 
 
  if ( !strcmp(card_name,"EIGENSOLVER_SPECIFICATIONS_DELIMITER_CARD") )            
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; } 
 
  if ( !strcmp(card_name,"EIGEN_INITIAL_SHIFTS_CARD") )           
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; } 
 
  if ( !strcmp(card_name,"EIGEN_INITIAL_VECTOR_WEIGHT_CARD") )            
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; } 
 
  if ( !strcmp(card_name,"EIGEN_MAXIMUM_ITERATIONS_CARD") )           
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; } 
 
  if ( !strcmp(card_name,"EIGEN_NUMBER_OF_FILTER_STEPS_CARD") )            
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; } 
 
  if ( !strcmp(card_name,"EIGEN_NUMBER_OF_MODES_CARD") )           
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; } 
 
  if ( !strcmp(card_name,"EIGEN_RECORD_MODES_CARD") )            
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; } 
  
  if ( !strcmp(card_name,"EIGEN_RECYCLE_CARD") )
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
   
  if ( !strcmp(card_name,"EIGEN_SIZE_OF_KRYLOV_SUBSPACE_CARD") )
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; } 
  
  if ( !strcmp(card_name,"EIGEN_TOLERANCE_CARD") )
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; } 
  
  if ( !strcmp(card_name,"EIGEN_WAVE_NUMBER_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; } 

  if ( !strcmp(card_name,"ELECTRICAL_CONDUCTIVITY_CARD") )  
  {  sprintf(msg, "Defaulting %s to CONSTANT 1. ",card_name); 
     declare_default_warning(msg); 
     mat_ptr->Elec_ConductivityModel = CONSTANT;	
     mat_ptr->electrical_conductivity = 1.;
     return; }    
          
  if ( !strcmp(card_name,"ELEMENT_MAPPING_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"END_OF_AUGMENTING_CONDITIONS_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"END_OF_BC_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"END_OF_CC_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }  
  
  if ( !strcmp(card_name,"END_OF_DATA_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }		  
  
  if ( !strcmp(card_name,"END_OF_DATA_SENS_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }		  
  
  if ( !strcmp(card_name,"END_OF_EQUATIONS_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"END_OF_FLUX_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }		  
  
  if ( !strcmp(card_name,"END_OF_FLUX_SENSITIVITY_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }		  
  
  if ( !strcmp(card_name,"END_OF_HC_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }  		
  
  if ( !strcmp(card_name,"END_OF_MAT_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"END_OF_PARTICLES_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }		  
  
  if ( !strcmp(card_name,"END_OF_TC_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }  
  
  if ( !strcmp(card_name,"END_OF_VOLUME_INT_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }		  
  
  if ( !strcmp(card_name,"END_TABLE_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"ENERGY_CONDUCTION_VECTORS_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }		  
  
  if ( !strcmp(card_name,"ENERGY_FLUXLINES_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }		  
  
  if ( !strcmp(card_name,"ERROR_ZZ_ELEM_SIZE_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }		  
  
  if ( !strcmp(card_name,"ERROR_ZZ_HEAT_FLUX_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }		  
  
  if ( !strcmp(card_name,"ERROR_ZZ_PRESSURE_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }		  
  
  if ( !strcmp(card_name,"ERROR_ZZ_VELOCITY_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }		  
  
  if ( !strcmp(card_name,"EXTERNAL_FIELD_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"FEM_FILE_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"FEM_FILE_SPECIFICATIONS_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"FILE_EQUALS_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"FILE_EQUALS_WITH_NAME_EQUALS_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"FILL_CONTOURS_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }		  
  
  if ( !strcmp(card_name,"FILL_SUBCYCLE_CARD") )  
  {  tran->exp_subcycle = 10;
     sprintf(msg, "Defaulting %s, tran->exp_subcycle = 10 ",card_name); 
     declare_default_warning(msg); return; }  
  
  if ( !strcmp(card_name,"FILTER_CONCENTRATION_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"FINAL_PARAMETER_VALUE_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }  
  
  if ( !strcmp(card_name,"FIRST_INVARIANT_OF_STRAIN_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }		  
  
  if ( !strcmp(card_name,"FIX_BC_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"FLUX_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }		  
  
  if ( !strcmp(card_name,"FLUX_SENSITIVITY_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }		  
  
  if ( !strcmp(card_name,"GAS_PHASE_DARCY_VELOCITY_IN_POROUS_MEDIA_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }		  
  
  if ( !strcmp(card_name,"BC_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"GENERAL_SPECIFICATIONS_DELIMITER_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }  
  
  if ( !strcmp(card_name,"GRID_PECLET_NUMBER_IN_POROUS_MEDIA_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }		  
  
  if ( !strcmp(card_name,"GUESS_FILE_CARD") )  
  {  sprintf(msg, " Defaulting %s to ARBITRARY",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"HC_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; } 
  
  if ( !strcmp(card_name,"HIGH_RATE_VISCOSITY_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"HUNTING_SPECIFICATIONS_DELIMITER_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }  
  
  if ( !strcmp(card_name,"INERTIA_COEFFICIENT_CARD") )  
  {  
    if( mat_ptr->PorousMediaType == POROUS_BRINKMAN )
    {
      sprintf(msg, "Meida type is not POROUS_BRINKMAN.  Ignoring Inertia Coefficient card and defaulting inertia coefficient to 0.0.");
      declare_default_warning(msg);  	  
      mat_ptr->Inertia_coefficient = 0.0;
    }
    return; 
  }  
  
  if ( !strcmp(card_name,"INITIALIZE_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }  
  
  if ( !strcmp(card_name,"INITIAL_GUESS_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"INITIAL_GUESS_OF_TP_PARAMETER_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }  
  
  if ( !strcmp(card_name,"INITIAL_PARAMETER_VALUE_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }  
  
  if ( !strcmp(card_name,"INITIAL_TIME_CARD") )  
  {  tran->init_time = 0.0;
     sprintf(msg, "Defaulting %s, tran->init_time = 0.0",card_name); 
     declare_default_warning(msg); return; } 
  
  if ( !strcmp(card_name,"JACOBIAN_REFORM_TIME_STRIDE_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }  
  
  if ( !strcmp(card_name,"LAGRANGIAN_CONVECTION_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }		  
  
  if ( !strcmp(card_name,"LAME_LAMBDA_CARD") )  
  {  sprintf(msg, "Defaulting %s, lame_lambda=1.0, lame_lambda_model=CONSTANT",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"LAME_MU_CARD") )  
  {  sprintf(msg, "Defaulting %s, lame_mu=1.0, lame_mu_model=CONSTANT",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"LEVEL_SET_CONTROL_WIDTH_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }  
  
  if ( !strcmp(card_name,"LEVEL_SET_INITIALIZATION_METHOD_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }  
  
  if ( !strcmp(card_name,"LEVEL_SET_INITIALIZATION_SURFACE_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }  
  
  if ( !strcmp(card_name,"LEVEL_SET_INTERFACE_TRACKING_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }                                		
  
  if ( !strcmp(card_name,"LEVEL_SET_RENORMALIZATION_FREQUENCY_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; } 
  
  if ( !strcmp(card_name,"LEVEL_SET_RENORMALIZATION_METHOD_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; } 
  
  if ( !strcmp(card_name,"LEVEL_SET_RENORMALIZATION_TOLERANCE_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }  
  
  if ( !strcmp(card_name,"LINEAR_STABILITY_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; } 
  
  if ( !strcmp(card_name,"LIQUID_CONSTITUTIVE_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"LIQUID_PHASE_DARCY_VELOCITY_IN_POROUS_MEDIA_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }		  
  
  if ( !strcmp(card_name,"LOCA_METHOD_CARD") ) 
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }  
  
  if ( !strcmp(card_name,"LOW_RATE_VISCOSITY_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"MASS_DIFFUSION_VECTORS_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }		  
  
  if ( !strcmp(card_name,"MASS_FLUXLINES_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }		  
  
  if ( !strcmp(card_name,"MATERIAL_FILE_NAME_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"MATERIAL_ID_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }  
  
  if ( !strcmp(card_name,"MATERIAL_PROPERTY_TAG_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }  
  
  if ( !strcmp(card_name,"MATERIAL_PROPERTY_TAG_SUBINDEX_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }  
  
  if ( !strcmp(card_name,"MATRIX_ABSOLUTE_THRESHOLD_CARD") ) 
  {  strcpy(Matrix_Absolute_Threshold,     "0");
     sprintf(msg, "Defaulting %s, Matrix_Absolute_Threshold = 0",card_name); 
     declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"MATRIX_AUXILIARY_VECTOR_CARD") )  
  {  strcpy(Matrix_Auxiliary_Vector,	"resid");
     sprintf(msg, "Defaulting %s, Matrix_Auxiliary_Vector = resid",card_name); 
     declare_default_warning(msg); return; } 
  
  if ( !strcmp(card_name,"MATRIX_BILU_THRESHOLD_CARD") )  
  {  strcpy(Matrix_BILU_Threshold,         "0");
     sprintf(msg, "Defaulting %s, Matrix_BILU_Threshold = 0",card_name); 
     declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"MATRIX_DROP_TOLERANCE_CARD") )  
  {  strcpy(Matrix_Drop_Tolerance,		"0");
     sprintf(msg, "Defaulting %s, Matrix_Drop_Tolerance = 0",card_name); 
     declare_default_warning(msg); return; }  
  
  if ( !strcmp(card_name,"MATRIX_FACTORIZATION_OVERLAP_CARD") )  
  {  strcpy(Matrix_Factor_Overlap,		"none");
     sprintf(msg, "Defaulting %s, Matrix_Factor_Overlap = none",card_name); 
     declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"MATRIX_FACTORIZATION_REUSE_CARD") )  
  {  strcpy(Matrix_Factorization_Reuse,	"recalc");
     sprintf(msg, "Defaulting %s, Matrix_Factorization_Reuse = recalc",card_name); 
     declare_default_warning(msg); return; }  
  
  if ( !strcmp(card_name,"MATRIX_FACTORIZATION_SAVE_CARD") )  		  
  {  strcpy(Matrix_Factorization_Save,	"0");
     sprintf(msg, "Defaulting %s, Matrix_Factorization_Save = 0",card_name); 
     declare_default_warning(msg); return; } 
   
  if ( !strcmp(card_name,"MATRIX_FILL_FACTOR_CARD") )    
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"MATRIX_GRAPH_FILLIN_CARD") )  
  {  strcpy(Matrix_Graph_Fillin,		"0"); /* Aztec 2 */
     sprintf(msg, "Defaulting %s, Matrix_Graph_Fillin = 0 (Aztec 2)",card_name); 
     declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"MATRIX_ILUT_FILL_FACTOR_CARD") )  		  
  {  strcpy(Matrix_ILUT_Fill_Factor,	"1");
     sprintf(msg, "Defaulting %s, Matrix_ILUT_Fill_Factor = 1",card_name); 
     declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"MATRIX_OUTPUT_TYPE_CARD") )    
  {  strcpy(Matrix_Output_Type,		"none");
     sprintf(msg, "Defaulting %s, Matrix_Output_Type = none",card_name); 
     declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"MATRIX_OVERLAP_TYPE_CARD") )  
  {  strcpy(Matrix_Overlap_Type,		"standard");
     sprintf(msg, "Defaulting %s, Matrix_Overlap_Type = standard",card_name); 
     declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"MATRIX_POLYNOMIAL_ORDER_CARD") )   
  {  strcpy(Matrix_Polynomial_Order,	"3");
     sprintf(msg, "Defaulting %s, Matrix_Polynomial_Order = 3",card_name); 
     declare_default_warning(msg); return; }  
  
  if ( !strcmp(card_name,"MATRIX_RELATIVE_THRESHOLD_CARD") )  
  {  strcpy(Matrix_Relative_Threshold,     "0");
     sprintf(msg, "Defaulting %s, Matrix_Relative_Threshold = 0",card_name); 
     declare_default_warning(msg); return; }  
  
  if ( !strcmp(card_name,"MATRIX_REORDER_CARD") )  
  {  strcpy(Matrix_Reorder,                "none");
     sprintf(msg, "Defaulting %s, Matrix_Reorder = none",card_name); 
     declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"MATRIX_RESIDUAL_NORM_TYPE_CARD") )    
  {  strcpy(Matrix_Residual_Norm_Type,	"r0");
     sprintf(msg, "Defaulting %s, Matrix_Residual_Norm_Type = r0",card_name); 
     declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"MATRIX_RILU_RELAX_FACTOR_CARD") )  		  
  {  strcpy(Matrix_RILU_Relax_Factor,	"1");
     sprintf(msg, "Defaulting %s, Matrix_RILU_Relax_Factor = 1",card_name); 
     declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"MATRIX_SCALING_CARD") )    
  {  strcpy(Matrix_Scaling,		"none");
     sprintf(msg, "Defaulting %s, Matrix_Scaling = none",card_name); 
     declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"MATRIX_STORAGE_FORMAT_CARD") )    
  {  strcpy(Matrix_Format,			"msr");
     sprintf(msg, "Defaulting %s, Matrix_Format = msr",card_name); 
     declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"MATRIX_SUBDOMAIN_SOLVER_CARD") )    
  {  strcpy(Matrix_Subdomain_Solver,	"ilut"); /* Aztec 2 */
     sprintf(msg, "Defaulting %s, Matrix_Subdomain_Solver = ilut (Aztec 2).",card_name); 
     declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"MAXIMUM_LINEAR_ITERATIONS_SOLVE_CARD") )  
  {  strcpy(Matrix_Maximum_Iterations,	"500");
     sprintf(msg, "Defaulting %s, Matrix_Maximum_Iterations = 500.",card_name); 
     declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"MAXIMUM_NUMBER_OF_PATH_STEPS_CARD") )    
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"MAXIMUM_NUMBER_OF_TIME_STEPS_CARD") )    
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"MAXIMUM_PATH_STEP_CARD") )    
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"MAXIMUM_PATH_VALUE_CARD") )    
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"MAXIMUM_TIME_CARD") )    
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"MAXIMUM_TIME_STEP_CARD") )    
  {  Delta_t_max = 1.e12;
     sprintf(msg, "Defaulting %s, Delta_t_max = 1.e12. ",card_name); 
     declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"MEAN_SHEAR_RATE_CARD") )  		  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"MECHANICAL_PROPERTIES_VISCOSITY_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"MECHANICAL_PROPERTIES_VISCOSITY_CARD") )  
  {  
    if( pd_glob[mn]->MeshMotion == TOTAL_ALE)
    {
      sprintf(msg, " Defaulting %s to NONLINEAR",card_name); 
      cr_glob[mn]->MeshFluxModel = NONLINEAR;
      declare_default_warning(msg); return; 
    }
  }  
  
  if ( !strcmp(card_name,"MESH_MOTION_CARD") )  
  {  
    for (i=0; i < number_of_materials_found; i++)
    {
      if( pd_glob[i]->MeshMotion == 0 ) 
      {
        pd_glob[i]->MeshMotion = ARBITRARY;
        cr_glob[i]->MeshMotion = ARBITRARY;
        sprintf(msg, " Defaulting %s for material %d to ARBITRARY.",card_name,i+1); declare_default_warning(msg);
      }
    }
    return; }
  
  if ( !strcmp(card_name,"MESH_STRAIN_TENSOR_CARD") )  		  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"MESH_STRESS_TENSOR_CARD") )  		  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"MINIMUM_PATH_STEP_CARD") )    
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"MINIMUM_TIME_STEP_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"MOBILITY_PARAMETER_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"MODIFIED_NEWTON_TOLERANCE_CARD") )    
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"MOVING_MESH_RESIDUALS_CARD") )  		  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"NAVIER_STOKES_RESIDUALS_CARD") )  		  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"NEWTON_CORRECTION_FACTOR_CARD") )    
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"NORMALIZED_CORRECTION_TOLERANCE_CARD") )    
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"NORMALIZED_RESIDUAL_TOLERANCE_CARD") )    
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"NORMAL_AND_TANGENT_VECTORS_CARD") )  		  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"NUMBER_OF_AUGMENTING_CONDITIONS_CARD") )  	
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"NUMBER_OF_BOUNDARY_CONDITIONS_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"NUMBER_OF_CONTINUATION_CONDITIONS_CARD") )    
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"NUMBER_OF_EQUATIONS_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"NUMBER_OF_HUNTING_CONDITIONS_CARD") )   		
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"NUMBER_OF_JACOBIAN_FILE DUMPS_CARD") )    
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"NUMBER_OF_MATERIALS_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"NUMBER_OF_NEWTONIAN_ITERATIONS_CARD") )    	
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"NUMBER_OF_TP_CONTINUATION_CONDITIONS_CARD") )    
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"ORTHOGONALIZATION_CARD") )    
  {  strcpy(Matrix_Orthogonalization,	"classic");
     sprintf(msg, "Defaulting %s, Matrix_Orthogonalization = classic",card_name); 
     declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"OUTPUT_FILE_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"OUTPUT_LEVEL_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"PARTICLE_CARD") )  		  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"PARTICLE_VELOCITY_DIVERGENCE_CARD") )  		  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"PATH_STEP_ERROR_CARD") )    
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"PATH_STEP_PARAMETER_CARD") )    
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }

  if ( !strcmp(card_name,"PLASTICITY_EQUATION_CARD") )           
  {  sprintf(msg, "Defaulting %s, ConstitutiveEquation = NO_MODEL",card_name); 
     evpl_glob[mn]->ConstitutiveEquation = NO_MODEL;
     declare_default_warning(msg); return; }
       
  if ( !strcmp(card_name,"POLYMER_CONSTITUTIVE_EQUATION_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"POLYMER_STRESS_FORMULATION_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"POLYMER_TIME_CONSTANT_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"POLYMER_VISCOSITY_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"POLYMER_WEIGHTING_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"POLYMER_WEIGHT_FUNCTION_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"POROUS_SATURATION_CARD") )  		  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"POST_PROCESSING_DATA_CARD") )  		  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"POST_PROCESSING_DATA_DELIMITER_CARD") )  		  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"POST_PROCESSING_DATA_SENSITIVITIES_DELIMITER_CARD") )  		  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"POST_PROCESSING_FLUXES_SECTION_DELIMITER_CARD") )  		  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"POST_PROCESSING_FLUX_SENSITIVITIES_SECTION_DELIMITER_CARD") )  		  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"POST_PROCESSING_PARTICLE_TRACKING_CALCULATIONS_DELIMITER_CARD") )  		  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"POST_PROCESSING_SPECIFICATION_SECTION_DELIMITER_CARD") )  		  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"POST_PROCESSING_VOLUMETRIC_INTEGRATION_DELIMITER_CARD") )  		  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"POWER_LAW_EXPONENT_EQUATION_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"PRECONDITIONER_CARD") )    
  {  strcpy(Matrix_Preconditioner,		"none");
     sprintf(msg, "Defaulting %s, Matrix_Preconditioner = none ",card_name);
     declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"PRESSURE_CONTOURS_CARD") )  		  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"PRESSURE_DATUM_CARD") )  
  {  upd->Pressure_Datum =  1.0132500000E6;
     sprintf(msg, "Defaulting %s, upd->Pressure_Datum =  1.0132500000E6.",card_name); 
     declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"PRESSURE_STABILIZATION_CARD") )   
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"PRESSURE_STABILIZATION_SCALING_CARD") )    
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"PRINTING_FREQUENCY_CARD") )    
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"PROBLEM_DESC_DELIMITER_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"REAL_SOLID_STRESS_TENSOR_CARD") )  		  	  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"RESIDUAL_RATIO_TOLERANCE_CARD") )  
  {  strcpy(Matrix_Convergence_Tolerance,	"1e-6");
     sprintf(msg, "Defaulting %s, Matrix_Convergence_Tolerance = 1e-6",card_name); 
     declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"SECOND_FREQUENCY_CARD") )    		  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"SECOND_FREQUENCY_TIME_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"SECOND_INVARIANT_OF_STRAIN_CARD") )  		  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"SIZE_OF_KRYLOV_SUBSPACE_CARD") )  
  {  strcpy(Matrix_Krylov_Subspace,	"30");
     sprintf(msg, "Defaulting %s, Matrix_Krylov_Subspace = 30",card_name); 
     declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"SOLID_CONSTITUTIVE_EQUATION_CARD") )  
  {  cr_glob[mn]->MeshFluxModel = NONLINEAR;
     if( pd_glob[mn]->MeshMotion == TOTAL_ALE)
     {  
       cr_glob[mn]->RealSolidFluxModel = NONLINEAR;
     }     
     sprintf(msg, "Defaulting %s, MeshFluxModel = NONLINEAR",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"SOLID_THERMAL_EXPANSION_CARD") )  
  { sprintf(msg, "Defaulting %s, thermal_expansion=0.0, thermal_expansion_model=CONSTANT",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"reference_temperature_card") )  
  { sprintf(msg, "Defaulting %s, solid_reference_temp=25.0, solid_reference_temp_model=CONSTANT",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"SOLN_FILE_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"SOLUTION_ALGORITHM_CARD") )    
  {  Linear_Solver = SPARSE13a;
     sprintf(msg, "Defaulting %s, Linear_Solver = lu",card_name); 
     declare_default_warning(msg); return; }
     
  if ( !strcmp(card_name,"SOLUTION_TEMPERATURE_CARD") )    
  {  
  }
     
  if ( !strcmp(card_name,"SOLVER_PRESSURE_DATUM_CARD") )    
  {  sprintf(msg, "Defaulting %s, pressure_datum_value = -1.0e12 & pressure_datum_element = -1.",card_name);
     pressure_datum_value = -1.0e12;  
     pressure_datum_element = -1;      
     declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"SOLVER_SPECIFICATIONS_DELIMITER_CARD") )    
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"STEP_CONTROL_AGGRESSIVENESS_CARD") )    
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"STREAMWISE_NORMAL_STRESS_CARD") )  		  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"STREAM_FUNCTION_CARD") )  		  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"STRESS_CONTOURS_CARD") )  		  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"STRESS_FREE_SOLVENT_VOL_FRAC_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"SUPG_VELOCITY_IN_POROUS_MEDIA_CARD") )  		  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"SURFACE_TENSION_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"SUSPENSION_SPECIES_NUMBER_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"TABLE_DATA_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"TC_CARD") )    
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"THERMAL_EXPONENT_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"THIRD_INVARIANT_OF_STRAIN_CARD") )  		  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"TIME_CONSTANT_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"TIME_DERIVATIVES_CARD") )  		  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"TIME_INTEGRATION_CARD") )    
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"TIME_INTEGRATION_SPECIFICATIONS_DELIMITER_CARD") )    
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"TIME_STEP_ERROR_CARD") )    
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"TIME_STEP_PARAMETER_CARD") )   
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"TOTAL_DENSITY_OF_SOLVENTS_IN_POROUS_MEDIA_CARD") )  		  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"TOTAL_VELOCITY_DIVERGENCE_CARD") )  		  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"TP_BC_DATA_FLOAT_TAG_CARD") )    
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"TP_BOUNDARY_CONDITION_ID_CARD") )    
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"TP_CONTINUATION_TYPE_CARD") )    
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"TP_MATERIAL_PROPERTY_TAG_SUBINDEX_CARD") )    
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"TP_PARAMETER_FINAL_VALUE_CARD") )    
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"TP_PARAMETER_MATERIAL_ID_CARD") )    
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"TP_PARAMETER_MATERIAL_PROPERTY_TAG_CARD") )    
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"UMF_IDIM_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"UMF_XDIM_CARD") )   
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"USER_DEFINED_POST_PROCESSING_CARD") )  	  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"VELOCITY_DIVERGENCE_CARD") )  		  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"VISCOPLASTIC_DEF_GRAD_TENSOR_CARD") )  		  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"VISCOSITY__CARD") )  		  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"VOLUME_INT_CARD") )  		  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"VORTICITY_VECTOR_CARD") )  		  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"WRITE_INITIAL_SOLUTION_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"WRITE_INTERMEDIATE_RESULTS_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"YEILD_EXPONENT_CARD") )  
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"YEILD_STRESS_CARD") )           
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); return; }
  
  if ( !strcmp(card_name,"EXAMPLE_CARD") )           
  {  sprintf(msg, "Defaulting %s, ",card_name); declare_default_warning(msg); 
     /* ADD THE DEFAULT PROCESSING FOR THE EXAMPLE CARD HERE. */
     return; 
  }  

  
sprintf(msg, "Trying to defualt unknown card %s",card_name); declare_default_warning(msg);

}

void consistency
  (
  int file_index
  )
{
  struct Elastic_Constitutive  *dum_ptr;
  
  if (file_index == 0)
  {
    /* Check consistency on Input file */
  }
  else
  {
    /* check consistency on Material file */
    
    /* ADD ANY COCNISTENCY CHECKS ASSOICATED WITH THE EXAMPLE CARD HERE.*/	
    
    if( pd_glob[mn]->MeshMotion == TOTAL_ALE)
    {
    
      /* First transfer already parsed real-properties to the appropriate
       * structure elements. Do this by swapping pointers. You need to maintain
       * though some of the elastic (elc) and constitutive relations (cr)
       * constants in both. 
       */

       cr_glob[mn]->RealSolidFluxModel = cr_glob[mn]->MeshFluxModel;
       dum_ptr = elc_rs_glob[mn];
       elc_rs_glob[mn] = elc_glob[mn];
       elc_glob[mn] = dum_ptr;
       elc_glob[mn]->Strss_fr_sol_vol_frac = elc_rs_glob[mn]->Strss_fr_sol_vol_frac;
       elc_glob[mn]->thermal_expansion = 0.;
       elc_glob[mn]->thermal_expansion_model = CONSTANT;
       elc_glob[mn]->solid_reference_temp_model = CONSTANT;	
       elc_glob[mn]->solid_reference_temp     = 25.; 
    }
    if(pd_glob[mn]->MassFluxModel == FICKIAN_CHARGED)
    {
       sprintf(msg, "The FICKIAN_CHARGED flux model is being used but no Solution Temperature card was specified."); 
       declare_error(msg);
    }
  }  
}

void initialize_input_file_values
  (
  )
{
/* this routine sets initial values associated with the input file.
   these values may be overwirtten later as input file cards are read in. 
*/
  int i;
  
  Num_Var_Init = 0;
  Num_Var_External = 0;
  efv->ev = T_NOTHING;
  ls = NULL;
  cont->sensvec_id = -1;   /* look here jjjinit */
  nCC = 0;
  nTC = 0;
  cont->sensvec_id = -1;
  ContType = cont->upType = -1;
  cont->upBCID = -1;
  cont->BegParameterValue = 0.0;
  cont->EndParameterValue = 1.0;
  Delta_s_max = 1.e+12; 
  loca_in->StepAggr = 0.0;
  loca_in->DpDs2 = 0.5;
  loca_in->DpDsHi = 1.0;
  loca_in->Texp = 0.0;
  loca_in->MaxTS = 0.0;
  loca_in->TPupBCID = -1;
  loca_in->TPGuess = 0.0;
  loca_in->TPFinal = 1.0;
  ERROR_ZZ_VEL_ELSIZE = ERROR_ZZ_Q_ELSIZE = ERROR_ZZ_P_ELSIZE = -1;
  for ( i=0; i<MAX_NUMBER_MATLS; i++)
    {
      pd_glob[i]->Num_Porous_Eqn  = 0;
    }
  
}

void initialize_material_file_values
  (
  int mn
  )
{
  int v;
  model_read = 1;
  mat_ptr = mp_glob[mn];
  if (pd_glob[mn]->MeshMotion == LAGRANGIAN || pd_glob[mn]->MeshMotion == TOTAL_ALE)
  {
    pd_glob[mn]->MeshInertia = 0; /* this pricessing from Convective Lagrangian Velocity card */
  }
  elc_glob[mn]->lame_mu = 1.;
  elc_glob[mn]->lame_lambda = 1.;
  elc_glob[mn]->lame_mu_model = CONSTANT;
  elc_glob[mn]->lame_lambda_model = CONSTANT;
  for ( v=0; v<MAX_CONC + MAX_VARIABLE_TYPES; v++)
    {
      elc_glob[mn]->d_lame_mu[v] = 0.;
      elc_glob[mn]->d_lame_lambda[v] = 0.;
    }
  elc_glob[mn]->thermal_expansion = 0.;
  elc_glob[mn]->solid_reference_temp = 25.;
  elc_glob[mn]->thermal_expansion_model = CONSTANT;
  elc_glob[mn]->solid_reference_temp = CONSTANT;

  elc_rs_glob[mn]->lame_mu = 1.;
  elc_rs_glob[mn]->lame_lambda = 1.;
  elc_rs_glob[mn]->lame_mu_model = CONSTANT;
  elc_rs_glob[mn]->lame_lambda_model = CONSTANT;
  for ( v=0; v<MAX_CONC + MAX_VARIABLE_TYPES; v++)
    {
      elc_rs_glob[mn]->d_lame_mu[v] = 0.;
      elc_rs_glob[mn]->d_lame_lambda[v] = 0.;
    }
  elc_rs_glob[mn]->thermal_expansion = 0.;
  elc_rs_glob[mn]->solid_reference_temp = 25.;
  elc_rs_glob[mn]->thermal_expansion_model = CONSTANT;
  elc_rs_glob[mn]->solid_reference_temp = CONSTANT;
  Num_Var_Init_Mat[mn]=0;
  number_of_species_source_cards_needed[mn] = mp_glob[mn]->Num_Species = pd_glob[mn]->Num_Species; 
  for (i=0; i<Num_BC; i++) 
  {
    if (BC_Types[i].BC_Name == VL_POLY_BC ||
	BC_Types[i].BC_Name == YFLUX_EQUIL_BC ||
	BC_Types[i].BC_Name == VL_EQUIL_PRXN_BC) 
    {
      read_bc_mp = i;
    }
  }
  pd_glob[mn]->MeshInertia = 0;
  elc_glob[mn]->v_mesh_sfs_model = 0;
}  /* end initialize_material_file_values */


int file_to_string 
  (
  FILE *file, char string[], 
  int max_number_chars_in
  )
  /******************************************************************
   *
   *  This routine reads a maximum of max_number_chars_in characters
   *  from a file and places them into a string.   Newline characters in
   *  the file are replaced with spaces in the string.
   *  If an error occurs, -1 is returned and an error message
   *  is written to standard error.  Upon successful
   *  completion, the number of characters read plus 1 for the
   *  null character at the end of the string is returned.
   *  After the routine, the file pointer is positioned after
   *  the terminating character.
   *
   *    Author:			Jim Simmons 9114
   *    Date:			7/3/00
   *    revised:		
   *
   *     Parameter list:
   *
   *    file   == pointer to file from which the characters are to be read.
   *              is positioned after the last character read.
   *    string == On output 'string' contains the characters read
   *  	          from the input stream, terminated by a null
   *              character.
   *    max_number_chars_in == The maximum numbers of character so be processed.
   *
   *    return == -1 if 1) EOF is encountered
   *                    2) MAX_CHAR_IN_INPUT characters are read before
   *                       an end of line or ch character is read.
   *           == strlen(string) + 1 if not an error condition
   *******************************************************************/
{	 
  int i = 0;
  int new_ch;

#ifdef linux 
{ /* do nothing */ }
#else
{
  yyclearin;
  yyrestart();
}
#endif

  bzero(mystring, 10000*sizeof(char));

  while ( (i < max_number_chars_in) && (new_ch = getc(file))
          && (new_ch != EOF) ) 
    {
    /*if (new_ch != '\n')*/ 
  	string[i++] = new_ch;
    /*else
  	string[i++] = ' ';*/
    }
  if (i == max_number_chars_in - 1)
    { 
    fprintf(stderr, "file_to_string:%d The maximum number of characters to be read is exceeded.\n", 
    max_number_chars_in);
    return(-1);
    } 
  string[i] = '\0';
  if(ProcID == 0) { /* Code to mpi broadcast the string goes here.  To be added later by geniuses; will work first time. */ }
  return (i+1);
}
        
int yyerror(const char *msg)
{ 
  yyclearin; 
  if(error_found_in_last_line == 0)
    {
      if (ProcID == 0) fprintf(parser_log, " <---<<< ERROR %d >>> ", error_number);
      error_found_in_last_line = 1;
      if(yytext[0] == '\n') 
      {
        if (ProcID == 0) 
        {
          fprintf(parser_log,"\n <<< ERROR %d Line %d: near \"%s\" >>>\n", error_number, line_number, debug_clue);
          fprintf(error_log,"  <<< ERROR %d Line %d: near \"%s\" >>>\n", error_number, line_number, debug_clue);
          fprintf(stdout,"\n <<< ERROR %d Line %d: near \"%s\" >>>", error_number, line_number, debug_clue);               
        }
        error_found_in_last_line = 0;
        line_number++;
        error_number++;
        if (ProcID == 0) fprintf(parser_log,"%d ",line_number);
      }
    } 
return 0; 
}

int declare_error
  (
  char *msg
  )
{
  if (ProcID == 0) 
  {
    fprintf(stdout, " <<< ERROR %d Line %d: %s >>>\n", error_number, line_number, msg);
    fprintf(error_log, " <<< ERROR %d Line %d: %s >>>\n", error_number, line_number, msg);
    fprintf(parser_log, "\n <<< ERROR %d Line %d: %s >>>", error_number, line_number, msg);
    /*fprintf(parser_log,"%d ",line_number );*/

  }
  error_number++;
  return 0; 
}

int declare_error_without_cr
  (
  char *msg
  )
{
  if (ProcID == 0) 
  {
    fprintf(parser_log, "<<< ERROR %d Line %d: %s >>>", (error_number), line_number, msg); 
    fprintf(stdout, " <<< ERROR %d: %s >>>\n", (error_number),  msg);
    fprintf(error_log, " <<< ERROR %d Line %d: %s >>>\n", (error_number), line_number, msg); 
  } 
  error_number++;
  line_number++;
  if(ProcID == 0)  fprintf(parser_log,"\n%d ",line_number );
  return 0; 
}

int declare_warning
  (
  char *msg
  )
{
  if (ProcID == 0)
  {
    fprintf(parser_log, "\n  <<< Warning %d Line %d: %s >>>", warning_number, line_number, msg);
    fprintf(stdout, " <<< Warning %d: %s >>>\n", warning_number, msg);  
    fprintf(error_log, " <<< Warning %d Line %d: %s >>>\n", warning_number, line_number, msg);  
  }
  warning_number++;
  /*line_number++;*/
  /*if (ProcID == 0) fprintf(parser_log,"%d ", line_number);*/
  return 0; 
}

int declare_default_warning
  (
  char *msg
  )
{
  if (ProcID == 0)
  {
    if (first_default)
    {
      fprintf(stdout, " <<< Warning %d: Some cards are being defualted.  See error.log or parser.log for details. >>>\n", warning_number);
      fprintf(parser_log, " <<< Warning %d: Some cards are being defualted.  See warnings below for details. >>>\n", warning_number);
      fprintf(error_log, " <<< Warning %d: Some cards are being defualted.  See warnings below for details. >>>\n", warning_number);
      first_default = FALSE; warning_number++;
    }
    fprintf(parser_log, "\n <<< Warning %d: %s >>>", warning_number, msg);
    fprintf(error_log, " <<< Warning %d: %s >>>\n", warning_number, msg); 
    fprintf(stdout, " <<< Warning %d: %s >>>\n", warning_number, msg);  
  }
  warning_number++;
  return 0; 
}

int setup_GD_table
  (
  )
{
  BC_Types[iBC].table = ( struct Data_Table * ) smalloc( sizeof( struct Data_Table ) );
  BC_Tables[num_BC_Tables] = BC_Types[iBC].table;
  BC_Types[iBC].table->interp_method = LINEAR;
  BC_Types[iBC].table->columns  = 2;  
  BC_Types[iBC].table->t_index[0] = BC_Types[iBC].BC_Data_Int[2];
  BC_Types[iBC].table->f_index = BC_Types[iBC].BC_Data_Int[0];    
   /* Find the names of the ordinate and the abscissa */
  for( k = 0; BC_Types[iBC].BC_Data_Int[0] != EQ_Name[k].Index && k < Num_EQ_Names; k++);
  {
    BC_Types[iBC].table->f_name = EQ_Name[k].name1;
  }
  if ( BC_Types[iBC].BC_Data_Int[2] == GD_TIME_TABLE )
  {
    strcpy( BC_Types[iBC].table->t_name[0], "time");
  }
  else
  {
    for( k = 0; BC_Types[iBC].BC_Data_Int[2] != Var_Name[k].Index && k < Num_Var_Names; k++);
    {
      strcpy( BC_Types[iBC].table->t_name[0], Var_Name[k].name1 );
    }
  } 
  BC_Types[iBC].table_index = num_BC_Tables++;  
  return(1);
}

int setup_bc_table 
  ( 
  int x_axis,
  int y_axis,
  int interpolation
  )
{
  BC_Types[iBC].table = ( struct Data_Table * ) smalloc( sizeof( struct Data_Table ) );
  BC_Types[iBC].table->columns  = columns_expected;
  switch(x_axis)
  {  
    case TIME_ORD:
      strcpy( BC_Types[iBC].table->t_name[0],"TIME");
      BC_Types[iBC].table->t_index[0] = -1;
      BC_Types[iBC].desc->sens[MESH_DISPLACEMENT1] = 1;      
      break;
    case X_ORD:
      strcpy( BC_Types[iBC].table->t_name[0],"X");
      BC_Types[iBC].table->t_index[0] = 0;
      BC_Types[iBC].desc->sens[MESH_DISPLACEMENT1] = 1;
      break;
    case Y_ORD:
      strcpy( BC_Types[iBC].table->t_name[0],"Y");
      BC_Types[iBC].table->t_index[0] = 1;
      BC_Types[iBC].desc->sens[MESH_DISPLACEMENT2] = 1;
      break;
    case Z_ORD:
      strcpy( BC_Types[iBC].table->t_name[0],"Z");
      BC_Types[iBC].table->t_index[0] = 2;
      BC_Types[iBC].desc->sens[MESH_DISPLACEMENT3] = 1;
      break;
    case XY_ORD:
      strcpy( BC_Types[iBC].table->t_name[0],"X");
      BC_Types[iBC].table->t_index[0] = 0;
      BC_Types[iBC].desc->sens[MESH_DISPLACEMENT1] = 1;
      strcpy( BC_Types[iBC].table->t_name[1],"Y");
      BC_Types[iBC].table->t_index[1] = 1;
      BC_Types[iBC].desc->sens[MESH_DISPLACEMENT2] = 1;
      break;
    case XZ_ORD:
      strcpy( BC_Types[iBC].table->t_name[0],"X");
      BC_Types[iBC].table->t_index[0] = 0;
      BC_Types[iBC].desc->sens[MESH_DISPLACEMENT1] = 1;
      strcpy( BC_Types[iBC].table->t_name[1],"Z");
      BC_Types[iBC].table->t_index[1] = 2;
      BC_Types[iBC].desc->sens[MESH_DISPLACEMENT3] = 1;
      break;
    case YZ_ORD:
      strcpy( BC_Types[iBC].table->t_name[0],"Y");
      BC_Types[iBC].table->t_index[0] = 1;
      BC_Types[iBC].desc->sens[MESH_DISPLACEMENT2] = 1;
      strcpy( BC_Types[iBC].table->t_name[1],"Z");
      BC_Types[iBC].table->t_index[1] = 2;
      BC_Types[iBC].desc->sens[MESH_DISPLACEMENT3] = 1;
      break;
    case YX_ORD:
      strcpy( BC_Types[iBC].table->t_name[0],"Y");
      BC_Types[iBC].table->t_index[0] = 1;
      BC_Types[iBC].desc->sens[MESH_DISPLACEMENT2] = 1;
      strcpy( BC_Types[iBC].table->t_name[1],"X");
      BC_Types[iBC].table->t_index[1] = 0;
      BC_Types[iBC].desc->sens[MESH_DISPLACEMENT1] = 1;
      break;
    case ZX_ORD:
      strcpy( BC_Types[iBC].table->t_name[0],"Z");
      BC_Types[iBC].table->t_index[0] = 2;
      BC_Types[iBC].desc->sens[MESH_DISPLACEMENT3] = 1;
      strcpy( BC_Types[iBC].table->t_name[1],"X");
      BC_Types[iBC].table->t_index[1] = 0;
      BC_Types[iBC].desc->sens[MESH_DISPLACEMENT1] = 1;
      break;
    case ZY_ORD:
      strcpy( BC_Types[iBC].table->t_name[0],"Z");
      BC_Types[iBC].table->t_index[0] = 2;
      BC_Types[iBC].desc->sens[MESH_DISPLACEMENT3] = 1;
      strcpy( BC_Types[iBC].table->t_name[1],"Y");
      BC_Types[iBC].table->t_index[1] = 1;
      BC_Types[iBC].desc->sens[MESH_DISPLACEMENT2] = 1;
      break;
  }

  switch(y_axis)
  {
    case VELOCITY1_ABS:
      BC_Types[iBC].table->f_name = "VELOCITY1";
      BC_Types[iBC].table->f_index = VELOCITY1;
      BC_Types[iBC].desc->equation = R_MOMENTUM1;
      BC_Types[iBC].desc->sens[VELOCITY1] = 1;
      break;
    case VELOCITY2_ABS:
      BC_Types[iBC].table->f_name = "VELOCITY2";
      BC_Types[iBC].table->f_index = VELOCITY2;
      BC_Types[iBC].desc->equation = R_MOMENTUM2;
      BC_Types[iBC].desc->sens[VELOCITY2] = 1;
      break;
    case VELOCITY3_ABS:
      BC_Types[iBC].table->f_name = "VELOCITY3";
      BC_Types[iBC].table->f_index = VELOCITY3;
      BC_Types[iBC].desc->equation = R_MOMENTUM3;
      BC_Types[iBC].desc->sens[VELOCITY3] = 1;
      break;
    case TEMPERATURE_ABS:
      BC_Types[iBC].table->f_name = "TEMPERATURE";
      BC_Types[iBC].table->f_index = TEMPERATURE;
      BC_Types[iBC].desc->equation = R_ENERGY;
      BC_Types[iBC].desc->sens[TEMPERATURE] = 1;
      break;
    case MASS_FRACTION_ABS:
      BC_Types[iBC].table->f_name = "MASS_FRACTION";
      BC_Types[iBC].table->f_index = MASS_FRACTION;
      BC_Types[iBC].desc->equation = R_MASS;
      BC_Types[iBC].desc->sens[MASS_FRACTION] = 1;
      break;
    case SPECIES_ABS:
      BC_Types[iBC].table->f_name = "MASS_FRACTION";
      BC_Types[iBC].table->f_index = MASS_FRACTION;
      BC_Types[iBC].desc->equation = R_MASS;
      BC_Types[iBC].desc->sens[MASS_FRACTION] = 1;
      /*if ( fscanf( ifp, "%d", &BC_Types[iBC].species_eq) != 1)
	{
	  sprintf (msg, "%s:\tError reading species number on TABLE BC \n", yo);
	  EH(-1,err_msg);
	}*/
      break;
    case MESH_DISPLACEMENT1_ABS:
      BC_Types[iBC].table->f_name = "MESH_DISPLACEMENT1";
      BC_Types[iBC].table->f_index = MESH_DISPLACEMENT1;
      BC_Types[iBC].desc->equation = R_MESH1;
      BC_Types[iBC].desc->sens[MESH_DISPLACEMENT1] = 1;
      break;
    case MESH_DISPLACEMENT2_ABS:
      BC_Types[iBC].table->f_name = "MESH_DISPLACEMENT2";
      BC_Types[iBC].table->f_index = MESH_DISPLACEMENT2;
      BC_Types[iBC].desc->equation = R_MESH2;
      BC_Types[iBC].desc->sens[MESH_DISPLACEMENT2] = 1;
      break;
    case MESH_DISPLACEMENT3_ABS:
      BC_Types[iBC].table->f_name = "MESH_DISPLACEMENT3";
      BC_Types[iBC].table->f_index = MESH_DISPLACEMENT3;
      BC_Types[iBC].desc->equation = R_MESH3;
      BC_Types[iBC].desc->sens[MESH_DISPLACEMENT3] = 1;
      break;
    case PRESSURE_ABS:
      BC_Types[iBC].table->f_name = "PRESSURE";
      BC_Types[iBC].table->f_index = PRESSURE;
      BC_Types[iBC].desc->equation = R_PRESSURE;
      BC_Types[iBC].desc->sens[PRESSURE] = 1;
      break;
    case SHEAR_RATE_ABS:
      BC_Types[iBC].table->f_name = "SHEAR_RATE";
      BC_Types[iBC].table->f_index = SHEAR_RATE;
      BC_Types[iBC].desc->equation = R_SHEAR_RATE;
      BC_Types[iBC].desc->sens[SHEAR_RATE] = 1;
      break;
    case S11_ABS:
      BC_Types[iBC].table->f_name = "S11";
      BC_Types[iBC].table->f_index = POLYMER_STRESS11;
      BC_Types[iBC].desc->equation = R_STRESS11;
      BC_Types[iBC].desc->sens[POLYMER_STRESS11] = 1;
      break;
    case S12_ABS:
      BC_Types[iBC].table->f_name = "S12";
      BC_Types[iBC].table->f_index = POLYMER_STRESS12;
      BC_Types[iBC].desc->equation = R_STRESS12;
      BC_Types[iBC].desc->sens[POLYMER_STRESS12] = 1;
      break;
    case S22_ABS:
      BC_Types[iBC].table->f_name = "S22";
      BC_Types[iBC].table->f_index = POLYMER_STRESS22;
      BC_Types[iBC].desc->equation = R_STRESS22;
      BC_Types[iBC].desc->sens[POLYMER_STRESS22] = 1;
      break;
    case S13_ABS:
      BC_Types[iBC].table->f_name = "S13";
      BC_Types[iBC].table->f_index = POLYMER_STRESS13;
      BC_Types[iBC].desc->equation = R_STRESS13;
      BC_Types[iBC].desc->sens[POLYMER_STRESS13] = 1;
      break;
    case S23_ABS:
      BC_Types[iBC].table->f_name = "S23";
      BC_Types[iBC].table->f_index = POLYMER_STRESS23;
      BC_Types[iBC].desc->equation = R_STRESS23;
      BC_Types[iBC].desc->sens[POLYMER_STRESS23] = 1;
      break;
    case S33_ABS:
      BC_Types[iBC].table->f_name = "S33";
      BC_Types[iBC].table->f_index = POLYMER_STRESS33;
      BC_Types[iBC].desc->equation = R_STRESS33;
      BC_Types[iBC].desc->sens[POLYMER_STRESS33] = 1;
      break;
    case S11_1_ABS:
      BC_Types[iBC].table->f_name = "S11_1";
      BC_Types[iBC].table->f_index = POLYMER_STRESS11_1;
      BC_Types[iBC].desc->equation = R_STRESS11_1;
      BC_Types[iBC].desc->sens[POLYMER_STRESS11_1] = 1;
      break;
    case S12_1_ABS:
      BC_Types[iBC].table->f_name = "S12_1";
      BC_Types[iBC].table->f_index = POLYMER_STRESS12_1;
      BC_Types[iBC].desc->equation = R_STRESS12_1;
      BC_Types[iBC].desc->sens[POLYMER_STRESS12_1] = 1;
      break;
    case S22_1_ABS:
      BC_Types[iBC].table->f_name = "S22_1";
      BC_Types[iBC].table->f_index = POLYMER_STRESS22_1;
      BC_Types[iBC].desc->equation = R_STRESS22_1;
      BC_Types[iBC].desc->sens[POLYMER_STRESS22_1] = 1;
      break;
    case S13_1_ABS:
      BC_Types[iBC].table->f_name = "S13_1";
      BC_Types[iBC].table->f_index = POLYMER_STRESS13_1;
      BC_Types[iBC].desc->equation = R_STRESS13_1;
      BC_Types[iBC].desc->sens[POLYMER_STRESS13_1] = 1;
      break;
    case S23_1_ABS:
      BC_Types[iBC].table->f_name = "S23_1";
      BC_Types[iBC].table->f_index = POLYMER_STRESS23_1;
      BC_Types[iBC].desc->equation = R_STRESS23_1;
      BC_Types[iBC].desc->sens[POLYMER_STRESS23_1] = 1;
      break;
    case S33_1_ABS:
      BC_Types[iBC].table->f_name = "S33_1";
      BC_Types[iBC].table->f_index = POLYMER_STRESS33_1;
      BC_Types[iBC].desc->equation = R_STRESS33_1;
      BC_Types[iBC].desc->sens[POLYMER_STRESS33_1] = 1;
      break;
    case S11_2_ABS:
      BC_Types[iBC].table->f_name = "S11_2";
      BC_Types[iBC].table->f_index = POLYMER_STRESS11_2;
      BC_Types[iBC].desc->equation = R_STRESS11_2;
      BC_Types[iBC].desc->sens[POLYMER_STRESS11_2] = 1;
      break;
   case S12_2_ABS:
      BC_Types[iBC].table->f_name = "S12_2";
      BC_Types[iBC].table->f_index = POLYMER_STRESS12_2;
      BC_Types[iBC].desc->equation = R_STRESS12_2;
      BC_Types[iBC].desc->sens[POLYMER_STRESS12_2] = 1;
      break;
   case S22_2_ABS:
      BC_Types[iBC].table->f_name = "S22_2";
      BC_Types[iBC].table->f_index = POLYMER_STRESS22_2;
      BC_Types[iBC].desc->equation = R_STRESS22_2;
      BC_Types[iBC].desc->sens[POLYMER_STRESS22_2] = 1;
      break;
   case S13_2_ABS:
      BC_Types[iBC].table->f_name = "S13_2";
      BC_Types[iBC].table->f_index = POLYMER_STRESS13_2;
      BC_Types[iBC].desc->equation = R_STRESS13_2;
      BC_Types[iBC].desc->sens[POLYMER_STRESS13_2] = 1;
      break;
   case S23_2_ABS:
      BC_Types[iBC].table->f_name = "S23_2";
      BC_Types[iBC].table->f_index = POLYMER_STRESS23_2;
      BC_Types[iBC].desc->equation = R_STRESS23_2;
      BC_Types[iBC].desc->sens[POLYMER_STRESS23_2] = 1;
      break;
   case S33_2_ABS:
      BC_Types[iBC].table->f_name = "S33_2";
      BC_Types[iBC].table->f_index = POLYMER_STRESS33_2;
      BC_Types[iBC].desc->equation = R_STRESS33_2;
      BC_Types[iBC].desc->sens[POLYMER_STRESS33_2] = 1;
      break;
   case S11_3_ABS:
      BC_Types[iBC].table->f_name = "S11_3";
      BC_Types[iBC].table->f_index = POLYMER_STRESS11_3;
      BC_Types[iBC].desc->equation = R_STRESS11_3;
      BC_Types[iBC].desc->sens[POLYMER_STRESS11_3] = 1;
      break;
   case S12_3_ABS:
      BC_Types[iBC].table->f_name = "S12_3";
      BC_Types[iBC].table->f_index = POLYMER_STRESS12_3;
      BC_Types[iBC].desc->equation = R_STRESS12_3;
      BC_Types[iBC].desc->sens[POLYMER_STRESS12_3] = 1;
      break;
   case S22_3_ABS:
      BC_Types[iBC].table->f_name = "S22_3";
      BC_Types[iBC].table->f_index = POLYMER_STRESS22_3;
      BC_Types[iBC].desc->equation = R_STRESS22_3;
      BC_Types[iBC].desc->sens[POLYMER_STRESS22_3] = 1;
      break;
   case S13_3_ABS:
      BC_Types[iBC].table->f_name = "S13_3";
      BC_Types[iBC].table->f_index = POLYMER_STRESS13_3;
      BC_Types[iBC].desc->equation = R_STRESS13_3;
      BC_Types[iBC].desc->sens[POLYMER_STRESS13_3] = 1;
      break;
   case S23_3_ABS:
      BC_Types[iBC].table->f_name = "S23_3";
      BC_Types[iBC].table->f_index = POLYMER_STRESS23_3;
      BC_Types[iBC].desc->equation = R_STRESS23_3;
      BC_Types[iBC].desc->sens[POLYMER_STRESS23_3] = 1;
      break;
   case S33_3_ABS:
      BC_Types[iBC].table->f_name = "S33_3";
      BC_Types[iBC].table->f_index = POLYMER_STRESS33_3;
      BC_Types[iBC].desc->equation = R_STRESS33_3;
      BC_Types[iBC].desc->sens[POLYMER_STRESS33_3] = 1;
      break;
   case S11_4_ABS:
      BC_Types[iBC].table->f_name = "S11_4";
      BC_Types[iBC].table->f_index = POLYMER_STRESS11_4;
      BC_Types[iBC].desc->equation = R_STRESS11_4;
      BC_Types[iBC].desc->sens[POLYMER_STRESS11_4] = 1;
      break;
   case S12_4_ABS:
      BC_Types[iBC].table->f_name = "S12_4";
      BC_Types[iBC].table->f_index = POLYMER_STRESS12_4;
      BC_Types[iBC].desc->equation = R_STRESS12_4;
      BC_Types[iBC].desc->sens[POLYMER_STRESS12_4] = 1;
      break;
   case S22_4_ABS:
      BC_Types[iBC].table->f_name = "S22_4";
      BC_Types[iBC].table->f_index = POLYMER_STRESS22_4;
      BC_Types[iBC].desc->equation = R_STRESS22_4;
      BC_Types[iBC].desc->sens[POLYMER_STRESS22_4] = 1;
      break;
   case S13_4_ABS:
      BC_Types[iBC].table->f_name = "S13_4";
      BC_Types[iBC].table->f_index = POLYMER_STRESS13_4;
      BC_Types[iBC].desc->equation = R_STRESS13_4;
      BC_Types[iBC].desc->sens[POLYMER_STRESS13_4] = 1;
      break;
   case S23_4_ABS:
      BC_Types[iBC].table->f_name = "S23_4";
      BC_Types[iBC].table->f_index = POLYMER_STRESS23_4;
      BC_Types[iBC].desc->equation = R_STRESS23_4;
      BC_Types[iBC].desc->sens[POLYMER_STRESS23_4] = 1;
      break;
   case S33_4_ABS:
      BC_Types[iBC].table->f_name = "S33_4";
      BC_Types[iBC].table->f_index = POLYMER_STRESS33_4;
      BC_Types[iBC].desc->equation = R_STRESS33_4;
      BC_Types[iBC].desc->sens[POLYMER_STRESS33_4] = 1;
      break;
   case S11_5_ABS:
      BC_Types[iBC].table->f_name = "S11_5";
      BC_Types[iBC].table->f_index = POLYMER_STRESS11_5;
      BC_Types[iBC].desc->equation = R_STRESS11_5;
      BC_Types[iBC].desc->sens[POLYMER_STRESS11_5] = 1;
      break;
   case S12_5_ABS:
      BC_Types[iBC].table->f_name = "S12_5";
      BC_Types[iBC].table->f_index = POLYMER_STRESS12_5;
      BC_Types[iBC].desc->equation = R_STRESS12_5;
      BC_Types[iBC].desc->sens[POLYMER_STRESS12_5] = 1;
      break;
   case S22_5_ABS:
      BC_Types[iBC].table->f_name = "S22_5";
      BC_Types[iBC].table->f_index = POLYMER_STRESS22_5;
      BC_Types[iBC].desc->equation = R_STRESS22_5;
      BC_Types[iBC].desc->sens[POLYMER_STRESS22_5] = 1;
      break;
   case S13_5_ABS:
      BC_Types[iBC].table->f_name = "S13_5";
      BC_Types[iBC].table->f_index = POLYMER_STRESS13_5;
      BC_Types[iBC].desc->equation = R_STRESS13_5;
      BC_Types[iBC].desc->sens[POLYMER_STRESS13_5] = 1;
      break;
   case S23_5_ABS:
      BC_Types[iBC].table->f_name = "S23_5";
      BC_Types[iBC].table->f_index = POLYMER_STRESS23_5;
      BC_Types[iBC].desc->equation = R_STRESS23_5;
      BC_Types[iBC].desc->sens[POLYMER_STRESS23_5] = 1;
      break;
   case S33_5_ABS:
      BC_Types[iBC].table->f_name = "S33_5";
      BC_Types[iBC].table->f_index = POLYMER_STRESS33_5;
      BC_Types[iBC].desc->equation = R_STRESS33_5;
      BC_Types[iBC].desc->sens[POLYMER_STRESS33_5] = 1;
      break;
   case S11_6_ABS:
      BC_Types[iBC].table->f_name = "S11_6";
      BC_Types[iBC].table->f_index = POLYMER_STRESS11_6;
      BC_Types[iBC].desc->equation = R_STRESS11_6;
      BC_Types[iBC].desc->sens[POLYMER_STRESS11_6] = 1;
      break;
   case S12_6_ABS:
      BC_Types[iBC].table->f_name = "S12_6";
      BC_Types[iBC].table->f_index = POLYMER_STRESS12_6;
      BC_Types[iBC].desc->equation = R_STRESS12_6;
      BC_Types[iBC].desc->sens[POLYMER_STRESS12_6] = 1;
      break;
   case S22_6_ABS:
      BC_Types[iBC].table->f_name = "S22_6";
      BC_Types[iBC].table->f_index = POLYMER_STRESS22_6;
      BC_Types[iBC].desc->equation = R_STRESS22_6;
      BC_Types[iBC].desc->sens[POLYMER_STRESS22_6] = 1;
      break;
   case S13_6_ABS:
      BC_Types[iBC].table->f_name = "S13_6";
      BC_Types[iBC].table->f_index = POLYMER_STRESS13_6;
      BC_Types[iBC].desc->equation = R_STRESS13_6;
      BC_Types[iBC].desc->sens[POLYMER_STRESS13_6] = 1;
      break;
   case S23_6_ABS:
      BC_Types[iBC].table->f_name = "S23_6";
      BC_Types[iBC].table->f_index = POLYMER_STRESS23_6;
      BC_Types[iBC].desc->equation = R_STRESS23_6;
      BC_Types[iBC].desc->sens[POLYMER_STRESS23_6] = 1;
      break;
   case S33_6_ABS:
      BC_Types[iBC].table->f_name = "S33_6";
      BC_Types[iBC].table->f_index = POLYMER_STRESS33_6;
      BC_Types[iBC].desc->equation = R_STRESS33_6;
      BC_Types[iBC].desc->sens[POLYMER_STRESS33_6] = 1;
      break;
   case S11_7_ABS:
      BC_Types[iBC].table->f_name = "S11_7";
      BC_Types[iBC].table->f_index = POLYMER_STRESS11_7;
      BC_Types[iBC].desc->equation = R_STRESS11_7;
      BC_Types[iBC].desc->sens[POLYMER_STRESS11_7] = 1;
      break;
   case S12_7_ABS:
      BC_Types[iBC].table->f_name = "S12_7";
      BC_Types[iBC].table->f_index = POLYMER_STRESS12_7;
      BC_Types[iBC].desc->equation = R_STRESS12_7;
      BC_Types[iBC].desc->sens[POLYMER_STRESS12_7] = 1;
      break;
   case S22_7_ABS:
      BC_Types[iBC].table->f_name = "S22_7";
      BC_Types[iBC].table->f_index = POLYMER_STRESS22_7;
      BC_Types[iBC].desc->equation = R_STRESS22_7;
      BC_Types[iBC].desc->sens[POLYMER_STRESS22_7] = 1;
      break;
   case S13_7_ABS:
      BC_Types[iBC].table->f_name = "S13_7";
      BC_Types[iBC].table->f_index = POLYMER_STRESS13_7;
      BC_Types[iBC].desc->equation = R_STRESS13_7;
      BC_Types[iBC].desc->sens[POLYMER_STRESS13_7] = 1;
      break;
   case S23_7_ABS:
      BC_Types[iBC].table->f_name = "S23_7";
      BC_Types[iBC].table->f_index = POLYMER_STRESS23_7;
      BC_Types[iBC].desc->equation = R_STRESS23_7;
      BC_Types[iBC].desc->sens[POLYMER_STRESS23_7] = 1;
      break;
   case S33_7_ABS:
      BC_Types[iBC].table->f_name = "S33_7";
      BC_Types[iBC].table->f_index = POLYMER_STRESS33_7;
      BC_Types[iBC].desc->equation = R_STRESS33_7;
      BC_Types[iBC].desc->sens[POLYMER_STRESS33_7] = 1;
      break;  
  }
  
  switch (interpolation)
  {
    case LINEAR_INT:
      BC_Types[iBC].table->interp_method = LINEAR;
      break;
    case QUADRATIC_INT:
      BC_Types[iBC].table->interp_method = QUADRATIC;
      break;
    case QUAD_GP_INT:
      BC_Types[iBC].table->interp_method = QUAD_GP;
      break;
    case  BIQUADRATIC_INT:
      BC_Types[iBC].table->interp_method = BIQUADRATIC;
      /*   if( BC_Types[iBC].BC_Name == TABLE_WICV_BC ) {BC_Types[iBC].table->columns = 5;}
         if( BC_Types[iBC].BC_Name == TABLE_WICS_BC ) {BC_Types[iBC].table->columns = 3;}*/ /*???? */
      break;
   }
BC_Types[iBC].table_index = num_BC_Tables++;  
return(1);
}

struct Data_Table *
setup_mp_table 
  ( 
  struct Data_Table * table,
  int    num_const,		/* this is the number of table columns */
  char * dependent_variable,
  char * first_independent_variable,
  char * second_independent_variable,
  char * third_independent_variable,
  int    interpolation_scheme,
  int 	 species_number  
  )
{
int	i,j;
char * independent_names[3];


/* Dependent Variable Porcessing */
table->f_name = dependent_variable;
table->f_index = *dependent_variable;

/* Number of Columns processing */

if( (num_const > 3) && (strcmp( table->f_name, "Saturation") != 0))
{
  sprintf( msg, "Multi DOF table lookup limited to bilinear");
  declare_error(msg);
}
if( (num_const > 4) && (strcmp( table->f_name, "Saturation") == 0))
{
  sprintf( msg, "Multi DOF table lookup limited to bilinear for Saturation");
  declare_error(msg);
}
table->columns = columns_expected= num_const;

/* Independent Variable Processing */

independent_names[0] = first_independent_variable;
independent_names[1] = second_independent_variable;
independent_names[2] = third_independent_variable;
for(i=0;i<table->columns-1;i++)
{
  if ( (strcmp( independent_names[i], "NA") == 0) )
  {
  if( (strcmp( independent_names[i], "TEMPERATURE") == 0) || (strcmp( independent_names[i], "T") == 0))
  {
    strcpy( table->t_name[i],"TEMPERATURE");
    table->t_index[i] = TEMPERATURE;
    for(j=0;j<i;j++)
    {
      if(strcmp(table->t_name[j],table->t_name[i])==0)
      {
        sprintf (msg, "Cannot set Temperature multiple times as Independent Variable.");
        declare_error(msg);	
      }
    }
  }
  else if ( (strcmp( independent_names[i], "MASS_FRACTION") == 0) || (strcmp( independent_names[i], "SPECIES") == 0) ||
		  ( strcmp( independent_names[i], "Y") == 0 ))
  {
    strcpy( table->t_name[i],"MASS_FRACTION");
    table->t_index[i] = MASS_FRACTION;
    if (species_number == NA)
    {
      declare_error("A species number is required for this independent variable.");
    }
    else
    {
      table->species_eq = species_number;
    }
  }
  else if( (strcmp( independent_names[i], "CAP_PRES") == 0) )
  {
    if( (strcmp( table->f_name, "Saturation") != 0) ) 
    {
      declare_error("Independent variable CAP_PRES can only be used with Saturation.");
    }
    strcpy( table->t_name[i],"CAP_PRES");
    table->t_index[i] = CAP_PRES;
    for(j=0;j<i;j++)
    {
      if(strcmp(table->t_name[j],table->t_name[i])==0)
      {
        declare_error("Cannot set CAP_PRES multiple times as Independent Variable.");	
      }
    }
  } 
  
  if(i == 0) { strcpy( table->t_name[i+1],"NULL"); table->t_index[i+1]=-1;}
  }  /* end if */
} /* end do for */

/* Interpolation scheme processing */
if( interpolation_scheme == LINEAR)
{
  if( (table->columns == 2) || (table->columns == 3 && (strcmp( table->f_name, "Saturation") == 0)))
  {
    table->interp_method = LINEAR;
  }
  else 
  {
    declare_error("Incorrect number of columns for material property table lookup.");
  }
}
else if( interpolation_scheme == BILINEAR)
{
  if( table->columns >= 3)
  {
    table->interp_method = BILINEAR;
  }
  else 
  {
    declare_error("Incorrect number of columns for material property table lookup.");
  }
}
mp_table_species_number = NA;
return(table);
}

int read_array_into_table
  (
  enum table_types Table_Type,
  struct Data_Table * table,
  int Num_Pnts /* Num_Pnts is the same as the numb_lines_read  */
  )
{
  double p,p2,p3,p4;
  int  i,j,k,ibegin,iend;
  char msg[100];
  int havent_declared_an_error = 1;
  
  
    table->tablelength = Num_Pnts;
    if( table->tablelength == 0 )
    {
      fprintf(parser_log, "<<< ERROR %i: Error reading tabular data. Can't find any points >>>", error_number++);
      line_number++;
      /*cards[END_TABLE_CARD]->times_read[0] = cards[END_TABLE_CARD]->times_read[0] + 1; jjj */
      if (ProcID == 0) fprintf(parser_log,"\n%d ",line_number );
      error_found_in_last_line = 0;
      return(-1);
    }
    if( table->tablelength == 1 ) 
    {
      fprintf(parser_log, "<<< ERROR %i: Error reading tabular data . Need more than 1 point >>>", error_number++);
      line_number++;
      /* cards[END_TABLE_CARD]->times_read[0] = cards[END_TABLE_CARD]->times_read[0] + 1; jjj */
      error_found_in_last_line = 0;
      if (ProcID == 0) fprintf(parser_log,"\n%d ",line_number );
      return(-1);
    }
    if( (table->interp_method == QUADRATIC || table->interp_method == BIQUADRATIC) && Num_Pnts%2 != 1)
    {
      fprintf(parser_log, "<<< ERROR %i: Need odd number of points for (bi)quadratic interpolation >>>", error_number++);
      line_number++;
      /* cards[END_TABLE_CARD]->times_read[0] = cards[END_TABLE_CARD]->times_read[0] + 1; jjj */
      error_found_in_last_line = 0;
      if (ProcID == 0) fprintf(parser_log,"\n%d ",line_number );
      return(-1);
    }
    if(table->interp_method == QUAD_GP && Num_Pnts%3 != 0)
    {
      fprintf(parser_log, "<<< ERROR %i: Need 3N points for quadratic gauss point interpolation >>>",error_number++);
      line_number++;
      /* cards[END_TABLE_CARD]->times_read[0] = cards[END_TABLE_CARD]->times_read[0] + 1; */
      error_found_in_last_line = 0;
      if (ProcID == 0) fprintf(parser_log,"\n%d ",line_number );
      return(-1);
    } 
    /* Assign memory */
    table->t = (double *) smalloc( sizeof(double)*Num_Pnts );
    table->t2 = (double *) smalloc( sizeof(double)*Num_Pnts );
        if(table->columns == 5) {
    table->f = (double *) smalloc( sizeof(double)*Num_Pnts*3 );}
        else if (table->columns == 3 && table->interp_method != BIQUADRATIC) {
    table->f = (double *) smalloc( sizeof(double)*Num_Pnts*2 );}
        else    {
    table->f = (double *) smalloc( sizeof(double)*Num_Pnts );}
    /* if !BIQUAD */
    /*   How about if we do some error checking here
    i.e. check for duplicate abscissa values in 1D tables  */
    if( table->columns == 2)
    {
      for (k=0;k<Num_Pnts;k++)
      {
        table->t[k] = table_array[k][0]; 
        table->f[k] = table_array[k][1]; 
      }
    } 
    if(table->columns ==3 && (table->interp_method == BIQUADRATIC || table->interp_method == BILINEAR) )
    {  
      for (k=0;k<Num_Pnts;k++)
      {
        table->t[k]  = table_array[k][0];
        table->t2[k] = table_array[k][1];
        table->f[k]  = table_array[k][2];
      }                 
    }
    if(table->columns ==3 && (table->interp_method != BIQUADRATIC && table->interp_method != BILINEAR) )
    {
      for (k=0;k<Num_Pnts;k++)
      {
        table->t[k]  = table_array[k][0];
        table->f[k]  = table_array[k][1];
        table->f[k+Num_Pnts]  = table_array[k][2];      /*?????*/
     }        
    } 
    if(table->columns ==5)
    {
      for (k=0;k<Num_Pnts;k++)
      {
        table->t[k]  = table_array[k][0];
        table->t2[k] = table_array[k][1];
        table->f[k]  = table_array[k][2];
        table->f[k+Num_Pnts] = table_array[k][3];
        table->f[k+2*Num_Pnts]  = table_array[k][4];      
      }     
    }
    /* 
    * Now sort the points from lowest to highest abscissa (Because people forget)
    */
    if(table->interp_method != BIQUADRATIC)  /* don't sort 2d bc tables */
    {
      for ( i = 1; i<Num_Pnts; i++)
      {
        p=table->t[i];
        p2=table->t2[i];
        p3=table->f[i];
        if(table->columns == 3 && table->interp_method != BILINEAR)
          {p4=table->f[i+Num_Pnts];}
        j=i-1;
        while (j>=0 && (p<table->t[j]) )
        {
          table->t[j+1]=table->t[j]; table->t[j]=p;
          table->t2[j+1]=table->t2[j]; table->t2[j]=p2;
          table->f[j+1]=table->f[j]; table->f[j]=p3;
          if(table->columns == 3 && table->interp_method != BILINEAR)
          { 
            table->f[j+1+Num_Pnts]=table->f[j+Num_Pnts];
            table->f[j+Num_Pnts]=p4;
          }
          if(( table->t[i] == table->t[j]) && table->columns == 2) 
          {
            if (havent_declared_an_error)
            {
              havent_declared_an_error = 0;
              fprintf(parser_log, "<<< ERROR %i: Multivalued function detected in table :%s >>>", error_number,table->f_name);
              error_number++;
            }  
          }
          j--;
        }
      }	      
      /* If 3 Columns then sort 2nd Column */
      if(table->columns == 3 && table->interp_method == BILINEAR)
      {
        ibegin=1;
        for(k=1;k<Num_Pnts;k++)
        {
          if((table->t[k-1] != table->t[k]) || (k==Num_Pnts-1) )
          {
            if(k==Num_Pnts-1)
	    {
	      iend=Num_Pnts;
	    }
	    else
	    {
	      iend=k;
	    }
      	    for(i=ibegin;i<iend; i++)
	    {
	      p2=table->t2[i];
	      p3=table->f[i];
	      j=i-1;
	      while (j>=ibegin-1 && (p2<table->t2[j]) )
	     {
	        table->t2[j+1]=table->t2[j]; table->t2[j]=p2;
	        table->f[j+1]=table->f[j]; table->f[j]=p3;
	        j--;
	      }
            }
	    ibegin=iend+1;
          }
        }
      }
    } /* How about if we do some error checking here, i.e., check for cuplicate abscissa values in 1D tables */
    if(table->interp_method != BIQUADRATIC && table->interp_method != BILINEAR)
    {
      for( i=1; i < Num_Pnts ; i++)
      {  
        if(table->t[i] == table->t[i-1])
        { 
            fprintf(parser_log, "<<< ERROR %i:  Multivalued function detected in table >>>", error_number);
            error_number++;
        }
      }
    }
    /*  determine number of grid points for 2d tables */
    if(table->interp_method == BIQUADRATIC)
    {
      double a, b, c, cosineC;
      for ( i=2 ; i < Num_Pnts ; i++ )
      {
        a = SQUARE(table->t[0]-table->t[i-1]) +
          SQUARE(table->t2[0]-table->t2[i-1]);
        b = SQUARE(table->t[0]-table->t[i]) +
          SQUARE(table->t2[0]-table->t2[i]);
        c = SQUARE(table->t[i]-table->t[i-1]) +
         SQUARE(table->t2[i]-table->t2[i-1]);
        cosineC = (a+b-c)/(2.*sqrt(a*b));
        if(fabs(cosineC) < 0.7)
        {
          table->ngrid = i;
          break;
        }
      }
      if(Num_Pnts%table->ngrid != 0)
      {
        fprintf(parser_log,"<<< ERROR %i: 2D table is not a rectangular grid3 (%d x %d != %d points) >>>",error_number,Num_Pnts/table->ngrid ,table->ngrid,Num_Pnts);
        error_number++;
      }  
    }
    line_number++;
    /* cards[END_TABLE_CARD]->times_read[0] = cards[END_TABLE_CARD]->times_read[0] + 1; jjj */
    if (ProcID == 0) fprintf(parser_log,"\n%d ",line_number );
    number_of_table_cards_read = 0;

  
  if ( table_type == BC_TABLE )
  {
  
  }
  
  if ( (table_type == CONDUCTIVITY_TABLE) || 
       (table_type == HEAT_CAPACITY_TABLE)|| 
       (table_type == VISCOSITY_TABLE) ||
       (table_type == DIFFUSIVITY_TABLE) ||
       (table_type == SPECIES_TABLE)
     )
  {
    num_MP_Tables++;
  } /* end of Material Property processing */  
  
  if (table_type == NO_TABLE)
  {
  
  } /* end of NO_TABLE processing */
  
  accept_table_data_cards = FALSE;
  table_type = NO_TABLE;
}

int common_BC_processing(char * BC_text_name, char * BC_set_type, int first_BC_ID)
{
  int k,l;
  int bc_found = 0;
  int bc_OK = 1;
  int description_index;
  if ( !accept_bcs_until_end_of_bc_card && (number_of_bcs_found == number_of_bcs_expected) ) return (0);
  for (k = 0; k < Num_BC_Names; k++)
  {
    if ( !strcmp(BC_text_name, BC_Desc[k].name1) || !strcmp(BC_text_name, BC_Desc[k].name2) )
    {
      bc_found = 1;
      description_index = k;
    }
  }
  if( !bc_found )
  {
    declare_error("Boundary condition not found in BC_Desc (mm_names.h file).");
    return(0);
  }
  /* Check that BC is appropriate to SS */
  if( !strcmp( BC_set_type, "SS") || !strcmp( BC_set_type, "SC") )
  {
    switch ( BC_Desc[description_index].method )
    {
      case   COLLOCATE_SURF:
      case    WEAK_INT_SURF:
      case  STRONG_INT_SURF:
      case  STRONG_INT_EDGE:
      case    WEAK_INT_EDGE:
      case   COLLOCATE_EDGE:
        break;
      default:
        declare_error("BC Consistency error detected, side set type should be NS or NC.");
	return(0);
    }
  }
  if( !strcmp( BC_set_type, "NS") || !strcmp( BC_set_type, "NC") )
  {
    switch ( BC_Desc[description_index].method )
    {
      case   DIRICHLET:
      case   SPECIAL:
        break;
     default:
       declare_error("BC Consistency error detected, side set type should be SS or SC.");
       return(0);
    }
  }
  /* for (l=0; l<=iBC; l++)
  {  
    if( ( BC_Desc[description_index].BC_Name == BC_Types[l].BC_Name) && (BC_Types[l].BC_ID == first_BC_ID) )
    {
      bc_OK = 0;
    }
  }
  if( !bc_OK )
  {
    declare_warning("Prior BC of this type has same side set.");
    return(0);
  } */
  if (number_of_bcs_found == 0)
  {
    BC_Types = alloc_struct_1(struct Boundary_Condition, 1);
  }
  else
  {
    BC_Types = realloc( BC_Types, (number_of_bcs_found + 1) * sizeof( struct Boundary_Condition ) );
  }
  iBC++;
  number_of_bcs_found++;
  initialize_Boundary_Condition(BC_Types + iBC);
  BC_Types[iBC].BC_Name       = BC_Desc[description_index].BC_Name;
  BC_Types[iBC].desc          = &BC_Desc[description_index];
  BC_Types[iBC].BC_Desc_index = k;
  if( (!strcmp("NC",BC_set_type)) || (!strcmp("NS",BC_set_type)) )
    {
      strncpy(BC_Types[iBC].Set_Type, "NS", 3);
    }
    else
    {
      strncpy(BC_Types[iBC].Set_Type, "SS", 3);
    }
    if( (!strcmp("NC",BC_set_type)) || (!strcmp("SC",BC_set_type)) )
    {
      cont->upBCID = iBC;
    }
    BC_Types[iBC].BC_ID = first_BC_ID;
    /* if ( !accept_bcs_until_end_of_bc_card ) number_of_bcs_expected--; */
    return(number_of_bcs_found);
}

int process_equation_and_variable_names(char * equation_name, char * variable_name)
{
  int k;
  int eq_found = 0;
  int var_found = 0;
  char msg[100];
  
  for(k=0; k<Num_EQ_Names; k++) /* equation names loop */
  {
    if ( !strcmp(equation_name, EQ_Name[k].name1) || !strcmp(equation_name, EQ_Name[k].name2) )
    {
      eq_found = 1;
      eq_index = EQ_Name[k].Index;
    }
  }
  for(k=0; k<Num_Var_Names; k++)  /* variable names loop */
  {
    if ( !strcmp(variable_name, Var_Name[k].name1) || !strcmp(variable_name, Var_Name[k].name2) )
    {
      var_found = 1;
      var_index = Var_Name[k].Index;
    }
  } 
  /* special case for GD_TIME card since GD_TIME names won't be in variable list */
  if ( !var_found)
  { 
    if ( !strcmp(variable_name, "LINEAR") )
    {
      BC_Types[iBC].BC_Data_Int[2] = GD_TIME_LIN;
      var_found = TRUE;
    } 
    else if ( !strcmp(variable_name, "EXPONENTILAL") )
    {
      BC_Types[iBC].BC_Data_Int[2] = GD_TIME_EXP;
      var_found = TRUE;
    } 
    else if ( !strcmp(variable_name, "SINUSOIDAL") )
    {
      BC_Types[iBC].BC_Data_Int[2] = GD_TIME_SIN;
      var_found = TRUE;
    } 
    else if ( !strcmp(variable_name, "TABLE") )
    {
      BC_Types[iBC].BC_Data_Int[2] = GD_TIME_TABLE;
      var_found = TRUE;
    } 
  }  
  if ( !var_found && !eq_found) {sprintf(msg, "Invalid equation (%s) and variable (%s) names.", equation_name, variable_name); declare_error(msg); return 0;}
  if ( !eq_found  ) {sprintf(msg, "Invalid equation name (%s).", equation_name); declare_error(msg); return 0;}
  if ( !var_found ) {sprintf(msg, "Invalid variable name (%s).", variable_name); declare_error(msg); return 0;}
  return 1;
}

int valid_fix_bc_variable_name ( char * variable_name)
{
  int k;
  int var_found = 0;
  char msg[100];
  for(k=0; k<MAX_VARIABLE_TYPES+1; k++) /* variable names loop */
  {
    if ( !strcmp(variable_name, Var_Name[k].name1) || !strcmp(variable_name, Var_Name[k].name2) )
    {
      BC_Types[iBC].BC_Data_Int[0] = Var_Name[k].Index;
      var_found = 1;
    }
  }
  if ( !var_found  ) {sprintf(msg, "Invalid variable name (%s).", variable_name); declare_error(msg); return 0;}
  return 1;
}

int allocate_fix_bc_description( )
{
  new_BC_Desc = ((struct BC_descriptions  **)
	 realloc(new_BC_Desc,  
        	 (num_new_BC_Desc+1) * 
	 sizeof(struct BC_descriptions *)));

  new_BC_Desc[num_new_BC_Desc] = 
	    alloc_BC_description(BC_Types[iBC].desc);
	  
  BC_Types[iBC].desc           = new_BC_Desc[num_new_BC_Desc];
  BC_Types[iBC].index_dad      = num_new_BC_Desc;
  num_new_BC_Desc++;
  BC_Types[iBC].desc->equation = BC_Types[iBC].BC_Data_Int[0];
  BC_Types[iBC].desc->sens[BC_Types[iBC].BC_Data_Int[0]] = 1;
  return(1);
}

int allocate_table_bcs_description( )
{
  new_BC_Desc = ((struct BC_descriptions  **)
	 realloc(new_BC_Desc,  
        	 (num_new_BC_Desc+1) * 
	 sizeof(struct BC_descriptions *)));

  new_BC_Desc[num_new_BC_Desc] = 
	    alloc_BC_description(BC_Types[iBC].desc);
	  
  BC_Types[iBC].desc           = new_BC_Desc[num_new_BC_Desc];
  BC_Types[iBC].index_dad      = num_new_BC_Desc++;  /* This is improtant to Phil */
  return(1);
}

int allocate_GD_bc_description( )
{	
  int var, eqn;
  /*
   * Need to Allocate Space for the BC Description for the
   * current boundary condition
   * 
   *  We do this by reallocating the pointer array to add an
   *  extra pointer onto the end of the current list. Then
   *  we malloc and initialize a new BC description structure
   *  and add it to the end of the list.
   */
  new_BC_Desc = (struct BC_descriptions  **)
	      realloc(new_BC_Desc,  
		      (num_new_BC_Desc+1) * sizeof(struct BC_descriptions *));
  new_BC_Desc[num_new_BC_Desc] = 
	      alloc_BC_description(BC_Types[iBC].desc);
  /*
   *  Now fill in the new BC_Types condition and BC_description
   *  structure with information about the boundary condition
   */
  BC_Types[iBC].desc           = new_BC_Desc[num_new_BC_Desc];
  BC_Types[iBC].index_dad      = num_new_BC_Desc;
  BC_Types[iBC].desc->equation = BC_Types[iBC].BC_Data_Int[0];
  var                          = BC_Types[iBC].BC_Data_Int[2];
  num_new_BC_Desc++;
  /*
   * Fill in an a default description of the sensitivity of this
   * boundary condition to the independent variables (i.e., what
   * has a nonzero Jacobian entry).
   *
   *  Notes:
   *    MESH_POSITION variables will point to MESH_DISPLACEMENT eqns.
   *    Time derivatives of regular variables will point to the
   *      associated equation.
   *  Question:
   *    What if you need more ?
   */
  if (var >= V_FIRST && var < V_LAST)
  {
    BC_Types[iBC].desc->sens[var] = 1;
  }
  else if (var >= MESH_POSITION1 && var <= MESH_POSITION3)
  {
    BC_Types[iBC].desc->sens[var - MESH_POSITION1 + MESH_DISPLACEMENT1] = 1;
  }
  else if (var >= D_VEL1_DT && var <= D_P_DT)
  {
    BC_Types[iBC].desc->sens[var - D_VEL1_DT] = 1;
  }
  /* Check to see if this condition implies rotation of mesh or momentum equations */
  eqn = BC_Types[iBC].BC_Data_Int[0];
  if ( (eqn >= R_MESH_NORMAL) && (eqn <= R_MESH_TANG2) )
  {
    BC_Types[iBC].desc->rotate = R_MESH1;
  }
  if ( (eqn >= R_MOM_NORMAL) && (eqn <= R_MOM_TANG2) )
  {
    BC_Types[iBC].desc->rotate = R_MOMENTUM1;
  }
  return(1);
}

int set_equation
  (
  char * card,
  int mn,    
  int eqnNum, 
  int Galerkin_wt, 
  int Variable_name, 
  int Interpolation_fnc, 
  float Mass_const, 
  float Advective_const, 
  float Boundary_const, 
  float Diffusive_const, 
  float Source_const, 
  float Porous_const
  )
/******************************************************************
 *
 * set_equation(
 *              equation name,
 *		material number, 
 *      	Galerkin_wt, 
 *		Variable_name, 
 *		Interpolation_fnc, 
 *      	Mass_const, 
 *		Advective_const, 
 *		Boundary_const, 
 *		Diffusive_const, 
 *      	Source_const, 
 *		Porous_const)
 *
 ******************************************************************/
{ 
  if (pd_glob[mn] != NULL)
  { 
  if( ( (number_of_equations_found[mn] < number_of_equations_expected[mn]) || accept_eqs_till_end_of_eq_card ) && accept_material_related_cards ) 
  {
  number_of_equations_found[mn]++;
  if (pd_glob[mn]->e[eqnNum] == T_SOMETHING) 
    {
    declare_error("Equation has already been activated.");
    } else {
    pd_glob[mn]->e[eqnNum] = T_SOMETHING;
    pd_glob[mn]->w[eqnNum] = Galerkin_wt;
    if (pd_glob[mn]->v[Variable_name]) 
      {
      if (ProcID == 0) 
      {
        fprintf(parser_log, "\n\n<<< ERROR %d Line %d: Variable has already been activated. >>>\n", error_number, line_number);
        fprintf(error_log, "\n\n<<< ERROR %d Line %d: Variable has already been activated. >>>\n", error_number, line_number);        
        fprintf(stdout, "\n\n<<< ERROR %d Line %d: Variable has already been activated. >>>\n", error_number, line_number);        
      }  
      error_number++;
      } else {      
      pd_glob[mn]->v[Variable_name] |= V_SOLNVECTOR;
      pd_glob[mn]->i[eqnNum] = Interpolation_fnc;
      if( Mass_const      != NA)  pd_glob[mn]->etm[eqnNum][(LOG2_MASS)]          = Mass_const;		
      if( Advective_const != NA)  pd_glob[mn]->etm[eqnNum][(LOG2_ADVECTION)]     = Advective_const;
      if( Boundary_const  != NA)  pd_glob[mn]->etm[eqnNum][(LOG2_BOUNDARY)]      = Boundary_const;
      if( Diffusive_const != NA)  pd_glob[mn]->etm[eqnNum][(LOG2_DIFFUSION)]     = Diffusive_const;
      if( Source_const    != NA)  pd_glob[mn]->etm[eqnNum][(LOG2_SOURCE)]        = Source_const;
      if( Porous_const    != NA)  pd_glob[mn]->etm[eqnNum][(LOG2_POROUS_BRINK)]  = Porous_const;	
      card_read(card, mn);
      }
     }
    } 
    else 
    {
      declare_warning("This equation card is being ignored.");	
    }
  }
  else
  {
    line_number++;
    if (ProcID == 0)  fprintf(parser_log,"\n%d ",line_number );
  }    
  /*return (eqnNum);*/
  return(1);
}

int set_variable
  (
  int varType, 
  int mn
  )
/******************************************************************
 *
 * set_variable()
 *
 *   Utility function used repeatedly as a kernal algorithm
 *   in rd_eq_specs.
 ******************************************************************/
{
  if (pd_glob[mn]->v[varType]) 
    {
    if (ProcID == 0) fprintf(parser_log, "\n\n<<< ERROR %d Line %d: Variable has already been activated. >>>\n", error_number, line_number);
    error_number++;
    }  else {
    pd_glob[mn]->v[varType] |= V_SOLNVECTOR;
    }
  return varType;
}

int create_new_material
  (
  char *name, 
  int Number_Matrl_Elem_Blks,
  int ID1,
  int ID2,
  int ID3
  )
{
  int i,j;
  mn++;
 /* pd_glob[mn] = (struct Problem_Description *) malloc(sizeof(struct Problem_Description)); */
  /*mp_glob[mn] = (struct Material_Properties *) malloc(sizeof(struct Material_Properties));*/

  mp_glob[mn]->Num_Matrl_Elem_Blk = Number_Matrl_Elem_Blks;
  mp_glob[mn]->Matrl_Elem_Blk_Ids = alloc_int_1(mp_glob[mn]->Num_Matrl_Elem_Blk, 0);

  for ( i=0; i<MAX_EQNS; i++ )
    {
    pd_glob[mn]->e[i]  = T_NOTHING;	/* No terms active in this equation. */
    pd_glob[mn]->v[i]  = V_NOTHING;	/* No active variables, either. */
    pd_glob[mn]->w[i]  = I_NOTHING;	/* No weighting function. */
    pd_glob[mn]->i[i]  = I_NOTHING;	/* Nothing to interpolate. */
    pd_glob[mn]->m[i]  = -1;		/* Map from i to R_EQNTYPE initially undefd. */
    for ( j=0; j<MAX_TERM_TYPES; j++)
      {
      pd_glob[mn]->etm[i][j] = 0.;	/* Zero out all terms. */
      }
    }
  /* Initialized Explicit_Fill to zero. Only set it if
   * we have a fill equation in the loop 
   */  
  Explicit_Fill = 0;
  strcpy(pd_glob[mn]->MaterialName, name);
  strcpy(mp_glob[mn]->Material_Name, name); 
  number_of_materials_found++;
  accept_material_related_cards = TRUE;
  switch (Number_Matrl_Elem_Blks)
  {
    case 1:
      mp_glob[mn]->Matrl_Elem_Blk_Ids[0] = ID1;
      break;
    case 2:
      mp_glob[mn]->Matrl_Elem_Blk_Ids[0] = ID1;
      mp_glob[mn]->Matrl_Elem_Blk_Ids[1] = ID2;
      break;
    case 3:
      mp_glob[mn]->Matrl_Elem_Blk_Ids[0] = ID1;
      mp_glob[mn]->Matrl_Elem_Blk_Ids[1] = ID2;      
      mp_glob[mn]->Matrl_Elem_Blk_Ids[2] = ID3;
      break;
  }
  return(mn);
}



void card_read
  (
  char * card, 
  int file_index
  )
{ 
  int i;
  for(i=0; i<num_cards; i++ )
  {
    if( !strcmp(card,cards[i].card_name) )
    {
    cards[i].times_read_per_file[file_index]++;
   /* fprintf(stdout,"\n %s file_index=%d mn=%d     %d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d", 
    	card, 
    	file_index,
    	mn,
    	cards[i].times_read_per_file[0],
    	cards[i].times_read_per_file[1],
    	cards[i].times_read_per_file[2],
    	cards[i].times_read_per_file[3],
    	cards[i].times_read_per_file[4],
    	cards[i].times_read_per_file[5],
    	cards[i].times_read_per_file[6],
    	cards[i].times_read_per_file[7],
    	cards[i].times_read_per_file[8],
    	cards[i].times_read_per_file[9],
    	cards[i].times_read_per_file[10],
    	cards[i].times_read_per_file[11]);*/
    }
  }
  floating_point_constant_list_index = -1; 
  /*line_number++;*/
  /*Uif (ProcID == 0)  fprintf(parser_log,"\n%d ",line_number ); */
} 

int times_card_read
  (
  char * card, 
  int file_index
  )
{ 
  int i;
  for(i=0; i<num_cards; i++ )
  {
    if( !strcmp(card,cards[i].card_name) )
    {
    return cards[i].times_read_per_file[file_index];
    }
    else
    {
      return (-1);
    }
  } 
} 

int common_AC_processing()
{  
  if ( !accept_acs_until_end_of_ac_card && (number_of_acs_found == number_of_acs_expected) ) return (0);
  if ( number_of_acs_found == 0)
  {
    augc = (struct AC_Information *) alloc_struct_1(struct AC_Information,1);
  }
  else
  {
    augc = realloc(augc, ( (number_of_acs_found+1) * sizeof (struct AC_Information) ) );
  }
  iAC++;
  number_of_acs_found++;
  augc[iAC].nAC = number_of_acs_found;
  augc[iAC].BCID = -99;
  augc[iAC].MTID = -99;
  augc[iAC].VOLID = -99;
  augc[iAC].MFID = -99;
  augc[iAC].iread=0;
  augc[iAC].len_AC = 0;
  if(AC_rd_file == TRUE) {augc[iAC].iread=1;}

  return(number_of_acs_found);
}

int common_CC_processing( int type, int first_integer, int second_integer, int third_integer, float first_float, float second_float )
{
  double range, vfloat;
  if ( !accept_ccs_until_end_of_cc_card && (number_of_ccs_found == number_of_ccs_expected) ) return (0);
  if ( number_of_ccs_found == 0)
  {
    cpcc = (struct Continuation_Conditions *) alloc_struct_1(struct Continuation_Conditions,1);
  }
  else
  {
    cpcc = realloc(cpcc, ( (number_of_ccs_found+1) * sizeof (struct Continuation_Conditions) ) );
  }
  iCC++;
  cpcc[iCC].nCC = nCC;
  cpcc[iCC].ratio = 1.0;
  cpcc[iCC].Type = cont->upType;
  cpcc[iCC].BCID = cont->upBCID;
  cpcc[iCC].DFID = cont->upDFID;
  cpcc[iCC].MTID = cont->upMTID;
  cpcc[iCC].MPID = cont->upMPID;
  cpcc[iCC].MDID = cont->upMDID;
  cpcc[iCC].Beg_CC_Value = cont->BegParameterValue;
  cpcc[iCC].End_CC_Value = cont->EndParameterValue;
  range = cpcc[iCC].End_CC_Value - cpcc[iCC].Beg_CC_Value;
  cpcc[iCC].Type = type;
  switch( type)
  {
    case 1: /* BC */
      cpcc[iCC].BCID = first_integer;
      cpcc[iCC].DFID = second_integer;
      /* fprintf(stdout, " %3d. BC: BCID=%3d DFID=%5d\n",iCC+1,
                cpcc[iCC].BCID, cpcc[iCC].DFID); */
      break;

     case 2: /* MT */
       cpcc[iCC].MTID = first_integer;
       cpcc[iCC].MPID = second_integer;
       /* fprintf(stdout, " %3d. MT: MTID=%3d MPID=%5d\n",iCC+1,
                cpcc[iCC].MTID, cpcc[iCC].MPID); */
       cpcc[iCC].MTID--;
       break;
  }    
  /*
  *  Third int input third_float indicates meaning of second float input "vfloat":
  *  0 - Value is same as continuation parameter, no floats needed
  *  1 - vfloat is End_CC_Value[iCC]; use to find range ratio d(val[iCC])/d(val[0])
  *  2 - vfloat is range ratio; use to find End_CC_Value[iCC]
  */
  if (third_integer == 0)
  {
    cpcc[iCC].ratio = 1.0;
    cpcc[iCC].Beg_CC_Value = cpcc[0].Beg_CC_Value;
    cpcc[iCC].End_CC_Value = cpcc[0].End_CC_Value;
  }
  else
  {
    switch (third_integer) 
    {
      case 1:
        cpcc[iCC].End_CC_Value = second_float;
        cpcc[iCC].Beg_CC_Value = first_float;
        if ( range != 0.0 )
        {
          cpcc[iCC].ratio = (vfloat - cpcc[iCC].Beg_CC_Value) / range;
        }
        else
        {
          declare_error("CC card is leading to a zero range value & a divide by zero in ratio calculation.");
        }
        number_of_ccs_found++;
        break;
         
      case 2:
        cpcc[iCC].ratio = second_float;
        cpcc[iCC].Beg_CC_Value = first_float;         
        cpcc[iCC].End_CC_Value = cpcc[iCC].Beg_CC_Value + vfloat * range;
        number_of_ccs_found++;
        break;
         
      default:
        sprintf(msg, " Third integer must be 0, 1, or 2."); 
        declare_error(msg);     
        break;
    }
  }  /* End of else block (third_integer == 0) */  
  return(number_of_ccs_found);
}

int common_TP_processing( int type, int first_integer, int second_integer, int third_integer, float first_float, float second_float )
{
  double range, vfloat;
  if ( !accept_tps_until_end_of_tp_card && (number_of_tps_found == number_of_tps_expected) ) return (0);
  if ( number_of_tps_found == 0)
  {
    tpcc = (struct Continuation_Conditions *) alloc_struct_1(struct Continuation_Conditions,1);
  }
  else
  {
    tpcc = realloc(cpcc, ( (number_of_tps_found+1) * sizeof (struct Continuation_Conditions) ) );
  }
  iTP++;
  tpcc[iTP].nCC = nTC;
  tpcc[iTP].ratio = 1.0;
  tpcc[iTP].Type = loca_in->TPupType;
  tpcc[iTP].BCID = loca_in->TPupBCID;
  tpcc[iTP].DFID = loca_in->TPupDFID;
  tpcc[iTP].MTID = loca_in->TPupMTID;
  tpcc[iTP].MPID = loca_in->TPupMPID;
  tpcc[iTP].MDID = loca_in->TPupMDID;
  tpcc[iTP].Beg_CC_Value = loca_in->TPGuess;
  tpcc[iTP].End_CC_Value = loca_in->TPFinal;
  range = tpcc[0].End_CC_Value - tpcc[0].Beg_CC_Value;
  tpcc[iTP].Type = type;
  switch( type)
  {
    case 1: /* BC */
      tpcc[iTP].BCID = first_integer;
      tpcc[iTP].DFID = second_integer;
      /* fprintf(stdout, " %3d. BC: BCID=%3d DFID=%5d\n",iTP+1,
                tpcc[iTP].BCID, tpcc[iTP].DFID); */
      break;

     case 2: /* MT */
       tpcc[iTP].MTID = first_integer;
       tpcc[iTP].MPID = second_integer;
       /* fprintf(stdout, " %3d. MT: MTID=%3d MPID=%5d\n",iTP+1,
                tpcc[iTP].MTID, tpcc[iTP].MPID); */
       tpcc[iTP].MTID--;
       break;
  }    
  /*
  *  Third int input third_float indicates meaning of second float input "vfloat":
  *  0 - Value is same as continuation parameter, no floats needed
  *  1 - vfloat is End_CC_Value[tp_index]; use to find range ratio d(val[iTP])/d(val[0])
  *  2 - vfloat is range ratio; use to find End_CC_Value[iTP]
  */
  if (third_integer == 0)
  {
    tpcc[iTP].ratio = 1.0;
    tpcc[iTP].Beg_CC_Value = tpcc[0].Beg_CC_Value;
    tpcc[iTP].End_CC_Value = tpcc[0].End_CC_Value;
  }
  else
  {
    switch (third_integer) 
    {
      case 1:
        tpcc[iTP].End_CC_Value = second_float;
        tpcc[iTP].Beg_CC_Value = first_float;
        if (range != 0.0)
        { 
          tpcc[iTP].ratio = (vfloat - tpcc[iTP].Beg_CC_Value) / range;
        }
        else
        {
          declare_error("TC card is leading to a zero value of range & a divide by zero in ratio calculation.");
        }        
        number_of_tps_found++;
        break;
         
      case 2:
        tpcc[iTP].ratio = second_float;
        tpcc[iTP].Beg_CC_Value = first_float; 
        tpcc[iTP].End_CC_Value = tpcc[iTP].Beg_CC_Value + vfloat * range;
        number_of_tps_found++;
        break;
         
      default:
        declare_error("Third integer must be 0, 1, or 2.");     
        break;
    }
  }  /* End of else block (third_integer == 0) */  
  return(number_of_tps_found);
}

int yywrap(const char *msg)
{
  return 1;
}
 
void summary_and_closes() 
{ 
  /* Error & Warning summary */
  if (ProcID == 0) 
  {
    if( error_number > 1)
      { /* there are errors */
        fprintf(stdout,"\nTOTAL Parse Warnings: %d \nTOTAL PARSE ERRORS: %d\nSee parser.log or error.log for details.\n\nParse failed.\n\n", warning_number-1, error_number-1);
        fprintf(error_log,"\nTOTAL Parse Warnings: %d \nTOTAL PARSE ERRORS: %d\nSee parser.log for details.\n\nParse failed.\n\n", warning_number-1, error_number-1);
        fprintf(parser_log,"\n\nTOTAL Parse Warnings: %d \nTOTAL PARSE ERRORS: %d\n\nParse failed.\n\n", warning_number-1, error_number-1);
      }
    else
      { /* no errors */
        fprintf(parser_log,"\n\nTOTAL Parse Warnings: %d \nTOTAL PARSE ERRORS: %d\n\nParser succeeded.\n\n", warning_number-1, error_number-1); 
        fprintf(error_log,"\nTOTAL Parse Warnings: %d \nTOTAL PARSE ERRORS: %d\n\nParser succeeded.\n\n", warning_number-1, error_number-1); 
        fprintf(stdout,"\nTOTAL Parse Warnings: %d \nTOTAL PARSE ERRORS: %d\nSee parser.log or error.log for details.\n\nParse succeeded.\n\n", warning_number-1, error_number-1);
      }
  /* Closes and memory deallocation */
  /* free(cards);*/
    fclose(parser_log);
    fprintf(error_log,"\n");
    fclose(error_log);
  }
}

int read_floats(  dbl **target_array, /* array of double constants 
				      * for user mat props 
                                      * target_array[species_no] is the
                                      * pointer to the data to be filled in. */
	       const int species_no) /* species number (zero if no species) */

{
  int np, num_const, i;
  int species_num=-1;
  num_const = floating_point_constant_list_index+1;
  floating_point_constant_list_index = -1;
  if (species_no >= 0) species_num = species_no;
  if (species_no < 0) species_num = 0;
  if (target_array == NULL)
  {
    sprintf(msg, " User Defined not allowed yet."); 
    declare_error(msg);    
    return(0);
  }
  if ( num_const  > 0 ) 
  {
    target_array[species_num] = alloc_dbl_1(num_const, 0.0);
    for(i = 0; i < num_const; i++) 
    {
      target_array[species_num][i] = floating_point_constant_list_array[i];
    }
  }
  /* else 
  { 
    declare_error("Invalid floating point list.");
    return(0);    
  } */
  return num_const;

}

int common_HC_processing( int type, int first_integer, int second_integer, int third_integer, 
      float first_float, float second_float, float third_float, float fourth_float,  float fifth_float)
{
  double range, range_0;
  
  if ( !accept_hcs_until_end_of_hc_card && (number_of_hcs_found == number_of_hcs_expected) ) return (0);
  if ( number_of_hcs_found == 0)
  {
    hunt = (struct HC_Information *) alloc_struct_1(struct HC_Information,1);
  }
  else
  {
    hunt = realloc(hunt, ( (number_of_hcs_found+1) * sizeof (struct HC_Information) ) );
  }
  iHC++;
  number_of_hcs_found++;
  hunt[iHC].nHC = iHC+1;
  hunt[iHC].Type = type;
  switch(type)
  {
    case 3: /* AC */
    case 1: /* BC */
      hunt[iHC].BCID = first_integer;
      hunt[iHC].DFID = second_integer;
      hunt[iHC].ramp = third_integer;
      hunt[iHC].BegParameterValue = first_float;
      hunt[iHC].EndParameterValue = second_float;
      hunt[iHC].Delta_s0 = third_float;
      hunt[iHC].Delta_s_min = fourth_float;
      hunt[iHC].Delta_s_max = fifth_float;
      break;
    case 2: /* MT */
      hunt[iHC].MTID = first_integer;
      hunt[iHC].MPID = second_integer;
      hunt[iHC].ramp = third_integer;
      hunt[iHC].BegParameterValue = first_float;
      hunt[iHC].EndParameterValue = second_float;
      hunt[iHC].Delta_s0 = third_float;
      hunt[iHC].Delta_s_min = fourth_float;
      hunt[iHC].Delta_s_max = fifth_float;
      hunt[iHC].MTID--;
      break;
  }   /* end of type switch */
  if (pd_glob[0]->Continuation == LOCA) /*This section is required for backward compatibility */
  {
    if ( iHC == 0)
    {
      cpcc = (struct Continuation_Conditions *) alloc_struct_1(struct Continuation_Conditions,1);
    }
    else
    {
      cpcc = realloc(cpcc, ( (iHC+1) * sizeof (struct Continuation_Conditions) ) );
    }    
    range_0 = hunt[0].EndParameterValue - hunt[0].BegParameterValue;
    if (hunt[0].ramp == 1)
    {
      hunt[0].Delta_s0 = range_0 / (double)(cont->MaxPathSteps - 1);
      loca_in->StepAggr = 0.0;
    }
    cpcc[iHC].nCC = hunt[iHC].nHC;
    cpcc[iHC].Type = hunt[iHC].Type;
    switch(hunt[iHC].Type) 
    {
      case 1:  /* BC */
        cpcc[iHC].BCID = hunt[iHC].BCID;
        cpcc[iHC].DFID = hunt[iHC].DFID;
        break;
      case 2:  /* MT */
        cpcc[iHC].MTID = hunt[iHC].MTID;
        cpcc[iHC].MPID = hunt[iHC].MPID;
         break;
    }
    cpcc[iHC].Beg_CC_Value = hunt[iHC].BegParameterValue;
    cpcc[iHC].End_CC_Value = hunt[iHC].EndParameterValue;
    range = hunt[iHC].EndParameterValue - hunt[iHC].BegParameterValue;
    if (range_0 != 0.0)
    {
      cpcc[iHC].ratio = range / range_0;
    }
    else
    {
      declare_error("Hunting condition card is leading to a zero range_0 value (End - Beg) & a divide by zero in ratio calculation.");
    }
    /* Fill in cont-> entries from hunt[0]. entries */
    hunt[0].Delta_s0 *= SGN(range_0);
    cont->upType = hunt[0].Type;
    cont->upBCID = hunt[0].BCID;
    cont->upDFID = hunt[0].DFID;
    cont->upMTID = hunt[0].MTID;
    cont->upMPID = hunt[0].MPID;
    cont->upMDID = hunt[0].MDID;
    cont->BegParameterValue = hunt[0].BegParameterValue;
    cont->EndParameterValue = hunt[0].EndParameterValue;
    cont->Delta_s0 = hunt[0].Delta_s0;
    cont->Delta_s_min = hunt[0].Delta_s_min;
    cont->Delta_s_max = hunt[0].Delta_s_max;	  
  }    
  return(number_of_hcs_found);
}

int common_error_element_size_processing ( float first_float, float second_float, float third_float, 
				float fourth_float,  float fifth_float, float sixth_float)
{
  if ( number_of_error_zz_elem_size_cards_found == 0)
  {
    pp_error_data = (struct Post_Processing_Error *) smalloc( sizeof(struct Post_Processing_Error) );
  }
  else
  {
    pp_error_data = realloc(pp_error_data, ( (number_of_error_zz_elem_size_cards_found+1) * sizeof (struct Post_Processing_Error) ) );
  }  
  pp_error_data[number_of_error_zz_elem_size_cards_found].error_params[0] = first_float;
  pp_error_data[number_of_error_zz_elem_size_cards_found].error_params[1] = second_float;
  pp_error_data[number_of_error_zz_elem_size_cards_found].error_params[2] = third_float;
  pp_error_data[number_of_error_zz_elem_size_cards_found].error_params[3] = fourth_float;
  pp_error_data[number_of_error_zz_elem_size_cards_found].error_params[4] = fifth_float;
  pp_error_data[number_of_error_zz_elem_size_cards_found].error_params[5] = sixth_float;
  number_of_error_zz_elem_size_cards_found++; 
  return(number_of_error_zz_elem_size_cards_found);
}

int common_pp_flux_processing ( char * flux_type_name, int flux_type, int ss_id, int blk_id, 
				int species_number,  char * flux_filenm, int profile_flag)				
{
  int i;
  if ( !accept_ppfs_until_end_of_ppf_card ) return (0);
  if ( nn_post_fluxes >= 20)
  {
    declare_error("Parser is limited to twenty post processing flux cards.");
    return (0);
  }
  if ( nn_post_fluxes == 0)
  {  
     pp_fluxes = (struct Post_Processing_Fluxes **) array_alloc(1, 20, sizeof(struct Post_Processing_Fluxes *));
     for(i = 0; i < 20; i++)      		
     {
       pp_fluxes[i] = (struct Post_Processing_Fluxes *) array_alloc(1, 1, sizeof(struct Post_Processing_Fluxes));
     }
  }
  strcpy(pp_fluxes[nn_post_fluxes]->flux_type_name, flux_type_name);
  pp_fluxes[nn_post_fluxes]->flux_type = flux_type;  
  pp_fluxes[nn_post_fluxes]->species_number = species_number;
  pp_fluxes[nn_post_fluxes]->ss_id = ss_id;
  pp_fluxes[nn_post_fluxes]->blk_id = blk_id;
  strcpy(pp_fluxes[nn_post_fluxes]->flux_filenm, flux_filenm);
  pp_fluxes[nn_post_fluxes]->profile_flag = profile_flag;

  nn_post_fluxes++;
  return(nn_post_fluxes);
}


void postprocess_EQ_cards ()
{
  int cv = 0;
  int ce, cem, cvm, n, neqn;
  PROBLEM_DESCRIPTION_STRUCT *pd_ptr = pd_glob[mn];
  

  /*
   *  Part 6 --
   * 
   *   Now figure out the mapping order, shape variables and
   *   the projection variable.
   *
   *   - pd_ptr->IntegrationMap variable is determined here
   *   - pd_ptr->ShapeVar       variable also
   *     pd_ptr->ProjectionVar
   */
  determine_ShapeVar(pd_ptr);
  determine_ProjectionVar(pd_ptr);
  
  /*
   *  Part 7 --
   * 
   *       In this section, we turn on the rest of the stress equations
   *       based on the number of modes specified and the base stress
   *       equation specified.
   */
   
  ce = R_STRESS11;
  cv = POLYMER_STRESS11;
  if ( pd_ptr->e[ce] )
  {
    cem = R_STRESS11_1;
    cvm = POLYMER_STRESS11_1;
    for(n=1; n<vn_glob[mn]->modes; n++)
    {
      /* set flag to turn on this stress mode */
      pd_ptr->e[cem]  = pd_ptr->e[ce];
      /* set weight function to be the same as mother mode,  R_STRESS11 */
      pd_ptr->w[cem] = pd_ptr->w[ce];
      /* make sure the variable associated with this equation is on */
      pd_ptr->v[cvm] = pd_ptr->v[cv];

      /* set interpolation function to be the same as mother mode,  R_STRESS11 */
      pd_ptr->i[cvm] = pd_ptr->i[cv];
      /* turn on the eqn term multipliers the same as mother's */
      pd_ptr->etm[cem][(LOG2_MASS)] = pd_ptr->etm[ce][(LOG2_MASS)];
      pd_ptr->etm[cem][(LOG2_ADVECTION)] = pd_ptr->etm[ce][(LOG2_ADVECTION)];
      pd_ptr->etm[cem][(LOG2_BOUNDARY)] = pd_ptr->etm[ce][(LOG2_BOUNDARY)];
      pd_ptr->etm[cem][(LOG2_DIFFUSION)] = pd_ptr->etm[ce][(LOG2_DIFFUSION)];
      pd_ptr->etm[cem][(LOG2_SOURCE)] = pd_ptr->etm[ce][(LOG2_SOURCE)];
	  
      pd_ptr->m[neqn] = cem;
	  
      /*
       * Setup actual equation to problem equation index arrays 
       */
      if( upd->ep[cem] == -1)
      {
	upd->ep[cem] = upd->Total_Num_EQ++;
      }

      if( upd->vp[cvm] == -1 )
      {
	upd->vp[cvm] = upd->Total_Num_Var++;
      }

      /* move on to the 11 component of the next stress tensor */
      cem += 6;
      cvm += 6;
      /* increment the number of equations to reflect the new stress modes */
      neqn++;
    }
  }

  ce = R_STRESS12;
  cv = POLYMER_STRESS12;
  if ( pd_ptr->e[ce] )
  {
    cem = R_STRESS12_1;
    cvm = POLYMER_STRESS12_1;
    for(n=1; n<vn_glob[mn]->modes; n++)
    {
      /* set flag to turn on this stress mode */
      pd_ptr->e[cem]  = pd_ptr->e[ce];
      /* set weight function to be the same as mother mode,  R_STRESS12 */
      pd_ptr->w[cem] = pd_ptr->w[ce];
      /* make sure the variable associated with this equation is on */
      pd_ptr->v[cvm] = pd_ptr->v[cv];

      /* set interpolation function to be the same as mother mode,  R_STRESS12 */
      pd_ptr->i[cvm] = pd_ptr->i[cv];
      /* turn on the eqn term multipliers the same as mother's */
      pd_ptr->etm[cem][(LOG2_MASS)] = pd_ptr->etm[ce][(LOG2_MASS)];
      pd_ptr->etm[cem][(LOG2_ADVECTION)] = pd_ptr->etm[ce][(LOG2_ADVECTION)];
      pd_ptr->etm[cem][(LOG2_BOUNDARY)] = pd_ptr->etm[ce][(LOG2_BOUNDARY)];
      pd_ptr->etm[cem][(LOG2_DIFFUSION)] = pd_ptr->etm[ce][(LOG2_DIFFUSION)];
      pd_ptr->etm[cem][(LOG2_SOURCE)] = pd_ptr->etm[ce][(LOG2_SOURCE)];
	  
      pd_ptr->m[neqn] = cem;
	  
      /*
       * Setup actual equation to problem equation index arrays 
       */
      if( upd->ep[cem] == -1)
      {
	upd->ep[cem] = upd->Total_Num_EQ++;
      }	 
      if( upd->vp[cvm] == -1 )
      {
	upd->vp[cvm] = upd->Total_Num_Var++;
      }
	  
      /* move on to the 12 component of the next stress tensor */
      cem += 6;
      cvm += 6;
      /* increment the number of equations to reflect the new stress modes */
      neqn++;
    }
  }

  ce = R_STRESS22;
  cv = POLYMER_STRESS22;
  if ( pd_ptr->e[ce] )
  {
    cem = R_STRESS22_1;
    cvm = POLYMER_STRESS22_1;
    for (n=1; n<vn_glob[mn]->modes; n++)
    {
      /* set flag to turn on this stress mode */
      pd_ptr->e[cem]  = pd_ptr->e[ce];
      /* set weight function to be the same as mother mode,  R_STRESS22 */
      pd_ptr->w[cem] = pd_ptr->w[ce];
      /* make sure the variable associated with this equation is on */
      pd_ptr->v[cvm] = pd_ptr->v[cv];

      /* set interpolation function to be the same as mother mode,
	 R_STRESS22 */
      pd_ptr->i[cvm] = pd_ptr->i[cv];
      /* turn on the eqn term multipliers the same as mother's */
      pd_ptr->etm[cem][(LOG2_MASS)] = pd_ptr->etm[ce][(LOG2_MASS)];
      pd_ptr->etm[cem][(LOG2_ADVECTION)] = pd_ptr->etm[ce][(LOG2_ADVECTION)];
      pd_ptr->etm[cem][(LOG2_BOUNDARY)] = pd_ptr->etm[ce][(LOG2_BOUNDARY)];
      pd_ptr->etm[cem][(LOG2_DIFFUSION)] = pd_ptr->etm[ce][(LOG2_DIFFUSION)];
      pd_ptr->etm[cem][(LOG2_SOURCE)] = pd_ptr->etm[ce][(LOG2_SOURCE)];
	  
      pd_ptr->m[neqn] = cem;
	  
      /*
       * Setup actual equation to problem equation index arrays 
       */
      if( upd->ep[cem] == -1)
      {
	upd->ep[cem] = upd->Total_Num_EQ++;
      }	 
      if( upd->vp[cvm] == -1 )
      {
	upd->vp[cvm] = upd->Total_Num_Var++;
      }
	  
      /* move on to the 22 component of the next stress tensor */
      cem += 6;
      cvm += 6;
      /* increment the number of equations to reflect the new stress modes */
      neqn++;
    }
  }

  ce = R_STRESS13;
  cv = POLYMER_STRESS13;
  if ( pd_ptr->e[ce] )
  {
    cem = R_STRESS13_1;
    cvm = POLYMER_STRESS13_1;
    for (n=1; n<vn_glob[mn]->modes; n++)
    {
      /* set flag to turn on this stress mode */
      pd_ptr->e[cem]  = pd_ptr->e[ce];
      /* set weight function to be the same as mother mode,  R_STRESS13 */
      pd_ptr->w[cem] = pd_ptr->w[ce];
      /* make sure the variable associated with this equation is on */
      pd_ptr->v[cvm] = pd_ptr->v[cv];

      /* set interpolation function to be the same as mother mode,  R_STRESS13 */
      pd_ptr->i[cvm] = pd_ptr->i[cv];
      /* turn on the eqn term multipliers the same as mother's */
      pd_ptr->etm[cem][(LOG2_MASS)] = pd_ptr->etm[ce][(LOG2_MASS)];
      pd_ptr->etm[cem][(LOG2_ADVECTION)] = pd_ptr->etm[ce][(LOG2_ADVECTION)];
      pd_ptr->etm[cem][(LOG2_BOUNDARY)] = pd_ptr->etm[ce][(LOG2_BOUNDARY)];
      pd_ptr->etm[cem][(LOG2_DIFFUSION)] = pd_ptr->etm[ce][(LOG2_DIFFUSION)];
      pd_ptr->etm[cem][(LOG2_SOURCE)] = pd_ptr->etm[ce][(LOG2_SOURCE)];
	  
      pd_ptr->m[neqn] = cem;
	  
      /*
       * Setup actual equation to problem equation index arrays 
       */
      if( upd->ep[cem] == -1)
      {
	upd->ep[cem] = upd->Total_Num_EQ++;
      }	 
      if( upd->vp[cvm] == -1 )
      {
	upd->vp[cvm] = upd->Total_Num_Var++;
      }
	  
      /* move on to the 13 component of the next stress tensor */
      cem += 6;
      cvm += 6;
      /* increment the number of equations to reflect the new stress modes */
      neqn++;
    }
  }

  ce = R_STRESS23;
  cv = POLYMER_STRESS23;
  if ( pd_ptr->e[ce] )
  {
    cem = R_STRESS23_1;
    cvm = POLYMER_STRESS23_1;
    for (n=1; n<vn_glob[mn]->modes; n++)
    {
      /* set flag to turn on this stress mode */
      pd_ptr->e[cem]  = pd_ptr->e[ce];
      /* set weight function to be the same as mother mode,  R_STRESS23 */
      pd_ptr->w[cem] = pd_ptr->w[ce];
      /* make sure the variable associated with this equation is on */
      pd_ptr->v[cvm] = pd_ptr->v[cv];

      /* set interpolation function to be the same as mother mode,  R_STRESS23 */
      pd_ptr->i[cvm] = pd_ptr->i[cv];
      /* turn on the eqn term multipliers the same as mother's */
      pd_ptr->etm[cem][(LOG2_MASS)] = pd_ptr->etm[ce][(LOG2_MASS)];
      pd_ptr->etm[cem][(LOG2_ADVECTION)] = pd_ptr->etm[ce][(LOG2_ADVECTION)];
      pd_ptr->etm[cem][(LOG2_BOUNDARY)] = pd_ptr->etm[ce][(LOG2_BOUNDARY)];
      pd_ptr->etm[cem][(LOG2_DIFFUSION)] = pd_ptr->etm[ce][(LOG2_DIFFUSION)];
      pd_ptr->etm[cem][(LOG2_SOURCE)] = pd_ptr->etm[ce][(LOG2_SOURCE)];
	  
      pd_ptr->m[neqn] = cem;
	  
      /*
       * Setup actual equation to problem equation index arrays 
       */
      if( upd->ep[cem] == -1)
      {
	upd->ep[cem] = upd->Total_Num_EQ++;
      }	 
      if( upd->vp[cvm] == -1 )
      {
	upd->vp[cvm] = upd->Total_Num_Var++;
      }
	  
      /* move on to the 23 component of the next stress tensor */
      cem += 6;
      cvm += 6;
      /* increment the number of equations to reflect the new stress modes */
      neqn++;
    }
  }

  ce = R_STRESS33;
  cv = POLYMER_STRESS33;
  if ( pd_ptr->e[ce] )
  {
    cem = R_STRESS33_1;
    cvm = POLYMER_STRESS33_1;
    for (n=1; n<vn_glob[mn]->modes; n++)
    {
      /* set flag to turn on this stress mode */
      pd_ptr->e[cem]  = pd_ptr->e[ce];
      /* set weight function to be the same as mother mode,  R_STRESS33 */
      pd_ptr->w[cem] = pd_ptr->w[ce];
      /* make sure the variable associated with this equation is on */
      pd_ptr->v[cvm] = pd_ptr->v[cv];

      /* set interpolation function to be the same as mother mode,  R_STRESS33 */
      pd_ptr->i[cvm] = pd_ptr->i[cv];
      /* turn on the eqn term multipliers the same as mother's */
      pd_ptr->etm[cem][(LOG2_MASS)] = pd_ptr->etm[ce][(LOG2_MASS)];
      pd_ptr->etm[cem][(LOG2_ADVECTION)] = pd_ptr->etm[ce][(LOG2_ADVECTION)];
      pd_ptr->etm[cem][(LOG2_BOUNDARY)] = pd_ptr->etm[ce][(LOG2_BOUNDARY)];
      pd_ptr->etm[cem][(LOG2_DIFFUSION)] = pd_ptr->etm[ce][(LOG2_DIFFUSION)];
      pd_ptr->etm[cem][(LOG2_SOURCE)] = pd_ptr->etm[ce][(LOG2_SOURCE)];
	  
      pd_ptr->m[neqn] = cem;
	  
      /*
       * Setup actual equation to problem equation index arrays 
       */
      if( upd->ep[cem] == -1)
      {
	upd->ep[cem] = upd->Total_Num_EQ++;
      }	
      if( upd->vp[cvm] == -1 )
      {
	upd->vp[cvm] = upd->Total_Num_Var++;
      }
 
      /* move on to the 33 component of the next stress tensor */
      cem += 6;
      cvm += 6;
      /* increment the number of equations to reflect the new stress modes */
      neqn++;
    }
  }

  /*
   *  Store the total number of equations that are active in this
   *   material in the Problem_Description structure
   */
  pd_ptr->Num_EQ = neqn;

  /*
   *  Check on tha bounds on the total number of equations permissible
   *  in a single material (NOTE: I don't know where the bounds comes
   *  from)
   */
  if ((neqn < 1 || neqn > 76)) {
    EH(-1,
       "rd_eq_specs: Too many (>76) or too few (<1) eqns. in this mat");
  }
  
  /* Consistency diagnostics for runs involving 3D stability of
   * 2D flows... */
  if(Linear_Stability == LSA_3D_OF_2D || Linear_Stability == LSA_3D_OF_2D_SAVE)
    {
      for(i = 0; i < upd->Num_Mat; i++)
	{
	  if(!pd_glob[i]->e[R_MOMENTUM3] || !pd_glob[i]->v[VELOCITY3])
	    {
	      fprintf(stderr, "\nR_MOMENTUM3/VELOCITY3 are required for a 3D stability of a 2D flow.\n");
	      fprintf(stderr, "They are missing in material %d (0-based)\n\n", i);
	      EH(-1, "missing equation for 3D stability of 2D flow");
	    }
	}
    }
} /* end postprocess_EQ_cards */


