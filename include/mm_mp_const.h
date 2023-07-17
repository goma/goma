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

#ifndef GOMA_MM_MP_CONST_H
#define GOMA_MM_MP_CONST_H

#ifndef MAX_NUMBER_MATLS
#define MAX_NUMBER_MATLS 9 /* maximum number of materials allowed */
#endif

#ifndef MAX_MODES
#define MAX_MODES 8 /* maximum number of viscoelastic modes allowed */
#endif

extern int Num_Var_Init_Mat[MAX_NUMBER_MATLS]; /* number of variables to overwrite with
                                          material-specific initialization */

#define DB_GOMA_MAT 1 /* Default location to obtain physical parameters */
#define DB_CHEMKIN_MAT                                      \
  4 /* Obtain physical parameters from the chemkin database \
     * for this material  */

/*
 *  Generic Model Types -> These can be used in all constituitive property
 *                         cases, unless otherwise overruled
 */
#define NO_MODEL -1 /* No model defined; good for initialization */
#define UNINITIALIZED_MODEL                            \
  0 /* Usually indicates that the model initialization \
     * process is not being handled with rigor     */
#define CHEMKIN_MODEL                                            \
  10005               /* Generic Chemkin model handled via calls \
                       * to the chemkin property database */
#define CONSTANT 1    /* Denotes CONSTANT property model */
#define RATIO    1001 /* Denotes RATIO property model for LS */
#define USER     2    /* Denotes USER-defined property model */
#define USER_GEN                                                    \
  3                     /* Denotes generalized User-defind property \
                             model with gradients */
#define TABLE        15 /* Denotes TABLE property model */
#define BILINEAR     6  /* Denotes option for table interpolation */
#define CAP_PRES     7
#define FAUX_PLASTIC 16
#define WAVELENGTH   17

/* Moment Property Models */
// growth rate
#define VISCOSITY_SCALED_GROWTH_RATE   4
#define VISCOSITY_PRESSURE_GROWTH_RATE 5

// coalescence
#define ADDITION_COALESCENCE               4
#define VISCOSITY_SCALED_COALESCENCE       5
#define VISCOSITY_BUBBLE_RATIO_COALESCENCE 6

/* Source term models */
#define BOUSS                                      \
  3 /* Boussinesq including rho*g*beta hydrostatic \
     * component */
#define VISC_DISS               4
#define JOULE                   5
#define SUSPENSION              6
#define BOUSSINESQ              7 /* Boussinesq as rho*g*beta*(T-Tref) */
#define BOUSS_JXB               8
#define SUSPEND                 9
#define FILL_SRC                10
#define VARIABLE_DENSITY        11   /* Drying of Polymeric Film */
#define LEVEL_SET               1212 /* density varies smoothly from - to + level set */
#define VE_LEVEL_SET            1213 /* density varies smoothly from - to + level set */
#define EHD_POLARIZATION        13   /* EHD Polarization force */
#define LS_QUADRATIC            14   /* variation on LEVEL_SET viscosity model */
#define ACOUSTIC                15   /* Acoustic energy density coupled to NS */
#define HS_FOAM                 16   /* fluorinert */
#define VISC_ACOUSTIC           17   /* heat generation by acoustics */
#define INGBER                  18
#define GRAV_VIBRATIONAL        19 /* momentum source for gravity + vibration */
#define MELT                    20 /* Lubrication source term model*/
#define EM_DISS                 21 /* heat generation by EM waves */
#define EM_VECTOR_DISS          22 /* heat generation by EM vector waves */
#define CONTINUUM_FLUID         23 /* Lubrication source term model*/
#define HS_FOAM_PBE             24
#define HS_FOAM_PMDI_10         25
#define VARIABLE_DENSITY_NO_GAS 26 /* Drying of Polymeric Film */

/* MMH */
/* #define  SUSPENSION_PM xxx Defiend below. */

/* Types of media */
/*#define  CONTINUOUS              3 */
/*#define  POROUS_SATURATED        4 */
/*#define  POROUS_UNSATURATED      5 */
/*#define  POROUS_TWO_PHASE        6 */
/*#define  POROUS_BRINKMAN         7 */

/* Types of media */
#define CONTINUOUS               0
#define POROUS_SATURATED         1
#define POROUS_UNSATURATED       2
#define POROUS_TWO_PHASE         3
#define POROUS_BRINKMAN          4
#define POROUS_SHELL_UNSATURATED 5

/* types of porosity relations */
#define DEFORM            3
#define EXTERNAL_FIELD    4
#define POROUS_CONST_INIT 5

/* types of drop patterns */
#define TFMP_SQUARE     300
#define TFMP_TRIANGULAR 301
#define TFMP_HEXAGONAL  302

// types of shell moment tensor calculation
#define SMT_SIMPLE   310
#define SMT_EXPANDED 311

// gap model types
// h = h0 - n dot d
#define GM_NDOTD 320
// h = r_web - r_roller
#define GM_RADIAL 321

// normal calculation methods
// use built-in fv->snormal
#define NCM_MAPPING 330
// for SIK_S normal of the web
#define NCM_PRIMITIVE_S_WEB 331
// for SIK_S normal of the roller
#define NCM_PRIMITIVE_S_ROLLER 332
// use built-in fv->normal
#define NCM_PRIMITIVE_XY 333

// shell integration kind
// integrate in 2d as in Scriven's ChEn 8104
// "Rudiments of Surface Geometry"
#define SIK_XY 340
// map 2d bar (shell) domain onto 1d element
#define SIK_S 341

/*
 * Options for k in potential equation
 * (electrical conductivity or permittivity)
 */
#define V_CONDUCTIVITY 0
#define V_PERMITTIVITY 1

// Electromagnetic
#define COMPLEX_CONSTANT 2
#define RADIAL_PML       3
#define REFRACTIVE_INDEX 4

/* Viscoelastic Constitutive equation parameters */

#define NOPOLYMER         0
#define GIESEKUS          6
#define WHITE_METZNER     7
#define OLDROYDB          8
#define PTT               9
#define SARAMITO_OLDROYDB 10
#define SARAMITO_GIESEKUS 11
#define SARAMITO_PTT      12
#define MODIFIED_JEFFREYS 13
#define ROLIE_POLY        14
#define ROLIE_POLY_FE     15
#define FENE_P            16
#define MODIFIED_WLF      39

// PTT forms
#define PTT_EXPONENTIAL 0
#define PTT_LINEAR      1

/**********************************************************************************/
/*  Density Models
 *
 *    This needs to be expanded
 *    However, all of the generic model types also apply here
 *      -1 NO_MODEL
 *       0 CHEMKIN_GENERIC
 *       1 CONSTANT
 *       2 USER
 *       3 USER_GEN
 */
#define DENSITY_FILL             29   /* same designation as FILL in rf_fem_const.h */
#define DENSITY_LEVEL_SET        1212 /* same designation as LEVEL_SET source model */
#define DENSITY_IDEAL_GAS        5    /* Ideal gas law. */
#define DENSITY_SUSPENSION       6    /* same designation as SUSPENSION IN mm_mp_const.h */
#define DENSITY_FOAM             7
#define DENSITY_FOAM_TIME        16
#define DENSITY_FOAM_TIME_TEMP   17
#define DENSITY_FOAM_CONC        18
#define DENSITY_CONST_PHASE_FUNC 15
#define DENSITY_CONSTANT_PMV                         \
  8 /* Density and concentration determined from     \
     * the assumption that the partial molar volumes \
     * are constant                                  \
     * -> appropriate for all species variable       \
     *    types, though aukward for mass fraction    \
     *    formulations.                              \
     */
#define DENSITY_CONSTANT_CONC                      \
  9 /* Density calculated from the assumption that \
     * the mixture concentration is constant       \
     * -> appropriate when the species variable    \
     *    type is mole fractions                   \
     */
#define DENSITY_CONSTANT_LAST_CONC 10
/* Density and mixture concentration calculated
 * assuming that the last species in the
 * phase has a constant concentration. All of the
 * other concentrations are allowed to freely
 * vary.
 * -> appropriate for dilute transport
 */
#define DENSITY_THERMEXP        13  /* Normal linear expansion for solids */
#define DENSITY_SUSPENSION_PM   14  /* special model for this system */
#define DENSITY_THERMAL_BATTERY 905 /* special density model for thermal bat work */
#define DENSITY_FOAM_PBE        34
#define DENSITY_FOAM_PBE_EQN    35
#define DENSITY_FOAM_PMDI_10    20
#define DENSITY_MOMENT_BASED    21
/**********************************************************************************/

/*
 *  Chemical Potential
 */
#define SSCHEMPOT_CONSTANT   1
#define SSCHEMPOT_POLYNOMIAL 5

#define PSCHEMPOT_PRESSURE_INDEPENDENT 1
#define PSCHEMPOT_IDEALGAS             5

#define CHEMPOT_IDEALSOLN   4
#define CHEMPOT_STOICHPHASE 5

/*
 * PBE Types
 */
#define PBE_R_11      0
#define PBE_N_PENTANE 1

/* #define FILL xxx already defined in rf_fem_const.h */

#define SOLVENT_POLYMER 100
#define REACTIVE_FOAM   101

/* MMH
 * I'm assuming that there is no chance of these overlapping with the
 * density model numbers.
 */
/* Generaralized Newtonian Constitutive equation parameters */

#define NEWTONIAN 3
#define POWER_LAW 4
#define CARREAU   5
/* #define SUSPENSION      6   key word defined in momentum source section   */
#define FILLED_EPOXY 7
#define CURE         8
#define THERMAL      9
/*
 *   HKM
 *         species source model, e.g. homogeneous reaction term
 *         this is also a heat source model and a viscosity model
 */
#define EPOXY                  10
#define SSM_EPOXY              10
#define BINGHAM                11
#define SSM_BINGHAM            11
#define EPOXY_DEA              12
#define SSM_EPOXY_DEA          12
#define CARREAU_SUSPENSION     13
#define SSM_CARREAU_SUSPENSION 13
#define SSM_BOND               14

#define SUSPENSION_PM                                                           \
  14                           /* Particle suspension, a la Yuri Buyevich       \
                                * suspension model.  It seemed best to define   \
                                * it here since this list is the longest.  This \
                                * isn't used for the actual c.r. b/c it is      \
                                * essentially Newtonian.  It is used for        \
                                * the source term, though.                      \
                                */
#define CARREAU_WLF         15 /*  Carreau viscosity with WLF temperature dependence */
#define POWERLAW_SUSPENSION 16

#define SSM_CHEMKIN_GAS 17 /* Chemkin Gas phase package call       */
#define SSM_CHEMKIN_LIQ 18 /* Chemkin Liquid phase package call    */
#define SSM_CHEMKIN_CPC 19 /* Chemkin Condensed phase package call */
#define FOAM            20 /* REF foam kinetics source model       */
#define CARREAU_WLF_CONC_PL                                            \
  21                        /*  Carreau viscosity with WLF temperature \
                               dependence and concentration shifting*/
#define HERSCHEL_BULKLEY 22 /* Herschel_bulkley model - power-law + yield stress */

#define BOND                 23 /* bond evolution structure model for viscosity */
#define CONST_PHASE_FUNCTION 24
#define CARREAU_WLF_CONC_EXP                    \
  25 /*  Carreau viscosity with WLF temperature \
        dependence and concentration shifting*/

#define FOAM_EPOXY     33
#define BINGHAM_WLF    27 /* Bingham WLF viscosity model */
#define SYLGARD        28 /* Sylgard viscosity model */
#define PRANDTL_MIXING 29 /* Shell Turbulent Viscosity Model */
#define BOND_SH        26 /* bond evolution structure model for viscosity with shear rate variable*/
#define BINGHAM_MIXED  30

#define FOAM_PBE_WATER 34
#define FOAM_PBE_OH    35
#define FOAM_PBE_CO2_L 36
#define FOAM_PBE_CO2_G 37
#define FOAM_PBE_BA_L  38
#define FOAM_PBE_BA_G  39

#define FOAM_PMDI_10         40
#define FOAM_PMDI_10_RXN     41
#define FOAM_PMDI_10_H2O     42
#define FOAM_PMDI_10_CO2     43
#define FOAM_PMDI_10_CO2_LIQ 44
#define FOAM_PMDI_10_CO2_GAS 45

#define MOMENT_CONSTANT_GROWTH            50
#define MOMENT_SIZE_DEPENDENT_COALESCENCE 51

/* Turbulent viscosity models for Reynolds Averaged NS */
#define TURBULENT_SA         52 /* Spallart Allmaras */
#define TURBULENT_SA_DYNAMIC 53 /* Spallart Allmaras */

/*
 *  Heat source modeling
 *
 */
#define HSM_EPOXY               10
#define HSM_BINGHAM             11
#define HSM_EPOXY_DEA           12
#define HSM_CARREAU_SUSPENSION  13
#define HSM_CARREAU_WLF         15 /*  Carreau viscosity with WLF temperature dependence */
#define HSM_POWERLAW_SUSPENSION 16
#define HSM_CK_GAS              17 /* Chemkin Gas phase package call       */
#define HSM_CK_LIQ              18 /* Chemkin Liquid phase package call    */
#define HSM_CK_CPC              19 /* Chemkin Condensed phase package call */
#define HSM_CARREAU_WLF_CONC    21
#define AVERAGE_CONTACT         22 /* Shell Energy Source Sliding Contact model */
#define LOCAL_CONTACT           23 /* Shell Energy Source Sliding Contact model */
#define LUBRICATION             24 /* Shell Energy Source Viscous Dissipation model */
#define CONSTANT_MELT           25 /* Shell Energy Source QCONV model */
#define CONSTANT_MELT_TURB      26 /* Shell Energy Source QCONV model */
#define LUBRICATION_FRICTION    27 /* Shell energy source VD model with solid-solid friction */
#define HSM_FOAM_PBE            34
/*
 *  Viscosity modeling
 *      -> HKM separated out the individual sections
 */
#define VISCM_EPOXY                10
#define VISCM_BINGHAM              11
#define VISCM_EPOXY_DEA            12
#define VISCM_CARREAU_SUSPENSION   13
#define VISCM_CARREAU_WLF          15 /*  Carreau viscosity with WLF temperature dependence */
#define VISCM_POWERLAW_SUSPENSION  16
#define VISCM_CK_GASMIXTAVG        17 /* Chemkin Gas phase mixture averaged      */
#define VISCM_CK_GASDIXONLEWIS     18 /* Chemkin Gas phase Multicomponent Dixon-Lewis  */
#define VISCM_CK_LIQ               19 /* Chemkin Liquid phase package call    */
#define VISCM_CK_CPC               20 /* Chemkin Condensed phase package call */
#define VISCM_CARREAU_WLF_CONC     21
#define VISCM_BOND                 23 /* bond evolution structure model for viscosity */
#define VISCM_CONST_PHASE_FUNCTION 24
#define VISCM_CARREAU_WLF_CONC_EXP              \
  25 /*  Carreau viscosity with WLF temperature \
  dependence and concentration shifting*/
#define VISCM_BOND_SH                                                \
  26 /* bond evolution structure model for viscosity with shear rate \
        variable*/

#define VISCM_FOAM_EPOXY     33
#define VISCM_BINGHAM_WLF    27 /* Bingham WLF viscosity model */
#define VISCM_SYLGARD        28 /* Sylgard viscosity model */
#define VISCM_PRANDTL_MIXING 29 /* Shell Turbulent Viscosity Model */
/*
 * Dilational Viscosity Model
 *
 */
//!               Assume that kappa = 2/3 mu, and term magically disappears (default)
#define DILVISCM_KAPPAWIPESMU 10
//!               Assume that kappa = constant
#define DILVISCM_KAPPACONSTANT 11
//!               Assume that kappa = R mu, where R is a constant, and borrow mu constitutive model
#define DILVISCM_KAPPAFIXEDRATIO 12
//!               Assume that kappa comes from bubble theory
#define DILVISCM_KAPPABUBBLES 13

/* types of diffusion coefficient relations */
#define POROUS               3
#define HYDRO                4
#define FREE_VOL             5
#define ANISOTROPIC          6
#define BISECTION            7
#define RZBISECTION          8
#define RICHARDSON_ZAKI      9
#define EXP_DECAY            10
#define GENERALIZED          11 /*used in conjunc. with Generalized Fickian */
#define GENERALIZED_FREE_VOL 12
#define SUSP_BAL             13
#define ARRHENIUS            14 /* for temperature-dependent S-M diffusivities, KSC */
#define SHOCK                15
#define PIECEWISE            16
#define CHAPMAN_GAS          17

/* Types of vapor or gas pressure relations */
#define KELVIN       3
#define FLAT         33
#define IDEAL_GAS    4
#define NON_VOLATILE 5
#define ANTOINE      6
#define RIEDEL       7

/* Types of saturation and permeability relations */
#define VAN_GENUCHTEN          3
#define SUM_TO_ONE             4
#define PSD_VOL                5
#define PSD_WEXP               6
#define PSD_SEXP               7
#define K_TENSOR               8
#define SOLIDIFICATION         9 /* permeability that slows down velocity for phase change  */
#define TANH                   10
#define TANH_HYST              11
#define KOZENY_CARMAN          12
#define SINK_MASS_PERM         13
#define ORTHOTROPIC            14
#define KC_TENSOR              15
#define SM_TENSOR              16
#define SHELL_CYLINDER_SQUARE  20
#define SHELL_TANH             21
#define TANH_EXTERNAL          22
#define VAN_GENUCHTEN_EXTERNAL 23
#define LEVER                  24
#define SATURATION             25
#define ATANH                  26
#define SINH                   27
#define VAN_GENUCHTEN_HYST     28
#define VAN_GENUCHTEN_HYST_EXT 29

/* Types of Flowing Liquid Viscosity Models */
#define MOLTEN_GLASS 3

/* Types of capillary pressure stress */
#define NO_CAP_STRESS     2
#define WETTING           3
#define PARTIALLY_WETTING 4
#define COMPRESSIBLE      5

/* Elastic Constitutive equation parameters */
#define LINEAR          3
#define NONLINEAR       4
#define INCOMP_3D       5
#define INCOMP_PSTRAIN  6
#define INCOMP_PSTRESS  7
#define HOOKEAN_PSTRAIN 8
#define HOOKEAN_PSTRESS 9

/* Viscoplastic consitutive equation params */
#define EVP_HYPER 10

/* Viscoelastic consitutive equation params */
#define KELVIN_VOIGT 41

/* Modulus parameters */
/*#define POWER_LAW    4  - defined rf_fem_const.h*/
#define CONTACT_LINE    5
#define SHEAR_HARDEN    6
#define EXPONENTIAL     7
#define DENSE_POWER_LAW 8
#define POISSON_RATIO   9
#define SHRINKAGE       10

/* Diffusion Constitutive equation parameters */
#define FICKIAN                3
#define DARCY                  4
#define DARCY_FICKIAN          5
#define NON_DIFFUSING          6
#define HYDRODYNAMIC           7
#define STEFAN_MAXWELL         8  /* Stefan-Maxwell diffusion of neutral species, KSC 7/98 */
#define STEFAN_MAXWELL_CHARGED 9  /* Stefan-Maxwell diffusion of charged species, KSC 9/98 */
#define DM_CK_GASMIXTUREAVG    10 /* Chemkin MixtureAveraged diffusivities */
#define DM_CK_GASMIXTUREAVG_VC                \
  11 /* Chemkin MixtureAveraged diffusivities \
      * with velocity correction thrown in */
#define DM_CK_GASDIXONLEWIS                \
  12 /* Chemkin Dixon-Lewis multicomponent \
      * diffusivities with velocity corrections */
#define DM_CK_LQSTEFAN_MAXWELL                   \
  13 /* Same as regular Stefan-Maxwell, but with \
      * coefficients coming from liquid chemkin  \
      * package */
#define DM_CK_LQSTEFAN_MAXWELL_CHARGED 14
/*
 * Same as regular Stefan-Maxwell_charged, but with
 * coefficients coming from liquid chemkin
 * package
 */
#define GENERALIZED_FICKIAN      15 /* generalized fickian ACS 4/00 */
#define STEFAN_MAXWELL_VOLUME    16 /* RSL 6/28/00 */
#define FICKIAN_CHARGED          17 /* Fickian diffusion of charged species, KSC 9/2000 */
#define DM_SUSPENSION_BALANCE    18 /* for Nott and Brady type models */
#define FICKIAN_CHARGED_X        19 /*  RSL 9/18/00  */
#define POWERLAW_DARCY_FICKIAN   20 /*PRS for P&G */
#define HYDRODYNAMIC_QTENSOR     21
#define HYDRODYNAMIC_QTENSOR_OLD 22
#define FICKIAN_SHELL            23 /* Shell version of Fickian diffusion equation */
/* surface tension laws */
#define DILATION       3
#define GIBBS_ISOTHERM 35

/* species only diffusion choices */
#define DIFF_OFF      0
#define DIFF_POSITIVE 10
#define DIFF_NEGATIVE 11

/* Species Time Integration choices */
#define STANDARD            0
#define TAYLOR_GALERKIN     1
#define TAYLOR_GALERKIN_EXP 2

/*Convective langrangian velocity models */
#define ROTATIONAL    2
#define ROTATIONAL_3D 25
#define OSC_LINEAR    252

/*Various thermophysical property models */
#define ENTHALPY     4
#define THERMAL_HEAT 45
#define FOAM_PBE     46

/*Electrode-kinetics Species Source model: KSC 10/13/98 */
#define ELECTRODE_KINETICS 904
#define ION_REACTIONS      908 /* RSL 3/13/01 */

/*Thermal-battery property model: KSC 3/1/99 */
#define THERMAL_BATTERY 905

/* Electrolyte-conductvity and current-density models for Charge-Species Transport: KSC 9/2000 */
#define ELECTRONEUTRALITY_SM      906
#define ELECTRONEUTRALITY_FICKIAN 907

/*Thermodynamic potential models: KSC 2/21/02 */
#define FeS2 908
#define LiSi 909

/* Atmospheric metal corrosion kinetic models: KSC 3/27/02 */
#define SOLID_DIFFUSION_SIMPLIFIED               910
#define SOLID_DIFFUSION_ELECTRONEUTRALITY        911
#define SOLID_DIFFUSION                          912
#define GAS_DIFFUSION                            913
#define FULL                                     914
#define ANNIHILATION_ELECTRONEUTRALITY           915
#define ANNIHILATION                             916
#define NET_CHARGE                               917 /* refer to F multiplied by sum of ci zi   */
#define SOLID_DIFFUSION_ELECTRONEUTRALITY_LINEAR 918
#define DEBYE_HUCKEL                             919 /* Debye-Huckel linearization approximation for potential */

/* Butler-Volmer kinetic and electro-osmotic transport models: KSC 4/20/2006 */
#define BUTLER_VOLMER  920
#define ELECTROOSMOTIC 921

/* Photocuring reaction model  */
#define PHOTO_CURING 934

/* KOH etching models */
#define ETCHING_KOH     935
#define ETCHING_KOH_EXT 936

/* Special material-related height function models */
#define CONSTANT_SPEED        1011
#define ROLL_ON               1012
#define ROLL                  1013
#define ROLLER                101301
#define LINEAR_TIME           10130
#define CONSTANT_SPEED_DEFORM 10131
#define CONSTANT_SPEED_MELT   10132
#define FLAT_GRAD_FLAT        10133
#define TANGENTIAL_ROTATE     10134
#define ROLL_ON_MELT          10135
#define POLY_TIME             10136
#define JOURNAL               10137
#define CAP_SQUEEZE           10140
#define SLIDER_POLY_TIME      10141
#define FLAT_GRAD_FLAT_MELT   10142
#define CIRCLE_MELT           10143
#define LOWER_DISTANCE        10144

/* Lubrication contact angle models */
#define DYNAMIC_CA        10201
#define DYNAMIC_LINEAR_CA 10202

/* Disjoining pressure model */
#define TWO_TERM        10150
#define TWO_TERM_EXT_CA 10152
#define ONE_TERM        10151

/* Film evaporation model */
#define CONC_POWER 10160

/* Diffusion coefficient model */
#define STOKES_EINSTEIN 10170

/* Special function models for structured porous shells */
#define MULTI_MODE 1014

/*

   CONSTANTS FOR MATERIAL PROPERTY TAGS

   NOT COMPLETE

   IAN GATES

*/

/*
 * General Model Constants
 */

#define TAGC_THERMAL_CONDUCTIVITY    1100
#define TAGC_ELECTRICAL_CONDUCTIVITY 1200
#define TAGC_PERMITTIVITY            1210
#define TAGC_VISCOSITY               1300
#define TAGC_SURFACE_TENSION         1400
#define TAGC_HEAT_CAPACITY           1500
#define TAGC_VOLUME_EXPANSION        1600
#define TAGC_DENSITY                 1700
#define TAGC_POROSITY                1800
#define TAGC_PERMEABILITY            1900
#define TAGC_REL_GAS_PERM            2000
#define TAGC_REL_LIQ_PERM            2100
#define TAGC_SATURATION              2200
#define TAGC_POROUS_COMPRESSIBILITY  2300
#define TAGC_MELTING_POINT_LIQUIDUS  2500
#define TAGC_MELTING_POINT_SOLIDUS   2600
#define TAGC_FLOWINGLIQUID_VISCOSITY 2700
#define TAGC_DIFFUSIVITY_0           2800
#define TAGC_DIFFUSIVITY_1           2801

/*
 *  Acoustic Model Constants
 */

#define TAGC_ACOUSTIC_WAVENUMBER 3000
#define TAGC_ACOUSTIC_IMPEDANCE  3010
#define TAGC_ACOUSTIC_ABSORPTION 3020
#define TAGC_REFRACTIVE_INDEX    3030
#define TAGC_LIGHT_ABSORPTION    3040
#define TAGC_EXTINCTION_INDEX    3050

/*
 * Generalized Newtonian Models:
 * Newtonian, Power Law, Carreau or Bingham(1,2,3)
 */

#define TAGC_MU0   4000
#define TAGC_NEXP  4100
#define TAGC_MUINF 4200
#define TAGC_LAM   4300
#define TAGC_AEXP  4400
#define TAGC_ATEXP 4500

/*  CARREAU_WLF  */

#define TAGC_WLFC2   4550
#define TAGC_REFTEMP 4560

/* these are for the BINGHAM yielding material model */

#define TAGC_TAU_Y 4600
#define TAGC_FEXP  4610

/* these are for SUSPENSION/FILLED_EPOXY models */

#define TAGC_MAXPACK    4800
#define TAGC_FICKDIFF_X 4810
#define TAGC_FICKDIFF_Y 4820

/* these are for Ryan's Qtensor model: */
#define TAGC_QTENSOR_EXTENSION_P 4830
#define TAGC_QTENSOR_NCT         4840

/* these are for CURE/EPOXY/FILLED_EPOXY models */

#define TAGC_GELPOINT 4900
#define TAGC_CUREAEXP 4910
#define TAGC_CUREBEXP 4920

/*
 * Viscoelastic Constitutive Equation: Giesekus Model
 * with proper coefficient choices it can become:
 * Maxwell Model
 * Oldroyd-B Model
 * White-Metzner Model
 * Leonov Model
 *
 * this struct contains the polymer viscosity
 * if it is shearthinning etc or NEWTONIAN
 */

/*
 *  this is the ID for the first ve mode
 *  subsequent mode IDs are incremented by 1
 *  TAGC_TIME_CONST(mode i) = TAGC_TIME_CONST +i
 */

#define TAGC_TIME_CONST             5000
#define TAGC_TIME_CONST1            5001
#define TAGC_TIME_CONST2            5002
#define TAGC_TIME_CONST3            5003
#define TAGC_TIME_CONST4            5004
#define TAGC_TIME_CONST5            5005
#define TAGC_TIME_CONST6            5006
#define TAGC_TIME_CONST7            5007
#define TAGC_WT_FUNC                5100
#define TAGC_ALPHA                  5200
#define TAGC_ALPHA1                 5201
#define TAGC_ALPHA2                 5202
#define TAGC_ALPHA3                 5203
#define TAGC_ALPHA4                 5204
#define TAGC_ALPHA5                 5205
#define TAGC_ALPHA6                 5206
#define TAGC_ALPHA7                 5207
#define TAGC_PTT_XI                 5300
#define TAGC_PTT_EPS                5400
#define TAGC_SHIFT_FUNC             5500
#define TAGC_SHIFT_FUNC1            5501
#define TAGC_POLYMER_YIELD_STRESS   5600
#define TAGC_POLYMER_YIELD_EXPONENT 5700

/*
 * Constants used in the Elasticity Constitutive Equations
 */

#define TAGC_LAME_MU                 6000
#define TAGC_LAME_MU_CONTACT_LINE_G0 6001
#define TAGC_LAME_MU_CONTACT_LINE_G1 6002
#define TAGC_LAME_MU_CONTACT_LINE_R0 6003
#define TAGC_LAME_LAMBDA             6100
#define TAGC_BEND_STIFFNESS          6110
#define TAGC_CONV_LAG_VELX           6201
#define TAGC_CONV_LAG_VELY           6202
#define TAGC_CONV_LAG_VELZ           6203
#define TAGC_CONV_LAG_ROTRATE        6221
#define TAGC_CONV_LAG_ROT_X0         6222
#define TAGC_CONV_LAG_ROT_Y0         6223
#define TAGC_CONV_LAG_ROT_Z0         6224

#define TAGC_RS_LAME_MU          6300
#define TAGC_RS_LAME_LAMBDA      6400
#define TAGC_RS_CONV_LAG_VELX    6301
#define TAGC_RS_CONV_LAG_VELY    6302
#define TAGC_RS_CONV_LAG_VELZ    6303
#define TAGC_RS_CONV_LAG_ROTRATE 6321
#define TAGC_RS_CONV_LAG_ROT_X0  6322
#define TAGC_RS_CONV_LAG_ROT_Y0  6323
#define TAGC_RS_CONV_LAG_ROT_Z0  6324

#define TAGC_POISSON               6600
#define TAGC_STRSS_FR_SOL_VOL_FRAC 6610

/*
 * Constants used for Source Term Models
 */

/*
 * General Navier-Stokes Source requires
 * three constants for vector components;
 * BOUSS_JXB and SUSPEND require 4th constant
 */

#define TAGC_NSS_A0 7000
#define TAGC_NSS_A1 7001
#define TAGC_NSS_A2 7002
#define TAGC_NSS_A3 7003

/*
 * Lubrication Constants:
 * heightU, heightL, veloU, veloL, dcaU, dcaL
 */

#define TAGC_LUB_HGT_U0 7010
#define TAGC_LUB_HGT_U1 7011
#define TAGC_LUB_HGT_U2 7012
#define TAGC_LUB_HGT_U3 7013
#define TAGC_LUB_HGT_U4 7014
#define TAGC_LUB_HGT_U5 7015
#define TAGC_LUB_HGT_U6 7016
#define TAGC_LUB_HGT_U7 7017

#define TAGC_LUB_HGT_L0 7018
#define TAGC_LUB_HGT_L1 7019
#define TAGC_LUB_HGT_L2 7020
#define TAGC_LUB_HGT_L3 7021
#define TAGC_LUB_HGT_L4 7022
#define TAGC_LUB_HGT_L5 7023
#define TAGC_LUB_HGT_L6 7024
#define TAGC_LUB_HGT_L7 7025

#define TAGC_U_LUB_VELO_U0 7026
#define TAGC_U_LUB_VELO_U1 7027
#define TAGC_U_LUB_VELO_U2 7028
#define TAGC_U_LUB_VELO_U3 7029
#define TAGC_U_LUB_VELO_U4 7030
#define TAGC_U_LUB_VELO_U5 7031

#define TAGC_U_LUB_VELO_L0 7032
#define TAGC_U_LUB_VELO_L1 7033
#define TAGC_U_LUB_VELO_L2 7034
#define TAGC_U_LUB_VELO_L3 7035
#define TAGC_U_LUB_VELO_L4 7036
#define TAGC_U_LUB_VELO_L5 7037

#define TAGC_LUB_VELO_U0 17026
#define TAGC_LUB_VELO_U1 17027
#define TAGC_LUB_VELO_U2 17028

#define TAGC_LUB_VELO_L0 17032
#define TAGC_LUB_VELO_L1 17033
#define TAGC_LUB_VELO_L2 17034

#define TAGC_LUB_DCA_U0 7038
#define TAGC_LUB_DCA_U1 7039
#define TAGC_LUB_DCA_U2 7040
#define TAGC_LUB_DCA_U3 7041
#define TAGC_LUB_DCA_L0 7042
#define TAGC_LUB_DCA_L1 7043
#define TAGC_LUB_DCA_L2 7044
#define TAGC_LUB_DCA_L3 7045

#define TAGC_LUB_SOURCE_0 7046
#define TAGC_LUB_SOURCE_1 7047
#define TAGC_LUB_SOURCE_2 7048

#define TAGC_HEAT_SOURCE_0       7050
#define TAGC_SPECIES_SOURCE_0_P0 70520
#define TAGC_SPECIES_SOURCE_0_P1 70521
#define TAGC_SPECIES_SOURCE_0_P2 70522
#define TAGC_SPECIES_SOURCE_0_P3 70523
#define TAGC_SPECIES_SOURCE_1_P0 70530
#define TAGC_SPECIES_SOURCE_1_P1 70531
#define TAGC_SPECIES_SOURCE_1_P2 70532
#define TAGC_SPECIES_SOURCE_1_P3 70533
#define TAGC_RST_FUNC_0          7060
#define TAGC_RST_FUNC_1          7061
#define TAGC_RST_FUNC_2          7062
#define TAGC_LATENT_HEAT_0       7070
#define TAGC_LATENT_HEAT_1       7071

/*  Problem Description Parameters   */

#define TAGC_ACOUSTIC_FREQ       8010
#define TAGC_PROCESS_TEMP        8011
#define TAGC_ACOUSTIC_WAVELENGTH 8012
#define TAGC_EM_FREQ             8013
#define TAGC_EM_WAVELENGTH       8014
/*
 * Thin film multiphase constants
 * */

#define TAGC_TFMP_REL_PERM_0  7100
#define TAGC_TFMP_REL_PERM_1  7101
#define TAGC_TFMP_REL_PERM_2  7102
#define TAGC_TFMP_REL_PERM_3  7103
#define TAGC_TFMP_DENSITY_0   7104
#define TAGC_TFMP_DENSITY_1   7105
#define TAGC_TFMP_DENSITY_2   7106
#define TAGC_TFMP_DENSITY_3   7107
#define TAGC_TFMP_VISCOSITY_0 7108
#define TAGC_TFMP_VISCOSITY_1 7109

#endif
