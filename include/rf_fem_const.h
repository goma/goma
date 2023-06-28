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
 *$Id: rf_fem_const.h,v 5.11 2010-05-17 20:33:27 sarober Exp $
 */

/* rf_fem_const.h
 *
 *	    Definition of globally occurring constants relating
 *	    to the details of the FEM problem.
 *
 *	    Author:
 *	    Date:           11/13/92
 *	    Revised:        Sat Mar 19 13:11:54 MST 1994 pasacki@sandia.gov
 *	    Revised:        1997/02/27 10:19 MST pasacki@sandia.gov
 */

/*
 * Use this symbol to conditionally include this include file if it hasn't
 * already been included...
 */

#ifndef GOMA_RF_FEM_CONST_H
#define GOMA_RF_FEM_CONST_H

/* Generic Logicals */
#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

/* Geometric parameters */
#define CARTESIAN         1
#define CYLINDRICAL       2
#define SPHERICAL         3
#define SWIRLING          4 /* axisymmetric 2 1/2 D */
#define POLAR             5
#define ELLIPTIC_CYLINDER 6
#define PROJECTED_CARTESIAN               \
  7 /* Used for 3D stability of 2D flows, \
     * for now.  This is 2D mesh, but 3D  \
     * velocities, in normal Cartesian    \
     * coordinates.  Lucky number 7! */
#define CARTESIAN_2pt5D 8

/* Integration mapping parameters */
#define ISOPARAMETRIC   1
#define SUBPARAMETRIC   2
#define SUPERPARAMETRIC 3

/* shell tangent computation selection is ISOPARAMETRIC or SEEDED */
#define SEEDED 2

/* Viscoelastic Constitutive equation weight functions */
#define GALERKIN    1
#define SUPG        2
#define SUPG_GP     3
#define SUPG_SHAKIB 4

#define SC_NONE       0
#define YZBETA_ONE    1
#define YZBETA_TWO    2
#define YZBETA_MIXED  3
#define YZBETA_CUSTOM 4
#define SC_DCDD       5
#define SC_YZBETA     6

/* Viscoelastic Constitutive equation formulation */
#define EVSS_G     1               /* Rajagopalan's formulation */
#define EVSS_F     2               /* Fortin's formulation */
#define EVSS_L     3               /* Level set solid-fluid formulation */
#define LOG_CONF   4               /* Log-conformation tensor formulation */
#define EVSS_GRADV 5               /* Fortin's formulation using GradV instead of G */
                                   /* for stress constitutive equations */
#define LOG_CONF_GRADV 6           /* Log-conformation tensor formulation using */
                                   /* GradV instead of G for stress constitutive equations */
#define LOG_CONF_TRANSIENT 7       /* Log-conformation tensor formulation using */
                                   /* lagged (explicit) terms for eigen-decomp parts */
#define LOG_CONF_TRANSIENT_GRADV 8 /* Log-conformation tensor formulation using */
                                   /* lagged (explicit) terms for eigen-decomp parts grad(v) form*/
#define SQRT_CONF 9                /* Log-conformation tensor formulation using */
                                   /* lagged (explicit) terms for eigen-decomp parts grad(v) form*/
#define CONF 10                    /* Log-conformation tensor formulation using */
/* Discontinuous Galerkin viscoelastic jacobian options */
#define EXPLICIT_DG 1
#define FULL_DG     2
#define SEGREGATED  3
#define LUMPED      4

/* Mesh Motion parameters */
#define ARBITRARY          1
#define LAGRANGIAN         2
#define TOTAL_ALE          3
#define DYNAMIC_LAGRANGIAN 4
#define DYNAMIC_TOTAL_ALE  5

/* Element quality metric parameters */
#define EQM_NONE  0
#define EQM_JAC   1
#define EQM_VOL   2
#define EQM_ANG   3
#define EQM_TRI   4
#define EQM_AVG   8
#define EQM_MIN   9
#define EQM_WTMIN 10

/* Fill Weight Function options */
#define FILL_WEIGHT_G           1
#define FILL_WEIGHT_TG          2
#define FILL_WEIGHT_SUPG        3
#define FILL_WEIGHT_GLS         4
#define FILL_WEIGHT_SC          5
#define FILL_WEIGHT_EXPLICIT    6
#define FILL_WEIGHT_SUPG_SHAKIB 7
#define FILL_WEIGHT_SUPG_GP     8

/* Fill Equation options */
#define FILL_EQN_ADVECT  1
#define FILL_EQN_EXT_V   2
#define FILL_EQN_EIKONAL 3

/* Fluid-structural interactions for shells */
#define FSI_MESH_BOTH           1
#define FSI_MESH_CONTINUUM      2
#define FSI_REALSOLID_CONTINUUM 3
#define FSI_MESH_SHELL          4
#define FSI_SHELL_ONLY          5
#define FSI_MESH_UNDEF          6
#define FSI_MESH_ONEWAY         7
#define FSI_SHELL_ONLY_MESH     8
#define FSI_SHELL_ONLY_UNDEF    9

/* Lubrication Viscosity Integration Options */
#define LUB_VISCINT_GAUSSIAN   25
#define LUB_VISCINT_ANALYTICAL 26
#define LUB_VISCINT_POWERLAW   27
#define MAX_LUB_NGP            5

/* Residence time kernel functions */
#define LINEAR_TIMETEMP      1110
#define EXPONENTIAL_TIMETEMP 1120
#define DROP_EVAP            1130

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/
/*
 * Variable Names of unknowns...
 *
 * Note: assignments are arbitrary except that:
 * 	(1) velocity1...velocity3 must be in sequence.
 *          i.e., VELOCITY1+1=VELOCITY2, etc.
 *          Ditto mesh displacement.
 *	(2) Must increase by 1 without gaps from lowest to highest
 *	    variable.
 */

#define VELOCITY1   0 /* aka "U1" */
#define VELOCITY2   1 /* aka "U2" */
#define VELOCITY3   2 /* aka "U3" */
#define TEMPERATURE 3 /* aka "T" */
#define MASS_FRACTION                                                               \
  4                          /* aka "Y" actually volume fraction (or gas and liquid \
                              *  pressures or liquid vol.frac. in porous media      \
                              *  -> Basically, this is a catchall for the species   \
                              *  unknowns (SHOULD ACTUALLY BE RENAMED,              \
                              *  SPECIES_UNKNOWNS!)                                 \
                              *                                                     \
                              * HKM NOTE-> This variable type will be going away,   \
                              *     to be replaced by the SPECIES_UNK_# variable    \
                              *     types described below.                          \
                              */
#define MESH_DISPLACEMENT1 5 /* aka "D1" */
#define MESH_DISPLACEMENT2 6 /* aka "D2" */
#define MESH_DISPLACEMENT3 7 /* aka "D3" */
#define SURFACE                                                         \
  8                /* aka "S"                                           \
                    *  HKM NOTE -> I don't see an application for this  \
                    *  variable type in its current form.               \
                    *  Therefore, it will either be going away or it    \
                    *  will be unrolled wrt surface species subvariable \
                    *  types (aka like the SPECIES_UNK_# variables)     \
                    */
#define PRESSURE 9 /* aka "P" */

/*
 *  Define names for tensor fields
 */
#define POLYMER_STRESS11 10 /* aka S11, 11 component of polymer stress tensor */
#define POLYMER_STRESS12 11 /* aka S12, 12 component of polymer stress tensor */
#define POLYMER_STRESS22 12 /* aka S22, 22 component of polymer stress tensor */
#define POLYMER_STRESS13 13 /* aka S13, 13 component of polymer stress tensor */
#define POLYMER_STRESS23 14 /* aka S23, 23 component of polymer stress tensor */
#define POLYMER_STRESS33 15 /* aka S33, 33 component of polymer stress tensor */

#define SOLID_DISPLACEMENT1 16 /* aka "d1'"  or real-solid displacement*/
#define SOLID_DISPLACEMENT2 17 /* aka "d2'"  or real-solid displacement*/
#define SOLID_DISPLACEMENT3 18 /* aka "d3'"  or real-solid displacement*/

#define VELOCITY_GRADIENT11 19 /* aka G11, 11 of the velocity gradient tensor */
#define VELOCITY_GRADIENT12 20 /* aka G12, 12 of the velocity gradient tensor */
#define VELOCITY_GRADIENT21 21 /* aka G21, 21 of the velocity gradient tensor */
#define VELOCITY_GRADIENT22 22 /* aka G22, 22 of the velocity gradient tensor */
#define VELOCITY_GRADIENT13 23 /* aka G13, 13 of the velocity gradient tensor */
#define VELOCITY_GRADIENT23 24 /* aka G23, 23 of the velocity gradient tensor */
#define VELOCITY_GRADIENT31 25 /* aka G31, 31 of the velocity gradient tensor */
#define VELOCITY_GRADIENT32 26 /* aka G32, 32 of the velocity gradient tensor */
#define VELOCITY_GRADIENT33 27 /* aka G33, 33 of the velocity gradient tensor */

#define VOLTAGE    28   /* aka "V" */
#define FILL       29   /* aka "F" */
#define LS         FILL /* That is FILL and LS share the same index */
#define SHEAR_RATE 30   /* aka "SH", shear rate from second invariant of ros tensor */

#define PVELOCITY1 31 /* aka "P_U", particle momentum */
#define PVELOCITY2 32 /* aka "P_V", particle momentum */
#define PVELOCITY3 33 /* aka "P_W", particle momentum */

/* Please leave these stress variables continuously numbered or
 * multimode viscoelasticity will never work!
 */
#define POLYMER_STRESS11_1 34 /* aka S11, 11 component of polymer stress tensor, mode 2 */
#define POLYMER_STRESS12_1 35 /* aka S12, 12 component of polymer stress tensor, mode 2 */
#define POLYMER_STRESS22_1 36 /* aka S22, 22 component of polymer stress tensor, mode 2 */
#define POLYMER_STRESS13_1 37 /* aka S13, 13 component of polymer stress tensor, mode 2 */
#define POLYMER_STRESS23_1 38 /* aka S23, 23 component of polymer stress tensor, mode 2 */
#define POLYMER_STRESS33_1 39 /* aka S33, 33 component of polymer stress tensor, mode 2 */

#define POLYMER_STRESS11_2 40 /* aka S11, 11 component of polymer stress tensor, mode 3 */
#define POLYMER_STRESS12_2 41 /* aka S12, 12 component of polymer stress tensor, mode 3 */
#define POLYMER_STRESS22_2 42 /* aka S22, 22 component of polymer stress tensor, mode 3 */
#define POLYMER_STRESS13_2 43 /* aka S13, 13 component of polymer stress tensor, mode 3 */
#define POLYMER_STRESS23_2 44 /* aka S23, 23 component of polymer stress tensor, mode 3 */
#define POLYMER_STRESS33_2 45 /* aka S33, 33 component of polymer stress tensor, mode 3 */

#define POLYMER_STRESS11_3 46 /* aka S11, 11 component of polymer stress tensor, mode 4 */
#define POLYMER_STRESS12_3 47 /* aka S12, 12 component of polymer stress tensor, mode 4 */
#define POLYMER_STRESS22_3 48 /* aka S22, 22 component of polymer stress tensor, mode 4 */
#define POLYMER_STRESS13_3 49 /* aka S13, 13 component of polymer stress tensor, mode 4 */
#define POLYMER_STRESS23_3 50 /* aka S23, 23 component of polymer stress tensor, mode 4 */
#define POLYMER_STRESS33_3 51 /* aka S33, 33 component of polymer stress tensor, mode 4 */

#define POLYMER_STRESS11_4 52 /* aka S11, 11 component of polymer stress tensor, mode 5 */
#define POLYMER_STRESS12_4 53 /* aka S12, 12 component of polymer stress tensor, mode 5 */
#define POLYMER_STRESS22_4 54 /* aka S22, 22 component of polymer stress tensor, mode 5 */
#define POLYMER_STRESS13_4 55 /* aka S13, 13 component of polymer stress tensor, mode 5 */
#define POLYMER_STRESS23_4 56 /* aka S23, 23 component of polymer stress tensor, mode 5 */
#define POLYMER_STRESS33_4 57 /* aka S33, 33 component of polymer stress tensor, mode 5 */

#define POLYMER_STRESS11_5 58 /* aka S11, 11 component of polymer stress tensor, mode 6 */
#define POLYMER_STRESS12_5 59 /* aka S12, 12 component of polymer stress tensor, mode 6 */
#define POLYMER_STRESS22_5 60 /* aka S22, 22 component of polymer stress tensor, mode 6 */
#define POLYMER_STRESS13_5 61 /* aka S13, 13 component of polymer stress tensor, mode 6 */
#define POLYMER_STRESS23_5 62 /* aka S23, 23 component of polymer stress tensor, mode 6 */
#define POLYMER_STRESS33_5 63 /* aka S33, 33 component of polymer stress tensor, mode 6 */

#define POLYMER_STRESS11_6 64 /* aka S11, 11 component of polymer stress tensor, mode 7 */
#define POLYMER_STRESS12_6 65 /* aka S12, 12 component of polymer stress tensor, mode 7 */
#define POLYMER_STRESS22_6 66 /* aka S22, 22 component of polymer stress tensor, mode 7 */
#define POLYMER_STRESS13_6 67 /* aka S13, 13 component of polymer stress tensor, mode 7 */
#define POLYMER_STRESS23_6 68 /* aka S23, 23 component of polymer stress tensor, mode 7 */
#define POLYMER_STRESS33_6 69 /* aka S33, 33 component of polymer stress tensor, mode 7 */

#define POLYMER_STRESS11_7 70 /* aka S11, 11 component of polymer stress tensor, mode 8 */
#define POLYMER_STRESS12_7 71 /* aka S12, 12 component of polymer stress tensor, mode 8 */
#define POLYMER_STRESS22_7 72 /* aka S22, 22 component of polymer stress tensor, mode 8 */
#define POLYMER_STRESS13_7 73 /* aka S13, 13 component of polymer stress tensor, mode 8 */
#define POLYMER_STRESS23_7 74 /* aka S23, 23 component of polymer stress tensor, mode 8 */
#define POLYMER_STRESS33_7 75 /* aka S33, 33 component of polymer stress tensor, mode 8 */

/*
 *  Following is an unrolling of the subvariable number for
 *  species concentrations in each domain. The subvariable number is now coded into
 *  the name of the species.
 *          A couple of specifications:
 *             The variable index should be contiguous and monatonically increasing
 *             with the coded subvar index.
 *             All subspecies in a material have contiguous variable types.
 *
 *  The species unknown variable types, SPECIES_UNK_#, may represent volume fractions,
 *  mass fractions, mole fractions, concentrations, or   pressures and/or liquid
 *  vol.frac. in porous media.
 *  -> Basically, this is a catchall for the species unknowns.
 */
#define SPECIES_UNK_0    76
#define SPECIES_UNK_1    77
#define SPECIES_UNK_2    78
#define SPECIES_UNK_3    79
#define SPECIES_UNK_4    80
#define SPECIES_UNK_5    81
#define SPECIES_UNK_6    82
#define SPECIES_UNK_7    83
#define SPECIES_UNK_8    84
#define SPECIES_UNK_9    85
#define SPECIES_UNK_10   86
#define SPECIES_UNK_11   87
#define SPECIES_UNK_12   88
#define SPECIES_UNK_13   89
#define SPECIES_UNK_14   90
#define SPECIES_UNK_15   91
#define SPECIES_UNK_16   92
#define SPECIES_UNK_17   93
#define SPECIES_UNK_18   94
#define SPECIES_UNK_19   95
#define SPECIES_UNK_20   96
#define SPECIES_UNK_21   97
#define SPECIES_UNK_22   98
#define SPECIES_UNK_23   99
#define SPECIES_UNK_24   100
#define SPECIES_UNK_25   101
#define SPECIES_UNK_26   102
#define SPECIES_UNK_27   103
#define SPECIES_UNK_28   104
#define SPECIES_UNK_29   105
#define SPECIES_UNK_LAST 105

/*
 *  MAX_SPECIES_UNK_NUM is the maximum number of species variables types
 *  coded into this file
 */
#define MAX_SPECIES_UNK_NUM 30

/*
 *  Subgrid volume fraction variables for phases. (Note these are
 *  thermodynamically distinct from species unknowns)
 *   (example -> precipated stoichiometrically-fixed salts in liquid
 *               solutions.
 *               These always have a mass fraction equal to one
 *               within the phase. However, the total amount of these
 *               phases (which you wouldn't want to resolve with the
 *               grid) would be specified by a nodal volume fraction.)
 *  These VOLF_PHASE_# variable types follow the same conventions as the
 *  SPECIES_UNK variable types
 */
#define VOLF_PHASE_0    106
#define VOLF_PHASE_1    107
#define VOLF_PHASE_2    108
#define VOLF_PHASE_3    109
#define VOLF_PHASE_4    110
#define VOLF_PHASE_LAST 110

/*
 *  MAX_VOLF_PHASE_NUM is the maximum number of subgrid volume fraction
 *  variables hardcoded into this file
 */
#define MAX_VOLF_PHASE_NUM 5

#define POR_LIQ_PRES   111 /* porous media solvent pressure */
#define POR_GAS_PRES   112 /* porous media gas pressure     */
#define POR_POROSITY   113 /* porous media porosity         */
#define POR_TEMP       114 /* porous media temperature      */
#define POR_SATURATION 115 /* porous media saturation       */
#define POR_LAST       115

/*
 *  MAX_POROUS_NUM is the maximum number of porous media
 *  variables hardcoded into this file
 */

#define MAX_POROUS_NUM 5

/*
 *  POR_SINK_MASS is a new porous-media-like variable that is used to
 *  model porous adsorbant media
 */
#define POR_SINK_MASS 116

/* MMH: This is the vorticity principle flow direction (not to be
 * confused with the vorticity vector), and the eigenvalue associated
 * with it.  AFAIK, this is only used for the QTENSOR Suspension
 * model...
 */
#define VORT_DIR1   117
#define VORT_DIR2   118
#define VORT_DIR3   119
#define VORT_LAMBDA 120

/*
 * TAB
 * This is the curvature of the level set field for computing
 * sharp surface tension effects and applying contact angles at the like
 */

#define CURVATURE 121

/*
 * These hold a generic lagrange multiplier field (vector or 3 scalars
 * for PDE-constrained problems.  cf. overlapping grid technique
 */
#define LAGR_MULT1 122
#define LAGR_MULT2 123
#define LAGR_MULT3 124

/*
 * This is a scalar related to the particle velocity fluctuations in
 * the suspension.
 */
#define BOND_EVOLUTION 125

#define SURF_CHARGE 126

#define EXT_VELOCITY                                             \
  127               /* Extension velocity in normal direction    \
                     * for level set simulations - it's a scalar \
                     */
#define EFIELD1 128 /* Electric field: E=grad(Voltage) */
#define EFIELD2 129
#define EFIELD3 130

/* MMH: The square of the norm of the potential field (grad V).
 * Needed for derivative requirements for dielectrophoresis
 * modeling. */
#define ENORM 131 /* ENORM = |E| */

#define NORMAL1 132 /* normal to LS field, component 1 */
#define NORMAL2 133 /* normal to LS field, component 2 */
#define NORMAL3 134 /* normal to LS field, component 3 */

#define SHELL_CURVATURE  135
#define SHELL_CURVATURE2 136
#define SHELL_TENSION    137
#define SHELL_X          138
#define SHELL_Y          139
#define SHELL_USER       140
#define PHASE1           141
#define PHASE2           142
#define PHASE3           143
#define PHASE4           144
#define PHASE5           145
#define SHELL_ANGLE1     146 /*  Smoothed normal orientation computed using shell elements, theta */
#define SHELL_ANGLE2     147 /*  Smoothed normal orientation computed using shell elements, phi */
#define SHELL_SURF_DIV_V 148 /*  (I-nn).div(V) ..From Surface Rheo Constitutive model */
#define SHELL_SURF_CURV \
  149 /*  2H = div_s(n), mean curvature ..From Surface Rheo Constitutive model */
#define N_DOT_CURL_V \
  150 /*  curl(V).n , higher-order quantity ..From Surface Rheo Constitutive model */
#define GRAD_S_V_DOT_N1 \
  151 /*  grad_s(v).n - x , higher-order quantity ..From Surface Rheo Constitutive model */
#define GRAD_S_V_DOT_N2 \
  152 /*  grad_s(v).n - y , higher-order quantity ..From Surface Rheo Constitutive model */
#define GRAD_S_V_DOT_N3 \
  153 /*  grad_s(v).n - z , higher-order quantity ..From Surface Rheo Constitutive model */
#define ACOUS_PREAL          154 /*  Acoustic harmonic pressure - real part    */
#define ACOUS_PIMAG          155 /*  Acoustic harmonic pressure - imag part    */
#define SHELL_DIFF_FLUX      156 /* sh_J - material diffusion on a shell surface */
#define SHELL_DIFF_CURVATURE 157 /* K = div_s_n used with shell normal vector unknowns */
#define SHELL_NORMAL1        158 /* X-component of shell normal vector unknown */
#define SHELL_NORMAL2        159 /* Y-component of shell normal vector unknown */
#define SHELL_NORMAL3        160 /* Z-component of shell normal vector unknown */
#define ACOUS_REYN_STRESS    161 /*  Acoustic Reynolds Stress    */
#define SHELL_BDYVELO        162 /*  Acoustic Boundary Velocity Squared    */
#define SHELL_LUBP           163 /*  Shell Lubrication Pressure    */
#define LUBP                 164 /*  Lubrication Pressure    */
#define SHELL_FILMP          165 /*  Lubrication pressure of a thin film */
#define SHELL_FILMH          166 /*  Film thickness */
#define SHELL_PARTC          167 /*  Particle concentration in the thin film */
#define SHELL_SAT_CLOSED     168 /*  Structured porous shell saturation - closed pores - SAR */
#define SHELL_PRESS_OPEN     169 /*  Structured porous shell saturation - open pores - SAR */
#define SHELL_TEMPERATURE    170 /*  Shell temperature */
#define SHELL_DELTAH         171 /*  Lubrication shell thickness change (melting problems)*/
#define SHELL_LUB_CURV       172 /*  Curvature from the level set field in a lubrication shell - SAR */
#define SHELL_SAT_GASN       173 /*  Structured porous shell saturation - gas compression - SAR */
#define SHELL_SHEAR_TOP \
  174 /*  Shear rate at top wall - lubrication shell for generalized Newtonian flow */
#define SHELL_SHEAR_BOT \
  175 /*  Shear rate at bottom wall - lubrication shell for generalized Newtonian flow */
#define SHELL_CROSS_SHEAR \
  176 /*  Cross stream shear stress - lubrication shell for generalized Newtonian flow */
#define MAX_STRAIN         177 /*  Time-history maximum von mises strain in a solid material */
#define CUR_STRAIN         178 /*  Current value of the von mises strain in a solid material */
#define LUBP_2             179 /*  Second Lubrication Pressure for multilayer problems  */
#define SHELL_PRESS_OPEN_2 180 /*  Second Structured porous shell pressure - open pores - SAR */
#define SHELL_LUB_CURV_2   181 /*  Curvature from the phase field in a lubrication_2 shell - PRS*/
#define LIGHT_INTP         182 /*  Light Intensity - Plus direction propagation-RBS*/
#define LIGHT_INTM         183 /*  Light Intensity - Minus direction propagation-RBS*/
#define LIGHT_INTD         184 /*  Light Intensity - Scattering Dispersion-RBS*/
#define TFMP_SAT           185 /*  Thin-Film Multi-Phase Saturation */
#define TFMP_PRES          186 /*  Thin-Film Multi-Phase Lubrication Pressure */
#define RESTIME            187 /*  Residence Time Function */
#define SHELL_SAT_1        188 /*  Porous shell layer 1 */
#define SHELL_SAT_2        189 /*  Porous shell layer 2 */
#define SHELL_SAT_3        190 /*  Porous shell layer 3 */
#define EM_E1_REAL         191 /*  EM wave variables */
#define EM_E2_REAL         192
#define EM_E3_REAL         193
#define EM_E1_IMAG         194 /*  EM wave variables */
#define EM_E2_IMAG         195
#define EM_E3_IMAG         196
#define EM_H1_REAL         197 /*  EM wave variables */
#define EM_H2_REAL         198
#define EM_H3_REAL         199
#define EM_H1_IMAG         200 /*  EM wave variables */
#define EM_H2_IMAG         201
#define EM_H3_IMAG         202
#define EM_CONT_REAL       203
#define EM_CONT_IMAG       204
#define MOMENT0            205
#define MOMENT1            206
#define MOMENT2            207
#define MOMENT3            208
#define DENSITY_EQN        209
#define PSTAR              210
#define USTAR              211
#define VSTAR              212
#define WSTAR              213
#define EDDY_NU            214
/*
 * define a variable to hold an external field which will be
 * held fixed in the problem but parametered by the basis functions
 * (aka "B-field" or something like that) */
#define EXTERNAL 0

/*
 * Define some special constants for input of boundary conditions
 * (not used in  setup and ordering of equations and sensitivities !!)
 *
 * requirements for numbering:
 *     All of the section below must be larger than V_LAST
 *     R_MOM_NORMAL  < R_MOM_TANG1  < R_MOM_TANG2       (and contiguous)
 *     R_MESH_NORMAL < R_MESH_TANG1 < R_MESH_TANG2      (and contiguous)
 *     MESH_POSITION1 < MESH_POSITION2 < MESH_POSITION3 (and contiguous)
 *     VEL_NORM < VEL_TANG1 < VEL_TANG2                 (and contiguous)
 *     MESH_NORM < MESH_TANG1 < MESH_TANG2              (and contiguous)
 *     D_VEL1_DT to D_P_DT -> monotonically contiguous and
 *                            numbered such that they correspond
 *                            to the associated variables at the
 *                            top of the variable list up to an
 *                            arbitrary additive offset.
 */
#define MESH_POSITION1  1147 /* aka "D1" + COOR */
#define MESH_POSITION2  1148 /* aka "D1" + COOR */
#define MESH_POSITION3  1149 /* aka "D1" + COOR */
#define R_MOM_NORMAL    1150 /* aka normal fluid momentum equation */
#define R_MOM_TANG1     1151 /* aka first tangent momentum equation */
#define R_MOM_TANG2     1152 /* aka second tangent momentum equation */
#define R_MESH_NORMAL   1155 /* aka normal mesh equation */
#define R_MESH_TANG1    1156 /* aka first tangent mesh equation */
#define R_MESH_TANG2    1157 /* aka second tangent mesh equation */
#define VEL_NORM        1150 /* aka normal velocity */
#define VEL_TANG1       1151 /* aka first normal velocity */
#define VEL_TANG2       1152 /* aka second normal velocity */
#define MESH_NORM       1155 /* aka mesh normal first component ?? */
#define MESH_TANG1      1156 /* aka first mesh tangent ?? */
#define MESH_TANG2      1157 /* aka second mesh tangent  ?? */
#define SOLID_POSITION1 1158 /* aka "D1" + COOR */
#define SOLID_POSITION2 1159 /* aka "D1" + COOR */
#define SOLID_POSITION3 1160 /* aka "D1" + COOR */
#define R_SOLID_NORMAL  1161 /* aka normal solid equation */
#define R_SOLID_TANG1   1162 /* aka first tangent solid equation */
#define R_SOLID_TANG2   1163 /* aka second tangent solid equation */
#define SOLID_NORM      1164 /* aka solid normal first component ?? */
#define SOLID_TANG1     1165 /* aka first solid tangent component */
#define SOLID_TANG2     1166 /* aka second solid tangent  component */
#define SPEED           1167 /* aka velocity magnitude */

#define R_ANYTHING -1 /* aka for bc's which can apply to any equation */

#define D_VEL1_DT 1260 /* aka "dU1/dt" */
#define D_VEL2_DT 1261 /* aka "dU2/dt" */
#define D_VEL3_DT 1262 /* aka "dU3/dt" */
#define D_T_DT    1263 /* aka "dT/dt" */
#define D_Y_DT    1264 /* aka "dY/dt" */
#define D_X1_DT   1265 /* aka "dD1/dt" */
#define D_X2_DT   1266 /* aka "dD2/dt" */
#define D_X3_DT   1267 /* aka "dD3/dt" */
#define D_S_DT    1268 /* aka "dS/dt" */
#define D_P_DT    1269 /* aka "dP/dt" */

/*
 *  Define some special definitions to provide explicitness as to what the species
 *  unknowns actually represent! Note, the definitions below are used to refine
 *  the definitions of the MASS_FRACTION and SPECIES_UNK_# variables types.
 *
 *    These are the acceptable values for the Species_Var_Type entries in the
 *  Problem_Description and Uniform_Problem_Description Structures.
 *
 *   Requirements for numbering:
 *       -> None that I know of
 */
#define SPECIES_MASS_FRACTION                     \
  2170 /* The Independent species unknown is      \
          mass fraction. The species equation has \
          units of mass/(cm^3 sec)  */
       /* The exodus unknowns have prefixes of Yk_ */
#define SPECIES_MOLE_FRACTION                     \
  2171 /* The Independent species unknown is      \
          mole fraction. The species equation has \
          units of mol/(cm^3 sec)  */
       /* The exodus unknowns have prefixes of Xk_ */
#define SPECIES_VOL_FRACTION                        \
  2172 /* The Independent species unknown is        \
          volume fraction. The species equation has \
          units of mol/(cm^3 sec)  */
       /*	The exodus unknowns have prefixes of Vk_ */
#define SPECIES_DENSITY                                 \
  2173 /* The independent species unknown is species    \
          density with units of  mass/cm^3 or mass/cm^2 \
          for surface species. The species equation     \
          has units of mass/(cm^3 sec) . */
       /*	The exodus unknowns have prefixes of Dk_ */
#define SPECIES_CONCENTRATION             \
  2174 /* Units of moles/cm^3 or mol/cm^2 \
          for surface species */
       /* The exodus unknowns have prefixes of Ck_ */
#define SPECIES_CAP_PRESSURE                                                      \
  2175                              /* Units of dynes/cm^2 -> this represents     \
                                       a species concentration when combined with \
                                       equilibrium data. The species equation has \
                                       units of mol/(cm^3 sec)  */
                                    /* The exodus unknowns have prefixes of Pk_ */
#define SPECIES_UNDEFINED_FORM 2176 /* Old form (default) here for compatibility */
                                    /* The exodus unknowns have prefixes of Y
                                       and consist of integer names starting with zero */

/*
 * V_FIRST and V_LAST:
 *      This facilitates the use of loops like...
 *	for (v=V_FIRST; v<V_LAST; v++)
 *
 *  HKM NOTE -> These types of loops will be taken out of all
 *              time-critical portions of the code (i.e., all
 *              residual and matrix fill routines)
 *              They are too slow, with the addition of all of
 *              the variable types created by unrolling the
 *              subvariable index.
 *
 * V_LAST is defined at the end of the #define R_* statements.
 *
 */
#define V_FIRST 0

/*
 *
 *     Equation Names, residuals...
 *
 *         MMH NOTE -> Should we add another vector eq for particle momentum?
 */
#define NUM_VECTOR_EQUATIONS 2
#define VECT_EQ_MOM          0
#define VECT_EQ_MESH         1

#define R_MOMENTUM1 0
#define R_MOMENTUM2 1
#define R_MOMENTUM3 2
#define R_ENERGY    3
#define R_MASS      4
#define R_MESH1     5
#define R_MESH2     6
#define R_MESH3     7
#define R_MASS_SURF 8
#define R_PRESSURE  9

#define R_STRESS11 10 /* aka S11 residual */
#define R_STRESS12 11 /* aka S12 residual */
#define R_STRESS22 12 /* aka S22 residual */
#define R_STRESS13 13 /* aka S13 residual */
#define R_STRESS23 14 /* aka S23 residual */
#define R_STRESS33 15 /* aka S33 residual */

#define R_SOLID1 16 /*AKA Real solid momentum residual, wrt d'*/
#define R_SOLID2 17 /*AKA Real solid momentum residual, wrt d'*/
#define R_SOLID3 18 /*AKA Real solid momentum residual, wrt d'*/

#define R_GRADIENT11 19 /* aka G11 residual */
#define R_GRADIENT12 20 /* aka G12 residual */
#define R_GRADIENT21 21 /* aka G21 residual */
#define R_GRADIENT22 22 /* aka G22 residual */
#define R_GRADIENT13 23 /* aka G13 residual */
#define R_GRADIENT23 24 /* aka G23 residual */
#define R_GRADIENT31 25 /* aka G31 residual */
#define R_GRADIENT32 26 /* aka G32 residual */
#define R_GRADIENT33 27 /* aka G33 residual */

#define R_POTENTIAL  28     /* aka V voltage potential residual */
#define R_FILL       29     /* aka F fill residual */
#define R_LEVEL_SET  R_FILL /* level set equation identical to fill */
#define R_SHEAR_RATE 30     /* aka SH shear rate from 2nd invariant residual */

#define R_PMOMENTUM1 31 /* aka Particle momentum residual */
#define R_PMOMENTUM2 32 /* aka Particle momentum residual */
#define R_PMOMENTUM3 33 /* aka Particle momentum residual */

/* Please leave these stress variables continuously numbered or
 * multimode viscoelasticity will never work!
 */
#define R_STRESS11_1 34 /* aka S11, 11 component of polymer stress tensor, mode 2 */
#define R_STRESS12_1 35 /* aka S12, 12 component of polymer stress tensor, mode 2 */
#define R_STRESS22_1 36 /* aka S22, 22 component of polymer stress tensor, mode 2 */
#define R_STRESS13_1 37 /* aka S13, 13 component of polymer stress tensor, mode 2 */
#define R_STRESS23_1 38 /* aka S23, 23 component of polymer stress tensor, mode 2 */
#define R_STRESS33_1 39 /* aka S33, 33 component of polymer stress tensor, mode 2 */

#define R_STRESS11_2 40 /* aka S11, 11 component of polymer stress tensor, mode 3 */
#define R_STRESS12_2 41 /* aka S12, 12 component of polymer stress tensor, mode 3 */
#define R_STRESS22_2 42 /* aka S22, 22 component of polymer stress tensor, mode 3 */
#define R_STRESS13_2 43 /* aka S13, 13 component of polymer stress tensor, mode 3 */
#define R_STRESS23_2 44 /* aka S23, 23 component of polymer stress tensor, mode 3 */
#define R_STRESS33_2 45 /* aka S33, 33 component of polymer stress tensor, mode 3 */

#define R_STRESS11_3 46 /* aka S11, 11 component of polymer stress tensor, mode 4 */
#define R_STRESS12_3 47 /* aka S12, 12 component of polymer stress tensor, mode 4 */
#define R_STRESS22_3 48 /* aka S22, 22 component of polymer stress tensor, mode 4 */
#define R_STRESS13_3 49 /* aka S13, 13 component of polymer stress tensor, mode 4 */
#define R_STRESS23_3 50 /* aka S23, 23 component of polymer stress tensor, mode 4 */
#define R_STRESS33_3 51 /* aka S33, 33 component of polymer stress tensor, mode 4 */

#define R_STRESS11_4 52 /* aka S11, 11 component of polymer stress tensor, mode 5 */
#define R_STRESS12_4 53 /* aka S12, 12 component of polymer stress tensor, mode 5 */
#define R_STRESS22_4 54 /* aka S22, 22 component of polymer stress tensor, mode 5 */
#define R_STRESS13_4 55 /* aka S13, 13 component of polymer stress tensor, mode 5 */
#define R_STRESS23_4 56 /* aka S23, 23 component of polymer stress tensor, mode 5 */
#define R_STRESS33_4 57 /* aka S33, 33 component of polymer stress tensor, mode 5 */

#define R_STRESS11_5 58 /* aka S11, 11 component of polymer stress tensor, mode 6 */
#define R_STRESS12_5 59 /* aka S12, 12 component of polymer stress tensor, mode 6 */
#define R_STRESS22_5 60 /* aka S22, 22 component of polymer stress tensor, mode 6 */
#define R_STRESS13_5 61 /* aka S13, 13 component of polymer stress tensor, mode 6 */
#define R_STRESS23_5 62 /* aka S23, 23 component of polymer stress tensor, mode 6 */
#define R_STRESS33_5 63 /* aka S33, 33 component of polymer stress tensor, mode 6 */

#define R_STRESS11_6 64 /* aka S11, 11 component of polymer stress tensor, mode 7 */
#define R_STRESS12_6 65 /* aka S12, 12 component of polymer stress tensor, mode 7 */
#define R_STRESS22_6 66 /* aka S22, 22 component of polymer stress tensor, mode 7 */
#define R_STRESS13_6 67 /* aka S13, 13 component of polymer stress tensor, mode 7 */
#define R_STRESS23_6 68 /* aka S23, 23 component of polymer stress tensor, mode 7 */
#define R_STRESS33_6 69 /* aka S33, 33 component of polymer stress tensor, mode 7 */

#define R_STRESS11_7 70 /* aka S11, 11 component of polymer stress tensor, mode 8 */
#define R_STRESS12_7 71 /* aka S12, 12 component of polymer stress tensor, mode 8 */
#define R_STRESS22_7 72 /* aka S22, 22 component of polymer stress tensor, mode 8 */
#define R_STRESS13_7 73 /* aka S13, 13 component of polymer stress tensor, mode 8 */
#define R_STRESS23_7 74 /* aka S23, 23 component of polymer stress tensor, mode 8 */
#define R_STRESS33_7 75 /* aka S33, 33 component of polymer stress tensor, mode 8 */

#define R_SPECIES_UNK_0    76
#define R_SPECIES_UNK_1    77
#define R_SPECIES_UNK_2    78
#define R_SPECIES_UNK_3    79
#define R_SPECIES_UNK_4    80
#define R_SPECIES_UNK_5    81
#define R_SPECIES_UNK_6    82
#define R_SPECIES_UNK_7    83
#define R_SPECIES_UNK_8    84
#define R_SPECIES_UNK_9    85
#define R_SPECIES_UNK_10   86
#define R_SPECIES_UNK_11   87
#define R_SPECIES_UNK_12   88
#define R_SPECIES_UNK_13   89
#define R_SPECIES_UNK_14   90
#define R_SPECIES_UNK_15   91
#define R_SPECIES_UNK_16   92
#define R_SPECIES_UNK_17   93
#define R_SPECIES_UNK_18   94
#define R_SPECIES_UNK_19   95
#define R_SPECIES_UNK_20   96
#define R_SPECIES_UNK_21   97
#define R_SPECIES_UNK_22   98
#define R_SPECIES_UNK_23   99
#define R_SPECIES_UNK_24   100
#define R_SPECIES_UNK_25   101
#define R_SPECIES_UNK_26   102
#define R_SPECIES_UNK_27   103
#define R_SPECIES_UNK_28   104
#define R_SPECIES_UNK_29   105
#define R_SPECIES_UNK_LAST 105

#define R_VOLF_PHASE_0    106
#define R_VOLF_PHASE_1    107
#define R_VOLF_PHASE_2    108
#define R_VOLF_PHASE_3    109
#define R_VOLF_PHASE_4    110
#define R_VOLF_PHASE_LAST 110

#define R_POR_LIQ_PRES   111
#define R_POR_GAS_PRES   112
#define R_POR_POROSITY   113
#define R_POR_SATURATION 114
#define R_POR_ENERGY     115
#define R_POR_LAST       115

#define R_POR_SINK_MASS 116

#define R_VORT_DIR1   117
#define R_VORT_DIR2   118
#define R_VORT_DIR3   119
#define R_VORT_LAMBDA 120

#define R_CURVATURE 121

#define R_LAGR_MULT1 122
#define R_LAGR_MULT2 123
#define R_LAGR_MULT3 124

/* This is a scalar equation related to the particle velocity
 * fluctuations in the suspension.
 */
#define R_BOND_EVOLUTION 125

#define R_SURF_CHARGE 126

#define R_EXT_VELOCITY                             \
  127 /* Extension velocity in normal direction    \
       * for level set simulations - it's a scalar \
       */

#define R_EFIELD1 128 /* Electric field: E=grad(Voltage) */
#define R_EFIELD2 129
#define R_EFIELD3 130

#define R_ENORM 131

#define R_NORMAL1 132 /*  Normal to LS field, component 1 */
#define R_NORMAL2 133 /*  Normal to LS field, component 2 */
#define R_NORMAL3 134 /*  Normal to LS field, component 3 */

/* Shell element variables */
#define R_SHELL_CURVATURE  135 /*  Curvature for structural shell elements*/
#define R_SHELL_CURVATURE2 136 /*  Second curvature for structural shell elements*/
#define R_SHELL_TENSION    137 /*  Tension variable for structural shell elements */
#define R_SHELL_X          138 /*  Shell X position equation */
#define R_SHELL_Y          139 /*  Shell Y position equation */
#define R_SHELL_USER       140 /*  Shell user equation */
#define R_PHASE1           141 /* First phase function */
#define R_PHASE2           142 /* Second phase function */
#define R_PHASE3           143 /* Third phase function */
#define R_PHASE4           144 /* Fourth phase function */
#define R_PHASE5           145 /* Fifth phase function */

/* more shell element vars */
#define R_SHELL_ANGLE1 146 /*  Smoothed normal orientation computed using shell elements, theta */
#define R_SHELL_ANGLE2 147 /*  Smoothed normal orientation computed using shell elements, phi */

/*Still more shell element vars */
#define R_SHELL_SURF_DIV_V 148 /*  (I-nn).div(V) ..From Surface Rheo Constitutive model */
#define R_SHELL_SURF_CURV \
  149 /*  2H = - div_s(n), mean curvature ..From Surface Rheo Constitutive model */
#define R_N_DOT_CURL_V \
  150 /*  curl(V).n , higher-order quantity ..From Surface Rheo Constitutive model */
#define R_GRAD_S_V_DOT_N1 \
  151 /*  grad_s(v).n - x , higher-order quantity ..From Surface Rheo Constitutive model */
#define R_GRAD_S_V_DOT_N2 \
  152 /*  grad_s(v).n - y , higher-order quantity ..From Surface Rheo Constitutive model */
#define R_GRAD_S_V_DOT_N3 \
  153 /*  grad_s(v).n - z , higher-order quantity ..From Surface Rheo Constitutive model */
#define R_ACOUS_PREAL          154 /*  Acoustic harmonic wave equation - real part */
#define R_ACOUS_PIMAG          155 /*  Acoustic harmonic wave equation - imaginary part */
#define R_SHELL_DIFF_FLUX      156 /* sh_J - material diffusion on a shell surface */
#define R_SHELL_DIFF_CURVATURE 157 /* K = div_s_n used with shell normal vector unknowns */
#define R_SHELL_NORMAL1        158 /* X-component of shell normal vector unknown */
#define R_SHELL_NORMAL2        159 /* Y-component of shell normal vector unknown */
#define R_SHELL_NORMAL3        160 /* Z-component of shell normal vector unknown */
#define R_ACOUS_REYN_STRESS    161 /*  Acoustic Reynolds Stress */
#define R_SHELL_BDYVELO        162 /*  Acoustic Streaming Velocity shell */
#define R_SHELL_LUBP           163 /*  Lubrication Pressure shell */
#define R_LUBP                 164 /*  Lubrication pressure Reynolds equation not on shell */
#define R_SHELL_FILMP          165 /*  Lubrication pressure of a thin film */
#define R_SHELL_FILMH          166 /*  Film thickness */
#define R_SHELL_PARTC          167 /*  Particles concentration */
#define R_SHELL_SAT_CLOSED     168 /*  Structured porous shell - Closed pores - SAR  */
#define R_SHELL_SAT_OPEN       169 /*  Structured porous shell - Open pores - SAR  */
#define R_SHELL_ENERGY         170 /*  Shell energy equation */
#define R_SHELL_DELTAH         171 /*  Shell gap evolution equation (e.g. melting problems) */
#define R_SHELL_LUB_CURV                                                                         \
  172                        /*  Curvature from the level set field in a lubrication shell - SAR \
                              */
#define R_SHELL_SAT_GASN 173 /*  Structured porous shell - Gas compression - SAR */
#define R_SHELL_SHEAR_TOP \
  174 /*  Shear rate at top wall - lubrication shell for generalized Newtonian flow */
#define R_SHELL_SHEAR_BOT \
  175 /*  Shear rate at bottom wall - lubrication shell for generalized Newtonian flow */
#define R_SHELL_CROSS_SHEAR \
  176 /*  Cross stream shear stress - lubrication shell for generalized Newtonian flow */
#define R_MAX_STRAIN       177 /*  Time-history maximum von mises strain in a solid material */
#define R_CUR_STRAIN       178 /*  Current value of the von mises strain in a solid material */
#define R_LUBP_2           179 /*  Second lubrication_pressure for multilayer flows */
#define R_SHELL_SAT_OPEN_2 180 /*  Second open pore shell lubrication saturation */
#define R_SHELL_LUB_CURV_2                                                                   \
  181                      /*  Curvature from the phase field in a lubrication_2 shell - PRS \
                            */
#define R_LIGHT_INTP   182 /*  Light Intensity - Plus direction propagation*/
#define R_LIGHT_INTM   183 /*  Light Intensity - Minus direction propagation*/
#define R_LIGHT_INTD   184 /*  Light Intensity - Scattering Dispersion*/
#define R_TFMP_MASS    185 /*  Thin-Film Multi-Phase Mass Equation */
#define R_TFMP_BOUND   186 /*  Thin-Film Multi-Phase Boundary Motion Equation */
#define R_RESTIME      187 /*  Resdience Time Function */
#define R_SHELL_SAT_1  188 /*  Porous shell layer 1 */
#define R_SHELL_SAT_2  189 /*  Porous shell layer 1 */
#define R_SHELL_SAT_3  190 /*  Porous shell layer 1 */
#define R_EM_E1_REAL   191 /*  EM wave variables */
#define R_EM_E2_REAL   192
#define R_EM_E3_REAL   193
#define R_EM_E1_IMAG   194 /*  EM wave variables */
#define R_EM_E2_IMAG   195
#define R_EM_E3_IMAG   196
#define R_EM_H1_REAL   197 /*  EM wave variables */
#define R_EM_H2_REAL   198
#define R_EM_H3_REAL   199
#define R_EM_H1_IMAG   200 /*  EM wave variables */
#define R_EM_H2_IMAG   201
#define R_EM_H3_IMAG   202
#define R_EM_CONT_REAL 203
#define R_EM_CONT_IMAG 204
#define R_MOMENT0      205
#define R_MOMENT1      206
#define R_MOMENT2      207
#define R_MOMENT3      208
#define R_DENSITY_EQN  209
#define R_PSTAR        210
#define R_USTAR        211
#define R_VSTAR        212
#define R_WSTAR        213
#define R_EDDY_NU      214
#define V_LAST         215

/* MMH
 * This is used for those parts of the code that want to ensure
 * a "real" equation is being calculated, not just a boundary
 * condition.  It used to be that if eqn <= R_GRADIENT33, eqn
 * was assumed to be a real equation.  Now you only have to change
 * it here when you add another new equation to the list.  It seems
 * to me that it should just be equal to V_LAST, so I will set it
 * to that.  If you want it to be something else, go ahead.
 *
 * NOTE: The R_POTENTIAL,R_FILL, and R_SHEAR_RATE were NOT considered
 * a "real" equation before this change.  Now they are.  I don't know
 * if this is going to break anything.
 */
#define LAST_REAL_EQ V_LAST

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

/* Selection of Time Integration Methods */

#define STEADY    0 /* Steady state solution method               */
#define TRANSIENT 1 /* Accurate transient solution method         */

/* Selection of Continuation Method Order */

#define ALC_NONE   -1 /* No continuation */
#define ALC_ZEROTH 0  /* Zero order   */
#define ALC_FIRST  1  /* First order  */
#define ALC_SECOND 2  /* Second order */
#define ALC_THIRD  3  /* Third order  */
#define ALC_FOURTH 4  /* Fourth order */
#define LOCA       9  /* Use LOCA     */

/* Selection of LOCA continuation algorithms */

#define PSEUDO          2 /* Pseudo transient method (not time accurate)*/
#define CONTINUATION    3 /* Continuation of steady state solutions     */
#define TP_CONTINUATION 4 /* Continuation of turning points             */
#define PF_CONTINUATION 7 /* Continuation of pitchfork bifurcations     */
#define HP_CONTINUATION 8 /* Continuation of Hopf bifurcations          */
#define LOCA_LSA_ONLY   9 /* Access ARPACK through LOCA w/o continuation*/

/* Selection of Hunting Continuation Method Order */

#define HUN_ZEROTH 100 /* Zero order   */
#define HUN_FIRST  101 /* First order  */

/* Selection augmenting condition type */

#define AC_USERBC        1
#define AC_USERMAT       2
#define AC_VOLUME        3
#define AC_FLUX          4
#define AC_LGRM          5
#define AC_LS_VEL        6
#define AC_ARC_LENGTH    7
#define AC_OVERLAP       8
#define AC_PERIODIC      9
#define AC_PF_CONSTRAINT 10
#define AC_POSITION      11 // Positional constraint
#define AC_FLUX_MAT      12
#define AC_ANGLE         13 // Angle constraint - 2D onl

/* Post Processing Volumetric Integration types - AMC Originally Aug 2013 */

#define PPVI_VERBOSE 1 //  This corresponds to the original output type
#define PPVI_CSV \
  2 //  This will output a properly labled csv file, useful for automated python interpretation

/* table interpolation types */

#define LINEAR       3
#define QUADRATIC    2
#define QUAD_GP      1
#define BIQUADRATIC  21
#define TRILINEAR    22
#define TRIQUADRATIC 23

/* maximum definitions */

/* If you haven't missed these since 8/97, then get rid of them! */

/* #define MAX_NODES          	2000000 */
/* #define MAX_ELEMS          	2000000*/

/* #define MAX_NEIGHBORS      	50 */

/*
 * New constructs permit easier overriding during compilation. Eg.,
 *
 *	cc -DMAX_CONC=39 file.c
 */

#ifndef MAX_VARIABLE_TYPES
#define MAX_VARIABLE_TYPES V_LAST
#endif

#ifndef MAX_INTERPOLATIONS
#define MAX_INTERPOLATIONS 5 /* on any given shape element */
#endif

#ifndef MAX_ELEMENT_SHAPES
#define MAX_ELEMENT_SHAPES 6 /* 1D, 2D, 3D together... */
#endif

#ifndef MAX_BASIS_FUNCTIONS
#define MAX_BASIS_FUNCTIONS 3
#endif

#ifndef MAX_ELEMENT_TYPES
#define MAX_ELEMENT_TYPES 3
#endif

#ifndef MAX_TERM_TYPES /* mass, advection, boundary, diffusion, etc. */
#define MAX_TERM_TYPES 10
#endif
/*
 *
 * PAS - "Use -DMAX_CONC=20 compiler flag instead; leave this default as low
 *        as reasonably possible to avoid showing our bloated slow
 *        implementation of concentration..."
 *
 */

#ifndef MAX_CONC
#define MAX_CONC 4 /* Max number bulk, surface species */
#endif

#ifndef MAX_PMV
#define MAX_PMV 4 /* Max number porous variables */
#endif

#ifndef MAX_POR_SHELL
#define MAX_POR_SHELL 3 /* Max number porous shell equations in a block */
#endif

#ifndef MAX_RXN
#define MAX_RXN 5 /* Max number of reactions */
#endif

#ifndef MAX_EXTERNAL_FIELD
#define MAX_EXTERNAL_FIELD 7 /* Maximum number of external fields */
#endif

#define MAX_EQNS MAX_VARIABLE_TYPES

#ifndef MAX_GEOMETRY_PARMS
#define MAX_GEOMETRY_PARMS 8 /* Used for "mp" structs. */
#endif

#ifndef MAX_ELEMENT_INDICES_RELATED
#define MAX_ELEMENT_INDICES_RELATED 2
#endif

#ifndef MAX_NUM_MATRICES
#define MAX_NUM_MATRICES \
  11 /* Maximum number of matrices to be solved in segregated solver fashion */
#endif

#define MAX_MOMENTS 4

/*
 * Magic numbers for adaptive time step selection -- how much acceleration,
 * deceleration, and maximum time step size.
 */

#ifndef TIME_STEP_ALPHA
#define TIME_STEP_ALPHA (5.0)
#endif

#ifndef TIME_STEP_BETA
#define TIME_STEP_BETA (1.5)
#endif

#ifndef TIME_STEP_GROWTH_CAP
#define TIME_STEP_GROWTH_CAP (1.5)
#endif

/*
 * This moves here from el_elm.h. The maximum number of degrees of freedom
 * for a variable within an element. The value of 12 is chosen because it
 * is the lowest number that satifies a typical bound found for 2D
 * problems using Lagrangian  biquadratic interpolation on
 * quadrilateral elements with an extra 3 dof
 * for phase-jumping accounting at elements on cross-phase boundaries added
 * in.  This number  needs to be made larger if, for example, you want
 * to use 27 dof/elem  interpolations in 3D. HEX/8 elements should be
 * handled just fine with the default number, though.
 */

#ifndef MDE
#define MDE                                         \
  12 /* maximum number of unknowns for a particular \
      * variable type per element */
#endif

#ifndef MAX_DOF_PER_NODE
#define MAX_DOF_PER_NODE                     \
  3 /* Maximum number of unknowns for a      \
     * single variable type at a single node \
     */
#endif

#ifndef MAX_DOFNAME
#define MAX_DOFNAME 80 /* Max string length dof description. */
#endif

#ifndef MAX_MATLNAME
#define MAX_MATLNAME 85 /* Ought to suffice in most cases. */
#endif

#ifndef MAX_PROB_VAR
#define MAX_PROB_VAR                     \
  14 /* maximum number of variable types \
      * in any one problem  */
#endif

#ifndef MAX_PROB_EQN
#define MAX_PROB_EQN                               \
  MAX_PROB_VAR /* maximum number of equation types \
                * in any one problem */
#endif

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/
#endif
