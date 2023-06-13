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
 *$Id: rf_bc_const.h,v 5.16 2010-04-07 22:27:00 prschun Exp $
 */

/*
 * Here's a RECIPE for adding new boundary conditions so you don't have any
 * excuses not to add new ones.  The changes should be made in at least
 * four files (rf_bc_const.h, mm_names.h, mm_input.c, and bc_[method].c)
 * for some boundary conditions you may want to make additions elsewhere also.
 * One example of extra additions is in el_exoII_io.c where the ss_dup_list
 * is created - you may want to adapt the logic for rotation of new conditions.
 * (note that these lines are repeated at each place where you need to
 * make changes):
 *  Boundary Condition  Step 1: add Macro Substitution for your new boundary
 *                              condition in rf_bc_const.h - this is how you
 *      rf_bc_const.h           will refer to the new boundary condition
 *                              throughout the code.  Make sure the integer
 *                              you choose is unique from all the other BC
 *                              types.
 *
 *  Boundary Condition  Step 2: add a description of your boundary condition
 *                              to the BC_Desc structure array in mm_names.h.
 *      mm_names.h              This structure includes information about the
 *                              type of boundary condition, which equation it
 *                              applies to, what variables it is sensitive to,
 *                              whether to rotate mesh or momentum equations,
 *                              etc.  It is very important that you fill out
 *                              this structure carefully, otherwise the code
 *                              won't know what to do.
 *
 *  Boundary Condition  Step 3: add your BC case to the correct format listing
 *                              for reading the necessary arguments from the input
 *      mm_input.c              file in mm_input.c.
 *
 *  Boundary Condition  Step 3a: some BC's may require changes in mm_bc.c
 *                              similar to what is done for PLANE, GEOM, CA,
 *                              and all dirichlet conditions
 *
 *  Boundary Condition  Step 4: Add a function call (and a function) in the
 *                              correct routine for evaluating your boundary
 *      bc_colloc.c             condition.  This will probably be in bc_colloc.c
 *      bc_integ.c              for collocated conditions or bc_integ.c for
 *                              strong or weak integrated conditions.
 *
 *  Boundary Condition  Step 5: use and enjoy your new boundary condition
 *
 */

#ifndef GOMA_RF_BC_CONST_H
#define GOMA_RF_BC_CONST_H

#include "dpi.h"          /* To know about Dpi type. */
#include "exo_struct.h"   /* To know about Exo_DB type. */
#include "rf_fem_const.h" /* To know about MAX_VARIABLE_TYPES */

#ifndef EXTERN
#define EXTERN extern
#endif

/* lengths used in boundary condition structures and I/O */

/*
 * New specification format permits easier compile line overriding...
 */

#ifndef MAX_BC_KEYWORD_LENGTH
#define MAX_BC_KEYWORD_LENGTH 30
#endif

#ifndef MAX_CS_KEYWORD_LENGTH
#define MAX_CS_KEYWORD_LENGTH 50 /* CS = Coordinate System */
#endif

#ifndef MAX_BC_INT_DATA
#define MAX_BC_INT_DATA 5
#endif

#ifndef MAX_BC_FLOAT_DATA
#define MAX_BC_FLOAT_DATA 16
#endif

#ifndef MAX_NODES_PER_SIDE
#define MAX_NODES_PER_SIDE 9
#endif

#ifndef MAX_BC_PER_SIDE
#define MAX_BC_PER_SIDE 64
#endif

/* #define MAX_NUM_SS_DUPLICATIONS 2000 */

#ifndef MAX_SS_PER_NODE
#define MAX_SS_PER_NODE 50
#endif

#ifndef MAX_MAT_PER_SS
#define MAX_MAT_PER_SS 20
#endif

#ifndef MAX_BC_TABLES
#define MAX_BC_TABLES 30
#endif

#ifndef MAX_MP_TABLES
#define MAX_MP_TABLES 20
#endif

#ifndef MAX_EXT_TABLES
#define MAX_EXT_TABLES 2
#endif

#ifndef MAX_AC_TABLES
#define MAX_AC_TABLES 2
#endif

/*
 * The hyperbolic VOF fill equation keeps careful track of inlet boundaries.
 */

#ifndef MAX_INLET_BC
#define MAX_INLET_BC 10
#endif

#ifndef MAX_NUM_NS_NO_XFEM
#define MAX_NUM_NS_NO_XFEM 10
#endif

/* penalty parameter for applying dirichlet or distinguishing conditions */
#ifndef BIG_PENALTY
#define BIG_PENALTY 1.e+12
#endif

#ifndef LITTLE_PENALTY
#define LITTLE_PENALTY 1.e+6
#endif

#ifndef DIRICHLET_PENALTY
#define DIRICHLET_PENALTY 1.
#endif

/*
 * Set some threshhold values for the number of nodesets, number of sidesets
 * and number of boundary conditions. If any of these threshholds is exceeded,
 * the 2D problem is deemed sufficiently complex that the output diagnostics
 * are diverted to a file, also named here.
 */
#ifndef DUP_THRESHHOLD_NODESETS
#define DUP_THRESHHOLD_NODESETS 5
#endif

#ifndef DUP_THRESHHOLD_SIDESETS
#define DUP_THRESHHOLD_SIDESETS 5
#endif

#ifndef DUP_THRESHHOLD_NUMBCS
#define DUP_THRESHHOLD_NUMBCS 7
#endif

#ifndef DUP_THRESHHOLD_FILENAME
#define DUP_THRESHHOLD_FILENAME "BCdup.txt"
#endif

#ifndef BC_3D_INFO_FILENAME
#define BC_3D_INFO_FILENAME "bc3D_output.txt"
#endif

/*   max number of CA cards +1  */

#ifndef MAX_CA
#define MAX_CA 10
#endif

/*
 *       define catagories for methods of boundary condition application
 *       these have been recently changed (10/31/96 RAC) to account for
 *       a wider variety of options in 3D
 */
#define WEAK_INT_SURF      1
#define WEAK_INT_EDGE      2
#define WEAK_POINT         3
#define WEAK_SHIFT         4
#define DIRICHLET          5
#define STRONG_INT_SURF    6
#define STRONG_INT_EDGE    7
#define COLLOCATE_SURF     8
#define COLLOCATE_EDGE     9
#define COLLOCATE_POINT    10
#define COLLOCATE          11
#define SPECIAL            12
#define CONTACT_SURF       13
#define EMBEDDED_SURF      14
#define WEAK_SHELL_GRAD    15
#define STRONG_SHELL_GRAD  16
#define WEAK_SHARP_INT     17
#define STRONG_SHARP_INT   18
#define LS_SPECIAL         19
#define STRONG_INT_NEDELEC 20
#define WEAK_INT_NEDELEC   21

/* define some other catagories */
#define STRESS 6 /* Six components in each mode */
#define VECTOR 3
#define SCALAR 1
#define STRESS 6

#define NO_ROT -1

/* define the applicability of the BC for different phases */
#define SINGLE_PHASE 1
#define CROSS_PHASE                           \
  2 /* for continuous variables but across    \
     * boundaries of conjugate problems, with \
     * different variables in each phase */
#define CROSS_PHASE_DISCONTINUOUS                  \
  3 /* for multivalued variables across boundaries \
     * of different phases */

/*
 * Possible Values for the DVI_Indexing_Type field
 */
/*   DVI_SINGLE_PHASE_DB
 *      These boundary conditions are simple boundary conditions applied to one
 *      variable specified by the equation field in the BC_description structure.
 *      carried out on domain boundaries. Integrations are only carried out on
 *      one side of an element face. DV_Indexing_Matid contains the current material
 *      ID. This is the default for SINGLE_PHASE boundary conditions.
 */
#define DVI_SINGLE_PHASE_DB SINGLE_PHASE
/*
 *   DVI_CROSS_PHASE_CONJGATE
 *      These are boundary conditions applied on material interfaces which
 *      have different variable types active in each material. The boundary
 *      condition involves a strongly integrated (or collocated) boundary
 *      condition involving variables that aren't active on both sides of the
 *      interface. DVI_Indexing_MatID refers to the material which owns the
 *      variables whose dof is receiving this boundary condition. Note, these
 *      boundary conditions may have additional complications due to the fact
 *      that the basis function may not be active on both sides of the
 *      interface. This is the default for CROSS_PHASE boundary conditions.
 */
#define DVI_CROSS_PHASE_CONJUGATE CROSS_PHASE
/*
 *   DVI_SIDTIE
 *      These are boundary conditions that are applied on boundaries with
 *      discontinuous variables. These are the strongly integrated
 *      Dirichlet boundary conditions that compliment the DVI_DVVSIG
 *      boundary conditions. Note, usually a DVI_SIDTIE and a DVI_DVVSIG
 *      boundary condition are paried up with each other. The DVI_Intexing_MatID
 *      variable contains the material ID on whose variable this boundary
 *      condition is applied. Of course, the paired DVI_SIDTIE and DVI_DVVSIG
 *      boundary conditions shouldn't have the same values of DVI_Indexing_Matid.
 *      This type is the default for CROSS_PHASE_DISCONTINUOUS boundary
 *      conditions.
 */
#define DVI_SIDTIE CROSS_PHASE_DISCONTINUOUS
/*
 *   DVI_DVVSIG
 *      These are discontinuous-variable volumetric surface integra galerking
 *      boundary conditions. These are WIC bondary conditions applied on
 *      internal boundaries where the surface flux integral from one phase,
 *      specified by DVI_Indexing_MatID,s replaced with an expression
 *      involving surface and volume integrals from the other phase across
 *      the interface. Thus, the volume contribution from the other phase
 *      is added into the residual expression for the variable from the
 *      phase, specified by DVI_Indexing_MatID variable. DVI_Indexing_MatID
 *      will be set to the MatID with the lowest index value. There may
 *      or may not be an additional surface integral arising out of this
 *      boundary condition
 */
#define DVI_DVVSIG 11
/*
 *   DVI_VSIG
 *      These are volumetric surface integral galerking boundary conditions
 *      applied to a variable type which isn't multi valued at the
 *      interface in contrast to the DVI_DVVSIG case above. They are WIC
 *      boundary conditions applied on internal boundaries where
 *      the surface flux integral from one phase, specified by the
 *      DVI_Indexing_MatID is replaced by an expression involving surface
 *      and volume integrals from the other phase across the interface.
 *      Thus, the volume contribution from the other phase is added
 *      into the residual expression for the variable from the phase,
 *      specifed by DVI_Indexing_MatID variable. DVI_Indexing_MatID will
 *      be set to the MatID with the lowest index value. There may or may
 *      not be an additional surface integral arising out of this
 *      boundary condition.
 */
#define DVI_VSIG 12
/*
 *   DVI_SID
 *      These are strongly integrated Dirichlet conditions that are not
 *      used as tie boundary conditions. Only one phase is involved.
 *      DV_Indexing_MatID contains the index of that phase.
 */
#define DVI_SID 13
/*
 *   DVI_MULTI_PHASE_SINGLE
 *      Boundary conditions which involve the addition of a surface
 *      integral to each side of an internal boundary. However, the
 *      variable is continuous across the interface DV_Indexing_MatID
 *      refers to the "defining" side of the interface. In other
 *      words, it identifies on side of the interface for situations
 *      where different surface integral boundary conditions are being
 *      applied. While the applied variable is continuous across the
 *      interface, other variables at the same interface may or may
 *      not be continuous at that interface. When evaluating integrals
 *      on each side of the interface, variables pertaining to the
 *      current material are used. KINEMATIC_DISC is an example of
 *      this case.
 */
#define DVI_MULTI_PHASE_SINGLE 14
/*
 *  DVI_MULTI_PHASE_VD
 *      Boundary conditions which involve the addition of a surface
 *      integral to each side of an internal boundary between adjacent
 *      materials. The variable is not continuous across the interface.
 *      DVI_Indexing_MatID refers to the "defining" side of the interface.
 *      In other words, it identifies one side of the interface situations
 *      where different surface integral boundary conditions are being applied.
 */
#define DVI_MULTI_PHASE_VD 15

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

/*
 *  Boundary Condition  Step 1: add Macro Substitution for your new boundary
 *                              condition in rf_bc_const.h - this is how you
 *      rf_bc_const.h           will refer to the new boundary condition
 *                              throughout the code.  Make sure the integer
 *                              you choose is unique from all the other BC
 *                              types.
 *  This should be done below:
 */

/* New section here to start adding keywords to boundary conditions */

#define B_USER 1 /* indicates Boundary condition user models*/
#define CIRCLE 2 /* just one model of CA_EDGE_OR_FIX bc     */

/*
 * General Boundary Conditions that may apply to any equation
 * (equation type is specified in the input file
 *  as well as variable type)
 */
#define FIX_BC             200000
#define PERIODIC_BC        210000
#define GD_CONST_BC        300000
#define GD_LINEAR_BC       310000
#define GD_INVERSE_BC      310005
#define GD_PARAB_BC        320000
#define GD_PARAB_OFFSET_BC 320050
#define GD_CIRC_BC         330000
#define GD_POLYN_BC        340000
#define GD_TABLE_BC        350000
#define GD_TIME_BC         400000
#define GD_TIME_LIN        400001
#define GD_TIME_EXP        400002
#define GD_TIME_SIN        400003
#define GD_TIME_TABLE      400004
#define GD_TIME_MAX        400005

/* velocity constants in cartesian form */

#define U_BC            1
#define UVARY_BC        3
#define U_PARABOLA_BC   35
#define UUSER_BC        4
#define UUSER_COLLOC_BC 5

#define PU_BC 6

#define V_BC            10
#define VVARY_BC        30
#define V_PARABOLA_BC   305
#define VUSER_BC        40
#define VUSER_COLLOC_BC 50

#define PV_BC 60

#define USTAR_BC 989991
#define VSTAR_BC 989992
#define WSTAR_BC 989993

#define W_BC 100

#define WVARY_BC        300
#define W_PARABOLA_BC   3005
#define WUSER_BC        400
#define WUSER_COLLOC_BC 500

#define PW_BC 600

#define DX_USER_BC 710
#define DY_USER_BC 720
#define DZ_USER_BC 730

#define DX_USER_NODE_BC 711
#define DY_USER_NODE_BC 721
#define DZ_USER_NODE_BC 731

/* velocity constants in normal/tangential form */

#define N_BC      1000
#define NSIDE_BC  2000
#define NVARY_BC  3000
#define NREACT_BC 4000
/* #define TNRMLSIDE_BC 5000 */
/* #define TNRMLVARY_BC 6000 */
#define SSIDE_BC  20000
#define SVARY_BC  30000
#define SREACT_BC 40000
/* #define TSHRSIDE_BC  50000 */
/* #define TSHRVARY_BC  60000 */

/* polymer stress */

#define S11_BC              1001
#define S12_BC              1002
#define S13_BC              1003
#define S22_BC              1004
#define S23_BC              1005
#define S33_BC              1006
#define U_VES11_PARABOLA_BC 10035
#define U_VES12_PARABOLA_BC 10036
#define U_VES22_PARABOLA_BC 10037
#define U_VES13_PARABOLA_BC 10038
#define U_VES23_PARABOLA_BC 10039
#define U_VES33_PARABOLA_BC 10040

#define S11_1_BC              11001
#define S12_1_BC              11002
#define S13_1_BC              11003
#define S22_1_BC              11004
#define S23_1_BC              11005
#define S33_1_BC              11006
#define U_VES11_1_PARABOLA_BC 11035
#define U_VES12_1_PARABOLA_BC 11036
#define U_VES22_1_PARABOLA_BC 11037
#define U_VES13_1_PARABOLA_BC 11038
#define U_VES23_1_PARABOLA_BC 11039
#define U_VES33_1_PARABOLA_BC 11040

#define S11_2_BC              21001
#define S12_2_BC              21002
#define S13_2_BC              21003
#define S22_2_BC              21004
#define S23_2_BC              21005
#define S33_2_BC              21006
#define U_VES11_2_PARABOLA_BC 21035
#define U_VES12_2_PARABOLA_BC 21036
#define U_VES22_2_PARABOLA_BC 21037
#define U_VES13_2_PARABOLA_BC 21038
#define U_VES23_2_PARABOLA_BC 21039
#define U_VES33_2_PARABOLA_BC 21040

#define S11_3_BC              31001
#define S12_3_BC              31002
#define S13_3_BC              31003
#define S22_3_BC              31004
#define S23_3_BC              31005
#define S33_3_BC              31006
#define U_VES11_3_PARABOLA_BC 31035
#define U_VES12_3_PARABOLA_BC 31036
#define U_VES22_3_PARABOLA_BC 31037
#define U_VES13_3_PARABOLA_BC 31038
#define U_VES23_3_PARABOLA_BC 31039
#define U_VES33_3_PARABOLA_BC 31040

#define S11_4_BC              41001
#define S12_4_BC              41002
#define S13_4_BC              41003
#define S22_4_BC              41004
#define S23_4_BC              41005
#define S33_4_BC              41006
#define U_VES11_4_PARABOLA_BC 41035
#define U_VES12_4_PARABOLA_BC 41036
#define U_VES22_4_PARABOLA_BC 41037
#define U_VES13_4_PARABOLA_BC 41038
#define U_VES23_4_PARABOLA_BC 41039
#define U_VES33_4_PARABOLA_BC 41040

#define S11_5_BC              51001
#define S12_5_BC              51002
#define S13_5_BC              51003
#define S22_5_BC              51004
#define S23_5_BC              51005
#define S33_5_BC              51006
#define U_VES11_5_PARABOLA_BC 51035
#define U_VES12_5_PARABOLA_BC 51036
#define U_VES22_5_PARABOLA_BC 51037
#define U_VES13_5_PARABOLA_BC 51038
#define U_VES23_5_PARABOLA_BC 51039
#define U_VES33_5_PARABOLA_BC 51040

#define S11_6_BC              61001
#define S12_6_BC              61002
#define S13_6_BC              61003
#define S22_6_BC              61004
#define S23_6_BC              61005
#define S33_6_BC              61006
#define U_VES11_6_PARABOLA_BC 61035
#define U_VES12_6_PARABOLA_BC 61036
#define U_VES22_6_PARABOLA_BC 61037
#define U_VES13_6_PARABOLA_BC 61038
#define U_VES23_6_PARABOLA_BC 61039
#define U_VES33_6_PARABOLA_BC 61040

#define S11_7_BC              71001
#define S12_7_BC              71002
#define S13_7_BC              71003
#define S22_7_BC              71004
#define S23_7_BC              71005
#define S33_7_BC              71006
#define U_VES11_7_PARABOLA_BC 71035
#define U_VES12_7_PARABOLA_BC 71036
#define U_VES22_7_PARABOLA_BC 71037
#define U_VES13_7_PARABOLA_BC 71038
#define U_VES23_7_PARABOLA_BC 71039
#define U_VES33_7_PARABOLA_BC 71040

#define STRESS_DEVELOPED_BC 81001

/* velocity gradient */

#define G11_BC 1007
#define G12_BC 1008
#define G13_BC 1009
#define G21_BC 1010
#define G22_BC 1011
#define G23_BC 1012
#define G31_BC 1013
#define G32_BC 1014
#define G33_BC 1015

/* bond evolution */

#define NN_BC 1016

/* extension velocity */

#define EXT_V_BC 1017

/* voltage */

#define VOLT_USER_BC        1018
#define VOLT_BC             1020
#define CURRENT_BC          1021
#define CURRENT_USER_BC     1022
#define CURRENT_BV_BC       1023
#define CURRENT_BV2_BC      1024 /* RSL 1/15/01 */
#define CURRENT_NI_BC       1025 /* RSL 3/9/01 */
#define QS_BC               1026
#define CURRENT_SIC_BC      1027
#define CURRENT_ORR_BC      1028
#define CURRENT_HOR_BC      1029
#define CURRENT_USER_SIC_BC 1037

/* fill */

#define F_BC               1030
#define F_DIODE_BC         1031
#define FILL_INLET_BC      1040
#define FILL_CA_BC         1041
#define STRONG_FILL_CA_BC  1042
#define WETTING_TENSION_BC 1043
#define H_BC               1044
#define LS_INLET_BC        1045
#define F1_BC              1046
#define F2_BC              1047
#define F3_BC              1048
#define F4_BC              1049
#define F5_BC              1050
#define APR_BC             1051
#define API_BC             1052
#define INTP_BC            1053
#define INTM_BC            1054
#define INTD_BC            1055

#define EM_E1R_BC                1056
#define EM_E2R_BC                1057
#define EM_E3R_BC                1058
#define EM_E1I_BC                1059
#define EM_E2I_BC                1060
#define EM_E3I_BC                1061
#define EM_H1R_BC                1062
#define EM_H2R_BC                1063
#define EM_H3R_BC                1064
#define EM_H1I_BC                1065
#define EM_H2I_BC                1066
#define EM_H3I_BC                1067
#define EM_ER_FARFIELD_DIRECT_BC 1068
#define EM_EI_FARFIELD_DIRECT_BC 1069
#define EM_HR_FARFIELD_DIRECT_BC 1070
#define EM_HI_FARFIELD_DIRECT_BC 1071

#define EM_ER_FREE_BC 1072
#define EM_EI_FREE_BC 1073
#define EM_HR_FREE_BC 1074
#define EM_HI_FREE_BC 1075

#define EM_ER_SOMMERFELD_BC 1084
#define EM_EI_SOMMERFELD_BC 1085
#define EM_HR_SOMMERFELD_BC 1086
#define EM_HI_SOMMERFELD_BC 1087

#define E_ER_PLANEWAVE_BC 1200
#define E_EI_PLANEWAVE_BC 1201

#define EM_CONT_REAL_BC 1202
#define EM_CONT_IMAG_BC 1203

#define E_ER_2D_BC 1204
#define E_EI_2D_BC 1205

#define EM_MMS_SIDE_BC       1206
#define EM_MMS_SIDE_IMAG_BC  1207
#define EM_ABSORBING_REAL_BC 1208
#define EM_ABSORBING_IMAG_BC 1209

#define E_ER_FARFIELD_BC 1210
#define E_EI_FARFIELD_BC 1211

#define EM_ER_MMS_BC 1212

#define EM_FARFIELD_REAL_NED_BC 1213
#define EM_FARFIELD_IMAG_NED_BC 1214

/* pressure */

#define P_BC     100000
#define PSPG_BC  100001
#define PSTAR_BC 989994

/* thermal constants */

#define T_BC                 1000000
#define T_USER_BC            3000000
#define QSIDE_BC             4000000
#define QSIDE_DIR_BC         4000001
#define Q_VELO_SLIP_BC       4100000
#define QVARY_BC             5000000
#define QUSER_BC             5600000
#define Q_LASER_WELD_BC      5600001
#define Q_VAPOR_BC           5600002
#define Q_RAIL_BC            5600003
#define QCONV_BC             6000000
#define QRAD_BC              7000000
#define QRAD_REPULSE_ROLL_BC 7100000
#define QREACT_BC            8000000
#define QNOBC_BC             8100000
#define T_MELT_BC            9000000
#define QSIDE_LS_BC          4000002

/* Acoustic wave bcs	*/

#define APR_PLANE_TRANS_BC 9100000
#define API_PLANE_TRANS_BC 9200000
#define APR_VELOCITY_BC    9130000
#define API_VELOCITY_BC    9140000
#define APR_NOBC_BC        9300000
#define API_NOBC_BC        9400000
#define LIGHTP_TRANS_BC    9500000
#define LIGHTM_TRANS_BC    9600000
#define LIGHTD_TRANS_BC    9700000
#define LIGHTP_JUMP_BC     9510000
#define LIGHTM_JUMP_BC     9610000
#define LIGHTP_JUMP_2_BC   9510001
#define LIGHTM_JUMP_2_BC   9610001

/* species unknown variables */

#define Y_BC                    10000000
#define YUSER_BC                11000000
#define SURFACE_CHARGE_BC       12000000
#define FICK_CHRGD_SURF_GRAD_BC 13000001
#define YFLUX_BC                20000000
#define YFLUX_DISC_RXN_BC       20000001
#define YTOTALFLUX_CONST_BC     20500000 /*  RSL 6/7/00  */
#define YFLUX_CONST_BC          21000000
#define YFLUX_USER_BC           22000000
#define YFLUX_SUS_BC            23000000
#define YFLUX_EQUIL_BC          24000000
#define YFLUX_BV_BC             24100000
#define YFLUX_HOR_BC            24400000
#define YFLUX_ORR_BC            24500000
#define YFLUX_H2O_ANODE_BC      24600000
#define YFLUX_H2O_CATHODE_BC    24900000
#define YFLUX_SULFIDATION_BC    24110000
#define YFLUX_ALLOY_BC          24200000
#define YFLUX_BV2_BC            24130000 /* RSL 3/9/01 */
#define YFLUX_NI_BC             24120000 /* RSL 3/9/01 */
#define YFLUX_ETCH_BC           24111111
#define RAOULT                  24300000
#define FLORY                   24700000
#define FLORY_CC                24800000
#define YREACT_BC               30000000

/* Shear rate conditions */
#define SH_BC 310000001

/* porous media boundary conditions */
#define POROUS_FLUX_BC              25000000
#define POROUS_CONV_BC              25100000
#define POROUS_PRESSURE_BC          25300000
#define POROUS_PRESSURE_LUB_BC      26000002
#define DARCY_CONTINUOUS_BC         25400000
#define DARCY_LUB_BC                25400001
#define POROUS_LIQ_PRESSURE_BC      25500000
#define POROUS_LIQ_FLUX_CONST_BC    25600000
#define POROUS_GAS_PRESSURE_BC      25700000
#define POROUS_GAS_FLUX_CONST_BC    25800000
#define POROUS_LIQ_PRESSURE_FILL_BC 25900000
#define POR_LIQ_FLUX_FILL_BC        25900500
#define POROUS_TEMP_BC              26000000
#define P_LIQ_USER_BC               26000001
#define POROUS_SINK_BC              26000003

/* real solid displacement */
#define DX_RS_BC 40000001
#define DY_RS_BC 70000001
#define DZ_RS_BC 100000001

/* mesh displacement */

#define DX_BC         40000000
#define DXVARY_BC     60000000
#define DY_BC         70000000
#define DYVARY_BC     90000000
#define DZ_BC         100000000
#define DZVARY_BC     300000000
#define DX_NOTHING_BC 40100000
#define DY_NOTHING_BC 40200000
#define DZ_NOTHING_BC 40300000

#define FORCE_BC                45000000
#define FORCE_SIC_BC            45020000
#define FORCE_USER_BC           45100000
#define FORCE_USER_SIC_BC       45120000
#define FORCE_RS_BC             45200000
#define FORCE_USER_RS_BC        45400000
#define REP_FORCE_BC            45500000
#define REP_FORCE_RS_BC         45600000
#define ATTR_FORCE_BC           45550000
#define ATTR_FORCE_RS_BC        45650000
#define REP_FORCE_ROLL_BC       45555000
#define REP_FORCE_ROLL_RS_BC    45655000
#define NORM_FORCE_BC           46000000
#define NORM_FORCE_RS_BC        46100000
#define FRICTION_BC             46500000
#define FRICTION_RS_BC          46600000
#define FRICTION_ACOUSTIC_BC    46620000
#define FRICTION_ACOUSTIC_RS_BC 46660000
#define REP_FORCE_SHU_BC        46680000
#define REP_FORCE_SHU_SIC_BC    46690000
#define SLOPE_BC                55000000
#define SLOPEX_BC               55100000
#define SLOPEY_BC               55200000
#define SLOPEZ_BC               55300000

#define DNORMALX_BC   400000000
#define DNORMALY_BC   500000000
#define DNORMALZ_BC   600000000
#define DISTNG_BC     710000000
#define DXDISTNG_BC   700000000
#define DYDISTNG_BC   800000000
#define DZDISTNG_BC   900000000
#define SPLINEX_BC    910000000
#define SPLINEY_BC    920000000
#define SPLINEZ_BC    930000000
#define SPLINE_BC     931000000
#define SPLINEX_RS_BC 911000000
#define SPLINEY_RS_BC 921000000
#define SPLINEZ_RS_BC 932000000
#define SPLINE_RS_BC  931100000
/* replace SPLINE boundary conditions by GEOM */
#define GEOMX_BC                   910000000
#define GEOMY_BC                   920000000
#define GEOMZ_BC                   930000000
#define GEOM_BC                    931000000
#define PLANEX_BC                  940000000
#define PLANEY_BC                  950000000
#define PLANEZ_BC                  960000000
#define PLANE_BC                   961000000
#define FILLET_BC                  961123400
#define DOUBLE_RAD_BC              961123500
#define FEATURE_ROLLON_BC          961223400
#define ROLL_FLUID_BC              961124500
#define TENSION_SHEET_BC           96210200
#define MOVING_PLANE_BC            96110000
#define MOVING_PLANE_ETCH_BC       96115000
#define SM_PLANE_BC                961200000 /* Solid Model PLANE BC */
#define MESH_CONSTRAINT_BC         961300000
#define KINEMATIC_BC               962000000
#define LUB_KINEMATIC_BC           962050000
#define KIN_LEAK_BC                962100000
#define KIN_ELECTRODEPOSITION_BC   962200000 /* RSL 5/24/02 */
#define KINEMATIC_EDGE_BC          962400000
#define KIN_CHEM_BC                962500000
#define KIN_LEAK_HEAT_BC           962600000
#define CAPILLARY_BC               963000000
#define CAPILLARY_SHEAR_VISC_BC    963000002
#define CAPILLARY_TABLE_BC         963000007
#define ELEC_TRACTION_BC           963000001 /* Include Maxwell Stress in CAPILLARY_BC */
#define ELEC_TRACTION_SOLID_BC     963000003 /* Include Maxwell Stress in solid */
#define CAP_REPULSE_BC             963100000
#define CAP_RECOIL_PRESS_BC        963110000
#define CAP_REPULSE_ROLL_BC        963120000
#define CAP_REPULSE_USER_BC        963130000
#define CAP_REPULSE_TABLE_BC       963140000
#define PRESSURE_USER_BC           963200000
#define FLOW_PRESSURE_BC           963300000
#define FLOW_HYDROSTATIC_BC        963400000
#define FLOW_PRESSURE_VAR_BC       963400001
#define FLOW_PRESS_USER_BC         963500000
#define HYDROSTATIC_SYMM_BC        963600000
#define FLOW_STRESSNOBC_BC         963700000
#define FLOW_GRADV_BC              963800000
#define FLOW_GRADV_SIC_BC          963800001
#define FLOW_GRADV_T_BC            963800002
#define SHEET_ENDSLOPE_BC          963900000
#define CA_BC                      964000000
#define CA_MOMENTUM_BC             964000008
#define CA_EDGE_BC                 964000001
#define CA_EDGE_INT_BC             964000002
#define CA_OR_FIX_BC               964000003
#define VAR_CA_EDGE_BC             964000004
#define VAR_CA_USER_BC             964000005
#define CA_EDGE_CURVE_BC           964000006
#define CA_EDGE_CURVE_INT_BC       964000007
#define SURFTANG_BC                964100000
#define SURFTANG_SCALAR_BC         964110000
#define CAP_ENDFORCE_BC            964100010
#define CAP_ENDFORCE_SCALAR_BC     964110010
#define SURFTANG_EDGE_BC           964100001
#define SURFTANG_SCALAR_EDGE_BC    964100002
#define VELO_NORMAL_BC             964200000
#define VELO_NORMAL_EDGE_INT_BC    964200001
#define VELO_NORMAL_EDGE_BC        964200002
#define VELO_NORM_COLLOC_BC        964200003
#define VELO_TANG1_COLLOC_BC       964200010
#define VELO_TANG2_COLLOC_BC       964200011
#define VELO_NORMAL_SOLID_BC       964200004
#define VELO_NORMAL_LS_BC          964200005
#define VELO_NORMAL_LS_COLLOC_BC   964200006
#define VELO_NORMAL_LS_PETROV_BC   964200007
#define VELO_NORMAL_LUB_BC         964200008
#define VNORM_LEAK_BC              964210000
#define VNORM_ELECTRODEPOSITION_BC 964220000 /* RSL 5/30/02 */
#define VELO_TANGENT_BC            964300000
#define VELO_TANGENT_EDGE_INT_BC   964300001
#define VELO_TANGENT_EDGE_BC       964300002
#define VELO_TANGENT_SOLID_BC      964300003
#define VELO_TANGENT_3D_BC         964300004
#define VELO_TANGENT_USER_BC       964300005
#define VELO_TANGENT_LS_BC         964300006
#define ZERO_VELO_TANGENT_3D_BC    964300007
#define VELO_SLIP_BC               964400000
#define VELO_SLIP_ROT_BC           964410000
#define VELO_SLIP_FILL_BC          964420000
#define VELO_SLIP_LEVEL_BC         964421000
#define VELO_SLIP_LEVEL_SIC_BC     964421010
#define VELO_SLIP_LS_ROT_BC        964421328
#define VELO_SLIP_LS_HEAVISIDE_BC  964421009
#define VELO_SLIP_LS_ORIENTED_BC   964421011
#define VELO_SLIP_SOLID_BC         964430000
#define VELO_SLIP_ROT_FILL_BC      964435000
#define VELO_SLIP_EK_BC            964440000
#define VELO_SLIP_POWER_BC         964400001
#define VELO_SLIP_POWER_CARD_BC    964400002

#define VELO_EK_3D_BC          964450000
#define MOVING_CA_BC           964500000
#define CA_EDGE_OR_FIX_BC      964500001
#define LS_ATTACH_BC           964500002
#define VELO_THETA_TPL_BC      964500003
#define VELO_THETA_HOFFMAN_BC  964500004
#define VELO_THETA_COX_BC      964500005
#define WETTING_SPEED_LIN_BC   964500006
#define VELO_THETA_SHIK_BC     964500008
#define LINEAR_WETTING_SIC_BC  964500009
#define VELO_STREAMING_BC      964600000
#define HYSTERESIS_WETTING_BC  964700000
#define AIR_FILM_BC            964800000
#define AIR_FILM_ROT_BC        964810000
#define VELO_SLIP_FLUID_BC     964900000
#define VELO_SLIP_ROT_FLUID_BC 964910000

#define EDDY_NU_BC 966666666

/* Structural Shells */
#define SH_K_BC             970000000
#define SH_TENS_BC          970000001
#define SH_X_BC             970000002
#define SH_Y_BC             970000003
#define SH_FLUID_STRESS_BC  970000004
#define SH_SLOPE_X_BC       970000005
#define SH_SLOPE_Y_BC       970000006
#define SHEAR_TO_SHELL_BC   970000009
#define SH_USER_BC          970000011
#define SH_LUBP_BC          970000013
#define SH_LUBP_SOLID_BC    970000014
#define SH_LUBP_SOLID_RS_BC 970000015
#define SH_S11_WEAK_BC      970000016
#define SH_S22_WEAK_BC      970000017

#define SH_SDET_BC       970000020
#define SH_MESH2_WEAK_BC 970000021

/* Shell variables that are not structural shells */
#define SH_GAMMA1_BC            980000001
#define SH_GAMMA1_DERIV_SYMM_BC 980000002
#define SH_GAMMA2_DERIV_SYMM_BC 980000003
#define DVZDR_ZERO_BC           980000004

/* Boundary conditions that apply at interface between two
   different materials */
#define FLUID_SOLID_BC          47000000
#define SOLID_FLUID_BC          48000000
#define FLUID_SOLID_RS_BC       47000100
#define SOLID_FLUID_RS_BC       48000100
#define SOLID_FLUID_CONTACT_BC  48000120
#define FLUID_SOLID_CONTACT_BC  48000121
#define LAGRANGE_NO_SLIP_BC     48000123
#define BAAIJENS_SOLID_FLUID_BC 48000124
#define BAAIJENS_FLUID_SOLID_BC 48000125
#define SOLID_LAGRANGE_MULT_BC  48000126

#define KIN_DISPLACEMENT_BC          48000110
#define KIN_DISPLACEMENT_RS_BC       48000115
#define KIN_DISPLACEMENT_PETROV_BC   48000210
#define KIN_DISPLACEMENT_COLLOC_BC   48000310
#define SHELL_SURFACE_CHARGE_BC      48000135
#define SHELL_SURFACE_CHARGE_SIC_BC  48000145
#define SURFACE_ELECTRIC_FIELD_BC    48000188
#define SURFACE_ACOUSTIC_VELOCITY_BC 48000193
#define POTENTIAL_NOBC_BC            48000199
#define SHELL_DIFF_KINEMATIC_BC      48000211
#define SURFACE_USER_SHELL_BC        48000250
#define SURFACE_LUBRICATION_BC       48000260
#define NO_SLIP_BC                   47000002
#define NO_SLIP_RS_BC                47000022
#define P_EQUIL_BC                   47000003
#define VL_EQUIL_BC                  47000004
#define VL_POLY_BC                   47000035
#define VOLUME                       47000037
#define MASS                         47000039
#define VP_EQUIL_BC                  47000005
#define VN_POROUS_BC                 47000006
#define POROUS_GAS_BC                47000007
#define CONT_TANG_VEL_BC             47000008
#define CONT_NORM_VEL_BC             47000009
#define RAOULTS_LAW_BC               47000010
#define KINEMATIC_SPECIES_BC         47000011
#define CONTACT_RESISTANCE_BC        47000012
#define DISCONTINUOUS_VELO_BC        47000013
#define EVAPORATION                  470000131
#define DISSOLUTION                  470000132
#define KINEMATIC_DISC_BC            47000014
#define VELO_NORMAL_DISC_BC          47000015
#define LATENT_HEAT_BC               47000016
#define LATENT_HEAT_INTERNAL_BC      47000017
#define SOLID_LIQUID                 47000018
#define LIQUID_VAPOR                 47000019
#define KINEMATIC_COLLOC_BC          47000020
#define KINEMATIC_PETROV_BC          47000021
#define Y_DISCONTINUOUS_BC           47000023
#define HEAT_OF_RXN_BC               47000024
#define T_CONTACT_RESIS_BC           47000025
#define T_CONTACT_RESIS_2_BC         47000026

/*
 * HKM Chemkin boundary conditions along surfaces
 *       -> Not sure if bit masks are used so will keep with the 47 theme!
 *             The "10" digit will represent equation unknowns
 *             The "1" digit will  represent variations of the BC
 * ACS DO not use the numbers between 47000999 and 47900000 unless
 * you have consulted with HKM.  This invokes new_way in
 * bc_integ.c.
 */
#define NEW_WAY_LOW_BC    47000999
#define MF_STEFANFLOW_BC  47001000
#define VL_EQUIL_PRXN_BC  47001040
#define IS_EQUIL_PRXN_BC  47001041
#define SF_CHEM_BC        47001030
#define VN_STEFANFLOW_BC  47001001
#define KINEMATIC_SC_BC   47001011
#define SDC_SURFRXN_BC    47001045
#define SDC_KIN_CHEM_BC   47001042
#define SDC_HEATRXN_BC    47001031
#define SDC_STEFANFLOW_BC 47001002
#define SDC_KIN_SF_BC     47001043
#define SDC_KIN_SFV_BC    47001044
#define NEW_WAY_HIGH_BC   47900000

/* level set-based bc */
#define LS_NO_SLIP_BC         47002000
#define LS_SOLID_FLUID_BC     47002001
#define LS_CAPILLARY_BC       47002002
#define LS_CAP_CURVE_BC       47002003
#define LS_CAP_DIV_N_BC       47002004
#define LS_CAP_DIV_S_N_BC     47002005
#define LS_ADC_BC             47002009
#define LS_ADC_OLD_BC         47002007
#define LS_CAP_HYSING_BC      47002010
#define LS_CAP_DENNER_DIFF_BC 47002011
#define LS_STRESS_JUMP_BC     47002012

/* surface normal dirichlet bc's */
#define N1_BC 47002100
#define N2_BC 47002101
#define N3_BC 47002102

/* lagrange multiplier dirichlet bc's */
#define LM1_BC 47002103
#define LM2_BC 47002104
#define LM3_BC 47002105

/* level set XFEM bcs */
#define LS_T_BC               47002301
#define LS_Q_BC               47002302
#define LS_VELO_NORMAL_BC     47002303
#define LS_VELO_TANGENT_BC    47002304
#define LS_FLOW_PRESSURE_BC   47002305
#define LS_CAPILLARY_GHOST_BC 47002306
#define LS_CONT_T_BC          47002307
#define LS_CONT_FLUX_BC       47002308
#define LS_U_BC               47002309
#define LS_V_BC               47002310
#define LS_W_BC               47002311
#define LS_CONT_TRACTION_BC   47002312
#define LS_CONT_VEL_BC        47002313
#define SHARP_CA_2D_BC        47002314
#define LS_QLASER_BC          47002315
#define LS_RECOIL_PRESSURE_BC 47002316
#define LS_EXTV_FLUID_SIC_BC  47002317
#define LS_EIK_KINEMATIC_BC   47002318
#define LS_EXTV_KINEMATIC_BC  47002319
#define LS_EIK_KIN_LEAK_BC    47002320
#define LS_EXTV_KIN_LEAK_BC   47002321
#define LS_QRAD_BC            47002322
#define LS_QVAPOR_BC          47002323
#define LS_YFLUX_BC           47002324
#define LS_EXTV_LATENT_BC     47002360
#define LS_LATENT_HEAT_BC     47002361
#define LS_ACOUSTIC_SOURCE_BC 47002370
#define LS_WALL_ANGLE_BC      49002371

#define SHARP_WETLIN_VELOCITY_BC  48002325
#define SHARP_BLAKE_VELOCITY_BC   48002326
#define SHARP_HOFFMAN_VELOCITY_BC 48002327
#define SHARP_COX_VELOCITY_BC     48002328
#define WETTING_SPEED_BLAKE_BC    48002329
#define WETTING_SPEED_HOFFMAN_BC  48002330
#define WETTING_SPEED_COX_BC      48002331
#define BLAKE_DIRICHLET_BC        48002332
#define SHARP_SHIK_VELOCITY_BC    48002341
#define WETTING_SPEED_SHIK_BC     48002342
#define HOFFMAN_DIRICHLET_BC      48002343
#define COX_DIRICHLET_BC          48002344
#define SHIK_DIRICHLET_BC         48002345
#define BLAKE_DIRICH_ROLL_BC      48002432
#define HOFFMAN_DIRICH_ROLL_BC    48002443
#define COX_DIRICH_ROLL_BC        48002444
#define SHIK_DIRICH_ROLL_BC       48002445
#define PF_CAPILLARY_BC           48003000

#define TABLE_BC      55500001
#define TABLE_WICV_BC 55500002
#define TABLE_WICS_BC 55500003

/* Boundary conditions associated with Lagrange multipliers */

#define LGR_FLOWRATE_BC 66700001

/* Boundary conditions on curvature equations */

#define H_FREE_BC  777000001
#define LS_CA_H_BC 777000002

/* boundary conditions on bulk lubrication equations */
#define LUB_PRESS_BC             777000003
#define SHELL_FILMP_BC           777000004
#define SHELL_FILMH_BC           777000005
#define SHELL_PARTC_BC           777000006
#define SHELL_GRAD_FP_BC         777000007
#define SHELL_GRAD_FP_NOBC_BC    777000008
#define SHELL_GRAD_FH_BC         777000009
#define SHELL_GRAD_FH_NOBC_BC    777000010
#define SHELL_GRAD_PC_BC         777000011
#define SHELL_GRAD_PC_NOBC_BC    777000012
#define SHELL_TEMP_BC            777000013
#define LUBP_SH_FP_MATCH_BC      777000014
#define LUBP_SH_FP_FLUX_BC       777000015
#define GRAD_LUB_PRESS_BC        777000016
#define SHELL_FLOW_DEVELOPED_BC  777000017
#define SHELL_OPEN_PRESS_BC      777000018
#define SHELL_GRAD_TEMP_BC       777000019
#define SHELL_GRAD_TEMP_NOBC_BC  777000020
#define SH_P_OPEN_USER_BC        777000021
#define LUB_PRESS_2_BC           777000022
#define SHELL_OPEN_PRESS_2_BC    777000023
#define LUB_STATIC_BC            777000024
#define LUB_PRESS_HYDROSTATIC_BC 777000025
#define GRAD_LUBP_NOBC_BC        777000026
#define SHELL_LUB_WALL_BC        777000027

#define SHELL_TFMP_PRES_BC           777000030
#define SHELL_TFMP_SAT_BC            777000031
#define SHELL_TFMP_GRAD_S_BC         777000034
#define SHELL_TFMP_FREE_LIQ_BC       777000041
#define SHELL_TFMP_NUM_DIFF_BC       777000042
#define SHELL_TFMP_AVG_PLATE_VELO_BC 777000043
#define SHELL_TFMP_FREE_GAS_BC       777000044

#define SHELL_SAT_1_BC 777000050
#define SHELL_SAT_2_BC 777000051
#define SHELL_SAT_3_BC 777000052

#define RESTIME_BC                   788000030
#define RESTIME_NOBC_BC              788000031
#define RESTIME_GRADSIC_BC           788000032
#define SHELL_LUBRICATION_OUTFLOW_BC 777000053

/* Vectors used for rotations */
#define ROT_NONE -1
#define ROT_N    -2
#define ROT_N2   -3
#define ROT_N3   -4
#define ROT_T    -5
#define ROT_T1   -6
#define ROT_T2   -7
#define ROT_B    -8
#define ROT_S    -9
#define ROT_X    -10
#define ROT_Y    -11
#define ROT_Z    -12
#define ROT_GD   -13

#define ROT_SEED         1
#define ROT_BASIS        2
#define ROT_BASIS_ONCE   3
#define ROT_BASIS_RESEED 4

/*
 * BC_descriptions Structure:
 *
 *  structure that provides a generic description of boundary
 *  conditions. A pointer to this structure is part of the
 *  Boundary_Condition structure.
 *
 *  Definition of Fields:
 *
 *   name1:    string for name of this bc in input deck
 *   name2:    alternate string for name of this bc in input deck
 *   method:   Descriptor of the method by which this condition is applied
 *   BC_Name:  integer which corresponds to this bc name (listed above)
 *   equation: equation type to which this condition is applied
 *   rotate:   flag to indicate if the corresponding equations should be rotated
 *             to the normal, tangent1, tangent2 coordinate system
 *   vector:   flag to indicate if this is a vector condition
 *   sens:     Flags to indicate which variable types this bc is sensitive to
 *   i_apply:  Defines the applicability of the BC for different phases
 *   DV_Index_Default: Default methodology for application of the variable
 *             in the case of discontinuous variables at an interface.
 *
 *     NOTE: It is important to put in all of the variables for each BC
 *           (or set them to default values) otherwise, the
 *           initialization of this array will be screwed up.
 */
struct BC_descriptions {
  char *name1;
  char *name2;
  int method;
  int BC_Name;
  int equation;
  int vector;
  int rotate;
  int sens[MAX_VARIABLE_TYPES];
  int i_apply;
  int DV_Index_Default;
};
/*
 *  BC_Desc: This is an array of initialized BC_descriptions structures
 *  initialized in the mm_names.h file
 */
extern struct BC_descriptions BC_Desc[];
extern int Num_BC_Names;
/*
 *  num_new_BC_Desc: This is an array of BC_Descriptions read in
 *                   from the input deck.
 */
extern int num_new_BC_Desc;
extern struct BC_descriptions **new_BC_Desc;

/*
 *  Equation_Names structure:
 *      This structure is used to keep lists for translating between
 *      character strings and numerical indecise for equation numbers,
 *      variable numbers, and postprocessing conditions
 *
 */
struct Equation_Names {
  char *name1;
  char *name2;
  int Index;
};
typedef struct Equation_Names EQUATION_NAMES_STRUCT;

extern EQUATION_NAMES_STRUCT EQ_Name[];
extern EQUATION_NAMES_STRUCT Var_Name[];
extern EQUATION_NAMES_STRUCT Post_Var_Name[];
extern EQUATION_NAMES_STRUCT Exo_Var_Names[];
extern EQUATION_NAMES_STRUCT Var_Units[];
extern int Num_EQ_Names;
extern int Num_Var_Names;
extern int Num_Post_Var_Names;
extern int Num_Exo_Var_Names;
extern int Num_Var_Units;

struct Data_Table {

  char t_name[3][132]; /* pointer to string with name of variable parameter */
  int t_index[3];      /* integer that serves to identify the t array (may or may not be used)*/
  int columns;         /* the number of columns in table */
  int interp_method;   /* the method used to interpolate between data points */
  int tablelength;     /* the number of pairs in the table */
  char *f_name;        /* pointer to string with name of dependent function */
  int f_index;         /* integer that serves to identify the f array (e.g. f_index = var ) */
  double *t;           /* pointer to array of variable data points (abscissa)*/
  double *t2;          /* pointer to second array of data points (2D abscissa)*/
  double *t3;          /* pointer to third array of data points (3D abscissa)*/
  double *f;           /* pointer to array of function data points (ordinate)
                        *     So f[i] = F( t[i] ) */
  double slope[3];
  int species_eq;
  int ngrid;     /* for 2d tables, the number of grid points in direction 2 */
  int ngrid2;    /* for 3d tables, the number of grid points in directions 1&2 */
  double yscale; /* Scaling value for the y axis */
  double Emin;   /* Minimum modulus value (for FAUX_PLASTICITY */
};

extern int num_BC_Tables;
extern struct Data_Table *BC_Tables[];
extern int num_MP_Tables;
extern struct Data_Table *MP_Tables[];
extern int num_ext_Tables;
extern struct Data_Table *ext_Tables[];
extern int num_AC_Tables;
extern struct Data_Table *AC_Tables[];

/* structure to hold b.c. information entered from the input file
 *-------------------------------------------------------------------------------
 *
 *   The input file should contain boundary condition information in the
 *   following format:
 *
 *     ------------------------------------------------------------
 *             Boundary Condition Specifications
 *     ------------------------------------------------------------
 *     Number of BC	       	 = 1
 *     BC		         = BC_Name Set_Type BC_ID BC_Data
 *
 *   where BC_Name is a string (maximum length of 'MAX_BC_KEYWORD_LENGTH'
 *   which matches one of the types defined in 'rf_bc_const.h', Set_Type is a 2
 *   character string, either 'NS' for node sets or 'SS' for side sets, BC_ID is
 *   the corresponding node-set or side-set ID value from the mesh and BC_Data is
 *   either integer or floating-point data, maximum in number either
 *   'MAX_BC_INT_DATA' or 'MAX_BC_FLOAT_DATA', respectively.  An example of a
 *   Dirichlet-type boundary condition on temperature set at a node-set with ID 2
 *   and a value of 3.14 is:
 *
 *     BC				 = T NS 2 3.14
 *
 *   The following b.c. types have a single floating-point data value:  U, USIDE,
 *   V, VSIDE, W, WSIDE, N, NSIDE, TNRMLSIDE, S, SSIDE, TSHRSIDE, P, T, TSIDE,
 *   QSIDE, DX, DXSIDE, DY, DYSIDE, DZ, DZSIDE , DXDISTNG, DYDISTNG, DZDISTNG.
 *
 *   RAC - adding extra entries to hold key information about each boundary condition
 *
 *   PRS - Added BC_relax to hold relaxation parameter for all dirichlet
 *         conditions if desired
 *
 *   pas - added indexing and sizing variables to assist in interprocessor
 *         communication of pointers to struct BC_description's, static or
 *         dynamically allocated, as well as any u_BC data of whatever length.
 */

struct Boundary_Condition {
  int BC_Name; /* Primary id of the bc -> matches one of the names
                * above */
  char Set_Type[3];
  int BC_ID;             /* Exodus ss or ns ID for the boundary condition */
  int BC_ID2;            /* 2nd Exodus ID needed for EDGE and VERTEX BC's */
  int BC_ID3;            /* 3rd Exodus ID needed for VERTEX BC's */
  int Set_Index;         /* This is the ss or ns index value that matches
                          * the first BC_ID value */
  int BC_matrl_index_1;  /* Material index for nodes on the first side
                            of the boundary  */
  int BC_matrl_index_2;  /* Material index for nodes on the second
                          * side of the boundary -> Pertinent for
                          * CROSS_PHASE and CROSS_PHASE_DISCONTINUOUS
                          * boundary conditions, or EDGE or VERTEX
                          * bc's */
  int BC_matrl_index_3;  /* Pertinent for EDGE and VERTEX BC's */
  int BC_matrl_index_4;  /* Pertinent for VERTEX BC's */
  int Internal_Boundary; /* Boolean. It is true if any part of the interface
                          * is an internal boundary, i.e., has part of
                          * the computational domain on both sides of it. */
  int BC_EBID_Apply;     /* If this boundary condition is restricted to
                          * be applied to variables from
                          * a particular element block,
                          * this entry will contain the element block
                          * ID of that element block. If not, this
                          * entry will be equal to -1.
                          */
  int BC_Data_Int[MAX_BC_INT_DATA];
  dbl BC_Data_Float[MAX_BC_FLOAT_DATA];
  int len_u_BC; /* number of elements in the user constant
                   list (0 most of the time) */
  int max_DFlt;
  int Storage_ID; /* ID of the quadature point storage for this bc
                   * Must be positive for it to exist. zero means
                   * that it does not yet exist
                   */
  double *u_BC;
  struct BC_descriptions *desc;
  int BC_Desc_index; /* transportable thru pseudo structure  */
  int index_dad;     /*
                      * "Dad, what do you do at work?"
                      *
                      * "Child, I work hard as an index substituting
                      *  as a pointer to local dynamic and statically
                      *  allocated memory."
                      *
                      *  index_dad = -1  --> refers to a static
                      *                      BC_Desc, use
                      *                      BC_Desc_index to
                      *                      say which one.
                      *
                      *  index_dad > -1  --> index of a dynamically
                      *                      allocated BC_description
                      *                      that also must be
                      *                      communicated cross proc.
                      *			Still use BC_Desc_index
                      *                      to point to template
                      *                      static BC_Desc since
                      *                      new one has pointers...
                      */
  int species_eq;
  dbl BC_relax; /* relaxation parameter for this dirichlet
                 *  conditions if desired */
  struct Data_Table *table;
  int table_index;
  int DV_Indexing_Type;     /* This field describes to Goma which equation
                             * to apply this boundary condition on when
                             * the equation and i_apply field in the
                             * BC_description structure isn't enough  */
  int DV_Indexing_MatID;    /* Used in conjunction with DV_Indexing_Type
                             * to describe the material id to apply
                             * the current boundary condition on.  */
  int BC_Memory_Allocation; /* Any temporary memory associated with this bc?
                             *  0 - no
                             *  1 - Yes, hanging off of node info structure
                             *  2 - Yes, hanging off of side  */
  int matrix;
  int equation;
};
typedef struct Boundary_Condition BOUNDARY_CONDITION_STRUCT;

#define CDIM 3
/*
 * Structure to hold input information about rotation conditions
 */
struct Rotation_Specs {
  int eq_type;
  int type;
  int ss_id[CDIM];
  int BC_Type[CDIM];
  int BC_SS[CDIM];
  int BC_id[CDIM];
  struct BC_descriptions *BC_desc[CDIM];
  int BC_desc_index[CDIM]; /* help teleport previous ptr cross procs */
  int method;
  int ROTATE; /* flag indicating whether any rotation
                 happens at this node */
  double seed[CDIM];
  int node; /* special integer which stores node number
               for vertex conditions */
  int *elems;
  int num_elem;
  int ss_ptr;
};

/* Structure for processing side boundary conditions:
 *
 *    id_side           = id of the side of the element on which the bc is to be
 *                        applied (1 to 6)
 *    num_nodes_on_side = number of nodes on the bc side
 *    local_node_id[]   = list of local node id for nodes which belong to the
 *		         element side
 * 			 (statically assigned , but real length is equal to
 *		          num_nodes_on_side)
 *    *next_side_bc     = Pointer to the next side of the current element
 *			 that needs an integral boundary condition applied.
 *			 A NULL pointer indicates that there are no more sides.
 *
 */

struct elem_side_bc_struct {
  int ielem;             /* element number */
  int id_side;           /* face index */
  int num_nodes_on_side; /* number of nodes on element face */
  int BC_applied;        /* Index of BC applied to this face */
  int *local_node_id;
  int *local_elem_node_id;
  int *BC_input_id;
  int Num_BC;
  int MatID_List[2];
  int Num_MatID;
  struct QP_Storage **Side_QP_Storage;
  struct elem_side_bc_struct *next_side_bc;
};
typedef struct elem_side_bc_struct ELEM_SIDE_BC_STRUCT;

/*
 *  Structure for processing edge boundary conditions:
 */
struct elem_edge_bc_struct {
  int ielem;
  int id_edge;
  int num_nodes_on_edge;
  int ipin;   /*Keeps track of whether pinned or not */
  int shared; /* TRUE if edged is shared by two elements */
  int BC_applied;
  int local_node_id[MAX_NODES_PER_SIDE];
  int edge_elem_node_id[MAX_NODES_PER_SIDE];
  int BC_input_id[MAX_BC_PER_SIDE];
  struct elem_side_bc_struct *elem_side_bc_1; /* side_bc of primary side */
  struct elem_side_bc_struct *elem_side_bc_2; /* side_bc of secondary side */
  struct elem_edge_bc_struct *next_edge_bc;
};

/* Structure for processing side set bc's on mesh motion:
 *    The idea here is to aid in placing distinguishing condtions on
 *    nodes belonging to two side sets (2D) or  two or three side sets (3D)
 *    so that penalty conditions are not mistakenly put on the same
 *    equation.  By default, the condition corresponding to
 *    the first ss_list entry at a given node in node_list will be assigned
 *    to the x-eqn, the second to the y-eqn., and so on.
 *
 *    node_list = nodes which are contained in more than one side set
 *      ss_list = list of side set id's attached to each node in node list.
 *
 */

/*
 * new array to store list of nodes that exist on SS
 * - for geometrical purposes
 */
extern int *boundary_node_list;
extern int **ss_on_boundary_node_list;
extern int num_boundary_nodes;
extern int **BC_dup_nodes;
extern int ****BC_dup_list;
extern int *BC_dup_ptr;
extern int *ss_to_blks[MAX_MAT_PER_SS + 1];
extern int dup_blks_list[MAX_MAT_PER_SS + 1];
extern int *SS_Internal_Boundary;
extern int **mesh_rotate_node;
extern int **mesh_rotate_ss;
extern int *num_mesh_rotate;
extern int **mom_rotate_node;
extern int **mom_rotate_ss;
extern int *num_mom_rotate;
extern int PRESSURE_DATUM;          /* flag to determine if a pressure datum is set */
extern int pressure_datum_element;  /* element in which the pressure datum is set */
extern double pressure_datum_value; /* value of the pressure datum */

#undef CDIM

#endif
